# -*- coding: utf-8 -*-
"""
SPECT attenuation images from a CT DICOM series (Python port of the SPECT branch of loadCorrections.m, HU_to_mu.m and
HU_rescale_coeff.m).
"""

import os

import numpy as np

# Energies (keV) and bilinear HU-to-mu coefficients, values from https://doi.org/10.1016/j.apradiso.2008.01.002
_ENERGIES = np.array([75., 80., 140., 159., 167., 171., 245., 364., 511.])
_INTERCEPT = 0.1 * np.array([0.16, 0.16, 0.15, 0.14, 0.13, 0.14, 0.12, 0.09, 0.09])
_SLOPE1 = 0.1 * 1e-4 * np.array([1.66, 1.66, 1.52, 1.38, 1.36, 1.39, 1.25, 0.941, 0.96])
_SLOPE2 = 0.1 * 1e-4 * np.array([1.48, 1.48, 1.14, 1.26, 1.25, 1.05, 0.861, 0.559, 0.573])


def HUtoMu(CTvol, keV):
    """
    Converts CT values (HU) to linear attenuation coefficients with a bilinear model (port of HU_to_mu.m).

    Parameters
    ----------
    CTvol : array_like or scalar
        CT values in Hounsfield units.
    keV : float
        Photon energy in keV, 75 to 511. The coefficients are linearly interpolated between the tabulated energies.

    Returns
    -------
    ndarray
        Linear attenuation coefficients in 1/mm (float64), same shape as CTvol.
    """
    keV = float(keV)
    if not (_ENERGIES[0] <= keV <= _ENERGIES[-1]):
        raise ValueError('options.keV = %g is outside the supported range of %g to %g keV!' % (keV, _ENERGIES[0], _ENERGIES[-1]))
    slope1 = np.interp(keV, _ENERGIES, _SLOPE1)
    slope2 = np.interp(keV, _ENERGIES, _SLOPE2)
    intercept = np.interp(keV, _ENERGIES, _INTERCEPT)
    CTvol = np.asarray(CTvol, dtype=np.float64)
    return np.where(CTvol < 0, slope1 * CTvol + intercept, slope2 * CTvol + intercept)


def _readCTSeries(fpathCT):
    """
    Reads all DICOM files (.dcm and .ima) of a folder. The slices are sorted by the z coordinate of ImagePositionPatient
    (ascending) and converted to HU with their own RescaleSlope and RescaleIntercept.

    Each slice uses its own coefficients (MATLAB's HU_rescale_coeff matches the files to the slices by ImagePositionPatient).

    Returns
    -------
    CTvol : ndarray
        (Rows, Columns, nSlices) float64 volume in HU.
    pixelSpacing : ndarray
        DICOM PixelSpacing (row spacing, column spacing) in mm.
    positions : ndarray
        (nSlices, 3) sorted ImagePositionPatient.
    """
    try:
        from pydicom import dcmread
    except ModuleNotFoundError:
        raise ModuleNotFoundError('pydicom package not found! DICOM files cannot be loaded. You can install pydicom package with "pip install pydicom".') from None
    if not os.path.isdir(fpathCT):
        raise ValueError('The CT DICOM folder "%s" does not exist!' % fpathCT)
    files = sorted(f for f in os.listdir(fpathCT) if os.path.splitext(f)[1].lower() in ('.dcm', '.ima')
                   and os.path.isfile(os.path.join(fpathCT, f)))
    if len(files) == 0:
        raise ValueError('No .dcm or .ima files found in "%s"!' % fpathCT)
    slices = []
    for f in files:
        ds = dcmread(os.path.join(fpathCT, f))
        pos = np.array([float(v) for v in ds.ImagePositionPatient])
        slope = float(getattr(ds, 'RescaleSlope', 1.))
        intercept = float(getattr(ds, 'RescaleIntercept', 0.))
        img = ds.pixel_array.astype(np.float64) * slope + intercept
        slices.append((pos, img, np.array([float(v) for v in ds.PixelSpacing])))
    slices.sort(key=lambda s: s[0][2])
    shape = slices[0][1].shape
    for s in slices:
        if s[1].shape != shape:
            raise ValueError('All CT slices need to have the same size!')
    CTvol = np.stack([s[1] for s in slices], axis=2)
    positions = np.array([s[0] for s in slices])
    return CTvol, slices[0][2], positions


def _askDirectory(title, optionName):
    """Asks for a folder with a dialog. optionName is the option that gives the folder without a dialog."""
    from ..reconstruction.prepass import _hasDisplay
    message = 'Cannot show the folder dialog "%s" because no display is available! Set %s instead.' % (title, optionName)
    if not _hasDisplay():
        raise ValueError(message)
    try:
        import tkinter as tk
        from tkinter.filedialog import askdirectory
    except ImportError:
        raise ValueError(message) from None
    try:
        root = tk.Tk()
    except tk.TclError:
        raise ValueError(message) from None
    root.withdraw()
    nimi = askdirectory(title=title)
    root.destroy()
    if len(nimi) == 0:
        raise ValueError('No folder selected!')
    return nimi


def _interpAxis(arr, axis, c):
    """
    Linear interpolation along one axis at the fractional (0-based) indices c. Returns the interpolated array and the mask of
    indices between the outermost voxel centres (0 to n-1). Like imwarp with linear interpolation, points outside the range of
    the voxel centres (including the outer half voxel) are filled by the caller.
    """
    n = arr.shape[axis]
    inside = (c >= -1e-9) & (c <= n - 1 + 1e-9)
    cc = np.clip(c, 0., n - 1.)
    if n == 1:
        return np.take(arr, np.zeros(len(c), dtype=np.int64), axis=axis), inside
    i0 = np.minimum(np.floor(cc).astype(np.int64), n - 2)
    w = cc - i0
    shp = [1] * arr.ndim
    shp[axis] = len(c)
    w = w.reshape(shp)
    return np.take(arr, i0, axis=axis) * (1. - w) + np.take(arr, i0 + 1, axis=axis) * w, inside


def loadSPECTCTAttenuation(options):
    """
    Forms the SPECT attenuation image from a CT DICOM series (port of the SPECT branch of loadCorrections.m). The CT is
    converted to HU, clamped to -1000, converted to linear attenuation coefficients (1/mm) at options.keV and resampled
    (trilinear, air outside the CT) onto the SPECT image grid options.refSPECT.

    Uses options.fpathCT (folder dialog if empty), options.keV (photon energy in keV, required) and options.refSPECT
    (dict with ImageSize, XWorldLimits, YWorldLimits and ZWorldLimits, as set by loadProSpectaData; if None, a centred
    grid from FOVa_x, FOVa_y, axial_fov and Nx, Ny, Nz is used).

    Returns
    -------
    ndarray
        (Nx, Ny, Nz) attenuation image in 1/mm (float64).
    """
    if getattr(options, 'keV', None) is None:
        raise ValueError('options.keV (the photon energy in keV) needs to be set in order to form the attenuation image from CT!')
    if len(getattr(options, 'fpathCT', '')) == 0:
        options.fpathCT = _askDirectory('Select the folder containing attenuation DICOM files', 'options.fpathCT')

    CTvol, pixelSpacing, pos = _readCTSeries(options.fpathCT)
    nSlices = CTvol.shape[2]
    if nSlices < 2:
        raise ValueError('The CT series needs at least two slices!')
    # DICOM PixelSpacing is (row spacing, column spacing); columns run along x and rows along y (as in imref3d)
    dx = pixelSpacing[1]  # Column spacing
    dy = pixelSpacing[0]  # Row spacing
    dz = abs(pos[1, 2] - pos[0, 2])
    nRowsCT = CTvol.shape[0]
    nColsCT = CTvol.shape[1]
    xLimits = pos[0, 0] + np.array([0., dx * nColsCT]) - 0.5 * dx
    yLimits = pos[0, 1] + np.array([0., dy * nRowsCT]) - 0.5 * dy
    zLimits = np.array([pos[:, 2].min(), pos[:, 2].max()]) + dz * np.array([-0.5, 0.5])

    CTvol[CTvol < -1000.] = -1000.
    MUvol = HUtoMu(CTvol, options.keV)
    muAir = float(HUtoMu(-1000., options.keV))

    ref = getattr(options, 'refSPECT', None)
    if ref is None:
        print('options.refSPECT not set, using the centred FOV as the SPECT image grid.')
        Nx0 = int(np.asarray(options.Nx).ravel()[0])
        Ny0 = int(np.asarray(options.Ny).ravel()[0])
        Nz0 = int(np.asarray(options.Nz).ravel()[0])
        fx = float(np.asarray(options.FOVa_x).ravel()[0])
        fy = float(np.asarray(options.FOVa_y).ravel()[0])
        fz = float(np.asarray(options.axial_fov).ravel()[0])
        ref = {'ImageSize': (Nx0, Ny0, Nz0), 'XWorldLimits': (-0.5 * fx, 0.5 * fx), 'YWorldLimits': (-0.5 * fy, 0.5 * fy),
               'ZWorldLimits': (-0.5 * fz, 0.5 * fz)}
    sz = [int(v) for v in ref['ImageSize']]
    xL = np.asarray(ref['XWorldLimits'], dtype=np.float64)
    yL = np.asarray(ref['YWorldLimits'], dtype=np.float64)
    zL = np.asarray(ref['ZWorldLimits'], dtype=np.float64)

    # imref3d convention: array dimension 1 (rows) <-> y, dimension 2 (columns) <-> x, dimension 3 <-> z.
    # World coordinates of the output voxel centres
    yo = yL[0] + (np.arange(sz[0]) + 0.5) * (yL[1] - yL[0]) / sz[0]
    xo = xL[0] + (np.arange(sz[1]) + 0.5) * (xL[1] - xL[0]) / sz[1]
    zo = zL[0] + (np.arange(sz[2]) + 0.5) * (zL[1] - zL[0]) / sz[2]
    # Fractional 0-based CT indices (CT rows <-> y with yLimits, columns <-> x with xLimits, slices <-> z with zLimits)
    cy = (yo - yLimits[0]) / ((yLimits[1] - yLimits[0]) / CTvol.shape[0]) - 0.5
    cx = (xo - xLimits[0]) / ((xLimits[1] - xLimits[0]) / CTvol.shape[1]) - 0.5
    cz = (zo - zLimits[0]) / ((zLimits[1] - zLimits[0]) / CTvol.shape[2]) - 0.5

    out, my = _interpAxis(MUvol, 0, cy)
    out, mx = _interpAxis(out, 1, cx)
    out, mz = _interpAxis(out, 2, cz)
    mask = my[:, None, None] & mx[None, :, None] & mz[None, None, :]
    return np.where(mask, out, muAir)
