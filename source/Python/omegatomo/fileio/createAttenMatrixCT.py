# -*- coding: utf-8 -*-
"""
Conversion of DICOM CT images to 511 keV attenuation images (Python ports of create_atten_matrix_CT.m and attenuationCT_to_511.m).
"""

import os

import numpy as np

# Effective (broad-beam) water attenuation at 511 keV in 1/cm, used without scatter correction
_WATER_MU_EFFECTIVE = 0.093
# Narrow-beam water attenuation at 511 keV in 1/cm, used with scatter correction
_WATER_MU_NARROW = 0.096


def attenuationCTTo511(ctImage, KVP, name=None, waterMu=None):
    """
    Scales CT values (HU) to 511 keV attenuation coefficients with a bilinear (piecewise linear) model. The X-ray tube
    potentials 80, 100, 120 and 140 kVp are tabulated, other values are interpolated (with a warning).

    By default water (0 HU) is mapped to 0.093 cm^-1. This is an effective (broad-beam) value that partly compensates for the
    scattered events that remain in the data when no scatter correction is applied. When scatter correction is used, the
    narrow-beam value 0.096 cm^-1 should be used instead. waterMu scales the whole converted curve by waterMu / 0.093 (air
    stays at 0).

    Parameters
    ----------
    ctImage : ndarray
        CT values in HU, any shape (the shape is kept).
    KVP : float
        The CT kVp value (e.g. 100).
    name : str, optional
        Examination name. If given, the attenuation coefficients are saved in the file
        511keV_attenuation_coefficients_for_<name>.mat (variable attenuation_factors).
    waterMu : float, optional
        The 511 keV attenuation coefficient of water in cm^-1. Default 0.093.

    Returns
    -------
    ndarray
        511 keV attenuation coefficients in cm^-1 (float64), same shape as ctImage.
    """
    if waterMu is None:
        waterMu = _WATER_MU_EFFECTIVE
    else:
        try:
            valid = np.ndim(waterMu) == 0 and bool(np.isfinite(float(waterMu))) and float(waterMu) > 0
        except (TypeError, ValueError):
            valid = False
        if not valid:
            raise ValueError('waterMu must be a positive scalar (attenuation coefficient of water in cm^-1)')
        waterMu = float(waterMu)
    # Values from: Accuracy of CT-based attenuation correction in PET/CT bone imaging
    if KVP == 100:
        a = np.array([9.3, 4., 0.5]) * 1e-5
        b = np.array([0.093, 0.093, 0.128])
    elif KVP == 80:
        a = np.array([9.3, 3.28, 0.41]) * 1e-5
        b = np.array([0.093, 0.093, 0.122])
    elif KVP == 120:
        a = np.array([9.3, 4.71, 0.589]) * 1e-5
        b = np.array([0.093, 0.093, 0.134])
    elif KVP == 140:
        a = np.array([9.3, 5.59, 0.698]) * 1e-5
        b = np.array([0.093, 0.093, 0.142])
    else:
        print('Unsupported kVp, interpolating initial values')
        from scipy.interpolate import CubicSpline
        aa = np.array([[9.3, 3.28, 0.41], [9.3, 4., 0.5], [9.3, 4.71, 0.589], [9.3, 5.59, 0.698]]) * 1e-5
        bb = np.array([[0.093, 0.093, 0.122], [0.093, 0.093, 0.128], [0.093, 0.093, 0.134], [0.093, 0.093, 0.142]])
        c = np.array([80., 100., 120., 140.])
        a = np.zeros(3)
        b = np.zeros(3)
        for kk in range(3):
            a[kk] = CubicSpline(c, aa[:, kk], bc_type='not-a-knot')(float(KVP))
            b[kk] = CubicSpline(c, bb[:, kk], bc_type='not-a-knot')(float(KVP))
    # Intersection points of the piecewise linear curve (b[2] is deliberately unused, as in MATLAB)
    xk = np.array([-1000., 0., 1000., 3000.])
    yk = np.array([b[0] - 1000. * a[0], b[1], b[1] + a[1] * 1000., b[1] + a[1] * 1000. + a[2] * 2000.])
    ct = np.asarray(ctImage, dtype=np.float64)
    attenuationFactors = np.interp(ct, xk, yk)
    # Linear extrapolation outside the intersection points with the slope of the end segments
    slopeLow = (yk[1] - yk[0]) / (xk[1] - xk[0])
    slopeHigh = (yk[3] - yk[2]) / (xk[3] - xk[2])
    low = ct < xk[0]
    high = ct > xk[3]
    attenuationFactors = np.where(low, yk[0] + (ct - xk[0]) * slopeLow, attenuationFactors)
    attenuationFactors = np.where(high, yk[3] + (ct - xk[3]) * slopeHigh, attenuationFactors)
    # Scale the whole curve to the requested water attenuation (the tabulated values correspond to 0.093 cm^-1 for water)
    if waterMu != _WATER_MU_EFFECTIVE:
        attenuationFactors = attenuationFactors * (waterMu / _WATER_MU_EFFECTIVE)
    if name is not None and len(str(name)) > 0:
        from scipy.io import savemat
        savemat('511keV_attenuation_coefficients_for_' + str(name) + '.mat', {'attenuation_factors': attenuationFactors})
    return attenuationFactors


def _readCTSliceHU(fname):
    """
    Reads a DICOM CT slice and converts it to HU: HU = RescaleSlope * value + RescaleIntercept. If either field is missing or
    empty, HU = value - 1000 is used and hasRescale is False. Values below -1000 HU (air) are set to -1000.
    """
    import pydicom
    ds = pydicom.dcmread(fname)
    slope = ds.get('RescaleSlope', None)
    intercept = ds.get('RescaleIntercept', None)
    hasRescale = slope is not None and intercept is not None and str(slope) != '' and str(intercept) != ''
    pixels = np.asarray(ds.pixel_array, dtype=np.float64)
    if hasRescale:
        hu = pixels * float(slope) + float(intercept)
    else:
        hu = pixels - 1000.
    # Values below air (e.g. -1024 HU outside the reconstruction circle in Siemens CT) are not physical and would give
    # negative attenuation coefficients in attenuationCTTo511, which extrapolates linearly
    hu = np.maximum(hu, -1000.)
    return hu, hasRescale


def _overlapMatrix(nIn, psIn, nOut, psOut):
    """
    Area-weighting matrix for 1D resampling. Returns a dense nOut x nIn matrix M, where M[j, k] is the length of the overlap
    between output pixel j and input pixel k divided by psOut. Both grids are centered at the same point. The rows sum to the
    fraction of the output pixel covered by the input grid (1 inside, less than 1 at the edges, 0 outside).
    """
    nIn = int(nIn)
    nOut = int(nOut)
    edgesIn = (2. * np.arange(nIn + 1) - nIn) * psIn / 2.
    edgesOut = (2. * np.arange(nOut + 1) - nOut) * psOut / 2.
    M = np.zeros((nOut, nIn))
    for j in range(nOut):
        # Candidate input pixels (one extra on both sides, the overlap is clipped to zero below). kLo and kHi are 1-based
        kLo = max(1, int(np.floor((edgesOut[j] - edgesIn[0]) / psIn)))
        kHi = min(nIn, int(np.ceil((edgesOut[j + 1] - edgesIn[0]) / psIn)) + 1)
        if kHi < kLo:
            continue
        # 0-based input pixels: the left edge of pixel k is edgesIn[k] and the right edge edgesIn[k + 1]
        k = np.arange(kLo - 1, kHi)
        ov = np.maximum(0., np.minimum(edgesOut[j + 1], edgesIn[k + 1]) - np.maximum(edgesOut[j], edgesIn[k]))
        M[j, k] = ov / psOut
    return M


def _cubicKernel(x):
    """The cubic convolution kernel of the MATLAB imresize 'bicubic' method (a = -0.5)."""
    ax = np.abs(x)
    ax2 = ax * ax
    ax3 = ax2 * ax
    return (1.5 * ax3 - 2.5 * ax2 + 1.) * (ax <= 1.) + (-0.5 * ax3 + 2.5 * ax2 - 4. * ax + 2.) * ((1. < ax) & (ax <= 2.))


def _imresizeSliceAxis(volume, nOut):
    """
    Resizes the third axis of a 3D array to nOut samples as MATLAB imresize does with its defaults (bicubic kernel, antialiasing
    on, symmetric boundary), applied along one dimension.
    """
    volume = np.asarray(volume, dtype=np.float64)
    nIn = volume.shape[2]
    scale = nOut / nIn
    if scale < 1.:
        kernelWidth = 4. / scale
        h = lambda x: scale * _cubicKernel(scale * x)
    else:
        kernelWidth = 4.
        h = _cubicKernel
    # One-based output index j
    j = np.arange(1, nOut + 1, dtype=np.float64)
    u = j / scale + 0.5 * (1. - 1. / scale)
    left = np.floor(u - kernelWidth / 2.)
    P = int(np.ceil(kernelWidth)) + 2
    indices = left[:, None] + np.arange(P)[None, :]
    weights = h(u[:, None] - indices)
    weights = weights / np.sum(weights, axis=1, keepdims=True)
    # Mirror the indices at the borders (one-based indices into [1:in, in:-1:1])
    aux = np.concatenate((np.arange(1, nIn + 1), np.arange(nIn, 0, -1)))
    indices = aux[np.mod(indices.astype(np.int64) - 1, 2 * nIn)] - 1
    R = np.zeros((nOut, nIn))
    rows = np.repeat(np.arange(nOut), P)
    np.add.at(R, (rows, indices.ravel()), weights.ravel())
    return np.tensordot(volume, R, axes=(2, 1))


def createAttenMatrixCT(options, KVP=None, name=None):
    """
    Loads DICOM CT images (.ima or .dcm files) and scales the CT (HU) values to the corresponding 511 keV attenuation
    coefficients with attenuationCTTo511. It is assumed that the PET FOV is centered in the center of the CT FOV.
    The CT images are resampled to the PET image grid (Nx x Ny, FOVa_x x FOVa_y) with exact area weighting. If the PET FOV is
    larger than the CT FOV, the part outside the CT is filled with air. The calculated attenuation coefficients are also saved
    on disk if a name for the examination is provided.

    If options.scatter_correction is True, the narrow-beam water attenuation coefficient 0.096 cm^-1 is used in the conversion
    (the scattered events are modelled separately). Otherwise the effective value 0.093 cm^-1 is used, which partly compensates
    for the scatter remaining in the data. See attenuationCTTo511 for details.

    The stored pixel values are converted to HU slice by slice with the DICOM RescaleSlope and RescaleIntercept
    (HU = slope * value + intercept). If a file contains no rescale information, the fixed offset 1000 is subtracted instead
    (HU = value - 1000) and a warning is printed once. HU values below -1000 (air) are set to -1000.

    The folder options.attenuation_datafile must contain a single CT series (all .ima and .dcm files of the folder are read, in
    the order of their names). The pydicom package is required.

    Parameters
    ----------
    options : projectorClass
        The options (Nx, Ny, Nz, axial_fov, FOVa_x, FOVa_y, attenuation_datafile and, optionally, scatter_correction and
        verbose are used). FOVa_y is optional, FOVa_x is used if it is missing or empty.
    KVP : float, optional
        The CT kVp value (e.g. 100). If omitted, the KVP of the DICOM header is used.
    name : str, optional
        Examination name. If given, the 511 keV attenuation coefficients (cm^-1) are saved in a MAT-file named with it.

    Returns
    -------
    ndarray
        The 511 keV attenuation image in 1/mm, float64, shape (Nx, Ny, Nz) in Fortran order. The first dimension corresponds to
        the DICOM image rows and the third to the sorted slices, so it can be input as options.vaimennus (with
        options.attIncm = False).
    """
    try:
        import pydicom
    except ModuleNotFoundError:
        raise ModuleNotFoundError('pydicom package not found! DICOM CT images cannot be loaded. You can install pydicom package with "pip install pydicom".') from None
    Nx = int(np.atleast_1d(options.Nx)[0])
    Ny = int(np.atleast_1d(options.Ny)[0])
    Nz = int(np.atleast_1d(options.Nz)[0])
    axialFOV = float(np.atleast_1d(options.axial_fov)[0])
    FOVaX = float(np.atleast_1d(options.FOVa_x)[0])
    FOVaY = getattr(options, 'FOVa_y', None)
    if FOVaY is None or np.size(FOVaY) == 0:
        FOVaY = FOVaX
    else:
        FOVaY = float(np.atleast_1d(FOVaY)[0])
    fpath = getattr(options, 'attenuation_datafile', '')
    if fpath is None or len(str(fpath)) == 0 or not os.path.isdir(str(fpath)):
        raise ValueError('options.attenuation_datafile must be the folder of the DICOM CT images (.ima or .dcm files)!')
    fpath = str(fpath)
    fnames = sorted(f for f in os.listdir(fpath) if f.lower().endswith(('.ima', '.dcm')) and os.path.isfile(os.path.join(fpath, f)))
    numfids = len(fnames)
    if numfids == 0:
        raise ValueError('No .ima or .dcm files found in the folder ' + fpath)
    info = pydicom.dcmread(os.path.join(fpath, fnames[0]), stop_before_pixels=True)
    H = int(info.Rows)
    W = int(info.Columns)
    if KVP is None:
        if info.get('KVP', None) is None:
            raise ValueError('The DICOM images contain no KVP, input it with the KVP parameter!')
        KVP = float(info.KVP)
    attenImage = np.zeros((H, W, numfids))
    noRescale = False
    for k in range(numfids):
        hu, hasRescale = _readCTSliceHU(os.path.join(fpath, fnames[k]))
        if hu.shape != (H, W):
            raise ValueError('The CT slice ' + fnames[k] + ' has a different size than the first slice! The folder must contain a single CT series.')
        attenImage[:, :, k] = hu
        if not hasRescale and not noRescale:
            noRescale = True
            print('DICOM RescaleSlope/RescaleIntercept not found, using the offset -1000 to convert the CT values to HU (file ' + fnames[k] + ')')
    pixSpace = np.atleast_1d(np.asarray(info.PixelSpacing, dtype=np.float64))
    sliceThick = float(info.SliceThickness)
    Nslices = Nz
    sliceLkmAlku = int(np.floor(((numfids * sliceThick - axialFOV) / 2.) / sliceThick))
    if sliceLkmAlku < 0:
        raise ValueError('The axial length of the CT (%d slices of %g mm) is shorter than options.axial_fov (%g mm)!' % (numfids, sliceThick, axialFOV))
    if numfids - 2 * sliceLkmAlku < 1:
        raise ValueError('No CT slices left inside options.axial_fov!')
    attenImage = attenImage[:, :, sliceLkmAlku:numfids - sliceLkmAlku]
    # Axial scaling to Nz slices
    resized = _imresizeSliceAxis(attenImage, Nslices)
    del attenImage
    # Transaxial resampling to the PET image grid by exact area weighting. The first dimension of the CT image is mapped to Nx
    # pixels covering FOVa_x, the second to Ny pixels covering FOVa_y. Both grids are centered at the image center. The parts
    # of the PET FOV outside the CT image are filled with air (-1000 HU).
    if pixSpace.size == 1:
        pixSpace = np.array([pixSpace[0], pixSpace[0]])
    Mx = _overlapMatrix(H, pixSpace[0], Nx, FOVaX / Nx)
    My = _overlapMatrix(W, pixSpace[1], Ny, FOVaY / Ny)
    coverage = np.minimum(np.maximum(np.outer(Mx.sum(axis=1), My.sum(axis=1)), 0.), 1.)
    airPart = (1. - coverage) * (-1000.)
    reresized = np.tensordot(Mx, resized, axes=(1, 0))
    reresized = np.tensordot(reresized, My, axes=(1, 1))
    del resized
    reresized = np.transpose(reresized, (0, 2, 1)) + airPart[:, :, None]
    # The axial interpolation may undershoot slightly below air
    reresized = np.maximum(reresized, -1000.)
    # Narrow-beam water attenuation (0.096 cm^-1) with scatter correction, effective value (0.093 cm^-1) otherwise
    if getattr(options, 'scatter_correction', False):
        waterMu = _WATER_MU_NARROW
        if getattr(options, 'verbose', 0) > 0:
            print('Scatter correction selected: using narrow-beam water attenuation 0.096 cm^-1 for the CT conversion')
    else:
        waterMu = _WATER_MU_EFFECTIVE
    attenuationFactors = attenuationCTTo511(reresized, KVP, name, waterMu)
    attenuationFactors = attenuationFactors / 10.
    return np.asfortranarray(attenuationFactors)
