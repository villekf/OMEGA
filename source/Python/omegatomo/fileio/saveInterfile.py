# -*- coding: utf-8 -*-
"""
Python port of saveInterfile.m. Saves an image in the Interfile format.
"""

import os
import time

import numpy as np

from ._saveUtil import _get, _first, _num2str, _resolveType, _castImage, _voxelSizeOrNone, _nTimeSteps, _totalTime, _asAtLeast3D, _columnMajorBytes


def _scannerParameters(options):
    """Returns the lines of the scanner parameters block of the header (an empty list if the scanner is not described by the
    options)."""
    if options is None or _get(options, 'CT', False) or _get(options, 'SPECT', False):
        return []
    try:
        if _get(options, 'rings') is not None:
            nRings = _first(_get(options, 'rings'))
            perRing = _first(_get(options, 'det_w_pseudo'))
            diameter = _first(_get(options, 'diameter'))
            axialBlocks = _first(_get(options, 'linear_multip'))
            crystAxial = _first(_get(options, 'cryst_per_block_axial'))
            cryst = _first(_get(options, 'cryst_per_block'))
            if crystAxial is None:
                crystAxial = cryst
            angles = _first(_get(options, 'Nang'))
            blocks = _first(_get(options, 'blocks_per_ring'))
        elif _get(options, ('nRings', 'machine_rings')) is not None:
            nRings = _first(_get(options, ('nRings', 'machine_rings')))
            perRing = _first(_get(options, ('detectorsPerRing', 'detectors_per_ring')))
            diameter = _first(_get(options, ('boreDiameter', 'bore_diameter')))
            axialBlocks = _first(_get(options, ('axialBlocks', 'axial_blocks')))
            cryst = _first(_get(options, ('crystalsPerBlock', 'crystals_per_block')))
            crystAxial = cryst
            angles = _first(_get(options, ('nSinogramAngles', 'Nang')))
            blocks = _first(_get(options, ('blocksPerRing', 'blocks_per_ring')))
        else:
            return []
        values = [nRings, perRing, diameter, axialBlocks, crystAxial, cryst, angles, blocks]
        if any(v is None for v in values) or nRings <= 0 or perRing <= 0:
            return []
        doi = _first(_get(options, 'DOI'))
        name = _get(options, ('scannerName', 'machine_name'), '')
    except (TypeError, ValueError, IndexError, KeyError):
        return []
    lines = ['Scanner parameters := ']
    lines.append('Scanner type := ' + str(name))
    lines.append('Number of rings := ' + _num2str(nRings))
    lines.append('Number of detectors per ring := ' + _num2str(perRing))
    lines.append('Inner ring diameter (cm) := ' + _num2str(diameter / 10))
    if doi is not None:
        lines.append('Average depth of interaction (cm) := ' + _num2str(doi / 10))
    else:
        lines.append('Average depth of interaction (cm) := 0')
    lines.append('View offset (degrees) := 0')
    lines.append('Maximum number of non-arc-corrected bins := ' + _num2str(angles))
    lines.append('Default number of arc-corrected bins := ' + _num2str(angles))
    lines.append('Number of blocks per bucket in transaxial direction := ' + _num2str(blocks))
    lines.append('Number of blocks per bucket in axial direction := ' + _num2str(axialBlocks))
    lines.append('Number of crystals per block in axial direction := ' + _num2str(crystAxial))
    lines.append('Number of crystals per block in transaxial direction := ' + _num2str(cryst))
    lines.append('Number of detector layers := 1')
    lines.append('end scanner parameters := ')
    return lines


def _isOne(value):
    value = np.asarray(value)
    return value.size > 0 and bool(np.all(value == 1))


def _formatExtreme(value):
    s = '%e' % float(value)
    return s if s.startswith('-') else '+' + s


def _interfileHeader(dataFile, shape, maxmin, reko, dt, options, voxel, sizeMode):
    """Builds the lines of the Interfile header (same keys and values as writeInterfileHeader in saveInterfile.m)."""
    L = []
    L.append('!INTERFILE := ')
    L.append('%comment := ')
    L.append('!originating system := OMEGA')
    L.append('%SMS-MI header name space := ')
    L.append('%SMS-MI version number := ')
    L.append('')
    L.append('!GENERAL DATA := ')
    L.append('%sinogram header file := ')
    L.append('%sinogram data file := ')
    L.append('!name of data file := ' + dataFile)
    L.append('')
    L.append('!GENERAL IMAGE DATA := ')
    nT = _nTimeSteps(options)
    if nT is not None and nT > 1:
        L.append('!type of data := Dynamic')
        L.append('!total number of images := ' + str(nT * shape[2]))
    else:
        L.append('!type of data := Static')
        L.append('!total number of images := ' + str(shape[2]))
    L.append('!number of slices := ' + str(shape[2]))
    now = time.localtime()
    L.append('%study date (yyyy:mm:dd) := ' + time.strftime('%Y:%m:%d', now))
    L.append('%study time (hh:mm:ss GMT+00:00) := ' + time.strftime('%H:%M:%S', now))
    L.append('isotope name := ')
    L.append('isotope gamma halflife (sec) := ')
    L.append('isotope branching factor := ')
    L.append('radiopharmaceutical := ')
    L.append('%tracer injection date (yyyy:mm:dd) := ')
    L.append('%tracer injection time (hh:mm:ss GMT+00:00) := ')
    L.append('tracer activity at time of injection (Bq) := ')
    L.append('relative time of tracer injection (sec) := ')
    L.append('injected volume (ml) := ')
    L.append('!imagedata byte order := LITTLEENDIAN')
    L.append('%patient orientation := ')
    L.append('%image orientation := ')
    L.append('!PET data type := image')
    nBytes = str(dt.itemsize)
    if dt.kind == 'f':
        L.append('!number of bytes per pixel := ' + nBytes)
        L.append('!number format := ' + ('float' if dt.itemsize == 4 else 'double'))
    elif dt.kind == 'i':
        L.append('!number format := signed integer')
        L.append('!number of bytes per pixel := ' + nBytes)
    else:
        L.append('!number format := unsigned integer')
        L.append('!number of bytes per pixel := ' + nBytes)
    is4D = len(shape) > 3 and shape[3] > 1
    L.append('number of dimensions := ' + ('4' if is4D else '3'))
    L.append('matrix axis label[1] := x')
    L.append('matrix axis label[2] := y')
    L.append('matrix axis label[3] := z')
    if is4D:
        L.append('matrix axis label[4] := t')
    L.append('!matrix size [1] := ' + str(shape[0]))
    L.append('!matrix size [2] := ' + str(shape[1]))
    L.append('!matrix size [3] := ' + str(shape[2]))
    if is4D:
        L.append('!matrix size [4] := ' + str(shape[3]))
    if sizeMode == 'size':
        L.append('scaling factor (mm/pixel) [1] := +' + '%e' % voxel[0])
        L.append('scaling factor (mm/pixel) [2] := +' + '%e' % voxel[1])
        L.append('scaling factor (mm/pixel) [3] := +' + '%e' % voxel[2])
        L.append('slice thickness (pixels) := +' + '%e' % voxel[2])
    else:
        L.append('scaling factor (mm/pixel) [1] := 1')
        L.append('scaling factor (mm/pixel) [2] := 1')
        L.append('scaling factor (mm/pixel) [3] := 1')
        L.append('slice thickness (pixels) := 1')
    L.append('horizontal bed translation := stepped')
    L.append('start horizontal bed position (mm) := ')
    L.append('end horizontal bed position (mm) := ')
    L.append('start vertical bed position (mm) := ')
    L.append('%reconstruction diameter (mm) := ')
    L.append('quantification units := ')
    L.append('%scanner quantification factor (Bq*s/ECAT counts) := ')
    L.append('%decay correction := ')
    L.append('%decay correction reference date (yyyy:mm:dd) := ')
    L.append('%decay correction reference time (hh:mm:ss GMT+00:00) := ')
    L.append('slice orientation := ')
    L.append('method of reconstruction := ' + ('' if reko is None else str(reko)))
    L.append('%PSF axial sigma (mm) := ')
    L.append('%PSF axial cutoff := ')
    L.append('%number of radial PSF coefficients := ')
    L.append('%PSF radial left coefficient (bins) [1] := ')
    L.append('%PSF radial left coefficient (1/mm) [2] := ')
    L.append('%PSF radial left coefficient (1/mm^2) [3] := ')
    L.append('%PSF radial left coefficient (1/mm^3) [4] := ')
    L.append('%PSF radial left coefficient (1/mm^4) [5] := ')
    L.append('%PSF radial right coefficient (bins) [1] := ')
    L.append('%PSF radial right coefficient (1/mm) [2] := ')
    L.append('%PSF radial right coefficient (1/mm^2) [3] := ')
    L.append('%PSF radial right coefficient (1/mm^3) [4] := ')
    L.append('%PSF radial right coefficient (1/mm^4) [5] := ')
    L.append('%PSF radial cutoff := ')
    L.append('%gantry offset (mm) [1] := 0')
    L.append('%gantry offset (mm) [2] := 0')
    L.append('%gantry offset (mm) [3] := 0')
    L.append('%gantry offset pitch (degrees) := ')
    L.append('%gantry offset yaw (degrees) := ')
    L.append('%gantry offset roll (degrees) := ')
    niter = _first(_get(options, ('Niter',)))
    L.append('number of iterations := ' + ('' if niter is None else _num2str(niter)))
    subsets = _first(_get(options, ('subsets',)))
    L.append('%number of subsets := ' + ('' if subsets is None else _num2str(subsets)))
    L.append('filter name := ')
    L.append('%xy-filter (mm) := ')
    L.append('%z-filter (mm) := ')
    L.append('%filter order := ')
    L.append('%image zoom := ')
    L.append('%x-offset (mm) := ')
    L.append('%y-offset (mm) := ')
    corrections = ''
    if _isOne(_get(options, ('normalization', 'normalization_correction', 'normalizationCorrection'), 0)):
        corrections += 'normalization,'
    if _isOne(_get(options, ('randoms', 'randoms_correction', 'randomsCorrection'), 0)):
        corrections += 'randoms,'
    if _isOne(_get(options, ('attenuation', 'attenuation_correction', 'attenuationCorrection'), 0)):
        corrections += 'attenuation,'
    if _isOne(_get(options, ('scatter', 'scatter_correction', 'scatterCorrection'), 0)):
        corrections += 'scatter'
    L.append('applied corrections := ' + corrections)
    L.append('method of attenuation correction := image')
    L.append('%CT coverage := ')
    L.append('method of scatter correction := ')
    L.append('%method of random correction := delayed coincidences')
    L.append('%TOF mashing factor := ')
    L.append('number of energy windows := ')
    L.append('%energy window lower level (keV) [1] := ')
    L.append('%energy window upper level (keV) [1] := ')
    L.append('%coincidence window width (ns) := ')
    L.append('')
    L.append('!IMAGE DATA DESCRIPTION := ')
    totalTime = _totalTime(options)
    if nT is not None:
        L.append('!total number of data sets := ' + str(nT))
        if totalTime is not None and np.isfinite(totalTime):
            L.append('!image duration (sec) := ' + _num2str(totalTime / nT))
        else:
            L.append('!image duration (sec) := 0')
    else:
        L.append('!total number of data sets := 1')
        L.append('!image duration (sec) := 0')
    start = _first(_get(options, ('start',)))
    if start is not None:
        L.append('!image relative start time (sec) := ' + _num2str(start))
    else:
        L.append('!image relative start time (sec) := 0')
    L.append('total prompts := 0')
    L.append('%total randoms := 0')
    L.append('%total net trues := 0')
    L.append('%GIM loss fraction := ')
    L.append('%PDR loss fraction := ')
    L.append('%total uncorrected singles rate := ')
    L.append('%image slope := ')
    L.append('%image intercept := ')
    L.append('maximum pixel count := ' + _formatExtreme(maxmin[0]))
    L.append('minimum pixel count := ' + _formatExtreme(maxmin[1]))
    L.append('')
    L.append('%SUPPLEMENTARY ATTRIBUTES := ')
    L.extend(_scannerParameters(options))
    return L


def saveInterfile(filename, img, reko=None, dtype=None, options=None, pixelSize=None, windowsFormat=False):
    """
    Saves the input image in Interfile format (Python port of saveInterfile.m). The image data are saved in column-major order
    (the first dimension, x, varies the fastest) in little endian byte order to [filename].i33 and the header to
    [filename].h33. By default the header uses the UNIX text format (LF).

    Parameters
    ----------
    filename : str
        Name (and optionally path) of the image and header files, without the file type.
    img : ndarray
        The 3D or 4D image (x, y, z[, t]). Vector images are not supported.
    reko : str, optional
        Name of the reconstruction algorithm, saved to the header (naming only).
    dtype : str or numpy dtype, optional
        Data type of the saved image, e.g. 'float32', 'uint16' (MATLAB names such as 'single' are also accepted). By default
        the type of the input image. The image is cast to this type before saving.
    options : options object or dict, optional
        The OMEGA options object (or a dict with the same field names, or the image properties of MATLAB). Used to save e.g.
        the voxel size (from FOVa_x/Nx, FOVa_y/Ny, axial_fov/Nz, i.e. the field of view divided by the size of the image),
        the number of time steps, the total time and the scanner information.
    pixelSize : sequence of float, optional
        Voxel size (mm) in x, y and z if no options are available. Overrides the voxel size from the options. If neither is
        given, the voxel size is not saved (1 mm).
    windowsFormat : bool, optional
        Use the Windows text format (CR LF) in the header. Default False (LF).
    """
    img = _asAtLeast3D(img)
    dt = _resolveType(dtype, img)
    img = _castImage(img, dt)
    voxel = _voxelSizeOrNone(img.shape, options, pixelSize)
    if voxel is not None:
        sizeMode = 'size'
    else:
        sizeMode = 'unit'
    finite = img[np.isfinite(img)] if img.dtype.kind == 'f' else img
    if finite.size == 0:
        maxmin = (0., 0.)
    else:
        maxmin = (finite.max(), finite.min())
    folder, base = os.path.split(filename)
    with open(filename + '.i33', 'wb') as f:
        f.write(_columnMajorBytes(img, dt, '<'))
    lines = _interfileHeader(base + '.i33', img.shape, maxmin, reko, dt, options, voxel, sizeMode)
    eol = '\r\n' if windowsFormat else '\n'
    with open(filename + '.h33', 'wb') as f:
        f.write((eol.join(lines) + eol).encode('ascii', errors='replace'))
