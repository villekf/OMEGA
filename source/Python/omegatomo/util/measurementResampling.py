# -*- coding: utf-8 -*-
"""
Automatic measurement resampling for the SPECT rotation-based (type 6) projector.

The type-6 projector works on a projection grid derived from the reconstruction volume. This module ports
the SPECT measurement-window preparation from MATLAB's projectorClass.m, including center-matched
resampling and the full-detector window needed when the object or extended FOV is off-centre.
"""
import numpy as np

from .matlabRound import matlabRound


def _cubicKernel(x):
    """The cubic convolution kernel (a = -0.5) used by MATLAB's imresize for the 'bicubic' method."""
    absx = np.abs(x)
    absx2 = absx ** 2
    absx3 = absx ** 3
    return (1.5 * absx3 - 2.5 * absx2 + 1.) * (absx <= 1.) + (-0.5 * absx3 + 2.5 * absx2 - 4. * absx + 2.) * ((absx > 1.) & (absx <= 2.))


def _imresizeWeights(inLength, outLength):
    """
    Dense (outLength x inLength) matrix of the 1-D resampling weights of MATLAB's imresize with the
    default bicubic kernel and antialiasing. Follows imresize's contributions(): the kernel is widened
    by 1 / scale when shrinking, each row is normalized to sum to one and out-of-range input indices
    are mirrored (symmetric padding).
    """
    scale = float(outLength) / float(inLength)
    kernelWidth = 4.
    if scale < 1.:
        h = lambda x: scale * _cubicKernel(scale * x)
        kernelWidth = kernelWidth / scale
    else:
        h = _cubicKernel
    x = np.arange(1, outLength + 1, dtype=np.float64)
    # Input-space coordinates; 0.5 in the output space maps to 0.5 in the input space
    u = x / scale + 0.5 * (1. - 1. / scale)
    left = np.floor(u - kernelWidth / 2.)
    P = int(np.ceil(kernelWidth)) + 2
    indices = left[:, None] + np.arange(P, dtype=np.float64)[None, :]
    weights = h(u[:, None] - indices)
    weights = weights / np.sum(weights, axis=1, keepdims=True)
    # Mirror the out-of-range indices (1-based indexing as in MATLAB)
    aux = np.concatenate((np.arange(1, inLength + 1), np.arange(inLength, 0, -1)))
    indices = aux[np.mod(indices - 1., 2. * inLength).astype(np.int64)] - 1
    matrix = np.zeros((outLength, inLength), dtype=np.float64)
    rows = np.broadcast_to(np.arange(outLength)[:, None], indices.shape)
    np.add.at(matrix, (rows, indices), weights)
    return matrix


def _imresizeAxis(data, outLength, axis):
    """Resamples one axis of a float64 array with the imresize bicubic weights."""
    if data.shape[axis] == outLength:
        return data
    matrix = _imresizeWeights(data.shape[axis], outLength)
    return np.moveaxis(np.tensordot(matrix, data, axes=([1], [axis])), 0, axis)


def _imresizeBicubic(data, outRows, outCols):
    """
    Port of MATLAB's imresize(data, [outRows, outCols]) (default bicubic method with antialiasing) for the
    first two dimensions of an N-D array; all further dimensions are resampled independently of each other.
    The computation is done in double precision. Floating point inputs keep their data type, other
    inputs (integers, logical) are returned as float32.
    """
    data = np.asarray(data)
    outDtype = data.dtype if data.dtype.kind == 'f' else np.dtype(np.float32)
    result = data.astype(np.float64)
    result = _imresizeAxis(result, int(outRows), 0)
    result = _imresizeAxis(result, int(outCols), 1)
    return result.astype(outDtype)


def _resizeCentered(data, length, axis):
    """
    Pads with zeros or crops one axis to the given length on both sides. As in MATLAB's
    resize(..., Side="both"), floor(|difference| / 2) elements are added/removed in the beginning and the
    rest at the end.
    """
    current = data.shape[axis]
    if current == length:
        return data
    index = [slice(None)] * data.ndim
    if current > length:
        start = (current - length) // 2
        index[axis] = slice(start, start + length)
        return data[tuple(index)]
    shape = list(data.shape)
    shape[axis] = length
    padded = np.zeros(shape, dtype=data.dtype)
    start = (length - current) // 2
    index[axis] = slice(start, start + current)
    padded[tuple(index)] = data
    return padded


def _resampleDetectorAxis(data, axis, nOut, ratio):
    """Bicubic resampling used by MATLAB's resampleDetectorAxis, with zero outside the detector."""
    data = np.asarray(data)
    nIn, nOut, ratio = int(data.shape[axis]), int(nOut), float(ratio)
    if nIn < 1 or nOut < 1 or ratio <= 0. or not np.isfinite(ratio):
        raise ValueError('Type-6 detector resampling needs positive finite dimensions and pitch ratios.')
    if nIn == nOut and ratio == 1.:
        return data
    if ratio > 1.:
        kernel = lambda x: (1. / ratio) * _cubicKernel(x / ratio)
        kernelWidth = 4. * ratio
    else:
        kernel = _cubicKernel
        kernelWidth = 4.
    # Match resampleDetectorAxis.m: normalize weights before dropping samples
    # outside the detector, so those samples contribute zeros rather than a
    # mirrored edge value.
    j = np.arange(1, nOut + 1, dtype=np.float64)
    u = (j - (nOut + 1.) / 2.) * ratio + (nIn + 1.) / 2.
    left = np.floor(u - kernelWidth / 2.)
    count = int(np.ceil(kernelWidth)) + 2
    indices = left[:, None] + np.arange(count, dtype=np.float64)[None, :]
    weights = kernel(u[:, None] - indices)
    weights /= np.sum(weights, axis=1, keepdims=True)
    valid = (indices >= 1.) & (indices <= nIn)
    matrix = np.zeros((nOut, nIn), dtype=np.float64)
    rows, cols = np.nonzero(valid)
    np.add.at(matrix, (rows, indices[valid].astype(np.int64) - 1), weights[valid])
    result = np.moveaxis(np.tensordot(matrix, data.astype(np.float64, copy=False), axes=([1], [axis])), 0, axis)
    return result.astype(data.dtype if data.dtype.kind == 'f' else np.float32, copy=False)


def _nearestIntegerRatio(ratio):
    """Return MATLAB's rounded physical size and whether the input ratio was exactly integral."""
    nearest = float(matlabRound(ratio))
    exact = float(ratio) == nearest
    if not exact and abs(float(ratio) - nearest) > 1e-9 * max(1., abs(float(ratio))):
        nearest = float(ratio)
    return nearest, exact


def _resampleArray(data, endRows, endCols, outRows, outCols):
    """
    Resample one type-6 measurement array to an output lattice of outRows x outCols.
    Integer detector/FOV ratios with symmetric padding/cropping use the historic
    centered path; other ratios use the centered detector-axis resampler.
    """
    result = np.asarray(data)
    rows, exactRows = _nearestIntegerRatio(endRows)
    cols, exactCols = _nearestIntegerRatio(endCols)
    if exactRows and (result.shape[0] - int(rows)) % 2 == 0:
        result = _resizeCentered(result, int(rows), 0)
    else:
        result = _resampleDetectorAxis(result, 0, outRows, rows / float(outRows))
    if exactCols and (result.shape[1] - int(cols)) % 2 == 0:
        result = _resizeCentered(result, int(cols), 1)
    else:
        result = _resampleDetectorAxis(result, 1, outCols, cols / float(outCols))
    if result.shape[0] != outRows or result.shape[1] != outCols:
        result = _imresizeBicubic(result, outRows, outCols)
    elif result.dtype.kind != 'f':
        result = result.astype(np.float32)
    if np.may_share_memory(result, data):
        result = result.copy()
    return np.asfortranarray(result)


def _firstArray(value):
    """The first non-empty ndarray found in a (nested) list/tuple, or None."""
    if isinstance(value, np.ndarray):
        return value if value.size > 0 else None
    if isinstance(value, (list, tuple)):
        for v in value:
            found = _firstArray(v)
            if found is not None:
                return found
    return None


def _resampleType6Measurements(options):
    """
    Resamples the SPECT measurements to the grid of the type-6 (rotation-based) projector, as the MATLAB
    projectorClass.m does. Centered data are mapped to the fine-volume Nx x Nz grid. For an off-centre
    object/eFOV shift (or explicit cropProjections), the full detector window is resampled onto a canvas
    centered on the scanner axis, then cropped to that window. This preserves the detector's absolute
    location while aligning the data to the type-6 voxel lattice.

    The following are resampled when they have the measurement grid (first two dimensions equal to that of
    options.SinM): options.SinM, options.ScatterC (a single array or the DEW/TEW energy window list, also
    per-timestep lists) and options.SinDelayed. Dynamic data (a list of frames or a 4-D array) is supported.
    Integer data is returned as float32 (resampled counts are not integers).

    This must run before parseInputs(). The built-in reconstruction may call it after addProjector() when it
    loads measurements from fpath, so the detector-size dependent fields are refreshed here too. The
    measurement arrays are rebound, never written in place (the arrays can be shared with the caller).

    Parameters
    ----------
    options : class object
        The (private copy of the) OMEGA options. Modified: SinM, ScatterC, SinDelayed, nRowsD, nColsD.
    """
    if getattr(options, '_type6MeasurementsResampled', False):
        return
    sinM = options.SinM
    reference = _firstArray(sinM)
    if reference is None:
        # Nothing to resample (e.g. the measurements are loaded later).
        return
    if getattr(options, 'listmode', False):
        return
    if (isinstance(sinM, np.ndarray) and sinM.ndim == 1 and
            np.asarray(getattr(options, 'x', np.empty(0))).size == sinM.size * 6):
        # Listmode SPECT uses SinM as event weights and x as one source/detector
        # coordinate tuple per event, not as a flattened detector image.
        return

    def scalar(name):
        values = np.asarray(getattr(options, name)).reshape(-1)
        if values.size == 0:
            raise ValueError('options.' + name + ' is empty, it is required for the type-6 measurement resampling!')
        return float(values[0])

    sinM = options.SinM
    reference = _firstArray(sinM)
    if reference is None:
        # Nothing to resample (e.g. the measurements are loaded later)
        return
    dPitchX = scalar('dPitchX')
    dPitchY = scalar('dPitchY')
    # The same pixel size defaults as in addProjector
    dPitch = float(np.asarray(getattr(options, 'dPitch', 0.)).reshape(-1)[0]) if np.size(getattr(options, 'dPitch', 0.)) > 0 else 0.
    if dPitchX == 0. and dPitch > 0.:
        dPitchX = dPitch
    if dPitchY == 0. and dPitch > 0.:
        dPitchY = dPitch
    if dPitchY == 0. and dPitchX > 0.:
        dPitchY = dPitchX
    FOVx = scalar('FOVa_x')
    axialFOV = scalar('axial_fov')
    if min(dPitchX, dPitchY, FOVx, axialFOV) <= 0.:
        raise ValueError('The type-6 measurement resampling requires positive options.dPitchX (or dPitch), options.dPitchY, '
                         'options.FOVa_x and options.axial_fov!')
    nx = int(scalar('Nx'))
    nz = int(scalar('Nz'))
    if reference.ndim >= 2:
        measRows, measCols = reference.shape[0], reference.shape[1]
    else:
        measRows, measCols = int(options.nRowsD), int(options.nColsD)
    if measRows < 1 or measCols < 1:
        raise ValueError('options.SinM must be at least two-dimensional or options.nRowsD and options.nColsD must be set!')

    offsets = np.asarray([getattr(options, 'oOffsetX', 0.), getattr(options, 'oOffsetY', 0.),
                          getattr(options, 'oOffsetZ', 0.)], dtype=np.float64)
    shifts = np.asarray(getattr(options, 'eFOVShift', (0., 0., 0.)), dtype=np.float64).reshape(-1)
    offCentre = bool(np.any(offsets != 0.) or np.any(shifts != 0.))
    cropSetting = getattr(options, 'cropProjections', False)
    if cropSetting is None:
        tolX = 0.5 * FOVx / nx
        tolZ = 0.5 * axialFOV / nz
        crop = (measRows * dPitchX < FOVx - tolX or
                measCols * dPitchY < axialFOV - tolZ)
    else:
        crop = bool(cropSetting)
    if offCentre:
        # Padding to the FOV symmetrically would place an offset object on the
        # detector center. Keep the full measured detector window instead.
        crop = True

    if crop:
        targetRows = measRows * dPitchX / (FOVx / nx)
        targetCols = measCols * dPitchY / (axialFOV / nz)
        newRows = int(matlabRound(targetRows))
        newCols = int(matlabRound(targetCols))
        if (newRows - nx) % 2:
            newRows += 1 if targetRows >= newRows else -1
        if (newCols - nz) % 2:
            newCols += 1 if targetCols >= newCols else -1
        if not offCentre:
            newRows = min(newRows, nx)
            newCols = min(newCols, nz)
        newRows, newCols = max(1, newRows), max(1, newCols)
        canvasRows, canvasCols = max(nx, newRows), max(nz, newCols)
        rowStart, colStart = (canvasRows - newRows) // 2, (canvasCols - newCols) // 2
        if canvasRows == nx:
            endRows = FOVx / dPitchX
        else:
            endRows = float(matlabRound(canvasRows * (FOVx / nx) / dPitchX))
        if canvasCols == nz:
            endCols = axialFOV / dPitchY
        else:
            endCols = float(matlabRound(canvasCols * (axialFOV / nz) / dPitchY))
        outputRows, outputCols = canvasRows, canvasCols
    else:
        newRows, newCols = nx, nz
        rowStart = colStart = 0
        endRows, endCols = FOVx / dPitchX, axialFOV / dPitchY
        outputRows, outputCols = nx, nz

    def matches(array):
        return isinstance(array, np.ndarray) and array.ndim >= 2 and array.shape[0] == measRows and array.shape[1] == measCols

    def resample(value, isSinM):
        if isinstance(value, list):
            return [resample(v, isSinM) for v in value]
        if isinstance(value, tuple):
            return tuple(resample(v, isSinM) for v in value)
        if isinstance(value, np.ndarray) and value.size > 0:
            if isSinM and value.ndim == 1:
                if value.size % (measRows * measCols) != 0:
                    raise ValueError('The number of elements in options.SinM is not a multiple of nRowsD * nColsD!')
                value = np.reshape(value, (measRows, measCols, -1), order='F')
            if matches(value):
                value = _resampleArray(value, endRows, endCols, outputRows, outputCols)
                if crop:
                    value = value[rowStart:rowStart + newRows, colStart:colStart + newCols, ...]
                return np.asfortranarray(value)
        return value

    if getattr(options, 'verbose', 1) > 0:
        suffix = ' (off-centre detector window)' if offCentre else ''
        print('Type-6 SPECT: resampling measurements from %d x %d to %d x %d pixels%s' %
              (measRows, measCols, newRows, newCols, suffix))
    options.SinM = resample(options.SinM, True)
    options.ScatterC = resample(options.ScatterC, False)
    if hasattr(options, 'SinDelayed'):
        options.SinDelayed = resample(options.SinDelayed, False)
    options.cropProjections = crop
    options.nRowsD = newRows
    options.nColsD = newCols
    options._type6MeasurementsResampled = True
    # Keep any detector dimensions that addProjector() initialized in sync
    # with the newly resampled arrays before native transfer or parseInputs().
    nProjections = int(getattr(options, 'nProjections', 0))
    for name, value in (('measurement_nRowsD', newRows), ('measurement_nColsD', newCols),
                        ('Ndist', newRows), ('Nang', newCols),
                        ('det_per_ring', newRows * nProjections),
                        ('det_w_pseudo', newRows * nProjections),
                        ('cryst_per_block', newRows * newCols)):
        if hasattr(options, name):
            setattr(options, name, value)
