# -*- coding: utf-8 -*-
"""
Automatic measurement resampling for the SPECT rotation-based (type 6) projector.

The type-6 projector works on a projection grid of Nx x Nz pixels (rows x columns), i.e. on the
transaxial/axial grid of the reconstruction volume. This module ports the sinogram resizing and
resampling block of MATLAB's projectorClass.m: every projection frame is first padded with zeros or
cropped (both sides, centered) to the number of detector pixels that cover the image FOV
(FOVa_x / dPitchX rows and axial_fov / dPitchY columns), then resampled with imresize (bicubic,
antialiasing) to Nx x Nz pixels. options.nRowsD and options.nColsD are set to Nx and Nz.
"""
import numpy as np


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


def _resampleArray(data, endRows, endCols, nx, nz):
    """
    The type-6 measurement resampling of one array (rows x columns x ...): pad/crop to endRows x endCols,
    then bicubic imresize to nx x nz. Returns a new array (never a view of the input).
    """
    result = _resizeCentered(data, endRows, 0)
    result = _resizeCentered(result, endCols, 1)
    if result.shape[0] != nx or result.shape[1] != nz or result.dtype.kind != 'f':
        result = _imresizeBicubic(result, nx, nz)
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
    projectorClass.m does. Each projection frame (nRowsD x nColsD) is padded with zeros/cropped (centered) to
    round(FOVa_x / dPitchX) x round(axial_fov / dPitchY) pixels and then resampled (imresize, bicubic with
    antialiasing) to Nx x Nz pixels. options.nRowsD and options.nColsD are set to Nx and Nz.

    The following are resampled when they have the measurement grid (first two dimensions equal to that of
    options.SinM): options.SinM, options.ScatterC (a single array or the DEW/TEW energy window list, also
    per-timestep lists) and options.SinDelayed. Dynamic data (a list of frames or a 4-D array) is supported.
    Integer data is returned as float32 (resampled counts are not integers).

    This must be called before options.addProjector(), since that derives the detector-size dependent values
    from nRowsD and nColsD, and before options.flipImageZ is applied to the projections, as in MATLAB.
    The measurement arrays are rebound, never written in place (the arrays can be shared with the caller).

    Parameters
    ----------
    options : class object
        The (private copy of the) OMEGA options. Modified: SinM, ScatterC, SinDelayed, nRowsD, nColsD.
    """
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
    # Desired number of detector pixels (rows, columns) that cover the image FOV
    endRows = max(1, int(round(FOVx / dPitchX)))
    endCols = max(1, int(round(axialFOV / dPitchY)))
    if reference.ndim >= 2:
        measRows, measCols = reference.shape[0], reference.shape[1]
    else:
        measRows, measCols = int(options.nRowsD), int(options.nColsD)
    if measRows < 1 or measCols < 1:
        raise ValueError('options.SinM must be at least two-dimensional or options.nRowsD and options.nColsD must be set!')

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
                return _resampleArray(value, endRows, endCols, nx, nz)
        return value

    if getattr(options, 'verbose', 1) > 0:
        print('Type-6 SPECT: resampling the measurements from %d x %d to %d x %d pixels (FOV-matched grid %d x %d)' % (measRows, measCols, nx, nz, endRows, endCols))
    options.SinM = resample(options.SinM, True)
    options.ScatterC = resample(options.ScatterC, False)
    if hasattr(options, 'SinDelayed'):
        options.SinDelayed = resample(options.SinDelayed, False)
    options.nRowsD = nx
    options.nColsD = nz
