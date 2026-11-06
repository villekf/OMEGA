# -*- coding: utf-8 -*-
"""
Private helpers shared by saveImage, saveInterfile, saveMetaImage and scaleImages.
"""

import numpy as np

# MATLAB class names (and a few aliases) that are accepted as the output data type
_MATLAB_TYPES = {'single': 'float32', 'double': 'float64', 'char': 'int8', 'uchar': 'uint8', 'logical': 'uint8'}

_SUPPORTED_TYPES = ('float32', 'float64', 'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64')


def _get(obj, names, default=None):
    """Returns the first of the attributes/keys `names` that exists in `obj` (an options object or a dict, e.g. the image
    properties), `default` if none exist."""
    if obj is None:
        return default
    if isinstance(names, str):
        names = (names,)
    for name in names:
        if isinstance(obj, dict):
            if name in obj:
                return obj[name]
        elif hasattr(obj, name):
            return getattr(obj, name)
    return default


def _first(x):
    """Returns the first element of x as a Python float (x can be a scalar or an array, e.g. a multi-resolution value).
    Returns None if x is None or empty."""
    if x is None:
        return None
    a = np.asarray(x).ravel()
    if a.size == 0:
        return None
    return float(a[0])


def _num2str(x):
    """Formats a real scalar in the same way as MATLAB's num2str (integer values without decimals, other values with at least
    five significant digits)."""
    x = float(x)
    if np.isnan(x):
        return 'NaN'
    if np.isinf(x):
        return 'Inf' if x > 0 else '-Inf'
    if x == np.floor(x):
        return '%d' % x
    digits = min(max(int(np.floor(np.log10(abs(x)))) + 5, 5), 16)
    return '%.*g' % (digits, x)


def _resolveType(dtype, img):
    """Returns the numpy dtype in which the image is saved. dtype can be a numpy dtype, a type name ('float32', 'uint16', ...)
    or a MATLAB class name ('single', 'double', ...). By default the type of the image is used."""
    if dtype is None:
        dt = np.dtype(img.dtype)
    else:
        if isinstance(dtype, str) and dtype in _MATLAB_TYPES:
            dtype = _MATLAB_TYPES[dtype]
        dt = np.dtype(dtype)
    if dt.kind == 'b':
        dt = np.dtype('uint8')
    elif dt == np.dtype('float16'):
        dt = np.dtype('float32')
    dt = dt.newbyteorder('=')
    if dt.name not in _SUPPORTED_TYPES:
        raise ValueError('Unsupported data type ' + str(dt) + '. Supported types are ' + ', '.join(_SUPPORTED_TYPES))
    return dt


def _castImage(img, dt):
    """Casts the image to the data type dt in the same way as MATLAB's cast: floating point values are rounded (halves away from
    zero) and values outside the range of an integer type saturate."""
    img = np.asarray(img)
    if img.dtype.kind == 'b':
        img = img.astype(np.uint8)
    if img.dtype == dt:
        return img
    if dt.kind in 'iu':
        info = np.iinfo(dt)
        if img.dtype.kind == 'f':
            a = np.where(np.isnan(img), 0., img).astype(np.float64)
            a = np.sign(a) * np.floor(np.abs(a) + 0.5)
            a = np.clip(a, float(info.min), float(info.max))
            return a.astype(dt)
        srcInfo = np.iinfo(img.dtype)
        lo = max(info.min, srcInfo.min)
        hi = min(info.max, srcInfo.max)
        return np.clip(img, lo, hi).astype(dt)
    return img.astype(dt)


def _voxelSizeOrNone(shape, options=None, pixelSize=None):
    """Returns the voxel size (dx, dy, dz) in mm. The explicit pixelSize has the highest priority, after that the field of view
    of the options (or the image properties) divided by the size of the image. None if neither is available."""
    if pixelSize is not None:
        p = np.asarray(pixelSize, dtype=np.float64).ravel()
        if p.size < 2:
            raise ValueError('pixelSize needs to have at least two elements (x and y)')
        if p.size == 2:
            p = np.append(p, 1.)
        return float(p[0]), float(p[1]), float(p[2])
    if options is not None:
        fx = _first(_get(options, ('FOVa_x', 'FOVx', 'FOV_x')))
        fy = _first(_get(options, ('FOVa_y', 'FOVy', 'FOV_y')))
        fz = _first(_get(options, ('axial_fov', 'axialFOV', 'axial_FOV')))
        if fx is not None and fy is not None and fz is not None and fx > 0 and fy > 0 and fz > 0:
            return fx / shape[0], fy / shape[1], fz / shape[2]
    return None


def _voxelSize(shape, options=None, pixelSize=None):
    """Same as _voxelSizeOrNone, but 1 mm voxels if the voxel size is not available."""
    voxel = _voxelSizeOrNone(shape, options, pixelSize)
    return (1., 1., 1.) if voxel is None else voxel


def _nTimeSteps(options):
    """The number of time steps of the options (partitions) or of the image properties (nTimeSteps), None if not available."""
    n = _get(options, ('nTimeSteps', 'n_time_steps', 'partitions'))
    if n is None:
        return None
    n = np.asarray(n)
    if n.ndim > 0:
        return max(int(n.size), 1)
    return int(n)


def _totalTime(options):
    """Total time (s) of the measurement, None if not available."""
    return _first(_get(options, ('totalTime', 'total_time', 'tot_time')))


def _lastIteration(img):
    """Returns the image that is saved from a 3D, 4D or 5D input: 4D input is X x Y x Z x iterations and the last iteration is
    used, 5D input is X x Y x Z x iterations x time and the result is X x Y x Z x time."""
    img = np.asarray(img)
    if img.ndim > 5:
        raise ValueError('The input image can have at most 5 dimensions')
    if img.ndim < 2:
        raise ValueError('The size of the input image has to be at least XxYxZ')
    if img.ndim == 2:
        img = img[:, :, np.newaxis]
    elif img.ndim >= 4:
        img = img[:, :, :, -1, ...]
    while img.ndim > 3 and img.shape[-1] == 1:
        img = img[..., 0]
    return img


def _asAtLeast3D(img):
    img = np.asarray(img)
    if img.ndim < 2 or img.ndim > 4:
        raise ValueError('The input image has to be 2D, 3D or 4D')
    if img.ndim == 2:
        img = img[:, :, np.newaxis]
    return img


def _columnMajorBytes(img, dt, byteOrder='<'):
    """The image data as raw bytes in column-major (MATLAB) order, converted to the data type dt and byte order ('<' or '>')."""
    return np.asarray(img, dtype=dt.newbyteorder(byteOrder)).tobytes(order='F')
