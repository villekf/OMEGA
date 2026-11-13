# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 13:50:49 2024

Copyright (C) 2024-2026 Ville-Veikko Wettenhovi

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
"""
import copy
import ctypes
import numpy as np
import warnings
from ..projector.proj import _copyValue

def _spect_listmode_sensitivity_weights(options, n_views):
    supplied = np.asarray(getattr(options, 'sensitivityViewWeights', np.empty(0)))
    if supplied.size:
        weights = np.asarray(supplied, dtype=np.float32)
        if weights.shape != (n_views, int(options.Nt)) or not np.all(np.isfinite(weights)) or np.any(weights < 0):
            raise ValueError('sensitivityViewWeights must be a finite, non-negative nViews-by-Nt array.')
        return np.asfortranarray(weights)
    if int(options.Nt) == 1:
        return np.ones((n_views, 1), dtype=np.float32, order='F')

    frame_index = np.asarray(getattr(options, 'temporalBinIndex', np.empty(0))).reshape(-1)
    if frame_index.size == n_views and np.all(np.isfinite(frame_index)) and np.all(frame_index == np.floor(frame_index)) and np.all((frame_index >= 0) & (frame_index < int(options.Nt))):
        weights = np.zeros((n_views, int(options.Nt)), dtype=np.float32, order='F')
        weights[np.arange(n_views), frame_index.astype(np.int64)] = 1
        return weights

    capture_start = np.asarray(getattr(options, 'measurementStartMs', np.empty(0)), dtype=np.float64).reshape(-1)
    capture_end = np.asarray(getattr(options, 'measurementEndMs', np.empty(0)), dtype=np.float64).reshape(-1)
    frame_start = np.asarray(getattr(options, 'dynamicPartitionStartMs', np.empty(0)), dtype=np.float64).reshape(-1)
    frame_end = np.asarray(getattr(options, 'dynamicPartitionEndMs', np.empty(0)), dtype=np.float64).reshape(-1)
    if capture_start.size == n_views and capture_end.size == n_views and frame_start.size == int(options.Nt) and frame_end.size == int(options.Nt):
        duration = capture_end - capture_start
        if np.all(np.isfinite(capture_start)) and np.all(np.isfinite(capture_end)) and np.all(duration >= 0) and np.all(np.isfinite(frame_start)) and np.all(np.isfinite(frame_end)):
            overlap = np.maximum(0, np.minimum(capture_end[:, None], frame_end[None, :]) - np.maximum(capture_start[:, None], frame_start[None, :]))
            weights = np.zeros((n_views, int(options.Nt)), dtype=np.float32, order='F')
            positive = duration > 0
            weights[positive, :] = (overlap[positive, :] / duration[positive, None]).astype(np.float32)
            zero = ~positive
            if np.any(zero):
                if frame_index.size != n_views:
                    raise ValueError('Zero-duration SPECT views require temporalBinIndex for dynamic sensitivity.')
                for timestep in range(int(options.Nt)):
                    weights[zero & (frame_index == timestep), timestep] = 1
            return weights
    raise ValueError('Dynamic listmode SPECT sensitivity requires view timing, temporalBinIndex, or explicit sensitivityViewWeights.')

def _saved_image_count(options):
    if options.save_iter:
        return int(options.Niter) + 1
    if options.saveNIter.size > 0:
        return int(options.saveNIter.size) + 1
    return 1

def _reshape_image_output(output, spatialShape, savedImageCount, timeFrameCount):
    # The C++ writer stores [voxel, timestep, save slot] in Fortran order,
    # matching the MATLAB MEX dimensions [Nx, Ny, Nz, Nt, saves].
    shape = (*map(int, spatialShape), int(timeFrameCount), int(savedImageCount))
    expectedSize = int(np.prod(shape, dtype=np.int64))
    if output.size != expectedSize:
        raise ValueError(f"Unexpected reconstruction output size: got {output.size} values, expected {expectedSize}")
    output = output.reshape(shape, order='F')
    if timeFrameCount == 1:
        output = output[..., 0, :]
    elif savedImageCount == 1:
        output = output[..., 0]
    return output

def _reshape_multiresolution_output(output, options):
    nVolumes = int(options.nMultiVolumes) + 1
    volSizes = (options.Nx[:nVolumes].astype(np.uint64) * options.Ny[:nVolumes].astype(np.uint64) * options.Nz[:nVolumes].astype(np.uint64)).astype(np.int64)
    savedImageCount = _saved_image_count(options)
    volumes = []
    offset = 0
    for ii, volSize in enumerate(volSizes):
        outputSize = int(volSize) * savedImageCount * int(options.Nt)
        nextOffset = offset + outputSize
        volumes.append(_reshape_image_output(
            output[offset:nextOffset],
            (options.Nx[ii], options.Ny[ii], options.Nz[ii]),
            savedImageCount,
            options.Nt,
        ))
        offset = nextOffset

    if offset != output.size:
        raise ValueError(f"Unexpected multi-resolution output size: consumed {offset} values from {output.size}")
    return volumes

# ctypes scalar type -> NumPy dtype, used by _as_ptr() below to validate/coerce
# the arrays pointed to by POINTER(...) struct fields.
_CTYPE_TO_NUMPY_DTYPE = {
    ctypes.c_uint8: np.uint8,
    ctypes.c_int8: np.int8,
    ctypes.c_uint16: np.uint16,
    ctypes.c_int16: np.int16,
    ctypes.c_uint32: np.uint32,
    ctypes.c_int32: np.int32,
    ctypes.c_uint64: np.uint64,
    ctypes.c_int64: np.int64,
    ctypes.c_float: np.float32,
    ctypes.c_double: np.float64,
    ctypes.c_bool: np.bool_,
}

# Scalar struct fields whose C-struct name differs from the options attribute
# name they are sourced from.
_SCALAR_NAME_OVERRIDES = {
    'POCS': 'ASD_POCS',
    'T': '_TVweight',
}

# Pointer struct fields whose C-struct name differs from the options
# attribute name they are sourced from.
_POINTER_NAME_OVERRIDES = {
    'atten': 'vaimennus',
    'norm': 'normalization',
    'pituus': 'nMeas',
    'gaussPSF': 'gaussK',
    'saveNiter': 'saveNIter',
    'randoms': 'SinDelayed',
    'offsetVal': 'OffsetLimit',
    'kerroin4': 'kerroin',
    'filter': 'filter0',
    'TV_ref': 'TV_referenceImage',
    'NLM_ref': 'NLM_referenceImage',
    'RDP_ref': 'RDP_referenceImage',
    'trIndices': 'trIndex',
    'axIndices': 'axIndex',
    'detectorVector': 'DetectorVector',
}

# Pointer struct fields set through bespoke logic in transferData() (the
# type-6 and ODRT PSF lookup buffers), rather than through
# _POINTER_NAME_OVERRIDES and _as_ptr().
_POINTER_SPECIAL_FIELDS = {'blurPlanes', 'blurPlanes2', 'gFilter', 'gFSize', 'gFilterODRT',
                           'blurPlanesF', 'blurPlanes2F', 'CDRFCoeffs6', 'CDRFOnTheFly', 'cdrfModel6', 'bpRot6', 'radiusPerProj', 'CORtoDetectorSurface'}


# Large measurement-sized data fields that _privateOptionsCopy() shares with the
# caller instead of deep-copying (memory): for these only the (nested) list
# containers are copied, ndarrays are shared. Every write to them in the
# reconstruction pipeline therefore has to be out-of-place (rebind, never write
# into the array) or take a private copy first, unless _sharesCallerMemory() says
# that the array no longer shares memory with the caller's arrays.
_SHARED_DATA_FIELDS = ('SinM', 'SinDelayed', 'ScatterC', 'normalization', 'corrVector', 'trIndex', 'axIndex', 'TOFIndices')
# Listmode event coordinate arrays, shared in the same way but only for
# listmode/index-based data (see _optionsAreListmode()). For sinogram and CT data
# they are derived/rescaled in place by the coordinate code and are copied.
_SHARED_LISTMODE_FIELDS = ('x', 'z', 'uV')


def _copyContainers(value):
    """
    Like _copyValue(), but only the list/tuple/dict containers are copied; the
    ndarrays inside them are shared with the original.
    """
    if isinstance(value, list):
        return [_copyContainers(v) for v in value]
    if isinstance(value, tuple):
        items = [_copyContainers(v) for v in value]
        if hasattr(value, '_fields'):
            return type(value)(*items)
        return tuple(items)
    if isinstance(value, dict):
        return {k: _copyContainers(v) for k, v in value.items()}
    return value


def _collectArrays(value, out):
    """Appends all non-empty ndarrays found in (nested) lists/tuples/dicts to `out`."""
    if isinstance(value, np.ndarray):
        if value.size > 0:
            out.append(value)
    elif isinstance(value, (list, tuple)):
        for v in value:
            _collectArrays(v, out)
    elif isinstance(value, dict):
        for v in value.values():
            _collectArrays(v, out)


def _sharesCallerMemory(options, array):
    """
    True if `array` may share memory with an array that _privateOptionsCopy()
    shared with the caller, i.e. if it must not be written in place. Conservative
    (np.may_share_memory only compares the memory bounds). Always False when the
    options were not copied (keepOptionsUnchanged = False) or for non-arrays.
    """
    if not isinstance(array, np.ndarray):
        return False
    return any(np.may_share_memory(array, s) for s in getattr(options, '_callerSharedArrays', ()))


def _optionsAreListmode(options):
    """
    True if the data in `options` is going to be handled as listmode/index-based
    data. options.listmode itself is only derived inside addProjector(), i.e. it
    is not yet set when the options are copied, so the same detection that
    addProjector() uses is evaluated here. CT and SPECT data never count (their
    x/z/uV are coordinates that are rescaled in place). Any failure to decide
    means "not listmode", which just results in the safe full copy.
    """
    if options.CT or options.SPECT:
        return False
    try:
        return bool(options.listmode) or bool(options.useIndexBasedReconstruction) or bool(options._isListModeCandidate())
    except Exception:
        return False


def _privateOptionsCopy(options):
    """
    Returns a private working copy of the options object.

    MATLAB passes the options struct by value, so a reconstruction never changes
    the caller's options. The Python pipeline (addProjector, parseInputs,
    prepassPhase, loadCorrections, linearizeData, the CT coordinate code, ...)
    modifies the options object in place and not idempotently (Nx/Ny/Nz become
    arrays, dL is rescaled, angles are offset, uV is flipped and scaled, the
    measurement data is linearized and reordered by subset, ...). Calling
    reconstructions_main twice with the same object would therefore start the
    second call from the already-modified values. Working on a copy keeps the
    caller's object exactly as the user set it.

    All instance attributes are copied recursively (see _copyValue()), except
    the large measurement/index data in _SHARED_DATA_FIELDS (and the listmode
    coordinates in _SHARED_LISTMODE_FIELDS for listmode data), which are shared
    with the caller. (The mutable class-level defaults, e.g. precondTypeMeas, are
    already instance attributes: projectorClass.__init__ gives every object its
    own copy.) The ctypes parameter struct is created anew. The shared arrays are
    recorded in _callerSharedArrays.
    Not used if options.keepOptionsUnchanged is False.
    """
    work = copy.copy(options)
    work.param = type(options.param)()
    shared = set(_SHARED_DATA_FIELDS)
    if _optionsAreListmode(options):
        shared.update(_SHARED_LISTMODE_FIELDS)
    sharedArrays = []
    for name, value in list(vars(options).items()):
        if name == 'param':
            continue
        if name in shared:
            _collectArrays(value, sharedArrays)
            setattr(work, name, _copyContainers(value))
        else:
            setattr(work, name, _copyValue(value))
    # Used by _sharesCallerMemory(): lets the pipeline write in place to arrays it owns
    work._callerSharedArrays = tuple(sharedArrays)
    return work


def _as_ptr(options, attr, ctype):
    """
    Returns a ctypes POINTER(ctype) to the array stored in options.<attr>.

    The pointed-to array is guaranteed to have the NumPy dtype matching
    `ctype` and to be contiguous, WITHOUT changing an already C- or
    F-contiguous array's memory order: OMEGA routinely hands this function
    Fortran-ordered multi-dimensional arrays (e.g. weighted_weights, the CT
    uV/z coordinate arrays, masks, gaussK), and the C++ side expects that
    same F order, so silently forcing C order here would reorder the memory
    and hand the native code different data (a bug that a 1-D-array-only
    regression check cannot see: a 1-D array is both C- and F-contiguous).

    If options.<attr> already has the right dtype and is C- or F-contiguous,
    it is used unchanged. If only its dtype is wrong, it is converted with
    `astype(dtype, order='K')`, which preserves its existing C/F memory
    order. If it is neither C- nor F-contiguous (a real edge case; today
    this silently produced a garbage pointer), it is copied to Fortran
    order (OMEGA's convention) with a warning naming the attribute. Either
    way, the (possibly new) array is stored back onto options.<attr> so it
    stays alive -- and the pointer stays valid -- for as long as options is
    used.
    """
    dtype = _CTYPE_TO_NUMPY_DTYPE[ctype]
    arr = getattr(options, attr)
    if arr.dtype != dtype or not (arr.flags['C_CONTIGUOUS'] or arr.flags['F_CONTIGUOUS']):
        if arr.flags['C_CONTIGUOUS'] or arr.flags['F_CONTIGUOUS']:
            arr = arr.astype(dtype, order='K')
        else:
            print('_as_ptr: options.%s is neither C- nor F-contiguous; copying to Fortran order' % attr)
            arr = np.asfortranarray(arr, dtype=dtype)
        setattr(options, attr, arr)
    return arr.ctypes.data_as(ctypes.POINTER(ctype))


def transferData(options):
    """
    Transfers the Python variables to the corresponding C-struct
    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Returns
    -------
    None.

    """
    pure_odrt = (bool(getattr(options, 'SPECT', False)) and
                 getattr(options, 'projector_type', None) in (2, 12, 21, 22))
    odrt_custom = pure_odrt and bool(getattr(options, 'gFilterCustom', False))
    if bool(getattr(options, 'gFilterCustom', False)) and not pure_odrt:
        raise ValueError('Custom SPECT ODRT gFilter is supported only for SPECT projector types 2, 12, 21, and 22.')
    if odrt_custom:
        raw_filter = np.asarray(options.gFilter)
        if not np.issubdtype(raw_filter.dtype, np.number) or np.iscomplexobj(raw_filter):
            raise ValueError('ODRT gFilter must be a real numeric array.')
        try:
            with np.errstate(over='ignore', invalid='ignore'):
                odrt_filter = np.asarray(raw_filter, dtype=np.float32)
                odrt_spacing = np.asarray(options.gFilterSpacing, dtype=np.float32).reshape(-1)
        except (TypeError, ValueError) as exc:
            raise ValueError('ODRT gFilter and gFilterSpacing must be representable as float32 arrays.') from exc
        if odrt_filter.ndim == 2:
            odrt_filter = odrt_filter[:, :, None]
        if odrt_filter.ndim != 3 or min(odrt_filter.shape) < 1:
            raise ValueError('ODRT gFilter must be a non-empty 2-D or 3-D array with axes (u, v, depth).')
        if any(int(dim) > np.iinfo(np.uint32).max for dim in odrt_filter.shape):
            raise ValueError('ODRT gFilter dimensions must fit in uint32.')
        if (not np.all(np.isfinite(odrt_filter)) or np.any(odrt_filter < 0.) or
                not np.any(odrt_filter > 0.)):
            raise ValueError('ODRT gFilter must contain finite, non-negative weights and at least one positive value.')
        if (odrt_spacing.size != 3 or not np.all(np.isfinite(odrt_spacing)) or
                np.any(odrt_spacing <= 0.)):
            raise ValueError('gFilterSpacing must contain three positive finite values (du, dv, dd).')
        options.gFilterODRTData = np.asfortranarray(odrt_filter)
        options.gFilterSpacing = np.ascontiguousarray(odrt_spacing, dtype=np.float32)
        options.gFilterNu, options.gFilterNv, options.gFilterNd = map(int, odrt_filter.shape)
    else:
        options.gFilterODRTData = np.empty(0, dtype=np.float32)
        options.gFilterNu = options.gFilterNv = options.gFilterNd = 0
    # Loaders may use None for an absent optional correction. The native
    # interface always receives a pointer and an element count, so normalize
    # that representation to an empty array before taking either.
    if options.normalization is None:
        options.normalization = np.empty(0, dtype=np.float32)
    # Optional user-defined local/work-group (block) size. Accepts a scalar or a sequence of up to 3
    # values; missing/negative entries keep the built-in defaults.
    localSize = options.local_size
    if np.isscalar(localSize):
        localSize = [localSize]
    else:
        localSize = list(np.asarray(localSize).ravel())
    localSize = (localSize + [-1, -1, -1])[:3]

    if isinstance(options.inffi, np.ndarray):
        inffiVal = options.inffi.item()
    else:
        inffiVal = options.inffi

    # sizeX: number of listmode sensitivity-image source coordinates when computing
    # the sensitivity image in listmode, otherwise the number of x-coordinates.

    if options.listmode and options.compute_sensitivity_image:
        sizeXVal = options.uV.size
    else:
        sizeXVal = options.x.size

    if options.SinDelayed.dtype != np.float32:
        options.SinDelayed = options.SinDelayed.astype(np.float32)
    # A TOF-resolved randoms/scatter term has the size of the TOF measurements (the native library checks this as well)
    if not options.randoms_correction:
        options.randomsTOF = False
    if options.randomsTOF and options.SinDelayed.size != options.SinM.size:
        raise ValueError('The TOF-resolved randoms/scatter term (randomsTOF) has {} elements, but the TOF measurements have {}'.format(
            options.SinDelayed.size, options.SinM.size))

    detectorFrames = getattr(options, 'DetectorVectorFrames', None)
    if options.SPECT and isinstance(detectorFrames, list):
        options.DetectorVector = np.ascontiguousarray(np.concatenate(detectorFrames), dtype=np.uint32)

    attenuation_dims = np.zeros((2, 3), dtype=np.uint32)
    attenuation_spacing = np.zeros((2, 3), dtype=np.float32)
    attenuation_origin = np.zeros((2, 3), dtype=np.float32)
    if hasattr(options, 'imageAttenuationGridDims'):
        dims = np.asarray(options.imageAttenuationGridDims, dtype=np.uint32).reshape(-1, 3)
        spacing = np.asarray(options.imageAttenuationGridSpacing, dtype=np.float32).reshape(-1, 3)
        origin = np.asarray(options.imageAttenuationGridOrigin, dtype=np.float32).reshape(-1, 3)
        attenuation_dims[:dims.shape[0]] = dims
        attenuation_spacing[:spacing.shape[0]] = spacing
        attenuation_origin[:origin.shape[0]] = origin

    # Scalar struct fields with bespoke values (array sizes, conditional
    # selection, explicit type coercion, etc.) rather than a plain
    # options.<field name> (or renamed) lookup.
    scalarValueOverrides = {
        'localSizeX': int(localSize[0]),
        'localSizeY': int(localSize[1]),
        'localSizeZ': int(localSize[2]),
        'regEveryIter': int(options.regEveryIter),
        'inffi': inffiVal,
        'mDim': options.SinM.size // options.Nt,
        'nIterSaved': options.saveNIter.size,
        'sizeScat': options.corrVector.size,
        'sizeRandoms': options.SinDelayed.size,
        'eFOV': options.eFOVIndices.size,
        'sizeX': sizeXVal,
        'sizeZ': options.z.size,
        'sizeAtten': options.vaimennus.size,
        'sizeNorm': options.normalization.size,
        'sizeNormSens': options.normSens.size,
        'sizePSF': options.gaussK.size,
        'sizeXYind': options.xy_index.size,
        'sizeZind': options.z_index.size,
        'xCenterSize': options.x_center.size,
        'yCenterSize': options.y_center.size,
        'zCenterSize': options.z_center.size,
        'sizeV': options.V.size,
        'measElem': options.SinM.size,
        'sizeZSens': options.zSens.size,
        'sizeDetectorVector': np.size(options.DetectorVector),
        'sizeSensitivityViewWeights': options.sensitivityViewWeights.size,
        'imageAttenuationGridDims': tuple(attenuation_dims.ravel()),
        'imageAttenuationGridSpacing': tuple(attenuation_spacing.ravel()),
        'imageAttenuationGridOrigin': tuple(attenuation_origin.ravel()),
        'gFilterODRTNu': int(options.gFilterNu),
        'gFilterODRTNv': int(options.gFilterNv),
        'gFilterODRTNd': int(options.gFilterNd),
        'gFilterODRTDu': float(options.gFilterSpacing[0]) if odrt_custom else 1.0,
        'gFilterODRTDv': float(options.gFilterSpacing[1]) if odrt_custom else 1.0,
        'gFilterODRTDd': float(options.gFilterSpacing[2]) if odrt_custom else 1.0,
        'gFilterODRTCustom': int(odrt_custom),
    }

    setFields = set()
    for name, ctype in options.param._fields_:
        if name in _POINTER_SPECIAL_FIELDS:
            continue
        if issubclass(ctype, ctypes._Pointer):
            attr = _POINTER_NAME_OVERRIDES.get(name, name)
            setattr(options.param, name, _as_ptr(options, attr, ctype._type_))
        elif issubclass(ctype, ctypes.Array):
            value = scalarValueOverrides.get(name, getattr(options, name, ()))
            setattr(options.param, name, ctype(*value))
        else:
            if name in scalarValueOverrides:
                value = scalarValueOverrides[name]
            else:
                attr = _SCALAR_NAME_OVERRIDES.get(name, name)
                value = getattr(options, attr)
            setattr(options.param, name, ctype(value))
        setFields.add(name)

    # The native type-6 (rotation-dependent PSF blurring) branch uses
    # precomputed per-plane blur indices and the volume-0 filter/size.
    if options.projector_type in (6, 16, 26, 61, 62, 66):
        blur_planes = options.blurPlanes[0] if isinstance(options.blurPlanes, list) else options.blurPlanes
        blur_planes2 = options.blurPlanes2[0] if isinstance(options.blurPlanes2, list) else options.blurPlanes2
        g_filter = options.gFilter[0] if isinstance(options.gFilter, list) else options.gFilter
        options.param.blurPlanes = blur_planes.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
        options.param.blurPlanes2 = blur_planes2.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
        # The native type-6 branch uses the volume-0 filter.
        options.param.gFilter = g_filter.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        options.gFSize = np.array(g_filter.shape if g_filter.size else (0, 0, 0), dtype=np.uint64)
        # Sub-voxel geometry and analytic CDRF of the main volume; the arrays are kept on options so they stay alive
        if np.size(options.Nx) > 1 and not options.CDRFOnTheFly:
            raise ValueError('The native multi-resolution projector type 6 requires CDRFOnTheFly = True '
                             '(a manual gFilter/sigma cannot be used)')
        options.blurPlanesF = np.ascontiguousarray(np.asarray(options.blurPlanesF, dtype=np.float32).reshape(-1, np.asarray(options.blurPlanesF).shape[-1]))
        options.blurPlanes2F = np.ascontiguousarray(np.asarray(options.blurPlanes2F, dtype=np.float32).reshape(-1, np.asarray(options.blurPlanes2F).shape[-1]))
        options.CDRFCoeffs6 = np.ascontiguousarray(np.asarray(options.CDRFCoeffs6, dtype=np.float32).reshape(-1))
        options.param.blurPlanesF = options.blurPlanesF[0].ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        options.param.blurPlanes2F = options.blurPlanes2F[0].ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        options.param.CDRFCoeffs6 = options.CDRFCoeffs6.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        options.radiusPerProj = np.ascontiguousarray(options.radiusPerProj, dtype=np.float32).reshape(-1)
        options.param.radiusPerProj = options.radiusPerProj.ctypes.data_as(ctypes.POINTER(ctypes.c_float)) if options.radiusPerProj.size else None
        options.param.CORtoDetectorSurface = float(options.CORtoDetectorSurface)
        options.param.CDRFOnTheFly = ctypes.c_bool(bool(options.CDRFOnTheFly)).value
        # The model of the collimator response, 'gaussian' or 'square'; the latter is a property of the on-the-fly CDRF
        from omegatomo.projector.proj import cdrfModelIndex
        cdrfModel = cdrfModelIndex(options.CDRFModel)
        if cdrfModel != 0 and not options.CDRFOnTheFly:
            raise ValueError("options.CDRFModel = '%s' requires the on-the-fly collimator-detector response (CDRFOnTheFly = True, no manual gFilter/sigma)"
                             % options.CDRFModel)
        options.param.cdrfModel6 = ctypes.c_uint32(cdrfModel).value
        # How the backprojection rotates a view back into the image, 'adjoint' or 'fast'
        from omegatomo.projector.proj import bpRotationIndex
        options.param.bpRot6 = ctypes.c_uint32(bpRotationIndex(options.type6BPRotation)).value
    else:
        options.param.blurPlanes = None
        options.param.blurPlanes2 = None
        options.param.gFilter = None
        options.param.blurPlanesF = None
        options.param.blurPlanes2F = None
        options.param.CDRFCoeffs6 = None
        options.param.CDRFOnTheFly = False
        options.param.cdrfModel6 = 0
        options.param.bpRot6 = 0
        options.gFSize = np.zeros(3, dtype=np.uint64)
    options.param.gFSize = options.gFSize.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64))
    if odrt_custom:
        options.param.gFilterODRT = options.gFilterODRTData.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    else:
        options.param.gFilterODRT = None
    setFields.update(_POINTER_SPECIAL_FIELDS)

    missing = [name for name, _ in options.param._fields_ if name not in setFields]
    if missing:
        raise RuntimeError('transferData: struct field(s) left unset: %s' % missing)


def _appliedCorrectionsFlag(value):
    """
    MATLAB's `if value` for a logical field of the appliedCorrections struct:
    true only for a non-empty value whose elements are all nonzero.
    """
    if value is None:
        return False
    value = np.asarray(value)
    return value.size > 0 and bool(np.all(value))


def _appliedCorrectionsText(value):
    """
    A char field of the appliedCorrections struct as a str (an empty field,
    which pymatreader returns as an empty array, gives '').
    """
    if value is None:
        return ''
    if isinstance(value, str):
        return value
    value = np.asarray(value)
    if value.size == 0 or value.dtype.kind not in ('U', 'S', 'O'):
        return ''
    return ''.join(str(v) for v in value.ravel())


def _resolveAppliedCorrections(options, appliedCorrections, loadRaw):
    """
    Port of the appliedCorrections block of MATLAB's loadInputData.m (corrections
    not applied during the reconstruction). A measurement file can store the
    precorrected measurements SinM together with the struct appliedCorrections
    that tells which corrections SinM has already been precorrected with
    (.normalization, .randoms and .scatter; the text of the last two names the
    applied options, e.g. 'randoms correction with variance reduction and
    smoothing'). From that and the selected corrections this decides whether SinM
    can be used as it is or the raw measurements (raw_SinM) have to be
    precorrected here, and which corrections are still left to apply: it updates
    options.normalization_correction, options.randoms_correction and
    options.scatter_correction accordingly (a correction that is already in SinM is
    turned off). The normalization is only re-enabled (when the raw data replaces
    SinM that was precorrected with it) if the user has selected it. The fields mashing and globalFactor are not used (as in
    MATLAB's loadInputData.m). The field gapFilling is not used here either, the sinogram gap filling that follows from it is done
    in _selectPrecorrectedMeasurement.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.
    appliedCorrections : dict
        The appliedCorrections struct of the file.
    loadRaw : bool
        True if the raw data is already used (SinM was not in the file).

    Returns
    -------
    bool
        True if the raw measurements (raw_SinM) have to be used.
    """
    normalizationApplied = _appliedCorrectionsFlag(appliedCorrections.get('normalization'))
    randomsApplied = _appliedCorrectionsText(appliedCorrections.get('randoms'))
    scatterApplied = _appliedCorrectionsText(appliedCorrections.get('scatter'))
    normalization_correction = bool(options.normalization_correction)
    # The normalization is only (re-)enabled below if the user has selected it
    normalizationSelected = bool(options.normalization_correction)
    # Normalization
    if normalizationApplied and not options.normalization_correction:
        print('Warning: Normalization correction not selected, but data precorrected with normalization! Precorrecting without normalization.')
        loadRaw = True
    if normalizationApplied:
        normalization_correction = False
    # Randoms: check if the data is precorrected with something that has not been selected
    randoms_correction = bool(options.randoms_correction)
    if 'variance reduction' in randomsApplied and not options.variance_reduction and options.randoms_correction:
        print('Warning: Randoms variance correction not selected, but data precorrected with randoms with applied variance reduction! Precorrecting without variance reduction.')
        loadRaw = True
        if normalizationApplied and not normalization_correction and normalizationSelected:
            normalization_correction = True
    elif 'smoothing' in randomsApplied and not options.randoms_smoothing and options.randoms_correction:
        print('Warning: Randoms smoothing not selected, but data precorrected with randoms with applied smoothing! Precorrecting without smoothing.')
        loadRaw = True
        if normalizationApplied and not normalization_correction and normalizationSelected:
            normalization_correction = True
    elif 'randoms correction' not in randomsApplied and not loadRaw and randoms_correction:
        # Randoms not in the precorrected data: precorrect the raw data
        loadRaw = True
        if normalizationApplied and not normalization_correction and normalizationSelected:
            normalization_correction = True
    elif ('randoms correction' in randomsApplied and randoms_correction
          and ((options.variance_reduction and 'variance reduction' not in randomsApplied)
               or (options.randoms_smoothing and 'smoothing' not in randomsApplied))):
        print('Warning: Randoms corrections selected, but data not precorrected with selected options! Precorrecting with selected options.')
        loadRaw = True
        if normalizationApplied and not normalization_correction and normalizationSelected:
            normalization_correction = True
    elif 'randoms correction' in randomsApplied and not loadRaw:
        randoms_correction = False
    # Same as above, but for scatter
    if 'variance reduction' in scatterApplied and not options.scatter_variance_reduction and options.scatter_correction:
        print('Warning: Scatter variance correction not selected, but data precorrected with scatter with applied variance reduction! Precorrecting without variance reduction.')
        if not randoms_correction and options.randoms_correction:
            randoms_correction = True
        if normalizationApplied and not normalization_correction and normalizationSelected:
            normalization_correction = True
        loadRaw = True
    elif 'smoothing' in scatterApplied and not options.scatter_smoothing and options.scatter_correction:
        print('Warning: Scatter smoothing not selected, but data precorrected with scatter with applied smoothing! Precorrecting without smoothing.')
        if not randoms_correction and options.randoms_correction:
            randoms_correction = True
        if normalizationApplied and not normalization_correction and normalizationSelected:
            normalization_correction = True
        loadRaw = True
    elif 'scatter correction' not in scatterApplied and not loadRaw and options.scatter_correction:
        if not randoms_correction and options.randoms_correction:
            randoms_correction = True
        if normalizationApplied and not normalization_correction and normalizationSelected:
            normalization_correction = True
        loadRaw = True
    elif 'scatter correction' in scatterApplied and not loadRaw:
        options.scatter_correction = False
    options.randoms_correction = randoms_correction
    options.normalization_correction = normalization_correction
    return loadRaw


def _gapFillingSelected(options):
    """
    True if sinogram gap filling is selected and applicable: options.fill_sinogram_gaps is set and the scanner has pseudo
    detectors (det_per_ring < det_w_pseudo), as in MATLAB's loadInputData.m. Not applicable to CT and SPECT.
    """
    if not options.fill_sinogram_gaps or options.CT or options.SPECT:
        return False
    return int(np.asarray(options.det_per_ring).ravel()[0]) < int(np.asarray(options.det_w_pseudo).ravel()[0])


def _gapFillMeasurement(options, sinM):
    """
    Sinogram gap filling (omegatomo.util.gapFilling) of a measurement array (or a list with one array per time frame) that was
    loaded from a file. The filled data are new arrays (float32), the input, which can be shared with the caller, is not modified.
    A flat array is treated as (Ndist, Nang, remaining dimensions) in Fortran order.
    """
    from omegatomo.util.gapFilling import gapFilling
    if isinstance(sinM, (list, tuple)):
        filled = list()
        gaps = None
        for llo, frame in enumerate(sinM):
            frame, gaps = gapFilling(options, _gapFillInput(options, frame), llo + 1, gaps)
            filled.append(frame)
        return filled
    sinM, _ = gapFilling(options, _gapFillInput(options, sinM))
    return sinM


def _gapFillInput(options, sinM):
    # A flat measurement array as (Ndist, Nang, remaining dimensions), the layout that gapFilling needs
    sinM = np.asarray(sinM)
    if sinM.ndim == 1:
        nSinogram = int(options.Ndist) * int(options.Nang)
        if nSinogram == 0 or sinM.size % nSinogram != 0:
            raise ValueError('The number of elements of the measurement data (%d) is not divisible by Ndist * Nang (%d), gap filling is not possible!' % (sinM.size, nSinogram))
        sinM = np.reshape(sinM, (int(options.Ndist), int(options.Nang), -1), order='F')
    return sinM


def _selectPrecorrectedMeasurement(options, getKey, appliedCorrections=None):
    """
    Selects either the corrected ('SinM') or the raw ('raw_SinM') measurement
    data key, matching the precorrection semantics shared by the MAT and NPZ
    measurement-file loaders. `getKey(key)` returns the array for that key,
    or raises KeyError.

    If the measurement file contains the struct appliedCorrections (a dict,
    as read from a MAT-file) and the data has not been explicitly requested to be
    precorrected here (options.precorrect), the choice and the corrections that are
    still applied follow MATLAB's loadInputData.m, see _resolveAppliedCorrections.
    Without appliedCorrections the stored SinM is assumed to be precorrected with
    everything that has been selected.

    Sinogram gap filling (options.fill_sinogram_gaps with pseudo detectors, as in loadInputData.m): the raw measurements are
    never gap filled in the file and are always filled here. The corrected SinM is filled here only if appliedCorrections is
    available and says that it has not been gap filled (appliedCorrections.gapFilling is false). Without appliedCorrections the
    stored SinM is used as it is.
    """
    useRaw = True
    sinM = None
    if (options.randoms_correction or options.scatter_correction or options.normalization_correction) and not options.corrections_during_reconstruction:
        if not options.precorrect:
            loadRaw = False
            try:
                sinM = getKey('SinM')
            except KeyError:
                options.precorrect = True
                loadRaw = True
            if isinstance(appliedCorrections, dict):
                loadRaw = _resolveAppliedCorrections(options, appliedCorrections, loadRaw)
                if loadRaw or options.randoms_correction or options.scatter_correction or options.normalization_correction:
                    # OMEGA has to apply the corrections that are not in SinM
                    options.precorrect = True
            useRaw = loadRaw
    if useRaw:
        sinM = getKey('raw_SinM')
    if _gapFillingSelected(options):
        if useRaw or (isinstance(appliedCorrections, dict) and not _appliedCorrectionsFlag(appliedCorrections.get('gapFilling'))):
            # New arrays are created: the measurement arrays can be shared with the caller
            sinM = _gapFillMeasurement(options, sinM)
    return sinM


def _nativeLibPath(libdir, name):
    """
    Returns the full path to the native reconstruction library `name` inside
    `libdir`, using the platform-appropriate shared-library extension
    ('.dll' on Windows, '.so' everywhere else).
    """
    import os
    ext = '.dll' if os.name == 'nt' else '.so'
    return str(os.path.join(libdir, name + ext))


def _delayedVariableName(options, names, allowFirst):
    """
    Selects the variable that holds the randoms (delayed coincidences) from the variable names `names` of a file, as in MATLAB's
    loadDelayedData.m: raw_SinDelayed, then SinDelayed (sinogram data), then delayed_coincidences (raw data), then, if allowFirst is
    set, the first variable. Returns None if no variable is selected.
    """
    names = list(names)
    if 'raw_SinDelayed' in names and not options.use_raw_data:
        return 'raw_SinDelayed'
    if 'SinDelayed' in names and not options.use_raw_data:
        return 'SinDelayed'
    if 'delayed_coincidences' in names and options.use_raw_data:
        return 'delayed_coincidences'
    if allowFirst and len(names) > 0:
        return names[0]
    return None


def _loadDelayedMeasurement(options, getKey, names):
    """
    Loads the randoms data from an already-open measurement file, matching the
    semantics shared by the MAT and NPZ measurement-file loaders. `getKey(key)`
    returns the array for that key and `names` are the variable names of the
    file. The first variable is not used as a fallback here since it is the
    measurement data itself (see _delayedVariableName for the used variables).
    """
    if options.randoms_correction and not options.reconstruct_scatter and not options.reconstruct_trues and options.SinDelayed.size < 1:
        name = _delayedVariableName(options, names, False)
        if name is None:
            print('Randoms correction selected but no randoms data found. The randoms data should be saved as SinDelayed')
        else:
            options.SinDelayed = getKey(name)


def _frameArray(options, frame, nDims):
    """
    A frame of the randoms/scatter data of a dynamic reconstruction as an array with nDims dimensions (Ndist, Nang, planes[, TOF
    bins]). A flat frame is reshaped to (Ndist, Nang, planes) and missing trailing dimensions are added as singleton dimensions.
    """
    frame = np.asarray(frame)
    if frame.ndim == 1:
        nSinogram = int(options.Ndist) * int(options.Nang)
        if nSinogram == 0 or frame.size % nSinogram != 0:
            raise ValueError('The number of elements of a randoms/scatter frame (%d) is not divisible by Ndist * Nang (%d)!' % (frame.size, nSinogram))
        frame = np.reshape(frame, (int(options.Ndist), int(options.Nang), -1), order='F')
    if frame.ndim > nDims:
        raise ValueError('A randoms/scatter frame has %d dimensions, but at most %d are supported!' % (frame.ndim, nDims))
    while frame.ndim < nDims:
        frame = frame[..., np.newaxis]
    return frame


def _expandFrameCorrections(options):
    """
    Dynamic data (Nt > 1): converts the randoms (options.SinDelayed) and the scatter (options.ScatterC) to the layouts of the
    reconstruction, the time frame is the last dimension: SinDelayed (Ndist, Nang, planes, Nt) and ScatterC (Ndist, Nang, planes,
    TOF bins, Nt). A list with one array for each frame is stacked, while a single-frame array (one sinogram for all the frames) is
    used for every frame, as in MATLAB's loadInputData.m and loadCorrections.m. If there are fewer randoms frames than time frames,
    the first one is used for all (with a warning). The data are rebound, never written in place, since the arrays can be the
    caller's arrays. A list of SPECT scatter data (ScatterC) holds the energy windows, not the frames, and is not touched. Dynamic
    list-mode and index-based data are not touched either.
    """
    nFrames = int(options.Nt)
    if nFrames <= 1 or options.listmode or options.useIndexBasedReconstruction or options.reconstruct_scatter or options.reconstruct_trues:
        return
    nSinogram = int(options.Ndist) * int(options.Nang)
    if options.randoms_correction and not options.randomsTOF:
        delayed = options.SinDelayed
        if isinstance(delayed, (list, tuple)):
            if len(delayed) > nFrames:
                raise ValueError('The number of delayed coincidence sinograms (%d) is larger than the number of time points (%d)!' % (len(delayed), nFrames))
            if len(delayed) > 0 and np.size(delayed[0]) > 1:
                frames = list(delayed)
                if len(frames) < nFrames:
                    print('Warning: The number of delayed coincidence sinograms is less than the number of time points. Using the first one')
                    frames = [frames[0]] * nFrames
                options.SinDelayed = np.stack([_frameArray(options, frame, 3) for frame in frames], axis=3)
        elif np.size(delayed) > 1:
            delayed = np.asarray(delayed)
            single = False
            if delayed.ndim == 1:
                single = delayed.size in (nSinogram * int(options.NSinos), nSinogram * int(options.TotSinos))
            elif delayed.ndim <= 3:
                single = True
            elif delayed.ndim == 4:
                single = delayed.shape[3] != nFrames
            if single:
                print('Warning: The number of delayed coincidence sinograms is less than the number of time points. Using the first one')
                if delayed.ndim == 4:
                    delayed = delayed[..., 0]
                frame = _frameArray(options, delayed, 3)
                options.SinDelayed = np.stack([frame] * nFrames, axis=3)
    if options.scatter_correction and not options.SPECT:
        scatter = options.ScatterC
        if isinstance(scatter, (list, tuple)):
            if len(scatter) > 0 and np.size(scatter[0]) > 1:
                if len(scatter) != nFrames:
                    raise ValueError('The number of scatter sinograms (%d) is different from the number of time points (%d)!' % (len(scatter), nFrames))
                options.ScatterC = np.stack([_frameArray(options, frame, 4) for frame in scatter], axis=4)
        elif np.size(scatter) > 1:
            scatter = np.asarray(scatter)
            tofBins = max(int(options.TOF_bins_used), 1)
            single = False
            if scatter.ndim == 1:
                single = scatter.size in {nSinogram * int(nSinos) * bins for nSinos in (options.NSinos, options.TotSinos) for bins in {1, tofBins}}
                if single:
                    # A flat single-frame scatter sinogram: TOF bins are the slowest dimension if there are any
                    if scatter.size == nSinogram * int(options.NSinos) * tofBins or scatter.size == nSinogram * int(options.TotSinos) * tofBins:
                        scatter = np.reshape(scatter, (int(options.Ndist), int(options.Nang), -1, tofBins), order='F')
                    else:
                        scatter = np.reshape(scatter, (int(options.Ndist), int(options.Nang), -1, 1), order='F')
            elif scatter.ndim <= 3:
                single = True
            elif scatter.ndim == 4:
                if scatter.shape[3] == nFrames and nFrames != tofBins:
                    # Non-TOF frames (Ndist, Nang, planes, Nt)
                    options.ScatterC = scatter[:, :, :, np.newaxis, :]
                else:
                    single = True
            elif scatter.ndim == 5:
                single = scatter.shape[4] != nFrames
                if single:
                    scatter = scatter[..., 0]
            if single:
                # One scatter sinogram shared by all the frames
                frame = _frameArray(options, scatter, 4)
                options.ScatterC = np.stack([frame] * nFrames, axis=4)



def reconstructions_mainCT(options):
    """
    This function simply does certain CT-specific adjustments before calling
    the main built-in reconstruction function
    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Returns
    -------
    pz : NumPy Array
        the reconstructed image volume.
    FPOutputP : NumPy Array
        The (optional) forward projections.
    residual : NumPy Array
        the (optional) residual/primal-dual gap.
    """
    options.CT = True
    return reconstructions_main(options)

def reconstructions_mainSPECT(options):
    """
    This function simply does certain SPECT-specific adjustments before calling
    the main built-in reconstruction function
    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Returns
    -------
    pz : NumPy Array
        the reconstructed image volume.
    FPOutputP : NumPy Array
        The (optional) forward projections.
    """
    options.SPECT = True
    return reconstructions_main(options)

def reconstructions_main(options):
    """
    The main built-in reconstruction function
    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Returns
    -------
    pz : NumPy Array
        the reconstructed image volume.
    FPOutputP : NumPy Array
        The (optional) forward projections.
    residual : NumPy Array
        the (optional) residual/primal-dual gap.

    Notes
    -----
    After the reconstruction the settings used are stored as the dictionary options.image_properties (the equivalent of the
    second output of MATLAB's reconstructions_main), see imageProperties(). The dictionary is stored in the object that was given
    to this function (and in the private working copy, if options.keepOptionsUnchanged is True). It can be given as the options
    input of saveInterfile/saveMetaImage/saveImage.
    """
    import time
    import os
    from .prepass import prepassPhase
    from .prepass import parseInputs
    from .prepass import loadCorrections
    from .prepass import _prepare_image_domain_attenuation
    tic = time.perf_counter()
    # Work on a private copy (MATLAB value semantics): the pipeline below changes the
    # options in place and not idempotently, so repeated calls with the same object
    # must always start from the user's original values (opt out: keepOptionsUnchanged = False).
    callerOptions = options
    if getattr(options, 'keepOptionsUnchanged', True):
        options = _privateOptionsCopy(options)
    if options.CT and np.size(options.flat) > 1:
        # An array-valued flat field (e.g. one value for each projection or detector pixel) is reduced to its maximum, as in MATLAB
        options.flat = np.max(np.asarray(options.flat))
    if options.SPECT and int(options.projector_type) == 6:
        # As in MATLAB's projectorClass: pad/crop and resample the SPECT measurements to the Nx x Nz grid of the rotation-based projector.
        # Done before addProjector(), which derives the detector-size dependent values and flips the projections.
        from omegatomo.util.measurementResampling import _resampleType6Measurements
        _resampleType6Measurements(options)
    options.addProjector()
    # A7: per-frame lists / one frame of the randoms and scatter data to the layouts of the dynamic reconstruction
    _expandFrameCorrections(options)
    # These projector types are accepted by the projector class (custom operators), but the native library has no
    # forward/backward projector mapping for them and would silently fall back to Siddon (MATLAB errors)
    if int(options.projector_type) in (16, 26, 61, 62, 66):
        raise ValueError('Projector type ' + str(options.projector_type) + ' is only available through the projector operators, '
                         'not the built-in reconstruction (reconstructions_main)!')
    print('Preparing for reconstruction...')
    if not options.builtin:
        raise ValueError('No reconstruction method selected, aborting.')
    if np.size(options.weights) > 0:
        options.empty_weight = False
    fname, suffix = os.path.splitext(options.fpath)
    if isinstance(options.SinM, list):
        if len(options.SinM) == 0:
            sinoSize = 0
        else:
            sinoSize = options.SinM[0].size
    else:
        sinoSize = options.SinM.size
    if sinoSize < 1 and (len(options.fpath) == 0 or len(suffix) == 0):
        import tkinter as tk
        from tkinter.filedialog import askopenfilename
        root = tk.Tk()
        root.withdraw()
        options.fpath = askopenfilename(title='Select measurement datafile',filetypes=(('NPY, NPZ and MAT files','*.mat *.npy *.npz'),('All','*.*')))
        if len(options.fpath) == 0:
            raise ValueError('No file selected')
        fname, suffix = os.path.splitext(options.fpath)
    if sinoSize < 1 and suffix == '.mat':
        from pymatreader import read_mat
        try:
            var = read_mat(options.fpath)
        except OSError:
            print('File not found, please select the measurement data file')
            import tkinter as tk
            from tkinter.filedialog import askopenfilename
            root = tk.Tk()
            root.withdraw()
            options.fpath = askopenfilename(title='Select measurement datafile',filetypes=(('MAT files','*.mat'),('All','*.*')))
            if len(options.fpath) == 0:
                raise ValueError('No file selected')
            var = read_mat(options.fpath)
        matGet = lambda key: np.array(var[key], order='F')
        if options.reconstruct_trues or options.reconstruct_scatter:
            options.SinM = matGet('SinTrues') if options.reconstruct_trues else matGet('SinScatter')
            if options.normalization_correction and not options.corrections_during_reconstruction:
                # As in MATLAB's loadCorrections.m, the trues/scatter data is precorrected
                # with the normalization (multiplied by it) instead of applying the
                # normalization during the reconstruction
                options.precorrect = True
            if _gapFillingSelected(options):
                if options.verbose:
                    print('Performing sinogram gap filling on ' + ('trues' if options.reconstruct_trues else 'scatter') + ' data')
                options.SinM = _gapFillMeasurement(options, options.SinM)
        else:
            options.SinM = _selectPrecorrectedMeasurement(options, matGet, var.get('appliedCorrections'))
        _loadDelayedMeasurement(options, matGet, var.keys())
    elif sinoSize < 1 and suffix == '.npy':
        options.SinM = np.load(options.fpath)
    elif sinoSize < 1 and suffix == '.npz':
        varList = np.load(options.fpath)
        npzGet = lambda key: varList[key]
        options.SinM = _selectPrecorrectedMeasurement(options, npzGet)
        _loadDelayedMeasurement(options, npzGet, varList.files)
    elif not options.corrections_during_reconstruction and not options.precorrect and (options.randoms_correction or options.scatter_correction or options.normalization_correction):
        if options.precorrect is None:
            options.precorrect = True
            print('Corrections selected and options.corrections_during_reconstruction = False: the input measurement data is precorrected (as in MATLAB). Set options.precorrect = False if the data has already been corrected.')
        else:
            print('Corrections selected and measurement data found, but options.precorrect = False was set explicitly. The input measurement data WILL NOT BE PRECORRECTED!!!!!! If you wish to have OMEGA-based precorrection, make sure options.precorrect = True')
    if options.precorrect is None:
        options.precorrect = False
    if options.randoms_correction and not options.reconstruct_scatter and not options.reconstruct_trues and options.SinDelayed.size < 1:
        import tkinter as tk
        from tkinter.filedialog import askopenfilename
        root = tk.Tk()
        root.withdraw()
        fpath = askopenfilename(title='Select randoms datafile',filetypes=(('NPY, NPZ and MAT files','*.mat *.npy *.npz'),('All','*.*')))
        if len(fpath) == 0:
            print('No file selected, disabling randoms correction')
            options.randoms_correction = False
        else:
            _, fsuffix = os.path.splitext(fpath)
            if fsuffix == '.mat' and options.randoms_correction:
                from pymatreader import read_mat
                var = read_mat(fpath)
                # As in MATLAB's loadDelayedData.m: raw_SinDelayed, SinDelayed, delayed_coincidences or otherwise the first variable
                name = _delayedVariableName(options, var.keys(), True)
                if name is None:
                    print('Randoms correction selected but no randoms data found. The randoms data should be saved as SinDelayed. Disabling randoms correction')
                    options.randoms_correction = False
                else:
                    options.SinDelayed = np.array(var[name],order='F')
            elif fsuffix == '.npy':
                options.SinDelayed = np.load(fpath)
            elif fsuffix == '.npz':
                varList = np.load(fpath)
                name = _delayedVariableName(options, varList.files, True)
                if name is None:
                    print('Randoms correction selected but no randoms data found. The randoms data should be saved as SinDelayed. Disabling randoms correction')
                else:
                    options.SinDelayed = varList[name]
    # The randoms/scatter data that were loaded from the files above
    _expandFrameCorrections(options)
    if options.TOF and options.TOF_bins_used == 1:
        options.TOF_bins = options.TOF_bins_used
        if isinstance(options.SinM, list):
            # Same list-stacking hazard as the SinM reshape below: np.sum on the raw
            # list would first stack frames along a new leading axis (via np.asarray),
            # shifting axis=3 off the intended TOF_bins axis onto nProjections instead.
            # Sum each frame individually so the TOF axis is reduced correctly per frame.
            options.SinM = [np.sum(np.asarray(frame), axis=3) for frame in options.SinM]
        else:
            options.SinM = np.sum(options.SinM, axis=3)
        options.TOF = False
    loadCorrections(options)
    # if options.normalization_correction and options.corrections_during_reconstruction == True:
    #     normdir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', '..', '..', '..', 'mat-files')) + "/"
    #     if os.path.exists(normdir):
    #         var = read_mat(options.fpath)
    if options.CT and options.flat <= 0 and not options.usingLinearizedData:
        if (options.normalization_correction and isinstance(options.SinM, np.ndarray) and options.normalization.size > 1
                and options.normalization.size == options.SinM.size):
            # CT normalization kept for the reconstruction (a gain of each measurement): the flat is the maximum of the normalized
            # raw data. If the data has been divided by the normalization (precorrection), the maximum of the data is the same value
            print('No flat value input! Using the maximum value of the normalized data, max(y / n), as the flat value. Alternatively, input the flat value into options.flat')
            ctNormalization = options.normalization.ravel('F').astype(np.float64)
            ctNormalizationValid = ctNormalization > 0
            if not np.any(ctNormalizationValid):
                raise ValueError('CT normalization has no positive values')
            options.flat = np.max(options.SinM.ravel('F').astype(np.float64)[ctNormalizationValid] / ctNormalization[ctNormalizationValid]).astype(dtype=np.float32)
        else:
            print('No flat value input! Using the maximum value as the flat value. Alternatively, input the flat value into options.flat')
            options.flat = np.max(options.SinM).astype(dtype=np.float32)
    if options.Nt <= 1 and not options.CT and not options.SPECT and not options.SinM.size == options.Ndist * options.Nang * options.TotSinos * options.Nt * options.TOF_bins_used and options.listmode == 0:
        raise ValueError('The number of elements in the input data does not match the input number of angles, radial distances and total number of sinograms multiplied together!')
    if not options.usingLinearizedData and (options.LSQR or options.CGLS or options.FISTA or options.FISTAL1 or options.PDHG or options.PDHGL1 or options.PDDY or options.FDK or options.SART or options.ASD_POCS or options.BB) and not options.largeDim and options.CT:
        from .prepass import linearizeData
        linearizeData(options)
        options.usingLinearizedData = True
    if options.useParkerWeights:
        from omegatomo.util.parkerWeights import ParkerWeights
        ParkerWeights(options)
    # Offset-detector FDK modifies measurements and geometry only during
    # the native reconstruction call; restore the user's options in finally.
    offsetFDKState = None
    if options.CT and options.FDK and options.offsetCorrection:
        from .offsetfdk import offsetFDKPreprocess
        offsetFDKState = offsetFDKPreprocess(options)
    try:
        if not options.listmode:
            if isinstance(options.SinM, list):
                # Preserve each dynamic frame as one contiguous Fortran-ordered block.
                options.SinM = np.concatenate([np.asarray(frame).ravel(order='F') for frame in options.SinM])
            options.SinM = np.reshape(
                options.SinM,
                (int(options.nRowsD), int(options.nColsD), options.nProjections, options.TOF_bins, options.Nt),
                order='F',
            )
        elif options.listmode and options.compute_sensitivity_image and options.SPECT:
            from omegatomo.projector.detcoord import getCoordinatesSPECT
            x_sensitivity, z_sensitivity = getCoordinatesSPECT(options)
            options.uV = np.float32(np.asfortranarray(x_sensitivity))
            options.zSens = np.float32(np.asfortranarray(z_sensitivity))
            options.sensitivityViewWeights = _spect_listmode_sensitivity_weights(
                options, options.uV.size // 6
            )
        elif options.listmode and options.compute_sensitivity_image:
            if (hasattr(options, 'xSens') and np.size(options.xSens) > 0
                    and hasattr(options, 'zSens') and np.size(options.zSens) > 0):
                options.uV = np.float32(np.asfortranarray(options.xSens))
                options.z = np.float32(np.asfortranarray(options.zSens))
                options.det_per_ring = options.uV.size // 2 // options.nLayers
                options.rings = options.z.size // options.nLayers
            else:
                from omegatomo.projector.detcoord import getCoordinates
                options.use_raw_data = True
                x, y, z = getCoordinates(options)
                options.use_raw_data = False
                options.uV = np.float32(x)
                options.z = np.float32(z)
        if (options.quad or options.FMH or options.L or options.weighted_mean or options.Huber or options.GGMRF) and options.MAP:
            if hasattr(options, 'weights') and np.size(options.weights) > 0:
                weights_flat = np.array(options.weights).flatten()
                expected_length = ((options.Ndx * 2 + 1) *
                    (options.Ndy * 2 + 1) *
                    (options.Ndz * 2 + 1))
                if len(weights_flat) < expected_length:
                    raise ValueError(
                        f'Weights vector is too short, needs to be {expected_length} in length'
                    )
                elif len(weights_flat) > expected_length:
                    raise ValueError(
                        f'Weights vector is too long, needs to be {expected_length} in length'
                    )
                middle_index = int(np.ceil(expected_length / 2)) - 1
                if not np.isinf(weights_flat[middle_index]):
                    weights_flat[middle_index] = np.inf
                options.weights = weights_flat
            else:
                options.empty_weight = True
        parseInputs(options, True)
        if (options.SPECT and options.attenuation_correction and options.CT_attenuation
                and options.useMultiResolutionVolumes and options.nMultiVolumes > 0):
            _prepare_image_domain_attenuation(options)
        # Each list entry is one timestep's frame, which can still be a multi-dimensional
        # (Ndist x Nang x NSinos(x TOF_bins)) array here: parseInputs only flattens per-frame
        # data when subsets > 1 (and subsetType > 0); with subsets == 1 (or largeDim) the
        # per-partition loop above never runs, so np.concatenate would otherwise join the
        # frames along axis 0 of their native (row-major-read, Fortran-ordered) shape instead
        # of stacking whole frames end to end -- interleaving timesteps instead of keeping
        # each one contiguous (the C++ side offsets into SinM/SinDelayed by that array's
        # per-timestep element count times the timestep index, i.e. expects frame 1 fully
        # before frame 2). Raveling each frame with order='F' first (a no-op if already 1-D,
        # as it is once the subsets > 1 path above already ran) guarantees frame-contiguous
        # memory order regardless of which path produced the list.
        if isinstance(options.SinM, list):
            options.SinM = np.concatenate([np.asarray(frame).ravel(order='F') for frame in options.SinM])
        if isinstance(options.SinDelayed, list):
            options.SinDelayed = np.concatenate([np.asarray(frame).ravel(order='F') for frame in options.SinDelayed])
        if not options.CT and (not options.LSQR and not options.CGLS):
            negativeMask = options.SinM < 0
            if negativeMask.any():
                if _sharesCallerMemory(options, options.SinM):
                    # SinM can be the caller's array: out-of-place
                    options.SinM = np.where(negativeMask, options.SinM.dtype.type(0), options.SinM)
                else:
                    options.SinM[negativeMask] = 0
        if options.FDK:
            options.precondTypeMeas[1] = True
        prepassPhase(options)
        options.tau = 2.5
        if options.use_32bit_atomics and options.use_64bit_atomics:
            options.use_64bit_atomics = False
        if options.use_64bit_atomics and (options.useCPU or options.useCUDA):
            options.use_64bit_atomics = False
        if options.use_32bit_atomics and (options.useCPU or options.useCUDA):
            options.use_32bit_atomics = False
        if options.storeMultiResolution:
            output = np.zeros(int(np.sum(options.N) * options.Nt), dtype=np.float32, order = 'F')
        elif options.useMultiResolutionVolumes:
            output = np.zeros(options.NxOrig * options.NyOrig * options.NzOrig * options.Nt, dtype=np.float32, order = 'F')
        else:
            output = np.zeros(options.Nx[0].item() * options.Ny[0].item() * options.Nz[0].item() * options.Nt, dtype=np.float32, order = 'F')
        if options.saveNIter.size > 0:
            output = np.tile(output, options.saveNIter.size + 1)
        elif options.save_iter:
            output = np.tile(output, options.Niter + 1)
        if options.storeFP:
            FPOutput = np.zeros(options.SinM.size * options.Niter, dtype=np.float32, order = 'F')
        else:
            FPOutput = np.empty(0, dtype=np.float32)
        if options.storeResidual:
            residual = np.zeros(options.Niter * options.subsets, dtype=np.float32)
        else:
            residual = np.zeros(1, dtype=np.float32)
        fPath = os.path.dirname( __file__ )
        if os.path.exists(os.path.join(fPath, '..', 'util', 'usingPyPi.py')):
            libdir = os.path.join(os.path.abspath(os.path.join(fPath, '..')), "libs")
        else:
            libdir = os.path.abspath(os.path.join(fPath, '..', '..'))
        from omegatomo.util.paths import opencl_header_dir
        options.headerDir = opencl_header_dir()
        # The native code reads one initial value block per (multi-resolution) volume from x0; an
        # undersized x0 would make it read past the end of the array
        if not options.largeDim and not options.FDK and options.x0.size < int(np.sum(options.N)):
            raise ValueError(f'Initial value x0 has {options.x0.size} elements but the reconstruction volume(s) require {int(np.sum(options.N))}!')
        # Per-frame (time-varying) normalization or measurement-domain attenuation, i.e. more than one frame of measurements, is not supported
        # by the library: it builds them (and the sensitivity images) from the first frame only, so every frame would be reconstructed with the
        # corrections of the first frame. Static (single-frame) arrays are shared by all the frames and image-based attenuation can differ
        # between the frames. The projector operators (forwardProject/backwardProject) do support per-frame corrections, so this is checked
        # here and not in OMEGAErrorCheck
        if options.Nt > 1 and not options.listmode and not options.largeDim:
            frameMeas = int(options.nTotMeas[options.subsets]) - int(options.nTotMeas[0])
            headNorm = options.SPECT and int(options.normZ) == int(options.nHeads)
            if options.normalization_correction and np.size(options.normalization) > frameMeas and not headNorm:
                raise ValueError('Per-frame (time-varying) normalization is not supported by reconstructions_main with dynamic data: every frame would use '
                                 'the normalization of the first frame. Supply a single-frame normalization, which is shared by all the frames, or use the '
                                 'projector operators (forwardProject/backwardProject), which support per-frame corrections.')
            if options.attenuation_correction and not options.CT_attenuation and np.size(options.vaimennus) > frameMeas:
                raise ValueError('Per-frame (time-varying) measurement-domain attenuation is not supported by reconstructions_main with dynamic data: every '
                                 'frame would use the attenuation of the first frame. Supply a single-frame attenuation, which is shared by all the frames, '
                                 'an image-based attenuation (which can differ between the frames), or use the projector operators '
                                 '(forwardProject/backwardProject), which support per-frame corrections.')
        transferData(options)
        inStr = options.headerDir.encode('utf-8')
        # point_ptr = ctypes.pointer(options.param)
        # Same frame-contiguity concern as the earlier SinM/SinDelayed concatenation above
        # (this list form can reappear here for the Nt <= 1 case via prepassPhase's
        # single_frame_list re-wrapping); ravel each frame with order='F' before joining them.
        if isinstance(options.SinM, list):
            options.SinM = np.concatenate([np.asarray(frame).ravel(order='F') for frame in options.SinM])
        keep_int = options.SinM.dtype in (np.uint16, np.uint8) and (options.largeDim or not options.loadTOF)
        if not keep_int and options.SinM.dtype != np.float32:
            options.SinM = options.SinM.astype(np.float32)
        if options.SinM.ndim > 1:
            options.SinM = options.SinM.ravel('F')
        if options.useCUDA:
            if options.SinM.dtype == 'uint16':
                libN = 'CUDA_matrixfree_uint16_lib'
            elif options.SinM.dtype == 'uint8':
                libN = 'CUDA_matrixfree_uint8_lib'
            else:
                libN = 'CUDA_matrixfree_lib'
            libname = _nativeLibPath(libdir, libN)
        elif options.useCPU:
            if options.SinM.dtype != np.float32:
                options.SinM = options.SinM.astype(np.float32)
            libname = _nativeLibPath(libdir, 'CPU_matrixfree_lib')
        else:
            if options.SinM.dtype == 'uint16':
                libN = 'OpenCL_matrixfree_uint16_lib'
            elif options.SinM.dtype == 'uint8':
                libN = 'OpenCL_matrixfree_uint8_lib'
            else:
                libN = 'OpenCL_matrixfree_lib'
            libname = _nativeLibPath(libdir, libN)
        residualP = residual.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        if options.SinM.dtype == 'uint16':
            SinoP = options.SinM.ctypes.data_as(ctypes.POINTER(ctypes.c_uint16))
        elif options.SinM.dtype == 'uint8':
            SinoP = options.SinM.ctypes.data_as(ctypes.POINTER(ctypes.c_uint8))
        else:
            SinoP = options.SinM.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        outputP = output.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        FPOutputP = FPOutput.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        from omegatomo.util.dllpath import addDLLDirectories
        addDLLDirectories()
        c_lib = ctypes.CDLL(libname)
        status = c_lib.omegaMain(options.param, ctypes.c_char_p(inStr), SinoP, outputP, FPOutputP, residualP)
    finally:
        # Put back the user-visible options (nRowsD, x, SinM, offsetCorrection, ...) that the offset FDK preprocessing changed. This covers
        # everything from the preprocessing to the end of the library call, also when any of the steps in between raises.
        if offsetFDKState is not None:
            from .offsetfdk import offsetFDKRestore
            offsetFDKRestore(options, offsetFDKState)
    if status != 0:
        raise RuntimeError(f'Native reconstruction failed with status {status}; see the backend diagnostics above.')
    try:
        if options.useMultiResolutionVolumes and options.storeMultiResolution:
            output = _reshape_multiresolution_output(output, options)
        elif options.useMultiResolutionVolumes:
            output = _reshape_image_output(
                output,
                (options.NxOrig, options.NyOrig, options.NzOrig),
                _saved_image_count(options),
                options.Nt,
            )
        else:
            output = _reshape_image_output(
                output,
                (options.Nx[0], options.Ny[0], options.Nz[0]),
                _saved_image_count(options),
                options.Nt,
            )
        # Port of reconstructions_main.m:722-731: crop the extended FOV back off
        # the output. This only makes sense once `output` has actually been
        # reshaped above (the first three axes are spatial in every branch
        # except the useMultiResolutionVolumes one, which MATLAB also skips
        # here) -- if any of the reshapes above raised, we never reach this
        # point and the (still flat) output is returned unmodified, exactly as
        # MATLAB's crop would be meaningless on an unreshaped array.
        if options.useEFOV and not options.useMultiResolutionVolumes:
            if options.transaxialEFOV:
                # MATLAB crops both x and y with the same nTrans computed from Nx.
                nTrans = int((int(options.Nx[0]) - int(options.NxOrig)) // 2)
                if nTrans > 0:
                    output = output[nTrans:-nTrans, nTrans:-nTrans, ...]
            if options.axialEFOV:
                nAxial = int((int(options.Nz[0]) - int(options.NzOrig)) // 2)
                if nAxial > 0:
                    output = output[:, :, nAxial:-nAxial, ...]
        if options.subsets == 1 and options.storeFP:
            FPOutput = FPOutput.reshape((options.nRowsD, options.nColsD, options.nProjections, options.TOF_bins, options.Niter), order = 'F')
        elif options.storeFP and options.subsets > 1 and options.subsetType >= 8 and not options.FDK:
            # Port of reconstructions_main.m:732-741 for the subsets > 1 (subset
            # type >= 8) case. device_to_host() in reconstructionAF.h fills the
            # flat FPptr/FPOutput buffer iteration-major, subset-minor (for ii in
            # 0..subsets*Niter-1: kk = ii // subsets is the iteration, jj = ii %
            # subsets is the subset), each block being that subset's forward
            # projection with TOF bins packed into the trailing dimension. Split
            # it into that same Niter*subsets list of chunks and reshape each to
            # (nRowsD, nColsD, meas*TOF_bins), mirroring MATLAB's fp{uu} cell.
            fpList = []
            offset = 0
            nMeas = np.asarray(options.nMeas).reshape(-1)
            for _ in range(int(options.Niter)):
                for ii in range(int(options.subsets)):
                    meas = int(nMeas[ii + 1] - nMeas[ii])
                    block = int(options.nRowsD) * int(options.nColsD) * meas * int(options.TOF_bins)
                    fpList.append(FPOutput[offset : offset + block].reshape(
                        (options.nRowsD, options.nColsD, meas * options.TOF_bins), order='F'))
                    offset += block
            FPOutput = fpList
    except Exception as e:
        # Keep the reconstruction even if the output dimensions do not match
        warnings.warn(f'Could not reshape the reconstruction output ({e}); returning the unreshaped (flat) arrays instead.')
    try:
        from .imageProperties import imageProperties
        options.image_properties = imageProperties(options)
        callerOptions.image_properties = options.image_properties
    except Exception as e:
        warnings.warn(f'Could not create options.image_properties ({e}).')
    toc = time.perf_counter()
    if options.verbose > 0:
        print(f"Reconstruction took {toc - tic:0.4f} seconds")
    if options.storeResidual:
        return output, FPOutput, residual
    else:
        return output, FPOutput
