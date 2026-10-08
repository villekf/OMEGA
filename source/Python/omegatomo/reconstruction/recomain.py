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
import ctypes
import numpy as np
import warnings

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
    'T': 'B',
    'POCS': 'ASD_POCS',
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
_POINTER_SPECIAL_FIELDS = {'blurPlanes', 'blurPlanes2', 'gFilter', 'gFSize', 'gFilterODRT'}


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
        'eFOV': options.eFOVIndices.size,
        'sizeX': sizeXVal,
        'sizeZ': options.z.size,
        'sizeAtten': options.vaimennus.size,
        'sizeNorm': options.normalization.size,
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
        options.param.blurPlanes = options.blurPlanes[0].ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
        options.param.blurPlanes2 = options.blurPlanes2[0].ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
        # The native type-6 branch uses the volume-0 filter.
        options.param.gFilter = options.gFilter[0].ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        options.gFSize = np.array(options.gFilter[0].shape, dtype=np.uint64)
    else:
        options.param.blurPlanes = None
        options.param.blurPlanes2 = None
        options.param.gFilter = None
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


def _selectPrecorrectedMeasurement(options, getKey):
    """
    Selects either the corrected ('SinM') or the raw ('raw_SinM') measurement
    data key, matching the precorrection semantics shared by the MAT and NPZ
    measurement-file loaders. `getKey(key)` returns the array for that key,
    or raises KeyError.
    """
    if (options.randoms_correction or options.scatter_correction or options.normalization_correction) and not options.corrections_during_reconstruction:
        if not options.precorrect:
            try:
                return getKey('SinM')
            except KeyError:
                options.precorrect = True
                return getKey('raw_SinM')
        else:
            return getKey('raw_SinM')
    else:
        return getKey('raw_SinM')


def _nativeLibPath(libdir, name):
    """
    Returns the full path to the native reconstruction library `name` inside
    `libdir`, using the platform-appropriate shared-library extension
    ('.dll' on Windows, '.so' everywhere else).
    """
    import os
    ext = '.dll' if os.name == 'nt' else '.so'
    return str(os.path.join(libdir, name + ext))


def _loadDelayedMeasurement(options, getKey):
    """
    Loads the randoms ('SinDelayed') data from an already-open measurement
    file, matching the semantics shared by the MAT and NPZ measurement-file
    loaders. `getKey(key)` returns the array for that key, or raises
    KeyError.
    """
    if options.randoms_correction and not options.reconstruct_scatter and not options.reconstruct_trues and options.SinDelayed.size < 1:
        try:
            options.SinDelayed = getKey('SinDelayed')
        except KeyError:
            print('Randoms correction selected but no randoms data found. The randoms data should be saved as SinDelayed')


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
    """
    import time
    import os
    from .prepass import prepassPhase
    from .prepass import parseInputs
    from .prepass import loadCorrections
    tic = time.perf_counter()
    options.addProjector()
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
        if options.reconstruct_trues:
            options.SinM = matGet('SinTrues')
        elif options.reconstruct_scatter:
            options.SinM = matGet('SinScatter')
        else:
            options.SinM = _selectPrecorrectedMeasurement(options, matGet)
        _loadDelayedMeasurement(options, matGet)
    elif sinoSize < 1 and suffix == '.npy':
        options.SinM = np.load(options.fpath)
    elif sinoSize < 1 and suffix == '.npz':
        varList = np.load(options.fpath)
        npzGet = lambda key: varList[key]
        options.SinM = _selectPrecorrectedMeasurement(options, npzGet)
        _loadDelayedMeasurement(options, npzGet)
    elif not options.corrections_during_reconstruction and not options.precorrect and (options.randoms_correction or options.scatter_correction or options.normalization_correction):
        print('Corrections selected and measurement data found. The input measurement data WILL NOT BE PRECORRECTED!!!!!! If you wish to have OMEGA-based precorrection, make sure options.precorrect = True')
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
                try:
                    options.SinDelayed = np.array(var["SinDelayed"],order='F')
                except KeyError:
                    print('Randoms correction selected but no randoms data found. The randoms data should be saved as SinDelayed. Disabling randoms correction')
                    options.randoms_correction = False
            elif fsuffix == '.npy':
                options.SinDelayed = np.load(fpath)
            elif fsuffix == '.npz':
                varList = np.load(fpath)
                try:
                    options.SinDelayed = varList['SinDelayed']
                except KeyError:
                    print('Randoms correction selected but no randoms data found. The randoms data should be saved as SinDelayed. Disabling randoms correction')
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
    if not options.listmode:
        if isinstance(options.SinM, list):
            # Frames must be contiguous (frame 1, then frame 2, ...), like MATLAB's
            # cell2mat(cellfun(@(x) x(:), ...)); np.reshape on a list would stack the
            # frames along a new leading axis and interleave them in F order.
            options.SinM = np.concatenate([np.asarray(frame).ravel(order='F') for frame in options.SinM])
        options.SinM = np.reshape(options.SinM, (int(options.nRowsD), int(options.nColsD), options.nProjections, options.TOF_bins, options.Nt), order='F')
    elif options.listmode and options.compute_sensitivity_image and options.SPECT:
        from omegatomo.projector.detcoord import getCoordinatesSPECT
        x_sensitivity, z_sensitivity = getCoordinatesSPECT(options)
        options.uV = np.float32(np.asfortranarray(x_sensitivity))
        options.zSens = np.float32(np.asfortranarray(z_sensitivity))
        options.sensitivityViewWeights = _spect_listmode_sensitivity_weights(
            options, options.uV.size // 6
        )
    elif options.listmode and options.compute_sensitivity_image:
        if hasattr(options, 'xSens') and np.size(options.xSens) > 0 and hasattr(options, 'zSens') and np.size(options.zSens) > 0:
            options.uV = np.float32(np.asfortranarray(options.xSens))
            options.z = np.float32(np.asfortranarray(options.zSens))
            options.det_per_ring = options.uV.size // 2
            options.rings = options.z.size
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
        options.SinM[options.SinM < 0] = 0
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
    if status != 0:
        raise RuntimeError(f'Native reconstruction failed with status {status}; see the backend diagnostics above.')
    try:
        # Number of saved image volumes per timestep: options.saveNIter/save_iter
        # cause the native code to store one volume per requested iteration
        # (plus the initial estimate), in addition to the per-timestep volumes.
        if options.saveNIter.size > 0:
            numSaves = int(options.saveNIter.size) + 1
        elif options.save_iter:
            numSaves = int(options.Niter) + 1
        else:
            numSaves = 1
        if options.useMultiResolutionVolumes and not options.storeMultiResolution:
            output = output.reshape((options.NxOrig, options.NyOrig, options.NzOrig, -1), order = 'F')
        elif not options.storeMultiResolution and options.Nt == 1:
            output = output.reshape((options.Nx[0], options.Ny[0], options.Nz[0], -1), order = 'F')
        elif numSaves > 1:
            # Memory layout (fastest to slowest): spatial voxels, timestep, save index
            output = output.reshape((options.Nx[0], options.Ny[0], options.Nz[0], options.Nt, numSaves), order = 'F')
        else:
            output = output.reshape((options.Nx[0], options.Ny[0], options.Nz[0], options.Nt), order = 'F')
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
    toc = time.perf_counter()
    if options.verbose > 0:
        print(f"Reconstruction took {toc - tic:0.4f} seconds")
    if options.storeResidual:
        return output, FPOutput, residual
    else:
        return output, FPOutput
