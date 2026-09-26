# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 13:57:46 2024

Copyright (C) 2024-2025 Ville-Veikko Wettenhovi

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
import os
import numpy as np

def _load_array_file(path, mat_key_index='auto'):
    """
    Loads an array stored in a .npy, .npz, or .mat file.

    Parameters
    ----------
    path : str
        Full path to the .npy/.npz/.mat file.
    mat_key_index : 'auto', 'last', or int, optional
        Selects which variable to use when reading a .mat file (ignored for
        .npy/.npz files, which always use the first stored array):
            'auto' - Mirrors the historical attenuation-loading convention:
                     if the first key in the loaded dict is '__header__' the
                     4th key (index 3) is used, otherwise the 1st key
                     (index 0) is used.
            'last' - Uses the last variable in the loaded dict
                     (list(var)[-1]), as used by the various reference
                     image loaders.
            int    - Uses that literal index into the loaded dict's keys.
        The default is 'auto'.

    Returns
    -------
    NumPy array, or None if pymatreader is required but not installed.

    """
    ext = os.path.splitext(path)[1].lower()
    if ext == '.npy':
        return np.load(path, allow_pickle=True)
    elif ext == '.npz':
        apu = np.load(path, allow_pickle=True)
        variables = list(apu.keys())
        return apu[variables[0]]
    elif ext == '.mat':
        try:
            from pymatreader import read_mat
        except ModuleNotFoundError:
            print('pymatreader package not found! Mat-files cannot be loaded. You can install pymatreader package with "pip install pymatreader".')
            return None
        var = read_mat(path)
        keys = list(var)
        if mat_key_index == 'last':
            return np.array(var[keys[-1]])
        elif mat_key_index == 'auto':
            if keys[0] == '__header__':
                return np.array(var[keys[3]])
            else:
                return np.array(var[keys[0]])
        else:
            return np.array(var[keys[mat_key_index]])
    else:
        raise ValueError('Unsupported datatype!')

def _load_reference_image(options, value, resize=True, emptyMessage=None,
                           squareCheck=None, squareOrder='C',
                           squareErrorMsg='Reference image has to be square',
                           squareResizeGuardNz1=False, squareUsesItem=True,
                           elseCheckNdim3=True, elseCheckShape0=True,
                           elseResizeGuardNz1=False,
                           finalize=True, castBeforeAsfortran=False, doAsfortran=True):
    """
    Loads a reference/anatomical-weighting image (from a .npy/.npz/.mat path
    if `value` is a string) and, unless `resize` is False, reshapes a flat/
    column input into a cubic (koko_apu, koko_apu, Nz) volume and/or resizes
    a 3D input to match the current reconstruction size (options.Nx/Ny/Nz).

    This factors out the "load if str, maybe resize to (Nx,Ny,Nz), ravel F
    float32" pattern shared by TVPrepass, APLSPrepass, NLMPrepass, the RDP
    reference image and the IEM referenceImage loading in prepassPhase.
    Every one of those call sites has slightly different quirks (reshape
    order, whether Nx/Ny/Nz are read with `.item()`, whether the resize is
    guarded by Nz > 1, whether shape[0] is checked, whether there even is a
    flat/column-reshape branch, ...) that were written independently; these
    are preserved exactly via the keyword arguments below rather than
    unified, since some may be latent bugs.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.
    value : str or NumPy array
        The reference image, or a path to load it from.
    resize : bool, optional
        If False, only the string-loading step is performed (used by RDP,
        which never resizes its reference image). The default is True.
    emptyMessage : str or None, optional
        If not None and `value` is an empty string, raises
        ValueError(emptyMessage). If None (the IEM referenceImage site), an
        empty string is instead passed straight to _load_array_file, which
        raises its own "Unsupported datatype!" ValueError.
    squareCheck : None, 'shape1', or 'ndim1_or_shape1', optional
        Selects how a flat/column input is detected for the square-reshape
        branch (TV: 'shape1', APLS: 'ndim1_or_shape1'). None (the default)
        disables this branch entirely, as for NLM, RDP and IEM.
    finalize : bool, optional
        If True (default), performs the final conversion (optionally
        np.asfortranarray, then ravel('F').astype(float32)). If False, the
        possibly reshaped/resized array is returned as-is; used by
        TVPrepass, which still needs the 3D array for its own min/max
        normalization and the TVtype == 1 assembleS() call before raveling
        it itself.
    castBeforeAsfortran : bool, optional
        If True, casts to float32 before np.asfortranarray (TV, APLS). If
        False (the default), the cast happens only in the final
        ravel().astype() step (NLM, RDP, IEM).
    doAsfortran : bool, optional
        If False, skips the np.asfortranarray step entirely (RDP). The
        default is True.

    Returns
    -------
    NumPy array
        The loaded (and, depending on the parameters, reshaped/resized/
        finalized) reference image.

    """
    if isinstance(value, str):
        if emptyMessage is not None and len(value) == 0:
            raise ValueError(emptyMessage)
        value = _load_array_file(value, mat_key_index='last')
    if resize:
        if squareCheck == 'shape1':
            isFlat = value.shape[1] == 1
        elif squareCheck == 'ndim1_or_shape1':
            isFlat = value.ndim == 1 or value.shape[1] == 1
        else:
            isFlat = False
        if isFlat:
            Nz0 = options.Nz[0].item() if squareUsesItem else options.Nz[0]
            Nx0 = options.Nx[0].item() if squareUsesItem else options.Nx[0]
            Ny0 = options.Ny[0].item() if squareUsesItem else options.Ny[0]
            koko_apu = np.sqrt(np.size(value) / Nz0)
            if np.floor(koko_apu) != koko_apu:
                raise ValueError(squareErrorMsg)
            koko_apu = int(koko_apu)
            value = value.reshape((koko_apu, koko_apu, Nz0), order=squareOrder)
            if koko_apu != Nx0 or value.shape[2] != Nz0:
                if not squareResizeGuardNz1 or Nz0 > 1:
                    from skimage.transform import resize as resizeImage
                    print('Resizing reference image')
                    value = resizeImage(value, (Nx0, Ny0, Nz0))
        else:
            Nx0 = options.Nx[0].item()
            Ny0 = options.Ny[0].item()
            Nz0 = options.Nz[0].item()
            if elseCheckNdim3:
                mismatch = value.ndim == 3 and (
                    (elseCheckShape0 and value.shape[0] != Nx0) or
                    value.shape[1] != Ny0 or value.shape[2] != Nz0
                )
            else:
                mismatch = value.shape[1] != Ny0 or value.shape[2] != Nz0
            if mismatch and (not elseResizeGuardNz1 or Nz0 > 1):
                from skimage.transform import resize as resizeImage
                print('Resizing reference image')
                value = resizeImage(value, (Nx0, Ny0, Nz0))
    if finalize:
        if castBeforeAsfortran:
            value = value.astype(dtype=np.float32)
        if doAsfortran:
            value = np.asfortranarray(value)
        value = value.ravel('F').astype(dtype=np.float32)
    return value

def _compute_lambda_vals(niter, subsets, stochastic):
    """
    Vectorized form of the per-iteration relaxation parameter loop used for
    BSREM/RAMLA/MBSREM/MRAMLA/ROSEM(_MAP)/PKMA/SPS/SART/ASD_POCS/SAGA.
    """
    i = np.arange(niter, dtype=np.float64)
    if stochastic:
        return 1. / (0.4 / subsets * i + 1.)
    else:
        return 1. / (i / 20. + 1.)

def _pkma_relaxation_values(niter, subsets, rho, delta):
    """
    Vectorized form of the (kk, ll) double loop used for PKMA-style momentum
    coefficients (alpha_PKMA, alphaPrecond, thetaCP). The original loops
    walk oo = 0 .. niter*subsets - 1 in row-major (kk outer, ll inner) order,
    computing 1 + rho*oo / (oo + delta) at each oo; this returns that same
    sequence directly as a function of oo.
    """
    oo = np.arange(niter * subsets, dtype=np.float64)
    return 1. + (rho * oo) / (oo + delta)

def linearizeData(options):
    """
    This function linearizes the input measurement data.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Returns
    -------
    None.

    """
    options.SinM = np.log(options.flat / options.SinM.astype(dtype=np.float32))

def loadCorrections(options):
    """
    This function loads all the corrections related data. It can also perform
    some precorrection steps such as sinogram interpolation or arc correction.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Raises
    ------
    ValueError
        If files are not found.

    Returns
    -------
    None.

    """
    import os
    from omegatomo.util.matlabRound import matlabRound
    normalization_shape = np.asarray(options.normalization).shape
    options.normZ = int(normalization_shape[2]) if options.SPECT and len(normalization_shape) == 3 else 1
    normalization_indexed_stack = bool(options.SPECT and int(options.normZ) == int(options.nHeads))

    def expand_detector_stack(values, target_shape, timestep=0):
        """Expand a detector-head stack for corrections applied before projection kernels."""
        compact = np.asarray(values).reshape((int(options.nRowsD), int(options.nColsD), int(options.nHeads)), order='F')
        detector_vector = np.asarray(options.DetectorVector, dtype=np.uint32).reshape(-1)
        projection_counts = np.asarray(
            getattr(options, 'nProjectionsPerFrame', [options.nProjections]), dtype=np.int64
        ).reshape(-1)
        if projection_counts.size == int(options.Nt):
            offsets = np.concatenate(([0], np.cumsum(projection_counts, dtype=np.int64)))
            detector_vector = detector_vector[int(offsets[timestep]) : int(offsets[timestep + 1])]
        frame_stride = int(options.nRowsD) * int(options.nColsD)
        if len(target_shape) == 1:
            if int(target_shape[0]) % frame_stride != 0:
                raise ValueError('The detector-indexed normalization cannot be expanded to the correction data shape.')
            n_projections = int(target_shape[0]) // frame_stride
        else:
            n_projections = int(target_shape[2])
        if detector_vector.size < n_projections:
            raise ValueError('DetectorVector is shorter than the projection data used for normalization correction.')
        expanded = compact[:, :, detector_vector[:n_projections]]
        if len(target_shape) == 1:
            return expanded.ravel(order='F')
        while expanded.ndim < len(target_shape):
            expanded = np.expand_dims(expanded, axis=-1)
        return expanded

    if options.attenuation_correction == 1:
        if options.vaimennus.size == 0:
            if len(options.attenuation_datafile) > 0 and os.path.splitext(options.attenuation_datafile)[1].lower() == '.mhd':
                try:
                    from SimpleITK import ReadImage as loadMetaImage
                    from SimpleITK import GetArrayFromImage
                    metaImage = loadMetaImage(options.attenuation_datafile)
                    options.vaimennus = GetArrayFromImage(metaImage)
                    options.vaimennus = np.asfortranarray(np.transpose(options.vaimennus, (2, 1, 0)))
                    apu = np.array(list(metaImage.GetSpacing()))
                    if options.CT_attenuation:
                        if matlabRound(apu[0].item()*100.)/100. > matlabRound(options.FOVa_x[0].item() / (options.Nx[0].item())*100.)/100. or matlabRound(apu[0].item()*100)/100 < matlabRound(options.FOVa_x[0].item() / (options.Nx[0].item())*100.)/100.:
                            options.vaimennus = options.vaimennus * (apu[0].item() / (options.FOVa_x[0].item() / (options.Nx[0].item())))
                except ModuleNotFoundError:
                    print('SimpleITK package not found! MetaImages cannot be loaded. You can install SimpleITK package with "pip install SimpleITK".')
            elif len(options.attenuation_datafile) > 0 and os.path.splitext(options.attenuation_datafile)[1].lower() == '.mat':
                apu = _load_array_file(options.attenuation_datafile, mat_key_index='auto')
                if apu is not None:
                    options.vaimennus = apu.astype(np.float32)
            elif len(options.attenuation_datafile) > 0 and os.path.splitext(options.attenuation_datafile)[1].lower() in ('.npy', '.npz'):
                options.vaimennus = _load_array_file(options.attenuation_datafile)
            else:
                import tkinter as tk
                from tkinter.filedialog import askopenfilename
                root = tk.Tk()
                root.withdraw()
                nimi = askopenfilename(title='Select attenuation datafile',filetypes=(('MHD, NPY, NPZ and MAT files','*.mhd *.mat *.npy *.npz'),('All','*.*')))
                if len(nimi) == 0:
                    raise ValueError("No file selected!")
                nimiExt = os.path.splitext(nimi)[1].lower()
                if nimiExt == '.mhd':
                    try:
                        from SimpleITK import ReadImage as loadMetaImage
                        from SimpleITK import GetArrayFromImage
                        from SimpleITK import ReadImage as loadMetaImage
                        from SimpleITK import GetArrayFromImage
                        metaImage = loadMetaImage(nimi)
                        options.vaimennus = GetArrayFromImage(metaImage)
                        apu = np.array(list(metaImage.GetSpacing()))
                        if options.CT_attenuation:
                            if matlabRound(apu[0].item()*100.)/100. > matlabRound(options.FOVa_x[0].item() / (options.Nx[0].item())*100.)/100. or matlabRound(apu[0].item()*100)/100 < matlabRound(options.FOVa_x[0].item() / (options.Nx[0].item())*100.)/100.:
                                options.vaimennus = options.vaimennus * (apu[0].item() / (options.FOVa_x[0].item() / (options.Nx[0].item())))
                    except ModuleNotFoundError:
                        print('SimpleITK package not found! MetaImages cannot be loaded. You can install SimpleITK package with "pip install SimpleITK".')
                elif nimiExt == '.mat':
                    apu = _load_array_file(nimi, mat_key_index='auto')
                    if apu is not None:
                        options.vaimennus = apu.astype(np.float32)
                elif nimiExt in ('.npy', '.npz'):
                    options.vaimennus = _load_array_file(nimi)
                else:
                    raise ValueError('Unsupported datatype!')
        if options.CT_attenuation:
            if options.vaimennus.ndim == 1:
                size_mismatch = options.vaimennus.shape[0] != options.N[0]
            else:
                size_mismatch = (not options.vaimennus.shape[0] == options.Nx[0] or not options.vaimennus.shape[1] == options.Ny[0].item() or not options.vaimennus.shape[2] == options.Nz[0].item())
            if size_mismatch:
                if options.vaimennus.shape[0] != options.N[0]:
                    print('Error: Attenuation data is of different size than the reconstructed image. Attempting resize!')
                    if options.vaimennus.ndim == 1:
                        raise ValueError('The attenuation image should be a 3D volume in order for the resize to work properly!')
                    from scipy.ndimage import zoom
                    options.vaimennus = zoom(options.vaimennus, (options.Nx[0] / options.vaimennus.shape[0], options.Ny[0] / options.vaimennus.shape[1], options.Nz[0] / options.vaimennus.shape[2]))
                    if (not options.vaimennus.shape[0] == options.Nx[0] or not options.vaimennus.shape[1] == options.Ny[0].item() or not options.vaimennus.shape[2] == options.Nz[0].item()) and not options.vaimennus.size == options.N[0]:
                        raise ValueError('Error: Attenuation data is of different size than the reconstructed image. Automatic resize failed.')
            if options.rotateAttImage != 0:
                atn = np.reshape(options.vaimennus, (options.Nx[0].item(), options.Ny[0].item(), options.Nz[0].item()))
                atn = np.rot90(atn,options.rotateAttImage)
                options.vaimennus = atn
            if options.flipAttImageXY:
                atn = np.reshape(options.vaimennus, (options.Nx[0].item(), options.Ny[0].item(), options.Nz[0].item()))
                atn = np.fliplr(atn)
                options.vaimennus = atn
            if options.flipAttImageZ:
                atn = np.reshape(options.vaimennus, (options.Nx[0].item(), options.Ny[0].item(), options.Nz[0].item()))
                atn = np.flip(atn,2)
                options.vaimennus = atn
            if options.attIncm:
                options.vaimennus /= 10.
        options.vaimennus = np.asfortranarray(options.vaimennus)
        options.vaimennus = options.vaimennus.ravel('F').astype(dtype=np.float32)
    if options.normalization_correction:
        normalizationFromFile = False
        if options.normalization.size == 0:
            normalizationFromFile = True
            normdir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', '..', '..', '..', 'mat-files')) + "/" +  options.machine_name + '_normalization_' + str(options.Ndist) + 'x' + str(options.Nang) + '_span' + str(options.span) + '.mat'
            if os.path.exists(normdir):
                try:
                    from pymatreader import read_mat
                    var = read_mat(normdir)
                    options.normalization = np.array(var["normalization"],order='F')
                except ModuleNotFoundError:
                    print('pymatreader package not found! Mat-files cannot be loaded. You can install pymatreader package with "pip install pymatreader".')
            else:
                normdir = os.path.join(os.path.dirname(options.fpath), options.machine_name + '_normalization_' + str(options.Ndist) + 'x' + str(options.Nang) + '_span' + str(options.span) + '.mat')
                if os.path.exists(normdir):
                    try:
                        from pymatreader import read_mat
                        var = read_mat(normdir)
                        options.normalization = np.array(var["normalization"],order='F')
                    except ModuleNotFoundError:
                        print('pymatreader package not found! Mat-files cannot be loaded. You can install pymatreader package with "pip install pymatreader".')
                else:
                    import tkinter as tk
                    from tkinter.filedialog import askopenfilename
                    root = tk.Tk()
                    root.withdraw()
                    nimi = askopenfilename(title='Select normalization datafile',filetypes=(('NRM, NPY, NPZ and MAT files','*.nrm *.mat *.npy *.npz'),('All','*.*')))
                    if len(nimi) == 0:
                        raise ValueError("No file selected!")
                    nimiExt = os.path.splitext(nimi)[1].lower()
                    if nimiExt == '.nrm':
                        options.normalization = np.fromfile(nimi, dtype=np.float32)
                        if options.normalization.size != options.Ndist * options.Nang * options.TotSinos and not options.use_raw_data:
                            raise ValueError('Size mismatch between the current data and the normalization data file')
                    elif nimiExt == '.mat':
                        from pymatreader import read_mat
                        var = read_mat(nimi)
                        options.normalization = np.array(var["normalization"])
                    elif nimiExt in ('.npy', '.npz'):
                        options.normalization = _load_array_file(nimi)
                    else:
                        raise ValueError('Unsupported datatype!')
            normalization_shape = np.asarray(options.normalization).shape
            options.normZ = int(normalization_shape[2]) if options.SPECT and len(normalization_shape) == 3 else 1
            normalization_indexed_stack = bool(options.SPECT and int(options.normZ) == int(options.nHeads))
            options.normalization = 1. / options.normalization.ravel('F').astype(dtype=np.float32)
            if not options.use_raw_data and options.NSinos != options.TotSinos and not normalization_indexed_stack:
                options.normalization = options.normalization[0 : options.Ndist * options.Nang * options.NSinos]
        if normalization_indexed_stack and options.normalization.size != int(options.nRowsD) * int(options.nColsD) * int(options.nHeads):
            raise ValueError('Detector-indexed normalization must contain one detector image for each detector head.')
        if not options.corrections_during_reconstruction and options.precorrect:
            if normalization_indexed_stack:
                if isinstance(options.SinM, list):
                    options.SinM = [
                        frame.astype(np.float32) / expand_detector_stack(options.normalization, frame.shape, timestep)
                        for timestep, frame in enumerate(options.SinM)
                    ]
                    normalization_for_data = None
                else:
                    normalization_for_data = expand_detector_stack(options.normalization, options.SinM.shape)
            else:
                normalization_for_data = np.reshape(options.normalization, options.SinM.shape, order='F')
            if normalization_for_data is not None:
                # File-loaded normalization was inverted on load, so dividing multiplies
                # by the raw coefficients. A prefilled PET normalization is multiplied
                # directly (as in MATLAB); SPECT divides (as in MATLAB).
                if not options.SPECT and not normalizationFromFile:
                    options.SinM = options.SinM.astype(np.float32) * normalization_for_data
                else:
                    options.SinM = options.SinM.astype(np.float32) / normalization_for_data
            options.normalization_correction = False
        else:
            options.normalization = options.normalization.ravel('F').astype(dtype=np.float32)
    # Other SPECT corrections: DEW/TEW scatter estimation from one (DEW) or
    # two (TEW) energy windows in options.ScatterC. Mirrors MATLAB's
    # loadCorrections.m "Other SPECT corrections" block. options.ScatterC is
    # here a list of 1 (DEW) or 2 (TEW) windows; each window is either a
    # plain array (static data) or a list of per-timestep arrays (dynamic
    # data), mirroring MATLAB's ScatterC{w}/ScatterC{w}{timestep} cells.
    if (options.SPECT and options.scatter_correction and isinstance(options.ScatterC, list)
            and options.SinDelayed.size <= 1 and options.subtract_scatter):  # From 10.1371/journal.pone.0269542
        nWindows = len(options.ScatterC)
        if nWindows not in (1, 2):
            raise ValueError('options.ScatterC must contain either one (DEW) or two (TEW) scatter energy windows for SPECT scatter correction!')
        if nWindows == 2:
            kLower = np.diff(np.asarray(options.eWin, dtype=np.float64)) / np.diff(np.asarray(options.eWinL, dtype=np.float64))
            kUpper = np.diff(np.asarray(options.eWin, dtype=np.float64)) / np.diff(np.asarray(options.eWinU, dtype=np.float64))

        def _scatterWindow(w, timestep):
            data = options.ScatterC[w]
            if timestep is not None and isinstance(data, list):
                data = data[timestep]
            return np.squeeze(data)

        def _combinedScatterEstimate(timestep=None):
            if nWindows == 1:  # DEW
                k = 1.
                return k * _scatterWindow(0, timestep)
            else:  # TEW
                return 0.5 * (kLower * _scatterWindow(0, timestep) + kUpper * _scatterWindow(1, timestep))

        useSingle = (options.implementation == 2 or options.implementation == 3
                     or options.implementation == 5 or options.useSingles)
        if isinstance(options.SinM, list):  # Dynamic data (SinM is a list of size options.Nt)
            if not options.corrections_during_reconstruction:
                for timestep in range(options.Nt):
                    options.SinM[timestep] = options.SinM[timestep] - _combinedScatterEstimate(timestep)
                options.scatter_correction = False
            else:
                options.SinDelayed = np.stack(
                    [np.asfortranarray(_combinedScatterEstimate(timestep)) for timestep in range(options.Nt)], axis=-1)
                if useSingle:
                    options.SinDelayed = options.SinDelayed.astype(np.float32)
                options.scatter_correction = False
                # Matches MATLAB loadCorrections.m's final flag-normalization
                # (~line 1074): during-reconstruction randoms_correction is
                # only enabled when options.ordinaryPoisson is set.
                options.randoms_correction = options.ordinaryPoisson
        else:  # Static data (SinM is not a list)
            scatterEstimate = _combinedScatterEstimate()
            if not options.corrections_during_reconstruction:
                options.SinM = options.SinM - scatterEstimate
                options.scatter_correction = False
            else:
                options.SinDelayed = np.asfortranarray(scatterEstimate)
                if useSingle:
                    options.SinDelayed = options.SinDelayed.astype(np.float32)
                options.scatter_correction = False
                # Matches MATLAB loadCorrections.m's final flag-normalization
                # (~line 1074): during-reconstruction randoms_correction is
                # only enabled when options.ordinaryPoisson is set.
                options.randoms_correction = options.ordinaryPoisson
        # ScatterC has now been fully consumed into SinM/SinDelayed above, so
        # it is reset to an empty array to avoid the generic (plain-array)
        # scatter-correction handling below re-processing it.
        options.ScatterC = np.empty(0, dtype=np.float32)
    if options.scatter_correction and options.normalization_correction and options.normalize_scatter and options.corrections_during_reconstruction:
        if normalization_indexed_stack:
            if isinstance(options.ScatterC, list):
                options.ScatterC = [
                    frame / expand_detector_stack(options.normalization, frame.shape, timestep)
                    for timestep, frame in enumerate(options.ScatterC)
                ]
            else:
                options.ScatterC /= expand_detector_stack(options.normalization, options.ScatterC.shape)
        else:
            options.ScatterC /= options.normalization
    if options.randoms_correction and options.ordinaryPoisson and options.variance_reduction:
        from omegatomo.util.Randoms_variance_reduction import Randoms_variance_reduction
        options.SinDelayed = Randoms_variance_reduction(options.SinDelayed, options)
    if options.randoms_correction and options.ordinaryPoisson and options.randoms_smoothing:
        from omegatomo.util.smoothing import randoms_smoothing
        options.SinDelayed = randoms_smoothing(options.SinDelayed, options)
    if options.scatter_correction and options.ordinaryPoisson and options.scatter_smoothing:
        from omegatomo.util.smoothing import randoms_smoothing
        options.ScatterC = randoms_smoothing(options.ScatterC, options)
    if (options.randoms_correction or options.scatter_correction) and not options.ordinaryPoisson:
        if options.randoms_correction and options.SinDelayed.size > 0 and options.randoms_smoothing:
            from omegatomo.util.smoothing import randoms_smoothing
            options.SinDelayed = randoms_smoothing(options.SinDelayed, options)
        if options.randoms_correction and options.SinDelayed.size > 0 and options.variance_reduction:
            from omegatomo.util.Randoms_variance_reduction import Randoms_variance_reduction
            options.SinDelayed = Randoms_variance_reduction(options.SinDelayed, options)
        if options.scatter_correction and options.ScatterC.size > 0 and options.scatter_smoothing:
            from omegatomo.util.smoothing import randoms_smoothing
            options.ScatterC = randoms_smoothing(options.ScatterC, options)
        if options.precorrect:
            if options.scatter_correction and options.ScatterC.size > 0 and options.SinDelayed.size > 0 and options.randoms_correction and options.subtract_scatter:
                options.SinM = options.SinM.astype(np.float32) - options.SinDelayed.astype(np.float32) - options.ScatterC.astype(np.float32)
            elif options.scatter_correction and options.ScatterC.size > 0 and not options.randoms_correction and options.subtract_scatter:
                options.SinM = options.SinM.astype(np.float32) - options.ScatterC.astype(np.float32) 
            elif options.SinDelayed.size > 0 and options.randoms_correction:
                options.SinM = options.SinM.astype(np.float32) - options.SinDelayed.astype(np.float32) 
    if options.scatter_correction and options.ScatterC.size > 0 and not options.subtract_scatter:
        options.additionalCorrection = True
        options.corrVector = options.ScatterC
    elif options.scatter_correction and options.ScatterC.size > 0 and options.SinDelayed.size <= 1 and options.subtract_scatter and options.ordinaryPoisson:
        options.SinDelayed = np.asfortranarray(options.ScatterC.astype(np.float32))
        options.scatter_correction = False
        # Matches MATLAB loadCorrections.m's final flag-normalization (~line 1078): during-
        # reconstruction randoms_correction is only enabled when SinDelayed ends up with more
        # than one element (and ordinaryPoisson, already guaranteed true by this branch).
        options.randoms_correction = options.SinDelayed.size > 1
    elif options.scatter_correction and options.ScatterC.size > 0 and options.SinDelayed.size > 1 and options.subtract_scatter and options.ordinaryPoisson:
        options.SinDelayed += np.asfortranarray(options.ScatterC.astype(np.float32))
    if options.arc_correction:
        from omegatomo.util.arcCorrection import arc_correction
        x, y, options = arc_correction(options, True)
    if options.sampling > 1:
        if options.corrections_during_reconstruction:
            from omegatomo.util.sampling import interpolateSinog
            if options.normalization_correction and not normalization_indexed_stack:
                options.normalization = interpolateSinog(options.normalization, options.sampling, options.Ndist, options.Nang, options.sampling_interpolation_method)
            if options.randoms_correction:
                options.SinDelayed = interpolateSinog(options.SinDelayed, options.sampling, options.Ndist, options.Nang, options.sampling_interpolation_method)
            if options.additionalCorrection:
                options.corrVector = interpolateSinog(options.corrVector, options.sampling, options.Ndist, options.Nang, options.sampling_interpolation_method)
        from omegatomo.util.sampling import increaseSampling
        x, y, options = increaseSampling(options, None, None, True)
                
        
        
def parseInputs(options, mDataFound = False):
    """
    This function parses the input measurement data such that the elements
    are correctly divided between the subsets. Also does the same for
    corrections if they are applied during reconstruction.
    
    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.
    mDataFound : bool, optional
        If True, measurement data was found to be present already. The division
        into subsets is thus skipped if no data has been input before this 
        step. The default is False.

    Returns
    -------
    None.

    """
    normalization_indexed_stack = bool(options.SPECT and int(options.normZ) == int(options.nHeads))
    single_frame_list = (
        options.Nt <= 1 and isinstance(options.SinM, list)
    )
    if single_frame_list:
        options.SinM = options.SinM[0]
        if isinstance(options.index, list):
            options.index = options.index[0]

    if options.subsets > 1 and options.subsetType > 0:
        if mDataFound and not options.largeDim:
            if options.Nt > 1:
                for ff in range(1, options.Nt + 1):
                    if not options.use_raw_data:
                        if isinstance(options.SinM, list):
                            temp = options.SinM[ff - 1]
                        else:
                            if options.listmode == 0:
                                if options.TOF:
                                    temp = options.SinM[:,:,:,:,ff - 1]
                                else:
                                    if options.SinM.ndim == 4:
                                        temp = options.SinM[:,:,:,ff - 1]
                                    else:
                                        temp = np.squeeze(options.SinM[:,:,:,:,ff - 1])
                                if options.NSinos != options.TotSinos:
                                    temp = temp[:, :, :options.NSinos, :]
                            else:
                                temp = options.SinM[:,ff - 1]
                    # else:
                    #     temp = np.single(np.full(options.SinM[ff - 1]))
            
                    if options.TOF and options.listmode == 0:
                        if options.subsetType >= 8:
                            temp = temp[:, :, options.index, :]
                        else:
                            koko = (temp.shape[0], temp.shape[1], temp.shape[2], temp.shape[3])
                            temp = np.reshape(temp, (temp.size // options.TOF_bins, options.TOF_bins),order='F')
                            temp = temp[options.index, :]
                            temp = np.reshape(temp, koko, order='F')
                    else:
                        if options.subsetType >= 8:
                            if isinstance(options.index, list):
                                temp = temp[:, :, options.index[ff - 1]]
                            else:
                                temp = temp[:, :, options.index]
                        else:
                            if temp.ndim == 3:
                                koko = (temp.shape[0], temp.shape[1], temp.shape[2])
                            else:
                                koko = temp.shape[0]
                            temp = temp.ravel(order='F')
                            if isinstance(options.index, list):
                                temp = temp[options.index[ff - 1]]
                            else:
                                temp = temp[options.index]
                            temp = np.reshape(temp, koko, order='F')
                    if isinstance(options.SinM, list):
                        options.SinM[ff - 1] = temp
                    else:
                        if options.TOF and options.listmode == 0:
                            options.SinM[:,:,:,:,ff - 1] = temp.ravel(order='F')
                        else:
                            if options.SinM.ndim == 4:
                                options.SinM[:,:,:,ff - 1] = temp
                            else:
                                options.SinM[:,:,:,0,ff - 1] = temp
            else:
                if not options.use_raw_data and options.listmode == 0:
                    if options.NSinos != options.TotSinos:
                        options.SinM = options.SinM[:, :, :options.NSinos, :]
                # else:
                #     options.SinM = np.single(np.full(options.SinM[0]))
            
                if options.subsetType >= 8:
                    options.SinM = np.reshape(options.SinM, (options.nRowsD, options.nColsD, options.nProjections, options.TOF_bins), order='F')
            
                if options.TOF and options.listmode == 0:
                    if options.subsetType >= 8:
                        options.SinM = options.SinM[:, :, options.index, :]
                    else:
                        options.SinM = np.reshape(options.SinM, (options.SinM.size // options.TOF_bins, options.TOF_bins), order='F')
                        options.SinM = options.SinM[options.index, :]
                else:
                    if options.subsetType >= 8:
                        options.SinM = options.SinM[:, :,options.index]
                    elif options.subsetType > 0:
                        options.SinM = options.SinM.ravel(order='F')
                        options.SinM = options.SinM[options.index]
        if options.normalization_correction and options.corrections_during_reconstruction:
            if options.Nt > 1 and not normalization_indexed_stack:
                # Mirror the per-frame handling already used for SinM/SinDelayed/ScatterC
                # above: options.normalization is a flat, frame-major-concatenated array (one
                # Ndist*Nang*TotSinos block per frame -- see init.py's [Nt][subsets] device
                # buffer construction, which slices this same flat array via
                # options.nTotMeas[timestep*subsets+subset : ...+1]) that must be truncated to
                # NSinos (if reduced) and subset-selected PER FRAME, not once for the whole
                # array. Without this loop (pre-existing bug, not previously reachable/tested
                # since nothing built genuinely per-frame normalization before), a single
                # frame's worth of (subset-reordered) data would silently replace the entire
                # array, and every other frame's device-buffer slice ends up empty or reads
                # out of bounds (observed as a native OpenCL kernel crash).
                single_frame_len = options.Ndist * options.Nang * options.TotSinos
                out_frames = []
                for ff in range(1, options.Nt + 1):
                    frame = options.normalization[(ff - 1) * single_frame_len: ff * single_frame_len]
                    frame = np.reshape(frame, (options.Ndist, options.Nang, -1), order='F')
                    if not options.use_raw_data and options.NSinos != options.TotSinos:
                        frame = frame[:, :, :options.NSinos]
                    idx = options.index[ff - 1] if isinstance(options.index, list) else options.index
                    if options.subsetType >= 8:
                        frame = frame[:, :, idx]
                    else:
                        frame = frame.ravel(order='F')[idx]
                    out_frames.append(np.asarray(frame).ravel(order='F'))
                options.normalization = np.concatenate(out_frames).astype(dtype=np.float32)
            else:
                if not normalization_indexed_stack and not options.use_raw_data and options.NSinos != options.TotSinos:
                    options.normalization = options.normalization[:options.NSinos * options.Ndist * options.Nang]
                if normalization_indexed_stack:
                    options.normalization = options.normalization.ravel(order='F').astype(dtype=np.float32)
                elif options.subsetType >= 8:
                    options.normalization = np.reshape(options.normalization, (options.Ndist, options.Nang, -1),order='F')
                    options.normalization = options.normalization[:, :, options.index]
                    options.normalization = options.normalization.ravel(order='F').astype(dtype=np.float32)
                else:
                    options.normalization = options.normalization[options.index]
        
        if options.additionalCorrection and hasattr(options, 'corrVector') and options.corrVector.size > 0:
            if options.subsetType >= 8:
                options.corrVector = np.reshape(options.corrVector, (options.Ndist, options.Nang, options.nProjections, -1), order='F')
                options.corrVector = options.corrVector[:, :, options.index, :]
                options.corrVector = options.corrVector.ravel(order='F')
            else:
                options.corrVector = np.reshape(options.corrVector, (options.Ndist * options.Nang * options.nProjections, -1), order='F')
                options.corrVector = options.corrVector[options.index,:]
                options.corrVector = options.corrVector.ravel(order='F').astype(dtype=np.float32)
        
        if (options.randoms_correction
                and options.corrections_during_reconstruction 
                and not options.reconstruct_trues and not options.reconstruct_scatter) and not options.largeDim:
            
            if options.SinDelayed.size > 1:
                if options.Nt > 1:
                    for ff in range(1, options.Nt + 1):
                        if not options.use_raw_data:
                            temp = options.SinDelayed[:,:,:,ff - 1]
                            if options.NSinos != options.TotSinos:
                                temp = temp[:, :, :options.NSinos]
                        # else:
                        #     temp = np.single(np.full(options.SinDelayed[ff - 1]))
            
                        if options.subsetType >= 8:
                            temp = np.reshape(temp, (options.Ndist, options.Nang, -1),order='F')
                            temp = temp[:, :, options.index]
                            temp = temp.ravel(order='F')
                        else:
                            temp = temp.ravel(order='F')
                            temp = temp[options.index]
                        options.SinDelayed[:,:,:,ff - 1] = np.reshape(temp, (options.Ndist, options.Nang, -1),order='F')
                    options.SinDelayed = options.SinDelayed.ravel(order='F').astype(dtype=np.float32)
            
                else:
                    # if isinstance(options.SinDelayed, list):
                    #     if options.subsetType >= 8:
                    #         options.SinDelayed[0] = np.reshape(options.SinDelayed[0], (options.Ndist, options.Nang, -1))
                    #         options.SinDelayed[0] = options.SinDelayed[0][:, :, options.index]
                    #         options.SinDelayed[0] = options.SinDelayed[0].ravel()
                    #     else:
                    #         options.SinDelayed[0] = options.SinDelayed[0][options.index]
                    # else:
                    if options.subsetType >= 8:
                        options.SinDelayed = np.reshape(options.SinDelayed, (options.Ndist, options.Nang, -1))
                        options.SinDelayed = options.SinDelayed[:, :, options.index]
                        options.SinDelayed = options.SinDelayed.ravel(order='F').astype(dtype=np.float32)
                    else:
                        options.SinDelayed = options.SinDelayed.ravel(order='F').astype(dtype=np.float32)
                        options.SinDelayed = options.SinDelayed[options.index]
        
        if (options.scatter_correction and options.corrections_during_reconstruction 
                and not options.reconstruct_trues and not options.reconstruct_scatter):
            if not options.largeDim:
                if options.Nt > 1: #and isinstance(options.ScatterC, list) and len(options.ScatterC) > 1:
                    for ff in range(1, options.Nt + 1):
                        if not options.use_raw_data:
                            temp = options.ScatterC[:,:,:,:,ff - 1]
                            if options.NSinos != options.TotSinos:
                                temp = temp[:, :, :options.NSinos,:]
                        # else:
                        #     temp = np.single(np.full(options.ScatterC[ff - 1]))
            
                        if options.subsetType >= 8:
                            temp = temp[:, :, options.index,:]
                        else:
                            temp = temp.ravel(order='F')
                            temp = temp[options.index]
                            temp = np.reshape(temp, (options.nRowsD, options.nColsD, options.nProjections, -1), order='F')
                        options.ScatterC[:,:,:,:,ff - 1] = temp
                    options.ScatterC = options.ScatterC.ravel(order='F').astype(dtype=np.float32)
            
                else:
                    # if isinstance(options.ScatterC, list):
                    #     if options.subsetType >= 8:
                    #         options.ScatterC[0] = np.reshape(options.ScatterC[0], (options.Ndist, options.Nang, -1))
                    #         options.ScatterC[0] = options.ScatterC[0][:, :, options.index]
                    #         options.ScatterC[0] = options.ScatterC[0].ravel()
                    #     else:
                    #         options.ScatterC = options.ScatterC[0][options.index]
                    # else:
                    if options.subsetType >= 8:
                        options.ScatterC = np.reshape(options.ScatterC, (options.Ndist, options.Nang, options.nProjections, -1), order='F')
                        options.ScatterC = options.ScatterC[:, :, options.index,:]
                        options.ScatterC = options.ScatterC.ravel(order='F').astype(dtype=np.float32)
                    else:
                        options.ScatterC = options.ScatterC.ravel(order='F').astype(dtype=np.float32)
                        options.ScatterC = options.ScatterC[options.index]
                if options.randoms_correction == 1 and options.SinDelayed.size == options.ScatterC.size:
                    options.SinDelayed = options.SinDelayed + options.ScatterC
                else:
                    options.SinDelayed = options.ScatterC
            else:
                if options.randoms_correction == 1 and options.SinDelayed.size == options.ScatterC.size:
                    options.SinDelayed = options.SinDelayed + options.ScatterC
                else:
                    options.SinDelayed = options.ScatterC
                
        
        if options.attenuation_correction and not options.CT_attenuation:
            if options.Nt > 1:
                # Measurement-domain attenuation is per-timestep in C++/MATLAB for dynamic
                # data (see init.py's [Nt][subsets] d_atten construction, sliced from this
                # same flat array via options.nTotMeas); mirror the per-frame handling above
                # (normalization_correction, SinM/SinDelayed/ScatterC) instead of subset-
                # selecting the whole Nt-frame-concatenated array with a single frame's worth
                # of indices (pre-existing bug -- see the normalization_correction branch
                # above for the full explanation of the failure mode this avoids).
                single_frame_len = options.Ndist * options.Nang * options.TotSinos
                out_frames = []
                for ff in range(1, options.Nt + 1):
                    frame = options.vaimennus[(ff - 1) * single_frame_len: ff * single_frame_len]
                    idx = options.index[ff - 1] if isinstance(options.index, list) else options.index
                    if options.subsetType >= 8:
                        frame = np.reshape(frame, (options.Ndist, options.Nang, -1), order='F')
                        frame = frame[:, :, idx]
                    else:
                        frame = frame[idx]
                    out_frames.append(np.asarray(frame).ravel(order='F'))
                options.vaimennus = np.concatenate(out_frames)
            elif options.subsetType >= 8:
                options.vaimennus = np.reshape(options.vaimennus, (options.Ndist, options.Nang, -1), order='F')
                options.vaimennus = options.vaimennus[:, :, options.index]
                options.vaimennus = options.vaimennus.ravel(order='F')
            else:
                options.vaimennus = options.vaimennus.ravel(order='F')
                options.vaimennus = options.vaimennus[options.index]
        if options.scatter_correction and not options.corrections_during_reconstruction:
            options.scatter_correction = False
        if options.randoms_correction and not options.corrections_during_reconstruction:
            options.randoms_correction = False
        if options.useMaskFP and options.maskFPZ > 1 and options.maskFPZ != options.nHeads and options.subsetType >= 8:
            options.maskFP = options.maskFP[:,:,options.index]
            
    
    if options.Nt <= 1 and mDataFound and not options.largeDim and options.loadTOF:
        if single_frame_list:
            options.SinM = [np.asfortranarray(options.SinM).astype(dtype=np.float32)]
            options.index = [options.index]
        else:
            options.SinM = np.asfortranarray(options.SinM)
            options.SinM = options.SinM.ravel(order='F').astype(dtype=np.float32)
    if single_frame_list and not isinstance(options.index, list):
        options.index = [options.index]


def TVPrepass(options):
    """
    Performs some possible TV-related prepass computations. These include
    loading of the anatomical weighting data and weighting values for TV type 1

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Raises
    ------
    ValueError
        If files are not found or the size is different.

    Returns
    -------
    S : NumPy array
        Weighting coefficients for TV type 1 when using anatomical weighting.

    """
    def assembleS(alkuarvo,T,Ny,Nx,Nz):
        f = -np.diff(alkuarvo, axis=1)
        f = np.concatenate((f, np.zeros((Nx, 1, Nz),order='F',dtype=np.float32)), axis=1)
        f = f.ravel('F')
        g = -np.diff(alkuarvo, axis=0)
        g = np.concatenate((g, np.zeros((1, Ny, Nz),order='F',dtype=np.float32)), axis=0)
        g = g.ravel('F')
        h = -np.diff(alkuarvo, axis=2)
        h = np.concatenate((h, np.zeros((Nx, Ny, 1),order='F',dtype=np.float32)), axis=2)
        h = h.ravel('F')

        gradvec = np.vstack((f, g, h))

        gradnorm = np.linalg.norm(gradvec, axis=0)

        gamma = np.exp(-gradnorm ** 2 / (T ** 2))

        # Construct the matrix S. Vectorized form of the per-voxel loop:
        #   B = I - (1 - gamma) * outer(nu, nu) where gradnorm > 0, else I
        #   S[3*ll:3*ll+3, :] = B
        L = np.size(gradnorm)
        eye3 = np.eye(3)
        B = np.broadcast_to(eye3, (L, 3, 3)).copy()
        mask = gradnorm > 0
        if np.any(mask):
            nu = gradvec[:, mask] / gradnorm[mask]
            outer = np.einsum('il,jl->lij', nu, nu)
            B[mask] = eye3 - (1. - gamma[mask])[:, None, None] * outer
        S = np.zeros((Nx * Ny * Nz * 3, 3),order='F',dtype=np.float32)
        S[:, :] = B.reshape(L * 3, 3)
        return S
    if options.TV_use_anatomical:
        options.TV_referenceImage = _load_reference_image(
            options, options.TV_referenceImage,
            emptyMessage='TV with anatomical weighting selected, but no reference image provided!',
            squareCheck='shape1', squareOrder='F',
            squareErrorMsg='Reference image has to be square',
            squareResizeGuardNz1=True, squareUsesItem=True,
            elseCheckNdim3=False, elseResizeGuardNz1=True,
            finalize=False,
        )
        options.TV_referenceImage = options.TV_referenceImage.astype(dtype=np.float32)
        options.TV_referenceImage = options.TV_referenceImage - np.min(options.TV_referenceImage)
        options.TV_referenceImage = options.TV_referenceImage / np.max(options.TV_referenceImage)
        if options.TVtype == 1:
            options.TV_referenceImage = options.TV_referenceImage.reshape((options.Nx[0].item(), options.Ny[0].item(), options.Nz[0].item()),order='F')
            S = assembleS(options.TV_referenceImage, options.B, options.Ny[0].item(), options.Nx[0].item(), options.Nz[0].item())
            S = S.astype(dtype=np.float32)
            s1 = S[0::3, 0]
            s2 = S[0::3, 1]
            s3 = S[0::3, 2]
            s4 = S[1::3, 0]
            s5 = S[1::3, 1]
            s6 = S[1::3, 2]
            s7 = S[2::3, 0]
            s8 = S[2::3, 1]
            s9 = S[2::3, 2]
            options.s = np.asfortranarray(np.concatenate((s1.flatten(), s2.flatten(), s3.flatten(), s4.flatten(), s5.flatten(), s6.flatten(), s7.flatten(), s8.flatten(), s9.flatten())))
        options.TV_referenceImage = np.asfortranarray(options.TV_referenceImage)
        options.TV_referenceImage = options.TV_referenceImage.ravel('F').astype(dtype=np.float32)
        
def APLSPrepass(options):
    """
    Loads the anatomical reference image for APLS.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Raises
    ------
    ValueError
        If files are not found or the size is different.

    Returns
    -------
    None.

    """
    options.APLS_ref_image = _load_reference_image(
        options, options.APLS_ref_image,
        emptyMessage='APLS selected, but no reference image provided!',
        squareCheck='ndim1_or_shape1', squareOrder='F',
        squareErrorMsg='Reference image has to be 2D/3D if different size than reconstruction size!',
        squareResizeGuardNz1=False, squareUsesItem=False,
        elseCheckNdim3=True, elseCheckShape0=True, elseResizeGuardNz1=False,
        finalize=True, castBeforeAsfortran=True, doAsfortran=True,
    )

def computeWeights(options, GGMRF):
    """
    Computes distance-based weights for various regularization methods. A 
    special case for GGMRF is included. The weighting is based on the distance
    of the "center" voxel to the other voxels. The number of weights depends on
    the size of the neighborhood.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.
    GGMRF : bool
        If True, GGMRF is selected.

    Returns
    -------
    None.

    """
    distX = options.FOVa_x[0] / options.Nx[0]
    distY = options.FOVa_y[0] / options.Ny[0]
    distZ = options.axial_fov[0] / options.Nz[0]

    if np.size(options.weights) == 0:
        # Offset vectors, each running from +N down to -N (matches the
        # element order the original loop-based implementation produced).
        xr = np.arange(options.Ndx, -options.Ndx-1, -1) * distX
        yr = np.arange(options.Ndy, -options.Ndy-1, -1) * distY
        zr = np.arange(options.Ndz, -options.Ndz-1, -1) * distZ
        if GGMRF:
            # GGMRF-style ordering: z varies fastest, then y, then x slowest.
            Zg, Yg, Xg = np.meshgrid(zr, yr, xr, indexing='ij')
            if options.Ndx == 0 or options.Nx[0] == 1:
                dist = np.sqrt(Zg**2 + Yg**2)
            else:
                dist = np.sqrt(Zg**2 + Yg**2 + Xg**2)
        else:
            # Default ordering: x varies fastest, then y, then z slowest.
            Xg, Yg, Zg = np.meshgrid(xr, yr, zr, indexing='ij')
            if options.Ndz == 0 or options.Nz[0] == 1:
                dist = np.sqrt(Xg**2 + Yg**2)
            else:
                dist = np.sqrt(Xg**2 + Yg**2 + Zg**2)
        options.weights = 1.0 / dist.flatten(order='F')
        options.weights = options.weights.astype(dtype=np.float32)
        
def quadWeights(options, isEmpty):
    """
    Normalizes the weights. If the weights are manually input, no normalization
    is performed.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.
    isEmpty : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    if isEmpty:
        non_inf_weights_sum = np.sum(options.weights[~np.isinf(options.weights)])
        options.weights_quad = options.weights / non_inf_weights_sum
        if not options.GGMRF:
            half_len = np.size(options.weights_quad) // 2
            options.weights_quad = np.concatenate((options.weights_quad[:half_len], options.weights_quad[half_len + 1:]))
    else:
        options.weights_quad = options.weights
    # if not options.GGMRF:
    options.weights_quad = options.weights_quad[~np.isinf(options.weights_quad)]
    options.weights_quad = options.weights_quad.astype(dtype=np.float32)
        
def huberWeights(options):
    """
    Normalizes Huber prior weights.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Returns
    -------
    None.

    """
    if np.size(options.weights_huber) == 0:
        non_inf_weights_sum = np.sum(options.weights[np.isfinite(options.weights)])
        options.weights_huber = options.weights / non_inf_weights_sum
        half_len = np.size(options.weights_huber) // 2
        options.weights_huber = np.concatenate((options.weights_huber[:half_len], options.weights_huber[half_len + 1:]))
    options.weights_huber = options.weights_huber[~np.isinf(options.weights_huber)]
    options.weights_huber = options.weights_huber.astype(dtype=np.float32)
    
def weightedWeights(options):
    """
    Special weighting for weighted mean.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Returns
    -------
    None.

    """
    if np.size(options.weighted_weights) == 0:
        distX = options.FOVa_x / float(options.Nx[0].item())
        kerroin = np.sqrt(2.) * distX
        options.weighted_weights = kerroin * options.weights
        options.weighted_weights[np.isinf(options.weighted_weights)] = options.weighted_center_weight
        options.weighted_weights /= np.sum(options.weighted_weights)
    options.weighted_weights = np.reshape(options.weighted_weights, (options.Ndx * 2 + 1, options.Ndy * 2 + 1, options.Ndz * 2 + 1),order='F').astype(dtype=np.float32)

def NLMPrepass(options):
    """
    Computes the Gaussian weights for NLM and loads the reference image if
    selected.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Raises
    ------
    ValueError
        If files are not found.

    Returns
    -------
    gaussK : NumPy array
        Gaussian weights.

    """
    def gaussianKernel(x, y, z, sigma_x, sigma_y, sigma_z = 0):
        gaussK = np.exp(-(np.add.outer(np.add.outer(x**2 / (2*sigma_x**2), y**2 / (2*sigma_y**2)), z**2 / (2*sigma_z**2))))
        return gaussK
    g_x = np.linspace(-options.Nlx, options.Nlx, 2 * options.Nlx + 1, dtype=np.float32)
    g_y = np.linspace(-options.Nly, options.Nly, 2 * options.Nly + 1, dtype=np.float32)
    g_z = np.linspace(-options.Nlz, options.Nlz, 2 * options.Nlz + 1, dtype=np.float32)
    gaussian = gaussianKernel(g_x, g_y, g_z, options.NLM_gauss, options.NLM_gauss, options.NLM_gauss)
    options.gaussianNLM = gaussian.flatten('F').astype(dtype=np.float32)
    if options.NLM_use_anatomical:
        options.NLM_referenceImage = _load_reference_image(
            options, options.NLM_referenceImage,
            emptyMessage='NLM with anatomical weighting selected, but no reference image provided!',
            squareCheck=None,
            elseCheckNdim3=True, elseCheckShape0=True, elseResizeGuardNz1=False,
            finalize=True, castBeforeAsfortran=False, doAsfortran=True,
        )

def prepassPhase(options):
    """
    Computes various preprocessing phases, such as computing weights, loading 
    reference images, computing relaxation parameters, and making sure that
    many of the input variables are correctly formatted.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Raises
    ------
    ValueError
        Incorrect input variables or files are not found.

    Returns
    -------
    None.

    """
    from .rampfilt import rampFilt
    from omegatomo.util.matlabRound import matlabRound
    options.Nf = options.nRowsD
    if not isinstance(options.tauCP, np.ndarray):
        options.tauCP = np.array(options.tauCP, dtype=np.float32, ndmin=1)
    if not isinstance(options.sigmaCP, np.ndarray):
        options.sigmaCP = np.array(options.sigmaCP, dtype=np.float32, ndmin=1)
    if not isinstance(options.sigma2CP, np.ndarray):
        options.sigma2CP = np.array(options.sigma2CP, dtype=np.float32, ndmin=1)
    if not isinstance(options.tauCPFilt, np.ndarray):
        options.tauCPFilt = np.array(options.tauCPFilt, dtype=np.float32, ndmin=1)
    if not isinstance(options.thetaCP, np.ndarray):
        options.thetaCP = np.array(options.thetaCP, dtype=np.float32, ndmin=1)
    if not isinstance(options.alpha_PKMA, np.ndarray):
        options.alpha_PKMA = np.array(options.alpha_PKMA, dtype=np.float32, ndmin=1)
    if options.precondTypeImage[2]:
        options.referenceImage = _load_reference_image(
            options, options.referenceImage,
            emptyMessage=None,
            squareCheck=None,
            elseCheckNdim3=True, elseCheckShape0=True, elseResizeGuardNz1=False,
            finalize=True, castBeforeAsfortran=False, doAsfortran=True,
        )
        if np.size(options.referenceImage) == int(matlabRound((options.NxFull - options.NxOrig) * options.multiResolutionScale)) * \
            int(matlabRound((options.NyFull - options.NyOrig) * options.multiResolutionScale)) * \
            int(matlabRound((options.NzFull - options.NzOrig) * options.multiResolutionScale)):
            skip = True
        else:
            skip = False
        
        if not skip and np.size(options.referenceImage) != options.NxFull * options.NyFull * options.NzFull:
            raise ValueError('The size of the reference image does not match the reconstructed image!')
        
        if not skip and options.nMultiVolumes > 0:
            options.referenceImage = options.referenceImage.reshape(options.NxFull, options.NyFull, options.NzFull, order = 'F')
            from scipy.ndimage import zoom
            apu = zoom(options.referenceImage, options.multiResolutionScale, order=1).astype(np.float32)

        if not skip:
            if options.nMultiVolumes == 6:
                options.referenceImage = options.referenceImage[
                    (options.referenceImage.shape[0] - options.NxOrig) // 2 :
                    (options.referenceImage.shape[0] - options.NxOrig) // 2 + options.NxOrig,
                    (options.referenceImage.shape[1] - options.NyOrig) // 2 :
                    (options.referenceImage.shape[1] - options.NyOrig) // 2 + options.NyOrig,
                    (options.referenceImage.shape[2] - options.NzOrig) // 2 :
                    (options.referenceImage.shape[2] - options.NzOrig) // 2 + options.NzOrig
                ]
                apu1 = apu[options.Nx[3].item() : options.Nx[3].item() + options.Nx[1].item(), options.Ny[5].item() : options.Ny[5].item() + options.Ny[1].item(),  : options.Nz[1].item()]
                options.referenceImage = np.concatenate((options.referenceImage.ravel(), apu1.ravel()))
                apu1 = apu[options.Nx[4].item() : options.Nx[4].item() + options.Nx[2].item(),options.Ny[6].item() : options.Ny[6].item() + options.Ny[2].item(),-options.Nz[1].item() : ]
                options.referenceImage = np.concatenate((options.referenceImage.ravel(), apu1.ravel()))
                apu1 = apu[ : options.Nx[3].item(),:,:]
                options.referenceImage = np.concatenate((options.referenceImage.ravel(), apu1.ravel()))
                apu1 = apu[options.Nx[4].item() + options.Nx[2].item() :, :, : ]
                options.referenceImage = np.concatenate((options.referenceImage.ravel(), apu1.ravel()))
                apu1 = apu[options.Nx[3].item() : options.Nx[3].item() + options.Nx[1].item(), : options.Ny[5].item(), : options.Nz[3].item() ]
                options.referenceImage = np.concatenate((options.referenceImage.ravel(), apu1.ravel()))
                apu1 = apu[ options.Nx[4].item() : options.Nx[4].item() + options.Nx[2].item(), options.Ny[6].item() + options.Ny[2].item() :, : options.Nz[4].item() ]
                options.referenceImage = np.concatenate((options.referenceImage.ravel(), apu1.ravel()))
            elif options.nMultiVolumes == 4:
                options.referenceImage = options.referenceImage[
                    (options.referenceImage.shape[0] - options.NxOrig) // 2 :
                    (options.referenceImage.shape[0] - options.NxOrig) // 2 + options.NxOrig,
                    (options.referenceImage.shape[1] - options.NyOrig) // 2 :
                    (options.referenceImage.shape[1] - options.NyOrig) // 2 + options.NyOrig,
                    : ]
                apu1 = apu[ : options.Nx[1].item(), :, :]
                options.referenceImage = np.concatenate((options.referenceImage.ravel(), apu1.ravel()))
                apu1 = apu[ options.Nx[1].item() + options.Nx[0].item() :, :,:                ]
                options.referenceImage = np.concatenate((options.referenceImage.ravel(), apu1.ravel()))
                apu1 = apu[ options.Nx[1].item() : options.Nx[1].item() + options.Nx[3].item(), : options.Ny[3].item(), : ]
                options.referenceImage = np.concatenate((options.referenceImage.ravel(), apu1.ravel()))
                apu1 = apu[ options.Nx[1].item() : options.Nx[1].item() + options.Nx[4].item(), options.Ny[3].item() + options.Ny[0].item() :, : ]
                options.referenceImage = np.concatenate((options.referenceImage.ravel(), apu1.ravel()))
            elif options.nMultiVolumes == 2:
                options.referenceImage = options.referenceImage[ :, :, (options.referenceImage.shape[2] - options.NzOrig) // 2 : (options.referenceImage.shape[2] - options.NzOrig) // 2 + options.NzOrig ]
                apu1 = apu[:, :,  : options.Nz[1].item()]
                options.referenceImage = np.concatenate((options.referenceImage.ravel(), apu1.ravel()))
                apu1 = apu[:, :, -options.Nz[2].item():]
                options.referenceImage = np.concatenate((options.referenceImage.ravel(), apu1.ravel()))
        
    # Check if any of the regularization options are selected
    if (options.MRP or options.quad or options.Huber or options.TV or options.FMH or options.L or options.weighted_mean or options.APLS or options.BSREM
        or options.RAMLA or options.MBSREM or options.MRAMLA or options.ROSEM or options.DRAMA or options.ROSEM_MAP or options.ECOSEM or options.SART or options.ASD_POCS 
        or options.COSEM or options.ACOSEM or options.AD or np.any(options.OSL_COSEM) or options.NLM or options.OSL_RBI or options.RBI or options.PKMA or options.SAGA
        or options.RDP or options.SPS or options.ProxNLM or options.GGMRF):
    
        # Compute and/or load necessary variables for the TV regularization
        if options.TV and options.MAP:
            TVPrepass(options)
    
        # Load necessary variables for the APLS regularization
        if options.APLS and options.MAP:
            APLSPrepass(options)
    
        if options.U == 0:
            if options.CT:
                options.U = 10.
            else:
                options.U = 10000.
    
        # Lambda values (relaxation parameters)
        if (options.BSREM or options.RAMLA or options.MBSREM or options.MRAMLA or options.ROSEM_MAP or options.ROSEM or options.PKMA or options.SPS or options.SART or options.ASD_POCS or options.SAGA) and (np.size(options.lambdaN) == 0 or np.sum(options.lambdaN) == 0.):
            options.lambdaN = _compute_lambda_vals(options.Niter, options.subsets, options.stochasticSubsetSelection).astype(np.float32)
            if options.CT and not options.SART and not options.ASD_POCS:
                options.lambdaN = options.lambdaN / 10000.
        elif (options.BSREM or options.RAMLA or options.MBSREM or options.MRAMLA or options.ROSEM_MAP or options.ROSEM or options.PKMA or options.SPS or options.SART or options.ASD_POCS or options.SAGA):
            if np.size(options.lambdaN) < options.Niter:
                print('Warning: The number of relaxation values must be at least the number of iterations times the number of subsets! Computing custom relaxation values.')
                options.lambdaN = _compute_lambda_vals(options.Niter, options.subsets, options.stochasticSubsetSelection).astype(np.float32)
                if options.CT and not options.SART and not options.ASD_POCS:
                    options.lambdaN = options.lambdaN / 10000.
            elif np.size(options.lambdaN) > options.Niter:
                print('Warning: The number of relaxation values is more than the number of iterations. Later values are ignored!')
    
        if options.DRAMA:
            # r(i, j) = i * subsets + j + 1 (the loop's running counter starting at 1);
            # the pre-loop lam_drama[0, 0] assignment is always overwritten by the
            # i = j = 0 iteration below (r = 1 there too), so it is folded in directly.
            options.lam_drama = np.zeros((options.Niter, options.subsets),order='F',dtype=np.float32)
            r_vals = np.arange(1, options.Niter * options.subsets + 1, dtype=np.float64).reshape((options.Niter, options.subsets))
            options.lam_drama[:, :] = options.beta_drama / (options.alpha_drama * options.beta0_drama + r_vals)

        if options.PKMA and (np.size(options.alpha_PKMA) < options.Niter * options.subsets or np.sum(options.alpha_PKMA) == 0.):
            if np.size(options.alpha_PKMA) < options.Niter * options.subsets:
                print('Warning: The number of PKMA alpha (momentum) values must be at least the number of iterations times the number of subsets! Computing custom alpha values.')
                options.alpha_PKMA = np.zeros(options.Niter * options.subsets, dtype=np.float32)
                options.alpha_PKMA[:] = _pkma_relaxation_values(options.Niter, options.subsets, options.rho_PKMA, options.delta_PKMA)
        elif options.PKMA:
            if np.size(options.alpha_PKMA) > options.Niter * options.subsets:
                print('Warning: The number of PKMA alpha (momentum) values is higher than the total number of iterations times subsets. The final values will be ignored.')
    
        # Compute the weights
        if (options.quad or options.L or options.FMH or options.weighted_mean or options.MRP or (options.TV and options.TVtype == 3 and options.TV_use_anatomical) or options.Huber or options.RDP or options.GGMRF or options.hyperbolic) and options.MAP:
            if options.quad or options.L or options.FMH or options.weighted_mean or (options.TV and options.TVtype == 3 and options.TV_use_anatomical) or options.Huber or options.RDP or options.GGMRF or options.hyperbolic:
                if options.GGMRF:
                    computeWeights(options, True)
                else:
                    computeWeights(options, False)
            # These values are needed in order to vectorize the calculation of
            # certain priors
            # Specifies the indices of the center pixel and its neighborhood
            if (options.L or options.FMH):
                raise ValueError('L-filter and FMH-filter are not yet implemented!')
                # options = computeOffsets(options)
            # else:
            #     if options.MRP:
            #         options.medx = options.Ndx * 2 + 1
            #         options.medy = options.Ndy * 2 + 1
            #         options.medz = options.Ndz * 2 + 1
            if options.quad or (options.TV and options.TVtype == 3) or options.GGMRF or options.hyperbolic or (options.RDP and options.RDPIncludeCorners):
                quadWeights(options, options.empty_weight)
            if options.Huber:
                huberWeights(options)
            # if options.RDP:
            #     options = RDPWeights(options)
            if options.L and np.size(options.a_L) == 0:
                raise ValueError('L-filter and FMH-filter are not yet implemented!')
                # options.a_L = lfilter_weights(options.Ndx, options.Ndy, options.Ndz, dx, dy, dz, options.oneD_weights)
            if options.FMH:
                raise ValueError('L-filter and FMH-filter are not yet implemented!')
                # options = fmhWeights(options)
            if (options.FMH or options.quad or options.Huber) and options.implementation == 2:
                options.weights = options.weights.astype(np.float32)
                options.inffi = np.where(np.isinf(options.weights))[0]
                if len(options.inffi) == 0:
                    options.inffi = options.weights.size // 2
            if options.weighted_mean:
                weightedWeights(options)
            if options.RDP and options.RDPIncludeCorners and options.RDP_use_anatomical:
                options.RDP_referenceImage = _load_reference_image(
                    options, options.RDP_referenceImage,
                    emptyMessage='RDP with anatomical weighting selected, but no reference image provided!',
                    resize=False,
                    finalize=True, castBeforeAsfortran=False, doAsfortran=False,
                )
            if options.verbose:
                print('Prepass phase for MRP, quadratic prior, L-filter, FMH, RDP and weighted mean completed')
        if (options.NLM and options.MAP):
            NLMPrepass(options)
    if options.PDHG or options.PDHGKL or options.PDHGL1 or options.PDDY:
        if not isinstance(options.thetaCP, np.ndarray):
            options.thetaCP = np.array(options.thetaCP, dtype=np.float32, ndmin=1)
        if np.size(options.thetaCP) != options.subsets * options.Niter:
            if np.size(options.thetaCP) > 1:
                raise ValueError('The number of elements in options.thetaCP has to be either one or options.subsets * options.Niter!')
            options.thetaCP = np.tile(options.thetaCP, (options.subsets * options.Niter, 1)).astype(dtype=np.float32)
    
    if (options.PKMA or options.MBSREM or options.SPS) and (options.ProxTV or options.TGV):
        if not isinstance(options.thetaCP, np.ndarray):
            options.thetaCP = np.array(options.thetaCP, dtype=np.float32, ndmin=1)
        if np.size(options.thetaCP) != options.subsets * options.Niter and np.size(options.alpha_PKMA) != options.subsets * options.Niter:
            options.thetaCP = np.zeros((options.Niter * options.subsets, 1), order='F', dtype=np.float32)
            options.thetaCP[:, 0] = _pkma_relaxation_values(options.Niter, options.subsets, options.rho_PKMA, options.delta_PKMA)
        else:
            options.thetaCP = options.alpha_PKMA.astype(dtype=np.float32)
        
    
    if options.PDHG or options.PDHGKL or options.PDHGL1 or options.ProxTV or options.TGV or options.FISTA or options.FISTAL1 or options.PDDY:
        if np.size(options.tauCP) < options.nMultiVolumes + 1:
            options.tauCP = np.repeat(options.tauCP, options.nMultiVolumes + 1).astype(dtype=np.float32)
        if np.size(options.sigmaCP) < options.nMultiVolumes + 1:
            options.sigmaCP = np.repeat(options.sigmaCP, options.nMultiVolumes + 1).astype(dtype=np.float32)
        if np.size(options.sigma2CP) < options.nMultiVolumes + 1:
            options.sigma2CP = np.repeat(options.sigma2CP, options.nMultiVolumes + 1).astype(dtype=np.float32)
        if np.size(options.tauCPFilt) < options.nMultiVolumes + 1:
            options.tauCPFilt = np.repeat(options.tauCPFilt, options.nMultiVolumes + 1).astype(dtype=np.float32)
        # if options.implementation == 1 or options.implementation == 4 or options.implementation == 5:
        #     if options.filteringIterations > 0 and options.precondTypeMeas[1]:
        #         apu = options.tauCP.copy()
        #         options.tauCP = options.tauCPFilt.copy()
        #         options.tauCPFilt = apu.copy()
    
    if options.precondTypeImage[5]:
        options.Nf = (2 ** np.ceil(np.log2(options.Nx[0].item()))).astype(dtype=np.uint32).item()
        options.filterIm = rampFilt(options.Nf, options.filterWindow, options.cutoffFrequency, options.normalFilterSigma, True)
        options.filterIm = options.filterIm.astype(dtype=np.float32)

    if options.precondTypeImage[3]:
        if np.size(options.alphaPrecond) < options.Niter * options.subsets:
            print('Warning: The number of alpha (momentum) values must be at least the number of iterations times the number of subsets! Computing custom alpha values.')
            options.alphaPrecond = np.zeros(options.Niter * options.subsets, dtype=np.float32)
            options.alphaPrecond[:] = _pkma_relaxation_values(options.Niter, options.subsets, options.rho_PKMA, options.delta_PKMA)
    
    if options.precondTypeMeas[1]:
        if options.subsets > 1 and options.subsetType == 5:
            options.Nf = 2 ** np.ceil(np.log2(options.nColsD))
        else:
            options.Nf = 2 ** np.ceil(np.log2(options.nRowsD))
        options.Nf = options.Nf.astype(dtype=np.uint32).item()
        options.filter0 = rampFilt(options.Nf, options.filterWindow, options.cutoffFrequency, options.normalFilterSigma)
        options.filter0[0] = 1e-6
        options.filter0 = options.filter0.astype(dtype=np.float32)
        if isinstance(options.sigmaCP, np.ndarray):
            options.Ffilter = np.fft.ifft(options.filter0) * options.sigmaCP[0].item()
        else:
            options.Ffilter = np.fft.ifft(options.filter0) * options.sigmaCP
        if (options.PDHG or options.PDHGL1 or options.PDHGKL or options.CV or options.PDDY) and options.TGV:
            options.Ffilter = options.Ffilter * options.beta
        options.Ffilter[0] = options.Ffilter[0] + 1
        options.Ffilter = np.real(np.fft.fft(options.Ffilter)).astype(dtype=np.float32)
        if options.subsets > 1 and options.subsetType == 5:
            options.filter2 = np.fft.ifft(rampFilt(options.nColsD, options.filterWindow, options.cutoffFrequency, options.normalFilterSigma))
        else:
            options.filter2 = np.fft.ifft(rampFilt(options.nRowsD, options.filterWindow, options.cutoffFrequency, options.normalFilterSigma))
        options.filter2 = np.real(options.filter2).astype(dtype=np.float32)
    if options.FDK and options.CT and options.useFDKWeights:
        options.sourceToCRot = options.sourceToDetector
    if isinstance(options.referenceImage, str):
        options.referenceImage = np.empty(0, dtype=np.float32)
    if isinstance(options.APLS_ref_image, str):
        options.APLS_ref_image = np.empty(0, dtype=np.float32)
    if isinstance(options.NLM_referenceImage, str):
        options.NLM_referenceImage = np.empty(0, dtype=np.float32)
    if isinstance(options.RDP_referenceImage, str):
        options.RDP_referenceImage = np.empty(0, dtype=np.float32)
    if isinstance(options.TV_referenceImage, str):
        options.TV_referenceImage = np.empty(0, dtype=np.float32)
    if (options.PKMA or options.MBSREM or options.SPS or options.RAMLA or options.BSREM or options.ROSEM or options.ROSEM_MAP or options.MRAMLA or options.SAGA) and (options.precondTypeMeas[1] or options.precondTypeImage[5]):
        if (options.lambdaFiltered.size != options.lambdaN.size and options.filteringIterations < options.subsets * options.Niter) or options.lambdaFiltered.size == 0:
            options.lambdaFiltered = options.lambdaN
            