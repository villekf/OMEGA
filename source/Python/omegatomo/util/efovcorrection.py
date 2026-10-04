# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 14:08:02 2024

@author: Ville-Veikko Wettenhovi
"""
def _round_away_from_zero(x):
    import numpy as np
    x = np.asarray(x)
    return np.sign(x) * np.floor(np.abs(x) + 0.5)

import warnings

import numpy as np

from .matlabRound import matlabRound


def _getOpt(options, name):
    """Return options.<name> if it exists and is not None, else None.

    This mirrors MATLAB's isfield(options, name) semantics for the Python
    projectorClass, whose attributes are always present but default to
    None when "not specified" is meant.
    """
    if hasattr(options, name):
        val = getattr(options, name)
        if val is not None:
            return val
    return None


def CTEFOVCorrection(options, extrapLengthTransaxial = None, extrapLengthAxial = None, eFOVLengthTransaxial = None, eFOVLengthAxial = None):
    """Extended FOV correction for (CB)CT data.

    Python port of source/m-files/CTEFOVCorrection.m, which is the
    reference implementation. The MATLAB version determines the
    transaxial/axial extrapolation and extended-FOV lengths from the
    parameters given in ``options`` (falling back to
    ``options.extrapLength`` / ``options.eFOVLength`` and finally to
    hard-coded defaults). For the axial extended FOV specifically, the
    length is instead derived automatically from the cone-beam geometry
    (``options.sourceToDetector``, ``options.sourceToCRot``,
    ``options.nColsD``, ``options.dPitchY``) whenever
    ``eFOVLengthAxial`` was not explicitly given and the geometry is
    available (``sourceToDetector > sourceToCRot``) -- this geometry-based
    override happens regardless of whether ``eFOVLength`` was set;
    ``eFOVLength``/``eFOVLengthAxial`` are only used as the *default* value
    of the axial extension when no usable geometry is present, or when the
    axial extension is explicitly given.

    The keyword arguments below mirror MATLAB's isfield-based field
    precedence: when given explicitly (not None) they take precedence
    exactly like an explicitly-set MATLAB options field. Otherwise the
    same-named ``options`` attribute is used if present (not None), and
    finally the hard-coded default described below.

    Parameters
    ----------
    options : projectorClass
        The main options/parameter object. ``options.SinM``,
        ``options.flat``, ``options.nProjections`` etc. must already be
        set when ``options.useExtrapolation`` is True, and
        ``options.Nx``/``Ny``/``Nz``/``FOVa_x``/``FOVa_y``/``axial_fov``
        (and the cone-beam geometry, for the automatic axial case) must
        already be set when ``options.useEFOV`` is True.
    extrapLengthTransaxial : float, optional
        Transaxial extrapolation length (fraction of the detector size,
        per side). Overrides ``options.extrapLengthTransaxial`` and
        ``options.extrapLength``. Default (if nothing else is set) is 0.25.
    extrapLengthAxial : float, optional
        Axial extrapolation length (fraction of the detector size, per
        side). Overrides ``options.extrapLengthAxial`` and
        ``options.extrapLength``. Default (if nothing else is set) is 0.25.
    eFOVLengthTransaxial : float, optional
        Transaxial extended FOV length (fraction of Nx/Ny, per side).
        Overrides ``options.eFOVLengthTransaxial`` and
        ``options.eFOVLength``. Default (if nothing else is set) is 0.4.
    eFOVLengthAxial : float, optional
        Axial extended FOV length (fraction of Nz, per side). Overrides
        ``options.eFOVLengthAxial`` and ``options.eFOVLength`` and disables
        the automatic geometry-based derivation described above. Default
        (if nothing else is set and no usable cone-beam geometry is
        present) is 0.3.

    Returns
    -------
    options : projectorClass
        The (mutated) options object.
    """
    if not hasattr(options, 'useExtrapolation') or options.useExtrapolation is None:
        options.useExtrapolation = False
    if not hasattr(options, 'scatter_correction') or options.scatter_correction is None:
        options.scatter_correction = False
    if not hasattr(options, 'useEFOV') or options.useEFOV is None:
        options.useEFOV = False
    if not hasattr(options, 'useInpaint') or options.useInpaint is None:
        options.useInpaint = False

    if options.useEFOV:
        size = getattr(options, 'eFOVSize', None)
        if size is None or not np.any(size):
            warnings.warn('Legacy EFOV parameters detected. Please use options.eFOVSize and options.eFOVShift instead of legacy EFOV lengths and flags. Converting legacy parameters for this reconstruction.', UserWarning, stacklevel=2)
            transaxial = getattr(options, 'transaxialEFOV', False)
            axial = getattr(options, 'axialEFOV', False)
            if not (transaxial or axial):
                warnings.warn('Neither transaxial nor axial extended FOV selected! Defaulting to axial EFOV!', UserWarning, stacklevel=2)
                axial = True
            def legacy_length(argument, name):
                if argument is not None:
                    return argument
                value = _getOpt(options, name)
                if value is None:
                    value = _getOpt(options, 'eFOVLength')
                return 0.4 if value is None else value
            size = np.zeros(3, dtype=np.float64)
            if transaxial:
                length = legacy_length(eFOVLengthTransaxial, 'eFOVLengthTransaxial')
                n_transaxial = np.floor(options.Nx * length) * 2
                size[0] = options.FOVa_x * (1 + n_transaxial / options.Nx)
                size[1] = options.FOVa_y * (1 + n_transaxial / options.Ny)
            if axial:
                length = legacy_length(eFOVLengthAxial, 'eFOVLengthAxial')
                explicit_axial = (eFOVLengthAxial is not None or
                                  _getOpt(options, 'eFOVLengthAxial') is not None)
                if (not explicit_axial and _getOpt(options, 'sourceToDetector') is not None and
                        _getOpt(options, 'sourceToCRot') is not None and
                        options.sourceToDetector > options.sourceToCRot):
                    distance = options.sourceToCRot + options.FOVa_x / 2.0
                    angle = options.sourceToDetector / (options.nColsD * options.dPitchY / 2.0)
                    length = max(0.0, (distance / angle - options.axial_fov / 2.0) / options.axial_fov)
                n_axial = np.floor(options.Nz * length) * 2
                size[2] = options.axial_fov * (1 + n_axial / options.Nz)
        options.eFOVSize = np.array(size, dtype=np.float64, copy=True)
        shift = getattr(options, 'eFOVShift', None)
        options.eFOVShift = (np.zeros(3, dtype=np.float64) if shift is None
                             else np.array(shift, dtype=np.float64, copy=True))

    # --- Resolve extrapolation/eFOV lengths, mirroring MATLAB's isfield
    #     fallback chain: explicit keyword argument > options.<name>Axial /
    #     options.<name>Transaxial > options.<name> (base) > hard-coded
    #     default. ---
    extrapLengthBase = _getOpt(options, 'extrapLength')
    extrapLengthBase = extrapLengthBase if extrapLengthBase is not None else .25

    if extrapLengthAxial is None:
        extrapLengthAxial = _getOpt(options, 'extrapLengthAxial')
        if extrapLengthAxial is None:
            extrapLengthAxial = extrapLengthBase

    if extrapLengthTransaxial is None:
        extrapLengthTransaxial = _getOpt(options, 'extrapLengthTransaxial')
        if extrapLengthTransaxial is None:
            extrapLengthTransaxial = extrapLengthBase

    eFOVLengthBase = _getOpt(options, 'eFOVLength')
    eFOVLengthBase = eFOVLengthBase if eFOVLengthBase is not None else .4

    if eFOVLengthTransaxial is None:
        eFOVLengthTransaxial = _getOpt(options, 'eFOVLengthTransaxial')
        if eFOVLengthTransaxial is None:
            eFOVLengthTransaxial = eFOVLengthBase

    # eFOVLengthAxial needs special handling: unlike every other length
    # above, MATLAB overrides it with a geometry-derived value whenever it
    # was NOT explicitly given (field missing), regardless of whether
    # eFOVLength was set. Track whether it was explicitly given so that the
    # geometry override below can be gated on that alone.
    eFOVLengthAxialExplicit = (eFOVLengthAxial is not None) or (_getOpt(options, 'eFOVLengthAxial') is not None)
    if eFOVLengthAxial is None:
        eFOVLengthAxial = _getOpt(options, 'eFOVLengthAxial')
        if eFOVLengthAxial is None:
            eFOVLengthAxial = eFOVLengthBase

    weighting = options.useExtrapolationWeighting
    offset = options.offsetCorrection

    # --- Normalize the transaxial/axial EFOV and extrapolation direction
    #     flags, warning (like MATLAB) if useEFOV/useExtrapolation is set
    #     but neither direction was selected. ---
    if options.useEFOV:
        if not options.transaxialEFOV and not options.axialEFOV:
            warnings.warn('Neither transaxial nor axial extended FOV selected! Defaulting to axial EFOV!')
            options.transaxialEFOV = False
            options.axialEFOV = True
    if options.useExtrapolation:
        if not options.transaxialExtrapolation and not options.axialExtrapolation:
            warnings.warn('Neither transaxial nor axial extrapolation selected! Defaulting to axial extrapolation!')
            options.transaxialExtrapolation = False
            options.axialExtrapolation = True


    if options.useExtrapolation:
        print('Extrapolating the projections')
        PnTr = int(np.floor(options.SinM.shape[0] * extrapLengthTransaxial))
        PnAx = int(np.floor(options.SinM.shape[1] * extrapLengthAxial))
        if options.useInpaint:
            # useInpaint is an experimental, unofficial MATLAB-only feature
            # (inpaint_nans-based projection extrapolation); it is
            # intentionally not ported to Python.
            raise ValueError('Projection inpainting (useInpaint) is not supported in Python')
        if options.transaxialExtrapolation:
            if offset:
                size1 = options.SinM.shape[0] + PnTr
            else:
                size1 = options.SinM.shape[0] + PnTr * 2
        else:
            size1 = options.SinM.shape[0]
        if options.axialExtrapolation:
            size2 = options.SinM.shape[1] + PnAx * 2
        else:
            size2 = options.SinM.shape[1]
        erotus1 = size1 - options.SinM.shape[0]
        erotus2 = size2 - options.SinM.shape[1]
        newProj = np.zeros((size1, size2, options.SinM.shape[2]), dtype=options.SinM.dtype)
        if offset:
            newProj[erotus1: options.SinM.shape[0] + erotus1, erotus2 // 2: options.SinM.shape[1] + erotus2 // 2, :] = options.SinM
        else:
            newProj[erotus1 // 2: options.SinM.shape[0] + erotus1 // 2, erotus2 // 2: options.SinM.shape[1] + erotus2 // 2, :] = options.SinM
        if options.transaxialExtrapolation:
            if offset:
                apu = np.tile(options.SinM[0,:,:], (erotus1, 1, 1)) + 1e-10
            else:
                apu = np.tile(options.SinM[0,:,:], (erotus1 // 2, 1, 1)) + 1e-10
            if weighting:
                apu = np.log(np.single(options.flat) / apu)
                pituus = int(matlabRound(apu.shape[0] / (6/6)))
                pituus2 = apu.shape[0] - pituus
                if pituus2 == 0:
                    apu = apu * (np.log(np.linspace(1, np.exp(1), pituus)).reshape((-1, 1, 1)) + 1e-10)
                else:
                    apu = apu * (np.concatenate((np.zeros((pituus2), dtype=apu.dtype), np.log(np.linspace(1, np.exp(1), pituus)))).reshape((-1, 1, 1)) + 1e-10)
                apu = np.single(options.flat) / np.exp(apu)

            if not offset:
                newProj[0: erotus1 // 2, erotus2 // 2 : options.SinM.shape[1] + erotus2 // 2, :] = apu
                apu = np.tile(options.SinM[-1,:,:], (erotus1 // 2, 1, 1)) + 1e-10
                if weighting:
                    apu = np.log(np.single(options.flat) / apu)
                    if pituus2 == 0:
                        apu = apu * (np.log(np.linspace(np.exp(1), 1, pituus)).reshape(-1, 1, 1) + 1e-10)
                    else:
                        apu = apu * (np.concatenate((np.log(np.linspace(np.exp(1), 1, pituus)), np.zeros((pituus2), dtype=apu.dtype))).reshape(-1, 1, 1) + 1e-10)
                    apu = np.single(options.flat) / np.exp(apu)
                newProj[options.SinM.shape[0] + erotus1 // 2 : , erotus2 // 2 : options.SinM.shape[1] + erotus2 // 2, :] = apu
            else:
                newProj[0: erotus1, erotus2 // 2 : options.SinM.shape[1] + erotus2 // 2, :] = apu
        if options.axialExtrapolation:
            apu = np.tile(newProj[:,erotus2 // 2, :].reshape(newProj.shape[0], 1, newProj.shape[2]), (1, erotus2 // 2, 1)) + 1e-10
            if weighting:
                apu = np.log(np.single(options.flat) / apu)
                apu = apu * (np.log(np.linspace(1, np.exp(1), apu.shape[1])).reshape(1, -1, 1) + 1e-10)
                apu = np.single(options.flat) / np.exp(apu)
            newProj[:, : erotus2 // 2, :] = apu
            apu = np.tile(newProj[:,options.SinM.shape[1] + erotus2 // 2 - 1, :].reshape(newProj.shape[0], 1, newProj.shape[2]), (1, erotus2 // 2, 1)) + 1e-10
            if weighting:
                apu = np.log(np.single(options.flat) / apu)
                apu = apu * (np.log(np.linspace(np.exp(1), 1, apu.shape[1])).reshape(1, -1, 1) + 1e-10)
                apu = np.single(options.flat) / np.exp(apu)
            newProj[:, options.SinM.shape[1] + erotus2 // 2 : , :] = apu
        options.SinM = newProj
        if options.scatter_correction and options.corrections_during_reconstruction:
            newProj = np.zeros((size1, size2, options.ScatterC.shape[2]), dtype=options.ScatterC.dtype)
            newProj[erotus1 // 2 : options.ScatterC.shape[0] + erotus1 // 2, erotus2 // 2 : options.ScatterC.shape[1] + erotus2 // 2,:] = options.ScatterC
            if options.transaxialExtrapolation:
                apu = np.tile(np.reshape(options.ScatterC[0,:,:], (1, options.ScatterC.shape[1], options.ScatterC.shape[2])), (erotus1 // 2, 1, 1))
                if weighting:
                    apu = np.log(np.single(options.flat) / apu)
                    pituus = int(matlabRound(apu.shape[0] / (6/6)))
                    pituus2 = apu.shape[0] - pituus
                    apu = apu * np.log(np.linspace(1, np.exp(1), pituus)).reshape(-1, 1, 1)
                    apu = np.single(options.flat) / np.exp(apu)
                newProj[0: erotus1 // 2, erotus2 // 2 : options.ScatterC.shape[1] + erotus2 // 2, :] = apu
                apu = np.tile(options.ScatterC[-1,:,:], (erotus1 // 2, 1, 1))
                if weighting:
                    apu = np.log(np.single(options.flat) / apu)
                    apu = apu * np.log(np.linspace(np.exp(1), 1, pituus)).reshape(-1, 1, 1)
                    apu = np.single(options.flat) / np.exp(apu)
                newProj[options.ScatterC.shape[0] + erotus1 // 2 : , erotus2 // 2 : options.ScatterC.shape[1] + erotus2 // 2, :] = apu
            if options.axialExtrapolation:
                apu = np.tile(np.reshape(newProj[:,erotus2 // 2, :], (newProj.shape[0], 1, newProj.shape[2])), (1, erotus2 // 2, 1))
                if weighting:
                    apu = np.log(np.single(options.flat) / apu)
                    apu = apu * np.reshape(np.log(np.linspace(1, np.exp(1), apu.shape[1])), (1, -1, 1))
                    apu = np.single(options.flat) / np.exp(apu)
                newProj[:, 0: erotus2 // 2, :] = apu
                apu = np.tile(np.reshape(newProj[:,options.ScatterC.shape[1] + erotus2 // 2, :], (newProj.shape[0], 1, newProj.shape[2])), (1, erotus2 // 2, 1))
                if weighting:
                    apu = np.log(np.single(options.flat) / apu)
                    apu = apu * np.reshape(np.log(np.linspace(np.exp(1), 1, apu.shape[1])), (1, -1, 1))
                    apu = np.single(options.flat) / np.exp(apu)
                newProj[:, options.ScatterC.shape[1] + erotus2 // 2 : , :] = apu
            options.ScatterC = newProj
        options.nRowsDOrig = options.nRowsD
        options.nColsDOrig = options.nColsD
        options.nRowsD = options.SinM.shape[0]
        options.nColsD = options.SinM.shape[1]
    
    if options.useEFOV:
        options.axialEFOV = False # Check if axial EFOV is inside FOV (after shift). If is outside (in both directions), set axialEFOV to true
        FOVmin_z = -options.axial_fov / 2.0
        FOVmax_z =  options.axial_fov / 2.0
        eFOVmin_z = -options.eFOVSize[2] / 2.0 + options.eFOVShift[2]
        eFOVmax_z =  options.eFOVSize[2] / 2.0 + options.eFOVShift[2]
        if FOVmin_z < eFOVmin_z or FOVmax_z > eFOVmax_z: # FOV not entirely inside eFOV
            print('The high-resolution FOV is not entirely inside the extended FOV in z-direction. No extension will be performed in the axial direction.')
            options.eFOVShift[2] = 0
            options.eFOVSize[2] = options.axial_fov
        else:
            options.axialEFOV = True
            
        options.transaxialEFOV = False
        FOVmin_x = -options.FOVa_x / 2.0
        FOVmax_x =  options.FOVa_x / 2.0
        eFOVmin_x = -options.eFOVSize[0] / 2.0 + options.eFOVShift[0]
        eFOVmax_x =  options.eFOVSize[0] / 2.0 + options.eFOVShift[0]

        FOVmin_y = -options.FOVa_y / 2.0
        FOVmax_y =  options.FOVa_y / 2.0
        eFOVmin_y = -options.eFOVSize[1] / 2.0 + options.eFOVShift[1]
        eFOVmax_y =  options.eFOVSize[1] / 2.0 + options.eFOVShift[1]

        if (FOVmin_x < eFOVmin_x or FOVmax_x > eFOVmax_x or
            FOVmin_y < eFOVmin_y or FOVmax_y > eFOVmax_y):
            print('Warning: The high-resolution FOV is not entirely inside the extended FOV in xy-direction. '
                'No extension will be performed in the transaxial direction.')
            options.eFOVShift[0] = 0
            options.eFOVSize[0] = options.FOVa_x
            options.eFOVShift[1] = 0
            options.eFOVSize[1] = options.FOVa_y
        else:
            options.transaxialEFOV = True

        if not (options.axialEFOV or options.transaxialEFOV):
            options.useEFOV = False
            print('Warning: FOV extension is not performed; turning off options.useEFOV')

    if options.useEFOV:
        print('Extending the FOV')
        options.FOVxOrig = options.FOVa_x
        options.FOVyOrig = options.FOVa_y
        options.axialFOVOrig = options.axial_fov
        options.NxOrig = options.Nx
        options.NyOrig = options.Ny
        options.NzOrig = options.Nz
        
        if options.transaxialEFOV:
            options.FOVa_x = options.eFOVSize[0]
            options.FOVa_y = options.eFOVSize[1]
            options.Nx = int(np.ceil(options.Nx * options.FOVa_x / options.FOVxOrig))
            options.Ny = int(np.ceil(options.Ny * options.FOVa_y / options.FOVyOrig))

        if options.axialEFOV:
            options.axial_fov = options.eFOVSize[2]
            options.Nz = int(np.ceil(options.Nz * options.axial_fov / options.axialFOVOrig))

        dx = options.FOVa_x / options.Nx
        dy = options.FOVa_y / options.Ny
        dz = options.axial_fov / options.Nz
        
        options.eFOVShift_Nx = int(_round_away_from_zero(options.eFOVShift[0] / dx))
        options.eFOVShift_Ny = int(_round_away_from_zero(options.eFOVShift[1] / dy))
        options.eFOVShift_Nz = int(_round_away_from_zero(options.eFOVShift[2] / dz))
    # Multiresolution total FOV size. numel(FOVa_x) is 1 unless this
    # function is (re-)run after setUpCorrections has already split the
    # volume into multi-resolution sub-volumes (in which case FOVa_x/FOVa_y
    # /axial_fov hold 3, 5 or 7 per-volume entries); mirrors
    # CTEFOVCorrection.m.
    FOVa_xArr = np.atleast_1d(options.FOVa_x)
    FOVa_yArr = np.atleast_1d(options.FOVa_y)
    axialFovArr = np.atleast_1d(options.axial_fov)
    if FOVa_xArr.size == 1:  # No eFOV
        FOV = np.array([FOVa_xArr[0], FOVa_yArr[0], axialFovArr[0]])
    elif FOVa_xArr.size == 3:  # Axial eFOV only
        FOV = np.array([FOVa_xArr[0], FOVa_yArr[0], np.sum(axialFovArr)])
    elif FOVa_xArr.size == 5:  # Transaxial eFOV only
        FOV = np.array([
            FOVa_xArr[0] + FOVa_xArr[1] + FOVa_xArr[2],
            FOVa_yArr[0] + FOVa_yArr[3] + FOVa_yArr[4],
            axialFovArr[0]
        ])
    elif FOVa_xArr.size == 7:  # Axial + transaxial eFOV
        FOV = np.array([
            FOVa_xArr[0] + FOVa_xArr[3] + FOVa_xArr[4],
            FOVa_yArr[0] + FOVa_yArr[5] + FOVa_yArr[6],
            axialFovArr[0] + axialFovArr[1] + axialFovArr[2]
        ])
    else:
        FOV = np.array([FOVa_xArr[0], FOVa_yArr[0], axialFovArr[0]])

    if options.ellipseParametersDerived or options.ellipseRadiusX == 0 or options.ellipseRadiusY == 0 or options.ellipseRadiusZ == 0:
        options.ellipseRadiusX = FOV[0] / 2
        options.ellipseRadiusY = FOV[1] / 2
        options.ellipseRadiusZ = FOV[2] / 2
        options.ellipseParametersDerived = True
    if not options.ellipseCenterOffsetApplied:
        options.ellipseCenterX += options.oOffsetX
        options.ellipseCenterY += options.oOffsetY
        options.ellipseCenterZ += options.oOffsetZ
        options.ellipseCenterOffsetApplied = True

    return options
