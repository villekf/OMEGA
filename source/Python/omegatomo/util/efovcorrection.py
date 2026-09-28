# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 14:08:02 2024

@author: Ville-Veikko Wettenhovi
"""
def _round_away_from_zero(x):
    import numpy as np
    x = np.asarray(x)
    return np.sign(x) * np.floor(np.abs(x) + 0.5)

def CTEFOVCorrection(options, extrapLengthTransaxial = None, extrapLengthAxial = None, eFOVLengthTransaxial = None, eFOVLengthAxial = None):
    import numpy as np
    if options.useEFOV:
        import warnings
        size = getattr(options, 'eFOVSize', None)
        # projectorClass supplies an all-zero size even for legacy callers.
        # Any nonzero explicit geometry takes precedence over legacy settings.
        if size is None or not np.any(size):
            warnings.warn('Legacy EFOV parameters detected. Please use options.eFOVSize '
                          'and options.eFOVShift instead of legacy EFOV lengths and flags. '
                          'Converting legacy parameters for this reconstruction.',
                          UserWarning, stacklevel=2)
            transaxial = getattr(options, 'transaxialEFOV', False)
            axial = getattr(options, 'axialEFOV', False)
            if not (transaxial or axial):
                warnings.warn('Neither transaxial nor axial extended FOV selected! '
                              'Defaulting to axial EFOV!', UserWarning, stacklevel=2)
                axial = True

            def legacy_length(argument, name):
                if argument is not None:
                    return argument
                value = getattr(options, name, None)
                if value is None:
                    value = getattr(options, 'eFOVLength', None)
                return 0.4 if value is None else value

            size = np.zeros(3, dtype=np.float64)
            if transaxial:
                length = legacy_length(eFOVLengthTransaxial, 'eFOVLengthTransaxial')
                n_transaxial = np.floor(options.Nx * length) * 2
                size[0] = options.FOVa_x * (1 + n_transaxial / options.Nx)
                size[1] = options.FOVa_y * (1 + n_transaxial / options.Ny)
            if axial:
                length = legacy_length(eFOVLengthAxial, 'eFOVLengthAxial')
                n_axial = np.floor(options.Nz * length) * 2
                size[2] = options.axial_fov * (1 + n_axial / options.Nz)
        options.eFOVSize = np.array(size, dtype=np.float64, copy=True)
        shift = getattr(options, 'eFOVShift', None)
        options.eFOVShift = (np.zeros(3, dtype=np.float64) if shift is None
                             else np.array(shift, dtype=np.float64, copy=True))
    if options.useExtrapolation:
        print('Extrapolating the projections')
        if extrapLengthTransaxial == None:
            if hasattr(options,'extrapLength'):
                if options.extrapLength == None:
                    PnTr = int(np.floor(options.SinM.shape[0] * 0.3))
                else:
                    PnTr = int(np.floor(options.SinM.shape[0] * options.extrapLength))
            else:
                PnTr = int(np.floor(options.SinM.shape[0] * 0.3))
        else:
            PnTr = int(np.floor(options.SinM.shape[0] * extrapLengthTransaxial))
        if extrapLengthAxial == None:
            if hasattr(options,'extrapLength'):
                if options.extrapLength == None:
                    PnAx = int(np.floor(options.SinM.shape[1] * 0.25))
                else:
                    PnAx = int(np.floor(options.SinM.shape[1] * options.extrapLength))
            else:
                PnAx = int(np.floor(options.SinM.shape[1] * 0.25))
        else:
            PnAx = int(np.floor(options.SinM.shape[1] * extrapLengthAxial))
        if options.transaxialExtrapolation:
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
        newProj[erotus1 // 2 : options.SinM.shape[0] + erotus1 // 2, erotus2 // 2 : options.SinM.shape[1] + erotus2 // 2, :] = options.SinM
        if options.transaxialExtrapolation:
            apu = np.tile(options.SinM[0,:,:], (erotus1 // 2, 1, 1)) + 1e-10
            if options.useExtrapolationWeighting:
                apu = np.log(np.single(options.flat) / apu)
                pituus = int(np.round(apu.shape[0] / (6/6)))
                pituus2 = apu.shape[0] - pituus
                if pituus2 == 0:
                    apu = apu * (np.log(np.linspace(1, np.exp(1), pituus)).reshape((-1, 1, 1)) + 1e-10)
                else:
                    apu = apu * (np.concatenate((np.zeros((pituus2), dtype=apu.dtype), np.log(np.linspace(1, np.exp(1), pituus)))).reshape((-1, 1, 1)) + 1e-10)
                apu = np.single(options.flat) / np.exp(apu)
            newProj[0: erotus1 // 2, erotus2 // 2 : options.SinM.shape[1] + erotus2 // 2, :] = apu
            apu = np.tile(options.SinM[-1,:,:], (erotus1 // 2, 1, 1)) + 1e-10
            if options.useExtrapolationWeighting:
                apu = np.log(np.single(options.flat) / apu)
                if pituus2 == 0:
                    apu = apu * (np.log(np.linspace(np.exp(1), 1, pituus)).reshape(-1, 1, 1) + 1e-10)
                else:
                    apu = apu * (np.concatenate((np.log(np.linspace(np.exp(1), 1, pituus))), np.zeros((pituus2), dtype=apu.dtype)).reshape(-1, 1, 1) + 1e-10)
                apu = np.single(options.flat) / np.exp(apu)
            newProj[options.SinM.shape[0] + erotus1 // 2 : , erotus2 // 2 : options.SinM.shape[1] + erotus2 // 2, :] = apu
        if options.axialExtrapolation:
            apu = np.tile(newProj[:,erotus2 // 2, :].reshape(newProj.shape[0], 1, newProj.shape[2]), (1, erotus2 // 2, 1)) + 1e-10
            if options.useExtrapolationWeighting:
                apu = np.log(np.single(options.flat) / apu)
                apu = apu * (np.log(np.linspace(1, np.exp(1), apu.shape[1])).reshape(1, -1, 1) + 1e-10)
                apu = np.single(options.flat) / np.exp(apu)
            newProj[:, : erotus2 // 2, :] = apu
            apu = np.tile(newProj[:,options.SinM.shape[1] + erotus2 // 2 - 1, :].reshape(newProj.shape[0], 1, newProj.shape[2]), (1, erotus2 // 2, 1)) + 1e-10
            if options.useExtrapolationWeighting:
                apu = np.log(np.single(options.flat) / apu)
                apu = apu * (np.log(np.linspace(np.exp(1), 1, apu.shape[1])).reshape(1, -1, 1) + 1e-10)
                apu = np.single(options.flat) / np.exp(apu)
            newProj[:, options.SinM.shape[1] + erotus2 // 2 : , :] = apu
        options.SinM = newProj
        options.nRowsDOrig = options.nRowsD
        options.nColsDOrig = options.nColsD
        options.nRowsD = options.SinM.shape[0]
        options.nColsD = options.SinM.shape[1]
        if options.scatter_correction and options.corrections_during_reconstruction:
            newProj = np.zeros((size1, size2, options.ScatterC.shape[2]), dtype=options.ScatterC.dtype)
            newProj[erotus1 // 2 : options.ScatterC.shape[0] + erotus1 // 2, erotus2 // 2 : options.ScatterC.shape[1] + erotus2 // 2,:] = options.ScatterC
            if options.transaxialExtrapolation:
                apu = np.tile(np.reshape(options.ScatterC[0,:,:], (1, options.ScatterC.shape[1], options.ScatterC.shape[2])), (erotus1 // 2, 1, 1))
                if options.useExtrapolationWeighting:
                    apu = np.log(np.single(options.flat) / apu)
                    pituus = int(np.round(apu.shape[0] / (6/6)))
                    pituus2 = apu.shape[0] - pituus
                    apu = apu * np.log(np.linspace(1, np.exp(1), pituus)).reshape(-1, 1, 1)
                    apu = np.single(options.flat) / np.exp(apu)
                newProj[0: erotus1 // 2, erotus2 // 2 : options.ScatterC.shape[1] + erotus2 // 2, :] = apu
                apu = np.tile(options.ScatterC[-1,:,:], (erotus1 // 2, 1, 1))
                if options.useExtrapolationWeighting:
                    apu = np.log(np.single(options.flat) / apu)
                    apu = apu * np.log(np.linspace(np.exp(1), 1, pituus)).reshape(-1, 1, 1)
                    apu = np.single(options.flat) / np.exp(apu)
                newProj[options.ScatterC.shape[0] + erotus1 // 2 : , erotus2 // 2 : options.ScatterC.shape[1] + erotus2 // 2, :] = apu
            if options.axialExtrapolation:
                apu = np.tile(np.reshape(newProj[:,erotus2 // 2, :], (newProj.shape[0], 1, newProj.shape[2])), (1, erotus2 // 2, 1))
                if options.useExtrapolationWeighting:
                    apu = np.log(np.single(options.flat) / apu)
                    apu = apu * np.reshape(np.log(np.linspace(1, np.exp(1), apu.shape[1])), (1, -1, 1))
                    apu = np.single(options.flat) / np.exp(apu)
                newProj[:, 0: erotus2 // 2, :] = apu
                apu = np.tile(np.reshape(newProj[:,options.ScatterC.shape[1] + erotus2 // 2, :], (newProj.shape[0], 1, newProj.shape[2])), (1, erotus2 // 2, 1))
                if options.useExtrapolationWeighting:
                    apu = np.log(np.single(options.flat) / apu)
                    apu = apu * np.reshape(np.log(np.linspace(np.exp(1), 1, apu.shape[1])), (1, -1, 1))
                    apu = np.single(options.flat) / np.exp(apu)
                newProj[:, options.ScatterC.shape[1] + erotus2 // 2 : , :] = apu
            options.ScatterC = newProj
    
    if options.useEFOV:
        options.axialEFOV = False # Check if axial EFOV is inside FOV (after shift). If is outside (in both directions), set axialEFOV to true
        FOVmin_z = -options.axial_fov / 2.0
        FOVmax_z =  options.axial_fov / 2.0
        print(options.eFOVShift)
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
