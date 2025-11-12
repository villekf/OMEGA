# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 14:08:02 2024

@author: Ville-Veikko Wettenhovi
"""

def CTEFOVCorrection(options, extrapLengthTransaxial = None, extrapLengthAxial = None, eFOVLengthTransaxial = None, eFOVLengthAxial = None):
    import numpy as np
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
        print('Extending the FOV')
        if not(options.transaxialEFOV) and not(options.axialEFOV):
            print('Neither transaxial nor axial EFOV selected, but EFOV itself is selected! Setting axial EFOV to True!')
            options.axialEFOV = True
        if options.transaxialEFOV:
            if eFOVLengthTransaxial == None:
                if hasattr(options,'eFOVLength'):
                    if options.eFOVLength == None:
                        nTransaxial = int(np.floor(options.Nx * 0.4)) * 2
                    else:
                        nTransaxial = int(np.floor(options.Nx * options.eFOVLength)) * 2
                else:
                    nTransaxial = int(np.floor(options.Nx * 0.4)) * 2
            else:
                nTransaxial = int(np.floor(options.Nx * eFOVLengthTransaxial)) * 2
            options.NxOrig = options.Nx
            options.NyOrig = options.Ny
            options.Nx += nTransaxial
            options.Ny += nTransaxial
            options.FOVxOrig = options.FOVa_x
            options.FOVyOrig = options.FOVa_y
            options.FOVa_x += options.FOVa_x / options.NxOrig * nTransaxial
            options.FOVa_y += options.FOVa_y / options.NyOrig * nTransaxial
        else:
            options.FOVxOrig = options.FOVa_x
            options.FOVyOrig = options.FOVa_y
            options.NxOrig = options.Nx
            options.NyOrig = options.Ny
        if options.axialEFOV:
            if eFOVLengthAxial == None:
                if hasattr(options,'eFOVLength'):
                    if options.eFOVLength == None:
                        if options.sourceToDetector > options.sourceToCRot:
                            pituus = options.sourceToDetector - options.sourceToCRot + options.FOVa_x / 2.
                            angle = (options.sourceToDetector / (options.nColsD * options.dPitchY))
                            eFOVLengthAxial = ((pituus / angle - options.axial_fov / 2.) / (options.axial_fov / 2.)) / 2.
                            nAxial = int(np.floor(options.Nz * eFOVLengthAxial)) * 2
                        else:
                            nAxial = int(np.floor(options.Nz * 0.3)) * 2
                    else:
                        nAxial = int(np.floor(options.Nz * options.eFOVLength)) * 2
                else:
                    if options.sourceToDetector > options.sourceToCRot:
                        pituus = options.sourceToDetector - options.sourceToCRot + options.FOVa_x / 2.
                        angle = (options.sourceToDetector / (options.nColsD * options.dPitchY))
                        eFOVLengthAxial = ((pituus / angle - options.axial_fov / 2.) / (options.axial_fov / 2.)) / 2.
                        nAxial = int(np.floor(options.Nz * eFOVLengthAxial)) * 2
                    else:
                        nAxial = int(np.floor(options.Nz * 0.3)) * 2
            else:
                nAxial = int(np.floor(options.Nz * eFOVLengthAxial)) * 2
            options.NzOrig = options.Nz
            options.Nz += nAxial
            options.axialFOVOrig = options.axial_fov
            options.axial_fov += options.axial_fov / options.NzOrig * nAxial
        else:
            options.axialFOVOrig = options.axial_fov
            options.NzOrig = options.Nz