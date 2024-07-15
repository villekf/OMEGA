# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 13:57:46 2024

@author: Ville-Veikko Wettenhovi
"""

import numpy as np
import os
import tkinter as tk
from tkinter.filedialog import askopenfilename
from .rampfilt import rampFilt

def linearizeData(options):
    options.SinM = np.log(options.flat / options.SinM.astype(dtype=np.float32))

def loadCorrections(options):
    if options.attenuation_correction == 1:
        if options.vaimennus.size == 0:
            if len(options.attenuation_datafile) > 0 and options.attenuation_datafile[len(options.attenuation_datafile)-3:len(options.attenuation_datafile)+1:1] == 'mhd':
                try:
                    from SimpleITK import ReadImage as loadMetaImage
                    from SimpleITK import GetArrayFromImage
                    metaImage = loadMetaImage(options.attenuation_datafile)
                    options.vaimennus = GetArrayFromImage(metaImage)
                    options.vaimennus = np.asfortranarray(np.transpose(options.vaimennus, (2, 1, 0)))
                    apu = np.array(list(metaImage.GetSpacing()))
                    if options.CT_attenuation:
                        if round(apu[0].item()*100.)/100. > round(options.FOVa_x[0].item() / (options.Nx[0].item())*100.)/100. or round(apu[0].item()*100)/100 < round(options.FOVa_x[0].item() / (options.Nx[0].item())*100.)/100.:
                            options.vaimennus = options.vaimennus * (apu[0].item() / (options.FOVa_x[0].item() / (options.Nx[0].item())))
                except ModuleNotFoundError:
                    print('SimpleITK package not found! MetaImages cannot be loaded. You can install SimpleITK package with "pip install SimpleITK".')
            elif len(options.attenuation_datafile) > 0 and  options.attenuation_datafile[len(options.attenuation_datafile)-3:len(options.attenuation_datafile)+1:1] == 'mat':
                try:
                    from pymatreader import read_mat
                    var = read_mat(options.attenuation_datafile)
                    options.vaimennus = np.array(next(iter(var.items()))[1])
                except ModuleNotFoundError:
                    print('pymatreader package not found! Mat-files cannot be loaded. You can install pymatreader package with "pip install pymatreader".')
            elif len(options.attenuation_datafile) > 0 and  (options.attenuation_datafile[len(options.attenuation_datafile)-3:len(options.attenuation_datafile)+1:1] == 'npy' or options.attenuation_datafile[len(options.attenuation_datafile)-3:len(options.attenuation_datafile)+1:1] == 'npz'):
                apu = np.load(options.attenuation_datafile, allow_pickle=True)
                variables = list(apu.keys())
                options.vaimennus = apu[variables[0]]
            else:
                root = tk.Tk()
                root.withdraw()
                nimi = askopenfilename(title='Select attenuation datafile',filetypes=(('MHD, NPY, NPZ and MAT files','*.mhd *.mat *.npy *.npz'),('All','*.*')))
                if len(nimi) == 0:
                    raise ValueError("No file selected!")
                if nimi[len(nimi)-3:len(nimi)+1:1] == 'mhd':
                    try:
                        from SimpleITK import ReadImage as loadMetaImage
                        from SimpleITK import GetArrayFromImage
                        from SimpleITK import ReadImage as loadMetaImage
                        from SimpleITK import GetArrayFromImage
                        metaImage = loadMetaImage(options.attenuation_datafile)
                        options.vaimennus = GetArrayFromImage(metaImage)
                        apu = np.array(list(metaImage.GetSpacing()))
                        if options.CT_attenuation:
                            if round(apu[0].item()*100.)/100. > round(options.FOVa_x[0].item() / (options.Nx[0].item())*100.)/100. or round(apu[0].item()*100)/100 < round(options.FOVa_x[0].item() / (options.Nx[0].item())*100.)/100.:
                                options.vaimennus = options.vaimennus * (apu[0].item() / (options.FOVa_x[0].item() / (options.Nx[0].item())))
                    except ModuleNotFoundError:
                        print('SimpleITK package not found! MetaImages cannot be loaded. You can install SimpleITK package with "pip install SimpleITK".')
                elif nimi[len(nimi)-3:len(nimi)+1:1] == 'mat':
                    try:
                        from pymatreader import read_mat
                        var = read_mat(options.attenuation_datafile)
                        options.vaimennus = np.array(next(iter(var.items()))[1])
                    except ModuleNotFoundError:
                        print('pymatreader package not found! Mat-files cannot be loaded. You can install pymatreader package with "pip install pymatreader".')
                elif (nimi[len(nimi)-3:len(nimi)+1:1] == 'npy' or nimi[len(nimi)-3:len(nimi)+1:1] == 'npz'):
                    apu = np.load(nimi, allow_pickle=True)
                    variables = list(apu.keys())
                    options.vaimennus = apu[variables[0]]
                else:
                    ValueError('Unsupported datatype!')
        if options.CT_attenuation:
            if not(options.vaimennus.shape[0] == options.Nx[0]) or not(options.vaimennus.shape[1] == options.Ny[0].item()) or not(options.vaimennus.shape[2] == options.Nz[0].item()):
                if options.vaimennus.shape[0] != options.N[0]:
                    ValueError('Error: Attenuation data is of different size than the reconstructed image')
            if options.rotateAttImage != 0:
                atn = np.reshape(options.vaimennus, (options.Nx[0].item(), options.Ny[0].item(), options.Nz[0].item()))
                atn = np.rot90(atn,options.rotateAttImage);
                options.vaimennus = atn
            if options.flipAttImageXY:
                atn = np.reshape(options.vaimennus, (options.Nx[0].item(), options.Ny[0].item(), options.Nz[0].item()))
                atn = np.fliplr(atn);
                options.vaimennus = atn
            if options.flipAttImageZ:
                atn = np.reshape(options.vaimennus, (options.Nx[0].item(), options.Ny[0].item(), options.Nz[0].item()))
                atn = np.flip(atn,2);
                options.vaimennus = atn
        options.vaimennus = options.vaimennus.ravel('F').astype(dtype=np.float32)
    if (options.normalization_correction and options.corrections_during_reconstruction):
        if options.normalization.size == 0:
            normdir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', '..', '..', '..', 'mat-files')) + "/" +  options.machine_name + '_normalization_' + str(options.Ndist) + 'x' + str(options.Nang) + '_span' + str(options.span) + '.mat'
            if os.path.exists(normdir):
                try:
                    from pymatreader import read_mat
                    var = read_mat(normdir)
                    options.normalization = np.array(var["normalization"],order='F')
                except ModuleNotFoundError:
                    print('pymatreader package not found! Mat-files cannot be loaded. You can install pymatreader package with "pip install pymatreader".')
            else:
                root = tk.Tk()
                root.withdraw()
                nimi = askopenfilename(title='Select normalization datafile',filetypes=(('NRM, NPY, NPZ and MAT files','*.nrm *.mat *.npy *.npz'),('All','*.*')))
                if len(nimi) == 0:
                    raise ValueError("No file selected!")
                if nimi[len(nimi)-3:len(nimi)+1:1] == 'nrm':
                    options.normalization = np.fromfile(nimi, dtype=np.float32)
                    if options.normalization.size != options.Ndist * options.Nang * options.TotSinos and ~options.use_raw_data:
                        ValueError('Size mismatch between the current data and the normalization data file')
                elif nimi[len(nimi)-3:len(nimi)+1:1] == 'mat':
                    var = read_mat(nimi)
                    options.normalization = np.array(var["normalization"])
                elif (nimi[len(nimi)-3:len(nimi)+1:1] == 'npy' or nimi[len(nimi)-3:len(nimi)+1:1] == 'npz'):
                    apu = np.load(nimi, allow_pickle=True)
                    variables = list(apu.keys())
                    options.normalization = apu[variables[0]]
                else:
                    ValueError('Unsupported datatype!')
            options.normalization = 1. / options.normalization.ravel('F').astype(dtype=np.float32)
            if ~options.use_raw_data and options.NSinos != options.TotSinos:
                options.normalization = options.normalization[0 : options.Ndist * options.Nang * options.NSinos]
        options.normalization = options.normalization.ravel('F').astype(dtype=np.float32)
    elif (options.normalization_correction and not options.corrections_during_reconstruction):
        options.normalization_correction = False
        
def parseInputs(options, mDataFound = False):
    if options.subsets > 1:
        if mDataFound and not options.largeDim:
            if options.Nt > 1:
                for ff in range(1, options.Nt + 1):
                    if not options.use_raw_data:
                        temp = options.SinM[:,:,:,:,ff - 1]
                        if options.NSinos != options.TotSinos:
                            temp = temp[:, :, :options.NSinos, :]
                    # else:
                    #     temp = np.single(np.full(options.SinM[ff - 1]))
            
                    if options.TOF:
                        if options.subsetType >= 8:
                            temp = temp[:, :, options.index, :]
                        else:
                            temp = np.reshape(temp, (temp.size // options.TOF_bins, options.TOF_bins),order='F')
                            temp = temp[options.index, :]
                    else:
                        if options.subsetType >= 8:
                            temp = temp[:, :, options.index]
                        else:
                            temp = temp.ravel(order='F')
                            temp = temp[options.index]
                            temp = np.reshape(temp, (options.nRowsD, options.nColsD, options.nProjections), order='F')
                    options.SinM[:,:,:,:,ff - 1] = temp
            else:
                if not options.use_raw_data:
                    if options.NSinos != options.TotSinos:
                        options.SinM = options.SinM[:, :, :options.NSinos, :]
                # else:
                #     options.SinM = np.single(np.full(options.SinM[0]))
            
                if options.subsetType >= 8:
                    options.SinM = np.reshape(options.SinM, (options.nRowsD, options.nColsD, options.nProjections, options.TOF_bins), order='F')
            
                if options.TOF:
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
            if not options.use_raw_data and options.NSinos != options.TotSinos:
                options.normalization = options.normalization[:options.NSinos * options.Ndist * options.Nang]
            if options.subsetType >= 8:
                options.normalization = np.reshape(options.normalization, (options.Ndist, options.Nang, -1),order='F')
                options.normalization = options.normalization[:, :, options.index]
                options.normalization = options.normalization.ravel(order='F').astype(dtype=np.float32)
            else:
                options.normalization = options.normalization[options.index]
        
        if options.additionalCorrection and hasattr(options, 'corrVector'):
            if options.subsetType >= 8:
                # if isinstance(options.corrVector, list):
                #     for kk in range(len(options.corrVector)):
                #         options.corrVector[kk] = np.reshape(options.corrVector[kk], (options.Ndist, options.Nang, -1))
                #         options.corrVector[kk] = options.corrVector[kk][:, :, index, :]
                #         options.corrVector[kk] = options.corrVector[kk].ravel()
                # else:
                options.corrVector = np.reshape(options.corrVector, (options.Ndist, options.Nang, options.nProjections, -1), order='F')
                options.corrVector = options.corrVector[:, :, options.index, :]
                options.corrVector = options.corrVector.ravel(order='F')
            else:
                # if isinstance(options.corrVector, list):
                #     for kk in range(len(options.corrVector)):
                #         options.corrVector[kk] = options.corrVector[kk][index]
                # else:
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
            if options.subsetType >= 8:
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
            
    
    if mDataFound and not options.largeDim:
        options.SinM = options.SinM.ravel(order='F').astype(dtype=np.float32)


def TVPrepass(options):
    from skimage.transform import resize #scikit-image
    def assembleS(alkuarvo,T,Ny,Nx,Nz):
        S = np.zeros((Nx * Ny * Nz * 3, 3),order='F',dtype=np.float32)
        # Compute the weighting gamma and tensor S
        f = -np.diff(alkuarvo, axis=1, append=0)
        f = np.concatenate((f, np.zeros((Nx, 1, Nz),order='F',dtype=np.float32)), axis=1)
        f = f.reshape(-1)
        g = -np.diff(alkuarvo, axis=0, append=0)
        g = np.concatenate((g, np.zeros((1, Ny, Nz),order='F',dtype=np.float32)), axis=0)
        g = g.reshape(-1)
        h = -np.diff(alkuarvo, axis=2, append=0)
        h = np.concatenate((h, np.zeros((Nx, Ny, 1),order='F',dtype=np.float32)), axis=2)
        h = h.reshape(-1)
        
        gradvec = np.vstack((f, g, h))
        
        gradnorm = np.linalg.norm(gradvec, axis=0)
        
        gamma = np.exp(-gradnorm ** 2 / (T ** 2))
        
        # Construct the matrix S
        for ll in range(np.size(gradnorm)):
            if gradnorm[ll] > 0:
                nu = gradvec[:, ll] / gradnorm[ll]
                B = np.eye(3) - (1. - gamma[ll]) * np.outer(nu, nu)
            else:
                B = np.eye(3)
            S[3 * ll:3 * ll + 3, :] = B
        return S
    if options.TV_use_anatomical:
        if isinstance(options.TV_referenceImage, str):
            if len(options.TV_referenceImage) == 0:
                raise ValueError('TV with anatomical weighting selected, but no reference image provided!')
            if options.TV_referenceImage[len(options.TV_referenceImage)-3:len(options.TV_referenceImage)+1:1] == 'mat':
                try:
                    from pymatreader import read_mat
                    apu = read_mat(options.TV_referenceImage)
                    options.TV_referenceImage = np.array(next(iter(apu.items()))[1])
                except ModuleNotFoundError:
                    print('pymatreader package not found! Mat-files cannot be loaded. You can install pymatreader package with "pip install pymatreader".')
            else:
                apu = np.load(options.TV_referenceImage, allow_pickle=True)
                variables = list(apu.keys())
                options.TV_referenceImage = apu[variables[0]]
        if options.TV_referenceImage.shape[1] == 1:
            koko_apu = np.sqrt(np.size(options.TV_referenceImage) / options.Nz[0].item())
            if np.floor(koko_apu) != koko_apu:
                raise ValueError('Reference image has to be square')
            else:
                options.TV_referenceImage = options.TV_referenceImage.reshape((koko_apu, koko_apu, options.Nz[0].item()))
                if koko_apu != options.Nx[0].item() or options.TV_referenceImage.shape[2] != options.Nz[0].item():
                    if options.Nz[0].item() > 1:
                        options.TV_referenceImage = resize(options.TV_referenceImage, (options.Nx[0].item(), options.Ny[0].item(), options.Nz[0].item()))
        else:
            if options.TV_referenceImage.shape[1] != options.Ny[0].item() or options.TV_referenceImage.shape[2] != options.Nz[0].item():
                if options.Nz[0].item() > 1:
                    options.TV_referenceImage = resize(options.TV_referenceImage, (options.Nx[0].item(), options.Ny[0].item(), options.Nz[0].item()))
        options.TV_referenceImage = options.TV_referenceImage.astype(dtype=np.float32)
        options.TV_referenceImage = options.TV_referenceImage - np.min(options.TV_referenceImage)
        options.TV_referenceImage = options.TV_referenceImage / np.max(options.TV_referenceImage)
        if options.TVtype == 1:
            options.TV_referenceImage = options.TV_referenceImage.reshape((koko_apu, koko_apu, options.Nz[0].item()),order='F')
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
            options.s = np.concatenate((s1.flatten(), s2.flatten(), s3.flatten(), s4.flatten(), s5.flatten(), s6.flatten(), s7.flatten(), s8.flatten(), s9.flatten()))
        options.TV_referenceImage = options.TV_referenceImage.ravel('F').astype(dtype=np.float32)
        
def APLSPrepass(options):
    from skimage.transform import resize #scikit-image
    if isinstance(options.APLS_ref_image, str):
        if len(options.APLS_ref_image) == 0:
            raise ValueError('APLS selected, but no reference image provided!')
        if options.APLS_ref_image[len(options.APLS_ref_image)-3:len(options.APLS_ref_image)+1:1] == 'mat':
            try:
                from pymatreader import read_mat
                apu = read_mat(options.APLS_ref_image)
                options.APLS_ref_image = np.array(next(iter(apu.items()))[1])
            except ModuleNotFoundError:
                print('pymatreader package not found! Mat-files cannot be loaded. You can install pymatreader package with "pip install pymatreader".')
        else:
            apu = np.load(options.APLS_ref_image, allow_pickle=True)
            variables = list(apu.keys())
            options.APLS_ref_image = apu[variables[0]]
    if options.APLS_ref_image.shape[1] == 1:
        koko_apu = np.sqrt(np.size(options.APLS_ref_image) / options.Nz[0])
        if koko_apu != np.floor(koko_apu):
            raise ValueError('Reference image has to be square')
        else:
            options.APLS_ref_image = options.APLS_ref_image.reshape((koko_apu, koko_apu, options.Nz[0]))
            if koko_apu != options.Nx[0] or options.APLS_ref_image.shape[2] != options.Nz[0]:
                options.APLS_ref_image = resize(options.APLS_ref_image, (options.Nx[0], options.Ny[0], options.Nz[0]))
    else:
        if options.APLS_ref_image.shape[1] != options.Ny or options.APLS_ref_image.shape[2] != options.Nz:
            options.APLS_ref_image = resize(options.APLS_ref_image, (options.Nx[0], options.Ny[0], options.Nz[0]))
    options.APLS_ref_image = options.APLS_ref_image.astype(dtype=np.float32)
    options.APLS_ref_image = options.APLS_ref_image.ravel('F')

def computeWeights(options, GGMRF):
    distX = options.FOVa_x[0] / options.Nx[0]
    distY = options.FOVa_y[0] / options.Ny[0]
    distZ = options.axial_fov[0] / options.Nz[0]
    
    if np.size(options.weights) == 0:
        options.weights = np.zeros(((options.Ndx * 2 + 1) * (options.Ndy * 2 + 1) * (options.Ndz * 2 + 1)),order='F',dtype=np.float32)
        cc = np.zeros((options.Ndy * 2 + 1) * (options.Ndx * 2 + 1),order='F',dtype=np.float32)
        lt = 0
        if GGMRF == True:
            for jj in range(options.Ndx, -options.Ndx-1, -1):
                lt += 1
                ll = 0
                for kk in range(options.Ndy, -options.Ndy-1, -1):
                    ll += 1
                    if options.Ndx == 0 or options.Nx[0] == 1:
                        apu = np.column_stack(((np.arange(options.Ndz, -options.Ndz-1, -1) * distZ), (np.repeat(kk, options.Ndy*2+1) * distY)))
                    else:
                        if options.Ndx != options.Ndz:
                            apu = np.column_stack(((np.arange(options.Ndz, -options.Ndz-1, -1) * distZ), (np.repeat(kk, options.Ndy*2+1) * distY), 
                                                   (np.concatenate((np.zeros(options.Ndz-options.Ndx), np.repeat(jj, options.Ndx*2+1) * distX, np.zeros(options.Ndx-options.Ndx))))))
                        else:
                            apu = np.column_stack(((np.arange(options.Ndz, -options.Ndz-1, -1) * distZ), (np.repeat(kk, options.Ndy*2+1) * distY), (np.repeat(jj, options.Ndx*2+1) * distX)))
                    edist = np.sqrt(np.sum(apu**2, axis=1))
                    cc[(options.Ndy*2+1)*(ll-1):(options.Ndy*2+1)*ll] = edist
                    options.weights[(options.Ndz*2+1)*(options.Ndy*2+1)*(lt-1):(options.Ndz*2+1)*(options.Ndy*2+1)*lt] = cc
        else:
            for jj in range(options.Ndz, -options.Ndz-1, -1):
                lt += 1
                ll = 0
                for kk in range(options.Ndy, -options.Ndy-1, -1):
                    ll += 1
                    if options.Ndz == 0 or options.Nz[0].item() == 1:
                        apu = np.column_stack(((np.arange(options.Ndx, -options.Ndx-1, -1) * distX), (np.repeat(kk, options.Ndy*2+1) * distY)))
                    else:
                        if options.Ndz != options.Ndx:
                            apu = np.column_stack(((np.arange(options.Ndx, -options.Ndx-1, -1) * distX), (np.repeat(kk, options.Ndy*2+1) * distY), 
                                                   (np.concatenate((np.zeros(options.Ndx-options.Ndz), np.repeat(jj, options.Ndz*2+1) * distZ, np.zeros(options.Ndx-options.Ndz))))))
                        else:
                            apu = np.column_stack(((np.arange(options.Ndx, -options.Ndx-1, -1) * distX), (np.repeat(kk, options.Ndy*2+1) * distY), (np.repeat(jj, options.Ndz*2+1) * distZ)))
                    edist = np.sqrt(np.sum(apu**2, axis=1))
                    cc[(options.Ndy*2+1)*(ll-1):(options.Ndy*2+1)*ll] = edist
                    options.weights[(options.Ndx*2+1)*(options.Ndy*2+1)*(lt-1):(options.Ndx*2+1)*(options.Ndy*2+1)*lt] = cc
        options.weights = 1.0 / options.weights
        options.weights = options.weights.astype(dtype=np.float32)
        
def quadWeights(options, isEmpty):
    if isEmpty:
        non_inf_weights_sum = np.nansum(options.weights)
        options.weights_quad = options.weights / non_inf_weights_sum
        if not options.GGMRF:
            half_len = np.size(options.weights_quad) // 2
            options.weights_quad = np.concatenate((options.weights_quad[:half_len], options.weights_quad[half_len + 1:]))
    else:
        options.weights_quad = options.weights
    if not options.GGMRF:
        options.weights_quad = options.weights_quad[~np.isinf(options.weights_quad)]
    options.weights_quad = options.weights_quad.astype(dtype=np.float32)
        
def huberWeights(options):
    if np.size(options.weights_huber) == 0:
        non_inf_weights_sum = np.nansum(options.weights)
        options.weights_huber = options.weights / non_inf_weights_sum
        half_len = np.size(options.weights_huber) // 2
        options.weights_huber = np.concatenate((options.weights_huber[:half_len], options.weights_huber[half_len + 1:]))
    options.weights_huber = options.weights_huber[~np.isinf(options.weights_huber)]
    options.weights_huber = options.weights_huber.astype(dtype=np.float32)
    
def weightedWeights(options):
    if np.size(options.weighted_weights) == 0:
        distX = options.FOVa_x / float(options.Nx[0].item())
        kerroin = np.sqrt(2.) * distX
        options.weighted_weights = kerroin * options.weights
        options.weighted_weights[np.isinf(options.weighted_weights)] = options.weighted_center_weight
        options.weighted_weights /= np.sum(options.weighted_weights)
    options.weighted_weights = np.reshape(options.weighted_weights, (options.Ndx * 2 + 1, options.Ndy * 2 + 1, options.Ndz * 2 + 1),order='F').astype(dtype=np.float32)

def NLMPrepass(options):
    def gaussianKernel(x, y, z, sigma_x, sigma_y, sigma_z = 0):
        gaussK = np.exp(-(np.add.outer(np.add.outer(x**2 / (2*sigma_x**2), y**2 / (2*sigma_y**2)), z**2 / (2*sigma_z**2))))
        return gaussK
    g_x = np.linspace(-options.Nlx, options.Nlx, 2 * options.Nlx + 1, dtype=np.float32)
    g_y = np.linspace(-options.Nly, options.Nly, 2 * options.Nly + 1, dtype=np.float32)
    g_z = np.linspace(-options.Nlz, options.Nlz, 2 * options.Nlz + 1, dtype=np.float32)
    gaussian = gaussianKernel(g_x, g_y, g_z, options.NLM_gauss, options.NLM_gauss, options.NLM_gauss)
    options.gaussianNLM = gaussian.flatten('F').astype(dtype=np.float32)
    if options.NLM_use_anatomical:
        if isinstance(options.NLM_referenceImage, str):
            if len(options.NLM_referenceImage) == 0:
                raise ValueError('NLM with anatomical weighting selected, but no reference image provided!')
            if options.NLM_referenceImage[len(options.NLM_referenceImage)-3:len(options.NLM_referenceImage)+1:1] == 'mat':
                try:
                    from pymatreader import read_mat
                    apu = read_mat(options.NLM_referenceImage)
                    options.NLM_referenceImage = np.array(next(iter(apu.items()))[1])
                except ModuleNotFoundError:
                    print('pymatreader package not found! Mat-files cannot be loaded. You can install pymatreader package with "pip install pymatreader".')
            else:
                apu = np.load(options.NLM_referenceImage, allow_pickle=True)
                variables = list(apu.keys())
                options.NLM_referenceImage = apu[variables[0]]
        options.NLM_referenceImage = options.NLM_referenceImage.ravel('F').astype(dtype=np.float32)

def prepassPhase(options):
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
        if isinstance(options.referenceImage, str):
            if options.referenceImage[len(options.referenceImage)-3:len(options.referenceImage)+1:1] == 'mat':
                try:
                    from pymatreader import read_mat
                    apu = read_mat(options.referenceImage)
                    options.referenceImage = np.array(next(iter(apu.items()))[1])
                except ModuleNotFoundError:
                    print('pymatreader package not found! Mat-files cannot be loaded. You can install pymatreader package with "pip install pymatreader".')
            else:
                apu = np.load(options.referenceImage, allow_pickle=True)
                if isinstance(list, apu):
                    variables = list(apu.keys())
                    options.referenceImage = apu[variables[0]]
                else:
                    options.referenceImage = apu;
        options.referenceImage = options.referenceImage.ravel('F').astype(dtype=np.float32)
        if np.size(options.referenceImage) == round((options.NxFull - options.NxOrig) * options.multiResolutionScale) * \
            round((options.NyFull - options.NyOrig) * options.multiResolutionScale) * \
            round((options.NzFull - options.NzOrig) * options.multiResolutionScale):
            skip = True
        else:
            skip = False
        
        if not skip and np.size(options.referenceImage) != options.NxFull * options.NyFull * options.NzFull:
            raise ValueError('The size of the reference image does not match the reconstructed image!')
        
        if not skip and options.nMultiVolumes > 0:
            options.referenceImage = options.referenceImage.reshape(options.NxFull, options.NyFull, options.NzFull, order = 'F')
            
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
        or options.RAMLA or options.MBSREM or options.MRAMLA or options.ROSEM or options.DRAMA or options.ROSEM_MAP or options.ECOSEM
        or options.COSEM or options.ACOSEM or options.AD or np.any(options.OSL_COSEM) or options.NLM or options.OSL_RBI or options.RBI or options.PKMA
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
        if (options.BSREM or options.RAMLA or options.MBSREM or options.MRAMLA or options.ROSEM_MAP or options.ROSEM or options.PKMA or options.SPS) and np.size(options.lambdaN) == 0:
            lambda_vals = np.zeros(options.Niter, dtype=np.float32)
            for i in range(options.Niter):
                lambda_vals[i] = 1. / (i / 20. + 1.)
            options.lambdaN = lambda_vals
            if options.CT == True:
                options.lambdaN = options.lambdaN / 10000.
        elif (options.BSREM or options.RAMLA or options.MBSREM or options.MRAMLA or options.ROSEM_MAP or options.ROSEM or options.param.PKMA or options.param.SPS):
            if np.size(options.lambdaN) < options.Niter:
                raise ValueError('The number of relaxation values needs to be at least equal to the number of iterations!')
            elif np.size(options.lambdaN) > options.Niter:
                print('Warning: The number of relaxation values is more than the number of iterations. Later values are ignored!')
    
        if options.DRAMA:
            options.lam_drama = np.zeros((options.Niter, options.subsets),order='F',dtype=np.float32)
            options.lam_drama[0, 0] = options.beta_drama / (options.alpha_drama * options.beta0_drama)
            r = 1
            for i in range(options.Niter):
                for j in range(options.subsets):
                    options.lam_drama[i, j] = options.beta_drama / (options.alpha_drama * options.beta0_drama + r)
                    r += 1
            
        if options.PKMA and np.size(options.alpha_PKMA) < options.Niter * options.subsets:
            if np.size(options.alpha_PKMA) < options.Niter * options.subsets:
                print('Warning: The number of PKMA alpha (momentum) values must be at least the number of iterations times the number of subsets! Computing custom alpha values.')
            options.alpha_PKMA = np.zeros(options.Niter * options.subsets, dtype=np.float32)
            oo = 0
            for kk in range(options.Niter):
                for ll in range(options.subsets):
                    options.alpha_PKMA[oo] = 1. + (options.rho_PKMA * (kk * options.subsets + ll)) / (kk * options.subsets + ll + options.delta_PKMA)
                    oo += 1
            # if options.CT == True:
            #     options.alpha_PKMA = options.alpha_PKMA / 10000.
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
            if options.quad or (options.TV and options.TVtype == 3) or options.GGMRF or options.hyperbolic:
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
                    options.inffi = np.floor(options.weights.size / 2)[0]
            if options.weighted_mean:
                weightedWeights(options)
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
            oo = 0
            for kk in range(1, options.Niter + 1):
                for ll in range(options.subsets):
                    options.thetaCP[oo] = 1. + (options.rho_PKMA * ((kk - 1) * options.subsets + ll)) / ((kk - 1) * options.subsets + ll + options.delta_PKMA)
                    oo += 1
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
        options.Nf = 2 ** np.ceil(np.log2(options.Nx[0].item()))
        options.filterIm = rampFilt(options.Nf, options.filterWindow, options.cutoffFrequency, options.normalFilterSigma, True)
        options.filterIm = options.filterIm.astype(dtype=np.float32)
    
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
    if isinstance(options.TV_referenceImage, str):
        options.TV_referenceImage = np.empty(0, dtype=np.float32)