# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 13:17:22 2025
"""

def initProjector(self):
    import arrayfire as af
    import numpy as np
    import os
    from omegatomo.reconstruction.prepass import prepassPhase
    from omegatomo.reconstruction.prepass import parseInputs
    from omegatomo.reconstruction.prepass import loadCorrections
    if self.useAF:
        # import arrayfire as af
        if af.get_active_backend() != 'opencl' and not self.useCUDA:
            af.set_backend('opencl')

        af.device.set_device(self.deviceNum)
    if self.useTorch and self.useAF:
        raise ValueError('Arrayfire and PyTorch cannot be used at the same time! Select only one!')
    if self.useTorch:
        import torch
        torch.cuda.init()
    if self.useTorch and not self.useCUDA:
        raise ValueError('PyTorch does not work with OpenCL! You can still use OpenCL manually with PyTorch, but you have to manually transfer the OpenCL data first to host (NumPy array) and then to Torch (or vice versa, i.e. Torch --> NumPy --> OpenCL)')
    if self.useCuPy and self.useCUDA:
        import cupy as cp
    
    self.NVOXELS = 8
    self.TH = 100000000000.
    self.TH32 = 100000.
    self.NVOXELS5 = 1
    self.NVOXELSFP = 8
    if np.size(self.weights) > 0:
        self.empty_weight = False
    if self.TOF_bins_used == 0:
        self.TOF_bins_used = 1
    mDataFound = self.SinM.size > 0
    loadCorrections(self)
    parseInputs(self, mDataFound)
    prepassPhase(self)
    if self.listmode > 0 and self.subsets > 1 and self.subsetType > 0:
        if self.useIndexBasedReconstruction:
            self.trIndex = self.trIndex[:,self.index]
            self.axIndex = self.axIndex[:,self.index]
        # else:
        #     self.x = np.reshape(self.x, [6, -1], order='F')
        #     self.x = self.x[:,self.index]
        #     self.x = self.x.ravel('F')
    if self.useIndexBasedReconstruction and self.listmode > 0:
        self.trIndex = self.trIndex.ravel('F')
        self.axIndex = self.axIndex.ravel('F')

    if self.projector_type in [1, 11, 14, 15, 12, 13]:
        self.FPType = 1
    elif self.projector_type in [2, 21, 22, 23, 24, 25]:
        self.FPType = 2
    elif self.projector_type in [3, 31, 32, 33, 34, 35]:
        self.FPType = 3
    elif self.projector_type in [4, 41, 42, 43, 44, 45]:
        self.FPType = 4
    elif self.projector_type in [5, 51, 52, 53, 54, 55]:
        self.FPType = 5
    elif self.projector_type == 6:
        self.FPType = 6
    else:
        raise ValueError('Invalid forward projector!')
    if self.projector_type in [1, 11, 21, 31, 41, 51]:
        self.BPType = 1
    elif self.projector_type in [2, 12, 22, 32, 42, 52]:
        self.BPType = 2
    elif self.projector_type in [3, 13, 23, 33, 43, 53]:
        self.BPType = 3
    elif self.projector_type in [4, 14, 24, 34, 44, 54]:
        self.BPType = 4
    elif self.projector_type in [5, 15, 25, 35, 45, 55]:
        self.BPType = 5
    elif self.projector_type == 6:
        self.BPType = 6
    else:
        raise ValueError('Invalid backprojector!')
    # if self.useAF == False and (self.FPType == 5 or self.BPType == 5):
    #     raise ValueError('Branchless distance-driven (projector type 5) can only be used with Arrayfire!')
    if (self.useAF == False and self.useCuPy == False) and self.projector_type == 6:
        raise ValueError('Projector type 6 can only be used with Arrayfire (OpenCL) or CuPy (CUDA)!')
        
    if not self.projector_type == 6:
        headerDir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', '..', '..', 'opencl')) + "/"
        with open(headerDir + 'general_opencl_functions.h', encoding="utf8") as f:
            hlines = f.read()
        if self.FPType in [1, 2, 3]:
            with open(headerDir + 'projectorType123.cl', encoding="utf8") as f:
                linesFP = f.read()
        elif self.FPType in [4]:
            with open(headerDir + 'projectorType4.cl', encoding="utf8") as f:
                linesFP = f.read()
        elif self.FPType in [5]:
            with open(headerDir + 'projectorType5.cl', encoding="utf8") as f:
                linesFP = f.read()
        if self.BPType in [1, 2, 3]:
            with open(headerDir + 'projectorType123.cl', encoding="utf8") as f:
                linesBP = f.read()
        elif self.BPType in [4]:
            with open(headerDir + 'projectorType4.cl', encoding="utf8") as f:
                linesBP = f.read()
        elif self.BPType in [5]:
            with open(headerDir + 'projectorType5.cl', encoding="utf8") as f:
                linesBP = f.read()
        globalSize = [None] * self.subsets
        # self.mSize = [None] * self.subsets
        for i in range(self.subsets):
            if (self.FPType == 5):
                globalSize[i] = (self.nRowsD, (self.nColsD + self.NVOXELSFP - 1) // self.NVOXELSFP, self.nProjSubset[i].item())
                localSize = (16, 16, 1)
                erotus = (localSize[0] - (globalSize[i][0] % localSize[0]), localSize[1] - (globalSize[i][1] % localSize[1]), 0)
                globalSize[i] = (self.nRowsD + erotus[0], (self.nColsD + self.NVOXELSFP - 1) // self.NVOXELSFP + erotus[1], self.nProjSubset[i].item())
            elif ((self.CT or self.SPECT or self.PET) and self.listmode == 0):
                globalSize[i] = (self.nRowsD, self.nColsD, self.nProjSubset[i].item())
                localSize = (16, 16, 1)
                erotus = (localSize[0] - (globalSize[i][0] % localSize[0]), localSize[1] - (globalSize[i][1] % localSize[1]), 0)
                globalSize[i] = (self.nRowsD + erotus[0], self.nColsD + erotus[1], self.nProjSubset[i].item())
            else:
                globalSize[i] = (self.nMeasSubset[i].item(), 1, 1)
                localSize = (128, 1, 1)
                erotus = (localSize[0] - (globalSize[i][0] % localSize[0]), localSize[1] - (globalSize[i][1] % localSize[1]), 0)
                globalSize[i] = (self.nMeasSubset[i].item() + erotus[0], 1, 1)
        self.globalSizeFP = globalSize.copy()
        self.localSizeFP = localSize + tuple()
        self.erotusBP = [0] * (self.nMultiVolumes + 1) * 2
        localSize = (16, 16, 1)
        for ii in range(self.nMultiVolumes + 1):
            apu = [self.Nx[ii].item() % localSize[0], self.Ny[ii].item() % localSize[1], 0]
            if apu[0] > 0:
                self.erotusBP[ii * 2] = localSize[0] - apu[0]
            if apu[1] > 0:
                self.erotusBP[ii * 2 + 1] = localSize[1] - apu[1]
        
        if self.BPType in [1, 2, 3] or (self.BPType == 4 and not self.CT):
            globalSize = [[None] * (self.nMultiVolumes + 1)] * self.subsets
            for i in range(self.subsets):
                for ii in range(self.nMultiVolumes + 1):
                    globalSize[i][ii] = self.globalSizeFP[i]
            self.localSizeBP = self.localSizeFP + tuple()
        else:
            globalSize = [[None] * (self.nMultiVolumes + 1)] * self.subsets
            for i in range(self.subsets):
                for ii in range(self.nMultiVolumes + 1):
                    if self.BPType == 4:
                        globalSize[i][ii] = (self.Nx[ii].item() + self.erotusBP[ii * 2], self.Ny[ii].item() + self.erotusBP[ii * 2 + 1], (self.Nz[ii].item()  + self.NVOXELS - 1) // self.NVOXELS)
                    elif self.BPType == 5:
                        if self.pitch:
                            globalSize[i][ii] = (self.Nx[ii].item() + self.erotusBP[ii * 2], self.Ny[ii].item() + self.erotusBP[ii * 2 + 1], self.Nz[ii].item())
                        else:
                            globalSize[i][ii] = (self.Nx[ii].item() + self.erotusBP[ii * 2], self.Ny[ii].item() + self.erotusBP[ii * 2 + 1], (self.Nz[ii].item()  + self.NVOXELS5 - 1) // self.NVOXELS5)
                    else:
                        globalSize[i][ii] = (self.Nx[ii].item() + self.erotusBP[ii * 2], self.Ny[ii].item() + self.erotusBP[ii * 2 + 1], self.Nz[ii].item())
            self.localSizeBP = localSize + tuple()
        self.globalSizeBP = globalSize.copy()
                            
        
        self.Nxy = self.Nx[0].item() * self.Ny[0].item()
        if self.useCUDA:
            if self.useImages and not self.useCuPy:
                self.useImages = False
            if self.use_64bit_atomics or self.use_32bit_atomics:
                self.use_64bit_atomics = False
                self.use_32bit_atomics = False
            bOpt = ('-DCUDA','-DPYTHON',)
        else:
            bOpt =('-cl-single-precision-constant -DOPENCL',)
        if self.useMAD:
            if self.useCUDA:
                bOpt += ('--use_fast_math','-DUSEMAD',)
            else:
                bOpt += (' -cl-fast-relaxed-math -DUSEMAD',)
        if self.useImages:
            bOpt += ('-DUSEIMAGES',)
        if (self.FPType == 2 or self.BPType == 2 or self.FPType == 3 or self.BPType == 3):
            if self.orthTransaxial:
                bOpt += ('-DCRYSTXY',)
            if self.orthAxial:
                bOpt += ('-DCRYSTZ',)
            with open(headerDir + 'opencl_functions_orth3D.h') as f:
                hlines2 = f.read()
            if self.FPType in [2, 3]:
                linesFP = hlines + hlines2 + linesFP
            else:
                linesFP = hlines + linesFP
            if self.BPType in [2, 3]:
                linesBP = hlines + hlines2 + linesBP
            else:
                linesBP = hlines + linesBP
        else:
            linesFP = hlines + linesFP
            linesBP = hlines + linesBP
        if self.FPType == 3 or self.BPType == 3:
            bOpt += ('-DVOL',)
        if self.useMaskFP:
            bOpt += ('-DMASKFP',)
        if self.useMaskBP:
            bOpt += ('-DMASKBP',)
        if self.useTotLength and not self.SPECT:
            bOpt += ('-DTOTLENGTH',)
        if self.OffsetLimit.size > 0:
            bOpt += ('-DOFFSET',)
        if self.attenuation_correction and self.CTAttenuation:
            bOpt += ('-DATN',)
        elif self.attenuation_correction and not self.CTAttenuation:
            bOpt += ('-DATNM',)
        if self.normalization_correction:
            bOpt += ('-DNORM',)
        if self.additionalCorrection:
            bOpt += ('-DSCATTER',)
        if self.randoms_correction:
            bOpt += ('-DRANDOMS',)
        if self.nLayers > 1:
            if self.useIndexBasedReconstruction:
                bOpt += ('-DNLAYERS=' + str(self.nLayers),)
            else:
                bOpt += ('-DNLAYERS=' + str(self.nProjections // (self.nLayers * self.nLayers)),)
        if self.TOF:
            bOpt += ('-DTOF',)
        if self.CT:
            bOpt += ('-DCT',)
        elif self.SPECT:
            if self.useCUDA:
                if self.useCuPy:
                    bOpt += ('-DSPECT',)
                else:
                    bOpt += ('-DSPECT', )
            else:
                bOpt += (' -DSPECT',)
        elif self.PET:
            bOpt += ('-DPET',)

        bOpt += ('-DNBINS=' + str(self.TOF_bins_used),)
        if self.listmode:
            bOpt += ('-DLISTMODE',)
        if self.listmode > 0 and self.useIndexBasedReconstruction:
            bOpt += ('-DINDEXBASED',)
        if self.listmode > 0 and ~self.useIndexBasedReconstruction:
            bOpt += ('-DUSEGLOBAL',)
        if (self.FPType == 1 or self.BPType == 1 or self.FPType == 4 or self.BPType == 4) and self.n_rays_transaxial * self.n_rays_axial > 1:
            bOpt += ('-DN_RAYS=' + str(self.n_rays_transaxial * self.n_rays_axial),)
            bOpt += ('-DN_RAYS2D=' + str(self.n_rays_transaxial),)
            bOpt += ('-DN_RAYS3D=' + str(self.n_rays_axial),)
        if self.pitch:
            bOpt += ('-DPITCH',)
        if (((self.subsets > 1 and (self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7))) and not self.CT and not self.SPECT and not self.PET and self.listmode == 0):
            bOpt += ('-DSUBSETS',)
        if self.subsets > 1 and self.listmode == 0:
            bOpt += ('-DSTYPE=' + str(self.subsetType),'-DNSUBSETS=' + str(self.subsets),)
        
        bOptFP = bOpt + ('-DFP',)
        if self.localSizeFP[1] > 1:
            bOptFP += ('-DLOCAL_SIZE=' + str(self.localSizeFP[0]),'-DLOCAL_SIZE2=' + str(self.localSizeFP[1]),)
        else:
            bOptFP += ('-DLOCAL_SIZE=' + str(self.localSizeFP[0]),'-DLOCAL_SIZE2=' + str(1),)
        if self.FPType in [1, 2, 3]:
            bOptFP += ('-DSIDDON',)
            if self.FPType in [2, 3]:
                bOptFP += ('-DORTH',)
            if self.FPType == 3:
                bOptFP += ('-DVOL',)
            if self.use_64bit_atomics:
                bOptFP += ('-DCAST=long',)
            elif self.use_32bit_atomics:
                bOptFP += ('-DCAST=int',)
            else:
                bOptFP += ('-DCAST=float',)
        elif self.FPType == 4:
            bOptFP += ('-DPTYPE4','-DNVOXELS=' + str(self.NVOXELS),)
            if not self.CT:
                if self.use_64bit_atomics:
                    bOptFP += ('-DCAST=long',)
                elif self.use_32bit_atomics:
                    bOptFP += ('-DCAST=int',)
                else:
                    bOptFP += ('-DCAST=float',)
        elif self.FPType == 5:
            bOptFP += ('-DPROJ5','-DNVOXELSFP=' + str(self.NVOXELSFP),)
            if self.meanFP:
                bOptFP += ('-DMEANDISTANCEFP',)
        
        bOptBP = bOpt + ('-DBP',)
        if self.localSizeBP[1] > 1:
            bOptBP += ('-DLOCAL_SIZE=' + str(self.localSizeBP[0]),'-DLOCAL_SIZE2=' + str(self.localSizeBP[1]),)
        else:
            bOptBP += ('-DLOCAL_SIZE=' + str(self.localSizeBP[0]),'-DLOCAL_SIZE2=' + str(1),)
        if self.BPType in [1, 2, 3]:
            bOptBP += ('-DSIDDON',)
            if self.BPType in [2, 3]:
                bOptBP += ('-DORTH',)
            if self.BPType == 3:
                bOptBP += ('-DVOL',)
            bOptBP += ('-DATOMICF',)
            if self.use_64bit_atomics:
                bOptBP += ('-DATOMIC','-DCAST=long','-DTH=' + str(self.TH),)
            elif self.use_32bit_atomics:
                bOptBP += (' -DATOMIC32',' -DCAST=int',' -DTH=' + str(self.TH32),)
            else:
                bOptBP += ('-DCAST=float',)
        elif self.BPType == 4:
            bOptBP += ('-DPTYPE4','-DNVOXELS=' + str(self.NVOXELS),)
            if not self.CT:
                bOptBP += ('-DATOMICF',)
                if self.use_64bit_atomics:
                    bOptBP += ('-DATOMIC','-DCAST=long','-DTH=' + str(self.TH),)
                elif self.use_32bit_atomics:
                    bOptBP += (' -DATOMIC32',' -DCAST=int',' -DTH=' + str(self.TH32),)
                else:
                    bOptBP += ('-DCAST=float',)
        elif self.BPType == 5:
            bOptBP += ('-DPROJ5','-DNVOXELS5=' + str(self.NVOXELS5),)
            if self.meanBP:
                bOptBP += ('-DMEANDISTANCEBP',)
    else:
        if self.useCUDA:
            if self.useTorch:
                self.gFilter = np.ascontiguousarray(self.gFilter)
                self.gFilter = np.transpose(self.gFilter, (1, 0, 2, 3))
                self.d_gFilter = torch.tensor(self.gFilter, device='cuda')
                self.angles = np.degrees(self.angles)
                self.d_gFilter = self.d_gFilter.permute(2, 0, 1, 3).unsqueeze(1)
        else:
            self.d_gFilter = af.interop.np_to_af_array(self.gFilter)
        self.uu = 0
    
    if self.useCUDA:
        if (self.BPType == 5 or self.FPType == 5 or self.FPType == 4) and not self.useCuPy:
            raise ValueError('Unsupported projector for CUDA! Only projector types 1, 11, 12, 13, 2, 21, 23, 3, 31, 32, 14 and 6 are supported')
        self.no_norm = 1
        self.mSize = self.nRowsD * self.nColsD * self.nProjections
        self.d_d = [None] * (self.nMultiVolumes + 1)
        self.d_b = [None] * (self.nMultiVolumes + 1)
        self.d_bmax = [None] * (self.nMultiVolumes + 1)
        self.d_Nxyz = [None] * (self.nMultiVolumes + 1)
        self.dSize = [None] * (self.nMultiVolumes + 1)
        self.d_Scale = [None] * (self.nMultiVolumes + 1)
        self.d_Scale4 = [None] * (self.nMultiVolumes + 1)
        self.d_x = [None] * self.subsets
        self.d_z = [None] * self.subsets
        if self.projector_type != 6:
            if self.useCuPy:
                # if self.FPType == 5:
                #     raise ValueError('Not yet supported')
                self.d_Sens = cp.empty(shape=(1,1), dtype=cp.float32)
                if (self.listmode == 0 and not (self.CT or self.SPECT)) or self.useIndexBasedReconstruction:
                    self.d_x[0] = cp.asarray(self.x.ravel())
                elif (self.CT or self.SPECT) and self.listmode == 0:
                    apu = self.x.ravel()
                    for i in range(self.subsets):
                        self.d_x[i] = cp.asarray(apu[self.nMeas[i] * 6 : self.nMeas[i + 1] * 6])
                elif self.listamode > 0 and not self.useIndexBasedReconstruction:
                    apu = self.x.ravel()
                    for i in range(self.subsets):
                        if self.loadTOF:
                            self.d_x[i] = cp.asarray(apu[self.nMeas[i] * 6 : self.nMeas[i + 1] * 6])
                if ((self.CT or self.SPECT) and self.listmode == 0):
                    if self.pitch:
                        kerroin = 6
                    else:
                        kerroin = 2
                    apu = self.z.ravel()
                    for i in range(self.subsets):
                        self.d_z[i] = cp.asarray(apu[self.nMeas[i] * kerroin : self.nMeas[i + 1] * kerroin])
                else:
                    if (self.PET and self.listmode == 0):
                        if self.nLayers > 1:
                            kerroin = 3
                        else:
                            kerroin = 2
                        apu = self.z.ravel()
                        for i in range(self.subsets):
                            self.d_z[i] = cp.asarray(apu[self.nMeas[i] * kerroin : self.nMeas[i + 1] * kerroin])
                    elif self.listmode == 0 or (self.listmode > 0 and self.useIndexBasedReconstruction):
                        self.d_z[0] = cp.asarray(self.z.ravel())
                    else:
                        for i in range(self.subsets):
                            self.d_z[i] = cp.asarray(np.zeros(1,dtype=np.float32))
                if (self.attenuation_correction and not self.CTAttenuation):
                    self.d_atten = [None] * self.subsets
                    for i in range(self.subsets):
                        self.d_atten[i] = cp.asarray(self.vaimennus[self.nTotMeas[i].item() : self.nTotMeas[i + 1].item()])
                elif (self.attenuation_correction and self.CTAttenuation):
                    if not self.useImages:
                        self.d_atten = cp.asarray(self.vaimennus)
                    else:
                        chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
                        array = cp.cuda.texture.CUDAarray(chl, self.Nx[0].item(), self.Ny[0].item(), self.Nz[0].item())
                        array.copy_from(self.vaimennus.reshape((self.Nz[0].item(), self.Ny[0].item(), self.Nx[0].item())))
                        res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                        tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=0)
                        self.d_atten = cp.cuda.texture.TextureObject(res, tdes)
                if self.useMaskFP:
                    if not self.useImages:
                        self.d_maskFP = cp.asarray(self.maskFP)
                    else:
                        chl = cp.cuda.texture.ChannelFormatDescriptor(8,0,0,0, cp.cuda.runtime.cudaChannelFormatKindUnsigned)
                        array = cp.cuda.texture.CUDAarray(chl, self.nRowsD, self.nColsD)
                        array.copy_from(self.maskFP.reshape((self.nColsD, self.nRowsD)))
                        res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                        tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=0)
                        self.d_maskFP = cp.cuda.texture.TextureObject(res, tdes)
                if self.useMaskBP:
                    if not self.useImages:
                        self.d_maskBP = cp.asarray(self.maskBP)
                    else:
                        chl = cp.cuda.texture.ChannelFormatDescriptor(8,0,0,0, cp.cuda.runtime.cudaChannelFormatKindUnsigned)
                        array = cp.cuda.texture.CUDAarray(chl, self.Nx[0].item(), self.Ny[0].item())
                        array.copy_from(self.maskFP.reshape((self.Ny[0].item(), self.Nx[0].item())))
                        res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                        tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=0)
                        self.d_maskFP = cp.cuda.texture.TextureObject(res, tdes)
                if self.TOF:
                    self.d_TOFCenter = cp.asarray(self.TOFCenter)
                if self.SPECT:
                    self.d_rayShiftsDetector = cp.asarray(self.rayShiftsDetector)
                    self.d_rayShiftsSource = cp.asarray(self.rayShiftsSource)
                if (self.BPType == 2 or self.BPType == 3 or self.FPType == 2 or self.FPType == 3):
                    self.d_V = cp.asarray(self.V)
                if (self.normalization_correction):
                    self.d_norm = [None] * self.subsets
                    for i in range(self.subsets):
                        self.d_norm[i] = cp.asarray(self.normalization[self.nTotMeas[i].item() : self.nTotMeas[i + 1].item()])
                if (self.additionalCorrection):
                    self.d_corr = [None] * self.subsets
                    for i in range(self.subsets):
                        self.d_corr[i] = cp.asarray(self.corrVector[self.nTotMeas[i].item() : self.nTotMeas[i + 1].item()])
                if (self.listmode != 1 and ((not self.CT and not self.SPECT and not self.PET) and (self.subsets > 1 and (self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7)))):
                    self.d_zindex = [None] * self.subsets
                    self.d_xyindex = [None] * self.subsets
                    for i in range(self.subsets):
                        self.d_xyindex[i] = cp.asarray(self.xy_index[self.nMeas[i] : self.nMeas[i + 1]])
                        self.d_zindex[i] = cp.asarray(self.z_index[self.nMeas[i] : self.nMeas[i + 1]])
                if (self.listmode > 0 and self.useIndexBasedReconstruction):
                    self.d_trIndex = [None] * self.subsets
                    self.d_axIndex = [None] * self.subsets
                    for i in range(self.subsets):
                        if self.loadTOF:
                            self.d_trIndex[i] = cp.asarray(self.trIndex[self.nMeas[i] * 2 : self.nMeas[i + 1] * 2])
                            self.d_axIndex[i] = cp.asarray(self.axIndex[self.nMeas[i] * 2 : self.nMeas[i + 1] * 2])
                if self.OffsetLimit.size > 0 and ((self.BPType == 4 and self.CT) or self.BPType == 5):
                    self.d_T = [None] * self.subsets
                    for i in range(self.subsets):
                        self.d_T[i] = cp.asarray(self.offsetLimit[self.nMeas[i].item() : self.nMeas[i + 1].item()])
                mod = cp.RawModule(code=linesFP, options=bOptFP)
                # import sys
                # mod.compile(log_stream=sys.stdout)
                if self.FPType in [1, 2, 3]:
                    self.knlF = mod.get_function('projectorType123')
                elif self.FPType == 4:
                    self.knlF = mod.get_function('projectorType4Forward')
                elif self.FPType == 5:
                    self.knlF = mod.get_function('projectorType5Forward')
                mod = cp.RawModule(code=linesBP, options=bOptBP)
                if self.BPType in [1, 2, 3]:
                    self.knlB = mod.get_function('projectorType123')
                elif self.BPType == 4 and not self.CT:
                    self.knlB = mod.get_function('projectorType4Forward')
                elif self.BPType == 4 and self.CT:
                    self.knlB = mod.get_function('projectorType4Backward')
                elif self.BPType == 5:
                    self.knlB = mod.get_function('projectorType5Backward')
                
                if self.use_psf:
                    with open(headerDir + 'auxKernels.cl', encoding="utf8") as f:
                        lines = f.read()
                    lines = hlines + lines
                    bOpt += ('-DCAST=float','-DPSF','-DLOCAL_SIZE=' + str(localSize[0]),'-DLOCAL_SIZE2=' + str(localSize[1]),)
                    mod = cp.RawModule(code=lines, options=bOpt)
                    self.knlPSF = mod.get_function('Convolution3D_f')
                    self.d_gaussPSF = cp.asarray(self.gaussK.ravel('F'))
                    
                if self.FPType in [1, 2, 3]:
                    self.kIndF = (cp.float32(self.global_factor), cp.float32(self.epps), cp.uint32(self.nRowsD), cp.uint32(self.det_per_ring), cp.float32(self.sigma_x),)
                    if self.SPECT:
                        self.kIndF += (self.d_rayShiftsDetector, self.d_rayShiftsSource, cp.float32(self.coneOfResponseStdCoeffA), cp.float32(self.coneOfResponseStdCoeffB), cp.float32(self.coneOfResponseStdCoeffC),)
                    self.kIndF += (cp.float32(self.dPitchX),cp.float32(self.dPitchY),)
                elif self.FPType == 4:
                    self.kIndF = (cp.uint32(self.nRowsD), cp.uint32(self.nColsD), cp.float32(self.dPitchX),cp.float32(self.dPitchY),cp.float32(self.dL),cp.float32(self.global_factor),)
                elif self.FPType == 5:
                    self.kIndF = (cp.uint32(self.nRowsD), cp.uint32(self.nColsD), cp.float32(self.dPitchX),cp.float32(self.dPitchY),)
                if self.FPType in [2,3]:
                    if self.FPType == 2:
                        self.kIndF += (cp.float32(self.tube_width_z),)
                    else:
                        self.kIndF += (cp.float32(self.tube_radius),)
                    self.kIndF += (cp.float32(self.bmin), cp.float32(self.bmax), cp.float32(self.Vmax),)
                if self.useMaskFP:
                    self.kIndF += (self.d_maskFP,)
                if self.FPType in [1, 2, 3]:
                    if self.TOF:
                        self.kIndF += (self.d_TOFCenter, )
                    if self.FPType in [2, 3]:
                        self.kIndF += (self.d_V, )
                    self.kIndF += (cp.uint32(self.nColsD),)
                if self.FPType == 4 and self.TOF:
                    self.kIndF += (self.d_TOFCenter, )
                    self.kIndF += (cp.float32(self.sigma_x), )
                if self.attenuation_correction and self.CTAttenuation and self.FPType in [1, 2, 3, 4]:
                    self.kIndF += (self.d_atten,)
                    
                
                if self.BPType == 4 or self.BPType == 5:
                    self.kIndB = (cp.uint32(self.nRowsD), cp.uint32(self.nColsD), cp.float32(self.dPitchX),cp.float32(self.dPitchY),)
                if self.BPType == 4 and not self.CT:
                    self.kIndB += (cp.float32(self.dL),)
                    self.kIndB += (cp.float32(self.global_factor),)
                if self.BPType in [1, 2, 3]:
                    self.kIndB = (cp.float32(self.global_factor), cp.float32(self.epps), cp.uint32(self.nRowsD), cp.uint32(self.det_per_ring), cp.float32(self.sigma_x),)
                    if self.SPECT:
                        self.kIndB += (self.d_rayShiftsDetector, self.d_rayShiftsSource, cp.float32(self.coneOfResponseStdCoeffA), cp.float32(self.coneOfResponseStdCoeffB), cp.float32(self.coneOfResponseStdCoeffC),)
                    self.kIndB += (cp.float32(self.dPitchX),cp.float32(self.dPitchY),)
                    if self.BPType in [2, 3]:
                        if self.BPType == 2:
                            self.kIndB  += (cp.float32(self.tube_width_z),)
                        else:
                            self.kIndB  += (cp.float32(self.tube_radius),)
                        self.kIndB += (cp.float32(self.bmin),)
                        self.kIndB += (cp.float32(self.bmax),)
                        self.kIndB += (cp.float32(self.Vmax),)
                    if self.useMaskFP:
                        self.kIndB += (self.d_maskFP,)
                if self.useMaskBP:
                    self.kIndB += (self.d_maskBP,)
                if self.BPType in [1, 2, 3]:
                    if self.TOF:
                        self.kIndB += (self.d_TOFCenter,)
                    if self.BPType in [2, 3]:
                        self.kIndB += (self.d_V,)
                    self.kIndB += (cp.uint32(self.nColsD),)
                if self.BPType == 4 and not self.CT and self.TOF:
                    self.kIndB += (self.d_TOFCenter,)
                    self.kIndB += (cp.float32(self.sigma_x),)
                if self.attenuation_correction and self.CTAttenuation and self.BPType in [1, 2, 3, 4] and not self.CT:
                    self.kIndB += (self.d_atten,)
            else:
                import pycuda as cuda
                from pycuda.compiler import SourceModule
                import pycuda.autoinit
                self.d_Sens = cuda.gpuarray.empty(shape=(1,1), dtype=np.float32)
                for k in range(self.nMultiVolumes + 1):
                    self.d_d[k] = cuda.gpuarray.vec.make_float3(self.dx[k].item(), self.dy[k].item(), self.dz[k].item())
                    self.d_b[k] = cuda.gpuarray.vec.make_float3(self.bx[k].item(), self.by[k].item(), self.bz[k].item())
                    self.d_bmax[k] = cuda.gpuarray.vec.make_float3(self.bx[k].item() + self.Nx[k].item() * self.dx[k].item(), self.by[k].item() + self.Ny[k].item() * self.dy[k].item(), self.bz[k].item() + self.Nz[k].item() * self.dz[k].item())
                    self.d_Nxyz[k] = cuda.gpuarray.vec.make_uint3(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item())
                    if (self.FPType == 4 or self.FPType == 5 or self.BPType == 4 or self.BPType == 5):
                        self.d_Scale4[k] = cuda.gpuarray.vec.make_float3(self.dScaleX4[k].item(), self.dScaleY4[k].item(), self.dScaleZ4[k].item())
                        if self.FPType == 5 or self.BPType == 5:
                            self.dSize[k] = cuda.gpuarray.vec.make_float2(self.dSizeX[k].item(), self.dSizeY[k].item())
                            self.d_Scale[k] = cuda.gpuarray.vec.make_float3(self.dScaleX[k].item(), self.dScaleY[k].item(), self.dScaleZ[k].item())
                            if k == 0:
                                self.dSizeBP = cuda.gpuarray.vec.make_float2(self.dSizeXBP, self.dSizeZBP)
                self.d_dPitch = cuda.gpuarray.vec.make_float2(self.dPitchX, self.dPitchY)
                if (self.listmode == 0 and not (self.CT or self.SPECT)) or self.useIndexBasedReconstruction:
                    self.d_x[0] = cuda.gpuarray.to_gpu(self.x.ravel())
                elif (self.CT or self.SPECT) and self.listmode == 0:
                    apu = self.x.ravel()
                    for i in range(self.subsets):
                        self.d_x[i] = cuda.gpuarray.to_gpu(apu[self.nMeas[i] * 6 : self.nMeas[i + 1] * 6])
                elif self.listamode > 0 and not self.useIndexBasedReconstruction:
                    apu = self.x.ravel()
                    for i in range(self.subsets):
                        if (i == 0 or self.loadTOF):
                            self.d_x[i] = cuda.gpuarray.to_gpu(apu[self.nMeas[i] * 6 : self.nMeas[i + 1] * 6])
                if ((self.CT or self.SPECT) and self.listmode == 0):
                    if self.pitch:
                        kerroin = 6
                    else:
                        kerroin = 2
                    apu = self.z.ravel()
                    for i in range(self.subsets):
                        self.d_z[i] = cuda.gpuarray.to_gpu(apu[self.nMeas[i] * kerroin : self.nMeas[i + 1] * kerroin])
                else:
                    if (self.PET and self.listmode == 0):
                        if self.nLayers > 1:
                            kerroin = 3
                        else:
                            kerroin = 2
                        apu = self.z.ravel()
                        for i in range(self.subsets):
                            self.d_z[i] = cuda.gpuarray.to_gpu(apu[self.nMeas[i] * kerroin : self.nMeas[i + 1] * kerroin])
                    elif self.listmode == 0:
                        self.d_z[0] = cuda.gpuarray.to_gpu(self.z.ravel())
                if (self.attenuation_correction and not self.CTAttenuation):
                    self.d_atten = [None] * self.subsets
                    for i in range(self.subsets):
                        self.d_atten[i] = cuda.gpuarray.to_gpu(self.vaimennus[self.nTotMeas[i].item() : self.nTotMeas[i + 1].item()])
                elif (self.attenuation_correction and self.CTAttenuation):
                    self.d_atten = cuda.gpuarray.to_gpu(self.vaimennus)
                if self.useMaskFP:
                    self.d_maskFP = cuda.gpuarray.to_gpu(self.maskFP)
                if self.useMaskBP:
                    self.d_maskBP = cuda.gpuarray.to_gpu(self.maskBP)
                if self.SPECT:
                    self.d_rayShiftsDetector = cuda.gpuarray.to_gpu(self.rayShiftsDetector)
                    self.d_rayShiftsSource = cuda.gpuarray.to_gpu(self.rayShiftsSource)
                if self.TOF:
                    self.d_TOFCenter = cuda.gpuarray.to_gpu(self.TOFCenter)
                if (self.BPType == 2 or self.BPType == 3 or self.FPType == 2 or self.FPType == 3):
                    self.d_V = cuda.gpuarray.to_gpu(self.V)
                if (self.normalization_correction):
                    self.d_norm = [None] * self.subsets
                    for i in range(self.subsets):
                        self.d_norm[i] = cuda.gpuarray.to_gpu(self.normalization[self.nTotMeas[i].item() : self.nTotMeas[i + 1].item()])
                if (self.additionalCorrection):
                    self.d_corr = [None] * self.subsets
                    for i in range(self.subsets):
                        self.d_corr[i] = cuda.gpuarray.to_gpu(self.corrVector[self.nTotMeas[i].item() : self.nTotMeas[i + 1].item()])
                if (self.listmode != 1 and ((not self.CT and not self.SPECT and not self.PET) and (self.subsets > 1 and (self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7)))):
                    self.d_zindex = [None] * self.subsets
                    self.d_xyindex = [None] * self.subsets
                    for i in range(self.subsets):
                        self.d_xyindex[i] = cuda.gpuarray.to_gpu(self.xy_index[self.nMeas[i] : self.nMeas[i + 1]])
                        self.d_zindex[i] = cuda.gpuarray.to_gpu(self.z_index[self.nMeas[i] : self.nMeas[i + 1]])
                if self.OffsetLimit.size > 0 and ((self.BPType == 4 and self.CT) or self.BPType == 5):
                    self.d_T = [None] * self.subsets
                    for i in range(self.subsets):
                        self.d_T[i] = cuda.gpuarray.to_gpu(self.offsetLimit[self.nMeas[i].item() : self.nMeas[i + 1].item()])
                try:
                    mod = SourceModule(linesFP, options=list(bOptFP), no_extern_c=True)
                except cuda.driver.CompileError as e:
                    print("Compilation error:", e)
                    raise ValueError("Forward projection compilation failed!")
                if self.FPType in [1, 2, 3]:
                    self.knlF = mod.get_function('projectorType123')
                elif self.FPType == 4:
                    self.knlF = mod.get_function('projectorType4Forward')
                elif self.FPType == 5:
                    self.knlF = mod.get_function('projectorType5Forward')
                try:
                    mod = SourceModule(linesBP, options=list(bOptBP), no_extern_c=True)
                except cuda.driver.CompileError as e:
                    print("Compilation error:", e)
                    raise ValueError("Backprojection compilation failed!")
                if self.BPType in [1, 2, 3]:
                    self.knlB = mod.get_function('projectorType123')
                elif self.BPType == 4 and not self.CT:
                    self.knlB = mod.get_function('projectorType4Forward')
                elif self.BPType == 4 and self.CT:
                    self.knlB = mod.get_function('projectorType4Backward')
                elif self.BPType == 5:
                    self.knlB = mod.get_function('projectorType5Backward')
                self.kIndF = (np.float32(self.global_factor), np.float32(self.epps), np.uint32(self.nRowsD), np.uint32(self.det_per_ring), np.float32(self.sigma_x), )
                if self.SPECT:
                    self.kIndF += (self.d_rayShiftsDetector.gpudata, self.d_rayShiftsSource.gpudata, np.float32(self.coneOfResponseStdCoeffA), np.float32(self.coneOfResponseStdCoeffB), np.float32(self.coneOfResponseStdCoeffC), )
                self.kIndF += (np.float32(self.d_dPitch['x'].item()),np.float32(self.d_dPitch['y'].item()), )
                
                if self.use_psf:
                    with open(headerDir + 'auxKernels.cl', encoding="utf8") as f:
                        lines = f.read()
                    lines = hlines + lines
                    bOpt += ('-DCAST=float','-DPSF','-DLOCAL_SIZE=' + str(localSize[0]),'-DLOCAL_SIZE2=' + str(localSize[1]),)
                    mod = SourceModule(lines, options=list(bOpt), no_extern_c=True)
                    self.knlPSF = mod.get_function('Convolution3D_f')
                    self.d_gaussPSF = cuda.gpuarray.to_gpu(self.gaussK.ravel('F'))
                
                if self.FPType in [2,3]:
                    if self.FPType == 2:
                        self.kIndF += (np.float32(self.tube_width_z),)
                    else:
                        self.kIndF += (np.float32(self.tube_radius),)
                    self.kIndF += (np.float32(self.bmin), np.float32(self.bmax), np.float32(self.Vmax),)
                if self.useMaskFP:
                    self.kIndF += (self.d_maskFP.gpudata,)
                if self.FPType in [1, 2, 3]:
                    if self.TOF:
                        self.kIndF += (self.d_TOFCenter.gpudata, )
                    if self.FPType in [2, 3]:
                        self.kIndF += (self.d_V.gpudata, )
                    self.kIndF += (np.uint32(self.nColsD),)
                if self.attenuation_correction and self.CTAttenuation and self.FPType in [1, 2, 3, 4]:
                    self.kIndF += (self.d_atten.gpudata,)
                    
                
                if self.BPType == 4 or self.BPType == 5:
                    self.kIndB = (np.uint32(self.nRowsD), np.uint32(self.nColsD), np.float32(self.d_dPitch['x'].item()),np.float32(self.d_dPitch['y'].item()),)
                if self.BPType == 4 and not self.CT:
                    self.kIndB += (np.float32(self.dL),)
                    self.kIndB += (np.float32(self.global_factor),)
                if self.BPType in [1, 2, 3]:
                    self.kIndB = (np.float32(self.global_factor), np.float32(self.epps), np.uint32(self.nRowsD), np.uint32(self.det_per_ring), np.float32(self.sigma_x), )
                    if self.SPECT:
                        self.kIndB += (self.d_rayShiftsDetector.gpudata, self.d_rayShiftsSource.gpudata, np.float32(self.coneOfResponseStdCoeffA), np.float32(self.coneOfResponseStdCoeffB), np.float32(self.coneOfResponseStdCoeffC), )
                    self.kIndB += (np.float32(self.d_dPitch['x'].item()),np.float32(self.d_dPitch['y'].item()), )
                    if self.BPType in [2, 3]:
                        if self.BPType == 2:
                            self.kIndB  += (np.float32(self.tube_width_z),)
                        else:
                            self.kIndB  += (np.float32(self.tube_radius),)
                        self.kIndB += (np.float32(self.bmin),)
                        self.kIndB += (np.float32(self.bmax),)
                        self.kIndB += (np.float32(self.Vmax),)
                    if self.useMaskFP:
                        self.kIndB += (self.d_maskFP.gpudata,)
                if self.useMaskBP:
                    self.kIndB += (self.d_maskBP.gpudata,)
                if self.BPType in [1, 2, 3]:
                    if self.TOF:
                        self.kIndB += (self.d_TOFCenter.gpudata,)
                    if self.BPType in [2, 3]:
                        self.kIndB += (self.d_V.gpudata,)
                    self.kIndB += (np.uint32(self.nColsD),)
                if self.BPType == 4 and not self.CT and self.TOF:
                    self.kIndB += (self.d_TOFCenter.gpudata,)
                    self.kIndB += (np.float32(self.sigma_x),)
                if self.attenuation_correction and self.CTAttenuation and self.BPType in [1, 2, 3, 4] and not self.CT:
                    self.kIndB += (self.d_atten.gpudata,)
    else:
        if self.projector_type != 6:
            import pyopencl as cl
            
            if self.useAF:
                ctx = af.opencl.get_context(retain=True)
                self.clctx = cl.Context.from_int_ptr(ctx)
                q = af.opencl.get_queue(True)
                self.queue = cl.CommandQueue.from_int_ptr(q)
            else:
                platforms = cl.get_platforms()
                dList = platforms[self.platform].get_devices()
                dList = [dList[self.deviceNum]]
                self.clctx = cl.Context(devices=dList)
                self.queue = cl.CommandQueue(self.clctx)
            
            self.no_norm = 1
            self.mSize = self.nRowsD * self.nColsD * self.nProjections
            # length = [None] * self.subsets
            # for i in range(self.subsets):
            #     if self.subsetType >= 8:
            self.d_Sens = cl.array.empty(self.queue, shape=(1,1), dtype=cl.cltypes.float)
            
            
            self.d_d = [None] * (self.nMultiVolumes + 1)
            self.d_b = [None] * (self.nMultiVolumes + 1)
            self.d_bmax = [None] * (self.nMultiVolumes + 1)
            self.d_Nxyz = [None] * (self.nMultiVolumes + 1)
            self.dSize = [None] * (self.nMultiVolumes + 1)
            self.d_Scale = [None] * (self.nMultiVolumes + 1)
            self.d_Scale4 = [None] * (self.nMultiVolumes + 1)
            for k in range(self.nMultiVolumes + 1):
                self.d_d[k] = cl.cltypes.make_float3(self.dx[k].item(), self.dy[k].item(), self.dz[k].item())
                self.d_b[k] = cl.cltypes.make_float3(self.bx[k].item(), self.by[k].item(), self.bz[k].item())
                self.d_bmax[k] = cl.cltypes.make_float3(self.bx[k].item() + self.Nx[k].item() * self.dx[k].item(), self.by[k].item() + self.Ny[k].item() * self.dy[k].item(), self.bz[k].item() + self.Nz[k].item() * self.dz[k].item())
                self.d_Nxyz[k] = cl.cltypes.make_uint3(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item())
                if (self.FPType == 4 or self.FPType == 5 or self.BPType == 4 or self.BPType == 5):
                    self.d_Scale4[k] = cl.cltypes.make_float3(self.dScaleX4[k].item(), self.dScaleY4[k].item(), self.dScaleZ4[k].item())
                    if self.FPType == 5 or self.BPType == 5:
                        self.dSize[k] = cl.cltypes.make_float2(self.dSizeX[k].item(), self.dSizeY[k].item())
                        self.d_Scale[k] = cl.cltypes.make_float3(self.dScaleX[k].item(), self.dScaleY[k].item(), self.dScaleZ[k].item())
                        if k == 0:
                            self.dSizeBP = cl.cltypes.make_float2(self.dSizeXBP, self.dSizeZBP)
            self.d_dPitch = cl.cltypes.make_float2(self.dPitchX, self.dPitchY)
            self.d_x = [None] * self.subsets
            self.d_z = [None] * self.subsets
            if (self.listmode == 0 and not (self.CT or self.SPECT)) or self.useIndexBasedReconstruction:
                self.d_x[0] = cl.array.to_device(self.queue, self.x.ravel())
            elif (self.CT or self.SPECT) and self.listmode == 0:
                apu = self.x.ravel()
                for i in range(self.subsets):
                    self.d_x[i] = cl.array.to_device(self.queue, apu[self.nMeas[i] * 6 : self.nMeas[i + 1] * 6])
            elif self.listmode > 0 and not self.useIndexBasedReconstruction:
                apu = self.x.ravel()
                for i in range(self.subsets):
                    if self.loadTOF:
                        self.d_x[i] = cl.array.to_device(self.queue, apu[self.nMeas[i] * 6 : self.nMeas[i + 1] * 6])
            if ((self.CT or self.SPECT) and self.listmode == 0):
                if self.pitch:
                    kerroin = 6
                else:
                    kerroin = 2
                apu = self.z.ravel()
                for i in range(self.subsets):
                    self.d_z[i] = cl.array.to_device(self.queue, apu[self.nMeas[i] * kerroin : self.nMeas[i + 1] * kerroin])
            else:
                if (self.PET and self.listmode == 0):
                    if self.nLayers > 1:
                        kerroin = 3
                    else:
                        kerroin = 2
                    apu = self.z.ravel()
                    for i in range(self.subsets):
                        self.d_z[i] = cl.array.to_device(self.queue, apu[self.nMeas[i] * kerroin : self.nMeas[i + 1] * kerroin])
                elif self.listmode == 0 or (self.listmode > 0 and self.useIndexBasedReconstruction):
                    self.d_z[0] = cl.array.to_device(self.queue, self.z.ravel())
                else:
                    for i in range(self.subsets):
                        self.d_z[i] = cl.array.to_device(self.queue, np.zeros(1,dtype=np.float32))
            if (self.attenuation_correction and not self.CTAttenuation):
                self.d_atten = [None] * self.subsets
                for i in range(self.subsets):
                    self.d_atten[i] = cl.array.to_device(self.queue, self.vaimennus[self.nTotMeas[i].item() : self.nTotMeas[i + 1].item()])
            elif (self.attenuation_correction and self.CTAttenuation):
                imformat = cl.ImageFormat(cl.channel_order.A, cl.channel_type.FLOAT)
                self.d_atten = cl.Image(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, imformat, hostbuf=self.vaimennus, shape=(self.Nx[0].item(), self.Ny[0].item(), self.Nz[0].item()))
                # self.d_atten = cl.image_from_array(self.clctx, np.reshape(self.vaimennus, (self.Nx[0].item(), self.Ny[0].item(), self.Nz[0].item()), order='F'))
            if self.SPECT:
                self.d_rayShiftsDetector = cl.array.to_device(self.queue, self.rayShiftsDetector)
                self.d_rayShiftsSource = cl.array.to_device(self.queue, self.rayShiftsSource)
            if self.useMaskFP:
                self.d_maskFP = cl.image_from_array(self.clctx, self.maskFP)
            if self.useMaskBP:
                self.d_maskBP = cl.image_from_array(self.clctx, self.maskBP)
            if self.TOF:
                self.d_TOFCenter = cl.array.to_device(self.queue, self.TOFCenter)
            if (self.BPType == 2 or self.BPType == 3 or self.FPType == 2 or self.FPType == 3):
                self.d_V = cl.array.to_device(self.queue, self.V)
            if (self.normalization_correction):
                self.d_norm = [None] * self.subsets
                for i in range(self.subsets):
                    self.d_norm[i] = cl.array.to_device(self.queue, self.normalization[self.nTotMeas[i].item() : self.nTotMeas[i + 1].item()])
            if (self.additionalCorrection):
                self.d_corr = [None] * self.subsets
                for i in range(self.subsets):
                    self.d_corr[i] = cl.array.to_device(self.queue, self.corrVector[self.nTotMeas[i].item() : self.nTotMeas[i + 1].item()])
            if (self.listmode != 1 and ((not self.CT and not self.SPECT and not self.PET) and (self.subsets > 1 and (self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7)))):
                self.d_zindex = [None] * self.subsets
                self.d_xyindex = [None] * self.subsets
                for i in range(self.subsets):
                    self.d_xyindex[i] = cl.array.to_device(self.queue, self.xy_index[self.nMeas[i] : self.nMeas[i + 1]])
                    self.d_zindex[i] = cl.array.to_device(self.queue, self.z_index[self.nMeas[i] : self.nMeas[i + 1]])
            if (self.listmode > 0 and self.useIndexBasedReconstruction):
                self.d_trIndex = [None] * self.subsets
                self.d_axIndex = [None] * self.subsets
                for i in range(self.subsets):
                    if self.loadTOF:
                        self.d_trIndex[i] = cl.array.to_device(self.queue, self.trIndex[self.nMeas[i] * 2 : self.nMeas[i + 1] * 2])
                        self.d_axIndex[i] = cl.array.to_device(self.queue, self.axIndex[self.nMeas[i] * 2 : self.nMeas[i + 1] * 2])
            if self.OffsetLimit.size > 0 and ((self.BPType == 4 and self.CT) or self.BPType == 5):
                self.d_T = [None] * self.subsets
                for i in range(self.subsets):
                    self.d_T[i] = cl.array.to_device(self.queue, self.offsetLimit[self.nMeas[i].item() : self.nMeas[i + 1].item()])
            # d_Sens = cl.Buffer(clctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=Sens)
            # d_x = cl.Buffer(self.clctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.x)
            # z = cl.Buffer(clctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.z)
            prg = cl.Program(self.clctx, linesFP).build(' '.join(bOptFP))
            if self.FPType in [1, 2, 3]:
                self.knlF = prg.projectorType123
            elif self.FPType == 4:
                self.knlF = prg.projectorType4Forward
            elif self.FPType == 5:
                self.knlF = prg.projectorType5Forward
            prg = cl.Program(self.clctx, linesBP).build(' '.join(bOptBP))
            if self.BPType in [1, 2, 3]:
                self.knlB = prg.projectorType123
            elif self.BPType == 4 and not self.CT:
                self.knlB = prg.projectorType4Forward
            elif self.BPType == 4 and self.CT:
                self.knlB = prg.projectorType4Backward
            elif self.BPType == 5:
                self.knlB = prg.projectorType5Backward
            
            if self.use_psf:
                with open(headerDir + 'auxKernels.cl', encoding="utf8") as f:
                    lines = f.read()
                lines = hlines + lines
                bOpt +=(' -DCAST=float',' -DPSF',' -DLOCAL_SIZE=' + str(localSize[0]), ' -DLOCAL_SIZE2=' + str(localSize[1]),)
                prg = cl.Program(self.clctx, lines).build(' '.join(bOpt))
                self.knlPSF = prg.Convolution3D_f
                self.d_gaussPSF = cl.array.to_device(self.queue, self.gaussK.ravel('F'))
                
            self.kIndF = 0
            if self.FPType == 4 or self.FPType == 5:
                self.knlF.set_arg(self.kIndF, (cl.cltypes.uint)(self.nRowsD))
                self.kIndF += 1
                self.knlF.set_arg(self.kIndF, (cl.cltypes.uint)(self.nColsD))
                self.kIndF += 1
                self.knlF.set_arg(self.kIndF, self.d_dPitch)
                self.kIndF += 1
            if self.FPType == 4:
                self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.dL))
                self.kIndF += 1
                self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.global_factor))
                self.kIndF += 1
            if self.FPType in [1, 2, 3]:
                self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.global_factor))
                self.kIndF += 1
                self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.epps))
                self.kIndF += 1
                self.knlF.set_arg(self.kIndF, (cl.cltypes.uint)(self.nRowsD))
                self.kIndF += 1
                self.knlF.set_arg(self.kIndF, (cl.cltypes.uint)(self.det_per_ring))
                self.kIndF += 1
                self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.sigma_x))
                self.kIndF += 1
                if self.SPECT:
                    self.knlF.set_arg(self.kIndF, self.d_rayShiftsDetector.data)
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, self.d_rayShiftsSource.data)
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.coneOfResponseStdCoeffA))
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.coneOfResponseStdCoeffB))
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.coneOfResponseStdCoeffC))
                    self.kIndF += 1
                self.knlF.set_arg(self.kIndF, self.d_dPitch)
                self.kIndF += 1
                if self.FPType in [2, 3]:
                    if self.FPType == 2:
                        self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.tube_width_z))
                        self.kIndF += 1
                    else:
                        self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.tube_radius))
                        self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.bmin))
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.bmax))
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.Vmax))
                    self.kIndF += 1
            if self.useMaskFP:
                self.knlF.set_arg(self.kIndF, self.d_maskFP)
                self.kIndF += 1
            if self.FPType in [1, 2, 3]:
                if self.TOF:
                    self.knlF.set_arg(self.kIndF, self.d_TOFCenter.data)
                    self.kIndF += 1
                if self.FPType in [2, 3]:
                    self.knlF.set_arg(self.kIndF, self.d_V.data)
                    self.kIndF += 1
                self.knlF.set_arg(self.kIndF, (cl.cltypes.uint)(self.nColsD))
                self.kIndF += 1
            if self.FPType == 4 and not self.CT and self.TOF:
                self.knlF.set_arg(self.kIndF, self.d_TOFCenter.data)
                self.kIndF += 1
                self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.sigma_x))
                self.kIndF += 1
            if self.attenuation_correction and self.CTAttenuation and self.FPType in [1, 2, 3, 4]:
                self.knlF.set_arg(self.kIndF, self.d_atten)
                self.kIndF += 1
                    
                
            
            self.kIndB = 0
            if self.BPType == 4 or self.BPType == 5:
                self.knlB.set_arg(self.kIndB, (cl.cltypes.uint)(self.nRowsD))
                self.kIndB += 1
                self.knlB.set_arg(self.kIndB, (cl.cltypes.uint)(self.nColsD))
                self.kIndB += 1
                self.knlB.set_arg(self.kIndB, self.d_dPitch)
                self.kIndB += 1
            if self.BPType == 4 and not self.CT:
                self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.dL))
                self.kIndB += 1
                self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.global_factor))
                self.kIndB += 1
            if self.BPType in [1, 2, 3]:
                self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.global_factor))
                self.kIndB += 1
                self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.epps))
                self.kIndB += 1
                self.knlB.set_arg(self.kIndB, (cl.cltypes.uint)(self.nRowsD))
                self.kIndB += 1
                self.knlB.set_arg(self.kIndB, (cl.cltypes.uint)(self.det_per_ring))
                self.kIndB += 1
                self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.sigma_x))
                self.kIndB += 1
                if self.SPECT:
                    self.knlB.set_arg(self.kIndB, self.d_rayShiftsDetector.data)
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, self.d_rayShiftsSource.data)
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.coneOfResponseStdCoeffA))
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.coneOfResponseStdCoeffB))
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.coneOfResponseStdCoeffC))
                    self.kIndB += 1
                self.knlB.set_arg(self.kIndB, self.d_dPitch)
                self.kIndB += 1
                if self.BPType in [2, 3]:
                    if self.BPType == 2:
                        self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.tube_width_z))
                        self.kIndB += 1
                    else:
                        self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.tube_radius))
                        self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.bmin))
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.bmax))
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.Vmax))
                    self.kIndB += 1
                if self.useMaskFP:
                    self.knlB.set_arg(self.kIndB, self.d_maskFP)
                    self.kIndB += 1
            if self.useMaskBP:
                self.knlB.set_arg(self.kIndB, self.d_maskBP)
                self.kIndB += 1
            if self.BPType in [1, 2, 3]:
                if self.TOF:
                    self.knlB.set_arg(self.kIndB, self.d_TOFCenter.data)
                    self.kIndB += 1
                if self.BPType in [2, 3]:
                    self.knlB.set_arg(self.kIndB, self.d_V.data)
                    self.kIndB += 1
                self.knlB.set_arg(self.kIndB, (cl.cltypes.uint)(self.nColsD))
                self.kIndB += 1
            if self.BPType == 4 and not self.CT and self.TOF:
                self.knlB.set_arg(self.kIndB, self.d_TOFCenter.data)
                self.kIndB += 1
                self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.sigma_x))
                self.kIndB += 1
            if self.attenuation_correction and self.CTAttenuation and self.BPType in [1, 2, 3, 4] and not self.CT:
                self.knlB.set_arg(self.kIndB, self.d_atten)
                self.kIndB += 1