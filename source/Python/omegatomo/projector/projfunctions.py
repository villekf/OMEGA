# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 13:25:14 2025
"""

import numpy as np

def conv3D(self, f, ii = 0):
    globalSize = (self.Nx[ii].item() + self.erotusBP[ii * 2], self.Ny[ii].item() + self.erotusBP[ii * 2 + 1], self.Nz[ii].item())
    kInd = 0
    if self.useCUDA:
        if self.useTorch:
            import torch
        if self.useCuPy:
            import cupy as cp
        else:
            import pycuda as cuda
    else:
        import pyopencl as cl
    if self.useAF:
        import arrayfire as af
        if isinstance(f, af.array.Array):
            ptr = f.raw_ptr()
            f = cl.MemoryObject.from_int_ptr(ptr)
            self.knlPSF.set_arg(kInd, f)
        else:
            self.knlPSF.set_arg(kInd, f.data)
    else:
        if self.useCUDA:
            if self.useTorch:
                output = torch.zeros(self.N[ii].item(), dtype=torch.float32, device='cuda')
            elif self.useCuPy:
                output = cp.zeros(self.N[ii].item(), dtype=cp.float32)
            else:
                output = cuda.gpuarray.zeros(self.N[ii].item(), dtype=np.float32)
        else:
            output = cl.array.zeros(self.queue, self.N[ii].item(), dtype=cl.cltypes.float)
        if not self.useCUDA:
            self.knlPSF.set_arg(kInd, f.data)
    if self.useCUDA:
        if self.useCuPy:
            fD = cp.asarray(f)
            outputD = cp.asarray(output)
            self.knlPSF((globalSize[0] // 16, globalSize[1] // 16, globalSize[2] // 16), (16,16,1),(fD, outputD, self.d_gaussPSF, cp.int32(self.g_dim_x), cp.int32(self.g_dim_y), cp.int32(self.g_dim_z)))
        else:
            if self.useTorch:
                class Holder(cuda.driver.PointerHolderBase):
                    def __init__(self, t):
                        super(Holder, self).__init__()
                        self.t = t
                        self.gpudata = t.data_ptr()
                    def get_pointer(self):
                        return self.t.data_ptr()
                fD = Holder(f)
                outputD = Holder(output)
                self.knlPSF(fD, outputD, self.d_gaussPSF.gpudata, np.int32(self.g_dim_x), np.int32(self.g_dim_y), np.int32(self.g_dim_z), block=(16,16,1), grid=(globalSize[0], globalSize[1], globalSize[2]))
            else:
                self.knlPSF(f.gpudata, output.gpudata, self.d_gaussPSF.gpudata, np.int32(self.g_dim_x), np.int32(self.g_dim_y), np.int32(self.g_dim_z), block=(16,16,1), grid=(globalSize[0], globalSize[1], globalSize[2]))
        if self.useTorch:
            torch.cuda.synchronize()
    else:
        kInd += 1
        if self.useAF:
            output = af.data.constant(0., self.N[ii].item(), dtype=af.Dtype.f32)
            outPtr = cl.MemoryObject.from_int_ptr(output.raw_ptr())
            self.knlPSF.set_arg(kInd, outPtr)
        else:
            self.knlPSF.set_arg(kInd, output.data)
        kInd += 1
        self.knlPSF.set_arg(kInd, self.d_gaussPSF.data)
        kInd += 1
        self.knlPSF.set_arg(kInd, (cl.cltypes.int)(self.g_dim_x))
        kInd += 1
        self.knlPSF.set_arg(kInd, (cl.cltypes.int)(self.g_dim_y))
        kInd += 1
        self.knlPSF.set_arg(kInd, (cl.cltypes.int)(self.g_dim_z))
        cl.enqueue_nd_range_kernel(self.queue, self.knlPSF, globalSize, (16, 16, 1))
        self.queue.finish()
    if self.useAF:
        af.device.unlock_array(output)
    return output

def forwardProjection(self, f, subset = -1):
    if subset == -1:
        subset = self.subset
    if self.projector_type == 6:
        if not self.useCUDA:
            import arrayfire as af
            u1 = np.sum(self.nProjSubset[:subset])
            y = af.data.constant(0., self.nRowsD, self.nColsD, self.nProjSubset[subset].item())
            for ii in range(self.nMultiVolumes + 1):
                if isinstance(f,list):
                    apuArr = af.data.moddims(f[ii], self.Nx[ii].item(), self.Ny[ii].item(), self.Nz[ii].item())
                else:
                    apuArr = af.data.moddims(f, self.Nx[ii].item(), self.Ny[ii].item(), self.Nz[ii].item())
                for kk in range(self.nProjSubset[subset].item()):
                    kuvaRot = af.image.rotate(apuArr, (180-self.angles[u1].item())*np.pi/180, method=af.INTERP.BILINEAR) # [128, 128, 96]
                    kuvaRot = af.data.reorder(kuvaRot, 2, 1, 0) # [96, 128, 128]
                    kuvaRot = af.signal.convolve2(kuvaRot, self.d_gFilter[:, :, :, u1]) # [96, 128, 128]
                    kuvaRot = kuvaRot[:, :, self.blurPlanes[u1].item():] # [96, 128, 90]
                    kuvaRot = af.data.reorder(kuvaRot, 2, 1, 0) # [90, 128, 96]
                    kuvaRot = af.data.reorder(af.algorithm.sum(kuvaRot, 0), 1, 2, 0) # [1, 128, 96] -->  [128, 96]
                    y[:, :, kk] += kuvaRot
                    u1 += 1
            y = af.flat(y)
            af.eval(y)
        else:
            if not self.useCuPy:
                raise ValueError('SPECT is only supported on CUDA with CuPy!')
            else:
                import torch
                from torchvision.transforms.functional import rotate
                from torchvision.transforms import InterpolationMode
                import torch.nn.functional as F
                u1 = np.sum(self.nProjSubset[:subset])
                y = torch.zeros((self.nProjSubset[subset].item(), self.nColsD, self.nRowsD), dtype=torch.float32).cuda()
                for ii in range(self.nMultiVolumes + 1):
                    if isinstance(f,list):
                        apuArr = torch.reshape(f[ii], (self.Nz[ii].item(), self.Ny[ii].item(), self.Nx[ii].item()))
                    else:
                        apuArr = torch.reshape(f, (self.Nz[ii].item(), self.Ny[ii].item(), self.Nx[ii].item()))
                    for kk in range(self.nProjSubset[subset].item()):
                        kuvaRot = rotate(apuArr, (180-self.angles[u1].item())*np.pi/180, InterpolationMode.BILINEAR)
                        kuvaRot = torch.permute(kuvaRot, (2, 1, 0))
                        kuvaRot = kuvaRot.unsqueeze(0)
                        kuvaRot = F.conv2d(kuvaRot, self.d_gFilter[:, :, :, :, u1], padding=(self.d_gFilter.shape[2] // 2, self.d_gFilter.shape[3] // 2), groups=kuvaRot.shape[1])
                        kuvaRot = kuvaRot.squeeze(0)
                        kuvaRot = kuvaRot[self.blurPlanes[u1].item():, :, :]
                        kuvaRot = torch.sum(kuvaRot, 0)
                        kuvaRot = torch.permute(kuvaRot, (1, 0))
                        y[kk, :, :] += kuvaRot
                        u1 += 1
                y = y.ravel()
    else:
        if self.useCUDA:
            if self.useCuPy:
                import cupy as cp
                if not self.loadTOF:
                    if self.useIndexBasedReconstruction and self.listmode > 0:
                        self.d_trIndex[0] = cp.asarray(self.trIndex[self.nMeas[subset] * 2 : self.nMeas[subset + 1] * 2])
                        self.d_axIndex[0] = cp.asarray(self.axIndex[self.nMeas[subset] * 2 : self.nMeas[subset + 1] * 2])
                    elif self.listmode > 0:
                        apu = self.x.ravel()
                        self.d_x[0] = cp.asarray(apu[self.nMeas[subset] * 6 : self.nMeas[subset + 1] * 6])
                if self.useTorch:
                    import torch
                    if self.subsetType > 7 or self.subsets == 1:
                        y = torch.zeros(self.nRowsD * self.nColsD * self.nProjSubset[subset].item(), dtype=torch.float32, device='cuda')
                    else:
                        y = torch.zeros(self.nMeasSubset[subset].item(), dtype=torch.float32, device='cuda')
                    yD = cp.asarray(y)
                else:
                    if self.subsetType > 7 or self.subsets == 1:
                        y = cp.zeros(self.nRowsD * self.nColsD * self.nProjSubset[subset].item(), dtype=cp.float32)
                    else:
                        y = cp.zeros(self.nMeasSubset[subset].item(), dtype=cp.float32)
                for k in range(self.nMultiVolumes + 1):
                    if isinstance(f,list):
                        if self.use_psf:
                            f[k] = self.computeConvolution(f[k])
                        if self.useTorch:
                            fD = cp.asarray(f[k])
                    else:
                        if self.use_psf:
                            f = self.computeConvolution(f)
                        if self.useTorch:
                            fD = cp.asarray(f)
                    if self.FPType == 5:
                        intIm = cp.zeros((self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item()), dtype=cp.float32, order='F')
                        if self.useTorch:
                            intIm[1:,1:,:] = cp.transpose(fD.reshape((self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()), order='F'), (1, 2, 0))
                        else:
                            if isinstance(f,list):
                                intIm[1:,1:,:] = cp.transpose(f[k].reshape((self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()), order='F'), (1, 2, 0))
                            else:
                                intIm[1:,1:,:] = cp.transpose(f.reshape((self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()), order='F'), (1, 2, 0))
                        intIm = intIm.cumsum(0)
                        intIm = intIm.cumsum(1)
                        intIm = intIm.ravel('F')
                        chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
                        array = cp.cuda.texture.CUDAarray(chl, self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item())
                        array.copy_from(intIm.reshape((self.Nx[k].item(), self.Nz[k].item() + 1, self.Ny[k].item() + 1)))
                        res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                        tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                filterMode=cp.cuda.runtime.cudaFilterModeLinear, normalizedCoords=1)
                        ff2 = cp.cuda.texture.TextureObject(res, tdes)
                        intIm = cp.zeros((self.Nx[k].item() + 1, self.Nz[k].item() + 1, self.Ny[k].item()), dtype=cp.float32, order='F')
                        if self.useTorch:
                            intIm[1:,1:,:] = cp.transpose(fD.reshape((self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()), order='F'), (0, 2, 1))
                        else:
                            if isinstance(f,list):
                                intIm[1:,1:,:] = cp.transpose(f[k].reshape((self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()), order='F'), (0, 2, 1))
                            else:
                                intIm[1:,1:,:] = cp.transpose(f.reshape((self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()), order='F'), (0, 2, 1))
                        intIm = intIm.cumsum(0)
                        intIm = intIm.cumsum(1)
                        intIm = intIm.ravel('F')
                        array2 = cp.cuda.texture.CUDAarray(chl, self.Nx[k].item() + 1, self.Nz[k].item() + 1, self.Ny[k].item())
                        array2.copy_from(intIm.reshape((self.Ny[k].item(), self.Nz[k].item() + 1, self.Nx[k].item() + 1)))
                        res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array2)
                        tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                filterMode=cp.cuda.runtime.cudaFilterModeLinear, normalizedCoords=1)
                        ff = cp.cuda.texture.TextureObject(res, tdes)
                    kIndLoc = self.kIndF
                    if self.FPType == 1 or self.FPType == 2 or self.FPType == 3 or self.FPType == 4:
                        if (self.attenuation_correction and not self.CTAttenuation):
                            kIndLoc += (self.d_atten[subset],)
                    if self.FPType == 5 or self.FPType == 4:
                        kIndLoc += (cp.uint32(self.Nx[k].item()),)
                        kIndLoc += (cp.uint32(self.Ny[k].item()),)
                        kIndLoc += (cp.uint32(self.Nz[k].item()),)
                        kIndLoc += (cp.float32(self.bx[k].item()),)
                        kIndLoc += (cp.float32(self.by[k].item()),)
                        kIndLoc += (cp.float32(self.bz[k].item()),)
                        if self.FPType == 5:
                            kIndLoc += (cp.float32(self.dSizeX[k].item()),)
                            kIndLoc += (cp.float32(self.dSizeY[k].item()),)
                            kIndLoc += (cp.float32(self.dx[k].item()),)
                            kIndLoc += (cp.float32(self.dy[k].item()),)
                            kIndLoc += (cp.float32(self.dz[k].item()),)
                            kIndLoc += (cp.float32(self.dScaleX[k].item()),)
                            kIndLoc += (cp.float32(self.dScaleY[k].item()),)
                            kIndLoc += (cp.float32(self.dScaleZ[k].item()),)
                        else:
                            kIndLoc += (cp.float32(self.bx[k].item() + self.Nx[k].item() * self.dx[k].item()),)
                            kIndLoc += (cp.float32(self.by[k].item() + self.Ny[k].item() * self.dy[k].item()),)
                            kIndLoc += (cp.float32(self.bz[k].item() + self.Nz[k].item() * self.dz[k].item()),)
                            kIndLoc += (cp.float32(self.dScaleX4[k].item()),)
                            kIndLoc += (cp.float32(self.dScaleY4[k].item()),)
                            kIndLoc += (cp.float32(self.dScaleZ4[k].item()),)
                    if self.FPType == 4:
                        if isinstance(f,list):
                            chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
                            array = cp.cuda.texture.CUDAarray(chl, self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item())
                            if self.useTorch:
                                array.copy_from(fD.reshape((self.Nz[k].item(), self.Ny[k].item(), self.Nx[k].item())))
                            else:
                                array.copy_from(f[k].reshape((self.Nz[k].item(), self.Ny[k].item(), self.Nx[k].item())))
                            res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                            tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                    filterMode=cp.cuda.runtime.cudaFilterModeLinear, normalizedCoords=1)
                            ff = cp.cuda.texture.TextureObject(res, tdes)
                            kIndLoc += (ff,)
                        else:
                            chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
                            array = cp.cuda.texture.CUDAarray(chl, self.Nx[0].item(), self.Ny[0].item(), self.Nz[0].item())
                            if self.useTorch:
                                array.copy_from(fD.reshape((self.Nz[0].item(), self.Ny[0].item(), self.Nx[0].item())))
                            else:
                                array.copy_from(f.reshape((self.Nz[0].item(), self.Ny[0].item(), self.Nx[0].item())))
                            res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                            tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                    filterMode=cp.cuda.runtime.cudaFilterModeLinear, normalizedCoords=1)
                            ff = cp.cuda.texture.TextureObject(res, tdes)
                            kIndLoc += (ff,)
                        if self.useTorch:
                            kIndLoc += (yD,)
                        else:
                            kIndLoc += (y,)
                        if (self.listmode == 0 and not self.CT):
                            kIndLoc += (self.d_x[0],)
                        else:
                            kIndLoc += (self.d_x[subset], )
                        if (self.CT or self.PET or self.listmode > 0):
                            kIndLoc += (self.d_z[subset],)
                        else:
                            kIndLoc += (self.d_z[0],)
                        kIndLoc += (cp.int64(self.nProjSubset[subset].item()),)
                        if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                            kIndLoc += (self.d_xyindex[subset],)
                            kIndLoc += (self.d_zindex[subset],)
                        if (self.normalization_correction):
                            kIndLoc += (self.d_norm[subset],)
                        if (self.additionalCorrection):
                            kIndLoc += (self.d_corr[subset],)
                        kIndLoc += (cp.uint8(self.no_norm),)
                        kIndLoc += (cp.uint64(self.nMeasSubset[subset].item()),)
                        kIndLoc += (cp.uint32(subset),)
                        kIndLoc += (cp.int32(k),)
                    elif self.FPType == 5:
                        kIndLoc += (self.d_x[subset], )
                        kIndLoc += (self.d_z[subset],)
                        # self.knlF.set_arg(kIndLoc, d_im)
                        # kIndLoc += 1
                        # self.knlF.set_arg(kIndLoc, d_imInt)
                        kIndLoc += (ff,)
                        kIndLoc += (ff2,)
                        if self.useTorch:
                            kIndLoc += (yD,)
                        else:
                            kIndLoc += (y,)
                        kIndLoc += (cp.int64(self.nProjSubset[subset].item()),)
                        # if self.meanFP:
                    elif self.FPType in [1, 2, 3]:
                        if (self.CT or self.PET or self.SPECT) and self.listmode == 0:
                            kIndLoc += (cp.int64(self.nProjSubset[subset].item()),)
                        if (((self.listmode == 0 and not (self.CT or self.SPECT)) or self.useIndexBasedReconstruction)) or (not self.loadTOF and self.listmode > 0):
                            kIndLoc += (self.d_x[0],)
                        else:
                            kIndLoc += (self.d_x[subset], )
                        if (self.CT or self.PET or self.SPECT or (self.listmode > 0 and not self.useIndexBasedReconstruction)):
                            kIndLoc += (self.d_z[subset],)
                        else:
                            kIndLoc += (self.d_z[0],)
                        if (self.normalization_correction):
                            kIndLoc += (self.d_norm[subset],)
                        if (self.additionalCorrection):
                            kIndLoc += (self.d_corr[subset],)
                        kIndLoc += (self.d_Sens,)
                        kIndLoc += (cp.uint32(self.Nx[k].item()),)
                        kIndLoc += (cp.uint32(self.Ny[k].item()),)
                        kIndLoc += (cp.uint32(self.Nz[k].item()),)
                        kIndLoc += (cp.float32(self.dx[k].item()),)
                        kIndLoc += (cp.float32(self.dy[k].item()),)
                        kIndLoc += (cp.float32(self.dz[k].item()),)
                        kIndLoc += (cp.float32(self.bx[k].item()),)
                        kIndLoc += (cp.float32(self.by[k].item()),)
                        kIndLoc += (cp.float32(self.bz[k].item()),)
                        kIndLoc += (cp.float32(self.bx[k].item() + self.Nx[k].item() * self.dx[k].item()),)
                        kIndLoc += (cp.float32(self.by[k].item() + self.Ny[k].item() * self.dy[k].item()),)
                        kIndLoc += (cp.float32(self.bz[k].item() + self.Nz[k].item() * self.dz[k].item()),)
                        if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                            kIndLoc += (self.d_xyindex[subset],)
                            kIndLoc += (self.d_zindex[subset],)
                        if self.useIndexBasedReconstruction and self.listmode > 0:
                            if not self.loadTOF:
                                kIndLoc += (self.d_trIndex[0],)
                                kIndLoc += (self.d_axIndex[0],)
                            else:
                                kIndLoc += (self.d_trIndex[subset],)
                                kIndLoc += (self.d_axIndex[subset],)
                        # if self.useTorch:
                        #     kIndLoc += (fD,)
                        # else:
                        if isinstance(f,list):
                            if self.useImages:
                                chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
                                array = cp.cuda.texture.CUDAarray(chl, self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item())
                                if self.useTorch:
                                    array.copy_from(fD.reshape((self.Nz[k].item(), self.Ny[k].item(), self.Nx[k].item())))
                                else:
                                    array.copy_from(f[k].reshape((self.Nz[k].item(), self.Ny[k].item(), self.Nx[k].item())))
                                res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                                tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                        filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=0)
                                ff = cp.cuda.texture.TextureObject(res, tdes)
                                kIndLoc += (ff,)
                            else:
                                if self.useTorch:
                                    kIndLoc += (fD,)
                                else:
                                    kIndLoc += (f[k],)
                        else:
                            if self.useImages:
                                chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
                                array = cp.cuda.texture.CUDAarray(chl, self.Nx[0].item(), self.Ny[0].item(), self.Nz[0].item())
                                if self.useTorch:
                                    array.copy_from(fD.reshape((self.Nz[0].item(), self.Ny[0].item(), self.Nx[0].item())))
                                else:
                                    array.copy_from(f.reshape((self.Nz[0].item(), self.Ny[0].item(), self.Nx[0].item())))
                                res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                                tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                        filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=0)
                                ff = cp.cuda.texture.TextureObject(res, tdes)
                                kIndLoc += (ff,)
                            else:
                                if self.useTorch:
                                    kIndLoc += (fD,)
                                else:
                                    kIndLoc += (f,)
                        if self.useTorch:
                            kIndLoc += (yD,)
                        else:
                            kIndLoc += (y,)
                        kIndLoc += (cp.uint8(self.no_norm),)
                        kIndLoc += (cp.uint64(self.nMeasSubset[subset].item()),)
                        kIndLoc += (cp.uint32(subset),)
                        kIndLoc += (cp.int32(k),)
                    self.knlF((self.globalSizeFP[subset][0] // self.localSizeFP[0], self.globalSizeFP[subset][1] // self.localSizeFP[1], self.globalSizeFP[subset][2]), (self.localSizeFP[0], self.localSizeFP[1], 1),kIndLoc)
            else:
                import pycuda as cuda
                if self.useTorch:
                    import torch
                    class Holder(cuda.driver.PointerHolderBase):
                        def __init__(self, t):
                            super(Holder, self).__init__()
                            self.t = t
                            self.gpudata = t.data_ptr()
                        def get_pointer(self):
                            return self.t.data_ptr()
                    if self.subsetType > 7 or self.subsets == 1:
                        y = torch.zeros(self.nRowsD * self.nColsD * self.nProjSubset[subset].item(), dtype=torch.float32, device='cuda')
                    else:
                        y = torch.zeros(self.nMeasSubset[subset].item(), dtype=torch.float32, device='cuda')
                    yD = Holder(y)
                else:
                    if self.subsetType > 7 or self.subsets == 1:
                        y = cuda.gpuarray.zeros(self.nRowsD * self.nColsD * self.nProjSubset[subset].item(), dtype=np.float32)
                    else:
                        y = cuda.gpuarray.zeros(self.nMeasSubset[subset].item(), dtype=np.float32)
                for k in range(self.nMultiVolumes + 1):
                    if isinstance(f,list):
                        if self.use_psf:
                            f[k] = self.computeConvolution(f[k])
                    else:
                        if self.use_psf:
                            f = self.computeConvolution(f)
                    if self.useTorch:
                        fD = Holder(f)
                    kIndLoc = self.kIndF
                    if self.FPType == 1 or self.FPType == 2 or self.FPType == 3 or self.FPType == 4:
                        if (self.attenuation_correction and not self.CTAttenuation):
                            kIndLoc += (self.d_atten[subset].gpudata,)
                    if self.FPType in [1, 2, 3]:
                        if (self.CT or self.PET or self.SPECT) and self.listmode == 0:
                            kIndLoc += (np.int64(self.nProjSubset[subset].item()),)
                        if (self.listmode == 0 and not (self.CT or self.SPECT)):
                            kIndLoc += (self.d_x[0].gpudata,)
                        else:
                            kIndLoc += (self.d_x[subset].gpudata, )
                        if (self.CT or self.PET or self.SPECT or self.listmode > 0):
                            kIndLoc += (self.d_z[subset].gpudata,)
                        else:
                            kIndLoc += (self.d_z[0].gpudata,)
                        if (self.normalization_correction):
                            kIndLoc += (self.d_norm[subset].gpudata,)
                        elif (self.additionalCorrection):
                            kIndLoc += (self.d_corr[subset].gpudata,)
                        kIndLoc += (self.d_Sens.gpudata,)
                        kIndLoc += (np.uint32(self.d_Nxyz[k]['x'].item()),)
                        kIndLoc += (np.uint32(self.d_Nxyz[k]['y'].item()),)
                        kIndLoc += (np.uint32(self.d_Nxyz[k]['z'].item()),)
                        kIndLoc += (np.float32(self.d_d[k]['x'].item()),)
                        kIndLoc += (np.float32(self.d_d[k]['y'].item()),)
                        kIndLoc += (np.float32(self.d_d[k]['z'].item()),)
                        kIndLoc += (np.float32(self.d_b[k]['x'].item()),)
                        kIndLoc += (np.float32(self.d_b[k]['y'].item()),)
                        kIndLoc += (np.float32(self.d_b[k]['z'].item()),)
                        kIndLoc += (np.float32(self.d_bmax[k]['x'].item()),)
                        kIndLoc += (np.float32(self.d_bmax[k]['y'].item()),)
                        kIndLoc += (np.float32(self.d_bmax[k]['z'].item()),)
                        if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                            kIndLoc += (self.d_xyindex[subset].gpudata,)
                            kIndLoc += (self.d_zindex[subset].gpudata,)
                        if self.useTorch:
                            kIndLoc += (fD,)
                        else:
                            if isinstance(f,list):
                                kIndLoc += (f[k].gpudata,)
                            else:
                                kIndLoc += (f.gpudata,)
                        if self.useTorch:
                            kIndLoc += (yD,)
                        else:
                            kIndLoc += (y.gpudata,)
                        kIndLoc += (np.uint8(self.no_norm),)
                        kIndLoc += (np.uint64(self.nMeasSubset[subset].item()),)
                        kIndLoc += (np.uint32(subset),)
                        kIndLoc += (np.int32(k),)
                        self.knlF(*kIndLoc, grid=(self.globalSizeFP[subset][0] // self.localSizeFP[0], self.globalSizeFP[subset][1] // self.localSizeFP[1], self.globalSizeFP[subset][2]), block=(self.localSizeFP[0], self.localSizeFP[1], 1))
            if self.useTorch:
                torch.cuda.synchronize()
            #     if self.useAF:
            #         if isinstance(f,list):
            #             af.device.unlock_array(f[k])
            #         else:
            #             af.device.unlock_array(f)
            # if self.useAF:
            #     af.device.unlock_array(y)
        else:
            import pyopencl as cl
            if not self.loadTOF:
                if self.useIndexBasedReconstruction and self.listmode > 0:
                    self.d_trIndex[0] = cl.array.to_device(self.queue, self.trIndex[self.nMeas[subset] * 2 : self.nMeas[subset + 1] * 2])
                    self.d_axIndex[0] = cl.array.to_device(self.queue, self.axIndex[self.nMeas[subset] * 2 : self.nMeas[subset + 1] * 2])
                elif self.listmode > 0:
                    apu = self.x.ravel()
                    self.d_x[0] = cl.array.to_device(self.queue, apu[self.nMeas[subset] * 6 : self.nMeas[subset + 1] * 6])
            if self.useAF:
                import arrayfire as af
                if self.subsetType > 7 or self.subsets == 1:
                    y = af.data.constant(0., self.nRowsD * self.nColsD * self.nProjSubset[subset].item())
                else:
                    y = af.data.constant(0., self.nMeasSubset[subset].item())
                yPtr = y.raw_ptr()
                yD = cl.MemoryObject.from_int_ptr(yPtr)
            else:
                if self.subsetType > 7 or self.subsets == 1:
                    y = cl.array.zeros(self.queue, self.nRowsD * self.nColsD * self.nProjSubset[subset].item(), dtype=cl.cltypes.float)
                else:
                    y = cl.array.zeros(self.queue, self.nMeasSubset[subset].item(), dtype=cl.cltypes.float)
            imformat = cl.ImageFormat(cl.channel_order.A, cl.channel_type.FLOAT)
            mf = cl.mem_flags
            for k in range(self.nMultiVolumes + 1):
                if self.FPType < 5:
                    d_im = cl.Image(self.clctx, mf.READ_ONLY, imformat, shape=(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()))
                else:
                    d_imInt = cl.Image(self.clctx, mf.READ_ONLY, imformat, shape=(self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item()))
                    d_im = cl.Image(self.clctx, mf.READ_ONLY, imformat, shape=(self.Nx[k].item() + 1, self.Nz[k].item() + 1, self.Ny[k].item()))
                if isinstance(f,list):
                    if self.use_psf:
                        f[k] = self.computeConvolution(f[k])
                    if self.useAF:
                        if self.FPType < 5:
                            fPtr = f[k].raw_ptr()
                            fD = cl.MemoryObject.from_int_ptr(fPtr)
                            cl.enqueue_copy(self.queue, d_im, fD, offset=(0), origin=(0,0,0), region=(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()));
                            af.device.unlock_array(f[k])
                        else:
                            intIm = af.data.constant(0., self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item())
                            if self.meanFP:
                                im = af.reorder(af.moddims(f[k], self.Nx[k].item(), d1=self.Ny[k].item(), d2=self.Nz[k].item()), d0=1, d1=2, d2=0)
                                d_meanFP = af.data.constant(0., self.Nx[k].item() + self.Ny[k].item())
                                d_meanFP[0:self.Nx[k].item()] = af.flat(af.mean(af.mean(im, dim=0), dim=1))
                                im -= af.tile(d_meanFP[0:self.Nx[k].item()], d0=im.shape[0], d1=im.shape[1], d3=1)
                                intIm[1:,1:,:] = af.sat(im)
                                im.eval()
                            else:
                                intIm[1:,1:,:] = af.sat(af.reorder(af.moddims(f[k], self.Nx[k].item(), d1=self.Ny[k].item(), d2=self.Nz[k].item()), d0=1, d1=2, d2=0))
                            intIm.eval()
                            intIm = af.flat(intIm)
                            fPtr = intIm.raw_ptr()
                            fD = cl.MemoryObject.from_int_ptr(fPtr)
                            cl.enqueue_copy(self.queue, d_imInt, fD, offset=(0), origin=(0,0,0), region=(self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item()));
                            af.device.unlock_array(intIm)
                            intIm = af.data.constant(0., self.Nx[k].item() + 1, self.Nz[k].item() + 1, self.Ny[k].item())
                            intIm[1:,1:,:] = af.sat(af.reorder(af.moddims(f[k], self.Nx[k].item(), d1=self.Ny[k].item(), d2=self.Nz[k].item()), d0=0, d1=2, d2=1))
                            intIm.eval()
                            intIm = af.flat(intIm)
                            fPtr = intIm.raw_ptr()
                            fD = cl.MemoryObject.from_int_ptr(fPtr)
                            cl.enqueue_copy(self.queue, d_im, fD, offset=(0), origin=(0,0,0), region=(self.Nx[k].item() + 1, self.Nz[k].item() + 1, self.Ny[k].item()));
                            af.device.unlock_array(intIm)
                    else:
                        cl.enqueue_copy(self.queue, d_im, f[k].data, offset=(0), origin=(0,0,0), region=(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()));
                else:
                    if self.use_psf:
                        f = self.computeConvolution(f)
                    if self.useAF:
                        if self.FPType < 5:
                            fPtr = f.raw_ptr()
                            fD = cl.MemoryObject.from_int_ptr(fPtr)
                            cl.enqueue_copy(self.queue, d_im, fD, offset=(0), origin=(0,0,0), region=(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()));
                            af.device.unlock_array(f)
                        else:
                            intIm = af.data.constant(0., self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item())
                            if self.meanFP:
                                im = af.reorder(af.moddims(f, self.Nx[k].item(), d1=self.Ny[k].item(), d2=self.Nz[k].item()), d0=1, d1=2, d2=0)
                                d_meanFP = af.data.constant(0., self.Nx[k].item() + self.Ny[k].item())
                                d_meanFP[0:self.Nx[k].item()] = af.flat(af.mean(af.mean(im, dim=0), dim=1))
                                im -= af.tile(d_meanFP[0:self.Nx[k].item()], d0=im.shape[0], d1=im.shape[1], d3=1)
                                intIm[1:,1:,:] = af.sat(im)
                                af.eval(im)
                            else:
                                intIm[1:,1:,:] = af.sat(af.reorder(af.moddims(f, self.Nx[k].item(), d1=self.Ny[k].item(), d2=self.Nz[k].item()), d0=1, d1=2, d2=0))
                            af.eval(intIm)
                            intIm = af.flat(intIm)
                            fPtr = intIm.raw_ptr()
                            fD = cl.MemoryObject.from_int_ptr(fPtr)
                            cl.enqueue_copy(self.queue, d_imInt, fD, offset=(0), origin=(0,0,0), region=(self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item()));
                            af.device.unlock_array(intIm)
                            intIm = af.data.constant(0., self.Nx[k].item() + 1, self.Nz[k].item() + 1, self.Ny[k].item())
                            intIm[1:,1:,:] = af.sat(af.reorder(af.moddims(f, self.Nx[k].item(), d1=self.Ny[k].item(), d2=self.Nz[k].item()), d0=0, d1=2, d2=1))
                            af.eval(intIm)
                            intIm = af.flat(intIm)
                            af.sync()
                            fPtr = intIm.raw_ptr()
                            fD = cl.MemoryObject.from_int_ptr(fPtr)
                            cl.enqueue_copy(self.queue, d_im, fD, offset=(0), origin=(0,0,0), region=(self.Nx[k].item() + 1, self.Nz[k].item() + 1, self.Ny[k].item()));
                            af.device.unlock_array(intIm)
                    else:
                        cl.enqueue_copy(self.queue, d_im, f.data, offset=(0), origin=(0,0,0), region=(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()));
                kIndLoc = self.kIndF
                if self.FPType == 1 or self.FPType == 2 or self.FPType == 3 or self.FPType == 4:
                    if (self.attenuation_correction and not self.CTAttenuation):
                        self.knlF.set_arg(kIndLoc, self.d_atten[subset].data)
                        kIndLoc += 1
                    # elif self.attenuation_correction and self.CTAttenuation:
                    #     self.knlF.set_arg(kIndLoc, self.d_atten.data)
                    #     kIndLoc += 1
                if self.FPType == 5 or self.FPType == 4:
                    self.knlF.set_arg(kIndLoc, self.d_Nxyz[k])
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, self.d_b[k])
                    kIndLoc += 1
                    if self.FPType == 5:
                        self.knlF.set_arg(kIndLoc, self.dSize[k])
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, self.d_d[k])
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, self.d_Scale[k])
                        kIndLoc += 1
                    else:
                        self.knlF.set_arg(kIndLoc, self.d_bmax[k])
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, self.d_Scale4[k])
                        kIndLoc += 1
                if self.FPType == 4:
                    self.knlF.set_arg(kIndLoc, d_im)
                    kIndLoc += 1
                    if self.useAF:
                        self.knlF.set_arg(kIndLoc, yD)
                    else:
                        self.knlF.set_arg(kIndLoc, y.data)
                    kIndLoc += 1
                    if (self.listmode == 0 and not self.CT):
                        self.knlF.set_arg(kIndLoc, self.d_x[0].data)
                    else:
                        self.knlF.set_arg(kIndLoc, self.d_x[subset].data)
                    kIndLoc += 1
                    if (self.CT or self.PET or self.listmode > 0):
                        self.knlF.set_arg(kIndLoc, self.d_z[subset].data)
                    else:
                        self.knlF.set_arg(kIndLoc, self.d_z[0].data)
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, (cl.cltypes.long)(self.nProjSubset[subset].item()))
                    kIndLoc += 1
                    if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                        self.knlF.set_arg(kIndLoc, self.d_xyindex[subset].data)
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, self.d_zindex[subset].data)
                        kIndLoc += 1
                    if (self.normalization_correction):
                        self.knlF.set_arg(kIndLoc, self.d_norm[subset].data)
                        kIndLoc += 1
                    elif (self.additionalCorrection):
                        self.knlF.set_arg(kIndLoc, self.d_corr[subset].data)
                        kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, (cl.cltypes.uchar)(self.no_norm))
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, (cl.cltypes.ulong)(self.nMeasSubset[subset].item()))
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, (cl.cltypes.uint)(subset))
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, (cl.cltypes.int)(k))
                elif self.FPType == 5:
                    self.knlF.set_arg(kIndLoc, self.d_x[subset].data)
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, self.d_z[subset].data)
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, d_im)
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, d_imInt)
                    kIndLoc += 1
                    if self.useAF:
                        self.knlF.set_arg(kIndLoc, yD)
                    else:
                        self.knlF.set_arg(kIndLoc, y.data)
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, (cl.cltypes.long)(self.nProjSubset[subset].item()))
                    # if self.meanFP:
                elif self.FPType in [1, 2, 3]:
                    if (self.CT or self.PET or self.SPECT) and self.listmode == 0:
                        self.knlF.set_arg(kIndLoc, (cl.cltypes.long)(self.nProjSubset[subset].item()))
                        kIndLoc += 1
                    if ((self.listmode == 0 or self.useIndexBasedReconstruction) and not (self.CT or self.SPECT)) or (not self.loadTOF and self.listmode > 0):
                        self.knlF.set_arg(kIndLoc, self.d_x[0].data)
                    else:
                        self.knlF.set_arg(kIndLoc, self.d_x[subset].data)
                    kIndLoc += 1
                    if (self.CT or self.PET or self.SPECT or (self.listmode > 0 and not self.useIndexBasedReconstruction)):
                        self.knlF.set_arg(kIndLoc, self.d_z[subset].data)
                    else:
                        self.knlF.set_arg(kIndLoc, self.d_z[0].data)
                    kIndLoc += 1
                    if (self.normalization_correction):
                        self.knlF.set_arg(kIndLoc, self.d_norm[subset].data)
                        kIndLoc += 1
                    elif (self.additionalCorrection):
                        self.knlF.set_arg(kIndLoc, self.d_corr[subset].data)
                        kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, self.d_Sens.data)
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, self.d_Nxyz[k])
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, self.d_d[k])
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, self.d_b[k])
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, self.d_bmax[k])
                    kIndLoc += 1
                    if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                        self.knlF.set_arg(kIndLoc, self.d_xyindex[subset].data)
                        kIndLoc += 1
                        self.knlF.set_arg(kIndLoc, self.d_zindex[subset].data)
                        kIndLoc += 1
                    if self.useIndexBasedReconstruction and self.listmode > 0:
                        if not self.loadTOF:
                            self.knlF.set_arg(kIndLoc, self.d_trIndex[0].data)
                            kIndLoc += 1
                            self.knlF.set_arg(kIndLoc, self.d_axIndex[0].data)
                            kIndLoc += 1
                        else:
                            self.knlF.set_arg(kIndLoc, self.d_trIndex[subset].data)
                            kIndLoc += 1
                            self.knlF.set_arg(kIndLoc, self.d_axIndex[subset].data)
                            kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, d_im)
                    kIndLoc += 1
                    if self.useAF:
                        self.knlF.set_arg(kIndLoc, yD)
                    else:
                        self.knlF.set_arg(kIndLoc, y.data)
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, (cl.cltypes.uchar)(self.no_norm))
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, (cl.cltypes.ulong)(self.nMeasSubset[subset].item()))
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, (cl.cltypes.uint)(subset))
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, (cl.cltypes.int)(k))
                cl.enqueue_nd_range_kernel(self.queue, self.knlF, self.globalSizeFP[subset], self.localSizeFP)
                self.queue.finish()
        if self.useAF:
            af.device.unlock_array(y)
    return y

def backwardProjection(self, y, subset = -1):
    if subset == -1:
        subset = self.subset
    if self.nMultiVolumes > 0:
        f = [None] * (self.nMultiVolumes + 1)
    if self.projector_type == 6:
        if not self.useCUDA:
            import arrayfire as af
            fProj = af.data.moddims(y, self.nRowsD, d1=self.nColsD, d2=self.nProjSubset[subset].item())
            for ii in range(self.nMultiVolumes + 1):
                u1 = np.sum(self.nProjSubset[:subset])
                if self.nMultiVolumes > 0:
                    f[ii] = af.data.constant(0, self.Nx[ii].item() * self.Ny[ii].item() * self.Nz[ii].item(), d1 = self.nProjSubset[subset].item())
                else:
                    f = af.data.constant(0, self.Nx[ii].item() * self.Ny[ii].item() * self.Nz[ii].item(), d1 = self.nProjSubset[subset].item())
                for kk in range(self.nProjSubset[subset].item()):
                    fApu = af.data.constant(0, self.Nx[ii].item(), d1=self.Ny[ii].item(), d2=self.Nz[ii].item()) # [128, 128, 96]
                    kuvaRot = fProj[:,:,kk] # [128, 96]
                    kuvaRot = af.data.reorder(kuvaRot, 1, 0, 2) # [96, 128]
                    kuvaRot = af.signal.convolve2(kuvaRot, self.d_gFilter[:, :, :, u1]) # [96, 128, 128]
                    kuvaRot = kuvaRot[:,:, self.blurPlanes[u1].item():] # [96, 128, 90]
                    kuvaRot = af.data.reorder(kuvaRot, 2, 1, 0) # [90, 128, 96]
                    fApu[self.blurPlanes[u1].item():,:,:] = kuvaRot
                    fApu = af.image.rotate(fApu, (180+self.angles[u1].item())*np.pi/180, method=af.INTERP.BILINEAR)
                    if isinstance(f, list):
                        f[ii][:, kk] = af.flat(fApu)
                    else:
                        f[:, kk] = af.flat(fApu)
                    u1 += 1
                if isinstance(f, list):
                    f[ii] = af.sum(f[ii], 1)
                else:
                    f = af.sum(f, 1)
        else:
            if not self.useCuPy:
                raise ValueError('SPECT is only supported on CUDA with CuPy!')
            else:
                import cupy as cp
                if self.useTorch:
                    import torch
                    from torchvision.transforms.functional import rotate
                    from torchvision.transforms import InterpolationMode
                    import torch.nn.functional as F
                    fProj = torch.reshape(y, (self.nProjSubset[subset].item(), self.nColsD, self.nRowsD))
                for ii in range(self.nMultiVolumes + 1):
                    u1 = np.sum(self.nProjSubset[:subset])
                    if self.nMultiVolumes > 0:
                        f[ii] = torch.zeros((self.nProjSubset[subset].item(), self.Nx[ii].item() * self.Ny[ii].item() * self.Nz[ii].item()), dtype=torch.float32).cuda()
                    else:
                        f = torch.zeros((self.nProjSubset[subset].item(), self.Nx[ii].item() * self.Ny[ii].item() * self.Nz[ii].item()), dtype=torch.float32).cuda()
                    for kk in range(self.nProjSubset[subset].item()):
                        fApu = torch.zeros((self.Nz[ii].item(), self.Ny[ii].item(), self.Nx[ii].item()), dtype=torch.float32).cuda()
                        kuvaRot = fProj[kk,:,:]
                        kuvaRot = torch.permute(kuvaRot, (1, 0)).unsqueeze(0).unsqueeze(0)
                        kuvaRot = F.conv2d(kuvaRot, self.d_gFilter[:, :, :, :, u1], padding=(self.d_gFilter.shape[2] // 2, self.d_gFilter.shape[3] // 2))
                        kuvaRot = kuvaRot.squeeze(0)
                        kuvaRot = kuvaRot[self.blurPlanes[u1].item():, :,:]
                        kuvaRot = torch.permute(kuvaRot, (2, 1, 0))
                        # kuvaRot = cp.transpose(kuvaRot, (2, 1, 0))
                        fApu[:,:,self.blurPlanes[u1].item():] = kuvaRot
                        fApu = rotate(fApu, (180+self.angles[u1].item())*np.pi/180, InterpolationMode.BILINEAR)
                        if isinstance(f, list):
                            f[ii][kk,:] = torch.ravel(fApu)
                        else:
                            f[kk,:] = torch.ravel(fApu)
                        u1 += 1
                    if isinstance(f, list):
                        f[ii] = torch.sum(f[ii], 0),
                    else:
                        f = torch.sum(f, 0)
    else:
        if self.useCUDA:
            if self.useCuPy:
                import cupy as cp
                for k in range(self.nMultiVolumes + 1):
                    if self.useTorch:
                        import torch
                        if self.nMultiVolumes > 0:
                            f[k] = torch.zeros(self.N[k].item(), dtype=torch.float32, device='cuda')
                            fD = cp.asarray(f[k])
                        else:
                            f = torch.zeros(self.N[k].item(), dtype=torch.float32, device='cuda')
                            fD = cp.asarray(f)
                        yD = cp.asarray(y)
                    else:
                        if self.nMultiVolumes > 0:
                            f[k] = cp.zeros(self.N[k].item(), dtype=cp.float32)
                        else:
                            f = cp.zeros(self.N[k].item(), dtype=cp.float32)
                    if self.BPType == 5:
                        # y1 = cp.asarray(np.load('testi.npy'))
                        yy = cp.zeros((self.nRowsD+1,self.nColsD+1,self.nProjSubset[subset].item()), dtype=cp.float32, order='F')
                        if self.useTorch:
                            yy[1:,1:,:] = yD.reshape((self.nRowsD,self.nColsD,self.nProjSubset[subset].item()), order='F')
                        else:
                            yy[1:,1:,:] = y.reshape((self.nRowsD,self.nColsD,self.nProjSubset[subset].item()), order='F')
                        yy = yy.cumsum(0)
                        yy = yy.cumsum(1)
                        yy = yy.ravel(order='F')
                    kIndLoc = self.kIndB
                    if self.BPType in [1, 2, 3]:
                        if (self.attenuation_correction and not self.CTAttenuation):
                            kIndLoc += (self.d_atten[subset],)
                        if (self.CT or self.PET or self.SPECT) and self.listmode == 0:
                            kIndLoc += ((self.nProjSubset[subset].item()),)
                        if ((self.listmode == 0 or self.useIndexBasedReconstruction) and not (self.CT or self.SPECT)) or (not self.loadTOF and self.listmode > 0):
                            kIndLoc += (self.d_x[0],)
                        else:
                            kIndLoc += (self.d_x[subset],)
                        if (self.CT or self.PET or self.SPECT or (self.listmode > 0 and not self.useIndexBasedReconstruction)):
                            kIndLoc += (self.d_z[subset],)
                        else:
                            kIndLoc += (self.d_z[0],)
                        if (self.normalization_correction):
                            kIndLoc += (self.d_norm[subset],)
                        if (self.additionalCorrection):
                            kIndLoc += (self.d_corr[subset],)
                        kIndLoc += (self.d_Sens,)
                        kIndLoc += (cp.uint32(self.Nx[k].item()),)
                        kIndLoc += (cp.uint32(self.Ny[k].item()),)
                        kIndLoc += (cp.uint32(self.Nz[k].item()),)
                        kIndLoc += (cp.float32(self.dx[k].item()),)
                        kIndLoc += (cp.float32(self.dy[k].item()),)
                        kIndLoc += (cp.float32(self.dz[k].item()),)
                        kIndLoc += (cp.float32(self.bx[k].item()),)
                        kIndLoc += (cp.float32(self.by[k].item()),)
                        kIndLoc += (cp.float32(self.bz[k].item()),)
                        kIndLoc += (cp.float32(self.bx[k].item() + self.Nx[k].item() * self.dx[k].item()),)
                        kIndLoc += (cp.float32(self.by[k].item() + self.Ny[k].item() * self.dy[k].item()),)
                        kIndLoc += (cp.float32(self.bz[k].item() + self.Nz[k].item() * self.dz[k].item()),)
                        if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                            kIndLoc += (self.d_xyindex[subset],)
                            kIndLoc += (self.d_zindex[subset],)
                        if self.useIndexBasedReconstruction and self.listmode > 0:
                            if not self.loadTOF:
                                kIndLoc += (self.d_trIndex[0],)
                                kIndLoc += (self.d_axIndex[0],)
                            else:
                                kIndLoc += (self.d_trIndex[subset],)
                                kIndLoc += (self.d_axIndex[subset],)
                        if self.useTorch:
                            kIndLoc += (yD,)
                        else:
                            kIndLoc += (y,)
                        if self.useTorch:
                            kIndLoc += (fD,)
                        else:
                            if self.nMultiVolumes > 0:
                                kIndLoc += (f[k],)
                            else:
                                kIndLoc += (f,)
                        kIndLoc += (cp.uint8(self.no_norm),)
                        kIndLoc += (cp.uint64(self.nMeasSubset[subset].item()),)
                        kIndLoc += (cp.uint32(subset),)
                        kIndLoc += (cp.int32(k),)
                    else:
                        if self.CT:
                            if self.OffsetLimit.size > 0:
                                kIndLoc += (self.d_T[subset],)
                            if self.BPType == 5 or self.BPType == 4:
                                kIndLoc += (cp.uint32(self.Nx[k].item()),)
                                kIndLoc += (cp.uint32(self.Ny[k].item()),)
                                kIndLoc += (cp.uint32(self.Nz[k].item()),)
                                kIndLoc += (cp.float32(self.bx[k].item()),)
                                kIndLoc += (cp.float32(self.by[k].item()),)
                                kIndLoc += (cp.float32(self.bz[k].item()),)
                                kIndLoc += (cp.float32(self.dx[k].item()),)
                                kIndLoc += (cp.float32(self.dy[k].item()),)
                                kIndLoc += (cp.float32(self.dz[k].item()),)
                                if self.BPType == 5:
                                    kIndLoc += (cp.float32(self.dScaleX[k].item()),)
                                    kIndLoc += (cp.float32(self.dScaleY[k].item()),)
                                    kIndLoc += (cp.float32(self.dScaleZ[k].item()),)
                                    kIndLoc += (cp.float32(self.dSizeXBP),)
                                    kIndLoc += (cp.float32(self.dSizeZBP),)
                                else:
                                    kIndLoc += (cp.float32(self.kerroin[k].item()),)
                            if self.BPType == 4:
                                # if self.useTorch:
                                #     kIndLoc += (yD,)
                                #     kIndLoc += (fD,)
                                # else:
                                if self.useImages:
                                    chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
                                    array = cp.cuda.texture.CUDAarray(chl, self.nRowsD, self.nColsD, self.nProjSubset[subset].item())
                                    if self.useTorch:
                                        array.copy_from(yD.reshape((self.nProjSubset[subset].item(), self.nColsD, self.nRowsD)))
                                    else:
                                        array.copy_from(y.reshape((self.nProjSubset[subset].item(), self.nColsD, self.nRowsD)))
                                    res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                                    tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                            filterMode=cp.cuda.runtime.cudaFilterModeLinear, normalizedCoords=1)
                                    yy = cp.cuda.texture.TextureObject(res, tdes)
                                    kIndLoc += (yy,)
                                else:
                                    if self.useTorch:
                                        kIndLoc += (yD,)
                                    else:
                                        kIndLoc += (y,)
                                if self.useTorch:
                                    kIndLoc += (fD,)
                                else:
                                    if isinstance(f, list):
                                        kIndLoc += (f[k],)
                                    else:
                                        kIndLoc += (f,)
                                kIndLoc += (self.d_x[subset],)
                                kIndLoc += (self.d_z[subset],)
                                kIndLoc += (self.d_Sens,)
                            else:
                                kIndLoc += (self.d_x[subset],)
                                kIndLoc += (self.d_z[subset],)
                                # if self.useTorch:
                                #     kIndLoc += (yD,)
                                #     kIndLoc += (fD,)
                                # else:
                                if self.useImages:
                                    chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
                                    array = cp.cuda.texture.CUDAarray(chl, self.nRowsD + 1, self.nColsD + 1, self.nProjSubset[subset].item())
                                    array.copy_from(yy.reshape((self.nProjSubset[subset].item(), self.nColsD + 1, self.nRowsD + 1)))
                                    res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                                    tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                            filterMode=cp.cuda.runtime.cudaFilterModeLinear, normalizedCoords=1)
                                    yy = cp.cuda.texture.TextureObject(res, tdes)
                                    kIndLoc += (yy,)
                                else:
                                    if self.useTorch:
                                        kIndLoc += (yD,)
                                    else:
                                        kIndLoc += (y,)
                                if self.useTorch:
                                    kIndLoc += (fD,)
                                else:
                                    if isinstance(f, list):
                                        kIndLoc += (f[k],)
                                    else:
                                        kIndLoc += (f,)
                                kIndLoc += (self.d_Sens,)
                                # if self.meanBP:
                                #     kIndLoc += (dMeanBP)
                        else:
                            kIndLoc += (cp.uint32(self.Nx[k].item()),)
                            kIndLoc += (cp.uint32(self.Ny[k].item()),)
                            kIndLoc += (cp.uint32(self.Nz[k].item()),)
                            kIndLoc += (cp.float32(self.bx[k].item()),)
                            kIndLoc += (cp.float32(self.by[k].item()),)
                            kIndLoc += (cp.float32(self.bz[k].item()),)
                            kIndLoc += (cp.float32(self.bx[k].item() + self.Nx[k].item() * self.dx[k].item()),)
                            kIndLoc += (cp.float32(self.by[k].item() + self.Ny[k].item() * self.dy[k].item()),)
                            kIndLoc += (cp.float32(self.bz[k].item() + self.Nz[k].item() * self.dz[k].item()),)
                            kIndLoc += (cp.float32(self.dScaleX4[k].item()),)
                            kIndLoc += (cp.float32(self.dScaleY4[k].item()),)
                            kIndLoc += (cp.float32(self.dScaleZ4[k].item()),)
                            # if self.useTorch:
                            #     kIndLoc += (yD,)
                            #     kIndLoc += (fD,)
                            # else:
                            if self.useImages:
                                chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
                                array = cp.cuda.texture.CUDAarray(chl, self.nRowsD, self.nColsD, self.nProjSubset[subset].item())
                                if self.useTorch:
                                    array.copy_from(yD.reshape((self.nProjSubset[subset].item(), self.nColsD, self.nRowsD)))
                                else:
                                    array.copy_from(y.reshape((self.nProjSubset[subset].item(), self.nColsD, self.nRowsD)))
                                res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                                tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                        filterMode=cp.cuda.runtime.cudaFilterModeLinear, normalizedCoords=1)
                                yy = cp.cuda.texture.TextureObject(res, tdes)
                                kIndLoc += (yy,)
                            else:
                                if self.useTorch:
                                    kIndLoc += (yD,)
                                else:
                                    kIndLoc += (y,)
                            if self.useTorch:
                                kIndLoc += (fD,)
                            else:
                                if isinstance(f, list):
                                    kIndLoc += (f[k],)
                                else:
                                    kIndLoc += (f,)
                            if self.listmode == 0 and not self.CT:
                                kIndLoc += (self.d_x[0],)
                            else:
                                kIndLoc += (self.d_x[subset],)
                            if (self.CT or self.PET or self.listmode > 0):
                                kIndLoc += (self.d_z[subset],)
                            else:
                                kIndLoc += (self.d_z[0],)
                            kIndLoc += (cp.int64(self.nProjSubset[subset].item()),)
                            if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                                kIndLoc += (self.d_xyindex[subset],)
                                kIndLoc += (self.d_zindex[subset],)
                            if (self.normalization_correction):
                                kIndLoc += (self.d_norm[subset],)
                            elif (self.additionalCorrection):
                                kIndLoc += (self.d_corr[subset],)
                            kIndLoc += (self.d_Sens,)
                        kIndLoc += (cp.uint8(self.no_norm),)
                        if self.CT:
                            kIndLoc += (cp.int64(self.nProjSubset[subset].item()),)
                        else:
                            kIndLoc += (cp.uint64(self.nMeasSubset[subset].item()),)
                            kIndLoc += (cp.uint32(subset),)
                        kIndLoc += (cp.int32(k),)
                    self.knlB((self.globalSizeBP[subset][k][0] // self.localSizeBP[0], self.globalSizeBP[subset][k][1] // self.localSizeBP[1], self.globalSizeBP[subset][k][2]), (self.localSizeBP[0], self.localSizeBP[1], 1), kIndLoc)
            else:
                import pycuda as cuda
                if self.useTorch:
                    import torch
                    class Holder(cuda.driver.PointerHolderBase):
                        def __init__(self, t):
                            super(Holder, self).__init__()
                            self.t = t
                            self.gpudata = t.data_ptr()
                        def get_pointer(self):
                            return self.t.data_ptr()
                    yD = Holder(y)
                for k in range(self.nMultiVolumes + 1):
                    if self.useTorch:
                        if self.nMultiVolumes > 0:
                            f[k] = torch.zeros(self.N[k].item(), dtype=torch.float32, device='cuda')
                            fD = Holder(f[k])
                        else:
                            f = torch.zeros(self.N[k].item(), dtype=torch.float32, device='cuda')
                            fD = Holder(f)
                    else:
                        if self.nMultiVolumes > 0:
                            f[k] = cuda.gpuarray.zeros(self.N[k].item(), dtype=np.float32)
                        else:
                            f = cuda.gpuarray.zeros(self.N[k].item(), dtype=np.float32)
                    kIndLoc = self.kIndB
                    if self.BPType in [1, 2, 3]:
                        if (self.attenuation_correction and not self.CTAttenuation):
                            kIndLoc += (self.d_atten[subset].gpudata,)
                        if (self.CT or self.PET or self.SPECT) and self.listmode == 0:
                            kIndLoc += (np.int64(self.nProjSubset[subset].item()),)
                        if (self.listmode == 0 and not (self.CT or self.SPECT)):
                            kIndLoc += (self.d_x[0].gpudata,)
                        else:
                            kIndLoc += (self.d_x[subset].gpudata,)
                        if (self.CT or self.PET or self.SPECT or self.listmode > 0):
                            kIndLoc += (self.d_z[subset].gpudata,)
                        else:
                            kIndLoc += (self.d_z[0].gpudata,)
                        if (self.normalization_correction):
                            kIndLoc += (self.d_norm[subset].gpudata,)
                        if (self.additionalCorrection):
                            kIndLoc += (self.d_corr[subset].gpudata,)
                        kIndLoc += (self.d_Sens.gpudata,)
                        kIndLoc += (np.uint32(self.d_Nxyz[k]['x'].item()),)
                        kIndLoc += (np.uint32(self.d_Nxyz[k]['y'].item()),)
                        kIndLoc += (np.uint32(self.d_Nxyz[k]['z'].item()),)
                        kIndLoc += (np.float32(self.d_d[k]['x'].item()),)
                        kIndLoc += (np.float32(self.d_d[k]['y'].item()),)
                        kIndLoc += (np.float32(self.d_d[k]['z'].item()),)
                        kIndLoc += (np.float32(self.d_b[k]['x'].item()),)
                        kIndLoc += (np.float32(self.d_b[k]['y'].item()),)
                        kIndLoc += (np.float32(self.d_b[k]['z'].item()),)
                        kIndLoc += (np.float32(self.d_bmax[k]['x'].item()),)
                        kIndLoc += (np.float32(self.d_bmax[k]['y'].item()),)
                        kIndLoc += (np.float32(self.d_bmax[k]['z'].item()),)
                        if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                            kIndLoc += (self.d_xyindex[subset].gpudata,)
                            kIndLoc += (self.d_zindex[subset].gpudata,)
                        if self.useTorch:
                            kIndLoc += (yD,)
                            kIndLoc += (fD,)
                        else:
                            kIndLoc += (y.gpudata,)
                            if self.nMultiVolumes > 0:
                                kIndLoc += (f[k].gpudata,)
                            else:
                                kIndLoc += (f.gpudata,)
                        kIndLoc += (np.uint8(self.no_norm),)
                        kIndLoc += (np.uint64(self.nMeasSubset[subset].item()),)
                        kIndLoc += (np.uint32(subset),)
                        kIndLoc += (np.int32(k),)
                    else:
                        if self.CT:
                            if self.OffsetLimit.size > 0:
                                kIndLoc += (self.d_T[subset].gpudata,)
                            if self.BPType == 5 or self.BPType == 4:
                                kIndLoc += (np.uint32(self.d_Nxyz[k]['x'].item()),)
                                kIndLoc += (np.uint32(self.d_Nxyz[k]['y'].item()),)
                                kIndLoc += (np.uint32(self.d_Nxyz[k]['z'].item()),)
                                kIndLoc += (np.float32(self.d_b[k]['x'].item()),)
                                kIndLoc += (np.float32(self.d_b[k]['y'].item()),)
                                kIndLoc += (np.float32(self.d_b[k]['z'].item()),)
                                kIndLoc += (np.float32(self.d_d[k]['x'].item()),)
                                kIndLoc += (np.float32(self.d_d[k]['y'].item()),)
                                kIndLoc += (np.float32(self.d_d[k]['z'].item()),)
                                if self.BPType == 5:
                                    kIndLoc += (np.float32(self.dScaleX[k].item()),)
                                    kIndLoc += (np.float32(self.dScaleY[k].item()),)
                                    kIndLoc += (np.float32(self.dScaleZ[k].item()),)
                                    kIndLoc += (np.float32(self.dSizeXBP),)
                                    kIndLoc += (np.float32(self.dSizeZBP),)
                                else:
                                    kIndLoc += (np.float32(self.kerroin[k].item()),)
                            if self.BPType == 4:
                                if self.useTorch:
                                    kIndLoc += (yD,)
                                    kIndLoc += (fD,)
                                else:
                                    kIndLoc += (y.gpudata,)
                                    if isinstance(f, list):
                                        kIndLoc += (f[k].gpudata,)
                                    else:
                                        kIndLoc += (f.gpudata,)
                                kIndLoc += (self.d_x[subset].gpudata,)
                                kIndLoc += (self.d_z[subset].gpudata,)
                                kIndLoc += (self.d_Sens.gpudata,)
                            else:
                                kIndLoc += (self.d_x[subset].gpudata,)
                                kIndLoc += (self.d_z[subset].gpudata,)
                                if self.useTorch:
                                    kIndLoc += (yD,)
                                    kIndLoc += (fD,)
                                else:
                                    kIndLoc += (y.gpudata,)
                                    if isinstance(f, list):
                                        kIndLoc += (f[k].gpudata,)
                                    else:
                                        kIndLoc += (f.gpudata,)
                                kIndLoc += (self.d_Sens.gpudata,)
                                # if self.meanBP:
                                #     kIndLoc += (dMeanBP)
                        else:
                            kIndLoc += (np.uint32(self.d_Nxyz[k]['x'].item()),)
                            kIndLoc += (np.uint32(self.d_Nxyz[k]['y'].item()),)
                            kIndLoc += (np.uint32(self.d_Nxyz[k]['z'].item()),)
                            kIndLoc += (np.float32(self.d_b[k]['x'].item()),)
                            kIndLoc += (np.float32(self.d_b[k]['y'].item()),)
                            kIndLoc += (np.float32(self.d_b[k]['z'].item()),)
                            kIndLoc += (np.float32(self.d_bmax[k]['x'].item()),)
                            kIndLoc += (np.float32(self.d_bmax[k]['y'].item()),)
                            kIndLoc += (np.float32(self.d_bmax[k]['z'].item()),)
                            kIndLoc += (np.float32(self.d_Scale4[k]['x'].item()),)
                            kIndLoc += (np.float32(self.d_Scale4[k]['y'].item()),)
                            kIndLoc += (np.float32(self.d_Scale4[k]['z'].item()),)
                            if self.useTorch:
                                kIndLoc += (yD,)
                                kIndLoc += (fD,)
                            else:
                                kIndLoc += (y.gpudata,)
                                if isinstance(f, list):
                                    kIndLoc += (f[k].gpudata,)
                                else:
                                    kIndLoc += (f.gpudata,)
                            if self.listmode == 0 and not self.CT:
                                kIndLoc += (self.d_x[0].gpudata,)
                            else:
                                kIndLoc += (self.d_x[subset].gpudata,)
                            if (self.CT or self.PET or self.listmode > 0):
                                kIndLoc += (self.d_z[subset].gpudata,)
                            else:
                                kIndLoc += (self.d_z[0].gpudata,)
                            kIndLoc += (np.int64(self.nProjSubset[subset].item()),)
                            if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                                kIndLoc += (self.d_xyindex[subset].gpudata,)
                                kIndLoc += (self.d_zindex[subset].gpudata,)
                            if (self.normalization_correction):
                                kIndLoc += (self.d_norm[subset].gpudata,)
                            elif (self.additionalCorrection):
                                kIndLoc += (self.d_corr[subset].gpudata,)
                            kIndLoc += (self.d_Sens.gpudata,)
                        kIndLoc += (np.uint8(self.no_norm),)
                        if self.CT:
                            kIndLoc += (np.int64(self.nProjSubset[subset].item()),)
                        else:
                            kIndLoc += (np.uint64(self.nMeasSubset[subset].item()),)
                            kIndLoc += (np.uint32(subset),)
                        kIndLoc += (np.int32(k),)
                    self.knlB(*kIndLoc, grid=(self.globalSizeBP[subset][k][0] // self.localSizeBP[0], self.globalSizeBP[subset][k][1] // self.localSizeBP[1], self.globalSizeBP[subset][k][2]), block=(self.localSizeBP[0], self.localSizeBP[1], 1))
                # if self.useAF:
                #     if self.nMultiVolumes > 0:
                #         af.device.unlock_array(f[k])
                #     else:
                #         af.device.unlock_array(f)
                #     af.device.unlock_array(y)
                
                if self.use_psf:
                    if self.nMultiVolumes > 0:
                        f[k] = self.computeConvolution(f[k])
                    else:
                        f = self.computeConvolution(f)
            if self.useTorch:
                torch.cuda.synchronize()
        else:
            import pyopencl as cl
            if self.useAF:
                import arrayfire as af
                cltype = af.Dtype.f32
                if self.use_64bit_atomics:
                    cltype = af.Dtype.u64
                elif self.use_32bit_atomics:
                    cltype = af.Dtype.u32
                yPtr = y.raw_ptr()
                yD = cl.MemoryObject.from_int_ptr(yPtr)
            else:
                cltype = cl.cltypes.float
                if self.use_64bit_atomics:
                    cltype = cl.cltypes.ulong
                elif self.use_32bit_atomics:
                    cltype = cl.cltypes.uint
            if self.CT and self.BPType in [4,5]:
                imformat = cl.ImageFormat(cl.channel_order.A, cl.channel_type.FLOAT)
                if self.BPType < 5:
                    d_im = cl.Image(self.clctx, cl.mem_flags.READ_ONLY, imformat, shape=(self.nRowsD, self.nColsD, self.nProjSubset[subset].item()))
                    if self.useAF:
                        cl.enqueue_copy(self.queue, d_im, yD, offset=(0), origin=(0,0,0), region=(self.nRowsD, self.nColsD, self.nProjSubset[subset].item()));
                    else:
                        cl.enqueue_copy(self.queue, d_im, y.data, offset=(0), origin=(0,0,0), region=(self.nRowsD, self.nColsD, self.nProjSubset[subset].item()));
                else:
                    d_im = cl.Image(self.clctx, cl.mem_flags.READ_ONLY, imformat, shape=(self.nRowsD + 1, self.nColsD + 1, self.nProjSubset[subset].item()))
                    y = af.moddims(y, self.nRowsD, d1=self.nColsD, d2=self.nProjSubset[subset].item())
                    if self.meanBP:
                        d_meanBP = af.mean(af.mean(y, dim=0), dim=1)
                        y -= af.tile(d_meanBP, self.nRowsD, d1=self.nColsD)
                        mPtr = d_meanBP.raw_ptr()
                        dMeanBP = cl.MemoryObject.from_int_ptr(mPtr)
                    y = af.sat(y)
                    y = af.join(0, af.data.constant(0, 1, d1=y.shape[1], d2=y.shape[2]), y)
                    y = af.flat(af.join(1, af.data.constant(0, y.shape[0], 1, y.shape[2]), y))
                    yPtr = y.raw_ptr()
                    yD = cl.MemoryObject.from_int_ptr(yPtr)
                    cl.enqueue_copy(self.queue, d_im, yD, offset=(0), origin=(0,0,0), region=(self.nRowsD + 1, self.nColsD + 1, self.nProjSubset[subset].item()));
            
            for k in range(self.nMultiVolumes + 1):
                if self.useAF:
                    if self.nMultiVolumes > 0:
                        f[k] = af.data.constant(0, self.N[k].item(), dtype=cltype)
                        fPtr = f[k].raw_ptr()
                    else:
                        f = af.data.constant(0, self.N[k].item(), dtype=cltype)
                        fPtr = f.raw_ptr()
                    fD = cl.MemoryObject.from_int_ptr(fPtr)
                else:
                    if self.nMultiVolumes > 0:
                        f[k] = cl.array.zeros(self.queue, self.N[k].item(), dtype=cltype)
                    else:
                        f = cl.array.zeros(self.queue, self.N[k].item(), dtype=cltype)
                kIndLoc = self.kIndB
                if self.BPType in [1, 2, 3]:
                    if (self.attenuation_correction and not self.CTAttenuation):
                        self.knlB.set_arg(kIndLoc, self.d_atten[subset].data)
                        kIndLoc += 1
                    if (self.CT or self.PET or self.SPECT) and self.listmode == 0:
                        self.knlB.set_arg(kIndLoc, (cl.cltypes.long)(self.nProjSubset[subset].item()))
                        kIndLoc += 1
                    if ((self.listmode == 0 or self.useIndexBasedReconstruction) and not (self.CT or self.SPECT)) or (not self.loadTOF and self.listmode > 0):
                        self.knlB.set_arg(kIndLoc, self.d_x[0].data)
                    else:
                        self.knlB.set_arg(kIndLoc, self.d_x[subset].data)
                    kIndLoc += 1
                    if (self.CT or self.PET or self.SPECT or (self.listmode > 0 and not self.useIndexBasedReconstruction)):
                        self.knlB.set_arg(kIndLoc, self.d_z[subset].data)
                    else:
                        self.knlB.set_arg(kIndLoc, self.d_z[0].data)
                    kIndLoc += 1
                    if (self.normalization_correction):
                        self.knlB.set_arg(kIndLoc, self.d_norm[subset].data)
                        kIndLoc += 1
                    if (self.additionalCorrection):
                        self.knlB.set_arg(kIndLoc, self.d_corr[subset].data)
                        kIndLoc += 1
                    self.knlB.set_arg(kIndLoc, self.d_Sens.data)
                    kIndLoc += 1
                    self.knlB.set_arg(kIndLoc, self.d_Nxyz[k])
                    kIndLoc += 1
                    self.knlB.set_arg(kIndLoc, self.d_d[k])
                    kIndLoc += 1
                    self.knlB.set_arg(kIndLoc, self.d_b[k])
                    kIndLoc += 1
                    self.knlB.set_arg(kIndLoc, self.d_bmax[k])
                    kIndLoc += 1
                    if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                        self.knlB.set_arg(kIndLoc, self.d_xyindex[subset].data)
                        kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, self.d_zindex[subset].data)
                        kIndLoc += 1
                    if self.useIndexBasedReconstruction and self.listmode > 0:
                        if not self.loadTOF:
                            self.knlB.set_arg(kIndLoc, self.d_trIndex[0].data)
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, self.d_axIndex[0].data)
                            kIndLoc += 1
                        else:
                            self.knlB.set_arg(kIndLoc, self.d_trIndex[subset].data)
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, self.d_axIndex[subset].data)
                            kIndLoc += 1
                    if self.useAF:
                        self.knlB.set_arg(kIndLoc, yD)
                        kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, fD)
                    else:
                        self.knlB.set_arg(kIndLoc, y.data)
                        kIndLoc += 1
                        if self.nMultiVolumes > 0:
                            self.knlB.set_arg(kIndLoc, f[k].data)
                        else:
                            self.knlB.set_arg(kIndLoc, f.data)
                    kIndLoc += 1
                    self.knlB.set_arg(kIndLoc, (cl.cltypes.uchar)(self.no_norm))
                    kIndLoc += 1
                    self.knlB.set_arg(kIndLoc, (cl.cltypes.ulong)(self.nMeasSubset[subset].item()))
                    kIndLoc += 1
                    self.knlB.set_arg(kIndLoc, (cl.cltypes.uint)(subset))
                    kIndLoc += 1
                    self.knlB.set_arg(kIndLoc, (cl.cltypes.int)(k))
                else:
                    if self.CT:
                        if self.OffsetLimit.size > 0:
                            self.knlB.set_arg(kIndLoc, self.d_T[subset].data)
                            kIndLoc += 1
                        if self.BPType == 5 or self.BPType == 4:
                            self.knlB.set_arg(kIndLoc, self.d_Nxyz[k])
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, self.d_b[k])
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, self.d_d[k])
                            kIndLoc += 1
                            if self.BPType == 5:
                                self.knlB.set_arg(kIndLoc, self.d_Scale[k])
                                kIndLoc += 1
                                self.knlB.set_arg(kIndLoc, self.dSizeBP)
                                kIndLoc += 1
                            else:
                                self.knlB.set_arg(kIndLoc, (cl.cltypes.float)(self.kerroin[k].item()))
                                kIndLoc += 1
                        if self.BPType == 4:
                            self.knlB.set_arg(kIndLoc, d_im)
                            kIndLoc += 1
                            if self.useAF:
                                self.knlB.set_arg(kIndLoc, fD)
                            else:
                                if isinstance(f, list):
                                    self.knlB.set_arg(kIndLoc, f[k].data)
                                else:
                                    self.knlB.set_arg(kIndLoc, f.data)
                            kIndLoc += 1
                            if not self.loadTOF and self.listmode > 0:
                                self.knlB.set_arg(kIndLoc, self.d_x[0].data)
                            else:
                                self.knlB.set_arg(kIndLoc, self.d_x[subset].data)
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, self.d_z[subset].data)
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, self.d_Sens.data)
                            kIndLoc += 1
                        else:
                            if not self.loadTOF and self.listmode > 0:
                                self.knlB.set_arg(kIndLoc, self.d_x[0].data)
                            else:
                                self.knlB.set_arg(kIndLoc, self.d_x[subset].data)
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, self.d_z[subset].data)
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, d_im)
                            kIndLoc += 1
                            if self.useAF:
                                self.knlB.set_arg(kIndLoc, fD)
                            else:
                                if isinstance(f, list):
                                    self.knlB.set_arg(kIndLoc, f[k].data)
                                else:
                                    self.knlB.set_arg(kIndLoc, f.data)
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, self.d_Sens.data)
                            kIndLoc += 1
                            if self.meanBP:
                                self.knlB.set_arg(kIndLoc, dMeanBP)
                                kIndLoc += 1
                    else:
                        self.knlB.set_arg(kIndLoc, self.d_Nxyz[k])
                        kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, self.d_b[k])
                        kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, self.d_bmax[k])
                        kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, self.d_Scale4[k])
                        kIndLoc += 1
                        if self.useAF:
                            self.knlB.set_arg(kIndLoc, yD)
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, fD)
                        else:
                            self.knlB.set_arg(kIndLoc, y.data)
                            kIndLoc += 1
                            if isinstance(f, list):
                                self.knlB.set_arg(kIndLoc, f[k].data)
                            else:
                                self.knlB.set_arg(kIndLoc, f.data)
                        kIndLoc += 1
                        if ((self.listmode == 0 or self.useIndexBasedReconstruction) and not self.CT) or (not self.loadTOF and self.listmode > 0):
                            self.knlB.set_arg(kIndLoc, self.d_x[0].data)
                        else:
                            self.knlB.set_arg(kIndLoc, self.d_x[subset].data)
                        kIndLoc += 1
                        if (self.CT or self.PET or (self.listmode > 0 and not self.useIndexBasedReconstruction)):
                            self.knlB.set_arg(kIndLoc, self.d_z[subset].data)
                        else:
                            self.knlB.set_arg(kIndLoc, self.d_z[0].data)
                        kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, (cl.cltypes.ulong)(self.nProjSubset[subset].item()))
                        kIndLoc += 1
                        if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                            self.knlB.set_arg(kIndLoc, self.d_xyindex[subset].data)
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, self.d_zindex[subset].data)
                            kIndLoc += 1
                        if (self.normalization_correction):
                            self.knlB.set_arg(kIndLoc, self.d_norm[subset].data)
                            kIndLoc += 1
                        elif (self.additionalCorrection):
                            self.knlB.set_arg(kIndLoc, self.d_corr[subset].data)
                            kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, self.d_Sens.data)
                        kIndLoc += 1
                    self.knlB.set_arg(kIndLoc, (cl.cltypes.uchar)(self.no_norm))
                    kIndLoc += 1
                    if self.CT:
                        self.knlB.set_arg(kIndLoc, (cl.cltypes.ulong)(self.nProjSubset[subset].item()))
                        kIndLoc += 1
                    else:
                        self.knlB.set_arg(kIndLoc, (cl.cltypes.ulong)(self.nMeasSubset[subset].item()))
                        kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, (cl.cltypes.uint)(subset))
                        kIndLoc += 1
                    self.knlB.set_arg(kIndLoc, (cl.cltypes.int)(k))
                            
                cl.enqueue_nd_range_kernel(self.queue, self.knlB, self.globalSizeBP[subset][k], self.localSizeBP)
                self.queue.finish()
                if self.useAF:
                    if self.nMultiVolumes > 0:
                        af.device.unlock_array(f[k])
                    else:
                        af.device.unlock_array(f)
                    af.device.unlock_array(y)
                    if self.use_64bit_atomics:
                        if self.nMultiVolumes > 0:
                            f[k] = f[k].as_type(af.Dtype.f32) / self.TH
                        else:
                            f = f.as_type(af.Dtype.f32) / self.TH
                    elif self.use_32bit_atomics:
                        if self.nMultiVolumes > 0:
                            f[k] = f[k].as_type(af.Dtype.f32) / self.TH32
                        else:
                            f = f.as_type(af.Dtype.f32) / self.TH32
                else:
                    if self.use_64bit_atomics:
                        if self.nMultiVolumes > 0:
                            f[k] = f[k].astype(cl.cltypes.float) / self.TH
                        else:
                            f = f.astype(cl.cltypes.float) / self.TH
                    elif self.use_32bit_atomics:
                        if self.nMultiVolumes > 0:
                            f[k] = f[k].astype(cl.cltypes.float) / self.TH32
                        else:
                            f = f.astype(cl.cltypes.float) / self.TH32
                if self.use_psf:
                    if self.nMultiVolumes > 0:
                        f[k] = self.computeConvolution(f[k])
                    else:
                        f = self.computeConvolution(f)
    return f