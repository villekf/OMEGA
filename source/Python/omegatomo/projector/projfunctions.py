# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 13:25:14 2025
"""

import numpy as np

def _mask_fp_resource(self, subset):
    if self.SPECT and self.maskFPZ == self.nHeads:
        return self.d_maskFP
    if self.maskFPZ > 1:
        return self.d_maskFP[subset]
    return self.d_maskFP


def conv3D(self, f, ii = 0):
    if getattr(self, "useMetal", False):
        raise NotImplementedError("The separate PSF convolution kernel is not yet wired to the Metal/MPS bridge.")
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
            if self.useTorch:
                fD = cp.asarray(f)
                outputD = cp.asarray(output)
                self.knlPSF((globalSize[0] // 16, globalSize[1] // 16, globalSize[2] // 1), (16,16,1),(fD, outputD, self.d_gaussPSF, cp.int32(self.g_dim_x), cp.int32(self.g_dim_y), cp.int32(self.g_dim_z)))
            else:
                self.knlPSF((globalSize[0] // 16, globalSize[1] // 16, globalSize[2] // 1), (16,16,1),(f, output, self.d_gaussPSF, cp.int32(self.g_dim_x), cp.int32(self.g_dim_y), cp.int32(self.g_dim_z)))
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

def forwardProjection(self, f, subset: int = -1, timestep: int = -1):
    if subset == -1:
        subset = self.subset
    if timestep == -1:
        timestep = self.timestep
    timestep = int(timestep)
    subset = int(subset)
    if self.useMetal:
        from omegatomo.projector.mps_backend import forward_projection_mps
        return forward_projection_mps(self, f, subset, timestep)
    if self.useCUDA and not self.useCuPy:
        raise ValueError('PyCUDA is no longer supported. Please use CuPy.')
    volumes = 0
    if self.projector_type in (6, 66):
        volume_count = int(self.nMultiVolumes) + 1
        inputs = list(f) if isinstance(f, (list, tuple)) else [f] * volume_count
        if len(inputs) != volume_count:
            raise ValueError(f'Expected {volume_count} volume inputs, got {len(inputs)}')
        rows = int(getattr(self, 'measurement_nRowsD', self.nRowsD))
        cols = int(getattr(self, 'measurement_nColsD', self.nColsD))
        if not self.useCUDA:
            import arrayfire as af
            ops = _type6_arrayfire_ops(self)
            for volume, image in enumerate(inputs):
                count = image.elements()
                if count != int(np.asarray(self.N).reshape(-1)[volume]):
                    raise ValueError(f'Volume {volume} has {count} elements; expected {self.N[volume]}')
            output = af.data.constant(0., rows * cols * int(self.nProjSubset[timestep, subset]))
            for volume, image in enumerate(inputs):
                partial = af.data.constant(0., output.elements())
                type6_forward(self, image, partial, volume, subset, timestep, ops=ops)
                ops.add_to(output, partial)
            af.eval(output)
            return output
        else:
            import torch
            ops = _type6_torch_ops(self)
            for volume, image in enumerate(inputs):
                image = image.contiguous()
                count = image.numel()
                if count != int(np.asarray(self.N).reshape(-1)[volume]):
                    raise ValueError(f'Volume {volume} has {count} elements; expected {self.N[volume]}')
                inputs[volume] = image
            output = torch.zeros(
                rows * cols * int(self.nProjSubset[timestep, subset]),
                dtype=torch.float32,
                device=inputs[0].device,
            )
            for volume, image in enumerate(inputs):
                partial = torch.zeros_like(output)
                type6_forward(self, image, partial, volume, subset, timestep, ops=ops)
                output += partial
            return output
    else:
        if self.nMultiVolumes > 0 and not(isinstance(f,list)):
            volumes = self.nMultiVolumes
            self.nMultiVolumes = 0
        if self.useCUDA:
            if self.useCuPy:
                import cupy as cp
                if not self.loadTOF:
                    if self.useIndexBasedReconstruction and self.listmode > 0:
                        self.d_trIndex[0] = cp.asarray(self.trIndex[self.nMeas[timestep * self.subsets + subset] * 2 : self.nMeas[timestep * self.subsets + subset + 1] * 2])
                        self.d_axIndex[0] = cp.asarray(self.axIndex[self.nMeas[timestep * self.subsets + subset] * 2 : self.nMeas[timestep * self.subsets + subset + 1] * 2])
                    elif self.listmode > 0:
                        apu = self.x.ravel()
                        self.d_x[timestep][0] = cp.asarray(apu[self.nMeas[timestep * self.subsets + subset] * 6 : self.nMeas[timestep * self.subsets + subset + 1] * 6])
                if self.useTorch:
                    import torch
                    if self.subsetType > 7 or self.subsets == 1:
                        y = torch.zeros(self.nRowsD * self.nColsD * self.nProjSubset[timestep, subset].item(), dtype=torch.float32, device='cuda')
                    else:
                        y = torch.zeros(self.nMeasSubset[timestep, subset].item(), dtype=torch.float32, device='cuda')
                    yD = cp.asarray(y)
                else:
                    if self.subsetType > 7 or self.subsets == 1:
                        y = cp.zeros(self.nRowsD * self.nColsD * self.nProjSubset[timestep, subset].item(), dtype=cp.float32)
                    else:
                        y = cp.zeros(self.nMeasSubset[timestep, subset].item(), dtype=cp.float32)
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
                            kIndLoc += (self.d_x[timestep][0],)
                        else:
                            kIndLoc += (self.d_x[timestep][subset], )
                        if (self.CT or self.PET or self.listmode > 0):
                            kIndLoc += (self.d_z[timestep][subset],)
                        else:
                            kIndLoc += (self.d_z[timestep][0],)
                        if self.useMaskFP:
                            kIndLoc += (_mask_fp_resource(self, subset),)
                        kIndLoc += (cp.int64(self.nProjSubset[timestep, subset].item()),)
                        if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                            kIndLoc += (self.d_xyindex[subset],)
                            kIndLoc += (self.d_zindex[subset],)
                        if (self.normalization_correction):
                            kIndLoc += (self.d_norm[subset],)
                        if (self.additionalCorrection):
                            kIndLoc += (self.d_corr[subset],)
                        kIndLoc += (cp.uint8(self.no_norm),)
                        kIndLoc += (cp.uint64(self.nMeasSubset[timestep, subset].item()),)
                        kIndLoc += (cp.uint32(subset),)
                        kIndLoc += (cp.int32(k),)
                    elif self.FPType == 5:
                        kIndLoc += (self.d_x[timestep][subset], )
                        kIndLoc += (self.d_z[timestep][subset],)
                        # self.knlF.set_arg(kIndLoc, d_im)
                        # kIndLoc += 1
                        # self.knlF.set_arg(kIndLoc, d_imInt)
                        kIndLoc += (ff,)
                        kIndLoc += (ff2,)
                        if self.useTorch:
                            kIndLoc += (yD,)
                        else:
                            kIndLoc += (y,)
                        if self.useMaskFP:
                            kIndLoc += (_mask_fp_resource(self, subset),)
                        kIndLoc += (cp.int64(self.nProjSubset[timestep, subset].item()),)
                        # if self.meanFP:
                    elif self.FPType in [1, 2, 3]:
                        if self.useMaskFP:
                            kIndLoc += (_mask_fp_resource(self, subset),)
                        if (self.CT or self.PET or self.SPECT) and self.listmode == 0:
                            kIndLoc += (cp.int64(self.nProjSubset[timestep, subset].item()),)
                        if (((self.listmode == 0 and not (self.CT or self.SPECT)) or self.useIndexBasedReconstruction)) or (not self.loadTOF and self.listmode > 0):
                            kIndLoc += (self.d_x[timestep][0],)
                        else:
                            kIndLoc += (self.d_x[timestep][subset], )
                        if (self.CT or self.PET or self.SPECT or (self.listmode > 0 and not self.useIndexBasedReconstruction)):
                            kIndLoc += (self.d_z[timestep][subset],)
                        else:
                            kIndLoc += (self.d_z[timestep][0],)
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
                                    apuArray = fD.reshape((self.Nx[0].item(), self.Ny[0].item(), self.Nz[0].item()), order='F')
                                    apuArray = np.transpose(apuArray, (2, 1, 0))
                                    array.copy_from(apuArray)
                                    # array.copy_from(fD.reshape((self.Nz[0].item(), self.Ny[0].item(), self.Nx[0].item())))
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
                        if self.SPECT:
                            kIndLoc += (self.d_detectorVector[timestep][subset],)
                        kIndLoc += (cp.uint8(self.no_norm),)
                        kIndLoc += (cp.uint64(self.nMeasSubset[timestep, subset].item()),)
                        kIndLoc += (cp.uint32(subset),)
                        kIndLoc += (cp.int32(k),)
                    self.knlF((self.globalSizeFP[timestep][subset][0] // self.localSizeFP[0], self.globalSizeFP[timestep][subset][1] // self.localSizeFP[1], self.globalSizeFP[timestep][subset][2]), (self.localSizeFP[0], self.localSizeFP[1], 1),kIndLoc)
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
            from pyopencl.version import VERSION
            if not self.loadTOF:
                if self.useIndexBasedReconstruction and self.listmode > 0:
                    self.d_trIndex[0] = cl.array.to_device(self.queue, self.trIndex[self.nMeas[timestep * self.subsets + subset] * 2 : self.nMeas[timestep * self.subsets + subset + 1] * 2])
                    self.d_axIndex[0] = cl.array.to_device(self.queue, self.axIndex[self.nMeas[timestep * self.subsets + subset] * 2 : self.nMeas[timestep * self.subsets + subset + 1] * 2])
                elif self.listmode > 0:
                    apu = self.x.ravel()
                    self.d_x[timestep][0] = cl.array.to_device(self.queue, apu[self.nMeas[timestep * self.subsets + subset] * 6 : self.nMeas[timestep * self.subsets + subset + 1] * 6])
            if self.useAF:
                import arrayfire as af
                if self.subsetType > 7 or self.subsets == 1:
                    y = af.data.constant(0., self.nRowsD * self.nColsD * self.nProjSubset[timestep, subset].item())
                else:
                    y = af.data.constant(0., self.nMeasSubset[timestep, subset].item())
                yPtr = y.raw_ptr()
                yD = cl.MemoryObject.from_int_ptr(yPtr)
            else:
                if self.subsetType > 7 or self.subsets == 1:
                    y = cl.array.zeros(self.queue, self.nRowsD * self.nColsD * self.nProjSubset[timestep, subset].item(), dtype=cl.cltypes.float)
                else:
                    y = cl.array.zeros(self.queue, self.nMeasSubset[timestep, subset].item(), dtype=cl.cltypes.float)
            imformat = cl.ImageFormat(cl.channel_order.A, cl.channel_type.FLOAT)
            mf = cl.mem_flags
            for k in range(self.nMultiVolumes + 1):
                if self.useImages:
                    if self.FPType < 5:
                        if VERSION[0] > 2024 or (VERSION[0] == 2024 and VERSION[1] > 2):
                            d_im = cl.create_image(self.clctx, mf.READ_ONLY, imformat, shape=(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()))
                        else:
                            d_im = cl.Image(self.clctx, mf.READ_ONLY, imformat, shape=(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()))
                    else:
                        if VERSION[0] > 2024 or (VERSION[0] == 2024 and VERSION[1] > 2):
                            d_imInt = cl.create_image(self.clctx, mf.READ_ONLY, imformat, shape=(self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item()))
                            d_im = cl.create_image(self.clctx, mf.READ_ONLY, imformat, shape=(self.Nx[k].item() + 1, self.Nz[k].item() + 1, self.Ny[k].item()))
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
                                    af.eval(im)
                                else:
                                    intIm[1:,1:,:] = af.sat(af.reorder(af.moddims(f[k], self.Nx[k].item(), d1=self.Ny[k].item(), d2=self.Nz[k].item()), d0=1, d1=2, d2=0))
                                af.eval(intIm)
                                intIm = af.flat(intIm)
                                fPtr = intIm.raw_ptr()
                                fD = cl.MemoryObject.from_int_ptr(fPtr)
                                cl.enqueue_copy(self.queue, d_imInt, fD, offset=(0), origin=(0,0,0), region=(self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item()));
                                af.device.unlock_array(intIm)
                                intIm = af.data.constant(0., self.Nx[k].item() + 1, self.Nz[k].item() + 1, self.Ny[k].item())
                                intIm[1:,1:,:] = af.sat(af.reorder(af.moddims(f[k], self.Nx[k].item(), d1=self.Ny[k].item(), d2=self.Nz[k].item()), d0=0, d1=2, d2=1))
                                af.eval(intIm)
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
                else:
                    if self.useAF:
                        if isinstance(f,list):
                            fPtr = f[k].raw_ptr()
                        else:
                            fPtr = f.raw_ptr()
                        d_im = cl.MemoryObject.from_int_ptr(fPtr)
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
                    if not self.useImages:
                        raise ValueError('Projector type 4 forward projection only works with images!')
                    self.knlF.set_arg(kIndLoc, d_im)
                    kIndLoc += 1
                    if self.useAF:
                        self.knlF.set_arg(kIndLoc, yD)
                    else:
                        self.knlF.set_arg(kIndLoc, y.data)
                    kIndLoc += 1
                    if (self.listmode == 0 and not self.CT):
                        self.knlF.set_arg(kIndLoc, self.d_x[timestep][0].data)
                    else:
                        self.knlF.set_arg(kIndLoc, self.d_x[timestep][subset].data)
                    kIndLoc += 1
                    if (self.CT or self.PET or self.listmode > 0):
                        self.knlF.set_arg(kIndLoc, self.d_z[timestep][subset].data)
                    else:
                        self.knlF.set_arg(kIndLoc, self.d_z[timestep][0].data)
                    kIndLoc += 1
                    if self.useMaskFP:
                        self.knlF.set_arg(kIndLoc, _mask_fp_resource(self, subset))
                        kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, (cl.cltypes.long)(self.nProjSubset[timestep, subset].item()))
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
                    self.knlF.set_arg(kIndLoc, (cl.cltypes.ulong)(self.nMeasSubset[timestep, subset].item()))
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, (cl.cltypes.uint)(subset))
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, (cl.cltypes.int)(k))
                elif self.FPType == 5:
                    self.knlF.set_arg(kIndLoc, self.d_x[timestep][subset].data)
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, self.d_z[timestep][subset].data)
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
                    if self.useMaskFP:
                        self.knlF.set_arg(kIndLoc, _mask_fp_resource(self, subset))
                        kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, (cl.cltypes.long)(self.nProjSubset[timestep, subset].item()))
                    # if self.meanFP:
                elif self.FPType in [1, 2, 3]:
                    if self.useMaskFP:
                        self.knlF.set_arg(kIndLoc, _mask_fp_resource(self, subset))
                        kIndLoc += 1
                    if (self.CT or self.PET or self.SPECT) and self.listmode == 0:
                        self.knlF.set_arg(kIndLoc, (cl.cltypes.long)(self.nProjSubset[timestep, subset].item()))
                        kIndLoc += 1
                    if ((self.listmode == 0 or self.useIndexBasedReconstruction) and not (self.CT or self.SPECT)) or (not self.loadTOF and self.listmode > 0):
                        self.knlF.set_arg(kIndLoc, self.d_x[timestep][0].data)
                    else:
                        self.knlF.set_arg(kIndLoc, self.d_x[timestep][subset].data)
                    kIndLoc += 1
                    if (self.CT or self.PET or self.SPECT or (self.listmode > 0 and not self.useIndexBasedReconstruction)):
                        self.knlF.set_arg(kIndLoc, self.d_z[timestep][subset].data)
                    else:
                        self.knlF.set_arg(kIndLoc, self.d_z[timestep][0].data)
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
                    if not self.useImages and not self.useAF:
                        self.knlF.set_arg(kIndLoc, f.data)
                    else:
                        self.knlF.set_arg(kIndLoc, d_im)
                    kIndLoc += 1
                    if self.useAF:
                        self.knlF.set_arg(kIndLoc, yD)
                    else:
                        self.knlF.set_arg(kIndLoc, y.data)
                    kIndLoc += 1
                    if self.SPECT:
                        self.knlF.set_arg(kIndLoc, self.d_detectorVector[timestep][subset].data)
                        kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, (cl.cltypes.uchar)(self.no_norm))
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, (cl.cltypes.ulong)(self.nMeasSubset[timestep, subset].item()))
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, (cl.cltypes.uint)(subset))
                    kIndLoc += 1
                    self.knlF.set_arg(kIndLoc, (cl.cltypes.int)(k))
                cl.enqueue_nd_range_kernel(self.queue, self.knlF, self.globalSizeFP[timestep][subset], self.localSizeFP)
                self.queue.finish()
        if volumes > 0 and not(isinstance(f,list)):
            self.nMultiVolumes = volumes
        if self.useAF:
            af.device.unlock_array(y)
            if not self.useImages:
                af.device.unlock_array(f)
    return y

def backwardProjection(self, y, subset = -1, timestep = -1):
    if subset == -1:
        subset = self.subset
    if timestep == -1:
        timestep = self.timestep
    timestep = int(timestep)
    subset = int(subset)
    if self.useMetal:
        from omegatomo.projector.mps_backend import backward_projection_mps
        return backward_projection_mps(self, y, subset, timestep)
    if self.useCUDA and not self.useCuPy:
        raise ValueError('PyCUDA is no longer supported. Please use CuPy.')
    volumes = 0
    if self.projector_type in (6, 66):
        rows = int(getattr(self, 'measurement_nRowsD', self.nRowsD))
        cols = int(getattr(self, 'measurement_nColsD', self.nColsD))
        expected = rows * cols * int(self.nProjSubset[timestep, subset])
        if not self.useCUDA:
            import arrayfire as af
            ops = _type6_arrayfire_ops(self)
            count = y.elements()
            if count != expected:
                raise ValueError(f'Backprojection input has {count} elements; expected {expected}')
            outputs = []
            for volume in range(int(self.nMultiVolumes) + 1):
                output = af.data.constant(0., int(np.asarray(self.N).reshape(-1)[volume]))
                type6_backward(self, y, output, volume, subset, timestep, ops=ops)
                af.eval(output)
                outputs.append(output)
            return outputs[0] if int(self.nMultiVolumes) == 0 else outputs
        else:
            import torch
            ops = _type6_torch_ops(self)
            y = y.contiguous()
            count = y.numel()
            if count != expected:
                raise ValueError(f'Backprojection input has {count} elements; expected {expected}')
            outputs = []
            for volume in range(int(self.nMultiVolumes) + 1):
                output = torch.zeros(
                    int(np.asarray(self.N).reshape(-1)[volume]),
                    dtype=torch.float32,
                    device=y.device,
                )
                type6_backward(self, y, output, volume, subset, timestep, ops=ops)
                outputs.append(output)
            return outputs[0] if int(self.nMultiVolumes) == 0 else outputs
    else:
        if self.nMultiVolumes > 0:
            f = [None] * (self.nMultiVolumes + 1)
        if self.nMultiVolumes > 0 and not(isinstance(f,list)):
            volumes = self.nMultiVolumes
            self.nMultiVolumes = 0
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
                        yy = cp.zeros((self.nRowsD+1,self.nColsD+1,self.nProjSubset[timestep, subset].item()), dtype=cp.float32, order='F')
                        if self.useTorch:
                            yy[1:,1:,:] = yD.reshape((self.nRowsD,self.nColsD,self.nProjSubset[timestep, subset].item()), order='F')
                        else:
                            yy[1:,1:,:] = y.reshape((self.nRowsD,self.nColsD,self.nProjSubset[timestep, subset].item()), order='F')
                        yy = yy.cumsum(0)
                        yy = yy.cumsum(1)
                        yy = yy.ravel(order='F')
                    kIndLoc = self.kIndB
                    if self.BPType in [1, 2, 3]:
                        if (self.attenuation_correction and not self.CTAttenuation):
                            kIndLoc += (self.d_atten[subset],)
                        if self.useMaskFP:
                            kIndLoc += (_mask_fp_resource(self, subset),)
                        if self.useMaskBP:
                            kIndLoc += (self.d_maskBP,)
                        if (self.CT or self.PET or self.SPECT) and self.listmode == 0:
                            kIndLoc += ((self.nProjSubset[timestep, subset].item()),)
                        if ((self.listmode == 0 or self.useIndexBasedReconstruction) and not (self.CT or self.SPECT)) or (not self.loadTOF and self.listmode > 0):
                            kIndLoc += (self.d_x[timestep][0],)
                        else:
                            kIndLoc += (self.d_x[timestep][subset],)
                        if (self.CT or self.PET or self.SPECT or (self.listmode > 0 and not self.useIndexBasedReconstruction)):
                            kIndLoc += (self.d_z[timestep][subset],)
                        else:
                            kIndLoc += (self.d_z[timestep][0],)
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
                        if self.SPECT:
                            kIndLoc += (self.d_detectorVector[timestep][subset],)
                        kIndLoc += (cp.uint8(self.no_norm),)
                        kIndLoc += (cp.uint64(self.nMeasSubset[timestep, subset].item()),)
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
                                    array = cp.cuda.texture.CUDAarray(chl, self.nRowsD, self.nColsD, self.nProjSubset[timestep, subset].item())
                                    if self.useTorch:
                                        array.copy_from(yD.reshape((self.nProjSubset[timestep, subset].item(), self.nColsD, self.nRowsD)))
                                    else:
                                        array.copy_from(y.reshape((self.nProjSubset[timestep, subset].item(), self.nColsD, self.nRowsD)))
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
                                kIndLoc += (self.d_x[timestep][subset],)
                                kIndLoc += (self.d_z[timestep][subset],)
                                kIndLoc += (self.d_Sens,)
                            else:
                                kIndLoc += (self.d_x[timestep][subset],)
                                kIndLoc += (self.d_z[timestep][subset],)
                                # Precomputed geometry; only present when the kernel was built with -DGEOM5
                                if self.listmode == 0:
                                    kIndLoc += (self.d_geom5[timestep][subset],)
                                # if self.useTorch:
                                #     kIndLoc += (yD,)
                                #     kIndLoc += (fD,)
                                # else:
                                if self.useImages:
                                    chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
                                    array = cp.cuda.texture.CUDAarray(chl, self.nRowsD + 1, self.nColsD + 1, self.nProjSubset[timestep, subset].item())
                                    array.copy_from(yy.reshape((self.nProjSubset[timestep, subset].item(), self.nColsD + 1, self.nRowsD + 1)))
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
                                kIndLoc += (self.d_x[timestep][0],)
                            else:
                                kIndLoc += (self.d_x[timestep][subset],)
                            if (self.CT or self.PET or self.listmode > 0):
                                kIndLoc += (self.d_z[timestep][subset],)
                            else:
                                kIndLoc += (self.d_z[timestep][0],)
                            if self.useMaskFP:
                                kIndLoc += (_mask_fp_resource(self, subset),)
                            if self.useMaskBP:
                                kIndLoc += (self.d_maskBP,)
                            kIndLoc += (cp.int64(self.nProjSubset[timestep, subset].item()),)
                            if ((self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0):
                                kIndLoc += (self.d_xyindex[subset],)
                                kIndLoc += (self.d_zindex[subset],)
                            if (self.normalization_correction):
                                kIndLoc += (self.d_norm[subset],)
                            if (self.additionalCorrection):
                                kIndLoc += (self.d_corr[subset],)
                            kIndLoc += (self.d_Sens,)
                        kIndLoc += (cp.uint8(self.no_norm),)
                        if self.CT:
                            if self.useMaskBP:
                                kIndLoc += (self.d_maskBP,)
                            kIndLoc += (cp.int64(self.nProjSubset[timestep, subset].item()),)
                        else:
                            kIndLoc += (cp.uint64(self.nMeasSubset[timestep, subset].item()),)
                            kIndLoc += (cp.uint32(subset),)
                        kIndLoc += (cp.int32(k),)
                    self.knlB((self.globalSizeBP[timestep][subset][k][0] // self.localSizeBP[0], self.globalSizeBP[timestep][subset][k][1] // self.localSizeBP[1], self.globalSizeBP[timestep][subset][k][2]), (self.localSizeBP[0], self.localSizeBP[1], 1), kIndLoc)
            if self.useTorch:
                torch.cuda.synchronize()
        else:
            import pyopencl as cl
            from pyopencl.version import VERSION
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
                    if VERSION[0] > 2024 or (VERSION[0] == 2024 and VERSION[1] > 2):
                        d_im = cl.create_image(self.clctx, cl.mem_flags.READ_ONLY, imformat, shape=(self.nRowsD, self.nColsD, self.nProjSubset[timestep, subset].item()))
                    else:
                        d_im = cl.Image(self.clctx, cl.mem_flags.READ_ONLY, imformat, shape=(self.nRowsD, self.nColsD, self.nProjSubset[timestep, subset].item()))
                    if self.useAF:
                        cl.enqueue_copy(self.queue, d_im, yD, offset=(0), origin=(0,0,0), region=(self.nRowsD, self.nColsD, self.nProjSubset[timestep, subset].item()));
                    else:
                        cl.enqueue_copy(self.queue, d_im, y.data, offset=(0), origin=(0,0,0), region=(self.nRowsD, self.nColsD, self.nProjSubset[timestep, subset].item()));
                else:
                    if VERSION[0] > 2024 or (VERSION[0] == 2024 and VERSION[1] > 2):
                        d_im = cl.create_image(self.clctx, cl.mem_flags.READ_ONLY, imformat, shape=(self.nRowsD + 1, self.nColsD + 1, self.nProjSubset[timestep, subset].item()))
                    else:
                        d_im = cl.Image(self.clctx, cl.mem_flags.READ_ONLY, imformat, shape=(self.nRowsD + 1, self.nColsD + 1, self.nProjSubset[timestep, subset].item()))
                    y = af.moddims(y, self.nRowsD, d1=self.nColsD, d2=self.nProjSubset[timestep, subset].item())
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
                    cl.enqueue_copy(self.queue, d_im, yD, offset=(0), origin=(0,0,0), region=(self.nRowsD + 1, self.nColsD + 1, self.nProjSubset[timestep, subset].item()));
            
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
                    if self.useMaskFP:
                        self.knlB.set_arg(kIndLoc, _mask_fp_resource(self, subset))
                        kIndLoc += 1
                    if self.useMaskBP:
                        self.knlB.set_arg(kIndLoc, self.d_maskBP)
                        kIndLoc += 1
                    if (self.CT or self.PET or self.SPECT) and self.listmode == 0:
                        self.knlB.set_arg(kIndLoc, (cl.cltypes.long)(self.nProjSubset[timestep, subset].item()))
                        kIndLoc += 1
                    if ((self.listmode == 0 or self.useIndexBasedReconstruction) and not (self.CT or self.SPECT)) or (not self.loadTOF and self.listmode > 0):
                        self.knlB.set_arg(kIndLoc, self.d_x[timestep][0].data)
                    else:
                        self.knlB.set_arg(kIndLoc, self.d_x[timestep][subset].data)
                    kIndLoc += 1
                    if (self.CT or self.PET or self.SPECT or (self.listmode > 0 and not self.useIndexBasedReconstruction)):
                        self.knlB.set_arg(kIndLoc, self.d_z[timestep][subset].data)
                    else:
                        self.knlB.set_arg(kIndLoc, self.d_z[timestep][0].data)
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
                    if self.SPECT:
                        self.knlB.set_arg(kIndLoc, self.d_detectorVector[timestep][subset].data)
                        kIndLoc += 1
                    self.knlB.set_arg(kIndLoc, (cl.cltypes.uchar)(self.no_norm))
                    kIndLoc += 1
                    self.knlB.set_arg(kIndLoc, (cl.cltypes.ulong)(self.nMeasSubset[timestep, subset].item()))
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
                                self.knlB.set_arg(kIndLoc, self.d_x[timestep][0].data)
                            else:
                                self.knlB.set_arg(kIndLoc, self.d_x[timestep][subset].data)
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, self.d_z[timestep][subset].data)
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, self.d_Sens.data)
                            kIndLoc += 1
                        else:
                            if not self.loadTOF and self.listmode > 0:
                                self.knlB.set_arg(kIndLoc, self.d_x[timestep][0].data)
                            else:
                                self.knlB.set_arg(kIndLoc, self.d_x[timestep][subset].data)
                            kIndLoc += 1
                            self.knlB.set_arg(kIndLoc, self.d_z[timestep][subset].data)
                            kIndLoc += 1
                            # Precomputed geometry; only present when the kernel was built with -DGEOM5
                            if self.listmode == 0:
                                self.knlB.set_arg(kIndLoc, self.d_geom5[timestep][subset].data)
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
                            self.knlB.set_arg(kIndLoc, self.d_x[timestep][0].data)
                        else:
                            self.knlB.set_arg(kIndLoc, self.d_x[timestep][subset].data)
                        kIndLoc += 1
                        if (self.CT or self.PET or (self.listmode > 0 and not self.useIndexBasedReconstruction)):
                            self.knlB.set_arg(kIndLoc, self.d_z[timestep][subset].data)
                        else:
                            self.knlB.set_arg(kIndLoc, self.d_z[timestep][0].data)
                        kIndLoc += 1
                        if self.useMaskFP:
                            self.knlB.set_arg(kIndLoc, _mask_fp_resource(self, subset))
                            kIndLoc += 1
                        if self.useMaskBP:
                            self.knlB.set_arg(kIndLoc, self.d_maskBP)
                            kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, (cl.cltypes.ulong)(self.nProjSubset[timestep, subset].item()))
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
                        if self.useMaskBP:
                            self.knlB.set_arg(kIndLoc, self.d_maskBP)
                            kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, (cl.cltypes.ulong)(self.nProjSubset[timestep, subset].item()))
                        kIndLoc += 1
                    else:
                        self.knlB.set_arg(kIndLoc, (cl.cltypes.ulong)(self.nMeasSubset[timestep, subset].item()))
                        kIndLoc += 1
                        self.knlB.set_arg(kIndLoc, (cl.cltypes.uint)(subset))
                        kIndLoc += 1
                    self.knlB.set_arg(kIndLoc, (cl.cltypes.int)(k))
                            
                cl.enqueue_nd_range_kernel(self.queue, self.knlB, self.globalSizeBP[timestep][subset][k], self.localSizeBP)
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
        if not(isinstance(f,list)) and volumes > 0:
            self.nMultiVolumes = volumes
        if self.use_psf:
            if self.nMultiVolumes > 0:
                f[k] = self.computeConvolution(f[k])
            else:
                f = self.computeConvolution(f)
    return f


from typing import Any

def _proj6_center_match(value: Any, target_rows: int, target_cols: int, *, torch_mode: bool) -> Any:
    """Center crop or zero pad detector rows and columns."""
    if torch_mode:
        import torch

        current_cols, current_rows = int(value.shape[-2]), int(value.shape[-1])
        if current_rows > target_rows:
            before = (current_rows - target_rows) // 2
            value = value[..., before:before + target_rows]
        elif current_rows < target_rows:
            total = target_rows - current_rows
            before = total // 2
            value = torch.cat((
                torch.zeros((*value.shape[:-1], before), dtype=value.dtype, device=value.device),
                value,
                torch.zeros((*value.shape[:-1], total - before), dtype=value.dtype, device=value.device),
            ), dim=-1)
        if current_cols > target_cols:
            before = (current_cols - target_cols) // 2
            value = value[..., before:before + target_cols, :]
        elif current_cols < target_cols:
            total = target_cols - current_cols
            before = total // 2
            value = torch.cat((
                torch.zeros((*value.shape[:-2], before, value.shape[-1]), dtype=value.dtype, device=value.device),
                value,
                torch.zeros((*value.shape[:-2], total - before, value.shape[-1]), dtype=value.dtype, device=value.device),
            ), dim=-2)
        return value
    out = value
    for axis, wanted in enumerate((target_rows, target_cols)):
        delta = wanted - out.shape[axis]
        if delta < 0:
            before = (-delta) // 2
            index = [slice(None)] * out.ndim
            index[axis] = slice(before, before + wanted)
            out = out[tuple(index)]
        elif delta > 0:
            before = delta // 2
            pads = [(0, 0)] * out.ndim
            pads[axis] = (before, delta - before)
            out = np.pad(out, pads)
    return out


def _proj6_resize_numpy(value: Any, source: tuple[int, int], target: tuple[int, int], physical: tuple[int, int], mode: str, direction: str, expected_frames: int | None, strict: bool) -> Any:
    if isinstance(value, (list, tuple)):
        return type(value)(_proj6_resize_numpy(item, source, target, physical, mode, direction, None, strict) for item in value)
    arr = np.asarray(value)
    if arr.size == 0 or arr.ndim == 0:
        return value
    flat = arr.ndim == 1
    source_pixels = int(np.prod(source))
    target_pixels = int(np.prod(target))
    if flat:
        if expected_frames is not None and arr.size == target_pixels * expected_frames:
            return value
        if arr.size % source_pixels:
            if strict:
                raise ValueError(f"projection vector has {arr.size} values; it is not a whole {source} detector stack")
            return value
        arr = arr.reshape((*source, arr.size // source_pixels), order='F')
    elif arr.ndim < 2:
        return value
    elif tuple(arr.shape[:2]) != source:
        if strict:
            raise ValueError(f"projection array has detector shape {tuple(arr.shape[:2])}; expected {source}")
        return value
    order = 0 if mode == 'mask' else 1
    from skimage.transform import resize
    if direction == 'measurement_to_image':
        out = _proj6_center_match(arr, physical[0], physical[1], torch_mode=False)
        if tuple(out.shape[:2]) != target:
            out = resize(out, target + tuple(out.shape[2:]), order=order, mode='reflect', anti_aliasing=order != 0, preserve_range=True)
    else:
        out = arr
        if tuple(out.shape[:2]) != physical:
            out = resize(out, physical + tuple(out.shape[2:]), order=order, mode='reflect', anti_aliasing=order != 0, preserve_range=True)
        out = _proj6_center_match(out, target[0], target[1], torch_mode=False)
    out = out.astype(arr.dtype, copy=False)
    return np.asfortranarray(out).ravel(order='F') if flat else out


def _proj6_resize_torch(value: Any, source: tuple[int, int], target: tuple[int, int], physical: tuple[int, int], mode: str, direction: str, expected_frames: int | None, strict: bool) -> Any:
    """Resize flat MPS/CPU torch projections without a host round trip."""
    import torch
    import torch.nn.functional as functional

    flat = value.ndim == 1
    if not flat and value.ndim == 3:
        out = value
    elif not flat:
        if strict:
            raise ValueError('custom-operator projection tensors must be one-dimensional or (views, cols, rows)')
        return value
    else:
        source_pixels = int(np.prod(source))
        target_pixels = int(np.prod(target))
        if expected_frames is not None and value.numel() == target_pixels * expected_frames:
            return value
        if value.numel() % source_pixels:
            if strict:
                raise ValueError(f"projection tensor has {value.numel()} values; it is not a whole {source} detector stack")
            return value
        out = value.reshape(value.numel() // source_pixels, source[1], source[0])
    if tuple(int(x) for x in out.shape[-2:]) != (source[1], source[0]):
        if strict:
            raise ValueError(f"projection tensor has detector shape {tuple(out.shape[-2:])}; expected {(source[1], source[0])}")
        return value
    interp = 'nearest' if mode == 'mask' else 'bilinear'
    def interpolate(array: Any, size: tuple[int, int]) -> Any:
        kwargs = {'size': size, 'mode': interp}
        if interp == 'bilinear':
            kwargs['align_corners'] = False
        return functional.interpolate(array.unsqueeze(1), **kwargs).squeeze(1)
    if direction == 'measurement_to_image':
        out = _proj6_center_match(out, physical[0], physical[1], torch_mode=True)
        if tuple(int(x) for x in out.shape[-2:]) != (target[1], target[0]):
            out = interpolate(out, (target[1], target[0]))
    else:
        if tuple(int(x) for x in out.shape[-2:]) != (physical[1], physical[0]):
            out = interpolate(out, (physical[1], physical[0]))
        out = _proj6_center_match(out, target[0], target[1], torch_mode=True)
    return out.reshape(-1) if flat else out


# Resize, resample, pad and/or crop measurement-domain arrays to match resolution and size of image domain for SPECT rotation projector.
def resample_resize_proj6(value: Any, options: Any, volume: int = 0, *, direction: str = 'measurement_to_image', mode: str = 'continuous', expected_frames: int | None = None, strict: bool = False) -> Any:
    """Transform a type-6 projection between its public and internal grids.

    The projector options are the metadata: public detector dimensions remain
    ``measurement_nRowsD``/``measurement_nColsD`` (or ``nRowsD``/``nColsD``),
    while the rotation projector grid is ``Nx[volume]`` by ``Nz[volume]``.
    """
    if mode not in {'continuous', 'mask'}:
        raise ValueError("mode must be 'continuous' or 'mask'")
    if direction not in {'measurement_to_image', 'image_to_measurement', 'measurement'}:
        raise ValueError("direction must be 'measurement_to_image', 'image_to_measurement', or 'measurement'")
    if direction == 'measurement':
        return value
    measurement = (int(getattr(options, 'measurement_nRowsD', options.nRowsD)), int(getattr(options, 'measurement_nColsD', options.nColsD)))
    def dimension(name: str) -> int:
        values = np.asarray(getattr(options, name)).reshape(-1)
        return int(values[min(volume, values.size - 1)])

    image = (dimension('Nx'), dimension('Nz'))
    try:
        def scalar(name: str) -> float:
            values = np.asarray(getattr(options, name)).reshape(-1)
            if values.size == 0:
                raise ValueError
            return float(values[min(volume, values.size - 1)])

        pitch_x, pitch_y = scalar('dPitchX'), scalar('dPitchY')
        fov_x, fov_z = scalar('FOVa_x'), scalar('axial_fov')
        if min(pitch_x, pitch_y, fov_x, fov_z) <= 0.0:
            raise ValueError
        physical = (max(1, int(round(fov_x / pitch_x))), max(1, int(round(fov_z / pitch_y))))
    except (AttributeError, TypeError, ValueError):
        physical = image
    source, target = (measurement, image) if direction == 'measurement_to_image' else (image, measurement)
    try:
        import torch
    except ImportError:
        torch = None
    if torch is not None and isinstance(value, torch.Tensor):
        return _proj6_resize_torch(value, source, target, physical, mode, direction, expected_frames, strict)
    return _proj6_resize_numpy(value, source, target, physical, mode, direction, expected_frames, strict)


def _type6_volume_view_value(self: Any, name: str, volume: int, view: int) -> int | float:
    """Read a custom type-6 geometry item without mixing volume and view axes."""
    values = np.asarray(getattr(self, name))
    if values.ndim == 2 and values.size:
        if volume >= values.shape[0] or view >= values.shape[1]:
            raise IndexError(f'{name} has no geometry for volume {volume}, view {view}')
        return values[volume, view].item()
    values = values.reshape(-1)
    if view >= values.size:
        raise IndexError(f'{name} has no geometry for view {view}')
    return values[view].item()


def _type6_volume_kernel(self: Any, volume: int) -> Any:
    """Select the CDRF sampled on the requested volume's pixel grid."""
    filters = self.gFilter
    if isinstance(filters, (list, tuple)) and len(filters):
        if volume >= len(filters):
            raise IndexError(f'gFilter has no filter for volume {volume}')
        return filters[volume]
    return filters


def _type6_path_weight(self: Any, volume: int, view: int, nx: int) -> float:
    """Physical x-segment weight for one type-6 volume/view contribution."""
    dx = float(np.asarray(self.dx).reshape(-1)[volume])
    if getattr(self, 'useTotLength', True):
        lengths = np.asarray(getattr(self, 'type6TotalLength', np.empty(0))).reshape(-1)
        if view >= lengths.size or lengths[view] <= 0.:
            raise ValueError(f'Type-6 total ray length is missing or invalid for view {view}')
        return dx / float(lengths[view])
    # Preserve the legacy, local-FOV normalization when explicitly requested.
    retained = max(1, nx - max(int(_type6_volume_view_value(self, 'blurPlanes', volume, view)), 0))
    return 1. / retained


def _type6_torch_resource(self: Any, name: str, timestep: int, subset: int, *, fallback: Any, dtype: Any, device: Any) -> Any:
    """Return a type-6 correction buffer on the current Torch device."""
    import torch

    value = getattr(self, name, None)
    if isinstance(value, list):
        try:
            value = value[timestep][subset]
        except (IndexError, TypeError):
            value = None
    if value is None:
        value = fallback
    if isinstance(value, torch.Tensor):
        return value.to(device=device, dtype=dtype)
    return torch.as_tensor(np.ascontiguousarray(np.asarray(value, dtype=np.float32)), device=device, dtype=dtype)


def _type6_torch_detector_vector(self: Any, timestep: int, subset: int, device: Any) -> Any:
    """Return detector-head IDs in the already subset-ordered view order."""
    import torch

    value = getattr(self, 'd_detectorVector', None)
    if isinstance(value, list):
        value = value[timestep][subset]
    if value is None:
        frames = getattr(self, 'DetectorVectorFrames', None)
        if isinstance(frames, list) and len(frames) == int(self.Nt):
            frame = np.asarray(frames[timestep], dtype=np.uint32).reshape(-1)
            offsets = np.concatenate(([0], np.cumsum(self.nProjSubset[timestep], dtype=np.int64)))
            value = frame[int(offsets[subset]):int(offsets[subset + 1])]
        else:
            start = int(np.sum(self.nProjSubset[:timestep, :]) + np.sum(self.nProjSubset[timestep, :subset]))
            stop = start + int(self.nProjSubset[timestep, subset])
            value = np.asarray(self.DetectorVector, dtype=np.uint32).reshape(-1)[start:stop]
    return _type6_torch_resource(
        self, '_unused_detector_vector', timestep, subset,
        fallback=value, dtype=torch.long, device=device,
    )


def _type6_torch_measurement_weights(self: Any, data: Any, projections: int, timestep: int, subset: int) -> Any:
    """Apply type-6 measurement-domain corrections on any Torch device."""
    import torch

    rows = int(getattr(self, 'measurement_nRowsD', self.nRowsD))
    cols = int(getattr(self, 'measurement_nColsD', self.nColsD))
    pixels = rows * cols
    device = data.device
    detector = _type6_torch_detector_vector(self, timestep, subset, device)
    weights = torch.ones_like(data)

    def expand(values: Any) -> Any | None:
        if values.numel() == 0:
            return None
        values = values.to(dtype=torch.float32, device=device)
        if values.numel() == pixels:
            return values.reshape(cols, rows).unsqueeze(0).expand(projections, -1, -1)
        if values.numel() == int(self.nHeads) * pixels:
            return values.reshape(int(self.nHeads), cols, rows).index_select(0, detector)
        return values.reshape(projections, cols, rows)

    if self.useMaskFP:
        value = _type6_torch_resource(
            self, 'd_maskFP', timestep, subset,
            fallback=getattr(self, 'maskFP', np.empty(0)), dtype=torch.float32, device=device,
        )
        value = expand(value)
        if value is not None:
            weights *= value
    if self.attenuation_correction and not self.CTAttenuation:
        value = _type6_torch_resource(
            self, 'd_attenuation', timestep, subset,
            fallback=getattr(self, 'vaimennus', np.empty(0)), dtype=torch.float32, device=device,
        )
        value = expand(value)
        if value is not None:
            weights *= value
    if self.normalization_correction:
        value = _type6_torch_resource(
            self, 'd_norm', timestep, subset,
            fallback=getattr(self, 'normalization', np.empty(0)), dtype=torch.float32, device=device,
        )
        value = expand(value)
        if value is not None:
            weights *= value
    return weights


def _type6_arrayfire_ops(self: Any):
    """Collect ArrayFire operations for the common type-6 projector.

    The shared implementation uses ``(view, z, x)`` measurements and
    ``(z, y, x)`` images.  ArrayFire's established type-6 layout is instead
    ``(x, y, z)`` / ``(row, column, view)``, so only these callables carry
    that layout conversion; the FP and BP physics remain backend-neutral.
    """
    import arrayfire as af
    from types import SimpleNamespace

    def host_resource(name: str, timestep: int, subset: int, fallback: Any) -> np.ndarray:
        # OpenCL initialization uploads ``d_*`` resources as OpenCL buffers
        # or images.  Keep this ArrayFire path on its corresponding host
        # correction arrays, then upload the composed result through AF.
        return np.asarray(fallback, dtype=np.float32).ravel(order='F')

    def from_vector(values: np.ndarray, *shape: int) -> Any:
        return af.data.moddims(af.interop.np_to_af_array(np.asfortranarray(values)), *shape)

    def canonical_image(values: Any, shape: tuple[int, int, int]) -> Any:
        nz, ny, nx = shape
        flat = np.asarray(values, dtype=np.float32).reshape(-1)
        if flat.size != nz * ny * nx:
            raise ValueError(f'Expected an image resource with {nz * ny * nx} elements, got {flat.size}')
        return af.data.reorder(from_vector(flat, nx, ny, nz), 2, 1, 0)

    def measurement_weights(data: Any, projections: int, timestep: int, subset: int) -> Any:
        rows = int(getattr(self, 'measurement_nRowsD', self.nRowsD))
        cols = int(getattr(self, 'measurement_nColsD', self.nColsD))
        pixels = rows * cols
        weights = np.ones((projections, cols, rows), dtype=np.float32)
        frames = getattr(self, 'DetectorVectorFrames', None)
        if isinstance(frames, list) and len(frames) == int(self.Nt):
            offsets = np.concatenate(([0], np.cumsum(self.nProjSubset[timestep], dtype=np.int64)))
            detector = np.asarray(frames[timestep], dtype=np.uint32).reshape(-1)[int(offsets[subset]):int(offsets[subset + 1])]
        else:
            start = int(np.sum(self.nProjSubset[:timestep, :]) + np.sum(self.nProjSubset[timestep, :subset]))
            detector = np.asarray(self.DetectorVector, dtype=np.uint32).reshape(-1)[start:start + projections]
        detector = np.asarray(detector, dtype=np.int64).reshape(-1)

        def expand(name: str, fallback: Any) -> np.ndarray | None:
            values = host_resource(name, timestep, subset, fallback)
            if values.size == 0:
                return None
            if values.size == pixels:
                return np.broadcast_to(values.reshape(cols, rows), weights.shape)
            if values.size == int(self.nHeads) * pixels:
                return values.reshape(int(self.nHeads), cols, rows)[detector]
            if values.size != projections * pixels:
                index = timestep * int(self.subsets) + subset
                start = int(self.nTotMeas[index])
                stop = int(self.nTotMeas[index + 1])
                if stop - start == projections * pixels:
                    values = values[start:stop]
            return values.reshape(projections, cols, rows)

        if self.useMaskFP:
            value = expand('d_maskFP', getattr(self, 'maskFP', np.empty(0)))
            if value is not None:
                weights *= value
        if self.attenuation_correction and not self.CTAttenuation:
            value = expand('d_attenuation', getattr(self, 'vaimennus', np.empty(0)))
            if value is not None:
                weights *= value
        if self.normalization_correction:
            value = expand('d_norm', getattr(self, 'normalization', np.empty(0)))
            if value is not None:
                weights *= value
        return from_vector(weights.ravel(order='F'), projections, cols, rows)

    def shift_axis_one(value: Any, amount: int) -> Any:
        shifted = af.data.constant(0., int(value.shape[0]), int(value.shape[1]))
        width = int(value.shape[1])
        if amount >= 0:
            if amount < width:
                shifted[:, amount:] = value[:, :width - amount]
        elif -amount < width:
            shifted[:, :width + amount] = value[:, -amount:]
        return shifted

    def shift_kernel(kernel: Any, plane: int, nx: int) -> Any:
        shifted = af.data.constant(0., int(kernel.shape[0]), int(kernel.shape[1]), int(kernel.shape[2]))
        depth = int(kernel.shape[2])
        if plane >= 0:
            if plane < depth:
                shifted[:, :, plane:] = kernel[:, :, :depth - plane]
        elif -plane < depth:
            shifted[:, :, :depth + plane] = kernel[:, :, -plane:]
        return shifted[:, :, :nx]

    def shift_projection(value: Any, pixels: float) -> Any:
        whole = int(np.floor(pixels))
        fraction = pixels - whole
        shifted = shift_axis_one(value, whole)
        return shifted if fraction == 0.0 else shifted * (1.0 - fraction) + shift_axis_one(value, whole + 1) * fraction

    def shift_image_y(value: Any, pixels: int) -> Any:
        shifted = af.data.constant(0., int(value.shape[0]), int(value.shape[1]), int(value.shape[2]))
        width = int(value.shape[1])
        if pixels >= 0:
            if pixels < width:
                shifted[:, pixels:, :] = value[:, :width - pixels, :]
        elif -pixels < width:
            shifted[:, :width + pixels, :] = value[:, -pixels:, :]
        return shifted

    def resize(value: Any, shape: tuple[int, int]) -> Any:
        return af.image.resize(value, int(shape[0]), int(shape[1]), method=af.INTERP.BILINEAR)

    def resize_stack(value: Any, shape: tuple[int, int]) -> Any:
        # AF resize acts on dimensions zero and one.  Present the detector
        # stack as (x, z, view), then restore the shared (view, z, x) form.
        resized = af.image.resize(
            af.data.reorder(value, 2, 1, 0), int(shape[1]), int(shape[0]),
            method=af.INTERP.BILINEAR,
        )
        return af.data.reorder(resized, 2, 1, 0)

    def center_match(value: Any, rows: int, cols: int) -> Any:
        views, current_cols, current_rows = (int(value.shape[0]), int(value.shape[1]), int(value.shape[2]))
        if current_rows != rows:
            matched = af.data.constant(0., views, current_cols, rows)
            length = min(current_rows, rows)
            source = (current_rows - length) // 2
            destination = (rows - length) // 2
            matched[:, :, destination:destination + length] = value[:, :, source:source + length]
            value = matched
        if current_cols != cols:
            matched = af.data.constant(0., views, cols, int(value.shape[2]))
            length = min(current_cols, cols)
            source = (current_cols - length) // 2
            destination = (cols - length) // 2
            matched[:, destination:destination + length, :] = value[:, source:source + length, :]
            value = matched
        return value

    def physical_grid(volume: int) -> tuple[int, int]:
        try:
            def scalar(name: str) -> float:
                values = np.asarray(getattr(self, name)).reshape(-1)
                if values.size == 0:
                    raise ValueError
                return float(values[min(volume, values.size - 1)])

            pitch_x, pitch_y = scalar('dPitchX'), scalar('dPitchY')
            fov_x, fov_z = scalar('FOVa_x'), scalar('axial_fov')
            if min(pitch_x, pitch_y, fov_x, fov_z) <= 0.0:
                raise ValueError
            return max(1, int(round(fov_x / pitch_x))), max(1, int(round(fov_z / pitch_y)))
        except (AttributeError, TypeError, ValueError):
            return int(self.Nx[volume]), int(self.Nz[volume])

    def resample(value: Any, volume: int, direction: str) -> Any:
        rows = int(getattr(self, 'measurement_nRowsD', self.nRowsD))
        cols = int(getattr(self, 'measurement_nColsD', self.nColsD))
        nx, nz = int(self.Nx[volume]), int(self.Nz[volume])
        physical_rows, physical_cols = physical_grid(volume)
        if direction == 'measurement_to_image':
            value = center_match(value, physical_rows, physical_cols)
            if int(value.shape[1]) != nz or int(value.shape[2]) != nx:
                value = resize_stack(value, (nz, nx))
            return value
        if direction == 'image_to_measurement':
            if int(value.shape[1]) != physical_cols or int(value.shape[2]) != physical_rows:
                value = resize_stack(value, (physical_cols, physical_rows))
            return center_match(value, rows, cols)
        raise ValueError(f'Unknown type-6 resampling direction: {direction}')

    def attenuation(value: Any, mu: Any, dx: float) -> Any:
        accumulate = getattr(af, 'accum', None)
        if accumulate is None:
            accumulate = af.algorithm.accum
        return value * af.exp(accumulate(mu, 2) * -dx)

    def apply_bp_mask(image: Any, volume: int, timestep: int, subset: int) -> Any:
        if not self.useMaskBP:
            return image
        nz, ny, nx = int(self.Nz[volume]), int(self.Ny[volume]), int(self.Nx[volume])
        values = host_resource('d_maskBP', timestep, subset, getattr(self, 'maskBP', np.empty(0)))
        if values.size == 0:
            return image
        if values.size == nx * ny:
            mask = np.broadcast_to(values.reshape(ny, nx), (nz, ny, nx))
        else:
            mask_z = int(self.maskBPZ)
            mask = values.reshape(mask_z, ny, nx)
            if mask_z != nz:
                mask = mask[np.arange(nz) * mask_z // nz]
        return image * canonical_image(mask, (nz, ny, nx))

    return SimpleNamespace(
        zeros=lambda shape, like: af.data.constant(0., *shape),
        reshape_image=lambda value, shape: af.data.reorder(af.data.moddims(value, shape[2], shape[1], shape[0]), 2, 1, 0),
        reshape_measurement=lambda value, shape: af.data.reorder(af.data.moddims(value, shape[2], shape[1], shape[0]), 2, 1, 0),
        flatten_image=lambda value: af.flat(af.data.reorder(value, 2, 1, 0)),
        flatten_measurement=lambda value: af.flat(af.data.reorder(value, 2, 1, 0)),
        rotate=lambda value, angle: af.data.reorder(
            af.image.rotate(af.data.reorder(value, 2, 1, 0), -angle * np.pi / 180.0, method=af.INTERP.BILINEAR), 2, 1, 0),
        attenuation=attenuation,
        shift_kernel=shift_kernel,
        blur=lambda value, kernel: af.data.reorder(
            af.signal.convolve2(af.data.reorder(value, 1, 0, 2), kernel), 1, 0, 2),
        resize=resize,
        shift_projection=shift_projection,
        shift_image_y=shift_image_y,
        sum_x=lambda value: af.data.moddims(af.algorithm.sum(value, 2), int(value.shape[0]), int(value.shape[1])),
        smear=lambda value, nx: af.tile(value, 1, 1, nx),
        set_view=lambda stack, index, value: stack.__setitem__((index, slice(None), slice(None)), value),
        add_to=lambda target, value: target.__setitem__(slice(None), target + value),
        multiply=lambda left, right: left * right,
        resample=resample,
        attenuation_image=lambda volume, timestep, subset, like: canonical_image(
            host_resource('d_attenuation_image', timestep, subset, getattr(self, 'vaimennus', np.empty(0))),
            (int(self.Nz[volume]), int(self.Ny[volume]), int(self.Nx[volume])),
        ),
        kernel=lambda volume, timestep, subset, like: (
            getattr(self, 'd_gFilter', [])[volume]
            if isinstance(getattr(self, 'd_gFilter', None), list)
            and volume < len(getattr(self, 'd_gFilter', []))
            else getattr(self, 'd_gFilter', af.interop.np_to_af_array(_type6_volume_kernel(self, volume)))
        ),
        weights=measurement_weights,
        apply_bp_mask=apply_bp_mask,
        device=lambda value: None,
        synchronize=lambda device: af.sync(),
    )


def _type6_torch_ops(self: Any):
    """Backend operations for the common type-6 Torch implementation.

    The same callables operate on MPS and CUDA tensors; the tensor's device
    selects the backend.  ArrayFire keeps its own callables at its entry point.
    """
    import torch
    import torch.nn.functional as functional
    from types import SimpleNamespace
    from torchvision.transforms import InterpolationMode
    from torchvision.transforms.functional import affine, rotate

    def shift_kernel(kernel: Any, plane: int, nx: int) -> Any:
        shifted = torch.zeros_like(kernel)
        depth = int(kernel.shape[2])
        if plane >= 0:
            if plane < depth:
                shifted[:, :, plane:] = kernel[:, :, :depth - plane]
        elif -plane < depth:
            shifted[:, :, :depth + plane] = kernel[:, :, -plane:]
        return shifted[:, :, :nx]

    def blur(value: Any, kernel: Any) -> Any:
        kernel = kernel.permute(2, 0, 1).contiguous().unsqueeze(1)
        return functional.conv2d(
            value.permute(2, 1, 0).unsqueeze(0),
            kernel,
            padding=(int(kernel.shape[2]) // 2, int(kernel.shape[3]) // 2),
            groups=int(kernel.shape[0]),
        ).squeeze(0).permute(2, 1, 0)

    def shift_projection(value: Any, pixels: float) -> Any:
        if pixels == 0.0:
            return value
        return affine(
            value.unsqueeze(0), angle=0.0, translate=[pixels, 0.0], scale=1.0,
            shear=[0.0, 0.0], interpolation=InterpolationMode.BILINEAR, fill=0.0,
        ).squeeze(0)

    def shift_image_y(value: Any, pixels: int) -> Any:
        shifted = torch.zeros_like(value)
        width = int(value.shape[1])
        if pixels >= 0:
            if pixels < width:
                shifted[:, pixels:, :] = value[:, :width - pixels, :]
        elif -pixels < width:
            shifted[:, :width + pixels, :] = value[:, -pixels:, :]
        return shifted

    def apply_bp_mask(image: Any, volume: int, timestep: int, subset: int) -> Any:
        if not self.useMaskBP:
            return image
        nz, ny, nx = (int(self.Nz[volume]), int(self.Ny[volume]), int(self.Nx[volume]))
        mask = _type6_torch_resource(
            self, 'd_maskBP', timestep, subset,
            fallback=np.asarray(getattr(self, 'maskBP', np.empty(0))).ravel(order='F'),
            dtype=torch.float32, device=image.device,
        )
        if mask.numel() == 0:
            return image
        if mask.numel() == nx * ny:
            mask = mask.reshape(ny, nx).unsqueeze(0).expand(nz, -1, -1)
        else:
            mask = mask.reshape(int(self.maskBPZ), ny, nx)
            if int(self.maskBPZ) != nz:
                mask = functional.interpolate(
                    mask.unsqueeze(0).unsqueeze(0), size=(nz, ny, nx), mode='nearest'
                ).squeeze(0).squeeze(0)
        return image * mask

    return SimpleNamespace(
        zeros=lambda shape, like: torch.zeros(shape, dtype=like.dtype, device=like.device),
        reshape_image=lambda value, shape: value.reshape(shape),
        reshape_measurement=lambda value, shape: value.reshape(shape),
        flatten_image=lambda value: value.reshape(-1),
        flatten_measurement=lambda value: value.reshape(-1),
        rotate=lambda value, angle: rotate(value, angle, interpolation=InterpolationMode.BILINEAR, fill=0.0),
        attenuation=lambda value, mu, dx: value * torch.exp(-torch.cumsum(mu, dim=2) * dx),
        shift_kernel=shift_kernel,
        blur=blur,
        resize=lambda value, shape: functional.interpolate(
            value.unsqueeze(0).unsqueeze(0), size=shape, mode='bilinear', align_corners=False,
        ).squeeze(0).squeeze(0),
        shift_projection=shift_projection,
        shift_image_y=shift_image_y,
        sum_x=lambda value: value.sum(dim=2),
        smear=lambda value, nx: value.unsqueeze(2).expand(-1, -1, nx),
        set_view=lambda stack, index, value: stack.__setitem__(index, value),
        add_to=lambda target, value: target.add_(value),
        multiply=lambda left, right: left * right,
        resample=lambda value, volume, direction: resample_resize_proj6(value, self, volume, direction=direction),
        attenuation_image=lambda volume, timestep, subset, like: _type6_torch_resource(
            self, 'd_attenuation_image', timestep, subset,
            fallback=getattr(self, 'vaimennus', np.empty(0)), dtype=like.dtype, device=like.device,
        ),
        kernel=lambda volume, timestep, subset, like: (
            getattr(self, 'd_gFilter', [])[volume].to(device=like.device, dtype=like.dtype)
            if isinstance(getattr(self, 'd_gFilter', None), list)
            and volume < len(getattr(self, 'd_gFilter', []))
            else _type6_torch_resource(
                self, 'd_gFilter', timestep, subset,
                fallback=_type6_volume_kernel(self, volume), dtype=like.dtype, device=like.device,
            )
        ),
        weights=lambda data, projections, timestep, subset: _type6_torch_measurement_weights(
            self, data, projections, timestep, subset,
        ),
        apply_bp_mask=apply_bp_mask,
        device=lambda value: value.device,
        synchronize=lambda device: torch.mps.synchronize() if device.type == 'mps' else None,
    )


def type6_forward(self: Any, image: Any, output: Any, volume: int, subset: int, timestep: int, *, ops: Any | None = None) -> None:
    """SPECT rotation projector shared logic."""
    ops = _type6_torch_ops(self) if ops is None else ops
    projections = int(self.nProjSubset[timestep][subset])
    start = int(np.sum(self.nProjSubset[:timestep, :]) + np.sum(self.nProjSubset[timestep, :subset]))
    nx, ny, nz = int(self.Nx[volume]), int(self.Ny[volume]), int(self.Nz[volume])
    source = ops.reshape_image(image, (nz, ny, nx))
    attenuation = None
    if self.attenuation_correction and self.CTAttenuation and volume == 0:
        attenuation = ops.reshape_image(ops.attenuation_image(volume, timestep, subset, image), (nz, ny, nx))
    kernel = ops.kernel(volume, timestep, subset, image)
    image_views = ops.zeros((projections, nz, nx), image)
    for local_view in range(projections):
        view = start + local_view
        angle = float(self.swivelAngles[view])
        rotated = ops.rotate(source, angle)
        rotated = ops.shift_image_y(rotated, -int(_type6_volume_view_value(self, 'blurPlanes2', volume, view)))
        if attenuation is not None:
            attenuation_rotated = ops.rotate(attenuation, angle)
            attenuation_rotated = ops.shift_image_y(
                attenuation_rotated, -int(_type6_volume_view_value(self, 'blurPlanes2', volume, view))
            )
            rotated = ops.attenuation(rotated, attenuation_rotated, float(self.dx[volume]))
        depth_shift = int(_type6_volume_view_value(self, 'blurPlanes', volume, view))
        shifted_kernel = ops.shift_kernel(kernel, depth_shift, nx)
        blurred = ops.blur(rotated, shifted_kernel)
        projection = ops.sum_x(blurred) * _type6_path_weight(self, volume, view, nx)
        if int(projection.shape[0]) != nz or int(projection.shape[1]) != nx:
            projection = ops.resize(projection, (nz, nx))
        ops.set_view(image_views, local_view, projection)
        if (local_view + 1) % 16 == 0:
            ops.synchronize(ops.device(image))
    measurement = ops.resample(image_views, volume, 'image_to_measurement')
    measurement = ops.multiply(measurement, ops.weights(measurement, projections, timestep, subset))
    ops.add_to(output, ops.flatten_measurement(measurement))


def type6_backward(self: Any, y: Any, output: Any, volume: int, subset: int, timestep: int, *, ops: Any | None = None) -> None:
    """SPECT rotation projector shared logic."""
    ops = _type6_torch_ops(self) if ops is None else ops
    projections = int(self.nProjSubset[timestep][subset])
    rows = int(getattr(self, 'measurement_nRowsD', self.nRowsD))
    cols = int(getattr(self, 'measurement_nColsD', self.nColsD))
    data = ops.reshape_measurement(y, (projections, cols, rows))
    data = ops.multiply(data, ops.weights(data, projections, timestep, subset))
    data = ops.resample(data, volume, 'measurement_to_image')
    start = int(np.sum(self.nProjSubset[:timestep, :]) + np.sum(self.nProjSubset[timestep, :subset]))
    nx, ny, nz = int(self.Nx[volume]), int(self.Ny[volume]), int(self.Nz[volume])
    attenuation = None
    if self.attenuation_correction and self.CTAttenuation and volume == 0:
        attenuation = ops.reshape_image(ops.attenuation_image(volume, timestep, subset, y), (nz, ny, nx))
    kernel = ops.kernel(volume, timestep, subset, y)
    image = ops.zeros((nz, ny, nx), y)
    for local_view in range(projections):
        view = start + local_view
        angle = float(self.swivelAngles[view])
        projection = data[local_view]
        if int(projection.shape[0]) != nz or int(projection.shape[1]) != ny:
            projection = ops.resize(projection, (nz, ny))
        depth_shift = int(_type6_volume_view_value(self, 'blurPlanes', volume, view))
        smeared = ops.smear(projection, nx) * _type6_path_weight(self, volume, view, nx)
        if attenuation is not None:
            attenuation_rotated = ops.rotate(attenuation, angle)
            attenuation_rotated = ops.shift_image_y(
                attenuation_rotated, -int(_type6_volume_view_value(self, 'blurPlanes2', volume, view))
            )
            smeared = ops.attenuation(smeared, attenuation_rotated, float(self.dx[volume]))
        rotated = ops.blur(smeared, ops.shift_kernel(kernel, depth_shift, nx))
        rotated = ops.shift_image_y(rotated, int(_type6_volume_view_value(self, 'blurPlanes2', volume, view)))
        ops.add_to(image, ops.rotate(rotated, -angle))
        if (local_view + 1) % 16 == 0:
            ops.synchronize(ops.device(y))
    if self.useMaskBP:
        image = ops.apply_bp_mask(image, volume, timestep, subset)
    ops.add_to(output, ops.flatten_image(image))
