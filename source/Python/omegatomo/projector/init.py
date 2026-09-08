# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 13:17:22 2025
"""
import numpy as np


def _coordinate_slice(self, name, timestep, subset, stride):
    """Return one frame/subset coordinate slice in kernel order."""
    frames = getattr(self, name + 'Frames', None)
    if isinstance(frames, list) and len(frames) == int(self.Nt):
        frame = np.asarray(frames[timestep], dtype=np.float32).ravel(order='F')
        offsets = np.concatenate(([0], np.cumsum(self.nProjSubset[timestep], dtype=np.int64)))
        start = int(offsets[subset]) * stride
        stop = int(offsets[subset + 1]) * stride
        return frame[start:stop]
    flat = np.asarray(getattr(self, name), dtype=np.float32).ravel(order='F')
    q = timestep * int(self.subsets) + subset
    start = int(self.nMeas[q]) * stride
    stop = int(self.nMeas[q + 1]) * stride
    return flat[start:stop]


def _full_coordinate_frame(self, name, timestep):
    frames = getattr(self, name + 'Frames', None)
    if isinstance(frames, list) and len(frames) == int(self.Nt):
        return np.asarray(frames[timestep], dtype=np.float32).ravel(order='F')
    return np.asarray(getattr(self, name), dtype=np.float32).ravel(order='F')


def _initialize_coordinate_buffers(self, upload):
    """Create the canonical ``[timestep][subset]`` geometry buffers."""
    self.d_x = [[None] * self.subsets for _ in range(self.Nt)]
    self.d_z = [[None] * self.subsets for _ in range(self.Nt)]
    z_stride = 6 if self.pitch else (3 if self.PET and getattr(self, 'nLayers', 0) > 1 else 2)
    subset_geometry = (self.CT or self.SPECT) and self.listmode == 0
    pet_geometry = self.PET and self.listmode == 0
    for timestep in range(self.Nt):
        if subset_geometry or (self.listmode > 0 and not self.useIndexBasedReconstruction and self.loadTOF):
            for subset in range(self.subsets):
                self.d_x[timestep][subset] = upload(
                    _coordinate_slice(self, 'x', timestep, subset, 6)
                )
        else:
            self.d_x[timestep][0] = upload(_full_coordinate_frame(self, 'x', timestep))

        if subset_geometry or pet_geometry:
            for subset in range(self.subsets):
                self.d_z[timestep][subset] = upload(
                    _coordinate_slice(self, 'z', timestep, subset, z_stride)
                )
        elif self.listmode == 0 or (self.listmode > 0 and self.useIndexBasedReconstruction):
            self.d_z[timestep][0] = upload(_full_coordinate_frame(self, 'z', timestep))
        else:
            self.d_z[timestep][0] = upload(np.zeros(1, dtype=np.float32))


def _initialize_detector_vector_buffers(self, upload, empty=None):
    """Create detector-head-index buffers aligned with each frame/subset geometry slice."""
    self.d_detectorVector = [[empty] * self.subsets for _ in range(self.Nt)]
    if not self.SPECT:
        return
    frames = getattr(self, 'DetectorVectorFrames', None)
    if not isinstance(frames, list) or len(frames) != self.Nt:
        frames = [np.asarray(self.DetectorVector, dtype=np.uint32).reshape(-1)] * self.Nt
    for timestep in range(self.Nt):
        frame = np.asarray(frames[timestep], dtype=np.uint32).reshape(-1)
        offsets = np.concatenate(([0], np.cumsum(self.nProjSubset[timestep], dtype=np.int64)))
        if frame.size != int(offsets[-1]):
            raise ValueError('DetectorVector does not match the reordered projections in a SPECT timeframe.')
        for subset in range(self.subsets):
            self.d_detectorVector[timestep][subset] = upload(
                frame[int(offsets[subset]) : int(offsets[subset + 1])]
            )

def computeGeom5(x, uv, nRowsD, nColsD, dPitchY, pitch):
    """
    Precompute the per-projection geometry used by the branchless distance-driven
    backprojection (projectorType5.cl built with -DGEOM5). Returns 16 floats per
    projection: s (3), d3 (3), normX (3), normY (3), crossP (3), upperPart (1).
    """
    import numpy as np
    nProj = x.size // 6
    xs = np.reshape(x, (nProj, 6)).astype(np.float32, copy=False)
    s = xs[:, 0:3]
    d = xs[:, 3:6]
    indX = np.float32(nRowsD) / np.float32(2.)
    indY = np.float32(nColsD) / np.float32(2.)
    if pitch:
        uvs = np.reshape(uv, (nProj, 6)).astype(np.float32, copy=False)
        apuX = uvs[:, 0:3] * indX
        apuY = uvs[:, 3:6] * indY
    else:
        uvs = np.reshape(uv, (nProj, 2)).astype(np.float32, copy=False)
        apuX = np.zeros((nProj, 3), dtype=np.float32)
        apuX[:, 0] = uvs[:, 0] * indX
        apuX[:, 1] = uvs[:, 1] * indX
        apuY = np.zeros((nProj, 3), dtype=np.float32)
        apuY[:, 2] = indY * np.float32(dPitchY)
    d3 = d - apuX - apuY
    d2 = apuX - apuY
    normX = apuX / np.linalg.norm(apuX, axis=1, keepdims=True)
    normY = apuY / np.linalg.norm(apuY, axis=1, keepdims=True)
    crossP = np.cross(d2, d3 - d)
    upperPart = np.sum(crossP * (s - d), axis=1, keepdims=True)
    geom = np.hstack((s, d3, normX, normY, crossP, upperPart)).astype(np.float32)
    return np.ascontiguousarray(geom).ravel()

def initProjector(self):
    self.CTAttenuation = self.CT_attenuation # TODO: consistent CT_attenuation vs CTAttenuation?
    try:
        import arrayfire as af
    except ModuleNotFoundError:
        if self.useAF:
            print('ArrayFire selected, but not found. Aborting.')
            return
    self.projectorInitialized = True
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
        if self.useMetal:
            if self.useCUDA:
                raise ValueError('Select either Metal/MPS or CUDA, not both.')
            if not torch.backends.mps.is_available():
                raise RuntimeError('PyTorch MPS is not available on this machine.')
            self.useCuPy = False
        else:
            torch.cuda.init()
            self.useCuPy = True
    if self.useTorch and not self.useCUDA and not self.useMetal:
        raise ValueError('PyTorch with the OMEGA OpenCL backend would require host staging. Select CUDA or the native Metal/MPS backend.')
    if self.useCuPy and self.useCUDA:
        import cupy as cp
    elif self.useCuPy and not self.useCUDA:
        print('CuPy can only be used when useCUDA is True. Setting useCUDA to True!')
        self.useCUDA = True
        import cupy as cp
    if self.useCuPy and self.useCUDA:
        def cupyROCm():
            try:
                return bool(cp.cuda.runtime.is_hip)
            except AttributeError:
                pass
            try:
                import io
                import contextlib
                buf = io.StringIO()
                with contextlib.redirect_stdout(buf):
                    cp.show_config()
                lower = buf.getvalue().lower()
                return "rocm" in lower or "hip" in lower
            except Exception:
                return False
    if not self.useCUDA and not self.useMetal:
        import pyopencl as cl
        from pyopencl.version import VERSION
        
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
    
    self.NVOXELS = 8
    self.TH = 100000000000.
    self.TH32 = 100000.
    self.NVOXELS5 = 1
    self.NVOXELSFP = 8
    if np.size(self.weights) > 0:
        self.empty_weight = False
    if self.TOF_bins_used == 0:
        self.TOF_bins_used = 1
    mDataFound = (
        any(np.asarray(frame).size > 0 for frame in self.SinM)
        if isinstance(self.SinM, list)
        else self.SinM.size > 0
    )
    loadCorrections(self)
    parseInputs(self, mDataFound)
    prepassPhase(self)
    # if self.listmode > 0 and self.subsets > 1 and self.subsetType > 0:
    #     if self.useIndexBasedReconstruction:
    #         self.trIndex = self.trIndex[:,self.index]
    #         self.axIndex = self.axIndex[:,self.index]
        # else:
        #     self.x = np.reshape(self.x, [6, -1], order='F')
        #     self.x = self.x[:,self.index]
        #     self.x = self.x.ravel('F')
    if self.useIndexBasedReconstruction and self.listmode > 0:
        self.trIndex = self.trIndex.ravel('F')
        self.axIndex = self.axIndex.ravel('F')

    if self.projector_type in [1, 11, 14, 15, 12, 13, 16]:
        self.FPType = 1
    elif self.projector_type in [2, 21, 22, 23, 24, 25, 26]:
        self.FPType = 2
    elif self.projector_type in [3, 31, 32, 33, 34, 35]:
        self.FPType = 3
    elif self.projector_type in [4, 41, 42, 43, 44, 45]:
        self.FPType = 4
    elif self.projector_type in [5, 51, 52, 53, 54, 55]:
        self.FPType = 5
    elif self.projector_type in [6, 61, 62, 66]:
        self.FPType = 6
    else:
        raise ValueError('Invalid forward projector!')
    if self.projector_type in [1, 11, 21, 31, 41, 51, 61]:
        self.BPType = 1
    elif self.projector_type in [2, 12, 22, 32, 42, 52, 62]:
        self.BPType = 2
    elif self.projector_type in [3, 13, 23, 33, 43, 53]:
        self.BPType = 3
    elif self.projector_type in [4, 14, 24, 34, 44, 54]:
        self.BPType = 4
    elif self.projector_type in [5, 15, 25, 35, 45, 55]:
        self.BPType = 5
    elif self.projector_type in [6, 16, 26, 66]:
        self.BPType = 6
    else:
        raise ValueError('Invalid backprojector!')
    # CuPy does not support the texture API (cupy.cuda.texture) on ROCm/HIP; creating a CUDA
    # array fails at runtime with hipErrorUnknown. Fall back to buffers where the kernels
    # support them, otherwise raise an error.
    if self.useCuPy and self.useCUDA and cupyROCm():
        if self.FPType in [4, 5] or self.BPType == 5:
            raise ValueError('Forward projector types 4 and 5 and backprojector type 5 require texture support, which CuPy does not provide on ROCm/HIP. Use forward projector types 1-3 and/or backprojector types 1-4 instead.')
        if self.useImages:
            print('CuPy does not support textures on ROCm/HIP. Setting useImages to False, buffers will be used instead!')
            self.useImages = False
    # if self.useAF == False and (self.FPType == 5 or self.BPType == 5):
    #     raise ValueError('Branchless distance-driven (projector type 5) can only be used with Arrayfire!')
    if (self.useAF == False and self.useCuPy == False and not self.useMetal) and self.projector_type in (6, 66):
        raise ValueError('Projector type 6 can only be used with Arrayfire (OpenCL), CuPy (CUDA), or PyTorch MPS!')
    if self.projector_type in (16, 26, 61, 62) and not self.useMetal:
        raise ValueError('Hybrid projector types 16, 26, 61, and 62 are supported only by the PyTorch MPS custom-operator path!')
        
    if self.FPType != 6 or self.BPType != 6:
        fPath = os.path.dirname( __file__ )
        if os.path.exists(os.path.join(fPath, '..', 'util', 'usingPyPi.py')):
            headerDir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'opencl')) + "/"
        else:
            headerDir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', '..', '..', 'opencl')) + "/"
        with open(headerDir + 'general_opencl_functions.h', encoding="utf8") as f:
            hlines = f.read()
        linesFP = None
        linesBP = None
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
        frame_count = int(getattr(self, 'Nt', 1))
        globalSize = [[None] * self.subsets for _ in range(frame_count)]
        # self.mSize = [None] * self.subsets
        for timestep in range(frame_count):
            for i in range(self.subsets):
                n_proj = int(self.nProjSubset[timestep, i])
                n_meas = int(self.nMeasSubset[timestep, i])
                if (self.FPType == 5):
                    globalSize[timestep][i] = (self.nRowsD, (self.nColsD + self.NVOXELSFP - 1) // self.NVOXELSFP, n_proj)
                    localSize = (16, 16, 1)
                    erotus = (localSize[0] - (globalSize[timestep][i][0] % localSize[0]), localSize[1] - (globalSize[timestep][i][1] % localSize[1]), 0)
                    globalSize[timestep][i] = (self.nRowsD + erotus[0], (self.nColsD + self.NVOXELSFP - 1) // self.NVOXELSFP + erotus[1], n_proj)
                elif ((self.CT or self.SPECT or self.PET) and self.listmode == 0):
                    globalSize[timestep][i] = (self.nRowsD, self.nColsD, n_proj)
                    localSize = (16, 16, 1)
                    erotus = (localSize[0] - (globalSize[timestep][i][0] % localSize[0]), localSize[1] - (globalSize[timestep][i][1] % localSize[1]), 0)
                    globalSize[timestep][i] = (self.nRowsD + erotus[0], self.nColsD + erotus[1], n_proj)
                else:
                    globalSize[timestep][i] = (n_meas, 1, 1)
                    localSize = (128, 1, 1)
                    erotus = (localSize[0] - (globalSize[timestep][i][0] % localSize[0]), localSize[1] - (globalSize[timestep][i][1] % localSize[1]), 0)
                    globalSize[timestep][i] = (n_meas + erotus[0], 1, 1)
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
            globalSize = [
                [[None] * (self.nMultiVolumes + 1) for _ in range(self.subsets)]
                for _ in range(frame_count)
            ]
            for timestep in range(frame_count):
                for subset in range(self.subsets):
                    for ii in range(self.nMultiVolumes + 1):
                        globalSize[timestep][subset][ii] = self.globalSizeFP[timestep][subset]
            self.localSizeBP = self.localSizeFP + tuple()
        else:
            globalSize = [
                [[None] * (self.nMultiVolumes + 1) for _ in range(self.subsets)]
                for _ in range(frame_count)
            ]
            for timestep in range(frame_count):
                for subset in range(self.subsets):
                    for ii in range(self.nMultiVolumes + 1):
                        if self.BPType == 4:
                            globalSize[timestep][subset][ii] = (self.Nx[ii].item() + self.erotusBP[ii * 2], self.Ny[ii].item() + self.erotusBP[ii * 2 + 1], (self.Nz[ii].item()  + self.NVOXELS - 1) // self.NVOXELS)
                        elif self.BPType == 5:
                            if self.pitch:
                                globalSize[timestep][subset][ii] = (self.Nx[ii].item() + self.erotusBP[ii * 2], self.Ny[ii].item() + self.erotusBP[ii * 2 + 1], self.Nz[ii].item())
                            else:
                                globalSize[timestep][subset][ii] = (self.Nx[ii].item() + self.erotusBP[ii * 2], self.Ny[ii].item() + self.erotusBP[ii * 2 + 1], (self.Nz[ii].item()  + self.NVOXELS5 - 1) // self.NVOXELS5)
                        else:
                            globalSize[timestep][subset][ii] = (self.Nx[ii].item() + self.erotusBP[ii * 2], self.Ny[ii].item() + self.erotusBP[ii * 2 + 1], self.Nz[ii].item())
            self.localSizeBP = localSize + tuple()
        self.globalSizeBP = globalSize.copy()
                            
        
        self.Nxy = self.Nx[0].item() * self.Ny[0].item()
        vendor = ''
        if self.useMetal:
            # The Metal branch compiles FP/BP source once below. The
            # existing bOpt logic remains the single source of specialization.
            bOpt = ('-DMETAL',)
            self.use_64bit_atomics = False
            self.use_32bit_atomics = False
        elif self.useCUDA:
            if self.use_64bit_atomics or self.use_32bit_atomics:
                self.use_64bit_atomics = False
                self.use_32bit_atomics = False
            if cupyROCm():
                bOpt = ('-DHIP','-DCUPY_HIP_FINITE_ELLIPSE_POWER','-DPYTHON',)
            else:
                bOpt = ('-DCUDA','-DPYTHON',)
        else:
            bOpt =('-cl-single-precision-constant -DOPENCL',)
            import pyopencl as cl
            device = self.clctx.get_info(cl.context_info.DEVICES)
            vendor = device[0].get_info(cl.device_info.VENDOR)
            ext = device[0].get_info(cl.device_info.EXTENSIONS)
            if vendor == 'NVIDIA Corporation':
                self.use_64bit_atomics = False
                self.use_32bit_atomics = False
                bOpt += ('-DNVIDIA',)
            elif vendor == 'Advanced Micro Devices, Inc.':
                self.use_64bit_atomics = False
                self.use_32bit_atomics = False
                bOpt += ('-DAMD',)
            elif ext.find('cl_ext_float_atomics') >= 0:
                self.use_64bit_atomics = False
                self.use_32bit_atomics = False
                bOpt += ('-DINTEL',)
        if self.useMAD:
            if self.useMetal:
                bOpt += ('-DUSEMAD',)
            elif self.useCUDA and cupyROCm():
                bOpt += ('-ffast-math','-DUSEMAD',)
            elif self.useCUDA:
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
            elif linesFP is not None:
                linesFP = hlines + linesFP
            if self.BPType in [2, 3]:
                linesBP = hlines + hlines2 + linesBP
            elif linesBP is not None:
                linesBP = hlines + linesBP
        else:
            linesFP = hlines + linesFP if linesFP is not None else None
            linesBP = hlines + linesBP if linesBP is not None else None
        # if self.FPType == 3 or self.BPType == 3:
        #     bOpt += ('-DVOL',)
        if self.useMaskFP:
            bOpt += ('-DMASKFP',)
            if self.maskFPZ > 1:
                bOpt += ('-DMASKFP3D',)
        if self.SPECT and self.maskFPZ == self.nHeads:
            bOpt += ('-DMASKFPBYDETECTOR',)
        if self.normalization_correction and self.SPECT and self.normZ == self.nHeads:
            bOpt += ('-DNORMBYDETECTOR',)
        if self.useMaskBP:
            bOpt += ('-DMASKBP',)
            if self.maskBPZ > 1:
                bOpt += ('-DMASKBP3D',)
        if self.useTotLength:# and not self.SPECT:
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
        if self.listmode > 0 and ~self.useIndexBasedReconstruction and not vendor == 'NVIDIA Corporation':
            bOpt += ('-DUSEGLOBAL',)
        elif not self.useMetal:
            bOpt += ('-DUSEGLOBAL',)
        if (((self.FPType == 1 or self.BPType == 1 or self.FPType == 4 or self.BPType == 4) and self.n_rays_transaxial * self.n_rays_axial > 1) or self.SPECT):
            bOpt += ('-DN_RAYS=' + str(self.n_rays_transaxial * self.n_rays_axial),)
            bOpt += ('-DN_RAYS2D=' + str(self.n_rays_transaxial),)
            bOpt += ('-DN_RAYS3D=' + str(self.n_rays_axial),)
        if self.pitch:
            bOpt += ('-DPITCH',)
        if (((self.subsets > 1 and (self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7))) and not self.CT and not self.SPECT and not self.PET and self.listmode == 0):
            bOpt += ('-DSUBSETS',)
        if self.subsets > 1 and self.listmode == 0:
            bOpt += ('-DSTYPE=' + str(self.subsetType),'-DNSUBSETS=' + str(self.subsets),)
        
        # hipRTC force-includes hiprtc_runtime.h, which uses FP as a type name, so a command-line
        # -DFP breaks that header. ROCm builds pass -DOMEGA_FP instead; general_opencl_functions.h
        # maps it back to FP after the hipRTC prelude.
        if self.useCUDA and cupyROCm():
            bOptFP = bOpt + ('-DOMEGA_FP',)
        else:
            bOptFP = bOpt + ('-DFP',)
        if self.localSizeFP[1] > 1:
            bOptFP += ('-DLOCAL_SIZE=' + str(self.localSizeFP[0]),'-DLOCAL_SIZE2=' + str(self.localSizeFP[1]),)
        else:
            bOptFP += ('-DLOCAL_SIZE=' + str(self.localSizeFP[0]),'-DLOCAL_SIZE2=' + str(1),)
        if self.FPType in [1, 2, 3]:
            bOptFP += ('-DSIDDON',)
            bOptFP += ('-DATOMICF',)
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
            if self.CT:
                bOptBP += ('-DBP4','-DNVOXELS=' + str(self.NVOXELS),)
            else:
                bOptBP += ('-DPTYPE4','-DNVOXELS=' + str(self.NVOXELS),)
                bOptBP += ('-DATOMICF',)
                if self.use_64bit_atomics:
                    bOptBP += ('-DATOMIC','-DCAST=long','-DTH=' + str(self.TH),)
                elif self.use_32bit_atomics:
                    bOptBP += (' -DATOMIC32',' -DCAST=int',' -DTH=' + str(self.TH32),)
                else:
                    bOptBP += ('-DCAST=float',)
        elif self.BPType == 5:
            bOptBP += ('-DPROJ5','-DNVOXELS5=' + str(self.NVOXELS5),)
            # Use the per-projection geometry precomputed with computeGeom5 (stored in self.d_geom5)
            if self.CT and self.listmode == 0:
                bOptBP += ('-DGEOM5',)
            if self.meanBP:
                bOptBP += ('-DMEANDISTANCEBP',)
    else:
        headerDir = None
        linesFP = None
        linesBP = None
        bOptFP = ()
        bOptBP = ()
        if self.useCUDA:
            if self.useTorch:
                #self.gFilter = np.ascontiguousarray(self.gFilter)
                #self.gFilter = np.transpose(self.gFilter, (1, 0, 2))
                if isinstance(self.gFilter, (list, tuple)) and len(self.gFilter):
                    self.d_gFilter = [torch.tensor(value, device='cuda') for value in self.gFilter]
                else:
                    self.d_gFilter = torch.tensor(self.gFilter, device='cuda')
                #self.d_gFilter = self.d_gFilter.permute(2, 0, 1).unsqueeze(1)
        elif not self.useMetal:
            if isinstance(self.gFilter, (list, tuple)) and len(self.gFilter):
                self.d_gFilter = [af.interop.np_to_af_array(value) for value in self.gFilter]
            else:
                self.d_gFilter = af.interop.np_to_af_array(self.gFilter)
        self.uu = 0
    
    if self.useMetal:
        from omegatomo.projector.mps_backend import init_mps_projector
        init_mps_projector(
            self,
            source_root=headerDir,
            source_fp=linesFP,
            source_bp=linesBP,
            options_fp=bOptFP,
            options_bp=bOptBP,
        )
        return

    if self.useCUDA:
        self.no_norm = 1
        self.mSize = self.nRowsD * self.nColsD * self.nProjections
        self.d_d = [None] * (self.nMultiVolumes + 1)
        self.d_b = [None] * (self.nMultiVolumes + 1)
        self.d_bmax = [None] * (self.nMultiVolumes + 1)
        self.d_Nxyz = [None] * (self.nMultiVolumes + 1)
        self.dSize = [None] * (self.nMultiVolumes + 1)
        self.d_Scale = [None] * (self.nMultiVolumes + 1)
        self.d_Scale4 = [None] * (self.nMultiVolumes + 1)
        if self.FPType != 6 or self.BPType != 6:
            if self.useCuPy:
                # if self.FPType == 5:
                #     raise ValueError('Not yet supported')
                self.d_Sens = cp.empty(shape=(1,1), dtype=cp.float32)
                _initialize_coordinate_buffers(self, lambda value: cp.asarray(value))
                # Precomputed per-projection geometry for the BDD backprojection (see -DGEOM5 in projectorType5.cl)
                if self.BPType == 5 and self.CT and self.listmode == 0:
                    self.d_geom5 = [[None] * self.subsets for _ in range(self.Nt)]
                    kerroin = 6 if self.pitch else 2
                    for timestep in range(self.Nt):
                        for subset in range(self.subsets):
                            geom = computeGeom5(
                                _coordinate_slice(self, 'x', timestep, subset, 6),
                                _coordinate_slice(self, 'z', timestep, subset, kerroin),
                                self.nRowsD,
                                self.nColsD,
                                self.dPitchY,
                                self.pitch,
                            )
                            self.d_geom5[timestep][subset] = cp.asarray(geom)
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
                        if self.BPType == 4 and not self.CT:
                            tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                    filterMode=cp.cuda.runtime.cudaFilterModeLinear, normalizedCoords=1)
                        else:
                            tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                    filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=0)
                        self.d_atten = cp.cuda.texture.TextureObject(res, tdes)
                if self.useMaskFP:
                    if not self.useImages:
                        self.d_maskFP = cp.asarray(self.maskFP)
                    else:
                        chl = cp.cuda.texture.ChannelFormatDescriptor(8,0,0,0, cp.cuda.runtime.cudaChannelFormatKindUnsigned)
                        self.maskFP = self.maskFP.ravel('F')
                        if self.SPECT and self.maskFPZ > 1 and self.maskFPZ == self.nHeads:
                            array = cp.cuda.texture.CUDAarray(chl, self.nRowsD, self.nColsD, self.nHeads)
                            array.copy_from(self.maskFP.reshape((self.nHeads, self.nColsD, self.nRowsD)))
                            res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                            tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp),
                                                                    filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=0)
                            self.d_maskFP = cp.cuda.texture.TextureObject(res, tdes)
                        elif self.maskFPZ > 1:
                            self.d_maskFP = [None] * self.subsets
                            for i in range(self.subsets):
                                array = cp.cuda.texture.CUDAarray(chl, self.nRowsD, self.nColsD, self.nMeas[i])
                                self.maskFP = self.maskFP.reshape((self.nMeas[i], self.nColsD, self.nRowsD))
                                array.copy_from(self.maskFP[self.nMeas[i] : self.nMeas[i + 1],:,:])
                                # array.copy_from(self.maskFP[:,:,self.nMeas[i] : self.nMeas[i + 1]].reshape((self.nMeas[i], self.nColsD, self.nRowsD)))
                                res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                                tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                        filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=0)
                                self.d_maskFP[i] = cp.cuda.texture.TextureObject(res, tdes)
                        else:
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
                        self.maskBP = self.maskBP.ravel('F')
                        if self.maskBPZ > 1:
                            array = cp.cuda.texture.CUDAarray(chl, self.Nx[0].item(), self.Ny[0].item(), self.maskBPZ)
                            array.copy_from(self.maskBP.reshape((self.maskBPZ, self.Ny[0].item(), self.Nx[0].item())))
                        else:
                            array = cp.cuda.texture.CUDAarray(chl, self.Nx[0].item(), self.Ny[0].item())
                            array.copy_from(self.maskBP.reshape((self.Ny[0].item(), self.Nx[0].item())))
                        res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                        if self.BPType == 4 and not self.CT:
                            tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                    filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=1)
                        else:
                            tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp), 
                                                                    filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=0)
                        self.d_maskBP = cp.cuda.texture.TextureObject(res, tdes)
                if self.TOF:
                    self.d_TOFCenter = cp.asarray(self.TOFCenter)
                _initialize_detector_vector_buffers(self, lambda value: cp.asarray(value))
                if self.SPECT:
                    self.d_rayShiftsDetector = cp.asarray(self.rayShiftsDetector)
                    self.d_rayShiftsSource = cp.asarray(self.rayShiftsSource)
                if (self.BPType == 2 or self.BPType == 3 or self.FPType == 2 or self.FPType == 3):
                    self.d_V = cp.asarray(self.V)
                if (self.normalization_correction):
                    self.d_norm = [None] * self.subsets
                    for i in range(self.subsets):
                        if self.SPECT and self.normZ == self.nHeads:
                            self.d_norm[i] = cp.asarray(self.normalization)
                        else:
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
                        self.d_T[i] = cp.asarray(self.OffsetLimit[self.nMeas[i].item() : self.nMeas[i + 1].item()])
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

                # ``inf`` is the public box-support value.  hipRTC receives
                # a finite sentinel instead, because ``isinf`` is unreliable
                # for a runtime scalar when HIP fast-math is enabled.
                ellipse_power_kernel = self.ellipsePower
                if cupyROCm() and np.isinf(ellipse_power_kernel):
                    ellipse_power_kernel = np.finfo(np.float32).max
                    
                if self.FPType in [1, 2, 3]:
                    self.kIndF = (cp.float32(self.global_factor), cp.float32(self.epps), cp.uint32(self.nRowsD), cp.uint32(self.det_per_ring), cp.float32(self.sigma_x),)
                    if self.SPECT:
                        self.kIndF += (self.d_rayShiftsDetector, self.d_rayShiftsSource, cp.float32(self.coneOfResponseStdCoeffA), cp.float32(self.coneOfResponseStdCoeffB), cp.float32(self.coneOfResponseStdCoeffC), cp.float32(self.ellipseCenterX), cp.float32(self.ellipseCenterY), cp.float32(self.ellipseCenterZ), cp.float32(self.ellipseRadiusX), cp.float32(self.ellipseRadiusY), cp.float32(self.ellipseRadiusZ), cp.float32(ellipse_power_kernel),)
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
                        self.kIndB += (self.d_rayShiftsDetector, self.d_rayShiftsSource, cp.float32(self.coneOfResponseStdCoeffA), cp.float32(self.coneOfResponseStdCoeffB), cp.float32(self.coneOfResponseStdCoeffC), cp.float32(self.ellipseCenterX), cp.float32(self.ellipseCenterY), cp.float32(self.ellipseCenterZ), cp.float32(self.ellipseRadiusX), cp.float32(self.ellipseRadiusY), cp.float32(self.ellipseRadiusZ), cp.float32(ellipse_power_kernel),)
                    self.kIndB += (cp.float32(self.dPitchX),cp.float32(self.dPitchY),)
                    if self.BPType in [2, 3]:
                        if self.BPType == 2:
                            self.kIndB  += (cp.float32(self.tube_width_z),)
                        else:
                            self.kIndB  += (cp.float32(self.tube_radius),)
                        self.kIndB += (cp.float32(self.bmin),)
                        self.kIndB += (cp.float32(self.bmax),)
                        self.kIndB += (cp.float32(self.Vmax),)
                    # if self.useMaskFP:
                        # self.kIndB += (self.d_maskFP,)
                # if self.useMaskBP:
                    # self.kIndB += (self.d_maskBP,)
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
                raise ValueError('Unsupported selection. Note that PyCUDA is no longer supported!')
    else:
        if self.FPType != 6 or self.BPType != 6:
            
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
            _initialize_coordinate_buffers(
                self,
                lambda value: cl.array.to_device(self.queue, value),
            )
            # Precomputed per-projection geometry for the BDD backprojection (see -DGEOM5 in projectorType5.cl)
            if self.BPType == 5 and self.CT and self.listmode == 0:
                self.d_geom5 = [[None] * self.subsets for _ in range(self.Nt)]
                kerroin = 6 if self.pitch else 2
                for timestep in range(self.Nt):
                    for subset in range(self.subsets):
                        geom = computeGeom5(
                            _coordinate_slice(self, 'x', timestep, subset, 6),
                            _coordinate_slice(self, 'z', timestep, subset, kerroin),
                            self.nRowsD,
                            self.nColsD,
                            self.dPitchY,
                            self.pitch,
                        )
                        self.d_geom5[timestep][subset] = cl.array.to_device(self.queue, geom)
            if (self.attenuation_correction and not self.CTAttenuation):
                self.d_atten = [None] * self.subsets
                for i in range(self.subsets):
                    self.d_atten[i] = cl.array.to_device(self.queue, self.vaimennus[self.nTotMeas[i].item() : self.nTotMeas[i + 1].item()])
            elif (self.attenuation_correction and self.CTAttenuation):
                if self.useImages:
                    imformat = cl.ImageFormat(cl.channel_order.A, cl.channel_type.FLOAT)
                    if VERSION[0] > 2024 or (VERSION[0] == 2024 and VERSION[1] > 2):
                        self.d_atten = cl.create_image(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, imformat, hostbuf=self.vaimennus, shape=(self.Nx[0].item(), self.Ny[0].item(), self.Nz[0].item()))
                    else:
                        self.d_atten = cl.Image(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, imformat, hostbuf=self.vaimennus, shape=(self.Nx[0].item(), self.Ny[0].item(), self.Nz[0].item()))
                else:
                    self.d_atten = cl.array.to_device(self.queue, self.vaimennus)
                # self.d_atten = cl.image_from_array(self.clctx, np.reshape(self.vaimennus, (self.Nx[0].item(), self.Ny[0].item(), self.Nz[0].item()), order='F'))
            _initialize_detector_vector_buffers(self, lambda value: cl.array.to_device(self.queue, value))
            if self.SPECT:
                self.d_rayShiftsDetector = cl.array.to_device(self.queue, self.rayShiftsDetector)
                self.d_rayShiftsSource = cl.array.to_device(self.queue, self.rayShiftsSource)
            if self.useMaskFP:
                imformat = cl.ImageFormat(cl.channel_order.A, cl.channel_type.UNSIGNED_INT8)
                if self.SPECT and self.maskFPZ > 1 and self.maskFPZ == self.nHeads:
                    if VERSION[0] > 2024 or (VERSION[0] == 2024 and VERSION[1] > 2):
                        self.d_maskFP = cl.create_image(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, imformat, hostbuf=self.maskFP, shape=(self.nRowsD, self.nColsD, self.nHeads))
                    else:
                        self.d_maskFP = cl.Image(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, imformat, hostbuf=self.maskFP, shape=(self.nRowsD, self.nColsD, self.nHeads))
                elif self.maskFPZ > 1:
                    self.d_maskFP = [None] * self.subsets
                    for i in range(self.subsets):
                        if VERSION[0] > 2024 or (VERSION[0] == 2024 and VERSION[1] > 2):
                            self.d_maskFP[i] = cl.create_image(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, imformat, hostbuf=self.maskFP, shape=(self.nRowsD, self.nColsD, self.nMeas[i]))
                        else:
                            self.d_maskFP[i] = cl.Image(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, imformat, hostbuf=self.maskFP, shape=(self.nRowsD, self.nColsD, self.nMeas[i]))
                else:
                    if VERSION[0] > 2024 or (VERSION[0] == 2024 and VERSION[1] > 2):
                        self.d_maskFP = cl.create_image(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, imformat, hostbuf=self.maskFP, shape=(self.nRowsD, self.nColsD))
                    else:
                        self.d_maskFP = cl.Image(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, imformat, hostbuf=self.maskFP, shape=(self.nRowsD, self.nColsD))
                # self.d_maskFP = cl.image_from_array(self.clctx, np.ascontiguousarray(self.maskFP))
            if self.useMaskBP:
                imformat = cl.ImageFormat(cl.channel_order.A, cl.channel_type.UNSIGNED_INT8)
                if self.maskBPZ > 1:
                    if VERSION[0] > 2024 or (VERSION[0] == 2024 and VERSION[1] > 2):
                        self.d_maskBP = cl.create_image(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, imformat, hostbuf=self.maskBP, shape=(self.Nx[0].item(), self.Ny[0].item(), self.maskBPZ))
                    else:
                        self.d_maskBP = cl.Image(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, imformat, hostbuf=self.maskBP, shape=(self.Nx[0].item(), self.Ny[0].item(), self.maskBPZ))
                else:
                    if VERSION[0] > 2024 or (VERSION[0] == 2024 and VERSION[1] > 2):
                        self.d_maskBP = cl.create_image(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, imformat, hostbuf=self.maskBP, shape=(self.Nx[0].item(), self.Ny[0].item()))
                    else:
                        self.d_maskBP = cl.Image(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, imformat, hostbuf=self.maskBP, shape=(self.Nx[0].item(), self.Ny[0].item()))
                # self.d_maskBP = cl.image_from_array(self.clctx, np.ascontiguousarray(self.maskBP))
            if self.TOF:
                self.d_TOFCenter = cl.array.to_device(self.queue, self.TOFCenter)
            if (self.BPType == 2 or self.BPType == 3 or self.FPType == 2 or self.FPType == 3):
                self.d_V = cl.array.to_device(self.queue, self.V)
            if (self.normalization_correction):
                self.d_norm = [None] * self.subsets
                for i in range(self.subsets):
                    if self.SPECT and self.normZ == self.nHeads:
                        self.d_norm[i] = cl.array.to_device(self.queue, self.normalization)
                    else:
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
                    self.d_T[i] = cl.array.to_device(self.queue, self.OffsetLimit[self.nMeas[i].item() : self.nMeas[i + 1].item()])
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
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.ellipseCenterX))
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.ellipseCenterY))
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.ellipseCenterZ))
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.ellipseRadiusX))
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.ellipseRadiusY))
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.ellipseRadiusZ))
                    self.kIndF += 1
                    self.knlF.set_arg(self.kIndF, (cl.cltypes.float)(self.ellipsePower))
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
            # if self.useMaskFP:
            #     self.knlF.set_arg(self.kIndF, self.d_maskFP)
            #     self.kIndF += 1
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
                if self.useImages:
                    self.knlF.set_arg(self.kIndF, self.d_atten)
                else:
                    self.knlF.set_arg(self.kIndF, self.d_atten.data)
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
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.ellipseCenterX))
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.ellipseCenterY))
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.ellipseCenterZ))
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.ellipseRadiusX))
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.ellipseRadiusY))
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.ellipseRadiusZ))
                    self.kIndB += 1
                    self.knlB.set_arg(self.kIndB, (cl.cltypes.float)(self.ellipsePower))
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
                if self.useImages:
                    self.knlB.set_arg(self.kIndB, self.d_atten)
                else:
                    self.knlB.set_arg(self.kIndB, self.d_atten.data)
                self.kIndB += 1
            # if self.BPType in [1, 2, 3] and self.useMaskFP:
            #     self.knlB.set_arg(self.kIndB, self.d_maskFP)
            #     self.kIndB += 1
            # if self.BPType in [1, 2, 3] and self.useMaskBP:
            #     self.knlB.set_arg(self.kIndB, self.d_maskBP)
            #     self.kIndB += 1
