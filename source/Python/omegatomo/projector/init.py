# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 13:17:22 2025
"""
import numpy as np


def _kernel_ellipse_power(value):
    # Kernels test for box support (ellipsePower = inf) with a finite threshold,
    # since isinf() is unreliable under fast-math
    value = float(value)
    return float(np.finfo(np.float32).max) if not np.isfinite(value) else value


def _build_kIndF(self):
    """The constant FP kernel-argument prefix (everything set once at init
    time, before any per-subset geometry/output arguments), shared by the
    OpenCL/AF ``set_arg`` chain and the CuPy argument tuple -- both used to
    hand-duplicate this exact sequence. See projfunctions._KernelArgs."""
    from omegatomo.projector.projfunctions import _KernelArgs
    ellipse_power_kernel = _kernel_ellipse_power(self.ellipsePower)
    a = _KernelArgs(self)
    if self.FPType in (1, 2, 3):
        a.f32(self.global_factor).f32(self.epps).u32(self.nRowsD).u32(self.det_per_ring).f32(self.sigma_x)
        if self.SPECT:
            a.buf(self.d_rayShiftsDetector).buf(self.d_rayShiftsSource)
            a.f32(self.coneOfResponseStdCoeffA).f32(self.coneOfResponseStdCoeffB).f32(self.coneOfResponseStdCoeffC)
            a.vec3f(self.ellipseCenterX, self.ellipseCenterY, self.ellipseCenterZ)
            a.vec3f(self.ellipseRadiusX, self.ellipseRadiusY, self.ellipseRadiusZ)
            a.f32(ellipse_power_kernel)
        a.vec2f(self.dPitchX, self.dPitchY)
        if self.FPType in (2, 3):
            if self.FPType == 2:
                a.f32(self.tube_width_z)
            else:
                a.f32(self.tube_radius)
            a.f32(self.bmin).f32(self.bmax).f32(self.Vmax)
    elif self.FPType == 4:
        a.u32(self.nRowsD).u32(self.nColsD).vec2f(self.dPitchX, self.dPitchY).f32(self.dL).f32(self.global_factor)
    elif self.FPType == 5:
        a.u32(self.nRowsD).u32(self.nColsD).vec2f(self.dPitchX, self.dPitchY)
    if self.FPType in (1, 2, 3):
        if self.TOF:
            a.buf(self.d_TOFCenter)
        if self.FPType in (2, 3):
            a.buf(self.d_V)
        a.u32(self.nColsD)
    # NOTE (pre-existing, preserved): the CuPy source this replaces never
    # had the "not self.CT" guard the OpenCL/AF source had here -- kept
    # per-backend exactly as before rather than silently unified, since TOF
    # is not currently exercised together with CT by any known config and
    # this is not the place to change behaviour.
    if self.FPType == 4 and self.TOF and (self.useCUDA or not self.CT):
        a.buf(self.d_TOFCenter)
        a.f32(self.sigma_x)
    if self.attenuation_correction and self.CTAttenuation and self.FPType in (1, 2, 3, 4):
        if self.useImages or self.useCUDA:
            a.img(self.d_atten)
        else:
            a.buf(self.d_atten)
    return a


def _build_kIndB(self):
    """The constant BP kernel-argument prefix -- see _build_kIndF."""
    from omegatomo.projector.projfunctions import _KernelArgs
    ellipse_power_kernel = _kernel_ellipse_power(self.ellipsePower)
    a = _KernelArgs(self)
    if self.BPType in (4, 5):
        a.u32(self.nRowsD).u32(self.nColsD).vec2f(self.dPitchX, self.dPitchY)
        if self.BPType == 4 and not self.CT:
            a.f32(self.dL).f32(self.global_factor)
    elif self.BPType in (1, 2, 3):
        a.f32(self.global_factor).f32(self.epps).u32(self.nRowsD).u32(self.det_per_ring).f32(self.sigma_x)
        if self.SPECT:
            a.buf(self.d_rayShiftsDetector).buf(self.d_rayShiftsSource)
            a.f32(self.coneOfResponseStdCoeffA).f32(self.coneOfResponseStdCoeffB).f32(self.coneOfResponseStdCoeffC)
            a.vec3f(self.ellipseCenterX, self.ellipseCenterY, self.ellipseCenterZ)
            a.vec3f(self.ellipseRadiusX, self.ellipseRadiusY, self.ellipseRadiusZ)
            a.f32(ellipse_power_kernel)
        a.vec2f(self.dPitchX, self.dPitchY)
        if self.BPType in (2, 3):
            if self.BPType == 2:
                a.f32(self.tube_width_z)
            else:
                a.f32(self.tube_radius)
            a.f32(self.bmin).f32(self.bmax).f32(self.Vmax)
    if self.BPType in (1, 2, 3):
        if self.TOF:
            a.buf(self.d_TOFCenter)
        if self.BPType in (2, 3):
            a.buf(self.d_V)
        a.u32(self.nColsD)
    if self.BPType == 4 and not self.CT and self.TOF:
        a.buf(self.d_TOFCenter)
        a.f32(self.sigma_x)
    if self.attenuation_correction and self.CTAttenuation and self.BPType in (1, 2, 3, 4) and not self.CT:
        if self.useImages or self.useCUDA:
            a.img(self.d_atten)
        else:
            a.buf(self.d_atten)
    return a


def _coordinate_slice(self, name, timestep, subset, stride):
    """Return one frame/subset coordinate slice in kernel order."""
    frames = getattr(self, name + 'Frames', None)
    if isinstance(frames, list) and len(frames) == int(self.Nt):
        frame = np.asarray(frames[timestep], dtype=np.float32).ravel(order='F')
        offsets = np.concatenate(([0], np.cumsum(self.nProjSubset[timestep], dtype=np.int64)))
        start = int(offsets[subset]) * stride
        stop = int(offsets[subset + 1]) * stride
        return frame[start:stop]
    # CT stores one projection per row; SPECT stores one per column.
    order = 'C' if self.CT and self.listmode == 0 else 'F'
    flat = np.asarray(getattr(self, name), dtype=np.float32).ravel(order=order)
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

        if self.SPECT and self.listmode > 0 and not self.useIndexBasedReconstruction:
            z_values = np.asarray(self.z, dtype=np.float32).ravel(order='F')
            for subset in range(self.subsets):
                index = timestep * self.subsets + subset
                start = int(self.nMeas[index]) * 5
                stop = int(self.nMeas[index + 1]) * 5
                self.d_z[timestep][subset] = upload(z_values[start:stop])
        elif subset_geometry or pet_geometry:
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
    if not self.SPECT or self.listmode > 0:
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

# ---------------------------------------------------------------------------
# FPType/BPType table
#
# options.projector_type is either a single digit N (FP=BP=N) or a two-digit
# number "FB" where F is the forward-projector type and B is the
# backprojector type -- but not every F/B combination exists (e.g. there is
# no FP3/BP6), so this is an explicit table, not a pure divmod. It is a
# direct transcription of the two membership-list if/elif chains this
# replaces; verify_fptype_bptype.py checks the two are identical (values
# AND the two raises) over range(0, 70).
# ---------------------------------------------------------------------------
_FPTYPE_GROUPS = {
    1: (1, 11, 12, 13, 14, 15, 16),
    2: (2, 21, 22, 23, 24, 25, 26),
    3: (3, 31, 32, 33, 34, 35),
    4: (4, 41, 42, 43, 44, 45),
    5: (5, 51, 52, 53, 54, 55),
    6: (6, 61, 62, 66),
}
_BPTYPE_GROUPS = {
    1: (1, 11, 21, 31, 41, 51, 61),
    2: (2, 12, 22, 32, 42, 52, 62),
    3: (3, 13, 23, 33, 43, 53),
    4: (4, 14, 24, 34, 44, 54),
    5: (5, 15, 25, 35, 45, 55),
    6: (6, 16, 26, 66),
}
FPTYPE_TABLE = {pt: fp for fp, members in _FPTYPE_GROUPS.items() for pt in members}
BPTYPE_TABLE = {pt: bp for bp, members in _BPTYPE_GROUPS.items() for pt in members}


def _read(directory, name, encoding='utf8'):
    """Read one OpenCL/CUDA/HIP kernel header/source file out of
    `directory` (as returned by omegatomo.util.paths.opencl_header_dir(),
    already ending in '/'). `encoding=None` reproduces the one call site
    (opencl_functions_orth3D.h) that historically opened without an
    explicit encoding."""
    with open(directory + name, encoding=encoding) as f:
        return f.read()


def _cl_image(clctx, flags, imformat, shape, hostbuf):
    """One READ_ONLY/COPY_HOST_PTR OpenCL image, hiding the pyopencl
    VERSION branch (cl.create_image() was added as the replacement for
    the old cl.Image() constructor call convention in newer pyopencl)."""
    import pyopencl as cl
    from pyopencl.version import VERSION
    if VERSION[0] > 2024 or (VERSION[0] == 2024 and VERSION[1] > 2):
        return cl.create_image(clctx, flags, imformat, hostbuf=hostbuf, shape=shape)
    else:
        return cl.Image(clctx, flags, imformat, hostbuf=hostbuf, shape=shape)


def _cupy_texture(data_c_order, shape_xyz, *, linear, normalized, uint8=False, address='clamp'):
    """One CuPy CUDA texture (ChannelFormatDescriptor + CUDAarray +
    ResourceDescriptor + TextureDescriptor + TextureObject), preserving
    each call site's own filterMode/normalizedCoords/channel-format
    choice via the linear/normalized/uint8 flags -- these are NOT the
    same at every site (e.g. the attenuation texture's BPType==4-and-not-CT
    branch uses filterMode=Linear while the mask textures' equivalent
    branch uses filterMode=Point), so callers must pass their own exact
    flags rather than relying on a default."""
    import cupy as cp
    if uint8:
        chl = cp.cuda.texture.ChannelFormatDescriptor(8, 0, 0, 0, cp.cuda.runtime.cudaChannelFormatKindUnsigned)
    else:
        chl = cp.cuda.texture.ChannelFormatDescriptor(32, 0, 0, 0, cp.cuda.runtime.cudaChannelFormatKindFloat)
    array = cp.cuda.texture.CUDAarray(chl, *shape_xyz)
    array.copy_from(data_c_order)
    res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
    addr = (cp.cuda.runtime.cudaAddressModeClamp if address == 'clamp'
            else cp.cuda.runtime.cudaAddressModeBorder)
    filt = cp.cuda.runtime.cudaFilterModeLinear if linear else cp.cuda.runtime.cudaFilterModePoint
    tdes = cp.cuda.texture.TextureDescriptor(addressModes=(addr, addr, addr),
                                              filterMode=filt, normalizedCoords=(1 if normalized else 0))
    return cp.cuda.texture.TextureObject(res, tdes)


def initProjector(self):
    if self.useMetal and not self.useTorch:
        raise ValueError('The Metal/MPS projector requires useTorch=True.')
    # CTAttenuation (the internal mirror of the user-facing CT_attenuation option) is derived once,
    # in addProjector() (proj.py), not here -- see the comment there. addProjector() always runs
    # before initProjector(), so self.CTAttenuation is already set by this point.
    if self.useAF:
        try:
            import arrayfire as af
        except (ImportError, OSError, RuntimeError):
            print('ArrayFire selected, but not found. Aborting.')
            return
    self.projectorInitialized = True
    import numpy as np
    from omegatomo.reconstruction.prepass import prepassPhase
    from omegatomo.reconstruction.prepass import parseInputs
    from omegatomo.reconstruction.prepass import loadCorrections
    if self.useAF:
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

    if self.projector_type not in FPTYPE_TABLE:
        raise ValueError('Invalid forward projector!')
    self.FPType = FPTYPE_TABLE[self.projector_type]
    if self.projector_type not in BPTYPE_TABLE:
        raise ValueError('Invalid backprojector!')
    self.BPType = BPTYPE_TABLE[self.projector_type]
    # CuPy does not support the texture API (cupy.cuda.texture) on ROCm/HIP; creating a CUDA
    # array fails at runtime with hipErrorUnknown. Fall back to buffers where the kernels
    # support them, otherwise raise an error.
    if self.useCuPy and self.useCUDA and cupyROCm():
        if self.FPType in [4, 5] or self.BPType == 5:
            raise ValueError('Forward projector types 4 and 5 and backprojector type 5 require texture support, which CuPy does not provide on ROCm/HIP. Use forward projector types 1-3 and/or backprojector types 1-4 instead.')
        if self.useImages:
            print('CuPy does not support textures on ROCm/HIP. Setting useImages to False, buffers will be used instead!')
            self.useImages = False
    # ProjectorClass.h:979-984 forces -DUSEIMAGES (silently, no message) whenever FPType is 4 or
    # 5 or BPType is 5 -- those kernels have no non-texture code path -- on every backend. The
    # CuPy argument builders (_build_kIndF/_build_kIndB, `self.useImages or self.useCUDA`) also
    # already always hand CUDA a texture for attenuation regardless of useImages. Force
    # useImages the same way here so every useImages-gated choice (the -DUSEIMAGES compile flag
    # below, and every image-vs-buffer argument builder) stays in agreement -- otherwise a
    # user-requested useImages=False could compile a kernel that expects textures while Python
    # still hands it buffers (or vice versa). Skip CuPy-on-ROCm, which cannot use textures at
    # all and was already handled (with its own error/fallback) just above.
    rocm_no_textures = self.useCuPy and self.useCUDA and cupyROCm()
    if not rocm_no_textures and (self.FPType in (4, 5) or self.BPType == 5 or self.useCUDA):
        self.useImages = True
    # if self.useAF == False and (self.FPType == 5 or self.BPType == 5):
    #     raise ValueError('Branchless distance-driven (projector type 5) can only be used with Arrayfire!')
    if (self.useAF == False and self.useCuPy == False and not self.useMetal) and self.projector_type in (6, 66):
        raise ValueError('Projector type 6 can only be used with Arrayfire (OpenCL), CuPy (CUDA), or PyTorch MPS!')
    if self.projector_type in (6, 66) and self.useCUDA and not self.useTorch:
        raise ValueError('Projector type 6 on CUDA requires useTorch=True (PyTorch tensors)!')
    if self.projector_type in (16, 26, 61, 62) and not self.useMetal:
        raise ValueError('Hybrid projector types 16, 26, 61, and 62 are supported only by the PyTorch MPS custom-operator path!')
        
    if self.FPType != 6 or self.BPType != 6:
        from omegatomo.util.paths import opencl_header_dir
        headerDir = opencl_header_dir()
        hlines = _read(headerDir, 'general_opencl_functions.h')
        linesFP = None
        linesBP = None
        if self.FPType in [1, 2, 3]:
            linesFP = _read(headerDir, 'projectorType123.cl')
        elif self.FPType in [4]:
            linesFP = _read(headerDir, 'projectorType4.cl')
        elif self.FPType in [5]:
            linesFP = _read(headerDir, 'projectorType5.cl')
        if self.BPType in [1, 2, 3]:
            linesBP = _read(headerDir, 'projectorType123.cl')
        elif self.BPType in [4]:
            linesBP = _read(headerDir, 'projectorType4.cl')
        elif self.BPType in [5]:
            linesBP = _read(headerDir, 'projectorType5.cl')
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
                bOpt = ('-DHIP','-DPYTHON',)
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
            hlines2 = _read(headerDir, 'opencl_functions_orth3D.h', encoding=None)
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
        if self.listmode > 0 or not self.useMetal:
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

    if (self.useMAD and self.useCUDA and self.useCuPy and cupyROCm()
            and (self.BPType in (1, 2, 3) or (self.BPType == 4 and not self.CT))):
        bOptBP += ("-munsafe-fp-atomics",)

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
                # Backend upload adapter (Stage 1: plain host->device array
                # upload only; image-vs-texture resources go through
                # _cl_image/_cupy_texture above instead, since the CuPy and
                # OpenCL image APIs differ too much to share one call).
                upload = cp.asarray
                self.d_Sens = cp.empty(shape=(1,1), dtype=cp.float32)
                _initialize_coordinate_buffers(self, upload)
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
                            self.d_geom5[timestep][subset] = upload(geom)
                if (self.attenuation_correction and not self.CTAttenuation):
                    # Measurement-domain attenuation is per-timestep in C++/MATLAB (dynamic data
                    # can have frame-dependent correction factors); build [Nt][subsets], mirroring
                    # mps_backend.py's d_attenuation, instead of a flat [subsets] list that only
                    # ever captured frame 0 via nTotMeas[0:subsets].
                    self.d_atten = [[None] * self.subsets for _ in range(self.Nt)]
                    for timestep in range(self.Nt):
                        for i in range(self.subsets):
                            index = timestep * self.subsets + i
                            self.d_atten[timestep][i] = upload(self.vaimennus[self.nTotMeas[index].item() : self.nTotMeas[index + 1].item()])
                elif (self.attenuation_correction and self.CTAttenuation):
                    if not self.useImages:
                        self.d_atten = upload(self.vaimennus)
                    else:
                        linear_norm = (self.BPType == 4 and not self.CT)
                        self.d_atten = _cupy_texture(
                            self.vaimennus.reshape((self.Nz[0].item(), self.Ny[0].item(), self.Nx[0].item())),
                            (self.Nx[0].item(), self.Ny[0].item(), self.Nz[0].item()),
                            linear=linear_norm, normalized=linear_norm)
                if self.useMaskFP:
                    if not self.useImages:
                        self.d_maskFP = upload(self.maskFP)
                    else:
                        self.maskFP = self.maskFP.ravel('F')
                        if self.SPECT and self.maskFPZ > 1 and self.maskFPZ == self.nHeads:
                            self.d_maskFP = _cupy_texture(
                                self.maskFP.reshape((self.nHeads, self.nColsD, self.nRowsD)),
                                (self.nRowsD, self.nColsD, self.nHeads),
                                linear=False, normalized=False, uint8=True)
                        elif self.maskFPZ > 1:
                            self.d_maskFP = [None] * self.subsets
                            # maskFP has already been reordered into subset-contiguous projection
                            # order by prepass.parseInputs (only done for subsetType >= 8); slice it
                            # per subset using the per-subset projection counts rather than the
                            # cumulative nMeas offsets (which are start boundaries, not per-subset
                            # depths, and previously caused every subset to reread from the start of
                            # the full mask array with a mismatched depth).
                            maskFP3 = self.maskFP.reshape((self.nRowsD, self.nColsD, -1), order='F')
                            offsets = np.concatenate(([0], np.cumsum(np.asarray(self.nProjSubset[0]))))
                            for i in range(self.subsets):
                                start = int(offsets[i])
                                stop = int(offsets[i + 1])
                                depth = stop - start
                                subMask = np.ascontiguousarray(np.transpose(maskFP3[:, :, start:stop], (2, 1, 0)))
                                self.d_maskFP[i] = _cupy_texture(
                                    subMask, (self.nRowsD, self.nColsD, depth),
                                    linear=False, normalized=False, uint8=True)
                        else:
                            self.d_maskFP = _cupy_texture(
                                self.maskFP.reshape((self.nColsD, self.nRowsD)),
                                (self.nRowsD, self.nColsD),
                                linear=False, normalized=False, uint8=True)
                if self.useMaskBP:
                    if not self.useImages:
                        self.d_maskBP = upload(self.maskBP)
                    else:
                        self.maskBP = self.maskBP.ravel('F')
                        normalized_bp = (self.BPType == 4 and not self.CT)
                        if self.maskBPZ > 1:
                            self.d_maskBP = _cupy_texture(
                                self.maskBP.reshape((self.maskBPZ, self.Ny[0].item(), self.Nx[0].item())),
                                (self.Nx[0].item(), self.Ny[0].item(), self.maskBPZ),
                                linear=False, normalized=normalized_bp, uint8=True)
                        else:
                            self.d_maskBP = _cupy_texture(
                                self.maskBP.reshape((self.Ny[0].item(), self.Nx[0].item())),
                                (self.Nx[0].item(), self.Ny[0].item()),
                                linear=False, normalized=normalized_bp, uint8=True)
                if self.TOF:
                    self.d_TOFCenter = upload(self.TOFCenter)
                _initialize_detector_vector_buffers(self, upload)
                if self.SPECT:
                    self.d_rayShiftsDetector = upload(self.rayShiftsDetector)
                    self.d_rayShiftsSource = upload(self.rayShiftsSource)
                if (self.BPType == 2 or self.BPType == 3 or self.FPType == 2 or self.FPType == 3):
                    self.d_V = upload(self.V)
                if (self.normalization_correction):
                    # Normalization is per-timestep in C++/MATLAB for dynamic data; build
                    # [Nt][subsets] (mirrors mps_backend.py's d_norm) instead of a flat
                    # [subsets] list that only ever captured frame 0.
                    self.d_norm = [[None] * self.subsets for _ in range(self.Nt)]
                    # A static (single-frame) normalization is shared by every timestep, as in C++: slice it
                    # with the frame-0 offsets. It is recognised by holding exactly frame 0's measurements
                    # (after prepass permutation); a frame-concatenated one holds all Nt frames.
                    normStatic = (self.Nt > 1 and np.size(self.normalization) == self.nTotMeas[self.subsets].item())
                    for timestep in range(self.Nt):
                        for i in range(self.subsets):
                            index = i if normStatic else timestep * self.subsets + i
                            if self.SPECT and self.normZ == self.nHeads:
                                self.d_norm[timestep][i] = upload(self.normalization)
                            else:
                                self.d_norm[timestep][i] = upload(self.normalization[self.nTotMeas[index].item() : self.nTotMeas[index + 1].item()])
                if (self.additionalCorrection):
                    # Additional corrections (randoms/scatter) are per-timestep in C++/MATLAB
                    # for dynamic data; build [Nt][subsets] (mirrors mps_backend.py's d_scatter).
                    self.d_corr = [[None] * self.subsets for _ in range(self.Nt)]
                    for timestep in range(self.Nt):
                        for i in range(self.subsets):
                            index = timestep * self.subsets + i
                            self.d_corr[timestep][i] = upload(self.corrVector[self.nTotMeas[index].item() : self.nTotMeas[index + 1].item()])
                if (self.listmode != 1 and ((not self.CT and not self.SPECT and not self.PET) and (self.subsets > 1 and (self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7)))):
                    self.d_zindex = [None] * self.subsets
                    self.d_xyindex = [None] * self.subsets
                    for i in range(self.subsets):
                        self.d_xyindex[i] = upload(self.xy_index[self.nMeas[i] : self.nMeas[i + 1]])
                        self.d_zindex[i] = upload(self.z_index[self.nMeas[i] : self.nMeas[i + 1]])
                if (self.listmode > 0 and self.useIndexBasedReconstruction):
                    self.d_trIndex = [None] * self.subsets
                    self.d_axIndex = [None] * self.subsets
                    for i in range(self.subsets):
                        if self.loadTOF:
                            self.d_trIndex[i] = upload(self.trIndex[self.nMeas[i] * 2 : self.nMeas[i + 1] * 2])
                            self.d_axIndex[i] = upload(self.axIndex[self.nMeas[i] * 2 : self.nMeas[i + 1] * 2])
                if self.OffsetLimit.size > 0 and ((self.BPType == 4 and self.CT) or self.BPType == 5):
                    self.d_T = [None] * self.subsets
                    for i in range(self.subsets):
                        self.d_T[i] = upload(self.OffsetLimit[self.nMeas[i].item() : self.nMeas[i + 1].item()])
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
                    lines = _read(headerDir, 'auxKernels.cl')
                    lines = hlines + lines
                    bOpt += ('-DCAST=float','-DPSF','-DLOCAL_SIZE=' + str(localSize[0]),'-DLOCAL_SIZE2=' + str(localSize[1]),)
                    mod = cp.RawModule(code=lines, options=bOpt)
                    self.knlPSF = mod.get_function('Convolution3D_f')
                    self.d_gaussPSF = upload(self.gaussK.ravel('F'))

                # ``inf`` is the public box-support value.  Kernels receive a
                # finite sentinel instead on every backend, because ``isinf``
                # is unreliable under fast-math.  (Consumed inside
                # _build_kIndF/_build_kIndB below.)

                self.kIndF = _build_kIndF(self).as_tuple()
                self.kIndB = _build_kIndB(self).as_tuple()
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
            # Backend upload adapter (see the matching CuPy `upload` above).
            def upload(value):
                return cl.array.to_device(self.queue, value)
            _initialize_coordinate_buffers(self, upload)
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
                        self.d_geom5[timestep][subset] = upload(geom)
            if (self.attenuation_correction and not self.CTAttenuation):
                # Measurement-domain attenuation is per-timestep in C++/MATLAB (dynamic data can
                # have frame-dependent correction factors); build [Nt][subsets], mirroring
                # mps_backend.py's d_attenuation, instead of a flat [subsets] list that only ever
                # captured frame 0 via nTotMeas[0:subsets].
                self.d_atten = [[None] * self.subsets for _ in range(self.Nt)]
                for timestep in range(self.Nt):
                    for i in range(self.subsets):
                        index = timestep * self.subsets + i
                        self.d_atten[timestep][i] = upload(self.vaimennus[self.nTotMeas[index].item() : self.nTotMeas[index + 1].item()])
            elif (self.attenuation_correction and self.CTAttenuation):
                if self.useImages:
                    imformat = cl.ImageFormat(cl.channel_order.A, cl.channel_type.FLOAT)
                    self.d_atten = _cl_image(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, imformat,
                                              (self.Nx[0].item(), self.Ny[0].item(), self.Nz[0].item()), self.vaimennus)
                else:
                    self.d_atten = upload(self.vaimennus)
                # self.d_atten = cl.image_from_array(self.clctx, np.reshape(self.vaimennus, (self.Nx[0].item(), self.Ny[0].item(), self.Nz[0].item()), order='F'))
            _initialize_detector_vector_buffers(self, upload)
            if self.SPECT:
                self.d_rayShiftsDetector = upload(self.rayShiftsDetector)
                self.d_rayShiftsSource = upload(self.rayShiftsSource)
            if self.useMaskFP:
                imformat = cl.ImageFormat(cl.channel_order.A, cl.channel_type.UNSIGNED_INT8)
                if self.SPECT and self.maskFPZ > 1 and self.maskFPZ == self.nHeads:
                    self.d_maskFP = _cl_image(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, imformat,
                                               (self.nRowsD, self.nColsD, self.nHeads), self.maskFP)
                elif self.maskFPZ > 1:
                    self.d_maskFP = [None] * self.subsets
                    # maskFP has already been reordered into subset-contiguous projection order by
                    # prepass.parseInputs (only done for subsetType >= 8); slice it per subset using
                    # the per-subset projection counts rather than the cumulative nMeas offsets
                    # (which are start boundaries, not per-subset depths, and previously caused
                    # every subset image to reread from the start of the full mask array with a
                    # mismatched depth).
                    maskFP3 = np.asfortranarray(self.maskFP).reshape((self.nRowsD, self.nColsD, -1), order='F')
                    offsets = np.concatenate(([0], np.cumsum(np.asarray(self.nProjSubset[0]))))
                    for i in range(self.subsets):
                        start = int(offsets[i])
                        stop = int(offsets[i + 1])
                        depth = stop - start
                        subMask = np.asfortranarray(maskFP3[:, :, start:stop])
                        self.d_maskFP[i] = _cl_image(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, imformat,
                                                      (self.nRowsD, self.nColsD, depth), subMask)
                else:
                    self.d_maskFP = _cl_image(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, imformat,
                                               (self.nRowsD, self.nColsD), self.maskFP)
                # self.d_maskFP = cl.image_from_array(self.clctx, np.ascontiguousarray(self.maskFP))
            if self.useMaskBP:
                imformat = cl.ImageFormat(cl.channel_order.A, cl.channel_type.UNSIGNED_INT8)
                if self.maskBPZ > 1:
                    self.d_maskBP = _cl_image(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, imformat,
                                               (self.Nx[0].item(), self.Ny[0].item(), self.maskBPZ), self.maskBP)
                else:
                    self.d_maskBP = _cl_image(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, imformat,
                                               (self.Nx[0].item(), self.Ny[0].item()), self.maskBP)
                # self.d_maskBP = cl.image_from_array(self.clctx, np.ascontiguousarray(self.maskBP))
            if self.TOF:
                self.d_TOFCenter = upload(self.TOFCenter)
            if (self.BPType == 2 or self.BPType == 3 or self.FPType == 2 or self.FPType == 3):
                self.d_V = upload(self.V)
            if (self.normalization_correction):
                # Normalization is per-timestep in C++/MATLAB for dynamic data; build
                # [Nt][subsets] (mirrors mps_backend.py's d_norm) instead of a flat [subsets]
                # list that only ever captured frame 0.
                self.d_norm = [[None] * self.subsets for _ in range(self.Nt)]
                # A static (single-frame) normalization is shared by every timestep, as in C++: slice it
                # with the frame-0 offsets. It is recognised by holding exactly frame 0's measurements
                # (after prepass permutation); a frame-concatenated one holds all Nt frames.
                normStatic = (self.Nt > 1 and np.size(self.normalization) == self.nTotMeas[self.subsets].item())
                for timestep in range(self.Nt):
                    for i in range(self.subsets):
                        index = i if normStatic else timestep * self.subsets + i
                        if self.SPECT and self.normZ == self.nHeads:
                            self.d_norm[timestep][i] = upload(self.normalization)
                        else:
                            self.d_norm[timestep][i] = upload(self.normalization[self.nTotMeas[index].item() : self.nTotMeas[index + 1].item()])
            if (self.additionalCorrection):
                # Additional corrections (randoms/scatter) are per-timestep in C++/MATLAB for
                # dynamic data; build [Nt][subsets] (mirrors mps_backend.py's d_scatter).
                self.d_corr = [[None] * self.subsets for _ in range(self.Nt)]
                for timestep in range(self.Nt):
                    for i in range(self.subsets):
                        index = timestep * self.subsets + i
                        self.d_corr[timestep][i] = upload(self.corrVector[self.nTotMeas[index].item() : self.nTotMeas[index + 1].item()])
            if (self.listmode != 1 and ((not self.CT and not self.SPECT and not self.PET) and (self.subsets > 1 and (self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7)))):
                self.d_zindex = [None] * self.subsets
                self.d_xyindex = [None] * self.subsets
                for i in range(self.subsets):
                    self.d_xyindex[i] = upload(self.xy_index[self.nMeas[i] : self.nMeas[i + 1]])
                    self.d_zindex[i] = upload(self.z_index[self.nMeas[i] : self.nMeas[i + 1]])
            if (self.listmode > 0 and self.useIndexBasedReconstruction):
                self.d_trIndex = [None] * self.subsets
                self.d_axIndex = [None] * self.subsets
                for i in range(self.subsets):
                    if self.loadTOF:
                        self.d_trIndex[i] = upload(self.trIndex[self.nMeas[i] * 2 : self.nMeas[i + 1] * 2])
                        self.d_axIndex[i] = upload(self.axIndex[self.nMeas[i] * 2 : self.nMeas[i + 1] * 2])
            if self.OffsetLimit.size > 0 and ((self.BPType == 4 and self.CT) or self.BPType == 5):
                self.d_T = [None] * self.subsets
                for i in range(self.subsets):
                    self.d_T[i] = upload(self.OffsetLimit[self.nMeas[i].item() : self.nMeas[i + 1].item()])
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
                lines = _read(headerDir, 'auxKernels.cl')
                lines = hlines + lines
                bOpt +=(' -DCAST=float',' -DPSF',' -DLOCAL_SIZE=' + str(localSize[0]), ' -DLOCAL_SIZE2=' + str(localSize[1]),)
                prg = cl.Program(self.clctx, lines).build(' '.join(bOpt))
                self.knlPSF = prg.Convolution3D_f
                self.d_gaussPSF = upload(self.gaussK.ravel('F'))
                
            fp_args = _build_kIndF(self)
            self.kIndF = fp_args.apply_opencl(self.knlF, 0)

            bp_args = _build_kIndB(self)
            self.kIndB = bp_args.apply_opencl(self.knlB, 0)
