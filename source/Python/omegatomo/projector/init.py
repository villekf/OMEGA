# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 13:17:22 2025
"""
import re
import numpy as np


def _kernel_ellipse_power(value):
    # Kernels test for box support (ellipsePower = inf) with a finite threshold,
    # since isinf() is unreliable under fast-math
    value = float(value)
    return float(np.finfo(np.float32).max) if not np.isfinite(value) else value


def _kernel_crystal_pitch(self):
    """The float2 handed to the kernels as crystalSize / dPitch.

    The kernels use .x as the transaxial and .y as the axial extent (multi-ray spread,
    detector pitch). For PET that is the crystal size (cr_p, cr_pz); dPitchX/dPitchY
    hold (cr_pz, cr_p) there for historical reasons and must not be passed as such.
    CT and SPECT use the detector pitch (dPitchX, dPitchY).
    """
    if self.CT or self.SPECT:
        return self.dPitchX, self.dPitchY
    return self.cr_p, self.cr_pz


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
        a.vec2f(*_kernel_crystal_pitch(self))
        if self.FPType in (2, 3):
            if self.FPType == 2:
                a.f32(self.tube_width_z)
            else:
                a.f32(self.tube_radius)
            a.f32(self.bmin).f32(self.bmax).f32(self.Vmax)
    elif self.FPType == 4:
        a.u32(self.nRowsD).u32(self.nColsD).vec2f(*_kernel_crystal_pitch(self))
        if self.CT and self.useHelical:
            a.f32(self.helicalRadius)
        a.f32(self.dL).f32(self.global_factor)
    elif self.FPType == 5:
        a.u32(self.nRowsD).u32(self.nColsD).vec2f(*_kernel_crystal_pitch(self))
        if self.CT and self.useHelical:
            a.f32(self.helicalRadius)
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
        # The type 4 kernel (projectorType4.cl) takes TOFCenter, sigma_x and d_epps as consecutive arguments (as in ProjectorClass.h)
        a.buf(self.d_TOFCenter)
        a.f32(self.sigma_x)
        a.f32(self.epps)
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
        a.u32(self.nRowsD).u32(self.nColsD).vec2f(*_kernel_crystal_pitch(self))
        if self.CT and self.useHelical:
            a.f32(self.helicalRadius)
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
        a.vec2f(*_kernel_crystal_pitch(self))
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
        a.f32(self.epps)
    if self.attenuation_correction and self.CTAttenuation and self.BPType in (1, 2, 3, 4) and not self.CT:
        if self.useImages or self.useCUDA:
            a.img(self.d_atten)
        else:
            a.buf(self.d_atten)
    return a


def _kernel_ravel(arr, default_order, components=0):
    """Flatten a coordinate array so that the per-ray components are consecutive.

    The kernels read the components of one ray (4 for PET xy, 6 for CT/SPECT/list-mode, 1-3 for z)
    as consecutive floats, i.e. the component axis has to vary fastest. PET and CT arrays are
    stored as (rays, components), SPECT and list-mode arrays as (components, rays). The component
    axis is the one equal to `components` if given, otherwise the smaller one. Square or otherwise
    ambiguous arrays use `default_order`.
    """
    arr = np.asarray(arr, dtype=np.float32)
    order = default_order
    if arr.ndim == 2:
        if components and (arr.shape[0] == components) != (arr.shape[1] == components):
            order = 'F' if arr.shape[0] == components else 'C'
        elif arr.shape[1] < arr.shape[0]:
            order = 'C'
        elif arr.shape[0] < arr.shape[1]:
            order = 'F'
    return arr.ravel(order=order)


def _coordinate_slice(self, name, timestep, subset, stride):
    """Return one frame/subset coordinate slice in kernel order."""
    frames = getattr(self, name + 'Frames', None)
    if isinstance(frames, list) and len(frames) == int(self.Nt):
        frame = _kernel_ravel(frames[timestep], 'F', stride)
        offsets = np.concatenate(([0], np.cumsum(self.nProjSubset[timestep], dtype=np.int64)))
        start = int(offsets[subset]) * stride
        stop = int(offsets[subset + 1]) * stride
        return frame[start:stop]
    # CT stores one projection per row; SPECT stores one per column.
    order = 'C' if self.CT and self.listmode == 0 else 'F'
    flat = _kernel_ravel(getattr(self, name), order, stride)
    q = timestep * int(self.subsets) + subset
    # Non-list-mode geometry without per-frame arrays (e.g. PET sinograms) is the same in every frame
    # and is stored once, i.e. for one frame only, so the subset offsets are counted from the start of
    # the frame. nMeas holds the cumulative counts over all frames, which is the right offset only for
    # the geometry that is stored for all frames (list-mode coordinates).
    frameStart = int(self.nMeas[timestep * int(self.subsets)]) if self.listmode == 0 else 0
    start = (int(self.nMeas[q]) - frameStart) * stride
    stop = (int(self.nMeas[q + 1]) - frameStart) * stride
    return flat[start:stop]


def _full_coordinate_frame(self, name, timestep):
    frames = getattr(self, name + 'Frames', None)
    if isinstance(frames, list) and len(frames) == int(self.Nt):
        return _kernel_ravel(frames[timestep], 'F')
    return _kernel_ravel(getattr(self, name), 'F')


def _initialize_coordinate_buffers(self, upload):
    """Create the canonical ``[timestep][subset]`` geometry buffers."""
    self.d_x = [[None] * self.subsets for _ in range(self.Nt)]
    self.d_z = [[None] * self.subsets for _ in range(self.Nt)]
    # Curved helical CT stores a single value (the projection angle) per projection in z
    if self.CT and self.useHelical:
        z_stride = 1
    else:
        # PET sinograms of a multi-layer scanner carry the layer combination as the third value of each sinogram
        z_stride = 6 if self.pitch else (3 if self.PET and self.listmode == 0 and self.nLayers > 1 else 2)
    subset_geometry = (self.CT or self.SPECT) and self.listmode == 0
    pet_geometry = self.PET and self.listmode == 0
    for timestep in range(self.Nt):
        if subset_geometry or (self.listmode > 0 and not self.useIndexBasedReconstruction and self.loadTOF):
            for subset in range(self.subsets):
                self.d_x[timestep][subset] = upload(_coordinate_slice(self, 'x', timestep, subset, 6))
        else:
            self.d_x[timestep][0] = upload(_full_coordinate_frame(self, 'x', timestep))

        if subset_geometry or pet_geometry:
            for subset in range(self.subsets):
                self.d_z[timestep][subset] = upload(_coordinate_slice(self, 'z', timestep, subset, z_stride))
        elif self.listmode == 0 or (self.listmode > 0 and self.useIndexBasedReconstruction):
            self.d_z[timestep][0] = upload(_full_coordinate_frame(self, 'z', timestep))
        else:
            self.d_z[timestep][0] = upload(np.zeros(1, dtype=np.float32))


def _tof_index_slice(self, timestep, subset):
    """The TOF bin indices (uint8, one per event) of one frame/subset of list-mode data. The events are ordered like the event
    coordinates and the index-based detector indices, i.e. the subset is the range [nMeas[q], nMeas[q + 1]) of the events,
    q = timestep * subsets + subset."""
    q = timestep * int(self.subsets) + subset
    start = int(self.nMeas[q])
    stop = int(self.nMeas[q + 1])
    tof = np.asarray(self.TOFIndices).ravel(order='F')
    if tof.size < stop:
        raise ValueError('The number of TOF indices (%d) is smaller than the number of events (%d)!' % (tof.size, stop))
    values = tof[start:stop]
    if self.TOF_bins_used == self.TOF_bins and values.size > 0 and (int(values.min()) < 0 or int(values.max()) >= self.TOF_bins):
        raise ValueError('TOF indices must be in the range 0 ... TOF_bins - 1 (%d)!' % (self.TOF_bins - 1))
    return np.ascontiguousarray(values, dtype=np.uint8)


def _initialize_tof_index_buffers(self, upload):
    """Create the canonical ``[timestep][subset]`` TOF index buffers of list-mode data (the uint8 TOFIndex kernel argument that
    selects the one TOF bin of each event, present when the kernels are built with both -DLISTMODE and -DTOF). Same layout as
    the event coordinates: with loadTOF every subset is uploaded here, otherwise only the entry [timestep][0] is used and
    projfunctions._load_listmode_subset replaces it with the subset that is about to be processed."""
    self.d_TOFIndex = [[None] * self.subsets for _ in range(self.Nt)]
    if self.listmode > 0 and self.TOF and self.loadTOF:
        for timestep in range(self.Nt):
            for subset in range(self.subsets):
                self.d_TOFIndex[timestep][subset] = upload(_tof_index_slice(self, timestep, subset))


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


def computeProjectionWindow(x, angles, r, nColsD, dPitchY, bx, by, bz, dx, dy, dz, Nx, Ny, Nz, nVoxels):
    """
    For each axial voxel block (nVoxels voxels in the axial direction) of a volume, the range of projections that can reach the
    curved helical detector (-DPROJWINDOW in projectorType4.cl and projectorType5.cl). Every projection is tested, so the range is
    correct for any projection order, though it only stays narrow when the source height is monotonic. Returns 2 int32 values per
    block: the first projection and one past the last projection (both 0 if no projection sees the block).
    """
    nProj = x.size // 6
    xs = np.reshape(x, (nProj, 6)).astype(np.float64, copy=False)
    a = np.asarray(angles, dtype=np.float64).ravel()
    r = float(r)
    sX = xs[:, 0]
    sY = xs[:, 1]
    sZ = xs[:, 2]
    # Arc centre. The source is inside the detector cylinder, so a ray from the source to a point of the image box has to travel at
    # least r - |s - C| in the transaxial plane, relative to the largest transaxial distance from the source to the box
    cX = xs[:, 3] - r * np.cos(a)
    cY = xs[:, 4] + r * np.sin(a)
    dSC = np.hypot(sX - cX, sY - cY)
    maxDist = np.zeros(nProj, dtype=np.float64)
    for cornerX in (bx, bx + Nx * dx):
        for cornerY in (by, by + Ny * dy):
            maxDist = np.maximum(maxDist, np.hypot(cornerX - sX, cornerY - sY))
    tMin = np.zeros(nProj, dtype=np.float64)
    np.divide(r - dSC, maxDist, out=tMin, where=maxDist > 0)
    # Axial extent of the detector with one row of margin on each side
    halfWidth = (nColsD / 2. + 1.) * dPitchY
    zMin = xs[:, 5] - halfWidth
    zMax = xs[:, 5] + halfWidth
    nBlocks = (int(Nz) + int(nVoxels) - 1) // int(nVoxels)
    zEnd = bz + Nz * dz
    window = np.zeros(2 * nBlocks, dtype=np.int32)
    for b in range(nBlocks):
        zLow = bz + b * nVoxels * dz
        zHigh = min(zLow + nVoxels * dz, zEnd)
        # A ray hits the cylinder at least tMin of the way (in the axial direction) from the source to the block
        invisible = (((zLow > sZ) & (sZ + tMin * (zLow - sZ) > zMax))
                     | ((zHigh < sZ) & (sZ + tMin * (zHigh - sZ) < zMin)))
        visible = np.flatnonzero(~invisible)
        if visible.size > 0:
            window[2 * b] = visible[0]
            window[2 * b + 1] = visible[-1] + 1
    return window


def _uses_projection_window(self):
    """True if and only if the curved helical BP of type 4 or 5 only visits the projections in the window of each axial voxel block."""
    return bool(self.CT and self.useHelical and not self.useParallelBeam and not self.useMetal and self.listmode == 0
                and self.subsetType != 12 and self.BPType in (4, 5))


def _initialize_projection_window_buffers(self, upload):
    """Create the ``[timestep][subset][volume]`` projection window buffers (see computeProjectionWindow)."""
    self.d_projWindow = [[[None] * (self.nMultiVolumes + 1) for _ in range(self.subsets)] for _ in range(self.Nt)]
    nVoxels = self.NVOXELS5 if self.BPType == 5 else self.NVOXELSHELICAL
    for timestep in range(self.Nt):
        for subset in range(self.subsets):
            xs = _coordinate_slice(self, 'x', timestep, subset, 6)
            angles = _coordinate_slice(self, 'z', timestep, subset, 1)
            for ii in range(self.nMultiVolumes + 1):
                window = computeProjectionWindow(
                    xs,
                    angles,
                    self.helicalRadius,
                    self.nColsD,
                    self.dPitchY,
                    self.bx[ii].item(),
                    self.by[ii].item(),
                    self.bz[ii].item(),
                    self.dx[ii].item(),
                    self.dy[ii].item(),
                    self.dz[ii].item(),
                    self.Nx[ii].item(),
                    self.Ny[ii].item(),
                    self.Nz[ii].item(),
                    nVoxels,
                )
                self.d_projWindow[timestep][subset][ii] = upload(window)

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


def _userLocalSize(self, default):
    """The local (work-group/block) size of the projector kernels: the default overridden per dimension by the positive
    entries of the user-supplied local_size option (the same semantics as in the C++ library; negative entries keep the default)."""
    user = np.atleast_1d(np.asarray(self.local_size)).ravel()
    out = list(default)
    for k in range(min(3, user.size)):
        if user[k] > 0:
            out[k] = int(user[k])
    return tuple(out)


def _isWave64(vendor=None, gfxipMajor=None, archName=None):
    """Whether the GPU is a wave64 AMD GPU (GCN/CDNA). Same rule as isWave64HipArch/isWave64OpenCL in ProjectorClass.h.

    HIP/ROCm (archName, the gcnArchName of the device): the part before the first ':' (the rest are feature suffixes such as
    gfx90a:sramecc+:xnack-) is gfx + 3 characters (gfx803, gfx90a, gfx942) for wave64 and gfx + 4 characters (gfx1030, gfx1100) for wave32.
    Anything unparsable is wave32.
    OpenCL (vendor and gfxipMajor, CL_DEVICE_GFXIP_MAJOR_AMD, None if the query is not available): AMD with a major version below 10.
    Everything else (CUDA, Metal, NVIDIA/Intel OpenCL, CPUs) is wave32."""
    if archName is not None:
        if isinstance(archName, bytes):
            archName = archName.decode('ascii', 'ignore')
        match = re.fullmatch(r'gfx([0-9a-fA-F]+)', str(archName).split(':')[0].rstrip('\0'))
        return match is not None and len(match.group(1)) == 3
    if vendor == 'Advanced Micro Devices, Inc.' and gfxipMajor is not None:
        try:
            return 0 < int(gfxipMajor) < 10
        except (TypeError, ValueError):
            return False
    return False


def _listmodeLocalSize(vendor=None, gfxipMajor=None, archName=None, deviceLimit=None):
    """The default 1D local (work-group/block) size of list-mode (and index-based) data. The events are unsorted, so 128 is a poor size:
    32 is the best on wave32 devices and 512 on the wave64 AMD GPUs (see _isWave64). Capped to the device limit (the maximum work-group size,
    None if unknown). Same rule as listmodeLocalSize in ProjectorClass.h."""
    size = 512 if _isWave64(vendor, gfxipMajor, archName) else 32
    try:
        if deviceLimit is not None and 0 < int(deviceLimit) < size:
            size = int(deviceLimit)
    except (TypeError, ValueError):
        pass
    return size


def _deviceListmodeLocalSize(self, rocm=False, cp=None):
    """_listmodeLocalSize for the device in use: the CuPy device on ROCm/HIP (rocm and the cupy module cp), the PyOpenCL device in
    self.clctx, 32 for everything else (CUDA, Metal/MPS, PyTorch) or if the device properties cannot be read."""
    if self.useMetal:
        return _listmodeLocalSize()
    try:
        if self.useCUDA:
            if not rocm or cp is None:
                return _listmodeLocalSize()
            props = cp.cuda.runtime.getDeviceProperties(cp.cuda.runtime.getDevice())
            return _listmodeLocalSize(archName=props.get('gcnArchName'), deviceLimit=props.get('maxThreadsPerBlock'))
        import pyopencl as cl
        device = self.clctx.get_info(cl.context_info.DEVICES)[0]
        vendor = device.get_info(cl.device_info.VENDOR)
        gfxipMajor = None
        # The query is an AMD extension, not every PyOpenCL version has the attribute
        if vendor == 'Advanced Micro Devices, Inc.' and 'cl_amd_device_attribute_query' in device.get_info(cl.device_info.EXTENSIONS) \
                and hasattr(cl.device_info, 'GFXIP_MAJOR_AMD'):
            gfxipMajor = device.get_info(cl.device_info.GFXIP_MAJOR_AMD)
        return _listmodeLocalSize(vendor=vendor, gfxipMajor=gfxipMajor, deviceLimit=device.max_work_group_size)
    except Exception:
        return _listmodeLocalSize()


def _orth_clip_needed(self):
    """Whether the projector type 2/3 kernels are built with the end point check of the tube (-DORTHCLIP, see orthogonalHelper3D in
    opencl_functions_orth3D.h). The tube of a ray ends at its end points, which only matters when a voxel beyond an end point can be
    reached by the tube, i.e. when a detector/source is inside the image.

    Rule (same as orthClipNeeded in precomp.h): enabled if and only if any detector/source coordinate point has its x and y inside the image box
    widened by m, and the z range of the coordinates overlaps the box widened by m in z. m is the tube reach (type 2: tube_width_z,
    type 3: max(bmax, tube_radius)) plus the largest voxel dimension, the box is the union of the volumes. Voxel centers are inside the
    box, so a voxel in the tube is closer than the reach to the line; for detectors around the image the line beyond an end point is
    outside the (convex) widened box and no voxel there is in the tube. The rule is conservative (it errs on the side of enabling the
    check), non-finite coordinates (multi-layer padding) are ignored. CT always enables it, SPECT never (the kernel skips it)."""
    if self.SPECT:
        return False
    if self.CT:
        return True
    reach = 0.0
    if self.FPType == 2 or self.BPType == 2:
        reach = max(reach, float(self.tube_width_z))
    if self.FPType == 3 or self.BPType == 3:
        reach = max(reach, float(self.tube_radius), float(self.bmax))
    nVol = int(self.nMultiVolumes) + 1
    lo = [np.inf] * 3
    hi = [-np.inf] * 3
    voxel = 0.0
    for k in range(nVol):
        for a, (b, d, n) in enumerate(((self.bx, self.dx, self.Nx), (self.by, self.dy, self.Ny), (self.bz, self.dz, self.Nz))):
            bk, dk, nk = float(np.ravel(b)[k]), float(np.ravel(d)[k]), float(np.ravel(n)[k])
            lo[a] = min(lo[a], bk)
            hi[a] = max(hi[a], bk + nk * dk)
            voxel = max(voxel, dk)
    m = reach + voxel
    lo = [v - m for v in lo]
    hi = [v + m for v in hi]

    def frames(name):
        f = getattr(self, name + 'Frames', None)
        if isinstance(f, list) and len(f) > 0:
            return f
        return [getattr(self, name)]

    xy = []
    z = []
    if self.listmode > 0 and not self.useIndexBasedReconstruction:
        # List-mode events: x, y, z of both end points of each event
        for f in frames('x'):
            p = _kernel_ravel(f, 'F', 6)
            p = p[:p.size // 6 * 6].reshape(-1, 6)
            xy += [p[:, 0:2], p[:, 3:5]]
            z += [p[:, 2], p[:, 5]]
    else:
        for f in frames('x'):
            p = _kernel_ravel(f, 'F')
            xy.append(p[:p.size // 2 * 2].reshape(-1, 2))
        for f in frames('z'):
            if self.listmode == 0 and not self.use_raw_data:
                # PET sinograms: source and detector z of each sinogram, with multiple layers the third value is the layer combination
                stride = 3 if self.nLayers > 1 else 2
                p = _kernel_ravel(f, 'F', stride)
                z.append(p[:p.size // stride * stride].reshape(-1, stride)[:, 0:2].ravel())
            else:
                z.append(_kernel_ravel(f, 'F'))
    xy = np.concatenate(xy) if len(xy) > 0 else np.empty((0, 2), dtype=np.float32)
    xy = xy[np.isfinite(xy).all(axis=1)]
    if xy.shape[0] == 0:
        # No coordinates available: nothing can be excluded
        return True
    inside = bool(np.any((xy[:, 0] >= lo[0]) & (xy[:, 0] <= hi[0]) & (xy[:, 1] >= lo[1]) & (xy[:, 1] <= hi[1])))
    if not inside:
        return False
    z = np.concatenate(z) if len(z) > 0 else np.empty(0, dtype=np.float32)
    z = z[np.isfinite(z)]
    if z.size == 0:
        # No z coordinates: the z range cannot exclude anything
        return True
    return bool(z.min() <= hi[2] and z.max() >= lo[2])


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
    if self.useCUDA and not self.useCuPy and not self.useMetal:
        # The operator interface (forwardProject/backwardProject) has no PyCUDA path, CuPy is the CUDA array library
        print('The CUDA operator interface requires CuPy. Setting useCuPy to True!')
        self.useCuPy = True
    if getattr(self, 'useRowMajor', None) and (self.useAF or self.useMetal):
        print('useRowMajor is not supported with ArrayFire or Metal/MPS, the operator output stays flat and Fortran-ordered.')
    if self.useCuPy and self.useCUDA:
        import cupy as cp
    elif self.useCuPy and not self.useCUDA:
        print('CuPy can only be used when useCUDA is True. Setting useCUDA to True!')
        self.useCUDA = True
        import cupy as cp
    # Always defined (also when CuPy is off), it is called unconditionally further below on the CUDA path
    def cupyROCm():
        if not (self.useCuPy and self.useCUDA):
            return False
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
        import pyopencl.array

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
    # Number of z-voxels each thread processes in the curved helical BP (mirrors NVOXELSHELICAL in the C++ headers)
    self.NVOXELSHELICAL = 1
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
        # List-mode (and index-based) data: 32, or 512 on wave64 AMD GPUs, instead of 128. Other data using the 1D launch keeps 128
        defaultLocalSize1D = 128
        if self.listmode > 0:
            cupyOnCUDA = self.useCuPy and self.useCUDA
            defaultLocalSize1D = _deviceListmodeLocalSize(self, cupyOnCUDA and cupyROCm(), cp if cupyOnCUDA else None)
        globalSize = [[None] * self.subsets for _ in range(frame_count)]
        # self.mSize = [None] * self.subsets
        for timestep in range(frame_count):
            for i in range(self.subsets):
                n_proj = int(self.nProjSubset[timestep, i])
                n_meas = int(self.nMeasSubset[timestep, i])
                if (self.FPType == 5):
                    globalSize[timestep][i] = (self.nRowsD, (self.nColsD + self.NVOXELSFP - 1) // self.NVOXELSFP, n_proj)
                    localSize = _userLocalSize(self, (16, 16, 1))
                    erotus = (localSize[0] - (globalSize[timestep][i][0] % localSize[0]), localSize[1] - (globalSize[timestep][i][1] % localSize[1]), 0)
                    globalSize[timestep][i] = (self.nRowsD + erotus[0], (self.nColsD + self.NVOXELSFP - 1) // self.NVOXELSFP + erotus[1], n_proj)
                elif ((self.CT or self.SPECT or self.PET) and self.listmode == 0):
                    globalSize[timestep][i] = (self.nRowsD, self.nColsD, n_proj)
                    # The multi-bin TOF kernels need many registers per thread, smaller blocks (16 x 8 instead of 16 x 16) give
                    # a higher occupancy (measured on an NVIDIA RTX 4090 with CUDA and OpenCL: up to 23 % faster, never slower). Same rule as in ProjectorClass.h
                    localSize = _userLocalSize(self, (16, 8, 1) if (self.TOF and self.TOF_bins_used > 1) else (16, 16, 1))
                    erotus = (localSize[0] - (globalSize[timestep][i][0] % localSize[0]), localSize[1] - (globalSize[timestep][i][1] % localSize[1]), 0)
                    globalSize[timestep][i] = (self.nRowsD + erotus[0], self.nColsD + erotus[1], n_proj)
                else:
                    globalSize[timestep][i] = (n_meas, 1, 1)
                    localSize = _userLocalSize(self, (defaultLocalSize1D, 1, 1))
                    erotus = (localSize[0] - (globalSize[timestep][i][0] % localSize[0]), localSize[1] - (globalSize[timestep][i][1] % localSize[1]), 0)
                    globalSize[timestep][i] = (n_meas + erotus[0], 1, 1)
        self.globalSizeFP = globalSize.copy()
        self.localSizeFP = localSize + tuple()
        self.erotusBP = [0] * (self.nMultiVolumes + 1) * 2
        localSize = _userLocalSize(self, (16, 16, 1))
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
                            nVoxelsZ = self.NVOXELSHELICAL if self.useHelical else self.NVOXELS
                            globalSize[timestep][subset][ii] = (self.Nx[ii].item() + self.erotusBP[ii * 2], self.Ny[ii].item() + self.erotusBP[ii * 2 + 1], (self.Nz[ii].item()  + nVoxelsZ - 1) // nVoxelsZ)
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
        # Fast math and the projector type 2/3 discontinuities: see the note at DIVIDE in general_opencl_functions.h
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
        orthClip = False
        if (self.FPType == 2 or self.BPType == 2 or self.FPType == 3 or self.BPType == 3):
            # End point check of the tube, only needed when a detector is inside the image (see _orth_clip_needed)
            orthClip = _orth_clip_needed(self)
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
            # The number of layers. For sinograms the layer combination of each sinogram is stored in z
            bOpt += ('-DNLAYERS=' + str(self.nLayers),)
        if self.TOF:
            bOpt += ('-DTOF',)
        if self.CT:
            bOpt += ('-DCT',)
            # Curved (helical) detector; the parallel-beam define takes precedence as in the C++ host code
            if self.useHelical and not self.useParallelBeam:
                bOpt += ('-DHELICAL',)
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
                if orthClip:
                    bOptFP += ('-DORTHCLIP',)
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
            if self.CT and self.useHelical:
                bOptFP += ('-DNVOXELSHELICAL=' + str(self.NVOXELSHELICAL),)
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
                if orthClip:
                    bOptBP += ('-DORTHCLIP',)
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
                if self.useHelical:
                    bOptBP += ('-DNVOXELSHELICAL=' + str(self.NVOXELSHELICAL),)
                # The helical backprojection only visits the projections in the window of its axial block
                if _uses_projection_window(self):
                    bOptBP += ('-DPROJWINDOW',)
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
            if self.CT and self.listmode == 0 and not self.useHelical:
                bOptBP += ('-DGEOM5',)
            # The helical backprojection only visits the projections in the window of its axial block
            if _uses_projection_window(self):
                bOptBP += ('-DPROJWINDOW',)
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
                _initialize_tof_index_buffers(self, upload)
                # Precomputed per-projection geometry for the BDD backprojection (see -DGEOM5 in projectorType5.cl)
                if self.BPType == 5 and self.CT and self.listmode == 0 and not self.useHelical:
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
                # Projection windows of the helical backprojection (see -DPROJWINDOW in projectorType4.cl and projectorType5.cl)
                if _uses_projection_window(self):
                    _initialize_projection_window_buffers(self, upload)
                if (self.attenuation_correction and not self.CTAttenuation):
                    # Measurement-domain attenuation: build [Nt][subsets], mirroring mps_backend.py's
                    # d_attenuation. As in C++ (a single d_atten[kk] buffer shared by all timesteps), a
                    # static (single-frame) attenuation is shared by every timestep and sliced with the
                    # frame-0 offsets. It is recognised by holding exactly frame 0's measurements (after
                    # prepass permutation); a frame-concatenated one (Nt blocks, a superset of C++) is
                    # sliced per frame.
                    self.d_atten = [[None] * self.subsets for _ in range(self.Nt)]
                    attenStatic = (self.Nt > 1 and np.size(self.vaimennus) == self.nTotMeas[self.subsets].item())
                    for timestep in range(self.Nt):
                        for i in range(self.subsets):
                            index = i if attenStatic else timestep * self.subsets + i
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
                    # Normalization: build [Nt][subsets] (mirrors mps_backend.py's d_norm).
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
                    # Additional corrections (scatter/corrVector): build [Nt][subsets] (mirrors
                    # mps_backend.py's d_scatter). As in C++, a corrVector holding more than one frame
                    # is per-timestep (frame-concatenated, sliced per frame); a single frame's worth
                    # (exactly frame 0's measurements after prepass permutation) is shared by every
                    # timestep and sliced with the frame-0 offsets.
                    self.d_corr = [[None] * self.subsets for _ in range(self.Nt)]
                    corrStatic = (self.Nt > 1 and np.size(self.corrVector) == self.nTotMeas[self.subsets].item())
                    for timestep in range(self.Nt):
                        for i in range(self.subsets):
                            index = i if corrStatic else timestep * self.subsets + i
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
            _initialize_tof_index_buffers(self, upload)
            # Precomputed per-projection geometry for the BDD backprojection (see -DGEOM5 in projectorType5.cl)
            if self.BPType == 5 and self.CT and self.listmode == 0 and not self.useHelical:
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
            # Projection windows of the helical backprojection (see -DPROJWINDOW in projectorType4.cl and projectorType5.cl)
            if _uses_projection_window(self):
                _initialize_projection_window_buffers(self, upload)
            if (self.attenuation_correction and not self.CTAttenuation):
                # Measurement-domain attenuation: build [Nt][subsets], mirroring mps_backend.py's
                # d_attenuation. As in C++ (a single d_atten[kk] buffer shared by all timesteps), a
                # static (single-frame) attenuation is shared by every timestep and sliced with the
                # frame-0 offsets. It is recognised by holding exactly frame 0's measurements (after
                # prepass permutation); a frame-concatenated one (Nt blocks, a superset of C++) is
                # sliced per frame.
                self.d_atten = [[None] * self.subsets for _ in range(self.Nt)]
                attenStatic = (self.Nt > 1 and np.size(self.vaimennus) == self.nTotMeas[self.subsets].item())
                for timestep in range(self.Nt):
                    for i in range(self.subsets):
                        index = i if attenStatic else timestep * self.subsets + i
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
                # Normalization: build [Nt][subsets] (mirrors mps_backend.py's d_norm).
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
                # Additional corrections (scatter/corrVector): build [Nt][subsets] (mirrors
                # mps_backend.py's d_scatter). As in C++, a corrVector holding more than one frame
                # is per-timestep (frame-concatenated, sliced per frame); a single frame's worth
                # (exactly frame 0's measurements after prepass permutation) is shared by every
                # timestep and sliced with the frame-0 offsets.
                self.d_corr = [[None] * self.subsets for _ in range(self.Nt)]
                corrStatic = (self.Nt > 1 and np.size(self.corrVector) == self.nTotMeas[self.subsets].item())
                for timestep in range(self.Nt):
                    for i in range(self.subsets):
                        index = i if corrStatic else timestep * self.subsets + i
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
