# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 13:25:14 2025
"""
from __future__ import annotations

import numpy as np
class _KernelArgs:
    """Backend-agnostic kernel-argument sequence builder.

    A projector family (FP types 1-3, FP type 4, FP type 5, BP types 1-3,
    BP type 4 non-CT, BP type 4 CT, BP type 5 CT) has exactly one argument
    ORDER, but three different call conventions to express it in:

      * plain PyOpenCL: ``knl.set_arg(i, ...)`` from some starting index
      * ArrayFire (still an OpenCL kernel under the hood): same as above,
        except device buffers are ArrayFire arrays that must be unwrapped
        to a ``cl.MemoryObject`` via ``x.raw_ptr()`` first
      * CuPy (CUDA/HIP, built with ``-DPYTHON``): one Python tuple passed
        straight to the compiled ``RawKernel``/``RawModule`` function

    and, on the CUDA side only, every OpenCL vector-typed argument
    (float3/uint3/float2/int3, matching the C++ kernel signatures in
    opencl/) is split into its scalar components -- see the vecNx methods.

    This class collects ONE sequence of logical arguments (built by
    calling one of the methods below, in kernel-argument order, once per
    argument) and only decides HOW to hand it to the kernel at the end,
    via ``apply_opencl`` (plain PyOpenCL / ArrayFire) or ``as_tuple``
    (CuPy). Call sites still perform their own backend-specific resource
    creation (image/texture upload, output allocation, the AF unlock
    calls, the 64/32-bit atomic conversion) exactly as before -- this
    class only replaces the parallel OpenCL-``set_arg``-chain /
    CuPy-tuple-concatenation duplication.
    """

    __slots__ = ('_useCUDA', '_useAF', '_useTorch', '_cp', '_cl', '_items')

    def __init__(self, owner):
        self._useCUDA = bool(owner.useCUDA)
        self._useAF = bool(owner.useAF)
        self._useTorch = bool(owner.useTorch)
        self._items = []
        if self._useCUDA:
            import cupy as cp
            self._cp = cp
            self._cl = None
        else:
            import pyopencl as cl
            self._cl = cl
            self._cp = None

    # -- device buffers / images ------------------------------------------------
    def buf(self, x):
        """A plain, ALREADY-init-time-uploaded device array (d_x/d_z/d_norm/
        d_corr/d_Sens/d_atten-as-buffer/d_xyindex/... -- every one of these
        is a ``pyopencl.array.Array`` under BOTH plain PyOpenCL and AF (AF
        only swaps context/queue creation at init time, never how these
        constant buffers themselves get uploaded), and a plain CuPy array
        under CuPy -- so OpenCL and AF share the exact same ``x.data``
        handling here.

        The genuinely ArrayFire-*native* per-call resources (the live f/y
        images/measurement vectors) are NOT passed through this method:
        unwrapping those needs a one-shot, explicitly-unlocked
        ``x.raw_ptr()`` -> ``cl.MemoryObject.from_int_ptr`` call, and the
        original code always performs that ONCE per forward/backward call
        (not once per multi-volume loop iteration or per argument), then
        reuses the resulting MemoryObject/CuPy-view for every argument that
        needs it. Callers resolve those the same way here (unchanged) and
        pass the resolved value through ``img()`` instead, which preserves
        the exact original raw_ptr()/lock-event count.

        CuPy+Torch -> ``cp.asarray(x)`` (zero-copy view), matching every
        existing ``fD = cp.asarray(f)`` / ``yD = cp.asarray(y)`` call site
        this replaces where the source is a live Torch tensor rather than
        an init-time-uploaded constant.
        """
        if self._useCUDA:
            self._items.append(self._cp.asarray(x) if self._useTorch else x)
        else:
            self._items.append(x.data)
        return self

    def img(self, obj):
        """An image/texture object, or any resource the caller has already
        resolved to its final backend-native kernel-argument form (an
        OpenCL ``cl.Image``, an AF-buffer-wrapped ``cl.MemoryObject`` from
        a prior ``x.raw_ptr()``, a CuPy ``TextureObject``, an already
        ``cp.asarray``-viewed Torch tensor, ...): appended as-is on every
        backend, no further unwrapping."""
        self._items.append(obj)
        return self

    # -- scalars ------------------------------------------------------------
    def f32(self, v):
        self._items.append(self._cp.float32(v) if self._useCUDA else self._cl.cltypes.float(v))
        return self

    def u32(self, v):
        self._items.append(self._cp.uint32(v) if self._useCUDA else self._cl.cltypes.uint(v))
        return self

    def i32(self, v):
        self._items.append(self._cp.int32(v) if self._useCUDA else self._cl.cltypes.int(v))
        return self

    def i64(self, v):
        self._items.append(self._cp.int64(v) if self._useCUDA else self._cl.cltypes.long(v))
        return self

    def u64(self, v):
        self._items.append(self._cp.uint64(v) if self._useCUDA else self._cl.cltypes.ulong(v))
        return self

    def u8(self, v):
        self._items.append(self._cp.uint8(v) if self._useCUDA else self._cl.cltypes.uchar(v))
        return self

    # -- vector-typed scalars -------------------------------------------------
    # OpenCL: one cl.cltypes.make_* argument. CUDA (-DPYTHON): the kernel
    # signature instead takes the components as separate scalars, so these
    # append 2-3 items there.
    def vec3f(self, a, b, c):
        if self._useCUDA:
            self._items.append(self._cp.float32(a))
            self._items.append(self._cp.float32(b))
            self._items.append(self._cp.float32(c))
        else:
            self._items.append(self._cl.cltypes.make_float3(a, b, c))
        return self

    def vec3u(self, a, b, c):
        if self._useCUDA:
            self._items.append(self._cp.uint32(a))
            self._items.append(self._cp.uint32(b))
            self._items.append(self._cp.uint32(c))
        else:
            self._items.append(self._cl.cltypes.make_uint3(a, b, c))
        return self

    def vec3i(self, a, b, c):
        if self._useCUDA:
            self._items.append(self._cp.int32(a))
            self._items.append(self._cp.int32(b))
            self._items.append(self._cp.int32(c))
        else:
            self._items.append(self._cl.cltypes.make_int3(a, b, c))
        return self

    def vec2f(self, a, b):
        if self._useCUDA:
            self._items.append(self._cp.float32(a))
            self._items.append(self._cp.float32(b))
        else:
            self._items.append(self._cl.cltypes.make_float2(a, b))
        return self

    # -- backend handoff ------------------------------------------------------
    def seed(self, prefix):
        """Prepend an already-built CuPy prefix tuple (self.kIndF/self.kIndB,
        built once at init time by _build_kIndF/_build_kIndB) before adding
        the per-call arguments. OpenCL/AF instead resume ``set_arg`` from
        the int index _build_kIndF/_build_kIndB already applied at init
        time (passed as ``start_index`` to ``apply_opencl``), so they never
        need this."""
        self._items.extend(prefix)
        return self

    def apply_opencl(self, knl, start_index):
        """set_arg every collected item onto `knl` starting at `start_index`;
        returns the next free index (mirrors the existing kIndLoc += 1 chains)."""
        idx = start_index
        for item in self._items:
            knl.set_arg(idx, item)
            idx += 1
        return idx

    def as_tuple(self):
        """The collected items as a plain tuple, for a CuPy kernel call."""
        return tuple(self._items)


def _mask_fp_resource(self, subset):
    if self.SPECT and self.maskFPZ == self.nHeads:
        return self.d_maskFP
    if self.maskFPZ > 1:
        return self.d_maskFP[subset]
    return self.d_maskFP


def _geometry_buffer(buffers, timestep, subset):
    """Select the per-subset d_x/d_z geometry buffer, falling back to the
    shared index-0 buffer whenever init.py (_initialize_coordinate_buffers)
    did not build a separate entry for this subset. This is the single
    source of truth for that choice -- callers must not re-derive the
    CT/PET/SPECT/listmode condition by hand, since that has repeatedly
    drifted out of sync with the buffer-population logic in init.py."""
    per = buffers[timestep]
    return per[subset] if per[subset] is not None else per[0]


def _load_listmode_subset(self, timestep, subset):
    """List-mode data without loadTOF: only one subset at a time is kept on the device (the entry [timestep][0] of d_x and
    d_TOFIndex, the entry [0] of d_trIndex/d_axIndex). Upload the events of the (timestep, subset) that is about to be
    projected. With loadTOF all subsets were uploaded in initProjector and nothing is done here. Used by both
    forwardProjection and backwardProjection (the backprojection needs the same events as the forward projection)."""
    if self.loadTOF or self.listmode <= 0:
        return
    if self.useCUDA:
        import cupy as cp
        upload = cp.asarray
    else:
        import pyopencl as cl
        import pyopencl.array
        def upload(value):
            return cl.array.to_device(self.queue, value)
    q = timestep * self.subsets + subset
    if self.useIndexBasedReconstruction:
        self.d_trIndex[0] = upload(self.trIndex[self.nMeas[q] * 2 : self.nMeas[q + 1] * 2])
        self.d_axIndex[0] = upload(self.axIndex[self.nMeas[q] * 2 : self.nMeas[q + 1] * 2])
    else:
        apu = self.x.ravel()
        self.d_x[timestep][0] = upload(apu[self.nMeas[q] * 6 : self.nMeas[q + 1] * 6])
    if self.TOF:
        from omegatomo.projector.init import _tof_index_slice
        self.d_TOFIndex[timestep][0] = upload(_tof_index_slice(self, timestep, subset))


def _tof_output_bins(self):
    """Number of TOF bins the FP output measurement vector must be widened by.

    The non-listmode TOF kernels (projectorType123.cl, projectorType4.cl) write
    NBINS values per LOR at idx + to * m_size (NBINS == self.TOF_bins_used, the
    -DNBINS the kernel was compiled with), so the output buffer must have
    m_size * TOF_bins_used elements or the kernel writes/reads out of bounds.
    Listmode TOF instead writes a single value per event (only the TOFid-selected
    bin), so m_size alone is already correct there -- no extra factor.
    """
    return int(self.TOF_bins_used) if (self.TOF and self.listmode == 0) else 1


def _expected_measurement_length(self, timestep, subset):
    """The number of elements the FP output / BP input measurement vector must
    have for this (timestep, subset): the base per-projection size (subsetType
    > 7 or subsets == 1 uses the full nRowsD*nColsD*nProjSubset image, other
    subset types use the raw nMeasSubset LOR count) times the TOF bin factor
    from _tof_output_bins. Single source of truth for both the FP allocation
    size and the BP input-length validation, so they cannot drift apart."""
    if self.listmode == 0 and (self.subsetType > 7 or self.subsets == 1):
        base = int(self.nRowsD) * int(self.nColsD) * int(self.nProjSubset[timestep, subset].item())
    else:
        base = int(self.nMeasSubset[timestep, subset].item())
    return base * _tof_output_bins(self)


def _element_count(x):
    """Cheap host-side element count for an AF/torch/CuPy/PyOpenCL array --
    no device sync, no host copy."""
    if hasattr(x, 'elements'):
        return int(x.elements())
    if hasattr(x, 'numel'):
        return int(x.numel())
    return int(x.size)


def _validate_forward_input(self, f):
    """Validate the FP input image element count(s) against N[k] before any
    device work happens, so a caller-side size mistake raises a clear
    ValueError instead of letting a kernel read/write out of bounds."""
    N = np.asarray(self.N).reshape(-1)
    if isinstance(f, list):
        volume_count = int(self.nMultiVolumes) + 1
        if len(f) != volume_count:
            raise ValueError(f'Expected {volume_count} volume inputs, got {len(f)}')
        for k, image in enumerate(f):
            count = _element_count(image)
            expected = int(N[k])
            if count != expected:
                raise ValueError(f'Volume {k} has {count} elements; expected {expected}')
    else:
        count = _element_count(f)
        expected = int(N[0])
        if count != expected:
            raise ValueError(f'Input image has {count} elements; expected {expected}')


def _validate_backward_input(self, y, timestep, subset):
    """Validate the BP input measurement length before any device work
    happens, so a caller-side size mistake raises a clear ValueError instead
    of letting a kernel read/write out of bounds."""
    count = _element_count(y)
    expected = _expected_measurement_length(self, timestep, subset)
    if count != expected:
        raise ValueError(f'Backprojection input has {count} elements; expected {expected}')


def _append_fp123_args(self, args, timestep, subset, k, f_arg, y_arg):
    """FP types 1-3 per-call kernel-argument sequence, shared by every
    backend (replaces the parallel CuPy-tuple / OpenCL-and-AF-set_arg
    chains). Must be called with `args` already seeded with the caller's
    LOCAL kIndLoc on CuPy (args.seed(kIndLoc)) -- not self.kIndF directly:
    kIndLoc already carries self.kIndF plus the FPType-1-4-shared non-CT
    attenuation arg the caller conditionally appends before dispatching on
    FPType. On OpenCL/AF that same prefix was already applied to the
    kernel via set_arg (by _build_kIndF at init time, plus the caller's own
    matching conditional set_arg immediately before this call), so `args`
    starts fresh there and the caller resumes from its local kIndLoc (the
    int index) via apply_opencl.

    `f_arg`/`y_arg` are the input image / output measurement kernel
    argument, ALREADY resolved to their final backend-native form by the
    caller (a resolved cl.Image/cl.Buffer/MemoryObject for OpenCL/AF, or
    the already torch-or-not resolved CuPy array/texture for CuPy) --
    exactly as before, this function only appends them, it never creates
    or unwraps them itself."""
    if self.useMaskFP:
        args.img(_mask_fp_resource(self, subset))
    if (self.CT or self.PET or self.SPECT) and self.listmode == 0:
        args.i64(self.nProjSubset[timestep, subset].item())
    args.buf(_geometry_buffer(self.d_x, timestep, subset))
    args.buf(_geometry_buffer(self.d_z, timestep, subset))
    if self.normalization_correction:
        args.buf(self.d_norm[timestep][subset])
    if self.additionalCorrection:
        args.buf(self.d_corr[timestep][subset])
    args.buf(self.d_Sens)
    args.vec3u(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item())
    args.vec3f(self.dx[k].item(), self.dy[k].item(), self.dz[k].item())
    args.vec3f(self.bx[k].item(), self.by[k].item(), self.bz[k].item())
    args.vec3f(self.bx[k].item() + self.Nx[k].item() * self.dx[k].item(),
               self.by[k].item() + self.Ny[k].item() * self.dy[k].item(),
               self.bz[k].item() + self.Nz[k].item() * self.dz[k].item())
    if (self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0:
        args.buf(self.d_xyindex[subset])
        args.buf(self.d_zindex[subset])
    if self.useIndexBasedReconstruction and self.listmode > 0:
        if not self.loadTOF:
            args.buf(self.d_trIndex[0])
            args.buf(self.d_axIndex[0])
        else:
            args.buf(self.d_trIndex[subset])
            args.buf(self.d_axIndex[subset])
    if self.listmode > 0 and self.TOF:
        args.buf(_geometry_buffer(self.d_TOFIndex, timestep, subset))
    args.img(f_arg)
    args.img(y_arg)
    if self.SPECT:
        args.buf(self.d_detectorVector[timestep][subset])
    args.u8(self.no_norm)
    args.u64(self.nMeasSubset[timestep, subset].item())
    args.u32(subset)
    args.i32(k)
    if (getattr(self, 'imageAttenuationIsMultiResolution', False)
            and self.attenuation_correction and self.CTAttenuation):
        grid = 0 if int(k) == 0 else 1
        dims = self.imageAttenuationGridDims[grid]
        spacing = self.imageAttenuationGridSpacing[grid]
        origin = self.imageAttenuationGridOrigin[grid]
        args.u32(dims[0]).u32(dims[1]).u32(dims[2])
        args.f32(spacing[0]).f32(spacing[1]).f32(spacing[2])
        args.f32(origin[0]).f32(origin[1]).f32(origin[2])
    return args


def _append_fp5_args(self, args, timestep, subset, k, im_arg, imInt_arg, y_arg, meanV_arg=None):
    """FP type 5 per-call kernel-argument TAIL, shared by every backend.

    Does NOT include the Nx/Ny/Nz/bx/by/bz/dSizeX/dSizeY/dx/dy/dz/
    dScaleX/Y/Z prefix (shared, unchanged, with FPType 4 by a single
    caller-side block run before dispatch -- see _append_fp4_args), and
    does NOT build im_arg/imInt_arg/y_arg itself: those are the summed-area-
    table (af.sat/cumsum/NumPy-cumsum) image textures, built by backend-
    specific code -- only ALREADY-resolved values are appended here.

    Unlike BP type 5's meanBP (a single trailing arg that the caller can
    tack on after this function returns, since it is the LAST kernel
    argument -- see _append_bp5_ct_args), meanFP's d_meanV sits BETWEEN
    d_forw (y_arg) and the optional maskFP in the kernel signature
    (projectorType5.cl's projectorType5Forward: ... d_forw BUF5,
    #ifdef MEANDISTANCEFP d_meanV BUF6, #ifdef MASKFP maskFP TEX7, ...),
    so it cannot be appended after this call -- it must be threaded
    through as a parameter instead. `meanV_arg` is the already-built,
    already-resolved mean-value buffer/array (only read when
    self.meanFP is set); callers not using meanFP omit it."""
    args.buf(_geometry_buffer(self.d_x, timestep, subset))
    args.buf(_geometry_buffer(self.d_z, timestep, subset))
    args.img(im_arg)
    args.img(imInt_arg)
    args.img(y_arg)
    if self.meanFP:
        args.img(meanV_arg)
    if self.useMaskFP:
        args.img(_mask_fp_resource(self, subset))
    args.i64(self.nProjSubset[timestep, subset].item())
    return args


def _fp5_integral_images_numpy(image_host, Nx, Ny, Nz, meanFP):
    """Host (NumPy) computation of the FP type 5 XZ-/YZ-plane summed-area-
    table (integral) images, for the plain-PyOpenCL backend (which has no
    ArrayFire af.sat/af.mean to build these on-device). Mirrors the C++
    reference (functions.hpp updateInputs()) and the CuPy branch of
    forwardProjection() exactly:

      * the flat, Fortran-ordered (Nx,Ny,Nz) image is reordered into two
        different (slice, plane-row, plane-col) layouts,
      * each layout optionally has its per-slice mean subtracted
        (meanFP) -- the C++/AF/CuPy analogue of af.mean(af.mean(im,0),1),
      * each is cumulatively summed along its first two axes (equivalent
        to af.sat of the unpadded array), then zero-padded with a
        leading row and column (matching intIm[1:,1:,:] = af.sat(...)).

    Returns (im_XZ, im_YZ, meanFP_host):
      * im_XZ has shape (Nx+1, Nz+1, Ny) -- the kernel's first image
        argument (d_image_os / d_IImageY / "XZ-plane"),
      * im_YZ has shape (Ny+1, Nz+1, Nx) -- the kernel's second image
        argument (d_image_os_int / d_IImageX / "YZ-plane"),
      * meanFP_host is a float32 array of Nx+Ny elements (first Nx =
        im_YZ's per-x-slice means, next Ny = im_XZ's per-y-slice means,
        matching vec.meanFP's layout in functions.hpp), or None when
        meanFP is False."""
    vol = np.asarray(image_host).reshape((Nx, Ny, Nz), order='F').astype(np.float32, copy=False)
    meanFP_host = np.zeros(Nx + Ny, dtype=np.float32) if meanFP else None

    im_yz = np.transpose(vol, (1, 2, 0))  # (Ny, Nz, Nx)
    if meanFP:
        meanFP_host[0:Nx] = im_yz.mean(axis=(0, 1)).astype(np.float32)
        im_yz = im_yz - meanFP_host[0:Nx].reshape((1, 1, Nx))
    im_yz = np.cumsum(im_yz, axis=0)
    im_yz = np.cumsum(im_yz, axis=1)
    padded_yz = np.zeros((Ny + 1, Nz + 1, Nx), dtype=np.float32, order='F')
    padded_yz[1:, 1:, :] = im_yz

    im_xz = np.transpose(vol, (0, 2, 1))  # (Nx, Nz, Ny)
    if meanFP:
        meanFP_host[Nx:Nx + Ny] = im_xz.mean(axis=(0, 1)).astype(np.float32)
        im_xz = im_xz - meanFP_host[Nx:Nx + Ny].reshape((1, 1, Ny))
    im_xz = np.cumsum(im_xz, axis=0)
    im_xz = np.cumsum(im_xz, axis=1)
    padded_xz = np.zeros((Nx + 1, Nz + 1, Ny), dtype=np.float32, order='F')
    padded_xz[1:, 1:, :] = im_xz

    return padded_xz, padded_yz, meanFP_host


def _append_bp5_ct_args(self, args, timestep, subset, k, f_arg, y_arg):
    """BP type 5 (CT) per-call kernel-argument sequence, shared by every
    backend. Does NOT include a trailing meanBP arg: on every backend
    (OpenCL/AF and CuPy) the caller appends the d_meanBP/dMeanBP argument
    itself, right after this call, guarded by 'if self.meanBP:' -- kept as
    a few lines of caller code rather than folded into this function.

    `y_arg` is the already-built summed-area-table (af.sat/cumsum) image
    texture -- built by code that must stay byte-identical and cannot be
    exercised by the AF variant of this repo's test harness -- appended
    here as-is exactly like before."""
    args.buf(_geometry_buffer(self.d_x, timestep, subset))
    args.buf(_geometry_buffer(self.d_z, timestep, subset))
    if self.listmode == 0:
        args.buf(self.d_geom5[timestep][subset])
    args.img(y_arg)
    args.img(f_arg)
    args.buf(self.d_Sens)
    return args


def _append_listmode_index_args(self, args, timestep, subset):
    """The list-mode arguments of the type 4 kernels that follow the optional subset index arguments (xyindex/zindex):
    the index-based detector indices (INDEXBASED && LISTMODE: trIndex, axIndex) and the TOF bin of each event
    (LISTMODE && TOF: TOFIndex), in the order of projectorType4.cl. (Types 1-3 append the same arguments inline.)"""
    if self.useIndexBasedReconstruction and self.listmode > 0:
        if not self.loadTOF:
            args.buf(self.d_trIndex[0])
            args.buf(self.d_axIndex[0])
        else:
            args.buf(self.d_trIndex[subset])
            args.buf(self.d_axIndex[subset])
    if self.listmode > 0 and self.TOF:
        args.buf(_geometry_buffer(self.d_TOFIndex, timestep, subset))
    return args


def _append_fp4_args(self, args, timestep, subset, k, f_arg, y_arg):
    """FP type 4 per-call kernel-argument TAIL, shared by every backend.

    Does NOT include the Nx/Ny/Nz/bx/by/bz/bmax/dScaleX4Y4Z4 prefix: that
    segment is shared, unchanged, between FPType 4 and FPType 5 by a single
    'if self.FPType == 5 or self.FPType == 4:' block the caller runs before
    dispatching on FPType (left untouched here, since unifying it would
    also have to touch the FPType-5-only ArrayFire code the task says must
    stay byte-identical and cannot be tested). `args` must be seeded from
    the caller's local kIndLoc (which already carries that prefix, plus
    self.kIndF and the FPType-1-4-shared non-CT attenuation arg) on CuPy;
    on OpenCL/AF the caller has already applied that same prefix via
    set_arg and resumes from its local kIndLoc (the int index).
    """
    args.img(f_arg)
    args.img(y_arg)
    args.buf(_geometry_buffer(self.d_x, timestep, subset))
    args.buf(_geometry_buffer(self.d_z, timestep, subset))
    if self.useMaskFP:
        args.img(_mask_fp_resource(self, subset))
    args.i64(self.nProjSubset[timestep, subset].item())
    if (self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0:
        args.buf(self.d_xyindex[subset])
        args.buf(self.d_zindex[subset])
    _append_listmode_index_args(self, args, timestep, subset)
    if self.normalization_correction:
        args.buf(self.d_norm[timestep][subset])
    if self.additionalCorrection:
        args.buf(self.d_corr[timestep][subset])
    args.u8(self.no_norm)
    args.u64(self.nMeasSubset[timestep, subset].item())
    args.u32(subset)
    args.i32(k)
    return args


def _append_bp123_args(self, args, timestep, subset, k, f_arg, y_arg, proj_count_arg):
    """BP types 1-3 per-call kernel-argument sequence, shared by every
    backend (replaces the parallel CuPy-tuple / OpenCL-and-AF-set_arg
    chains). Unlike FP1-3, the whole sequence (including the non-CT
    attenuation arg) lives inside one shared 'if self.BPType in [1,2,3]:'
    block on both backends, so `args` can be seeded straight from
    self.kIndB (via args.seed(kIndLoc), kIndLoc == self.kIndB unchanged at
    the call site).

    `f_arg`/`y_arg` are the output image / input measurement kernel
    argument, already resolved to their final backend-native form by the
    caller exactly as before. `proj_count_arg` is the
    "(CT or PET or SPECT) and listmode==0" projection-count argument,
    ALSO pre-resolved by the caller: OpenCL/AF wrap it as cl.cltypes.long
    and CuPy wraps it as cp.int64 -- both match the kernel's
    `const LONG d_nProjections` (projectorType123.cl)."""
    if (getattr(self, 'imageAttenuationIsMultiResolution', False)
            and self.attenuation_correction and self.CTAttenuation):
        attenuation = self.d_imageAttenuation[timestep][k]
        if self.useImages:
            args.img(attenuation)
        else:
            args.buf(attenuation)
    elif self.attenuation_correction and not self.CTAttenuation:
        args.buf(self.d_atten[timestep][subset])
    if self.useMaskFP:
        args.img(_mask_fp_resource(self, subset))
    if self.useMaskBP:
        args.img(self.d_maskBP)
    if (self.CT or self.PET or self.SPECT) and self.listmode == 0:
        args.img(proj_count_arg)
    args.buf(_geometry_buffer(self.d_x, timestep, subset))
    args.buf(_geometry_buffer(self.d_z, timestep, subset))
    if self.normalization_correction:
        args.buf(self.d_norm[timestep][subset])
    if self.additionalCorrection:
        args.buf(self.d_corr[timestep][subset])
    args.buf(self.d_Sens)
    args.vec3u(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item())
    args.vec3f(self.dx[k].item(), self.dy[k].item(), self.dz[k].item())
    args.vec3f(self.bx[k].item(), self.by[k].item(), self.bz[k].item())
    args.vec3f(self.bx[k].item() + self.Nx[k].item() * self.dx[k].item(),
               self.by[k].item() + self.Ny[k].item() * self.dy[k].item(),
               self.bz[k].item() + self.Nz[k].item() * self.dz[k].item())
    if (self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0:
        args.buf(self.d_xyindex[subset])
        args.buf(self.d_zindex[subset])
    if self.useIndexBasedReconstruction and self.listmode > 0:
        if not self.loadTOF:
            args.buf(self.d_trIndex[0])
            args.buf(self.d_axIndex[0])
        else:
            args.buf(self.d_trIndex[subset])
            args.buf(self.d_axIndex[subset])
    if self.listmode > 0 and self.TOF:
        args.buf(_geometry_buffer(self.d_TOFIndex, timestep, subset))
    args.img(y_arg)
    args.img(f_arg)
    if self.SPECT:
        args.buf(self.d_detectorVector[timestep][subset])
    args.u8(self.no_norm)
    args.u64(self.nMeasSubset[timestep, subset].item())
    args.u32(subset)
    args.i32(k)
    if (getattr(self, 'imageAttenuationIsMultiResolution', False)
            and self.attenuation_correction and self.CTAttenuation):
        grid = 0 if int(k) == 0 else 1
        dims = self.imageAttenuationGridDims[grid]
        spacing = self.imageAttenuationGridSpacing[grid]
        origin = self.imageAttenuationGridOrigin[grid]
        args.u32(dims[0]).u32(dims[1]).u32(dims[2])
        args.f32(spacing[0]).f32(spacing[1]).f32(spacing[2])
        args.f32(origin[0]).f32(origin[1]).f32(origin[2])
    return args


def _append_bp4_ct_args(self, args, timestep, subset, k, f_arg, y_arg):
    """BP type 4 (CT) per-call kernel-argument MIDDLE segment, shared by
    every backend. Deliberately tiny: the OffsetLimit/Nx/Ny/Nz/bx/by/bz/
    dx/dy/dz/kerroin prefix is shared, unchanged, between BPType 4 and
    BPType 5 by code the caller runs before dispatching on BPType (left
    untouched, since unifying it would also have to touch the BPType-5-only
    ArrayFire code that must stay byte-identical and cannot be tested), and
    the no_norm/mask/projection-count/k tail is shared between BPType 4-CT,
    5-CT AND the non-CT branch by more caller code that runs after (also
    left untouched, for the same reason). `args` is seeded from -- and its
    result read back into -- the caller's local kIndLoc, which already
    carries everything before and will go on to carry everything after.

    `f_arg`/`y_arg` are the output image / input measurement kernel
    argument, already resolved to their final backend-native form by the
    caller exactly as before."""
    args.img(y_arg)
    args.img(f_arg)
    args.buf(_geometry_buffer(self.d_x, timestep, subset))
    args.buf(_geometry_buffer(self.d_z, timestep, subset))
    args.buf(self.d_Sens)
    return args


def _append_bp4_nonct_args(self, args, timestep, subset, k, f_arg, y_arg, proj_count_arg):
    """BP type 4 (non-CT) per-call kernel-argument sequence (minus the
    shared no_norm/nMeasSubset/subset/k tail -- see _append_bp4_ct_args'
    docstring; the same shared-tail caller code follows both), shared by
    every backend.

    `proj_count_arg` is the projection-count argument, pre-resolved by the
    caller: the pre-existing CuPy source wraps it as cp.int64 while
    OpenCL/AF wrap the SAME value as cl.cltypes.ulong (a signed/unsigned
    mismatch -- preserved verbatim rather than unified, since this file
    must not change kernel-argument types)."""
    args.vec3u(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item())
    args.vec3f(self.bx[k].item(), self.by[k].item(), self.bz[k].item())
    args.vec3f(self.bx[k].item() + self.Nx[k].item() * self.dx[k].item(),
               self.by[k].item() + self.Ny[k].item() * self.dy[k].item(),
               self.bz[k].item() + self.Nz[k].item() * self.dz[k].item())
    args.vec3f(self.dScaleX4[k].item(), self.dScaleY4[k].item(), self.dScaleZ4[k].item())
    args.img(y_arg)
    args.img(f_arg)
    args.buf(_geometry_buffer(self.d_x, timestep, subset))
    args.buf(_geometry_buffer(self.d_z, timestep, subset))
    if self.useMaskFP:
        args.img(_mask_fp_resource(self, subset))
    if self.useMaskBP:
        args.img(self.d_maskBP)
    args.img(proj_count_arg)
    if (self.subsetType == 3 or self.subsetType == 6 or self.subsetType == 7) and self.subsets > 1 and self.listmode == 0:
        args.buf(self.d_xyindex[subset])
        args.buf(self.d_zindex[subset])
    _append_listmode_index_args(self, args, timestep, subset)
    if self.normalization_correction:
        args.buf(self.d_norm[timestep][subset])
    if self.additionalCorrection:
        args.buf(self.d_corr[timestep][subset])
    args.buf(self.d_Sens)
    return args


def _applyOffsetSqrtFP(self, y, subset):
    """Weights a forward projection (or measurement data) y of the given subset with the square root of the offset redundancy weights,
    sqrt(W), see offsetSqrtWeights in detcoord.py. The backprojection kernels apply the same weights on the other side (-DOFFSETSQRT), so
    the projector as a whole is the linear sqrt(W) * A. The weights are applied once to the final output, after all the volumes have been
    projected. Returns y as is unless the square root offset weighting is in use; otherwise returns the weighted array (y itself is not modified)."""
    if not getattr(self, '_offsetSqrtActive', False):
        return y
    w = self.d_offsetSqrtW[subset]
    nR = int(self.nRowsD)
    nC = int(self.nColsD)
    # The data are ordered with the detector row as the fastest index, then the column and then the projection
    if self.useMetal or (self.useCUDA and self.useTorch):
        nP = w.numel() // nR
        return (y.reshape(nP, nC, nR) * w.reshape(nP, 1, nR)).reshape(y.shape)
    elif self.useCUDA:
        nP = w.size // nR
        return (y.reshape(nP, nC, nR) * w.reshape(nP, 1, nR)).reshape(y.shape)
    elif self.useAF:
        import arrayfire as af
        nP = w.elements() // nR
        return af.flat(af.moddims(y, nR, nC, nP) * af.tile(af.moddims(w, nR, 1, nP), 1, nC, 1))
    else:
        nP = w.size // nR
        yW = y.copy()
        self._knlOffsetSqrt(yW, w, np.uint32(nR), np.uint32(nR * nC))
        return yW


def conv3D(self, f, ii = 0):
    if self.useTorch:
        import torch
        import torch.nn.functional as functional
        nx, ny, nz = (int(self.Nx[ii]), int(self.Ny[ii]), int(self.Nz[ii]))
        rx, ry, rz = (int(self.g_dim_x), int(self.g_dim_y), int(self.g_dim_z))
        image = f.reshape(1, 1, nz, ny, nx)
        # The native convolution reflects at the voxel boundary, including the edge voxel.
        for axis, size, radius in ((2, nz, rz), (3, ny, ry), (4, nx, rx)):
            if radius:
                indices = torch.arange(-radius, size + radius, device=f.device) % (2 * size)
                indices = torch.where(indices < size, indices, 2 * size - indices - 1).long()
                image = image.index_select(axis, indices)
        kernel = np.asarray(self.gaussK, dtype=np.float32).ravel(order='F').copy()
        kernel = torch.as_tensor(kernel, device=f.device).reshape(1, 1, 2 * rz + 1, 2 * ry + 1, 2 * rx + 1)
        return functional.conv3d(image, kernel).reshape(-1)
    globalSize = (self.Nx[ii].item() + self.erotusBP[ii * 2], self.Ny[ii].item() + self.erotusBP[ii * 2 + 1], self.Nz[ii].item())
    if self.useCUDA:
        if self.useTorch:
            import torch
        if self.useCuPy:
            import cupy as cp
        else:
            raise ValueError('PyCUDA is no longer supported. Please use CuPy.')
    else:
        import pyopencl as cl
        import pyopencl.array
    args = _KernelArgs(self)
    if self.useAF:
        import arrayfire as af
        if isinstance(f, af.array.Array):
            ptr = f.raw_ptr()
            f = cl.MemoryObject.from_int_ptr(ptr)
            args.img(f)
        else:
            args.buf(f)
    else:
        if self.useCUDA:
            if self.useTorch:
                output = torch.zeros(self.N[ii].item(), dtype=torch.float32, device='cuda')
            else:
                output = cp.zeros(self.N[ii].item(), dtype=cp.float32)
        else:
            output = cl.array.zeros(self.queue, self.N[ii].item(), dtype=cl.cltypes.float)
        if not self.useCUDA:
            args.buf(f)
    if self.useCUDA:
        # Convolution3D_f (auxKernels.cl) now splits its trailing int3 N (image dimensions)
        # into three separate ints under -DPYTHON, matching every other vector-typed kernel
        # argument in this codebase (see RDPKernel/TVKernel/NLM for the same pattern).
        if self.useTorch:
            fD = cp.asarray(f)
            outputD = cp.asarray(output)
        else:
            fD = f
            outputD = output
        args.img(fD).img(outputD).img(self.d_gaussPSF)
        args.i32(self.g_dim_x).i32(self.g_dim_y).i32(self.g_dim_z)
        args.vec3i(self.Nx[ii].item(), self.Ny[ii].item(), self.Nz[ii].item())
        self.knlPSF((globalSize[0] // 16, globalSize[1] // 16, globalSize[2] // 1), (16, 16, 1), args.as_tuple())
        if self.useTorch:
            torch.cuda.synchronize()
    else:
        if self.useAF:
            output = af.data.constant(0., self.N[ii].item(), dtype=af.Dtype.f32)
            outPtr = cl.MemoryObject.from_int_ptr(output.raw_ptr())
            args.img(outPtr)
        else:
            args.buf(output)
        args.buf(self.d_gaussPSF)
        args.i32(self.g_dim_x).i32(self.g_dim_y).i32(self.g_dim_z)
        args.vec3i(self.Nx[ii].item(), self.Ny[ii].item(), self.Nz[ii].item())
        args.apply_opencl(self.knlPSF, 0)
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
    if isinstance(f, list):
        # Work on a shallow copy so PSF convolution below (which replaces f[k] in place)
        # never mutates the caller's own list.
        f = list(f)
    if self.useMetal:
        from omegatomo.projector.mps_backend import forward_projection_mps
        return _applyOffsetSqrtFP(self, forward_projection_mps(self, f, subset, timestep), subset)
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
                inputs[volume] = self.computeConvolution(image, volume) if self.use_psf else image
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
        _validate_forward_input(self, f)
        if self.nMultiVolumes > 0 and not(isinstance(f,list)):
            volumes = self.nMultiVolumes
            self.nMultiVolumes = 0
        if self.useCUDA:
            if self.useCuPy:
                import cupy as cp
                _load_listmode_subset(self, timestep, subset)
                measLen = _expected_measurement_length(self, timestep, subset)
                if self.useTorch:
                    import torch
                    y = torch.zeros(measLen, dtype=torch.float32, device='cuda')
                    yD = cp.asarray(y)
                else:
                    y = cp.zeros(measLen, dtype=cp.float32)
                for k in range(self.nMultiVolumes + 1):
                    if isinstance(f,list):
                        if self.use_psf:
                            f[k] = self.computeConvolution(f[k], k)
                        if self.useTorch:
                            fD = cp.asarray(f[k])
                    else:
                        if self.use_psf:
                            f = self.computeConvolution(f)
                        if self.useTorch:
                            fD = cp.asarray(f)
                    if self.FPType == 5:
                        if self.useTorch:
                            volSrc = fD
                        elif isinstance(f,list):
                            volSrc = f[k]
                        else:
                            volSrc = f
                        vol5 = volSrc.reshape((self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()), order='F')
                        dMeanFP = cp.zeros(self.Nx[k].item() + self.Ny[k].item(), dtype=cp.float32) if self.meanFP else None
                        intIm = cp.zeros((self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item()), dtype=cp.float32, order='F')
                        im_yz = cp.transpose(vol5, (1, 2, 0))
                        if self.meanFP:
                            dMeanFP[0:self.Nx[k].item()] = cp.mean(im_yz, axis=(0, 1)).astype(cp.float32)
                            im_yz = im_yz - dMeanFP[0:self.Nx[k].item()].reshape((1, 1, self.Nx[k].item()))
                        intIm[1:,1:,:] = im_yz
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
                        im_xz = cp.transpose(vol5, (0, 2, 1))
                        if self.meanFP:
                            dMeanFP[self.Nx[k].item():self.Nx[k].item() + self.Ny[k].item()] = cp.mean(im_xz, axis=(0, 1)).astype(cp.float32)
                            im_xz = im_xz - dMeanFP[self.Nx[k].item():self.Nx[k].item() + self.Ny[k].item()].reshape((1, 1, self.Ny[k].item()))
                        intIm[1:,1:,:] = im_xz
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
                        if (getattr(self, 'imageAttenuationIsMultiResolution', False)
                                and self.attenuation_correction and self.CTAttenuation):
                            kIndLoc += (self.d_imageAttenuation[timestep][k],)
                        elif (self.attenuation_correction and not self.CTAttenuation):
                            kIndLoc += (self.d_atten[timestep][subset],)
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
                        f_arg = ff
                        y_arg = yD if self.useTorch else y
                        args = _KernelArgs(self).seed(kIndLoc)
                        _append_fp4_args(self, args, timestep, subset, k, f_arg, y_arg)
                        kIndLoc = args.as_tuple()
                    elif self.FPType == 5:
                        y_arg = yD if self.useTorch else y
                        args = _KernelArgs(self).seed(kIndLoc)
                        _append_fp5_args(self, args, timestep, subset, k, ff, ff2, y_arg, meanV_arg=dMeanFP if self.meanFP else None)
                        kIndLoc = args.as_tuple()
                    elif self.FPType in [1, 2, 3]:
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
                                f_arg = ff
                            else:
                                f_arg = fD if self.useTorch else f[k]
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
                                f_arg = ff
                            else:
                                f_arg = fD if self.useTorch else f
                        y_arg = yD if self.useTorch else y
                        # Seed from the LOCAL kIndLoc, not self.kIndF directly: kIndLoc
                        # already carries self.kIndF plus the FPType-1-4-shared
                        # non-CT-attenuation arg conditionally appended above.
                        args = _KernelArgs(self).seed(kIndLoc)
                        _append_fp123_args(self, args, timestep, subset, k, f_arg, y_arg)
                        kIndLoc = args.as_tuple()
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
            import pyopencl.array
            from pyopencl.version import VERSION
            _load_listmode_subset(self, timestep, subset)
            measLen = _expected_measurement_length(self, timestep, subset)
            if self.useAF:
                import arrayfire as af
                y = af.data.constant(0., measLen)
                yPtr = y.raw_ptr()
                yD = cl.MemoryObject.from_int_ptr(yPtr)
            else:
                y = cl.array.zeros(self.queue, measLen, dtype=cl.cltypes.float)
            imformat = cl.ImageFormat(cl.channel_order.A, cl.channel_type.FLOAT)
            mf = cl.mem_flags
            for k in range(self.nMultiVolumes + 1):
                # FP type 5's -DMEANDISTANCEFP d_meanV kernel argument (only used when
                # self.meanFP is set); built below alongside the SAT images and consumed
                # right after the FPType == 5 _append_fp5_args() call further down.
                meanV_arg = None
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
                            f[k] = self.computeConvolution(f[k], k)
                        if self.useAF:
                            if self.FPType < 5:
                                fPtr = f[k].raw_ptr()
                                fD = cl.MemoryObject.from_int_ptr(fPtr)
                                cl.enqueue_copy(self.queue, d_im, fD, offset=(0), origin=(0,0,0), region=(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()));
                                af.device.unlock_array(f[k])
                            else:
                                d_meanFP = af.data.constant(0., self.Nx[k].item() + self.Ny[k].item()) if self.meanFP else None
                                intIm = af.data.constant(0., self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item())
                                if self.meanFP:
                                    im = af.reorder(af.moddims(f[k], self.Nx[k].item(), d1=self.Ny[k].item(), d2=self.Nz[k].item()), d0=1, d1=2, d2=0)
                                    mFPxy = af.mean(af.mean(im, dim=0), dim=1)
                                    d_meanFP[0:self.Nx[k].item()] = af.flat(mFPxy)
                                    im -= af.tile(mFPxy, d0=im.shape[0], d1=im.shape[1])
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
                                # Second (XZ-plane) integral image: mirrors functions.hpp
                                # updateInputs()'s second af::sat block exactly, including the
                                # meanFP subtraction using d_meanFP's LAST Ny entries (the first
                                # image above used the first Nx entries) -- previously missing
                                # here entirely, so meanFP only ever touched the first image.
                                if self.meanFP:
                                    im = af.reorder(af.moddims(f[k], self.Nx[k].item(), d1=self.Ny[k].item(), d2=self.Nz[k].item()), d0=0, d1=2, d2=1)
                                    mFPyz = af.mean(af.mean(im, dim=0), dim=1)
                                    d_meanFP[self.Nx[k].item():self.Nx[k].item() + self.Ny[k].item()] = af.flat(mFPyz)
                                    im -= af.tile(mFPyz, d0=im.shape[0], d1=im.shape[1])
                                    intIm[1:,1:,:] = af.sat(im)
                                    af.eval(im)
                                else:
                                    intIm[1:,1:,:] = af.sat(af.reorder(af.moddims(f[k], self.Nx[k].item(), d1=self.Ny[k].item(), d2=self.Nz[k].item()), d0=0, d1=2, d2=1))
                                af.eval(intIm)
                                intIm = af.flat(intIm)
                                fPtr = intIm.raw_ptr()
                                fD = cl.MemoryObject.from_int_ptr(fPtr)
                                cl.enqueue_copy(self.queue, d_im, fD, offset=(0), origin=(0,0,0), region=(self.Nx[k].item() + 1, self.Nz[k].item() + 1, self.Ny[k].item()));
                                af.device.unlock_array(intIm)
                                if self.meanFP:
                                    # d_meanFP stays LOCKED (unlike the SAT images above, which get
                                    # copied into their own cl.Image and can be unlocked right away):
                                    # it is used directly as the kernel's d_meanV argument, so it must
                                    # remain valid through the kernel dispatch and is only unlocked
                                    # after self.queue.finish() below.
                                    meanV_arg = cl.MemoryObject.from_int_ptr(d_meanFP.raw_ptr())
                        else:
                            if self.FPType < 5:
                                cl.enqueue_copy(self.queue, d_im, f[k].data, offset=(0), origin=(0,0,0), region=(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()));
                            else:
                                # Plain PyOpenCL has no ArrayFire af.sat/af.mean, so build both
                                # summed-area-table images (and, when meanFP, the mean array) on
                                # the host with NumPy -- previously this branch did a single
                                # wrong-shape copy into d_im and never touched d_imInt at all,
                                # leaving it uninitialized device memory read by the kernel.
                                im_xz, im_yz, meanFP_host = _fp5_integral_images_numpy(f[k].get(), self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item(), self.meanFP)
                                cl.enqueue_copy(self.queue, d_imInt, im_yz, origin=(0, 0, 0), region=(self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item()));
                                cl.enqueue_copy(self.queue, d_im, im_xz, origin=(0, 0, 0), region=(self.Nx[k].item() + 1, self.Nz[k].item() + 1, self.Ny[k].item()));
                                if self.meanFP:
                                    meanV_arg = cl.Buffer(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=meanFP_host)
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
                                d_meanFP = af.data.constant(0., self.Nx[k].item() + self.Ny[k].item()) if self.meanFP else None
                                intIm = af.data.constant(0., self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item())
                                if self.meanFP:
                                    im = af.reorder(af.moddims(f, self.Nx[k].item(), d1=self.Ny[k].item(), d2=self.Nz[k].item()), d0=1, d1=2, d2=0)
                                    mFPxy = af.mean(af.mean(im, dim=0), dim=1)
                                    d_meanFP[0:self.Nx[k].item()] = af.flat(mFPxy)
                                    im -= af.tile(mFPxy, d0=im.shape[0], d1=im.shape[1])
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
                                if self.meanFP:
                                    im = af.reorder(af.moddims(f, self.Nx[k].item(), d1=self.Ny[k].item(), d2=self.Nz[k].item()), d0=0, d1=2, d2=1)
                                    mFPyz = af.mean(af.mean(im, dim=0), dim=1)
                                    d_meanFP[self.Nx[k].item():self.Nx[k].item() + self.Ny[k].item()] = af.flat(mFPyz)
                                    im -= af.tile(mFPyz, d0=im.shape[0], d1=im.shape[1])
                                    intIm[1:,1:,:] = af.sat(im)
                                    af.eval(im)
                                else:
                                    intIm[1:,1:,:] = af.sat(af.reorder(af.moddims(f, self.Nx[k].item(), d1=self.Ny[k].item(), d2=self.Nz[k].item()), d0=0, d1=2, d2=1))
                                af.eval(intIm)
                                intIm = af.flat(intIm)
                                af.sync()
                                fPtr = intIm.raw_ptr()
                                fD = cl.MemoryObject.from_int_ptr(fPtr)
                                cl.enqueue_copy(self.queue, d_im, fD, offset=(0), origin=(0,0,0), region=(self.Nx[k].item() + 1, self.Nz[k].item() + 1, self.Ny[k].item()));
                                af.device.unlock_array(intIm)
                                if self.meanFP:
                                    meanV_arg = cl.MemoryObject.from_int_ptr(d_meanFP.raw_ptr())
                        else:
                            if self.FPType < 5:
                                cl.enqueue_copy(self.queue, d_im, f.data, offset=(0), origin=(0,0,0), region=(self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item()));
                            else:
                                im_xz, im_yz, meanFP_host = _fp5_integral_images_numpy(f.get(), self.Nx[k].item(), self.Ny[k].item(), self.Nz[k].item(), self.meanFP)
                                cl.enqueue_copy(self.queue, d_imInt, im_yz, origin=(0, 0, 0), region=(self.Ny[k].item() + 1, self.Nz[k].item() + 1, self.Nx[k].item()));
                                cl.enqueue_copy(self.queue, d_im, im_xz, origin=(0, 0, 0), region=(self.Nx[k].item() + 1, self.Nz[k].item() + 1, self.Ny[k].item()));
                                if self.meanFP:
                                    meanV_arg = cl.Buffer(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=meanFP_host)
                else:
                    if self.useAF:
                        if isinstance(f,list):
                            fPtr = f[k].raw_ptr()
                        else:
                            fPtr = f.raw_ptr()
                        d_im = cl.MemoryObject.from_int_ptr(fPtr)
                kIndLoc = self.kIndF
                if self.FPType == 1 or self.FPType == 2 or self.FPType == 3 or self.FPType == 4:
                    if (getattr(self, 'imageAttenuationIsMultiResolution', False)
                            and self.attenuation_correction and self.CTAttenuation):
                        resource = self.d_imageAttenuation[timestep][k]
                        self.knlF.set_arg(kIndLoc, resource if self.useImages else resource.data)
                        kIndLoc += 1
                    elif (self.attenuation_correction and not self.CTAttenuation):
                        self.knlF.set_arg(kIndLoc, self.d_atten[timestep][subset].data)
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
                    y_arg = yD if self.useAF else y.data
                    args = _KernelArgs(self)
                    _append_fp4_args(self, args, timestep, subset, k, d_im, y_arg)
                    kIndLoc = args.apply_opencl(self.knlF, kIndLoc)
                elif self.FPType == 5:
                    y_arg = yD if self.useAF else y.data
                    args = _KernelArgs(self)
                    _append_fp5_args(self, args, timestep, subset, k, d_im, d_imInt, y_arg, meanV_arg=meanV_arg)
                    kIndLoc = args.apply_opencl(self.knlF, kIndLoc)
                elif self.FPType in [1, 2, 3]:
                    if not self.useImages and not self.useAF:
                        f_arg = f[k].data if isinstance(f, list) else f.data
                    else:
                        f_arg = d_im
                    y_arg = yD if self.useAF else y.data
                    args = _KernelArgs(self)
                    _append_fp123_args(self, args, timestep, subset, k, f_arg, y_arg)
                    kIndLoc = args.apply_opencl(self.knlF, kIndLoc)
                cl.enqueue_nd_range_kernel(self.queue, self.knlF, self.globalSizeFP[timestep][subset], self.localSizeFP)
                self.queue.finish()
                if self.useAF and meanV_arg is not None:
                    # d_meanFP (see the FPType == 5 image-building block above) was kept
                    # locked through the kernel dispatch since meanV_arg wraps its raw
                    # device pointer directly (no intermediate cl.Image copy) -- unlock it
                    # now that the kernel has finished reading it.
                    af.device.unlock_array(d_meanFP)
        if volumes > 0 and not(isinstance(f,list)):
            self.nMultiVolumes = volumes
        if self.useAF:
            af.device.unlock_array(y)
            if not self.useImages:
                af.device.unlock_array(f)
    return _applyOffsetSqrtFP(self, y, subset)

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
                outputs.append(self.computeConvolution(output, volume) if self.use_psf else output)
            return outputs[0] if int(self.nMultiVolumes) == 0 else outputs
    else:
        _validate_backward_input(self, y, timestep, subset)
        _load_listmode_subset(self, timestep, subset)
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
                        nProjBP5 = self.nProjSubset[timestep, subset].item()
                        yy = cp.zeros((self.nRowsD+1,self.nColsD+1,nProjBP5), dtype=cp.float32, order='F')
                        if self.useTorch:
                            yUnpadded = yD.reshape((self.nRowsD,self.nColsD,nProjBP5), order='F')
                        else:
                            yUnpadded = y.reshape((self.nRowsD,self.nColsD,nProjBP5), order='F')
                        if self.meanBP:
                            # Mirrors the C++ reference (functions.hpp computeIntegralImage) and the
                            # plain-PyOpenCL NumPy fallback below: per-projection mean over the
                            # detector rows/columns, subtracted before the summed-area-table cumsum.
                            dMeanBP = cp.mean(yUnpadded, axis=(0, 1)).astype(cp.float32)
                            yUnpadded = yUnpadded - dMeanBP.reshape((1, 1, nProjBP5))
                        yy[1:,1:,:] = yUnpadded
                        yy = yy.cumsum(0)
                        yy = yy.cumsum(1)
                        yy = yy.ravel(order='F')
                    kIndLoc = self.kIndB
                    if self.BPType in [1, 2, 3]:
                        y_arg = yD if self.useTorch else y
                        if self.useTorch:
                            f_arg = fD
                        else:
                            f_arg = f[k] if self.nMultiVolumes > 0 else f
                        # Kernel parameter is `const LONG d_nProjections` (projectorType123.cl);
                        # explicit cp.int64 to match every other LONG-typed arg in this launch
                        # tuple -- a bare Python int is not guaranteed to push the same 8-byte
                        # kernel-argument encoding CuPy's RawKernel uses for an explicit int64.
                        proj_count_arg = cp.int64(self.nProjSubset[timestep, subset].item())
                        args = _KernelArgs(self).seed(kIndLoc)
                        _append_bp123_args(self, args, timestep, subset, k, f_arg, y_arg, proj_count_arg)
                        kIndLoc = args.as_tuple()
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
                                    y_arg = yy
                                else:
                                    y_arg = yD if self.useTorch else y
                                if self.useTorch:
                                    f_arg = fD
                                else:
                                    f_arg = f[k] if isinstance(f, list) else f
                                args = _KernelArgs(self).seed(kIndLoc)
                                _append_bp4_ct_args(self, args, timestep, subset, k, f_arg, y_arg)
                                kIndLoc = args.as_tuple()
                            else:
                                if self.useImages:
                                    chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
                                    array = cp.cuda.texture.CUDAarray(chl, self.nRowsD + 1, self.nColsD + 1, self.nProjSubset[timestep, subset].item())
                                    array.copy_from(yy.reshape((self.nProjSubset[timestep, subset].item(), self.nColsD + 1, self.nRowsD + 1)))
                                    res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
                                    tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp,cp.cuda.runtime.cudaAddressModeClamp),
                                                                            filterMode=cp.cuda.runtime.cudaFilterModeLinear, normalizedCoords=1)
                                    yy = cp.cuda.texture.TextureObject(res, tdes)
                                    y_arg = yy
                                else:
                                    y_arg = yD if self.useTorch else y
                                if self.useTorch:
                                    f_arg = fD
                                else:
                                    f_arg = f[k] if isinstance(f, list) else f
                                args = _KernelArgs(self).seed(kIndLoc)
                                _append_bp5_ct_args(self, args, timestep, subset, k, f_arg, y_arg)
                                kIndLoc = args.as_tuple()
                                if self.meanBP:
                                    kIndLoc += (dMeanBP,)
                        else:
                            y_arg = yD if self.useTorch else y
                            if self.useTorch:
                                f_arg = fD
                            else:
                                f_arg = f[k] if isinstance(f, list) else f
                            # NOTE: preserved as cp.int64 (pre-existing CuPy source),
                            # matching the cl.cltypes.ulong OpenCL/AF uses for the same
                            # logical value -- see _append_bp4_nonct_args.
                            proj_count_arg = cp.int64(self.nProjSubset[timestep, subset].item())
                            args = _KernelArgs(self).seed(kIndLoc)
                            _append_bp4_nonct_args(self, args, timestep, subset, k, f_arg, y_arg, proj_count_arg)
                            kIndLoc = args.as_tuple()
                        kIndLoc += (cp.uint8(self.no_norm),)
                        if self.CT:
                            if self.useMaskBP:
                                kIndLoc += (self.d_maskBP,)
                            kIndLoc += (cp.int64(self.nProjSubset[timestep, subset].item()),)
                        else:
                            kIndLoc += (cp.uint64(self.nMeasSubset[timestep, subset].item()),)
                            kIndLoc += (cp.uint32(subset),)
                        kIndLoc += (cp.int32(k),)
                    self.knlB(
                        (
                            self.globalSizeBP[timestep][subset][k][0] // self.localSizeBP[0],
                            self.globalSizeBP[timestep][subset][k][1] // self.localSizeBP[1],
                            self.globalSizeBP[timestep][subset][k][2],
                        ),
                        (self.localSizeBP[0], self.localSizeBP[1], 1),
                        kIndLoc,
                    )
            if self.useTorch:
                torch.cuda.synchronize()
        else:
            import pyopencl as cl
            import pyopencl.array
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
                    if self.useAF:
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
                    else:
                        # Non-ArrayFire (plain PyOpenCL) path: reproduce the AF computation above
                        # on the host with NumPy -- reshape into (rows, cols, projections), optional
                        # meanBP subtraction, 2-D cumulative sum (summed-area table) over the row and
                        # column axes, then zero-pad one leading row and column before uploading.
                        nProj = int(self.nProjSubset[timestep, subset].item())
                        y_host = np.asarray(y.get()).reshape((self.nRowsD, self.nColsD, nProj), order='F').astype(np.float32)
                        if self.meanBP:
                            meanBP_host = np.mean(y_host, axis=(0, 1)).astype(np.float32)
                            y_host = y_host - meanBP_host.reshape((1, 1, nProj))
                            dMeanBP = cl.Buffer(self.clctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=meanBP_host)
                        y_host = np.cumsum(y_host, axis=0)
                        y_host = np.cumsum(y_host, axis=1)
                        padded = np.zeros((self.nRowsD + 1, self.nColsD + 1, nProj), dtype=np.float32, order='F')
                        padded[1:, 1:, :] = y_host
                        cl.enqueue_copy(self.queue, d_im, padded, origin=(0, 0, 0), region=(self.nRowsD + 1, self.nColsD + 1, nProj));
            
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
                    if self.useAF:
                        y_arg = yD
                        f_arg = fD
                    else:
                        y_arg = y.data
                        f_arg = f[k].data if self.nMultiVolumes > 0 else f.data
                    proj_count_arg = (cl.cltypes.long)(self.nProjSubset[timestep, subset].item())
                    args = _KernelArgs(self)
                    _append_bp123_args(self, args, timestep, subset, k, f_arg, y_arg, proj_count_arg)
                    kIndLoc = args.apply_opencl(self.knlB, kIndLoc)
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
                            if self.useAF:
                                f_arg = fD
                            else:
                                f_arg = f[k].data if isinstance(f, list) else f.data
                            args = _KernelArgs(self)
                            _append_bp4_ct_args(self, args, timestep, subset, k, f_arg, d_im)
                            kIndLoc = args.apply_opencl(self.knlB, kIndLoc)
                        else:
                            if self.useAF:
                                f_arg = fD
                            else:
                                f_arg = f[k].data if isinstance(f, list) else f.data
                            args = _KernelArgs(self)
                            _append_bp5_ct_args(self, args, timestep, subset, k, f_arg, d_im)
                            kIndLoc = args.apply_opencl(self.knlB, kIndLoc)
                            if self.meanBP:
                                self.knlB.set_arg(kIndLoc, dMeanBP)
                                kIndLoc += 1
                    else:
                        if self.useAF:
                            y_arg = yD
                            f_arg = fD
                        else:
                            y_arg = y.data
                            f_arg = f[k].data if isinstance(f, list) else f.data
                        proj_count_arg = (cl.cltypes.ulong)(self.nProjSubset[timestep, subset].item())
                        args = _KernelArgs(self)
                        _append_bp4_nonct_args(self, args, timestep, subset, k, f_arg, y_arg, proj_count_arg)
                        kIndLoc = args.apply_opencl(self.knlB, kIndLoc)
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
                for k in range(self.nMultiVolumes + 1):
                    f[k] = self.computeConvolution(f[k], k)
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


def _proj6_resize_numpy(value: Any, source: tuple[int, int], target: tuple[int, int], physical: tuple[int, int], mode: str, direction: str, expected_frames: int | None, strict: bool, exact: tuple[float | None, float | None] | None = None) -> Any:
    if isinstance(value, (list, tuple)):
        return type(value)(_proj6_resize_numpy(item, source, target, physical, mode, direction, None, strict, exact) for item in value)
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
    original = arr
    if exact is not None:
        # The FOV is not a whole number of detector pixels or the detector and the FOV differ by an odd number of pixels (see
        # exactAxisRatios): resample these axes directly between the detector and the image lattices. The other axes (if any) are
        # padded/trimmed and resized as always below, the resampled axes are already at their target size
        from omegatomo.util.sinogramresample import resampleAxis
        out = arr
        for axis in range(2):
            if exact[axis] is None:
                continue
            ratio = exact[axis] / target[axis] if direction == 'measurement_to_image' else source[axis] / exact[axis]
            out = resampleAxis(out, axis, target[axis], ratio, mode)
        if all(item is not None for item in exact):
            out = out.astype(original.dtype, copy=False)
            return np.asfortranarray(out).ravel(order='F') if flat else out
        physical = tuple(target[axis] if exact[axis] is not None else physical[axis] for axis in range(2))
        arr = out
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
    out = out.astype(original.dtype, copy=False)
    return np.asfortranarray(out).ravel(order='F') if flat else out


def _proj6_resize_torch(value: Any, source: tuple[int, int], target: tuple[int, int], physical: tuple[int, int], mode: str, direction: str, expected_frames: int | None, strict: bool, exact: tuple[float | None, float | None] | None = None) -> Any:
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
    if exact is not None:
        # The FOV is not a whole number of detector pixels or the detector and the FOV differ by an odd number of pixels (see
        # exactAxisRatios): resample these axes directly between the detector and the image lattices. The other axes (if any) are
        # padded/trimmed and resized as always below, the resampled axes are already at their target size
        from omegatomo.util.sinogramresample import axisMatrix
        def matrix(axis: int) -> Any:
            ratio = exact[axis] / target[axis] if direction == 'measurement_to_image' else source[axis] / exact[axis]
            return torch.as_tensor(axisMatrix(source[axis], target[axis], ratio, mode), dtype=out.dtype, device=out.device)
        if exact[1] is not None:
            out = torch.matmul(matrix(1), out)  # columns (the second to last axis)
        if exact[0] is not None:
            out = torch.matmul(out, matrix(0).transpose(0, 1))  # rows (the last axis)
        if all(item is not None for item in exact):
            return out.reshape(-1) if flat else out
        physical = tuple(target[axis] if exact[axis] is not None else physical[axis] for axis in range(2))
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
def resample_resize_proj6(value: Any, options: Any, volume: int = 0, *, direction: str = 'measurement_to_image', mode: str = 'continuous', expected_frames: int | None = None, strict: bool = False, grid: tuple[int, int] | None = None) -> Any:
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

    # Type-6 sums the depth axis Nx and projects onto the transaxial Ny by
    # axial Nz grid. The old Nx usage only happened to work for square volumes.
    image = (dimension('Ny'), dimension('Nz')) if grid is None else tuple(map(int, grid))
    try:
        def scalar(name: str) -> float:
            values = np.asarray(getattr(options, name)).reshape(-1)
            if values.size == 0:
                raise ValueError
            return float(values[min(volume, values.size - 1)])

        pitch_x, pitch_y = scalar('dPitchX'), scalar('dPitchY')
        fov_y, fov_z = scalar('FOVa_y'), scalar('axial_fov')
        dy = scalar('dy')
        dz = scalar('dz')
        if min(pitch_x, pitch_y, fov_y, fov_z, dy, dz) <= 0.0:
            raise ValueError
        if grid is None:
            physical = (max(1, int(round(fov_y / pitch_x))), max(1, int(round(fov_z / pitch_y))))
            exact_fov = (fov_y / pitch_x, fov_z / pitch_y)
        else:
            exact_fov = (image[0] * dy / pitch_x, image[1] * dz / pitch_y)
            physical = (max(1, int(round(exact_fov[0]))), max(1, int(round(exact_fov[1]))))
        # A FOV that is not a whole number of detector pixels (relative tolerance 1e-9), or that differs from the detector by an odd
        # number of pixels, is resampled exactly (the symmetric padding/trimming is then impossible), see sinogramresample.py
        from omegatomo.util.sinogramresample import exactAxisRatios
        exact = exactAxisRatios(exact_fov, measurement, physical)
    except (AttributeError, TypeError, ValueError):
        physical = image
        exact = None
    source, target = (measurement, image) if direction == 'measurement_to_image' else (image, measurement)
    try:
        import torch
    except ImportError:
        torch = None
    if torch is not None and isinstance(value, torch.Tensor):
        return _proj6_resize_torch(value, source, target, physical, mode, direction, expected_frames, strict, exact)
    return _proj6_resize_numpy(value, source, target, physical, mode, direction, expected_frames, strict, exact)


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
    """Type-6 ray normalization, shared with the native rotation projector."""
    del view, nx
    dx = np.asarray(self.dx, dtype=np.float64).reshape(-1)
    dy = np.asarray(self.dy, dtype=np.float64).reshape(-1)
    dz = np.asarray(self.dz, dtype=np.float64).reshape(-1)
    return float((dx[volume] / dx[0]) / (dy[0] * dz[0]))


def _type6_axis_offset(self: Any, volume: int) -> tuple[float, float]:
    """Gantry-axis location in pixel coordinates relative to a volume's centre."""
    n = int(self.Nx[volume]), int(self.Ny[volume])
    spacing = np.asarray(self.dx, dtype=np.float64).reshape(-1), np.asarray(self.dy, dtype=np.float64).reshape(-1)
    bx = np.asarray(self.bx, dtype=np.float64).reshape(-1)
    by = np.asarray(self.by, dtype=np.float64).reshape(-1)
    return (-float(bx[volume]) / spacing[0][volume] - n[0] / 2.,
            -float(by[volume]) / spacing[1][volume] - n[1] / 2.)


def _type6_panel_shift(self: Any, volume: int, view: int, grid_ny: int) -> float:
    """Native type-6 detector-panel shift in the requested volume's voxels."""
    values = np.asarray(getattr(self, 'blurPlanes2F', np.empty(0)), dtype=np.float64).reshape(-1)
    if not values.size:
        return 0.
    dx = np.asarray(self.dx, dtype=np.float64).reshape(-1)
    shift = float(values[view]) * float(dx[0]) / float(dx[volume])

    # Native gridType6 adds half a row only when its detector placement is on
    # the unresampled lattice and the centered detector window starts halfway
    # between rows. Keep the same parity correction for single-resolution grids.
    pad_lo, pad_hi = _type6_padding_z(self, volume)
    ox, oy = _type6_axis_offset(self, 0)
    dz = np.asarray(self.dz, dtype=np.float64).reshape(-1)
    bz = np.asarray(self.bz, dtype=np.float64).reshape(-1)
    centered = (abs(ox) < 1e-3 and abs(oy) < 1e-3 and
                abs(float(bz[0] + int(self.Nz[0]) * dz[0] * .5)) < 1e-3 * float(dz[0]))
    grid_nz = int(self.Nz[volume]) + pad_lo + pad_hi
    resample = (volume > 0 or not centered or pad_lo != pad_hi or
                int(self.nRowsD) > grid_ny or int(self.nColsD) > grid_nz)
    if not resample and (grid_ny - int(self.nRowsD)) % 2:
        shift += .5
    return shift


def _type6_rotation_grid(self: Any, volume: int) -> tuple[int, int]:
    """Native type-6 axis-centred grid dimensions (depth Nx, transverse Ny)."""
    nx, ny = int(self.Nx[volume]), int(self.Ny[volume])
    if volume == 0 and int(getattr(self, 'nMultiVolumes', 0)) == 0 and np.allclose(
        np.asarray(self.bx).reshape(-1)[0] + nx * float(np.asarray(self.dx).reshape(-1)[0]) / 2., 0., atol=1e-3 * float(np.asarray(self.dx).reshape(-1)[0])
    ) and np.allclose(
        np.asarray(self.by).reshape(-1)[0] + ny * float(np.asarray(self.dy).reshape(-1)[0]) / 2., 0., atol=1e-3 * float(np.asarray(self.dy).reshape(-1)[0])
    ):
        return nx, ny
    dx = np.asarray(self.dx, dtype=np.float64).reshape(-1)
    dy = np.asarray(self.dy, dtype=np.float64).reshape(-1)
    bx = np.asarray(self.bx, dtype=np.float64).reshape(-1)
    by = np.asarray(self.by, dtype=np.float64).reshape(-1)
    x0, x1 = bx[volume], bx[volume] + nx * dx[volume]
    y0, y1 = by[volume], by[volume] + ny * dy[volume]
    radius = float(np.sqrt(max(x0 * x0, x1 * x1) + max(y0 * y0, y1 * y1)))
    radius_y = radius
    dx0 = float(dx[0])
    dy0 = float(dy[0])
    if .5 * float(getattr(self, 'nRowsD', 0)) * dy0 > radius:
        shifts = np.asarray(getattr(self, 'blurPlanes2F', np.empty(0)), dtype=np.float64)
        if shifts.size:
            radius_y += float(np.max(np.abs(shifts))) * dx0

    def fit(r: float, spacing: float, parity: int) -> int:
        size = int(np.ceil(2. * r / spacing))
        if size % 2 != parity:
            size += 1
        return max(size, 1)

    return fit(radius, float(dx[volume]), nx % 2 if volume == 0 else 0), fit(
        radius_y, float(dy[volume]), ny % 2 if volume == 0 else 0
    )


def _type6_padding_z(self: Any, volume: int) -> tuple[int, int]:
    """Native CDRF spill padding beyond a multi-resolution volume's axial ends."""
    if int(getattr(self, 'nMultiVolumes', 0)) <= 0:
        return 0, 0
    nx_grid, _ = _type6_rotation_grid(self, volume)
    dx = np.asarray(self.dx, dtype=np.float64).reshape(-1)
    dz = np.asarray(self.dz, dtype=np.float64).reshape(-1)
    bz = np.asarray(self.bz, dtype=np.float64).reshape(-1)
    dist_values = np.asarray(self.blurPlanesF, dtype=np.float64).reshape(-1)
    min_dist = float(np.min(dist_values))
    collimator = (min_dist - .5 * float(self.Nx[0])) * dx[0]
    dist = (collimator + .5 * nx_grid * dx[volume]) / dx[volume]
    z_max = max((nx_grid - .5 - dist) * dx[volume], 0.)
    coeff = np.asarray(self.CDRFCoeffs6, dtype=np.float64).reshape(-1)
    az, bz_cdrf, intrinsic = map(float, (coeff[0], coeff[1], coeff[4]))
    sigma_coeff = 1. / (2. * np.sqrt(2. * np.log(2.)))
    extent = 0.
    for zmm in (0., z_max):
        rg = max(az * zmm + bz_cdrf, 0.)
        if str(getattr(self, 'CDRFModel', 'gaussian')).lower() == 'square':
            value = rg / dz[volume] + 5.2565 * intrinsic / dz[volume] * sigma_coeff
        else:
            value = 5.2565 * np.sqrt(rg * rg + intrinsic * intrinsic) / dz[volume] * sigma_coeff
        extent = max(extent, value)
    radius = int(np.ceil(extent))
    pad_max = ((radius + 7) // 8) * 8
    cols = int(getattr(self, 'measurement_nColsD', self.nColsD))
    z_col_lo = -.5 * cols * dz[0] + .5 * dz[0]
    z_col_hi = .5 * cols * dz[0] - .5 * dz[0]
    pos_lo = (z_col_lo - bz[volume]) / dz[volume] - .5
    pos_hi = (z_col_hi - bz[volume]) / dz[volume] - .5
    need_lo = -int(np.floor(pos_lo + 1e-4))
    need_hi = int(np.ceil(pos_hi - 1e-4)) - (int(self.Nz[volume]) - 1)
    return max(min(pad_max, need_lo), 0), max(min(pad_max, need_hi), 0)


def _type6_rotation_variance(self: Any, volume: int, view: int, grid: tuple[int, int]) -> float:
    """Native rotVarType6 tent-interpolation variance for one view and volume."""
    nx, ny = int(self.Nx[volume]), int(self.Ny[volume])
    out_nx, out_ny = grid
    angle = np.deg2rad(float(self.swivelAngles[view]))
    cosa, sina = float(np.cos(angle)), float(np.sin(angle))
    panel = np.asarray(getattr(self, 'blurPlanes2F', np.zeros(1)), dtype=np.float64).reshape(-1)
    panel_shift = (float(panel[min(view, panel.size - 1)]) if panel.size else 0.)
    dx = np.asarray(self.dx, dtype=np.float64).reshape(-1)
    dy = np.asarray(self.dy, dtype=np.float64).reshape(-1)
    panel_shift *= dx[0] / dx[volume]
    axis_x, axis_y = _type6_axis_offset(self, volume)
    off_x, off_y = -sina * panel_shift + axis_x, cosa * panel_shift + axis_y
    hh = abs(cosa) + abs(sina)
    csx, csy = .5 * nx - .5, .5 * ny - .5
    cox, coy = .5 * out_nx - .5, .5 * out_ny - .5
    aligned = abs(sina) < 1e-6 or abs(cosa) < 1e-6
    n_voxels = nx * ny
    if aligned:
        samples = [(nx // 2, ny // 2)]
    elif n_voxels <= 1024:
        samples = [(i % nx, i // nx) for i in range(n_voxels)]
    else:
        def halton(index: int, base: int) -> float:
            factor, result = 1., 0.
            while index:
                factor /= base
                result += factor * (index % base)
                index //= base
            return result
        samples = [(min(int(halton(i + 1, 2) * nx), nx - 1),
                    min(int(halton(i + 1, 3) * ny), ny - 1)) for i in range(1024)]
    numerator = denominator = 0.
    for ix, iy in samples:
        ux, uy = ix - csx, iy - csy
        qx = cosa * (ux - off_x) + sina * (uy - off_y)
        qy = -sina * (ux - off_x) + cosa * (uy - off_y)
        jx0, jx1 = int(np.ceil(qx + cox - hh)), int(np.floor(qx + cox + hh))
        jy0, jy1 = int(np.ceil(qy + coy - hh)), int(np.floor(qy + coy + hh))
        for jy in range(jy0, jy1 + 1):
            ay = jy - coy
            for jx in range(jx0, jx1 + 1):
                ax = jx - cox
                sx = cosa * ax - sina * ay + off_x
                sy = sina * ax + cosa * ay + off_y
                weight = max(1. - abs(sx - ux), 0.) * max(1. - abs(sy - uy), 0.)
                numerator += weight * (ay - qy) ** 2
                denominator += weight
    variance = numerator / denominator if denominator else 1. / 6.
    if volume > 0:
        n_det = int(getattr(self, 'nRowsD', 0))
        a = float(dy[0] / dy[volume])
        b = .5 * out_ny + (.5 - .5 * n_det) * a - .5
        num = den = 0.
        for row in range(n_det):
            t = a * row + b
            low = np.floor(t)
            for offset in (0., 1.):
                d = t - (low + offset)
                weight = max(1. - abs(d), 0.)
                num += weight * d * d
                den += weight
        if den:
            variance += num / den
    return float(variance)


def _type6_on_the_fly_cdrf(self: Any, volume: int, view: int) -> np.ndarray:
    """Native type-6 per-view CDRF sampled on the rotated grid."""
    model = str(getattr(self, 'CDRFModel', 'gaussian')).strip().lower()
    if model not in ('gaussian', 'square'):
        raise ValueError(f"Unknown type-6 CDRF model: {model!r}.")
    grid_nx, grid_ny = _type6_rotation_grid(self, volume)
    grid_nz = int(self.Nz[volume]) + sum(_type6_padding_z(self, volume))
    dx = np.asarray(self.dx, dtype=np.float64).reshape(-1)
    dy = np.asarray(self.dy, dtype=np.float64).reshape(-1)
    dz = np.asarray(self.dz, dtype=np.float64).reshape(-1)
    dist_float = np.asarray(self.blurPlanesF, dtype=np.float64).reshape(-1)
    collimator = (dist_float[view] - .5 * float(self.Nx[0])) * dx[0]
    distance = (collimator + .5 * grid_nx * dx[volume]) / dx[volume]
    first_plane = int(np.clip(np.ceil(distance - .5), 0, grid_nx - 1))
    first_plane = (first_plane // 16) * 16
    coeff = np.asarray(self.CDRFCoeffs6, dtype=np.float64).reshape(-1)
    az, bz, axy, bxy, intrinsic = map(float, coeff[:5])
    sigma_coeff = 1. / (2. * np.sqrt(2. * np.log(2.)))
    import os
    if os.environ.get('TYPE6_ROTCOMP', '1') != '0':
        comp = _type6_rotation_variance(self, volume, view, (grid_nx, grid_ny))
    else:
        comp = sigma_coeff * sigma_coeff
    planes = np.arange(grid_nx, dtype=np.float64)
    depth_mm = (planes + .5 - distance) * dx[volume]
    rg_z = np.maximum(az * depth_mm + bz, 0.)
    rg_xy = np.maximum(axy * depth_mm + bxy, 0.)
    fwhm_z = np.sqrt(rg_z * rg_z + intrinsic * intrinsic) / dz[volume]
    fwhm_xy = np.sqrt(rg_xy * rg_xy + intrinsic * intrinsic) / dy[volume]
    sigma_z = fwhm_z * sigma_coeff
    sigma_xy = np.sqrt(np.maximum((fwhm_xy * sigma_coeff) ** 2 - comp,
                                  1e-16 * sigma_coeff * sigma_coeff))
    # psfHelperType6 chooses a common radius from the two depth endpoints,
    # rounds it up for stable work sizes, and clamps it to the rotated grid.
    endpoint_depths = ((first_plane + .5 - distance) * dx[volume],
                       (grid_nx - .5 - distance) * dx[volume])
    extent = 0.
    for zmm in endpoint_depths:
        rz = max(az * zmm + bz, 0.)
        rxy = max(axy * zmm + bxy, 0.)
        if model == 'square':
            extent_z = (rz + 5.2565 * intrinsic * sigma_coeff) / dz[volume]
            extent_xy = (rxy + 5.2565 * intrinsic * sigma_coeff) / dy[volume]
        else:
            extent_z = 5.2565 * np.hypot(rz, intrinsic) / dz[volume] * sigma_coeff
            extent_xy = 5.2565 * np.hypot(rxy, intrinsic) / dy[volume] * sigma_coeff
        extent = max(extent, extent_z, extent_xy)
    radius = int(np.ceil(extent))
    radius = ((radius + 7) // 8) * 8
    radius = min(radius, max(grid_ny, grid_nz))
    taps = np.arange(-radius, radius + 1, dtype=np.float64)
    visible = (planes >= first_plane) & (depth_mm > 0.)
    if model == 'square':
        from scipy.special import erfc

        def square_response(rg, spacing, compensation):
            variance = (intrinsic / spacing * sigma_coeff) ** 2 - compensation
            sigma = np.sqrt(max(variance, 1e-6))
            half_base = np.sqrt(np.maximum((rg / spacing) ** 2 - 6. * max(1e-6 - variance, 0.), 0.))
            u = np.abs(taps[:, None])
            r = half_base[None, :]
            a = r / sigma
            z2 = np.minimum((u / sigma) ** 2, 100.)
            a2 = a * a
            he2 = z2 - 1.
            he4 = (z2 - 6.) * z2 + 3.
            he6 = ((z2 - 15.) * z2 + 45.) * z2 - 15.
            he8 = (((z2 - 28.) * z2 + 210.) * z2 - 420.) * z2 + 105.
            series = np.exp(-.5 * z2) / (np.sqrt(2. * np.pi) * sigma) * (
                1. + a2 * (he2 / 12. + a2 * (he4 / 360. + a2 * (he6 / 20160. + a2 * he8 / 1814400.))))
            series = np.where((u / sigma) ** 2 <= 100., series, 0.)

            def tail(x):
                z = x / sigma
                return sigma * np.exp(-.5 * z * z) / np.sqrt(2. * np.pi) - .5 * x * erfc(z / np.sqrt(2.))

            closed = (tail(u + r) - 2. * tail(u) + tail(np.abs(u - r)) + np.maximum(r - u, 0.))
            closed /= np.maximum(r * r, 1e-30)
            response = np.maximum(np.where(a < .5, series, closed), 0.)
            cutoff = np.minimum(np.ceil(half_base + 5.2565 * sigma).astype(np.int64), radius)
            return response, cutoff

        axial, rz = square_response(rg_z, dz[volume], 0.)
        trans, rxy = square_response(rg_xy, dy[volume], comp)
    else:
        rz = np.minimum(np.ceil(5.2565 * sigma_z).astype(np.int64), radius)
        rxy = np.minimum(np.ceil(5.2565 * sigma_xy).astype(np.int64), radius)
        positive_sigma_z = sigma_z > 0.
        axial = np.empty((taps.size, sigma_z.size), dtype=np.float64)
        axial[:, positive_sigma_z] = np.exp(
            -(taps[:, None] ** 2) / (2. * sigma_z[None, positive_sigma_z] ** 2)
        )
        # A zero Gaussian width is the delta response: all weight belongs at
        # the center tap. Evaluating 0 / 0 first would leave NaNs after masking.
        axial[:, ~positive_sigma_z] = (taps[:, None] == 0.)
        trans = np.exp(-(taps[:, None] ** 2) / (2. * sigma_xy[None, :] ** 2))
    axial *= (np.abs(taps[:, None]) <= rz[None, :]) & visible[None, :]
    trans *= (np.abs(taps[:, None]) <= rxy[None, :]) & visible[None, :]
    axial /= np.maximum(axial.sum(axis=0, keepdims=True), 1e-30)
    trans /= np.maximum(trans.sum(axis=0, keepdims=True), 1e-30)
    # Torch's (z,y,x) array is convolved as (x,y,z): row taps are
    # transaxial and column taps are axial.
    return np.asfortranarray(trans[:, None, :] * axial[None, :, :], dtype=np.float32)


def _type6_attenuation_map_extent(self: Any) -> tuple[float, float, float, int, int, int]:
    """Union grid of the type-6 attenuation map, at the main voxel spacing."""
    nvol = int(self.nMultiVolumes) + 1
    origins = np.column_stack([np.asarray(getattr(self, name), dtype=np.float32).reshape(-1)[:nvol]
                               for name in ('bx', 'by', 'bz')])
    counts = np.column_stack([np.asarray(getattr(self, name), dtype=np.float32).reshape(-1)[:nvol]
                              for name in ('Nx', 'Ny', 'Nz')])
    spacing = np.column_stack([np.asarray(getattr(self, name), dtype=np.float32).reshape(-1)[:nvol]
                               for name in ('dx', 'dy', 'dz')])
    low = np.min(origins, axis=0)
    high = np.max(origins + counts * spacing, axis=0)
    dims = np.ceil(((high - low) / spacing[0]).astype(np.float64) - 1e-5).astype(np.int64)
    return float(low[0]), float(low[1]), float(low[2]), int(dims[0]), int(dims[1]), int(dims[2])


def _type6_resample_attenuation_numpy(attenMap: np.ndarray, self: Any, volume: int,
                                       bxMap: float, byMap: float, bzMap: float,
                                       NxMap: int, NyMap: int, NzMap: int) -> np.ndarray:
    """Box-average-resamples `attenMap` (shape (NzMap, NyMap, NxMap), the dense map covering
    the union bounding box from _type6_attenuation_map_extent) onto multi-resolution volume
    `volume`'s own (Nz[volume], Ny[volume], Nx[volume]) grid, using PHYSICAL position
    (bx/by/bz) and voxel size (dx/dy/dz) -- not assuming any volume is centred on the origin.

    Implemented as a mean-preserving box average via a 3D summed-area table (SAT): each
    target voxel's exact physical box is mapped into fractional source-map index space, the
    (zero-padded) SAT is trilinearly sampled at the box's 8 corners, and the box's mean is the
    standard inclusion-exclusion box-sum divided by its volume (in source-voxel units). This
    is preferred over plain trilinear point-sampling because multi-resolution side volumes are
    normally coarser than the map (multiResolutionScale < 1) -- a point sample would alias the
    many source voxels each target voxel actually covers -- and the same formula reduces to
    (numerically) the same result as point-sampling when a target voxel is no coarser than the
    source, so one formula covers both directions. Edge-clamping of the box corners only
    guards against floating-point spill a fraction of a voxel past the map's own edge -- the
    volumes are built to exactly tile the map's extent. Mirrors
    source/cpp/functions.hpp's resampleAttenuationType6Volume exactly, in pure NumPy so the
    identical algorithm can be shared by every Python backend (Torch, ArrayFire)."""
    nx, ny, nz = int(self.Nx[volume]), int(self.Ny[volume]), int(self.Nz[volume])
    if (volume == 0 and nx == NxMap and ny == NyMap and nz == NzMap
            and np.isclose(float(self.bx[0]), bxMap) and np.isclose(float(self.by[0]), byMap)
            and np.isclose(float(self.bz[0]), bzMap)):
        return attenMap  # Identity: no multi-resolution split, the map already IS volume 0's own grid

    sat = np.cumsum(np.cumsum(np.cumsum(attenMap.astype(np.float64), axis=0), axis=1), axis=2)
    sat = np.pad(sat, ((1, 0), (1, 0), (1, 0)))

    dx0, dy0, dz0 = float(self.dx[0]), float(self.dy[0]), float(self.dz[0])
    bxV, byV, bzV = float(self.bx[volume]), float(self.by[volume]), float(self.bz[volume])
    dxV, dyV, dzV = float(self.dx[volume]), float(self.dy[volume]), float(self.dz[volume])
    zEdges = np.clip((bzV + np.arange(nz + 1, dtype=np.float64) * dzV - bzMap) / dz0, 0., float(NzMap))
    yEdges = np.clip((byV + np.arange(ny + 1, dtype=np.float64) * dyV - byMap) / dy0, 0., float(NyMap))
    xEdges = np.clip((bxV + np.arange(nx + 1, dtype=np.float64) * dxV - bxMap) / dx0, 0., float(NxMap))

    z0 = np.clip(np.floor(zEdges).astype(np.int64), 0, NzMap); z1 = np.minimum(z0 + 1, NzMap)
    y0 = np.clip(np.floor(yEdges).astype(np.int64), 0, NyMap); y1 = np.minimum(y0 + 1, NyMap)
    x0 = np.clip(np.floor(xEdges).astype(np.int64), 0, NxMap); x1 = np.minimum(x0 + 1, NxMap)
    fz = (zEdges - z0).reshape(-1, 1, 1)
    fy = (yEdges - y0).reshape(1, -1, 1)
    fx = (xEdges - x0).reshape(1, 1, -1)

    def gather(zi: np.ndarray, yi: np.ndarray, xi: np.ndarray) -> np.ndarray:
        return sat[np.ix_(zi, yi, xi)]

    c000, c001 = gather(z0, y0, x0), gather(z0, y0, x1)
    c010, c011 = gather(z0, y1, x0), gather(z0, y1, x1)
    c100, c101 = gather(z1, y0, x0), gather(z1, y0, x1)
    c110, c111 = gather(z1, y1, x0), gather(z1, y1, x1)
    c00 = c000 * (1. - fx) + c001 * fx
    c01 = c010 * (1. - fx) + c011 * fx
    c10 = c100 * (1. - fx) + c101 * fx
    c11 = c110 * (1. - fx) + c111 * fx
    c0 = c00 * (1. - fy) + c01 * fy
    c1 = c10 * (1. - fy) + c11 * fy
    corners = c0 * (1. - fz) + c1 * fz  # shape (nz+1, ny+1, nx+1)

    boxSum = (
          corners[1:, 1:, 1:] - corners[:-1, 1:, 1:] - corners[1:, :-1, 1:] - corners[1:, 1:, :-1]
        + corners[:-1, :-1, 1:] + corners[:-1, 1:, :-1] + corners[1:, :-1, :-1] - corners[:-1, :-1, :-1]
    )
    dZe = np.diff(zEdges).reshape(-1, 1, 1)
    dYe = np.diff(yEdges).reshape(1, -1, 1)
    dXe = np.diff(xEdges).reshape(1, 1, -1)
    # Native `resampleAttenuationType6Host(..., clampedVolu=false)` divides by the
    # full target voxel volume even when a target box extends beyond the union map.
    # The out-of-map part therefore contributes zero attenuation. Dividing by the
    # clipped overlap volume incorrectly extends boundary attenuation into those
    # uncovered regions, especially for the larger side-volume source lattices.
    volu_nominal = max((float(dxV) / dx0) * (float(dyV) / dy0) * (float(dzV) / dz0), 1e-6)
    volu = volu_nominal
    return (boxSum / volu).astype(np.float32)


def _type6_attenuation_source_numpy(self: Any, attenMap: np.ndarray, volume: int, *,
                                   bz_override: float | None = None,
                                   nz_override: int | None = None) -> tuple[np.ndarray, int]:
    """Resample the union attenuation map onto native type-6's axis-centred source lattice.

    The returned source is `(Nz, srcNy, srcNx)` and `attNx` is the x extent after rotation.
    Keeping the source lattice separate from the emission volume is essential: attenuation
    from every volume along the detector ray contributes before the emission volume is cropped.
    """
    from types import SimpleNamespace

    bxMap, byMap, bzMap, NxMap, NyMap, NzMap = _type6_attenuation_map_extent(self)
    expected = (NzMap, NyMap, NxMap)
    attenMap = np.asarray(attenMap, dtype=np.float32)
    if attenMap.shape != expected:
        raise ValueError(f'Expected type-6 union attenuation map shape {expected}, got {attenMap.shape}')

    gridNx, gridNy = _type6_rotation_grid(self, volume)
    axisX = axisY = 0.0  # Native axisType6 currently returns the origin.
    runion = 0.0
    for v in range(int(self.nMultiVolumes) + 1):
        x0 = float(self.bx[v]) - axisX
        x1 = x0 + int(self.Nx[v]) * float(self.dx[v])
        y0 = float(self.by[v]) - axisY
        y1 = y0 + int(self.Ny[v]) * float(self.dy[v])
        runion = max(runion, float(np.sqrt(max(x0*x0, x1*x1) + max(y0*y0, y1*y1))))
    minNx = int(np.ceil(2.0 * runion / float(self.dx[volume])))
    minNy = int(np.ceil(2.0 * runion / float(self.dy[volume])))
    attNx = minNx + ((minNx - gridNx) % 2)
    srcNx = attNx
    srcNy = minNy + ((minNy - gridNy) % 2)

    # Reuse the shared SAT box-average resampler with an axis-centred target lattice.
    source_bz = float(self.bz[volume]) if bz_override is None else float(bz_override)
    source_nz = int(self.Nz[volume]) if nz_override is None else int(nz_override)
    # The resampler reads index 0 as the map lattice and `volume` as the target.
    # Keep those separate: using a single target entry made coarse side volumes
    # silently reinterpret the union map at their coarse voxel spacing.
    proxy = SimpleNamespace(
        Nx=[NxMap, srcNx], Ny=[NyMap, srcNy], Nz=[NzMap, source_nz],
        dx=[float(self.dx[0]), float(self.dx[volume])],
        dy=[float(self.dy[0]), float(self.dy[volume])],
        dz=[float(self.dz[0]), float(self.dz[volume])],
        bx=[bxMap, axisX - .5 * srcNx * float(self.dx[volume])],
        by=[byMap, axisY - .5 * srcNy * float(self.dy[volume])],
        bz=[bzMap, source_bz],
    )
    source = _type6_resample_attenuation_numpy(
        attenMap, proxy, 1, bxMap, byMap, bzMap, NxMap, NyMap, NzMap)
    return source, attNx


def _type6_attenuation_group(self: Any, volume: int) -> tuple[list[int], int, float, int, int, list[int]]:
    """Return the native attenuation group and its shared rotation/source grid."""
    nvol = int(self.nMultiVolumes) + 1
    grids = [_type6_rotation_grid(self, v) for v in range(nvol)]
    bx = np.asarray(self.bx, dtype=np.float64).reshape(-1)
    by = np.asarray(self.by, dtype=np.float64).reshape(-1)
    bz = np.asarray(self.bz, dtype=np.float64).reshape(-1)
    dx = np.asarray(self.dx, dtype=np.float64).reshape(-1)
    dy = np.asarray(self.dy, dtype=np.float64).reshape(-1)
    dz = np.asarray(self.dz, dtype=np.float64).reshape(-1)
    nx = np.asarray(self.Nx, dtype=np.int64).reshape(-1)
    ny = np.asarray(self.Ny, dtype=np.int64).reshape(-1)
    nz = np.asarray(self.Nz, dtype=np.int64).reshape(-1)
    runion = 0.0
    for v in range(nvol):
        x0, x1 = float(bx[v]), float(bx[v] + nx[v] * dx[v])
        y0, y1 = float(by[v]), float(by[v] + ny[v] * dy[v])
        runion = max(runion, np.sqrt(max(x0*x0, x1*x1) + max(y0*y0, y1*y1)))

    specs = []
    for v, (grid_x, grid_y) in enumerate(grids):
        min_nx = int(np.ceil(2.0 * runion / dx[v]))
        min_ny = int(np.ceil(2.0 * runion / dy[v]))
        att_nx = min_nx + ((min_nx - grid_x) % 2)
        src_ny = min_ny + ((min_ny - grid_y) % 2)
        specs.append((int(att_nx), int(att_nx), int(src_ny)))

    groups: list[list[int]] = []
    for v in range(nvol):
        found = None
        for members in groups:
            first = members[0]
            if (dx[v] != dx[first] or dy[v] != dy[first] or dz[v] != dz[first]
                    or specs[v] != specs[first]
                    or (grids[v][0] - grids[first][0]) % 2
                    or (grids[v][1] - grids[first][1]) % 2):
                continue
            offset = (bz[v] - bz[first]) / dz[v]
            if abs(offset - np.round(offset)) < 1e-3:
                found = members
                break
        if found is None:
            found = []
            groups.append(found)
        found.append(v)

    members = next(group for group in groups if volume in group)
    group_bz = float(min(bz[v] for v in members))
    z_offsets = [0] * nvol
    group_nz = 0
    rows = 0
    for v in members:
        z_offsets[v] = int(np.rint((bz[v] - group_bz) / dz[v]))
        group_nz = max(group_nz, z_offsets[v] + int(nz[v]))
        rows = max(rows, grids[v][1])
    return members, rows, group_bz, group_nz, specs[members[0]][0], z_offsets


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


def _type6_numpy_central_ray_weights(self: Any, projections: int, timestep: int, subset: int) -> np.ndarray:
    """Central-ray voxel-volume/chord diagonal in (view, detector-column, row) order."""
    rows = int(getattr(self, 'measurement_nRowsD', self.nRowsD))
    cols = int(getattr(self, 'measurement_nColsD', self.nColsD))
    row_grid = (np.arange(rows, dtype=np.float64) - .5 * (rows - 1))[:, None] * float(self.dPitchX)
    col_grid = (np.arange(cols, dtype=np.float64) - .5 * (cols - 1))[None, :] * float(self.dPitchY)
    dx0, dy0, dz0 = (float(np.asarray(getattr(self, key)).reshape(-1)[0]) for key in ('dx', 'dy', 'dz'))
    voxel_volume = dx0 * dy0 * dz0
    centers = (float(self.ellipseCenterX), float(self.ellipseCenterY), float(self.ellipseCenterZ))
    radii = (float(self.ellipseRadiusX), float(self.ellipseRadiusY), float(self.ellipseRadiusZ))
    swivel = np.asarray(self.swivelAngles, dtype=np.float64).reshape(-1)
    angles = np.asarray(self.angles, dtype=np.float64).reshape(-1)
    radii_view = np.asarray(self.radiusPerProj, dtype=np.float64).reshape(-1)
    # `distFloat2` is the native C++ Weighting name. Python options expose the
    # same per-view detector-panel shift as `blurPlanes2F` (in main-volume voxels).
    panel_values = np.asarray(getattr(self, 'blurPlanes2F', np.zeros_like(swivel)), dtype=np.float64)
    if panel_values.size == 0:
        panel_values = np.zeros_like(swivel)
    start = int(np.sum(self.nProjSubset[:timestep, :]) + np.sum(self.nProjSubset[timestep, :subset]))
    result = np.zeros((projections, cols, rows), dtype=np.float32)
    for local in range(projections):
        view = start + local
        theta = np.deg2rad(swivel[view])
        nx, ny = np.cos(theta), np.sin(theta)
        if panel_values.ndim == 2:
            offset_voxels = panel_values[0, view]
        else:
            offset_voxels = panel_values.reshape(-1)[view]
        offset = float(offset_voxels) * dx0
        qx, qy = -ny * (row_grid + offset), nx * (row_grid + offset)
        valid = np.broadcast_to(np.abs(col_grid - centers[2]) <= radii[2], (rows, cols)).copy()
        source_n = radii_view[view] * np.cos(np.deg2rad(angles[view] - swivel[view])) + float(self.CORtoDetectorSurface)
        if np.isinf(float(self.ellipsePower)):
            lo, hi = np.full((rows, cols), -np.inf), np.full((rows, cols), np.inf)
            for direction, coord, center, extent in ((nx, qx, centers[0], radii[0]), (ny, qy, centers[1], radii[1])):
                if abs(direction) < 1e-8:
                    valid &= (coord >= center - extent) & (coord <= center + extent)
                else:
                    t0, t1 = (center - extent - coord) / direction, (center + extent - coord) / direction
                    lo, hi = np.maximum(lo, np.minimum(t0, t1)), np.minimum(hi, np.maximum(t0, t1))
            chord = np.maximum(hi - np.maximum(lo, source_n), 0.)
        elif float(self.ellipsePower) == 2.:
            inv_x2, inv_y2 = 1. / radii[0] ** 2, 1. / radii[1] ** 2
            aa = nx * nx * inv_x2 + ny * ny * inv_y2
            px, py = qx - centers[0], qy - centers[1]
            bb = 2. * (px * nx * inv_x2 + py * ny * inv_y2)
            cc = px * px * inv_x2 + py * py * inv_y2 - 1.
            disc = bb * bb - 4. * aa * cc
            root = np.sqrt(np.maximum(disc, 0.))
            exit_t, entry_t = (-bb + root) / (2. * aa), (-bb - root) / (2. * aa)
            chord = np.maximum(exit_t - np.maximum(entry_t, source_n), 0.)
            valid &= disc > 1e-12
        else:
            raise ValueError('Type-6 path normalization requires box or elliptical SPECT support')
        result[local] = np.where(valid & (chord > 1e-6), voxel_volume / np.maximum(chord, 1e-6), 0.).T
    return result


def _type6_torch_measurement_weights(self: Any, data: Any, projections: int, timestep: int, subset: int) -> Any:
    """Apply type-6 measurement-domain corrections on any Torch device."""
    import torch

    rows = int(getattr(self, 'measurement_nRowsD', self.nRowsD))
    cols = int(getattr(self, 'measurement_nColsD', self.nColsD))
    pixels = rows * cols
    device = data.device
    detector = _type6_torch_detector_vector(self, timestep, subset, device)
    weights = torch.ones_like(data)

    weights *= torch.as_tensor(
        _type6_numpy_central_ray_weights(self, projections, timestep, subset),
        dtype=data.dtype, device=device,
    )

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
    if int(getattr(self, 'scatter_correction', 0)) == 1:
        value = _type6_torch_resource(
            self, 'd_scatter', timestep, subset,
            fallback=getattr(self, 'corrVector', np.empty(0)), dtype=torch.float32, device=device,
        )
        value = expand(value)
        if value is not None:
            weights *= value
    global_factor = float(getattr(self, 'global_factor', 1.))
    if global_factor != 1.:
        weights *= global_factor
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

    _weights_unset = object()  # sentinel distinguishing "attribute missing" from any real id()
    # Ops are constructed once per FP/BP call; cache the host-side SAT resample
    # for every view in that call without risking stale attenuation across calls.
    attenuation_source_cache: dict[tuple, tuple[np.ndarray, int]] = {}

    def _weights_cache_validity(timestep: int, subset: int) -> tuple:
        # measurement_weights only ever reads self.maskFP / self.vaimennus /
        # self.normalization as its host-side sources (host_resource() above
        # ignores its 'name' argument and always uses the caller-supplied
        # fallback, which is one of exactly these three attributes) -- so
        # the composed AF weights array for a given (timestep, subset) stays
        # valid exactly as long as none of these three has been reassigned
        # to a different object (mutating one of them in place would be
        # invisible to id() and is not a case that arises in this codebase).
        return (
            id(getattr(self, 'maskFP', _weights_unset)),
            id(getattr(self, 'vaimennus', _weights_unset)),
            id(getattr(self, 'normalization', _weights_unset)),
            id(getattr(self, 'corrVector', _weights_unset)),
            float(getattr(self, 'global_factor', 1.)),
        )

    def measurement_weights(data: Any, projections: int, timestep: int, subset: int) -> Any:
        cache = getattr(self, '_type6_af_weights_cache', None)
        if cache is None:
            cache = {}
            self._type6_af_weights_cache = cache
        cache_key = (timestep, subset)
        validity = _weights_cache_validity(timestep, subset)
        cached = cache.get(cache_key)
        if cached is not None and cached[0] == validity:
            return cached[1]

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
        if int(getattr(self, 'scatter_correction', 0)) == 1:
            value = expand('d_scatter', getattr(self, 'corrVector', np.empty(0)))
            if value is not None:
                weights *= value
        weights *= float(getattr(self, 'global_factor', 1.))
        # The native type-6 projector applies the central-ray diagonal in
        # addition to the user supplied normalization array.  Keep the AF
        # path on the same detector/view grid as the native and Torch paths.
        weights *= _type6_numpy_central_ray_weights(self, projections, timestep, subset)
        result = from_vector(weights.ravel(order='F'), projections, cols, rows)
        cache[cache_key] = (validity, result)
        return result

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
        if int(value.shape[0]) == 1 and len(value.shape) > 2:
            value = af.data.moddims(value, int(value.shape[1]), int(value.shape[2]))
        return af.image.resize(value, scale=None, odim0=int(shape[0]), odim1=int(shape[1]), method=af.INTERP.BILINEAR)

    def resize_stack(value: Any, shape: tuple[int, int]) -> Any:
        # AF resize acts on dimensions zero and one.  Present the detector
        # stack as (x, z, view), then restore the shared (view, z, x) form.
        resized = af.image.resize(
            af.data.reorder(value, 2, 1, 0), scale=None, odim0=int(shape[1]), odim1=int(shape[0]),
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

    def physical_grid(volume: int, grid: tuple[int, int] | None = None) -> tuple[int, int]:
        if grid is not None:
            try:
                pitch_x = float(np.asarray(self.dPitchX).reshape(-1)[0])
                pitch_y = float(np.asarray(self.dPitchY).reshape(-1)[0])
                dy = float(np.asarray(self.dy).reshape(-1)[min(volume, np.asarray(self.dy).size - 1)])
                dz = float(np.asarray(self.dz).reshape(-1)[min(volume, np.asarray(self.dz).size - 1)])
                return max(1, int(round(int(grid[0]) * dy / pitch_x))), max(1, int(round(int(grid[1]) * dz / pitch_y)))
            except (AttributeError, TypeError, ValueError, IndexError, ZeroDivisionError):
                return int(grid[0]), int(grid[1])
        try:
            def scalar(name: str) -> float:
                values = np.asarray(getattr(self, name)).reshape(-1)
                if values.size == 0:
                    raise ValueError
                return float(values[min(volume, values.size - 1)])

            pitch_x, pitch_y = scalar('dPitchX'), scalar('dPitchY')
            fov_y, fov_z = scalar('FOVa_y'), scalar('axial_fov')
            if min(pitch_x, pitch_y, fov_y, fov_z) <= 0.0:
                raise ValueError
            return max(1, int(round(fov_y / pitch_x))), max(1, int(round(fov_z / pitch_y)))
        except (AttributeError, TypeError, ValueError):
            return int(self.Ny[volume]), int(self.Nz[volume])

    def exact_grid(volume: int, grid: tuple[int, int] | None = None) -> tuple[float | None, float | None] | None:
        # The FOV in detector pixels (x, z) for the axes that are resampled directly: not a whole number (relative tolerance 1e-9) or
        # an odd difference to the detector size (None for the other axes, None if there are no such axes), see sinogramresample.py
        try:
            def scalar(name: str) -> float:
                values = np.asarray(getattr(self, name)).reshape(-1)
                if values.size == 0:
                    raise ValueError
                return float(values[min(volume, values.size - 1)])

            pitch_x, pitch_y = scalar('dPitchX'), scalar('dPitchY')
            fov_y, fov_z = scalar('FOVa_y'), scalar('axial_fov')
            if min(pitch_x, pitch_y, fov_y, fov_z) <= 0.0:
                raise ValueError
            from omegatomo.util.sinogramresample import exactAxisRatios
            measurement = (int(getattr(self, 'measurement_nRowsD', self.nRowsD)), int(getattr(self, 'measurement_nColsD', self.nColsD)))
            physical = physical_grid(volume, grid)
            if grid is None:
                exact_fov = (fov_y / pitch_x, fov_z / pitch_y)
            else:
                dy, dz = scalar('dy'), scalar('dz')
                exact_fov = (int(grid[0]) * dy / pitch_x, int(grid[1]) * dz / pitch_y)
            return exactAxisRatios(exact_fov, measurement, physical)
        except (AttributeError, TypeError, ValueError):
            return None

    def apply_axis_matrix(value: Any, matrix: np.ndarray, dim: int) -> Any:
        # Multiplies the dimension dim (1 = columns, 2 = rows) of the (view, column, row) stack by the (out, in) matrix
        weights = af.interop.np_to_af_array(np.asfortranarray(matrix.astype(np.float32)))
        views, current_cols, current_rows = (int(value.shape[0]), int(value.shape[1]), int(value.shape[2]))
        if dim == 2:
            moved = af.data.moddims(af.data.reorder(value, 2, 0, 1), current_rows, views * current_cols)
            moved = af.data.moddims(af.blas.matmul(weights, moved), int(matrix.shape[0]), views, current_cols)
            return af.data.reorder(moved, 1, 2, 0)
        moved = af.data.moddims(af.data.reorder(value, 1, 0, 2), current_cols, views * current_rows)
        moved = af.data.moddims(af.blas.matmul(weights, moved), int(matrix.shape[0]), views, current_rows)
        return af.data.reorder(moved, 1, 0, 2)

    def detector_axis_matrices(volume: int, grid: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
        """Native tent-sampling maps from rotated image grid (Ny,Nz) to detector."""
        rows = int(getattr(self, 'measurement_nRowsD', self.nRowsD))
        cols = int(getattr(self, 'measurement_nColsD', self.nColsD))
        dy = np.asarray(self.dy, dtype=np.float64).reshape(-1)
        dz = np.asarray(self.dz, dtype=np.float64).reshape(-1)
        bz = np.asarray(self.bz, dtype=np.float64).reshape(-1)
        pad_lo, _ = _type6_padding_z(self, volume)
        trans_a = float(dy[0] / dy[volume])
        trans_b = .5 * int(grid[0]) + (.5 - .5 * rows) * trans_a - .5
        axial_a = float(dz[0] / dz[volume])
        axial_b = (-float(bz[volume]) / dz[volume]
                   + (.5 - .5 * cols) * axial_a - .5 + pad_lo)
        det_rows = np.arange(rows, dtype=np.float64)[:, None]
        src_rows = np.arange(int(grid[0]), dtype=np.float64)[None, :]
        det_cols = np.arange(cols, dtype=np.float64)[:, None]
        src_cols = np.arange(int(grid[1]), dtype=np.float64)[None, :]
        trans = np.maximum(1. - np.abs(det_rows * trans_a + trans_b - src_rows), 0.).astype(np.float32)
        axial = np.maximum(1. - np.abs(det_cols * axial_a + axial_b - src_cols), 0.).astype(np.float32)
        return trans, axial

    def exact_resample(value: Any, volume: int, direction: str, exact: tuple[float, float], grid: tuple[int, int] | None = None) -> Any:
        from omegatomo.util.sinogramresample import axisMatrix
        rows = int(getattr(self, 'measurement_nRowsD', self.nRowsD))
        cols = int(getattr(self, 'measurement_nColsD', self.nColsD))
        nx, nz = physical_grid(volume, grid)
        current_cols, current_rows = int(value.shape[1]), int(value.shape[2])
        if direction == 'measurement_to_image':
            if exact[0] is not None:
                value = apply_axis_matrix(value, axisMatrix(current_rows, nx, exact[0] / nx), 2)
            if exact[1] is not None:
                value = apply_axis_matrix(value, axisMatrix(current_cols, nz, exact[1] / nz), 1)
        else:
            if exact[0] is not None:
                value = apply_axis_matrix(value, axisMatrix(current_rows, rows, current_rows / exact[0]), 2)
            if exact[1] is not None:
                value = apply_axis_matrix(value, axisMatrix(current_cols, cols, current_cols / exact[1]), 1)
        return value

    def resample(value: Any, volume: int, direction: str, grid: tuple[int, int] | None = None) -> Any:
        if grid is not None:
            trans, axial = detector_axis_matrices(volume, grid)
            if direction == 'image_to_measurement':
                value = apply_axis_matrix(value, trans, 2)
                return apply_axis_matrix(value, axial, 1)
            if direction == 'measurement_to_image':
                value = apply_axis_matrix(value, axial.T, 1)
                return apply_axis_matrix(value, trans.T, 2)
            raise ValueError(f'Unknown type-6 resampling direction: {direction}')
        rows = int(getattr(self, 'measurement_nRowsD', self.nRowsD))
        cols = int(getattr(self, 'measurement_nColsD', self.nColsD))
        nx, nz = physical_grid(volume, grid)
        exact = exact_grid(volume, grid)
        physical_rows, physical_cols = physical_grid(volume, grid)
        if exact is not None and direction in ('measurement_to_image', 'image_to_measurement'):
            # The axes in exact are resampled directly and are at their target size afterwards; the others are padded/trimmed and
            # resized below as always
            value = exact_resample(value, volume, direction, exact, grid)
            if all(item is not None for item in exact):
                return value
            if exact[0] is not None:
                physical_rows = nx if direction == 'measurement_to_image' else rows
            if exact[1] is not None:
                physical_cols = nz if direction == 'measurement_to_image' else cols
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

    def resize_axis_matrix(source: int, target: int) -> np.ndarray:
        """Return ArrayFire resize's exact 1D linear map (target x source)."""
        basis = af.interop.np_to_af_array(np.asfortranarray(np.eye(source, dtype=np.float32)))
        resized = af.image.resize(basis, scale=None, odim0=target, odim1=source, method=af.INTERP.BILINEAR)
        return np.asarray(resized.to_ndarray(), dtype=np.float32)

    def center_matrix(source: int, target: int) -> np.ndarray:
        matrix = np.zeros((target, source), dtype=np.float32)
        length = min(source, target)
        src = (source - length) // 2
        dst = (target - length) // 2
        matrix[dst:dst + length, src:src + length] = np.eye(length, dtype=np.float32)
        return matrix

    def resample_adjoint(value: Any, volume: int, grid: tuple[int, int] | None = None) -> Any:
        """Transpose the actual image-to-measurement map used by FP."""
        if grid is not None:
            trans, axial = detector_axis_matrices(volume, grid)
            value = apply_axis_matrix(value, axial.T, 1)
            return apply_axis_matrix(value, trans.T, 2)
        from omegatomo.util.sinogramresample import axisMatrix

        rows = int(getattr(self, 'measurement_nRowsD', self.nRowsD))
        cols = int(getattr(self, 'measurement_nColsD', self.nColsD))
        physical_rows, physical_cols = physical_grid(volume, grid)
        image_rows, image_cols = (
            (int(grid[0]), int(grid[1])) if grid is not None
            else (int(self.Ny[volume]), int(self.Nz[volume]))
        )
        exact = exact_grid(volume, grid)
        if exact is not None:
            if exact[0] is None:
                row_map = center_matrix(physical_rows, rows) @ resize_axis_matrix(image_rows, physical_rows)
            else:
                row_map = axisMatrix(image_rows, rows, image_rows / exact[0]).astype(np.float32)
            if exact[1] is None:
                col_map = center_matrix(physical_cols, cols) @ resize_axis_matrix(image_cols, physical_cols)
            else:
                col_map = axisMatrix(image_cols, cols, image_cols / exact[1]).astype(np.float32)
        else:
            row_map = center_matrix(physical_rows, rows) @ resize_axis_matrix(image_rows, physical_rows)
            col_map = center_matrix(physical_cols, cols) @ resize_axis_matrix(image_cols, physical_cols)
        # axisMatrix in exact_grid is (transverse rows, axial columns).
        value = apply_axis_matrix(value, row_map.T, 2)
        value = apply_axis_matrix(value, col_map.T, 1)
        return value

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

    def attenuation_image(volume: int, timestep: int, subset: int, like: Any) -> Any:
        # 'vaimennus'/'d_attenuation_image' covers the union bounding box of every
        # multi-resolution volume (see _type6_attenuation_map_extent); resample it onto
        # `volume`'s own grid (identity when nMultiVolumes == 0) before handing it to
        # canonical_image exactly as before -- host_resource's Fortran ravel makes x the
        # fastest-varying axis, matching a plain 'C'-order reshape into (Nz, Ny, Nx) here.
        bxMap, byMap, bzMap, NxMap, NyMap, NzMap = _type6_attenuation_map_extent(self)
        flat = host_resource('d_attenuation_image', timestep, subset, getattr(self, 'vaimennus', np.empty(0)))
        if flat.size != NxMap * NyMap * NzMap:
            raise ValueError(
                f'Expected the type-6 attenuation map to cover the union of every multi-resolution '
                f'volume ({NxMap}x{NyMap}x{NzMap} = {NxMap * NyMap * NzMap} elements), got {flat.size}')
        attenMap = flat.reshape((NzMap, NyMap, NxMap), order='C')
        resampled = _type6_resample_attenuation_numpy(attenMap, self, volume, bxMap, byMap, bzMap, NxMap, NyMap, NzMap)
        return canonical_image(
            resampled.ravel(order='C'), (int(self.Nz[volume]), int(self.Ny[volume]), int(self.Nx[volume])))

    def attenuation_factors(volume: int, timestep: int, subset: int, view: int, like: Any) -> Any:
        """Build the native axis-centred attenuation ray factor for one view.

        The union attenuation source can be wider than any emission volume. Rotate it
        onto the volume's padded detector grid, shift it by the native panel offset,
        then integrate along x before cropping to the emission rotation grid.
        """
        bxMap, byMap, bzMap, NxMap, NyMap, NzMap = _type6_attenuation_map_extent(self)
        flat = host_resource('d_attenuation_image', timestep, subset,
                             getattr(self, 'vaimennus', np.empty(0)))
        if flat.size != NxMap * NyMap * NzMap:
            raise ValueError(
                f'Expected the type-6 attenuation map to cover the union of every multi-resolution '
                f'volume ({NxMap}x{NyMap}x{NzMap} = {NxMap * NyMap * NzMap} elements), got {flat.size}')
        atten_map = flat.reshape((NzMap, NyMap, NxMap), order='C')
        members, group_rows, group_bz, group_nz, att_nx, z_offsets = _type6_attenuation_group(self, volume)
        representative = members[0]
        source_key = (timestep, subset, id(getattr(self, 'vaimennus', None)),
                      tuple(members), group_bz, group_nz)
        cached_source = attenuation_source_cache.get(source_key)
        if cached_source is None:
            source, source_att_nx = _type6_attenuation_source_numpy(
                self, atten_map, representative, bz_override=group_bz, nz_override=group_nz)
            if source_att_nx != att_nx:
                raise RuntimeError('Type-6 attenuation group grid changed while building its source lattice')
            attenuation_source_cache[source_key] = (source, source_att_nx)
        else:
            source, source_att_nx = cached_source
        src_nz, src_ny, src_nx = source.shape
        grid_nx, grid_ny = _type6_rotation_grid(self, volume)
        dx = float(np.asarray(self.dx).reshape(-1)[volume])
        dy = float(np.asarray(self.dy).reshape(-1)[volume])
        angle = float(np.asarray(self.swivelAngles).reshape(-1)[view])
        radians = np.deg2rad(angle)
        cosine, sine = np.cos(radians), np.sin(radians)

        panel_values = np.asarray(getattr(self, 'blurPlanes2F', np.empty(0)))
        panel_shift_vox = (float(_type6_volume_view_value(self, 'blurPlanes2F', 0, view))
                           if panel_values.size else 0.)
        panel_shift = panel_shift_vox * float(np.asarray(self.dx).reshape(-1)[0]) / dx
        offset_x, offset_y = -sine * panel_shift, cosine * panel_shift

        # Match native/Torch affine_grid sampling with align_corners=False.
        # All members of an attenuation group share a source z lattice and
        # row count; each volume crops its own centered detector-grid rows.
        yy, xx = np.mgrid[:group_rows, :att_nx].astype(np.float64)
        xxn = 2. * (xx + .5) / att_nx - 1.
        yyn = 2. * (yy + .5) / group_rows - 1.
        sx = ((cosine * att_nx / src_nx) * xxn
              - (sine * dy * group_rows / (dx * src_nx)) * yyn
              + 2. * offset_x / src_nx + 1.) * src_nx * .5 - .5
        sy = ((sine * dx * att_nx / (dy * src_ny)) * xxn
              + (cosine * group_rows / src_ny) * yyn
              + 2. * offset_y / src_ny + 1.) * src_ny * .5 - .5
        x0, y0 = np.floor(sx).astype(np.int64), np.floor(sy).astype(np.int64)
        fx, fy = sx - x0, sy - y0
        rotated = np.zeros((src_nz, group_rows, att_nx), dtype=np.float32)
        for oy, wy in ((0, 1. - fy), (1, fy)):
            iy = y0 + oy
            for ox, wx in ((0, 1. - fx), (1, fx)):
                ix = x0 + ox
                valid = (iy >= 0) & (iy < src_ny) & (ix >= 0) & (ix < src_nx)
                if np.any(valid):
                    rotated += (source[:, np.clip(iy, 0, src_ny - 1),
                                        np.clip(ix, 0, src_nx - 1)]
                                * (wy * wx * valid).astype(np.float32)[None, :, :])
        cumulative = np.cumsum(rotated, axis=2, dtype=np.float32)
        import os
        try:
            attenuation_mode = int(os.environ.get('TYPE6_ATTCONV', '2')[:1])
        except (TypeError, ValueError):
            attenuation_mode = 2
        if attenuation_mode not in (0, 1, 2):
            attenuation_mode = 2
        if attenuation_mode == 0:
            factors_group = np.exp(-cumulative * dx)
        else:
            previous = np.zeros_like(cumulative)
            previous[:, :, 1:] = cumulative[:, :, :-1]
            if attenuation_mode == 1:
                factors_group = np.exp(-dx * (cumulative + previous) * .5)
            else:
                interval = dx * (cumulative - previous)
                average = np.ones_like(interval)
                nonzero = np.abs(interval) > 1e-6
                average[nonzero] = -np.expm1(-interval[nonzero]) / interval[nonzero]
                small = ~nonzero
                average[small] = 1. - interval[small] * .5 + interval[small] ** 2 / 6.
                factors_group = np.exp(-dx * previous) * average
        pad_lo, _ = _type6_padding_z(self, volume)
        row_off = (group_rows - grid_ny) // 2
        xoff = (att_nx - grid_nx) // 2
        z_group = z_offsets[volume] + np.arange(_type6_padding_z(self, volume)[0] +
                                                  int(self.Nz[volume]) +
                                                  _type6_padding_z(self, volume)[1]) - pad_lo
        y_group = row_off + np.arange(grid_ny)
        x_group = xoff + np.arange(grid_nx)
        valid_z = (z_group >= 0) & (z_group < src_nz)
        valid_y = (y_group >= 0) & (y_group < group_rows)
        valid_x = (x_group >= 0) & (x_group < att_nx)
        factors = np.ones((z_group.size, grid_ny, grid_nx), dtype=np.float32)
        if np.any(valid_z) and np.any(valid_y) and np.any(valid_x):
            factors[np.ix_(valid_z, valid_y, valid_x)] = factors_group[
                np.ix_(z_group[valid_z], y_group[valid_y], x_group[valid_x])]
        return canonical_image(factors.ravel(order='C'), factors.shape)

    def rotation_matrix(volume: int, angle: float, view: int):
        """Sparse bilinear forward rotation for one volume and view.

        This host sparse matrix is deliberately shared by FP and BP: BP uses
        its transpose, which preserves the exact zero-padded bilinear weights
        at clipped boundaries. ArrayFire's gather-only rotate cannot provide
        that transpose. The fallback is CPU-bound but remains exact.
        """
        from scipy.sparse import coo_matrix

        nx, ny = int(self.Nx[volume]), int(self.Ny[volume])
        out_nx, out_ny = _type6_rotation_grid(self, volume)
        dx = float(np.asarray(self.dx).reshape(-1)[volume])
        dy = float(np.asarray(self.dy).reshape(-1)[volume])
        ox, oy = _type6_axis_offset(self, volume)
        if np.asarray(getattr(self, 'blurPlanes2F', np.empty(0))).size:
            panel_shift = _type6_panel_shift(self, volume, view, out_ny)
            panel_angle = np.deg2rad(float(angle))
            ox += -float(np.sin(panel_angle)) * panel_shift
            oy += float(np.cos(panel_angle)) * panel_shift
        key = (volume, float(angle), nx, ny, out_nx, out_ny, dx, dy, ox, oy)
        cache = getattr(self, '_type6_af_rotation_cache', None)
        if cache is None:
            cache = {}
            self._type6_af_rotation_cache = cache
        if key in cache:
            matrix = cache.pop(key)
            cache[key] = matrix
            return matrix

        yy, xx = np.mgrid[:out_ny, :out_nx].astype(np.float64)
        xxn = (2. * (xx + .5) / out_nx) - 1.
        yyn = (2. * (yy + .5) / out_ny) - 1.
        # Native rotateHelperType6 receives -swivelAngles and samples
        # src_x = cos(a)*x - sin(a)*y, src_y = sin(a)*x + cos(a)*y.
        # This sparse matrix is the output-to-input sampling map, so its
        # angle must have the same sign as that native sampling map.
        rad = np.deg2rad(-float(angle))
        c, s = np.cos(rad), np.sin(rad)
        sx = ((c * out_nx / nx) * xxn + (s * dy * out_ny / (dx * nx)) * yyn + 2. * ox / nx + 1.) * nx * .5 - .5
        sy = ((-s * dx * out_nx / (dy * ny)) * xxn + (c * out_ny / ny) * yyn + 2. * oy / ny + 1.) * ny * .5 - .5
        x0, y0 = np.floor(sx).astype(np.int64), np.floor(sy).astype(np.int64)
        fx, fy = sx - x0, sy - y0
        out_row = (yy.astype(np.int64) * out_nx + xx.astype(np.int64)).ravel()
        rows, cols, vals = [], [], []
        for oy_i, wy in ((0, 1. - fy), (1, fy)):
            iy = y0 + oy_i
            for ox_i, wx in ((0, 1. - fx), (1, fx)):
                ix = x0 + ox_i
                weight = (wy * wx).ravel()
                valid = (iy >= 0) & (iy < ny) & (ix >= 0) & (ix < nx) & (weight.reshape(out_ny, out_nx) != 0.)
                if np.any(valid):
                    rows.append(out_row[valid.ravel()])
                    cols.append((iy * nx + ix).ravel()[valid.ravel()])
                    vals.append(weight[valid.ravel()].astype(np.float32, copy=False))
        if vals:
            matrix = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                                shape=(out_ny * out_nx, ny * nx), dtype=np.float32).tocsr()
        else:
            # A shifted volume can lie wholly outside this view's rotation grid.
            matrix = coo_matrix((out_ny * out_nx, ny * nx), dtype=np.float32).tocsr()
        cache[key] = matrix
        while len(cache) > 4:
            cache.pop(next(iter(cache)))
        return matrix

    def apply_rotation(value: Any, volume: int, angle: float, view: int, *, adjoint: bool) -> Any:
        nx, ny, nz = (int(self.Nx[volume]), int(self.Ny[volume]), int(self.Nz[volume]))
        out_nx, out_ny = _type6_rotation_grid(self, volume)
        pad_lo, pad_hi = _type6_padding_z(self, volume)
        out_nz = nz + pad_lo + pad_hi
        canonical = af.data.reorder(value, 2, 1, 0).to_ndarray()
        matrix = rotation_matrix(volume, angle, view)
        if adjoint:
            matrix = matrix.transpose().tocsr()
            src_x, src_y, dst_x, dst_y = out_nx, out_ny, nx, ny
            depth_offset = pad_lo
            output_z = nz
        else:
            src_x, src_y, dst_x, dst_y = nx, ny, out_nx, out_ny
            depth_offset = 0
            output_z = out_nz
        result = np.zeros((dst_x, dst_y, output_z), dtype=np.float32, order='F')
        for z in range(nz):
            plane = np.asarray(canonical[:, :, z + depth_offset]).T.reshape(-1, order='C') if adjoint else np.asarray(canonical[:, :, z]).T.reshape(-1, order='C')
            projected = matrix @ plane
            result[:, :, z + (0 if adjoint else pad_lo)] = projected.reshape((dst_y, dst_x), order='C').T
        af_value = af.interop.np_to_af_array(np.asfortranarray(result))
        return af.data.reorder(af_value, 2, 1, 0)

    def expand_kernel(kernel: Any, volume: int, nx: int) -> Any:
        local_nx = int(self.Nx[volume])
        if nx == local_nx:
            return kernel
        dx = float(np.asarray(self.dx).reshape(-1)[volume])
        bx = float(np.asarray(self.bx).reshape(-1)[volume])
        start = int(round((bx + .5 * nx * dx) / dx))
        values = kernel.to_ndarray()
        expanded = np.zeros((values.shape[0], values.shape[1], nx), dtype=np.float32, order='F')
        lo, hi = max(start, 0), min(start + values.shape[2], nx)
        if hi > lo:
            expanded[:, :, lo:hi] = values[:, :, lo - start:hi - start]
        return af.interop.np_to_af_array(np.asfortranarray(expanded))

    def blur(value: Any, kernel: Any, *, adjoint: bool = False) -> Any:
        # ArrayFire's convolve2 is mathematical convolution, while the common
        # type-6 forward operator uses cross-correlation. Flip for FP and use
        # the unflipped kernel for its exact transpose.
        if not adjoint:
            kernel = af.flip(af.flip(kernel, 0), 1)
        return af.data.reorder(
            af.signal.convolve2(af.data.reorder(value, 1, 0, 2), kernel), 1, 0, 2)

    def stack_views(values: list[Any]) -> Any:
        if not values:
            raise ValueError('Cannot stack an empty type-6 view list')
        if len(values) == 1:
            joined = af.data.moddims(values[0], int(values[0].shape[0]), int(values[0].shape[1]), 1)
        else:
            joined = values[0]
            for value in values[1:]:
                joined = af.join(2, joined, value)
        return af.data.reorder(joined, 2, 0, 1)

    def kernel(volume: int, timestep: int, subset: int, view: int, like: Any) -> Any:
        if bool(getattr(self, 'CDRFOnTheFly', False)):
            return af.interop.np_to_af_array(np.asfortranarray(_type6_on_the_fly_cdrf(self, volume, view)))
        filters = getattr(self, 'd_gFilter', None)
        if isinstance(filters, list) and volume < len(filters):
            return filters[volume]
        return af.interop.np_to_af_array(np.asfortranarray(_type6_volume_kernel(self, volume)))

    return SimpleNamespace(
        zeros=lambda shape, like: af.data.constant(0., *shape),
        reshape_image=lambda value, shape: af.data.reorder(af.data.moddims(value, shape[2], shape[1], shape[0]), 2, 1, 0),
        reshape_measurement=lambda value, shape: af.data.reorder(af.data.moddims(value, shape[2], shape[1], shape[0]), 2, 1, 0),
        flatten_image=lambda value: af.flat(af.data.reorder(value, 2, 1, 0)),
        flatten_measurement=lambda value: af.flat(af.data.reorder(value, 2, 1, 0)),
        stack_views=stack_views,
        rotation_grid=lambda volume: _type6_rotation_grid(self, volume),
        axial_grid=lambda volume: (int(self.Nz[volume]) + sum(_type6_padding_z(self, volume)[:2]), *_type6_padding_z(self, volume)),
        expand_kernel=expand_kernel,
        rotate_volume=lambda value, angle, volume, view: apply_rotation(value, volume, angle, view, adjoint=False),
        rotate_volume_adjoint=lambda value, angle, volume, view: apply_rotation(value, volume, angle, view, adjoint=True),
        rotate=lambda value, angle: af.data.reorder(
            af.image.rotate(af.data.reorder(value, 2, 1, 0), -angle * np.pi / 180.0, method=af.INTERP.BILINEAR), 2, 1, 0),
        attenuation=attenuation,
        shift_kernel=shift_kernel,
        blur=blur,
        resize=resize,
        shift_projection=shift_projection,
        shift_image_y=shift_image_y,
        sum_x=lambda value: af.data.moddims(af.algorithm.sum(value, 2), int(value.shape[0]), int(value.shape[1])),
        smear=lambda value, nx: af.tile(value, 1, 1, nx),
        set_view=lambda stack, index, value: stack.__setitem__((index, slice(None), slice(None)), value),
        add_to=lambda target, value: target.__setitem__(slice(None), target + value),
        multiply=lambda left, right: left * right,
        resample=lambda value, volume, direction, grid=None: resample(value, volume, direction, grid=grid),
        resample_adjoint=lambda value, volume, grid=None: resample_adjoint(value, volume, grid=grid),
        attenuation_image=attenuation_image,
        attenuation_factors=attenuation_factors,
        kernel=kernel,
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
    attenuation_source_cache = {}

    def shift_kernel(kernel: Any, plane: int, nx: int) -> Any:
        shifted = torch.zeros_like(kernel)
        depth = int(kernel.shape[2])
        if plane >= 0:
            if plane < depth:
                shifted[:, :, plane:] = kernel[:, :, :depth - plane]
        elif -plane < depth:
            shifted[:, :, :depth + plane] = kernel[:, :, -plane:]
        return shifted[:, :, :nx]

    def expand_kernel(kernel: Any, volume: int, nx: int) -> Any:
        local_nx = int(self.Nx[volume])
        if nx == local_nx:
            return kernel
        dx = float(np.asarray(self.dx).reshape(-1)[volume])
        bx = float(np.asarray(self.bx).reshape(-1)[volume])
        start = int(round((bx + .5 * nx * dx) / dx))
        result = torch.zeros((*kernel.shape[:2], nx), dtype=kernel.dtype, device=kernel.device)
        lo, hi = max(start, 0), min(start + int(kernel.shape[2]), nx)
        if hi > lo:
            result[:, :, lo:hi] = kernel[:, :, lo - start:hi - start]
        return result

    def blur(value: Any, kernel: Any, *, adjoint: bool = False) -> Any:
        kernel = kernel.permute(2, 0, 1).contiguous()
        if adjoint:
            kernel = kernel.flip((-2, -1))
        kernel = kernel.unsqueeze(1)
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

    def rotate_volume(value: Any, angle: float, volume: int, view: int) -> Any:
        grid_nx, grid_ny = _type6_rotation_grid(self, volume)
        nz, in_ny, in_nx = (int(v) for v in value.shape)
        dx = float(np.asarray(self.dx).reshape(-1)[volume])
        dy = float(np.asarray(self.dy).reshape(-1)[volume])
        ox, oy = _type6_axis_offset(self, volume)
        if np.asarray(getattr(self, 'blurPlanes2F', np.empty(0))).size:
            panel_shift = _type6_panel_shift(self, volume, view, grid_ny)
            panel_angle = float(angle) * np.pi / 180.
            ox += -np.sin(panel_angle) * panel_shift
            oy += np.cos(panel_angle) * panel_shift
        # affine_grid maps output coordinates back to the source. Match the
        # native type-6 sampler's src_x = c*x - s*y, src_y = s*x + c*y.
        radians = -float(angle) * np.pi / 180.
        cosine, sine = float(np.cos(radians)), float(np.sin(radians))
        theta = torch.tensor([[[
            cosine * grid_nx / in_nx, sine * dy * grid_ny / (dx * in_nx), ox * 2. / in_nx,
        ], [
            -sine * dx * grid_nx / (dy * in_ny), cosine * grid_ny / in_ny, oy * 2. / in_ny,
        ]]], dtype=value.dtype, device=value.device)
        grid = functional.affine_grid(theta.expand(nz, -1, -1), (nz, 1, grid_ny, grid_nx), align_corners=False)
        rotated = functional.grid_sample(value.unsqueeze(1), grid, mode='bilinear', padding_mode='zeros', align_corners=False).squeeze(1)
        pad_lo, pad_hi = _type6_padding_z(self, volume)
        return functional.pad(rotated, (0, 0, 0, 0, pad_lo, pad_hi))

    def rotate_volume_adjoint(value: Any, angle: float, volume: int, view: int) -> Any:
        # The inverse angle is not the transpose of bilinear sampling at clipped
        # boundaries. Autograd computes the exact scatter transpose of this affine.
        with torch.enable_grad():
            source = torch.zeros(
                (int(self.Nz[volume]), int(self.Ny[volume]), int(self.Nx[volume])),
                dtype=value.dtype, device=value.device, requires_grad=True,
            )
            transformed = rotate_volume(source, angle, volume, view)
            return torch.autograd.grad(transformed, source, grad_outputs=value,
                                       create_graph=False, retain_graph=False)[0]

    def resample_adjoint(value: Any, volume: int, grid: tuple[int, int]) -> Any:
        axial, trans = detector_matrices(value, volume, grid)
        return torch.matmul(torch.matmul(axial.transpose(0, 1), value), trans)

    def detector_matrices(like: Any, volume: int, grid: tuple[int, int]) -> tuple[Any, Any]:
        # Native mapType6/detectorWeightsType6: detector bin d samples lattice
        # coordinate A*d+B with linear tent weights and zero support outside.
        rows = int(getattr(self, 'measurement_nRowsD', self.nRowsD))
        cols = int(getattr(self, 'measurement_nColsD', self.nColsD))
        dx = np.asarray(self.dx, dtype=np.float64).reshape(-1)
        dy = np.asarray(self.dy, dtype=np.float64).reshape(-1)
        dz = np.asarray(self.dz, dtype=np.float64).reshape(-1)
        bz = np.asarray(self.bz, dtype=np.float64).reshape(-1)
        pad_lo, _ = _type6_padding_z(self, volume)
        axial_a = float(dz[0] / dz[volume])
        axial_b = -float(bz[volume]) / dz[volume] + (.5 - .5 * cols) * axial_a - .5 + pad_lo
        trans_a = float(dy[0] / dy[volume])
        trans_b = .5 * grid[0] + (.5 - .5 * rows) * trans_a - .5
        det_z = torch.arange(cols, dtype=like.dtype, device=like.device)[:, None]
        src_z = torch.arange(grid[1], dtype=like.dtype, device=like.device)[None, :]
        det_y = torch.arange(rows, dtype=like.dtype, device=like.device)[:, None]
        src_y = torch.arange(grid[0], dtype=like.dtype, device=like.device)[None, :]
        axial = torch.clamp(1. - torch.abs(det_z * axial_a + axial_b - src_z), min=0.)
        trans = torch.clamp(1. - torch.abs(det_y * trans_a + trans_b - src_y), min=0.)
        return axial, trans

    def resample_grid(value: Any, volume: int, direction: str, grid: tuple[int, int]) -> Any:
        axial, trans = detector_matrices(value, volume, grid)
        if direction == 'image_to_measurement':
            return torch.matmul(torch.matmul(axial, value), trans.transpose(0, 1))
        if direction == 'measurement_to_image':
            return torch.matmul(torch.matmul(axial.transpose(0, 1), value), trans)
        raise ValueError(f'Unknown type-6 resampling direction: {direction}')

    def axial_grid(volume: int) -> tuple[int, int, int]:
        lo, hi = _type6_padding_z(self, volume)
        return int(self.Nz[volume]) + lo + hi, lo, hi

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

    def attenuation_map(timestep: int, subset: int, like: Any) -> tuple[Any, tuple[float, float, float, int, int, int]]:
        extent = _type6_attenuation_map_extent(self)
        _, _, _, nx, ny, nz = extent
        raw = _type6_torch_resource(
            self, 'd_attenuation_image', timestep, subset,
            fallback=getattr(self, 'vaimennus', np.empty(0)), dtype=like.dtype, device=like.device,
        )
        frame_size = nx * ny * nz
        frames = int(getattr(self, 'Nt', 1))
        if raw.numel() == frame_size * frames and frames > 1:
            raw = raw[..., timestep] if raw.ndim == 4 else raw.reshape(-1)[timestep * frame_size:(timestep + 1) * frame_size]
        elif raw.numel() != frame_size:
            raise ValueError(f'Expected one or {frames} type-6 attenuation frames of {frame_size} elements, got {raw.numel()}')
        if raw.ndim >= 3:
            return raw.reshape(nx, ny, nz).permute(2, 1, 0).contiguous(), extent
        return raw.reshape(nz, ny, nx), extent

    def attenuation_image(volume: int, timestep: int, subset: int, like: Any) -> Any:
        mapTensor, (bxMap, byMap, bzMap, NxMap, NyMap, NzMap) = attenuation_map(timestep, subset, like)
        if volume == 0 and NxMap == int(self.Nx[0]) and NyMap == int(self.Ny[0]) and NzMap == int(self.Nz[0]):
            return mapTensor
        host_map = mapTensor.detach().to(device='cpu').to(dtype=torch.float64)
        resampled = _type6_resample_attenuation_numpy(
            host_map.numpy(), self, volume, bxMap, byMap, bzMap, NxMap, NyMap, NzMap)
        return torch.as_tensor(resampled, dtype=like.dtype, device=like.device)

    def attenuation_factors(volume: int, timestep: int, subset: int, view: int, like: Any) -> Any:
        attenMap, (bxMap, byMap, bzMap, NxMap, NyMap, NzMap) = attenuation_map(timestep, subset, like)
        members, group_rows, group_bz, group_nz, attNx, z_offsets = _type6_attenuation_group(self, volume)
        representative = members[0]
        cache_key = (int(timestep), int(subset), representative, str(like.device), like.dtype)
        source = attenuation_source_cache.get(cache_key)
        if source is None:
            host_map = attenMap.detach().to(device='cpu').to(dtype=torch.float64)
            source_np, source_attNx = _type6_attenuation_source_numpy(
                self, host_map.numpy(), representative,
                bz_override=group_bz, nz_override=group_nz)
            if source_attNx != attNx:
                raise RuntimeError('Type-6 attenuation group grid changed while building its source lattice')
            source = torch.as_tensor(source_np, dtype=like.dtype, device=like.device)
            attenuation_source_cache[cache_key] = source
        gridNx, gridNy = _type6_rotation_grid(self, volume)
        gridNz, _, _ = axial_grid(volume)
        srcNz, srcNy, srcNx = (int(v) for v in source.shape)
        dx, dy = float(self.dx[volume]), float(self.dy[volume])
        radians = float(self.swivelAngles[view]) * np.pi / 180.
        cosine, sine = float(np.cos(radians)), float(np.sin(radians))
        panel_values = np.asarray(getattr(self, 'blurPlanes2F', np.empty(0)))
        panel_shift_vox = (float(_type6_volume_view_value(self, 'blurPlanes2F', 0, view))
                           if panel_values.size else 0.)
        panel_shift = panel_shift_vox * float(self.dx[0]) / dx
        offset_x, offset_y = -sine * panel_shift, cosine * panel_shift
        theta = torch.tensor([[[(cosine * attNx / srcNx),
                               (-sine * dy * group_rows / (dx * srcNx)), 2. * offset_x / srcNx],
                              [(sine * dx * attNx / (dy * srcNy)),
                               (cosine * group_rows / srcNy), 2. * offset_y / srcNy]]], dtype=source.dtype, device=source.device)
        sample_grid = functional.affine_grid(
            theta.expand(srcNz, -1, -1), (srcNz, 1, group_rows, attNx), align_corners=False)
        rotated = functional.grid_sample(
            source.unsqueeze(1), sample_grid, mode='bilinear', padding_mode='zeros', align_corners=False
        ).squeeze(1)
        accumulated = torch.cumsum(rotated, dim=2)
        previous = torch.cat((torch.zeros_like(accumulated[:, :, :1]), accumulated[:, :, :-1]), dim=2)
        dx_path = float(self.dx[volume])
        import os
        try:
            attenuation_mode = int(os.environ.get('TYPE6_ATTCONV', '2')[:1])
        except (TypeError, ValueError):
            attenuation_mode = 2
        if attenuation_mode == 0:
            factors_group = torch.exp(-dx_path * accumulated)
        elif attenuation_mode == 1:
            factors_group = torch.exp(-0.5 * dx_path * (previous + accumulated))
        else:
            optical_depth = dx_path * (accumulated - previous)
            safe_depth = torch.where(torch.abs(optical_depth) < 1e-4, torch.ones_like(optical_depth), optical_depth)
            voxel_average = -torch.expm1(-safe_depth) / safe_depth
            small = torch.abs(optical_depth) < 1e-4
            voxel_average = torch.where(
                small, 1. - optical_depth * (0.5 - optical_depth * (1. / 6. - optical_depth / 24.)),
                voxel_average,
            )
            factors_group = torch.exp(-dx_path * previous) * voxel_average
        pad_lo, _ = _type6_padding_z(self, volume)
        row_off = (group_rows - gridNy) // 2
        xoff = (attNx - gridNx) // 2
        z_group = z_offsets[volume] + torch.arange(gridNz, device=like.device) - pad_lo
        y_group = row_off + torch.arange(gridNy, device=like.device)
        x_group = xoff + torch.arange(gridNx, device=like.device)
        valid_z = (z_group >= 0) & (z_group < srcNz)
        z_safe = z_group.clamp(0, max(srcNz - 1, 0))
        y_safe = y_group.clamp(0, max(group_rows - 1, 0))
        x_safe = x_group.clamp(0, max(attNx - 1, 0))
        factors = factors_group[z_safe[:, None, None], y_safe[None, :, None], x_safe[None, None, :]]
        valid = (valid_z[:, None, None]
                 & (y_group[None, :, None] >= 0) & (y_group[None, :, None] < group_rows)
                 & (x_group[None, None, :] >= 0) & (x_group[None, None, :] < attNx))
        return torch.where(valid, factors, torch.ones_like(factors))

    def kernel(volume: int, timestep: int, subset: int, view: int, like: Any) -> Any:
        if bool(getattr(self, 'CDRFOnTheFly', False)):
            return torch.as_tensor(_type6_on_the_fly_cdrf(self, volume, view),
                                   dtype=like.dtype, device=like.device)
        return (
            getattr(self, 'd_gFilter', [])[volume].to(device=like.device, dtype=like.dtype)
            if isinstance(getattr(self, 'd_gFilter', None), list)
            and volume < len(getattr(self, 'd_gFilter', []))
            else _type6_torch_resource(
                self, 'd_gFilter', timestep, subset,
                fallback=_type6_volume_kernel(self, volume), dtype=like.dtype, device=like.device,
            )
        )

    return SimpleNamespace(
        zeros=lambda shape, like: torch.zeros(shape, dtype=like.dtype, device=like.device),
        reshape_image=lambda value, shape: value.reshape(shape),
        reshape_measurement=lambda value, shape: value.reshape(shape),
        flatten_image=lambda value: value.reshape(-1),
        flatten_measurement=lambda value: value.reshape(-1),
        stack_views=lambda values: torch.stack(values, dim=0),
        rotation_grid=lambda volume: _type6_rotation_grid(self, volume),
        axial_grid=axial_grid,
        expand_kernel=expand_kernel,
        rotate=lambda value, angle: rotate(value, angle, interpolation=InterpolationMode.BILINEAR, fill=0.0),
        rotate_volume=rotate_volume,
        rotate_volume_adjoint=rotate_volume_adjoint,
        resample_adjoint=resample_adjoint,
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
        resample=lambda value, volume, direction, grid=None: (
            resample_grid(value, volume, direction, grid) if grid is not None
            else resample_resize_proj6(value, self, volume, direction=direction)
        ),
        attenuation_image=attenuation_image,
        attenuation_factors=attenuation_factors,
        kernel=kernel,
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
    grid_nx, grid_ny = ops.rotation_grid(volume)
    grid_nz, _, _ = ops.axial_grid(volume)
    source = ops.reshape_image(image, (nz, ny, nx))
    use_attenuation = bool(self.attenuation_correction and self.CTAttenuation)
    image_views = []
    for local_view in range(projections):
        view = start + local_view
        angle = float(self.swivelAngles[view])
        # Match native type-6: the continuous panel translation is folded into
        # the rotation sampling coordinates, scaled from the main-volume grid.
        rotated = ops.rotate_volume(source, angle, volume, view)
        if use_attenuation:
            rotated = rotated * ops.attenuation_factors(volume, timestep, subset, view, image)
        depth_shift = int(_type6_volume_view_value(self, 'blurPlanes', volume, view))
        kernel = ops.kernel(volume, timestep, subset, view, image)
        shifted_kernel = (ops.shift_kernel(kernel, 0, grid_nx)
                          if bool(getattr(self, 'CDRFOnTheFly', False))
                          else ops.expand_kernel(ops.shift_kernel(kernel, depth_shift, nx), volume, grid_nx))
        blurred = ops.blur(rotated, shifted_kernel)
        projection = ops.sum_x(blurred) * _type6_path_weight(self, volume, view, nx)
        if int(projection.shape[0]) != grid_nz or int(projection.shape[1]) != grid_ny:
            projection = ops.resize(projection, (grid_nz, grid_ny))
        image_views.append(projection)
        if (local_view + 1) % 16 == 0:
            ops.synchronize(ops.device(image))
    measurement = ops.resample(ops.stack_views(image_views), volume, 'image_to_measurement', grid=(grid_ny, grid_nz))
    measurement = ops.multiply(measurement, ops.weights(measurement, projections, timestep, subset))
    ops.add_to(output, ops.flatten_measurement(measurement))


def type6_backward(self: Any, y: Any, output: Any, volume: int, subset: int, timestep: int, *, ops: Any | None = None) -> None:
    """SPECT rotation projector shared logic."""
    if str(getattr(self, 'type6BPRotation', 'adjoint')).lower() != 'adjoint':
        raise ValueError(
            "The Python type-6 operator currently supports type6BPRotation='adjoint' only; "
            "the 'fast' gather is not implemented by this backend."
        )
    ops = _type6_torch_ops(self) if ops is None else ops
    projections = int(self.nProjSubset[timestep][subset])
    rows = int(getattr(self, 'measurement_nRowsD', self.nRowsD))
    cols = int(getattr(self, 'measurement_nColsD', self.nColsD))
    data = ops.reshape_measurement(y, (projections, cols, rows))
    data = ops.multiply(data, ops.weights(data, projections, timestep, subset))
    start = int(np.sum(self.nProjSubset[:timestep, :]) + np.sum(self.nProjSubset[timestep, :subset]))
    nx, ny, nz = int(self.Nx[volume]), int(self.Ny[volume]), int(self.Nz[volume])
    grid_nx, grid_ny = ops.rotation_grid(volume)
    grid_nz, _, _ = ops.axial_grid(volume)
    data = ops.resample_adjoint(data, volume, grid=(grid_ny, grid_nz))
    use_attenuation = bool(self.attenuation_correction and self.CTAttenuation)
    image = ops.zeros((nz, ny, nx), y)
    for local_view in range(projections):
        view = start + local_view
        angle = float(self.swivelAngles[view])
        projection = data[local_view]
        if int(projection.shape[0]) != grid_nz or int(projection.shape[1]) != grid_ny:
            projection = ops.resize(projection, (grid_nz, grid_ny))
        depth_shift = int(_type6_volume_view_value(self, 'blurPlanes', volume, view))
        kernel = ops.kernel(volume, timestep, subset, view, y)
        smeared = ops.smear(projection, grid_nx) * _type6_path_weight(self, volume, view, nx)
        # The forward projection attenuates before blurring (see type6_forward and the C++
        # reference backProjectionType6Angle in functions.hpp), and the two steps do not
        # commute, so the adjoint has to blur first and attenuate second, using the same
        # rotated/shifted attenuation map as before.
        shifted_kernel = (ops.shift_kernel(kernel, 0, grid_nx)
                          if bool(getattr(self, 'CDRFOnTheFly', False))
                          else ops.expand_kernel(ops.shift_kernel(kernel, depth_shift, nx), volume, grid_nx))
        blurred = ops.blur(smeared, shifted_kernel, adjoint=True)
        if use_attenuation:
            blurred = blurred * ops.attenuation_factors(volume, timestep, subset, view, y)
        ops.add_to(image, ops.rotate_volume_adjoint(blurred, angle, volume, view))
        if (local_view + 1) % 16 == 0:
            ops.synchronize(ops.device(y))
    if self.useMaskBP:
        image = ops.apply_bp_mask(image, volume, timestep, subset)
    ops.add_to(output, ops.flatten_image(image))
