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
    if self.subsetType > 7 or self.subsets == 1:
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
    args.img(f_arg)
    args.img(y_arg)
    if self.SPECT:
        args.buf(self.d_detectorVector[timestep][subset])
    args.u8(self.no_norm)
    args.u64(self.nMeasSubset[timestep, subset].item())
    args.u32(subset)
    args.i32(k)
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
    if self.attenuation_correction and not self.CTAttenuation:
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
    args.img(y_arg)
    args.img(f_arg)
    if self.SPECT:
        args.buf(self.d_detectorVector[timestep][subset])
    args.u8(self.no_norm)
    args.u64(self.nMeasSubset[timestep, subset].item())
    args.u32(subset)
    args.i32(k)
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
    if self.normalization_correction:
        args.buf(self.d_norm[timestep][subset])
    if self.additionalCorrection:
        args.buf(self.d_corr[timestep][subset])
    args.buf(self.d_Sens)
    return args


def conv3D(self, f, ii = 0):
    if getattr(self, "useMetal", False):
        raise NotImplementedError("The separate PSF convolution kernel is not yet wired to the Metal/MPS bridge.")
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
        _validate_forward_input(self, f)
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
                        if (self.attenuation_correction and not self.CTAttenuation):
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
            from pyopencl.version import VERSION
            if not self.loadTOF:
                if self.useIndexBasedReconstruction and self.listmode > 0:
                    self.d_trIndex[0] = cl.array.to_device(self.queue, self.trIndex[self.nMeas[timestep * self.subsets + subset] * 2 : self.nMeas[timestep * self.subsets + subset + 1] * 2])
                    self.d_axIndex[0] = cl.array.to_device(self.queue, self.axIndex[self.nMeas[timestep * self.subsets + subset] * 2 : self.nMeas[timestep * self.subsets + subset + 1] * 2])
                elif self.listmode > 0:
                    apu = self.x.ravel()
                    self.d_x[timestep][0] = cl.array.to_device(self.queue, apu[self.nMeas[timestep * self.subsets + subset] * 6 : self.nMeas[timestep * self.subsets + subset + 1] * 6])
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
                    if (self.attenuation_correction and not self.CTAttenuation):
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
                        f_arg = f.data
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
        _validate_backward_input(self, y, timestep, subset)
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


# NOTE: multi-resolution volumes (nMultiVolumes > 0) are not yet supported/validated for
# projector type 6 on this branch, independent of the attenuation resampling below -- type 6's
# rotation step rotates each volume about its own local grid centre and never reads bx/by, so
# laterally-offset (transaxial) side volumes are misplaced, and the side-volume normalization
# (~(1/multiResolutionScale)^2) is not accounted for; the full fix lives in a separate branch.
# Keep this attenuation resampling regardless (it works in physical coordinates, so it becomes
# correct once volume placement is), but do not treat multi-volume type-6 results as validated.


def _type6_attenuation_map_extent(self: Any) -> tuple[float, float, float, int, int, int]:
    """Union bounding box, at volume 0's own (finest) voxel size, of every type-6
    multi-resolution volume (0..nMultiVolumes) -- the grid the host-provided attenuation
    map ('vaimennus') is expected to cover, exactly like source/cpp/functions.hpp's
    buildType6AttenuationVolume. Positions come from bx/by/bz as already stored (NOT
    assumed centred on the origin), so this keeps working once a shifted main volume is
    supported. With nMultiVolumes == 0 this reduces to volume 0's own grid exactly."""
    nVol = int(self.nMultiVolumes) + 1
    bx = np.asarray(self.bx, dtype=np.float64).reshape(-1)[:nVol]
    by = np.asarray(self.by, dtype=np.float64).reshape(-1)[:nVol]
    bz = np.asarray(self.bz, dtype=np.float64).reshape(-1)[:nVol]
    dx = np.asarray(self.dx, dtype=np.float64).reshape(-1)[:nVol]
    dy = np.asarray(self.dy, dtype=np.float64).reshape(-1)[:nVol]
    dz = np.asarray(self.dz, dtype=np.float64).reshape(-1)[:nVol]
    Nx = np.asarray(self.Nx, dtype=np.float64).reshape(-1)[:nVol]
    Ny = np.asarray(self.Ny, dtype=np.float64).reshape(-1)[:nVol]
    Nz = np.asarray(self.Nz, dtype=np.float64).reshape(-1)[:nVol]
    bxMap, byMap, bzMap = float(np.min(bx)), float(np.min(by)), float(np.min(bz))
    bxHi = float(np.max(bx + Nx * dx))
    byHi = float(np.max(by + Ny * dy))
    bzHi = float(np.max(bz + Nz * dz))
    NxMap = int(round((bxHi - bxMap) / float(dx[0])))
    NyMap = int(round((byHi - byMap) / float(dy[0])))
    NzMap = int(round((bzHi - bzMap) / float(dz[0])))
    return bxMap, byMap, bzMap, NxMap, NyMap, NzMap


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
    if volume == 0 and nx == NxMap and ny == NyMap and nz == NzMap:
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
    volu = np.maximum(dZe * dYe * dXe, 1e-6)
    return (boxSum / volu).astype(np.float32)


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

    _weights_unset = object()  # sentinel distinguishing "attribute missing" from any real id()

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
        attenuation_image=attenuation_image,
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

    def attenuation_image(volume: int, timestep: int, subset: int, like: Any) -> Any:
        # 'vaimennus'/'d_attenuation_image' covers the union bounding box of every
        # multi-resolution volume (see _type6_attenuation_map_extent); resample it onto
        # `volume`'s own grid (identity when nMultiVolumes == 0, and bit-identical to the
        # pre-existing single-reshape behavior in that case) before returning, so the outer
        # reshape_image(..., (nz, ny, nx)) call downstream is left a no-op as before.
        bxMap, byMap, bzMap, NxMap, NyMap, NzMap = _type6_attenuation_map_extent(self)
        raw = _type6_torch_resource(
            self, 'd_attenuation_image', timestep, subset,
            fallback=getattr(self, 'vaimennus', np.empty(0)), dtype=like.dtype, device=like.device,
        )
        if raw.numel() != NxMap * NyMap * NzMap:
            raise ValueError(
                f'Expected the type-6 attenuation map to cover the union of every multi-resolution '
                f'volume ({NxMap}x{NyMap}x{NzMap} = {NxMap * NyMap * NzMap} elements), got {raw.numel()}')
        if raw.dim() >= 3:
            # A plain, not-pre-flattened (Nx, Ny, Nz)-shaped array (this codebase's usual
            # convention for a directly-assigned vaimennus, e.g. np.asfortranarray(mu) with
            # mu.shape == (Nx, Ny, Nz)): permute into (Nz, Ny, Nx) explicitly. A same-total
            # reshape() alone would silently mis-interpret the axes whenever the map is not
            # cubic (NxMap != NzMap), which the union bounding box normally is not.
            mapTensor = raw.reshape(NxMap, NyMap, NzMap).permute(2, 1, 0).contiguous()
        else:
            # Flat (1-D): x fastest-varying, matching this file's other host_resource /
            # canonical_image Fortran-ravel convention, so a plain C-order reshape into
            # (Nz, Ny, Nx) makes Nx (the last axis) the fastest-varying, as required.
            mapTensor = raw.reshape(NzMap, NyMap, NxMap)
        if volume == 0 and NxMap == int(self.Nx[0]) and NyMap == int(self.Ny[0]) and NzMap == int(self.Nz[0]):
            return mapTensor
        resampled = _type6_resample_attenuation_numpy(
            mapTensor.detach().to('cpu', dtype=torch.float64).numpy(), self, volume, bxMap, byMap, bzMap, NxMap, NyMap, NzMap)
        return torch.as_tensor(resampled, dtype=like.dtype, device=like.device)

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
        attenuation_image=attenuation_image,
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
    if self.attenuation_correction and self.CTAttenuation:
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
    if self.attenuation_correction and self.CTAttenuation:
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
        # The forward projection attenuates before blurring (see type6_forward and the C++
        # reference backProjectionType6Angle in functions.hpp), and the two steps do not
        # commute, so the adjoint has to blur first and attenuate second, using the same
        # rotated/shifted attenuation map as before.
        blurred = ops.blur(smeared, ops.shift_kernel(kernel, depth_shift, nx))
        if attenuation is not None:
            attenuation_rotated = ops.rotate(attenuation, angle)
            attenuation_rotated = ops.shift_image_y(
                attenuation_rotated, -int(_type6_volume_view_value(self, 'blurPlanes2', volume, view))
            )
            blurred = ops.attenuation(blurred, attenuation_rotated, float(self.dx[volume]))
        rotated = ops.shift_image_y(blurred, int(_type6_volume_view_value(self, 'blurPlanes2', volume, view)))
        ops.add_to(image, ops.rotate(rotated, -angle))
        if (local_view + 1) % 16 == 0:
            ops.synchronize(ops.device(y))
    if self.useMaskBP:
        image = ops.apply_bp_mask(image, volume, timestep, subset)
    ops.add_to(output, ops.flatten_image(image))
