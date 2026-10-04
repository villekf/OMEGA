# -*- coding: utf-8 -*-

import numpy as np

# Kernel source (general_opencl_functions.h + auxKernels.cl) is identical for
# every call in a process; read it from disk only once.
_KERNEL_SOURCE_CACHE = {}

# Compiled program/function cache, keyed by (rType, context identity, bOpt
# tuple, kernel name). Recompiling on every call is the main cost of these
# standalone prior functions when used inside an iterative reconstruction
# loop, so the compiled kernel is reused whenever the same backend context,
# build options (which also encode e.g. NLType/SW/PW/useAdaptive/useRef/
# Lange defines) and kernel name are requested again.
_PROGRAM_CACHE = {}


def _kernel_source():
    if 'lines' not in _KERNEL_SOURCE_CACHE:
        from omegatomo.util.paths import opencl_header_dir
        headerDir = opencl_header_dir()
        with open(headerDir + 'general_opencl_functions.h', encoding="utf8") as f:
            hlines = f.read()
        with open(headerDir + 'auxKernels.cl', encoding="utf8") as f:
            lines = f.read()
        _KERNEL_SOURCE_CACHE['lines'] = hlines + lines
    return _KERNEL_SOURCE_CACHE['lines']


def _backend_setup(rType, clctx, queue):
    """Imports the backend modules required for rType, derives (clctx, queue)
    and the base build options, and reports whether this is a HIP CuPy
    context (which the standalone prior functions cannot support, since they
    require CUDA texture support)."""
    isHip = False
    if rType == 0:
        import pyopencl as cl
        import arrayfire as af
        ctx = af.opencl.get_context(retain=True)
        clctx = cl.Context.from_int_ptr(ctx)
        q = af.opencl.get_queue(True)
        queue = cl.CommandQueue.from_int_ptr(q)
        bOpt = ('-cl-single-precision-constant', '-DOPENCL', '-DCAST=float',)
    elif rType == 1 or rType == 2:
        import cupy as cp
        isHip = bool(getattr(cp.cuda.runtime, 'is_hip', False))
        bOpt = ('-DHIP', '-DPYTHON',) if isHip else ('-DCUDA', '-DPYTHON',)
    elif rType == 3:
        bOpt = ('-cl-single-precision-constant', '-DOPENCL', '-DCAST=float',)
    else:
        raise ValueError('Unsupported rType: %r' % (rType,))
    return clctx, queue, bOpt, isHip


def _global_size(Nx, Ny, Nz):
    localSize = (16, 16, 1)
    apu = [Nx % localSize[0], Ny % localSize[1], 0]
    erotus = [0] * 2
    if apu[0] > 0:
        erotus[0] = localSize[0] - apu[0]
    if apu[1] > 0:
        erotus[1] = localSize[1] - apu[1]
    globalSize = [Nx + erotus[0], Ny + erotus[1], Nz]
    return globalSize, localSize


def _compiled_kernel(rType, name, bOpt, clctx=None):
    if rType == 0 or rType == 3:
        import pyopencl as cl
        ctxKey = int(clctx.int_ptr)
    else:
        import cupy as cp
        # No cheap, stable context handle is exposed by CuPy the way
        # PyOpenCL exposes clctx.int_ptr; the current device index is used as
        # the context-identity proxy instead (one context per device here).
        ctxKey = int(cp.cuda.runtime.getDevice())
    key = (rType, ctxKey, bOpt, name)
    knl = _PROGRAM_CACHE.get(key)
    if knl is not None:
        return knl
    lines = _kernel_source()
    if rType == 0 or rType == 3:
        prg = cl.Program(clctx, lines).build(bOpt)
        knl = getattr(prg, name)
    else:
        mod = cp.RawModule(code=lines, options=bOpt)
        knl = mod.get_function(name)
    _PROGRAM_CACHE[key] = knl
    return knl


def _output_array(rType, n, queue):
    """Allocates the zero-initialized output array and returns (f, fArg),
    where fArg is the ready-to-use kernel argument for f (the raw OpenCL
    memory object for ArrayFire/torch, or the array itself for plain
    PyOpenCL/CuPy)."""
    if rType == 0:
        import pyopencl as cl
        import arrayfire as af
        f = af.data.constant(0, n, dtype=af.Dtype.f32)
        fPtr = f.raw_ptr()
        fArg = cl.MemoryObject.from_int_ptr(fPtr)
    elif rType == 3:
        import pyopencl as cl
        f = cl.array.zeros(queue, n, dtype=cl.cltypes.float)
        fArg = f.data
    elif rType == 2:
        import cupy as cp
        import torch
        f = torch.zeros(n, dtype=torch.float32, device='cuda')
        fArg = cp.asarray(f)
    else:
        import cupy as cp
        f = cp.zeros(n, dtype=cp.float32)
        fArg = f
    return f, fArg


def _image_input(rType, im, Nx, Ny, Nz, clctx, queue):
    """Builds the backend image/texture for `im`, copies the data in, and
    returns the ready-to-use kernel argument. Handles the ArrayFire
    raw-pointer unlock (OpenCL) and the torch->CuPy conversion (CUDA/HIP)."""
    if rType == 0 or rType == 3:
        import pyopencl as cl
        imformat = cl.ImageFormat(cl.channel_order.A, cl.channel_type.FLOAT)
        mf = cl.mem_flags
        d_im = cl.Image(clctx, mf.READ_ONLY, imformat, shape=(Nx, Ny, Nz))
        if rType == 0:
            import arrayfire as af
            imPtr = im.raw_ptr()
            imD = cl.MemoryObject.from_int_ptr(imPtr)
            cl.enqueue_copy(queue, d_im, imD, offset=(0), origin=(0, 0, 0), region=(Nx, Ny, Nz))
            af.device.unlock_array(im)
        else:
            cl.enqueue_copy(queue, d_im, im.data, offset=(0), origin=(0, 0, 0), region=(Nx, Ny, Nz))
        return d_im
    else:
        import cupy as cp
        chl = cp.cuda.texture.ChannelFormatDescriptor(32, 0, 0, 0, cp.cuda.runtime.cudaChannelFormatKindFloat)
        array = cp.cuda.texture.CUDAarray(chl, Nx, Ny, Nz)
        if rType == 2:
            imD = cp.asarray(im)
            array.copy_from(imD.reshape((Nz, Ny, Nx)))
        else:
            array.copy_from(im.reshape((Nz, Ny, Nx)))
        res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
        tdes = cp.cuda.texture.TextureDescriptor(
            addressModes=(cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp, cp.cuda.runtime.cudaAddressModeClamp),
            filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=0)
        return cp.cuda.texture.TextureObject(res, tdes)


def _to_device_1d(rType, queue, arr_np):
    """Uploads a 1-D NumPy array and returns (deviceArray, kernelArg)."""
    if rType == 0 or rType == 3:
        import pyopencl as cl
        d = cl.array.to_device(queue, arr_np)
        return d, d.data
    else:
        import cupy as cp
        d = cp.asarray(arr_np)
        return d, d


def _dims_args(rType, Nx, Ny, Nz):
    """The Nxyz-dims kernel argument pair (used twice, back to back, by every
    kernel here): a single shared int3 for OpenCL, or six separate int32
    scalars for CUDA/HIP (matching each kernel's parameter list)."""
    if rType == 0 or rType == 3:
        import pyopencl as cl
        d = cl.cltypes.make_int3(Nx, Ny, Nz)
        return [d, d]
    else:
        import cupy as cp
        return [cp.int32(Nx), cp.int32(Ny), cp.int32(Nz), cp.int32(Nx), cp.int32(Ny), cp.int32(Nz)]


def _scalar(rType, value):
    if rType == 0 or rType == 3:
        import pyopencl as cl
        return cl.cltypes.float(value)
    else:
        import cupy as cp
        return cp.float32(value)


def _launch(rType, queue, knl, f, fArg, args, globalSize, localSize):
    """Sets fArg as kernel argument 0 followed by `args` in order, launches
    the kernel, synchronizes, and (ArrayFire only) unlocks the output array."""
    if rType == 0 or rType == 3:
        import pyopencl as cl
        knl.set_arg(0, fArg)
        for i, a in enumerate(args, start=1):
            knl.set_arg(i, a)
        cl.enqueue_nd_range_kernel(queue, knl, globalSize, localSize)
        queue.finish()
        if rType == 0:
            import arrayfire as af
            af.device.unlock_array(f)
    else:
        gridSize = (globalSize[0] // localSize[0], globalSize[1] // localSize[1], globalSize[2])
        blockSize = (localSize[0], localSize[1], 1)
        knl(gridSize, blockSize, (fArg,) + tuple(args))
        if rType == 2:
            import torch
            torch.cuda.synchronize()
    return f


def _rdp_corner_weights(Ndx, Ndy, Ndz, Nz, dx, dy, dz):
    """Distance-based neighborhood weights for the RDPCORNERS kernel variant.

    This is a direct port of reconstruction/prepass.py's computeWeights(GGMRF=False)
    followed by quadWeights(isEmpty=True), which is exactly what the full OMEGA
    pipeline computes into options.weights_quad when options.RDP and
    options.RDPIncludeCorners are both set (see prepass.py's "Compute the weights"
    block). Reusing that logic verbatim (rather than re-deriving it) keeps the
    standalone function's neighbor ordering consistent with the already-verified
    native RDPCORNERS pipeline. Returns a flat float32 array of length
    (2*Ndx+1)*(2*Ndy+1)*(2*Ndz+1) - 1 (the center/self entry removed), matching the
    CONSTANT float* weight argument of RDPKernel(RDPCORNERS).
    """
    distX, distY, distZ = float(dx), float(dy), float(dz)
    # Offset vectors, each running from +N down to -N (matches the element
    # order the original loop-based implementation produced). Ordering:
    # x varies fastest, then y, then z slowest (see prepass.py computeWeights,
    # GGMRF=False branch, which this mirrors).
    xr = np.arange(Ndx, -Ndx - 1, -1) * distX
    yr = np.arange(Ndy, -Ndy - 1, -1) * distY
    zr = np.arange(Ndz, -Ndz - 1, -1) * distZ
    Xg, Yg, Zg = np.meshgrid(xr, yr, zr, indexing='ij')
    if Ndz == 0 or Nz == 1:
        dist = np.sqrt(Xg ** 2 + Yg ** 2)
    else:
        dist = np.sqrt(Xg ** 2 + Yg ** 2 + Zg ** 2)
    weights = dist.flatten(order='F')
    with np.errstate(divide='ignore'):
        # The center (self) entry has distance 0 -> 1/0 = inf by design; it is
        # dropped below, matching computeWeights.m/computeWeights (prepass.py).
        weights = 1.0 / weights
    finite_sum = np.sum(weights[~np.isinf(weights)])
    weights_quad = weights / finite_sum
    half_len = weights_quad.size // 2
    weights_quad = np.concatenate((weights_quad[:half_len], weights_quad[half_len + 1:]))
    weights_quad = weights_quad[~np.isinf(weights_quad)]
    return weights_quad.astype(np.float32)


def RDP(im, Nx, Ny, Nz, gamma, beta, rType = 0, clctx = -1, queue = -1, includeCorners = False,
        Ndx = 1, Ndy = 1, Ndz = 1, dx = 1., dy = 1., dz = 1., weights = None):
    """
    Relative difference prior
    This is a standalone function for computing relative difference prior.
    Supports ArrayFire arrays, PyOpenCL arrays, CuPy arrays or PyTorch tensors
    as the input. Use rType to specify the input type.

    Args:
        im: The image where the regularization should be applied. This should be a
        vector and column-major!

        Nx/Ny/Nz: Number of voxels in x/y/z-direction. Each should be either a scalar
        or NumPy array

        gamma: The adjustable value for RDP. Scalar float. Lower values smooth the
        image, while larger make it sharper.

        beta: Regularization paramerer/hyperparameter. Scalar float.

        rType: Reconstruction type. 0 = ArrayFire (OpenCL), 1 = CuPy, 2 = PyTorch,
        3 = PyOpenCL. Default is 0.

        clctx: Only used by PyOpenCL, omit otherwise. The PyOpenCL context value.

        queue: Only used by PyOpenCL, omit otherwise. The PyOpenCL command queue value.

        includeCorners: If True, uses the full (2*Ndx+1)x(2*Ndy+1)x(2*Ndz+1) voxel
        neighborhood (corners included) instead of only the 6 face-adjacent voxels,
        matching options.RDPIncludeCorners in the full OMEGA pipeline. Default is False
        (unchanged behavior, 6-neighbor RDP).

        Ndx/Ndy/Ndz: Neighborhood size (number of voxels) in x/y/z-direction. Only used
        when includeCorners is True. Default is 1 for each (a 3x3x3 neighborhood).

        dx/dy/dz: Voxel size (physical pitch) in x/y/z-direction. Only used when
        includeCorners is True, to weight neighbors by inverse Euclidean distance.
        Default is 1. for each (isotropic voxels). Only the relative values matter.

        weights: Optional precomputed neighborhood weights (e.g. an existing
        options.weights_quad from the full pipeline) to use instead of computing them
        from Ndx/Ndy/Ndz/dx/dy/dz. Only used when includeCorners is True. Should be a
        flat NumPy array of length (2*Ndx+1)*(2*Ndy+1)*(2*Ndz+1) - 1 (center excluded).

    Returns:
        f: The gradient of the RDP. Vector of the same type as the input im.
    """
    if type(Nx) == np.ndarray:
        Nx = Nx.item()
    if type(Ny) == np.ndarray:
        Ny = Ny.item()
    if type(Nz) == np.ndarray:
        Nz = Nz.item()

    clctx, queue, bOpt, isHip = _backend_setup(rType, clctx, queue)
    if isHip:
        raise ValueError('RDP standalone function requires CUDA texture support, which CuPy does not provide on ROCm/HIP. Use forward/backward projector types 1-4 with the full OMEGA reconstruction pipeline instead, or use a non-ROCm CuPy build.')

    bOpt += ('-DRDP', '-DUSEIMAGES', '-DLOCAL_SIZE=16', '-DLOCAL_SIZE2=16',)
    if includeCorners:
        bOpt += ('-DRDPCORNERS', '-DSWINDOWX=' + str(Ndx), '-DSWINDOWY=' + str(Ndy), '-DSWINDOWZ=' + str(Ndz),)

    epps = 1e-8
    globalSize, localSize = _global_size(Nx, Ny, Nz)

    if rType == 0 or rType == 3:
        queue.finish()
    knl = _compiled_kernel(rType, 'RDPKernel', bOpt, clctx)
    f, fArg = _output_array(rType, Nx * Ny * Nz, queue)
    d_im = _image_input(rType, im, Nx, Ny, Nz, clctx, queue)
    args = [d_im] + _dims_args(rType, Nx, Ny, Nz) + [
        _scalar(rType, gamma), _scalar(rType, epps), _scalar(rType, beta)]
    if includeCorners:
        w = np.asarray(weights, dtype=np.float32) if weights is not None else _rdp_corner_weights(Ndx, Ndy, Ndz, Nz, dx, dy, dz)
        w = np.ascontiguousarray(w, dtype=np.float32)
        expected = (Ndx * 2 + 1) * (Ndy * 2 + 1) * (Ndz * 2 + 1) - 1
        if w.size != expected:
            raise ValueError('weights must contain (2*Ndx+1)*(2*Ndy+1)*(2*Ndz+1) - 1 = %d entries for the requested RDPCORNERS neighborhood, got %d' % (expected, w.size))
        _, weightArg = _to_device_1d(rType, queue, w)
        args.append(weightArg)
    f = _launch(rType, queue, knl, f, fArg, args, globalSize, localSize)
    return f


def NLReg(im, Nx, Ny, Nz, h, beta, SW = (1, 1, 1), PW = (1, 1, 1), rType = 0, clctx = -1, queue = -1, NLType = 0, STD = 1., gamma = 10., phi = 10., useAdaptive = False, adaptiveConstant = 5e-6,
       GGMRFpqc = (2., 1.5, 0.001), refIm = []):
    """
    Non-local regularization methods
    This is a standalone function for computing non-local regularization. Supported
    non-local methods are: non-local means (NLM), non-local TV (NLTV), non-local
    relative difference (NLRD), NLM filtering, non-local Lange (NLLange), NL filtering
    with Lange and non-local GGMRF. NLM is used by default.
    Supports ArrayFire arrays, PyOpenCL arrays, CuPy arrays or PyTorch tensors
    as the input. Use rType to specify the input type.

    Args:
        im: The image where the regularization should be applied. This should be a
        vector and column-major!

        Nx/Ny/Nz: Number of voxels in x/y/z-direction. Each should be either a scalar
        or NumPy array

        h: The filter parameter for the non-local methods. Higher values smooth
        the image, while lower values make it sharper.

        beta: Regularization paramerer/hyperparameter. Scalar float.

        rType: Reconstruction type. 0 = ArrayFire (OpenCL), 1 = CuPy, 2 = PyTorch,
        3 = PyOpenCL. Default is 0.

        NLType: The regularization type. 0 = NLM, 1 = NLTV, 2 = NLM filtered, 3 =
        NLRD, 4 = NL Lange, 5 = NL filtered with Lange, 6 = NLGGMRF, and
        7 = NL Geman-McClure.

        SW: The search window (neighborhood) size. A tuple that that contains the number
        of voxels included for each dimension. Default is (1, 1, 1) which corresponds
        to a search window of size 3x3x3. The dimension is thus always * 2 + 1.

        PW: The patch window size. Otherwise identical to SW in function. It is
        recommended to keep this small. Default is (1, 1, 1).

        STD: The standard deviation for the Gaussian weighted Euclidian distance.
        This is used to weight the patch values based on the distance. Default is 1.
        Higher values give more emphasis to the voxels further from the center,
        while smaller values emphasize the center voxel.

        gamma: The adjustable value for NLRD. Scalar float. Only required by NLRD.
        Lower values smooth the image, while larger make it sharper. Also used as the
        delta of NL Geman-McClure (NLType 7), where differences much larger than it
        are ignored by the prior.

        phi: The adjustable value for NLLange. Scalar float. Only required by NLLange.

        useAdaptive: Use the adaptive weighting, based on the mean of the patch.
        Default is False. If you use this, you should also input adaptiveConstant.

        adaptiveConstant: Used by the adaptive weighting. This is the additive part
        of the adaptive method. While h controls the strength of the mean part,
        this value affects the whole image and thus large values can completely
        overwrite the effect of the adaptive weighting. Scalar float

        GGMRFpqc: The p, q, and c values for NLGGMRF. Only needed if NLGGMRF is selected.
        The input is a tuple that should contain all three values. p is first, q second,
        and c last.

        refIm: Optional reference image used in the patch window computations. This
        has to be in the same format and have the same size as im.

        clctx: Only used by PyOpenCL, omit otherwise. The PyOpenCL context value.

        queue: Only used by PyOpenCL, omit otherwise. The PyOpenCL command queue value.

    Returns:
        f: The gradient of the selected NL regularization. Vector of the same type as the
        input im.
    """
    if type(Nx) == np.ndarray:
        Nx = Nx.item()
    if type(Ny) == np.ndarray:
        Ny = Ny.item()
    if type(Nz) == np.ndarray:
        Nz = Nz.item()

    clctx, queue, bOpt, isHip = _backend_setup(rType, clctx, queue)
    if isHip:
        raise ValueError('NLReg standalone function requires CUDA texture support, which CuPy does not provide on ROCm/HIP. Use forward/backward projector types 1-4 with the full OMEGA reconstruction pipeline instead, or use a non-ROCm CuPy build.')

    bOpt += ('-DNLM_', '-DUSEIMAGES', '-DLOCAL_SIZE=16', '-DLOCAL_SIZE2=16', '-DNLTYPE=' + str(NLType), '-DSWINDOWX=' + str(SW[0]),
             '-DSWINDOWY=' + str(SW[1]), '-DSWINDOWZ=' + str(SW[2]), '-DPWINDOWX=' + str(PW[0]), '-DPWINDOWY=' + str(PW[1]),
             '-DPWINDOWZ=' + str(PW[2]),)
    if (useAdaptive):
        bOpt += ('-DNLMADAPTIVE',)
    useRef = type(im) == type(refIm)
    if useRef:
        bOpt += ('-DNLMREF',)

    epps = 1e-8
    globalSize, localSize = _global_size(Nx, Ny, Nz)

    x = np.linspace(-PW[0], PW[0], 2 * PW[0] + 1, dtype=np.float32)
    y = np.linspace(-PW[1], PW[1], 2 * PW[1] + 1, dtype=np.float32)
    z = np.linspace(-PW[2], PW[2], 2 * PW[2] + 1, dtype=np.float32)
    gaussK = np.exp(-(np.add.outer(np.add.outer(x**2 / (2*STD**2), y**2 / (2*STD**2)), z**2 / (2*STD**2))))
    gaussK = gaussK.flatten('F').astype(dtype=np.float32)
    if NLType == 4 or NLType == 5:
        gamma = phi

    if rType == 0 or rType == 3:
        queue.finish()
    knl = _compiled_kernel(rType, 'NLM', bOpt, clctx)
    f, fArg = _output_array(rType, Nx * Ny * Nz, queue)
    _, gaussArg = _to_device_1d(rType, queue, gaussK)
    d_im = _image_input(rType, im, Nx, Ny, Nz, clctx, queue)
    args = [d_im, gaussArg] + _dims_args(rType, Nx, Ny, Nz) + [
        _scalar(rType, h * h), _scalar(rType, epps), _scalar(rType, beta)]
    if NLType >= 3:
        args.append(_scalar(rType, gamma))
    if NLType == 6:
        args.append(_scalar(rType, GGMRFpqc[0]))
        args.append(_scalar(rType, GGMRFpqc[1]))
        args.append(_scalar(rType, GGMRFpqc[2]))
    if useAdaptive:
        args.append(_scalar(rType, adaptiveConstant))
    if useRef:
        d_refIm = _image_input(rType, refIm, Nx, Ny, Nz, clctx, queue)
        args.append(d_refIm)
    f = _launch(rType, queue, knl, f, fArg, args, globalSize, localSize)
    return f

def TV(im, Nx, Ny, Nz, beta, sValue = 1e-4, rType = 0, clctx = -1, queue = -1, Lange = False, sigma = 10.):
    """
    Total variation prior
    This is a standalone function for computing total variation prior. This is the
    gradient version of the prior and is thus not differentiable without additional
    "smoothing" parameter.
    Supports ArrayFire arrays, PyOpenCL arrays, CuPy arrays or PyTorch tensors
    as the input. Use rType to specify the input type.

    Args:
        im: The image where the regularization should be applied. This should be a
        vector and column-major!

        Nx/Ny/Nz: Number of voxels in x/y/z-direction. Each should be either a scalar
        or NumPy array

        beta: Regularization paramerer/hyperparameter. Scalar float.

        sValue: The smoothing value to guarantee differentiability. Default value is
        1e-4. Scalar float.

        rType: Reconstruction type. 0 = ArrayFire (OpenCL), 1 = CuPy, 2 = PyTorch,
        3 = PyOpenCL. Default is 0.

        clctx: Only used by PyOpenCL, omit otherwise. The PyOpenCL context value.

        queue: Only used by PyOpenCL, omit otherwise. The PyOpenCL command queue value.

        Lange: If True, computes the Lange prior instead of TV. Default is False.

        sigma: Adjustable parameter for the Lange prior. Scalar float. Default value
        is 10. Only used with Lange prior.

    Returns:
        f: The gradient of the TV prior. Vector of the same type as the input im.
    """
    if type(Nx) == np.ndarray:
        Nx = Nx.item()
    if type(Ny) == np.ndarray:
        Ny = Ny.item()
    if type(Nz) == np.ndarray:
        Nz = Nz.item()

    clctx, queue, bOpt, isHip = _backend_setup(rType, clctx, queue)
    if isHip:
        raise ValueError('TV standalone function requires CUDA texture support, which CuPy does not provide on ROCm/HIP. Use forward/backward projector types 1-4 with the full OMEGA reconstruction pipeline instead, or use a non-ROCm CuPy build.')

    bOpt += ('-DTVGRAD', '-DUSEIMAGES', '-DLOCAL_SIZE=16', '-DLOCAL_SIZE2=16',)
    if Lange:
        bOpt += ('-DSATV',)

    globalSize, localSize = _global_size(Nx, Ny, Nz)

    if rType == 0 or rType == 3:
        queue.finish()
    knl = _compiled_kernel(rType, 'TVKernel', bOpt, clctx)
    f, fArg = _output_array(rType, Nx * Ny * Nz, queue)
    d_im = _image_input(rType, im, Nx, Ny, Nz, clctx, queue)
    args = [d_im] + _dims_args(rType, Nx, Ny, Nz) + [
        _scalar(rType, sigma), _scalar(rType, sValue), _scalar(rType, beta)]
    f = _launch(rType, queue, knl, f, fArg, args, globalSize, localSize)
    return f
