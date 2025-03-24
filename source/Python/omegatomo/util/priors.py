# -*- coding: utf-8 -*-

def RDP(im, Nx, Ny, Nz, gamma, beta, rType = 0, clctx = -1, queue = -1):
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
    
    Returns:
        f: The gradient of the RDP. Vector of the same type as the input im.
    """
    if rType == 0:
        import pyopencl as cl
        import arrayfire as af
    elif rType == 1:
        import cupy as cp
    elif rType == 2:
        import cupy as cp
        import torch
    elif rType == 3:
        import pyopencl as cl
    import os
    import numpy as np
    
    if type(Nx) == np.ndarray:
        Nx = Nx.item()
    if type(Ny) == np.ndarray:
        Ny = Ny.item()
    if type(Nz) == np.ndarray:
        Nz = Nz.item()
        
    if rType == 0:
        ctx = af.opencl.get_context(retain=True)
        clctx = cl.Context.from_int_ptr(ctx)
        q = af.opencl.get_queue(True)
        queue = cl.CommandQueue.from_int_ptr(q)
        bOpt = ('-cl-single-precision-constant', '-DOPENCL', '-DCAST=float',)
    elif rType == 1 or rType == 2:
        bOpt = ('-DCUDA', '-DPYTHON',)
    elif rType == 3:
        bOpt = ('-cl-single-precision-constant', '-DOPENCL', '-DCAST=float',)
    
    bOpt += ('-DRDP', '-DUSEIMAGES', '-DLOCAL_SIZE=16', '-DLOCAL_SIZE2=16',)
    
    epps = 1e-8
    localSize = (16, 16, 1)
    apu = [Nx % localSize[0], Ny % localSize[1], 0]
    erotus = [0] * 2
    if apu[0] > 0:
        erotus[0] = localSize[0] - apu[0]
    if apu[1] > 0:
        erotus[1] = localSize[1] - apu[1]
    
    globalSize = [None] * 3
    globalSize[0] = Nx + erotus[0]
    globalSize[1] = Ny + erotus[1]
    globalSize[2] = Nz
    
    headerDir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', '..', '..', 'opencl')) + "/"
    with open(headerDir + 'general_opencl_functions.h', encoding="utf8") as f:
        hlines = f.read()
    with open(headerDir + 'auxKernels.cl', encoding="utf8") as f:
        lines = f.read()
    lines = hlines + lines
    if rType == 0 or rType == 3:
        queue.finish()
        d_Nxyz = cl.cltypes.make_int3(Nx, Ny, Nz)
        prg = cl.Program(clctx, lines).build(bOpt)
        knlRDP = prg.RDPKernel
        if rType == 0:
            imPtr = im.raw_ptr()
            imD = cl.MemoryObject.from_int_ptr(imPtr)
            f = af.data.constant(0, Nx * Ny * Nz, dtype=af.Dtype.f32)
            fPtr = f.raw_ptr()
            fD = cl.MemoryObject.from_int_ptr(fPtr)
        else:
            f = cl.array.zeros(queue, Nx * Ny * Nz, dtype=cl.cltypes.float)
        imformat = cl.ImageFormat(cl.channel_order.A, cl.channel_type.FLOAT)
        mf = cl.mem_flags
        d_im = cl.Image(clctx, mf.READ_ONLY, imformat, shape=(Nx, Ny, Nz))
        if rType == 0:
            cl.enqueue_copy(queue, d_im, imD, offset=(0), origin=(0,0,0), region=(Nx, Ny, Nz))
            af.device.unlock_array(im)
        else:
            cl.enqueue_copy(queue, d_im, im.data, offset=(0), origin=(0,0,0), region=(Nx, Ny, Nz))
        kIndLoc = 0
        if rType == 0:
            knlRDP.set_arg(kIndLoc, fD)
        else:
            knlRDP.set_arg(kIndLoc, f.data)
        kIndLoc += 1
        knlRDP.set_arg(kIndLoc, d_im)
        kIndLoc += 1
        knlRDP.set_arg(kIndLoc, d_Nxyz)
        kIndLoc += 1
        knlRDP.set_arg(kIndLoc, d_Nxyz)
        kIndLoc += 1
        knlRDP.set_arg(kIndLoc, (cl.cltypes.float)(gamma))
        kIndLoc += 1
        knlRDP.set_arg(kIndLoc, (cl.cltypes.float)(epps))
        kIndLoc += 1
        knlRDP.set_arg(kIndLoc, (cl.cltypes.float)(beta))
        cl.enqueue_nd_range_kernel(queue, knlRDP, globalSize, localSize)
        queue.finish()
        if rType == 0:
            af.device.unlock_array(f)
    else:
        mod = cp.RawModule(code=lines, options=bOpt)
        knlRDP = mod.get_function('RDPKernel')
        if rType == 2:
            imD = cp.asarray(im)
            f = torch.zeros(Nx * Ny * Nz, dtype=torch.float32, device='cuda')
            fD = cp.asarray(f)
        else:
            f = cp.zeros(Nx * Ny * Nz, dtype=cp.float32)
        chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
        array = cp.cuda.texture.CUDAarray(chl, Nx, Ny, Nz)
        if rType == 2:
            array.copy_from(imD.reshape((Nz, Ny, Nx)))
        else:
            array.copy_from(im.reshape((Nz, Ny, Nx)))
        res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
        tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeBorder, cp.cuda.runtime.cudaAddressModeBorder,cp.cuda.runtime.cudaAddressModeBorder), 
                                                filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=0)
        d_im = cp.cuda.texture.TextureObject(res, tdes)
        if rType == 2:
            kIndLoc = (fD,)
        else:
            kIndLoc = (f,)
        kIndLoc += (d_im,)
        # kIndLoc += (cp.array([cp.int32(Nx), cp.int32(Ny), cp.int32(Nz)], dtype=cp.int32),)
        # kIndLoc += (cp.array([cp.int32(Nx), cp.int32(Ny), cp.int32(Nz)], dtype=cp.int32),)
        kIndLoc += (cp.int32(Nx),)
        kIndLoc += (cp.int32(Ny),)
        kIndLoc += (cp.int32(Nz),)
        kIndLoc += (cp.int32(Nx),)
        kIndLoc += (cp.int32(Ny),)
        kIndLoc += (cp.int32(Nz),)
        kIndLoc += (cp.float32(gamma),)
        kIndLoc += (cp.float32(epps),)
        kIndLoc += (cp.float32(beta),)
        knlRDP((globalSize[0] // localSize[0], globalSize[1] // localSize[1], globalSize[2]), (localSize[0], localSize[1], 1), kIndLoc)
        if rType == 2:
            torch.cuda.synchronize()
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
        NLRD, 4 = NL Lange, 5 = NL filtered with Lange, and 6 = NLGGMRF.
        
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
        Lower values smooth the image, while larger make it sharper.
        
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
    if rType == 0:
        import pyopencl as cl
        import arrayfire as af
    elif rType == 1:
        import cupy as cp
    elif rType == 2:
        import cupy as cp
        import torch
    elif rType == 3:
        import pyopencl as cl
    import os
    import numpy as np
    
    if type(Nx) == np.ndarray:
        Nx = Nx.item()
    if type(Ny) == np.ndarray:
        Ny = Ny.item()
    if type(Nz) == np.ndarray:
        Nz = Nz.item()
        
    if rType == 0:
        ctx = af.opencl.get_context(retain=True)
        clctx = cl.Context.from_int_ptr(ctx)
        q = af.opencl.get_queue(True)
        queue = cl.CommandQueue.from_int_ptr(q)
        bOpt = ('-cl-single-precision-constant', '-DOPENCL', '-DCAST=float',)
    elif rType == 1 or rType == 2:
        bOpt = ('-DCUDA', '-DPYTHON',)
    elif rType == 3:
        bOpt = ('-cl-single-precision-constant', '-DOPENCL', '-DCAST=float',)
    
    bOpt += ('-DNLM_', '-DUSEIMAGES', '-DLOCAL_SIZE=16', '-DLOCAL_SIZE2=16', '-DNLTYPE=' + str(NLType), '-DSWINDOWX=' + str(SW[0]), 
             '-DSWINDOWY=' + str(SW[1]), '-DSWINDOWZ=' + str(SW[2]), '-DPWINDOWX=' + str(PW[0]), '-DPWINDOWY=' + str(PW[1]), 
             '-DPWINDOWZ=' + str(PW[2]),)
    if (useAdaptive):
        bOpt += ('-DNLMADAPTIVE',)
    if type(im) == type(refIm):
        bOpt += ('-DNLMREF',)
        useRef = True
    else:
        useRef = False
    
    epps = 1e-8
    localSize = (16, 16, 1)
    apu = [Nx % localSize[0], Ny % localSize[1], 0]
    erotus = [0] * 2
    if apu[0] > 0:
        erotus[0] = localSize[0] - apu[0]
    if apu[1] > 0:
        erotus[1] = localSize[1] - apu[1]
    
    globalSize = [None] * 3
    globalSize[0] = Nx + erotus[0]
    globalSize[1] = Ny + erotus[1]
    globalSize[2] = Nz
    
    headerDir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', '..', '..', 'opencl')) + "/"
    with open(headerDir + 'general_opencl_functions.h', encoding="utf8") as f:
        hlines = f.read()
    with open(headerDir + 'auxKernels.cl', encoding="utf8") as f:
        lines = f.read()
    lines = hlines + lines
    x = np.linspace(-PW[0], PW[0], 2 * PW[0] + 1, dtype=np.float32)
    y = np.linspace(-PW[1], PW[1], 2 * PW[1] + 1, dtype=np.float32)
    z = np.linspace(-PW[2], PW[2], 2 * PW[2] + 1, dtype=np.float32)
    gaussK = np.exp(-(np.add.outer(np.add.outer(x**2 / (2*STD**2), y**2 / (2*STD**2)), z**2 / (2*STD**2))))
    gaussK = gaussK.flatten('F').astype(dtype=np.float32)
    if NLType == 4 or NLType == 5:
        gamma = phi
    if rType == 0 or rType == 3:
        queue.finish()
        d_Nxyz = cl.cltypes.make_int3(Nx, Ny, Nz)
        prg = cl.Program(clctx, lines).build(bOpt)
        knlNLM = prg.NLM
        if rType == 0:
            imPtr = im.raw_ptr()
            imD = cl.MemoryObject.from_int_ptr(imPtr)
            f = af.data.constant(0, Nx * Ny * Nz, dtype=af.Dtype.f32)
            fPtr = f.raw_ptr()
            fD = cl.MemoryObject.from_int_ptr(fPtr)
            if useRef:
                RefimPtr = refIm.raw_ptr()
                RefimD = cl.MemoryObject.from_int_ptr(RefimPtr)
        else:
            f = cl.array.zeros(queue, Nx * Ny * Nz, dtype=cl.cltypes.float)
        d_gaussian = cl.array.to_device(queue, gaussK)
        imformat = cl.ImageFormat(cl.channel_order.A, cl.channel_type.FLOAT)
        mf = cl.mem_flags
        d_im = cl.Image(clctx, mf.READ_ONLY, imformat, shape=(Nx, Ny, Nz))
        if useRef:
            d_refIm = cl.Image(clctx, mf.READ_ONLY, imformat, shape=(Nx, Ny, Nz))
        if rType == 0:
            cl.enqueue_copy(queue, d_im, imD, offset=(0), origin=(0,0,0), region=(Nx, Ny, Nz))
            af.device.unlock_array(im)
            if useRef:
                cl.enqueue_copy(queue, d_refIm, RefimD, offset=(0), origin=(0,0,0), region=(Nx, Ny, Nz))
                af.device.unlock_array(refIm)
        else:
            cl.enqueue_copy(queue, d_im, im.data, offset=(0), origin=(0,0,0), region=(Nx, Ny, Nz))
            if useRef:
                cl.enqueue_copy(queue, d_refIm, refIm.data, offset=(0), origin=(0,0,0), region=(Nx, Ny, Nz))
        kIndLoc = 0
        if rType == 0:
            knlNLM.set_arg(kIndLoc, fD)
        else:
            knlNLM.set_arg(kIndLoc, f.data)
        kIndLoc += 1
        knlNLM.set_arg(kIndLoc, d_im)
        kIndLoc += 1
        knlNLM.set_arg(kIndLoc, d_gaussian.data)
        kIndLoc += 1
        knlNLM.set_arg(kIndLoc, d_Nxyz)
        kIndLoc += 1
        knlNLM.set_arg(kIndLoc, d_Nxyz)
        kIndLoc += 1
        knlNLM.set_arg(kIndLoc, (cl.cltypes.float)(h))
        kIndLoc += 1
        knlNLM.set_arg(kIndLoc, (cl.cltypes.float)(epps))
        kIndLoc += 1
        knlNLM.set_arg(kIndLoc, (cl.cltypes.float)(beta))
        if NLType >= 3:
            kIndLoc += 1
            knlNLM.set_arg(kIndLoc, (cl.cltypes.float)(gamma))
        if NLType == 6:
            kIndLoc += 1
            knlNLM.set_arg(kIndLoc, (cl.cltypes.float)(GGMRFpqc[0]))
            kIndLoc += 1
            knlNLM.set_arg(kIndLoc, (cl.cltypes.float)(GGMRFpqc[1]))
            kIndLoc += 1
            knlNLM.set_arg(kIndLoc, (cl.cltypes.float)(GGMRFpqc[2]))
        if useAdaptive:
            kIndLoc += 1
            knlNLM.set_arg(kIndLoc, (cl.cltypes.float)(adaptiveConstant))
        if useRef:
            kIndLoc += 1
            knlNLM.set_arg(kIndLoc, d_refIm)
        cl.enqueue_nd_range_kernel(queue, knlNLM, globalSize, localSize)
        queue.finish()
        if rType == 0:
            af.device.unlock_array(f)
    else:
        mod = cp.RawModule(code=lines, options=bOpt)
        knlNLM = mod.get_function('NLM')
        if rType == 2:
            imD = cp.asarray(im)
            f = torch.zeros(Nx * Ny * Nz, dtype=torch.float32, device='cuda')
            fD = cp.asarray(f)
            if useRef:
                RefimD = cp.asarray(refIm)
        else:
            f = cp.zeros(Nx * Ny * Nz, dtype=cp.float32)
        d_gaussian = cp.asarray(gaussK)
        chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
        array = cp.cuda.texture.CUDAarray(chl, Nx, Ny, Nz)
        if useRef:
            arrayRef = cp.cuda.texture.CUDAarray(chl, Nx, Ny, Nz)
        if rType == 2:
            array.copy_from(imD.reshape((Nz, Ny, Nx)))
            if useRef:
                arrayRef.copy_from(RefimD.reshape((Nz, Ny, Nx)))
        else:
            array.copy_from(im.reshape((Nz, Ny, Nx)))
            if useRef:
                arrayRef.copy_from(refIm.reshape((Nz, Ny, Nx)))
        res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
        tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeBorder, cp.cuda.runtime.cudaAddressModeBorder,cp.cuda.runtime.cudaAddressModeBorder), 
                                                filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=0)
        d_im = cp.cuda.texture.TextureObject(res, tdes)
        if useRef:
            resRef = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=arrayRef)
            d_refIm = cp.cuda.texture.TextureObject(resRef, tdes)
        if rType == 2:
            kIndLoc = (fD,)
        else:
            kIndLoc = (f,)
        kIndLoc += (d_im,)
        kIndLoc += (d_gaussian,)
        # kIndLoc += (cp.array([cp.int32(Nx), cp.int32(Ny), cp.int32(Nz)], dtype=cp.int32),)
        # kIndLoc += (cp.array([cp.int32(Nx), cp.int32(Ny), cp.int32(Nz)], dtype=cp.int32),)
        kIndLoc += (cp.int32(Nx),)
        kIndLoc += (cp.int32(Ny),)
        kIndLoc += (cp.int32(Nz),)
        kIndLoc += (cp.int32(Nx),)
        kIndLoc += (cp.int32(Ny),)
        kIndLoc += (cp.int32(Nz),)
        kIndLoc += (cp.float32(h),)
        kIndLoc += (cp.float32(epps),)
        kIndLoc += (cp.float32(beta),)
        if NLType >= 3:
            kIndLoc += (cp.float32(gamma),)
        if NLType == 6:
            kIndLoc += (cp.float32(GGMRFpqc[0]),)
            kIndLoc += (cp.float32(GGMRFpqc[1]),)
            kIndLoc += (cp.float32(GGMRFpqc[2]),)
        if useAdaptive:
            kIndLoc += (cp.float32(adaptiveConstant),)
        if useRef:
            kIndLoc += (d_refIm,)
        knlNLM((globalSize[0] // localSize[0], globalSize[1] // localSize[1], globalSize[2]), (localSize[0], localSize[1], 1), kIndLoc)
        if rType == 2:
            torch.cuda.synchronize()
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
    if rType == 0:
        import pyopencl as cl
        import arrayfire as af
    elif rType == 1:
        import cupy as cp
    elif rType == 2:
        import cupy as cp
        import torch
    elif rType == 3:
        import pyopencl as cl
    import os
    import numpy as np
    
    if type(Nx) == np.ndarray:
        Nx = Nx.item()
    if type(Ny) == np.ndarray:
        Ny = Ny.item()
    if type(Nz) == np.ndarray:
        Nz = Nz.item()
        
    if rType == 0:
        ctx = af.opencl.get_context(retain=True)
        clctx = cl.Context.from_int_ptr(ctx)
        q = af.opencl.get_queue(True)
        queue = cl.CommandQueue.from_int_ptr(q)
        bOpt = ('-cl-single-precision-constant', '-DOPENCL', '-DCAST=float',)
    elif rType == 1 or rType == 2:
        bOpt = ('-DCUDA', '-DPYTHON',)
    elif rType == 3:
        bOpt = ('-cl-single-precision-constant', '-DOPENCL', '-DCAST=float',)
    
    bOpt += ('-DTVGRAD', '-DUSEIMAGES', '-DLOCAL_SIZE=16', '-DLOCAL_SIZE2=16',)
    
    if Lange:
        bOpt += ('-DSATV',)
    
    localSize = (16, 16, 1)
    apu = [Nx % localSize[0], Ny % localSize[1], 0]
    erotus = [0] * 2
    if apu[0] > 0:
        erotus[0] = localSize[0] - apu[0]
    if apu[1] > 0:
        erotus[1] = localSize[1] - apu[1]
    
    globalSize = [None] * 3
    globalSize[0] = Nx + erotus[0]
    globalSize[1] = Ny + erotus[1]
    globalSize[2] = Nz
    
    headerDir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', '..', '..', 'opencl')) + "/"
    with open(headerDir + 'general_opencl_functions.h', encoding="utf8") as f:
        hlines = f.read()
    with open(headerDir + 'auxKernels.cl', encoding="utf8") as f:
        lines = f.read()
    lines = hlines + lines
    if rType == 0 or rType == 3:
        queue.finish()
        d_Nxyz = cl.cltypes.make_int3(Nx, Ny, Nz)
        prg = cl.Program(clctx, lines).build(bOpt)
        knlTV = prg.TVKernel
        if rType == 0:
            imPtr = im.raw_ptr()
            imD = cl.MemoryObject.from_int_ptr(imPtr)
            f = af.data.constant(0, Nx * Ny * Nz, dtype=af.Dtype.f32)
            fPtr = f.raw_ptr()
            fD = cl.MemoryObject.from_int_ptr(fPtr)
        else:
            f = cl.array.zeros(queue, Nx * Ny * Nz, dtype=cl.cltypes.float)
        imformat = cl.ImageFormat(cl.channel_order.A, cl.channel_type.FLOAT)
        mf = cl.mem_flags
        d_im = cl.Image(clctx, mf.READ_ONLY, imformat, shape=(Nx, Ny, Nz))
        if rType == 0:
            cl.enqueue_copy(queue, d_im, imD, offset=(0), origin=(0,0,0), region=(Nx, Ny, Nz))
            af.device.unlock_array(im)
        else:
            cl.enqueue_copy(queue, d_im, im.data, offset=(0), origin=(0,0,0), region=(Nx, Ny, Nz))
        kIndLoc = 0
        if rType == 0:
            knlTV.set_arg(kIndLoc, fD)
        else:
            knlTV.set_arg(kIndLoc, f.data)
        kIndLoc += 1
        knlTV.set_arg(kIndLoc, d_im)
        kIndLoc += 1
        knlTV.set_arg(kIndLoc, d_Nxyz)
        kIndLoc += 1
        knlTV.set_arg(kIndLoc, d_Nxyz)
        kIndLoc += 1
        knlTV.set_arg(kIndLoc, (cl.cltypes.float)(sigma))
        kIndLoc += 1
        knlTV.set_arg(kIndLoc, (cl.cltypes.float)(sValue))
        kIndLoc += 1
        knlTV.set_arg(kIndLoc, (cl.cltypes.float)(beta))
        cl.enqueue_nd_range_kernel(queue, knlTV, globalSize, localSize)
        queue.finish()
        if rType == 0:
            af.device.unlock_array(f)
    else:
        mod = cp.RawModule(code=lines, options=bOpt)
        knlTV = mod.get_function('TVKernel')
        if rType == 2:
            imD = cp.asarray(im)
            f = torch.zeros(Nx * Ny * Nz, dtype=torch.float32, device='cuda')
            fD = cp.asarray(f)
        else:
            f = cp.zeros(Nx * Ny * Nz, dtype=cp.float32)
        chl = cp.cuda.texture.ChannelFormatDescriptor(32,0,0,0, cp.cuda.runtime.cudaChannelFormatKindFloat)
        array = cp.cuda.texture.CUDAarray(chl, Nx, Ny, Nz)
        if rType == 2:
            array.copy_from(imD.reshape((Nz, Ny, Nx)))
        else:
            array.copy_from(im.reshape((Nz, Ny, Nx)))
        res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=array)
        tdes= cp.cuda.texture.TextureDescriptor(addressModes=(cp.cuda.runtime.cudaAddressModeBorder, cp.cuda.runtime.cudaAddressModeBorder,cp.cuda.runtime.cudaAddressModeBorder), 
                                                filterMode=cp.cuda.runtime.cudaFilterModePoint, normalizedCoords=0)
        d_im = cp.cuda.texture.TextureObject(res, tdes)
        if rType == 2:
            kIndLoc = (fD,)
        else:
            kIndLoc = (f,)
        kIndLoc += (d_im,)
        # kIndLoc += (cp.array([cp.int32(Nx), cp.int32(Ny), cp.int32(Nz)], dtype=cp.int32),)
        # kIndLoc += (cp.array([cp.int32(Nx), cp.int32(Ny), cp.int32(Nz)], dtype=cp.int32),)
        kIndLoc += (cp.int32(Nx),)
        kIndLoc += (cp.int32(Ny),)
        kIndLoc += (cp.int32(Nz),)
        kIndLoc += (cp.int32(Nx),)
        kIndLoc += (cp.int32(Ny),)
        kIndLoc += (cp.int32(Nz),)
        kIndLoc += (cp.float32(sigma),)
        kIndLoc += (cp.float32(sValue),)
        kIndLoc += (cp.float32(beta),)
        knlTV((globalSize[0] // localSize[0], globalSize[1] // localSize[1], globalSize[2]), (localSize[0], localSize[1], 1), kIndLoc)
        if rType == 2:
            torch.cuda.synchronize()
    return f