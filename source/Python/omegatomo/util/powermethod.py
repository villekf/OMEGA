# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 17:44:17 2024

@author: Ville-Veikko Wettenhovi
"""

from omegatomo.util.measprecond import applyMeasPreconditioning

def powerMethod(A):
    if A.useAF:
        import arrayfire as af
    elif A.useTorch:
        import torch
    elif A.useCUDA:
        if A.useCuPy:
            import cupy as cp
        else:
            import pycuda as cuda
            from pycuda import cumath
    else:
        import pyopencl as cl
        from pyopencl import clmath
    L = [None] * (A.nMultiVolumes + 1)
    if A.nMultiVolumes > 0:
        x = [None] * (A.nMultiVolumes + 1)
    for i in range(A.nMultiVolumes + 1):
        if A.nMultiVolumes > 0:
            if A.useAF:
                x[i] = af.abs(af.randn(A.N[i].item()))
                x[i] = x[i] / af.norm(x[i])
            else:
                if A.useTorch:
                    x[i] = torch.randn(A.N[i].item(), dtype=torch.float32, device='cuda').abs()
                    x[i] = x[i] / torch.norm(x[i])
                else:
                    if A.useCuPy:
                        rng = cp.random.default_rng()
                        x[i] = cp.abs(rng.standard_normal(A.N[i].item(), dtype=cp.float32))
                        x[i] = x[i] / cp.sqrt(cp.dot(x[i], x[i]))
                    else:
                        import numpy as np
                        rng = np.random.default_rng()
                        if A.useCUDA:
                            x[i] = cuda.gpuarray.to_gpu(np.abs(rng.standard_normal(A.N[i].item(), dtype=np.float32)))
                            x[i] = x[i] / cumath.sqrt(cuda.gpuarray.dot(x[i], x[i]))
                        else:
                            x[i] = cl.array.to_device(A.queue, np.abs(rng.standard_normal(A.N[i].item(), dtype=np.float32)))
                            x[i] = x[i] / clmath.sqrt(cl.array.dot(x[i], x[i]))
        else:
            if A.useAF:
                x = af.abs(af.randn(A.N[0].item()))
                x = x / af.norm(x)
            else:
                if A.useTorch:
                    x = torch.randn(A.N[0].item(), dtype=torch.float32, device='cuda').abs()
                    x = x / torch.norm(x)
                else:
                    if A.useCuPy:
                        rng = cp.random.default_rng()
                        x = cp.abs(rng.standard_normal(A.N[0].item(), dtype=cp.float32))
                        x = x / cp.sqrt(cp.dot(x, x))
                    else:
                        import numpy as np
                        rng = np.random.default_rng()
                        if A.useCUDA:
                            x = cuda.gpuarray.to_gpu(np.abs(rng.standard_normal(A.N[0].item(), dtype=np.float32)))
                            x = x / cumath.sqrt(cuda.gpuarray.dot(x, x))
                        else:
                            x = cl.array.to_device(A.queue, np.abs(rng.standard_normal(A.N[0].item(), dtype=np.float32)))
                            x = x / clmath.sqrt(cl.array.dot(x, x))
    for k in range(A.powerIterations):
        x2 = A * x
        if A.useAF or A.useTorch:
            x2 = applyMeasPreconditioning(A, x2)
        x2 = A.T() * x2
        if A.nMultiVolumes > 0:
            for i in range(A.nMultiVolumes + 1):
                if A.useAF:
                    L[i] = (((af.dot(x[i], x2[i])) / (af.dot(x[i], x[i])) * A.subsets).to_ndarray()).item()
                    x[i] = x2[i] / af.norm(x2[i])
                else:
                    if A.useTorch:
                        LL = torch.dot(x[i], x2[i]) / torch.dot(x[i], x[i]) * A.subsets
                        L[i] = LL.cpu().numpy().item()
                        x[i] = x2[i] / torch.norm(x2[i])
                    else:
                        if A.useCUDA:
                            if A.useCuPy:
                                L[i] = (((cp.dot(x[i], x2[i])) / (cp.dot(x[i], x[i])) * A.subsets).get()).item()
                                x[i] = x2[i] / cp.sqrt(cp.dot(x2[i], x2[i]))
                            else:
                                L[i] = (((cuda.gpuarray.dot(x[i], x2[i])) / (cuda.gpuarray.dot(x[i], x[i])) * A.subsets).get()).item()
                                x[i] = x2[i] / cumath.sqrt(cuda.gpuarray.dot(x2[i], x2[i]))
                        else:
                            L[i] = (((cl.array.dot(x[i], x2[i])) / (cl.array.dot(x[i], x[i])) * A.subsets).get(A.queue)).item()
                            x[i] = x2[i] / clmath.sqrt(cl.array.dot(x2[i], x2[i]))
                if A.verbose > 0:
                    print('Largest eigenvalue at iteration ' + str(k) + ' and in volume ' + str(i) + ' is ' + str(L[i]))
        else:
            if A.useAF:
                L = (((af.dot(x, x2)) / (af.dot(x, x)) * A.subsets).to_ndarray()).item()
                x = x2 / af.norm(x2)
            else:
                if A.useTorch:
                    LL = torch.dot(x, x2) / torch.dot(x, x) * A.subsets
                    L = LL.cpu().numpy().item()
                    x = x2 / torch.norm(x2)
                else:
                    if A.useCUDA:
                        if A.useCuPy:
                            L = (((cp.dot(x, x2)) / (cp.dot(x, x)) * A.subsets).get()).item()
                            x = x2 / cp.sqrt(cp.dot(x2, x2))
                        else:
                            L = (((cuda.gpuarray.dot(x, x2)) / (cuda.gpuarray.dot(x, x)) * A.subsets).get()).item()
                            x = x2 / cuda.cumath.sqrt(cuda.gpuarray.dot(x2, x2))
                    else:
                        L = (((cl.array.dot(x, x2)) / (cl.array.dot(x, x)) * A.subsets).get(A.queue)).item()
                        x = x2 / clmath.sqrt(cl.array.dot(x2, x2))
            if A.verbose > 0:
                print('Largest eigenvalue at iteration ' + str(k) + ' is ' + str(L))
    if A.nMultiVolumes == 0:
        L = 1. / L
    else:
        for i in range(A.nMultiVolumes + 1):
            L[i] = 1. / L[i]
    if A.useAF:
        af.device_gc()
    return L