# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 17:44:17 2024

@author: Ville-Veikko Wettenhovi
"""


class _PowerMethodOps:
    """Tiny array-namespace adapter selected once (per powerMethod call) from
    A.useAF/useTorch/useCuPy/else (plain PyOpenCL), so the backend branching
    that used to be repeated at every random-init/dot/norm/divide site in
    powerMethod only happens here.

    randn_abs(n): abs() of a standard-normal random vector of length n, using
    this instance's own (optionally seeded) random state, matching what each
    backend previously did at every draw site.

    dot(a, b): the backend's dot product, converted to a plain Python float
    (the immediate host conversion needed for the L eigenvalue estimate;
    unlike the original per-backend inline code, the numerator/denominator
    division itself then happens in Python rather than on-device, so results
    can differ from the original by float32-vs-float64 rounding, per the
    refactor's brief).

    norm(a): the backend's own norm of a, returned AS THE BACKEND NORMALLY
    RETURNS IT (a device scalar for torch/CuPy/PyOpenCL, a Python float for
    ArrayFire) -- this preserves the original on-device division used for
    the actual vector normalization step.

    scale(a, s): a / s, backend-agnostic (every array type here supports
    dividing by a Python float or a device 0-d scalar).
    """

    def __init__(self, A, seed=None):
        self.A = A
        self.af = None
        self.torch = None
        self.cp = None
        self.cl = None
        self.clmath = None
        if A.useAF:
            import arrayfire as af
            self.backend = 'af'
            self.af = af
            if seed is not None:
                af.set_seed(seed)
        elif A.useTorch:
            import torch
            self.backend = 'torch'
            self.torch = torch
            self._device = 'mps' if getattr(A, 'useMetal', False) else 'cuda'
            self._gen = None
            if seed is not None:
                self._gen = torch.Generator(device=self._device)
                self._gen.manual_seed(seed)
        elif A.useCUDA:
            import cupy as cp
            self.backend = 'cupy'
            self.cp = cp
            self._rng = cp.random.default_rng(seed)
        else:
            import numpy as np
            import pyopencl as cl
            from pyopencl import clmath
            self.backend = 'cl'
            self._np = np
            self.cl = cl
            self.clmath = clmath
            self._rng = np.random.default_rng(seed)

    def randn_abs(self, n):
        if self.backend == 'af':
            return self.af.abs(self.af.randn(n))
        elif self.backend == 'torch':
            torch = self.torch
            if self._gen is not None:
                return torch.randn(n, dtype=torch.float32, device=self._device, generator=self._gen).abs()
            return torch.randn(n, dtype=torch.float32, device=self._device).abs()
        elif self.backend == 'cupy':
            return self.cp.abs(self._rng.standard_normal(n, dtype=self.cp.float32))
        else:
            return self.cl.array.to_device(self.A.queue, self._np.abs(self._rng.standard_normal(n, dtype=self._np.float32)))

    def dot(self, a, b):
        if self.backend == 'af':
            return (self.af.dot(a, b).to_ndarray()).item()
        elif self.backend == 'torch':
            return self.torch.dot(a, b).cpu().numpy().item()
        elif self.backend == 'cupy':
            return (self.cp.dot(a, b).get()).item()
        else:
            return (self.cl.array.dot(a, b).get(self.A.queue)).item()

    def norm(self, a):
        if self.backend == 'af':
            return self.af.norm(a)
        elif self.backend == 'torch':
            return self.torch.norm(a)
        elif self.backend == 'cupy':
            return self.cp.sqrt(self.cp.dot(a, a))
        else:
            return self.clmath.sqrt(self.cl.array.dot(a, a))

    def scale(self, a, s):
        return a / s


def powerMethod(A, seed=None):
    """
    Power method for estimating the largest eigenvalue of A^T * A (used to
    derive a step-size normalization constant for several algorithms).

    Args:
        A: The (initialized or not-yet-initialized) projectorClass instance.

        seed: Optional RNG seed for the random starting vector(s). Default
        None reproduces the previous, unseeded behaviour. Passing a seed
        makes the method's randomness (and hence its result) reproducible,
        which is useful for testing.

    Returns:
        L: The estimated largest eigenvalue (or a list of one per volume, for
        multi-resolution/multi-volume reconstructions).
    """
    if not(A.projectorInitialized):
        A.initProj()
    from omegatomo.util.measprecond import applyMeasPreconditioning
    ops = _PowerMethodOps(A, seed=seed)

    L = [None] * (A.nMultiVolumes + 1)
    if A.nMultiVolumes > 0:
        x = [None] * (A.nMultiVolumes + 1)
    for i in range(A.nMultiVolumes + 1):
        if A.nMultiVolumes > 0:
            x[i] = ops.randn_abs(A.N[i].item())
            x[i] = ops.scale(x[i], ops.norm(x[i]))
        else:
            x = ops.randn_abs(A.N[0].item())
            x = ops.scale(x, ops.norm(x))
    if A.nMultiVolumes > 0:
        i = 0
        for k in range(A.powerIterations):
            if A.useTorch and getattr(A, 'useMetal', False):
                x2 = A * [x[0], *[ops.torch.zeros_like(v) for v in x[1:]]]
            else:
                x2 = A * x[0]
            x2 = applyMeasPreconditioning(A, x2, subIter=0)
            x2 = A.T() * x2
            L[i] = ops.dot(x[i], x2[i]) / ops.dot(x[i], x[i]) * A.subsets
            x[i] = ops.scale(x2[i], ops.norm(x2[i]))
            if A.verbose > 0:
                print('Largest eigenvalue at iteration ' + str(k) + ' in the main volume is ' + str(L[i]))
    for k in range(A.powerIterations):
        if A.nMultiVolumes > 0:
            x2 = A * x
            x2 = applyMeasPreconditioning(A, x2, subIter=0)
            x2 = A.T() * x2
            for i in range(A.nMultiVolumes + 1):
                if i > 0:
                    L[i] = ops.dot(x[i], x2[i]) / ops.dot(x[i], x[i]) * A.subsets
                x[i] = ops.scale(x2[i], ops.norm(x2[i]))
                if A.verbose > 0 and i > 0:
                    print('Largest eigenvalue at iteration ' + str(k) + ' and in volume ' + str(i) + ' is ' + str(L[i]))
        else:
            x2 = A * x
            x2 = applyMeasPreconditioning(A, x2, subIter=0)
            x2 = A.T() * x2
            L = ops.dot(x, x2) / ops.dot(x, x) * A.subsets
            x = ops.scale(x2, ops.norm(x2))
            if A.verbose > 0:
                print('Largest eigenvalue at iteration ' + str(k) + ' is ' + str(L))
    if A.nMultiVolumes == 0:
        L = 1. / L
    else:
        for i in range(A.nMultiVolumes + 1):
            L[i] = 1. / L[i]
    if A.useAF:
        ops.af.device_gc()
    return L
