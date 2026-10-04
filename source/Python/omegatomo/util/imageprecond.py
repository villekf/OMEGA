# -*- coding: utf-8 -*-
"""
Image-domain preconditioning for the custom-algorithm / powerMethod path.

Mirrors MATLAB's applyImagePreconditioning.m. See that file (and
options.precondTypeImage) for the authoritative definition of each of the
seven preconditioner types (index i here == MATLAB's precondTypeImage(i+1)):

    0 : Diagonal normalization preconditioner,  input = input / D
    1 : EM preconditioner,                      input = input * (im / D)
    2 : IEM preconditioner,                     input = input * (max(im, max(VAL, ref)) / D)
    3 : Momentum-like preconditioner,           input = input * alphaPrecond(kk)
    4 : Gradient-based preconditioner           NOT IMPLEMENTED (see below)
    5 : Filtering-based preconditioner,         input = filtering2D(filterIm, input, Nf)
    6 : Curvature-based preconditioner          NOT IMPLEMENTED (see below)

Types 4 and 6 depend on inputs MATLAB computes internally
(applyImagePreconditioning.m -> gradientPreconditioner.m for options.gradF,
and reconstructions_main.m's SPS-style E/Sino curvature backprojection for
options.dP) that have no Python equivalent yet, and are not invented here;
enabling either type raises NotImplementedError naming exactly what is
missing (see the two functions below for the proposed scope of a real
implementation).

Like applyMeasPreconditioning.m/applyMeasPreconditioning (measprecond.py),
this intentionally omits the MATLAB function's `options.verbose >= 3` disp()
diagnostics -- this runs inside the power-iteration hot loop.

Copyright (C) 2024-2025 Ville-Veikko Wettenhovi

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
"""

VAL = 0.00001


def _maximum(options, a, b):
    """Backend-dispatched elementwise max(a, b). `b` may be a plain Python
    scalar or an array of the same backend as `a`."""
    if options.useAF:
        import arrayfire as af
        if not hasattr(b, 'dims'):
            b = float(b)
            return af.maxof(a, af.constant(b, *a.dims(), dtype=a.dtype()))
        return af.maxof(a, b)
    elif options.useTorch:
        import torch
        if isinstance(b, (int, float)):
            return torch.clamp(a, min=b)
        return torch.maximum(a, b)
    elif options.useCuPy:
        import cupy as cp
        return cp.maximum(a, b)
    else:
        import numpy as np
        import pyopencl.array as clarray
        if isinstance(b, (int, float)):
            # pyopencl.array.maximum(array, python_scalar) generates an
            # elementwise kernel whose fmax_nanprop(a, b) macro call is
            # ambiguous to overload resolution on at least some OpenCL C
            # compilers (observed: Intel's CPU OpenCL runtime) when `b` is a
            # bare scalar rather than a same-shaped array -- broadcast the
            # scalar to a device array first to sidestep that kernel
            # entirely (plain elementwise +/*, unambiguous).
            b = a * np.float32(0) + np.float32(b)
        return clarray.maximum(a, b)


def _filter2D_af(options, var, Nx, Ny, Nz, Nf):
    import arrayfire as af
    if not hasattr(options, 'filterImG'):
        options.filterImG = af.interop.np_to_af_array(options.filterIm.astype('float32', copy=False))
    filterG = options.filterImG
    im = af.moddims(var, Nx, d1=Ny, d2=Nz)
    temp = af.fft2(im, Nf, Nf)
    tiled = af.tile(filterG, 1, 1, d2=temp.shape[2])
    temp = temp * tiled
    af.ifft2_inplace(temp)
    out = af.real(temp[0:Nx, 0:Ny, :])
    return af.flat(out)


def _filter2D_torch(options, var, Nx, Ny, Nz, Nf):
    import torch
    if not hasattr(options, 'filterImG'):
        options.filterImG = torch.as_tensor(options.filterIm, dtype=var.dtype if var.is_floating_point() else torch.float32,
                                             device=var.device)
    filterT = options.filterImG
    # var is F-order-flattened (Nx fastest); reshape(Nz, Ny, Nx) (C-order) + permute
    # gives a tensor indexed [x, y, z] equal to the Fortran-order volume (see
    # the same trick used for CuPy/torch texture uploads in util/priors.py).
    im = var.reshape(Nz, Ny, Nx).permute(2, 1, 0)
    temp = torch.fft.fft2(im, s=(Nf, Nf), dim=(0, 1))
    temp = temp * filterT.unsqueeze(-1)
    temp = torch.fft.ifft2(temp, dim=(0, 1))
    out = temp[:Nx, :Ny, :].real
    return out.permute(2, 1, 0).reshape(-1)


def _filter2D_cupy(options, var, Nx, Ny, Nz, Nf):
    import cupy as cp
    filterC = cp.asarray(options.filterIm)
    im = cp.reshape(var, (Nx, Ny, Nz), order='F')
    temp = cp.fft.fft2(im, s=(Nf, Nf), axes=(0, 1))
    temp = temp * filterC[:, :, None]
    temp = cp.fft.ifft2(temp, axes=(0, 1))
    out = cp.real(temp[:Nx, :Ny, :])
    return cp.ravel(out, order='F')


def _filter2D_cl(options, var, Nx, Ny, Nz, Nf):
    import numpy as np
    import pyopencl as cl
    queue = var.queue
    varNp = var.get()
    im = np.reshape(varNp, (Nx, Ny, Nz), order='F')
    temp = np.fft.fft2(im, s=(Nf, Nf), axes=(0, 1))
    temp = temp * options.filterIm[:, :, None]
    temp = np.fft.ifft2(temp, axes=(0, 1))
    out = np.real(temp[:Nx, :Ny, :])
    outFlat = np.ravel(out, order='F').astype(np.float32)
    return cl.array.to_device(queue, outFlat)


def applyImagePreconditioning(options, input, im, kk, ii=0):
    """
    Computes the image-based preconditioning for the input data.

    Mirrors applyImagePreconditioning.m. Called from powerMethod at the same
    point MATLAB's powerMethod.m calls it: right after the backprojection
    step, using the pre-update image estimate as `im`.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.
    input : arrayfire array, torch tensor, cupy array or pyopencl array
        The (just backprojected) image-domain data to be preconditioned.
    im : same type as `input`
        The current image estimate (x before this power-iteration update),
        matching MATLAB's `im` argument.
    kk : int
        One-based power-iteration index, matching applyImagePreconditioning.m's
        `kk` (i.e. the caller should pass the loop counter + 1, not a
        zero-based index).
    ii : int, optional
        Zero-based volume index (0 = main volume). Used to index options.D
        when it is given as a list/tuple (one per volume), mirroring
        MATLAB's `iscell(options.D)` handling, and to look up
        options.Nx[ii]/Ny[ii]/Nz[ii] for the type-5 filtering reshape.
        ii + 1 is also used to build the 'referenceImageN' attribute name
        for multi-volume IEM, mirroring applyImagePreconditioning.m's
        dynamic fieldname ('referenceImage', 'referenceImage2', ...).
        options.alphaPrecond and options.filterIm are used as single
        shared arrays regardless of ii, matching how MATLAB's
        prepass_phase.m computes them (not per-volume).

    Returns
    -------
    input : same type as the input `input`
        The preconditioned data.

    """
    if options.precondTypeImage[4].item() and kk >= int(options.gradInitIter) + 1:
        raise NotImplementedError(
            "Image preconditioner type 4 (gradient-based preconditioner, "
            "options.precondTypeImage[4]) is not implemented in the Python "
            "powerMethod/custom-algorithm path. MATLAB computes this via "
            "gradientPreconditioner.m, which needs: (1) a finite-difference "
            "image gradient (computeGradient.m, controlled by options.derivType) "
            "with no Python equivalent, and (2) options.gradF, cached per "
            "volume and refreshed once between options.gradInitIter+1 and "
            "options.gradLastIter+1, clamped to [options.gradV1, options.gradV2]. "
            "A Python implementation would need a computeGradient-equivalent "
            "(finite differences matching options.derivType) plus the "
            "gradF caching/clamping logic above before this can be enabled."
        )

    if options.precondTypeImage[3].item():
        input = input * options.alphaPrecond[kk - 1]

    if options.precondTypeImage[0].item() or options.precondTypeImage[1].item() or options.precondTypeImage[2].item():
        if not hasattr(options, 'D') or options.D is None:
            raise ValueError(
                "Image preconditioner types 0/1/2 (diagonal/EM/IEM, "
                "options.precondTypeImage[0:3]) require options.D (the image "
                "sensitivity image), which is not set. Compute it the same way "
                "MATLAB's reconstructions_main.m does for the built-in algorithms: "
                "D = 1 + sum over subsets of A.T() applied to a ones-vector of "
                "that subset's measurement length (optionally PSF-convolved), "
                "divided by the number of subsets, with any exact-zero entries "
                "replaced by 1. For multi-volume reconstructions, options.D may "
                "be a list/tuple with one such image per volume, mirroring "
                "MATLAB's iscell(options.D)."
            )
        D = options.D[ii] if isinstance(options.D, (list, tuple)) else options.D
        if options.precondTypeImage[0].item():
            input = input / D
        elif options.precondTypeImage[1].item():
            input = input * (im / D)
        elif options.precondTypeImage[2].item():
            fieldname = 'referenceImage' if ii == 0 else 'referenceImage' + str(ii + 1)
            ref = getattr(options, fieldname, None)
            if ref is None or isinstance(ref, str):
                raise ValueError(
                    f"IEM preconditioner (type 2, options.precondTypeImage[2]) requires "
                    f"options.{fieldname} to be set to a reference image array matching "
                    f"the shape/backend of the current image estimate."
                )
            input = input * (_maximum(options, im, _maximum(options, ref, VAL)) / D)

    if options.precondTypeImage[6].item():
        raise NotImplementedError(
            "Image preconditioner type 6 (curvature-based preconditioner, "
            "options.precondTypeImage[6]) is not implemented in the Python "
            "powerMethod/custom-algorithm path. MATLAB computes options.dP in "
            "reconstructions_main.m as an SPS-style curvature term: "
            "E = 1 + sum over subsets of forwardProject(ones), then an "
            "SPS/MBSREM-style curvature formula combining E with the "
            "measurement data (Sino), the flat-field (CT) or randoms/scatter "
            "correction (PET/SPECT) terms, then options.dP = subsets / "
            "backwardProject(E) with NaN/Inf entries replaced by 1. A Python "
            "implementation would need that same forward/backward pass plus "
            "the CT- and emission-specific curvature formulas from "
            "reconstructions_main.m (lines ~288-356) before this can be enabled."
        )

    if options.precondTypeImage[5].item() and kk <= int(options.filteringIterations):
        Nx = int(options.Nx[ii])
        Ny = int(options.Ny[ii])
        Nz = int(options.Nz[ii])
        Nf = int(options.Nf)
        if options.useAF:
            input = _filter2D_af(options, input, Nx, Ny, Nz, Nf)
        elif options.useTorch:
            input = _filter2D_torch(options, input, Nx, Ny, Nz, Nf)
        elif options.useCuPy:
            input = _filter2D_cupy(options, input, Nx, Ny, Nz, Nf)
        else:
            input = _filter2D_cl(options, input, Nx, Ny, Nz, Nf)

    return input
