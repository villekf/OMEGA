# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 17:45:47 2024

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

def _apply_filter_af(options, var, filterAttr, filterSource, divide):
    """Shared reshape / FFT(axis 0, n=Nf) / filter-multiply-or-divide /
    IFFT / crop / flatten used by both applyMeasPreconditioning and
    circulantInverse, ArrayFire backend. `filterAttr` is the options
    attribute the (cached) device filter is stored under ('filterG' or
    'FilterG'), `filterSource` the NumPy attribute it is built from ('filter0'
    or 'Ffilter'), and `divide` selects multiply (False) vs divide (True)."""
    import arrayfire as af
    if not hasattr(options, filterAttr):
        setattr(options, filterAttr, af.interop.np_to_af_array(getattr(options, filterSource)))
    filterG = getattr(options, filterAttr)
    if (options.subsets > 1 and (options.subsetType == 5 or options.subsetType == 4)):
        if options.subsetType == 4:
            var = af.moddims(var, options.nRowsD, d1=var.elements() // options.nRowsD)
        else:
            var = af.moddims(var, options.nColsD, d1=var.elements() // options.nColsD)
    else:
        var = af.moddims(var, options.nRowsD, d1=options.nColsD, d2=var.elements() // (options.nRowsD * options.nColsD))
    temp = af.fft(var, options.Nf)
    tiled = af.tile(filterG, 1, d1=temp.shape[1], d2=temp.shape[2])
    if divide:
        temp /= tiled
    else:
        temp = temp * tiled
    af.eval(temp)
    af.ifft_inplace(temp)
    return af.flat(af.real(temp[:var.shape[0], :, :]))


def _apply_filter_torch(options, var, filterAttr, filterSource, divide):
    """Same as _apply_filter_af, torch backend. NOTE: keeps the torch path's
    existing reshape semantics exactly (reshape(-1, width)); this is not
    "fixed" here, only preserved."""
    import torch
    if not hasattr(options, filterAttr):
        setattr(options, filterAttr, torch.as_tensor(getattr(options, filterSource), dtype=var.dtype, device=var.device).reshape(-1))
    filterT = getattr(options, filterAttr)
    width = options.nColsD if options.subsets > 1 and options.subsetType == 5 else options.nRowsD
    var = var.reshape(-1, width)
    temp = torch.fft.fft(var, n=options.Nf, dim=-1)
    if divide:
        temp /= filterT
    else:
        temp *= filterT
    temp = torch.fft.ifft(temp, dim=-1)
    return temp.real[:, :width].reshape(-1).contiguous()


def _apply_filter_cupy(options, var, filterAttr, filterSource, divide):
    """Same as _apply_filter_af, CuPy backend."""
    import cupy as cp
    if not hasattr(options, filterAttr):
        setattr(options, filterAttr, cp.asarray(getattr(options, filterSource)))
    filterC = getattr(options, filterAttr)
    if options.subsets > 1 and (options.subsetType == 5 or options.subsetType == 4):
        if options.subsetType == 4:
            var = cp.reshape(var, (options.nRowsD, var.size // options.nRowsD), order='F')
        else:
            var = cp.reshape(var, (options.nColsD, var.size // options.nColsD), order='F')
    else:
        var = cp.reshape(var, (options.nRowsD, options.nColsD, var.size // (options.nRowsD * options.nColsD)), order='F')
    temp = cp.fft.fft(var, n=options.Nf, axis=0)
    filterShape = (-1,) + (1,) * (var.ndim - 1)
    if divide:
        temp /= filterC.reshape(filterShape, order='F')
    else:
        temp *= filterC.reshape(filterShape, order='F')
    temp = cp.fft.ifft(temp, axis=0)
    return cp.real(temp[:var.shape[0], ...]).ravel(order='F')


def _apply_filter_cl(options, var, filterAttr, filterSource, divide):
    """Same as _apply_filter_af/_apply_filter_cupy, plain PyOpenCL backend.
    PyOpenCL has no built-in array FFT, so the (host-cached) NumPy filter
    from `filterSource` is applied via a NumPy FFT round-trip: the device
    array is pulled to the host, filtered exactly like the CuPy path, and
    pushed back to the same command queue."""
    import numpy as np
    import pyopencl as cl
    filterNp = getattr(options, filterSource)
    queue = var.queue
    varNp = var.get()
    if options.subsets > 1 and (options.subsetType == 5 or options.subsetType == 4):
        if options.subsetType == 4:
            varNp = np.reshape(varNp, (options.nRowsD, varNp.size // options.nRowsD), order='F')
        else:
            varNp = np.reshape(varNp, (options.nColsD, varNp.size // options.nColsD), order='F')
    else:
        varNp = np.reshape(varNp, (options.nRowsD, options.nColsD, varNp.size // (options.nRowsD * options.nColsD)), order='F')
    temp = np.fft.fft(varNp, n=options.Nf, axis=0)
    filterShape = (-1,) + (1,) * (varNp.ndim - 1)
    if divide:
        temp /= filterNp.reshape(filterShape, order='F')
    else:
        temp *= filterNp.reshape(filterShape, order='F')
    temp = np.fft.ifft(temp, axis=0)
    out = np.real(temp[:varNp.shape[0], ...]).ravel(order='F').astype(np.float32)
    return cl.array.to_device(queue, out)


def applyMeasPreconditioning(options, var, subIter=None):
    """
    Computes the measurement-based preconditioning for the input data.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.
    var : arrayfire array, torch tensor, cupy array or pyopencl array
        The input data that is filtered.
    subIter : int, optional
        Zero-based subset/sub-iteration index used to index options.M when
        the diagonal (1 / (A1)) preconditioner (precondTypeMeas[0]) is
        enabled. Required only in that case, mirroring
        applyMeasPreconditioning.m's `subIter` argument (there 1-based).

    Returns
    -------
    var : arrayfire array, torch tensor, cupy array or pyopencl array
        The filtered input data.

    """
    if (options.precondTypeMeas[0].item() or options.precondTypeMeas[1].item()):
        if options.precondTypeMeas[1].item():
            if options.useAF:
                var = _apply_filter_af(options, var, 'filterG', 'filter0', divide=False)
            elif options.useTorch:
                var = _apply_filter_torch(options, var, 'filterG', 'filter0', divide=False)
            elif options.useCuPy:
                var = _apply_filter_cupy(options, var, 'filterG', 'filter0', divide=False)
            else:
                var = _apply_filter_cl(options, var, 'filterG', 'filter0', divide=False)
        if options.precondTypeMeas[0].item():
            if subIter is None:
                raise ValueError(
                    "applyMeasPreconditioning: 'subIter' must be provided when "
                    "options.precondTypeMeas[0] (diagonal 1 / (A1) preconditioner) is enabled."
                )
            var = var / options.M[subIter]

    return var

def circulantInverse(options, var):
    """
    Computes the circulant inverse for PDHG. Applies only when using the
    filtering-based preconditioner (above).

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.
    var : arrayfire array or torch tensor
        The partially computed dual estimate of the PDHG.

    Returns
    -------
    var : arrayfire array or torch tensor
        The fully computed dual estimate.

    """
    if options.useAF:
        var = _apply_filter_af(options, var, 'FilterG', 'Ffilter', divide=True)
    elif options.useTorch:
        var = _apply_filter_torch(options, var, 'FilterG', 'Ffilter', divide=True)
    elif options.useCuPy:
        var = _apply_filter_cupy(options, var, 'FilterG', 'Ffilter', divide=True)
    return var
