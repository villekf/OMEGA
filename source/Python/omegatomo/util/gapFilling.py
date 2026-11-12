# -*- coding: utf-8 -*-
"""
Sinogram gap filling for PET scanners with pseudo detectors.

Python port of gapFilling.m and applyGapFilling.m (MATLAB). If pseudo
detectors are used (e.g. Biograph mCT), these detectors don't actually receive
any measurements. This causes zero measurements when no mashing and reduced
number of measurements with mashing. The gap filling interpolation attempts to
fix these gaps.

The method 'fillmissing' (MATLAB's fillmissing, applied separately along the
radial and the angular dimension of each sinogram and averaged) is supported.
The method 'inpaint_nans' (File Exchange function) and 'fillmissing2' are not.
"""

import numpy as np


# Interpolation methods of MATLAB's fillmissing that are supported
_FILLMISSING_METHODS = ('linear', 'nearest', 'previous', 'next', 'pchip', 'spline', 'makima')


def _scalar(value):
    # First element of a scalar or array-like parameter as a Python number
    return np.asarray(value).ravel()[0].item()


def _gapMask(options):
    """
    Computes the location of the sinogram gaps (the bins that are affected by
    the pseudo detectors). Port of the mask computation in gapFilling.m, using
    the LOR to sinogram bin indexing of sinogram_coordinates_2D.m.

    Parameters
    ----------
    options : class object
        OMEGA class object. Uses det_w_pseudo, Nang, Ndist, cryst_per_block,
        ndist_side, normalization_correction and
        corrections_during_reconstruction.

    Returns
    -------
    gaps : NumPy array (bool)
        Array of size (Ndist, Nang), True for the bins that are filled.

    """
    det = int(_scalar(options.det_w_pseudo))
    Nang = int(_scalar(options.Nang))
    Ndist = int(_scalar(options.Ndist))
    cryst = int(_scalar(options.cryst_per_block))
    nDistSide = int(_scalar(options.ndist_side))
    if det % (2 * Nang) != 0 or det < 2 * Nang:
        raise ValueError('The number of detectors with pseudo detectors (' + str(det) + ') divided by two must be a multiple of the number of ' +
                         'sinogram angles (' + str(Nang) + ') for gap filling!')
    mashing = det // Nang // 2
    # The mask is computed at the full angular resolution (before mashing)
    NangFull = Nang * mashing
    half = det // 2

    # Detector pair indices (1-based, ordered by the first, then the second detector). Same as formDetectorIndices.m
    ya, xa = np.triu_indices(det)
    L = np.column_stack((ya + 1, xa + 1))
    # sinogram_coordinates_2D.m works with zero-based indices (xa = max, ya = min)
    # Angle
    j = ((xa + ya + half + 1) % det) // 2
    b = j + half
    # Distance
    i = np.abs(xa - ya - half)
    neg = (ya < j) | (b < xa)
    i[neg] = -i[neg]
    # Accepted LORs (distances that are within the predefined value)
    if Ndist % 2 == 0:
        accepted = (i <= (Ndist // 2 + min(0, nDistSide))) & (i >= (-(Ndist // 2) + max(0, nDistSide)))
    else:
        accepted = (i <= Ndist // 2) & (i >= -(Ndist // 2))
    j = j // (half // NangFull)
    i = i[accepted]
    j = j[accepted]
    L = L[accepted, :]
    # MATLAB one-based bin numbers
    if np.min(i) <= 0:
        i = i + np.abs(np.min(i)) + 1
    j = j + 1
    if np.max(i) > Ndist or np.max(j) > NangFull:
        raise ValueError('The sinogram bin indices exceed the sinogram size (Ndist = ' + str(Ndist) + ', Nang = ' + str(NangFull) + ')!')
    lin = (i - 1) + Ndist * (j - 1)

    # Pseudo detector locations (the detector number is divisible by cryst_per_block + 1)
    nBins = Ndist * NangFull
    ff = np.bincount(lin, weights=(L[:, 0] % (cryst + 1) != 0).astype(np.float64), minlength=nBins)
    ff2 = np.bincount(lin, weights=(L[:, 1] % (cryst + 1) != 0).astype(np.float64), minlength=nBins)
    uusf = ff + ff2
    uusf[uusf == 1] = 0
    uusf[uusf > 0] = 1
    apu = (uusf == 0).astype(np.float64).reshape((Ndist, NangFull), order='F')
    if mashing > 1:
        apu = np.mean(apu.reshape((Ndist, Nang, mashing)), axis=2)
    normCorr = bool(getattr(options, 'normalization_correction', False))
    corrRecon = bool(getattr(options, 'corrections_during_reconstruction', True))
    if normCorr and not corrRecon and mashing > 1:
        gaps = apu >= 1
    else:
        gaps = apu > 0
    return gaps


def _interpolateColumn(kx, y, mx, method):
    """
    Interpolates the missing samples mx from the known samples (kx, y) of one
    column. y has shape (len(kx), S) with S sinograms. Returns (values, valid)
    where valid marks the entries of mx that were filled (MATLAB's fillmissing
    leaves e.g. the leading entries of 'previous' and all entries of
    two-point methods with a single known sample as NaN).
    """
    nk = kx.size
    nm = mx.size
    out = np.full((nm, y.shape[1]), np.nan, dtype=np.float64)
    valid = np.zeros(nm, dtype=bool)
    if nk == 0:
        return out, valid
    r = np.searchsorted(kx, mx)
    if method == 'nearest':
        left = np.clip(r - 1, 0, nk - 1)
        right = np.clip(r, 0, nk - 1)
        # Ties are resolved towards the later sample (as MATLAB does)
        useRight = (kx[right] - mx) <= (mx - kx[left])
        useRight[r >= nk] = False
        useRight[r < 1] = True
        idx = np.where(useRight, right, left)
        out[:] = y[idx, :]
        valid[:] = True
    elif method == 'previous':
        valid = r >= 1
        out[valid, :] = y[r[valid] - 1, :]
    elif method == 'next':
        valid = r < nk
        out[valid, :] = y[r[valid], :]
    elif nk >= 2:
        # Interpolation with extrapolation at the ends (MATLAB's fillmissing default EndValues = 'extrap')
        if method == 'linear':
            lo = np.clip(r - 1, 0, nk - 2)
            t = ((mx - kx[lo]) / (kx[lo + 1] - kx[lo]))[:, None]
            out[:] = y[lo, :] * (1. - t) + y[lo + 1, :] * t
        elif method == 'pchip':
            from scipy.interpolate import PchipInterpolator
            out[:] = PchipInterpolator(kx, y, axis=0, extrapolate=True)(mx)
        elif method == 'spline':
            from scipy.interpolate import CubicSpline
            out[:] = CubicSpline(kx, y, axis=0, bc_type='not-a-knot', extrapolate=True)(mx)
        else:
            from scipy.interpolate import Akima1DInterpolator
            out[:] = Akima1DInterpolator(kx, y, axis=0, method='makima', extrapolate=True)(mx)
        valid[:] = True
    return out, valid


def _fillMissingDim0(A, method):
    """
    Equivalent of MATLAB's fillmissing(A, method) applied to each column (the
    first dimension) of each (n, m) slice of the (n, m, S) float array A.
    A has NaN at the missing entries and is not modified. The NaN pattern is
    assumed to be the same in all slices S (the caller guarantees it).
    """
    out = A.copy()
    miss = np.isnan(A[:, :, 0])
    for jj in range(A.shape[1]):
        col = miss[:, jj]
        if not col.any():
            continue
        kx = np.flatnonzero(~col)
        mx = np.flatnonzero(col)
        vals, valid = _interpolateColumn(kx, A[kx, jj, :].astype(np.float64), mx, method)
        out[mx[valid], jj, :] = vals[valid, :]
    return out


def _fillMissing(A, method):
    """
    fillmissing along the first dimension for the (n, m, S) array A. If the
    NaN pattern differs between the slices, each slice is processed separately.
    """
    miss = np.isnan(A)
    if A.shape[2] == 1 or np.array_equal(miss, np.broadcast_to(miss[:, :, :1], miss.shape)):
        return _fillMissingDim0(A, method)
    out = np.empty_like(A)
    for kk in range(A.shape[2]):
        out[:, :, kk:kk + 1] = _fillMissingDim0(A[:, :, kk:kk + 1], method)
    return out


def gapFilling(options, Sin, llo = 1, gaps = None):
    """
    Performs gap filling to pseudo detector sinograms (gapFilling.m).

    Parameters
    ----------
    options : class object
        OMEGA class object. Uses fill_sinogram_gaps related options
        gap_filling_method ('fillmissing') and
        interpolation_method_fillmissing ('linear', 'nearest', 'previous',
        'next', 'pchip', 'spline' or 'makima'), as well as the sinogram and
        scanner properties (det_w_pseudo, Nang, Ndist, cryst_per_block,
        ndist_side).
    Sin : NumPy array
        The original sinogram. The first two dimensions are the radial
        (Ndist) and the angular (Nang) dimensions, the remaining dimensions
        (sinograms, TOF bins, time steps, ...) are all processed
        independently.
    llo : int, optional
        Current time step (1 is the first). The gap locations are computed
        for the first time step. The default is 1.
    gaps : NumPy array (bool), optional
        Location of the gaps to be filled (output of a previous call, used
        for subsequent time steps).

    Returns
    -------
    Sin : NumPy array
        The gap-filled sinogram (float32).
    gaps : NumPy array (bool)
        The location of the sinogram gaps for subsequent time steps.

    """
    import warnings
    if options.verbose:
        print('Starting sinogram gap filling. This can take several minutes.')
    if int(_scalar(getattr(options, 'nLayers', 1))) > 1:
        raise NotImplementedError('Sinogram gap filling is not supported with multiple crystal layers!')
    if llo == 1 or gaps is None or np.size(gaps) == 0:
        gaps = _gapMask(options)
    method = options.gap_filling_method
    if method == 'fillmissing':
        interp = options.interpolation_method_fillmissing
        if interp not in _FILLMISSING_METHODS:
            raise ValueError('Unsupported interpolation method for fillmissing: ' + str(interp) + '! Supported methods: ' + ', '.join(_FILLMISSING_METHODS))
        Sin = np.asarray(Sin)
        if Sin.ndim < 2 or Sin.shape[0] != gaps.shape[0] or Sin.shape[1] != gaps.shape[1]:
            raise ValueError('Size mismatch between the sinogram ' + str(Sin.shape) + ' and the gap locations ' + str(gaps.shape) + '!')
        shape = Sin.shape
        # (Ndist, Nang, all the other dimensions)
        A = np.array(Sin, dtype=np.float32, order='F').reshape((shape[0], shape[1], -1), order='F')
        A = np.ascontiguousarray(A)
        A[gaps, :] = np.nan
        # Interpolation along the radial dimension and the angular dimension, the results are averaged
        jelppi1 = _fillMissing(A, interp).astype(np.float32)
        jelppi2 = np.swapaxes(_fillMissing(np.ascontiguousarray(np.swapaxes(A, 0, 1)), interp), 0, 1).astype(np.float32)
        A = (jelppi1 + jelppi2) / np.float32(2.)
        A[A < 0] = 0
        Sin = np.asfortranarray(A.reshape(shape, order='F'))
    elif method == 'inpaint_nans':
        raise NotImplementedError('The gap filling method inpaint_nans is not available in Python. inpaint_nans is a MATLAB File Exchange function. ' +
                                  'Use options.gap_filling_method = "fillmissing" instead.')
    elif method == 'fillmissing2':
        raise NotImplementedError('The gap filling method fillmissing2 (2-D interpolation) is not available in Python. ' +
                                  'Use options.gap_filling_method = "fillmissing" instead.')
    else:
        warnings.warn('Unsupported gap filling method! No gap filling was performed!')
    if options.verbose:
        print('Gap filling completed')
    return Sin, gaps


def applyGapFilling(options):
    """
    Performs sinogram gap filling for the input data options.SinM (applyGapFilling.m). SinM can be a
    single sinogram array or a list with one array per time step.

    Parameters
    ----------
    options : class object
        OMEGA class object.

    Returns
    -------
    options : class object
        The same object with the gap-filled options.SinM.

    """
    gaps = None
    if isinstance(options.SinM, (list, tuple)):
        SinM = list()
        for llo in range(len(options.SinM)):
            Sin, gaps = gapFilling(options, options.SinM[llo], llo + 1, gaps)
            SinM.append(Sin)
        options.SinM = SinM
    else:
        options.SinM, gaps = gapFilling(options, options.SinM, 1, gaps)
    return options
