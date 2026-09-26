# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 11:40:44 2025

@author: Ville-Veikko Wettenhovi
"""

import numpy as np

def arc_correction(options, interpolate_sinogram):
    import copy
    import warnings
    from omegatomo.projector.detcoord import detectorCoordinates, sinogramCoordinates2D
    from omegatomo.util.matlabRound import matlabRound

    # Original (non-arc-corrected) LOR coordinates, computed with the TRUE
    # options (flip_image as given). These are what the non-arc projector
    # would use, and are also what the interpolation (section C) maps from.
    xp_o, yp_o = detectorCoordinates(options)
    x_o, y_o = sinogramCoordinates2D(options, xp_o, yp_o)

    # flip_image breaks the quadrant mirror-fill below (it assumes an
    # unflipped, quadrant-symmetric detector layout), so the arc-corrected
    # geometry (new_xp/new_yp, J, eka, vika, rotations) is built from a local
    # copy of options with flip_image forced to False. The returned x_o, y_o
    # above (used for interpolation) still use the true options.
    optionsGeom = copy.copy(options)
    optionsGeom.flip_image = False
    xp, yp = detectorCoordinates(optionsGeom)

    # Detector radius: the ideal circle used for the arc intersection has to
    # match the radius the detectors are actually placed at by
    # detectorCoordinates/computeCoordinates (diameter/2 + DOI), otherwise
    # the radial FOV shrinks by DOI.
    DOI = float(options.DOI) if np.size(options.DOI) > 0 else 0.0
    D = options.diameter + 2 * DOI

    new_xp = np.zeros_like(xp)
    new_yp = np.zeros_like(yp)

    # Detector angles
    if options.blocks_per_ring % 4 == 0:
        l_angles = np.linspace(0, 90, options.blocks_per_ring // 4 + 1)
    else:
        l_angles = np.linspace(0, 180, options.blocks_per_ring // 2 + 1)
        l_angles = l_angles[:int(np.ceil(options.blocks_per_ring / 4))]
    ll = 0  # Python index

    # cryst_per_block may be an array (multi-layer PET); the shift below
    # only ever needs the first layer's value.
    if isinstance(options.cryst_per_block, np.ndarray):
        cryst_per_block_val = options.cryst_per_block[0].item()
    else:
        cryst_per_block_val = options.cryst_per_block

    # Shift to zero angle. MATLAB rounds offangle away from zero and floors
    # the (pseudo-adjusted) crystals-per-block/2 term separately, which can
    # differ from a single combined int() truncation whenever cryst_per_block
    # is odd.
    if options.det_w_pseudo > options.det_per_ring:
        shift_val = int(matlabRound(options.offangle) + np.floor((cryst_per_block_val + 1) / 2))
    else:
        shift_val = int(matlabRound(options.offangle) + np.floor(cryst_per_block_val / 2))
    xp = np.roll(xp, -shift_val)
    yp = np.roll(yp, -shift_val)

    # xp -= options.diameter / 2
    # yp -= options.diameter / 2

    quarter_len = len(xp) // 4
    # Vectorized over kk: val_range is a fixed-width sliding window whose
    # start shifts by 1 with kk, so it is built as one (quarter_len, window)
    # gather-index matrix instead of being recomputed every iteration.
    # angle_ref (= l_angles[ll]) and the blocks_per_ring%4 branch below are
    # both invariant across kk (ll is never reassigned in the loop, and both
    # branches compute the same ind1), matching the original per-kk results
    # exactly.
    kk_idx = np.arange(quarter_len)
    r = int(options.det_w_pseudo / 2)
    if options.det_w_pseudo > options.det_per_ring:
        half_width = options.cryst_per_block + 1
    else:
        half_width = options.cryst_per_block
    offsets = np.arange(-half_width, half_width + 1)
    val_range_matrix = kk_idx[:, None] + r + offsets[None, :]

    angle_ref = l_angles[ll]
    dx = xp[kk_idx][:, None] - xp[val_range_matrix]
    dy = yp[kk_idx][:, None] - yp[val_range_matrix]
    angle1 = matlabRound(np.degrees(np.arctan2(dy, dx)) * 1e4) / 1e4

    diffs = np.abs(np.abs(angle1) - angle_ref)
    ind1 = np.argmin(diffs, axis=1)

    sel_idx = val_range_matrix[kk_idx, ind1]
    x2 = xp[sel_idx]
    y2 = yp[sel_idx]

    p = np.column_stack((xp[kk_idx], yp[kk_idx]))
    q = np.column_stack((x2, y2))
    d = q - p
    d_dot = np.sum(d * d, axis=1)
    p_dot = np.sum(p * p, axis=1)
    pd_dot = np.sum(2 * p * d, axis=1)
    rad2 = (D / 2) ** 2

    sqrt_term = pd_dot ** 2 - 4 * d_dot * (p_dot - rad2)
    l_param = (-pd_dot - np.sqrt(sqrt_term)) / (2 * d_dot)

    lx = p + l_param[:, None] * d
    new_xp[:quarter_len] = lx[:, 0]
    new_yp[:quarter_len] = lx[:, 1]

    # Mirror fill
    new_xp[:quarter_len] += D / 2
    new_yp[:quarter_len] += D / 2

    new_yp[quarter_len:2*quarter_len] = np.flip(new_yp[:quarter_len])
    diffi = np.diff(np.concatenate([np.flip(new_yp[:quarter_len]), [D/2]]))
    diffi = D/2 + np.cumsum(np.flip(diffi))
    new_yp[2*quarter_len:] = np.concatenate([diffi, np.flip(diffi)])

    diffi = np.diff(np.concatenate([new_xp[:quarter_len], [D/2]]))
    diffi = D/2 + np.cumsum(np.flip(diffi))
    new_xp[quarter_len:2*quarter_len] = diffi
    new_xp[2*quarter_len:] = np.flip(new_xp[:2*quarter_len])

    # Undo shift
    xp = np.roll(new_xp, shift_val)
    yp = np.roll(new_yp, shift_val)

    x, y = sinogramCoordinates2D(optionsGeom, xp, yp)

    # Reshape LOR coordinates
    xx1 = x[:, 0].reshape(options.Ndist, options.Nang, order='F')
    xx2 = x[:, 1].reshape(options.Ndist, options.Nang, order='F')
    yy1 = y[:, 0].reshape(options.Ndist, options.Nang, order='F')
    yy2 = y[:, 1].reshape(options.Ndist, options.Nang, order='F')

    # Compute LOR angles (degrees), shift to avoid negatives, find J
    angle = np.degrees(np.arctan((yy1 - yy2) / (xx1 - xx2))) + 90
    angle[angle == 180] = 0
    J = np.argmin(np.mean(angle, axis=0))

    # Create arc-corrected coordinates from perpendicular LOR
    eka = xx1[0, J]
    vika = xx1[-1, J]
    alkux = np.linspace(eka, vika, options.Ndist)

    alkuy = np.sqrt((D / 2)**2 - (alkux - D / 2)**2) + D / 2
    alku = np.asfortranarray(np.vstack([alkux, alkuy]))
    alku2 = alku - D / 2

    alku = np.asfortranarray(np.vstack([alkux, np.abs(alkuy - D)]))
    alku1 = alku - D / 2

    # Compute angles
    angles = np.linspace(0, 180, options.Nang + 1)[:-1]
    angles = np.roll(angles, J + 1)
    angles = angles.reshape(1, 1, -1)

    # Rotation matrices
    cos_a = np.cos(np.radians(angles)).squeeze()
    sin_a = np.sin(np.radians(angles)).squeeze()
    rot_matrix = np.array([[cos_a, -sin_a], [sin_a, cos_a]], order='F')  # shape: (2, 2, Nang)

    # Rotate all points with each matrix (broadcasting).
    # out[n, l, a] = sum_b rot_matrix[a, b, n] * alku[b, l], i.e. exactly
    # (rot_matrix[:, :, n] @ alku).T[l, a] for every n -- a batched version
    # of the per-ii 2x2 matmul below, reshaped (Nang, Ndist, 2) -> (Nang *
    # Ndist, 2) in the same n-major, l-minor order the loop filled.
    new_xy1 = (np.einsum('abn,bl->nla', rot_matrix, alku1) + D / 2).reshape(-1, 2)
    new_xy2 = (np.einsum('abn,bl->nla', rot_matrix, alku2) + D / 2).reshape(-1, 2)
    new_xy1 = np.asfortranarray(new_xy1)
    new_xy2 = np.asfortranarray(new_xy2)

    # Multiply each matrix with alku1 and alku2 (2 x Ndist vectors)
    # Result will be shape (options.Nang, 2) after transpose
    # new_xy1 = np.array([(R @ alku1).T for R in rot_matrix]) + D / 2
    # new_xy2 = np.array([(R @ alku2).T for R in rot_matrix]) + D / 2
    # Update x and y arrays
    x = np.column_stack([new_xy1[:, 0], new_xy2[:, 0]])
    y = np.column_stack([new_xy1[:, 1], new_xy2[:, 1]])

    # Flip half of the LORs to match symmetry (MATLAB-style)
    xx1 = x[:, 0].reshape(options.Ndist, options.Nang, order='F')
    xx2 = x[:, 1].reshape(options.Ndist, options.Nang, order='F')
    yy1 = y[:, 0].reshape(options.Ndist, options.Nang, order='F')
    yy2 = y[:, 1].reshape(options.Ndist, options.Nang, order='F')

    apu = xx1[:, :J+1].copy()
    xx1[:, :J+1] = np.flipud(xx2[:, :J+1])
    xx2[:, :J+1] = np.flipud(apu)

    apu = yy1[:, :J+1].copy()
    yy1[:, :J+1] = np.flipud(yy2[:, :J+1])
    yy2[:, :J+1] = np.flipud(apu)

    # Flatten back to 2D coordinates
    x[:, 0] = xx1.ravel('F')
    x[:, 1] = xx2.ravel('F')
    y[:, 0] = yy1.ravel('F')
    y[:, 1] = yy2.ravel('F')

    x -= D / 2
    y -= D / 2

    # Equidistance check: within each angular column, the signed radial
    # distance between adjacent LORs should be (near) constant. A failure
    # usually means the geometry construction above could not be made
    # self-consistent for this detector configuration.
    # Unfolded signed distance (not _lorParams, which folds the angle into
    # [0, pi) and flips the sign of s along with it): a column whose LOR
    # angle lies exactly at the 0/pi fold would otherwise have its sign
    # flip mid-column, producing a false equidistance warning (seen for
    # the oddcpb geometry). Endpoint order (columns 0/1 of x, y) is
    # consistent within each column.
    s_check = (x[:, 0] * y[:, 1] - x[:, 1] * y[:, 0]) / np.hypot(x[:, 1] - x[:, 0], y[:, 1] - y[:, 0])
    sMat = s_check.reshape(options.Ndist, options.Nang, order='F')
    dS = np.diff(sMat, axis=0)
    if dS.size > 0:
        spreadPerCol = np.max(dS, axis=0) - np.min(dS, axis=0)
        medAbs = np.median(np.abs(dS))
        if medAbs > 0 and np.max(spreadPerCol) > 1e-3 * medAbs:
            warnings.warn('Arc correction failed to make all the LORs equidistant.')

    if interpolate_sinogram:
        W = arcInterpWeights(x_o, y_o, x, y, options.Ndist, options.Nang, options.arc_interpolation)
        if W is None:
            u_o, v_o, u_n, v_n = _lorUV(x_o, y_o, x, y, options.Ndist, options.Nang)
        else:
            u_o = v_o = u_n = v_n = None

        def _apply(A):
            if isinstance(A, list):
                return [_apply(a) for a in A]
            if W is not None:
                return _applyArc(A, W, options.Ndist, options.Nang)
            return _applyArcFallback(A, u_o, v_o, u_n, v_n, options.Ndist, options.Nang, options.arc_interpolation)

        def _hasData(A):
            if isinstance(A, list):
                return len(A) > 0
            return hasattr(A, 'size') and A.size > 1

        options.SinM = _apply(options.SinM)

        # Also interpolate any corrections that are still stored as separate
        # sinogram-shaped arrays and applied during reconstruction, so they
        # stay in the same (arc-corrected) LOR geometry as SinM/x/y. When a
        # correction has already been baked into SinM directly (precorrect,
        # not corrections_during_reconstruction), it is left untouched here.
        if options.normalization_correction and options.corrections_during_reconstruction and _hasData(options.normalization):
            options.normalization = _apply(options.normalization)
        if _hasData(options.SinDelayed):
            options.SinDelayed = _apply(options.SinDelayed)
        if options.additionalCorrection and _hasData(options.corrVector):
            options.corrVector = _apply(options.corrVector)
        if options.scatter_correction and _hasData(options.ScatterC):
            options.ScatterC = _apply(options.ScatterC)

        if options.verbose:
            print("Arc correction interpolation complete")

    return x, y, options


def arcInterpWeights(x_o, y_o, x, y, Ndist, Nang, method='linear'):
    """
    Builds the sparse interpolation-weight matrix that maps sinogram data in
    the original (non-arc-corrected) LOR geometry onto the arc-corrected LOR
    geometry: new_values = W @ old_values.

    Parameters
    ----------
    x_o : NumPy array
        x-coordinates (N x 2) of the original LORs.
    y_o : NumPy array
        y-coordinates (N x 2) of the original LORs.
    x : NumPy array
        x-coordinates (N x 2) of the arc-corrected LORs.
    y : NumPy array
        y-coordinates (N x 2) of the arc-corrected LORs.
    Ndist : int
        Number of radial distances.
    Nang : int
        Number of angles.
    method : str, optional
        Interpolation method. 'linear' and 'nearest' are supported directly
        (a sparse matrix is returned); anything else returns None and the
        caller has to fall back to a slower, generic interpolation (e.g.
        scipy.interpolate.griddata). The default is 'linear'.

    Returns
    -------
    scipy.sparse.csr_matrix or None
        Sparse (N x N) matrix with rows = new LORs, columns = original LORs.
        Every row sums to 1. None if ``method`` is not 'linear' or 'nearest'.

    """
    if method not in ('linear', 'nearest'):
        return None

    from scipy.sparse import coo_matrix
    from scipy.spatial import Delaunay, cKDTree

    N = x_o.shape[0]
    u_o, v_o, u_n, v_n = _lorUV(x_o, y_o, x, y, Ndist, Nang)

    P, srcIdx = _periodicExtend(u_o, v_o, Nang)
    uniquePts, invIdx, counts, starts, sortedSrc = _uniqueGroups(P, srcIdx)

    Q = np.column_stack((u_n, v_n))

    if method == 'linear':
        tri = Delaunay(uniquePts)
        simplex = tri.find_simplex(Q)
        inside = simplex >= 0

        rows_list = []
        cols_list = []
        data_list = []

        if np.any(inside):
            simpIn = simplex[inside]
            Tinv = tri.transform[simpIn, :2]
            r = Q[inside] - tri.transform[simpIn, 2]
            b = np.einsum('ijk,ik->ij', Tinv, r)
            bary = np.column_stack((b, 1.0 - b.sum(axis=1)))
            # Clip tiny negative barycentric weights (numerical noise) to 0
            bary = np.where((bary < 0) & (bary > -1e-12), 0.0, bary)
            verts = tri.simplices[simpIn]
            rowIdxIn = np.nonzero(inside)[0]
            for k in range(3):
                rows_list.append(rowIdxIn)
                cols_list.append(verts[:, k])
                data_list.append(bary[:, k])

        outside = ~inside
        if np.any(outside):
            tree = cKDTree(uniquePts)
            _, nn = tree.query(Q[outside])
            rows_list.append(np.nonzero(outside)[0])
            cols_list.append(nn)
            data_list.append(np.ones(nn.shape[0]))

        rowsU = np.concatenate(rows_list)
        colsU = np.concatenate(cols_list)
        dataU = np.concatenate(data_list)
    else:  # 'nearest'
        tree = cKDTree(uniquePts)
        _, nn = tree.query(Q)
        rowsU = np.arange(Q.shape[0])
        colsU = nn
        dataU = np.ones(nn.shape[0])

    # Expand each (row, uniquePoint, weight) triple into its duplicate
    # original columns (periodic copies that collapsed onto the same (u,v)
    # point), splitting the weight equally between them -- the inverse of
    # griddata's duplicate averaging.
    cnt = counts[colsU]
    maxCnt = int(cnt.max())
    offsets = np.arange(maxCnt)
    posGrid = starts[colsU][:, None] + offsets[None, :]
    maskGrid = offsets[None, :] < cnt[:, None]
    flatPos = posGrid[maskGrid]
    finalCols = sortedSrc[flatPos]
    finalRows = np.repeat(rowsU, cnt)
    finalData = np.repeat(dataU, cnt) / np.repeat(cnt, cnt)

    W = coo_matrix((finalData, (finalRows, finalCols)), shape=(Q.shape[0], N)).tocsr()
    W.sum_duplicates()
    return W


def _lorParams(xArr, yArr):
    """
    Canonical (angle, signed distance) LOR parameters, folded into the
    [0, pi) / consistent-sign convention used for both the original and the
    arc-corrected LOR sets.
    """
    x1 = xArr[:, 0]
    y1 = yArr[:, 0]
    x2 = xArr[:, 1]
    y2 = yArr[:, 1]
    dx = x2 - x1
    dy = y2 - y1
    th = np.arctan2(dy, dx)
    s = (x1 * y2 - x2 * y1) / np.hypot(dx, dy)
    neg = th < 0
    th = np.where(neg, th + np.pi, th)
    s = np.where(neg, -s, s)
    ge = th >= np.pi
    th = np.where(ge, th - np.pi, th)
    s = np.where(ge, -s, s)
    return th, s


def _lorUV(x_o, y_o, x, y, Ndist, Nang):
    """
    Converts the canonical LOR parameters of the original and arc-corrected
    sets to (angular bin, radial bin) units, using the radial bin size
    computed from the ORIGINAL set for both.
    """
    th_o, s_o = _lorParams(x_o, y_o)
    th_n, s_n = _lorParams(x, y)
    ds = (np.max(s_o) - np.min(s_o)) / (Ndist - 1)
    u_o = th_o / (np.pi / Nang)
    v_o = s_o / ds
    u_n = th_n / (np.pi / Nang)
    v_n = s_n / ds
    return u_o, v_o, u_n, v_n


def _periodicExtend(u_o, v_o, Nang):
    """
    Periodic extension of the original (u, v) points (angle is pi-periodic
    with a sign flip on v), needed so the triangulation/nearest-neighbour
    search below sees continuous data across the angle wrap-around.
    """
    N = u_o.shape[0]
    P = np.column_stack((
        np.concatenate((u_o, u_o - Nang, u_o + Nang)),
        np.concatenate((v_o, -v_o, -v_o))
    ))
    srcIdx = np.tile(np.arange(N), 3)
    return P, srcIdx


def _uniqueGroups(P, srcIdx):
    """
    Rounds the periodically-extended points to remove floating-point noise
    and groups the (possibly duplicated) points, returning everything needed
    to both triangulate the unique points and later map back to the original
    (non-extended) column indices.
    """
    Pr = np.round(P * 1e6) / 1e6
    uniquePts, invIdx = np.unique(Pr, axis=0, return_inverse=True)
    invIdx = np.asarray(invIdx).ravel()
    order = np.argsort(invIdx, kind='stable')
    counts = np.bincount(invIdx, minlength=uniquePts.shape[0])
    starts = np.concatenate(([0], np.cumsum(counts)[:-1]))
    sortedSrc = srcIdx[order]
    return uniquePts, invIdx, counts, starts, sortedSrc


def _applyArc(A, W, Ndist, Nang):
    """
    Applies the sparse interpolation matrix W to A. A can be of any shape
    whose total number of elements is a positive multiple of Ndist*Nang
    (trailing dimensions -- slices, TOF bins, dynamic frames, ...-- are
    handled generically); otherwise A is returned unchanged.
    """
    A = np.asarray(A)
    n = Ndist * Nang
    if A.size <= 0 or A.size % n != 0:
        return A
    origShape = A.shape
    origDtype = A.dtype
    flat = A.reshape((n, -1), order='F')
    result = W @ flat.astype(np.float64)
    outDtype = np.float32 if (np.issubdtype(origDtype, np.integer) or origDtype == np.float32) else origDtype
    result = result.astype(outDtype)
    return result.reshape(origShape, order='F')


def _applyArcFallback(A, u_o, v_o, u_n, v_n, Ndist, Nang, method):
    """
    Generic (slow) fallback used when arcInterpWeights returns None, i.e.
    for interpolation methods other than 'linear'/'nearest' (e.g. 'natural',
    'cubic', 'v4'). Mirrors _applyArc's shape handling.
    """
    A = np.asarray(A)
    n = Ndist * Nang
    if A.size <= 0 or A.size % n != 0:
        return A
    origShape = A.shape
    origDtype = A.dtype
    flat = A.reshape((n, -1), order='F')
    result = _griddataFallbackApply(flat, u_o, v_o, u_n, v_n, Ndist, Nang, method)
    outDtype = np.float32 if (np.issubdtype(origDtype, np.integer) or origDtype == np.float32) else origDtype
    result = result.astype(outDtype)
    return result.reshape(origShape, order='F')


def _griddataFallbackApply(flat, u_o, v_o, u_n, v_n, Ndist, Nang, method):
    """
    Per-column scipy.interpolate.griddata fallback on the same periodic
    (u, v) points used by arcInterpWeights. The point set (positions) is
    built once and reused for every column; only the (duplicate-averaged)
    values change per column.
    """
    from scipy.interpolate import griddata
    P, _ = _periodicExtend(u_o, v_o, Nang)
    Pr = np.round(P * 1e6) / 1e6
    uniquePts, invIdx = np.unique(Pr, axis=0, return_inverse=True)
    invIdx = np.asarray(invIdx).ravel()
    nUnique = uniquePts.shape[0]
    counts = np.bincount(invIdx, minlength=nUnique)

    Q = np.column_stack((u_n, v_n))
    result = np.empty((Q.shape[0], flat.shape[1]), dtype=np.float64)
    for c in range(flat.shape[1]):
        vals3 = np.tile(flat[:, c].astype(np.float64), 3)
        sumVal = np.bincount(invIdx, weights=vals3, minlength=nUnique)
        avgVal = sumVal / counts
        interp = griddata(uniquePts, avgVal, Q, method=method, fill_value=0.0)
        interp = np.nan_to_num(interp, nan=0.0)
        result[:, c] = interp
    return result
