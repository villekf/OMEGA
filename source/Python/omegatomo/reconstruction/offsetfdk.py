# -*- coding: utf-8 -*-
"""
Offset-detector (half-fan) CBCT FDK preprocessing: redundancy weighting with conjugate completion.

The FDK reconstruction of an offset detector has to weight the redundant (conjugate) rays BEFORE the ramp filtering and the filtered
projection needs data beyond the physical edge of the detector on the missing side. Weighting inside the backprojection kernel (after the
filtering) leaves a bright central artifact. The preprocessing implemented here works on the line integrals:

1. the central ray of every projection is computed exactly from the geometry,
2. the missing side of the detector is completed with virtual detector rows. Each virtual pixel gets the value of the measured conjugate ray
   (theta + pi, -rho) of the same detector column, linearly interpolated in the (theta, rho) space of the measured rays,
3. the physical and the virtual rows are multiplied with a smooth redundancy weight w(s) with the half-width of the transition region
   Delta (options.offsetFDKDelta, wider than the physical overlap) and with 2 (the weights of the conjugate rays sum to one, but the
   full-scan FDK of the backprojector expects two),
4. the detector is padded with zeros after the virtual rows (the ramp filter needs a zero tail) and its centre is shifted so that the physical rows
   keep their exact world positions.

The library FDK is then called with offsetCorrection = False.
"""
import numpy as np

# Number of the last reliable physical rows (counting from the missing side) that the hold value of the virtual rows outside the interpolation
# hull is computed from: rows N-8 ... N-4 for the missing side being the high-index side
_HOLD_FIRST = 8
_HOLD_LAST = 4
# Margin (pixels) kept between the transition region and the far edge of the detector and number of extra rows in the interpolation window
_EDGE_MARGIN = 2.
_WINDOW_MARGIN = 16
# Zero rows beyond the virtual rows (the zero tail of the ramp filter)
_ZERO_TAIL = 256
# The part of the periodic theta axis that is duplicated over the seam (radians)
_THETA_WRAP = 0.5

# The attributes of the options that offsetFDKPreprocess changes and offsetFDKRestore puts back
_SAVED_ATTRIBUTES = ('nRowsD', 'Ndist', 'det_per_ring', 'det_w_pseudo', 'cryst_per_block', 'detectors', 'nTotMeas', 'nMeasSubset', 'dSizeXBP',
                     'x', 'offsetCorrection', 'OffsetLimit', 'SinM')


def _centralRay(x, u, nRows):
    """Position of the central ray (the ray through the rotation axis) of every projection in pixels from the row-0 edge of the detector.

    x : (K, 6) source and detector centre coordinates, u : (K, >=2) the detector row direction vector with the length of one pixel.
    Same closed form as OffsetLimit (divided by the pixel size) in setUpCorrections."""
    s = x[:, 0:2]
    d = x[:, 3:5]
    t = (s[:, 1] * d[:, 0] - s[:, 0] * d[:, 1]) / (s[:, 0] * u[:, 1] - s[:, 1] * u[:, 0])
    return t + nRows / 2.


def _weightProfile(nRows, T, delta):
    """The redundancy weight w_Delta(s) (before the factor of two) at the pixel centres of the first nRows rows, shape (nRows, K) float64.

    s = (i + 0.5) - T is the signed distance from the central ray in pixels, positive towards the missing side. w = 1 for s < -Delta, 0 for
    s > Delta and 0.5 * (1 - sin(pi * s / (2 * Delta))) in between, i.e. w(s) + w(-s) = 1."""
    s = (np.arange(nRows, dtype=np.float64)[:, None] + 0.5) - np.asarray(T, dtype=np.float64)[None, :]
    return np.where(s < -delta, 1., np.where(s > delta, 0., 0.5 * (1. - np.sin(np.pi * s / (2. * delta)))))


def _rayParameters(x, u, rows, nRows):
    """The transaxial ray parameters of the rays through the centres of the detector rows `rows` (indices, may be beyond the physical detector):
    the direction angle theta = atan2 of the ray direction and the signed distance rho = s_x e_y - s_y e_x of the ray from the rotation
    axis. Both are (K, len(rows))."""
    t = np.asarray(rows, dtype=np.float64)[None, :] + 0.5 - nRows / 2.
    qx = x[:, 3, None] + t * u[:, 0, None]
    qy = x[:, 4, None] + t * u[:, 1, None]
    ex = qx - x[:, 0, None]
    ey = qy - x[:, 1, None]
    nn = np.hypot(ex, ey)
    ex = ex / nn
    ey = ey / nn
    return np.arctan2(ey, ex), x[:, 0, None] * ey - x[:, 1, None] * ex


def _interpolationMatrix(x, u, nRows, nVirtual, nWindow, pitch, nProj):
    """Sparse linear interpolation from the measured rays to the conjugate rays of the virtual pixels.

    The measured rays are those of the last nWindow rows of all the projections, the queries are the conjugate rays (theta + pi, -rho) of
    the nVirtual virtual rows N, ..., N + nVirtual - 1. The interpolation is the piecewise linear one on the Delaunay triangulation of the
    measured rays in the (theta / (2 pi / K), rho / (pitch * SOD / SDD)) plane, with theta periodic.

    Returns S : (K * nVirtual, K * nWindow) scipy.sparse CSR (float32) with the row index k * nVirtual + e (projection k, virtual row e) and
    the column index k * nWindow + r (projection k, measured row N - nWindow + r), and the boolean (K, nVirtual) array of the queries that
    are outside the hull of the measured rays (their rows of S are empty)."""
    from scipy.spatial import Delaunay
    from scipy.sparse import coo_matrix
    sod = np.hypot(x[:, 0], x[:, 1]).mean()
    sdd = np.linalg.norm(x[:, 0:3] - x[:, 3:6], axis=1).mean()
    scTheta = 2. * np.pi / nProj
    scRho = pitch * sod / sdd
    th, rho = _rayParameters(x, u, np.arange(nRows - nWindow, nRows), nRows)
    thF = th.ravel()
    rhoF = rho.ravel()
    nPts = thF.size
    lo = thF < -np.pi + _THETA_WRAP
    hi = thF > np.pi - _THETA_WRAP
    src = np.concatenate((np.arange(nPts), np.nonzero(lo)[0], np.nonzero(hi)[0]))
    pts = np.column_stack((np.concatenate((thF, thF[lo] + 2. * np.pi, thF[hi] - 2. * np.pi)) / scTheta,
                           np.concatenate((rhoF, rhoF[lo], rhoF[hi])) / scRho))
    # Coincident points are not allowed in the triangulation (the first occurrence is kept)
    _, first = np.unique(pts, axis=0, return_index=True)
    first = np.sort(first)
    pts = pts[first]
    src = src[first]
    tri = Delaunay(pts)
    thq, rhoq = _rayParameters(x, u, np.arange(nRows, nRows + nVirtual), nRows)
    thc = np.mod(thq, 2. * np.pi) - np.pi
    q = np.column_stack((thc.ravel() / scTheta, (-rhoq).ravel() / scRho))
    simplex = tri.find_simplex(q)
    valid = simplex >= 0
    ids = np.nonzero(valid)[0]
    sv = simplex[ids]
    tr = tri.transform[sv]
    b = np.einsum('nij,nj->ni', tr[:, :2, :], q[ids] - tr[:, 2, :])
    bary = np.column_stack((b, 1. - b.sum(axis=1)))
    good = np.all(np.isfinite(bary), axis=1)
    ids = ids[good]
    bary = bary[good]
    verts = tri.simplices[sv[good]]
    inside = np.zeros(q.shape[0], dtype=bool)
    inside[ids] = True
    S = coo_matrix((bary.ravel().astype(np.float32), (np.repeat(ids, 3), src[verts].ravel())), shape=(q.shape[0], nPts)).tocsr()
    return S, ~inside.reshape(nProj, nVirtual)


def _padHighSide(y, x, u, pitch, deltaMm, pad):
    """Weighting, conjugate completion and padding for the missing side being the high-index side of the detector.

    y : (N, M, K) line integrals, x : (K, 6) float64, u : (K, 3) float64 (pixel step of the detector row direction), pitch (mm).
    Returns the padded and weighted data (N + P, M, K) float32 in Fortran order, P, and a dictionary of the details. The centre of the padded
    detector is to be shifted by P / 2 * u (see _shiftDetector)."""
    N, M, K = y.shape
    T = _centralRay(x, u, N)
    if np.any(T <= 0.) or np.any(T >= N):
        raise ValueError('The central ray is outside the detector for some projections')
    TT = N - T
    delta = deltaMm / pitch
    if delta < TT.max():
        delta = float(TT.max())
        raised = True
    else:
        raised = False
    deltaMax = float(T.min()) - _EDGE_MARGIN
    if delta > deltaMax:
        if deltaMax < TT.max():
            raise ValueError('The detector is too close to being centred for the offset FDK preprocessing')
        delta = deltaMax
        clamped = True
    else:
        clamped = False
    nVirtual = int(np.ceil(delta - TT.min())) + 2
    P = min(N, nVirtual + _ZERO_TAIL) if pad is None else int(pad)
    if P < 1:
        raise ValueError('The padding must be at least one row')
    if nVirtual > P:
        nVirtual = P
    nWindow = int(min(N, np.ceil(2. * delta) + _WINDOW_MARGIN))
    S, outside = _interpolationMatrix(x, u, N, nVirtual, nWindow, pitch, K)
    w = 2. * _weightProfile(N + nVirtual, T, delta)
    ypad = np.zeros((N + P, M, K), dtype=np.float32, order='F')
    ypad[:N] = y * w[:N, None, :].astype(np.float32)
    # The hold value of the queries outside the hull: the mean of the last reliable physical rows
    h0 = max(0, N - _HOLD_FIRST)
    h1 = max(h0 + 1, N - _HOLD_LAST + 1)
    wv = w[N:N + nVirtual, None, :].astype(np.float32)
    chunk = 64
    for j0 in range(0, M, chunk):
        j1 = min(M, j0 + chunk)
        vals = np.ascontiguousarray(np.transpose(y[N - nWindow:N, j0:j1, :], (2, 0, 1))).reshape(K * nWindow, j1 - j0)
        out = (S @ vals).reshape(K, nVirtual, j1 - j0)
        hold = np.asarray(y[h0:h1, j0:j1, :], dtype=np.float32).mean(axis=0).T
        out = np.where(outside[:, :, None], hold[:, None, :], out)
        ypad[N:N + nVirtual, j0:j1, :] = np.transpose(out, (1, 2, 0)) * wv
    info = dict(P=P, nVirtual=nVirtual, delta=delta, deltaMm=delta * pitch, raised=raised, clamped=clamped, nWindow=nWindow,
                outsideFraction=float(outside.mean()), T=T)
    return ypad, P, info


def _shiftDetector(x, u, shift):
    """Moves the detector centre of the (K, 6) coordinates x by shift * u (shift in pixels, u (K, 3)). The result has the dtype and the memory order of x."""
    xn = np.array(x, copy=True, order='K')
    view = xn.reshape(-1, 6)
    if not np.may_share_memory(view, xn):
        raise RuntimeError('The detector coordinates cannot be shifted in place')
    view[:, 3:6] = (view[:, 3:6].astype(np.float64) + shift * u).astype(xn.dtype)
    return xn


def _offsetFDKCore(y, x, u, pitch, deltaMm=30., pad=None):
    """Pure array version of the offset FDK preprocessing (no options object).

    y : (N, M, K) line integrals (detector rows fastest), x : (K, 6) [sx sy sz dx dy dz] per projection, u : (K, 2 or 3) the detector row direction
    vector scaled by the pixel size (the options.z of the projector), pitch : the pixel size (mm), deltaMm : the half-width of the transition
    region on the detector (mm), pad : the number of padded rows (None = automatic).

    Returns ypad (N + P, M, K) float32 Fortran-ordered, the new detector coordinates (the dtype and memory order of x) and a dictionary of details."""
    y = np.asarray(y)
    if y.ndim != 3:
        raise ValueError('The projection data must be a three-dimensional array')
    N, M, K = y.shape
    x64 = np.asarray(x, dtype=np.float64).reshape(-1, 6)
    if x64.shape[0] != K:
        raise ValueError('The number of the source/detector coordinate rows does not match the number of projections')
    u = np.asarray(u, dtype=np.float64).reshape(K, -1)
    u3 = np.zeros((K, 3))
    u3[:, :min(3, u.shape[1])] = u[:, :min(3, u.shape[1])]
    T = _centralRay(x64, u3, N)
    high = bool(np.mean(T) > N / 2.)
    if high:
        ypad, P, info = _padHighSide(y, x64, u3, pitch, deltaMm, pad)
        xNew = _shiftDetector(x, u3, P / 2.)
    else:
        # The low-index side is missing: mirror the detector (reverse the rows and the row direction), process it as the high-index
        # side and mirror back. The physical rows then start at the row P of the padded detector.
        yf = y[::-1]
        uf = -u3
        ypadF, P, info = _padHighSide(yf, x64, uf, pitch, deltaMm, pad)
        ypad = np.asfortranarray(ypadF[::-1])
        xNew = _shiftDetector(x, uf, P / 2.)
        info['T'] = N - info['T']
    info['high'] = high
    info['N'] = N
    return ypad, xNew, info


def offsetFDKPreprocess(options):
    """Prepares the offset-detector CT data of an FDK reconstruction (see the module documentation).

    Applies only when options.CT, options.FDK and options.offsetCorrection are set. The measurement data options.SinM (line integrals, except
    for largeDim) are weighted, completed and padded, the detector centre options.x is shifted, options.nRowsD and the other quantities that
    depend on it are updated and options.offsetCorrection is switched off so that the in-kernel weights are not applied on top of these.
    Call this after the projector has been set up (options.addProjector()) and the data have been linearised, before anything else uses nRowsD.

    Returns a state dictionary for offsetFDKRestore, or None if the preprocessing was not applied (the original offset correction is then in effect)."""
    if not (options.CT and options.FDK and options.offsetCorrection) or options.listmode or options.useHelical:
        return None
    N = int(options.nRowsD)
    M = int(options.nColsD)
    K = int(options.nProjections)
    pitch = float(options.dPitchX)
    if int(options.Nt) > 1:
        print('Warning: Offset correction (FDK): the redundancy weighting with conjugate completion is not supported with dynamic data, using the in-kernel weights.')
        return None
    x = np.asarray(options.x)
    z = np.asarray(options.z)
    if x.size != 6 * K or z.size not in (2 * K, 6 * K):
        print('Warning: Offset correction (FDK): unsupported detector coordinates, using the in-kernel weights.')
        return None
    SinM = options.SinM
    if isinstance(SinM, list):
        if len(SinM) != 1:
            print('Warning: Offset correction (FDK): unsupported measurement data layout, using the in-kernel weights.')
            return None
        SinM = SinM[0]
    if SinM.size != N * M * K:
        print('Warning: Offset correction (FDK): the measurement data size does not match the detector, using the in-kernel weights.')
        return None
    if np.size(options.corrVector) > 1 or np.size(options.normalization) > 1:
        print('Warning: Offset correction (FDK): corrections during the reconstruction are not supported with the redundancy weighting '
              'with conjugate completion, using the in-kernel weights.')
        return None
    zz = z.reshape(K, -1)
    u = np.zeros((K, 3))
    nu = 3 if zz.shape[1] >= 6 else 2
    u[:, :nu] = zz[:, :nu]
    y = np.reshape(SinM, (N, M, K), order='F')
    largeDim = bool(options.largeDim)
    # The largeDim FDK linearises raw (not yet linearised) data in the library (log(flat / measurement)), so such data are linearised here for
    # the preprocessing and converted back to the measurement domain below. Already linearised data (usingLinearizedData) stay as they are
    if largeDim and not options.usingLinearizedData:
        y = np.log(np.float32(options.flat) / y.astype(np.float32))
    pad = getattr(options, 'offsetFDKPad', None)
    try:
        ypad, xNew, info = _offsetFDKCore(y, x, u, pitch, float(options.offsetFDKDelta), pad)
    except ValueError as e:
        print('Warning: Offset correction (FDK): %s, using the in-kernel weights.' % str(e))
        return None
    if info['raised']:
        print('Warning: options.offsetFDKDelta (%g mm) is smaller than the physical overlap, using %g mm.' % (float(options.offsetFDKDelta), info['deltaMm']))
    if info['clamped']:
        print('Warning: options.offsetFDKDelta is too large for the detector, using %g mm.' % info['deltaMm'])
    # Back to the measurement domain for the library, unless the input was already linearised (the library then does not linearise it either)
    if largeDim and not options.usingLinearizedData:
        ypad = (np.float32(options.flat) * np.exp(-ypad)).astype(np.float32)
        ypad = np.asfortranarray(ypad)
    P = info['P']
    state = {name: getattr(options, name) for name in _SAVED_ATTRIBUTES if hasattr(options, name)}
    options.SinM = ypad
    newRows = N + P
    options.nRowsD = newRows
    options.Ndist = newRows
    options.det_per_ring = newRows * K
    options.det_w_pseudo = options.det_per_ring
    options.cryst_per_block = M * newRows
    if hasattr(options, 'detectors') and options.rings > 0:
        options.detectors = options.det_w_pseudo * options.rings
    kerroin = M * newRows
    options.nTotMeas = options.nMeas * kerroin
    options.nMeasSubset = np.zeros_like(options.nMeasSubset)
    options.nMeasSubset[:, :] = options.nMeasPerFrameSubset * kerroin
    if options.projector_type in (5, 15, 45, 54, 51):
        options.dSizeXBP = (newRows + 1) * options.dPitchY
    options.x = xNew
    options.offsetCorrection = False
    options.OffsetLimit = np.empty(0, dtype=np.float32)
    if options.verbose >= 1:
        print('Offset correction (FDK): redundancy weights with conjugate completion (half-width %.1f mm), detector padded by %d rows'
              % (info['deltaMm'], P))
    return state


def offsetFDKRestore(options, state):
    """Puts back the options that offsetFDKPreprocess changed (state is its return value)."""
    if state is None:
        return
    for name, value in state.items():
        setattr(options, name, value)
