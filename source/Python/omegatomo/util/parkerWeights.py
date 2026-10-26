# -*- coding: utf-8 -*-
import numpy as np
import warnings

def _S(x):
    """Smooth S-function transition."""
    y = np.zeros_like(x)
    m1 = (x > -0.5) & (x < 0.5)
    m2 = (x >= 0.5)
    y[m1] = 0.5 * (1 + np.sin(np.pi * x[m1]))
    y[m2] = 1
    return y

def _projectionGeometry(options, nProj):
    """
    In-plane (x/y) geometry of every projection as set by setCTCoordinates: the source position, the detector centre and the
    unit vector of the detector row (u) direction. Returns None if the coordinates are not available in the expected form
    (for example helical CT), in which case only the projection angles can be used.
    """
    try:
        x = np.asarray(options.x, dtype=np.float64)
        z = np.asarray(options.z, dtype=np.float64)
    except (AttributeError, TypeError, ValueError):
        return None
    if x.ndim != 2 or z.ndim != 2 or x.shape[1] != 6 or x.shape[0] < nProj or z.shape[0] < nProj or z.shape[1] < 2:
        return None
    sx, sy, dx, dy = x[:nProj, 0], x[:nProj, 1], x[:nProj, 3], x[:nProj, 4]
    un = np.hypot(z[:nProj, 0], z[:nProj, 1])
    if np.any(un <= 0) or np.any(np.hypot(sx, sy) <= 0):
        return None
    return sx, sy, dx, dy, z[:nProj, 0] / un, z[:nProj, 1] / un

def ParkerWeights(options):
    """
    Computes Parker weights as defined in DOI: 10.1118/1.1450132

    The weights are computed from the actual source and detector positions of every projection, i.e. the source polar angle,
    the rotation direction and the fan angle of each detector pixel (including detector offsets, mirrored detector axes and
    arbitrary start angles/arcs that cross 0 or 2*pi) are derived from options.x/options.z (set by setCTCoordinates).
    If the coordinates are not available, the built-in geometry is assumed and the projection angles are used.

    The angular quadrature is also included in the weights: the FDK backprojection integrates with the constant pi/nProjections
    (the full scan step 2*pi/nProjections times the 1/2 of the 360 degree redundancy), and the Parker weights of two conjugate
    rays sum to 1 instead of the 1/2 + 1/2 of a full scan. Each weight is therefore multiplied by dBeta_k * nProjections / pi,
    where dBeta_k is the actual angular step of projection k. For a full 2*pi scan dBeta_k = 2*pi/nProjections and the factor is 2.

    Parameters:
    options (struct): Struct containing all scan parameters and projection data.

    Returns:
    options: The modified options struct where options.SinM has been weighted.
    """
    # Required inputs
    betaIn = np.asarray(options.angles, dtype=np.float64).flatten()
    DSD = options.sourceToDetector
    nU = int(options.nRowsD)
    du = options.dPitchX
    nProj = len(betaIn)

    # Optional parameters with defaults
    if hasattr(options, 'ParkerWeight') and options.ParkerWeight > 0:
        q = options.ParkerWeight
    else:
        q = 0.25

    # Input validation
    if not (np.isscalar(DSD) and np.isscalar(nU) and np.isscalar(du)):
        raise ValueError('sourceToDetector, nRowsD, and dPitchX must be scalars.')

    if q <= 0 or q > 1:
        raise ValueError('options.ParkerWeight must satisfy 0 < q <= 1.')

    if nProj < 2:
        raise ValueError('options.angles must contain at least two projection angles.')

    geom = _projectionGeometry(options, nProj)
    if geom is not None:
        # Source polar angle (unwrapped), rotation direction and fan angle of every pixel from the actual geometry
        sx, sy, dx, dy, ux, uy = geom
        beta = np.unwrap(np.arctan2(sy, sx))
        ccw = beta[-1] >= beta[0]
        sad = np.hypot(sx, sy)
        # Unit vector from the source towards the axis of rotation and its counter-clockwise normal
        cx, cy = -sx / sad, -sy / sad
        j = np.arange(nU) - (nU - 1) / 2
        rx = dx[:, None] + j[None, :] * du * ux[:, None] - sx[:, None]
        ry = dy[:, None] + j[None, :] * du * uy[:, None] - sy[:, None]
        # Fan angle, positive counter-clockwise from the central ray
        alpha = np.arctan2(-rx * cy[:, None] + ry * cx[:, None], rx * cx[:, None] + ry * cy[:, None])
        if not ccw:
            # Clockwise scan: mirror both the source angle and the fan angle
            alpha = -alpha
    else:
        # Built-in geometry from the angles. Detector offset: any sign, scalar or one value per projection
        beta = np.unwrap(betaIn)
        ccw = beta[-1] >= beta[0]
        detOffset = 0.
        if hasattr(options, 'detOffsetRow') and np.size(options.detOffsetRow) > 0:
            off = np.asarray(options.detOffsetRow, dtype=np.float64).reshape(-1)
            if np.any(off != 0):
                detOffset = off[0] if off.size == 1 else np.resize(off, nProj).reshape(-1, 1)
        # Flat-panel detector coordinate u (the built-in detector u axis points in the direction of the counter-clockwise fan angle)
        u = (np.arange(nU) - (nU - 1) / 2) * du + detOffset
        alpha = np.arctan2(np.broadcast_to(u, (nProj, nU)), DSD)
        if not ccw:
            alpha = -alpha

    # Use actual maximum fan half-angle from detector edges
    delta = np.max(np.abs(alpha))

    # Relative source angle along the scan
    betaRel = np.abs(beta - beta[0])
    scanRange = betaRel[-1]

    # Short-scan validity
    minShortScan = np.pi + 2 * delta
    if scanRange + 1e-8 < minShortScan:
        raise ValueError(f'Angular range is too short. Need at least pi + 2*delta = {minShortScan:.6f} rad, got {scanRange:.6f} rad.')

    if scanRange > 2 * np.pi + 1e-8:
        warnings.warn('Angular range exceeds 2*pi. Weights are intended for short/over-scan up to 2*pi.')

    epsilon = max(scanRange - (np.pi + 2 * delta), 0)

    # Angular quadrature, see above
    dBeta = np.abs(np.gradient(beta))
    angFactor = dBeta * nProj / np.pi

    # Compute weights. If options.SinM can be the caller's array (see _privateOptionsCopy())
    # the weighted data is written to a new array, otherwise in place.
    from omegatomo.reconstruction.recomain import _sharesCallerMemory
    sino = options.SinM
    flat = sino.ndim != 3
    if flat:
        sino = np.reshape(sino, (nU, -1, nProj), order='F')
    weightedSinM = np.empty_like(sino) if _sharesCallerMemory(options, sino) else sino
    for iu in range(nU):
        g = alpha[:, iu]

        B = epsilon + 2 * delta - 2 * g
        b = q * B
        B2 = epsilon + 2 * delta + 2 * g
        b2 = q * B2

        if np.any(b <= 0) or np.any(b2 <= 0):
            raise ValueError('Encountered nonpositive transition length b. Check geometry.')

        x1 = betaRel / b - 0.5
        x2 = (betaRel - B) / b + 0.5
        x3 = (betaRel - np.pi + 2 * g) / b2 - 0.5
        x4 = (betaRel - np.pi - 2 * delta - epsilon) / b2 + 0.5

        w2 = 0.5 * (_S(x1) + _S(x2) - _S(x3) - _S(x4))
        w2 = np.clip(w2, 0, 1) * angFactor
        weightedSinM[iu, :, :] = sino[iu, :, :] * w2.reshape(1, -1)
    if flat:
        weightedSinM = np.reshape(weightedSinM, options.SinM.shape, order='F')
    options.SinM = weightedSinM

    # return options
