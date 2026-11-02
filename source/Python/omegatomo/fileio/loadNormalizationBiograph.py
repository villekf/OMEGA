# -*- coding: utf-8 -*-
"""
Loader of the normalization component files of the Siemens Biograph Vision.
"""

import os
import re
import warnings

import numpy as np

_NR = 520  # Radial bins of the scanner
_NS = 798  # Crystals per ring, including the gaps
_R = 80  # Rings
_MCT_MAX_FLOATS = 1255709  # Largest mCT normalization file (in floats)
_VISION_DEFAULT_OFFSETS = (0, 13312000, 13355680, 13611040)
_VISION_MIN_BYTES = 13636640


def _readNormHeader(fpath):
    """
    Reads the component offsets and sizes from the Interfile header that is embedded in a
    .ptd file (after the data) or stored in the sidecar <file>.hdr. Returns a list of
    (offset in bytes, dimensions) for each component or an empty list if there is no header.
    """
    fsz = os.path.getsize(fpath)
    nTail = min(fsz, 4 * 2**20)
    with open(fpath, 'rb') as f:
        f.seek(fsz - nTail)
        buf = f.read(nTail)
        pos = buf.find(b'!INTERFILE')
        if pos < 0 and fsz > nTail:
            f.seek(0)
            buf = f.read()
            pos = buf.find(b'!INTERFILE')
    txt = ''
    if pos >= 0:
        txt = buf[pos:].decode('latin-1')
        # The header ends where the binary (DICOM) data continues
        z = txt.find('\x00')
        if z >= 0:
            txt = txt[:z]
    else:
        base = os.path.splitext(fpath)[0]
        for hdr in (fpath + '.hdr', base + '.hdr'):
            if os.path.isfile(hdr):
                with open(hdr, 'r', encoding='latin-1') as f:
                    txt = f.read()
                break
    if len(txt) == 0:
        return []
    offsets = re.findall(r'data offset in bytes\[(\d+)\]\s*:=\s*(\d+)', txt)
    sizes = re.findall(r'matrix size\[(\d+)\]\s*:=\s*\{([^}]*)\}', txt)
    if len(offsets) == 0 or len(offsets) != len(sizes):
        warnings.warn('The normalization header could not be parsed, using the default layout')
        return []
    off = {int(k): int(v) for k, v in offsets}
    dim = {int(k): tuple(int(x) for x in v.replace(',', ' ').split()) for k, v in sizes}
    if set(off) != set(dim):
        warnings.warn('The normalization header could not be parsed, using the default layout')
        return []
    return [(off[k], dim[k]) for k in sorted(off)]


def _readComponent(f, offset, dims):
    """Reads one float32 component (the first dimension varies fastest)"""
    n = int(np.prod(dims))
    f.seek(offset)
    x = np.fromfile(f, dtype='<f4', count=n)
    if x.size != n:
        raise ValueError('The normalization file ended before the end of the normalization component')
    return x.astype(np.float64).reshape(dims, order='F')


def _siemensBin(s1, s2):
    """
    Siemens (CTI) sinogram bin (t, v) of the crystal pairs (s1, s2) of a ring of 798 crystals,
    oriented such that s1 = (v + t // 2) % 798 and s2 = (v - (t + 1) // 2 + 399) % 798. The
    pair is swapped if needed. Returns t, v and the (possibly swapped) crystals s1 and s2.
    """
    s1 = np.asarray(s1, dtype=np.int64)
    s2 = np.asarray(s2, dtype=np.int64)
    t = (s1 - s2 + 399) % 798
    t = np.where(t >= 399, t - 798, t)
    v = ((s1 + s2 - 399 + (t % 2)) % 798) // 2
    bad = (v + t // 2) % 798 != s1
    if np.any(bad):
        s1o = np.where(bad, s2, s1)
        s2o = np.where(bad, s1, s2)
        t = (s1o - s2o + 399) % 798
        t = np.where(t >= 399, t - 798, t)
        v = ((s1o + s2o - 399 + (t % 2)) % 798) // 2
        s1, s2 = s1o, s2o
    return t, v, s1, s2


def _sinogramPlane(rn1, rn2, span, rings, seg):
    """Zero-based sinogram plane of the ring numbers rn1 and rn2 (the same rule as in saveSinogram)"""
    rn1 = np.asarray(rn1, dtype=np.int64)
    rn2 = np.asarray(rn2, dtype=np.int64)
    if span <= 1:
        return rn2 * rings + rn1
    h = span // 2
    erotus = rn1 - rn2
    summa = rn1 + rn2
    far = np.abs(erotus) > h
    k = (np.abs(erotus) + h) // span
    out = summa.copy()
    if np.any(far):
        kf = k[far]
        if int(kf.max()) * 2 > seg.size:
            raise ValueError('The segment table has too few segments for the ring difference')
        off = np.where(erotus[far] < 0, seg[(kf - 1) * 2], seg[(kf - 1) * 2 + 1])
        out[far] = summa[far] - (h * (2 * kf - 1) + kf) + off
    return out


def _lorLengthFactor(options, W, x1, x2, n):
    """
    Transaxial LOR length factor D^(-n) of the detector pairs (x1, x2) (zero-based), where D is
    the distance between the two detectors, from the detector coordinates of OMEGA
    (detectorCoordinates). Pseudo detectors have finite coordinates; non-finite factors are
    set to zero.
    """
    from omegatomo.projector.detcoord import detectorCoordinates
    for name in ('cr_p', 'diameter', 'blocks_per_ring', 'cryst_per_block'):
        v = np.asarray(getattr(options, name, None), dtype=np.float64)
        if v.size == 0 or np.any(~(v > 0)):
            raise ValueError('The Vision LOR length factor needs the detector geometry, but options.' + name +
                             ' is missing or not positive. Set the detector geometry or set '
                             'options.visionLORLengthExponent = 0 to disable the factor.')
    x, y = detectorCoordinates(options)
    x = np.asarray(x, dtype=np.float64).ravel()
    y = np.asarray(y, dtype=np.float64).ravel()
    if x.size != W:
        raise ValueError('The number of transaxial detector coordinates (' + str(x.size) + ') differs from det_w_pseudo (' +
                         str(W) + '). Set the detector geometry correctly or set options.visionLORLengthExponent = 0 '
                         'to disable the LOR length factor.')
    with np.errstate(divide='ignore', invalid='ignore'):
        h = np.hypot(x[x1] - x[x2], y[x1] - y[x2]) ** (-n)
    h[~np.isfinite(h)] = 0.
    return h


def loadNormalizationBiograph(fpath, options, returnEfficiency=False, returnScale=False):
    """
    Loads the normalization component file (.n or the .ptd calibration file) of the Siemens
    Biograph Vision and converts it into sinogram normalization coefficients with the
    sinogram dimensions of OMEGA.

    The component layout is read from the Interfile header that follows the data in a .ptd
    file or, for a plain .n file, from the sidecar header <file>.hdr. Without a header the
    standard Vision layout is assumed.

    The efficiency of every LOR is
        eff = E(s1, ring1) * E(s2, ring2) * I(p, r) / (G(r, P) * A(P)),
    where E are the crystal efficiencies, I the crystal interference, G the geometric effects
    and A the axial effects of the span-1 plane P of the scanner. The efficiencies of all LORs
    that are binned into the same sinogram bin by OMEGA (the same binning as in saveSinogram,
    including the mashing, span and the sinogram size) are summed, i.e. the multiplicity of the
    bin is included, and the sum is normalized such that the mean efficiency of the LORs is one.
    The component 5 of the file (dead-time parameters) depends on the count rate and is
    ignored; a global dead-time factor can be given with options.global_correction_factor.

    The efficiency of every LOR is also multiplied by D^(-n), where D is the transaxial (x-y)
    distance between the two detectors of the LOR and n = options.visionLORLengthExponent
    (default 2, 0 = no length factor). The factor is empirical: it was validated against the
    Siemens reconstruction of measured Vision data, where the radial RMS error decreased from
    4.66 to 2.15 % with n = 2. The detector positions are those of
    omegatomo.projector.detcoord.detectorCoordinates, so the detector geometry of the options
    (diameter, DOI, cr_p, blocks_per_ring, cryst_per_block, transaxialGaps, firstBlockAngle,
    flip_image and offangle) has to be set.

    Requirements: 80 rings, 760 detectors per ring without (det_w_pseudo = 760) or 798 with
    (det_w_pseudo = 798, one pseudo detector at the end of each block of 20 crystals) pseudo
    detectors. Ndist, Nang, span, ring_difference, segment_table and ndist_side can be freely
    selected.

    LORs whose radial distance is outside of the 520 radial bins of the scanner (possible
    with the 760 crystal sinogram and a large Ndist) are given the crystal interference one
    and the geometric effect of the nearest edge bin. A warning with the number of such LORs
    and their share of the total efficiency is shown.

    Parameters
    ----------
    fpath : str
        Path to the normalization file.
    options : class object
        OMEGA class object with the scanner properties (Ndist, Nang, rings, det_per_ring,
        det_w_pseudo, cryst_per_block, span, ring_difference, segment_table, ndist_side and
        visionLORLengthExponent).
    returnEfficiency : bool, optional
        If True, also the summed efficiency of each sinogram bin is returned.
    returnScale : bool, optional
        If True, also the scale that was divided out when the mean efficiency of the LORs was
        normalized to one is returned (last): the mean over the LORs of
        E1 * E2 * I / (G * A) * (Dref / D)^n, with Dref = options.diameter + 2 * options.DOI (the
        effective detector ring diameter in mm) and n = options.visionLORLengthExponent. With
        n = 0 it is the mean of the plain efficiencies. It is used by quantificationFactorBiograph
        to convert the images to Bq/mL.

    Returns
    -------
    normalization : ndarray, float32, Fortran-order [Ndist, Nang, TotSinos]
        The normalization coefficients. They MULTIPLY the measured sinogram
        (normalized = measured * normalization), and 1 / normalization is the efficiency
        (options.normalization during the reconstruction). Bins without any efficiency (e.g.
        only pseudo detector LORs) are zero and their efficiency should be taken as zero,
        not infinity.
    efficiency : ndarray, float32, Fortran-order [Ndist, Nang, TotSinos]
        Only if returnEfficiency is True. The summed efficiency of each sinogram bin, i.e.
        normalization = 1 / efficiency (zero where the efficiency is zero).
    meanEfficiency : float
        Only if returnScale is True (returned last). The scale described above.

    Raises
    ------
    NotImplementedError
        If the file is not a Vision normalization file.
    """
    comp = _readNormHeader(fpath)
    if len(comp) > 0:
        isVision = len(comp) >= 4 and tuple(comp[2][1]) == (_NS, _R)
    else:
        isVision = os.path.getsize(fpath) / 4 > _MCT_MAX_FLOATS
    if not isVision:
        raise NotImplementedError('Only the Siemens Biograph Vision normalization files are supported in Python')

    Nd = int(options.Ndist)
    Nang = int(options.Nang)
    rings = int(options.rings)
    W = int(options.det_w_pseudo)
    span = int(options.span)
    rd = int(options.ring_difference)
    nSide = getattr(options, 'ndist_side', 1)
    nSide = 1 if nSide is None or np.asarray(nSide).size == 0 else int(np.asarray(nSide).ravel()[0])
    nLen = getattr(options, 'visionLORLengthExponent', 2.)
    nLen = 2. if nLen is None or np.asarray(nLen).size == 0 else nLen
    if np.asarray(nLen).size != 1 or not np.isreal(nLen) or not np.isfinite(float(np.asarray(nLen).ravel()[0])):
        raise ValueError('options.visionLORLengthExponent has to be a finite real scalar')
    nLen = float(np.asarray(nLen).ravel()[0])
    NR, NS, R = _NR, _NS, _R
    if rings != R:
        raise ValueError('The Vision normalization requires 80 crystal rings')
    if W == int(options.det_per_ring):
        if W != 760 or int(options.cryst_per_block) != 20:
            raise ValueError('The Vision normalization requires 760 detectors per ring (38 blocks of 20 crystals)')
    elif W != NS:
        raise ValueError('The Vision normalization requires det_w_pseudo of 760 (no pseudo detectors) or 798 (one pseudo detector per block)')
    if span <= 1:
        totSinos = rings * rings
        seg = np.empty(0, dtype=np.int64)
    else:
        seg = np.cumsum(np.asarray(options.segment_table, dtype=np.int64).ravel())
        totSinos = int(seg[-1])
    if rd > rings - 1:
        raise ValueError('Ring difference cannot be larger than the number of rings - 1')

    # Read the components
    if len(comp) == 0:
        if os.path.getsize(fpath) < _VISION_MIN_BYTES:
            raise ValueError('Unrecognized normalization file layout')
        off = _VISION_DEFAULT_OFFSETS
    else:
        dimsC = ((NR, 6400), (21, NR), (NS, R), (6400,))
        for kk in range(4):
            dims = tuple(comp[kk][1])
            if int(np.prod(dims)) != int(np.prod(dimsC[kk])) or (kk < 3 and dims != dimsC[kk]):
                raise ValueError('Unexpected size of the Vision normalization component ' + str(kk + 1))
        off = [comp[kk][0] for kk in range(4)]
    with open(fpath, 'rb') as f:
        G = _readComponent(f, off[0], (NR, 6400))
        I = _readComponent(f, off[1], (21, NR))
        E = _readComponent(f, off[2], (NS, R))
        A = _readComponent(f, off[3], (6400,))

    # Siemens span-1 planes: blocks of [80, 79, 79, 78, 78, ..., 1, 1] planes
    sizes = np.concatenate(([R], np.repeat(np.arange(R - 1, 0, -1), 2)))
    starts = np.concatenate(([0], np.cumsum(sizes[:-1])))
    RA, RB = np.meshgrid(np.arange(R), np.arange(R), indexing='ij')
    dd = np.abs(RA - RB)
    blk = np.maximum(2 * dd - 1, 0)
    Pmat = starts[blk] + np.minimum(RA, RB)
    # The blocks 2d-1 and 2d should be identical
    dmax = 0.
    for d in range(1, R):
        n = R - d
        ia = starts[2 * d - 1] + np.arange(n)
        ib = starts[2 * d] + np.arange(n)
        dmax = max(dmax, np.max(np.abs(A[ia] - A[ib])) / np.max(np.abs(A[ia])),
                   np.max(np.abs(G[:, ia] - G[:, ib])) / np.max(np.abs(G[:, ia])))
    if dmax > 1e-6:
        warnings.warn('The axial or geometric effects of the plane blocks 2d-1 and 2d differ by up to ' + str(dmax) + ' (relative), only the block 2d-1 is used')

    # Transaxial LOR table: all detector pairs x1 < x2, binned as in saveSinogram
    x1, x2 = np.triu_indices(W, 1)
    x1 = x1.astype(np.int64)
    x2 = x2.astype(np.int64)
    half = W // 2
    j = ((x2 + x1 + half + 1) % W) // 2
    i = np.abs(x2 - x1 - half)
    i = np.where((x1 < j) | (j + half < x2), -i, i)
    swap = ((j * 2) < -i) | (i <= (j - half) * 2)
    if Nd % 2 == 0:
        if nSide == 0:
            nSide = 1
        acc = (i <= (Nd // 2 + min(0, nSide))) & (i >= (-(Nd // 2) + max(0, nSide)))
        ib = i + Nd // 2 - max(0, nSide)
    else:
        acc = (i <= Nd // 2) & (i >= -(Nd // 2))
        ib = i + Nd // 2
    jb = j // (half // Nang)
    if np.any(acc & (jb >= Nang)):
        raise ValueError('Nang is not compatible with the number of detectors per ring')
    # Siemens crystals of the OMEGA detectors
    x = np.arange(W, dtype=np.int64)
    if W == NS:
        sx = (x + 1) % NS
    else:
        sx = x + x // int(options.cryst_per_block) + 1
    s1 = sx[x1]
    s2 = sx[x2]
    # Transaxial LOR length factor D^(-n)
    if nLen != 0:
        hLen = _lorLengthFactor(options, W, x1, x2, nLen)
    else:
        hLen = np.ones(x1.shape, dtype=np.float64)
    # Pseudo detectors (gap crystals) have no efficiency
    acc = acc & (s1 % 21 != 0) & (s2 % 21 != 0)
    swap = swap[acc]
    s1 = s1[acc]
    s2 = s2[acc]
    hLen = hLen[acc]
    NdNa = Nd * Nang
    bin0 = (ib + jb * Nd)[acc]
    # Orientation-corrected Siemens bin
    t, v, _, _ = _siemensBin(s1, s2)
    r = t + 260
    p = v % 21
    oor = (r < 0) | (r > NR - 1)
    rc = np.clip(r, 0, NR - 1)
    Ifac = np.ones(r.shape, dtype=np.float64)
    Ifac[~oor] = I[p[~oor], r[~oor]]
    # LORs with the swap flag first
    ord_ = np.argsort(~swap, kind='stable')
    s1 = s1[ord_]
    s2 = s2[ord_]
    bin0 = bin0[ord_]
    rc = rc[ord_]
    Ifac = Ifac[ord_] * hLen[ord_]
    oor = oor[ord_]
    nS = int(np.sum(swap))
    nTr = bin0.size
    oorIdx = np.flatnonzero(oor)
    bS = bin0[:nS]
    bN = bin0[nS:]

    E1 = np.ascontiguousarray(E[s1, :].T)
    E2 = np.ascontiguousarray(E[s2, :].T)
    with np.errstate(divide='ignore', invalid='ignore'):
        Winv = 1. / (G * A[np.newaxis, :])
    Winv[~np.isfinite(Winv)] = 0.
    Winv = np.ascontiguousarray(Winv.T)  # [plane, radial bin]

    # Sinogram planes of the ring pairs: the detector 1 is at x1 and the detector 2 at x2
    # (ring numbers are swapped for the swapped LORs and then again since x2 > x1)
    okPair = dd <= rd
    # (the ring pairs beyond the ring difference are not binned and are set to a valid pair)
    planeS = _sinogramPlane(np.where(okPair, RA, 0), np.where(okPair, RB, 0), span, rings, seg)
    planeN = _sinogramPlane(np.where(okPair, RB, 0), np.where(okPair, RA, 0), span, rings, seg)
    if np.any(okPair & ((planeS < 0) | (planeS >= totSinos) | (planeN < 0) | (planeN >= totSinos))):
        raise ValueError('The segment table does not match the number of rings and the ring difference')
    nRP = int(np.sum(okPair))

    S = np.zeros((totSinos, NdNa), dtype=np.float32)
    oorSum = 0.
    for ra in range(R):
        E1a = E1[ra]
        for rb in range(R):
            if not okPair[ra, rb]:
                continue
            w = Winv[Pmat[ra, rb]]
            e = E1a * E2[rb] * Ifac * w[rc]
            if nS > 0:
                S[planeS[ra, rb]] += np.bincount(bS, weights=e[:nS], minlength=NdNa).astype(np.float32)
            if nS < nTr:
                S[planeN[ra, rb]] += np.bincount(bN, weights=e[nS:], minlength=NdNa).astype(np.float32)
            if oorIdx.size > 0:
                oorSum += float(e[oorIdx].sum())

    # Mean efficiency of the LORs is one
    totEff = float(S.sum(dtype=np.float64))
    if totEff <= 0:
        raise ValueError('The Vision normalization has no non-zero efficiencies')
    eff = S / np.float32(totEff / (nTr * nRP))
    # Mean efficiency with the length factor referred to the effective ring diameter,
    # (Dref / D)^n instead of D^(-n) (a missing DOI is zero)
    meanEff = totEff / (nTr * nRP)
    if nLen != 0:
        DOI = getattr(options, 'DOI', None)
        DOI = 0. if DOI is None or np.asarray(DOI).size == 0 else float(np.asarray(DOI).ravel()[0])
        meanEff *= (float(np.asarray(options.diameter).ravel()[0]) + 2. * DOI) ** nLen
    if oorIdx.size > 0:
        warnings.warn('%d LORs (%.3g %% of all LORs) are outside of the %d radial bins of the Vision normalization; '
                      'their crystal interference is set to one and the geometric effect of the edge bin is used. '
                      'They carry %.3g %% of the total efficiency.' % (oorIdx.size * nRP, 100. * oorIdx.size / nTr, NR, 100. * oorSum / totEff))
    coef = np.zeros_like(eff)
    nz = eff > 0
    coef[nz] = np.float32(1.) / eff[nz]
    coef = np.asfortranarray(coef.reshape((totSinos, Nang, Nd)).transpose(2, 1, 0))
    out = (coef,)
    if returnEfficiency:
        out += (np.asfortranarray(eff.reshape((totSinos, Nang, Nd)).transpose(2, 1, 0)),)
    if returnScale:
        out += (meanEff,)
    return out if len(out) > 1 else coef
