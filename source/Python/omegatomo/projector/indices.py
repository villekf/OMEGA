# -*- coding: utf-8 -*-
"""
Copyright (C) 2024-2026 Ville-Veikko Wettenhovi

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
import numpy as np


def _prime_factor_order(n, tyyppi):
    """
    Computes the prime-factor (mixed-radix) permutation of range(n) used by
    subsetType 11. Equivalent to computing, for every projection index rr,
    a mixed-radix digit expansion of rr using the prime factors of n (from
    smallest to largest place value) and re-reading those digits with the
    place values reversed.

    Parameters
    ----------
    n : int
        Number of projections to order (must be > 1).
    tyyppi : NumPy dtype
        Integer dtype used for the returned index array (np.uint32 for
        non-CT, np.uint64 for CT).

    Returns
    -------
    indices : NumPy array (tyyppi)
        The prime-factor ordering of range(n).
    """
    factors = []
    value = n
    divisor = 2
    while value > 1:
        while value % divisor == 0:
            factors.append(divisor)
            value //= divisor
        divisor += 1
    n2 = np.array(factors)
    nn = np.flipud(np.cumprod(np.flipud(n2[1:])))
    N = np.size(factors)
    p = np.zeros((N, n), dtype=tyyppi)
    for ll in range(1, N + 1):
        p1 = np.arange(factors[ll - 1])
        if ll == 1:
            p1 = np.tile(p1, n // np.size(p1))
        else:
            p1 = np.repeat(p1, np.prod(n2[:ll - 1]))
            p1 = np.tile(p1, n // np.size(p1))
        p[ll - 1, :] = p1
    # Vectorized form of: for rr in range(n): indices[rr] = sum(p[:len(nn), rr] * nn) + p[-1, rr]
    indices = np.zeros(n, dtype=tyyppi)
    indices[:] = np.sum(p[:len(nn), :] * nn[:, None], axis=0) + p[-1, :]
    return indices


def _goldenAngleProjectionOrder(angles, subsets, tyyppi):
    """
    Golden angle ordering of projection images (subsetType 10). Port of
    index_maker.m:279-317. The projection whose angle is closest to the current
    golden angle position is picked repeatedly (and removed from the pool),
    and the picks are dealt to the subsets one after another. The first
    subsets - 1 subsets get floor(n / subsets) projections each and the last
    one ceil(n / subsets) (as in MATLAB, this loses a projection if n is not
    divisible by the number of subsets).

    Parameters
    ----------
    angles : NumPy array
        Projection angles.
    subsets : int
        Number of subsets (> 1).
    tyyppi : NumPy dtype
        Integer dtype used for the returned index arrays.

    Returns
    -------
    index : list of NumPy arrays (tyyppi)
        The (0-based) projection indices of each subset.
    """
    angles = np.abs(np.asarray(angles).ravel())
    if not np.issubdtype(angles.dtype, np.floating):
        angles = angles.astype(np.float64)
    dt = angles.dtype.type
    n = angles.size
    angles = angles - np.min(angles)
    remaining = np.arange(n)
    maksimi = np.max(angles)
    ga = dt(2.39996322972865332) * (maksimi / dt(2. * np.pi))
    angle = dt(0.)
    index = []
    for kk in range(subsets):
        if kk < subsets - 1:
            count = n // subsets
        else:
            count = -(-n // subsets)
        ind = np.zeros(count, dtype=tyyppi)
        for ii in range(count):
            I = np.argmin(np.abs(angles - angle))
            ind[ii] = remaining[I]
            angles = np.delete(angles, I)
            remaining = np.delete(remaining, I)
            angle = angle + ga
            if angle > maksimi:
                angle = angle - maksimi
        index.append(ind)
    return index


def _sinogramLORAngles(options):
    """
    Computes the angle (in degrees, 0 to 180) of the LOR of every sinogram bin
    from the detector coordinates. Common part of subset_angles.m (non-raw
    branch, lines 37-51) and goldenAngleSubsets.m. options.Ndist must already
    include the possible sampling increase (indexMaker does this).

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Returns
    -------
    kulma : NumPy array, shape (Ndist * Nang,)
        The unrounded angle of each sinogram bin in column-major (Ndist
        fastest) order.
    """
    from omegatomo.projector.detcoord import detectorCoordinates, sinogramCoordinates2D
    xp, yp = detectorCoordinates(options)
    x, y = sinogramCoordinates2D(options, xp, yp)
    # Lengths of the triangle sides
    jokux = np.abs(x[:, 0] - x[:, 1])
    jokuy = np.abs(y[:, 0] - y[:, 1])
    # Indices on the other side of the circle (angle over 90 degrees)
    ind = ((x[:, 0] > x[:, 1]) & (y[:, 0] > y[:, 1])) | ((x[:, 1] > x[:, 0]) & (y[:, 1] > y[:, 0]))
    # LORs that are vertical (or horizontal) up to numerical precision are
    # treated as exactly vertical (0 degrees) or horizontal (90 degrees), so
    # that rounding errors (e.g. single vs. double precision) cannot move them
    # to the other side of the circle
    vertical = jokux <= 1e-6 * np.max(jokuy)
    horizontal = jokuy <= 1e-6 * np.max(jokux)
    ind = ind & ~(vertical | horizontal)
    with np.errstate(divide='ignore', invalid='ignore'):
        kulma = np.degrees(np.arctan(jokux / jokuy))
    kulma[horizontal] = 90
    kulma[vertical] = 0
    kulma[ind] = kulma[ind] + 90
    return kulma


def _subsetAngles(options, tyyppi):
    """
    Angle-based subsets for sinogram data (subsetType 6). Port of the non-raw
    branch of subset_angles.m. The sinogram bins are grouped by the angle of
    their LOR, n_angles degrees per subset.

    2D data (only the direct planes are used: no oblique segments or
    NSinos < TotSinos): 180 / n_angles subsets, the subset j contains the bins
    with a (rounded) angle in [j * n_angles, (j + 1) * n_angles) of all
    (direct) planes.

    3D data: 360 / n_angles subsets. The oblique segments of all segment pairs
    belong to the subset that contains their full angle (the angle of the first
    member of a pair, angle + 180 for the second member). The direct planes are
    split by the plane parity: the bins with the angle range j of the even
    planes (0, 2, ...) go to the subset j and those of the odd planes to the
    subset j + 180 / n_angles, i.e. every subset spans the whole axial range.

    Order inside a subset: the direct planes first, then the oblique segments
    in increasing segment order; every block is plane-major (a plane after the
    other), the bins of a plane in ascending sinogram order.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data. Requires
        n_angles, segment_table, NSinos and TotSinos. options.Ndist must
        already include the possible sampling increase.
    tyyppi : NumPy dtype
        Integer dtype used for the returned index array.

    Returns
    -------
    index : NumPy array (tyyppi)
        The (0-based) sinogram indices, subset-major.
    pituus : NumPy array (int64)
        Number of indices in each subset.
    """
    if not hasattr(options, 'n_angles'):
        raise ValueError('Subset type 6 requires options.n_angles (the number of angles combined per subset) to be set!')
    nAngles = options.n_angles
    maara = int(round(180 / nAngles))
    if maara < 1 or abs(180 / nAngles - maara) > 1e-9:
        raise ValueError('options.n_angles has to divide 180 evenly!')
    segmentTable = np.asarray(options.segment_table).ravel().astype(np.int64)
    sinoLength = int(options.Nang) * int(options.Ndist)
    from omegatomo.util.matlabRound import matlabRound
    # The angles are first snapped to 0.001 degrees, so that LORs with an angle that is an exact half
    # degree (up to numerical precision) are rounded identically in single and double precision
    kulma = matlabRound(matlabRound(_sinogramLORAngles(options).astype(np.float64) * 1e3) / 1e3)
    # Combine similar angles together
    kulma[kulma == 180] = 179
    nPairs = segmentTable.size // 2
    use3D = nPairs > 0 and not (int(options.NSinos) < int(options.TotSinos))
    if use3D:
        nSubsets = 2 * maara
        nDirect = int(segmentTable[0])
        # First plane and number of planes of the first member of each segment pair
        firstPlane = [int(np.sum(segmentTable[:(rd + 1) * 2 - 1])) for rd in range(nPairs)]
        nPlanes = [int(segmentTable[1 + rd * 2]) for rd in range(nPairs)]
    else:
        nSubsets = maara
        nDirect = int(options.NSinos)
    index = []
    pituus = []
    for kk in range(nSubsets):
        j = kk % maara
        par = kk // maara
        bins = np.flatnonzero((kulma >= j * nAngles) & (kulma < (j + 1) * nAngles)).astype(np.int64)
        if use3D:
            planes = [np.arange(par, nDirect, 2, dtype=np.int64)]
            for rd in range(nPairs):
                planes.append(firstPlane[rd] + nPlanes[rd] * par + np.arange(nPlanes[rd], dtype=np.int64))
            planes = np.concatenate(planes)
        else:
            planes = np.arange(nDirect, dtype=np.int64)
        ind = (planes[:, None] * sinoLength + bins[None, :]).ravel()
        index.append(ind)
        pituus.append(ind.size)
    index = np.concatenate(index).astype(tyyppi)
    return index, np.asarray(pituus, dtype=np.int64)


def _goldenAngleSubsets(options, tyyppi):
    """
    Golden angle sampling of the sinogram angles (subsetType 7, sinogram data).
    Port of goldenAngleSubsets.m. Whole sinogram columns (all radial bins and
    all sinograms of one angle) are dealt to the subsets in the order of the
    golden angle (180 / golden ratio degrees) steps.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.
        options.Ndist must already include the possible sampling increase.
    tyyppi : NumPy dtype
        Integer dtype used for the returned index array.

    Returns
    -------
    index : NumPy array (tyyppi)
        The (0-based) sinogram indices, subset-major.
    pituus : NumPy array (int64)
        Number of indices in each subset.
    """
    Ndist = int(options.Ndist)
    Nang = int(options.Nang)
    NSinos = int(options.NSinos)
    subsets = int(options.subsets)
    kulma = _sinogramLORAngles(options).reshape((Ndist, Nang), order='F')
    ga = 180. / ((1. + np.sqrt(5.)) / 2.)
    koko = -(-Nang // subsets)
    if Nang % subsets > 0:
        last = Nang - (subsets - 1) * koko
    else:
        last = koko
    # The column where the angle wraps around from 180 to 0 degrees gets the largest angle
    loc = np.abs(np.diff(kulma, axis=0)[0, :])
    loc = np.flatnonzero(loc > np.mean(loc))
    if loc.size > 0:
        kulma[:, loc[0]] = np.nanmax(kulma[:, loc[0]])
    keskiarvo = np.mean(kulma, axis=0, dtype=np.float64).astype(kulma.dtype)
    origKeskiarvo = keskiarvo.copy()
    angle0 = 0.
    used = set()
    indices = []
    for kk in range(subsets):
        if kk < subsets - 1:
            count = koko
        else:
            count = last
        cols = []
        for ll in range(count):
            closestIndex = np.argmin(np.abs(keskiarvo - angle0))
            # Equal averages of different columns are taken in order
            closestIndex2 = [c for c in np.flatnonzero(origKeskiarvo == keskiarvo[closestIndex]) if c not in used][0]
            used.add(closestIndex2)
            cols.append(closestIndex2)
            keskiarvo = np.delete(keskiarvo, closestIndex)
            angle0 = angle0 + ga
            if angle0 >= 180:
                angle0 = angle0 - 180
        indices.append(cols)
    index1 = np.arange(Ndist * Nang * NSinos, dtype=np.int64).reshape((Ndist, Nang, NSinos), order='F')
    index = []
    pituus = np.zeros(subsets, dtype=np.int64)
    for i in range(subsets):
        index3 = [index1[:, c, :].ravel(order='F') for c in indices[i]]
        index3 = np.concatenate(index3) if len(index3) > 0 else np.empty(0, dtype=np.int64)
        index.append(index3)
        pituus[i] = index3.size
    return np.concatenate(index).astype(tyyppi), pituus


def _projection_subset_indices(nProjections, subsets, subsetType, generator, tyyppi, angles=None):
    """
    Splits `nProjections` projections (either the full, static acquisition or
    a single dynamic timeframe) into `subsets` subsets for subsetType 8, 9, 10,
    11, 12, or any other value, which falls back to a single subset
    containing every projection in order (this is the branch actually taken
    whenever subsets == 1, regardless of the requested subsetType).

    Parameters
    ----------
    nProjections : int
        Number of projections to split.
    subsets : int
        Number of subsets.
    subsetType : int
        The subset selection method.
    generator : numpy.random.Generator or None
        Used (and advanced) only for subsetType 9 with subsets > 1.
    tyyppi : NumPy dtype
        Integer dtype used for the returned index array.
    angles : NumPy array or None
        Projection angles. Required (only) for subsetType 10 with subsets > 1.

    Returns
    -------
    index : NumPy array (tyyppi)
        Concatenated (subset-major) index array.
    lengths : NumPy array (int64), shape (subsets,)
        Number of indices in each subset. Only lengths[0] is populated by the
        fallback (non 8/9/10/11) branch; the caller is responsible for any
        CT/PET/SPECT-specific override of that single value.
    """
    sProjections = nProjections // subsets
    modi = np.mod(nProjections, subsets)
    uu = (modi > 0).astype(np.int64)
    ind1 = 0
    ind2 = sProjections + uu
    temp_ind = []
    lengths = np.zeros(subsets, dtype=np.int64)
    if subsetType == 8 and subsets > 1:
        for kk in range(subsets):
            index1 = np.arange(kk, nProjections, subsets).astype(tyyppi)
            temp_ind.append(index1)
            lengths[kk] = np.size(index1)
    elif subsetType == 9 and subsets > 1:
        apu = generator.permutation(nProjections).astype(tyyppi)
        for kk in range(subsets):
            index1 = apu[ind1:ind2]
            temp_ind.append(index1)
            lengths[kk] = np.size(index1)
            modi = modi - 1
            if modi <= 0:
                uu = 0
            ind1 = ind2
            ind2 = ind2 + sProjections + uu
    elif subsetType == 10 and subsets > 1:
        if angles is None or np.size(angles) == 0:
            raise ValueError('Subset type 10 is only available with CT data and when the projections angles (options.angles) are supplied!')
        temp_ind = _goldenAngleProjectionOrder(angles, subsets, tyyppi)
        for kk in range(subsets):
            lengths[kk] = np.size(temp_ind[kk])
    elif subsetType == 11 and subsets > 1:
        indices = _prime_factor_order(nProjections, tyyppi)
        for kk in range(subsets):
            index1 = indices[ind1:ind2]
            temp_ind.append(index1)
            lengths[kk] = np.size(index1)
            modi = modi - 1
            if modi <= 0:
                uu = 0
            ind1 = ind2
            ind2 = ind2 + sProjections + uu
    elif subsetType == 12 and subsets > 1:
        if nProjections % 2 != 0:
            raise ValueError('Subset type 12 requires an even number of projections!')
        half = nProjections // 2
        for kk in range(subsets):
            # Every nth projection from the first half, each followed by its
            # counterpart (half a rotation later) from the second half
            first = np.arange(kk, half, subsets)
            index1 = np.column_stack((first, first + half)).ravel().astype(tyyppi)
            temp_ind.append(index1)
            lengths[kk] = np.size(index1)
    else:
        index1 = np.arange(0, nProjections).astype(tyyppi)
        temp_ind.append(index1)
        lengths[0] = np.size(index1)
    index = np.concatenate(temp_ind) if temp_ind else np.empty(0, dtype=tyyppi)
    return index, lengths


def indexMaker(options):
    """
    This function is used to generate the subset indices for the selected 
    subsetType. Subset types 6 and 7 are supported only for sinogram
    data (not raw/list-mode data) and type 6 can change options.subsets (the
    number of subsets is then 180 / n_angles for 2D data or 360 / n_angles for
    3D data).

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Raises
    ------
    ValueError
        Invalid subset type.

    Returns
    -------
    None.

    """
    subsets = options.subsets
    if options.sampling > 1:
        options.Ndist = options.Ndist * options.sampling
    Ndist = options.Ndist
    Nang = options.Nang
    NSinos = options.NSinos
    options.nMeas = np.array(options.Ndist * options.Nang * options.NSinos, dtype=np.int64)
    if options.CT:
        tyyppi = np.uint64
    else:
        tyyppi = np.uint32
    dynamicProjectionSubsets = (
        options.listmode == 0 and options.Nt > 1 and
        ((subsets > 1 and options.subsetType in [8, 9, 10, 11, 12]) or subsets == 1) and
        hasattr(options, 'nProjectionsPerFrame') and
        np.size(options.nProjectionsPerFrame) == options.Nt
    )
    dynamicListmodeSubsets = (
        options.listmode > 0 and options.Nt > 1 and
        np.size(getattr(options, 'listmodeIndices', [])) == options.Nt and
        options.subsetType in (0, 1, 3)
    )
    if dynamicProjectionSubsets:
        projectionCounts = np.asarray(options.nProjectionsPerFrame, dtype=float).reshape(-1)
        if np.any((projectionCounts < 1) | (projectionCounts != np.fix(projectionCounts))):
            raise ValueError('The number of projection images in every timeframe must be a positive integer!')
        projectionCounts = projectionCounts.astype(np.int64)
        options.projectionFrameOffsets = np.concatenate(
            (np.asarray([0], dtype=np.int64), np.cumsum(projectionCounts, dtype=np.int64))
        )
        options.index = [None] * options.Nt
        pituus = np.zeros((subsets, options.Nt), dtype=np.int64)
        generator = None
        if options.subsetType == 9:
            if options.seed < 0:
                generator = np.random.default_rng()
            else:
                generator = np.random.default_rng(options.seed)
        if subsets > 1 and options.subsetType == 10:
            raise ValueError('Subset type 10 is not supported with dynamic projection subsets (timeframes with a separate number of projections)!')
        usedFallbackArm = not (subsets > 1 and options.subsetType in (8, 9, 10, 11, 12))
        for tt in range(options.Nt):
            nProjections = int(projectionCounts[tt])
            index1, lengths = _projection_subset_indices(nProjections, subsets, options.subsetType, generator, tyyppi)
            if usedFallbackArm and not (options.CT or options.PET or options.SPECT):
                lengths[0] = options.Nang * options.Ndist * options.NSinos
            options.index[tt] = index1
            pituus[:, tt] = lengths
        options.nMeas = pituus.flatten(order='F')
    elif dynamicListmodeSubsets:
        eventCounts = np.asarray(options.listmodeIndices, dtype=np.int64).reshape(-1)
        if np.any(eventCounts < 0):
            raise ValueError('Listmode event counts must be non-negative.')
        options.index = [None] * options.Nt
        pituus = np.zeros((subsets, options.Nt), dtype=np.int64)
        generator = np.random.default_rng(options.seed) if options.subsetType == 3 and options.seed >= 0 else np.random.default_rng()
        for tt, eventCount in enumerate(eventCounts):
            if options.subsetType == 1:
                frameIndices = [np.arange(kk, eventCount, subsets, dtype=tyyppi) for kk in range(subsets)]
                options.index[tt] = np.concatenate(frameIndices)
                for kk, indices in enumerate(frameIndices):
                    pituus[kk, tt] = indices.size
            elif options.subsetType == 3:
                frameIndices = np.array_split(generator.permutation(int(eventCount)).astype(tyyppi), subsets)
                options.index[tt] = np.concatenate(frameIndices)
                for kk, indices in enumerate(frameIndices):
                    pituus[kk, tt] = indices.size
            else:
                val = int(eventCount) // subsets
                pituus[:-1, tt] = val
                pituus[-1, tt] = int(eventCount) - val * (subsets - 1)
                options.index[tt] = np.zeros(1, dtype=tyyppi)
        options.nMeas = pituus.flatten(order='F')

    if subsets > 1 and options.subsetType < 8 and not dynamicListmodeSubsets:
        totalLength = Ndist*Nang*NSinos
        options.index = np.empty(0, dtype=tyyppi)
        options.nMeas = np.zeros((subsets, 1), dtype = np.int64)
        if options.subsetType == 4 and options.use_raw_data == 0:
            for i in range(subsets):
                osa = np.size(np.arange(i, Nang, subsets))
                index1 = np.tile(np.arange(i*Ndist, (i + 1) * Ndist, 1), (osa*NSinos)).astype(tyyppi)
                index1 = index1 + np.repeat(np.arange(0, (osa)*NSinos, 1) * Ndist*subsets,Ndist).astype(tyyppi)
                if Nang % subsets > 0:
                    erotus = Nang - osa * subsets
                    index1 = (np.int64(index1) + np.int64(np.repeat((np.arange(0, NSinos, 1)) * Ndist * erotus, Ndist * osa))).astype(tyyppi)
                options.index = np.append(options.index, index1)
                options.nMeas[i] = np.size(index1)
        elif options.subsetType == 5:
            apu = np.tile((np.tile(np.arange(0, Nang * Ndist, Ndist, dtype=tyyppi).T[None,:], (Ndist, 1)) + np.arange(0, Ndist, dtype=tyyppi)[:,None])[:,:,None], (1, 1, NSinos)) + np.transpose(np.arange(0, Nang * Ndist * NSinos, Nang * Ndist, dtype=tyyppi)[:,None,None],(2, 1, 0))
            apu = np.transpose(apu, (0, 2, 1)).reshape(Ndist * NSinos, -1, order='F')
            for i in range(subsets):
                index1 = apu[i : : subsets, :]
                index1 = index1.T
                index1 = index1.ravel('F')
                options.index = np.append(options.index, index1)
                options.nMeas[i] = np.size(index1)
            # Take every nth (column) measurement
        elif options.subsetType == 1:
            if options.listmode > 0 and options.Nt > 1:
                options.index = [None] * options.Nt
                pituus = np.zeros((subsets, options.Nt), dtype=np.int64)
                for j in range(options.Nt):
                    tot_length = options.listmodeIndices[j]
                    temp_ind = []
                    for i in range(subsets):
                        index1 = np.arange(i, tot_length, subsets, dtype=tyyppi)
                        temp_ind.append(index1)
                        pituus[i, j] = len(index1)
                    options.index[j] = np.concatenate(temp_ind)
                options.nMeas = pituus.flatten(order='F')
            else:
                for i in range(subsets):
                    index1 = np.arange(i, totalLength, subsets).astype(tyyppi)
                    [I,J,K] = ind2sub(index1, np.array((Nang, Ndist, NSinos)))
                    index1 = sub2ind(np.array((Ndist, Nang, NSinos)), J, I,K)
                    options.index = np.append(options.index, index1)
                    options.nMeas[i] = np.size(index1)
            # Take every nth (row) measurement
            # Every nth measurements
        elif options.subsetType == 2:
            for i in range(subsets):
                index1 = np.arange(i, totalLength, subsets).astype(tyyppi)
                options.index = np.append(options.index, index1)
                options.nMeas[i] = np.size(index1)
            # Pick the measurements randomly
        elif options.subsetType == 3:
            port = totalLength // subsets
            if options.seed < 0:
                generator = np.random.default_rng()
            else:
                generator = np.random.default_rng(options.seed)
            apu = generator.permutation(totalLength).astype(tyyppi)
            for i in range(subsets):
                if i == subsets - 1:
                    index1 = apu[port * i :]
                else:
                    index1 = apu[port * i : port * (i + 1)]
                options.index = np.append(options.index, index1)
                options.nMeas[i] = np.size(index1)
            # Pick the subsets based on the angles of the LORs
        elif options.subsetType == 6:
            if options.use_raw_data:
                raise ValueError('Subset type 6 is not supported with raw data in Python! Use sinogram data.')
            options.index, lengths = _subsetAngles(options, tyyppi)
            subsets = lengths.size
            if options.NSinos < options.TotSinos:
                # Only the direct planes are used (180 / n_angles subsets)
                if options.NSinos != int(np.asarray(options.Nz).ravel()[0]):
                    raise ValueError('Number of sinograms with subset type 6 has to be either ' + str(int(np.asarray(options.Nz).ravel()[0])) + ' or ' + str(options.TotSinos))
            options.nMeas = lengths.astype(np.int64).reshape((subsets, 1))
            # Use golden angle sampling
        elif options.subsetType == 7:
            if options.use_raw_data:
                raise ValueError('Subset type 7 is not supported with raw data in Python! Use sinogram data.')
            options.index, lengths = _goldenAngleSubsets(options, tyyppi)
            options.nMeas = lengths.astype(np.int64).reshape((subsets, 1))
        elif options.subsetType == 0:
            if options.listmode and options.Nt > 1:
                pituus = np.zeros((subsets, options.Nt), dtype=np.int64)
                for j in range(options.Nt):
                    totLength  = options.listmodeIndices[j]
                    val = totLength  // subsets
                    if totLength  % subsets > 0:
                        valEnd = totLength  - val * (subsets - 1)
                    else:
                        valEnd = val
                    for i in range(subsets - 1):
                        pituus[i, j] = val
                    pituus[subsets - 1, j] = valEnd
                options.nMeas = pituus.flatten(order='F')
                options.index = np.zeros(1,dtype=tyyppi)
            else:
                if options.listmode:
                    val = totalLength // subsets
                    if totalLength % subsets > 0:
                        valEnd = totalLength - val * (subsets - 1)
                    else:
                        valEnd = val
                else:
                    val = options.nProjections // subsets
                    if options.nProjections % subsets > 0:
                        valEnd = options.nProjections - val * (subsets - 1)
                    else:
                        valEnd = val
                options.nMeas = np.full(subsets - 1, val, dtype=np.int64)
                options.nMeas = np.append(options.nMeas, valEnd)
                options.index = np.zeros(1,dtype=tyyppi)
    elif ((subsets > 1 and options.subsetType in [8, 9, 10, 11, 12]) or subsets == 1) and not dynamicProjectionSubsets:
        generator = None
        if options.subsetType == 9 and subsets > 1:
            if options.seed < 0:
                generator = np.random.default_rng()
            else:
                generator = np.random.default_rng(options.seed)
        options.index, lengths = _projection_subset_indices(options.nProjections, subsets, options.subsetType, generator, tyyppi, options.angles)
        options.nMeas = lengths.reshape((subsets, 1))
        if options.subsets == 1:
            if options.listmode > 0 and options.Nt > 1 and np.size(getattr(options, 'listmodeIndices', [])) == options.Nt:
                options.nMeas = np.asarray(options.listmodeIndices, dtype=np.int64).reshape(-1)
            elif (options.CT or options.PET or options.SPECT) and options.listmode == 0:
                options.nMeas[0] = options.NSinos
            elif options.listmode == 1:
                if options.useIndexBasedReconstruction:
                    options.nMeas[0] = options.trIndex.size // 2
                else:
                    options.nMeas[0] = options.x.size // 6
            else:
                options.nMeas[0] = options.Nang * options.Ndist * options.NSinos
    elif options.subsetType > 12:
        raise ValueError('Invalid subset type!')
    if options.listmode == 0 and options.Nt > 1 and not dynamicProjectionSubsets:
        options.nMeas = np.tile(options.nMeas, options.Nt)
    nMeas = np.asarray(options.nMeas, dtype=np.int64)
    if nMeas.size == subsets * options.Nt:
        if nMeas.ndim == 2 and nMeas.shape == (subsets, options.Nt):
            options.nMeasPerFrameSubset = nMeas.T.copy()
        else:
            options.nMeasPerFrameSubset = nMeas.reshape((subsets, options.Nt), order='F').T
    else:
        options.nMeasPerFrameSubset = nMeas.reshape((subsets, 1), order='F').T
    if options.sampling > 1:
        options.Ndist = int(options.Ndist / options.sampling)
    options.subsets = subsets
    
def ind2sub(index, dims):
    """
    Transforms linear indices to subscripts.

    Parameters
    ----------
    index : NumPy array
        Linear indices.
    dims : Tuple
        Size for x, y and z dimensions.

    Returns
    -------
    I : NumPy array
        First dimension subscripts.
    J : NumPy array
        Second dimension subscripts.
    K : NumPy array
        Third dimension subscripts.

    """
    K = index // (dims[0] * dims[1])
    J = (index - K * (dims[0] * dims[1])) // dims[0]
    I = np.mod(index - K * (dims[0] * dims[1]), dims[0])
    return (I, J, K)
    
def sub2ind(dims, I, J, K):
    """
    Same as above, but the other way around.

    Parameters
    ----------
    dims : Tuple
        Size for x, y and z dimensions.
    I : NumPy array
        First dimension subscripts.
    J : NumPy array
        Second dimension subscripts.
    K : NumPy array
        Third dimension subscripts.

    Returns
    -------
    index : NumPy array
        Linear indices.

    """
    index = K * dims[0] * dims[1] + I + J * dims[0]
    return index


def formSubsetIndices(options):
    """
    Handles some of the subset indexing. Mainly reorganizes the index-based
    reconstruction indices to match the subsets but also the UV-direction 
    vectors for CT. Also, computes the "special" indices for subset type 3.

    Parameters
    ----------
    options : class object
        OMEGA class object used to contain all the necessary data.

    Returns
    -------
    None.

    """
    if options.sampling > 1:
        options.Ndist = options.Ndist * options.sampling
    if options.listmode > 0 and options.subsetType < 8 and options.subsets > 1:
        if options.subsetType > 0:
            if options.useIndexBasedReconstruction:
                if options.Nt > 1:
                    if isinstance(options.trIndex, list):
                        for kk in range(options.Nt):
                            if options.trIndex[kk].shape[0] != 2:
                                options.trIndex[kk] = options.trIndex[kk].reshape(2, -1, order='F')
                            options.trIndex[kk] = options.trIndex[kk][:, options.index[kk]]
                            if options.axIndex[kk].shape[0] != 2:
                                options.axIndex[kk] = options.axIndex[kk].reshape(2, -1, order='F')
                            options.axIndex[kk] = options.axIndex[kk][:, options.index[kk]]
                    else:
                        if isinstance(options.index, list):
                           options.index = np.concatenate(options.index)
                        options.trIndex = options.trIndex[:, options.index, :]
                        options.axIndex = options.axIndex[:, options.index, :]
                else:
                    if options.trIndex.shape[0] != 2:
                        options.trIndex = np.reshape(options.trIndex, (2, -1), order='F')
                    options.trIndex = options.trIndex[:,options.index]
                    if options.axIndex.shape[0] != 2:
                        options.axIndex = np.reshape(options.axIndex, (2, -1), order='F')
                    options.axIndex = options.axIndex[:,options.index]
            else:
                if options.Nt > 1:
                    if isinstance(options.x, list):
                        for kk in range(options.Nt):
                            if options.x[kk].shape[0] != 6:
                                options.x[kk] = options.x[kk].reshape(6, -1, order='F')
                            options.x[kk] = options.x[kk][:, options.index[kk]]
                            if options.SPECT and isinstance(options.z, list):
                                if options.z[kk].shape[0] != 5:
                                    options.z[kk] = options.z[kk].reshape(5, -1, order='F')
                                options.z[kk] = options.z[kk][:, options.index[kk]]
                    else:
                        if isinstance(options.index, list):
                            options.index = np.concatenate(options.index)
                        options.x = options.x[:, options.index, :]
                        if options.SPECT and not isinstance(options.z, list):
                            z = np.asarray(options.z)
                            if z.ndim == 1:
                                z = z.reshape(5, -1, order='F')
                            options.z = z[:, options.index]
                else:
                    if not options.x.shape[0] == 6:
                        options.x = np.reshape(options.x, (6, options.x.size // 6), order='F')
                    options.x = options.x[:,options.index]
                    if options.SPECT:
                        z = np.asarray(options.z)
                        if z.ndim == 1:
                            z = z.reshape(5, -1, order='F')
                        options.z = z[:, options.index]
            if options.TOF_bins_used > 1:
                if options.Nt > 1:
                    if isinstance(options.TOFIndices, list):
                        for kk in range(options.Nt):
                            options.TOFIndices[kk] = options.TOFIndices[kk][options.index[kk]]
                    else:
                        if isinstance(options.index, list):
                            options.index = np.concatenate(options.index)
                        options.TOFIndices = options.TOFIndices[options.index]
                else:
                    options.TOFIndices = options.TOFIndices[options.index]
        if not isinstance(options.x, list):
            options.x = options.x.ravel('F')
        if options.useIndexBasedReconstruction and options.listmode > 0 and isinstance(options.trIndex, list):
            for kk in range(options.Nt):
                options.trIndex[kk] = options.trIndex[kk].flatten(order='F')
                options.axIndex[kk] = options.axIndex[kk].flatten(order='F')
            options.trIndex = np.concatenate(options.trIndex)
            options.axIndex = np.concatenate(options.axIndex)
            if options.TOF_bins > 1 and isinstance(options.TOFIndices, list):
                options.TOFIndices = np.concatenate(options.TOFIndices)
        elif not options.useIndexBasedReconstruction and options.listmode > 0 and isinstance(options.x, list):
            for kk in range(options.Nt):
                options.x[kk] = options.x[kk].flatten(order='F')
                if options.SPECT and isinstance(options.z, list):
                    options.z[kk] = options.z[kk].flatten(order='F')
            options.x = np.concatenate(options.x)
            if options.SPECT and isinstance(options.z, list):
                options.z = np.concatenate(options.z)
            if options.TOF_bins > 1 and isinstance(options.TOFIndices, list):
                options.TOFIndices = np.concatenate(options.TOFIndices)
        options.xy_index = np.empty(0, dtype=np.uint32)
        options.z_index = np.empty(0, dtype=np.uint16)
    elif options.listmode > 0 and (options.subsetType in [8, 9]) and options.subsets > 1:
        if options.Nt > 1:
            options.z = options.z[options.index,:,:]
            if options.uV.shape[0] == options.nProjections:
                options.uV = options.uV[options.index,:,:]
        else:
            options.z = options.z[options.index,:]
            if options.uV.shape[0] == options.nProjections:
                options.uV = options.uV[options.index,:]
        options.xy_index = np.empty(0, dtype=np.uint32)
        options.z_index = np.empty(0, dtype=np.uint16)
    elif options.use_raw_data == 0 and (options.subsets > 1 and options.subsetType in [3, 6, 7]):
        options.xy_index = np.arange(0, options.Nang * options.Ndist, dtype=np.uint32)
        if options.span > 1:
            xy_index2 = np.tile(np.arange(0, options.Nang * options.Ndist, dtype=np.uint32), (options.NSinos - (options.rings * 2 - 1)))
            options.xy_index = np.concatenate((np.tile(options.xy_index,  (options.rings * 2 - 1)), xy_index2))
        else:
            xy_index2 = np.tile(np.arange(0, options.Nang * options.Ndist, dtype=np.uint32), (options.NSinos - options.rings))
            options.xy_index = np.concatenate((np.tile(options.xy_index, (options.rings)), xy_index2))
        options.z_index = np.arange(0, options.NSinos, dtype=np.uint16)
        options.z_index = np.repeat(options.z_index, (options.Nang * options.Ndist))
        options.z_index = options.z_index[options.index]

    if options.listmode > 0 and not options.useIndexBasedReconstruction and isinstance(options.x, list):
        options.x = np.concatenate([np.asarray(value).ravel(order='F') for value in options.x])
        if options.SPECT and isinstance(options.z, list):
            options.z = np.concatenate([np.asarray(value).ravel(order='F') for value in options.z])
    
        options.xy_index = options.xy_index[options.index]
    else:
        options.xy_index = np.empty(0, dtype=np.uint32)
        options.z_index = np.empty(0, dtype=np.uint16)
    if options.sampling > 1:
        options.Ndist = int(options.Ndist / options.sampling)
