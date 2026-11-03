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

The side volumes of the multi-resolution reconstruction have exactly the coarse voxel size, so their lengths are
rounded up to a whole number of coarse voxels: the extended FOV never shrinks and grows by less than one coarse voxel on
each outer face (see setUpCorrections of the projector class). The whole-grid inputs (initial image, BP mask, type 6
attenuation map) are given for the requested extended FOV, NxFull x NyFull x NzFull voxels of the main volume, and are
padded here with the enlargement (options.multiResPadLow/High, in voxels of the main volume).
"""
import numpy as np


def sideVoxels(length, d):
    """The number of coarse voxels (size d) of a side volume of the length `length`: rounded up, so that the extended FOV
    never shrinks. The tolerance keeps a length that is a whole number of voxels (up to rounding error) from getting one
    voxel more."""
    return int(np.ceil(float(length) / float(d) - 1e-6))


def setMultiResPadding(options, enlargeLow, enlargeHigh, FOVRequested, dxMain, dM):
    """Stores the enlargement of the extended FOV (mm, low and high side of x, y and z) as the padding (voxels of the main
    volume) of the whole-grid inputs, whether the padding is exact (the requested FOV and the enlargement are whole numbers of
    voxels of the main volume) and the size of the union of the volumes, and prints the enlargement."""
    enlargeLow = np.asarray(enlargeLow, dtype=np.float64)
    enlargeHigh = np.asarray(enlargeHigh, dtype=np.float64)
    dxMain = np.asarray(dxMain, dtype=np.float64)
    gridFull = np.asarray([int(np.asarray(options.NxFull).reshape(-1)[0]), int(np.asarray(options.NyFull).reshape(-1)[0]),
                           int(np.asarray(options.NzFull).reshape(-1)[0])], dtype=np.float64)
    options.multiResPadLow = enlargeLow / dxMain
    options.multiResPadHigh = enlargeHigh / dxMain
    options.multiResPadExact = bool(np.all(np.abs(options.multiResPadLow - np.round(options.multiResPadLow)) < 1e-3)
                                    and np.all(np.abs(options.multiResPadHigh - np.round(options.multiResPadHigh)) < 1e-3)
                                    and np.all(np.abs(gridFull * dxMain - np.asarray(FOVRequested, dtype=np.float64)) < 1e-3 * dxMain))
    options.multiResUnionDims = np.floor(gridFull + options.multiResPadLow + options.multiResPadHigh + 0.5).astype(np.int64)
    dM = np.asarray(dM, dtype=np.float64)
    if options.verbose > 0 and np.any(np.concatenate((enlargeLow, enlargeHigh)) > 1e-4 * np.tile(dM, 2)):
        print(f"The extended FOV was enlarged to a whole number of voxels ({dM[0]:.6g} x {dM[1]:.6g} x {dM[2]:.6g} mm) in the side volumes. "
              f"The enlargement (mm) in the low and high side of x is [{enlargeLow[0]:.6g}, {enlargeHigh[0]:.6g}], "
              f"y [{enlargeLow[1]:.6g}, {enlargeHigh[1]:.6g}] and z [{enlargeLow[2]:.6g}, {enlargeHigh[2]:.6g}]")


def padMultiResVolume(image, options, fillZero=False):
    """Extends a whole-grid 3D volume (NxFull x NyFull x NzFull) by the enlargement of the extended FOV, with zeros if
    fillZero, otherwise with the values of the edge voxels. Other sizes (and an enlargement that is not a whole number of
    voxels) give the volume as it is."""
    if not getattr(options, 'multiResPadExact', False):
        return image
    lo = np.round(np.asarray(options.multiResPadLow)).astype(np.int64)
    hi = np.round(np.asarray(options.multiResPadHigh)).astype(np.int64)
    dims = tuple(int(np.asarray(getattr(options, axis)).reshape(-1)[0]) for axis in ('NxFull', 'NyFull', 'NzFull'))
    if image.ndim != 3 or tuple(image.shape) != dims or (not lo.any() and not hi.any()):
        return image
    pad = [(int(a), int(b)) for a, b in zip(lo, hi)]
    if fillZero:
        return np.pad(image, pad, mode='constant', constant_values=0)
    return np.pad(image, pad, mode='edge')


def padType6Attenuation(options, attenuation):
    """Pads the type 6 attenuation map (x, y, z[, frame]) of the requested extended FOV by the enlargement of the extended FOV:
    the type 6 projector takes the map on the union bounding box of the volumes. A map that already is of the size of the union
    is returned as it is.

    The padding is air (mu = 0) transaxially (x and y), since the body is inside the transaxial FOV, but axially (z) the edge
    slice of the map is replicated: the body continues beyond the axial FOV, and an air layer would make the lowest coarse slice
    of the side volume an average of air and tissue, i.e. the activity near the axial edge would be under-attenuated (up to a
    few percent 20 mm from the edge and far more at the edge itself). With no axial enlargement the map is not changed."""
    dims = getattr(options, 'multiResUnionDims', None)
    if dims is None:
        return attenuation
    oldDims = tuple(int(np.asarray(getattr(options, axis)).reshape(-1)[0]) for axis in ('NxFull', 'NyFull', 'NzFull'))
    newDims = tuple(int(v) for v in np.asarray(dims).reshape(-1))
    if oldDims == newDims:
        return attenuation
    lo = np.maximum(np.floor(np.asarray(options.multiResPadLow) + 0.5).astype(np.int64), 0)
    hi = np.maximum(np.asarray(newDims) - np.asarray(oldDims) - lo, 0)
    attenuation = np.asarray(attenuation)
    nOld = int(np.prod(oldDims))
    if attenuation.size == 0 or attenuation.size % nOld != 0:
        return attenuation
    nFrames = attenuation.size // nOld
    isFlat = attenuation.ndim == 1
    volume = attenuation.reshape(oldDims + (nFrames,), order='F')
    # Air in x and y, the nearest edge slice in z (the identity index when there is no axial padding)
    padded = np.pad(volume, [(int(lo[0]), int(hi[0])), (int(lo[1]), int(hi[1])), (0, 0), (0, 0)], mode='constant', constant_values=0)
    padded = np.take(padded, np.clip(np.arange(newDims[2]) - int(lo[2]), 0, oldDims[2] - 1), axis=2)
    if isFlat:
        return padded.ravel(order='F')
    return padded.reshape(newDims + (nFrames,), order='F') if nFrames > 1 or attenuation.ndim == 4 else padded[..., 0]


def prepareType6Attenuation(options, attenuation):
    """The type 6 attenuation map of the multi-resolution reconstruction for the union bounding box of the volumes. A map that
    is of the size of the union, or of the requested extended FOV (NxFull x NyFull x NzFull, padded with air by the enlargement,
    see padType6Attenuation), is used as it is; a 3D or 4D (x, y, z[, frame]) map of any other size is first resized to the
    requested extended FOV (the whole grid, not the main volume). A flattened map cannot be resized and is used as it is."""
    attenuation = np.asarray(attenuation)
    dims = getattr(options, 'multiResUnionDims', None)
    if dims is None:
        return attenuation
    wholeGrid = tuple(int(np.asarray(getattr(options, axis)).reshape(-1)[0]) for axis in ('NxFull', 'NyFull', 'NzFull'))
    unionGrid = tuple(int(v) for v in np.asarray(dims).reshape(-1))
    if attenuation.ndim >= 3 and tuple(attenuation.shape[:3]) not in (wholeGrid, unionGrid):
        from scipy.ndimage import zoom
        print('Error: Attenuation data is of different size than the extended FOV (all the volumes). Attempting resize!')
        factors = tuple(w / s for w, s in zip(wholeGrid, attenuation.shape[:3])) + (1.,) * (attenuation.ndim - 3)
        attenuation = zoom(attenuation, factors).astype(attenuation.dtype, copy=False)
    return padType6Attenuation(options, attenuation)
