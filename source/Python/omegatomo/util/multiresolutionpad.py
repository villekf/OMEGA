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
each outer face (see setUpCorrections of the projector class). The whole-grid inputs (initial image and BP mask) are
given for the requested extended FOV, NxFull x NyFull x NzFull voxels of the main volume, and are padded here with the
enlargement (options.multiResPadLow/High, in voxels of the main volume). Image-domain attenuation is prepared by the
shared attenuation routine in reconstruction/prepass.py.
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
