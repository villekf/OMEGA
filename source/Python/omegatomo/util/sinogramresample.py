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

Resampling of the detector images of the SPECT rotation-based projector (type 6) when the FOV is not a whole number of detector
pixels (FOV / detector pitch is not an integer).

The detector image is otherwise padded with zeros or trimmed symmetrically to the (integer) number of detector pixels that span the
FOV and then resized to the voxel lattice of the image. A fraction of a pixel cannot be padded (and rounding the number of pixels
scales and shifts the image), so a ratio that is not an integer (up to the relative tolerance RATIO_TOLERANCE) is resampled in
one step directly between the two lattices, whose centers coincide: the output pixel j (0-based, of nOut) is at the input index
u(j) = (j + 0.5 - nOut / 2) * ratio + nIn / 2 - 0.5, ratio = output pixel pitch / input pixel pitch. This is the same interpolation
as that of the resize of the integer case (skimage.transform.resize: linear, order 1, with the gaussian antialiasing of sigma =
(ratio - 1) / 2 when shrinking, or nearest for masks), so the values of the pixels and not the sums of the counts are preserved, and
the input is zero outside of the array (like the zero padding). It is a weighted mean of the input values: the weights of a pixel sum
to one (to less than one next to the edge of the array).

An integer number of FOV pixels is padded or trimmed as always, except that an axis where the difference between the detector size and
the number of FOV pixels is odd is resampled too: padding or trimming a odd number of pixels symmetrically is impossible, the extra
pixel goes to the end of the axis and the image is shifted by half a detector pixel, see exactAxisRatios.
"""
import numpy as np

RATIO_TOLERANCE = 1e-9


def snapRatio(ratio):
    """(ratio, True) with the ratio rounded to an integer if it is an integer up to the relative RATIO_TOLERANCE, else
    (ratio, False)."""
    ratio = float(ratio)
    nearest = float(np.round(ratio))
    if abs(ratio - nearest) <= RATIO_TOLERANCE * max(1.0, abs(ratio)):
        return nearest, True
    return ratio, False


def exactAxisRatios(ratios, measurement, physical):
    """The axes of the detector image (rows, columns) that have to be resampled directly between the detector and the image lattices.
    ratios are the FOV sizes in detector pixels (FOV / pitch, rows and columns), measurement the detector size (rows, columns) and
    physical the (integer) number of FOV pixels. Returns None if all axes can be padded/trimmed (and then resized) as always, else
    a tuple of the (snapped) ratio of each axis, or None for an axis that is padded/trimmed as always.

    An axis is resampled if its ratio is not an integer (up to RATIO_TOLERANCE), in which case both axes are resampled, or if the
    ratio is an integer but the difference of the detector size and the integer number of FOV pixels is odd: the padding/trimming
    would then be asymmetric, i.e. shift the image by half a detector pixel. An even difference is centered exactly."""
    snapped = [snapRatio(ratio) for ratio in ratios]
    if not all(flag for _, flag in snapped):
        return tuple(value for value, _ in snapped)
    exact = tuple(value if (int(measurement[axis]) - int(physical[axis])) % 2 else None for axis, (value, _) in enumerate(snapped))
    return None if all(value is None for value in exact) else exact


def axisMatrix(nIn, nOut, ratio, mode='continuous'):
    """The (nOut, nIn) matrix of the 1D resampling of nIn pixels to nOut pixels of the pitch `ratio` times the pitch of the input,
    centered on the center of the input. mode 'mask' is the nearest neighbour, otherwise linear interpolation with antialiasing
    when ratio > 1. The input is zero outside of the array."""
    nIn, nOut, ratio = int(nIn), int(nOut), float(ratio)
    matrix = np.zeros((nOut, nIn), dtype=np.float64)
    u = (np.arange(nOut) + 0.5 - nOut / 2.0) * ratio + nIn / 2.0 - 0.5
    if mode == 'mask':
        nearest = np.floor(u + 0.5).astype(np.int64)
        ok = (nearest >= 0) & (nearest < nIn)
        matrix[np.nonzero(ok)[0], nearest[ok]] = 1.0
        return matrix
    sigma = max(0.0, (ratio - 1.0) / 2.0)
    if sigma > 0.0:
        radius = int(4.0 * sigma + 0.5)
        offsets = np.arange(-radius, radius + 1)
        kernel = np.exp(-0.5 * (offsets / sigma) ** 2)
        kernel /= kernel.sum()
    else:
        offsets = np.zeros(1, dtype=np.int64)
        kernel = np.ones(1)
    base = np.floor(u).astype(np.int64)
    fraction = u - base
    for neighbour, weight in ((0, 1.0 - fraction), (1, fraction)):
        for offset, gauss in zip(offsets, kernel):
            index = base + neighbour + offset
            ok = (index >= 0) & (index < nIn)
            np.add.at(matrix, (np.nonzero(ok)[0], index[ok]), (weight * gauss)[ok])
    return matrix


def resampleAxis(array, axis, nOut, ratio, mode='continuous'):
    """Resamples the axis `axis` of the numpy array to nOut pixels, see axisMatrix."""
    array = np.asarray(array)
    matrix = axisMatrix(array.shape[axis], nOut, ratio, mode)
    moved = np.moveaxis(array.astype(np.float64, copy=False), axis, 0)
    result = np.tensordot(matrix, moved, axes=(1, 0))
    return np.moveaxis(result, 0, axis).astype(array.dtype if np.issubdtype(array.dtype, np.floating) else np.float32, copy=False)
