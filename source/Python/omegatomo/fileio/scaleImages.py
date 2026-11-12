# -*- coding: utf-8 -*-
"""
Python port of scaleImages.m. Scales images to Bq/mL.
"""

import numpy as np

from ._saveUtil import _get, _first, _totalTime


def scaleImages(img, options=None, totalTime=None, voxelSize=None):
    """
    Scales the input image(s) in the same way as scaleImages.m: the image is divided by the voxel volume (dx * dy * dz,
    with the voxel size in mm) times the total time of the measurement (s). By default OMEGA outputs images with the number of
    counts per voxel. The total time has to be finite. If the options were created with an infinite total time, input the
    correct total time with totalTime.

    Parameters
    ----------
    img : ndarray or list of ndarray
        The input image(s), any ND image or a list of images (None elements are kept as is).
    options : options object or dict, optional
        The options (or the image properties). Total time (tot_time) and the field of view and image sizes (FOVa_x/Nx,
        FOVa_y/Ny, axial_fov/Nz) are taken from this. For compatibility with scaleImages(pz, time) of MATLAB, a number can also
        be input here, it is then used as the total time.
    totalTime : float, optional
        Total time of the measurement (s). Overrides the total time of the options.
    voxelSize : sequence of 3 floats, optional
        Distance between adjacent voxels (dx, dy, dz) in mm. Overrides the voxel size from the options.

    Returns
    -------
    ndarray or list of ndarray
        The scaled image(s). Floating point images keep their type.
    """
    if options is not None and (np.isscalar(options) or isinstance(options, np.ndarray)):
        if totalTime is None:
            totalTime = float(options)
        options = None
    if totalTime is None:
        totalTime = _totalTime(options)
    if totalTime is None:
        raise ValueError('Total time not found, input it with scaleImages(img, totalTime=...) or give the options')
    if not np.isfinite(totalTime):
        raise ValueError('Total time is infinite, input finite total time with scaleImages(img, totalTime=total_time)')
    if voxelSize is not None:
        v = np.asarray(voxelSize, dtype=np.float64).ravel()
        if v.size != 3:
            raise ValueError('voxelSize has to contain three values (dx, dy, dz)')
        dx, dy, dz = float(v[0]), float(v[1]), float(v[2])
    else:
        fx = _first(_get(options, ('FOVa_x', 'FOVx', 'FOV_x')))
        fy = _first(_get(options, ('FOVa_y', 'FOVy', 'FOV_y')))
        fz = _first(_get(options, ('axial_fov', 'axialFOV', 'axial_FOV')))
        nx = _first(_get(options, 'Nx'))
        ny = _first(_get(options, 'Ny'))
        nz = _first(_get(options, 'Nz'))
        if any(x is None for x in (fx, fy, fz, nx, ny, nz)):
            raise ValueError('Voxel size not found, input it with voxelSize=(dx, dy, dz) or give the options')
        dx = fx / nx
        dy = fy / ny
        dz = fz / nz
    scale = dx * dy * dz * float(totalTime)
    if isinstance(img, (list, tuple)):
        return [None if x is None else _scale(x, scale) for x in img]
    return _scale(img, scale)


def _scale(x, scale):
    x = np.asarray(x)
    if x.dtype.kind != 'f':
        x = x.astype(np.float64)
    return x / x.dtype.type(scale)
