# -*- coding: utf-8 -*-
"""
Python port of scaleImages.m. Scales images to Bq/mL.
"""

import numpy as np

from ._saveUtil import _get, _first


def scaleImages(img, options=None, totalTime=None, voxelSize=None, frameDurations=None, timeAxis=None):
    """
    Scales the input image(s) to Bq/mL in the same way as scaleImages.m. By default OMEGA outputs images with the number of
    counts per voxel. The image is divided by the voxel volume in mL (dx * dy * dz with the voxel size in mm, i.e. 1 mL =
    1000 mm^3) and by the duration (s) of the time window the counts were collected in.

    Static data: the duration is the length of the time window, min(end, tot_time) - start, which is the total time of the
    measurement when the whole measurement is reconstructed. It has to be finite. If the options were created with an infinite
    total time, input the correct total time with totalTime.

    Dynamic data: each time step is divided by its own duration. The durations are taken, in order of priority, from
    frameDurations, from options.partitions when it is a vector of time-step durations (as it is after the data has been
    loaded) or, when the number of time steps is known (options.Nt, options.nTimeSteps or a scalar options.partitions), by
    dividing the time window equally. The time axis of the images is axis 3 (0-based) of an (Nx, Ny, Nz, Nt[, saves]) image
    as output by OMEGA, or axis 4 of a MATLAB-style (Nx, Ny, Nz, iterations, Nt) image (use timeAxis to override).

    Parameters
    ----------
    img : ndarray or list of ndarray
        The input image(s), any ND image or a list of images (None elements are kept as is).
    options : options object or dict, optional
        The options (or the image properties). Total time (tot_time), start and end times, time steps and the field of view
        and image sizes (FOVa_x/Nx, FOVa_y/Ny, axial_fov/Nz) are taken from this. For compatibility with
        scaleImages(pz, time) of MATLAB, a number can also be input here, it is then used as the length of the time window.
    totalTime : float, optional
        Total time of the measurement (s). Overrides the total time of the options.
    voxelSize : sequence of 3 floats, optional
        Distance between adjacent voxels (dx, dy, dz) in mm. Overrides the voxel size from the options.
    frameDurations : sequence of floats, optional
        Duration (s) of each time step, one value per time step. Overrides the durations from the options.
    timeAxis : int, optional
        The (0-based) time axis of the images when there is more than one time step. By default 3, or 4 for 5D images whose
        axis 3 is not the time axis but axis 4 matches the number of time steps.

    Returns
    -------
    ndarray or list of ndarray
        The scaled image(s). Floating point images keep their type.
    """
    if options is not None and (np.isscalar(options) or isinstance(options, np.ndarray)):
        if totalTime is None:
            totalTime = float(options)
        options = None
        optionsTime = False
    else:
        optionsTime = True
    # Voxel volume (mm^3)
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
    volume = dx * dy * dz

    durations = _frameDurations(options if optionsTime else None, totalTime, frameDurations)
    if isinstance(img, (list, tuple)):
        return [None if x is None else _scale(x, volume, durations, timeAxis) for x in img]
    return _scale(img, volume, durations, timeAxis)


def _frameDurations(options, totalTime, frameDurations):
    """Returns the duration(s) (s) of the time step(s) as a float64 array (one element for static data)."""
    if frameDurations is not None:
        d = np.asarray(frameDurations, dtype=np.float64).ravel()
    else:
        steps = _get(options, ('partitions', 'nTimeSteps'))
        if steps is not None and np.asarray(steps).size > 1:
            # A vector of time-step durations
            d = np.asarray(steps, dtype=np.float64).ravel()
        else:
            nSteps = None
            if steps is not None and np.asarray(steps).size == 1:
                nSteps = int(np.asarray(steps).ravel()[0])
            if nSteps is None:
                nSteps = _get(options, 'Nt')
                nSteps = None if nSteps is None else int(np.asarray(nSteps).ravel()[0])
            nSteps = 1 if nSteps is None or nSteps < 1 else nSteps
            tot = totalTime if totalTime is not None else _first(_get(options, ('totalTime', 'total_time', 'tot_time')))
            if tot is None:
                raise ValueError('Total time not found, input it with scaleImages(img, totalTime=...) or give the options')
            start = _first(_get(options, ('startTime', 'start')))
            end = _first(_get(options, ('endTime', 'end')))
            start = 0. if start is None else start
            end = np.inf if end is None else end
            windowEnd = min(end, float(tot))
            if not np.isfinite(windowEnd):
                raise ValueError('Total time is infinite, input finite total time with scaleImages(img, totalTime=total_time)')
            d = np.full(nSteps, (windowEnd - start) / nSteps, dtype=np.float64)
    if not np.all(np.isfinite(d)) or np.any(d <= 0):
        raise ValueError('The durations of the time steps have to be positive and finite')
    return d


def _scale(x, volume, durations, timeAxis):
    """Scales the image x to Bq/mL, volume is the voxel volume in mm^3 and durations the duration(s) of the time step(s) in s."""
    x = np.asarray(x)
    if x.dtype.kind != 'f':
        x = x.astype(np.float64)
    if durations.size == 1:
        return x / x.dtype.type(volume * durations[0] / 1000.)
    nt = durations.size
    if timeAxis is None:
        if x.ndim >= 4 and x.shape[3] == nt:
            timeAxis = 3
        elif x.ndim == 5 and x.shape[4] == nt:
            timeAxis = 4
        else:
            raise ValueError('The number of time steps (' + str(nt) + ' durations) does not match the image size ' + str(x.shape) +
                             '. The time axis has to be axis 3 (Nx, Ny, Nz, Nt[, saves]) or axis 4 (Nx, Ny, Nz, iterations, Nt).')
    elif x.ndim <= timeAxis or x.shape[timeAxis] != nt:
        raise ValueError('The size of axis ' + str(timeAxis) + ' of the image ' + str(x.shape) + ' does not match the number of time steps ' + str(nt))
    shape = [1] * x.ndim
    shape[timeAxis] = nt
    return x / (volume * durations / 1000.).reshape(shape).astype(x.dtype)
