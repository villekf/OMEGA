# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 14:36:12 2025

@author: Ville-Veikko Wettenhovi
"""

import numpy as np

def increaseSampling(options, x=None, y=None, interpolate_sinogram=False):
    from omegatomo.projector.detcoord import detectorCoordinates, sinogramCoordinates2D
    # If x is empty, get from detector coordinates
    if x is None or len(x) == 0:
        x, y = detectorCoordinates(options)
        x, y = sinogramCoordinates2D(options, x, y)

    xx1 = x[:, 0].reshape((options.Ndist, options.Nang), order='F')
    xx2 = x[:, 1].reshape((options.Ndist, options.Nang), order='F')
    yy1 = y[:, 0].reshape((options.Ndist, options.Nang), order='F')
    yy2 = y[:, 1].reshape((options.Ndist, options.Nang), order='F')

    # MATLAB's interp2 based interpolation, see _interpolateRows
    x1 = _interpolateRows(xx1, options.sampling, options.sampling_interpolation_method)
    x2 = _interpolateRows(xx2, options.sampling, options.sampling_interpolation_method)
    y1 = _interpolateRows(yy1, options.sampling, options.sampling_interpolation_method)
    y2 = _interpolateRows(yy2, options.sampling, options.sampling_interpolation_method)
    
    x = np.asfortranarray(np.column_stack((x1.ravel('F'), x2.ravel('F'))))
    y = np.asfortranarray(np.column_stack((y1.ravel('F'), y2.ravel('F'))))

    if interpolate_sinogram:
        if options.SinM.ndim == 1:
            options.SinM = options.SinM.reshape((options.Ndist, options.Nang, -1), order='F')

        options.SinM = interpolateSinog(options.SinM, options.sampling, options.Ndist, options.Nang, options.sampling_interpolation_method)
        options.Ndist *= options.sampling
        options.nRowsD *= options.sampling

        if options.verbose:
            print(f"Sinogram sampling increased by {options.sampling}x")

    return x, y, options


def _interpolateRows(data, sampling, method):
    """Interpolates a 2D (Ndist, Nang) array to Ndist * sampling rows, following MATLAB's interpolateSinog."""
    from scipy.interpolate import interpn
    if method == 'makima':
        raise ValueError("Sampling interpolation method 'makima' is not available in Python. Use 'pchip' or 'spline' instead.")
    spline = method == 'spline'
    # MATLAB's spline is a cubic spline
    scipyMethod = 'cubic' if spline else method
    data = np.asarray(data, dtype=np.float64)
    Ndist, Nang = data.shape
    xGrid = np.arange(1, Nang + 1)
    yq = np.arange(1, Ndist * sampling + 1)
    if spline:
        # No mirrored row, the last sampling - 1 query rows are extrapolated
        yGrid = np.arange(1, Ndist * sampling + 1, sampling)
        fillValue = None
    else:
        # Append a mirrored row and extend the grid by one sample. The extra last query row is not requested.
        data = np.concatenate((data, data[-2:-1, :]), axis=0)
        yGrid = np.arange(1, Ndist * sampling + 2, sampling)
        fillValue = np.nan
    Yq, Xq = np.meshgrid(yq, xGrid, indexing='ij')
    query_points = np.stack([Yq.ravel(), Xq.ravel()], axis=-1)
    out = interpn((yGrid, xGrid), data, query_points, method=scipyMethod, bounds_error=False, fill_value=fillValue)
    return np.asfortranarray(out.reshape(Yq.shape))


def interpolateSinog(SinM, sampling, Ndist, Nang, sampling_interpolation_method):
    if SinM.ndim < 3 or SinM.shape[2] == 1:
        SinM = np.reshape(SinM, (Ndist, Nang, -1), order='F')
    if SinM.ndim == 3:
        SinM = np.reshape(SinM, (SinM.shape[0], SinM.shape[1], SinM.shape[2], 1), order='F')
    
    SinM_uus = np.zeros((SinM.shape[0] * sampling, Nang, SinM.shape[2], SinM.shape[3]), dtype=np.float32, order='F')

    for uu in range(SinM.shape[3]):
        for kk in range(SinM.shape[2]):
            sinogram_slice = SinM[:, :, kk, uu].astype(np.float64)

            interpolated = _interpolateRows(sinogram_slice, sampling, sampling_interpolation_method)

            SinM_uus[:, :, kk, uu] = interpolated

    return SinM_uus.astype(np.float32).ravel('F')