# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 14:36:12 2025

@author: Ville-Veikko Wettenhovi
"""

import numpy as np
from scipy.interpolate import interp2d

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

    x_vals = np.arange(1, options.Nang + 1)
    y_vals = np.arange(1, options.Ndist * options.sampling + 1, options.sampling)
    points = (y_vals, x_vals)
    
    yq_vals = np.arange(1, options.Ndist * options.sampling + 1)
    xq_vals = np.arange(1, options.Nang + 1)
    Yq, Xq = np.meshgrid(yq_vals, xq_vals, indexing='ij')
    query_points = np.stack([Yq.ravel(), Xq.ravel()], axis=-1)
    
    # Interpolation
    def interp_surface(data):
        from scipy.interpolate import interpn
        return np.asfortranarray((interpn(points, data.astype(np.float64), query_points, method=options.sampling_interpolation_method, bounds_error=False).reshape(Yq.shape)))
    
    x1 = interp_surface(xx1)
    x2 = interp_surface(xx2)
    y1 = interp_surface(yy1)
    y2 = interp_surface(yy2)
    
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


def interpolateSinog(SinM, sampling, Ndist, Nang, sampling_interpolation_method):
    from scipy.interpolate import interpn
    if SinM.ndim == 1 or SinM.shape[2] == 1:
        SinM = np.reshape(SinM, (Ndist, Nang, -1), order='F')
    if SinM.ndim == 3:
        SinM = np.reshape(SinM, (SinM.shape[0], SinM.shape[1], SinM.shape[2], 1), order='F')
    
    y = np.arange(1, Ndist * sampling + 1, sampling)
    x = np.arange(1, Nang + 1)
    grid = (np.asfortranarray(y), np.asfortranarray(x))

    yq = np.arange(1, Ndist * sampling + 1)
    xq = np.arange(1, Nang + 1)
    Xq, Yq = np.meshgrid(xq, yq, indexing='xy')
    Xq = Xq.T
    Yq = Yq.T
    query_points = np.asfortranarray(np.column_stack((Yq.ravel(), Xq.ravel())))
    
    SinM_uus = np.zeros((SinM.shape[0] * sampling, Nang, SinM.shape[2], SinM.shape[3]), dtype=np.float32, order='F')

    for uu in range(SinM.shape[3]):
        for kk in range(SinM.shape[2]):
            sinogram_slice = SinM[:, :, kk, uu].astype(np.float64)

            interpolated = interpn(
                grid,
                sinogram_slice,
                query_points,
                method=sampling_interpolation_method,
                bounds_error=False
            )
            interpolated = interpolated.reshape(Ndist * sampling, Nang, order='F')

            SinM_uus[:, :, kk, uu] = np.asfortranarray(interpolated)

    return SinM_uus.astype(np.float32).ravel('F')