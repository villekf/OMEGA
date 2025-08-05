# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 17:37:24 2025

@author: Ville-Veikko Wettenhovi
"""

import numpy as np

def Randoms_variance_reduction(Randoms, options):
    """Applies 3D fan-sum variance reduction to input randoms."""
    from omegatomo.projector.detcoord import detectorCoordinates, sinogramCoordinates2D, sinogramCoordinates3D
    if options.verbose > 0:
        print("Starting Randoms variance reduction")

    # Get coordinates and convert mm to cm
    z = sinogramCoordinates3D(options) / 10.0
    detectors_x, detectors_y = detectorCoordinates(options)
    x, y = sinogramCoordinates2D(options, detectors_x, detectors_y)
    x, y = x / 10.0, y / 10.0
    detectors_x = detectors_x / 10.0
    detectors_y = detectors_y / 10.0

    detectors_ring = options.detectors // options.rings

    if Randoms.ndim == 3 and Randoms.shape[2] == 1:
        Randoms = Randoms.reshape((options.Ndist, options.Nang, options.TotSinos))

    sino_amount = int(np.sum(options.segment_table))
    coeff_matrix = np.zeros_like(Randoms, dtype=np.float32)
    det_num = np.zeros((options.Ndist, options.Nang, 2), dtype=np.int32, order='F')

    # Normalize z to start from zero
    z = z - np.min(z)
    ring = (np.round(z / z[1, 0]) + 1).astype(np.int32)

    # Determine each LOR's detector index
    for u in range(options.Ndist * options.Nang):
        i = u % options.Ndist
        j = u // options.Ndist

        for d in range(len(detectors_x)):
            if np.linalg.norm([x[u, 0] - detectors_x[d], y[u, 0] - detectors_y[d]]) < 1e-3:
                det_num[i, j, 0] = d + 1
                break
        for d in range(len(detectors_x)):
            if np.linalg.norm([x[u, 1] - detectors_x[d], y[u, 1] - detectors_y[d]]) < 1e-3:
                det_num[i, j, 1] = d + 1
                break

    # Vectorized construction
    testi1 = det_num[:, :, 0].flatten('F')
    testi2 = det_num[:, :, 1].flatten('F')

    ring1 = ring[:, 0].reshape(1, -1)
    ring2 = ring[:, 1].reshape(1, -1)

    testi1 = (np.reshape(testi1, (-1, 1)) + ((ring1 - 1) * detectors_ring)).astype(np.int32) - 1
    testi2 = (np.reshape(testi2, (-1, 1)) + ((ring2 - 1) * detectors_ring)).astype(np.int32) - 1

    randoms_flat = Randoms.flatten('F')
    size = max(testi1.flatten('F').max(), testi2.flatten('F').max()) + 1

    # Sum values into bins (like accumarray)
    randoms_det = np.zeros(size, dtype=np.float32)
    hits_det = np.zeros(size, dtype=np.float32)
    
    # Bin summation
    randoms_det1 = np.bincount(testi1.flatten('F'), weights=randoms_flat, minlength=size)
    randoms_det2 = np.bincount(testi2.flatten('F'), weights=randoms_flat, minlength=size)
    hits_det1 = np.bincount(testi1.flatten('F'), minlength=size)
    hits_det2 = np.bincount(testi2.flatten('F'), minlength=size)
    
    # Combine results
    randoms_det[:size] = randoms_det1 + randoms_det2
    hits_det[:size] = hits_det1 + hits_det2

    with np.errstate(divide='ignore', invalid='ignore'):
        randoms_det = np.divide(randoms_det, hits_det, out=np.zeros_like(randoms_det), where=hits_det > 0)

    # Compute mean inverse
    mean_det = np.nanmean(randoms_det[hits_det > 0])
    coeffs_detectors = np.zeros_like(randoms_det)
    coeffs_detectors[hits_det > 0] = mean_det / randoms_det[hits_det > 0]

    for k in range(sino_amount):
        r1 = ring[k, 0]
        r2 = ring[k, 1]
        for i in range(options.Ndist):
            for j in range(options.Nang):
                d1 = det_num[i, j, 0] + (r1 - 1) * detectors_ring - 1
                d2 = det_num[i, j, 1] + (r2 - 1) * detectors_ring - 1
                coeff_matrix[i, j, k] = coeffs_detectors[d1] * coeffs_detectors[d2]

    New_randoms = Randoms.astype(dtype=np.float32) * coeff_matrix

    # Scale back total counts
    total_original = np.reshape(np.sum(np.sum(Randoms, axis=0), axis = 0), (1, 1, -1))
    total_new = np.reshape(np.sum(np.sum(New_randoms, axis=0), axis = 0), (1, 1, -1))
    New_randoms *= total_original / total_new

    if options.verbose > 0:
        print("Randoms variance reduction completed")

    return New_randoms