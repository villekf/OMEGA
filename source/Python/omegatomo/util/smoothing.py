# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 16:27:43 2025

@author: Ville-Veikko Wettenhovi
"""

import numpy as np

def randoms_smoothing(randoms: np.ndarray, options):
    from omegatomo.util.padding import padding
    from scipy.ndimage import convolve
    """
    Performs a moving mean smoothing on randoms or scatter data.

    Parameters:
        randoms (np.ndarray): Input randoms/scatter data (2D or 3D).

    Returns:
        np.ndarray: Smoothed data.
    """
    if options.verbose > 0:
        print("Beginning randoms/scatter smoothing")

    # The 7x7 moving-mean below is a spatial (Ndist x Nang) smoothing and is
    # only meaningful when `randoms` genuinely has that 2-D sinogram/detector
    # structure (e.g. a real sinogram, or index-based reconstruction fed a
    # sinogram-shaped custom input). Pure list-mode data -- one randoms/scatter
    # value per event, with no Ndist/Nang extent -- has no spatial neighbors to
    # average over, so averaging consecutive entries would just mix unrelated
    # events together. OMEGA's own listmode/useIndexBasedReconstruction setup
    # collapses Ndist and Nang to 1 for that case (see proj.py's listmode
    # branches), so checking the first two axes distinguishes the two cases
    # without needing to inspect options.listmode/useIndexBasedReconstruction
    # directly (both a raw event-list and an index-based-without-sinogram input
    # end up with the same degenerate Ndist=Nang=1 shape). This is
    # intentionally conservative: skip smoothing (with a warning) whenever
    # spatial structure cannot be established, rather than risk smoothing
    # across unrelated list-mode events.
    randoms = np.asarray(randoms)
    if randoms.ndim < 2 or randoms.shape[0] <= 1 or randoms.shape[1] <= 1:
        print("Warning: randoms/scatter smoothing was requested, but the input data does not "
              "have a detectable 2-D (Ndist x Nang) sinogram structure -- it looks like flat "
              "list-mode/event-based data instead. Skipping smoothing for this array, since a "
              "spatial moving-mean average is not meaningful for purely event-based data.")
        return randoms

    Ndx, Ndy, Ndz = 7, 7, 0

    if Ndz == 0:
        kernel = np.ones((Ndx, Ndy, 1)) / (Ndx * Ndy)
    else:
        kernel = np.ones((Ndx, Ndy, Ndz)) / (Ndx * Ndy * Ndz)

    pad_size = [Ndx // 2, Ndy // 2, Ndz // 2]
    padded = padding(randoms, pad_size, mode='symmetric').astype(np.float32)
    smoothed = convolve(padded, kernel, mode='constant', cval=0.0)

    # Trim the padding
    smoothed = smoothed[pad_size[0]:smoothed.shape[0] - pad_size[0], pad_size[1]:smoothed.shape[1] - pad_size[1], pad_size[2]:smoothed.shape[2] - pad_size[2]]

    if options.verbose > 0:
        print("Smoothing complete")

    return smoothed