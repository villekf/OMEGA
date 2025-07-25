# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 16:27:43 2025

@author: Ville-Veikko Wettenhovi
"""

from scipy.ndimage import convolve
from omegatomo.util.padding import padding
import numpy as np

def randoms_smoothing(randoms: np.ndarray, options):
    """
    Performs a moving mean smoothing on randoms or scatter data.

    Parameters:
        randoms (np.ndarray): Input randoms/scatter data (2D or 3D).

    Returns:
        np.ndarray: Smoothed data.
    """
    if options.verbose > 0:
        print("Beginning randoms/scatter smoothing")

    Ndx, Ndy, Ndz = 2, 2, 0

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