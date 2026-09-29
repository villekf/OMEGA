# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 16:26:38 2025

@author: Ville-Veikko Wettenhovi
"""

import numpy as np

def padding(A: np.ndarray, sizeP, mode='symmetric'):
    """
    Pads the input array symmetrically or with zeros.

    Parameters:
        A (np.ndarray): Input 2D or 3D array.
        sizeP (list or tuple): Amount of padding in each direction. 
                               E.g., [2, 2] means padding 2 on each side.
        mode (str): 'symmetric' (default) or 'zeros'.

    Returns:
        np.ndarray: Padded array.
    """
    sizeP = list(sizeP)
    if len(sizeP) < 3:
        sizeP += [0] * (3 - len(sizeP))  # Ensure 3D compatibility

    if mode == 'symmetric':
        # padding.m's symmetric branch (lines 38-39) pads dim1 by sizeP(2)
        # and dim2 by sizeP(1), i.e. the axes are swapped.
        pad_width = [
            (sizeP[1], sizeP[1]),  # pad rows (y)
            (sizeP[0], sizeP[0]),  # pad cols (x)
        ]
    elif mode == 'zeros':
        # padding.m's zeros branch (line 46) pads dim1 by sizeP(1) and dim2
        # by sizeP(2), i.e. NOT swapped (unlike the symmetric branch above).
        pad_width = [
            (sizeP[0], sizeP[0]),  # pad rows (y)
            (sizeP[1], sizeP[1]),  # pad cols (x)
        ]
    else:
        raise ValueError("Mode must be 'symmetric' or 'zeros'")

    if A.ndim >= 3:
        pad_width.append((sizeP[2], sizeP[2]))  # pad depth (z)
    if A.ndim == 4:
        pad_width.append((0, 0))  # do not pad batch/channel dimension

    if mode == 'symmetric':
        A = np.pad(A, pad_width, mode='symmetric')
    else:
        A = np.pad(A, pad_width, mode='constant', constant_values=0)

    return A