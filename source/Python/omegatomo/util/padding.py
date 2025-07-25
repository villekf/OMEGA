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

    pad_width = [
        (sizeP[1], sizeP[1]),  # pad rows (y)
        (sizeP[0], sizeP[0]),  # pad cols (x)
    ]

    if A.ndim >= 3:
        pad_width.append((sizeP[2], sizeP[2]))  # pad depth (z)
    if A.ndim == 4:
        pad_width.append((0, 0))  # do not pad batch/channel dimension

    if mode == 'symmetric':
        A = np.pad(A, pad_width, mode='symmetric')
    elif mode == 'zeros':
        A = np.pad(A, pad_width, mode='constant', constant_values=0)
    else:
        raise ValueError("Mode must be 'symmetric' or 'zeros'")
    
    return A