# -*- coding: utf-8 -*-
"""
@author: Ville-Veikko Wettenhovi
"""

import numpy as np


def matlabRound(x):
    """Round half away from zero, matching MATLAB's round().

    MATLAB's round() rounds ties away from zero, e.g. round(4.5) == 5 and
    round(-4.5) == -5. Python's built-in round() and numpy.round()/
    numpy.around() instead round ties to the nearest even value (banker's
    rounding), so round(4.5) == 4. Several routines translated from the
    MATLAB m-files (e.g. the multi-resolution/extended-FOV volume sizing
    in setUpCorrections.m: NzM2 = round((axial_fov - axialFOVOrig) / 2 /
    dzM) * 2) rely on the away-from-zero tie-breaking behaviour of MATLAB's
    round(), so a direct substitution of round()/np.round() can silently
    produce a different result whenever the argument is exactly at a .5
    tie.

    Parameters
    ----------
    x : scalar (int/float) or array_like
        Value(s) to round.

    Returns
    -------
    float
        If ``x`` is a scalar (Python int/float or a 0-d/scalar NumPy
        value), a plain Python float is returned, mirroring how MATLAB's
        round() is used in the m-files (its result is typically assigned
        to a numeric variable and cast to an integer type, e.g. int() or
        uint32(), by the caller).
    numpy.ndarray
        If ``x`` is array_like (including a NumPy array), an ndarray of
        the same shape with dtype float64 is returned.
    """
    arr = np.asarray(x, dtype=np.float64)
    rounded = np.sign(arr) * np.floor(np.abs(arr) + 0.5)
    if np.ndim(x) == 0:
        return float(rounded)
    return rounded
