# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 2025

@author: Ville-Veikko Wettenhovi
"""

import os


def opencl_header_dir():
    """
    Returns the directory containing the OpenCL/CUDA/HIP kernel header and
    source files (e.g. general_opencl_functions.h, auxKernels.cl), with a
    trailing '/'.

    This mirrors the directory-selection logic used in
    projector/init.py: a PyPI install ships the kernels under
    omegatomo/opencl (detected via the presence of
    omegatomo/util/usingPyPi.py), while a source checkout keeps them
    under source/opencl, three directories above omegatomo/util.

    Returns
    -------
    str
        Absolute path to the kernel directory, with a trailing '/'.
    """
    fPath = os.path.dirname(__file__)
    if os.path.exists(os.path.join(fPath, 'usingPyPi.py')):
        headerDir = os.path.abspath(os.path.join(fPath, '..', 'opencl')) + "/"
    else:
        headerDir = os.path.abspath(os.path.join(fPath, '..', '..', '..', 'opencl')) + "/"
    return headerDir
