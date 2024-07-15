# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 16:58:41 2024

@author: Ville-Veikko Wettenhovi
"""

# import numpy as np
# import ctypes
# import math
# import tkinter as tk
# from tkinter.filedialog import askopenfilename
# import os
# from skimage.transform import resize #scikit-image
# from skimage.transform import downscale_local_mean
# from scipy import interpolate
# from skimage.measure import block_reduce
# from scipy.integrate import quad
# from scipy.special import ellipk, ellipe
# try:
#     import numpy_groupies as npg
# except ModuleNotFoundError:
#     print('numpy_groupies package not found! PET support will be limited to only custom detector/source coordinates! Install numpy_groupies with "pip install numpy_groupies" to enable full support for PET data.')
# try:
#     from pymatreader import read_mat
# except ModuleNotFoundError:
#     print('pymatreader not found! MAT-files cannot be loaded! If you want to use MAT-files, install pymatreader with "pip install pymatreader"')
# try:
#     from SimpleITK import ReadImage as loadMetaImage
#     from SimpleITK import GetArrayFromImage
# except ModuleNotFoundError:
#     print('SimpleITK not found! MetaImage-files cannot be loaded! If you want to use MetaImage-files, install SimpleITK with "pip install SimpleITK"')
# try:
#     import arrayfire as af
# except ModuleNotFoundError:
#     print('ArrayFire package not found! ArrayFire features are not supported. You can install ArrayFire package with "pip install arrayfire".')

from .projector import proj


    
__version__ = "2.0.0"