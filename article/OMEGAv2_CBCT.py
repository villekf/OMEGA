# -*- coding: utf-8 -*-

from CBCT_example import CBCTExample

# Used data available from: https://doi.org/10.5281/zenodo.12722386

# Set the path (to the folder) of the above mat-file to here:
path = ''

# type 0 = PDHG with filtering
# type 1 = PDHG with filtering and NLRDP regularization
# type 2 = PKMA with RDP regularization
iType = 0

# The selected GPU device
GPUDevice = 0

# If your GPU doesn't have enough memory to contain all the data set the
# below variable to true. It divides the image into segments that are
# reconstructed separately. This slows down the computations, but allows
# them to be computed even on GPUs with very limited amount of memory
# 4 GB of memory is recommended even when this is set as true, 12 GB
# otherwise
largeDim = True
z, t = CBCTExample(path, iType, GPUDevice, largeDim)