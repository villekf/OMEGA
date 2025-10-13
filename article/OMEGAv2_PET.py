# -*- coding: utf-8 -*-

from PET_TOF_example import PETTOFExample

# Data available from: https://doi.org/10.5281/zenodo.17185907

# Set the path of the above mat-file to here:
path = ''

# iType 0 = OSEM with improved version of Siddon ray-based projector
# iType 1 = OSEM with volume of intersection ray-based projector
# iType 2 = list-mode PET, PKMA with relative difference prior regularization and with improved version of Siddon ray-based projector
iType = 0

# The selected GPU device
GPUDevice = 0

z, t = PETTOFExample(path, iType, GPUDevice)