# -*- coding: utf-8 -*-
"""
## Python codes for SPECT custom algorithm reconstruction
This example contains a simplified example for custom algorithm
reconstruction using projection SPECT data. Currently the support for
some of the additional features is limited. The default configuration
uses MLEM, but OSEM with subsets is also available.

Note that custom algorithm refers to your own algorithms and not the built-in 
algorithms. This example merely has the MLEM/OSEM algorithm shown as an example.
The forward and/or backward projections of OMEGA are utilized for the computation
of these algorithms.

This example uses PyTorch and CuPy and thus requires CUDA (and CuPy and PyTorch)!
"""
import numpy as np
from omegatomo import proj
import matplotlib as plt
from pydicom import dcmread

options = proj.projectorClass()

###########################################################################
###########################################################################
###########################################################################
########################### SCANNER PROPERTIES ############################
###########################################################################
###########################################################################
###########################################################################

# Header file location
dcm = dcmread('/path/to/DICOM')


# Load projection images
options.SinM = dcm.pixel_array
options.SinM = options.SinM.transpose((2, 1, 0))

options.SinM = options.SinM[:,64//4:128-64//4,:]

### Crystal thickness (mm)
options.cr_p = 9.525

### Transaxial FOV size (mm), this is the length of the x [horizontal] side
# of the FOV
# Note that with SPECT data using projector_type = 6, this is not exactly
# used as the FOV size but rather as the value used to compute the voxel
# size
options.FOVa_x = 4.7952*128

### Transaxial FOV size (mm), this is the length of the y [vertical] side
# of the FOV
options.FOVa_y = options.FOVa_x

### Axial FOV (mm)
# This is unused if projector_type = 6. Cubic voxels are always assumed!
options.axial_fov = 4.7952*128

# Number of rows in a projection image
options.nRowsD = options.SinM.shape[0]

# Number of columns in a projection image
options.nColsD = options.SinM.shape[1]

# Number of projections
options.nProjections = options.SinM.shape[2]

# Number of detector heads
options.nHeads = 2

startAngle1 = (float)(dcm.DetectorInformationSequence[0].StartAngle)
startAngle2 = (float)(dcm.DetectorInformationSequence[1].StartAngle)

angleIncrement = (float)(dcm.RotationInformationSequence[0].AngularStep)

# Rotation angles
options.angles = np.concatenate((np.arange(startAngle2, startAngle2 + angleIncrement * (options.nProjections / options.nHeads), angleIncrement),np.arange(startAngle1, startAngle1 + angleIncrement * (options.nProjections / options.nHeads - 1), angleIncrement)))

# Radial distance of the detector panel from the COR
options.radiusPerProj = np.concatenate((dcm.DetectorInformationSequence[0].RadialPosition, dcm.DetectorInformationSequence[1].RadialPosition))

### Scanner name
# Used for naming purposes (measurement data)
options.machine_name = 'Prospecta'
 
###########################################################################
###########################################################################
###########################################################################
###########################################################################
 

###########################################################################
###########################################################################
###########################################################################
########################### IMAGE PROPERTIES ##############################
###########################################################################
###########################################################################
###########################################################################
 
### Reconstructed image pixel count (X-direction)
# NOTE: Non-square image sizes (X- and Y-direction) may not work
options.Nx = 128

### Y-direction
options.Ny = 128

### Z-direction (number of slices) (axial)
options.Nz = 96

### Flip the image (in vertical direction)?
options.flip_image = True

### How much is the image rotated?
# You need to run the precompute phase again if you modify this
# NOTE: The rotation is done in the detector space (before reconstruction).
# This current setting is for systems whose detector blocks start from the
# right hand side when viewing the device from front.
# Positive values perform the rotation in clockwise direction
options.offangle = (3*np.pi)/2
 
###########################################################################
###########################################################################
###########################################################################
###########################################################################
 
 
 

###########################################################################
###########################################################################
###########################################################################
########################## COLLIMATOR PROPERTIES ##########################
###########################################################################
###########################################################################
###########################################################################

### Collimator-detector response function (CDRF)
# You can either input the (Gaussian) PSF filter, or the standard
# deviations for both transaxial and axial directions or simply the
# collimator parameters (see below) for an analytic solution for round (and
# hexagonal) holes (this may be unoptimal).

# If you want to compute the CDRF analytically, input the following values:
# Collimator hole length (mm)
options.colL = 24.05
# Collimator hole radius
options.colR = 1.11/2
# Distance from collimator to the detector
options.colD = 0
# Intrinsic resolution
options.iR = 3.8

# If you have the standard deviations for transaxial (XY) and axial (Z)
# directions, you can input them here instead of the above values (the
# dimensions need to be options.nProjections x options.Nx):
# Transaxial standard deviation
# options.sigmaXY = np.tile(0, (options.nProjection, options.Nx))
# Axial standard deviation
# options.sigmaZ = np.tile(0, (options.nProjection, options.Nx))

# Lastly, you can input the filter for the CDRF directly. This should be
# filterSizeXY x filterSizeZ x options.nProjections:
# options.gFilter = np.ones((10,10,options.nProjections), dtype=np.float32)
 
###########################################################################
###########################################################################
###########################################################################
###########################################################################



###########################################################################
###########################################################################
###########################################################################
######################## RECONSTRUCTION PROPERTIES ########################
###########################################################################
###########################################################################
###########################################################################
 
############################### PROJECTOR #################################
### Type of projector to use for the geometric matrix
# 6 = Rotation-based projector
# See the documentation for more information:
# https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
options.projector_type = 6
 
######################### RECONSTRUCTION SETTINGS #########################
### Number of iterations (all reconstruction methods)
options.Niter = 1

### Number of subsets (all excluding MLEM and subset_type = 6)
options.subsets = 1

### Subset type (n = subsets)
# 8 = Use every nth projection image
# 9 = Randomly select the full projection images
# 11 = Use prime factor sampling to select the full projection images
# Most of the time subset_type 8 is sufficient.
options.subsetType = 8

# Required for SPECT data
options.SPECT = True


options.addProjector()

# Assumes that PyTorch tensors are used as input to either forward or backward projections
options.useTorch = True

# Required for PyTorch
options.useCUDA = True

# Uses CuPy instead of PyCUDA (recommended)
options.useCuPy = True

# Compute forward projection with options * f
# Compute backprojection with options.T() * y

options.initProj()
import torch
d_f = torch.tensor(options.x0, device='cuda')
m = options.SinM.ravel('F')


"""
MLEM
"""
d_m = torch.tensor(m, device='cuda')
Sens = options.T() * torch.ones(d_m.numel(), dtype=torch.float32, device='cuda')
for it in range(options.Niter):
    fp = options * d_f
    bp = options.T() * (d_m / fp)
    d_f = d_f / Sens * bp
    
    

# """
# OSEM
# """
# d_m = [None] * options.subsets
# for k in range(options.subsets):
#     d_m[k] = torch.tensor(m[options.nTotMeas[k].item() : options.nTotMeas[k + 1].item()], device='cuda')
# for it in range(options.Niter):
#     for k in range(options.subsets):
#         # This is necessary when using subsets
#         # Alternative, call options.forwardProject(d_f, k) to use forward projection
#         # options.backwardProject(m, k) for backprojection
#         options.subset = k
#         fp = options * d_f
#         Sens = options.T() * torch.ones(d_m[k].numel(), dtype=torch.float32, device='cuda')
#         Sens[Sens <= 0] = options.epps
#         bp = options.T() * (d_m[k] / fp)
#         d_f = d_f / Sens * bp
    
f_np = d_f.cpu().numpy()
f_np = np.reshape(f_np, (options.Nx[0].item(), options.Ny[0].item(), options.Nz[0].item()), order='F')
plt.pyplot.imshow(f_np[:,:,20], vmin=0)