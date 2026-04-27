# -*- coding: utf-8 -*-
"""
Python codes for CBCT reconstruction
This example file computes the Figure 4 of the OMEGA V2 article with FDK.
DOI will be added later.
You need to provide the folder path to the input data
Planmeca_VisoG7_100kV_80mAs_500proj_kneePhantom.mat!
Used data available from: https://doi.org/10.5281/zenodo.12722386
"""
import numpy as np
import os
from omegatomo.projector import proj
from omegatomo.util import CTEFOVCorrection
from omegatomo.reconstruction import reconstructions_mainCT
from omegatomo.util.checkCUDA import checkCUDA
import matplotlib.pyplot as plt

# Set the path (to the folder) of the above mat-file to here:
path = ''

options = proj.projectorClass()

###########################################################################
###########################################################################
###########################################################################
########################### SCANNER PROPERTIES ############################
###########################################################################
###########################################################################
###########################################################################

### Name of current datafile/examination
# This is used for (potential) naming purposes only
options.name = 'Planmeca_CT_data'

### Show status messages
# These are e.g. time elapsed on various functions and what steps have been
# completed. It is recommended to keep this at 1 or 2. With value of 2, 
# you get more detailed timing information. Maximum is 3. Minimum is 0.
options.verbose = 1


###########################################################################
###########################################################################
###########################################################################
###########################################################################

# Load input data
fpath = os.path.join(path, 'Planmeca_VisoG7_100kV_80mAs_500proj_kneePhantom.mat')
from pymatreader import read_mat
var = read_mat(fpath)

# Flat field corrected projections
options.SinM = var['proj']

# Number of projections
options.nProjections = var['nProjections']

# Field of view
options.FOVa_x = var['FOV'][0]
options.FOVa_y = var['FOV'][1]
options.axial_fov = var['FOV'][2]

# Number of rows and columns in a single projection image
options.nRowsD = var['nRowsD']
options.nColsD = var['nColsD']

# Object offset values from the origin, i.e. how much is the origin of the
# FOV shifted
options.oOffsetX = var['oOffset'][0]
options.oOffsetY = var['oOffset'][1]
options.oOffsetZ = var['oOffset'][2]

# Flat value
# If omitted, will use the maximum value from the input
options.flat = var['flatValue']

# Distance to center of rotation and detector
options.sourceToCRot = np.float32(var['sourceToCRot'])
options.sourceToDetector = np.float32(var['sourceToDetector'])

# Detector pixel size
options.dPitchX = np.float32(var['dPitch'][0])
options.dPitchY = np.float32(var['dPitch'][1])

# Projection angles
options.angles = np.float32(var['projAngles'])

# Rotation of the detector panel
options.pitchRoll = np.float32(var['panelRot'])

# Coordinates for the source and center of the detector
# As source-detector pairs
options.x = var['xCoord']
options.y = var['yCoord']
options.z = var['zCoord']

del var


###########################################################################
###########################################################################
###########################################################################
########################### IMAGE PROPERTIES ##############################
###########################################################################
###########################################################################
###########################################################################

### Reconstructed image pixel count (X/row-direction)
options.Nx = 801

### Y/column-direction
options.Ny = 801

### Z-direction (number of slices) (axial)
options.Nz = 668

# Use these two to rotate/flip the final image
### Flip the image (in column direction)?
options.flip_image = True

### How much is the image rotated (radians)?
# The angle (in radians) on how much the image is rotated BEFORE
# reconstruction, i.e. the rotation is performed in the detector space.
options.offangle = (3*np.pi)/2

###########################################################################
###########################################################################
###########################################################################
###########################################################################

### Use projection extrapolation
# If True, extrapolates the projection data. You can select below whether
# this extrapolation is done only in the axial or transaxial directions, or
# both. Default extrapolation length is 20% of the original length, for
# both sides. For example if axial extrapolation is enabled, then the left
# and right regions of the projection get 20% increase in size. This value
# can be adjusted in CTEFOVCorrection. The values are scaled to air with
# the use of logarithmic scaling.
options.useExtrapolation = False

### Use extended FOV
# Similar to above, but expands the FOV. The benefit of expanding the FOV
# this way is to enable to the use of multi-resolution reconstruction or
# computation of the priors/regularization only in the original FOV. The
# default extension is 40% per side
options.useEFOV = False

# Use transaxial extended FOV (this is off by default)
options.transaxialEFOV = False

# Use axial extended FOV (this is on by default. If both this and
# transaxialEFOV are False but useEFOV is True, the axial EFOV will be
# turned on)
options.axialEFOV = True

# Same as above, but for extrapolation. Same default behavior exists.
options.transaxialExtrapolation = True

# Same as above, but for extrapolation. Same default behavior exists.
options.axialExtrapolation = True

# Setting this to True uses multi-resolution reconstruction when using
# extended FOV. Only applies to extended FOV!
options.useMultiResolutionVolumes = True

# This is the scale value for the multi-resolution volumes. The original
# voxel size is divided by this value and then used as the voxel size for
# the multi-resolution volumes. Default is 4 times the original voxel size.
# This means that the multi-resolution regions have larger voxel sizes if
# this is < 1, i.e. 1/4 = 4 times the original voxel size.
options.multiResolutionScale = 1/4

# Performs the extrapolation and adjusts the image size accordingly
CTEFOVCorrection(options)

# Use offset-correction
# If you use offset imaging, i.e. the center of rotation is not in the
# origin but rather a circle around the origin, you can enable automatic
# offset weighting by setting this to True.
options.offsetCorrection = False


###########################################################################
###########################################################################
###########################################################################
######################## RECONSTRUCTION PROPERTIES ########################
###########################################################################
###########################################################################
###########################################################################

############################# IMPLEMENTATIONS #############################
### OpenCL/CUDA device used 
# NOTE: Use 
# from omegatomo.util.devinfo import deviceInfo
# deviceInfo(True)
# to determine the device nubmers
options.deviceNum = 0

### Use CUDA
# Selecting this to True will use CUDA kernels/code instead of OpenCL. This
# only works if the CUDA code was successfully built. This is recommended
# if you have CUDA-capable device.
options.useCUDA = checkCUDA(options.deviceNum)

############################### PROJECTOR #################################
### Type of projector to use for the geometric matrix
# 1 = Improved/accelerated Siddon's algorithm
# 4 = Interpolation-based projector (ray- and voxel-based)
# 5 = Branchless distance-driven projector
# NOTE: You can mix and match most of the projectors. I.e. 45 will use
# interpolation-based projector for forward projection while branchless
# distance-driven is used for backprojection
# NOTE 2: The below additional options apply also in hybrid cases as long
# as the other projector is the corresponding projector.
# See the doc for more information:
# https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
options.projector_type = 4

### Interpolation length (projector type = 4 only)
# This specifies the length after which the interpolation takes place. This
# value will be multiplied by the voxel size which means that 1 means that
# the interpolation length corresponds to a single voxel (transaxial)
# length. Larger values lead to faster computation but at the cost of
# accuracy. Recommended values are between [0.5 1].
options.dL = 1

###########################################################################




###########################################################################
###########################################################################
###########################################################################
######################## RECONSTRUCTION ALGORITHMS ########################
###########################################################################
###########################################################################
###########################################################################
# Reconstruction algorithms to use (choose only one algorithm

### Feldkamp-Davis-Kress (FDK)
options.FDK = True


############################# FDK PROPERTIES #############################
# Use Parker weights
# By default, FDK does NOT use Parker weights. However, if you are
# reconstructing a scan with less than 2pi coverage, it is HIGHLY
# recommended to set the Parker weights true as they are here. Note that
# you can also set an optional Parker weight value which affects the Parker
# weights as described in https://doi.org/10.1118/1.1450132
options.useParkerWeights = True
# Specify the Parker weight below. The default value is 0.25 as in the
# above article, but 1 gives the "original" Parker weights.
options.ParkerWeight = 0.25

###########################################################################
###########################################################################
###########################################################################
###########################################################################

import time
tic = time.perf_counter()
pz, fp = reconstructions_mainCT(options)
toc = time.perf_counter()
print(f"{toc - tic:0.4f} seconds")

z = np.int16(pz[:,:,:] * 55000) - 1000

# plt.imshow(z[:,:,50], vmin=-1000, vmax=2000)
# plt.show()
from omegatomo.util.volume3Dviewer import volume3Dviewer
volume3Dviewer(z, [-1000, 2000])