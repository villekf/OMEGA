# -*- coding: utf-8 -*-
"""
Python codes for CBCT reconstruction
This example file computes the Figure 4 (a) of the OMEGA V2 article. 
DOI will be added later.
You need to provide the folder path to the input data
Planmeca_VisoG7_100kV_80mAs_500proj_kneePhantom.mat!
Used data available from: https://doi.org/10.5281/zenodo.12722386
"""
import numpy as np
from omegatomo.projector import proj
from omegatomo.util import CTEFOVCorrection
from omegatomo.reconstruction import reconstructions_mainCT
from omegatomo.util import checkCUDA
import matplotlib as plt

# Set the path (to the folder) of the above mat-file to here:
path = ''

options = proj.projectorClass()

# Set this to true, if your GPU has less than 12 GB of memory
largeDim = False

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
options.verbose = 2


###########################################################################
###########################################################################
###########################################################################
###########################################################################

# Load input data
fpath = path + '\Planmeca_VisoG7_100kV_80mAs_500proj_kneePhantom.mat'
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
options.useEFOV = True

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


######################### RECONSTRUCTION SETTINGS #########################
### Number of iterations (all reconstruction methods)
options.Niter = 8
### Save specific intermediate iterations
# You can specify the intermediate iterations you wish to save here. Note
# that this uses zero-based indexing, i.e. 0 is the first iteration (not
# the initial value). By default only the last iteration is saved. Only
# full iterations (epochs) can be saved.
options.saveNIter = np.empty(0)
# Alternatively you can save ALL intermediate iterations by setting the
# below to True and uncommenting it. As above, only full iterations
# (epochs) are saved.
# options.save_iter = False

### Number of subsets
options.subsets = 20

### Subset type (n = subsets)
# 1 = Every nth (column) measurement is taken
# 2 = Every nth (row) measurement is taken (e.g. if subsets = 3, then
# first subset has measurements 1, 4, 7, etc., second 2, 5, 8, etc.) 
# 3 = Measurements are selected randomly
# 4 = (Sinogram only) Take every nth column in the sinogram
# 5 = (Sinogram only) Take every nth row in the sinogram
# 6 = Sort the LORs according to their angle with positive X-axis, combine
# n_angles together and have 180/n_angles subsets for 2D slices and
# 360/n_angles for 3D, see docs for more information:
# https://omega-doc.readthedocs.io/en/latest/algorithms.html#type-6
# 7 = Form the subsets by using golden angle sampling
# 8 = Use every nth sinogram
# 9 = Randomly select the full sinograms
# 10 = Use golden angle sampling to select the subsets (not recommended for
# PET)
# 11 = Use prime factor sampling to select the full sinograms
# Most of the time subsetType 8 is sufficient.
options.subsetType = 8

### Initial value for the reconstruction
options.x0 = np.ones((options.Nx, options.Ny, options.Nz), dtype=np.float32) * 1e-5

###########################################################################




###########################################################################
###########################################################################
###########################################################################
######################## RECONSTRUCTION ALGORITHMS ########################
###########################################################################
###########################################################################
###########################################################################
# Reconstruction algorithms to use (choose only one algorithm and
# optionally one prior)

############################### MAP-METHODS ###############################
# These algorithms can utilize any of the selected priors, though only one
# prior can be used at a time
### Primal-dual hybrid gradient (PDHG)
options.PDHG = True;


############################ ENFORCE POSITIVITY ###########################
### Applies to PDHG, PDHGL1, PDDY, FISTA, FISTAL1, MBSREM, MRAMLA, PKMA
# Enforces positivity in the estimate after each iteration
options.enforcePositivity = True


############################# PDHG PROPERTIES #############################
# Primal value
# If left zero, or empty, it will be automatically computed
# Note that if you change any of the model parameters, i.e. image volume
# size, number of projections or use binning, this needs to be recomputed
# or scaled accordingly!
# The computed largest eigenvalue is printed if verbose > 0. This can be 
# used as the below value as long as one is divided by it. For example, 
# if "Largest eigenvalue for volume 0 is 100" then options.tauCP should be 
# 1/100 (if you use filtering-based preconditioner this is the "without 
# filtering" value)
# if you have a multi-resolution situation, you should input the values
# for each volume or use zero/empty
options.tauCP = np.array([1/7540.685059,1/117180.726562,1/92894.101562],dtype=np.float32)
# Primal value for filtered iterations, applicable only if
# options.precondTypeMeas[1] = True. As with above, automatically computed
# if left zero or empty. Same restrictions apply here as above.
# Use the "Largest eigenvalue for volume 0 with filtering" value here!
# if you have a multi-resolution situation, you should input the values
# for each volume or use zero/empty
options.tauCPFilt = np.array([1/127.564476,1/1606.458252,1/1479.921875],dtype=np.float32)
# Dual value. Recommended to set at 1.
options.sigmaCP = 1
# Next estimate update variable, recommended to keep at 1.
options.thetaCP = 1

# Use adaptive update of the primal and dual variables
# Currently two methods available
# Setting this to 1 or 2 uses an adaptive update for both the primal and 
# dual variables.
# Can lead to unstable behavior when using with multi-resolution
# Minimal to none use with filtering-based preconditioner
options.PDAdaptiveType = 0


############################# PRECONDITIONERS #############################
### Applies to PDHG, PDHGL1, PDHGKL, PKMA, MBSREM, MRAMLA, PDDY, FISTA and
### FISTAL1
# Measurement-based preconditioners
# Default values are False
# precondTypeMeas(0) = Diagonal normalization preconditioner (1 / (A1))
# precondTypeMeas(1) = Filtering-based preconditioner
options.precondTypeMeas[1] = True

# Number of filtering iterations
# Applies to both precondTypeMeas(1) and precondTypeImage(5)
options.filteringIterations = 80

###########################################################################
###########################################################################
###########################################################################
###########################################################################

options.largeDim = largeDim

import time
tic = time.perf_counter()
pz, fp = reconstructions_mainCT(options)
toc = time.perf_counter()
print(f"{toc - tic:0.4f} seconds")

z = np.int16(pz[:,:,:] * 55000) - 1000

plt.pyplot.imshow(z[:,:,50], vmin=-1000, vmax=2000)
