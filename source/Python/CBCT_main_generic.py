# -*- coding: utf-8 -*-
"""
# Python codes for CBCT reconstruction for the Planmeca data
This example contains a simplified example for Planmeca CBCT data.
Simplified refers to less number of selectable options. You can always
add other options manually. This example uses PDHG with subsets and
filtering. Other algorithms that have been left here are PKMA and PDDY.
Regularization includes TV, NLM, RDP, GGMRF and the non-local variations.
You can manually add any of the built-in priors if desired.

This example also showcases CBCT cases where you have the source and center
of detector panel coordinates, as well as additional rotation of the panel.
Instead of the rotation of the panel, you can also input the direction vectors.

Example data available from: https://doi.org/10.5281/zenodo.12722386
"""
import numpy as np
from omegatomo import proj
from omegatomo.util import CTEFOVCorrection
from omegatomo.reconstruction import reconstructions_mainCT
import matplotlib as plt

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

### Compute only the reconstructions
# If this file is run with this set to True, then the data load and
# sinogram formation steps are always skipped. Precomputation step is
# only performed if precompute_lor = True and precompute_all = True
# (below). Normalization coefficients are not computed even if selected.
options.only_reconstructions = False

### Show status messages
# These are e.g. time elapsed on various functions and what steps have been
# completed. It is recommended to keep this 1.
options.verbose = 1


###########################################################################
###########################################################################
###########################################################################
###########################################################################

# Load input data
fpath = 'Planmeca_VisoG7_100kV_80mAs_500proj_kneePhantom.mat'
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
options.sourceToCRot = var['sourceToCRot']
options.sourceToDetector = var['sourceToDetector']

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

# NOTE: If you want to reduce the number of projections, you need to do
# this manually as outlined below:
# lasku = 1
# options.SinM = options.SinM[:,:,::lasku]
# options.angles = options.angles[::lasku]
# options.pitchRoll = options.pitchRoll[::lasku,:]
# options.x = options.x[::lasku,:]
# options.y = options.y[::lasku,:]
# options.z = options.z[::lasku,:]
# options.nProjections = options.angles.size



###########################################################################
###########################################################################
###########################################################################
########################### IMAGE PROPERTIES ##############################
###########################################################################
###########################################################################
###########################################################################
### Reconstructed image pixel size (X-direction)
options.Nx = 801

### Y-direction
options.Ny = 801

### Z-direction (number of slices) (axial)
options.Nz = 668

# Use these two to rotate/flip the final image
### Flip the image (in vertical direction)?
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
options.useExtrapolation = False

### Use extended FOV
options.useEFOV = False

# Use transaxial extended FOV (this is off by default)
options.transaxialEFOV = True

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
# the multi-resolution volumes. Default is 1/4 of the original voxel size.
options.multiResolutionScale = 1/4

# Performs the extrapolation and adjusts the image size accordingly
CTEFOVCorrection(options)

# Use offset-correction
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
# only works if the CUDA code was successfully built. Recommended only for
# Siddon as the orthogonal/volume-based ray tracer are slower in CUDA.
options.useCUDA = False

### Use CPU
# Selecting this to True will use CPU-based code instead of OpenCL or CUDA.
options.useCPU = False

############################### PROJECTOR #################################
### Type of projector to use for the geometric matrix
# 1 = Improved/accelerated Siddon's algorithm
# 2 = Orthogonal distance based ray tracer
# 3 = Volume of intersection based ray tracer
# 4 = Interpolation-based projector (ray- and voxel-based)
# 5 = Branchless distance-driven projector
# NOTE: You can mix and match most of the projectors. I.e. 41 will use
# interpolation-based projector for forward projection while improved
# Siddon is used for backprojection.
# See the doc for more information:
# https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
options.projector_type = 4

### Use mask
# The mask needs to be a binary mask (uint8 or logical) where 1 means that
# the pixel is included while 0 means it is skipped. Separate masks can be
# used for both forward and backward projection and either one or both can
# be utilized at the same time. E.g. if only backprojection mask is input,
# then only the voxels which have 1 in the mask are reconstructed.
# Currently the masks need to be a 2D image that is applied identically at
# each slice.
# Forward projection mask
# If nonempty, the mask will be applied. If empty, or completely omitted, no
# mask will be considered.
# options.maskFP = True(options.nRowsD,options.nColsD)
# Backprojection mask
# If nonempty, the mask will be applied. If empty, or completely omitted, no
# mask will be considered.
# Create a circle that fills the FOV:
# columns_in_image, rows_in_image = np.meshgrid(np.arange(1, options.Nx + 1), np.arange(1, options.Ny + 1))
# centerX = options.Nx / 2
# centerY = options.Ny / 2
# radius = options.Nx / 2
# options.maskBP = ((rows_in_image - centerY)**2 + (columns_in_image - centerX)**2 <= radius**2).astype(np.uint8)

### Interpolation length (projector type = 4 only)
# This specifies the length after which the interpolation takes place. This
# value will be multiplied by the voxel size which means that 1 means that
# the interpolation length corresponds to a single voxel (transaxial)
# length. Larger values lead to faster computation but at the cost of
# accuracy. Recommended values are between [0.5 1].
options.dL = 0.5


######################### RECONSTRUCTION SETTINGS #########################
### Number of iterations (all reconstruction methods)
options.Niter = 10
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

### Number of subsets (all excluding MLEM and subset_type = 5)
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
# 360/n_angles for 3D, see GitHub wiki for more information:
# https://github.com/villekf/OMEGA/wiki/Function-help#reconstruction-settings
# 7 = Form the subsets by using golden angle sampling
# 8 = Use every nth sinogram
# 9 = Randomly select the full sinograms
# 10 = Use golden angle sampling to select the subsets (not recommended for
# PET)
# 11 = Use prime factor sampling to select the full sinograms
# Most of the time subset_type 8 is sufficient.
options.subsetType = 8

### Initial value for the reconstruction
options.x0 = np.ones((options.Nx, options.Ny, options.Nz), dtype=np.float32) * 1e-4

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

############################# NON-REGULARIZED #############################

### LSQR
options.LSQR = False

### Conjugate Gradient Least-squares (CGLS)
options.CGLS = False

### Feldkamp-Davis-Kress (FDK)
options.FDK = False


############################### MAP-METHODS ###############################
# These algorithms can utilize any of the selected priors, though only one
# prior can be used at a time

### PKMA
options.PKMA = False

### Primal-dual hybrid gradient (PDHG)
options.PDHG = True

### Primal-dual hybrid gradient (PDHG) with L1 minimization
options.PDHGL1 = False

### Primal-dual Davis-Yin (PDDY)
options.PDDY = False


################################# PRIORS ##################################
### Total Generalized Variation (TGV) prior
options.TGV = False

### Proximal TV
options.ProxTV = False

### Non-local Means (NLM) prior
options.NLM = False

### Relative difference prior
options.RDP = False

### Generalized Gaussian Markov random field (GGMRF) prior
options.GGMRF = False


############################ ENFORCE POSITIVITY ###########################
### Applies to PDHG, PDHGL1, PDDY, FISTA, FISTAL1, MBSREM, MRAMLA, PKMA
# Enforces positivity in the estimate after each iteration
options.enforcePositivity = True


############################# PDHG PROPERTIES #############################
# Primal value
# If left zero, or empty, it will be automatically computed
options.tauCP = 0
# Primal value for filtered iterations, applicable only if
# options.precondTypeMeas[2] = True. As with above, automatically computed
# if left zero or empty.
options.tauCPFilt = 0
# Dual value. Recommended to set at 1.
options.sigmaCP = 1
# Next estimate update variable
options.thetaCP = 1

# Use adaptive update of the primal and dual variables
# Currently only one method available
# Setting this to 1 uses an adaptive update for both the primal and dual
# variables.
# Can lead to unstable behavior when using with multi-resolution!
# Minimal to none benefit with filtering-based preconditioner
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
options.filteringIterations = 100
 

############################# PKMA PROPERTIES #############################
### Relaxation parameter for PKMA
# If a scalar (or an empty) value is used, then the relaxation parameter is
# computed automatically as lambda(i) = (1 / ((i - 1)/20 + 1)) / 10000,
# where i is the iteration number. The input number thus has no effect.
# If, on the other hand, a vector is input then the input lambda values are
# used as is without any modifications (the length has to be at least the
# number of iterations).
options.lambdaN = np.zeros(1, dtype=np.float32)

### Step size (alpha) parameter for PKMA
# If a scalar (or an empty) value is used, then the alpha parameter is
# computed automatically as alpha_PKMA(oo) = 1 + (options.rho_PKMA *((i -
# 1) * options.subsets + ll)) / ((i - 1) * options.subsets + ll +
# options.delta_PKMA), where i is the iteration number and l the subset
# number. The input number thus has no effect. options.rho_PKMA and
# options.delta_PKMA are defined below.
# If, on the other hand, a vector is input then the input alpha values are
# used as is without any modifications (the length has to be at least the
# number of iterations * number of subsets).
options.alpha_PKMA = np.zeros(1, dtype=np.float32)

### rho_PKMA
# This value is ignored if a vector input is used with alpha_PKMA
options.rho_PKMA = 0.95

### delta_PKMA
# This value is ignored if a vector input is used with alpha_PKMA
options.delta_PKMA = 100.


######################### REGULARIZATION PARAMETER ########################
### The regularization parameter for ALL regularization methods (priors)
# 50-100 is good starting point for NLM
# !1 is good starting region for RDP and NLRD
options.beta = 1


######################### NEIGHBORHOOD PROPERTIES #########################
### How many neighboring pixels are considered
# With MRP, QP, L, FMH, NLM and weighted mean
# E.g. if Ndx = 1, Ndy = 1, Ndz = 0, then you have 3x3 square area where
# the pixels are taken into account (I.e. (Ndx*2+1)x(Ndy*2+1)x(Ndz*2+1)
# area).
# NOTE: Currently Ndx and Ndy must be identical.
# For NLM this is often called the "search window".
options.Ndx = 1
options.Ndy = 1
options.Ndz = 1


############################## TGV PROPERTIES #############################
### TGV weights
# First part
options.alpha0TGV = 1
# Second part (symmetrized derivative)
options.alpha1TGV = 2


############################## NLM PROPERTIES #############################
### Filter parameter
# Higher values smooth the image, smaller values make it sharper
options.sigma = 6.00e-3

### Patch radius
options.Nlx = 1
options.Nly = 1
options.Nlz = 1

### Standard deviation of the Gaussian filter
options.NLM_gauss = 0.75

# By default, the original NLM is used. You can, however, use another
# potential function by selecting ONE of the options below.
### Use Non-local total variation (NLTV)
# If selected, will overwrite regular NLM regularization as well as the
# below MRP version.
options.NLTV = False

### Use Non-local relative difference (NLRD)
options.NLRD = True

### Use Non-local Lange prior (NLLange)
options.NLLange = False

# Tuning parameter for Lange
options.SATVPhi = 5

### Use Non-local GGMRF (NLGGMRF)
options.NLGGMRF = False

### Use MRP algorithm (without normalization)
# I.e. gradient = im - NLM_filtered(im)
options.NLM_MRP = False


############################## RDP PROPERTIES #############################
### Edge weighting factor
# Note that this affects NLRD as well
options.RDP_gamma = 10


############################# GGMRF PROPERTIES ############################
### GGMRF parameters
# See the original article for details
# These affect the NLGGMRF as well
options.GGMRF_p = 1.5
options.GGMRF_q = 1
options.GGMRF_c = 5

###########################################################################
###########################################################################
###########################################################################
###########################################################################


# Store the intermediate forward projections. Unlike image estimates, this
# also stores subiteration results.
options.storeFP = False


###########################################################################
###########################################################################
###########################################################################
############################### DEVICE INFO ###############################
###########################################################################
###########################################################################
###########################################################################

# Uncomment the below lines and run it to determine the available device
# numbers:
# from omegatomo.util.devinfo import deviceInfo
# deviceInfo(True)

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

plt.pyplot.imshow(pz[:,:,120])