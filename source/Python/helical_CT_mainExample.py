""" Python code for curved helical CT data
This file shows an example reconstruction of curved helical CT data.
Flat panel helical data can be used as CBCT data, but with varying
z-coordinate. This example is thus specifically designed for curved
detectors with cylindrical shape. The example data is from: 
https://aapm.app.box.com/s/eaw4jddb53keg1bptavvvd1sf4x3pe9h/folder/144226105715
Any projection data from above should work, but only L506 was tested
For curved data, setting options.useHelical = True is mandatory
"""
import numpy as np
from omegatomo.projector import proj
from omegatomo.fileio import loadDICOMCTPD
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


### Binning
# The level of binning used for the raw data. For example binning of 2
# reduces the size of the projections by two from both dimensions (e.g.
# 2048x3072 becomes 1024x1536).
# Note: At the moment, binning is not supported with curved helical data
options.binning = 1

### Name of current datafile/examination
# This is used for (potential) naming purposes only
options.name = 'helical_CT_data'

### Compute only the reconstructions
# If this file is run with this set to True, then the data load is always
# skipped 
options.only_reconstructions = False

### Show status messages
# completed. It is recommed to keep this at 1 or 2. With value of 2, 
# you get more detailed timing information. Maximum is 3. Minimum is 0.
options.verbose = 2


###########################################################################
###########################################################################
###########################################################################
###########################################################################

# Put the folder of the DICOM projections here. Only one dataset should be
# in this folder
path = '/path/to/DICOM-CT-PD_FD'

if ~options.only_reconstructions:
    [proj, vars] = loadDICOMCTPD(path)
#else:
    # load presaved data, for example, here


proj[proj < 0] = 0.

# Rotation angles
options.angles = vars['angles']
# The projection images
options.SinM = proj
# Source-detector coordinates for each projection
options.x = np.asfortranarray(np.column_stack((vars['xs'],vars['xd'])))
options.y = np.asfortranarray(np.column_stack((vars['ys'],vars['yd'])))
options.z = np.asfortranarray(np.column_stack((vars['zs'],vars['zd'])))
# This is the radius of the circle formed by the detector array
# Usually the distance of the source to the detector
options.helicalRadius = vars['r']
# Total number of projections
options.nProjections = vars['nProjections']
# Detector pixel sizes
options.dPitchX = vars['dPitchX']
options.dPitchY = vars['dPitchY']
# Source to detector and center of rotation distances
options.sourceToDetector = vars['r']
options.sourceToCRot = vars['sourceToCRot']

# Size of the field-of-view (FOV)
# Transaxial FOV
options.FOVa_x = 500.
# Axial FOV
options.axial_fov = 500.

# Number of rows and columns in the detector array (projection)
options.nRowsD = proj.shape[0]
options.nColsD = proj.shape[1]

# This MUST be set as True when using already linearized data
options.usingLinearizedData = True

del proj

# NOTE: If you want to reduce the number of projections, you need to do
# this manually as outlined below:
# lasku = 1
# options.SinM = options.SinM(:,:,1:lasku:options.nProjections)
# options.angles = options.angles(1:lasku:options.nProjections)
# options.x = options.x(1:lasku:options.nProjections,:)
# options.y = options.y(1:lasku:options.nProjections,:)
# options.z = options.z(1:lasku:options.nProjections,:)
# options.nProjections = numel(options.angles)

# These have to be currently adjusted manually
# The offset of the center of the FOV from the origin of the coordinate
# system
# options.oOffsetX = 0
# options.oOffsetY = 0
options.oOffsetZ = 300.

###########################################################################
###########################################################################
###########################################################################
########################### IMAGE PROPERTIES ##############################
###########################################################################
###########################################################################
###########################################################################

### Reconstructed image voxel count (X/row-direction)
options.Nx = 512

### Y/column-direction
options.Ny = 512

### Z-direction (number of slices) (axial)
options.Nz = 256

# Use these two to rotate/flip the final image
### Flip the image (in column direction)?
options.flip_image = False

### How much is the image rotated (radians)?
# The angle (in radians) on how much the image is rotated BEFORE
# reconstruction, i.e. the rotation is performed in the detector space.
# Positive values perform the rotation in counter-clockwise direction
# Note: Rotation is not yet supported with helical data
options.offangle = (0.*np.pi)/2.

###########################################################################
###########################################################################
###########################################################################
###########################################################################

### Use exted FOV
options.useEFOV = False

# Use axial exted FOV (this is on by default. If both this and
# transaxialEFOV are False but useEFOV is True, the axial EFOV will be
# turned on)
options.axialEFOV = True

# Setting this to True uses multi-resolution reconstruction when using
# exted FOV. Only applies to exted FOV!
options.useMultiResolutionVolumes = True

# This is the scale value for the multi-resolution volumes. The original
# voxel size is divided by this value and then used as the voxel size for
# the multi-resolution volumes. Default is 4 times the original voxel size.
options.multiResolutionScale = 1/4

# Performs the extrapolation and adjusts the image size accordingly
CTEFOVCorrection(options)


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
# only works if the CUDA code was successfully built. This is recommed
# if you have CUDA-capable device.
options.useCUDA = False

### Use CPU
# Selecting this to True will use CPU-based code instead of OpenCL or CUDA.
# Not recommended, even OpenCL with CPU should be used before this.
options.useCPU = False

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
# NOTE 3: Helical CT currently only supports projector_type 4!
options.projector_type = 4

### Use mask
# The mask needs to be a binary mask (uint8 or logical) where 1 means that
# the pixel is included while 0 means it is skipped. Separate masks can be
# used for both forward and backward projection and either one or both can
# be utilized at the same time. E.g. if only backprojection mask is input,
# then only the voxels which have 1 in the mask are reconstructed.
# The mask can be either a 2D image that is applied identically to each slice
# or a 3D mask that is applied as-is
# Forward projection mask
# If nonempty, the mask will be applied. If empty, or completely omitted, no
# mask will be considered.
# options.maskFP = True(options.nRowsD,options.nColsD)
# Backprojection mask
# If nonempty, the mask will be applied. If empty, or completely omitted, no
# mask will be considered.
# Create a circle that fills the FOV:
# [columnsInImage, rowsInImage] = meshgrid(1:options.Nx, 1:options.Ny)
# centerX = options.Nx/2
# centerY = options.Ny/2
# radius = options.Nx/2
# options.maskBP = uint8((rowsInImage - centerY).^2 ...
#     + (columnsInImage - centerX).^2 <= radius.^2)

### Interpolation length (projector type = 4 only)
# This specifies the length after which the interpolation takes place. This
# value will be multiplied by the voxel size which means that 1 means that
# the interpolation length corresponds to a single voxel (transaxial)
# length. Larger values lead to faster computation but at the cost of
# accuracy. Recommed values are between [0.5 1].
options.dL = 1


######################### RECONSTRUCTION SETTINGS #########################
### Number of iterations (all reconstruction methods)
options.Niter = 3
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
# Note that, at the moment, it is recommed to use a large number of
# subsets with helical CT data
options.subsets = 30

### Subset type (n = subsets)
# 8 = Use every nth projection image
# 9 = Randomly select the projection images
# 10 = Use golden angle sampling to select the subsets (not recommed for
# PET)
# 11 = Use prime factor sampling to select the projection images
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

### Preconditioner Krasnoselskii-Mann algorithm (PKMA)
# Supported by implementations 1, 2, 4, and 5
options.PKMA = False

### Primal-dual hybrid gradient (PDHG)
# Supported by implementations 1, 2, 4, and 5
options.PDHG = True

### Primal-dual hybrid gradient (PDHG) with L1 minimization
# Supported by implementations 1, 2, 4, and 5
options.PDHGL1 = False

### Primal-dual Davis-Yin (PDDY)
# Supported by implementation 2
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
options.tauCP = 0
# Primal value for filtered iterations, applicable only if
# options.precondTypeMeas[2] = True. As with above, automatically computed
# if left zero or empty. Same restrictions apply here as above.
# Use the "Largest eigenvalue for volume 0 with filtering" value here!
# if you have a multi-resolution situation, you should input the values
# for each volume or use zero/empty
options.tauCPFilt = 0
# Dual value. Recommed to set at 1.
options.sigmaCP = 1
# Next estimate update variable
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
# precondTypeMeas(1) = Diagonal normalization preconditioner (1 / (A1))
# precondTypeMeas(2) = Filtering-based preconditioner
if options.PDHG or options.PDHGL1 or options.PDDY:
    options.precondTypeMeas[1] = True


# Number of filtering iterations
# Applies to both precondTypeMeas(2) and precondTypeImage(6)
options.filteringIterations = 30

# Image-based preconditioners
# Setting options.precondTypeImage(2) = True when using PKMA, MRAMLA or
# MBSREM is recommed
# precondTypeImage(1) = Diagonal normalization preconditioner (division with
# the sensitivity image 1 / (A^T1), A is the system matrix) 
# precondTypeImage(2) = EM preconditioner (f / (A^T1), where f is the current
# estimate) 
# precondTypeImage(3) = IEM preconditioner (max(n, fhat, f)/ (A^T1), where
# fhat is an estimate of the final image and n is a small positive number) 
# precondTypeImage(4) = Momentum-like preconditioner (basically a step size
# inclusion) 
# precondTypeImage(5) = Gradient-based preconditioner (Uses the normalized
# divergence (sum of the gradient) of the current estimate) 
# precondTypeImage(6) = Filtering-based preconditioner
# precondTypeImage(7) = Curvature-based preconditioner
if options.PKMA:
    options.precondTypeImage[1] = True



############################# PKMA PROPERTIES #############################
### Relaxation parameter for PKMA
# If a scalar (or an empty) value is used, then the relaxation parameter is
# computed automatically as lambda(i) = (1 / ((i - 1)/20 + 1)) / 10000,
# where i is the iteration number. The input number thus has no effect.
# If, on the other hand, a vector is input then the input lambda values are
# used as is without any modifications (the length has to be at least the
# number of iterations).
options.lambdaN = np.zeros(0, dtype=np.float32)

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
options.alpha_PKMA = 0

### rho_PKMA
# This value is ignored if a vector input is used with alpha_PKMA
options.rho_PKMA = 0.95

### delta_PKMA
# This value is ignored if a vector input is used with alpha_PKMA
options.delta_PKMA = 100


######################### REGULARIZATION PARAMETER ########################
### The regularization parameter for ALL regularization methods (priors)
# 50-100 is good starting point for NLM
# ~1 is good starting region for RDP and NLRD
options.beta = 50


######################### NEIGHBORHOOD PROPERTIES #########################
### How many neighboring pixels are considered
# With MRP, QP, L, FMH, NLM and weighted mean
# E.g. if Ndx = 1, Ndy = 1, Ndz = 0, then you have 3x3 square area where
# the pixels are taken into account (I.e. (Ndx*2+1)x(Ndy*2+1)x(Ndz*2+1)
# area).
# NOTE: Currently Ndx and Ndy must be identical.
# For NLM this is often called the "search window".
options.Ndx = 2
options.Ndy = 2
options.Ndz = 1


############################## TGV PROPERTIES #############################
### TGV weights
# First part
options.alpha0TGV = 0.3
# Second part (symmetrized derivative)
options.alpha1TGV = 0.05


############################## NLM PROPERTIES #############################
### Filter parameter
# Higher values smooth the image, smaller values make it sharper
options.sigma = 1.0e-2

### Patch radius
options.Nlx = 1
options.Nly = 1
options.Nlz = 1

### Standard deviation of the Gaussian filter
options.NLM_gauss = 2

### Adaptive NL methods
options.NLAdaptive = False

### Summed constant for adaptive NL
options.NLAdaptiveConstant = 2.0e-5

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
options.SATVPhi = 10

### Use Non-local GGMRF (NLGGMRF)
options.NLGGMRF = False

### Use MRP algorithm (without normalization)
# I.e. gradient = im - NLM_filtered(im)
options.NLM_MRP = False


############################## RDP PROPERTIES #############################
### Edge weighting factor
# Note that this affects NLRD as well
# Higher values sharpen the image, lower values smooth it
options.RDP_gamma = 1

############################# GGMRF PROPERTIES ############################
### GGMRF parameters
# See the original article for details
# These affect the NLGGMRF as well
options.GGMRF_p = 2
options.GGMRF_q = 1.15
options.GGMRF_c = .0005

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
########################### OPENCL DEVICE INFO ############################
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

# Use helical CT
# This must be set to True to use curved helical CT data
options.useHelical = True


import time
tic = time.perf_counter()
pz, fp = reconstructions_mainCT(options)
toc = time.perf_counter()
print(f"{toc - tic:0.4f} seconds")

plt.pyplot.imshow(pz[:,:,150])