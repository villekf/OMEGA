# -*- coding: utf-8 -*-
"""
# Python codes for CT reconstruction using Skyscan µCT projection images
This example uses Skyscan µCT projection images. Furthermore, this example
showcases the use of the "scalable" reconstruction for high-dimensional
cases. This means that data of dozens of gigabytes can be successfully
reconstructed on GPUs that cannot store the projection data and the final
reconstructed image volume. The caveat is that functionality is limited.
Only some of the algorithms are supported, such as PDHG, FDK and PKMA.
Furthermore, only filtering-based preconditioner is supported.
Multi-resolution reconstruction is not supported. The data and image are
divided into options.subsets number of segments where only one segment is
present at the GPU at a time. This means that the more subsets you use,
the less memory will be used on the GPU side. The intermediate data will
be stored in host (CPU) so high physical memory amount is recommended.
Only subset types 1 and 8 are supported, though 8 should be used with CT
data. Furthermore, only projector_type = 4 is supported. Only some of the
regularization methods are supported, such as RDP, TV, NLM and GGMRF.

You can use https://doi.org/10.5281/zenodo.12744181 as example data
"""

import numpy as np
from omegatomo import proj
from omegatomo.io import loadSkyscanData
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

# This setting determines whether the high-dimensional scalable
# reconstruction is used (if set as True). Otherwise, the regular
# reconstruction is performed.
# NOTE: Currently the high-dimensional reconstructions are scaled
# differently than the regular ones
options.largeDim = True

### Binning
# The level of binning used for the raw data. For example binning of 2
# reduces the size of the projections by two from both dimensions (e.g.
# 2048x3072 becomes 1024x1536).
options.binning = 1

### Name of current datafile/examination
# This is used for naming purposes only
options.name = 'Skyscan_jyva_data'

### Compute only the reconstructions
# If this file is run with this set to True, then the data load and
# sinogram formation steps are always skipped. Precomputation step is
# only performed if precompute_lor = True and precompute_all = True
# (below). Normalization coefficients are not computed even if selected.
options.only_reconstructions = False

### Show status messages
# These are e.g. time elapsed on various functions and what steps have been
# completed. It is recommended to keep this 1. Maximum value of 3 is supported.
options.verbose = 1

### Transaxial FOV size (mm), this is the length of the x (horizontal) side
# of the FOV
options.FOVa_x = 13

### Transaxial FOV size (mm), this is the length of the y (vertical) side
# of the FOV
options.FOVa_y = options.FOVa_x

### Axial FOV (mm)
options.axial_fov = 7.5


###########################################################################
###########################################################################
###########################################################################
###########################################################################

# If left blank, you will be prompted for the log-file. Alternatively,
# input the full path to the desired log-file here
options.fpath = '/path/to/jyvat.log'


loadSkyscanData(options)
### Projection angles (degree or radian)
# The angles corresponding to the projections
options.angles = -np.arange(0, 0.2*options.nProjections, 0.2, dtype=np.float32)
# NOTE: If you want to reduce the number of projections, you need to do
# this manually as outlined below:
# options.SinM = options.SinM(:,:,1:4:options.nProjections)
# options.angles = options.angles(1:4:numel(options.angles))
# options.nProjections = numel(options.angles)




###########################################################################
###########################################################################
###########################################################################
########################### IMAGE PROPERTIES ##############################
###########################################################################
###########################################################################
###########################################################################

### Reconstructed image pixel size (X-direction)
options.Nx = 1000 * 2

### Y-direction
options.Ny = 1000 * 2

### Z-direction (number of slices) (axial)
options.Nz = 481 * 2

### Flip the image (in vertical direction)?
options.flip_image = False

### How much is the image rotated (radians)?
# The angle (in radians) on how much the image is rotated BEFORE
# reconstruction, i.e. the rotation is performed in the detector space.
options.offangle = (0*np.pi)/2

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

############################# IMPLEMENTATIONS #############################
### Device (GPU) used
# In implementation 2 this determines the device used image reconstruction.
# NOTE: Use the following lines to determine the numbers:
# from omegatomo.util.devinfo import deviceInfo
# deviceInfo(True)
options.use_device = 0

### Use CUDA
# Selecting this to True will use CUDA kernels/code instead of OpenCL. This
# only works if the CUDA code was successfully built.
options.use_CUDA = False


############################### PROJECTOR #################################
### Type of projector to use for the geometric matrix
# 4 = Interpolation-based projector (ray- and voxel-based)
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
# The default mask covers a circle that fills the squre FOV
columns_in_image, rows_in_image = np.meshgrid(np.arange(1, options.Nx + 1), np.arange(1, options.Ny + 1))
centerX = options.Nx / 2
centerY = options.Ny / 2
radius = options.Nx/2
options.maskBP = ((rows_in_image - centerY)**2 + (columns_in_image - centerX)**2 <= radius**2).astype(np.uint8)

######################### RECONSTRUCTION SETTINGS #########################
### Number of iterations (all reconstruction methods except FDK)
# Note that using FDK automatically sets this to 1
options.Niter = 4

### Number of subsets
# Note that with high-dimensional data this is required for FDK as well.
# As mentioned above, for high-dimensional data this controls the amount of
# memory required by the GPU. More subsets, less memory, but using too many
# subsets can lead to reduced performance.
options.subsets = 20

### Subset type (n = subsets)
# 8 = Use every nth projection image
# For high-dimensional data, do not change this!
options.subsetType = 8

### Initial value for the reconstruction
# Should not be used with high-dimensional reconstruction
if not options.largeDim:
    options.x0 = np.ones((options.Nx, options.Ny, options.Nz),dtype=np.float32) * 1e-4

###########################################################################




###########################################################################
###########################################################################
###########################################################################
######################## RECONSTRUCTION ALGORITHMS ########################
###########################################################################
###########################################################################
###########################################################################
# Reconstruction algorithms to use (you can choose one)

### Feldkamp-Davis-Kress (FDK)
options.FDK = True

### Preconditioned Krasnoselskii-Mann algorithm (PKMA)
options.PKMA = False

### Primal-dual hybrid gradient (PDHG)
options.PDHG = False

### Primal-dual hybrid gradient (PDHG) with L1 minimization
options.PDHGL1 = False


################################# PRIORS ##################################
# These priors do not work with FDK!
# Note that currently regularization will not work optimally with
# high-dimensional data. What this means is that you'll get a "flashing
# effect" in certain slices.
### Proximal TV
options.ProxTV = False

### Non-local Means (NLM) prior
options.NLM = False

### Relative difference prior
options.RDP = False

### Generalized Gaussian Markov random field (GGMRF) prior
options.GGMRF = False


######################### REGULARIZATION PARAMETER ########################
### The regularization parameter for ALL regularization methods (priors)
# ~.1 is good starting region for RDP and NLRD
options.beta = .1


############################ ENFORCE POSITIVITY ###########################
### Applies to PDHG, PDHGL1, PKMA
# Enforces positivity in the estimate after each iteration
options.enforcePositivity = True


#################### MEASUREMENT-DOMAIN PRECONDITIONERS ###################
# 1 = Filtering-based preconditioner
# Note: At the moment it is not recommended to use filtering with PKMA, use
# at your own risk!
options.precondTypeMeas[1] = True

### Filtering-steps
# The number of filtering steps for image- and measurement-domain
# filtering. Includes subsets/sub-iterations as well
# E.g. if you have 10 filtering steps, 3 iterations, and 5 subsets, the
# first 2 iterations will use filtering
# There are no restrictions on when to stop the filtering
options.filteringIterations = 200


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
options.tauCP = 0
# Primal value for filtered iterations, applicable only if
# options.precondTypeMeas[2] = True. As with above, automatically computed
# if left zero or empty. Same restrictions apply here as above.
# Use the "Largest eigenvalue for volume 0 with filtering" value here!
options.tauCPFilt = 0
# Dual value. Recommended to set at 1.
options.sigmaCP = 1
# Next estimate update variable, recommended to keep at 1.
options.thetaCP = 1

############################# PKMA PROPERTIES #############################
### Relaxation parameter for PKMA
# If a scalar (or an empty) value is used, then the relaxation parameter is
# computed automatically as lambda(i) = (1 / ((i - 1)/20 + 1)) / 10000,
# where i is the iteration number. The input number thus has no effect.
# If, on the other hand, a vector is input then the input lambda values are
# used as is without any modifications (the length has to be at least the
# number of iterations).
options.lambdaN = 0

# If the reconstruction doesn't work or there are holes or bright spots, 
# then this relaxation value is too high. Reduce the last value (default 
# is 8.1e-7) in such cases and try again.
options.lambdaN = np.zeros(options.Niter, dtype=np.float32)
for i in range(options.Niter):
    if np.sum(options.precondTypeMeas) == 0:
        options.lambdaN[i] = 1 / ((1/3500) * i + 1) / 1 * 8.1e-7
    else:
        options.lambdaN[i] = 1 / ((1/3500) * i + 1) / 1 * 5.1e-3

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
options.NLRD = False

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


###########################################################################
###########################################################################
###########################################################################
########################### OPENCL DEVICE INFO ############################
###########################################################################
###########################################################################
###########################################################################

# Uncomment the below lines and run them to determine the available device
# numbers
# from omegatomo.util.devinfo import deviceInfo
# deviceInfo(True)
###########################################################################
###########################################################################
###########################################################################
###########################################################################



import time
tic = time.perf_counter()
pz = reconstructions_mainCT(options)
toc = time.perf_counter()
print(f"{toc - tic:0.4f} seconds")
plt.pyplot.imshow(pz[:,:,500], vmin=0)