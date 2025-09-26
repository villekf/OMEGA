# -*- coding: utf-8 -*-
"""
Python codes for SPECT reconstruction from DICOM data
"""
import numpy as np
from omegatomo.projector import proj
from omegatomo.reconstruction import reconstructions_mainSPECT
from omegatomo.fileio.loadProSpectaData import loadProSpectaData

# Initialize projector class for reconstruction
options = proj.projectorClass()

options.fpath = ''

###########################################################################
###########################################################################
###########################################################################
########################### SCANNER PROPERTIES ############################
###########################################################################
###########################################################################
###########################################################################

### Crystal thickness (mm)
options.cr_p = 9.525

### Crystal width (mm)
options.dPitchX = 4.7952
options.dPitchY = 4.7952

### Scanner name
# Used for naming purposes (measurement data)
options.machine_name = 'Prospecta'
 
###########################################################################
###########################################################################
###########################################################################
########################### IMAGE PROPERTIES ##############################
###########################################################################
###########################################################################
###########################################################################
 
### Reconstructed image pixel count
# NOTE: Non-square image sizes (X- and Y-direction) may not work
options.Nx = 128; # X-direction
options.Ny = 128; # Y-direction
options.Nz = 128; # Z-direction (number of axial slices)

### FOV size [mm]
# NOTE: Non-cubical voxels may not work
options.FOVa_x = options.dPitchX*128; # [mm], x-axis of FOV (transaxial)
options.FOVa_y = options.dPitchX*128; # [mm], y-axis of FOV (transaxial)
options.axial_fov = options.dPitchY*128; # [mm], z-axis of FOV (axial)

### Flip the image?
options.flipImageX = False
options.flipImageY = False
options.flipImageZ = False

### Use back projection mask?
options.useMaskBP = False
options.maskBP = np.ones((options.Nx, options.Ny, options.Nz))

### Attenuation correction
options.attenuation_correction = False
options.vaimennus = np.zeros((options.Nx, options.Ny, options.Nz))

### How much is the image rotated in degrees?
# NOTE: The rotation is done in the detector space (before reconstruction).
# Positive values perform the rotation in counterclockwise direction
options.offangle = 0
 
###########################################################################
###########################################################################
###########################################################################
###########################################################################
 
loadProSpectaData(options)

###########################################################################
###########################################################################
###########################################################################
########################## COLLIMATOR PROPERTIES ##########################
###########################################################################
###########################################################################
###########################################################################

### Collimator-detector response function (CDRF)
# For projector types 1, 2 and 6 you can either input either:
# 1. the collimator parameters (default) for an analytic solution for round (and hexagonal) holes (this may be unoptimal),
# 2. the standard deviations for both transaxial and axial directions or
# 3. the (Gaussian) PSF filter
# 4. the shifts of each ray traced 

# NOTE: For projector type 1 the CDRF is determined by
# options.rayShiftsDetector and options.rayShiftsSource defined in option 
# 4. These can also be calculated automatically when collimator parameters
# (1.) are input.

# NOTE: With projector_type == 2 (orthogonal distance projector), only the
# collimator parameters below are used for CDR calculation i.e. the
# collimator hole is assumed to be a circle. Thus only 1. below is
# supported with projector_type == 2

# 1. The collimator parameters (projector types 1, 2 and 6)
# Collimator hole length (mm)
options.colL = 24.05
# Collimator hole radius (mm)
options.colR = 1.11/2
# Distance from collimator to the detector (mm)
options.colD = 0
# Intrinsic resolution (mm)
options.iR = 3.8

# 2. If you have the standard deviations for transaxial (XY) and axial (Z)
# directions, you can input them here instead of the above values The
# dimensions need to be options.nProjections x options.Nx. Only for
# projector type 6.
# options.sigmaZ = np.tile(1, options.nProjections, options.Nx)
# options.sigmaXY = np.tile(1, options.nProjections, options.Nx)

# 3. You can input the filter for the CDRF directly. This should be of the
# size filterSizeXY x filterSizeZ x options.nProjections. Only for
# projector type 6.
# options.gFilter = np.ones(1, 1, options.Nx, options.nProjections)

# 4. For the Siddon ray tracer, the CDRF is defined by shifting the rays to
# the shape of the collimator hole. The values of rayShiftsDetector and
# rayShiftsSource should be in the range of [-1, 1]. For example,
# dx1=dy1=-1 below would have the end of the ray shifted to the corner of a
# detector element.
options.nRays = 1; # Number of rays traced per detector element
# options.rayShiftsDetector = np.zeros(2*options.nRays, 1) # The relative shifts (dx1, dy1, dx2, dy2, ...) at the collimator-detector interface
# options.rayShiftsSource = np.zeros(2*options.nRays, 1) # The relative shifts (dx1, dy1, dx2, dy2, ...) at the other end of the collimator

###########################################################################
###########################################################################
###########################################################################
############################# MISC PROPERTIES #############################
###########################################################################
###########################################################################
###########################################################################

### Name of current datafile/examination
# This is used to name the saved measurement data and also load it in
# future sessions.
options.name = 'spect_example'

### Show status messages
# These are e.g. time elapsed on various functions and what steps have been
# completed. It is recommended to keep this 1.  Maximum value of 3 is
# supported.
options.verbose = 1

###########################################################################
###########################################################################
###########################################################################
######################## RECONSTRUCTION PROPERTIES ########################
###########################################################################
###########################################################################
###########################################################################
 
############################# IMPLEMENTATIONS #############################
### OpenCL/CUDA device used 
# NOTE: to obtain the device numbers uncomment the following two lines.
# from omegatomo.util.devinfo import deviceInfo
# deviceInfo(True)
options.deviceNum = 0

### Use CUDA
# Selecting this to True will use CUDA kernels/code instead of OpenCL. This
# only works if the CUDA code was successfully built. 
options.useCUDA = False

### Use CPU
# Selecting this to True will use CPU-based code instead of OpenCL or CUDA.
options.useCPU = False
 
############################### PROJECTOR #################################
### Type of projector to use for the geometric matrix
# 1 = (Improved) Siddon ray-based projector
# 2 = Orthogonal distance ray tracing
# 6 = Rotation-based projector
# See the documentation on some details on the projectors:
# https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
options.projector_type = 1

### Use images instead of buffers? For rotation-based projector this
# implies hardware texture interpolation, which typically has 8 bit 
# precision. With buffers, software interpolation with 32 bit floats is
# used.
options.useImages = True

######################### RECONSTRUCTION SETTINGS #########################
### Number of iterations (all reconstruction methods)
options.Niter = 1
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
options.subsets = 8

### Subset type (n = subsets)
# 8 = Use every nth projection image (recommended for projector_type = 6)
# 9 = Randomly select the projection images
# 10 = Use golden angle sampling to select the subsets (not recommended for PET)
# 11 = Use prime factor sampling to select the projection images
options.subsetType = 8

### Initial value for the reconstruction
options.x0 = np.ones((options.Nx, options.Ny, options.Nz), dtype=np.float32)

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
### Ordered Subsets Expectation Maximization (OSEM) OR Maximum-Likelihood
### Expectation Maximization (MLEM) (if subsets = 1)
# Supported by all implementations
options.OSEM = True

### Modified Row-Action Maximum Likelihood Algorithm (MRAMLA)
# Supported by implementations 2, and 4
options.MRAMLA = False

### Row-Action Maximum Likelihood Algorithm (RAMLA)
# Supported by implementations 2, and 4
options.RAMLA = False

### Relaxed Ordered Subsets Expectation Maximization (ROSEM)
# Supported by implementations 2, and 4
options.ROSEM = False

### LSQR
# Supported by implementations 2, and 4
options.LSQR = False

### Conjugate Gradient Least-squares (CGLS)
# Supported by implementations 2, and 4
options.CGLS = False

### Rescaled Block Iterative Expectation Maximization (RBI-EM)
# Supported by implementations 2, and 4
options.RBI = False

### Dynamic RAMLA (DRAMA)
# Supported by implementations 2, and 4
options.DRAMA = False

### Complete data OSEM (COSEM)
# Supported by implementations 2, and 4
options.COSEM = False

### Enhanced COSEM (ECOSEM)
# Supported by implementations 2, and 4
options.ECOSEM = False

### Accelerated COSEM (ACOSEM)
# Supported by implementations 2, and 4
options.ACOSEM = False
 
 
############################### MAP-METHODS ###############################
# These algorithms can utilize any of the selected priors, though only one
# prior can be used at a time

### One-Step Late OSEM (OSL-OSEM)
# Supported by implementations 2, and 4
options.OSL_OSEM = False

### Modified BSREM (MBSREM)
# Supported by implementations 2, and 4
options.MBSREM = False

### Block Sequential Regularized Expectation Maximization (BSREM)
# Supported by implementations 2, and 4
options.BSREM = False

### ROSEM-MAP
# Supported by implementations 2, and 4
options.ROSEM_MAP = False

### RBI-OSL
# Supported by implementations 2, and 4
options.OSL_RBI = False

### (A)COSEM-OSL
# 0/False = No COSEM-OSL, 1/True = ACOSEM-OSL, 2 = COSEM-OSL
# Supported by implementations 2, and 4
options.OSL_COSEM = False

### Preconditioner Krasnoselskii-Mann algorithm (PKMA)
# Supported by implementations 2, and 4
options.PKMA = False

### Primal-dual hybrid gradient (PDHG)
# Supported by implementations 2, and 4
options.PDHG = False

### Primal-dual hybrid gradient (PDHG) with L1 minimization
# Supported by implementations 2, and 4
options.PDHGL1 = False

### Primal-dual hybrid gradient (PDHG) with Kullback-Leibler minimization
# Supported by implementations 2, and 4
options.PDHGKL = False

### Primal-dual Davis-Yin (PDDY)
# Supported by implementation 2
options.PDDY = False
 
 
################################# PRIORS ##################################
### Median Root Prior (MRP)
options.MRP = False

### Quadratic Prior (QP)
options.quad = False

### Huber Prior (QP)
options.Huber = False

### L-filter prior
options.L = False

### Finite impulse response (FIR) Median Hybrid (FMH) prior
options.FMH = False

### Weighted mean prior
options.weighted_mean = False

### Total Variation (TV) prior
options.TV = False

### Anisotropic Diffusion Median Root Prior (ADMRP)
options.AD = False

### Asymmetric Parallel Level Set (APLS) prior
options.APLS = False

### Hyperbolic prior
options.hyperbolic = False

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
 
 
############################ ACOSEM PROPERTIES ############################
### Acceleration parameter for ACOSEM (1 equals COSEM)
options.h = 2


########################## RELAXATION PARAMETER ###########################
### Relaxation parameter for MRAMLA, RAMLA, ROSEM, BSREM, MBSREM and PKMA
# Use scalar if you want it to decrease as
# lambda / ((current_iteration - 1)/20 + 1). Use vector (length = Niter) if
# you want your own relaxation parameters. Use empty array or zero if you
# want to OMEGA to compute the relaxation parameter using the above formula
# with lamda = 1. Note that current_iteration is one-based, i.e. it starts
# at 1.
options.lambdaN = np.zeros(1, dtype=np.float32)
 

######################## MRAMLA & MBSREM PROPERTIES #######################
### Upper bound for MRAMLA/MBSREM (use 0 for default (computed) value)
options.U = 0
 

############################# PKMA PROPERTIES #############################
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
options.delta_PKMA = 1

 
############################ DRAMA PROPERTIES #############################
### Beta_0 value
options.beta0_drama = 0.1
### Beta value
options.beta_drama = 1
### Alpha value
options.alpha_drama = 0.1

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
# Can lead to unstable behavior with using multi-resolution
# Minimal to none use with filtering-based preconditioner
options.PDAdaptiveType = 0

############################# PRECONDITIONERS #############################
### Applies to PDHG, PDHGL1, PDHGKL, PKMA, MBSREM, MRAMLA, PDDY, FISTA and
### FISTAL1
# Measurement-based preconditioners
# precondTypeMeas(0) = Diagonal normalization preconditioner (1 / (A1))
# precondTypeMeas(1) = Filtering-based preconditioner
options.precondTypeMeas[0] = False
options.precondTypeMeas[1] = False

# Image-based preconditioners
# Setting options.precondTypeImage(1) = true when using PKMA, MRAMLA or
# MBSREM is recommended
# precondTypeImage(0) = Diagonal normalization preconditioner (division with
# the sensitivity image 1 / (A^T1), A is the system matrix) 
# precondTypeImage(1) = EM preconditioner (f / (A^T1), where f is the current
# estimate) 
# precondTypeImage(2) = IEM preconditioner (max(n, fhat, f)/ (A^T1), where
# fhat is an estimate of the final image and n is a small positive number) 
# precondTypeImage(3) = Momentum-like preconditioner (basically a step size
# inclusion) 
# precondTypeImage(4) = Gradient-based preconditioner (Uses the normalized
# divergence (sum of the gradient) of the current estimate) 
# precondTypeImage(5) = Filtering-based preconditioner
# precondTypeImage(6) = Curvature-based preconditioner
options.precondTypeImage[0] = False
options.precondTypeImage[1] = False
options.precondTypeImage[2] = False
options.precondTypeImage[3] = False
options.precondTypeImage[4] = False
options.precondTypeImage[5] = False
options.precondTypeImage[6] = False

# Reference image for precondTypeImage(3). Can be either a mat-file or a
# variable
options.referenceImage = ''

# Momentum parameter for precondTypeImage(4)
# Set the desired momentum parameters to the following variable (note that
# the length should be options.Niter * options.subsets): 
# options.alphaPrecond = []
# Otherwise set the following parameters:
options.rhoPrecond = options.rho_PKMA
options.delta1Precond = options.delta_PKMA

# Parameters for precondTypeImage(5)
# See the article for details
options.gradV1 = 1.5
options.gradV2 = 2
# Note that these include subiterations (options.Niter * options.subsets)
options.gradInitIter = 1
options.gradLastIter = 100

# Number of filtering iterations
# Applies to both precondTypeMeas(2) and precondTypeImage(6)
options.filteringIterations = 100


######################### REGULARIZATION PARAMETER ########################
### The regularization parameter for ALL regularization methods (priors)
options.beta = 1
 
 
######################### NEIGHBORHOOD PROPERTIES #########################
### How many neighboring pixels are considered 
# With MRP, QP, L, FMH, NLM, GGMRF and weighted mean
# E.g. if Ndx = 1, Ndy = 1, Ndz = 0, then you have 3x3 square area where
# the pixels are taken into account (I.e. (Ndx*2+1)x(Ndy*2+1)x(Ndz*2+1)
# area).
# NOTE: Currently Ndx and Ndy must be identical.
# For NLM this is often called the "search window".
options.Ndx = 1
options.Ndy = 1
options.Ndz = 1
 
 
############################## QP PROPERTIES ##############################
### Pixel weights for quadratic prior
# The number of pixels need to be the amount of neighboring pixels,
# e.g. if the above Nd values are all 1, then 27 weights need to be
# included where the center pixel (if Nd values are 1, element 14) should
# be Inf. Size is (Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1). If left empty then
# they will be calculated by the algorithm and are based on the distance of
# the voxels from the center.
options.weights = np.empty(0)
 
 
############################## HP PROPERTIES ##############################
### Delta parameter for Huber prior
# Upper and lower bounds for the prior
options.huber_delta = 5

### Pixel weights for Huber prior
# Same rules apply as with quadratic prior weights.
# If left empty then they will be calculated by the algorithm and are based
# on the distance of the voxels from the center.
options.weights_huber = np.empty(0)
 
 
########################### L-FILTER PROPERTIES ###########################
### Weighting factors for the L-filter pixels
# Otherwise the same as in quadratic prior, but center pixel is not Inf.
# If left empty then they will be calculated by the algorithm such that the
# weights resemble a Laplace distribution.
options.a_L = np.empty(0)

### If the weighting factors are set empty, then this option will determine
# whether the computed weights follow a 1D weighting scheme (True) or 2D 
# (False).
# See the docs for more information:
# https://omega-doc.readthedocs.io/en/latest/algorithms.html#l-filter
options.oneD_weights = False
 
 
############################## FMH PROPERTIES #############################
### Pixel weights for FMH
# The matrix size needs to be [Ndx*2+1, 4] if Nz = 1 or Ndz = 0, or
# [Ndx*2+1, 13] otherwise.
# The center pixel weight should be in the middle of the weight matrix.
# If the sum of each column is > 1, then the weights will be normalized
# such that the sum = 1.
# If left empty then they will be calculated by the algorithm such that the
# weights follow the same pattern as in the original article.
options.fmh_weights = np.empty(0)

### Weighting value for the center pixel
# Default value is 4, which was used in the original article.
# NOTE: This option is ignored if you provide your own weights.
options.fmh_center_weight = 4
 
 
######################### WEIGHTED MEAN PROPERTIES ########################
### Mean type
# 1 = Arithmetic mean, 2 = Harmonic mean, 3 = Geometric mean
options.mean_type = 1

### Pixel weights for weighted mean
# The number of pixels needs to be the amount of neighboring pixels,
# e.g. if the above Ndx/y/z values are all 1, then 27 weights need to be
# included. Size is (Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1). If left empty then
# they will be calculated by the algorithm such that the weights are
# dependent on the distance from the center pixel to the neighboring
# pixels.
options.weighted_weights = np.empty(0)

### Center pixel weight for weighted mean.
# NOTE: This option is ignored if you provide your own weights.
options.weighted_center_weight = 4
 
 
############################## TV PROPERTIES ##############################
### "Smoothing" parameter
# Also used to prevent zero values in square root.
options.TVsmoothing = 1e-5

### Whether to use an anatomical reference/weighting image with the TV
options.TV_use_anatomical = False

### If the TV_use_anatomical value is set to True, specify filename for the
# reference image here (same rules apply as with attenuation correction
# above). Alternatively you can specifiy the variable that holds the
# reference image.
options.TV_reference_image = 'reference_image.mat'

### Three different TV methods are available.
# Value can be 1, 2, 3, 4 or 6.
# Type 3 is not recommended!
# Types 1 and 2 are the same if anatomical prior is not included
# Type 3 uses the same weights as quadratic prior
# Type 4 is the Lange prior, does not support anatomic weighting.
# Type 6 is a weighted TV, does not support anatomic weighting.
# See the docs for more information:
# https://omega-doc.readthedocs.io/en/latest/algorithms.html#tv
options.TVtype = 1

### Weighting parameters for the TV prior. 
# Applicable only if use_anatomical = True. T-value is specific to the used
# TVtype, e.g. for type 1 it is the edge threshold parameter. See the wiki
# for more details:
# https://omega-doc.readthedocs.io/en/latest/algorithms.html#tv
options.T = 0.5

### C is the weight for the original image in type 3 and is ignored with
# other types
options.C = 1

### Tuning parameter for TV and APLS
options.tau = 1e-8

### Tuning parameter for Lange function in SATV (type 4) or weight factor
### for weighted TV (type 6)
# Setting this to 0 gives regular anisotropic TV with type 4
options.SATVPhi = 0.2
 
 
############################# ADMRP PROPERTIES ############################
### Time step variable for AD (implementation 2 only)
options.TimeStepAD = 0.0625

### Conductivity/connectivity for AD (edge threshold)
options.KAD = 2

### Number of iterations for AD filter
# NOTE: This refers to the AD smoothing part, not the actual reconstruction
# phase.
options.NiterAD = 10

### Flux/conduction type for AD filter
# 1 = Exponential
# 2 = Quadratic
options.FluxType = 1

### Diffusion type for AD (implementation 2 only)
# 1 = Gradient
# 2 = Modified curvature
options.DiffusionType = 1
 
 
############################# APLS PROPERTIES #############################
### Scaling parameter (eta)
# See the wiki for details:
# https://omega-doc.readthedocs.io/en/latest/algorithms.html#tv
options.eta = 1e-5

### "Smoothing" parameter (beta)
# Also used to prevent zero values in square root.
options.APLSsmoothing = 1e-5

### Specify filename for the reference image here (same rules apply as with
# attenuation correction above). As before, this can also be a variable
# instead.
# NOTE: For APSL, the reference image is required.
options.APLS_ref_image = 'reference_image.mat'


########################## HYPERBOLIC PROPERTIES ##########################
### Edge weighting factor
options.hyperbolicDelta = 800.
 
 
############################## TGV PROPERTIES #############################
### TGV weights
# First part
options.alpha0TGV = 1
# Second part (symmetrized derivative)
options.alpha1TGV = 2
 
 
############################## NLM PROPERTIES #############################
### Filter parameter
options.sigma = 10

### Patch radius
options.Nlx = 1
options.Nly = 1
options.Nlz = 1

### Standard deviation of the Gaussian filter
options.NLM_gauss = 1

# Search window radius is controlled by Ndx, Ndy and Ndz parameters
# Use anatomical reference image for the patches
options.NLM_use_anatomical = False

### Specify filename for the reference image here (same rules apply as with
# attenuation correction above)
options.NLM_reference_image = 'reference_image.mat'

# Note that only one of the below options for NLM can be selected!
### Use Non-local total variation (NLTV)
# If selected, will overwrite regular NLM regularization as well as the
# below MRP version.
options.NLTV = False

### Use MRP algorithm (without normalization)
# I.e. gradient = im - NLM_filtered(im)
options.NLM_MRP = False

### Use non-local relative difference prior (NLRD)
options.NLRD = False

### Use non-local GGMRF (NLGGMRF)
options.NLGGMRF = False


############################## RDP PROPERTIES #############################
### Edge weighting factor
options.RDP_gamma = 10


############################# GGMRF PROPERTIES ############################
### GGMRF parameters
# See the original article for details
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
# numbers:
# from omegatomo.util.devinfo import deviceInfo
# deviceInfo(True)
###########################################################################
###########################################################################
###########################################################################
###########################################################################
 
## Reconstructions
import time
tic = time.perf_counter()
pz, fp = reconstructions_mainSPECT(options)
toc = time.perf_counter()
print(f"{toc - tic:0.4f} seconds")

# Plot
from matplotlib import pyplot as plt
plt.imshow(pz[:,:,48], vmin=0)
plt.show()
#from omegatomo.util.volume3Dviewer import volume3Dviewer
#volume3Dviewer(pz)