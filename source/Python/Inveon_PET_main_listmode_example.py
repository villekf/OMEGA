# -*- coding: utf-8 -*-
"""
## Python codes for PET reconstruction using Inveon PET LST or GATE output
# This main-file provides an example on how to obtain & perform list-mode
# (event-by-event) reconstruction. 
# For the input measurement data, you can use the open preclinical PET data
# available from: https://doi.org/10.5281/zenodo.3528056
"""
import numpy as np
from omegatomo.projector import proj
from omegatomo.reconstruction import reconstructions_main
from omegatomo.fileio.loadInveon import loadInveonData

options = proj.projectorClass()


###########################################################################
###########################################################################
###########################################################################
########################### SCANNER PROPERTIES ############################
###########################################################################
###########################################################################
###########################################################################

### R-sectors/modules/blocks/buckets in transaxial direction
options.blocks_per_ring = (16)

### R-sectors/modules/blocks/buckets in axial direction (i.e. number of physical
### scanner/crystal rings) 
# Multiplying this with the below cryst_per_block should equal the total
# number of crystal rings. 
options.linear_multip = (4)

### Number of detectors on the side of R-sector/block/module (transaxial
### direction)
# (e.g. 13 if 13x13, 20 if 20x10)
options.cryst_per_block = (20)

### Number of detectors on the side of R-sector/block/module (axial
### direction)
# (e.g. 13 if 13x13, 10 if 20x10)
options.cryst_per_block_axial = 20

### Crystal pitch/size in x- and y-directions (transaxial) (mm)
options.cr_p = 1.59

### Crystal pitch/size in z-direction (axial) (mm)
options.cr_pz = 1.59

### Ring diameter (distance between perpendicular detectors) (mm)
options.diameter = 161

# Note that non-square transaxial FOV sizes should work, but might not work
# always. Square transaxial FOV is thus recommended.
### Transaxial FOV size (mm), this is the length of the x (vertical) side
# of the FOV
options.FOVa_x = 100

### Transaxial FOV size (mm), this is the length of the y (horizontal) side
# of the FOV
options.FOVa_y = options.FOVa_x

# The above recommendation doesn't apply to axial FOV, i.e. this can be
# different from the transaxial FOV size(s). 
### Axial FOV (mm)
options.axial_fov = 127

### Number of pseudo rings between physical rings (use 0 or np.empty(0, dtype=np.float32) if none)
# NOTE: Inveon has no pseudo detectors/rings
options.pseudot = np.empty(0, dtype=np.float32)

### Number of gaps between rings
# These should be the rings after which there is a gap
# For example, if there are a total of 6 rings with two gaps and the gaps
# are after ring number 2 and 4. Then options.ringGaps = [2,4]
# Note that the below options.rings should not include the gaps, though by
# default it should be correctly computed.
options.ringGaps = np.empty(0, dtype=np.float32)

### Number of detectors per ring (without pseudo detectors)
options.det_per_ring = options.blocks_per_ring*options.cryst_per_block

### Number of detectors per ring (with pseudo detectors)
# NOTE: Inveon has no pseudo detectors/rings
# If you have a single pseudo detector per block, use the commented line
options.det_w_pseudo = options.blocks_per_ring*(options.cryst_per_block)
# options.det_w_pseudo = options.blocks_per_ring*(options.cryst_per_block + 1)

### Number of crystal rings
options.rings = options.linear_multip * options.cryst_per_block

### Number of detectors
options.detectors = options.det_per_ring*options.rings

### Scanner name
# Used for naming purposes (measurement data)
options.machine_name = 'Inveon'

###########################################################################




###########################################################################
###########################################################################
######################## GATE SPECIFIC SETTINGS ###########################
###########################################################################
###########################################################################

### Obtain Trues (True coincidences)
# If this is set to True then, in addition to the normal coincidences
# (prompts), Trues are also obtained and saved.
options.obtain_Trues = True

### Reconstruct the True coincidences
# If this is set to True, then the True coincidences will be used for
# reconstruction
# NOTE: If both this and reconstruct_scatter are set, then the Trues are
# reconstructed, but not the scatter.
options.reconstruct_Trues = False

### Obtain scattered coincidences
# If this is set to True, then scattered coincidences are saved separately.
# These events are not used for scatter correction though, but a separate
# scatter sinogram/raw data matrix will be created. The scatter elements
# included can be selected below.
options.store_scatter = True

### What scatter components are included in the scatter part
# (1 means that component is included, 0 means it is not included in the
# scatter data) 
# First: Compton scattering in the phantom, second: Compton scattering in
# the detector, third: Rayleigh scattering in the phantom, fourth: Rayleigh
# scattering in the detector.
# If store_scatter is set to True, at least one value has to be 1. E.g. [1
# 0 1 0] will save Compton and Rayleigh scattering in the phantom. 
# NOTE: LMF will always include only Compton scattering in the phantom,
# regardless of the choice below (as long as scatter is selected).
options.scatter_components = np.array([True, True, False, False])

### Reconstruct the scattered coincidences
# If this is set to True, then the scattered coincidences will be used for
# reconstruction
# NOTE: If both this and reconstruct_Trues are set, then the Trues are
# reconstructed, but not the scatter.
options.reconstruct_scatter = False

### Obtain (True) random coincidences.
# If this is set to True then coincidence events that are genuine random
# events are stored separately.
# These events are not used for randoms correction (see the
# Corrections-section for delayed coincidence window randoms correction),
# but a separate randoms sinogram/raw data matrix will be created.
options.store_randoms = True

### Obtain source coordinates (used in forming the "True" image).
# If this is set to True, then the "True" decay image is also saved during
# data load, i.e. the locations where the decay has occurred.
# If any of the above settings are set to True, then the True images are
# also obtained for them. E.g. if store_scatter = True, then an image
# showing the locations and number of counts of where the scattered events
# originated will be saved in a mat-file. Scatter and Trues contain
# coincidence events while randoms contain singles.
# NOTE: If you use LMF data, the source images are not considered reliable.
options.source = False
 
###########################################################################
 
 
 
 
###########################################################################
###########################################################################
###########################################################################
######################## ROOT DATA FORMAT SETTINGS ########################
###########################################################################
###########################################################################
###########################################################################
 
### Is ROOT data used
options.use_root = False
 
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
 
# Note that non-square transaxial image sizes can be unreliable just as the
# non-square transaxial FOV, but they should, generally, work
### Reconstructed image pixel count (X/row-direction)
options.Nx = 128

### Y/column-direction
options.Ny = 128

### Z-direction (number of slices) (axial)
options.Nz = options.rings*2 - 1

### Flip the image (in horizontal direction)?
options.flip_image = False

### How much is the image rotated?
# You need to run the precompute phase again if you modify this
# NOTE: The rotation is done in the detector space (before reconstruction).
# This current setting is for scanner list-mode data or sinogram data.
# Positive values perform the rotation in clockwise direction
options.offangle = options.det_w_pseudo * (2/4)
 
###########################################################################
###########################################################################
###########################################################################
###########################################################################
 
 
 
 
###########################################################################
###########################################################################
###########################################################################
############################# CORRECTIONS #################################
###########################################################################
###########################################################################
###########################################################################
# List-mode reconstruction only supports automatic attenuation correction.
# Randoms and/or scatter has to be done manually. Normalization is not
# supported.

######################### Attenuation correction ##########################
### Image-based attenuation correction
# Include attenuation correction from images (e.g. CT-images) (for this you
# need attenuation images of each slice correctly rotated and scaled for
# 511 keV). For CT-images you can use attenuationCT_to_511 or
# create_atten_matrix_CT functions. If the below attenuation_datafile is
# empty, UMAPs created with Inveon AW can be used (if the below
# CT_attenuation is set to True) or, alternatively, use the produces .atn
# files instead (CT_attenuation is False). 
# E.g. if this is set to True, CT_attenuation = True and
# attenuation_datafile = '', then the user will be automatically prompted
# for the UMAP-files and they will be automatically saved as a mat-file
# with the filename saved in the attenuation_datafile field. Alternatively,
# if this is set to True, CT_attenuation = False and  attenuation_datafile
# = '', then the user will be automatically prompted for the .atn-files
# from which the attenuation images will be automatically created and
# saved (filename is saved in the attenuation_datafile field).
options.attenuation_correction = True

### CT-image attenuation
# Use CT-images (UMAP-image) for the attenuation. If set to False, uses the
# .atn-files instead (if above attenuation is set to True). 
options.CT_attenuation = False

### Attenuation image data file
# Specify the path (if not in MATLAB path) and filename.
# NOTE: the attenuation data must be the only variable in the file and
# have the dimensions of the final reconstructed image. Previously
# saved attenuation images can be used here.
options.attenuation_datafile = ''


############################ Global corrections ###########################
### Global correction factor
# This correction factor will be applied (if nonzero) to all LORs equally.
# This can be e.g. dead time correction factor.
options.global_correction_factor = np.empty(0, dtype=np.float32)
 

#################### Corrections during reconstruction ####################
# If set to True, all the corrections are performed during the
# reconstruction step, otherwise the corrections are performed to the
# sinogram/raw data before reconstruction. I.e. this can be considered as
# e.g. normalization weighted reconstruction if normalization correction is
# applied.
# NOTE: Attenuation correction is always performed during reconstruction
# regardless of the choice here.
options.corrections_during_reconstruction = True
 
###########################################################################
###########################################################################
###########################################################################
###########################################################################
 
 
 
 
###########################################################################
###########################################################################
###########################################################################
####################### DYNAMIC IMAGING PROPERTIES ########################
###########################################################################
###########################################################################
###########################################################################
 
### Total time of the measurement (s)
# Use inf if you want the whole examination (static measurement only)
options.tot_time = np.inf

### Number of time points/dynamic frames (if a static measurement, use 1)
options.partitions = 1

### Start time (s) (all measurements BEFORE this will be ignored)
options.start = 0

### End time (s) (all measurements AFTER this will be ignored)
# Use inf if you want to the end of the examination (static measurement
# only)
options.end = options.tot_time
 
###########################################################################
###########################################################################
###########################################################################
###########################################################################
 
 
 
 
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
options.name = 'open_PET_data_list'

### Location of the datafile
# If no files are located in the path provided below, then the current
# folder is also checked. If no files are detected there either, an error
# is thrown.
# NOTE: for .lst or .scn files the user will be prompted for their
# locations and as such this path is ignored.
options.fpath = 'C:\\path\\to\\GATE\\output\\'

### Compute only the reconstructions
# If this file is run with this set to True, then the data load and
# sinogram formation steps are always skipped. Precomputation step is
# only performed if precompute_lor = True and precompute_all = True
# (below). Normalization coefficients are not computed even if selected.
options.only_reconstructions = False

### Show status messages
# These are e.g. time elapsed on various functions and what steps have been
# completed. It is recommended to keep this 1. This can be at most 3.
options.verbose = 1
 
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

### OpenCL/CUDA device used 
options.deviceNum = 0

### Use 64-bit integer atomic functions
# If True, then 64-bit integer atomic functions (atomic add) will be used
# if they are supported by the selected device.
# Setting this to True will make computations faster on GPUs that support
# the functions, but might make results slightly less reliable due to
# floating point rounding. Recommended for GPUs.
options.use_64bit_atomics = True

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
# 4 = Interpolation-based projector
# NOTE: You can mix and match most of the projectors. I.e. 41 will use
# interpolation-based projector for forward projection while improved
# Siddon is used for backprojection.
# See the documentation for more information:
# https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
options.projector_type = 1

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
# options.maskFP = np.ones((options.nRowsD,options.nColsD),dtype=np.uint8)
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

### Use point spread function (PSF) blurring
# Applies PSF blurring through convolution to the image space. This is the
# same as multiplying the geometric matrix with an image blurring matrix.
options.use_psf = False

# FWHM of the Gaussian used in PSF blurring in all three dimensions
# options.FWHM = [options.cr_p options.cr_p options.cr_pz]
options.FWHM = np.array([options.cr_p, options.cr_p, options.cr_pz])

# Orthogonal ray tracer (projector_type = 2) only
### The 2D (XY) width of the "strip/tube" where the orthogonal distances are
# included. If tube_width_z below is non-zero, then this value is ignored.
options.tube_width_xy = options.cr_p

# Orthogonal ray tracer (projector_type = 2) only
### The 3D (Z) width of the "tube" where the orthogonal distances are
# included. If set to 0, then the 2D orthogonal ray tracer is used. If this
# value is non-zero then the above value is IGNORED.
options.tube_width_z = options.cr_pz

# Volume ray tracer (projector_type = 3) only
### Radius of the tube-of-response (cylinder)
# The radius of the cylinder that approximates the tube-of-response.
options.tube_radius = np.sqrt(2) * (options.cr_pz / 2)

# Volume ray tracer (projector_type = 3 only)
### Relative size of the voxel (sphere)
# In volume ray tracer, the voxels are modeled as spheres. This value
# specifies the relative radius of the sphere such that with 1 the sphere
# is just large enoough to encompass an entire cubic voxel, i.e. the
# corners of the cubic voxel intersect with the sphere shell. Larger values
# create larger spheres, while smaller values create smaller spheres.
options.voxel_radius = 1

# Siddon (projector_type = 1 only)
### Number of rays
# Number of rays used per detector if projector_type = 1 (i.e. Improved
# Siddon is used) and precompute_lor = False. I.e. when using precomputed
# LOR data, only 1 rays is always used.
# Number of rays in transaxial direction
options.n_rays_transaxial = 1
# Number of rays in axial direction
options.n_rays_axial = 1
 
######################### RECONSTRUCTION SETTINGS #########################
### Number of iterations (all reconstruction methods)
options.Niter = 1

### Save specific intermediate iterations
# You can specify the intermediate iterations you wish to save here. Note
# that this uses zero-based indexing, i.e. 0 is the first iteration (not
# the initial value). By default only the last iteration is saved.
options.saveNIter = np.empty(0, dtype=np.float32)
# Alternatively you can save ALL intermediate iterations by setting the
# below to True and uncommenting it
# options.save_iter = False

### Number of subsets (all excluding MLEM and subset_type = 6)
options.subsets = 8

### Subset type (n = subsets)
# 1 = Every nth (column) measurement is taken
# 2 = Every nth (row) measurement is taken (e.g. if subsets = 3, then
# first subset has measurements 1, 4, 7, etc., second 2, 5, 8, etc.) 
# 3 = Measurements are selected randomly
# 4 = (Sinogram only) Take every nth column in the sinogram
# 5 = (Sinogram only) Take every nth row in the sinogram
# 8 = Use every nth sinogram
# 9 = Randomly select the full sinograms
# 11 = Use prime factor sampling to select the full sinograms
# Most of the time subset_type 1 is sufficient.
options.subsetType = 1

### Initial value for the reconstruction
options.x0 = np.ones((options.Nx, options.Ny, options.Nz), dtype=np.float32)
 

###########################################################################
###########################################################################
###########################################################################
######################## RECONSTRUCTION ALGORITHMS ########################
###########################################################################
###########################################################################
###########################################################################
# Reconstruction algorithms to use (choose only one algorithm and
# optionally one prior)
 
############################### ML-METHODS ################################
### Ordered Subsets Expectation Maximization (OSEM) OR Maximum-Likelihood
### Expectation Maximization (MLEM) (if subsets = 1)
# Supported by all implementations
options.OSEM = True

### Modified Row-Action Maximum Likelihood Algorithm (MRAMLA)
# Supported by implementations 1, 2, 4, and 5
options.MRAMLA = False

### Row-Action Maximum Likelihood Algorithm (RAMLA)
# Supported by implementations 1, 2, 4, and 5
options.RAMLA = False

### Relaxed Ordered Subsets Expectation Maximization (ROSEM)
# Supported by implementations 1, 2, 4, and 5
options.ROSEM = False

### Rescaled Block Iterative Expectation Maximization (RBI-EM)
# Supported by implementations 1, 2, 4, and 5
options.RBI = False

### Dynamic RAMLA (DRAMA)
# Supported by implementations 1, 2, 4, and 5
options.DRAMA = False

### Complete data OSEM (COSEM)
# Supported by implementations 1, 2, 4, and 5
options.COSEM = False

### Enhanced COSEM (ECOSEM)
# Supported by implementations 1, 2, 4, and 5
options.ECOSEM = False

### Accelerated COSEM (ACOSEM)
# Supported by implementations 1, 2, 4, and 5
options.ACOSEM = False

### FISTA
# Supported by implementations 1, 2, 4, and 5
options.FISTA = False

### FISTA with L1 regularization (FISTAL1)
# Supported by implementations 1, 2, 4, and 5
options.FISTAL1 = False

### LSQR
# Supported by implementations 1, 2, 4, and 5
options.LSQR = False

### CGLS
# Supported by implementations 1, 2, 4, and 5
options.CGLS = False
 
 
############################### MAP-METHODS ###############################
# Any algorithm selected here will utilize any of the priors selected below
# this. Note that only one algorithm and prior combination is allowed! You
# can also use most of these algorithms without priors (such as PKMA or
# PDHG).
### One-Step Late MLEM (OSL-MLEM)
# Supported by implementations 1, 2, 4, and 5
options.OSL_MLEM = False

### One-Step Late OSEM (OSL-OSEM)
# Supported by implementations 1, 2, 4, and 5
options.OSL_OSEM = False

### Modified BSREM (MBSREM)
# Supported by implementations 1, 2, 4, and 5
options.MBSREM = False

### Block Sequential Regularized Expectation Maximization (BSREM)
# Supported by implementations 1, 2, 4, and 5
options.BSREM = False

### ROSEM-MAP
# Supported by implementations 1, 2, 4, and 5
options.ROSEM_MAP = False

### RBI-OSL
# Supported by implementations 1, 2, 4, and 5
options.OSL_RBI = False

### (A)COSEM-OSL
# 0/False = No COSEM-OSL, 1/True = ACOSEM-OSL, 2 = COSEM-OSL
# Supported by implementations 1, 2, 4, and 5
options.OSL_COSEM = False

### Preconditioner Krasnoselskii-Mann algorithm (PKMA)
# Supported by implementations 1, 2, 4, and 5
options.PKMA = False

### Primal-dual hybrid gradient (PDHG)
# Supported by implementations 1, 2, 4, and 5
options.PDHG = False

### Primal-dual hybrid gradient (PDHG) with L1 minimization
# Supported by implementations 1, 2, 4, and 5
options.PDHGL1 = False

### Primal-dual hybrid gradient (PDHG) with KUllback-Leibler minimization
# Supported by implementations 1, 2, 4, and 5
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
# options.alphaPrecond = np.empty(0, dtype=np.float32)
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
options.weights = np.empty(0, dtype=np.float32)
 
 
############################## HP PROPERTIES ##############################
### Delta parameter for Huber prior
# Upper and lower bounds for the prior
options.huber_delta = 5

### Pixel weights for Huber prior
# Same rules apply as with quadratic prior weights.
# If left empty then they will be calculated by the algorithm and are based
# on the distance of the voxels from the center.
options.weights_huber = np.empty(0, dtype=np.float32)
 
 
########################### L-FILTER PROPERTIES ###########################
### Weighting factors for the L-filter pixels
# Otherwise the same as in quadratic prior, but center pixel is not Inf.
# If left empty then they will be calculated by the algorithm such that the
# weights resemble a Laplace distribution.
options.a_L = np.empty(0, dtype=np.float32)

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
options.fmh_weights = np.empty(0, dtype=np.float32)

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
options.weighted_weights = np.empty(0, dtype=np.float32)

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
########################## DEPTH OF INTERACTION ###########################
###########################################################################
# Uncomment the below value to set a depth of interaction (mm)
# NOTE: By default this is set to 0, i.e. it is assumed that all the
# interactions occur at the surface of the detector crystals. What this
# value changes is the depth of where the interactions are assumed to
# occur, i.e. it only changes the detector coordinates such that the
# transaxial coordinates are "deeper" in the crystal.
# options.DOI = 4.584
###########################################################################

# Loads the detector coordinates for each coincidence event

# You can either use the index-based reconstruction or the detector coordinates

# Use index-based reconstruction
# If true, requires options.trIndex and options.axIndex variables
# These should include the transaxial and axial, respectively, detector
# coordinate indices. Two indices are required per measurement, i.e.
# source-detector or detector-detector pairs. The indexing has to be
# zero-based! The transaxial coordinates should be stored in options.x and
# axial coordinates in options.z. The indices should correspond to the
# coordinates in each. Note that options.x should have both the x- and
# y-coordinates while options.z should have only the z-coordinates. You can
# also include randoms by inputting them as negative measurements. The
# indices are used in the same order as measurements.
options.useIndexBasedReconstruction = True

# For list-mode data, the core component you need are the detector
# coordinates for each event. options.x should be 6xNumberOfEvents, where
# there first three rows correspond to the x/y/z-coordinates of the first
# detector and the next three the x/y/z-coordinates for the second
# detector. options.x can also be a vector of size 6 * NumberOfEvents.
# The above applies when using Fourier data ordering, if you use the default
# C ordering, then the dimensions need to be NumberOfEventsx6. 
# Alternatively, you can input the coordinates separately in options.x, 
# options.y and options.z such that x-coordinates for source + detector are in
# options.x, etc. These can be either Fourier or C order and will be correctly
# ordered later. However, this consumes more intermediate memory.
# Note: List-mode reconstruction can be much more memory intensive than
# regular reconstruction
# Note: Origin is assumed to be at the center. If this is not the case, you
# can shift it with options.oOffsetX, options.oOffsetY and options.oOffsetZ
# Note: The True in the function call extracts the coordinates instead of 
# sinograms

options.randoms_correction = True
if not options.only_reconstructions:
    if options.useIndexBasedReconstruction:
        Sino, SinoD, options.trIndex, options.axIndex, DtrIndices, DaxIndices = loadInveonData(options)
        from omegatomo.projector.detcoord import detectorCoordinates
        x, y = detectorCoordinates(options)
        options.x = np.asfortranarray(np.concatenate((np.reshape(x.T, (1, -1)), np.reshape(y.T, (1,-1))),axis=0))
        z_length = options.linear_multip * options.cryst_per_block_axial * options.cr_pz
        z1 = np.linspace(-(z_length / 2 - options.cr_pz/2), z_length / 2 - options.cr_pz/2, options.rings,dtype=np.float32)
        options.z = z1
    else:
        Sino, SinoD, options.x, rand, temp1, temp2 = loadInveonData(options, True)

## Reconstructions
# List-mode reconstruction is performed for the selected algorithms. As
# long as the detector coordinates are in the options-struct, the
# list-mode/custom detector reconstruction is performed automatically. The
# detector coordinates can also correspond to the LORs in sinograms and the
# input measurement data should then be in sinogram format (matrix or
# vector).

# The measurement data
# For list-mode this is all ones as each coincidence is handled
# individually
if options.randoms_correction == False:
    if options.useIndexBasedReconstruction:
        options.SinM = np.ones((options.trIndex.shape[1]),dtype=np.float32)
    else:
        options.SinM = np.ones((options.x.shape[1]),dtype=np.float32)

# Note that if you want to add randoms, you need to add these manually:
if options.randoms_correction == True:
    if options.useIndexBasedReconstruction:
        options.SinM = np.concatenate((np.ones((options.trIndex.shape[1]),dtype=np.float32),-np.ones((DtrIndices.shape[1]),dtype=np.float32)))
        options.trIndex = np.concatenate((options.trIndex, DtrIndices),axis=1)
        options.axIndex = np.concatenate((options.axIndex, DaxIndices),axis=1)
    else:
        options.SinM = np.concatenate((np.ones((options.x.shape[1]),dtype=np.float32),-np.ones((rand.shape[1]),dtype=np.float32)))
        options.x = np.concatenate((options.x, rand),axis=1)
        del rand
# Randoms don't have to be situated last in the coordinate vector or in the
# "measurement" vector, as long as the coordinates and the "measurement"
# vector have the the same type of events in the same index.

# Since the measurement vector contains all ones, but does not contain
# every possible LOR (and also contains some duplicates) the sensitivity
# image is computed for the every possible LOR. Having this value set to
# True, the sensitivity image is computed for all applicable LORs. This
# requires correct values in the scanner properties, mainly the number of
# detectors per ring, number of rings and possible pseudo gaps. If this is
# set to False, then the sensitivity image is computed for the input
# detector coordinates only.
# NOTE: If you use your own detector coordinates, you can include the
# necessary coordinates for sensitivity image, that are missing from the
# measurements, by adding zero measurements to the corresponding
# coordinates. Coordinates that appear multiple times, however, will cause
# the sensitivity image values to be incorrect in these coordinates as the
# values are added multiple times. 
options.compute_sensitivity_image = True

# In case the data takes a lot of memory and there isn't enough memory on
# the GPU to store all the measurement data, set the below value to false
options.loadTOF = True

# In case you're running out of device (GPU) memory, you can uncomment the
# below line which only loads the measurement and coordinate data for the
# current subset
# options.loadTOF = False
# Note that the above only helps if you use subsets

pz, fp = reconstructions_main(options)

import matplotlib as plt

plt.pyplot.imshow(pz[:,:,100])