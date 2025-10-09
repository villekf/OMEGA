# -*- coding: utf-8 -*-
## MATLAB/Octave code for reconstruction of any ray-tracing compatible data
# This example showcases how to use data that is not in (standard) PET,
# SPECT or CT format. Essentially the only things that are needed are the
# FOV size, the number of voxels in each direction and the source/detector
# coordinates. The example data is a cylindrical PET data but in reality it
# could be anything. This is a simplified example.

# NOTE: This example has no error checking!

import numpy as np
from omegatomo.projector import proj
import matplotlib as plt
from omegatomo.reconstruction import reconstructions_main

options = proj.projectorClass()


### Transaxial FOV size (mm), this is the length of the x (vertical/row) side
# of the FOV
options.FOVa_x = 300

### Transaxial FOV size (mm), this is the length of the y (horizontal/column) side
# of the FOV
options.FOVa_y = options.FOVa_x

### Axial FOV (mm)
options.axial_fov = np.floor(76.8 - 2.4/10)

# Note: Origin is assumed to be at the center. If this is not the case, you
# can shift it with options.oOffsetX, options.oOffsetY and options.oOffsetZ
# That is row, column and slice directions
# options.oOffsetX = 0
# options.oOffsetY = 0
# options.oOffsetZ = 0


###########################################################################
###########################################################################
###########################################################################
########################### IMAGE PROPERTIES ##############################
###########################################################################
###########################################################################
###########################################################################

### Reconstructed image pixel count (X/row-direction)
# NOTE: Non-square image sizes (X- and Y-direction) may not work
options.Nx = 256

### Y/column-direction
options.Ny = 256

### Z-direction (number of slices) (axial)
options.Nz = 63

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

### Show status messages
# These are e.g. time elapsed on various functions and what steps have been
# completed. It is recommended to keep this at 1 or 2. With value of 2, 
# you get more detailed timing information. Maximum is 3, minimum 0.
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

############################# IMPLEMENTATIONS #############################
### OpenCL/CUDA device used 
# NOTE: Use 
# from omegatomo.util.devinfo import deviceInfo
# deviceInfo(True)
# to determine the device nubmers
options.deviceNum = 0

### Use 64-bit integer atomic functions
# If True, then 64-bit integer atomic functions (atomic add) will be used
# if they are supported by the selected device.
# Setting this to True will make computations faster on GPUs that support
# the functions, but might make results slightly less reliable due to
# floating point rounding. Recommended for OpenCL GPUs. Not recommended for
# CUDA. Doesn't apply for CPU.
options.use_64bit_atomics = True

### Use 32-bit integer atomic functions
# If True, then 32-bit integer atomic functions (atomic add) will be used.
# This is even faster than the above 64-bit atomics version, but will also
# have significantly higher reduction in numerical/floating point accuracy.
# This should be about 20-30# faster than the above 64-bit version, but
# might lead to integer overflow if you have a high count measurement
# (thousands of coincidences per sinogram bin). Use this only if speed is
# of utmost importance. 32-bit atomics take precedence over 64-bit ones,
# i.e. if options.use_32bit_atomics = true then the 64-bit version will be 
# always set as false.
options.use_32bit_atomics = False

### Use CUDA
# Selecting this to True will use CUDA kernels/code instead of OpenCL. This
# only works if the CUDA code was successfully built. This is recommended
# if you have CUDA-capable device.
options.useCUDA = False

### Use CPU
# Selecting this to True will use CPU-based code instead of OpenCL or CUDA.
# Not recommended, even OpenCL with CPU should be used before this.
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
# NOTE 2: The below additional options apply also in hybrid cases as long
# as the other projector is the corresponding projector.
# See the documentation for more information:
# https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
options.projector_type = 1

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
# value will be multiplied by the voxel size which means that the
# interpolation length of 1 corresponds to a single voxel (transaxial) 
# length. Larger values lead to faster computation but at the cost of
# accuracy. Recommended values are between [0.5 1], though values up to 2
# should be fine.
options.dL = 0.5

### Use point spread function (PSF) blurring
# Applies PSF blurring through convolution to the image space. This is the
# same as multiplying the geometric matrix with an image blurring matrix.
options.use_psf = False

# FWHM (mm) of the Gaussian used in PSF blurring in all three dimensions (X/Y/Z)
options.FWHM = np.array([2.4, 2.4, 2.4])

# Orthogonal ray tracer (projector_type = 2 only)
### The 2D (XY) width of the "strip/tube" where the orthogonal distances are
# included. If tube_width_z below is non-zero, then this value is ignored.
options.tube_width_xy = 2.4

# Orthogonal ray tracer (projector_type = 2 only)
### The 3D (Z) width (mm) of the "tube" where the orthogonal distances are
# included. If set to 0, then the 2D orthogonal ray tracer is used. If this
# value is non-zero then the above value is IGNORED.
# If you want the projector to be a tube, use this, if you want it to be 
# strip, use the above
# This slows down the reconstruction, but makes it more accurate
options.tube_width_z = 2.4

# Volume ray tracer (projector_type = 3) only
### Radius (mm) of the tube-of-response (cylinder)
# The radius of the cylinder that approximates the tube-of-response.
options.tube_radius = np.sqrt(2) * (2.4 / 2)

# Volume ray tracer (projector_type = 3 only)
### Relative size of the voxel (sphere)
# In volume ray tracer, the voxels are modeled as spheres. This value
# specifies the relative radius of the sphere such that with 1 the sphere
# is just large enough to encompass an entire cubic voxel, i.e. the
# corners of the cubic voxel intersect with the sphere shell. Larger values
# create larger spheres, while smaller values create smaller spheres.
options.voxel_radius = 1

# projector_type = 1 and 4 only
### Number of rays
# Number of rays used per detector if projector_type = 1 (i.e. Improved
# Siddon is used) or projector_type = 4 (interpolation).
# The total number of rays per detector is the multiplication of the two
# below values!
# Number of rays in transaxial (row) direction
options.n_rays_transaxial = 1;
# Number of rays in axial (column) direction
options.n_rays_axial = 1;

######################### RECONSTRUCTION SETTINGS #########################
### Number of iterations (all reconstruction methods)
options.Niter = 5

### Save specific intermediate iterations
# You can specify the intermediate iterations you wish to save here. Note
# that this uses zero-based indexing, i.e. 0 is the first iteration (not
# the initial value). By default only the last iteration is saved.
# Note: Subiterations cannot be saved!
options.saveNIter = np.empty(0)
# Alternatively you can save ALL intermediate iterations by setting the
# below to True and uncommenting it
# Note: Only affects full iterations (epochs)
# options.save_iter = False

### Number of subsets (excluding subset_type = 6 and algorithms that do not
### support subsets)
options.subsets = 8

### Subset type (n = subsets)
# 0 = Measurements are divided into n segments
# 1 = Every nth measurement is taken
# 3 = Measurements are selected randomly
options.subsetType = 3

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
# Note: OSEM, MRAMLA, RAMLA, ROSEM, RBI, COSEM, ACOSEM, ECOSEM are not
# recommended for non-Poisson data
options.OSEM = False

### Modified Row-Action Maximum Likelihood Algorithm (MRAMLA)
options.MRAMLA = False

### Row-Action Maximum Likelihood Algorithm (RAMLA)
options.RAMLA = False

### Relaxed Ordered Subsets Expectation Maximization (ROSEM)
options.ROSEM = False

### Rescaled Block Iterative Expectation Maximization (RBI-EM)
options.RBI = False

### Dynamic RAMLA (DRAMA)
options.DRAMA = False

### Complete data OSEM (COSEM)
options.COSEM = False

### Enhanced COSEM (ECOSEM)
options.ECOSEM = False

### Accelerated COSEM (ACOSEM)
options.ACOSEM = False

### FISTA
options.FISTA = False

### FISTA with L1 regularization (FISTAL1)
options.FISTAL1 = False

### LSQR
options.LSQR = False

### CGLS
options.CGLS = False


############################### MAP-METHODS ###############################
# Any algorithm selected here will utilize any of the priors selected below
# this. Note that only one algorithm and prior combination is allowed! You
# can also use most of these algorithms without priors (such as PKMA or
# PDHG).
### Modified BSREM (MBSREM)
options.MBSREM = False

### Block Sequential Regularized Expectation Maximization (BSREM)
options.BSREM = False

### Preconditioned Krasnoselskii-Mann algorithm (PKMA)
options.PKMA = False

### Primal-dual hybrid gradient (PDHG)
options.PDHG = True

### Primal-dual hybrid gradient (PDHG) with L1 minimization
options.PDHGL1 = False

### Primal-dual hybrid gradient (PDHG) with Kullback-Leibler minimization
options.PDHGKL = False

### Primal-dual Davis-Yin (PDDY)
options.PDDY = False

# You can input other reconstruction parameters and priors as in the other
# examples

# Input the measurement data
from pymatreader import read_mat
var = read_mat('Cylindrical_PET_example_cylpet_example_new_sinograms_combined_static_200x168x703_span3.mat')
options.SinM = np.float32(var['raw_SinM'])

# Set the coordinates for EACH measurement
# The format for options.x variable is
# [sourceX1;sourceY1;sourceZ1;detectorX1;detectorY1;detectorZ1;sourceX2;sourceY2;sourceZ2;detectorX2...]
# where the number 1 refers to the first measurement, 2 to the second
# measurement, etc. X/Y/Z refers to the X/Y/Z (Cartesian) coordinates. This
# means that there should be SIX (6) coordinates for EACH measurement. Note
# this should be COLUMN-MAJOR, which means that column values are read
# first. As such, if you wish to use a matrix (which is fine) use
# 6xNumberOfMeasurements matrix with Fortran-ordering. Vector format is recommended though.
var = read_mat('cylpet_example_det_coord.mat')

x = var['x']
z_det = var['z_det']

sizeX = x.shape[1]
x = np.tile(x, (1, z_det.shape[1]))

z = np.repeat(z_det, sizeX, axis=1)

options.x = np.asfortranarray(np.vstack([x[0, :], x[1, :], z[0, :], x[2, :], x[3, :], z[1, :]]))

# Under normal situations, you would want to enable the "CT" mode which
# uses the length of intersection rather than the probability
# options.CT = True

## Reconstructions


import time
tic = time.perf_counter()
pz, fp = reconstructions_main(options)
toc = time.perf_counter()
print(f"{toc - tic:0.4f} seconds")
plt.pyplot.imshow(pz[:,:,30], vmin=0)