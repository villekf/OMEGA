# -*- coding: utf-8 -*-
## MATLAB/Octave code for custom algorithm reconstruction of any ray-tracing compatible data
# This example showcases how to use data that is not in (standard) PET,
# SPECT or CT format. Essentially the only things that are needed are the
# FOV size, the number of voxels in each direction and the source/detector
# coordinates. The example data is a cylindrical PET data but in reality it
# could be anything. This is a simplified example and computes the
# reconstructions by using the built-in class object rather than the
# built-in algorithms. This example can be modified to compute your own
# algorithm.

# NOTE: This example has no error checking!

import numpy as np
from omegatomo.projector import proj
import matplotlib as plt

options = proj.projectorClass()


### Transaxial FOV size (mm), this is the length of the x (vertical) side
# of the FOV
options.FOVa_x = 150

### Transaxial FOV size (mm), this is the length of the y (horizontal) side
# of the FOV
options.FOVa_y = options.FOVa_x

### Axial FOV (mm)
options.axial_fov = np.floor(76.8 - 2.4/10)

# Note: Origin is assumed to be at the center. If this is not the case, you
# can shift it with options.oOffsetX, options.oOffsetY and options.oOffsetZ
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

### Use CUDA
# Selecting this to True will use CUDA kernels/code instead of OpenCL. This
# only works if the CUDA code was successfully built.
options.useCUDA = False

### Use CPU
# Selecting this to True will use CPU-based code instead of OpenCL or CUDA.
# Some features are not supported by CPU such as projector_type 4 and 5.
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
# each slice/sinogram/projection.
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

# FWHM of the Gaussian used in PSF blurring in all three dimensions (X/Y/Z)
options.FWHM = np.array([2.4, 2.4, 2.4])

# Orthogonal ray tracer (projector_type = 2 only)
### The 2D (XY) width of the "strip/tube" where the orthogonal distances are
# included. If tube_width_z below is non-zero, then this value is ignored.
options.tube_width_xy = 2.4

# Orthogonal ray tracer (projector_type = 2 only)
### The 3D (Z) width of the "tube" where the orthogonal distances are
# included. If set to 0, then the 2D orthogonal ray tracer is used. If this
# value is non-zero then the above value is IGNORED.
options.tube_width_z = 2.4

# Volume ray tracer (projector_type = 3 only)
### Radius of the tube-of-response (cylinder)
# The radius of the cylinder that approximates the tube-of-response.
# Default uses circle size that is just large enough to fit one detector
# crystal
options.tube_radius = np.sqrt(2) * (2.4 / 2)

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
# Siddon is used).
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
options.subsetType = 1

### Initial value for the reconstruction
options.x0 = np.ones((options.Nx, options.Ny, options.Nz), dtype=np.float32)

# Set the coordinates for EACH measurement
# The format for options.x variable is
# [sourceX1;sourceY1;sourceZ1;detectorX1;detectorY1;detectorZ1;sourceX2;sourceY2;sourceZ2;detectorX2...]
# where the number 1 refers to the first measurement, 2 to the second
# measurement, etc. X/Y/Z refers to the X/Y/Z (Cartesian) coordinates. This
# means that there should be SIX (6) coordinates for EACH measurement. Note
# this should be COLUMN-MAJOR, which means that column values are read
# first. As such, if you wish to use a matrix (which is fine) use
# 6xNumberOfMeasurements matrix with Fortran-ordering. Vector format is recommended though.
from pymatreader import read_mat
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

## Class example (OSEM)

# Here is an example of how to obtain the same results as above by using a
# specific MATLAB class. This is a bit more simplified from above and also
# allows more easily to use other properties files (such as
# Inveon_PET_main.m). PSF blurring will be performed automatically if it
# has been selected.

# Load data
var = read_mat('Cylindrical_PET_example_cylpet_example_new_sinograms_combined_static_200x168x703_span3.mat')
# When using custom coordinates, it is important to store the measurements
# in options.SinM BEFORE constructing the class object
options.SinM = np.float32(var['raw_SinM'])


# Construct the forward and backward projections object (you need to rerun
# this if you make any changes to the system):
options.addProjector()

# Assumes that Arrayfire arrays are used as input to either forward or backward projections
options.useAF = True

# Compute forward projection with options * f
# Compute backprojection with options.T() * y
# Load the input data into options.SinM before this step to automatically handle subsets
options.initProj()

import arrayfire as af
d_f = af.interop.np_to_af_array(options.x0)
m = options.SinM.ravel('F')


# """
# MLEM
# """
# d_m = af.interop.np_to_af_array(m)
# Sens = options.T() * af.constant(1, d_m.elements())
# for it in range(options.Niter):
#     fp = options * d_f
#     bp = options.T() * (d_m / fp)
#     d_f = d_f / Sens * bp
    
    

"""
OSEM
"""
d_m = [None] * options.subsets
for k in range(options.subsets):
    d_m[k] = af.interop.np_to_af_array(m[options.nTotMeas[k].item() : options.nTotMeas[k + 1].item()])
for it in range(options.Niter):
    for k in range(options.subsets):
        # This is necessary when using subsets
        # Alternative, call options.forwardProject(d_f, k) to use forward projection
        # options.backwardProject(m, k) for backprojection
        options.subset = k
        fp = options * d_f
        Sens = options.T() * af.constant(1, d_m[k].elements())
        Sens[Sens <= 0] = options.epps
        bp = options.T() * (d_m[k] / fp)
        d_f = d_f / Sens * bp
        af.eval(d_f)
    
f_np = d_f.to_ndarray()
f_np = np.reshape(f_np, (options.Nx[0].item(), options.Ny[0].item(), options.Nz[0].item()), order='F')
plt.pyplot.imshow(f_np[:,:,20], vmin=0)

af.sync()
af.device_gc()