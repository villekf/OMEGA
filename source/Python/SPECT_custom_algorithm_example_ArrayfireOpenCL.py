# -*- coding: utf-8 -*-
"""
## Python codes for SPECT custom algorithm reconstruction
This example contains a simplified example for custom algorithm
reconstruction using projection SPECT data. Currently the support for
some of the additional features is limited.

Note that custom algorithm refers to your own algorithms and not the built-in 
algorithms. The forward and/or backward projections of OMEGA are utilized for the computation
of these algorithms.

This example uses Arrayfire with PyOpenCL and thus requires OpenCL (with PyOpenCL and Arrayfire)!
"""

import numpy as np
from omegatomo import proj
import matplotlib.pyplot as plt
import arrayfire as af

options = proj.projectorClass()

# Required for SPECT data
options.SPECT = True

# Assumes that Arrayfire arrays are used as input to either forward or backward projections
options.useAF = True

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
options.crXY = 4.7952

### Scanner name
# Used for naming purposes (measurement data)
options.machine_name = 'Example'
 
###########################################################################
###########################################################################
###########################################################################
########################### IMAGE PROPERTIES ##############################
###########################################################################
###########################################################################
###########################################################################
 
### Reconstructed image pixel count
# NOTE: Non-square image sizes (X- and Y-direction) may not work
options.Nx = 128 # X-direction
options.Ny = 128 # Y-direction
options.Nz = 128 # Z-direction (number of axial slices)

### FOV size [mm]
# NOTE: Non-cubical voxels may not work
options.FOVa_x = options.crXY*128 # [mm], x-axis of FOV (transaxial)
options.FOVa_y = options.crXY*128 # [mm], y-axis of FOV (transaxial)
options.axial_fov = options.crXY*128 # [mm], z-axis of FOV (axial)

### Flip the image?
options.flipImageX = False
options.flipImageY = False
options.flipImageZ = False

### Use back projection mask?
options.useMaskBP = False

### Attenuation correction
options.attenuation_correction = False

# Linear attenuation coefficients, size and dimensions should match FOV
options.vaimennus = np.zeros((options.Nx, options.Ny, options.Nz), dtype=np.float32)

### How much is the image rotated in degrees?
# NOTE: The rotation is done in the detector space (before reconstruction).
# Positive values perform the rotation in clockwise direction
options.offangle = 0

###########################################################################
###########################################################################
###########################################################################
########################### PROJECTION DATA ###############################
###########################################################################
###########################################################################
###########################################################################

# Gantry angles
options.angles = np.array([0])

# Detector swivel angles
options.swivelAngles = options.angles+180

# Distance between detector surface and FOV centre (origin)
options.radiusPerProj = 48*options.crXY*np.ones_like(options.angles); 

# Initial value for the forward projection example
x0 = np.zeros((options.Nx, options.Ny, options.Nz), dtype=np.float32)
x0[63, 63, 63] = 1


# Projection images for backward projection example
y0 = np.zeros((128, 128, len(options.angles)), dtype=np.float32)
y0[63, 63, :] = 1

# Number of rows in a projection image
options.nRowsD = y0.shape[0]

# Number of columns in a projection image
options.nColsD = y0.shape[1]

# Number of projections
options.nProjections = y0.shape[2]

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
# options.sigmaZ = np.ones((options.nProjections, options.Nx), dtype=np.float32)
# options.sigmaXY = np.ones((options.nProjections, options.Nx), dtype=np.float32) # Transaxial standard deviation, (projector type 6)

# 3. You can input the filter for the CDRF directly. This should be of the
# size filterSizeXY x filterSizeZ x options.nProjections. Only for
# projector type 6.
# options.gFilter = np.ones((10, 10, options.nProjections), dtype=np.float32)

# 4. For the Siddon ray tracer, the CDRF is defined by shifting the rays to
# the shape of the collimator hole. The below example is for random
# (uniform distribution) rays with one square collimator hole at the centre
# of each detector element. For 1 ray, the ray perpendicular to the
# detector element.
options.nRays = 1 # Number of rays traced per detector element
# options.rayShiftsDetector = options.colR*(2*np.random.rand(2 * options.nRays, 1).astype(np.float32)-1)/options.crXY # The relative shifts (dx1, dy1, dx2, dy2, ...) at the collimator-detector interface
# options.rayShiftsSource = options.colR*(2*np.random.rand(2 * options.nRays, 1).astype(np.float32)-1)/options.crXY # The relative shifts (dx1, dy1, dx2, dy2, ...) at the other end of the collimator
 
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
options.verbose = 3

###########################################################################
###########################################################################
###########################################################################
######################## RECONSTRUCTION PROPERTIES ########################
###########################################################################
###########################################################################
###########################################################################
 
############################### PROJECTOR #################################
### Type of projector to use for the geometric matrix
# 1 = (Improved) Siddon ray-based projector
# 2 = Orthogonal distance ray tracing
# 6 = Rotation-based projector
# See the documentation on some details on the projectors:
# https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
# NOTE: with rotation-based projector, the sinogram must be resized and
# resampled to match FOV XZ-plane size and resolution.
options.projector_type = 6

### Use images instead of buffers? For rotation-based projector this
# implies hardware texture interpolation, which typically has 8 bit 
# precision. With buffers, software interpolation with 32 bit floats is
# used.
options.useImages = True

###########################################################################
###########################################################################
###########################################################################
###########################################################################

# Initialize projector
options.addProjector()
options.initProj()
A = options

d_x0 = af.interop.np_to_af_array(x0.ravel('F')) # Input as Arrayfire array
d_y0 = af.interop.np_to_af_array(y0.ravel('F')) # Input as Arrayfire array

# Compute forward projection with A * x
d_y = A * d_x0

y = d_y.to_ndarray()
y = np.reshape(y, (options.nRowsD, options.nColsD, options.nProjections), order='F')
plt.subplot(1, 2, 1)
plt.imshow(y[:,:,0], vmin=0) # Plot first sinogram
plt.title('First projection image', fontsize=6)

# Compute backprojection with A.T() * y
d_x = A.T() * d_y0

x = d_x.to_ndarray()
x = np.reshape(x, (options.Nx[0], options.Ny[0], options.Nz[0]), order='F')
plt.subplot(1, 2, 2)
plt.imshow(x[:, :, 63], vmin=0) # Plot cross-section of backprojection
plt.title('Cross-section of back projection', fontsize=6)

plt.show()

af.sync()
af.device_gc()