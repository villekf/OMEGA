# -*- coding: utf-8 -*-
"""
## Python codes for SPECT custom algorithm reconstruction
This example contains a simplified example for custom algorithm
reconstruction using projection SPECT data. Currently the support for
some of the additional features is limited.

Note that custom algorithm refers to your own algorithms and not the built-in 
algorithms. The forward and/or backward projections of OMEGA are utilized for the computation
of these algorithms.

This example uses PyTorch and CuPy and thus requires CUDA (with CuPy and PyTorch)!
"""

import numpy as np
from omegatomo import proj
import torch

options = proj.projectorClass()

# Required for SPECT data
options.SPECT = True

# Assumes that PyTorch tensors are used as input to either forward or backward projections
options.useTorch = True

# Required for PyTorch
options.useCUDA = True

# Uses CuPy instead of PyCUDA (recommended)
options.useCuPy = True

# Use Pro.specta DICOM file?
useProSpectaData = True
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
# Positive values perform the rotation in counterclockwise direction
options.offangle = 0

###########################################################################
###########################################################################
###########################################################################
########################### PROJECTION DATA ###############################
###########################################################################
###########################################################################
###########################################################################

if useProSpectaData:
    from omegatomo.io.loadProSpectaData import loadProSpectaData
    loadProSpectaData(options)
else:
    # Gantry angles
    options.angles = np.array([0])

    # Detector swivel angles
    options.swivelAngles = options.angles+180

    # Distance between detector surface and FOV centre (origin)
    options.radiusPerProj = 48*options.crXY*np.ones_like(options.angles); 

    # Projection images for backward projection example
    options.SinM = np.zeros((128, 128, len(options.angles)), dtype=np.float32)
    options.SinM[63, 63, :] = 1

    # Number of rows in a projection image
    options.nRowsD = options.SinM.shape[0]

    # Number of columns in a projection image
    options.nColsD = options.SinM.shape[1]

    # Number of projections
    options.nProjections = options.SinM.shape[2]

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
options.verbose = 1

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
options.projector_type = 1

### Use images instead of buffers? For rotation-based projector this
# implies hardware texture interpolation, which typically has 8 bit 
# precision. With buffers, software interpolation with 32 bit floats is
# used.
options.useImages = True

###########################################################################
###########################################################################
###########################################################################
###########################################################################

options.Niter = 5 # Number of iterations
options.subsets = 8 # 1 for MLEM
# 8: every nth projection (n=options.subsets)
options.subsetType = 8 

# Initialize projector
options.addProjector()
options.initProj()

### MLEM/OSEM
m = options.SinM.ravel('F') # Measurements ()
d_m = [None] * options.subsets
d_f = torch.tensor(options.x0,device='cuda') # Transfer initial value to GPU (default = array of ones)
for k in range(options.subsets): # Split data to subsets
    d_m[k] = torch.tensor(m[options.nTotMeas[k].item() : options.nTotMeas[k + 1].item()], device='cuda')
for it in range(options.Niter):
    for k in range(options.subsets):
        # This is necessary when using subsets
        # Alternative, call options.forwardProject(d_f, k) to use forward projection
        # options.backwardProject(m, k) for backprojection
        options.subset = k
        fp = options * d_f
        Sens = options.T() * torch.ones(d_m[k].numel(), dtype=torch.float32, device='cuda')
        Sens[Sens <= 0] = options.epps
        bp = options.T() * (d_m[k] / fp)
        d_f = d_f / Sens * bp
        d_f = torch.clamp(d_f, min=1e-6)
    print(f'{"ML" if options.subsets==1 else "OS"}EM iteration {it+1}/{options.Niter} finished', end=f'{"\r" if it!=options.Niter-1 else "\n"}')

# Back to CPU
f = d_f.cpu().numpy()
f = np.reshape(f, (options.Nx[0], options.Ny[0], options.Nz[0]), order='F')

# Plot
from matplotlib import pyplot as plt
plt.imshow(pz[:,:,48], vmin=0)
plt.show()
#from omegatomo.util.volume3Dviewer import volume3Dviewer
#volume3Dviewer(f)
