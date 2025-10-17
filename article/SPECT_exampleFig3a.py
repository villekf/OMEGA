# -*- coding: utf-8 -*-
"""
Python codes for SPECT reconstruction from projection images

This example outlines the reconstruction of SPECT data. In this case the
data is Siemens Pro.specta projection data available at DOI
10.5281/zenodo.17315440
"""
import numpy as np
from omegatomo.projector import proj
from omegatomo.reconstruction import reconstructions_mainSPECT
from pymatreader import read_mat

# Initialize projector class for reconstruction
options = proj.projectorClass()

options.fpath = '' # Path to .mat file

###########################################################################
###########################################################################
###########################################################################
############################### LOAD DATA #################################
###########################################################################
###########################################################################
###########################################################################

data = read_mat(options.fpath)
options.SinM = np.array(data['projection_data'])
options.angles = np.array(data['angular_position']).squeeze()
options.radiusPerProj = np.array(data['radial_position']).squeeze()
options.nRowsD = options.SinM.shape[0]
options.nColsD = options.SinM.shape[1]
options.nProjections = options.SinM.shape[2]
energy_window = data.get("energy_window", None)
pixel_spacing = np.array(data["pixel_spacing"]).squeeze()
detector_thickness = float(np.squeeze(data["detector_thickness"]))

###########################################################################
###########################################################################
###########################################################################
########################### SCANNER PROPERTIES ############################
###########################################################################
###########################################################################
###########################################################################

### Crystal thickness (mm)
options.cr_p = detector_thickness

### Crystal width (mm)
options.dPitchX = float(pixel_spacing[0])
options.dPitchY = float(pixel_spacing[1])

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
options.Nx = 64; # X-direction
options.Ny = 64; # Y-direction
options.Nz = 128; # Z-direction (number of axial slices)

### FOV size [mm]
# NOTE: Non-cubical voxels may not work
options.FOVa_x = options.dPitchX*64; # [mm], x-axis of FOV (transaxial)
options.FOVa_y = options.dPitchX*64; # [mm], y-axis of FOV (transaxial)
options.axial_fov = options.dPitchY*128; # [mm], z-axis of FOV (axial)

### Flip the image?
options.flipImageX = False
options.flipImageY = False
options.flipImageZ = False

### Use back projection mask?
options.useMaskBP = False
options.maskBP = np.ones((options.Nx, options.Ny, options.Nz))

### How much is the image rotated in degrees?
# NOTE: The rotation is done in the detector space (before reconstruction).
# Positive values perform the rotation in counterclockwise direction
options.offangle = 0

###########################################################################
###########################################################################
###########################################################################
############################## CORRECTIONS ################################
###########################################################################
###########################################################################
###########################################################################

######################### Attenuation correction ##########################
# Currently scaling and resampling is not supported for the attenuation map.
options.attenuation_correction = False

######################### Normalization correction ########################
# If set to true, normalization correction is applied to either the
# projection data or in the image reconstruction by using predefined
# normalization coefficients.
options.normalization_correction = False
options.normalization = np.ndarray([])

############################ Scatter correction ###########################
# Uses linear interpolation between scatter windows. options.ScatterC{1} 
# contains the lower scatter window and options.ScatterC{2} contains the 
# upper scatter window (sizes equal options.SinM).
# See for example: 10.1371/journal.pone.0269542
options.scatter_correction = False
options.ScatterC = np.ndarray([])
options.eWin = np.array(energy_window).squeeze() # Main energy window: [lowerLimit upperLimit]
options.eWinL = None  # Lower energy window: [lowerLimit upperLimit]
options.eWinU = None  # Upper energy window: [lowerLimit upperLimit]

########################### Resolution recovery ##########################
### Collimator-detector response function (CDRF)
# For projector types 2 and 6 you can either input either:
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
#
# 1. The collimator parameters (projector types 1, 2 and 6)
# Collimator hole length (mm)
options.colL = float(np.squeeze(data["collimator_thickness"]))
# Collimator hole radius (mm)
options.colR = float(np.squeeze(data["collimator_hole_radius"]))
# Distance from collimator to the detector (mm)
options.colD = 0.0
# Intrinsic resolution (mm)
options.iR = float(np.squeeze(data["detector_intrinsic_resolution"]))
# Focal distance (XY)
options.colFxy = np.inf
# Focal distance (Z)
options.colFz = np.inf

# 2. If you have the standard deviations for transaxial (XY) and axial (Z)
# directions, you can input them here instead of the above values The
# dimensions need to be options.nProjections x options.Nx. Only for
# projector type 6.
# options.sigmaZ = np.ones((options.nProjections, options.Nx), dtype=np.float32)
# options.sigmaXY = np.ones((options.nProjections, options.Nx), dtype=np.float32)

# 3. You can input the filter for the CDRF directly. This should be of the
# size filterSizeXY x filterSizeZ. Only for
# projector type 6.
# options.gFilter = np.ones((1, 1, options.Nx), dtype=np.float32)

# 4. For the Siddon ray tracer, the CDRF is defined by shifting the rays to
# the shape of the collimator hole. The values of rayShiftsDetector and
# rayShiftsSource represent [shift1XY, shift1Z, shift2XY, ...] in mm. Size
# should be 2*nRays x nColsD x nRowsD x nProjections. If not input, values
# are calculated automatically.
options.nRays = 1  # Number of rays traced per detector element
# options.rayShiftsDetector = np.zeros((2*options.nRays, options.nColsD, options.nRowsD, options.nProjections));
# options.rayShiftsSource = np.zeros((2*options.nRays, options.nColsD, options.nRowsD, options.nProjections));

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
options.useImages = False

######################### RECONSTRUCTION SETTINGS #########################
### Number of iterations (all reconstruction methods)
options.Niter = 5
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
######################## RECONSTRUCTION ALGORITHMS ########################
###########################################################################
###########################################################################
###########################################################################

# See examples in main-files folder for more algorithms.
options.OSEM = True


## Reconstructions
import time
tic = time.perf_counter()
pz, fp = reconstructions_mainSPECT(options)
toc = time.perf_counter()
print(f"{toc - tic:0.4f} seconds")

# Plot
from matplotlib import pyplot as plt
plt.imshow(pz[:,:,39], vmin=0)
plt.show()
#from omegatomo.util.volume3Dviewer import volume3Dviewer
#volume3Dviewer(pz)