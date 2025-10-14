# -*- coding: utf-8 -*-
"""
## Python codes for SPECT custom algorithm reconstruction
This example contains a simplified example for custom algorithm
reconstruction using projection SPECT data. In this case the
data is Siemens Pro.specta projection data available at DOI
10.5281/zenodo.17315440. Currently the support for
some of the additional features is limited.

Note that custom algorithm refers to your own algorithms and not the built-in 
algorithms. The forward and/or backward projections of OMEGA are utilized for the computation
of these algorithms.

This example uses PyTorch and CuPy and thus requires CUDA (with CuPy and PyTorch)!
"""

import numpy as np
from omegatomo.projector import proj
import torch
import h5py

options = proj.projectorClass()

# Required for SPECT data
options.SPECT = True

# Assumes that PyTorch tensors are used as input to either forward or backward projections
options.useTorch = True

# Required for PyTorch
options.useCUDA = True

# Uses CuPy instead of PyCUDA (recommended)
options.useCuPy = True

options.fpath = '' # Path to .mat file

###########################################################################
###########################################################################
###########################################################################
############################### LOAD DATA #################################
###########################################################################
###########################################################################
###########################################################################

# ---------------------------
# Helpers for MATLAB v7.3 HDF5
# ---------------------------
def _read(obj):
    if isinstance(obj, h5py.Dataset):
        data = obj[()]
        try:
            return np.array(data)
        except Exception:
            return data
    elif isinstance(obj, h5py.Group):
        out = {}
        for k, v in obj.items():
            out[k] = _read(v)
        return out
    else:
        return obj

def load_mat73(path):
    out = {}
    with h5py.File(path, 'r') as f:
        for k in f.keys():
            try:
                out[k] = _read(f[k])
            except Exception:
                pass
    return out

mat = load_mat73(options.fpath)
options.SinM = np.array(mat['projection_data'])
options.SinM = np.transpose(options.SinM, (2, 1, 0)) # Fortran vs C order
options.angles = np.array(mat['angular_position']).squeeze()
options.radiusPerProj = np.array(mat['radial_position']).squeeze()
options.nRowsD = options.SinM.shape[0]
options.nColsD = options.SinM.shape[1]
options.nProjections = options.SinM.shape[2]
energy_window = mat.get("energy_window", None)
pixel_spacing = np.array(mat["pixel_spacing"]).squeeze()
detector_thickness = float(np.squeeze(mat["detector_thickness"]))

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
options.colL = float(np.squeeze(mat["collimator_thickness"]))
# Collimator hole radius (mm)
options.colR = float(np.squeeze(mat["collimator_hole_radius"]))
# Distance from collimator to the detector (mm)
options.colD = 0.0
# Intrinsic resolution (mm)
options.iR = float(np.squeeze(mat["detector_intrinsic_resolution"]))
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
# options.rayShiftsDetector = [];
# options.rayShiftsSource = [];
 
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

# This has to be True if you want to use the filtering-based preconditioner
options.PDHG = False

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
plt.imshow(f[:,:,39], vmin=0)
plt.show()
#from omegatomo.util.volume3Dviewer import volume3Dviewer
#volume3Dviewer(f)
