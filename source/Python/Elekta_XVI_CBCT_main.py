# -*- coding: utf-8 -*-
""" Python codes for CBCT reconstruction using Elekta XVI projection data
This example reconstructs an Elekta XVI cone-beam CT scan from the raw
projection export of the XVI software: one Heimann HIS file (.his) per
projection, 512x512 pixels, and the _Frames.xml file that contains the
gantry angle and the detector offsets (UCentre, VCentre) of each
projection. The scan is a short scan (about 200 degrees, clockwise or
counter-clockwise), so Parker weights are used. The data and the geometry
are read with loadElektaXVIData, which sets all the scanner properties
(angles, source distances, detector size, pitch and the per projection
detector offsets) and the global air value (flat). The projection images
are intensities; OMEGA computes the line integrals log(flat / intensity)
itself.
NOTE: This script uses implementation 2 (OpenCL/CUDA, GPU) by default.
"""

import numpy as np
from omegatomo.projector import proj
from omegatomo.fileio.loadElektaXVIData import loadElektaXVIData
from omegatomo.reconstruction import reconstructions_mainCT
import matplotlib.pyplot as plt

options = proj.projectorClass()

###########################################################################
###########################################################################
############################ SCANNER PROPERTIES ###########################
###########################################################################
###########################################################################

### Folder of the XVI projections
# The folder that contains the .his files and the _Frames.xml file. If left
# empty, you will be prompted to select the folder.
options.fpath = ''

### Binning
# The projections are binned (mean) from both dimensions if this is > 1.
# Binning of 2 turns the 512x512 projections into 256x256 projections
# (1.6 mm pixels).
options.binning = 1

### Geometry
# The loader uses the nominal XVI geometry: source to center of rotation
# 1000 mm, source to detector 1536 mm, raw detector pixel 0.8 mm, and the
# UCentre/VCentre values of the _Frames.xml file as the detector offsets
# (mm). To use other values, set options.sourceToCRot,
# options.sourceToDetector or options.dPitchX/options.dPitchY (raw pitch
# times binning) before loading the data. options.elektaUCentreScale and
# options.elektaVCentreScale (default 1) scale the UCentre/VCentre offsets.

### Compute only the reconstructions
# If this file is run with this set to True, then the data load will be
# skipped if the options.SinM variable exists
options.only_reconstructions = False

### Load the data
# Sets SinM (intensity), flat, angles, nRowsD, nColsD, nProjections,
# sourceToCRot, sourceToDetector, dPitchX, dPitchY, detOffsetRow and
# detOffsetCol. The global air value (flat) is the 99.9th percentile of all
# the intensities. It can be changed with options.elektaFlatPercentile, or
# by setting options.flat before the load.
loadElektaXVIData(options)

### Name of current datafile/examination
# This is used for naming purposes only
options.name = 'Elekta_XVI_CBCT'

### Show status messages
# These are e.g. time elapsed on various functions and what steps have been
# completed. It is recommended to keep this at 1 or 2. With value of 2,
# you get more detailed timing information. Maximum is 3. Minimum is 0.
options.verbose = 1

### Transaxial FOV size (mm)
# The detector covers about 266 mm at the center of rotation
# (256 pixels * 1.6 mm * 1000 / 1536) so 256 mm FOV is used.
options.FOVa_x = 256
options.FOVa_y = options.FOVa_x

### Axial FOV (mm)
# The axial coverage is similar to the transaxial one, but the cone beam
# artifacts increase towards the axial ends of the volume. A smaller FOV can
# be used if only the central part is of interest.
options.axial_fov = 256

###########################################################################
###########################################################################
############################ IMAGE PROPERTIES #############################
###########################################################################
###########################################################################

### Reconstructed image pixel count (X/row-direction)
options.Nx = 256

### Y/column-direction
options.Ny = 256

### Z-direction (number of slices) (axial)
options.Nz = 256

### Flip the image (in column direction)?
options.flip_image = False

### How much is the image rotated (radians)?
# The loader maps the XVI gantry angles to OMEGA's angles, no additional
# rotation is needed
options.offangle = 0.

###########################################################################
###########################################################################
######################## RECONSTRUCTION PROPERTIES ########################
###########################################################################
###########################################################################

### Reconstruction type
# False = FDK, True = iterative PDHG with RDP
iterative = True

### OpenCL/CUDA device used
options.deviceNum = 0

### Use CUDA
# Selecting this to True will use CUDA kernels/code instead of OpenCL. This
# only works if the CUDA code was successfully built. This is recommended
# if you have CUDA-capable device.
options.useCUDA = False

### Type of projector to use for the geometric matrix
# 1 = Improved/accelerated Siddon's algorithm
# 4 = Interpolation-based projector (ray- and voxel-based)
# 5 = Branchless distance-driven projector
options.projector_type = 4

### Interpolation length (projector type = 4 only)
# This specifies the length after which the interpolation takes place. This
# value will be multiplied by the voxel size which means that 1 means that
# the interpolation length corresponds to a single voxel (transaxial)
# length. Larger values lead to faster computation but at the cost of
# accuracy. Recommended values are between [0.5 1].
options.dL = 0.5

### Parker weights
# The scan is a short scan, so Parker weights are required with FDK. The
# weights are computed from the actual geometry (including the per
# projection detector offsets) and they include the angular scaling of FDK,
# so they must NOT be used with the iterative algorithms.
options.useParkerWeights = not iterative

### Number of iterations (iterative methods)
options.Niter = 5

### Number of subsets (iterative methods)
options.subsets = 10

### Subset type (n = subsets)
# 8 = Use every nth projection image
# 9 = Randomly select the full projection images
# 10 = Use golden angle sampling to select the subsets
# 11 = Use prime factor sampling to select the full projection images
options.subsetType = 8
if not iterative:
    options.subsets = 1

### Initial value for the reconstruction
options.x0 = np.ones((options.Nx, options.Ny, options.Nz), dtype=np.float32) * 1e-4

############################### ALGORITHMS ################################
### Feldkamp-Davis-Kress (FDK)
options.FDK = not iterative

### Primal-dual hybrid gradient (PDHG)
options.PDHG = iterative

### Relative difference prior (RDP) for PDHG
options.RDP = iterative

### The regularization parameter for the prior
options.beta = 3.

### Edge weighting factor of RDP
# Higher values sharpen the image, smaller values make it smoother
options.RDP_gamma = 10.

### Enforce positivity (PDHG)
options.enforcePositivity = True

### Measurement-based preconditioners (PDHG)
# precondTypeMeas[0] = Diagonal normalization preconditioner (1 / (A1))
# precondTypeMeas[1] = Filtering-based preconditioner
options.precondTypeMeas[1] = True

### Window of the filter (FDK and the filtering-based preconditioner)
# 'hamming' (default), 'hann', 'blackman', ... (see the docs)
options.filterWindow = 'hamming'

###########################################################################
###########################################################################
###########################################################################
###########################################################################

import time
tic = time.perf_counter()
# pz is the reconstructed attenuation image (1/mm)
pz, fp = reconstructions_mainCT(options)
toc = time.perf_counter()
print(f"{toc - tic:0.4f} seconds")

# Visualization. Soft tissue is about 0.02 1/mm at XVI energies
plt.imshow(np.squeeze(pz[:,:,options.Nz // 2]), vmin=0, vmax=0.04)

#from omegatomo.util.volume3Dviewer import volume3Dviewer
#volume3Dviewer(pz)
