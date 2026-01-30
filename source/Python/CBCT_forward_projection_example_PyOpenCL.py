# -*- coding: utf-8 -*-
"""
This example shows how to do forward projections in OMEGA using (CB)CT data
"""
import numpy as np
from omegatomo import proj
import matplotlib as plt

A = proj.projectorClass()


# Load input data
fpath = 'Planmeca_VisoG7_100kV_80mAs_500proj_kneePhantom.mat'
from pymatreader import read_mat
var = read_mat(fpath)

# Flat field corrected projections
A.SinM = var['proj']


# Number of projections
A.nProjections = var['nProjections']

# Field of view
A.FOVa_x = var['FOV'][0]
A.FOVa_y = var['FOV'][1]
A.axial_fov = var['FOV'][2]

# Number of rows and columns in a single projection image
A.nRowsD = var['nRowsD']
A.nColsD = var['nColsD']

# Object offset values from the origin, i.e. how much is the origin of the
# FOV shifted
A.oOffsetX = var['oOffset'][0]
A.oOffsetY = var['oOffset'][1]
A.oOffsetZ = var['oOffset'][2]

# Flat value
A.flat = var['flatValue']

# Distance to center of rotation and detector
A.sourceToCRot = var['sourceToCRot']
A.sourceToDetector = var['sourceToDetector']

# Detector pixel size
A.dPitchX = var['dPitch'][0]
A.dPitchY = var['dPitch'][1]

# Projection angles
A.angles = np.float32(var['projAngles'])

# Rotation of the detector panel
A.pitchRoll = np.float32(var['panelRot'])

# Coordinates for the source and center of the detector
A.x = var['xCoord']
A.y = var['yCoord']
A.z = var['zCoord']

del var

# This example reduces the resolution of the detector panel by two
# This is done by first reducing the number of pixels by half and then increasing
# the size of the pixels by two
# nRowsD and nColsD affects the number of pixels, while dPitch values affect the
# size (in mm) of these pixels
A.nRowsD //= 2
A.nColsD //= 2
A.dPitchX *= 2.
A.dPitchY *= 2.

# Increasing the number of projections requires more coordinates to be input, 
# but reducing the number of projections simply means removing elements from
# A.x, A.y, A.z, A.angles and A.pitchRoll
# The below example halves the number of projections
A.nProjections //= 2
A.x = A.x[::2,:]
A.y = A.y[::2,:]
A.z = A.z[::2,:]
A.angles = A.angles[::2]
A.pitchRoll = A.pitchRoll[::2,:]

### Input image pixel size (X-direction)
A.Nx = 801

### Y-direction
A.Ny = 801

### Z-direction (number of slices) (axial)
A.Nz = 668

# Rotate/flip the image, i.e. the detector coordinates
A.flip_image = True
A.offangle = (3.*np.pi)/2.
# Computation device
A.deviceNum = 0
# Projector, see https://omega-doc.readthedocs.io/en/latest/selectingprojector.html for details
A.projector_type = 4
# Interpolation length, 1 means the length of one voxel
A.dL = 1

A.Niter = 1
# Number of subsets. If you use non-subset algorithm, you need to manually set this to 1
A.subsets = 1
A.subsetType = 8

# Apply offset correction
A.offsetCorrection = False

A.verbose = 1

# Needed for CT data
A.CT = True

# The above should work in all cases, but the below part should be adjusted if
# you use CuPy, Arrayfire or PyTorch instead!
# If True, uses CUDA
# options.useCUDA = True

# If True, assumes that CuPy arrays are the input for forward and backward projections, unless useTorch = True
# Requires useCUDA = True
# Default is False
# NOTE: OMEGA is column-major, i.e. it uses Fortran ordering. It is recommended to use Fortran-ordered CuPy arrays if you need to use multi-dimensional arrays
# options.useCuPy = True

# You'll need to adjust all the variables before this point!
# Modifying the geometries after this point most likely won't work
A.addProjector()


A.initProj()
# Modify these to use CuPy if you work with that instead!
import pyopencl as cl
d_f = cl.array.to_device(A.queue, A.x0)

# The forward projection
fp = A * d_f

# Convert into NumPy array
f_np = fp.get(A.queue)
f_np = np.reshape(f_np, (A.nRowsD, A.nColsD, -1), order='F')
plt.pyplot.imshow(f_np[:,:,200], vmin=0)
