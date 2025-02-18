# -*- coding: utf-8 -*-
"""
This example shows how to do custom algorithm reconstructions with CBCT data. This uses given source and center of detector panel coordinates,
as well as additional rotation by the panel. The algorithm is always primal-dual hybrid gradient, but there are several variations of it:
with or without subsets and with or without multi-resolution reconstruction. For the multi-resolution cases, be sure the turn extended FOV
on before trying them. Likewise, turn the multi-resolution off when using non-multi-resolution cases. Default is with subsets but without
multi-resolution. Filtering-based preconditioner is on.

For CUDA, there are no recommendations on which version to use. This version is pure CuPy, but PyTorch version is available too.

Example data available from: https://doi.org/10.5281/zenodo.12722386
"""
import numpy as np
from omegatomo import proj
from omegatomo.util import CTEFOVCorrection
from omegatomo.reconstruction.prepass import linearizeData
from omegatomo.util.powermethod import powerMethod
import matplotlib as plt

options = proj.projectorClass()


# Load input data
fpath = 'Planmeca_VisoG7_100kV_80mAs_500proj_kneePhantom.mat'
from pymatreader import read_mat
var = read_mat(fpath)

# Flat field corrected projections
options.SinM = var['proj']


# Number of projections
options.nProjections = var['nProjections']

# Field of view
options.FOVa_x = var['FOV'][0]
options.FOVa_y = var['FOV'][1]
options.axial_fov = var['FOV'][2]

# Number of rows and columns in a single projection image
options.nRowsD = var['nRowsD']
options.nColsD = var['nColsD']

# Object offset values from the origin, i.e. how much is the origin of the
# FOV shifted
options.oOffsetX = var['oOffset'][0]
options.oOffsetY = var['oOffset'][1]
options.oOffsetZ = var['oOffset'][2]

# Flat value
options.flat = var['flatValue']

# Distance to center of rotation and detector
options.sourceToCRot = var['sourceToCRot']
options.sourceToDetector = var['sourceToDetector']

# Detector pixel size
options.dPitchX = np.float32(var['dPitch'][0])
options.dPitchY = np.float32(var['dPitch'][1])

# Projection angles
options.angles = np.float32(var['projAngles'])

# Rotation of the detector panel
options.pitchRoll = np.float32(var['panelRot'])

# Coordinates for the source and center of the detector
options.x = var['xCoord']
options.y = var['yCoord']
options.z = var['zCoord']

del var


### Reconstructed image pixel size (X-direction)
options.Nx = 801

### Y-direction
options.Ny = 801

### Z-direction (number of slices) (axial)
options.Nz = 668


# Extrapolation and extended FOV
# You will have to make sure that extended FOV is on when using multi-resolution reconstruction below
# Furthermore, make sure multi-resolution is off otherwise (except when extended FOV is off as there is no multi-resolution reconstruction then)
options.useExtrapolation = False
options.useEFOV = False
options.transaxialEFOV = False
options.axialEFOV = True
options.transaxialExtrapolation = False
options.axialExtrapolation = True

options.useMultiResolutionVolumes = True

CTEFOVCorrection(options)

# Rotate/flip the image, i.e. the detector coordinates
options.flip_image = True
options.offangle = (3.*np.pi)/2.
# Computation device (this has no effect on CUDA)
options.deviceNum = 0
# Projector, see https://omega-doc.readthedocs.io/en/latest/selectingprojector.html for details
options.projector_type = 14
# Interpolation length, 1 means the length of one voxel
options.dL = 1

options.Niter = 1
# Number of subsets. If you use non-subset algorithm, you need to manually set this to 1
options.subsets = 20
options.subsetType = 8

options.offsetCorrection = False
# Should have no effect at the moment
options.enforcePositivity = True

# Unused at the moment
options.beta = .1
options.NLMsigma = 2.00e-3
options.RDP_gamma = 10.

options.verbose = 1

# Number of power method iterations
options.powerIterations = 10

# Needed for CT data
options.CT = True

# If True, uses CUDA
options.useCUDA = True

# If True, assumes that CuPy arrays are the input for forward and backward projections, unless useTorch = True
# Requires useCUDA = True
# Default is False
# NOTE: OMEGA is column-major, i.e. it uses Fortran ordering. It is recommended to use Fortran-ordered CuPy arrays if you need to use multi-dimensional arrays
options.useCuPy = True

options.addProjector()

linearizeData(options)


# Compute forward projection with options * f
# Compute backprojection with options.T() * y

options.initProj()
import cupy as cp
d_f = cp.asarray(options.x0)
m = options.SinM.ravel('F')
L = powerMethod(options)
sigma = 1.
theta = 1.
tau = L
# tau = 1. / 20.
# tau = 1. / 7535.
import time
tic = time.perf_counter()
    

# """
# PDHG without subsets or multi-resolution
# """
# d_m = cp.asarray(m)
# p = cp.zeros(m.size, dtype=cp.float32)
# f = cp.zeros(options.x0.size, dtype=cp.float32)
# for k in range(options.Niter):
#     apu = options * d_f
#     p = (p + sigma * (apu - d_m)) / (1 + sigma)
#     apu = options.T() * p
#     fPrev = f.copy()
#     f = f - tau * apu
#     f[f < 0] = 0
#     d_f = f + theta * (f - fPrev)
#     print('Iteration ' + str(k))


"""
PDHG with subsets but without multi-resolution
"""
d_m = [None] * options.subsets
for k in range(options.subsets):
    d_m[k] = cp.asarray(m[options.nTotMeas[k].item() : options.nTotMeas[k + 1].item()])
p = [None] * options.subsets
for k in range(options.subsets):
    p[k] = cp.zeros(d_m[k].size, dtype=cp.float32)
g1 = cp.zeros(d_f.size, dtype=cp.float32)
for k in range(options.Niter):
    for i in range(options.subsets):
        options.subset = i
        apu = options * d_f
        pl = (p[i] + sigma * (apu - d_m[i])) / (1. + sigma)
        dg = options.T() * (pl - p[i])
        p[i] = pl
        g1 = g1 + dg
        g = g1 + (theta * options.subsets) * dg
        d_f = d_f - tau * g
        d_f[d_f < 0] = 0
        print('Sub-iteration ' + str(i))
    print('Iteration ' + str(k))
    
# """
# PDHG without subsets but with multi-resolution
# """
# d_m = cp.asarray(m)
# p = cp.zeros(m.size, dtype=cp.float32)
# f = [None] * (options.nMultiVolumes + 1)
# d_f = [None] * (options.nMultiVolumes + 1)
# NN = 0
# NN2 = 0
# for i in range((options.nMultiVolumes + 1)):
#     f[i] = cp.zeros(options.N[i].item(), dtype=cp.float32)
#     NN += options.N[i].item()
#     d_f[i] = cp.asarray(options.x0[NN2 : NN])
#     NN2 += options.N[i].item()
# for k in range(options.Niter):
#     apu = options * d_f
#     p = (p + sigma * (apu - d_m)) / (1 + sigma)
#     apu = options.T() * p
#     for i in range((options.nMultiVolumes + 1)):
#         fPrev = f[i]
#         f[i] = f[i] - tau[i] * apu[i]
#         f[i][f[i] < 0] = 0
#         d_f[i] = f[i] + theta * (f[i] - fPrev)
#     print('Iteration ' + str(k))


# """
# PDHG with subsets and multi-resolution
# """
# d_m = [None] * options.subsets
# for k in range(options.subsets):
#     d_m[k] = cp.asarray(m[options.nTotMeas[k].item() : options.nTotMeas[k + 1].item()])
# p = [None] * options.subsets
# for k in range(options.subsets):
#     p[k] = cp.zeros(d_m[k].size, dtype=cp.float32)
# g1 = [None] * (options.nMultiVolumes + 1)
# d_f = [None] * (options.nMultiVolumes + 1)
# NN = 0
# NN2 = 0
# for i in range((options.nMultiVolumes + 1)):
#     g1[i] = cp.zeros(options.N[i].item(), dtype=cp.float32)
#     NN += options.N[i].item()
#     d_f[i] = cp.asarray(options.x0[NN2 : NN])
#     NN2 += options.N[i].item()
# for k in range(options.Niter):
#     for i in range(options.subsets):
#         options.subset = i
#         apu = options * d_f
#         pl = (p[i] + sigma * (apu - d_m[i])) / (1 + sigma)
#         dg = options.T() * (pl - p[i])
#         p[i] = pl.copy()
#         for j in range((options.nMultiVolumes + 1)):
#             g1[j] = g1[j] + dg[j]
#             g = g1[j] + (theta * options.subsets) * dg[j]
#             d_f[j] = d_f[j] - tau[j] * g
#             d_f[j][d_f[j] < 0] = 0
#     print('Iteration ' + str(k))

toc = time.perf_counter()
print(f"{toc - tic:0.4f} seconds")
if isinstance(d_f, list):
    f_np = d_f[0].get()
else:
    f_np = d_f.get()
f_np = cp.reshape(f_np, (options.Nx[0].item(), options.Ny[0].item(), options.Nz[0].item()), order='F')
plt.pyplot.imshow(f_np[:,:,200], vmin=0)