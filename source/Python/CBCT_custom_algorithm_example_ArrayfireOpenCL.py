# -*- coding: utf-8 -*-
"""
This example shows how to do custom algorithm reconstructions with CBCT data. This uses given source and center of detector panel coordinates,
as well as additional rotation by the panel. The algorithm is always primal-dual hybrid gradient, but there are several variations of it:
with or without subsets and with or without multi-resolution reconstruction. For the multi-resolution cases, be sure the turn extended FOV
on before trying them. Likewise, turn the multi-resolution off when using non-multi-resolution cases. Default is with subsets but without
multi-resolution. Filtering-based preconditioner is on.

This example is the recommended way for OpenCL reconstructions, although a PyOpenCL version exists too.

Example data available from: https://doi.org/10.5281/zenodo.12722386
"""
import numpy as np
from omegatomo.projector import proj
from omegatomo.util import CTEFOVCorrection
from omegatomo.reconstruction.prepass import linearizeData
from omegatomo.util.powermethod import powerMethod
from omegatomo.util.measprecond import applyMeasPreconditioning, circulantInverse
from omegatomo.util.priors import RDP
from omegatomo.util.priors import NLReg
from omegatomo.util.priors import TV
import matplotlib as plt

# Note that the name can be anything
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


### Reconstructed image pixel count (X/row-direction)
options.Nx = 801

### Y/column-direction
options.Ny = 801

### Z-direction (number of slices) (axial)
options.Nz = 668

# Extrapolation and extended FOV
# You will have to make sure that extended FOV is on when using multi-resolution reconstruction below
# Furthermore, make sure multi-resolution is off otherwise (except when extended FOV is off as there is no multi-resolution reconstruction then)

### Use projection extrapolation
# If True, extrapolates the projection data. You can select below whether
# this extrapolation is done only in the axial or transaxial directions, or
# both. Default extrapolation length is 20% of the original length, for
# both sides. For example if axial extrapolation is enabled, then the left
# and right regions of the projection get 20% increase in size. This value
# can be adjusted in CTEFOVCorrection. The values are scaled to air with
# the use of logarithmic scaling.
options.useExtrapolation = False

### Use extended FOV
# Similar to above, but expands the FOV. The benefit of expanding the FOV
# this way is to enable to the use of multi-resolution reconstruction or
# computation of the priors/regularization only in the original FOV. The
# default extension is 40% per side (see below)
options.useEFOV = False

# Use transaxial extended FOV
options.transaxialEFOV = True

# Use axial extended FOV (this is on by default. If both this and
# transaxialEFOV are False but useEFOV is True, the axial EFOV will be
# turned on)
options.axialEFOV = True

options.transaxialExtrapolation = False
options.axialExtrapolation = True

# Setting this to True uses multi-resolution reconstruction when using
# extended FOV. Only applies to extended FOV!
options.useMultiResolutionVolumes = False

# This is the scale value for the multi-resolution volumes. The original
# voxel size is divided by this value and then used as the voxel size for
# the multi-resolution volumes. Default is 4 times the original voxel size.
# This means that the multi-resolution regions have larger voxel sizes if
# this is < 1, i.e. 1/4 = 4 times the original voxel size.
options.multiResolutionScale = 1/4

CTEFOVCorrection(options)

### Flip the image (in column direction)?
options.flip_image = True

### How much is the image rotated (radians)?
# The angle (in radians) on how much the image is rotated BEFORE
# reconstruction, i.e. the rotation is performed in the detector space.
# Positive values perform the rotation in counter-clockwise direction
options.offangle = (3.*np.pi)/2.

# Computation device
options.deviceNum = 0
# Projector, see https://omega-doc.readthedocs.io/en/latest/selectingprojector.html for details
options.projector_type = 4
# Interpolation length, 1 means the length of one voxel
options.dL = 1

options.Niter = 10
# Number of subsets. If you use non-subset algorithm, you need to manually set this to 1
options.subsets = 20
options.subsetType = 8

# This has to be True if you want to use the filtering-based preconditioner
options.PDHG = True

# Use offset-correction
# If you use offset imaging, i.e. the center of rotation is not in the
# origin but rather a circle around the origin, you can enable automatic
# offset weighting by setting this to True.
options.offsetCorrection = False

# Regularization parameter
options.beta = .1
# Filter parameter for non-local regularization
options.NLMsigma = 4.00e-3
# Edge preservation parameter for RDP and NLRDP
options.RDP_gamma = 10.

### Show status messages
# These are e.g. time elapsed on various functions and what steps have been
# completed. It is recommended to keep this at 1 or 2. With value of 2, 
# you get more detailed timing information. Maximum is 3. Minimum is 0.
options.verbose = 1

options.precondTypeImage[6] = False
options.precondTypeImage[1] = False

# Filtering-based preconditioner
options.precondTypeMeas[1] = True

# At the moment this has no effect. Filtering is either always on or off
options.filteringIterations = 200

# Number of power method iterations
options.powerIterations = 10

# Needed for CT data
options.CT = True

options.addProjector()

linearizeData(options)

# Assumes that Arrayfire arrays are used as input to either forward or backward projections
options.useAF = True

# Compute forward projection with options * f
# Compute backprojection with options.T() * y

options.initProj()
# fp = options * d_f
# bp = options.T() * d_m


import arrayfire as af
d_f = af.interop.np_to_af_array(options.x0)
m = options.SinM.ravel('F')
L = powerMethod(options)
sigma = 1.
theta = 1.
tau = L
import time
tic = time.perf_counter()
    

# """
# PDHG without subsets or multi-resolution
# """
# d_m = af.interop.np_to_af_array(m)
# p = af.constant(0., m.size)
# f = af.constant(0., options.x0.size)
# for k in range(options.Niter):
#     # tic = time.perf_counter()
#     apu = options * d_f
#     if options.precondTypeMeas[1]:
#         apu -= d_m
#         af.eval(apu)
#         apu = applyMeasPreconditioning(options, apu)
#         p = (p + sigma * apu)
#         p = circulantInverse(options, p)
#     else:
#         p = (p + sigma * (apu - d_m)) / (1 + sigma)
#     apu = options.T() * p
#     fPrev = f.copy()
#     f = f - tau * apu
#     f[f < 0] = 0
#     af.eval(f)
#     d_f = f + theta * (f - fPrev)
#     af.device_gc()
#     print('Iteration ' + str(k))

"""
PDHG with subsets but without multi-resolution
"""
d_m = [None] * options.subsets
for k in range(options.subsets):
    d_m[k] = af.interop.np_to_af_array(m[options.nTotMeas[k].item() : options.nTotMeas[k + 1].item()])
p = [None] * options.subsets
for k in range(options.subsets):
    p[k] = af.constant(0., d_m[k].elements())
g1 = af.constant(0., d_f.elements())
for k in range(options.Niter):
    for i in range(options.subsets):
        options.subset = i
        apu = options * d_f
        if options.precondTypeMeas[1]:
            apu = applyMeasPreconditioning(options, apu - d_m[i])
            pl = p[i] + sigma * apu
            pl = circulantInverse(options, pl)
        else:
            pl = (p[i] + sigma * (apu - d_m[i])) / (1 + sigma)
        dg = options.T() * (pl - p[i])
        p[i] = pl
        g1 = g1 + dg
        af.eval(g1)
        g = g1 + (theta * options.subsets) * dg
		# Regularized versions are enabled by uncommenting one regularizer and the update equation
        # grad = RDP(d_f, options.Nx, options.Ny, options.Nz, options.RDP_gamma, options.beta)
        # grad = NLReg(d_f, options.Nx, options.Ny, options.Nz, options.NLMsigma, options.beta, NLType = 3, useAdaptive=True, adaptiveConstant=5e-6)
        # grad = TV(d_f, options.Nx, options.Ny, options.Nz, options.beta)
        # d_f = d_f - tau * (g + grad)
        d_f = d_f - tau * g
        d_f[d_f < 0] = 0.
        af.eval(d_f)
        af.device_gc()
        print('Sub-iteration ' + str(i))
    print('Iteration ' + str(k))
    
# """
# PDHG without subsets but with multi-resolution
# """
# d_m = af.interop.np_to_af_array(m)
# p = af.constant(0., m.size)
# f = [None] * (options.nMultiVolumes + 1)
# d_f = [None] * (options.nMultiVolumes + 1)
# NN = 0
# NN2 = 0
# for i in range((options.nMultiVolumes + 1)):
#     f[i] = af.constant(0., options.N[i].item())
#     NN += options.N[i].item()
#     d_f[i] = af.interop.np_to_af_array(options.x0[NN2 : NN])
#     NN2 += options.N[i].item()
# for k in range(options.Niter):
#     apu = options * d_f
#     if options.precondTypeMeas[1]:
#         apu -= d_m
#         af.eval(apu)
#         apu = applyMeasPreconditioning(options, apu)
#         p = (p + sigma * apu)
#         p = circulantInverse(options, p)
#     else:
#         p = (p + sigma * (apu - d_m)) / (1 + sigma)
#     af.eval(p)
#     apu = options.T() * p
#     for i in range((options.nMultiVolumes + 1)):
#         fPrev = f[i]
#         if i == 0:
#             f[i] = f[i] - tau[i] * apu[i]
#         else:
#             f[i] = f[i] - tau[6] / 8 * apu[i]
#         f[i][f[i] < 0] = 0.
#         af.eval(f[i])
#         d_f[i] = f[i] + theta * (f[i] - fPrev)
#     print('Iteration ' + str(k))

# """
# PDHG with subsets and multi-resolution
# """
# d_m = [None] * options.subsets
# for k in range(options.subsets):
#     d_m[k] = af.interop.np_to_af_array(m[options.nTotMeas[k].item() : options.nTotMeas[k + 1].item()])
# p = [None] * options.subsets
# for k in range(options.subsets):
#     p[k] = af.constant(0., d_m[k].elements())
# g1 = [None] * (options.nMultiVolumes + 1)
# d_f = [None] * (options.nMultiVolumes + 1)
# NN = 0
# NN2 = 0
# for i in range((options.nMultiVolumes + 1)):
#     g1[i] = af.constant(0., options.N[i].item())
#     NN += options.N[i].item()
#     d_f[i] = af.interop.np_to_af_array(options.x0[NN2 : NN])
#     NN2 += options.N[i].item()
# for k in range(options.Niter):
#     for i in range(options.subsets):
#         options.subset = i
#         apu = options * d_f
#         if options.precondTypeMeas[1]:
#             apu = applyMeasPreconditioning(options, apu - d_m[i])
#             pl = p[i] + sigma * apu
#             pl = circulantInverse(options, pl)
#         else:
#             pl = (p[i] + sigma * (apu - d_m[i])) / (1 + sigma)
#         dg = options.T() * (pl - p[i])
#         p[i] = pl.copy()
#         for j in range((options.nMultiVolumes + 1)):
#             g1[j] = g1[j] + dg[j]
#             af.eval(g1[j])
#             g = g1[j] + (theta * options.subsets) * dg[j]
#             if (j > 0):
#                 d_f[j] = d_f[j] - tau[6] / 8 * g
#             else:
#                 d_f[j] = d_f[j] - tau[0] * g
#             d_f[j][d_f[j] < 0] = 0.
#             af.eval(d_f[j])
#         af.sync()
#         af.device_gc()
#         print('Sub-iteration ' + str(i))
#     print('Iteration ' + str(k))

toc = time.perf_counter()
leike = 0
print(f"{toc - tic:0.4f} seconds")
if isinstance(d_f, list):
    f_np = d_f[0].to_ndarray()
else:
    f_np = d_f.to_ndarray()
f_np = np.reshape(f_np, (options.Nx[0].item(), options.Ny[0].item(), options.Nz[0].item()), order='F')



z = np.int16(f_np * 55000) - 1000

plt.pyplot.imshow(z[:,:,420], vmin=-1000, vmax=2000)

af.sync()
af.device_gc()