# -*- coding: utf-8 -*-
"""
This example shows how to do custom algorithm reconstructions with CBCT data. This uses given source and center of detector panel coordinates,
as well as additional rotation by the panel. The algorithm is always primal-dual hybrid gradient, but there are several variations of it:
with or without subsets and with or without multi-resolution reconstruction. For the multi-resolution cases, be sure the turn extended FOV
on before trying them. Likewise, turn the multi-resolution off when using non-multi-resolution cases. Default is with subsets but without
multi-resolution. Filtering-based preconditioner is on.

This example is the recommended way for OpenCL reconstructions, although a PyOpenCL version exists too.
"""
import numpy as np
from omegatomo import proj
from omegatomo.util import CTEFOVCorrection
from omegatomo.reconstruction.prepass import linearizeData
from omegatomo.util.powermethod import powerMethod
from omegatomo.util.measprecond import applyMeasPreconditioning, circulantInverse
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
options.dPitchX = var['dPitch'][0]
options.dPitchY = var['dPitch'][1]

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
options.transaxialEFOV = True
options.axialEFOV = True
options.transaxialExtrapolation = False
options.axialExtrapolation = True

options.useMultiResolutionVolumes = False

options.multiResolutionScale = 1/4

CTEFOVCorrection(options)

# Rotate/flip the image, i.e. the detector coordinates
options.flip_image = True
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

options.offsetCorrection = False
# Should have no effect at the moment
options.enforcePositivity = True

# Unused at the moment
options.beta = .1
options.NLMsigma = 2.00e-3
options.RDP_gamma = 10.

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
# # u2 = p.to_ndarray()
# # u2 = np.reshape(u2, (options.nRowsD, options.nColsD, -1), order='F')
# # plt.pyplot.imshow(u2[:,:,260])
# f_np2 = d_f.to_ndarray()
# f_np2 = np.reshape(f_np2, (options.Nx[0].item(), options.Ny[0].item(), options.Nz[0].item()), order='F')
# plt.pyplot.imshow(f_np2[:,:,200])

obj = np.zeros(options.Niter*options.subsets,dtype=np.float32)
"""
PDHG with subsets but without multi-resolution
"""
uu = 0
# FilterG = af.interop.np_to_af_array(options.filter0)
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
        d_f = d_f - tau * g
        d_f[d_f < 0] = 0.
        af.eval(d_f)
        ar = (options * d_f - d_m[i])
        var = applyMeasPreconditioning(options, ar)
        # var = af.moddims(ar, options.nRowsD, d1=options.nColsD, d2=ar.elements() // (options.nRowsD * options.nColsD));
        # temp = af.fft(var, options.Nf)
        # temp /= af.tile(FilterG, 1, d1=temp.shape[1], d2=temp.shape[2])
        # af.eval(temp)
        # af.ifft_inplace(temp)
        # var = af.flat(af.real(temp[:var.shape[0], :, :]))
        ar = af.blas.matmulTN(ar, var) * .5
        # ar = .5 * af.blas.matmulTN(ar, ar)
        obj[uu] = ar.to_ndarray()
        uu += 1
        af.device_gc()
        print('Sub-iteration ' + str(i))
    print('Iteration ' + str(k))
    
# from skimage.transform import resize #scikit-image
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
# # u1 = p.to_ndarray()
# # u1 = np.reshape(u1, (options.nRowsD, options.nColsD, -1), order='F')
# # plt.pyplot.imshow(u1[:,:,260])
# if isinstance(d_f, list):
#     f_np = np.zeros((681, 681, 401), dtype=np.float32, order='F')
#     var = 1
#     var2 = -1
#     for k in range(len(d_f)):
#         ff = d_f[k].to_ndarray()
#         ff = np.reshape(ff, (options.Nx[k].item(), options.Ny[k].item(), options.Nz[k].item()), order='F')
#         if k == 0:
#             f_np[int(options.Nx[3].item() / options.multiResolutionScale) - var : int(options.Nx[3].item() / options.multiResolutionScale) - var + options.Nx[0].item(), 
#                   int(options.Ny[5].item() / options.multiResolutionScale) - var : int(options.Ny[5].item() / options.multiResolutionScale) - var + options.Ny[0].item(), 
#                   int(options.Nz[1].item() / options.multiResolutionScale) : int(options.Nz[1].item() / options.multiResolutionScale) + options.Nz[0].item()] = ff
#         if k == 1:
#             if options.multiResolutionScale < 1:
#                 ff = resize(ff, (options.Nx[0].item(), options.Ny[0].item(), int(options.Nz[1].item() / options.multiResolutionScale)))
#             f_np[int(options.Nx[3].item() / options.multiResolutionScale) - var : int(options.Nx[3].item() / options.multiResolutionScale) - var + options.Nx[0].item(), 
#                  int(options.Ny[5].item() / options.multiResolutionScale) - var: int(options.Ny[5].item() / options.multiResolutionScale) - var + options.Ny[0].item(),  
#                  : int(options.Nz[1].item() / options.multiResolutionScale)] = ff
#         if k == 2:
#             if options.multiResolutionScale < 1:
#                 ff = resize(ff, (options.Nx[0].item(), options.Ny[0].item(), int(options.Nz[2].item() / options.multiResolutionScale)))
#             f_np[int(options.Nx[3].item() / options.multiResolutionScale) - var : int(options.Nx[3].item() / options.multiResolutionScale) - var + options.Nx[0].item(), 
#                  int(options.Ny[5].item() / options.multiResolutionScale) - var : int(options.Ny[5].item() / options.multiResolutionScale) - var + options.Ny[0].item(), 
#                  -int(options.Nz[1].item() / options.multiResolutionScale):] = ff
#         if k == 3:
#             if options.multiResolutionScale < 1:
#                 ff = resize(ff, (int(options.Nx[k].item() / options.multiResolutionScale) - var, 681, 401))
#             f_np[ : int(options.Nx[k].item() / options.multiResolutionScale) - var, :, :] = ff
#         if k == 4:
#             if options.multiResolutionScale < 1:
#                 ff = resize(ff, (int(options.Nx[k].item() / options.multiResolutionScale) - var, 681, 401))
#             f_np[int(options.Nx[k].item() / options.multiResolutionScale) - var + options.Nx[0].item() :, :, :] = ff
#         if k == 5:
#             if options.multiResolutionScale < 1:
#                 ff = resize(ff, (int(options.Nx[k].item() / options.multiResolutionScale) + var2, int(options.Ny[k].item() / options.multiResolutionScale) - var, 401))
#             f_np[int(options.Nx[3].item() / options.multiResolutionScale) - var : int(options.Nx[3].item() / options.multiResolutionScale) - var + options.Nx[0].item(), 
#                  : int(options.Ny[5].item() / options.multiResolutionScale) - var,  :] = ff
#         if k == 6:
#             if options.multiResolutionScale < 1:
#                 ff = resize(ff, (int(options.Nx[k].item() / options.multiResolutionScale) + var2, int(options.Ny[k].item() / options.multiResolutionScale) - var, 401))
#             f_np[int(options.Nx[3].item() / options.multiResolutionScale) - var : int(options.Nx[3].item() / options.multiResolutionScale) - var + options.Nx[0].item(), 
#                  int(options.Ny[5].item() / options.multiResolutionScale) - var + options.Ny[0].item() :,  :] = ff
#     plt.pyplot.imshow(f_np[:,:,305])

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
    # f_np1 = apu[leike].to_ndarray()
    # f_np1 = np.reshape(f_np1, (options.Nx[leike].item(), options.Ny[leike].item(), options.Nz[leike].item()), order='F')
    # plt.pyplot.imshow(f_np1[:,:,2])
else:
    f_np = d_f.to_ndarray()
    # f_np2 = apu.to_ndarray()
    # f_np2 = np.reshape(f_np2, (options.Nx[0].item(), options.Ny[0].item(), options.Nz[0].item()), order='F')
    # plt.pyplot.imshow(f_np2[:,:,90], vmin=-10, vmax=5)
f_np = np.reshape(f_np, (options.Nx[0].item(), options.Ny[0].item(), options.Nz[0].item()), order='F')



z = np.int16(f_np * 55000) - 1000

# apu1 = options * d_f
# u1 = apu1.to_ndarray()
# u1 = np.reshape(u1, (options.nRowsD, options.nColsD, -1), order='F')
# plt.pyplot.imshow(u1[:,:,500])
# apu = options * d_f
# u = apu.to_ndarray()
# u = np.reshape(u, (options.nRowsD, options.nColsD, -1), order='F')
# plt.pyplot.imshow(u[:,:,260], vmin=0.030643528, vmax=0.031839356)

plt.pyplot.imshow(z[:,:,420], vmin=-1000, vmax=2000)

# apu1 = options.T() * d_m
# if isinstance(apu1, list):
#     f_np = np.zeros((681, 681, 401), dtype=np.float32, order='F')
#     for k in range(len(apu1)):
#         ff = apu1[k].to_ndarray()
#         ff = np.reshape(ff, (options.Nx[k].item(), options.Ny[k].item(), options.Nz[k].item()), order='F')
#         if k == 0:
#             f_np[options.Nx[3].item() : options.Nx[3].item() + options.Nx[0].item(), options.Ny[5].item() : options.Ny[5].item() + options.Ny[0].item(), options.Nz[1].item() : options.Nz[1].item() + options.Nz[0].item()] = ff
#         if k == 1:
#             f_np[options.Nx[3].item() : options.Nx[3].item() + options.Nx[0].item(), options.Ny[5].item() : options.Ny[5].item() + options.Ny[0].item(),  : options.Nz[1].item()] = ff
#         if k == 2:
#             f_np[options.Nx[3].item() : options.Nx[3].item() + options.Nx[0].item(), options.Ny[5].item() : options.Ny[5].item() + options.Ny[0].item(), options.Nz[1].item() + options.Nz[0].item() :] = ff
#         if k == 3:
#             f_np[ : options.Nx[3].item(), :, :] = ff
#         if k == 4:
#             f_np[options.Nx[3].item() + options.Nx[0].item() :, :, :] = ff
#         if k == 5:
#             f_np[options.Nx[3].item() : options.Nx[3].item() + options.Nx[0].item(), : options.Ny[5].item(),  :] = ff
#         if k == 6:
#             f_np[options.Nx[3].item() : options.Nx[3].item() + options.Nx[0].item(), options.Ny[5].item() + options.Ny[0].item() : ,  :] = ff
#     plt.pyplot.imshow(f_np[:,:,200])
# else:
#     f_np2 = apu1.to_ndarray()
#     f_np2 = np.reshape(f_np2, (options.Nx[0].item(), options.Ny[0].item(), options.Nz[0].item()), order='F')
#     plt.pyplot.imshow(f_np2[:,:,5])
    

af.sync()
af.device_gc()