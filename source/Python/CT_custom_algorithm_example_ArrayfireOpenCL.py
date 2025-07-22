# -*- coding: utf-8 -*-
"""
## Python codes for CBCT custom algorithm reconstruction
This example contains a simplified example for custom algorithm
reconstruction using TIFF projection CBCT data. Currently the support for
some of the additional features is limited. The default configuration
uses PDHG without multi-resolution reconstruction, but with subsets. PDHG
should also work with multi-resolution reconstruction and without subsets.

Note that custom algorithm refers to your own algorithms and not the built-in 
algorithms. This example merely has the PDHG algorithm shown as an example.
The forward and/or backward projections of OMEGA are utilized for the computation
of these algorithms.

This example uses Arrayfire with PyOpenCL and thus requires OpenCl (and PyOpenCL and Arrayfire)!
"""
import numpy as np
from omegatomo.projector import proj
from omegatomo.fileio.loadProjectionImages import loadProjectionImages
from omegatomo.util import CTEFOVCorrection
from omegatomo.reconstruction.prepass import linearizeData
from omegatomo.util.powermethod import powerMethod
from omegatomo.util.measprecond import applyMeasPreconditioning, circulantInverse
import matplotlib as plt

options = proj.projectorClass()

###########################################################################
###########################################################################
###########################################################################
########################### SCANNER PROPERTIES ############################
###########################################################################
###########################################################################
###########################################################################

# This setting determines whether the high-dimensional scalable
# reconstruction is used (if set as True). Otherwise, the regular
# reconstruction is performed.
options.largeDim = False

### Binning
# The level of binning used for the raw data. For example binning of 2
# reduces the size of the projections by two from both dimensions (e.g.
# 2048x3072 becomes 1024x1536).
options.binning = 4

### Number of detector pixels (horizontal)
# The number of detector pixels in the detector panel (horizontal
# direction)
# NOTE: if you use binning, this value has to use the final binned
# dimensions
options.nColsD = 2368//options.binning

### Number of detector pixels (vertical)
# The number of detector pixels in the detector panel (vertical
# direction)
# NOTE: if you use binning, this value has to use the final binned
# dimensions
options.nRowsD = 2240//options.binning

### Number of projections
# Total number of projections used
options.nProjections = 721

### Projection angles (degree or radian)
# The angles corresponding to the projections
options.angles = -np.linspace(0, 360, options.nProjections, dtype=np.float32)

### Detector pixel pitch/size (mm)
# The size of the detector/distance between adjacent detectors
# NOTE: if you use binning, this value has to use the final binned
# dimensions
options.dPitch = 0.05*options.binning

### Source to detector distance (mm)
# The orthogonal distance from the source to the detector panel
options.sourceToDetector = 553.74

### Source to center of rotation distance (mm)
# The distance from the source to the center of rotation/object/origin
options.sourceToCRot = 210.66

### Name of current datafile/examination
# This is used for naming purposes only
options.name = 'Walnut3DCT_data'

### Compute only the reconstructions
# If this file is run with this set to True, then the data load will be
# skipped if the options.SinM variable exists
options.only_reconstructions = False

### Show status messages
# These are e.g. time elapsed on various functions and what steps have been
# completed. It is recommended to keep this 1. Maximum value of 3 is supported.
options.verbose = 1

### Transaxial FOV size (mm), this is the length of the x (horizontal) side
# of the FOV
options.FOVa_x = 40.1

### Transaxial FOV size (mm), this is the length of the y (vertical) side
# of the FOV
options.FOVa_y = options.FOVa_x

### Axial FOV (mm)
options.axial_fov = 40

### Source row offset (mm)
# The center of rotation is not exactly in the origin. With this parameter
# the source location can be offset by the specifed amount (row direction).
# This has a similar effect as circulary shifting the projection images.
# Use vector values if these vary with each projection (this is untested at
# the moment).
# NOTE: The default value has been obtained experimentally and is not based
# on any known value.
options.sourceOffsetRow = -0.16

### Source column offset (mm)
# Same as above, but for column direction.
options.sourceOffsetCol = 0

### Detector panel row offset (mm)
# Same as above, but the offset value for the detector panel.
options.detOffsetRow = 0

### Detector panel column offset (mm)
# Same as above, but for column direction.
options.detOffsetCol = 0

### Bed offset (mm)
# The offset values for multi-bed step-and-shoot examinations. Each bed
# position should have its own offset value.
options.bedOffset = np.empty(0, dtype=np.float32)

### Pitch/roll angles for the detector panel (radian)
# Sometimes the detector panel is slightly rotated in all three directions.
# The main rotation is included in the above options.angles, but the other
# two directions can be included here. pitchRoll should be column vector,
# where the first column corresponds to the rotation in the XY-plane and
# the second to rotation in the ZY-plane. 
options.pitchRoll = np.empty(0, dtype=np.float32)

### Direction vectors (normalized)
# This one is optional, but you can also input straight the direction
# vectors for all dimensions. The number of required dimensions depends on
# the axes where rotation occurs. If pitchRoll would be empty, i.e.
# rotation is only in the ZX-plane (angles) then only two dimensions are
# needed, one for X- and one for Y-direction (row direction). Z-direction
# (column direction) is handled automatically. If the other two rotations
# are included, then six dimensions are needed. The first three are the
# direction vectors for the placement in the row-direction. I.e. they are
# used to determine the current detector pixel location in the
# row-direction. The latter three are for the detector pixel location in
# the column direction. The dimensions should be such that the number of
# rows for uV should be either 2 or 6, while the number of columns should
# equal the number of projections. Note that these values have to be
# normalized values as they are later multiplied with the size of the
# detector pixel. If you input the above pitchRoll, these are computed
# automatically or if there is no rotation other than the general rotation
# of angles, then this is also generally not required.
options.uV = np.empty(0, dtype=np.float32)


###########################################################################
###########################################################################
###########################################################################
###########################################################################

# You can input the first projection below or leave it blank in which case
# you will be prompted to select the first projection image
options.fpath = '/path/to/20201111_walnut_0001.tif'

if ~options.only_reconstructions or not hasattr(options,'SinM'):
    options.SinM = loadProjectionImages(options.nProjections,options.binning,options.fpath)
    options.SinM = np.transpose(options.SinM, (1, 0, 2))
# NOTE: If you want to reduce the number of projections, you need to do
# this manually as outlined below:
# options.SinM = options.SinM(:,:,1:4:options.nProjections)
# options.angles = options.angles(1:4:numel(options.angles))
# options.nProjections = numel(options.angles)

# Flat value
# Needed for both linearized and Poisson-based data
options.flat = np.max(options.SinM)


###########################################################################
###########################################################################
###########################################################################
########################### IMAGE PROPERTIES ##############################
###########################################################################
###########################################################################
###########################################################################

# Note that non-square transaxial image sizes can be unreliable just as the
# non-square transaxial FOV, but they should, generally, work
### Reconstructed image pixel size (X-direction)
options.Nx = 280*2

### Y-direction
options.Ny = 280*2

# The above, again, doesn't apply to axial direction
### Z-direction (number of slices) (axial)
options.Nz = 280*2

### Flip the image (in vertical direction)?
options.flip_image = False

### How much is the image rotated (radians)?
# The angle (in radians) on how much the image is rotated BEFORE
# reconstruction, i.e. the rotation is performed in the detector space.
options.offangle = (2*np.pi)/2

###########################################################################
###########################################################################
###########################################################################
###########################################################################


### Use projection extrapolation
# If True, extrapolates the projection data. You can select below whether
# this extrapolation is done only in the axial or transaxial directions, or
# both. Default extrapolation length is 20# of the original length, for
# both sides. For example if axial extrapolation is enabled, then the left
# and right regions of the projection get 20# increase in size. This value
# can be adjusted in CTEFOVCorrection. The values are scaled to air with
# the use of logarithmic scaling.
options.useExtrapolation = False

### Use extended FOV
# Similar to above, but expands the FOV. The benefit of expanding the FOV
# this way is to enable to the use of multi-resolution reconstruction or
# computation of the priors/regularization only in the original FOV. The
# default extension is 40# per side.
options.useEFOV = False

# Use transaxial extended FOV (this is off by default)
options.transaxialEFOV = False

# Use axial extended FOV (this is on by default. If both this and
# transaxialEFOV are False but useEFOV is True, the axial EFOV will be
# turned on)
options.axialEFOV = False

# Same as above, but for extrapolation. Same default behavior exists.
options.transaxialExtrapolation = False

# Same as above, but for extrapolation. Same default behavior exists.
options.axialExtrapolation = False

# Setting this to True uses multi-resolution reconstruction when using
# extended FOV. Only applies to extended FOV!
options.useMultiResolutionVolumes = True

# This is the scale value for the multi-resolution volumes. The original
# voxel size is divided by this value and then used as the voxel size for
# the multi-resolution volumes. Default is 1/4 of the original voxel size.
# This means that the multi-resolution regions have smaller voxel sizes if
# this is < 1.
options.multiResolutionScale = 1/4

# Performs the extrapolation and adjusts the image size accordingly
CTEFOVCorrection(options)

# Use offset-correction
# If you use offset imaging, i.e. the center of rotation is not in the
# origin but rather a circle around the origin, you can enable automatic
# offset weighting by setting this to True.
options.offsetCorrection = False


############################### PROJECTOR #################################
### Type of projector to use for the geometric matrix
# 1 = Improved/accelerated Siddon's algorithm
# 2 = Orthogonal distance based ray tracer (not recommended in CT)
# 3 = Volume of intersection based ray tracer (not recommended in CT)
# 4 = Interpolation-based projector (ray- and voxel-based)
# 5 = Branchless distance-driven projector
# See the doc for more information:
# https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
options.projector_type = 4

### Use mask
# The mask needs to be a binary mask (uint8 or logical) where 1 means that
# the pixel is included while 0 means it is skipped. Separate masks can be
# used for both forward and backward projection and either one or both can
# be utilized at the same time. E.g. if only backprojection mask is input,
# then only the voxels which have 1 in the mask are reconstructed.
# Currently the masks need to be a 2D image that is applied identically at
# each slice.
# Forward projection mask
# If nonempty, the mask will be applied. If empty, or completely omitted, no
# mask will be considered.
# options.maskFP = True(options.nRowsD,options.nColsD)
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
# value will be multiplied by the voxel size which means that 1 means that
# the interpolation length corresponds to a single voxel (transaxial)
# length. Larger values lead to faster computation but at the cost of
# accuracy. Recommended values are between [0.5 1].
options.dL = 0.5


######################### RECONSTRUCTION SETTINGS #########################
### Number of iterations (all reconstruction methods)
options.Niter = 2

### Number of subsets (all excluding MLEM and subset_type = 5)
# Note that with high-dimensional data this is required for FDK as well.
# For high-dimensional data this controls the amount of memory required by 
# the GPU. More subsets, less memory, but using too many subsets can lead 
# to reduced performance.
options.subsets = 10

### Subset type (n = subsets)
# 1 = Every nth (column) measurement is taken
# 2 = Every nth (row) measurement is taken (e.g. if subsets = 3, then
# first subset has measurements 1, 4, 7, etc., second 2, 5, 8, etc.) 
# 3 = Measurements are selected randomly
# 4 = (Sinogram only) Take every nth column in the projection image
# 5 = (Sinogram only) Take every nth row in the projection image
# 8 = Use every nth projection image
# 9 = Randomly select the full projection images
# 11 = Use prime factor sampling to select the full projection images
# Most of the time subset_type 8 is sufficient.
options.subsetType = 8

###########################################################################

# This has to be True if you want to use the filtering-based preconditioner
options.PDHG = True

# Image-based preconditioners
# Setting options.precondTypeImage(1) = true when using PKMA, MRAMLA or
# MBSREM is recommended
# precondTypeImage(0) = Diagonal normalization preconditioner (division with
# the sensitivity image 1 / (A^T1), A is the system matrix) 
# precondTypeImage(1) = EM preconditioner (f / (A^T1), where f is the current
# estimate) 
# precondTypeImage(2) = IEM preconditioner (max(n, fhat, f)/ (A^T1), where
# fhat is an estimate of the final image and n is a small positive number) 
# precondTypeImage(3) = Momentum-like preconditioner (basically a step size
# inclusion) 
# precondTypeImage(4) = Gradient-based preconditioner (Uses the normalized
# divergence (sum of the gradient) of the current estimate) 
# precondTypeImage(5) = Filtering-based preconditioner
# precondTypeImage(6) = Curvature-based preconditioner
options.precondTypeImage[0] = False
options.precondTypeImage[1] = False
options.precondTypeImage[2] = False
options.precondTypeImage[3] = False
options.precondTypeImage[4] = False
options.precondTypeImage[5] = False
options.precondTypeImage[6] = False

# Measurement-based preconditioners
# precondTypeMeas(0) = Diagonal normalization preconditioner (1 / (A1))
# precondTypeMeas(1) = Filtering-based preconditioner
options.precondTypeMeas[0] = False
options.precondTypeMeas[1] = True

# At the moment this has no effect. Filtering is either always on or off
options.filteringIterations = 100

# Number of power method iterations
options.powerIterations = 10

# Needed for CT data
options.CT = True

### OpenCL/CUDA device used
options.deviceNum = 0

options.addProjector()

linearizeData(options)

# Assumes that Arrayfire arrays are used as input to either forward or backward projections
options.useAF = True

# Compute forward projection with options * f
# Compute backprojection with options.T() * y

options.initProj()
import arrayfire as af
d_f = af.interop.np_to_af_array(options.x0)
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
#     p = (p + sigma * (apu - d_m)) / (1 + sigma)
#     af.eval(p)
#     apu = options.T() * p
#     for i in range((options.nMultiVolumes + 1)):
#         fPrev = f[i]
#         f[i] = f[i] - tau[i] * apu[i]
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
#         pl = (p[i] + sigma * (apu - d_m[i])) / (1 + sigma)
#         dg = options.T() * (pl - p[i])
#         p[i] = pl.copy()
#         for j in range((options.nMultiVolumes + 1)):
#             g1[j] = g1[j] + dg[j]
#             af.eval(g1[j])
#             g = g1[j] + (theta * options.subsets) * dg[j]
#             d_f[j] = d_f[j] - tau[j] * g
#             d_f[j][d_f[j] < 0] = 0.
#             af.eval(d_f[j])
#         af.sync()
#         af.device_gc()
#         print('Sub-iteration ' + str(i))
#     print('Iteration ' + str(k))

toc = time.perf_counter()
print(f"{toc - tic:0.4f} seconds")
if isinstance(d_f, list):
    f_np = d_f[0].to_ndarray()
else:
    f_np = d_f.to_ndarray()
f_np = np.reshape(f_np, (options.Nx[0].item(), options.Ny[0].item(), options.Nz[0].item()), order='F')
plt.pyplot.imshow(f_np[:,:,200], vmin=0)

af.sync()
af.device_gc()