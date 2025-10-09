# -*- coding: utf-8 -*-
"""
## Python codes for PET custom algorithm reconstruction
This example contains a simplified example for custom algorithm
reconstruction using sinogram PET data. Currently the support for
some of the additional features is limited. The default configuration
uses OSEM, but MLEM without subsets is also available.

Note that custom algorithm refers to your own algorithms and not the built-in 
algorithms. This example merely has the MLEM/OSEM algorithm shown as an example.
The forward and/or backward projections of OMEGA are utilized for the computation
of these algorithms.

This example uses Arrayfire with PyOpenCL and thus requires OpenCl (and PyOpenCL and Arrayfire)!
"""
import numpy as np
from omegatomo.projector import proj
import matplotlib as plt

# Note that the name can be anything
options = proj.projectorClass()

# Load input data
fpath = 'Cylindrical_PET_example_cylpet_example_new_sinograms_combined_static_200x168x703_span3.mat'
from pymatreader import read_mat
var = read_mat(fpath)

# Measurement data needs to be in options.SinM (or rather in the SinM of 
# whatever you name the proj.projectorClass())
options.SinM = var['raw_SinM']

###########################################################################
###########################################################################
###########################################################################
########################### SCANNER PROPERTIES ############################
###########################################################################
###########################################################################
###########################################################################
 
### R-sectors/modules/blocks/buckets in transaxial direction
options.blocks_per_ring = (42)

### R-sectors/modules/blocks/buckets in axial direction (i.e. number of physical
### scanner/crystal rings) 
# Multiplying this with the below cryst_per_block_axial should equal the total
# number of crystal rings. 
options.linear_multip = (4)

### R-sectors/modules/blocks/buckets in transaxial direction
# Required only if larger than 1
options.transaxial_multip = 1

### Number of detectors on the side of R-sector/block/module (transaxial
### direction)
# (e.g. 13 if 13x13, 20 if 20x10)
options.cryst_per_block = (8)

### Number of detectors on the side of R-sector/block/module (axial
### direction)
# (e.g. 13 if 13x13, 10 if 20x10)
options.cryst_per_block_axial = 8

### Crystal pitch/size in x- and y-directions (transaxial) (mm)
options.cr_p = 2.4

### Crystal pitch/size in z-direction (axial) (mm)
options.cr_pz = 2.4

### Ring diameter (distance between perpendicular detectors) (mm)
options.diameter = 130*2

# Note that non-square transaxial FOV sizes should work, but might not work
# always. Square transaxial FOV is thus recommended.
### Transaxial FOV size (mm), this is the length of the x (vertical/row) side
# of the FOV
options.FOVa_x = 151

### Transaxial FOV size (mm), this is the length of the y (horizontal/column) side
# of the FOV
options.FOVa_y = options.FOVa_x

# The above recommendation doesn't apply to axial FOV, i.e. this can be
# different from the transaxial FOV size(s).
### Axial FOV (mm)
options.axial_fov = np.floor(76.8 - options.cr_pz/10)

### Number of pseudo rings between physical rings (use 0 or np.empty(0, dtype=np.float32) if none)
options.pseudot = np.empty(0, dtype=np.float32)

### Ring gaps (mm)
# Each ring is assumed to contain options.cryst_per_block_axial crystals
# Input the gap between each of these rings here, for every gap
# If there are no gaps, leave this empty or zero
# If the gap values are the same, you need to repeat the value for each gap
options.ringGaps = np.empty(0, dtype=np.float32)

### Number of detectors per crystal ring (without pseudo detectors)
options.det_per_ring = options.blocks_per_ring * options.cryst_per_block * options.transaxial_multip

### Number of detectors per crystal ring (with pseudo detectors)
# If your scanner has a single pseudo detector on each (transaxial) side of
# the crystal block then simply add +1 inside the parenthesis (or uncomment
# the one below).
options.det_w_pseudo = options.blocks_per_ring*(options.cryst_per_block)
#options.det_w_pseudo = options.blocks_per_ring*(options.cryst_per_block + 1)

### Number of crystal rings
options.rings = options.linear_multip * options.cryst_per_block_axial

### Total number of detectors
options.detectors = options.det_per_ring*options.rings

### Scanner name
# Used for naming purposes (measurement data)
options.machine_name = 'Cylindrical_PET_example'

# Note: Origin is assumed to be at the center. If this is not the case, you
# can shift it with options.oOffsetX, options.oOffsetY and options.oOffsetZ
# That is row, column and slice directions
# options.oOffsetX = 0;
# options.oOffsetY = 0;
# options.oOffsetZ = 0;
 
###########################################################################
###########################################################################
###########################################################################
###########################################################################


###########################################################################
###########################################################################
###########################################################################
########################### IMAGE PROPERTIES ##############################
###########################################################################
###########################################################################
###########################################################################

# Note that non-square transaxial image sizes can be unreliable just as the
# non-square transaxial FOV, but they should, generally, work
### Reconstructed image pixel count (X/row-direction)
options.Nx = 128

### Y/column-direction
options.Ny = 128

### Z-direction (number of slices) (axial)
options.Nz = options.rings*2-1

### Flip the image (in column direction)?
options.flip_image = False

### How much is the image rotated?
# NOTE: The rotation is done in the detector space (before reconstruction).
# This current setting is for systems whose detector blocks start from the
# right hand side when viewing the scanner from front.
# Positive values perform the rotation in clockwise direction
# The units are crystals, i.e. if the value is 1, the rotation is done by
# rotating the coordinates equaling to one crystal pitch		
options.offangle = options.det_w_pseudo * (3/4)
 
###########################################################################
###########################################################################
###########################################################################
###########################################################################
 
 
 

###########################################################################
###########################################################################
###########################################################################
########################### SINOGRAM PROPERTIES ###########################
###########################################################################
###########################################################################
###########################################################################
 
### Span factor/axial compression
options.span = 3

### Maximum ring difference
options.ring_difference = options.rings - 1

### Number of radial positions (views) in sinogram
# You should primarily use the same number as the scanner uses.
# However, if that information is not available you can use ndist_max
# function to determine potential values (see help ndist_max for usage).
# This is the ROW dimension, i.e. the number of rows in the sinogram
options.Ndist = 200

### Number of angles (tangential positions) in sinogram 
# This is the final amount after possible mashing, maximum allowed is the
# number of detectors per ring/2.
# This is the COLUMN dimension, i.e. the number of columns in the sinogram
options.Nang = options.det_per_ring//2

### Specify the amount of sinograms contained on each segment
# (this should total the total number of sinograms).
# Currently this is computed automatically, but you can also manually
# specify the segment sizes.
options.segment_table = np.concatenate((np.array(options.rings*2-1,ndmin=1), np.arange(options.rings*2-1 - (options.span + 1), max(options.Nz - options.ring_difference*2, options.rings - options.ring_difference), -options.span*2)))
options.segment_table = np.insert(np.repeat(options.segment_table[1:], 2), 0, options.segment_table[0])

### Total number of sinograms
options.TotSinos = np.sum(options.segment_table)

### Number of sinograms used in reconstruction
# The first NSinos sinograms will be used for the image reconstruction.
options.NSinos = options.TotSinos
 
###########################################################################
###########################################################################
###########################################################################
###########################################################################


 
######################## Normalization correction #########################
### Apply normalization correction
# If set to True, normalization correction is applied to either data
# formation or in the image reconstruction by using precomputed 
# normalization coefficients.
# Insert the normalization data into options.normalization or insert it
# when you are prompted for it
options.normalization_correction = True


######################### Attenuation correction ##########################
### Image-based attenuation correction
# Include attenuation correction from images (e.g. CT-images) (for this you
# need attenuation images of each slice correctly rotated and scaled for
# 511 keV). 
options.attenuation_correction = True

### Rotate the attenuation image before correction
# Rotates the attenuation image N * 90 degrees where N is the number
# specified below. Positive values are clockwise, negative
# counter-clockwise.
options.rotateAttImage = 1

### Attenuation image data file
# Specify the path and filename.
# NOTE: the attenuation data must be the only variable in the file and
# have the dimensions of the final reconstructed image.
# If no file is specified here, the user will be prompted to select one
# Alternatively, input the attenuation data into options.vaimennus
# NOTE: For GATE data, the MuMap actor output can be used here
options.attenuation_datafile = '/path/to/cylpet_example_atn1-MuMap.mhd'



###########################################################################
###########################################################################
###########################################################################
######################## RECONSTRUCTION PROPERTIES ########################
###########################################################################
###########################################################################
###########################################################################
 
############################# IMPLEMENTATIONS #############################
### Device used 
# Uncomment the below lines and run them to determine the available device
# numbers:
# from omegatomo.util.devinfo import deviceInfo
# deviceInfo(True)
options.deviceNum = 0

### Use 64-bit integer atomic functions
# If True, then 64-bit integer atomic functions (atomic add) will be used
# if they are supported by the selected device.
# Setting this to True will make computations faster on GPUs that support
# the functions, but might make results slightly less reliable due to
# floating point rounding. Recommended for GPUs.
# Note: This should be used only with OpenCL
options.use_64bit_atomics = True

### Use 32-bit integer atomic functions
# If True, then 32-bit integer atomic functions (atomic add) will be used.
# This is even faster than the above 64-bit atomics version, but will also
# have significantly higher reduction in numerical/floating point accuracy.
# This should be about 20-30# faster than the above 64-bit version, but
# might lead to integer overflow if you have a high count measurement
# (thousands of coincidences per sinogram bin). Use this only if speed is
# of utmost importance. 32-bit atomics take precedence over 64-bit ones,
# i.e. if options.use_32bit_atomics = true then the 64-bit version will be 
# always set as false.
options.use_32bit_atomics = False
 
############################### PROJECTOR #################################
### Type of projector to use for the geometric matrix
# 1 = Improved/accelerated Siddon's algorithm
# 2 = Orthogonal distance based ray tracer
# 3 = Volume of intersection based ray tracer
# 4 = Interpolation-based projector
# NOTE: You can mix and match most of the projectors. I.e. 41 will use
# interpolation-based projector for forward projection while improved
# Siddon is used for backprojection.
# NOTE 2: The below additional options apply also in hybrid cases as long
# as the other projector is the corresponding projector.
# See the doc for more information:
# https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
options.projector_type = 1

### Use point spread function (PSF) blurring
# Applies PSF blurring through convolution to the image space. This is the
# same as multiplying the geometric matrix with an image blurring matrix.
options.use_psf = True

# FWHM (mm) of the Gaussian used in PSF blurring in all three dimensions
options.FWHM = np.array([options.cr_p, options.cr_p, options.cr_pz])

# Orthogonal ray tracer (projector_type = 2) only
### The 2D (XY) width (mm) of the "strip/tube" where the orthogonal distances are
# included. If tube_width_z below is non-zero, then this value is ignored.
options.tube_width_xy = options.cr_p

# Orthogonal ray tracer (projector_type = 2) only
### The 3D (Z) width (mm) of the "tube" where the orthogonal distances are
# included. If set to 0, then the 2D orthogonal ray tracer is used. If this
# value is non-zero then the above value is IGNORED.
# If you want the projector to be a tube, use this, if you want it to be 
# strip, use the above
# This slows down the reconstruction, but makes it more accurate
options.tube_width_z = options.cr_pz

# Volume ray tracer (projector_type = 3) only
### Radius (mm) of the tube-of-response (cylinder)
# The radius of the cylinder that approximates the tube-of-response.
options.tube_radius = np.sqrt(2) * (options.cr_pz / 2)

# Volume ray tracer (projector_type = 3 only)
### Relative size of the voxel (sphere)
# In volume ray tracer, the voxels are modeled as spheres. This value
# specifies the relative radius of the sphere such that with 1 the sphere
# is just large enoough to encompass an entire cubic voxel, i.e. the
# corners of the cubic voxel intersect with the sphere shell. Larger values
# create larger spheres, while smaller values create smaller spheres.
options.voxel_radius = 1

# projector_type = 1 and 4 only
### Number of rays
# Number of rays used per detector if projector_type = 1 (i.e. Improved
# Siddon is used) or projector_type = 4 (interpolation).
# The total number of rays per detector is the multiplication of the two
# below values!
# Number of rays in transaxial (row) direction
options.n_rays_transaxial = 1;
# Number of rays in axial (column) direction
options.n_rays_axial = 1;

### Interpolation length (projector type = 4 only)
# This specifies the length after which the interpolation takes place. This
# value will be multiplied by the voxel size which means that 1 means that
# the interpolation length corresponds to a single voxel (transaxial)
# length. Larger values lead to faster computation but at the cost of
# accuracy. Recommended values are between [0.5 1].
options.dL = 0.5
 
######################### RECONSTRUCTION SETTINGS #########################
### Number of iterations (all reconstruction methods)
options.Niter = 10

### Number of subsets
options.subsets = 8

### Subset type (n = subsets)
# 1 = Every nth (column) measurement is taken
# 2 = Every nth (row) measurement is taken (e.g. if subsets = 3, then
# first subset has measurements 1, 4, 7, etc., second 2, 5, 8, etc.) 
# 3 = Measurements are selected randomly
# 4 = (Sinogram only) Take every nth column in the sinogram
# 5 = (Sinogram only) Take every nth row in the sinogram
options.subsetType = 1



options.addProjector()

# Assumes that Arrayfire arrays are used as input to either forward or backward projections
options.useAF = True

# Compute forward projection with options * f
# Compute backprojection with options.T() * y

options.initProj()
import arrayfire as af
d_f = af.interop.np_to_af_array(options.x0)
m = options.SinM.ravel('F')

# if you use MLEM or similar non-subset algorithms, be sure to put subsets to 1

# """
# MLEM
# """
# d_m = af.interop.np_to_af_array(m)
# Sens = options.T() * af.constant(1, d_m.elements())
# for it in range(options.Niter):
#     fp = options * d_f
#     bp = options.T() * (d_m / fp)
#     d_f = d_f / Sens * bp
    
    

"""
OSEM
"""
d_m = [None] * options.subsets
for k in range(options.subsets):
    d_m[k] = af.interop.np_to_af_array(m[options.nTotMeas[k].item() : options.nTotMeas[k + 1].item()])
for it in range(options.Niter):
    for k in range(options.subsets):
        # This is necessary when using subsets
        # Alternative, call options.forwardProject(d_f, k) to use forward projection
        # options.backwardProject(m, k) for backprojection
        options.subset = k
        fp = options * d_f
        Sens = options.T() * af.constant(1, d_m[k].elements())
        Sens[Sens <= 0] = options.epps
        bp = options.T() * (d_m[k] / fp)
        d_f = d_f / Sens * bp
        af.eval(d_f)
    
f_np = d_f.to_ndarray()
f_np = np.reshape(f_np, (options.Nx[0].item(), options.Ny[0].item(), options.Nz[0].item()), order='F')
plt.pyplot.imshow(f_np[:,:,20], vmin=0)

af.sync()
af.device_gc()