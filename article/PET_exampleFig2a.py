"""
Python codes for PET reconstruction
This example file computes the Figure 2 (a) of the OMEGA V2 article. 
DOI will be added later.
Used data available from: https://doi.org/10.5281/zenodo.17185907
"""
import numpy as np
from omegatomo.projector import proj
from omegatomo.reconstruction import reconstructions_main
from omegatomo.util.checkCUDA import checkCUDA

# Set the folder containing the above input data here:
path = ''

options = proj.projectorClass()


###########################################################################
###########################################################################
###########################################################################
########################### SCANNER PROPERTIES ############################
###########################################################################
###########################################################################
###########################################################################

### R-sectors/modules/blocks/buckets in transaxial direction
options.blocks_per_ring = (38)

### R-sectors/modules/blocks/buckets in axial direction (i.e. number of physical
### scanner/crystal rings) 
# Multiplying this with the below cryst_per_block_axial should equal the total
# number of crystal rings. 
options.linear_multip = (8)

### Number of detectors on the side of R-sector/block/module (transaxial
### direction)
# (e.g. 13 if 13x13, 20 if 20x10)
options.cryst_per_block = (20)

### Number of detectors on the side of R-sector/block/module (axial
### direction)
# (e.g. 13 if 13x13, 10 if 20x10)
options.cryst_per_block_axial = 10

### Crystal pitch/size in x- and y-directions (transaxial) (mm)
options.cr_p = 3.2

### Crystal pitch/size in z-direction (axial) (mm)
options.cr_pz = 3.2

### Ring diameter (distance between perpendicular detectors) (mm)
options.diameter = 780.

# Note that non-square transaxial FOV sizes should work, but might not work
# always. Square transaxial FOV is thus recommended.
### Transaxial FOV size (mm), this is the length of the x (vertical/row) side
# of the FOV
options.FOVa_x = (704./np.sqrt(2.))

### Transaxial FOV size (mm), this is the length of the y (horizontal/column) side
# of the FOV
options.FOVa_y = options.FOVa_x

# The above recommendation doesn't apply to axial FOV, i.e. this can be
# different from the transaxial FOV size(s). 
### Axial FOV (mm)
options.axial_fov = 261.

### Ring gaps (mm)
# Each ring is assumed to contain options.cryst_per_block_axial crystals
# Input the gap between each of these rings here, for every gap
# If there are no gaps, leave this empty or zero
# If the gap values are the same, you need to repeat the value for each gap
options.ringGaps = np.repeat(0.7142857, 7)

### Number of detectors per crystal ring (without pseudo detectors)
options.det_per_ring = options.blocks_per_ring*options.cryst_per_block

### Number of detectors per crystal ring (with pseudo detectors)
# NOTE: Vision normally has one pseudo detector per block, but it is
# omitted here by default
options.det_w_pseudo = options.blocks_per_ring*(options.cryst_per_block)

### Number of crystal rings
options.rings = options.linear_multip * options.cryst_per_block_axial

### Total number of detectors
options.detectors = options.det_per_ring*options.rings

### Scanner name
# Used for naming purposes (measurement data)
options.machine_name = 'Siemens_Vision'

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
options.Nx = 220

### Y/column-direction
options.Ny = 220

### Z-direction (number of slices) (axial)
options.Nz = 88

### Flip the image (in column direction)?
options.flip_image = True

### How much is the image rotated?
# You need to run the precompute phase again if you modify this
# NOTE: The rotation is done in the detector space (before reconstruction).
# This current setting is for scanner list-mode data or sinogram data.
# Positive values perform the rotation in clockwise direction
# The unit is the detector element, i.e. the below 390 means that the
# rotation is done by the size of 390 detector elements
options.offangle = 390.
 
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
options.span = 1

### Maximum ring difference
options.ring_difference = options.rings - 1

### Number of radial positions (views) in sinogram
# You should primarily use the same number as the scanner uses.
# However, if that information is not available you can use ndist_max
# function to determine potential values (see help ndist_max for usage).
options.Ndist = 520

### Number of angles (tangential positions) in sinogram 
# This is the final amount after possible mashing, maximum allowed is the
# number of detectors per ring/2.
options.Nang = options.det_per_ring//8

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

### If Ndist value is even, take one extra out of the negative side (+1) or
# from the positive side (-1). E.g. if Ndist = 200, then with +1 the
# interval is [-99,100] and with -1 [-100,99]. This varies from scanner to
# scanner. If you see a slight shift in the sinograms when comparing with
# the scanner sinograms then use the other option here.
options.ndist_side = 1
 
###########################################################################
###########################################################################
###########################################################################
###########################################################################
 
 
 
 
###########################################################################
###########################################################################
###########################################################################
############################# CORRECTIONS #################################
###########################################################################
###########################################################################
###########################################################################

########################### Randoms correction ############################
# If set to true, performs randoms correction during reconstruction or 
# performs precorrection, depending on the selections below
# If you are using your own data, the randoms
# data can be input either manually into options.SinDelayed or input when
# prompted (has to be stored in a mat-file beforehand!)
options.randoms_correction = True
 
######################### Attenuation correction ##########################
### Image-based attenuation correction
# Include attenuation correction from images (e.g. CT-images) (for this you
# need attenuation images of each slice correctly rotated and scaled for
# 511 keV) or from attenuation sinograms. Note that all the attenuation
# data has to be correctly scaled before reconstruction.
# You can either use the path below to input the data or manually input
# the attenuation data into options.vaimennus
options.attenuation_correction = True

### Attenuation image data file
# Specify the path and filename.
# NOTE: the attenuation data must be the only variable in the file and
# have the dimensions of the final reconstructed image.
# If no file is specified here, the user will be prompted to select one
# if options.vaimennus is empty (or does not exist)
# Alternatively, just input the attenuation data into options.vaimennus
options.attenuation_datafile = path + '/511keV_attenuation_coefficients_for_VisionDerenzoAttenuation.mat'


#################### Corrections during reconstruction ####################
# If set to True, all the corrections are performed during the
# reconstruction step, otherwise the corrections are performed to the
# sinogram/raw data before reconstruction. I.e. this can be considered as
# e.g. normalization weighted reconstruction if normalization correction is
# applied.
# NOTE: Attenuation correction is always performed during reconstruction
# regardless of the choice here.
# If you have manually precorrected the data, do not put those corrections
# to true that have already been applied! Otherwise, the data will be 
# precorrected twice. This obviously only applies when this is set to false
options.corrections_during_reconstruction = True
 
###########################################################################
###########################################################################
###########################################################################
###########################################################################





###########################################################################
###########################################################################
###########################################################################
############################# TOF PROPERTIES ##############################
###########################################################################
###########################################################################
###########################################################################

### Total number of TOF bins
options.TOF_bins = 33
# options.TOF_bins = 1

### Length of each TOF bin (s)
# The time length of each TOF bin in seconds
# This multiplied with the number of bins total the entire time frame that
# the TOF data contains. For example with 10 bins of size 400 ps all time
# differences of at most 4 ns will be included in the TOF data. The
# multiplied value should be, at most, the size of the coincidence window.
options.TOF_width = 143e-12

### TOF offset (s)
# If your TOF bins are not centered on zero (center of FOV) you can specify
# the offset value here.
options.TOF_offset = 4.3333e-12

### FWHM of the temporal noise/FWHM of the TOF data
# This parameter has two properties. The first one applies to any TOF data
# that is saved by OMEGA (GATE, Inveon/Biograph list-mode), the second only
# to GATE data.
# Firstly this specifies the FWHM of TOF data used for file naming and
# loading purposes. This value is included in the filename when data is
# imported/saved and also used when that same data is later loaded. 
# Secondly, this is the FWHM of the ADDED temporal noise to the time
# differences. If you are using GATE data and have set a custom temporal
# blurring in GATE then you should set to this zero if you wish to use the
# same temporal resolution. If no custom temporal blurring was applied then
# use this value to control the accuracy of the TOF data. For example if
# you want to have TOF data with 500 ps FWHM then set this value to
# 500e-12. 
options.TOF_noise_FWHM = 214e-12

### FWHM of the TOF data
# Applies to ALL data.
# This value specifies the TOF accuracy during the reconstruction process
# and thus can be different from above. If you are using GATE data with
# temporal blurring, you need to multiply that FWHM with sqrt(2) here.
options.TOF_FWHM = 214e-12

### Number of TOF bins used in reconstruction
# Number of TOF bins used during reconstruction phase.
# NOTE: Currently supports only either all bins specified by
# options.TOF_bins or 1 (or 0) which converts the TOF data into non-TOF
# data during reconstruction phase.
options.TOF_bins_used = options.TOF_bins
 
###########################################################################
###########################################################################
###########################################################################
###########################################################################
 
 
 
 
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
options.name = 'Vision_derenzo'

### Show status messages
# These are e.g. time elapsed on various functions and what steps have been
# completed. It is recommended to keep this at 1 or 2. With value of 2, 
# you get more detailed timing information. Maximum is 3, minimum 0.
options.verbose = 2
 
###########################################################################
###########################################################################
###########################################################################
###########################################################################
 
 
 
 
###########################################################################
###########################################################################
###########################################################################
######################## RECONSTRUCTION PROPERTIES ########################
###########################################################################
###########################################################################
###########################################################################
### OpenCL/CUDA device used 
# Uncomment the below lines and run them to determine the available device
# numbers:
# from omegatomo.util.devinfo import deviceInfo
# deviceInfo(True)
options.deviceNum = 0

### Use CUDA
# Selecting this to True will use CUDA kernels/code instead of OpenCL. This
# only works if the CUDA code was successfully built. This is recommended
# if you have CUDA-capable device.
options.useCUDA = checkCUDA(options.deviceNum)

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
if not options.useCUDA:
    options.use_32bit_atomics = True
    
 
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
# See the documentation for more information:
# https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
options.projector_type = 1

### Interpolation length (projector type = 4 only)
# This specifies the length after which the interpolation takes place. This
# value will be multiplied by the voxel size which means that 1 means that
# the interpolation length corresponds to a single voxel (transaxial)
# length. Larger values lead to faster computation but at the cost of
# accuracy. Recommended values are between [0.5 1].
options.dL = 1

### Use point spread function (PSF) blurring
# Applies PSF blurring through convolution to the image space. This is the
# same as multiplying the geometric matrix with an image blurring matrix.
options.use_psf = True

# FWHM (mm) of the Gaussian used in PSF blurring in all three dimensions
options.FWHM = np.array([options.cr_p, options.cr_p, options.cr_pz])
 
######################### RECONSTRUCTION SETTINGS #########################
### Number of iterations (all reconstruction methods)
options.Niter = 1

### Save specific intermediate iterations
# You can specify the intermediate iterations you wish to save here. Note
# that this uses zero-based indexing, i.e. 0 is the first iteration (not
# the initial value). By default only the last iteration is saved. Only
# full iterations (epochs) can be saved.
options.saveNIter = np.empty(0, dtype=np.float32)
# Alternatively you can save ALL intermediate iterations by setting the
# below to True and uncommenting it. As above, only full iterations
# (epochs) are saved.
# options.save_iter = False

### Number of subsets
options.subsets = 30

### Subset type (n = subsets)
# 1 = Every nth (column) measurement is taken
# 2 = Every nth (row) measurement is taken (e.g. if subsets = 3, then
# first subset has measurements 1, 4, 7, etc., second 2, 5, 8, etc.) 
# 3 = Measurements are selected randomly
# 4 = (Sinogram only) Take every nth column in the sinogram
# 5 = (Sinogram only) Take every nth row in the sinogram
# 8 = Use every nth sinogram
# 9 = Randomly select the full sinograms
# 11 = Use prime factor sampling to select the full sinograms
# Most of the time subsetType 1 or 4 is sufficient.
options.subsetType = 1

### Initial value for the reconstruction
options.x0 = np.ones((options.Nx, options.Ny, options.Nz), dtype=np.float32)
 

###########################################################################
###########################################################################
###########################################################################
######################## RECONSTRUCTION ALGORITHMS ########################
###########################################################################
###########################################################################
###########################################################################
# Reconstruction algorithms to use (choose only one algorithm and
# optionally one prior)
 
############################### ML-METHODS ################################
### Ordered Subsets Expectation Maximization (OSEM) OR Maximum-Likelihood
### Expectation Maximization (MLEM) (if subsets = 1)
options.OSEM = True

###########################################################################
###########################################################################
###########################################################################
###########################################################################

options.epps = 1e-3

# Loads the measurement data
sinogram = np.fromfile(path + '/SiemensVision_DerenzoPhantom_TOF214psFWHM_33bins_sinogram_520x95x6400_span1.bin', dtype=np.uint8)
options.SinM = sinogram.reshape((options.Ndist, options.Nang, options.rings**2, options.TOF_bins), order='F')
del sinogram
# Load randoms correction data
fpath = path + '/SiemensVision_DerenzoPhantom_TOF214psFWHM_33bins_listmode.mat'
from pymatreader import read_mat
var = read_mat(fpath, ['SinDelayed'])
options.SinDelayed = var['SinDelayed']

del var
    
# Only load the current subset into the device (GPU) memory
options.loadTOF = False

## Reconstructions
import time
tic = time.perf_counter()
pz, fp = reconstructions_main(options)
toc = time.perf_counter()
t = toc - tic
print(f"Reconstruction process took {t:0.4f} seconds")

columns_in_image, rows_in_image = np.meshgrid(np.arange(1, options.Nx + 1), np.arange(1, options.Ny + 1))
centerX = options.Nx / 2
centerY = options.Ny / 2
radius = options.Nx / 2
mask = ((rows_in_image - centerY)**2 + (columns_in_image - centerX)**2 <= radius**2)
mask = np.dstack([mask]*pz.shape[2])

pz[~mask] = 0.

import matplotlib as plt

plt.pyplot.imshow(pz[:,:,44])