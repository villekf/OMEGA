# -*- coding: utf-8 -*-
"""
## Python code for GATE PET reconstruction using ROOT or precomputed input
# Note that this file contains much less adjustable parameters than the other 
# examples. All omitted parameters will use default values. For the list of all
# adjustable parameters see main_PET_full.py file.

# You can use the same example data as with MATLAB/Octave version. mat-files
# are also supported: 
"""
import numpy as np
from omegatomo.projector import proj
from omegatomo.reconstruction import reconstructions_main

options = proj.projectorClass()

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
# Multiplying this with the below cryst_per_block should equal the total
# number of crystal rings. options.
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

### Transaxial FOV size (mm), this is the length of the x (horizontal) side
# of the FOV
options.FOVa_x = 151

### Transaxial FOV size (mm), this is the length of the y (vertical) side
# of the FOV
options.FOVa_y = options.FOVa_x

### Axial FOV (mm)
options.axial_fov = np.floor(76.8 - options.cr_pz/10)

### Number of pseudo rings between physical rings (use 0 or [] if none)
options.pseudot = 0

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
 
###########################################################################
###########################################################################
###########################################################################
###########################################################################
 
 
 

###########################################################################
###########################################################################
###########################################################################
######################## ROOT DATA FORMAT SETTINGS ########################
###########################################################################
###########################################################################
###########################################################################
 
### Is ROOT data loaded
# If True, loads the ROOT data
# If False, will use saved preloaded data (such as mat or npz)
options.use_root = True
 
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
 
### Reconstructed image pixel count (X-direction)
# NOTE: Non-square image sizes (X- and Y-direction) may not work
options.Nx = 128

### Y-direction
options.Ny = 128

### Z-direction (number of slices) (axial)
options.Nz = options.rings*2-1

### Flip the image (in vertical direction)?
options.flip_image = False

### How much is the image rotated?
# You need to run the precompute phase again if you modify this
# NOTE: The rotation is done in the detector space (before reconstruction).
# This current setting is for systems whose detector blocks start from the
# right hand side when viewing the device from front.
# Positive values perform the rotation in clockwise direction
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
# You should primarily use the same number as the device uses.
# However, if that information is not available you can use ndist_max
# function to determine potential values (see help ndist_max for usage).
options.Ndist = 200

### Number of angles (tangential positions) in sinogram 
# This is the final amount after possible mashing, maximum allowed is the
# number of detectors per ring/2.
options.Nang = options.det_per_ring//2

### Specify the amount of sinograms contained on each segment
# (this should total the total number of sinograms).
# Currently this is computed automatically, but you can also manually
# specify the segment sizes.
options.segment_table = np.concatenate((np.array(options.rings*2-1,ndmin=1), np.arange(options.rings*2-1 - (options.span + 1), max(options.Nz - options.ring_difference*2, options.rings - options.ring_difference), -options.span*2)))
options.segment_table = np.insert(np.repeat(options.segment_table[1:], 2), 0, options.segment_table[0])

### Total number of sinograms
options.TotSinos = np.sum(options.segment_table)
 
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
# If set to True, stores the delayed coincidences during data load and
# later corrects for randoms during the data formation/load or during
# reconstruction. Delayes need to be stored in GATE data for this to work.
options.randoms_correction = False
 
############################ Scatter correction ###########################
# If set to True, will prompt the user to load the scatter sinogram/raw
# data. Corrects for scatter during data formation/load or during
# reconstruction.
# NOTE: Scatter data is not created by this software and as such must be
# provided by the user. Previously created scatter sinogram/raw data matrix
# can be used though.
options.scatter_correction = False
 
######################### Attenuation correction ##########################
### Image-based attenuation correction
# Include attenuation correction from images (e.g. CT-images) (for this you
# need attenuation images of each slice correctly rotated and scaled for
# 511 keV). For CT-images you can use attenuationCT_to_511 or
# create_atten_matrix_CT functions.
options.attenuation_correction = True

### Rotate the attenuation image before correction
# Rotates the attenuation image N * 90 degrees where N is the number
# specified below. Positive values are clocwise, negative
# counter-clockwise.
options.rotateAttImage = 1

### Attenuation image data file
# Specify the path (if not in MATLAB path) and filename.
# NOTE: the attenuation data must be the only variable in the file and
# have the dimensions of the final reconstructed image.
# If no file is specified here, the user will be prompted to select one
options.attenuation_datafile = '/path/to/cylpet_example_atn1-MuMap.mhd'
 
######################## Normalization correction #########################
### Apply normalization correction
# If set to True, normalization correction is applied to either data
# formation or in the image reconstruction by using precomputed 
# normalization coefficients. I.e. once you have computed the normalization
# coefficients, turn above compute_normalization to False and set this to
# True.
options.normalization_correction = True

### Use user-made normalization
# Use either a .mat or .nrm file containing the normalization coefficients
# for normalization correction if normalization_correction is also set to
# True. 
# User will be prompted for the location of the file either during sinogram
# formation or before image reconstruction (see below).
# NOTE: If you have previously computed normalization coefficients with
# OMEGA, you do not need to set this to True. The normalization
# coefficients for the specified scanner will be automatically loaded. Use
# this only if you want to use normalization coefficients computed outside
# of OMEGA.
options.use_user_normalization = False
 
#################### Corrections during reconstruction ####################
# If set to True, all the corrections are performed during the
# reconstruction step, otherwise the corrections are performed to the
# sinogram/raw data before reconstruction. I.e. this can be considered as
# e.g. normalization weighted reconstruction if normalization correction is
# applied.
# NOTE: Attenuation correction is always performed during reconstruction
# regardless of the choice here.
options.corrections_during_reconstruction = True
 
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
options.name = 'cylpet_example_new'

### Folder for the data (.dat ASCII, .ccs LMF, .root ROOT) files
# If no files are located in the path provided below, then the current
# folder is also checked. If no files are detected there either, an error
# is thrown.
options.fpath = '/path/to/'

### Show status messages
# These are e.g. time elapsed on various functions and what steps have been
# completed. It is recommended to keep this 1.
options.verbose = 1
 
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
 
############################# IMPLEMENTATIONS #############################
### Device used 
# Uncomment the below lines and run them to determine the available device
# numbers:
# from omegatomo.util.devinfo import deviceInfo
# deviceInfo(True)
options.deviceNum = 0
 
############################### PROJECTOR #################################
### Type of projector to use for the geometric matrix
# 0 = Regular Siddon's algorithm (only available with implementation 1 and
# when precomputed_lor = False) NOT RECOMMENDED.
# 1 = Improved/accelerated Siddon's algorithm
# 2 = Orthogonal distance based ray tracer
# 3 = Volume of intersection based ray tracer
# See the docs for more information:
# https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
options.projector_type = 1

### Use point spread function (PSF) blurring
# Applies PSF blurring through convolution to the image space. This is the
# same as multiplying the geometric matrix with an image blurring matrix.
options.use_psf = False

# FWHM of the Gaussian used in PSF blurring in all three dimensions
# options.FWHM = [options.cr_p options.cr_p options.cr_pz]
options.FWHM = np.array([options.cr_p, options.cr_p, options.cr_pz])
 
######################### RECONSTRUCTION SETTINGS #########################
### Number of iterations (all reconstruction methods)
options.Niter = 4

### Number of subsets
options.subsets = 8

###########################################################################
 
 
 

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
# Supported by all implementations
options.OSEM = True

### Row-Action Maximum Likelihood Algorithm (RAMLA)
# Supported by implementations 1, 2, 4, and 5
options.RAMLA = False

### Accelerated COSEM (ACOSEM)
# Supported by implementations 1, 2, 4, and 5
options.ACOSEM = False
 
 
############################### MAP-METHODS ###############################
# Any algorithm selected here will utilize any of the priors selected below
# this. Note that only one algorithm and prior combination is allowed! You
# can also use most of these algorithms without priors (such as PKMA or
# PDHG).
### Modified BSREM (MBSREM)
# Supported by implementations 1, 2, 4, and 5
options.MBSREM = False

### Block Sequential Regularized Expectation Maximization (BSREM)
# Supported by implementations 1, 2, 4, and 5
options.BSREM = False

### Preconditioner Krasnoselskii-Mann algorithm (PKMA)
# Supported by implementations 1, 2, 4, and 5
options.PKMA = False

### Primal-dual hybrid gradient (PDHG) with KUllback-Leibler minimization
# Supported by implementations 1, 2, 4, and 5
options.PDHGKL = False


 
 
################################# PRIORS ##################################
### Median Root Prior (MRP)
options.MRP = False

### Quadratic Prior (QP)
options.quad = False

### Huber Prior (QP)
options.Huber = False

### Non-local Means (NLM) prior
options.NLM = False

### Relative difference prior
options.RDP = False
 
 
############################ ACOSEM PROPERTIES ############################
### Acceleration parameter for ACOSEM (1 equals COSEM)
options.h = 2


######################### REGULARIZATION PARAMETER ########################
### The regularization parameter for ALL regularization methods (priors)
options.beta = 1
 
 
######################### NEIGHBORHOOD PROPERTIES #########################
### How many neighboring pixels are considered 
# With MRP, QP, L, FMH, NLM, GGMRF and weighted mean
# E.g. if Ndx = 1, Ndy = 1, Ndz = 0, then you have 3x3 square area where
# the pixels are taken into account (I.e. (Ndx*2+1)x(Ndy*2+1)x(Ndz*2+1)
# area).
# NOTE: Currently Ndx and Ndy must be identical.
# For NLM this is often called the "search window".
options.Ndx = 1
options.Ndy = 1
options.Ndz = 1
 
 
############################## HP PROPERTIES ##############################
### Delta parameter for Huber prior
# Upper and lower bounds for the prior
options.huber_delta = 5
 
 
############################## NLM PROPERTIES #############################
### Filter parameter
options.sigma = 10

### Patch radius
options.Nlx = 1
options.Nly = 1
options.Nlz = 1

### Standard deviation of the Gaussian filter
options.NLM_gauss = 1


############################## RDP PROPERTIES #############################
### Edge weighting factor
options.RDP_gamma = 10
 
###########################################################################
###########################################################################
###########################################################################
###########################################################################


# Load ROOT data
if options.use_root:
    from omegatomo.fileio import loadROOT
    # Sino = uncorrected sinogram
    # SinoT = Trues sinogram
    # SinoC = Scattered events sinogram
    # SinoR = Random event sinogram (these are the true randoms)
    # SinoD = Delayed coincidences
    # Fcoord = Coordinates for each event
    # FDcoord = Coordinates for each delayed event
    Sino, SinoT, SinoC, SinoR, SinoD, Fcoord, FDcoord = loadROOT(options)
    if options.reconstruct_trues:
        options.SinM = SinoT
    else:
        options.SinM = Sino
    np.save(options.machine_name + options.name +  ' _sinograms_combined_static_' + str(options.Ndist) + str(options.Nang) + '_' + str(options.NSinos) + '_span' + str(options.span), Sino)
    

# pz, the reconstructed image
# fp, the forward projections (if selected)
pz, fp = reconstructions_main(options)

import matplotlib as plt

plt.pyplot.imshow(pz[:,:,20])
