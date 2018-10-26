%% MATLAB codes for GATE PET reconstruction using ASCII or LMF output

clear

%%%%%%%%%%%%%%%%%%%% Specify the below machine properties %%%%%%%%%%%%%%%%%%%%
% Number of real rings (possible pseudo rings are not included here)
options.rings = (32);
% R-sectors in transaxial direction
options.blocks_per_ring = (42);
% R-sectors/modules in axial directions
options.linear_multip = (4);
% number of detectors on the side of R-sector(block) (e.g. 13 if 13x13)
options.cryst_per_block = (8);
% crystal pitch in x- and y-directions (mm)
options.cr_p = 2.4;
% crystal pitch in z-direction (mm)
options.cr_pz = 2.4;
% diameter of the bore (distance between perpendicular detectors) in mm
options.diameter = 130*2;
% Transaxial FOV size (mm), assuming square FOV (this is the length of one
% side of the square FOV)
options.FOVa = 151;
% Axial FOV (mm)
options.axial_fov = floor(76.8 - options.cr_pz/2);
% Ring numbers of pseudo rings (use empty vector if none), use MATLAB
% numbering (starting from 1)
options.pseudot = [];
% Number of detectors per ring (without pseudo detectors)
options.det_per_ring = options.blocks_per_ring*options.cryst_per_block;
% Number of detectors per ring (with pseudo detectors)
options.det_w_pseudo = options.blocks_per_ring*(options.cryst_per_block);
% number of detectors
options.detectors = options.det_per_ring*options.rings;
% machine name
options.machine_name = 'Cylindrical_PET_example';





%%%%%%%%%%%%%%%%%%%% ASCII data format settings %%%%%%%%%%%%%%%%%%%%
% Is ASCII data used? (Only one data type can be used at a time)
options.use_ASCII = true;
% Columns corresponding to the R-sectors (or blocks) (see GATE wiki on ASCII file
% format)
% Corresponds to (from GATE Wiki): Column 12 to 17 : volume IDs* (first single)
options.rsector_ind1 = (6);
options.rsector_ind2 = (16);
% Module indices (columns)
% if no modules present use 0 (i.e. ECAT geometry)
options.module_ind1 = (7);
options.module_ind2 = (17);
% Crystal indices (columns)
options.crs_ind1 = (9);
options.crs_ind2 = (19);
% Time index (column), use either the first single or second
% Corresponds to (from GATE Wiki): Column 7 : Time stamp (first single)
options.time_index = (4);
% index of the source location of the first single (X-dimension, first column of the three)
% Corresponds to (from GATE Wiki): Column 4 to 6 : XYZ position of the source in world referential (first single)
% Use 0 or [] to skip this step
options.source_index1 = 1;
% index of the source location of the second single (X-dimension, first column of the three)
% Corresponds to (from GATE Wiki): Column 27 to 29 : XYZ position of the source in world referential (second single)
% Use 0 or [] to skip this step
options.source_index2 = 11;





%%%%%%%%%%%%%%%%%%%% LFM data format settings %%%%%%%%%%%%%%%%%%%%
% Is LMF data used? (Only one data type can be used at a time)
options.use_LMF = false;
% How many bytes are in the LMF header part?
options.header_bytes = (16);
% How many bytes at each event packet?
options.data_bytes = (8 + 2 + 6);
% How many bits dedicated to R-sectors?
options.R_bits = 7;
% How many bits dedicated to modules?
options.M_bits = 2;
% How many bits dedicated to submodules?
options.S_bits = 1;
% How many bits dedicated to crystals?
options.C_bits = 6;
% How many bits dedicated to layers?
options.L_bits = 0;
% What is the coincidence window (in seconds)?
options.coincidence_window = 120e-9;
% Obtain source coordinates? (used in forming the "true" image)
options.source = true;
% What is the clock time step? (see the .cch files)
% If e.g. 1 ps then use 1e12, if 1 ns use 1e9, etc.
options.clock_time_step = 1e12;





%%%%%%%%%%%%%%%%%%%% Root data format settings %%%%%%%%%%%%%%%%%%%%
% Is root data used? (Only one data type can be used at a time)
options.use_root = false;
% Obtain source coordinates? (used in forming the "true" image)
options.source = false;





%%%%%%%%%%%%%%%%%%%% Image properties %%%%%%%%%%%%%%%%%%%%
% reconstructed image pixel size (X-direction)
options.Nx = 128;
% Y-direction
options.Ny = 128;
% Z-direction (number of slices)
options.Nz = options.rings*2-1;





%%%%%%%%%%%%%%%%%%%% Sinogram properties %%%%%%%%%%%%%%%%%%%%
% span factor/axial compression
options.span = 3;
% maximum ring difference
options.ring_difference = 29;
% number of angles in sinogram (this is the final amount after possible
% mashing), maximum allowed is the number of detectors per ring
options.Nang = options.det_per_ring/2;
% number of distances in sinogram
options.Ndist = 200;
% for sinograms, specify the amount of sinograms contained on each segment
% (this should total the total number of sinograms)
options.segment_table = [options.Nz, options.Nz - (options.span + 1):-options.span*2:options.span];
options.segment_table = [options.segment_table(1), repelem(options.segment_table(2:end),2)];
% is the data arc corrected? (NOTE: arc correction is not applied to the
% data even if it is not arc corrected)
% NOTE: At the moment this option does nothing
options.arc_corrected = false;
% Total number of sinograms
options.TotSinos = sum(options.segment_table);
% number of sinograms used in reconstruction
options.NSinos = options.TotSinos;
% If Ndist value is even, take one extra out of the negative side (+1) or
% from the positive side (-1). E.g. if Ndist = 200, then with +1 the
% interval is [-99,100] and with -1 [-100,99].
options.ndist_side = 1;
% How much is the sinogram "rotated"?
options.offangle = options.det_w_pseudo * (3/4);





%%%%%%%%%%%%%%%%%%%% Dynamic imaging properties %%%%%%%%%%%%%%%%%%%%
options.tot_time = 1800; % Total time of the measurement (s)
options.partitions = 1;% how many time points (if a static measurement, use 1)





%%%%%%%%%%%%%%%%%%%% Misc properties %%%%%%%%%%%%%%%%%%%%
% name of current datafile/examination
options.name = 'cylpet_example';

% precompute data (this should be done when using data from a certain
% machine the first time)
% This includes e.g. detector coordinates, sinogram coordinates, etc.
options.precompute = true;

% path: folder for the data (.dat, ASCII) files (must include / at end [or
% \ on Windows])
if isunix % Unix
    options.fpath = '/path/to/GATE/output/';
elseif ispc % Windows
    options.fpath = 'C:\path\to\GATE\output\';
end

% Form only sinograms (no reconstructions)
options.only_sinos = false;

% Precompute the observation matrix for the reconstruction (this might require a
% lot of memory), if false then the observation matrix is calculated on the
% fly (slower)
% NOTE: Supports only MLEM reconstruction
options.precompute_obs_matrix = false;

% Compute only the reconstructions (this option overwrites the precompute
% option)
options.only_reconstructions = false;

% Use raw list mode data
% This means that the data is used as is without any sinogramming and thus
% without any "compression"
options.use_raw_data = false;

% Use precomputed geometrical matrix information
% During the precompute-phase the number of pixels each LOR traverse is
% counted (this phase requires the above precompute-option to true). These
% are saved and later on used during the reconstruction process (compulsory
% for OpenCL reconstruction). Once the precompute-phase has been done once
% for the specific sinogram and image resolution, it is not necessary to do
% it again
options.precompute_lor = true;

% Set this option to true if you want to precompute all possible
% combinations in one go (i.e. raw data, precomputed LORs, sinogram format)
% Requires for precompute option to be true to take effect
options.precompute_all = false;

% Show status messages
options.verbose = true;





%%%%%%%%%%%%%%%%%%%% Reconstruction properties %%%%%%%%%%%%%%%%%%%%
% Reconstruction method used
% 1 = Reconstructions in MATLAB
% 2 = Everything done in OpenCL (CPU or GPU)
% 3 = Reconstructions done sequentially (matrix-free method)
% 4 = Matrix-free reconstruction with OpenCL (parallel)
options.reconstruction_method = 3;
% Device used (this is applicable to methods 2 and 4)
% In methods 2 and 4 this determines the device used for both system matrix
% formation and image reconstruction
% NOTE: On OpenCL the first devices are your GPUs (if you have only 1 GPU,
% then device number 0 is your GPU and 1 your CPU)
options.use_device = 0;
% Reconstruction algorithms to use (you can choose several):
% NOTE: MLEM requires precomputed observation matrix or a matrix-free
% method
% Maximum-Likelihood Expectation Maximization (MLEM)
% Supported by methods 1, 2 and 4
options.mlem = true;
% Ordered Subsets Expectation Maximization (OSEM)
% Supported by all methods
options.osem = false;
% Modified Row-Action Maximum Likelihood Algorithm (MRAMLA, modified BSREM)
% Supported by method 1 only
options.mramla = false;
% Row-Action Maximum Likelihood Algorithm (RAMLA)
% Supported by method 1 only
options.ramla = false;
% Enhanced COSEM (ECOSEM)
% Supported by method 1 only
options.ecosem = false;
% Complete data OSEM (COSEM)
% Supported by method 1 only
options.cosem = false;
% Accelerated COSEM (ACOSEM)
% Supported by method 1 only
options.acosem = false;
% Median root prior (MRP) with One Step Late (OSL) algorithm
% Supported by method 1 only
options.mrp_osl = false;
% Median root prior with Block Sequential Regularized Expectation
% Maximization (BSREM)
% Supported by method 1 only
options.mrp_bsrem = false;
% Quadratic prior with OSL
% Supported by method 1 only
options.quad_osl = false;
% Quadratic prior with BSREM
% Supported by method 1 only
options.quad_bsrem = false;
% L-filter with OSL
% Supported by method 1 only
options.L_osl = false;
% L-filter with BSREM
% Supported by method 1 only
options.L_bsrem = false;
% Finite impulse response (FIR) Median Hybrid (FMH) with OSL
% Supported by method 1 only
options.FMH_osl = false;
% Weighted mean with OSL
% Supported by method 1 only
options.weighted_mean_osl = false;
% number of iterations
options.Niter = 12;
% number of subsets (all excluding MLEM)
options.subsets = 1;
% Initial value for the reconstruction
options.x0 = ones(options.Nx, options.Ny, options.Nz);
% epsilon value (small value to prevent division by zero)
options.epps = 1e-8;
% Include attenuation correction (for this you need attenuation images of
% each slice correctly rotated and scaled for 511 keV)
options.attenuation_correction = false;
% Attenuation image data file (specify the path (if not in MATLAB path) and filename)
% NOTE: the attenuation data has to be the only variable in the file and have the
% dimensions of the final reconstructed image
options.attenuation_datafile = 'GATE_xcat_vaimennus_129.mat';
% Use Shuffle? (recommended)
% NOTE: Applies only when using raw data
% Download from: https://se.mathworks.com/matlabcentral/fileexchange/27076-shuffle
options.use_Shuffle = true;
% Use fast sparse? (recommended)
% Download from: https://github.com/stefanengblom/stenglib
% NOTE: This applies only to method 1 when precompute_lor is false
options.use_fsparse = true;
% accleration parameter for ACOSEM
options.h = 1;
% relaxation parameter for RAMLA
options.b0 = 1;
% Upper bound for MRAMLA
options.U = 0;
% regularization parameter for MRP with OSL
options.beta_mrp_osl = 0.05;
% regularization parameter for MRP with BSREM
options.beta_mrp_bsrem = 0.05;
% distance of where the median is taken on MRP
options.medx = 3;
options.medy = 3;
options.medz = 3;
% regularization parameter for quadratic prior with OSL
options.beta_quad_osl = 0.05;
% regularization parameter for quadratic prior with BSREM
options.beta_quad_bsrem = 0.05;
% how many neighboring pixels are taken into account (with quadratic, MRP, L, FMH and weighted)
% NOTE: At the moment Ndx and Ndy has to be identical
options.Ndx = 1;
options.Ndy = 1;
options.Ndz = 1;
% pixel weights for quadratic prior
% the number of pixels need to be the amount of neighboring pixels
% e.g. if the above Nd values are all 1, then 27 weights need to be included
% where the center pixel (if Nd values are 1, element 14) should be Inf
% Size is (Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)
% If left empty then they will be calculated by the algorithm
options.weights = [];
% regularization parameter for L-filter with OSL
options.beta_L_osl = 0.05;
% regularization parameter for L-filter with BSREM
options.beta_L_bsrem = 0.05;
% weighting factors for the L-filter pixels
% Otherwise the same as in quadratic prior, but center pixel is not Inf
options.a_L = [];
% pixel weights for FMH
% The size needs to be [3, floor((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)/2)]
options.fmh_weights = [];
% additional weighting factor for the center row of fmh_weights
% if 1 then the maximum value of all weights are used as the center row weight
% In other words it multiplies the maximum value of all other weights
options.fmh_center_weight = 1;
% regularization parameter for FMH with OSL
options.beta_fmh_osl = 0.05;
% regularization parameter for weighted mean with OSL
options.beta_weighted_osl = 0.05;
% pixel weights for weighted mean
% the number of pixels need to be the amount of neighboring pixels
% e.g. if the above Nd values are all 1, then 27 weights need to be included
% where the center pixel (if Nd values are 1, element 14) should be Inf
% Size is (Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)
% If left empty then they will be calculated by the algorithm
options.weighted_weights = [];
% center pixel weight for weighted mean
% if 1 then the maximum value of all weights are used as the center weight
% In other words it multiplies the maximum value of all other weights
options.weighted_center_weight = 1;

%% Precompute the necessary data

if options.precompute && options.only_reconstructions == false
    precompute_data(options);
end

%% Load the ASCII/LMF coindicence data

if options.only_reconstructions == false && options.only_sinos == false
    options.coincidences = load_GATE_data(options);
end

%% Form the sinograms

if options.only_reconstructions == false && options.use_raw_data == false
    options.SinM = form_sinograms(options);
end

%% Reconstructions

if options.only_sinos == false
    
    tStart = tic;
    pz = reconstructions_main(options);
    tElapsed = toc(tStart);
    disp(['Reconstruction process took ' num2str(tElapsed) ' seconds'])
    
end

save([options.name '_reconstruction_' num2str(options.subsets) 'subsets_' num2str(options.Niter) 'iterations.mat'], 'pz');