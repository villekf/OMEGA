%% MATLAB codes for GATE PET reconstruction using ASCII or LMF output

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Machine properties %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R-sectors/blocks in transaxial direction
options.blocks_per_ring = (42);
% R-sectors/modules/blocks in axial direction (i.e. number of physical
% machine rings)
options.linear_multip = (4);
% number of detectors on the side of R-sector/block/module (e.g. 13 if 13x13)
options.cryst_per_block = (8);
% crystal pitch in x- and y-directions (mm)
options.cr_p = 2.4;
% crystal pitch in z-direction (mm)
options.cr_pz = 2.4;
% ring diameter (distance between perpendicular detectors) in mm
options.diameter = 130*2;
% Transaxial FOV size (mm), this is the length of the x (horizontal) side
% of the FOV
options.FOVa_x = 151;
% Transaxial FOV size (mm), this is the length of the y (vertical) side
% of the FOV
options.FOVa_y = options.FOVa_x;
% Axial FOV (mm)
options.axial_fov = floor(76.8 - options.cr_pz/2);
% Number of pseudo rings between physical rings (use 0 or [] if none)
options.pseudot = [];
% Number of detectors per crystal ring (without pseudo detectors)
options.det_per_ring = options.blocks_per_ring*options.cryst_per_block;
% Number of detectors per crystal ring (with pseudo detectors)
% If your scanner has a single pseudo detector on each side of the crystal
% block then simply add +1 inside the parenthesis
options.det_w_pseudo = options.blocks_per_ring*(options.cryst_per_block);
% number of crystal rings
options.rings = options.linear_multip * options.cryst_per_block;
% number of detectors
options.detectors = options.det_per_ring*options.rings;
% machine name
options.machine_name = 'Cylindrical_PET_example';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% GATE specific settings %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain trues (true coincidences)
% If this is set to true then, in addition to the normal coincidences,
% trues are also obtained and saved
options.obtain_trues = true;
% Reconstruct the true coincidences
% If this is set to true, then the true coincidences will be used for
% reconstruction
% NOTE: If both this and reconstruct_trues are set, then the trues are
% reconstructed, but not the scatter
options.reconstruct_trues = false;
% Obtain scattered coincidences
% If this is set to true, then scattered coincidences are saved separately
% These events are not used for scatter correction though, but a separate
% scatter sinogram will be created if sinograms are created with this
% option set to true
options.store_scatter = true;
% What scatter components are included in the scatter part(1 means that
% component is included, 0 means it is not included in the scatter data)
% First: Compton scattering in the phantom, second: Compton scattering in
% the detector, third: Rayleigh scattering in the phantom, fourth: Rayleigh
% scattering in the detector
% If store_scatter is set to true, at least one value has to be 1
% NOTE: LMF will always include only Compton scattering in the phantom,
% regardless of the choice below (as long as scatter is selected)
options.scatter_components = [1 0 1 0];
% Reconstruct the scattered coincidences
% If this is set to true, then the scattered coincidences will be used for
% reconstruction
% NOTE: If both this and reconstruct_trues are set, then the trues are
% reconstructed, but not the scatter
options.reconstruct_scatter = false;
% Obtain (true) random coincidences.
% If this is set to true then coincidence events that are actually random
% events are stored separately.
% These events are not used for randoms correction (see the
% Corrections-section for delayed coincidence window randoms correction),
% but a separate randoms sinogram will be created if sinograms are created
% with this option set to true.
options.store_randoms = true;
% Obtain source coordinates? (used in forming the "true" image).
% If this is set to true, then the "true" decay image is also saved during
% data load.
% If any of the above settings are set to true, then the true images are
% also obtained for them. E.g. if store_scatter = true, then an image
% showing the locations and number of counts of where the scattered events
% originated will be saved in a mat-file. Scatter and trues contain
% coincidence events while randoms contain singles.
% NOTE: If you use LMF data, the source images are not considered reliable
options.source = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% ASCII data format settings %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Is ASCII data used? (Only one data type can be used at a time)
options.use_ASCII = true;
% Columns corresponding to the R-sectors (or blocks) (see GATE wiki on ASCII file
% format)
% Corresponds to (from GATE Wiki): Column 12 to 17 : volume IDs* (first single)
options.rsector_ind1 = (7);
options.rsector_ind2 = (22);
% Module indices (columns)
% if no modules present use 0 (i.e. ECAT geometry)
options.module_ind1 = (8);
options.module_ind2 = (23);
% Crystal indices (columns)
options.crs_ind1 = (10);
options.crs_ind2 = (25);
% Time index (column), use either the first single or second
% Corresponds to (from GATE Wiki): Column 7 : Time stamp (first single)
options.time_index = (5);
% index of the source location of the first single (X-dimension, first column of the three)
% Corresponds to (from GATE Wiki): Column 4 to 6 : XYZ position of the source in world referential (first single)
% Use 0 or [] to skip this step
options.source_index1 = 2;
% index of the source location of the second single (X-dimension, first column of the three)
% Corresponds to (from GATE Wiki): Column 27 to 29 : XYZ position of the source in world referential (second single)
% Use 0 or [] to skip this step
options.source_index2 = 17;
% event ID index of first single
% Corresponds to (from GATE Wiki): Column 2 : ID of the event (first single)
% This is only applicable if any of the GATE specific settings are set to
% true
options.event_index1 = 1;
% event ID index of second single
% Corresponds to (from GATE Wiki): Column 25 : ID of the event (second single)
% This is only applicable if any of the GATE specific settings are set to
% true
options.event_index2 = 16;
% Index for number of Compton interaction in phantom/detector, first single
% Corresponds to (from GATE Wiki): Column 18 : Number of Compton interactions in phantoms before reaching the detector (first single)
% This is only applicable if store_scatter = true
% If Compton interactions in phantom are NOT selected, but detector
% interactions are, then this corresponds to the column for Compton scatter
% in detector
% If Compton interactions are NOT selected, but Rayleigh are then this
% corresponds to the column index of Rayleigh scatter interactions
options.scatter_index1 = 12;
% Index for number of Compton interaction in phantom, second single
% Corresponds to (from GATE Wiki): Column 41 : Number of Compton interactions in phantoms before reaching the detector (second single)
% This is only applicable if store_scatter = true
% Same rules apply as above
options.scatter_index2 = 27;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% LFM data format settings %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Is LMF data used? (Only one data type can be used at a time)
options.use_LMF = false;
% How many bytes are in the LMF header part?
options.header_bytes = (16);
% How many bytes in each event packet?
options.data_bytes = (8 + 2 + 6 + 4 + 1);
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
options.coincidence_window = 10e-9;
% Obtain source coordinates? (used in forming the "true" image)
options.source = true;
% What is the clock time step? (see the .cch files)
% If e.g. 1 ps then use 1e-12, if 1 ns use 1e-9, etc.
options.clock_time_step = 1e-12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Root data format settings %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Is root data used? (Only one data type can be used at a time)
options.use_root = false;
% Obtain source coordinates? (used in forming the "true" image)
options.root_source = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Image properties %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reconstructed image pixel size (X-direction)
options.Nx = 128;
% Y-direction
options.Ny = 128;
% Z-direction (number of slices)
options.Nz = options.rings*2-1;
% Flip the image (in vertical direction)?
options.flip_image = false;
% How much is the image rotated?
% You need to run the precompute phase again if you modify this
% NOTE: The rotation is done in the detector space (before reconstruction)
options.offangle = options.det_w_pseudo * (3/4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Sinogram properties %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% span factor/axial compression
options.span = 3;
% maximum ring difference
options.ring_difference = 29;
% number of angles in sinogram (this is the final amount after possible
% mashing), maximum allowed is the number of detectors per ring/2
options.Nang = options.det_per_ring/2;
% number of angular positions in sinogram
options.Ndist = 200;
% for sinograms, specify the amount of sinograms contained on each segment
% (this should total the total number of sinograms)
options.segment_table = [options.Nz, options.Nz - (options.span + 1):-options.span*2:options.span];
options.segment_table = [options.segment_table(1), repelem(options.segment_table(2:end),2)];
% Total number of sinograms
options.TotSinos = sum(options.segment_table);
% number of sinograms used in reconstruction
options.NSinos = options.TotSinos;
% If Ndist value is even, take one extra out of the negative side (+1) or
% from the positive side (-1). E.g. if Ndist = 200, then with +1 the
% interval is [-99,100] and with -1 [-100,99].
options.ndist_side = 1;
% Store also the raw sinogram without any corrections applied (gap filling,
% randoms correction, etc.)
% This is applicable only if any of the corrections are selected
options.store_raw_sinogram = false;
% Fill the gaps caused by pseudo detectors?
% NOTE: Applicable only if options.pseudot > 0
options.fill_sinogram_gaps = true;
% Which method used to fill the gaps?
% Either MATLAB's built-in fillmissing or inpaint_nans from file exchange
% For inpaint_nans see: https://se.mathworks.com/matlabcentral/fileexchange/4551-inpaint_nans
% See wiki for more details
options.gap_filling_method = 'fillmissing';
% interpolation method used with fillmissing
% Possible methods are those listed under method-section in fillmissing
options.interpolation_method_fillmissing = 'linear';
% interpolation method used with inpaint_nans
% See inpaint_nans.m for details
options.interpolation_method_inpaint = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Corrections %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Randoms correction
% If set to true, stores the delayed coincidences during data load and
% later corrects for randoms during sinogram formation (if sinogram data)
% or during data load (if raw list-mode data)
options.randoms_correction = false;
% Scatter correction
% If set to true, will prompt the user to load the scatter sinogram/raw
% data. Corrects for scatter during sinogram formation (if sinogram data)
% or during data load (if raw list-mode data)
% NOTE: Scatter data is not created by this software and as such has to be
% provided by the user
options.scatter_correction = false;
% Attenuation correction
% Blank scan
% If set to true, the blank data will be used to form the blank
% sinogram/raw list-mode data file
% This will prompt the folder where the blank data is located
% This step should be done only once for each blank file
options.blank = false;
% Transmission scan
% If set to true, form the transmission sinogram/raw list-mode data file
% This will prompt the folder where the transmission data is located
options.transmission_scan = false;
% Image-based attenuation correction
% Include attenuation correction from images (e.g. CT-images) (for this you
% need attenuation images of each slice correctly rotated and scaled for
% 511 keV) 
options.attenuation_correction = false;
% Attenuation image data file (specify the path (if not in MATLAB path) and filename)
% NOTE: the attenuation data has to be the only variable in the file and have the
% dimensions of the final reconstructed image
options.attenuation_datafile = '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Dynamic imaging properties %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total time of the measurement (s)
options.tot_time = 1800;
% how many time points/dynamic frames (if a static measurement, use 1)
options.partitions = 1;
% Start time (s) (all measurements before this will be ignored)
options.start = 0;
% End time (s) (all measurements after this will be ignored)
options.end = options.tot_time;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Misc properties %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% name of current datafile/examination
options.name = 'cylpet_example';
% precompute data (this should be done when using data from a certain
% machine the first time)
% This includes e.g. detector coordinates, sinogram coordinates, etc.
options.precompute = true;
% folder for the data (.dat ASCII, .ccs LMF, .root ROOT) files
if isunix % Unix
    options.fpath = '/path/to/GATE/output/';
elseif ispc % Windows
%     options.fpath = 'C:\path\to\GATE\output\';
    options.fpath = 'C:\Users\villewe\OneDrive - University of Eastern Finland\MAT-tiedostot\PET\GATE\example\Uus';
%     options.fpath = 'I:\Väikkäri\MAT-tiedostot\PET\example\Uus\Arc';
end
% Form only sinograms (no reconstructions)
options.only_sinos = false;
% Precompute the observation matrix for the reconstruction (this might require a
% lot of memory), if false then the observation matrix is calculated on the
% fly
% NOTE: Supports only MLEM reconstruction
options.precompute_obs_matrix = false;
% Compute only the reconstructions (this option overwrites the precompute
% option, but not precompute_all)
options.only_reconstructions = true;
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
% Precompute variables needed by MBSREM, RBI and COSEM (and its variants)
% If set to true, some variables needed by the above algorithms are
% computed during the precompute phase
% It is recommended to have this true if any of the above methods are ever
% used, but they can also be computed on-the-fly
options.precompute_MBSREM = true;
% Set this option to true if you want to precompute all possible
% combinations in one go (i.e. raw data, precomputed LORs, sinogram format)
% Requires for precompute option to be true to take effect
% Setting this option to true also causes the precompute phase to be
% performed even if only_reconstructions is true
% Additionally, setting this to true also causes the data load phase to
% save both the data formatted for sinogram formation as well as the raw
% list-mode data
options.precompute_all = true;
% Show status messages
options.verbose = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Reconstruction method %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction method used
% 1 = Reconstructions in MATLAB
% 2 = Matrix-free reconstruction with OpenCL (parallel) Recommended method
% 3 = Reconstructions done sequentially (matrix-free method)
% 4 = Everything done in OpenCL (CPU or GPU) using matrices
options.reconstruction_method = 2;
% Device used (this is applicable to methods 2 and 4)
% In methods 2 and 4 this determines the device used for both system matrix
% formation and image reconstruction
% NOTE: On OpenCL the first devices are your GPUs (if you have only 1 GPU,
% then device number 0 is your GPU and 1 your CPU), you can use
% Arrayfire_OpenCL_device_info() to determine the device numbers
options.use_device = 1;
% NOTE: if you switch devices then you need to run the below line
% (uncommented) as well:
% clear mex
% Force the (re)building of OpenCL binaries
% If set to true, the OpenCL binaries are rebuilt even if they have been
% previously built
% Use this once if you update your drivers or there are changes made to the
% .cl-file
options.force_build = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Reconstruction algorithms %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction algorithms to use (you can choose several)
% NOTE: MLEM requires precomputed observation matrix or a matrix-free
% method

%%%%%%%%%%%% ML-Methods %%%%%%%%%%%%
% These are non-regularized versions
% Maximum-Likelihood Expectation Maximization (MLEM)
% Supported by methods 1, 2 and 3
options.mlem = false;
% Ordered Subsets Expectation Maximization (OSEM)
% Supported by all methods
options.osem = true;
% Modified Row-Action Maximum Likelihood Algorithm (MRAMLA)
% Supported by methods 1 and 2
options.mramla = true;
% Row-Action Maximum Likelihood Algorithm (RAMLA)
% Supported by methods 1 and 2
options.ramla = false;
% Relaxed Ordered Subsets Expectation Maximization (ROSEM)
% Supported by methods 1 and 2
options.rosem = true;
% Rescaled Block Iterative Expectation Maximization (RBI-EM)
% Supported by methods 1 and 2
options.rbi = false;
% Dynamic RAMLA (DRAMA)
% Supported by methods 1 and 2
options.drama = false;
% Complete data OSEM (COSEM)
% Supported by methods 1 and 2
options.cosem = false;
% Enhanced COSEM (ECOSEM)
% Supported by methods 1 and 2
options.ecosem = false;
% Accelerated COSEM (ACOSEM)
% Supported by methods 1 and 2
options.acosem = false;


%%%%%%%%%%%% MAP-Methods %%%%%%%%%%%%
% Any algorithm selected here will utilize all the priors selected below.
% For example, if OSL-OSEM is set to true and MRP and Quad are set to true,
% then OSL-OSEM estimates will be computed for both MRP and Quadratic
% prior.
% One-Step Late MLEM (OSL-MLEM)
% Supported by method 2 only
options.OSL_MLEM = false;
% One-Step Late OSEM (OSL-OSEM)
% Supported by methods 1 and 2
options.OSL_OSEM = true;
% Modified BSREM (MBSREM)
% Supported by methods 1 and 2
options.MBSREM = false;
% Block Sequantial Regularized Expectaction Maximation (BSREM)
% Supported by methods 1 and 2
options.BSREM = false;
% ROSEM-MAP
% Supported by methods 1 and 2
options.ROSEM_MAP = false;
% RBI-MAP
% Supported by methods 1 and 2
options.RBI_MAP = false;
% (A)COSEM-MAP (OSL)
% 0/false = No COSEM-MAP, 1/true = ACOSEM-MAP, 2 = COSEM-MAP
% Supported by methods 1 and 2
options.COSEM_MAP = false;


%%%%%%%%%%%% Priors %%%%%%%%%%%%
% Median Root Prior (MRP)
options.MRP = true;
% Quadratic prior (QP)
options.quad = false;
% L-filter prior
options.L = false;
% Finite impulse response (FIR) Median Hybrid (FMH) prior
options.FMH = false;
% Weighted mean prior
options.weighted_mean = false;
% Total Variation (TV) prior
options.TV = false;
% Anisotropic Diffusion (AD) prior
options.AD = false;
% Asymmetric Parallel Level Set (APLS) prior
options.APLS = false;
% Total Generalized Variation (TGV) prior
options.TGV = false;
% Non-local Means (NLM) prior
options.NLM = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Reconstruction properties %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of iterations (all reconstruction methods)
options.Niter = 1;
% Number of subsets (all excluding MLEM and subset_type = 5)
options.subsets = 8;
% Subset type (n = subsets)
% 1 = Every nth measurement is taken (e.g. if subsets = 3, then first 
% subset has measurements 1, 4, 7, etc., second 2, 5, 8, etc.) 
% 2 = Measurements are selected randomly
% 3 = (Sinogram only) Take every nth column in the sinogram
% 4 = (Sinogram only) Take every nth row in the sinogram
% 5 = Sort the LORs according to their angle with positive X-axis, combine
% n_angles together and have 180/n_angles subsets for 2D slices and
% 360/n_angles for 3D, see GitHub wiki for more information
options.subset_type = 1;
% How many angles are combined in subset_type = 5
% E.g. there are 180 angles, in n_angles = 2, then angles 0 and 1 are
% combined, 2 and 3, etc.
options.n_angles = 2;
% Initial value for the reconstruction
options.x0 = ones(options.Nx, options.Ny, options.Nz);
% options.x0 = load('initial_value.mat','alkuarvo');
% options.x0 = options.x0.alkuarvo;
% epsilon value (small value to prevent division by zero)
options.epps = 1e-8;
% Include attenuation correction (for this you need attenuation images of
% each slice correctly rotated and scaled for 511 keV)
options.attenuation_correction = false;
% Attenuation image data file (specify the path (if not in MATLAB path) and filename)
% NOTE: the attenuation data has to be the only variable in the file and have the
% dimensions of the final reconstructed image
options.attenuation_datafile = '';
% Use Shuffle? (recommended)
% NOTE: Applies only when using subset_type = 2
% Download from: https://se.mathworks.com/matlabcentral/fileexchange/27076-shuffle
options.use_Shuffle = true;
% Use fast sparse? (recommended)
% Download from: https://github.com/stefanengblom/stenglib
% NOTE: This applies only to method 1 when precompute_lor is false
options.use_fsparse = true;
% Skip the normalization phase in MRP, FMH, L-filter, AD and/or weighted
% mean
% E.g. if set to true the MRP prior is (x - median(x))
% E.g. if set to false the MRP prior is (x - median(x)) / median(x)
options.med_no_norm = false;


%%%%%%%%%%%%%%%%%%%% ACOSEM properties %%%%%%%%%%%%%%%%%%%%
% Accleration parameter for ACOSEM (1 equals COSEM)
options.h = 2;

%%%%%%%%%%%%%%%%%%%% MRAMLA and MBSREM properties %%%%%%%%%%%%%%%%%%%%
% Relaxation parameter for MRAMLA and MBSREM
% Use scalar if you want it to decrease as lambda0_mbsrem/current_iteration
% Use vector (length = Niter) if you want your own relaxation parameters
options.lambda0_mbsrem = 0.2;
% Upper bound for MRAMLA/MBSREM (use 0 for default value)
options.U = 0;

%%%%%%%%%%%%%%%%%%%% RAMLA and BSREM properties %%%%%%%%%%%%%%%%%%%%
% Relaxation parameter for RAMLA and BSREM
% Use scalar if you want it to decrease as lambda0/current_iteration
% Use vector (length = Niter) if you want your own relaxation parameters
options.lambda0 = 0.1;

%%%%%%%%%%%%%%%%%%%% ROSEM and ROSEM-MAP properties %%%%%%%%%%%%%%%%%%%%
% Relaxation parameter for ROSEM and ROSEM-MAP
% Use scalar if you want it to decrease as lambda0_rosem/current_iteration
% Use vector (length = Niter) if you want your own relaxation parameters
options.lambda0_rosem = 0.5;


%%%%%%%%%%%%%%%%%%%% DRAMA properties %%%%%%%%%%%%%%%%%%%%
% Beta_0 value
options.beta0_drama = 0.1;
% Beta value
options.beta_drama = 1;
% Alpha value
options.alpha_drama = 0.1;


%%%%%%%%%%%%%%%%%%%% Weighting properties %%%%%%%%%%%%%%%%%%%%
% how many neighboring pixels are taken into account (with MRP, quadratic,
% L, FMH, NLM and weighted mean)
% E.g. if Ndx = 1, Ndy = 1, Ndz = 0, then you have 3x3 square area where
% the pixels are taken into account (I.e. (Ndx*2+1)x(Ndy*2+1)x(Ndz*2+1)
% area)
% NOTE: At the moment Ndx and Ndy has to be identical
options.Ndx = 1;
options.Ndy = 1;
options.Ndz = 1;


%%%%%%%%%%%%%%%%%%%% MRP properties %%%%%%%%%%%%%%%%%%%%
% Regularization parameter for MRP with OSL-OSEM
options.beta_mrp_osem = 0.1;%0.3;
% Regularization parameter for MRP with OSL-MLEM
options.beta_mrp_mlem = 1.5;
% Regularization parameter for MRP with MBSREM
options.beta_mrp_mbsrem = 0.3;
% Regularization parameter for MRP with BSREM
options.beta_mrp_bsrem = 0.1;
% Regularization parameter for MRP with ROSEM
options.beta_mrp_rosem = 2;
% Regularization parameter for MRP with RBI
options.beta_mrp_rbi = 0.1;
% Regularization parameter for MRP with OSL-(A)COSEM
options.beta_mrp_cosem = 1;


%%%%%%%%%%%%%%%%%%%% Quadratic prior properties %%%%%%%%%%%%%%%%%%%%
% Regularization parameter for quadratic prior with OSL-OSEM
options.beta_quad_osem = 0.01;%0.1;
% Regularization parameter for quadratic prior with OSL-MLEM
options.beta_quad_mlem = 0.1;
% Regularization parameter for quadratic prior with MBSREM
options.beta_quad_mbsrem = 0.05;
% Regularization parameter for quadratic prior with BSREM
options.beta_quad_bsrem = 0.03;
% Regularization parameter for quadratic prior with ROSEM
options.beta_quad_rosem = 0.1;
% Regularization parameter for quadratic prior with RBI
options.beta_quad_rbi = 0.05;
% Regularization parameter for quadratic prior (OSL-(A)COSEM)
options.beta_quad_cosem = 0.01;
% pixel weights for quadratic prior
% the number of pixels need to be the amount of neighboring pixels
% e.g. if the above Nd values are all 1, then 27 weights need to be included
% where the center pixel (if Nd values are 1, element 14) should be Inf
% Size is (Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)
% If left empty then they will be calculated by the algorithm
options.weights = [];


%%%%%%%%%%%%%%%%%%%% L-filter properties %%%%%%%%%%%%%%%%%%%%
% Regularization parameter for L-filter with OSL-OSEM
options.beta_L_osem = 0.1;
% Regularization parameter for L-filter with OSL-MLEM
options.beta_L_mlem = 0.1;
% Regularization parameter for L-filter with MBSREM
options.beta_L_mbsrem = 0.1;
% Regularization parameter for L-filter with BSREM
options.beta_L_bsrem = 0.03;
% Regularization parameter for L-filter with ROSEM
options.beta_L_rosem = 3;
% Regularization parameter for L-filter with RBI
options.beta_L_rbi = 0.09;
% Regularization parameter for L-filter (OSL-(A)COSEM)
options.beta_L_cosem = 0.1;
% weighting factors for the L-filter pixels
% Otherwise the same as in quadratic prior, but center pixel is not Inf
% If left empty then they will be calculated by the algorithm such that the
% weights resemble a Laplace distribution
options.a_L = [];
% If the weighting factors are set empty, then this option will determine
% whether the computed weights follow a 1D weighting scheme (true) or 2D 
% (false)
% See the wiki for more information
options.oneD_weights = false;


%%%%%%%%%%%%%%%%%%%% FMH properties %%%%%%%%%%%%%%%%%%%%
% Regularization parameter for FMH with OSL-OSEM
options.beta_fmh_osem = 0.1;
% Regularization parameter for FMH with OSL-MLEM
options.beta_fmh_mlem = 0.1;
% Regularization parameter for FMH with MBSREM
options.beta_fmh_mbsrem = 0.6;
% Regularization parameter for FMH with BSREM
options.beta_fmh_bsrem = 5;
% Regularization parameter for FMH with ROSEM
options.beta_fmh_rosem = 8;
% Regularization parameter for FMH with RBI
options.beta_fmh_rbi = 0.5;
% Regularization parameter for FMH (OSL-(A)COSEM)
options.beta_fmh_cosem = 0.1;
% pixel weights for FMH
% The size needs to be [Ndx*2+1, 4] if Nz = 1 or Ndz = 0, or [Ndx*2+1, 13]
% otherwise
% The center pixel weight should be in the middle
% If the sum of each column is > 1, then the weights will be normalized
% such that the sum = 1
% If left empty then they will be calculated by the algorithm such that the
% weights follow the same pattern as in the original article
options.fmh_weights = [];
% Weighting value for the center pixel
% Default value is 4, which was used in the original article
% NOTE: This option is ignored if you provide your own weights
options.fmh_center_weight = 4;


%%%%%%%%%%%%%%%%%%%% Weighted mean properties %%%%%%%%%%%%%%%%%%%%
% Mean type
% 1 = Arithmetic mean, 2 = Harmonic mean, 3 = Geometric mean
options.mean_type = 1;
% Regularization parameter for weighted mean with OSL-OSEM
options.beta_weighted_osem = 0.02;%0.2;
% Regularization parameter for weighted mean with OSL-MLEM
options.beta_weighted_mlem = 0.1;
% Regularization parameter for weighted mean with MBSREM
options.beta_weighted_mbsrem = 0.1;
% Regularization parameter for weighted mean with BSREM
options.beta_weighted_bsrem = 5;
% Regularization parameter for weighted mean with ROSEM
options.beta_weighted_rosem = 3;
% Regularization parameter for weighted mean with RBI
options.beta_weighted_rbi = 0.04;
% Regularization parameter for weighted mean (OSL-(A)COSEM)
options.beta_weighted_cosem = 0.2;
% pixel weights for weighted mean
% the number of pixels need to be the amount of neighboring pixels
% e.g. if the above Nd values are all 1, then 27 weights need to be included
% Size is (Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)
% If left empty then they will be calculated by the algorithm such that the
% weights are dependent on the distance from the center pixel to the
% neighboring pixels
options.weighted_weights = [];
% center pixel weight for weighted mean
% NOTE: This option is ignored if you provide your own weights
options.weighted_center_weight = 4;


%%%%%%%%%%%%%%%%%%%% TV properties %%%%%%%%%%%%%%%%%%%%
% Regularization parameter for TV OSL-OSEM
options.beta_TV_osem = 0.3;
% Regularization parameter for TV with OSL-MLEM
options.beta_TV_mlem = 0.1;
% Regularization parameter for TV with MBSREM
options.beta_TV_mbsrem = 0.1;
% Regularization parameter for TV with BSREM
options.beta_TV_bsrem = 0.05;
% Regularization parameter for TV with ROSEM
options.beta_TV_rosem = 0.7;
% Regularization parameter for TV with RBI
options.beta_TV_rbi = 0.02;
% Regularization parameter for TV (OSL-(A)COSEM)
options.beta_TV_cosem = 0.03;
% "Smoothing" parameter
% Also used to prevent zero values in square root
options.TVsmoothing = 1e-1;
% Whether to use an anatomical reference/weighting image with the TV
options.TV_use_anatomical = false;
% Three different TV methods are available.
% Value can be 1, 2 or 3.
% Types 1 and 2 are the same if anatomical prior is not included
% Type 3 uses the same weights as quadratic prior
% See the wiki for more information
options.TVtype = 3;
% If the TV_use_anatomical value is set to true, specify filename for the
% reference image here (same rules apply as with attenuation correction
% above)
options.TV_reference_image = 'reference_image.mat';
% Weighting parameters for the TV prior. Applicable only if
% use_anatomical = true
% T-value is specific to the used TVtype, e.g. for type 1 it is the edge
% threshold parameter. See the wiki for more details
options.T = 0.1;
% C is the weight for the original image in type 3 and is ignored with
% other types
options.C = 1;
% Tuning parameter for TV and APLS
options.tau = 1e-8;


%%%%%%%%%%%%%%%%%%%% AD properties %%%%%%%%%%%%%%%%%%%%
% Regularization parameter for AD with OSL-OSEM
options.beta_ad_osem = 0.1;
% Regularization parameter for AD with OSL-MLEM
options.beta_ad_mlem = 0.1;
% Regularization parameter for AD with MBSREM
options.beta_ad_mbsrem = 0.3;
% Regularization parameter for AD with BSREM
options.beta_ad_bsrem = 0.2;
% Regularization parameter for AD with ROSEM
options.beta_ad_rosem = 3;
% Regularization parameter for AD with RBI
options.beta_ad_rbi = 0.05;
% Regularization parameter for AD with (OSL-(A)COSEM)
options.beta_ad_cosem = 0.1;
% Time step variable for AD (method 2 only)
options.TimeStepAD = 0.0625;
% Conductivitiy/connectivity for AD (edge threshold)
options.KAD = 1;
% Number of iterations for AD filter
% NOTE: This refers to the AD smoothing part, not the actual reconstruction
% phase
options.NiterAD = 5;
% Flux/conduction type for AD filter
% 1 = Exponential
% 2 = Quadratic
options.FluxType = 1;
% Diffusion type for AD (method 2 only)
% 1 = Gradient
% 2 = Modified curvature
options.DiffusionType = 1;


%%%%%%%%%%%%%%%%%%%% APLS properties %%%%%%%%%%%%%%%%%%%%
% Regularization parameter for APLS with OSL-OSEM
options.beta_APLS_osem = 0.01;
% Regularization parameter for APLS with OSL-MLEM
options.beta_APLS_mlem = 0.1;
% Regularization parameter for APLS with MBSREM
options.beta_APLS_mbsrem = 0.1;
% Regularization parameter for APLS with BSREM
options.beta_APLS_bsrem = 0.005;
% Regularization parameter for APLS with ROSEM
options.beta_APLS_rosem = 0.1;
% Regularization parameter for APLS with RBI
options.beta_APLS_rbi = 0.1;
% Regularization parameter for APLS (OSL-(A)COSEM)
options.beta_APLS_cosem = 0.01;
% "Smoothing" parameter 1 (eta)
% Also used to prevent zero values in square root
options.eta = 1e-5;
% "Smoothing" parameter 2 (beta)
% Also used to prevent zero values in square root
options.APLSsmoothing = 1e-5;
% Specify filename for the reference image here (same rules apply as with
% attenuation correction above)
% NOTE: For APSL, the prior is required
options.APLS_reference_image = 'reference_image.mat';


%%%%%%%%%%%%%%%%%%%% TGV properties %%%%%%%%%%%%%%%%%%%%
% Regularization parameter for TGV with OSL-OSEM
options.beta_TGV_osem = 0.05;
% Regularization parameter for TGV with OSL-MLEM
options.beta_TGV_mlem = 0.1;
% Regularization parameter for TGV with MBSREM
options.beta_TGV_mbsrem = 0.1;
% Regularization parameter for TGV with BSREM
options.beta_TGV_bsrem = 1;
% Regularization parameter for TGV with ROSEM
options.beta_TGV_rosem = 0.25;
% Regularization parameter for TGV with RBI
options.beta_TGV_rbi = 0.1;
% Regularization parameter for TGV (OSL-(A)COSEM)
options.beta_TGV_cosem = 0.05;
% "Smoothing" parameter 1 (eta)
% Also used to prevent zero values in square root
options.alphaTGV = 2;
options.betaTGV = 1;
options.NiterTGV = 30;


%%%%%%%%%%%%%%%%%%%% NLM properties %%%%%%%%%%%%%%%%%%%%
% EXPERIMENTAL FEATURE
% Regularization parameter for NLM with OSL-OSEM
options.beta_NLM_osem = 0.025;
% Regularization parameter for NLM with MBSREM
options.beta_NLM_mbsrem = 0.05;
% Regularization parameter for NLM with BSREM
options.beta_NLM_bsrem = 0.01;
% Regularization parameter for NLM with ROSEM
options.beta_NLM_rosem = 0.1;
% Regularization parameter for NLM with RBI
options.beta_NLM_rbi = 0.01;
% Regularization parameter for NLM (OSL-(A)COSEM)
options.beta_NLM_cosem = 0.01;
% Filter parameter
options.sigma = 0.01;
% Patch radius
options.Nlx = 1;
options.Nly = 1;
options.Nlz = 0;
% Search window radius is controlled by Ndx, Ndy and Ndz parameters
% Use anatomical reference image for the patches
options.NLM_use_anatomical = true;
% Specify filename for the reference image here (same rules apply as with
% attenuation correction above)
options.NLM_reference_image = 'reference_image.mat';
% Use MRP algorithm (without normalization)
% I.e. gradient = im - NLM_filtered(im)
options.NLM_MRP = true;


%%%%%%%%%%%%%%%%%%%% Custom prior properties %%%%%%%%%%%%%%%%%%%%
% Regularization parameter for custom prior with OSL-OSEM
options.beta_custom_osem = 0.1;%0.3;
% Regularization parameter for custom prior with OSL-MLEM
options.beta_custom_mlem = 1.5;
% Regularization parameter for custom prior with MBSREM
options.beta_custom_mbsrem = 0.3;
% Regularization parameter for custom prior with BSREM
options.beta_custom_bsrem = 0.1;
% Regularization parameter for custom prior with ROSEM
options.beta_custom_rosem = 2;
% Regularization parameter for custom prior with RBI
options.beta_custom_rbi = 0.1;
% Regularization parameter for custom prior with OSL-(A)COSEM
options.beta_custom_cosem = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% OpenCL device info %%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uncomment the below line and run it to determine the available devices
% and their respective numbers
% OpenCL_device_info();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Load scatter data %%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load scatter data (if applicable)
if options.scatter_correction && ~options.only_reconstructions || options.scatter_correction && options.use_raw_data
    [scatter_file, s_fpath] = uigetfile('*.mat','Select scatter correction data');
    
    FileName = fullfile(s_fpath, scatter_file);
    storedStructure = load(FileName);
    variables = fields(storedStructure);
    
    options.ScatterC = storedStructure.(variables{1});
    clear scatter_file s_fpath FileName storedStructure variables
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Error checking %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic error checking is done here
options = OMEGA_error_check(options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Precompute the necessary data

if options.precompute && options.only_reconstructions == false && options.precompute_all == false || options.precompute_all && options.precompute
    precompute_data(options);
end

%% Load data

if options.only_reconstructions == false
    options = load_custom_data(options);
end

%% Reconstructions

[options, pz] = custom_prior_prepass(options);


for t = 1 : options.partitions
    
    for iter = 1 : options.Niter
        
        for osa_iter = 1 : options.subsets
            %%% Your custom prior here
            %%% Replace custom_prior with your own function and uncomment the
            %%% methods that have been selected above
            %%% OSL-OSEM
            % options.grad_OSEM = custom_prior(options.im_vectors.custom_OSEM(:,iter));
            options.grad_OSEM = MRP(options.im_vectors.custom_OSEM(:,iter), options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
            %%% OSL-MLEM (method 2 only)
            % options.grad_MLEM = custom_prior(options.im_vectors.custom_MLEM(:,iter));
            %%% MBSREM
            % options.grad_MBSREM = custom_prior(options.im_vectors.custom_MBSREM(:,iter));
            %%% BSREM
            % options.grad_BSREM = custom_prior(options.im_vectors.custom_BSREM(:,iter));
            %%% ROSEM-MAP
            % options.grad_ROSEM = custom_prior(options.im_vectors.custom_ROSEM(:,iter));
            %%% RBI-MAP
            % options.grad_RBI = custom_prior(options.im_vectors.custom_RBI(:,iter));
            %%% OSL-COSEM
            % options.grad_COSEM = custom_prior(options.im_vectors.custom_COSEM(:,iter));
            options = custom_prior_reconstruction(options, t, iter, osa_iter);
        end
        options = init_next_iter(options, iter);
    end
    [options, pz] = save_custom_prior_iterations(options, t, pz);
end

% save([options.name '_reconstruction_' num2str(options.subsets) 'subsets_' num2str(options.Niter) 'iterations.mat'], 'pz');