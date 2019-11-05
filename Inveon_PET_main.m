%% MATLAB codes for PET reconstruction using Inveon PET LST or GATE output

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% MACHINE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Blocks in transaxial direction
options.blocks_per_ring = (16);
%%% Blocks in axial direction (i.e.  number of physical machine rings)
options.linear_multip = (4);
%%% number of detectors on the side of R-sector(block) (e.g. 13 if 13x13)
options.cryst_per_block = (20);
%%% crystal pitch in x- and y-directions (mm)
options.cr_p = 1.63;
%%% crystal pitch in z-direction (mm)
options.cr_pz = 1.592;
%%% ring diameter (distance between perpendicular detectors) in mm
options.diameter = 161;
%%% Transaxial FOV size (mm), this is the length of the x (horizontal) side
% of the FOV
options.FOVa_x = 100;
%%% Transaxial FOV size (mm), this is the length of the y (vertical) side
% of the FOV
options.FOVa_y = options.FOVa_x;
%%% Axial FOV (mm)
options.axial_fov = 127;
%%% Number of pseudo rings between physical rings (use 0 or [] if none)
options.pseudot = [];
%%% Number of detectors per ring (without pseudo detectors)
options.det_per_ring = options.blocks_per_ring*options.cryst_per_block;
%%% Number of detectors per ring (with pseudo detectors)
options.det_w_pseudo = options.blocks_per_ring*(options.cryst_per_block);
%%% number of crystal rings
options.rings = options.linear_multip * options.cryst_per_block;
%%% number of detectors
options.detectors = options.det_per_ring*options.rings;
%%% machine name
options.machine_name = 'Inveon';
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INVEON DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Use Inveon data
% 0/false = Don't use Inveon machine data (use GATE data instead)
% 1 = Use list-mode data (.lst)
% 2 = Use machine created sinogram data (.scn)
options.use_machine = 1;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% GATE SPECIFIC SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Obtain trues (true coincidences)
% If this is set to true then, in addition to the normal coincidences,
% trues are also obtained and saved
options.obtain_trues = false;
%%% Reconstruct the true coincidences
% If this is set to true, then the true coincidences will be used for
% reconstruction
% NOTE: If both this and reconstruct_scatter are set, then the trues are
% reconstructed, but not the scatter.
options.reconstruct_trues = false;
%%% Obtain scattered coincidences
% If this is set to true, then scattered coincidences are saved separately
% These events are not used for scatter correction though, but a separate
% scatter sinogram/raw data matrix will be created.
options.store_scatter = false;
%%% What scatter components are included in the scatter part
% (1 means that component is included, 0 means it is not included in the
% scatter data) 
% First: Compton scattering in the phantom, second: Compton scattering in
% the detector, third: Rayleigh scattering in the phantom, fourth: Rayleigh
% scattering in the detector
% If store_scatter is set to true, at least one value has to be 1
% NOTE: LMF will always include only Compton scattering in the phantom,
% regardless of the choice below (as long as scatter is selected)
options.scatter_components = [1 1 1 1];
%%% Reconstruct the scattered coincidences
% If this is set to true, then the scattered coincidences will be used for
% reconstruction
% NOTE: If both this and reconstruct_trues are set, then the trues are
% reconstructed, but not the scatter.
options.reconstruct_scatter = false;
%%% Obtain (true) random coincidences.
% If this is set to true then coincidence events that are genuine random
% events are stored separately.
% These events are not used for randoms correction (see the
% Corrections-section for delayed coincidence window randoms correction),
% but a separate randoms sinogram/raw data matrix will be created.
options.store_randoms = false;
%%% Obtain source coordinates (used in forming the "true" image).
% If this is set to true, then the "true" decay image is also saved during
% data load, i.e. the locations where the decay has occurred.
% If any of the above settings are set to true, then the true images are
% also obtained for them. E.g. if store_scatter = true, then an image
% showing the locations and number of counts of where the scattered events
% originated will be saved in a mat-file. Scatter and trues contain
% coincidence events while randoms contain singles.
% NOTE: If you use LMF data, the source images are not considered reliable.
options.source = false;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% ASCII DATA FORMAT SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% Is ASCII data used (Only one data type can be used at a time)
options.use_ASCII = true;
% Copy-paste the ASCII coincidence mask used in your macro inside the
% brackets below. If no coincidence mask is used, use an empty array ([]).
options.coincidence_mask = [0 1 0 1 1 1 1 0 0 0 0 1 1 1 1 1 0 0 0 1 0 1 1 1 1 0 0 0 0 1 1 1 1 1 0 0];
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% LMF DATA FORMAT SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% Is LMF data used (Only one data type can be used at a time)
options.use_LMF = false;
%%% How many bytes are in the LMF header part?
options.header_bytes = (16);
% How many bytes in each event packet?
options.data_bytes = (8 + 2 + 6 + 4 + 1);
% How many bits dedicated to R-sectors?
options.R_bits = 4;
% How many bits dedicated to modules?
options.M_bits = 2;
% How many bits dedicated to submodules?
options.S_bits = 1;
% How many bits dedicated to crystals?
options.C_bits = 9;
% How many bits dedicated to layers?
options.L_bits = 0;
% What is the coincidence window (in seconds)?
options.coincidence_window = 10e-9;
% Obtain source coordinates? (used in forming the "true" image)
options.source = true;
% What is the clock time step? (see the .cch files)
% If e.g. 1 ps then use 1e-12, if 1 ns use 1e-9, etc.
options.clock_time_step = 1e-12;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% ROOT DATA FORMAT SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% Is ROOT data used (Only one data type can be used at a time)
% NOTE: On Windows ROOT works only with 32-bit MATLAB, but has not been
% tested with it
options.use_root = false;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% IMAGE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% Reconstructed image pixel size (X-direction)
options.Nx = 128;
%%% Y-direction
options.Ny = 128;
%%% Z-direction (number of slices)
options.Nz = options.rings*2 - 1;
%%% Flip the image (in vertical direction)?
options.flip_image = false;
%%% How much is the image rotated?
% You need to run the precompute phase again if you modify this
% NOTE: The rotation is done in the detector space (before reconstruction)
options.offangle = options.det_w_pseudo * (2/4) - options.cryst_per_block/2;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SINOGRAM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% Span factor/axial compression
options.span = 3;
% maximum ring difference
options.ring_difference = options.rings - 1;
% number of angular positions in sinogram
options.Ndist = 128;
% number of angles in sinogram (this is the final amount after possible
% mashing), maximum allowed is the number of detectors per ring/2
options.Nang = 160;
% for sinograms, specify the amount of sinograms contained on each segment
% (this should total the total number of sinograms)
options.segment_table = [options.Nz, options.Nz - (options.span + 1):-options.span*2:options.span];
if verLessThan('matlab','8.5')
    options.segment_table = [options.segment_table(1), repeat_elem(options.segment_table(2:end),2,1)];
else
    options.segment_table = [options.segment_table(1), repelem(options.segment_table(2:end),2)];
end
% Total number of sinograms
options.TotSinos = sum(options.segment_table);
% number of sinograms used in reconstruction
options.NSinos = options.TotSinos;
% If Ndist value is even, take one extra out of the negative side (+1) or
% from the positive side (-1). E.g. if Ndist = 200, then with +1 the
% interval is [-99,100] and with -1 [-100,99].
options.ndist_side = -1;
% Fill the gaps caused by pseudo detectors?
% NOTE: Applicable only if options.pseudot > 0
options.fill_sinogram_gaps = false;
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
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CORRECTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%% Randoms correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If set to true, stores the delayed coincidences during data load and
% later corrects for randoms during the data formation/load or during
% reconstruction
options.randoms_correction = false;

%%% Variance reduction
% If set to true, variance reduction will be performed to delayed
% coincidence (randoms corrections) data if randoms correction is selected
options.variance_reduction = false;

%%% Randoms smoothing
% If set to true, applies a 8x8 moving mean smoothing to the delayed
% coincidence data. This is applied on all cases (randoms correction data
% is smoothed before subtraction of before reconstruction)
% NOTE: Mean window size can be adjusted by modifying the randoms_smoothing
% function
options.randoms_smoothing = false;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Scatter correction %%%%%%%%%%%%%%%%%%%%%%%%%%%
% If set to true, will prompt the user to load the scatter sinogram/raw
% data. Corrects for scatter during data formation/load or during
% reconstruction
% NOTE: Scatter data is not created by this software and as such must be
% provided by the user
options.scatter_correction = false;

%%% Scatter normalization
% If set to true, normalizes the scatter coincidences data before
% reconstruction. This applies only if the below
% options.corrections_during_reconstruction = true, otherwise it will have
% no effect (scatter correction is applied before the sinogram/raw data is
% normalized).
options.normalize_scatter = false;

%%% Scatter smoothing
% If set to true, applies a 8x8 moving mean smoothing to the scattered
% coincidences data. This is applied on all cases (scatter correction data
% is smoothed before subtraction of before reconstruction)
% NOTE: Mean window size can be adjusted by modifying the randoms_smoothing
% function
options.scatter_smoothing = false;
 
%%%%%%%%%%%%%%%%%%%%%%%%% Attenuation correction %%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Image-based attenuation correction
% Include attenuation correction from images (e.g. CT-images) (for this you
% need attenuation images of each slice correctly rotated and scaled for
% 511 keV) 
options.attenuation_correction = false;
%%% Attenuation image data file
% specify the path (if not in MATLAB path) and filename
% NOTE: the attenuation data must be the only variable in the file and
% have the dimensions of the final reconstructed image
options.attenuation_datafile = '';
 
%%%%%%%%%%%%%%%%%%%%%%%% Normalization correction %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute the normalization coefficients
% If set to true, then the normalization coefficients are computed after
% the measurement data has been loaded
options.compute_normalization = false;
% Normalization correction components to include (1 means that the
% component is included, 0 that it is not included)
% First: Axial geometric correction 
% Second: Detector efficiency correction, use 1 for fan-sum algorithm (both
% sinogram and raw list-mode data) or 2 for SPC (only raw list-mode data)
% Third: Block profile correction
% Fourth: Transaxial geometric correction (NOT recommended when using
% normalization data that does not encompass the entire FOV)
options.normalization_options = [1 1 1 1];
% If a cylinder that is smaller than the FOV was used for the normalization
% measurement, specify the radius of this cylinder (cm) otherwise use an
% empty array or inf.
options.normalization_phantom_radius = 3.5;
% Apply scatter correction to normalization cylinder
% If cylinder is used for normalization correction, applies also the
% scatter correction. Requires the above cylinder radius.
% NOTE: Applicable only to sinogram data
options.normalization_scatter_correction = true;
 
%%% Apply normalization correction
% If set to true, normalization correction is applied.
options.normalization_correction = false;
%%% Use user-made normalization
% Use either a .mat or .nrm file containing the normalization coefficients
% for normalization correction if normalization_correction is also set to
% true. 
% User will be prompted for the location of the file either during sinogram
% formation or before image reconstruction (see below)
% NOTE: If you have previously computed normalization coefficients with
% OMEGA, you do not need to set this to true. The normalization
% coefficients for the specified machine will be automatically loaded. Use
% this only if you want to use normalization coefficients computed outside
% of OMEGA.
options.use_user_normalization = false;
 
%%%%%%%%%%%%%%%%%%%% Corrections during reconstruction %%%%%%%%%%%%%%%%%%%%
% If set to true, all the corrections are performed during the
% reconstruction step, otherwise the corrections are performed to the
% sinogram/raw data before reconstruction.
% NOTE: Attenuation correction is always performed during reconstruction
% regardless of the choice below
options.corrections_during_reconstruction = false;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% DYNAMIC IMAGING PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% Total time of the measurement (s)
% Use inf if you want the whole examination
options.tot_time = inf;
%%% Number of time points/dynamic frames (if a static measurement, use 1)
options.partitions = 1;
%%% Start time (s) (all measurements before this will be ignored)
options.start = 0;
%%% End time (s) (all measurements after this will be ignored)
options.end = options.tot_time;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MISC PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% Name of current datafile/examination
options.name = 'open_PET_data';
%%% Precompute data 
% This should be done when using data from a certain machine the first time
% as it can speed up reconstruction. Especially recommended for raw
% list-mode data.
options.precompute = true;
%%% Folder for the data (.dat ASCII, .ccs LMF, .root ROOT) files
% If no files are located in the path provided below, then the current
% folder is also checked. If no files are detected there either, an error
% is thrown.
% NOTE: for .lst or .scn files the user will be prompted for their
% locations and as such this path is ignored
if ispc % Windows
    options.fpath = 'C:\path\to\GATE\output\';
else % Unix/Mac
    options.fpath = '/path/to/GATE/output/';
end
%%% Form only sinograms (no reconstructions)
options.only_sinos = false;
%%% Precompute the observation matrix for the reconstruction (this might require a
% lot of memory), if false then the observation matrix is calculated on the
% fly
% NOTE: Supports only MLEM reconstruction and is an experimental feature
options.precompute_obs_matrix = false;
%%% Compute only the reconstructions
% If this file is run with this set to true, then the data load and
% sinogram formation steps are always skipped. Precomputation step is
% only performed if precompute_lor = true and precompute_all = true.
options.only_reconstructions = false;
%%% Use raw list mode data
% This means that the data is used as is without any sinogramming and thus
% without any "compression"
options.use_raw_data = false;
%%% Use precomputed geometrical matrix information
% During the precompute-phase the number of pixels each LOR traverse is
% counted (this phase requires the above precompute-option to true). These
% are saved and later on used during the reconstruction process. Once the
% precompute-phase has been done once for the specific sinogram and image
% resolution, it is not necessary to do it again. Recommended for raw
% list-mode data.
options.precompute_lor = false;
%%% Precompute all
% Set this option to true if you want to precompute all possible
% combinations in one go (i.e. raw data, precomputed LORs, sinogram format)
% Requires for precompute option to be true to take effect
% Setting this option to true also causes the precompute phase to be
% performed even if only_reconstructions is true
options.precompute_all = false;
% Show status messages
options.verbose = true;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPLEMENTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reconstruction implementation used
% 1 = Reconstructions in MATLAB (projector in a MEX-file)
% 2 = Matrix-free reconstruction with OpenCL/ArrayFire (Recommended)
% 3 = Multi-GPU/device matrix-free OpenCL (OSEM & MLEM only)
% 4 = Matrix-free reconstruction with OpenMP (parallel) or C++ (sequential) 
% (OSEM & MLEM only)
options.implementation = 4;
% Implementations 2 and 3 ONLY
%%% Device used (this is applicable to implementation 2), or platform used
% (implementation 3)
% In implementation 2 this determines the device used for both system
% matrix formation and image reconstruction
% NOTE: Use ArrayFire_OpenCL_device_info() to determine the device numbers
% In implementation 3, this determines the platform from where the
% device(s) are taken
% NOTE: Use OpenCL_device_info() to determine the platform numbers and
% their respective devices
% NOTE: if you switch devices then you need to run the below line
% (uncommented) as well:
% clear mex
options.use_device = 0;
%%% Use 64-bit integer atomic functions
% If true, then 64-bit integer atomic functions (atomic add) will be used
% if they are supported by the selected device
% Setting this to true will make computations faster on GPUs that support
% the functions, but might make results slightly less reliable
options.use_64bit_atomics = true;
% Implementation 2 ONLY
%%% Force the (re)building of OpenCL binaries 
% If set to true, the OpenCL binaries are rebuilt even if they have been
% previously built
% Use this once if you update your drivers or there are changes made to the
% .cl-files
options.force_build = false;
% Implementation 3 ONLY
%%% How many times more measurements/LORs are in the GPU part (applicable if
% heterogenous computing is used) 
% Alternatively, set this to 0 to use only a single device on the specific
% platform (the one with the highest memory count will be used)
options.cpu_to_gpu_factor = 2.5;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROJECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Type of projector to use for the geometric matrix
% 0 = Regular Siddon's algorithm (only available with implementation 1 and
% when precomputed_lor = false.
% 1 = Improved/accelerated Siddon's algorithm
% 2 = Orthogonal distance based ray tracer
options.projector_type = 1;
% Orthogonal ray tracer only
%%% The 2D (XY) width of the "strip/tube" where the orthogonal distances are
% included. If the tube_width_z is non-zero, then this value is ignored.
options.tube_width_xy = options.cr_p;
% Orthogonal ray tracer only
%%% The 3D (Z) width of the "tube" where the orthogonal distances are
% included. If set to 0, then the 2D orthogonal ray tracer is used. If this
% value is non-zero then the above value is IGNORED.
options.tube_width_z = options.cr_pz;
% 3D Orthogonal ray tracer only
%%% Accuracy factor
% This value controls the approximations in the 3D orthogonal ray tracer.
% Higher values lead to more accurate results, but slower computational
% speeds. Default value is 5 and is a compromise between accuracy and
% speed. Values above the X/Y pixel count have no effect.
options.accuracy_factor = 5;
%%% Number of rays
% Number of rays used if projector_type = 1 (i.e. Improved Siddon is used)
options.n_rays = 5;
 
%%%%%%%%%%%%%%%%%%%%%%%%% RECNSTRUCTION SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Number of iterations (all reconstruction methods)
options.Niter = 4;
%%% Number of subsets (all excluding MLEM and subset_type = 5)
options.subsets = 8;
%%% Subset type (n = subsets)
% 1 = Every nth (column) measurement is taken
% 2 = Every nth (row) measurement is taken (e.g. if subsets = 3, then
% first subset has measurements 1, 4, 7, etc., second 2, 5, 8, etc.) 
% 3 = Measurements are selected randomly
% 4 = (Sinogram only) Take every nth column in the sinogram
% 5 = (Sinogram only) Take every nth row in the sinogram
% 6 = Sort the LORs according to their angle with positive X-axis, combine
% n_angles together and have 180/n_angles subsets for 2D slices and
% 360/n_angles for 3D, see GitHub wiki for more information
% 7 = Form the subsets by using golden angle sampling
options.subset_type = 1;
%%% How many angles are combined in subset_type = 6
% E.g. there are 180 angles, in n_angles = 2, then angles 0 and 1 are
% combined to the same subset, 2 and 3, etc.
options.n_angles = 2;
%%% Initial value for the reconstruction
options.x0 = ones(options.Nx, options.Ny, options.Nz);
%%% Epsilon value 
% small value to prevent division by zero
options.epps = 1e-8;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MISC SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Use Shuffle (recommended)
% NOTE: Applies only when using subset_type = 3; accelerates the subset
% formation and uses less memory. Not included in OMEGA, needs to be
% manually downloaded and installed.
% Download from: 
% https://se.mathworks.com/matlabcentral/fileexchange/27076-shuffle
options.use_Shuffle = true;
%%% Use fast sparse
% Not included in OMEGA, needs to be manually downloaded and installed.
% Download from: https://github.com/stefanengblom/stenglib
% NOTE: This applies only to implementation 1 when precompute_lor is false.
options.use_fsparse = true;
%%% Skip the normalization phase in MRP, FMH, L-filter, ADMRP and/or weighted
% mean
% E.g. if set to true the MRP prior is (x - median(x))
% E.g. if set to false the MRP prior is (x - median(x)) / median(x)
options.med_no_norm = false;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION ALGORITHMS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction algorithms to use (you can choose several)
% NOTE: MLEM requires precomputed observation matrix or a matrix-free
% method
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ML-METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are non-regularized versions
%%% Maximum-Likelihood Expectation Maximization (MLEM)
% Supported by all implementations
options.mlem = false;
%%% Ordered Subsets Expectation Maximization (OSEM)
% Supported by all implementations
options.osem = true;
%%% Modified Row-Action Maximum Likelihood Algorithm (MRAMLA)
% Supported by implementations 1 and 2
options.mramla = false;
%%% Row-Action Maximum Likelihood Algorithm (RAMLA)
% Supported by implementations 1, 2 and 4
options.ramla = false;
%%% Relaxed Ordered Subsets Expectation Maximization (ROSEM)
% Supported by implementations 1, 2 and 4
options.rosem = false;
%%% Rescaled Block Iterative Expectation Maximization (RBI-EM)
% Supported by implementations 1 and 2
options.rbi = false;
%%% Dynamic RAMLA (DRAMA)
% Supported by implementations 1 and 2
options.drama = false;
%%% Complete data OSEM (COSEM)
% Supported by implementations 1 and 2
options.cosem = false;
%%% Enhanced COSEM (ECOSEM)
% Supported by implementations 1 and 2
options.ecosem = false;
%%% Accelerated COSEM (ACOSEM)
% Supported by implementations 1 and 2
options.acosem = false;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAP-METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Any algorithm selected here will utilize all the priors selected below.
% For example, if OSL-OSEM is set to true and MRP and Quad are set to true,
% then OSL-OSEM estimates will be computed for both MRP and Quadratic
% prior.
%%% One-Step Late MLEM (OSL-MLEM)
% Supported by implementations 2 and 4
options.OSL_MLEM = false;
%%% One-Step Late OSEM (OSL-OSEM)
% Supported by implementations 1, 2 and 4
options.OSL_OSEM = false;
%%% Modified BSREM (MBSREM)
% Supported by implementations 1 and 2
options.MBSREM = false;
%%% Block Sequential Regularized Expectation Maximation (BSREM)
% Supported by implementations 1, 2 and 4
options.BSREM = false;
%%% ROSEM-MAP
% Supported by implementations 1, 2 and 4
options.ROSEM_MAP = false;
%%% RBI-MAP
% Supported by implementations 1 and 2
options.RBI_MAP = false;
%%% (A)COSEM-MAP (OSL)
% 0/false = No COSEM-MAP, 1/true = ACOSEM-MAP, 2 = COSEM-MAP
% Supported by implementations 1 and 2
options.COSEM_MAP = false;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Median Root Prior (MRP)
options.MRP = false;
%%% Quadratic Prior (QP)
options.quad = false;
%%% L-filter prior
options.L = false;
%%% Finite impulse response (FIR) Median Hybrid (FMH) prior
options.FMH = false;
%%% Weighted mean prior
options.weighted_mean = false;
%%% Total Variation (TV) prior
options.TV = false;
%%% Anisotropic Diffusion Median Root Prior (ADMRP)
options.AD = false;
%%% Asymmetric Parallel Level Set (APLS) prior
options.APLS = false;
%%% Total Generalized Variation (TGV) prior
options.TGV = false;
%%% Non-local Means (NLM) prior (implementations 1 and 4 only)
options.NLM = false;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% ACOSEM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accleration parameter for ACOSEM (1 equals COSEM)
options.h = 2;

%%%%%%%%%%%%%%%%%%%%%%%% MRAMLA & MBSREM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%
%%% Relaxation parameter for MRAMLA and MBSREM
% Use scalar if you want it to decrease as lambda0_mbsrem/current_iteration
% Use vector (length = Niter) if you want your own relaxation parameters
options.lambda0_mbsrem = 0.2;
% Upper bound for MRAMLA/MBSREM (use 0 for default value)
options.U = 0;

%%%%%%%%%%%%%%%%%%%%%%%%% RAMLA & BSREM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%
%%% Relaxation parameter for RAMLA and BSREM
% Use scalar if you want it to decrease as lambda0/current_iteration
% Use vector (length = Niter) if you want your own relaxation parameters
options.lambda0 = 0.1;

%%%%%%%%%%%%%%%%%%%%%%% ROSEM & ROSEM-MAP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%
%%% Relaxation parameter for ROSEM and ROSEM-MAP
% Use scalar if you want it to decrease as lambda0_rosem/current_iteration
% Use vector (length = Niter) if you want your own relaxation parameters
options.lambda0_rosem = 0.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%% DRAMA PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Beta_0 value
options.beta0_drama = 0.1;
% Beta value
options.beta_drama = 1;
% Alpha value
options.alpha_drama = 0.1;


%%%%%%%%%%%%%%%%%%%%%%%%% NEIGHBORHOOD PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%
%%% How many neighboring pixels are considered 
% L, FMH and weighted mean)
% E.g. if Ndx = 1, Ndy = 1, Ndz = 0, then you have 3x3 square area where
% the pixels are taken into account (I.e. (Ndx*2+1)x(Ndy*2+1)x(Ndz*2+1)
% area)
% NOTE: At the moment Ndx and Ndy has to be identical
options.Ndx = 1;
options.Ndy = 1;
options.Ndz = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MRP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for MRP with OSL-OSEM
options.beta_mrp_osem = 0.1;%0.3;
% regularization parameter for MRP with OSL-MLEM
options.beta_mrp_mlem = 1.5;
% regularization parameter for MRP with MBSREM
options.beta_mrp_mbsrem = 0.3;
% regularization parameter for MRP with BSREM
options.beta_mrp_bsrem = 0.1;
% regularization parameter for MRP with ROSEM
options.beta_mrp_rosem = 2;
% regularization parameter for MRP with RBI
options.beta_mrp_rbi = 0.1;
% Regularization parameter for MRP with OSL-(A)COSEM
options.beta_mrp_cosem = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for quadratic prior with OSL-OSEM
options.beta_quad_osem = 0.01;%0.1;
% regularization parameter for quadratic prior with OSL-MLEM
options.beta_quad_mlem = 0.1;
% regularization parameter for quadratic prior with MBSREM
options.beta_quad_mbsrem = 0.05;
% regularization parameter for quadratic prior with BSREM
options.beta_quad_bsrem = 0.03;
% regularization parameter for quadratic prior with ROSEM
options.beta_quad_rosem = 0.1;
% regularization parameter for quadratic prior with RBI
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


%%%%%%%%%%%%%%%%%%%%%%%%%%% L-FILTER PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for L-filter with OSL-OSEM
options.beta_L_osem = 0.1;
% regularization parameter for L-filter with OSL-MLEM
options.beta_L_mlem = 0.1;
% regularization parameter for L-filter with MBSREM
options.beta_L_mbsrem = 0.1;
% regularization parameter for L-filter with BSREM
options.beta_L_bsrem = 0.03;
% regularization parameter for L-filter with ROSEM
options.beta_L_rosem = 3;
% regularization parameter for L-filter with RBI
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
% regularization parameter for FMH with OSL-OSEM
options.beta_fmh_osem = 0.1;
% regularization parameter for FMH with OSL-MLEM
options.beta_fmh_mlem = 0.1;
% regularization parameter for FMH with MBSREM
options.beta_fmh_mbsrem = 0.6;
% regularization parameter for FMH with BSREM
options.beta_fmh_bsrem = 5;
% regularization parameter for FMH with ROSEM
options.beta_fmh_rosem = 8;
% regularization parameter for FMH with RBI
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
% regularization parameter for weighted mean with OSL-OSEM
options.beta_weighted_osem = 0.02;%0.2;
% regularization parameter for weighted mean with OSL-MLEM
options.beta_weighted_mlem = 0.1;
% regularization parameter for weighted mean with MBSREM
options.beta_weighted_mbsrem = 0.1;
% regularization parameter for weighted mean with BSREM
options.beta_weighted_bsrem = 5;
% regularization parameter for weighted mean with ROSEM
options.beta_weighted_rosem = 3;
% regularization parameter for weighted mean with RBI
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
options.beta_TV_osem = 0.03;
% regularization parameter for TV with OSL-MLEM
options.beta_TV_mlem = 0.1;
% regularization parameter for TV with MBSREM
options.beta_TV_mbsrem = 0.1;
% regularization parameter for TV with BSREM
options.beta_TV_bsrem = 0.05;
% regularization parameter for TV with ROSEM
options.beta_TV_rosem = 0.7;
% regularization parameter for TV with RBI
options.beta_TV_rbi = 0.02;
% Regularization parameter for TV (OSL-(A)COSEM)
options.beta_TV_cosem = 0.03;
% "Smoothing" parameter
% Also used to prevent zero values in square root
options.TVsmoothing = 1e-1;
% Whether to use an anatomical reference/weighting image with the TV
options.TV_use_anatomical = false;
% If the TV_use_anatomical is set to true, you can specify the type of
% anatomical TV that is used. Three different types are available, see the
% wiki and readme for algorithm details
% Value can be 1, 2 or 3
options.TVtype = 1;
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
% regularization parameter for AD with OSL-MLEM
options.beta_ad_mlem = 0.1;
% regularization parameter for AD with MBSREM
options.beta_ad_mbsrem = 0.3;
% regularization parameter for AD with BSREM
options.beta_ad_bsrem = 0.2;
% regularization parameter for AD with ROSEM
options.beta_ad_rosem = 3;
% regularization parameter for AD with RBI
options.beta_ad_rbi = 0.05;
% regularization parameter for AD with (OSL-(A)COSEM)
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
% regularization parameter for APLS with OSL-MLEM
options.beta_APLS_mlem = 0.1;
% regularization parameter for APLS with MBSREM
options.beta_APLS_mbsrem = 0.1;
% regularization parameter for APLS with BSREM
options.beta_APLS_bsrem = 0.005;
% regularization parameter for APLS with ROSEM
options.beta_APLS_rosem = 0.1;
% regularization parameter for APLS with RBI
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
% regularization parameter for TGV with OSL-MLEM
options.beta_TGV_mlem = 0.1;
% regularization parameter for TGV with MBSREM
options.beta_TGV_mbsrem = 0.1;
% regularization parameter for TGV with BSREM
options.beta_TGV_bsrem = 1;
% regularization parameter for TGV with ROSEM
options.beta_TGV_rosem = 0.25;
% regularization parameter for TGV with RBI
options.beta_TGV_rbi = 0.1;
% Regularization parameter for TGV (OSL-(A)COSEM)
options.beta_TGV_cosem = 0.05;
% "Smoothing" parameter 1 (eta)
% Also used to prevent zero values in square root
options.alphaTGV = 2;
options.betaTGV = 1;
options.NiterTGV = 30;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NLM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for NLM with OSL-OSEM
options.beta_NLM_osem = 0.001;
%%% Regularization parameter for NLM with MBSREM
options.beta_NLM_mbsrem = 0.05;
%%% Regularization parameter for NLM with BSREM
options.beta_NLM_bsrem = 0.01;
%%% Regularization parameter for NLM with ROSEM
options.beta_NLM_rosem = 0.1;
%%% Regularization parameter for NLM with RBI
options.beta_NLM_rbi = 0.01;
%%% Regularization parameter for NLM (OSL-(A)COSEM)
options.beta_NLM_cosem = 0.01;
%%% Filter parameter
options.sigma = 0.01;
%%% Patch radius
options.Nlx = 2;
options.Nly = 2;
options.Nlz = 1;
%%% Standard deviation of the Gaussian filter
options.NLM_gauss = 1;
% Search window radius is controlled by Ndx, Ndy and Ndz parameters
% Use anatomical reference image for the patches
options.NLM_use_anatomical = false;
%%% Specify filename for the reference image here (same rules apply as with
% attenuation correction above)
options.NLM_reference_image = 'reference_image.mat';
%%% Use Non-local total variation (NLTV)
% If selected, will overwrite regular NLM regularization as well as the
% below MRP version
options.NLTV = false;
%%% Use MRP algorithm (without normalization)
% I.e. gradient = im - NLM_filtered(im)
options.NLM_MRP = false;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%% OpenCL device info %%%%%%%%%%%%%%%%%%%% 
% Uncomment the below line and run it to determine the available devices
% and their respective numbers
% OpenCL_device_info()







%%%%%%%%%%%%%%%%%%%% Load scatter data %%%%%%%%%%%%%%%%%%%% 
% Load scatter data (if applicable)
if options.scatter_correction && ~options.only_reconstructions || options.scatter_correction && options.use_raw_data
    options = loadScatterData(options);
end





%%%%%%%%%%%%%%%%%%%% Error checking %%%%%%%%%%%%%%%%%%%%
% Basic error checking is done here
options = OMEGA_error_check(options);


%% Precompute data if necessary

if options.precompute && options.only_reconstructions == false && options.precompute_lor || options.precompute_all && options.precompute && options.precompute_lor
    precompute_data(options);
end

%% Load the list-mode coindicence data

if options.only_reconstructions == false && options.use_machine ~= 2 || options.precompute_all && ~options.only_reconstructions && options.use_machine ~= 2
    options.coincidences = load_data(options);
end

%% Form the sinograms

if options.only_reconstructions == false && (options.use_raw_data == false || options.use_machine == 2) || options.precompute_all && ~options.only_reconstructions
    options.SinM = form_sinograms(options);
end


%% Form the attenuation correction image from Inveon atn-file

if options.attenuation_correction && ~options.only_reconstructions
    options.attenuation_datafile = attenuation_correction_factors(options);
end

%% Compute the normalization coefficients

if options.compute_normalization && ~options.only_reconstructions
    [norm_matrix, options.SinM, axial_geom_coeffs, cross_coefs_save, block_profile_matrix, coeff_matrix, det_coeffs, tr_geom_matrix] = normalization_coefficients(options);
end

%% Reconstructions

if options.only_sinos == false
    
    tStart = tic;
    pz = reconstructions_main(options);
    tElapsed = toc(tStart);
    disp(['Reconstruction process took ' num2str(tElapsed) ' seconds'])
    
end

save([options.name '_reconstruction_' num2str(options.subsets) 'subsets_' num2str(options.Niter) 'iterations_' ...
    num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '_NLMTV.mat'], 'pz');