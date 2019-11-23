%% MATLAB codes for PET reconstruction using sinogram input from any machine

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% MACHINE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% R-sectors/blocks in transaxial direction
options.blocks_per_ring = (42);
%%% R-sectors/modules/blocks in axial direction (i.e. number of physical
% machine rings). Multiplying this with the below cryst_per_block should
% equal the number of crystal rings
options.linear_multip = (4);
%%% number of detectors on the side of R-sector/block/module (transaxial
%%% direction)
% (e.g. 13 if 13x13, 20 if 20x10)
options.cryst_per_block = (8);
%%% crystal pitch in x- and y-directions (mm)
options.cr_p = 2.4;
%%% crystal pitch in z-direction (mm)
options.cr_pz = 2.4;
%%% ring diameter (distance between perpendicular detectors) in mm
options.diameter = 130*2;
%%% Transaxial FOV size (mm), this is the length of the x (horizontal) side
% of the FOV
options.FOVa_x = 151;
%%% Transaxial FOV size (mm), this is the length of the y (vertical) side
% of the FOV
options.FOVa_y = options.FOVa_x;
%%% Axial FOV (mm)
options.axial_fov = floor(76.8 - options.cr_pz/10);
%%% Number of pseudo rings between physical rings (use 0 or [] if none)
options.pseudot = [];
%%% Number of detectors per crystal ring (without pseudo detectors)
options.det_per_ring = options.blocks_per_ring*options.cryst_per_block;
%%% Number of detectors per crystal ring (with pseudo detectors)
% If your scanner has a single pseudo detector on each side of the crystal
% block then simply add +1 inside the parenthesis
options.det_w_pseudo = options.blocks_per_ring*(options.cryst_per_block);
%%% Number of crystal rings
options.rings = options.linear_multip * options.cryst_per_block;
%%% Total number of detectors
options.detectors = options.det_per_ring*options.rings;
%%% Machine name
options.machine_name = 'Cylindrical_PET_example';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% IMAGE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Reconstructed image pixel count (X-direction)
% NOTE: Non-square image sizes (X- and Y-direction) may not work
options.Nx = 128;
%%% Y-direction
options.Ny = 128;
%%% Z-direction (number of slices)
options.Nz = options.rings*2-1;
%%% Flip the image (in vertical direction)?
options.flip_image = false;
%%% How much is the image rotated?
% You need to run the precompute phase again if you modify this
% NOTE: The rotation is done in the detector space (before reconstruction)
options.offangle = options.det_w_pseudo * (3/4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SINOGRAM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Span factor/axial compression
options.span = 3;
%%% Maximum ring difference
options.ring_difference = 31;
%%% Number of angles (tangential positions) in sinogram
% This is the final amount after possible mashing, maximum allowed is the
% number of detectors per ring/2
options.Nang = options.det_per_ring/2;
%%% Number of angular positions (views) in sinogram
options.Ndist = 200;
%%% Specify the amount of sinograms contained on each segment
% (this should total the total number of sinograms)
options.segment_table = [options.Nz, options.Nz - (options.span + 1):-options.span*2:max(options.Nz - options.ring_difference*2, options.span)];
if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
    options.segment_table = [options.segment_table(1), repeat_elem(options.segment_table(2:end),2,1)];
else
    options.segment_table = [options.segment_table(1), repelem(options.segment_table(2:end),2)];
end
%%% Total number of sinograms
options.TotSinos = sum(options.segment_table);
%%% Number of sinograms used in reconstruction
options.NSinos = options.TotSinos;
%%% If Ndist value is even, take one extra out of the negative side (+1) or
% from the positive side (-1). E.g. if Ndist = 200, then with +1 the
% interval is [-99,100] and with -1 [-100,99].
options.ndist_side = 1;
%%% Increase the sampling rate of the sinogram
% Increasing this interpolates additional rows to the sinogram
% Can be used to prevent aliasing artifacts
% NOTE: Has to be either 1 or divisible by two
options.sampling = 1;
%%% Interpolation method used for sampling rate increase
% All the methods are available that are supported by interp1
options.sampling_interpolation_method = 'linear';
%%% Fill the gaps caused by pseudo detectors?
% NOTE: Applicable only if options.pseudot > 0
options.fill_sinogram_gaps = false;
%%% Which method used to fill the gaps?
% Either MATLAB's built-in fillmissing or inpaint_nans from file exchange
% For inpaint_nans see:
% https://se.mathworks.com/matlabcentral/fileexchange/4551-inpaint_nans
% See wiki for more details
options.gap_filling_method = 'fillmissing';
%%% Interpolation method used with fillmissing
% Possible methods are those listed under method-section in fillmissing
options.interpolation_method_fillmissing = 'linear';
%%% Interpolation method used with inpaint_nans
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
% If set to true, then the variance reduction will be performed to delayed
% coincidence (randoms corrections) data if is selected
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
options.normalization_options = [1 1 1 0];
% If a cylinder that is smaller than the FOV was used for the normalization
% measurement, specify the radius of this cylinder (cm) otherwise use an
% empty array or inf.
options.normalization_phantom_radius = inf;
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
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Arc correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Apply arc correction
% NOTE: Arc correction is an experimental feature. It is currently
% relatively slow and supports only sinogram data. Generally it is not
% recommended to use arc correction (Inveon data is an exception).
% Uses parallel computing toolbox if it is available (parfor)
options.arc_correction = false;
%%% Arc correction interpolation method
% The interpolation method used to interpolate the arc corrected sinogram.
% Available methods are those supported by scatteredInterpolant and
% griddata. If an interpolation method is used which is not supported by
% scatteredInterpolant then griddata will be used instead
% NOTE: griddata is used if scatteredInterpolant is not found
options.arc_interpolation = 'linear';

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
% Use inf if you want the whole examination (static measurement only)
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
options.name = 'cylpet_example';
%%% Precompute data
% This should be done when using data from a certain machine the first time
% as it can speed up reconstruction. Especially recommended for raw
% list-mode data.
options.precompute = false;
%%% Compute only the system matrix (no reconstructions)
options.only_system_matrix = false;
%%% Compute single reconstruction
% Compute e.g. OSEM, ROSEM, OSL-MRP by using their corresponding functions
options.single_reconstructions = false;
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
%%% Show status messages
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
% heterogeneous computing is used)
% Alternatively, set this to 0 to use only a single device on the specific
% platform (the one with the highest memory count will be used)
options.cpu_to_gpu_factor = 2.4;

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
options.n_rays = 1;

%%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%
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
options.use_Shuffle = false;
%%% Use fast sparse
% Not included in OMEGA, needs to be manually downloaded and installed.
% Download from: https://github.com/stefanengblom/stenglib
% NOTE: This applies only to implementation 1 when precompute_lor is false.
options.use_fsparse = false;
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
% Supported by implementations 1, 2, 3, and 5
options.mlem = false;
%%% Ordered Subsets Expectation Maximization (OSEM)
% Supported by all implementations
options.osem = true;
%%% Modified Row-Action Maximum Likelihood Algorithm (MRAMLA)
% Supported by implementations 1 and 2
options.mramla = false;
%%% Row-Action Maximum Likelihood Algorithm (RAMLA)
% Supported by implementations 1 and 2
options.ramla = false;
%%% Relaxed Ordered Subsets Expectation Maximization (ROSEM)
% Supported by implementations 1 and 2
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
% Supported by implementation 2 only
options.OSL_MLEM = false;
%%% One-Step Late OSEM (OSL-OSEM)
% Supported by implementations 1 and 2
options.OSL_OSEM = false;
%%% Modified BSREM (MBSREM)
% Supported by implementations 1 and 2
options.MBSREM = false;
%%% Block Sequential Regularized Expectation Maximization (BSREM)
% Supported by implementations 1 and 2
options.BSREM = false;
%%% ROSEM-MAP
% Supported by implementations 1 and 2
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
%%% Non-local Means (NLM) prior (implementation 1 only)
options.NLM = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%% ACOSEM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Acceleration parameter for ACOSEM (1 equals COSEM)
options.h = 2;

%%%%%%%%%%%%%%%%%%%%%%%% MRAMLA & MBSREM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%
%%% Relaxation parameter for MRAMLA and MBSREM
% Use scalar if you want it to decrease as lambda0_mbsrem/current_iteration
% Use vector (length = Niter) if you want your own relaxation parameters
options.lambda0_mbsrem = 0.2;
%%% Upper bound for MRAMLA/MBSREM (use 0 for default value)
options.U = 0;

%%%%%%%%%%%%%%%%%%%%%%%%% RAMLA & BSREM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%
%%% Relaxation parameter for RAMLA and BSREM
% Use scalar if you want it to decrease as lambda0/current_iteration
% Use vector (length = Niter) if you want your own relaxation parameters
options.lambda0 = 0.2;

%%%%%%%%%%%%%%%%%%%%%%% ROSEM & ROSEM-MAP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%
%%% Relaxation parameter for ROSEM and ROSEM-MAP
% Use scalar if you want it to decrease as lambda0_rosem/current_iteration
% Use vector (length = Niter) if you want your own relaxation parameters
options.lambda0_rosem = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%% DRAMA PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Beta_0 value
options.beta0_drama = 0.1;
%%% Beta value
options.beta_drama = 1;
%%% Alpha value
options.alpha_drama = 0.1;


%%%%%%%%%%%%%%%%%%%%%%%%% NEIGHBORHOOD PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%
%%% How many neighboring pixels are considered
% With MRP, QP, L, FMH, NLM and weighted mean
% E.g. if Ndx = 1, Ndy = 1, Ndz = 0, then you have 3x3 square area where
% the pixels are taken into account (I.e. (Ndx*2+1)x(Ndy*2+1)x(Ndz*2+1)
% area)
% NOTE: Currently Ndx and Ndy must be identical
options.Ndx = 1;
options.Ndy = 1;
options.Ndz = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MRP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for MRP with OSL-OSEM
options.beta_mrp_osem = 0.1;%0.3;
%%% Regularization parameter for MRP with OSL-MLEM
options.beta_mrp_mlem = 1.5;
%%% Regularization parameter for MRP with MBSREM
options.beta_mrp_mbsrem = 0.3;
%%% Regularization parameter for MRP with BSREM
options.beta_mrp_bsrem = 0.1;
%%% Regularization parameter for MRP with ROSEM
options.beta_mrp_rosem = 2;
%%% Regularization parameter for MRP with RBI
options.beta_mrp_rbi = 0.1;
%%% Regularization parameter for MRP with OSL-(A)COSEM
options.beta_mrp_cosem = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for quadratic prior with OSL-OSEM
options.beta_quad_osem = 0.01;%0.1;
%%% Regularization parameter for quadratic prior with OSL-MLEM
options.beta_quad_mlem = 0.1;
%%% Regularization parameter for quadratic prior with MBSREM
options.beta_quad_mbsrem = 0.05;
%%% Regularization parameter for quadratic prior with BSREM
options.beta_quad_bsrem = 0.03;
%%% Regularization parameter for quadratic prior with ROSEM
options.beta_quad_rosem = 0.1;
%%% Regularization parameter for quadratic prior with RBI
options.beta_quad_rbi = 0.05;
%%% Regularization parameter for quadratic prior (OSL-(A)COSEM)
options.beta_quad_cosem = 0.01;
%%% Pixel weights for quadratic prior
% the number of pixels need to be the amount of neighboring pixels
% e.g. if the above Nd values are all 1, then 27 weights need to be included
% where the center pixel (if Nd values are 1, element 14) should be Inf
% Size is (Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)
% If left empty then they will be calculated by the algorithm
options.weights = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%% L-FILTER PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for L-filter with OSL-OSEM
options.beta_L_osem = 0.1;
%%% Regularization parameter for L-filter with OSL-MLEM
options.beta_L_mlem = 0.1;
%%% Regularization parameter for L-filter with MBSREM
options.beta_L_mbsrem = 0.1;
%%% Regularization parameter for L-filter with BSREM
options.beta_L_bsrem = 0.03;
%%% Regularization parameter for L-filter with ROSEM
options.beta_L_rosem = 3;
%%% Regularization parameter for L-filter with RBI
options.beta_L_rbi = 0.09;
%%% Regularization parameter for L-filter (OSL-(A)COSEM)
options.beta_L_cosem = 0.1;
%%% Weighting factors for the L-filter pixels
% Otherwise the same as in quadratic prior, but center pixel is not Inf
% If left empty then they will be calculated by the algorithm such that the
% weights resemble a Laplace distribution
options.a_L = [];
%%% If the weighting factors are set empty, then this option will determine
% whether the computed weights follow a 1D weighting scheme (true) or 2D
% (false)
% See the wiki for more information
options.oneD_weights = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FMH PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for FMH with OSL-OSEM
options.beta_fmh_osem = 0.1;
%%% Regularization parameter for FMH with OSL-MLEM
options.beta_fmh_mlem = 0.1;
%%% Regularization parameter for FMH with MBSREM
options.beta_fmh_mbsrem = 0.6;
%%% Regularization parameter for FMH with BSREM
options.beta_fmh_bsrem = 5;
%%% Regularization parameter for FMH with ROSEM
options.beta_fmh_rosem = 8;
%%% Regularization parameter for FMH with RBI
options.beta_fmh_rbi = 0.5;
%%% Regularization parameter for FMH (OSL-(A)COSEM)
options.beta_fmh_cosem = 0.1;
%%% Pixel weights for FMH
% The matrix size needs to be [Ndx*2+1, 4] if Nz = 1 or Ndz = 0, or
% [Ndx*2+1, 13] otherwise
% The center pixel weight should be in the middle
% If the sum of each column is > 1, then the weights will be normalized
% such that the sum = 1
% If left empty then they will be calculated by the algorithm such that the
% weights follow the same pattern as in the original article
options.fmh_weights = [];
%%% Weighting value for the center pixel
% Default value is 4, which was used in the original article
% NOTE: This option is ignored if you provide your own weights
options.fmh_center_weight = 4;


%%%%%%%%%%%%%%%%%%%%%%%%% WEIGHTED MEAN PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%
%%% Mean type
% 1 = Arithmetic mean, 2 = Harmonic mean, 3 = Geometric mean
options.mean_type = 1;
%%% Regularization parameter for weighted mean with OSL-OSEM
options.beta_weighted_osem = 0.2;%0.2;
%%% Regularization parameter for weighted mean with OSL-MLEM
options.beta_weighted_mlem = 0.1;
%%% Regularization parameter for weighted mean with MBSREM
options.beta_weighted_mbsrem = 0.1;
%%% Regularization parameter for weighted mean with BSREM
options.beta_weighted_bsrem = 5;
%%% Regularization parameter for weighted mean with ROSEM
options.beta_weighted_rosem = 3;
%%% Regularization parameter for weighted mean with RBI
options.beta_weighted_rbi = 0.04;
%%% Regularization parameter for weighted mean (OSL-(A)COSEM)
options.beta_weighted_cosem = 0.2;
%%% Pixel weights for weighted mean
% the number of pixels need to be the amount of neighboring pixels
% e.g. if the above Ndx/y/z values are all 1, then 27 weights need to be included
% Size is (Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)
% If left empty then they will be calculated by the algorithm such that the
% weights are dependent on the distance from the center pixel to the
% neighboring pixels
options.weighted_weights = [];
%%% Center pixel weight for weighted mean
% NOTE: This option is ignored if you provide your own weights
options.weighted_center_weight = 4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TV PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for TV OSL-OSEM
options.beta_TV_osem = 0.1;
%%% Regularization parameter for TV with OSL-MLEM
options.beta_TV_mlem = 0.1;
%%% Regularization parameter for TV with MBSREM
options.beta_TV_mbsrem = 0.01;
%%% Regularization parameter for TV with BSREM
options.beta_TV_bsrem = 0.05;
%%% Regularization parameter for TV with ROSEM
options.beta_TV_rosem = 0.07;
%%% Regularization parameter for TV with RBI
options.beta_TV_rbi = 0.002;
%%% Regularization parameter for TV (OSL-(A)COSEM)
options.beta_TV_cosem = 0.003;
%%% "Smoothing" parameter
% Also used to prevent zero values in square root
options.TVsmoothing = 1e-1;
%%% Whether to use an anatomical reference/weighting image with the TV
options.TV_use_anatomical = false;
%%% If the TV_use_anatomical value is set to true, specify filename for the
% reference image here (same rules apply as with attenuation correction
% above)
options.TV_reference_image = 'reference_image.mat';
%%% Three different TV methods are available.
% Value can be 1, 2 or 3.
% Types 1 and 2 are the same if anatomical prior is not included
% Type 3 uses the same weights as quadratic prior
% See the wiki for more information
options.TVtype = 3;
%%% Weighting parameters for the TV prior. Applicable only if
% use_anatomical = true
% T-value is specific to the used TVtype, e.g. for type 1 it is the edge
% threshold parameter. See the wiki for more details
options.T = 0.1;
%%% C is the weight for the original image in type 3 and is ignored with
% other types
options.C = 1;
%%% Tuning parameter for TV and APLS
options.tau = 1e-8;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADMRP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for AD with OSL-OSEM
options.beta_ad_osem = 0.1;
%%% Regularization parameter for AD with OSL-MLEM
options.beta_ad_mlem = 0.1;
%%% Regularization parameter for AD with MBSREM
options.beta_ad_mbsrem = 0.3;
%%% Regularization parameter for AD with BSREM
options.beta_ad_bsrem = 0.2;
%%% Regularization parameter for AD with ROSEM
options.beta_ad_rosem = 0.0003;
%%% Regularization parameter for AD with RBI
options.beta_ad_rbi = 0.05;
%%% Regularization parameter for AD with (OSL-(A)COSEM)
options.beta_ad_cosem = 0.1;
%%% Time step variable for AD (implementation 2 only)
options.TimeStepAD = 0.0625;
%%% Conductivity/connectivity for AD (edge threshold)
options.KAD = 2;
%%% Number of iterations for AD filter
% NOTE: This refers to the AD smoothing part, not the actual reconstruction
% phase
options.NiterAD = 1;
%%% Flux/conduction type for AD filter
% 1 = Exponential
% 2 = Quadratic
options.FluxType = 1;
%%% Diffusion type for AD (implementation 2 only)
% 1 = Gradient
% 2 = Modified curvature
options.DiffusionType = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% APLS PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for APLS with OSL-OSEM
options.beta_APLS_osem = 0.01;
%%% Regularization parameter for APLS with OSL-MLEM
options.beta_APLS_mlem = 0.1;
%%% Regularization parameter for APLS with MBSREM
options.beta_APLS_mbsrem = 0.1;
%%% Regularization parameter for APLS with BSREM
options.beta_APLS_bsrem = 0.005;
%%% Regularization parameter for APLS with ROSEM
options.beta_APLS_rosem = 0.1;
%%% Regularization parameter for APLS with RBI
options.beta_APLS_rbi = 0.1;
%%% Regularization parameter for APLS (OSL-(A)COSEM)
options.beta_APLS_cosem = 0.01;
%%% Scaling parameter (eta)
% See the wiki for details
options.eta = 1e-5;
%%% "Smoothing" parameter (beta)
% Also used to prevent zero values in square root
options.APLSsmoothing = 1e-5;
%%% Specify filename for the reference image here (same rules apply as with
% attenuation correction above)
% NOTE: For APSL, the prior is required
options.APLS_reference_image = 'reference_image.mat';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TGV PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for TGV with OSL-OSEM
options.beta_TGV_osem = 0.1;
%%% Regularization parameter for TGV with OSL-MLEM
options.beta_TGV_mlem = 0.1;
%%% Regularization parameter for TGV with MBSREM
options.beta_TGV_mbsrem = 0.1;
%%% Regularization parameter for TGV with BSREM
options.beta_TGV_bsrem = 1;
%%% Regularization parameter for TGV with ROSEM
options.beta_TGV_rosem = 0.25;
%%% Regularization parameter for TGV with RBI
options.beta_TGV_rbi = 0.1;
%%% Regularization parameter for TGV (OSL-(A)COSEM)
options.beta_TGV_cosem = 0.05;
%%% TGV weights
% First part
options.betaTGV = 1;
% Second part (symmetrized derivative)
options.alphaTGV = 2;
%%% Number of TGV iterations
options.NiterTGV = 30;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NLM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% EXPERIMENTAL FEATURE %%%%
%%% Regularization parameter for NLM with OSL-OSEM
options.beta_NLM_osem = 0.025;
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
options.Nlx = 1;
options.Nly = 1;
options.Nlz = 0;
% Search window radius is controlled by Ndx, Ndy and Ndz parameters
% Use anatomical reference image for the patches
options.NLM_use_anatomical = true;
%%% Specify filename for the reference image here (same rules apply as with
% attenuation correction above)
options.NLM_reference_image = 'reference_image.mat';
%%% Use MRP algorithm (without normalization)
% I.e. gradient = im - NLM_filtered(im)
options.NLM_MRP = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% OPENCL DEVICE INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Implementation 3
% Uncomment the below line and run it to determine the available platforms,
% their respective numbers and device numbers
% OpenCL_device_info();

%%% Implementation 2
% Uncomment the below line and run it to determine the available device
% numbers
% ArrayFire_OpenCL_device_info();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD SCATTER DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load user scatter data
% Load scatter data (if applicable)
if options.scatter_correction && ~options.only_reconstructions || options.scatter_correction && options.use_raw_data
    options = loadScatterData(options);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ERROR CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = set_GATE_variables(options);
% Basic error checking is done here
options = OMEGA_error_check(options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






if options.only_system_matrix == false
    % Load the measurement data
    options = loadMeasurementData(options);
else
    if options.use_raw_data == false
        options.SinM = false(options.Ndist, options.Nang, options.NSinos);
    else
        options.coincidences = {false(options.detectors ^2/2 + options.detectors/2)};
    end
end


%% Precompute the necessary data

if options.precompute && options.only_reconstructions == false
    precompute_data(options);
end

%% Reconstructions

if options.only_system_matrix == false
    
    tStart = tic;
    pz = reconstructions_main(options);
    tElapsed = toc(tStart);
    disp(['Reconstruction process took ' num2str(tElapsed) ' seconds'])
    
end

% save([options.name '_reconstruction_' num2str(options.subsets) 'subsets_' num2str(options.Niter) 'iterations.mat'], 'pz');

%% System matrix formation
% System matrix formation is ALWAYS computed with implementation 1,
% regardless of the choices made above

if options.only_system_matrix || options.single_reconstructions
    
    % This is necessary to produce the system matrix
    options = custom_prior_prepass(options, options.single_reconstructions);
    
    if options.single_reconstructions
        im_vectors = form_image_vectors(options, options.Nx * options.Ny * options.Nz);
    end
    
    %%
    for llo = 1 : partitions
        SinD = 0;
        if options.single_reconstructions
            if iscell(options.SinM)
                Sino = options.SinM{llo};
            else
                Sino = options.SinM;
            end
            
            Sino = Sino(:);
            
            if issparse(Sino)
                Sino = (full(Sino));
            end
            if options.randoms_correction
                if iscell(options.SinDelayed)
                    SinD = double(options.SinDelayed{llo}(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                else
                    SinD = double(options.SinDelayed(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                end
                if issparse(SinD)
                    SinD = (full(SinD));
                end
                SinD = SinD(:);
            end
        end
        for iter = 1 : options.Niter
            for osa_iter = 1 : options.subsets
                %%%%%%%%%%%%%%%% Compute the system matrix %%%%%%%%%%%%%%%%
                A = observation_matrix_formation(options, osa_iter);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Use A here (y = Ax)
                %%% Separate reconstruction algorithms
                if options.single_reconstructions
                    % Sensitivity image/normalization constant
                    if options.is_transposed
                        Summ = full(sum(A,2));
                    else
                        Summ = full(sum(A,1))';
                    end
                    Summ(Summ == 0) = options.epps;
                    uu = double(Sino(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1)));
                    % Compute OSEM
                    if options.osem || options.ecosem
                        if verbose
                            tStart = tic;
                        end
                        im_vectors.OSEM_apu = OSEM_im(im_vectors.OSEM_apu, A, options.epps, uu, Summ, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['OSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute MRAMLA
                    if options.mramla
                        if verbose
                            tStart = tic;
                        end
                        im_vectors.MRAMLA_apu = MBSREM(im_vectors.MRAMLA_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, SinD, randoms_correction, ...
                            options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['MRAMLA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute RAMLA
                    if options.ramla
                        if verbose
                            tStart = tic;
                        end
                        im_vectors.RAMLA_apu = BSREM_subiter(im_vectors.RAMLA_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed);
                        if any(im_vectors.RAMLA_apu < 0)
                            error('Negative values in RAMLA, lower lambda value!')
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['RAMLA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute ROSEM
                    if options.rosem
                        if verbose
                            tStart = tic;
                        end
                        im_vectors.ROSEM_apu = ROSEM_subiter(im_vectors.ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute RBI
                    if options.rbi
                        if verbose
                            tStart = tic;
                        end
                        im_vectors.RBI_apu = RBI_subiter(im_vectors.RBI_apu, A, uu, options.epps, Summ, 0, 0, options.D, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['RBI sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute DRAMA
                    if options.drama
                        if verbose
                            tStart = tic;
                        end
                        im_vectors.DRAMA_apu = DRAMA_subiter(im_vectors.DRAMA_apu, options.lam_drama, options.epps, iter, Summ, osa_iter, A, uu, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['DRAMA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute COSEM
                    if options.cosem || options.ecosem
                        if verbose
                            tStart = tic;
                        end
                        [im_vectors.COSEM_apu, options.C_co] = COSEM_im(im_vectors.COSEM_apu, A, options.epps, uu, options.C_co, options.D, osa_iter, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['COSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute ECOSEM
                    if options.ecosem
                        if verbose
                            tStart = tic;
                        end
                        im_vectors.ECOSEM_apu = ECOSEM_im(im_vectors.ECOSEM_apu, options.epps, options.D, im_vectors.COSEM_apu, im_vectors.OSEM_apu);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ECOSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute ACOSEM
                    if options.acosem
                        if verbose
                            tStart = tic;
                        end
                        [im_vectors.ACOSEM_apu, options.C_aco] = ACOSEM_im(im_vectors.ACOSEM_apu, A, options.epps, uu, options.C_aco, options.D, options.h, osa_iter, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ACOSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute OSL with MRP
                    if options.MRP && options.OSL_OSEM
                        if verbose
                            tStart = tic;
                        end
                        med = MRP(im_vectors.MRP_OSL_apu, options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
                        im_vectors.MRP_OSL_apu = OSL_OSEM(im_vectors.MRP_OSL_apu, Summ, options.beta_mrp_osem, med, options.epps, A, uu, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute MBSREM with MRP
                    if options.MRP && options.MBSREM
                        if verbose
                            tStart = tic;
                        end
                        med = MRP(im_vectors.MRP_MBSREM_apu, options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
                        im_vectors.MRP_MBSREM_apu = MBSREM(im_vectors.MRP_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, SinD, randoms_correction, options.is_transposed, ...
                            options.beta_mrp_mbsrem, med);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['MBSREM MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute BSREM with MRP
                    if options.MRP && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1 && options.ramla
                            im_vectors.MRP_BSREM_apu = im_vectors.RAMLA_apu;
                        else
                            im_vectors.MRP_BSREM_apu = BSREM_subiter(im_vectors.MRP_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed);
                        end
                        if any(im_vectors.MRP_BSREM_apu < 0)
                            error('Negative values in BSREM, lower lambda value!')
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute ROSEM-MAP with MRP
                    if options.MRP && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1 && options.rosem
                            im_vectors.MRP_ROSEM_apu = im_vectors.ROSEM_apu;
                        else
                            im_vectors.MRP_ROSEM_apu = ROSEM_subiter(im_vectors.MRP_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, options.is_transposed);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM-MAP MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute RBI-MAP with MRP
                    if options.MRP && options.RBI_MAP
                        if verbose
                            tStart = tic;
                        end
                        med = MRP(im_vectors.MRP_RBI_apu, options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
                        im_vectors.MRP_RBI_apu = RBI_subiter(im_vectors.MRP_RBI_apu, A, uu, options.epps, Summ, options.beta_mrp_rbi, med, options.D, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['RBI-MAP MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute COSEM-MAP with MRP
                    if options.MRP && any(options.COSEM_MAP)
                        if verbose
                            tStart = tic;
                        end
                        med = MRP(im_vectors.MRP_COSEM_apu, options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
                        if options.COSEM_MAP == 1
                            [im_vectors.MRP_COSEM_apu, options.C_osl] = COSEM_OSL(im_vectors.MRP_COSEM_apu, options.D, options.beta_mrp_cosem, med, options.epps, A, uu, ...
                                options.C_osl, options.h, options.COSEM_MAP, osa_iter, SinD, options.is_transposed);
                        else
                            [im_vectors.MRP_COSEM_apu, options.C_osl] = COSEM_OSL(im_vectors.MRP_COSEM_apu, options.D, options.beta_mrp_cosem, med, options.epps, A, uu, ...
                                options.C_osl, 0, options.COSEM_MAP, osa_iter, SinD, options.is_transposed);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['COSEM-MAP MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute OSL with Quadratic prior
                    if options.quad && options.OSL_OSEM
                        if verbose
                            tStart = tic;
                        end
                        med = Quadratic_prior(im_vectors.Quad_OSL_apu, options.tr_offsets, options.weights, options.weights_quad, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz);
                        im_vectors.Quad_OSL_apu = OSL_OSEM(im_vectors.Quad_OSL_apu, Summ, options.beta_quad_osem, med, options.epps, A, uu, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['OSL Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute MBSREM with Quadratic prior
                    if options.quad && options.MBSREM
                        if verbose
                            tStart = tic;
                        end
                        med = Quadratic_prior(im_vectors.Quad_MBSREM_apu, options.tr_offsets, options.weights, options.weights_quad, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz);
                        im_vectors.Quad_MBSREM_apu = MBSREM(im_vectors.Quad_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, SinD, randoms_correction, options.is_transposed, ...
                            options.beta_quad_mbsrem, med);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['MBSREM Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute BSREM with Quadratic prior
                    if options.quad && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1 && options.ramla
                            im_vectors.Quad_BSREM_apu = im_vectors.RAMLA_apu;
                        else
                            im_vectors.Quad_BSREM_apu = BSREM_subiter(im_vectors.Quad_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed);
                        end
                        if any(im_vectors.Quad_BSREM_apu < 0)
                            error('Negative values in BSREM, lower lambda value')
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute ROSEM-MAP with Quadratic prior
                    if options.quad && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1 && options.rosem
                            im_vectors.Quad_ROSEM_apu = im_vectors.ROSEM_apu;
                        else
                            im_vectors.Quad_ROSEM_apu = ROSEM_subiter(im_vectors.Quad_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, options.is_transposed);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM-MAP Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute RBI-MAP with Quadratic prior
                    if options.quad && options.RBI_MAP
                        if verbose
                            tStart = tic;
                        end
                        med = Quadratic_prior(im_vectors.Quad_RBI_apu, options.tr_offsets, options.weights, options.weights_quad, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz);
                        im_vectors.Quad_RBI_apu = RBI_subiter(im_vectors.Quad_RBI_apu, A, uu, options.epps, Summ, options.beta_quad_rbi, med, options.D, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['RBI-MAP Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute COSEM-MAP with Quadratic prior
                    if options.quad && any(options.COSEM_MAP)
                        if verbose
                            tStart = tic;
                        end
                        med = Quadratic_prior(im_vectors.Quad_COSEM_apu, options.tr_offsets, options.weights, options.weights_quad, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz);
                        if options.COSEM_MAP == 1
                            [im_vectors.Quad_COSEM_apu, options.C_osl] = COSEM_OSL(im_vectors.Quad_COSEM_apu, options.D, options.beta_quad_cosem, med, options.epps, A, uu, ...
                                options.C_osl, options.h, options.COSEM_MAP, osa_iter, SinD, options.is_transposed);
                        else
                            [im_vectors.Quad_COSEM_apu, options.C_osl] = COSEM_OSL(im_vectors.Quad_COSEM_apu, options.D, options.beta_quad_cosem, med, options.epps, A, uu, ...
                                options.C_osl, 0, options.COSEM_MAP, osa_iter, SinD, options.is_transposed);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['COSEM-MAP Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute OSL with L-filter prior
                    if options.L && options.OSL_OSEM
                        if verbose
                            tStart = tic;
                        end
                        med = L_filter(im_vectors.L_OSL_apu, options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, options.epps, options.med_no_norm);
                        im_vectors.L_OSL_apu = OSL_OSEM(im_vectors.L_OSL_apu, Summ, options.beta_L_osem, med, options.epps, A, uu, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['OSL L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute MBSREM with L-filter prior
                    if options.L && options.MBSREM
                        if verbose
                            tStart = tic;
                        end
                        med = L_filter(im_vectors.L_MBSREM_apu, options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, options.epps, options.med_no_norm);
                        im_vectors.L_MBSREM_apu = MBSREM(im_vectors.L_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, SinD, randoms_correction, options.is_transposed, ...
                            options.beta_L_mbsrem, med);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['MBSREM L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute BSREM with L-filter prior
                    if options.L && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1 && options.ramla
                            im_vectors.L_BSREM_apu = im_vectors.RAMLA_apu;
                        else
                            im_vectors.L_BSREM_apu = BSREM_subiter(im_vectors.L_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed);
                        end
                        if any(im_vectors.L_BSREM_apu < 0)
                            error('Negative values in BSREM, lower lambda value')
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute ROSEM-MAP with L-filter prior
                    if options.L && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1 && options.rosem
                            im_vectors.L_ROSEM_apu = im_vectors.ROSEM_apu;
                        else
                            im_vectors.L_ROSEM_apu = ROSEM_subiter(im_vectors.L_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, options.is_transposed);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM-MAP L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute RBI-MAP with L-filter prior
                    if options.L && options.RBI_MAP
                        if verbose
                            tStart = tic;
                        end
                        med = L_filter(im_vectors.L_RBI_apu, options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, options.epps, options.med_no_norm);
                        im_vectors.L_RBI_apu = RBI_subiter(im_vectors.L_RBI_apu, A, uu, options.epps, Summ, options.beta_L_rbi, med, options.D, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['RBI-MAP L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute COSEM-MAP with L-filter prior
                    if options.L && any(options.COSEM_MAP)
                        if verbose
                            tStart = tic;
                        end
                        med = L_filter(im_vectors.L_COSEM_apu, options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, options.epps, options.med_no_norm);
                        if options.COSEM_MAP == 1
                            [im_vectors.L_COSEM_apu, options.C_osl] = COSEM_OSL(im_vectors.L_COSEM_apu, options.D, options.beta_L_cosem, med, options.epps, A, uu, options.C_osl, ...
                                options.h, options.COSEM_MAP, osa_iter, SinD, options.is_transposed);
                        else
                            [im_vectors.L_COSEM_apu, options.C_osl] = COSEM_OSL(im_vectors.L_COSEM_apu, options.D, options.beta_L_cosem, med, options.epps, A, uu, options.C_osl, 0, ...
                                options.COSEM_MAP, osa_iter, SinD, options.is_transposed);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['COSEM-MAP L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute OSL with FMH prior
                    if options.FMH && options.OSL_OSEM
                        if verbose
                            tStart = tic;
                        end
                        med = FMH(im_vectors.FMH_OSL_apu, options.tr_offsets, options.fmh_weights, options.weights, options.Nx, options.Ny, options.Nz, N, Ndx, Ndy, Ndz, options.epps, ...
                            options.med_no_norm);
                        im_vectors.FMH_OSL_apu = OSL_OSEM(im_vectors.FMH_OSL_apu, Summ, options.beta_fmh_osem, med, options.epps, A, uu, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['OSL FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.FMH && options.MBSREM
                        if verbose
                            tStart = tic;
                        end
                        med = FMH(im_vectors.FMH_MBSREM_apu, options.tr_offsets, options.fmh_weights, options.weights, options.Nx, options.Ny, options.Nz, N, Ndx, Ndy, Ndz, options.epps, ...
                            options.med_no_norm);
                        im_vectors.FMH_MBSREM_apu = MBSREM(im_vectors.FMH_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, SinD, randoms_correction, options.is_transposed, ...
                            options.beta_fmh_mbsrem, med);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['MBSREM FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.FMH && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1 && options.ramla
                            im_vectors.FMH_BSREM_apu = im_vectors.RAMLA_apu;
                        else
                            im_vectors.FMH_BSREM_apu = BSREM_subiter(im_vectors.FMH_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed);
                        end
                        if any(im_vectors.FMH_BSREM_apu < 0)
                            error('Negative values in BSREM, lower lambda value')
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.FMH && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1 && options.rosem
                            im_vectors.FMH_ROSEM_apu = im_vectors.ROSEM_apu;
                        else
                            im_vectors.FMH_ROSEM_apu = ROSEM_subiter(im_vectors.FMH_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, options.is_transposed);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM-MAP FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.FMH && options.RBI_MAP
                        if verbose
                            tStart = tic;
                        end
                        med = FMH(im_vectors.FMH_RBI_apu, options.tr_offsets, options.fmh_weights, options.weights, options.Nx, options.Ny, options.Nz, N, Ndx, Ndy, Ndz, options.epps, ...
                            options.med_no_norm);
                        im_vectors.FMH_RBI_apu = RBI_subiter(im_vectors.FMH_RBI_apu, A, uu, options.epps, Summ, options.beta_fmh_rbi, med, options.D, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['RBI-MAP FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.FMH && any(options.COSEM_MAP)
                        if verbose
                            tStart = tic;
                        end
                        med = FMH(im_vectors.FMH_COSEM_apu, options.tr_offsets, options.fmh_weights, options.weights, options.Nx, options.Ny, options.Nz, N, Ndx, Ndy, Ndz, options.epps, ...
                            options.med_no_norm);
                        if options.COSEM_MAP == 1
                            [im_vectors.FMH_COSEM_apu, options.C_osl] = COSEM_OSL(im_vectors.FMH_COSEM_apu, options.D, options.beta_fmh_cosem, med, options.epps, A, uu, ...
                                options.C_osl, options.h, options.COSEM_MAP, osa_iter, SinD, options.is_transposed);
                        else
                            [im_vectors.FMH_COSEM_apu, options.C_osl] = COSEM_OSL(im_vectors.FMH_COSEM_apu, options.D, options.beta_fmh_cosem, med, options.epps, A, uu, ...
                                options.C_osl, 0, options.COSEM_MAP, osa_iter, SinD, options.is_transposed);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['COSEM-MAP FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute OSL with weighted mean prior
                    if options.weighted_mean && options.OSL_OSEM
                        if verbose
                            tStart = tic;
                        end
                        med = Weighted_mean(im_vectors.Weighted_OSL_apu, options.tr_offsets, options.weighted_weights, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, ...
                            options.mean_type, options.epps, options.w_sum, options.med_no_norm);
                        im_vectors.Weighted_OSL_apu = OSL_OSEM(im_vectors.Weighted_OSL_apu, Summ, options.beta_weighted_osem, med, options.epps, A, uu, SinD, ...
                            options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['OSL weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.weighted_mean && options.MBSREM
                        if verbose
                            tStart = tic;
                        end
                        med = Weighted_mean(im_vectors.Weighted_MBSREM_apu, options.tr_offsets, options.weighted_weights, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, ...
                            options.mean_type, options.epps, options.w_sum, options.med_no_norm);
                        im_vectors.Weighted_MBSREM_apu = MBSREM(im_vectors.Weighted_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, ...
                            options.lam_mbsrem, iter, SinD, randoms_correction, options.is_transposed, options.beta_weighted_mbsrem, med);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['MBSREM weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.weighted_mean && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1 && options.ramla
                            im_vectors.Weighted_BSREM_apu = im_vectors.RAMLA_apu;
                        else
                            im_vectors.Weighted_BSREM_apu = BSREM_subiter(im_vectors.Weighted_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed);
                        end
                        if any(im_vectors.Weighted_BSREM_apu < 0)
                            error('Negative values in BSREM, lower lambda value')
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.weighted_mean && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1 && options.rosem
                            im_vectors.Weighted_ROSEM_apu = im_vectors.ROSEM_apu;
                        else
                            im_vectors.Weighted_ROSEM_apu = ROSEM_subiter(im_vectors.Weighted_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, ...
                                options.is_transposed);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM-MAP weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.weighted_mean && options.RBI_MAP
                        if verbose
                            tStart = tic;
                        end
                        med = Weighted_mean(im_vectors.Weighted_RBI_apu, options.tr_offsets, options.weighted_weights, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, ...
                            options.mean_type, options.epps, options.w_sum, options.med_no_norm);
                        im_vectors.Weighted_RBI_apu = RBI_subiter(im_vectors.Weighted_RBI_apu, A, uu, options.epps, Summ, options.beta_weighted_rbi, ...
                            med, options.D, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['RBI-MAP weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.weighted_mean && any(options.COSEM_MAP)
                        if verbose
                            tStart = tic;
                        end
                        med = Weighted_mean(im_vectors.Weighted_COSEM_apu, options.tr_offsets, options.weighted_weights, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, ...
                            options.mean_type, options.epps, options.w_sum, options.med_no_norm);
                        if options.COSEM_MAP == 1
                            [im_vectors.Weighted_COSEM_apu, options.C_osl] = COSEM_OSL(im_vectors.Weighted_COSEM_apu, options.D, options.beta_weighted_cosem, ...
                                med, options.epps, A, uu, options.C_osl, options.h, options.COSEM_MAP, osa_iter, SinD, options.is_transposed);
                        else
                            [im_vectors.Weighted_COSEM_apu, options.C_osl] = COSEM_OSL(im_vectors.Weighted_COSEM_apu, options.D, options.beta_weighted_cosem, ...
                                med, options.epps, A, uu, options.C_osl, 0, options.COSEM_MAP, osa_iter, SinD, options.is_transposed);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['COSEM-MAP weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute OSL with TV prior
                    if options.TV && options.OSL_OSEM
                        if verbose
                            tStart = tic;
                        end
                        grad = TVpriorFinal(im_vectors.TV_OSL_apu, options.TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, ...
                            options.tr_offsets);
                        im_vectors.TV_OSL_apu = OSL_OSEM(im_vectors.TV_OSL_apu, Summ, options.beta_TV_osem, grad, options.epps, A, uu, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['OSL TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.TV && options.MBSREM
                        if verbose
                            tStart = tic;
                        end
                        grad = TVpriorFinal(im_vectors.TV_MBSREM_apu, options.TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, ...
                            options.tr_offsets);
                        im_vectors.TV_MBSREM_apu = MBSREM(im_vectors.TV_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                            iter, SinD, randoms_correction, options.is_transposed, options.beta_TV_mbsrem, grad);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['MBSREM TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.TV && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1 && options.ramla
                            im_vectors.TV_BSREM_apu = im_vectors.RAMLA_apu;
                        else
                            im_vectors.TV_BSREM_apu = BSREM_subiter(im_vectors.TV_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed);
                        end
                        if any(im_vectors.TV_BSREM_apu < 0)
                            error('Negative values in BSREM, lower lambda value')
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.TV && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1 && options.rosem
                            im_vectors.TV_ROSEM_apu = im_vectors.ROSEM_apu;
                        else
                            im_vectors.TV_ROSEM_apu = ROSEM_subiter(im_vectors.TV_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, options.is_transposed);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM-MAP TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.TV && options.RBI_MAP
                        if verbose
                            tStart = tic;
                        end
                        grad = TVpriorFinal(im_vectors.TV_RBI_apu, options.TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, ...
                            options.tr_offsets);
                        im_vectors.TV_RBI_apu = RBI_subiter(im_vectors.TV_RBI_apu, A, uu, options.epps, Summ, options.beta_TV_rbi, grad, options.D, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['RBI-MAP TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.TV && any(options.COSEM_MAP)
                        if verbose
                            tStart = tic;
                        end
                        grad = TVpriorFinal(im_vectors.TV_COSEM_apu, options.TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, ...
                            options.tr_offsets);
                        if options.COSEM_MAP == 1
                            [im_vectors.TV_COSEM_apu, options.C_osl] = COSEM_OSL(im_vectors.TV_COSEM_apu, options.D, options.beta_TV_cosem, grad, options.epps, A, uu, ...
                                options.C_osl, options.h, options.COSEM_MAP, osa_iter, SinD, options.is_transposed);
                        else
                            [im_vectors.TV_COSEM_apu, options.C_osl] = COSEM_OSL(im_vectors.TV_COSEM_apu, options.D, options.beta_TV_cosem, grad, options.epps, A, uu, ...
                                options.C_osl, 0, options.COSEM_MAP, osa_iter, SinD, options.is_transposed);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['COSEM-MAP TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute OSL with MRP-AD prior
                    if options.AD && options.OSL_OSEM
                        if verbose
                            tStart = tic;
                        end
                        if osa_iter > 1
                            med = AD(im_vectors.AD_OSL_apu, options.FluxType, options.Nx, options.Ny, options.Nz, options);
                            im_vectors.AD_OSL_apu = OSL_OSEM(im_vectors.AD_OSL_apu, Summ, options.beta_ad_osem, med, options.epps, A, uu, SinD, options.is_transposed);
                        else
                            im_vectors.AD_OSL_apu = OSEM_im(im_vectors.AD_OSL_apu, A, options.epps, uu, Summ);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['OSL AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.AD && options.MBSREM
                        if verbose
                            tStart = tic;
                        end
                        med = AD(im_vectors.AD_MBSREM_apu, options.FluxType, options.Nx, options.Ny, options.Nz, options);
                        im_vectors.AD_MBSREM_apu = MBSREM(im_vectors.AD_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                            iter, SinD, randoms_correction, options.is_transposed, options.beta_ad_mbsrem, med);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['MBSREM AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.AD && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1 && options.ramla
                            im_vectors.AD_BSREM_apu = im_vectors.RAMLA_apu;
                        else
                            im_vectors.AD_BSREM_apu = BSREM_subiter(im_vectors.AD_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed);
                        end
                        if any(im_vectors.AD_BSREM_apu < 0)
                            error('Negative values in BSREM, lower lambda value')
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.AD && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1 && options.rosem
                            im_vectors.AD_ROSEM_apu = im_vectors.ROSEM_apu;
                        else
                            im_vectors.AD_ROSEM_apu = ROSEM_subiter(im_vectors.AD_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, options.is_transposed);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM-MAP AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.AD && options.RBI_MAP
                        if verbose
                            tStart = tic;
                        end
                        med = AD(im_vectors.AD_RBI_apu, options.FluxType, options.Nx, options.Ny, options.Nz, options);
                        im_vectors.AD_RBI_apu = RBI_subiter(im_vectors.AD_RBI_apu, A, uu, options.epps, Summ, options.beta_ad_rbi, med, options.D, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['RBI-MAP AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.AD && any(options.COSEM_MAP)
                        if verbose
                            tStart = tic;
                        end
                        med = AD(im_vectors.AD_COSEM_apu, options.FluxType, options.Nx, options.Ny, options.Nz, options);
                        if options.COSEM_MAP == 1
                            [im_vectors.AD_COSEM_apu, options.C_osl] = COSEM_OSL(im_vectors.AD_COSEM_apu, options.D, options.beta_ad_cosem, med, options.epps, A, uu, ...
                                options.C_osl, options.h, options.COSEM_MAP, osa_iter, SinD, options.is_transposed);
                        else
                            [im_vectors.AD_COSEM_apu, options.C_osl] = COSEM_OSL(im_vectors.AD_COSEM_apu, options.D, options.beta_ad_cosem, med, options.epps, A, uu, ...
                                options.C_osl, 0, options.COSEM_MAP, osa_iter, SinD, options.is_transposed);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['COSEM-MAP AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute OSL with APLS prior
                    if options.APLS && options.OSL_OSEM
                        if verbose
                            tStart = tic;
                        end
                        grad = TVpriorFinal(im_vectors.APLS_OSL_apu, [], options.Nx, options.Ny, options.Nz, true, options, 4);
                        im_vectors.APLS_OSL_apu = OSL_OSEM(im_vectors.APLS_OSL_apu, Summ, options.beta_APLS_osem, grad, options.epps, A, uu, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['OSL APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.APLS && options.MBSREM
                        if verbose
                            tStart = tic;
                        end
                        grad = TVpriorFinal(im_vectors.APLS_MBSREM_apu, [], options.Nx, options.Ny, options.Nz, true, options, 4);
                        im_vectors.APLS_MBSREM_apu = MBSREM(im_vectors.APLS_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                            iter, SinD, randoms_correction, options.is_transposed, options.beta_APLS_mbsrem, grad);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['MBSREM APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.APLS && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1 && options.ramla
                            im_vectors.APLS_BSREM_apu = im_vectors.RAMLA_apu;
                        else
                            im_vectors.APLS_BSREM_apu = BSREM_subiter(im_vectors.APLS_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed);
                        end
                        if any(im_vectors.APLS_BSREM_apu < 0)
                            error('Negative values in BSREM, lower lambda value')
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.APLS && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1 && options.rosem
                            im_vectors.APLS_ROSEM_apu = im_vectors.ROSEM_apu;
                        else
                            im_vectors.APLS_ROSEM_apu = ROSEM_subiter(im_vectors.APLS_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, options.is_transposed);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM-MAP APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.APLS && options.RBI_MAP
                        if verbose
                            tStart = tic;
                        end
                        grad = TVpriorFinal(im_vectors.APLS_RBI_apu, [], options.Nx, options.Ny, options.Nz, true, options, 4);
                        im_vectors.APLS_RBI_apu = RBI_subiter(im_vectors.APLS_RBI_apu, A, uu, options.epps, Summ, SinD, options.beta_APLS_rbi, grad, options.D, ...
                            options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['RBI-MAP APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.APLS && any(options.COSEM_MAP)
                        if verbose
                            tStart = tic;
                        end
                        grad = TVpriorFinal(im_vectors.APLS_COSEM_apu, [], options.Nx, options.Ny, options.Nz, true, options, 4);
                        if options.COSEM_MAP == 1
                            [im_vectors.APLS_COSEM_apu, options.C_osl] = COSEM_OSL(im_vectors.APLS_COSEM_apu, options.D, options.beta_APLS_cosem, grad, A, uu, ...
                                options.epps, options.C_osl, options.h, options.COSEM_MAP, osa_iter, SinD, options.is_transposed);
                        else
                            [im_vectors.APLS_COSEM_apu, options.C_osl] = COSEM_OSL(im_vectors.APLS_COSEM_apu, options.D, options.beta_APLS_cosem, grad, A, uu, ...
                                options.epps, options.C_osl, 0, options.COSEM_MAP, osa_iter, SinD, options.is_transposed);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['COSEM-MAP APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute OSL with TGV prior
                    if options.TGV && options.OSL_OSEM
                        if verbose
                            tStart = tic;
                        end
                        grad = TGV(im_vectors.TGV_OSL_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
                        im_vectors.TGV_OSL_apu = OSL_OSEM(im_vectors.TGV_OSL_apu, Summ, options.beta_TGV_osem, grad, options.epps, A, uu, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['OSL TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.TGV && options.MBSREM
                        if verbose
                            tStart = tic;
                        end
                        grad = TGV(im_vectors.TGV_MBSREM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
                        im_vectors.TGV_MBSREM_apu = MBSREM(im_vectors.TGV_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                            iter, SinD, randoms_correction, options.is_transposed, options.beta_TGV_mbsrem, grad);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['MBSREM TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.TGV && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1 && options.ramla
                            im_vectors.TGV_BSREM_apu = im_vectors.RAMLA_apu;
                        else
                            im_vectors.TGV_BSREM_apu = BSREM_subiter(im_vectors.TGV_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed);
                        end
                        if any(im_vectors.TGV_BSREM_apu < 0)
                            error('Negative values in BSREM, lower lambda value')
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.TGV && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1 && options.rosem
                            im_vectors.TGV_ROSEM_apu = im_vectors.ROSEM_apu;
                        else
                            im_vectors.TGV_ROSEM_apu = ROSEM_subiter(im_vectors.TGV_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, options.is_transposed);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM-MAP TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.TGV && options.RBI_MAP
                        if verbose
                            tStart = tic;
                        end
                        grad = TGV(im_vectors.TGV_RBI_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
                        im_vectors.TGV_RBI_apu = RBI_subiter(im_vectors.TGV_RBI_apu, A, uu, options.epps, Summ, options.beta_TGV_rbi, grad, options.D, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['RBI-MAP TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.TGV && any(options.COSEM_MAP)
                        if verbose
                            tStart = tic;
                        end
                        grad = TGV(im_vectors.TGV_COSEM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
                        if options.COSEM_MAP == 1
                            [im_vectors.TGV_COSEM_apu, options.C_osl] = COSEM_OSL(im_vectors.TGV_COSEM_apu, options.D, options.beta_TGV_cosem, grad, A, uu, ...
                                options.epps, options.C_osl, options.h, options.COSEM_MAP, osa_iter, SinD, options.is_transposed);
                        else
                            [im_vectors.TGV_COSEM_apu, options.C_osl] = COSEM_OSL(im_vectors.TGV_COSEM_apu, options.D, options.beta_TGV_cosem, grad, A, uu, ...
                                options.epps, options.C_osl, 0, options.COSEM_MAP, osa_iter, SinD, options.is_transposed);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['COSEM-MAP TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    % Compute OSL with NLM prior
                    if options.NLM && options.OSL_OSEM
                        if verbose
                            tStart = tic;
                        end
                        med = NLM(im_vectors.NLM_OSL_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                            options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
                        im_vectors.NLM_OSL_apu = OSL_OSEM(im_vectors.NLM_OSL_apu, Summ, options.beta_NLM_osem, med, options.epps, A, uu, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['OSL NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.NLM && options.MBSREM
                        if verbose
                            tStart = tic;
                        end
                        med = NLM(im_vectors.NLM_MBSREM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                            options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
                        im_vectors.NLM_MBSREM_apu = MBSREM(im_vectors.NLM_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                            iter, SinD, randoms_correction, options.is_transposed, options.beta_NLM_mbsrem, med);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['MBSREM NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.NLM && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1 && options.ramla
                            im_vectors.NLM_BSREM_apu = im_vectors.RAMLA_apu;
                        else
                            im_vectors.NLM_BSREM_apu = BSREM_subiter(im_vectors.NLM_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed);
                        end
                        if any(im_vectors.NLM_BSREM_apu < 0)
                            error('Negative values in BSREM, lower lambda value')
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.NLM && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1 && options.rosem
                            im_vectors.NLM_ROSEM_apu = im_vectors.ROSEM_apu;
                        else
                            im_vectors.NLM_ROSEM_apu = ROSEM_subiter(im_vectors.NLM_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, options.is_transposed);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM-MAP NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.NLM && options.RBI_MAP
                        if verbose
                            tStart = tic;
                        end
                        med = NLM(im_vectors.NLM_RBI_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                            options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
                        im_vectors.NLM_RBI_apu = RBI_subiter(im_vectors.NLM_RBI_apu, A, uu, options.epps, Summ, options.beta_NLM_rbi, med, options.D, SinD, options.is_transposed);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['RBI-MAP NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                    if options.NLM && any(options.COSEM_MAP)
                        if verbose
                            tStart = tic;
                        end
                        med = NLM(im_vectors.NLM_COSEM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                            options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
                        if options.COSEM_MAP == 1
                            [im_vectors.NLM_COSEM_apu, options.C_osl] = COSEM_OSL(im_vectors.NLM_COSEM_apu, options.D, options.beta_NLM_cosem, med, A, uu, ...
                                options.epps, options.C_osl, options.h, options.COSEM_MAP, osa_iter, SinD, options.is_transposed);
                        else
                            [im_vectors.NLM_COSEM_apu, options.C_osl] = COSEM_OSL(im_vectors.NLM_COSEM_apu, options.D, options.beta_NLM_cosem, med, A, uu, ...
                                options.epps, options.C_osl, 0, options.COSEM_MAP, osa_iter, SinD, options.is_transposed);
                        end
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['COSEM-MAP NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                    end
                end
            end
            if options.MRP && options.BSREM
                if verbose
                    tStart = tic;
                end
                med = MRP(im_vectors.MRP_BSREM_apu, options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
                im_vectors.MRP_BSREM_apu = BSREM_iter(im_vectors.MRP_BSREM_apu, options.lam, iter, options.beta_mrp_bsrem, med, options.epps);
                if verbose
                    tElapsed = toc(tStart);
                    disp(['BSREM MRP iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                else
                    disp(['BSREM MRP iteration ' num2str(iter) ' finished'])
                end
            end
            
            if options.MRP && options.ROSEM_MAP
                if verbose
                    tStart = tic;
                end
                med = MRP(im_vectors.MRP_ROSEM_apu, options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
                im_vectors.MRP_ROSEM_apu = BSREM_iter(im_vectors.MRP_ROSEM_apu, options.lam_rosem, iter, options.beta_mrp_rosem, med, options.epps);
                if verbose
                    tElapsed = toc(tStart);
                    disp(['ROSEM MRP iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                else
                    disp(['ROSEM MRP iteration ' num2str(iter) ' finished'])
                end
            end
            
            if options.quad && options.BSREM
                if verbose
                    tStart = tic;
                end
                med = Quadratic_prior(im_vectors.Quad_BSREM_apu, options.tr_offsets, options.weights, options.weights_quad, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz);
                im_vectors.Quad_BSREM_apu = BSREM_iter(im_vectors.Quad_BSREM_apu, options.lam, iter, options.beta_quad_bsrem, med, options.epps);
                if verbose
                    tElapsed = toc(tStart);
                    disp(['BSREM quadratic iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                else
                    disp(['BSREM quadratic iteration ' num2str(iter) ' finished'])
                end
            end
            
            if options.quad && options.ROSEM_MAP
                if verbose
                    tStart = tic;
                end
                med = Quadratic_prior(im_vectors.Quad_ROSEM_apu, options.tr_offsets, options.weights, options.weights_quad, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz);
                im_vectors.Quad_ROSEM_apu = BSREM_iter(im_vectors.Quad_ROSEM_apu, options.lam_rosem, iter, options.beta_quad_rosem, med, options.epps);
                if verbose
                    tElapsed = toc(tStart);
                    disp(['ROSEM quadratic iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                else
                    disp(['ROSEM quadratic iteration ' num2str(iter) ' finished'])
                end
            end
            
            if options.L && options.BSREM
                if verbose
                    tStart = tic;
                end
                med = L_filter(im_vectors.L_BSREM_apu, options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, options.epps, options.med_no_norm);
                im_vectors.L_BSREM_apu = BSREM_iter(im_vectors.L_BSREM_apu, options.lam, iter, options.beta_L_bsrem, med, options.epps);
                if verbose
                    tElapsed = toc(tStart);
                    disp(['BSREM L-filter iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                else
                    disp(['BSREM L-filter iteration ' num2str(iter) ' finished'])
                end
            end
            
            if options.L && options.ROSEM_MAP
                if verbose
                    tStart = tic;
                end
                med = L_filter(im_vectors.L_ROSEM_apu, options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, options.epps, options.med_no_norm);
                im_vectors.L_ROSEM_apu = BSREM_iter(im_vectors.L_ROSEM_apu, options.lam_rosem, iter, options.beta_L_rosem, med, options.epps);
                if verbose
                    tElapsed = toc(tStart);
                    disp(['ROSEM L-filter iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                else
                    disp(['ROSEM L-filter iteration ' num2str(iter) ' finished'])
                end
            end
            
            if options.FMH && options.BSREM
                if verbose
                    tStart = tic;
                end
                med = FMH(im_vectors.FMH_BSREM_apu, options.tr_offsets, options.fmh_weights, options.weights, options.Nx, options.Ny, options.Nz, N, Ndx, Ndy, Ndz, options.epps, ...
                    options.med_no_norm);
                im_vectors.FMH_BSREM_apu = BSREM_iter(im_vectors.FMH_BSREM_apu, options.lam, iter, options.beta_fmh_bsrem, med, options.epps);
                if verbose
                    tElapsed = toc(tStart);
                    disp(['BSREM FMH iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                else
                    disp(['BSREM FMH iteration ' num2str(iter) ' finished'])
                end
            end
            
            if options.FMH && options.ROSEM_MAP
                if verbose
                    tStart = tic;
                end
                med = FMH(im_vectors.FMH_ROSEM_apu, options.tr_offsets, options.fmh_weights, options.weights, options.Nx, options.Ny, options.Nz, N, Ndx, Ndy, Ndz, options.epps, ...
                    options.med_no_norm);
                im_vectors.FMH_ROSEM_apu = BSREM_iter(im_vectors.FMH_ROSEM_apu, options.lam_rosem, iter, options.beta_fmh_rosem, med, options.epps);
                if verbose
                    tElapsed = toc(tStart);
                    disp(['ROSEM FMH iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                else
                    disp(['ROSEM FMH iteration ' num2str(iter) ' finished'])
                end
            end
            
            if options.weighted_mean && options.BSREM
                if verbose
                    tStart = tic;
                end
                med = Weighted_mean(im_vectors.Weighted_BSREM_apu, options.tr_offsets, options.weighted_weights, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, ...
                    options.mean_type, options.epps, options.w_sum, options.med_no_norm);
                im_vectors.Weighted_BSREM_apu = BSREM_iter(im_vectors.Weighted_BSREM_apu, options.lam, iter, options.beta_weighted_bsrem, ...
                    med, options.epps);
                if verbose
                    tElapsed = toc(tStart);
                    disp(['BSREM weighted mean iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                else
                    disp(['BSREM weighted mean iteration ' num2str(iter) ' finished'])
                end
            end
            
            if options.weighted_mean && options.ROSEM_MAP
                if verbose
                    tStart = tic;
                end
                med = Weighted_mean(im_vectors.Weighted_ROSEM_apu, options.tr_offsets, options.weighted_weights, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, ...
                    options.mean_type, options.epps, options.w_sum, options.med_no_norm);
                im_vectors.Weighted_ROSEM_apu = BSREM_iter(im_vectors.Weighted_ROSEM_apu, options.lam_rosem, iter, options.beta_weighted_rosem, ...
                    med, options.epps);
                if verbose
                    tElapsed = toc(tStart);
                    disp(['ROSEM weighted mean iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                else
                    disp(['ROSEM weighted mean iteration ' num2str(iter) ' finished'])
                end
            end
            
            if options.TV && options.BSREM
                if verbose
                    tStart = tic;
                end
                grad = TVpriorFinal(im_vectors.TV_BSREM_apu, options.TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, options.tr_offsets);
                im_vectors.TV_BSREM_apu = BSREM_iter(im_vectors.TV_BSREM_apu, options.lam, iter, options.beta_TV_bsrem, grad, options.epps);
                if verbose
                    tElapsed = toc(tStart);
                    disp(['BSREM TV iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                else
                    disp(['BSREM TV iteration ' num2str(iter) ' finished'])
                end
            end
            
            if options.TV && options.ROSEM_MAP
                if verbose
                    tStart = tic;
                end
                grad = TVpriorFinal(im_vectors.TV_ROSEM_apu, options.TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, options.tr_offsets);
                im_vectors.TV_ROSEM_apu = BSREM_iter(im_vectors.TV_ROSEM_apu, options.lam_rosem, iter, options.beta_TV_rosem, grad, options.epps);
                if verbose
                    tElapsed = toc(tStart);
                    disp(['ROSEM TV iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                else
                    disp(['ROSEM TV iteration ' num2str(iter) ' finished'])
                end
            end
            
            if options.AD && options.BSREM
                if verbose
                    tStart = tic;
                end
                med = AD(im_vectors.AD_BSREM_apu, options.FluxType, options.Nx, options.Ny, options.Nz, options);
                im_vectors.AD_BSREM_apu = BSREM_iter(im_vectors.AD_BSREM_apu, options.lam, iter, options.beta_ad_bsrem, med, options.epps);
                if verbose
                    tElapsed = toc(tStart);
                    disp(['BSREM AD iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                else
                    disp(['BSREM AD iteration ' num2str(iter) ' finished'])
                end
            end
            if options.AD && options.ROSEM_MAP
                if verbose
                    tStart = tic;
                end
                med = AD(im_vectors.AD_ROSEM_apu, options.FluxType, options.Nx, options.Ny, options.Nz, options);
                im_vectors.AD_ROSEM_apu = BSREM_iter(im_vectors.AD_ROSEM_apu, options.lam_rosem, iter, options.beta_ad_rosem, med, options.epps);
                if verbose
                    tElapsed = toc(tStart);
                    disp(['ROSEM AD iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                else
                    disp(['ROSEM AD iteration ' num2str(iter) ' finished'])
                end
            end
            
            if options.APLS && options.BSREM
                if verbose
                    tStart = tic;
                end
                grad = TVpriorFinal(im_vectors.APLS_BSREM_apu, 0, options.Nx, options.Ny, options.Nz, true, options, 4);
                im_vectors.APLS_BSREM_apu = BSREM_iter(im_vectors.APLS_BSREM_apu, options.lam, iter, options.beta_APLS_bsrem, grad, options.epps);
                if verbose
                    tElapsed = toc(tStart);
                    disp(['BSREM APLS iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                else
                    disp(['BSREM APLS iteration ' num2str(iter) ' finished'])
                end
            end
            if options.APLS && options.ROSEM_MAP
                if verbose
                    tStart = tic;
                end
                grad = TVpriorFinal(im_vectors.APLS_ROSEM_apu, 0, options.Nx, options.Ny, options.Nz, true, options, 4);
                im_vectors.APLS_ROSEM_apu = BSREM_iter(im_vectors.APLS_ROSEM_apu, options.lam_rosem, iter, options.beta_APLS_rosem, grad, ...
                    options.epps);
                if verbose
                    tElapsed = toc(tStart);
                    disp(['ROSEM APLS iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                else
                    disp(['ROSEM APLS iteration ' num2str(iter) ' finished'])
                end
            end
            
            if options.TGV && options.BSREM
                if verbose
                    tStart = tic;
                end
                grad = TGV(im_vectors.TGV_BSREM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
                im_vectors.TGV_BSREM_apu = BSREM_iter(im_vectors.TGV_BSREM_apu, options.lam, iter, options.beta_TGV_bsrem, grad, options.epps);
                if verbose
                    tElapsed = toc(tStart);
                    disp(['BSREM TGV iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                else
                    disp(['BSREM TGV iteration ' num2str(iter) ' finished'])
                end
            end
            if options.TGV && options.ROSEM_MAP
                if verbose
                    tStart = tic;
                end
                grad = TGV(im_vectors.TGV_ROSEM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
                im_vectors.TGV_ROSEM_apu = BSREM_iter(im_vectors.TGV_ROSEM_apu, options.lam_rosem, iter, options.beta_TGV_rosem, grad, options.epps);
                if verbose
                    tElapsed = toc(tStart);
                    disp(['ROSEM TGV iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                else
                    disp(['ROSEM TGV iteration ' num2str(iter) ' finished'])
                end
            end
            
            if options.NLM && options.BSREM
                if verbose
                    tStart = tic;
                end
                med = NLM(im_vectors.NLM_BSREM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                    options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
                im_vectors.NLM_BSREM_apu = BSREM_iter(im_vectors.NLM_BSREM_apu, options.lam_rosem, iter, options.beta_NLM_bsrem, med, options.epps);
                if verbose
                    tElapsed = toc(tStart);
                    disp(['BSREM NLM iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                else
                    disp(['BSREM NLM iteration ' num2str(iter) ' finished'])
                end
            end
            if options.NLM && options.ROSEM_MAP
                if verbose
                    tStart = tic;
                end
                med = NLM(im_vectors.NLM_ROSEM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                    options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
                im_vectors.NLM_ROSEM_apu = BSREM_iter(im_vectors.NLM_ROSEM_apu, options.lam_rosem, iter, options.beta_NLM_rosem, med, options.epps);
                if verbose
                    tElapsed = toc(tStart);
                    disp(['ROSEM NLM iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                else
                    disp(['ROSEM NLM iteration ' num2str(iter) ' finished'])
                end
            end
            
            disp(['Iteration ' num2str(iter) ' finished'])
        end
    end
    
end