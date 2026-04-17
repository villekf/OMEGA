%% MATLAB/Octave codes for dynamic PET reconstruction using Inveon PET LST
% This file provides an example of dynamic reconstruction in PET, using
% either sinogram data or listmode data.
% This example uses BOTH spatial AND temporal regularization (RDP for
% spatial and TV for temporal). You can disable RDP by setting options.RDP
% to false below and the temporal TV by setting options.temporalTV = false.
% By default, this example uses sinogram data, but you can use listmode
% data by modifying line 1022 variable listmode to true.
% For the input measurement data, you can use the open preclinical PET data
% available from: https://doi.org/10.5281/zenodo.3528056
% Documentation: https://omega-doc.readthedocs.io/en/latest/dynamic.html
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SCANNER PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% R-sectors/modules/blocks/buckets in transaxial direction
options.blocks_per_ring = (16);

%%% R-sectors/modules/blocks/buckets in axial direction (i.e. number of physical
%%% scanner/crystal rings)
% Multiplying this with the below cryst_per_block_axial should equal the
% total number of crystal rings.
options.linear_multip = (4);

%%% Number of detectors on the side of R-sector/block/module (transaxial
%%% direction)
% (e.g. 13 if 13x13, 20 if 20x10)
options.cryst_per_block = (20);

%%% Number of detectors on the side of R-sector/block/module (axial
%%% direction)
% (e.g. 13 if 13x13, 10 if 20x10)
options.cryst_per_block_axial = 20;

%%% Crystal pitch/size in x- and y-directions (transaxial) (mm)
options.cr_p = 1.59;

%%% Crystal pitch/size in z-direction (axial) (mm)
options.cr_pz = 1.59;

%%% Ring diameter (distance between perpendicular detectors) (mm)
options.diameter = 161;

%%% Transaxial FOV size (mm), this is the length of the x (vertical/row) side
% of the FOV
options.FOVa_x = 100;

%%% Transaxial FOV size (mm), this is the length of the y (horizontal/column) side
% of the FOV
options.FOVa_y = options.FOVa_x;

%%% Axial FOV (mm)
options.axial_fov = 127;

%%% Ring gaps (mm)
% Each ring is assumed to contain options.cryst_per_block_axial crystals
% Input the gap between each of these rings here, for every gap
% If there are no gaps, leave this empty or zero
% If the gap values are the same, you need to repeat the value for each gap
options.ringGaps = [];

%%% Number of detectors per crystal ring (without pseudo detectors)
options.det_per_ring = options.blocks_per_ring*options.cryst_per_block;

%%% Number of crystal rings
options.rings = options.linear_multip * options.cryst_per_block;

%%% Total number of detectors
options.detectors = options.det_per_ring*options.rings;

%%% Scanner name
% Used for naming purposes (measurement data)
options.machine_name = 'Inveon';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% NOTE! The below sinogram properties are only needed for sinogram data,
%%% not for listmode. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SINOGRAM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Span factor/axial compression
options.span = 3;

%%% Maximum ring difference
options.ring_difference = options.rings - 1;

%%% Number of radial positions (views) in sinogram
% You should primarily use the same number as the scanner uses.
% However, if that information is not available you can use ndist_max
% function to determine potential values (see help ndist_max for usage).
% This is the ROW dimension, i.e. the number of rows in the sinogram
options.Ndist = 128;

%%% Number of angles (tangential positions) in sinogram
% This is the final amount after possible mashing, maximum allowed is the
% number of detectors per ring/2.
% This is the COLUMN dimension, i.e. the number of columns in the sinogram
options.Nang = 160;

%%% Specify the amount of sinograms contained on each segment
% (this should total the total number of sinograms).
% Currently this is computed automatically, but you can also manually
% specify the segment sizes.
options.segment_table = [options.rings*2-1, options.rings*2-1 - (options.span + 1):-options.span*2:options.rings - options.ring_difference];
if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
    options.segment_table = [options.segment_table(1), repeat_elem(options.segment_table(2:end),2,1)];
else
    options.segment_table = [options.segment_table(1), repelem(options.segment_table(2:end),2)];
end

%%% Total number of sinograms
options.TotSinos = sum(options.segment_table);

%%% Number of sinograms used in reconstruction
% The first NSinos sinograms will be used for the image reconstruction.
% This is the number of slices in the sinogram
options.NSinos = options.TotSinos;

%%% If Ndist value is even, take one extra out of the negative side (+1) or
% from the positive side (-1). E.g. if Ndist = 200, then with +1 the
% interval is [-99,100] and with -1 [-100,99]. This varies from scanner to
% scanner. If you see a slight shift in the sinograms when comparing with
% the scanner sinograms then use the other option here. For Inveon, this
% should be -1.
options.ndist_side = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INVEON DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Use Inveon data
% 0/false = Don't use Inveon scanner data or use user-input data
% 1 = Use list-mode data (.lst)
% 2 = Use scanner created sinogram data (.scn)
options.use_machine = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% IMAGE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note that non-square transaxial image sizes can be unreliable just as the
% non-square transaxial FOV, but they should, generally, work
%%% Reconstructed image pixel count (X/row-direction)
options.Nx = 128;

%%% Y/column-direction
options.Ny = 128;

%%% Z-direction (number of slices) (axial)
options.Nz = options.rings*2 - 1;

%%% Flip the image (in column direction)?
options.flip_image = false;

%%% How much is the image rotated?
% NOTE: The rotation is done in the detector space (before reconstruction).
% This current setting is for systems whose detector blocks start from the
% right hand side when viewing the scanner from front.
% Positive values perform the rotation in clockwise direction
% The units are crystals, i.e. if the value is 1, the rotation is done by
% rotating the coordinates equaling to one crystal pitch
options.offangle = options.det_per_ring * (2/4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CORRECTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List-mode reconstruction only supports automatic attenuation correction.
% Randoms and/or scatter has to be done manually. Normalization is not yet
% supported, but you can perform it manually. For sinogram, all corrections
% are supported.

%%%%%%%%%%%%%%%%%%%%%%%%% Attenuation correction %%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Image-based attenuation correction
% Include attenuation correction from images (e.g. CT-images) (for this you
% need attenuation images of each slice correctly rotated and scaled for
% 511 keV). 
% You can either use the path below to input the data or manually input
% the attenuation data into options.vaimennus
% For dynamic data, the attenuation can either be static (same for all
% timesteps) or dynamic (different for each timestep). For the latter, make
% sure that the image is of size Nx*Ny*Nz*numberOfTimesteps
options.attenuation_correction = true;

%%% CT-image attenuation
% Use CT-images for the attenuation. If set to false, uses the
% .atn-files instead (if above attenuation is set to true).
options.CT_attenuation = true;

%%% Attenuation image data file
% Specify the path (if not in MATLAB path) and filename.
% NOTE: the attenuation data must be the only variable in the file and
% should have the dimensions of the final reconstructed image. Previously
% saved attenuation images can be used here.
options.attenuation_datafile = 'Inveon_attenuation_coefficients_for_open_PET_data_128x128x159.mat';
 

%%%%%%%%%%%%%%%%%%%%%%%% Normalization correction %%%%%%%%%%%%%%%%%%%%%%%%%
% Sinogram data only!
%%% Apply normalization correction
% If set to true, normalization correction is applied in either data
% formation or in the image reconstruction by using precomputed 
% normalization coefficients. I.e. once you have computed the normalization
% coefficients, turn above compute_normalization to false and set this to
% true. Alternatively, input your own normalization data (see below)
options.normalization_correction = false;

%%% Use user-made normalization
% Use either a .mat or .nrm file containing the normalization coefficients
% for normalization correction, or input the normalization data into 
% options.normalization if normalization_correction is also set to true
% User will be prompted for the location of the file either during sinogram
% formation or before image reconstruction (see below).
% NOTE: If you have previously computed normalization coefficients with
% OMEGA, you do not need to set this to true. The normalization
% coefficients for the specified scanner will be automatically loaded. Use
% this only if you want to use normalization coefficients computed outside
% of OMEGA.
% NOTE: Supports .nrm files created by the Inveon AW. 
% Use the Inveon nrm-file here
options.use_user_normalization = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Global corrections %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Global correction factor
% This correction factor will be applied (if nonzero) to all LORs equally.
% This can be e.g. dead time correction factor.
options.global_correction_factor = [];


%%%%%%%%%%%%%%%%%%%% Corrections during reconstruction %%%%%%%%%%%%%%%%%%%%
% If set to true, all the corrections are performed during the
% reconstruction step, otherwise the corrections are performed to the
% sinogram/raw data before reconstruction (i.e. precorrected). I.e. this
% can be considered as e.g. normalization weighted reconstruction if
% normalization correction is applied.
% NOTE: Attenuation correction is always performed during reconstruction
% regardless of the choice here.
options.corrections_during_reconstruction = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% DYNAMIC IMAGING PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Total time of the measurement (s)
% Use inf if you want the whole examination (static measurement only)
% Note that this value is only used when LOADING data from GATE or Inveon
% files using OMEGA's built-in functions
% In this example we use a 30 min measurement
options.tot_time = 60 * 30;

%%% Number of time points/dynamic frames (if a static measurement, use 1)
%%% or alternatively the size of the time window for each dynamic step (in
%%% seconds). If you use a scalar value that is bigger than 1, then the
%%% dynamic time windows will use a constant width. If you want to use
%%% custom time windows you can use them e.g. with options.partitions =
%%% [30;30;60;120]; where each element is the width of the time window (in
%%% seconds). Note that the sum should in this case equal the end time
%%% minus the start time.
% NOTE: The above applies ONLY when using OMEGA to load the data. If you
% use your own data, this should be the number of time steps!
% In this example case, we first specify a 30 second timestep at the
% beginning (there's very little activity between 0-20 seconds) and then
% have 27 timesteps that are all 10 seconds long. After this, five
% timesteps of 300 seconds. There's no particular reason for this specific
% number of timesteps here and you can use any (nonzero) length.
options.partitions = [30, repelem(10,27), repelem(300,5)];

%%% Start time (s) (all measurements BEFORE this will be ignored)
% Note that this value is only used when LOADING data from GATE or Inveon
% files using OMEGA's built-in functions
options.start = 0;

%%% End time (s) (all measurements AFTER this will be ignored)
% Use inf if you want to the end of the examination (static measurement
% only)
% Note that this value is only used when LOADING data from GATE or Inveon
% files using OMEGA's built-in functions
options.end = options.tot_time;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MISC PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Name of current datafile/examination
% This is used to name the saved measurement data and also load it in
% future sessions.
options.name = 'dynamic_PET_data';

%%% Folder for the data (.dat ASCII, .root ROOT) files
% If no files are located in the path provided below, then the current
% folder is also checked. If no files are detected there either, an error
% is thrown.
% NOTE: for .lst or .scn files the user will be prompted for their
% locations and as such this path is ignored.
if ispc % Windows
    options.fpath = 'C:\path\to\GATE\output\';
else % Unix/Mac
    options.fpath = '/path/to/GATE/output/';
end

%%% Compute only the reconstructions
% If this file is run with this set to true, then the data load and
% sinogram formation steps are always skipped. Normalization coefficients
% are not computed even if selected.
options.only_reconstructions = false;

%%% Show status messages
% These are e.g. time elapsed on various functions and what steps have been
% completed. It is recommended to keep this at 1 or 2. With value of 2,
% you get more detailed timing information. Maximum is 3, minimum 0.
options.verbose = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPLEMENTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reconstruction implementation used
% 1 = Reconstructions in MATLAB (projector in a MEX-file), uses matrices.
% (Slow and memory intensive)
% 2 = Matrix-free reconstruction with OpenCL/CUDA (Recommended)
% (Requires ArrayFire).
% 3 = Multi-GPU/device matrix-free OpenCL (OSEM & MLEM only).
% 4 = Matrix-free reconstruction with OpenMP (CPU, parallel), standard C++
% 5 = Matrix-free reconstruction with OpenCL (parallel)
% See the docs for more information:
% https://omega-doc.readthedocs.io/en/latest/implementation.html
options.implementation = 2;

% Applies to implementations 3 and 5 ONLY
%%% OpenCL platform used
% NOTE: Use OpenCL_device_info() to determine the platform numbers and
% their respective devices with implementations 3 or 5.
options.platform = 0;

% Applies to implementations 2, 3 and 5 ONLY
%%% OpenCL/CUDA device used
% NOTE: Use ArrayFire_OpenCL_device_info() to determine the device numbers
% with implementation 2.
% NOTE: Use OpenCL_device_info() to determine the platform numbers and
% their respective devices with implementations 3 or 5.
% NOTE: The device numbers might be different between implementation 2 and
% implementations 3 and 5
% NOTE: if you switch devices then you might need to run the below line
% (uncommented) as well:
% clear mex
options.use_device = 0;

% Implementation 2 ONLY
%%% Use CUDA
% Selecting this to true will use CUDA kernels/code instead of OpenCL. This
% only works if the CUDA code was successfully built. This is recommended
% if you have CUDA-capable device.
options.use_CUDA = false;

% Implementation 2 ONLY
%%% Use CPU
% Selecting this to true will use CPU-based code instead of OpenCL or CUDA.
% Not recommended, even OpenCL with CPU should be used before this.
options.use_CPU = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROJECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Type of projector to use for the geometric matrix
% 1 = Improved/accelerated Siddon's algorithm
% 2 = Orthogonal distance based ray tracer
% 3 = Volume of intersection based ray tracer
% 4 = Interpolation-based projector
% NOTE: You can mix and match most of the projectors. I.e. 41 will use
% interpolation-based projector for forward projection while improved
% Siddon is used for backprojection.
% NOTE 2: The below additional options apply also in hybrid cases as long
% as the other projector is the corresponding projector.
% See the documentation for more information:
% https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
options.projector_type = 1;

%%% Use mask
% The mask needs to be a binary mask (uint8 or logical) where 1 means that
% the pixel is included while 0 means it is skipped. Separate masks can be
% used for both forward and backward projection and either one or both can
% be utilized at the same time. E.g. if only backprojection mask is input,
% then only the voxels which have 1 in the mask are reconstructed.
% The mask can be either a 2D image that is applied identically to each slice
% or a 3D mask that is applied as-is
% Forward projection mask
% If nonempty, the mask will be applied. If empty, or completely omitted, no
% mask will be considered.
% options.maskFP = true(options.Ndist,options.Nang);
% Backprojection mask
% If nonempty, the mask will be applied. If empty, or completely omitted, no
% mask will be considered.
% Create a circle that fills the FOV:
% [columnsInImage, rowsInImage] = meshgrid(1:options.Nx, 1:options.Ny);
% centerX = options.Nx/2;
% centerY = options.Ny/2;
% radius = options.Nx/2;
% options.maskBP = uint8((rowsInImage - centerY).^2 ...
%     + (columnsInImage - centerX).^2 <= radius.^2);

%%% Interpolation length (projector type = 4 only)
% This specifies the length after which the interpolation takes place. This
% value will be multiplied by the voxel size which means that 1 is
% the interpolation length corresponding to a single voxel (transaxial)
% length. Larger values lead to faster computation but at the cost of
% accuracy. Recommended values are between [0.5 1].
options.dL = 0.5;

%%% Use point spread function (PSF) blurring
% Applies PSF blurring through convolution to the image space. This is the
% same as multiplying the geometric matrix with an image blurring matrix.
options.use_psf = true;

% FWHM (mm) of the Gaussian used in PSF blurring in all three dimensions
options.FWHM = [options.cr_p options.cr_p options.cr_pz];

% Orthogonal ray tracer (projector_type = 2 only)
%%% The 2D (XY) width (mm) of the "strip/tube" where the orthogonal distances are
% included. If tube_width_z below is non-zero, then this value is ignored.
options.tube_width_xy = options.cr_p;

% Orthogonal ray tracer (projector_type = 2 only)
%%% The 3D (Z) width (mm) of the "tube" where the orthogonal distances are
% included. If set to 0, then the 2D orthogonal ray tracer is used. If this
% value is non-zero then the above value is IGNORED.
% If you want the projector to be a tube, use this, if you want it to be
% strip, use the above
% This slows down the reconstruction, but makes it more accurate
options.tube_width_z = options.cr_pz;

% Volume ray tracer (projector_type = 3 only)
%%% Radius of the tube-of-response (cylinder)
% The radius (mm) of the cylinder that approximates the tube-of-response.
options.tube_radius = sqrt(2) * (options.cr_pz / 2);

% Volume ray tracer (projector_type = 3 only)
%%% Relative size of the voxel (sphere)
% In volume ray tracer, the voxels are modeled as spheres. This value
% specifies the relative radius of the sphere such that with 1 the sphere
% is just large enough to encompass an entire cubic voxel, i.e. the
% corners of the cubic voxel intersect with the sphere shell. Larger values
% create larger spheres, while smaller values create smaller spheres.
options.voxel_radius = 1;

% projector_type = 1 and 4 only
%%% Number of rays
% Number of rays used per detector if projector_type = 1 (i.e. Improved
% Siddon is used) or projector_type = 4 (interpolation).
% The total number of rays per detector is the multiplication of the two
% below values!
% Number of rays in transaxial (row) direction
options.n_rays_transaxial = 1;
% Number of rays in axial (column) direction
options.n_rays_axial = 1;

%%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Number of iterations (all reconstruction methods)
options.Niter = 2;

%%% Save specific intermediate iterations
% You can specify the intermediate iterations you wish to save here. Note
% that this uses zero-based indexing, i.e. 0 is the first iteration (not
% the initial value). By default only the last iteration is saved.
options.saveNIter = [];
% Alternatively you can save ALL intermediate iterations by setting the
% below to true and uncommenting it
% options.save_iter = false;

%%% Number of subsets (excluding subset_type = 6)
options.subsets = 8;

%%% Subset type (n = subsets)
% Note that you can use more subset types with sinogram data, but this
% example covers both sinogram and listmode usage. Listmode data can only
% use subset types 0 and 1 with dynamic data.
% 0 = Measurements are divided into n segments
% 1 = Every nth measurement is taken
options.subset_type = 1;

%%% Initial value for the reconstruction
options.x0 = ones(options.Nx, options.Ny, options.Nz);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION ALGORITHMS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction algorithms to use (choose only one algorithm and
% optionally one prior)

%%% In this example we use PKMA, in order to use both spatial and temporal
%%% regularization. You can freely deselect one or both of them.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ML-METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are non-regularized versions
%%% Ordered Subsets Expectation Maximization (OSEM) OR Maximum-Likelihood
%%% Expectation Maximization (MLEM) (if subsets = 1)
% Supported by all implementations
options.OSEM = false;

%%% Accelerated COSEM (ACOSEM)
% Supported by implementations 1, 2, 4, and 5
options.ACOSEM = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAP-METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Any algorithm selected here will utilize any of the priors selected below
% this. Note that only one algorithm and prior combination is allowed! You
% can also use most of these algorithms without priors (such as PKMA or
% PDHG).
%%% One-Step Late OSEM (OSL-OSEM) or MLEM (if subsets = 1)
% Supported by implementations 1, 2, 4, and 5
options.OSL_OSEM = false;

%%% Modified BSREM (MBSREM)
% Supported by implementations 1, 2, 4, and 5
options.MBSREM = false;

%%% Block Sequential Regularized Expectation Maximization (BSREM)
% Supported by implementations 1, 2, 4, and 5
options.BSREM = false;

%%% Preconditioned Krasnoselskii-Mann algorithm (PKMA)
% Supported by implementations 1, 2, 4, and 5
options.PKMA = true;

%%% Primal-dual hybrid gradient (PDHG)
% Supported by implementations 1, 2, 4, and 5
options.PDHG = false;

%%% Primal-dual hybrid gradient (PDHG) with Kullback-Leibler minimization
% Supported by implementations 1, 2, 4, and 5
options.PDHGKL = false;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPATIAL PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Median Root Prior (MRP)
options.MRP = false;

%%% Quadratic Prior (QP)
options.quad = false;

%%% Huber Prior (QP)
options.Huber = false;

%%% Total Variation (TV) prior
options.TV = false;

%%% Non-local Means (NLM) prior
options.NLM = false;

%%% Relative difference prior
% In this case, we use the RDP regularization as the spatial regularization
options.RDP = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEMPORAL PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First derivative smoothness prior
options.temporal_smoothness = false;

% Temporal total variations prior
% Temporal TV is used in this example and recommended in general
options.temporalTV = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%% ENFORCE POSITIVITY %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Applies to PDHG, PDHGL1, PDDY, FISTA, FISTAL1, MBSREM, MRAMLA, PKMA
% Enforces positivity in the estimate after each iteration
options.enforcePositivity = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%% ACOSEM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Acceleration parameter for ACOSEM (1 equals COSEM)
options.h = 2;


%%%%%%%%%%%%%%%%%%%%%%%%%% RELAXATION PARAMETER %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Relaxation parameter for MRAMLA, RAMLA, ROSEM, BSREM, MBSREM and PKMA
% Use scalar if you want it to decrease as
% lambda / ((current_iteration - 1)/20 + 1). Use vector (length = Niter) if
% you want your own relaxation parameters. Use empty array or zero if you
% want to OMEGA to compute the relaxation parameter using the above formula
% with lambda = 1. Note that current_iteration is one-based, i.e. it starts
% at 1.
options.lambda = 1;


%%%%%%%%%%%%%%%%%%%%%%%% MRAMLA & MBSREM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%
%%% Upper bound for MRAMLA/MBSREM (use 0 for default (computed) value)
options.U = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PKMA PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step size (alpha) parameter for PKMA
% If a scalar (or an empty) value is used, then the alpha parameter is
% computed automatically as alpha_PKMA(oo) = 1 + (options.rho_PKMA *((i -
% 1) * options.subsets + ll)) / ((i - 1) * options.subsets + ll +
% options.delta_PKMA), where i is the iteration number and l the subset
% number. The input number thus has no effect. options.rho_PKMA and
% options.delta_PKMA are defined below.
% If, on the other hand, a vector is input then the input alpha values are
% used as is without any modifications (the length has to be at least the
% number of iterations * number of subsets).
options.alpha_PKMA = 0;

%%% rho_PKMA
% This value is ignored if a vector input is used with alpha_PKMA
options.rho_PKMA = 0.95;

%%% delta_PKMA
% This value is ignored if a vector input is used with alpha_PKMA
options.delta_PKMA = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PDHG PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Primal value
% If left zero, or empty, it will be automatically computed
% Note that if you change any of the model parameters, i.e. image volume
% size, number of projections or use binning, this needs to be recomputed
% or scaled accordingly!
% The computed largest eigenvalue is printed if verbose > 0. This can be
% used as the below value as long as one is divided by it. For example,
% if "Largest eigenvalue for volume 0 is 100" then options.tauCP should be
% 1/100 (if you use filtering-based preconditioner this is the "without
% filtering" value)
% if you have a multi-resolution situation, you should input the values
% for each volume or use zero/empty
options.tauCP = 0;
% Primal value for filtered iterations, applicable only if
% options.precondTypeMeas[2] = true. As with above, automatically computed
% if left zero or empty. Same restrictions apply here as above.
% Use the "Largest eigenvalue for volume 0 with filtering" value here!
% if you have a multi-resolution situation, you should input the values
% for each volume or use zero/empty
options.tauCPFilt = 0;
% Dual value. Recommended to set at 1.
options.sigmaCP = 1;
% Next estimate update variable
options.thetaCP = 1;
% Dual value for TV and/or TGV. For faster convergence, set this to higher
% than 1.
options.sigma2CP = 1;

% Use adaptive update of the primal and dual variables
% Currently two methods available
% Setting this to 1 or 2 uses an adaptive update for both the primal and
% dual variables.
% Can lead to unstable behavior when using with multi-resolution
% Minimal to none use with filtering-based preconditioner
options.PDAdaptiveType = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRECONDITIONERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Applies to PDHG, PDHGL1, PDHGKL, PKMA, MBSREM, MRAMLA, PDDY, FISTA,
%%% FISTAL1, and SAGA
% Measurement-based preconditioners
% precondTypeMeas(1) = Diagonal normalization preconditioner (1 / (A1))
% precondTypeMeas(2) = Filtering-based preconditioner
options.precondTypeMeas = [false;false];

% Image-based preconditioners
% Setting options.precondTypeImage(2) = true when using PKMA, MRAMLA or
% MBSREM is recommended
% precondTypeImage(1) = Diagonal normalization preconditioner (division with
% the sensitivity image 1 / (A^T1), A is the system matrix)
% precondTypeImage(2) = EM preconditioner (f / (A^T1), where f is the current
% estimate)
% precondTypeImage(3) = IEM preconditioner (max(n, fhat, f)/ (A^T1), where
% fhat is an estimate of the final image and n is a small positive number)
% precondTypeImage(4) = Momentum-like preconditioner (basically a step size
% inclusion)
% precondTypeImage(5) = Gradient-based preconditioner (Uses the normalized
% divergence (sum of the gradient) of the current estimate)
% precondTypeImage(6) = Filtering-based preconditioner
% precondTypeImage(7) = Curvature-based preconditioner
options.precondTypeImage = [false;false;false;false;false;false;false];
if options.PKMA || options.MBSREM
    options.precondTypeImage(2) = true;
end

% Reference image for precondTypeImage(3). Can be either a mat-file or a
% variable
options.referenceImage = '';

% Momentum parameter for precondTypeImage(4)
% Set the desired momentum parameters to the following variable (note that
% the length should be options.Niter * options.subsets):
% options.alphaPrecond = [];
% Otherwise set the following parameters:
options.rhoPrecond = options.rho_PKMA;
options.delta1Precond = options.delta_PKMA;

% Parameters for precondTypeImage(5)
% See the article for details:
% https://omega-doc.readthedocs.io/en/latest/algorithms.html#gradient-based-preconditioner
options.gradV1 = 1.5;
options.gradV2 = 2;
% Note that these include subiterations (options.Niter * options.subsets)
% The first iteration where to start the gradient computation
options.gradInitIter = 1;
% Last iteration of the gradient computation
options.gradLastIter = 100;

% Number of filtering iterations
% Applies to both precondTypeMeas(2) and precondTypeImage(6)
% The filtering is applies to this many (sub)iterations
% Note that this include subiterations (options.Niter * options.subsets)
options.filteringIterations = 100;


%%%%%%%%%%%%%%%%%%%%% SPATIAL REGULARIZATION PARAMETER %%%%%%%%%%%%%%%%%%%%
%%% The regularization parameter for ALL spatial regularization methods
%%% (priors)
% Adjusted for sinogram data
options.beta = 0.1; % 0.1 is good for listmode


%%%%%%%%%%%%%%%%%%%% TEMPORAL REGULARIZATION PARAMETER %%%%%%%%%%%%%%%%%%%%
%%% The regularization parameter for ALL temporal regularization methods
%%% (priors)
% Adjusted for sinogram data
options.beta_temporal = 0.01; % 0.1 is good for listmode


%%%%%%%%%%%%%%%%%%%%%%%%% NEIGHBORHOOD PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%
%%% How many neighboring pixels are considered
% With MRP, QP, L, FMH, NLM, (RDP), GGMRF and weighted mean
% E.g. if Ndx = 1, Ndy = 1, Ndz = 0, then you have 3x3 square area where
% the pixels are taken into account (I.e. (Ndx*2+1)x(Ndy*2+1)x(Ndz*2+1)
% area).
% NOTE: Currently Ndx and Ndy must be identical.
% For NLM this is often called the "search window".
options.Ndx = 2;
options.Ndy = 2;
options.Ndz = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pixel weights for quadratic prior
% The number of pixels need to be the amount of neighboring pixels,
% e.g. if the above Nd values are all 1, then 27 weights need to be
% included where the center pixel (if Nd values are 1, element 14) should
% be Inf. Size is (Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1). If left empty then
% they will be calculated by the algorithm and are based on the distance of
% the voxels from the center.
options.weights = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Delta parameter for Huber prior
% Upper and lower bounds for the prior
options.huber_delta = 5;

%%% Pixel weights for Huber prior
% Same rules apply as with quadratic prior weights.
% If left empty then they will be calculated by the algorithm and are based
% on the distance of the voxels from the center.
options.weights_huber = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TV PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% "Smoothing" parameter
% Also used to prevent zero values in square root.
options.TVsmoothing = 1e-5;

%%% Whether to use an anatomical reference/weighting image with the TV
options.TV_use_anatomical = false;

%%% If the TV_use_anatomical value is set to true, specify filename for the
% reference image here (same rules apply as with attenuation correction
% above). Alternatively you can specify the variable that holds the
% reference image.
options.TV_reference_image = 'reference_image.mat';

%%% Five different TV methods are available.
% Value can be 1, 2, 3, 4 or 6.
% Type 3 is not recommended!
% Types 1 and 2 are the same if anatomical prior is not included
% Type 3 uses the same weights as quadratic prior
% Type 4 is the Lange prior, does not support anatomic weighting.
% Type 6 is a weighted TV, does not support anatomic weighting.
% See the docs for more information:
% https://omega-doc.readthedocs.io/en/latest/algorithms.html#tv
options.TVtype = 1;

%%% Weighting parameters for the TV prior.
% Applicable only if use_anatomical = true. T-value is specific to the used
% TVtype, e.g. for type 1 it is the edge threshold parameter. See the docs
% for more details:
% https://omega-doc.readthedocs.io/en/latest/algorithms.html#tv
options.T = 0.5;

%%% C is the weight for the original image in type 3 and is ignored with
% other types
options.C = 1;

%%% Tuning parameter for TV and APLS
options.tau = 1e-8;

%%% Tuning parameter for Lange function in SATV (type 4) or weight factor
%%% for weighted TV (type 6)
% Setting this to 0 gives regular anisotropic TV with type 4
% This affects also non-local Lange
options.SATVPhi = 0.2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NLM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filter parameter
% Higher values smooth the image, smaller values make it sharper
options.sigma = 10;

%%% Patch radius
% Works exactly the same as the neighborhood size
options.Nlx = 1;
options.Nly = 1;
options.Nlz = 1;

%%% Standard deviation of the Gaussian-weighted Euclidean norm
options.NLM_gauss = 2;

% Search window radius is controlled by Ndx, Ndy and Ndz parameters
% Use anatomical reference image for the patches
options.NLM_use_anatomical = false;

%%% Specify filename for the reference image here (same rules apply as with
% attenuation correction above) or the variable containing the reference image
options.NLM_reference_image = 'reference_image.mat';

% Note that only one of the below options for NLM can be selected!
% If all the below ones are false, regular NLM is used!
%%% Use Non-local total variation (NLTV)
options.NLTV = false;

%%% Use Non-local Lange prior (NLLange)
options.NLLange = false;

%%% Use MRP algorithm (without normalization)
% I.e. gradient = im - NLM_filtered(im)
options.NLM_MRP = false;

%%% Use non-local relative difference prior (NLRD)
options.NLRD = false;

%%% Use non-local GGMRF (NLGGMRF)
options.NLGGMRF = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RDP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Edge weighting factor
% Higher values sharpen the image, smaller values make it smoother
% Note that this affects NLRD as well
options.RDP_gamma = 10;

% If true, includes also the "diagonal" corners in the neighborhood in RDP
% By default, only the sides which the current voxel shares a side are
% included
% See https://omega-doc.readthedocs.io/en/latest/algorithms.html#rdp for
% details
% Default is false
options.RDPIncludeCorners = false;

% Applies only if the above RDPIncludeCorners is true
% Use anatomical reference image weighting for RDP
options.RDP_use_anatomical = false;

% Set the file containing the reference image or the variable of the reference
% image here
options.RDP_reference_image = '';


%%%%%%%%%%%%%%%%%%%%%%%%%% TEMPORAL TV PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Smoothing parameter for temporal TV
% Prevents potential division by zero
options.temporalTVsmoothing = 1e-6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% OPENCL DEVICE INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Implementation 2
% Uncomment the below line and run it to determine the available device
% numbers
% ArrayFire_OpenCL_device_info();

%%% Implementation 3
% Uncomment the below line and run it to determine the available platforms,
% their respective numbers and device numbers
% OpenCL_device_info();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ERROR CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic error checking is done here
options = OMEGA_error_check(options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% DEPTH OF INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uncomment the below value to set a depth of interaction (mm)
% NOTE: By default this is set to 0, i.e. it is assumed that all the
% interactions occur at the surface of the detector crystals. What this
% value changes is the depth of where the interactions are assumed to
% occur, i.e. it only changes the detector coordinates such that the
% transaxial coordinates are "deeper" in the crystal.
% options.DOI = 4.584;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load the LST coincidence data
% Loads the detector coordinates or indices for each coincidence event

% Whether sinogram (false) or listmode (true) data is used
% Default is sinogram
listmode = false;

% You can either use the index-based reconstruction or the detector
% coordinates when using listmode data
% Index-based is on by default

% Use index-based reconstruction
% If true, requires options.trIndex and options.axIndex variables
% These should include the transaxial and axial, respectively, detector
% coordinate indices. Two indices are required per measurement, i.e.
% source-detector or detector-detector pairs. The indexing has to be
% zero-based! The transaxial coordinates should be stored in options.x and
% axial coordinates in options.z. The indices should correspond to the
% coordinates in each. Note that options.x should have both the x- and
% y-coordinates while options.z should have only the z-coordinates. You can
% also include randoms by inputting them as negative measurements. The
% indices are used in the same order as measurements.
options.useIndexBasedReconstruction = listmode;

% For list-mode data, when not using the index-based reconstruction defined
% above the core component you need are the detector coordinates for each
% event. options.x should be 6xNumberOfEvents, where the first three rows
% correspond to the x/y/z-coordinates of the first detector and the next
% three the x/y/z-coordinates for the second detector. options.x can also
% be a vector of size 6 * NumberOfEvents. 
% Note: List-mode reconstruction can be much more memory intensive than
% regular reconstruction

% Note: Origin is assumed to be at the center. If this is not the case, you
% can shift it with options.oOffsetX, options.oOffsetY and options.oOffsetZ

if ~options.only_reconstructions && options.use_machine ~= 2
    if listmode
        if ~options.useIndexBasedReconstruction
            % Save the coordinates of the detected coincidences to work space:
            % Replace options.x, etc. with your own coordinates to use user input
            % data.
            [~, ~, ~, ~, ~, options.x] = load_data(options);
        else
            % Index-based reconstruction
            % Load the detector indices, or LOR numbers
            [~, ~, ~, ~, ~, options.trIndex, options.axIndex] = load_data(options);
        end
    else
        % Sinogram creation
        options.SinM = load_data(options);
    end
end

if listmode && options.useIndexBasedReconstruction
    % Coordinates for the detector pairs
    % Transaxial
    [x, y] = detector_coordinates(options);
    % Input the transaxial coordinates to options.x
    options.x = [x';y'];
    % Axial coordinates
    z_length = double(options.linear_multip .* options.cryst_per_block_axial(1) * options.cr_pz);
    z1 = linspace(-(z_length / 2 - options.cr_pz/2), z_length / 2 - options.cr_pz/2, options.rings)';
    % Input the axial coordinates to options.z
    options.z = z1;
end

%% Reconstructions
% When using list-mode reconstruction, the detector coordinates can also
% correspond to the LORs in sinograms and the input measurement data should
% then be in sinogram format (matrix or vector).

% The measurement data
% For list-mode this is all ones as each coincidence is handled
% individually
% For sinogram data, we load the data earlier
% If you want to use your own data, load it here and input it to
% options.SinM
if listmode
    for t = 1 : numel(options.partitions)
        if options.useIndexBasedReconstruction == true
            options.SinM{t} = ones(size(options.trIndex{t},2),1,'uint8');
        else
            options.SinM{t} = ones(size(options.x{t},2),1,'uint8');
        end
    end
end

% The below applies only to list-mode reconstructions:
% Since the measurement vector contains all ones, but does not contain
% every possible LOR (and also contains some duplicates) the sensitivity
% image is computed for the every measured LOR. Having this value set to
% true, the sensitivity image is computed for all applicable LORs. This
% requires correct values in the scanner properties, mainly the number of
% detectors per ring, number of rings and possible (pseudo) gaps. If this is
% set to false, then the sensitivity image is computed for the input
% detector coordinates only.
% NOTE: If you use your own detector coordinates, you can include the
% necessary coordinates for sensitivity image, that are missing from the
% measurements, by adding zero measurements to the corresponding
% coordinates. Coordinates that appear multiple times, however, will cause
% the sensitivity image values to be incorrect in these coordinates as the
% values are added multiple times.
% NOTE: If you use sinogram data with index-based recon, setting this to
% false should be OK.
options.compute_sensitivity_image = true;

% In case the data takes a lot of memory and there isn't enough memory on
% the GPU to store all the measurement data, set the below value to false
% If false, loads only the current subset and timestep to the GPU at a time
options.loadTOF = true;

if options.only_sinos == false

    tStart = tic;
    pz = reconstructions_main(options);
    tElapsed = toc(tStart);
    disp(['Reconstruction process took ' num2str(tElapsed) ' seconds'])

    pz = squeeze(pz);
    % You can visualize dynamic images as well
    % The upper slider is for slices while the lower one is for timesteps
    volume3Dviewer(pz, [], [0 0 1])

end
