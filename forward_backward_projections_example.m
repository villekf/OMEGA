%% An example file demonstrating the use of the forward and backprojection functions

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% MACHINE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% R-sectors/blocks in transaxial direction
options.blocks_per_ring = (42);
%%% R-sectors/modules/blocks in axial direction (i.e. number of physical
% machine rings) 
options.linear_multip = (4);
%%% number of detectors on the side of R-sector/block/module
% (e.g. 13 if 13x13)
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
%%% number of crystal rings
options.rings = options.linear_multip * options.cryst_per_block;
%%% number of detectors
options.detectors = options.det_per_ring*options.rings;
%%% machine name
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
options.ring_difference = options.rings - 1;
%%% Number of angles (tangential positions) in sinogram 
% This is the final amount after possible mashing, maximum allowed is the
% number of detectors per ring/2
options.Nang = options.det_per_ring/2;
%%% Number of angular positions (views) in sinogram
options.Ndist = 200;
%%% Specify the amount of sinograms contained on each segment
% (this should total the total number of sinograms)
options.segment_table = [options.Nz, options.Nz - (options.span + 1):-options.span*2:options.span];
if verLessThan('matlab','8.5')
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
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CORRECTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Randoms correction %%%%%%%%%%%%%%%%%%%%%%%%%%%
% If set to true, stores the delayed coincidences during data load and
% later corrects for randoms during sinogram formation (if sinogram data)
% or during data load (if raw list-mode data)
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
% data. Corrects for scatter during sinogram formation (if sinogram data)
% or during data load (if raw list-mode data)
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

%%%%%%%%%%%%%%%%%%%%%%%%%% Attenuation correction %%%%%%%%%%%%%%%%%%%%%%%%%
% Image-based attenuation correction
% Include attenuation correction from images (e.g. CT-images) (for this you
% need attenuation images of each slice correctly rotated and scaled for
% 511 keV) 
options.attenuation_correction = false;
%%% Attenuation image data file
% specify the path (if not in MATLAB path) and filename
% NOTE: the attenuation data must be the only variable in the file and
% have the dimensions of the final reconstructed image
options.attenuation_datafile = '';

%%%%%%%%%%%%%%%%%%%%%%%%% Normalization correction %%%%%%%%%%%%%%%%%%%%%%%%
% If set to true, normalization correction is applied.
options.normalization_correction = false;
%%% Use user-made normalization
% Use either a .mat or .nrm file containing the normalization coefficients
% for normalization correction if normalization_correction is also set to
% true. 
% User will be prompted for the location of the file either during sinogram
% formation or before image reconstruction (see below)
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
options.name = 'cylpet_example';
%%% Precompute data 
% This should be done when using data from a certain machine the first time
% as it can speed up reconstruction. Especially recommended for raw
% list-mode data.
options.precompute = false;
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
% NOTE: Only the below ones are available in forward/backward projections
% 1 = Reconstructions in MATLAB (projector in a MEX-file)
% 3 = Multi-GPU/device matrix-free OpenCL (OSEM & MLEM only)
options.implementation = 3;
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
% Implementation 3 ONLY
%%% How many times more measurements/LORs are in the GPU part (applicable if
% heterogenous computing is used) 
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
options.use_Shuffle = false;
%%% Use fast sparse
% Not included in OMEGA, needs to be manually downloaded and installed.
% Download from: https://github.com/stefanengblom/stenglib
% NOTE: This applies only to implementation 1 when precompute_lor is false.
options.use_fsparse = false;
 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Precompute the necessary data

if options.precompute
    precompute_data(options);
end

%% Forward and backward projections (and normalization constant)

% Get the subset indices for the current subset type
% Index contains the subset indices (which measurements belong to each
% subset)
% n_meas the number of measurements in each subset
[index, n_meas, options.subsets] = index_maker(options.Nx, options.Ny, options.Nz, options.subsets, options.use_raw_data, ...
    options.machine_name, options, options.Nang, options.Ndist, options.TotSinos, options.NSinos);

nn = [0;cumsum(n_meas)];
if iscell(index)
    index = cell2mat(index);
end

f = ones(options.Nx*options.Ny*options.Nz,1);

load('Cylindrical_PET_example_cylpet_example_sinograms_combined_static_200x168x703_span3.mat','raw_SinM')
if options.implementation == 1
    raw_SinM = double(raw_SinM(index));
else
    raw_SinM = single(raw_SinM(index));
end

% Computes the OSEM-estimates by using the forward and backprojections (not
% recommended for an actual OSEM computation due to reduced performance)
for iter = 1 : options.Niter
    for osa_iter = 1 : options.subsets
        [y, options] = forward_project(options, index(nn(osa_iter) + 1:nn(osa_iter+1)), n_meas(osa_iter), f, [nn(osa_iter) + 1 , nn(osa_iter+1)]);
        [x, norm] = backproject(options, index(nn(osa_iter) + 1:nn(osa_iter+1)), n_meas(osa_iter), raw_SinM(nn(osa_iter) + 1:nn(osa_iter+1)) ./ (y + options.epps), ...
            nn(osa_iter) + 1:nn(osa_iter+1));
        f = (f./(norm + options.epps)).*(x + options.epps);
    end
end
f = reshape(f, options.Nx,options.Ny,options.Nz);
