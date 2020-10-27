%% An example file demonstrating the use of the forward and backprojection functions

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% MACHINE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% R-sectors/modules/blocks/buckets in transaxial direction
options.blocks_per_ring = (42);

%%% R-sectors/modules/blocks/buckets in axial direction (i.e. number of physical
%%% machine/crystal rings) 
% Multiplying this with the below cryst_per_block should equal the total
% number of crystal rings. 
options.linear_multip = (4);

%%% R-sectors/modules/blocks/buckets in transaxial direction
% Required only if larger than 1
options.transaxial_multip = 1;

%%% Number of detectors on the side of R-sector/block/module (transaxial
%%% direction)
% (e.g. 13 if 13x13, 20 if 20x10)
options.cryst_per_block = (8);

%%% Number of detectors on the side of R-sector/block/module (axial
%%% direction)
% (e.g. 13 if 13x13, 10 if 20x10)
options.cryst_per_block_axial = 8;

%%% Crystal pitch/size in x- and y-directions (transaxial) (mm)
options.cr_p = 2.4;

%%% Crystal pitch/size in z-direction (axial) (mm)
options.cr_pz = 2.4;

%%% Ring diameter (distance between perpendicular detectors) (mm)
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
options.det_per_ring = options.blocks_per_ring * options.cryst_per_block * options.transaxial_multip;

%%% Number of detectors per crystal ring (with pseudo detectors)
% If your scanner has a single pseudo detector on each (transaxial) side of
% the crystal block then simply add +1 inside the parenthesis (or uncomment
% the one below).
options.det_w_pseudo = options.blocks_per_ring * options.cryst_per_block * options.transaxial_multip;

%%% Number of crystal rings
options.rings = options.linear_multip * options.cryst_per_block_axial;

%%% Total number of detectors
options.detectors = options.det_per_ring*options.rings;

%%% Machine name
% Used for naming purposes (measurement data)
options.machine_name = 'Cylindrical_PET_example';
 
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
 
%%% Reconstructed image pixel count (X-direction)
% NOTE: Non-square image sizes (X- and Y-direction) may not work
options.Nx = 128;

%%% Y-direction
options.Ny = 128;

%%% Z-direction (number of slices) (axial)
options.Nz = options.rings*2-1;

%%% Flip the image (in vertical direction)?
options.flip_image = false;

%%% How much is the image rotated?
% You need to run the precompute phase again if you modify this
% NOTE: The rotation is done in the detector space (before reconstruction).
% This current setting is for systems whose detector blocks start from the
% right hand side when viewing the device from front.
options.offangle = options.det_w_pseudo * (3/4);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% You should primarily use the same number as the device uses.
% However, if that information is not available you can use ndist_max
% function to determine potential values (see help ndist_max for usage).
options.Ndist = 200;

%%% Number of angles (tangential positions) in sinogram 
% This is the final amount after possible mashing, maximum allowed is the
% number of detectors per ring/2.
options.Nang = options.det_per_ring/2;

%%% Specify the amount of sinograms contained on each segment
% (this should total the total number of sinograms).
% Currently this is computed automatically, but you can also manually
% specify the segment sizes.
options.segment_table = [options.Nz, options.Nz - (options.span + 1):-options.span*2:max(options.Nz - options.ring_difference*2, options.rings - options.ring_difference)];
if exist('OCTAVE_VERSION','builtin') == 0 && exist('repelem', 'builtin') == 0
    options.segment_table = [options.segment_table(1), repeat_elem(options.segment_table(2:end),2,1)];
else
    options.segment_table = [options.segment_table(1), repelem(options.segment_table(2:end),2)];
end

%%% Total number of sinograms
options.TotSinos = sum(options.segment_table);

%%% Number of sinograms used in reconstruction
% The first NSinos sinograms will be used for the image reconstruction.
options.NSinos = options.TotSinos;

%%% If Ndist value is even, take one extra out of the negative side (+1) or
% from the positive side (-1). E.g. if Ndist = 200, then with +1 the
% interval is [-99,100] and with -1 [-100,99]. This varies from device to
% device. If you see a slight shift in the sinograms when comparing with
% the machine sinograms then use the other option here.
options.ndist_side = 1;

%%% Increase the sampling rate of the sinogram
% Increasing this interpolates additional rows to the sinogram.
% Can be used to prevent aliasing artifacts.
% NOTE: Has to be either 1 or divisible by two.
options.sampling = 1;

%%% Interpolation method used for sampling rate increase
% All the methods are available that are supported by interp1 (see help
% interp1).
options.sampling_interpolation_method = 'linear';

%%% Fill the gaps caused by pseudo detectors?
% NOTE: Applicable only if options.pseudot > 0
options.fill_sinogram_gaps = false;

%%% Which method used to fill the gaps?
% Either MATLAB's built-in fillmissing or inpaint_nans from file exchange.
% For inpaint_nans see: 
% https://se.mathworks.com/matlabcentral/fileexchange/4551-inpaint_nans
% See wiki for more details.
% NOTE: GNU Octave does not support fillmissing.
options.gap_filling_method = 'fillmissing';

%%% Interpolation method used with fillmissing
% Possible methods are those listed under method-section in fillmissing.
options.interpolation_method_fillmissing = 'linear';

%%% Interpolation method used with inpaint_nans
% See inpaint_nans.m for details.
options.interpolation_method_inpaint = 0;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% RAW DATA PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Use raw data
% This means that the data is used as is without any sinogramming and thus
% without any "compression". Measurement data is stored as diagonal matrix
% containing the counts on every LOR available. Raw data can be visualized
% with visualizeRawData function.
options.use_raw_data = false;

%%% Store raw data
% If the above use_raw_data is set to false, you can still save the raw
% data during data load by setting this to true.
options.store_raw_data = false;
 
%%% Maximum ring difference in raw data
options.ring_difference_raw = options.rings;

%%% Increase the sampling rate of the raw data
% Increasing this interpolates additional rows and columns to the raw data.
% Can be used to prevent aliasing artifacts.
% NOTE: Has to be either 1 or divisible by two.
options.sampling_raw = 1;

%%% Interpolation method used for sampling rate increase
% All the methods are available that are supported by interp2 (see help
% interp2).
options.sampling_interpolation_method_raw = 'linear';
 
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
 
%%%%%%%%%%%%%%%%%%%%%%%%%%% Randoms correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If set to true, stores the delayed coincidences during data load and
% later corrects for randoms during the data formation/load or during
% reconstruction. Delayes need to be stored in GATE data for this to work.
options.randoms_correction = false;

%%% Variance reduction
% If set to true, variance reduction will be performed to delayed
% coincidence (randoms corrections) data if randoms correction is selected.
options.variance_reduction = false;

%%% Randoms smoothing
% If set to true, applies a 8x8 moving mean smoothing to the delayed
% coincidence data. This is applied on all cases (i.e. randoms correction
% data is smoothed before subtraction of before reconstruction).
% NOTE: Mean window size can be adjusted by modifying the randoms_smoothing
% function.
options.randoms_smoothing = false;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Scatter correction %%%%%%%%%%%%%%%%%%%%%%%%%%%
% If set to true, will prompt the user to load the scatter sinogram/raw
% data. Corrects for scatter during data formation/load or during
% reconstruction.
% NOTE: Scatter data is not created by this software and as such must be
% provided by the user. Previously created scatter sinogram/raw data matrix
% can be used though.
options.scatter_correction = false;

%%% Variance reduction
% If set to true, variance reduction will be performed to scatter data if
% scatter correction is selected.
options.scatter_variance_reduction = false;

%%% Scatter normalization
% If set to true, normalizes the scatter coincidences data during data
% formation or before reconstruction. If set to false, the scatter data is
% subtracted from the sinogram before normalization (and when
% options.corrections_during_reconstruction = false).
options.normalize_scatter = false;

%%% Scatter smoothing
% If set to true, applies a 8x8 moving mean smoothing to the scattered
% coincidences data. This is applied on all cases (i.e. scatter correction
% data is smoothed before subtraction of before reconstruction).
% NOTE: Mean window size can be adjusted by modifying the randoms_smoothing
% function.
options.scatter_smoothing = false;
 
%%% Subtract scatter
% If set to true, the scattered coincidences are subtracted from the
% sinogram or the forward projection. If set to false, the scattered
% coincidences are multiplied with the sinogram or included to the system
% matrix. The latter choices are applied if
% options.corrections_during_reconstruction = true.
options.subtract_scatter = true;
 
%%%%%%%%%%%%%%%%%%%%%%%%% Attenuation correction %%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Image-based attenuation correction
% Include attenuation correction from images (e.g. CT-images) (for this you
% need attenuation images of each slice correctly rotated and scaled for
% 511 keV). For CT-images you can use attenuationCT_to_511 or
% create_atten_matrix_CT functions.
options.attenuation_correction = false;

%%% Attenuation image data file
% Specify the path (if not in MATLAB path) and filename.
% NOTE: the attenuation data must be the only variable in the file and
% have the dimensions of the final reconstructed image.
options.attenuation_datafile = '';
 
%%%%%%%%%%%%%%%%%%%%%%%% Normalization correction %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute the normalization coefficients
% If set to true, then the normalization coefficients are computed after
% the measurement data has been loaded.
options.compute_normalization = false;

% Normalization correction components to include (1 means that the
% component is included, 0 that it is not included)
% First: Axial geometric correction 
% Second: Detector efficiency correction, use 1 for fan-sum algorithm (both
% sinogram and raw list-mode data) or 2 for SPC (only raw list-mode data)
% Third: Block profile correction
% Fourth: Transaxial geometric correction (NOT recommended when using
% normalization data that does not encompass the entire FOV)
% E.g. [1 1 0 0] computes normalization correction for axial geometric
% effects and detector efficiency, the latter by using fan-sum.
options.normalization_options = [1 1 1 1];

% If a cylinder, that is smaller than the FOV, was used for the normalization
% measurement, specify the radius of this cylinder (cm). Otherwise use an
% empty array or inf.
options.normalization_phantom_radius = inf;

% If the above radius is smaller than the FOV and attenuation has been
% included in the data, then the normalization data can be corrected for
% attenuation. Specify the attenuation coefficient (1/cm) here if you wish
% to include attenuation. Leave empty ([]) if no attenuation should be
% included. If the above radius is inf, this value is ignored.
options.normalization_attenuation = [];

% Apply scatter correction to normalization cylinder
% If cylinder is used for normalization correction, applies also scatter
% correction. Requires the above cylinder radius. 
% NOTE: Applicable only to sinogram data,
options.normalization_scatter_correction = false;
 
%%% Apply normalization correction
% If set to true, normalization correction is applied to either data
% formation or in the image reconstruction by using precomputed 
% normalization coefficients. I.e. once you have computed the normalization
% coefficients, turn above compute_normalization to false and set this to
% true.
options.normalization_correction = true;

%%% Use user-made normalization
% Use either a .mat or .nrm file containing the normalization coefficients
% for normalization correction if normalization_correction is also set to
% true. 
% User will be prompted for the location of the file either during sinogram
% formation or before image reconstruction (see below).
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
% Uses parallel computing toolbox if it is available (parfor).
options.arc_correction = false;

%%% Arc correction interpolation method
% The interpolation method used to interpolate the arc corrected sinogram.
% Available methods are those supported by scatteredInterpolant and
% griddata. If an interpolation method is used which is not supported by
% scatteredInterpolant then griddata will be used instead.
% NOTE: griddata is used if scatteredInterpolant is not found.
% NOTE: GNU Octave supports only griddata.
options.arc_interpolation = 'linear';

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Global corrections %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Global correction factor
% This correction factor will be applied (if nonzero) to all LORs equally.
% This can be e.g. dead time correction factor.
options.global_correction_factor = [];
 
%%%%%%%%%%%%%%%%%%%% Corrections during reconstruction %%%%%%%%%%%%%%%%%%%%
% If set to true, all the corrections are performed during the
% reconstruction step, otherwise the corrections are performed to the
% sinogram/raw data before reconstruction. I.e. this can be considered as
% e.g. normalization weighted reconstruction if normalization correction is
% applied.
% NOTE: Attenuation correction is always performed during reconstruction
% regardless of the choice here.
options.corrections_during_reconstruction = false;
 
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
options.tot_time = inf;

%%% Number of time points/dynamic frames (if a static measurement, use 1)
options.partitions = 1;

%%% Start time (s) (all measurements BEFORE this will be ignored)
options.start = 0;

%%% End time (s) (all measurements AFTER this will be ignored)
% Use inf if you want to the end of the examination (static measurement
% only)
options.end = options.tot_time;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TOF PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Total number of TOF bins
options.TOF_bins = 1;

%%% Length of each TOF bin (s)
% The time length of each TOF bin in seconds
% This multiplied with the number of bins total the entire time frame that
% the TOF data contains. For example with 10 bins of size 400 ps all time
% differences of at most 4 ns will be included in the TOF data. The
% multiplied value should be, at most, the size of the coincidence window.
options.TOF_width = 100e-12;

%%% TOF offset (s)
% If your TOF bins are not centered on zero (center of FOV) you can specify
% the offset value here.
options.TOF_offset = 0;

%%% FWHM of the temporal noise/FWHM of the TOF data (s)
% This parameter has two properties. The first one applies to any TOF data
% that is saved by OMEGA (GATE, Inveon/Biograph list-mode), the second only
% to GATE data.
% Firstly this specifies the FWHM of TOF data used for file naming and
% loading purposes. This value is included in the filename when data is
% imported/saved and also used when that same data is later loaded. 
% Secondly, this is the FWHM of the ADDED temporal noise to the time
% differences. If you are using GATE data and have set a custom temporal
% blurring in GATE then you should set to this zero if you wish to use the
% same temporal resolution. If no custom temporal blurring was applied then
% use this value to control the accuracy of the TOF data. For example if
% you want to have TOF data with 500 ps FWHM then set this value to
% 500e-12. 
options.TOF_noise_FWHM = 200e-12;

%%% FWHM of the TOF data (s)
% Applies to ALL data.
% This value specifies the TOF accuracy during the reconstruction process
% and thus can be different from above. If you are using GATE data with
% temporal blurring, you need to multiply that FWHM with sqrt(2) here.
options.TOF_FWHM = 210e-12;

%%% Number of TOF bins used in reconstruction
% Number of TOF bins used during reconstruction phase.
% NOTE: Currently supports only either all bins specified by
% options.TOF_bins or 1 (or 0) which converts the TOF data into non-TOF
% data during reconstruction phase.
options.TOF_bins_used = options.TOF_bins;
 
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
options.name = 'cylpet_example';

%%% Precompute data 
% This should be done when using data from a certain machine the first time
% as it can speed up reconstruction. Especially recommended for raw
% list-mode data. Not mandatory and the precomputed data is only used if 
% the below precompute_lor is set to true. If using solely implementation 
% 1, this is HIGHLY recommended. 
options.precompute = false;

%%% Use precomputed geometrical matrix information
% During the precompute-phase the number of voxels each LOR traverse is
% counted (this phase requires the above precompute-option to true). These
% are saved and later on used during the reconstruction process. Once the
% precompute-phase has been done once for the specific sinogram/raw data
% and image resolution, it is not necessary to do it again. Recommended for
% raw list-mode data, but speeds up sinogram reconstruction as well. HIGHLY
% recommended when using implementation 1.
options.precompute_lor = false;

%%% Precompute all
% Set this option to true if you want to precompute all possible
% combinations in one go (i.e. precomputed LORs with both sinogram format
% and raw data). Requires for precompute option to be true to take effect.
% Setting this option to true also causes the precompute phase to be
% performed even if only_reconstructions is true (and precompute_lor =
% true).
options.precompute_all = false;

%%% Show status messages
% These are e.g. time elapsed on various functions and what steps have been
% completed. It is recommended to keep this true.
options.verbose = true;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
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
% 3 = Multi-GPU/device matrix-free OpenCL.
% 4 = Matrix-free reconstruction with OpenMP (parallel), standard C++.
% See the wiki for more information: 
% https://github.com/villekf/OMEGA/wiki/Useful-information#selecting-the-correct-implementation
% NOTE: Forward and/or backward projections are ONLY supported with
% implementations 1, 3 and 4.
options.implementation = 3;

% Applies to implementation 3 ONLY
%%% Platform used
% In implementation 3, this determines the platform from where the
% device(s) are taken.
% NOTE: Use OpenCL_device_info() to determine the platform numbers and
% their respective devices.
% NOTE: if you switch platforms then you need to run the below line
% (uncommented) as well:
% clear mex
options.use_device = 0;

% Applies to implementation 3 ONLY
%%% Use 64-bit integer atomic functions
% If true, then 64-bit integer atomic functions (atomic add) will be used
% if they are supported by the selected device.
% Setting this to true will make computations faster on GPUs that support
% the functions, but might make results slightly less reliable due to
% floating point rounding. Recommended for GPUs.
options.use_64bit_atomics = true;

% Implementation 3 ONLY
%%% How many times more measurements/LORs are in the GPU part (applicable if
% heterogeneous computing (CPU + GPU) is used).
% Alternatively, set this to 0 to use only a single device on the specific
% platform (the one with the highest memory count will be used).
options.cpu_to_gpu_factor = 1;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROJECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Type of projector to use for the geometric matrix
% 0 = Regular Siddon's algorithm (only available with implementation 1 and
% when precomputed_lor = false) NOT RECOMMENDED.
% 1 = Improved/accelerated Siddon's algorithm
% 2 = Orthogonal distance based ray tracer
% 3 = Volume of intersection based ray tracer
% See the wiki for more information:
% https://github.com/villekf/OMEGA/wiki/Useful-information#selecting-the-projector
options.projector_type = 1;

%%% Use point spread function (PSF) blurring
% Applies PSF blurring through convolution to the image space. This is the
% same as multiplying the geometric matrix with an image blurring matrix.
options.use_psf = false;

% FWHM of the Gaussian used in PSF blurring in all three dimensions
options.FWHM = [options.cr_p options.cr_p options.cr_pz];

% Use deblurring phase (PSF ONLY)
% If enabled, a deblurring phase is performed once the reconstruction has
% completed. This step is performed for all iterations (deblurred estimates
% are NOT used in the reconstruction phase). This is used ONLY when PSF
% blurring is used.
options.deblurring = true;
% Number of deblurring iterations
% How many iterations of the deblurring step is performed
options.deblur_iterations = 10;

% Orthogonal ray tracer (projector_type = 2) only
%%% The 2D (XY) width of the "strip/tube" where the orthogonal distances are
% included. If tube_width_z below is non-zero, then this value is ignored.
options.tube_width_xy = options.cr_p;

% Orthogonal ray tracer (projector_type = 2) only
%%% The 3D (Z) width of the "tube" where the orthogonal distances are
% included. If set to 0, then the 2D orthogonal ray tracer is used. If this
% value is non-zero then the above value is IGNORED.
options.tube_width_z = options.cr_pz;

% Volume ray tracer (projector_type = 3) only
%%% Radius of the tube-of-response (cylinder)
% The radius of the cylinder that approximates the tube-of-response.
options.tube_radius = sqrt(2) * (options.cr_pz / 2);

% Volume ray tracer (projector_type = 3) only
%%% Relative size of the voxel (sphere)
% In volume ray tracer, the voxels are modeled as spheres. This value
% specifies the relative radius of the sphere such that with 1 the sphere
% is just large enoough to encompass an entire cubic voxel, i.e. the
% corners of the cubic voxel intersect with the sphere shell. Larger values
% create larger spheres, while smaller values create smaller spheres.
options.voxel_radius = 1;

% Siddon (projector_type = 1) only
%%% Number of rays
% Number of rays used per detector if projector_type = 1 (i.e. Improved
% Siddon is used) and precompute_lor = false. I.e. when using precomputed
% LOR data, only 1 rays is always used.
% Number of rays in transaxial direction
options.n_rays_transaxial = 1;
% Number of rays in axial direction
options.n_rays_axial = 1;

% Orthogonal and volume ray tracers (projector_type = 2 and 3) only
% Implementation 3 ONLY
%%% Apply acceleration
% If true, then intermediate results are saved in memory. If you run out
% memory, you should set this to false. Does not apply to improved Siddon.
% Applies only if using implementation 3. Can speed up computations by
% around 30% if set to true.
options.apply_acceleration = true;
 
%%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Number of iterations
options.Niter = 4;

%%% Number of subsets (all excluding MLEM and subset_type = 6)
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
% 360/n_angles for 3D, see GitHub wiki for more information:
% https://github.com/villekf/OMEGA/wiki/Function-help#reconstruction-settings
% 7 = Form the subsets by using golden angle sampling
options.subset_type = 1;

%%% How many angles are combined in subset_type = 6
% E.g. there are 180 angles, in n_angles = 2, then angles 0 and 1 are
% combined to the same subset, 2 and 3, etc.
options.n_angles = 2;

%%% Initial value for the reconstruction
options.x0 = ones(options.Nx, options.Ny, options.Nz);

%%% Epsilon value 
% A small value to prevent division by zero and square root of zero. Should
% not be smaller than eps.
options.epps = 1e-8;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MISC SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Use Shuffle (recommended)
% NOTE: Applies only when using subset_type = 3. 
% Accelerates the subset formation and uses less memory. Not included in
% OMEGA, needs to be manually downloaded and installed.
% Download from: 
% https://se.mathworks.com/matlabcentral/fileexchange/27076-shuffle
options.use_Shuffle = false;

%%% Use fast sparse
% Not included in OMEGA, needs to be manually downloaded and installed.
% Download from: https://github.com/stefanengblom/stenglib
% NOTE: This applies only to implementation 1 when precompute_lor is false.
% Enabling this will automatically cause the function to be called if it is
% found on MATLAB path. 
% NOTE: Suggested only for MATLAB 2019b and earlier.
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% CUSTOM DETECTOR COORDINATES %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load your custom detector coordinates and replace the below examples. In
% the example below x contains the detector coordinates for the x-direction
% and so on. For sinogram data, the dimensions need to be the following:
% size(x) = [(number of elements in a SINGLE sinogram) 2]
% size(y) = [(number of elements in a SINGLE sinogram) 2]
% size(z_det) = [(total number of sinograms) 2]
% E.g. each column represents the detector coordinates of one of the
% detectors in a LOR.
% For raw data:
% size(x) = [(detectors per ring) 1]
% size(y) = [(detectors per ring) 1]
% size(z_det) = [(total number of rings) 1]
% OR
% size(x) = [(total number of measurements) 2]
% size(y) = [(total number of measurements) 2]
% size(z_det) = [(total number of measurements) 2]
% For raw data, as long as it is formatted as in OMEGA, you only need each
% detector coordinate once.

% options.x = x;
% options.y = y;
% options.z_det = z;

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
% options.DOI = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 



%% Precompute the necessary data

if options.precompute
    precompute_data(options);
end

%% Class example (OSEM/MLEM)

% Here is an example of how to obtain the same results as above by using a
% specific MATLAB class. This is a bit more simplified from above and also
% allows more easily to use other properties files (such as
% Inveon_PET_main.m). PSF blurring will be performed automatically if it
% has been selected.

% Construct the forward and backward projections object (you need to rerun
% this if you make any changes to the system): 
A = forwardBackwardProject(options);
% Load the measurement data
load('Cylindrical_PET_example_cylpet_example_sinograms_combined_static_200x168x703_span3.mat','raw_SinM')
if options.implementation == 1
    raw_SinM = double(raw_SinM(A.index));
else
    raw_SinM = single(raw_SinM(A.index));
end

% Compute the OSEM-estimate for the specified number of iterations and
% subsets (use 1 subset if you want to use a method that doesn't use
% subsets):
f = options.x0(:);
for iter = 1 : options.Niter
    for osa_iter = 1 : options.subsets
        % The result is stored in y
        y = forwardProject(A, f, osa_iter);
        % The result is stored in x
        if iter == 1
            % Sensitivity image can be computed during the first iteration
            % Computed ONLY if the second output is present
            [x, A] = backwardProject(A, raw_SinM(A.nn(osa_iter) + 1:A.nn(osa_iter+1)) ./ (y + options.epps), osa_iter);
        else
            x = backwardProject(A, raw_SinM(A.nn(osa_iter) + 1:A.nn(osa_iter+1)) ./ (y + options.epps), osa_iter);
        end
        f = (f ./ (A.sens(:,osa_iter) + options.epps)) .* (x + options.epps);
    end
end
% PSF deblurring phase
if options.use_psf && options.deblurring
    f = deblur(f, options, gaussK, options.Nx, options.Ny, options.Nz);
end
ff = reshape(f, options.Nx,options.Ny,options.Nz);
f_osem = ff;


%% Conjugate-gradient Least-squares (CGLS) example

% This example demonstrates the "operator overloading" available in the
% forwardBackwardProject class. This means that both the forwardProject and
% backwardProject can be replaced with simple multiplication *. The
% dimension of the input vector determines whether forward or backward
% projection is used.

% Subsets should be set to 1 if no subsets are used
options.subsets = 1;

% Create the class object
A = forwardBackwardProject(options);
% Load measurement data
load('Cylindrical_PET_example_cylpet_example_sinograms_combined_static_200x168x703_span3.mat','raw_SinM')
if options.implementation == 1
    raw_SinM = double(raw_SinM(A.index));
else
    raw_SinM = single(raw_SinM(A.index));
end
f = options.x0(:);
r = raw_SinM;
s = A * raw_SinM;
p = s;
gamma = norm(s)^2;
for iter = 1 : options.Niter
    q = A * p;
    alpha = gamma / norm(q)^2;
    f = f + alpha * p;
    r = r - alpha * q;
    s = A * r; % This is equivalent to A' * r
    gamma_ = norm(s)^2;
    beta = gamma_ / gamma;
    p = s + beta * p;
    gamma = gamma_;
end
% PSF deblurring phase
if options.use_psf && options.deblurring
    f = deblur(f, options, gaussK, options.Nx, options.Ny, options.Nz);
end
ff = reshape(f, options.Nx,options.Ny,options.Nz);

f_osem = ff;

%% LSQR example

options.subsets = 1;

A = forwardBackwardProject(options);
load('Cylindrical_PET_example_cylpet_example_sinograms_combined_static_200x168x703_span3.mat','raw_SinM')
if options.implementation == 1
    raw_SinM = double(raw_SinM(A.index));
else
    raw_SinM = single(raw_SinM(A.index));
end
f = zeros(options.Nx*options.Ny*options.Nz,1);
u = raw_SinM;
beta = 1 / norm(u);
u = u * beta;
vh = backwardProject(A, u);
alpha = 1 / norm(vh);
vh = vh * alpha;
w = vh;
phi = beta;
rho = alpha;
for iter = 1 : options.Niter
    u = forwardProject(A, vh) - alpha * u;
    beta = 1 / norm(u);
    u = u * beta;
    vh = backwardProject(A, u) - beta * vh;
    alpha = 1 / norm(vh);
    vh = vh * alpha;
    rho_ = sqrt(rho^2 + beta^2);
    c = rho / rho_;
    s = beta / rho_;
    theta = s * alpha;
    rho = -c * alpha;
    phi_ = c * phi;
    phi = s* phi;
    f = f + (phi_ / rho_)*w;
    w = vh + (theta/rho_)*w;
end
% PSF deblurring phase
if options.use_psf && options.deblurring
    f = deblur(f, options, gaussK, options.Nx, options.Ny, options.Nz);
end
ff = reshape(f, options.Nx,options.Ny,options.Nz);

f_osem = ff;

%% SA-WLS example (includes TOF example)

options.subsets = 1;

A = forwardBackwardProject(options);
load('Cylindrical_PET_example_cylpet_example_sinograms_combined_static_200x168x703_span3.mat','raw_SinM')
if options.TOF_bins > 1
    raw_SinM = reshape(raw_SinM, numel(raw_SinM) / options.TOF_bins, options.TOF_bins);
else
    raw_SinM = raw_SinM(:);
end
if options.implementation == 1
    raw_SinM = double(raw_SinM(A.index,:));
else
    raw_SinM = single(raw_SinM(A.index,:));
end
% Sensitivity image
[f, A] = backwardProject(A, raw_SinM);
y = A * f;
for iter = 1 : options.Niter
    % Brackets are needed here
    % It is not necessary to transpose the matrix when doing backprojection
    f_ = f .* sqrt((1./A.sens) .* (A * (raw_SinM.^2 ./ (y.^2 + 1e-4))));
    y = y + (A * (f_ - f));
    f = f_;
end

ff = reshape(f, options.Nx,options.Ny,options.Nz);

f_osem = ff;

%% System matrix (OSEM) example

% This example demonstrates the use of the forwardBackwardProject class to
% obtain the system matrix itself (either the full matrix or a subset). If
% you are using precomputed data (i.e. options.precompute_lor = true), then
% the system matrix will be the TRANSPOSE of the system matrix. This
% example assumes that options.precompute_lor = true, which is the
% recommended method for system matrix creation.
% NOTE: Unlike with the forward and backward projections, the PSF blurring
% has to be performed manually here.
options.subsets = 16;

A = forwardBackwardProject(options);
load('Cylindrical_PET_example_cylpet_example_sinograms_combined_static_200x168x703_span3.mat','raw_SinM')
if options.implementation == 1
    raw_SinM = double(raw_SinM(A.index));
else
    raw_SinM = single(raw_SinM(A.index));
end
f = options.x0(:);
% Form the full system matrix with this:
% sysMat = formMatrix(A);
for iter = 1 : options.Niter
    for osa_iter = 1 : options.subsets
        % Form a subset of the system matrix with this:
        sysMat = formMatrix(A, osa_iter);
        
        if options.use_psf
            f = computeConvolution(f, options, options.Nx, options.Ny, options.Nz, gaussK);
        end
        
        y = sysMat' * f;
        Summ = full(sum(sysMat,2));
        if options.use_psf
            Summ = computeConvolution(Summ, options, options.Nx, options.Ny, options.Nz, gaussK);
        end
        x = sysMat * (raw_SinM ./ y);
        if options.use_psf
            x = computeConvolution(x, options, options.Nx, options.Ny, options.Nz, gaussK);
        end
        f = (f ./ (Summ + options.epps)) .* (x + options.epps);
    end
end
ff = reshape(f, options.Nx,options.Ny,options.Nz);

f_osem = ff;

%% ATP-WLS example, using * to compute forward projection and '* backward projection

options.subsets = 1;

A = forwardBackwardProject(options);
load('Cylindrical_PET_example_cylpet_example_sinograms_combined_static_200x168x703_span3.mat','raw_SinM')
if options.implementation == 1
    raw_SinM = double(raw_SinM(A.index));
else
    raw_SinM = single(raw_SinM(A.index));
end
f = options.x0(:);
NN = sum(raw_SinM(:));
tau = 0.1 / NN;
for iter = 1 : options.Niter
    % The backprojection can be alternatively computed with A' * f
    f = f .* ((A' * (raw_SinM.^2 ./ (A * f).^2) + 2 * tau * NN) ./ (1 + 2 * tau * sum(f)));
end

ff = reshape(f, options.Nx,options.Ny,options.Nz);

f_osem = ff;


%% EM-PKMA with TV regularization
% This example demonstrates the use of included priors, in this case the TV
% prior

options.subsets = 1;

A = forwardBackwardProject(options);
load('Cylindrical_PET_example_cylpet_example_sinograms_combined_static_200x168x703_span3.mat','raw_SinM')
if options.implementation == 1
    raw_SinM = double(raw_SinM(A.index));
else
    raw_SinM = single(raw_SinM(A.index));
end
[f, A] = backwardProject(A, raw_SinM);
A.sens(A.sens == 0) = 1;
beta = 0.5;
lambda = 0.2;
N = options.Nx*options.Ny*options.Nz;
delta = 0.1;
rho = 0.5;
% Such prepass functions exist for all the priors that need one, see the
% GitHub wiki or the included documentation for more information (Computing
% the forward and or backward projections)
options = TVPrepass(options);


for iter = 1 : options.Niter
    S = (f + options.epps) ./ A.sens;
    f_ = f;
    % Compute the gradient of the prior
    grad = TVpriorFinal(f,options.TVdata,options.Nx, options.Ny, options.Nz, false, options, 1);
    ff = A' * (ones(length(raw_SinM),1,'single') - raw_SinM ./ (A * f));
    f = f - beta .* S .* (ff + lambda * grad);
    f(f < 0) = 0;
    alpha = 1 + (rho * iter) / (iter + delta);
    f = (1 - alpha) .* f_ + alpha .* f;
end

ff = reshape(f, options.Nx,options.Ny,options.Nz);

f_osem = ff;

%% EM-PKMA with NLM regularization
% This example demonstrates the use of included priors, in this case the NLM
% prior

options.subsets = 1;

A = forwardBackwardProject(options);
load('Cylindrical_PET_example_cylpet_example_sinograms_combined_static_200x168x703_span3.mat','raw_SinM')
if options.implementation == 1
    raw_SinM = double(raw_SinM(A.index));
else
    raw_SinM = single(raw_SinM(A.index));
end
[f, A] = backwardProject(A, raw_SinM);
A.sens(A.sens == 0) = 1;
beta = 0.5;
lambda = 0.2;
N = options.Nx*options.Ny*options.Nz;
delta = 0.1;
rho = 0.5;
% NLM parameters need to be specified
options.Ndx = 2;
options.Ndy = 2;
options.Ndz = 2;
options.sigma = 10;
options.Nlx = 2;
options.Nly = 2;
options.Nlz = 2;
options.NLM_gauss = 1;
options.NLM_use_anatomical = false;
options.NLM_reference_image = 'reference_image.mat';
options.NLTV = false;
options.NLM_MRP = false;
% Such prepass functions exist for all the priors that need one, see the
% GitHub wiki or the included documentation for more information (Computing
% the forward and or backward projections)
options = NLMPrepass(options);


for iter = 1 : options.Niter
    S = (f + options.epps) ./ A.sens;
    f_ = f;
    % Compute the gradient of the prior
    grad = NLM(f, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
        options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
    ff = A' * (ones(length(raw_SinM),1,'single') - raw_SinM ./ (A * f));
    f = f - beta .* S .* (ff + lambda * grad);
    f(f < 0) = 0;
    alpha = 1 + (rho * iter) / (iter + delta);
    f = (1 - alpha) .* f_ + alpha .* f;
end

ff = reshape(f, options.Nx,options.Ny,options.Nz);

f_osem = ff;

%% ADMM with TV

options.subsets = 1;

A = forwardBackwardProject(options);
load('Cylindrical_PET_example_cylpet_example_sinograms_combined_static_200x168x703_span3.mat','raw_SinM')
if options.implementation == 1
    raw_SinM = double(raw_SinM(A.index));
else
    raw_SinM = single(raw_SinM(A.index));
end
% Compute the sensitivity image A.sens
[f, A] = backwardProject(A, raw_SinM);
options = TVPrepass(options);
beta = 4;
mu = 0.001;
d = f;
u = f;

for iter = 1 : 20
    grad = TVpriorFinal(u,options.TVdata,options.Nx, options.Ny, options.Nz, false, options, 1);
    gradPhi = mu .* (u - f + d) + beta * grad;
    grad2 = TVpriorFinal(gradPhi,options.TVdata,options.Nx, options.Ny, options.Nz, false, options, 1);
    gradPhi2 = mu * gradPhi + grad2;
    alpha = norm(gradPhi)^2 / (gradPhi' * gradPhi2);
    u = u - alpha * gradPhi;
    w = A.sens - mu * (u - d);
    v = f .* (A * (raw_SinM ./ (A * f)));
    f = (-w + sqrt(w.^2 + 4*mu * v))/ (2*mu);
    d = d + f - u;
end

ff = reshape(f, options.Nx,options.Ny,options.Nz);

f_osem = ff;

%% Forward and backward projections (and sensitivity image)
% This example uses the older, function-based, way of computing the forward
% and backward projections

% Get the subset indices for the current subset type
% Index contains the subset indices (which measurements belong to each
% subset)
% n_meas is the number of measurements in each subset
[index, n_meas, options.subsets] = index_maker(options.Nx, options.Ny, options.Nz, options.subsets, options.use_raw_data, ...
    options.machine_name, options, options.Nang, options.Ndist, options.TotSinos, options.NSinos);

nn = [0;cumsum(n_meas)];
if iscell(index)
    index = cell2mat(index);
end

[gaussK, options] = PSFKernel(options);

f = options.x0(:);

% Load your own data here
% Either use the same raw_SinM variable name or rename all instances of
% raw_SinM
load('Cylindrical_PET_example_cylpet_example_sinograms_combined_static_200x168x703_span3.mat','raw_SinM')
if options.implementation == 1
    raw_SinM = double(raw_SinM(index));
else
    raw_SinM = single(raw_SinM(index));
end

sens = zeros(size(f,1),options.subsets);
% Computes the OSEM-estimates by using the forward and backprojections (not
% recommended for an actual OSEM computation due to reduced performance)
for iter = 1 : options.Niter
    for osa_iter = 1 : options.subsets
        % Compute PSF blurring
        if options.use_psf
            f = computeConvolution(f, options, options.Nx, options.Ny, options.Nz, gaussK);
        end
        % Compute the forward projection (y = A*f)
        % Only those indices that have non-zero measurements (counts) are
        % included, i.e. the "matrix-vector multiplication" A*f is
        % performed ONLY for the indices where raw_SinM is non-zero.
        % To include ALL indices replace raw_SinM(nn(osa_iter) + 1:nn(osa_iter+1))
        % with a vector of ones, e.g. as in the commented function call
        % below. The measurements themselves are not used in the
        % computation of forward projection.
        [y, options] = forward_project(options, index(nn(osa_iter) + 1:nn(osa_iter+1)), n_meas(osa_iter), f, [nn(osa_iter) + 1 , nn(osa_iter+1)], ...
            iter + osa_iter - 1);
        %         [y, options] = forward_project(options, index(nn(osa_iter) + 1:nn(osa_iter+1)), n_meas(osa_iter), f, [nn(osa_iter) + 1 , nn(osa_iter+1)], ...
        %             iter + osa_iter - 1, ones(length(nn(osa_iter) + 1:nn(osa_iter+1)),1,'single'));
        
        % Compute the backprojection (x = A'*b, where in this case b =
        % raw_SinM(nn(osa_iter) + 1:nn(osa_iter+1)) ./ (y + options.epps))
        if iter == 1
            % Compute the sensitivity image on the first iteration.
            [x, norm] = backproject(options, index(nn(osa_iter) + 1:nn(osa_iter+1)), n_meas(osa_iter), raw_SinM(nn(osa_iter) + 1:nn(osa_iter+1)) ./ (y + options.epps), ...
                [nn(osa_iter) + 1,nn(osa_iter+1)], iter + osa_iter - 1, raw_SinM(nn(osa_iter) + 1:nn(osa_iter+1)));
            sens(:,osa_iter) = norm;
            if options.use_psf
                sens(:,osa_iter) = computeConvolution(sens(:,osa_iter), options, options.Nx, options.Ny, options.Nz, gaussK);
            end
            % Compute only x. E.g. x = backproject(...) to compute only the
            % backprojection without the sensitivity image.
            %             x = backproject(options, index(nn(osa_iter) + 1:nn(osa_iter+1)), n_meas(osa_iter), raw_SinM(nn(osa_iter) + 1:nn(osa_iter+1)) ./ (y + options.epps), ...
            %                 nn(osa_iter) + 1:nn(osa_iter+1), iter + osa_iter - 1, raw_SinM(nn(osa_iter) + 1:nn(osa_iter+1)));
        else
            [x] = backproject(options, index(nn(osa_iter) + 1:nn(osa_iter+1)), n_meas(osa_iter), raw_SinM(nn(osa_iter) + 1:nn(osa_iter+1)) ./ (y + options.epps), ...
                [nn(osa_iter) + 1,nn(osa_iter+1)], iter + osa_iter - 1, raw_SinM(nn(osa_iter) + 1:nn(osa_iter+1)));
        end
        % Compute PSF blurring
        if options.use_psf
            x = computeConvolution(x, options, options.Nx, options.Ny, options.Nz, gaussK);
        end
        f = (f ./ (sens(:,osa_iter) + options.epps)) .* (x + options.epps);
    end
end
% PSF deblurring phase
if options.use_psf && options.deblurring
    f = deblur(f, options, gaussK, options.Nx, options.Ny, options.Nz);
end
ff = reshape(f, options.Nx,options.Ny,options.Nz);

f_osem = ff;