%% MATLAB codes for generic PET reconstruction using Inveon list-mode data
% This example show cases a "generic" PET example. The example data is the
% Inveon preclinical list-mode PET data
% (https://zenodo.org/records/3528056), but any sinogram PET data can be
% used. https://doi.org/10.5281/zenodo.4646897 works also
% Note that https://doi.org/10.5281/zenodo.4646897 has better attenuation
% correction

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
% Multiplying this with the below cryst_per_block should equal the total
% number of crystal rings. 
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
options.cr_p = 1.63;

%%% Crystal pitch/size in z-direction (axial) (mm)
options.cr_pz = 1.592;

%%% Ring diameter (distance between perpendicular detectors) (mm)
options.diameter = 161.08;

% Note that non-square transaxial FOV sizes should work, but might not work
% always. Square transaxial FOV is thus recommended.
%%% Transaxial FOV size (mm), this is the length of the x (vertical) side
% of the FOV
options.FOVa_x = 100;

%%% Transaxial FOV size (mm), this is the length of the y (horizontal) side
% of the FOV
options.FOVa_y = options.FOVa_x;

% The above recommendation doesn't apply to axial FOV, i.e. this can be
% different from the transaxial FOV size(s). 
%%% Axial FOV (mm)
options.axial_fov = 127;

%%% Number of pseudo rings between physical rings (use 0 or [] if none)
% NOTE: Inveon has no pseudo detectors/rings
options.pseudot = [];

%%% Number of gaps between rings
% These should be the rings after which there is a gap
% For example, if there are a total of 6 rings with two gaps and the gaps
% are after ring number 2 and 4. Then options.ringGaps = [2,4];
% Note that the below options.rings should not include the gaps, though by
% default it should be correctly computed.
options.ringGaps = [];

%%% Number of detectors per ring (without pseudo detectors)
options.det_per_ring = options.blocks_per_ring*options.cryst_per_block;

%%% Number of detectors per ring (with pseudo detectors)
% NOTE: Inveon has no pseudo detectors/rings
% If you have a single pseudo detector per block, use the commented line
options.det_w_pseudo = options.blocks_per_ring*(options.cryst_per_block);
% options.det_w_pseudo = options.blocks_per_ring*(options.cryst_per_block + 1);

%%% Number of crystal rings
options.rings = options.linear_multip * options.cryst_per_block;

%%% Number of detectors
options.detectors = options.det_per_ring*options.rings;

%%% Scanner name
% Used for naming purposes (measurement data)
options.machine_name = 'Inveon';


% Note: Origin is assumed to be at the center. If this is not the case, you
% can shift it with options.oOffsetX, options.oOffsetY and options.oOffsetZ
% options.oOffsetX = 0;
% options.oOffsetY = 0;
% options.oOffsetZ = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INVEON DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note that if you use your own data, this step is not required. Instead,
% you should load your sinogram data directly into options.SinM variable.
% See sinogram properties below for the correct orientation
%%% Use Inveon data
% 0/false = Don't use Inveon list-mode data
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
options.Nx = 256;

%%% Y/column-direction
options.Ny = 256;

% The above, again, doesn't apply to axial direction
%%% Z-direction (number of slices) (axial)
options.Nz = options.rings*2 - 1;

%%% Flip the image (in horizontal direction)?
options.flip_image = false;

%%% How much is the image rotated?
% NOTE: The rotation is done in the detector space (before reconstruction).
% This current setting is for scanner list-mode data or sinogram data.
% Positive values perform the rotation in clockwise direction
options.offangle = options.det_w_pseudo * (2/4) - options.cryst_per_block/2;
 
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
% This is the ROW dimension, i.e. the number of rows in the sinogram
options.Ndist = 128;

%%% Number of angles (tangential positions) in sinogram 
% This is the final amount after possible mashing, maximum allowed is the
% number of detectors per ring/2.
% This is the COLUMN dimension, i.e. the number of clumns in the sinogram
options.Nang = 160;

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
% This is the number of slices in the sinogram
options.NSinos = options.TotSinos;

%%% If Ndist value is even, take one extra out of the negative side (+1) or
% from the positive side (-1). E.g. if Ndist = 200, then with +1 the
% interval is [-99,100] and with -1 [-100,99]. This varies from device to
% device. If you see a slight shift in the sinograms when comparing with
% the scanner sinograms then use the other option here. For Inveon, this
% should be -1.
options.ndist_side = -1;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CORRECTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that generally precorrection is not available with custom PET data,
% only in the built-in cases such as the Inveon PET data used in this
% example. Thus, when using your own PET data, it is recommended to have
% options.corrections_during_reconstruction = true if you want to use
% corrections. 
% NOTE: When using custom data, the additional corrections below, such as
% variance reduction or smoothing, might not work properly!
 
%%%%%%%%%%%%%%%%%%%%%%%%%%% Randoms correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If true, uses the input randoms correction data. The randoms need to
% input into options.SinDelayed. If you use the Inveon data of this
% example, setting this to true will automatically load the data.
options.randoms_correction = true;

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
% If set to true, will prompt the user to load the scatter sinogram.
% Corrects for scatter during data formation/load or during
% reconstruction.
% NOTE: Scatter data is not created by this software and as such must be
% provided by the user. Previously created scatter sinogram/raw data matrix
% can be used though. Scatter sinograms created by Inveon AW are supported.
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


%%%%%%%%%%%%%%%%%%%%%%%%% Attenuation correction %%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Image-based attenuation correction
% Include attenuation correction from images (e.g. CT-images) (for this you
% need attenuation images of each slice correctly rotated and scaled for
% 511 keV). For CT-images you can use attenuationCT_to_511 or
% create_atten_matrix_CT functions. If the below attenuation_datafile is
% empty, UMAPs created with Inveon AW can be used (if the below
% CT_attenuation is set to true) or, alternatively, use the produced .atn
% files instead (CT_attenuation is false). 
% E.g. if this is set to true, CT_attenuation = true and
% attenuation_datafile = '', then the user will be automatically prompted
% for the UMAP-files and they will be automatically saved as a mat-file
% with the filename saved in the attenuation_datafile field. Alternatively,
% if this is set to true, CT_attenuation = false and attenuation_datafile
% = '', then the user will be automatically prompted for the .atn-files
% from which the attenuation images will be automatically created and
% saved (filename is saved in the attenuation_datafile field).
options.attenuation_correction = false;

%%% CT-image attenuation
% Use CT-images (UMAP-image) for the attenuation. If set to false, uses the
% .atn-files instead (if above attenuation is set to true). 
options.CT_attenuation = false;

%%% Attenuation image data file
% Specify the path (if not in MATLAB path) and filename.
% NOTE: the attenuation data must be the only variable in the file and
% have the dimensions of the final reconstructed image. Previously
% saved attenuation images can be used here.
options.attenuation_datafile = '';

%%% Rotate the attenuation image before correction
% Rotates the attenuation image N * 90 degrees where N is the number
% specified below. Positive values are clockwise, negative
% counter-clockwise.
options.rotateAttImage = 0;
 

%%%%%%%%%%%%%%%%%%%%%%%% Normalization correction %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute the normalization coefficients
% If set to true, then the normalization coefficients are computed after
% the measurement data has been loaded.
% You can use this to compute normalization coefficients for the selected
% scanner, assuming that the input data is measurement data fit for
% normalization correction.
options.compute_normalization = false;

% Normalization correction components to include (1 means that the
% component is included, 0 that it is not included)
% First: Axial geometric correction 
% Second: Detector efficiency correction
% Third: Block profile correction
% Fourth: Transaxial geometric correction (NOT recommended when using
% normalization data that does not encompass the entire FOV)
% E.g. [1 1 0 0] computes normalization correction for axial geometric
% effects and detector efficiency.
options.normalization_options = [1 1 1 1];

% If a cylinder, that is smaller than the FOV, was used for the normalization
% measurement, specify the radius of this cylinder (cm). Otherwise use an
% empty array or inf. Inveon uses a small cylinder for normalization
% correction (default radius is 3 cm).
options.normalization_phantom_radius = 3;

% If the above radius is smaller than the FOV and attenuation has been
% included in the data, then the normalization data can be corrected for
% attenuation. Specify the attenuation coefficient (1/cm) here if you wish
% to include attenuation. Leave empty ([]) if no attenuation should be
% included. If the above radius is inf, this value is ignored.
options.normalization_attenuation = 0.4371;

% Apply scatter correction to normalization cylinder
% If cylinder is used for normalization correction, applies also scatter
% correction. Requires the above cylinder radius. 
% NOTE: Applicable only to sinogram data.
options.normalization_scatter_correction = true;
 
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
% coefficients for the specified scanner will be automatically loaded. Use
% this only if you want to use normalization coefficients computed outside
% of OMEGA.
% NOTE: Supports .nrm files created by the Inveon AW. 
options.use_user_normalization = true;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Arc correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Apply arc correction
% NOTE: Arc correction is an experimental feature. It is currently
% relatively slow. Generally it is not recommended to use arc correction
% (Inveon data is an exception). Uses parallel computing toolbox if it is
% available (parfor). 
% NOTE: For Inveon data, arc correction reduces aliasing artifacts when
% using improved Siddon without PSF. 
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
% applied or ordinary Poisson reconstruction.
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
options.name = 'open_PET_data';
% options.name = 'NEMA_S6_data';

%%% Folder for the data file(s)
% If no files are located in the path provided below, then the user will
% prompted for the data .lst or .scn for Inveon, .mat, .npy or .npz
% otherwise.
if ispc % Windows
    options.fpath = 'C:\path\to\GATE\output\';
else % Unix/Mac
    options.fpath = '/path/to/GATE/output/';
end

%%% Form only sinograms (no reconstructions)
% If this is set to true, running this file will only produce the
% measurement data matrices (sinograms and raw data). Also computes the
% normalization coefficients if they have been selected.
options.only_sinos = false;

%%% Do not perform data load/sinogram creation
% Setting this to true will always skip both the data load and sinogram
% creation when running this file. Sinogram corrections are applied if
% selected. Overwrites the above only_sinos variable. Normalization
% correction coefficients are computed if selected.
options.no_data_load = false;

%%% Compute only the reconstructions
% If this file is run with this set to true, then the data load and
% sinogram formation steps are always skipped. Precomputation step is
% only performed if precompute_lor = true and precompute_all = true
% (below). Normalization coefficients are not computed even if selected.
options.only_reconstructions = false;

%%% Show status messages
% These are e.g. time elapsed on various functions and what steps have been
% completed. It is recommended to keep this 1. This can be at most 3.
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
% 2 = Matrix-free reconstruction with OpenCL/ArrayFire (Recommended)
% (Requires ArrayFire. Compiles with MinGW ONLY when ArrayFire was compiled
% with MinGW as well (cannot use the prebuilt binaries)).
% 3 = Multi-GPU/device matrix-free OpenCL (OSEM & MLEM only).
% 4 = Matrix-free reconstruction with OpenMP (parallel), standard C++
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
% NOTE: if you switch devices then you need to run the below line
% (uncommented) as well:
% clear mex
options.use_device = 0;

% Applies to implementations 2, 3 and 5 ONLY
%%% Use 64-bit integer atomic functions
% If true, then 64-bit integer atomic functions (atomic add) will be used
% if they are supported by the selected device.
% Setting this to true will make computations faster on GPUs that support
% the functions, but might make results slightly less reliable due to
% floating point rounding. Recommended for GPUs.
options.use_64bit_atomics = true;

% Implementation 2 ONLY
%%% Use CUDA
% Selecting this to true will use CUDA kernels/code instead of OpenCL. This
% only works if the CUDA code was successfully built. Recommended only for
% Siddon as the orthogonal/volume-based ray tracer are slower in CUDA.
options.use_CUDA = false;

% Implementation 2 ONLY
%%% Use CPU
% Selecting this to true will use CPU-based code instead of OpenCL or CUDA.
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
% See the documentation for more information:
% https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
options.projector_type = 1;

%%% Use mask
% The mask needs to be a binary mask (uint8 or logical) where 1 means that
% the pixel is included while 0 means it is skipped. Separate masks can be
% used for both forward and backward projection and either one or both can
% be utilized at the same time. E.g. if only backprojection mask is input,
% then only the voxels which have 1 in the mask are reconstructed.
% Currently the masks need to be a 2D image that is applied identically at
% each slice.
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
% value will be multiplied by the voxel size which means that 1 means that
% the interpolation length corresponds to a single voxel (transaxial)
% length. Larger values lead to faster computation but at the cost of
% accuracy. Recommended values are between [0.5 1].
options.dL = 0.5;

%%% Use point spread function (PSF) blurring
% Applies PSF blurring through convolution to the image space. This is the
% same as multiplying the geometric matrix with an image blurring matrix.
options.use_psf = true;

% FWHM of the Gaussian used in PSF blurring in all three dimensions
options.FWHM = [options.cr_p options.cr_p options.cr_pz];

% Orthogonal ray tracer (projector_type = 2 only)
%%% The 2D (XY) width of the "strip/tube" where the orthogonal distances are
% included. If tube_width_z below is non-zero, then this value is ignored.
options.tube_width_xy = options.cr_p;

% Orthogonal ray tracer (projector_type = 2 only)
%%% The 3D (Z) width of the "tube" where the orthogonal distances are
% included. If set to 0, then the 2D orthogonal ray tracer is used. If this
% value is non-zero then the above value is IGNORED.
options.tube_width_z = options.cr_pz;

% Volume ray tracer (projector_type = 3 only)
%%% Radius of the tube-of-response (cylinder)
% The radius of the cylinder that approximates the tube-of-response.
options.tube_radius = sqrt(2) * (options.cr_pz / 2);

% Volume ray tracer (projector_type = 3 only)
%%% Relative size of the voxel (sphere)
% In volume ray tracer, the voxels are modeled as spheres. This value
% specifies the relative radius of the sphere such that with 1 the sphere
% is just large enoough to encompass an entire cubic voxel, i.e. the
% corners of the cubic voxel intersect with the sphere shell. Larger values
% create larger spheres, while smaller values create smaller spheres.
options.voxel_radius = 1;

% Siddon (projector_type = 1 only)
%%% Number of rays
% Number of rays used per detector if projector_type = 1 (i.e. Improved
% Siddon is used) and precompute_lor = false. I.e. when using precomputed
% LOR data, only 1 rays is always used.
% Number of rays in transaxial direction
options.n_rays_transaxial = 1;
% Number of rays in axial direction
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
% 360/n_angles for 3D, see docs for more information:
% https://omega-doc.readthedocs.io/en/latest/algorithms.html#type-6
% 7 = Form the subsets by using golden angle sampling
% 8 = Use every nth sinogram
% 9 = Randomly select the full sinograms
% 10 = Use golden angle sampling to select the subsets (not recommended for
% PET)
% 11 = Use prime factor sampling to select the full sinograms
% Most of the time subset_type 1 or 4 is sufficient.
options.subset_type = 1;

%%% How many angles are combined in subset_type = 6
% E.g. there are 180 angles, in n_angles = 2, then angles 0 and 1 are
% combined to the same subset, 2 and 3, etc.
options.n_angles = 90;

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
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ML-METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ordered Subsets Expectation Maximization (OSEM) OR Maximum-Likelihood
%%% Expectation Maximization (MLEM) (if subsets = 1)
% Supported by all implementations
options.OSEM = true;

%%% Modified Row-Action Maximum Likelihood Algorithm (MRAMLA)
% Supported by implementations 1, 2, 4, and 5
options.MRAMLA = false;

%%% Row-Action Maximum Likelihood Algorithm (RAMLA)
% Supported by implementations 1, 2, 4, and 5
options.RAMLA = false;

%%% Relaxed Ordered Subsets Expectation Maximization (ROSEM)
% Supported by implementations 1, 2, 4, and 5
options.ROSEM = false;

%%% Rescaled Block Iterative Expectation Maximization (RBI-EM)
% Supported by implementations 1, 2, 4, and 5
options.RBI = false;

%%% Dynamic RAMLA (DRAMA)
% Supported by implementations 1, 2, 4, and 5
options.DRAMA = false;

%%% Complete data OSEM (COSEM)
% Supported by implementations 1, 2, 4, and 5
options.COSEM = false;

%%% Enhanced COSEM (ECOSEM)
% Supported by implementations 1, 2, 4, and 5
options.ECOSEM = false;

%%% Accelerated COSEM (ACOSEM)
% Supported by implementations 1, 2, 4, and 5
options.ACOSEM = false;

%%% FISTA
% Supported by implementations 1, 2, 4, and 5
options.FISTA = false;

%%% FISTA with L1 regularization (FISTAL1)
% Supported by implementations 1, 2, 4, and 5
options.FISTAL1 = false;

%%% LSQR
% Supported by implementations 1, 2, 4, and 5
options.LSQR = false;

%%% CGLS
% Supported by implementations 1, 2, 4, and 5
options.CGLS = false;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAP-METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Any algorithm selected here will utilize any of the priors selected below
% this. Note that only one algorithm and prior combination is allowed! You
% can also use most of these algorithms without priors (such as PKMA or
% PDHG).
%%% One-Step Late MLEM (OSL-MLEM)
% Supported by implementations 1, 2, 4, and 5
options.OSL_MLEM = false;

%%% One-Step Late OSEM (OSL-OSEM)
% Supported by implementations 1, 2, 4, and 5
options.OSL_OSEM = false;

%%% Modified BSREM (MBSREM)
% Supported by implementations 1, 2, 4, and 5
options.MBSREM = false;

%%% Block Sequential Regularized Expectation Maximization (BSREM)
% Supported by implementations 1, 2, 4, and 5
options.BSREM = false;

%%% ROSEM-MAP
% Supported by implementations 1, 2, 4, and 5
options.ROSEM_MAP = false;

%%% RBI-OSL
% Supported by implementations 1, 2, 4, and 5
options.OSL_RBI = false;

%%% (A)COSEM-OSL
% 0/false = No COSEM-OSL, 1/true = ACOSEM-OSL, 2 = COSEM-OSL
% Supported by implementations 1, 2, 4, and 5
options.OSL_COSEM = false;

%%% Preconditioner Krasnoselskii-Mann algorithm (PKMA)
% Supported by implementations 1, 2, 4, and 5
options.PKMA = false;

%%% Primal-dual hybrid gradient (PDHG)
% Supported by implementations 1, 2, 4, and 5
options.PDHG = false;

%%% Primal-dual hybrid gradient (PDHG) with L1 minimization
% Supported by implementations 1, 2, 4, and 5
options.PDHGL1 = false;

%%% Primal-dual hybrid gradient (PDHG) with KUllback-Leibler minimization
% Supported by implementations 1, 2, 4, and 5
options.PDHGKL = false;

%%% Primal-dual Davis-Yin (PDDY)
% Supported by implementation 2
options.PDDY = false;


 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Median Root Prior (MRP)
options.MRP = false;

%%% Quadratic Prior (QP)
options.quad = false;

%%% Huber Prior (QP)
options.Huber = false;

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

%%% Hyperbolic prior
options.hyperbolic = false;

%%% Total Generalized Variation (TGV) prior
options.TGV = false;

%%% Non-local Means (NLM) prior
options.NLM = false;

%%% Relative difference prior
options.RDP = false;

%%% Generalized Gaussian Markov random field (GGMRF) prior
options.GGMRF = false;


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
% with lamda = 1. Note that current_iteration is one-based, i.e. it starts
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

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% DRAMA PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Beta_0 value
options.beta0_drama = 0.1;
%%% Beta value
options.beta_drama = 1;
%%% Alpha value
options.alpha_drama = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PDHG PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Primal value
% If left zero, or empty, it will be automatically computed
options.tauCP = 0;
% Primal value for filtered iterations, applicable only if
% options.precondTypeMeas[2] = true. As with above, automatically computed
% if left zero or empty.
options.tauCPFilt = 0;
% Dual value. Recommended to set at 1.
options.sigmaCP = 1;
% Next estimate update variable
options.thetaCP = 1;

% Use adaptive update of the primal and dual variables
% Currently only one method available
% Setting this to 1 uses an adaptive update for both the primal and dual
% variables.
% Can lead to unstable behavior with using multi-resolution
% Minimal to none use with filtering-based preconditioner
options.PDAdaptiveType = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRECONDITIONERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Applies to PDHG, PDHGL1, PDHGKL, PKMA, MBSREM, MRAMLA, PDDY, FISTA and
%%% FISTAL1
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
% See the article for details
options.gradV1 = 1.5;
options.gradV2 = 2;
% Note that these include subiterations (options.Niter * options.subsets)
options.gradInitIter = 1;
options.gradLastIter = 100;

% Number of filtering iterations
% Applies to both precondTypeMeas(2) and precondTypeImage(6)
options.filteringIterations = 100;


%%%%%%%%%%%%%%%%%%%%%%%%% REGULARIZATION PARAMETER %%%%%%%%%%%%%%%%%%%%%%%%
%%% The regularization parameter for ALL regularization methods (priors)
options.beta = 1;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%% NEIGHBORHOOD PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%
%%% How many neighboring pixels are considered 
% With MRP, QP, L, FMH, NLM, GGMRF and weighted mean
% E.g. if Ndx = 1, Ndy = 1, Ndz = 0, then you have 3x3 square area where
% the pixels are taken into account (I.e. (Ndx*2+1)x(Ndy*2+1)x(Ndz*2+1)
% area).
% NOTE: Currently Ndx and Ndy must be identical.
% For NLM this is often called the "search window".
options.Ndx = 1;
options.Ndy = 1;
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
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%% L-FILTER PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Weighting factors for the L-filter pixels
% Otherwise the same as in quadratic prior, but center pixel is not Inf.
% If left empty then they will be calculated by the algorithm such that the
% weights resemble a Laplace distribution.
options.a_L = [];

%%% If the weighting factors are set empty, then this option will determine
% whether the computed weights follow a 1D weighting scheme (true) or 2D 
% (false).
% See the docs for more information:
% https://omega-doc.readthedocs.io/en/latest/algorithms.html#l-filter
options.oneD_weights = false;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FMH PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pixel weights for FMH
% The matrix size needs to be [Ndx*2+1, 4] if Nz = 1 or Ndz = 0, or
% [Ndx*2+1, 13] otherwise.
% The center pixel weight should be in the middle of the weight matrix.
% If the sum of each column is > 1, then the weights will be normalized
% such that the sum = 1.
% If left empty then they will be calculated by the algorithm such that the
% weights follow the same pattern as in the original article.
options.fmh_weights = [];

%%% Weighting value for the center pixel
% Default value is 4, which was used in the original article.
% NOTE: This option is ignored if you provide your own weights.
options.fmh_center_weight = 4;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%% WEIGHTED MEAN PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%
%%% Mean type
% 1 = Arithmetic mean, 2 = Harmonic mean, 3 = Geometric mean
options.mean_type = 1;

%%% Pixel weights for weighted mean
% The number of pixels needs to be the amount of neighboring pixels,
% e.g. if the above Ndx/y/z values are all 1, then 27 weights need to be
% included. Size is (Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1). If left empty then
% they will be calculated by the algorithm such that the weights are
% dependent on the distance from the center pixel to the neighboring
% pixels.
options.weighted_weights = [];

%%% Center pixel weight for weighted mean.
% NOTE: This option is ignored if you provide your own weights.
options.weighted_center_weight = 4;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TV PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% "Smoothing" parameter
% Also used to prevent zero values in square root.
options.TVsmoothing = 1e-5;

%%% Whether to use an anatomical reference/weighting image with the TV
options.TV_use_anatomical = false;

%%% If the TV_use_anatomical value is set to true, specify filename for the
% reference image here (same rules apply as with attenuation correction
% above). Alternatively you can specifiy the variable that holds the
% reference image.
options.TV_reference_image = 'reference_image.mat';

%%% Three different TV methods are available.
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
% TVtype, e.g. for type 1 it is the edge threshold parameter. See the wiki
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
options.SATVPhi = 0.2;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADMRP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Time step variable for AD (implementation 2 only)
options.TimeStepAD = 0.0625;

%%% Conductivity/connectivity for AD (edge threshold)
options.KAD = 2;

%%% Number of iterations for AD filter
% NOTE: This refers to the AD smoothing part, not the actual reconstruction
% phase.
options.NiterAD = 10;

%%% Flux/conduction type for AD filter
% 1 = Exponential
% 2 = Quadratic
options.FluxType = 1;

%%% Diffusion type for AD (implementation 2 only)
% 1 = Gradient
% 2 = Modified curvature
options.DiffusionType = 1;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% APLS PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Scaling parameter (eta)
% See the wiki for details:
% https://omega-doc.readthedocs.io/en/latest/algorithms.html#tv
options.eta = 1e-5;

%%% "Smoothing" parameter (beta)
% Also used to prevent zero values in square root.
options.APLSsmoothing = 1e-5;

%%% Specify filename for the reference image here (same rules apply as with
% attenuation correction above). As before, this can also be a variable
% instead.
% NOTE: For APSL, the reference image is required.
options.APLS_reference_image = 'reference_image.mat';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HYPERBOLIC PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Edge weighting factor
options.hyperbolicDelta = 800;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TGV PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TGV weights
% First part
options.alpha0TGV = 1;
% Second part (symmetrized derivative)
options.alpha1TGV = 2;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NLM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filter parameter
options.sigma = 10;

%%% Patch radius
options.Nlx = 1;
options.Nly = 1;
options.Nlz = 1;

%%% Standard deviation of the Gaussian filter
options.NLM_gauss = 1;

% Search window radius is controlled by Ndx, Ndy and Ndz parameters
% Use anatomical reference image for the patches
options.NLM_use_anatomical = false;

%%% Specify filename for the reference image here (same rules apply as with
% attenuation correction above)
options.NLM_reference_image = 'reference_image.mat';

% Note that only one of the below options for NLM can be selected!
%%% Use Non-local total variation (NLTV)
% If selected, will overwrite regular NLM regularization as well as the
% below MRP version.
options.NLTV = false;

%%% Use MRP algorithm (without normalization)
% I.e. gradient = im - NLM_filtered(im)
options.NLM_MRP = false;

%%% Use non-local relative difference prior (NLRD)
options.NLRD = false;

%%% Use non-local GGMRF (NLGGMRF)
options.NLGGMRF = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RDP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Edge weighting factor
options.RDP_gamma = 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GGMRF PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GGMRF parameters
% See the original article for details
options.GGMRF_p = 1.5;
options.GGMRF_q = 1;
options.GGMRF_c = 5;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Store the intermediate forward projections. Unlike image estimates, this
% also stores subiteration results.
options.storeFP = false;
 
 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD SCATTER DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% Load user scatter data
% Load scatter data (if applicable)
if ~isfield(options,'ScatterC') && ((options.scatter_correction && ~options.only_reconstructions && ~options.corrections_during_reconstruction) ...
        || (options.scatter_correction && options.use_raw_data && ~options.corrections_during_reconstruction) ...
        || (options.scatter_correction && options.corrections_during_reconstruction && ~options.only_sinos))
  options = loadScatterData(options);
end
 
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
options.DOI = 4.584;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load the LST coincidence data

if ~options.only_reconstructions && ~options.no_data_load && options.use_machine ~= 2
    options.SinM = load_data(options);
end


%% Form the attenuation correction image from Inveon atn-file or UMAP-file

if options.attenuation_correction && isempty(options.attenuation_datafile)
    options.attenuation_datafile = attenuation_correction_factors(options);
end

%% Compute the normalization coefficients, if selected

if options.compute_normalization && ~options.only_reconstructions
    [norm_matrix, options.SinM, axial_geom_coeffs, axial_block_profile, block_profile_matrix, det_eff_coeffs, tr_geom_matrix] = normalization_coefficients(options);
    
    if options.use_raw_data
        options.coincidences{1} = options.SinM;
        options = rmfield(options, 'SinM');
    end
end

%% Reconstructions

if options.only_sinos == false
    
    tStart = tic;
    % pz contains the 3D or 4D image
    % recPar contains certain reconstruction parameters in a struct that can be
    % easily included with the reconstruction
    % options is the modified options struct (some modifications might be
    % applied before reconstruction)
    % fp are the stored forward projections, if options.storeFP = true
    [pz, recPar, ~, fp] = reconstructions_main(options);
    tElapsed = toc(tStart);
    disp(['Reconstruction process took ' num2str(tElapsed) ' seconds'])

	% save([options.name '_reconstruction_' num2str(options.subsets) 'subsets_' num2str(options.Niter) 'iterations_' ...
	%     num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '.mat'], 'pz');


    volume3Dviewer(pz, [0], [0 0 1])
    
end

%% Scale the images

% Use this function to scale the reconstructed images to Bq/mL (default
% values are number of counts / voxel). Finite total time is required.
% pz = scaleImages(pz);

% If infinite time was used, input the total time of the examination as the
% second input
% pz = scaleImages(pz, total_time);