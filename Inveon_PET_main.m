%% MATLAB codes for PET reconstruction using Inveon PET LST or GATE output

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% MACHINE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% R-sectors/modules/blocks/buckets in transaxial direction
options.blocks_per_ring = (16);

%%% R-sectors/modules/blocks/buckets in axial direction (i.e. number of physical
%%% machine/crystal rings) 
% Multiplying this with the below cryst_per_block should equal the total
% number of crystal rings. 
options.linear_multip = (4);

%%% Number of detectors on the side of R-sector/block/module (transaxial
%%% direction)
% (e.g. 13 if 13x13, 20 if 20x10)
options.cryst_per_block = (20);

%%% Crystal pitch/size in x- and y-directions (transaxial) (mm)
options.cr_p = 1.63;

%%% Crystal pitch/size in z-direction (axial) (mm)
options.cr_pz = 1.592;

%%% Ring diameter (distance between perpendicular detectors) (mm)
options.diameter = 161.08;

%%% Transaxial FOV size (mm), this is the length of the x (horizontal) side
% of the FOV
options.FOVa_x = 100;

%%% Transaxial FOV size (mm), this is the length of the y (vertical) side
% of the FOV
options.FOVa_y = options.FOVa_x;

%%% Axial FOV (mm)
options.axial_fov = 127;

%%% Number of pseudo rings between physical rings (use 0 or [] if none)
% NOTE: Inveon has no pseudo detectors/rings
options.pseudot = [];

%%% Number of detectors per ring (without pseudo detectors)
options.det_per_ring = options.blocks_per_ring*options.cryst_per_block;

%%% Number of detectors per ring (with pseudo detectors)
% NOTE: Inveon has no pseudo detectors/rings
options.det_w_pseudo = options.blocks_per_ring*(options.cryst_per_block);

%%% Number of crystal rings
options.rings = options.linear_multip * options.cryst_per_block;

%%% Number of detectors
options.detectors = options.det_per_ring*options.rings;

%%% Machine name
% Used for naming purposes (measurement data)
options.machine_name = 'Inveon';

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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% GATE SPECIFIC SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Obtain trues (true coincidences)
% If this is set to true then, in addition to the normal coincidences
% (prompts), trues are also obtained and saved.
options.obtain_trues = false;

%%% Reconstruct the true coincidences
% If this is set to true, then the true coincidences will be used for
% reconstruction
% NOTE: If both this and reconstruct_scatter are set, then the trues are
% reconstructed, but not the scatter.
options.reconstruct_trues = false;

%%% Obtain scattered coincidences
% If this is set to true, then scattered coincidences are saved separately.
% These events are not used for scatter correction though, but a separate
% scatter sinogram/raw data matrix will be created. The scatter elements
% included can be selected below.
options.store_scatter = false;

%%% What scatter components are included in the scatter part
% (1 means that component is included, 0 means it is not included in the
% scatter data) 
% First: Compton scattering in the phantom, second: Compton scattering in
% the detector, third: Rayleigh scattering in the phantom, fourth: Rayleigh
% scattering in the detector.
% If store_scatter is set to true, at least one value has to be 1. E.g. [1
% 0 1 0] will save Compton and Rayleigh scattering in the phantom. 
% NOTE: LMF will always include only Compton scattering in the phantom,
% regardless of the choice below (as long as scatter is selected).
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% Is LMF data used (Only one data type can be used at a time)
options.use_LMF = false;

% The variables below are used only if the above variable is set to true
%%% How many bytes are in the LMF header part?
options.header_bytes = (16);

%%% How many bytes in each event packet?
% From left to right: time + detector ID + source XYZ position + event ID +
% N Compton in phantom
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% Is ROOT data used (Only one data type can be used at a time)
% NOTE: On Windows ROOT works only with 32-bit MATLAB/Octave, but has not
% been tested with it.
% NOTE 2: If you are using MATLAB R2018b or earlier, ROOT will eventually
% cause MATLAB to crash. This can be circumvent by running MATLAB with
% matlab -nojvm. 2019a and up are unaffected, GNU Octave is unaffected.
options.use_root = false;
 
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
 
%%% Reconstructed image pixel size (X-direction)
options.Nx = 128;

%%% Y-direction
options.Ny = 128;

%%% Z-direction (number of slices) (axial)
options.Nz = options.rings*2 - 1;

%%% Flip the image (in vertical direction)?
options.flip_image = false;

%%% How much is the image rotated?
% You need to run the precompute phase again if you modify this
% NOTE: The rotation is done in the detector space (before reconstruction).
% This current setting is for machine list-mode data or sinogram data.
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
options.Ndist = 128;

%%% Number of angles (tangential positions) in sinogram 
% This is the final amount after possible mashing, maximum allowed is the
% number of detectors per ring/2.
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
options.NSinos = options.TotSinos;

%%% If Ndist value is even, take one extra out of the negative side (+1) or
% from the positive side (-1). E.g. if Ndist = 200, then with +1 the
% interval is [-99,100] and with -1 [-100,99]. This varies from device to
% device. If you see a slight shift in the sinograms when comparing with
% the machine sinograms then use the other option here. For Inveon, this
% should be -1.
options.ndist_side = -1;

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
 
%%% Maximum ring difference
options.ring_difference_raw = options.rings;

%%% Increase the sampling rate of the raw data
% Increasing this interpolates additional rows and columns to the raw data
% Can be used to prevent aliasing artifacts
% NOTE: Has to be either 1 or divisible by two
options.sampling_raw = 1;

%%% Interpolation method used for sampling rate increase
% All the methods are available that are supported by interp2
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
options.randoms_correction = true;

%%% Variance reduction
% If set to true, variance reduction will be performed to delayed
% coincidence (randoms corrections) data if randoms correction is selected
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
% create_atten_matrix_CT functions. If the below attenuation_datafile is
% empty, UMAPs created with Inveon AW can be used (if the below
% CT_attenuation is set to true) or, alternatively, use the produces .atn
% files instead (CT_attenuation is false). 
% E.g. if this is set to true, CT_attenuation = true and
% attenuation_datafile = '', then the user will be automatically prompted
% for the UMAP-files and they will be automatically saved as a mat-file
% with the filename saved in the attenuation_datafile field. Alternatively,
% if this is set to true, CT_attenuation = false and  attenuation_datafile
% = '', then the user will be automatically prompted for the .atn-files
% from which the attenuation images will be automatically created and
% saved (filename is saved in the attenuation_datafile field).
options.attenuation_correction = true;

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
% coefficients for the specified machine will be automatically loaded. Use
% this only if you want to use normalization coefficients computed outside
% of OMEGA.
% NOTE: Supports .nrm files created by the Inveon AW. 
options.use_user_normalization = true;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Arc correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Apply arc correction
% NOTE: Arc correction is an experimental feature. It is currently
% relatively slow and supports only sinogram data. Generally it is not
% recommended to use arc correction (Inveon data is an exception).
% Uses parallel computing toolbox if it is available (parfor).
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
% applied.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MISC PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% Name of current datafile/examination
% This is used to name the saved measurement data and also load it in
% future sessions.
options.name = 'open_PET_data';

%%% Precompute data 
% This should be done when using data from a certain machine the first time
% as it can speed up reconstruction. Especially recommended for raw
% list-mode data. Not mandatory and the precomputed data is only used if 
% the below precompute_lor is set to true. If using solely implementation 
% 1, this is HIGHLY recommended. 
options.precompute = true;

%%% Folder for the data (.dat ASCII, .ccs LMF, .root ROOT) files
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

%%% Form only sinograms and raw data matrix (no reconstructions)
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

%%% Precompute the observation/system matrix for the reconstruction (this
%%% might require a lot of memory), if false then the observation matrix is
%%% calculated on the fly. 
% NOTE: Supports only MLEM reconstruction and is an experimental feature.
% This applies ONLY when using implementation 1.
options.precompute_obs_matrix = false;

%%% Compute only the reconstructions
% If this file is run with this set to true, then the data load and
% sinogram formation steps are always skipped. Precomputation step is
% only performed if precompute_lor = true and precompute_all = true
% (below). Normalization coefficients are not computed even if selected.
options.only_reconstructions = false;

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
% (Supports only one algorithm at a time)
% See the wiki for more information: 
% https://github.com/villekf/OMEGA/wiki/Useful-information#selecting-the-correct-implementation
options.implementation = 4;

% Applies to implementations 2 and 3 ONLY
%%% Device used (this is applicable to implementation 2), or platform used
% (implementation 3)
% In implementation 2 this determines the device used for both system
% matrix formation and image reconstruction.
% NOTE: Use ArrayFire_OpenCL_device_info() to determine the device numbers.
% In implementation 3, this determines the platform from where the
% device(s) are taken.
% NOTE: Use OpenCL_device_info() to determine the platform numbers and
% their respective devices.
% NOTE: if you switch devices then you need to run the below line
% (uncommented) as well:
% clear mex
% NOTE: Using implementation 3 after using 2 might fail until MATLAB is
% restarted.
options.use_device = 0;

% Applies to implementations 2 and 3 ONLY
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

% Implementation 3 ONLY
%%% How many times more measurements/LORs are in the GPU part (applicable if
% heterogenous computing is used) 
% Alternatively, set this to 0 to use only a single device on the specific
% platform (the one with the highest memory count will be used)
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
options.use_psf = true;

% FWHM of the Gaussian used in PSF blurring in all three dimensions
options.FWHM = [options.cr_p options.cr_p options.cr_pz];

% Use deblurring phase
% If enabled, a deblurring phase is performed once the reconstruction has
% completed. This step is performed for all iterations (deblurred estimates
% are NOT used in the reconstruction phase). This is used ONLY when PSF
% blurring is used.
options.deblurring = false;
% Number of deblurring iterations
% How many iterations of the deblurring step is performed
options.deblur_iterations = 10;

% Orthogonal ray tracer (projector_type = 2) only
%%% The 2.5D (XY) width of the "strip/tube" where the orthogonal distances are
% included. If the tube_width_z is non-zero, then this value is ignored.
options.tube_width_xy = options.cr_p;

% Orthogonal ray tracer (projector_type = 2) only
%%% The 3D (Z) width of the "tube" where the orthogonal distances are
% included. If set to 0, then the 2.5D orthogonal ray tracer is used. If this
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
% Implementations 2 and 3 ONLY
%%% Apply acceleration
% If true, then intermediate results are saved in memory. If you run out
% memory, you should set this to false. Does not apply to improved Siddon.
% Applies only if using either implementation 2 or 3. Can speed up
% computations by around 30% if set to true.
options.apply_acceleration = true;
 

%%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Number of iterations (all reconstruction methods)
options.Niter = 10;
% Save ALL iterations
% Set this to false if you do not want to save all the intermediate
% iterations, but only the very last one.
options.save_iter = true;

%%% Number of subsets (all excluding MLEM and subset_type = 5)
options.subsets = 16;

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
% NOTE: Applies only when using subset_type = 3; accelerates the subset
% formation and uses less memory. Not included in OMEGA, needs to be
% manually downloaded and installed. Enabling this will automatically cause
% the function to be called if it is found on MATLAB path.
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

%%% Skip the normalization phase in MRP, FMH, L-filter, ADMRP and
%%% weighted mean
% E.g. if set to true the MRP prior is (x - median(x))
% E.g. if set to false the MRP prior is (x - median(x)) / median(x)
% The published MRP uses the one that is obtained when this is set to
% false, however, you might get better results with true. I.e. this should
% be set to false if you wish to use the original prior implementation.
options.med_no_norm = false;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION ALGORITHMS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction algorithms to use (you can choose several)
% NOTE: MLEM requires precomputed observation matrix or a matrix-free
% method
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ML-METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are non-regularized versions
%%% Maximum-Likelihood Expectation Maximization (MLEM)
% Supported by all implementations
% For implementation 1 requires the precomputed observation/system matrix
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
% Supported by implementations 1, 2 and 4
options.drama = false;

%%% Complete data OSEM (COSEM)
% Supported by implementations 1, 2 and 4
options.cosem = false;

%%% Enhanced COSEM (ECOSEM)
% Supported by implementations 1, 2 and 4
options.ecosem = false;

%%% Accelerated COSEM (ACOSEM)
% Supported by implementations 1, 2 and 4
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

%%% RBI-OSL
% Supported by implementations 1, 2 and 4
options.RBI_OSL = false;

%%% (A)COSEM-OSL
% 0/false = No COSEM-OSL, 1/true = ACOSEM-OSL, 2 = COSEM-OSL
% Supported by implementations 1, 2 and 4
options.COSEM_OSL = false;
 
 
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

%%% Total Generalized Variation (TGV) prior
options.TGV = false;

%%% Non-local Means (NLM) prior
options.NLM = false;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% ACOSEM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Acceleration parameter for ACOSEM (1 equals COSEM)
options.h = 2;
 

%%%%%%%%%%%%%%%%%%%%%%%% MRAMLA & MBSREM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%
%%% Relaxation parameter for MRAMLA and MBSREM
% Use scalar if you want it to decrease as
% lambda0_mbsrem/current_iteration_number. Use vector (length = Niter) if
% you want your own relaxation parameters.
options.lambda0_mbsrem = 0.2;

%%% Upper bound for MRAMLA/MBSREM (use 0 for default (computed) value)
options.U = 0;
 

%%%%%%%%%%%%%%%%%%%%%%%%% RAMLA & BSREM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%
%%% Relaxation parameter for RAMLA and BSREM
% Use scalar if you want it to decrease as
% lambda0/current_iteration_number. Use vector (length = Niter) if you want
% your own relaxation parameters.
options.lambda0 = 0.2;
 

%%%%%%%%%%%%%%%%%%%%%%% ROSEM & ROSEM-MAP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%
%%% Relaxation parameter for ROSEM and ROSEM-MAP
% Use scalar if you want it to decrease as
% lambda0_rosem/current_iteration_number. Use vector (length = Niter) if
% you want your own relaxation parameters.
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
% area).
% NOTE: Currently Ndx and Ndy must be identical.
% For NLM this is often called the "search window".
options.Ndx = 5;
options.Ndy = 5;
options.Ndz = 5;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MRP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for MRP with OSL-OSEM
options.beta_mrp_osem = 0.1;
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
options.beta_quad_osem = 0.01;
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
% The number of pixels need to be the amount of neighboring pixels,
% e.g. if the above Nd values are all 1, then 27 weights need to be
% included where the center pixel (if Nd values are 1, element 14) should
% be Inf. Size is (Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1). If left empty then
% they will be calculated by the algorithm and are based on the distance of
% the voxels from the center.
options.weights = [];
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for Huber prior with OSL-OSEM
options.beta_huber_osem = 0.01;
%%% Regularization parameter for Huber prior with OSL-MLEM
options.beta_huber_mlem = 0.1;
%%% Regularization parameter for Huber prior with MBSREM
options.beta_huber_mbsrem = 0.05;
%%% Regularization parameter for Huber prior with BSREM
options.beta_huber_bsrem = 0.03;
%%% Regularization parameter for Huber prior with ROSEM
options.beta_huber_rosem = 0.1;
%%% Regularization parameter for Huber prior with RBI
options.beta_huber_rbi = 0.05;
%%% Regularization parameter for Huber prior (OSL-(A)COSEM)
options.beta_huber_cosem = 0.01;

%%% Delta parameter for Huber prior
% Upper and lower bounds for the prior
options.huber_delta = 5;

%%% Pixel weights for Huber prior
% Same rules apply as with quadratic prior weights.
% If left empty then they will be calculated by the algorithm and are based
% on the distance of the voxels from the center.
options.weights_huber = [];
 
 
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
% Otherwise the same as in quadratic prior, but center pixel is not Inf.
% If left empty then they will be calculated by the algorithm such that the
% weights resemble a Laplace distribution.
options.a_L = [];

%%% If the weighting factors are set empty, then this option will determine
% whether the computed weights follow a 1D weighting scheme (true) or 2D 
% (false).
% See the wiki for more information:
% https://github.com/villekf/OMEGA/wiki/Function-help#reconstruction-algorithms
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

%%% Regularization parameter for weighted mean with OSL-OSEM
options.beta_weighted_osem = 0.2;
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
%%% Regularization parameter for TV OSL-OSEM
options.beta_TV_osem = 0.01;
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
% Also used to prevent zero values in square root.
options.TVsmoothing = 1e-1;

%%% Whether to use an anatomical reference/weighting image with the TV
options.TV_use_anatomical = true;

%%% If the TV_use_anatomical value is set to true, specify filename for the
% reference image here (same rules apply as with attenuation correction
% above).
options.TV_reference_image = 'reference_image.mat';

%%% Three different TV methods are available.
% Value can be 1, 2 or 3.
% Types 1 and 2 are the same if anatomical prior is not included
% Type 3 uses the same weights as quadratic prior
% See the wiki for more information:
% https://github.com/villekf/OMEGA/wiki/Function-help#reconstruction-algorithms
options.TVtype = 1;

%%% Weighting parameters for the TV prior. 
% Applicable only if use_anatomical = true. T-value is specific to the used
% TVtype, e.g. for type 1 it is the edge threshold parameter. See the wiki
% for more details:
% https://github.com/villekf/OMEGA/wiki/Function-help#reconstruction-algorithms
options.T = 0.5;

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
% phase.
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
% See the wiki for details:
% https://github.com/villekf/OMEGA/wiki/Function-help#reconstruction-algorithms
options.eta = 1e-5;

%%% "Smoothing" parameter (beta)
% Also used to prevent zero values in square root.
options.APLSsmoothing = 1e-5;

%%% Specify filename for the reference image here (same rules apply as with
% attenuation correction above)
% NOTE: For APSL, the prior is required.
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
%%% Regularization parameter for NLM with OSL-OSEM
options.beta_NLM_osem = 0.01;
%%% Regularization parameter for NLM with OSL-MLEM
options.beta_NLM_mlem = 0.01;
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
options.sigma = 10;

%%% Patch radius
options.Nlx = 2;
options.Nly = 2;
options.Nlz = 2;

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
% below MRP version.
options.NLTV = false;

%%% Use MRP algorithm (without normalization)
% I.e. gradient = im - NLM_filtered(im)
options.NLM_MRP = false;
 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD SCATTER DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% Load user scatter data
% Load scatter data (if applicable)
if (options.scatter_correction && ~options.only_reconstructions && ~options.corrections_during_reconstruction) ...
        || (options.scatter_correction && options.use_raw_data && ~options.corrections_during_reconstruction) ...
        || (options.scatter_correction && options.corrections_during_reconstruction)
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
% size(x) = [(number of elements in a SINLGE sinogram) 2]
% size(y) = [(number of elements in a SINLGE sinogram) 2]
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
% options.DOI = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Precompute data if necessary

if (options.precompute && options.only_reconstructions == false && ~options.only_sinos) || (options.precompute_all && options.precompute)
    precompute_data(options);
end

%% Load the ASCII/LMF/ROOT/LST coincidence data

if ~options.only_reconstructions && ~options.no_data_load && options.use_machine ~= 2
    options.coincidences = load_data(options);
    % Save the delayed coincidences, trues, scatter and randoms to
    % workspace (GATE only):
    % [options.coincidences, delayed_coincidences, true_coincidences, scattered_coincidences, random_coincidences] = load_data(options);
    % Save the coordinates of the detected coincidences to work space (GATE only):
    % [options.coincidences, ~, ~, ~, ~, x, y, z] = load_data(options);
end

%% Form the sinograms and perform corrections if selected

%%% Forms the sinogram from the raw data.
%%% Should be used only to create sinograms of different size as
%%% previously, for corrections use the below version.
% if ~options.only_reconstructions && ~options.use_raw_data && ~isfield(options,'coincidences') && ~options.no_data_load
%    options.SinM = form_sinograms(options);
% end

%%% Use this function to perform only selected corrections to the raw
%%% uncorrected sinogram. This overwrites any previous corrections applied
%%% to the sinogram. Using this prevents the need to completely compute the
%%% sinogram again. The uncorrected sinogram remains untouched.
% if ~options.corrections_during_reconstruction && ~options.use_raw_data && (options.randoms_correction || options.scatter_correction || options.normalization_correction) ...
%         && ~options.only_sinos
%     options.SinM = form_sinograms(options, true);
% end


%% Form the attenuation correction image from Inveon atn-file or UMAP-file

if options.attenuation_correction && ~options.only_reconstructions && isempty(options.attenuation_datafile)
    options.attenuation_datafile = attenuation_correction_factors(options);
end

%% Compute the normalization coefficients

if options.compute_normalization && ~options.only_reconstructions
    [norm_matrix, options.SinM, axial_geom_coeffs, axial_block_profile, block_profile_matrix, det_eff_coeffs, tr_geom_matrix] = normalization_coefficients(options);
    
    if options.use_raw_data
        options.coincidences{1} = options.SinM;
        options = rmfield(options, 'SinM');
    end
end

%% Reconstructions

if ~options.use_raw_data && isfield(options,'coincidences')
    options = rmfield(options, 'coincidences');
end

if options.only_sinos == false
    
    tStart = tic;
    pz = reconstructions_main(options);
    tElapsed = toc(tStart);
    disp(['Reconstruction process took ' num2str(tElapsed) ' seconds'])
    
end

% save([options.name '_reconstruction_' num2str(options.subsets) 'subsets_' num2str(options.Niter) 'iterations_' ...
%     num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '.mat'], 'pz');