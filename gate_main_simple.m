%% MATLAB codes for GATE PET reconstruction using ASCII or LMF output

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
%%% Machine name
options.machine_name = 'Cylindrical_PET_example';
 
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
options.obtain_trues = true;
%%% Obtain scattered coincidences
% If this is set to true, then scattered coincidences are saved separately
% These events are not used for scatter correction though, but a separate
% scatter sinogram/raw data matrix will be created.
options.store_scatter = true;
%%% Reconstruct the true coincidences
% If this is set to true, then the true coincidences will be used for
% reconstruction
% NOTE: If both this and reconstruct_scatter are set, then the trues are
% reconstructed, but not the scatter.
options.reconstruct_trues = false;
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
options.store_randoms = true;
%%% Obtain source coordinates (used in forming the "true" image).
% If this is set to true, then the "true" decay image is also saved during
% data load, i.e. the locations where the decay has occurred.
% If any of the above settings are set to true, then the true images are
% also obtained for them. E.g. if store_scatter = true, then an image
% showing the locations and number of counts of where the scattered events
% originated will be saved in a mat-file. Scatter and trues contain
% coincidence events while randoms contain singles.
% NOTE: If you use LMF data, the source images are not considered reliable.
options.source = true;
 
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
% The variables below are used only if the above variable is set to true
%%% How many bytes are in the LMF header part?
options.header_bytes = (16);
%%% How many bytes in each event packet?
% From left to right: time + detector ID + source XYZ position + event ID +
% N Compton in phantom
options.data_bytes = (8 + 2 + 6 + 4 + 1);
%%% How many bits dedicated to R-sectors?
options.R_bits = 7;
%%% How many bits dedicated to modules?
options.M_bits = 2;
%%% How many bits dedicated to submodules?
options.S_bits = 1;
%%% How many bits dedicated to crystals?
options.C_bits = 6;
%%% How many bits dedicated to layers?
options.L_bits = 0;
%%% What is the coincidence window (in seconds)?
options.coincidence_window = 10e-9;
%%% What is the clock time step? (see the .cch files)
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
% NOTE 2: If you are using MATLAB R2018b or earlier, ROOT will eventually
% cause MATLAB to crash. This can be circumvent by running MATLAB with
% matlab -nojvm
options.use_root = false;
 
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
options.Nz = 63;
%%% Flip the image (in vertical direction)?
options.flip_image = false;
%%% How much is the image rotated?
% You need to run the precompute phase again if you modify this
% NOTE: The rotation is done in the detector space (before reconstruction)
options.offangle = 0;
 
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
%%% Number of angular positions (views) in sinogram
% You should primarily use the same number as the device uses
% However, if that information is not available you can use ndist_max
% function to determine potential values (see help ndist_max for usage)
options.Ndist = 200;
%%% Number of angles (tangential positions) in sinogram 
% This is the final amount after possible mashing, maximum allowed is the
% number of detectors per ring/2
options.Nang = 168;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MISC PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% Name of current datafile/examination
options.name = 'cylpet_example';
%%% Folder for the data (.dat ASCII, .ccs LMF, .root ROOT) files
% If no files are located in the path provided below, then the current
% folder is also checked. If no files are detected there either, an error
% is thrown.
if ispc % Windows
    options.fpath = 'C:\path\to\GATE\output\';
else % Unix/MAC
    options.fpath = '/path/to/GATE/output/';
end
%%% Form only sinograms and raw data matrix (no reconstructions)
options.only_sinos = false;
%%% Compute only the reconstructions
% If this file is run with this set to true, then the data load and
% sinogram formation steps are always skipped. Use this only after the data
% has been imported and saved
options.only_reconstructions = false;
%%% Use raw list-mode data
% This means that the data is used as is without any sinogramming and thus
% without any "compression"
options.use_raw_data = false;
%%% Show status messages
options.verbose = true;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% Number of iterations
options.Niter = 4;
%%% Number of subsets
options.subsets = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up necessary parameters that are not adjusted in the simple version
options = set_up_parameters_simple(options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ERROR CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic error checking is done here
options = OMEGA_error_check(options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load the ASCII/LMF/ROOT coindicence data

if options.only_reconstructions == false
    options.coincidences = load_data(options);
end

%% Form the sinograms

if options.only_reconstructions == false && options.use_raw_data == false
    options.SinM = form_sinograms(options);
end

%% Reconstructions

if ~options.only_sinos
    tStart = tic;
    f_osem = reconstructions_main_simple(options);
    tElapsed = toc(tStart);
    disp(['Reconstruction process took ' num2str(tElapsed) ' seconds'])
end
