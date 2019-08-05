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
% machine name
options.machine_name = 'Cylindrical_PET_example';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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
% What is the clock time step? (see the .cch files)
% If e.g. 1 ps then use 1e-12, if 1 ns use 1e-9, etc.
options.clock_time_step = 1e-12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Root data format settings %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Is root data used? (Only one data type can be used at a time)
options.use_root = false;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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
% If no files are located in the path provided below, then the current
% folder is also checked
if ispc % Windows
%     options.fpath = 'C:\path\to\GATE\output\';
    options.fpath = 'C:\Users\villewe\OneDrive - University of Eastern Finland\MAT-tiedostot\PET\GATE\example\Uus';
%     options.fpath = 'I:\Väikkäri\MAT-tiedostot\PET\example\Uus';
else % Unix/MAC
    options.fpath = '/path/to/GATE/output/';
end
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
% Show status messages
options.verbose = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Error checking %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic error checking is done here
options = OMEGA_error_check(options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Precompute the necessary data

if options.precompute && options.only_reconstructions == false
    precompute_data(options);
end

%% Load the ASCII/LMF coindicence data

if options.only_reconstructions == false
    options.coincidences = load_data(options);
end

%% Form the sinograms

if options.only_reconstructions == false && options.use_raw_data == false
    options.SinM = form_sinograms(options);
end

%% Reconstructions
    
tStart = tic;
f_osem = reconstructions_main_simple(options);
tElapsed = toc(tStart);
disp(['Reconstruction process took ' num2str(tElapsed) ' seconds'])


% save([options.name '_reconstruction_' num2str(options.subsets) 'subsets_' num2str(options.Niter) 'iterations.mat'], 'pz');