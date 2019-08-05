%% MATLAB codes for PET reconstruction using sinogram input from any machine

clear

%%%%%%%%%%%%%%%%%%%% Specify the below machine properties %%%%%%%%%%%%%%%%%%%%
% Blocks in transaxial direction
options.blocks_per_ring = (42);
% Blocks in axial direction (i.e.  number of physical machine rings)
options.linear_multip = (4);
% number of detectors on the side of block (e.g. 13 if 13x13)
options.cryst_per_block = (8);
% crystal pitch in x- and y-directions (mm)
options.cr_p = 2.4;
% crystal pitch in z-direction (mm)
options.cr_pz = 2.4;
% ring diameter (distance between perpendicular detectors) in mm
options.diameter = 130*2;
% Transaxial FOV size (mm), assuming square FOV (this is the length of one
% side of the square FOV)
options.FOVa = 151;
% Axial FOV (mm)
options.axial_fov = round(76.8 - options.cr_pz/2);
% Number of pseudo rings between physical rings (use 0 or [] if none)
options.pseudot = [];
% Number of detectors per ring (without pseudo detectors)
options.det_per_ring = options.blocks_per_ring*options.cryst_per_block;
% Number of detectors per ring (with pseudo detectors)
options.det_w_pseudo = options.blocks_per_ring*(options.cryst_per_block);
% number of crystal rings
options.rings = options.linear_multip * options.cryst_per_block;
% number of detectors
options.detectors = options.det_per_ring*options.rings;
% machine name
options.machine_name = 'Cylindrical_PET_example';





%%%%%%%%%%%%%%%%%%%% Image properties %%%%%%%%%%%%%%%%%%%%
% reconstructed image pixel size (X-direction)
options.Nx = 128;
% Y-direction
options.Ny = 128;
% Z-direction (number of slices)
options.Nz = options.rings*2-1;





%%%%%%%%%%%%%%%%%%%% Sinogram properties %%%%%%%%%%%%%%%%%%%%
% span factor/axial compression
options.span = 3;
% maximum ring difference
options.ring_difference = 29;
% number of angles in sinogram (this is the final amount after possible
% mashing), maximum allowed is the number of detectors per ring
options.Nang = options.det_per_ring/2;
% number of distances in sinogram
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
% interval is [-99,100].
% This option has almost no effect, default value can be used unless you
% know the value the device uses
options.ndist_side = 1;





%%%%%%%%%%%%%%%%%%%% Dynamic imaging properties %%%%%%%%%%%%%%%%%%%%
options.tot_time = 1800; % Total time of the measurement (s)
options.partitions = 1;% how many time points (if a static measurement, use 1)





%%%%%%%%%%%%%%%%%%%% Misc properties %%%%%%%%%%%%%%%%%%%%
% name of current datafile/examination
options.name = 'cylpet_example';

% precompute data (this should be done when using data from a certain
% machine the first time)
% This includes e.g. detector coordinates, sinogram coordinates, etc.
options.precompute = true;

% Precompute the observation matrix for the reconstruction (this might require a
% lot of memory), if false then the observation matrix is calculated on the
% fly (slower)
% NOTE: Supports only MLEM reconstruction
options.precompute_obs_matrix = false;

% Compute only the reconstructions (this option overwrites the precompute
% option)
options.only_reconstructions = false;

% Compute all the algorithms individually
options.single_reconstructions = false;

% Use raw list mode data
% This means that the data is used as is without any sinogramming and thus
% without any "compression"
options.use_raw_data = false;

% Compute only the system matrix (this option overwrites the only
% reconstructions option)
options.only_system_matrix = false;

% Use precomputed geometrical matrix information
% During the precompute-phase the number of pixels each LOR traverse is
% counted (this phase requires the above precompute-option to true). These
% are saved and later on used during the reconstruction process (compulsory
% for OpenCL reconstruction). Once the precompute-phase has been done once
% for the specific sinogram and image resolution, it is not necessary to do
% it again
options.precompute_lor = true;

% Set this option to true if you want to precompute all possible
% combinations in one go (i.e. raw data, precomputed LORs, sinogram format)
% Requires for precompute option to be true to take effect
options.precompute_all = false;

% Show status messages
options.verbose = true;





%%%%%%%%%%%%%%%%%%%% Reconstruction properties %%%%%%%%%%%%%%%%%%%%
% Reconstruction method used
% 1 = Reconstructions in MATLAB
% 2 = Everything done in OpenCL (CPU or GPU)
% 3 = Reconstructions done sequentially (matrix-free method)
% 4 = Matrix-free reconstruction with OpenCL (parallel)
options.reconstruction_method = 1;
% Device used (this is applicable to methods 2 and 4)
% In methods 2 and 4 this determines the device used for both system matrix
% formation and image reconstruction
% NOTE: On OpenCL the first devices are your GPUs (if you have only 1 GPU,
% then device number 0 is your GPU and 1 your CPU)
options.use_device = 0;
% Reconstruction algorithms to use (you can choose several):
% NOTE: MLEM requires precomputed observation matrix or a matrix-free
% method
% Maximum-Likelihood Expectation Maximization (MLEM)
% Supported by methods 1, 2 and 4
options.mlem = false;
% Ordered Subsets Expectation Maximization (OSEM)
% Supported by all methods
options.osem = true;
% Modified Row-Action Maximum Likelihood Algorithm (MRAMLA, modified BSREM)
% Supported by method 1 only
options.mramla = true;
% Row-Action Maximum Likelihood Algorithm (RAMLA)
% Supported by method 1 only
options.ramla = true;
% Enhanced Convergent OSEM (ECOSEM)
% Supported by method 1 only
options.ecosem = true;
% Complete data OSEM (COSEM)
% Supported by method 1 only
options.cosem = true;
% Accelerated COSEM (ACOSEM)
% Supported by method 1 only
options.acosem = true;
% Median root prior (MRP) with One Step Late (OSL) algorithm
% Supported by method 1 only
options.mrp_osl = true;
% Median root prior with Block Sequential Regularized Expectation
% Supported by method 1 only
% Maximization (BSREM)
options.mrp_bsrem = true;
% Quadratic prior with OSL
% Supported by method 1 only
options.quad_osl = true;
% Quadratic prior with BSREM
% Supported by method 1 only
options.quad_bsrem = true;
% L-filter with OSL
% Supported by method 1 only
options.L_osl = true;
% L-filter with BSREM
% Supported by method 1 only
options.L_bsrem = true;
% Finite impulse response (FIR) Median Hybrid (FMH) with OSL
% Supported by method 1 only
options.FMH_osl = true;
% Weighted mean with OSL
% Supported by method 1 only
options.weighted_mean_osl = true;
% number of iterations
options.Niter = 2;
% number of subsets (all excluding MLEM, use 1 if MLEM is selected)
options.subsets = 8;
% Initial value for the reconstruction
options.x0 = ones(options.Nx, options.Ny, options.Nz);
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
% NOTE: Applies only when using raw data
% Download from: https://se.mathworks.com/matlabcentral/fileexchange/27076-shuffle
options.use_Shuffle = true;
% use fast sparse? (recommended)
% Download from: https://github.com/stefanengblom/stenglib
% NOTE: This applies only to method 1 when precompute_lor is false
options.use_fsparse = false;
% accleration parameter for ACOSEM
options.h = 2;
% relaxation parameter for RAMLA
options.b0 = 1;
% Upper bound for MRAMLA
options.U = 0;
% regularization parameter for MRP with OSL
options.beta_mrp_osl = 0.05;
% regularization parameter for MRP with BSREM
options.beta_mrp_bsrem = 0.05;
% distance of where the median is taken on MRP
options.medx = 3;
options.medy = 3;
options.medz = 3;
% regularization parameter for quadratic prior with OSL
options.beta_quad_osl = 0.005;
% regularization parameter for quadratic prior with BSREM
options.beta_quad_bsrem = 0.05;
% how many neighboring pixels are taken into account (with quadratic, MRP, L, FMH and weighted)
% NOTE: At the moment Ndx and Ndy has to be identical
options.Ndx = 1;
options.Ndy = 1;
options.Ndz = 1;
% pixel weights for quadratic prior
% the number of pixels need to be the amount of neighboring pixels
% e.g. if the above Nd values are all 1, then 27 weights need to be included
% where the center pixel (if Nd values are 1, element 14) should be Inf
% Size is (Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)
% If left empty then they will be calculated by the algorithm
options.weights = [];
% regularization parameter for L-filter with OSL
options.beta_L_osl = 0.05;
% regularization parameter for L-filter with BSREM
options.beta_L_bsrem = 0.05;
% weighting factors for the L-filter pixels
% Otherwise the same as in quadratic prior, but center pixel is not Inf
options.a_L = [];
% pixel weights for FMH
% The size needs to be [3, floor((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)/2)]
options.fmh_weights = [];
% additional weighting factor for the center row of fmh_weights
% if 1 then the maximum value of all weights are used as the center row weight
% In other words it multiplies the maximum value of all other weights
options.fmh_center_weight = 1;
% regularization parameter for FMH with OSL
options.beta_fmh_osl = 0.05;
% regularization parameter for weighted mean with OSL
options.beta_weighted_osl = 0.05;
% pixel weights for weighted mean
% the number of pixels need to be the amount of neighboring pixels
% e.g. if the above Nd values are all 1, then 27 weights need to be included
% where the center pixel (if Nd values are 1, element 14) should be Inf
% Size is (Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)
% If left empty then they will be calculated by the algorithm
options.weighted_weights = [];
% center pixel weight for weighted mean
% if 1 then the maximum value of all weights are used as the center weight
% In other words it multiplies the maximum value of all other weights
options.weighted_center_weight = 1;


if options.only_system_matrix == false
    % Load the sinogram/raw data
    [options.file, options.fpath] = uigetfile('*.mat','Select PET Sinogram or list-mode data');
    
    FileName = fullfile(options.fpath, options.file);
    
    storedStructure = load(FileName);
    variables = fields(storedStructure);
    
    if options.use_raw_data == false
        options.SinM = storedStructure.(variables{1});
    else
        options.coincidences = storedStructure.(variables{1});
    end
    
end


%% Precompute the necessary data

if options.precompute && options.only_reconstructions == false
    precompute_data_nongate(options);
end

%% Reconstructions

if options.only_system_matrix == false
    
    if options.use_raw_data
        load([options.machine_name '_detector_locations_' num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '_raw.mat'],'discard')
        if iscell(options.coincidences)
            if numel(options.coincidences{1}) == numel(discard)
                options.SinM = options.coincidences(discard);
            end
        else
            if numel(options.coincidences{1}) == numel(discard)
                options.coincidences = options.coincidences(discard);
            end
        end
        clear discard
    end
    
    tStart = tic;
    pz = reconstructions_main(options);
    tElapsed = toc(tStart);
    disp(['Reconstruction process took ' num2str(tElapsed) ' seconds'])
    
end

save([options.name '_reconstruction_' num2str(options.subsets) 'subsets_' num2str(options.Niter) 'iterations.mat'], 'pz');

%% System matrix formation

if options.only_system_matrix
    %%
    LL = [];
    index = [];
    pituus = [];
    lor = [];
    
    if options.use_raw_data == false && options.subsets > 1
        if options.precompute_lor || options.reconstruction_method == 3
            load([options.machine_name '_lor_pixel_count_' num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '_sino_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'lor','discard')
            if length(discard) ~= options.TotSinos*options.Nang*options.Ndist
                error('Error: Size mismatch between sinogram and LORs to be removed')
            end
            if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                discard = discard(1:options.NSinos*options.Nang*options.Ndist);
            end
            ind_apu = uint32(find(discard));
            port = ceil((options.Nang-options.subsets+1)/options.subsets);
            over = options.Nang - port*options.subsets;
            index = cell(options.subsets,1);
            pituus = zeros(options.subsets, 1, 'uint32');
            for i=1:options.subsets
                if over>0
                    index1 = uint32(sort(sub2ind([options.Nang options.Ndist options.NSinos],repmat(repelem(i:options.subsets:(port + 1)*options.subsets,options.Ndist)',options.NSinos,1),repmat((1:options.Ndist)',(port+1)*options.NSinos,1),repelem((1:options.NSinos)',options.Ndist*(port+1),1))));
                    over = over - 1;
                else
                    index1 = uint32(sort(sub2ind([options.Nang options.Ndist options.NSinos],repmat(repelem(i:options.subsets:port*options.subsets,options.Ndist)',options.NSinos,1),repmat((1:options.Ndist)',port*options.NSinos,1),repelem((1:options.NSinos)',options.Ndist*port,1))));
                end
                index{i} = index1(ismember(index1, ind_apu));
                pituus(i) = int32(length(index{i}));
            end
            index = cell2mat(index);
            index = index(ismember(index, ind_apu));
            clear index1 ind_apu
        else
            port = ceil((options.Nang-options.subsets+1)/options.subsets);
            over = options.Nang - port*options.subsets;
            index = cell(options.subsets,1);
            pituus = zeros(options.subsets, 1, 'uint32');
            for i=1:options.subsets
                if over>0
                    index1 = uint32(sort(sub2ind([options.Nang options.Ndist options.NSinos],repmat(repelem(i:options.subsets:(port + 1)*options.subsets,options.Ndist)',options.NSinos,1),repmat((1:options.Ndist)',(port+1)*options.NSinos,1),repelem((1:options.NSinos)',options.Ndist*(port+1),1))));
                    over = over - 1;
                else
                    index1 = uint32(sort(sub2ind([options.Nang options.Ndist options.NSinos],repmat(repelem(i:options.subsets:port*options.subsets,options.Ndist)',options.NSinos,1),repmat((1:options.Ndist)',port*options.NSinos,1),repelem((1:options.NSinos)',options.Ndist*port,1))));
                end
                index{i} = uint32(index1);
                pituus(i) = int32(length(index1));
            end
            clear index1
        end
    elseif options.subsets > 1
        % for raw list-mode data, take the options.subsets randomly
        % last subset has all the spare indices
        if options.precompute_lor || options.reconstruction_method == 3 || options.reconstruction_method == 2
            load([options.machine_name '_detector_locations_' num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '_raw.mat'],'LL','lor')
            indices = uint32(length(LL));
            index = cell(options.subsets, 1);
            port = uint32(floor(length(LL)/options.subsets));
            if options.use_Shuffle
                apu = Shuffle(indices(end), 'index')';
            else
                apu = uint32(randperm(indices(end)))';
            end
            pituus = zeros(options.subsets, 1, 'uint32');
            for i = 1 : options.subsets
                if i == options.subsets
                    index{i} = apu(port*(i-1)+1:end);
                else
                    index{i} = apu(port*(i-1)+1:(port*(i)));
                end
                pituus(i) = int32(length(index{i}));
            end
            clear apu
        else
            load([options.machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'],'LL')
            indices = uint32(length(LL));
            index = cell(options.subsets, 1);
            port = uint32(floor(length(LL)/options.subsets));
            if options.use_Shuffle
                apu = Shuffle(indices(end), 'index')';
            else
                apu = uint32(randperm(indices(end)))';
            end
            for i = 1 : options.subsets
                if i == options.subsets
                    index{i} = apu(port*(i-1)+1:end);
                else
                    index{i} = apu(port*(i-1)+1:(port*(i)));
                end
            end
            clear apu
        end
    end
    
    %%
    for kk = 1 : options.subsets
        A = observation_matrix_formation_nongate(options, kk, index, LL, pituus, lor);
        % Use A here (y = Ax)
        % index-vector contains the measurement numbers used in subset kk
        
    end
    
end

%% Separate reconstruction algorithms

if options.single_reconstructions
    
    % OSEM
    x_osem = OSEM(options);
    
    % RAMLA
    x_ramla = RAMLA(options);
    
    % MRAMLA
    x_mramla = MRAMLA(options);
    
    % COSEM
    x_cosem = COSEM(options);
    
    % ECOSEM
    x_ecosem = ECOSEM(options);
    
    % ACOSEM
    x_acosem = ACOSEM(options);
    
    % OSL MRP
    x_oslmrp = OSL_MRP(options);
    
    % BSREM MRP
    x_bsremmrp = BSREM_MRP(options);
    
    % OSL Quadratic
    x_oslquad = OSL_quad(options);
    
    % BSREM Quadratic
    x_bsremquad = BSREM_quad(options);
    
    % OSL L-filter
    x_oslL = OSL_L(options);
    
    % BSREM L-filter
    x_bsremL = BSREM_L(options);
    
    % OSL FMH
    x_oslFMH = OSL_FMH(options);
    
    % OSL weighted mean
    x_oslwmean = OSL_weightedmean(options);
    
end