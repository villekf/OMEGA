function options = setMissingValues(options)
%SETMISSINGVALUES Sets default values for variables that are missing
%   Utility function

if ~isfield(options, 'CT')
    options.CT = false;
end
if ~isfield(options, 'use_CPU')
    options.use_CPU = false;
end
if ~isfield(options, 'projector_type')
    if options.CT && ~options.use_CPU
        options.projector_type = 4;
    else
        options.projector_type = 11;
    end
end
if ~isfield(options, 'PET')
    options.PET = false;
end
if ~isfield(options, 'SPECT')
    options.SPECT = false;
end
if ~isfield(options, 'useSingles')
    options.useSingles = true;
end
if ~isfield(options, 'largeDim')
    options.largeDim = false;
end
if ~isfield(options, 'loadTOF')
    options.loadTOF = true;
end
if ~isfield(options, 'saveSens')
    options.saveSens = true;
end
if ~isfield(options, 'storeResidual')
    options.storeResidual = false;
end
if ~isfield(options, 'sourceToCRot')
    options.sourceToCRot = 0;
end
if ~isfield(options, 'sourceToDetector')
    options.sourceToDetector = 1;
end
if ~isfield(options,'use_raw_data')
    options.use_raw_data = false;
end
if ~isfield(options,'attenuation_phase')
    options.attenuation_phase = false;
end
if ~isfield(options,'rotateAttImage')
    options.rotateAttImage = 0;
end
if ~isfield(options, 'store_raw_data')
    options.store_raw_data = options.use_raw_data;
end
if ~isfield(options,'cryst_per_block')
    options.cryst_per_block = 0;
end
if ~isfield(options,'linear_multip')
    options.linear_multip = 1;
end
if ~isfield(options, 'dPitch')
    options.dPitch = 0;
end
if ~isfield(options, 'dPitchX') && isfield(options, 'dPitch')
    options.dPitchX = options.dPitch;
end
if ~isfield(options, 'dPitchY') && isfield(options, 'dPitch')
    options.dPitchY = options.dPitch;
end
if ~isfield(options, 'dPitchX')
    options.dPitchX = 0;
end
if ~isfield(options, 'dPitchY')
    if options.dPitchX > 0
        options.dPitchY = options.dPitchX;
    else
        options.dPitchY = 0;
    end
end
if ~isfield(options, 'cr_p')
    if options.dPitchX > 0
        options.cr_p = options.dPitchX;
    else
        options.cr_p = 0;
    end
end
if ~isfield(options, 'cr_pz')
    if options.dPitchY > 0
        options.cr_pz = options.dPitchY;
    else
        if options.cr_p > 0
            options.cr_pz = options.cr_p;
        else
            options.cr_pz = 0;
        end
    end
end
if ~isfield(options, 'axial_fov')
    if options.dPitchY > 0
        options.axial_fov = options.dPitchY;
    else
        options.axial_fov = options.cr_pz;
    end
end
if ~isfield(options, 'custom')
    options.custom = false;
end
if ~isfield(options, 'usingLinearizedData')
    options.usingLinearizedData = false;
end
if ~isfield(options, 'no_data_load')
    options.no_data_load = false;
end
if ~isfield(options, 'TOF_bins') || options.TOF_bins == 0
    options.TOF_bins = 1;
end
if ~isfield(options, 'TOF_bins_used') || options.TOF_bins_used == 0
    if options.TOF_bins > 1
        options.TOF_bins_used = options.TOF_bins;
    else
        options.TOF_bins_used = 1;
    end
end
if ~isfield(options, 'TOF_FWHM')
    options.TOF_FWHM = 0;
end
if ~isfield(options,'TOF_width')
    options.TOF_width = 0;
end
if ~isfield(options, 'cryst_per_block_axial')
    options.cryst_per_block_axial = options.cryst_per_block;
end
if ~isfield(options, 'transaxial_multip')
    options.transaxial_multip = 1;
end
if ~isfield(options,'use_machine')
    options.use_machine = 0;
end
if ~isfield(options,'machine_name')
    options.machine_name = '';
end
if ~isfield(options,'use_LMF')
    options.use_LMF = false;
end
if ~isfield(options, 'use_ASCII')
    options = set_GATE_variables(options);
end
if ~isfield(options, 'only_sinos')
    options.only_sinos = false;
end
if ~isfield(options, 'only_reconstructions')
    options.only_reconstructions = false;
end
if ~isfield(options, 'custom')
    options.custom = false;
end
recA = recNames(0);
for kk = 1 : numel(recA)
    if ~isfield(options, recA{kk})
        options.(recA{kk}) = false;
    end
end
recP = recNames(1);
for kk = 1 : numel(recP)
    if ~isfield(options, recP{kk})
        options.(recP{kk}) = false;
    end
end
if isfield(options, 'verticalOffset')
    options.sourceOffsetCol = options.verticalOffset;
end
if ~isfield(options, 'sourceOffsetCol')
    options.sourceOffsetCol = 0;
end
if isfield(options, 'horizontalOffset')
    options.sourceOffsetRow = options.horizontalOffset;
end
if ~isfield(options, 'sourceOffsetRow')
    options.sourceOffsetRow = 0;
end
if ~isfield(options, 'horizontalOffset')
    options.horizontalOffset = options.sourceOffsetRow;
end
if ~isfield(options, 'verticalOffset')
    options.verticalOffset = options.sourceOffsetCol;
end
if ~isfield(options, 'subsets')
    options.subsets = 1;
end
if ~isfield(options, 'subset_type')
    if isfield(options, 'subsetType')
        options.subset_type = options.subsetType;
    else
        options.subset_type = 8;
    end
end
if ~isfield(options, 'useMaskFP')
    options.useMaskFP = false;
end
if ~isfield(options, 'useMaskBP')
    options.useMaskBP = false;
end
if ~isfield(options, 'usePriorMask')
    options.usePriorMask = false;
end
if ~isfield(options, 'offsetCorrection')
    options.offsetCorrection = false;
end
if ~isfield(options,'bedOffset')
    if ~options.SPECT
        options.bedOffset = [];
    else
        options.bedOffset = 0;
    end
end
if isfield(options,'uCenter')
    options.detOffsetRow = options.uCenter;
end
if isfield(options,'vCenter')
    options.detOffsetCol = options.vCenter;
end
if ~isfield(options,'detOffsetRow')
    options.detOffsetRow = [];
end
if ~isfield(options,'detOffsetCol')
    options.detOffsetCol = [];
end
if ~isfield(options,'uCenter')
    options.detOffsetRow = options.detOffsetRow;
end
if ~isfield(options,'vCenter')
    options.detOffsetCol = options.detOffsetCol;
end
if ~isfield(options, 'pitchRoll')
    options.pitchRoll = [];
end
if ~isfield(options,'nBed')
    options.nBed = 1;
end
if ~isfield(options,'flip_image')
    options.flip_image = false;
end
if ~isfield(options, 'flipImageX')
    options.flipImageX = false;
end
if ~isfield(options, 'flipImageY')
    options.flipImageY = false;
end
if ~isfield(options, 'flipImageZ')
    options.flipImageZ = false;
end
if ~isfield(options,'offangle')
    options.offangle = 0;
end
if ~isfield(options, 'oOffsetX')
    options.oOffsetX = 0;
end
if ~isfield(options, 'oOffsetY')
    options.oOffsetY = 0;
end
if ~isfield(options, 'oOffsetZ')
    options.oOffsetZ = 0;
end
if ~isfield(options, 'tube_width_z')
    options.tube_width_z = 0;
end
if ~isfield(options, 'tube_width_xy')
    options.tube_width_xy = 0;
end
if ~isfield(options, 'use_psf')
    options.use_psf = false;
end
if ~isfield(options, 'FWHM')
    options.FWHM = [options.dPitchX options.dPitchX options.dPitchY];
end
if ~isfield(options, 'save_iter')
    options.save_iter = false;
end
if ~isfield(options, 'apply_acceleration')
    options.apply_acceleration = false;
end
if ~isfield(options, 'deblurring')
    options.deblurring = false;
end
if ~isfield(options,'use_device')
    options.use_device = uint32(0);
end
if ~isfield(options,'implementation')
    options.implementation = 2;
end
if ~isfield(options,'platform')
    options.platform = 0;
end
if ~isfield(options, 'use_64bit_atomics')
    options.use_64bit_atomics = ~ismac;
end
if ~isfield(options, 'use_CUDA')
    if ismac
        options.use_CUDA = false;
    else
        options.use_CUDA = checkCUDA(options.use_device);
    end
end
if ~isfield(options, 'nRays')
    options.nRays = 1;
end
if ~isfield(options, 'n_rays_transaxial')
    if options.SPECT
        options.n_rays_transaxial = options.nRays;
    else
        options.n_rays_transaxial = 1;
    end
end
if ~isfield(options, 'n_rays_axial')
    options.n_rays_axial = 1;
end
if ~isfield(options, 'cpu_to_gpu_factor')
    options.cpu_to_gpu_factor = 1;
end
if ~isfield(options,'meanFP')
    options.meanFP = false;
end
if ~isfield(options,'meanBP')
    options.meanBP = false;
end
if ~isfield(options,'useFDKWeights')
    options.useFDKWeights = true;
end
if ~isfield(options, 'Niter')
    options.Niter = 4;
end
if ~isfield(options, 'filteringIterations') || isempty(options.filteringIterations)
    options.filteringIterations = 0;
end
if ~isfield(options, 'saveNIter') || options.save_iter
    options.saveNIter = uint32([]);
end
if options.Niter - 1 <= max(options.saveNIter(:))
    loc = find(options.saveNIter >= options.Niter - 1, 1, 'first');
    if loc > 1
        options.saveNIter = options.saveNIter(1 : loc - 1);
    else
        options.saveNIter = uint32([]);
    end
end
if ~isa(options.saveNIter,'uint32')
    options.saveNIter = uint32(options.saveNIter);
end
if ~isfield(options,'multiResolutionScale')
    options.multiResolutionScale = .25;
end
if ~isfield(options,'nLayers')
    options.nLayers = 1;
end
if ~isfield(options,'tauCP') || (isfield(options,'tauCP') && isempty(options.tauCP))
    options.tauCP = 0;
end
if ~isfield(options,'thetaCP') || (isfield(options,'thetaCP') && isempty(options.thetaCP))
    options.thetaCP = 1;
end
if ~isfield(options,'sigmaCP') || (isfield(options,'sigmaCP') && isempty(options.sigmaCP))
    options.sigmaCP = 1;
end
if ~isfield(options, 'sigma2CP')
    options.sigma2CP = options.sigmaCP;
end
if ~isfield(options, 'tauCPFilt') || isempty(options.tauCPFilt)
    options.tauCPFilt = 0;
end
if ~isfield(options,'powerIterations') || (isfield(options,'powerIterations') && isempty(options.powerIterations))
    options.powerIterations = 20;
end
if ~isfield(options,'derivativeType')
    options.derivativeType = 0;
end
if ~isfield(options,'enforcePositivity')
    options.enforcePositivity = true;
end
if ~isfield(options,'precondTypeImage')
    options.precondTypeImage = [false;false;false;false;false;false;false];
end
if numel(options.precondTypeImage) < 7
    options.precondTypeImage = [options.precondTypeImage;false(7-numel(options.precondTypeImage),1)];
end
if ~isfield(options,'precondTypeMeas')
    options.precondTypeMeas = [false;false];
end
if numel(options.precondTypeMeas) < 2
    options.precondTypeMeas = [options.precondTypeMeas;false(2-numel(options.precondTypeMeas),1)];
end
if ~isfield(options, 'gradV1')
    options.gradV1 = 0.5;
end
if ~isfield(options, 'gradV2')
    options.gradV1 = 2.5;
end
if ~isfield(options, 'gradInitIter')
    options.gradInitIter = options.subsets;
end
if ~isfield(options, 'gradLastIter')
    options.gradLastIter = options.gradInitIter;
end
if ~isfield(options, 'filterWindow')
    options.filterWindow = 'hamming';
end
if ~isfield(options, 'cutoffFrequency')
    options.cutoffFrequency = 1;
end
if ~isfield(options, 'normalFilterSigma')
    options.normalFilterSigma = 0.25;
end
if ~isfield(options, 'useZeroPadding')
    options.useZeroPadding = false;
end
if ~isfield(options, 'use_binary')
    options.use_binary = false;
end
if ~isfield(options, 'pitch')
    options.pitch = false;
end
if ~isfield(options,'TOF_bins')
    options.TOF_bins = 1;
end
if ~isfield(options,'TOF_width')
    options.TOF_width = 0;
end
if ~isfield(options,'compute_sensitivity_image')
    options.compute_sensitivity_image = false;
end
if ~isfield(options,'listmode')
    options.listmode = false;
end
if ~isfield(options,'useIndexBasedReconstruction')
    options.useIndexBasedReconstruction = false;
end
if ~isfield(options,'sampling_raw')
    options.sampling_raw = 1;
end
if ~isfield(options,'nProjections')
    options.nProjections = 1;
end
if ~isfield(options,'oOffsetX')
    options.oOffsetX = 0;
end
if ~isfield(options,'oOffsetY')
    options.oOffsetY = 0;
end
if ~isfield(options,'oOffsetZ')
    options.oOffsetZ = 0;
end
if ~isfield(options, 'FOVa_y')
    if isfield(options, 'FOVa_x')
        options.FOVa_y = options.FOVa_x;
    else
        options.FOVa_y = 0;
    end
end
if ~isfield(options, 'FOVa_x')
    options.FOVa_x = 0;
end
if ~isfield(options, 'Nz')
    options.Nz = 1;
end
if ~isfield(options, 'x0')
    options.x0 = ones(options.Nx, options.Ny, options.Nz) * 1e-5;
end
if ~isfield(options, 'dL')
    options.dL = 0;
end
if ~isfield(options, 'epps')
    options.epps = 1e-5;
end
if ~isfield(options, 'use_Shuffle')
    options.use_Shuffle = false;
end
if ~isfield(options, 'use_fsparse')
    options.use_fsparse = false;
end
if ~isfield(options, 'med_no_norm')
    options.med_no_norm = false;
end
if ~isfield(options, 'errorChecking')
    options.errorChecking = false;
end
if ~isfield(options, 'blocks_per_ring')
    options.blocks_per_ring = 1;
end
if ~isfield(options, 'diameter')
    options.diameter = 1;
end
if ~isfield(options, 'span')
    options.span = 3;
end
if ~isfield(options, 'binning')
    options.binning = 1;
end
if ~isfield(options, 'tube_radius')
    options.tube_radius = sqrt(2) * (options.cr_pz / 2);
end
if ~isfield(options, 'voxel_radius')
    options.voxel_radius = 1;
end
if ~isfield(options, 'use_32bit_atomics')
    options.use_32bit_atomics = false;
end
if ~isfield(options, 'precompute')
    options.precompute = false;
end
if ~isfield(options, 'precompute_obs_matrix')
    options.precompute_obs_matrix = false;
end
if ~isfield(options, 'precompute_lor')
    options.precompute_lor = false;
end
if ~isfield(options, 'precompute_all')
    options.precompute_all = false;
end
if ~isfield(options, 'verbose')
    options.verbose = 1;
end
if ~isfield(options, 'tot_time')
    options.tot_time = inf;
end
if ~isfield(options, 'partitions')
    options.partitions = 1;
end
if ~isfield(options, 'start')
    options.start = 0;
end
if ~isfield(options, 'end')
    options.end = options.tot_time;
end
if ~isfield(options, 'randoms_correction')
    options.randoms_correction = false;
end
if ~isfield(options, 'variance_reduction')
    options.variance_reduction = false;
end
if ~isfield(options, 'randoms_smoothing')
    options.randoms_smoothing = false;
end
if ~isfield(options, 'scatter_correction')
    options.scatter_correction = false;
end
if ~isfield(options, 'scatter_variance_reduction')
    options.scatter_variance_reduction = false;
end
if ~isfield(options, 'normalize_scatter')
    options.normalize_scatter = false;
end
if ~isfield(options, 'scatter_smoothing')
    options.scatter_smoothing = false;
end
if ~isfield(options, 'subtract_scatter')
    options.subtract_scatter = true;
end
if ~isfield(options, 'scatter')
    options.scatter = false;
end
if ~isfield(options, 'additionalCorrection')
    options.additionalCorrection = false;
end
if ~isfield(options, 'ScatterC')
    options.ScatterC = [];
end
if ~isfield(options, 'attenuation_correction')
    options.attenuation_correction = false;
end
if ~isfield(options, 'CT_attenuation')
    options.CT_attenuation = true;
end
if ~isfield(options, 'dualLayerSubmodule')
    options.dualLayerSubmodule = false;
end
if ~isfield(options, 'attenuation_datafile')
    options.attenuation_datafile = '';
end
if ~isfield(options, 'vaimennus')
    if options.useSingles
        options.vaimennus = single([]);
    else
        options.vaimennus = [];
    end
end
if ~isfield(options, 'SinM')
    options.SinM = [];
end
if ~isfield(options, 'SinDelayed')
    options.SinDelayed = [];
end
if ~isfield(options, 'compute_normalization')
    options.compute_normalization = false;
end
if ~isfield(options, 'normalization_options')
    options.normalization_options = [1 1 1 1];
end
if ~isfield(options, 'normalization_phantom_radius')
    options.normalization_phantom_radius = inf;
end
if ~isfield(options, 'normalization_attenuation')
    options.normalization_attenuation = [];
end
if ~isfield(options, 'normalization_scatter_correction')
    options.normalization_scatter_correction = false;
end
if ~isfield(options, 'normalization_correction')
    options.normalization_correction = false;
end
if ~isfield(options, 'use_user_normalization')
    options.use_user_normalization = false;
end
if ~isfield(options, 'normalization')
    if options.useSingles
        options.normalization = single([]);
    else
        options.normalization = [];
    end
end
if ~isfield(options, 'arc_correction')
    options.arc_correction = false;
end
if ~isfield(options, 'arc_interpolation')
    options.arc_interpolation = 'linear';
end
if ~isfield(options, 'global_correction_factor')
    options.global_correction_factor = [];
end
if ~isfield(options, 'corrections_during_reconstruction')
    options.corrections_during_reconstruction = true;
end
if ~isfield(options, 'ordinaryPoisson')
    options.ordinaryPoisson = options.corrections_during_reconstruction;
end
if ~isfield(options, 'ndist_side')
    options.ndist_side = 1;
end
if ~isfield(options, 'sampling')
    options.sampling = 1;
end
if ~isfield(options, 'sampling_interpolation_method')
    options.sampling_interpolation_method = 'linear';
end
if ~isfield(options, 'fill_sinogram_gaps')
    options.fill_sinogram_gaps = false;
end
if ~isfield(options, 'obtain_trues')
    options.obtain_trues = false;
end
if ~isfield(options, 'reconstruct_trues')
    options.reconstruct_trues = false;
end
if ~isfield(options, 'store_scatter')
    options.store_scatter = false;
end
if ~isfield(options, 'scatter_components')
    options.scatter_components = [1 1 0 0];
end
if ~isfield(options, 'reconstruct_scatter')
    options.reconstruct_scatter = false;
end
if ~isfield(options, 'store_randoms')
    options.store_randoms = false;
end
if ~isfield(options, 'source')
    options.source = false;
end
if ~isfield(options, 'pseudot')
    options.pseudot = [];
end
if ~isfield(options,'Ndist')
    options.Ndist = 400;
end
if ~isfield(options, 'rings')
    if isfield(options, 'cryst_per_block')
        if isfield(options, 'linear_multip') && isfield(options, 'cryst_per_block')
            options.rings = options.linear_multip * options.cryst_per_block(end);
        elseif isfield(options, 'cryst_per_block')
            options.rings = options.cryst_per_block(end);
        end
    else
        options.rings = 1;
    end
end
if ~isfield(options,'ring_difference')
    options.ring_difference = options.rings - 1;
end
if ~isfield(options, 'ring_difference_raw')
    options.ring_difference_raw = options.rings - 1;
end
if ~isfield(options, 'det_per_ring')
    if isfield(options, 'blocks_per_ring') && isfield(options, 'cryst_per_block')
        options.det_per_ring = options.blocks_per_ring*options.cryst_per_block(end);
    else
        options.det_per_ring = options.Nang * 2 * options.Ndist;
    end
end
if ~isfield(options, 'det_w_pseudo')
    options.det_w_pseudo = options.det_per_ring(end);
end
if ~isfield(options, 'detectors')
    options.detectors = options.det_w_pseudo*options.rings;
end
if ~isfield(options,'Nang')
    if isfield(options,'det_w_pseudo') && options.det_w_pseudo > 0 
        options.Nang = options.det_w_pseudo / 2;
    else
        options.Nang = 1;
    end
end
if ~isfield(options, 'segment_table') && isfield(options, 'span') && options.span > 0
    if options.span == 1
        options.segment_table = options.rings^2;
    else
        options.segment_table = [options.Nz, options.Nz - (options.span + 1):-options.span*2:max(options.Nz - options.ring_difference*2, options.rings - options.ring_difference)];
        if exist('OCTAVE_VERSION','builtin') == 0 && exist('repelem', 'builtin') == 0
            options.segment_table = [options.segment_table(1), repeat_elem(options.segment_table(2:end),2,1)];
        else
            options.segment_table = [options.segment_table(1), repelem(options.segment_table(2:end),2)];
        end
    end
end
if ~isfield(options,'NSinos')
    if isfield(options, 'TotSinos')
        options.NSinos = options.TotSinos;
    elseif isfield(options, 'segment_table')
        options.NSinos = sum(options.segment_table);
    else
        options.NSinos = 1;
    end
end
if ~isfield(options,'TotSinos')
    options.TotSinos = options.NSinos;
end
if ~isfield(options, 'h')
    options.h = 2;
end
if ~isfield(options, 'U')
    options.U = 10000;
end
if ~isfield(options, 'Ndx')
    options.Ndx = 1;
end
if ~isfield(options, 'Ndy')
    options.Ndy = 1;
end
if ~isfield(options, 'Ndz')
    options.Ndz = 1;
end
if ~isfield(options, 'weights')
    options.weights = [];
end
if ~isfield(options, 'weights_huber')
    options.weights_huber = [];
end
if ~isfield(options, 'a_L')
    options.a_L = [];
end
if ~isfield(options, 'oneD_weights')
    options.oneD_weights = false;
end
if ~isfield(options, 'fmh_weights')
    options.fmh_weights = [];
end
if ~isfield(options, 'fmh_center_weight')
    options.fmh_center_weight = 4;
end
if ~isfield(options, 'mean_type')
    options.mean_type = 4;
end
if ~isfield(options, 'weighted_weights')
    options.weighted_weights = [];
end
if ~isfield(options, 'weighted_center_weight')
    options.weighted_center_weight = 2;
end
if ~isfield(options, 'TVsmoothing')
    options.TVsmoothing = 1e-2;
end
if ~isfield(options, 'tau')
    options.tau = 1e-8;
end
if ~isfield(options, 'TV_use_anatomical')
    options.TV_use_anatomical = false;
end
if ~isfield(options, 'TVtype')
    options.TVtype = 1;
end
if ~isfield(options, 'T')
    options.T = 0.01;
end
if ~isfield(options, 'C')
    options.C = 1;
end
if ~isfield(options, 'FluxType')
    options.FluxType = 1;
end
if ~isfield(options, 'DiffusionType')
    options.DiffusionType = 2;
end
if ~isfield(options, 'eta')
    options.eta = 1e-5;
end
if ~isfield(options, 'APLSsmoothing')
    options.APLSsmoothing = 1e-5;
end
if ~isfield(options, 'Nlx')
    options.Nlx = 1;
end
if ~isfield(options, 'Nly')
    options.Nly = 1;
end
if ~isfield(options, 'Nlz')
    options.Nlz = 1;
end
if ~isfield(options, 'RDP_gamma')
    options.RDP_gamma = 1;
end
if ~isfield(options, 'RDPIncludeCorners')
    options.RDPIncludeCorners = false;
end
if ~isfield(options, 'RDP_use_anatomical')
    options.RDP_use_anatomical = false;
end
if ~isfield(options, 'NLM_use_anatomical')
    options.NLM_use_anatomical = false;
end
if ~isfield(options, 'NLTV')
    options.NLTV = false;
end
if ~isfield(options, 'NLRD')
    options.NLRD = false;
end
if ~isfield(options, 'NLLange')
    options.NLLange = false;
end
if ~isfield(options, 'NLGGMRF')
    options.NLGGMRF = false;
end
if ~isfield(options, 'NLM_MRP')
    options.NLM_MRP = false;
end
if ~isfield(options, 'NLAdaptive')
    options.NLAdaptive = false;
end
if ~isfield(options, 'NLAdaptiveConstant')
    options.NLAdaptiveConstant = 1e-5;
end
if ~isfield(options, 'sigma')
    options.sigma = 1;
end
if ~isfield(options, 'weights_RDP')
    options.weights_RDP = [];
end
if ~isfield(options, 'name')
    options.name = [];
end
if ~isfield(options, 'segment_table')
    options.segment_table = [];
end
if ~isfield(options, 'lor_a')
    options.lor_a = [];
end
if ~isfield(options, 'lor_orth')
    options.lor_orth = [];
end
if ~isfield(options, 'alpha0TGV')
    options.alpha0TGV = 0;
end
if ~isfield(options, 'alpha1TGV')
    options.alpha1TGV = 0;
end
if ~isfield(options, 'useL2Ball')
    options.useL2Ball = true;
end
if ~isfield(options, 'useMAD')
    options.useMAD = true;
end
if ~isfield(options, 'useImages')
    if options.projector_type == 6
        options.useImages = false;
    else
        options.useImages = true;
    end
end
if ~isfield(options, 'useEFOV')
    options.useEFOV = false;
end
if ~isfield(options, 'useExtrapolation')
    options.useExtrapolation = false;
end
if ~isfield(options, 'flat')
    options.flat = 0;
end
if ~isfield(options, 'eFOVIndices')
    options.eFOVIndices = [];
end
if ~isfield(options, 'NxOrig')
    options.NxOrig = options.Nx;
end
if ~isfield(options, 'NyOrig')
    options.NyOrig = options.Ny;
end
if ~isfield(options, 'NzOrig')
    options.NzOrig = options.Nz;
end
if ~isfield(options, 'NxPrior')
    options.NxPrior = options.NxOrig;
end
if ~isfield(options, 'NyPrior')
    options.NyPrior = options.NyOrig;
end
if ~isfield(options, 'NzPrior')
    options.NzPrior = options.NzOrig;
end
if ~isfield(options, 'use2DTGV')
    options.use2DTGV = false;
end
if ~isfield(options, 'rho_PKMA')
    options.rho_PKMA = .45;
end
if ~isfield(options, 'delta_PKMA')
    options.delta_PKMA = 100;
end
if ~isfield(options, 'useMultiResolutionVolumes')
    options.useMultiResolutionVolumes = false;
end
if ~isfield(options, 'nMultiVolumes')
    options.nMultiVolumes = 0;
end
if ~isfield(options,'relaxationScaling')
    options.relaxationScaling = false;
end
if ~isfield(options,'computeRelaxationParameters')
    options.computeRelaxationParameters = false;
end
if ~isfield(options, 'lambda')
    options.lambda = [];
end
if ~isfield(options, 'PDAdaptiveType')
    options.PDAdaptiveType = 0;
end
if ~isfield(options, 'storeFP')
    options.storeFP = false;
end
if isfield(options,'xSize') && ~isfield(options,'nRowsD')
    options.nRowsD = options.xSize;
end
if isfield(options,'ySize') && ~isfield(options,'nColsD')
    options.nColsD = options.ySize;
end
if isfield(options,'Ndist') && ~isfield(options,'nRowsD')
    options.nRowsD = options.Ndist;
end
if isfield(options,'Nang') && ~isfield(options,'nColsD')
    options.nColsD = options.Nang;
end
if ~isfield(options,'nColsD')
    options.nColsD = 1;
end
if ~isfield(options, 'derivType')
    options.derivType = 0;
end
if ~isfield(options, 'Nf')
    options.Nf = 0;
end
if ~isfield(options,'hyperbolicDelta')
    options.hyperbolicDelta = 1;
end
if ~isfield(options,'POCS_NgradIter')
    options.POCS_NgradIter = 20;
end
if ~isfield(options,'POCS_alpha')
    options.POCS_alpha = 0.2;
end
if ~isfield(options,'POCS_rMax')
    options.POCS_rMax = .95;
end
if ~isfield(options,'POCS_alphaRed')
    options.POCS_alphaRed = .95;
end
if ~isfield(options,'POCSepps')
    options.POCSepps = 1e-4;
end
if ~isfield(options, 'FISTA_acceleration')
    options.FISTA_acceleration = false;
end
if ~isfield(options, 'FISTAType')
    options.FISTAType = 0;
end
if ~isfield(options,'maskFPZ')
    options.maskFPZ = 1;
end
if ~isfield(options,'maskBPZ')
    options.maskBPZ = 1;
end
if ~isfield(options, 'stochasticSubsetSelection')
    options.stochasticSubsetSelection = false;
end
if ~isfield(options, 'seed')
    options.seed = -1;
end
if ~isfield(options, 'useParallelBeam')
    options.useParallelBeam = false;
end
if ~isfield(options, 'useTotLength')
    options.useTotLength = true;
end
if ~isfield(options, 'coneOfResponseStdCoeffA')
    options.coneOfResponseStdCoeffA = 0.1;
end
if ~isfield(options, 'coneOfResponseStdCoeffB')
    options.coneOfResponseStdCoeffB = 0.1;
end
if ~isfield(options, 'coneOfResponseStdCoeffC')
    options.coneOfResponseStdCoeffC = 0.1;
end
if ~isfield(options, 'colL') % SPECT collimator hole length
    options.colL = 1;
end
if ~isfield(options, 'colD') % SPECT collimator distance from detector surface
    options.colD = 1;
end
if ~isfield(options, 'colR') % SPECT collimator hole radius
    options.colR = 1;
end
if ~isfield(options, 'colFxy') % SPECT collimator focal distance
    options.colFxy = Inf; % Parallel-hole
end
if ~isfield(options, 'colFz') % SPECT collimator focal distance
    options.colFz = Inf; % Parallel-hole
end
if ~isfield(options, 'crXY') % SPECT detector crystal size
    options.crXY = 1;
end
if ~isfield(options, 'rayShiftsDetector')
    options.rayShiftsDetector = [];
end
if ~isfield(options, 'rayShiftsSource')
    options.rayShiftsSource = [];
end
if ~isfield(options, 'iR')
    options.iR = 1;
end
if ~isfield(options, 'CORtoDetectorSurface') % SPECT detector surface distance from detector centre of rotation
    options.CORtoDetectorSurface = 0;
end
if ~isfield(options, 'swivelAngles')
    options.swivelAngles = [];
end
if ~isfield(options, 'useHalf') % When true, computations are performed with FP16 precision instead of FP32
    options.useHalf = 0;
end