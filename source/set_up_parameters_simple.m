function options = set_up_parameters_simple(options)
%set_up_parameters_simple Set up missing parameters to simple
%reconstruction
if ~isfield(options, 'cryst_per_block_axial')
    options.cryst_per_block_axial = options.cryst_per_block;
end
if ~isfield(options, 'transaxial_multip')
    options.transaxial_multip = 1;
end
options.pseudot = [];
options.det_per_ring = options.blocks_per_ring*options.cryst_per_block * options.transaxial_multip;
options.det_w_pseudo = options.blocks_per_ring*(options.cryst_per_block) * options.transaxial_multip;
options.rings = options.linear_multip * options.cryst_per_block_axial;
options.detectors = options.det_per_ring*options.rings;
options.sampling_raw = 1;
options.segment_table = [options.Nz, options.Nz - (options.span + 1):-options.span*2:max(options.Nz - options.ring_difference*2, options.rings - options.ring_difference)];
if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
    options.segment_table = [options.segment_table(1), repeat_elem(options.segment_table(2:end),2,1)];
else
    options.segment_table = [options.segment_table(1), repelem(options.segment_table(2:end),2)];
end
options.TotSinos = sum(options.segment_table);
options.time_index = 0;
options.NSinos = options.TotSinos;
options.ndist_side = -1;
options.sampling = 1;
options.fill_sinogram_gaps = false;
options.gap_filling_method = 'fillmissing';
options.interpolation_method_fillmissing = 'linear';
options.interpolation_method_inpaint = 0;
options.use_raw_data = false;
options.store_raw_data = false;
options.randoms_correction = false;
options.variance_reduction = false;
options.randoms_smoothing = false;
options.scatter_correction = false;
options.normalize_scatter = false;
options.scatter_smoothing = false;
options.subtract_scatter = true;
options.attenuation_correction = false;
options.attenuation_datafile = '';
options.compute_normalization = false;
options.normalization_options = [1 1 1 0];
options.normalization_phantom_radius = inf;
options.normalization_scatter_correction = false;
options.normalization_correction = false;
options.use_user_normalization = false;
options.arc_correction = false;
options.corrections_during_reconstruction = false;
options.tot_time = inf;
options.partitions = 1;
options.start = 0;
options.end = options.tot_time;
options.precompute_obs_matrix = false;
options.precompute_lor = false;
options.precompute_all = false;
options.implementation = 4;
options.use_device = 0;
options.use_64bit_atomics = false;
options.force_build = false;
options.use_CUDA = false;
options.cpu_to_gpu_factor = 0;
options.projector_type = 1;
options.deblurring = false;
options.deblur_iterations = 0;
options.tube_width_xy = options.cr_p;
options.tube_width_z = 0;
options.tube_radius = sqrt(2) * (options.cr_pz / 2);
options.voxel_radius = 1;
options.n_rays_transaxial = 1;
options.n_rays_axial = 1;
options.apply_acceleration = false;
options.save_iter = false;
options.subset_type = 2;
options.n_angles = 2;
options.x0 = ones(options.Nx, options.Ny, options.Nz);
options.epps = 1e-8;
options.use_Shuffle = false;
options.use_fsparse = false;
options.med_no_norm = false;
options.mlem = false;
options.osem = true;
options.mramla = false;
options.ramla = false;
options.rosem = false;
options.rbi = false;
options.drama = false;
options.cosem = false;
options.ecosem = false;
options.acosem = false;
options.OSL_MLEM = false;
options.OSL_OSEM = false;
options.MBSREM = false;
options.BSREM = false;
options.ROSEM_MAP = false;
options.RBI_OSL = false;
options.COSEM_OSL = false;
options.MRP = false;
options.quad = false;
options.Huber = false;
options.L = false;
options.FMH = false;
options.weighted_mean = false;
options.TV = false;
options.AD = false;
options.APLS = false;
options.TGV = false;
options.NLM = false;
options.h = 2;
options.lambda0_mbsrem = 0.2;
options.U = 0;
options.lambda0 = 0.2;
options.lambda0_rosem = 1;
options.beta0_drama = 0.1;
options.beta_drama = 1;
options.alpha_drama = 0.1;
options.Ndx = 1;
options.Ndy = 1;
options.Ndz = 1;
options.beta_mrp_osem = 0.1;
options.beta_mrp_mlem = 1.5;
options.beta_mrp_mbsrem = 0.3;
options.beta_mrp_bsrem = 0.1;
options.beta_mrp_rosem = 2;
options.beta_mrp_rbi = 0.1;
options.beta_mrp_cosem = 1;
options.beta_quad_osem = 0.01;
options.beta_quad_mlem = 0.1;
options.beta_quad_mbsrem = 0.05;
options.beta_quad_bsrem = 0.03;
options.beta_quad_rosem = 0.1;
options.beta_quad_rbi = 0.05;
options.beta_quad_cosem = 0.01;
options.weights = [];
options.beta_huber_osem = 0.01;
options.beta_huber_mlem = 0.1;
options.beta_huber_mbsrem = 0.05;
options.beta_huber_bsrem = 0.03;
options.beta_huber_rosem = 0.1;
options.beta_huber_rbi = 0.05;
options.beta_huber_cosem = 0.01;
options.huber_delta = 5;
options.weights_huber = [];
options.beta_L_osem = 0.1;
options.beta_L_mlem = 0.1;
options.beta_L_mbsrem = 0.1;
options.beta_L_bsrem = 0.03;
options.beta_L_rosem = 3;
options.beta_L_rbi = 0.09;
options.beta_L_cosem = 0.1;
options.a_L = [];
options.oneD_weights = false;
options.beta_fmh_osem = 0.1;
options.beta_fmh_mlem = 0.1;
options.beta_fmh_mbsrem = 0.6;
options.beta_fmh_bsrem = 5;
options.beta_fmh_rosem = 8;
options.beta_fmh_rbi = 0.5;
options.beta_fmh_cosem = 0.1;
options.fmh_weights = [];
options.fmh_center_weight = 4;
options.mean_type = 1;
options.beta_weighted_osem = 0.02;
options.beta_weighted_mlem = 0.1;
options.beta_weighted_mbsrem = 0.1;
options.beta_weighted_bsrem = 5;
options.beta_weighted_rosem = 3;
options.beta_weighted_rbi = 0.04;
options.beta_weighted_cosem = 0.2;
options.weighted_weights = [];
options.weighted_center_weight = 4;
options.beta_TV_osem = 0.03;
options.beta_TV_mlem = 0.1;
options.beta_TV_mbsrem = 0.01;
options.beta_TV_bsrem = 0.05;
options.beta_TV_rosem = 0.07;
options.beta_TV_rbi = 0.002;
options.beta_TV_cosem = 0.003;
options.TVsmoothing = 1e-1;
options.TV_use_anatomical = false;
options.TVtype = 3;
options.TV_reference_image = '';
options.T = 0.1;
options.C = 1;
options.tau = 1e-8;
options.beta_ad_osem = 0.1;
options.beta_ad_mlem = 0.1;
options.beta_ad_mbsrem = 0.3;
options.beta_ad_bsrem = 0.2;
options.beta_ad_rosem = 0.0003;
options.beta_ad_rbi = 0.05;
options.beta_ad_cosem = 0.1;
options.TimeStepAD = 0.0625;
options.KAD = 2;
options.NiterAD = 1;
options.FluxType = 1;
options.DiffusionType = 1;
options.beta_APLS_osem = 0.01;
options.beta_APLS_mlem = 0.1;
options.beta_APLS_mbsrem = 0.1;
options.beta_APLS_bsrem = 0.005;
options.beta_APLS_rosem = 0.1;
options.beta_APLS_rbi = 0.1;
options.beta_APLS_cosem = 0.01;
options.eta = 1e-5;
options.APLSsmoothing = 1e-5;
options.APLS_reference_image = '';
options.beta_TGV_osem = 0.05;
options.beta_TGV_mlem = 0.1;
options.beta_TGV_mbsrem = 0.1;
options.beta_TGV_bsrem = 1;
options.beta_TGV_rosem = 0.25;
options.beta_TGV_rbi = 0.1;
options.beta_TGV_cosem = 0.05;
options.alphaTGV = 2;
options.betaTGV = 1;
options.NiterTGV = 30;
options.beta_NLM_osem = 0.025;
options.beta_NLM_mlem = 0.1;
options.beta_NLM_mbsrem = 0.05;
options.beta_NLM_bsrem = 0.01;
options.beta_NLM_rosem = 0.1;
options.beta_NLM_rbi = 0.01;
options.beta_NLM_cosem = 0.01;
options.sigma = 0.01;
options.Nlx = 1;
options.Nly = 1;
options.Nlz = 0;
options.NLM_use_anatomical = false;
options.NLM_reference_image = '';
options.NLM_MRP = false;
options.simple = true;
end

