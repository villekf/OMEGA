function options = init_next_iter(options, iter)
% Initialize the next iteration in the custom prior reconstruction

if options.osem
    options.im_vectors.OSEM(:, iter + 1) = options.im_vectors.OSEM_apu;
end

if options.mlem
    options.im_vectors.MLEM(:, iter + 1) = options.im_vectors.MLEM_apu;
end

if options.mramla
    options.im_vectors.MRAMLA(:, iter + 1) = options.im_vectors.MRAMLA_apu;
end

if options.ramla
    options.im_vectors.RAMLA(:, iter + 1) = options.im_vectors.RAMLA_apu;
end

if options.rosem
    options.im_vectors.ROSEM(:, iter + 1) = options.im_vectors.ROSEM_apu;
end

if options.rbi
    options.im_vectors.RBI(:, iter + 1) = options.im_vectors.RBI_apu;
end

if options.drama
    options.im_vectors.DRAMA(:, iter + 1) = options.im_vectors.DRAMA_apu;
end

if options.cosem
    options.im_vectors.COSEM(:, iter + 1) = options.im_vectors.COSEM_apu;
end

if options.ecosem
    options.im_vectors.ECOSEM(:, iter + 1) = options.im_vectors.ECOSEM_apu;
end

if options.acosem
    options.im_vectors.ACOSEM(:, iter + 1) = options.im_vectors.ACOSEM_apu;
end

if options.MRP && options.OSL_OSEM
    options.im_vectors.MRP_OSL(:, iter + 1) = options.im_vectors.MRP_OSL_apu;
end
if options.MRP && options.OSL_MLEM
    options.im_vectors.MRP_MLEM(:, iter + 1) = options.im_vectors.MRP_MLEM_apu;
end
if options.MRP && options.MBSREM
    options.im_vectors.MRP_MBSREM(:, iter + 1) = options.im_vectors.MRP_MBSREM_apu;
end

if options.MRP && options.BSREM
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    med = MRP(options.im_vectors.MRP_BSREM_apu, options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
    options.im_vectors.MRP_BSREM(:,iter+1) = BSREM_iter(options.im_vectors.MRP_BSREM_apu, options.lam, iter, options.beta_mrp_bsrem, med, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['BSREM MRP iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['BSREM MRP iteration ' num2str(iter) ' finished'])
    end
    else
    end
end

if options.MRP && options.ROSEM_MAP
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    med = MRP(options.im_vectors.MRP_ROSEM_apu, options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
    options.im_vectors.MRP_ROSEM(:,iter+1) = BSREM_iter(options.im_vectors.MRP_ROSEM_apu, options.lam_rosem, iter, options.beta_mrp_rosem, med, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['ROSEM MRP iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['ROSEM MRP iteration ' num2str(iter) ' finished'])
    end
    else
    end
end

if options.MRP && options.RBI_MAP
    options.im_vectors.MRP_RBI(:, iter + 1) = options.im_vectors.MRP_RBI_apu;
end

if options.MRP && any(options.COSEM_MAP)
    options.im_vectors.MRP_COSEM(:, iter + 1) = options.im_vectors.MRP_COSEM_apu;
end

if options.quad && options.OSL_OSEM
    options.im_vectors.Quad_OSL(:, iter + 1) = options.im_vectors.Quad_OSL_apu;
end
if options.quad && options.OSL_MLEM
    options.im_vectors.Quad_MLEM(:, iter + 1) = options.im_vectors.Quad_MLEM_apu;
end
if options.quad && options.MBSREM
    options.im_vectors.Quad_MBSREM(:, iter + 1) = options.im_vectors.Quad_MBSREM_apu;
end

if options.quad && options.BSREM
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    med = Quadratic_prior(options.im_vectors.Quad_BSREM_apu, options.tr_offsets, options.weights, options.weights_quad, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz);
    options.im_vectors.Quad_BSREM(:,iter+1) = BSREM_iter(options.im_vectors.Quad_BSREM_apu, options.lam, iter, options.beta_quad_bsrem, med, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['BSREM quadratic iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['BSREM quadratic iteration ' num2str(iter) ' finished'])
    end
    else
    end
end

if options.quad && options.ROSEM_MAP
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    med = Quadratic_prior(options.im_vectors.Quad_ROSEM_apu, options.tr_offsets, options.weights, options.weights_quad, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz);
    options.im_vectors.Quad_ROSEM(:,iter+1) = BSREM_iter(options.im_vectors.Quad_ROSEM_apu, options.lam_rosem, iter, options.beta_quad_rosem, med, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['ROSEM quadratic iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['ROSEM quadratic iteration ' num2str(iter) ' finished'])
    end
    else
    end
end

if options.quad && options.RBI_MAP
    options.im_vectors.Quad_RBI(:, iter + 1) = options.im_vectors.Quad_RBI_apu;
end

if options.quad && any(options.COSEM_MAP)
    options.im_vectors.Quad_COSEM(:, iter + 1) = options.im_vectors.Quad_COSEM_apu;
end

if options.L && options.OSL_OSEM
    options.im_vectors.L_OSL(:, iter + 1) = options.im_vectors.L_OSL_apu;
end
if options.L && options.OSL_MLEM
    options.im_vectors.L_MLEM(:, iter + 1) = options.im_vectors.L_MLEM_apu;
end
if options.L && options.MBSREM
    options.im_vectors.L_MBSREM(:, iter + 1) = options.im_vectors.L_MBSREM_apu;
end

if options.L && options.BSREM
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    med = L_filter(options.im_vectors.L_BSREM_apu, options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, options.epps, options.med_no_norm);
    options.im_vectors.L_BSREM(:,iter+1) = BSREM_iter(options.im_vectors.L_BSREM_apu, options.lam, iter, options.beta_L_bsrem, med, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['BSREM L-filter iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['BSREM L-filter iteration ' num2str(iter) ' finished'])
    end
    else
    end
end

if options.L && options.ROSEM_MAP
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    med = L_filter(options.im_vectors.L_ROSEM_apu, options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, options.epps, options.med_no_norm);
    options.im_vectors.L_ROSEM(:,iter+1) = BSREM_iter(options.im_vectors.L_ROSEM_apu, options.lam_rosem, iter, options.beta_L_rosem, med, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['ROSEM L-filter iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['ROSEM L-filter iteration ' num2str(iter) ' finished'])
    end
    else
    end
end

if options.L && options.RBI_MAP
    options.im_vectors.L_RBI(:, iter + 1) = options.im_vectors.L_RBI_apu;
end

if options.L && any(options.COSEM_MAP)
    options.im_vectors.L_COSEM(:, iter + 1) = options.im_vectors.L_COSEM_apu;
end

if options.FMH && options.OSL_OSEM
    options.im_vectors.FMH_OSL(:, iter + 1) = options.im_vectors.FMH_OSL_apu;
end
if options.FMH && options.OSL_MLEM
    options.im_vectors.FMH_MLEM(:, iter + 1) = options.im_vectors.FMH_MLEM_apu;
end
if options.FMH && options.MBSREM
    options.im_vectors.FMH_MBSREM(:, iter + 1) = options.im_vectors.FMH_MBSREM_apu;
end

if options.FMH && options.BSREM
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    med = FMH(options.im_vectors.FMH_BSREM_apu, options.tr_offsets, options.fmh_weights, options.weights, options.Nx, options.Ny, options.Nz, N, options.Ndx, options.Ndy, options.Ndz, options.epps, options.med_no_norm);
    options.im_vectors.FMH_BSREM(:,iter+1) = BSREM_iter(options.im_vectors.FMH_BSREM_apu, options.lam, iter, options.beta_fmh_bsrem, med, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['BSREM FMH iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['BSREM FMH iteration ' num2str(iter) ' finished'])
    end
    else
    end
end

if options.FMH && options.ROSEM_MAP
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    med = FMH(options.im_vectors.FMH_ROSEM_apu, options.tr_offsets, options.fmh_weights, options.weights, options.Nx, options.Ny, options.Nz, N, options.Ndx, options.Ndy, options.Ndz, options.epps, options.med_no_norm);
    options.im_vectors.FMH_ROSEM(:,iter+1) = BSREM_iter(options.im_vectors.FMH_ROSEM_apu, options.lam_rosem, iter, options.beta_fmh_rosem, med, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['ROSEM FMH iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['ROSEM FMH iteration ' num2str(iter) ' finished'])
    end
    else
    end
end

if options.FMH && options.RBI_MAP
    options.im_vectors.FMH_RBI(:, iter + 1) = options.im_vectors.FMH_RBI_apu;
end

if options.FMH && any(options.COSEM_MAP)
    options.im_vectors.FMH_COSEM(:, iter + 1) = options.im_vectors.FMH_COSEM_apu;
end

if options.weighted_mean && options.OSL_OSEM
    options.im_vectors.Weighted_OSL(:, iter + 1) = options.im_vectors.Weighted_OSL_apu;
end
if options.weighted_mean && options.OSL_MLEM
    options.im_vectors.Weighted_MLEM(:, iter + 1) = options.im_vectors.Weighted_MLEM_apu;
end
if options.weighted_mean && options.MBSREM
    options.im_vectors.Weighted_MBSREM(:, iter + 1) = options.im_vectors.Weighted_MBSREM_apu;
end

if options.weighted_mean && options.BSREM
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    med = Weighted_mean(options.im_vectors.Weighted_BSREM_apu, options.tr_offsets, options.weighted_weights, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, options.mean_type, options.epps, options.w_sum, options.med_no_norm);
    options.im_vectors.Weighted_BSREM(:,iter+1) = BSREM_iter(options.im_vectors.Weighted_BSREM_apu, options.lam, iter, options.beta_weighted_bsrem, med, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['BSREM weighted mean iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['BSREM weighted mean iteration ' num2str(iter) ' finished'])
    end
    else
    end
end

if options.weighted_mean && options.ROSEM_MAP
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    med = Weighted_mean(options.im_vectors.Weighted_ROSEM_apu, options.tr_offsets, options.weighted_weights, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, options.mean_type, options.epps, options.w_sum, options.med_no_norm);
    options.im_vectors.Weighted_ROSEM(:,iter+1) = BSREM_iter(options.im_vectors.Weighted_ROSEM_apu, options.lam_rosem, iter, options.beta_weighted_rosem, med, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['ROSEM weighted mean iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['ROSEM weighted mean iteration ' num2str(iter) ' finished'])
    end
    else
    end
end

if options.weighted_mean && options.RBI_MAP
    options.im_vectors.Weighted_RBI(:, iter + 1) = options.im_vectors.Weighted_RBI_apu;
end

if options.weighted_mean && any(options.COSEM_MAP)
    options.im_vectors.Weighted_COSEM(:, iter + 1) = options.im_vectors.Weighted_COSEM_apu;
end

if options.TV && options.OSL_OSEM
    options.im_vectors.TV_OSL(:, iter + 1) = options.im_vectors.TV_OSL_apu;
end
if options.TV && options.OSL_MLEM
    options.im_vectors.TV_MLEM(:, iter + 1) = options.im_vectors.TV_MLEM_apu;
end
if options.TV && options.MBSREM
    options.im_vectors.TV_MBSREM(:, iter + 1) = options.im_vectors.TV_MBSREM_apu;
end

if options.TV && options.BSREM
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    grad = TVpriorFinal(options.im_vectors.TV_BSREM_apu, options.TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, options.tr_offsets);
    options.im_vectors.TV_BSREM(:, iter + 1) = BSREM_iter(options.im_vectors.TV_BSREM_apu, options.lam, iter, options.beta_TV_bsrem, grad, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['BSREM TV iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['BSREM TV iteration ' num2str(iter) ' finished'])
    end
    else
    end
end

if options.TV && options.ROSEM_MAP
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    grad = TVpriorFinal(options.im_vectors.TV_ROSEM_apu, options.TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, options.tr_offsets);
    options.im_vectors.TV_ROSEM(:, iter + 1) = BSREM_iter(options.im_vectors.TV_ROSEM_apu, options.lam_rosem, iter, options.beta_TV_rosem, grad, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['ROSEM TV iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['ROSEM TV iteration ' num2str(iter) ' finished'])
    end
    else
    end
end

if options.TV && options.RBI_MAP
    options.im_vectors.TV_RBI(:, iter + 1) = options.im_vectors.TV_RBI_apu;
end

if options.TV && any(options.COSEM_MAP)
    options.im_vectors.TV_COSEM(:, iter + 1) = options.im_vectors.TV_COSEM_apu;
end

if options.AD && options.OSL_OSEM
    options.im_vectors.AD_OSL(:, iter + 1) = options.im_vectors.AD_OSL_apu;
end
if options.AD && options.OSL_MLEM
    options.im_vectors.AD_MLEM(:, iter + 1) = options.im_vectors.AD_MLEM_apu;
end
if options.AD && options.MBSREM
    options.im_vectors.AD_MBSREM(:, iter + 1) = options.im_vectors.AD_MBSREM_apu;
end

if options.AD && options.BSREM
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    med = AD(options.im_vectors.AD_BSREM_apu, options.FluxType, options.Nx, options.Ny, options.Nz, options);
    options.im_vectors.AD_BSREM(:,iter+1) = BSREM_iter(options.im_vectors.AD_BSREM_apu, options.lam, iter, options.beta_ad_bsrem, med, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['BSREM AD iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['BSREM AD iteration ' num2str(iter) ' finished'])
    end
    else
    end
end
if options.AD && options.ROSEM_MAP
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    med = AD(options.im_vectors.AD_ROSEM_apu, options.FluxType, options.Nx, options.Ny, options.Nz, options);
    options.im_vectors.AD_ROSEM(:, iter + 1) = BSREM_iter(options.im_vectors.AD_ROSEM_apu, options.lam_rosem, iter, options.beta_ad_rosem, med, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['ROSEM AD iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['ROSEM AD iteration ' num2str(iter) ' finished'])
    end
    else
    end
end
if options.AD && options.RBI_MAP
    options.im_vectors.AD_RBI(:, iter + 1) = options.im_vectors.AD_RBI_apu;
end
if options.AD && any(options.COSEM_MAP)
    options.im_vectors.AD_COSEM(:, iter + 1) = options.im_vectors.AD_COSEM_apu;
end

if options.APLS && options.OSL_OSEM
    options.im_vectors.APLS_OSL(:, iter + 1) = options.im_vectors.APLS_OSL_apu;
end
if options.APLS && options.OSL_MLEM
    options.im_vectors.APLS_MLEM(:, iter + 1) = options.im_vectors.APLS_MLEM_apu;
end
if options.APLS && options.MBSREM
    options.im_vectors.APLS_MBSREM(:, iter + 1) = options.im_vectors.APLS_MBSREM_apu;
end
if options.APLS && options.BSREM
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    grad = TVpriorFinal(options.im_vectors.APLS_BSREM(:, iter), 0, options.Nx, options.Ny, options.Nz, true, options, 4);
    options.im_vectors.APLS_BSREM(:, iter + 1) = BSREM_iter(options.im_vectors.APLS_BSREM_apu, options.lam, iter, options.beta_APLS_bsrem, grad, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['BSREM APLS iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['BSREM APLS iteration ' num2str(iter) ' finished'])
    end
    else
    end
end
if options.APLS && options.ROSEM_MAP
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    grad = TVpriorFinal(options.im_vectors.APLS_ROSEM(:, iter), 0, options.Nx, options.Ny, options.Nz, true, options, 4);
    options.im_vectors.APLS_ROSEM(:, iter + 1) = BSREM_iter(options.im_vectors.APLS_ROSEM_apu, options.lam_rosem, iter, options.beta_APLS_rosem, grad, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['ROSEM APLS iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['ROSEM APLS iteration ' num2str(iter) ' finished'])
    end
    else
    end
end
if options.APLS && options.RBI_MAP
    options.im_vectors.APLS_RBI(:, iter + 1) = options.im_vectors.APLS_RBI_apu;
end
if options.APLS && any(options.COSEM_MAP)
    options.im_vectors.APLS_COSEM(:, iter + 1) = options.im_vectors.APLS_COSEM_apu;
end

if options.TGV && options.OSL_OSEM
    options.im_vectors.TGV_OSL(:, iter + 1) = options.im_vectors.TGV_OSL_apu;
end
if options.TGV && options.OSL_MLEM
    options.im_vectors.TGV_MLEM(:, iter + 1) = options.im_vectors.TGV_MLEM_apu;
end
if options.TGV && options.MBSREM
    options.im_vectors.TGV_MBSREM(:, iter + 1) = options.im_vectors.TGV_MBSREM_apu;
end
if options.TGV && options.BSREM
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    grad = TGV(options.im_vectors.TGV_BSREM(:, iter),options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
    options.im_vectors.TGV_BSREM(:, iter + 1) = BSREM_iter(options.im_vectors.TGV_BSREM_apu, options.lam, iter, options.beta_TGV_bsrem, grad, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['BSREM TGV iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['BSREM TGV iteration ' num2str(iter) ' finished'])
    end
    else
    end
end
if options.TGV && options.ROSEM_MAP
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    grad = TGV(options.im_vectors.TGV_ROSEM(:, iter),options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
    options.im_vectors.TGV_ROSEM(:, iter + 1) = BSREM_iter(options.im_vectors.TGV_ROSEM_apu, options.lam_rosem, iter, options.beta_TGV_rosem, grad, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['ROSEM TGV iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['ROSEM TGV iteration ' num2str(iter) ' finished'])
    end
    else
    end
end
if options.TGV && options.RBI_MAP
    options.im_vectors.TGV_RBI(:, iter + 1) = options.im_vectors.TGV_RBI_apu;
end
if options.TGV && any(options.COSEM_MAP)
    options.im_vectors.TGV_COSEM(:, iter + 1) = options.im_vectors.TGV_COSEM_apu;
end
if options.NLM && options.OSL_OSEM
    options.im_vectors.NLM_OSL(:, iter + 1) = options.im_vectors.NLM_OSL_apu;
end
if options.NLM && options.MBSREM
    options.im_vectors.NLM_MBSREM(:, iter + 1) = options.im_vectors.NLM_MBSREM_apu;
end
if options.NLM && options.BSREM
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    med = NLM(options.im_vectors.NLM_BSREM_apu, options.options.Ndx, options.options.Ndy, options.options.Ndz, options.Nlx, options.Nly, options.Nlz, options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
    options.im_vectors.NLM_BSREM(:, iter + 1) = BSREM_iter(options.im_vectors.NLM_BSREM_apu, options.lam_rosem, iter, options.beta_NLM_bsrem, med, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['BSREM NLM iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['BSREM NLM iteration ' num2str(iter) ' finished'])
    end
    else
    end
end
if options.NLM && options.ROSEM_MAP
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    med = NLM(options.im_vectors.NLM_ROSEM_apu, options.options.Ndx, options.options.Ndy, options.options.Ndz, options.Nlx, options.Nly, options.Nlz, options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
    options.im_vectors.NLM_ROSEM(:, iter + 1) = BSREM_iter(options.im_vectors.NLM_ROSEM_apu, options.lam_rosem, iter, options.beta_NLM_rosem, med, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['ROSEM NLM iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['ROSEM NLM iteration ' num2str(iter) ' finished'])
    end
    else
    end
end
if options.NLM && options.RBI_MAP
    options.im_vectors.NLM_RBI(:, iter + 1) = options.im_vectors.NLM_RBI_apu;
end
if options.NLM && any(options.COSEM_MAP)
    options.im_vectors.NLM_COSEM(:, iter + 1) = options.im_vectors.NLM_COSEM_apu;
end

if options.OSL_OSEM && options.custom
    options.im_vectors.custom_OSL(:, iter + 1) = options.im_vectors.custom_OSL_apu;
end
if options.OSL_MLEM && options.custom
    options.im_vectors.custom_MLEM(:, iter + 1) = options.im_vectors.custom_MLEM_apu;
end
if options.MBSREM && options.custom
    options.im_vectors.custom_MBSREM(:, iter + 1) = options.im_vectors.custom_MBSREM_apu;
end
if options.BSREM && options.custom
    if options.implementation == 1
        if verbose
            tStart = tic;
        end
        options.im_vectors.custom_BSREM(:, iter + 1) = BSREM_iter(options.im_vectors.custom_BSREM_apu, options.lam_rosem, iter, options.beta_custom_bsrem, options.grad_BSREM, options.epps);
        if verbose
            tElapsed = toc(tStart);
            disp(['BSREM custom prior iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
        else
            disp(['BSREM custom prior iteration ' num2str(iter) ' finished'])
        end
    else
        
    end
end
if options.ROSEM_MAP && options.custom
    if options.implementation == 1
    if verbose
        tStart = tic;
    end
    options.im_vectors.custom_ROSEM(:, iter + 1) = BSREM_iter(options.im_vectors.custom_ROSEM_apu, options.lam_rosem, iter, options.beta_custom_rosem, options.grad_ROSEM, options.epps);
    if verbose
        tElapsed = toc(tStart);
        disp(['ROSEM custom prior iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp(['ROSEM custom prior iteration ' num2str(iter) ' finished'])
    end
    else
    end
end
if options.RBI_MAP && options.custom
    options.im_vectors.custom_RBI(:, iter + 1) = options.im_vectors.custom_RBI_apu;
end
if any(options.COSEM_MAP) && options.custom
    options.im_vectors.custom_COSEM(:, iter + 1) = options.im_vectors.custom_COSEM_apu;
end
disp(['Iteration ' num2str(iter) ' finished'])