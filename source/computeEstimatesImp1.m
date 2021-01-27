function [im_vectors,C_co,C_aco,C_osl] = computeEstimatesImp1(im_vectors, options, A, uu, Summ, SinD, is_transposed, gaussK, iter, osa_iter, C_co, C_aco,C_osl,...
    randoms_correction, N, Ndx, Ndy, Ndz, D)
%COMPUTEESTIMATESIMP1 Computes the subset estimates for implementation 1
%   Utility function
% Compute OSEM
if options.osem || options.ecosem || options.attenuation_phase
    if options.verbose
        tStart = tic;
    end
    im_vectors.OSEM_apu = OSEM_im(im_vectors.OSEM_apu, A, options.epps, uu, Summ, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['OSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute MRAMLA
if options.mramla
    if options.verbose
        tStart = tic;
    end
    im_vectors.MRAMLA_apu = MBSREM(im_vectors.MRAMLA_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, ...
        SinD, randoms_correction, is_transposed, [], [], options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['MRAMLA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute RAMLA
if options.ramla
    if options.verbose
        tStart = tic;
    end
    im_vectors.RAMLA_apu = BSREM_subiter(im_vectors.RAMLA_apu, options.lam, options.epps, iter, A, uu, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    %                             if any(im_vectors.RAMLA_apu < 0)
    %                                 warning('Negative values in RAMLA, lower lambda value!')
    %                             end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['RAMLA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute ROSEM
if options.rosem
    if options.verbose
        tStart = tic;
    end
    im_vectors.ROSEM_apu = ROSEM_subiter(im_vectors.ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['ROSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute RBI
if options.rbi
    if options.verbose
        tStart = tic;
    end
    im_vectors.RBI_apu = RBI_subiter(im_vectors.RBI_apu, A, uu, options.epps, Summ, D, SinD, is_transposed, [], [], options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['RBI sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute DRAMA
if options.drama
    if options.verbose
        tStart = tic;
    end
    im_vectors.DRAMA_apu = DRAMA_subiter(im_vectors.DRAMA_apu, options.lam_drama, options.epps, iter, Summ, osa_iter, A, uu, SinD, is_transposed, options, ...
        options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['DRAMA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute COSEM
if options.cosem || options.ecosem
    if options.verbose
        tStart = tic;
    end
    [im_vectors.COSEM_apu, C_co] = COSEM_im(im_vectors.COSEM_apu, A, options.epps, uu, C_co, D, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['COSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute ECOSEM
if options.ecosem
    if options.verbose
        tStart = tic;
    end
    im_vectors.ECOSEM_apu = ECOSEM_im(im_vectors.ECOSEM_apu, options.epps, D, im_vectors.COSEM_apu, im_vectors.OSEM_apu);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['ECOSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute ACOSEM
if options.acosem
    if options.verbose
        tStart = tic;
    end
    [im_vectors.ACOSEM_apu, C_aco] = ACOSEM_im(im_vectors.ACOSEM_apu, A, options.epps, uu, C_aco, D, options.h, osa_iter, SinD, is_transposed, options, ...
        options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['ACOSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute OSL with MRP
if options.MRP && options.OSL_OSEM
    if options.verbose
        tStart = tic;
    end
    med = MRP(im_vectors.MRP_OSL_apu, options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
    im_vectors.MRP_OSL_apu = OSEM_im(im_vectors.MRP_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_mrp_osem, med, options.epps), SinD, is_transposed, ...
        options, options.Nx, options.Ny, options.Nz, gaussK);
    %                             im_vectors.MRP_OSL_apu = OSL_OSEM(im_vectors.MRP_OSL_apu, Summ, options.beta_mrp_osem, med, options.epps, A, uu, SinD, is_transposed);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['OSEM-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute MBSREM with MRP
if options.MRP && options.MBSREM
    if options.verbose
        tStart = tic;
    end
    med = MRP(im_vectors.MRP_MBSREM_apu, options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
    im_vectors.MRP_MBSREM_apu = MBSREM(im_vectors.MRP_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
        iter, SinD, randoms_correction, is_transposed, options.beta_mrp_mbsrem, med, options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['MBSREM MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute BSREM with MRP
if options.MRP && options.BSREM
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.ramla
        im_vectors.MRP_BSREM_apu = im_vectors.RAMLA_apu;
    else
        im_vectors.MRP_BSREM_apu = BSREM_subiter(im_vectors.MRP_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, is_transposed, options, ...
            options.Nx, options.Ny, options.Nz, gaussK);
    end
    %                             if any(im_vectors.MRP_BSREM_apu < 0)
    %                                 warning('Negative values in BSREM, lower lambda value!')
    %                             end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['BSREM MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute ROSEM-MAP with MRP
if options.MRP && options.ROSEM_MAP
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.rosem
        im_vectors.MRP_ROSEM_apu = im_vectors.ROSEM_apu;
    else
        im_vectors.MRP_ROSEM_apu = ROSEM_subiter(im_vectors.MRP_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, is_transposed, ...
            options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['ROSEM-MAP MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute RBI-OSL with MRP
if options.MRP && options.RBI_OSL
    if options.verbose
        tStart = tic;
    end
    med = MRP(im_vectors.MRP_RBI_apu, options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
    im_vectors.MRP_RBI_apu = RBI_subiter(im_vectors.MRP_RBI_apu, A, uu, options.epps, Summ, D, SinD, is_transposed, options.beta_mrp_rbi, med, options, ...
        options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['RBI-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute COSEM-OSL with MRP
if options.MRP && any(options.COSEM_OSL)
    if options.verbose
        tStart = tic;
    end
    med = MRP(im_vectors.MRP_COSEM_apu, options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
    if options.COSEM_OSL == 1
        [im_vectors.MRP_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.MRP_COSEM_apu, D, options.beta_mrp_cosem, med, A, uu, ...
            options.epps, C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    else
        [im_vectors.MRP_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.MRP_COSEM_apu, D, options.beta_mrp_cosem, med, A, uu, ...
            options.epps, C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['COSEM-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute OSL with Quadratic prior
if options.quad && options.OSL_OSEM
    if options.verbose
        tStart = tic;
    end
    med = Quadratic_prior(im_vectors.Quad_OSL_apu, options.weights_quad, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, options);
    im_vectors.Quad_OSL_apu = OSEM_im(im_vectors.Quad_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_quad_osem, med, options.epps), SinD, is_transposed, ...
        options, options.Nx, options.Ny, options.Nz, gaussK);
    %                             im_vectors.Quad_OSL_apu = OSL_OSEM(im_vectors.Quad_OSL_apu, Summ, options.beta_quad_osem, med, options.epps, A, uu, SinD, is_transposed);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['OSEM-OSL Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute MBSREM with Quadratic prior
if options.quad && options.MBSREM
    if options.verbose
        tStart = tic;
    end
    med = Quadratic_prior(im_vectors.Quad_MBSREM_apu, options.weights_quad, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz);
    im_vectors.Quad_MBSREM_apu = MBSREM(im_vectors.Quad_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
        iter, SinD, randoms_correction, is_transposed, options.beta_quad_mbsrem, med, options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['MBSREM Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute BSREM with Quadratic prior
if options.quad && options.BSREM
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.ramla
        im_vectors.Quad_BSREM_apu = im_vectors.RAMLA_apu;
    else
        im_vectors.Quad_BSREM_apu = BSREM_subiter(im_vectors.Quad_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, is_transposed, options, ...
            options.Nx, options.Ny, options.Nz, gaussK);
    end
    %                             if any(im_vectors.Quad_BSREM_apu < 0)
    %                                 warning('Negative values in BSREM, it is recommended to lower lambda value')
    %                             end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['BSREM Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute ROSEM-MAP with Quadratic prior
if options.quad && options.ROSEM_MAP
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.rosem
        im_vectors.Quad_ROSEM_apu = im_vectors.ROSEM_apu;
    else
        im_vectors.Quad_ROSEM_apu = ROSEM_subiter(im_vectors.Quad_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, is_transposed, ...
            options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['ROSEM-MAP Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute RBI-OSL with Quadratic prior
if options.quad && options.RBI_OSL
    if options.verbose
        tStart = tic;
    end
    med = Quadratic_prior(im_vectors.Quad_RBI_apu, options.weights_quad, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz);
    im_vectors.Quad_RBI_apu = RBI_subiter(im_vectors.Quad_RBI_apu, A, uu, options.epps, Summ, D, SinD, is_transposed, options.beta_quad_rbi, med, ...
        options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['RBI-OSL Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute COSEM-OSL with Quadratic prior
if options.quad && any(options.COSEM_OSL)
    if options.verbose
        tStart = tic;
    end
    med = Quadratic_prior(im_vectors.Quad_COSEM_apu, options.weights_quad, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz);
    if options.COSEM_OSL == 1
        [im_vectors.Quad_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.Quad_COSEM_apu, D, options.beta_quad_cosem, med, A, uu, ...
            options.epps, C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    else
        [im_vectors.Quad_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.Quad_COSEM_apu, D, options.beta_quad_cosem, med, A, uu, ...
            options.epps, C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['COSEM-OSL Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute OSL with Huber prior
if options.Huber && options.OSL_OSEM
    if options.verbose
        tStart = tic;
    end
    med = Huber_prior(im_vectors.Huber_OSL_apu, options.weights_huber, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, options.huber_delta);
    im_vectors.Huber_OSL_apu = OSEM_im(im_vectors.Huber_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_huber_osem, med, options.epps), SinD, is_transposed, ...
        options, options.Nx, options.Ny, options.Nz, gaussK);
    %                             im_vectors.Huber_OSL_apu = OSL_OSEM(im_vectors.Huber_OSL_apu, Summ, options.beta_huber_osem, med, options.epps, A, uu, SinD, is_transposed);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['OSEM-OSL Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute MBSREM with Huber prior
if options.Huber && options.MBSREM
    if options.verbose
        tStart = tic;
    end
    med = Huber_prior(im_vectors.Huber_MBSREM_apu, options.weights_huber, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, options.huber_delta);
    im_vectors.Huber_MBSREM_apu = MBSREM(im_vectors.Huber_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
        iter, SinD, randoms_correction, is_transposed, options.beta_huber_mbsrem, med, options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['MBSREM Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute BSREM with Huber prior
if options.Huber && options.BSREM
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.ramla
        im_vectors.Huber_BSREM_apu = im_vectors.RAMLA_apu;
    else
        im_vectors.Huber_BSREM_apu = BSREM_subiter(im_vectors.Huber_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, is_transposed, options, ...
            options.Nx, options.Ny, options.Nz, gaussK);
    end
    %                             if any(im_vectors.Huber_BSREM_apu < 0)
    %                                 warning('Negative values in BSREM, it is recommended to lower lambda value')
    %                             end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['BSREM Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute ROSEM-MAP with Huber prior
if options.Huber && options.ROSEM_MAP
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.rosem
        im_vectors.Huber_ROSEM_apu = im_vectors.ROSEM_apu;
    else
        im_vectors.Huber_ROSEM_apu = ROSEM_subiter(im_vectors.Huber_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, is_transposed, ...
            options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['ROSEM-MAP Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute RBI-OSL with Huber prior
if options.Huber && options.RBI_OSL
    if options.verbose
        tStart = tic;
    end
    med = Huber_prior(im_vectors.Huber_RBI_apu, options.weights_huber, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, options.huber_delta);
    im_vectors.Huber_RBI_apu = RBI_subiter(im_vectors.Huber_RBI_apu, A, uu, options.epps, Summ, D, SinD, is_transposed, options.beta_huber_rbi, med, ...
        options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['RBI-OSL Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute COSEM-OSL with Huber prior
if options.Huber && any(options.COSEM_OSL)
    if options.verbose
        tStart = tic;
    end
    med = Huber_prior(im_vectors.Huber_COSEM_apu, options.weights_huber, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, options.huber_delta);
    if options.COSEM_OSL == 1
        [im_vectors.Huber_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.Huber_COSEM_apu, D, options.beta_huber_cosem, med, A, uu, ...
            options.epps, C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    else
        [im_vectors.Huber_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.Huber_COSEM_apu, D, options.beta_huber_cosem, med, A, uu, ...
            options.epps, C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['COSEM-OSL Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute OSL with L-filter prior
if options.L && options.OSL_OSEM
    if options.verbose
        tStart = tic;
    end
    med = L_filter(im_vectors.L_OSL_apu, options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, options.epps, options.med_no_norm);
    im_vectors.L_OSL_apu = OSEM_im(im_vectors.L_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_L_osem, med, options.epps), SinD, is_transposed, options, ...
        options.Nx, options.Ny, options.Nz, gaussK);
    %                             im_vectors.L_OSL_apu = OSL_OSEM(im_vectors.L_OSL_apu, Summ, options.beta_L_osem, med, options.epps, A, uu, SinD, is_transposed);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['OSEM-OSL L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute MBSREM with L-filter prior
if options.L && options.MBSREM
    if options.verbose
        tStart = tic;
    end
    med = L_filter(im_vectors.L_MBSREM_apu, options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, options.epps, options.med_no_norm);
    im_vectors.L_MBSREM_apu = MBSREM(im_vectors.L_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
        iter, SinD, randoms_correction, is_transposed, options.beta_L_mbsrem, med, options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['MBSREM L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute BSREM with L-filter prior
if options.L && options.BSREM
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.ramla
        im_vectors.L_BSREM_apu = im_vectors.RAMLA_apu;
    else
        im_vectors.L_BSREM_apu = BSREM_subiter(im_vectors.L_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    %                             if any(im_vectors.L_BSREM_apu < 0)
    %                                 warning('Negative values in BSREM, it is recommended to lower lambda value')
    %                             end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['BSREM L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute ROSEM-MAP with L-filter prior
if options.L && options.ROSEM_MAP
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.rosem
        im_vectors.L_ROSEM_apu = im_vectors.ROSEM_apu;
    else
        im_vectors.L_ROSEM_apu = ROSEM_subiter(im_vectors.L_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, is_transposed, options, ...
            options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['ROSEM-MAP L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute RBI-OSL with L-filter prior
if options.L && options.RBI_OSL
    if options.verbose
        tStart = tic;
    end
    med = L_filter(im_vectors.L_RBI_apu, options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, options.epps, options.med_no_norm);
    im_vectors.L_RBI_apu = RBI_subiter(im_vectors.L_RBI_apu, A, uu, options.epps, Summ, D, SinD, is_transposed, options.beta_L_rbi, med, options, ...
        options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['RBI-OSL L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute COSEM-OSL with L-filter prior
if options.L && any(options.COSEM_OSL)
    if options.verbose
        tStart = tic;
    end
    med = L_filter(im_vectors.L_COSEM_apu, options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, options.epps, options.med_no_norm);
    if options.COSEM_OSL == 1
        [im_vectors.L_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.L_COSEM_apu, D, options.beta_L_cosem, med, A, uu, options.epps, C_osl, ...
            options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    else
        [im_vectors.L_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.L_COSEM_apu, D, options.beta_L_cosem, med, A, uu, options.epps, C_osl, 0, ...
            options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['COSEM-OSL L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute OSL with FMH prior
if options.FMH && options.OSL_OSEM
    if options.verbose
        tStart = tic;
    end
    med = FMH(im_vectors.FMH_OSL_apu, options.tr_offsets, options.fmh_weights, options.Nx, options.Ny, options.Nz, N, Ndx, Ndy, Ndz, options.epps, ...
        options.med_no_norm);
    im_vectors.FMH_OSL_apu = OSEM_im(im_vectors.FMH_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_fmh_osem, med, options.epps), SinD, is_transposed, options, ...
        options.Nx, options.Ny, options.Nz, gaussK);
    %                             im_vectors.FMH_OSL_apu = OSL_OSEM(im_vectors.FMH_OSL_apu, Summ, options.beta_fmh_osem, med, options.epps, A, uu, SinD, is_transposed);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['OSEM-OSL FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.FMH && options.MBSREM
    if options.verbose
        tStart = tic;
    end
    med = FMH(im_vectors.FMH_MBSREM_apu, options.tr_offsets, options.fmh_weights, options.Nx, options.Ny, options.Nz, N, Ndx, Ndy, Ndz, options.epps, ...
        options.med_no_norm);
    im_vectors.FMH_MBSREM_apu = MBSREM(im_vectors.FMH_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
        iter, SinD, randoms_correction, is_transposed, options.beta_fmh_mbsrem, med, options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['MBSREM FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.FMH && options.BSREM
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.ramla
        im_vectors.FMH_BSREM_apu = im_vectors.RAMLA_apu;
    else
        im_vectors.FMH_BSREM_apu = BSREM_subiter(im_vectors.FMH_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, is_transposed, options, ...
            options.Nx, options.Ny, options.Nz, gaussK);
    end
    %                             if any(im_vectors.FMH_BSREM_apu < 0)
    %                                 warning('Negative values in BSREM, it is recommended to lower lambda value')
    %                             end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['BSREM FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.FMH && options.ROSEM_MAP
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.rosem
        im_vectors.FMH_ROSEM_apu = im_vectors.ROSEM_apu;
    else
        im_vectors.FMH_ROSEM_apu = ROSEM_subiter(im_vectors.FMH_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, is_transposed, ...
            options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['ROSEM-MAP FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.FMH && options.RBI_OSL
    if options.verbose
        tStart = tic;
    end
    med = FMH(im_vectors.FMH_RBI_apu, options.tr_offsets, options.fmh_weights, options.Nx, options.Ny, options.Nz, N, Ndx, Ndy, Ndz, options.epps, ...
        options.med_no_norm);
    im_vectors.FMH_RBI_apu = RBI_subiter(im_vectors.FMH_RBI_apu, A, uu, options.epps, Summ, D, SinD, is_transposed, options.beta_fmh_rbi, med, options, ...
        options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['RBI-OSL FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.FMH && any(options.COSEM_OSL)
    if options.verbose
        tStart = tic;
    end
    med = FMH(im_vectors.FMH_COSEM_apu, options.tr_offsets, options.fmh_weights, options.Nx, options.Ny, options.Nz, N, Ndx, Ndy, Ndz, options.epps, ...
        options.med_no_norm);
    if options.COSEM_OSL == 1
        [im_vectors.FMH_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.FMH_COSEM_apu, D, options.beta_fmh_cosem, med, A, uu, ...
            options.epps, C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    else
        [im_vectors.FMH_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.FMH_COSEM_apu, D, options.beta_fmh_cosem, med, A, uu, ...
            options.epps, C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['COSEM-OSL FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute OSL with weighted mean prior
if options.weighted_mean && options.OSL_OSEM
    if options.verbose
        tStart = tic;
    end
    med = Weighted_mean(im_vectors.Weighted_OSL_apu, options.weighted_weights, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, ...
        options.mean_type, options.epps, options.med_no_norm);
    im_vectors.Weighted_OSL_apu = OSEM_im(im_vectors.Weighted_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_weighted_osem, med, options.epps), SinD, ...
        is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    %                             im_vectors.Weighted_OSL_apu = OSL_OSEM(im_vectors.Weighted_OSL_apu, Summ, options.beta_weighted_osem, med, options.epps, A, uu, SinD, ...
    %                                 is_transposed);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['OSEM-OSL weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.weighted_mean && options.MBSREM
    if options.verbose
        tStart = tic;
    end
    med = Weighted_mean(im_vectors.Weighted_MBSREM_apu, options.weighted_weights, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, ...
        options.mean_type, options.epps, options.med_no_norm);
    im_vectors.Weighted_MBSREM_apu = MBSREM(im_vectors.Weighted_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, ...
        options.lam_mbsrem, iter, SinD, randoms_correction, is_transposed, options.beta_weighted_mbsrem, med, options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['MBSREM weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.weighted_mean && options.BSREM
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.ramla
        im_vectors.Weighted_BSREM_apu = im_vectors.RAMLA_apu;
    else
        im_vectors.Weighted_BSREM_apu = BSREM_subiter(im_vectors.Weighted_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, is_transposed, ...
            options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    %                             if any(im_vectors.Weighted_BSREM_apu < 0)
    %                                 warning('Negative values in BSREM, it is recommended to lower lambda value')
    %                             end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['BSREM weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.weighted_mean && options.ROSEM_MAP
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.rosem
        im_vectors.Weighted_ROSEM_apu = im_vectors.ROSEM_apu;
    else
        im_vectors.Weighted_ROSEM_apu = ROSEM_subiter(im_vectors.Weighted_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, ...
            is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['ROSEM-MAP weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.weighted_mean && options.RBI_OSL
    if options.verbose
        tStart = tic;
    end
    med = Weighted_mean(im_vectors.Weighted_RBI_apu, options.weighted_weights, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, ...
        options.mean_type, options.epps, options.med_no_norm);
    im_vectors.Weighted_RBI_apu = RBI_subiter(im_vectors.Weighted_RBI_apu, A, uu, options.epps, Summ, ...
        D, SinD, is_transposed, options.beta_weighted_rbi, med, options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['RBI-OSL weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.weighted_mean && any(options.COSEM_OSL)
    if options.verbose
        tStart = tic;
    end
    med = Weighted_mean(im_vectors.Weighted_COSEM_apu, options.weighted_weights, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, ...
        options.mean_type, options.epps, options.med_no_norm);
    if options.COSEM_OSL == 1
        [im_vectors.Weighted_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.Weighted_COSEM_apu, D, options.beta_weighted_cosem, ...
            med, A, uu, options.epps, C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    else
        [im_vectors.Weighted_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.Weighted_COSEM_apu, D, options.beta_weighted_cosem, ...
            med, A, uu, options.epps, C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['COSEM-OSL weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute OSL with TV prior
if options.TV && options.OSL_OSEM
    if options.verbose
        tStart = tic;
    end
    grad = TVpriorFinal(im_vectors.TV_OSL_apu, options.TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, ...
        options.tr_offsets);
    im_vectors.TV_OSL_apu = OSEM_im(im_vectors.TV_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_TV_osem, grad, options.epps), SinD, is_transposed, ...
        options, options.Nx, options.Ny, options.Nz, gaussK);
    %                             im_vectors.TV_OSL_apu = OSL_OSEM(im_vectors.TV_OSL_apu, Summ, options.beta_TV_osem, grad, options.epps, A, uu, SinD, is_transposed);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['OSEM-OSL TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.TV && options.MBSREM
    if options.verbose
        tStart = tic;
    end
    grad = TVpriorFinal(im_vectors.TV_MBSREM_apu, options.TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, ...
        options.tr_offsets);
    im_vectors.TV_MBSREM_apu = MBSREM(im_vectors.TV_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
        iter, SinD, randoms_correction, is_transposed, options.beta_TV_mbsrem, grad, options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['MBSREM TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.TV && options.BSREM
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.ramla
        im_vectors.TV_BSREM_apu = im_vectors.RAMLA_apu;
    else
        im_vectors.TV_BSREM_apu = BSREM_subiter(im_vectors.TV_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, is_transposed, options, ...
            options.Nx, options.Ny, options.Nz, gaussK);
    end
    %                             if any(im_vectors.TV_BSREM_apu < 0)
    %                                 warning('Negative values in BSREM, it is recommended to lower lambda value')
    %                             end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['BSREM TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.TV && options.ROSEM_MAP
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.rosem
        im_vectors.TV_ROSEM_apu = im_vectors.ROSEM_apu;
    else
        im_vectors.TV_ROSEM_apu = ROSEM_subiter(im_vectors.TV_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, is_transposed, ...
            options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['ROSEM-MAP TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.TV && options.RBI_OSL
    if options.verbose
        tStart = tic;
    end
    grad = TVpriorFinal(im_vectors.TV_RBI_apu, options.TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, ...
        options.tr_offsets);
    im_vectors.TV_RBI_apu = RBI_subiter(im_vectors.TV_RBI_apu, A, uu, options.epps, Summ, D, SinD, is_transposed, options.beta_TV_rbi, grad, options, ...
        options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['RBI-OSL TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.TV && any(options.COSEM_OSL)
    if options.verbose
        tStart = tic;
    end
    grad = TVpriorFinal(im_vectors.TV_COSEM_apu, options.TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, ...
        options.tr_offsets);
    if options.COSEM_OSL == 1
        [im_vectors.TV_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.TV_COSEM_apu, D, options.beta_TV_cosem, grad, A, uu, ...
            options.epps, C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    else
        [im_vectors.TV_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.TV_COSEM_apu, D, options.beta_TV_cosem, grad, A, uu, ...
            options.epps, C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['COSEM-OSL TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute OSL with MRP-AD prior
if options.AD && options.OSL_OSEM
    if options.verbose
        tStart = tic;
    end
    if osa_iter > 1
        med = AD(im_vectors.AD_OSL_apu, options.FluxType, options.Nx, options.Ny, options.Nz, options);
        im_vectors.AD_OSL_apu = OSEM_im(im_vectors.AD_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_ad_osem, med, options.epps), SinD, is_transposed, ...
            options, options.Nx, options.Ny, options.Nz, gaussK);
        %                                 im_vectors.AD_OSL_apu = OSL_OSEM(im_vectors.AD_OSL_apu, Summ, options.beta_ad_osem, med, options.epps, A, uu, SinD, is_transposed);
    else
        im_vectors.AD_OSL_apu = OSEM_im(im_vectors.AD_OSL_apu, A, options.epps, uu, Summ, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['OSEM-OSL AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.AD && options.MBSREM
    if options.verbose
        tStart = tic;
    end
    if osa_iter > 1
        med = AD(im_vectors.AD_MBSREM_apu, options.FluxType, options.Nx, options.Ny, options.Nz, options);
        im_vectors.AD_MBSREM_apu = MBSREM(im_vectors.AD_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
            iter, SinD, randoms_correction, is_transposed, options.beta_ad_mbsrem, med, options, options.Nx, options.Ny, options.Nz, gaussK);
    else
        im_vectors.AD_MBSREM_apu = MBSREM(im_vectors.AD_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, ...
            SinD, randoms_correction, is_transposed, [], [], options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['MBSREM AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.AD && options.BSREM
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.ramla
        im_vectors.AD_BSREM_apu = im_vectors.RAMLA_apu;
    else
        im_vectors.AD_BSREM_apu = BSREM_subiter(im_vectors.AD_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    %                             if any(im_vectors.AD_BSREM_apu < 0)
    %                                 warning('Negative values in BSREM, it is recommended to lower lambda value')
    %                             end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['BSREM AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.AD && options.ROSEM_MAP
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.rosem
        im_vectors.AD_ROSEM_apu = im_vectors.ROSEM_apu;
    else
        im_vectors.AD_ROSEM_apu = ROSEM_subiter(im_vectors.AD_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, is_transposed, options, ...
            options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['ROSEM-MAP AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.AD && options.RBI_OSL
    if options.verbose
        tStart = tic;
    end
    med = AD(im_vectors.AD_RBI_apu, options.FluxType, options.Nx, options.Ny, options.Nz, options);
    im_vectors.AD_RBI_apu = RBI_subiter(im_vectors.AD_RBI_apu, A, uu, options.epps, Summ, D, SinD, is_transposed, options.beta_ad_rbi, med, options, ...
        options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['RBI-OSL AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.AD && any(options.COSEM_OSL)
    if options.verbose
        tStart = tic;
    end
    med = AD(im_vectors.AD_COSEM_apu, options.FluxType, options.Nx, options.Ny, options.Nz, options);
    if options.COSEM_OSL == 1
        [im_vectors.AD_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.AD_COSEM_apu, D, options.beta_ad_cosem, med, A, uu, ...
            options.epps, C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    else
        [im_vectors.AD_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.AD_COSEM_apu, D, options.beta_ad_cosem, med, A, uu, ...
            options.epps, C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['COSEM-OSL AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute OSL with APLS prior
if options.APLS && options.OSL_OSEM
    if options.verbose
        tStart = tic;
    end
    grad = TVpriorFinal(im_vectors.APLS_OSL_apu, [], options.Nx, options.Ny, options.Nz, true, options, 5);
    im_vectors.APLS_OSL_apu = OSEM_im(im_vectors.APLS_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_APLS_osem, grad, options.epps), SinD, is_transposed, options, ...
        options.Nx, options.Ny, options.Nz, gaussK);
    %                             im_vectors.APLS_OSL_apu = OSL_OSEM(im_vectors.APLS_OSL_apu, Summ, options.beta_APLS_osem, grad, options.epps, A, uu, SinD, is_transposed);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['OSEM-OSL APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.APLS && options.MBSREM
    if options.verbose
        tStart = tic;
    end
    grad = TVpriorFinal(im_vectors.APLS_MBSREM_apu, [], options.Nx, options.Ny, options.Nz, true, options, 5);
    im_vectors.APLS_MBSREM_apu = MBSREM(im_vectors.APLS_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
        iter, SinD, randoms_correction, is_transposed, options.beta_APLS_mbsrem, grad, options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['MBSREM APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.APLS && options.BSREM
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.ramla
        im_vectors.APLS_BSREM_apu = im_vectors.RAMLA_apu;
    else
        im_vectors.APLS_BSREM_apu = BSREM_subiter(im_vectors.APLS_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, is_transposed, options, ...
            options.Nx, options.Ny, options.Nz, gaussK);
    end
    %                             if any(im_vectors.APLS_BSREM_apu < 0)
    %                                 warning('Negative values in BSREM, it is recommended to lower lambda value')
    %                             end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['BSREM APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.APLS && options.ROSEM_MAP
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.rosem
        im_vectors.APLS_ROSEM_apu = im_vectors.ROSEM_apu;
    else
        im_vectors.APLS_ROSEM_apu = ROSEM_subiter(im_vectors.APLS_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, is_transposed, ...
            options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['ROSEM-MAP APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.APLS && options.RBI_OSL
    if options.verbose
        tStart = tic;
    end
    grad = TVpriorFinal(im_vectors.APLS_RBI_apu, [], options.Nx, options.Ny, options.Nz, true, options, 5);
    im_vectors.APLS_RBI_apu = RBI_subiter(im_vectors.APLS_RBI_apu, A, uu, options.epps, Summ, D, SinD, ...
        is_transposed, options.beta_APLS_rbi, grad, options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['RBI-OSL APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.APLS && any(options.COSEM_OSL)
    if options.verbose
        tStart = tic;
    end
    grad = TVpriorFinal(im_vectors.APLS_COSEM_apu, [], options.Nx, options.Ny, options.Nz, true, options, 5);
    if options.COSEM_OSL == 1
        [im_vectors.APLS_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.APLS_COSEM_apu, D, options.beta_APLS_cosem, grad, A, uu, ...
            options.epps, C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    else
        [im_vectors.APLS_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.APLS_COSEM_apu, D, options.beta_APLS_cosem, grad, A, uu, ...
            options.epps, C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['COSEM-OSL APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute OSL with TGV prior
if options.TGV && options.OSL_OSEM
    if options.verbose
        tStart = tic;
    end
    grad = TGV(im_vectors.TGV_OSL_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
    im_vectors.TGV_OSL_apu = OSEM_im(im_vectors.TGV_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_TGV_osem, grad, options.epps), SinD, is_transposed, options, ...
        options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['OSEM-OSL TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.TGV && options.MBSREM
    if options.verbose
        tStart = tic;
    end
    grad = TGV(im_vectors.TGV_MBSREM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
    im_vectors.TGV_MBSREM_apu = MBSREM(im_vectors.TGV_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
        iter, SinD, randoms_correction, is_transposed, options.beta_TGV_mbsrem, grad, options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['MBSREM TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.TGV && options.BSREM
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.ramla
        im_vectors.TGV_BSREM_apu = im_vectors.RAMLA_apu;
    else
        im_vectors.TGV_BSREM_apu = BSREM_subiter(im_vectors.TGV_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, is_transposed, options, ...
            options.Nx, options.Ny, options.Nz, gaussK);
    end
    %                             if any(im_vectors.TGV_BSREM_apu < 0)
    %                                 warning('Negative values in BSREM, it is recommended to lower lambda value')
    %                             end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['BSREM TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.TGV && options.ROSEM_MAP
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.rosem
        im_vectors.TGV_ROSEM_apu = im_vectors.ROSEM_apu;
    else
        im_vectors.TGV_ROSEM_apu = ROSEM_subiter(im_vectors.TGV_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, is_transposed, options, ...
            options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['ROSEM-MAP TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.TGV && options.RBI_OSL
    if options.verbose
        tStart = tic;
    end
    grad = TGV(im_vectors.TGV_RBI_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
    im_vectors.TGV_RBI_apu = RBI_subiter(im_vectors.TGV_RBI_apu, A, uu, options.epps, Summ, D, SinD, is_transposed, options.beta_TGV_rbi, grad, options, ...
        options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['RBI-OSL TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.TGV && any(options.COSEM_OSL)
    if options.verbose
        tStart = tic;
    end
    grad = TGV(im_vectors.TGV_COSEM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
    if options.COSEM_OSL == 1
        [im_vectors.TGV_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.TGV_COSEM_apu, D, options.beta_TGV_cosem, grad, A, uu, ...
            options.epps, C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    else
        [im_vectors.TGV_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.TGV_COSEM_apu, D, options.beta_TGV_cosem, grad, A, uu, ...
            options.epps, C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['COSEM-OSL TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute OSL with NLM prior
if options.NLM && options.OSL_OSEM
    if options.verbose
        tStart = tic;
    end
    med = NLM(im_vectors.NLM_OSL_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
        options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
    im_vectors.NLM_OSL_apu = OSEM_im(im_vectors.NLM_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_NLM_osem, med, options.epps), SinD, is_transposed, options, ...
        options.Nx, options.Ny, options.Nz, gaussK);
    %                             im_vectors.NLM_OSL_apu = OSL_OSEM(im_vectors.NLM_OSL_apu, Summ, options.beta_NLM_osem, med, options.epps, A, uu, SinD, is_transposed);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['OSEM-OSL NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.NLM && options.MBSREM
    if options.verbose
        tStart = tic;
    end
    med = NLM(im_vectors.NLM_MBSREM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
        options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
    im_vectors.NLM_MBSREM_apu = MBSREM(im_vectors.NLM_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
        iter, SinD, randoms_correction, is_transposed, options.beta_NLM_mbsrem, med, options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['MBSREM NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.NLM && options.BSREM
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.ramla
        im_vectors.NLM_BSREM_apu = im_vectors.RAMLA_apu;
    else
        im_vectors.NLM_BSREM_apu = BSREM_subiter(im_vectors.NLM_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, is_transposed, options, ...
            options.Nx, options.Ny, options.Nz, gaussK);
    end
    %                             if any(im_vectors.NLM_BSREM_apu < 0)
    %                                 warning('Negative values in BSREM, it is recommended to lower lambda value')
    %                             end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['BSREM NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.NLM && options.ROSEM_MAP
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.rosem
        im_vectors.NLM_ROSEM_apu = im_vectors.ROSEM_apu;
    else
        im_vectors.NLM_ROSEM_apu = ROSEM_subiter(im_vectors.NLM_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, is_transposed, ...
            options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['ROSEM-MAP NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.NLM && options.RBI_OSL
    if options.verbose
        tStart = tic;
    end
    med = NLM(im_vectors.NLM_RBI_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
        options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
    im_vectors.NLM_RBI_apu = RBI_subiter(im_vectors.NLM_RBI_apu, A, uu, options.epps, Summ, D, SinD, is_transposed, options.beta_NLM_rbi, med, options, ...
        options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['RBI-OSL NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.NLM && any(options.COSEM_OSL)
    if options.verbose
        tStart = tic;
    end
    med = NLM(im_vectors.NLM_COSEM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
        options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
    if options.COSEM_OSL == 1
        [im_vectors.NLM_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.NLM_COSEM_apu, D, options.beta_NLM_cosem, med, A, uu, ...
            options.epps, C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    else
        [im_vectors.NLM_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.NLM_COSEM_apu, D, options.beta_NLM_cosem, med, A, uu, ...
            options.epps, C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['COSEM-OSL NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end

if options.OSL_OSEM && options.custom
    if options.verbose
        tStart = tic;
    end
    im_vectors.custom_OSL_apu = OSEM_im(im_vectors.custom_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_custom_osem, options.grad_OSEM, options.epps), SinD, ...
        is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['OSL custom prior sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.MBSREM && options.custom
    if options.verbose
        tStart = tic;
    end
    im_vectors.custom_MBSREM_apu = MBSREM(im_vectors.custom_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, ...
        options.lam_mbsrem, iter, SinD, randoms_correction, is_transposed, options.beta_custom_mbsrem, options.grad_MBSREM, options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['MBSREM custom prior sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.BSREM && options.custom
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.ramla
        im_vectors.custom_BSREM_apu = im_vectors.RAMLA_apu;
    else
        im_vectors.custom_BSREM_apu = BSREM_subiter(im_vectors.custom_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, is_transposed, ...
            options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    %         if any(im_vectors.custom_BSREM(:,iter) < 0)
    %             error('Negative values in BSREM, lower options.lambda value')
    %         end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['BSREM custom prior sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.ROSEM_MAP && options.custom
    if options.verbose
        tStart = tic;
    end
    if iter == 1 && options.rosem
        im_vectors.custom_ROSEM_apu = im_vectors.ROSEM_apu;
    else
        im_vectors.custom_ROSEM_apu = ROSEM_subiter(im_vectors.custom_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, ...
            is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['ROSEM-MAP custom prior sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if options.RBI_OSL && options.custom
    if options.verbose
        tStart = tic;
    end
    im_vectors.custom_RBI_apu = RBI_subiter(im_vectors.custom_RBI_apu, A, uu, options.epps, Summ, D, SinD, is_transposed, ...
        options.beta_custom_rbi, options.grad_RBI, options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['RBI-OSL custom prior sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
if any(options.COSEM_OSL) && options.custom
    if options.verbose
        tStart = tic;
    end
    if options.COSEM_OSL == 1
        [im_vectors.custom_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.custom_COSEM_apu, D, options.beta_custom_cosem, options.grad_COSEM, ...
            A, uu, options.epps, C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    else
        [im_vectors.custom_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.custom_COSEM_apu, D, options.beta_custom_cosem, options.grad_COSEM, ...
            A, uu, options.epps, C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    end
    if options.verbose
        tElapsed = toc(tStart);
        disp(['COSEM-OSL custom prior sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
fn = fieldnames(im_vectors);
for kk = 2 : 2 : numel(fn)
    im_vectors.(fn{kk})(im_vectors.(fn{kk}) < 0) = options.epps;
end
end

