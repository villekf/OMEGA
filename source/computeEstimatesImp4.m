function [im_vectors,C_co,C_aco,C_osl] = computeEstimatesImp4(im_vectors, options, rhs, uu, f_Summ, SinD, gaussK, iter, osa_iter, C_co, C_aco,C_osl,...
    randoms_correction, N, Ndx, Ndy, Ndz, D, tStart, epps, use_raw_data, pituus,normalization_correction,...
    Nx, Ny, Nz, dx, dy, dz, bx, by, bz, x, y, z_det, xx, yy, size_x, NSinos, NSlices, zmax, attenuation_correction, pseudot, det_per_ring, ...
    TOF, TOFSize, sigma_x, TOFCenter, dec, nCores, L_input, lor_a_input, xy_index_input, z_index_input,  ...
    x_center, y_center, z_center, bmin, bmax, Vmax, V, scatter_input, norm_input, dc_z)
%COMPUTEESTIMATESIMP4 Computes the subset estimates for implementation 4
%   Utility function


varList = recNames(4);
varPrior = recNames(1);
varPriorName = recNames(10);

apu = false(numel(varList),1);
for kk = 1 : numel(varList)
    if options.(varList{kk})
        apu(kk) = true;
    end
end
varList = varList(apu);

apu = false(numel(varPrior),1);
for kk = 1 : numel(varPrior)
    if options.(varPrior{kk}) && options.(varList{1})
        apu(kk) = true;
    end
end
varPrior = varPrior(apu);
varPriorName = varPriorName(apu);

% for kk = 1 : numel(varList)
%     for ll = 1 : numel(varPrior)
% if options.verbose
%     tStart = tic;
% end
if numel(varPrior) > 0
    if strcmp(varPrior{1},'MRP') && ~strcmp(varList{1},'BSREM') && ~strcmp(varList{1},'ROSEM_MAP')
        grad = MRP(im_vectors.OSEM_apu, options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, ...
            options.tr_offsets, options.med_no_norm);
    elseif strcmp(varPrior{1},'quad') && ~strcmp(varList{1},'BSREM') && ~strcmp(varList{1},'ROSEM_MAP')
        grad = Quadratic_prior(im_vectors.OSEM_apu, options.weights_quad, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz);
    elseif strcmp(varPrior{1},'Huber') && ~strcmp(varList{1},'BSREM') && ~strcmp(varList{1},'ROSEM_MAP')
        grad = Huber_prior(im_vectors.OSEM_apu, options.weights_huber, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, options.huber_delta);
    elseif strcmp(varPrior{1},'FMH') && ~strcmp(varList{1},'BSREM') && ~strcmp(varList{1},'ROSEM_MAP')
        grad = FMH(im_vectors.OSEM_apu, options.tr_offsets, options.fmh_weights, options.Nx, options.Ny, options.Nz, N, Ndx, Ndy, Ndz, options.epps, ...
            options.med_no_norm);
    elseif strcmp(varPrior{1},'L') && ~strcmp(varList{1},'BSREM') && ~strcmp(varList{1},'ROSEM_MAP')
        grad = L_filter(im_vectors.OSEM_apu, options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, ...
            options.epps, options.med_no_norm);
    elseif strcmp(varPrior{1},'weighted_mean') && ~strcmp(varList{1},'BSREM') && ~strcmp(varList{1},'ROSEM_MAP')
        grad = Weighted_mean(im_vectors.OSEM_apu, options.weighted_weights, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, ...
            options.mean_type, options.epps, options.med_no_norm);
    elseif strcmp(varPrior{1},'TV') && ~strcmp(varList{1},'BSREM') && ~strcmp(varList{1},'ROSEM_MAP')
        grad = TVpriorFinal(im_vectors.OSEM_apu, options.TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, ...
            options.tr_offsets);
    elseif strcmp(varPrior{1},'AD') && ~strcmp(varList{1},'BSREM') && ~strcmp(varList{1},'ROSEM_MAP')
        if osa_iter > 1
            grad = AD(im_vectors.OSEM_apu, options.FluxType, options.Nx, options.Ny, options.Nz, options);
        else
            grad = zeros(options.Nx * options.Ny * options.Nz, 1);
        end
    elseif strcmp(varPrior{1},'APLS') && ~strcmp(varList{1},'BSREM') && ~strcmp(varList{1},'ROSEM_MAP')
        grad = TVpriorFinal(im_vectors.OSEM_apu, [], options.Nx, options.Ny, options.Nz, true, options, 5);
    elseif strcmp(varPrior{1},'TGV') && ~strcmp(varList{1},'BSREM') && ~strcmp(varList{1},'ROSEM_MAP')
        grad = TGV(im_vectors.OSEM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
    elseif strcmp(varPrior{1},'NLM') && ~strcmp(varList{1},'BSREM') && ~strcmp(varList{1},'ROSEM_MAP')
        grad = NLM(im_vectors.OSEM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
            options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
    elseif strcmp(varPrior{1},'RDP') && ~strcmp(varList{1},'BSREM') && ~strcmp(varList{1},'ROSEM_MAP')
        grad = RDP(im_vectors.OSEM_apu, options.weights_RDP, options.RDP_gamma, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, options.tr_offsets);
    end
end

if strcmp(varList{1},'OSEM')
    im_vectors.OSEM_apu = OSEM_im(im_vectors.OSEM_apu, rhs, f_Summ(:,osa_iter));
elseif strcmp(varList{1},'RAMLA')
    im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, f_Summ(:,osa_iter), rhs);
elseif strcmp(varList{1},'MRAMLA')
    error('MRAMLA is not supported when using implementation 4')
elseif strcmp(varList{1},'ROSEM')
    im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, options.lam_ROSEM, iter, f_Summ(:,osa_iter), epps, rhs);
elseif strcmp(varList{1},'RBI')
    im_vectors.OSEM_apu = RBI_subiter(im_vectors.OSEM_apu, f_Summ(:,osa_iter), rhs, [], [], D, [], []);
elseif strcmp(varList{1},'DRAMA')
    im_vectors.OSEM_apu = DRAMA_subiter(im_vectors.OSEM_apu, options.lam_drama, epps, iter, f_Summ(:,osa_iter), osa_iter, rhs);
elseif strcmp(varList{1},'COSEM')
    [im_vectors.OSEM_apu, C_co] = COSEM_im(im_vectors.OSEM_apu, rhs, C_co, D, osa_iter, [], [], [], []);
elseif strcmp(varList{1},'ECOSEM')
    no_norm_ecosem = true;
    if options.use_psf
        OSEM_apu = computeConvolution(im_vectors.COSEM_apu, options, Nx, Ny, Nz, gaussK);
    else
        OSEM_apu = im_vectors.COSEM_apu;
    end
    [~,rhs] = computeImplementation4(options,use_raw_data,randoms_correction, pituus,osa_iter, normalization_correction,...
        Nx, Ny, Nz, dx, dy, dz, bx, by, bz, x, y, z_det, xx, yy, size_x, NSinos, NSlices, zmax, attenuation_correction, pseudot, det_per_ring, ...
        TOF, TOFSize, sigma_x, TOFCenter, dec, nCores, L_input, lor_a_input, xy_index_input, z_index_input, epps, uu, OSEM_apu, no_norm_ecosem, ...
        x_center, y_center, z_center, bmin, bmax, Vmax, V, scatter_input, norm_input, SinD, dc_z);
    [im_vectors.COSEM_apu, C_co] = COSEM_im(im_vectors.COSEM_apu, rhs, C_co, D, osa_iter, [], [], [], []);
    
    if options.use_psf
        OSEM_apu = computeConvolution(im_vectors.OSEM_apu2, options, Nx, Ny, Nz, gaussK);
    else
        OSEM_apu = im_vectors.OSEM_apu2;
    end
    [~,rhs] = computeImplementation4(options,use_raw_data,randoms_correction, pituus,osa_iter, normalization_correction,...
        Nx, Ny, Nz, dx, dy, dz, bx, by, bz, x, y, z_det, xx, yy, size_x, NSinos, NSlices, zmax, attenuation_correction, pseudot, det_per_ring, ...
        TOF, TOFSize, sigma_x, TOFCenter, dec, nCores, L_input, lor_a_input, xy_index_input, z_index_input, epps, uu, OSEM_apu, no_norm_ecosem, ...
        x_center, y_center, z_center, bmin, bmax, Vmax, V, scatter_input, norm_input, SinD, dc_z);
    im_vectors.OSEM_apu2 = OSEM_im(im_vectors.OSEM_apu2, rhs, f_Summ(:,osa_iter));
    im_vectors.OSEM_apu = ECOSEM_im(im_vectors.OSEM_apu, epps, D, im_vectors.COSEM_apu, im_vectors.OSEM_apu2);
elseif strcmp(varList{1},'ACOSEM')
    [im_vectors.OSEM_apu, C_aco] = ACOSEM_im(im_vectors.OSEM_apu, rhs, C_aco, D, osa_iter, options.h, [], [], [], []);
    im_vectors = ACOSEM_coefficient(im_vectors, options, Nx, Ny, Nz, dx, dz, dy, bx, by, bz, z_det, x, y, yy, xx, NSinos, NSlices, size_x, ...
        zmax, norm_input, SinD, pituus, osa_iter, attenuation_correction, normalization_correction, randoms_correction, options.scatter, scatter_input, lor_a_input, ...
        xy_index_input, z_index_input, L_input, pseudot, det_per_ring, use_raw_data, epps, uu, dc_z, x_center, y_center, z_center, dec, bmin, ...
        bmax, Vmax, V, gaussK, TOF, TOFSize, sigma_x, TOFCenter, options.listmode);
elseif strcmp(varList{1},'OSL_OSEM')
    im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.(['beta_' varPrior{1} '_' varList{1}]), grad, epps, rhs);
elseif strcmp(varList{1},'BSREM')
    im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, f_Summ(:,osa_iter), rhs);
elseif strcmp(varList{1},'MBSREM')
    error('MBSREM is not supported when using implementation 4')
elseif strcmp(varList{1},'ROSEM_MAP')
    im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, options.lam_ROSEM, iter, f_Summ(:,osa_iter), epps, rhs);
elseif strcmp(varList{1},'OSL_RBI')
    im_vectors.OSEM_apu = RBI_subiter(im_vectors.OSEM_apu, f_Summ(:,osa_iter), rhs, [], [], D, [], [], options.(['beta_' varPrior{1} '_' varList{1}]), grad);
elseif strcmp(varList{1},'OSL_COSEM')
    [im_vectors.OSEM_apu, C_osl] = COSEM_OSL(im_vectors.OSEM_apu, D, options.(['beta_' varPrior{1} '_' varList{1}]), grad, rhs, osa_iter, options.h, ...
        C_osl, options.COSEM_OSL, [], [], [], []);
    if options.COSEM_OSL == 1
        im_vectors = ACOSEM_coefficient(im_vectors, options, Nx, Ny, Nz, dx, dz, dy, bx, by, bz, z_det, x, y, yy, xx, NSinos, NSlices, size_x, ...
            zmax, norm_input, SinD, pituus, osa_iter, attenuation_correction, normalization_correction, randoms_correction, options.scatter, scatter_input, lor_a_input, ...
            xy_index_input, z_index_input, L_input, pseudot, det_per_ring, use_raw_data, epps, uu, dc_z, x_center, y_center, z_center, dec, bmin, ...
            bmax, Vmax, V, gaussK, TOF, TOFSize, sigma_x, TOFCenter, options.listmode);
    end
elseif strcmp(varList{1},'PKMA')
    im_vectors.OSEM_apu = PKMA(im_vectors.OSEM_apu, options, rhs, f_Summ, options.lam_PKMA, options.alpha_PKMA, options.sigma_PKMA, options.epps, grad, ...
        options.(['beta_' varPrior{1} '_' varList{1}]), iter, osa_iter);
end
if options.verbose
    apu = varList{1};
    apu = strrep(apu, '_', '-');
    tElapsed = toc(tStart);
    if numel(varPrior) > 0
        disp([apu ' ' varPriorName{1} ' sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    else
        disp([apu ' sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
%     end
% end

im_vectors.OSEM_apu(im_vectors.OSEM_apu < epps) = epps;
% fn = fieldnames(im_vectors);
% for kk = 2 : 2 : numel(fn)
%     im_vectors.(fn{1})(im_vectors.(fn{1}) < epps) = epps;
% end
end

