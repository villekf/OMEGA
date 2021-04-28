function [im_vectors,C_co,C_aco,C_osl] = computeEstimatesImp1(im_vectors, options, A, uu, Summ, SinD, is_transposed, gaussK, iter, osa_iter, C_co, C_aco,C_osl,...
    randoms_correction, N, Ndx, Ndy, Ndz, D)
%COMPUTEESTIMATESIMP1 Computes the subset estimates for implementation 1
%   Utility function
% Compute OSEM
if options.OSEM || options.ECOSEM || options.attenuation_phase
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
if options.MRAMLA
    if options.verbose
        tStart = tic;
    end
    im_vectors.MRAMLA_apu = MBSREM(im_vectors.MRAMLA_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_MBSREM, iter, ...
        SinD, randoms_correction, is_transposed, [], [], options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['MRAMLA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute RAMLA
if options.RAMLA
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
if options.ROSEM
    if options.verbose
        tStart = tic;
    end
    im_vectors.ROSEM_apu = ROSEM_subiter(im_vectors.ROSEM_apu, options.lam_ROSEM, iter, Summ, options.epps, A, uu, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
    if options.verbose
        tElapsed = toc(tStart);
        disp(['ROSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
    end
end
% Compute RBI
if options.RBI
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
if options.DRAMA
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
if options.COSEM || options.ECOSEM
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
if options.ECOSEM
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
if options.ACOSEM
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

varMAP = recNames(2);
varPrior = recNames(1);

apu = false(numel(varMAP),1);
for kk = 1 : numel(varMAP)
    if options.(varMAP{kk})
        apu(kk) = true;
    end
end
varMAP = varMAP(apu);

apu = false(numel(varPrior),1);
for kk = 1 : numel(varPrior)
    if options.(varPrior{kk})
        apu(kk) = true;
    end
end
varPrior = varPrior(apu);
varPriorName = recNames(10);
varPriorName = varPriorName(apu);

varapu = strcat(varMAP,'_apu');

for kk = 1 : numel(varMAP)
    for ll = 1 : numel(varPrior)
        if options.verbose
            tStart = tic;
        end
        if strcmp(varPrior{ll},'MRP') && ~strcmp(varMAP{kk},'BSREM') && ~strcmp(varMAP{kk},'ROSEM_MAP')
            grad = MRP(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, ...
                options.tr_offsets, options.med_no_norm);
        elseif strcmp(varPrior{ll},'quad') && ~strcmp(varMAP{kk},'BSREM') && ~strcmp(varMAP{kk},'ROSEM_MAP')
            grad = Quadratic_prior(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.weights_quad, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz);
        elseif strcmp(varPrior{ll},'Huber') && ~strcmp(varMAP{kk},'BSREM') && ~strcmp(varMAP{kk},'ROSEM_MAP')
            grad = Huber_prior(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.weights_huber, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, options.huber_delta);
        elseif strcmp(varPrior{ll},'FMH') && ~strcmp(varMAP{kk},'BSREM') && ~strcmp(varMAP{kk},'ROSEM_MAP')
            grad = FMH(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.tr_offsets, options.fmh_weights, options.Nx, options.Ny, options.Nz, N, Ndx, Ndy, Ndz, options.epps, ...
                options.med_no_norm);
        elseif strcmp(varPrior{ll},'L') && ~strcmp(varMAP{kk},'BSREM') && ~strcmp(varMAP{kk},'ROSEM_MAP')
            grad = L_filter(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, ...
                options.epps, options.med_no_norm);
        elseif strcmp(varPrior{ll},'weighted_mean') && ~strcmp(varMAP{kk},'BSREM') && ~strcmp(varMAP{kk},'ROSEM_MAP')
            grad = Weighted_mean(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.weighted_weights, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, ...
                options.mean_type, options.epps, options.med_no_norm);
        elseif strcmp(varPrior{ll},'TV') && ~strcmp(varMAP{kk},'BSREM') && ~strcmp(varMAP{kk},'ROSEM_MAP')
            grad = TVpriorFinal(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, ...
                options.tr_offsets);
        elseif strcmp(varPrior{ll},'AD') && ~strcmp(varMAP{kk},'BSREM') && ~strcmp(varMAP{kk},'ROSEM_MAP')
            if osa_iter > 1
                grad = AD(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.FluxType, options.Nx, options.Ny, options.Nz, options);
            else
                grad = zeros(options.Nx * options.Ny * options.Nz, 1);
            end
        elseif strcmp(varPrior{ll},'APLS') && ~strcmp(varMAP{kk},'BSREM') && ~strcmp(varMAP{kk},'ROSEM_MAP')
            grad = TVpriorFinal(im_vectors.([varPrior{ll} '_' varapu{kk}]), [], options.Nx, options.Ny, options.Nz, true, options, 5);
        elseif strcmp(varPrior{ll},'TGV') && ~strcmp(varMAP{kk},'BSREM') && ~strcmp(varMAP{kk},'ROSEM_MAP')
            grad = TGV(im_vectors.([varPrior{ll} '_' varapu{kk}]),options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
        elseif strcmp(varPrior{ll},'NLM') && ~strcmp(varMAP{kk},'BSREM') && ~strcmp(varMAP{kk},'ROSEM_MAP')
            grad = NLM(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
        elseif strcmp(varPrior{ll},'RDP') && ~strcmp(varMAP{kk},'BSREM') && ~strcmp(varMAP{kk},'ROSEM_MAP')
            grad = RDP(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.weights_RDP, options.RDP_gamma, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, options.tr_offsets);
        end
            
        if strcmp(varMAP{kk},'OSL_OSEM')
            im_vectors.([varPrior{ll} '_' varapu{kk}]) = OSEM_im(im_vectors.([varPrior{ll} '_' varapu{kk}]), A, options.epps, uu, ...
                OSL(Summ, options.(['beta_' varPrior{ll} '_' varMAP{kk}]), grad, options.epps), SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
        elseif strcmp(varMAP{kk},'BSREM')
            im_vectors.([varPrior{ll} '_' varapu{kk}]) = BSREM_subiter(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.lam, options.epps, iter, A, uu, SinD, ...
                is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
        elseif strcmp(varMAP{kk},'MBSREM')
            im_vectors.([varPrior{ll} '_' varapu{kk}]) = MBSREM(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.U, options.pj3, A, options.epps, uu, ...
                options.epsilon_mramla, options.lam_MBSREM, iter, SinD, randoms_correction, is_transposed, options.(['beta_' varPrior{ll} '_' varMAP{kk}]), grad, ...
                options, options.Nx, options.Ny, options.Nz, gaussK);
        elseif strcmp(varMAP{kk},'ROSEM_MAP')
            im_vectors.([varPrior{ll} '_' varapu{kk}]) = ROSEM_subiter(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.lam_ROSEM, iter, Summ, options.epps, A, ...
                uu, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
        elseif strcmp(varMAP{kk},'OSL_RBI')
            im_vectors.([varPrior{ll} '_' varapu{kk}]) = RBI_subiter(im_vectors.([varPrior{ll} '_' varapu{kk}]), A, uu, options.epps, Summ, D, SinD, is_transposed, ...
                options.(['beta_' varPrior{ll} '_' varMAP{kk}]), grad, options, options.Nx, options.Ny, options.Nz, gaussK);
        elseif strcmp(varMAP{kk},'OSL_COSEM')
            if options.OSL_COSEM == 1
                [im_vectors.([varPrior{ll} '_' varapu{kk}]), C_osl] = COSEM_OSL(im_vectors.([varPrior{ll} '_' varapu{kk}]), D, options.(['beta_' varPrior{ll} '_' varMAP{kk}]), ...
                    grad, A, uu, options.epps, C_osl, options.h, options.OSL_COSEM, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
            else
                [im_vectors.([varPrior{ll} '_' varapu{kk}]), C_osl] = COSEM_OSL(im_vectors.([varPrior{ll} '_' varapu{kk}]), D, options.(['beta_' varPrior{ll} '_' varMAP{kk}]), ...
                    grad, A, uu, options.epps, C_osl, 0, options.OSL_COSEM, osa_iter, SinD, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
            end
        elseif strcmp(varMAP{kk},'PKMA')
            im_vectors.([varPrior{ll} '_' varapu{kk}]) = PKMA(im_vectors.([varPrior{ll} '_' varapu{kk}]), options, A, Summ, options.lam_PKMA, options.alpha_PKMA, ...
                options.sigma_PKMA, options.epps, grad, options.(['beta_' varPrior{ll} '_' varMAP{kk}]), iter, osa_iter, uu, is_transposed, SinD, options.Nx, options.Ny, ...
                options.Nz, gaussK);
        end
        if options.verbose
            apu = varMAP{kk};
            apu = strrep(apu, '_', '-');
            tElapsed = toc(tStart);
            disp([apu ' ' varPriorName{ll} ' sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
end
fn = fieldnames(im_vectors);
fn = fn(~cellfun('isempty',strfind(fn,'apu')));
for kk = 1 : numel(fn)
    im_vectors.(fn{kk})(im_vectors.(fn{kk}) < options.epps) = options.epps;
end
end

