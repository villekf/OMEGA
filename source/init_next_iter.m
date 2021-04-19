function im_vectors = init_next_iter(im_vectors, options, iter)
% Initialize the next iteration in the custom prior reconstruction

if options.save_iter
    iter_n = iter + 1;
else
    iter_n = 1;
end
if ~isfield(options,'N')
    options.N = options.Nx * options.Ny * options.Nz;
end

if options.OSEM
    im_vectors.OSEM(:, iter_n) = im_vectors.OSEM_apu;
end

if options.MLEM
    im_vectors.MLEM(:, iter_n) = im_vectors.MLEM_apu;
end

if options.MRAMLA
    im_vectors.MRAMLA(:, iter_n) = im_vectors.MRAMLA_apu;
end

if options.RAMLA
    im_vectors.RAMLA(:, iter_n) = im_vectors.RAMLA_apu;
end

if options.ROSEM
    im_vectors.ROSEM(:, iter_n) = im_vectors.ROSEM_apu;
end

if options.RBI
    im_vectors.RBI(:, iter_n) = im_vectors.RBI_apu;
end

if options.DRAMA
    im_vectors.DRAMA(:, iter_n) = im_vectors.DRAMA_apu;
end

if options.COSEM
    im_vectors.COSEM(:, iter_n) = im_vectors.COSEM_apu;
end

if options.ECOSEM
    im_vectors.ECOSEM(:, iter_n) = im_vectors.ECOSEM_apu;
end

if options.ACOSEM
    im_vectors.ACOSEM(:, iter_n) = im_vectors.ACOSEM_apu;
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

varapu = strcat(varMAP,'_apu');

for kk = 1 : numel(varMAP)
    for ll = 1 : numel(varPrior)
        if options.verbose
            tStart = tic;
        end
        if strcmp(varPrior{ll},'MRP') && (strcmp(varMAP{kk},'BSREM') || strcmp(varMAP{kk},'ROSEM_MAP'))
            grad = MRP(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, ...
                options.tr_offsets, options.med_no_norm);
        elseif strcmp(varPrior{ll},'quad') && (strcmp(varMAP{kk},'BSREM') || strcmp(varMAP{kk},'ROSEM_MAP'))
            grad = Quadratic_prior(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.weights_quad, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz);
        elseif strcmp(varPrior{ll},'Huber') && (strcmp(varMAP{kk},'BSREM') || strcmp(varMAP{kk},'ROSEM_MAP'))
            grad = Huber_prior(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.weights_huber, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, options.huber_delta);
        elseif strcmp(varPrior{ll},'FMH') && (strcmp(varMAP{kk},'BSREM') || strcmp(varMAP{kk},'ROSEM_MAP'))
            grad = FMH(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.tr_offsets, options.fmh_weights, options.Nx, options.Ny, options.Nz, options.N, options.Ndx, options.Ndy, options.Ndz, options.epps, ...
                options.med_no_norm);
        elseif strcmp(varPrior{ll},'L') && (strcmp(varMAP{kk},'BSREM') || strcmp(varMAP{kk},'ROSEM_MAP'))
            grad = L_filter(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, ...
                options.epps, options.med_no_norm);
        elseif strcmp(varPrior{ll},'weighted_mean') && (strcmp(varMAP{kk},'BSREM') || strcmp(varMAP{kk},'ROSEM_MAP'))
            grad = Weighted_mean(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.weighted_weights, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, ...
                options.mean_type, options.epps, options.med_no_norm);
        elseif strcmp(varPrior{ll},'TV') && (strcmp(varMAP{kk},'BSREM') || strcmp(varMAP{kk},'ROSEM_MAP'))
            grad = TVpriorFinal(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, ...
                options.tr_offsets);
        elseif strcmp(varPrior{ll},'AD') && (strcmp(varMAP{kk},'BSREM') || strcmp(varMAP{kk},'ROSEM_MAP'))
            grad = AD(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.FluxType, options.Nx, options.Ny, options.Nz, options);
        elseif strcmp(varPrior{ll},'APLS') && (strcmp(varMAP{kk},'BSREM') || strcmp(varMAP{kk},'ROSEM_MAP'))
            grad = TVpriorFinal(im_vectors.([varPrior{ll} '_' varapu{kk}]), [], options.Nx, options.Ny, options.Nz, true, options, 5);
        elseif strcmp(varPrior{ll},'TGV') && (strcmp(varMAP{kk},'BSREM') || strcmp(varMAP{kk},'ROSEM_MAP'))
            grad = TGV(im_vectors.([varPrior{ll} '_' varapu{kk}]),options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
        elseif strcmp(varPrior{ll},'NLM') && (strcmp(varMAP{kk},'BSREM') || strcmp(varMAP{kk},'ROSEM_MAP'))
            grad = NLM(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
        end
            
        if strcmp(varMAP{kk},'BSREM')
            im_vectors.([varPrior{ll} '_' varMAP{kk}])(:,iter_n) = BSREM_iter(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.lam, iter, options.(['beta_' varPrior{ll} '_' varMAP{kk}]), grad, options.epps);
            im_vectors.([varPrior{ll} '_' varapu{kk}]) = im_vectors.([varPrior{ll} '_' varMAP{kk}])(:,iter_n);
        elseif strcmp(varMAP{kk},'ROSEM_MAP')
            im_vectors.([varPrior{ll} '_' varMAP{kk}])(:,iter_n) = BSREM_iter(im_vectors.([varPrior{ll} '_' varapu{kk}]), options.lam_ROSEM, iter, options.(['beta_' varPrior{ll} '_' varMAP{kk}]), grad, options.epps);
            im_vectors.([varPrior{ll} '_' varapu{kk}]) = im_vectors.([varPrior{ll} '_' varMAP{kk}])(:,iter_n);
        else
            im_vectors.([varPrior{ll} '_' varMAP{kk}])(:,iter_n) = im_vectors.([varPrior{ll} '_' varapu{kk}]);
        end
        if options.verbose
            apu = varMAP{kk};
            apu = strrep(apu, '_', '-');
            tElapsed = toc(tStart);
            disp([apu ' ' varPrior{ll} ' iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
end

disp(['Iteration ' num2str(iter) ' finished'])