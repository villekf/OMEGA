function [im_vectors] = computeEstimatesImp4Iter(im_vectors, options, gaussK, iter, iter_n, N, Ndx, Ndy, Ndz, epps, Nx, Ny, Nz, osem, tStart, varargin)
%COMPUTEESTIMATESIMP1 Computes the subset estimates for implementation 1
%   Utility function

if osem
    
    % Compute OS-based algorithms
    
    varList = recNames(4);
    varML = recNames(6);
    
    varApu = 'OSEM';
else
    
    varList = recNames(3);
    varML = recNames(5);
    
    varApu = 'MLEM';
end

varPrior = recNames(1);
apu = false(numel(varList),1);
for kk = 1 : numel(varList)
    if options.(varList{kk})
        apu(kk) = true;
    end
end
varList = varList(apu);
ind = strcmp(varML,varList);
if any(ind)
    varMAP = [];
else
    varMAP = varList;
end

apu = false(numel(varPrior),1);
if numel(varMAP) > 0
    for kk = 1 : numel(varPrior)
        if options.(varPrior{kk}) && options.(varMAP{1})
            apu(kk) = true;
        end
    end
end
varPrior = varPrior(apu);
varPriorName = recNames(10);
varPriorName = varPriorName(apu);

if numel(varPrior) > 0
    varTot = [varPrior{1} '_' varList{1}];
else
    varTot = varList{1};
end

% for kk = 1 : numel(varList)
%     for ll = 1 : numel(varPrior)
% if options.verbose
%     tStart = tic;
% end
if numel(varPrior) > 0 && (strcmp(varList{1},'BSREM') || strcmp(varList{1},'ROSEM_MAP') || ~osem)
    if strcmp(varPrior{1},'MRP')
        grad = MRP(im_vectors.([varApu '_apu']), options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, ...
            options.tr_offsets, options.med_no_norm);
    elseif strcmp(varPrior{1},'quad')
        grad = Quadratic_prior(im_vectors.([varApu '_apu']), options.weights_quad, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz);
    elseif strcmp(varPrior{1},'Huber')
        grad = Huber_prior(im_vectors.([varApu '_apu']), options.weights_huber, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, options.huber_delta);
    elseif strcmp(varPrior{1},'FMH')
        grad = FMH(im_vectors.([varApu '_apu']), options.tr_offsets, options.fmh_weights, options.Nx, options.Ny, options.Nz, N, Ndx, Ndy, Ndz, options.epps, ...
            options.med_no_norm);
    elseif strcmp(varPrior{1},'L')
        grad = L_filter(im_vectors.([varApu '_apu']), options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, ...
            options.epps, options.med_no_norm);
    elseif strcmp(varPrior{1},'weighted_mean')
        grad = Weighted_mean(im_vectors.([varApu '_apu']), options.weighted_weights, options.Nx, options.Ny, options.Nz, Ndx, Ndy, Ndz, ...
            options.mean_type, options.epps, options.med_no_norm);
    elseif strcmp(varPrior{1},'TV')
        grad = TVpriorFinal(im_vectors.([varApu '_apu']), options.TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, ...
            options.tr_offsets);
    elseif strcmp(varPrior{1},'AD')
        grad = AD(im_vectors.([varApu '_apu']), options.FluxType, options.Nx, options.Ny, options.Nz, options);
    elseif strcmp(varPrior{1},'APLS')
        grad = TVpriorFinal(im_vectors.([varApu '_apu']), [], options.Nx, options.Ny, options.Nz, true, options, 5);
    elseif strcmp(varPrior{1},'TGV')
        grad = TGV(im_vectors.([varApu '_apu']),options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
    elseif strcmp(varPrior{1},'NLM')
        grad = NLM(im_vectors.([varApu '_apu']), options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
            options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
    end
end

if strcmp(varList{1},'BSREM')
    im_vectors.(varTot)(:,iter_n) = BSREM_iter(im_vectors.([varApu '_apu']), options.lam, iter, options.(['beta_' varTot]), grad, options.epps);
    im_vectors.([varApu '_apu']) = im_vectors.(varTot)(:,iter_n);
elseif strcmp(varList{1},'ROSEM_MAP')
    im_vectors.(varTot)(:,iter_n) = BSREM_iter(im_vectors.([varApu '_apu']), options.lam_ROSEM, iter, options.(['beta_' varTot]), grad, options.epps);
    im_vectors.([varApu '_apu']) = im_vectors.(varTot)(:,iter_n);
elseif ~osem
    % Non-subset-based algorithms
    if numel(varPrior) > 0
%         if strcmp(varList{1},'OSL_MLEM')
            im_vectors.([varApu '_apu']) = OSL_OSEM(im_vectors.([varApu '_apu']), varargin{1}, options.(['beta_' varTot]), grad, epps, varargin{2});
            im_vectors.(varTot)(:, iter_n) = im_vectors.([varApu '_apu']);
%         end
    else
%         if strcmp(varList{1},'MLEM')
            im_vectors.([varApu '_apu']) = MLEM_im(im_vectors.([varApu '_apu']), varargin{1}, epps, varargin{2});
            im_vectors.(varTot)(:, iter_n) = im_vectors.([varApu '_apu']);
%         end
    end
else
    im_vectors.(varTot)(:,iter_n) = im_vectors.([varApu '_apu']);
end
if options.verbose
    apu = varList{1};
    apu = strrep(apu, '_', '-');
    if numel(varPrior) > 0
        apu2 = varPriorName{1};
    else
        apu2 = '';
    end
    tElapsed = toc(tStart);
    disp([apu ' ' apu2 ' iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
end

end