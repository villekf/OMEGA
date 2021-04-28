function [options, D, C_co, C_aco, C_osl, Amin, E] = prepass_phase(options, pituus, index, SinM, pseudot, x, y, xx, yy, z_det, dz, dx, dy, bz, bx, by, NSlices, zmax, size_x, block1, blocks,...
    normalization_correction, randoms_correction, xy_index, z_index, lor_a, lor_orth, summa, LL, is_transposed, x_center, y_center, z_center, ind_size, gaussK, bmin, bmax, Vmax, V)
%PREPASS_PHASE Prepass step for various priors and algorithms
% Computes the necessary variables (e.g. weights) for certain
% algorithms/priors if they have been selected. Also converts various
% values to single precision if implementation 2 has been selected.

if options.sampling > 1 && ~options.use_raw_data && ~options.precompute_lor
    options.Ndist = options.Ndist * options.sampling;
end
Nx = options.Nx;
Ny = options.Ny;
Nz = options.Nz;
epps = options.epps;
attenuation_correction = options.attenuation_correction;
NSinos = options.NSinos;
det_per_ring = options.det_per_ring;
use_raw_data = options.use_raw_data;
verbose = options.verbose;
Ny = uint32(Ny);
Nx = uint32(Nx);
Nz = uint32(Nz);
N = (Nx)*(Ny)*(Nz);
det_per_ring = uint32(det_per_ring);
D = [];
C_co = [];
C_aco = [];
C_osl = [];
Amin = [];
E = [];

if exist('feature','builtin') == 5
    nCores = uint32(feature('numcores'));
else
    nCores = uint32(1);
end
if ~isfield(options,'listmode')
    options.listmode = false;
end
if ~isfield(options,'PKMA')
    options.PKMA = false;
end
if ~isfield(options,'CT')
    options.CT = false;
end

if (options.MRP || options.quad || options.Huber || options.TV ||options. FMH || options.L || options.weighted_mean || options.APLS || options.BSREM ...
        || options.RAMLA || options.MBSREM || options.MRAMLA || options.ROSEM || options.DRAMA || options.ROSEM_MAP || options.ECOSEM ...
        || options.COSEM || options.ACOSEM || options.AD || any(options.OSL_COSEM) || options.NLM || options.OSL_RBI || options.RBI || options.PKMA || options.RDP)
    
    % Compute and/or load necessary variables for the TV regularization
    if options.TV && options.MAP
        options = TVPrepass(options);
    end
    
    if options.TV && options.MAP && options.implementation == 2
        options.alphaTGV = single(options.alphaTGV);
        options.betaTGV = single(options.betaTGV);
        options.NiterTGV = uint32(options.NiterTGV);
    end
    
    % Load necessary variables for the APLS regularization
    if options.APLS && options.MAP
        options = APLSPrepass(options);
    end
    
    % Compute the necessary variables for MRAMLA, RBI-OSL and/or various
    % COSEM algorithms
    % E.g. for COSEM compute the complete data matrix, for RBI-OSL compute
    % the sum of all the rows of the system matrix
    if ((options.MRAMLA || options.MBSREM || options.OSL_RBI || options.RBI) && options.MBSREM_prepass || options.ECOSEM || options.COSEM ...
            || options.ACOSEM || any(options.OSL_COSEM) || options.PKMA)  && options.implementation == 1
        
        if options.ACOSEM
            C_aco = zeros(double(N), options.subsets);
        end
        if options.COSEM || options.ECOSEM
            C_co = zeros(double(N), options.subsets);
        end
        if any(options.OSL_COSEM)
            C_osl = zeros(double(N), options.subsets);
        end
        if options.ACOSEM || options.COSEM || options.ECOSEM || any(options.OSL_COSEM)
            if options.use_psf
                im_apu = computeConvolution(options.x0(:), options, Nx, Ny, Nz, gaussK);
            else
                im_apu = options.x0(:);
            end
            if options.ACOSEM || options.OSL_COSEM == 1
                im = power(options.x0(:), 1/options.h);
            end
        end
        if options.MRAMLA || options.MBSREM
            if options.precompute_lor == false
                Amin = zeros(options.Nang*options.Ndist*options.NSinos,1);
            else
                Amin = zeros(pituus(end),1);
            end
        end
        
        if ~use_raw_data
            if isempty(pseudot)
                pseudot = uint32(0);
            end
        end
        
        D = zeros(N,1);
        if options.precompute_lor
            if normalization_correction || options.attenuation_correction
                E = zeros(length(lor_a),1);
            else
                E = ones(length(lor_a),1);
            end
        else
            if normalization_correction || options.attenuation_correction
                E = zeros(options.Nang*options.Ndist*options.NSinos,1);
            else
                E = ones(options.Nang*options.Ndist*options.NSinos,1);
            end
        end
        if options.implementation == 1 && options.precompute_lor == false
            iij = double(0:Nx);
            jji = double(0:Ny);
            kkj = double(0:Nz);
        else
            iij = 0;
            jji = 0;
            kkj = 0;
        end
        
        if verbose
            disp('Prepass phase for MRAMLA, RBI, PKMA, COSEM, ACOSEM and ECOSEM started')
        end
        if iscell(SinM)
            Sino = SinM{1};
        else
            Sino = SinM;
        end
        
        Sino = Sino(:);
        
        if issparse(Sino)
            Sino = (full(Sino));
        end
        for osa_iter = 1 : options.subsets
            if randoms_correction
                if iscell(options.SinDelayed)
                    SinD = double(options.SinDelayed{1}(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                else
                    SinD = double(options.SinDelayed(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                end
                if issparse(SinD)
                    SinD = (full(SinD));
                end
                SinD = SinD(:);
            else
                SinD = 0;
            end
            if normalization_correction
                norm_input = options.normalization(pituus(osa_iter)+1:pituus(osa_iter + 1));
            else
                norm_input = 0;
            end
            if options.scatter_correction && ~options.subtract_scatter
                scatter_input = double(options.ScatterC(pituus(osa_iter)+1:pituus(osa_iter + 1)));
            else
                scatter_input = 0;
            end
            koko = pituus(osa_iter + 1) - pituus(osa_iter);
            [A] = computeImplementation1(options,use_raw_data,randoms_correction, pituus,osa_iter, normalization_correction,...
                Nx, Ny, Nz, dx, dy, dz, bx, by, bz, x, y, z_det, xx, yy, size_x, NSinos, NSlices, zmax, attenuation_correction, pseudot, det_per_ring, ...
                false, 0, 0, uint32(0), nCores, ind_size, block1, blocks, index, iij, jji, kkj, LL, N, summa, lor_a, xy_index, z_index, ...
                x_center, y_center, z_center, bmin, bmax, Vmax, V, lor_orth, gaussK, is_transposed, scatter_input, norm_input, SinD, koko);
            if is_transposed
                % Sensitivity image
                D = D + A * ones(size(A,2),1,'double');
                % Required for MRAMLA/MBSREM epsilon value
                if (normalization_correction || options.attenuation_correction) && (options.MRAMLA || options.MBSREM)
                    E(pituus(osa_iter)+1:pituus(osa_iter + 1)) = full(sum(A,1))';
                end
            else
                D = D + full(sum(A,1))';
                if normalization_correction || options.attenuation_correction && (options.MRAMLA || options.MBSREM)
                    E(pituus(osa_iter)+1:pituus(osa_iter + 1)) = full(sum(A,2))';
                end
            end
            if options.ECOSEM || options.COSEM || options.ACOSEM || any(options.OSL_COSEM)
                uu = double(Sino(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                if is_transposed
                    FP = A' * im_apu + epps + SinD;
                    if options.CT
                        FP = exp(FP);
                    end
                    RHS = A * (uu ./ FP);
                else
                    FP = A * im_apu + epps + SinD;
                    if options.CT
                        FP = exp(FP);
                    end
                    RHS = A' * (uu ./ FP);
                end
                if options.use_psf
                    RHS = computeConvolution(RHS, options, Nx, Ny, Nz, gaussK);
                end
            end
            if options.COSEM || options.ECOSEM
                if osa_iter > 1
                    if options.verbose
                        tic
                    end
                    C_co(:,osa_iter) = options.x0(:) .* RHS;
                    if options.verbose
                        disp(['COSEM complete data calculation took ' num2str(toc) ' seconds'])
                    end
                end
            end
            if options.ACOSEM
                if osa_iter > 1
                    if options.verbose
                        tic
                    end
                    C_aco(:,osa_iter) = im .* RHS;
                    if options.verbose
                        disp(['ACOSEM complete data calculation took ' num2str(toc) ' seconds'])
                    end
                end
            end
            if any(options.OSL_COSEM)
                if options.OSL_COSEM == 2
                    if osa_iter > 1
                        C_osl(:,osa_iter) = options.x0(:) .* RHS;
                    end
                else
                    if osa_iter > 1
                        C_osl(:,osa_iter) = im .* RHS;
                    end
                end
            end
            % Required for upper bound of MRAMLA/MBSREM
            if options.MBSREM_prepass && options.U == 0 && (options.MBSREM || options.MRAMLA)
                %%%% This particular piece of code was taken from:
                %%%% https://se.mathworks.com/matlabcentral/answers/35309-max-min-of-sparse-matrices
                if is_transposed
                    [~,m] = size(A);
                    rowMin = nan(m, 1);
                    [~,I,S] = find(A);
                else
                    [m,~] = size(A);
                    rowMin = nan(m, 1);
                    [I,~,S] = find(A);
                end
                I = I(S>1e-10);
                S = S(S>1e-10);
                [I,K] = sort(I);
                S = S(K);
                markers = [find([1; diff(I)]); numel(I)+1];
                iRows = uint32(I(markers(1:end-1)));
                for i = 1:numel(iRows)
                    s = S(markers(i):(markers(i+1)-1));
                    rowMin(iRows(i)) = min(s);
                end
                rowMin(isnan(rowMin)) = epps;
                if options.precompute_lor == false
                    if iscell(index)
                        Amin(index{osa_iter}) = rowMin;
                    else
                        Amin(index(pituus(osa_iter)+1:pituus(osa_iter + 1))) = rowMin;
                    end
                else
                    Amin(pituus(osa_iter)+1:pituus(osa_iter + 1)) = rowMin;
                end
                clear I K S markers rowMin s iRows
            end
            clear A
        end
        %             D = sum(pj,2);
        if verbose
            disp('Prepass phase for COSEM, ACOSEM and ECOSEM completed')
        end
    end
    
    % Lambda values (relaxation parameters)
    if (options.BSREM || options.RAMLA) && length(options.lambda0) == 1
        lam = zeros(options.Niter,1);
        lam(1) = options.lambda0;
        %             orig_lam = lam;
        %             if lam(1) > 1/max(max(pj))
        %                 lam(1) = min(min(pj));
        %             end
        for i=2:options.Niter
            %                 lam(i) = 0.5*lam(i-1);
            lam(i) = lam(1)/i;
            %                 lam(i) = lam(1)/1.01;
        end
        if options.implementation == 2
            options.lam = single(lam);
        else
            options.lam = lam;
        end
    elseif (options.BSREM || options.RAMLA)
        if numel(options.lambda0) < options.Niter
            error('The number of relaxation values needs to be at least equal to the number of iterations')
        end
        if options.implementation == 2
            options.lam = single(options.lambda0);
        else
            options.lam = double(options.lambda0);
        end
    end
    if (options.MBSREM || options.MRAMLA) && length(options.lambda0_MBSREM) == 1
        lam_MBSREM = zeros(options.Niter,1);
        lam_MBSREM(1) = options.lambda0_MBSREM;
        for i=2:options.Niter
            lam_MBSREM(i) = lam_MBSREM(1)/(i);
        end
        if options.implementation == 2
            options.lam_MBSREM = single(lam_MBSREM);
        else
            options.lam_MBSREM = lam_MBSREM;
        end
    elseif (options.MBSREM || options.MRAMLA)
        if numel(options.lambda0_MBSREM) < options.Niter
            error('The number of relaxation values needs to be at least equal to the number of iterations')
        end
        if options.implementation == 2
            options.lam_MBSREM = single(options.lambda0_MBSREM);
        else
            options.lam_MBSREM = double(options.lambda0_MBSREM);
        end
    end
    if (options.ROSEM_MAP || options.ROSEM) && length(options.lambda0_ROSEM) == 1
        lam_ROSEM = zeros(options.Niter,1);
        lam_ROSEM(1) = options.lambda0_ROSEM;
        for i=2:options.Niter
            lam_ROSEM(i) = lam_ROSEM(1)/i;
        end
        if options.implementation == 2
            options.lam_ROSEM = single(lam_ROSEM);
        else
            options.lam_ROSEM = lam_ROSEM;
        end
    elseif (options.ROSEM_MAP || options.ROSEM)
        if numel(options.lambda0_ROSEM) < options.Niter
            error('The number of relaxation values needs to be at least equal to the number of iterations')
        end
        if options.implementation == 2
            options.lam_ROSEM = single(options.lambda0_ROSEM);
        else
            options.lam_ROSEM = double(options.lambda0_ROSEM);
        end
    end
    if (options.PKMA) && length(options.lambda0_PKMA) == 1
%         warning('PKMA requires the relaxation parameter to be a vector, using the default values')
        lam_PKMA = zeros(options.Niter,1);
        for i = 1 : options.Niter
            lam_PKMA(i) = 1 / ((i - 1)/12 + 1);
        end
        if options.implementation == 2
            options.lam_PKMA = single(lam_PKMA);
        else
            options.lam_PKMA = lam_PKMA;
        end
    elseif (options.PKMA)
        if numel(options.lambda0_PKMA) < options.Niter
            error('The number of PKMA relaxation values needs to be at least equal to the number of iterations')
        end
        if options.implementation == 2
            options.lam_PKMA = single(options.lambda0_PKMA);
        else
            options.lam_PKMA = double(options.lambda0_PKMA);
        end
    end
    if options.DRAMA
        lam_drama = zeros(options.Niter,options.subsets);
        lam_drama(1,1) = options.beta_drama/(options.alpha_drama*options.beta0_drama);
        r = 1;
        for i=1:options.Niter
            for j = 1 : options.subsets
                lam_drama(i,j) = options.beta_drama/(options.alpha_drama*options.beta0_drama + r);
                r = r + 1;
            end
        end
        if options.implementation == 2
            options.lam_drama = single(lam_drama);
        else
            options.lam_drama = lam_drama;
        end
    end
    if (options.PKMA) && (~isfield(options,'alpha_PKMA') || numel(options.alpha_PKMA) < options.Niter * options.subsets)
        options.alpha_PKMA = zeros(options.Niter * options.subsets,1);
        oo = 1;
        for kk = 1 : options.Niter
            for ll = 0 : options.subsets - 1
                options.alpha_PKMA(oo) = 1 + (options.rho_PKMA *((kk - 1) * options.subsets + ll)) / ((kk - 1) * options.subsets + ll + options.delta_PKMA);
                oo = oo + 1;
            end
        end
        if options.implementation == 2
            options.alpha_PKMA = single(options.alpha_PKMA);
        end
    elseif (options.PKMA)
        if numel(options.alpha_PKMA) < options.Niter
            error('The number of PKMA alpha values must be at least the number of iterations!')
        end
        if options.implementation == 2
            options.alpha_PKMA = single(options.alpha_PKMA);
        end
    end
    if (options.PKMA) && (~isfield(options,'sigma_PKMA') || numel(options.sigma_PKMA) < options.Niter * options.subsets)
        options.alpha_PKMA = zeros(options.Niter * options.subsets,1);
        oo = 1;
        for kk = 1 : options.Niter
            for ll = 0 : options.subsets - 1
                options.alpha_PKMA(oo) = sqrt(((kk - 1) * options.subsets + ll)) / ((kk - 1) * options.subsets + ll + options.delta_PKMA);
                oo = oo + 1;
            end
        end
        if options.implementation == 2
            options.alpha_PKMA = single(options.alpha_PKMA);
        end
        if ~isfield(options,'sigma_PKMA') || sum(options.sigma_PKMA) == 0
            options.sigma_PKMA = ones(size(options.alpha_PKMA),class(options.alpha_PKMA));
        else
            options.sigma_PKMA = 1 - options.alpha_PKMA;
        end
    elseif (options.PKMA)
        if numel(options.sigma_PKMA) < options.Niter
            error('The number of PKMA sigma values must be at least the number of iterations!')
        end
        if options.implementation == 2
            options.sigma_PKMA = single(options.sigma_PKMA);
        end
    end
    % Compute the weights
    if (options.quad || options.L || options.FMH || options.weighted_mean || options.MRP || (options.TV && options.TVtype == 3) || options.Huber || options.RDP) && options.MAP
        if options.quad || options.L || options.FMH || options.weighted_mean || (options.TV && options.TVtype == 3) || options.Huber || options.RDP
            options = computeWeights(options);
        end
        % These values are needed in order to vectorize the calculation of
        % certain priors
        % Specifies the indices of the center pixel and its neighborhood
        if (options.MRP && ((options.implementation == 2 && options.use_CUDA) || ~license('test', 'image_toolbox'))) || options.L || options.FMH || ...
                (options.TV && options.TVtype == 3) || (options.RDP)
            options = computeOffsets(options);
        else
            options.tr_offsets = uint32(0);
            if options.implementation == 2 && isfield(options,'Ndx')
                options.Ndx = uint32(options.Ndx);
                options.Ndy = uint32(options.Ndy);
                options.Ndz = uint32(options.Ndz);
            end
            if options.MRP
                options.medx = options.Ndx*2 + 1;
                options.medy = options.Ndy*2 + 1;
                options.medz = options.Ndz*2 + 1;
            end
        end
        if options.quad || (options.TV && options.TVtype == 3)
            options = quadWeights(options, options.empty_weight);
        end
        if options.Huber
            options = huberWeights(options);
        end
        if options.RDP
            options = RDPWeights(options);
        end
        if options.L
            if isempty(options.a_L)
                options.a_L = lfilter_weights(options.Ndx, options.Ndy, options.Ndz, dx, dy, dz, options.oneD_weights);
            end
            if options.implementation == 2
                options.a_L = single(options.a_L);
            end
        end
        if options.FMH
            options = fmhWeights(options);
        end
        if (options.FMH || options.quad || options.Huber || options.RDP) && options.implementation == 2
            options.weights = single(options.weights);
            options.inffi = uint32(find(isinf(options.weights)) - 1);
            if isempty(options.inffi)
                options.inffi = uint32(floor(numel(options.weights) / 2));
            end
        end
        if options.weighted_mean
            options = weightedWeights(options);
        end
        if verbose
            disp('Prepass phase for MRP, quadratic prior, L-filter, FMH, RDP and weighted mean completed')
        end
    end
    if options.AD && options.MAP
        if options.implementation == 2
            options.NiterAD = uint32(options.NiterAD);
            options.KAD = single(options.KAD);
            options.TimeStepAD = single(options.TimeStepAD);
            options.FluxType = uint32(options.FluxType);
            options.DiffusionType = uint32(options.DiffusionType);
        else
            if options.FluxType == 1
                options.FluxType = 'exponential';
            elseif options.FluxType == 2
                options.FluxType = 'quadratic';
            end
        end
    end
    if options.NLM && options.MAP
        options = NLMPrepass(options);
    end
end
if options.sampling > 1 && ~options.use_raw_data && ~options.precompute_lor
    options.Ndist = options.Ndist / options.sampling;
end