function [options] = prepass_phase(options, dz, dx, dy)
%PREPASS_PHASE Prepass step for various priors and algorithms
% Computes the necessary variables (e.g. weights) for certain
% algorithms/priors if they have been selected. Also converts various
% values to single precision if implementation 2 has been selected.

if options.sampling > 1 && ~options.use_raw_data && ~options.precompute_lor
    options.Ndist = options.Ndist * options.sampling;
end
verbose = options.verbose;

if ~isfield(options,'listmode')
    options.listmode = false;
end
if ~isfield(options,'CT')
    options.CT = false;
end
options.Nf = options.nRowsD;
if options.precondTypeImage(3)
    if ischar(options.referenceImage)
        apu = load(options.referenceImage);
        variables = fieldnames(apu);
        options.referenceImage = apu.(variables{1});
    end
    if numel(options.referenceImage) == round((options.NxFull - options.NxOrig) * options.multiResolutionScale) * round((options.NyFull - options.NyOrig) * options.multiResolutionScale) * round((options.NzFull - options.NzOrig) * options.multiResolutionScale)
        skip = true;
    else
        false;
    end
    if skip == false && numel(options.referenceImage) ~= options.NxFull * options.NyFull * options.NzFull
        error('The size of the reference image does not match the reconstructed image!')
    end
    if (options.implementation == 2 || options.implementation == 5 || options.useSingles) && ~isa(options.referenceImage, 'single')
        options.referenceImage = single(options.referenceImage);
    end
    if skip == false && options.nMultiVolumes > 0
        options.referenceImage = reshape(options.referenceImage, options.NxFull, options.NyFull, options.NzFull);
    end
    if skip == false && options.nMultiVolumes == 6
        if exist('imresize3','file') == 2
            apu = imresize3(options.referenceImage,options.multiResolutionScale);
        else
            apu = imresize(options.referenceImage,options.multiResolutionScale);
            [XX, YY, ZZ] = meshgrid(1:size(apu,1), 1:size(apu,1), 1:1/options.multiResolutionScale:options.Nz);
            apu = interp3(apu, XX, YY, ZZ);
        end
        options.referenceImage = single(options.referenceImage(1 + (size(options.referenceImage,1) - options.NxOrig) / 2 : ...
            (size(options.referenceImage,1) - options.NxOrig) / 2 + options.NxOrig, ...
            1 + (size(options.referenceImage,2) - options.NyOrig) / 2 : ...
            (size(options.referenceImage,2) - options.NyOrig) / 2 + options.NyOrig,...
            1 + (size(options.referenceImage,3) - options.NzOrig) / 2 : ...
            (size(options.referenceImage,3) - options.NzOrig) / 2 + options.NzOrig));
        apu1 = apu(options.Nx(4) + 1 : options.Nx(4) + options.Nx(2), options.Ny(6) + 1 : options.Ny(6) + options.Ny(2), ...
            1 : options.Nz(2));
        options.referenceImage = [options.referenceImage(:);apu1(:)];
        apu1 = apu(options.Nx(5) + 1 : options.Nx(5) + options.Nx(3), options.Ny(7) + 1 : options.Ny(7) + options.Ny(3), ...
            end - options.Nz(2) + 1:end);
        options.referenceImage = [options.referenceImage(:);apu1(:)];
        apu1 = apu(1 : options.Nx(4), :, :);
        options.referenceImage = [options.referenceImage(:);apu1(:)];
        apu1 = apu(options.Nx(5) + options.Nx(3) + 1 : end, :, :);
        options.referenceImage = [options.referenceImage(:);apu1(:)];
        apu1 = apu(options.Nx(4) + 1 : options.Nx(4) + options.Nx(2), 1:options.Ny(6), ...
            1 : options.Nz(4));
        options.referenceImage = [options.referenceImage(:);apu1(:)];
        apu1 = apu(options.Nx(5) + 1 : options.Nx(5) + options.Nx(3), 1 + options.Ny(7) + options.Ny(3) : end, ...
            1 : options.Nz(5));
        options.referenceImage = [options.referenceImage(:);apu1(:)];
    elseif skip == false && options.nMultiVolumes == 4
        if exist('imresize3','file') == 2
            apu = single(imresize3(options.referenceImage,options.multiResolutionScale));
        else
            apu = single(imresize(options.referenceImage,options.multiResolutionScale));
            [XX, YY, ZZ] = meshgrid(1:size(apu,1), 1:size(apu,1), 1:1/options.multiResolutionScale:options.Nz);
            apu = interp3(apu, XX, YY, ZZ);
        end
        options.referenceImage = single(options.referenceImage(1 + (size(options.referenceImage,1) - options.NxOrig) / 2 : ...
            (size(options.referenceImage,1) - options.NxOrig) / 2 + options.NxOrig, ...
            1 + (size(options.referenceImage,2) - options.NyOrig) / 2 : ...
            (size(options.referenceImage,2) - options.NyOrig) / 2 + options.NyOrig, :));
        apu1 = apu(1 : options.Nx(2), :, :);
        options.referenceImage = [options.referenceImage(:);apu1(:)];
        apu1 = apu(options.Nx(2) + options.Nx(1) + 1 : end, :, :);
        options.referenceImage = [options.referenceImage(:);apu1(:)];
        apu1 = apu(options.Nx(2) + 1 : options.Nx(2) + options.Nx(4), 1:options.Ny(4), :);
        options.referenceImage = [options.referenceImage(:);apu1(:)];
        apu1 = apu(options.Nx(2) + 1 : options.Nx(2) + options.Nx(5), 1 + options.Ny(4) + options.Ny(1) : end, :);
        options.referenceImage = [options.referenceImage(:);apu1(:)];
    elseif skip == false && options.nMultiVolumes == 2
        if exist('imresize3','file') == 2
            apu = single(imresize3(options.referenceImage,options.multiResolutionScale));
        else
            apu = single(imresize(options.referenceImage,options.multiResolutionScale));
            [XX, YY, ZZ] = meshgrid(1:size(apu,1), 1:size(apu,1), 1:1/options.multiResolutionScale:options.Nz);
            apu = interp3(apu, XX, YY, ZZ);
        end
        options.referenceImage = single(options.referenceImage(:, :, ...
            1 + (size(options.referenceImage,3) - options.NzOrig) / 2 : ...
            (size(options.referenceImage,3) - options.NzOrig) / 2 + options.NzOrig));
        apu1 = apu(:, :, 1 : options.Nz(2));
        options.referenceImage = [options.referenceImage(:);apu1(:)];
        apu1 = apu(:, :, end - options.Nz(3) + 1:end);
        options.referenceImage = [options.referenceImage(:);apu1(:)];
    end
end

if (options.MRP || options.quad || options.Huber || options.TV ||options. FMH || options.L || options.weighted_mean || options.APLS || options.BSREM ...
        || options.RAMLA || options.MBSREM || options.MRAMLA || options.ROSEM || options.DRAMA || options.ROSEM_MAP || options.ECOSEM || options.SART || options.ASD_POCS ...
        || options.COSEM || options.ACOSEM || options.AD || any(options.OSL_COSEM) || options.NLM || options.OSL_RBI || options.RBI || options.PKMA...
        || options.RDP || options.SPS || options.ProxNLM || options.GGMRF ||options.hyperbolic)

    % Compute and/or load necessary variables for the TV regularization
    if options.TV && options.MAP
        options = TVPrepass(options);
    end

    %     if options.TV && options.MAP && options.implementation == 2
    %         options.alphaTGV = single(options.alphaTGV);
    %         options.betaTGV = single(options.betaTGV);
    %         options.NiterTGV = uint32(options.NiterTGV);
    %     end

    % Load necessary variables for the APLS regularization
    if options.APLS && options.MAP
        options = APLSPrepass(options);
    end

    if options.U == 0
        if options.CT
            options.U = 10;
        else
            options.U = 10000;
        end
    end

    % Lambda values (relaxation parameters)
    if (options.BSREM || options.RAMLA || options.MBSREM || options.MRAMLA || options.ROSEM_MAP || options.ROSEM || options.PKMA || options.SPS || options.SART || options.ASD_POCS) && (~isfield(options,'lambda') || isempty(options.lambda) || sum(options.lambda) == 0)
        lambda = zeros(options.Niter,1);
        for i = 1 : options.Niter
            lambda(i) = 1 / ((i - 1)/20 + 1);
        end
        if options.CT && ~options.SART && ~options.ASD_POCS
            lambda = lambda / 10000;
        end
        if options.implementation == 2 || options.useSingles || options.implementation == 5
            options.lambda = single(lambda);
        else
            options.lambda = lambda;
        end
    elseif (options.BSREM || options.RAMLA || options.MBSREM || options.MRAMLA || options.ROSEM_MAP || options.ROSEM || options.PKMA || options.SPS || options.SART || options.ASD_POCS) && numel(options.lambda) == 1 && options.Niter > 1
        lambdaT = zeros(options.Niter,1);
        for i = 1 : options.Niter
            lambdaT(i) = options.lambda / ((i - 1)/20 + 1);
        end
        if options.CT && ~options.SART && ~options.ASD_POCS
            lambdaT = lambdaT / 10000;
        end
        if options.implementation == 2 || options.useSingles || options.implementation == 5
            options.lambda = single(lambdaT);
        else
            options.lambda = lambdaT;
        end
    elseif (options.BSREM || options.RAMLA || options.MBSREM || options.MRAMLA || options.ROSEM_MAP || options.ROSEM || options.PKMA || options.SPS || options.SART || options.ASD_POCS)
        if numel(options.lambda) < options.Niter
            error('The number of relaxation values needs to be at least equal to the number of iterations!')
        end
        if numel(options.lambda) > options.Niter
            warning('The number of relaxation values is more than the number of iterations. Later values are ignored!')
        end
        if options.implementation == 2 || options.useSingles || options.implementation == 5
            options.lambda = single(options.lambda);
        else
            options.lambda = double(options.lambda);
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
        if options.implementation == 2 || options.useSingles || options.implementation == 5
            options.lam_drama = single(lam_drama);
        else
            options.lam_drama = lam_drama;
        end
    end
    if (options.PKMA) && (~isfield(options,'alpha_PKMA') || numel(options.alpha_PKMA) < options.Niter * options.subsets)
        if isfield(options,'alpha_PKMA') && numel(options.alpha_PKMA) < options.Niter * options.subsets
            warning('The number of PKMA alpha (momentum) values must be at least the number of iterations times the number of subsets! Computing custom alpha values.')
        end
        options.alpha_PKMA = zeros(options.Niter * options.subsets,1);
        oo = 1;
        for kk = 1 : options.Niter
            for ll = 0 : options.subsets - 1
                options.alpha_PKMA(oo) = 1 + (options.rho_PKMA *((kk - 1) * options.subsets + ll)) / ((kk - 1) * options.subsets + ll + options.delta_PKMA);
                oo = oo + 1;
            end
        end
        if options.implementation == 2 || options.useSingles || options.implementation == 5
            options.alpha_PKMA = single(options.alpha_PKMA);
        end
    elseif (options.PKMA)
        if numel(options.alpha_PKMA) > options.Niter * options.subsets
            warning('The number of PKMA alpha (momentum) values is higher than the total number of iterations times subsets. The final values will be ignored.')
        end
        if options.implementation == 2 || options.useSingles || options.implementation == 5
            options.alpha_PKMA = single(options.alpha_PKMA);
        end
    end
    % Compute the weights
    if (options.quad || options.L || options.FMH || options.weighted_mean || options.MRP || (options.TV && options.TVtype == 3 && options.TV_use_anatomical) || options.Huber || options.RDP || options.GGMRF ||options.hyperbolic) && options.MAP
        if options.quad || options.L || options.FMH || options.weighted_mean || (options.TV && options.TVtype == 3 && options.TV_use_anatomical) || options.Huber || options.RDP || options.GGMRF ||options.hyperbolic
            if options.GGMRF
                options = computeWeights(options, true);
            else
                options = computeWeights(options, false);
            end
        end
        % These values are needed in order to vectorize the calculation of
        % certain priors
        % Specifies the indices of the center pixel and its neighborhood
        if (options.MRP && (options.implementation ~= 2 && ~license('test', 'image_toolbox'))) || options.L || options.FMH || ...
                (options.TV && options.TVtype == 3 && options.TV_use_anatomical && options.implementation ~= 2) || (options.RDP && options.implementation ~= 2)
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
        if options.quad || options.hyperbolic || options.GGMRF
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
            if options.implementation == 2 || options.useSingles || options.implementation == 5
                options.a_L = single(options.a_L);
            end
        end
        if options.FMH
            options = fmhWeights(options);
        end
        if (options.FMH || options.quad || options.Huber) && options.implementation == 2
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
    if (options.NLM && options.MAP) || options.ProxNLM
        options = NLMPrepass(options);
    end
end
if options.precondTypeImage(4)
    if ~isfield(options,'alphaPrecond') || isempty(options.alphaPrecond)
        if isfield(options,'alpha_PKMA') && numel(options.alpha_PKMA) == options.Niter * options.subsets
            options.alphaPrecond = options.alpha_PKMA;
        end
        if numel(options.alphaPrecond) < options.Niter * options.subsets
            options.alphaPrecond = zeros(options.Niter * options.subsets,1);
            if ~isfield(options,'rhoPrecond')
                if isfield(options,'rho_PKMA')
                    options.rhoPrecond = options.rho_PKMA;
                else
                    error('rhoPrecond value is missing!')
                end
            end
            if ~isfield(options,'delta1Precond')
                if isfield(options,'delta_PKMA')
                    options.delta1Precond = options.delta_PKMA;
                else
                    error('delta1Precond value is missing!')
                end
            end
            % if ~isfield(options,'delta2Precond')
            %     if isfield(options,'delta2_PKMA')
            %         options.delta2Precond = options.delta2_PKMA;
            %     else
            %         error('delta2Precond value is missing!')
            %     end
            % end
            oo = 1;
            for kk = 0 : options.Niter - 1
                for ll = 1 : options.subsets
                    options.alphaPrecond(oo) = (options.rhoPrecond * (kk * options.subsets + ll - 1)) / (kk * options.subsets + ll - 1 + options.delta1Precond);
                    oo = oo + 1;
                end
            end
        end
    end
    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
        options.alphaPrecond = single(options.alphaPrecond);
    end
end
if options.PDHG || options.PDHGKL || options.PDHGL1 || options.PDDY
    if numel(options.thetaCP) < options.subsets * options.Niter
        if numel(options.thetaCP) > 1
            error('The number of elements in options.thetaCP has to be either one or options.subsets * options.Niter!')
        end
        options.thetaCP = repmat(options.thetaCP,options.subsets * options.Niter,1);
        if options.implementation == 2 || options.useSingles || options.implementation == 5
            options.thetaCP = single(options.thetaCP);
        end
    end
end
if (options.PKMA || options.MBSREM || options.SPS) && (options.ProxTV || options.TGV || options.ProxRDP)
    if ~isfield(options,'alpha_PKMA') && (~isfield(options,'thetaCP') || numel(options.thetaCP) < options.subsets * options.Niter)
        options.thetaCP = zeros(options.Niter * options.subsets,1);
        oo = 1;
        for kk = 1 : options.Niter
            for ll = 0 : options.subsets - 1
                options.thetaCP(oo) = 1 + (options.rho_PKMA *((kk - 1) * options.subsets + ll)) / ((kk - 1) * options.subsets + ll + options.delta_PKMA);
                oo = oo + 1;
            end
        end
    elseif isfield(options,'alpha_PKMA') && (~isfield(options,'thetaCP') || numel(options.thetaCP) < options.subsets * options.Niter)
        options.thetaCP = single(options.alpha_PKMA);
    end
end
if options.PDHG || options.PDHGKL || options.PDHGL1 || options.ProxTV || options.TGV || options.FISTA || options.FISTAL1 || options.ProxRDP || options.PDDY
    if numel(options.tauCP) < options.nMultiVolumes + 1
        options.tauCP = repelem(options.tauCP, options.nMultiVolumes + 1);
    end
    if numel(options.sigmaCP) < options.nMultiVolumes + 1
        options.sigmaCP = repelem(options.sigmaCP, options.nMultiVolumes + 1);
    end
    if numel(options.sigma2CP) < options.nMultiVolumes + 1
        options.sigma2CP = repelem(options.sigma2CP, options.nMultiVolumes + 1);
    end
    if numel(options.tauCPFilt) < options.nMultiVolumes + 1
        options.tauCPFilt = repelem(options.tauCPFilt, options.nMultiVolumes + 1);
    end
    if options.implementation == 1 || options.implementation == 4 || options.implementation == 5
        if options.filteringIterations > 0 && options.precondTypeMeas(2)
            apu = options.tauCP;
            options.tauCP = options.tauCPFilt;
            options.tauCPFilt = apu;
        end
    end
end

if options.precondTypeImage(6)
    %     if options.useZeroPadding
    options.Nf = double(2^nextpow2(options.Nx(1)));
    %     end
    if options.implementation == 2 || options.useSingles || options.implementation == 5
        options.filterIm = single(rampFilt(options.Nf,options.filterWindow, options.cutoffFrequency, options.normalFilterSigma, true));
    else
        options.filterIm = double(rampFilt(options.Nf,options.filterWindow, options.cutoffFrequency, options.normalFilterSigma, true));
    end
end
if options.precondTypeMeas(2)
    %     if options.useZeroPadding
    if options.subsets > 1 && options.subset_type == 5
        options.Nf = 2^nextpow2(options.nColsD);
    else
        options.Nf = 2^nextpow2(options.nRowsD);
    end
    %     end
    options.filter = rampFilt(options.Nf,options.filterWindow, options.cutoffFrequency, options.normalFilterSigma);
    options.filter(1) = 1e-6;
    options.Ffilter = ifft(options.filter) .* options.sigmaCP(1);
    if (options.PDHG || options.PDHGL1 || options.PDHGKL || options.CV || options.PDDY) && options.TGV
        options.Ffilter = options.Ffilter * options.beta;
    end
    options.Ffilter(1) = options.Ffilter(1) + 1;
    options.Ffilter = real(fft(options.Ffilter));
    if options.subsets > 1 && options.subset_type == 5
        options.filter2 = ifft(rampFilt(options.nColsD,options.filterWindow, options.cutoffFrequency, options.normalFilterSigma));
    else
        options.filter2 = ifft(rampFilt(options.nRowsD,options.filterWindow, options.cutoffFrequency, options.normalFilterSigma));
    end
    %     options.filter = 1 ./ options.filter;
end
if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
    options.thetaCP = single(options.thetaCP);
    options.tauCP = single(options.tauCP);
    options.tauCPFilt = single(options.tauCPFilt);
    options.sigmaCP = single(options.sigmaCP);
    options.sigma2CP = single(options.sigma2CP);
    if options.precondTypeMeas(2)
        options.filter = single(options.filter);
        options.filter2 = single(options.filter2);
        options.Ffilter = single(options.Ffilter);
    end
end
if options.sampling > 1 && ~options.use_raw_data && ~options.precompute_lor
    options.Ndist = options.Ndist / options.sampling;
end
if (options.PKMA || options.MBSREM || options.SPS || options.RAMLA || options.BSREM || options.ROSEM || options.ROSEM_MAP || options.MRAMLA) && (options.precondTypeMeas(2) || options.precondTypeImage(6))
    if ~isfield(options,'lambdaFiltered')
        options.lambdaFiltered = options.lambda;
    end
    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
        options.lambdaFiltered = single(options.lambdaFiltered);
    end
end
if options.FDK && options.CT && options.useFDKWeights
    % if options.largeDim
    %     [X, Y] = meshgrid(0:options.nColsD-1, 0:options.nRowsD-1);
    %     X = single(X - options.nColsD / 2);
    %     Y = single(Y - options.nRowsD / 2);
    %     for kk = 1 : options.nProjections
    %         if size(options.uV,1) == 2
    %             grid_points = options.x(4:6,kk)' + X(:) .* [options.uV(1,kk);options.dPitchX;options.dPitchX]' + Y(:) .* [options.dPitchX,options.uV(2,kk),options.dPitchX];
    %         else
    %             grid_points = options.x(4:6,kk)' + X(:) .* options.uV(1:3,kk)' + Y(:) .* options.uV(4:6,kk)';
    %         end
    %         options.FDKWeights = sqrt(sum((options.x(4:6,kk) - options.x(1:3,kk)).^2,2)) ./ squeeze(sqrt(sum((grid_points - options.x(1:3,kk)').^2, 2)));
    %         options.SinM(:,:,kk) = options.SinM(:,:,kk) .* reshape(options.FDKWeights, size(options.SinM,1), size(options.SinM,2));
    %     end
    % else
    %     [X, Y] = meshgrid(0:options.nColsD-1, 0:options.nRowsD-1);
    %     X = single(X - options.nColsD / 2);
    %     Y = single(Y - options.nRowsD / 2);
    %     grid_points = permute(options.x(4:6,:)',[3 2 1]) + X(:) .* permute(options.uV(1:3,:)',[3 2 1]) + Y(:) .* permute(options.uV(4:6,:)',[3 2 1]);
    %     options.FDKWeights = sqrt(sum((permute(options.x(4:6,:),[3 2 1]) - permute(options.x(1:3,:),[3 2 1])).^2,3)) ./ squeeze(sqrt(sum((grid_points - permute(options.x(1:3,:)',[3 2 1])).^2, 2)));
    %     options.FDKWeights = options.FDKWeights(:);
    %     options.SinM = options.SinM .* reshape(options.FDKWeights, size(options.SinM)) ./ 2;
    % end
    options.angles = single(options.angles);
    options.sourceToCRot = options.sourceToDetector;
    % options.sourceToCRot = mean(norm(options.x(1:3)));
    % distances = reshape(distances, options.nRowsD, options.nColsD, options.nProjections);
end