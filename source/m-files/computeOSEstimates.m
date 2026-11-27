function [im_vectors, options] = computeOSEstimates(im_vectors, options, iter, osa_iter, uu, corrim, varargin)
%COMPUTEOSESTIMATES Subiteration update of the m-file algorithms (computeOSEstimates in subiterStep.h)
%   iter is the outer iteration and osa_iter the subset (sub-iteration)
%   number, both 1-based. The optional seventh input is the randoms/scatter
%   data of the current subset (used by ACOSEM and OSL-COSEM type 1)

% Global 1-based subiteration number, i.e. the C++ counter osa_iter + subsets * iter plus one
kk = osa_iter + options.param.subsets * (iter - 1);

% Image-based filter switch, done once for all volumes before the updates (C++: 0-based subiteration == filterIter)
if (options.param.precondTypeImage(6) && options.param.filteringIterations > 0 && kk - 1 == options.param.filteringIterations)
    if (options.param.verbose >= 3)
        disp("Image-based filter iterations complete.");
    end
    if (options.param.CPType || options.param.FISTA || options.param.FISTAL1 || options.param.TGV || options.param.ProxTV)
        for ii = 1 : options.param.nMultiVolumes + 1
            if (options.param.sigmaCP(ii) == 1)
                options.param.tauCP(ii) = options.param.tauCP2(ii);
            elseif (options.param.sigmaCP(ii) == options.param.tauCP(ii))
                options.param.tauCP(ii) = options.param.tauCP2(ii);
                options.param.sigmaCP(ii) = options.param.tauCP2(ii);
            else
                options.param.sigmaCP(ii) = options.param.tauCP2(ii);
            end
        end
    end
    % lambda holds the filter-phase relaxation parameters and lambdaFiltered the ones after the filtering (swapped in
    % prepass_phase as in the C++ code)
    if (options.param.MRAMLA || options.param.MBSREM || options.param.SPS || options.param.RAMLA || options.param.BSREM || options.param.ROSEM || options.param.ROSEM_MAP || options.param.PKMA)
        options.param.lambda = options.param.lambdaFiltered;
    end
    options.param.precondTypeImage(6) = false;
end

for ii = 1 : options.param.nMultiVolumes + 1
    if (options.param.PDHG || options.param.PDHGKL || options.param.PDHGL1 || options.param.PDDY)
        [im_vectors] = PDHG1(options.param, im_vectors, kk, ii);
    end
    % PDDY (Salim, Condat, Mishchenko and Richtarik, JOTA 2022), one subset. p = uCP is the auxiliary primal variable and
    % the image estimate s = recApu was used for the forward projection (dual update). Here rhs = A^T y+ (PDHG1 does
    % nothing without subsets). The point x = p - tau * A^T y+ is where the prior gradient is evaluated, so it is placed
    % into recApu. PDHG2 then forms s+ = prox(x - tau * (A^T y+ + grad R(x))) and p+ = p + s+ - x
    if (options.param.PDDY)
        if (options.param.verbose == 3)
            disp("Computing PDDY step");
        end
        if iscell(im_vectors.recApu)
            im_vectors.recApu{ii} = im_vectors.uCP{ii} - options.param.tauCP(ii) .* im_vectors.rhs{ii};
        else
            im_vectors.recApu = im_vectors.uCP - options.param.tauCP(ii) .* im_vectors.rhs;
        end
    end
end

% Subset scaling of the prior strength. These methods combine one subset's data gradient (about 1/S of the full gradient)
% with the prior on every subiteration, so beta/S per subiteration gives the objective L + beta*R independently of the
% number of subsets. PDHG (the dual sums over the subsets), OSL-COSEM (full sensitivity) and BSREM/ROSEM-MAP (prior once
% per iteration) are not scaled
if (options.param.OSL_OSEM || options.param.OSL_RBI || options.param.PKMA || options.param.MBSREM || options.param.MRAMLA || options.param.SPS)
    subsetScale = 1 / options.param.subsets;
else
    subsetScale = 1;
end
% The prior gradient is recomputed only every regEveryIter-th (sub)iteration and the last computed gradient
% (im_vectors.dU) is used in between. The (0-based) counter runs from 0 to subsets * Niter - 1; the first and last are
% always computed. BSREM and ROSEM-MAP apply the prior once per iteration (init_next_iter)
if isfield(options.param, 'regEveryIter') && ~isempty(options.param.regEveryIter)
    regEveryIter = double(options.param.regEveryIter);
else
    regEveryIter = 1;
end
regCounter = kk - 1;
regCounterMax = options.param.subsets * options.param.Niter - 1;
computeReg = regEveryIter <= 1 || regCounter == 0 || regCounter == regCounterMax || mod(regCounter, regEveryIter) == 0;
if (options.param.BSREM || options.param.ROSEM_MAP)
    dU = [];
elseif computeReg || ~isfield(im_vectors, 'dU')
    % Prior gradient of 1-subset PDHG is evaluated at the primal iterate uCP, not at the extrapolated image (MATLAB has no
    % dual-split prox priors here, ProxTV errors in applyPrior). PDDY evaluates it at its own point x (recApu, see above)
    priorAtU = (options.param.PDHG || options.param.PDHGKL || options.param.PDHGL1) && ~options.param.PDDY && options.param.subsets == 1 && isfield(im_vectors, 'uCP');
    if priorAtU
        priorIm = im_vectors.uCP;
    else
        priorIm = im_vectors.recApu;
    end
    % The global subiteration number is passed (the AD prior is skipped on the very first subiteration)
    if iscell(priorIm)
        dU = applyPrior(priorIm{1}, options.param, options.param.beta * subsetScale, kk);
    else
        dU = applyPrior(priorIm, options.param, options.param.beta * subsetScale, kk);
    end
    im_vectors.dU = dU;
else
    dU = im_vectors.dU;
end

% Sum of the additive randoms/scatter term of the subset for the ACOSEM weight (the randoms are repeated for every TOF bin,
% see computeForwardStep)
randSum = [];
if ((options.param.ACOSEM || options.param.OSL_COSEM == 1) && ~options.param.CT && options.param.randoms_correction && ~isempty(varargin) && ~isempty(varargin{1}))
    randSum = sum(varargin{1}(:));
    if options.param.TOF && numel(varargin{1}) < numel(uu)
        randSum = randSum * options.param.TOF_bins;
    end
end

for ii = 1 : options.param.nMultiVolumes + 1
    if iscell(im_vectors.Sens)
        if size(im_vectors.Sens{ii},2) > 1
            Sens = im_vectors.Sens{ii}(:, osa_iter);
        else
            Sens = im_vectors.Sens{ii};
        end
    else
        if size(im_vectors.Sens,2) > 1
            Sens = im_vectors.Sens(:, osa_iter);
        else
            Sens = im_vectors.Sens;
        end
    end
    if iscell(im_vectors.recApu)
        im = im_vectors.recApu{ii};
        rhs = im_vectors.rhs{ii};
        if isfield(options.param, 'D')
            D = options.param.D{ii};
        end
    else
        im = im_vectors.recApu;
        rhs = im_vectors.rhs;
        if isfield(options.param, 'D')
            D = options.param.D;
        end
    end
    % For list-mode data the -1 term of the Poisson gradient is replaced by the sensitivity image D
    % (computeForwardStep forms only BP(y / ybar) over the events), as in subiterStep.h
    if options.param.listmode > 0 && isfield(options.param, 'D')
        if options.param.PKMA
            rhs = D - rhs;
        elseif options.param.RAMLA || options.param.MRAMLA || options.param.BSREM || options.param.MBSREM || options.param.SPS ...
                || options.param.DRAMA || options.param.RBI || options.param.OSL_RBI
            rhs = rhs - D;
        end
    end
    % The prior is applied to the main volume only, the extra (multi-resolution) volumes use the non-regularized update
    if ii == 1
        dUi = dU;
    else
        dUi = [];
    end

    % Compute the (matrix free) algorithms
    % Ordered Subsets Expectation Maximization (OSEM)
    if (options.param.OSEM || options.param.ECOSEM)
        if (options.param.verbose >= 3)
            disp("Computing OSEM/ECOSEM");
        end
        if (options.param.ECOSEM)
            OSEMApu = EM(im, Sens, rhs + options.param.epps);
        else
            if (options.param.CT)
                im = EM(im, Sens, options.param.flat * rhs);
            else
                im = EM(im, Sens, rhs + options.param.epps);
            end
        end
    end

    % Modfied Row-action Maximum Likelihood (MRAMLA)
    if (options.param.MRAMLA)
        if (options.param.verbose >= 3)
            disp("Computing MRAMLA");
        end
        [im, options.param] = MBSREM(im, rhs, options.param, options.param.U, options.param.lambda, iter, osa_iter, ii);
    end

    % Row-action Maximum Likelihood (RAMLA)
    if (options.param.RAMLA)
        if (options.param.verbose >= 3)
            disp("Computing RAMLA");
        end
        im = BSREM(im, rhs, options.param.lambda, iter, options.param.epps, options.param.enforcePositivity);
    end

    % Relaxed OSEM (ROSEM)
    if (options.param.ROSEM)
        if (options.param.verbose >= 3)
            disp("Computing ROSEM");
        end
        if (options.param.CT)
            im = ROSEM(im, Sens, options.param.flat * rhs, options.param.lambda, iter);
        else
            im = ROSEM(im, Sens, rhs, options.param.lambda, iter);
        end
    end

    % Rescaled Block Iterative EM (RBI)
    if (options.param.RBI)
        if (options.param.verbose >= 3)
            disp("Computing RBI");
        end
        im = RBI(im, Sens, rhs, D);
    end

    % Dynamic RAMLA
    if (options.param.DRAMA)
        if (options.param.verbose >= 3)
            disp("Computing DRAMA");
        end
        im = DRAMA(im, Sens, rhs, options.param.lam_drama, iter, osa_iter, options.param.subsets);
    end

    % Complete data OSEM
    if (options.param.COSEM || options.param.ECOSEM)
        if (options.param.verbose >= 3)
            disp("Computing COSEM/ECOSEM");
        end
        % D is the sensitivity image divided by the number of subsets, the COSEM family needs the full sensitivity
        if (options.param.ECOSEM)
            [COSEMApu, im_vectors.C_co] = COSEM(im, rhs, im_vectors.C_co, options.param.subsets .* D, options.param.h, 2, osa_iter);
        else
            [im, im_vectors.C_co] = COSEM(im, rhs, im_vectors.C_co, options.param.subsets .* D, options.param.h, 2, osa_iter);
        end
    end

    % Enhanced COSEM
    if (options.param.ECOSEM)
        if (options.param.verbose >= 3)
            disp("Computing ECOSEM");
        end
        im = ECOSEM(im, D, OSEMApu, COSEMApu, options.param.epps);
    end

    % Accelerated COSEM
    if (options.param.ACOSEM)
        if (options.param.verbose >= 3)
            disp("Computing ACOSEM");
        end
        [im, im_vectors.C_co] = COSEM(im, rhs, im_vectors.C_co, options.param.subsets .* D, options.param.h, 1, osa_iter);
        ACOSEM_rhs = computeACOSEMWeight(options, im, uu, osa_iter, corrim, ii, randSum);
        im = im .* ACOSEM_rhs;
    end
    if (options.param.OSL_OSEM)
        if (options.param.verbose >= 3)
            disp("Computing OSL-OSEM");
        end
        % The denominator (sensitivity image + prior gradient) is limited to at least 10 % of the sensitivity image (and at
        % least epps where the sensitivity is ~0), so that the update cannot change sign or divide by zero with a strong
        % prior. Without a prior this is identical to OSEM
        if isempty(dUi)
            denom = Sens;
        else
            denom = max(Sens + dUi, max(0.1 * Sens, options.param.epps));
        end
        if (options.param.CT)
            im = EM(im, denom, options.param.flat * rhs);
        else
            im = EM(im, denom, rhs + options.param.epps);
        end
    elseif (options.param.BSREM)
        if (options.param.verbose >= 3)
            disp("Computing BSREM");
        end
        im = BSREM(im, rhs, options.param.lambda, iter, options.param.epps, options.param.enforcePositivity);
    elseif (options.param.MBSREM)
        if (options.param.verbose >= 3)
            disp("Computing MBSREM");
        end
        % MBSREM/SPS use an ascent-direction rhs (update is im + lambda * rhs), so the prior gradient is subtracted
        if isempty(dUi)
            [im, options.param] = MBSREM(im, rhs, options.param, options.param.U, options.param.lambda, iter, osa_iter, ii);
        else
            [im, options.param] = MBSREM(im, rhs - dUi, options.param, options.param.U, options.param.lambda, iter, osa_iter, ii);
        end
    elseif (options.param.ROSEM_MAP)
        if (options.param.verbose >= 3)
            disp("Computing ROSEM-MAP");
        end
        if (options.param.CT)
            im = ROSEM(im, Sens, options.param.flat * rhs, options.param.lambda, iter);
        else
            im = ROSEM(im, Sens, rhs, options.param.lambda, iter);
        end
    elseif (options.param.OSL_RBI)
        if (options.param.verbose >= 3)
            disp("Computing OSL-RBI");
        end
        if isempty(dUi)
            im = RBI(im, Sens, rhs, D);
        else
            im = RBI(im, Sens, rhs, D, options.param.beta, dUi, options.param.epps);
        end
    elseif (options.param.OSL_COSEM > 0)
        if (options.param.verbose >= 3)
            disp(['Computing OSL-COSEM ' num2str(options.param.OSL_COSEM)]);
        end
        % D is the sensitivity image divided by the number of subsets, the COSEM family needs the full sensitivity. The
        % denominator with the prior gradient is limited to at least 10 % of the sensitivity and at least epps
        if isempty(dUi)
            [im, im_vectors.C_co] = COSEM(im, rhs, im_vectors.C_co, options.param.subsets .* D, options.param.h, options.param.OSL_COSEM, osa_iter);
        else
            denom = max(options.param.subsets .* D + dUi, max(0.1 * options.param.subsets .* D, options.param.epps));
            [im, im_vectors.C_co] = COSEM(im, rhs, im_vectors.C_co, denom, options.param.h, options.param.OSL_COSEM, osa_iter);
        end
        if (options.param.OSL_COSEM == 1)
            ACOSEM_rhs = computeACOSEMWeight(options, im, uu, osa_iter, corrim, ii, randSum);
            im = im .* ACOSEM_rhs;
        end
    elseif (options.param.PKMA)
        if (options.param.verbose >= 3)
            disp("Computing PKMA");
        end
        if isempty(dUi)
            [im, options.param] = PKMA(im, rhs, options.param, iter, osa_iter, ii);
        else
            [im, options.param] = PKMA(im, rhs + dUi, options.param, iter, osa_iter, ii);
        end
    elseif (options.param.SPS)
        if (options.param.verbose >= 3)
            disp("Computing SPS");
        end
        if isempty(dUi)
            [im, options.param] = SPS(im, rhs, options.param.U, options.param.lambda, iter, osa_iter, ii, options.param);
        else
            [im, options.param] = SPS(im, rhs - dUi, options.param.U, options.param.lambda, iter, osa_iter, ii, options.param);
        end
    elseif (options.param.LSQR)
        if (options.param.verbose >= 3)
            disp("Computing LSQR");
        end
        [im_vectors,options.param] = LSQRim(options.param, iter, im_vectors, ii);
    elseif (options.param.CGLS)
        if (options.param.verbose >= 3)
            disp("Computing CGLS");
        end
        [im_vectors,options.param] = CGLS(options.param, iter, im_vectors, ii);
    elseif (options.param.PDHG || options.param.PDHGKL || options.param.PDHGL1 || options.param.PDDY)
        if (options.param.verbose >= 3)
            disp("Computing PDHG/PDHGKL/PDHGL1/PDDY");
        end
        if isempty(dUi)
            [im, options.param, im_vectors] = PDHG2(im, rhs, im_vectors, options.param, iter, osa_iter, ii);
        else
            [im, options.param, im_vectors] = PDHG2(im, rhs + dUi, im_vectors, options.param, iter, osa_iter, ii);
        end
    elseif (options.param.FISTA)
        if (options.param.verbose >= 3)
            disp("Computing FISTA");
        end
        % FISTA/FISTAL1 are non-MAP algorithms (checkAlgorithmSelection disables the priors), so no prior gradient is added
        [im_vectors, im, options.param] = FISTA(im, rhs, options.param, im_vectors, iter, osa_iter, ii);
    elseif (options.param.FISTAL1)
        if (options.param.verbose >= 3)
            disp("Computing FISTAL1");
        end
        [im_vectors, im, options.param] = FISTAL1(im, rhs, options.param, im_vectors, options.param.beta, iter, osa_iter, ii);
    end
    if ~options.param.LSQR && ~options.param.CGLS
        if iscell(im_vectors.recApu)
            im_vectors.recApu{ii} = im;
        else
            im_vectors.recApu = im;
        end
    end
end

% Enforce positivity after each subiteration when subsets are used (reconstructionAF.h). Without subsets the clamp is
% done after the iteration (after the BSREM/ROSEM-MAP step and the saving of the iterate)
if (options.param.enforcePositivity && options.param.subsets > 1 && ~options.param.CPType)
    if iscell(im_vectors.recApu)
        for ii = 1 : options.param.nMultiVolumes + 1
            im_vectors.recApu{ii}(im_vectors.recApu{ii} < options.param.epps) = options.param.epps;
        end
    else
        im_vectors.recApu(im_vectors.recApu < options.param.epps) = options.param.epps;
    end
end
end
