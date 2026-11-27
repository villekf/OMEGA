function [input, options, im_vector, y] = computeForwardStep(options, y, input, subIter, iter, length, ii, im_vector, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
kk = subIter + (iter - 1) * options.subsets;
if options.randoms_correction && ~isempty(varargin)
    rand = varargin{1};
end
% TOF-resolved additive term (randomsTOF) has a value for each TOF bin, otherwise the same values are used for each TOF bin
% (list-mode data has one value per event, in the TOF bin of the event, for both the data and the additive term)
randomsTOF = isfield(options, 'randomsTOF') && options.randomsTOF;
tileRandoms = options.TOF && ~randomsTOF && options.listmode == 0;
% LSQR and CGLS are linear least-squares methods and need a linear operator: the epps only guards the divisions of the
% Poisson methods. The epps is added before the randoms/scatter (and before the modified log-likelihood indices below)
if ~options.CT && ~(options.LSQR || options.CGLS)
    input = input + options.epps;
end
if (options.precondTypeMeas(2) && options.filteringIterations > 0 && kk == options.filteringIterations + 1)
    if (options.verbose >= 3)
        disp("Filter iterations complete. Switching tau/sigma-values");
    end
    if (options.CPType || options.FISTA || options.FISTAL1 || options.TGV || options.ProxTV)
        % The loop variable must not overwrite the volume index ii of this function
        for vv = 1 : options.nMultiVolumes + 1
            if (options.sigmaCP(vv) == 1)
                options.tauCP(vv) = options.tauCPFilt(vv);
            elseif (options.sigmaCP(vv) == options.tauCP(vv))
                options.tauCP(vv) = options.tauCPFilt(vv);
                options.sigmaCP(vv) = options.tauCPFilt(vv);
            else
                options.sigmaCP(vv) = options.tauCPFilt(vv);
            end
            if (options.CPType)
                options.sigma2CP(vv) = options.sigmaCP(vv);
            end
            if (options.PDAdaptiveType == 1)
                options.alphaCP(vv) = 1;
            end
        end
        % LCP2 exists only if the step sizes were computed with the power method
        if isfield(options, 'LCP2')
            options.LCP = options.LCP2;
        end
    end
    % options.lambda holds the filter-phase values and options.lambdaFiltered the values after the filtering (swapped in
    % prepass_phase), as in the native code
    if (options.MRAMLA || options.MBSREM || options.SPS || options.RAMLA || options.BSREM || options.ROSEM || options.ROSEM_MAP || options.PKMA)
        options.lambda = options.lambdaFiltered;
    end
    options.precondTypeMeas(2) = false;
end
% Ahn & Fessler's modified log-likelihood applies to all LORs with zero additive term (also without randoms/scatter
% correction). The indices must be computed before the randoms/scatter are added to the forward projection
indS = false;
if ((options.MBSREM || options.MRAMLA || options.SPS) && ~options.CT)
    if ~options.randoms_correction
        indeksit = y > 0 & input <= options.epsilon_mramla;
        indS = any(indeksit);
    elseif ~all(rand > 0)
        if tileRandoms
            indeksit = y > 0 & repmat(rand == 0, options.TOF_bins, 1) & input <= options.epsilon_mramla;
        else
            indeksit = y > 0 & rand == 0 & input <= options.epsilon_mramla;
        end
        indS = any(indeksit);
    end
end
if (options.randoms_correction)
    if (~options.CT)
        if (options.verbose >= 3)
            disp("Adding randoms/scatter data to forward projection");
        end
        if tileRandoms
            input = input + repmat(rand, options.TOF_bins, 1);
        else
            input = input + rand;
        end
    end
end
if options.CPType && options.subsets > 1
    im_vector.p0CP = im_vector.pCP{subIter};
end
if (options.ACOSEM || options.OSL_COSEM > 0 || options.OSEM || options.COSEM || options.ECOSEM || options.ROSEM || options.OSL_OSEM || options.ROSEM_MAP)
    if (options.verbose >= 3)
        disp("Computing (OSL-)(A)(E)(C)(R)OSEM");
    end
    if (options.CT)
        input(input < 0) = 0;
        if (options.verbose >= 3)
            disp("CT mode");
        end
        if (options.verbose >= 3 && options.randoms_correction)
            disp("Adding scatter data to forward projection");
        end
        if (options.OSEM || options.ROSEM || options.OSL_OSEM || options.ROSEM_MAP)
            input = exp(-input);
        else
            input = exp(-input) ./ y;
        end
    else
        if (options.verbose >= 3)
            disp("PET mode");
        end
        % if (options.randoms_correction)
        %     input = y ./ (input + rand);
        % else
        % input(input == 0) = 1;
            input = y ./ (input);
        % end
    end
elseif (options.RAMLA || options.BSREM || options.RBI || options.OSL_RBI || options.DRAMA)
    if (options.verbose >= 3)
        disp("Computing RAMLA/BSREM/RBI/DRAMA");
    end
    if (options.CT)
        input(input < 0) = 0;
        if (options.verbose >= 3)
            disp("CT mode");
        end
        if (options.verbose >= 3 && options.randoms_correction)
            disp("Adding scatter data to forward projection");
        end
        if (options.randoms_correction)
            input = exp(-input) .* options.flat - (y .* exp(-input)) ./ (exp(-input) + rand);
        else
            input = exp(-input) .* options.flat - y;
        end
    else
        if (options.verbose >= 3)
            disp("PET mode");
        end
        % List-mode data: the -1 term is replaced by the sensitivity image D (subtracted after the backprojection)
        if options.listmode > 0
            input = y ./ (input);
        else
            input = y ./ (input) - 1;
        end
    end
elseif (options.PKMA)
    if (options.verbose >= 3)
        disp("Computing PKMA");
    end
    if (options.CT)
        input(input < 0) = 0;
        if (options.verbose >= 3)
            disp("CT mode");
        end
        if (options.verbose >= 3 && options.randoms_correction)
            disp("Adding scatter data to forward projection");
        end
        if (options.randoms_correction)
            input = (y .* exp(-input)) ./ (exp(-input) + rand) - exp(-input) .* options.flat;
        else
            input = y - exp(-input) .* options.flat;
        end
    else
        if (options.verbose >= 3)
            disp("PET mode");
        end
        % List-mode data: the 1 term is replaced by the sensitivity image D (the backprojection is subtracted from D)
        if options.listmode > 0
            input = y ./ (input);
        else
            input = 1 - y ./ (input);
        end
    end
    input = applyMeasPreconditioning(options, input, length, subIter);
elseif (options.MBSREM || options.MRAMLA || options.SPS)
    if (options.verbose >= 3)
        disp("Computing MBSREM/MRAMLA/SPS");
    end
    if (options.CT)
        input(input < 0) = 0;
        if (options.verbose >= 3)
            disp("CT mode");
        end
        if (options.verbose >= 3 && options.randoms_correction)
            disp("Adding scatter data to forward projection");
        end
        if (options.randoms_correction)
            input = exp(-input) .* options.flat - (y .* exp(-input)) ./ (exp(-input) + rand);
        else
            input = exp(-input) .* options.flat - y;
        end
    else
        if (options.verbose >= 3)
            disp("PET mode");
        end
        if (indS)
            input(indeksit) = y(indeksit) ./ options.epsilon_mramla - 1 - (y(indeksit) ./ (options.epsilon_mramla .* options.epsilon_mramla)) .* (input(indeksit) - options.epsilon_mramla);
            input(~indeksit) = y(~indeksit) ./ (input(~indeksit) + options.epps) - 1;
            % The formulas above include the -1 term. For list-mode data it is replaced by the sensitivity image D
            % (subtracted after the backprojection), so remove it here
            if options.listmode > 0
                input = input + 1;
            end
        elseif options.listmode > 0
            input = y ./ (input);
        else
            input = y ./ (input) - 1;
        end
    end
    input = applyMeasPreconditioning(options, input, length, subIter);
elseif (options.LSQR)
    if (options.verbose >= 3)
        disp("Computing LSQR");
    end
    input = input - options.alphaLSQR .* y;
    options.betaLSQR = norm(input);
    input = input ./ options.betaLSQR;
    y = input;
elseif (options.CGLS)
    if (options.verbose >= 3)
        disp("Computing CGLS");
    end
    normi = norm(input);
    options.alphaCGLS = options.gammaCGLS ./ (normi .* normi);
    input = im_vector.rCGLS - options.alphaCGLS .* input;
    im_vector.rCGLS = input;
elseif (options.PDHG || options.PDDY)
    if (options.verbose >= 3)
        disp("Computing PDHG");
    end
    res = input - y;
    res = applyMeasPreconditioning(options, res, length, subIter);
    if (options.precondTypeMeas(2) && kk <= options.filteringIterations)
        if (options.verbose >= 3)
            disp("Computing inverse with circulant matrix");
        end
        input = (im_vector.pCP{subIter} + options.sigmaCP(ii) .* res);
        if options.subset_type == 4
            input = reshape(input, options.nRowsD, length / options.nRowsD);
        else
            input = reshape(input, options.nRowsD, options.nColsD, length / (options.nRowsD * options.nColsD));
        end
        if (options.PDAdaptiveType == 1)
            options.Ffilter = ifft(options.filter) .* options.sigmaCP(ii);
            options.Ffilter(1) = options.Ffilter(1) + 1;
            options.Ffilter = real(fft(options.Ffilter));
        end
        input = filteringInverse(options.Ffilter, input, options.Nf);
    else
        % if (options.ProxTGV)
        %     if (options.verbose >= 3)
        %         disp("Computing Proximal TGV");
        %     end
        %     input = (im_vector.pCP{subIter} + options.sigmaCP(ii) .* res) ./ (1 + options.sigmaCP(ii) .* options.betaReg);
        % else
            if (options.verbose >= 3)
                disp("Computing PDHG");
            end
            % Diagonal preconditioning (Pock and Chambolle 2011) uses Sigma = sigma / M. res has already been divided by M,
            % so the numerator contains Sigma * (Ax - y). The prox of the conjugate of ||z - y||^2 / 2 is
            % (v - Sigma y) / (1 + Sigma) elementwise, i.e. the denominator must also use Sigma
            if (options.precondTypeMeas(1))
                input = (im_vector.pCP{subIter} + options.sigmaCP(ii) .* res) ./ (1 + options.sigmaCP(ii) ./ options.M{subIter});
            else
                input = (im_vector.pCP{subIter} + options.sigmaCP(ii) .* res) ./ (1 + options.sigmaCP(ii));
            end
        % end
    end
    im_vector.pCP{subIter} = input;
elseif (options.PDHGKL)
    if (options.verbose >= 3)
        disp("Computing PDHG with KL");
    end
    % if options.randoms_correction
    %     input = input + rand;
    % end
    if (options.precondTypeMeas(1) || options.precondTypeMeas(2))
        if (options.precondTypeMeas(2))
            apu1 = y;
            apu1 = applyMeasPreconditioning(options, apu1, length, subIter);
            input = applyMeasPreconditioning(options, input, length, subIter);
            % The filtered data and forward projection can be negative, which the KL prox does not allow
            apu1(apu1 < 0) = 0;
            input(input < 0) = 0;
            input = .5 * (1 + im_vector.pCP{subIter} + options.sigmaCP(ii) .* input - sqrt((im_vector.pCP{subIter} + options.sigmaCP(ii) .* input - 1).^2 + 4 .* options.sigmaCP(ii) .* apu1));
        else
            if (options.verbose >= 3)
                disp("Applying diagonal normalization preconditioner (1 / (A1)), type 0");
            end
            input = .5 * (1 + im_vector.pCP{subIter} + options.sigmaCP(ii) .* input ./ options.M{subIter} - sqrt((im_vector.pCP{subIter} + options.sigmaCP(ii) .* input ./ options.M{subIter} - 1).^2 + 4 * options.sigmaCP(ii) .* y ./ options.M{subIter}));
        end
    else
        input = .5 * (1 + im_vector.pCP{subIter} + options.sigmaCP(ii) .* input - sqrt((im_vector.pCP{subIter} + options.sigmaCP(ii) .* input - 1).^2 + 4 .* options.sigmaCP(ii) .* y));
    end
    im_vector.pCP{subIter} = input;
elseif (options.PDHGL1)
    if (options.verbose >= 3)
        disp("Computing CPL1/TVL1/TGVL1");
    end
    res = input - y;
    res = applyMeasPreconditioning(options, res, length, subIter);
    input = (im_vector.pCP{subIter} + options.sigmaCP(ii) .* res);
    input = input ./ max(1, abs(input));
    im_vector.pCP{subIter} = input;
end
if (options.CPType && options.subsets > 1)
    if (options.verbose >= 3)
        disp("Computing PDHG with subsets");
    end
    input = input - im_vector.p0CP;
    % im_vector.p0CP{subIter} = im_vector.pCP{subIter};
end
if (options.FISTA || options.FISTAL1)
    if (options.verbose >= 3)
        disp("Computing FISTA/L1");
    end
    input = input - y;
    input = applyMeasPreconditioning(options, input, length, subIter);
end
% As in the native code, only infinite values are replaced (NaN values are kept)
input(isinf(input)) = options.epps;
end