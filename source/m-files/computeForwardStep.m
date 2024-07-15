function [input, options, im_vector, y] = computeForwardStep(options, y, input, subIter, iter, length, ii, im_vector, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
kk = subIter + (iter - 1) * options.subsets;
if options.randoms_correction && ~isempty(varargin)
    rand = varargin{1};
end
if (options.precondTypeMeas(2) && options.filteringIterations > 0 && kk == options.filteringIterations + 1)
    if (options.verbose >= 3)
        disp("Filter iterations complete. Switching tau/sigma-values");
    end
    if (options.CPType || options.FISTA || options.FISTAL1 || options.ProxTGV || options.ProxTV)
        for ii = 1 : options.nMultiVolumes + 1
            if (options.sigmaCP(ii) == 1)
                options.tauCP(ii) = options.tauCPFilt(ii);
            elseif (options.sigmaCP(ii) == options.tauCP(ii))
                options.tauCP(ii) = options.tauCPFilt(ii);
                options.sigmaCP(ii) = options.tauCPFilt(ii);
            else
                options.sigmaCP(ii) = options.tauCPFilt(ii);
            end
            if (options.CPType)
                options.sigma2CP(ii) = options.sigmaCP(ii);
            end
            if (options.PDAdaptiveType == 1)
                options.alphaCP(ii) = 1;
            end
        end
        options.LCP = options.LCP2;
    end
    if (options.MRAMLA || options.MBSREM || options.SPS || options.RAMLA || options.BSREM || options.ROSEM || options.ROSEM_MAP || options.PKMA)
        options.lambda = options.lambdaFiltered;
    end
    options.precondTypeMeas(2) = false;
end
if (options.randoms_correction)
    if ((options.MBSREM || options.MRAMLA || options.SPS) && ~options.CT && ~all(rand > 0))
        indeksit = y > 0 && rand == 0 && input <= options.epsilon_mramla;
        indS = any(indeksit);
    end
    if (~options.CT)
        if (options.verbose >= 3)
            disp("Adding randoms/scatter data to forward projection");
        end
        if (options.TOF)
            for to = 1 : options.TOF_bins
                input(numel(rand) * (to - 1) + 1 : numel(rand) * to) = input(numel(rand) * (to - 1) + 1 : numel(rand) * to) + rand;
            end
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
        % if (options.randoms_correction)
        %     input = y ./ (input + rand) - 1;
        % else
            input = y ./ (input) - 1;
            input(isinf(input)) = 0;
        % end
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
        % if (options.randoms_correction)
        %     input = 1 - y ./ (input + rand);
        % else
            input = 1 - y ./ (input);
        % end
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
        if (options.randoms_correction && indS)
            input(indeksit) = y(indeksit) ./ options.epsilon_mramla - 1 - (y(indeksit) ./ (options.epsilon_mramla .* options.epsilon_mramla)) .* (input(indeksit) - options.epsilon_mramla);
            input(~indeksit) = y(~indeksit) ./ (input(~indeksit)) - 1;
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
            input = (im_vector.pCP{subIter} + options.sigmaCP(ii) .* res) ./ (1 + options.sigmaCP(ii));
        % end
    end
    im_vector.pCP{subIter} = input;
elseif (options.PDHGKL)
    if (options.verbose >= 3)
        disp("Computing PDHG/CV with KL");
    end
    % if options.randoms_correction
    %     input = input + rand;
    % end
    if (options.precondTypeMeas(1) || options.precondTypeMeas(2))
        if (options.precondTypeMeas(2))
            apu1 = y;
            apu1 = applyMeasPreconditioning(options, apu1, length, subIter);
            input = applyMeasPreconditioning(options, input, length, subIter);
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
input(isnan(input)) = options.epps;
input(isinf(input)) = options.epps;
end