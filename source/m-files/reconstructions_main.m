function [pz,varargout] = reconstructions_main(inputStruct,varargin)
%% Main reconstruction file
% This function is used to compute various reconstructions with the selected method. Can be used with any sinogram or
% raw data.
%
% OUTPUT:
%   pz = The reconstructed image. Can be a 4D matrix if N number of
%   iterations are saved of if a dynamic reconstruction is performed.
%   image_properties = Reconstruction specific variables that are initially
%   stored in the options-struct (optional)
%   options = The modified input options-struct (optional)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2018-2024 Ville-Veikko Wettenhovi, Samuli Summala, Jonna Kangasniemi
%
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program. If not, see
% <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(inputStruct,'use_raw_data') && ~inputStruct.use_raw_data && isfield(inputStruct,'coincidences')
    inputStruct = rmfield(inputStruct, 'coincidences');
end

if nargin > 1 && ~isempty(varargin) && ~isempty(varargin{1})
    tyyppi = varargin{1};
else
    tyyppi = 0;
end

inputStruct.empty_weight = false;

inputStruct.listmode = false;
% tStart = 0;
% tStart_iter = 0;


if tyyppi == 0
    disp('Preparing for reconstruction...')
end

% Crate the class object
options = projectorClass(inputStruct);

if numel(options.param.partitions) > 1
    partitions = numel(options.param.partitions);
else
    partitions = options.param.partitions;
end

if options.param.listmode && partitions > 1
    error('Dynamic reconstruction is not yet supported with listmode/custom detector data')
end

if tyyppi < 2
    [options.param, RandProp, ScatterProp] = loadInputData(options.param);
end

if (options.param.CT && isfield(options.param,'flat') && options.param.flat == 0) || (options.param.CT && ~isfield(options.param,'flat'))
    warning('No flat value input! Using the maximum value as the flat value. Alternatively, input the flat value into options.flat')
    if options.param.useSingles
        options.param.flat = single(max(options.param.SinM(:)));
    else
        options.param.flat = double(max(options.param.SinM(:)));
    end
end


if options.param.L && options.param.MAP
    if ~isempty(options.param.a_L)
        if length(options.param.a_L(:)) < ((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1))
            error(['Weights vector options.param.a_L is too short, needs to be ' num2str(((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1))) ' in length'])
        elseif length(options.param.a_L(:)) > ((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1))
            error(['Weights vector options.param.a_L is too large, needs to be ' num2str(((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1))) ' in length'])
        end
    end
end
if (options.param.quad || options.param.FMH || options.param.L || options.param.weighted_mean || options.param.Huber || options.param.GGMRF) && options.param.MAP
    if ~isempty(options.param.weights)
        if length(options.param.weights(:)) < ((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1))
            error(['Weights vector is too short, needs to be ' num2str(((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1))) ' in length'])
        elseif length(options.param.weights(:)) > ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
            error(['Weights vector is too long, needs to be ' num2str(((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1))) ' in length'])
        end
        if ~isinf(options.param.weights(ceil(((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1)/2))))
            options.param.weights(ceil(((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1)/2))) = Inf;
        end
    else
        options.param.empty_weight = true;
    end
end
if options.param.Huber && options.param.MAP
    if ~isempty(options.param.weights_huber)
        if length(options.param.weights_huber(:)) < ((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1))
            error(['Huber weights vector is too short, needs to be ' num2str(((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1))) ' in length'])
        elseif length(options.param.weights_huber(:)) > ((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1))
            error(['Huber weights vector is too long, needs to be ' num2str(((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1))) ' in length'])
        end
        if ~isinf(options.param.weights_huber(ceil(((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1)/2))))
            options.param.weights_huber(ceil(((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1)/2))) = Inf;
        end
    end
end
if options.param.FMH && options.param.MAP
    if ~isempty(options.param.fmh_weights)
        if options.param.Nz == 1 || options.param.Ndz == 0
            if length(options.param.fmh_weights(:)) < (4*(options.param.Ndx*2+1))
                error(['Weights vector options.param.fmh_weights is too short, needs to be [' num2str(options.param.Ndx*2+1) ', 4] in size'])
            elseif length(options.param.fmh_weights(:)) > (4*(options.param.Ndx*2+1))
                error(['Weights vector options.param.fmh_weights is too large, needs to be [' num2str(options.param.Ndx*2+1) ', 4] in size'])
            end
        else
            if length(options.param.fmh_weights(:)) < (13*(options.param.Ndx*2+1))
                error(['Weights vector options.param.fmh_weights is too short, needs to be [' num2str(options.param.Ndx*2+1) ', 13] in size'])
            elseif length(options.param.fmh_weights(:)) > (13*(options.param.Ndx*2+1))
                error(['Weights vector options.param.fmh_weights is too long, needs to be [' num2str(options.param.Ndx*2+1) ', 13] in size'])
            end
        end
    end
end
if options.param.weighted_mean && options.param.MAP
    if ~isempty(options.param.weighted_weights)
        if length(options.param.weighted_weights(:)) < ((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1))
            error(['Weights vector options.param.weighted_weights is too short, needs to be ' num2str(((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1))) ' in length'])
        elseif length(options.param.weighted_weights(:)) > ((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1))
            error(['Weights vector options.param.weighted_weights is too long, needs to be ' num2str(((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1))) ' in length'])
        end
        if ~isinf(options.param.weighted_weights(ceil(((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1)/2))))
            options.param.weighted_weights(ceil(((options.param.Ndx*2+1) * (options.param.Ndy*2+1) * (options.param.Ndz*2+1)/2))) = Inf;
        end
    end
end

if partitions == 1 && options.param.nLayers == 1
    if ~options.param.CT && ~options.param.SPECT && numel(options.param.SinM) ~= options.param.Ndist * options.param.Nang * options.param.TotSinos * options.param.TOF_bins_used && options.param.listmode == 0 && ~options.param.attenuation_phase
        error('The number of elements in the input data does not match the input number of angles, radial distances and total number of sinograms multiplied together!')
    end
end

linAlg = recNames(5);

if ~options.param.usingLinearizedData
    for kk = 1 : numel(linAlg)
        if options.param.CT && options.param.(linAlg{kk}) && ~options.param.largeDim
            options.param = linearizeData(options.param);
            break
        end
    end
end

if ~options.param.listmode
    % Load correction data
    [options.param.normalization_correction, options.param.randoms_correction, options.param] = loadCorrections(options.param, RandProp, ScatterProp);
else
    % list-mode does not support corrections
    options.param.normalization_correction = false;
    options.param.randoms_correction = false;
    options.param.additionalCorrection = false;
    [options.param.normalization_correction, options.param.randoms_correction, options.param] = loadCorrections(options.param, RandProp, ScatterProp);
end

options.param = parseInputData(options.param, options.index);

% Remove negative values
if ~options.param.CT && (~options.param.LSQR && ~options.param.CGLS)
    if iscell(options.param.SinM)
        for llo = 1 : partitions
            options.param.SinM{llo}(options.param.SinM{llo} < 0) = 0;
        end
    else
        options.param.SinM(options.param.SinM < 0) = 0;
    end
end
if options.param.FDK
    options.param.precondTypeMeas(2) = true;
end
% Perform various prepass steps, if necessary
% These include computing weights, and loading anatomic reference images
[options.param] = prepass_phase(options.param, options.param.dz, options.param.dx, options.param.dy);

%%

% Compute the reconstructions
disp('Starting image reconstruction')

% Implementations 1, 4 and 5
if options.param.implementation ~= 2 && options.param.implementation ~= 3

    if options.param.projector_type > 3 && options.param.implementation == 4 && options.param.projector_type < 6
        error('Implementation 4 supports only projector types 1-3 and 6!')
    end

    if options.param.OSEM || options.param.ECOSEM || options.param.ROSEM || options.param.RBI || options.param.OSL_RBI || options.param.DRAMA || ...
            options.param.OSL_OSEM || options.param.ROSEM_MAP
        noSensIm = false;
    else
        noSensIm = true;
    end

    if (options.param.implementation == 1 || options.param.implementation == 4 || options.param.implementation == 5) || tyyppi == 1
        im_vectors = form_image_vectors(options.param, options.param.N, noSensIm);
    end

    if (options.param.RBI || options.param.OSL_RBI || options.param.COSEM || options.param.ACOSEM || options.param.OSL_COSEM > 0 || options.param.ECOSEM || (options.param.precondTypeImage(1) || options.param.precondTypeImage(2) || options.param.precondTypeImage(3)))
        computeD = true;
    else
        computeD = false;
    end

    if (options.param.precondTypeMeas(1))
        computeM = true;
    else
        computeM = false;
    end

    if (options.param.PDHG || options.param.PDHGKL || options.param.PDHGL1 || options.param.CV || options.param.PDDY)
        options.param.CPType = true;
    else
        options.param.CPType = false;
    end
    if options.param.storeFP
        fp = zeros(options.param.nRowsD, options.param.nColsD, options.param.nProjections, options.param.Niter, options.param.cType);
    end
    if options.param.PDAdaptiveType == 1
        options.param.alphaCP = ones(options.param.nMultiVolumes + 1, 1, options.param.cType);
    end

    % Loop through all time steps
    for llo = 1 : partitions
        if iscell(options.param.SinM)
            Sino = options.param.SinM{llo};
        else
            Sino = options.param.SinM;
        end
        %         koko = numel(Sino);
        %         TOFSize = koko / options.param.TOF_bins;
        if options.param.additionalCorrection
            if iscell(options.param.corrVector)
                if length(options.param.corrVector) > 1
                    corrVector = cast(options.param.corrVector{llo}, options.param.cType);
                else
                    corrVector = cast(options.param.corrVector{1}, options.param.cType);
                end
            else
                corrVector = cast(options.param.corrVector, options.param.cType);
            end
        else
            if options.param.implementation == 1 || (options.param.implementation == 4 && ~options.param.useSingles)
                corrVector = 0;
            else
                corrVector = single(0);
            end
        end
        if options.param.randoms_correction
            if iscell(options.param.SinDelayed)
                SinD = cast(options.param.SinDelayed{llo}, options.param.cType);
            else
                SinD = cast(options.param.SinDelayed, options.param.cType);
                options.param = rmfield(options.param, 'SinDelayed');
            end
            SinD = SinD(:);
        else
            SinD = cast([], options.param.cType);
        end

        Sino = Sino(:);

        if issparse(Sino)
            Sino = (full(Sino));
        end
        Sino = cast(Sino, options.param.cType);
        if options.param.normalization_correction
            options.param.normalization = cast(options.param.normalization, options.param.cType);
        end

        if llo > 1
            im_vectors = form_image_vectors(options.param, options.param.N, noSensIm, llo, im_vectors);
        end

        if (options.param.MBSREM || options.param.MRAMLA || options.param.SPS || options.param.precondTypeImage(7))
            if options.param.randoms_correction && ~options.param.CT && (options.param.MBSREM || options.param.MRAMLA)
                if iscell(options.param.SinDelayed)
                    allTrue = options.param.SinDelayed{llo} > 0;
                else
                    allTrue = options.param.SinDelayed > 0;
                end
            else
                allTrue = false;
            end
            if ((options.param.precondTypeImage(7) || options.param.SPS) || (~options.param.CT && ((options.param.randoms_correction && ~all(allTrue)) || options.param.randoms_correction == 0)))
                recApu = cell(options.param.nMultiVolumes + 1, 1);
                for ii = 1 : options.param.nMultiVolumes + 1
                    recApu{ii} = ones(options.param.N(ii), 1, options.param.cType);
                end
                E = ones(numel(Sino), 1, options.param.cType);
                for ll = 1 : options.param.subsets
                    if options.param.implementation == 1
                        A = formMatrix(options, ll);
                        E(options.index(options.nMeas(1) + 1 : options.nMeas(2))) = E(options.index(options.nMeas(1) + 1 : options.nMeas(2))) + A' * recApu{1};
                    else
                        E(options.index(options.nMeas(1) + 1 : options.nMeas(2))) = E(options.index(options.nMeas(1) + 1 : options.nMeas(2))) + forwardProject(options, recApu, ll, corrVector);
                    end
                end
                if (~options.param.CT && ((options.param.randoms_correction && ~all(allTrue)) || options.param.randoms_correction == 0) && (options.param.MBSREM || options.param.MRAMLA || options.param.SPS))
                    if (options.param.verbose >= 3)
                        disp("Computing epsilon value for MBSREM/MRAMLA");
                    end
                    options.param.epsilon_mramla = MBSREM_epsilon(Sino, options.param.epps, options.param.randoms_correction, SinD, E, options.param.TOF, options.param.TOF_bins, options.param.CT);
                    if (options.param.verbose >= 3)
                        disp("Epsilon value for MBSREM/MRAMLA computed");
                    end
                end
                if (options.param.precondTypeImage(7) || options.param.SPS)
                    if (options.param.verbose >= 3)
                        disp("Starting computation of dP for preconditioner type 6");
                    end
                    if (options.param.CT)
                        apu6 = Sino ./ options.param.flat;
                        if (options.param.randoms_correction)
                            EE = (options.param.flat .* apu6 - (Sino .* SinD .* apu6) ./ ((apu6 + SinD) * (apu6 + SinD)));
                        else
                            EE = (options.param.flat .* apu6);
                        end
                        E = E .* EE;
                    else
                        EE = zeros(size(Sino), options.param.cType);
                        EE(Sino > 0) = E(Sino > 0) ./ Sino(Sino > 0);
                        E = EE;
                    end
                    if options.param.implementation == 1 && options.param.nMultiVolumes == 0
                        options.param.dP = A * E;
                        if options.param.use_psf
                            options.param.dP = computeConvolution(options.param.dP, options.param, options.param.Nx(1), options.param.Ny(1), options.param.Nz(1), options.param.gaussK);
                        end
                    else
                        options.param.dP = backwardProject(options, E, 0, corrVector);
                    end
                    if iscell(options.param.dP)
                        for ii = 1 : options.param.nMultiVolumes + 1
                            options.param.dP{ii} = (options.param.subsets) ./ options.param.dP{ii};
                            options.param.dP{ii}(isnan(options.param.dP{ii})) = 1;
                            options.param.dP{ii}(isinf(options.param.dP{ii})) = 1;
                        end
                    else
                        options.param.dP = (options.param.subsets) ./ options.param.dP;
                        options.param.dP(isnan(options.param.dP)) = 1;
                        options.param.dP(isinf(options.param.dP)) = 1;
                    end
                    if (options.param.verbose >= 3)
                        disp("dP computation finished");
                    end
                end
                if options.param.implementation == 1
                    clear A
                end
            end
        end


        if llo == 1
            % Compute sensitivity image for the whole measurement domain
            if (computeD || options.param.compute_sensitivity_image)
                if (options.param.verbose >= 3)
                    disp("Starting computation of sensitivity image (D)");
                end
                options.param.D = ones(options.param.N(1), 1, options.param.cType);
                psf = options.param.use_psf;
                options.param.use_psf = false;
                for ll = 1 : options.param.subsets
                    oneInput = ones(options.nMeasSubset(ll), 1, options.param.cType);
                    if options.param.implementation == 1 && options.param.nMultiVolumes == 0
                        A = formMatrix(options, ll);
                        apu = A * oneInput;
                        options.param.D = options.param.D + apu;
                    else
                        options.param.D = options.param.D + backwardProject(options, oneInput, ll, corrVector);
                    end
                end
                options.param.use_psf = psf;
                if options.param.use_psf
                    options.param.D = computeConvolution(options.param.D, options.param, options.param.Nx(1), options.param.Ny(1), options.param.Nz(1), options.param.gaussK);
                end
                if iscell(options.param.D)
                    for ii = 1 : options.param.nMultiVolumes + 1
                        if (options.param.subsets > 1)
                            options.param.D{ii} = options.param.D{ii} ./ (options.param.subsets);
                        end
                        options.param.D{ii}(options.param.D{ii} == 0) = 1;
                    end
                else
                    if (options.param.subsets > 1)
                        options.param.D = options.param.D ./ (options.param.subsets);
                    end
                    options.param.D(options.param.D == 0) = 1;
                end
                if options.param.compute_sensitivity_image && ~computeD
                    im_vectors.Sens = repmat(options.param.D, 1, options.param.subsets);
                    options.param = rmfield(options.param, 'D');
                    noSensIm = 1;
                end
                if options.param.compute_sensitivity_image
                    options.param.compute_sensitivity_image = false;
                end
                if (options.param.verbose >= 3)
                    disp("Sensitivity image D computed");
                end
            end
            % Compute the measurement sensitivity image for each subset
            if (computeM)
                if (options.param.verbose >= 3)
                    disp("Starting computation of measurement sensitivity image (M)");
                end
                options.param.M = cell(options.param.subsets, 1);
                oneInput = cell(options.param.nMultiVolumes + 1, 1);
                for ii = 1 : options.param.nMultiVolumes + 1
                    oneInput{ii} = ones(options.param.N(ii), 1, options.param.cType);
                end
                for ll = 1 : options.param.subsets
                    if options.implementation == 1
                        if options.param.use_psf
                            oneInput{1} = computeConvolution(oneInput{1}, options.param, options.param.Nx(1), options.param.Ny(1), options.param.Nz(1), options.param.gaussK);
                        end
                        A = formMatrix(options, ll);
                        options.param.M{ll} = A' * oneInput{1};
                    else
                        options.param.M{ll} = forwardProject(options, oneInput, ll, corrVector);
                    end
                    options.param.M{ll}(options.param.M{ll} == 0) = 1;
                end
                if (options.param.verbose >= 3)
                    disp("M computation finished");
                end
            end
            % Use power method to compute the tau value, if necessary
            if ((options.param.CPType || options.param.FISTA || options.param.FISTAL1) && (isempty(options.param.tauCP) || options.param.tauCP(1) == 0))
                s = powerMethod(options);
                options.param.LCP = s;
                if options.param.precondTypeMeas(2)
                    options.param.tauCP = 1 ./ s;
                    apuVar = options.param.precondTypeMeas(2);
                    options.param.precondTypeMeas(2) = false;
                    s = powerMethod(options);
                    options.param.tauCPFilt = 1 ./ s;
                    options.param.precondTypeMeas(2) = apuVar;
                    options.param.LCP2 = s;
                else
                    options.param.tauCP = 1 ./ s;
                end
                options.param.sigmaCP = ones(size(options.param.tauCP), options.param.cType);
            end
            if (options.param.verbose >= 2 && (options.param.CPType || options.param.FISTA || options.param.FISTAL1))
                if options.param.precondTypeMeas(2)
                    disp(['Primal step size (tau) is ' num2str(options.param.tauCP(1)) ' with filtering']);
                    disp(['Primal step size (tau) is ' num2str(options.param.tauCPFilt(1)) ' without filtering']);
                else
                    disp(['Primal step size (tau) is ' num2str(options.param.tauCP(1))]);
                end
                disp(['Dual step size (sigma) is ' num2str(options.param.sigmaCP(1))]);
            end
        end

        %%
        for iter = 1 : options.param.Niter
            % if (options.param.save_iter || numel(options.param.saveNIter) > 0)
            % curIter = iter;

            for osa_iter = 1 : options.param.subsets
                if (~options.param.LSQR && ~options.param.CGLS) || iter == 1
                    [uu, rand] = splitMeas(options.param, Sino, options.nMeas, osa_iter, SinD);
                end
                [im_vectors, uu, options] = initializationStep(options, uu, im_vectors, options.nMeasSubset(osa_iter), iter, osa_iter);
                if options.param.verbose > 1
                    disp('Starting forward projection')
                end
                if options.param.implementation == 1
                    if options.param.subsets > 1 || iter == 1
                        if options.param.subsets > 1
                            A = formMatrix(options, osa_iter, 1, corrVector);
                        else
                            A = formMatrix(options, 0, 1, corrVector);
                        end
                    end
                    if iscell(im_vectors.recApu)
                        if options.param.use_psf
                            apu = computeConvolution(im_vectors.recApu{1}, options.param, options.param.Nx(1), options.param.Ny(1), options.param.Nz(1), options.param.gaussK);
                            outputFP = A' * apu;
                            clear apu
                        else
                            outputFP = A' * im_vectors.recApu{1};
                        end
                    else
                        if options.param.use_psf
                            apu = computeConvolution(im_vectors.recApu, options.param, options.param.Nx(1), options.param.Ny(1), options.param.Nz(1), options.param.gaussK);
                            outputFP = A' * apu;
                            clear apu
                        else
                            outputFP = A' * im_vectors.recApu;
                        end
                    end
                else
                    outputFP = forwardProject(options, im_vectors.recApu, osa_iter, corrVector);
                end
                if options.param.verbose > 1
                    disp('Forward projection complete')
                end
                % if options.param.listmode
                %    if iter == 1 && osa_iter == 1
                %       lor_a = lor_pixel_count_prepass(options.param, false);
                %       for lla = 1 : subsets
                %          summa(lla) = uint64(sum(lor_a(pituus(osa_iter)*koko+1:pituus(osa_iter + 1)*kerroin)));
                %       end
                %    end
                % end
                [outputFP, options.param, im_vectors, uu] = computeForwardStep(options.param, uu, outputFP, osa_iter, iter, options.nMeasSubset(osa_iter), 1, im_vectors, rand);
                if options.param.storeFP
                    if options.param.subset_type < 8 && options.param.subsets > 1
                        fp(options.index(options.nMeas(osa_iter) + 1 : options.nMeas(osa_iter + 1)) + options.param.totMeas * (iter - 1)) = outputFP;
                    elseif options.param.subsets == 1
                        fp(:,:,:,iter) = reshape(outputFP, options.param.nRowsD, options.param.nColsD, []);
                    else
                        fp(:, :, options.index(options.nMeas(osa_iter) + 1 : options.nMeas(osa_iter + 1)), iter) = reshape(outputFP, options.param.nRowsD, options.param.nColsD, []);
                    end
                end
                if (options.param.CT && (options.param.ACOSEM || options.param.OSL_COSEM > 0 || options.param.OSEM || options.param.COSEM || options.param.ECOSEM || ...
                        options.param.ROSEM || options.param.OSL_OSEM || options.param.ROSEM_MAP))
                    noSensIm = 1;
                    if (options.param.randoms_correction)
                        OSEMapu = (uu .* outputFP) ./ (outputFP + rand);
                    elseif (iter == 1  || (~options.param.saveSens && noSensIm == 0))
                        OSEMapu = uu;
                    end
                end
                if (options.param.CT && (options.param.ACOSEM || options.param.OSL_COSEM > 0 || options.param.OSEM || options.param.COSEM || options.param.ECOSEM || ...
                        options.param.ROSEM || options.param.OSL_OSEM || options.param.ROSEM_MAP))
                    if (options.param.randoms_correction || iter == 1 || (~options.param.saveSens && noSensIm == 0))

                        if options.param.implementation == 1 && options.param.nMultiVolumes == 0
                            if ~options.param.saveSens
                                im_vectors.Sens = A * OSEMapu;
                                if options.param.use_psf
                                    im_vectors.Sens = computeConvolution(im_vectors.Sens, options.param, options.param.Nx(1), options.param.Ny(1), options.param.Nz(1), options.param.gaussK);
                                end
                                im_vectors.Sens(im_vectors.Sens < options.param.epps, 1) = options.param.epps;
                            else
                                im_vectors.Sens(:,osa_iter) = A * OSEMapu;
                                if options.param.use_psf
                                    im_vectors.Sens(:,osa_iter) = computeConvolution(im_vectors.Sens(:,osa_iter), options.param, options.param.Nx(1), options.param.Ny(1), options.param.Nz(1), options.param.gaussK);
                                end
                                im_vectors.Sens(im_vectors.Sens(:,osa_iter) < options.param.epps, osa_iter) = options.param.epps;
                            end
                        else
                            im_vectors.rhs = backwardProject(options, OSEMapu, osa_iter, corrVector);
                            if ~options.param.saveSens || options.param.randoms_correction
                                im_vectors.Sens = im_vectors.rhs;
                                if iscell(im_vectors.rhs)
                                    for ii = 1 : options.param.nMultiVolumes + 1
                                        im_vectors.Sens{ii}(im_vectors.Sens{ii}(:,1) < options.param.epps, 1) = options.param.epps;
                                    end
                                else
                                    im_vectors.Sens(im_vectors.Sens < options.param.epps, 1) = options.param.epps;
                                end
                            else
                                if iscell(im_vectors.rhs)
                                    for ii = 1 : options.param.nMultiVolumes + 1
                                        im_vectors.Sens{ii}(:,osa_iter) = im_vectors.rhs{ii};
                                        im_vectors.Sens{ii}(im_vectors.Sens{ii}(:,osa_iter) < options.param.epps, osa_iter) = options.param.epps;
                                    end
                                else
                                    im_vectors.Sens(:,osa_iter) = im_vectors.rhs;
                                    im_vectors.Sens(im_vectors.Sens(:,osa_iter) < options.param.epps, osa_iter) = options.param.epps;
                                end
                            end
                        end
                    end
                end
                if noSensIm
                    if options.param.implementation == 1 && options.param.nMultiVolumes == 0
                        im_vectors.rhs = A * outputFP;
                        if options.param.use_psf
                            im_vectors.rhs = computeConvolution(im_vectors.rhs, options.param, options.param.Nx(1), options.param.Ny(1), options.param.Nz(1), options.param.gaussK);
                        end
                    else
                        im_vectors.rhs = backwardProject(options, outputFP, osa_iter, corrVector);
                    end
                else
                    if options.param.implementation == 1 && options.param.nMultiVolumes == 0
                        im_vectors.rhs = A * outputFP;
                        if options.param.use_psf
                            im_vectors.rhs = computeConvolution(im_vectors.rhs, options.param, options.param.Nx(1), options.param.Ny(1), options.param.Nz(1), options.param.gaussK);
                        end
                        if ~options.param.saveSens
                            im_vectors.Sens(:,1) = sum(A, 2);
                            if options.param.use_psf
                                im_vectors.Sens(:,1) = computeConvolution(im_vectors.Sens(:,1), options.param, options.param.Nx(1), options.param.Ny(1), options.param.Nz(1), options.param.gaussK);
                            end
                            im_vectors.Sens(im_vectors.Sens(:,1) < options.param.epps, 1) = options.param.epps;
                        else
                            im_vectors.Sens(:,osa_iter) = sum(A, 2);
                            if options.param.use_psf
                                im_vectors.Sens(:,osa_iter) = computeConvolution(im_vectors.Sens(:,osa_iter), options.param, options.param.Nx(1), options.param.Ny(1), options.param.Nz(1), options.param.gaussK);
                            end
                            im_vectors.Sens(im_vectors.Sens(:,osa_iter) < options.param.epps, osa_iter) = options.param.epps;
                        end
                    else
                        if ~options.param.saveSens
                            [im_vectors.rhs, im_vectors.Sens] = backwardProject(options, outputFP, osa_iter, corrVector);
                            if iscell(im_vectors.Sens)
                                for ii = 1 : options.param.nMultiVolumes + 1
                                    im_vectors.Sens{ii}(im_vectors.Sens{ii} < options.param.epps, 1) = options.param.epps;
                                end
                            else
                                im_vectors.Sens(im_vectors.Sens < options.param.epps, 1) = options.param.epps;
                            end
                        else
                            [im_vectors.rhs, Summ] = backwardProject(options, outputFP, osa_iter, corrVector);
                            if iscell(Summ)
                                for ii = 1 : options.param.nMultiVolumes + 1
                                    Summ{ii}(Summ{ii} < options.param.epps, 1) = options.param.epps;
                                    im_vectors.Sens{ii}(:, osa_iter) = Summ{ii};
                                end
                            else
                                Summ(Summ < options.param.epps, 1) = options.param.epps;
                                im_vectors.Sens(:, osa_iter) = Summ;
                            end
                        end
                    end
                end
                if options.param.subsets > 1
                    clear A
                end
                [im_vectors, options] = computeOSEstimates(im_vectors, options, iter, osa_iter, uu, corrVector);
                if options.param.verbose > 0 && options.param.subsets > 1
                    disp(['Sub-iteration ' num2str(osa_iter) ' completed'])
                end
            end
            im_vectors = init_next_iter(im_vectors, options.param, iter, llo);
            if options.param.use_psf && options.param.deblurring
                im_vectors = computeDeblur(im_vectors, options.param, iter, options.param.gaussK, options.param.Nx(1), options.param.Ny(1), options.param.Nz(1), llo);
            end
        end
        if partitions > 1 && options.param.verbose > 0
            disp(['Reconstructions for timestep ' num2str(llo) ' completed'])
        end
        if options.param.saveSens
            noSensIm = true;
        end
    end
    im_vectors = reshape_vectors(im_vectors, options.param);
    pz = images_to_cell(im_vectors);
    %% Implementation 2
    % OpenCL matrix free
    % Uses ArrayFire libraries
    % Only C++ code (no pure MATLAB implementation)
    % Supports all features as implementation 1 except for NLM
elseif options.param.implementation == 2 || options.param.implementation == 3
    %%

    options.param = double_to_single(options.param);
    % if partitions == 1
    % end
    if options.param.additionalCorrection
        if iscell(options.param.corrVector) && ~isa(options.param.corrVector{1}, 'single')
            for kk = 1 : length(options.param.corrVector)
                options.param.corrVector{kk} = single(full(options.param.corrVector{kk}));
            end
            options.param.corrVector = cell2mat(options.param.corrVector);
        elseif ~isa(options.param.corrVector, 'single')
            options.param.corrVector = single(options.param.corrVector);
        end
    else
        options.param.corrVector = single(0);
    end
    if ~options.param.largeDim
        if iscell(options.param.SinM) && issparse(options.param.SinM{1})
            for kk = 1 : length(options.param.SinM)
                options.param.SinM{kk} = single(full(options.param.SinM{kk}));
            end
        elseif issparse(options.param.SinM)
            options.param.SinM = single(full(options.param.SinM));
        end
        if iscell(options.param.SinM)
            options.param.SinM = cell2mat(options.param.SinM);
        end
        if ~isa(options.param.SinM, 'single')
            options.param.SinM = single(options.param.SinM);
        end
        if (options.param.randoms_correction || options.param.scatter_correction) && options.param.corrections_during_reconstruction
            options.param.randSize = uint64(diff(options.nMeas));
            if iscell(options.param.SinDelayed) && issparse(options.param.SinDelayed{1})
                for kk = 1 : length(options.param.SinDelayed)
                    options.param.SinDelayed{kk} = single(full(options.param.SinDelayed{kk}));
                end
            end
            if iscell(options.param.SinDelayed)
                options.param.SinDelayed = single(full(options.param.SinDelayed));
            end
            if ~isa(options.param.SinDelayed, 'single')
                options.param.SinDelayed = single(options.param.SinDelayed);
            end
        else
            options.param.randSize = uint64(1);
            options.param.SinDelayed = single(0);
        end
    end
    if ~isfield(options.param, 'normalization')
        options.param.normalization = single(0);
    end
    options.param.scatter = options.param.additionalCorrection;
    [pz, fp, res] = computeImplementation23(options.param, options.x, options.z, options.nMeas);
else
    error('Unsupported reconstruction method.');
end


if options.param.useEFOV && ~options.param.useMultiResolutionVolumes
    if options.param.transaxialEFOV
        nTrans = (options.param.Nx - options.param.NxOrig)/2;
        pz = pz(1 + nTrans:end-nTrans,1 + nTrans:end-nTrans,:,:,:);
    end
    if options.param.axialEFOV
        nAxial = (options.param.Nz - options.param.NzOrig)/2;
        pz = pz(:,:,1 + nAxial:end-nAxial,:,:,:);
    end
end
if options.param.storeFP && (options.param.subsets == 1 || options.param.subset_type >= 8) && options.param.implementation == 2 && ~options.param.FDK
    uu = 1;
    for kk = 1 : options.param.Niter
        for ii = 1 : options.param.subsets
            meas = options.nMeas(ii + 1) - options.nMeas(ii);
            fp{uu} = reshape(fp{uu}, options.param.nRowsD, options.param.nColsD, meas * options.param.TOF_bins);
            uu = uu + 1;
        end
    end
end
% Save various image properties, e.g. matrix size, sinogram dimensions, FOV
% size, regularization parameters, etc.
if nargout > 1
    varargout{1} = save_image_properties(options.param);
end
% if nargin > 2 && ~isempty(varargin{2}) && varargin{2}
%     apu = cell(2,1);
%     apu{1} = pz;
%     apu{2} = save_image_properties(options.param);
%     pz = apu;
% end
if nargout > 2
    varargout{2} = options;
end
if nargout > 3 && options.param.storeFP
    varargout{3} = fp;
elseif nargout > 3
    varargout{3} = [];
end
if nargout > 4 && options.param.storeResidual
    varargout{4} = res;
elseif nargout > 3
    varargout{4} = [];
end

end
