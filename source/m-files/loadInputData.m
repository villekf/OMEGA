function [options, varargout] = loadInputData(options)
%LOADINPUTDATA Loads the correct previously saved input data
%   Loads all the necessary input data. This includes measurement data,
%   potentially randoms and scatter. Also precorrects data if it has not
%   been precorrected and precorrection has been selected.
% Load the measurement data if it does not exist in options.SinM
% Raw data
if options.listmode
    if nargout > 1
        varargout{1} = [];
        varargout{2} = [];
    end
    return
end
if numel(options.partitions) > 1
    partitions = numel(options.partitions);
elseif isempty(options.partitions)
    partitions = 1;
else
    partitions = options.partitions;
end
if options.use_raw_data
    RandProp.smoothing = false;
    RandProp.variance_reduction = false;
    ScatterProp.smoothing = false;
    ScatterProp.variance_reduction = false;
    ScatterProp.normalization = false;
    % Static case
    if partitions == 1
        load_string = [options.machine_name '_measurements_' options.name '_static_raw'];
        if options.use_ASCII && options.use_machine == 0
            load_string =  [load_string '_ASCII.mat'];
        elseif options.use_LMF && options.use_machine == 0
            load_string =  [load_string '_LMF.mat'];
        elseif options.use_root && options.use_machine == 0
            load_string =  [load_string '_root.mat'];
        else
            load_string =  [load_string '_listmode.mat'];
        end
    else
        % Dynamic case
        load_string = [options.machine_name '_measurements_' options.name '_' num2str(partitions) 'timepoints_for_total_of_' ...
            num2str(options.tot_time) 's_raw'];
        load_string2 = [options.machine_name '_measurements_' options.name '_' num2str(partitions) 'timepoints_for_total_of_ ' ...
            num2str(options.tot_time) 's_raw'];
        if options.use_ASCII && options.use_machine == 0
            if exist([load_string '_ASCII.mat'], 'file') == 0
                load_string = load_string2;
            end
            load_string =  [load_string '_ASCII.mat'];
        elseif options.use_LMF && options.use_machine == 0
            if exist([load_string '_LMF.mat'], 'file') == 0
                load_string = load_string2;
            end
            load_string =  [load_string '_LMF.mat'];
        elseif options.use_root && options.use_machine == 0
            if exist([load_string '_root.mat'], 'file') == 0
                load_string = load_string2;
            end
            load_string =  [load_string '_root.mat'];
        else
            if exist([load_string '_listmode.mat'], 'file') == 0
                load_string = load_string2;
            end
            load_string =  [load_string '_listmode.mat'];
        end
    end
    % Load the actual data and put it into options.SinM
    if options.reconstruct_trues == false && ~options.reconstruct_scatter && (isfield(options, 'coincidences') == 0 ...
            || numel(options.coincidences) == 0) && options.use_machine < 2
        options.SinM = loadStructFromFile(load_string, 'coincidences');
    elseif ~options.reconstruct_trues && ~options.reconstruct_scatter && isfield(options, 'coincidences')
        options.SinM = options.coincidences;
    elseif options.reconstruct_trues
        % Load Trues
        options.SinM = loadStructFromFile(load_string, 'true_coincidences');
    elseif options.reconstruct_scatter
        % Load scattered coincidences
        options.SinM = loadStructFromFile(load_string, 'scattered_coincidences');
    end
    % Perform corrections if needed
    if options.randoms_correction && ~options.reconstruct_trues && ~options.reconstruct_scatter
        if ((options.use_ASCII || options.use_LMF || options.use_root) && options.use_machine == 0) || options.use_machine == 1
            options.SinDelayed = loadStructFromFile(load_string, 'delayed_coincidences');
            if ~isfield(options,'SinDelayed')
                warning('No randoms correction data detected, disabling randoms correction!')
                options.randoms_correction = false;
            elseif iscell(options.SinDelayed)
                if numel(options.SinDelayed{1}) == 1
                    warning('No randoms correction data detected, disabling randoms correction!')
                    options.randoms_correction = false;
                end
            else
                if numel(options.SinDelayed) == 1
                    warning('No randoms correction data detected, disabling randoms correction!')
                    options.randoms_correction = false;
                end
            end
        else
            if ~isfield(options,'SinDelayed') || numel(options.SinDelayed) == 0
                options = loadDelayedData(options);
            end
        end
    end
    if options.scatter_correction && ~options.corrections_during_reconstruction
        if ~isfield(options,'ScatterC') || numel(options.ScatterC) == 0
            options = loadScatterData(options);
        end
    end
    clear coincidences true_coincidences delayed_coincidences
    options = rmfield(options, 'coincidences');
    % Sinogram data
else
    RandProp.smoothing = false;
    RandProp.variance_reduction = false;
    ScatterProp.smoothing = false;
    ScatterProp.variance_reduction = false;
    ScatterProp.normalization = false;
    % Static case
    if partitions == 1
        % TOF
        if options.TOF
            load_string = [options.machine_name '_' options.name '_TOFsinograms_combined_static_' num2str(options.Ndist) ...
                'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span)];
            load_string_TOF = ['_' num2str(options.TOF_bins) 'bins_' num2str(options.TOF_width*1e12) 'psBinSize_' ...
                num2str(options.TOF_noise_FWHM*1e12) 'psFWHM'];
            load_string = [load_string load_string_TOF];
        else
            load_string = [options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) ...
                'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span)];
        end
        if options.use_machine == 0
            sinoFile = [load_string '.mat'];
        elseif options.use_machine == 1
            sinoFile = [load_string '_listmode.mat'];
        elseif options.use_machine == 2
            sinoFile = [load_string '_machine_sinogram.mat'];
        elseif options.use_machine == 3
            sinoFile = [load_string '_listmode_sinogram.mat'];
        end
    else
        % Dynamic case
        % TOF
        if options.TOF
            load_string = [options.machine_name '_' options.name '_TOFsinograms_combined_' num2str(poptions.artitions) ...
                'timepoints_for_total_of_' num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) ...
                'x' num2str(options.TotSinos) '_span' num2str(options.span)];
            load_string_TOF = ['_' num2str(options.TOF_bins) 'bins_' num2str(options.TOF_width*1e12) 'psBinSize_' ...
                num2str(options.TOF_noise_FWHM*1e12) 'psFWHM'];
            load_string = [load_string load_string_TOF];
            load_string2 = load_string;
        else
            load_string = [options.machine_name '_' options.name '_sinograms_combined_' num2str(partitions) 'timepoints_for_total_of_' ...
                num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                num2str(options.span)];
            load_string2 = [options.machine_name '_' options.name '_sinograms_combined_' num2str(partitions) 'timepoints_for_total_of_ ' ...
                num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                num2str(options.span)];
        end
        if options.use_machine == 0
            sinoFile = [load_string '.mat'];
            if exist(sinoFile, 'file') == 0
                sinoFile = [load_string2 '.mat'];
            end
        elseif options.use_machine == 1
            sinoFile = [load_string '_listmode.mat'];
            if exist(sinoFile, 'file') == 0
                sinoFile = [load_string2 '_listmode.mat'];
            end
        elseif options.use_machine == 2
            sinoFile = [load_string '_machine_sinogram.mat'];
            if exist(sinoFile, 'file') == 0
                sinoFile = [load_string2 '_machine_sinogram.mat'];
            end
        elseif options.use_machine == 3
            sinoFile = [load_string '_listmode_sinogram.mat'];
            if exist(sinoFile, 'file') == 0
                sinoFile = [load_string2 '_listmode_sinogram.mat'];
            end
        end
    end
    % Put the data into options.SinM
    if ~options.reconstruct_trues && ~options.reconstruct_scatter
        loadRaw = false;
        if isfield(options,'SinM') == 0 || numel(options.SinM) == 0
            if ((options.randoms_correction || options.scatter_correction || options.normalization_correction) && ~options.corrections_during_reconstruction) ...
                    || options.fill_sinogram_gaps
                try
                    [options.SinM, appliedCorrections] = loadStructFromFile(sinoFile, 'SinM','appliedCorrections');
                catch
                    [options.SinM, appliedCorrections] = loadStructFromFile(sinoFile, 'raw_SinM','appliedCorrections');
                    loadRaw = true;
                end
            else
                options.SinM = loadStructFromFile(sinoFile, 'raw_SinM');
                appliedCorrections = [];
            end
        else
            try
                appliedCorrections = loadStructFromFile(sinoFile,'appliedCorrections');
            catch
                appliedCorrections = [];
            end
        end
        % Apply corrections
        % Corrections not applied during reconstruction
        if ~options.corrections_during_reconstruction && ~isempty(appliedCorrections)
            normalization_correction = options.normalization_correction;
            % Normalization
            if appliedCorrections.normalization && ~options.normalization_correction
                warning('Normalization correction not selected, but data precorrected with normalization! Precorrecting without normalization.')
                [options.SinM, appliedCorrections] = loadStructFromFile(sinoFile, 'raw_SinM','appliedCorrections');
                loadRaw = true;
            end
            if appliedCorrections.normalization
                normalization_correction = false;
            end
            % Randoms
            randoms_correction = options.randoms_correction;
            % Check if data is precorrected with something that has not
            % been selected
            % Variance reduction
            if ~isempty(strfind(appliedCorrections.randoms,'variance reduction')) && ~options.variance_reduction && options.randoms_correction
                warning('Randoms variance correction not selected, but data precorrected with randoms with applied variance reduction! Precorrecting without variance reduction.')
                [options.SinM, appliedCorrections] = loadStructFromFile(sinoFile, 'raw_SinM','appliedCorrections');
                loadRaw = true;
                if appliedCorrections.normalization && ~normalization_correction
                    normalization_correction = true;
                end
                % Randoms smoothing
            elseif ~isempty(strfind(appliedCorrections.randoms,'smoothing')) && ~options.randoms_smoothing && options.randoms_correction
                warning('Randoms smoothing not selected, but data precorrected with randoms with applied smoothing! Precorrecting without smoothing.')
                [options.SinM, appliedCorrections] = loadStructFromFile(sinoFile, 'raw_SinM','appliedCorrections');
                loadRaw = true;
                if appliedCorrections.normalization && ~normalization_correction
                    normalization_correction = true;
                end
                % Simple randoms correction
            elseif isempty(strfind(appliedCorrections.randoms,'randoms correction')) && ~loadRaw && randoms_correction
                [options.SinM, appliedCorrections] = loadStructFromFile(sinoFile, 'raw_SinM','appliedCorrections');
                loadRaw = true;
                if appliedCorrections.normalization && ~normalization_correction
                    normalization_correction = true;
                end
                % If data has not been precorrected before, precorrect
                % it now
            elseif ~isempty(strfind(appliedCorrections.randoms,'randoms correction')) && randoms_correction && ...
                    ((options.variance_reduction && isempty(strfind(appliedCorrections.randoms,'variance reduction'))) ...
                    || (options.randoms_smoothing && isempty(strfind(appliedCorrections.randoms,'smoothing'))))
                warning('Randoms corrections selected, but data not precorrected with selected options! Precorrecting with selected options.')
                [options.SinM, appliedCorrections] = loadStructFromFile(sinoFile, 'raw_SinM','appliedCorrections');
                loadRaw = true;
                if appliedCorrections.normalization && ~normalization_correction
                    normalization_correction = true;
                end
            elseif ~isempty(strfind(appliedCorrections.randoms,'randoms correction')) && ~loadRaw
                randoms_correction = false;
            end
            % Same as above, but for scatter
            if ~isempty(strfind(appliedCorrections.scatter,'variance reduction')) && ~options.scatter_variance_reduction && options.scatter_correction
                warning('Scatter variance correction not selected, but data precorrected with scatter with applied variance reduction! Precorrecting without variance reduction.')
                [options.SinM, appliedCorrections] = loadStructFromFile(sinoFile, 'raw_SinM','appliedCorrections');
                if ~randoms_correction && options.randoms_correction
                    randoms_correction = true;
                end
                if appliedCorrections.normalization && ~normalization_correction
                    normalization_correction = true;
                end
                loadRaw = true;
            elseif ~isempty(strfind(appliedCorrections.scatter,'smoothing')) && ~options.scatter_smoothing && options.scatter_correction
                warning('Scatter smoothing not selected, but data precorrected with scatter with applied smoothing! Precorrecting without smoothing.')
                [options.SinM, appliedCorrections] = loadStructFromFile(sinoFile, 'raw_SinM','appliedCorrections');
                if ~randoms_correction && options.randoms_correction
                    randoms_correction = true;
                end
                if appliedCorrections.normalization && ~normalization_correction
                    normalization_correction = true;
                end
                loadRaw = true;
            elseif isempty(strfind(appliedCorrections.scatter,'scatter correction')) && ~loadRaw && options.scatter_correction
                [options.SinM, appliedCorrections] = loadStructFromFile(sinoFile, 'raw_SinM','appliedCorrections');
                if ~randoms_correction && options.randoms_correction
                    randoms_correction = true;
                end
                if appliedCorrections.normalization && ~normalization_correction
                    normalization_correction = true;
                end
                loadRaw = true;
            elseif ~isempty(strfind(appliedCorrections.scatter, 'scatter correction')) && ~loadRaw
                options.scatter_correction = false;
            end
            options.randoms_correction = randoms_correction;
            options.normalization_correction = normalization_correction;
            % Sinogram gap filling
            if (~appliedCorrections.gapFilling || loadRaw) && options.fill_sinogram_gaps
                appliedCorrections.gapFilling = true;
                options = applyGapFilling(options);
            end
            % Corrections are applied during reconstruction
        elseif options.corrections_during_reconstruction && ~isempty(appliedCorrections)
            % Gap filling
            if (appliedCorrections.normalization || ~isempty(appliedCorrections.randoms) || ~isempty(appliedCorrections.scatter) || ~appliedCorrections.gapFilling) ...
                    && options.fill_sinogram_gaps && options.det_per_ring < options.det_w_pseudo
                options.SinM = loadStructFromFile(sinoFile, 'raw_SinM');
                appliedCorrections = [];
                appliedCorrections.gapFilling = true;
                options = applyGapFilling(options);
            end
        end
        % Trues
    elseif options.reconstruct_trues && options.use_machine == 0
        if partitions == 1 && isfield(options, 'SinTrues') == 0
            options.SinM = loadStructFromFile(sinoFile,'SinTrues');
        elseif isfield(options, 'SinTrues') == 0
            options.SinM = loadStructFromFile(sinoFile, 'SinTrues');
        end
        if options.fill_sinogram_gaps && options.det_per_ring < options.det_w_pseudo
            if options.verbose
                disp('Performing sinogram gap filling on trues data')
            end
            options = applyGapFilling(options);
        end
        % Scatter
    elseif options.reconstruct_scatter && options.use_machine == 0
        if partitions == 1 && isfield(options, 'SinScatter') == 0
            options.SinM = loadStructFromFile(sinoFile,'SinScatter');
        elseif isfield(options, 'SinScatter') == 0
            options.SinM = loadStructFromFile(sinoFile, 'SinScatter');
        end
        if options.fill_sinogram_gaps && options.det_per_ring < options.det_w_pseudo
            if options.verbose
                disp('Performing sinogram gap filling on scatter data')
            end
            options = applyGapFilling(options);
        end
    end
    % Load randoms data
    if options.randoms_correction && ~options.reconstruct_scatter && ~options.reconstruct_trues && (~isfield(options,'SinDelayed') || numel(options.SinDelayed) == 0)
        if options.use_machine == 0 || options.use_machine == 1 || options.use_machine == 3
            try
                [options.SinDelayed,RandProp] = loadStructFromFile(sinoFile,'SinDelayed','RandProp');
            catch
                options = loadDelayedData(options);
            end
        else
            options = loadDelayedData(options);
        end
    end
    if partitions == 1 && options.randoms_correction && ~options.reconstruct_scatter && ~options.reconstruct_trues
        if iscell(options.SinDelayed)
            if numel(options.SinDelayed{1}) <= 1
                warning('No randoms correction data detected, disabling randoms correction!')
                options.randoms_correction = false;
            end
        else
            if numel(options.SinDelayed) <= 1
                warning('No randoms correction data detected, disabling randoms correction!')
                options.randoms_correction = false;
            end
        end
    elseif options.randoms_correction && ~options.reconstruct_scatter && ~options.reconstruct_trues
        if length(options.SinDelayed) < partitions && iscell(options.SinDelayed)
            warning('The number of delayed coincidence sinograms is less than the number of time points. Using the first one')
            temp = options.SinDelayed;
            options.SinDelayed = cell(partitions,1);
            if sum(size(temp{1})) > 1
                if size(temp{1},1) ~= size(options.Nang)
                    temp{1} = permute(temp{1},[2 1 3]);
                end
            end
            for kk = 1 : partitions
                options.SinDelayed{kk} = temp{1};
            end
        elseif partitions > 1 && ~iscell(options.SinDelayed)
            warning('The number of delayed coincidence sinograms is less than the number of time points. Using the first one')
            temp = options.SinDelayed;
            options.SinDelayed = cell(partitions,1);
            if sum(size(temp)) > 1
                if size(temp,1) ~= size(options.Nang)
                    temp = permute(temp,[2 1 3]);
                end
            end
            for kk = 1 : partitions
                options.SinDelayed{kk} = temp;
            end
        else
            if iscell(options.SinDelayed)
                for kk = 1 : length(options.SinDelayed)
                    if sum(size(options.SinDelayed{kk})) > 1
                        if size(options.SinDelayed{kk},1) ~= size(options.Nang)
                            options.SinDelayed{kk} = permute(options.SinDelayed{kk},[2 1 3]);
                        end
                    end
                end
            else
                if sum(size(options.SinDelayed)) > 1
                    if size(options.SinDelayed,1) ~= size(options.Nang)
                        options.SinDelayed = permute(options.SinDelayed,[2 1 3]);
                    end
                end
            end
        end
        if iscell(options.SinDelayed)
            if numel(options.SinDelayed{1}) <= 1
                warning('No randoms correction data detected, disabling randoms correction!')
                options.randoms_correction = false;
            end
        else
            if numel(options.SinDelayed) <= 1
                warning('No randoms correction data detected, disabling randoms correction!')
                options.randoms_correction = false;
            end
        end
    end
end
if options.TOF && options.TOF_bins_used == 1
    options.TOF_bins = options.TOF_bins_used;
    if iscell(options.SinM)
        for kk = 1 : numel(options.SinM)
            options.SinM{kk} = sum(options.SinM{kk}, 4,'native');
        end
    else
        options.SinM = sum(options.SinM, 4,'native');
    end
    options.TOF = false;
end
if nargout > 1
    varargout{1} = RandProp;
    varargout{2} = ScatterProp;
end