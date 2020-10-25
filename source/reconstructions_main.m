function pz = reconstructions_main(options)
%% Main reconstruction file
% This function is used to compute various reconstructions with the
% selected method. Can be used with any sinogram or raw data.
%
% OUTPUT:
%   pz = A cell matrix containing output from each of the selected
%   algorithms and/or priors. E.g. if OSEM and ROSEM are selected pz{2}
%   contains the OSEM estimates and pz{5} the ROSEM estimates.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ville-Veikko Wettenhovi, Samuli Summala
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(options,'use_machine')
    options.use_machine = 0;
end

if ~isfield(options, 'custom')
    options.custom = false;
end

if ~isfield(options,'attenuation_phase')
    options.attenuation_phase = false;
end

if ~isfield(options,'save_iter')
    options.save_iter = true;
end

if ~options.use_raw_data && isfield(options,'coincidences')
    options = rmfield(options, 'coincidences');
end

if ~isfield(options,'TOF_bins')
    options.TOF_bins = 1;
end
if ~isfield(options,'TOF_width')
    options.TOF_width = 0;
end

TOF = options.TOF_bins > 1 && options.projector_type == 1;

folder = fileparts(which('reconstructions_main.m'));
folder = strrep(folder, 'source','mat-files/');
folder = strrep(folder, '\','/');

if options.span == 1
    options.TotSinos = options.rings^2;
    options.NSinos = options.TotSinos;
end

disp('Preparing for reconstruction')

% Load the measurement data if it does not exist in options.SinM
% Raw data
if options.use_raw_data
    RandProp.smoothing = false;
    RandProp.variance_reduction = false;
    ScatterProp.smoothing = false;
    if options.partitions == 1
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
        load_string = [options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_' ...
            num2str(options.tot_time) 's_raw'];
        load_string2 = [options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
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
    if options.reconstruct_trues == false && ~options.reconstruct_scatter && (isfield(options, 'coincidences') == 0 || options.precompute_all) && options.use_machine < 2
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
            if iscell(options.SinDelayed)
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
            if isfield(options,'SinDelayed')
                if iscell(options.SinM) && iscell(options.SinDelayed)
                    for kk = 1 : length(options.SinM)
                        if options.variance_reduction
                            options.SinDelayed{kk} = Randoms_variance_reduction(single(options.SinDelayed{kk}), options);
                            RandProp.variance_reduction = true;
                        end
                        if options.randoms_smoothing
                            options.SinDelayed{kk} = randoms_smoothing(options.SinDelayed{kk}, options);
                            RandProp.smoothing = true;
                        end
                        if ~options.corrections_during_reconstruction
                            options.SinM{kk} = options.SinM{kk} - single(options.SinDelayed{kk} / options.TOF_bins);
                        end
                        if ~options.corrections_during_reconstruction
                            options.SinM{kk}(options.SinM{kk} < 0) = 0;
                        end
                    end
                else
                    if options.variance_reduction
                        options.SinDelayed = Randoms_variance_reduction(single(options.SinDelayed), options);
                        RandProp.variance_reduction = true;
                    end
                    if options.randoms_smoothing
                        options.SinDelayed = randoms_smoothing(options.SinDelayed, options);
                        RandProp.smoothing = true;
                    end
                    if ~options.corrections_during_reconstruction
                        options.SinM = options.SinM - single(options.SinDelayed / options.TOF_bins);
                    end
                    if ~options.corrections_during_reconstruction
                        options.SinM(options.SinM < 0) = 0;
                    end
                end
                if ~options.corrections_during_reconstruction
                    options = rmfield(options,'SinDelayed');
                end
            else
                disp('Delayed coincidences not found, randoms correction not performed')
                options.SinDelayed = 0;
                options.randoms_correction = false;
            end
        else
            if isfield(options,'SinDelayed')
                options = loadDelayedData(options);
                if iscell(options.SinDelayed)
                    if iscell(options.SinM)
                        for kk = 1 : length(options.SinM)
                            if numel(options.SinDelayed{kk}) ~= numel(options.SinM{kk})
                                error('Size mismatch between randoms correction data and measurement data')
                            end
                            if options.variance_reduction
                                options.SinDelayed{kk} = Randoms_variance_reduction(single(options.SinDelayed{kk}), options);
                                RandProp.variance_reduction = true;
                            end
                            if options.randoms_smoothing
                                options.SinDelayed{kk} = randoms_smoothing(double(options.SinDelayed{kk}), options);
                                RandProp.smoothing = true;
                            end
                            if ~options.corrections_during_reconstruction
                                options.SinM{kk} = options.SinM{kk} - single(options.SinDelayed{kk} / options.TOF_bins);
                            end
                            if ~options.corrections_during_reconstruction
                                options.SinM{kk}(options.SinM{kk} < 0) = 0;
                            end
                        end
                    else
                        if numel(options.SinDelayed{1}) ~= numel(options.SinM)
                            error('Size mismatch between randoms correction data and measurement data')
                        end
                        options.SinM = options.SinM - single(options.SinDelayed{1}(:) / options.TOF_bins);
                        options.SinM(options.SinM < 0) = 0;
                    end
                else
                    if iscell(options.SinM)
                        for kk = 1 : length(options.SinM)
                            if numel(options.SinDelayed) ~= numel(options.SinM{kk})
                                error('Size mismatch between randoms correction data and measurement data')
                            end
                            if options.variance_reduction
                                options.SinDelayed = Randoms_variance_reduction(single(options.SinDelayed), options);
                                RandProp.variance_reduction = true;
                            end
                            if options.randoms_smoothing
                                options.SinDelayed = randoms_smoothing(single(options.SinDelayed), options);
                                RandProp.smoothing = true;
                            end
                            if ~options.corrections_during_reconstruction
                                options.SinM{kk} = (options.SinM{kk}) - single(options.SinDelayed(:) / options.TOF_bins);
                            end
                            if ~options.corrections_during_reconstruction
                                options.SinM{kk}(options.SinM{kk} < 0) = 0;
                            end
                        end
                    else
                        if numel(options.SinDelayed) ~= numel(options.SinM)
                            error('Size mismatch between randoms correction data and measurement data')
                        end
                        if options.variance_reduction
                            options.SinDelayed = Randoms_variance_reduction(single(options.SinDelayed), options);
                            RandProp.variance_reduction = true;
                        end
                        if options.randoms_smoothing
                            options.SinDelayed = randoms_smoothing(options.SinDelayed, options);
                            RandProp.smoothing = true;
                        end
                        if ~options.corrections_during_reconstruction
                            options.SinM = options.SinM - single(options.SinDelayed(:) / options.TOF_bins);
                        end
                        if ~options.corrections_during_reconstruction
                            options.SinM(options.SinM < 0) = 0;
                        end
                    end
                end
            end
        end
    end
    if options.scatter_correction && ~options.corrections_during_reconstruction
        if ~isfield(options,'ScatterC')
            options = loadScatterData(options);
        end
        if iscell(options.ScatterC)
            if iscell(options.SinM)
                for kk = 1 : length(options.SinM)
                    if numel(options.ScatterC{kk}) ~= numel(options.SinM{kk})
                        error('Size mismatch between options.scatter correction data and measurement data')
                    end
                    if options.scatter_variance_reduction
                        options.ScatterC{kk} = Randoms_variance_reduction(single(options.ScatterC{kk}), options);
                    end
                    if options.scatter_smoothing
                        options.ScatterC{kk} = randoms_smoothing(options.ScatterC{kk}, options);
                        if options.subtract_scatter
                            options.SinM{kk} = full(options.SinM{kk}) - single(options.ScatterC{kk}(:));
                        else
                            options.SinM{kk} = full(options.SinM{kk}) .* single(options.ScatterC{kk}(:));
                        end
                        ScatterProp.smoothing = true;
                    else
                        if options.subtract_scatter
                            options.SinM{kk} = full(options.SinM{kk}) - single(options.ScatterC{kk}(:));
                        else
                            options.SinM{kk} = full(options.SinM{kk}) .* single(options.ScatterC{kk}(:));
                        end
                    end
                    options.SinM{kk}(options.SinM{kk} < 0) = 0;
                end
            else
                if numel(options.ScatterC{1}) ~= numel(options.SinM)
                    error('Size mismatch between options.scatter correction data and measurement data')
                end
                if options.scatter_variance_reduction
                    options.ScatterC{1} = Randoms_variance_reduction(single(options.ScatterC{1}), options);
                end
                if options.scatter_smoothing
                    options.ScatterC{1} = randoms_smoothing(single(options.ScatterC{1}), options);
                    if options.subtract_scatter
                        options.SinM = full(options.SinM) - single(options.ScatterC{1}(:));
                    else
                        options.SinM = full(options.SinM) .* single(options.ScatterC{1}(:));
                    end
                    ScatterProp.smoothing = true;
                else
                    if options.subtract_scatter
                        options.SinM = (options.SinM) - single(options.ScatterC{1}(:));
                    else
                        options.SinM = (options.SinM) .* single(options.ScatterC{1}(:));
                    end
                end
                options.SinM(options.SinM < 0) = 0;
            end
        else
            if iscell(options.SinM)
                for kk = 1 : length(options.SinM)
                    if numel(options.ScatterC) ~= numel(options.SinM{kk})
                        error('Size mismatch between options.scatter correction data and measurement data')
                    end
                    if options.scatter_variance_reduction
                        options.ScatterC = Randoms_variance_reduction(single(options.ScatterC), options);
                    end
                    if options.scatter_smoothing
                        options.ScatterC = randoms_smoothing(options.ScatterC, options);
                        if options.subtract_scatter
                            options.SinM{kk} = full(options.SinM{kk}) - single(options.ScatterC(:));
                        else
                            options.SinM{kk} = full(options.SinM{kk}) .* single(options.ScatterC(:));
                        end
                        ScatterProp.smoothing = true;
                    else
                        if options.subtract_scatter
                            options.SinM{kk} = (options.SinM{kk}) - single(options.ScatterC(:));
                        else
                            options.SinM{kk} = (options.SinM{kk}) .* single(options.ScatterC(:));
                        end
                    end
                    options.SinM{kk}(options.SinM{kk} < 0) = 0;
                end
            else
                if numel(options.ScatterC) ~= numel(options.SinM)
                    error('Size mismatch between options.scatter correction data and measurement data')
                end
                if options.scatter_variance_reduction
                    options.ScatterC = Randoms_variance_reduction(single(options.ScatterC), options);
                end
                if options.scatter_smoothing
                    options.ScatterC = randoms_smoothing(options.ScatterC, options);
                    if options.subtract_scatter
                        options.SinM = full(options.SinM) - single(options.ScatterC(:));
                    else
                        options.SinM = full(options.SinM) .* single(options.ScatterC(:));
                    end
                    ScatterProp.smoothing = true;
                else
                    if options.subtract_scatter
                        options.SinM = options.SinM - single(options.ScatterC(:));
                    else
                        options.SinM = options.SinM .* single(options.ScatterC(:));
                    end
                end
                options.SinM(options.SinM < 0) = 0;
            end
        end
    end
    if options.normalization_correction && ~options.corrections_during_reconstruction
        if options.use_user_normalization
            [file, fpath] = uigetfile({'*.mat'},'Select normalization datafile');
            if isequal(file, 0)
                error('No file was selected')
            end
            if any(strfind(file, '.nrm'))
                error('Inveon normalization data cannot be used with raw list-mode data')
            else
                data = load([fpath file]);
                variables = fieldnames(data);
                normalization = data.(variables{1});
                clear data
                if numel(normalization) ~= sum(1:options.detectors)
                    error('Size mismatch between the current data and the normalization data file')
                end
            end
        else
            norm_file = [folder options.machine_name '_normalization_listmode.mat'];
            if exist(norm_file, 'file') == 2
                normalization = loadStructFromFile(norm_file,'normalization');
            else
                error('Normalization correction selected, but no normalization data found')
            end
        end
        if iscell(options.SinM)
            for kk = 1 : length(options.SinM)
                options.SinM{kk} = options.SinM{kk} .* single(full(normalization));
            end
        else
            options.SinM = options.SinM .* single(full(normalization));
        end
    end
    
    clear coincidences options.coincidences true_coincidences delayed_coincidences
    % Sinogram data
else
    RandProp.smoothing = false;
    RandProp.variance_reduction = false;
    ScatterProp.smoothing = false;
    ScatterProp.variance_reduction = false;
    ScatterProp.normalization = false;
    if options.partitions == 1
        if TOF
            load_string = [options.machine_name '_' options.name '_TOFsinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span)];
            load_string_TOF = ['_' num2str(options.TOF_bins) 'bins_' num2str(options.TOF_width*1e12) 'psBinSize_' num2str(options.TOF_noise_FWHM*1e12) 'psFWHM'];
            load_string = [load_string load_string_TOF];
        else
            load_string = [options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                num2str(options.TotSinos) '_span' num2str(options.span)];
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
        if TOF
            load_string = [options.machine_name '_' options.name '_TOFsinograms_combined_' num2str(poptions.artitions) 'timepoints_for_total_of_' num2str(options.tot_time) 's_' ...
                num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span)];
            load_string_TOF = ['_' num2str(options.TOF_bins) 'bins_' num2str(options.TOF_width*1e12) 'psBinSize_' num2str(options.TOF_noise_FWHM*1e12) 'psFWHM'];
            load_string = [load_string load_string_TOF];
            load_string2 = load_string;
        else
            load_string = [options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_' ...
                num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                num2str(options.span)];
            load_string2 = [options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
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
    if ~options.reconstruct_trues && ~options.reconstruct_scatter
        if isfield(options,'SinM') == 0
            if ((options.randoms_correction || options.scatter_correction || options.normalization_correction) && ~options.corrections_during_reconstruction) ...
                    || options.fill_sinogram_gaps
                try
                    [options.SinM, appliedCorrections] = loadStructFromFile(sinoFile, 'SinM','appliedCorrections');
                catch
                    [options.SinM, appliedCorrections] = loadStructFromFile(sinoFile, 'raw_SinM','appliedCorrections');
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
                %                 appliedCorrections.normalization = false;
                %                 appliedCorrections.randoms = false;
                %                 appliedCorrections.scatter = false;
                %                 appliedCorrections.gapFilling = false;
            end
        end
        if ~options.corrections_during_reconstruction && ~isempty(appliedCorrections)
            if ~appliedCorrections.normalization && options.normalization_correction
                error('Normalization correction selected, but no normalization correction applied to the sinogram')
            end
            if (strcmp(appliedCorrections.randoms, 'randoms correction') || strcmp(appliedCorrections.randoms, 'randoms correction with smoothing')) ...
                    && options.variance_reduction
                error('Randoms variance correction selected, but no variance reduction has been performed')
            end
            if (strcmp(appliedCorrections.randoms, 'randoms correction') || strcmp(appliedCorrections.randoms, 'randoms correction with variance reduction')) ...
                    && options.randoms_smoothing
                error('Randoms smoothing selected, but no smoothing has been performed')
            end
            if (strcmp(appliedCorrections.scatter, 'scatter correction') || strcmp(appliedCorrections.scatter, 'scatter correction with smoothing')) ...
                    && options.scatter_variance_reduction
                error('Scatter variance correction selected, but no variance reduction has been performed')
            end
            if (strcmp(appliedCorrections.scatter, 'scatter correction') || strcmp(appliedCorrections.scatter, 'scatter correction with variance reduction')) ...
                    && options.scatter_smoothing
                error('Scatter smoothing selected, but no smoothing has been performed')
            end
            if ~isempty(strfind(appliedCorrections.scatter,'normalized scatter correction')) && options.normalized_scatter
                error('Normalized scatter correction selected, but no normalization has been performed')
            end
            if ~appliedCorrections.gapFilling && options.fill_sinogram_gaps
                appliedCorrections.gapFilling = true;
                [~, ~, xp, yp] = detector_coordinates(options);
                for llo = 1 : options.partitions
                    if llo == 1
                        gaps = [];
                    end
                    if options.partitions > 1
                        Sin = options.SinM{llo};
                    else
                        Sin = options.SinM;
                    end
                    [Sin, gaps] = gapFilling(options, Sin, xp, yp, llo, gaps);
                    if options.partitions > 1
                        options.SinM{llo} = Sin;
                    else
                        options.SinM = Sin;
                    end
                end
                clear Sin
            end
        elseif options.corrections_during_reconstruction && ~isempty(appliedCorrections)
            if (appliedCorrections.normalization || ~isempty(appliedCorrections.randoms) || ~isempty(appliedCorrections.scatter) || ~appliedCorrections.gapFilling) ...
                    && options.fill_sinogram_gaps && options.det_per_ring < options.det_w_pseudo
                options.SinM = loadStructFromFile(sinoFile, 'raw_SinM');
                appliedCorrections = [];
                appliedCorrections.gapFilling = true;
                [~, ~, xp, yp] = detector_coordinates(options);
                for llo = 1 : options.partitions
                    if llo == 1
                        gaps = [];
                    end
                    if options.partitions > 1
                        Sin = options.SinM{llo};
                    else
                        Sin = options.SinM;
                    end
                    [Sin, gaps] = gapFilling(options, Sin, xp, yp, llo, gaps);
                    if options.partitions > 1
                        options.SinM{llo} = Sin;
                    else
                        options.SinM = Sin;
                    end
                end
                clear Sin
            end
        end
    elseif options.reconstruct_trues && options.use_machine == 0
        if options.partitions == 1 && isfield(options, 'SinTrues') == 0
            options.SinM = loadStructFromFile(sinoFile,'SinTrues');
        elseif isfield(options, 'SinTrues') == 0
            options.SinM = loadStructFromFile(sinoFile, 'SinTrues');
        end
        if options.fill_sinogram_gaps && options.det_per_ring < options.det_w_pseudo
            if options.verbose
                disp('Performing sinogram gap filling on trues data')
            end
            [~, ~, xp, yp] = detector_coordinates(options);
            for llo = 1 : options.partitions
                if llo == 1
                    gaps = [];
                end
                if options.partitions > 1
                    Sin = options.SinM{llo};
                else
                    Sin = options.SinM;
                end
                [Sin, gaps] = gapFilling(options, Sin, xp, yp, llo, gaps);
                if options.partitions > 1
                    options.SinM{llo} = Sin;
                else
                    options.SinM = Sin;
                end
            end
            clear Sin
        end
    elseif options.reconstruct_scatter && options.use_machine == 0
        if options.partitions == 1 && isfield(options, 'SinScatter') == 0
            options.SinM = loadStructFromFile(sinoFile,'SinScatter');
        elseif isfield(options, 'SinScatter') == 0
            options.SinM = loadStructFromFile(sinoFile, 'SinScatter');
        end
        if options.fill_sinogram_gaps && options.det_per_ring < options.det_w_pseudo
            if options.verbose
                disp('Performing sinogram gap filling on scatter data')
            end
            [~, ~, xp, yp] = detector_coordinates(options);
            for llo = 1 : options.partitions
                if llo == 1
                    gaps = [];
                end
                if options.partitions > 1
                    Sin = options.SinM{llo};
                else
                    Sin = options.SinM;
                end
                [Sin, gaps] = gapFilling(options, Sin, xp, yp, llo, gaps);
                if options.partitions > 1
                    options.SinM{llo} = Sin;
                else
                    options.SinM = Sin;
                end
            end
            clear Sin
        end
    end
    if options.partitions == 1 && options.randoms_correction && options.corrections_during_reconstruction && ~options.reconstruct_scatter && ~options.reconstruct_trues
        if options.use_machine == 0
            [options.SinDelayed,RandProp] = loadStructFromFile(sinoFile,'SinDelayed','RandProp');
        elseif options.use_machine == 1 || options.use_machine == 3
            [options.SinDelayed,RandProp] = loadStructFromFile(sinoFile,'SinDelayed','RandProp');
        else
            options = loadDelayedData(options);
        end
        if iscell(options.SinDelayed)
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
    elseif options.randoms_correction && options.corrections_during_reconstruction && ~options.reconstruct_scatter && ~options.reconstruct_trues
        if options.use_machine == 0
            [options.SinDelayed,RandProp] = loadStructFromFile(sinoFile, 'SinDelayed','RandProp');
        elseif options.use_machine == 1 || options.use_machine == 3
            [options.SinDelayed,RandProp] = loadStructFromFile(sinoFile, 'SinDelayed','RandProp');
        else
            options = loadDelayedData(options);
            if length(options.SinDelayed) < options.partitions && iscell(options.SinDelayed)
                warning('The number of delayed coincidence sinograms is less than the number of time points. Using the first one')
                temp = options.SinDelayed;
                options.SinDelayed = cell(options.partitions,1);
                if sum(size(temp{1})) > 1
                    if size(temp{1},1) ~= size(options.Nang)
                        temp{1} = permute(temp{1},[2 1 3]);
                    end
                end
                for kk = 1 : options.partitions
                    options.SinDelayed{kk} = temp{1};
                end
            elseif options.partitions > 1
                warning('The number of delayed coincidence sinograms is less than the number of time points. Using the first one')
                temp = options.SinDelayed;
                options.SinDelayed = cell(options.partitions,1);
                if sum(size(temp)) > 1
                    if size(temp,1) ~= size(options.Nang)
                        temp = permute(temp,[2 1 3]);
                    end
                end
                for kk = 1 : options.partitions
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
        end
        if iscell(options.SinDelayed)
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
    end
end

rekot = reko_maker(options);

Nx = options.Nx;
Ny = options.Ny;
Nz = options.Nz;
Niter = options.Niter;
subsets = options.subsets;
epps = options.epps;
precompute_obs_matrix = options.precompute_obs_matrix;
attenuation_correction = options.attenuation_correction;
diameter = options.diameter;
FOVax = options.FOVa_x;
FOVay = options.FOVa_y;
axial_fov = options.axial_fov;
NSinos = options.NSinos;
pseudot = uint32(options.pseudot);
rings = options.rings;
if options.use_raw_data && isfield(options,'x')
    det_per_ring = numel(options.x);
else
    det_per_ring = options.det_per_ring;
end
if options.use_raw_data
    rings = rings - sum(options.pseudot);
    options.rings = rings;
    options.detectors = det_per_ring * rings;
end
machine_name = options.machine_name;
Nang = options.Nang;
Ndist = options.Ndist;
use_raw_data = options.use_raw_data;
TotSinos = options.TotSinos;
partitions = options.partitions;
verbose = options.verbose;
device = uint32(options.use_device);
options.empty_weight = false;
options.MBSREM_prepass = true;
options.tr_offsets = 0;
list_mode_format = false;
options.MAP = (options.OSL_MLEM || options.OSL_OSEM || options.BSREM || options.MBSREM || options.ROSEM_MAP || options.RBI_OSL...
    || any(options.COSEM_OSL));

MLEM_bool = options.OSL_MLEM || options.mlem;
OS_bool = options.osem || options.rosem || options.ramla || options.OSL_OSEM || options.BSREM || options.ROSEM_MAP || options.rbi || options.drama ...
    || options.cosem || options.ecosem || options.acosem || options.RBI_OSL || any(options.COSEM_OSL);

pz = cell(length(rekot),partitions);

N = Nx * Ny * Nz;

temp = pseudot;
if ~isempty(temp) && temp > 0
    for kk = uint32(1) : temp
        pseudot(kk) = uint32(options.cryst_per_block + 1) * kk;
    end
elseif temp == 0
    pseudot = [];
end



if options.implementation == 3
    if options.osem && options.mlem
        if subsets == 1
            disp(['Both OSEM and MLEM set for method ' num2str(options.implementation) ', using MLEM'])
            options.osem = false;
        else
            disp(['Both OSEM and MLEM set for method ' num2str(options.implementation) ', using OSEM'])
            options.mlem = false;
        end
    end
end

if (options.quad || options.FMH || options.L || options.weighted_mean || options.MRP || options.TV || options.Huber) && options.MAP
    Ndx = options.Ndx;
    Ndy = options.Ndy;
    Ndz = options.Ndz;
end
if options.L && options.MAP
    if ~isempty(options.a_L)
        if length(options.a_L(:)) < ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
            error(['Weights vector options.a_L is too small, needs to be ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) ' in length'])
        elseif length(options.a_L(:)) > ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
            error(['Weights vector options.a_L is too large, needs to be ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) ' in length'])
        end
    end
end
if (options.quad || options.FMH || options.L || options.weighted_mean || options.Huber) && options.MAP
    if ~isempty(options.weights)
        if length(options.weights(:)) < ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
            error(['Weights vector is too small, needs to be ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) ' in length'])
        elseif length(options.weights(:)) > ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
            error(['Weights vector is too large, needs to be ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) ' in length'])
        end
        if ~isinf(options.weights(ceil(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)/2))))
            options.weights(ceil(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)/2))) = Inf;
        end
    else
        options.empty_weight = true;
    end
end
if options.Huber && options.MAP
    if ~isempty(options.weights_huber)
        if length(options.weights_huber(:)) < ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
            error(['Huber weights vector is too small, needs to be ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) ' in length'])
        elseif length(options.weights_huber(:)) > ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
            error(['Huber weights vector is too large, needs to be ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) ' in length'])
        end
        if ~isinf(options.weights_huber(ceil(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)/2))))
            options.weights_huber(ceil(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)/2))) = Inf;
        end
    end
end
if options.FMH && options.MAP
    if ~isempty(options.fmh_weights)
        if Nz == 1 || Ndz == 0
            if length(options.fmh_weights(:)) < (4*(Ndx*2+1))
                error(['Weights vector options.fmh_weights is too small, needs to be [' num2str(Ndx*2+1) ', 4] in size'])
            elseif length(options.fmh_weights(:)) > (4*(Ndx*2+1))
                error(['Weights vector options.fmh_weights is too large, needs to be [' num2str(Ndx*2+1) ', 4] in size'])
            end
        else
            if length(options.fmh_weights(:)) < (13*(Ndx*2+1))
                error(['Weights vector options.fmh_weights is too small, needs to be [' num2str(Ndx*2+1) ', 13] in size'])
            elseif length(options.fmh_weights(:)) > (13*(Ndx*2+1))
                error(['Weights vector options.fmh_weights is too large, needs to be [' num2str(Ndx*2+1) ', 13] in size'])
            end
        end
    end
end
if options.weighted_mean && options.MAP
    if ~isempty(options.weighted_weights)
        if length(options.weighted_weights(:)) < ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
            error(['Weights vector options.weighted_weights is too small, needs to be ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) ' in length'])
        elseif length(options.weighted_weights(:)) > ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
            error(['Weights vector options.weighted_weights is too large, needs to be ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) ' in length'])
        end
        if ~isinf(options.weighted_weights(ceil(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)/2))))
            options.weighted_weights(ceil(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)/2))) = Inf;
        end
    end
end

if options.implementation == 1 || options.implementation == 4
    im_vectors = form_image_vectors(options, N);
end

if options.precompute_lor
    is_transposed = true;
else
    is_transposed = false;
end

if (options.implementation == 3 || options.implementation == 1 ) && ~any(rekot(cellfun('isempty',strfind(algorithms_char(),'MLEM'))))
    subsets = 1;
end

if TOF
    if options.TOF_bins_used ~= options.TOF_bins
        if iscell(options.SinM)
            for kk = 1 : options.partitions
                options.SinM{kk} = sum(options.SinM{kk},4,'native');
            end
        else
            options.SinM = sum(options.SinM,4,'native');
        end
        sigma_x = 0;
        TOFCenter = 0;
        TOF = false;
        options.TOF_bins = 1;
    else
        c = 2.99792458e11;
        sigma_x = (c*options.TOF_FWHM/2) / (2 * sqrt(2 * log(2)));
        edges_user = linspace(-options.TOF_width * options.TOF_bins/2, options.TOF_width * options.TOF_bins / 2, options.TOF_bins + 1);
        edges_user = edges_user(1:end-1) + options.TOF_width/2; % the most probable value where annihilation occured
        TOFCenter = zeros(size(edges_user));
        TOFCenter(1) = edges_user(ceil(length(edges_user)/2));
        TOFCenter(2:2:end) = edges_user(ceil(length(edges_user)/2) + 1:end);
        TOFCenter(3:2:end) = edges_user(ceil(length(edges_user)/2) - 1: -1 : 1);
        if isfield(options, 'TOF_offset') && options.TOF_offset > 0
            TOFCenter = TOFCenter + options.TOF_offset;
        end
        TOFCenter = -TOFCenter * c / 2;
    end
else
    sigma_x = 0;
    TOFCenter = 0;
end
if options.implementation == 2 || options.implementation == 3
    sigma_x = single(sigma_x);
    TOFCenter = single(TOFCenter);
end


%%

% Compute the indices for the subsets used.
% For Sinogram data, six different methods to select the subsets are
% available. For raw list-mode data, three methods are available.
if (isfield(options,'x') && isfield(options,'y') && (isfield(options,'z') || isfield(options,'z_det'))) && numel(options.x) / 2 == numel(options.SinM)
    %     index = uint32(1:numel(options.SinM))';
    index = 0;
    det_per_ring = numel(options.SinM);
    pituus = floor(det_per_ring / subsets);
    pituus = int64([repmat(pituus,subsets - 1,1);det_per_ring - pituus*(subsets - 1)]);
    list_mode_format = true;
    if options.implementation ~= 4
        error('List-mode reconstruction with custom detectors is currently supported only with implementation 4.')
    end
else
    [index, pituus, subsets] = index_maker(Nx, Ny, Nz, subsets, use_raw_data, machine_name, options, Nang, Ndist, TotSinos, NSinos);
end

%%

% Diameter of the PET-device (bore) (mm)
R = double(diameter);
% Transaxial FOV (x-direction, horizontal) (mm)
FOVax = double(FOVax);
% Transaxial FOV (y-direction, vertical) (mm)
FOVay = double(FOVay);
% Axial FOV (mm)
axial_fow = double(axial_fov);
% Number of rings
blocks = uint32(rings - 1);
% First ring
block1=uint32(0);

NSinos = uint32(NSinos);
NSlices = uint32(Nz);

% Coordinates of the detectors
[x, y, z_det, options] = get_coordinates(options, blocks, pseudot);

if ~list_mode_format
    % Load correction data
    [normalization_correction, randoms_correction, options] = set_up_corrections(options, blocks, RandProp, ScatterProp);
else
    normalization_correction = false;
    randoms_correction = false;
end

if options.use_raw_data
    size_x = uint32(numel(x));
else
    size_x = uint32(options.Nang*options.Ndist);
    if options.sampling > 1 && ~options.precompute_lor
        size_x = size_x * options.sampling;
    end
end

if subsets > 1
    pituus = [int64(0);cumsum(pituus)];
    if iscell(index)
        index = cell2mat(index);
    end
end

% Compute the necessary indices required for subsets (e.g. the index of
% the detector coordinates for the current LOR)
if ~list_mode_format
    [options, lor_a, xy_index, z_index, LL, summa, pituus, options.SinM, lor_orth] = form_subset_indices(options, pituus, subsets, index, size_x, y, z_det, rings, false, TOF, ...
        options.SinM);
else
    LL = uint16(0);
    xy_index = uint32(0);
    z_index = uint16(0);
    lor_orth = uint16(0);
    summa = zeros(subsets, 1, 'uint64');
end
if ~options.precompute_lor
    lor_a = uint16(0);
end

if use_raw_data
    if isempty(pseudot)
        pseudot = uint32(1e5);
    else
        pseudot = pseudot - 1;
    end
end


% Pixel boundaries
etaisyys_x = (R - FOVax)/2;
etaisyys_y = (R - FOVay)/2;
if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
    zz = linspace(single(0), single(axial_fow), Nz + 1);
    xx = single(linspace(etaisyys_x, R - etaisyys_x, Nx + 1));
    yy = single(linspace(etaisyys_y, R - etaisyys_y, Ny + 1));
else
    zz = linspace(double(0), double(axial_fow), Nz + 1);
    xx = double(linspace(etaisyys_x, R - etaisyys_x, Nx + 1));
    yy = double(linspace(etaisyys_y, R - etaisyys_y, Ny + 1));
end
%     zz = zz(2*block1 + 1 : 2*blocks + 2);

% Distance of adjacent pixels
dx = diff(xx(1:2));
dy = diff(yy(1:2));
dz = diff(zz(1:2));

% Distance of image from the origin
bx = xx(1);
by = yy(1);
bz = zz(1);

% Number of pixels
Ny = uint32(Ny);
Nx = uint32(Nx);
Nz = uint32(Nz);

N = (Nx)*(Ny)*(Nz);
det_per_ring = uint32(det_per_ring);

% How much memory is preallocated
% Implementation 1, no precompute
if ~use_raw_data
    ind_size = uint32(NSinos / subsets * (det_per_ring) * Nx * (Ny));
else
    ind_size = uint32((numel(x))^2 / subsets * Nx * (Ny));
end

zmax = max(max(z_det));
if zmax==0
    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
        zmax = single(1);
    else
        zmax = double(1);
    end
end

% Coordinates of the centers of the voxels
if options.projector_type == 2 || options.projector_type == 3
    x_center = xx(1 : end - 1)' + dx/2;
    y_center = yy(1 : end - 1)' + dy/2;
    z_center = zz(1 : end - 1)' + dz/2;
    temppi = min([options.FOVa_x / options.Nx, options.axial_fov / options.Nz]);
    if options.tube_width_z > 0
        temppi = max([1,round(options.tube_width_z / temppi)]);
    else
        temppi = max([1,round(options.tube_width_xy / temppi)]);
    end
    temppi = temppi * temppi * 4;
    if options.apply_acceleration || options.implementation == 4
        if options.tube_width_z == 0
            dec = uint32(sqrt(options.Nx^2 + options.Ny^2) * temppi);
        else
            dec = uint32(sqrt(options.Nx^2 + options.Ny^2 + options.Nz^2) * temppi);
        end
    else
        dec = uint32(0);
    end
elseif (options.projector_type == 1 && TOF)
    x_center = xx(1);
    y_center = yy(1);
    z_center = zz(1);
    if (options.apply_acceleration || options.implementation == 4) && options.n_rays_transaxial * options.n_rays_axial == 1
        dec = uint32(sqrt(options.Nx^2 + options.Ny^2 + options.Nz^2) * 2);
    else
        dec = uint32(0);
    end
else
    x_center = xx(1);
    y_center = yy(1);
    z_center = zz(1);
    dec = uint32(0);
end

% Compute the various volumes available for the spherical voxels in
% volume-based ray tracer
if options.projector_type == 3
    dp = max([dx,dy,dz]);
    options.voxel_radius = sqrt(2) * options.voxel_radius * (dp / 2);
    bmax = options.tube_radius + options.voxel_radius;
    b = linspace(0, bmax, 10000)';
    b(options.tube_radius > (b + options.voxel_radius)) = [];
    b = unique(round(b*10^3)/10^3);
    V = volumeIntersection(options.tube_radius, options.voxel_radius, b);
    diffis = [diff(V);0];
    b = b(diffis <= 0);
    V = V(diffis <= 0);
    Vmax = (4*pi)/3*options.voxel_radius^3;
    bmin = min(b);
else
    V = 0;
    Vmax = 0;
    bmin = 0;
    bmax = 0;
end
if options.implementation == 2 || options.implementation == 3
    V = single(V);
    Vmax = single(Vmax);
    bmin = single(bmin);
    bmax = single(bmax);
end

if options.implementation == 1
    iij = double(0:Nx);
    jji = double(0:Ny);
    kkj = double(0:Nz);
end

%% This computes a whole observation matrix and uses it to compute the MLEM (no on-the-fly calculations)
% NOTE: Only attenuation correction is supported
% This section is largely untested
if precompute_obs_matrix && options.implementation == 1
    
    for llo = 1 : partitions
        
        
        if options.mlem
            if TOF
                error('TOF is not supported when using precomputed system matrix')
            end
            if iscell(options.SinM)
                Sino = options.SinM{llo};
            else
                Sino = options.SinM;
                clear options.SinM
            end
            
            Sino = Sino(:);
            
            if issparse(Sino)
                Sino = (full(Sino));
            end
            
            if llo == 1
                options.pituus = pituus;
                options.lor_a = lor_a;
                options.lor_orth = lor_orth;
                options.index = index;
                options.xy_index = xy_index;
                options.LL = LL;
                options.z_index = z_index;
                options.x_center = x_center;
                options.y_center = y_center;
                options.z_center = z_center;
                options.Nx = Nx;
                options.Ny = Ny;
                options.Nz = Nz;
                options.N = N;
                options.zmax = zmax;
                options.summa = summa;
                options.normalization_correction = normalization_correction;
                %             options.randoms_correction = false;
                %             options.scatter_correction = false;
                options.det_per_ring = det_per_ring;
                options.block1 = block1;
                options.blocks = blocks;
                options.NSinos = NSinos;
                options.NSlices = NSlices;
                options.z_det = z_det;
                options.x = x;
                options.y = y;
                options.size_x = size_x;
                options.xx = xx;
                options.yy = yy;
                options.pseudot = pseudot;
                options.ind_size = ind_size;
                options.dx = dx;
                options.dy = dy;
                options.dz = dz;
                options.bx = bx;
                options.by = by;
                options.bz = bz;
                options.bmin = 0;
                options.bmax = 0;
                options.Vmax = 0;
                options.V = 0;
                
                A = observation_matrix_formation(options);
                D = full(sum(A,1))';
                D(D < epps) = epps;
            end
            
            MLEM = ones(N,Niter);
            MLEM(:,1) = options.x0(:);
            Sino = double(Sino);
            
            for iter=1:Niter
                if options.save_iter
                    iter_n = iter + 1;
                else
                    iter_n = 1;
                end
                if options.mlem
                    MLEM(:, iter_n) = MLEM_im(MLEM(:,iter), D, epps, A, Sino, false);
                    disp(['MLEM iteration ' num2str(iter) ' finished'])
                end
                if ~any(rekot(cellfun('isempty',strfind(algorithms_char(),'MLEM'))))
                    warning('Only MLEM is supported with precomputed observation matrix')
                end
            end
            MLEM = reshape(MLEM,Nx,Ny,Nz,Niter+1);
        else
            error('Only MLEM is supported with precomputed observation matrix')
        end
        
        pz{1, llo} = MLEM;
        
    end
    
else
    %% This part is used when the observation matrix is calculated on-the-fly
    
    % Multi-ray Siddon
    % Compute the multi-ray coordinates
    if options.implementation > 1 && (options.n_rays_transaxial > 1 || options.n_rays_axial > 1) && ~options.precompute_lor && options.projector_type == 1
        [x,y] = getMultirayCoordinates(options);
    end
    
    if length(pituus) == 1
        pituus = [int64(0);pituus];
    end
    
    % Remove negative values
    if iscell(options.SinM)
        for llo = 1 : partitions
            options.SinM{llo}(options.SinM{llo} < 0) = 0;
        end
    else
        options.SinM(options.SinM < 0) = 0;
    end
    
    % Compute PSF kernel
    if options.use_psf
        g_pituus_x = ceil(2*(options.FWHM(1) / (2 * sqrt(2 * log(2)))) / dx);
        g_pituus_y = ceil(2*(options.FWHM(2) / (2 * sqrt(2 * log(2)))) / dy);
        g_pituus_z = ceil(2*(options.FWHM(3) / (2 * sqrt(2 * log(2)))) / dz);
        g_x = linspace(-g_pituus_x * dx, g_pituus_x * dx, 2*g_pituus_x + 1)';
        g_y = linspace(-g_pituus_y * dy, g_pituus_y * dy, 2*g_pituus_y + 1);
        g_z = zeros(1,1,g_pituus_z*2+1);
        g_z(1,1,:) = linspace(-g_pituus_z * dz, g_pituus_z * dz, 2*g_pituus_z + 1);
        gaussK = gaussianKernel(g_x, g_y, g_z, options.FWHM(1) / (2 * sqrt(2 * log(2))), options.FWHM(2) / (2 * sqrt(2 * log(2))), options.FWHM(3) / (2 * sqrt(2 * log(2))));
        
        options.g_dim_x = uint32(g_pituus_x);
        options.g_dim_y = uint32(g_pituus_y);
        options.g_dim_z = uint32(g_pituus_z);
        if options.implementation == 2 || options.implementation == 3
            gaussK = single(gaussK(:));
        end
    else
        gaussK = single(0);
    end
    
    % Perform various prepass steps, if necessary
    % These include computing weights, matrices required by some algorithms
    % (COSEM, etc.) and loading anatomic reference images
    [options, D, C_co, C_aco, C_osl, Amin, E] = prepass_phase(options, pituus, index, options.SinM, pseudot, x, y, xx, yy, z_det, dz, dx, dy, bz, bx, by, NSlices, ...
        zmax, size_x, block1, blocks, normalization_correction, randoms_correction, xy_index, z_index, lor_a, lor_orth, summa, LL, is_transposed, x_center, y_center, ...
        z_center, ind_size, gaussK, bmin, bmax, Vmax, V);
    
    if options.use_psf && ((options.mramla || options.MBSREM || options.RBI_OSL || options.rbi) && options.MBSREM_prepass || options.ecosem || options.cosem ...
            || options.acosem || any(options.COSEM_OSL)) && options.implementation == 1
        D = computeConvolution(D, options, Nx, Ny, Nz, gaussK);
    end
    %%
    
    % Compute the reconstructions
    disp('Starting image reconstruction')
    
    % Implementations 1, 4 and 5
    if options.implementation ~= 2 && options.implementation ~= 3
        
        % Loop through all time steps
        for llo = 1 : partitions
            if iscell(options.SinM)
                Sino = options.SinM{llo};
            else
                Sino = options.SinM;
                clear options.SinM
            end
            if options.scatter_correction && ~options.subtract_scatter && iscell(options.ScatterC)
                if iscell(options.ScatterC)
                    options.ScatterC = double(options.ScatterC{1});
                else
                    options.ScatterC = double(options.ScatterC);
                end
            end
            
            Sino = Sino(:);
            
            if issparse(Sino)
                Sino = (full(Sino));
            end
            
            % Implementation 1
            if options.implementation == 1
                % Upper bound for MRAMLA
                if options.MBSREM || options.mramla
                    if options.U == 0 || isempty(options.U)
                        
                        options.U = max(double(Sino)./Amin);
                    end
                end
                % Compute the epsilon value for MRAMLA
                if options.MBSREM || options.mramla
                    if iscell(options.SinDelayed)
                        options.epsilon_mramla = MBSREM_epsilon(Sino, options.epps, randoms_correction, options.SinDelayed{llo}, E);
                    else
                        options.epsilon_mramla = MBSREM_epsilon(Sino, options.epps, randoms_correction, options.SinDelayed, E);
                    end
                end
                
                if ~use_raw_data
                    if isempty(pseudot)
                        pseudot = uint32(0);
                    end
                end
                
                %%
                for iter = 1 : Niter
                    
                    for osa_iter = 1 : subsets
                        
                        % Randoms correction during reconstruction
                        if randoms_correction
                            if iscell(options.SinDelayed)
                                SinD = double(options.SinDelayed{llo}(pituus(osa_iter)+1:pituus(osa_iter + 1)));
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
                        % No precomputation done
                        % This is sequential (non-parallel) code
                        % SLOW
                        % Supports PURE MATLAB computation
                        if options.precompute_lor == false
                            if use_raw_data == false
                                TOFSize = int64(pituus(osa_iter + 1) - pituus(osa_iter));
                                % Siddon
                                if options.projector_type == 1 || options.projector_type == 0
                                    if exist('OCTAVE_VERSION','builtin') == 0 && exist('projector_mex','file') == 3
                                        [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                            zmax, options.vaimennus, norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, ...
                                            randoms_correction, options.scatter, scatter_input, options.global_correction_factor, uint16(0), uint32(0), uint32(0), NSinos, uint16(0), pseudot, det_per_ring, ...
                                            TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                            use_raw_data, uint32(2), ind_size, block1, blocks, index(pituus(osa_iter)+1:pituus(osa_iter + 1)), ...
                                            uint32(options.projector_type), iij, jji, kkj);
                                    elseif exist('OCTAVE_VERSION','builtin') == 5 && exist('projector_oct','file') == 3
                                        [ lor, indices, alkiot] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                            zmax, options.vaimennus, norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, ...
                                            randoms_correction, options.scatter, scatter_input, options.global_correction_factor, uint16(0), uint32(0), uint32(0), NSinos, uint16(0), pseudot, det_per_ring, ...
                                            TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                            use_raw_data, uint32(2), ind_size, block1, blocks, index(pituus(osa_iter)+1:pituus(osa_iter + 1)), ...
                                            uint32(options.projector_type), iij, jji, kkj);
                                    else
                                        % The below lines allow for pure MATLAB
                                        % implemention, i.e. no MEX-files will be
                                        % used. Currently the below function
                                        % uses parfor-loops (requires parallel
                                        % computing toolbox).
                                        % NOTE: The below method is not
                                        % recommended since it is much slower
                                        % method.
                                        [ lor, indices, alkiot, ~] = improved_siddon_atten( int32(Ny), int32(Nx), int32(Nz), dx, dz, by, bx, bz, z_det, x, y, yy, xx, ...
                                            NSinos, NSlices, options.vaimennus, index(pituus(osa_iter)+1:pituus(osa_iter + 1)), pituus(osa_iter + 1) - pituus(osa_iter), ...
                                            attenuation_correction);
                                        alkiot = cell2mat(alkiot);
                                        indices = cell2mat(indices) - 1;
                                        lor = lor(:,2);
                                    end
                                    % Orthogonal distance based
                                elseif options.projector_type == 2
                                    error('Unsupported projector type')
                                else
                                    error('Unsupported projector type')
                                end
                            else
                                L = LL(pituus(osa_iter) * 2 + 1 : pituus(osa_iter + 1) * 2);
                                TOFSize = int64(size(L,1));
                                if options.projector_type == 1
                                    if exist('OCTAVE_VERSION','builtin') == 0
                                        [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                            zmax, options.vaimennus, norm_input, SinD, uint32(0), attenuation_correction, normalization_correction, ...
                                            randoms_correction, options.scatter, scatter_input, options.global_correction_factor, uint16(0), uint32(0), uint32(0), NSinos, L, pseudot, det_per_ring, ...
                                            TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                            use_raw_data, uint32(2), ind_size, block1, blocks, uint32(0), uint32(options.projector_type), iij, jji, kkj);
                                    elseif exist('OCTAVE_VERSION','builtin') == 5
                                        [ lor, indices, alkiot] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                            zmax, options.vaimennus, norm_input, SinD, uint32(0), attenuation_correction, normalization_correction, ...
                                            randoms_correction, options.scatter, scatter_input, options.global_correction_factor, uint16(0), uint32(0), uint32(0), NSinos, L, pseudot, det_per_ring, ...
                                            TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                            use_raw_data, uint32(2), ind_size, block1, blocks, uint32(0), uint32(options.projector_type), iij, jji, kkj);
                                    end
                                elseif options.projector_type == 2
                                    error('Unsupported projector type')
                                else
                                    error('Unsupported projector type')
                                end
                            end
                            if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                                lor = repeat_elem(uint32(1:length(lor))',uint32(lor));
                            elseif exist('OCTAVE_VERSION','builtin') == 5
                                lor = repelem(uint32(1:length(lor)),uint32(lor));
                            else
                                lor = repelem(uint32(1:length(lor)),uint32(lor))';
                            end
                            uu = double(Sino(pituus(osa_iter) + 1 : pituus(osa_iter + 1)));
                            
                            A_length = length(uu);
                            if verbose
                                tStart = tic;
                            end
                            % Form the sparse matrix
                            if exist('OCTAVE_VERSION','builtin') == 0 && ~verLessThan('matlab','9.8')
                                indices = uint32(indices) + 1;
                                A = sparse(lor,indices,double(alkiot), A_length, double(N));
                            elseif options.use_fsparse && exist('fsparse','file') == 3
                                indices = int32(indices) + 1;
                                A = fsparse(int32(lor),(indices),double(alkiot),[A_length double(N) length(alkiot)]);
                            elseif options.use_fsparse && exist('fsparse','file') == 0
                                warning('options.fsparse set to true, but no FSparse mex-file found. Using regular sparse')
                                indices = double(indices) + 1;
                                A = sparse(double(lor),(indices),double(alkiot), A_length, double(N));
                            else
                                indices = double(indices) + 1;
                                A = sparse(double(lor),(indices),double(alkiot), A_length, double(N));
                            end
                            clear indices alkiot lor
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['Sparse matrix formation took ' num2str(tElapsed) ' seconds'])
                            end
                            % Precomputation performed
                            % Parallel
                            % Faster
                            % Only C++ code (no pure MATLAB implementation)
                        else
                            if use_raw_data
                                L_input = LL(pituus(osa_iter) * 2 + 1 : pituus(osa_iter + 1) * 2);
                                xy_index_input = uint32(0);
                                z_index_input = uint16(0);
                                TOFSize = int64(size(L_input,1));
                            else
                                L_input = uint16(0);
                                xy_index_input = xy_index(pituus(osa_iter)+1:pituus(osa_iter + 1));
                                z_index_input = z_index(pituus(osa_iter)+1:pituus(osa_iter + 1));
                                TOFSize = int64(size(xy_index_input,1));
                            end
                            if options.projector_type == 2 || options.projector_type == 3
                                lor2 = [uint64(0); cumsum(uint64(lor_orth(pituus(osa_iter)+1:pituus(osa_iter + 1))))];
                            else
                                lor2 = [uint64(0); cumsum(uint64(lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1))))];
                            end
                            if exist('OCTAVE_VERSION','builtin') == 0
                                [A, ll] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                    norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction,...
                                    randoms_correction, options.scatter, scatter_input, options.global_correction_factor, lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1)), xy_index_input, ...
                                    z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                    TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, use_raw_data, uint32(0), lor2, summa(osa_iter), ...
                                    options.attenuation_phase, uint32(options.projector_type), options.tube_width_xy, x_center, y_center, z_center, ...
                                    options.tube_width_z, bmin, bmax, Vmax, V);
                            elseif exist('OCTAVE_VERSION','builtin') == 5
                                [A, ll] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                    norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction,...
                                    randoms_correction, options.scatter, scatter_input, options.global_correction_factor, lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1)), xy_index_input, ...
                                    z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                    TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, use_raw_data, uint32(0), lor2, summa(osa_iter), ...
                                    options.attenuation_phase, uint32(options.projector_type), options.tube_width_xy, x_center, y_center, z_center, ...
                                    options.tube_width_z, bmin, bmax, Vmax, V);
                            end
                            uu = double(Sino(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                            if options.attenuation_phase
                                uu = uu ./ ll;
                            end
                            clear lor2
                        end
                        % Sensitivity image
                        if is_transposed
                            Summ = full(sum(A,2));
                        else
                            Summ = full(sum(A,1))';
                        end
                        Summ(Summ < options.epps) = options.epps;
                        if options.use_psf
                            Summ = computeConvolution(Summ, options, Nx, Ny, Nz, gaussK);
                        end
                        % Compute OSEM
                        if options.osem || options.ecosem || options.attenuation_phase
                            if verbose
                                tStart = tic;
                            end
                            im_vectors.OSEM_apu = OSEM_im(im_vectors.OSEM_apu, A, epps, uu, Summ, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute MRAMLA
                        if options.mramla
                            if verbose
                                tStart = tic;
                            end
                            im_vectors.MRAMLA_apu = MBSREM(im_vectors.MRAMLA_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, ...
                                SinD, randoms_correction, is_transposed, [], [], options, Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['MRAMLA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute RAMLA
                        if options.ramla
                            if verbose
                                tStart = tic;
                            end
                            im_vectors.RAMLA_apu = BSREM_subiter(im_vectors.RAMLA_apu, options.lam, epps, iter, A, uu, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            if any(im_vectors.RAMLA_apu < 0)
                                warning('Negative values in RAMLA, lower lambda value!')
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RAMLA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute ROSEM
                        if options.rosem
                            if verbose
                                tStart = tic;
                            end
                            im_vectors.ROSEM_apu = ROSEM_subiter(im_vectors.ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute RBI
                        if options.rbi
                            if verbose
                                tStart = tic;
                            end
                            im_vectors.RBI_apu = RBI_subiter(im_vectors.RBI_apu, A, uu, epps, Summ, D, SinD, is_transposed, [], [], options, Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute DRAMA
                        if options.drama
                            if verbose
                                tStart = tic;
                            end
                            im_vectors.DRAMA_apu = DRAMA_subiter(im_vectors.DRAMA_apu, options.lam_drama, epps, iter, Summ, osa_iter, A, uu, SinD, is_transposed, options, ...
                                Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['DRAMA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute COSEM
                        if options.cosem || options.ecosem
                            if verbose
                                tStart = tic;
                            end
                            [im_vectors.COSEM_apu, C_co] = COSEM_im(im_vectors.COSEM_apu, A, epps, uu, C_co, D, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute ECOSEM
                        if options.ecosem
                            if verbose
                                tStart = tic;
                            end
                            im_vectors.ECOSEM_apu = ECOSEM_im(im_vectors.ECOSEM_apu, epps, D, im_vectors.COSEM_apu, im_vectors.OSEM_apu);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ECOSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute ACOSEM
                        if options.acosem
                            if verbose
                                tStart = tic;
                            end
                            [im_vectors.ACOSEM_apu, C_aco] = ACOSEM_im(im_vectors.ACOSEM_apu, A, epps, uu, C_aco, D, options.h, osa_iter, SinD, is_transposed, options, ...
                                Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ACOSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute OSL with MRP
                        if options.MRP && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            med = MRP(im_vectors.MRP_OSL_apu, options.medx, options.medy, options.medz, Nx, Ny, Nz, epps, options.tr_offsets, options.med_no_norm);
                            im_vectors.MRP_OSL_apu = OSEM_im(im_vectors.MRP_OSL_apu, A, epps, uu, OSL(Summ, options.beta_mrp_osem, med, epps), SinD, is_transposed, ...
                                options, Nx, Ny, Nz, gaussK);
                            %                             im_vectors.MRP_OSL_apu = OSL_OSEM(im_vectors.MRP_OSL_apu, Summ, options.beta_mrp_osem, med, epps, A, uu, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute MBSREM with MRP
                        if options.MRP && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            med = MRP(im_vectors.MRP_MBSREM_apu, options.medx, options.medy, options.medz, Nx, Ny, Nz, epps, options.tr_offsets, options.med_no_norm);
                            im_vectors.MRP_MBSREM_apu = MBSREM(im_vectors.MRP_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                                iter, SinD, randoms_correction, is_transposed, options.beta_mrp_mbsrem, med, options, Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['MBSREM MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute BSREM with MRP
                        if options.MRP && options.BSREM
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.ramla
                                im_vectors.MRP_BSREM_apu = im_vectors.RAMLA_apu;
                            else
                                im_vectors.MRP_BSREM_apu = BSREM_subiter(im_vectors.MRP_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed, options, ...
                                    Nx, Ny, Nz, gaussK);
                            end
                            if any(im_vectors.MRP_BSREM_apu < 0)
                                warning('Negative values in BSREM, lower lambda value!')
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['BSREM MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute ROSEM-MAP with MRP
                        if options.MRP && options.ROSEM_MAP
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.rosem
                                im_vectors.MRP_ROSEM_apu = im_vectors.ROSEM_apu;
                            else
                                im_vectors.MRP_ROSEM_apu = ROSEM_subiter(im_vectors.MRP_ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed, ...
                                    options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute RBI-OSL with MRP
                        if options.MRP && options.RBI_OSL
                            if verbose
                                tStart = tic;
                            end
                            med = MRP(im_vectors.MRP_RBI_apu, options.medx, options.medy, options.medz, Nx, Ny, Nz, epps, options.tr_offsets, options.med_no_norm);
                            im_vectors.MRP_RBI_apu = RBI_subiter(im_vectors.MRP_RBI_apu, A, uu, epps, Summ, D, SinD, is_transposed, options.beta_mrp_rbi, med, options, ...
                                Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute COSEM-OSL with MRP
                        if options.MRP && any(options.COSEM_OSL)
                            if verbose
                                tStart = tic;
                            end
                            med = MRP(im_vectors.MRP_COSEM_apu, options.medx, options.medy, options.medz, Nx, Ny, Nz, epps, options.tr_offsets, options.med_no_norm);
                            if options.COSEM_OSL == 1
                                [im_vectors.MRP_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.MRP_COSEM_apu, D, options.beta_mrp_cosem, med, epps, A, uu, ...
                                    C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            else
                                [im_vectors.MRP_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.MRP_COSEM_apu, D, options.beta_mrp_cosem, med, epps, A, uu, ...
                                    C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute OSL with Quadratic prior
                        if options.quad && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            med = Quadratic_prior(im_vectors.Quad_OSL_apu, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz, options);
                            im_vectors.Quad_OSL_apu = OSEM_im(im_vectors.Quad_OSL_apu, A, epps, uu, OSL(Summ, options.beta_quad_osem, med, epps), SinD, is_transposed, ...
                                options, Nx, Ny, Nz, gaussK);
                            %                             im_vectors.Quad_OSL_apu = OSL_OSEM(im_vectors.Quad_OSL_apu, Summ, options.beta_quad_osem, med, epps, A, uu, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute MBSREM with Quadratic prior
                        if options.quad && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            med = Quadratic_prior(im_vectors.Quad_MBSREM_apu, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                            im_vectors.Quad_MBSREM_apu = MBSREM(im_vectors.Quad_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                                iter, SinD, randoms_correction, is_transposed, options.beta_quad_mbsrem, med, options, Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['MBSREM Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute BSREM with Quadratic prior
                        if options.quad && options.BSREM
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.ramla
                                im_vectors.Quad_BSREM_apu = im_vectors.RAMLA_apu;
                            else
                                im_vectors.Quad_BSREM_apu = BSREM_subiter(im_vectors.Quad_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed, options, ...
                                    Nx, Ny, Nz, gaussK);
                            end
                            if any(im_vectors.Quad_BSREM_apu < 0)
                                warning('Negative values in BSREM, it is recommended to lower lambda value')
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['BSREM Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute ROSEM-MAP with Quadratic prior
                        if options.quad && options.ROSEM_MAP
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.rosem
                                im_vectors.Quad_ROSEM_apu = im_vectors.ROSEM_apu;
                            else
                                im_vectors.Quad_ROSEM_apu = ROSEM_subiter(im_vectors.Quad_ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed, ...
                                    options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute RBI-OSL with Quadratic prior
                        if options.quad && options.RBI_OSL
                            if verbose
                                tStart = tic;
                            end
                            med = Quadratic_prior(im_vectors.Quad_RBI_apu, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                            im_vectors.Quad_RBI_apu = RBI_subiter(im_vectors.Quad_RBI_apu, A, uu, epps, Summ, D, SinD, is_transposed, options.beta_quad_rbi, med, ...
                                options, Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-OSL Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute COSEM-OSL with Quadratic prior
                        if options.quad && any(options.COSEM_OSL)
                            if verbose
                                tStart = tic;
                            end
                            med = Quadratic_prior(im_vectors.Quad_COSEM_apu, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                            if options.COSEM_OSL == 1
                                [im_vectors.Quad_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.Quad_COSEM_apu, D, options.beta_quad_cosem, med, epps, A, uu, ...
                                    C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            else
                                [im_vectors.Quad_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.Quad_COSEM_apu, D, options.beta_quad_cosem, med, epps, A, uu, ...
                                    C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-OSL Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute OSL with Huber prior
                        if options.Huber && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            med = Huber_prior(im_vectors.Huber_OSL_apu, options.weights_huber, Nx, Ny, Nz, Ndx, Ndy, Ndz, options.huber_delta);
                            im_vectors.Huber_OSL_apu = OSEM_im(im_vectors.Huber_OSL_apu, A, epps, uu, OSL(Summ, options.beta_huber_osem, med, epps), SinD, is_transposed, ...
                                options, Nx, Ny, Nz, gaussK);
                            %                             im_vectors.Huber_OSL_apu = OSL_OSEM(im_vectors.Huber_OSL_apu, Summ, options.beta_huber_osem, med, epps, A, uu, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute MBSREM with Huber prior
                        if options.Huber && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            med = Huber_prior(im_vectors.Huber_MBSREM_apu, options.weights_huber, Nx, Ny, Nz, Ndx, Ndy, Ndz, options.huber_delta);
                            im_vectors.Huber_MBSREM_apu = MBSREM(im_vectors.Huber_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                                iter, SinD, randoms_correction, is_transposed, options.beta_huber_mbsrem, med, options, Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['MBSREM Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute BSREM with Huber prior
                        if options.Huber && options.BSREM
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.ramla
                                im_vectors.Huber_BSREM_apu = im_vectors.RAMLA_apu;
                            else
                                im_vectors.Huber_BSREM_apu = BSREM_subiter(im_vectors.Huber_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed, options, ...
                                    Nx, Ny, Nz, gaussK);
                            end
                            if any(im_vectors.Huber_BSREM_apu < 0)
                                warning('Negative values in BSREM, it is recommended to lower lambda value')
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['BSREM Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute ROSEM-MAP with Huber prior
                        if options.Huber && options.ROSEM_MAP
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.rosem
                                im_vectors.Huber_ROSEM_apu = im_vectors.ROSEM_apu;
                            else
                                im_vectors.Huber_ROSEM_apu = ROSEM_subiter(im_vectors.Huber_ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed, ...
                                    options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute RBI-OSL with Huber prior
                        if options.Huber && options.RBI_OSL
                            if verbose
                                tStart = tic;
                            end
                            med = Huber_prior(im_vectors.Huber_RBI_apu, options.weights_huber, Nx, Ny, Nz, Ndx, Ndy, Ndz, options.huber_delta);
                            im_vectors.Huber_RBI_apu = RBI_subiter(im_vectors.Huber_RBI_apu, A, uu, epps, Summ, D, SinD, is_transposed, options.beta_huber_rbi, med, ...
                                options, Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-OSL Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute COSEM-OSL with Huber prior
                        if options.Huber && any(options.COSEM_OSL)
                            if verbose
                                tStart = tic;
                            end
                            med = Huber_prior(im_vectors.Huber_COSEM_apu, options.weights_huber, Nx, Ny, Nz, Ndx, Ndy, Ndz, options.huber_delta);
                            if options.COSEM_OSL == 1
                                [im_vectors.Huber_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.Huber_COSEM_apu, D, options.beta_huber_cosem, med, epps, A, uu, ...
                                    C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            else
                                [im_vectors.Huber_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.Huber_COSEM_apu, D, options.beta_huber_cosem, med, epps, A, uu, ...
                                    C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-OSL Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute OSL with L-filter prior
                        if options.L && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            med = L_filter(im_vectors.L_OSL_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            im_vectors.L_OSL_apu = OSEM_im(im_vectors.L_OSL_apu, A, epps, uu, OSL(Summ, options.beta_L_osem, med, epps), SinD, is_transposed, options, ...
                                Nx, Ny, Nz, gaussK);
                            %                             im_vectors.L_OSL_apu = OSL_OSEM(im_vectors.L_OSL_apu, Summ, options.beta_L_osem, med, epps, A, uu, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute MBSREM with L-filter prior
                        if options.L && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            med = L_filter(im_vectors.L_MBSREM_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            im_vectors.L_MBSREM_apu = MBSREM(im_vectors.L_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                                iter, SinD, randoms_correction, is_transposed, options.beta_L_mbsrem, med, options, Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['MBSREM L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute BSREM with L-filter prior
                        if options.L && options.BSREM
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.ramla
                                im_vectors.L_BSREM_apu = im_vectors.RAMLA_apu;
                            else
                                im_vectors.L_BSREM_apu = BSREM_subiter(im_vectors.L_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            end
                            if any(im_vectors.L_BSREM_apu < 0)
                                warning('Negative values in BSREM, it is recommended to lower lambda value')
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['BSREM L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute ROSEM-MAP with L-filter prior
                        if options.L && options.ROSEM_MAP
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.rosem
                                im_vectors.L_ROSEM_apu = im_vectors.ROSEM_apu;
                            else
                                im_vectors.L_ROSEM_apu = ROSEM_subiter(im_vectors.L_ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed, options, ...
                                    Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute RBI-OSL with L-filter prior
                        if options.L && options.RBI_OSL
                            if verbose
                                tStart = tic;
                            end
                            med = L_filter(im_vectors.L_RBI_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            im_vectors.L_RBI_apu = RBI_subiter(im_vectors.L_RBI_apu, A, uu, epps, Summ, D, SinD, is_transposed, options.beta_L_rbi, med, options, ...
                                Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-OSL L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute COSEM-OSL with L-filter prior
                        if options.L && any(options.COSEM_OSL)
                            if verbose
                                tStart = tic;
                            end
                            med = L_filter(im_vectors.L_COSEM_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            if options.COSEM_OSL == 1
                                [im_vectors.L_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.L_COSEM_apu, D, options.beta_L_cosem, med, epps, A, uu, C_osl, ...
                                    options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            else
                                [im_vectors.L_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.L_COSEM_apu, D, options.beta_L_cosem, med, epps, A, uu, C_osl, 0, ...
                                    options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-OSL L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute OSL with FMH prior
                        if options.FMH && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            med = FMH(im_vectors.FMH_OSL_apu, options.tr_offsets, options.fmh_weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            im_vectors.FMH_OSL_apu = OSEM_im(im_vectors.FMH_OSL_apu, A, epps, uu, OSL(Summ, options.beta_fmh_osem, med, epps), SinD, is_transposed, options, ...
                                Nx, Ny, Nz, gaussK);
                            %                             im_vectors.FMH_OSL_apu = OSL_OSEM(im_vectors.FMH_OSL_apu, Summ, options.beta_fmh_osem, med, epps, A, uu, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.FMH && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            med = FMH(im_vectors.FMH_MBSREM_apu, options.tr_offsets, options.fmh_weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            im_vectors.FMH_MBSREM_apu = MBSREM(im_vectors.FMH_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                                iter, SinD, randoms_correction, is_transposed, options.beta_fmh_mbsrem, med, options, Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['MBSREM FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.FMH && options.BSREM
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.ramla
                                im_vectors.FMH_BSREM_apu = im_vectors.RAMLA_apu;
                            else
                                im_vectors.FMH_BSREM_apu = BSREM_subiter(im_vectors.FMH_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed, options, ...
                                    Nx, Ny, Nz, gaussK);
                            end
                            if any(im_vectors.FMH_BSREM_apu < 0)
                                warning('Negative values in BSREM, it is recommended to lower lambda value')
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['BSREM FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.FMH && options.ROSEM_MAP
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.rosem
                                im_vectors.FMH_ROSEM_apu = im_vectors.ROSEM_apu;
                            else
                                im_vectors.FMH_ROSEM_apu = ROSEM_subiter(im_vectors.FMH_ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed, ...
                                    options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.FMH && options.RBI_OSL
                            if verbose
                                tStart = tic;
                            end
                            med = FMH(im_vectors.FMH_RBI_apu, options.tr_offsets, options.fmh_weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            im_vectors.FMH_RBI_apu = RBI_subiter(im_vectors.FMH_RBI_apu, A, uu, epps, Summ, D, SinD, is_transposed, options.beta_fmh_rbi, med, options, ...
                                Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-OSL FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.FMH && any(options.COSEM_OSL)
                            if verbose
                                tStart = tic;
                            end
                            med = FMH(im_vectors.FMH_COSEM_apu, options.tr_offsets, options.fmh_weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            if options.COSEM_OSL == 1
                                [im_vectors.FMH_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.FMH_COSEM_apu, D, options.beta_fmh_cosem, med, epps, A, uu, ...
                                    C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            else
                                [im_vectors.FMH_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.FMH_COSEM_apu, D, options.beta_fmh_cosem, med, epps, A, uu, ...
                                    C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-OSL FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute OSL with weighted mean prior
                        if options.weighted_mean && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            med = Weighted_mean(im_vectors.Weighted_OSL_apu, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.med_no_norm);
                            im_vectors.Weighted_OSL_apu = OSEM_im(im_vectors.Weighted_OSL_apu, A, epps, uu, OSL(Summ, options.beta_weighted_osem, med, epps), SinD, ...
                                is_transposed, options, Nx, Ny, Nz, gaussK);
                            %                             im_vectors.Weighted_OSL_apu = OSL_OSEM(im_vectors.Weighted_OSL_apu, Summ, options.beta_weighted_osem, med, epps, A, uu, SinD, ...
                            %                                 is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.weighted_mean && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            med = Weighted_mean(im_vectors.Weighted_MBSREM_apu, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.med_no_norm);
                            im_vectors.Weighted_MBSREM_apu = MBSREM(im_vectors.Weighted_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, ...
                                options.lam_mbsrem, iter, SinD, randoms_correction, is_transposed, options.beta_weighted_mbsrem, med, options, Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['MBSREM weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.weighted_mean && options.BSREM
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.ramla
                                im_vectors.Weighted_BSREM_apu = im_vectors.RAMLA_apu;
                            else
                                im_vectors.Weighted_BSREM_apu = BSREM_subiter(im_vectors.Weighted_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed, ...
                                    options, Nx, Ny, Nz, gaussK);
                            end
                            if any(im_vectors.Weighted_BSREM_apu < 0)
                                warning('Negative values in BSREM, it is recommended to lower lambda value')
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['BSREM weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.weighted_mean && options.ROSEM_MAP
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.rosem
                                im_vectors.Weighted_ROSEM_apu = im_vectors.ROSEM_apu;
                            else
                                im_vectors.Weighted_ROSEM_apu = ROSEM_subiter(im_vectors.Weighted_ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, ...
                                    is_transposed, options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.weighted_mean && options.RBI_OSL
                            if verbose
                                tStart = tic;
                            end
                            med = Weighted_mean(im_vectors.Weighted_RBI_apu, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.med_no_norm);
                            im_vectors.Weighted_RBI_apu = RBI_subiter(im_vectors.Weighted_RBI_apu, A, uu, epps, Summ, ...
                                D, SinD, is_transposed, options.beta_weighted_rbi, med, options, Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-OSL weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.weighted_mean && any(options.COSEM_OSL)
                            if verbose
                                tStart = tic;
                            end
                            med = Weighted_mean(im_vectors.Weighted_COSEM_apu, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.med_no_norm);
                            if options.COSEM_OSL == 1
                                [im_vectors.Weighted_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.Weighted_COSEM_apu, D, options.beta_weighted_cosem, ...
                                    med, epps, A, uu, C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            else
                                [im_vectors.Weighted_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.Weighted_COSEM_apu, D, options.beta_weighted_cosem, ...
                                    med, epps, A, uu, C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-OSL weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute OSL with TV prior
                        if options.TV && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            grad = TVpriorFinal(im_vectors.TV_OSL_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, ...
                                options.tr_offsets);
                            im_vectors.TV_OSL_apu = OSEM_im(im_vectors.TV_OSL_apu, A, epps, uu, OSL(Summ, options.beta_TV_osem, grad, epps), SinD, is_transposed, ...
                                options, Nx, Ny, Nz, gaussK);
                            %                             im_vectors.TV_OSL_apu = OSL_OSEM(im_vectors.TV_OSL_apu, Summ, options.beta_TV_osem, grad, epps, A, uu, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.TV && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            grad = TVpriorFinal(im_vectors.TV_MBSREM_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, ...
                                options.tr_offsets);
                            im_vectors.TV_MBSREM_apu = MBSREM(im_vectors.TV_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                                iter, SinD, randoms_correction, is_transposed, options.beta_TV_mbsrem, grad, options, Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['MBSREM TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.TV && options.BSREM
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.ramla
                                im_vectors.TV_BSREM_apu = im_vectors.RAMLA_apu;
                            else
                                im_vectors.TV_BSREM_apu = BSREM_subiter(im_vectors.TV_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed, options, ...
                                    Nx, Ny, Nz, gaussK);
                            end
                            if any(im_vectors.TV_BSREM_apu < 0)
                                warning('Negative values in BSREM, it is recommended to lower lambda value')
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['BSREM TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.TV && options.ROSEM_MAP
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.rosem
                                im_vectors.TV_ROSEM_apu = im_vectors.ROSEM_apu;
                            else
                                im_vectors.TV_ROSEM_apu = ROSEM_subiter(im_vectors.TV_ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed, ...
                                    options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.TV && options.RBI_OSL
                            if verbose
                                tStart = tic;
                            end
                            grad = TVpriorFinal(im_vectors.TV_RBI_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, ...
                                options.tr_offsets);
                            im_vectors.TV_RBI_apu = RBI_subiter(im_vectors.TV_RBI_apu, A, uu, epps, Summ, D, SinD, is_transposed, options.beta_TV_rbi, grad, options, ...
                                Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-OSL TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.TV && any(options.COSEM_OSL)
                            if verbose
                                tStart = tic;
                            end
                            grad = TVpriorFinal(im_vectors.TV_COSEM_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, ...
                                options.tr_offsets);
                            if options.COSEM_OSL == 1
                                [im_vectors.TV_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.TV_COSEM_apu, D, options.beta_TV_cosem, grad, epps, A, uu, ...
                                    C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            else
                                [im_vectors.TV_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.TV_COSEM_apu, D, options.beta_TV_cosem, grad, epps, A, uu, ...
                                    C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-OSL TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute OSL with MRP-AD prior
                        if options.AD && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            if osa_iter > 1
                                med = AD(im_vectors.AD_OSL_apu, options.FluxType, Nx, Ny, Nz, options);
                                im_vectors.AD_OSL_apu = OSEM_im(im_vectors.AD_OSL_apu, A, epps, uu, OSL(Summ, options.beta_ad_osem, med, epps), SinD, is_transposed, ...
                                    options, Nx, Ny, Nz, gaussK);
                                %                                 im_vectors.AD_OSL_apu = OSL_OSEM(im_vectors.AD_OSL_apu, Summ, options.beta_ad_osem, med, epps, A, uu, SinD, is_transposed);
                            else
                                im_vectors.AD_OSL_apu = OSEM_im(im_vectors.AD_OSL_apu, A, epps, uu, Summ, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.AD && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            if osa_iter > 1
                                med = AD(im_vectors.AD_MBSREM_apu, options.FluxType, Nx, Ny, Nz, options);
                                im_vectors.AD_MBSREM_apu = MBSREM(im_vectors.AD_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                                    iter, SinD, randoms_correction, is_transposed, options.beta_ad_mbsrem, med, options, Nx, Ny, Nz, gaussK);
                            else
                                im_vectors.AD_MBSREM_apu = MBSREM(im_vectors.AD_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, ...
                                    SinD, randoms_correction, is_transposed, [], [], options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['MBSREM AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.AD && options.BSREM
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.ramla
                                im_vectors.AD_BSREM_apu = im_vectors.RAMLA_apu;
                            else
                                im_vectors.AD_BSREM_apu = BSREM_subiter(im_vectors.AD_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            end
                            if any(im_vectors.AD_BSREM_apu < 0)
                                warning('Negative values in BSREM, it is recommended to lower lambda value')
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['BSREM AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.AD && options.ROSEM_MAP
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.rosem
                                im_vectors.AD_ROSEM_apu = im_vectors.ROSEM_apu;
                            else
                                im_vectors.AD_ROSEM_apu = ROSEM_subiter(im_vectors.AD_ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed, options, ...
                                    Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.AD && options.RBI_OSL
                            if verbose
                                tStart = tic;
                            end
                            med = AD(im_vectors.AD_RBI_apu, options.FluxType, Nx, Ny, Nz, options);
                            im_vectors.AD_RBI_apu = RBI_subiter(im_vectors.AD_RBI_apu, A, uu, epps, Summ, D, SinD, is_transposed, options.beta_ad_rbi, med, options, ...
                                Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-OSL AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.AD && any(options.COSEM_OSL)
                            if verbose
                                tStart = tic;
                            end
                            med = AD(im_vectors.AD_COSEM_apu, options.FluxType, Nx, Ny, Nz, options);
                            if options.COSEM_OSL == 1
                                [im_vectors.AD_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.AD_COSEM_apu, D, options.beta_ad_cosem, med, epps, A, uu, ...
                                    C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            else
                                [im_vectors.AD_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.AD_COSEM_apu, D, options.beta_ad_cosem, med, epps, A, uu, ...
                                    C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-OSL AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute OSL with APLS prior
                        if options.APLS && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            grad = TVpriorFinal(im_vectors.APLS_OSL_apu, [], Nx, Ny, Nz, true, options, 5);
                            im_vectors.APLS_OSL_apu = OSEM_im(im_vectors.APLS_OSL_apu, A, epps, uu, OSL(Summ, options.beta_APLS_osem, grad, epps), SinD, is_transposed, options, ...
                                Nx, Ny, Nz, gaussK);
                            %                             im_vectors.APLS_OSL_apu = OSL_OSEM(im_vectors.APLS_OSL_apu, Summ, options.beta_APLS_osem, grad, epps, A, uu, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.APLS && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            grad = TVpriorFinal(im_vectors.APLS_MBSREM_apu, [], Nx, Ny, Nz, true, options, 5);
                            im_vectors.APLS_MBSREM_apu = MBSREM(im_vectors.APLS_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                                iter, SinD, randoms_correction, is_transposed, options.beta_APLS_mbsrem, grad, options, Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['MBSREM APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.APLS && options.BSREM
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.ramla
                                im_vectors.APLS_BSREM_apu = im_vectors.RAMLA_apu;
                            else
                                im_vectors.APLS_BSREM_apu = BSREM_subiter(im_vectors.APLS_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed, options, ...
                                    Nx, Ny, Nz, gaussK);
                            end
                            if any(im_vectors.APLS_BSREM_apu < 0)
                                warning('Negative values in BSREM, it is recommended to lower lambda value')
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['BSREM APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.APLS && options.ROSEM_MAP
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.rosem
                                im_vectors.APLS_ROSEM_apu = im_vectors.ROSEM_apu;
                            else
                                im_vectors.APLS_ROSEM_apu = ROSEM_subiter(im_vectors.APLS_ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed, ...
                                    options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.APLS && options.RBI_OSL
                            if verbose
                                tStart = tic;
                            end
                            grad = TVpriorFinal(im_vectors.APLS_RBI_apu, [], Nx, Ny, Nz, true, options, 5);
                            im_vectors.APLS_RBI_apu = RBI_subiter(im_vectors.APLS_RBI_apu, A, uu, epps, Summ, SinD, D, ...
                                is_transposed, options.beta_APLS_rbi, grad, options, Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-OSL APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.APLS && any(options.COSEM_OSL)
                            if verbose
                                tStart = tic;
                            end
                            grad = TVpriorFinal(im_vectors.APLS_COSEM_apu, [], Nx, Ny, Nz, true, options, 5);
                            if options.COSEM_OSL == 1
                                [im_vectors.APLS_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.APLS_COSEM_apu, D, options.beta_APLS_cosem, grad, A, uu, ...
                                    epps, C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            else
                                [im_vectors.APLS_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.APLS_COSEM_apu, D, options.beta_APLS_cosem, grad, A, uu, ...
                                    epps, C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-OSL APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute OSL with TGV prior
                        if options.TGV && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            grad = TGV(im_vectors.TGV_OSL_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                            im_vectors.TGV_OSL_apu = OSEM_im(im_vectors.TGV_OSL_apu, A, epps, uu, OSL(Summ, options.beta_TGV_osem, grad, epps), SinD, is_transposed, options, ...
                                Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.TGV && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            grad = TGV(im_vectors.TGV_MBSREM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                            im_vectors.TGV_MBSREM_apu = MBSREM(im_vectors.TGV_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                                iter, SinD, randoms_correction, is_transposed, options.beta_TGV_mbsrem, grad, options, Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['MBSREM TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.TGV && options.BSREM
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.ramla
                                im_vectors.TGV_BSREM_apu = im_vectors.RAMLA_apu;
                            else
                                im_vectors.TGV_BSREM_apu = BSREM_subiter(im_vectors.TGV_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed, options, ...
                                    Nx, Ny, Nz, gaussK);
                            end
                            if any(im_vectors.TGV_BSREM_apu < 0)
                                warning('Negative values in BSREM, it is recommended to lower lambda value')
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['BSREM TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.TGV && options.ROSEM_MAP
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.rosem
                                im_vectors.TGV_ROSEM_apu = im_vectors.ROSEM_apu;
                            else
                                im_vectors.TGV_ROSEM_apu = ROSEM_subiter(im_vectors.TGV_ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed, options, ...
                                    Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.TGV && options.RBI_OSL
                            if verbose
                                tStart = tic;
                            end
                            grad = TGV(im_vectors.TGV_RBI_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                            im_vectors.TGV_RBI_apu = RBI_subiter(im_vectors.TGV_RBI_apu, A, uu, epps, Summ, D, SinD, is_transposed, options.beta_TGV_rbi, grad, options, ...
                                Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-OSL TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.TGV && any(options.COSEM_OSL)
                            if verbose
                                tStart = tic;
                            end
                            grad = TGV(im_vectors.TGV_COSEM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                            if options.COSEM_OSL == 1
                                [im_vectors.TGV_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.TGV_COSEM_apu, D, options.beta_TGV_cosem, grad, A, uu, ...
                                    epps, C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            else
                                [im_vectors.TGV_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.TGV_COSEM_apu, D, options.beta_TGV_cosem, grad, A, uu, ...
                                    epps, C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-OSL TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute OSL with NLM prior
                        if options.NLM && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            med = NLM(im_vectors.NLM_OSL_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                options.sigma, epps, Nx, Ny, Nz, options);
                            im_vectors.NLM_OSL_apu = OSEM_im(im_vectors.NLM_OSL_apu, A, epps, uu, OSL(Summ, options.beta_NLM_osem, med, epps), SinD, is_transposed, options, ...
                                Nx, Ny, Nz, gaussK);
                            %                             im_vectors.NLM_OSL_apu = OSL_OSEM(im_vectors.NLM_OSL_apu, Summ, options.beta_NLM_osem, med, epps, A, uu, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.NLM && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            med = NLM(im_vectors.NLM_MBSREM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                options.sigma, epps, Nx, Ny, Nz, options);
                            im_vectors.NLM_MBSREM_apu = MBSREM(im_vectors.NLM_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                                iter, SinD, randoms_correction, is_transposed, options.beta_NLM_mbsrem, med, options, Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['MBSREM NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.NLM && options.BSREM
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.ramla
                                im_vectors.NLM_BSREM_apu = im_vectors.RAMLA_apu;
                            else
                                im_vectors.NLM_BSREM_apu = BSREM_subiter(im_vectors.NLM_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed, options, ...
                                    Nx, Ny, Nz, gaussK);
                            end
                            if any(im_vectors.NLM_BSREM_apu < 0)
                                warning('Negative values in BSREM, it is recommended to lower lambda value')
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['BSREM NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.NLM && options.ROSEM_MAP
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.rosem
                                im_vectors.NLM_ROSEM_apu = im_vectors.ROSEM_apu;
                            else
                                im_vectors.NLM_ROSEM_apu = ROSEM_subiter(im_vectors.NLM_ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed, ...
                                    options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.NLM && options.RBI_OSL
                            if verbose
                                tStart = tic;
                            end
                            med = NLM(im_vectors.NLM_RBI_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                options.sigma, epps, Nx, Ny, Nz, options);
                            im_vectors.NLM_RBI_apu = RBI_subiter(im_vectors.NLM_RBI_apu, A, uu, epps, Summ, D, SinD, is_transposed, options.beta_NLM_rbi, med, options, ...
                                Nx, Ny, Nz, gaussK);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-OSL NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.NLM && any(options.COSEM_OSL)
                            if verbose
                                tStart = tic;
                            end
                            med = NLM(im_vectors.NLM_COSEM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                options.sigma, epps, Nx, Ny, Nz, options);
                            if options.COSEM_OSL == 1
                                [im_vectors.NLM_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.NLM_COSEM_apu, D, options.beta_NLM_cosem, med, A, uu, ...
                                    epps, C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            else
                                [im_vectors.NLM_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.NLM_COSEM_apu, D, options.beta_NLM_cosem, med, A, uu, ...
                                    epps, C_osl, 0, options.COSEM_OSL, osa_iter, SinD, is_transposed, options, Nx, Ny, Nz, gaussK);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-OSL NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        clear A
                        fn = fieldnames(im_vectors);
                        for kk = 2 : 2 : numel(fn)
                            im_vectors.(fn{kk})(im_vectors.(fn{kk}) < 0) = epps;
                        end
                    end
                    if options.save_iter
                        iter_n = iter + 1;
                    else
                        iter_n = 1;
                    end
                    if options.osem
                        im_vectors.OSEM(:, iter_n) = im_vectors.OSEM_apu;
                    end
                    
                    if options.mramla
                        im_vectors.MRAMLA(:, iter_n) = im_vectors.MRAMLA_apu;
                    end
                    
                    if options.ramla
                        im_vectors.RAMLA(:, iter_n) = im_vectors.RAMLA_apu;
                    end
                    
                    if options.rosem
                        im_vectors.ROSEM(:, iter_n) = im_vectors.ROSEM_apu;
                    end
                    
                    if options.rbi
                        im_vectors.RBI(:, iter_n) = im_vectors.RBI_apu;
                    end
                    
                    if options.drama
                        im_vectors.DRAMA(:, iter_n) = im_vectors.DRAMA_apu;
                    end
                    
                    if options.cosem
                        im_vectors.COSEM(:, iter_n) = im_vectors.COSEM_apu;
                    end
                    
                    if options.ecosem
                        im_vectors.ECOSEM(:, iter_n) = im_vectors.ECOSEM_apu;
                    end
                    
                    if options.acosem
                        im_vectors.ACOSEM(:, iter_n) = im_vectors.ACOSEM_apu;
                    end
                    
                    if options.MRP && options.OSL_OSEM
                        im_vectors.MRP_OSL(:, iter_n) = im_vectors.MRP_OSL_apu;
                    end
                    if options.MRP && options.MBSREM
                        im_vectors.MRP_MBSREM(:, iter_n) = im_vectors.MRP_MBSREM_apu;
                    end
                    
                    if options.MRP && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        med = MRP(im_vectors.MRP_BSREM_apu, options.medx, options.medy, options.medz, Nx, Ny, Nz, epps, options.tr_offsets, options.med_no_norm);
                        im_vectors.MRP_BSREM(:, iter_n) = BSREM_iter(im_vectors.MRP_BSREM_apu, options.lam, iter, options.beta_mrp_bsrem, med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM MRP iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['BSREM MRP iteration ' num2str(iter) ' finished'])
                        end
                    end
                    
                    if options.MRP && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        med = MRP(im_vectors.MRP_ROSEM_apu, options.medx, options.medy, options.medz, Nx, Ny, Nz, epps, options.tr_offsets, options.med_no_norm);
                        im_vectors.MRP_ROSEM(:, iter_n) = BSREM_iter(im_vectors.MRP_ROSEM_apu, options.lam_rosem, iter, options.beta_mrp_rosem, med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM MRP iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM MRP iteration ' num2str(iter) ' finished'])
                        end
                    end
                    
                    if options.MRP && options.RBI_OSL
                        im_vectors.MRP_RBI(:, iter_n) = im_vectors.MRP_RBI_apu;
                    end
                    
                    if options.MRP && any(options.COSEM_OSL)
                        im_vectors.MRP_COSEM(:, iter_n) = im_vectors.MRP_COSEM_apu;
                    end
                    
                    if options.quad && options.OSL_OSEM
                        im_vectors.Quad_OSL(:, iter_n) = im_vectors.Quad_OSL_apu;
                    end
                    if options.quad && options.MBSREM
                        im_vectors.Quad_MBSREM(:, iter_n) = im_vectors.Quad_MBSREM_apu;
                    end
                    
                    if options.quad && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        med = Quadratic_prior(im_vectors.Quad_BSREM_apu, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                        im_vectors.Quad_BSREM(:, iter_n) = BSREM_iter(im_vectors.Quad_BSREM_apu, options.lam, iter, options.beta_quad_bsrem, med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM quadratic iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['BSREM quadratic iteration ' num2str(iter) ' finished'])
                        end
                    end
                    
                    if options.quad && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        med = Quadratic_prior(im_vectors.Quad_ROSEM_apu, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                        im_vectors.Quad_ROSEM(:, iter_n) = BSREM_iter(im_vectors.Quad_ROSEM_apu, options.lam_rosem, iter, options.beta_quad_rosem, med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM quadratic iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM quadratic iteration ' num2str(iter) ' finished'])
                        end
                    end
                    
                    if options.quad && options.RBI_OSL
                        im_vectors.Quad_RBI(:, iter_n) = im_vectors.Quad_RBI_apu;
                    end
                    
                    if options.quad && any(options.COSEM_OSL)
                        im_vectors.Quad_COSEM(:, iter_n) = im_vectors.Quad_COSEM_apu;
                    end
                    
                    if options.Huber && options.OSL_OSEM
                        im_vectors.Huber_OSL(:, iter_n) = im_vectors.Huber_OSL_apu;
                    end
                    if options.Huber && options.MBSREM
                        im_vectors.Huber_MBSREM(:, iter_n) = im_vectors.Huber_MBSREM_apu;
                    end
                    
                    if options.Huber && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        med = Huber_prior(im_vectors.Huber_BSREM_apu, options.weights_huber, Nx, Ny, Nz, Ndx, Ndy, Ndz, options.huber_delta);
                        im_vectors.Huber_BSREM(:, iter_n) = BSREM_iter(im_vectors.Huber_BSREM_apu, options.lam, iter, options.beta_huber_bsrem, med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM Huber iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['BSREM Huber iteration ' num2str(iter) ' finished'])
                        end
                    end
                    
                    if options.Huber && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        med = Huber_prior(im_vectors.Huber_ROSEM_apu, options.weights_huber, Nx, Ny, Nz, Ndx, Ndy, Ndz, options.huber_delta);
                        im_vectors.Huber_ROSEM(:, iter_n) = BSREM_iter(im_vectors.Huber_ROSEM_apu, options.lam_rosem, iter, options.beta_huber_rosem, med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM Huber iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM Huber iteration ' num2str(iter) ' finished'])
                        end
                    end
                    
                    if options.Huber && options.RBI_OSL
                        im_vectors.Huber_RBI(:, iter_n) = im_vectors.Huber_RBI_apu;
                    end
                    
                    if options.Huber && any(options.COSEM_OSL)
                        im_vectors.Huber_COSEM(:, iter_n) = im_vectors.Huber_COSEM_apu;
                    end
                    
                    if options.L && options.OSL_OSEM
                        im_vectors.L_OSL(:, iter_n) = im_vectors.L_OSL_apu;
                    end
                    if options.L && options.MBSREM
                        im_vectors.L_MBSREM(:, iter_n) = im_vectors.L_MBSREM_apu;
                    end
                    
                    if options.L && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        med = L_filter(im_vectors.L_BSREM_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                        im_vectors.L_BSREM(:, iter_n) = BSREM_iter(im_vectors.L_BSREM_apu, options.lam, iter, options.beta_L_bsrem, med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM L-filter iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['BSREM L-filter iteration ' num2str(iter) ' finished'])
                        end
                    end
                    
                    if options.L && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        med = L_filter(im_vectors.L_ROSEM_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                        im_vectors.L_ROSEM(:, iter_n) = BSREM_iter(im_vectors.L_ROSEM_apu, options.lam_rosem, iter, options.beta_L_rosem, med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM L-filter iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM L-filter iteration ' num2str(iter) ' finished'])
                        end
                    end
                    
                    if options.L && options.RBI_OSL
                        im_vectors.L_RBI(:, iter_n) = im_vectors.L_RBI_apu;
                    end
                    
                    if options.L && any(options.COSEM_OSL)
                        im_vectors.L_COSEM(:, iter_n) = im_vectors.L_COSEM_apu;
                    end
                    
                    if options.FMH && options.OSL_OSEM
                        im_vectors.FMH_OSL(:, iter_n) = im_vectors.FMH_OSL_apu;
                    end
                    if options.FMH && options.MBSREM
                        im_vectors.FMH_MBSREM(:, iter_n) = im_vectors.FMH_MBSREM_apu;
                    end
                    
                    if options.FMH && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        med = FMH(im_vectors.FMH_BSREM_apu, options.tr_offsets, options.fmh_weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                            options.med_no_norm);
                        im_vectors.FMH_BSREM(:, iter_n) = BSREM_iter(im_vectors.FMH_BSREM_apu, options.lam, iter, options.beta_fmh_bsrem, med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM FMH iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['BSREM FMH iteration ' num2str(iter) ' finished'])
                        end
                    end
                    
                    if options.FMH && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        med = FMH(im_vectors.FMH_ROSEM_apu, options.tr_offsets, options.fmh_weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                            options.med_no_norm);
                        im_vectors.FMH_ROSEM(:, iter_n) = BSREM_iter(im_vectors.FMH_ROSEM_apu, options.lam_rosem, iter, options.beta_fmh_rosem, med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM FMH iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM FMH iteration ' num2str(iter) ' finished'])
                        end
                    end
                    
                    if options.FMH && options.RBI_OSL
                        im_vectors.FMH_RBI(:, iter_n) = im_vectors.FMH_RBI_apu;
                    end
                    
                    if options.FMH && any(options.COSEM_OSL)
                        im_vectors.FMH_COSEM(:, iter_n) = im_vectors.FMH_COSEM_apu;
                    end
                    
                    if options.weighted_mean && options.OSL_OSEM
                        im_vectors.Weighted_OSL(:, iter_n) = im_vectors.Weighted_OSL_apu;
                    end
                    if options.weighted_mean && options.MBSREM
                        im_vectors.Weighted_MBSREM(:, iter_n) = im_vectors.Weighted_MBSREM_apu;
                    end
                    
                    if options.weighted_mean && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        med = Weighted_mean(im_vectors.Weighted_BSREM_apu, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                            options.mean_type, epps, options.med_no_norm);
                        im_vectors.Weighted_BSREM(:, iter_n) = BSREM_iter(im_vectors.Weighted_BSREM_apu, options.lam, iter, options.beta_weighted_bsrem, ...
                            med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM weighted mean iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['BSREM weighted mean iteration ' num2str(iter) ' finished'])
                        end
                    end
                    
                    if options.weighted_mean && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        med = Weighted_mean(im_vectors.Weighted_ROSEM_apu, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                            options.mean_type, epps, options.med_no_norm);
                        im_vectors.Weighted_ROSEM(:, iter_n) = BSREM_iter(im_vectors.Weighted_ROSEM_apu, options.lam_rosem, iter, options.beta_weighted_rosem, ...
                            med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM weighted mean iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM weighted mean iteration ' num2str(iter) ' finished'])
                        end
                    end
                    
                    if options.weighted_mean && options.RBI_OSL
                        im_vectors.Weighted_RBI(:, iter_n) = im_vectors.Weighted_RBI_apu;
                    end
                    
                    if options.weighted_mean && any(options.COSEM_OSL)
                        im_vectors.Weighted_COSEM(:, iter_n) = im_vectors.Weighted_COSEM_apu;
                    end
                    
                    if options.TV && options.OSL_OSEM
                        im_vectors.TV_OSL(:, iter_n) = im_vectors.TV_OSL_apu;
                    end
                    if options.TV && options.MBSREM
                        im_vectors.TV_MBSREM(:, iter_n) = im_vectors.TV_MBSREM_apu;
                    end
                    
                    if options.TV && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        grad = TVpriorFinal(im_vectors.TV_BSREM_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, options.tr_offsets);
                        im_vectors.TV_BSREM(:, iter_n) = BSREM_iter(im_vectors.TV_BSREM_apu, options.lam, iter, options.beta_TV_bsrem, grad, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM TV iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['BSREM TV iteration ' num2str(iter) ' finished'])
                        end
                    end
                    
                    if options.TV && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        grad = TVpriorFinal(im_vectors.TV_ROSEM_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, options.tr_offsets);
                        im_vectors.TV_ROSEM(:, iter_n) = BSREM_iter(im_vectors.TV_ROSEM_apu, options.lam_rosem, iter, options.beta_TV_rosem, grad, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM TV iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM TV iteration ' num2str(iter) ' finished'])
                        end
                    end
                    
                    if options.TV && options.RBI_OSL
                        im_vectors.TV_RBI(:, iter_n) = im_vectors.TV_RBI_apu;
                    end
                    
                    if options.TV && any(options.COSEM_OSL)
                        im_vectors.TV_COSEM(:, iter_n) = im_vectors.TV_COSEM_apu;
                    end
                    
                    if options.AD && options.OSL_OSEM
                        im_vectors.AD_OSL(:, iter_n) = im_vectors.AD_OSL_apu;
                    end
                    if options.AD && options.MBSREM
                        im_vectors.AD_MBSREM(:, iter_n) = im_vectors.AD_MBSREM_apu;
                    end
                    
                    if options.AD && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        med = AD(im_vectors.AD_BSREM_apu, options.FluxType, Nx, Ny, Nz, options);
                        im_vectors.AD_BSREM(:, iter_n) = BSREM_iter(im_vectors.AD_BSREM_apu, options.lam, iter, options.beta_ad_bsrem, med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM AD iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['BSREM AD iteration ' num2str(iter) ' finished'])
                        end
                    end
                    if options.AD && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        med = AD(im_vectors.AD_ROSEM_apu, options.FluxType, Nx, Ny, Nz, options);
                        im_vectors.AD_ROSEM(:, iter_n) = BSREM_iter(im_vectors.AD_ROSEM_apu, options.lam_rosem, iter, options.beta_ad_rosem, med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM AD iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM AD iteration ' num2str(iter) ' finished'])
                        end
                    end
                    if options.AD && options.RBI_OSL
                        im_vectors.AD_RBI(:, iter_n) = im_vectors.AD_RBI_apu;
                    end
                    if options.AD && any(options.COSEM_OSL)
                        im_vectors.AD_COSEM(:, iter_n) = im_vectors.AD_COSEM_apu;
                    end
                    if options.APLS && options.OSL_OSEM
                        im_vectors.APLS_OSL(:, iter_n) = im_vectors.APLS_OSL_apu;
                    end
                    if options.APLS && options.MBSREM
                        im_vectors.APLS_MBSREM(:, iter_n) = im_vectors.APLS_MBSREM_apu;
                    end
                    if options.APLS && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        grad = TVpriorFinal(im_vectors.APLS_BSREM_apu, 0, Nx, Ny, Nz, true, options, 5);
                        im_vectors.APLS_BSREM(:, iter_n) = BSREM_iter(im_vectors.APLS_BSREM_apu, options.lam, iter, options.beta_APLS_bsrem, grad, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM APLS iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['BSREM APLS iteration ' num2str(iter) ' finished'])
                        end
                    end
                    if options.APLS && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        grad = TVpriorFinal(im_vectors.APLS_ROSEM_apu, 0, Nx, Ny, Nz, true, options, 5);
                        im_vectors.APLS_ROSEM(:, iter_n) = BSREM_iter(im_vectors.APLS_ROSEM_apu, options.lam_rosem, iter, options.beta_APLS_rosem, grad, ...
                            epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM APLS iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM APLS iteration ' num2str(iter) ' finished'])
                        end
                    end
                    if options.APLS && options.RBI_OSL
                        im_vectors.APLS_RBI(:, iter_n) = im_vectors.APLS_RBI_apu;
                    end
                    if options.APLS && any(options.COSEM_OSL)
                        im_vectors.APLS_COSEM(:, iter_n) = im_vectors.APLS_COSEM_apu;
                    end
                    if options.TGV && options.OSL_OSEM
                        im_vectors.TGV_OSL(:, iter_n) = im_vectors.TGV_OSL_apu;
                    end
                    if options.TGV && options.MBSREM
                        im_vectors.TGV_MBSREM(:, iter_n) = im_vectors.TGV_MBSREM_apu;
                    end
                    if options.TGV && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        grad = TGV(im_vectors.TGV_BSREM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                        im_vectors.TGV_BSREM(:, iter_n) = BSREM_iter(im_vectors.TGV_BSREM_apu, options.lam, iter, options.beta_TGV_bsrem, grad, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM TGV iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['BSREM TGV iteration ' num2str(iter) ' finished'])
                        end
                    end
                    if options.TGV && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        grad = TGV(im_vectors.TGV_ROSEM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                        im_vectors.TGV_ROSEM(:, iter_n) = BSREM_iter(im_vectors.TGV_ROSEM_apu, options.lam_rosem, iter, options.beta_TGV_rosem, grad, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM TGV iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM TGV iteration ' num2str(iter) ' finished'])
                        end
                    end
                    if options.TGV && options.RBI_OSL
                        im_vectors.TGV_RBI(:, iter_n) = im_vectors.TGV_RBI_apu;
                    end
                    if options.TGV && any(options.COSEM_OSL)
                        im_vectors.TGV_COSEM(:, iter_n) = im_vectors.TGV_COSEM_apu;
                    end
                    if options.NLM && options.OSL_OSEM
                        im_vectors.NLM_OSL(:, iter_n) = im_vectors.NLM_OSL_apu;
                    end
                    if options.NLM && options.MBSREM
                        im_vectors.NLM_MBSREM(:, iter_n) = im_vectors.NLM_MBSREM_apu;
                    end
                    if options.NLM && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        med = NLM(im_vectors.NLM_BSREM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                            options.sigma, epps, Nx, Ny, Nz, options);
                        im_vectors.NLM_BSREM(:, iter_n) = BSREM_iter(im_vectors.NLM_BSREM_apu, options.lam_rosem, iter, options.beta_NLM_bsrem, med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['BSREM NLM iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['BSREM NLM iteration ' num2str(iter) ' finished'])
                        end
                    end
                    if options.NLM && options.ROSEM_MAP
                        if verbose
                            tStart = tic;
                        end
                        med = NLM(im_vectors.NLM_ROSEM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                            options.sigma, epps, Nx, Ny, Nz, options);
                        im_vectors.NLM_ROSEM(:, iter_n) = BSREM_iter(im_vectors.NLM_ROSEM_apu, options.lam_rosem, iter, options.beta_NLM_rosem, med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM NLM iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM NLM iteration ' num2str(iter) ' finished'])
                        end
                    end
                    if options.NLM && options.RBI_OSL
                        im_vectors.NLM_RBI(:, iter_n) = im_vectors.NLM_RBI_apu;
                    end
                    if options.NLM && any(options.COSEM_OSL)
                        im_vectors.NLM_COSEM(:, iter_n) = im_vectors.NLM_COSEM_apu;
                    end
                    if options.use_psf && options.deblurring
                        im_vectors = computeDeblur(im_vectors, options, iter, gaussK, Nx, Ny, Nz);
                    end
                    disp(['Iteration ' num2str(iter) ' finished'])
                end
                %% Implementation 5 (unsupported and untested)
            elseif options.implementation == 5
                if use_raw_data
                    xy_index = uint32(0);
                    z_index = uint32(0);
                else
                    if isempty(pseudot)
                        pseudot = uint32(0);
                    end
                    LL = uint16(0);
                end
                %                 if use_raw_data == false
                kernel_path = which('siddon_kernel.cl');
                [im_vectors.OSEM] = improved_Siddon_openCL( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, ...
                    single(NSlices), size_x, zmax, NSinos, options.vaimennus, pituus, uint32(attenuation_correction), uint32(Niter), uint32(subsets), rekot, ...
                    single(epps), single(full(Sino)), single(options.x0(:)), lor_a, summa, xy_index, z_index, LL, pseudot, det_per_ring, ...
                    uint8(use_raw_data), options.verbose, device);
                %% Implementation 4
                % Matrix-free method
                % Parallel
                % Only C++ code (no pure MATLAB implementation)
                % Does not support MRAMLA or MBSREM
                % Only one algorihtm/prior at a time
            elseif options.implementation == 4
                if llo == 1
                    no_norm = false;
                end
                if ~use_raw_data
                    if isempty(pseudot)
                        pseudot = uint32(0);
                    end
                end
                options.n_rays_transaxial = uint16(options.n_rays_transaxial);
                options.n_rays_axial = uint16(options.n_rays_axial);
                if options.rings > 1
                    dc_z = z_det(2,1) - z_det(1,1);
                else
                    dc_z = options.cr_pz;
                end
                if options.cosem || options.ecosem || options.acosem || options.RBI_OSL || options.rbi || any(options.COSEM_OSL)
                    if llo == 1
                        f_Summ = zeros(Nx*Ny*Nz,subsets);
                    end
                    D = zeros(Nx*Ny*Nz, 1);
                    if options.cosem || options.ecosem || options.COSEM_OSL == 2
                        options.h = 1;
                    end
                    if options.cosem || options.ecosem
                        C_co = zeros(Nx*Ny*Nz,subsets);
                    elseif options.acosem
                        C_aco = zeros(Nx*Ny*Nz,subsets);
                    elseif any(options.COSEM_OSL)
                        C_osl = zeros(Nx*Ny*Nz,subsets);
                    end
                    if options.ecosem
                        im_vectors.COSEM_apu = im_vectors.OSEM_apu;
                        im_vectors.OSEM_apu2 = im_vectors.OSEM_apu;
                    end
                    if options.use_psf
                        OSEM_apu = computeConvolution(im_vectors.OSEM_apu, options, Nx, Ny, Nz, gaussK);
                    else
                        OSEM_apu = im_vectors.OSEM_apu;
                    end
                    for osa_iter = 1 : subsets
                        pituusS = pituus;
                        if randoms_correction
                            if iscell(options.SinDelayed)
                                SinD = double(options.SinDelayed{llo}(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                            else
                                SinD = double(options.SinDelayed(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                            end
                            if issparse(SinD)
                                SinD = (full(SinD));
                            end
                            SinD = SinD(:);
                            if TOF
                                SinD = SinD / options.TOF_bins;
                            end
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
                        uu = double(full(Sino(pituusS(osa_iter)+1:pituusS(osa_iter + 1))));
                        if use_raw_data
                            if ~list_mode_format
                                L_input = LL(pituus(osa_iter) * 2 + 1 : pituus(osa_iter + 1) * 2);
                            else
                                L_input = LL;
                            end
                            xy_index_input = uint32(0);
                            z_index_input = uint32(0);
                            TOFSize = int64(size(L_input,1));
                        else
                            L_input = uint16(0);
                            xy_index_input = xy_index(pituus(osa_iter)+1:pituus(osa_iter + 1));
                            z_index_input = z_index(pituus(osa_iter)+1:pituus(osa_iter + 1));
                            TOFSize = int64(numel(xy_index_input));
                        end
                        if options.precompute_lor
                            lor_a_input = lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1));
                        else
                            lor_a_input = uint16(0);
                        end
                        
                        if options.projector_type == 1
                            if exist('OCTAVE_VERSION','builtin') == 0
                                [Summ, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                    norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                    options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                    TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                    (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, false, ...
                                    list_mode_format, options.n_rays_transaxial, options.n_rays_axial, dc_z);
                            elseif exist('OCTAVE_VERSION','builtin') == 5
                                [Summ, rhs] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                    norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                    options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                    TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                    (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, false, ...
                                    list_mode_format, options.n_rays_transaxial, options.n_rays_axial, dc_z);
                            end
                        elseif options.projector_type == 2
                            if exist('OCTAVE_VERSION','builtin') == 0
                                [Summ, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                    norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                    options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                    TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                    (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, false, ...
                                    list_mode_format, options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z);
                            elseif exist('OCTAVE_VERSION','builtin') == 5
                                [Summ, rhs] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                    norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                    options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                    TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                    (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, false, ...
                                    list_mode_format, options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z);
                            end
                        elseif options.projector_type == 3
                            if exist('OCTAVE_VERSION','builtin') == 0
                                [Summ, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                    norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                    options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                    TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                    (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, false, ...
                                    list_mode_format, x_center, y_center, z_center, bmin, bmax, Vmax, V);
                            elseif exist('OCTAVE_VERSION','builtin') == 5
                                [Summ, rhs] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                    norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                    options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                    TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                    (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, false, ...
                                    list_mode_format, x_center, y_center, z_center, bmin, bmax, Vmax, V);
                            end
                        else
                            error('Unsupported projector')
                        end
                        if options.cosem || options.ecosem
                            C_co(:, osa_iter) = im_vectors.OSEM_apu .* rhs;
                        elseif options.acosem
                            C_aco(:, osa_iter) = im_vectors.OSEM_apu.^(1/options.h) .* rhs;
                        elseif options.COSEM_OSL == 1
                            C_osl(:, osa_iter) = im_vectors.OSEM_apu.^(1/options.h) .* rhs;
                        elseif options.COSEM_OSL == 2
                            C_osl(:, osa_iter) = im_vectors.OSEM_apu .* rhs;
                        end
                        if llo == 1
                            D = D + Summ;
                            if options.use_psf
                                Summ = computeConvolution(Summ, options, Nx, Ny, Nz, gaussK);
                            end
                            Summ(Summ < epps) = epps;
                            f_Summ(:,osa_iter) = Summ;
                        end
                    end
                    if options.use_psf
                        D = computeConvolution(D, options, Nx, Ny, Nz, gaussK);
                    end
                    D(D < epps) = epps;
                end
                for iter = 1 : Niter
                    if OS_bool
                        if verbose
                            tStart_iter = tic;
                        end
                        if iter == 1 && llo == 1 && ~no_norm
                            f_Summ = ones(Nx*Ny*Nz,subsets);
                        end
                        for osa_iter = 1 : subsets
                            if randoms_correction
                                if iscell(options.SinDelayed)
                                    SinD = double(options.SinDelayed{llo}(pituus(osa_iter)+1:pituus(osa_iter + 1)));
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
                            if verbose
                                tStart = tic;
                            end
                            if use_raw_data
                                if ~list_mode_format
                                    L_input = LL(pituus(osa_iter) * 2 + 1 : pituus(osa_iter + 1) * 2);
                                else
                                    L_input = LL;
                                end
                                xy_index_input = uint32(0);
                                z_index_input = uint32(0);
                                TOFSize = int64(size(L_input,1));
                                fullSize = size(LL,1);
                            else
                                L_input = uint16(0);
                                xy_index_input = xy_index(pituus(osa_iter)+1:pituus(osa_iter + 1));
                                z_index_input = z_index(pituus(osa_iter)+1:pituus(osa_iter + 1));
                                TOFSize = int64(numel(xy_index_input));
                                fullSize = length(xy_index);
                            end
                            if options.precompute_lor
                                lor_a_input = lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1));
                            else
                                lor_a_input = uint16(0);
                            end
                            if TOF
                                uu = zeros(TOFSize * options.TOF_bins, 1);
                                for dd = 1 : options.TOF_bins
                                    uu(1 + TOFSize * (dd - 1) : TOFSize * dd) = double(full(Sino(pituus(osa_iter) + 1 + fullSize * (dd - 1) : pituus(osa_iter + 1) + fullSize * (dd - 1))));
                                end
                            else
                                uu = double(full(Sino(pituus(osa_iter)+1:pituus(osa_iter + 1))));
                            end
                            uu(isnan(uu)) = 0;
                            uu(isinf(uu)) = 0;
                            
                            if options.use_psf
                                OSEM_apu = computeConvolution(im_vectors.OSEM_apu, options, Nx, Ny, Nz, gaussK);
                            else
                                OSEM_apu = im_vectors.OSEM_apu;
                            end
                            if options.projector_type == 1
                                if exist('OCTAVE_VERSION','builtin') == 0
                                    [Summ, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                        norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                        options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                        TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                        (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, false, ...
                                        list_mode_format, options.n_rays_transaxial, options.n_rays_axial, dc_z);
                                elseif exist('OCTAVE_VERSION','builtin') == 5
                                    [Summ, rhs] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                        norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                        options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                        TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                        (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, false, ...
                                        list_mode_format, options.n_rays_transaxial, options.n_rays_axial, dc_z);
                                end
                            elseif options.projector_type == 2
                                if exist('OCTAVE_VERSION','builtin') == 0
                                    [Summ, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                        norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                        options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                        TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                        (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, false, ...
                                        list_mode_format, options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z);
                                elseif exist('OCTAVE_VERSION','builtin') == 5
                                    [Summ, rhs] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                        norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                        options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                        TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                        (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, false, ...
                                        list_mode_format, options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z);
                                end
                            elseif options.projector_type == 3
                                if exist('OCTAVE_VERSION','builtin') == 0
                                    [Summ, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                        norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                        options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                        TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                        (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, false, ...
                                        list_mode_format, x_center, y_center, z_center, bmin, bmax, Vmax, V);
                                elseif exist('OCTAVE_VERSION','builtin') == 5
                                    [Summ, rhs] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                        norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                        options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                        TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                        (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, false, ...
                                        list_mode_format, x_center, y_center, z_center, bmin, bmax, Vmax, V);
                                end
                            else
                                error('Unsupported projector')
                            end
                            
                            if options.use_psf
                                rhs = computeConvolution(rhs, options, Nx, Ny, Nz, gaussK);
                            end
                            rhs(rhs < epps) = epps;
                            if iter == 1 && llo == 1 && ~no_norm
                                if options.use_psf
                                    Summ = computeConvolution(Summ, options, Nx, Ny, Nz, gaussK);
                                end
                                Summ(Summ < epps) = epps;
                                f_Summ(:,osa_iter) = Summ;
                                if list_mode_format && sum(abs(Summ-rhs)) < 1e-2
                                    f_Summ(:,osa_iter) = max(f_Summ(:,osa_iter));
                                end
                            end
                            if options.osem
                                im_vectors.OSEM_apu = OSEM_im(im_vectors.OSEM_apu, rhs, f_Summ(:,osa_iter));
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.ramla
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, f_Summ(:,osa_iter), rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    warning(['Negative values in RAMLA, lower lambda value! lambda <= ' num2str(min(1./f_Summ(:,osa_iter)))])
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['RAMLA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.mramla
                                error('MRAMLA is not supported when using implementation 4')
                            elseif options.rosem
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, options.lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute RBI
                            elseif options.rbi
                                im_vectors.OSEM_apu = RBI_subiter(im_vectors.OSEM_apu, f_Summ(:,osa_iter), rhs, [], [], D, [], []);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['RBI sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute DRAMA
                            elseif options.drama
                                im_vectors.OSEM_apu = DRAMA_subiter(im_vectors.OSEM_apu, options.lam_drama, epps, iter, f_Summ(:,osa_iter), osa_iter, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['DRAMA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute COSEM
                            elseif options.cosem
                                [im_vectors.OSEM_apu, C_co] = COSEM_im(im_vectors.OSEM_apu, rhs, C_co, D, osa_iter, [], [], [], []);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['COSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute ECOSEM
                            elseif options.ecosem
                                no_norm_ecosem = true;
                                if options.use_psf
                                    OSEM_apu = computeConvolution(im_vectors.COSEM_apu, options, Nx, Ny, Nz, gaussK);
                                else
                                    OSEM_apu = im_vectors.COSEM_apu;
                                end
                                if options.projector_type == 1
                                    if exist('OCTAVE_VERSION','builtin') == 0
                                        [~, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                            norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                            options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                            TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                            (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm_ecosem, options.precompute_lor, false, ...
                                            list_mode_format, options.n_rays_transaxial, options.n_rays_axial, dc_z);
                                    elseif exist('OCTAVE_VERSION','builtin') == 5
                                        [~, rhs] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                            norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                            options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                            TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                            (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm_ecosem, options.precompute_lor, false, ...
                                            list_mode_format, options.n_rays_transaxial, options.n_rays_axial, dc_z);
                                    end
                                elseif options.projector_type == 2
                                    if exist('OCTAVE_VERSION','builtin') == 0
                                        [~, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                            norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                            options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                            TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                            (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm_ecosem, options.precompute_lor, false, ...
                                            list_mode_format, options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, dec);
                                    elseif exist('OCTAVE_VERSION','builtin') == 5
                                        [~, rhs] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                            norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                            options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                            TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                            (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm_ecosem, options.precompute_lor, false, ...
                                            list_mode_format, options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, dec);
                                    end
                                elseif options.projector_type == 3
                                    if exist('OCTAVE_VERSION','builtin') == 0
                                        [~, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                            norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                            options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                            TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                            (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm_ecosem, options.precompute_lor, false, ...
                                            list_mode_format, x_center, y_center, z_center, dec, bmin, bmax, Vmax, V);
                                    elseif exist('OCTAVE_VERSION','builtin') == 5
                                        [~, rhs] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                            norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                            options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                            TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                            (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm_ecosem, options.precompute_lor, false, ...
                                            list_mode_format, x_center, y_center, z_center, dec, bmin, bmax, Vmax, V);
                                    end
                                else
                                    error('Unsupported projector')
                                end
                                [im_vectors.COSEM_apu, C_co] = COSEM_im(im_vectors.COSEM_apu, rhs, C_co, D, osa_iter, [], [], [], []);
                                
                                if options.use_psf
                                    OSEM_apu = computeConvolution(im_vectors.OSEM_apu2, options, Nx, Ny, Nz, gaussK);
                                else
                                    OSEM_apu = im_vectors.OSEM_apu2;
                                end
                                if options.projector_type == 1
                                    if exist('OCTAVE_VERSION','builtin') == 0
                                        [~, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                            norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                            options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                            TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                            (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm_ecosem, options.precompute_lor, false, ...
                                            list_mode_format, options.n_rays_transaxial, options.n_rays_axial, dc_z);
                                    elseif exist('OCTAVE_VERSION','builtin') == 5
                                        [~, rhs] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                            norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                            options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                            TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                            (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm_ecosem, options.precompute_lor, false, ...
                                            list_mode_format, options.n_rays_transaxial, options.n_rays_axial, dc_z);
                                    end
                                elseif options.projector_type == 2
                                    if exist('OCTAVE_VERSION','builtin') == 0
                                        [~, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                            norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                            options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                            TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                            (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm_ecosem, options.precompute_lor, false, ...
                                            list_mode_format, options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, dec);
                                    elseif exist('OCTAVE_VERSION','builtin') == 5
                                        [~, rhs] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                            norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                            options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                            TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                            (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm_ecosem, options.precompute_lor, false, ...
                                            list_mode_format, options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, dec);
                                    end
                                elseif options.projector_type == 3
                                    if exist('OCTAVE_VERSION','builtin') == 0
                                        [~, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                            norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                            options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                            TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                            (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm_ecosem, options.precompute_lor, false, ...
                                            list_mode_format, x_center, y_center, z_center, dec, bmin, bmax, Vmax, V);
                                    elseif exist('OCTAVE_VERSION','builtin') == 5
                                        [~, rhs] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                            norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                            options.scatter, scatter_input, options.global_correction_factor, lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                            TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                            (use_raw_data), uint32(1), epps, uu, OSEM_apu, uint32(options.projector_type), no_norm_ecosem, options.precompute_lor, false, ...
                                            list_mode_format, x_center, y_center, z_center, dec, bmin, bmax, Vmax, V);
                                    end
                                else
                                    error('Unsupported projector')
                                end
                                im_vectors.OSEM_apu2 = OSEM_im(im_vectors.OSEM_apu2, rhs, f_Summ(:,osa_iter));
                                im_vectors.OSEM_apu = ECOSEM_im(im_vectors.OSEM_apu, epps, D, im_vectors.COSEM_apu, im_vectors.OSEM_apu2);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ECOSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute ACOSEM
                            elseif options.acosem
                                [im_vectors.OSEM_apu, C_aco] = ACOSEM_im(im_vectors.OSEM_apu, rhs, C_aco, D, osa_iter, options.h, [], [], [], []);
                                im_vectors = ACOSEM_coefficient(im_vectors, options, Nx, Ny, Nz, dx, dz, dy, bx, by, bz, z_det, x, y, yy, xx, NSinos, NSlices, size_x, ...
                                    zmax, norm_input, SinD, pituus, osa_iter, attenuation_correction, normalization_correction, randoms_correction, options.scatter, scatter_input, lor_a_input, ...
                                    xy_index_input, z_index_input, L_input, pseudot, det_per_ring, use_raw_data, epps, uu, dc_z, x_center, y_center, z_center, dec, bmin, ...
                                    bmax, Vmax, V, gaussK, TOF, TOFSize, sigma_x, TOFCenter, list_mode_format);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ACOSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.MRP && options.OSL_OSEM
                                med = MRP(im_vectors.OSEM_apu, options.medx, options.medy, options.medz, Nx, Ny, Nz, epps, options.tr_offsets, options.med_no_norm);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_mrp_osem, med, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSEM-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.MRP && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    warning('Negative values in BSREM, lower lambda value!')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.MRP && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, options.lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM-MAP MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute RBI-OSL with MRP
                            elseif options.MRP && options.RBI_OSL
                                if verbose
                                    tStart = tic;
                                end
                                med = MRP(im_vectors.OSEM_apu, options.medx, options.medy, options.medz, Nx, Ny, Nz, epps, options.tr_offsets, options.med_no_norm);
                                im_vectors.OSEM_apu = RBI_subiter(im_vectors.OSEM_apu, f_Summ(:,osa_iter), rhs, [], [], D, [], [], options.beta_mrp_rbi, med);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['RBI-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute COSEM-OSL with MRP
                            elseif options.MRP && any(options.COSEM_OSL)
                                if verbose
                                    tStart = tic;
                                end
                                med = MRP(im_vectors.OSEM_apu, options.medx, options.medy, options.medz, Nx, Ny, Nz, epps, options.tr_offsets, options.med_no_norm);
                                [im_vectors.OSEM_apu, C_osl] = COSEM_OSL(im_vectors.OSEM_apu, D, options.beta_mrp_cosem, med, rhs, osa_iter, options.h, ...
                                    C_osl, options.COSEM_OSL, [], [], [], []);
                                if options.COSEM_OSL == 1
                                    im_vectors = ACOSEM_coefficient(im_vectors, options, Nx, Ny, Nz, dx, dz, dy, bx, by, bz, z_det, x, y, yy, xx, NSinos, NSlices, size_x, ...
                                        zmax, norm_input, SinD, pituus, osa_iter, attenuation_correction, normalization_correction, randoms_correction, options.scatter, scatter_input, lor_a_input, ...
                                        xy_index_input, z_index_input, L_input, pseudot, det_per_ring, use_raw_data, epps, uu, dc_z, x_center, y_center, z_center, dec, bmin, ...
                                        bmax, Vmax, V, gaussK, TOF, TOFSize, sigma_x, TOFCenter, list_mode_format);
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['COSEM-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.quad && options.OSL_OSEM
                                med = Quadratic_prior(im_vectors.OSEM_apu, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_quad_osem, med, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSEM-OSL Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.quad && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    warning('Negative values in BSREM, it is recommended to lower lambda value')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.quad && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, options.lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM-MAP Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute RBI-OSL with quad
                            elseif options.quad && options.RBI_OSL
                                if verbose
                                    tStart = tic;
                                end
                                med = Quadratic_prior(im_vectors.OSEM_apu, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                                im_vectors.OSEM_apu = RBI_subiter(im_vectors.OSEM_apu, f_Summ(:,osa_iter), rhs, [], [], D, [], [], options.beta_quad_rbi, med);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['RBI-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute COSEM-OSL with quad
                            elseif options.quad && any(options.COSEM_OSL)
                                if verbose
                                    tStart = tic;
                                end
                                med = Quadratic_prior(im_vectors.OSEM_apu, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                                [im_vectors.OSEM_apu, C_osl] = COSEM_OSL(im_vectors.OSEM_apu, D, options.beta_quad_cosem, med, rhs, osa_iter, options.h, ...
                                    C_osl, options.COSEM_OSL, [], [], [], []);
                                if options.COSEM_OSL == 1
                                    im_vectors = ACOSEM_coefficient(im_vectors, options, Nx, Ny, Nz, dx, dz, dy, bx, by, bz, z_det, x, y, yy, xx, NSinos, NSlices, size_x, ...
                                        zmax, norm_input, SinD, pituus, osa_iter, attenuation_correction, normalization_correction, randoms_correction, options.scatter, scatter_input, lor_a_input, ...
                                        xy_index_input, z_index_input, L_input, pseudot, det_per_ring, use_raw_data, epps, uu, dc_z, x_center, y_center, z_center, dec, bmin, ...
                                        bmax, Vmax, V, gaussK, TOF, TOFSize, sigma_x, TOFCenter, list_mode_format);
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['COSEM-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.Huber && options.OSL_OSEM
                                med = Huber_prior(im_vectors.OSEM_apu, options.weights_huber, Nx, Ny, Nz, Ndx, Ndy, Ndz, options.huber_delta);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_huber_osem, med, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSEM-OSL Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.Huber && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    warning('Negative values in BSREM, it is recommended to lower lambda value')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.Huber && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, options.lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM-MAP Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute RBI-OSL with Huber
                            elseif options.Huber && options.RBI_OSL
                                if verbose
                                    tStart = tic;
                                end
                                med = Huber_prior(im_vectors.OSEM_apu, options.weights_huber, Nx, Ny, Nz, Ndx, Ndy, Ndz, options.huber_delta);
                                im_vectors.OSEM_apu = RBI_subiter(im_vectors.OSEM_apu, f_Summ(:,osa_iter), rhs, [], [], D, [], [], options.beta_huber_rbi, med);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['RBI-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute COSEM-OSL with Huber
                            elseif options.Huber && any(options.COSEM_OSL)
                                if verbose
                                    tStart = tic;
                                end
                                med = Huber_prior(im_vectors.OSEM_apu, options.weights_huber, Nx, Ny, Nz, Ndx, Ndy, Ndz, options.huber_delta);
                                [im_vectors.OSEM_apu, C_osl] = COSEM_OSL(im_vectors.OSEM_apu, D, options.beta_huber_cosem, med, rhs, osa_iter, options.h, ...
                                    C_osl, options.COSEM_OSL, [], [], [], []);
                                if options.COSEM_OSL == 1
                                    im_vectors = ACOSEM_coefficient(im_vectors, options, Nx, Ny, Nz, dx, dz, dy, bx, by, bz, z_det, x, y, yy, xx, NSinos, NSlices, size_x, ...
                                        zmax, norm_input, SinD, pituus, osa_iter, attenuation_correction, normalization_correction, randoms_correction, options.scatter, scatter_input, lor_a_input, ...
                                        xy_index_input, z_index_input, L_input, pseudot, det_per_ring, use_raw_data, epps, uu, dc_z, x_center, y_center, z_center, dec, bmin, ...
                                        bmax, Vmax, V, gaussK, TOF, TOFSize, sigma_x, TOFCenter, list_mode_format);
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['COSEM-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.L && options.OSL_OSEM
                                med = L_filter(im_vectors.OSEM_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_L_osem, med, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSEM-OSL L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.L && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    warning('Negative values in BSREM, it is recommended to lower lambda value')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.L && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, options.lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM-MAP L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute RBI-OSL with L
                            elseif options.L && options.RBI_OSL
                                if verbose
                                    tStart = tic;
                                end
                                med = L_filter(im_vectors.OSEM_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                                im_vectors.OSEM_apu = RBI_subiter(im_vectors.OSEM_apu, f_Summ(:,osa_iter), rhs, [], [], D, [], [], options.beta_L_rbi, med);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['RBI-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute COSEM-OSL with L
                            elseif options.L && any(options.COSEM_OSL)
                                if verbose
                                    tStart = tic;
                                end
                                med = L_filter(im_vectors.OSEM_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                                [im_vectors.OSEM_apu, C_osl] = COSEM_OSL(im_vectors.OSEM_apu, D, options.beta_L_cosem, med, rhs, osa_iter, options.h, ...
                                    C_osl, options.COSEM_OSL, [], [], [], []);
                                if options.COSEM_OSL == 1
                                    im_vectors = ACOSEM_coefficient(im_vectors, options, Nx, Ny, Nz, dx, dz, dy, bx, by, bz, z_det, x, y, yy, xx, NSinos, NSlices, size_x, ...
                                        zmax, norm_input, SinD, pituus, osa_iter, attenuation_correction, normalization_correction, randoms_correction, options.scatter, scatter_input, lor_a_input, ...
                                        xy_index_input, z_index_input, L_input, pseudot, det_per_ring, use_raw_data, epps, uu, dc_z, x_center, y_center, z_center, dec, bmin, ...
                                        bmax, Vmax, V, gaussK, TOF, TOFSize, sigma_x, TOFCenter, list_mode_format);
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['COSEM-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.FMH && options.OSL_OSEM
                                med = FMH(im_vectors.OSEM_apu, options.tr_offsets, options.fmh_weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                    options.med_no_norm);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_fmh_osem, med, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSEM-OSL FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.FMH && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    warning('Negative values in BSREM, it is recommended to lower lambda value')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.FMH && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, options.lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM-MAP FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute RBI-OSL with FMH
                            elseif options.FMH && options.RBI_OSL
                                if verbose
                                    tStart = tic;
                                end
                                med = FMH(im_vectors.OSEM_apu, options.tr_offsets, options.fmh_weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                    options.med_no_norm);
                                im_vectors.OSEM_apu = RBI_subiter(im_vectors.OSEM_apu, f_Summ(:,osa_iter), rhs, [], [], D, [], [], options.beta_fmh_rbi, med);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['RBI-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute COSEM-OSL with FMH
                            elseif options.FMH && any(options.COSEM_OSL)
                                if verbose
                                    tStart = tic;
                                end
                                med = FMH(im_vectors.OSEM_apu, options.tr_offsets, options.fmh_weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                    options.med_no_norm);
                                [im_vectors.OSEM_apu, C_osl] = COSEM_OSL(im_vectors.OSEM_apu, D, options.beta_fmh_cosem, med, rhs, osa_iter, options.h, ...
                                    C_osl, options.COSEM_OSL, [], [], [], []);
                                if options.COSEM_OSL == 1
                                    im_vectors = ACOSEM_coefficient(im_vectors, options, Nx, Ny, Nz, dx, dz, dy, bx, by, bz, z_det, x, y, yy, xx, NSinos, NSlices, size_x, ...
                                        zmax, norm_input, SinD, pituus, osa_iter, attenuation_correction, normalization_correction, randoms_correction, options.scatter, scatter_input, lor_a_input, ...
                                        xy_index_input, z_index_input, L_input, pseudot, det_per_ring, use_raw_data, epps, uu, dc_z, x_center, y_center, z_center, dec, bmin, ...
                                        bmax, Vmax, V, gaussK, TOF, TOFSize, sigma_x, TOFCenter, list_mode_format);
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['COSEM-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.weighted_mean && options.OSL_OSEM
                                med = Weighted_mean(im_vectors.OSEM_apu, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                    options.mean_type, epps, options.med_no_norm);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_weighted_osem, med, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSEM-OSL weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.weighted_mean && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    warning('Negative values in BSREM, it is recommended to lower lambda value')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.weighted_mean && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, options.lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM-MAP weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute RBI-OSL with weighted mean
                            elseif options.weighted_mean && options.RBI_OSL
                                if verbose
                                    tStart = tic;
                                end
                                med = Weighted_mean(im_vectors.OSEM_apu, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                    options.mean_type, epps, options.med_no_norm);
                                im_vectors.OSEM_apu = RBI_subiter(im_vectors.OSEM_apu, f_Summ(:,osa_iter), rhs, [], [], D, [], [], options.beta_weighted_rbi, med);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['RBI-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute COSEM-OSL with weighted mean
                            elseif options.weighted_mean && any(options.COSEM_OSL)
                                if verbose
                                    tStart = tic;
                                end
                                med = Weighted_mean(im_vectors.OSEM_apu, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                    options.mean_type, epps, options.med_no_norm);
                                [im_vectors.OSEM_apu, C_osl] = COSEM_OSL(im_vectors.OSEM_apu, D, options.beta_weighted_cosem, med, rhs, osa_iter, options.h, ...
                                    C_osl, options.COSEM_OSL, [], [], [], []);
                                if options.COSEM_OSL == 1
                                    im_vectors = ACOSEM_coefficient(im_vectors, options, Nx, Ny, Nz, dx, dz, dy, bx, by, bz, z_det, x, y, yy, xx, NSinos, NSlices, size_x, ...
                                        zmax, norm_input, SinD, pituus, osa_iter, attenuation_correction, normalization_correction, randoms_correction, options.scatter, scatter_input, lor_a_input, ...
                                        xy_index_input, z_index_input, L_input, pseudot, det_per_ring, use_raw_data, epps, uu, dc_z, x_center, y_center, z_center, dec, bmin, ...
                                        bmax, Vmax, V, gaussK, TOF, TOFSize, sigma_x, TOFCenter, list_mode_format);
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['COSEM-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.TV && options.OSL_OSEM
                                grad = TVpriorFinal(im_vectors.OSEM_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, ...
                                    options.tr_offsets);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_TV_osem, grad, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSEM-OSL TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.TV && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    warning('Negative values in BSREM, it is recommended to lower lambda value')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.TV && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, options.lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM-MAP TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute RBI-OSL with TV
                            elseif options.TV && options.RBI_OSL
                                if verbose
                                    tStart = tic;
                                end
                                grad = TVpriorFinal(im_vectors.OSEM_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, ...
                                    options.tr_offsets);
                                im_vectors.OSEM_apu = RBI_subiter(im_vectors.OSEM_apu, f_Summ(:,osa_iter), rhs, [], [], D, [], [], options.beta_TV_rbi, grad);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['RBI-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute COSEM-OSL with TV
                            elseif options.TV && any(options.COSEM_OSL)
                                if verbose
                                    tStart = tic;
                                end
                                grad = TVpriorFinal(im_vectors.OSEM_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, ...
                                    options.tr_offsets);
                                [im_vectors.OSEM_apu, C_osl] = COSEM_OSL(im_vectors.OSEM_apu, D, options.beta_TV_cosem, grad, rhs, osa_iter, options.h, ...
                                    C_osl, options.COSEM_OSL, [], [], [], []);
                                if options.COSEM_OSL == 1
                                    im_vectors = ACOSEM_coefficient(im_vectors, options, Nx, Ny, Nz, dx, dz, dy, bx, by, bz, z_det, x, y, yy, xx, NSinos, NSlices, size_x, ...
                                        zmax, norm_input, SinD, pituus, osa_iter, attenuation_correction, normalization_correction, randoms_correction, options.scatter, scatter_input, lor_a_input, ...
                                        xy_index_input, z_index_input, L_input, pseudot, det_per_ring, use_raw_data, epps, uu, dc_z, x_center, y_center, z_center, dec, bmin, ...
                                        bmax, Vmax, V, gaussK, TOF, TOFSize, sigma_x, TOFCenter, list_mode_format);
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['COSEM-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.AD && options.OSL_OSEM
                                if osa_iter > 1
                                    med = AD(im_vectors.OSEM_apu, options.FluxType, Nx, Ny, Nz, options);
                                    im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_ad_osem, med, epps, rhs);
                                else
                                    im_vectors.OSEM_apu = OSEM_im(im_vectors.OSEM_apu, rhs, f_Summ(:,osa_iter), epps);
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSEM-OSL AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.AD && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    warning('Negative values in BSREM, it is recommended to lower lambda value')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.AD && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, options.lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM-MAP AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute RBI-OSL with AD
                            elseif options.AD && options.RBI_OSL
                                if verbose
                                    tStart = tic;
                                end
                                if osa_iter > 1
                                    med = AD(im_vectors.OSEM_apu, options.FluxType, Nx, Ny, Nz, options);
                                    im_vectors.OSEM_apu = RBI_subiter(im_vectors.OSEM_apu, f_Summ(:,osa_iter), rhs, [], [], D, [], [], options.beta_ad_rbi, med);
                                else
                                    im_vectors.OSEM_apu = RBI_subiter(im_vectors.OSEM_apu, f_Summ(:,osa_iter), rhs, [], [], D, [], []);
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['RBI-OSL AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute COSEM-OSL with AD
                            elseif options.AD && any(options.COSEM_OSL)
                                if verbose
                                    tStart = tic;
                                end
                                if osa_iter > 1
                                    med = AD(im_vectors.OSEM_apu, options.FluxType, Nx, Ny, Nz, options);
                                    [im_vectors.OSEM_apu, C_osl] = COSEM_OSL(im_vectors.OSEM_apu, D, options.beta_mrp_cosem, med, rhs, osa_iter, options.h, ...
                                        C_osl, options.COSEM_OSL, [], [], [], []);
                                    if options.COSEM_OSL == 1
                                        im_vectors = ACOSEM_coefficient(im_vectors, options, Nx, Ny, Nz, dx, dz, dy, bx, by, bz, z_det, x, y, yy, xx, NSinos, NSlices, size_x, ...
                                            zmax, norm_input, SinD, pituus, osa_iter, attenuation_correction, normalization_correction, randoms_correction, options.scatter, scatter_input, lor_a_input, ...
                                            xy_index_input, z_index_input, L_input, pseudot, det_per_ring, use_raw_data, epps, uu, dc_z, x_center, y_center, z_center, dec, bmin, ...
                                            bmax, Vmax, V, gaussK, TOF, TOFSize, sigma_x, TOFCenter, list_mode_format);
                                    end
                                else
                                    if options.COSEM_OSL == 1
                                        [im_vectors.OSEM_apu, C_osl] = ACOSEM_im(im_vectors.OSEM_apu, rhs, C_osl, D, osa_iter, options.h, [], [], [], []);
                                        im_vectors = ACOSEM_coefficient(im_vectors, options, Nx, Ny, Nz, dx, dz, dy, bx, by, bz, z_det, x, y, yy, xx, NSinos, NSlices, size_x, ...
                                            zmax, norm_input, SinD, pituus, osa_iter, attenuation_correction, normalization_correction, randoms_correction, options.scatter, scatter_input, lor_a_input, ...
                                            xy_index_input, z_index_input, L_input, pseudot, det_per_ring, use_raw_data, epps, uu, dc_z, x_center, y_center, z_center, dec, bmin, ...
                                            bmax, Vmax, V, gaussK, TOF, TOFSize, sigma_x, TOFCenter, list_mode_format);
                                    else
                                        [im_vectors.OSEM_apu, C_osl] = COSEM_im(im_vectors.OSEM_apu, rhs, C_osl, D, osa_iter, [], [], [], []);
                                    end
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['COSEM-OSL AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.APLS && options.OSL_OSEM
                                grad = TVpriorFinal(im_vectors.OSEM_apu, [], Nx, Ny, Nz, true, options, 5);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_APLS_osem, grad, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSEM-OSL APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.APLS && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    warning('Negative values in BSREM, it is recommended to lower lambda value')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.APLS && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, options.lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM-MAP APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute RBI-OSL with APLS
                            elseif options.APLS && options.RBI_OSL
                                if verbose
                                    tStart = tic;
                                end
                                grad = TVpriorFinal(im_vectors.OSEM_apu, [], Nx, Ny, Nz, true, options, 5);
                                im_vectors.OSEM_apu = RBI_subiter(im_vectors.OSEM_apu, f_Summ(:,osa_iter), rhs, [], [], D, [], [], options.beta_APLS_rbi, grad);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['RBI-OSL APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute COSEM-OSL with APLS
                            elseif options.APLS && any(options.COSEM_OSL)
                                if verbose
                                    tStart = tic;
                                end
                                grad = TVpriorFinal(im_vectors.OSEM_apu, [], Nx, Ny, Nz, true, options, 5);
                                [im_vectors.OSEM_apu, C_osl] = COSEM_OSL(im_vectors.OSEM_apu, D, options.beta_APLS_cosem, grad, rhs, osa_iter, options.h, ...
                                    C_osl, options.COSEM_OSL, [], [], [], []);
                                if options.COSEM_OSL == 1
                                    im_vectors = ACOSEM_coefficient(im_vectors, options, Nx, Ny, Nz, dx, dz, dy, bx, by, bz, z_det, x, y, yy, xx, NSinos, NSlices, size_x, ...
                                        zmax, norm_input, SinD, pituus, osa_iter, attenuation_correction, normalization_correction, randoms_correction, options.scatter, scatter_input, lor_a_input, ...
                                        xy_index_input, z_index_input, L_input, pseudot, det_per_ring, use_raw_data, epps, uu, dc_z, x_center, y_center, z_center, dec, bmin, ...
                                        bmax, Vmax, V, gaussK, TOF, TOFSize, sigma_x, TOFCenter, list_mode_format);
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['COSEM-OSL APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.TGV && options.OSL_OSEM
                                grad = TGV(im_vectors.OSEM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_TGV_osem, grad, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSEM-OSL TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.TGV && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    warning('Negative values in BSREM, it is recommended to lower lambda value')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.TGV && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, options.lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM-MAP TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute RBI-OSL with TGV
                            elseif options.TGV && options.RBI_OSL
                                if verbose
                                    tStart = tic;
                                end
                                grad = TGV(im_vectors.OSEM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                                im_vectors.OSEM_apu = RBI_subiter(im_vectors.OSEM_apu, f_Summ(:,osa_iter), rhs, [], [], D, [], [], options.beta_TGV_rbi, grad);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['RBI-OSL TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute COSEM-OSL with TGV
                            elseif options.TGV && any(options.COSEM_OSL)
                                if verbose
                                    tStart = tic;
                                end
                                grad = TGV(im_vectors.OSEM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                                [im_vectors.OSEM_apu, C_osl] = COSEM_OSL(im_vectors.OSEM_apu, D, options.beta_TGV_cosem, grad, rhs, osa_iter, options.h, ...
                                    C_osl, options.COSEM_OSL, [], [], [], []);
                                if options.COSEM_OSL == 1
                                    im_vectors = ACOSEM_coefficient(im_vectors, options, Nx, Ny, Nz, dx, dz, dy, bx, by, bz, z_det, x, y, yy, xx, NSinos, NSlices, size_x, ...
                                        zmax, norm_input, SinD, pituus, osa_iter, attenuation_correction, normalization_correction, randoms_correction, options.scatter, scatter_input, lor_a_input, ...
                                        xy_index_input, z_index_input, L_input, pseudot, det_per_ring, use_raw_data, epps, uu, dc_z, x_center, y_center, z_center, dec, bmin, ...
                                        bmax, Vmax, V, gaussK, TOF, TOFSize, sigma_x, TOFCenter, list_mode_format);
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['COSEM-OSL TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.NLM && options.OSL_OSEM
                                med = NLM(im_vectors.OSEM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                    options.sigma, epps, Nx, Ny, Nz, options);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_NLM_osem, med, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSEM-OSL NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.NLM && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    warning('Negative values in BSREM, it is recommended to lower lambda value')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.NLM && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, options.lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM-MAP NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute RBI-OSL with NLM
                            elseif options.NLM && options.RBI_OSL
                                if verbose
                                    tStart = tic;
                                end
                                med = NLM(im_vectors.OSEM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                    options.sigma, epps, Nx, Ny, Nz, options);
                                im_vectors.OSEM_apu = RBI_subiter(im_vectors.OSEM_apu, f_Summ(:,osa_iter), rhs, [], [], D, [], [], options.beta_NLM_rbi, med);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['RBI-OSL NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                                % Compute COSEM-OSL with NLM
                            elseif options.NLM && any(options.COSEM_OSL)
                                if verbose
                                    tStart = tic;
                                end
                                med = NLM(im_vectors.OSEM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                    options.sigma, epps, Nx, Ny, Nz, options);
                                [im_vectors.OSEM_apu, C_osl] = COSEM_OSL(im_vectors.OSEM_apu, D, options.beta_NLM_cosem, med, rhs, osa_iter, options.h, ...
                                    C_osl, options.COSEM_OSL, [], [], [], []);
                                if options.COSEM_OSL == 1
                                    im_vectors = ACOSEM_coefficient(im_vectors, options, Nx, Ny, Nz, dx, dz, dy, bx, by, bz, z_det, x, y, yy, xx, NSinos, NSlices, size_x, ...
                                        zmax, norm_input, SinD, pituus, osa_iter, attenuation_correction, normalization_correction, randoms_correction, options.scatter, scatter_input, lor_a_input, ...
                                        xy_index_input, z_index_input, L_input, pseudot, det_per_ring, use_raw_data, epps, uu, dc_z, x_center, y_center, z_center, dec, bmin, ...
                                        bmax, Vmax, V, gaussK, TOF, TOFSize, sigma_x, TOFCenter, list_mode_format);
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['COSEM-OSL NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            end
                            
                            clear Summ rhs
                            im_vectors.OSEM_apu(im_vectors.OSEM_apu < epps) = epps;
                        end
                        if options.save_iter
                            iter_n = iter + 1;
                        else
                            iter_n = 1;
                        end
                        if options.osem
                            im_vectors.OSEM(:,iter_n) = im_vectors.OSEM_apu;
                        elseif options.ramla
                            im_vectors.RAMLA(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.rosem
                            im_vectors.ROSEM(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.drama
                            im_vectors.DRAMA(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.rbi
                            im_vectors.RBI(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.cosem
                            im_vectors.COSEM(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.ecosem
                            im_vectors.ECOSEM(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.acosem
                            im_vectors.ACOSEM(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.MRP && options.OSL_OSEM
                            im_vectors.MRP_OSL(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.MRP && options.BSREM
                            med = MRP(im_vectors.OSEM_apu, options.medx, options.medy, options.medz, Nx, Ny, Nz, epps, options.tr_offsets, options.med_no_norm);
                            im_vectors.MRP_BSREM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam, iter, SinD, randoms_correction, is_transposed, ...
                                options.beta_mrp_bsrem, med, epps);
                            im_vectors.OSEM_apu = im_vectors.MRP_BSREM(:, iter_n);
                        elseif options.MRP && options.ROSEM_MAP
                            med = MRP(im_vectors.OSEM_apu, options.medx, options.medy, options.medz, Nx, Ny, Nz, epps, options.tr_offsets, options.med_no_norm);
                            im_vectors.MRP_ROSEM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, ...
                                options.beta_mrp_rosem, med, epps);
                            im_vectors.OSEM_apu = im_vectors.MRP_ROSEM(:, iter_n);
                        elseif options.MRP && options.RBI_OSL
                            im_vectors.MRP_RBI(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.MRP && any(options.COSEM_OSL)
                            im_vectors.MRP_COSEM(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.quad && options.OSL_OSEM
                            im_vectors.Quad_OSL(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.quad && options.BSREM
                            med = Quadratic_prior(im_vectors.OSEM_apu, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                            im_vectors.Quad_BSREM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam, iter, SinD, randoms_correction, is_transposed, options.beta_quad_bsrem, ...
                                med, epps);
                            im_vectors.OSEM_apu = im_vectors.Quad_BSREM(:, iter_n);
                        elseif options.quad && options.ROSEM_MAP
                            med = Quadratic_prior(im_vectors.OSEM_apu, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                            im_vectors.Quad_ROSEM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, ...
                                options.beta_quad_rosem, med, epps);
                            im_vectors.OSEM_apu = im_vectors.Quad_ROSEM(:, iter_n);
                        elseif options.quad && options.RBI_OSL
                            im_vectors.Quad_RBI(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.quad && any(options.COSEM_OSL)
                            im_vectors.Quad_COSEM(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.Huber && options.OSL_OSEM
                            im_vectors.Huber_OSL(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.Huber && options.BSREM
                            med = Huberratic_prior(im_vectors.OSEM_apu, options.weights_huber, Nx, Ny, Nz, Ndx, Ndy, Ndz, options.huber_delta);
                            im_vectors.Huber_BSREM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam, iter, SinD, randoms_correction, is_transposed, ...
                                options.beta_huber_bsrem, med, epps);
                            im_vectors.OSEM_apu = im_vectors.Huber_BSREM(:, iter_n);
                        elseif options.Huber && options.ROSEM_MAP
                            med = Huberratic_prior(im_vectors.OSEM_apu, options.weights_huber, Nx, Ny, Nz, Ndx, Ndy, Ndz, options.huber_delta);
                            im_vectors.Huber_ROSEM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, ...
                                options.beta_huber_rosem, med, epps);
                            im_vectors.OSEM_apu = im_vectors.Huber_ROSEM(:, iter_n);
                        elseif options.Huber && options.RBI_OSL
                            im_vectors.Huber_RBI(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.Huber && any(options.COSEM_OSL)
                            im_vectors.Huber_COSEM(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.L && options.OSL_OSEM
                            im_vectors.L_OSL(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.L && options.BSREM
                            med = L_filter(im_vectors.OSEM_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            im_vectors.L_BSREM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam, iter, SinD, randoms_correction, is_transposed, options.beta_L_bsrem, ...
                                med, epps);
                            im_vectors.OSEM_apu = im_vectors.L_BSREM(:, iter_n);
                        elseif options.L && options.ROSEM_MAP
                            med = L_filter(im_vectors.OSEM_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            im_vectors.L_ROSEM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, ...
                                options.beta_L_rosem, med, epps);
                            im_vectors.OSEM_apu = im_vectors.L_ROSEM(:, iter_n);
                        elseif options.L && options.RBI_OSL
                            im_vectors.L_RBI(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.L && any(options.COSEM_OSL)
                            im_vectors.L_COSEM(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.FMH && options.OSL_OSEM
                            im_vectors.FMH_OSL(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.FMH && options.BSREM
                            med = FMH(im_vectors.OSEM_apu, options.tr_offsets, options.fmh_weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            im_vectors.FMH_BSREM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam, iter, SinD, randoms_correction, is_transposed, options.beta_fmh_bsrem, ...
                                med, epps);
                            im_vectors.OSEM_apu = im_vectors.FMH_BSREM(:, iter_n);
                        elseif options.FMH && options.ROSEM_MAP
                            med = FMH(im_vectors.OSEM_apu, options.tr_offsets, options.fmh_weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            im_vectors.FMH_ROSEM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, ...
                                options.beta_fmh_rosem, med, epps);
                            im_vectors.OSEM_apu = im_vectors.FMH_ROSEM(:, iter_n);
                        elseif options.FMH && options.RBI_OSL
                            im_vectors.FMH_RBI(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.FMH && any(options.COSEM_OSL)
                            im_vectors.FMH_COSEM(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.weighted_mean && options.OSL_OSEM
                            im_vectors.Weighted_OSL(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.weighted_mean && options.BSREM
                            med = Weighted_mean(im_vectors.OSEM_apu, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.med_no_norm);
                            im_vectors.Weighted_BSREM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam, iter, SinD, randoms_correction, is_transposed, ...
                                options.beta_weighted_bsrem, med, epps);
                            im_vectors.OSEM_apu = im_vectors.Weighted_BSREM(:, iter_n);
                        elseif options.weighted_mean && options.ROSEM_MAP
                            med = Weighted_mean(im_vectors.OSEM_apu, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.med_no_norm);
                            im_vectors.Weighted_ROSEM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, ...
                                options.beta_weighted_rosem, ...
                                med, epps);
                            im_vectors.OSEM_apu = im_vectors.Weighted_ROSEM(:, iter_n);
                        elseif options.weighted_mean && options.RBI_OSL
                            im_vectors.Weighted_RBI(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.weighted_mean && any(options.COSEM_OSL)
                            im_vectors.Weighted_COSEM(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.TV && options.OSL_OSEM
                            im_vectors.TV_OSL(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.TV && options.BSREM
                            grad = TVpriorFinal(im_vectors.OSEM_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, options.tr_offsets);
                            im_vectors.TV_BSREM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam, iter, SinD, randoms_correction, is_transposed, options.beta_TV_bsrem, ...
                                grad, epps);
                            im_vectors.OSEM_apu = im_vectors.TV_BSREM(:, iter_n);
                        elseif options.TV && options.ROSEM_MAP
                            grad = TVpriorFinal(im_vectors.OSEM_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, options.tr_offsets);
                            im_vectors.TV_ROSEM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, ...
                                options.beta_TV_rosem, grad, epps);
                            im_vectors.OSEM_apu = im_vectors.TV_ROSEM(:, iter_n);
                        elseif options.TV && options.RBI_OSL
                            im_vectors.TV_RBI(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.TV && any(options.COSEM_OSL)
                            im_vectors.TV_COSEM(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.AD && options.OSL_OSEM
                            im_vectors.AD_OSL(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.AD && options.BSREM
                            med = AD(im_vectors.OSEM_apu, options.FluxType, Nx, Ny, Nz, options);
                            im_vectors.AD_BSREM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam, iter, SinD, randoms_correction, is_transposed, options.beta_ad_bsrem, med, epps);
                            im_vectors.OSEM_apu = im_vectors.AD_BSREM(:, iter_n);
                        elseif options.AD && options.ROSEM_MAP
                            med = AD(im_vectors.OSEM_apu, options.FluxType, Nx, Ny, Nz, options);
                            im_vectors.AD_ROSEM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, ...
                                options.beta_ad_rosem, med, epps);
                            im_vectors.OSEM_apu = im_vectors.AD_ROSEM(:, iter_n);
                        elseif options.AD && options.RBI_OSL
                            im_vectors.AD_RBI(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.AD && any(options.COSEM_OSL)
                            im_vectors.AD_COSEM(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.APLS && options.OSL_OSEM
                            im_vectors.APLS_OSL(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.APLS && options.BSREM
                            grad = TVpriorFinal(im_vectors.OSEM_apu, 0, Nx, Ny, Nz, true, options, 5);
                            im_vectors.APLS_BSREM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam, iter, SinD, randoms_correction, is_transposed, options.beta_APLS_bsrem, ...
                                grad, epps);
                            im_vectors.OSEM_apu = im_vectors.APLS_BSREM(:, iter_n);
                        elseif options.APLS && options.ROSEM_MAP
                            grad = TVpriorFinal(im_vectors.OSEM_apu, 0, Nx, Ny, Nz, true, options, 5);
                            im_vectors.APLS_ROSEM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, ...
                                options.beta_APLS_rosem, grad, epps);
                            im_vectors.OSEM_apu = im_vectors.APLS_ROSEM(:, iter_n);
                        elseif options.APLS && options.RBI_OSL
                            im_vectors.APLS_RBI(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.APLS && any(options.COSEM_OSL)
                            im_vectors.APLS_COSEM(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.TGV && options.OSL_OSEM
                            im_vectors.TGV_OSL(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.TGV && options.BSREM
                            grad = TGV(im_vectors.OSEM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                            im_vectors.TGV_BSREM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam, iter, SinD, randoms_correction, is_transposed, options.beta_TGV_bsrem, ...
                                grad, epps);
                            im_vectors.OSEM_apu = im_vectors.TGV_BSREM(:, iter_n);
                        elseif options.TGV && options.ROSEM_MAP
                            grad = TGV(im_vectors.OSEM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                            im_vectors.TGV_ROSEM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, ...
                                options.beta_TGV_rosem, grad, epps);
                            im_vectors.OSEM_apu = im_vectors.TGV_ROSEM(:, iter_n);
                        elseif options.TGV && options.RBI_OSL
                            im_vectors.TGV_RBI(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.TGV && any(options.COSEM_OSL)
                            im_vectors.TGV_COSEM(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.NLM && options.OSL_OSEM
                            im_vectors.NLM_OSL(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.NLM && options.BSREM
                            med = NLM(im_vectors.OSEM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                options.sigma, epps, Nx, Ny, Nz, options);
                            im_vectors.NLM_BSREM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, ...
                                options.beta_NLM_bsrem, med, epps);
                            im_vectors.OSEM_apu = im_vectors.NLM_BSREM(:, iter_n);
                        elseif options.NLM && options.ROSEM_MAP
                            med = NLM(im_vectors.OSEM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                options.sigma, epps, Nx, Ny, Nz, options);
                            im_vectors.NLM_ROSEM(:, iter_n) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, ...
                                options.beta_NLM_rosem, med, epps);
                            im_vectors.OSEM_apu = im_vectors.NLM_ROSEM(:, iter_n);
                        elseif options.NLM && options.RBI_OSL
                            im_vectors.NLM_RBI(:, iter_n) = im_vectors.OSEM_apu;
                        elseif options.NLM && any(options.COSEM_OSL)
                            im_vectors.NLM_COSEM(:, iter_n) = im_vectors.OSEM_apu;
                        end
                        if options.use_psf && options.deblurring && (options.save_iter || (~options.save_iter && iter == options.Niter))
                            im_vectors = computeDeblur(im_vectors, options, iter, gaussK, Nx, Ny, Nz);
                        end
                        if verbose
                            tElapsed = toc(tStart_iter);
                            disp(['OS iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        end
                        no_norm = true;
                    end
                    % Are MLEM-methods used
                    if MLEM_bool
                        if verbose
                            tStart = tic;
                        end
%                         TOFOffset = int64(0);
                        if iter == 1 && llo == 1
                            if no_norm && OS_bool
                                f_Summ_ml = sum(f_Summ,2);
                            else
                                f_Summ_ml = zeros(Nx*Ny*Nz,1);
                            end
                        end
                        if randoms_correction
                            if iscell(options.SinDelayed)
                                SinD = double(options.SinDelayed{llo});
                            else
                                SinD = double(options.SinDelayed);
                            end
                            if issparse(SinD)
                                SinD = (full(SinD));
                            end
                            SinD = SinD(:);
                            if TOF
                                SinD = SinD / options.TOF_bins;
                            end
                        else
                            SinD = 0;
                        end
                        if use_raw_data
                            xy_index = uint32(0);
                            z_index = uint32(0);
                            TOFSize = int64(size(LL,1));
                        else
                            LL = uint16(0);
                            TOFSize = int64(numel(xy_index));
                        end
                        if ~options.precompute_lor
                            lor_a = uint16(0);
                        end
                        if options.use_psf
                            MLEM_apu = computeConvolution(im_vectors.MLEM_apu, options, Nx, Ny, Nz, gaussK);
                        else
                            MLEM_apu = im_vectors.MLEM_apu;
                        end
                        
                        if options.projector_type == 1
                            if exist('OCTAVE_VERSION','builtin') == 0
                                [Summ, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                    options.normalization, SinD, TOFSize, attenuation_correction, normalization_correction, randoms_correction,...
                                    options.scatter, options.ScatterC, options.global_correction_factor, lor_a, xy_index, z_index, NSinos, LL, pseudot, det_per_ring, ...
                                    TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                    (use_raw_data), uint32(1), epps, double(Sino), MLEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, false, ...
                                    list_mode_format, options.n_rays_transaxial, options.n_rays_axial, dc_z);
                            elseif exist('OCTAVE_VERSION','builtin') == 5
                                [Summ, rhs] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                    options.normalization, SinD, pituus(end) - pituus(1), attenuation_correction, normalization_correction, randoms_correction,...
                                    options.scatter, options.ScatterC, options.global_correction_factor, lor_a, xy_index, z_index, NSinos, LL, pseudot, det_per_ring, ...
                                    TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, ...
                                    (use_raw_data), uint32(1), epps, double(Sino), MLEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, false, ...
                                    list_mode_format, options.n_rays_transaxial, options.n_rays_axial, dc_z);
                            end
                        elseif options.projector_type == 2
                            if exist('OCTAVE_VERSION','builtin') == 0
                                [Summ, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                    options.normalization, SinD, pituus(end) - pituus(1), attenuation_correction, normalization_correction, randoms_correction, ...
                                    options.scatter, options.ScatterC, options.global_correction_factor, lor_a, xy_index, z_index, NSinos, LL, pseudot, det_per_ring, ...
                                    TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, (use_raw_data), uint32(1), ...
                                    epps, double(Sino), MLEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, false, list_mode_format, ...
                                    options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z);
                            elseif exist('OCTAVE_VERSION','builtin') == 5
                                [Summ, rhs] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                    options.normalization, SinD, pituus(end) - pituus(1), attenuation_correction, normalization_correction, randoms_correction, ...
                                    options.scatter, options.ScatterC, options.global_correction_factor, lor_a, xy_index, z_index, NSinos, LL, pseudot, det_per_ring, ...
                                    TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, (use_raw_data), uint32(1), ...
                                    epps, double(Sino), MLEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, false, list_mode_format, ...
                                    options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z);
                            end
                        elseif options.projector_type == 3
                            if exist('OCTAVE_VERSION','builtin') == 0
                                [Summ, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                    options.normalization, SinD, pituus(end) - pituus(1), attenuation_correction, normalization_correction, randoms_correction, ...
                                    options.scatter, options.ScatterC, options.global_correction_factor, lor_a, xy_index, z_index, NSinos, LL, pseudot, det_per_ring, ...
                                    TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, (use_raw_data), uint32(1), ...
                                    epps, double(Sino), MLEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, false, list_mode_format, x_center, ...
                                    y_center, z_center, bmin, bmax, Vmax, V);
                            elseif exist('OCTAVE_VERSION','builtin') == 5
                                [Summ, rhs] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                    options.normalization, SinD, pituus(end) - pituus(1), attenuation_correction, normalization_correction, randoms_correction, ...
                                    options.scatter, options.ScatterC, options.global_correction_factor, lor_a, xy_index, z_index, NSinos, LL, pseudot, det_per_ring, ...
                                    TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, (use_raw_data), uint32(1), ...
                                    epps, double(Sino), MLEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, false, list_mode_format, x_center, ...
                                    y_center, z_center, bmin, bmax, Vmax, V);
                            end
                        else
                            error('Unsupported projector')
                        end
                        
                        if iter == 1 && ~OS_bool && llo == 1
                            if options.use_psf
                                Summ = computeConvolution(Summ, options, Nx, Ny, Nz, gaussK);
                            end
                            Summ(Summ < epps) = epps;
                            f_Summ_ml = Summ;
                        end
                        if options.use_psf
                            rhs = computeConvolution(rhs, options, Nx, Ny, Nz, gaussK);
                        end
                        if options.save_iter
                            iter_n = iter + 1;
                        else
                            iter_n = 1;
                        end
                        if options.mlem
                            im_vectors.MLEM_apu = MLEM_im(im_vectors.MLEM_apu, f_Summ_ml, epps, rhs);
                            im_vectors.MLEM(:, iter_n) = im_vectors.MLEM_apu;
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['MLEM iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.MRP && options.OSL_MLEM
                            med = MRP(im_vectors.MLEM_apu, options.medx, options.medy, options.medz, Nx, Ny, Nz, epps, options.tr_offsets, options.med_no_norm);
                            im_vectors.MLEM_apu = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_mrp_osem, med, epps, rhs);
                            im_vectors.MRP_MLEM(:, iter_n) = im_vectors.MLEM_apu;
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL MRP iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.quad && options.OSL_MLEM
                            med = Quadratic_prior(im_vectors.MLEM_apu, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                            im_vectors.MLEM_apu = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_quad_osem, med, epps, rhs);
                            im_vectors.Quad_MLEM(:, iter_n) = im_vectors.MLEM_apu;
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL Quadratic iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.Huber && options.OSL_MLEM
                            med = Huber_prior(im_vectors.MLEM_apu, options.weights_huber, Nx, Ny, Nz, Ndx, Ndy, Ndz, options.huber_delta);
                            im_vectors.MLEM_apu = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_huber_osem, med, epps, rhs);
                            im_vectors.Huber_MLEM(:, iter_n) = im_vectors.MLEM_apu;
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL Huber iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.L && options.OSL_MLEM
                            med = L_filter(im_vectors.MLEM_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            im_vectors.MLEM_apu = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_L_osem, med, epps, rhs);
                            im_vectors.L_MLEM(:, iter_n) = im_vectors.MLEM_apu;
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL L-filter iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.FMH && options.OSL_MLEM
                            med = FMH(im_vectors.MLEM_apu, options.tr_offsets, options.fmh_weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            im_vectors.MLEM_apu = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_fmh_osem, med, epps, rhs);
                            im_vectors.FMH_MLEM(:, iter_n) = im_vectors.MLEM_apu;
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL FMH iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.weighted_mean && options.OSL_MLEM
                            med = Weighted_mean(im_vectors.MLEM_apu, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.med_no_norm);
                            im_vectors.MLEM_apu = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_weighted_osem, med, rhs);
                            im_vectors.Weighted_MLEM(:, iter_n) = im_vectors.MLEM_apu;
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL weighted mean iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.TV && options.OSL_MLEM
                            grad = TVpriorFinal(im_vectors.MLEM_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, ...
                                options.tr_offsets);
                            im_vectors.MLEM_apu = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_TV_osem, grad, rhs);
                            im_vectors.TV_MLEM(:, iter_n) = im_vectors.MLEM_apu;
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL TV sub-iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.AD && options.OSL_MLEM
                            if iter > 1
                                med = AD(im_vectors.MLEM_apu, options.FluxType, Nx, Ny, Nz, options);
                                im_vectors.MLEM_apu = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_ad_osem, med, epps, rhs);
                            else
                                im_vectors.MLEM_apu = OSEM_im(im_vectors.MLEM_apu, rhs, f_Summ_ml(:,osa_iter), epps);
                            end
                            im_vectors.AD_MLEM(:, iter_n) = im_vectors.MLEM_apu;
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL AD iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.APLS && options.OSL_MLEM
                            grad = TVpriorFinal(im_vectors.MLEM_apu, [], Nx, Ny, Nz, true, options, 5);
                            im_vectors.MLEM_apu = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_APLS_osem, grad, rhs);
                            im_vectors.APLS_MLEM(:, iter_n) = im_vectors.MLEM_apu;
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL APLS iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.TGV && options.OSL_MLEM
                            grad = TGV(im_vectors.MLEM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                            im_vectors.MLEM_apu = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_TGV_osem, grad, rhs);
                            im_vectors.TGV_MLEM(:, iter_n) = im_vectors.MLEM_apu;
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL TGV iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.NLM && options.OSL_MLEM
                            med = NLM(im_vectors.MLEM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                options.sigma, epps, Nx, Ny, Nz, options);
                            im_vectors.MLEM_apu = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_NLM_osem, med, epps, rhs);
                            im_vectors.NLM_MLEM(:, iter_n) = im_vectors.MLEM_apu;
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM-OSL NLM iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.use_psf && options.deblurring && (options.save_iter || (~options.save_iter && iter == options.Niter))
                            im_vectors = computeDeblurMLEM(im_vectors, options, iter, gaussK, Nx, Ny, Nz);
                        end
                        no_norm = true;
                    end
                end
            end
            pause(0.1)
            if options.implementation ~= 2 && options.implementation ~= 3
                im_vectors = reshape_vectors(im_vectors, options);
            end
            pz = images_to_cell(im_vectors, llo, pz, options);
            if partitions > 1 && options.verbose
                disp(['Reconstructions for timestep ' num2str(llo) ' completed'])
            end
            if partitions > 1
                im_vectors = form_image_vectors(options, N);
            end
        end
        %% Implementation 2
        % OpenCL matrix free
        % Uses ArrayFire libraries
        % Only C++ code (no pure MATLAB implementation)
        % Supports all features as implementation 1 except for NLM
    elseif options.implementation == 2
        %%
        options = double_to_single(options);
        
        if partitions == 1
            if ~iscell(options.SinM)
                options.SinM = {options.SinM};
            end
        end
        if options.scatter_correction && ~options.subtract_scatter
            if ~iscell(options.ScatterC)
                options.ScatterC = {single(options.ScatterC)};
            elseif ~isa(options.ScatterC{1}, 'single')
                for kk = 1 : length(options.ScatterC)
                    options.ScatterC{kk} = single(full(options.ScatterC{kk}));
                end
            end
        else
            options.ScatterC = {single(0)};
        end
        if issparse(options.SinM{1})
            for kk = 1 : length(options.SinM)
                options.SinM{kk} = single(full(options.SinM{kk}));
            end
        elseif ~isa(options.SinM{1}, 'single')
            for kk = 1 : length(options.SinM)
                options.SinM{kk} = single(options.SinM{kk});
            end
        end
        if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
            if ~iscell(options.SinDelayed)
                options.SinDelayed = {options.SinDelayed};
            end
            if issparse(options.SinDelayed{1})
                for kk = 1 : length(options.SinDelayed)
                    options.SinDelayed{kk} = single(full(options.SinDelayed{kk}));
                end
            elseif ~isa(options.SinDelayed{1}, 'single')
                for kk = 1 : length(options.SinDelayed)
                    options.SinDelayed{kk} = single(options.SinDelayed{kk});
                end
            end
        end
        if use_raw_data
            xy_index = uint32(0);
            z_index = uint16(0);
            TOFSize = int64(size(LL,1));
        else
            if isempty(pseudot)
                pseudot = uint32(100000);
            end
            LL = uint16(0);
            TOFSize = int64(size(xy_index,1));
        end
        if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction ...
                && ~options.reconstruct_trues && ~options.reconstruct_scatter
            randoms = uint32(1);
        else
            randoms = uint32(0);
        end
        n_rays = uint16(options.n_rays_transaxial);
        n_rays3D = uint16(options.n_rays_axial);
        if n_rays * n_rays3D > 1
            dc_z = single(z_det(2,1) - z_det(1,1));
        else
            dc_z = single(options.cr_pz);
        end
        
        
        tube_width_xy = single(options.tube_width_xy);
        crystal_size_z = single(options.tube_width_z);
        if (options.projector_type == 1 && (options.precompute_lor || (n_rays + n_rays3D) <= 2)) || options.projector_type == 2 || options.projector_type == 3
            kernel_file = 'multidevice_kernel.cl';
            kernel_path = which(kernel_file);
            kernel_path = strrep(kernel_path, '\', '/');
            kernel_path = strrep(kernel_path, '.cl', '');
            filename = 'OMEGA_matrix_free_OpenCL_binary_device';
            header_directory = strrep(kernel_path,'multidevice_kernel','');
        elseif options.projector_type == 1 && ~options.precompute_lor
            kernel_file = 'multidevice_siddon_no_precomp.cl';
            kernel_path = which(kernel_file);
            kernel_path = strrep(kernel_path, '\', '/');
            kernel_path = strrep(kernel_path, '.cl', '');
            filename = 'OMEGA_matrix_free_OpenCL_binary_device';
            header_directory = strrep(kernel_path,'multidevice_siddon_no_precomp','');
        else
            error('Invalid projector for OpenCL')
        end
        filename = [header_directory, filename];
        if options.use_CUDA
            header_directory = strcat('"', header_directory);
            header_directory = strcat(header_directory,'"');
        end
        joku = algorithms_char();
        % n_rekos = uint32(sum(rekot(~contains(joku,'MLEM'))));
        n_rekos = uint32(sum(rekot(cellfun('isempty',strfind(joku,'MLEM')))));
        n_rekos_mlem = uint32(sum(rekot(~cellfun('isempty',strfind(joku,'MLEM')))));
        reko_type = zeros(length(rekot),1,'uint8');
        reko_type(~cellfun('isempty',strfind(joku,'MBSREM'))) = 1;
        reko_type(~cellfun('isempty',strfind(joku,'MRAMLA'))) = 1;
        reko_type(~cellfun('isempty',strfind(joku,'COSEM'))) = 2;
        reko_type(~cellfun('isempty',strfind(joku,'ACOSEM'))) = 3;
        reko_type(~cellfun('isempty',strfind(joku,'OSL-COSEM')) & options.COSEM_OSL == 1) = 2;
        reko_type(~cellfun('isempty',strfind(joku,'OSL-COSEM')) & options.COSEM_OSL == 2) = 3;
        ind = cellfun('isempty',strfind(joku,'MLEM'));
        reko_type = reko_type(rekot & ind);
        joku = joku(rekot & ind);
        if options.ecosem
            if options.cosem && options.osem
                reko_type(~cellfun('isempty',strfind(joku,'ECOSEM'))) = [];
            elseif options.cosem && ~options.osem
                reko_type(~cellfun('isempty',strfind(joku,'ECOSEM'))) = [];
                reko_type = [0;reko_type];
            elseif ~options.cosem && ~options.osem
                reko_type = [0;reko_type];
            end
        end
        reko_type_mlem = zeros(n_rekos_mlem,1,'uint8');
        % reko_type(contains(joku,'MBSREM')) = 1;
        % reko_type(contains(joku,'MRAMLA')) = 1;
        % reko_type(contains(joku,'COSEM')) = 2;
        % reko_type(contains(joku,'ACOSEM')) = 3;
        % reko_type(contains(joku,'OSL-COSEM') & options.COSEM_OSL == 1) = 2;
        % reko_type(contains(joku,'OSL-COSEM') & options.COSEM_OSL == 2) = 3;
        tic
        if ~options.use_CUDA
            [pz] = OpenCL_matrixfree( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end), NSinos, single(NSlices), size_x, zmax, NSinos, ...
                options.verbose, LL, pseudot, det_per_ring, TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), int32(dec), device, uint8(use_raw_data), ...
                filename, uint32(0), options.use_psf, header_directory, options.vaimennus, options.normalization, pituus, uint32(attenuation_correction), ...
                uint32(normalization_correction), uint32(Niter), uint32(subsets), uint8(rekot), single(epps), lor_a, xy_index, z_index, any(n_rekos), tube_width_xy, ...
                crystal_size_z, x_center, y_center, z_center, options.SinDelayed, randoms, uint32(options.projector_type), options.precompute_lor, n_rays, n_rays3D, ...
                dc_z, options, options.SinM, uint32(partitions), logical(options.use_64bit_atomics), n_rekos, n_rekos_mlem, reko_type, reko_type_mlem, ...
                options.global_correction_factor, bmin, bmax, Vmax, V, gaussK);
        else
            header_directory = strrep(header_directory,'"','');
            [pz] = CUDA_matrixfree( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end), NSinos, single(NSlices), size_x, zmax, NSinos, ...
                options.verbose, LL, pseudot, det_per_ring, TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), int32(dec), device, uint8(use_raw_data), ...
                filename, uint32(0), options.use_psf, header_directory, options.vaimennus, options.normalization, pituus, uint32(attenuation_correction), ...
                uint32(normalization_correction), uint32(Niter), uint32(subsets), uint8(rekot), single(epps), lor_a, xy_index, z_index, any(n_rekos), tube_width_xy, ...
                crystal_size_z, x_center, y_center, z_center, options.SinDelayed, randoms, uint32(options.projector_type), options.precompute_lor, n_rays, n_rays3D, ...
                dc_z, options, options.SinM, uint32(partitions), logical(options.use_64bit_atomics), n_rekos, n_rekos_mlem, reko_type, reko_type_mlem, ...
                options.global_correction_factor, bmin, bmax, Vmax, V, gaussK);
        end
        toc
        
        pz(end,:) = [];
        %% Implementation 3
        % OpenCL matrix free with multiple device support
        % Uses pure OpenCL
        % Only C++ code (no pure MATLAB implementation)
        % Supports only MLEM and OSEM
    elseif options.implementation == 3
        options = double_to_single(options);
        
        if partitions == 1
            if ~iscell(options.SinM)
                options.SinM = {options.SinM};
            end
        end
        if options.scatter_correction && ~options.subtract_scatter
            if ~iscell(options.ScatterC)
                options.ScatterC = {single(options.ScatterC)};
            elseif ~isa(options.ScatterC{1}, 'single')
                for kk = 1 : length(options.ScatterC)
                    options.ScatterC{kk} = single(full(options.ScatterC{kk}));
                end
            end
        else
            options.ScatterC = {single(0)};
        end
        if issparse(options.SinM{1})
            for kk = 1 : length(options.SinM)
                options.SinM{kk} = single(full(options.SinM{kk}));
            end
        elseif ~isa(options.SinM{1}, 'single')
            for kk = 1 : length(options.SinM)
                options.SinM{kk} = single(options.SinM{kk});
            end
        end
        if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
            if ~iscell(options.SinDelayed)
                options.SinDelayed = {options.SinDelayed};
            end
            if issparse(options.SinDelayed{1})
                for kk = 1 : length(options.SinDelayed)
                    options.SinDelayed{kk} = single(full(options.SinDelayed{kk}));
                end
            elseif ~isa(options.SinDelayed{1}, 'single')
                for kk = 1 : length(options.SinDelayed)
                    options.SinDelayed{kk} = single(options.SinDelayed{kk});
                end
            end
        end
        if use_raw_data
            xy_index = uint32(0);
            z_index = uint16(0);
            TOFSize = int64(size(LL,1));
        else
            if isempty(pseudot)
                pseudot = uint32(100000);
            end
            LL = uint16(0);
            TOFSize = int64(size(xy_index,1));
        end
        if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction ...
                && ~options.reconstruct_trues && ~options.reconstruct_scatter
            randoms = uint32(1);
        else
            randoms = uint32(0);
        end
        n_rays = uint16(options.n_rays_transaxial);
        n_rays3D = uint16(options.n_rays_axial);
        dc_z = single(z_det(2,1) - z_det(1,1));
        
        tube_width_xy = single(options.tube_width_xy);
        crystal_size_z = single(options.tube_width_z);
        if (options.projector_type == 1 && (options.precompute_lor || (n_rays + n_rays3D) <= 2)) || options.projector_type == 2 || options.projector_type == 3
            kernel_file = 'multidevice_kernel.cl';
            kernel_path = which(kernel_file);
            kernel_path = strrep(kernel_path, '\', '/');
            kernel_path = strrep(kernel_path, '.cl', '');
            filename = 'OMEGA_matrix_free_OpenCL_binary_device';
            header_directory = strrep(kernel_path,'multidevice_kernel','');
        elseif options.projector_type == 1 && ~options.precompute_lor
            kernel_file = 'multidevice_siddon_no_precomp.cl';
            kernel_path = which(kernel_file);
            kernel_path = strrep(kernel_path, '\', '/');
            kernel_path = strrep(kernel_path, '.cl', '');
            filename = 'OMEGA_matrix_free_OpenCL_binary_device';
            header_directory = strrep(kernel_path,'multidevice_siddon_no_precomp','');
        else
            error('Invalid projector for OpenCL')
        end
        filename = [header_directory, filename];
%         header_directory = strcat('-I "', header_directory);
%         header_directory = strcat(header_directory,'"');
        tic
        [tz] = OpenCL_matrixfree_multi_gpu( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end), single(NSlices), size_x, zmax, options.verbose, ...
            LL, pseudot, det_per_ring, TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), int32(dec), uint32(options.use_device), filename, uint8(use_raw_data), single(options.cpu_to_gpu_factor), uint32(1), header_directory, ...
            options.vaimennus, options.normalization, pituus, uint32(attenuation_correction), uint32(normalization_correction), lor_a, xy_index, z_index, tube_width_xy, ...
            crystal_size_z, x_center, y_center, z_center, options.SinDelayed, randoms, uint32(options.projector_type), options.precompute_lor, ...
            n_rays, n_rays3D, dc_z, options.SinM, logical(options.use_64bit_atomics), NSinos, NSinos, uint32(Niter), uint32(subsets), uint8(rekot), ...
            single(epps), uint32(partitions), options.osem, options.use_psf, options.global_correction_factor, bmin, bmax, Vmax, V, gaussK, options);
        toc
        
        
        for ll = 1 : options.partitions
            apu = tz{1,ll};
            if options.save_iter
                apu(:,1) = options.x0;
                if options.mlem
                    pz{1,ll} = reshape(apu, Nx, Ny, Nz, Niter + 1);
                elseif options.osem
                    pz{2,ll} = reshape(apu, Nx, Ny, Nz, Niter + 1);
                end
            else
                if options.mlem
                    pz{1,ll} = reshape(apu, Nx, Ny, Nz);
                elseif options.osem
                    pz{2,ll} = reshape(apu, Nx, Ny, Nz);
                end
            end
        end
        clear tz
    else
        error('Unsupported reconstruction method.');
    end
end

% Save various image properties, e.g. matrix size, sinogram dimensions, FOV
% size, regularization parameters, etc.
pz = save_image_properties(options, pz, subsets);

end
