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
% Copyright (C) 2019  Ville-Veikko Wettenhovi, Samuli Summala
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

if ~isfield(options,'attenuation_phase')
    options.attenuation_phase = false;
end

folder = fileparts(which('reconstructions_main.m'));
folder = strrep(folder, 'source','mat-files/');
folder = strrep(folder, '\','/');

disp('Preparing for reconstruction')

% Load the measurement data if it does not exist in options.options.SinM
if options.use_raw_data
    if options.reconstruct_trues == false && ~options.reconstruct_scatter && (isfield(options, 'coincidences') == 0 || options.precompute_all) && options.use_machine < 2
        if options.partitions == 1
            if options.use_ASCII && options.use_machine == 0
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_static_raw_ASCII.mat'], 'coincidences');
            elseif options.use_LMF && options.use_machine == 0
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_static_raw_LMF.mat'], 'coincidences');
            elseif options.use_root && options.use_machine == 0
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_static_raw_root.mat'], 'coincidences');
            else
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_static_raw_listmode.mat'], 'coincidences');
            end
        else
            if options.use_ASCII && options.use_machine == 0
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_ASCII.mat'], 'coincidences');
            elseif options.use_LMF && options.use_machine == 0
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_LMF.mat'], 'coincidences');
            elseif options.use_root && options.use_machine == 0
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_root.mat'], 'coincidences');
            else
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_listmode.mat'], 'coincidences');
            end
        end
        %         options.SinM = coincidences;
    elseif ~options.reconstruct_trues && ~options.reconstruct_scatter && isfield(options, 'coincidences')
        options.SinM = options.coincidences;
    elseif options.reconstruct_trues
        % Load Trues
        if options.partitions == 1
            if options.use_ASCII
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_static_raw_ASCII.mat'], 'true_coincidences');
            elseif options.use_LMF
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_static_raw_LMF.mat'], 'true_coincidences');
            elseif options.use_root
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_static_raw_root.mat'], 'true_coincidences');
            else
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_static_raw_real.mat'], 'true_coincidences');
            end
        else
            if options.use_ASCII
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_ASCII.mat'], 'true_coincidences');
            elseif options.use_LMF
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_LMF.mat'], 'true_coincidences');
            elseif options.use_root
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_root.mat'], 'true_coincidences');
            else
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_real.mat'], 'true_coincidences');
            end
        end
    elseif options.reconstruct_scatter
        % Load scattered coincidences
        if options.partitions == 1
            if options.use_ASCII
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_static_raw_ASCII.mat'], 'scattered_coincidences');
            elseif options.use_LMF
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_static_raw_LMF.mat'], 'scattered_coincidences');
            elseif options.use_root
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_static_raw_root.mat'], 'scattered_coincidences');
            else
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_static_raw_real.mat'], 'scattered_coincidences');
            end
        else
            if options.use_ASCII
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_ASCII.mat'], 'scattered_coincidences');
            elseif options.use_LMF
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_LMF.mat'], 'scattered_coincidences');
            elseif options.use_root
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_root.mat'], 'scattered_coincidences');
            else
                options.SinM = loadStructFromFile([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_real.mat'], 'scattered_coincidences');
            end
        end
    end
    % Perform corrections if needed
    if options.randoms_correction && ~options.reconstruct_trues && ~options.reconstruct_scatter
        if (options.use_ASCII || options.use_LMF || options.use_root) && options.use_machine == 0
            if options.partitions == 1
                if options.use_ASCII && options.use_machine == 0
                    options.SinDelayed = loadStructFromFile([options.machine_name '_measurements_' options.name '_static_raw_ASCII.mat'], 'delayed_coincidences');
                elseif options.use_LMF && options.use_machine == 0
                    options.SinDelayed = loadStructFromFile([options.machine_name '_measurements_' options.name '_static_raw_LMF.mat'], 'delayed_coincidences');
                elseif options.use_root && options.use_machine == 0
                    options.SinDelayed = loadStructFromFile([options.machine_name '_measurements_' options.name '_static_raw_root.mat'], 'delayed_coincidences');
                else
                    options.SinDelayed = loadStructFromFile([options.machine_name '_measurements_' options.name '_static_raw_listmode.mat'], 'delayed_coincidences');
                end
            else
                if options.use_ASCII && options.use_machine == 0
                    options.SinDelayed = loadStructFromFile([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                        num2str(options.tot_time) 's_raw_ASCII.mat'], 'delayed_coincidences');
                elseif options.use_LMF && options.use_machine == 0
                    options.SinDelayed = loadStructFromFile([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                        num2str(options.tot_time) 's_raw_LMF.mat'], 'delayed_coincidences');
                elseif options.use_root && options.use_machine == 0
                    options.SinDelayed = loadStructFromFile([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                        num2str(options.tot_time) 's_raw_root.mat'], 'delayed_coincidences');
                else
                    options.SinDelayed = loadStructFromFile([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                        num2str(options.tot_time) 's_raw_listmode.mat'], 'delayed_coincidences');
                end
            end
            if isfield(options,'SinDelayed')
                %                 if ~options.corrections_during_reconstruction
                if iscell(options.SinM) && iscell(options.SinDelayed)
                    for kk = 1 : length(options.SinM)
                        if options.variance_reduction
                            options.SinDelayed{kk} = Randoms_variance_reduction(double(options.SinDelayed{kk}), options);
                        end
                        if options.randoms_smoothing
                            options.SinDelayed{kk} = randoms_smoothing(options.SinDelayed{kk}, options);
                            if ~options.corrections_during_reconstruction
                                options.SinM{kk} = full(options.SinM{kk}) - options.SinDelayed{kk};
                            end
                        elseif ~options.corrections_during_reconstruction
                            options.SinM{kk} = options.SinM{kk} - options.SinDelayed{kk};
                        end
                        if ~options.corrections_during_reconstruction
                            options.SinM{kk}(options.SinM{kk} < 0) = 0;
                        end
                    end
                else
                    if options.variance_reduction
                        options.SinDelayed = Randoms_variance_reduction(double(options.SinDelayed), options);
                    end
                    if options.randoms_smoothing
                        options.SinDelayed = randoms_smoothing(options.SinDelayed, options);
                        if ~options.corrections_during_reconstruction
                            options.SinM = full(options.SinM) - options.SinDelayed;
                        end
                    elseif ~options.corrections_during_reconstruction
                        options.SinM = options.SinM - options.SinDelayed;
                    end
                    if ~options.corrections_during_reconstruction
                        options.SinM(options.SinM < 0) = 0;
                    end
                end
                if ~options.corrections_during_reconstruction
                    options = rmfield(options,'SinDelayed');
                end
                %                 end
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
                                options.SinDelayed{kk} = Randoms_variance_reduction(double(options.SinDelayed{kk}), options);
                            end
                            if options.randoms_smoothing
                                options.SinDelayed{kk} = randoms_smoothing(double(options.SinDelayed{kk}), options);
                                if ~options.corrections_during_reconstruction
                                    options.SinM{kk} = full(options.SinM{kk}) - options.SinDelayed{kk};
                                end
                            elseif ~options.corrections_during_reconstruction
                                options.SinM{kk} = options.SinM{kk} - options.SinDelayed{kk};
                            end
                            if ~options.corrections_during_reconstruction
                                options.SinM{kk}(options.SinM{kk} < 0) = 0;
                            end
                        end
                    else
                        if numel(options.SinDelayed{1}) ~= numel(options.SinM)
                            error('Size mismatch between randoms correction data and measurement data')
                        end
                        options.SinM = options.SinM - double(options.SinDelayed{1}(:));
                        options.SinM(options.SinM < 0) = 0;
                    end
                else
                    if iscell(options.SinM)
                        for kk = 1 : length(options.SinM)
                            if numel(options.SinDelayed) ~= numel(options.SinM{kk})
                                error('Size mismatch between randoms correction data and measurement data')
                            end
                            if options.variance_reduction
                                options.SinDelayed = Randoms_variance_reduction(double(options.SinDelayed), options);
                            end
                            if options.randoms_smoothing
                                options.SinDelayed = randoms_smoothing(double(options.SinDelayed), options);
                                if ~options.corrections_during_reconstruction
                                    options.SinM{kk} = double(full(options.SinM{kk})) - options.SinDelayed;
                                end
                            elseif ~options.corrections_during_reconstruction
                                options.SinM{kk} = (options.SinM{kk}) - double(options.SinDelayed(:));
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
                            options.SinDelayed = Randoms_variance_reduction(double(options.SinDelayed), options);
                        end
                        if options.randoms_smoothing
                            options.SinDelayed = randoms_smoothing(options.SinDelayed, options);
                            if ~options.corrections_during_reconstruction
                                options.SinM = full(options.SinM) - options.SinDelayed;
                            end
                        elseif ~options.corrections_during_reconstruction
                            options.SinM = options.SinM - double(options.SinDelayed(:));
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
                        error('Size mismatch between scatter correction data and measurement data')
                    end
                    %                     if options.variance_reduction
                    %                         options.ScatterC{kk} = Randoms_variance_reduction(double(options.ScatterC{kk}), options);
                    %                     end
                    if options.randoms_smoothing
                        options.ScatterC{kk} = randoms_smoothing(options.ScatterC{kk}, options);
                        options.SinM{kk} = full(options.SinM{kk}) - options.ScatterC{kk};
                    else
                        options.SinM{kk} = options.SinM{kk} - double(options.ScatterC{kk}(:));
                    end
                    options.SinM{kk}(options.SinM{kk} < 0) = 0;
                end
            else
                if numel(options.ScatterC{1}) ~= numel(options.SinM)
                    error('Size mismatch between scatter correction data and measurement data')
                end
                if options.variance_reduction
                    options.ScatterC{1} = Randoms_variance_reduction(double(options.ScatterC{1}), options);
                end
                if options.randoms_smoothing
                    options.ScatterC{1} = randoms_smoothing(options.ScatterC{1}, options);
                    options.SinM = full(options.SinM) - options.ScatterC{1};
                else
                    options.SinM = options.SinM - double(options.ScatterC{1}(:));
                end
                options.SinM(options.SinM < 0) = 0;
            end
        else
            if iscell(options.SinM)
                for kk = 1 : length(options.SinM)
                    if numel(options.ScatterC) ~= numel(options.SinM{kk})
                        error('Size mismatch between scatter correction data and measurement data')
                    end
                    %                     if options.variance_reduction
                    %                         options.ScatterC = Randoms_variance_reduction(double(options.ScatterC), options);
                    %                     end
                    if options.randoms_smoothing
                        options.ScatterC = randoms_smoothing(options.ScatterC, options);
                        options.SinM{kk} = full(options.SinM{kk}) - options.ScatterC;
                    else
                        options.SinM{kk} = options.SinM{kk} - double(options.ScatterC(:));
                    end
                    options.SinM{kk} = options.SinM{kk} - double(options.ScatterC(:));
                    options.SinM{kk}(options.SinM{kk} < 0) = 0;
                end
            else
                if numel(options.ScatterC) ~= numel(options.SinM)
                    error('Size mismatch between scatter correction data and measurement data')
                end
                %                 if options.variance_reduction
                %                     options.ScatterC = Randoms_variance_reduction(double(options.ScatterC), options);
                %                 end
                if options.randoms_smoothing
                    options.ScatterC = randoms_smoothing(options.ScatterC, options);
                    options.SinM = full(options.SinM) - options.ScatterC;
                else
                    options.SinM = options.SinM - double(options.ScatterC(:));
                end
                options.SinM = options.SinM - double(options.ScatterC(:));
                options.SinM(options.SinM < 0) = 0;
            end
        end
    end
    if options.normalization_correction && ~options.corrections_during_reconstruction
        if options.use_user_normalization
            [file, ~] = uigetfile({'*.mat'},'Select normalization datafile');
            if isequal(file, 0)
                error('No file was selected')
            end
            if any(strfind(file, '.nrm'))
                error('Inveon normalization data cannot be used with raw list-mode data')
            else
                data = load(file);
                variables = fields(data);
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
                options.SinM{kk} = options.SinM{kk} .* full(normalization);
            end
        else
            options.SinM = options.SinM .* full(normalization);
        end
    end
    
    clear coincidences options.coincidences true_coincidences delayed_coincidences
    % Sinogram data
else
    if (~options.reconstruct_trues && ~options.reconstruct_scatter) || options.use_machine > 0
        if options.partitions == 1 && isfield(options, 'SinM') == 0
            if (options.randoms_correction || options.scatter_correction || options.normalization_correction) && ~options.corrections_during_reconstruction && ...
                    options.use_machine < 2
                if options.use_machine == 0
                    options.SinM = loadStructFromFile([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                        num2str(options.TotSinos) '_span' num2str(options.span) '.mat'],'SinM');
                elseif  options.use_machine == 1
                    options.SinM = loadStructFromFile([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                        num2str(options.TotSinos) '_span' num2str(options.span) '_listmode.mat'],'SinM');
                end
            elseif options.use_machine < 2
                if options.use_machine == 0
                    options.SinM = loadStructFromFile([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                        num2str(options.TotSinos) '_span' num2str(options.span) '.mat'],'raw_SinM');
                elseif  options.use_machine == 1
                    options.SinM = loadStructFromFile([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                        num2str(options.TotSinos) '_span' num2str(options.span) '_listmode.mat'],'raw_SinM');
                end
            else
                options.SinM = loadStructFromFile([options.machine_name '_' options.name '_sinogram_original_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                    num2str(options.TotSinos) '_span' num2str(options.span) '_machine_sinogram.mat'],'SinM');
            end
        elseif isfield(options, 'SinM') == 0
            if (options.randoms_correction || options.scatter_correction || options.normalization_correction) && ~options.corrections_during_reconstruction ...
                    && options.use_machine < 2
                if options.use_machine == 0
                    options.SinM = loadStructFromFile([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                        num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                        num2str(options.span) '.mat'], 'SinM');
                elseif  options.use_machine == 1
                    options.SinM = loadStructFromFile([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                        num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                        num2str(options.span) '_listmode.mat'], 'SinM');
                end
            elseif options.use_machine < 2
                if options.use_machine == 0
                    options.SinM = loadStructFromFile([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ '
                        num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                        num2str(options.span) '.mat'], 'raw_SinM');
                elseif  options.use_machine == 1
                    options.SinM = loadStructFromFile([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                        num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                        num2str(options.span) '_listmode.mat'], 'raw_SinM');
                end
                %                 options.SinM = raw_SinM;
                %                 clear raw_SinM
            else
                options.SinM = loadStructFromFile([options.machine_name '_' options.name '_sinograms_original_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                    num2str(options.span) '_machine_sinogram.mat'], 'SinM');
            end
        else
        end
    elseif options.reconstruct_trues && options.use_machine == 0
        if options.partitions == 1 && isfield(options, 'SinTrues') == 0
            options.SinM = loadStructFromFile([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                num2str(options.TotSinos) '_span' num2str(options.span) '.mat'],'SinTrues');
        elseif isfield(options, 'SinTrues') == 0
            options.SinM = loadStructFromFile([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                num2str(options.span) '.mat'], 'SinTrues');
        end
    elseif options.reconstruct_scatter && options.use_machine == 0
        if options.partitions == 1 && isfield(options, 'SinScatter') == 0
            options.SinM = loadStructFromFile([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                num2str(options.TotSinos) '_span' num2str(options.span) '.mat'],'SinScatter');
        elseif isfield(options, 'SinScatter') == 0
            options.SinM = loadStructFromFile([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                num2str(options.span) '.mat'], 'SinScatter');
        end
    end
    if options.partitions == 1 && options.randoms_correction && options.corrections_during_reconstruction
        if options.use_machine == 0
            options.SinDelayed = loadStructFromFile([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                num2str(options.TotSinos) '_span' num2str(options.span) '.mat'],'SinDelayed');
        elseif  options.use_machine == 1
            options.SinDelayed = loadStructFromFile([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                num2str(options.TotSinos) '_span' num2str(options.span) '_listmode.mat'],'SinDelayed');
        else
            options = loadDelayedData(options);
        end
    elseif options.randoms_correction && options.corrections_during_reconstruction
        if options.use_machine == 0
            options.SinDelayed = loadStructFromFile([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                num2str(options.span) '.mat'], 'SinDelayed');
        elseif  options.use_machine == 1
            options.SinDelayed = loadStructFromFile([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                num2str(options.span) '_listmode.mat'], 'SinDelayed');
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
    end
end

rekot = reko_maker(options);

% OS = any([options.osem; options.mramla; options.ramla; options.rosem; options.rbi; options.cosem;...
%     options.ecosem; options.acosem; options.OSL_OSEM; options.MBSREM; options.BSREM;...
%     options.ROSEM_MAP; options.RBI_MAP; options.COSEM_MAP]);

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
det_per_ring = options.det_per_ring;
machine_name = options.machine_name;
Nang = options.Nang;
Ndist = options.Ndist;
% cr_pz = options.cr_pz;
use_raw_data = options.use_raw_data;
TotSinos = options.TotSinos;
% attenuation_datafile = options.attenuation_datafile;
partitions = options.partitions;
verbose = options.verbose;
device = uint32(options.use_device);
options.empty_weight = false;
options.MBSREM_prepass = true;
options.tr_offsets = 0;
% options.U = [];
% options.weights = [];
% options.a_L = [];
% options.fmh_weights = [];
% options.weighted_weights = [];
% options.weights_quad = [];
options.MAP = (options.OSL_MLEM || options.OSL_OSEM || options.BSREM || options.MBSREM || options.ROSEM_MAP || options.RBI_MAP...
    || any(options.COSEM_MAP));

MLEM_bool = options.OSL_MLEM || options.mlem;
OS_bool = options.osem || options.rosem || options.ramla || options.OSL_OSEM || options.BSREM || options.ROSEM_MAP;

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

if (options.quad || options.FMH || options.L || options.weighted_mean || options.MRP || options.TV) && options.MAP
    Ndx = options.Ndx;
    Ndy = options.Ndy;
    Ndz = options.Ndz;
end
if options.L && options.MAP
    %     options.a_L = options.a_L;
    if ~isempty(options.a_L)
        if length(options.a_L(:)) < ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
            error(['Weights vector options.a_L is too small, needs to be ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) ' in length'])
        elseif length(options.a_L(:)) > ((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))
            error(['Weights vector options.a_L is too large, needs to be ' num2str(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1))) ' in length'])
        end
    end
end
if (options.quad || options.FMH || options.L || options.weighted_mean) && options.MAP
    %     options.weights = options.weights;
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
if options.FMH && options.MAP
    %     options.fmh_weights = options.fmh_weights;
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
    %     options.weighted_weights = options.weighted_weights;
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
    %     pituus = [pituus(1); pituus(end)];
end

%%

% Compute the indices for the subsets used.
% For Sinogram data, six different methods to select the subsets are
% available. For raw list-mode data, three methods are available.
[index, pituus, subsets] = index_maker(Nx, Ny, Nz, subsets, use_raw_data, machine_name, options, Nang, Ndist, TotSinos, NSinos);



%%

% Diameter of the PET-device (bore) (mm)
R=double(diameter);
% Transaxial FOV (x-direction, horizontal) (mm)
FOVax=double(FOVax);
% Transaxial FOV (y-direction, vertical) (mm)
FOVay=double(FOVay);
% Axial FOV (mm)
axial_fow = double(axial_fov);
% Number of rings
blocks=uint32(rings + length(pseudot) - 1);
% First ring
block1=uint32(0);

NSinos = uint32(NSinos);
NSlices = uint32(Nz);
TotSinos = uint32(TotSinos);

% Coordinates of the detectors
[x, y, z] = get_coordinates(options, blocks);

% Load correction data
[normalization_correction, randoms_correction, options] = set_up_corrections(options, blocks);

if min(z(:)) == 0
    z = z + (axial_fow - max(max(z)))/2;
end

if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
    x = single(x);
    y = single(y);
    z_det = single(z);
else
    x = double(x);
    y = double(y);
    z_det = double(z);
end
clear z


size_x = uint32(size(x,1));

if subsets > 1
    pituus = [uint32(0);cumsum(pituus)];
    if iscell(index)
        index = cell2mat(index);
    end
end

% Compute the necessary indices required for subsets (e.g. the index of
% the detector coordinates for the current LOR)
[options, lor_a, xy_index, z_index, LL, summa, pituus, options.SinM, lor_orth] = form_subset_indices(options, pituus, subsets, index, size_x, y, z_det, rings, false, options.SinM);
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


% Pixels
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
    ind_size = uint32((det_per_ring)^2 / subsets * Nx * (Ny));
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
if options.projector_type == 2
    x_center = xx(1 : end - 1)' + dx/2;
    y_center = yy(1 : end - 1)' + dy/2;
    if options.tube_width_z > 0
        z_center = zz(1 : end - 1)' + dz/2;
    else
        z_center = zz(1);
    end
else
    x_center = xx(1);
    y_center = yy(1);
    z_center = zz(1);
end

%% This computes a whole observation matrix and uses it to compute the MLEM (no on-the-fly calculations)
% NOTE: Only attenuation correction is supported
% This section is largely untested
if precompute_obs_matrix && options.implementation == 1
    
    for llo = 1 : partitions
        
        
        if options.mlem
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

                A = observation_matrix_formation(options);
                D = full(sum(A,1))';
                D(D <= 0) = epps;
            end
            
            MLEM = ones(N,Niter);
            MLEM(:,1) = options.x0(:);
            Sino = double(Sino);
            
            for iter=1:Niter
                if options.mlem
                    MLEM(:,iter+1) = MLEM_im(MLEM(:,iter), D, epps, A, Sino, false);
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
    
    
    [options, D, C_co, C_aco, C_osl, Amin, E] = prepass_phase(options, pituus, index, options.SinM, pseudot, x, y, xx, yy, z_det, dz, dx, dy, bz, bx, by, NSlices, zmax, size_x, block1, blocks,...
        normalization_correction, randoms_correction, xy_index, z_index, lor_a, lor_orth, summa, LL, is_transposed, x_center, y_center, z_center);
    
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
            
            Sino = Sino(:);
            
            if issparse(Sino)
                Sino = (full(Sino));
            end
            
            % Implementation 1
            if options.implementation == 1
                
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
                        % No precomputation done
                        % This is sequential (non-parallel) code
                        % SLOW
                        % Supports PURE MATLAB computation
                        if options.precompute_lor == false
                            if use_raw_data == false
                                % Siddon
                                if options.projector_type == 1 || options.projector_type == 0
                                    if exist('projector_mex','file') == 3
                                        [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                            zmax, options.vaimennus, options.normalization, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, ...
                                            randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, uint16(0), pseudot, det_per_ring, options.verbose, ...
                                            use_raw_data, uint32(2), ind_size, block1, blocks, index(pituus(osa_iter)+1:pituus(osa_iter + 1)), uint32(options.projector_type));
                                    else
                                        % The below lines allow for pure MATLAB
                                        % implemention, i.e. no MEX-files will be
                                        % used. Currently the below function
                                        % uses parfor-loops (requires parallel
                                        % computing toolbox).
                                        % NOTE: The below method is not
                                        % recommended since it is much slower
                                        % method.
                                        [ lor, indices, alkiot, discard] = improved_siddon_atten( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, yy, xx, ...
                                            NSinos, NSlices, options.vaimennus, index(pituus(osa_iter)+1:pituus(osa_iter + 1)), pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction);
                                        alkiot = cell2mat(alkiot);
                                        indices = indices(discard);
                                        indices = cell2mat(indices) - 1;
                                    end
                                    % Orthogonal distance based
                                elseif options.projector_type == 2
                                    [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                        zmax, options.vaimennus, options.normalization, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, ...
                                        randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, uint16(0), pseudot, det_per_ring, options.verbose, ...
                                        use_raw_data, uint32(2), ind_size, block1, blocks, index(pituus(osa_iter)+1:pituus(osa_iter + 1)), uint32(options.projector_type), ...
                                        options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, int32(options.accuracy_factor));
                                else
                                    error('Unsupported projector type')
                                end
                            else
                                L = LL(pituus(osa_iter) * 2 + 1 : pituus(osa_iter + 1) * 2);
                                if options.projector_type == 1
                                    [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                        zmax, options.vaimennus, options.normalization, SinD, uint32(0), attenuation_correction, normalization_correction, ...
                                        randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, L, pseudot, det_per_ring, options.verbose, ...
                                        use_raw_data, uint32(2), ind_size, block1, blocks, uint32(0), uint32(options.projector_type));
                                elseif options.projector_type == 2
                                    [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                        zmax, options.vaimennus, options.normalization, SinD, uint32(0), attenuation_correction, normalization_correction, ...
                                        randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, L, pseudot, det_per_ring, options.verbose, ...
                                        use_raw_data, uint32(2), ind_size, block1, blocks, uint32(0), uint32(options.projector_type), options.tube_width_xy, ...
                                        x_center, y_center, z_center, options.tube_width_z, int32(options.accuracy_factor));
                                else
                                    error('Unsupported projector type')
                                end
                            end
                            if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                                lor = repeat_elem(int32(1:length(lor))',int32(lor));
                            elseif exist('OCTAVE_VERSION','builtin') == 5
                                lor = repelem(int32(1:length(lor)),int32(lor));
                            else
                                lor = repelem(int32(1:length(lor)),int32(lor))';
                            end
                            uu = double(Sino(pituus(osa_iter) + 1 : pituus(osa_iter + 1)));
                            
                            A_length = length(uu);
                            indices = int32(indices) + 1;
                            if verbose
                                tStart = tic;
                            end
                            if options.use_fsparse && exist('fsparse','file') == 3
                                A = fsparse(lor,int32(indices),double(alkiot),[A_length double(N) length(alkiot)]);
                            elseif options.use_fsparse && exist('fsparse','file') == 0
                                warning('options.fsparse set to true, but no FSparse mex-file found. Using regular sparse')
                                A = sparse(double(lor),double(indices),double(alkiot), A_length, double(options.N));
                            else
                                A = sparse(double(lor),double(indices),double(alkiot), A_length, double(N));
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
                            else
                                L_input = uint16(0);
                                xy_index_input = xy_index(pituus(osa_iter)+1:pituus(osa_iter + 1));
                                z_index_input = z_index(pituus(osa_iter)+1:pituus(osa_iter + 1));
                            end
                            if options.projector_type == 2
                                lor2 = [uint64(0); cumsum(uint64(lor_orth(pituus(osa_iter)+1:pituus(osa_iter + 1))))];
                            else
                                lor2 = [uint64(0); cumsum(uint64(lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1))))];
                            end
                            [A, ll] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                options.normalization, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction,...
                                randoms_correction, lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1)), xy_index_input, z_index_input, NSinos, L_input, pseudot, ...
                                det_per_ring, options.verbose, use_raw_data, uint32(0), lor2, summa(osa_iter), options.attenuation_phase, uint32(options.projector_type), ...
                                options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, int32(options.accuracy_factor));
                            uu = double(Sino(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                            if options.attenuation_phase
                                uu = uu ./ ll;
                            end
                            clear lor2
                        end
                        if is_transposed
                            Summ = full(sum(A,2));
                        else
                            Summ = full(sum(A,1))';
                        end
                        Summ(Summ == 0) = options.epps;
                        % Compute OSEM
                        if options.osem || options.ecosem || options.attenuation_phase
                            if verbose
                                tStart = tic;
                            end
                            im_vectors.OSEM_apu = OSEM_im(im_vectors.OSEM_apu, A, epps, uu, Summ, SinD, is_transposed);
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
                            im_vectors.MRAMLA_apu = MBSREM(im_vectors.MRAMLA_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, SinD, randoms_correction, ...
                                is_transposed);
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
                            im_vectors.RAMLA_apu = BSREM_subiter(im_vectors.RAMLA_apu, options.lam, epps, iter, A, uu, SinD, is_transposed);
                            if any(im_vectors.RAMLA_apu < 0)
                                error('Negative values in RAMLA, lower lambda value!')
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
                            im_vectors.ROSEM_apu = ROSEM_subiter(im_vectors.ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed);
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
                            im_vectors.RBI_apu = RBI_subiter(im_vectors.RBI_apu, A, uu, epps, Summ, 0, 0, D, SinD, is_transposed);
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
                            im_vectors.DRAMA_apu = DRAMA_subiter(im_vectors.DRAMA_apu, options.lam_drama, epps, iter, Summ, osa_iter, A, uu, SinD, is_transposed);
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
                            [im_vectors.COSEM_apu, C_co] = COSEM_im(im_vectors.COSEM_apu, A, epps, uu, C_co, D, osa_iter, SinD, is_transposed);
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
                            [im_vectors.ACOSEM_apu, C_aco] = ACOSEM_im(im_vectors.ACOSEM_apu, A, epps, uu, C_aco, D, options.h, osa_iter, SinD, is_transposed);
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
                            im_vectors.MRP_OSL_apu = OSL_OSEM(im_vectors.MRP_OSL_apu, Summ, options.beta_mrp_osem, med, epps, A, uu, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute MBSREM with MRP
                        if options.MRP && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            med = MRP(im_vectors.MRP_MBSREM_apu, options.medx, options.medy, options.medz, Nx, Ny, Nz, epps, options.tr_offsets, options.med_no_norm);
                            im_vectors.MRP_MBSREM_apu = MBSREM(im_vectors.MRP_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, SinD, randoms_correction, is_transposed, ...
                                options.beta_mrp_mbsrem, med);
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
                                im_vectors.MRP_BSREM_apu = BSREM_subiter(im_vectors.MRP_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed);
                            end
                            if any(im_vectors.MRP_BSREM_apu < 0)
                                error('Negative values in BSREM, lower lambda value!')
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
                                im_vectors.MRP_ROSEM_apu = ROSEM_subiter(im_vectors.MRP_ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute RBI-MAP with MRP
                        if options.MRP && options.RBI_MAP
                            if verbose
                                tStart = tic;
                            end
                            med = MRP(im_vectors.MRP_RBI_apu, options.medx, options.medy, options.medz, Nx, Ny, Nz, epps, options.tr_offsets, options.med_no_norm);
                            im_vectors.MRP_RBI_apu = RBI_subiter(im_vectors.MRP_RBI_apu, A, uu, epps, Summ, options.beta_mrp_rbi, med, D, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-MAP MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute COSEM-MAP with MRP
                        if options.MRP && any(options.COSEM_MAP)
                            if verbose
                                tStart = tic;
                            end
                            med = MRP(im_vectors.MRP_COSEM_apu, options.medx, options.medy, options.medz, Nx, Ny, Nz, epps, options.tr_offsets, options.med_no_norm);
                            if options.COSEM_MAP == 1
                                [im_vectors.MRP_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.MRP_COSEM_apu, D, options.beta_mrp_cosem, med, epps, A, uu, ...
                                    C_osl, options.h, options.COSEM_MAP, osa_iter, SinD, is_transposed);
                            else
                                [im_vectors.MRP_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.MRP_COSEM_apu, D, options.beta_mrp_cosem, med, epps, A, uu, ...
                                    C_osl, 0, options.COSEM_MAP, osa_iter, SinD, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-MAP MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute OSL with Quadratic prior
                        if options.quad && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            med = Quadratic_prior(im_vectors.Quad_OSL_apu, options.tr_offsets, options.weights, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                            im_vectors.Quad_OSL_apu = OSL_OSEM(im_vectors.Quad_OSL_apu, Summ, options.beta_quad_osem, med, epps, A, uu, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute MBSREM with Quadratic prior
                        if options.quad && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            med = Quadratic_prior(im_vectors.Quad_MBSREM_apu, options.tr_offsets, options.weights, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                            im_vectors.Quad_MBSREM_apu = MBSREM(im_vectors.Quad_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, SinD, randoms_correction, is_transposed, ...
                                options.beta_quad_mbsrem, med);
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
                                im_vectors.Quad_BSREM_apu = BSREM_subiter(im_vectors.Quad_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed);
                            end
                            if any(im_vectors.Quad_BSREM_apu < 0)
                                error('Negative values in BSREM, lower lambda value')
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
                                im_vectors.Quad_ROSEM_apu = ROSEM_subiter(im_vectors.Quad_ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute RBI-MAP with Quadratic prior
                        if options.quad && options.RBI_MAP
                            if verbose
                                tStart = tic;
                            end
                            med = Quadratic_prior(im_vectors.Quad_RBI_apu, options.tr_offsets, options.weights, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                            im_vectors.Quad_RBI_apu = RBI_subiter(im_vectors.Quad_RBI_apu, A, uu, epps, Summ, options.beta_quad_rbi, med, D, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-MAP Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute COSEM-MAP with Quadratic prior
                        if options.quad && any(options.COSEM_MAP)
                            if verbose
                                tStart = tic;
                            end
                            med = Quadratic_prior(im_vectors.Quad_COSEM_apu, options.tr_offsets, options.weights, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                            if options.COSEM_MAP == 1
                                [im_vectors.Quad_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.Quad_COSEM_apu, D, options.beta_quad_cosem, med, epps, A, uu, ...
                                    C_osl, options.h, options.COSEM_MAP, osa_iter, SinD, is_transposed);
                            else
                                [im_vectors.Quad_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.Quad_COSEM_apu, D, options.beta_quad_cosem, med, epps, A, uu, ...
                                    C_osl, 0, options.COSEM_MAP, osa_iter, SinD, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-MAP Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute OSL with L-filter prior
                        if options.L && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            med = L_filter(im_vectors.L_OSL_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            im_vectors.L_OSL_apu = OSL_OSEM(im_vectors.L_OSL_apu, Summ, options.beta_L_osem, med, epps, A, uu, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute MBSREM with L-filter prior
                        if options.L && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            med = L_filter(im_vectors.L_MBSREM_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            im_vectors.L_MBSREM_apu = MBSREM(im_vectors.L_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, SinD, randoms_correction, is_transposed, ...
                                options.beta_L_mbsrem, med);
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
                                im_vectors.L_BSREM_apu = BSREM_subiter(im_vectors.L_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed);
                            end
                            if any(im_vectors.L_BSREM_apu < 0)
                                error('Negative values in BSREM, lower lambda value')
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
                                im_vectors.L_ROSEM_apu = ROSEM_subiter(im_vectors.L_ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute RBI-MAP with L-filter prior
                        if options.L && options.RBI_MAP
                            if verbose
                                tStart = tic;
                            end
                            med = L_filter(im_vectors.L_RBI_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            im_vectors.L_RBI_apu = RBI_subiter(im_vectors.L_RBI_apu, A, uu, epps, Summ, options.beta_L_rbi, med, D, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-MAP L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute COSEM-MAP with L-filter prior
                        if options.L && any(options.COSEM_MAP)
                            if verbose
                                tStart = tic;
                            end
                            med = L_filter(im_vectors.L_COSEM_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            if options.COSEM_MAP == 1
                                [im_vectors.L_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.L_COSEM_apu, D, options.beta_L_cosem, med, epps, A, uu, C_osl, ...
                                    options.h, options.COSEM_MAP, osa_iter, SinD, is_transposed);
                            else
                                [im_vectors.L_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.L_COSEM_apu, D, options.beta_L_cosem, med, epps, A, uu, C_osl, 0, ...
                                    options.COSEM_MAP, osa_iter, SinD, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-MAP L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute OSL with FMH prior
                        if options.FMH && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            med = FMH(im_vectors.FMH_OSL_apu, options.tr_offsets, options.fmh_weights, options.weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            im_vectors.FMH_OSL_apu = OSL_OSEM(im_vectors.FMH_OSL_apu, Summ, options.beta_fmh_osem, med, epps, A, uu, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.FMH && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            med = FMH(im_vectors.FMH_MBSREM_apu, options.tr_offsets, options.fmh_weights, options.weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            im_vectors.FMH_MBSREM_apu = MBSREM(im_vectors.FMH_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, SinD, randoms_correction, is_transposed, ...
                                options.beta_fmh_mbsrem, med);
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
                                im_vectors.FMH_BSREM_apu = BSREM_subiter(im_vectors.FMH_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed);
                            end
                            if any(im_vectors.FMH_BSREM_apu < 0)
                                error('Negative values in BSREM, lower lambda value')
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
                                im_vectors.FMH_ROSEM_apu = ROSEM_subiter(im_vectors.FMH_ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.FMH && options.RBI_MAP
                            if verbose
                                tStart = tic;
                            end
                            med = FMH(im_vectors.FMH_RBI_apu, options.tr_offsets, options.fmh_weights, options.weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            im_vectors.FMH_RBI_apu = RBI_subiter(im_vectors.FMH_RBI_apu, A, uu, epps, Summ, options.beta_fmh_rbi, med, D, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-MAP FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.FMH && any(options.COSEM_MAP)
                            if verbose
                                tStart = tic;
                            end
                            med = FMH(im_vectors.FMH_COSEM_apu, options.tr_offsets, options.fmh_weights, options.weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            if options.COSEM_MAP == 1
                                [im_vectors.FMH_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.FMH_COSEM_apu, D, options.beta_fmh_cosem, med, epps, A, uu, ...
                                    C_osl, options.h, options.COSEM_MAP, osa_iter, SinD, is_transposed);
                            else
                                [im_vectors.FMH_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.FMH_COSEM_apu, D, options.beta_fmh_cosem, med, epps, A, uu, ...
                                    C_osl, 0, options.COSEM_MAP, osa_iter, SinD, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-MAP FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute OSL with weighted mean prior
                        if options.weighted_mean && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            med = Weighted_mean(im_vectors.Weighted_OSL_apu, options.tr_offsets, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.w_sum, options.med_no_norm);
                            im_vectors.Weighted_OSL_apu = OSL_OSEM(im_vectors.Weighted_OSL_apu, Summ, options.beta_weighted_osem, med, epps, A, uu, SinD, ...
                                is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.weighted_mean && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            med = Weighted_mean(im_vectors.Weighted_MBSREM_apu, options.tr_offsets, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.w_sum, options.med_no_norm);
                            im_vectors.Weighted_MBSREM_apu = MBSREM(im_vectors.Weighted_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, ...
                                options.lam_mbsrem, iter, SinD, randoms_correction, is_transposed, options.beta_weighted_mbsrem, med);
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
                                im_vectors.Weighted_BSREM_apu = BSREM_subiter(im_vectors.Weighted_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed);
                            end
                            if any(im_vectors.Weighted_BSREM_apu < 0)
                                error('Negative values in BSREM, lower lambda value')
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
                                    is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.weighted_mean && options.RBI_MAP
                            if verbose
                                tStart = tic;
                            end
                            med = Weighted_mean(im_vectors.Weighted_RBI_apu, options.tr_offsets, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.w_sum, options.med_no_norm);
                            im_vectors.Weighted_RBI_apu = RBI_subiter(im_vectors.Weighted_RBI_apu, A, uu, epps, Summ, options.beta_weighted_rbi, ...
                                med, D, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-MAP weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.weighted_mean && any(options.COSEM_MAP)
                            if verbose
                                tStart = tic;
                            end
                            med = Weighted_mean(im_vectors.Weighted_COSEM_apu, options.tr_offsets, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.w_sum, options.med_no_norm);
                            if options.COSEM_MAP == 1
                                [im_vectors.Weighted_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.Weighted_COSEM_apu, D, options.beta_weighted_cosem, ...
                                    med, epps, A, uu, C_osl, options.h, options.COSEM_MAP, osa_iter, SinD, is_transposed);
                            else
                                [im_vectors.Weighted_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.Weighted_COSEM_apu, D, options.beta_weighted_cosem, ...
                                    med, epps, A, uu, C_osl, 0, options.COSEM_MAP, osa_iter, SinD, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-MAP weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute OSL with TV prior
                        if options.TV && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            grad = TVpriorFinal(im_vectors.TV_OSL_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, ...
                                options.tr_offsets);
                            im_vectors.TV_OSL_apu = OSL_OSEM(im_vectors.TV_OSL_apu, Summ, options.beta_TV_osem, grad, epps, A, uu, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.TV && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            grad = TVpriorFinal(im_vectors.TV_MBSREM_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, ...
                                options.tr_offsets);
                            im_vectors.TV_MBSREM_apu = MBSREM(im_vectors.TV_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                                iter, SinD, randoms_correction, is_transposed, options.beta_TV_mbsrem, grad);
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
                                im_vectors.TV_BSREM_apu = BSREM_subiter(im_vectors.TV_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed);
                            end
                            if any(im_vectors.TV_BSREM_apu < 0)
                                error('Negative values in BSREM, lower lambda value')
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
                                im_vectors.TV_ROSEM_apu = ROSEM_subiter(im_vectors.TV_ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.TV && options.RBI_MAP
                            if verbose
                                tStart = tic;
                            end
                            grad = TVpriorFinal(im_vectors.TV_RBI_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, ...
                                options.tr_offsets);
                            im_vectors.TV_RBI_apu = RBI_subiter(im_vectors.TV_RBI_apu, A, uu, epps, Summ, options.beta_TV_rbi, grad, D, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-MAP TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.TV && any(options.COSEM_MAP)
                            if verbose
                                tStart = tic;
                            end
                            grad = TVpriorFinal(im_vectors.TV_COSEM_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, ...
                                options.tr_offsets);
                            if options.COSEM_MAP == 1
                                [im_vectors.TV_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.TV_COSEM_apu, D, options.beta_TV_cosem, grad, epps, A, uu, ...
                                    C_osl, options.h, options.COSEM_MAP, osa_iter, SinD, is_transposed);
                            else
                                [im_vectors.TV_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.TV_COSEM_apu, D, options.beta_TV_cosem, grad, epps, A, uu, ...
                                    C_osl, 0, options.COSEM_MAP, osa_iter, SinD, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-MAP TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute OSL with MRP-AD prior
                        if options.AD && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            if osa_iter > 1
                                med = AD(im_vectors.AD_OSL_apu, options.FluxType, Nx, Ny, Nz, options);
                                im_vectors.AD_OSL_apu = OSL_OSEM(im_vectors.AD_OSL_apu, Summ, options.beta_ad_osem, med, epps, A, uu, SinD, is_transposed);
                            else
                                im_vectors.AD_OSL_apu = OSEM_im(im_vectors.AD_OSL_apu, A, epps, uu, Summ);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.AD && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            med = AD(im_vectors.AD_MBSREM_apu, options.FluxType, Nx, Ny, Nz, options);
                            im_vectors.AD_MBSREM_apu = MBSREM(im_vectors.AD_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                                iter, SinD, randoms_correction, is_transposed, options.beta_ad_mbsrem, med);
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
                                im_vectors.AD_BSREM_apu = BSREM_subiter(im_vectors.AD_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed);
                            end
                            if any(im_vectors.AD_BSREM_apu < 0)
                                error('Negative values in BSREM, lower lambda value')
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
                                im_vectors.AD_ROSEM_apu = ROSEM_subiter(im_vectors.AD_ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.AD && options.RBI_MAP
                            if verbose
                                tStart = tic;
                            end
                            med = AD(im_vectors.AD_RBI_apu, options.FluxType, Nx, Ny, Nz, options);
                            im_vectors.AD_RBI_apu = RBI_subiter(im_vectors.AD_RBI_apu, A, uu, epps, Summ, options.beta_ad_rbi, med, D, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-MAP AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.AD && any(options.COSEM_MAP)
                            if verbose
                                tStart = tic;
                            end
                            med = AD(im_vectors.AD_COSEM_apu, options.FluxType, Nx, Ny, Nz, options);
                            if options.COSEM_MAP == 1
                                [im_vectors.AD_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.AD_COSEM_apu, D, options.beta_ad_cosem, med, epps, A, uu, ...
                                    C_osl, options.h, options.COSEM_MAP, osa_iter, SinD, is_transposed);
                            else
                                [im_vectors.AD_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.AD_COSEM_apu, D, options.beta_ad_cosem, med, epps, A, uu, ...
                                    C_osl, 0, options.COSEM_MAP, osa_iter, SinD, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-MAP AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute OSL with APLS prior
                        if options.APLS && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            grad = TVpriorFinal(im_vectors.APLS_OSL_apu, [], Nx, Ny, Nz, true, options, 4);
                            im_vectors.APLS_OSL_apu = OSL_OSEM(im_vectors.APLS_OSL_apu, Summ, options.beta_APLS_osem, grad, epps, A, uu, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.APLS && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            grad = TVpriorFinal(im_vectors.APLS_MBSREM_apu, [], Nx, Ny, Nz, true, options, 4);
                            im_vectors.APLS_MBSREM_apu = MBSREM(im_vectors.APLS_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                                iter, SinD, randoms_correction, is_transposed, options.beta_APLS_mbsrem, grad);
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
                                im_vectors.APLS_BSREM_apu = BSREM_subiter(im_vectors.APLS_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed);
                            end
                            if any(im_vectors.APLS_BSREM_apu < 0)
                                error('Negative values in BSREM, lower lambda value')
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
                                im_vectors.APLS_ROSEM_apu = ROSEM_subiter(im_vectors.APLS_ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.APLS && options.RBI_MAP
                            if verbose
                                tStart = tic;
                            end
                            grad = TVpriorFinal(im_vectors.APLS_RBI_apu, [], Nx, Ny, Nz, true, options, 4);
                            im_vectors.APLS_RBI_apu = RBI_subiter(im_vectors.APLS_RBI_apu, A, uu, epps, Summ, SinD, options.beta_APLS_rbi, grad, D, ...
                                is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-MAP APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.APLS && any(options.COSEM_MAP)
                            if verbose
                                tStart = tic;
                            end
                            grad = TVpriorFinal(im_vectors.APLS_COSEM_apu, [], Nx, Ny, Nz, true, options, 4);
                            if options.COSEM_MAP == 1
                                [im_vectors.APLS_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.APLS_COSEM_apu, D, options.beta_APLS_cosem, grad, A, uu, ...
                                    epps, C_osl, options.h, options.COSEM_MAP, osa_iter, SinD, is_transposed);
                            else
                                [im_vectors.APLS_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.APLS_COSEM_apu, D, options.beta_APLS_cosem, grad, A, uu, ...
                                    epps, C_osl, 0, options.COSEM_MAP, osa_iter, SinD, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-MAP APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute OSL with TGV prior
                        if options.TGV && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            grad = TGV(im_vectors.TGV_OSL_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                            im_vectors.TGV_OSL_apu = OSL_OSEM(im_vectors.TGV_OSL_apu, Summ, options.beta_TGV_osem, grad, epps, A, uu, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.TGV && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            grad = TGV(im_vectors.TGV_MBSREM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                            im_vectors.TGV_MBSREM_apu = MBSREM(im_vectors.TGV_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                                iter, SinD, randoms_correction, is_transposed, options.beta_TGV_mbsrem, grad);
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
                                im_vectors.TGV_BSREM_apu = BSREM_subiter(im_vectors.TGV_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed);
                            end
                            if any(im_vectors.TGV_BSREM_apu < 0)
                                error('Negative values in BSREM, lower lambda value')
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
                                im_vectors.TGV_ROSEM_apu = ROSEM_subiter(im_vectors.TGV_ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.TGV && options.RBI_MAP
                            if verbose
                                tStart = tic;
                            end
                            grad = TGV(im_vectors.TGV_RBI_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                            im_vectors.TGV_RBI_apu = RBI_subiter(im_vectors.TGV_RBI_apu, A, uu, epps, Summ, options.beta_TGV_rbi, grad, D, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-MAP TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.TGV && any(options.COSEM_MAP)
                            if verbose
                                tStart = tic;
                            end
                            grad = TGV(im_vectors.TGV_COSEM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                            if options.COSEM_MAP == 1
                                [im_vectors.TGV_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.TGV_COSEM_apu, D, options.beta_TGV_cosem, grad, A, uu, ...
                                    epps, C_osl, options.h, options.COSEM_MAP, osa_iter, SinD, is_transposed);
                            else
                                [im_vectors.TGV_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.TGV_COSEM_apu, D, options.beta_TGV_cosem, grad, A, uu, ...
                                    epps, C_osl, 0, options.COSEM_MAP, osa_iter, SinD, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-MAP TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        % Compute OSL with NLM prior
                        if options.NLM && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            med = NLM(im_vectors.NLM_OSL_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                options.sigma, epps, Nx, Ny, Nz, options);
                            im_vectors.NLM_OSL_apu = OSL_OSEM(im_vectors.NLM_OSL_apu, Summ, options.beta_NLM_osem, med, epps, A, uu, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.NLM && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            med = NLM(im_vectors.NLM_MBSREM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                options.sigma, epps, Nx, Ny, Nz, options);
                            im_vectors.NLM_MBSREM_apu = MBSREM(im_vectors.NLM_MBSREM_apu, options.U, options.pj3, A, epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                                iter, SinD, randoms_correction, is_transposed, options.beta_NLM_mbsrem, med);
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
                                im_vectors.NLM_BSREM_apu = BSREM_subiter(im_vectors.NLM_BSREM_apu, options.lam, epps, iter, A, uu, SinD, is_transposed);
                            end
                            if any(im_vectors.NLM_BSREM_apu < 0)
                                error('Negative values in BSREM, lower lambda value')
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
                                im_vectors.NLM_ROSEM_apu = ROSEM_subiter(im_vectors.NLM_ROSEM_apu, options.lam_rosem, iter, Summ, epps, A, uu, SinD, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.NLM && options.RBI_MAP
                            if verbose
                                tStart = tic;
                            end
                            med = NLM(im_vectors.NLM_RBI_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                options.sigma, epps, Nx, Ny, Nz, options);
                            im_vectors.NLM_RBI_apu = RBI_subiter(im_vectors.NLM_RBI_apu, A, uu, epps, Summ, options.beta_NLM_rbi, med, D, SinD, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-MAP NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.NLM && any(options.COSEM_MAP)
                            if verbose
                                tStart = tic;
                            end
                            med = NLM(im_vectors.NLM_COSEM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                options.sigma, epps, Nx, Ny, Nz, options);
                            if options.COSEM_MAP == 1
                                [im_vectors.NLM_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.NLM_COSEM_apu, D, options.beta_NLM_cosem, med, A, uu, ...
                                    epps, C_osl, options.h, options.COSEM_MAP, osa_iter, SinD, is_transposed);
                            else
                                [im_vectors.NLM_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.NLM_COSEM_apu, D, options.beta_NLM_cosem, med, A, uu, ...
                                    epps, C_osl, 0, options.COSEM_MAP, osa_iter, SinD, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-MAP NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        clear A
                        fn = fieldnames(im_vectors);
                        for kk = 2 : 2 : numel(fn)
                            im_vectors.(fn{kk})(im_vectors.(fn{kk}) < 0) = epps;
                        end
                    end
                    if options.osem
                        im_vectors.OSEM(:, iter + 1) = im_vectors.OSEM_apu;
                    end
                    
                    if options.mramla
                        im_vectors.MRAMLA(:, iter + 1) = im_vectors.MRAMLA_apu;
                    end
                    
                    if options.ramla
                        im_vectors.RAMLA(:, iter + 1) = im_vectors.RAMLA_apu;
                    end
                    
                    if options.rosem
                        im_vectors.ROSEM(:, iter + 1) = im_vectors.ROSEM_apu;
                    end
                    
                    if options.rbi
                        im_vectors.RBI(:, iter + 1) = im_vectors.RBI_apu;
                    end
                    
                    if options.drama
                        im_vectors.DRAMA(:, iter + 1) = im_vectors.DRAMA_apu;
                    end
                    
                    if options.cosem
                        im_vectors.COSEM(:, iter + 1) = im_vectors.COSEM_apu;
                    end
                    
                    if options.ecosem
                        im_vectors.ECOSEM(:, iter + 1) = im_vectors.ECOSEM_apu;
                    end
                    
                    if options.acosem
                        im_vectors.ACOSEM(:, iter + 1) = im_vectors.ACOSEM_apu;
                    end
                    
                    if options.MRP && options.OSL_OSEM
                        im_vectors.MRP_OSL(:, iter + 1) = im_vectors.MRP_OSL_apu;
                    end
                    if options.MRP && options.MBSREM
                        im_vectors.MRP_MBSREM(:, iter + 1) = im_vectors.MRP_MBSREM_apu;
                    end
                    
                    if options.MRP && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        med = MRP(im_vectors.MRP_BSREM_apu, options.medx, options.medy, options.medz, Nx, Ny, Nz, epps, options.tr_offsets, options.med_no_norm);
                        im_vectors.MRP_BSREM(:,iter+1) = BSREM_iter(im_vectors.MRP_BSREM_apu, options.lam, iter, options.beta_mrp_bsrem, med, epps);
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
                        im_vectors.MRP_ROSEM(:,iter+1) = BSREM_iter(im_vectors.MRP_ROSEM_apu, options.lam_rosem, iter, options.beta_mrp_rosem, med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM MRP iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM MRP iteration ' num2str(iter) ' finished'])
                        end
                    end
                    
                    if options.MRP && options.RBI_MAP
                        im_vectors.MRP_RBI(:, iter + 1) = im_vectors.MRP_RBI_apu;
                    end
                    
                    if options.MRP && any(options.COSEM_MAP)
                        im_vectors.MRP_COSEM(:, iter + 1) = im_vectors.MRP_COSEM_apu;
                    end
                    
                    if options.quad && options.OSL_OSEM
                        im_vectors.Quad_OSL(:, iter + 1) = im_vectors.Quad_OSL_apu;
                    end
                    if options.quad && options.MBSREM
                        im_vectors.Quad_MBSREM(:, iter + 1) = im_vectors.Quad_MBSREM_apu;
                    end
                    
                    if options.quad && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        med = Quadratic_prior(im_vectors.Quad_BSREM_apu, options.tr_offsets, options.weights, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                        im_vectors.Quad_BSREM(:,iter+1) = BSREM_iter(im_vectors.Quad_BSREM_apu, options.lam, iter, options.beta_quad_bsrem, med, epps);
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
                        med = Quadratic_prior(im_vectors.Quad_ROSEM_apu, options.tr_offsets, options.weights, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                        im_vectors.Quad_ROSEM(:,iter+1) = BSREM_iter(im_vectors.Quad_ROSEM_apu, options.lam_rosem, iter, options.beta_quad_rosem, med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM quadratic iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM quadratic iteration ' num2str(iter) ' finished'])
                        end
                    end
                    
                    if options.quad && options.RBI_MAP
                        im_vectors.Quad_RBI(:, iter + 1) = im_vectors.Quad_RBI_apu;
                    end
                    
                    if options.quad && any(options.COSEM_MAP)
                        im_vectors.Quad_COSEM(:, iter + 1) = im_vectors.Quad_COSEM_apu;
                    end
                    
                    if options.L && options.OSL_OSEM
                        im_vectors.L_OSL(:, iter + 1) = im_vectors.L_OSL_apu;
                    end
                    if options.L && options.MBSREM
                        im_vectors.L_MBSREM(:, iter + 1) = im_vectors.L_MBSREM_apu;
                    end
                    
                    if options.L && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        med = L_filter(im_vectors.L_BSREM_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                        im_vectors.L_BSREM(:,iter+1) = BSREM_iter(im_vectors.L_BSREM_apu, options.lam, iter, options.beta_L_bsrem, med, epps);
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
                        im_vectors.L_ROSEM(:,iter+1) = BSREM_iter(im_vectors.L_ROSEM_apu, options.lam_rosem, iter, options.beta_L_rosem, med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM L-filter iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM L-filter iteration ' num2str(iter) ' finished'])
                        end
                    end
                    
                    if options.L && options.RBI_MAP
                        im_vectors.L_RBI(:, iter + 1) = im_vectors.L_RBI_apu;
                    end
                    
                    if options.L && any(options.COSEM_MAP)
                        im_vectors.L_COSEM(:, iter + 1) = im_vectors.L_COSEM_apu;
                    end
                    
                    if options.FMH && options.OSL_OSEM
                        im_vectors.FMH_OSL(:, iter + 1) = im_vectors.FMH_OSL_apu;
                    end
                    if options.FMH && options.MBSREM
                        im_vectors.FMH_MBSREM(:, iter + 1) = im_vectors.FMH_MBSREM_apu;
                    end
                    
                    if options.FMH && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        med = FMH(im_vectors.FMH_BSREM_apu, options.tr_offsets, options.fmh_weights, options.weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                            options.med_no_norm);
                        im_vectors.FMH_BSREM(:,iter+1) = BSREM_iter(im_vectors.FMH_BSREM_apu, options.lam, iter, options.beta_fmh_bsrem, med, epps);
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
                        med = FMH(im_vectors.FMH_ROSEM_apu, options.tr_offsets, options.fmh_weights, options.weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                            options.med_no_norm);
                        im_vectors.FMH_ROSEM(:,iter+1) = BSREM_iter(im_vectors.FMH_ROSEM_apu, options.lam_rosem, iter, options.beta_fmh_rosem, med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM FMH iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM FMH iteration ' num2str(iter) ' finished'])
                        end
                    end
                    
                    if options.FMH && options.RBI_MAP
                        im_vectors.FMH_RBI(:, iter + 1) = im_vectors.FMH_RBI_apu;
                    end
                    
                    if options.FMH && any(options.COSEM_MAP)
                        im_vectors.FMH_COSEM(:, iter + 1) = im_vectors.FMH_COSEM_apu;
                    end
                    
                    if options.weighted_mean && options.OSL_OSEM
                        im_vectors.Weighted_OSL(:, iter + 1) = im_vectors.Weighted_OSL_apu;
                    end
                    if options.weighted_mean && options.MBSREM
                        im_vectors.Weighted_MBSREM(:, iter + 1) = im_vectors.Weighted_MBSREM_apu;
                    end
                    
                    if options.weighted_mean && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        med = Weighted_mean(im_vectors.Weighted_BSREM_apu, options.tr_offsets, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                            options.mean_type, epps, options.w_sum, options.med_no_norm);
                        im_vectors.Weighted_BSREM(:,iter+1) = BSREM_iter(im_vectors.Weighted_BSREM_apu, options.lam, iter, options.beta_weighted_bsrem, ...
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
                        med = Weighted_mean(im_vectors.Weighted_ROSEM_apu, options.tr_offsets, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                            options.mean_type, epps, options.w_sum, options.med_no_norm);
                        im_vectors.Weighted_ROSEM(:,iter+1) = BSREM_iter(im_vectors.Weighted_ROSEM_apu, options.lam_rosem, iter, options.beta_weighted_rosem, ...
                            med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM weighted mean iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM weighted mean iteration ' num2str(iter) ' finished'])
                        end
                    end
                    
                    if options.weighted_mean && options.RBI_MAP
                        im_vectors.Weighted_RBI(:, iter + 1) = im_vectors.Weighted_RBI_apu;
                    end
                    
                    if options.weighted_mean && any(options.COSEM_MAP)
                        im_vectors.Weighted_COSEM(:, iter + 1) = im_vectors.Weighted_COSEM_apu;
                    end
                    
                    if options.TV && options.OSL_OSEM
                        im_vectors.TV_OSL(:, iter + 1) = im_vectors.TV_OSL_apu;
                    end
                    if options.TV && options.MBSREM
                        im_vectors.TV_MBSREM(:, iter + 1) = im_vectors.TV_MBSREM_apu;
                    end
                    
                    if options.TV && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        grad = TVpriorFinal(im_vectors.TV_BSREM_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, options.tr_offsets);
                        im_vectors.TV_BSREM(:, iter + 1) = BSREM_iter(im_vectors.TV_BSREM_apu, options.lam, iter, options.beta_TV_bsrem, grad, epps);
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
                        im_vectors.TV_ROSEM(:, iter + 1) = BSREM_iter(im_vectors.TV_ROSEM_apu, options.lam_rosem, iter, options.beta_TV_rosem, grad, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM TV iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM TV iteration ' num2str(iter) ' finished'])
                        end
                    end
                    
                    if options.TV && options.RBI_MAP
                        im_vectors.TV_RBI(:, iter + 1) = im_vectors.TV_RBI_apu;
                    end
                    
                    if options.TV && any(options.COSEM_MAP)
                        im_vectors.TV_COSEM(:, iter + 1) = im_vectors.TV_COSEM_apu;
                    end
                    
                    if options.AD && options.OSL_OSEM
                        im_vectors.AD_OSL(:, iter + 1) = im_vectors.AD_OSL_apu;
                    end
                    if options.AD && options.MBSREM
                        im_vectors.AD_MBSREM(:, iter + 1) = im_vectors.AD_MBSREM_apu;
                    end
                    
                    if options.AD && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        med = AD(im_vectors.AD_BSREM_apu, options.FluxType, Nx, Ny, Nz, options);
                        im_vectors.AD_BSREM(:,iter+1) = BSREM_iter(im_vectors.AD_BSREM_apu, options.lam, iter, options.beta_ad_bsrem, med, epps);
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
                        im_vectors.AD_ROSEM(:, iter + 1) = BSREM_iter(im_vectors.AD_ROSEM_apu, options.lam_rosem, iter, options.beta_ad_rosem, med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM AD iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM AD iteration ' num2str(iter) ' finished'])
                        end
                    end
                    if options.AD && options.RBI_MAP
                        im_vectors.AD_RBI(:, iter + 1) = im_vectors.AD_RBI_apu;
                    end
                    if options.AD && any(options.COSEM_MAP)
                        im_vectors.AD_COSEM(:, iter + 1) = im_vectors.AD_COSEM_apu;
                    end
                    if options.APLS && options.OSL_OSEM
                        im_vectors.APLS_OSL(:, iter + 1) = im_vectors.APLS_OSL_apu;
                    end
                    if options.APLS && options.MBSREM
                        im_vectors.APLS_MBSREM(:, iter + 1) = im_vectors.APLS_MBSREM_apu;
                    end
                    if options.APLS && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        grad = TVpriorFinal(im_vectors.APLS_BSREM_apu, 0, Nx, Ny, Nz, true, options, 4);
                        im_vectors.APLS_BSREM(:, iter + 1) = BSREM_iter(im_vectors.APLS_BSREM_apu, options.lam, iter, options.beta_APLS_bsrem, grad, epps);
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
                        grad = TVpriorFinal(im_vectors.APLS_ROSEM_apu, 0, Nx, Ny, Nz, true, options, 4);
                        im_vectors.APLS_ROSEM(:, iter + 1) = BSREM_iter(im_vectors.APLS_ROSEM_apu, options.lam_rosem, iter, options.beta_APLS_rosem, grad, ...
                            epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM APLS iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM APLS iteration ' num2str(iter) ' finished'])
                        end
                    end
                    if options.APLS && options.RBI_MAP
                        im_vectors.APLS_RBI(:, iter + 1) = im_vectors.APLS_RBI_apu;
                    end
                    if options.APLS && any(options.COSEM_MAP)
                        im_vectors.APLS_COSEM(:, iter + 1) = im_vectors.APLS_COSEM_apu;
                    end
                    if options.TGV && options.OSL_OSEM
                        im_vectors.TGV_OSL(:, iter + 1) = im_vectors.TGV_OSL_apu;
                    end
                    if options.TGV && options.MBSREM
                        im_vectors.TGV_MBSREM(:, iter + 1) = im_vectors.TGV_MBSREM_apu;
                    end
                    if options.TGV && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        grad = TGV(im_vectors.TGV_BSREM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                        im_vectors.TGV_BSREM(:, iter + 1) = BSREM_iter(im_vectors.TGV_BSREM_apu, options.lam, iter, options.beta_TGV_bsrem, grad, epps);
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
                        im_vectors.TGV_ROSEM(:, iter + 1) = BSREM_iter(im_vectors.TGV_ROSEM_apu, options.lam_rosem, iter, options.beta_TGV_rosem, grad, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM TGV iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM TGV iteration ' num2str(iter) ' finished'])
                        end
                    end
                    if options.TGV && options.RBI_MAP
                        im_vectors.TGV_RBI(:, iter + 1) = im_vectors.TGV_RBI_apu;
                    end
                    if options.TGV && any(options.COSEM_MAP)
                        im_vectors.TGV_COSEM(:, iter + 1) = im_vectors.TGV_COSEM_apu;
                    end
                    if options.NLM && options.OSL_OSEM
                        im_vectors.NLM_OSL(:, iter + 1) = im_vectors.NLM_OSL_apu;
                    end
                    if options.NLM && options.MBSREM
                        im_vectors.NLM_MBSREM(:, iter + 1) = im_vectors.NLM_MBSREM_apu;
                    end
                    if options.NLM && options.BSREM
                        if verbose
                            tStart = tic;
                        end
                        med = NLM(im_vectors.NLM_BSREM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                            options.sigma, epps, Nx, Ny, Nz, options);
                        im_vectors.NLM_BSREM(:, iter + 1) = BSREM_iter(im_vectors.NLM_BSREM_apu, options.lam_rosem, iter, options.beta_NLM_bsrem, med, epps);
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
                        im_vectors.NLM_ROSEM(:, iter + 1) = BSREM_iter(im_vectors.NLM_ROSEM_apu, options.lam_rosem, iter, options.beta_NLM_rosem, med, epps);
                        if verbose
                            tElapsed = toc(tStart);
                            disp(['ROSEM NLM iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['ROSEM NLM iteration ' num2str(iter) ' finished'])
                        end
                    end
                    if options.NLM && options.RBI_MAP
                        im_vectors.NLM_RBI(:, iter + 1) = im_vectors.NLM_RBI_apu;
                    end
                    if options.NLM && any(options.COSEM_MAP)
                        im_vectors.NLM_COSEM(:, iter + 1) = im_vectors.NLM_COSEM_apu;
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
                % Supports only MLEM, OSEM, RAMLA and ROSEM (and their MAP
                % versions)
            elseif options.implementation == 4
                if llo == 1
                    no_norm = false;
                end
                if ~use_raw_data
                    if isempty(pseudot)
                        pseudot = uint32(0);
                    end
                end
                options.n_rays = uint16(options.n_rays);
                if options.n_rays > 1 && ~options.precompute_lor && options.projector_type == 1
                    if options.n_rays < 4 || options.n_rays > 5
                        error('Only ray counts of 1, 4 and 5 are supported')
                    end
                    if options.use_raw_data
                        [x, y] = detector_coordinates_multiray(options);
                        x = x(:,[2 1 3]);
                        y = y(:,[2 1 3]);
                    else
                        [x, y] = sinogram_coordinates_2D_multiray(options);
                        x = x(:,:,[2 1 3]);
                        y = y(:,:,[2 1 3]);
                    end
                    x = x(:);
                    y = y(:);
                end
                dc_z = z_det(2,1) - z_det(1,1);
                for iter = 1 : Niter
                    if OS_bool
                        if verbose
                            tStart_iter = tic;
                        end
                        if iter == 1 && llo == 1
                            f_Summ = zeros(Nx*Ny*Nz,subsets);
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
                            if verbose
                                tStart = tic;
                            end
                            uu = double(full(Sino(pituus(osa_iter)+1:pituus(osa_iter + 1))));
                            if use_raw_data
                                L_input = LL(pituus(osa_iter) * 2 + 1 : pituus(osa_iter + 1) * 2);
                                xy_index_input = uint32(0);
                                z_index_input = uint32(0);
                            else
                                L_input = uint16(0);
                                xy_index_input = xy_index(pituus(osa_iter)+1:pituus(osa_iter + 1));
                                z_index_input = z_index(pituus(osa_iter)+1:pituus(osa_iter + 1));
                            end
                            if options.precompute_lor
                                lor_a_input = lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1));
                            else
                                lor_a_input = uint16(0);
                            end
                            
                            [Summ, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                                options.normalization, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                lor_a_input, xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, options.verbose, (use_raw_data), uint32(1), epps, uu, ...
                                im_vectors.OSEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, options.tube_width_xy, x_center, y_center, z_center, ...
                                options.tube_width_z, int32(options.accuracy_factor), options.n_rays, dc_z);
                            
                            if iter == 1 && llo == 1
                                f_Summ(:,osa_iter) = Summ + epps;
                            end
                            if options.osem
                                im_vectors.OSEM_apu = OSEM_im(im_vectors.OSEM_apu, rhs, f_Summ(:,osa_iter), epps);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.ramla
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, f_Summ(:,osa_iter), rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error(['Negative values in RAMLA, lower lambda value! lambda <= ' num2str(min(1./f_Summ(:,osa_iter)))])
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['RAMLA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.rosem
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, options.lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.MRP && options.OSL_OSEM
                                med = MRP(im_vectors.OSEM_apu, options.medx, options.medy, options.medz, Nx, Ny, Nz, epps, options.tr_offsets, options.med_no_norm);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_mrp_osem, med, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.MRP && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error('Negative values in BSREM, lower lambda value!')
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
                            elseif options.quad && options.OSL_OSEM
                                med = Quadratic_prior(im_vectors.OSEM_apu, options.tr_offsets, options.weights, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_quad_osem, med, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSL Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.quad && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error('Negative values in BSREM, lower lambda value')
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
                            elseif options.L && options.OSL_OSEM
                                med = L_filter(im_vectors.OSEM_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_L_osem, med, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSL L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.L && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error('Negative values in BSREM, lower lambda value')
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
                            elseif options.FMH && options.OSL_OSEM
                                if verbose
                                    tStart = tic;
                                end
                                med = FMH(im_vectors.OSEM_apu, options.tr_offsets, options.fmh_weights, options.weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                    options.med_no_norm);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_fmh_osem, med, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSL FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.FMH && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error('Negative values in BSREM, lower lambda value')
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
                            elseif options.weighted_mean && options.OSL_OSEM
                                med = Weighted_mean(im_vectors.OSEM_apu, options.tr_offsets, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                    options.mean_type, epps, options.w_sum, options.med_no_norm);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_weighted_osem, med, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSL weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.weighted_mean && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error('Negative values in BSREM, lower lambda value')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.weighted_mean && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, options.lam_rosem, epps, iter, f_Summ(:,osa_iter), rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM-MAP weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.TV && options.OSL_OSEM
                                grad = TVpriorFinal(im_vectors.OSEM_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, ...
                                    options.tr_offsets);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_TV_osem, grad, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSL TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.TV && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error('Negative values in BSREM, lower lambda value')
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
                            elseif options.AD && options.OSL_OSEM
                                if osa_iter > 1
                                    med = AD(im_vectors.OSEM_apu, options.FluxType, Nx, Ny, Nz, options);
                                    im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_ad_osem, med, epps, rhs);
                                else
                                    im_vectors.OSEM_apu = OSEM_im(im_vectors.OSEM_apu, rhs, f_Summ(:,osa_iter), epps);
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSL AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.AD && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error('Negative values in BSREM, lower lambda value')
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
                            elseif options.APLS && options.OSL_OSEM
                                grad = TVpriorFinal(im_vectors.OSEM_apu, [], Nx, Ny, Nz, true, options, 4);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_APLS_osem, grad, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSL APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.APLS && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error('Negative values in BSREM, lower lambda value')
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
                            elseif options.TGV && options.OSL_OSEM
                                grad = TGV(im_vectors.OSEM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_TGV_osem, grad, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSL TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.TGV && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error('Negative values in BSREM, lower lambda value')
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
                            elseif options.NLM && options.OSL_OSEM
                                if verbose
                                    tStart = tic;
                                end
                                med = NLM(im_vectors.OSEM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                    options.sigma, epps, Nx, Ny, Nz, options);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_NLM_osem, med, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSL NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.NLM && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, options.lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error('Negative values in BSREM, lower lambda value')
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
                            end
                            
                            clear Summ rhs
                        end
                        if options.osem
                            im_vectors.OSEM(:,iter + 1) = im_vectors.OSEM_apu;
                        elseif options.ramla
                            im_vectors.RAMLA(:, iter + 1) = im_vectors.OSEM_apu;
                        elseif options.rosem
                            im_vectors.ROSEM(:, iter + 1) = im_vectors.OSEM_apu;
                        elseif options.MRP && options.OSL_OSEM
                            im_vectors.MRP_OSL(:, iter + 1) = im_vectors.OSEM_apu;
                        elseif options.MRP && options.BSREM
                            med = MRP(im_vectors.OSEM_apu, options.medx, options.medy, options.medz, Nx, Ny, Nz, epps, options.tr_offsets, options.med_no_norm);
                            im_vectors.MRP_BSREM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, options.lam, iter, SinD, randoms_correction, is_transposed, options.beta_mrp_bsrem, med, epps);
                        elseif options.MRP && options.ROSEM_MAP
                            med = MRP(im_vectors.OSEM_apu, options.medx, options.medy, options.medz, Nx, Ny, Nz, epps, options.tr_offsets, options.med_no_norm);
                            im_vectors.MRP_ROSEM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, options.beta_mrp_rosem, med, epps);
                        elseif options.quad && options.OSL_OSEM
                            im_vectors.Quad_OSL(:, iter + 1) = im_vectors.OSEM_apu;
                        elseif options.quad && options.BSREM
                            med = Quadratic_prior(im_vectors.OSEM_apu, options.tr_offsets, options.weights, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                            im_vectors.Quad_BSREM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, options.lam, iter, SinD, randoms_correction, is_transposed, options.beta_quad_bsrem, med, epps);
                        elseif options.quad && options.ROSEM_MAP
                            med = Quadratic_prior(im_vectors.OSEM_apu, options.tr_offsets, options.weights, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                            im_vectors.Quad_ROSEM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, options.beta_quad_rosem, med, epps);
                        elseif options.L && options.OSL_OSEM
                            im_vectors.L_OSL(:, iter + 1) = im_vectors.L_OSL_apu;
                        elseif options.L && options.BSREM
                            med = L_filter(im_vectors.OSEM_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            im_vectors.L_BSREM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, options.lam, iter, SinD, randoms_correction, is_transposed, options.beta_L_bsrem, med, epps);
                        elseif options.L && options.ROSEM_MAP
                            med = L_filter(im_vectors.OSEM_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            im_vectors.L_ROSEM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, options.beta_L_rosem, med, epps);
                        elseif options.FMH && options.OSL_OSEM
                            im_vectors.FMH_OSL(:, iter + 1) = im_vectors.OSEM_apu;
                        elseif options.FMH && options.BSREM
                            med = FMH(im_vectors.OSEM_apu, options.tr_offsets, options.fmh_weights, options.weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            im_vectors.FMH_BSREM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, options.lam, iter, SinD, randoms_correction, is_transposed, options.beta_fmh_bsrem, med, epps);
                        elseif options.FMH && options.ROSEM_MAP
                            med = FMH(im_vectors.OSEM_apu, options.tr_offsets, options.fmh_weights, options.weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            im_vectors.FMH_ROSEM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, options.beta_fmh_rosem, med, epps);
                        elseif options.weighted_mean && options.OSL_OSEM
                            im_vectors.Weighted_OSL(:, iter + 1) = im_vectors.Weighted_OSL_apu;
                        elseif options.weighted_mean && options.BSREM
                            med = Weighted_mean(im_vectors.OSEM_apu, options.tr_offsets, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.w_sum, options.med_no_norm);
                            im_vectors.Weighted_BSREM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, options.lam, iter, SinD, randoms_correction, is_transposed, options.beta_weighted_bsrem, ...
                                med, epps);
                        elseif options.weighted_mean && options.ROSEM_MAP
                            med = Weighted_mean(im_vectors.OSEM_apu, options.tr_offsets, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.w_sum, options.med_no_norm);
                            im_vectors.Weighted_ROSEM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, options.beta_weighted_rosem, ...
                                med, epps);
                        elseif options.TV && options.OSL_OSEM
                            im_vectors.TV_OSL(:, iter + 1) = im_vectors.OSEM_apu;
                        elseif options.TV && options.BSREM
                            grad = TVpriorFinal(im_vectors.OSEM_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, options.tr_offsets);
                            im_vectors.TV_BSREM(:, iter + 1) = BSREM_iter(im_vectors.OSEM_apu, options.lam, iter, SinD, randoms_correction, is_transposed, options.beta_TV_bsrem, grad, epps);
                        elseif options.TV && options.ROSEM_MAP
                            grad = TVpriorFinal(im_vectors.OSEM_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, options.tr_offsets);
                            im_vectors.TV_ROSEM(:, iter + 1) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, options.beta_TV_rosem, grad, epps);
                        elseif options.AD && options.OSL_OSEM
                            im_vectors.AD_OSL(:, iter + 1) = im_vectors.AD_OSL_apu;
                        elseif options.AD && options.BSREM
                            med = AD(im_vectors.OSEM_apu, options.FluxType, Nx, Ny, Nz, options);
                            im_vectors.AD_BSREM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, options.lam, iter, SinD, randoms_correction, is_transposed, options.beta_ad_bsrem, med, epps);
                        elseif options.AD && options.ROSEM_MAP
                            med = AD(im_vectors.OSEM_apu, options.FluxType, Nx, Ny, Nz, options);
                            im_vectors.AD_ROSEM(:, iter + 1) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, options.beta_ad_rosem, med, epps);
                        elseif options.APLS && options.OSL_OSEM
                            im_vectors.APLS_OSL(:, iter + 1) = im_vectors.APLS_OSL_apu;
                        elseif options.APLS && options.BSREM
                            grad = TVpriorFinal(im_vectors.OSEM_apu, 0, Nx, Ny, Nz, true, options, 4);
                            im_vectors.APLS_BSREM(:, iter + 1) = BSREM_iter(im_vectors.OSEM_apu, options.lam, iter, SinD, randoms_correction, is_transposed, options.beta_APLS_bsrem, grad, epps);
                        elseif options.APLS && options.ROSEM_MAP
                            grad = TVpriorFinal(im_vectors.OSEM_apu, 0, Nx, Ny, Nz, true, options, 4);
                            im_vectors.APLS_ROSEM(:, iter + 1) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, options.beta_APLS_rosem, grad, ...
                                epps);
                        elseif options.TGV && options.OSL_OSEM
                            im_vectors.TGV_OSL(:, iter + 1) = im_vectors.OSEM_apu;
                        elseif options.TGV && options.BSREM
                            grad = TGV(im_vectors.OSEM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                            im_vectors.TGV_BSREM(:, iter + 1) = BSREM_iter(im_vectors.OSEM_apu, options.lam, iter, SinD, randoms_correction, is_transposed, options.beta_TGV_bsrem, grad, epps);
                        elseif options.TGV && options.ROSEM_MAP
                            grad = TGV(im_vectors.OSEM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                            im_vectors.TGV_ROSEM(:, iter + 1) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, options.beta_TGV_rosem, grad, epps);
                        elseif options.NLM && options.OSL_OSEM
                            im_vectors.NLM_OSL(:, iter + 1) = im_vectors.NLM_OSL_apu;
                        elseif options.NLM && options.BSREM
                            med = NLM(im_vectors.OSEM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                options.sigma, epps, Nx, Ny, Nz, options);
                            im_vectors.NLM_BSREM(:, iter + 1) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, options.beta_NLM_bsrem, med, epps);
                        elseif options.NLM && options.ROSEM_MAP
                            med = NLM(im_vectors.OSEM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                options.sigma, epps, Nx, Ny, Nz, options);
                            im_vectors.NLM_ROSEM(:, iter + 1) = BSREM_iter(im_vectors.OSEM_apu, options.lam_rosem, iter, SinD, randoms_correction, is_transposed, options.beta_NLM_rosem, med, epps);
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
                        else
                            SinD = 0;
                        end
                        if use_raw_data
                            xy_index = uint32(0);
                            z_index = uint32(0);
                        else
                            LL = uint16(0);
                        end
                        if ~options.precompute_lor
                            lor_a = uint16(0);
                        end
                        
                        [Summ, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                            options.normalization, SinD, pituus(end) - pituus(1), attenuation_correction, normalization_correction, randoms_correction, lor_a, xy_index, ...
                            z_index, NSinos, LL, pseudot, det_per_ring, options.verbose, (use_raw_data), uint32(1), epps, double(Sino), im_vectors.MLEM_apu, ...
                            uint32(options.projector_type), no_norm, options.precompute_lor, options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, ...
                            int32(options.accuracy_factor), options.n_rays, dc_z);
                        
                        if iter == 1 && (Niter > 1 || options.partitions > 1) && ~OS_bool && llo == 1
                            f_Summ_ml = Summ;
                        end
                        if options.mlem
                            im_vectors.MLEM(:,iter + 1) = MLEM_im(im_vectors.MLEM_apu, f_Summ_ml, epps, rhs);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['MLEM iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.MRP && options.OSL_MLEM
                            med = MRP(im_vectors.MLEM_apu, options.medx, options.medy, options.medz, Nx, Ny, Nz, epps, options.tr_offsets, options.med_no_norm);
                            im_vectors.MRP_MLEM(:,iter + 1) = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_mrp_osem, med, epps, rhs);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL MRP iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.quad && options.OSL_MLEM
                            med = Quadratic_prior(im_vectors.MLEM_apu, options.tr_offsets, options.weights, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                            im_vectors.Quad_MLEM(:, iter + 1) = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_quad_osem, med, epps, rhs);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL Quadratic iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.L && options.OSL_MLEM
                            med = L_filter(im_vectors.MLEM_apu, options.tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            im_vectors.L_MLEM(:, iter + 1) = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_L_osem, med, epps, rhs);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL L-filter iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.FMH && options.OSL_MLEM
                            med = FMH(im_vectors.MLEM_apu, options.tr_offsets, options.fmh_weights, options.weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            im_vectors.FMH_MLEM(:, iter + 1) = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_fmh_osem, med, epps, rhs);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL FMH iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.weighted_mean && options.OSL_MLEM
                            med = Weighted_mean(im_vectors.MLEM_apu, options.tr_offsets, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.w_sum, options.med_no_norm);
                            im_vectors.Weighted_MLEM(:, iter + 1) = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_weighted_osem, med, rhs);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL weighted mean iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.TV && options.OSL_MLEM
                            grad = TVpriorFinal(im_vectors.MLEM_apu, options.TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, ...
                                options.tr_offsets);
                            im_vectors.TV_MLEM(:, iter + 1) = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_TV_osem, grad, rhs);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL TV sub-iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.AD && options.OSL_MLEM
                            if iter > 1
                                med = AD(im_vectors.MLEM_apu, options.FluxType, Nx, Ny, Nz, options);
                                im_vectors.AD_MLEM(:, iter + 1) = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_ad_osem, med, epps, rhs);
                            else
                                im_vectors.AD_MLEM(:, iter + 1) = OSEM_im(im_vectors.MLEM_apu, rhs, f_Summ_ml(:,osa_iter), epps);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL AD iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.APLS && options.OSL_MLEM
                            grad = TVpriorFinal(im_vectors.MLEM_apu, [], Nx, Ny, Nz, true, options, 4);
                            im_vectors.APLS_MLEM(:, iter + 1) = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_APLS_osem, grad, rhs);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL APLS iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.TGV && options.OSL_MLEM
                            grad = TGV(im_vectors.MLEM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                            im_vectors.TGV_MLEM(:, iter + 1) = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_TGV_osem, grad, rhs);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL TGV iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.NLM && options.OSL_MLEM
                            med = NLM(im_vectors.MLEM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                options.sigma, epps, Nx, Ny, Nz, options);
                            im_vectors.NLM_MLEM(:, iter + 1) = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_NLM_osem, med, epps, rhs);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL NLM iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        no_norm = true;
                    end
                end
            end
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
        else
            if isempty(pseudot)
                pseudot = uint32(100000);
            end
            LL = uint16(0);
        end
        if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
            randoms = uint32(1);
        else
            randoms = uint32(0);
        end
        n_rays = uint16(options.n_rays);
        if options.n_rays > 1 && ~options.precompute_lor && options.projector_type == 1
            if options.n_rays > 5
                error('Only ray counts of 1 to 5 are supported')
            end
            if options.use_raw_data
                [x, y] = detector_coordinates_multiray(options);
                x = x(:,[2 1 3]);
                y = y(:,[2 1 3]);
            else
                [x, y] = sinogram_coordinates_2D_multiray(options);
                x = x(:,:,[2 1 3]);
                y = y(:,:,[2 1 3]);
            end
            x = single(x(:));
            y = single(y(:));
        end
        dc_z = single(z_det(2,1) - z_det(1,1));
        
        
        tube_width_xy = single(options.tube_width_xy);
        crystal_size_z = single(options.tube_width_z);
        if options.projector_type == 1
            kernel_file = 'siddon_kernel_matrixfree_GPU_noprecomp.cl';
            kernel_path = which(kernel_file);
            kernel_path = strrep(kernel_path, '\', '/');
            kernel_path = strrep(kernel_path, '.cl', '');
            %         kernel_file = 'siddon_kernel_matrixfree_GPU - Copy (2).cl';
            filename = 'OMEGA_matrix_free_OpenCL_binary_siddon_device';
            header_directory = strrep(kernel_path,'siddon_kernel_matrixfree_GPU_noprecomp','');
        elseif options.projector_type == 2
            kernel_file = 'orthogonal_kernel_matrixfree_noprecomp.cl';
            kernel_path = which(kernel_file);
            kernel_path = strrep(kernel_path, '\', '/');
            kernel_path = strrep(kernel_path, '.cl', '');
            filename = 'OMEGA_matrix_free_OpenCL_binary_orth_device';
            header_directory = strrep(kernel_path,'orthogonal_kernel_matrixfree_noprecomp','');
        end
        filename = [header_directory, filename];
        header_directory = strcat('-I "', header_directory);
        header_directory = strcat(header_directory,'"');
        joku = algorithms_char();
        %         n_rekos = uint32(sum(rekot(~contains(joku,'MLEM'))));
        n_rekos = uint32(sum(rekot(cellfun('isempty',strfind(joku,'MLEM')))));
        n_rekos_mlem = uint32(sum(rekot(~cellfun('isempty',strfind(joku,'MLEM')))));
        reko_type = zeros(length(rekot),1,'uint8');
        reko_type(~cellfun('isempty',strfind(joku,'MBSREM'))) = 1;
        reko_type(~cellfun('isempty',strfind(joku,'MRAMLA'))) = 1;
        reko_type(~cellfun('isempty',strfind(joku,'COSEM'))) = 2;
        reko_type(~cellfun('isempty',strfind(joku,'ACOSEM'))) = 3;
        reko_type(~cellfun('isempty',strfind(joku,'OSL-COSEM')) & options.COSEM_MAP == 1) = 2;
        reko_type(~cellfun('isempty',strfind(joku,'OSL-COSEM')) & options.COSEM_MAP == 2) = 3;
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
        %         reko_type(contains(joku,'MBSREM')) = 1;
        %         reko_type(contains(joku,'MRAMLA')) = 1;
        %         reko_type(contains(joku,'COSEM')) = 2;
        %         reko_type(contains(joku,'ACOSEM')) = 3;
        %         reko_type(contains(joku,'OSL-COSEM') & options.COSEM_MAP == 1) = 2;
        %         reko_type(contains(joku,'OSL-COSEM') & options.COSEM_MAP == 2) = 3;
        tic
        [pz] = OpenCL_matrixfree( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end) , NSinos, single(NSlices), size_x, zmax, NSinos, ...
            options.verbose, LL, pseudot, det_per_ring, device, uint8(use_raw_data), filename, uint32(0), options.force_build, header_directory, options.vaimennus, options.normalization, ...
            pituus, uint32(attenuation_correction), uint32(normalization_correction), uint32(Niter), uint32(subsets), uint8(rekot), single(epps), lor_a, xy_index, ...
            z_index, any(n_rekos), tube_width_xy, crystal_size_z, x_center, y_center, z_center, options.SinDelayed, randoms, uint32(options.projector_type), options.precompute_lor, ...
            int32(options.accuracy_factor), n_rays, dc_z, options, options.SinM, uint32(partitions), logical(options.use_64bit_atomics), n_rekos, n_rekos_mlem, reko_type);
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
        else
            if isempty(pseudot)
                pseudot = uint32(100000);
            end
            LL = uint16(0);
        end
        if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
            randoms = uint32(1);
        else
            randoms = uint32(0);
        end
        n_rays = uint16(options.n_rays);
        if options.n_rays > 1 && ~options.precompute_lor && options.projector_type == 1
            if options.n_rays > 5
                error('Only ray counts of 1 to 5 are supported')
            end
            if options.use_raw_data
                [x, y] = detector_coordinates_multiray(options);
                x = x(:,[2 1 3]);
                y = y(:,[2 1 3]);
            else
                [x, y] = sinogram_coordinates_2D_multiray(options);
                x = x(:,:,[2 1 3]);
                y = y(:,:,[2 1 3]);
            end
            x = single(x(:));
            y = single(y(:));
        end
        dc_z = single(z_det(2,1) - z_det(1,1));
        
        tube_width_xy = single(options.tube_width_xy);
        crystal_size_z = single(options.tube_width_z);
        if options.projector_type == 1 && options.precompute_lor
            kernel_file = 'multidevice_siddon.cl';
            kernel_path = which(kernel_file);
            kernel_path = strrep(kernel_path, '\', '/');
            kernel_path = strrep(kernel_path, '.cl', '');
            filename = 'OMEGA_matrix_free_OpenCL_binary_device';
            header_directory = strrep(kernel_path,'multidevice_siddon','');
        elseif options.projector_type == 2 && options.precompute_lor
            filename = 'OMEGA_matrix_free_orthogonal_OpenCL_binary_device';
            kernel_file = 'multidevice_orth.cl';
            kernel_path = which(kernel_file);
            kernel_path = strrep(kernel_path, '\', '/');
            kernel_path = strrep(kernel_path, '.cl', '');
            header_directory = strrep(kernel_path,'multidevice_orth','');
        elseif options.projector_type == 1 && ~options.precompute_lor
            kernel_file = 'multidevice_siddon_no_precomp.cl';
            kernel_path = which(kernel_file);
            kernel_path = strrep(kernel_path, '\', '/');
            kernel_path = strrep(kernel_path, '.cl', '');
            filename = 'OMEGA_matrix_free_OpenCL_binary_device';
            header_directory = strrep(kernel_path,'multidevice_siddon_no_precomp','');
        elseif options.projector_type == 2 && ~options.precompute_lor
            kernel_file = 'multidevice_orth_no_precomp.cl';
            kernel_path = which(kernel_file);
            kernel_path = strrep(kernel_path, '\', '/');
            kernel_path = strrep(kernel_path, '.cl', '');
            filename = 'OMEGA_matrix_free_OpenCL_binary_device';
            header_directory = strrep(kernel_path,'multidevice_orth_no_precomp','');
        else
            error('Invalid projector for OpenCL')
        end
        filename = [header_directory, filename];
        header_directory = strcat('-I "', header_directory);
        header_directory = strcat(header_directory,'"');
        tic
        [tz] = OpenCL_matrixfree_multi_gpu( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end), single(NSlices), size_x, zmax, options.verbose, ...
            LL, pseudot, det_per_ring, uint32(options.use_device), filename, uint8(use_raw_data), single(options.cpu_to_gpu_factor), uint32(1), header_directory, ...
            options.vaimennus, options.normalization, pituus, uint32(attenuation_correction), uint32(normalization_correction), lor_a, xy_index, z_index, tube_width_xy, ...
            crystal_size_z, x_center, y_center, z_center, options.SinDelayed, randoms, uint32(options.projector_type), options.precompute_lor, ...
            int32(options.accuracy_factor), n_rays, dc_z, NSinos, uint16(TotSinos), uint32(Niter), uint32(subsets), uint8(rekot), single(epps), options.SinM, uint32(partitions), ...
            options.osem, options.force_build, options, logical(options.use_64bit_atomics));
        toc
        
        for ll = 1 : options.partitions
            apu = tz{1,ll};
            apu(:,1) = options.x0;
            if options.mlem
                pz{1,ll} = reshape(apu, Nx, Ny, Nz, Niter + 1);
            elseif options.osem
                pz{2,ll} = reshape(apu, Nx, Ny, Nz, Niter + 1);
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
