function [pz,varargout] = reconstructions_main(options,varargin)
%% Main reconstruction file
% This function is used to compute various reconstructions with the
% selected method. Can be used with any sinogram or raw data.
%
% OUTPUT:
%   pz = A cell matrix containing output from each of the selected
%   algorithms and/or priors. E.g. if OSEM and ROSEM are selected pz{2}
%   contains the OSEM estimates and pz{5} the ROSEM estimates.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021 Ville-Veikko Wettenhovi, Samuli Summala
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

% Check if certain variables are present
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
if ~isfield(options,'compute_sensitivity_image')
    options.compute_sensitivity_image = false;
end
if ~isfield(options,'listmode')
    options.listmode = false;
end
if ~isfield(options,'sampling_raw')
    options.sampling_raw = 1;
end
if ~isfield(options,'CT')
    options.CT = false;
end

if nargin > 1
    tyyppi = varargin{1};
else
    tyyppi = 0;
end

options.listmode = false;
tStart = 0;

% Is TOF enabled?
TOF = options.TOF_bins > 1 && options.projector_type == 1;

folder = fileparts(which('reconstructions_main.m'));
folder = strrep(folder, 'source','mat-files/');
folder = strrep(folder, '\','/');

% Special requirements for span 1 case
if options.span == 1
    options.TotSinos = options.rings^2;
    options.NSinos = options.TotSinos;
end

options = convertOptions(options);

if tyyppi == 0
    disp('Preparing for reconstruction')
end

Nx = options.Nx;
Ny = options.Ny;
Nz = options.Nz;
N = Nx * Ny * Nz;
options.N = N;


var = recNames(2);
ll = 0;
kk = 1;
while ll == 0 && kk <= numel(var)
    ll = ll + options.(var{kk});
    kk = kk +1;
end
options.MAP = ll > 0;

% options.MAP = (options.OSL_MLEM || options.OSL_OSEM || options.BSREM || options.MBSREM || options.ROSEM_MAP || options.OSL_RBI...
%     || any(options.OSL_COSEM) || options.PKMA);

var = recNames(3);
ll = 0;
kk = 1;
while ll == 0 && kk <= numel(var)
    ll = ll + options.(var{kk});
    kk = kk +1;
end
MLEM_bool = ll > 0;

var = recNames(4);
ll = 0;
kk = 1;
while ll == 0 && kk <= numel(var)
    ll = ll + options.(var{kk});
    kk = kk +1;
end
OS_bool = ll > 0;

if tyyppi < 2
    % Load the measurement data if it does not exist in options.SinM
    % Raw data
    if options.use_raw_data
        RandProp.smoothing = false;
        RandProp.variance_reduction = false;
        ScatterProp.smoothing = false;
        ScatterProp.variance_reduction = false;
        ScatterProp.normalization = false;
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
                if ~isfield(options,'SinDelayed')
                    options = loadDelayedData(options);
                end
            end
        end
        if options.scatter_correction && ~options.corrections_during_reconstruction
            if ~isfield(options,'ScatterC')
                options = loadScatterData(options);
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
            loadRaw = false;
            if isfield(options,'SinM') == 0
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
            if ~options.corrections_during_reconstruction && ~isempty(appliedCorrections)
                normalization_correction = options.normalization_correction;
                if appliedCorrections.normalization && ~options.normalization_correction
                    warning('Normalization correction not selected, but data precorrected with normalization! Precorrecting without normalization.')
                    [options.SinM, appliedCorrections] = loadStructFromFile(sinoFile, 'raw_SinM','appliedCorrections');
                    loadRaw = true;
                end
                if appliedCorrections.normalization
                    normalization_correction = false;
                end
                randoms_correction = options.randoms_correction;
                if ~isempty(strfind(appliedCorrections.randoms,'variance reduction')) && ~options.variance_reduction && options.randoms_correction
                    warning('Randoms variance correction not selected, but data precorrected with randoms with applied variance reduction! Precorrecting without variance reduction.')
                    [options.SinM, appliedCorrections] = loadStructFromFile(sinoFile, 'raw_SinM','appliedCorrections');
                    loadRaw = true;
                    if appliedCorrections.normalization && ~normalization_correction
                        normalization_correction = true;
                    end
                elseif ~isempty(strfind(appliedCorrections.randoms,'smoothing')) && ~options.randoms_smoothing && options.randoms_correction
                    warning('Randoms smoothing not selected, but data precorrected with randoms with applied smoothing! Precorrecting without smoothing.')
                    [options.SinM, appliedCorrections] = loadStructFromFile(sinoFile, 'raw_SinM','appliedCorrections');
                    loadRaw = true;
                    if appliedCorrections.normalization && ~normalization_correction
                        normalization_correction = true;
                    end
                elseif isempty(strfind(appliedCorrections.randoms,'randoms correction')) && ~loadRaw && randoms_correction
                    [options.SinM, appliedCorrections] = loadStructFromFile(sinoFile, 'raw_SinM','appliedCorrections');
                    loadRaw = true;
                    if appliedCorrections.normalization && ~normalization_correction
                        normalization_correction = true;
                    end
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
                if (~appliedCorrections.gapFilling || loadRaw) && options.fill_sinogram_gaps
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
        if options.partitions == 1 && options.randoms_correction && ~options.reconstruct_scatter && ~options.reconstruct_trues
            if options.use_machine == 0 || options.use_machine == 1 || options.use_machine == 3
                try
                    [options.SinDelayed,RandProp] = loadStructFromFile(sinoFile,'SinDelayed','RandProp');
                catch
                    options = loadDelayedData(options);
                end
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
        elseif options.randoms_correction && ~options.reconstruct_scatter && ~options.reconstruct_trues
            if options.use_machine == 0 || options.use_machine == 1 || options.use_machine == 3
                try
                    [options.SinDelayed,RandProp] = loadStructFromFile(sinoFile,'SinDelayed','RandProp');
                catch
                    options = loadDelayedData(options);
                end
            else
                options = loadDelayedData(options);
            end
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
    
    if (options.quad || options.FMH || options.L || options.weighted_mean || options.MRP || options.TV || options.Huber) && options.MAP
        Ndx = options.Ndx;
        Ndy = options.Ndy;
        Ndz = options.Ndz;
    else
        Ndx = 0;
        Ndy = 0;
        Ndz = 0;
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
    
    
    if (options.implementation == 1 || options.implementation == 4) || tyyppi == 1
        im_vectors = form_image_vectors(options, N);
    end
end

% Vector to identify the reconstruction algorithms
rekot = reko_maker(options);

Niter = options.Niter;
subsets = options.subsets;
epps = options.epps;
if tyyppi == 0
    precompute_obs_matrix = options.precompute_obs_matrix;
    attenuation_correction = options.attenuation_correction;
end
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
if ~isfield(options, 'global_correction_factor') || isempty(options.global_correction_factor)
    if options.implementation == 2 || options.implementation == 3
        options.global_correction_factor = single(1);
    else
        options.global_correction_factor = 1;
    end
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


% MLEM_bool = options.OSL_MLEM || options.MLEM;
% OS_bool = options.OSEM || options.rosem || options.ramla || options.OSL_OSEM || options.BSREM || options.ROSEM_MAP || options.RBI || options.drama ...
%     || options.COSEM || options.ECOSEM || options.ACOSEM || options.OSL_RBI || any(options.OSL_COSEM) || options.PKMA;

pz = cell(length(rekot),partitions);


temp = pseudot;
if ~isempty(temp) && temp > 0
    for kk = uint32(1) : temp
        pseudot(kk) = uint32(options.cryst_per_block + 1) * kk;
    end
elseif temp == 0
    pseudot = [];
end



if options.implementation == 3 && tyyppi == 0
    if options.OSEM && options.MLEM
        if subsets == 1
            disp(['Both OSEM and MLEM set for method ' num2str(options.implementation) ', using MLEM'])
            options.OSEM = false;
        else
            disp(['Both OSEM and MLEM set for method ' num2str(options.implementation) ', using OSEM'])
            options.MLEM = false;
        end
    end
end

% if options.CT && (options.xSize == 1 || size(options.SinM,3) == 1 || numel(options.SinM) < 6340608) && isfield(options,'x') && numel(options.x)/2 ~= numel(options.SinM)
%     options.precompute_lor = false;
if options.CT && options.implementation == 1
    options.precompute_lor = true;
end

if options.precompute_lor
    is_transposed = true;
else
    is_transposed = false;
end

if (options.implementation == 3 || options.implementation == 1 ) && ~any(rekot(cellfun('isempty',strfind(algorithms_char(),'MLEM'))))
    subsets = 1;
end

% TOF parameters
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

% Whether list-mode or sinogram/raw data is used
if (isfield(options,'x') && isfield(options,'y') && (isfield(options,'z') || isfield(options,'z_det'))) && numel(options.x) / 2 == numel(options.SinM)
    index = 0;
    det_per_ring = numel(options.SinM);
    pituus = floor(det_per_ring / subsets);
    pituus = int64([repmat(pituus,subsets - 1,1); det_per_ring - pituus*(subsets - 1)]);
    list_mode_format = true;
    options.listmode = true;
    %     if options.implementation == 2 && options.use_CUDA
    %         error('CUDA support not enabled for list-mode data')
    %     end
    if options.implementation == 4 || options.implementation == 1
        use_raw_data = true;
        options.use_raw_data = use_raw_data;
    end
    if abs(min(options.x(:))) < abs(max(options.x(:))) / 2 && options.diameter == 0
        diameter = (min(options.x(:))) + (max(options.x(:)));
    elseif abs(min(options.x(:))) < abs(max(options.x(:))) / 2 && options.diameter > 0
        diameter = options.diameter;
    else
        diameter = 0;
    end
    if isfield(options,'z')
        if abs(min(options.z(:))) < abs(max(options.z(:))) / 2
            if min(options.z(:)) < 0
                Z = options.axial_fov - min(options.z(:)) * 2;
            elseif max(options.z(:)) > options.axial_fov
                Z = options.axial_fov + (max(options.z(:)) - options.axial_fov) * 2;
            else
                Z = options.axial_fov;
            end
        else
            Z = 0;
        end
    else
        if abs(min(options.z_det(:))) < abs(max(options.z_det(:))) / 2
            if min(options.z_det(:)) < 0
                Z = options.axial_fov - min(options.z_det(:)) * 2;
            elseif max(options.z_det(:)) > options.axial_fov
                Z = options.axial_fov + (max(options.z_det(:)) - options.axial_fov) * 2;
            else
                Z = options.axial_fov;
            end
        else
            Z = 0;
        end
    end
    options.Z = Z;
else
    % Compute the indices for the subsets used.
    % For Sinogram data, six different methods to select the subsets are
    % available. For data, three methods are available.
    options.listmode = false;
    [index, pituus, subsets, lor_a, lor_orth] = index_maker(Nx, Ny, Nz, subsets, use_raw_data, machine_name, options, Nang, Ndist, TotSinos, NSinos);
    Z = axial_fov;
    options.Z = Z;
end

%%

% Diameter of the PET-scanner (bore) (mm)
R = double(diameter);
% Transaxial FOV (x-direction, horizontal) (mm)
FOVax = double(FOVax);
% Transaxial FOV (y-direction, vertical) (mm)
FOVay = double(FOVay);
% Number of rings
blocks = uint32(rings - 1);
% First ring
block1=uint32(0);

NSinos = uint32(NSinos);
NSlices = uint32(Nz);

if ~list_mode_format
    % Load correction data
    [normalization_correction, randoms_correction, options] = set_up_corrections(options, blocks, RandProp, ScatterProp);
else
    % list-mode does not support corrections
    options.normalization_correction = false;
    options.randoms_correction = false;
    options.scatter_correction = false;
    [normalization_correction, randoms_correction, options] = set_up_corrections(options, blocks, RandProp, ScatterProp);
end

% Coordinates of the detectors
[x, y, z_det, options] = get_coordinates(options, blocks, pseudot);

if options.use_raw_data
    if list_mode_format
        size_x = uint32(numel(x) / 2);
    else
        size_x = uint32(numel(x));
    end
else
    if list_mode_format
        size_x = uint32(numel(x) / 2);
    else
        size_x = uint32(options.Nang*options.Ndist);
    end
    if options.sampling > 1 && ~options.precompute_lor
        size_x = size_x * options.sampling;
    end
end
if options.CT
    size_x = uint32(options.ySize);
    options.size_y = uint32(options.xSize);
    if options.listmode
        %         size_x = size_x * uint32(options.xSize * options.nProjections);
        size_x = uint32(numel(x) / 2);
    end
    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
        options.angles = single(options.angles);
        options.dPitch = single(options.dPitch);
        options.nProjections = int64(options.nProjections);
        options.xSize = uint32(options.xSize);
    end
else
    options.angles = 0;
    options.dPitch = 0;
    options.size_y = uint32(0);
    options.nProjections = int64(0);
    options.xSize = 0;
end

if subsets > 1
    pituus = [int64(0);int64(cumsum(pituus))];
    if iscell(index)
        index = cell2mat(index);
    end
end

% Compute the necessary indices required for subsets (e.g. the index of
% the detector coordinates for the current LOR)
if ~list_mode_format
    [options, lor_a, xy_index, z_index, LL, summa, pituus, options.SinM, lor_orth] = form_subset_indices(options, pituus, subsets, index, size_x, y, z_det, rings, false, TOF, ...
        options.SinM, lor_a, lor_orth);
else
    LL = uint16(0);
    xy_index = uint32(0);
    z_index = uint16(0);
    lor_orth = uint16(0);
    summa = zeros(subsets, 1, 'uint64');
end
if options.implementation > 1
    index = uint32(0);
end
if ~options.precompute_lor
    lor_a = uint16(0);
    lor_orth = uint16(0);
elseif options.precompute_lor && exist('lor_a','var') ~= 1
    [lor_a, ~, ~, lor_orth] = lor_pixel_count_prepass(options, false);
end

if use_raw_data
    if isempty(pseudot)
        pseudot = uint32(1e5);
    else
        pseudot = pseudot - 1;
    end
end

if options.CT
    R = 0;
    Z = options.dPitch * double(options.xSize) + z_det(options.nProjections);
end
[xx,yy,zz,dx,dy,dz,bx,by,bz] = computePixelSize(R, FOVax, FOVay, Z, axial_fov, Nx, Ny, Nz, options.implementation);

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
if options.CT
    zmax = z_det(options.nProjections) + options.dPitch * (double(options.xSize) - 1);
    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
        zmax = single(zmax);
    end
end

[x_center,y_center,z_center,dec] = computePixelCenters(xx,yy,zz,dx,dy,dz,TOF,options);

[V,Vmax,bmin,bmax] = computeVoxelVolumes(dx,dy,dz,options);

if options.implementation == 1
    iij = double(0:Nx);
    jji = double(0:Ny);
    kkj = double(0:Nz);
end

if exist('feature','builtin') == 5
    nCores = uint32(feature('numcores'));
else
    nCores = uint32(1);
end

if options.use_raw_data && options.sampling_raw > 1
    det_per_ring = det_per_ring * options.sampling_raw;
end

% Compute PSF kernel
[gaussK, options] = PSFKernel(options);

%% This computes a whole observation matrix and uses it to compute the MLEM (no on-the-fly calculations)
% NOTE: Only attenuation correction is supported
% This section is largely untested
if tyyppi == 0 && precompute_obs_matrix && options.implementation == 1
    
    for llo = 1 : partitions
        
        
        if options.MLEM
            if TOF
                error('TOF is not supported when using precomputed system matrix')
            end
            if iscell(options.SinM)
                Sino = options.SinM{llo};
            else
                Sino = options.SinM;
                clear options.SinM
            end
            
            Sino = double(Sino(:));
            
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
                options.gaussK = gaussK;
                options.is_transposed = is_transposed;
                
                A = observation_matrix_formation(options);
                if is_transposed
                    D = full(sum(A,2));
                else
                    D = full(sum(A,1))';
                end
                D(D < epps) = epps;
            end
            
            if options.save_iter
                MLEM = ones(N,Niter + 1);
            else
                MLEM = ones(N,1);
            end
            MLEM(:,1) = options.x0(:);
            Sino = double(Sino);
            
            for iter=1:Niter
                if options.save_iter
                    iter_n = iter + 1;
                else
                    iter_n = 1;
                end
                MLEM(:, iter_n) = MLEM_im(MLEM(:,iter_n), D, epps, A, Sino, is_transposed, options, options.Nx, options.Ny, options.Nz, gaussK);
                disp(['MLEM iteration ' num2str(iter) ' finished'])
                if ~any(rekot(cellfun('isempty',strfind(algorithms_char(),'MLEM'))))
                    warning('Only MLEM is supported with precomputed observation matrix')
                end
            end
            if options.save_iter
                MLEM = reshape(MLEM,Nx,Ny,Nz,Niter+1);
            else
                MLEM = reshape(MLEM,Nx,Ny,Nz);
            end
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
    
    % Perform various prepass steps, if necessary
    % These include computing weights, matrices required by some algorithms
    % (COSEM, etc.) and loading anatomic reference images
    [options, D, C_co, C_aco, C_osl, Amin, E] = prepass_phase(options, pituus, index, options.SinM, pseudot, x, y, xx, yy, z_det, dz, dx, dy, bz, bx, by, NSlices, ...
        zmax, size_x, block1, blocks, normalization_correction, randoms_correction, xy_index, z_index, lor_a, lor_orth, summa, LL, is_transposed, x_center, y_center, ...
        z_center, ind_size, gaussK, bmin, bmax, Vmax, V);
    
    if options.use_psf && ((options.MRAMLA || options.MBSREM || options.OSL_RBI || options.RBI || options.PKMA) && options.MBSREM_prepass || options.ECOSEM || options.COSEM ...
            || options.ACOSEM || any(options.OSL_COSEM)) && options.implementation == 1
        D = computeConvolution(D, options, Nx, Ny, Nz, gaussK);
    end
    % Sensitivity image for MRAMLA/MBSREM
    if (options.MBSREM || options.MRAMLA || options.PKMA) && options.implementation == 1
        options.pj3 = D/options.subsets;
    end
    if tyyppi > 0
        % Upper bound for MRAMLA
        if options.MBSREM || options.MRAMLA
            if options.U == 0 || isempty(options.U)
                if options.CT
                    if options.implementation == 1
                        if iscell(options.SinM)
                            options.U = max(-log(1./options.SinM{1})./Amin);
                        else
                            options.U = max(-log(1./options.SinM)./Amin);
                        end
                    else
                        if iscell(options.SinM)
                            options.U = max(-log(options.SinM{1})./Amin);
                        else
                            options.U = max(-log(options.SinM)./Amin);
                        end
                    end
                else
                    if iscell(options.SinM)
                        options.U = max(double(options.SinM{1}) ./ Amin);
                    else
                        options.U = max(double(options.SinM) ./ Amin);
                    end
                end
            end
        end
        % Compute the epsilon value for MRAMLA
        if options.MBSREM || options.MRAMLA
            if iscell(options.SinDelayed)
                if iscell(options.SinM)
                    options.epsilon_mramla = MBSREM_epsilon(options.SinM{1}, epps, randoms_correction, options.SinDelayed{1}, E);
                else
                    options.epsilon_mramla = MBSREM_epsilon(options.SinM, epps, randoms_correction, options.SinDelayed{1}, E);
                end
            else
                if iscell(options.SinM)
                    options.epsilon_mramla = MBSREM_epsilon(options.SinM{1}, epps, randoms_correction, options.SinDelayed, E);
                else
                    options.epsilon_mramla = MBSREM_epsilon(options.SinM, epps, randoms_correction, options.SinDelayed, E);
                end
            end
        end
        options.det_per_ring = det_per_ring;
        options.ind_size = ind_size;
        options.zmax = zmax;
        options.x_center = x_center;
        options.y_center = y_center;
        options.z_center = z_center;
        options.dec = dec;
        options.V = V;
        options.Vmax = Vmax;
        options.bmin = bmin;
        options.bmax = bmax;
        options.dx = dx;
        options.dy = dy;
        options.dz = dz;
        options.xx = xx;
        options.yy = yy;
        options.zz = zz;
        options.D = D;
        clear D
        options.C_co = C_co;
        clear C_co
        options.C_aco = C_aco;
        clear C_aco
        options.C_osl = C_osl;
        clear C_osl
        options.Amin = Amin;
        clear Amin
        options.E = E;
        clear E
        options.pituus = pituus;
        options.index = index;
        clear index
        options.pseudot = pseudot;
        options.x = x;
        options.y = y;
        options.z_det = z_det;
        options.bz = bz;
        options.bx = bx;
        options.by = by;
        options.NSlices = NSlices;
        options.size_x = size_x;
        options.block1 = block1;
        options.blocks = blocks;
        options.normalization_correction = normalization_correction;
        options.randoms_correction = randoms_correction;
        options.xy_index = xy_index;
        clear xy_index
        options.z_index = z_index;
         clear z_index
        options.lor_a = lor_a;
         clear lor_a
        options.lor_orth = lor_orth;
         clear lor_orth;
        options.summa = summa;
        options.LL = LL;
        clear LL
        options.is_transposed = is_transposed;
        options.TOFCenter = TOFCenter;
        options.sigma_x = sigma_x;
        options.im_vectors = im_vectors;
        clear im_vectors
        options.OS_bool = OS_bool;
        options.MLEM_bool = MLEM_bool;
        varargout{1} = options;
        return
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
            if options.scatter_correction && ~options.subtract_scatter
                if iscell(options.ScatterC)
                    if length(options.ScatterC) > 1
                        ScatterC = double(options.ScatterC{llo});
                    else
                        ScatterC = double(options.ScatterC{1});
                    end
                else
                    ScatterC = double(options.ScatterC);
                end
            else
                ScatterC = 0;
            end
            
            Sino = Sino(:);
            
            if issparse(Sino)
                Sino = (full(Sino));
            end
            
            % Implementation 1
            if options.implementation == 1
                % Upper bound for MRAMLA
                if options.MBSREM || options.MRAMLA
                    if options.U == 0 || isempty(options.U)
                        if options.CT
                            if options.implementation == 1
                                options.U = max(-log(1./Sino)./Amin);
                            else
                                options.U = max(-log(Sino)./Amin);
                            end
                        else
                            options.U = max(double(Sino)./Amin);
                        end
                    end
                end
                % Compute the epsilon value for MRAMLA
                if options.MBSREM || options.MRAMLA
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
                            norm_input = single(options.normalization(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                        else
                            norm_input = single(0);
                        end
                        if options.scatter_correction && ~options.subtract_scatter
                            scatter_input = double(ScatterC(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                        else
                            scatter_input = 0;
                        end
                        if options.listmode
                            if iter == 1 && osa_iter == 1
                                lor_a = lor_pixel_count_prepass(options, false);
                                for lla = 1 : subsets
                                    summa(lla) = uint64(sum(lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1))));
                                end
                            end
                        end
                        uu = double(Sino(pituus(osa_iter) + 1 : pituus(osa_iter + 1)));
                        koko = length(uu);
                        [A,ll, Summ] = computeImplementation1(options,use_raw_data,randoms_correction, pituus,osa_iter, normalization_correction,...
                            Nx, Ny, Nz, dx, dy, dz, bx, by, bz, x, y, z_det, xx, yy, size_x, NSinos, NSlices, zmax, attenuation_correction, pseudot, det_per_ring, ...
                            TOF, sigma_x, TOFCenter, dec, nCores, ind_size, block1, blocks, index, iij, jji, kkj, LL, N, summa, lor_a, xy_index, z_index, ...
                            x_center, y_center, z_center, bmin, bmax, Vmax, V, lor_orth, gaussK,is_transposed, scatter_input, norm_input, SinD, koko);
                        if options.attenuation_phase
                            uu = uu ./ ll;
                        end
                        [im_vectors,C_co,C_aco,C_osl] = computeEstimatesImp1(im_vectors, options, A, uu, Summ, SinD, is_transposed, gaussK, iter, osa_iter, C_co, C_aco,C_osl,...
                            randoms_correction, N, Ndx, Ndy, Ndz, D);
                        clear A
                    end
                    im_vectors = init_next_iter(im_vectors, options, iter);
                    if options.use_psf && options.deblurring
                        im_vectors = computeDeblur(im_vectors, options, iter, gaussK, Nx, Ny, Nz);
                    end
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
                
                if ~isfield(options, 'vaimennus')
                    options.vaimennus = 0;
                end
                if ~isfield(options, 'normalization')
                    options.normalization = 0;
                end
                if ~isfield(options, 'scatter')
                    options.scatter = false;
                end
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
                %%% SENSITIVITY IMAGE FOR LISTMODE DATA %%%
                if options.listmode && options.compute_sensitivity_image
                    %                     vec = repmat((1:options.det_per_ring)', options.rings,1);
                    %                     vec2 = repelem((1:options.rings)',options.det_per_ring);
                    %                     pituusD = int64(sum(vec .* vec2 * options.det_per_ring));
                    options.listmode = uint8(2);
                    [xd, yd] = detector_coordinates(options);
                    
                    z_length = double(rings + 1 + sum(options.pseudot)) * options.cr_pz;
                    zd = linspace(0, z_length, rings + 2 + sum(options.pseudot))';
                    if sum(options.pseudot) > 0
                        zd(options.pseudot) = [];
                    end
                    if min(zd(:)) == 0
                        zd = zd + (options.axial_fov - (options.rings + sum(options.pseudot)) * options.cr_pz)/2 + options.cr_pz/2;
                    end
                    zd = zd(1:end-2);
                    LL = form_detector_pairs_raw(options.rings, options.det_per_ring)';
                    LL = LL(:);
                    pituusD = int64(numel(LL)/2);
                    if options.projector_type == 1
                        if exist('OCTAVE_VERSION','builtin') == 0
                            [f_Summ, ~] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, zd, xd, yd, dy, yy, xx , NSinos, NSlices, size_x, max(zd(:)), options.vaimennus, ...
                                0, 0, pituusD, attenuation_correction, false, false,...
                                false, 0, options.global_correction_factor, uint16(0), uint32(0), uint16(0), NSinos, LL, pseudot, uint32(options.det_per_ring), ...
                                TOF, int64(0), sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, nCores, ...
                                (use_raw_data), uint32(1), options.listmode, epps, 0, im_vectors.OSEM_apu, uint32(options.projector_type), false, false, false, ...
                                options.n_rays_transaxial, options.n_rays_axial, dc_z);
                        elseif exist('OCTAVE_VERSION','builtin') == 5
                            [f_Summ, ~] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, zd, xd, yd, dy, yy, xx , NSinos, NSlices, size_x, max(zd(:)), options.vaimennus, ...
                                0, 0, pituusD, attenuation_correction, false, false,...
                                false, 0, options.global_correction_factor, uint16(0), uint32(0), uint16(0), NSinos, LL, pseudot, uint32(options.det_per_ring), ...
                                TOF, int64(0), sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, nCores, ...
                                (use_raw_data), uint32(1), options.listmode, epps, 0, im_vectors.OSEM_apu, uint32(options.projector_type), false, false, false, ...
                                options.n_rays_transaxial, options.n_rays_axial, dc_z);
                        end
                    else
                        error('Unsupported projector')
                    end
                    if options.use_psf
                        f_Summ = computeConvolution(f_Summ, options, Nx, Ny, Nz, gaussK);
                    end
                    f_Summ(f_Summ < epps) = epps;
                    f_Summ = repmat(f_Summ, 1, subsets);
                    no_norm = true;
                    options.listmode = uint8(1);
                    LL = uint16(0);
                end
                %%% PREPASS PHASE %%%
                if options.COSEM || options.ECOSEM || options.ACOSEM || options.OSL_RBI || options.RBI || options.PKMA || any(options.OSL_COSEM)
                    if llo == 1 && ~options.compute_sensitivity_image
                        f_Summ = zeros(Nx*Ny*Nz,subsets);
                    end
                    D = zeros(Nx*Ny*Nz, 1);
                    if options.COSEM || options.ECOSEM || options.OSL_COSEM == 2
                        options.h = 1;
                    end
                    if options.COSEM || options.ECOSEM
                        C_co = zeros(Nx*Ny*Nz,subsets);
                    elseif options.ACOSEM
                        C_aco = zeros(Nx*Ny*Nz,subsets);
                    elseif any(options.OSL_COSEM)
                        C_osl = zeros(Nx*Ny*Nz,subsets);
                    end
                    if options.ECOSEM
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
                                SinD = single(full(options.SinDelayed{llo}(pituus(osa_iter)+1:pituus(osa_iter + 1))));
                            else
                                SinD = single(full(options.SinDelayed(pituus(osa_iter)+1:pituus(osa_iter + 1))));
                            end
                            SinD = SinD(:);
                            if TOF
                                SinD = SinD / options.TOF_bins;
                            end
                        else
                            SinD = 0;
                        end
                        if normalization_correction
                            norm_input = single(options.normalization(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                        else
                            norm_input = 0;
                        end
                        if options.scatter_correction && ~options.subtract_scatter
                            scatter_input = double(ScatterC(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                        else
                            scatter_input = 0;
                        end
                        uu = single(full(Sino(pituusS(osa_iter)+1:pituusS(osa_iter + 1))));
                        if use_raw_data
                            if ~list_mode_format
                                L_input = LL(pituus(osa_iter) * 2 + 1 : pituus(osa_iter + 1) * 2);
                            else
                                L_input = LL;
                                apux = x;
                                apuy = y;
                                apuz = z_det;
                                x = reshape(x, numel(x)/2,2);
                                y = reshape(y, numel(y)/2,2);
                                z_det = reshape(z_det, numel(z_det)/2,2);
                                x = x(pituus(osa_iter) + 1 : pituus(osa_iter + 1),:);
                                y = y(pituus(osa_iter) + 1 : pituus(osa_iter + 1),:);
                                z_det = z_det(pituus(osa_iter) + 1 : pituus(osa_iter + 1),:);
                                x = x(:);
                                y = y(:);
                                z_det = z_det(:);
                                det_per_ring = uint32(numel(x)/2);
                            end
                            xy_index_input = uint32(0);
                            z_index_input = uint32(0);
                            TOFSize = int64(size(L_input,1));
                        else
                            L_input = uint16(0);
                            if ~list_mode_format
                                xy_index_input = xy_index(pituus(osa_iter)+1:pituus(osa_iter + 1));
                                z_index_input = z_index(pituus(osa_iter)+1:pituus(osa_iter + 1));
                            else
                                xy_index_input = xy_index;
                                z_index_input = z_index;
                            end
                            TOFSize = int64(numel(xy_index_input));
                        end
                        if options.precompute_lor
                            lor_a_input = lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1));
                        else
                            lor_a_input = uint16(0);
                        end
                        [Summ,rhs] = computeImplementation4(options,use_raw_data,randoms_correction, pituus,osa_iter, normalization_correction,...
                            Nx, Ny, Nz, dx, dy, dz, bx, by, bz, x, y, z_det, xx, yy, size_x, NSinos, NSlices, zmax, attenuation_correction, pseudot, det_per_ring, ...
                            TOF, TOFSize, sigma_x, TOFCenter, dec, nCores, L_input, lor_a_input, xy_index_input, z_index_input, epps, uu, OSEM_apu, no_norm, ...
                            x_center, y_center, z_center, bmin, bmax, Vmax, V, scatter_input, norm_input, SinD, dc_z);
                        
                        if list_mode_format
                            x = apux;
                            y = apuy;
                            z_det = apuz;
                        end
                        if options.COSEM || options.ECOSEM
                            C_co(:, osa_iter) = im_vectors.OSEM_apu .* rhs;
                        elseif options.ACOSEM
                            C_aco(:, osa_iter) = im_vectors.OSEM_apu.^(1/options.h) .* rhs;
                        elseif options.OSL_COSEM == 1
                            C_osl(:, osa_iter) = im_vectors.OSEM_apu.^(1/options.h) .* rhs;
                        elseif options.OSL_COSEM == 2
                            C_osl(:, osa_iter) = im_vectors.OSEM_apu .* rhs;
                        end
                        if llo == 1 && ~options.compute_sensitivity_image
                            D = D + Summ;
                            if options.use_psf
                                Summ = computeConvolution(Summ, options, Nx, Ny, Nz, gaussK);
                            end
                            Summ(Summ < epps) = epps;
                            f_Summ(:,osa_iter) = Summ;
                        end
                    end
                    if ~options.compute_sensitivity_image
                        if options.use_psf
                            D = computeConvolution(D, options, Nx, Ny, Nz, gaussK);
                        end
                        D(D < epps) = epps;
                    else
                        D = sum(f_Summ,2);
                    end
                    if  options.PKMA
                        options.pj3 = D / options.subsets;
                    end
                end
                %%% RECONSTRUCTION PHASE %%%
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
                                    SinD = single(full(options.SinDelayed{llo}(pituus(osa_iter)+1:pituus(osa_iter + 1))));
                                else
                                    SinD = single(full(options.SinDelayed(pituus(osa_iter)+1:pituus(osa_iter + 1))));
                                end
                                SinD = SinD(:);
                            else
                                SinD = 0;
                            end
                            if normalization_correction
                                norm_input = single(options.normalization(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                            else
                                norm_input = 0;
                            end
                            if options.scatter_correction && ~options.subtract_scatter
                                scatter_input = double(ScatterC(pituus(osa_iter)+1:pituus(osa_iter + 1)));
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
                                    apux = x;
                                    apuy = y;
                                    apuz = z_det;
                                    x = reshape(x, numel(x)/2,2);
                                    y = reshape(y, numel(y)/2,2);
                                    z_det = reshape(z_det, numel(z_det)/2,2);
                                    x = x(pituus(osa_iter) + 1 : pituus(osa_iter + 1),:);
                                    y = y(pituus(osa_iter) + 1 : pituus(osa_iter + 1),:);
                                    z_det = z_det(pituus(osa_iter) + 1 : pituus(osa_iter + 1),:);
                                    x = x(:);
                                    y = y(:);
                                    z_det = z_det(:);
                                    det_per_ring = uint32(numel(x)/2);
                                end
                                xy_index_input = uint32(0);
                                z_index_input = uint32(0);
                                TOFSize = int64(size(L_input,1));
                                fullSize = size(LL,1);
                            else
                                L_input = uint16(0);
                                if options.CT && options.subsets == 1
                                    xy_index_input = xy_index;
                                    z_index_input = z_index;
                                else
                                    xy_index_input = xy_index(pituus(osa_iter)+1:pituus(osa_iter + 1));
                                    z_index_input = z_index(pituus(osa_iter)+1:pituus(osa_iter + 1));
                                end
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
                                    uu(1 + TOFSize * (dd - 1) : TOFSize * dd) = single(full(Sino(pituus(osa_iter) + 1 + fullSize * (dd - 1) : pituus(osa_iter + 1) + fullSize * (dd - 1))));
                                end
                            else
                                uu = single(full(Sino(pituus(osa_iter)+1:pituus(osa_iter + 1))));
                            end
                            uu(isnan(uu)) = 0;
                            uu(isinf(uu)) = 0;
                            
                            if options.use_psf
                                OSEM_apu = computeConvolution(im_vectors.OSEM_apu, options, Nx, Ny, Nz, gaussK);
                            else
                                OSEM_apu = im_vectors.OSEM_apu;
                            end
                            [Summ,rhs] = computeImplementation4(options,use_raw_data,randoms_correction, pituus,osa_iter, normalization_correction,...
                                Nx, Ny, Nz, dx, dy, dz, bx, by, bz, x, y, z_det, xx, yy, size_x, NSinos, NSlices, zmax, attenuation_correction, pseudot, det_per_ring, ...
                                TOF, TOFSize, sigma_x, TOFCenter, dec, nCores, L_input, lor_a_input, xy_index_input, z_index_input, epps, uu, OSEM_apu, no_norm, ...
                                x_center, y_center, z_center, bmin, bmax, Vmax, V, scatter_input, norm_input, SinD, dc_z);
                            if list_mode_format
                                x = apux;
                                y = apuy;
                                z_det = apuz;
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
                            end
                            [im_vectors,C_co,C_aco,C_osl] = computeEstimatesImp4(im_vectors, options, rhs, uu, f_Summ, SinD, gaussK, iter, osa_iter, C_co, C_aco,C_osl,...
                                randoms_correction, N, Ndx, Ndy, Ndz, D, tStart, epps, use_raw_data, pituus,normalization_correction,...
                                Nx, Ny, Nz, dx, dy, dz, bx, by, bz, x, y, z_det, xx, yy, size_x, NSinos, NSlices, zmax, attenuation_correction, pseudot, det_per_ring, ...
                                TOF, TOFSize, sigma_x, TOFCenter, dec, nCores, L_input, lor_a_input, xy_index_input, z_index_input,  ...
                                x_center, y_center, z_center, bmin, bmax, Vmax, V, scatter_input, norm_input, dc_z);
                            
                            clear Summ rhs
                        end
                        if options.save_iter
                            iter_n = iter + 1;
                        else
                            iter_n = 1;
                        end
                        im_vectors = computeEstimatesImp4Iter(im_vectors, options, gaussK, iter, iter_n, N, Ndx, Ndy, Ndz, epps, Nx, Ny, Nz, true, tStart_iter);
                        no_norm = true;
                    end
                    % Are MLEM-methods used?
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
                                SinD = single(full(options.SinDelayed{llo}));
                            else
                                SinD = single(full(options.SinDelayed));
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
                        
                        [Summ,rhs] = computeImplementation4(options,use_raw_data,randoms_correction, pituus, 1, normalization_correction,...
                            Nx, Ny, Nz, dx, dy, dz, bx, by, bz, x, y, z_det, xx, yy, size_x, NSinos, NSlices, zmax, attenuation_correction, pseudot, det_per_ring, ...
                            TOF, TOFSize, sigma_x, TOFCenter, dec, nCores, LL, lor_a, xy_index, z_index, epps, single(full(Sino)), MLEM_apu, no_norm, ...
                            x_center, y_center, z_center, bmin, bmax, Vmax, V, ScatterC, single(options.normalization), SinD, dc_z);
                        
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
                        im_vectors = computeEstimatesImp4Iter(im_vectors, options, gaussK, iter, iter_n, N, Ndx, Ndy, Ndz, epps, Nx, Ny, Nz, ...
                            false, tStart, f_Summ_ml, rhs);
                        no_norm = true;
                    end
                    if options.use_psf && options.deblurring && (options.save_iter || (~options.save_iter && iter == options.Niter))
                        im_vectors = computeDeblur(im_vectors, options, iter, gaussK, Nx, Ny, Nz);
                    end
                end
            end
            pause(0.01)
            if options.implementation ~= 2 && options.implementation ~= 3
                im_vectors = reshape_vectors(im_vectors, options);
            end
            pz = images_to_cell(im_vectors, llo, pz, options, rekot);
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
    elseif options.implementation == 2 || options.implementation == 3
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
            options.randSize = uint64(diff(pituus));
            if ~iscell(options.SinDelayed)
                options.SinDelayed = {single(full(options.SinDelayed))};
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
        else
            options.randSize = uint64(1);
            options.SinDelayed = {single(0)};
        end
        if ~isfield(options, 'normalization')
            options.normalization = single(0);
        end
        if ~isfield(options, 'scatter')
            options.scatter = false;
        end
        [pz] = computeImplementation23(options, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx, NSinos, NSlices, size_x, zmax, ...
            LL, pseudot, det_per_ring, TOF, sigma_x, TOFCenter, dec, device, use_raw_data, options.normalization, pituus, attenuation_correction, ...
            normalization_correction, Niter, subsets, epps, lor_a, xy_index, z_index, x_center, y_center, z_center, options.SinDelayed, ...
            options.SinM, bmin, bmax, Vmax, V, gaussK, 0, rekot, pz);
    else
        error('Unsupported reconstruction method.');
    end
end

% Save various image properties, e.g. matrix size, sinogram dimensions, FOV
% size, regularization parameters, etc.
pz = save_image_properties(options, pz, subsets);

end
