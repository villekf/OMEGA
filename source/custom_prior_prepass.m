function [options, varargout] = custom_prior_prepass(options, varargin)
%% Prepass phase for the custom prior file
% Computes all the necessary variables needed for the reconstruction
% process
%

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

if nargin > 1
    custom = false;
else
    custom = true;
end

if ~isfield(options,'use_Inveon')
    options.use_Inveon = 0;
end

if ~isfield(options,'save_iter')
    options.save_iter = true;
end


if ~isfield(options,'SinM')
    options.SinM = [];
end

if ~isfield(options,'TOF_bins') || options.TOF_bins == 0
    options.TOF_bins = 1;
end
if ~isfield(options,'TOF_width') || options.TOF_bins == 0
    options.TOF_width = 0;
end
TOF = options.TOF_bins > 1;

if options.precompute_lor == false && options.implementation == 3
    error('precompute_lor must be set to true if using method 3')
end

% if (options.MBSREM) && options.implementation == 2
%     error('BSREM, MBSREM and ROSEM-MAP do not work with implementation 2. Use implementation 1 instead')
% end

Nx = options.Nx;
Ny = options.Ny;
Nz = options.Nz;
subsets = options.subsets;
% attenuation_correction = options.attenuation_correction;
diameter = options.diameter;
FOVax = options.FOVa_x;
FOVay = options.FOVa_y;
axial_fov = options.axial_fov;
NSinos = options.NSinos;
pseudot = int32(options.pseudot);
rings = options.rings;
det_per_ring = options.det_per_ring;
machine_name = options.machine_name;
Nang = options.Nang;
Ndist = options.Ndist;
% cr_pz = options.cr_pz;
use_raw_data = options.use_raw_data;
TotSinos = options.TotSinos;
% attenuation_datafile = options.attenuation_datafile;
% partitions = options.partitions;
% verbose = options.verbose;

options.N = Nx * Ny * Nz;
% options.U = [];
% options.weights = [];
% options.a_L = [];
% options.fmh_weights = [];
% options.weighted_weights = [];
% options.weights_quad = [];

options.MAP = (options.OSL_MLEM || options.OSL_OSEM || options.BSREM || options.MBSREM || options.ROSEM_MAP || options.OSL_RBI || any(options.OSL_COSEM));
options.empty_weight = false;
if options.MBSREM || options.MRAMLA || options.RBI || options.OSL_RBI || options.COSEM || options.ECOSEM || options.ACOSEM || any(options.OSL_COSEM)
    options.MBSREM_prepass = true;
else
    options.MBSREM_prepass = false;
end

% if custom
options.rekot = reko_maker(options);
pz = cell(length(options.rekot),options.partitions);
% end

temp = pseudot;
if ~isempty(temp) && temp > 0
    for kk = uint32(1) : temp
        pseudot(kk) = uint32(options.cryst_per_block + 1) * kk;
    end
elseif temp == 0
    pseudot = [];
end
options.pseudot = pseudot;

if options.precompute_lor
    options.is_transposed = true;
else
    options.is_transposed = false;
end

[gaussK, options] = PSFKernel(options);

if custom || (isfield(options,'quad') && isfield(options,'Huber') && isfield(options,'FMH') && isfield(options,'L') && ...
        isfield(options,'weighted_mean') && isfield(options,'MRP') && isfield(options,'MAP'))
    if (options.quad || options.Huber || options.FMH || options.L || options.weighted_mean || options.MRP) && options.MAP
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
    if (options.quad || options.Huber || options.FMH || options.L || options.weighted_mean) && options.MAP
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
    
    options.im_vectors = form_image_vectors(options, options.N);
    
    if options.use_raw_data
        options.SinM = options.coincidences;
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
        % Perform corrections if needed
        if options.randoms_correction && ~options.reconstruct_trues && ~options.reconstruct_scatter
            if (options.use_ASCII || options.use_LMF || options.use_root) && options.use_machine == 0
                try
                    load(load_string, 'delayed_coincidences')
                    if ~options.corrections_during_reconstruction
                        if iscell(options.SinM) && iscell(delayed_coincidences)
                            for kk = 1 : length(options.SinM)
                                options.SinM{kk} = options.SinM{kk} - delayed_coincidences{kk};
                                options.SinM{kk}(options.SinM{kk} < 0) = 0;
                            end
                        else
                            options.SinM = options.SinM - delayed_coincidences;
                            options.SinM(options.SinM < 0) = 0;
                        end
                    end
                catch
                    options = loadDelayedData(options);
                    if ~options.corrections_during_reconstruction
                        if iscell(options.SinM) && iscell(options.SinDelayed)
                            for kk = 1 : length(options.SinM)
                                options.SinM{kk} = options.SinM{kk} - options.SinDelayed{kk};
                                options.SinM{kk}(options.SinM{kk} < 0) = 0;
                            end
                        else
                            options.SinM = options.SinM - options.SinDelayed;
                            options.SinM(options.SinM < 0) = 0;
                        end
                    end
                end
            else
                if iscell(options.SinD)
                    if iscell(options.SinM)
                        for kk = 1 : length(options.SinM)
                            if numel(options.SinD{kk}) ~= numel(options.SinM{kk})
                                error('Size mismatch between randoms correction data and measurement data')
                            end
                            options.SinM{kk} = options.SinM{kk} - single(options.SinD{kk}(:));
                            options.SinM{kk}(options.SinM{kk} < 0) = 0;
                        end
                    else
                        if numel(options.SinD{1}) ~= numel(options.SinM)
                            error('Size mismatch between randoms correction data and measurement data')
                        end
                        options.SinM = options.SinM - single(options.SinD{1}(:));
                        options.SinM(options.SinM < 0) = 0;
                    end
                else
                    if iscell(options.SinM)
                        for kk = 1 : length(options.SinM)
                            if numel(options.SinD) ~= numel(options.SinM{kk})
                                error('Size mismatch between randoms correction data and measurement data')
                            end
                            options.SinM{kk} = options.SinM{kk} - single(options.SinD(:));
                            options.SinM{kk}(options.SinM{kk} < 0) = 0;
                        end
                    else
                        if numel(options.SinD) ~= numel(options.SinM)
                            error('Size mismatch between randoms correction data and measurement data')
                        end
                        options.SinM = options.SinM - single(options.SinD(:));
                        options.SinM(options.SinM < 0) = 0;
                    end
                end
            end
        end
        if options.scatter_correction && ~options.corrections_during_reconstruction
            if iscell(options.ScatterC)
                if iscell(options.SinM)
                    for kk = 1 : length(options.SinM)
                        if numel(options.ScatterC{kk}) ~= numel(options.SinM{kk})
                            error('Size mismatch between scatter correction data and measurement data')
                        end
                        options.SinM{kk} = options.SinM{kk} - single(options.ScatterC{kk}(:));
                        options.SinM{kk}(options.SinM{kk} < 0) = 0;
                    end
                else
                    if numel(options.ScatterC{1}) ~= numel(options.SinM)
                        error('Size mismatch between scatter correction data and measurement data')
                    end
                    options.SinM = options.SinM - single(options.ScatterC{1}(:));
                    options.SinM(options.SinM < 0) = 0;
                end
            else
                if iscell(options.SinM)
                    for kk = 1 : length(options.SinM)
                        if numel(options.ScatterC) ~= numel(options.SinM{kk})
                            error('Size mismatch between scatter correction data and measurement data')
                        end
                        options.SinM{kk} = options.SinM{kk} - single(options.ScatterC(:));
                        options.SinM{kk}(options.SinM{kk} < 0) = 0;
                    end
                else
                    if numel(options.ScatterC) ~= numel(options.SinM)
                        error('Size mismatch between scatter correction data and measurement data')
                    end
                    options.SinM = options.SinM - single(options.ScatterC(:));
                    options.SinM(options.SinM < 0) = 0;
                end
            end
        end
        
        clear coincidences options.coincidences true_coincidences delayed_coincidences
        % Sinogram data
    else
        if options.partitions == 1
            load_string = [options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                num2str(options.TotSinos) '_span' num2str(options.span)];
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
            load_string = [options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_' ...
                num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                num2str(options.span)];
            load_string2 = [options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                num2str(options.span)];
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
        if options.partitions == 1 && options.randoms_correction && options.corrections_during_reconstruction
            if options.use_machine == 0 || options.use_machine == 1
                try
                    options.SinDelayed = loadStructFromFile(sinoFile,'SinDelayed');
                catch
                    options = loadDelayedData(options);
                end
            else
                [dfile, dfpath] = uigetfile('*.mat','Select delayed coincidence datafile');
                if isequal(dfile, 0)
                    error('No file was selected')
                end
                data = load(fullfile(dfpath, dfile));
                variables = fieldnames(data);
                if length(variables) > 1
                    if (any(strcmp('SinDelayed',variables)))
                        options.SinDelayed = double(data.(variables{strcmp('SinDelayed',variables)}));
                    else
                        error('The provided delayed coincidence file contains more than one variable and none of them are named "SinDelayed".')
                    end
                else
                    options.SinDelayed = single(data.(variables{1}));
                end
                clear data variables
            end
        elseif options.randoms_correction && options.corrections_during_reconstruction
            if options.use_machine == 0 || options.use_machine == 1
                try
                    options.SinDelayed = loadStructFromFile(sinoFile,'SinDelayed');
                catch
                    options = loadDelayedData(options);
                end
            else
                [options.file, options.fpath] = uigetfile('*.mat','Select delayed coincidence datafile');
                if isequal(options.file, 0)
                    error('No file was selected')
                end
                data = load(fullfile(options.fpath, options.file));
                variables = fieldnames(data);
                if length(variables) > 1
                    if (any(strcmp('SinDelayed',variables)))
                        options.SinDelayed = single(data.(variables{strcmp('SinDelayed',variables)}));
                    else
                        error('The provided delayed coincidence file contains more than one variable and none of them are named "SinDelayed".')
                    end
                else
                    options.SinDelayed = single(data.(variables{1}));
                end
                clear data variables
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
        options.sigma_x = 0;
        options.TOFCenter = 0;
        TOF = false;
        options.TOF_bins = 1;
    else
        c = 2.99792458e11;
        options.sigma_x = (c*options.TOF_FWHM/2) / (2 * sqrt(2 * log(2)));
        edges_user = linspace(-options.TOF_width * options.TOF_bins/2, options.TOF_width * options.TOF_bins / 2, options.TOF_bins + 1);
        edges_user = edges_user(1:end-1) + options.TOF_width/2; % the most probable value where annihilation occured
        TOFCenter = zeros(size(edges_user));
        TOFCenter(1) = edges_user(ceil(length(edges_user)/2));
        TOFCenter(2:2:end) = edges_user(ceil(length(edges_user)/2) + 1:end);
        TOFCenter(3:2:end) = edges_user(ceil(length(edges_user)/2) - 1: -1 : 1);
        options.TOFCenter = TOFCenter * c / 2;
        if isfield(options, 'TOF_offset') && options.TOF_offset > 0
            options.TOFCenter = options.TOFCenter + options.TOF_offset;
        end
    end
else
    options.sigma_x = 0;
    options.TOFCenter = 0;
end
if options.implementation == 2 || options.implementation == 3
    options.sigma_x = single(options.sigma_x);
    options.TOFCenter = single(options.TOFCenter);
end

%% This part is used when the observation matrix is calculated on-the-fly
% Compute the indices for the subsets used.
% For Sinogram data, five different methods to select the subsets are
% available. For raw list-mode data, three methods are available.

if (isfield(options,'x') && isfield(options,'y') && (isfield(options,'z') || isfield(options,'z_det'))) && numel(options.x) / 2 == numel(options.SinM)
    %     index = uint32(1:numel(options.SinM))';
    options.index = 0;
    det_per_ring = numel(options.SinM);
    options.det_per_ring = det_per_ring;
    options.pituus = floor(det_per_ring / subsets);
    options.pituus = int64([repmat(options.pituus,subsets - 1,1); det_per_ring - options.pituus*(subsets - 1)]);
    options.listmode = true;
%     if options.implementation == 1
%         error('List-mode reconstruction with custom detectors is currently not supported with implementation 1.')
%     end
    if options.implementation == 2 && options.use_CUDA
        error('CUDA support not enabled for list-mode data')
    end
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
    options.listmode = false;
    [options.index, options.pituus, options.subsets] = index_maker(Nx, Ny, Nz, subsets, use_raw_data, machine_name, options, Nang, Ndist, TotSinos, NSinos);
    Z = axial_fov;
    options.Z = Z;
end


%%

% Diameter of the PET-device (bore) (mm)
R=double(diameter);
% Transaxial FOV (x-direction, horizontal) (mm)
FOVax=double(FOVax);
% Transaxial FOV (y-direction, vertical) (mm)
FOVay=double(FOVay);
% Axial FOV (mm)
axial_fov = double(axial_fov);
% Number of rings
options.blocks=int32(rings + length(pseudot) - 1);
% First ring
options.block1=int32(0);

options.NSinos = int32(NSinos);
options.NSlices = int32(Nz);
options.TotSinos = int32(TotSinos);

[x, y, z_det, options] = get_coordinates(options, options.blocks, options.pseudot);

if ~custom && ~isfield(options,'randoms_correction')
    options.randoms_correction = false;
end
if ~custom && ~isfield(options,'scatter_correction')
    options.scatter_correction = false;
end

[normalization_correction, randoms_correction, options] = set_up_corrections(options, [], []);
options.normalization_correction = normalization_correction;
options.randoms_correction = randoms_correction;

if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
    options.x = single(x);
    options.y = single(y);
    options.z_det = single(z_det);
else
    options.x = double(x);
    options.y = double(y);
    options.z_det = double(z_det);
end

if options.use_raw_data
    options.size_x = uint32(options.det_w_pseudo);
else
    options.size_x = uint32(options.Nang*options.Ndist);
    if options.sampling > 1 && ~options.precompute_lor
        options.size_x = options.size_x * options.sampling;
    end
end

if options.precompute_lor && options.subsets > 1 || options.implementation == 2  && options.subsets > 1 || options.implementation == 4 && options.subsets > 1
    if exist('OCTAVE_VERSION','builtin') == 5
        options.pituus = [int64(0);cumsum(options.pituus,'native')];
    else
        options.pituus = [0;cumsum(options.pituus)];
    end
    if iscell(options.index)
        options.index = cell2mat(options.index);
    end
    indeksi = options.index;
elseif options.implementation == 1
    if exist('OCTAVE_VERSION','builtin') == 5
        options.pituus = [int64(0);cumsum(options.pituus,'native')];
    else
        options.pituus = [0;cumsum(options.pituus)];
    end
    if iscell(options.index)
        options.index = cell2mat(options.index);
    end
    indeksi = options.index;
else
    if iscell(options.index)
        options.index = cell2mat(options.index);
    end
    indeksi = options.index;
end

if options.use_raw_data
    if isempty(options.pseudot)
        options.pseudot = int32(1e5);
    else
        options.pseudot = options.pseudot - 1;
    end
end

% for the precomputed version, index vectors are needed

[options, options.lor_a, options.xy_index, options.z_index, options.LL, options.summa, options.pituus, options.SinM, options.lor_orth] = ...
    form_subset_indices(options, options.pituus, subsets, indeksi, options.size_x, options.y, options.z_det, rings, false, TOF, options.SinM);
if ~options.precompute_lor
    options.lor_a = uint16(0);
    options.lor_orth = uint16(0);
end

[options.xx,options.yy,options.zz,options.dx,options.dy,options.dz,options.bx,options.by,options.bz] = computePixelSize(R, FOVax, FOVay, Z, axial_fov, Nx, Ny, Nz, options.implementation);

options.zz=options.zz(2*options.block1+1:2*options.blocks+2);

% Number of pixels
options.Ny=int32(Ny);
options.Nx=int32(Nx);
options.Nz=int32(Nz);

options.N=(Nx)*(Ny)*(Nz);
options.det_per_ring = int32(det_per_ring);

% How much memory is preallocated
if use_raw_data == false
    options.ind_size = int32(NSinos/options.subsets*(det_per_ring)* Nx * (Ny));
else
    options.ind_size = int32((det_per_ring)^2/options.subsets* Nx * (Ny));
end


options.zmax = max(max(options.z_det));
if options.zmax==0
    if options.implementation == 2 || options.implementation == 3
        options.zmax = single(1);
    else
        options.zmax = double(1);
    end
end
[options.x_center,options.y_center,options.z_center,options.dec] = computePixelCenters(options.xx,options.yy,options.zz,options.dx,options.dy,options.dz,TOF,options);

[options.V,options.Vmax,options.bmin,options.bmax] = computeVoxelVolumes(options.dx,options.dy,options.dz,options);
%%

[options, options.D, options.C_co, options.C_aco, options.C_osl, options.Amin, options.E] = ...
    prepass_phase(options, options.pituus, options.index, options.SinM, options.pseudot, options.x, options.y, options.xx, options.yy, options.z_det, options.dz, options.dx, options.dy, ...
    options.bz, options.bx, options.by, options.NSlices, options.zmax, options.size_x, options.block1, options.blocks, options.normalization_correction, options.randoms_correction, ...
    options.xy_index, options.z_index, options.lor_a, options.lor_orth, options.summa, options.LL, options.is_transposed, options.x_center, options.y_center, options.z_center, 0, ...
    gaussK, options.bmin, options.bmax, options.Vmax, options.V);


if (options.MBSREM || options.MRAMLA) && (options.implementation == 1 || options.implementation == 4)
    if iscell(options.SinDelayed)
        if iscell(options.SinM)
            options.epsilon_mramla = MBSREM_epsilon(options.SinM{1}, options.epps, options.randoms_correction, options.SinDelayed{1}, options.E);
        else
            options.epsilon_mramla = MBSREM_epsilon(options.SinM, options.epps, options.randoms_correction, options.SinDelayed{1}, options.E);
        end
    else
        if iscell(options.SinM)
            options.epsilon_mramla = MBSREM_epsilon(options.SinM{1}, options.epps, options.randoms_correction, options.SinDelayed, options.E);
        else
            options.epsilon_mramla = MBSREM_epsilon(options.SinM, options.epps, options.randoms_correction, options.SinDelayed, options.E);
        end
    end
    if options.U == 0 || isempty(options.U)
        if iscell(options.SinM)
            options.U = max(double(options.SinM{1}) ./ options.Amin);
        else
            options.U = max(double(options.SinM) ./ options.Amin);
        end
    end
end
% if custom
varargout{1} = pz;
% end
end