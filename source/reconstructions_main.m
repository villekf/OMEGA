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

if options.use_raw_data
    if options.reconstruct_trues == false && ~options.reconstruct_scatter && (isfield(options, 'coincidences') == 0 && ~exist('coincidences','var') ...
            || options.precompute_all) && options.use_machine < 2
        if options.partitions == 1
            if options.use_ASCII && options.use_machine == 0
                load([options.machine_name '_measurements_' options.name '_static_raw_ASCII.mat'], 'coincidences')
            elseif options.use_LMF && options.use_machine == 0
                load([options.machine_name '_measurements_' options.name '_static_raw_LMF.mat'], 'coincidences')
            elseif options.use_root && options.use_machine == 0
                load([options.machine_name '_measurements_' options.name '_static_raw_root.mat'], 'coincidences')
            else
                load([options.machine_name '_measurements_' options.name '_static_raw_listmode.mat'], 'coincidences')
            end
        else
            if options.use_ASCII && options.use_machine == 0
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_ASCII.mat'], 'coincidences')
            elseif options.use_LMF && options.use_machine == 0
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_LMF.mat'], 'coincidences')
            elseif options.use_root && options.use_machine == 0
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_root.mat'], 'coincidences')
            else
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_listmode.mat'], 'coincidences')
            end
        end
        SinM = coincidences;
    elseif ~options.reconstruct_trues && ~options.reconstruct_scatter && ~exist('coincidences','var') && isfield(options, 'coincidences')
        SinM = options.coincidences;
    elseif options.reconstruct_trues
        if options.partitions == 1
            if options.use_ASCII
                load([options.machine_name '_measurements_' options.name '_static_raw_ASCII.mat'], 'true_coincidences')
            elseif options.use_LMF
                load([options.machine_name '_measurements_' options.name '_static_raw_LMF.mat'], 'true_coincidences')
            elseif options.use_root
                load([options.machine_name '_measurements_' options.name '_static_raw_root.mat'], 'true_coincidences')
            else
                load([options.machine_name '_measurements_' options.name '_static_raw_real.mat'], 'true_coincidences')
            end
        else
            if options.use_ASCII
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_ASCII.mat'], 'true_coincidences')
            elseif options.use_LMF
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_LMF.mat'], 'true_coincidences')
            elseif options.use_root
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_root.mat'], 'true_coincidences')
            else
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_real.mat'], 'true_coincidences')
            end
        end
        SinM = true_coincidences;
    elseif options.reconstruct_scatter
        if options.partitions == 1
            if options.use_ASCII
                load([options.machine_name '_measurements_' options.name '_static_raw_ASCII.mat'], 'scattered_coincidences')
            elseif options.use_LMF
                load([options.machine_name '_measurements_' options.name '_static_raw_LMF.mat'], 'scattered_coincidences')
            elseif options.use_root
                load([options.machine_name '_measurements_' options.name '_static_raw_root.mat'], 'scattered_coincidences')
            else
                load([options.machine_name '_measurements_' options.name '_static_raw_real.mat'], 'scattered_coincidences')
            end
        else
            if options.use_ASCII
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_ASCII.mat'], 'scattered_coincidences')
            elseif options.use_LMF
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_LMF.mat'], 'scattered_coincidences')
            elseif options.use_root
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_root.mat'], 'scattered_coincidences')
            else
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_real.mat'], 'scattered_coincidences')
            end
        end
        SinM = scattered_coincidences;
    end
    if options.randoms_correction && ~options.reconstruct_trues && ~options.reconstruct_scatter
        if (options.use_ASCII || options.use_LMF || options.use_root) && options.use_machine == 0
            if options.partitions == 1
                if options.use_ASCII && options.use_machine == 0
                    load([options.machine_name '_measurements_' options.name '_static_raw_ASCII.mat'], 'delayed_coincidences')
                elseif options.use_LMF && options.use_machine == 0
                    load([options.machine_name '_measurements_' options.name '_static_raw_LMF.mat'], 'delayed_coincidences')
                elseif options.use_root && options.use_machine == 0
                    load([options.machine_name '_measurements_' options.name '_static_raw_root.mat'], 'delayed_coincidences')
                else
                    load([options.machine_name '_measurements_' options.name '_static_raw_listmode.mat'], 'delayed_coincidences')
                end
            else
                if options.use_ASCII && options.use_machine == 0
                    load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                        num2str(options.tot_time) 's_raw_ASCII.mat'], 'delayed_coincidences')
                elseif options.use_LMF && options.use_machine == 0
                    load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                        num2str(options.tot_time) 's_raw_LMF.mat'], 'delayed_coincidences')
                elseif options.use_root && options.use_machine == 0
                    load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                        num2str(options.tot_time) 's_raw_root.mat'], 'delayed_coincidences')
                else
                    load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                        num2str(options.tot_time) 's_raw_listmode.mat'], 'delayed_coincidences')
                end
            end
            if exist('delayed_coincidences','var')
                if ~options.corrections_during_reconstruction
                    if iscell(SinM) && iscell(delayed_coincidences)
                        for kk = 1 : length(SinM)
                            SinM{kk} = SinM{kk} - delayed_coincidences{kk};
                            SinM{kk}(SinM{kk} < 0) = 0;
                        end
                    else
                        SinM = SinM - delayed_coincidences;
                        SinM(SinM < 0) = 0;
                    end
                end
            else
                disp('Delayed coincidences not found, randoms correction not performed')
            end
        else
            if iscell(options.SinD)
                if iscell(SinM)
                    for kk = 1 : length(SinM)
                        if numel(options.SinD{kk}) ~= numel(SinM{kk})
                            error('Size mismatch between randoms correction data and measurement data')
                        end
                        SinM{kk} = SinM{kk} - double(options.SinD{kk}(:));
                        SinM{kk}(SinM{kk} < 0) = 0;
                    end
                else
                    if numel(options.SinD{1}) ~= numel(SinM)
                        error('Size mismatch between randoms correction data and measurement data')
                    end
                    SinM = SinM - double(options.SinD{1}(:));
                    SinM(SinM < 0) = 0;
                end
            else
                if iscell(SinM)
                    for kk = 1 : length(SinM)
                        if numel(options.SinD) ~= numel(SinM{kk})
                            error('Size mismatch between randoms correction data and measurement data')
                        end
                        SinM{kk} = SinM{kk} - double(options.SinD(:));
                        SinM{kk}(SinM{kk} < 0) = 0;
                    end
                else
                    if numel(options.SinD) ~= numel(SinM)
                        error('Size mismatch between randoms correction data and measurement data')
                    end
                    SinM = SinM - double(options.SinD(:));
                    SinM(SinM < 0) = 0;
                end
            end
        end
    end
    if options.scatter_correction
        if iscell(options.ScatterC)
            if iscell(SinM)
                for kk = 1 : length(SinM)
                    if numel(options.ScatterC{kk}) ~= numel(SinM{kk})
                        error('Size mismatch between scatter correction data and measurement data')
                    end
                    SinM{kk} = SinM{kk} - double(options.ScatterC{kk}(:));
                    SinM{kk}(SinM{kk} < 0) = 0;
                end
            else
                if numel(options.ScatterC{1}) ~= numel(SinM)
                    error('Size mismatch between scatter correction data and measurement data')
                end
                SinM = SinM - double(options.ScatterC{1}(:));
                SinM(SinM < 0) = 0;
            end
        else
            if iscell(SinM)
                for kk = 1 : length(SinM)
                    if numel(options.ScatterC) ~= numel(SinM{kk})
                        error('Size mismatch between scatter correction data and measurement data')
                    end
                    SinM{kk} = SinM{kk} - double(options.ScatterC(:));
                    SinM{kk}(SinM{kk} < 0) = 0;
                end
            else
                if numel(options.ScatterC) ~= numel(SinM)
                    error('Size mismatch between scatter correction data and measurement data')
                end
                SinM = SinM - double(options.ScatterC(:));
                SinM(SinM < 0) = 0;
            end
        end
    end
    
    clear coincidences options.coincidences true_coincidences delayed_coincidences
else
    if (~options.reconstruct_trues && ~options.reconstruct_scatter) || options.use_machine > 0
        if options.partitions == 1 && isfield(options, 'SinM') == 0
            if (options.randoms_correction || options.scatter_correction || options.normalization_correction) && ~options.corrections_during_reconstruction && options.use_machine < 2
                if options.use_machine == 0
                    load([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                        num2str(options.TotSinos) '_span' num2str(options.span) '.mat'],'SinM')
                elseif  options.use_machine == 1
                    load([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                        num2str(options.TotSinos) '_span' num2str(options.span) '_listmode.mat'],'SinM')
                end
            elseif options.use_machine < 2
                if options.use_machine == 0
                    load([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                        num2str(options.TotSinos) '_span' num2str(options.span) '.mat'],'raw_SinM')
                elseif  options.use_machine == 1
                    load([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                        num2str(options.TotSinos) '_span' num2str(options.span) '_listmode.mat'],'raw_SinM')
                end
                SinM = raw_SinM;
                clear raw_SinM
            else
                load([options.machine_name '_' options.name '_sinogram_original_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                    num2str(options.TotSinos) '_span' num2str(options.span) '_machine_sinogram.mat'],'SinM')
            end
        elseif isfield(options, 'SinM') == 0
            %             variableInfo = who('-file', [options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) ...
            %                 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
            %                 num2str(options.TotSinos) '_span' num2str(options.span) '.mat']);
            if (options.randoms_correction || options.scatter_correction || options.normalization_correction) && ~options.corrections_during_reconstruction && options.use_machine < 2
                if options.use_machine == 0
                    load([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                        num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                        num2str(options.span) '.mat'], 'SinM')
                elseif  options.use_machine == 1
                    load([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                        num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                        num2str(options.span) '_listmode.mat'], 'SinM')
                end
            elseif options.use_machine < 2
                if options.use_machine == 0
                    load([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ '
                        num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                        num2str(options.span) '.mat'], 'raw_SinM')
                elseif  options.use_machine == 1
                    load([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                        num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                        num2str(options.span) '_listmode.mat'], 'raw_SinM')
                end
                SinM = raw_SinM;
                clear raw_SinM
            else
                load([options.machine_name '_' options.name '_sinograms_original_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                    num2str(options.span) '_machine_sinogram.mat'], 'SinM')
            end
        else
            SinM = options.SinM;
            clear options.SinM
        end
    elseif options.reconstruct_trues && options.use_machine == 0
        if options.partitions == 1 && isfield(options, 'SinTrues') == 0
            load([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                num2str(options.TotSinos) '_span' num2str(options.span) '.mat'],'SinTrues')
        elseif isfield(options, 'SinTrues') == 0
            load([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                num2str(options.span) '.mat'], 'SinTrues')
        end
        SinM = SinTrues;
        clear SinTrues
    elseif options.reconstruct_scatter && options.use_machine == 0
        if options.partitions == 1 && isfield(options, 'SinScatter') == 0
            load([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                num2str(options.TotSinos) '_span' num2str(options.span) '.mat'],'SinScatter')
        elseif isfield(options, 'SinScatter') == 0
            load([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                num2str(options.span) '.mat'], 'SinScatter')
        end
        SinM = SinScatter;
        clear SinScatter
    end
    if options.partitions == 1 && options.randoms_correction && options.corrections_during_reconstruction
        if options.use_machine == 0
            load([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                num2str(options.TotSinos) '_span' num2str(options.span) '.mat'],'SinDelayed')
        elseif  options.use_machine == 1
            load([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                num2str(options.TotSinos) '_span' num2str(options.span) '_listmode.mat'],'SinDelayed')
        else
            [options.file, options.fpath] = uigetfile('*.mat','Select delayed coincidence datafile');
            if isequal(options.file, 0)
                error('No file was selected')
            end
            data = load(fullfile(options.fpath, options.file));
            variables = fields(data);
            if length(variables) > 1
                if (any(strcmp('SinDelayed',variables)))
                    SinDelayed = double(data.(variables{strcmp('SinDelayed',variables)}));
                else
                    error('The provided delayed coincidence file contains more than one variable and none of them are named "SinDelayed".')
                end
            else
                SinDelayed = double(data.(variables{1}));
            end
            clear data variables
        end
    elseif options.randoms_correction && options.corrections_during_reconstruction
        if options.use_machine == 0
            load([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                num2str(options.span) '.mat'], 'SinDelayed')
        elseif  options.use_machine == 1
            load([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                num2str(options.span) '_listmode.mat'], 'SinDelayed')
        else
            [options.file, options.fpath] = uigetfile('*.mat','Select delayed coincidence datafile');
            if isequal(options.file, 0)
                error('No file was selected')
            end
            data = load(fullfile(options.fpath, options.file));
            variables = fields(data);
            if length(variables) > 1
                if (any(strcmp('SinDelayed',variables)))
                    SinDelayed = double(data.(variables{strcmp('SinDelayed',variables)}));
                else
                    error('The provided delayed coincidence file contains more than one variable and none of them are named "SinDelayed".')
                end
            else
                SinDelayed = double(data.(variables{1}));
            end
            clear data variables
            if length(SinDelayed) < options.partitions && iscell(SinDelayed)
                warning('The number of delayed coincidence sinograms is less than the number of time points. Using the first one')
                temp = SinDelayed;
                SinDelayed = cell(options.partitions,1);
                if sum(size(temp{1})) > 1
                    if size(temp{1},1) ~= size(options.Nang)
                        temp{1} = permute(temp{1},[2 1 3]);
                    end
                end
                for kk = 1 : options.partitions
                    SinDelayed{kk} = temp{1};
                end
            elseif options.partitions > 1
                warning('The number of delayed coincidence sinograms is less than the number of time points. Using the first one')
                temp = SinDelayed;
                SinDelayed = cell(options.partitions,1);
                if sum(size(temp)) > 1
                    if size(temp,1) ~= size(options.Nang)
                        temp = permute(temp,[2 1 3]);
                    end
                end
                for kk = 1 : options.partitions
                    SinDelayed{kk} = temp;
                end
            else
                if iscell(SinDelayed)
                    for kk = 1 : length(SinDelayed)
                        if sum(size(SinDelayed{kk})) > 1
                            if size(SinDelayed{kk},1) ~= size(options.Nang)
                                SinDelayed{kk} = permute(SinDelayed{kk},[2 1 3]);
                            end
                        end
                    end
                else
                    if sum(size(SinDelayed)) > 1
                        if size(SinDelayed,1) ~= size(options.Nang)
                            SinDelayed = permute(SinDelayed,[2 1 3]);
                        end
                    end
                end
            end
        end
    end
end

% if options.precompute_lor == false && options.reconstruction_method == 4
%     error('precompute_lor must be set to true if using method 4')
% end

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
det_per_ring = options.det_per_ring;
machine_name = options.machine_name;
Nang = options.Nang;
Ndist = options.Ndist;
cr_pz = options.cr_pz;
use_raw_data = options.use_raw_data;
TotSinos = options.TotSinos;
attenuation_datafile = options.attenuation_datafile;
partitions = options.partitions;
verbose = options.verbose;
device = uint32(options.use_device);
empty_weight = false;
options.MBSREM_prepass = true;
tr_offsets = 0;
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

pz = cell(length(rekot) + 1,partitions);

N = Nx * Ny * Nz;

temp = pseudot;
if ~isempty(temp) && temp > 0
    for kk = uint32(1) : temp
        pseudot(kk) = uint32(options.cryst_per_block + 1) * kk;
    end
elseif temp == 0
    pseudot = [];
end



if options.reconstruction_method == 3
    if options.osem && options.mlem
        if subsets == 1
            disp(['Both OSEM and MLEM set for method ' num2str(options.reconstruction_method) ', using MLEM'])
            options.osem = false;
        else
            disp(['Both OSEM and MLEM set for method ' num2str(options.reconstruction_method) ', using OSEM'])
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
        empty_weight = true;
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

if options.reconstruction_method == 1 || options.reconstruction_method == 4
    im_vectors = form_image_vectors(options, N);
end

if options.precompute_lor
    is_transposed = true;
else
    is_transposed = false;
end

if ~any(sum(rekot(11:end))) && options.reconstruction_method == 3
    subsets = 1;
    %     pituus = [pituus(1); pituus(end)];
end

%%
% for llo = 1 : partitions
if iscell(SinM)
    Sino = SinM{1};
else
    Sino = SinM;
end

Sino = Sino(:);

if issparse(Sino)
    Sino = (full(Sino));
end

if use_raw_data == false && NSinos ~= TotSinos
    Sino = Sino(1:NSinos*Ndist*Nang);
end

%% This computes a whole observation matrix and uses it to compute the MLEM (no on-the-fly calculations)
if precompute_obs_matrix && options.reconstruction_method == 1
    
    for llo = 1 : partitions
        
        
        if options.mlem
            if iscell(SinM)
                Sino = SinM{llo};
            else
                Sino = SinM;
                clear SinM
            end
            
            Sino = Sino(:);
            
            if issparse(Sino)
                Sino = (full(Sino));
            end
            
            if use_raw_data == false && NSinos ~= TotSinos
                Sino = Sino(1:NSinos*Ndist*Nang);
            end
            if use_raw_data == false && options.precompute_lor
                lor_file = [folder machine_name '_lor_pixel_count_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_sino_' num2str(Ndist) 'x' ...
                    num2str(Nang) 'x' num2str(TotSinos) '.mat'];
                if exist(lor_file, 'file') == 2
                    load(lor_file,'lor')
                else
                    lor_pixel_count_prepass(options);
                    load(lor_file,'lor')
                end
                if NSinos ~= TotSinos
                    lor = lor(1:NSinos*Ndist*Nang);
                end
                Sino = Sino(lor > 0);
                clear lor
            end
            
            pituus = uint32(length(Sino));
            if llo == 1
                A = observation_matrix_formation(diameter, options.FOVa_x, axial_fov, rings, pseudot, Nx, Ny, Nz, det_per_ring, cr_pz, ...
                    options.use_fsparse, attenuation_correction, attenuation_datafile, options.precompute_lor, options.use_raw_data, pituus, ...
                    options, Nz);
                D = full(sum(A,2));
                D(D <= 0) = epps;
            end
            
            %             ll = ll + 1;
            
            MLEM = ones(N,Niter);
            MLEM(:,1) = options.x0(:);
            Sino = double(Sino);
            
            for iter=1:Niter
                if options.mlem
                    MLEM(:,iter+1) = MLEM_im(MLEM(:,iter), D, epps, A, Sino);
                    disp(['MLEM iteration ' num2str(iter) ' finished'])
                end
                if sum(rekot(2:end)) > 0
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
    
    % Compute the indices for the subsets used.
    % For Sinogram data, five different methods to select the subsets are
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
    
    if use_raw_data == false
        sino_file = [folder machine_name '_app_coordinates_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'];
        if exist(sino_file, 'file') == 2
            variableInfo = who('-file', sino_file);
            if ismember('x', variableInfo) && ismember('y', variableInfo)
                load(sino_file,'x', 'y');
            else
                sinogram_coordinates_2D(options);
                load(sino_file,'x', 'y');
            end
        else
            sinogram_coordinates_2D(options);
            load(sino_file,'x', 'y');
        end
        z = sinogram_coordinates_3D(options);
        if NSinos ~= TotSinos
            z = z(1:NSinos,:);
        end
    else
        %         load([machine_name '_detector_coordinates.mat'],'x','y');
        [x, y, ~, ~] = detector_coordinates(options);
        
        z_length = double(blocks) * cr_pz;
        z = linspace(0, z_length, blocks + 1);
    end
    if attenuation_correction
        data = load(attenuation_datafile);
        variables = fields(data);
        vaimennus = double(data.(variables{1}));
        if size(vaimennus,1) ~= Nx || size(vaimennus,2) ~= Ny || size(vaimennus,3) ~= Nz
            if size(vaimennus,1) ~= Nx*Ny*Nz
                error("Error: Attenuation data is of different size than the reconstructed image")
            end
        end
        if size(vaimennus,2) == 1
            vaimennus = vaimennus(:,:,2*block1+1:2*blocks+1);
        else
            vaimennus = vaimennus(2*block1+1:(2*blocks+1)*Nx*Ny);
        end
        vaimennus = vaimennus(:) / 10;
        if options.reconstruction_method == 2 || options.reconstruction_method == 3 || options.reconstruction_method == 5
            vaimennus = single(vaimennus);
        end
        clear data
    else
        if options.reconstruction_method == 2 || options.reconstruction_method == 3 || options.reconstruction_method == 5
            vaimennus = single(0);
        else
            vaimennus = 0;
        end
    end
    
    if (options.normalization_correction && options.corrections_during_reconstruction) && ~options.use_user_normalization
        normalization_correction = true;
        
    elseif options.normalization_correction && options.use_user_normalization && options.corrections_during_reconstruction
        normalization_correction = true;
        [options.file, options.fpath] = uigetfile({'*.nrm;*.mat'},'Select normalization datafile');
        if isequal(options.file, 0)
            error('No file was selected')
        end
        if any(strfind(options.file, '.nrm'))
            fid = fopen(options.file);
            normalization = fread(fid, inf, 'single=>single',0,'l');
            fclose(fid);
        else
            data = load(options.file);
            variables = fields(data);
            normalization = data.(variables{1});
            clear data
            if numel(normalization) ~= options.Ndist*options.Nang*options.TotSinos
                error('Size mismatch between the current data and the normalization data file')
            end
        end
        normalization = normalization(:);
        %         Sino = Sino .* normalization;
        if (options.reconstruction_method == 1 || options.reconstruction_method == 4)
            normalization = double(normalization);
            normalization = 1 ./ normalization;
        else
            normalization = single(1) ./ normalization;
        end
    else
        normalization_correction = false;
        if options.reconstruction_method == 2 || options.reconstruction_method == 3 || options.reconstruction_method == 5
            normalization = single(0);
        else
            normalization = 0;
        end
    end
    
    if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
        randoms_correction = true;
        if iscell(SinDelayed)
            for k = 1 : options.partitions
                SinM{kk} = SinM{kk} + SinDelayed{kk};
                SinDelayed{kk} = 2 * SinDelayed{kk};
                if options.scatter_correction
                    if sum(size(options.ScatterC)) > 1 && ~iscell(options.ScatterC)
                        if size(options.ScatterC,1) ~= size(options.Nang)
                            ScatterC = permute(options.ScatterC,[2 1 3]);
                        end
                    elseif iscell(options.ScatterC)
                        if sum(size(options.ScatterC{1})) > 1
                            if size(options.ScatterC{1},1) ~= size(options.Nang)
                                if length(options.ScatterC) > 1
                                    ScatterC = permute(options.ScatterC{kk},[2 1 3]);
                                else
                                    ScatterC = permute(options.ScatterC{1},[2 1 3]);
                                end
                            end
                        end
                    end
                    SinDelayed{kk} = SinDelayed{kk} + ScatterC(:);
                end
            end
        else
            SinDelayed = SinDelayed(:);
            if use_raw_data == false && NSinos ~= TotSinos
                SinDelayed = SinDelayed(1:NSinos*Ndist*Nang);
            end
            Sino = Sino + SinDelayed;
            SinDelayed = 2 * SinDelayed;
            if options.scatter_correction
                if sum(size(options.ScatterC)) > 1 && ~iscell(options.ScatterC)
                    if size(options.ScatterC,1) ~= size(options.Nang)
                        options.ScatterC = permute(options.ScatterC,[2 1 3]);
                    end
                elseif iscell(options.ScatterC)
                    if sum(size(options.ScatterC{1})) > 1
                        if size(options.ScatterC{1},1) ~= size(options.Nang)
                            options.ScatterC = permute(options.ScatterC{1},[2 1 3]);
                        end
                    end
                end
                SinDelayed = SinDelayed + options.ScatterC(:);
            end
        end
    else
        randoms_correction = false;
        if options.reconstruction_method == 2 || options.reconstruction_method == 3 || options.reconstruction_method == 5
            SinDelayed = {single(0)};
        else
            SinDelayed = {0};
        end
    end
    
    if min(min(z)) == 0
        z = z + (axial_fow - max(max(z)))/2;
    end
    
    if options.reconstruction_method == 2 || options.reconstruction_method == 3 || options.reconstruction_method == 5
        x=single(x);
        y=single(y);
        z_det = single(z);
    else
        x=double(x);
        y=double(y);
        z_det = double(z);
    end
    clear z
    
    
    size_x = uint32(size(x,1));
    
    if subsets > 1 && ((options.precompute_lor && options.reconstruction_method == 1) || options.reconstruction_method > 1)
        pituus = [0;cumsum(pituus)];
        if iscell(index)
            index = cell2mat(index);
        end
    end
    if ~options.precompute_lor
        lor_a = uint16(0);
    end
    
    %     if ~options.precompute_lor && options.reconstruction_method > 1
    %         index = cell2mat(index);
    %     end
    
    if use_raw_data
        if isempty(pseudot)
            pseudot = uint32(1e5);
        else
            pseudot = pseudot - 1;
        end
    end
    
    % for the precomputed version, index vectors are needed
    if use_raw_data == false && options.precompute_lor
        
        lor_file = [folder machine_name '_lor_pixel_count_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_sino_' num2str(Ndist) 'x' ...
            num2str(Nang) 'x' num2str(TotSinos) '.mat'];
        
        if exist(lor_file, 'file') == 2
            if options.reconstruction_method == 1 || options.reconstruction_method == 4
                variableInfo = who('-file', lor_file);
                if any(ismember('lor', variableInfo))
                    load(lor_file,'lor')
                else
                    lor_pixel_count_prepass(options);
                    load(lor_file,'lor')
                end
                if options.projector_type == 2 && options.reconstruction_method == 1
                    if any(ismember('lor_orth', variableInfo))
                        load(lor_file,'lor_orth')
                    else
                        lor_pixel_count_prepass(options);
                        load(lor_file,'lor_orth')
                    end
                    load(lor_file,'crystal_size_z')
                    load(lor_file,'crystal_size_xy')
                    if options.tube_width_z == 0
                        if crystal_size_xy ~= options.tube_width_xy
                            lor_pixel_count_prepass(options);
                            load(lor_file,'lor_orth')
                        end
                        if length(lor_orth) > length(Sino)
                            lor_orth = lor_orth(1:length(Sino));
                        end
                    elseif options.tube_width_z > 0
                        if crystal_size_z ~= options.tube_width_z
                            lor_pixel_count_prepass(options);
                            load(lor_file,'lor_orth')
                        end
                        if length(lor_orth) > length(Sino)
                            lor_orth = lor_orth(length(Sino)+1:end);
                        elseif length(lor_orth) == length(Sino)
                            lor_pixel_count_prepass(options);
                            load(lor_file,'lor_orth')
                            lor_orth = lor_orth(length(Sino)+1:end);
                        end
                    end
                end
            else
                variableInfo = who('-file', lor_file);
                if any(ismember('lor_opencl', variableInfo))
                    load(lor_file,'lor_opencl')
                else
                    lor_pixel_count_prepass(options);
                    load(lor_file,'lor_opencl')
                end
                lor = lor_opencl;
                clear lor_opencl
            end
        else
            lor_pixel_count_prepass(options);
            if options.reconstruction_method == 1 || options.reconstruction_method == 4
                load(lor_file,'lor')
                if options.projector_type == 2 && options.reconstruction_method == 1
                    load(lor_file,'lor_orth')
                    if options.tube_width_z == 0
                        if length(lor_orth) > length(Sino)
                            lor_orth = lor_orth(1:length(Sino));
                        end
                    elseif options.tube_width_z > 0
                        if length(lor_orth) > length(Sino)
                            lor_orth = lor_orth(length(Sino)+1:end);
                        end
                    end
                end
            else
                load(lor_file,'lor_opencl')
                lor = lor_opencl;
                clear lor_opencl
            end
        end
        if subsets > 1
            lor_a = (lor(index));
            if options.projector_type == 2 && options.reconstruction_method == 1
                lor_orth = (lor_orth(index));
            end
            Sino = Sino(index);
            if options.normalization_correction && options.corrections_during_reconstruction
                normalization = normalization(index);
            end
            if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
                if partitions > 1
                    for ff = 1 : partitions
                        temp = SinDelayed{ff};
                        if NSinos ~= TotSinos
                            temp = temp(:,:,1:NSinos);
                        end
                        temp = temp(index);
                        SinDelayed{ff} = temp;
                    end
                    clear temp
                else
                    SinDelayed = SinDelayed(index);
                end
            end
            if partitions > 1
                for ff = 1 : partitions
                    temp = SinM{ff};
                    if NSinos ~= TotSinos
                        temp = temp(:,:,1:NSinos);
                    end
                    temp = temp(index);
                    SinM{ff} = temp;
                end
                clear temp
            end
            clear lor
        else
            discard = lor > 0;
            if length(discard) ~= TotSinos*Ndist*Nang
                error('Error: Size mismatch between sinogram and LORs to be removed')
            end
            if use_raw_data == false && NSinos ~= TotSinos
                discard = discard(1:NSinos*Ndist*Nang);
            end
            lor_a = (lor(discard));
            if options.projector_type == 2 && options.reconstruction_method == 1
                if use_raw_data == false && NSinos ~= TotSinos
                    lor_orth = lor_orth(1:NSinos*Ndist*Nang);
                end
                lor_orth = (lor_orth(discard));
            end
            Sino = Sino(discard);
            if options.normalization_correction && options.corrections_during_reconstruction
                if use_raw_data == false && NSinos ~= TotSinos
                    normalization = normalization(1:NSinos*Ndist*Nang);
                end
                normalization = normalization(discard);
            end
            if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
                if partitions > 1
                    for ff = 1 : partitions
                        temp = SinDelayed{ff};
                        if NSinos ~= TotSinos
                            temp = temp(:,:,1:NSinos);
                        end
                        temp = temp(discard);
                        SinDelayed{ff} = temp;
                    end
                    clear temp
                else
                    SinDelayed = SinDelayed(discard);
                end
            end
            if partitions > 1
                for ff = 1 : partitions
                    temp = SinM{ff};
                    if NSinos ~= TotSinos
                        temp = temp(:,:,1:NSinos);
                    end
                    temp = temp(discard);
                    SinM{ff} = temp;
                end
                clear temp
            end
            clear lor
        end
        [~, I] = sort(y, 2);
        sy = size(y);
        I = sub2ind(sy, repmat((1:sy(1)).', 1, sy(2)), I);
        
        xy_index = uint32(I(:,1));
        xy_index2 = repmat(uint32(1:size_x)', NSinos - Nz, 1);
        xy_index = [repmat(xy_index, Nz, 1); xy_index2];
        [~, I] = sort(z_det, 2);
        sy = size(z_det);
        I = sub2ind(sy, repmat((1:sy(1)).', 1, sy(2)), I);
        
        z_index = uint16(I(:,1));
        if verLessThan('matlab','8.5')
            z_index = repeat_elem(z_index, size_x);
        else
            z_index = repelem(z_index, size_x);
        end
        if subsets > 1
            z_index = z_index(index);
        else
            z_index = (z_index(discard));
        end
        apu = z_index > NSinos;
        z_index = z_index - 1;
        
        if subsets > 1
            xy_index = xy_index(index);
        else
            xy_index = (xy_index(discard));
        end
        xy_index(apu) = xy_index(apu) + uint32(size_x);
        xy_index = xy_index - 1;
        
        summa = zeros(subsets, 1, 'uint64');
        
        if subsets > 1
            for kk = 1 : subsets
                if options.projector_type == 2 && options.reconstruction_method == 1
                    summa(kk) = sum(int64(lor_orth(pituus(kk)+1:pituus(kk+1))));
                else
                    summa(kk) = sum(int64(lor_a(pituus(kk)+1:pituus(kk+1))));
                end
            end
        else
            if options.projector_type == 2 && options.reconstruction_method == 1
                summa = uint64(sum(int64(lor_orth)));
            else
                summa = uint64(sum(int64(lor_a)));
            end
            pituus = uint32([0;length(Sino)]);
        end
        
        
        clear discard I yt xt xy_index2 index apu
    elseif use_raw_data && options.precompute_lor
        
        lor_file = [folder machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'];
        if exist(lor_file, 'file') == 2
            variableInfo = who('-file', lor_file);
            if options.reconstruction_method == 1 || options.reconstruction_method == 4
                if ismember('lor', variableInfo)
                    load(lor_file,'lor')
                else
                    lor_pixel_count_prepass(options);
                    load(lor_file,'lor')
                end
                if options.projector_type == 2 && options.reconstruction_method == 1
                    if any(ismember('lor_orth', variableInfo))
                        load(lor_file,'lor_orth')
                    else
                        lor_pixel_count_prepass(options);
                        load(lor_file,'lor_orth')
                    end
                    load(lor_file,'crystal_size_z')
                    load(lor_file,'crystal_size_xy')
                    if options.tube_width_z == 0
                        if crystal_size_xy ~= options.tube_width_xy
                            lor_pixel_count_prepass(options);
                            load(lor_file,'lor_orth')
                        end
                        if length(lor_orth) > length(Sino)
                            lor_orth = lor_orth(1:length(Sino));
                        end
                    elseif options.tube_width_z > 0
                        if crystal_size_z ~= options.tube_width_z
                            lor_pixel_count_prepass(options);
                            load(lor_file,'lor_orth')
                        end
                        if length(lor_orth) > length(Sino)
                            lor_orth = lor_orth(length(Sino)+1:end);
                        elseif length(lor_orth) == length(Sino)
                            lor_pixel_count_prepass(options);
                            load(lor_file,'lor_orth')
                            lor_orth = lor_orth(length(Sino)+1:end);
                        end
                    end
                end
            else
                if ismember('lor_opencl', variableInfo)
                    load(lor_file,'lor_opencl')
                else
                    lor_pixel_count_prepass(options);
                    load(lor_file,'lor_opencl')
                end
                lor = lor_opencl;
                clear lor_opencl
            end
        else
            lor_pixel_count_prepass(options);
            if options.reconstruction_method == 1 || options.reconstruction_method == 4
                load(lor_file,'lor')
                if options.projector_type == 2 && options.reconstruction_method == 1
                    load(lor_file,'lor_orth')
                    if options.tube_width_z == 0
                        if length(lor_orth) > length(Sino)
                            lor_orth = lor_orth(1:length(Sino));
                        end
                    elseif options.tube_width_z > 0
                        if length(lor_orth) > length(Sino)
                            lor_orth = lor_orth(length(Sino)+1:end);
                        end
                    end
                end
            else
                load(lor_file,'lor_opencl')
                lor = lor_opencl;
                clear lor_opencl
            end
        end
        discard = lor > 0;
        if subsets > 1
            if ~exist('LL','var')
                LL = form_detector_pairs_raw(rings, det_per_ring);
            end
            LL = LL(discard,:);
            Sino = Sino(discard);
            lor = lor(discard);
            if options.normalization_correction && options.corrections_during_reconstruction
                normalization = normalization(discard);
            end
            if options.projector_type == 2 && options.reconstruction_method == 1
                lor_orth = (lor_orth(discard));
            end
            if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
                if partitions > 1
                    for ff = 1 : partitions
                        temp = SinDelayed{ff};
                        temp = temp(discard);
                        SinDelayed{ff} = temp;
                    end
                    clear temp
                else
                    SinDelayed = SinDelayed(discard);
                end
            end
            LL = LL(index,:);
            lor_a = (lor(index));
            Sino = Sino(index);
            if partitions > 1
                for ff = 1 : partitions
                    temp = SinM{ff};
                    temp = temp(discard);
                    if subsets > 1
                        temp = temp(index);
                    end
                    SinM{ff} = temp;
                end
                clear temp
            end
            if options.normalization_correction && options.corrections_during_reconstruction
                normalization = normalization(index);
            end
            if options.projector_type == 2 && options.reconstruction_method == 1
                lor_orth = (lor_orth(index));
            end
            if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
                if partitions > 1
                    for ff = 1 : partitions
                        temp = SinDelayed{ff};
                        temp = temp(index);
                        SinDelayed{ff} = temp;
                    end
                    clear temp
                else
                    SinDelayed = SinDelayed(index);
                end
            end
            clear lor
        else
            if ~exist('LL','var')
                LL = form_detector_pairs_raw(rings, det_per_ring);
                LL = LL(discard,:);
            end
            Sino = Sino(discard);
            lor_a = lor(discard);
            if options.normalization_correction && options.corrections_during_reconstruction
                normalization = normalization(discard);
            end
            if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
                if partitions > 1
                    for ff = 1 : partitions
                        temp = SinDelayed{ff};
                        temp = temp(discard);
                        SinDelayed{ff} = temp;
                    end
                    clear temp
                else
                    SinDelayed = SinDelayed(discard);
                end
            end
            pituus = uint32([0;length(Sino)]);
            if options.projector_type == 2 && options.reconstruction_method == 1
                lor_orth = (lor_orth(discard));
            end
            clear lor
        end
        summa = zeros(subsets, 1, 'uint64');
        
        for kk = 1 : subsets
            apu = LL(pituus(kk) + 1 : pituus(kk + 1),:) - 1;
            apu2 = idivide(apu, uint16(det_per_ring));
            idx = apu2(:,1) == apu2(:,2);
            apu2 = apu(idx,:);
            ind = mod(apu2, uint16(det_per_ring)) + 1;
            yt = y(ind);
            y_i = yt(:,1) > yt(:,2);
            apu2(y_i,:) = fliplr(apu2(y_i,:));
            apu(idx,:) = apu2;
            LL(pituus(kk) + 1 : pituus(kk + 1),:) = apu + 1;
            if options.projector_type == 2 && options.reconstruction_method == 1
                summa(kk) = sum(int64(lor_orth(pituus(kk)+1:pituus(kk+1))));
            else
                summa(kk) = sum(int64(lor_a(pituus(kk)+1:pituus(kk+1))));
            end
        end
        
        clear apu apu2 idx ind yt y_i index discard
        
        LL = LL';
        LL = LL(:);
    elseif use_raw_data == false && ~options.precompute_lor && options.reconstruction_method > 1
        
        if subsets > 1
            Sino = Sino(index);
            if options.normalization_correction && options.corrections_during_reconstruction
                normalization = normalization(index);
            end
            if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
                if partitions > 1
                    for ff = 1 : partitions
                        temp = SinDelayed{ff};
                        if NSinos ~= TotSinos
                            temp = temp(:,:,1:NSinos);
                        end
                        temp = temp(index);
                        SinDelayed{ff} = temp;
                    end
                    clear temp
                else
                    SinDelayed = SinDelayed(index);
                end
            end
            if partitions > 1
                for ff = 1 : partitions
                    temp = SinM{ff};
                    if NSinos ~= TotSinos
                        temp = temp(:,:,1:NSinos);
                    end
                    temp = temp(index);
                    SinM{ff} = temp;
                end
                clear temp
            end
            clear lor
        end
        [~, I] = sort(y, 2);
        sy = size(y);
        I = sub2ind(sy, repmat((1:sy(1)).', 1, sy(2)), I);
        
        xy_index = uint32(I(:,1));
        xy_index2 = repmat(uint32(1:size_x)', NSinos - Nz, 1);
        xy_index = [repmat(xy_index, Nz, 1); xy_index2];
        [~, I] = sort(z_det, 2);
        sy = size(z_det);
        I = sub2ind(sy, repmat((1:sy(1)).', 1, sy(2)), I);
        
        z_index = uint16(I(:,1));
        if verLessThan('matlab','8.5')
            z_index = repeat_elem(z_index, size_x);
        else
            z_index = repelem(z_index, size_x);
        end
        if subsets > 1
            z_index = z_index(index);
        end
        apu = z_index > NSinos;
        z_index = z_index - 1;
        
        if subsets > 1
            xy_index = xy_index(index);
        end
        xy_index(apu) = xy_index(apu) + uint32(size_x);
        xy_index = xy_index - 1;
        
        summa = zeros(subsets, 1, 'uint64');
        
        if subsets > 1
        else
            pituus = uint32([0;length(Sino)]);
        end
        
        
        clear discard I yt xt xy_index2 index apu
    elseif use_raw_data && ~options.precompute_lor
        
        if subsets > 1
            if ~exist('LL','var')
                LL = form_detector_pairs_raw(rings, det_per_ring);
            end
            LL = LL(index,:);
            Sino = Sino(index);
            if options.normalization_correction && options.corrections_during_reconstruction
                normalization = normalization(index);
            end
            if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
                if partitions > 1
                    for ff = 1 : partitions
                        temp = SinDelayed{ff};
                        temp = temp(index);
                        SinDelayed{ff} = temp;
                    end
                    clear temp
                else
                    SinDelayed = SinDelayed(index);
                end
            end
            clear lor
        end
        if partitions > 1
            for ff = 1 : partitions
                temp = SinM{ff};
                if subsets > 1
                    temp = temp(index);
                end
                SinM{ff} = temp;
            end
            clear temp
        end
        summa = zeros(subsets, 1, 'uint64');
        
        for kk = 1 : subsets
            apu = LL(pituus(kk) + 1 : pituus(kk + 1),:) - 1;
            apu2 = idivide(apu, uint16(det_per_ring));
            idx = apu2(:,1) == apu2(:,2);
            apu2 = apu(idx,:);
            ind = mod(apu2, uint16(det_per_ring)) + 1;
            yt = y(ind);
            y_i = yt(:,1) > yt(:,2);
            apu2(y_i,:) = fliplr(apu2(y_i,:));
            apu(idx,:) = apu2;
            LL(pituus(kk) + 1 : pituus(kk + 1),:) = apu + 1;
        end
        
        clear apu apu2 idx ind yt y_i index discard
        
        LL = LL';
        LL = LL(:);
    end
    
    % Pixels
    etaisyys_x = (R - FOVax)/2;
    etaisyys_y = (R - FOVay)/2;
    if options.reconstruction_method == 2 || options.reconstruction_method == 3 || options.reconstruction_method == 5
        zz = linspace(single(0), single(axial_fow), Nz + 1);
        xx = single(linspace(etaisyys_x, R - etaisyys_x, Nx + 1));
        yy = single(linspace(etaisyys_y, R - etaisyys_y, Ny + 1));
    else
        zz = linspace(double(0), double(axial_fow), Nz + 1);
        xx = double(linspace(etaisyys_x, R - etaisyys_x, Nx + 1));
        yy = double(linspace(etaisyys_y, R - etaisyys_y, Ny + 1));
    end
    zz = zz(2*block1 + 1 : 2*blocks + 2);
    
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
    
    %     if options.reconstruction_method == 2 || options.reconstruction_method == 4
    %         iij=single(0:Nx);
    %         jji=single(0:Ny);
    %         kkj=single(0:Nz);
    %     else
    %         iij=double(0:Nx);
    %         jji=double(0:Ny);
    %         kkj=double(0:Nz);
    %     end
    
    N = (Nx)*(Ny)*(Nz);
    det_per_ring = uint32(det_per_ring);
    
    % How much memory is preallocated
    if ~use_raw_data
        ind_size = uint32(NSinos / subsets * (det_per_ring) * Nx * (Ny));
    else
        ind_size = uint32((det_per_ring)^2 / subsets * Nx * (Ny));
    end
    
    
    zmax = max(max(z_det));
    if zmax==0
        if options.reconstruction_method == 2 || options.reconstruction_method == 3 || options.reconstruction_method == 5
            zmax = single(1);
        else
            zmax = double(1);
        end
    end
    
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
    %%
    
    if (options.MRP || options.quad || options.TV ||options. FMH || options.L || options.weighted_mean || options.APLS || options.BSREM ...
            || options.ramla || options.MBSREM || options.mramla || options.rosem || options.drama || options.ROSEM_MAP || options.ecosem ...
            || options.cosem || options.acosem || options.AD || any(options.COSEM_MAP) || (options.NLM && options.NLM_use_anatomical))
        
        if options.TV && options.MAP
            if options.TV_use_anatomical
                apu = load(options.TV_reference_image);
                variables = fields(apu);
                alkuarvo = double(apu.(variables{1}));
                if size(alkuarvo,2) == 1
                    koko_apu = sqrt(length(alkuarvo)/double(Nz));
                    if floor(koko_apu) ~= koko_apu
                        error("Reference image has to be square")
                    else
                        alkuarvo = reshape(alkuarvo, koko_apu,koko_apu,Nz);
                        if koko_apu ~= Nx
                            alkuarvo = imresize(alkuarvo, Nx, Ny, Nz);
                        end
                    end
                else
                    if size(alkuarvo,2) ~= Nx
                        alkuarvo = imresize(alkuarvo, Nx, Ny, Nz);
                    end
                end
                alkuarvo = alkuarvo - min(min(min(alkuarvo)));
                alkuarvo = alkuarvo/max(max(max(alkuarvo)));
                if options.TVtype == 1
                    S = assembleS(alkuarvo,options.T,Ny,Nx,Nz);
                    if options.reconstruction_method == 2
                        TVdata.s1 = single(S(1:3:end,1));
                        TVdata.s2 = single(S(1:3:end,2));
                        TVdata.s3 = single(S(1:3:end,3));
                        TVdata.s4 = single(S(2:3:end,1));
                        TVdata.s5 = single(S(2:3:end,2));
                        TVdata.s6 = single(S(2:3:end,3));
                        TVdata.s7 = single(S(3:3:end,1));
                        TVdata.s8 = single(S(3:3:end,2));
                        TVdata.s9 = single(S(3:3:end,3));
                    else
                        TVdata.s1 = S(1:3:end,1);
                        TVdata.s2 = S(1:3:end,2);
                        TVdata.s3 = S(1:3:end,3);
                        TVdata.s4 = S(2:3:end,1);
                        TVdata.s5 = S(2:3:end,2);
                        TVdata.s6 = S(2:3:end,3);
                        TVdata.s7 = S(3:3:end,1);
                        TVdata.s8 = S(3:3:end,2);
                        TVdata.s9 = S(3:3:end,3);
                    end
                end
                if options.reconstruction_method == 2
                    TVdata.reference_image = single(alkuarvo);
                    TVdata.T = single(options.T);
                    TVdata.C = single(options.C);
                else
                    TVdata.reference_image = alkuarvo;
                    TVdata.T = options.T;
                    TVdata.C = options.C;
                end
                clear apu variables alkuarvo S
            end
            if options.reconstruction_method == 2
                options.tau = single(options.tau);
                TVdata.beta = single(options.TVsmoothing);
                options.TVdata = TVdata;
                clear TVdata;
            else
                TVdata.beta = options.TVsmoothing;
            end
        end
        
        if options.TV && options.MAP && options.reconstruction_method == 2
            options.alphaTGV = single(options.alphaTGV);
            options.betaTGV = single(options.betaTGV);
            options.NiterTGV = uint32(options.NiterTGV);
        end
        
        if options.APLS && options.MAP
            apu = load(options.im_vectors.APLS_reference_image);
            variables = fields(apu);
            alkuarvo = double(apu.(variables{1}));
            if size(alkuarvo,2) == 1
                koko_apu = sqrt(length(alkuarvo)/double(Nz));
                if floor(koko_apu) ~= koko_apu
                    error("Reference image has to be square")
                else
                    alkuarvo = reshape(alkuarvo, koko_apu,koko_apu,Nz);
                    if koko_apu ~= Nx
                        alkuarvo = imresize(alkuarvo, Nx, Ny, Nz);
                    end
                end
            else
                if size(alkuarvo,2) ~= Nx
                    alkuarvo = imresize(alkuarvo, Nx, Ny, Nz);
                end
            end
            alkuarvo = alkuarvo - min(min(min(alkuarvo)));
            alkuarvo = alkuarvo/max(max(max(alkuarvo)));
            if options.reconstruction_method == 2
                options.im_vectors.APLS_ref_image = single(alkuarvo);
                options.eta = single(options.eta);
                options.APLSsmoothing = single(options.APLSsmoothing);
            else
                options.im_vectors.APLS_ref_image = (alkuarvo);
            end
            clear alkuarvo apu variables
        end
        
        if ((options.mramla || options.MBSREM || options.rbi || options.RBI_MAP) && options.MBSREM_prepass || options.ecosem || options.cosem ...
                || options.acosem || any(options.COSEM_MAP))  && options.reconstruction_method == 1
            
            if options.acosem
                C_aco = zeros(double(N), subsets);
            end
            if options.cosem || options.ecosem
                C_co = zeros(double(N), subsets);
            end
            if any(options.COSEM_MAP)
                C_osl = zeros(double(N), subsets);
            end
            if options.mramla || options.MBSREM
                Amin = zeros(length(Sino));
            end
            
            if ~use_raw_data
                if isempty(pseudot)
                    pseudot = uint32(0);
                end
            end
            
            D = zeros(N,1);
            
            if verbose
                disp('Prepass phase for MRAMLA, COSEM, ACOSEM and ECOSEM started')
            end
            for osa_iter = 1 : subsets
                if options.precompute_lor == false
                    if use_raw_data == false
                        if options.projector_type == 1 || options.projector_type == 0
                            if exist('projector_mex','file') == 3
                                [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                    zmax, vaimennus, normalization, 0, pituus(osa_iter), attenuation_correction, normalization_correction, ...
                                    randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, uint16(0), pseudot, det_per_ring, options.verbose, ...
                                    use_raw_data, uint32(2), ind_size, block1, blocks, index{osa_iter}, uint32(options.projector_type));
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
                                    NSinos, NSlices, vaimennus, index{osa_iter}, pituus(osa_iter), attenuation_correction);
                                alkiot = cell2mat(alkiot);
                                indices = indices(discard);
                                indices = cell2mat(indices) - 1;
                            end
                        elseif options.projector_type == 2
                            [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                zmax, vaimennus, normalization, 0, pituus(osa_iter), attenuation_correction, normalization_correction, ...
                                randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, uint16(0), pseudot, det_per_ring, options.verbose, ...
                                use_raw_data, uint32(2), ind_size, block1, blocks, index{osa_iter}, uint32(options.projector_type), ...
                                options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z);
                        else
                            error('Unsupported projector type')
                        end
                    else
                        L = LL(index{osa_iter},:);
                        L = L';
                        L = L(:);
                        if options.projector_type == 1
                            [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                zmax, vaimennus, normalization, 0, uint32(0), attenuation_correction, normalization_correction, ...
                                randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, L, pseudot, det_per_ring, options.verbose, ...
                                use_raw_data, uint32(2), ind_size, block1, blocks, uint32(0), uint32(options.projector_type));
                        elseif options.projector_type == 2
                            [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                zmax, vaimennus, normalization, 0, uint32(0), attenuation_correction, normalization_correction, ...
                                randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, L, pseudot, det_per_ring, options.verbose, ...
                                use_raw_data, uint32(2), ind_size, block1, blocks, uint32(0), uint32(options.projector_type), options.tube_width_xy, ...
                                x_center, y_center, z_center, options.tube_width_z);
                        else
                            error('Unsupported projector type')
                        end
                    end
                    lor = reshape(lor,[],2);
                    if verLessThan('matlab','8.5')
                        lor = repeat_elem(uint32((lor(:,1))),lor(:,2));
                    else
                        lor = repelem(uint32((lor(:,1))),lor(:,2));
                    end
                    uu = double(Sino(index{osa_iter}));
                    
                    A_length = length(uu);
                    indices=indices + 1;
                    if verbose
                        tStart = tic;
                    end
                    if options.use_fsparse == false
                        A = sparse(double(lor),double(indices),double(alkiot), A_length, double(N));
                    else
                        A = fsparse(lor,indices,double(alkiot),[A_length double(N) length(alkiot)]);
                    end
                    clear indices alkiot lor
                    if verbose
                        tElapsed = toc(tStart);
                        disp(['Sparse matrix formation took ' num2str(tElapsed) ' seconds'])
                    end
                else
                    if options.projector_type == 2
                        lor2 = [0; cumsum(uint64(lor_orth(pituus(osa_iter)+1:pituus(osa_iter + 1))))];
                        if ~use_raw_data
                            [A, ~] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                zmax, vaimennus, normalization, 0, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, ...
                                randoms_correction, lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1)), xy_index(pituus(osa_iter)+1:pituus(osa_iter + 1)), ...
                                z_index(pituus(osa_iter)+1:pituus(osa_iter + 1)), NSinos, uint16(0), pseudot, det_per_ring, options.verbose, ...
                                use_raw_data, uint32(0), lor2, summa(osa_iter), options.attenuation_phase, uint32(options.projector_type), ...
                                options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, int32(options.accuracy_factor));
                        else
                            [A, ~] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                zmax, vaimennus, normalization, 0, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, ...
                                randoms_correction, lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1)), uint32(0), uint32(0), NSinos, ...
                                LL(pituus(osa_iter) * 2 + 1 : pituus(osa_iter + 1) * 2), pseudot, det_per_ring, options.verbose, use_raw_data, ...
                                uint32(0), lor2, summa(osa_iter), options.attenuation_phase, uint32(options.projector_type), options.tube_width_xy, ...
                                x_center, y_center, z_center, options.tube_width_z);
                        end
                    else
                        lor2 = [0; cumsum(uint64(lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1))))];
                        if ~use_raw_data
                            [A, ~] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                zmax, vaimennus, normalization, 0, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, ...
                                randoms_correction, lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1)), xy_index(pituus(osa_iter)+1:pituus(osa_iter + 1)), ...
                                z_index(pituus(osa_iter)+1:pituus(osa_iter + 1)), NSinos, uint16(0), pseudot, det_per_ring, options.verbose, ...
                                use_raw_data, uint32(0), lor2, summa(osa_iter), options.attenuation_phase, uint32(options.projector_type));
                        else
                            [A, ~] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                zmax, vaimennus, normalization, 0, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, ...
                                randoms_correction, lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1)), uint32(0), uint32(0), NSinos, ...
                                LL(pituus(osa_iter) * 2 + 1 : pituus(osa_iter + 1) * 2), pseudot, det_per_ring, options.verbose, ...
                                use_raw_data, uint32(0), lor2, summa(osa_iter), options.attenuation_phase, uint32(options.projector_type));
                        end
                    end
                    uu = double(Sino(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                    clear lor2
                end
                if is_transposed
                    D = D + A * ones(size(A,2),1,'double');
                else
                    D = D + full(sum(A,1))';
                end
                if options.ecosem || options.cosem || options.acosem || any(options.COSEM_MAP)
                    if options.precompute_lor == false
                        uu = double(Sino(index{osa_iter}));
                    else
                        uu = double(Sino(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                    end
                end
                if options.cosem || options.ecosem
                    if osa_iter > 1
                        if is_transposed
                            C_co(:,osa_iter) = full(sum(spdiags(options.x0(:),0,size(A,1),size(A,1)) * A * spdiags(uu ./ (A' * options.x0(:) + epps),...
                                0,size(A,2),size(A,2)),2));
                        else
                            C_co(:,osa_iter) = full(sum(spdiags(options.x0(:),0,size(A,2),size(A,2)) * A' * spdiags(uu ./ (A * options.x0(:) + epps),...
                                0,size(A,1),size(A,1)),2));
                        end
                    end
                end
                if options.acosem
                    if osa_iter > 1
                        if is_transposed
                            C_aco(:,osa_iter) = full(sum(spdiags(power(options.x0(:),1/options.h),0,size(A,1),size(A,1)) * A * ...
                                spdiags(uu ./ (A' * options.x0(:) + epps),0,size(A,2),size(A,2)),2));
                        else
                            C_aco(:,osa_iter) = full(sum(spdiags(power(options.x0(:),1/options.h),0,size(A,2),size(A,2)) * A' * ...
                                spdiags(uu ./ (A * options.x0(:) + epps),0,size(A,1),size(A,1)),2));
                        end
                    end
                end
                if any(options.COSEM_MAP)
                    if options.COSEM_MAP == 2
                        if osa_iter > 1
                            if is_transposed
                                C_osl(:,osa_iter) = full(sum(spdiags(options.x0(:),0,size(A,1),size(A,1)) * A * spdiags(uu ./ (A' * options.x0(:) + epps),...
                                    0,size(A,2),size(A,2)),2));
                            else
                                C_osl(:,osa_iter) = full(sum(spdiags(options.x0(:),0,size(A,2),size(A,2)) * A' * spdiags(uu ./ (A * options.x0(:) + epps),...
                                    0,size(A,1),size(A,1)),2));
                            end
                        end
                    else
                        if osa_iter > 1
                            if is_transposed
                                C_osl(:,osa_iter) = full(sum(spdiags(power(options.x0(:),1/options.h),0,size(A,1),size(A,1)) * A * ...
                                    spdiags(uu ./ (A' * options.x0(:) + epps),0,size(A,2),size(A,2)),2));
                            else
                                C_osl(:,osa_iter) = full(sum(spdiags(power(options.x0(:),1/options.h),0,size(A,2),size(A,2)) * A' * ...
                                    spdiags(uu ./ (A * options.x0(:) + epps),0,size(A,1),size(A,1)),2));
                            end
                        end
                    end
                end
                if options.MBSREM_prepass && options.U == 0 && (options.MBSREM || options.mramla)
                    %%%% This particular piece of code was taken from:
                    %%%% https://se.mathworks.com/matlabcentral/answers/35309-max-min-of-sparse-matrices
                    if is_transposed
                        [~,m] = size(A);
                        rowMin = nan(m, 1);
                        [~,I,S] = find(A);
                    else
                        [m,~] = size(A);
                        rowMin = nan(m, 1);
                        [I,~,S] = find(A);
                    end
                    I = I(S>1e-10);
                    S = S(S>1e-10);
                    [I,K] = sort(I);
                    S = S(K);
                    markers = [find([1; diff(I)]); numel(I)+1];
                    iRows = I(markers(1:end-1));
                    for i = 1:numel(iRows)
                        s = S(markers(i):(markers(i+1)-1));
                        rowMin(iRows(i)) = min(s);
                    end
                    rowMin(isnan(rowMin)) = epps;
                    if options.precompute_lor == false
                        Amin(index{osa_iter}) = rowMin;
                    else
                        Amin(pituus(osa_iter)+1:pituus(osa_iter + 1)) = rowMin;
                    end
                    clear I K S markers rowMin s iRows
                end
                clear A
            end
            %             D = sum(pj,2);
            if verbose
                disp('Prepass phase for COSEM, ACOSEM and ECOSEM completed')
            end
        end
        if (options.BSREM || options.ramla) && length(options.lambda0) == 1
            lam = zeros(Niter,1);
            lam(1) = options.lambda0;
%             orig_lam = lam;
            %             if lam(1) > 1/max(max(pj))
            %                 lam(1) = min(min(pj));
            %             end
            for i=2:Niter
                %                 lam(i) = 0.5*lam(i-1);
                lam(i) = lam(1)/i;
                %                 lam(i) = lam(1)/1.01;
            end
            if options.reconstruction_method == 2
                options.lam = single(lam);
            else
                options.lambda0 = lam;
            end
        elseif (options.BSREM || options.ramla) && options.reconstruction_method == 2
            options.lam = single(options.lam);
        end
        if (options.MBSREM || options.mramla) && length(options.lambda0_mbsrem) == 1
            lam_mbsrem = zeros(Niter,1);
            lam_mbsrem(1) = options.lambda0_mbsrem;
            for i=2:Niter
                lam_mbsrem(i) = lam_mbsrem(1)/sqrt(i);
            end
            if options.reconstruction_method == 2
                options.lam_mbsrem = single(lam_mbsrem);
            else
                options.lam_mbsrem = lam_mbsrem;
            end
        elseif (options.MBSREM || options.mramla) && options.reconstruction_method == 2
            options.lam_mbsrem = single(options.lam_mbsrem);
        end
        if (options.ROSEM_MAP || options.rosem) && length(options.lambda0_rosem) == 1
            lam_rosem = zeros(Niter,1);
            lam_rosem(1) = options.lambda0_rosem;
            for i=2:Niter
                lam_rosem(i) = lam_rosem(1)/i;
            end
            if options.reconstruction_method == 2
                options.lam_rosem = single(lam_rosem);
            else
                options.lam_rosem = lam_rosem;
            end
        elseif (options.MBSREM || options.mramla) && options.reconstruction_method == 2
            options.lambda0_rosem = single(options.lambda0_rosem);
        end
        if options.drama
            lam_drama = zeros(Niter,subsets);
            lam_drama(1,1) = options.beta_drama/(options.alpha_drama*options.beta0_drama);
            r = 1;
            for i=1:Niter
                for j = 1 : subsets
                    lam_drama(i,j) = options.beta_drama/(options.alpha_drama*options.beta0_drama + r);
                    r = r + 1;
                end
            end
            if options.reconstruction_method == 2
                options.lam_drama = single(lam_drama);
            else
                options.lam_drama = lam_drama;
            end
        end
        if (options.MBSREM || options.mramla) && options.reconstruction_method == 1
            pj3 = D/subsets;
        end
        if (options.quad || options.L || options.FMH || options.weighted_mean || options.MRP || (options.TV && options.TVtype == 3)) && options.MAP
            if options.quad || options.L || options.FMH || options.weighted_mean || (options.TV && options.TVtype == 3)
                distX = FOVax/double(Nx);
                distY = FOVay/double(Ny);
                distZ = (double(axial_fov)/double(Nz));
                if isempty(options.weights)
                    options.weights = zeros(((Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1)),1);
                    edist = zeros((Ndx*2+1),1);
                    cc = zeros((Ndy*2+1)*(Ndx*2+1),1);
                    lt = 0;
                    for jj = Ndz : -1 : -Ndz
                        lt = lt + 1;
                        ll = 0;
                        for kk = Ndy : -1 : -Ndy
                            ll = ll + 1;
                            if Ndz == 0 || Nz == 1
                                if verLessThan('matlab','8.5')
                                    apu = [((Ndx:-1:-Ndx) * distX)', (repeat_elem(kk,Ndy*2+1) * distY)'];
                                else
                                    apu = [((Ndx:-1:-Ndx) * distX)', (repelem(kk,Ndy*2+1) * distY)'];
                                end
                            else
                                if Ndz ~= Ndx
                                    if verLessThan('matlab','8.5')
                                        apu = [((Ndx:-1:-Ndx) * distX)', (repeat_elem(kk,Ndy*2+1) * distY)', [zeros(Ndx-Ndz,1),(repeat_elem(jj,Ndz*2+1) * distZ),...
                                            zeros(Ndx-Ndz,1)]'];
                                    else
                                        apu = [((Ndx:-1:-Ndx) * distX)', (repelem(kk,Ndy*2+1) * distY)', [zeros(Ndx-Ndz,1),(repelem(jj,Ndz*2+1) * distZ),...
                                            zeros(Ndx-Ndz,1)]'];
                                    end
                                else
                                    if verLessThan('matlab','8.5')
                                        apu = [((Ndx:-1:-Ndx) * distX)', (repeat_elem(kk,Ndy*2+1) * distY)', (repeat_elem(jj,Ndz*2+1) * distZ)'];
                                    else
                                        apu = [((Ndx:-1:-Ndx) * distX)', (repelem(kk,Ndy*2+1) * distY)', (repelem(jj,Ndz*2+1) * distZ)'];
                                    end
                                end
                            end
                            for ii = 1 : length(apu)
                                edist(ii) = sqrt(apu(ii,:)*apu(ii,:)');
                            end
                            cc((Ndy*2+1)*(ll-1)+1:(Ndy*2+1)*ll) = edist;
                        end
                        options.weights((Ndx*2+1) * (Ndy*2+1) * (lt - 1) + 1: (Ndx*2+1) * (Ndy*2+1) * lt) = cc;
                    end
                    options.weights = 1./options.weights;
                end
            end
            %             pz_pad = padding(reshape(options.x0(:),Nx,Ny,Nz),[Ndx Ndy Ndz]);
            s = [Nx + Ndx*2 Ny + Ndy*2 Nz + Ndz*2];
            N_pad = min(3, Ndx + Ndy + Ndz);
            [c1{1:N_pad}]=ndgrid(1:(Ndx*2+1));
            c2(1:N_pad)={Ndy+1};
            if Ndz > Ndx && Ndz > 1
                c1{1} = cat(3, c1{1}, zeros(size(c1{1},1), size(c1{1},2), Ndz));
                c1{2} = cat(3, c1{2}, zeros(size(c1{2},1), size(c1{2},2), Ndz));
                c1{3} = cat(3, c1{3}, zeros(size(c1{3},1), size(c1{3},2), Ndz));
                %                 apu2 = c1{2};
                %                 apu3 = c1{3};
                for kk = Ndz - 1 : - 1 : 0
                    %                     apu(:,:,end+1) = apu(:,:,end);
                    %                     apu2(:,:,end+1) = apu2(:,:,end);
                    %                     apu3(:,:,end+1) = apu3(:,:,end) + 1;
                    c1{1}(:,:,end-kk) = c1{1}(:,:,end - kk - 1);
                    c1{2}(:,:,end-kk) = c1{2}(:,:,end - kk - 1);
                    c1{3}(:,:,end-kk) = c1{3}(:,:,end - kk - 1) + 1;
                end
                %                 c1{1} = apu;
                %                 c1{2} = apu2;
                %                 c1{3} = apu3;
                c2(end) = {Ndz+1};
            elseif Ndz < Ndx && Ndz > 1
                %                 apu = c1{1};
                %                 apu2 = c1{2};
                %                 apu3 = c1{3};
                c1{1}(:,:,end-2*(Ndx-Ndz) + 1) = [];
                c1{2}(:,:,end-2*(Ndx-Ndz) + 1) = [];
                c1{3}(:,:,end-2*(Ndx-Ndz) + 1) = [];
                %                 c1{1} = apu;
                %                 c1{2} = apu2;
                %                 c1{3} = apu3;
                c2(end) = {Ndz+1};
            end
            offsets=sub2ind(s,c1{:}) - sub2ind(s,c2{:});
            if Nz == 1
                tr_ind = sub2ind([Nx+Ndx*2 Ny + Ndy*2],mod((1:N)'-1,Nx)+(Ndx + 1),mod(floor(((1:double(N))'-1)/double(Nx)),double(Ny))+(Ndy + 1));
            else
                tr_ind = sub2ind([Nx+Ndx*2 Ny+Ndy*2 Nz+Ndz*2], mod((1:N)' - 1, Nx) + (Ndx + 1), mod(floor(((1:double(N))' - 1)/double(Nx)), ...
                    double(Ny)) + (Ndy + 1), floor(((1:double(N))' - 1)/double(Nx * Ny)) + (Ndz + 1));
            end
            tr_offsets = uint32(bsxfun(@plus,tr_ind,offsets(:)'));
            if options.reconstruction_method == 2
                options.tr_offsets = tr_offsets - 1;
                options.Ndx = uint32(Ndx);
                options.Ndy = uint32(Ndy);
                options.Ndz = uint32(Ndz);
                clear tr_offsets
            end
            %             if (options.OSL_OSEM || options.OSL_MLEM) && options.quad
            %                 pz_pad_osl = pz_pad;
            %                 if options.reconstruction_method == 2
            %                     options.pz_pad_osl = single(pz_pad_osl);
            %                     clear pz_pad_osl
            %                 end
            %             end
            if options.quad || options.TV && options.TVtype == 3
                if empty_weight
                    options.weights_quad = options.weights/sum(options.weights(~isinf(options.weights)));
                    options.weights_quad = [options.weights_quad(1:floor(length(options.weights_quad) / 2)); ...
                        options.weights_quad(ceil(length(options.weights_quad)/2) + 1 : end)];
                else
                    options.weights_quad = options.weights;
                end
                if options.reconstruction_method == 2
                    options.weights_quad = single(options.weights_quad);
                    clear weights_quad
                end
            end
            if options.L
                if isempty(options.a_L)
                    options.a_L = lfilter_weights(Ndx, Ndy, Ndz, dx, dy, dz, options.oneD_weights);
                end
                if options.reconstruction_method == 2
                    options.a_L = single(options.a_L);
                    clear a_L
                end
                clear dd
            end
            %             if (options.OSL_OSEM || options.OSL_MLEM )&& options.L
            %                 pz_pad_L_osl  = pz_pad;
            %                 if options.reconstruction_method == 2
            %                     options.pz_pad_L_osl = single(pz_pad_L_osl);
            %                     clear pz_pad_L_osl
            %                 end
            %             end
            if options.FMH
                if isempty(options.fmh_weights)
                    kerroin = sqrt(2)*distX;
                    if Nz == 1 || Ndz == 0
                        options.fmh_weights = zeros(Ndx*2+1, 4);
                        for jjj = 1 : 4
                            lll = lll + 1;
                            apu = zeros(Ndx*2+1,1);
                            hhh = 0;
                            if jjj == 1 || jjj == 3
                                for iii = Ndx : -1 : -Ndx
                                    hhh = hhh + 1;
                                    if iii == 0
                                        apu(hhh) = options.fmh_center_weight;
                                    else
                                        apu(hhh) = kerroin/sqrt((distX*iii)^2+(distY*iii)^2);
                                    end
                                end
                            elseif jjj == 2
                                for iii = Ndx : -1 : -Ndx
                                    hhh = hhh + 1;
                                    if iii == 0
                                        apu(hhh) = options.fmh_center_weight;
                                    else
                                        apu(hhh) = kerroin/abs(distX*iii);
                                    end
                                end
                            elseif jjj == 4
                                for iii = Ndx : -1 : -Ndx
                                    hhh = hhh + 1;
                                    if iii == 0
                                        apu(hhh) = options.fmh_center_weight;
                                    else
                                        apu(hhh) = kerroin/abs(distY*iii);
                                    end
                                end
                            end
                            options.fmh_weights(:, jjj) = apu;
                        end
                    else
                        options.fmh_weights = zeros(max([Ndx*2+1,Ndz*2+1]), 13);
                        lll = 0;
                        for kkk = 1 : -1 : 0
                            for jjj = 1 : 9
                                lll = lll + 1;
                                if kkk == 1
                                    apu = zeros(Ndz*2+1,1);
                                    hhh = 0;
                                    if jjj == 1 || jjj == 3 || jjj == 7 || jjj == 9
                                        for iii = Ndz : -1 : -Ndz
                                            hhh = hhh + 1;
                                            if iii == 0
                                                apu(hhh) = options.fmh_center_weight;
                                            else
                                                apu(hhh) = kerroin/sqrt(sqrt((distZ*iii)^2+(distX*iii)^2)^2+(distY*iii)^2);
                                            end
                                        end
                                    elseif jjj == 2 || jjj == 8
                                        for iii = Ndz : -1 : -Ndz
                                            hhh = hhh + 1;
                                            if iii == 0
                                                apu(hhh) = options.fmh_center_weight;
                                            else
                                                apu(hhh) = kerroin/sqrt((distZ*iii)^2+(distX*iii)^2);
                                            end
                                        end
                                    elseif jjj == 4 || jjj == 6
                                        for iii = Ndz : -1 : -Ndz
                                            hhh = hhh + 1;
                                            if iii == 0
                                                apu(hhh) = options.fmh_center_weight;
                                            else
                                                apu(hhh) = kerroin/sqrt((distZ*iii)^2+(distY*iii)^2);
                                            end
                                        end
                                    elseif jjj == 5
                                        for iii = Ndz : -1 : -Ndz
                                            hhh = hhh + 1;
                                            if iii == 0
                                                apu(hhh) = options.fmh_center_weight;
                                            else
                                                apu(hhh) = kerroin/abs(distZ*iii);
                                            end
                                        end
                                    end
                                    options.fmh_weights(:, lll) = apu;
                                else
                                    apu = zeros(Ndx*2+1,1);
                                    hhh = 0;
                                    if jjj == 1 || jjj == 3
                                        for iii = Ndx : -1 : -Ndx
                                            hhh = hhh + 1;
                                            if iii == 0
                                                apu(hhh) = options.fmh_center_weight;
                                            else
                                                apu(hhh) = kerroin/sqrt((distX*iii)^2+(distY*iii)^2);
                                            end
                                        end
                                    elseif jjj == 2
                                        for iii = Ndx : -1 : -Ndx
                                            hhh = hhh + 1;
                                            if iii == 0
                                                apu(hhh) = options.fmh_center_weight;
                                            else
                                                apu(hhh) = kerroin/abs(distX*iii);
                                            end
                                        end
                                    elseif jjj == 4
                                        for iii = Ndx : -1 : -Ndx
                                            hhh = hhh + 1;
                                            if iii == 0
                                                apu(hhh) = options.fmh_center_weight;
                                            else
                                                apu(hhh) = kerroin/abs(distY*iii);
                                            end
                                        end
                                    else
                                        break
                                    end
                                    options.fmh_weights(:, lll) = apu;
                                end
                            end
                        end
                    end
                    options.fmh_weights = options.fmh_weights./sum(options.fmh_weights,1);
                end
                %                 if options.OSL_OSEM || options.OSL_MLEM
                %                     pz_pad_fmh = pz_pad;
                %                 end
                
                if options.reconstruction_method == 2
                    options.fmh_weights = single(options.fmh_weights);
                    clear fmh_weights pz_pad_fmh
                end
            end
            if (options.FMH || options.quad) && options.reconstruction_method == 2
                options.weights = single(options.weights);
                options.inffi = uint32(find(isinf(options.weights)) - 1);
            end
            if options.MRP
                medx = Ndx*2 + 1;
                medy = Ndy*2 + 1;
                medz = Ndz*2 + 1;
            end
            if options.weighted_mean
                if isempty(options.weighted_weights)
                    kerroin = sqrt(2)*distX;
                    options.weighted_weights = kerroin.*options.weights;
                    options.weighted_weights(isinf(options.weighted_weights)) = options.weighted_center_weight;
                    %                 options.weighted_weights = options.weighted_weights/sum(options.weighted_weights);
                end
                %                 if options.OSL_OSEM || options.OSL_MLEM
                %                     pz_pad_weighted = pz_pad;
                %                 end
                options.w_sum = sum(options.weighted_weights);
                if options.reconstruction_method == 2
                    options.weighted_weights = single(options.weighted_weights);
                    %                     if options.OSL_OSEM || options.OSL_MLEM
                    %                         options.pz_pad_weighted = single(pz_pad_weighted);
                    %                     end
                    clear weighted_weights pz_pad_weighted
                end
            end
            clear tr_ind offsets c1 c2 apu apu2 apu3 N_pad cc pz_pad
            if verbose
                disp('Prepass phase for MRP, quadratic prior, L-filter, FMH and weighted mean completed')
            end
        end
        if options.AD && options.MAP
            if options.reconstruction_method == 2
                options.NiterAD = uint32(options.NiterAD);
                options.KAD = single(options.KAD);
                options.TimeStepAD = single(options.TimeStepAD);
                options.FluxType = uint32(options.FluxType);
                options.DiffusionType = uint32(options.DiffusionType);
            else
                if options.FluxType == 1
                    FluxType = 'exponential';
                elseif options.FluxType == 2
                    FluxType = 'quadratic';
                end
            end
        end
        if options.MBSREM || options.mramla
            epsilon_mramla = MBSREM_epsilon(Sino, options.epps);
            if options.reconstruction_method == 2
                options.epsilon_mramla = single(epsilon_mramla);
                clear epsilon_mramla
            end
        end
        if (options.NLM && options.NLM_use_anatomical) && options.MAP
            apu = load(options.NLM_reference_image);
            variables = fields(apu);
            options.NLM_ref = double(apu.(variables{1}));
            options.NLM_ref = reshape(options.NLM_ref, Nx, Ny, Nz);
        end
    end
    
    %%
    
    
    disp('Starting image reconstruction')
    
    if options.reconstruction_method ~= 2 && options.reconstruction_method ~= 3
        
        
        for llo = 1 : partitions
            if partitions > 1
                if iscell(SinM)
                    Sino = SinM{llo};
                else
                    Sino = SinM;
                    clear SinM
                end
                
                Sino = Sino(:);
                
                if issparse(Sino)
                    Sino = (full(Sino));
                end
                
                if randoms_correction
                    SinD = SinDelayed{llo};
                    SinD = SinD(:);
                    if issparse(SinD)
                        SinD = (full(SinD));
                    end
                else
                    SinD = 0;
                end
            else
                if randoms_correction
                    SinD = SinDelayed;
                    SinD = SinD(:);
                    if issparse(SinD)
                        SinD = (full(SinD));
                    end
                else
                    SinD = 0;
                end
            end
            
            if options.reconstruction_method == 1
                
                if options.MBSREM || options.mramla
                    if options.U == 0 || isempty(options.U)
                        options.U = max(double(Sino)./Amin);
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
                        if options.precompute_lor == false
                            if use_raw_data == false
                                if options.projector_type == 1 || options.projector_type == 0
                                    if exist('projector_mex','file') == 3
                                        [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                            zmax, vaimennus, normalization, SinD, pituus(osa_iter), attenuation_correction, normalization_correction, ...
                                            randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, uint16(0), pseudot, det_per_ring, options.verbose, ...
                                            use_raw_data, uint32(2), ind_size, block1, blocks, index{osa_iter}, uint32(options.projector_type));
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
                                            NSinos, NSlices, vaimennus, index{osa_iter}, pituus(osa_iter), attenuation_correction);
                                        alkiot = cell2mat(alkiot);
                                        indices = indices(discard);
                                        indices = cell2mat(indices) - 1;
                                    end
                                elseif options.projector_type == 2
                                    [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                            zmax, vaimennus, normalization, SinD, pituus(osa_iter), attenuation_correction, normalization_correction, ...
                                            randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, uint16(0), pseudot, det_per_ring, options.verbose, ...
                                            use_raw_data, uint32(2), ind_size, block1, blocks, index{osa_iter}, uint32(options.projector_type), ...
                                            options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, int32(options.accuracy_factor));
                                else
                                    error('Unsupported projector type')
                                end
                            else
                                L = LL(index{osa_iter},:);
                                L = L';
                                L = L(:);
                                if options.projector_type == 1
                                    [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                            zmax, vaimennus, normalization, SinD, uint32(0), attenuation_correction, normalization_correction, ...
                                            randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, L, pseudot, det_per_ring, options.verbose, ...
                                            use_raw_data, uint32(2), ind_size, block1, blocks, uint32(0), uint32(options.projector_type));
                                elseif options.projector_type == 2
                                    [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                            zmax, vaimennus, normalization, SinD, uint32(0), attenuation_correction, normalization_correction, ...
                                            randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, L, pseudot, det_per_ring, options.verbose, ...
                                            use_raw_data, uint32(2), ind_size, block1, blocks, uint32(0), uint32(options.projector_type), options.tube_width_xy, ...
                                            x_center, y_center, z_center, options.tube_width_z, int32(options.accuracy_factor));
                                else
                                    error('Unsupported projector type')
                                end
                            end
                            %                             lor = reshape(lor,[],2);
                            %                             lor=repelem(uint32((lor(:,1))),lor(:,2));
                            if verLessThan('matlab','8.5')
                                lor = int32(repeat_elem(1:pituus(osa_iter),lor)');
                            else
                                lor = int32(repelem(1:pituus(osa_iter),lor)');
                            end
                            uu = double(Sino(index{osa_iter}));
                            
                            A_length = length(uu);
                            indices = int32(indices) + 1;
                            if verbose
                                tStart = tic;
                            end
                            if options.use_fsparse == false
                                A = sparse(double(lor),double(indices),double(alkiot), A_length, double(N));
                            else
                                A = fsparse(lor,indices,double(alkiot),[A_length double(N) length(alkiot)]);
                            end
                            clear indices alkiot lor
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['Sparse matrix formation took ' num2str(tElapsed) ' seconds'])
                            end
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
                                lor2 = [0; cumsum(uint64(lor_orth(pituus(osa_iter)+1:pituus(osa_iter + 1))))];
                            else
                                lor2 = [0; cumsum(uint64(lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1))))];
                            end
                            [A, ll] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, ...
                                size_x, zmax, vaimennus, normalization, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction,...
                                randoms_correction, lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1)), xy_index_input, ...
                                z_index_input, NSinos, L_input, pseudot, det_per_ring, options.verbose, ...
                                use_raw_data, uint32(0), lor2, summa(osa_iter), options.attenuation_phase, uint32(options.projector_type), ...
                                options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, int32(options.accuracy_factor));
                            uu = double(Sino(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                            if options.attenuation_phase
                                uu = uu ./ ll;
                            end
                            clear lor2
                        end
                        %                         osa_iter
                        if is_transposed
                            Summ = full(sum(A,2));
                        else
                            Summ = full(sum(A,1))';
                        end
                        Summ(Summ == 0) = options.epps;
                        if options.osem || options.ecosem || attenuation_phase
                            if verbose
                                tStart = tic;
                            end
                            im_vectors.OSEM_apu = OSEM_im(im_vectors.OSEM_apu, A, epps, uu, Summ, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.mramla
                            if verbose
                                tStart = tic;
                            end
                            im_vectors.MRAMLA_apu = MBSREM(im_vectors.MRAMLA_apu, options.U, pj3, A, epps, uu, epsilon_mramla, lam_mbsrem, iter, 0, 0, ...
                                is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['MRAMLA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.ramla
                            if verbose
                                tStart = tic;
                            end
                            im_vectors.RAMLA_apu = BSREM_subiter(im_vectors.RAMLA_apu, lam, epps, iter, A, uu, is_transposed);
                            if any(im_vectors.RAMLA_apu < 0)
                                error('Negative values in RAMLA, lower lambda value!')
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RAMLA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.rosem
                            if verbose
                                tStart = tic;
                            end
                            im_vectors.ROSEM_apu = ROSEM_subiter(im_vectors.ROSEM_apu, lam_rosem, iter, Summ, epps, A, uu, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.rbi
                            if verbose
                                tStart = tic;
                            end
                            im_vectors.RBI_apu = RBI_subiter(im_vectors.RBI_apu, A, uu, epps, Summ, 0, 0, D, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.drama
                            if verbose
                                tStart = tic;
                            end
                            im_vectors.DRAMA_apu = DRAMA_subiter(im_vectors.DRAMA_apu, lam_drama, epps, iter, Summ, osa_iter, A, uu, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['DRAMA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.cosem || options.ecosem
                            if verbose
                                tStart = tic;
                            end
                            [im_vectors.COSEM_apu, C_co] = COSEM_im(im_vectors.COSEM_apu, A, epps, uu, C_co, D, osa_iter, is_transposed);
                            %                             C_co(:,osa_iter) = Summ;
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
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
                        if options.acosem
                            if verbose
                                tStart = tic;
                            end
                            [im_vectors.ACOSEM_apu, C_aco] = ACOSEM_im(im_vectors.ACOSEM_apu, A, epps, uu, C_aco, D, options.h, osa_iter, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ACOSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.MRP && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            med = MRP(im_vectors.MRP_OSL_apu, medx, medy, medz, Nx, Ny, Nz, epps, tr_offsets, options.med_no_norm);
                            im_vectors.MRP_OSL_apu = OSL_OSEM(im_vectors.MRP_OSL_apu, Summ, options.beta_mrp_osem, med, epps, A, uu, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.MRP && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            med = MRP(im_vectors.MRP_MBSREM_apu, medx, medy, medz, Nx, Ny, Nz, epps, tr_offsets, options.med_no_norm);
                            im_vectors.MRP_MBSREM_apu = MBSREM(im_vectors.MRP_MBSREM_apu, options.U, pj3, A, epps, uu, epsilon_mramla, lam_mbsrem, iter, ...
                                options.beta_mrp_mbsrem, med, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['MBSREM MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.MRP && options.BSREM
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.ramla
                                im_vectors.MRP_BSREM_apu = im_vectors.RAMLA_apu;
                            else
                                im_vectors.MRP_BSREM_apu = BSREM_subiter(im_vectors.MRP_BSREM_apu, lam, epps, iter, A, uu, is_transposed);
                            end
                            if any(im_vectors.MRP_BSREM_apu < 0)
                                error('Negative values in BSREM, lower lambda value!')
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['BSREM MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.MRP && options.ROSEM_MAP
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.rosem
                                im_vectors.MRP_ROSEM_apu = im_vectors.ROSEM_apu;
                            else
                                im_vectors.MRP_ROSEM_apu = ROSEM_subiter(im_vectors.MRP_ROSEM_apu, lam_rosem, iter, Summ, epps, A, uu, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.MRP && options.RBI_MAP
                            if verbose
                                tStart = tic;
                            end
                            med = MRP(im_vectors.MRP_RBI_apu, medx, medy, medz, Nx, Ny, Nz, epps, tr_offsets, options.med_no_norm);
                            im_vectors.MRP_RBI_apu = RBI_subiter(im_vectors.MRP_RBI_apu, A, uu, epps, Summ, options.beta_mrp_rbi, med, D, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-MAP MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.MRP && any(options.COSEM_MAP)
                            if verbose
                                tStart = tic;
                            end
                            med = MRP(im_vectors.MRP_COSEM_apu, medx, medy, medz, Nx, Ny, Nz, epps, tr_offsets, options.med_no_norm);
                            if options.COSEM_MAP == 1
                                [im_vectors.MRP_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.MRP_COSEM_apu, D, options.beta_mrp_cosem, med, A, uu, epps, ...
                                    C_osl, options.h, options.COSEM_MAP, osa_iter, is_transposed);
                            else
                                [im_vectors.MRP_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.MRP_COSEM_apu, D, options.beta_mrp_cosem, med, A, uu, epps, ...
                                    C_osl, 0, options.COSEM_MAP, osa_iter, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-MAP MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.quad && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            med = Quadratic_prior(im_vectors.Quad_OSL_apu, tr_offsets, options.weights, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                            im_vectors.Quad_OSL_apu = OSL_OSEM(im_vectors.Quad_OSL_apu, Summ, options.beta_quad_osem, med, epps, A, uu, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.quad && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            med = Quadratic_prior(im_vectors.Quad_MBSREM_apu, tr_offsets, options.weights, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                            im_vectors.Quad_MBSREM_apu = MBSREM(im_vectors.Quad_MBSREM_apu, options.U, pj3, A, epps, uu, epsilon_mramla, lam_mbsrem, iter, ...
                                options.beta_quad_mbsrem, med, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['MBSREM Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.quad && options.BSREM
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.ramla
                                im_vectors.Quad_BSREM_apu = im_vectors.RAMLA_apu;
                            else
                                im_vectors.Quad_BSREM_apu = BSREM_subiter(im_vectors.Quad_BSREM_apu, lam, epps, iter, A, uu, is_transposed);
                            end
                            if any(im_vectors.Quad_BSREM_apu < 0)
                                error('Negative values in BSREM, lower lambda value')
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['BSREM Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.quad && options.ROSEM_MAP
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.rosem
                                im_vectors.Quad_ROSEM_apu = im_vectors.ROSEM_apu;
                            else
                                im_vectors.Quad_ROSEM_apu = ROSEM_subiter(im_vectors.Quad_ROSEM_apu, lam_rosem, iter, Summ, epps, A, uu, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.quad && options.RBI_MAP
                            if verbose
                                tStart = tic;
                            end
                            med = Quadratic_prior(im_vectors.Quad_RBI_apu, tr_offsets, options.weights, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                            im_vectors.Quad_RBI_apu = RBI_subiter(im_vectors.Quad_RBI_apu, A, uu, epps, Summ, options.beta_quad_rbi, med, D, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-MAP Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.quad && any(options.COSEM_MAP)
                            if verbose
                                tStart = tic;
                            end
                            med = Quadratic_prior(im_vectors.Quad_COSEM_apu, tr_offsets, options.weights, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                            if options.COSEM_MAP == 1
                                [im_vectors.Quad_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.Quad_COSEM_apu, D, options.beta_quad_cosem, med, A, uu, epps, ...
                                    C_osl, options.h, options.COSEM_MAP, osa_iter, is_transposed);
                            else
                                [im_vectors.Quad_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.Quad_COSEM_apu, D, options.beta_quad_cosem, med, A, uu, epps, ...
                                    C_osl, 0, options.COSEM_MAP, osa_iter, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-MAP Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.L && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            med = L_filter(im_vectors.L_OSL_apu, tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            im_vectors.L_OSL_apu = OSL_OSEM(im_vectors.L_OSL_apu, Summ, options.beta_L_osem, med, epps, A, uu, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.L && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            med = L_filter(im_vectors.L_MBSREM_apu, tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            im_vectors.L_MBSREM_apu = MBSREM(im_vectors.L_MBSREM_apu, options.U, pj3, A, epps, uu, epsilon_mramla, lam_mbsrem, iter, ...
                                options.beta_L_mbsrem, med, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['MBSREM L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.L && options.BSREM
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.ramla
                                im_vectors.L_BSREM_apu = im_vectors.RAMLA_apu;
                            else
                                im_vectors.L_BSREM_apu = BSREM_subiter(im_vectors.L_BSREM_apu, lam, epps, iter, A, uu, is_transposed);
                            end
                            if any(im_vectors.L_BSREM_apu < 0)
                                error('Negative values in BSREM, lower lambda value')
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['BSREM L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.L && options.ROSEM_MAP
                            if verbose
                                tStart = tic;
                            end
                            if iter == 1 && options.rosem
                                im_vectors.L_ROSEM_apu = im_vectors.ROSEM_apu;
                            else
                                im_vectors.L_ROSEM_apu = ROSEM_subiter(im_vectors.L_ROSEM_apu, lam_rosem, iter, Summ, epps, A, uu, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['ROSEM-MAP L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.L && options.RBI_MAP
                            if verbose
                                tStart = tic;
                            end
                            med = L_filter(im_vectors.L_RBI_apu, tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            im_vectors.L_RBI_apu = RBI_subiter(im_vectors.L_RBI_apu, A, uu, epps, Summ, options.beta_L_rbi, med, D, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-MAP L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.L && any(options.COSEM_MAP)
                            if verbose
                                tStart = tic;
                            end
                            med = L_filter(im_vectors.L_COSEM_apu, tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            if options.COSEM_MAP == 1
                                [im_vectors.L_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.L_COSEM_apu, D, options.beta_L_cosem, med, A, uu, epps, C_osl, ...
                                    options.h, options.COSEM_MAP, osa_iter, is_transposed);
                            else
                                [im_vectors.L_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.L_COSEM_apu, D, options.beta_L_cosem, med, A, uu, epps, C_osl, 0, ...
                                    options.COSEM_MAP, osa_iter, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-MAP L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.FMH && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            med = FMH(im_vectors.FMH_OSL_apu, tr_offsets, options.fmh_weights, options.weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            im_vectors.FMH_OSL_apu = OSL_OSEM(im_vectors.FMH_OSL_apu, Summ, options.beta_fmh_osem, med, epps, A, uu, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.FMH && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            med = FMH(im_vectors.FMH_MBSREM_apu, tr_offsets, options.fmh_weights, options.weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            im_vectors.FMH_MBSREM_apu = MBSREM(im_vectors.FMH_MBSREM_apu, options.U, pj3, A, epps, uu, epsilon_mramla, lam_mbsrem, iter, ...
                                options.beta_fmh_mbsrem, med, is_transposed);
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
                                im_vectors.FMH_BSREM_apu = BSREM_subiter(im_vectors.FMH_BSREM_apu, lam, epps, iter, A, uu, is_transposed);
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
                                im_vectors.FMH_ROSEM_apu = ROSEM_subiter(im_vectors.FMH_ROSEM_apu, lam_rosem, iter, Summ, epps, A, uu, is_transposed);
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
                            med = FMH(im_vectors.FMH_RBI_apu, tr_offsets, options.fmh_weights, options.weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            im_vectors.FMH_RBI_apu = RBI_subiter(im_vectors.FMH_RBI_apu, A, uu, epps, Summ, options.beta_fmh_rbi, med, D, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-MAP FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.FMH && any(options.COSEM_MAP)
                            if verbose
                                tStart = tic;
                            end
                            med = FMH(im_vectors.FMH_COSEM_apu, tr_offsets, options.fmh_weights, options.weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            if options.COSEM_MAP == 1
                                [im_vectors.FMH_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.FMH_COSEM_apu, D, options.beta_fmh_cosem, med, A, uu, epps, ...
                                    C_osl, options.h, options.COSEM_MAP, osa_iter, is_transposed);
                            else
                                [im_vectors.FMH_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.FMH_COSEM_apu, D, options.beta_fmh_cosem, med, A, uu, epps, ...
                                    C_osl, 0, options.COSEM_MAP, osa_iter, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-MAP FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.weighted_mean && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            med = Weighted_mean(im_vectors.Weighted_OSL_apu, tr_offsets, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.w_sum, options.med_no_norm);
                            im_vectors.Weighted_OSL_apu = OSL_OSEM(im_vectors.Weighted_OSL_apu, Summ, options.beta_weighted_osem, med, A, uu, epps,...
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
                            med = Weighted_mean(im_vectors.Weighted_MBSREM_apu, tr_offsets, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.w_sum, options.med_no_norm);
                            im_vectors.Weighted_MBSREM_apu = MBSREM(im_vectors.Weighted_MBSREM_apu, options.U, pj3, A, epps, uu, epsilon_mramla, ...
                                lam_mbsrem, iter, options.beta_weighted_mbsrem, med, is_transposed);
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
                                im_vectors.Weighted_BSREM_apu = BSREM_subiter(im_vectors.Weighted_BSREM_apu, lam, epps, iter, A, uu, is_transposed);
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
                                im_vectors.Weighted_ROSEM_apu = ROSEM_subiter(im_vectors.Weighted_ROSEM_apu, lam_rosem, A, uu, epps, iter, Summ, ...
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
                            med = Weighted_mean(im_vectors.Weighted_RBI_apu, tr_offsets, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.w_sum, options.med_no_norm);
                            im_vectors.Weighted_RBI_apu = RBI_subiter(im_vectors.Weighted_RBI_apu, A, uu, epps, Summ, options.beta_weighted_rbi, ...
                                med, D, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-MAP weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.weighted_mean && any(options.COSEM_MAP)
                            if verbose
                                tStart = tic;
                            end
                            med = Weighted_mean(im_vectors.Weighted_COSEM_apu, tr_offsets, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.w_sum, options.med_no_norm);
                            if options.COSEM_MAP == 1
                                [im_vectors.Weighted_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.Weighted_COSEM_apu, D, options.beta_weighted_cosem, ...
                                    med, A, uu, epps, C_osl, options.h, options.COSEM_MAP, osa_iter, is_transposed);
                            else
                                [im_vectors.Weighted_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.Weighted_COSEM_apu, D, options.beta_weighted_cosem, ...
                                    med, A, uu, epps, C_osl, 0, options.COSEM_MAP, osa_iter, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-MAP weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.TV && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            grad = TVpriorFinal(im_vectors.TV_OSL_apu, TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, ...
                                tr_offsets);
                            im_vectors.TV_OSL_apu = OSL_OSEM(im_vectors.TV_OSL_apu, Summ, options.beta_TV_osem, grad, A, uu, epps, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.TV && options.MBSREM
                            if verbose
                                tStart = tic;
                            end
                            grad = TVpriorFinal(im_vectors.TV_MBSREM_apu, TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, ...
                                tr_offsets);
                            im_vectors.TV_MBSREM_apu = MBSREM(im_vectors.TV_MBSREM_apu, options.U, pj3, A, epps, uu, epsilon_mramla, lam_mbsrem, ...
                                iter, options.beta_TV_mbsrem, grad, is_transposed);
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
                                im_vectors.TV_BSREM_apu = BSREM_subiter(im_vectors.TV_BSREM_apu, lam, epps, iter, A, uu, is_transposed);
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
                                im_vectors.TV_ROSEM_apu = ROSEM_subiter(im_vectors.TV_ROSEM_apu, lam_rosem, iter, Summ, epps, A, uu, is_transposed);
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
                            grad = TVpriorFinal(im_vectors.TV_RBI_apu, TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, ...
                                tr_offsets);
                            im_vectors.TV_RBI_apu = RBI_subiter(im_vectors.TV_RBI_apu, A, uu, epps, Summ, options.beta_TV_rbi, grad, D, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-MAP TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.TV && any(options.COSEM_MAP)
                            if verbose
                                tStart = tic;
                            end
                            grad = TVpriorFinal(im_vectors.TV_COSEM_apu, TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, ...
                                tr_offsets);
                            if options.COSEM_MAP == 1
                                [im_vectors.TV_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.TV_COSEM_apu, D, options.beta_TV_cosem, grad, A, uu, epps, ...
                                    C_aco, options.h, options.COSEM_MAP, osa_iter, is_transposed);
                            else
                                [im_vectors.TV_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.TV_COSEM_apu, D, options.beta_TV_cosem, grad, A, uu, epps, ...
                                    C_co, 0, options.COSEM_MAP, osa_iter, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-MAP TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.AD && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            if osa_iter > 1
                                med = AD(im_vectors.AD_OSL_apu, FluxType, Nx, Ny, Nz, options);
                                im_vectors.AD_OSL_apu = OSL_OSEM(im_vectors.AD_OSL_apu, Summ, options.beta_ad_osem, med, epps, A, uu, is_transposed);
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
                            med = AD(im_vectors.AD_MBSREM_apu, FluxType, Nx, Ny, Nz, options);
                            im_vectors.AD_MBSREM_apu = MBSREM(im_vectors.AD_MBSREM_apu, options.U, pj3, A, epps, uu, epsilon_mramla, lam_mbsrem, ...
                                iter, options.beta_ad_mbsrem, med, is_transposed);
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
                                im_vectors.AD_BSREM_apu = BSREM_subiter(im_vectors.AD_BSREM_apu, lam, epps, iter, A, uu, is_transposed);
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
                                im_vectors.AD_ROSEM_apu = ROSEM_subiter(im_vectors.AD_ROSEM_apu, lam_rosem, iter, Summ, epps, A, uu, is_transposed);
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
                            med = AD(im_vectors.AD_RBI_apu, FluxType, Nx, Ny, Nz, options);
                            im_vectors.AD_RBI_apu = RBI_subiter(im_vectors.AD_RBI_apu, A, uu, epps, Summ, options.beta_ad_rbi, med, D, is_transposed);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['RBI-MAP AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.AD && any(options.COSEM_MAP)
                            if verbose
                                tStart = tic;
                            end
                            med = AD(im_vectors.AD_COSEM_apu, FluxType, Nx, Ny, Nz, options);
                            if options.COSEM_MAP == 1
                                [im_vectors.AD_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.AD_COSEM_apu, D, options.beta_ad_cosem, med, A, uu, epps, ...
                                    C_osl, options.h, options.COSEM_MAP, osa_iter, is_transposed);
                            else
                                [im_vectors.AD_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.AD_COSEM_apu, D, options.beta_ad_cosem, med, A, uu, epps, ...
                                    C_osl, 0, options.COSEM_MAP, osa_iter, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-MAP AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.APLS && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            grad = TVpriorFinal(im_vectors.APLS_OSL_apu, [], Nx, Ny, Nz, true, options, 4);
                            im_vectors.APLS_OSL_apu = OSL_OSEM(im_vectors.APLS_OSL_apu, Summ, options.beta_APLS_osem, grad, A, uu, epps, is_transposed);
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
                            im_vectors.APLS_MBSREM_apu = MBSREM(im_vectors.APLS_MBSREM_apu, options.U, pj3, A, epps, uu, epsilon_mramla, lam_mbsrem, ...
                                iter, options.beta_APLS_mbsrem, grad, is_transposed);
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
                                im_vectors.APLS_BSREM_apu = BSREM_subiter(im_vectors.APLS_BSREM_apu, lam, epps, iter, A, uu, is_transposed);
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
                                im_vectors.APLS_ROSEM_apu = ROSEM_subiter(im_vectors.APLS_ROSEM_apu, lam_rosem, iter, Summ, epps, A, uu, is_transposed);
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
                            im_vectors.APLS_RBI_apu = RBI_subiter(im_vectors.APLS_RBI_apu, A, uu, epps, Summ, options.beta_APLS_rbi, grad, D, ...
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
                                    epps, C_aco, options.h, options.COSEM_MAP, osa_iter, is_transposed);
                            else
                                [im_vectors.APLS_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.APLS_COSEM_apu, D, options.beta_APLS_cosem, grad, A, uu, ...
                                    epps, C_co, 0, options.COSEM_MAP, osa_iter, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-MAP APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.TGV && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            grad = TGV(im_vectors.TGV_OSL_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                            im_vectors.TGV_OSL_apu = OSL_OSEM(im_vectors.TGV_OSL_apu, Summ, options.beta_TGV_osem, grad, A, uu, epps, is_transposed);
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
                            im_vectors.TGV_MBSREM_apu = MBSREM(im_vectors.TGV_MBSREM_apu, options.U, pj3, A, epps, uu, epsilon_mramla, lam_mbsrem, ...
                                iter, options.beta_TGV_mbsrem, grad, is_transposed);
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
                                im_vectors.TGV_BSREM_apu = BSREM_subiter(im_vectors.TGV_BSREM_apu, lam, epps, iter, A, uu, is_transposed);
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
                                im_vectors.TGV_ROSEM_apu = ROSEM_subiter(im_vectors.TGV_ROSEM_apu, lam_rosem, iter, Summ, epps, A, uu, is_transposed);
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
                            im_vectors.TGV_RBI_apu = RBI_subiter(im_vectors.TGV_RBI_apu, A, uu, epps, Summ, options.beta_TGV_rbi, grad, D, is_transposed);
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
                                    epps, C_aco, options.h, options.COSEM_MAP, osa_iter, is_transposed);
                            else
                                [im_vectors.TGV_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.TGV_COSEM_apu, D, options.beta_TGV_cosem, grad, A, uu, ...
                                    epps, C_co, 0, options.COSEM_MAP, osa_iter, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-MAP TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        if options.NLM && options.OSL_OSEM
                            if verbose
                                tStart = tic;
                            end
                            med = NLM(im_vectors.NLM_OSL_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                options.sigma, epps, Nx, Ny, Nz, options);
                            im_vectors.NLM_OSL_apu = OSL_OSEM(im_vectors.NLM_OSL_apu, Summ, options.beta_NLM_osem, med, epps, A, uu, is_transposed);
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
                            im_vectors.NLM_MBSREM_apu = MBSREM(im_vectors.NLM_MBSREM_apu, options.U, pj3, A, epps, uu, epsilon_mramla, lam_mbsrem, ...
                                iter, options.beta_NLM_mbsrem, med, is_transposed);
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
                                im_vectors.NLM_BSREM_apu = BSREM_subiter(im_vectors.NLM_BSREM_apu, lam, epps, iter, A, uu, is_transposed);
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
                                im_vectors.NLM_ROSEM_apu = ROSEM_subiter(im_vectors.NLM_ROSEM_apu, lam_rosem, iter, Summ, epps, A, uu, is_transposed);
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
                            im_vectors.NLM_RBI_apu = RBI_subiter(im_vectors.NLM_RBI_apu, A, uu, epps, Summ, options.beta_NLM_rbi, med, D, is_transposed);
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
                                    epps, C_aco, options.h, options.COSEM_MAP, osa_iter, is_transposed);
                            else
                                [im_vectors.NLM_COSEM_apu, C_osl] = COSEM_OSL(im_vectors.NLM_COSEM_apu, D, options.beta_NLM_cosem, med, A, uu, ...
                                    epps, C_co, 0, options.COSEM_MAP, osa_iter, is_transposed);
                            end
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['COSEM-MAP NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        end
                        clear A
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
                        med = MRP(im_vectors.MRP_BSREM_apu, medx, medy, medz, Nx, Ny, Nz, epps, tr_offsets, options.med_no_norm);
                        im_vectors.MRP_BSREM(:,iter+1) = BSREM_iter(im_vectors.MRP_BSREM_apu, lam, iter, options.beta_mrp_bsrem, med, epps);
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
                        med = MRP(im_vectors.MRP_ROSEM_apu, medx, medy, medz, Nx, Ny, Nz, epps, tr_offsets, options.med_no_norm);
                        im_vectors.MRP_ROSEM(:,iter+1) = BSREM_iter(im_vectors.MRP_ROSEM_apu, lam_rosem, iter, options.beta_mrp_rosem, med, epps);
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
                        med = Quadratic_prior(im_vectors.Quad_BSREM_apu, tr_offsets, options.weights, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                        im_vectors.Quad_BSREM(:,iter+1) = BSREM_iter(im_vectors.Quad_BSREM_apu, lam, iter, options.beta_quad_bsrem, med, epps);
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
                        med = Quadratic_prior(im_vectors.Quad_ROSEM_apu, tr_offsets, options.weights, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                        im_vectors.Quad_ROSEM(:,iter+1) = BSREM_iter(im_vectors.Quad_ROSEM_apu, lam_rosem, iter, options.beta_quad_rosem, med, epps);
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
                        med = L_filter(im_vectors.L_BSREM_apu, tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                        im_vectors.L_BSREM(:,iter+1) = BSREM_iter(im_vectors.L_BSREM_apu, lam, iter, options.beta_L_bsrem, med, epps);
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
                        med = L_filter(im_vectors.L_ROSEM_apu, tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                        im_vectors.L_ROSEM(:,iter+1) = BSREM_iter(im_vectors.L_ROSEM_apu, lam_rosem, iter, options.beta_L_rosem, med, epps);
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
                        med = FMH(im_vectors.FMH_BSREM_apu, tr_offsets, options.fmh_weights, options.weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                            options.med_no_norm);
                        im_vectors.FMH_BSREM(:,iter+1) = BSREM_iter(im_vectors.FMH_BSREM_apu, lam, iter, options.beta_fmh_bsrem, med, epps);
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
                        med = FMH(im_vectors.FMH_ROSEM_apu, tr_offsets, options.fmh_weights, options.weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                            options.med_no_norm);
                        im_vectors.FMH_ROSEM(:,iter+1) = BSREM_iter(im_vectors.FMH_ROSEM_apu, lam_rosem, iter, options.beta_fmh_rosem, med, epps);
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
                        med = Weighted_mean(im_vectors.Weighted_BSREM_apu, tr_offsets, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                            options.mean_type, epps, options.w_sum, options.med_no_norm);
                        im_vectors.Weighted_BSREM(:,iter+1) = BSREM_iter(im_vectors.Weighted_BSREM_apu, lam, iter, options.beta_weighted_bsrem, ...
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
                        med = Weighted_mean(im_vectors.Weighted_ROSEM_apu, tr_offsets, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                            options.mean_type, epps, options.w_sum, options.med_no_norm);
                        im_vectors.Weighted_ROSEM(:,iter+1) = BSREM_iter(im_vectors.Weighted_ROSEM_apu, lam_rosem, iter, options.beta_weighted_rosem, ...
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
                        grad = TVpriorFinal(im_vectors.TV_BSREM_apu, TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, tr_offsets);
                        im_vectors.TV_BSREM(:, iter + 1) = BSREM_iter(im_vectors.TV_BSREM_apu, lam, iter, options.beta_TV_bsrem, grad, epps);
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
                        grad = TVpriorFinal(im_vectors.TV_ROSEM_apu, TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, tr_offsets);
                        im_vectors.TV_ROSEM(:, iter + 1) = BSREM_iter(im_vectors.TV_ROSEM_apu, lam_rosem, iter, options.beta_TV_rosem, grad, epps);
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
                        med = AD(im_vectors.AD_BSREM_apu, FluxType, Nx, Ny, Nz, options);
                        im_vectors.AD_BSREM(:,iter+1) = BSREM_iter(im_vectors.AD_BSREM_apu, lam, iter, options.beta_ad_bsrem, med, epps);
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
                        med = AD(im_vectors.AD_ROSEM_apu, FluxType, Nx, Ny, Nz, options);
                        im_vectors.AD_ROSEM(:, iter + 1) = BSREM_iter(im_vectors.AD_ROSEM_apu, lam_rosem, iter, options.beta_ad_rosem, med, epps);
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
                        im_vectors.APLS_BSREM(:, iter + 1) = BSREM_iter(im_vectors.APLS_BSREM_apu, lam, iter, options.beta_APLS_bsrem, grad, epps);
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
                        im_vectors.APLS_ROSEM(:, iter + 1) = BSREM_iter(im_vectors.APLS_ROSEM_apu, lam_rosem, iter, options.beta_APLS_rosem, grad, ...
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
                        im_vectors.TGV_BSREM(:, iter + 1) = BSREM_iter(im_vectors.TGV_BSREM_apu, lam, iter, options.beta_TGV_bsrem, grad, epps);
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
                        im_vectors.TGV_ROSEM(:, iter + 1) = BSREM_iter(im_vectors.TGV_ROSEM_apu, lam_rosem, iter, options.beta_TGV_rosem, grad, epps);
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
                        im_vectors.NLM_BSREM(:, iter + 1) = BSREM_iter(im_vectors.NLM_BSREM_apu, lam_rosem, iter, options.beta_NLM_bsrem, med, epps);
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
                        im_vectors.NLM_ROSEM(:, iter + 1) = BSREM_iter(im_vectors.NLM_ROSEM_apu, lam_rosem, iter, options.beta_NLM_rosem, med, epps);
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
            elseif options.reconstruction_method == 5
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
                    single(NSlices), size_x, zmax, NSinos, vaimennus, pituus, uint32(attenuation_correction), uint32(Niter), uint32(subsets), rekot, ...
                    single(epps), single(full(Sino)), single(options.x0(:)), lor_a, summa, xy_index, z_index, LL, pseudot, det_per_ring, ...
                    uint8(use_raw_data), options.verbose, device);
            elseif options.reconstruction_method == 4
                no_norm = false;
                if ~use_raw_data
                    if isempty(pseudot)
                        pseudot = uint32(0);
                    end
                end
                for iter = 1 : Niter
                    if OS_bool
                        if verbose
                            tStart_iter = tic;
                        end
                        if iter == 1
                            f_Summ = zeros(Nx*Ny*Nz,subsets);
                        end
                        for osa_iter = 1 : subsets
                            if verbose
                                tStart = tic;
                            end
                            uu = double(full(Sino(pituus(osa_iter)+1:pituus(osa_iter + 1))));
                            if randoms_correction
                                dd = double(full(SinD(pituus(osa_iter)+1:pituus(osa_iter + 1))));
                            else
                                dd = 0;
                            end
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
                            
                            [Summ, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, ...
                                size_x, zmax, vaimennus, normalization, dd, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction,...
                                lor_a_input, xy_index_input, ...
                                z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
                                options.verbose, (use_raw_data), uint32(1), epps, uu, im_vectors.OSEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, ...
                                options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, int32(options.accuracy_factor));
                            
                            if iter == 1
                                f_Summ(:,osa_iter) = Summ + epps;
                            end
                            if options.osem
                                im_vectors.OSEM_apu = OSEM_im(im_vectors.OSEM_apu, rhs, f_Summ(:,osa_iter), epps);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.ramla
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, lam, epps, iter, f_Summ(:,osa_iter), rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error(['Negative values in RAMLA, lower lambda value! lambda <= ' num2str(min(1./f_Summ(:,osa_iter)))])
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['RAMLA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.rosem
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.MRP && options.OSL_OSEM
                                med = MRP(im_vectors.OSEM_apu, medx, medy, medz, Nx, Ny, Nz, epps, tr_offsets, options.med_no_norm);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_mrp_osem, med, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.MRP && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error('Negative values in BSREM, lower lambda value!')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.MRP && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM-MAP MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.quad && options.OSL_OSEM
                                med = Quadratic_prior(im_vectors.OSEM_apu, tr_offsets, options.weights, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_quad_osem, med, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSL Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.quad && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error('Negative values in BSREM, lower lambda value')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.quad && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM-MAP Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.L && options.OSL_OSEM
                                med = L_filter(im_vectors.OSEM_apu, tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_L_osem, med, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSL L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.L && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error('Negative values in BSREM, lower lambda value')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.L && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM-MAP L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.FMH && options.OSL_OSEM
                                if verbose
                                    tStart = tic;
                                end
                                med = FMH(im_vectors.OSEM_apu, tr_offsets, options.fmh_weights, options.weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                    options.med_no_norm);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_fmh_osem, med, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSL FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.FMH && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error('Negative values in BSREM, lower lambda value')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.FMH && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM-MAP FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.weighted_mean && options.OSL_OSEM
                                med = Weighted_mean(im_vectors.OSEM_apu, tr_offsets, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                    options.mean_type, epps, options.w_sum, options.med_no_norm);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_weighted_osem, med, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSL weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.weighted_mean && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error('Negative values in BSREM, lower lambda value')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.weighted_mean && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, lam_rosem, epps, iter, f_Summ(:,osa_iter), rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM-MAP weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.TV && options.OSL_OSEM
                                grad = TVpriorFinal(im_vectors.OSEM_apu, TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, ...
                                    tr_offsets);
                                im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_TV_osem, grad, epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSL TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.TV && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error('Negative values in BSREM, lower lambda value')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.TV && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['ROSEM-MAP TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.AD && options.OSL_OSEM
                                if osa_iter > 1
                                    med = AD(im_vectors.OSEM_apu, FluxType, Nx, Ny, Nz, options);
                                    im_vectors.OSEM_apu = OSL_OSEM(im_vectors.OSEM_apu, f_Summ(:,osa_iter), options.beta_ad_osem, med, epps, rhs);
                                else
                                    im_vectors.OSEM_apu = OSEM_im(im_vectors.OSEM_apu, rhs, f_Summ(:,osa_iter), epps);
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['OSL AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.AD && options.BSREM
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error('Negative values in BSREM, lower lambda value')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.AD && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
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
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error('Negative values in BSREM, lower lambda value')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.APLS && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
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
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error('Negative values in BSREM, lower lambda value')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.TGV && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
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
                                im_vectors.OSEM_apu = BSREM_subiter(im_vectors.OSEM_apu, lam, epps, iter, rhs);
                                if any(im_vectors.OSEM_apu < 0)
                                    error('Negative values in BSREM, lower lambda value')
                                end
                                if verbose
                                    tElapsed = toc(tStart);
                                    disp(['BSREM NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
                                end
                            elseif options.NLM && options.ROSEM_MAP
                                im_vectors.OSEM_apu = ROSEM_subiter(im_vectors.OSEM_apu, lam_rosem, iter, f_Summ(:,osa_iter), epps, rhs);
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
                            med = MRP(im_vectors.OSEM_apu, medx, medy, medz, Nx, Ny, Nz, epps, tr_offsets, options.med_no_norm);
                            im_vectors.MRP_BSREM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, lam, iter, options.beta_mrp_bsrem, med, epps);
                        elseif options.MRP && options.ROSEM_MAP
                            med = MRP(im_vectors.OSEM_apu, medx, medy, medz, Nx, Ny, Nz, epps, tr_offsets, options.med_no_norm);
                            im_vectors.MRP_ROSEM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, lam_rosem, iter, options.beta_mrp_rosem, med, epps);
                        elseif options.quad && options.OSL_OSEM
                            im_vectors.Quad_OSL(:, iter + 1) = im_vectors.OSEM_apu;
                        elseif options.quad && options.BSREM
                            med = Quadratic_prior(im_vectors.OSEM_apu, tr_offsets, options.weights, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                            im_vectors.Quad_BSREM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, lam, iter, options.beta_quad_bsrem, med, epps);
                        elseif options.quad && options.ROSEM_MAP
                            med = Quadratic_prior(im_vectors.OSEM_apu, tr_offsets, options.weights, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                            im_vectors.Quad_ROSEM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, lam_rosem, iter, options.beta_quad_rosem, med, epps);
                        elseif options.L && options.OSL_OSEM
                            im_vectors.L_OSL(:, iter + 1) = im_vectors.L_OSL_apu;
                        elseif options.L && options.BSREM
                            med = L_filter(im_vectors.OSEM_apu, tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            im_vectors.L_BSREM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, lam, iter, options.beta_L_bsrem, med, epps);
                        elseif options.L && options.ROSEM_MAP
                            med = L_filter(im_vectors.OSEM_apu, tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            im_vectors.L_ROSEM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, lam_rosem, iter, options.beta_L_rosem, med, epps);
                        elseif options.FMH && options.OSL_OSEM
                            im_vectors.FMH_OSL(:, iter + 1) = im_vectors.OSEM_apu;
                        elseif options.FMH && options.BSREM
                            med = FMH(im_vectors.OSEM_apu, tr_offsets, options.fmh_weights, options.weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            im_vectors.FMH_BSREM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, lam, iter, options.beta_fmh_bsrem, med, epps);
                        elseif options.FMH && options.ROSEM_MAP
                            med = FMH(im_vectors.OSEM_apu, tr_offsets, options.fmh_weights, options.weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            im_vectors.FMH_ROSEM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, lam_rosem, iter, options.beta_fmh_rosem, med, epps);
                        elseif options.weighted_mean && options.OSL_OSEM
                            im_vectors.Weighted_OSL(:, iter + 1) = im_vectors.Weighted_OSL_apu;
                        elseif options.weighted_mean && options.BSREM
                            med = Weighted_mean(im_vectors.OSEM_apu, tr_offsets, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.w_sum, options.med_no_norm);
                            im_vectors.Weighted_BSREM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, lam, iter, options.beta_weighted_bsrem, ...
                                med, epps);
                        elseif options.weighted_mean && options.ROSEM_MAP
                            med = Weighted_mean(im_vectors.OSEM_apu, tr_offsets, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.w_sum, options.med_no_norm);
                            im_vectors.Weighted_ROSEM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, lam_rosem, iter, options.beta_weighted_rosem, ...
                                med, epps);
                        elseif options.TV && options.OSL_OSEM
                            im_vectors.TV_OSL(:, iter + 1) = im_vectors.OSEM_apu;
                        elseif options.TV && options.BSREM
                            grad = TVpriorFinal(im_vectors.OSEM_apu, TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, tr_offsets);
                            im_vectors.TV_BSREM(:, iter + 1) = BSREM_iter(im_vectors.OSEM_apu, lam, iter, options.beta_TV_bsrem, grad, epps);
                        elseif options.TV && options.ROSEM_MAP
                            grad = TVpriorFinal(im_vectors.OSEM_apu, TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, tr_offsets);
                            im_vectors.TV_ROSEM(:, iter + 1) = BSREM_iter(im_vectors.OSEM_apu, lam_rosem, iter, options.beta_TV_rosem, grad, epps);
                        elseif options.AD && options.OSL_OSEM
                            im_vectors.AD_OSL(:, iter + 1) = im_vectors.AD_OSL_apu;
                        elseif options.AD && options.BSREM
                            med = AD(im_vectors.OSEM_apu, FluxType, Nx, Ny, Nz, options);
                            im_vectors.AD_BSREM(:,iter+1) = BSREM_iter(im_vectors.OSEM_apu, lam, iter, options.beta_ad_bsrem, med, epps);
                        elseif options.AD && options.ROSEM_MAP
                            med = AD(im_vectors.OSEM_apu, FluxType, Nx, Ny, Nz, options);
                            im_vectors.AD_ROSEM(:, iter + 1) = BSREM_iter(im_vectors.OSEM_apu, lam_rosem, iter, options.beta_ad_rosem, med, epps);
                        elseif options.APLS && options.OSL_OSEM
                            im_vectors.APLS_OSL(:, iter + 1) = im_vectors.APLS_OSL_apu;
                        elseif options.APLS && options.BSREM
                            grad = TVpriorFinal(im_vectors.OSEM_apu, 0, Nx, Ny, Nz, true, options, 4);
                            im_vectors.APLS_BSREM(:, iter + 1) = BSREM_iter(im_vectors.OSEM_apu, lam, iter, options.beta_APLS_bsrem, grad, epps);
                        elseif options.APLS && options.ROSEM_MAP
                            grad = TVpriorFinal(im_vectors.OSEM_apu, 0, Nx, Ny, Nz, true, options, 4);
                            im_vectors.APLS_ROSEM(:, iter + 1) = BSREM_iter(im_vectors.OSEM_apu, lam_rosem, iter, options.beta_APLS_rosem, grad, ...
                                epps);
                        elseif options.TGV && options.OSL_OSEM
                            im_vectors.TGV_OSL(:, iter + 1) = im_vectors.OSEM_apu;
                        elseif options.TGV && options.BSREM
                            grad = TGV(im_vectors.OSEM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                            im_vectors.TGV_BSREM(:, iter + 1) = BSREM_iter(im_vectors.OSEM_apu, lam, iter, options.beta_TGV_bsrem, grad, epps);
                        elseif options.TGV && options.ROSEM_MAP
                            grad = TGV(im_vectors.OSEM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, Nx, Ny, Nz);
                            im_vectors.TGV_ROSEM(:, iter + 1) = BSREM_iter(im_vectors.OSEM_apu, lam_rosem, iter, options.beta_TGV_rosem, grad, epps);
                        elseif options.NLM && options.OSL_OSEM
                            im_vectors.NLM_OSL(:, iter + 1) = im_vectors.NLM_OSL_apu;
                        elseif options.NLM && options.BSREM
                            med = NLM(im_vectors.OSEM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                options.sigma, epps, Nx, Ny, Nz, options);
                            im_vectors.NLM_BSREM(:, iter + 1) = BSREM_iter(im_vectors.OSEM_apu, lam_rosem, iter, options.beta_NLM_bsrem, med, epps);
                        elseif options.NLM && options.ROSEM_MAP
                            med = NLM(im_vectors.OSEM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
                                options.sigma, epps, Nx, Ny, Nz, options);
                            im_vectors.NLM_ROSEM(:, iter + 1) = BSREM_iter(im_vectors.OSEM_apu, lam_rosem, iter, options.beta_NLM_rosem, med, epps);
                        end
                        if verbose
                            tElapsed = toc(tStart_iter);
                            disp(['OS iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                        else
                            disp(['OS iteration ' num2str(iter) ' finished'])
                        end
                        no_norm = true;
                    end
                    if MLEM_bool
                        if verbose
                            tStart = tic;
                        end
                        if iter == 1
                            if no_norm && OS_bool
                                f_Summ_ml = sum(f_Summ,2);
                            else
                                f_Summ_ml = zeros(Nx*Ny*Nz,1);
                            end
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
                        
                        [Summ, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, ...
                            size_x, zmax, vaimennus, normalization, dd, pituus(end) - pituus(1), attenuation_correction, normalization_correction, randoms_correction, lor_a, xy_index, z_index, NSinos, ...
                            LL, pseudot, det_per_ring, options.verbose, (use_raw_data), uint32(1), epps, double(Sino), ...
                            im_vectors.MLEM_apu, uint32(options.projector_type), no_norm, options.precompute_lor, ...
                            options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, int32(options.accuracy_factor));
                        
                        if iter == 1 && Niter > 1 && ~OS_bool
                            f_Summ_ml = Summ;
                        end
                        if options.mlem
                            im_vectors.MLEM(:,iter + 1) = MLEM_im(im_vectors.MLEM_apu, f_Summ_ml, epps, rhs);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['MLEM iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            else
                                disp(['MLEM iteration ' num2str(iter) ' finished'])
                            end
                        elseif options.MRP && options.OSL_MLEM
                            med = MRP(im_vectors.MLEM_apu, medx, medy, medz, Nx, Ny, Nz, epps, tr_offsets, options.med_no_norm);
                            im_vectors.MRP_MLEM(:,iter + 1) = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_mrp_osem, med, epps, rhs);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL MRP iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.quad && options.OSL_MLEM
                            med = Quadratic_prior(im_vectors.MLEM_apu, tr_offsets, options.weights, options.weights_quad, Nx, Ny, Nz, Ndx, Ndy, Ndz);
                            im_vectors.Quad_MLEM(:, iter + 1) = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_quad_osem, med, epps, rhs);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL Quadratic iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.L && options.OSL_MLEM
                            med = L_filter(im_vectors.MLEM_apu, tr_offsets, options.a_L, Nx, Ny, Nz, Ndx, Ndy, Ndz, epps, options.med_no_norm);
                            im_vectors.L_MLEM(:, iter + 1) = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_L_osem, med, epps, rhs);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL L-filter iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.FMH && options.OSL_MLEM
                            med = FMH(im_vectors.MLEM_apu, tr_offsets, options.fmh_weights, options.weights, Nx, Ny, Nz, N, Ndx, Ndy, Ndz, epps, ...
                                options.med_no_norm);
                            im_vectors.FMH_MLEM(:, iter + 1) = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_fmh_osem, med, epps, rhs);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL FMH iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.weighted_mean && options.OSL_MLEM
                            med = Weighted_mean(im_vectors.MLEM_apu, tr_offsets, options.weighted_weights, Nx, Ny, Nz, Ndx, Ndy, Ndz, ...
                                options.mean_type, epps, options.w_sum, options.med_no_norm);
                            im_vectors.Weighted_MLEM(:, iter + 1) = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_weighted_osem, med, rhs);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL weighted mean iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.TV && options.OSL_MLEM
                            grad = TVpriorFinal(im_vectors.MLEM_apu, TVdata, Nx, Ny, Nz, options.TV_use_anatomical, options, options.TVtype, ...
                                tr_offsets);
                            im_vectors.TV_MLEM(:, iter + 1) = OSL_OSEM(im_vectors.MLEM_apu, f_Summ_ml(:,osa_iter), options.beta_TV_osem, grad, rhs);
                            if verbose
                                tElapsed = toc(tStart);
                                disp(['OSL TV sub-iteration ' num2str(iter) ' took ' num2str(tElapsed) ' seconds'])
                            end
                        elseif options.AD && options.OSL_MLEM
                            if iter > 1
                                med = AD(im_vectors.MLEM_apu, FluxType, Nx, Ny, Nz, options);
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
            if options.reconstruction_method ~= 2 && options.reconstruction_method ~= 3
                im_vectors = reshape_vectors(im_vectors, options);
            end
            pz = images_to_cell(im_vectors, llo, pz, options);
            if partitions > 1
                disp(['Reconstructions for timestep ' num2str(llo) ' completed'])
            end
        end
    elseif options.reconstruction_method == 2
        %%
        %         maksimi = zeros(subsets,1,'uint32');
        %         for osa_iter = 1 : subsets
        %             maksimi(osa_iter) = uint32(max(lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1))));
        %         end
        options = double_to_single(options);
        
        if partitions == 1
            clear SinM
            SinM{1} = single(full(Sino));
            clear Sino
        end
        if issparse(SinM{1})
            for kk = 1 : length(SinM)
                SinM{kk} = single(full(SinM{kk}));
            end
            clear Sino
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
        
        
        tube_width = single(options.tube_width_xy);
        %         if options.projector_type == 1
        kernel_file = 'siddon_kernel_matrixfree_GPU';
        %         kernel_file = 'siddon_kernel_matrixfree_GPU - Copy (2).cl';
        filename = 'OMEGA_matrix_free_OpenCL_binary_device';
        %             x_center = xx(end);
        %             y_center = yy(end);
        %         elseif options.projector_type == 2
        %             kernel_file = 'orthogonal_kernel_matrixfree.cl';
        %             filename = 'OMEGA_matrix_free_orthogonal_OpenCL_binary_device';
        % %             x_center = xx(1:end-1)'+dx/2;
        % %             y_center = yy(1:end-1)'+dy/2;
        %         end
        kernel_path = which(kernel_file);
        kernel_path = strrep(kernel_path, '\', '/');
        filename = [kernel_path(1:end-length(kernel_file)), filename];
        %         filename = [kernel_path(1:end-length(kernel_file) + 3), filename];
        tic
        [pz] = OpenCL_matrixfree( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end) , NSinos, ...
            single(NSlices), size_x, zmax, NSinos, options.verbose, LL, pseudot, det_per_ring, device, uint8(use_raw_data), ...
            filename, uint32(0), options.force_build, vaimennus, normalization, pituus, uint32(attenuation_correction), uint32(options.normalization_correction), uint32(Niter), uint32(subsets), ...
            uint8(rekot), single(epps), lor_a, xy_index, z_index, any(sum(rekot(11:end))), tube_width, x_center, y_center, SinDelayed, randoms, options, ...
            SinM, uint32(partitions), uint32(options.projector_type), int32(options.accuracy_factor));
        toc
    elseif options.reconstruction_method == 3
        options = double_to_single(options);
        
        if partitions == 1
            clear SinM
            SinM{1} = single(full(Sino));
            if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
                apu = single(full(SinDelayed));
                SinDelayed = cell(1,1);
                SinDelayed{1} = apu;
            end
            clear Sino apu
        end
        if issparse(SinM{1})
            for kk = 1 : length(SinM)
                SinM{kk} = single(full(SinM{kk}));
            end
            clear Sino
        end
        if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
            if issparse(SinDelayed{1})
                for kk = 1 : length(SinDelayed)
                    SinDelayed{kk} = single(full(SinDelayed{kk}));
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
        
        %         options.normalization_correction = false;
        %         SinM{1} = SinM{1} .* normalization;
        %         normalization = single(0);
        
        tube_width_xy = single(options.tube_width_xy);
        crystal_size_z = single(options.tube_width_z);
        if options.projector_type == 1 && options.precompute_lor
            %         kernel_file = 'siddon_kernel_matrixfree_GPU - Copy (2).cl';
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
        %         kernel_file = 'multidevice_kernels.cl';
        % %         kernel_file = 'siddon_kernel_matrixfree_GPU - Copy (2).cl';
        %         filename = 'OMEGA_matrix_free_OpenCL_binary_device';
        filename = [header_directory, filename];
        header_directory = strcat('-I "', header_directory);
        header_directory = strcat(header_directory,'"');
        tic
        [pz] = OpenCL_matrixfree_multi_gpu( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end), ...
            single(NSlices), size_x, zmax, options.verbose, LL, pseudot, det_per_ring, uint32(options.use_device), filename, uint8(use_raw_data), ...
            single(options.cpu_to_gpu_factor), uint32(1), header_directory, vaimennus, normalization, pituus, uint32(attenuation_correction), ...
            uint32(normalization_correction), lor_a, xy_index, z_index, tube_width_xy, crystal_size_z, x_center, y_center, z_center, SinDelayed, randoms, ...
            uint32(options.projector_type), options.precompute_lor, int32(options.accuracy_factor), NSinos, uint16(TotSinos), uint32(Niter), uint32(subsets), uint8(rekot), single(epps), SinM, ...
            uint32(partitions), options.osem, options.force_build, options);
        toc
        
        for kk = 1 : length(pz)
            if ~isempty(pz{kk})
                pz{kk} = reshape(pz{kk}, Nx, Ny, Nz, Niter + 1);
            end
        end
        pz{end + 1} = [];
    else
        error('Unsupported reconstruction method.');
    end
end

pz = save_image_properties(options, pz, subsets);
% clear mex

end
