function [options, pz] = custom_prior_prepass(options)
%% Prepass phase for the custom prior file
% Computes all the necessary variables needed for the reconstruction
% process
%

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

if ~isfield(options,'use_Inveon')
    options.use_Inveon = 0;
end

if options.precompute_lor == false && options.reconstruction_method == 3
    error('precompute_lor must be set to true if using method 3')
end

Nx = options.Nx;
Ny = options.Ny;
Nz = options.Nz;
Niter = options.Niter;
subsets = options.subsets;
epps = options.epps;
attenuation_correction = options.attenuation_correction;
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
cr_pz = options.cr_pz;
use_raw_data = options.use_raw_data;
TotSinos = options.TotSinos;
attenuation_datafile = options.attenuation_datafile;
partitions = options.partitions;
verbose = options.verbose;
empty_weight = false;
options.MBSREM_prepass = true;
% options.U = [];
% options.weights = [];
% options.a_L = [];
% options.fmh_weights = [];
% options.weighted_weights = [];
% options.weights_quad = [];
options.MAP = (options.OSL_MLEM || options.OSL_OSEM || options.BSREM || options.MBSREM || options.ROSEM_MAP || options.RBI_MAP || any(options.COSEM_MAP));

options.N = Nx * Ny * Nz;

options.rekot = reko_maker(options);
pz = cell(length(options.rekot) + 7,options.partitions);

temp = pseudot;
if ~isempty(temp) && temp > 0
    for kk = int32(1) : temp
        pseudot(kk) = int32(options.cryst_per_block + 1) * kk;
    end
end
options.pseudot = pseudot;

if (options.quad || options.FMH || options.L || options.weighted_mean || options.MRP) && options.MAP
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

options.im_vectors = form_image_vectors(options, options.N);
if options.reconstruction_method ~= 2
    if options.OSL_OSEM
        options.im_vectors.custom_OSEM = ones(options.N,options.Niter + 1);
        options.im_vectors.custom_OSEM(:,1) = options.x0(:);
    end
    if options.OSL_MLEM
        options.im_vectors.custom_MLEM = ones(options.N,options.Niter + 1);
        options.im_vectors.custom_MLEM(:,1) = options.x0(:);
    end
    if options.MBSREM
        options.im_vectors.custom_MBSREM = ones(options.N,options.Niter + 1);
        options.im_vectors.custom_MBSREM(:,1) = options.x0(:);
    end
    if options.BSREM
        options.im_vectors.custom_BSREM = ones(options.N,options.Niter + 1);
        options.im_vectors.custom_BSREM(:,1) = options.x0(:);
    end
    if options.ROSEM_MAP
        options.im_vectors.custom_ROSEM = ones(options.N,options.Niter + 1);
        options.im_vectors.custom_ROSEM(:,1) = options.x0(:);
    end
    if options.RBI_MAP
        options.im_vectors.custom_RBI = ones(options.N,options.Niter + 1);
        options.im_vectors.custom_RBI(:,1) = options.x0(:);
    end
    if any(options.COSEM_MAP)
        options.im_vectors.custom_COSEM = ones(options.N,options.Niter + 1);
        options.im_vectors.custom_COSEM(:,1) = options.x0(:);
    end
    
else
    
    if options.OSL_OSEM
        options.im_vectors.custom_OSEM = ones(options.N,options.Niter + 1,'single');
        options.im_vectors.custom_OSEM(:,1) = options.x0(:);
    end
    if options.OSL_MLEM
        options.im_vectors.custom_MLEM = ones(options.N,options.Niter + 1,'single');
        options.im_vectors.custom_MLEM(:,1) = options.x0(:);
    end
    if options.MBSREM
        options.im_vectors.custom_MBSREM = ones(options.N,options.Niter + 1,'single');
        options.im_vectors.custom_MBSREM(:,1) = options.x0(:);
    end
    if options.BSREM
        options.im_vectors.custom_BSREM = ones(options.N,options.Niter + 1,'single');
        options.im_vectors.custom_BSREM(:,1) = options.x0(:);
    end
    if options.ROSEM_MAP
        options.im_vectors.custom_ROSEM = ones(options.N,options.Niter + 1,'single');
        options.im_vectors.custom_ROSEM(:,1) = options.x0(:);
    end
    if options.RBI_MAP
        options.im_vectors.custom_RBI = ones(options.N,options.Niter + 1,'single');
        options.im_vectors.custom_RBI(:,1) = options.x0(:);
    end
    if any(options.COSEM_MAP)
        options.im_vectors.custom_COSEM = ones(options.N,options.Niter + 1,'single');
        options.im_vectors.custom_COSEM(:,1) = options.x0(:);
    end
end

if options.use_raw_data
    if options.reconstruct_trues == false && ~options.reconstruct_scatter && (isfield(options, 'coincidences') == 0 && ~exist('coincidences','var') || options.precompute_all || options.use_machine > 0)
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
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_raw_ASCII.mat'], 'coincidences')
            elseif options.use_LMF && options.use_machine == 0
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_raw_LMF.mat'], 'coincidences')
            elseif options.use_root && options.use_machine == 0
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_raw_root.mat'], 'coincidences')
            else
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_raw_listmode.mat'], 'coincidences')
            end
        end
        options.SinM = coincidences;
    elseif ~options.reconstruct_trues && ~options.reconstruct_scatter && ~exist('coincidences','var')
        options.SinM = options.coincidences;
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
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_raw_ASCII.mat'], 'true_coincidences')
            elseif options.use_LMF
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_raw_LMF.mat'], 'true_coincidences')
            elseif options.use_root
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_raw_root.mat'], 'true_coincidences')
            else
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_raw_real.mat'], 'true_coincidences')
            end
        end
        options.SinM = true_coincidences;
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
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_raw_ASCII.mat'], 'scattered_coincidences')
            elseif options.use_LMF
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_raw_LMF.mat'], 'scattered_coincidences')
            elseif options.use_root
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_raw_root.mat'], 'scattered_coincidences')
            else
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_raw_real.mat'], 'scattered_coincidences')
            end
        end
        options.SinM = scattered_coincidences;
    end
    if options.randoms_correction && ~options.reconstruct_trues && ~options.reconstruct_scatter
        if options.use_ASCII || options.use_LMF || options.use_root || options.use_machine > 0
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
                    load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_raw_ASCII.mat'], 'delayed_coincidences')
                elseif options.use_LMF && options.use_machine == 0
                    load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_raw_LMF.mat'], 'delayed_coincidences')
                elseif options.use_root && options.use_machine == 0
                    load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_raw_root.mat'], 'delayed_coincidences')
                else
                    load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_raw_listmode.mat'], 'delayed_coincidences')
                end
            end
            if exist('delayed_coincidences','var')
                if iscell(options.SinM) && iscell(delayed_coincidences)
                    for kk = 1 : length(options.SinM)
                        options.SinM{kk} = options.SinM{kk} - delayed_coincidences{kk};
                        options.SinM{kk}(options.SinM{kk} < 0) = 0;
                    end
                else
                    options.SinM = options.SinM - delayed_coincidences;
                    options.SinM(options.SinM < 0) = 0;
                end
            else
                disp('Delayed coincidences not found, randoms correction not performed')
            end
        else
            if iscell(options.SinD)
                if iscell(options.SinM)
                    for kk = 1 : length(options.SinM)
                        if numel(options.SinD{kk}) ~= numel(options.SinM{kk})
                            error('Size mismatch between randoms correction data and measurement data')
                        end
                        options.SinM{kk} = options.SinM{kk} - double(options.SinD{kk}(:));
                        options.SinM{kk}(options.SinM{kk} < 0) = 0;
                    end
                else
                    if numel(options.SinD{1}) ~= numel(options.SinM)
                        error('Size mismatch between randoms correction data and measurement data')
                    end
                    options.SinM = options.SinM - double(options.SinD{1}(:));
                    options.SinM(options.SinM < 0) = 0;
                end
            else
                if iscell(options.SinM)
                    for kk = 1 : length(options.SinM)
                        if numel(options.SinD) ~= numel(options.SinM{kk})
                            error('Size mismatch between randoms correction data and measurement data')
                        end
                        options.SinM{kk} = options.SinM{kk} - double(options.SinD(:));
                        options.SinM{kk}(options.SinM{kk} < 0) = 0;
                    end
                else
                    if numel(options.SinD) ~= numel(options.SinM)
                        error('Size mismatch between randoms correction data and measurement data')
                    end
                    options.SinM = options.SinM - double(options.SinD(:));
                    options.SinM(options.SinM < 0) = 0;
                end
            end
        end
    end
    if options.scatter_correction
        if iscell(options.ScatterC)
            if iscell(options.SinM)
                for kk = 1 : length(options.SinM)
                    if numel(options.ScatterC{kk}) ~= numel(options.SinM{kk})
                        error('Size mismatch between scatter correction data and measurement data')
                    end
                    options.SinM{kk} = options.SinM{kk} - double(options.ScatterC{kk}(:));
                    options.SinM{kk}(options.SinM{kk} < 0) = 0;
                end
            else
                if numel(options.ScatterC{1}) ~= numel(options.SinM)
                    error('Size mismatch between scatter correction data and measurement data')
                end
                options.SinM = options.SinM - double(options.ScatterC{1}(:));
                options.SinM(options.SinM < 0) = 0;
            end
        else
            if iscell(options.SinM)
                for kk = 1 : length(options.SinM)
                    if numel(options.ScatterC) ~= numel(options.SinM{kk})
                        error('Size mismatch between scatter correction data and measurement data')
                    end
                    options.SinM{kk} = options.SinM{kk} - double(options.ScatterC(:));
                    options.SinM{kk}(options.SinM{kk} < 0) = 0;
                end
            else
                if numel(options.ScatterC) ~= numel(options.SinM)
                    error('Size mismatch between scatter correction data and measurement data')
                end
                options.SinM = options.SinM - double(options.ScatterC(:));
                options.SinM(options.SinM < 0) = 0;
            end
        end
    end
    
    clear coincidences options.coincidences true_coincidences delayed_coincidences
else
    if (~options.reconstruct_trues && ~options.reconstruct_scatter) || options.use_machine > 0
        if options.partitions == 1 && (isfield(options, 'SinM') == 0 || options.precompute_all)
            if options.use_machine == 0
                variableInfo = who('-file', [options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                    num2str(options.span) '.mat']);
            elseif  options.use_machine == 1
                variableInfo = who('-file', [options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                    num2str(options.span) '_listmode.mat']);
            else
                variableInfo = who('-file', [options.machine_name '_' options.name '_sinogram_original_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                    num2str(options.span) '_machine_sinogram.mat']);
            end
            if options.randoms_correction || options.scatter_correction || (options.fill_sinogram_gaps && sum(options.pseudot) > 0) || ~ismember('raw_SinM', variableInfo)
                if options.use_machine == 0
                    load([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span) '.mat'],'SinM')
                elseif  options.use_machine == 1
                    load([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span) '_listmode.mat'],'SinM')
                else
                    load([options.machine_name '_' options.name '_sinogram_original_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span) '_machine_sinogram.mat'],'SinM')
                end
                options.SinM = SinM;
                clear SinM
            else
                if options.use_machine == 0
                    load([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span) '.mat'],'raw_SinM')
                elseif  options.use_machine == 1
                    load([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span) '_listmode.mat'],'raw_SinM')
                else
                    load([options.machine_name '_' options.name '_sinogram_original_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span) '_machine_sinogram.mat'],'raw_SinM')
                end
                options.SinM = raw_SinM;
                clear raw_SinM
            end
        elseif isfield(options, 'SinM') == 0 || options.precompute_all
            variableInfo = who('-file', [options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_' ...
                num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span) '.mat']);
            if options.randoms_correction || options.scatter_correction || (options.fill_sinogram_gaps && sum(options.pseudot) > 0) || ~ismember('raw_SinM', variableInfo)
                if options.use_machine == 0
                    load([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_' num2str(options.Ndist) ...
                        'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span) '.mat'], 'SinM')
                elseif  options.use_machine == 1
                    load([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_' num2str(options.Ndist) ...
                    'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span) '_listmode.mat'], 'SinM')
                else
                    load([options.machine_name '_' options.name '_sinograms_original_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_' num2str(options.Ndist) ...
                    'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span) '_machine_sinogram.mat'], 'SinM')
                end
                options.SinM = SinM;
                clear SinM
            else
                if options.use_machine == 0
                    load([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_' num2str(options.Ndist) ...
                        'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span) '.mat'], 'raw_SinM')
                elseif  options.use_machine == 1
                    load([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_' num2str(options.Ndist) ...
                    'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span) '_listmode.mat'], 'raw_SinM')
                else
                    load([options.machine_name '_' options.name '_sinograms_original_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_' num2str(options.Ndist) ...
                    'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span) '_machine_sinogram.mat'], 'raw_SinM')
                end
                options.SinM = raw_SinM;
                clear raw_SinM
            end
        else
            options.SinM = options.SinM;
            clear options.SinM
        end
    elseif options.reconstruct_trues && options.use_machine == 0
        if options.partitions == 1 && isfield(options, 'SinTrues') == 0
            load([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span) '.mat'],'SinTrues')
        elseif isfield(options, 'SinTrues') == 0
            load([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' ...
                num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span) '.mat'], 'SinTrues')
        end
        options.SinM = SinTrues;
        clear SinTrues
    elseif options.reconstruct_scatter && options.use_machine == 0
        if options.partitions == 1 && isfield(options, 'SinScatter') == 0
            load([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span) '.mat'],'SinScatter')
        elseif isfield(options, 'SinScatter') == 0
            load([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' ...
                num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' num2str(options.span) '.mat'], 'SinScatter')
        end
        options.SinM = SinScatter;
        clear SinScatter
    end
end


%%
if iscell(options.SinM)
    Sino = options.SinM{1};
else
    Sino = options.SinM;
end

Sino = Sino(:);

if issparse(Sino)
    Sino = (full(Sino));
end

if use_raw_data == false && NSinos ~= TotSinos
    Sino = Sino(1:NSinos*Nang*Ndist);
end

%% This computes a whole observation matrix and uses it to compute the MLEM (no on-the-fly calculations)
%% This part is used when the observation matrix is calculated on-the-fly

% Compute the indices for the subsets used.
% For Sinogram data, five different methods to select the subsets are
% available. For raw list-mode data, three methods are available.
[index, options.pituus, options.subsets, options.lor_a] = index_maker(Nx, Ny, Nz, options.subsets, use_raw_data, Nang, Ndist, TotSinos, NSinos, machine_name, options);



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
options.blocks=int32(rings + length(pseudot) - 1);
% First ring
options.block1=int32(0);

options.NSinos = int32(NSinos);
options.NSlices = int32(Nz);
options.TotSinos = int32(TotSinos);

if use_raw_data == false
    load([machine_name '_app_coordinates_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'x','y');
    load([machine_name '_3D_coordinates_span' num2str(options.span) '_ringdiff' num2str(options.ring_difference) '_' num2str(TotSinos) '.mat'],'z')
    if NSinos ~= TotSinos
        z = z(1:NSinos,:);
    end
else
    load([machine_name '_detector_coordinates.mat'],'x','y');
    if ~exist('LL','var')
        if ~exist('discard','var')
            load([machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'],'lor','discard')
        end
        options.LL = form_detector_pairs_raw(rings, det_per_ring);
        options.LL = options.LL(discard,:);
    end
    
    options.z_length = double(options.blocks) * cr_pz;
    z = linspace(0, z_length, options.blocks + 1);
    
    
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
    options.vaimennus = vaimennus(:);
    if options.reconstruction_method == 2 || options.reconstruction_method == 4
        options.vaimennus = single(options.vaimennus);
    end
    clear data
else
    if options.reconstruction_method == 2 || options.reconstruction_method == 4
        options.vaimennus = single(0);
    else
        options.vaimennus = 0;
    end
end

if min(min(z)) == 0
    z = z + (axial_fow - max(max(z)))/2;
end

%     x = circshift(x,length(x)/2);
%     y = circshift(y,length(y)/2);

if options.reconstruction_method == 2 || options.reconstruction_method == 4
    options.x=single(x);
    options.y=single(y);
    options.z_det = single(z);
else
    options.x=double(x);
    options.y=double(y);
    options.z_det = double(z);
end
clear z


options.size_x = int32(size(x,1));

if options.precompute_lor && options.subsets > 1 || options.reconstruction_method == 2  && options.subsets > 1 || options.reconstruction_method == 4 && options.subsets > 1
    options.pituus = [0;cumsum(options.pituus)];
    if iscell(index)
        index = cell2mat(index);
    end
end

if use_raw_data
    if isempty(options.pseudot)
        options.pseudot = int32(1e5);
    else
        options.pseudot = options.pseudot - 1;
    end
end

% for the precomputed version, index vectors are needed
if use_raw_data == false && (options.precompute_lor || options.reconstruction_method == 2 || options.reconstruction_method == 4)
    
    
    
    if subsets > 1
        if options.mramla || options.MBSREM || options.rbi || options.RBI_MAP || options.cosem || options.ecosem...
                || options.acosem || any(options.COSEM_MAP)
            if options.use_raw_data
                f_name = [machine_name '_MBSREM_precompute_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'];
            else
                f_name = [machine_name '_MBSREM_precompute_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_sino_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'];
            end
            if exist(f_name,'file') == 2
                load(f_name,'D','Amin');
                options.MBSREM_prepass = false;
                options.Amin = Amin(index);
                options.D = D;
                if options.reconstruction_method == 2
                    options.D = single(D);
                    if options.mramla || options.MBSREM
                        options.Amin = single(Amin);
                    end
                end
                if ~options.mramla && ~options.MBSREM
                    clear Amin
                end
                clear D Amin
            end
        end
        options.lor_a = (options.lor_a(index));
        Sino = Sino(index);
        if partitions > 1
            for ff = 1 : partitions
                temp = options.SinM{ff};
                if NSinos ~= TotSinos
                    temp = temp(:,:,1:NSinos);
                end
                temp = temp(index);
                options.SinM{ff} = temp;
            end
            clear temp
        else
            temp = options.SinM;
            if NSinos ~= TotSinos
                temp = temp(:,:,1:NSinos);
            end
            temp = temp(index);
            options.SinM = temp;
        end
    else
        load([machine_name '_lor_pixel_count_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_sino_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'lor','discard')
        if length(discard) ~= TotSinos*Nang*Ndist
            error('Error: Size mismatch between sinogram and LORs to be removed')
        end
        if use_raw_data == false && NSinos ~= TotSinos
            discard = discard(1:NSinos*Nang*Ndist);
        end
        options.lor_a = (lor(discard));
        Sino = Sino(discard);
        if partitions > 1
            for ff = 1 : partitions
                temp = options.SinM{ff};
                if NSinos ~= TotSinos
                    temp = temp(:,:,1:NSinos);
                end
                temp = temp(discard);
                options.SinM{ff} = temp;
            end
            clear temp
        else
            temp = options.SinM;
            if NSinos ~= TotSinos
                temp = temp(:,:,1:NSinos);
            end
            temp = temp(discard);
            options.SinM = temp;
        end
        clear lor
    end
    [~, I] = sort(y, 2);
    sy = size(y);
    I = sub2ind(sy, repmat((1:sy(1)).', 1, sy(2)), I);
    
    xy_index = uint32(I(:,1));
    xy_index2 = repmat(uint32(1:options.size_x)', NSinos - Nz, 1);
    xy_index = [repmat(xy_index, Nz, 1); xy_index2];
    [~, I] = sort(options.z_det, 2);
    sy = size(options.z_det);
    I = sub2ind(sy, repmat((1:sy(1)).', 1, sy(2)), I);
    
    z_index = uint32(I(:,1));
    z_index = repelem(z_index, options.size_x);
    if options.subsets > 1
        z_index = z_index(index);
    else
        z_index = (z_index(discard));
    end
    apu = z_index > NSinos;
    options.z_index = z_index - 1;
    
    if options.subsets > 1
        xy_index = xy_index(index);
    else
        xy_index = (xy_index(discard));
    end
    xy_index(apu) = xy_index(apu) + uint32(options.size_x);
    options.xy_index = xy_index - 1;
    
    options.summa = zeros(options.subsets, 1, 'uint64');
    
    if options.subsets > 1
        for kk = 1 : options.subsets
            options.summa(kk) = sum(int64(options.lor_a(options.pituus(kk)+1:options.pituus(kk+1))));
            %                 koko(kk) = -(pituus(kk))+pituus(kk+1);
        end
        %             vector_size = max(koko);
    else
        options.summa = uint64(sum(int64(options.lor_a)));
        options.pituus = int32([0;length(Sino)]);
    end
    
    
    clear discard I yt xt xy_index2 index apu xy_index z_index
elseif use_raw_data && (options.precompute_lor || options.reconstruction_method == 2  || options.reconstruction_method == 4)
    
    if options.subsets > 1
        if options.mramla || options.MBSREM || options.rbi || options.RBI_MAP || options.cosem || options.ecosem...
                || options.acosem || any(options.COSEM_MAP)
            if options.use_raw_data
                f_name = [machine_name '_MBSREM_precompute_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'];
            else
                f_name = [machine_name '_MBSREM_precompute_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_sino_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'];
            end
            if exist(f_name,'file') == 2
                load(f_name,'D','Amin');
                options.MBSREM_prepass = false;
                options.Amin = Amin(index);
                options.D = D;
                if options.reconstruction_method == 2
                    options.D = single(D);
                    if options.mramla || options.MBSREM
                        options.Amin = single(Amin);
                    end
                end
                if ~options.mramla && ~options.MBSREM
                    clear Amin
                end
                clear D Amin
            end
        end
        options.LL = options.LL(index,:);
        options.lor_a = (options.lor_a(index));
        Sino = Sino(index);
        if partitions > 1
            for ff = 1 : partitions
                temp = options.SinM{ff};
                if NSinos ~= TotSinos
                    temp = temp(:,:,1:NSinos);
                end
                temp = temp(index);
                options.SinM{ff} = temp;
            end
            clear temp
        end
    else
        if ~exist('LL','var')
            if ~exist('discard','var')
                load([machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'],'lor','discard')
            end
            options.LL = form_detector_pairs_raw(rings, det_per_ring);
            options.LL = options.LL(discard,:);
        end
        options.pituus = int32([0;length(Sino)]);
        options.lor_a = lor;
        clear lor discard
    end
    options.summa = zeros(options.subsets, 1, 'uint64');
    %             y = circshift(y, -length(y)/4);
    %             x = circshift(x, -length(x)/4);
    
    for kk = 1 : options.subsets
        apu = options.LL(options.pituus(kk) + 1 : options.pituus(kk + 1),:) - 1;
        apu2 = idivide(apu, uint16(det_per_ring));
        idx = apu2(:,1) == apu2(:,2);
        apu2 = apu(idx,:);
        ind = mod(apu2, uint16(det_per_ring)) + 1;
        yt = y(ind);
        y_i = yt(:,1) > yt(:,2);
        apu2(y_i,:) = fliplr(apu2(y_i,:));
        apu(idx,:) = apu2;
        options.LL(options.pituus(kk) + 1 : options.pituus(kk + 1),:) = apu + 1;
        options.summa(kk) = sum(int64(options.lor_a(options.pituus(kk)+1:options.pituus(kk+1))));
        %             koko(kk) = -(pituus(kk))+pituus(kk+1);
    end
    
    clear apu apu2 idx ind yt y_i index discard
    
    options.LL = options.LL';
    options.LL = options.LL(:);
end


% Pixels
etaisyys_x=(R-FOVax)/2;
etaisyys_y=(R-FOVay)/2;
if options.reconstruction_method == 2 || options.reconstruction_method == 4
    options.zz=linspace(single(0),single(axial_fow),Nz+1);
    options.xx = single(linspace(etaisyys_x,R-etaisyys_x,Nx+1));
    options.yy = single(linspace(etaisyys_y,R-etaisyys_y,Ny+1));
else
    options.zz=linspace(double(0),double(axial_fow),Nz+1);
    options.xx = double(linspace(etaisyys_x,R-etaisyys_x,Nx+1));
    options.yy = double(linspace(etaisyys_y,R-etaisyys_y,Ny+1));
end
options.zz=options.zz(2*options.block1+1:2*options.blocks+2);

% Distance of adjacent pixels
options.dx=diff(options.xx(1:2));
options.dy=diff(options.yy(1:2));
options.dz=diff(options.zz(1:2));

% Distance of image from the origin
options.bx=options.xx(1);
options.by=options.yy(1);
options.bz=options.zz(1);

% Number of pixels
options.Ny=int32(Ny);
options.Nx=int32(Nx);
options.Nz=int32(Nz);

%     if options.reconstruction_method == 2 || options.reconstruction_method == 4
%         iij=single(0:Nx);
%         jji=single(0:Ny);
%         kkj=single(0:Nz);
%     else
%         iij=double(0:Nx);
%         jji=double(0:Ny);
%         kkj=double(0:Nz);
%     end

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
    if options.reconstruction_method == 2 || options.reconstruction_method == 4
        options.zmax = single(1);
    else
        options.zmax = double(1);
    end
end
%%

if (options.MRP || options.quad || options.TV ||options. FMH || options.L || options.weighted_mean || options.APLS || options.BSREM || options.ramla || options.MBSREM || options.mramla ...
        || options.rosem || options.drama || options.ROSEM_MAP || options.ecosem || options.cosem || options.acosem || options.AD || any(options.COSEM_MAP)...
        || (options.NLM && options.NLM_use_anatomical))
    
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
                    options.TVdata.s1 = single(S(1:3:end,1));
                    options.TVdata.s2 = single(S(1:3:end,2));
                    options.TVdata.s3 = single(S(1:3:end,3));
                    options.TVdata.s4 = single(S(2:3:end,1));
                    options.TVdata.s5 = single(S(2:3:end,2));
                    options.TVdata.s6 = single(S(2:3:end,3));
                    options.TVdata.s7 = single(S(3:3:end,1));
                    options.TVdata.s8 = single(S(3:3:end,2));
                    options.TVdata.s9 = single(S(3:3:end,3));
                else
                    options.TVdata.s1 = S(1:3:end,1);
                    options.TVdata.s2 = S(1:3:end,2);
                    options.TVdata.s3 = S(1:3:end,3);
                    options.TVdata.s4 = S(2:3:end,1);
                    options.TVdata.s5 = S(2:3:end,2);
                    options.TVdata.s6 = S(2:3:end,3);
                    options.TVdata.s7 = S(3:3:end,1);
                    options.TVdata.s8 = S(3:3:end,2);
                    options.TVdata.s9 = S(3:3:end,3);
                end
            end
            if options.reconstruction_method == 2
                options.TVdata.reference_image = single(alkuarvo);
                options.TVdata.T = single(options.T);
                options.TVdata.C = single(options.C);
            else
                options.TVdata.reference_image = alkuarvo;
                options.TVdata.T = options.T;
                options.TVdata.C = options.C;
            end
            clear apu variables alkuarvo S
        end
        if options.reconstruction_method == 2
            options.tau = single(options.tau);
            options.TVdata.beta = single(options.TVsmoothing);
%             options.TVdata = TVdata;
%             clear TVdata;
        else
            options.TVdata.beta = options.TVsmoothing;
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
    
    if ((options.mramla || options.MBSREM || options.rbi || options.RBI_MAP) && options.MBSREM_prepass || options.ecosem || options.cosem || options.acosem || any(options.COSEM_MAP))  && options.reconstruction_method == 1
        
        if options.acosem
            options.C_aco = zeros(double(N), subsets);
        end
        if options.cosem || options.ecosem
            options.C_co = zeros(double(N), subsets);
        end
        if any(options.COSEM_MAP)
            options.C_osl = zeros(double(N), subsets);
        end
        
        if ~use_raw_data
            if isempty(options.pseudot)
                options.pseudot = int32(0);
            end
        end
        
        options.D = zeros(options.N,1);
        
        if verbose
            disp('Prepass phase for MRAMLA, COSEM, ACOSEM and ECOSEM started')
        end
        for osa_iter = 1 : subsets
            if options.precompute_lor == false
                if options.use_raw_data == false
                    [ lor, indices, alkiot] = improved_Siddon_algorithm( options.verbose, options.options.options.Ny, options.options.Nx, options.options.Nz, options.dx, options.dz, options.by, ...
                        options.bx, options.bz, options.z_det, options.x, options.y, options.dy, options.yy, options.xx , options.NSinos, options.NSlices, ...
                        options.size_x, options.zmax, options.NSinos, options.ind_size, options.vaimennus, options.index{osa_iter}, options.pituus(osa_iter), ...
                        options.attenuation_correction, options.use_raw_data, uint16(0), options.pseudot, options.block1, options.blocks, options.det_per_ring);
                else
                    L = options.LL(options.index{osa_iter},:);
                    L = L';
                    L = L(:);
                    [ lor, indices, alkiot] = improved_Siddon_algorithm(options.verbose, options.options.options.Ny, options.options.Nx, options.options.Nz, options.dx, options.dz, options.by, ...
                        options.bx, options.bz, options.z_det, options.x, options.y, options.dy, options.yy, options.xx , options.NSinos, options.NSlices, ...
                        options.size_x, options.zmax, options.NSinos, options.ind_size, options.vaimennus, uint32(0), int32(0), options.attenuation_correction, ...
                        options.use_raw_data, L, options.pseudot, options.block1, options.blocks, options.det_per_ring);
                end
                lor = reshape(lor,[],2);
                lor=repelem(int32((lor(:,1))),lor(:,2));
                A_length = length(Sino(index{osa_iter}));
                indices=indices + 1;
                if verbose
                    tStart = tic;
                end
                if use_raw_data == false
                    if options.use_fsparse == false
                        A = sparse(double(lor),double(indices),double(alkiot),A_length,double(N));
                    else
                        A = fsparse(lor,indices,double(alkiot),[A_length double(N) length(alkiot)]);
                    end
                else
                    if options.use_fsparse == false
                        A = sparse(double(lor),double(indices),double(alkiot),A_length,double(N));
                    else
                        A = fsparse(lor,indices,double(alkiot),[A_length double(N) length(alkiot)]);
                    end
                end
                clear indices alkiot lor
                if verbose
                    tElapsed = toc(tStart);
                    disp(['Sparse matrix formation took ' num2str(tElapsed) ' seconds'])
                end
                A = A';
            else
                lor2 = [0; cumsum(uint32(options.lor_a(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1))))];
                if options.use_raw_data == false
                    [A] = improved_Siddon_algorithm_array( options.options.options.Ny, options.options.Nx, options.options.Nz, options.dx, options.dz, options.by, options.bx, options.bz, options.z_det, ...
                        options.x, options.y, options.dy, options.yy, options.xx , options.NSinos, options.NSlices, options.size_x, options.zmax, options.vaimennus, ...
                        lor2, options.pituus(osa_iter + 1) - options.pituus(osa_iter), options.attenuation_correction, options.lor_a(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1)),...
                        options.summa(osa_iter), options.xy_index(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1)), options.z_index(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1)), options.NSinos, ...
                        uint16(0), options.pseudot, options.det_per_ring, options.verbose, options.use_raw_data);
                else
                    [A] = improved_Siddon_algorithm_array(options.options.options.Ny, options.options.Nx, options.options.Nz, options.dx, options.dz, options.by, options.bx, options.bz, options.z_det, ...
                        options.x, options.y, options.dy, options.yy, options.xx , options.NSinos, options.NSlices, options.size_x, options.zmax, options.vaimennus, ...
                        lor2, options.pituus(osa_iter + 1) - options.pituus(osa_iter), options.attenuation_correction, options.lor_a(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1)),...
                        options.summa(osa_iter), uint32(0), uint32(0), options.NSinos, options.LL(options.pituus(osa_iter) * 2 + 1 : options.pituus(osa_iter + 1) * 2), ...
                        options.pseudot, options.det_per_ring, options.verbose, options.use_raw_data);
                end
                clear lor2
            end
            options.D = options.D + A*ones(size(A,2),1,'double');
            %                 if subsets > 2000
            %                     pj(:,osa_iter) = single(A*ones(size(A,2),1,'double'));
            %                 else
            %                     pj(:,osa_iter) = A*ones(size(A,2),1,'double');
            %                 end
            if options.ecosem || options.cosem || options.acosem || any(options.COSEM_MAP)
                if options.precompute_lor == false
                    uu = double(Sino(index{osa_iter}));
                else
                    uu = double(Sino(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                end
            end
            if options.cosem || options.ecosem
                % C_co(:,osa_iter) = full(sum(spdiags(uu./(A'*options.x0(:)+epps),0,size(A,2),size(A,2))*(bsxfun(@times,A,options.x0(:))))');
                if osa_iter > 1
                    options.C_co(:,osa_iter) = full(sum(spdiags(options.x0(:),0,size(A,1),size(A,1))*A*spdiags(uu./(A'*options.x0(:)+epps),0,size(A,2),size(A,2)),2));
                    %                         B = A;
                    %                         B(B > 0) = 1;
                    %                         C_co(:,osa_iter) = full(sum(B*spdiags(uu,0,size(A,2),size(A,2)),2));
                end
            end
            if options.acosem
                % C_aco(:,osa_iter) = full(sum(spdiags(uu./(A'*options.x0(:)+epps),0,size(A,2),size(A,2))*(bsxfun(@times,A',(options.x0(:)')).^(1/options.h)))');
                if osa_iter > 1
                    options.C_aco(:,osa_iter) = full(sum(spdiags(power(options.x0(:),1/options.h),0,size(A,1),size(A,1))*A*spdiags(uu./(A'*options.x0(:)+epps),0,size(A,2),size(A,2)),2));
                end
            end
            if any(options.COSEM_MAP)
                if options.COSEM_MAP == 2
                    if osa_iter > 1
                        options.C_osl(:,osa_iter) = full(sum(spdiags(options.x0(:),0,size(A,1),size(A,1))*A*spdiags(uu./(A'*options.x0(:)+epps),0,size(A,2),size(A,2)),2));
                    end
                else
                    if osa_iter > 1
                        options.C_osl(:,osa_iter) = full(sum(spdiags(power(options.x0(:),1/options.h),0,size(A,1),size(A,1))*A*spdiags(uu./(A'*options.x0(:)+epps),0,size(A,2),size(A,2)),2));
                    end
                end
            end
            if options.MBSREM_prepass && options.U == 0 && (options.MBSREM || options.mramla)
                %%%% This particular piece of code was taken from:
                %%%% https://se.mathworks.com/matlabcentral/answers/35309-max-min-of-sparse-matrices
                [~,m] = size(A);
                rowMin = nan(m, 1);
                [~,I,S] = find(A);
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
                    options.Amin(index{osa_iter}) = rowMin;
                else
                    options.Amin(pituus(osa_iter)+1:pituus(osa_iter + 1)) = rowMin;
                end
                clear I K S markers rowMin s iRows
            end
            clear A
        end
        %             D = sum(pj,2);
        %             if rekot(8) || rekot(9)
        %                 B = sum(C_co,2);
        %             end
        if verbose
            disp('Prepass phase for COSEM, ACOSEM and ECOSEM completed')
        end
    end
    if (options.BSREM || options.ramla) && length(options.lambda0) == 1
        options.lam = zeros(Niter,1);
        options.lam(1) = options.lambda0;
        %             if lam(1) > 1/max(max(pj))
        %                 lam(1) = min(min(pj));
        %             end
        for i=2:Niter
            %                 lam(i) = 0.5*lam(i-1);
            options.lam(i) = options.lam(1)/i;
            %                 lam(i) = lam(1)/1.01;
        end
        if options.reconstruction_method == 2
            options.lam = single(options.lam);
        else
            options.lambda0 = options.lam;
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
        options.pj3 = options.D/options.subsets;
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
                            apu = [((Ndx:-1:-Ndx) * distX)', (repelem(kk,Ndy*2+1) * distY)'];
                        else
                            if Ndz ~= Ndx
                                apu = [((Ndx:-1:-Ndx) * distX)', (repelem(kk,Ndy*2+1) * distY)', [zeros(Ndx-Ndz,1),(repelem(jj,Ndz*2+1) * distZ),zeros(Ndx-Ndz,1)]'];
                            else
                                apu = [((Ndx:-1:-Ndx) * distX)', (repelem(kk,Ndy*2+1) * distY)', (repelem(jj,Ndz*2+1) * distZ)'];
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
            tr_ind = sub2ind([Nx+Ndx*2 Ny+Ndy*2],mod((1:N)'-1,Nx)+(Ndx + 1),mod(floor(((1:double(N))'-1)/double(Nx)),double(Ny))+(Ndy + 1));
        else
            tr_ind = sub2ind([Nx+Ndx*2 Ny+Ndy*2 Nz+Ndz*2],mod((1:options.N)'-1,Nx)+(Ndx + 1),mod(floor(((1:double(options.N))'-1)/double(Nx)),double(Ny))+(Ndy + 1),floor(((1:double(options.N))'-1)/double(Nx*Ny))+(Ndz+1));
        end
        options.tr_offsets = uint32(bsxfun(@plus,tr_ind,offsets(:)'));
        if options.reconstruction_method == 2
            options.tr_offsets = options.tr_offsets - 1;
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
                options.weights_quad = [options.weights_quad(1:floor(length(options.weights_quad)/2));options.weights_quad(ceil(length(options.weights_quad)/2) + 1:end)];
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
        if options.FMH || options.quad
            options.weights = single(options.weights);
            options.inffi = uint32(find(isinf(options.weights)) - 1);
        end
        if options.MRP
            options.medx = options.Ndx*2 + 1;
            options.medy = options.Ndy*2 + 1;
            options.medz = options.Ndz*2 + 1;
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
            options.FluxType = int32(options.FluxType);
            options.DiffusionType = int32(options.DiffusionType);
        else
            if options.FluxType == 1
                options.FluxType = 'exponential';
            elseif options.FluxType == 2
                options.FluxType = 'quadratic';
            end
        end
    end
    if options.MBSREM || options.mramla
        options.epsilon_mramla = MBSREM_epsilon(Sino, options.epps);
        if options.reconstruction_method == 2
            options.epsilon_mramla = single(options.epsilon_mramla);
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
end