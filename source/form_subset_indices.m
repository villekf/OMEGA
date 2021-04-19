function [options, lor_a, xy_index, z_index, LL, summa, pituus, varargout] = form_subset_indices(options, pituus, subsets, index, size_x, y, z_det, rings, fpbp, varargin)
%% FORM SUBSET INDICES
% This function is used to form the necessary index variables for the
% reconstruction depending on the number of subsets and whether sinogram or
% raw list-mode data is used.
% INPUTS:
%   options = Contains the machine and (optionally) sinogram specific
%   information (e.g. detectors per ring, image (matrix) size, FOV size).
%   pituus = The number of measurements (LORs) that each subset has.
%   subsets = The number of subsets used.
%   index = The indices corresponding to each subset.
%   size_x = Length of the detector coordinate matrix.
%   y = Detector coordinates in the y-axis.
%   z_det = Same as above, but for z-axis.
%   fbbp = Is the forward/backward projection code used.
%   TOF = Is TOF data used.
%   measurement_data = The measurement data (sinogram or raw list-mode
%   data). This input is optional and is only used if the measurement data
%   is also output.
%
% OUTPUTS:
%   options = Contains the machine and (optionally) sinogram specific
%   information (e.g. detectors per ring, image (matrix) size, FOV size).
%   lor_a = If options.precompute_lor = true, then this vector contains the
%   measurement rows (LORs) that are included in the reconstruction at each
%   subset. E.g. LORs outside the FOV are removed with it. If the value is
%   set to false, returns an empty array.
%   xy_index = The indices for the detector elements for each subset in the
%   x/y-axis. Each row gives the index of the detector element coordinate
%   for the corresponding measurement (LOR). Used only for sinogram data.
%   If raw data is used, returns an empty array.
%   z_index = Same as above, but for z-axis.
%   LL = The order of the detectors in the reconstruction for each subset.
%   Used only for raw data. If sinogram data is used, returns an empty
%   array.
%   summa = Total number of voxels intersected at each subset. This is
%   needed for implementation 1 in order to preallocate sufficiently large
%   output matrix. Only applicable if options.precompute_lor = true,
%   otherwise returns an array full of zeros.
%   pituus = The number of measurements (LORs) that each subset has.
%   measurement_data = The measurement data (sinogram or raw list-mode
%   data) sorted according to the subset indices. This output is optional.
%
% See also index_maker


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ville-Veikko Wettenhovi
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

folder = fileparts(which('form_subset_indices.m'));
folder = [folder(1:end-6), 'mat-files/'];
folder = strrep(folder, '\','/');

if ~isempty(varargin) && ~isempty(varargin{1})
    TOF = varargin{1};
else
    TOF = false;
end
if nargout > 7
    if nargin > 10 && ~isempty(varargin{2})
        varargout{1} = varargin{2};
    end
end
if nargin >= 12 && ~isempty(varargin{3}) && options.precompute_lor
    lor = varargin{3};
    if options.projector_type == 3 && ~isempty(varargin{4})
        lor_orth = varargin{4};
    end
end
if nargin >= 14 && ~isempty(varargin{5}) && options.precompute_lor
    storeMatrix = varargin{5};
else
    storeMatrix = false;
end
if nargout > 8
    varargout{2} = [];
end
if nargout > 9
    varargout{3} = false;
end
if options.sampling > 1 && ~options.use_raw_data && ~options.precompute_lor
    options.Ndist = options.Ndist * options.sampling;
end
if ~isfield(options,'sampling_raw')
    options.sampling_raw = 1;
end
if ~isfield(options,'CT')
    options.CT = false;
end
% for the precomputed version, index vectors are needed
if options.use_raw_data == false && options.precompute_lor
    
    lor_file = [folder options.machine_name '_lor_pixel_count_' num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '_sino_' num2str(options.Ndist) 'x' ...
        num2str(options.Nang) 'x' num2str(options.TotSinos) '.mat'];
    
    if exist(lor_file, 'file') == 2 && exist('lor','var') ~= 1 && ~options.CT
        if options.implementation == 1 || options.implementation == 4
            variableInfo = who('-file', lor_file);
            if any(ismember('lor', variableInfo))
                load(lor_file,'lor')
            else
                lor_pixel_count_prepass(options);
                load(lor_file,'lor')
            end
            if options.projector_type == 2 && options.implementation == 1
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
                    if length(lor_orth) > options.Nang*options.Ndist*options.NSinos
                        lor_orth = lor_orth(1:options.Nang*options.Ndist*options.NSinos);
                    end
                elseif options.tube_width_z > 0
                    if crystal_size_z ~= options.tube_width_z
                        lor_pixel_count_prepass(options);
                        load(lor_file,'lor_orth')
                    end
                    if length(lor_orth) > options.Nang*options.Ndist*options.TotSinos
                        lor_orth = lor_orth(options.Nang*options.Ndist*options.TotSinos+1:end);
                    elseif length(lor_orth) == options.Nang*options.Ndist*options.TotSinos
                        lor_pixel_count_prepass(options);
                        load(lor_file,'lor_orth')
                        lor_orth = lor_orth(options.Nang*options.Ndist*options.TotSinos+1:end);
                    end
                end
            end
            if options.projector_type == 3 && options.implementation == 1
                if any(ismember('lor_vol', variableInfo))
                    load(lor_file,'lor_vol')
                else
                    lor_pixel_count_prepass(options);
                    load(lor_file,'lor_vol')
                end
                lor_orth = lor_vol;
                clear lor_vol
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
    elseif exist('lor','var') ~= 1 && ~options.CT
        lor_pixel_count_prepass(options);
        if options.implementation == 1 || options.implementation == 4
            load(lor_file,'lor')
            if options.projector_type == 2 && options.implementation == 1
                load(lor_file,'lor_orth')
                if options.tube_width_z == 0
                    if length(lor_orth) > options.Nang*options.Ndist*options.NSinos
                        lor_orth = lor_orth(1:options.Nang*options.Ndist*options.NSinos);
                    end
                elseif options.tube_width_z > 0
                    if length(lor_orth) > options.Nang*options.Ndist*options.NSinos
                        lor_orth = lor_orth(options.Nang*options.Ndist*options.NSinos+1:end);
                    end
                end
            end
            if options.projector_type == 3 && options.implementation == 1
                load(lor_file,'lor_vol')
                lor_orth = lor_vol;
                clear lor_vol
                if length(lor_orth) > options.Nang*options.Ndist*options.NSinos
                    lor_orth = lor_orth(options.Nang*options.Ndist*options.NSinos+1:end);
                end
            end
        else
            load(lor_file,'lor_opencl')
            lor = lor_opencl;
            clear lor_opencl
        end
    end
    if (exist('lor','var') ~= 1 || isempty(lor)) && options.CT
        [lor, ~, ~, lor_orth] = lor_pixel_count_prepass(options, false);
    end
    if (subsets > 1 || fpbp)
        lor_a = (lor(index));
        if (options.projector_type == 2 || options.projector_type == 3) && options.implementation == 1
            lor_orth = (lor_orth(index));
        end
        if options.normalization_correction && options.corrections_during_reconstruction
            options.normalization = options.normalization(index);
        end
        if (options.randoms_correction || (options.scatter_correction  && options.subtract_scatter)) && options.corrections_during_reconstruction ...
                && ~options.reconstruct_trues && ~options.reconstruct_scatter
            if options.partitions > 1
                for ff = 1 : options.partitions
                    temp = options.SinDelayed{ff};
                    if options.NSinos ~= options.TotSinos
                        temp = temp(:,:,1:options.NSinos);
                    end
                    temp = temp(index);
                    options.SinDelayed{ff} = temp;
                end
                clear temp
            else
                if iscell(options.SinDelayed)
                    options.SinDelayed{1} = options.SinDelayed{1}(index);
                else
                    options.SinDelayed = options.SinDelayed(index);
                end
            end
        end
        if options.scatter_correction && options.corrections_during_reconstruction && ~options.subtract_scatter ...
                && ~options.reconstruct_trues && ~options.reconstruct_scatter
            if options.partitions > 1 && iscell(options.ScatterC) && length(options.ScatterC) > 1
                for ff = 1 : options.partitions
                    temp = options.ScatterC{ff};
                    if options.NSinos ~= options.TotSinos
                        temp = temp(:,:,1:options.NSinos);
                    end
                    temp = temp(index);
                    options.ScatterC{ff} = temp;
                end
                clear temp
            else
                if iscell(options.ScatterC)
                    options.ScatterC = options.ScatterC{1}(index);
                else
                    options.ScatterC = options.ScatterC(index);
                end
            end
        end
        if nargout >= 8 && ~isempty(varargin) && length(varargin) > 1 && ~isempty(varargin{2})
            if options.partitions > 1
                for ff = 1 : options.partitions
                    temp = varargin{2}{ff};
                    if options.NSinos ~= options.TotSinos
                        temp = temp(:,:,1:options.NSinos,:);
                    end
                    if TOF
                        temp = reshape(temp, numel(temp) / options.TOF_bins, options.TOF_bins);
                        temp = temp(index,:);
                    else
                        temp = temp(index);
                    end
                    varargout{1}{ff} = temp(:);
                end
                clear temp
            else
                if options.NSinos ~= options.TotSinos
                    varargin{2} = varargin{2}(:,:,1:options.NSinos,:);
                end
                if TOF
                    varargout{1} = reshape(varargin{2}, numel(varargin{2}) / options.TOF_bins, options.TOF_bins);
                    varargout{1} = varargout{1}(index,:);
                else
                    varargout{1} = varargin{2}(index);
                end
                varargout{1} = varargout{1}(:);
            end
        end
        clear lor
    else
        if storeMatrix
            discard = true(size(lor));
        else
            discard = lor > 0;
        end
        if length(discard) ~= options.TotSinos*options.Ndist*options.Nang
            error('Error: Size mismatch between sinogram and LORs to be removed')
        end
        if options.use_raw_data == false && options.NSinos ~= options.TotSinos
            discard = discard(1:options.NSinos*options.Ndist*options.Nang);
        end
        if nargout > 9
            varargout{3} = discard;
        end
        lor_a = (lor(discard));
        if (options.projector_type == 2 || options.projector_type == 3) && options.implementation == 1
            if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                lor_orth = lor_orth(1:options.NSinos*options.Ndist*options.Nang);
            end
            lor_orth = (lor_orth(discard));
        end
        if options.normalization_correction && options.corrections_during_reconstruction
            if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                options.normalization = options.normalization(1:options.NSinos*options.Ndist*options.Nang);
            end
            options.normalization = options.normalization(discard);
        end
        if (options.randoms_correction || (options.scatter_correction  && options.subtract_scatter)) && options.corrections_during_reconstruction ...
                && ~options.reconstruct_trues && ~options.reconstruct_scatter
            if options.partitions > 1
                for ff = 1 : options.partitions
                    temp = options.SinDelayed{ff};
                    if options.NSinos ~= options.TotSinos
                        temp = temp(:,:,1:options.NSinos);
                    end
                    temp = temp(discard);
                    options.SinDelayed{ff} = temp;
                end
                clear temp
            else
                if iscell(options.SinDelayed)
                    options.SinDelayed{1} = options.SinDelayed{1}(discard);
                else
                    options.SinDelayed = options.SinDelayed(discard);
                end
            end
        end
        if options.scatter_correction && options.corrections_during_reconstruction && ~options.subtract_scatter ...
                && ~options.reconstruct_trues && ~options.reconstruct_scatter
            if options.partitions > 1 && iscell(options.ScatterC) && length(options.ScatterC) > 1
                for ff = 1 : options.partitions
                    temp = options.ScatterC{ff};
                    if options.NSinos ~= options.TotSinos
                        temp = temp(:,:,1:options.NSinos);
                    end
                    temp = temp(discard);
                    options.ScatterC{ff} = temp;
                end
                clear temp
            else
                if iscell(options.ScatterC)
                    options.ScatterC = options.ScatterC{1}(discard);
                else
                    options.ScatterC = options.ScatterC(discard);
                end
            end
        end
        if nargout >= 8 && ~isempty(varargin) && length(varargin) > 1 && ~isempty(varargin{2})
            if options.partitions > 1
                for ff = 1 : options.partitions
                    temp = varargin{2}{ff};
                    if options.NSinos ~= options.TotSinos
                        temp = temp(:,:,1:options.NSinos,:);
                    end
                    if TOF
                        temp = reshape(temp, numel(temp) / options.TOF_bins, options.TOF_bins);
                        temp = temp(discard,:);
                    else
                        temp = temp(discard);
                    end
                    varargout{1}{ff} = temp(:);
                end
                clear temp
            else
                if options.NSinos ~= options.TotSinos
                    varargin{2} = varargin{2}(:,:,1:options.NSinos,:);
                end
                if TOF
                    varargin{2} = reshape(varargin{2}, numel(varargin{2}) / options.TOF_bins, options.TOF_bins);
                    varargin{2} = varargin{2}(discard,:);
                else
                    varargin{2} = varargin{2}(discard);
                end
                varargout{1} = varargin{2}(:);
            end
        end
        clear lor
    end
    if options.CT
        if (subsets > 1 || fpbp) || options.precompute_lor
            xy_index = repmat(uint32(1:size_x * options.xSize)', options.nProjections, 1);
            z_index = uint16(1:options.nProjections)';
            if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                z_index = repeat_elem(z_index, options.xSize * size_x);
            else
                z_index = repelem(z_index, options.xSize * size_x);
            end
        else
            xy_index = uint32(1);
            z_index = uint16(1);
        end
    else
        xy_index = uint32(1:size_x)';
        if options.span > 1
            xy_index2 = repmat(uint32(1:size_x)', options.NSinos - options.Nz, 1);
            xy_index = [repmat(xy_index, options.Nz, 1); xy_index2];
        else
            xy_index2 = repmat(uint32(1:size_x)', options.NSinos - options.rings, 1);
            xy_index = [repmat(xy_index, options.rings, 1); xy_index2];
        end
        z_index = uint16(1:size(z_det,1))';
        if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
            z_index = repeat_elem(z_index, size_x);
        else
            z_index = repelem(z_index, size_x);
        end
    end
    if (subsets > 1 || fpbp)
        z_index = z_index(index);
    elseif ~options.CT || (options.precompute_lor && options.CT)
        z_index = (z_index(discard));
    end
    z_index = z_index - 1;
    
    if (subsets > 1 || fpbp)
        xy_index = xy_index(index);
    elseif ~options.CT || (options.precompute_lor && options.CT)
        xy_index = (xy_index(discard));
    end
    xy_index = xy_index - 1;
    
    summa = zeros(subsets, 1, 'uint64');
    
    if (subsets > 1 && length(pituus) > 1 || fpbp)
        for kk = 1 : subsets
            if (options.projector_type == 2 || options.projector_type == 3) && options.implementation == 1
                summa(kk) = uint64(sum(uint64(lor_orth(pituus(kk)+1:pituus(kk+1)))));
            else
                summa(kk) = uint64(sum(uint64(lor_a(pituus(kk)+1:pituus(kk+1)))));
            end
        end
    else
        %         if TOF
        %             pituus = int64([0;pituus]);
        %             koko = int64(size_x * options.NSinos);
        %             for kk = 1 : subsets
        %                 alku = pituus(kk) + 1;
        %                 loppu = pituus(kk + 1);
        %                 bin = floor(alku / koko);
        %                 if alku - loppu > koko || bin < floor(loppu / koko)
        %                     alku = alku - koko * bin;
        %                     loppu = loppu - koko * floor(loppu / koko);
        %                     if (options.projector_type == 2 || options.projector_type == 3) && options.implementation == 1
        %                         summa(kk) = uint64(sum(uint64(lor_orth(alku:end))) + sum(uint64(lor_orth(1:loppu))));
        %                     else
        %                         summa(kk) = uint64(sum(uint64(lor_a(alku:end))) + sum(uint64(lor_a(1:loppu))));
        %                     end
        %                 else
        %                     alku = alku - koko * bin;
        %                     loppu = loppu - koko * floor(loppu / koko);
        %                     if (options.projector_type == 2 || options.projector_type == 3) && options.implementation == 1
        %                         summa(kk) = uint64(sum(uint64(lor_orth(alku:loppu))));
        %                     else
        %                         summa(kk) = uint64(sum(uint64(lor_a(alku:loppu))));
        %                     end
        %                 end
        %             end
        %         else
        if (options.projector_type == 2 || options.projector_type == 3) && options.implementation == 1
            summa = uint64(sum(int64(lor_orth)));
        else
            summa = uint64(sum(int64(lor_a)));
        end
        pituus = int64([0;length(lor_a)]);
        %         end
    end
    if (options.projector_type == 2 || options.projector_type == 3) && options.implementation == 1
        varargout{2} = lor_orth;
    end
    
    LL = [];
    clear discard I yt xt xy_index2 index apu
elseif options.use_raw_data && options.precompute_lor
    
    lor_file = [folder options.machine_name '_detector_locations_' num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '_raw.mat'];
    if exist(lor_file, 'file') == 2
        variableInfo = who('-file', lor_file);
        if options.implementation == 1 || options.implementation == 4
            if ismember('lor', variableInfo)
                load(lor_file,'lor')
            else
                lor_pixel_count_prepass(options);
                load(lor_file,'lor')
            end
            if options.projector_type == 2 && options.implementation == 1
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
                    if length(lor_orth) > length(lor)
                        lor_orth = lor_orth(1:length(lor));
                    end
                elseif options.tube_width_z > 0
                    if crystal_size_z ~= options.tube_width_z
                        lor_pixel_count_prepass(options);
                        load(lor_file,'lor_orth')
                    end
                    if length(lor_orth) > length(lor)
                        lor_orth = lor_orth(length(lor)+1:end);
                    elseif length(lor_orth) == length(lor)
                        lor_pixel_count_prepass(options);
                        load(lor_file,'lor_orth')
                        lor_orth = lor_orth(length(lor)+1:end);
                    end
                end
            end
            if options.projector_type == 3 && options.implementation == 1
                if any(ismember('lor_vol', variableInfo))
                    load(lor_file,'lor_vol')
                else
                    lor_pixel_count_prepass(options);
                    load(lor_file,'lor_vol')
                end
                lor_orth = lor_vol;
                clear lor_vol
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
        if options.implementation == 1 || options.implementation == 4
            load(lor_file,'lor')
            if options.projector_type == 2 && options.implementation == 1
                load(lor_file,'lor_orth')
                if options.tube_width_z == 0
                    if length(lor_orth) > length(lor)
                        lor_orth = lor_orth(1:length(lor));
                    end
                elseif options.tube_width_z > 0
                    if length(lor_orth) > length(lor)
                        lor_orth = lor_orth(length(lor)+1:end);
                    end
                end
            end
            if options.projector_type == 3 && options.implementation == 1
                load(lor_file,'lor_vol')
                lor_orth = lor_vol;
                clear lor_vol
            end
        else
            load(lor_file,'lor_opencl')
            lor = lor_opencl;
            clear lor_opencl
        end
    end
    if options.ring_difference_raw < options.rings
        testi3 = zeros(options.detectors,options.detectors,'uint16');
        testi3(tril(true(size(testi3)), 0)) = lor;
        for kk = options.rings : - 1 : options.ring_difference_raw + 1
            for ll = 1 : kk - options.ring_difference_raw
                testi3(1 + (kk - 1) * options.det_per_ring : kk * options.det_per_ring, 1 + (ll - 1) * options.det_per_ring : ll * options.det_per_ring) = ...
                    zeros(options.det_per_ring, options.det_per_ring, 'uint16');
            end
        end
        lor = testi3(tril(true(size(testi3)), 0));
    end
    discard = lor > 0;
    if nargout > 9
        varargout{3} = discard;
    end
    if (subsets > 1 || fpbp)
        if ~exist('LL','var')
            LL = form_detector_pairs_raw(rings, options.det_per_ring);
        end
        LL = LL(discard,:);
        lor = lor(discard);
        if options.normalization_correction && options.corrections_during_reconstruction
            options.normalization = options.normalization(discard);
        end
        if (options.projector_type == 2 || options.projector_type == 3) && options.implementation == 1
            lor_orth = (lor_orth(discard));
        end
        if (options.randoms_correction || (options.scatter_correction  && options.subtract_scatter)) && options.corrections_during_reconstruction ...
                && ~options.reconstruct_trues && ~options.reconstruct_scatter
            if options.partitions > 1
                for ff = 1 : options.partitions
                    temp = single(full(options.SinDelayed{ff}));
                    temp = temp(discard);
                    temp = temp(index);
                    options.SinDelayed{ff} = temp;
                end
                clear temp
            else
                if iscell(options.SinDelayed)
                    options.SinDelayed{1} = options.SinDelayed{1}(discard);
                    options.SinDelayed{1} = options.SinDelayed{1}(index);
                else
                    options.SinDelayed = options.SinDelayed(discard);
                    options.SinDelayed = options.SinDelayed(index);
                end
            end
        end
        if options.scatter_correction  && ~options.subtract_scatter && options.corrections_during_reconstruction ...
                && ~options.reconstruct_trues && ~options.reconstruct_scatter
            if options.partitions > 1 && iscell(options.ScatterC) && length(options.ScatterC) > 1
                for ff = 1 : options.partitions
                    temp = single(full(options.ScatterC{ff}));
                    temp = temp(discard);
                    temp = temp(index);
                    options.ScatterC{ff} = temp;
                end
                clear temp
            else
                if iscell(options.ScatterC)
                    options.ScatterC = options.ScatterC{1}(discard);
                    options.ScatterC = options.ScatterC{1}(index);
                else
                    options.ScatterC = options.ScatterC(discard);
                    options.ScatterC = options.ScatterC(index);
                end
            end
        end
        LL = LL(index,:);
        lor_a = (lor(index));
        if nargout >= 8 && ~isempty(varargin) && length(varargin) > 1 && ~isempty(varargin{2})
            if options.partitions > 1
                for ff = 1 : options.partitions
                    temp = single(full(varargin{2}{ff}));
                    temp = temp(discard);
                    if subsets > 1 || fpbp
                        if TOF
                            temp = reshape(temp, numel(temp) / options.TOF_bins, options.TOF_bins);
                            temp = temp(index,:);
                        else
                            temp = temp(index);
                        end
                    end
                    varargout{1}{ff} = temp(:);
                end
                clear temp
            else
                if iscell(varargin{2})
                    temp = single(full(varargin{2}{1}));
                    if TOF
                        temp = reshape(temp, numel(temp) / options.TOF_bins, options.TOF_bins);
                        temp = temp(discard,:);
                        temp = temp(index,:);
                    else
                        temp = temp(discard);
                        temp = temp(index);
                    end
                    varargout{1}{1} = temp(:);
                    clear temp
                    %                     varargout{1}{1} = varargin{2}{1}(discard);
                    %                     varargout{1}{1} = varargout{1}{1}(index);
                else
                    temp = single(full(varargin{2}));
                    if TOF
                        temp = reshape(temp, numel(temp) / options.TOF_bins, options.TOF_bins);
                        temp = temp(discard,:);
                        temp = temp(index,:);
                    else
                        temp = temp(discard);
                        temp = temp(index);
                    end
                    varargout{1} = temp(:);
                    clear temp
                    %                     varargout{1} = varargin{2}(discard);
                    %                     varargout{1} = varargout{1}(index);
                end
            end
        end
        if options.normalization_correction && options.corrections_during_reconstruction
            options.normalization = options.normalization(index);
        end
        if (options.projector_type == 2 || options.projector_type == 3) && options.implementation == 1
            lor_orth = (lor_orth(index));
        end
        clear lor
    else
        if ~exist('LL','var')
            LL = form_detector_pairs_raw(rings, options.det_per_ring);
            LL = LL(discard,:);
        end
        if nargout >= 8 && ~isempty(varargin) && length(varargin) > 1 && ~isempty(varargin{2})
            if options.partitions > 1
                for ff = 1 : options.partitions
                    temp = single(full(varargin{2}{ff}));
                    if TOF
                        temp = reshape(temp, numel(temp) / options.TOF_bins, options.TOF_bins);
                        temp = temp(discard, :);
                    else
                        temp = temp(discard);
                    end
                    varargout{1}{ff} = temp(:);
                end
                clear temp
            else
                if iscell(varargin{2})
                    apu = single(full(varargin{2}{1}));
                    if TOF
                        apu = reshape(apu, numel(apu) / options.TOF_bins, options.TOF_bins);
                        apu = apu(discard, :);
                        varargout{1}{1} = apu(:);
                    else
                        varargout{1}{1} = (apu(discard));
                    end
                    clear apu
                else
                    apu = single(full(varargin{2}));
                    if TOF
                        apu = reshape(apu, numel(apu) / options.TOF_bins, options.TOF_bins);
                        apu = apu(discard, :);
                        varargout{1} = apu(:);
                    else
                        varargout{1} = (apu(discard));
                    end
                    clear apu
                end
            end
        end
        lor_a = lor(discard);
        if options.normalization_correction && options.corrections_during_reconstruction
            options.normalization = options.normalization(discard);
        end
        if (options.randoms_correction || (options.scatter_correction  && options.subtract_scatter)) && options.corrections_during_reconstruction ...
                && ~options.reconstruct_trues && ~options.reconstruct_scatter
            if options.partitions > 1
                for ff = 1 : options.partitions
                    temp = single(full(options.SinDelayed{ff}));
                    temp = temp(discard);
                    options.SinDelayed{ff} = temp;
                end
                clear temp
            else
                if iscell(options.SinDelayed)
                    %                     options.SinDelayed{1} = options.SinDelayed{1}(discard);
                    temp = single(full(options.SinDelayed{1}));
                    temp = temp(discard);
                    options.SinDelayed{1} = temp;
                    clear temp
                else
                    %                     options.SinDelayed = options.SinDelayed(discard);
                    temp = single(full(options.SinDelayed));
                    temp = temp(discard);
                    options.SinDelayed = temp;
                    clear temp
                end
            end
        end
        if options.scatter_correction  && ~options.subtract_scatter && options.corrections_during_reconstruction ...
                && ~options.reconstruct_trues && ~options.reconstruct_scatter
            if options.partitions > 1 && iscell(options.ScatterC) && length(options.ScatterC) > 1
                for ff = 1 : options.partitions
                    temp = single(full(options.ScatterC{ff}));
                    temp = temp(discard);
                    options.ScatterC{ff} = temp;
                end
                clear temp
            else
                if iscell(options.ScatterC)
                    %                     options.ScatterC = options.ScatterC{1}(discard);
                    temp = single(full(options.ScatterC{1}));
                    temp = temp(discard);
                    options.ScatterC{1} = temp;
                    clear temp
                else
                    %                     options.ScatterC = options.ScatterC(discard);
                    temp = single(full(options.ScatterC));
                    temp = temp(discard);
                    options.ScatterC = temp;
                    clear temp
                end
            end
        end
        pituus = int64([0;length(lor_a)]);
        if (options.projector_type == 2 || options.projector_type == 3) && options.implementation == 1
            lor_orth = (lor_orth(discard));
        end
        clear lor
    end
    summa = zeros(subsets, 1, 'uint64');
    
    %     if ~TOF
    for kk = 1 : subsets
        apu = LL(pituus(kk) + 1 : pituus(kk + 1),:) - 1;
        apu2 = idivide(apu, uint16(options.det_per_ring));
        idx = apu2(:,1) == apu2(:,2);
        apu2 = apu(idx,:);
        ind = mod(apu2, uint16(options.det_per_ring)) + 1;
        yt = y(ind);
        y_i = yt(:,1) > yt(:,2);
        apu2(y_i,:) = fliplr(apu2(y_i,:));
        apu(idx,:) = apu2;
        LL(pituus(kk) + 1 : pituus(kk + 1),:) = apu + 1;
        if (options.projector_type == 2 || options.projector_type == 3) && options.implementation == 1
            summa(kk) = uint64(sum(int64(lor_orth(pituus(kk)+1:pituus(kk+1)))));
        else
            summa(kk) = uint64(sum(int64(lor_a(pituus(kk)+1:pituus(kk+1)))));
        end
    end
    %     else
    %         apu = LL - 1;
    %         apu2 = idivide(apu, uint16(options.det_per_ring));
    %         idx = apu2(:,1) == apu2(:,2);
    %         apu2 = apu(idx,:);
    %         ind = mod(apu2, uint16(options.det_per_ring)) + 1;
    %         yt = y(ind);
    %         y_i = yt(:,1) > yt(:,2);
    %         apu2(y_i,:) = fliplr(apu2(y_i,:));
    %         apu(idx,:) = apu2;
    %         LL = apu + 1;
    %
    %         koko = length(lor_a);
    %         for kk = 1 : subsets
    %             alku = pituus(kk) + 1;
    %             loppu = pituus(kk + 1);
    %             bin = floor(alku / koko);
    %             if alku - loppu > koko || bin < floor(loppu / koko)
    %                 alku = alku - koko * bin;
    %                 loppu = loppu - koko * floor(loppu / koko);
    %                 if (options.projector_type == 2 || options.projector_type == 3) && options.implementation == 1
    %                     summa(kk) = uint64(sum(uint64(lor_orth(alku:end))) + sum(uint64(lor_orth(1:loppu))));
    %                 else
    %                     summa(kk) = uint64(sum(uint64(lor_a(alku:end))) + sum(uint64(lor_a(1:loppu))));
    %                 end
    %             else
    %                 alku = alku - koko * bin;
    %                 loppu = loppu - koko * floor(loppu / koko);
    %                 if (options.projector_type == 2 || options.projector_type == 3) && options.implementation == 1
    %                     summa(kk) = uint64(sum(uint64(lor_orth(alku:loppu))));
    %                 else
    %                     summa(kk) = uint64(sum(uint64(lor_a(alku:loppu))));
    %                 end
    %             end
    %         end
    %     end
    
    clear apu apu2 idx ind yt y_i index discard
    
    LL = LL';
    LL = LL(:);
    xy_index =  [];
    z_index = [];
    if (options.projector_type == 2 || options.projector_type == 3) && options.implementation == 1
        varargout{2} = lor_orth;
    end
elseif options.use_raw_data == false && ~options.precompute_lor
    
    if (subsets > 1 || fpbp)
        if options.NSinos ~= options.TotSinos
            index = index(1:options.Ndist*options.Nang*options.NSinos);
        end
        if options.normalization_correction && options.corrections_during_reconstruction
            if options.NSinos ~= options.TotSinos
                options.normalization = options.normalization(index(1:options.Ndist*options.Nang*options.NSinos));
            else
                options.normalization = options.normalization(index);
            end
        end
        if (options.randoms_correction || (options.scatter_correction  && options.subtract_scatter)) && options.corrections_during_reconstruction ...
                && ~options.reconstruct_trues && ~options.reconstruct_scatter
            if options.partitions > 1
                for ff = 1 : options.partitions
                    temp = options.SinDelayed{ff};
                    if options.NSinos ~= options.TotSinos
                        temp = temp(:,:,1:options.NSinos);
                    end
                    temp = temp(index);
                    options.SinDelayed{ff} = temp;
                end
                clear temp
            else
                if iscell(options.SinDelayed)
                    options.SinDelayed{1} = options.SinDelayed{1}(index);
                else
                    options.SinDelayed = options.SinDelayed(index);
                end
            end
        end
        if options.scatter_correction  && ~options.subtract_scatter && options.corrections_during_reconstruction ...
                && ~options.reconstruct_trues && ~options.reconstruct_scatter
            if options.partitions > 1 && iscell(options.ScatterC) && length(options.ScatterC) > 1
                for ff = 1 : options.partitions
                    temp = options.ScatterC{ff};
                    if options.NSinos ~= options.TotSinos
                        temp = temp(:,:,1:options.NSinos);
                    end
                    temp = temp(index);
                    options.ScatterC{ff} = temp;
                end
                clear temp
            else
                if iscell(options.ScatterC)
                    options.ScatterC = options.ScatterC{1}(index);
                else
                    options.ScatterC = options.ScatterC(index);
                end
            end
        end
        if nargout >= 8 && ~isempty(varargin) && length(varargin) > 1 && ~isempty(varargin{2})
            if options.partitions > 1
                for ff = 1 : options.partitions
                    temp = varargin{2}{ff};
                    if options.NSinos ~= options.TotSinos
                        temp = temp(:,:,1:options.NSinos,:);
                        if TOF
                            temp = reshape(temp, numel(temp) / options.TOF_bins, options.TOF_bins);
                        end
                    else
                        if TOF
                            temp = reshape(temp, numel(temp) / options.TOF_bins, options.TOF_bins);
                        end
                    end
                    temp = temp(index);
                    varargout{1}{ff} = temp(:);
                end
                clear temp
            else
                if options.NSinos ~= options.TotSinos
                    varargin{2} = varargin{2}(:,:,1:options.NSinos,:);
                    if TOF
                        varargout{1} = reshape(varargin{2}, numel(varargin{2}) / options.TOF_bins, options.TOF_bins);
                        varargout{1} = varargout{1}(index,:);
                    else
                        varargout{1} = varargin{2}(index);
                    end
                else
                    if TOF
                        varargout{1} = reshape(varargin{2}, numel(varargin{2}) / options.TOF_bins, options.TOF_bins);
                        varargout{1} = varargout{1}(index,:);
                    else
                        varargout{1} = varargin{2}(index);
                    end
                end
                varargout{1} = varargout{1}(:);
            end
        end
    else
        if options.normalization_correction && options.corrections_during_reconstruction
            if options.NSinos ~= options.TotSinos
                options.normalization = options.normalization(1:options.Ndist*options.Nang*options.NSinos);
            end
            options.normalization = options.normalization(:);
        end
        if (options.randoms_correction || (options.scatter_correction  && options.subtract_scatter)) && options.corrections_during_reconstruction ...
                && ~options.reconstruct_trues && ~options.reconstruct_scatter
            if options.partitions > 1
                for ff = 1 : options.partitions
                    temp = options.SinDelayed{ff};
                    if options.NSinos ~= options.TotSinos
                        temp = temp(:,:,1:options.NSinos);
                    end
                    temp = temp(:);
                    options.SinDelayed{ff} = temp;
                end
                clear temp
            else
                if iscell(options.SinDelayed)
                    if options.NSinos ~= options.TotSinos
                        options.SinDelayed{1} = options.SinDelayed{1}(1:options.Ndist*options.Nang*options.NSinos);
                    end
                    options.SinDelayed{1} = options.SinDelayed{1}(:);
                else
                    if options.NSinos ~= options.TotSinos
                        options.SinDelayed = options.SinDelayed(1:options.Ndist*options.Nang*options.NSinos);
                    end
                    options.SinDelayed = options.SinDelayed(:);
                end
            end
        end
        if options.scatter_correction  && ~options.subtract_scatter && options.corrections_during_reconstruction ...
                && ~options.reconstruct_trues && ~options.reconstruct_scatter
            if options.partitions > 1 && iscell(options.ScatterC) && length(options.ScatterC) > 1
                for ff = 1 : options.partitions
                    temp = options.ScatterC{ff};
                    if options.NSinos ~= options.TotSinos
                        temp = temp(:,:,1:options.NSinos);
                    end
                    temp = temp(:);
                    options.ScatterC{ff} = temp;
                end
                clear temp
            else
                if iscell(options.ScatterC)
                    if options.NSinos ~= options.TotSinos
                        options.ScatterC = options.ScatterC{1}(1:options.Ndist*options.Nang*options.NSinos);
                    end
                    options.ScatterC = options.ScatterC{1}(:);
                else
                    if options.NSinos ~= options.TotSinos
                        options.ScatterC = options.ScatterC(1:options.Ndist*options.Nang*options.NSinos);
                    end
                    options.ScatterC = options.ScatterC(:);
                end
            end
        end
        if nargout >= 8 && ~isempty(varargin) && length(varargin) > 1 && ~isempty(varargin{2})
            if options.partitions > 1
                for ff = 1 : options.partitions
                    temp = varargin{2}{ff};
                    if options.NSinos ~= options.TotSinos
                        temp = temp(:,:,1:options.NSinos,:);
                    end
                    varargout{1}{ff} = temp(:);
                end
                clear temp
            else
                if options.NSinos ~= options.TotSinos
                    varargout{1} = varargin{2}(:,:,1:options.NSinos,:);
                end
                varargout{1} = varargout{1}(:);
            end
        end
    end
    if options.CT
        if (subsets > 1 || fpbp)
            xy_index = repmat(uint32(1:size_x * options.xSize)', options.nProjections, 1);
            z_index = uint16(1:options.nProjections)';
            if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                z_index = repeat_elem(z_index, options.xSize * size_x);
            else
                z_index = repelem(z_index, options.xSize * size_x);
            end
        else
            xy_index = uint32(1);
            z_index = uint16(1);
        end
    else
        xy_index = uint32(1:size_x)';
        if options.span > 1
            xy_index2 = repmat(uint32(1:size_x)', options.NSinos - options.Nz, 1);
            xy_index = [repmat(xy_index, options.Nz, 1); xy_index2];
        else
            xy_index2 = repmat(uint32(1:size_x)', options.NSinos - options.rings, 1);
            xy_index = [repmat(xy_index, options.rings, 1); xy_index2];
        end
        z_index = uint16(1:size(z_det,1))';
        if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
            z_index = repeat_elem(z_index, size_x);
        else
            z_index = repelem(z_index, size_x);
        end
    end
    if (subsets > 1 || fpbp)
        z_index = z_index(index);
    end
    z_index = z_index - 1;
    
    if (subsets > 1 || fpbp)
        xy_index = xy_index(index);
    end
    xy_index = xy_index - 1;
    
    summa = zeros(subsets, 1, 'uint64');
    
    if (subsets > 1 || fpbp)
    else
        pituus = int64([0;options.Nang*options.Ndist*options.NSinos]);
    end
    
    lor_a = [];
    LL = [];
    clear discard I yt xt xy_index2 index apu
elseif options.use_raw_data && ~options.precompute_lor
    det_per_ring = options.det_per_ring * options.sampling_raw;
    if ~exist('LL','var')
        LL = form_detector_pairs_raw(rings, det_per_ring);
    end
    if subsets > 1 || fpbp
        LL = LL(index,:);
        if options.normalization_correction && options.corrections_during_reconstruction
            options.normalization = options.normalization(index);
        end
        if (options.randoms_correction || (options.scatter_correction  && options.subtract_scatter)) && options.corrections_during_reconstruction ...
                && ~options.reconstruct_trues && ~options.reconstruct_scatter
            if options.partitions > 1
                for ff = 1 : options.partitions
                    temp = single(full(options.SinDelayed{ff}));
                    temp = temp(index);
                    options.SinDelayed{ff} = temp;
                end
                clear temp
            else
                if iscell(options.SinDelayed)
                    %                     options.SinDelayed{1} = options.SinDelayed{1}(index);
                    temp = single(full(options.SinDelayed{1}));
                    temp = temp(index);
                    options.SinDelayed{1} = temp;
                    clear temp
                else
                    %                     options.SinDelayed = options.SinDelayed(index);
                    temp = single(full(options.SinDelayed));
                    temp = temp(index);
                    options.SinDelayed = temp;
                    clear temp
                end
            end
        end
        if options.scatter_correction  && ~options.subtract_scatter && options.corrections_during_reconstruction ...
                && ~options.reconstruct_trues && ~options.reconstruct_scatter
            if options.partitions > 1 && iscell(options.ScatterC) && length(options.ScatterC) > 1
                for ff = 1 : options.partitions
                    temp = single(full(options.ScatterC{ff}));
                    temp = temp(index);
                    options.ScatterC{ff} = temp;
                end
                clear temp
            else
                if iscell(options.SinDelayed)
                    %                     options.ScatterC = options.ScatterC{1}(index);
                    temp = single(full(options.ScatterC{1}));
                    temp = temp(index);
                    options.ScatterC{1} = temp;
                    clear temp
                else
                    %                     options.ScatterC = options.ScatterC(index);
                    temp = single(full(options.ScatterC));
                    temp = temp(index);
                    options.ScatterC = temp;
                    clear temp
                end
            end
        end
        if nargout >= 8 && ~isempty(varargin) && length(varargin) > 1 && ~isempty(varargin{2})
            if options.partitions > 1
                for ff = 1 : options.partitions
                    temp = single(full(varargin{2}{ff}));
                    if subsets > 1 || fpbp
                        if TOF
                            temp = reshape(temp, numel(temp) / options.TOF_bins, options.TOF_bins);
                            temp = temp(index,:);
                        else
                            temp = temp(index);
                        end
                    end
                    varargout{1}{ff} = temp(:);
                end
                clear temp
            else
                if iscell(varargin{2})
                    apu = single(full(varargin{2}{1}));
                    if TOF
                        apu = reshape(apu, numel(apu) / options.TOF_bins, options.TOF_bins);
                        varargout{1}{1} = (apu(index,:));
                    else
                        varargout{1}{1} = (apu(index));
                    end
                    clear apu
                    varargout{1}{1} = varargout{1}{1}(:);
                else
                    apu = single(full(varargin{2}));
                    if TOF
                        apu = reshape(apu, numel(apu) / options.TOF_bins, options.TOF_bins);
                        varargout{1} = (apu(index,:));
                    else
                        varargout{1} = (apu(index));
                    end
                    clear apu
                    varargout{1} = varargout{1}(:);
                end
            end
        end
        clear lor
    end
    summa = zeros(subsets, 1, 'uint64');
    
    if subsets > 1
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
    else
        apu = LL - 1;
        apu2 = idivide(apu, uint16(det_per_ring));
        idx = apu2(:,1) == apu2(:,2);
        apu2 = apu(idx,:);
        ind = mod(apu2, uint16(det_per_ring)) + 1;
        yt = y(ind);
        y_i = yt(:,1) > yt(:,2);
        apu2(y_i,:) = fliplr(apu2(y_i,:));
        apu(idx,:) = apu2;
        LL = apu + 1;
    end
    
    clear apu apu2 idx ind yt y_i index discard
    
    LL = LL';
    LL = LL(:);
    lor_a = [];
    xy_index = [];
    z_index = [];
end
if options.sampling > 1 && ~options.use_raw_data && ~options.precompute_lor
    options.Ndist = options.Ndist / options.sampling;
end
