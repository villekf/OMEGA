function [index, pituus, subsets, varargout] = index_maker(Nx, Ny, Nz, subsets, use_raw_data, machine_name, options, varargin)
% INDEX_MAKER Form the subset indices for any of the subset types
%
% Example:
%   [index, pituus, subsets] = index_maker(Nx, Ny, Nz, subsets,
%   use_raw_data, machine_name, options, Nang, Ndist, TotSinos, NSinos)
% INPUTS:
%   Nx, Ny, Nz = Image (matrix) size in x-, y-, and z-directions
%   subsets = Number of subsets used
%   use_raw_data = True (or 1) if raw lis-mode data is used, false (or 0)
%   if sinogram data is used
%   machine_name = Name of the machine, in order to load precomputed data
%   options = Holds various misc parameters (subset type, no. rings,
%   detectors per ring, etc.)
%   Nang = The number of angles in the sinogram (necessary if sinogram data
%   is used)
%   Ndist = Number of angular positions in sinogram (necessary if sinogram
%   data is used)
%   TotSinos = Total number of sinograms in the sinogram (necessary if
%   sinogram data is used)
%   NSinos = Number of sinograms used (e.g. take only the 2D sinograms)
%   (necessary if sinogram data is used)
%
% OUTPUTS:
%   index = All indices in the order specified by the subset type
%   pituus = Number of indices per subset
%   subsets = If using subset_type = 5, this updates the number of subsets
%   used, otherwise returns the previously used value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2022 Ville-Veikko Wettenhovi
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


if nargin - 7 < 4 && ~use_raw_data
    error('Sinogram data selected, but not enough sinogram size information input')
end
if nargin - 7 >= 4 && ~isempty(varargin) && ~isempty(varargin{1}) && ~isempty(varargin{2}) && ~isempty(varargin{3}) && ~isempty(varargin{4})
    Nang = varargin{1};
    Ndist = varargin{2};
    TotSinos = varargin{3};
    NSinos = varargin{4};
    
    if options.sampling > 1 && ~options.use_raw_data && ~options.precompute_lor
        Ndist = Ndist * options.sampling;
    end
end
if nargin >= 12 && ~isempty(varargin) && ~isempty(varargin{5})
    storeMatrix = varargin{5};
else
    storeMatrix = false;
end

if ~isfield(options,'sampling_raw')
    options.sampling_raw = 1;
end
if ~isfield(options,'CT')
    options.CT = false;
end
folder = fileparts(which('index_maker.m'));
folder = [folder(1:end-6), 'mat-files/'];
folder = strrep(folder, '\','/');
% lor_a = 0;
if options.use_raw_data
    pituus = int64(options.detectors^2/2 + options.detectors/2);
else
    pituus = int64(Ndist * Nang * NSinos);
end
index = 0;
lor = [];
lor_orth = [];
if options.precompute_lor
    if options.CT
        [lor, ~, ~, lor_orth] = lor_pixel_count_prepass(options, false);
    elseif ~options.CT
        if ~options.use_raw_data
            lor_file = [folder machine_name '_lor_pixel_count_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_sino_' num2str(Ndist) 'x' ...
                num2str(Nang) 'x' num2str(TotSinos) '.mat'];
        else
            lor_file = [folder machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'];
        end
        if exist(lor_file, 'file') == 2
            if options.implementation == 1 || options.implementation == 4
                variableInfo = who('-file', lor_file);
                if ismember('lor', variableInfo)
                    load(lor_file,'lor')
                else
                    lor_pixel_count_prepass(options);
                    load(lor_file,'lor')
                end
            else
                variableInfo = who('-file', lor_file);
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
            else
                load(lor_file,'lor_opencl')
                lor = lor_opencl;
                clear lor_opencl
            end
        end
    end
end
if options.CT
    tyyppi = 'uint64';
else
    tyyppi = 'uint32';
end
% Subset data
if subsets > 1 && options.subset_type < 8
    if options.use_raw_data
        if options.precompute_lor
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
            totalLength = numel(lor);
        else
            det_per_ring = options.det_per_ring * options.sampling_raw;
            % Same as above, but for precompute_lor = false case
            totalLength = sum(1 : det_per_ring * options.rings);
            if options.ring_difference_raw < options.rings
                if options.sampling_raw > 1
                    error('Increased sampling cannot be used with smaller ring difference than the number of rings!')
                end
                testi = zeros(options.detectors,options.detectors,tyyppi);
                testi(tril(true(size(testi)), 0)) = (1:totalLength);
                for kk = options.rings : - 1 : options.ring_difference_raw + 1
                    for ll = 1 : kk - options.ring_difference_raw
                        testi(1 + (kk - 1) * options.det_per_ring : kk * options.det_per_ring, 1 + (ll - 1) * options.det_per_ring : ll * options.det_per_ring) = ...
                            zeros(options.det_per_ring, options.det_per_ring, tyyppi);
                    end
                end
                ind = testi(tril(true(size(testi)), 0));
            end
        end
    else
        totalLength = Ndist*Nang*NSinos;
    end
    index = cell(subsets,1);
    pituus = zeros(subsets, 1, 'int64');
    % Take every nth column from the sinogram
    if options.subset_type == 4 && ~options.use_raw_data
        for i=1:subsets
            osa = length(i:subsets:Nang);
            index1 = repmat(cast((i-1)*Ndist+1:i*Ndist, tyyppi)', osa*NSinos,1);
            if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                index1 = index1 + cast(repeat_elem((0:(osa)*NSinos-1)'*Ndist*subsets,Ndist), tyyppi);
            else
                index1 = index1 + cast(repelem((0:(osa)*NSinos-1)'*Ndist*subsets,Ndist), tyyppi);
            end
            index{i} = index1;
            pituus(i) = int64(length(index{i}));
        end
        % Take every nth row from the sinogram
    elseif options.subset_type == 5 && ~options.use_raw_data
        for i=1:subsets
            osa = length(i:subsets:Ndist);
            index1 = repmat(uint32((i-1)+1:Ndist:Nang*Ndist)', osa*NSinos,1);
            if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                index1 = index1 + repmat(cast(repeat_elem((0:subsets:Ndist-1)',Nang), tyyppi), NSinos,1) + repeat_elem(cast(0:Ndist*Nang:totalLength-1, tyyppi)',Nang*osa);
            else
                index1 = index1 + repmat(cast(repelem((0:subsets:Ndist-1)',Nang), tyyppi), NSinos,1) + repelem(cast(0:Ndist*Nang:totalLength-1, tyyppi)',Nang*osa);
            end
            index{i} = index1;
            pituus(i) = int64(length(index{i}));
        end
        % Take every nth (column) measurement
    elseif options.subset_type == 1 && ~options.use_raw_data
        for i=1:subsets
            index1 = cast(i:subsets:totalLength, tyyppi)';
            [I,J,K] = ind2sub([Nang Ndist NSinos], index1);
            index1 = cast(sub2ind([Ndist Nang NSinos], J, I,K), tyyppi);
            index{i} = index1;
            pituus(i) = int64(length(index{i}));
        end
        % Take every nth (row) measurement
        % Every nth measurements
    elseif options.subset_type == 2 || (options.subset_type ~= 3 && options.subset_type ~= 6 && options.use_raw_data)
        for i=1:subsets
            index1 = cast(i:subsets:totalLength, tyyppi)';
            if options.use_raw_data && ~options.precompute_lor && options.ring_difference_raw < options.rings
                index1 = index1(ismember(index1,ind));
            end
            index{i} = index1;
            pituus(i) = int64(length(index{i}));
        end
        % Pick the measurements randomly
    elseif options.subset_type == 3
        indices = cast(totalLength, tyyppi);
        port = cast(floor(totalLength/subsets), tyyppi);
        if options.use_Shuffle && exist('Shuffle','file') == 3
            apu = Shuffle(indices(end), 'index')';
        elseif options.use_Shuffle && exist('Shuffle','file') == 0
            warning('options.Shuffle was set to true, but no Shuffle mex-file found. Using randperm')
            apu = cast(randperm(indices(end)), tyyppi)';
        else
            apu = cast(randperm(indices(end)), tyyppi)';
        end
        for i = 1 : subsets
            if i == subsets
                index1 = cast(apu(port*(i-1)+1:end), tyyppi);
            else
                index1 = cast(apu(port*(i-1)+1:(port*(i))), tyyppi);
            end
            if options.use_raw_data && ~options.precompute_lor && options.ring_difference_raw < options.rings
                index1 = index1(ismember(index1,ind));
            end
            index{i} = index1;
            pituus(i) = int64(length(index{i}));
        end
        % Pick the subsets based on the angles of the LORs
    elseif options.subset_type == 6
        if options.use_raw_data
            LL = form_detector_pairs_raw(options.rings, options.det_per_ring);
            if options.precompute_lor
                % Use only LORs that go through the FOV
                LL = LL(lor > 0,:);
            end
        else
            LL = [];
        end
        [index, pituus] = subset_angles(options, LL);
        subsets = length(pituus);
        if ~options.use_raw_data
            if options.NSinos < options.TotSinos
                if options.NSinos == options.Nz
                    subsets = 180 / options.n_angles;
                else
                    error(['Number of sinograms with subset_type = 6 has to be either ' num2str(options.Nz) ' or ' num2str(options.TotSinos)]);
                end
            end
        end
        % Use golden angle sampling
    elseif options.subset_type == 7 && ~options.use_raw_data
        [index, pituus] = goldenAngleSubsets(options);
    end
        
        %         lor_a = lor;
        %         clear lor
        
    if options.precompute_lor
        if iscell(index)
            index = cell2mat(index);
        end
        if ~storeMatrix
            discard = lor > 0;
            if ~options.use_raw_data && length(discard) ~= TotSinos*Ndist*Nang
                error('Error: Size mismatch between sinogram and LORs to be removed')
            end
            if use_raw_data == false && NSinos ~= TotSinos
                discard = discard(1:totalLength);
            end
            if options.CT
                ind_apu = uint64(find(discard));
            else
                ind_apu = uint32(find(discard));
            end
            % Remove indices that do not go through the FOV
            joku = ismember(index, ind_apu);
            pituus2 = cumsum(pituus);
            for kk = 1 : length(pituus2)
                if kk == 1
                    pituus(kk) = int64(sum(joku(1:pituus2(kk))));
                else
                    pituus(kk) = int64(sum(joku(1+pituus2(kk-1) : pituus2(kk))));
                end
            end
            index = index(joku);
        else
            pituus2 = cumsum(pituus);
            for kk = 1 : length(pituus2)
                if kk == 1
                    pituus(kk) = int64(sum((1:pituus2(kk))));
                else
                    pituus(kk) = int64(sum((1+pituus2(kk-1) : pituus2(kk))));
                end
            end
        end
        %         if options.CT
        %             lor = lor(index);
        %         end
    end
elseif (subsets > 1 && (options.subset_type == 8 || options.subset_type == 9)) || (subsets == 1  && ~options.precompute_lor)
    sProjections = floor(options.nProjections / subsets);
    pituus = zeros(subsets, 1, 'int64');
    index = cell(subsets, 1);
    if options.subset_type == 8
        for kk = 1 : subsets
            index{kk} = (kk : subsets : options.nProjections)';
            pituus(kk) = numel(index{kk});
        end
    elseif options.subset_type == 9
        if options.use_Shuffle && exist('Shuffle','file') == 3
            apu = Shuffle(options.nProjections, 'index')';
        elseif options.use_Shuffle && exist('Shuffle','file') == 0
            warning('options.Shuffle was set to true, but no Shuffle mex-file found. Using randperm')
            apu = randperm(options.nProjections)';
        else
            apu = randperm(options.nProjections)';
        end
        for kk = 1 : subsets
            if kk == subsets
                index{kk} = apu(sProjections*(kk-1)+1:end);
            else
                index{kk} = apu(sProjections*(kk-1)+1:(sProjections*(kk)));
            end
            pituus(kk) = numel(index{kk});
        end
    else
        index{1} = (1 : options.nProjections)';
        pituus(1) = numel(index{1});
    end
    index = cell2mat(index);
%     pituus = int64(cumsum(pituus));
elseif subsets == 1 && options.precompute_lor
%     if options.CT
%         index = uint64(1 : numel(lor))';
%     else
%         index = uint32(1 : pituus)';
%     end
    index = lor > 0;
    pituus = int64(sum(index));
end
if ~iscell(index) && size(index,1) == 1 && ~options.precompute_lor
%     index = true(pituus,1);
    index = 0;
%     if options.CT
%         index = uint64(1 : pituus)';
%     else
%         index = uint32(1 : pituus)';
%     end
end
if nargout >= 4
    varargout{1} = lor;
end
if nargout >=5
    varargout{2} = lor_orth;
end