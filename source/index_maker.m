function [index, pituus, subsets] = index_maker(Nx, Ny, Nz, subsets, use_raw_data, machine_name, options, varargin)
% INDEX_MAKER Form the subset indices for any of the subset types
%
% Example:
%   [index, pituus, subsets] = index_maker(Nx, Ny, Nz, subsets,
%   use_raw_data, machine_name, options, Nang, Ndist, TotSinos, NSinos)
% INPUTS:
%   Nx, Ny, Nz = Image (matrix) size in x-, y-, and z-directions
%   subsets = Number of subsets used
%   use_raw_data = True (or 1) if raw list-mode data is used, false (or 0)
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

if ~isfield(options,'sampling_raw')
    options.sampling_raw = 1;
end
folder = fileparts(which('index_maker.m'));
folder = strrep(folder, 'source','mat-files/');
folder = strrep(folder, '\','/');
% lor_a = 0;
if options.use_raw_data
    pituus = int64(options.detectors^2/2 + options.detectors/2);
else
    pituus = int64(Ndist * Nang * NSinos);
end
index = 0;
% Sinogram data
if use_raw_data == false && subsets > 1
    index = cell(subsets,1);
    pituus = zeros(subsets, 1, 'int64');
    % Take every nth column from the sinogram
    if options.subset_type == 4
        for i=1:subsets
            osa = length(i:subsets:Nang);
            index1 = repmat(uint32((i-1)*Ndist+1:i*Ndist)', osa*NSinos,1);
            if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                index1 = index1 + uint32(repeat_elem((0:(osa)*NSinos-1)'*Ndist*subsets,Ndist));
            else
                index1 = index1 + uint32(repelem((0:(osa)*NSinos-1)'*Ndist*subsets,Ndist));
            end
            index{i} = index1;
            pituus(i) = int64(length(index{i}));
        end
        % Take every nth row from the sinogram
    elseif options.subset_type == 5
        for i=1:subsets
            osa = length(i:subsets:Ndist);
            index1 = repmat(uint32((i-1)+1:Ndist:Nang*Ndist)', osa*NSinos,1);
            if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                index1 = index1 + repmat(uint32(repeat_elem((0:subsets:Ndist-1)',Nang)), NSinos,1) + repeat_elem(uint32(0:Ndist*Nang:Ndist*Nang*NSinos-1)',Nang*osa);
            else
                index1 = index1 + repmat(uint32(repelem((0:subsets:Ndist-1)',Nang)), NSinos,1) + repelem(uint32(0:Ndist*Nang:Ndist*Nang*NSinos-1)',Nang*osa);
            end
            index{i} = index1;
            pituus(i) = int64(length(index{i}));
        end
        % Take every nth (column) measurement
    elseif options.subset_type == 1
        for i=1:subsets
            index1 = uint32(i:subsets:Ndist*Nang*NSinos)';
            [I,J,K] = ind2sub([Nang Ndist NSinos], index1);
            index1 = uint32(sub2ind([Ndist Nang NSinos], J, I,K));
            index{i} = index1;
            pituus(i) = int64(length(index{i}));
        end
        % Take every nth (row) measurement
    elseif options.subset_type == 2
        for i=1:subsets
            index1 = uint32(i:subsets:Ndist*Nang*NSinos)';
            index{i} = index1;
            pituus(i) = int64(length(index{i}));
        end
        % Pick the measurements randomly
    elseif options.subset_type == 3
        indices = uint32(Ndist*Nang*NSinos);
        port = uint32(floor(Ndist*Nang*NSinos/subsets));
        if options.use_Shuffle
            apu = Shuffle(indices(end), 'index')';
        else
            apu = uint32(randperm(indices(end)))';
        end
        for i = 1 : subsets
            if i == subsets
                index1 = uint32(apu(port*(i-1)+1:end));
            else
                index1 = uint32(apu(port*(i-1)+1:(port*(i))));
            end
            index{i} = index1;
            pituus(i) = int64(length(index{i}));
        end
        % Pick the subsets based on the angles of the LORs
    elseif options.subset_type == 6
        [index, pituus] = subset_angles(options);
        subsets = length(pituus);
        if options.NSinos < options.TotSinos
            if options.NSinos == options.Nz
                subsets = 180 / options.n_angles;
            else
                error(['Number of sinograms with subset_type = 6 has to be either ' num2str(options.Nz) ' or ' num2str(options.TotSinos)]);
            end
        end
        % Use golden angle sampling
    elseif options.subset_type == 7
        [index, pituus] = goldenAngleSubsets(options);
    end
    if options.precompute_lor
        lor_file = [folder machine_name '_lor_pixel_count_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_sino_' num2str(Ndist) 'x' ...
            num2str(Nang) 'x' num2str(TotSinos) '.mat'];
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
        
        discard = lor > 0;
        if length(discard) ~= TotSinos*Ndist*Nang
            error('Error: Size mismatch between sinogram and LORs to be removed')
        end
        if use_raw_data == false && NSinos ~= TotSinos
            discard = discard(1:NSinos*Ndist*Nang);
        end
        %         lor_a = lor;
        clear lor
        ind_apu = uint32(find(discard));
        
        if iscell(index)
            index = cell2mat(index);
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
    end
elseif subsets > 1
    % For raw data
    if options.precompute_lor
        lor_file = [folder machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'];
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
        LL = form_detector_pairs_raw(options.rings, options.det_per_ring);
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
        % Use only LORs that go through the FOV
        LL = LL(lor > 0,:);
        %         lor_a = lor(lor > 0);
        clear lor
        index = cell(subsets, 1);
        pituus = zeros(subsets, 1, 'int64');
        % Pick the measurements randomly
        if options.subset_type == 3
            indices = uint32(length(LL));
            port = uint32(floor(length(LL)/subsets));
            if options.use_Shuffle
                apu = Shuffle(indices(end), 'index')';
            else
                apu = uint32(randperm(indices(end)))';
            end
            for i = 1 : subsets
                if i == subsets
                    index1 = uint32(apu(port*(i-1)+1:end));
                else
                    index1 = uint32(apu(port*(i-1)+1:(port*(i))));
                end
                index{i} = index1;
                pituus(i) = int64(length(index{i}));
            end
            % Based on the angles of the LORs
        elseif options.subset_type == 6
            [index, pituus] = subset_angles(options, LL);
            subsets = length(pituus);
        else
            % Every nth measurements
            for i=1:subsets
                index1 = uint32(i:subsets:length(LL))';
                index{i} = index1;
                pituus(i) = int64(length(index{i}));
            end
        end
    else
        det_per_ring = options.det_per_ring * options.sampling_raw;
        % Same as above, but for precompute_lor = false case
        LL = form_detector_pairs_raw(options.rings, det_per_ring);
        index = cell(subsets, 1);
        pituus = zeros(subsets, 1, 'int64');
        if options.ring_difference_raw < options.rings
            if options.sampling_raw > 1
                error('Increasing sampling cannot be used with smaller ring difference than the number of rings!')
            end
            testi = zeros(options.detectors,options.detectors,'uint32');
            testi(tril(true(size(testi)), 0)) = uint32(1:length(LL));
            for kk = options.rings : - 1 : options.ring_difference_raw + 1
                for ll = 1 : kk - options.ring_difference_raw
                    testi(1 + (kk - 1) * options.det_per_ring : kk * options.det_per_ring, 1 + (ll - 1) * options.det_per_ring : ll * options.det_per_ring) = ...
                        zeros(options.det_per_ring, options.det_per_ring, 'uint32');
                end
            end
            ind = testi(tril(true(size(testi)), 0));
        end
        if options.subset_type == 3
            indices = uint32(length(LL));
            port = uint32(floor(length(LL)/subsets));
            if options.use_Shuffle && exist('Shuffle','file') == 3
                apu = Shuffle(indices(end), 'index')';
            elseif options.use_Shuffle && exist('Shuffle','file') == 0
                warning('options.Shuffle was set to true, but no Shuffle mex-file found. Using randperm')
                apu = uint32(randperm(indices(end)))';
            else
                apu = uint32(randperm(indices(end)))';
            end
            for i = 1 : subsets
                if i == subsets
                    index1 = uint32(apu(port*(i-1)+1:end));
                else
                    index1 = uint32(apu(port*(i-1)+1:(port*(i))));
                end
                if options.ring_difference_raw < options.rings
                    index1 = index1(ismember(index1,ind));
                end
                index{i} = index1;
                pituus(i) = int64(length(index{i}));
            end
        elseif options.subset_type == 6
            [index, pituus] = subset_angles(options, LL);
            if options.ring_difference_raw < options.rings
                index = index(ismember(index,ind));
            end
            subsets = length(pituus);
        else
            for i=1:subsets
                index1 = uint32(i:subsets:length(LL))';
                if options.ring_difference_raw < options.rings
                    index1 = index1(ismember(index1,ind));
                end
                index{i} = index1;
                pituus(i) = int64(length(index{i}));
            end
        end
    end
end
if ~iscell(index) && size(index,1) == 1 && ~options.precompute_lor
    index = uint32(1 : pituus)';
end