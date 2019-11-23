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
% Copyright (C) 2019  Ville-Veikko Wettenhovi
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
if nargin - 7 == 4
    Nang = varargin{1};
    Ndist = varargin{2};
    TotSinos = varargin{3};
    NSinos = varargin{4};
end

folder = fileparts(which('index_maker.m'));
folder = strrep(folder, 'source','mat-files/');
folder = strrep(folder, '\','/');
% lor_a = 0;
if options.use_raw_data
    pituus = options.detectors ^2/2 + options.detectors/2;
else
    pituus = Ndist * Nang * NSinos;
end
index = 0;
% Sinogram data
if use_raw_data == false && subsets > 1
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
        index = cell(subsets, 1);
        pituus = zeros(subsets, 1, 'uint32');
        % Take every nth column from the sinogram
        if options.subset_type == 4
            for i=1:subsets
                osa = length(i-1:subsets:Ndist);
                if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                    index1 = uint32(repmat(repmat((1:Ndist)', osa,1) + repeat_elem((i-1:subsets:Ndist)'*Ndist,Ndist), NSinos, 1) + repeat_elem(Ndist*Nang*(0:NSinos-1)',...
                        Ndist*osa));
                else
                    index1 = uint32(repmat(repmat((1:Ndist)', osa,1) + repelem((i-1:subsets:Ndist)'*Ndist,Ndist), NSinos, 1) + repelem(Ndist*Nang*(0:NSinos-1)',...
                        Ndist*osa));
                end
                index{i} = index1;
                pituus(i) = uint32(length(index{i}));
            end
        % Take every nth row from the sinogram
        elseif options.subset_type == 5
            for i=1:subsets
                osa = length(i-1:subsets:Nang);
                if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                    index1 = uint32(repmat(repmat((1:Ndist:Ndist*Nang)', osa,1) + repeat_elem((i-1:subsets:Nang)',Nang), NSinos, 1) + ...
                        repeat_elem(Ndist*Nang*(0:NSinos-1)',Nang*osa));
                else
                    index1 = uint32(repmat(repmat((1:Ndist:Ndist*Nang)', osa,1) + repelem((i-1:subsets:Nang)',Nang), NSinos, 1) + ...
                        repelem(Ndist*Nang*(0:NSinos-1)',Nang*osa));
                end
                index{i} = index1;
                pituus(i) = uint32(length(index{i}));
            end
        % Take every nth (column) measurement
        elseif options.subset_type == 1
            for i=1:subsets
                index1 = uint32(i:subsets:Ndist*Nang*NSinos)';
                [I,J,K] = ind2sub([Nang Ndist NSinos], index1);
                index1 = uint32(sub2ind([Ndist Nang NSinos], J, I,K));
                index{i} = index1;
                pituus(i) = uint32(length(index{i}));
            end
        % Take every nth (row) measurement
        elseif options.subset_type == 2
            for i=1:subsets
                index1 = uint32(i:subsets:Ndist*Nang*NSinos)';
                index{i} = index1;
                pituus(i) = uint32(length(index{i}));
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
                pituus(i) = uint32(length(index{i}));
            end
        % Pick the subsets based on the angles of the LORs
        elseif options.subset_type == 6
            [index, pituus] = subset_angles(options);
            subsets = length(pituus);
        % Use golden angle sampling
        elseif options.subset_type == 7
            [index, pituus] = goldenAngleSubsets(options);
        end
        if iscell(index)
            index = cell2mat(index);
        end
        % Remove indices that do not go through the FOV
        joku = ismember(index, ind_apu);
        pituus2 = cumsum(pituus);
        for kk = 1 : length(pituus2)
            if kk == 1
                pituus(kk) = sum(joku(1:pituus2(kk)));
            else
                pituus(kk) = sum(joku(1+pituus2(kk-1) : pituus2(kk)));
            end
        end
        index = index(joku);
    else
        % Same as above, but for precompute_lor = false case
        index = cell(subsets,1);
        pituus = zeros(subsets, 1, 'uint32');
        if options.subset_type == 4
            for i=1:subsets
                osa = length(i-1:subsets:Ndist);
                if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                    index1 = uint32(repmat(repmat((1:Ndist)', osa,1) + repeat_elem((i-1:subsets:Ndist)'*Ndist,Ndist), NSinos, 1) + repeat_elem(Ndist*Nang*(0:NSinos-1)',...
                        Ndist*osa));
                else
                    index1 = uint32(repmat(repmat((1:Ndist)', osa,1) + repelem((i-1:subsets:Ndist)'*Ndist,Ndist), NSinos, 1) + repelem(Ndist*Nang*(0:NSinos-1)',...
                        Ndist*osa));
                end
                index{i} = index1;
                pituus(i) = uint32(length(index{i}));
            end
        elseif options.subset_type == 5
            for i=1:subsets
                osa = length(i-1:subsets:Nang);
                if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                    index1 = uint32(repmat(repmat((1:Ndist:Ndist*Nang)', osa,1) + repeat_elem((i-1:subsets:Nang)',Nang), NSinos, 1) + ...
                        repeat_elem(Ndist*Nang*(0:NSinos-1)',Nang*osa));
                else
                    index1 = uint32(repmat(repmat((1:Ndist:Ndist*Nang)', osa,1) + repelem((i-1:subsets:Nang)',Nang), NSinos, 1) + ...
                        repelem(Ndist*Nang*(0:NSinos-1)',Nang*osa));
                end
                index{i} = index1;
                pituus(i) = uint32(length(index{i}));
            end
        elseif options.subset_type == 1
            for i=1:subsets
                index1 = uint32(i:subsets:Ndist*Nang*NSinos)';
                [I,J,K] = ind2sub([Nang Ndist NSinos], index1);
                index1 = uint32(sub2ind([Ndist Nang NSinos], J, I,K));
                index{i} = index1;
                pituus(i) = uint32(length(index{i}));
            end
        elseif options.subset_type == 2
            for i=1:subsets
                index1 = uint32(i:subsets:Ndist*Nang*NSinos)';
                index{i} = index1;
                pituus(i) = uint32(length(index{i}));
            end
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
                pituus(i) = uint32(length(index{i}));
            end
        % Pick the subsets based on the angles of the LORs
        elseif options.subset_type == 6
            [index, pituus] = subset_angles(options);
            subsets = length(pituus);
        % Use golden angle sampling
        elseif options.subset_type == 7
            [index, pituus] = goldenAngleSubsets(options);
        end
    end
elseif subsets > 1
    % For raw list-mode data
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
        % Use only LORs that go through the FOV
        LL = LL(lor > 0,:);
%         lor_a = lor(lor > 0);
        clear lor
        index = cell(subsets, 1);
        pituus = zeros(subsets, 1, 'uint32');
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
                pituus(i) = uint32(length(index{i}));
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
                pituus(i) = uint32(length(index{i}));
            end
        end
    else
        % Same as above, but for precompute_lor = false case
        LL = form_detector_pairs_raw(options.rings, options.det_per_ring);
        index = cell(subsets, 1);
        pituus = zeros(subsets, 1, 'uint32');
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
                index{i} = index1;
                pituus(i) = uint32(length(index{i}));
            end
        elseif options.subset_type == 6
            [index, pituus] = subset_angles(options, LL);
            subsets = length(pituus);
        else
            for i=1:subsets
                index1 = uint32(i:subsets:length(LL))';
                index{i} = index1;
                pituus(i) = uint32(length(index{i}));
            end
        end
    end
end