function [index, pituus, subsets] = index_maker(options)
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
% Copyright (C) 2020-2024 Ville-Veikko Wettenhovi
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


if ~isfield(options,'sampling_raw')
    options.sampling_raw = 1;
end
subsets = options.subsets;
Ndist = options.Ndist;
Nang = options.Nang;
NSinos = options.NSinos;
if options.use_raw_data
    pituus = int64(options.detectors^2/2 + options.detectors/2);
else
    pituus = int64(Ndist * Nang * NSinos);
end
index = 0;
if options.CT
    tyyppi = 'uint64';
else
    tyyppi = 'uint32';
end
% Subset data
if subsets > 1 && options.subset_type < 8
    if options.use_raw_data
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
        Ndist = options.det_per_ring;
        Nang = options.det_per_ring;
        NSinos = sum(options.ring_difference_raw : - 1 : 1);
    end
    totalLength = Ndist*Nang*NSinos;
    index = cell(subsets,1);
    pituus = zeros(subsets, 1, 'int64');
    % Take every nth column from the sinogram
    if options.subset_type == 4 && ~options.use_raw_data
        maksimi = numel(1 : subsets : Nang);
        for i=1:subsets
            koko = i:subsets:Nang;
            osa = length(koko);
            index1 = repmat(cast((i-1)*Ndist+1:i*Ndist, tyyppi)', osa*NSinos,1);
            if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                index1 = index1 + cast(repeat_elem((0:(osa)*NSinos-1)'*Ndist*subsets,Ndist), tyyppi);
            else
                index1 = index1 + cast(repelem((0:(osa)*NSinos-1)'*Ndist*subsets,Ndist), tyyppi);
            end
            if mod(Nang,subsets) > 0
                if osa < maksimi
                    erotus = osa - 1;
                else
                    erotus = mod(Nang,subsets) - subsets;
                end
                index1 = cast((int64(index1) + int64(repelem((0:NSinos-1)'*Ndist*erotus,Ndist*osa))), tyyppi);
                % ero = [zeros(Ndist*osa, 1, 'int64');repelem(Ndist*erotus, numel(index1) - Ndist*osa)'];
                % index1 = cast(int64(index1) + ero, tyyppi);
            end
            index{i} = index1;
            pituus(i) = int64(length(index{i}));
        end
        % Take every nth row from the sinogram
    elseif options.subset_type == 5
        apu = repmat(repmat(cast((1 : Ndist : Nang * Ndist), tyyppi), Ndist, 1) + cast((0 : Ndist - 1)', tyyppi), 1, 1, NSinos) + permute(cast((0 : Nang * Ndist : Nang * Ndist * (NSinos - 1))', tyyppi),[3 2 1]);
        apu = permute(apu, [1 3 2]);
        apu = reshape(apu, Ndist * NSinos, []);
        for i=1:subsets
            index1 = apu((i - 1) + 1 : subsets : end, :);
            index1 = index1';
            index1 = index1(:);
            index{i} = index1;
            pituus(i) = int64(length(index{i}));
        end
        % Take every nth (column) measurement
    elseif options.subset_type == 1
        for i=1:subsets
            index1 = cast(i:subsets:totalLength, tyyppi)';
            [I,J,K] = ind2sub([Nang Ndist NSinos], index1);
            index1 = cast(sub2ind([Ndist Nang NSinos], J, I,K), tyyppi);
            index{i} = index1;
            pituus(i) = int64(length(index{i}));
        end
        % Take every nth (row) measurement
        % Every nth measurements
    elseif options.subset_type == 2
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
    elseif options.subset_type == 0
        if options.listmode == 1
            val = floor(totalLength / subsets);
            if mod(totalLength, subsets) > 0
                valEnd = totalLength - val * (subsets - 1);
            else
                valEnd = val;
            end
        else
            val = floor(options.nProjections / subsets);
            if mod(options.nProjections, subsets) > 0
                valEnd = options.nProjections - val * (subsets - 1);
            else
                valEnd = val;
            end
        end
        for i = 1 : subsets - 1
            pituus(i) = val;
        end
        pituus(subsets) = valEnd;
        index = {0};
    end
elseif (subsets > 1 && (options.subset_type == 8 || options.subset_type == 9 || options.subset_type == 10 || options.subset_type == 11 || options.subset_type == 12)) || (subsets == 1)
    sProjections = floor(options.nProjections / subsets);
    pituus = zeros(subsets, 1, 'int64');
    index = cell(subsets, 1);
    modi = mod(options.nProjections, subsets);
    uu = double(modi > 0);
    ind1 = 1;
    ind2 = sProjections + uu;
    if options.subset_type == 8 && subsets > 1
        for kk = 1 : subsets
            index{kk} = (kk : subsets : options.nProjections)';
            pituus(kk) = numel(index{kk});
        end
    elseif options.subset_type == 9 && subsets > 1
        apu = randperm(options.nProjections)';
        for kk = 1 : subsets
            index{kk} = apu(ind1 : ind2);
            pituus(kk) = numel(index{kk});
            modi = modi - 1;
            if modi <= 0
                uu = 0;
            end
            ind1 = ind2 + 1;
            ind2 = ind2 + (floor(options.nProjections / subsets)) + uu;
        end
    elseif options.subset_type == 10 && subsets > 1
        if ~isfield(options,'angles')
            error('Subset type 10 is only available with CT data and when the projections angles (options.angles) are supplied!')
        end
        ga = 2.39996322972865332;
        options.angles = abs(options.angles);
        angles = options.angles - min(options.angles);
        anglesOrig = angles;
        maksimi = max(angles);
        ga = ga * (maksimi / (2*pi));
        angle = 0;
        for kk = 1 : subsets - 1
            ind = zeros(floor(numel(options.angles) / subsets),1);
            for ii = 1 : floor(numel(options.angles) / subsets)
                [~,I] = min(abs(angles-angle));
                II = find(ismember(anglesOrig,angles(I)));
                angles(I) = [];
                ind(ii) = II;
                angle = angle + ga;
                if angle > maksimi
                    angle = angle - maksimi;
                end
            end
            index{kk} = ind;
            pituus(kk) = numel(index{kk});
        end
        ind = zeros(ceil(numel(options.angles) / subsets),1);
        for ii = 1 : ceil(numel(options.angles) / subsets)
            [~,I] = min(abs(angles-angle));
            II = find(ismember(anglesOrig,angles(I)));
            angles(I) = [];
            ind(ii) = II;
            angle = angle + ga;
            if angle > maksimi
                angle = angle - maksimi;
            end
        end
        index{subsets} = ind;
        pituus(subsets) = numel(index{subsets});
    elseif options.subset_type == 11 && subsets > 1
        v = powerFactorSubsets(options.nProjections);
        for kk = 1 : subsets
            index{kk} = v(ind1 : ind2);
            pituus(kk) = numel(index{kk});
            modi = modi - 1;
            if modi <= 0
                uu = 0;
            end
            ind1 = ind2 + 1;
            ind2 = ind2 + (floor(options.nProjections / subsets)) + uu;
        end
    elseif options.subset_type == 12 && subsets > 1
        for kk = 1 : subsets
            index{kk} = [(kk : subsets : options.nProjections/2);(kk + options.nProjections/2 : subsets : options.nProjections)];
            index{kk} = index{kk}(:);
            pituus(kk) = numel(index{kk});
        end
    else
        if options.use_raw_data
            pituus = int64(options.detectors^2/2 + options.detectors/2);
            index{1} = uint32(1 : pituus)';
        else
            index{1} = (1 : options.nProjections)';
        end
        pituus(1) = numel(index{1});
    end
    if options.subsets == 1
        if (options.CT || options.PET || options.SPECT) && options.listmode == 0
            pituus(1) = options.NSinos;
        elseif options.listmode
            if options.useIndexBasedReconstruction
                pituus(1) = numel(options.trIndex) / 2;
            else
                pituus(1) = numel(options.x) / 6;
            end
            index = {0};
        else
            pituus(1) = options.Nang * options.Ndist * options.NSinos;
        end
    end
    index = cell2mat(index);
elseif options.subset_type > 11
    error('Invalid subset type!')
end
if ~iscell(index) && size(index,1) == 1 && ~options.precompute_lor
    index = 0;
end