function [options, lor_a, xy_index, z_index, LL, summa, index] = form_subset_indices(options, pituus, index, y, varargin)
%% FORM SUBSET INDICES
% This function is used to form the necessary index variables for the
% reconstruction depending on the number of subsets and whether sinogram or
% raw list-mode data is used.
% INPUTS:
%   options = Contains the machine and (optionally) sinogram specific
%   information (e.g. detectors per ring, image (matrix) size, FOV size).
%   pituus = The number of measurements (LORs) that each subset has.
%   options.subsets = The number of options.subsets used.
%   index = The indices corresponding to each subset.
%   size_x = Length of the detector coordinate matrix.
%   y = Detector coordinates in the y-axis.
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

lor_a = uint16([]);
LL = uint16([]);
if nargin >= 5 && ~isempty(varargin{1}) && options.precompute_lor
    lor_a = varargin{1};
end
if nargin >= 6 && ~isempty(varargin{2}) && options.use_raw_data
    LL = varargin{2};
end
if options.subset_type >= 8 || options.subsets == 1
    kerroin = options.nRowsD * options.nColsD;
else
    kerroin = 1;
end
% if nargin >= 14 && ~isempty(varargin{5}) && options.precompute_lor
%     storeMatrix = varargin{5};
% else
%     storeMatrix = false;
% end
if options.sampling > 1 && ~options.use_raw_data && ~options.precompute_lor
    options.Ndist = options.Ndist * options.sampling;
end
if ~isfield(options,'sampling_raw')
    options.sampling_raw = 1;
end
if ~isfield(options,'CT')
    options.CT = false;
end
if ~isfield(options,'SPECT')
    options.SPECT = false;
end
if ~options.precompute_lor
    if options.NSinos ~= options.TotSinos && options.listmode == 0
        index = index(1:options.Ndist*options.Nang*options.NSinos);
    end
end
if (options.subset_type < 8 && options.subsets > 1) || options.subsets == 1
    if options.use_raw_data 
%         LL = form_detector_pairs_raw(rings, options.det_per_ring);
        LL = LL(index,:);
    end
end
if options.subsets > 1 && options.precompute_lor
    if options.subset_type >= 8
        lor_a = reshape(lor_a, options.Ndist, options.Nang, []);
        lor_a = lor_a(:,:,index);
    else
        lor_a = (lor_a(index));
    end
end
if options.listmode > 0 && options.subset_type < 8 && options.subsets > 1
    if options.subset_type > 0
        if options.useIndexBasedReconstruction
            if size(options.trIndex, 1) ~= 2
                options.trIndex = reshape(options.trIndex, 2, []);
            end
            options.trIndex = options.trIndex(:,index);
            if size(options.axIndex, 1) ~= 2
                options.axIndex = reshape(options.axIndex, 2, []);
            end
            options.axIndex = options.axIndex(:,index);
        else
            if size(options.x, 1) ~= 6
                options.x = reshape(options.x, 6, []);
            end
            options.x = options.x(:,index);
        end
    end
    options.x = options.x(:);
    % options.y = options.y(index,:);
    % if isfield(options, 'z_det')
    %     options.z_det = options.z_det(index,:);
    % else
    %     options.z = options.z(index,:);
    % end
    xy_index = uint32(0);
    z_index = uint16(0);
elseif options.listmode > 0 && (options.subset_type == 8 || options.subset_type == 9) && options.subsets > 1
    if isfield(options, 'z_det')
        options.z_det = options.z_det(index,:);
    else
        options.z = options.z(index,:);
    end
    if isfield(options, 'uV') && size(options.uV,1) == options.nProjections
        options.uV = options.uV(index,:);
    end
    xy_index = uint32(0);
    z_index = uint16(0);
elseif ~options.use_raw_data && ((options.subsets > 1 && (options.subset_type == 3 || options.subset_type == 6 || options.subset_type == 7)))
    xy_index = uint32(1:options.Nang * options.Ndist)';
    if options.span > 1
        xy_index2 = repmat(uint32(1:options.Nang * options.Ndist)', options.NSinos - (options.rings * 2 - 1), 1);
        xy_index = [repmat(xy_index,  options.rings * 2 - 1, 1); xy_index2];
    else
        xy_index2 = repmat(uint32(1:options.Nang * options.Ndist)', options.NSinos - options.rings, 1);
        xy_index = [repmat(xy_index, options.rings, 1); xy_index2];
    end
    z_index = uint16(1:options.NSinos)';
    if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
        z_index = repeat_elem(z_index, options.Nang * options.Ndist);
    else
        z_index = repelem(z_index, options.Nang * options.Ndist);
    end
    z_index = z_index(index);
    z_index = z_index - 1;

    xy_index = xy_index(index);
    xy_index = xy_index - 1;
else
    xy_index = uint32(0);
    z_index = uint16(0);
end

summa = zeros(options.subsets, 1, 'uint64');

if ~options.use_raw_data
    if (options.subsets > 1 && length(pituus) > 1)
        if options.precompute_lor
            for kk = 1 : options.subsets
                if options.subset_type >= 8
                    summa(kk) = uint64(sum(sum(sum(uint64(lor_a(:,:,pituus(kk) + 1 : pituus(kk + 1)))))));
                else
                    summa(kk) = uint64(sum(uint64(lor_a(pituus(kk) * kerroin + 1 : pituus(kk + 1) * kerroin))));
                end
            end
        else
            summa = uint64(pituus);
        end
    else
        if options.precompute_lor
            summa = uint64(sum(int64(lor_a)));
            pituus = int64([0;length(lor_a)]);
        else
            if options.PET || options.CT
                pituus = int64([0;options.NSinos]);
            else
                pituus = int64([0;options.Nang*options.Ndist*options.NSinos]);
            end
        end
    end
end
if options.use_raw_data
    if options.subsets > 1
        for kk = 1 : options.subsets
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
            if options.precompute_lor
                    summa(kk) = uint64(sum(int64(lor_a(pituus(kk)+1:pituus(kk+1)))));
            end
        end
    else
        apu = LL - 1;
        apu2 = idivide(apu, uint16(options.det_per_ring));
        idx = apu2(:,1) == apu2(:,2);
        apu2 = apu(idx,:);
        ind = mod(apu2, uint16(options.det_per_ring)) + 1;
        yt = y(ind);
        y_i = yt(:,1) > yt(:,2);
        apu2(y_i,:) = fliplr(apu2(y_i,:));
        apu(idx,:) = apu2;
        LL = apu + 1;
        if options.precompute_lor
                summa = uint64(sum(int64(lor_a)));
        end
    end
    LL = LL';
    LL = LL(:);
end


if options.sampling > 1 && ~options.use_raw_data && ~options.precompute_lor
    options.Ndist = options.Ndist / options.sampling;
end
