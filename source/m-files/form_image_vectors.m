function im_vectors = form_image_vectors(options, N, varargin)
%FORM_IMAGE_VECTORS This function simply forms the applicable image vectors
%for each volume
%
%
%
% Copyright (C) 2019-2024 Ville-Veikko Wettenhovi
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

if options.save_iter
    Niter = options.Niter + 1;
else
    Niter = max(1, numel(options.saveNIter));
end
if (isfield(options,'useSingles') && options.useSingles) || options.implementation == 5
    type = 'single';
else
    type = 'double';
end
if nargin > 2 && ~isempty(varargin{1})
    noSensIm = varargin{1};
else
    noSensIm = true;
end
if numel(options.partitions) > 1
    partitions = numel(options.partitions);
else
    partitions = options.partitions;
end

im_vectors.recImage = cell(options.nMultiVolumes + 1, 1); 
im_vectors.recApu = cell(partitions, options.nMultiVolumes + 1); % Holds current estimate
im_vectors.rhs = cell(partitions, options.nMultiVolumes + 1); % Holds current backprojections
im_vectors.Sens = cell(partitions, options.nMultiVolumes + 1); % Holds sensitivity image
if ~noSensIm
    for timestep = 1:partitions
        for kk = 1 : options.nMultiVolumes + 1
            if options.saveSens
                im_vectors.Sens{timestep, kk} = zeros(N(kk), options.subsets, type);
            else
                im_vectors.Sens{timestep, kk} = zeros(N(kk), 1, type);
            end
        end
    end
end
for kk = 1 : options.nMultiVolumes + 1
    im_vectors.recImage{kk} = ones(N(kk), Niter, partitions, type);
end

if options.save_iter
    for tt = 1:partitions
        im_vectors.recImage{1}(:,1,tt) = cast(options.x0(:), type);
        if options.nMultiVolumes >= 2
            im_vectors.recImage{2}(:,1,tt) = cast(options.x1(:), type);
            im_vectors.recImage{3}(:,1,tt) = cast(options.x2(:), type);
        end
        if options.nMultiVolumes >= 4
            im_vectors.recImage{4}(:,1,tt) = cast(options.x3(:), type);
            im_vectors.recImage{5}(:,1,tt) = cast(options.x4(:), type);
        end
        if options.nMultiVolumes == 6
            im_vectors.recImage{6}(:,1,tt) = cast(options.x5(:), type);
            im_vectors.recImage{7}(:,1,tt) = cast(options.x6(:), type);
        end
    end
end

for timestep = 1:partitions
    im_vectors.recApu{timestep, 1} = cast(options.x0(:), type);
    if options.nMultiVolumes >= 2
        im_vectors.recApu{timestep, 2} = cast(options.x1(:), type);
        im_vectors.recApu{timestep, 3} = cast(options.x2(:), type);
    end
    if options.nMultiVolumes >= 4
        im_vectors.recApu{timestep, 4} = cast(options.x3(:), type);
        im_vectors.recApu{timestep, 5} = cast(options.x4(:), type);
    end
    if options.nMultiVolumes == 6
        im_vectors.recApu{timestep, 6} = cast(options.x5(:), type);
        im_vectors.recApu{timestep, 7} = cast(options.x6(:), type);
    end
end

% Special ECOSEM case guarantees that OSEM and COSEM are initialized even
% if they haven't been selected
if options.ECOSEM
    if ~options.OSEM
        im_vectors.OSEMApu = cast(options.x0(:), type);
    end
    if ~options.COSEM
        im_vectors.COSEMApu = cast(options.x0(:), type);
    end
end