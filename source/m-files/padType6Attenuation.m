function options = padType6Attenuation(options)
%PADTYPE6ATTENUATION Pads the type 6 attenuation map to the union of the multi-resolution volumes
%   The attenuation map of the SPECT projector type 6 with multi-resolution volumes covers the union bounding box of the volumes
%   (options.multiResUnionDims voxels of the main volume), i.e. the requested extended FOV (options.NxFull x NyFull x NzFull) and the
%   enlargement of it by the side volumes (their lengths are rounded up to a whole number of coarse voxels, see setUpCorrections). A
%   map of the requested extended FOV is padded here by the enlargement; a map that already is of the size of the union is left as it
%   is. A dynamic map is the frames one after the other, or a cell array of the maps.
%   The padding is air (mu = 0) transaxially (x and y), since the body is inside the transaxial FOV, but axially (z) the edge slice of
%   the map is replicated: the body continues beyond the axial FOV, and an air layer would make the lowest coarse slice of the side
%   volume an average of air and tissue, i.e. the activity near the axial edge would be under-attenuated (up to a few percent 20 mm
%   from the edge and far more at the edge itself). With no axial enlargement the map is not changed at all.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021-2025 Ville-Veikko Wettenhovi
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(options, 'multiResUnionDims') || (isfield(options, 'multiResAttenuationPadded') && options.multiResAttenuationPadded)
    return
end
oldDims = double([options.NxFull, options.NyFull, options.NzFull]);
newDims = double(options.multiResUnionDims(:).');
if isequal(oldDims, newDims)
    return
end
lo = max(floor(double(options.multiResPadLow(:).') + 0.5), 0);
hi = max(newDims - oldDims - lo, 0);
if isa(options.vaimennus, 'cell')
    for kk = 1 : numel(options.vaimennus)
        options.vaimennus{kk} = padOneAttenuation(options.vaimennus{kk}, oldDims, newDims, lo, hi);
    end
else
    options.vaimennus = padOneAttenuation(options.vaimennus, oldDims, newDims, lo, hi);
end
options.multiResAttenuationPadded = true;
end

function out = padOneAttenuation(vol, oldDims, newDims, lo, hi)
out = vol;
nVox = numel(vol);
if nVox == 0
    return
end
if mod(nVox, prod(oldDims)) ~= 0
    if mod(nVox, prod(newDims)) ~= 0
        warning('The size of the type 6 attenuation map (%d) is not that of the requested extended FOV (%d x %d x %d voxels) or of the volumes, it is used as it is.', nVox, oldDims(1), oldDims(2), oldDims(3))
    end
    return % Already the size of the union (or unknown)
end
isVec = isvector(vol);
nFrames = nVox / prod(oldDims);
vol = reshape(vol, [oldDims, nFrames]);
out = zeros([newDims, nFrames], 'like', vol);
% Air in x and y, the nearest edge slice in z (the identity index when there is no axial padding)
zIdx = min(max((1 : newDims(3)) - lo(3), 1), oldDims(3));
out(lo(1) + (1 : oldDims(1)), lo(2) + (1 : oldDims(2)), :, :) = vol(:, :, zIdx, :);
if isVec
    out = out(:);
end
end
