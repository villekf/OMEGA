function out = warpToReference(vol, refIn, refOut, fillValue)
%WARPTOREFERENCE Resample a volume from one imref3d grid to another
%   OUT = WARPTOREFERENCE(VOL, REFIN, REFOUT, FILLVALUE) resamples the
%   volume VOL, defined on the imref3d grid REFIN, onto the grid REFOUT
%   with trilinear interpolation (identity transformation). Voxels of REFOUT
%   outside REFIN are set to FILLVALUE.
%
%   If imwarp and affinetform3d are available they are used. Otherwise a
%   fallback based on interpn is used. In the fallback the volume is
%   plain linear interpolation is used inside the input voxel-centre range
%   and output voxels whose centres fall outside that range get FILLVALUE,
%   as imwarp does. The result is single if VOL is single and double
%   otherwise.
%
% See also imwarp, imref3d, interpn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2026 Ville-Veikko Wettenhovi
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
if exist('imwarp', 'file') == 2 && exist('affinetform3d') ~= 0
    out = imwarp(vol, refIn, affinetform3d(eye(4)), 'OutputView', refOut, 'FillValue', fillValue, 'InterpolationMethod', 'linear');
    return
end

% imref3d: dimension 1 (rows) runs along Y, dimension 2 (columns) along X
% and dimension 3 along Z
limIn = {refIn.YWorldLimits, refIn.XWorldLimits, refIn.ZWorldLimits};
limOut = {refOut.YWorldLimits, refOut.XWorldLimits, refOut.ZWorldLimits};
szIn = refIn.ImageSize;
szOut = refOut.ImageSize;

c = cell(1, 3);
for d = 1 : 3
    extIn = (limIn{d}(2) - limIn{d}(1)) / szIn(d);
    extOut = (limOut{d}(2) - limOut{d}(1)) / szOut(d);
    w = limOut{d}(1) + ((1 : szOut(d)) - 0.5) * extOut;
    % Intrinsic coordinates of the output voxel centres in the input grid
    c{d} = (w - limIn{d}(1)) / extIn + 0.5;
end

% Output voxels whose centres map outside the input voxel-centre range
% [1, n] get fillValue (extrapolation value of interpn)
isSingle = isa(vol, 'single');
[C1, C2, C3] = ndgrid(c{1}, c{2}, c{3});
out = interpn(double(vol), C1, C2, C3, 'linear', double(fillValue));
out = reshape(out, szOut);
if isSingle
    out = single(out);
end
end
