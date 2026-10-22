function outputVol = padMultiResVolume(inputVol, options, fillZero)
%PADMULTIRESVOLUME Extends a whole-grid volume by the enlargement of the extended FOV
%   The side volumes of the multi-resolution reconstruction have exactly the coarse voxel size, so their lengths are rounded up to
%   a whole number of coarse voxels and the extended FOV grows by less than one coarse voxel on each outer face (see setUpCorrections).
%   A whole-grid volume (options.NxFull x NyFull x NzFull voxels, the requested extended FOV) is extended here by the enlargement,
%   options.multiResPadLow/High (in voxels of the main volume), so that it covers the FOV of the volumes: with zeros if fillZero is
%   true, otherwise with the values of the edge voxels. A volume of another size, or an enlargement that is not a whole number of voxels,
%   gives the input volume as it is.
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

outputVol = inputVol;
if ~isfield(options, 'multiResPadExact') || ~options.multiResPadExact
    return
end
lo = round(options.multiResPadLow);
hi = round(options.multiResPadHigh);
sz = [size(inputVol, 1), size(inputVol, 2), size(inputVol, 3)];
if (~any(lo) && ~any(hi)) || ~isequal(sz, double([options.NxFull, options.NyFull, options.NzFull]))
    return
end
if fillZero
    outputVol = zeros(sz + lo + hi, 'like', inputVol);
    outputVol(lo(1) + (1 : sz(1)), lo(2) + (1 : sz(2)), lo(3) + (1 : sz(3))) = inputVol;
else
    outputVol = inputVol([ones(1, lo(1)), 1 : sz(1), sz(1) * ones(1, hi(1))], [ones(1, lo(2)), 1 : sz(2), sz(2) * ones(1, hi(2))], ...
        [ones(1, lo(3)), 1 : sz(3), sz(3) * ones(1, hi(3))]);
end
end
