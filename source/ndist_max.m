function distance = ndist_max(cr_p, diameter, cryst_per_block, blocks_per_ring, FOVa_x, varargin)
%NDIST_MAX Determine the optimal Ndist values
%   This function can be used to determine the maximum suggested Ndist
%   value as well as an Ndist value that creates a circle that is just
%   inside the FOV.
%
% Examples:
%   Ndist = ndist_max(cr_p, diameter, cryst_per_block, blocks_per_ring,
%   FOVa_x)
%   Ndist = ndist_max(cr_p, diameter, cryst_per_block, blocks_per_ring,
%   FOVa_x, pseudot)
%
% INPUTS:
%   cr_p = Crystal pitch in transaxial direction
%   diameter = Diameter of the system (bore)
%   cryst_per_block = Transaxial crystal count per block
%   blocks_per_ring = The number of blocks/buckets per ring
%   FOVa_x = FOV size in transaxial direction (only square FOVs are
%   supported)
%   pseudot = (optional) The number of pseudo rings. Should be used only
%   when using a system with pseudo detectors.
%   verbose = Turn of the messages; output only the distances.
%
% OUTPUTS:
%   dist = The orthogonal distance of a line joining two detectors and the
%   origin. The vector has been sorted from smallest distance to the
%   longest (with the last one being NaN, as both the detectors are
%   actually the same one).
%
% See also sinogram_coordinates_2D, detector_coordinates

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

options.cr_p = cr_p;
options.diameter = diameter;
options.cryst_per_block = cryst_per_block;
options.blocks_per_ring = blocks_per_ring;
options.offangle = 0;
options.flip_image = false;
if nargin >= 6 && ~isempty(varargin{1})
    options.pseudot = varargin{1};
else
    options.pseudot = 0;
end
if nargin >= 7 && ~isempty(varargin{2})
    verbose = varargin{2};
else
    verbose = true;
end
if nargin > 5 && ~isempty(options.pseudot) && options.pseudot > 0
    [~, ~, x, y] = detector_coordinates(options);
else
    [x, y] = detector_coordinates(options);
end
y1 = max(y);
x1 = x(y == y1);
x1 = x1(floor(length(x1)/2));
x0 = max(x)/2;
y0 = max(y)/2;

distance = sort((abs((y-y1)*x0 - (x - x1)*y0 + x.*y1 - y.*x1)./sqrt((y-y1).^2 + (x-x1).^2)));

dist = numel((distance(distance<(FOVa_x/2*sqrt(2)))));

if verbose
    disp(['Maximum suggested Ndist value is ' num2str(dist)])
end

dist = numel((distance(distance<=(FOVa_x/2))));

if verbose
    disp(['Circle inside the FOV is achieved with Ndist value of ' num2str(dist)])
end