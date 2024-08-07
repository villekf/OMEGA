function distance = ndist_max(options, varargin)
%NDIST_MAX Determine the optimal Ndist values
%   This function can be used to determine the maximum suggested Ndist
%   value as well as an Ndist value that creates a circle that is just
%   inside the FOV.
%
% Examples:
%   Ndist = ndist_max(options)
%   Ndist = ndist_max(options, pseudot)
%
% INPUTS:
%   options = The options struct created by any of main-files, only the
%   machine properties and image properties sections are required.
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

if nargin >= 2 && ~isempty(varargin{1})
    options.pseudot = varargin{1};
else
    options.pseudot = 0;
end
if nargin >= 3 && ~isempty(varargin{2})
    verbose = varargin{2};
else
    verbose = true;
end
if ~isfield(options,'implementation')
    options.implementation = 4;
end
options.nLayers = 1;
if nargin > 1 && ~isempty(options.pseudot) && options.pseudot > 0
    [~, ~, x, y] = detector_coordinates(options);
else
    [x, y] = detector_coordinates(options);
end
y1 = max(y);
x1 = x(y == y1);
if numel(x1) > 1
    x1 = x1(floor(length(x1)/2));
end
x0 = max(x)/2;
y0 = max(y)/2;

distance = sort((abs((y-y1)*x0 - (x - x1)*y0 + x.*y1 - y.*x1)./sqrt((y-y1).^2 + (x-x1).^2)));

dist = numel((distance(distance<(options.FOVa_x/sqrt(2)))));

if verbose
    disp(['Maximum suggested Ndist value is ' num2str(dist)])
end

dist = numel((distance(distance<=(options.FOVa_x/2))));

if verbose
    disp(['Circle inside the FOV is achieved with Ndist value of ' num2str(dist)])
end