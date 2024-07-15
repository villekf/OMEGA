function [output, varargout] = computeIntegralImage(input, varargin)
% COMPUTEINTEGRALIMAGE Computes the integral image for the input data
%   Outputs the x-, y- and z-coordinates for the current machine, depending
%   on whether sinogram or raw data is used
%
% EXAMPLES:
%   [x, y, z, options] = get_coordinates(options)
%   [x, y, z, options] = get_coordinates(options, rings, pseudot)
% INPUTS:
%   options = Machine properties, sinogram properties and whether raw data
%   is used are needed.
%   rings = Number of crystal rings (required only for raw data)
%   pseudot = The numbers of the pseudo rings (required only for raw data)
%
% See also sinogram_coordinates_2D, sinogram_coordinates_3D,
% detector_coordinates, getMultirayCoordinates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2023-2024 Ville-Veikko Wettenhovi
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

subtractMean = false;
if nargin >= 2
    subtractMean = varargin{1};
end
type = class(input);
input = double(input);
output = zeros(size(input,1) + 1, size(input,2) + 1, size(input,3));
if subtractMean
    meanVal = mean(mean(input, 1), 2);
    input = bsxfun(@minus, input, meanVal);
    if nargout == 1
        error('Mean subtracted, but not output. Set two output values when subtracting mean!')
    end
    varargout{1} = meanVal;
else
    if nargout >= 2
        varargout{1} = [];
    end
end
for kk = 1 : size(input,3)
    output(2:end,2:end,kk) = cumsum(cumsum(input(:,:,kk),1),2);
end
output = cast(output(:), type);