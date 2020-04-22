function [x,y] = getMultirayCoordinates(options)
%GETMULTIRAYCOORDINATES Computes the detector coordinates for multi-ray
%Siddon
%   Outputs the x- and y-coordinates for the current machine, depending
%   on whether sinogram or raw data is used. Multi-ray case.
%
% EXAMPLES:
%   [x,y] = getMultirayCoordinates(options)
% INPUTS:
%   options = Machine properties, sinogram properties and whether raw data
%   is used are needed.
%
% See also sinogram_coordinates_2D, sinogram_coordinates_3D,
% detector_coordinates, get_coordinates

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


if options.use_raw_data
    [x, y] = detector_coordinates_multiray(options);
    
else
    if options.arc_correction
        [~, ~, xp, yp] = detector_coordinates_multiray(options);
        new_x = zeros(options.Ndist * options.Nang, 2, size(xp,2));
        new_y = zeros(options.Ndist * options.Nang, 2, size(yp,2));
        for Nightwish = 1 : size(yp,2)
            [new_x(:,:,Nightwish), new_y(:,:,Nightwish), options] = arcCorrection(options, xp(:,Nightwish), yp(:,Nightwish), false);
        end
        x = new_x;
        y = new_y;
        clear new_x new_y
    else
        [x, y] = sinogram_coordinates_2D_multiray(options);
    end
    if options.sampling > 1
        new_x = zeros(options.Ndist * options.sampling * options.Nang, size(x,2), size(x,3));
        new_y = zeros(options.Ndist * options.sampling * options.Nang, size(y,2), size(y,3));
        for KalPa = 1 : size(y,3)
            [new_x(:,:,KalPa), new_y(:,:,KalPa), options] = increaseSampling(options, x(:,:,KalPa), y(:,:,KalPa), false);
        end
        x = new_x;
        y = new_y;
        clear new_x new_y
    end
end
% x = x(:);
% y = y(:);
if options.implementation == 2 || options.implementation == 3
    x = single(x);
    y = single(y);
end
end

