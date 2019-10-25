function [x, y, z] = get_coordinates(options, rings)
% GET_COORDINATES Get the coordinates for the current machine
%   Outputs the x-, y- and z-coordinates for the current machine, depending
%   on whether sinogram or raw data is used
%
% EXAMPLE:
%   [x, y, z] = get_coordinates(options, rings)
% INPUTS:
%   options = Machine properties, sinogram properties and whether raw data
%   is used are needed.
%   rings = Number of crystal rings
%
% See also sinogram_coordinates_2D, sinogram_coordinates_3D, detector_coordinates

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



if options.use_raw_data == false
    [x, y] = sinogram_coordinates_2D(options);
    z = sinogram_coordinates_3D(options);
    if options.NSinos ~= options.TotSinos
        z = z(1:options.NSinos,:);
    end
else
    [x, y] = detector_coordinates(options);
    
    z_length = double(rings) * options.cr_pz;
    z = linspace(0, z_length, rings + 1)';
end