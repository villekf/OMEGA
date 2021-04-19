function [x, y, z, options] = get_coordinates(options, varargin)
% GET_COORDINATES Get the coordinates for the current machine
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

if nargin >= 2
    rings = varargin{1};
else
    rings = options.rings;
end
if nargin >= 3
    pseudot = varargin{2};
else
    pseudot = [];
end
if nargin >= 4
    interpolateSinogram = varargin{3};
else
    interpolateSinogram = true;
end

if isfield(options,'x') && isfield(options,'y') && (isfield(options,'z') || isfield(options,'z_det'))
    x = options.x(:);
    y = options.y(:);
    if isfield(options,'z')
        z = options.z(:);
    elseif isfield(options,'z_det')
        z = options.z_det(:);
    else
        z = zeros(numel(x),1);
    end
    %     x = x + max(abs(x(:)));
    %     y = y + max(abs(y(:)));
    %     z = z + max(abs(z(:)));
else
    if options.use_raw_data == false
        [~, ~, xp, yp] = detector_coordinates(options);
        [x, y] = sinogram_coordinates_2D(options, xp, yp);
        
        if options.arc_correction && ~options.precompute_lor
            [~, ~, xp, yp] = detector_coordinates(options);
            [x, y, options] = arcCorrection(options, xp, yp, interpolateSinogram);
        end
        if options.sampling > 1 && ~options.precompute_lor
            [x, y, options] = increaseSampling(options, x, y, interpolateSinogram);
        end
        z = sinogram_coordinates_3D(options);
        if options.NSinos ~= options.TotSinos
            z = z(1:options.NSinos,:);
        end
    else
        if options.det_per_ring < options.det_w_pseudo
            options.offangle = options.offangle / options.det_w_pseudo;
            options.offangle = options.offangle * options.det_per_ring;
        end
        [x, y] = detector_coordinates(options);
        if options.sampling_raw > 1 && ~options.precompute_lor
            [x, y, options] = increaseSampling(options, x, y, interpolateSinogram);
        end
        
        z_length = double(rings + 1 + sum(options.pseudot)) * options.cr_pz;
        z = linspace(0, z_length, rings + 2 + sum(options.pseudot))';
        if sum(pseudot) > 0
            z(pseudot) = [];
        end
    end
    if min(z(:)) == 0
        z = z + (options.axial_fov - (options.rings + sum(options.pseudot)) * options.cr_pz)/2 + options.cr_pz/2;
    end
    
    if options.use_raw_data
        
        %         z = z + options.cr_pz/2;
        z(end) = [];
    end
end

if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
    x = single(x);
    y = single(y);
    z = single(z);
else
    x = double(x);
    y = double(y);
    z = double(z);
end