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
% Copyright (C) 2022 Ville-Veikko Wettenhovi
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
%     if min(z) < 0
%         z = z - min(z);
%     end
    %     x = x + max(abs(x(:)));
    %     y = y + max(abs(y(:)));
    %     z = z + max(abs(z(:)));
else
    if options.use_raw_data == false
        [~, ~, xp, yp] = detector_coordinates(options);
        if options.nLayers > 1
            koko = numel(xp)/2;
            x = zeros(options.Ndist * options.Nang * 4, 2);
            y = zeros(options.Ndist * options.Nang * 4, 2);
            for kk = 1 : options.nLayers
                if options.cryst_per_block(1) == options.cryst_per_block(2)
                    [x1, y1] = sinogram_coordinates_2D(options, xp(1 + (kk - 1) * koko : kk * koko), yp(1 + (kk - 1) * koko : kk * koko));
                else
                    if kk == 2
                        [x1, y1] = sinogram_coordinates_2D(options, xp(1 + (kk - 1) * koko : kk * koko), yp(1 + (kk - 1) * koko : kk * koko), options.nLayers);
                    else
                        [x1, y1] = sinogram_coordinates_2D(options, xp(1 + (kk - 1) * koko : kk * koko), yp(1 + (kk - 1) * koko : kk * koko));
                    end
                end
                x1(ismember(x1, [0 0], 'rows'),:) = repmat([inf inf], nnz(ismember(x1, [0 0], 'rows')), 1);
                y1(ismember(y1, [0 0], 'rows'),:) = repmat([inf inf], nnz(ismember(y1, [0 0], 'rows')), 1);
                x(1 + (kk - 1) * options.Ndist * options.Nang * 3: options.Ndist * options.Nang + (kk - 1) * options.Ndist * options.Nang * 3,:) = x1;
                y(1 + (kk - 1) * options.Ndist * options.Nang * 3: options.Ndist * options.Nang + (kk - 1) * options.Ndist * options.Nang * 3,:) = y1;
            end
            ind2 = 1 + options.Ndist * options.Nang * 3;
            ind1 = options.Ndist * options.Nang;
            x(1 + options.Ndist * options.Nang : options.Ndist * options.Nang * 2,:) = [x(ind2:end,1) x(1 : ind1, 2)];
            y(1 + options.Ndist * options.Nang : options.Ndist * options.Nang * 2,:) = [y(ind2:end,1) y(1 : ind1, 2)];
            x(1 + options.Ndist * options.Nang * 2 : options.Ndist * options.Nang * 3,:) = [x(1 : ind1, 1) x(ind2:end,2)];
            y(1 + options.Ndist * options.Nang * 2 : options.Ndist * options.Nang * 3,:) = [y(1 : ind1, 1) y(ind2:end,2)];
        else
            [x, y] = sinogram_coordinates_2D(options, xp, yp);
        end
        
        if options.arc_correction && ~options.precompute_lor
%             [~, ~, xp, yp] = detector_coordinates(options);
            [x, y, options] = arcCorrection(options, xp, yp, interpolateSinogram);
        end
        if options.sampling > 1 && ~options.precompute_lor
            [x, y, options] = increaseSampling(options, x, y, interpolateSinogram);
        end
        z = sinogram_coordinates_3D(options);
        if options.NSinos ~= options.TotSinos
            z = z(1:options.NSinos,:);
        end
        x = [x(:,1)';y(:,1)';x(:,2)';y(:,2)'];
        z = z';
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
        z = linspace(-(z_length / 2 - options.cr_pz/2), z_length / 2 - options.cr_pz/2, rings + 1 + sum(options.pseudot))';
        if sum(pseudot) > 0
            z(pseudot) = [];
        end
        x = [x';y'];
    end
%     if min(z(:)) == 0
%         z = z + (options.axial_fov - (options.rings + sum(options.pseudot)) * options.cr_pz)/2 + options.cr_pz/2;
%     end
    
%     if options.use_raw_data
%         
%         %         z = z + options.cr_pz/2;
%         z(end) = [];
%     end
%     z = z - z(end) / 2;
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