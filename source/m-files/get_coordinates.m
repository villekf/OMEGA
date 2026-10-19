function [x, y, z, options] = get_coordinates(options, varargin)
% GET_COORDINATES Get the coordinates for the current scanner
%   Outputs the x-, y- and z-coordinates for the current scanner, depending
%   on whether sinogram or raw data is used
%
% EXAMPLES:
%   [x, y, z, options] = get_coordinates(options)
%   [x, y, z, options] = get_coordinates(options, rings, pseudot)
% INPUTS:
%   options = Scanner properties, sinogram properties and whether raw data
%   is used are needed.
%   rings = Number of crystal rings (required only for raw data)
%   pseudot = The numbers of the pseudo rings (required only for raw data)
%
% See also sinogram_coordinates_2D, sinogram_coordinates_3D,
% detector_coordinates, getMultirayCoordinates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019-2024 Ville-Veikko Wettenhovi
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
    rings = varargin{1};
else
    rings = options.rings;
end
if nargin >= 3 && ~isempty(varargin{2})
    pseudot = varargin{2};
else
    pseudot = [];
end
if nargin >= 4 && ~isempty(varargin{3})
    interpolateSinogram = varargin{3};
else
    interpolateSinogram = false;
end

if (isfield(options,'x') && isfield(options,'y') && (isfield(options,'z') || isfield(options,'z_det'))) && ~options.listmode
    x = options.x(:);
    y = options.y(:);
    if isfield(options,'z')
        z = options.z(:);
    elseif isfield(options,'z_det')
        z = options.z_det(:);
    else
        z = zeros(numel(x),1);
    end
else
    if options.use_raw_data == false
        [~, ~, xp, yp] = detector_coordinates(options);
        if options.nLayers > 1
            for kk = 1 : options.nLayers^2
                if kk == 1
                    [x1, y1] = sinogram_coordinates_2D(options, xp, yp, [1, 1]);
                    z1 = sinogram_coordinates_3D(options, [1, 1]);
                elseif kk == 2
                    [x2, y2] = sinogram_coordinates_2D(options, xp, yp, [1, 2]);
                    z2 = sinogram_coordinates_3D(options, [1, 2]);
                elseif kk == 3
                    [x3, y3] = sinogram_coordinates_2D(options, xp, yp, [2, 1]);
                    z3 = sinogram_coordinates_3D(options, [2, 1]);
                elseif kk == 4
                    [x4, y4] = sinogram_coordinates_2D(options, xp, yp, [2, 2]);
                    z4 = sinogram_coordinates_3D(options, [2, 2]);
                end
            end
            x = [x1;x2;x3;x4];
            y = [y1;y2;y3;y4];
            z = [z1;z2;z3;z4];
            % The layer combination (block index) of each sinogram is stored as the third value of z so that it follows
            % the sinogram when the subsets are formed
            z = [z, [zeros(size(z1,1),1); ones(size(z2,1),1); 2 * ones(size(z3,1),1); 3 * ones(size(z4,1),1)]];
        else
            [x, y] = sinogram_coordinates_2D(options, xp, yp);
            z = sinogram_coordinates_3D(options);
            if options.NSinos ~= options.TotSinos
                z = z(1:options.NSinos,:);
            end
        end

        if options.arc_correction && ~options.precompute_lor
            [x, y, options] = arcCorrection(options, interpolateSinogram);
        end
        if options.sampling > 1 && ~options.precompute_lor
            [x, y, options] = increaseSampling(options, x, y, interpolateSinogram);
        end
        x = [x(:,1)';y(:,1)';x(:,2)';y(:,2)'];
        z = z';
    else
        if options.det_per_ring < options.det_w_pseudo
            options.offangle = options.offangle / options.det_w_pseudo;
            options.offangle = options.offangle * options.det_per_ring;
        end
        if options.nLayers > 1
            % Multi-layer: one block of det_w_pseudo detectors per layer (the loader's layerID * det_w_pseudo offset)
            [~, ~, xp, yp] = detector_coordinates(options);
            [x, y] = rawLayerDetectors(options, xp, yp);
        else
            [x, y] = detector_coordinates(options);
        end
        if options.sampling_raw > 1 && ~options.precompute_lor
            [x, y, options] = increaseSampling(options, x, y, interpolateSinogram);
        end

        z_length = double(rings + sum(options.pseudot)) * options.cr_pz;
        z = linspace(-(z_length / 2 - options.cr_pz/2), z_length / 2 - options.cr_pz/2, rings + sum(options.pseudot))';
        if sum(pseudot) > 0
            z(pseudot) = [];
        end
        % Gaps between the axial blocks
        if isfield(options, 'ringGaps') && numel(options.ringGaps) > 0 && sum(options.ringGaps) > 0
            gaps = repelem([0;cumsum(options.ringGaps(:))], options.cryst_per_block_axial(1), 1);
            z = z + gaps(1:numel(z)) - sum(options.ringGaps) / 2;
        end
        % Multi-layer: the rings of layer l are located at an offset of l * rings
        if options.nLayers > 1
            z = rawLayerRings(options, z, rings);
        end
        x = [x';y'];
    end
end
% cType is normally set by projectorClass, default to single when get_coordinates is called without it
if ~isfield(options, 'cType')
    options.cType = 'single';
end
x = cast(x, options.cType);
y = cast(y, options.cType);
z = cast(z, options.cType);

function [x, y] = rawLayerDetectors(options, xp, yp)
% Raw/index-based detector coordinates of a multi-layer scanner. The coordinates of layer l are stored in a block
% of det_w_pseudo entries starting at l * det_w_pseudo, i.e. the same offset the list-mode loader adds to the
% transaxial detector index of a layer. A block contains the real detectors of the layer in the order of the
% loader's detector index, any padding (Inf) is moved to the end of the block
S = options.det_w_pseudo(1);
nLayers = options.nLayers;
if numel(xp) < nLayers * S
    error('Raw/index-based multi-layer detector coordinates are not supported with pseudo detectors!')
end
x = Inf(nLayers * S, 1, class(xp));
y = Inf(nLayers * S, 1, class(yp));
for ll = 1 : nLayers
    if ll < nLayers
        ind = (ll - 1) * S + 1 : ll * S;
    else
        ind = (ll - 1) * S + 1 : numel(xp);
    end
    valid = isfinite(xp(ind));
    xl = xp(ind);
    yl = yp(ind);
    n = nnz(valid);
    if n > S
        error('Too many detectors in layer %d!', ll)
    end
    x((ll - 1) * S + (1 : n)) = xl(valid);
    y((ll - 1) * S + (1 : n)) = yl(valid);
end

function zL = rawLayerRings(options, z, rings)
% Raw/index-based ring coordinates of a multi-layer scanner. The coordinates of layer l are stored in a block of
% rings(1) entries starting at l * rings(1), i.e. the same offset the list-mode loader adds to the ring number of
% a layer. With layers that have a different number of axial crystals (cryst_per_block_axial) the rings of each
% layer are obtained from sinogram_coordinates_3D (which includes the gaps between the modules) and the real
% rings are followed by Inf padding up to rings(1). Otherwise all layers have the same axial positions
nLayers = options.nLayers;
if numel(options.cryst_per_block_axial) > 1 && options.cryst_per_block_axial(1) ~= options.cryst_per_block_axial(2)
    S = rings(1);
    % The first loppu rows of the second column of the span 1 coordinates are the axial positions of the rings of
    % the second LOR end (the layer given by the second element of layer)
    opts = options;
    opts.span = 1;
    zL = Inf(nLayers * S, 1);
    for ll = 1 : nLayers
        zl = sinogram_coordinates_3D(opts, [ll, ll]);
        zl = zl(1 : S, 2);
        zl = zl(isfinite(zl));
        if numel(zl) > S
            error('Too many rings in layer %d!', ll)
        end
        zL((ll - 1) * S + (1 : numel(zl))) = zl;
    end
else
    zL = repmat(z, nLayers, 1);
end