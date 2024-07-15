function [varargout] = detector_coordinates(options)
%% Compute raw detector coordinates
% This function computes the original detector coordinates for one ring,
% both non-pseudo and pseudo cases (the latter only done if pseudo
% detectors are present).
%
% OUTPUTS:
%   x = X detector coordinates
%   y = Y detector coordinates
%   xp = X detector coordinates with pseudo detector(s)
%   yp = Y detector coordinates with pseudo detector(s)
%
% See also sinogram_coordinates_2D, computeCoordinates

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

if nargout > 4
    error('Too many output arguments')
end
if ~isfield(options,'listmode')
    options.listmode = false;
end

cr_p = options.cr_p;
diameter = options.diameter;
cryst_per_block = options.cryst_per_block(1);
blocks_per_ring = options.blocks_per_ring;
x = 0;
y = 0;

if isfield(options, 'DOI')
    DOI = options.DOI;
else
    DOI = 0;
end

if isfield(options, 'transaxial_multip')
    transaxial_multip = options.transaxial_multip;
else
    transaxial_multip = 1;
end

% orig_diameter = diameter;
diameter = diameter + DOI * 2;


%% Pseudo detectors

% Same operations as above, but with pseudo detectors (usually +1 detector
% per block per ring)

% determine if pseudo detectors are present
if options.det_w_pseudo > options.det_per_ring
    
    [xp,yp] = computeCoordinates(blocks_per_ring, transaxial_multip,cryst_per_block, diameter, cr_p, true);
    x = xp;
    y = yp;
    
elseif options.nLayers > 1
    if options.cryst_per_block(1) == options.cryst_per_block(2)
        [x1,y1] = computeCoordinates(blocks_per_ring, transaxial_multip,cryst_per_block, diameter, cr_p, false);
        % orig_diameter = diameter + options.crystH(end) * 2;
        diameter = diameter + DOI * 2 + options.crystH(1) * 2;
        [x2,y2] = computeCoordinates(blocks_per_ring, transaxial_multip,cryst_per_block, diameter, cr_p, false);
    elseif options.cryst_per_block(1) > options.cryst_per_block(2)
        [x1,y1] = computeCoordinates(blocks_per_ring, transaxial_multip,options.cryst_per_block(1), diameter, cr_p, false);
        % orig_diameter = diameter + options.crystH(end) * 2;
        diameter = diameter + DOI * 2 + options.crystH(1) * 2;
        [x2,y2] = computeCoordinates(blocks_per_ring, transaxial_multip,options.cryst_per_block(end), diameter, cr_p, false);
        x = Inf(options.det_w_pseudo * options.nLayers,1);
        insert_indices = setdiff(1:length(x)/2, options.cryst_per_block(1):options.cryst_per_block(1):length(x)/2);
        x(insert_indices) = x2;
        y = Inf(options.det_w_pseudo * options.nLayers,1);
        y(insert_indices) = y2;
        x2 = x;
        y2 = y;
        koko = size(x2,1);
    else
        [x1,y1] = computeCoordinates(blocks_per_ring, transaxial_multip,options.cryst_per_block(1), diameter, cr_p, false);
        x = Inf(options.det_w_pseudo / 2 * options.nLayers,1);
        insert_indices = setdiff(1:length(x), options.cryst_per_block(end):options.cryst_per_block(end):length(x));
        x(insert_indices) = x1;
        y = Inf(options.det_w_pseudo / 2 * options.nLayers,1);
        y(insert_indices) = y1;
        x1 = x;
        y1 = y;
        % orig_diameter = diameter + options.crystH(end) * 2;
        diameter = diameter + DOI * 2 + options.crystH(1) * 2;
        [x2,y2] = computeCoordinates(blocks_per_ring, transaxial_multip,options.cryst_per_block(end), diameter, cr_p, false);
        koko = size(x1,1);
    end
    xp = [x1;x2];
    yp = [y1;y2];
    x = xp;
    x(any(isinf(x),2),:) = [];
    y = yp;
    y(any(isinf(y),2),:) = [];
else
    %% Non-pseudo detectors

    [x,y] = computeCoordinates(blocks_per_ring, transaxial_multip,cryst_per_block, diameter, cr_p, false);
    % If pseudo detectors are not present
    xp = x;
    yp = y;
end

% flip the image if so desired
if options.flip_image
    x = flip(x);
    xp = flip(xp);
    y = flip(y);
    yp = flip(yp);
end
% Rotate the image if desired
if options.offangle ~= 0
    if options.nLayers > 1
        x1 = xp(1 : koko,:);
        x1 = circshift(x1, round(options.offangle));
        x2 = xp(koko + 1 : end,:);
        x2 = circshift(x2, round(options.offangle));
        xp = [x1;x2];
        y1 = yp(1 : koko,:);
        y1 = circshift(y1, round(options.offangle));
        y2 = yp(koko + 1 : end,:);
        y2 = circshift(y2, round(options.offangle));
        yp = [y1;y2];
    else
        x = circshift(x, round(options.offangle));
        xp = circshift(xp, round(options.offangle));
        y = circshift(y, round(options.offangle));
        yp = circshift(yp, round(options.offangle));
    end
end

if (options.implementation == 2 || options.implementation == 3 || options.implementation == 5) && options.listmode > 0
    x = single(x);
    y = single(y);
    xp = single(xp);
    yp = single(yp);
end

if nargout >= 1
    varargout{1} = x;
end
if nargout >= 2
    varargout{2} = y;
end
if nargout >= 3
    varargout{3} = xp;
end
if nargout == 4
    varargout{4} = yp;
end

end