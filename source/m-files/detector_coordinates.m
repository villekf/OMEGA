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

% offangle for which the pseudo-detector block size warning has been given
% (the function is called several times per run)
persistent offangleWarned
if nargout > 4
    error('Too many output arguments')
end
if ~isfield(options,'listmode')
    options.listmode = false;
end
if ~isfield(options,'nLayers')
    options.nLayers = 1;
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
usePseudo = options.det_w_pseudo > options.det_per_ring;
if usePseudo && options.nLayers < 2
    
    [xp,yp] = computeCoordinates(blocks_per_ring, transaxial_multip,cryst_per_block, diameter, cr_p, true, gapsFor(options, cryst_per_block), fbaFor(options));
    x = xp;
    y = yp;
    
elseif options.nLayers > 1
    % Both layers are mapped into blocks of blockSize detectors as in the loader: the block size including the pseudo
    % detectors when they are used, otherwise the block size of the layer with the most crystals per block. The detectors
    % that a layer does not have are located at the end of each block (Inf)
    nBlocksRing = blocks_per_ring * transaxial_multip;
    if usePseudo
        blockSize = options.det_w_pseudo / nBlocksRing;
        if blockSize ~= floor(blockSize) || blockSize < max(options.cryst_per_block(1:2)) + 1
            error(['det_w_pseudo (' num2str(options.det_w_pseudo) ') has to be a multiple of the number of blocks per ring (' num2str(nBlocksRing) ') with at least one pseudo detector per block!'])
        end
    else
        blockSize = max(options.cryst_per_block(1:2));
    end
    [x1,y1] = computeCoordinates(blocks_per_ring, transaxial_multip,options.cryst_per_block(1), diameter, cr_p, usePseudo, gapsFor(options, options.cryst_per_block(1)), fbaFor(options));
    [x1,y1] = padLayerBlocks(x1, y1, options.cryst_per_block(1) + usePseudo, blockSize, nBlocksRing);
    % orig_diameter = diameter + options.crystH(end) * 2;
    diameter = diameter + DOI * 2 + options.crystH(1) * 2;
    [x2,y2] = computeCoordinates(blocks_per_ring, transaxial_multip,options.cryst_per_block(2), diameter, cr_p, usePseudo, gapsFor(options, options.cryst_per_block(2)), fbaFor(options));
    [x2,y2] = padLayerBlocks(x2, y2, options.cryst_per_block(2) + usePseudo, blockSize, nBlocksRing);
    koko = size(x1,1);
    xp = [x1;x2];
    yp = [y1;y2];
    x = xp;
    x(any(isinf(x),2),:) = [];
    y = yp;
    y(any(isinf(y),2),:) = [];
else
    %% Non-pseudo detectors

    [x,y] = computeCoordinates(blocks_per_ring, transaxial_multip,cryst_per_block, diameter, cr_p, false, gapsFor(options, cryst_per_block), fbaFor(options));
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
    if options.det_w_pseudo > options.det_per_ring
        % The pseudo-inclusive arrays are rotated together with their pseudo
        % gaps. Unless the shift is a whole number of pseudo-inclusive
        % blocks, the gaps no longer match the pseudo slots of the data.
        blockSizePseudo = round(options.det_w_pseudo / blocks_per_ring);
        if blockSizePseudo > 0 && mod(round(options.offangle), blockSizePseudo) ~= 0 && ~isequal(offangleWarned, options.offangle)
            offangleWarned = options.offangle;
            warning('offangle (%g) is not a multiple of the block size including pseudo detectors (%d). The detector coordinates are rotated together with the pseudo-detector gaps, so the gaps no longer align with the measured data and the geometry is distorted. Use a multiple of %d for exact results.', options.offangle, blockSizePseudo, blockSizePseudo)
        end
    end
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

function [x, y] = padLayerBlocks(x, y, slots, blockSize, nBlocks)
% Maps the detector coordinates of one layer into blocks of blockSize detectors (the same mapping as in the list-mode
% loader). The slots of each block that the layer does not have are set to Inf and located at the end of the block.
if slots >= blockSize
    return
end
xOut = Inf(blockSize, nBlocks, class(x));
yOut = Inf(blockSize, nBlocks, class(y));
xOut(1 : slots, :) = reshape(x, slots, nBlocks);
yOut(1 : slots, :) = reshape(y, slots, nBlocks);
x = xOut(:);
y = yOut(:);
end

function gaps = gapsFor(options, nCryst)
% Returns options.transaxialGaps if it applies to a block of nCryst crystals, otherwise [] (no gaps). A layer of a
% multi-layer scanner with a different number of crystals per block gets no gaps.
gaps = [];
if isfield(options, 'transaxialGaps') && ~isempty(options.transaxialGaps) && any(options.transaxialGaps(:))
    if numel(options.transaxialGaps) == nCryst - 1
        gaps = options.transaxialGaps;
    elseif ~isfield(options, 'nLayers') || options.nLayers < 2
        error(['options.transaxialGaps has to contain cryst_per_block - 1 = ' num2str(nCryst - 1) ' elements (one gap after each crystal of a block except the last) or be empty, but it has ' num2str(numel(options.transaxialGaps)) ' elements!'])
    end
end
end

function fba = fbaFor(options)
% Angle of the first block (degrees), 0 if not set
fba = 0;
if isfield(options, 'firstBlockAngle') && ~isempty(options.firstBlockAngle)
    fba = options.firstBlockAngle;
end
end
