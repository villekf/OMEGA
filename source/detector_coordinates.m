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
% See also sinogram_coordinates_2D

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

if nargout > 4
    error('Too many output arguments')
end

cr_p = options.cr_p;
diameter = options.diameter;
cryst_per_block = options.cryst_per_block;
blocks_per_ring = options.blocks_per_ring;

%% Non-pseudo detectors

% All angles, starting from 90 degrees
angle = 90:-360/blocks_per_ring:-270;

% Width of the topmost blocks (diameter)
widths = (cryst_per_block)*cr_p*cosd(angle(1:(blocks_per_ring/2 + 1)));
widthsy = (cryst_per_block)*cr_p*sind(angle(1:(blocks_per_ring/2)));

% Gap between adjacent blocks
% If negative, then the crystal size is smaller on the edges of the block
% than the crystal pitch
erotus = (diameter - sum(abs(widths)))/sum(cosd(angle(1:(blocks_per_ring/2 + 1))))/2;
erotusy = (diameter - sum(abs(widthsy)))/sum(abs(sind(angle(1:(blocks_per_ring/2)))))/2;

% starting points
alkupistex = diameter;
alkupistey = diameter/2 - (cryst_per_block/2 + 0.5)*cr_p;

ii = 1;
x = zeros(blocks_per_ring*cryst_per_block,1);
y = zeros(blocks_per_ring*cryst_per_block,1);
%%
% Compute the detector coordinates of each detector (crystal) in each block
% Only for the 1/4th of the ring
for blocks = 1:ceil(blocks_per_ring/4)
    for crystals = 1:cryst_per_block
        % Moving to the next block
        if blocks > 1 && crystals == 1
            x(ii) = alkupistex - (cr_p * 0.5)*cosd(angle(blocks)) - erotus * cosd(angle(blocks - 1)) - erotus * cosd(angle(blocks)) - (cr_p * 0.5)*cosd(angle(blocks - 1));
            y(ii) = alkupistey + (cr_p * 0.5)*sind(angle(blocks)) + erotusy * sind(angle(blocks - 1)) + erotusy * sind(angle(blocks)) + (cr_p * 0.5)*sind(angle(blocks - 1));
        else
            % While in the same block
            x(ii) = alkupistex - (cr_p)*cosd(angle(blocks));
            y(ii) = alkupistey + (cr_p)*sind(angle(blocks));
        end
        alkupistex = x(ii);
        alkupistey = y(ii);
        ii=ii+1;
    end
end
%%
% This section ensures symmetry of the coordinates
% Fills the rest of coordinates (3/4)

% If the blocks can be divided by 4, compute an extra block (i.e. the top
% block) first
if mod(blocks_per_ring,4) == 0
    blocks = blocks + 1;
    x(ii) = alkupistex - (cr_p * 0.5)*cosd(angle(blocks)) - erotus * cosd(angle(blocks - 1)) - erotus * cosd(angle(blocks)) - (cr_p * 0.5)*cosd(angle(blocks - 1));
    y(ii) = alkupistey + (cr_p * 0.5)*sind(angle(blocks)) + erotus * sind(angle(blocks - 1)) + erotus * sind(angle(blocks)) + (cr_p * 0.5)*sind(angle(blocks - 1));
    alkupistex = x(ii);
    alkupistey = y(ii);
    ii = ii+1;
    for ll = 2:cryst_per_block
        x(ii) = alkupistex - (cr_p)*cosd(angle(blocks));
        y(ii) = alkupistey + (cr_p)*sind(angle(blocks));
        alkupistex = x(ii);
        alkupistey = y(ii);
        ii=ii+1;
    end
    % Symmmetry
    x(ii:ii + (blocks_per_ring*cryst_per_block)/4 - 1) = -(flip(x(1:(blocks_per_ring*cryst_per_block)/4) - diameter));
    y(ii:ii + (blocks_per_ring*cryst_per_block)/4 - 1) = (flip(y(1:(blocks_per_ring*cryst_per_block)/4)));
    x(ii + (blocks_per_ring*cryst_per_block)/4:end) = flip(x(cryst_per_block + 1:(blocks_per_ring*cryst_per_block)/2));
    y(ii + (blocks_per_ring*cryst_per_block)/4:end) = -(y(cryst_per_block + 1:(blocks_per_ring*cryst_per_block)/2) - diameter);
else
    x(ii:ii + (blocks_per_ring*cryst_per_block)/4 + cryst_per_block/2-1) = -(flip(x(1:(blocks_per_ring*cryst_per_block)/4 + cryst_per_block/2)) - diameter);
    y(ii:ii + (blocks_per_ring*cryst_per_block)/4 + cryst_per_block/2-1) = (flip(y(1:(blocks_per_ring*cryst_per_block)/4 + cryst_per_block/2)));
    x(ii + (blocks_per_ring*cryst_per_block)/4 + cryst_per_block/2:end) = flip(x(cryst_per_block + 1:(blocks_per_ring*cryst_per_block)/2));
    y(ii + (blocks_per_ring*cryst_per_block)/4 + cryst_per_block/2:end) = -(y(cryst_per_block + 1:(blocks_per_ring*cryst_per_block)/2) - diameter);
end

%% Pseudo detectors

% Same operations as above, but with pseudo detectors (usually +1 detector
% per block per ring)

% determine if pseudo detectors are present
if options.det_w_pseudo > options.det_per_ring
    
    angle = linspace(90,-270,blocks_per_ring + 1)';
    % starting points
    alkupistex = diameter;
    alkupistey = diameter/2 - ((cryst_per_block + 1)/2 + 0.5)*cr_p;
    
    % Width of the topmost blocks (diameter)
    widths = (cryst_per_block + 1)*cr_p*cosd(angle(1:(blocks_per_ring/2 + 1)));
    
    % Gap
    erotus = (diameter - sum(abs(widths)))/sum(cosd(angle(1:(blocks_per_ring/2 + 1))))/2;
    
    ii = 1;
    xp = zeros(blocks_per_ring*(cryst_per_block + 1),1);
    yp = zeros(blocks_per_ring*(cryst_per_block + 1),1);
    
    % Compute the detector coordinates of each detector (crystal) in each block
    for blocks = 1:ceil(blocks_per_ring/4)
        for crystals = 1:cryst_per_block+1
            if blocks > 1 && crystals == 1
                xp(ii) = alkupistex - (cr_p * 0.5)*cosd(angle(blocks)) - erotus * cosd(angle(blocks - 1)) - erotus * cosd(angle(blocks)) - (cr_p * 0.5)*cosd(angle(blocks - 1));
                yp(ii) = alkupistey + (cr_p * 0.5)*sind(angle(blocks)) + erotus * sind(angle(blocks - 1)) + erotus * sind(angle(blocks)) + (cr_p * 0.5)*sind(angle(blocks - 1));
            else
                xp(ii) = alkupistex - (cr_p)*cosd(angle(blocks));
                yp(ii) = alkupistey + (cr_p)*sind(angle(blocks));
            end
            alkupistex = xp(ii);
            alkupistey = yp(ii);
            ii=ii+1;
        end
    end
    if mod(blocks_per_ring,4) == 0
        blocks = blocks + 1;
        xp(ii) = alkupistex - (cr_p * 0.5)*cosd(angle(blocks)) - erotus * cosd(angle(blocks - 1)) - erotus * cosd(angle(blocks)) - (cr_p * 0.5)*cosd(angle(blocks - 1));
        yp(ii) = alkupistey + (cr_p * 0.5)*sind(angle(blocks)) + erotus * sind(angle(blocks - 1)) + erotus * sind(angle(blocks)) + (cr_p * 0.5)*sind(angle(blocks - 1));
        alkupistex = xp(ii);
        alkupistey = yp(ii);
        ii = ii+1;
        for ll = 2:cryst_per_block+1
            xp(ii) = alkupistex - (cr_p)*cosd(angle(blocks));
            yp(ii) = alkupistey + (cr_p)*sind(angle(blocks));
            alkupistex = xp(ii);
            alkupistey = yp(ii);
            ii=ii+1;
        end
        xp(ii:ii+(blocks_per_ring*(cryst_per_block+1))/4 - 1) = -(flip(xp(1:(blocks_per_ring*(cryst_per_block+1))/4) - diameter));
        yp(ii:ii+(blocks_per_ring*(cryst_per_block+1))/4 - 1) = (flip(yp(1:(blocks_per_ring*(cryst_per_block+1))/4)));
        xp(ii+(blocks_per_ring*(cryst_per_block+1))/4:end) = flip(xp((cryst_per_block+1)+1:(blocks_per_ring*(cryst_per_block+1))/2));
        yp(ii+(blocks_per_ring*(cryst_per_block+1))/4:end) = -(yp((cryst_per_block+1)+1:(blocks_per_ring*(cryst_per_block+1))/2) - diameter);
    else
        xp(ii:ii+(blocks_per_ring*(cryst_per_block+1))/4+(cryst_per_block+1)/2-1) = -(flip(xp(1:(blocks_per_ring*(cryst_per_block+1))/4 +(cryst_per_block+1)/2)) - diameter);
        yp(ii:ii+(blocks_per_ring*(cryst_per_block+1))/4+(cryst_per_block+1)/2-1) = (flip(yp(1:(blocks_per_ring*(cryst_per_block+1))/4 +(cryst_per_block+1)/2)));
        xp(ii+(blocks_per_ring*(cryst_per_block+1))/4+(cryst_per_block+1)/2:end) = flip(xp((cryst_per_block+1)+1:(blocks_per_ring*(cryst_per_block+1))/2));
        yp(ii+(blocks_per_ring*(cryst_per_block+1))/4+(cryst_per_block+1)/2:end) = -(yp((cryst_per_block+1)+1:(blocks_per_ring*(cryst_per_block+1))/2) - diameter);
    end
    
    
else
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
    x = circshift(x, options.offangle);
    xp = circshift(xp, options.offangle);
    y = circshift(y, options.offangle);
    yp = circshift(yp, options.offangle);
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