function [x_final, y_final, xp_final, yp_final] = detector_coordinates_multiray(options)
%% Compute raw detector coordinates
% This function computes the original detector coordinates for one ring,
% both non-pseudo and pseudo cases (the latter only done if pseudo
% detectors are present).

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

x_final = zeros(blocks_per_ring*cryst_per_block,3);
y_final = zeros(blocks_per_ring*cryst_per_block,3);

hh = 1;

for kk = - 1 : 1 : 1
    
    cr = kk * cr_p / 3;
    
    % starting points
    alkupistex = diameter;
    % alkupistey = diameter/2 - (cryst_per_block/2 + 0.5)*cr_p;
    alkupistey = diameter/2 - (cryst_per_block/2 + 0.5)*cr_p + cr;
    
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
                x(ii) = alkupistex - (cr_p * 0.5)*cosd(angle(blocks)) - erotus * cosd(angle(blocks - 1)) - erotus * cosd(angle(blocks)) - (cr_p * 0.5)*cosd(angle(blocks - 1)) - cr*cosd(angle(blocks));
                y(ii) = alkupistey + (cr_p * 0.5)*sind(angle(blocks)) + erotusy * sind(angle(blocks - 1)) + erotusy * sind(angle(blocks)) + (cr_p * 0.5)*sind(angle(blocks - 1)) + cr*sind(angle(blocks));
            else
                % While in the same block
                x(ii) = alkupistex - (cr_p)*cosd(angle(blocks));
                y(ii) = alkupistey + (cr_p)*sind(angle(blocks));
            end
            alkupistex = x(ii);
            alkupistey = y(ii);
            ii = ii + 1;
        end
        alkupistex = x(ii - 1) + cr*cosd(angle(blocks));
        alkupistey = y(ii - 1) - cr*sind(angle(blocks));
    end
    if mod(blocks_per_ring,4) == 0
        blocks = blocks + 1;
        x(ii) = alkupistex - (cr_p * 0.5)*cosd(angle(blocks)) - erotus * cosd(angle(blocks - 1)) - erotus * cosd(angle(blocks)) - (cr_p * 0.5)*cosd(angle(blocks - 1)) - cr*cosd(angle(blocks));
        y(ii) = alkupistey + (cr_p * 0.5)*sind(angle(blocks)) + erotusy * sind(angle(blocks - 1)) + erotusy * sind(angle(blocks)) + (cr_p * 0.5)*sind(angle(blocks - 1)) + cr*sind(angle(blocks));
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
    end
    x_final(:,hh) = x;
    y_final(:,hh) = y;
    hh = hh + 1;
end
%%
% This section ensures symmetry of the coordinates
% Fills the rest of coordinates (3/4)

% If the blocks can be divided by 4, compute an extra block (i.e. the top
% block) first
split = [3 2 1];
x_apu = x_final;
y_apu = y_final;
for hh = 1 : 3
    if mod(blocks_per_ring,4) == 0
        % Symmmetry
            x_final(ii:ii + (blocks_per_ring*cryst_per_block)/4 - 1,hh) = -(flip(x_apu(1:(blocks_per_ring*cryst_per_block)/4,split(hh)) - diameter));
            y_final(ii:ii + (blocks_per_ring*cryst_per_block)/4 - 1,hh) = (flip(y_apu(1:(blocks_per_ring*cryst_per_block)/4,split(hh))));
%             x_final(ii + (blocks_per_ring*cryst_per_block)/4:end) = flip(x_final(cryst_per_block + 1:(blocks_per_ring*cryst_per_block)/2,hh));
%             y_final(ii + (blocks_per_ring*cryst_per_block)/4:end) = -(y_final(cryst_per_block + 1:(blocks_per_ring*cryst_per_block)/2,hh) - diameter);
    else
        x_final(ii:ii + (blocks_per_ring*cryst_per_block)/4 + cryst_per_block/2-1,hh) = -(flip(x_apu(1:(blocks_per_ring*cryst_per_block)/4 + cryst_per_block/2,split(hh))) - diameter);
        y_final(ii:ii + (blocks_per_ring*cryst_per_block)/4 + cryst_per_block/2-1,hh) = (flip(y_apu(1:(blocks_per_ring*cryst_per_block)/4 + cryst_per_block/2,split(hh))));
    end
end
x_apu = x_final;
y_apu = y_final;
for hh = 1 : 3
    if mod(blocks_per_ring,4) == 0
        x_final(ii + (blocks_per_ring*cryst_per_block)/4:end,hh) = flip(x_apu(cryst_per_block + 1:(blocks_per_ring*cryst_per_block)/2,split(hh)));
        y_final(ii + (blocks_per_ring*cryst_per_block)/4:end,hh) = -(y_apu(cryst_per_block + 1:(blocks_per_ring*cryst_per_block)/2,hh) - diameter);
    else
        x_final(ii + (blocks_per_ring*cryst_per_block)/4 + cryst_per_block/2:end,(hh)) = flip(x_apu(cryst_per_block + 1:(blocks_per_ring*cryst_per_block)/2,split(hh)));
        y_final(ii + (blocks_per_ring*cryst_per_block)/4 + cryst_per_block/2:end,(hh)) = -(y_apu(cryst_per_block + 1:(blocks_per_ring*cryst_per_block)/2,(hh)) - diameter);
    end
end

%% Pseudo detectors

% Same operations as above, but with pseudo detectors (usually +1 detector
% per block per ring)

% determine if pseudo detectors are present
if ~isempty(options.pseudot) && options.pseudot > 0
    
    angle = linspace(90,-270,blocks_per_ring + 1)';
    % starting points
    alkupistex = diameter;
    alkupistey = diameter/2 - ((cryst_per_block + 1)/2 + 0.5)*cr_p;
    
    % Width of the topmost blocks (diameter)
    widths = (cryst_per_block + 1)*cr_p*cosd(angle(1:(blocks_per_ring/2 + 1)));
    
    erotus = (diameter - sum(abs(widths)))/sum(cosd(angle(1:(blocks_per_ring/2 + 1))))/2;
    xp_final = zeros(blocks_per_ring*cryst_per_block,3);
    yp_final = zeros(blocks_per_ring*cryst_per_block,3);
    
    hh = 1;
    
    for kk = - 1 : 1 : 1
        
        cr = kk * cr_p / 3;
        
        ii = 1;
        xp = zeros(blocks_per_ring*(cryst_per_block + 1),1);
        yp = zeros(blocks_per_ring*(cryst_per_block + 1),1);
        
        % Compute the detector coordinates of each detector (crystal) in each block
        for blocks = 1:ceil(blocks_per_ring/4)
            for crystals = 1:cryst_per_block+1
                if blocks > 1 && crystals == 1
                    x(ii) = alkupistex - (cr_p * 0.5)*cosd(angle(blocks)) - erotus * cosd(angle(blocks - 1)) - erotus * cosd(angle(blocks)) - (cr_p * 0.5)*cosd(angle(blocks - 1)) - cr*cosd(angle(blocks));
                    y(ii) = alkupistey + (cr_p * 0.5)*sind(angle(blocks)) + erotusy * sind(angle(blocks - 1)) + erotusy * sind(angle(blocks)) + (cr_p * 0.5)*sind(angle(blocks - 1)) + cr*sind(angle(blocks));
                else
                    xp(ii) = alkupistex - (cr_p)*cosd(angle(blocks));
                    yp(ii) = alkupistey + (cr_p)*sind(angle(blocks));
                end
                alkupistex = xp(ii);
                alkupistey = yp(ii);
                ii = ii + 1;
            end
            alkupistex = x(ii - 1) + cr*cosd(angle(blocks));
            alkupistey = y(ii - 1) - cr*sind(angle(blocks));
        end
        if mod(blocks_per_ring,4) == 0
            blocks = blocks + 1;
            x(ii) = alkupistex - (cr_p * 0.5)*cosd(angle(blocks)) - erotus * cosd(angle(blocks - 1)) - erotus * cosd(angle(blocks)) - (cr_p * 0.5)*cosd(angle(blocks - 1)) - cr*cosd(angle(blocks));
            y(ii) = alkupistey + (cr_p * 0.5)*sind(angle(blocks)) + erotusy * sind(angle(blocks - 1)) + erotusy * sind(angle(blocks)) + (cr_p * 0.5)*sind(angle(blocks - 1)) + cr*sind(angle(blocks));
            alkupistex = xp(ii);
            alkupistey = yp(ii);
            ii = ii+1;
            for ll = 2:cryst_per_block+1
                xp(ii) = alkupistex - (cr_p)*cosd(angle(blocks));
                yp(ii) = alkupistey + (cr_p)*sind(angle(blocks));
                alkupistex = xp(ii);
                alkupistey = yp(ii);
                ii = ii + 1;
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
        xp_final(:,hh) = xp;
        yp_final(:,hh) = yp;
        hh = hh + 1;
    end
else
    % If pseudo detectors are not present
    xp_final = x_final;
    yp_final = y_final;
end

% flip the image if so desired
if options.flip_image
    x_final = flip(x_final);
    xp_final = flip(xp_final);
    y_final = flip(y_final);
    yp_final = flip(yp_final);
end
% Rotate the image if desired
if options.offangle ~= 0
    x_final = circshift(x_final, options.offangle);
    xp_final = circshift(xp_final, options.offangle);
    y_final = circshift(y_final, options.offangle);
    yp_final = circshift(yp_final, options.offangle);
end

% save([options.machine_name '_detector_coordinates.mat'], 'x', 'y', 'angle', 'xp', 'yp')

end