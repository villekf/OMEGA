function [x_final, y_final, xp_final, yp_final] = detector_coordinates_multiray(options)
%% Compute raw multiray detector coordinates
% This function computes the original detector coordinates for one ring,
% both non-pseudo and pseudo cases (the latter only done if pseudo
% detectors are present). This function is for the multiray case with
% n_transaxial rays.

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

cr_p = options.cr_p;
diameter = options.diameter;
cryst_per_block = options.cryst_per_block;
blocks_per_ring = options.blocks_per_ring;

if isfield(options, 'DOI')
    DOI = options.DOI;
else
    DOI = 0;
end

diameter = diameter + DOI * 2;

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

if erotus < 0
    cr_p = diameter/sum((cryst_per_block)*cosd(angle(1:(blocks_per_ring/2 + 1))));
    widths = (cryst_per_block)*cr_p*cosd(angle(1:(blocks_per_ring/2 + 1)));
    widthsy = (cryst_per_block)*cr_p*sind(angle(1:(blocks_per_ring/2)));
    erotus = (diameter - sum(abs(widths)))/sum(cosd(angle(1:(blocks_per_ring/2 + 1))))/2;
    erotusy = (diameter - sum(abs(widthsy)))/sum(abs(sind(angle(1:(blocks_per_ring/2)))))/2;
end



hh = 1;

% Is the ray count even or odd
if mod(options.n_rays_transaxial,2) == 1
    kk = - 1 * floor(options.n_rays_transaxial/2) : 1 : 1 * floor(options.n_rays_transaxial/2);
else
    kk = - 1 * options.n_rays_transaxial/2 : 1 : 1 * options.n_rays_transaxial/2;
    kk(kk == 0) = [];
end
cr_kerroin = cr_p / (options.n_rays_transaxial);
cr = cr_kerroin/2 : cr_kerroin : cr_p - cr_kerroin/2;
x_final = zeros(blocks_per_ring*cryst_per_block,length(kk));
y_final = zeros(blocks_per_ring*cryst_per_block,length(kk));
for oo = kk
    
    % starting points
    alkupistex = 0;
    alkupistey = diameter/2 - (cryst_per_block/2 + 1)*cr_p + cr(hh);
    
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
                x(ii) = alkupistex + erotus * cosd(angle(blocks - 1)) + erotus * cosd(angle(blocks)) + cr(hh)*cosd(angle(blocks));
                y(ii) = alkupistey + erotusy * sind(angle(blocks - 1)) + erotusy * sind(angle(blocks)) + cr(hh)*sind(angle(blocks));
            else
                % While in the same block
                x(ii) = alkupistex + (cr_p)*cosd(angle(blocks));
                y(ii) = alkupistey + (cr_p)*sind(angle(blocks));
            end
            alkupistex = x(ii);
            alkupistey = y(ii);
            ii = ii + 1;
        end
        alkupistex = alkupistex + (cr_p - cr(hh))*cosd(angle(blocks));
        alkupistey = alkupistey + (cr_p - cr(hh))*sind(angle(blocks));
    end
    if mod(blocks_per_ring,4) == 0
        blocks = blocks + 1;
        x(ii) = alkupistex + erotus * cosd(angle(blocks - 1)) + erotus * cosd(angle(blocks)) + cr(hh)*cosd(angle(blocks));
        y(ii) = alkupistey + erotusy * sind(angle(blocks - 1)) + erotusy * sind(angle(blocks)) + cr(hh)*sind(angle(blocks));
        alkupistex = x(ii);
        alkupistey = y(ii);
        ii = ii+1;
        for ll = 2:cryst_per_block
            x(ii) = alkupistex + (cr_p)*cosd(angle(blocks));
            y(ii) = alkupistey + (cr_p)*sind(angle(blocks));
            alkupistex = x(ii);
            alkupistey = y(ii);
            ii=ii+1;
        end
    end
    if mod(blocks_per_ring,4) == 0
        y(ii:ii + (blocks_per_ring*cryst_per_block)/4 - 1) = flip(y(1:(blocks_per_ring*cryst_per_block)/4));
    else
        y(ii:ii + (blocks_per_ring*cryst_per_block)/4 - 1 + cryst_per_block/2) = flip(y(1:(blocks_per_ring*cryst_per_block)/4 + cryst_per_block/2));
    end
    
    if mod(blocks_per_ring,4) == 0
        y(options.det_per_ring/2 + cryst_per_block + 1:options.det_per_ring/2 + cryst_per_block + options.det_per_ring/4) = flip(x(1:options.det_per_ring/4));
        y(options.det_per_ring/2 + options.det_per_ring/4 + 1:end) = flip(y(options.det_per_ring/2 + cryst_per_block + 1:options.det_per_ring/2 + cryst_per_block + options.det_per_ring/4));
    else
        y(options.det_per_ring/2 + cryst_per_block + 1 : end) = flipud(abs(options.diameter - y(cryst_per_block + 1 : options.det_per_ring/2)));
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
split = length(kk) : - 1 : 1;
if mod(blocks_per_ring,4) == 0
    for hh = 1 : length(kk)
        x_final(ii:ii + options.det_per_ring/4 - 1,hh) = flip(y_final(options.det_per_ring/4 + 1:options.det_per_ring/2,split(hh)));
        x_final(options.det_per_ring/2 + cryst_per_block + 1:options.det_per_ring/2 + options.det_per_ring/4  + cryst_per_block,hh) = (y_final(options.det_per_ring/4 + cryst_per_block + 1:options.det_per_ring/2 + cryst_per_block,hh));
        x_final(options.det_per_ring/2 + options.det_per_ring/4 + cryst_per_block + 1:end,hh) = y_final(cryst_per_block + options.det_per_ring/2 + 1:options.det_per_ring/2 + options.det_per_ring/4,split(hh));
    end
else
    x_jelppi = x_final;
    for hh = 1 : length(kk)
        x_final(ii:ii + options.det_per_ring/4 - 1 - cryst_per_block/2,hh) = ...
            abs(options.diameter - flip(x_jelppi(cryst_per_block + 1 : ii - 1,(hh))));
        x_final(ii + options.det_per_ring/4 - cryst_per_block/2 : ii + options.det_per_ring/4 + cryst_per_block/2 - 1,hh) = options.diameter;
    end
    x_jelppi = x_final;
    for hh = 1 : length(kk)
        x_final(ii + options.det_per_ring/4 + cryst_per_block/2 : end,hh) = flipud(x_jelppi(cryst_per_block + 1 : ii + options.det_per_ring/4 - 1 - cryst_per_block/2,(hh)));
    end
end

x_final = (circshift(x_final, size(x_final,1)/2));

%% Pseudo detectors

% Same operations as above, but with pseudo detectors (usually +1 detector
% per block per ring)

% determine if pseudo detectors are present
if options.det_w_pseudo > options.det_per_ring
    
    angle = linspace(90,-270,blocks_per_ring + 1)';
    
    
    % Width of the topmost blocks (diameter)
    widths = (cryst_per_block + 1)*cr_p*cosd(angle(1:(blocks_per_ring/2 + 1)));
    widthsy = (cryst_per_block + 1)*cr_p*sind(angle(1:(blocks_per_ring/2)));
    
    erotus = (diameter - sum(abs(widths)))/sum(cosd(angle(1:(blocks_per_ring/2 + 1))))/2;
    erotusy = (diameter - sum(abs(widthsy)))/sum(abs(sind(angle(1:(blocks_per_ring/2)))))/2;
    
    if erotus < 0
        cr_p = diameter/sum((cryst_per_block + 1)*cosd(angle(1:(blocks_per_ring/2 + 1))));
        widths = (cryst_per_block + 1)*cr_p*cosd(angle(1:(blocks_per_ring/2 + 1)));
        widthsy = (cryst_per_block + 1)*cr_p*sind(angle(1:(blocks_per_ring/2)));
        erotus = (diameter - sum(abs(widths)))/sum(cosd(angle(1:(blocks_per_ring/2 + 1))))/2;
        erotusy = (diameter - sum(abs(widthsy)))/sum(abs(sind(angle(1:(blocks_per_ring/2)))))/2;
    end
    xp_final = zeros(blocks_per_ring*cryst_per_block,length(kk));
    yp_final = zeros(blocks_per_ring*cryst_per_block,length(kk));
    
    
    hh = 1;
    
    if mod(options.n_rays_transaxial,2) == 1
        kk = - 1 * floor(options.n_rays_transaxial/2) : 1 : 1 * floor(options.n_rays_transaxial/2);
    else
        kk = - 1 * options.n_rays_transaxial/2 : 1 : 1 * options.n_rays_transaxial/2;
        kk(kk == 0) = [];
    end
    cr_kerroin = cr_p / (options.n_rays_transaxial);
    cr = cr_kerroin/2 : cr_kerroin : cr_p - cr_kerroin/2;
    
    for oo = kk
        
        % starting points
        alkupistex = 0;
        alkupistey = diameter/2 - ((cryst_per_block + 1)/2 + 1)*cr_p + cr(hh);
        
        ii = 1;
        xp = zeros(blocks_per_ring*(cryst_per_block + 1),1);
        yp = zeros(blocks_per_ring*(cryst_per_block + 1),1);
        
        % Compute the detector coordinates of each detector (crystal) in each block
        for blocks = 1:ceil(blocks_per_ring/4)
            for crystals = 1:cryst_per_block + 1
                % Moving to the next block
                if blocks > 1 && crystals == 1
                    xp(ii) = alkupistex + erotus * cosd(angle(blocks - 1)) + erotus * cosd(angle(blocks)) + cr(hh)*cosd(angle(blocks));
                    yp(ii) = alkupistey + erotusy * sind(angle(blocks - 1)) + erotusy * sind(angle(blocks)) + cr(hh)*sind(angle(blocks));
                else
                    % While in the same block
                    xp(ii) = alkupistex + (cr_p)*cosd(angle(blocks));
                    yp(ii) = alkupistey + (cr_p)*sind(angle(blocks));
                end
                alkupistex = xp(ii);
                alkupistey = yp(ii);
                ii = ii + 1;
            end
            alkupistex = alkupistex + (cr_p - cr(hh)) * cosd(angle(blocks));
            alkupistey = alkupistey + (cr_p - cr(hh)) * sind(angle(blocks));
        end
        if mod(blocks_per_ring,4) == 0
            blocks = blocks + 1;
            x(ii) = alkupistex + erotus * cosd(angle(blocks - 1)) + erotus * cosd(angle(blocks)) + cr(hh) * cosd(angle(blocks));
            yp(ii) = alkupistey + erotusy * sind(angle(blocks - 1)) + erotusy * sind(angle(blocks)) + cr(hh) * sind(angle(blocks));
            alkupistex = xp(ii);
            alkupistey = yp(ii);
            ii = ii + 1;
            for ll = 2 : cryst_per_block + 1
                xp(ii) = alkupistex + (cr_p)*cosd(angle(blocks));
                yp(ii) = alkupistey + (cr_p)*sind(angle(blocks));
                alkupistex = xp(ii);
                alkupistey = yp(ii);
                ii=ii+1;
            end
        end
        if mod(blocks_per_ring, 4) == 0
            yp(ii:ii + (blocks_per_ring * (cryst_per_block + 1))/4 - 1) = flip(yp(1:(blocks_per_ring*(cryst_per_block + 1))/4));
        else
            yp(ii:ii + (blocks_per_ring * (cryst_per_block + 1))/4 - 1 + (cryst_per_block + 1)/2) = flip(yp(1:(blocks_per_ring*(cryst_per_block + 1))/4 + (cryst_per_block + 1)/2));
        end
        
        if mod(blocks_per_ring, 4) == 0
            yp(options.det_w_pseudo/2 + (cryst_per_block + 1) + 1:options.det_w_pseudo/2 + (cryst_per_block + 1) + options.det_w_pseudo/4) = flip(xp(1:options.det_w_pseudo/4));
            yp(options.det_w_pseudo/2 + options.det_w_pseudo/4 + 1:end) = flip(yp(options.det_w_pseudo/2 + (cryst_per_block + 1) + 1:options.det_w_pseudo/2 + (cryst_per_block + 1) + options.det_w_pseudo/4));
        else
            yp(options.det_w_pseudo/2 + (cryst_per_block + 1) + 1 : end) = flipud(abs(options.diameter - yp((cryst_per_block + 1) + 1 : options.det_w_pseudo/2)));
        end
        xp_final(:,hh) = xp;
        yp_final(:,hh) = yp;
        hh = hh + 1;
    end
    split = length(kk) : - 1 : 1;
    if mod(blocks_per_ring, 4) == 0
        for hh = 1 : length(kk)
            xp_final(ii:ii + options.det_w_pseudo/4 - 1,hh) = flip(yp_final(options.det_w_pseudo/4 + 1:options.det_w_pseudo/2,split(hh)));
            xp_final(options.det_w_pseudo/2 + (cryst_per_block + 1) + 1:options.det_w_pseudo/2 + options.det_w_pseudo/4  + (cryst_per_block + 1),hh) = (yp_final(options.det_w_pseudo/4 + (cryst_per_block + 1) + 1:options.det_w_pseudo/2 + (cryst_per_block + 1),hh));
            xp_final(options.det_w_pseudo/2 + options.det_per_ring/4 + (cryst_per_block + 1) + 1:end,hh) = yp_final((cryst_per_block + 1) + options.det_w_pseudo/2 + 1:options.det_w_pseudo/2 + options.det_w_pseudo/4,split(hh));
        end
    else
        x_jelppi = xp_final;
        for hh = 1 : length(kk)
            xp_final(ii:ii + options.det_w_pseudo/4 - 1 - (cryst_per_block + 1)/2,hh) = ...
                abs(options.diameter - flip(x_jelppi((cryst_per_block + 1) + 1 : ii - 1,(hh))));
            xp_final(ii + options.det_w_pseudo/4 - (cryst_per_block + 1)/2 : ii + options.det_w_pseudo/4 + (cryst_per_block + 1)/2 - 1,hh) = options.diameter;
        end
        x_jelppi = xp_final;
        for hh = 1 : length(kk)
            xp_final(ii + options.det_w_pseudo/4 + (cryst_per_block + 1)/2 : end,hh) = flipud(x_jelppi((cryst_per_block + 1) + 1 : ii + options.det_w_pseudo/4 - 1 - (cryst_per_block + 1)/2,(hh)));
        end
    end
    
    xp_final = (circshift(xp_final, size(xp_final,1)/2));
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

end