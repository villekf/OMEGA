function [x,y] = computeCoordinates(blocks_per_ring, transaxial_multip,cryst_per_block, diameter, cr_p, usePseudo)
%COMPUTECOORDINATES Computes the coordinates for each detector on a single
%ring
%   This function computes the transaxial detector coordinates (x,y) for
%   all crystals on one ring. Can be used also for cases with pseudo
%   detectors.
%
% See also detector_coordinates

cryst_per_block_orig = cryst_per_block;
if usePseudo
    extraVar = 1;
    cryst_per_block = cryst_per_block + 1;
else
    extraVar = 0;
end
angle = (90:-360/(blocks_per_ring * transaxial_multip):-270)';

% Width of the topmost blocks (diameter)
widths = (cryst_per_block_orig) * cr_p * cosd(angle(1:((blocks_per_ring * transaxial_multip)/2 + 1)));
widthsy = (cryst_per_block_orig) * cr_p * sind(angle(1:((blocks_per_ring * transaxial_multip)/2)));

% Gap between adjacent blocks
% If negative, then the crystal size is smaller on the edges of the block
% than the crystal pitch
% diameter = orig_diameter + DOI * 2;
erotus = (diameter - sum(abs(widths)))/sum(cosd(angle(1:((blocks_per_ring * transaxial_multip)/2 + 1))))/2;
erotusy = (diameter - sum(abs(widthsy)))/sum(abs(sind(angle(1:((blocks_per_ring * transaxial_multip)/2)))))/2;

% starting points
alkupistex = (diameter/2);
alkupistey = - ((cryst_per_block_orig)/2 + 0.5) * cr_p;

ii = 1;
x = zeros((blocks_per_ring * transaxial_multip) * (cryst_per_block),1);
y = zeros((blocks_per_ring * transaxial_multip) * (cryst_per_block),1);
%%
% Compute the detector coordinates of each detector (crystal) in each block
% Only for the 1/4th of the ring
for blocks = 1:ceil((blocks_per_ring * transaxial_multip)/4)
    for crystals = 1:cryst_per_block
        % Moving to the next block
        if blocks > 1 && crystals == 1
            x(ii) = alkupistex - (cr_p * 0.5) * cosd(angle(blocks)) - erotus * cosd(angle(blocks - 1)) - erotus * cosd(angle(blocks)) - (cr_p * 0.5) * cosd(angle(blocks - 1));
            y(ii) = alkupistey + (cr_p * 0.5) * sind(angle(blocks)) + erotusy * sind(angle(blocks - 1)) + erotusy * sind(angle(blocks)) + (cr_p * 0.5) * sind(angle(blocks - 1));
        else
            % While in the same block
            x(ii) = alkupistex - (cr_p) * cosd(angle(blocks));
            y(ii) = alkupistey + (cr_p) * sind(angle(blocks));
        end
        if crystals == cryst_per_block
            alkupistex = x(ii) + (cr_p) * cosd(angle(blocks)) * extraVar;
            alkupistey = y(ii) - (cr_p) * sind(angle(blocks)) * extraVar;
        else
            alkupistex = x(ii);
            alkupistey = y(ii);
        end
        ii=ii+1;
    end
end
%%
% This section ensures symmetry of the coordinates
% Fills the rest of coordinates (3/4)

% If the blocks can be divided by 4, compute an extra block (i.e. the top
% block) first
if mod((blocks_per_ring * transaxial_multip),4) == 0
    blocks = blocks + 1;
    x(ii) = alkupistex - (cr_p * 0.5) * cosd(angle(blocks)) - erotus * cosd(angle(blocks - 1)) - erotus * cosd(angle(blocks)) - (cr_p * 0.5) * cosd(angle(blocks - 1));
    y(ii) = alkupistey + (cr_p * 0.5) * sind(angle(blocks)) + erotusy * sind(angle(blocks - 1)) + erotusy * sind(angle(blocks)) + (cr_p * 0.5) * sind(angle(blocks - 1));
    alkupistex = x(ii);
    alkupistey = y(ii);
    ii = ii + 1;
    for ll = 2:cryst_per_block
        x(ii) = alkupistex - (cr_p) * cosd(angle(blocks));
        y(ii) = alkupistey + (cr_p) * sind(angle(blocks));
        alkupistex = x(ii);
        alkupistey = y(ii);
        ii = ii + 1;
    end
    % Symmmetry
    if usePseudo
        ii = ii - 1;
        alkupistex = x(ii) + (cr_p) * cosd(angle(blocks)) * extraVar;
        alkupistey = y(ii) - (cr_p) * sind(angle(blocks)) * extraVar;
        alkublock = blocks + 1;
        ii = ii + 1;
        for blocks = alkublock:blocks_per_ring * transaxial_multip
            for crystals = 1:cryst_per_block
                % Moving to the next block
                if blocks > 1 && crystals == 1
                    x(ii) = alkupistex - (cr_p * 0.5) * cosd(angle(blocks)) - erotus * cosd(angle(blocks - 1)) - erotus * cosd(angle(blocks)) - (cr_p * 0.5)*cosd(angle(blocks - 1));
                    y(ii) = alkupistey + (cr_p * 0.5) * sind(angle(blocks)) + erotusy * sind(angle(blocks - 1)) + erotusy * sind(angle(blocks)) + (cr_p * 0.5)*sind(angle(blocks - 1));
                else
                    % While in the same block
                    x(ii) = alkupistex - (cr_p) * cosd(angle(blocks));
                    y(ii) = alkupistey + (cr_p) * sind(angle(blocks));
                end
                if crystals == cryst_per_block
                    alkupistex = x(ii) + (cr_p) * cosd(angle(blocks)) * extraVar;
                    alkupistey = y(ii) - (cr_p) * sind(angle(blocks)) * extraVar;
                else
                    alkupistex = x(ii);
                    alkupistey = y(ii);
                end
                ii = ii + 1;
            end
        end
    else
        x(ii:ii + ((blocks_per_ring * transaxial_multip)*cryst_per_block)/4 - 1) = -(flip(x(1:((blocks_per_ring * transaxial_multip)*cryst_per_block)/4)));
        y(ii:ii + ((blocks_per_ring * transaxial_multip)*cryst_per_block)/4 - 1) = (flip(y(1:((blocks_per_ring * transaxial_multip)*cryst_per_block)/4)));
        x(ii + ((blocks_per_ring * transaxial_multip)*cryst_per_block)/4:end) = flip(x(cryst_per_block + 1:((blocks_per_ring * transaxial_multip)*cryst_per_block)/2));
        y(ii + ((blocks_per_ring * transaxial_multip)*cryst_per_block)/4:end) = -(y(cryst_per_block + 1:((blocks_per_ring * transaxial_multip)*cryst_per_block)/2));
    end
else
    if usePseudo
        alkublock = blocks + 1;
        for blocks = alkublock:blocks_per_ring * transaxial_multip
            for crystals = 1:cryst_per_block
                % Moving to the next block
                if blocks > 1 && crystals == 1
                    x(ii) = alkupistex - (cr_p * 0.5) * cosd(angle(blocks)) - erotus * cosd(angle(blocks - 1)) - erotus * cosd(angle(blocks)) - (cr_p * 0.5)*cosd(angle(blocks - 1));
                    y(ii) = alkupistey + (cr_p * 0.5) * sind(angle(blocks)) + erotusy * sind(angle(blocks - 1)) + erotusy * sind(angle(blocks)) + (cr_p * 0.5)*sind(angle(blocks - 1));
                else
                    % While in the same block
                    x(ii) = alkupistex - (cr_p) * cosd(angle(blocks));
                    y(ii) = alkupistey + (cr_p) * sind(angle(blocks));
                end
                if crystals == cryst_per_block
                    alkupistex = x(ii) + (cr_p) * cosd(angle(blocks)) * extraVar;
                    alkupistey = y(ii) - (cr_p) * sind(angle(blocks)) * extraVar;
                else
                    alkupistex = x(ii);
                    alkupistey = y(ii);
                end
                ii = ii + 1;
            end
        end
    else
        x(ii:ii + ((blocks_per_ring * transaxial_multip)*cryst_per_block)/4 + cryst_per_block/2-1) = -(flip(x(1:((blocks_per_ring * transaxial_multip)*cryst_per_block)/4 + cryst_per_block/2)));
        y(ii:ii + ((blocks_per_ring * transaxial_multip)*cryst_per_block)/4 + cryst_per_block/2-1) = (flip(y(1:((blocks_per_ring * transaxial_multip)*cryst_per_block)/4 + cryst_per_block/2)));
        x(ii + ((blocks_per_ring * transaxial_multip)*cryst_per_block)/4 + cryst_per_block/2:end) = flip(x(cryst_per_block + 1:((blocks_per_ring * transaxial_multip)*cryst_per_block)/2));
        y(ii + ((blocks_per_ring * transaxial_multip)*cryst_per_block)/4 + cryst_per_block/2:end) = -(y(cryst_per_block + 1:((blocks_per_ring * transaxial_multip)*cryst_per_block)/2));
    end
end
end