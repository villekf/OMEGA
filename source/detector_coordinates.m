function detector_coordinates(options)
%% Compute raw detector coordinates
% This function computes the original detector coordinates for one ring,
% both non-pseudo and pseudo cases (the latter only done if pseudo
% detectors are present).

cr_p = options.cr_p;
diameter = options.diameter;
cryst_per_block = options.cryst_per_block;
blocks_per_ring = options.blocks_per_ring;

%% Non-pseudo detectors

% All angles, starting from 90 degrees
angle = linspace(90,-270,blocks_per_ring + 1)';
% starting points
alkupistex = diameter;
alkupistey=diameter/2-(cryst_per_block/2 + 0.5)*cr_p;

% determine the gaps between adjecent blocks
if mod(blocks_per_ring,2) == 0
    widths = (cryst_per_block)*cr_p*cosd(angle(2:blocks_per_ring/2));
    widths_y = (cryst_per_block)*cr_p*sind(angle(1:(blocks_per_ring/2)));
    gaps = diameter-sum(widths);
else
    widths = (cryst_per_block)*cr_p*cosd(angle(2:floor(blocks_per_ring/2) + 1));
    widths_y = (cryst_per_block)*cr_p*sind(angle(1:floor(blocks_per_ring/2)));
    gaps = diameter-sum(widths);
end
erotus_x = diameter - sum(abs(widths));
erotus_y = diameter - sum(abs(widths_y));

if erotus_x < 0 || erotus_y < 0
    osuudet = widths/((cryst_per_block-1)*cr_p);
    osuudet = osuudet/sum(osuudet);
    gap = osuudet*gaps;
    if mod(blocks_per_ring,4) > 0
        jump = [gap(1:floor(length(gap)/2));gap(ceil(length(gap)/2):end)];
    else
        jump = [gap(1:floor(length(gap)/2));gap(round(length(gap)/2));gap(ceil(length(gap)/2):end)];
    end
    jump = jump / (sum(jump)/gaps);
    jump = jump / (sum(jump)/gaps);
    jump = [jump;-jump];
    jumpy = [jump(ceil(blocks_per_ring/4):-1:1);-jump(blocks_per_ring/2:-1:ceil(blocks_per_ring/4)+1);jump(blocks_per_ring/2 + ceil(blocks_per_ring/4):-1:blocks_per_ring/2 + 1);-jump(end:-1:blocks_per_ring/2 + ceil(blocks_per_ring/4) + 1)];
end
ii=1;
x=zeros(blocks_per_ring*cryst_per_block,1);
y=zeros(blocks_per_ring*cryst_per_block,1);

% Compute the detector coordinates of each detector (crystal) in each block
for blocks=1:blocks_per_ring
    for crystals=1:cryst_per_block
        % Moving to the next block
        if erotus_x < 0 || erotus_y < 0
            if crystals==1 && blocks>1
                x(ii)=alkupistex-jump(blocks-1);
                y(ii)=alkupistey+jumpy(blocks-1);
                
            else
                x(ii)=alkupistex-(cr_p)*cosd(angle(blocks));
                y(ii)=alkupistey+(cr_p)*sind(angle(blocks));
            end
        else
            if blocks > 1 && crystals == 1
                x(ii) = alkupistex - (cr_p * 0.5)*cosd(angle(blocks));
                y(ii) = alkupistey + (cr_p * 0.5)*sind(angle(blocks));
            else
                x(ii)=alkupistex-(cr_p)*cosd(angle(blocks));
                y(ii)=alkupistey+(cr_p)*sind(angle(blocks));
            end
        end
        alkupistex=x(ii);
        alkupistey=y(ii);
        ii=ii+1;
    end
    if erotus_x >= 0 || erotus_y >= 0
        merkki = sign(cosd(angle(blocks)));
        if merkki == 0
            merkki = sign(cosd(angle(blocks+1)));
        end
        alkupistex=alkupistex - (cr_p * 0.5)*cosd(angle(blocks)) - merkki*erotus_x/(floor(blocks_per_ring/2));
        merkki = sign(sind(angle(blocks)));
        if merkki == 0
            merkki = sign(sind(angle(blocks+1)));
        end
        alkupistey=alkupistey + (cr_p * 0.5)*sind(angle(blocks)) +  merkki*erotus_y/(floor(blocks_per_ring/2));
    end
end

% Values usually do not reach exact zero
% We force small values to exactly zero
x(abs(x)<abs(min(x))*1.1) = round(x(abs(x)<abs(min(x))*1.1));
y(abs(y)<abs(min(y))*1.1) = round(y(abs(y)<abs(min(y))*1.1));
x(abs(x)==abs(max(x))) = round(x(abs(x)==abs(max(x))));
y(abs(y)==abs(max(y))) = round(y(abs(y)==abs(max(y))));


%% Pseudo detectors

% determine if pseudo detectors are present
if ~isempty(options.pseudot) %mod(options.det_per_ring, options.Nang) ~= 0
    
    angle = linspace(90,-270,blocks_per_ring + 1)';
    % starting points
    alkupistex = diameter;
    alkupistey=diameter/2-((cryst_per_block+1)/2 + 0.5)*cr_p;
    
    % determine the gaps between adjecent blocks
    if mod(blocks_per_ring,2) == 0
        widths = (cryst_per_block+1)*cr_p*cosd(angle(2:blocks_per_ring/2));
        widths_y = (cryst_per_block+1)*cr_p*sind(angle(1:(blocks_per_ring/2)));
        gaps = diameter-sum(widths);
    else
        widths = (cryst_per_block+1)*cr_p*cosd(angle(2:floor(blocks_per_ring/2) + 1));
        widths_y = (cryst_per_block+1)*cr_p*sind(angle(1:floor(blocks_per_ring/2)));
        gaps = diameter-sum(widths);
    end
    erotus_x = diameter - sum(abs(widths));
    erotus_y = diameter - sum(abs(widths_y));
    
    ii=1;
    xp=zeros(blocks_per_ring*(cryst_per_block+1),1);
    yp=zeros(blocks_per_ring*(cryst_per_block+1),1);
    
    % Compute the detector coordinates of each detector (crystal) in each block
    for blocks=1:blocks_per_ring
        for crystals=1:cryst_per_block+1
            % Moving to the next block
            if blocks > 1 && crystals == 1
                xp(ii) = alkupistex - (cr_p * 0.5)*cosd(angle(blocks));
                yp(ii) = alkupistey + (cr_p * 0.5)*sind(angle(blocks));
            else
                xp(ii)=alkupistex-(cr_p)*cosd(angle(blocks));
                yp(ii)=alkupistey+(cr_p)*sind(angle(blocks));
            end
            alkupistex=xp(ii);
            alkupistey=yp(ii);
            ii=ii+1;
        end
        merkki = sign(cosd(angle(blocks)));
        if merkki == 0
            merkki = sign(cosd(angle(blocks+1)));
        end
        alkupistex=alkupistex - (cr_p * 0.5)*cosd(angle(blocks)) - merkki*erotus_x/(floor(blocks_per_ring/2));
        merkki = sign(sind(angle(blocks)));
        if merkki == 0
            merkki = sign(sind(angle(blocks+1)));
        end
        alkupistey=alkupistey + (cr_p * 0.5)*sind(angle(blocks)) +  merkki*erotus_y/(floor(blocks_per_ring/2));
    end
    % Values usually do not reach exact zero
    % We force small values to exactly zero
    xp(abs(xp)<abs(min(xp))*1.1) = round(xp(abs(xp)<abs(min(xp))*1.1));
    yp(abs(yp)<abs(min(yp))*1.1) = round(yp(abs(yp)<abs(min(yp))*1.1));
    xp(abs(xp)==abs(max(xp))) = round(xp(abs(xp)==abs(max(xp))));
    yp(abs(yp)==abs(max(yp))) = round(yp(abs(yp)==abs(max(yp))));
    
else
    % If pseudo detectors are not present
    xp = x;
    yp = y;
end

save([options.machine_name '_detector_coordinates.mat'], 'x', 'y', 'angle', 'xp', 'yp')

end