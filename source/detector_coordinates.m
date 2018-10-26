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
% if mod(cryst_per_block,2) == 0
    alkupistey=diameter/2-(cryst_per_block/2 + 0.5)*cr_p;
% else
%     alkupistey=diameter/2-ceil(cryst_per_block/2)*cr_p;
% end
% determine the gaps between adjecent blocks
if mod(blocks_per_ring,2) == 0
    widths = (cryst_per_block)*cr_p*cosd(angle(2:blocks_per_ring/2));
    widths_y = (cryst_per_block)*cr_p*sind(angle(1:(blocks_per_ring/2)));
    gaps = diameter-sum(widths);
    %     jump = gaps/(blocks_per_ring/2);
else
    widths = (cryst_per_block)*cr_p*cosd(angle(2:floor(blocks_per_ring/2) + 1));
    widths_y = (cryst_per_block)*cr_p*sind(angle(1:floor(blocks_per_ring/2)));
    gaps = diameter-sum(widths);
    %     jump = gaps/(floor(blocks_per_ring/2) + 0.5);
end
erotus_x = diameter - sum(abs(widths));
erotus_y = diameter - sum(abs(widths_y));
% jumpy = jump;

% if mod(blocks_per_ring,2) == 0
%     widths(length(widths)/2) = [];
% end

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
    % jump = gaps/2;
    jumpy = [jump(ceil(blocks_per_ring/4):-1:1);-jump(blocks_per_ring/2:-1:ceil(blocks_per_ring/4)+1);jump(blocks_per_ring/2 + ceil(blocks_per_ring/4):-1:blocks_per_ring/2 + 1);-jump(end:-1:blocks_per_ring/2 + ceil(blocks_per_ring/4) + 1)];
end
ii=1;
x=zeros(blocks_per_ring*cryst_per_block,1);
y=zeros(blocks_per_ring*cryst_per_block,1);
% blocks = 1;
% alkupistex = alkupistex+(cr_p)*cosd(angle(blocks));
% alkupistey = alkupistey-(cr_p)*sind(angle(blocks));
% Compute the detector coordinates of each detector (crystal) in each block
for blocks=1:blocks_per_ring
    for crystals=1:cryst_per_block
        % Moving to the next block
%         if crystals==1 && blocks>1
%             x(ii)=alkupistex-jump(blocks-1);
%             y(ii)=alkupistey+jumpy(blocks-1);
%             %             if sign(cosd(angle(blocks))) ~= 0 && sign(cosd(angle(blocks))) ~= sign(cosd(angle(blocks - 1))) && blocks > 2
%             %                 jump = -jump;
%             %             end
%             %             if sign(sind(angle(blocks))) ~= 0 && sign(sind(angle(blocks))) ~= sign(sind(angle(blocks - 1)))
%             %                 jumpy = -jumpy;
%             %             end
%             %             if abs(cosd(angle(blocks))) == abs(cosd(angle(blocks - 1))) && sign(cosd(angle(blocks))) ~= sign(cosd(angle(blocks - 1)))
%             %             	x(ii)=alkupistex;
%             %             else
%             %                 x(ii)=alkupistex-jump;
%             %             end
%             %             if abs(sind(angle(blocks))) == abs(sind(angle(blocks - 1))) && sign(sind(angle(blocks))) ~= sign(sind(angle(blocks - 1)))
%             %                 y(ii)=alkupistey;
%             %             else
%             %                 y(ii)=alkupistey+jumpy;
%             %             end
%         else
            % While still in the same block
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
        %         end
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
% x(abs(x)<abs(min(x))*2) = x(abs(x)<abs(min(x))*2) - min(x);
% y(abs(y)<abs(min(y))*2) = y(abs(y)<abs(min(y))*2) - min(y);
x(abs(x)<abs(min(x))*1.1) = round(x(abs(x)<abs(min(x))*1.1));
y(abs(y)<abs(min(y))*1.1) = round(y(abs(y)<abs(min(y))*1.1));
x(abs(x)==abs(max(x))) = round(x(abs(x)==abs(max(x))));
y(abs(y)==abs(max(y))) = round(y(abs(y)==abs(max(y))));


%% Pseudo detectors

% determine if pseudo detectors are present
if ~isempty(options.pseudot) %mod(options.det_per_ring, options.Nang) ~= 0
    
    %     angle = linspace(90,-270,blocks_per_ring+1)';
    %     alkupistex=diameter;
    %     if mod(cryst_per_block,2) == 0
    %         alkupistey=diameter/2-(cryst_per_block/2 - 0.5)*cr_p;
    %     else
    %         alkupistey=diameter/2-floor(cryst_per_block/2)*cr_p;
    %     end
    %     widths = (cryst_per_block)*cr_p*cosd(angle(2:blocks_per_ring/2));
    %     gaps = diameter-sum(widths);
    %     osuudet = widths/(cryst_per_block*cr_p);
    %     osuudet = osuudet/sum(osuudet);
    %     gap = osuudet*gaps;
    %     jump = [gap(1:floor(length(gap)/2));gap(floor(length(gap)/2));gap(ceil(length(gap)/2):end)];
    %     jump = jump / (sum(jump)/gaps);
    %     jump = [jump;-jump];
    %     jumpy = [jump(blocks_per_ring/4:-1:1);-jump(blocks_per_ring/2:-1:blocks_per_ring/4+1);jump(blocks_per_ring/2 + blocks_per_ring/4:-1:blocks_per_ring/2 + 1);-jump(end:-1:blocks_per_ring/2 + blocks_per_ring/4 + 1)];
    %     ii=1;
    %     xp=zeros(blocks_per_ring*(cryst_per_block+1),1);
    %     yp=zeros(blocks_per_ring*(cryst_per_block+1),1);
    %     blocks = 1;
    %     alkupistex = alkupistex+(cr_p)*cosd(angle(blocks));
    %     alkupistey = alkupistey-(cr_p)*sind(angle(blocks));
    %     kk = 0;
    %     for blocks=1:blocks_per_ring
    %         for crystals=1:(cryst_per_block+1)
    %             if crystals==1 && blocks>1
    %                 xx = x(ii - kk);
    %                 yy = y(ii - kk);
    %                 xp(ii) = xx;
    %                 yp(ii) = yy;
    %             elseif crystals==(cryst_per_block+1)
    %                 xp(ii)=alkupistex-jump(blocks);
    %                 yp(ii)=alkupistey+jumpy(blocks);
    %                 kk = kk + 1;
    %                 xx = xp(ii);
    %                 yy = yp(ii);
    %             else
    %                 xx= x(ii - kk);
    %                 yy= y(ii - kk);
    %                 xp(ii) = xx;
    %                 yp(ii) = yy;
    %             end
    %             alkupistex=xx;
    %             alkupistey=yy;
    %             ii=ii+1;
    %         end
    %     end
    %     xp(xp<min(xp)*2) = xp(xp<min(xp)*2) - min(xp);
    %     yp(yp<min(yp)*2) = yp(yp<min(yp)*2) - min(yp);
    
    angle = linspace(90,-270,blocks_per_ring + 1)';
    % starting points
    alkupistex = diameter;
% if mod(cryst_per_block,2) == 0
    alkupistey=diameter/2-((cryst_per_block+1)/2 + 0.5)*cr_p;
% else
%     alkupistey=diameter/2-floor((cryst_per_block+1)/2)*cr_p;
% end
% determine the gaps between adjecent blocks
if mod(blocks_per_ring,2) == 0
    widths = (cryst_per_block+1)*cr_p*cosd(angle(2:blocks_per_ring/2));
    widths_y = (cryst_per_block+1)*cr_p*sind(angle(1:(blocks_per_ring/2)));
    gaps = diameter-sum(widths);
    %     jump = gaps/(blocks_per_ring/2);
else
    widths = (cryst_per_block+1)*cr_p*cosd(angle(2:floor(blocks_per_ring/2) + 1));
    widths_y = (cryst_per_block+1)*cr_p*sind(angle(1:floor(blocks_per_ring/2)));
    gaps = diameter-sum(widths);
    %     jump = gaps/(floor(blocks_per_ring/2) + 0.5);
end
erotus_x = diameter - sum(abs(widths));
erotus_y = diameter - sum(abs(widths_y));
% jumpy = jump;

% if mod(blocks_per_ring,2) == 0
%     widths(length(widths)/2) = [];
% end

if erotus_x < 0 || erotus_y < 0
    osuudet = widths/((cryst_per_block)*cr_p);
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
    % jump = gaps/2;
    jumpy = [jump(ceil(blocks_per_ring/4):-1:1);-jump(blocks_per_ring/2:-1:ceil(blocks_per_ring/4)+1);jump(blocks_per_ring/2 + ceil(blocks_per_ring/4):-1:blocks_per_ring/2 + 1);-jump(end:-1:blocks_per_ring/2 + ceil(blocks_per_ring/4) + 1)];
end
ii=1;
xp=zeros(blocks_per_ring*(cryst_per_block+1),1);
yp=zeros(blocks_per_ring*(cryst_per_block+1),1);
% blocks = 1;
% alkupistex = alkupistex+(cr_p)*cosd(angle(blocks));
% alkupistey = alkupistey-(cr_p)*sind(angle(blocks));
% Compute the detector coordinates of each detector (crystal) in each block
for blocks=1:blocks_per_ring
    for crystals=1:cryst_per_block+1
        % Moving to the next block
%         if crystals==1 && blocks>1
%             x(ii)=alkupistex-jump(blocks-1);
%             y(ii)=alkupistey+jumpy(blocks-1);
%             %             if sign(cosd(angle(blocks))) ~= 0 && sign(cosd(angle(blocks))) ~= sign(cosd(angle(blocks - 1))) && blocks > 2
%             %                 jump = -jump;
%             %             end
%             %             if sign(sind(angle(blocks))) ~= 0 && sign(sind(angle(blocks))) ~= sign(sind(angle(blocks - 1)))
%             %                 jumpy = -jumpy;
%             %             end
%             %             if abs(cosd(angle(blocks))) == abs(cosd(angle(blocks - 1))) && sign(cosd(angle(blocks))) ~= sign(cosd(angle(blocks - 1)))
%             %             	x(ii)=alkupistex;
%             %             else
%             %                 x(ii)=alkupistex-jump;
%             %             end
%             %             if abs(sind(angle(blocks))) == abs(sind(angle(blocks - 1))) && sign(sind(angle(blocks))) ~= sign(sind(angle(blocks - 1)))
%             %                 y(ii)=alkupistey;
%             %             else
%             %                 y(ii)=alkupistey+jumpy;
%             %             end
%         else
            % While still in the same block
%             if erotus_x < 0 || erotus_y < 0
%                 if crystals==1 && blocks>1
%                     xp(ii)=alkupistex-jump(blocks-1);
%                     yp(ii)=alkupistey+jumpy(blocks-1);
%                 else
%                     xp(ii)=alkupistex-(cr_p)*cosd(angle(blocks));
%                     yp(ii)=alkupistey+(cr_p)*sind(angle(blocks));
%                 end
%             else
                if blocks > 1 && crystals == 1
                    xp(ii) = alkupistex - (cr_p * 0.5)*cosd(angle(blocks));
                    yp(ii) = alkupistey + (cr_p * 0.5)*sind(angle(blocks));
                else
                    xp(ii)=alkupistex-(cr_p)*cosd(angle(blocks));
                    yp(ii)=alkupistey+(cr_p)*sind(angle(blocks));
                end
%             end
        %         end
        alkupistex=xp(ii);
        alkupistey=yp(ii);
        ii=ii+1;
    end
%     if erotus_x >= 0 || erotus_y >= 0
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
%     end
end    
    % Values usually do not reach exact zero
    % We force small values to exactly zero
    % xp(abs(xp)<abs(min(xp))*2) = xp(abs(xp)<abs(min(xp))*2) - min(xp);
    % yp(abs(yp)<abs(min(yp))*2) = yp(abs(yp)<abs(min(yp))*2) - min(yp);
%     xp(abs(xp)<abs(min(xp))*2) = round(xp(abs(xp)<abs(min(xp))*2), 6);
%     yp(abs(yp)<abs(min(yp))*2) = round(yp(abs(yp)<abs(min(yp))*2), 6);
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