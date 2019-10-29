function sinogram_coordinates(options)
%% Coordinates for the sinogram detectors
% This code is used to compute the coordinates for the detector coordinates
% in sinogram space. It also provides the indexing needed for the formation
% of the sinograms from the raw list-mode data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi, Samuli Summala
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

cr_pz = options.cr_pz;
cryst_per_block = options.cryst_per_block;
blocks_per_ring = options.blocks_per_ring;
det_per_ring = options.det_per_ring;
Nang = options.Nang;
Ndist = options.Ndist;
machine_name = options.machine_name;
Nz = options.Nz;
span = options.span;
ring_difference = options.ring_difference;

%% 2D coordinates

load([machine_name '_detector_coordinates.mat'], 'xp', 'yp')

det_w_pseudo = length(xp);

mashing = det_w_pseudo/Nang/2;


ll = cryst_per_block + (det_w_pseudo-det_per_ring)/blocks_per_ring;
if (det_w_pseudo-det_per_ring) > 0
    pseudo_d = (cryst_per_block+1:cryst_per_block+1:det_w_pseudo)';
else
    pseudo_d = [];
end

% Determine the sinogram indices for each of the LOR

% Form the detector vector pair
L = zeros(sum(1:det_w_pseudo),2,'int32');
ii = 2;
jh = int32(1);
for kk = int32(1) : (det_w_pseudo)
    if kk == ll && (det_w_pseudo-det_per_ring) > 0
        if kk < (det_w_pseudo)
            ll = ll + cryst_per_block + (det_w_pseudo-det_per_ring)/blocks_per_ring;
            ii = ii + 1;
        end
    else
        L(jh:(jh + (det_w_pseudo) - kk),:) = [repelem((kk), det_w_pseudo-(kk-1))', ((kk):det_w_pseudo)'];
    end
    jh = jh + (det_w_pseudo) -kk + 1;
end
if (det_w_pseudo-det_per_ring) > 0
    L(ismember(L(:,2),pseudo_d),:) = [];
end
L(L(:,1) == 0,:) = [];

L = L - 1;

xa = max(L,[],2);
ya = min(L,[],2);

j=idivide(mod(xa+ya+det_w_pseudo/2+1,det_w_pseudo),2);

b=j+det_w_pseudo/2;

i=abs(xa-ya-det_w_pseudo/2);
for kk = 1 : length(ya)
    if (ya(kk)<j(kk)) || (b(kk)<xa(kk))
        i(kk)=-i(kk);
    end
end

% The sinogram corners need to the swapped
temppi = j*2 < -i;
temppi2 = (i <= (j-det_w_pseudo/2)*2);

temppi3 = false(det_w_pseudo - numel(pseudo_d));
temppi3(tril(true(det_w_pseudo - numel(pseudo_d)))) = temppi;
temppi = logical(temppi3 + tril(temppi3,-1)');

temppi3 = false(det_w_pseudo - numel(pseudo_d));
temppi3(tril(true(det_w_pseudo - numel(pseudo_d)))) = temppi2;
temppi2 = logical(temppi3 + tril(temppi3,-1)');

apu1 = triu(temppi);
apu3 = tril(temppi);
apu2 = triu(temppi2);
apu4 = tril(temppi2);

save([num2str(machine_name) '_swap_corners_' num2str(Nang) '.mat'], 'apu1', 'apu2','apu3','apu4');

if mod(Ndist,2) == 0
    joku12 = (i <= (Ndist/2 + min(0,options.ndist_side)) & i >= (-Ndist/2 + max(0,options.ndist_side)));
else
    joku12 = (i <= Ndist/2 & i >= (-Ndist/2));
end

j = idivide(j,det_w_pseudo/2/Nang);

i = i(joku12);
j = j(joku12);
if min(i) <= 0
    i = i + abs(min(i)) + 1;
end
j = j + 1;

L = L(joku12,:);

if det_w_pseudo > det_per_ring%det_w_pseudo/2/Nang > 1
    
    LL = zeros(sum(1:det_w_pseudo),2,'int32');
    ii = 2;
    jh = int32(1);
    for kk = int32(1) : (det_w_pseudo)
        if kk == ll && (det_w_pseudo-det_per_ring) > 0
            if kk < (det_w_pseudo)
                ll = ll + cryst_per_block + (det_w_pseudo-det_per_ring)/blocks_per_ring;
                ii = ii + 1;
            end
        else
            LL(jh:(jh + (det_w_pseudo) - kk),:) = [repelem((kk), det_w_pseudo-(kk-1))', ((kk):det_w_pseudo)'];
        end
        jh = jh + (det_w_pseudo) -kk + 1;
    end
    LL(LL(:,1) == 0,:) = [];
    
    LL = LL - 1;
    
    xa = max(LL,[],2);
    ya = min(LL,[],2);
    
    jj=idivide(mod(xa+ya+det_w_pseudo/2+1,det_w_pseudo),2);
    
    bb=jj+det_w_pseudo/2;
    
    ii=abs(xa-ya-det_w_pseudo/2);
    for kk = 1 : length(ya)
        if (ya(kk)<jj(kk)) || (bb(kk)<xa(kk))
            ii(kk)=-ii(kk);
        end
    end
    
    if mod(Ndist,2) == 0
        joku121 = (ii <= (Ndist/2 + min(0,options.ndist_side)) & ii >= (-Ndist/2 + max(0,options.ndist_side)));
    else
        joku121 = (ii <= Ndist/2 & ii >= (-Ndist/2));
    end
    
    ii = ii(joku121);
    jj = jj(joku121);
    
    jj = idivide(jj,det_w_pseudo/2/Nang);
    
    if min(ii) <= 0
        ii = ii + abs(min(ii)) + 1;
    end
    jj = jj + 1;
    
    LL = LL(joku121,:);
    
    LL = LL +1;
    
    xx3 = xp(LL(:,1));
    yy3 = yp(LL(:,1));
    xx4 = xp(LL(:,2));
    yy4 = yp(LL(:,2));
end

L = L +1;

xx1 = xp(L(:,1));
yy1 = yp(L(:,1));
xx2 = xp(L(:,2));
yy2 = yp(L(:,2));

% If pseudo detectors present, locate the gaps in the sinogram
if det_w_pseudo > det_per_ring
    
    xx = accumarray([j i], xx1, [Nang Ndist],@mean, NaN);
    
    if mashing > 1
        [C,~,ic] = unique(double([j i]),'rows');
        nOccurances = histc(ic, 1:length(C));
        
        indeksi1 = nOccurances < max(nOccurances);
        gaps = sub2ind([size(xx,1) size(xx,2)], C(indeksi1,1),C(indeksi1,2));
    else
        gaps = find(isnan(xx));
    end
    
    save([num2str(machine_name) '_pseudo_gaps_' num2str(Nang) '.mat'], 'gaps');
end

%%

% If mashing is present
% This code is still slightly inaccurate
if det_w_pseudo/2/Nang > 1

    x = accumarray([j i], xx1, [Nang Ndist],@mean, NaN);
    y = accumarray([j i], yy1, [Nang Ndist],@mean, NaN);
    x2 = accumarray([j i], xx2, [Nang Ndist],@mean, NaN);
    y2 = accumarray([j i], yy2, [Nang Ndist],@mean, NaN);
    
    x_diff = [zeros(1,size(x,2));diff(x)];
    
    thresholdx = mean(abs(x_diff(:)),'omitnan')*2;
    
    xx = accumarray([j i], xx1, [Nang Ndist],@mean, NaN);
    
    nanx = find(isnan(xx));
    
    y_diff = [zeros(1,size(y,2));diff(y)];
    
    thresholdy = mean(abs(y_diff(:)),'omitnan')*2;
        
    for kk = 1 : size(x,1)
        for ll = 1 : size(x,2)
            if any(nanx == sub2ind(size(x),kk, ll))
%                 if ll == 17 && kk == 46
%                     return
%                 end
                if abs(diff(xx3(jj==kk & ii==ll))) >= thresholdx %abs(diff([x(kk+1,ll) x(kk-1,ll)])) > threshold
                    %                     return
                    apu1 = xx3(jj==kk & ii==ll);
                    apu2 = xx4(jj==kk & ii==ll);
                    x(kk,ll) = mean([apu1(1) apu2(2)]);
                    x2(kk,ll) = mean([apu1(2) apu2(1)]);
                elseif abs(diff(xx3(jj==kk & ii==ll))) < thresholdx
                    apu1 = xx3(jj==kk & ii==ll);
                    apu2 = xx4(jj==kk & ii==ll);
                    x(kk,ll) = mean([apu1(1) apu1(2)]);
                    x2(kk,ll) = mean([apu2(2) apu2(1)]);
                else
                    x(kk,ll) = mean([x(kk+1,ll) x(kk-1,ll)]);
                    x2(kk,ll) = mean([x2(kk+1,ll) x2(kk-1,ll)]);
                end
                if abs(diff(yy3(jj==kk & ii==ll))) >= thresholdy %abs(diff([y(kk+1,ll) y(kk-1,ll)])) > threshold
                    apu1 = yy3(jj==kk & ii==ll);
                    apu2 = yy4(jj==kk & ii==ll);
                    y(kk,ll) = mean([apu1(1) apu2(2)]);
                    y2(kk,ll) = mean([apu1(2) apu2(1)]);
                elseif abs(diff(yy3(jj==kk & ii==ll))) < thresholdy
                    apu1 = yy3(jj==kk & ii==ll);
                    apu2 = yy4(jj==kk & ii==ll);
                    y(kk,ll) = mean([apu1(1) apu1(2)]);
                    y2(kk,ll) = mean([apu2(2) apu2(1)]);
                else
                    y(kk,ll) = mean([y(kk+1,ll) y(kk-1,ll)]);
                    y2(kk,ll) = mean([y2(kk+1,ll) y2(kk-1,ll)]);
                end
            end
        end
    end
    
else
    % If no mashing is done
    x = accumarray([j i], xx1, [Nang Ndist],@mean, NaN);
    y = accumarray([j i], yy1, [Nang Ndist],@mean, NaN);
    x2 = accumarray([j i], xx2, [Nang Ndist],@mean, NaN);
    y2 = accumarray([j i], yy2, [Nang Ndist],@mean, NaN);
    
end


x = [x(:) x2(:)];
y = [y(:) y2(:)];
%%


save([machine_name '_app_coordinates_' num2str(Ndist) 'x' num2str(Nang) '.mat'],'x','y','i','j','joku12');



%% Compute the 3D coordinates

dif = cr_pz/2;
% Michelogram row/column indices for each segment
p = zeros(floor((ring_difference-ceil(span/2))/span) + 2,1);
for kk = 0 : floor((ring_difference-ceil(span/2))/span)
    p(kk+2) = ceil(span/2) + span*kk;
end

z = zeros(options.TotSinos,2);

z(1:options.Nz,:) = [(0:dif:(Nz-1)*dif)' (0:dif:(Nz-1)*dif)'];        % Parallel sinograms
pp = ((span+1)/4:span:options.rings)';

% Oblique sinograms
for t=2:ceil(length(options.segment_table)/2)
    ddd(1,:) = [0 dif*2*p(t)];
    for i=1:floor(span/2)-1
        ddd(2*i+1,:) = [dif*i dif*(2*p(t)+3*i)];
        ddd(2*i,:) = [dif*(i-1) dif*(2*p(t)+3*i-1)];
    end
    
    d1 = [(2*(pp(1)-1)*dif:dif:((Nz-1)-(2*(pp(t)-1)))*dif)' ((2*(pp(t)-1))*dif:dif:((Nz-1)-2*(pp(1)-1))*dif)'];
    d = [ddd;d1;((Nz-1)*dif-fliplr(flip(ddd)))];
    z(sum(options.segment_table(1:2*t-3))+1:sum(options.segment_table(1:2*t-1)),:) = [d;fliplr(d)];
    
end

save([machine_name '_3D_coordinates_span' num2str(span) '_ringdiff' num2str(ring_difference) '_' num2str(size(z,1)) '.mat'],'z')

end
