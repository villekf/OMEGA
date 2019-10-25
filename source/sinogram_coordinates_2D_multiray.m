function [varargout] = sinogram_coordinates_2D_multiray(options)
%% Coordinates for the sinogram detectors
% This code is used to compute the 2D coordinates for the detector
% coordinates in sinogram space. It also provides the indexing needed for
% the formation of the sinograms from the raw list-mode data
%
% OUTPUTS:
%   x = X-coordinates of each of the sinogram bin in one sinogram
%   y = Y-coordinates of each of the sinogram bin in one sinogram
%   i = Sinogram bin number (distance) for each measurement
%   j = Sinogram bin number (angle) for each measurement
%   accepted_lors = The indices of the LORs that are within the specified
%   distance value (Ndist)
%   swap = Indices of sinogram corners to be swapped
%   gaps = Location of pseudo detector gaps in the sinogram
%
% See also sinogram_coordinates_3D, detector_coordinates, form_sinograms

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

% if nargout > 7
%     error('Too many output arguments')
% end
cryst_per_block = options.cryst_per_block;
blocks_per_ring = options.blocks_per_ring;
det_per_ring = options.det_per_ring;
Nang = options.Nang;
Ndist = options.Ndist;
% machine_name = options.machine_name;

% folder = fileparts(which('sinogram_coordinates_2D.m'));
% folder = strrep(folder, 'source','mat-files/');
% folder = strrep(folder, '\','/');

%% 2D coordinates

% Detector coordinates
[~, ~, xp, yp] = detector_coordinates_multiray(options);

x_final = zeros(options.Nang*options.Ndist,2,3);
y_final = zeros(options.Nang*options.Ndist,2,3);

for uu = 1 : 3

det_w_pseudo = length(xp(:,uu));

% Determine whether mashing is applied
mashing = det_w_pseudo/Nang/2;

% Pseudo rings if present
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

% if nargout >= 6
%     temppi = j*2 < -i;
%     temppi2 = (i <= (j-det_w_pseudo/2)*2);
%     
%     temppi3 = false(det_w_pseudo - numel(pseudo_d));
%     temppi3(tril(true(det_w_pseudo - numel(pseudo_d)))) = temppi;
%     temppi = logical(temppi3 + tril(temppi3,-1)');
%     
%     temppi3 = false(det_w_pseudo - numel(pseudo_d));
%     temppi3(tril(true(det_w_pseudo - numel(pseudo_d)))) = temppi2;
%     temppi2 = logical(temppi3 + tril(temppi3,-1)');
%     
%     swap1 = triu(temppi);
%     swap3 = tril(temppi);
%     swap2 = triu(temppi2);
%     swap4 = tril(temppi2);
%     varargout{6} = cat(3, swap1, swap2, swap3, swap4);
% end

% save([folder machine_name '_app_coordinates_' num2str(Ndist) 'x' num2str(Nang) '.mat'], 'swap1', 'swap2','swap3','swap4');

% Determine the accepted LORs (distances that are within the predefined
% value)
if mod(Ndist,2) == 0
    accepted_lors = (i <= (Ndist/2 + min(0,options.ndist_side)) & i >= (-Ndist/2 + max(0,options.ndist_side)));
else
    accepted_lors = (i <= Ndist/2 & i >= (-Ndist/2));
end
% if nargout >= 5
%     varargout{5} = accepted_lors;
% end

j = idivide(j,det_w_pseudo/2/Nang);

i = i(accepted_lors);
j = j(accepted_lors);
if min(i) <= 0
    i = i + abs(min(i)) + 1;
end
j = j + 1;

L = L(accepted_lors,:);

L = L + 1;

xx1 = xp(L(:,1),uu);
yy1 = yp(L(:,1),uu);
xx2 = xp(L(:,2),uu);
yy2 = yp(L(:,2),uu);

if mashing > 1
    n = 2;
    xp = arrayfun(@(i) mean(xp(i:i+n-1),uu),1:n:det_w_pseudo-n+1)';
    yp = arrayfun(@(i) mean(yp(i:i+n-1),uu),1:n:det_w_pseudo-n+1)';
    det_w_pseudo = det_w_pseudo / mashing;
    det_per_ring = det_per_ring / mashing;
    
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
    
%     jj = idivide(jj,det_w_pseudo/2);
    
    if min(ii) <= 0
        ii = ii + abs(min(ii)) + 1;
    end
    jj = jj + 1;
    
    LL = LL(joku121,:);
    
    LL = LL +1;
    
    xx3 = xp(LL(:,1),uu);
    yy3 = yp(LL(:,1),uu);
    xx4 = xp(LL(:,2),uu);
    yy4 = yp(LL(:,2),uu);
end

% If pseudo detectors present, locate the gaps in the sinogram
% if det_w_pseudo > det_per_ring && nargout >= 6
%     
%     xx = accumarray([i j], xx1, [Ndist Nang],@mean, NaN);
%     
%     if mashing > 1
%         [C,~,ic] = unique(double([j i]),'rows');
%         nOccurances = histc(ic, 1:length(C));
%         
%         indeksi1 = nOccurances < max(nOccurances);
%         gaps = sub2ind([size(xx,1) size(xx,2)], C(indeksi1,1),C(indeksi1,2));
%     else
%         gaps = find(isnan(xx));
%     end
% %     if nargout >= 6
% %         varargout{6} = gaps;
% %     end
%     
% %     save([folder machine_name '_app_coordinates_' num2str(Ndist) 'x' num2str(Nang) '.mat'], 'gaps','-append');
% end

%%

% If mashing is present
if mashing > 1
    x = accumarray([ii jj], xx3, [Ndist Nang],@mean, NaN);
    y = accumarray([ii jj], yy3, [Ndist Nang],@mean, NaN);
    x2 = accumarray([ii jj], xx4, [Ndist Nang],@mean, NaN);
    y2 = accumarray([ii jj], yy4, [Ndist Nang],@mean, NaN);
else
    % If no mashing is done
    x = accumarray([i j], xx1, [Ndist Nang],@mean, NaN);
    y = accumarray([i j], yy1, [Ndist Nang],@mean, NaN);
    x2 = accumarray([i j], xx2, [Ndist Nang],@mean, NaN);
    y2 = accumarray([i j], yy2, [Ndist Nang],@mean, NaN);
    
end


x = [x(:) x2(:)];
y = [y(:) y2(:)];

x_final(:,:,uu) = x;
y_final(:,:,uu) = y;

end

if nargout >= 1
    varargout{1} = x_final;
end
if nargout >= 2
    varargout{2} = y_final;
end
% if nargout >= 3
%     varargout{3} = i;
% end
% if nargout >= 4
%     varargout{4} = j;
% end
%%
% save([folder machine_name '_app_coordinates_' num2str(Ndist) 'x' num2str(Nang) '.mat'],'x','y','i','j','accepted_lors','-append');