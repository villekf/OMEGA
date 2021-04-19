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

% cryst_per_block = options.cryst_per_block;
% blocks_per_ring = options.blocks_per_ring;
% det_per_ring = options.det_per_ring;
Nang = options.Nang;
Ndist = options.Ndist;
mashing = 1;
if Nang < options.det_w_pseudo/2
    mashing = (options.det_w_pseudo / Nang / 2);
    Nang = Nang * mashing;
end

%% 2D coordinates

% Detector coordinates
% [~, ~, xp, yp] = detector_coordinates_multiray(options);
[xp, yp] = detector_coordinates_multiray(options);

x_final = zeros(options.Nang*options.Ndist,2,options.n_rays_transaxial);
y_final = zeros(options.Nang*options.Ndist,2,options.n_rays_transaxial);

for uu = 1 : options.n_rays_transaxial
    
    det_w_pseudo = options.det_w_pseudo;
    
    % Pseudo rings if present
%     ll = cryst_per_block + (det_w_pseudo-det_per_ring)/blocks_per_ring;
%     if (det_w_pseudo-det_per_ring) > 0
%         pseudo_d = (cryst_per_block+1:cryst_per_block+1:det_w_pseudo)';
%     else
%         pseudo_d = [];
%     end
    
    % Determine the sinogram indices for each of the LOR
    
    % Form the detector vector pair
    L = zeros(sum(1:det_w_pseudo),2,'int32');
    L = reshape(L, numel(L)/2, 2);
    jh = int32(1);
    for kk = int32(1) : (det_w_pseudo)
%         if kk == ll && (det_w_pseudo-det_per_ring) > 0
%             if kk < (det_w_pseudo)
%                 ll = ll + cryst_per_block + (det_w_pseudo-det_per_ring)/blocks_per_ring;
%                 ii = ii + 1;
%             end
%         else
            L(jh:(jh + (det_w_pseudo) - kk),:) = [reshape(repelem((kk), det_w_pseudo-(kk-1)), [], 1), ((kk):det_w_pseudo)'];
%         end
        jh = jh + (det_w_pseudo) -kk + 1;
    end
%     if (det_w_pseudo-det_per_ring) > 0
%         L(ismember(L(:,2),pseudo_d),:) = [];
%     end
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
    
    % Determine the accepted LORs (distances that are within the predefined
    % value)
    if mod(Ndist,2) == 0
        accepted_lors = (i <= (Ndist/2 + min(0,options.ndist_side)) & i >= (-Ndist/2 + max(0,options.ndist_side)));
    else
        accepted_lors = (i <= Ndist/2 & i >= (-Ndist/2));
    end
    
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
    
    %%
    
    x = accumarray([i j], xx1, [Ndist Nang],@mean, NaN);
    y = accumarray([i j], yy1, [Ndist Nang],@mean, NaN);
    x2 = accumarray([i j], xx2, [Ndist Nang],@mean, NaN);
    y2 = accumarray([i j], yy2, [Ndist Nang],@mean, NaN);
    
    if mashing > 1
        xx1 = reshape(x,options.Ndist,Nang);
        xx2 = reshape(x2,options.Ndist,Nang);
        yy1 = reshape(y,options.Ndist,Nang);
        yy2 = reshape(y2,options.Ndist,Nang);
        
        testi = abs(diff(xx1));
        [I,J] = find(testi > options.cr_p*2);
        testi2 = abs(diff(J));
        ind1 = find(testi2 > mean(testi2)*2, 1, 'first');
        if xx1(1,J(ind1)) <= xx1(1,J(ind1) + 1)
            indices2 = J(ind1) + 1 : - 1 : 1;
            indices2 = [indices2(1);repelem(indices2(2:end),mashing)'];
        elseif xx1(1,J(ind1)) <= xx1(2,J(ind1))
            indices2 = J(ind1) : - 1 : 1;
            indices2 = [indices2(1);indices2(1);repelem(indices2(2:end),mashing)'];
        else
            indices2 = J(ind1) : - 1 : 1;
            indices2 = [indices2(1);repelem(indices2(2:end),mashing)'];
        end
        indices1 = false(size(xx1));
        for kk = 1 : I(1)
            indices1(kk,1:indices2(kk)) = true(indices2(kk),1);
        end
        testi = abs(diff(fliplr(xx1)));
        [I,J] = find(testi > options.cr_p*2);
        testi2 = abs(diff(J));
        ind1 = find(testi2 > mean(testi2)*2, 1, 'first');
        if xx1(1,Nang - J(ind1) + 1) >= xx1(2,Nang - J(ind1) + 1)
            indices2 = Nang - J(ind1) + 1 : Nang;
            indices2 = [indices2(1);indices2(1);repelem(indices2(2:end),mashing)'];
        elseif xx1(1,Nang - J(ind1) + 1) >= xx1(1,Nang - J(ind1))
            indices2 = Nang - J(ind1) + 1 : Nang;
            indices2 = [indices2(1);repelem(indices2(2:end),mashing)'];
        else
            indices2 = Nang - J(ind1) + 1 : Nang;
            indices2 = [indices2(1);repelem(indices2(2:end),mashing)'];
        end
        for kk = 1 : I(1)
            indices1(kk,indices2(kk):end) = true(length(indices1(kk,indices2(kk):end)),1);
        end
        temp = xx1(indices1);
        x(indices1) = xx2(indices1);
        x2(indices1) = temp;
        temp = yy1(indices1);
        y(indices1) = yy2(indices1);
        y2(indices1) = temp;
        
        x = cell2mat(arrayfun(@(i) mean(x(:,i:i+mashing-1),2),1:mashing:size(x,2)-mashing+1,'UniformOutput',false));
        x2 = cell2mat(arrayfun(@(i) mean(x2(:,i:i+mashing-1),2),1:mashing:size(x2,2)-mashing+1,'UniformOutput',false));
        y = cell2mat(arrayfun(@(i) mean(y(:,i:i+mashing-1),2),1:mashing:size(y,2)-mashing+1,'UniformOutput',false));
        y2 = cell2mat(arrayfun(@(i) mean(y2(:,i:i+mashing-1),2),1:mashing:size(y2,2)-mashing+1,'UniformOutput',false));
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