function [varargout] = sinogram_coordinates_2D(options, xp, yp)
%% Coordinates for the sinogram detectors
% This code is used to compute the 2D coordinates for the detector
% coordinates in sinogram space. It also provides the indexing needed for
% the formation of the sinograms from the raw list-mode data.
%
% INPUTS:
%   options = Machine and sinogram properties
%   x/y = Detector coordinates in x- and y-directions
%
% OUTPUTS:
%   x = X-coordinates of each sinogram bin in one sinogram
%   y = Y-coordinates of each sinogram bin in one sinogram
%   i = Sinogram bin number (distance) for each measurement
%   j = Sinogram bin number (angle) for each measurement
%   accepted_lors = The indices of the LORs that are within the specified
%   distance value (Ndist)
%   swap = Indices of sinogram corners to be swapped
%
% See also sinogram_coordinates_3D, detector_coordinates, form_sinograms


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. 
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details. 
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout > 6
    error('Too many output arguments')
end
det_w_pseudo = options.det_w_pseudo;
Nang = options.Nang;
Ndist = options.Ndist;


%% 2D coordinates

% Determine the sinogram indices for each of the LOR

% Form the detector vector pair
L = zeros(sum(1:det_w_pseudo),2,'int32');
jh = int32(1);
for kk = int32(1) : (det_w_pseudo)
    L(jh:(jh + (det_w_pseudo) - kk),:) = [repelem((kk), det_w_pseudo-(kk-1))', ((kk):det_w_pseudo)'];
    jh = jh + (det_w_pseudo) -kk + 1;
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

if nargout >= 6
    temppi = j*2 < -i;
    temppi2 = (i <= (j-det_w_pseudo/2)*2);
    
    temppi3 = false(det_w_pseudo);
    temppi3(tril(true(det_w_pseudo))) = temppi;
    temppi = logical(temppi3 + tril(temppi3,-1)');
    
    temppi3 = false(det_w_pseudo);
    temppi3(tril(true(det_w_pseudo))) = temppi2;
    temppi2 = logical(temppi3 + tril(temppi3,-1)');
    
    swap1 = triu(temppi);
    swap3 = tril(temppi);
    swap2 = triu(temppi2);
    swap4 = tril(temppi2);
    varargout{6} = cat(3, swap1, swap2, swap3, swap4);
end

% Determine the accepted LORs (distances that are within the predefined
% value)
if mod(Ndist,2) == 0
    accepted_lors = (i <= (Ndist/2 + min(0,options.ndist_side)) & i >= (-Ndist/2 + max(0,options.ndist_side)));
else
    accepted_lors = (i <= Ndist/2 & i >= (-Ndist/2));
end
if nargout >= 5
    varargout{5} = accepted_lors;
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

xx1 = xp(L(:,1));
yy1 = yp(L(:,1));
xx2 = xp(L(:,2));
yy2 = yp(L(:,2));

%%

x = accumarray([i j], xx1, [Ndist Nang],@mean, NaN);
y = accumarray([i j], yy1, [Ndist Nang],@mean, NaN);
x2 = accumarray([i j], xx2, [Ndist Nang],@mean, NaN);
y2 = accumarray([i j], yy2, [Ndist Nang],@mean, NaN);

if sum(isnan(x)) > 0
    x = fillmissing(x,'linear');
    y = fillmissing(y,'linear');
    x2 = fillmissing(x2,'linear');
    y2 = fillmissing(y2,'linear');
end


x = [x(:) x2(:)];
y = [y(:) y2(:)];

if nargout >= 1
    varargout{1} = x;
end
if nargout >= 2
    varargout{2} = y;
end
if nargout >= 3
    varargout{3} = i;
end
if nargout >= 4
    varargout{4} = j;
end