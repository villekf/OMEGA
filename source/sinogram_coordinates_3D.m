function z = sinogram_coordinates_3D(options)
%% Coordinates for the sinogram detectors
% This code is used to compute the z-coordinates for the detector
% coordinates in sinogram space. 
%
% OUTPUT:
%   z = Z-coordinates for the sinograms
%
% See also sinogram_coordinates_2D, detector_coordinates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Samuli Summala, Ville-Veikko Wettenhovi
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
Nz = options.Nz;
span = options.span;
ring_difference = options.ring_difference;

%% Compute the 3D coordinates

dif = cr_pz/2;
% dif = options.axial_fov/(options.Nz);
% Michelogram row/column indices for each segment
p = zeros(floor((ring_difference-ceil(span/2))/span) + 2,1);
for kk = 0 : floor((ring_difference-ceil(span/2))/span)
    p(kk+2) = ceil(span/2)*2 + (span*2)*kk;
end

z = zeros(options.NSinos,2);

% Parallel sinograms
if exist('OCTAVE_VERSION', 'builtin') == 0
    z(1:options.Nz,:) = [round((0:dif:(Nz-1)*dif)',5) round((0:dif:(Nz-1)*dif)',5)];
else
    z(1:options.Nz,:) = [round((0:dif:(Nz-1)*dif)'.*10e5)./10e5 round((0:dif:(Nz-1)*dif)'.*10e5)./10e5];
end

% Oblique sinograms
for t=2:ceil(length(options.segment_table)/2)
    % If the last segment is only partially filled, take less initial
    % sinograms
    if floor(options.segment_table(t*2-1)/2) < span-2
        maara = floor((options.segment_table(t*2-1) - 3)/4);
    else
        maara = floor(span/2) - 1;
    end
    % There's uneven amount of initial sinograms that do not yet follow the
    % floor(span/2) -> ceil(span/2) pattern
    ddd = zeros(maara*2+1,2);
    % The first element is present in all cases (span >= 3)
    ddd(1,:) = [0 dif*p(t)];
    for i=1:maara
        ddd(2*i,:) = [dif*(i-1) dif*(p(t)+i*2+(i-1))];
        ddd(2*i+1,:) = [dif*i dif*(p(t)+i*2+1+(i-1))];
    end
    % Mirror image
    if exist('OCTAVE_VERSION', 'builtin') == 0
        ddd2 = round((Nz-1)*dif-rot90(ddd,2), 5);
    else
        ddd2 = round((Nz-1)*dif-rot90(ddd,2).*10e5)./10e5;
    end
    % Sinograms in the middle of the intial sinograms
    % There'a always a minimum of 3 of these
    d1 = [(ddd(end,1):dif:ddd2(1,1)-2*dif)' (ddd(end,2)+dif*2:dif:ddd2(1,2))'];
    if exist('OCTAVE_VERSION', 'builtin') == 0
        d = round([ddd;d1;ddd2],5);
    else
        d = round([ddd;d1;ddd2].*10e5)./10e5;
    end
    z(sum(options.segment_table(1:2*t-3))+1:sum(options.segment_table(1:2*t-1)),:) = [d;fliplr(d)];
end

end
