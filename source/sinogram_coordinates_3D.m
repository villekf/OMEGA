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
% machine_name = options.machine_name;
Nz = options.Nz;
span = options.span;
ring_difference = options.ring_difference;

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

% save([machine_name '_3D_coordinates_span' num2str(span) '_ringdiff' num2str(ring_difference) '_' num2str(size(z,1)) '.mat'],'z')

end
