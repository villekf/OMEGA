function detector_pair_formation_raw(options)
%% Form valid detector pairs (LORs) for raw (list-mode) data
% This function forms the valid LORs (i.e. LORs that go through the FOV)
% for the raw list-mode data.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi
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

rings = options.rings;
det_per_ring = options.det_per_ring;
diameter = options.diameter;
FOVax = options.FOVa_x;
FOVay = options.FOVa_y;
cr_pz = options.cr_pz;
axial_fov = options.axial_fov;
Nx = options.Nx;
Ny = options.Ny;
Nz = options.Nz;
machine_name = options.machine_name;
pseudot = options.pseudot;

% Detector pair (LOR) vector

LL = form_detector_pairs_raw(rings, det_per_ring);

load([options.machine_name '_detector_coordinates.mat'], 'x', 'y')

% Parameters for (improved) Siddon algorithm

z=(rings+length(pseudot) - 1)*cr_pz;

z_det=linspace(0,(z),(rings)+length(pseudot));
if min(min(z_det)) == 0
    z_det = z_det + (axial_fov - max(max(z_det)))/2;
end

% The image (pixel) space
R = diameter;
etaisyys=(R-FOVax)/2;
xx=linspace(etaisyys,R-etaisyys,Nx+1);
etaisyys=(R-FOVay)/2;
yy=linspace(etaisyys,R-etaisyys,Ny+1);
zz=linspace(0,axial_fov,Nz+1);

% The distance of the pixels from each other
d=diff(xx(1:2));
dy=diff(yy(1:2));
dz=diff(zz(1:2));

% The distance of the image from the origo
bx=xx(1);
by=yy(1);
bz=zz(1);

L1 = LL';
L1=L1(:);

% Determine which LORs go throught the FOV
[discard, lor] = Siddon_algorithm_discard_L( int32(Nx), int32(1), int32(rings), int32(det_per_ring), d, dz, by, bx, bz, z_det, x, y, dy, yy, xx, L1, int32(Ny), int32(Nz), int32(pseudot));
clear mex L1

% Remove the LORs that do not go through the FOV
% LL = LL(discard,:);
lor = lor(discard);
save([machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'],'discard','lor','-v7.3')

end