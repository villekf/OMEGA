function detector_pair_formation_raw(options)
%% Form valid detector pairs (LORs) for raw (list-mode) data
% This function forms the valid LORs (i.e. LORs that go through the FOV)
% for the raw list-mode data.

rings = options.rings;
det_per_ring = options.det_per_ring;
diameter = options.diameter;
FOVa = options.FOVa;
cr_pz = options.cr_pz;
axial_fov = options.axial_fov;
Nx = options.Nx;
Ny = options.Ny;
Nz = options.Nz;
machine_name = options.machine_name;
pseudot = options.pseudot;

% Detector pair (LOR) vector

temp = uint16(rings*det_per_ring);

LL = (uint16(1):temp)';
LL = repelem(LL, flip(LL));

apu = zeros(size(LL),'uint16');

loc = uint32(1);
for kk = uint32(1) : uint32(temp)
    apu(loc: loc + (uint32(temp(end)) - kk)) = (uint16(kk) :  temp)';
    loc = loc + (uint32(temp(end)) - kk) + 1;
end

LL = [LL,apu];
clear apu

load([options.machine_name '_detector_coordinates.mat'], 'x', 'y')

% Parameters for (improved) Siddon algorithm

z=(rings+length(pseudot) - 1)*cr_pz;

z_det=linspace(0,(z),(rings)+length(pseudot));
if min(min(z_det)) == 0
    z_det = z_det + (axial_fov - max(max(z_det)))/2;
end

% The image (pixel) space
R = diameter;
etaisyys=(R-FOVa)/2;
xx=linspace(etaisyys,R-etaisyys,Nx+1);
yy=linspace(etaisyys,R-etaisyys,Ny+1);
zz=linspace(0,axial_fov,Nz+1);

% The distance of the pixels from each other
d=diff(xx(1:2));
dz=diff(zz(1:2));

% The distance of the image from the origo
bx=xx(1);
by=yy(1);
bz=zz(1);

iij=([0,Nx]);
jji=([0,Ny]);
kkj=([0,Nz]);

L1 = LL';
L1=L1(:);

% Determine which LORs go throught the FOV
[discard, lor] = Siddon_algorithm_discard_L( int32(Nx), int32(1), int32(rings), int32(det_per_ring), d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx, L1, int32(Ny), int32(Nz), int32(pseudot));
clear mex L1

% Remove the LORs that do not go through the FOV
LL = LL(discard,:);
lor = lor(discard);
save([machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'],'LL','discard','lor','-v7.3')

end