function lor_pixel_count_prepass(options)
%% Count the number of pixels each LOR traverses
% This function counts the number of pixels that each LOR traverses. This
% is needed for methods 2, 3 and 4 or for method 1 when precompute_lor has
% been selected.

rings = options.rings;
diameter = options.diameter;
FOVa = options.FOVa;
axial_fov = options.axial_fov;
Nx = options.Nx;
Ny = options.Ny;
Nz = options.Nz;
machine_name = options.machine_name;
pseudot = options.pseudot;
NSinos = options.NSinos;
TotSinos = options.TotSinos;

% PET-laitteen halkaisija (mm)
R=double(diameter);
% Transaxial FOV (mm)
FOVa=double(FOVa);
% Axial FOV (mm)
axial_fow = double(axial_fov);
% Kristallien lukum��r�
blocks=int32(rings + length(pseudot));
% ensimm�inen kristalli
block1=int32(1);
% pikseleiden lukum��r� x ja y suunnissa
pikselikoko=int32(Nx);

NSinos = int32(NSinos);
NSlices = int32(Nz);
TotSinos = int32(TotSinos);

load([machine_name '_app_coordinates_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'x','y');
load([machine_name '_3D_coordinates_span' num2str(options.span) '_ringdiff' num2str(options.ring_difference) '_' num2str(TotSinos) '.mat'],'z')

if min(min(z)) == 0
    z = z + (axial_fow - max(max(z)))/2;
end



%
x=double(x);
y=double(y);
z_det = double(z);



% Pikselit
etaisyys=(R-FOVa)/2;
xx=linspace(etaisyys,R-etaisyys,pikselikoko+1);
yy=linspace(etaisyys,R-etaisyys,pikselikoko+1);
% zz=linspace(double(0),double(axial_fow),Nz+1);
%         if use_raw_data == false
zz=linspace(double(0),double(axial_fow),Nz+1);
%         else
%             zz=linspace(double((z_length-axial_fow)/2),double(axial_fow + (z_length-axial_fow)/2),Nz+1);
%         end
zz=zz(2*block1-1:2*blocks);

% Pikselien et�isyys toisistaan
d=diff(xx(1:2));
dz=diff(zz(1:2));

% Kuvan et�isyys origosta
bx=xx(1);
by=yy(1);
bz=zz(1);

% Pikselien m��r�
Ny=int32(Ny);
Nx=int32(Nx);
Nz=int32(Nz);

iij=double(0:Nx);
jji=double(0:Ny);
kkj=double(0:Nz);

size_x = int32(size(x,1));

zmax = max(max(z_det));
if zmax==0
    zmax = double(1);
end

% Determine which LORs go through the FOV
[ discard, lor] = improved_Siddon_algorithm_discard( TotSinos, Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , NSinos, NSlices, size_x, zmax);
clear mex
% lor = lor(discard);

save([machine_name '_lor_pixel_count_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_sino_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'lor','discard','-v7.3')

end