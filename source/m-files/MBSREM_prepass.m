function MBSREM_prepass(options)
%% Count the number of pixels each LOR traverses
% This function counts the number of pixels that each LOR traverses. This
% is needed for methods 2, 3 and 4 or for method 1 when precompute_lor has
% been selected.

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
diameter = options.diameter;
FOVax = options.FOVa_x;
FOVay = options.FOVa_y;
axial_fov = options.axial_fov;
Nx = options.Nx;
Ny = options.Ny;
Nz = options.Nz;
machine_name = options.machine_name;
pseudot = int32(options.pseudot);
NSinos = options.NSinos;
TotSinos = options.TotSinos;
det_per_ring = options.det_per_ring;
% Nang = options.Nang;
% Ndist = options.Ndist;

blocks=int32(rings + sum(pseudot));
block1=int32(1);

temp = pseudot;
if ~isempty(temp) && temp > 0
    for kk = int32(1) : temp
        pseudot(kk) = int32(options.cryst_per_block + 1) * kk;
    end
end

% PET-laitteen halkaisija (mm)
R=double(diameter);
% Transaxial FOV (mm)
FOVax=double(FOVax);
FOVay=double(FOVay);
% Axial FOV (mm)
axial_fov = double(axial_fov);

NSinos = int32(NSinos);
NSlices = int32(Nz);
TotSinos = int32(TotSinos);

% Voxel counts
Ny=int32(Ny);
Nx=int32(Nx);
Nz=int32(Nz);

%%

if options.use_raw_data || options.precompute_all
    
    load([machine_name '_detector_coordinates.mat'],'x','y');
    LL = form_detector_pairs_raw(rings, det_per_ring);
    
    z_length = double(blocks) * options.cr_pz;
    zl = z_length - axial_fov;
    z = linspace(-zl/2, z_length-zl/2, blocks + 1);
    z = z(2:end) - options.cr_pz/2;
    
    if options.attenuation_correction
        data = load(options.attenuation_datafile);
        variables = fields(data);
        vaimennus = double(data.(variables{1}));
        if size(vaimennus,1) ~= Nx || size(vaimennus,2) ~= Ny || size(vaimennus,3) ~= Nz
            if size(vaimennus,1) ~= Nx*Ny*Nz
                error("Error: Attenuation data is of different size than the reconstructed image")
            end
        end
        if size(vaimennus,2) == 1
            vaimennus = vaimennus(:,:,2*block1+1:2*blocks+1);
        else
            vaimennus = vaimennus(2*block1-1:(2*blocks-1)*Nx*Ny);
        end
        vaimennus = vaimennus(:);
        clear data
    else
        vaimennus = 0;
    end
    
%     if min(min(z)) == 0
%         z = z + (axial_fov - max(max(z)))/2;
%     end
    
    x=double(x);
    y=double(y);
    z_det = double(z);
    clear z
    
%     if isempty(pseudot)
%         pseudot = int32(1e5);
%     else
%         pseudot = pseudot - 1;
%     end
    if sum(pseudot) == 0
        pseudot = [];
    end
    
    % for the precomputed version, index vectors are needed
    
    LL = LL';
    LL = LL(:);
    
    
    
    % Pikselit
    etaisyys=(R-FOVax)/2;
    xx=linspace(etaisyys,R-etaisyys,Nx+1);
    etaisyys=(R-FOVay)/2;
    yy=linspace(etaisyys,R-etaisyys,Ny+1);
    zz=linspace(double(0),double(axial_fov),Nz+1);
    zz=zz(2*block1-1:2*blocks);
    
    % Pikselien etäisyys toisistaan
    dx=diff(xx(1:2));
    dy=diff(yy(1:2));
    dz=diff(zz(1:2));
    
    % Kuvan etäisyys origosta
    bx=xx(1);
    by=yy(1);
    bz=zz(1);
    
    size_x = int32(size(x,1));
    
    zmax = max(max(z_det));
    if zmax==0
        zmax = double(1);
    end
    
    % Determine which LORs go through the FOV
    
    [discard, lor, D, Amin] = MBSREM_precompute( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, ...
        size_x, zmax, vaimennus, options.attenuation_correction, NSinos, LL, pseudot, int32(det_per_ring), ...
        options.verbose, true, block1, blocks);
    Amin = Amin(discard);
    lor = lor(discard);
    clear LL
    
    save([machine_name '_MBSREM_precompute_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'],'D','Amin','-v7.3')
    save([machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'],'discard','lor','-v7.3')
    clear discard Amin lor D
end
%%

if ~options.use_raw_data || options.precompute_all
    
    load([machine_name '_app_coordinates_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'x','y');
    load([machine_name '_3D_coordinates_span' num2str(options.span) '_ringdiff' num2str(options.ring_difference) '_' num2str(TotSinos) '.mat'],'z')
    
    
    if min(min(z)) == 0
        z = z + (axial_fov - max(max(z)))/2;
    end
    
    x=double(x);
    y=double(y);
    z_det = double(z);
    
    
    
    % Pikselit
    etaisyys=(R-FOVax)/2;
    xx=linspace(etaisyys,R-etaisyys,Nx+1);
    etaisyys=(R-FOVay)/2;
    yy=linspace(etaisyys,R-etaisyys,Ny+1);
    zz=linspace(double(0),double(axial_fov),Nz+1);
    zz=zz(2*block1-1:2*blocks);
    
    % Pikselien etäisyys toisistaan
    d=diff(xx(1:2));
    dy=diff(yy(1:2));
    dz=diff(zz(1:2));
    
    % Kuvan etäisyys origosta
    bx=xx(1);
    by=yy(1);
    bz=zz(1);
    
    size_x = int32(size(x,1));
    
    zmax = max(max(z_det));
    if zmax==0
        zmax = double(1);
    end
    
    if options.attenuation_correction
        data = load(options.attenuation_datafile);
        variables = fields(data);
        vaimennus = double(data.(variables{1}));
        if size(vaimennus,1) ~= Nx || size(vaimennus,2) ~= Ny || size(vaimennus,3) ~= Nz
            if size(vaimennus,1) ~= Nx*Ny*Nz
                error("Error: Attenuation data is of different size than the reconstructed image")
            end
        end
        if size(vaimennus,2) == 1
            vaimennus = vaimennus(:,:,2*block1-1:2*blocks-1);
        else
            vaimennus = vaimennus(2*block1-1:(2*blocks-1)*Nx*Ny);
        end
        vaimennus = vaimennus(:);
        clear data
    else
        vaimennus = 0;
    end
    
    % Determine which LORs go through the FOV
    
    [discard, lor, D, Amin] = MBSREM_precompute( Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, ...
        size_x, zmax, vaimennus, options.attenuation_correction, NSinos, uint16(0), pseudot, int32(det_per_ring), ...
        options.verbose, false, block1, blocks);
    Amin = single(Amin);
%     Amin = Amin(discard);
%     lor = lor(discard);
    
    save([machine_name '_MBSREM_precompute_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_sino_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'D','Amin','-v7.3')
    save([machine_name '_lor_pixel_count_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_sino_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'lor','discard','-v7.3')
end

end