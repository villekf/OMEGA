function [A] = observation_matrix_formation_nongate(options, current_subset, index, LL, pituus, lor_a)
%% PRECOMPUTED SYSTEM MATRIX
% Precomputes the system matrix for PET reconstruction to be used in
% equations of form y = Ax.
%
% Input arguments are options variables from main_nongate and the current
% subset (use current_subset = 1 and options.subsets = 1 to compute the
% entire system matrix)
%
% Required options inputs are: Nx, Ny, Nz, attenuation_correction, diameter,
% FOVa, axial_fov, NSinos, pseudot, rings, det_per_ring, machine_name,
% TotSinos, attenuation_datafile, verbose, precompute_lor, Ndist, Nang,
% subsets, span, ring difference
%
% Outputs the transposed system matrix for PET reconstruction.
%
% Requires prepass phase.

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


use_openCL = false;

Nx = options.Nx;
Ny = options.Ny;
Nz = options.Nz;
attenuation_correction = options.attenuation_correction;
diameter = options.diameter;
FOVa = options.FOVa;
axial_fov = options.axial_fov;
NSinos = options.NSinos;
pseudot = int32(options.pseudot);
rings = options.rings;
det_per_ring = options.det_per_ring;
machine_name = options.machine_name;
TotSinos = options.TotSinos;
attenuation_datafile = options.attenuation_datafile;
subsets = options.subsets;


% Diameter of the PET-device (bore) (mm)
R=double(diameter);
% Transaxial FOV (mm)
FOVa=double(FOVa);
% Axial FOV (mm)
axial_fow = double(axial_fov);
% Number of rings
blocks=int32(rings + length(pseudot) - 1);
% First ring
block1=int32(0);
% Pixel count in x- and y-directions
pikselikoko=int32(Nx);

NSinos = int32(NSinos);
NSlices = int32(Nz);
TotSinos = int32(TotSinos);

load([machine_name '_app_coordinates_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'x','y');
load([machine_name '_3D_coordinates_span' num2str(options.span) '_ringdiff' num2str(options.ring_difference) '_' num2str(TotSinos) '.mat'],'z')
if NSinos ~= TotSinos
    z = z(1:NSinos,:);
end

if attenuation_correction
    data = load(attenuation_datafile);
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
        vaimennus = vaimennus(2*block1+1:(2*blocks+1)*Nx*Ny);
    end
    vaimennus = vaimennus(:);
    clear data
else
    vaimennus = 0;
end

if min(min(z)) == 0
    z = z + (axial_fow - max(max(z)))/2;
end

x=double(x);
y=double(y);
z_det = double(z);

size_x = int32(size(x,1));

if options.precompute_lor && subsets > 1 || options.reconstruction_method == 3  && subsets > 1 || options.reconstruction_method == 2 && subsets > 1
    pituus = [0;cumsum(pituus)];
    if iscell(index)
        index = cell2mat(index);
    end
end

if options.use_raw_data
    if isempty(pseudot)
        pseudot = int32(1e5);
    else
        pseudot = pseudot - 1;
    end
end


% for the precomputed version, index vectors are needed
if options.use_raw_data == false && options.precompute_lor || options.use_raw_data == false && options.reconstruction_method == 3
    
    if subsets > 1
        load([machine_name '_lor_pixel_count_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_sino_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'lor')
        lor_a = (lor(index));
        clear lor
    else
        load([machine_name '_lor_pixel_count_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_sino_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'lor','discard')
        if length(discard) ~= options.TotSinos*options.Nang*options.Ndist
            error('Error: Size mismatch between sinogram and LORs to be removed')
        end
        if options.use_raw_data == false && options.NSinos ~= options.TotSinos
            discard = discard(1:options.NSinos*options.Nang*options.Ndist);
        end
        lor_a = (lor(discard));
        clear lor
    end
    [~, I] = sort(y, 2);
    sy = size(y);
    I = sub2ind(sy, repmat((1:sy(1)).', 1, sy(2)), I);
    
    xy_index = uint32(I(:,1));
    xy_index2 = repmat(uint32(1:size_x)', NSinos - Nz, 1);
    xy_index = [repmat(xy_index, Nz, 1); xy_index2];
    [~, I] = sort(z_det, 2);
    sy = size(z_det);
    I = sub2ind(sy, repmat((1:sy(1)).', 1, sy(2)), I);
    
    z_index = uint32(I(:,1));
    z_index = repelem(z_index, size_x);
    if subsets > 1
        z_index = z_index(index);
    else
        z_index = (z_index(discard));
    end
    apu = z_index > NSinos;
    z_index = z_index - 1;
    
    if subsets > 1
        xy_index = xy_index(index);
    else
        xy_index = (xy_index(discard));
    end
    xy_index(apu) = xy_index(apu) + uint32(size_x);
    xy_index = xy_index - 1;
    
    summa = zeros(subsets, 1, 'uint64');
    
    if subsets > 1
        for kk = 1 : subsets
            summa(kk) = sum(int64(lor_a(pituus(kk)+1:pituus(kk+1))));
%             koko(kk) = -(pituus(kk))+pituus(kk+1);
        end
    end
    
    clear discard I yt xt xy_index2 index apu
elseif options.use_raw_data && options.precompute_lor || options.use_raw_data && options.reconstruction_method == 3  || options.use_raw_data && options.reconstruction_method == 2
    
    LL = LL(index,:);
    lor_a = (lor_a(index));
    summa = zeros(subsets, 1, 'uint64');
%     y = circshift(y, -length(y)/4);
%     x = circshift(x, -length(x)/4);
    
    for kk = 1 : subsets
        apu = LL(pituus(kk) + 1 : pituus(kk + 1),:) - 1;
        apu2 = idivide(apu, uint16(det_per_ring));
        idx = apu2(:,1) == apu2(:,2);
        apu2 = apu(idx,:);
        ind = mod(apu2, uint16(det_per_ring)) + 1;
        yt = y(ind);
        y_i = yt(:,1) > yt(:,2);
        apu2(y_i,:) = fliplr(apu2(y_i,:));
        apu(idx,:) = apu2;
        LL(pituus(kk) + 1 : pituus(kk + 1),:) = apu + 1;
        summa(kk) = sum(int64(lor_a(pituus(kk)+1:pituus(kk+1))));
%         koko(kk) = -(pituus(kk))+pituus(kk+1);
    end
    
    clear apu apu2 idx ind yt y_i index discard
    
    LL = LL';
    LL = LL(:);
end


% Pikselit
etaisyys=(R-FOVa)/2;
zz=linspace(double(0),double(axial_fow),Nz+1);
xx = double(linspace(etaisyys,R-etaisyys,pikselikoko+1));
yy = double(linspace(etaisyys,R-etaisyys,pikselikoko+1));
zz=zz(2*block1+1:2*blocks+2);

% Pikselien et?isyys toisistaan
dx=diff(xx(1:2));
dy=diff(yy(1:2));
dz=diff(zz(1:2));

% Kuvan et?isyys origosta
bx=xx(1);
by=yy(1);
bz=zz(1);

% Pikselien m??r?
Ny=int32(Ny);
Nx=int32(Nx);
Nz=int32(Nz);


% iij=double(0:Nx);
% jji=double(0:Ny);
% kkj=double(0:Nz);

N=(Nx)*(Ny)*(Nz);
det_per_ring = int32(det_per_ring);

% Kuinka monelle alkiolle tallennetaan muistia etuk?teen
ind_size = int32(NSinos/subsets*(det_per_ring)* Nx * (Ny));


zmax = max(max(z_det));
if zmax==0
    zmax = double(1);
end

%%

if use_openCL
    kernel_path = which('system_matrix_kernel.cl');
    lor2 = [0; cumsum(uint32(lor_a(pituus(current_subset)+1:pituus(current_subset + 1))))];
    if options.use_raw_data == false
        [A] = system_matrix_openCL( kernel_path, Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , NSinos, double(NSlices), size_x, zmax, NSinos, ...
            vaimennus, (pituus(current_subset + 1) - pituus(current_subset)), int32(attenuation_correction), lor_a(pituus(current_subset)+1:pituus(current_subset + 1)), ...
            summa(current_subset), xy_index(pituus(current_subset)+1:pituus(current_subset + 1)), z_index(pituus(current_subset)+1:pituus(current_subset + 1)), lor2, ...
            options.verbose, int32(options.use_device));
    else
        
    end
else
    if options.precompute_lor == false
        if options.use_raw_data == false
%             [ lor, indices, alkiot] = improved_Siddon_algorithm( pikselikoko, Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , NSinos, NSlices, size_x, ...
%                 zmax, NSinos, ind_size, vaimennus, index{current_subset}, pituus(current_subset), attenuation_correction, options.verbose);
            [ lor, indices, alkiot] = improved_Siddon_algorithm( options.verbose, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx ,...
                NSinos, NSlices, size_x, zmax, NSinos, ind_size, vaimennus, index{current_subset}, pituus(current_subset), attenuation_correction, use_raw_data,...
                uint16(0), pseudot, block1, blocks, det_per_ring);
        else
            L = LL(index{current_subset},:);
            L = L';
            L = L(:);
%             [ lor, indices, alkiot] = improved_Siddon_algorithm_raw( pikselikoko, Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , NSinos, NSlices, ...
%                 size_x, zmax, NSinos, ind_size, vaimennus, L, pseudot, block1, blocks, det_per_ring, attenuation_correction, options.verbose);
            [ lor, indices, alkiot] = improved_Siddon_algorithm( options.verbose, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx ,...
                NSinos, NSlices, size_x, zmax, NSinos, ind_size, vaimennus, uint32(0), int32(0), attenuation_correction, use_raw_data,...
                L, pseudot, block1, blocks, det_per_ring);
        end
        lor = reshape(lor,[],2);
        lor=repelem(int32((lor(:,1))),lor(:,2));
        
        A_length = pituus(current_subset+1) - pituus(current_subset);
        indices=indices + 1;
        if verbose
            tStart = tic;
        end
        if options.use_fsparse == false
            A = sparse(double(lor),double(indices),double(alkiot), A_length, double(N));
        else
            A = fsparse(lor,indices,double(alkiot),[A_length double(N) length(alkiot)]);
        end
        A = A';
        clear indices alkiot lor
        if verbose
            tElapsed = toc(tStart);
            disp(['Sparse matrix formation took ' num2str(tElapsed) ' seconds'])
        end
    else
        lor2 = [0; cumsum(uint32(lor_a(pituus(current_subset)+1:pituus(current_subset + 1))))];
        if options.use_raw_data == false
%             [A] = improved_Siddon_algorithm_array( Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , NSinos, NSlices, size_x, zmax, vaimennus, lor2, ...
%                 pituus(current_subset + 1) - pituus(current_subset), attenuation_correction, lor_a(pituus(current_subset)+1:pituus(current_subset + 1)), ...
%                 summa(current_subset), xy_index(pituus(current_subset)+1:pituus(current_subset + 1)), z_index(pituus(current_subset)+1:pituus(current_subset + 1)), NSinos, ...
%                 options.verbose);
            [A] = improved_Siddon_algorithm_array( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, vaimennus, ...
                lor2, pituus(current_subset + 1) - pituus(current_subset), attenuation_correction, lor_a(pituus(current_subset)+1:pituus(current_subset + 1)), ...
                summa(current_subset), xy_index(pituus(current_subset)+1:pituus(current_subset + 1)), z_index(pituus(current_subset)+1:pituus(current_subset + 1)), ...
                NSinos, uint16(0), pseudot, det_per_ring, options.verbose, use_raw_data);
        else
%             [A] = improved_Siddon_algorithm_array_raw( Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , size_x, zmax, vaimennus, lor2, ...
%                 pituus(current_subset + 1) - pituus(current_subset), attenuation_correction, lor_a(pituus(current_subset)+1:pituus(current_subset + 1)), ...
%                 summa(current_subset), LL(pituus(current_subset) * 2 + 1 : pituus(current_subset + 1) * 2), pseudot, det_per_ring, options.verbose);
            [A] = improved_Siddon_algorithm_array( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, vaimennus, ...
                lor2, pituus(current_subset + 1) - pituus(current_subset), attenuation_correction, lor_a(pituus(current_subset)+1:pituus(current_subset + 1)), ...
                summa(current_subset), uint32(0), uint32(0), NSinos, LL(pituus(current_subset) * 2 + 1 : pituus(current_subset + 1) * 2), pseudot, ...
                det_per_ring, options.verbose, use_raw_data);
        end
        clear lor2
    end
end
end