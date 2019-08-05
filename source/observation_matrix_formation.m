function A = observation_matrix_formation(diameter, FOVa, axial_fov, rings, pseudot, Nx, Ny, Nz, det_per_ring, cr_pz, use_fsparse, attenuation_correction, attenuation_datafile, ...
    precompute_lor, use_raw_data, pituus, options, NSlices)
%% PRECOMPUTED SYSTEM MATRIX
% Precomputes the system matrix for PET reconstruction to be used in
% equations of form y = Ax.
%
% Input arguments are machine diameter, transaxial FOV, axial FOV, number
% of rings, possible pseudo rings, image size for X, Y and Z directions,
% number of detectors per ring, crystal pitch in z-direction, boolean value
% on whether fsparse is used, boolean value on whether attenuation
% correction is applied, filename for attenuation images, boolean value on
% whether precomputed LORs are used, boolean value on whether raw data is
% used, the total number of measurements.
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

% Diameter of the PET device/bore (mm)
R=double(diameter);
% Transaxial FOV (mm)
FOVa=double(FOVa);
% Axial FOV (mm)
axial_fow = double(axial_fov);
% Number of rings
blocks=int32(rings + length(pseudot) - 1);
% First ring
block1=int32(0);
% Number of pixels in x- and y-directions
pikselikoko=int32(Nx);

machine_name = options.machine_name;

if use_raw_data == false
    if precompute_lor
        load([machine_name '_lor_pixel_count_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_sino_' num2str(options.Ndist) 'x' num2str(options.Nang)  'x' num2str(options.TotSinos) '.mat'],'lor')
        if options.NSinos ~= options.TotSinos
            lor = lor(1:options.NSinos*options.Nang*options.Ndist);
        end
        summa = uint64(sum(lor));
        discard = lor > 0;
        lor_a = lor(discard);
        clear lor
    end
    load([machine_name '_app_coordinates_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'x','y');
    load([machine_name '_3D_coordinates_span' num2str(options.span) '_ringdiff' num2str(options.ring_difference) '_' num2str(options.TotSinos) '.mat'],'z')
    if options.NSinos ~= options.TotSinos
        z = z(1:options.NSinos,:);
    end
else
    load([machine_name '_detector_coordinates.mat'],'x','y');
    if precompute_lor
        load([machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'],'LL','lor')
        summa = uint64(sum(lor));
        discard = lor > 0;
        lor_a = lor(discard);
        clear lor
    else
        load([machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'],'LL')
    end
    z_length = double(blocks) * cr_pz;
    z = linspace(0, z_length, blocks + 1);
    
    
end


x=double(x);
y=double(y);
z_det = double(z);

size_x = int32(size(x,1));

if use_raw_data == false && options.precompute_lor || use_raw_data == false && options.reconstruction_method == 3
    
    [~, I] = sort(y, 2);
    sy = size(y);
    I = sub2ind(sy, repmat((1:sy(1)).', 1, sy(2)), I);
    
    xy_index = uint32(I(:,1));
    xy_index2 = repmat(uint32(1:size_x)', options.NSinos - Nz, 1);
    xy_index = [repmat(xy_index, Nz, 1); xy_index2];
    [~, I] = sort(z_det, 2);
    sy = size(z_det);
    I = sub2ind(sy, repmat((1:sy(1)).', 1, sy(2)), I);
    
    z_index = uint32(I(:,1));
    z_index = repelem(z_index, size_x);
    apu = z_index > options.NSinos;
    z_index = z_index - 1;
    
    xy_index(apu) = xy_index(apu) + uint32(size_x);
    xy_index = xy_index - 1;
    
    xy_index = xy_index(discard);
    z_index = z_index(discard);
    
    clear discard I yt xt xy_index2 index apu
elseif use_raw_data && options.precompute_lor || use_raw_data && options.reconstruction_method == 3
    
    apu = LL - 1;
    apu2 = idivide(apu, uint16(det_per_ring));
    idx = apu2(:,1) == apu2(:,2);
    apu2 = apu(idx,:);
    ind = mod(apu2, uint16(det_per_ring)) + 1;
    yt = y(ind);
    y_i = yt(:,1) > yt(:,2);
    apu2(y_i,:) = fliplr(apu2(y_i,:));
    apu(idx,:) = apu2;
    LL = apu + 1;
    
%     vector_size = max(koko);
    clear apu apu2 idx ind yt y_i index discard
    LL = LL';
    LL = LL(:);
end

% Attenuation
if attenuation_correction == true
    data = load(attenuation_datafile);
    variables = fields(data);
    vaimennus = data.(variables{1});
    if numel(vaimennus) ~= Nx*Ny*Nz
        error("Error: Attenuation data is of different size than the reconstructed image")
    elseif size(vaimennus,2) > 1 && size(vaimennus,3) > min(1, Nz - 1) && (size(vaimennus,2) ~= Ny || size(vaimennus,1) ~= Nx || size(vaimennus,3) ~= Nz)
        error("Error: Attenuation data is of different size than the reconstructed image")
    end
    vaimennus = vaimennus(:);
    vaimennus = double(vaimennus) ./ 10;
    clear data
else
    vaimennus = 1;
end

%%

% Pixels
etaisyys=(R-FOVa)/2;
xx=linspace(etaisyys,R-etaisyys,pikselikoko+1);
yy=linspace(etaisyys,R-etaisyys,pikselikoko+1);
zz=linspace(double(0),double(axial_fow),Nz+1);
zz=zz(2*block1+1:2*blocks+2);

% Distance between adjecent pixels
dx=diff(xx(1:2));
dy=diff(yy(1:2));
dz=diff(zz(1:2));

% Distance of the image from the origin
bx=xx(1);
by=yy(1);
bz=zz(1);

% Number of pixels
Ny=int32(Ny);
Nx=int32(Nx);
Nz=int32(Nz);

% iij=double(0:Nx);
% jji=double(0:Ny);
% kkj=double(0:Nz);

N=(Nx)*(Ny)*(Nz);
det_per_ring = int32(det_per_ring);

% Preallocation
ind_size = int32((det_per_ring^2)* Nx * (Ny/2));

size_x = int32(size(x,1));

zmax = max(max(z_det));
if zmax==0
    zmax = double(1);
end
if options.precompute_lor == false
    if use_raw_data == false
%         [ lor, indices, alkiot] = improved_Siddon_algorithm( pikselikoko, Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , NSinos, NSlices, size_x, zmax, NSinos, ...
%             ind_size, vaimennus, int32(1), pituus, attenuation_correction, options.verbose);
        [ lor, indices, alkiot] = improved_Siddon_algorithm( options.verbose, int32(Ny), int32(Nx), int32(Nz), dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx ,...
            int32(options.NSinos), int32(NSlices), size_x, zmax, int32(options.NSinos), ind_size, vaimennus, int32(1), pituus, attenuation_correction, use_raw_data,...
            uint16(0), pseudot, block1, blocks, det_per_ring);
    else
%         [ lor, indices, alkiot] = improved_Siddon_algorithm_raw( pikselikoko, Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , NSinos, NSlices, size_x, zmax, ...
%             NSinos, ind_size, vaimennus, LL, pseudot, block1, blocks, det_per_ring, attenuation_correction, options.verbose);
        [ lor, indices, alkiot] = improved_Siddon_algorithm( options.verbose, int32(Ny), int32(Nx), int32(Nz), dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx ,...
            int32(options.NSinos), int32(NSlices), size_x, zmax, int32(options.NSinos), ind_size, vaimennus, uint32(0), int32(0), attenuation_correction, use_raw_data,...
            L, pseudot, block1, blocks, det_per_ring);
    end
    lor = reshape(lor,[],2);
    lor=repelem(int32((lor(:,1))),lor(:,2));
    A_length = pituus;
    indices=indices + 1;
    if verbose
        tStart = tic;
    end
    if use_raw_data == false
        if use_fsparse == false
            A = sparse(double(lor),double(indices),double(alkiot),A_length,double(N));
        else
            A = fsparse(lor,indices,double(alkiot),[A_length double(N) length(alkiot)]);
        end
    else
        if use_fsparse == false
            A = sparse(double(lor),double(indices),double(alkiot),A_length,double(N));
        else
            A = fsparse(lor,indices,double(alkiot),[A_length double(N) length(alkiot)]);
        end
    end
    clear indices alkiot lor
    if verbose
        tElapsed = toc(tStart);
        disp(['Sparse matrix formation took ' num2str(tElapsed) ' seconds'])
    end
    A = A';
else
    lor2 = [0; cumsum(uint32(lor_a))];
    if use_raw_data == false
%         [A] = improved_Siddon_algorithm_array( Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , NSinos, NSlices, size_x, zmax, vaimennus, lor2, ...
%             pituus, attenuation_correction, lor_a, summa, xy_index, z_index, NSinos, options.verbose);
        [A, ~] = improved_Siddon_algorithm_array( int32(Ny), int32(Nx), int32(Nz), dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , int32(options.NSinos), int32(NSlices), size_x, zmax, vaimennus, ...
            lor2, pituus, attenuation_correction, lor_a, summa, xy_index, z_index, int32(options.NSinos), uint16(0), pseudot, det_per_ring, options.verbose, use_raw_data, false);
    else
%         [A] = improved_Siddon_algorithm_array_raw( Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , size_x, zmax, vaimennus, lor2, pituus, ...
%             attenuation_correction, lor_a, summa, LL, pseudot, det_per_ring, options.verbose);
        [A, ~] = improved_Siddon_algorithm_array( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , int32(NSinos), int32(NSlices), size_x, zmax, vaimennus, ...
            lor2, pituus, attenuation_correction, lor_a, summa, uint32(0), uint32(0), int32(options.NSinos), LL, pseudot, det_per_ring, options.verbose, use_raw_data, false);
    end
    clear lor2
end