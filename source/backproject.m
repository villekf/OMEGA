function [bp, varargout] = backproject(options, index, n_meas, rhs, nn, iternn, SinM)
%BACKPROJECT Calculates the backprojection
% Example:
%   bp = backproject(options, index, n_meas, rhs, nn)
%   [bp, norm] = backproject(options, index, n_meas, rhs, nn)
% INPUTS:
%   options = Contains the machine and (optionally) sinogram specific
%   information (e.g. detectors per ring, image (matrix) size, FOV size)
%   index = The indices (LORs) used to compute the system matrix (you can
%   use index_maker to produce the indices)
%   n_meas = Number of measurements used
%   rhs = The right hand side of the input (i.e. bp = A' * rhs)
%   nn = The interval from where the measurements are taken (current
%   subset)
%
% OUTPUTS:
%   bp = The backprojection (bp = A' * rhs)
%   norm = The (optional) normalization constant (norm = sum(A,1))
%
% See also index_maker, forward_project

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ville-Veikko Wettenhovi
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


if nargout == 1
    no_norm = true;
else
    no_norm = false;
end

if nargout > 2
    error('Too many output arguments')
end
if nargout == 0
    error('Too few output arguments')
end

% folder = fileparts(which('reconstructions_main.m'));
% folder = strrep(folder, 'source','mat-files/');
% folder = strrep(folder, '\','/');

Nx = options.Nx;
Ny = options.Ny;
Nz = options.Nz;
NSlices = uint32(Nz);
attenuation_correction = options.attenuation_correction;
FOVax = options.FOVa_x;
FOVay = options.FOVa_y;
rings = options.rings;
if options.use_raw_data && isfield(options,'x')
    det_per_ring = numel(options.x);
else
    det_per_ring = options.det_per_ring;
end
% machine_name = options.machine_name;
% Nang = options.Nang;
% Ndist = options.Ndist;
% cr_pz = options.cr_pz;
use_raw_data = options.use_raw_data;
% attenuation_datafile = options.attenuation_datafile;
pseudot = uint32(options.pseudot);
temp = pseudot;
if ~isempty(temp) && temp > 0
    for kk = int32(1) : temp
        pseudot(kk) = int32(options.cryst_per_block + 1) * kk;
    end
elseif temp == 0
    pseudot = [];
end
% Diameter of the PET-device (bore) (mm)
R=double(options.diameter);
% Transaxial FOV (x-direction, horizontal) (mm)
FOVax=double(FOVax);
% Transaxial FOV (y-direction, vertical) (mm)
FOVay=double(FOVay);
% Axial FOV (mm)
axial_fow = double(options.axial_fov);
% Number of rings
blocks=uint32(rings + length(pseudot) - 1);
% First ring
block1=uint32(0);

NSinos = uint32(options.NSinos);
% TotSinos = int32(options.TotSinos);

if iternn == 1 || (options.implementation > 1 && (options.n_rays_transaxial > 1 || options.n_rays_axial > 1) && ~options.precompute_lor && options.projector_type == 1)
    if (options.n_rays_transaxial > 1 || options.n_rays_axial > 1) && isfield(options,'x') && isfield(options,'y')
        options = rmfield(options, 'x');
        options = rmfield(options, 'y');
    end
    [x, y, z_det, options] = get_coordinates(options, blocks, pseudot);
    options.x = x;
    options.y = y;
    options.z_det = z_det;
else
    x = options.x;
    y = options.y;
    z_det = options.z_det;
end

[normalization_correction, randoms_correction, options] = set_up_corrections(options, blocks);

if randoms_correction
    if iscell(options.SinDelayed)
        if options.implementation == 1
            SinDelayed = options.SinDelayed{1}(nn(1) + 1 : nn(2));
        else
            SinDelayed{1} = options.SinDelayed{1}(nn(1) + 1 : nn(2));
        end
    else
        if options.implementation == 1
            SinDelayed = options.SinDelayed(nn(1) + 1 : nn(2));
        else
            SinDelayed{1} = options.SinDelayed(nn(1) + 1 : nn(2));
        end
    end
else
    if options.implementation == 1
        SinDelayed = 0;
    else
        SinDelayed{1} = single(0);
    end
end

if normalization_correction
    normalization = options.normalization(nn(1) + 1 : nn(2));
else
    if options.implementation == 1
        normalization = 0;
    else
        normalization = single(0);
    end
end

if options.scatter_correction && ~options.subtract_scatter
    if options.implementation == 1
        scatter_input = options.ScatterC(nn(1) : nn(2));
    else
        if iscell(options.SinDelayed)
            options.ScatterFB{1} = {single(options.ScatterC{1}(nn(1) : nn(2)))};
        else
            options.ScatterFB{1} = {single(options.ScatterC(nn(1) : nn(2)))};
        end
    end
else
    if options.implementation == 1
        scatter_input = 0;
    else
        options.ScatterFB{1} = single(0);
    end
end

if options.use_raw_data
    size_x = uint32(options.det_per_ring);
else
    size_x = uint32(options.Nang*options.Ndist);
end

if (options.precompute_lor  || options.implementation == 5 || options.implementation == 2 || options.implementation == 3)
    n_meas = [0;cumsum(n_meas)];
    if iscell(index)
        index = cell2mat(index);
    end
end

if use_raw_data
    if isempty(pseudot)
        pseudot = int32(1e5);
    else
        pseudot = pseudot - 1;
    end
end



[options, lor_a, xy_index, z_index, LL, summa, n_meas] = form_subset_indices(options, n_meas, 1, index, size_x, y, z_det, blocks, true);
if ~options.precompute_lor
    lor_a = uint16(0);
end


% Pixels
etaisyys_x=(R-FOVax)/2;
etaisyys_y=(R-FOVay)/2;
if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
    zz=linspace(single(0),single(axial_fow),Nz+1);
    xx = single(linspace(etaisyys_x,R-etaisyys_x,Nx+1));
    yy = single(linspace(etaisyys_y,R-etaisyys_y,Ny+1));
else
    zz=linspace(double(0),double(axial_fow),Nz+1);
    xx = double(linspace(etaisyys_x,R-etaisyys_x,Nx+1));
    yy = double(linspace(etaisyys_y,R-etaisyys_y,Ny+1));
end
zz=zz(2*block1+1:2*blocks);

% Distance of adjacent pixels
dx=diff(xx(1:2));
dy=diff(yy(1:2));
dz=diff(zz(1:2));

% Distance of image from the origin
bx=xx(1);
by=yy(1);
bz=zz(1);

% Number of pixels
Ny=uint32(Ny);
Nx=uint32(Nx);
Nz=uint32(Nz);

N=(Nx)*(Ny)*(Nz);
det_per_ring = uint32(det_per_ring);

% How much memory is preallocated
if use_raw_data == false
    ind_size = uint32(NSinos/8*(det_per_ring)* Nx * (Ny));
else
    ind_size = uint32((det_per_ring)^2/8* Nx * (Ny));
end

if ~options.precompute_lor
    lor_a = uint16(0);
end

zmax = max(max(z_det));
if zmax==0
    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
        zmax = single(1);
    else
        zmax = double(1);
    end
end

if options.projector_type == 2 || options.projector_type == 3
    x_center = xx(1 : end - 1)' + dx/2;
    y_center = yy(1 : end - 1)' + dy/2;
%     if options.tube_width_z > 0
        z_center = zz(1 : end - 1)' + dz/2;
%     else
%         z_center = zz(1);
%     end
    temppi = min([options.FOVa_x / options.Nx, options.axial_fov / options.Nz]);
    if options.tube_width_z > 0
        temppi = max([1,round(options.tube_width_z / temppi)]);
    else
        temppi = max([1,round(options.tube_width_xy / temppi)]);
    end
    temppi = temppi * temppi * 4;
    if options.apply_acceleration
        if options.tube_width_z == 0
            dec = uint32(sqrt(options.Nx^2 + options.Ny^2) * temppi);
        else
            dec = uint32(sqrt(options.Nx^2 + options.Ny^2 + options.Nz^2) * temppi);
        end
    else
        dec = uint32(0);
    end
else
    x_center = xx(1);
    y_center = yy(1);
    z_center = zz(1);
    if options.use_psf && options.apply_acceleration
        dec = uint32(sqrt(options.Nx^2 + options.Ny^2 + options.Nz^2) * ((ceil(options.cr_pz / dx) * 2 + 1)^2));
    else
        dec = uint32(0);
    end
end
if options.projector_type == 3
    voxel_radius = (sqrt(2) * options.voxel_radius * dx) / 2;
    bmax = options.tube_radius + voxel_radius;
    b = linspace(0, bmax, 10000)';
    b(options.tube_radius > (b + voxel_radius)) = [];
    b = unique(round(b*10^3)/10^3);
    V = volumeIntersection(options.tube_radius, voxel_radius, b);
    Vmax = (4*pi)/3*voxel_radius^3;
    bmin = min(b);
else
    V = 0;
    Vmax = 0;
    bmin = 0;
    bmax = 0;
end
if options.implementation == 2 || options.implementation == 3
    V = single(V);
    Vmax = single(Vmax);
    bmin = single(bmin);
    bmax = single(bmax);
end
% Multi-ray Siddon
if options.implementation > 1 && options.n_rays_transaxial > 1 && ~options.precompute_lor && options.projector_type == 1
    [x,y] = getMultirayCoordinates(options);
    options.x = x;
    options.y = y;
end

if options.implementation == 1
    if options.precompute_lor == false
        iij = double(0:Nx);
        jji = double(0:Ny);
        kkj = double(0:Nz);
        if use_raw_data == false
            if options.projector_type == 1 || options.projector_type == 0
                [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                    zmax, options.vaimennus, normalization, SinDelayed, n_meas(end), attenuation_correction, normalization_correction, ...
                    randoms_correction, options.scatter, scatter_input, options.global_correction_factor, uint16(0), uint32(0), uint32(0), NSinos, uint16(0), pseudot, det_per_ring, options.verbose, ...
                    use_raw_data, uint32(2), ind_size, block1, blocks, index, ...
                    uint32(options.projector_type), iij, jji, kkj);
            else
                error('Unsupported projector type')
            end
        else
%             L = LL(index,:);
            LL = LL';
            LL = LL(:);
            if options.projector_type == 1
                [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                    zmax, options.vaimennus, normalization, SinDelayed, uint32(0), attenuation_correction, normalization_correction, ...
                    randoms_correction, options.scatter, scatter_input, options.global_correction_factor, uint16(0), uint32(0), uint32(0), NSinos, LL, pseudot, det_per_ring, options.verbose, ...
                    use_raw_data, uint32(2), ind_size, block1, blocks, uint32(0), uint32(options.projector_type), iij, jji, kkj);
            else
                error('Unsupported projector type')
            end
        end
        if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
            lor = repeat_elem(uint32(1:length(lor))',uint32(lor));
        elseif exist('OCTAVE_VERSION','builtin') == 5
            lor = repelem(uint32(1:length(lor)),uint32(lor));
        else
            lor = repelem(uint32(1:length(lor)),uint32(lor))';
        end
        
        A_length = length(SinM);
        if options.verbose
            tStart = tic;
        end
        if exist('OCTAVE_VERSION','builtin') == 0 && ~verLessThan('matlab','9.8')
            indices = uint32(indices) + 1;
            A = sparse(lor,indices,double(alkiot), A_length, double(N));
        elseif options.use_fsparse && exist('fsparse','file') == 3
            indices = int32(indices) + 1;
            A = fsparse(int32(lor),(indices),double(alkiot),[A_length double(N) length(alkiot)]);
        elseif options.use_fsparse && exist('fsparse','file') == 0
            warning('options.fsparse set to true, but no FSparse mex-file found. Using regular sparse')
            indices = double(indices) + 1;
            A = sparse(double(lor),(indices),double(alkiot), A_length, double(N));
        else
            indices = double(indices) + 1;
            A = sparse(double(lor),(indices),double(alkiot), A_length, double(N));
        end
        clear indices alkiot lor
        if options.verbose
            tElapsed = toc(tStart);
            disp(['Sparse matrix formation took ' num2str(tElapsed) ' seconds'])
        end
    else
        if use_raw_data
            xy_index = uint32(0);
            z_index = uint16(0);
        else
            LL = uint16(0);
        end
        if options.projector_type == 2
            lor2 = [0; cumsum(uint64(lor_orth(n_meas(1)+1:n_meas(2))))];
        else
            lor2 = [0; cumsum(uint64(lor_a(n_meas(1)+1:n_meas(2))))];
        end
        [A, ~] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
            normalization, SinDelayed, n_meas(end), attenuation_correction, normalization_correction,...
            randoms_correction, options.scatter, scatter_input, options.global_correction_factor, lor_a, xy_index, z_index, NSinos, ...
            LL, pseudot, det_per_ring, options.verbose, use_raw_data, uint32(0), lor2, summa, false, ...
            uint32(options.projector_type), options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, int32(0), bmin, bmax, Vmax, V);
        clear lor2
    end
    
    if size(A,2) == size(rhs,1)
        if no_norm
            bp = A * rhs;
        else
            bp = A * rhs;
            varargout{1} = full(sum(A,2));
        end
    else
        if no_norm
            bp = A' * rhs;
        else
            bp = A' * rhs;
            varargout{1} = full(sum(A,1))';
        end
    end
    
else
%     options = double_to_single(options);
    if use_raw_data
        xy_index = uint32(0);
        z_index = uint32(0);
    else
        if isempty(pseudot)
            pseudot = int32(100000);
        end
        LL = uint16(0);
    end
    if ~iscell(SinM)
        SinM = {single(SinM)};
    end
    tube_width_xy = single(options.tube_width_xy);
    crystal_size_z = single(options.tube_width_z);
    n_rays = uint16(options.n_rays_transaxial);
    n_rays3D = uint16(options.n_rays_axial);
    dc_z = single(z_det(2,1) - z_det(1,1));
    if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
        randoms = uint32(1);
    else
        randoms = uint32(0);
    end
    if (options.projector_type == 1 && (options.precompute_lor || (n_rays + n_rays3D) <= 2)) || options.projector_type == 2 || options.projector_type == 3
        kernel_file = 'multidevice_kernel.cl';
        kernel_path = which(kernel_file);
        kernel_path = strrep(kernel_path, '\', '/');
        kernel_path = strrep(kernel_path, '.cl', '');
        filename = 'OMEGA_matrix_free_OpenCL_binary_device';
        header_directory = strrep(kernel_path,'multidevice_kernel','');
    elseif options.projector_type == 1 && ~options.precompute_lor
        kernel_file = 'multidevice_siddon_no_precomp.cl';
        kernel_path = which(kernel_file);
        kernel_path = strrep(kernel_path, '\', '/');
        kernel_path = strrep(kernel_path, '.cl', '');
        filename = 'OMEGA_matrix_free_OpenCL_binary_device';
        header_directory = strrep(kernel_path,'multidevice_siddon_no_precomp','');
    else
        error('Invalid projector for OpenCL')
    end

    filename = [header_directory, filename];
    header_directory = strcat('-I "', header_directory);
    header_directory = strcat(header_directory,'"');
    tic
    [output] = OpenCL_matrixfree_multi_gpu( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end), ...
        single(NSlices), size_x, zmax, options.verbose, LL, pseudot, det_per_ring, uint32(options.use_device), filename, uint8(use_raw_data), ...
        single(options.cpu_to_gpu_factor), uint32(0), header_directory, options.vaimennus, normalization, n_meas(end), uint32(attenuation_correction), ...
        uint32(normalization_correction), lor_a, xy_index, z_index, tube_width_xy, crystal_size_z, x_center, y_center, z_center, SinDelayed, randoms, ...
        uint32(options.projector_type), options.precompute_lor, int32(dec), n_rays, n_rays3D, dc_z, SinM, logical(options.use_64bit_atomics), rhs, no_norm, ...
        options.global_correction_factor, bmin, bmax, Vmax, V, options.use_psf, options);
    toc
    if isa(output{1},'uint64')
        bp = single(output{1}) / 100000000000;
    else
        bp = output{1};
    end
    if isa(output{2},'uint64')
        varargout{1} = single(output{2}) / 100000000000;
    else
        varargout{1} = output{2};
    end
end
if options.verbose
    disp('Backprojection done')
end