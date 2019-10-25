function [bp, varargout] = backproject(options, index, n_meas, rhs, nn)
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
det_per_ring = options.det_per_ring;
% machine_name = options.machine_name;
% Nang = options.Nang;
% Ndist = options.Ndist;
% cr_pz = options.cr_pz;
use_raw_data = options.use_raw_data;
% attenuation_datafile = options.attenuation_datafile;
pseudot = int32(options.pseudot);
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
blocks=int32(rings + length(pseudot) - 1);
% First ring
block1=int32(0);

NSinos = int32(options.NSinos);
% TotSinos = int32(options.TotSinos);

[x, y, z] = get_coordinates(options, blocks);

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

if min(min(z)) == 0
    z = z + (axial_fow - max(max(z)))/2;
end

if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
    x=single(x);
    y=single(y);
    z_det = single(z);
else
    x=double(x);
    y=double(y);
    z_det = double(z);
end
clear z


size_x = int32(size(x,1));

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
zz=zz(2*block1+1:2*blocks+2);

% Distance of adjacent pixels
dx=diff(xx(1:2));
dy=diff(yy(1:2));
dz=diff(zz(1:2));

% Distance of image from the origin
bx=xx(1);
by=yy(1);
bz=zz(1);

% Number of pixels
Ny=int32(Ny);
Nx=int32(Nx);
Nz=int32(Nz);

N=(Nx)*(Ny)*(Nz);
det_per_ring = int32(det_per_ring);

% How much memory is preallocated
if use_raw_data == false
    ind_size = int32(NSinos/8*(det_per_ring)* Nx * (Ny));
else
    ind_size = int32((det_per_ring)^2/8* Nx * (Ny));
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

if options.projector_type == 2
    x_center = xx(1 : end - 1)' + dx/2;
    y_center = yy(1 : end - 1)' + dy/2;
    if options.tube_width_z > 0
        z_center = zz(1 : end - 1)' + dz/2;
    else
        z_center = zz(1);
    end
else
    x_center = xx(1);
    y_center = yy(1);
    z_center = zz(1);
end

if options.implementation == 1
    if options.precompute_lor == false
        if use_raw_data == false
            if options.projector_type == 1 || options.projector_type == 0
                [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                    zmax, vaimennus, normalization, SinDelayed, n_meas(end), attenuation_correction, normalization_correction, ...
                    randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, uint16(0), pseudot, det_per_ring, options.verbose, ...
                    use_raw_data, uint32(2), ind_size, block1, blocks, index{osa_iter}, uint32(options.projector_type));
            elseif options.projector_type == 2
                [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                    zmax, options.vaimennus, normalization, SinDelayed, n_meas(end), attenuation_correction, normalization_correction, ...
                    randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, uint16(0), pseudot, det_per_ring, options.verbose, ...
                    use_raw_data, uint32(2), ind_size, block1, blocks, index{osa_iter}, uint32(options.projector_type), ...
                    options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, int32(options.accuracy_factor));
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
                    randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, LL, pseudot, det_per_ring, options.verbose, ...
                    use_raw_data, uint32(2), ind_size, block1, blocks, uint32(0), uint32(options.projector_type));
            elseif options.projector_type == 2
                [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                    zmax, options.vaimennus, normalization, SinDelayed, uint32(0), attenuation_correction, normalization_correction, ...
                    randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, LL, pseudot, det_per_ring, options.verbose, ...
                    use_raw_data, uint32(2), ind_size, block1, blocks, uint32(0), uint32(options.projector_type), options.tube_width_xy, ...
                    x_center, y_center, z_center, options.tube_width_z, int32(options.accuracy_factor));
            else
                error('Unsupported projector type')
            end
        end
        lor = reshape(lor,[],2);
        lor=repelem(int32((lor(:,1))),lor(:,2));
        
        A_length = length(rhs);
        indices=indices + 1;
        if verbose
            tStart = tic;
        end
        if options.use_fsparse == false
            A = sparse(double(lor),double(indices),double(alkiot), A_length, double(N));
        else
            A = fsparse(lor,indices,double(alkiot),[A_length double(N) length(alkiot)]);
        end
        clear indices alkiot lor
        if verbose
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
        [A, ~] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx, NSinos, NSlices, ...
            size_x, zmax, options.vaimennus, normalization, SinDelayed, n_meas(end), attenuation_correction, normalization_correction,...
            randoms_correction, lor_a, xy_index, ...
            z_index, NSinos, LL, pseudot, det_per_ring, options.verbose, ...
            use_raw_data, uint32(0), lor2, summa, false, uint32(options.projector_type), ...
            options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, int32(options.accuracy_factor));
        clear lor2
    end
    
    if size(A,2) == size(rhs,1)
        if no_norm
            bp = A * rhs;
        else
            bp = A * rhs;
            varargout{2} = full(sum(A,1));
        end
    else
        if no_norm
            bp = A' * rhs;
        else
            bp = A' * rhs;
            varargout{2} = full(sum(A,1));
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
    tube_width_xy = single(options.tube_width_xy);
    crystal_size_z = single(options.tube_width_z);
    n_rays = uint16(options.n_rays);
    dc_z = single(z_det(2,1) - z_det(1,1));
    if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
        randoms = uint32(1);
    else
        randoms = uint32(0);
    end
    if options.projector_type == 1 && options.precompute_lor
        kernel_file = 'multidevice_siddon_bpfp.cl';
        kernel_path = which(kernel_file);
        kernel_path = strrep(kernel_path, '\', '/');
        kernel_path = strrep(kernel_path, '.cl', '');
        filename = 'OMEGA_backproject_OpenCL_binary_device';
        filename = [kernel_path(1:end-length(kernel_file)), filename];
        header_directory = strrep(kernel_path,'multidevice_siddon_bpfp','');
    elseif options.projector_type == 2 && options.precompute_lor
        filename = 'OMEGA_matrix_free_orthogonal_OpenCL_binary_device';
        kernel_file = 'multidevice_orth_bpfp.cl';
        kernel_path = which(kernel_file);
        kernel_path = strrep(kernel_path, '\', '/');
        kernel_path = strrep(kernel_path, '.cl', '');
        header_directory = strrep(kernel_path,'multidevice_orth_bpfp','');
    elseif options.projector_type == 1 && ~options.precompute_lor
        kernel_file = 'multidevice_siddon_no_precomp_bpfp.cl';
        kernel_path = which(kernel_file);
        kernel_path = strrep(kernel_path, '\', '/');
        kernel_path = strrep(kernel_path, '.cl', '');
        filename = 'OMEGA_matrix_free_OpenCL_binary_device';
        header_directory = strrep(kernel_path,'multidevice_siddon_no_precomp_bpfp','');
    elseif options.projector_type == 2 && ~options.precompute_lor
        kernel_file = 'multidevice_orth_no_precomp_bpfp.cl';
        kernel_path = which(kernel_file);
        kernel_path = strrep(kernel_path, '\', '/');
        kernel_path = strrep(kernel_path, '.cl', '');
        filename = 'OMEGA_matrix_free_OpenCL_binary_device';
        header_directory = strrep(kernel_path,'multidevice_orth_no_precomp_bpfp','');
    else
        error('Invalid projector for OpenCL')
    end

    filename = [header_directory, filename];
    header_directory = strcat('-I "', header_directory);
    header_directory = strcat(header_directory,'"');
    tic
    [bp, norm] = OpenCL_matrixfree_multi_gpu( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end), ...
        single(NSlices), size_x, zmax, options.verbose, LL, pseudot, det_per_ring, uint32(options.use_device), filename, uint8(use_raw_data), ...
        single(options.cpu_to_gpu_factor), uint32(0), header_directory, options.vaimennus, normalization, n_meas(end), uint32(attenuation_correction), ...
        uint32(normalization_correction), lor_a, xy_index, z_index, tube_width_xy, crystal_size_z, x_center, y_center, z_center, SinDelayed, randoms, ...
        uint32(options.projector_type), options.precompute_lor, int32(options.accuracy_factor), n_rays, dc_z, rhs, no_norm);
    toc
    varargout{1} = norm;
end
if options.verbose
    disp('Backprojection done')
end