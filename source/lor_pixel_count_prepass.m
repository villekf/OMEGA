function [varargout] = lor_pixel_count_prepass(options)
%% Count the number of pixels each LOR traverses
% This function counts the number of pixels that each LOR traverses. This
% is needed for methods 2, 3 and 4 or for method 1 when precompute_lor has
% been selected.
% Separate codes for the sinogram and raw list-mode data.
%
% OUTPUTS:
%   lor = The number of voxels each LOR traverses (double precision)
%   lor_opencl = The number of voxels each LOR traverses (single precision)
%   lor_orth = The number of voxels each LOR traverses, orthogonal distance
%   based ray tracer. Stores the number of voxels separately for the 2D and
%   3D orthogonal ray tracer. 2D counts are stored in the earlier part, 3D
%   counts in the latter (3D applicable only if 3D orthogonal has been
%   selected) (double precision)

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

if nargout > 3
    error('Too many output arguments')
end
rings = options.rings;
diameter = options.diameter;
FOVax = options.FOVa_x;
FOVay = options.FOVa_y;
axial_fov = options.axial_fov;
Nx = options.Nx;
Ny = options.Ny;
Nz = options.Nz;
machine_name = options.machine_name;
pseudot = uint32(options.pseudot);
% NSinos = options.NSinos;
TotSinos = options.TotSinos;
det_per_ring = options.det_per_ring;

temp = pseudot;
if sum(temp) > 0
    for kk = uint32(1) : temp
        pseudot(kk) = uint32(options.cryst_per_block + 1) * kk;
    end
end

% PET-laitteen halkaisija (mm)
R=double(diameter);
% Transaxial FOV (mm)
FOVax=double(FOVax);
FOVay=double(FOVay);
% Axial FOV (mm)
axial_fov = double(axial_fov);
% Kristallien lukumäärä
blocks=uint32(rings + length(pseudot));
% ensimmäinen kristalli
block1=uint32(1);

% NSinos = uint32(NSinos);
NSlices = uint32(Nz);
TotSinos = uint32(TotSinos);

% Voxel counts
Ny=uint32(Ny);
Nx=uint32(Nx);
Nz=uint32(Nz);

folder = fileparts(which('lor_pixel_count_prepass.m'));
folder = strrep(folder, 'source','mat-files/');
folder = strrep(folder, '\','/');

%% Raw data and non-OpenCL methods
if (options.use_raw_data && (options.implementation == 1 || options.implementation == 4)) || options.precompute_all
    % Detector pair (LOR) vector
    LL = form_detector_pairs_raw(rings, det_per_ring);
    [x, y, ~, ~] = detector_coordinates(options);
    if sum(pseudot) == 0
        pseudot = [];
    end
    z_length = double(blocks) * options.cr_pz;
    zl = z_length - axial_fov;
    z = linspace(-zl/2, z_length-zl/2, blocks + 1);
    z = z(2:end) - options.cr_pz/2;
    blocks = uint32(rings);
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
    d = diff(xx(1:2));
    dy = diff(yy(1:2));
    dz = diff(zz(1:2));
    
    % Kuvan etäisyys origosta
    bx = xx(1);
    by = yy(1);
    bz = zz(1);
    
    size_x = uint32(size(x,1));
    
    zmax = max(max(z_det));
    if zmax==0
        zmax = double(1);
    end
    if options.projector_type == 1 && ~options.precompute_all
        type = uint32(0);
        x_center = xx(end);
        y_center = yy(end);
        z_center = zz(end);
    else
        if options.tube_width_z > 0 && options.implementation == 1
            type = uint32(2);
        elseif options.implementation == 1
            type = uint32(1);
        else
            type = uint(0);
        end
        x_center = xx(1:end-1)'+d/2;
        y_center = yy(1:end-1)'+dy/2;
        z_center = zz(1 : end - 1)' + dz/2;
    end
    
    % Determine which LORs go through the FOV
    LL = LL';
    LL = LL(:);
    [ lor, lor_orth] = improved_Siddon_algorithm_discard( TotSinos, Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, dy, yy, xx , TotSinos, NSlices, ...
        size_x, zmax, block1, blocks, uint32(det_per_ring), LL, pseudot, true, type, options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z);
    clear LL
    %     lor = lor(lor > 0);
    % save([folder machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'],'lor','-v7.3')
    file_string = [folder machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'];
    crystal_size_xy = options.tube_width_xy;
    crystal_size_z = options.tube_width_z;
    if exist(file_string,'file') == 2
        if options.projector_type == 1 && ~options.precompute_all
            save(file_string,'lor','-append')
        else
            variableInfo = who('-file', file_string);
            if any(strcmp(variableInfo,'lor_orth'))
                lor_file = matfile(file_string);
                if length(lor_file.lor_orth) > length(lor_orth)
                    lor_temp = load(file_string,'lor_orth');
                    lor_temp.lor_orth(1:length(lor_orth)) = lor_orth;
                    lor_orth = lor_temp.lor_orth;
                    clear lor_temp
                end
            end
            save(file_string,'lor','lor_orth','crystal_size_xy','crystal_size_z','-append')
        end
    else
        if options.projector_type == 1 && ~options.precompute_all
            save(file_string,'lor','-v7.3')
        else
            save(file_string,'lor','lor_orth','crystal_size_xy','crystal_size_z','-v7.3')
        end
    end
    if nargout >= 1
        varargout{1} = lor;
    end
    if nargout >= 3
        varargout{3} = lor_orth;
    end

end
%% Raw data, OpenCL-methods
if (options.use_raw_data && (options.implementation == 2 || options.implementation == 3 || options.implementation == 5)) ...
        || options.precompute_all
    if exist('OpenCL_matrixfree_multi_gpu','file') == 3
        % Detector pair (LOR) vector
        LL = form_detector_pairs_raw(rings, det_per_ring);
        [x, y, ~, ~] = detector_coordinates(options);
        if sum(pseudot) == 0
            pseudot = [];
        end
        z_length = single(blocks) * options.cr_pz;
        zl = z_length - axial_fov;
        z = linspace(-zl/2, z_length-zl/2, blocks + 1);
        z = z(2:end) - options.cr_pz/2;
        blocks = uint32(rings);
        if min(min(z)) == 0
            z = z + (axial_fov - max(max(z)))/2;
        end
        
        x=single(x);
        y=single(y);
        z_det = single(z);
        
        % Pikselit
        etaisyys = (R-FOVax)/2;
        xx = single(linspace(etaisyys,R-etaisyys,Nx+1));
        etaisyys = (R-FOVay)/2;
        yy = single(linspace(etaisyys,R-etaisyys,Ny+1));
        zz = linspace(single(0),single(axial_fov),Nz+1);
        zz = zz(2*block1-1:2*blocks);
        
        % Pikselien etäisyys toisistaan
        d = diff(xx(1:2));
        dy = diff(yy(1:2));
        dz = diff(zz(1:2));
        
        % Kuvan etäisyys origosta
        bx = xx(1);
        by = yy(1);
        bz = zz(1);
        
        size_x = uint32(size(x,1));
        
        zmax = max(max(z_det));
        if zmax==0
            zmax = single(1);
        end
        kernel_file = 'find_lors_opencl.cl';
        kernel_path = which(kernel_file);
        kernel_path = strrep(kernel_path, '\', '/');
        kernel_path = strrep(kernel_path, '.cl', '');
        header_directory = strrep(kernel_path,'find_lors_opencl','');
        filename = 'OMEGA_precompute_OpenCL_binary_device';
        filename = [header_directory, filename];
        header_directory = strcat('-I "', header_directory);
        header_directory = strcat(header_directory,'"');
        if isempty(pseudot)
            pseudot = uint32(100000);
        end
        
        
        LL = LL';
        LL = LL(:);
        [lor_opencl] = OpenCL_matrixfree_multi_gpu( kernel_path, Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end), ...
            single(NSlices), size_x, zmax, options.verbose, LL, pseudot, uint32(det_per_ring), uint32(options.use_device), filename, uint8(1), ...
            single(options.cpu_to_gpu_factor), uint32(2), header_directory, block1, blocks, uint32(options.NSinos), uint16(TotSinos));
        clear LL
        %     lor = lor(lor > 0);
        % save([folder machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'],'lor_opencl','-v7.3','-append')
        file_string = [folder machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'];
        if exist(file_string,'file') == 2
            save(file_string,'lor_opencl','-append')
        else
            save(file_string,'lor_opencl','-v7.3')
        end
        if nargout >= 2
            varargout{2} = lor_opencl;
        end
    else
        error('OpenCL mex-file not found. If you want to use OpenCL methods, build the OpenCL mex-files first')
    end
end
%% Sinogram data, non-OpenCL
if (~options.use_raw_data && (options.implementation == 1 || options.implementation == 4)) || options.precompute_all
    if exist([folder machine_name '_app_coordinates_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'], 'file') == 2
        load([folder machine_name '_app_coordinates_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'x','y');
    else
        sinogram_coordinates_2D(options);
        load([folder machine_name '_app_coordinates_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'x','y');
    end
    z = sinogram_coordinates_3D(options);
    
    
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
    
    size_x = uint32(size(x,1));
    
    zmax = max(max(z_det));
    if zmax==0
        zmax = double(1);
    end
    if options.projector_type == 1 && ~options.precompute_all
        type = uint32(0);
        x_center = xx(end);
        y_center = yy(end);
        z_center = zz(end);
    else
        if options.tube_width_z > 0 && options.implementation == 1
            type = uint32(2);
        elseif options.implementation == 1
            type = uint32(1);
        else
            type = uint32(0);
        end
        x_center = xx(1:end-1)'+d/2;
        y_center = yy(1:end-1)'+dy/2;
        z_center = zz(1 : end - 1)' + dz/2;
    end
    
    % Determine which LORs go through the FOV
    
    [ lor, lor_orth] = improved_Siddon_algorithm_discard( TotSinos, Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, dy, yy, xx , TotSinos, NSlices, ...
        size_x, zmax, block1, blocks, uint32(det_per_ring), uint16(0), [], false, type, options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z);
    
    file_string = [folder machine_name '_lor_pixel_count_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_sino_' num2str(options.Ndist) 'x' ...
        num2str(options.Nang) 'x' num2str(options.TotSinos) '.mat'];
    crystal_size_xy = options.tube_width_xy;
    crystal_size_z = options.tube_width_z;
    if exist(file_string,'file') == 2
        if options.projector_type == 1 && ~options.precompute_all
            save(file_string,'lor','-append')
        else
            variableInfo = who('-file', file_string);
            if any(strcmp(variableInfo,'lor_orth'))
                lor_file = matfile(file_string);
                if length(lor_file.lor_orth) > length(lor_orth)
                    lor_temp = load(file_string,'lor_orth');
                    lor_temp.lor_orth(1:length(lor_orth)) = lor_orth;
                    lor_orth = lor_temp.lor_orth;
                    clear lor_temp
                end
            end
            save(file_string,'lor','lor_orth','crystal_size_xy','crystal_size_z','-append')
        end
    else
        if options.projector_type == 1 && ~options.precompute_all
            save(file_string,'lor','-v7.3')
        else
            save(file_string,'lor','lor_orth','crystal_size_xy','crystal_size_z','-v7.3')
        end
    end
    if nargout >= 1
        varargout{1} = lor;
    end
    if nargout >= 3
        varargout{3} = lor_orth;
    end
end
%% Sinogram data, OpenCL
if (~options.use_raw_data && (options.implementation == 2 || options.implementation == 3 || options.implementation == 5)) ...
        || options.precompute_all
    if exist('OpenCL_matrixfree_multi_gpu','file') == 3
        if exist([folder machine_name '_app_coordinates_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'], 'file') == 2
            load([folder machine_name '_app_coordinates_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'x','y');
        else
            sinogram_coordinates_2D(options);
            load([folder machine_name '_app_coordinates_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'x','y');
        end
        z = sinogram_coordinates_3D(options);
        
        
        if min(min(z)) == 0
            z = z + (axial_fov - max(max(z)))/2;
        end
        
        x = single(x);
        y = single(y);
        z_det = single(z);
        
        
        
        % Pikselit
        etaisyys = (R-FOVax)/2;
        xx = single(linspace(etaisyys,R-etaisyys,Nx+1));
        etaisyys = (R-FOVay)/2;
        yy = single(linspace(etaisyys,R-etaisyys,Ny+1));
        zz = linspace(single(0),single(axial_fov),Nz+1);
        zz = zz(2*block1-1:2*blocks);
        
        % Pikselien etäisyys toisistaan
        d = diff(xx(1:2));
        dy = diff(yy(1:2));
        dz = diff(zz(1:2));
        
        % Kuvan etäisyys origosta
        bx = xx(1);
        by = yy(1);
        bz = zz(1);
        
        size_x = uint32(size(x,1));
        
        zmax = max(max(z_det));
        if zmax==0
            zmax = single(1);
        end
        kernel_file = 'find_lors_opencl.cl';
        kernel_path = which(kernel_file);
        kernel_path = strrep(kernel_path, '\', '/');
        kernel_path = strrep(kernel_path, '.cl', '');
        header_directory = strrep(kernel_path,'find_lors_opencl','');
        filename = 'OMEGA_precompute_OpenCL_binary_device';
        filename = [header_directory, filename];
        header_directory = strcat('-I "', header_directory);
        header_directory = strcat(header_directory,'"');
        if isempty(pseudot)
            pseudot = uint32(100000);
        end
        
        % Determine which LORs go through the FOV
        
        [lor_opencl] = OpenCL_matrixfree_multi_gpu( kernel_path, Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end), ...
            single(NSlices), size_x, zmax, options.verbose, uint16(0), pseudot, uint32(det_per_ring), uint32(options.use_device), filename, uint8(0),...
            single(options.cpu_to_gpu_factor), uint32(2), header_directory, block1, blocks, uint32(options.NSinos), uint16(TotSinos));
        
        file_string = [folder machine_name '_lor_pixel_count_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_sino_' num2str(options.Ndist) 'x' ...
            num2str(options.Nang) 'x' num2str(options.TotSinos) '.mat'];
        if exist(file_string,'file') == 2
            save(file_string,'lor_opencl','-append')
        else
            save(file_string,'lor_opencl','-v7.3')
        end
        if nargout >= 2
            varargout{2} = lor_opencl;
        end
    else
        error('OpenCL mex-file not found. If you want to use OpenCL methods, build the OpenCL mex-files first')
    end
end

end