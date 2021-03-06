function [varargout] = lor_pixel_count_prepass(options, varargin)
%% Count the number of voxels each LOR traverses
% This function counts the number of voxels that each LOR traverses, i.e.
% the total number of voxels along each LOR. This is needed for all
% implementations when precompute_lor = true. 
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
%   lor_vol = The number of voxels in each LOR when using the volume-based
%   ray tracer.

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

if nargout > 4
    error('Too many output arguments')
end
if nargin >= 2 && ~isempty(varargin)
    save_file = varargin{1};
else
    save_file = true;
end
if ~isfield(options,'listmode')
    options.listmode = false;
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
if options.use_raw_data
    rings = rings - sum(options.pseudot);
end

temp = pseudot;
if sum(temp) > 0
    for kk = uint32(1) : temp
        pseudot(kk) = uint32(options.cryst_per_block + 1) * kk;
    end
end
if ~isfield(options,'angles')
    options.angles = 0;
end
if ~isfield(options,'dPitch')
    options.dPitch = 0;
end
if ~isfield(options,'size_y')
    options.size_y = uint32(0);
end
if ~isfield(options,'nProjections')
    options.nProjections = int64(0);
end
if ~isfield(options,'CT')
    options.CT = false;
end

% Diameter of the device (mm)
R=double(diameter);
if isfield(options,'Z')
    Z = options.Z;
else
    Z = axial_fov;
end
% Transaxial FOV (mm)
FOVax=double(FOVax);
FOVay=double(FOVay);
% Axial FOV (mm)
axial_fov = double(axial_fov);
% Number of axial rings
blocks=uint32(rings - 1);
% First ring
block1=uint32(0);

% NSinos = uint32(NSinos);
NSlices = uint32(Nz);
TotSinos = uint32(TotSinos);

% Voxel counts
Ny=uint32(Ny);
Nx=uint32(Nx);
Nz=uint32(Nz);

if options.CT
    if isfield(options,'x') && isfield(options,'y')
        x = options.x;
        y = options.y;
        if isfield(options,'z')
            z = options.z;
        else
            z = options.z_det;
        end
    else
        [x,y,z] = CTDetectorCoordinates(options.angles,options.sourceToDetector,options.sourceToCRot,options.dPitch,options.xSize,...
            options.ySize,options.horizontalOffset,options.verticalOffset,options.bedOffset);
    end
    if numel(z)/2 > numel(options.angles)
        if size(options.angles,1) == 1
            options.angles = reshape(options.angles, [],1);
        end
        options.angles = repmat(options.angles,numel(z)/2/numel(options.angles),1);
    end
end

if options.CT
    R = 0;
    if options.listmode
        Z = options.dPitch * double(options.xSize) + z(numel(z)/2);
    else
        Z = options.dPitch * double(options.xSize) + z(options.nProjections);
    end
    if isfield(options,'Z')
        Z = options.Z;
    end
end
[xx,yy,zz,dx,dy,dz,bx,by,bz] = computePixelSize(R, FOVax, FOVay, Z, axial_fov, Nx, Ny, Nz, options.implementation);

folder = fileparts(which('lor_pixel_count_prepass.m'));
folder = [folder(1:end-6), 'mat-files/'];
folder = strrep(folder, '\','/');

if exist('feature','builtin') == 5
    nCores = uint32(feature('numcores'));
else
    nCores = uint32(1);
end

lor = [];
lor_orth = [];
lor_vol = [];
lor_opencl = [];

%% Raw data and non-OpenCL methods
if (options.use_raw_data && (options.implementation == 1 || options.implementation == 4)) || options.precompute_all
    % Detector coordinates
    [x, y, z] = get_coordinates(options, blocks, pseudot, false);
    % Detector pair (LOR) vector
    if ~save_file
        LL = zeros(numel(x),1,'uint16');
    else
        LL = form_detector_pairs_raw(rings, det_per_ring);
    end
    if sum(pseudot) == 0
        pseudot = [];
    end
    
%     blocks = uint32(rings);
    
    x = double(x(:));
    y = double(y(:));
    z_det = double(z(:));
%     zz=zz(2*block1-1:2*blocks);
    
    size_x = uint32(numel(x)) / 2;
    
    if isfield(options,'zmax')
        zmax = options.zmax;
    else
        zmax = max(max(z_det));
        if zmax==0
            zmax = double(1);
        end
    end
    % Voxel center coordinates are needed for orthogonal/distance-based ray
    % tracer 
    if options.projector_type == 1
        type = uint32(0);
        x_center = xx(end);
        y_center = yy(end);
        z_center = zz(end);
        V = 0;
        Vmax = 0;
        bmin = 0;
        bmax = 0;
    elseif options.projector_type == 2
        if options.tube_width_z > 0 && options.implementation == 1
            type = uint32(2);
        elseif options.implementation == 1 && (options.tube_width_z > 0 || options.tube_width_xy > 0)
            type = uint32(1);
        else
            type = uint(0);
        end
        x_center = xx(1:end-1)'+dx/2;
        y_center = yy(1:end-1)'+dy/2;
        z_center = zz(1 : end - 1)' + dz/2;
        V = 0;
        Vmax = 0;
        bmin = 0;
        bmax = 0;
    elseif options.projector_type == 3
        if options.implementation == 1
            type = uint32(3);
        else
            type = uint32(0);
        end
        x_center = xx(1:end-1)'+dx/2;
        y_center = yy(1:end-1)'+dy/2;
        z_center = zz(1 : end - 1)' + dz/2;
        options.voxel_radius = (sqrt(2) * options.voxel_radius * dx) / 2;
        bmax = options.tube_radius + options.voxel_radius;
        b = linspace(0, bmax, 10000)';
        b(options.tube_radius > (b + options.voxel_radius)) = [];
        b = unique(round(b*10^3)/10^3);
        V = volumeIntersection(options.tube_radius, options.voxel_radius, b);
        Vmax = (4*pi)/3*options.voxel_radius^3;
        bmin = min(b);
    end
    
    % Determine which LORs go through the FOV
    LL = LL';
    LL = LL(:);
    if exist('OCTAVE_VERSION','builtin') == 0
        [ lor, lor_orth, lor_vol] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx, TotSinos, NSlices, size_x, ...
            zmax, 0, 0, 0, uint32(0), false, false, false, false, 0, 0, uint16(0), uint32(0), uint16(0), TotSinos, LL, pseudot, uint32(det_per_ring), ...
            false, int64(0), 0, 0, int64(0), uint32(0), options.verbose, nCores, ...
            options.use_raw_data, uint32(3), options.listmode, block1, blocks, uint32(options.projector_type), options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, ...
            bmin, bmax, Vmax, V, type);
    elseif exist('OCTAVE_VERSION','builtin') == 5
        [ lor, lor_orth, lor_vol] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx, TotSinos, NSlices, size_x, ...
            zmax, 0, 0, 0, uint32(0), false, false, false, false, 0, 0, uint16(0), uint32(0), uint16(0), TotSinos, LL, pseudot, uint32(det_per_ring), ...
            false, int64(0), 0, 0, int64(0), uint32(0), options.verbose, nCores, ...
            options.use_raw_data, uint32(3), options.listmode, block1, blocks, uint32(options.projector_type), options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, ...
            bmin, bmax, Vmax, V, type);
    end
    clear LL
    % Save the data
    file_string = [folder machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'];
    crystal_size_xy = options.tube_width_xy;
    crystal_size_z = options.tube_width_z;
    if save_file
        if exist(file_string,'file') == 2
            if type == 0
                if exist('OCTAVE_VERSION', 'builtin') == 0
                    save(file_string,'lor','-append')
                else
                    save(file_string,'lor','-append','-v7')
                end
            elseif type == 1 || type == 2
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
                if exist('OCTAVE_VERSION', 'builtin') == 0
                    save(file_string,'lor','lor_orth','crystal_size_xy','crystal_size_z','-append')
                else
                    save(file_string,'lor','lor_orth','crystal_size_xy','crystal_size_z','-append','-v7')
                end
            elseif type == 3
                if exist('OCTAVE_VERSION', 'builtin') == 0
                    save(file_string,'lor','lor_vol','bmax','Vmax','-append')
                else
                    save(file_string,'lor','lor_vol','bmax','Vmax','-append','-v7')
                end
            end
        else
            if type == 0
                if exist('OCTAVE_VERSION', 'builtin') == 0
                    save(file_string,'lor','-v7.3')
                else
                    save(file_string,'lor','-v7')
                end
            elseif type == 1 || type == 2
                if exist('OCTAVE_VERSION', 'builtin') == 0
                    save(file_string,'lor','lor_orth','crystal_size_xy','crystal_size_z','-v7.3')
                else
                    save(file_string,'lor','lor_orth','crystal_size_xy','crystal_size_z','-v7')
                end
            elseif type == 3
                if exist('OCTAVE_VERSION', 'builtin') == 0
                    save(file_string,'lor','lor_vol','bmax','Vmax','-v7.3')
                else
                    save(file_string,'lor','lor_vol','bmax','Vmax','-v7')
                end
            end
        end
    end
    if nargout >= 1
        varargout{1} = lor;
    end
    if nargout >= 2
        varargout{2} = [];
    end
    if nargout >= 3
        varargout{3} = lor_orth;
    end
    if nargout >= 4
        varargout{4} = lor_vol;
    end

end
%% Raw data, OpenCL-methods
if (options.use_raw_data && (options.implementation == 2 || options.implementation == 3 || options.implementation == 5)) ...
        || options.precompute_all
    if exist('OpenCL_matrixfree_multi_gpu','file') == 3
        % Detector pair (LOR) vector
        LL = form_detector_pairs_raw(rings, det_per_ring);
        if sum(pseudot) == 0
            pseudot = [];
        end
        [x, y, z] = get_coordinates(options, blocks, pseudot, false);
        
%         blocks = uint32(rings);
        
        x=single(x);
        y=single(y);
        z_det = single(z);
        
        size_x = uint32(size(x,1));
        
        zmax = max(max(z_det));
        if zmax==0
            zmax = single(1);
        end
        kernel_file = 'multidevice_kernel.cl';
        kernel_path = which(kernel_file);
        kernel_path = strrep(kernel_path, '\', '/');
        kernel_path = strrep(kernel_path, '.cl', '');
        header_directory = strrep(kernel_path,'multidevice_kernel','');
        filename = 'OMEGA_precompute_OpenCL_binary_device';
        filename = [header_directory, filename];
%         header_directory = strcat('-I "', header_directory);
%         header_directory = strcat(header_directory,'"');
        if isempty(pseudot)
            pseudot = uint32(100000);
        end
        
        
        LL = LL';
        LL = LL(:);
        if options.implementation == 3
            [lor_opencl] = OpenCL_matrixfree_multi_gpu( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end), ...
                single(NSlices), size_x, zmax, options.verbose, LL, pseudot, uint32(det_per_ring), false, int64(0), 0, 0, int64(0), int32(0), uint32(options.use_device), filename, uint8(1), ...
                single(options.cpu_to_gpu_factor), uint32(2), header_directory, block1, blocks, uint32(options.NSinos), uint16(TotSinos));
        elseif options.implementation == 2
            lor_opencl = OpenCL_matrixfree( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end), ...
                TotSinos, single(NSlices), size_x, zmax, TotSinos, options.verbose, LL, pseudot, uint32(det_per_ring), false, int64(0), 0, 0, int64(0), int32(0), uint32(options.use_device), uint8(1), ...
                filename, uint32(2), false, header_directory, block1, blocks, uint32(options.NSinos), uint16(TotSinos));
        end
        clear LL
        %     lor = lor(lor > 0);
        % save([folder machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'],'lor_opencl','-v7.3','-append')
        file_string = [folder machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'];
        if save_file
            if exist(file_string,'file') == 2
                if exist('OCTAVE_VERSION', 'builtin') == 0
                    save(file_string,'lor_opencl','-append')
                else
                    save(file_string,'lor_opencl','-append','-v7')
                end
            else
                if exist('OCTAVE_VERSION', 'builtin') == 0
                    save(file_string,'lor_opencl','-v7.3')
                else
                    save(file_string,'lor_opencl','-v7')
                end
            end
            if nargout >= 2
                varargout{2} = lor_opencl;
            end
        end
    else
        error('OpenCL mex-file not found. If you want to use OpenCL methods, build the OpenCL mex-files first')
    end
end
%% Sinogram data, non-OpenCL
if (~options.use_raw_data && (options.implementation == 1 || options.implementation == 4)) || options.precompute_all
    if ~options.CT
        [x, y, z] = get_coordinates(options, options.rings - 1, [], false);
    end
    
    x = double(x(:));
    y = double(y(:));
    z_det = double(z(:));
    
    if options.CT
        if options.listmode
            size_x = uint32(numel(x)) / 2;
        else
            size_x = uint32(options.ySize);
        end
    else
        size_x = uint32(numel(x)) / 2;
    end
    if ~save_file
        TotSinos = uint32(1);
    end
    if options.CT && ~options.listmode
        TotSinos = uint32(options.nProjections * options.xSize);
    end
    
    if isfield(options,'zmax')
        zmax = options.zmax;
    else
        zmax = max(max(z_det));
        if zmax==0
            zmax = double(1);
        end
        if options.CT
            zmax = z_det(options.nProjections) + options.dPitch * (options.xSize - 1);
        end
    end
    if options.projector_type == 1
        type = uint32(0);
        x_center = xx(end);
        y_center = yy(end);
        z_center = zz(end);
        V = 0;
        Vmax = 0;
        bmin = 0;
        bmax = 0;
    elseif options.projector_type == 2
        if options.tube_width_z > 0 && options.implementation == 1
            type = uint32(2);
        elseif options.implementation == 1
            type = uint32(1);
        else
            type = uint32(0);
        end
        x_center = xx(1 : end - 1)' + dx/2;
        y_center = yy(1 : end - 1)' + dy/2;
        z_center = zz(1 : end - 1)' + dz/2;
        V = 0;
        Vmax = 0;
        bmin = 0;
        bmax = 0;
    elseif options.projector_type == 3
        if options.implementation == 1
            type = uint32(3);
        else
            type = uint32(0);
        end
        x_center = xx(1 : end - 1)' + dx/2;
        y_center = yy(1 : end - 1)' + dy/2;
        z_center = zz(1 : end - 1)' + dz/2;
        dp = max([dx,dy,dz]);
        options.voxel_radius = (sqrt(2) * options.voxel_radius * dp) / 2;
        bmax = options.tube_radius + options.voxel_radius;
        b = linspace(0, bmax, 10000)';
        b(options.tube_radius > (b + options.voxel_radius)) = [];
        b = unique(round(b*10^3)/10^3);
        V = volumeIntersection(options.tube_radius, options.voxel_radius, b);
        Vmax = (4*pi)/3*options.voxel_radius^3;
        bmin = min(b);
    end
    
    % Determine which LORs go through the FOV
    
    if options.CT
        if exist('OCTAVE_VERSION','builtin') == 0
            [ lor, lor_orth, lor_vol] = projector_mexCT( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx, TotSinos, NSlices, size_x, ...
                zmax, 0, 0, 0, uint32(0), false, false, false, false, 0, 0, uint16(0), uint32(0), uint16(0), TotSinos, uint16(0), pseudot, uint32(det_per_ring), ...
                false, int64(0), 0, 0, int64(0), uint32(0), options.verbose, nCores, ...
                options.use_raw_data, uint32(3), options.listmode, block1, blocks, uint32(options.projector_type), options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, ...
                bmin, bmax, Vmax, V, type, uint32(1), options.angles, uint32(options.xSize), options.dPitch, ...
                int64(numel(options.angles)));
        elseif exist('OCTAVE_VERSION','builtin') == 5
            [ lor, lor_orth, lor_vol] = projector_octCT( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx, TotSinos, NSlices, size_x, ...
                zmax, 0, 0, 0, uint32(0), false, false, false, false, 0, 0, uint16(0), uint32(0), uint16(0), TotSinos, uint16(0), pseudot, uint32(det_per_ring), ...
                false, int64(0), 0, 0, int64(0), uint32(0), options.verbose, nCores, ...
                options.use_raw_data, uint32(3), options.listmode, block1, blocks, uint32(options.projector_type), options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, ...
                bmin, bmax, Vmax, V, type, uint32(1), options.angles, uint32(options.xSize), options.dPitch, ...
                int64(numel(options.angles)));
        end
    else
        if exist('OCTAVE_VERSION','builtin') == 0
            [ lor, lor_orth, lor_vol] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx, TotSinos, NSlices, size_x, ...
                zmax, 0, 0, 0, uint32(0), false, false, false, false, 0, 0, uint16(0), uint32(0), uint16(0), TotSinos, uint16(0), pseudot, uint32(det_per_ring), ...
                false, int64(0), 0, 0, int64(0), uint32(0), options.verbose, nCores, ...
                options.use_raw_data, uint32(3), options.listmode, block1, blocks, uint32(options.projector_type), options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, ...
                bmin, bmax, Vmax, V, type);
        elseif exist('OCTAVE_VERSION','builtin') == 5
            [ lor, lor_orth, lor_vol] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx, TotSinos, NSlices, size_x, ...
                zmax, 0, 0, 0, uint32(0), false, false, false, false, 0, 0, uint16(0), uint32(0), uint16(0), TotSinos, uint16(0), pseudot, uint32(det_per_ring), ...
                false, int64(0), 0, 0, int64(0), uint32(0), options.verbose, nCores, ...
                options.use_raw_data, uint32(3), options.listmode, block1, blocks, uint32(options.projector_type), options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, ...
                bmin, bmax, Vmax, V, type);
        end
    end
    
    file_string = [folder machine_name '_lor_pixel_count_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_sino_' num2str(options.Ndist) 'x' ...
        num2str(options.Nang) 'x' num2str(options.TotSinos) '.mat'];
    crystal_size_xy = options.tube_width_xy;
    crystal_size_z = options.tube_width_z;
    if save_file
        if exist(file_string,'file') == 2
            if type == 0
                if exist('OCTAVE_VERSION', 'builtin') == 0
                    save(file_string,'lor','-append')
                else
                    save(file_string,'lor','-append','-v7')
                end
            elseif type == 1 || type == 2
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
                if exist('OCTAVE_VERSION', 'builtin') == 0
                    save(file_string,'lor','lor_orth','crystal_size_xy','crystal_size_z','-append')
                else
                    save(file_string,'lor','lor_orth','crystal_size_xy','crystal_size_z','-append','-v7')
                end
            elseif type == 3
                if exist('OCTAVE_VERSION', 'builtin') == 0
                    save(file_string,'lor','lor_vol','bmax','Vmax','-append')
                else
                    save(file_string,'lor','lor_vol','bmax','Vmax','-append','-v7')
                end
            end
        else
            if type == 0
                if exist('OCTAVE_VERSION', 'builtin') == 0
                    save(file_string,'lor','-v7.3')
                else
                    save(file_string,'lor','-v7')
                end
            elseif type == 1 || type == 2
                if exist('OCTAVE_VERSION', 'builtin') == 0
                    save(file_string,'lor','lor_orth','crystal_size_xy','crystal_size_z','-v7.3')
                else
                    save(file_string,'lor','lor_orth','crystal_size_xy','crystal_size_z','-v7')
                end
            elseif type == 3
                if exist('OCTAVE_VERSION', 'builtin') == 0
                    save(file_string,'lor','lor_vol','bmax','Vmax','-v7.3')
                else
                    save(file_string,'lor','lor_vol','bmax','Vmax','-v7')
                end
            end
        end
    end
    if nargout >= 1
        varargout{1} = lor;
    end
    if nargout >= 3
        varargout{3} = lor_orth;
    end
    if nargout >= 4
        varargout{4} = lor_vol;
    end
end
%% Sinogram data, OpenCL
if (~options.use_raw_data && (options.implementation == 2 || options.implementation == 3 || options.implementation == 5)) ...
        || options.precompute_all
    if exist('OpenCL_matrixfree_multi_gpu','file') == 3
        [x, y, z] = get_coordinates(options, blocks, [], false);
        
        x = single(x);
        y = single(y);
        z_det = single(z);
        
        size_x = uint32(size(x,1));
        
        zmax = max(max(z_det));
        if zmax==0
            zmax = single(1);
        end
        kernel_file = 'multidevice_kernel.cl';
        kernel_path = which(kernel_file);
        kernel_path = strrep(kernel_path, '\', '/');
        kernel_path = strrep(kernel_path, '.cl', '');
        header_directory = strrep(kernel_path,'multidevice_kernel','');
        filename = 'OMEGA_precompute_OpenCL_binary_device';
        filename = [header_directory, filename];
%         header_directory = strcat('-I "', header_directory);
%         header_directory = strcat(header_directory,'"');
        if isempty(pseudot)
            pseudot = uint32(100000);
        end
        
        % Determine which LORs go through the FOV
        
        if options.implementation == 3
            [lor_opencl] = OpenCL_matrixfree_multi_gpu( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end), ...
                single(NSlices), size_x, zmax, options.verbose, uint16(0), pseudot, uint32(det_per_ring), false, int64(0), 0, 0, int64(0), int32(0), uint32(options.use_device), filename, uint8(0),...
                single(options.cpu_to_gpu_factor), uint32(2), header_directory, block1, blocks, uint32(options.NSinos), uint16(TotSinos));
        elseif options.implementation == 2
            [lor_opencl] = OpenCL_matrixfree( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end), ...
                TotSinos, single(NSlices), size_x, zmax, TotSinos, options.verbose, uint16(0), pseudot, uint32(det_per_ring), false, int64(0), 0, 0, int64(0), int32(0), uint32(options.use_device), uint8(0),...
                filename, uint32(2), true, header_directory, block1, blocks, uint32(options.NSinos), uint16(TotSinos));
        end
        
        file_string = [folder machine_name '_lor_pixel_count_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_sino_' num2str(options.Ndist) 'x' ...
            num2str(options.Nang) 'x' num2str(options.TotSinos) '.mat'];
        if save_file
            if exist(file_string,'file') == 2
                if exist('OCTAVE_VERSION', 'builtin') == 0
                    save(file_string,'lor_opencl','-append')
                else
                    save(file_string,'lor_opencl','-append','-v7')
                end
            else
                if exist('OCTAVE_VERSION', 'builtin') == 0
                    save(file_string,'lor_opencl','-v7.3')
                else
                    save(file_string,'lor_opencl','-v7')
                end
            end
        end
        if nargout >= 2
            varargout{2} = lor_opencl;
        end
    else
        error('OpenCL mex-file not found. If you want to use OpenCL methods, build the OpenCL mex-files first')
    end
end

end