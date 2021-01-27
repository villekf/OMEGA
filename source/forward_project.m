function varargout = forward_project(options, index, n_meas, f, nn, iternn, varargin)
%FORWARD_PROJECT Calculates the forward projection
% Examples:
%   fp = forward_project(options, index, n_meas, f, nn)
%   [fp, options] = forward_project(options, index, n_meas, f, nn)
% INPUTS:
%   options = Contains the machine and (optionally) sinogram specific
%   information (e.g. detectors per ring, image (matrix) size, FOV size)
%   index = The indices (LORs) used to compute the system matrix (you can
%   use index_maker to produce the indices)
%   n_meas = Number of measurements used
%   f = The current estimate
%   nn = The interval from where the measurements are taken (current
%   subset)
%   iternn = Current subset and iteration numbers summed minus 1. E.g. this
%   is 1 if it is the very first iteration and subset
%   SinM = The measurements, used to skip zero measurements. If you want to
%   include also zero measurements, then input a vector of ones.
%
% OUTPUTS:
%   fp = The forward projection (fp = A * f)
%   options = If randoms, scatter or normalization correction is used, then
%   they are stored in the options variable at sub-iteration 1
%
% See also index_maker, backproject

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


if nargout > 2
    error('Too many output arguments')
end
if nargout == 0
    error('Too few output arguments')
end

if nargin >= 7 && ~isempty(varargin{1}) && islogical(varargin{1})
    luokka = varargin{1};
else
    luokka = false;
end
if nargin >= 8 && ~isempty(varargin{2}) && islogical(varargin{2}) && options.implementation == 1
    store_matrix = varargin{1};
else
    store_matrix = false;
end

% folder = fileparts(which('reconstructions_main.m'));
% folder = strrep(folder, 'source','mat-files/');
% folder = strrep(folder, '\','/');

if ~isfield(options,'TOF_bins')
    options.TOF_bins = 1;
end
TOF = options.TOF_bins > 1 && options.projector_type == 1 && options.implementation > 1;
if ~isfield(options,'simple')
    options.simple = false;
end
if ~isfield(options,'listmode')
    options.listmode = false;
end
if ~isfield(options,'compute_sensitivity_image')
    options.compute_sensitivity_image = false;
end

f = f(:);

rings = options.rings;
Nx = options.Nx;
Ny = options.Ny;
Nz = options.Nz;
if Nz > options.NSinos && ~options.simple
    Nz = options.NSinos;
    rings = ceil(Nz / 2 + 0.25);
end
NSlices = uint32(Nz);
attenuation_correction = options.attenuation_correction;
FOVax = options.FOVa_x;
FOVay = options.FOVa_y;
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
if ~options.listmode
    Z = options.axial_fov;
else
    if abs(min(options.x(:))) < abs(max(options.x(:))) / 2 && options.diameter == 0
        options.diameter = (min(options.x(:))) + (max(options.x(:)));
    end
    if abs(min(options.z_det(:))) < abs(max(options.z_det(:))) / 2
        if min(options.z_det(:)) < 0
            Z = options.axial_fov - min(options.z_det(:)) * 2;
        elseif max(options.z_det(:)) > options.axial_fov
            Z = options.axial_fov + (max(options.z_det(:)) - options.axial_fov) * 2;
        else
            Z = options.axial_fov;
        end
    else
        Z = 0;
    end
    if options.implementation == 4 || options.implementation == 1
        use_raw_data = true;
        options.use_raw_data = use_raw_data;
    end
end
% Number of rings
blocks=uint32(rings + length(pseudot) - 1);
% First ring
block1=uint32(0);

save_norm = false;
save_rand = false;
save_scat = false;

if exist('feature','builtin') == 5
    nCores = uint32(feature('numcores'));
else
    nCores = uint32(1);
end

NSinos = uint32(options.NSinos);
% TotSinos = int32(options.TotSinos);

if numel(f) ~= Nx*Ny*Nz && ~isempty(f)
    error('Estimate has different amount of elements than the image size')
end

if ~isa(n_meas,'int64')
    n_meas = int64(n_meas);
end

if TOF
    if options.TOF_FWHM < 0
        if iscell(options.SinM)
            for kk = 1 : options.partitions
                options.SinM{kk} = sum(options.SinM{kk},4,'native');
            end
        else
            options.SinM = sum(options.SinM,4,'native');
        end
        sigma_x = 0;
        TOFCenter = 0;
        TOF = false;
        options.TOF_bins = 1;
    else
        c = 2.99792458e11;
        sigma_x = (c*options.TOF_FWHM/2) / (2 * sqrt(2 * log(2)));
        edges_user = linspace(-options.TOF_width * options.TOF_bins/2, options.TOF_width * options.TOF_bins / 2, options.TOF_bins + 1);
        edges_user = edges_user(1:end-1) + options.TOF_width/2; % the most probable value where annihilation occured
        TOFCenter = zeros(size(edges_user));
        TOFCenter(1) = edges_user(ceil(length(edges_user)/2));
        TOFCenter(2:2:end) = edges_user(ceil(length(edges_user)/2) + 1:end);
        TOFCenter(3:2:end) = edges_user(ceil(length(edges_user)/2) - 1: -1 : 1);
        TOFCenter = TOFCenter * c / 2;
        if isfield(options, 'TOF_offset') && options.TOF_offset > 0
            TOFCenter = TOFCenter + options.TOF_offset;
        end
    end
else
    sigma_x = 0;
    TOFCenter = 0;
end
if options.implementation == 3 || options.implementation == 2
    sigma_x = single(sigma_x);
    TOFCenter = single(TOFCenter);
end

if ~luokka
    if iternn == 1 || (options.implementation > 1 && (options.n_rays_transaxial > 1 || options.n_rays_axial > 1) && ~options.precompute_lor && options.projector_type == 1)
        if (options.n_rays_transaxial > 1 || options.n_rays_axial > 1) && isfield(options,'x') && isfield(options,'y')
            options = rmfield(options, 'x');
            options = rmfield(options, 'y');
        end
        [x, y, z_det, options] = get_coordinates(options, blocks, pseudot);
        options.x = x(:);
        options.y = y(:);
        options.z_det = z_det(:);
    else
        x = options.x(:);
        y = options.y(:);
        z_det = options.z_det(:);
    end
    if abs(min(options.x(:))) < abs(max(options.x(:))) / 2 && options.diameter == 0
        R = (min(options.x(:))) + (max(options.x(:)));
    elseif abs(min(options.x(:))) < abs(max(options.x(:))) / 2 && options.diameter > 0
        R = options.diameter;
    else
        R = 0;
    end
    if isfield(options,'z')
        if abs(min(options.z(:))) < abs(max(options.z(:))) / 2
            if min(options.z(:)) < 0
                Z = options.axial_fov - min(options.z(:)) * 2;
            elseif max(options.z(:)) > options.axial_fov
                Z = options.axial_fov + (max(options.z(:)) - options.axial_fov) * 2;
            else
                Z = options.axial_fov;
            end
        else
            Z = 0;
        end
    else
        if abs(min(options.z_det(:))) < abs(max(options.z_det(:))) / 2
            if min(options.z_det(:)) < 0
                Z = options.axial_fov - min(options.z_det(:)) * 2;
            elseif max(options.z_det(:)) > options.axial_fov
                Z = options.axial_fov + (max(options.z_det(:)) - options.axial_fov) * 2;
            else
                Z = options.axial_fov;
            end
        else
            Z = 0;
        end
    end
    [xx,yy,zz,dx,dy,dz,bx,by,bz] = computePixelSize(R, FOVax, FOVay, Z, options.axial_fov, Nx, Ny, Nz, options.implementation);
    
    if ~isfield(options,'normalization')
        save_norm = true;
    end
    if ~isfield(options,'SinDelayed')
        save_rand = true;
    end
    if ~isfield(options,'ScatterC')
        save_scat = true;
    end
else
    if options.listmode
        x = options.x(nn(1):nn(2),:);
        y = options.y(nn(1):nn(2),:);
        z_det = options.z_det(nn(1):nn(2),:);
        x = x(:);
        y = y(:);
        z_det = z_det(:);
    else
        x = options.x(:);
        y = options.y(:);
        z_det = options.z_det(:);
    end
    dx = options.dx;
    dy = options.dy;
    dz = options.dz;
    bx = options.bx;
    by = options.by;
    bz = options.bz;
    xx = options.xx;
    yy = options.yy;
    zz = options.zz;
end
if use_raw_data && isfield(options,'x')
    if options.listmode
        det_per_ring = numel(x) / 2;
        options.det_per_ring = det_per_ring;
    else
        det_per_ring = numel(x);
    end
else
    det_per_ring = options.det_per_ring;
end

if isfield(options,'norm_full')
    options.normalization = options.norm_full;
end
if isfield(options,'rand_full')
    options.SinDelayed = options.rand_full;
end
if isfield(options,'scat_full')
    options.ScatterC = options.scat_full;
end

if ~luokka
    [normalization_correction, randoms_correction, options] = set_up_corrections(options, blocks);
    
    
    if save_norm
        options.norm_full = options.normalization;
    end
    if save_rand
        options.rand_full = options.SinDelayed;
    end
    if save_scat
        options.scat_full = options.ScatterC;
    end
else
    normalization_correction = options.normalization_correction;
    randoms_correction = options.randoms_correction;
end

if options.use_raw_data
    if options.listmode
        size_x = uint32(numel(x) / 2);
    else
        size_x = uint32(options.det_w_pseudo);
    end
else
    if options.listmode
        size_x = uint32(numel(x) / 2);
    else
        size_x = uint32(options.Nang*options.Ndist);
    end
    if isfield(options, 'sampling') && options.sampling > 1 && ~options.precompute_lor
        size_x = size_x * options.sampling;
    end
end

if (options.precompute_lor || options.implementation == 5 || options.implementation == 2 || options.implementation == 3)
    n_meas = [int64(0);int64(cumsum(n_meas))];
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

if ~luokka
    [options, lor_a, xy_index, z_index, LL, summa, n_meas, ~, lor_orth] = form_subset_indices(options, n_meas, 1, index, size_x, y, z_det, blocks, true, TOF);
    if ~options.precompute_lor
        lor_a = uint16(0);
        lor_orth = uint16(0);
    end
    if normalization_correction
        if options.implementation == 1 || options.implementation == 4
            normalization = double(options.normalization);
        else
            normalization = single(options.normalization);
        end
    else
        if options.implementation == 1 || options.implementation == 4
            normalization = 0;
        else
            normalization = single(0);
        end
    end
    
    if randoms_correction
        if iscell(options.SinDelayed)
            if options.implementation == 1 || options.implementation == 4
                SinDelayed = double(options.SinDelayed{1});
            else
                SinDelayed{1} = single(options.SinDelayed{1});
            end
        else
            if options.implementation == 1 || options.implementation == 4
                SinDelayed = double(options.SinDelayed);
            else
                SinDelayed{1} = single(options.SinDelayed);
            end
        end
    else
        if options.implementation == 1 || options.implementation == 4
            SinDelayed = 0;
        else
            SinDelayed{1} = single(0);
        end
    end
    
    if options.scatter_correction && ~options.subtract_scatter
        if options.implementation == 1 || options.implementation == 4
            scatter_input = double(options.ScatterC);
        else
            if iscell(options.ScatterFB)
                options.ScatterFB{1} = {single(options.ScatterC{1})};
            else
                options.ScatterFB{1} = {single(options.ScatterC)};
            end
        end
    else
        if options.implementation == 1 || options.implementation == 4
            scatter_input = 0;
        else
            options.ScatterFB{1} = single(0);
        end
    end
else
    if normalization_correction
        if options.implementation == 1 || options.implementation == 4
            normalization = double(options.normalization(nn(1) : nn(2)));
        else
            normalization = single(options.normalization(nn(1) : nn(2)));
        end
    else
        if options.implementation == 1 || options.implementation == 4
            normalization = 0;
        else
            normalization = single(0);
        end
    end
    
    if randoms_correction
        if iscell(options.SinDelayed)
            if options.implementation == 1 || options.implementation == 4
                SinDelayed = double(options.SinDelayed{1}(nn(1) : nn(2)));
            else
                SinDelayed{1} = single(options.SinDelayed{1}(nn(1) : nn(2)));
            end
        else
            if options.implementation == 1 || options.implementation == 4
                SinDelayed = double(options.SinDelayed(nn(1) : nn(2)));
            else
                SinDelayed{1} = single(options.SinDelayed(nn(1) : nn(2)));
            end
        end
    else
        if options.implementation == 1 || options.implementation == 4
            SinDelayed = 0;
        else
            SinDelayed{1} = single(0);
        end
    end
    
    if options.scatter_correction && ~options.subtract_scatter
        if options.implementation == 1 || options.implementation == 4
            scatter_input = options.ScatterC(nn(1) : nn(2));
        else
            if iscell(options.ScatterFB)
                options.ScatterFB{1} = {single(options.ScatterC{1}(nn(1) : nn(2)))};
            else
                options.ScatterFB{1} = {single(options.ScatterC(nn(1) : nn(2)))};
            end
        end
    else
        if options.implementation == 1 || options.implementation == 4
            scatter_input = 0;
        else
            options.ScatterFB{1} = single(0);
        end
    end
    if options.use_raw_data && ~options.listmode
        LL = options.LL(nn(1) : nn(2));
    elseif ~options.listmode
        xy_index = options.xy_index(nn(1) : nn(2));
        z_index = options.z_index(nn(1) : nn(2));
        LL = uint16(0);
    else
        LL = options.LL;
        xy_index = options.xy_index;
        z_index = options.z_index;
    end
    if options.implementation == 1 || (options.implementation == 4 && ~options.listmode)
        summa = options.summa(iternn);
    else
        summa = options.summa;
    end
    if options.precompute_lor
        if options.projector_type == 2 || options.projector_type == 3
            lor_orth = options.lor_orth(nn(1) : nn(2));
        else
            lor_orth = uint16(0);
        end
        lor_a = options.lor_a(nn(1) : nn(2));
    else
        lor_a = uint16(0);
        lor_orth = uint16(0);
    end
end

% Number of pixels
Ny=uint32(Ny);
Nx=uint32(Nx);
Nz=uint32(Nz);

N=(Nx)*(Ny)*(Nz);
det_per_ring = uint32(det_per_ring);

% How much memory is preallocated
if options.implementation == 1
    if use_raw_data == false
        ind_size = uint32(NSinos/8*(det_per_ring)* Nx * (Ny));
    else
        ind_size = uint32((det_per_ring)^2/8* Nx * (Ny));
    end
end


zmax = max(max(z_det));
if zmax==0
    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
        zmax = single(1);
    else
        zmax = double(1);
    end
end
if ~luokka
    [x_center,y_center,z_center,dec] = computePixelCenters(xx,yy,zz,dx,dy,dz,TOF,options);
    
    [V,Vmax,bmin,bmax] = computeVoxelVolumes(dx,dy,dz,options);
    % Multi-ray Siddon
    if options.implementation > 1 && options.n_rays_transaxial > 1 && ~options.precompute_lor && options.projector_type == 1
        [x,y] = getMultirayCoordinates(options);
        options.x = x;
        options.y = y;
    end
else
    x_center = options.x_center;
    y_center = options.y_center;
    z_center = options.z_center;
    dec = options.dec;
    V = options.V;
    Vmax = options.Vmax;
    bmin = options.bmin;
    bmax = options.bmax;
end

if options.implementation == 1
    if numel(n_meas) == 1
        n_meas = [int64(0);n_meas];
    end
    if options.precompute_lor
        is_transposed = true;
    else
        is_transposed = false;
    end
    iij = double(0:Nx);
    jji = double(0:Ny);
    kkj = double(0:Nz);
    if options.listmode
        temp_x = options.x;
        temp_y = options.y;
        temp_z = options.z_det;
        options.x = x;
        options.y = y;
        options.z_det = z_det;
        options.Z = Z;
        options.precompute_all = false;
        lor_a = lor_pixel_count_prepass(options, false);
        options.x = temp_x;
        options.y = temp_y;
        options.z_det = temp_z;
        summa = uint64(sum(lor_a));
    end
    options.attenuation_phase = false;
    [A,~] = computeImplementation1(options,use_raw_data,randoms_correction, n_meas,1, normalization_correction,...
        Nx, Ny, Nz, dx, dy, dz, bx, by, bz, x, y, z_det, xx, yy, size_x, NSinos, NSlices, zmax, attenuation_correction, pseudot, det_per_ring, ...
        TOF, sigma_x, TOFCenter, dec, nCores, ind_size, block1, blocks, index, iij, jji, kkj, LL, N, summa, lor_a, xy_index, z_index, ...
        x_center, y_center, z_center, bmin, bmax, Vmax, V, lor_orth, options.gaussK,is_transposed, scatter_input, normalization, SinDelayed, double(n_meas(end)));
    
    if store_matrix
        varargout{1} = A;
    else
        if size(A,2) ~= size(f,1)
            if randoms_correction
                varargout{1} = A' * f + SinDelayed;
            else
                varargout{1} = A' * f;
            end
        else
            if randoms_correction
                varargout{1} = A * f + SinDelayed;
            else
                varargout{1} = A * f;
            end
        end
    end
    if nargout >= 2
        varargout{2} = options;
    end
elseif options.implementation == 4
    f = double(f);
    no_norm = uint8(1);
    
    epps = 1e-8;
    list_mode_format = options.listmode;
    if options.rings > 1
        dc_z = z_det(2,1) - z_det(1,1);
    else
        dc_z = options.cr_pz;
    end
    if use_raw_data
        xy_index = uint32(0);
        z_index = uint32(0);
        TOFSize = int64(size(LL,1));
    else
        LL = uint16(0);
        TOFSize = int64(numel(xy_index));
    end
    uu = ones(n_meas(end) * options.TOF_bins,1);
    
    if options.projector_type == 1
        if exist('OCTAVE_VERSION','builtin') == 0
            [~, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                normalization, SinDelayed, n_meas(end), attenuation_correction, normalization_correction, randoms_correction,...
                options.scatter, scatter_input, options.global_correction_factor, lor_a, xy_index, z_index, NSinos, LL, pseudot, det_per_ring, ...
                TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, nCores, ...
                (use_raw_data), uint32(1), list_mode_format, epps, uu, f, uint32(options.projector_type), no_norm, options.precompute_lor, uint8(1), ...
                options.n_rays_transaxial, options.n_rays_axial, dc_z);
        elseif exist('OCTAVE_VERSION','builtin') == 5
            [~, rhs] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                normalization, SinDelayed, n_meas(end), attenuation_correction, normalization_correction, randoms_correction,...
                options.scatter, scatter_input, options.global_correction_factor, lor_a, xy_index, z_index, NSinos, LL, pseudot, det_per_ring, ...
                TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, nCores, ...
                (use_raw_data), uint32(1), list_mode_format, epps, uu, f, uint32(options.projector_type), no_norm, options.precompute_lor, uint8(1), ...
                options.n_rays_transaxial, options.n_rays_axial, dc_z);
        end
    elseif options.projector_type == 2
        if exist('OCTAVE_VERSION','builtin') == 0
            [~, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                normalization, SinDelayed, n_meas(end), attenuation_correction, normalization_correction, randoms_correction,...
                options.scatter, scatter_input, options.global_correction_factor, lor_a, xy_index, z_index, NSinos, LL, pseudot, det_per_ring, ...
                TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, nCores, ...
                (use_raw_data), uint32(1), list_mode_format, epps, uu, f, uint32(options.projector_type), no_norm, options.precompute_lor, uint8(1), ...
                options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z);
        elseif exist('OCTAVE_VERSION','builtin') == 5
            [~, rhs] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                normalization, SinDelayed, n_meas(end), attenuation_correction, normalization_correction, randoms_correction,...
                options.scatter, scatter_input, options.global_correction_factor, lor_a, xy_index, z_index, NSinos, LL, pseudot, det_per_ring, ...
                TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, nCores, ...
                (use_raw_data), uint32(1), list_mode_format, epps, uu, f, uint32(options.projector_type), no_norm, options.precompute_lor, uint8(1), ...
                options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z);
        end
    elseif options.projector_type == 3
        if exist('OCTAVE_VERSION','builtin') == 0
            [~, rhs] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                normalization, SinDelayed, n_meas(end), attenuation_correction, normalization_correction, randoms_correction,...
                options.scatter, scatter_input, options.global_correction_factor, lor_a, xy_index, z_index, NSinos, LL, pseudot, det_per_ring, ...
                TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, nCores, ...
                (use_raw_data), uint32(1), list_mode_format, epps, uu, f, uint32(options.projector_type), no_norm, options.precompute_lor, uint8(1), ...
                x_center, y_center, z_center, bmin, bmax, Vmax, V);
        elseif exist('OCTAVE_VERSION','builtin') == 5
            [~, rhs] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                normalization, SinDelayed, n_meas(end), attenuation_correction, normalization_correction, randoms_correction,...
                options.scatter, scatter_input, options.global_correction_factor, lor_a, xy_index, z_index, NSinos, LL, pseudot, det_per_ring, ...
                TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, nCores, ...
                (use_raw_data), uint32(1), list_mode_format, epps, uu, f, uint32(options.projector_type), no_norm, options.precompute_lor, uint8(1), ...
                x_center, y_center, z_center, bmin, bmax, Vmax, V);
        end
    else
        error('Unsupported projector')
    end
    
    varargout{1} = rhs;
    if nargout >= 2
        varargout{2} = options;
    end
elseif options.implementation == 3 || options.implementation == 2
    %     options = double_to_single(options);
    
    f = single(f);
    if use_raw_data
        xy_index = uint32(0);
        z_index = uint32(0);
        TOFSize = int64(size(LL,1));
    else
        if isempty(pseudot)
            pseudot = int32(100000);
        end
        LL = uint16(0);
        TOFSize = int64(size(xy_index,1));
    end
    %     if ~iscell(SinM)
    %         SinM = {single(SinM)};
    %     end
    SinM = {ones(n_meas(end) * options.TOF_bins,1,'single')};
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
%     header_directory = strcat('-I "', header_directory);
%     header_directory = strcat(header_directory,'"');
    
    if options.verbose
        tStart = tic;
    end
    [output] = OpenCL_matrixfree_multi_gpu( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end), ... 15
        single(NSlices), size_x, zmax, options.verbose, LL, pseudot, det_per_ring, TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), int32(dec), uint32(options.use_device), filename, uint8(use_raw_data), ...25
        single(options.cpu_to_gpu_factor), uint32(0), header_directory, options.vaimennus, normalization, n_meas(end), uint32(attenuation_correction), ...32
        uint32(normalization_correction), lor_a, xy_index, z_index, tube_width_xy, crystal_size_z, x_center, y_center, z_center, SinDelayed, randoms, ...43
        uint32(options.projector_type), options.precompute_lor, n_rays, n_rays3D, dc_z, SinM, logical(options.use_64bit_atomics), f, uint8(true), ...53
        options.global_correction_factor, bmin, bmax, Vmax, V, options.use_psf, options);
    if options.verbose
        toc(tStart)
    end
    if isa(output{1},'int64')
        varargout{1} = single(output{1}) / single(100000000000);
    else
        varargout{1} = output{1};
    end
    if nargout >= 2
        varargout{2} = options;
    end
else
    error('Only implementations 1 and 3 are available in forward/backward projection')
end
if options.verbose
    disp('Forward projection done')
end