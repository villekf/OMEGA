function varargout = forward_project(options, index, n_meas, f)
%FORWARD_PROJECT Calculates the forward projection
% Examples:
%   fp = forward_project(options, index, n_meas, f)
%   [fp,norm] = forward_project(options, index, n_meas, f)
% INPUTS:
%   options = Contains the machine and (optionally) sinogram specific
%   information (e.g. detectors per ring, image (matrix) size, FOV size)
%   index = The indices (LORs) used to compute the system matrix (you can
%   use index_maker to produce the indices)
%   n_meas = Number of measurements used
%   f = The current estimate
%
% OUTPUTS:
%   fp = The forward projection (fp = A * f)
%   norm = The (optional) normalization constant (norm = sum(A,1))
%   options = If randoms, scatter or normalization correction is used, then
%   they are stored in the options variable at sub-iteration 1
%
% See also index_maker, backproject

if nargout > 2
    error('Too many output arguments')
end
if nargout == 0
    error('Too few output arguments')
end

folder = fileparts(which('reconstructions_main.m'));
folder = strrep(folder, 'source','mat-files/');
folder = strrep(folder, '\','/');

f = f(:);

Nx = options.Nx;
Ny = options.Ny;
Nz = options.Nz;
NSlices = uint32(Nz);
attenuation_correction = options.attenuation_correction;
FOVax = options.FOVa_x;
FOVay = options.FOVa_y;
rings = options.rings;
det_per_ring = options.det_per_ring;
machine_name = options.machine_name;
Nang = options.Nang;
Ndist = options.Ndist;
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
TotSinos = int32(options.TotSinos);

if numel(f) ~= Nx*Ny*Nz
    error('Estimate has different amount of elements than the image size')
end


[x, y, z] = get_coordinates(options, blocks);

[normalization_correction, randoms_correction] = set_up_corrections(options, blocks);
% if attenuation_correction
%     data = load(attenuation_datafile);
%     variables = fields(data);
%     vaimennus = double(data.(variables{1}));
%     if size(vaimennus,1) ~= Nx || size(vaimennus,2) ~= Ny || size(vaimennus,3) ~= Nz
%         if size(vaimennus,1) ~= Nx*Ny*Nz
%             error("Error: Attenuation data is of different size than the reconstructed image")
%         end
%     end
%     if size(vaimennus,2) == 1
%         vaimennus = vaimennus(:,:,2*block1+1:2*blocks+1);
%     else
%         vaimennus = vaimennus(2*block1+1:(2*blocks+1)*Nx*Ny);
%     end
%     vaimennus = vaimennus(:) / 10;
%     if options.reconstruction_method == 2 || options.reconstruction_method == 3 || options.reconstruction_method == 5
%         vaimennus = single(vaimennus);
%     end
%     clear data
% else
%     if options.reconstruction_method == 2 || options.reconstruction_method == 3 || options.reconstruction_method == 5
%         vaimennus = single(0);
%     else
%         vaimennus = 0;
%     end
% end
% 
% if options.normalization_correction && options.use_user_normalization && options.corrections_during_reconstruction
%     normalization_correction = true;
%     if ~isfield(options,'normalization')
%         [options.file, options.fpath] = uigetfile({'*.nrm;*.mat'},'Select normalization datafile');
%         if isequal(options.file, 0)
%             error('No file was selected')
%         end
%         if any(strfind(options.file, '.nrm'))
%             fid = fopen(options.file);
%             normalization = fread(fid, inf, 'single=>single',0,'l');
%             fclose(fid);
%         else
%             data = load(options.file);
%             variables = fields(data);
%             normalization = data.(variables{1});
%             clear data
%             if numel(normalization) ~= options.Ndist*options.Nang*options.TotSinos
%                 error('Size mismatch between the current data and the normalization data file')
%             end
%         end
%         normalization = normalization(:);
%         options.normalization = normalization;
%     else
%         normalization = options.normalization;
%     end
%     if (options.reconstruction_method == 1 || options.reconstruction_method == 3)
%         normalization = double(normalization);
%         normalization = 1 ./ normalization;
%     else
%         normalization = single(1) ./ normalization;
%     end
%     normalization = normalization(index);
% else
%     normalization_correction = false;
%     if options.reconstruction_method == 2 || options.reconstruction_method == 3 || options.reconstruction_method == 5
%         normalization = single(0);
%     else
%         normalization = 0;
%     end
% end
% 
% if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
%     if options.randoms_correction && ~isfield(options,'SinDelayed')
%         randoms_correction = true;
%         [randoms_file, r_fpath] = uigetfile('*.mat','Select randoms correction data');
%         
%         FileName = fullfile(r_fpath, randoms_file);
%         storedStructure = load(FileName);
%         variables = fields(storedStructure);
%         
%         SinDelayed = storedStructure.(variables{1});
%         options.SinDelayed = SinDelayed;
%     else
%         SinDelayed = options.SinDelayed;
%     end
%     if options.scatter_correction && ~isfield(options,'ScatterC')
%         [scatter_file, s_fpath] = uigetfile('*.mat','Select scatter correction data');
%         
%         FileName = fullfile(s_fpath, scatter_file);
%         storedStructure = load(FileName);
%         variables = fields(storedStructure);
%         
%         ScatterC = storedStructure.(variables{1});
%         options.ScatterC = ScatterC;
%     else
%         ScatterC = options.ScatterC;
%     end
%     if iscell(SinDelayed)
%         if options.scatter_correction
%             if sum(size(ScatterC)) > 1 && ~iscell(ScatterC)
%                 if size(ScatterC,1) ~= size(options.Nang)
%                     ScatterC = permute(ScatterC,[2 1 3]);
%                 end
%             elseif iscell(ScatterC)
%                 if sum(size(ScatterC{1})) > 1
%                     if size(ScatterC{1},1) ~= size(options.Nang)
%                         if length(ScatterC) > 1
%                             ScatterC = permute(ScatterC{1},[2 1 3]);
%                         else
%                             ScatterC = permute(ScatterC{1},[2 1 3]);
%                         end
%                     end
%                 end
%             end
%             SinDelayed = SinDelayed{1} + ScatterC(:);
%         end
%     else
%         SinDelayed = SinDelayed(:);
%         if use_raw_data == false && NSinos ~= TotSinos
%             SinDelayed = SinDelayed(1:NSinos*Ndist*Nang);
%         end
%         if options.scatter_correction
%             if sum(size(ScatterC)) > 1 && ~iscell(ScatterC)
%                 if size(ScatterC,1) ~= size(options.Nang)
%                     ScatterC = permute(ScatterC,[2 1 3]);
%                 end
%             elseif iscell(ScatterC)
%                 if sum(size(ScatterC{1})) > 1
%                     if size(ScatterC{1},1) ~= size(options.Nang)
%                         ScatterC = permute(ScatterC{1},[2 1 3]);
%                     end
%                 end
%             end
%             SinDelayed = SinDelayed + ScatterC(:);
%         end
%     end
%     SinDelayed = SinDelayed(index);
% else
%     randoms_correction = false;
%     if options.reconstruction_method == 2 || options.reconstruction_method == 3 || options.reconstruction_method == 5
%         SinDelayed = {single(0)};
%     else
%         SinDelayed = 0;
%     end
% end

if min(min(z)) == 0
    z = z + (axial_fow - max(max(z)))/2;
end

if options.reconstruction_method == 2 || options.reconstruction_method == 3 || options.reconstruction_method == 5
    x=single(x);
    y=single(y);
    z_det = single(z);
else
    x=double(x);
    y=double(y);
    z_det = double(z);
end
clear z


if ~options.precompute_lor
    lor_a = uint16(0);
end


size_x = int32(size(x,1));

if (options.precompute_lor  || options.reconstruction_method == 5 || options.reconstruction_method == 2 || options.reconstruction_method == 3)
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

% for the precomputed version, index vectors are needed
if use_raw_data == false && options.precompute_lor
    
    
    lor_file = [folder machine_name '_lor_pixel_count_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_sino_' num2str(Ndist) 'x' ...
        num2str(Nang) 'x' num2str(TotSinos) '.mat'];
    
    if exist(lor_file, 'file') == 2
        if options.reconstruction_method == 1 || options.reconstruction_method == 3
            variableInfo = who('-file', lor_file);
            if any(ismember('lor', variableInfo))
                load(lor_file,'lor')
            else
                lor_pixel_count_prepass(options);
                load(lor_file,'lor')
            end
            if options.projector_type == 2 && options.reconstruction_method == 1
                if any(ismember('lor_orth', variableInfo))
                    load(lor_file,'lor_orth')
                else
                    lor_pixel_count_prepass(options);
                    load(lor_file,'lor_orth')
                end
                load(lor_file,'crystal_size_z')
                load(lor_file,'crystal_size_xy')
                if options.tube_width_z == 0
                    if crystal_size_xy ~= options.tube_width_xy
                        lor_pixel_count_prepass(options);
                        load(lor_file,'lor_orth')
                    end
                    if length(lor_orth) > options.Nang*options.Ndist*options.TotSinos
                        lor_orth = lor_orth(1:options.Nang*options.Ndist*options.TotSinos);
                    end
                elseif options.tube_width_z > 0
                    if crystal_size_z ~= options.tube_width_z
                        lor_pixel_count_prepass(options);
                        load(lor_file,'lor_orth')
                    end
                    if length(lor_orth) > options.Nang*options.Ndist*options.TotSinos
                        lor_orth = lor_orth(options.Nang*options.Ndist*options.TotSinos+1:end);
                    elseif length(lor_orth) == options.Nang*options.Ndist*options.TotSinos
                        lor_pixel_count_prepass(options);
                        load(lor_file,'lor_orth')
                        lor_orth = lor_orth(options.Nang*options.Ndist*options.TotSinos+1:end);
                    end
                end
            end
        else
            variableInfo = who('-file', lor_file);
            if any(ismember('lor_opencl', variableInfo))
                load(lor_file,'lor_opencl')
            else
                lor_pixel_count_prepass(options);
                load(lor_file,'lor_opencl')
            end
            lor = lor_opencl;
            clear lor_opencl
        end
    else
        lor_pixel_count_prepass(options);
        if options.reconstruction_method == 1 || options.reconstruction_method == 3
            load(lor_file,'lor')
            if options.projector_type == 2 && options.reconstruction_method == 1
                load(lor_file,'lor_orth')
                if options.tube_width_z == 0
                    if length(lor_orth) > options.Nang*options.Ndist*options.TotSinos
                        lor_orth = lor_orth(1:options.Nang*options.Ndist*options.TotSinos);
                    end
                elseif options.tube_width_z > 0
                    if length(lor_orth) > options.Nang*options.Ndist*options.TotSinos
                        lor_orth = lor_orth(options.Nang*options.Ndist*options.TotSinos+1:end);
                    end
                end
            end
        else
            load(lor_file,'lor_opencl')
            lor = lor_opencl;
            clear lor_opencl
        end
    end
    lor_a = (lor(index));
    clear lor
    if options.projector_type == 2 && options.reconstruction_method == 1
        lor_orth = (lor_orth(index));
    end
    if options.normalization_correction && options.corrections_during_reconstruction
        normalization = normalization(index);
    end
    if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
        SinDelayed = SinDelayed(index);
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
    
    z_index = uint16(I(:,1));
    if verLessThan('matlab','8.5')
        z_index = repeat_elem(z_index, size_x);
    else
        z_index = repelem(z_index, size_x);
    end
    z_index = z_index(index);
    apu = z_index > NSinos;
    z_index = z_index - 1;
    
    xy_index = xy_index(index);
    xy_index(apu) = xy_index(apu) + uint32(size_x);
    xy_index = xy_index - 1;
    
    summa = int64(sum(lor_a(n_meas(1)+1:n_meas(2))));
    
    
    clear discard I yt xt xy_index2 index apu
elseif use_raw_data && options.precompute_lor
    
    lor_file = [folder machine_name '_detector_locations_' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '_raw.mat'];
    if exist(lor_file, 'file') == 2
        variableInfo = who('-file', lor_file);
        if options.reconstruction_method == 1 || options.reconstruction_method == 3
            if ismember('lor', variableInfo)
                load(lor_file,'lor')
            else
                lor_pixel_count_prepass(options);
                load(lor_file,'lor')
            end
            if options.projector_type == 2 && options.reconstruction_method == 1
                if any(ismember('lor_orth', variableInfo))
                    load(lor_file,'lor_orth')
                else
                    lor_pixel_count_prepass(options);
                    load(lor_file,'lor_orth')
                end
                load(lor_file,'crystal_size_z')
                load(lor_file,'crystal_size_xy')
                if options.tube_width_z == 0
                    if crystal_size_xy ~= options.tube_width_xy
                        lor_pixel_count_prepass(options);
                        load(lor_file,'lor_orth')
                    end
                    if length(lor_orth) > options.Nang*options.Ndist*options.TotSinos
                        lor_orth = lor_orth(1:options.Nang*options.Ndist*options.TotSinos);
                    end
                elseif options.tube_width_z > 0
                    if crystal_size_z ~= options.tube_width_z
                        lor_pixel_count_prepass(options);
                        load(lor_file,'lor_orth')
                    end
                    if length(lor_orth) > options.Nang*options.Ndist*options.TotSinos
                        lor_orth = lor_orth(options.Nang*options.Ndist*options.TotSinos+1:end);
                    elseif length(lor_orth) == options.Nang*options.Ndist*options.TotSinos
                        lor_pixel_count_prepass(options);
                        load(lor_file,'lor_orth')
                        lor_orth = lor_orth(options.Nang*options.Ndist*options.TotSinos+1:end);
                    end
                end
            end
        else
            if ismember('lor_opencl', variableInfo)
                load(lor_file,'lor_opencl')
            else
                lor_pixel_count_prepass(options);
                load(lor_file,'lor_opencl')
            end
            lor = lor_opencl;
            clear lor_opencl
        end
    else
        lor_pixel_count_prepass(options);
        if options.reconstruction_method == 1 || options.reconstruction_method == 3
            load(lor_file,'lor')
            if options.projector_type == 2 && options.reconstruction_method == 1
                load(lor_file,'lor_orth')
                if options.tube_width_z == 0
                    if length(lor_orth) > options.Nang*options.Ndist*options.TotSinos
                        lor_orth = lor_orth(1:options.Nang*options.Ndist*options.TotSinos);
                    end
                elseif options.tube_width_z > 0
                    if length(lor_orth) > options.Nang*options.Ndist*options.TotSinos
                        lor_orth = lor_orth(options.Nang*options.Ndist*options.TotSinos+1:end);
                    end
                end
            end
        else
            load(lor_file,'lor_opencl')
            lor = lor_opencl;
            clear lor_opencl
        end
    end
    LL = form_detector_pairs_raw(rings, det_per_ring);
    LL = LL(index,:);
    lor_a = (lor(index));
    clear lor
    if options.projector_type == 2 && options.reconstruction_method == 1
        lor_orth = (lor_orth(index));
    end
    if options.normalization_correction && options.corrections_during_reconstruction
        normalization = normalization(index);
    end
    if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
        SinDelayed = SinDelayed(index);
    end
    kk = 1;
    apu = LL(n_meas(kk) + 1 : n_meas(kk + 1),:) - 1;
    apu2 = idivide(apu, uint16(det_per_ring));
    idx = apu2(:,1) == apu2(:,2);
    apu2 = apu(idx,:);
    ind = mod(apu2, uint16(det_per_ring)) + 1;
    yt = y(ind);
    y_i = yt(:,1) > yt(:,2);
    apu2(y_i,:) = fliplr(apu2(y_i,:));
    apu(idx,:) = apu2;
    LL(n_meas(kk) + 1 : n_meas(kk + 1),:) = apu + 1;
    summa = int64(sum(lor_a(n_meas(kk)+1:n_meas(kk+1))));
    
    clear apu apu2 idx ind yt y_i index
    
    LL = LL';
    LL = LL(:);    
elseif use_raw_data == false && ~options.precompute_lor && options.reconstruction_method > 1
        
    if options.normalization_correction && options.corrections_during_reconstruction
        normalization = normalization(index);
    end
    if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
        SinDelayed = SinDelayed(index);
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
    
    z_index = uint16(I(:,1));
    if verLessThan('matlab','8.5')
        z_index = repeat_elem(z_index, size_x);
    else
        z_index = repelem(z_index, size_x);
    end
    z_index = z_index(index);
    apu = z_index > NSinos;
    z_index = z_index - 1;
    
    xy_index = xy_index(index);
    xy_index(apu) = xy_index(apu) + uint32(size_x);
    xy_index = xy_index - 1;
    
    
    clear I yt xt xy_index2 index apu
elseif use_raw_data && ~options.precompute_lor
        
    if ~exist('LL','var')
        LL = form_detector_pairs_raw(rings, det_per_ring);
    end
    LL = LL(index,:);
    if options.normalization_correction && options.corrections_during_reconstruction
        normalization = normalization(index);
    end
    if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction
        SinDelayed = SinDelayed(index);
    end
    clear lor
    
    kk = 1;
    apu = LL(n_meas(kk) + 1 : n_meas(kk + 1),:) - 1;
    apu2 = idivide(apu, uint16(det_per_ring));
    idx = apu2(:,1) == apu2(:,2);
    apu2 = apu(idx,:);
    ind = mod(apu2, uint16(det_per_ring)) + 1;
    yt = y(ind);
    y_i = yt(:,1) > yt(:,2);
    apu2(y_i,:) = fliplr(apu2(y_i,:));
    apu(idx,:) = apu2;
    LL(n_meas(kk) + 1 : n_meas(kk + 1),:) = apu + 1;
    
    
    clear apu apu2 idx ind yt y_i index discard
    
    LL = LL';
    LL = LL(:);
end


% Pixels
etaisyys_x=(R-FOVax)/2;
etaisyys_y=(R-FOVay)/2;
if options.reconstruction_method == 2 || options.reconstruction_method == 3 || options.reconstruction_method == 5
    zz = linspace(single(0), single(axial_fow), Nz + 1);
    xx = single(linspace(etaisyys_x, R - etaisyys_x, Nx + 1));
    yy = single(linspace(etaisyys_y, R - etaisyys_y, Ny + 1));
else
    zz = linspace(double(0), double(axial_fow), Nz + 1);
    xx = double(linspace(etaisyys_x, R - etaisyys_x, Nx + 1));
    yy = double(linspace(etaisyys_y, R - etaisyys_y, Ny + 1));
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


zmax = max(max(z_det));
if zmax==0
    if options.reconstruction_method == 2 || options.reconstruction_method == 3 || options.reconstruction_method == 5
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

if options.reconstruction_method == 1
    if options.precompute_lor == false
        if use_raw_data == false
            if options.projector_type == 1 || options.projector_type == 0
                [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                    zmax, vaimennus, normalization, SinDelayed, pituus(osa_iter), attenuation_correction, normalization_correction, ...
                    randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, uint16(0), pseudot, det_per_ring, options.verbose, ...
                    use_raw_data, uint32(2), ind_size, block1, blocks, index{osa_iter}, uint32(options.projector_type));
            elseif options.projector_type == 2
                [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                    zmax, vaimennus, normalization, SinDelayed, pituus(osa_iter), attenuation_correction, normalization_correction, ...
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
                    zmax, vaimennus, normalization, SinDelayed, uint32(0), attenuation_correction, normalization_correction, ...
                    randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, LL, pseudot, det_per_ring, options.verbose, ...
                    use_raw_data, uint32(2), ind_size, block1, blocks, uint32(0), uint32(options.projector_type));
            elseif options.projector_type == 2
                [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                    zmax, vaimennus, normalization, SinDelayed, uint32(0), attenuation_correction, normalization_correction, ...
                    randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, LL, pseudot, det_per_ring, options.verbose, ...
                    use_raw_data, uint32(2), ind_size, block1, blocks, uint32(0), uint32(options.projector_type), options.tube_width_xy, ...
                    x_center, y_center, z_center, options.tube_width_z, int32(options.accuracy_factor));
            else
                error('Unsupported projector type')
            end
        end
        lor = reshape(lor,[],2);
        if verLessThan('matlab','8.5')
            lor=repeat_elem(int32((lor(:,1))),lor(:,2));
        else
            lor=repelem(int32((lor(:,1))),lor(:,2));
        end
        
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
            size_x, zmax, vaimennus, normalization, SinDelayed, n_meas(end), attenuation_correction, normalization_correction,...
            randoms_correction, lor_a, xy_index, ...
            z_index, NSinos, LL, pseudot, det_per_ring, options.verbose, ...
            use_raw_data, uint32(0), lor2, summa, false, uint32(options.projector_type), ...
            options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, int32(options.accuracy_factor));
        clear lor2
    end
    
    
    if size(A,2) ~= size(f,1)
        if no_norm
            varargout{1} = A' * f;
        else
            varargout{1} = A' * f;
        end
    else
        if no_norm
            varargout{1} = A * f;
        else
            varargout{1} = A * f;
        end
    end
    varargout{2} = options;
elseif options.reconstruction_method == 3
    %     options = double_to_single(options);
    
    f = single(f);
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
    [varargout{1},~] = OpenCL_matrixfree_multi_gpu( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end), ...
        single(NSlices), size_x, zmax, options.verbose, LL, pseudot, det_per_ring, uint32(options.use_device), filename, uint8(use_raw_data), ...
        single(options.cpu_to_gpu_factor), uint32(0), header_directory, vaimennus, normalization, n_meas(end), uint32(attenuation_correction), ...
        uint32(normalization_correction), lor_a, xy_index, z_index, tube_width_xy, crystal_size_z, x_center, y_center, z_center, SinDelayed, randoms, ...
        uint32(options.projector_type), options.precompute_lor, int32(options.accuracy_factor), f, false);
    toc
    varargout{2} = options;
else
    error('Only implementations 1 and 3 are available in forward/backward projection')
end