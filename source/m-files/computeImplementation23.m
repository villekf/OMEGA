function [pz, fp, res] = computeImplementation23(options, x, z_det, pituus)
%COMPUTEIMPLEMENTATION23 Computes the image reconstruction phase of
%implementations 2 or 3
%   Utility function

if ~isfield(options,'compute_sensitivity_image')
    options.compute_sensitivity_image = false;
end
fp = [];
res = [];
if ismac && options.use_64bit_atomics
    options.use_64bit_atomics = false;
end
if options.use_32bit_atomics && options.use_64bit_atomics
    options.use_64bit_atomics = false;
end
if options.use_64bit_atomics && (options.use_CPU || options.use_CUDA || options.projector_type == 6)
    options.use_64bit_atomics = false;
end
% if options.use_32bit_atomics && (options.use_CPU || options.use_CUDA || options.projector_type == 6)
%     options.use_32bit_atomics = false;
% end
options.use_device = uint32(options.use_device);
if options.use_raw_data
    TOFSize = int64(size(options.LL,1));
else
    TOFSize = int64(size(options.xy_index,1));
end
if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction ...
        && ~options.reconstruct_trues && ~options.reconstruct_scatter
    randoms = uint32(1);
else
    randoms = uint32(0);
end
n_rays = uint16(options.n_rays_transaxial);
n_rays3D = uint16(options.n_rays_axial);
if ~isfield(options, 'vaimennus')
    options.vaimennus = single(0);
end
if numel(options.partitions) > 1
    partitions = numel(options.partitions);
elseif isempty(options.partitions)
    partitions = 1;
else
    partitions = options.partitions;
end

% Dynamic SPECT: the C++ side reads the detector-head index of timestep tt / subset kk
% from DetectorVector at the offset pituus[kk + tt * subsets], i.e. it expects the vector
% concatenated over all Nt frames (as the Python interface does). OMEGA_error_check,
% setMissingValues and parseInputData only handle a single frame's worth
% (nProjections values, already permuted into subset order), which would make every
% timestep after the first read past the end of the array. The head layout is identical
% for every frame, so replicate the (already subset-ordered) frame Nt times.
if options.SPECT && partitions > 1 && isfield(options, 'DetectorVector') ...
        && ~isempty(options.DetectorVector) && numel(options.DetectorVector) == options.nProjections
    options.DetectorVector = repmat(uint32(options.DetectorVector(:)), partitions, 1);
end

% Each cell here holds one timestep's data, which may still be a multi-dimensional
% (Ndist x Nang x NSinos(x TOF_bins)) array rather than a column vector: parseInputData
% only columnizes per-frame cells when subsets > 1 (and subset_type > 0), so this branch
% -- reached whenever largeDim skipped reconstructions_main.m's own cell2mat, or whenever
% subsets == 1 -- must columnize each frame itself before cell2mat. Otherwise cell2mat
% concatenates the 1xNt cell along dimension 2 instead of stacking whole frames end to
% end, which interleaves timesteps instead of keeping each one contiguous (the C++ side
% offsets into each array by that array's per-timestep element count times the timestep
% index, i.e. expects frame 1 fully before frame 2).
if iscell(options.SinM)
    options.SinM = cell2mat(cellfun(@(x) x(:), options.SinM, 'UniformOutput', false));
end
if options.attenuation_correction
    if iscell(options.vaimennus)
        options.vaimennus = cell2mat(cellfun(@(x) x(:), options.vaimennus, 'UniformOutput', false));
    end
end
if options.randoms_correction
    if iscell(options.SinDelayed)
        options.SinDelayed = cell2mat(cellfun(@(x) x(:), options.SinDelayed, 'UniformOutput', false));
    end
end
options.currentSubset = 0;
options.currentTimestep = 0;

if options.orthAxial
    crystal_size_z = (options.tube_width_z);
else
    crystal_size_z = (options.tube_width_xy);
end
if options.projector_type == 1 || options.projector_type == 11 ...
        || options.projector_type == 2 || options.projector_type == 3 || options.projector_type == 22 || options.projector_type == 33 ...
        || options.projector_type == 13 || options.projector_type == 23 || options.projector_type == 12 || options.projector_type == 21 ...
        || options.projector_type == 32 || options.projector_type == 31
    kernel_file = 'projectorType123.cl';
    kernel_path = which(kernel_file);
    kernel_path = strrep(kernel_path, '\', '/');
    kernel_path = strrep(kernel_path, '.cl', '');
    header_directory = strrep(kernel_path,'projectorType123','');
elseif options.projector_type == 4 || options.projector_type == 41 || options.projector_type == 14 || options.projector_type == 45 ...
        || options.projector_type == 42 || options.projector_type == 43 || options.projector_type == 24 || options.projector_type == 34
    kernel_file = 'projectorType4.cl';
    kernel_path = which(kernel_file);
    kernel_path = strrep(kernel_path, '\', '/');
    kernel_path = strrep(kernel_path, '.cl', '');
    header_directory = strrep(kernel_path,'projectorType4','');
elseif options.projector_type == 5 || options.projector_type == 51 || options.projector_type == 15 || options.projector_type == 54 ...
         || options.projector_type == 52  || options.projector_type == 53  || options.projector_type == 35  || options.projector_type == 25
    kernel_file = 'projectorType5.cl';
    kernel_path = which(kernel_file);
    kernel_path = strrep(kernel_path, '\', '/');
    kernel_path = strrep(kernel_path, '.cl', '');
    header_directory = strrep(kernel_path,'projectorType5','');
elseif options.projector_type == 6
    kernel_file = 'projectorType123.cl';
    kernel_path = which(kernel_file);
    kernel_path = strrep(kernel_path, '\', '/');
    kernel_path = strrep(kernel_path, '.cl', '');
    header_directory = strrep(kernel_path,'projectorType123','');
elseif options.projector_type == 7 || options.projector_type == 71 || options.projector_type == 17 || options.projector_type == 74 ...
        || options.projector_type == 47 || options.projector_type == 75 || options.projector_type == 57
    kernel_file = 'projectorType7.cl';
    kernel_path = which(kernel_file);
    kernel_path = strrep(kernel_path, '\', '/');
    kernel_path = strrep(kernel_path, '.cl', '');
    header_directory = strrep(kernel_path,'projectorType7','');
else
    error('Invalid projector for OpenCL')
end
if options.listmode > 0 && options.compute_sensitivity_image && isfield(options,'xSens') && isfield(options,'zSens') && ~isempty(options.xSens) && ~isempty(options.zSens) && ~options.SPECT
    x = options.xSens;
    z_det = options.zSens;
	options.det_per_ring = numel(x) / 2;
	options.rings = numel(z_det);
elseif options.listmode > 0 && options.compute_sensitivity_image && ~options.useIndexBasedReconstruction
    options.use_raw_data = true;
    [x, ~, z_det, ~] = get_coordinates(options);
    options.use_raw_data = false;
end
if (~options.largeDim && ~(~options.loadTOF && ~isa(options.SinM,'single'))) || (options.largeDim && isa(options.SinM,'uint32')) || options.use_CPU
    options.SinM = single(options.SinM);
end
if options.implementation == 2
    if options.verbose
        tStart = tic;
    end
    if options.use_CUDA
        if isa(options.SinM,'uint16')
            [pz,fp, res] = CUDA_matrixfree_uint16( options.Nx, options.Ny, options.Nz, options.dx, options.dy, options.dz, options.bx, options.by, options.bz, ...
                z_det, x, options.nRowsD, options.verbose, options.LL, options.TOF, ... % 15
                TOFSize, options.sigma_x, options.TOFCenter, options.TOF_bins, options.use_device, options.use_raw_data, options.use_psf, header_directory, options.vaimennus, ... % 24
                options.normalization, pituus, options.attenuation_correction, options.normalization_correction, options.Niter, options.subsets, options.epps, options.xy_index, ...
                options.z_index, crystal_size_z, ... % 34
                options.x_center, options.y_center, options.z_center, options.SinDelayed, randoms, options.projector_type, n_rays, n_rays3D, ... % 42
                options, options.SinM, partitions, options.use_64bit_atomics, options.bmin, options.bmax, options.Vmax, options.V, options.gaussK); % 21
        elseif isa(options.SinM,'uint8')
            [pz,fp, res] = CUDA_matrixfree_uint8( options.Nx, options.Ny, options.Nz, options.dx, options.dy, options.dz, options.bx, options.by, options.bz, ...
                z_det, x, options.nRowsD, options.verbose, options.LL, options.TOF, ... % 15
                TOFSize, options.sigma_x, options.TOFCenter, options.TOF_bins, options.use_device, options.use_raw_data, options.use_psf, header_directory, options.vaimennus, ... % 24
                options.normalization, pituus, options.attenuation_correction, options.normalization_correction, options.Niter, options.subsets, options.epps, options.xy_index, ...
                options.z_index, crystal_size_z, ... % 34
                options.x_center, options.y_center, options.z_center, options.SinDelayed, randoms, options.projector_type, n_rays, n_rays3D, ... % 42
                options, options.SinM, partitions, options.use_64bit_atomics, options.bmin, options.bmax, options.Vmax, options.V, options.gaussK); % 21
        else
            [pz,fp, res] = CUDA_matrixfree( options.Nx, options.Ny, options.Nz, options.dx, options.dy, options.dz, options.bx, options.by, options.bz, ...
                z_det, x, options.nRowsD, options.verbose, options.LL, options.TOF, ... % 15
                TOFSize, options.sigma_x, options.TOFCenter, options.TOF_bins, options.use_device, options.use_raw_data, options.use_psf, header_directory, options.vaimennus, ... % 24
                options.normalization, pituus, options.attenuation_correction, options.normalization_correction, options.Niter, options.subsets, options.epps, options.xy_index, ...
                options.z_index, crystal_size_z, ... % 34
                options.x_center, options.y_center, options.z_center, options.SinDelayed, randoms, options.projector_type, n_rays, n_rays3D, ... % 42
                options, options.SinM, partitions, options.use_64bit_atomics, options.bmin, options.bmax, options.Vmax, options.V, options.gaussK); % 21
        end
    elseif options.use_CPU
        [pz,fp, res] = CPU_matrixfree( options.Nx, options.Ny, options.Nz, options.dx, options.dy, options.dz, options.bx, options.by, options.bz, ...
            z_det, x, options.nRowsD, options.verbose, options.LL, options.TOF, ... % 15
            TOFSize, options.sigma_x, options.TOFCenter, options.TOF_bins, options.use_device, options.use_raw_data, options.use_psf, header_directory, options.vaimennus, ... % 24
            options.normalization, pituus, options.attenuation_correction, options.normalization_correction, options.Niter, options.subsets, options.epps, options.xy_index, ...
            options.z_index, crystal_size_z, ... % 34
            options.x_center, options.y_center, options.z_center, options.SinDelayed, randoms, options.projector_type, n_rays, n_rays3D, ... % 42
            options, options.SinM, partitions, options.use_64bit_atomics, options.bmin, options.bmax, options.Vmax, options.V, options.gaussK); % 21
    else
        if isa(options.SinM,'uint16')
            [pz,fp, res] = OpenCL_matrixfree_uint16( options.Nx, options.Ny, options.Nz, options.dx, options.dy, options.dz, options.bx, options.by, options.bz, ...
                z_det, x, options.nRowsD, options.verbose, options.LL, options.TOF, ... % 15
                TOFSize, options.sigma_x, options.TOFCenter, options.TOF_bins, options.use_device, options.use_raw_data, options.use_psf, header_directory, options.vaimennus, ... % 24
                options.normalization, pituus, options.attenuation_correction, options.normalization_correction, options.Niter, options.subsets, options.epps, options.xy_index, ...
                options.z_index, crystal_size_z, ... % 34
                options.x_center, options.y_center, options.z_center, options.SinDelayed, randoms, options.projector_type, n_rays, n_rays3D, ... % 42
                options, options.SinM, partitions, options.use_64bit_atomics, options.bmin, options.bmax, options.Vmax, options.V, options.gaussK); % 21
        elseif isa(options.SinM,'uint8')
            [pz,fp, res] = OpenCL_matrixfree_uint8( options.Nx, options.Ny, options.Nz, options.dx, options.dy, options.dz, options.bx, options.by, options.bz, ...
                z_det, x, options.nRowsD, options.verbose, options.LL, options.TOF, ... % 15
                TOFSize, options.sigma_x, options.TOFCenter, options.TOF_bins, options.use_device, options.use_raw_data, options.use_psf, header_directory, options.vaimennus, ... % 24
                options.normalization, pituus, options.attenuation_correction, options.normalization_correction, options.Niter, options.subsets, options.epps, options.xy_index, ...
                options.z_index, crystal_size_z, ... % 34
                options.x_center, options.y_center, options.z_center, options.SinDelayed, randoms, options.projector_type, n_rays, n_rays3D, ... % 42
                options, options.SinM, partitions, options.use_64bit_atomics, options.bmin, options.bmax, options.Vmax, options.V, options.gaussK); % 21
        else
            [pz,fp, res] = OpenCL_matrixfree( options.Nx, options.Ny, options.Nz, options.dx, options.dy, options.dz, options.bx, options.by, options.bz, ...
                z_det, x, options.nRowsD, options.verbose, options.LL, options.TOF, ... % 15
                TOFSize, options.sigma_x, options.TOFCenter, options.TOF_bins, options.use_device, options.use_raw_data, options.use_psf, header_directory, options.vaimennus, ... % 24
                options.normalization, pituus, options.attenuation_correction, options.normalization_correction, options.Niter, options.subsets, options.epps, options.xy_index, ...
                options.z_index, crystal_size_z, ... % 34
                options.x_center, options.y_center, options.z_center, options.SinDelayed, randoms, options.projector_type, n_rays, n_rays3D, ... % 42
                options, options.SinM, partitions, options.use_64bit_atomics, options.bmin, options.bmax, options.Vmax, options.V, options.gaussK); % 21
        end
    end
    if options.verbose
        toc(tStart)
    end

    %     pz(end,:) = [];
else

    if options.verbose
        tStart = tic;
    end
    try
        [pz] = OpenCL_matrixfree_multi_gpu( options.Nx, options.Ny, options.Nz, options.dx, options.dy, options.dz, options.bx, options.by, options.bz, ...
            z_det, x, options.nRowsD, options.verbose, options.LL, options.TOF, ... % 15
            TOFSize, options.sigma_x, options.TOFCenter, options.TOF_bins, options.platform, options.use_raw_data, options.use_psf, header_directory, options.vaimennus, ... % 24
            options.normalization, pituus, options.attenuation_correction, options.normalization_correction, options.Niter, options.subsets, options.epps, options.xy_index, ...
            options.z_index, crystal_size_z, ... % 34
            options.x_center, options.y_center, options.z_center, options.SinDelayed, randoms, options.projector_type, n_rays, n_rays3D, ... % 42
            options, options.SinM, partitions, options.use_64bit_atomics, options.bmin, options.bmax, options.Vmax, options.V, options.gaussK, 0, 0); % 51
        if options.verbose
            toc(tStart)
        end
    catch ME
        error(['Error when running MEX-file: ' ME.message])
    end
    pz = reshape(pz, options.Nx(1), options.Ny(1), options.Nz(1));

    % if type == 0
    %     for ll = 1 : options.partitions
    %         apu = tz{2,ll};
    %         if options.save_iter
    %             apu(:,1) = options.x0;
    %             pz = reshape(apu, Nx, Ny, Nz, Niter + 1, options.partitions);
    %         else
    %             pz = reshape(apu, Nx, Ny, Nz, options.partitions);
    %         end
    %     end
    %     clear tz
    % end
end

end
