function [pz] = computeImplementation23(options, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx, NSinos, NSlices, size_x, zmax, ...
    LL, pseudot, det_per_ring, TOF, sigma_x, TOFCenter, dec, device, use_raw_data, normalization, pituus, attenuation_correction, ...
    normalization_correction, Niter, subsets, epps, lor_a, xy_index, z_index, x_center, y_center, z_center, SinDelayed, ...
    SinM, bmin, bmax, Vmax, V, gaussK, varargin)
%COMPUTEIMPLEMENTATION23 Computes the image reconstruction phase of
%implementations 2 or 3
%   Utility function

if ~isfield(options,'compute_sensitivity_image')
    options.compute_sensitivity_image = false;
end
if ~isfield(options,'listmode')
    options.listmode = false;
end
if ~isfield(options,'use_CUDA')
    options.use_CUDA = false;
end
if ~isfield(options,'use_32bit_atomics')
    options.use_32bit_atomics = false;
end
if isempty(varargin)
    type = 0;
else
    type = varargin{1};
end
if options.use_32bit_atomics && options.use_64bit_atomics
    options.use_32bit_atomics = false;
end
if type <= 0
    rekot = varargin{2};
    if numel(varargin) >= 3 && ~isempty(varargin{3})
        pz = varargin{3};
    else
        pz = {[]};
    end
end
if use_raw_data
    xy_index = uint32(0);
    z_index = uint16(0);
    TOFSize = int64(size(LL,1));
else
    if isempty(pseudot)
        pseudot = uint32(100000);
    end
    LL = uint16(0);
    TOFSize = int64(size(xy_index,1));
end
if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction ...
        && ~options.reconstruct_trues && ~options.reconstruct_scatter
    randoms = uint32(1);
else
    randoms = uint32(0);
end
n_rays = uint16(options.n_rays_transaxial);
n_rays3D = uint16(options.n_rays_axial);
if n_rays * n_rays3D > 1
    dc_z = single(z_det(2,1) - z_det(1,1));
else
    dc_z = single(options.cr_pz);
end
if ~isfield(options, 'vaimennus')
    options.vaimennus = single(0);
end


tube_width_xy = single(options.tube_width_xy);
crystal_size_z = single(options.tube_width_z);
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
if options.use_CUDA
    header_directory = strcat('-I"', header_directory);
    header_directory = strcat(header_directory,'"');
end
if options.implementation == 2
    joku = algorithms_char();
%     varList = recNames();
    varPrior = recNames(1);
    varMAP = recNames(2);
    varMAPML = recNames(7);
%     varMAPOS = recNames(8);
    varMLEM = recNames(5);
    varOS = recNames(6);
    options.nTot = uint32(sum(rekot));
    apu = false(numel(varPrior),1);
    for kk = 1 : numel(varPrior)
        if options.(varPrior{kk}) && ~strcmp(varPrior{kk},'custom')
            apu(kk) = true;
        end
    end
    varPrior = varPrior(apu);
    options.nPriors = uint32(numel(varPrior));
    apu = false(numel(varMAP),1);
    for kk = 1 : numel(varMAP)
        if options.(varMAP{kk})
            apu(kk) = true;
        end
    end
    varMAP = varMAP(apu);
    options.nMAP = uint32(numel(varMAP));
    apu = false(numel(varMAPML),1);
    for kk = 1 : numel(varMAPML)
        if options.(varMAPML{kk})
            apu(kk) = true;
        end
    end
    varMAPML = varMAPML(apu);
    options.nMAPML = uint32(numel(varMAPML));
    apu = false(numel(varMLEM),1);
    for kk = 1 : numel(varMLEM)
        if options.(varMLEM{kk})
            apu(kk) = true;
        end
    end
    varMLEM = varMLEM(apu);
    options.nMLEM = uint32(numel(varMLEM));
    apu = false(numel(varOS),1);
    for kk = 1 : numel(varOS)
        if options.(varOS{kk})
            apu(kk) = true;
        end
    end
    varOS = varOS(apu);
    options.nOS = uint32(numel(varOS));
    options.rekoList = uint32(find(rekot) - 1);
    % Add algorithms that use prior between iterations rather than
    % subiterations here
    options.mIt = cell(2,1);
    options.mIt{1} = int32(find(strcmp(varMAP,'BSREM')));
    options.mIt{2} = int32(find(strcmp(varMAP,'ROSEM_MAP')));
    for kk = 1 : numel(options.mIt)
        if isempty(options.mIt{kk})
            options.mIt{kk} = int32(-1);
        else
            options.mIt{kk} = int32(options.mIt{kk} - 1 - double(options.nMAPML));
        end
    end
    if options.nPriors > 0 && options.MAP
        options.varList = cell(options.nMAP * options.nPriors,1);
        uu = 1;
        for kk = 1 : options.nPriors
            for ll = 1 : options.nMAP
                options.varList{uu} = ['beta_' varPrior{kk} '_' varMAP{ll}];
                uu = uu + 1;
            end
        end
    end
end
if options.listmode && options.compute_sensitivity_image
    options.listmode = uint8(2);
    [xd, yd] = detector_coordinates(options);
    
    z_length = single(options.rings + 1 + sum(options.pseudot)) * options.cr_pz;
    zd = linspace(0, z_length, options.rings + 2 + sum(options.pseudot))';
    if sum(options.pseudot) > 0
        zd(options.pseudot) = [];
    end
    if min(zd(:)) == 0
        zd = zd + (options.axial_fov - (options.rings + sum(options.pseudot)) * options.cr_pz)/2 + options.cr_pz/2;
    end
    zd = single(zd(1:end-2));
    LL = form_detector_pairs_raw(options.rings, options.det_per_ring)';
    LL = LL(:);
    n_rekos = uint32(0);
    n_rekos_mlem = uint32(1);
    reko_type = uint8([]);
    reko_type_mlem = uint8(0);
    pituusD = int64([0;numel(LL)/2]);
    
    if options.implementation == 2
        if ~options.use_CUDA
            [pz] = OpenCL_matrixfree( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, zd, xd, yd, dy, yy(end), xx(end), NSinos, single(NSlices), size_x, max(zd(:)), NSinos, ...
                options.verbose, LL, pseudot, uint32(options.det_per_ring), TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), uint32(dec), device, uint8(true), ...
                filename, uint32(0), options.use_psf, header_directory, options.vaimennus, normalization, pituusD, uint32(attenuation_correction), ...
                uint32(normalization_correction), uint32(1), uint32(1), uint8(rekot), single(epps), lor_a, xy_index, z_index, any(n_rekos), tube_width_xy, ...
                crystal_size_z, x_center, y_center, z_center, SinDelayed, randoms, uint32(options.projector_type), options.precompute_lor, n_rays, n_rays3D, ...
                dc_z, options, SinM, uint32(1), logical(options.use_64bit_atomics), n_rekos, n_rekos_mlem, reko_type, reko_type_mlem, ...
                options.global_correction_factor, bmin, bmax, Vmax, V, gaussK);
        else
            header_directory = strrep(header_directory,'"','');
            [pz] = CUDA_matrixfree( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end), NSinos, single(NSlices), size_x, max(zd(:)), NSinos, ...
                options.verbose, LL, pseudot, uint32(options.det_per_ring), TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), uint32(dec), device, uint8(true), ...
                filename, uint32(0), options.use_psf, header_directory, options.vaimennus, normalization, pituus, uint32(attenuation_correction), ...
                uint32(normalization_correction), uint32(1), uint32(1), uint8(rekot), single(epps), lor_a, xy_index, z_index, any(n_rekos), tube_width_xy, ...
                crystal_size_z, x_center, y_center, z_center, SinDelayed, randoms, uint32(options.projector_type), options.precompute_lor, n_rays, n_rays3D, ...
                dc_z, options, SinM, uint32(1), logical(options.use_64bit_atomics), n_rekos, n_rekos_mlem, reko_type, reko_type_mlem, ...
                options.global_correction_factor, bmin, bmax, Vmax, V, gaussK);
        end
        options.listmode = uint8(1);
        options.Summ = pz{1}(:,:,:,end);
    else
        [tz] = OpenCL_matrixfree_multi_gpu( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, zd, xd, yd, dy, yy(end), xx(end), single(NSlices), size_x, max(zd(:)), options.verbose, ...
            LL, pseudot, uint32(options.det_per_ring), TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), uint32(dec), uint32(options.use_device), filename, uint8(true), single(options.cpu_to_gpu_factor), uint32(1), header_directory, ...
            options.vaimennus, normalization, pituusD, uint32(attenuation_correction), uint32(normalization_correction), lor_a, xy_index, z_index, tube_width_xy, ...
            crystal_size_z, x_center, y_center, z_center, SinDelayed, randoms, uint32(options.projector_type), options.precompute_lor, ...
            n_rays, n_rays3D, dc_z, SinM, logical(options.use_64bit_atomics), NSinos, uint16(NSinos), uint32(1), uint32(1), uint8(rekot), ...
            single(epps), uint32(1), options.OSEM, options.use_psf, options.global_correction_factor, bmin, bmax, Vmax, V, gaussK, options);
        options.Summ = tz{1}(:,end) / single(subsets);
        if (options.use_psf)
            options.Summ = computeConvolution(options.Summ, options, Nx, Ny, Nz, gaussK);
        end
    end
    LL = uint16(0);
end
if options.implementation == 2
    % n_rekos = uint32(sum(rekot(~contains(joku,'MLEM'))));
    n_rekos = uint32(sum(rekot(cellfun('isempty',strfind(joku,'MLEM')))));
    n_rekos_mlem = uint32(sum(rekot(~cellfun('isempty',strfind(joku,'MLEM')))));
    reko_type = zeros(length(rekot),1,'uint8');
    reko_type(~cellfun('isempty',strfind(joku,'MBSREM'))) = 1;
    reko_type(~cellfun('isempty',strfind(joku,'MRAMLA'))) = 1;
    reko_type(~cellfun('isempty',strfind(joku,'COSEM'))) = 2;
    reko_type(~cellfun('isempty',strfind(joku,'ACOSEM'))) = 3;
    reko_type(~cellfun('isempty',strfind(joku,'OSL-COSEM')) & options.OSL_COSEM == 1) = 2;
    reko_type(~cellfun('isempty',strfind(joku,'OSL-COSEM')) & options.OSL_COSEM == 2) = 3;
    ind = cellfun('isempty',strfind(joku,'MLEM'));
    reko_type = reko_type(rekot & ind);
    joku = joku(rekot & ind);
    if options.ECOSEM
        if options.COSEM && options.OSEM
            reko_type(~cellfun('isempty',strfind(joku,'ECOSEM'))) = [];
        elseif options.COSEM && ~options.OSEM
            reko_type(~cellfun('isempty',strfind(joku,'ECOSEM'))) = [];
            reko_type = [0;reko_type];
        elseif ~options.COSEM && ~options.OSEM
            reko_type = [0;reko_type];
        end
    end
    reko_type_mlem = zeros(n_rekos_mlem,1,'uint8');
    if type < 0
        rekot = [rekot; false; false; false; false];
        options.save_iter = false;
        options.use_psf = false;
    end
    % reko_type(contains(joku,'MBSREM')) = 1;
    % reko_type(contains(joku,'MRAMLA')) = 1;
    % reko_type(contains(joku,'COSEM')) = 2;
    % reko_type(contains(joku,'ACOSEM')) = 3;
    % reko_type(contains(joku,'OSL-COSEM') & options.COSEM_OSL == 1) = 2;
    % reko_type(contains(joku,'OSL-COSEM') & options.COSEM_OSL == 2) = 3;
    if options.verbose
        tStart = tic;
    end
    if ~options.use_CUDA
        [pz] = OpenCL_matrixfree( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end), NSinos, single(NSlices), size_x, zmax, NSinos, ...
            options.verbose, LL, pseudot, det_per_ring, TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), uint32(dec), device, uint8(use_raw_data), ...
            filename, uint32(0), options.use_psf, header_directory, options.vaimennus, normalization, pituus, uint32(attenuation_correction), ...
            uint32(normalization_correction), uint32(Niter), uint32(subsets), uint8(rekot), single(epps), lor_a, xy_index, z_index, any(n_rekos), tube_width_xy, ...
            crystal_size_z, x_center, y_center, z_center, SinDelayed, randoms, uint32(options.projector_type), options.precompute_lor, n_rays, n_rays3D, ...
            dc_z, options, SinM, uint32(options.partitions), logical(options.use_64bit_atomics), n_rekos, n_rekos_mlem, reko_type, reko_type_mlem, ...
            options.global_correction_factor, bmin, bmax, Vmax, V, gaussK);
    else
        header_directory = strrep(header_directory,'"','');
        [pz] = CUDA_matrixfree( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end), NSinos, single(NSlices), size_x, zmax, NSinos, ...
            options.verbose, LL, pseudot, det_per_ring, TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), uint32(dec), device, uint8(use_raw_data), ...
            filename, uint32(0), options.use_psf, header_directory, options.vaimennus, normalization, pituus, uint32(attenuation_correction), ...
            uint32(normalization_correction), uint32(Niter), uint32(subsets), uint8(rekot), single(epps), lor_a, xy_index, z_index, any(n_rekos), tube_width_xy, ...
            crystal_size_z, x_center, y_center, z_center, SinDelayed, randoms, uint32(options.projector_type), options.precompute_lor, n_rays, n_rays3D, ...
            dc_z, options, SinM, uint32(options.partitions), logical(options.use_64bit_atomics), n_rekos, n_rekos_mlem, reko_type, reko_type_mlem, ...
            options.global_correction_factor, bmin, bmax, Vmax, V, gaussK);
    end
    if options.verbose
        toc(tStart)
    end
    
    pz(end,:) = [];
else
    
    if options.verbose
        tStart = tic;
    end
    try
        if type <= 0
            [tz] = OpenCL_matrixfree_multi_gpu( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end), single(NSlices), size_x, zmax, options.verbose, ... 19
                LL, pseudot, det_per_ring, TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), uint32(dec), uint32(options.use_device), filename, uint8(use_raw_data), ... 31
                single(options.cpu_to_gpu_factor), uint32(1), header_directory, options.vaimennus, normalization, pituus, uint32(attenuation_correction), ... 38
                uint32(normalization_correction), lor_a, xy_index, z_index, tube_width_xy, crystal_size_z, x_center, y_center, z_center, SinDelayed, randoms, ...49
                uint32(options.projector_type), options.precompute_lor, n_rays, n_rays3D, dc_z, SinM, logical(options.use_64bit_atomics), NSinos, uint16(NSinos), uint32(Niter), ... 59
                uint32(subsets), uint8(rekot), single(epps), uint32(options.partitions), options.OSEM, options.use_psf, options.global_correction_factor, bmin, bmax, Vmax, ...
                V, gaussK, options);
        else
            [pz] = OpenCL_matrixfree_multi_gpu( kernel_path, Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy(end), xx(end), ... 15
                single(NSlices), size_x, zmax, options.verbose, LL, pseudot, det_per_ring, TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), uint32(dec), ... 28
                uint32(options.use_device), filename, uint8(use_raw_data), single(options.cpu_to_gpu_factor), uint32(0), header_directory, options.vaimennus, ... 35
                normalization, pituus(end), uint32(attenuation_correction), uint32(normalization_correction), lor_a, xy_index, z_index, tube_width_xy, crystal_size_z, ... 44
                x_center, y_center, z_center, SinDelayed, randoms, uint32(options.projector_type), options.precompute_lor, n_rays, n_rays3D, dc_z, SinM, ... 55
                logical(options.use_64bit_atomics), varargin{2}, uint8(varargin{3}), options.global_correction_factor, bmin, bmax, Vmax, V, options.use_psf, options);
        end
        if options.verbose
            toc(tStart)
        end
    catch ME
        error(['Error when running MEX-file: ' ME.message])
    end
    
    if type == 0
        for ll = 1 : options.partitions
            apu = tz{1,ll};
            if options.save_iter
                apu(:,1) = options.x0;
                if options.MLEM
                    pz{1,ll} = reshape(apu, Nx, Ny, Nz, Niter + 1);
                elseif options.OSEM
                    pz{2,ll} = reshape(apu, Nx, Ny, Nz, Niter + 1);
                end
            else
                if options.MLEM
                    pz{1,ll} = reshape(apu, Nx, Ny, Nz);
                elseif options.OSEM
                    pz{2,ll} = reshape(apu, Nx, Ny, Nz);
                end
            end
        end
        clear tz
    end
end

end

