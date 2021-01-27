function options = custom_prior_reconstruction(options, t, iter, osa_iter)
%% Main reconstruction file for the custom prior file
% This function is used to compute various reconstructions with the
% selected methods and the custom prior gradient.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ville-Veikko Wettenhovi, Samuli Summala
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

if iscell(options.SinM)
    Sino = options.SinM{t};
else
    Sino = options.SinM;
end

Sino = Sino(:);

if issparse(Sino)
    Sino = (full(Sino));
end

if ~isfield(options, 'attenuation_phase')
    options.attenuation_phase = false;
end
if ~isfield(options,'TOF_bins')
    options.TOF_bins = 1;
end
if ~isfield(options,'listmode')
    options.listmode = false;
end
if ~isfield(options,'compute_sensitivity_image')
    options.compute_sensitivity_image = false;
end
TOF = options.TOF_bins > 1 && options.projector_type == 1;

if t == 1 && osa_iter == 1 && iter == 1
    options.x00 = options.x0;
end

[gaussK, options] = PSFKernel(options);

if exist('feature','builtin') == 5
    nCores = uint32(feature('numcores'));
else
    nCores = uint32(1);
end

%%
if options.implementation == 1
    
    Nx = options.Nx;
    Ny = options.Ny;
    Nz = options.Nz;
    
    iij = double(0:options.Nx);
    jji = double(0:options.Ny);
    kkj = double(0:options.Nz);
    if ~options.use_raw_data
        if isempty(options.pseudot)
            options.pseudot = int32(0);
        end
    end
    if options.randoms_correction
        if iscell(options.SinDelayed)
            SinD = options.SinDelayed{llo}(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1));
        else
            SinD = options.SinDelayed(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1));
        end
        if issparse(SinD)
            SinD = (full(SinD));
        end
        SinD = SinD(:);
    else
        SinD = 0;
    end
    
    if options.normalization_correction
        normalization = options.normalization(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1));
    else
        normalization = 0;
    end
    
    if options.scatter_correction && ~options.subtract_scatter
        if options.implementation == 1
            scatter_input = options.ScatterC(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1));
        else
            if iscell(options.SinDelayed)
                options.ScatterFB{1} = {single(options.ScatterC{1}(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1)))};
            else
                options.ScatterFB{1} = {single(options.ScatterC(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1)))};
            end
        end
    else
        scatter_input = 0;
    end
    
    uu = double(Sino(options.pituus(osa_iter) + 1 : options.pituus(osa_iter + 1)));
    koko = length(uu);
    [A,~, Summ] = computeImplementation1(options,options.use_raw_data,options.randoms_correction, options.pituus,osa_iter, options.normalization_correction,...
        options.Nx, options.Ny, options.Nz, options.dx, options.dy, options.dz, options.bx, options.by, options.bz, options.x, options.y, options.z_det, options.xx, options.yy, options.size_x, options.NSinos, options.NSlices, options.zmax, options.attenuation_correction, options.pseudot, options.det_per_ring, ...
        TOF, options.sigma_x, options.TOFCenter, options.dec, nCores, options.ind_size, options.block1, options.blocks, options.index, iij, jji, kkj, options.LL, options.N, options.summa, options.lor_a, options.xy_index, options.z_index, ...
        options.x_center, options.y_center, options.z_center, options.bmin, options.bmax, options.Vmax, options.V, options.lor_orth, gaussK,options.is_transposed, scatter_input, normalization, SinD, koko);
    [options.im_vectors,options.C_co,options.C_aco,options.C_osl] = computeEstimatesImp1(options.im_vectors, options, A, uu, Summ, SinD, options.is_transposed, gaussK, iter, osa_iter, options.C_co, options.C_aco,options.C_osl,...
        options.randoms_correction, options.N, options.Ndx, options.Ndy, options.Ndz, options.D);
    clear A
else
    %%
    options = double_to_single(options);
    rekot = reko_maker(options);
    
    if t == 1 && iter == 1 && osa_iter == 1
        options.im_vectors = initialize_im_vectors(options.im_vectors, iter, options);
    end
    if t > 1 && osa_iter == 1 && iter == 1
        options.x0 = options.x00;
    else
        options.x0 = updateInitialValue(options.im_vectors, options);
    end
    if (options.rbi || options.MBSREM || options.RBI_OSL || options.COSEM_OSL > 0) && t == 1 && iter == 1 && osa_iter == 1
        options.D = zeros(options.N,1,'single');
    end
    
    if isfield(options, 'grad_OSEM') && options.OSL_OSEM
        options.grad_OSEM = single(options.grad_OSEM);
        options.beta_custom_osem = single(options.beta_custom_osem);
        options.custom_osl_apu = options.im_vectors.custom_OSL_apu;
    end
    if isfield(options, 'grad_MLEM') && options.OSL_MLEM
        options.grad_MLEM = single(options.grad_MLEM);
        options.beta_custom_mlem = single(options.beta_custom_mlem);
        options.custom_mlem_apu = options.im_vectors.custom_MLEM_apu;
    end
    if isfield(options, 'grad_BSREM') && options.BSREM
        options.grad_BSREM = single(options.grad_BSREM);
        options.beta_custom_bsrem = single(options.beta_custom_bsrem);
        options.custom_bsrem_apu = options.im_vectors.custom_BSREM_apu;
    end
    if isfield(options, 'grad_MBSREM') && options.MBSREM
        options.grad_MBSREM = single(options.grad_MBSREM);
        options.beta_custom_mbsrem = single(options.beta_custom_mbsrem);
        options.custom_mbsrem_apu = options.im_vectors.custom_MBSREM_apu;
    end
    if isfield(options, 'grad_ROSEM') && options.ROSEM
        options.grad_ROSEM = single(options.grad_ROSEM);
        options.beta_custom_rosem = single(options.beta_custom_rosem);
        options.custom_rosem_apu = options.im_vectors.custom_ROSEM_apu;
    end
    if isfield(options, 'grad_RBI') && options.RBI_OSL
        options.grad_RBI = single(options.grad_RBI);
        options.beta_custom_rbi = single(options.beta_custom_rbi);
        options.custom_rbi_apu = options.im_vectors.custom_RBI_apu;
    end
    if isfield(options, 'grad_COSEM') && any(options.COSEM_OSL)
        options.grad_COSEM = single(options.grad_COSEM);
        options.beta_custom_cosem = single(options.beta_custom_cosem);
        options.custom_cosem_apu = options.im_vectors.custom_COSEM_apu;
    end
    
    %     if (options.cosem || options.ecosem) && t == 1 && osa_iter == 1 && iter == 1
    %         options.C_co = zeros(options.N, options.subsets, 'single');
    %     end
    %     if options.acosem && t == 1 && osa_iter == 1 && iter == 1
    %         options.C_aco = zeros(options.N, options.subsets, 'single');
    %     end
    if any(options.COSEM_OSL) && t == 1 && osa_iter == 1 && iter == 1
        options.C_osl = zeros(options.N, options.subsets, 'single');
    end
    
    if options.partitions == 1
        if ~iscell(options.SinM)
            options.SinM = {options.SinM};
        end
        if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction ...
                && ~options.reconstruct_trues && ~options.reconstruct_scatter
            apu = single(full(options.SinDelayed));
            options.SinDelayed = cell(1,1);
            options.SinDelayed{1} = apu;
        end
        clear apu
    end
    if issparse(options.SinM{1})
        for kk = 1 : length(options.SinM)
            options.SinM{kk} = single(full(options.SinM{kk}));
        end
    elseif isa(options.SinM{1}, 'uint16')
        for kk = 1 : length(options.SinM)
            options.SinM{kk} = single(options.SinM{kk});
        end
    end
    if (options.randoms_correction || options.scatter_correction) && options.corrections_during_reconstruction ...
            && ~options.reconstruct_trues && ~options.reconstruct_scatter
        options.randSize = uint64(diff(pituus));
        if issparse(options.SinDelayed{1})
            for kk = 1 : length(options.SinDelayed)
                options.SinDelayed{kk} = single(full(options.SinDelayed{kk}));
            end
        end
    else
        options.randSize = uint64(1);
    end
    if options.use_raw_data
        options.xy_index = uint32(0);
        options.z_index = uint16(0);
        TOFSize = int64(size(options.LL,1));
    else
        if isempty(options.pseudot)
            options.pseudot = uint32(100000);
        end
        options.LL = uint16(0);
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
    if n_rays * n_rays3D > 1 && ~options.precompute_lor && options.projector_type == 1
        [x,y] = getMultirayCoordinates(options);
        options.x = single(x(:));
        options.y = single(y(:));
    end
    dc_z = single(options.z_det(2,1) - options.z_det(1,1));
    
    
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
        header_directory = strcat('-I "', header_directory);
        header_directory = strcat(header_directory,'"');
    end
    joku = algorithms_char();
    %         n_rekos = uint32(sum(rekot(~contains(joku,'MLEM'))));
    n_rekos = uint32(sum(rekot(cellfun('isempty',strfind(joku,'MLEM')))));
    n_rekos_mlem = uint32(sum(rekot(~cellfun('isempty',strfind(joku,'MLEM')))));
    reko_type = zeros(length(rekot),1,'uint8');
    reko_type(~cellfun('isempty',strfind(joku,'MBSREM'))) = 1;
    reko_type(~cellfun('isempty',strfind(joku,'MRAMLA'))) = 1;
    reko_type(~cellfun('isempty',strfind(joku,'COSEM'))) = 2;
    reko_type(~cellfun('isempty',strfind(joku,'ACOSEM'))) = 3;
    reko_type(~cellfun('isempty',strfind(joku,'OSL-COSEM')) & options.COSEM_OSL == 1) = 2;
    reko_type(~cellfun('isempty',strfind(joku,'OSL-COSEM')) & options.COSEM_OSL == 2) = 3;
    ind = cellfun('isempty',strfind(joku,'MLEM'));
    reko_type = reko_type(rekot & ind);
    joku = joku(rekot & ind);
    if options.ecosem
        if options.cosem && options.osem
            reko_type(~cellfun('isempty',strfind(joku,'ECOSEM'))) = [];
        elseif options.cosem && ~options.osem
            reko_type(~cellfun('isempty',strfind(joku,'ECOSEM'))) = [];
            reko_type = [0;reko_type];
        elseif ~options.cosem && ~options.osem
            reko_type = [0;reko_type];
        end
    end
    reko_type_mlem = zeros(n_rekos_mlem,1,'uint8');
    options.tt = uint32(t - 1);
    options.iter = uint32(iter - 1);
    options.osa_iter = uint32(osa_iter - 1);
    %         filename = [kernel_path(1:end-length(kernel_file) + 3), filename];
    rekot = [options.rekot; false; false; false; false];
    tic
    if ~options.use_CUDA
        [pz] = OpenCL_matrixfree( kernel_path, options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, options.bx, options.bz, options.z_det, options.x, ...
            options.y, options.dy, options.yy(end), options.xx(end) , options.NSinos, single(options.NSlices), options.size_x, options.zmax, options.NSinos, ...
            options.verbose, options.LL, options.pseudot, options.det_per_ring, TOF, TOFSize, options.sigma_x, options.TOFCenter, int64(options.TOF_bins), int32(options.dec), uint32(options.use_device), uint8(options.use_raw_data), filename, uint32(0), options.use_psf, ...
            header_directory, options.vaimennus, options.normalization, options.pituus, uint32(options.attenuation_correction), uint32(options.normalization_correction), ...
            uint32(options.Niter), uint32(options.subsets), uint8(rekot), single(options.epps), options.lor_a, options.xy_index, options.z_index, any(n_rekos), tube_width_xy, ...
            crystal_size_z, options.x_center, options.y_center, options.z_center, options.SinDelayed, randoms, uint32(options.projector_type), options.precompute_lor, ...
            n_rays, n_rays3D, dc_z, options, options.SinM, uint32(options.partitions), logical(options.use_64bit_atomics), n_rekos, n_rekos_mlem, reko_type, ...
            reko_type_mlem, options.global_correction_factor, options.bmin, options.bmax, options.Vmax, options.V, gaussK);
    else
        header_directory = strrep(header_directory,'"','');
        [pz] = CUDA_matrixfree( kernel_path, options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, options.bx,options. bz, options.z_det, options.x, ...
            options.y, options.dy, options.yy(end), options.xx(end) , options.NSinos, single(options.NSlices), options.size_x, options.zmax, options.NSinos, ...
            options.verbose, options.LL, options.pseudot, options.det_per_ring, TOF, TOFSize, options.sigma_x, options.TOFCenter, int64(options.TOF_bins), int32(options.dec), uint32(options.use_device), uint8(options.use_raw_data), filename, uint32(0), options.use_psf, ...
            header_directory, options.vaimennus, options.normalization, options.pituus, uint32(options.attenuation_correction), uint32(options.normalization_correction), ...
            uint32(options.Niter), uint32(options.subsets), uint8(rekot), single(options.epps), options.lor_a, options.xy_index, options.z_index, any(n_rekos), tube_width_xy, ...
            crystal_size_z, options.x_center, options.y_center, options.z_center, options.SinDelayed, randoms, uint32(options.projector_type), options.precompute_lor, ...
            n_rays, n_rays3D, dc_z, options, options.SinM, uint32(options.partitions), logical(options.use_64bit_atomics), n_rekos, n_rekos_mlem, reko_type, ...
            reko_type_mlem, options.global_correction_factor, options.bmin, options.bmax, options.Vmax, options.V, gaussK);
    end
    toc
    
    options.im_vectors = transfer_im_vectors(options.im_vectors, pz, options, iter);
    %     if (options.cosem || options.ecosem)
    %         options.C_co = pz{end-3};
    %     end
    %     if options.acosem
    %         options.C_aco = pz{end-2};
    %     end
    if any(options.COSEM_OSL) && osa_iter == 1 && iter == 1 && t == 1
        options.C_osl = pz{end-5};
    end
    if (options.mramla || options.MBSREM || options.RBI_OSL || options.rbi || options.cosem || options.ecosem...
            || options.acosem || any(options.COSEM_OSL)) && options.MBSREM_prepass && osa_iter == 1 && iter == 1 && t == 1
        options.D = pz{end-4}(:);
    end
    if options.MBSREM && osa_iter == 1 && iter == 1 && t == 1
        options.epsilon_mramla = pz{end-3};
        options.U = pz{end - 2};
    end
    pz{end} = 0;
    
end
