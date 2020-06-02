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

if t == 1 && osa_iter == 1 && iter == 1
    options.x00 = options.x0;
end

[gaussK, options] = PSFKernel(options);

%%
if options.implementation == 1
    
    Nx = options.Nx;
    Ny = options.Ny;
    Nz = options.Nz;

    if options.precompute_lor == false
        iij = double(0:options.Nx);
        jji = double(0:options.Ny);
        kkj = double(0:options.Nz);
    end
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
    
    if options.precompute_lor == false
        if options.use_raw_data == false
            if options.projector_type == 1 || options.projector_type == 0
                if exist('projector_mex','file') == 3
                [ lor, indices, alkiot] = projector_mex( options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, options.bx, options.bz, options.z_det, ...
                    options.x, options.y, options.dy, options.yy, options.xx , options.NSinos, options.NSlices, options.size_x, options.zmax, options.vaimennus, ...
                    options.normalization, SinD, options.pituus(osa_iter), options.attenuation_correction, options.normalization_correction, options.randoms_correction, ...
                    options.global_correction_factor, uint16(0), uint32(0), uint32(0), options.NSinos, uint16(0), options.pseudot, options.det_per_ring, options.verbose, ...
                    options.use_raw_data, uint32(2), options.ind_size, options.block1, options.blocks, options.index{osa_iter}, uint32(options.projector_type), iij, jji, kkj);
                else
                    % The below lines allow for pure MATLAB
                    % implemention, i.e. no MEX-files will be
                    % used. Currently the below function
                    % uses parfor-loops (requires parallel
                    % computing toolbox).
                    % NOTE: The below method is not
                    % recommended since it is much slower
                    % method.
                    [ lor, indices, alkiot, discard] = improved_siddon_atten( options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, options.bx, options.bz, ...
                        options.z_det, options.x, options.y, options.yy, options.xx, options.NSinos, options.NSlices, options.vaimennus, options.index{osa_iter}, ...
                        options.pituus(osa_iter), options.attenuation_correction);
                    alkiot = cell2mat(alkiot);
                    indices = indices(discard);
                    indices = cell2mat(indices) - 1;
                end
            else
                error('Unsupported projector type')
            end
        else
            L = options.LL(options.pituus(osa_iter) * 2 + 1 : options.pituus(osa_iter + 1) * 2);
            %                         L = L';
            %                         L = L(:);
            if options.projector_type == 1
            [ lor, indices, alkiot] = projector_mex( options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, options.bx, options.bz, options.z_det, options.x, ...
                options.y, options.dy, options.yy, options.xx , options.NSinos, options.NSlices, options.size_x, options.zmax, options.vaimennus, options.normalization, SinD, ...
                uint32(0), options.attenuation_correction, options.normalization_correction, options.randoms_correction, options.global_correction_factor, uint16(0), uint32(0), ...
                uint32(0), options.NSinos, L, options.pseudot, options.det_per_ring, options.verbose, options.use_raw_data, uint32(2), options.ind_size, options.block1, ...
                options.blocks, uint32(0), uint32(options.projector_type), iij, jji, kkj);
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
        uu = double(Sino(pituus(osa_iter) + 1 : pituus(osa_iter + 1)));
        
        A_length = length(uu);
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
        if options.use_raw_data && options.subsets > 1
            L_input = options.LL(options.pituus(osa_iter) * 2 + 1 : options.pituus(osa_iter + 1) * 2);
            xy_index_input = uint32(0);
            z_index_input = uint16(0);
        else
            L_input = uint16(0);
            xy_index_input = options.xy_index(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1));
            z_index_input = options.z_index(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1));
        end
        if (options.projector_type == 2 || options.projector_type == 3) && options.subsets > 1
            lor2 = [0; cumsum(uint64(options.lor_orth(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1))))];
        elseif (options.projector_type == 2 || options.projector_type == 3) && options.subsets == 1
            lor2 = [0; uint64(options.lor_orth)];
        elseif options.projector_type == 1 && options.subsets == 1
            lor2 = [0; uint64(options.lor_a)];
        else
            lor2 = [0; cumsum(uint64(options.lor_a(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1))))];
        end
        [A, ~] = projector_mex( options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, options.bx, options.bz, options.z_det, options.x, options.y, ...
            options.dy, options.yy, options.xx , options.NSinos, options.NSlices, options.size_x, options.zmax, options.vaimennus, options.normalization, SinD, ...
            options.pituus(osa_iter + 1) - options.pituus(osa_iter), options.attenuation_correction, options.normalization_correction, options.randoms_correction, ...
            options.global_correction_factor, options.lor_a(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1)), xy_index_input, z_index_input, options.NSinos, ...
            L_input, options.pseudot, options.det_per_ring, options.verbose, options.use_raw_data, uint32(0), lor2, options.summa(osa_iter), options.attenuation_phase, ...
            uint32(options.projector_type), options.tube_width_xy, options.x_center, options.y_center, options.z_center, options.tube_width_z, int32(0), options.bmin, ...
            options.bmax, options.Vmax, options.V);
        uu = double(Sino(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1)));
        clear lor2
    end
    
    if options.is_transposed
        Summ = full(sum(A,2));
    else
        Summ = full(sum(A,1))';
    end
    Summ(Summ <= 0) = options.epps;
    if options.use_psf
        Summ = computeConvolution(Summ, options, Nx, Ny, Nz, gaussK);
    end
    if options.osem || options.ecosem || options.attenuation_phase
        if options.verbose
            tStart = tic;
        end
        options.im_vectors.OSEM_apu = OSEM_im(options.im_vectors.OSEM_apu, A, options.epps, uu, Summ, SinD, options.is_transposed, options, ...
            Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute MRAMLA
    if options.mramla
        if options.verbose
            tStart = tic;
        end
        options.im_vectors.MRAMLA_apu = MBSREM(options.im_vectors.MRAMLA_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, ...
            options.lam_mbsrem, iter, SinD, options.randoms_correction, options.is_transposed, [], [], options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['MRAMLA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute RAMLA
    if options.ramla
        if options.verbose
            tStart = tic;
        end
        options.im_vectors.RAMLA_apu = BSREM_subiter(options.im_vectors.RAMLA_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed, ...
            options, Nx, Ny, Nz, gaussK);
        if any(options.im_vectors.RAMLA_apu < 0)
            warning('Negative values in RAMLA, lower options.lambda value!')
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RAMLA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute ROSEM
    if options.rosem
        if options.verbose
            tStart = tic;
        end
        options.im_vectors.ROSEM_apu = ROSEM_subiter(options.im_vectors.ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, ...
            options.is_transposed, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute RBI
    if options.rbi
        if options.verbose
            tStart = tic;
        end
        options.im_vectors.RBI_apu = RBI_subiter(options.im_vectors.RBI_apu, A, uu, options.epps, Summ, 0, 0, options.D, SinD, options.is_transposed, ...
            [], [], options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute DRAMA
    if options.drama
        if options.verbose
            tStart = tic;
        end
        options.im_vectors.DRAMA_apu = DRAMA_subiter(options.im_vectors.DRAMA_apu, options.lam_drama, options.epps, iter, Summ, osa_iter, A, uu, SinD, ...
            options.is_transposed, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['DRAMA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute COSEM
    if options.cosem || options.ecosem
        if options.verbose
            tStart = tic;
        end
        [options.im_vectors.COSEM_apu, options.C_co] = COSEM_im(options.im_vectors.COSEM_apu, A, options.epps, uu, options.C_co, options.D, osa_iter, SinD, ...
            options.is_transposed, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute ECOSEM
    if options.ecosem
        if options.verbose
            tStart = tic;
        end
        options.im_vectors.ECOSEM_apu = ECOSEM_im(options.im_vectors.ECOSEM_apu, options.epps, options.D, options.im_vectors.COSEM_apu, options.im_vectors.OSEM_apu);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ECOSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute ACOSEM
    if options.acosem
        if options.verbose
            tStart = tic;
        end
        [options.im_vectors.ACOSEM_apu, options.C_aco] = ACOSEM_im(options.im_vectors.ACOSEM_apu, A, options.epps, uu, options.C_aco, options.D, options.h, osa_iter, ...
            SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ACOSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute OSL with MRP
    if options.MRP && options.OSEM_im
        if options.verbose
            tStart = tic;
        end
        med = MRP(options.im_vectors.MRP_OSL_apu, options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
        options.im_vectors.MRP_OSL_apu = OSEM_im(im_vectors.MRP_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_mrp_osem, med, options.epps), SinD, ...
            options.is_transposed, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute MBSREM with MRP
    if options.MRP && options.MBSREM
        if options.verbose
            tStart = tic;
        end
        med = MRP(options.im_vectors.MRP_MBSREM_apu, options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
        options.im_vectors.MRP_MBSREM_apu = MBSREM(options.im_vectors.MRP_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
            iter, SinD, options.randoms_correction, options.is_transposed, options.beta_mrp_mbsrem, med, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['MBSREM MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute BSREM with MRP
    if options.MRP && options.BSREM
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.ramla
            options.im_vectors.MRP_BSREM_apu = options.im_vectors.RAMLA_apu;
        else
            options.im_vectors.MRP_BSREM_apu = BSREM_subiter(options.im_vectors.MRP_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed, options, ...
                Nx, Ny, Nz, gaussK);
        end
        if any(options.im_vectors.MRP_BSREM_apu < 0)
            error('Negative values in BSREM, lower options.lambda value!')
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['BSREM MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute ROSEM-MAP with MRP
    if options.MRP && options.ROSEM_MAP
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.rosem
            options.im_vectors.MRP_ROSEM_apu = options.im_vectors.ROSEM_apu;
        else
            options.im_vectors.MRP_ROSEM_apu = ROSEM_subiter(options.im_vectors.MRP_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, options.is_transposed, ...
                options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute RBI-OSL with MRP
    if options.MRP && options.RBI_OSL
        if options.verbose
            tStart = tic;
        end
        med = MRP(options.im_vectors.MRP_RBI_apu, options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
        options.im_vectors.MRP_RBI_apu = RBI_subiter(options.im_vectors.MRP_RBI_apu, A, uu, options.epps, Summ, options.beta_mrp_rbi, med, options.D, SinD, ...
            options.is_transposed, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute COSEM-OSL with MRP
    if options.MRP && any(options.COSEM_OSL)
        if options.verbose
            tStart = tic;
        end
        med = MRP(options.im_vectors.MRP_COSEM_apu, options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
        if options.COSEM_OSL == 1
            [options.im_vectors.MRP_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.MRP_COSEM_apu, options.D, options.beta_mrp_cosem, med, options.epps, A, uu, ...
                options.C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        else
            [options.im_vectors.MRP_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.MRP_COSEM_apu, options.D, options.beta_mrp_cosem, med, options.epps, A, uu, ...
                options.C_osl, 0, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute OSL with Quadratic prior
    if options.quad && options.OSEM_im
        if options.verbose
            tStart = tic;
        end
        med = Quadratic_prior(options.im_vectors.Quad_OSL_apu, options.tr_offsets, options.weights, options.weights_quad, options.Nx, options.Ny, options.Nz, ...
            options.Ndx, options.Ndy, options.Ndz);
        options.im_vectors.Quad_OSL_apu = OSEM_im(options.im_vectors.Quad_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_quad_osem, med, options.epps), SinD, ...
            options.is_transposed, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute MBSREM with Quadratic prior
    if options.quad && options.MBSREM
        if options.verbose
            tStart = tic;
        end
        med = Quadratic_prior(options.im_vectors.Quad_MBSREM_apu, options.tr_offsets, options.weights, options.weights_quad, options.Nx, options.Ny, options.Nz, ...
            options.Ndx, options.Ndy, options.Ndz);
        options.im_vectors.Quad_MBSREM_apu = MBSREM(options.im_vectors.Quad_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
            iter, SinD, options.randoms_correction, options.is_transposed, options.beta_quad_mbsrem, med, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['MBSREM Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute BSREM with Quadratic prior
    if options.quad && options.BSREM
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.ramla
            options.im_vectors.Quad_BSREM_apu = options.im_vectors.RAMLA_apu;
        else
            options.im_vectors.Quad_BSREM_apu = BSREM_subiter(options.im_vectors.Quad_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed, ...
                options, Nx, Ny, Nz, gaussK);
        end
        if any(options.im_vectors.Quad_BSREM_apu < 0)
            warning('Negative values in BSREM, lower options.lambda value')
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['BSREM Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute ROSEM-MAP with Quadratic prior
    if options.quad && options.ROSEM_MAP
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.rosem
            options.im_vectors.Quad_ROSEM_apu = options.im_vectors.ROSEM_apu;
        else
            options.im_vectors.Quad_ROSEM_apu = ROSEM_subiter(options.im_vectors.Quad_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, options.is_transposed, ...
                options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute RBI-OSL with Quadratic prior
    if options.quad && options.RBI_OSL
        if options.verbose
            tStart = tic;
        end
        med = Quadratic_prior(options.im_vectors.Quad_RBI_apu, options.tr_offsets, options.weights, options.weights_quad, options.Nx, options.Ny, options.Nz, ...
            options.Ndx, options.Ndy, options.Ndz);
        options.im_vectors.Quad_RBI_apu = RBI_subiter(options.im_vectors.Quad_RBI_apu, A, uu, options.epps, Summ, options.beta_quad_rbi, med, options.D, SinD, ...
            options.is_transposed, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-OSL Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute COSEM-OSL with Quadratic prior
    if options.quad && any(options.COSEM_OSL)
        if options.verbose
            tStart = tic;
        end
        med = Quadratic_prior(options.im_vectors.Quad_COSEM_apu, options.tr_offsets, options.weights, options.weights_quad, options.Nx, options.Ny, options.Nz, ...
            options.Ndx, options.Ndy, options.Ndz);
        if options.COSEM_OSL == 1
            [options.im_vectors.Quad_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.Quad_COSEM_apu, options.D, options.beta_quad_cosem, med, options.epps, A, uu, ...
                options.C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        else
            [options.im_vectors.Quad_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.Quad_COSEM_apu, options.D, options.beta_quad_cosem, med, options.epps, A, uu, ...
                options.C_osl, 0, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-OSL Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute OSL with Huber prior
    if options.Huber && options.OSEM_im
        if verbose
            tStart = tic;
        end
        med = Huber_prior(options.im_vectors.Huber_OSL_apu, options.tr_offsets, options.weights, options.weights_huber, options.Nx, options.Ny, options.Nz, ...
            options.Ndx, options.Ndy, options.Ndz, options.huber_delta);
        options.im_vectors.Huber_OSL_apu = OSEM_im(options.im_vectors.Huber_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_huber_osem, med, options.epps), SinD, ...
            options.is_transposed, options, Nx, Ny, Nz, gaussK);
        if verbose
            tElapsed = toc(tStart);
            disp(['OSL Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute MBSREM with Huber prior
    if options.Huber && options.MBSREM
        if verbose
            tStart = tic;
        end
        med = Huber_prior(options.im_vectors.Huber_MBSREM_apu, options.tr_offsets, options.weights, options.weights_huber, options.Nx, options.Ny, options.Nz, ...
            options.Ndx, options.Ndy, options.Ndz, options.huber_delta);
        options.im_vectors.Huber_MBSREM_apu = MBSREM(options.im_vectors.Huber_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, ...
            options.lam_mbsrem, iter, SinD, options.randoms_correction, options.is_transposed, options.beta_huber_mbsrem, med, options, Nx, Ny, Nz, gaussK);
        if verbose
            tElapsed = toc(tStart);
            disp(['MBSREM Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute BSREM with Huber prior
    if options.Huber && options.BSREM
        if verbose
            tStart = tic;
        end
        if iter == 1 && options.ramla
            options.im_vectors.Huber_BSREM_apu = options.im_vectors.RAMLA_apu;
        else
            options.im_vectors.Huber_BSREM_apu = BSREM_subiter(options.im_vectors.Huber_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed, ...
                options, Nx, Ny, Nz, gaussK);
        end
        if any(options.im_vectors.Huber_BSREM_apu < 0)
            warning('Negative values in BSREM, it is recommended to lower lambda value')
        end
        if verbose
            tElapsed = toc(tStart);
            disp(['BSREM Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute ROSEM-MAP with Huber prior
    if options.Huber && options.ROSEM_MAP
        if verbose
            tStart = tic;
        end
        if iter == 1 && options.rosem
            options.im_vectors.Huber_ROSEM_apu = options.im_vectors.ROSEM_apu;
        else
            options.im_vectors.Huber_ROSEM_apu = ROSEM_subiter(options.im_vectors.Huber_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, options.is_transposed, ...
                options, Nx, Ny, Nz, gaussK);
        end
        if verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute RBI-OSL with Huber prior
    if options.Huber && options.RBI_OSL
        if verbose
            tStart = tic;
        end
        med = Huber_prior(options.im_vectors.Huber_RBI_apu, options.tr_offsets, options.weights, options.weights_huber, options.Nx, options.Ny, options.Nz, ...
            options.Ndx, options.Ndy, options.Ndz, options.huber_delta);
        options.im_vectors.Huber_RBI_apu = RBI_subiter(options.im_vectors.Huber_RBI_apu, A, uu, options.epps, Summ, options.beta_huber_rbi, med, options.D, SinD, ...
            options.is_transposed, options, Nx, Ny, Nz, gaussK);
        if verbose
            tElapsed = toc(tStart);
            disp(['RBI-OSL Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute COSEM-OSL with Huber prior
    if options.Huber && any(options.COSEM_OSL)
        if verbose
            tStart = tic;
        end
        med = Huber_prior(options.im_vectors.Huber_COSEM_apu, options.tr_offsets, options.weights, options.weights_huber, options.Nx, options.Ny, options.Nz, ...
            options.Ndx, options.Ndy, options.Ndz, options.huber_delta);
        if options.COSEM_OSL == 1
            [options.im_vectors.Huber_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.Huber_COSEM_apu, options.D, options.beta_huber_cosem, med, options.epps, A, uu, ...
                options.C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        else
            [options.im_vectors.Huber_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.Huber_COSEM_apu, options.D, options.beta_huber_cosem, med, options.epps, A, uu, ...
                options.C_osl, 0, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        end
        if verbose
            tElapsed = toc(tStart);
            disp(['COSEM-OSL Huber sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute OSL with L-filter prior
    if options.L && options.OSEM_im
        if options.verbose
            tStart = tic;
        end
        med = L_filter(options.im_vectors.L_OSL_apu, options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, ...
            options.epps, options.med_no_norm);
        options.im_vectors.L_OSL_apu = OSEM_im(options.im_vectors.L_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_L_osem, med, options.epps), SinD, ...
            options.is_transposed, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute MBSREM with L-filter prior
    if options.L && options.MBSREM
        if options.verbose
            tStart = tic;
        end
        med = L_filter(options.im_vectors.L_MBSREM_apu, options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, ...
            options.epps, options.med_no_norm);
        options.im_vectors.L_MBSREM_apu = MBSREM(options.im_vectors.L_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
            iter, SinD, options.randoms_correction, options.is_transposed, options.beta_L_mbsrem, med, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['MBSREM L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute BSREM with L-filter prior
    if options.L && options.BSREM
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.ramla
            options.im_vectors.L_BSREM_apu = options.im_vectors.RAMLA_apu;
        else
            options.im_vectors.L_BSREM_apu = BSREM_subiter(options.im_vectors.L_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed, ...
                options, Nx, Ny, Nz, gaussK);
        end
        if any(options.im_vectors.L_BSREM_apu < 0)
            error('Negative values in BSREM, lower options.lambda value')
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['BSREM L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute ROSEM-MAP with L-filter prior
    if options.L && options.ROSEM_MAP
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.rosem
            options.im_vectors.L_ROSEM_apu = options.im_vectors.ROSEM_apu;
        else
            options.im_vectors.L_ROSEM_apu = ROSEM_subiter(options.im_vectors.L_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, options.is_transposed, ...
                options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute RBI-OSL with L-filter prior
    if options.L && options.RBI_OSL
        if options.verbose
            tStart = tic;
        end
        med = L_filter(options.im_vectors.L_RBI_apu, options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, options.epps, ...
            options.med_no_norm);
        options.im_vectors.L_RBI_apu = RBI_subiter(options.im_vectors.L_RBI_apu, A, uu, options.epps, Summ, options.beta_L_rbi, med, options.D, SinD, options.is_transposed, ...
            options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-OSL L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute COSEM-OSL with L-filter prior
    if options.L && any(options.COSEM_OSL)
        if options.verbose
            tStart = tic;
        end
        med = L_filter(options.im_vectors.L_COSEM_apu, options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, options.epps, ...
            options.med_no_norm);
        if options.COSEM_OSL == 1
            [options.im_vectors.L_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.L_COSEM_apu, options.D, options.beta_L_cosem, med, options.epps, A, uu, C_osl, ...
                options.h, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        else
            [options.im_vectors.L_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.L_COSEM_apu, options.D, options.beta_L_cosem, med, options.epps, A, uu, C_osl, 0, ...
                options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-OSL L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute OSL with FMH prior
    if options.FMH && options.OSEM_im
        if options.verbose
            tStart = tic;
        end
        med = FMH(options.im_vectors.FMH_OSL_apu, options.tr_offsets, options.fmh_weights, options.weights, options.Nx, options.Ny, options.Nz, N, ...
            options.Ndx, options.Ndy, options.Ndz, options.epps, options.med_no_norm);
        options.im_vectors.FMH_OSL_apu = OSEM_im(options.im_vectors.FMH_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_fmh_osem, med, options.epps), SinD, ...
            options.is_transposed, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.FMH && options.MBSREM
        if options.verbose
            tStart = tic;
        end
        med = FMH(options.im_vectors.FMH_MBSREM_apu, options.tr_offsets, options.fmh_weights, options.weights, options.Nx, options.Ny, options.Nz, N, ...
            options.Ndx, options.Ndy, options.Ndz, options.epps, options.med_no_norm);
        options.im_vectors.FMH_MBSREM_apu = MBSREM(options.im_vectors.FMH_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
            iter, SinD, options.randoms_correction, options.is_transposed, options.beta_fmh_mbsrem, med, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['MBSREM FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.FMH && options.BSREM
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.ramla
            options.im_vectors.FMH_BSREM_apu = options.im_vectors.RAMLA_apu;
        else
            options.im_vectors.FMH_BSREM_apu = BSREM_subiter(options.im_vectors.FMH_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed, ...
                options, Nx, Ny, Nz, gaussK);
        end
        if any(options.im_vectors.FMH_BSREM_apu < 0)
            error('Negative values in BSREM, lower options.lambda value')
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['BSREM FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.FMH && options.ROSEM_MAP
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.rosem
            options.im_vectors.FMH_ROSEM_apu = options.im_vectors.ROSEM_apu;
        else
            options.im_vectors.FMH_ROSEM_apu = ROSEM_subiter(options.im_vectors.FMH_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, options.is_transposed, ...
                options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.FMH && options.RBI_OSL
        if options.verbose
            tStart = tic;
        end
        med = FMH(options.im_vectors.FMH_RBI_apu, options.tr_offsets, options.fmh_weights, options.weights, options.Nx, options.Ny, options.Nz, N, ...
            options.Ndx, options.Ndy, options.Ndz, options.epps, options.med_no_norm);
        options.im_vectors.FMH_RBI_apu = RBI_subiter(options.im_vectors.FMH_RBI_apu, A, uu, options.epps, Summ, options.beta_fmh_rbi, med, options.D, SinD, options.is_transposed, ...
            options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-OSL FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.FMH && any(options.COSEM_OSL)
        if options.verbose
            tStart = tic;
        end
        med = FMH(options.im_vectors.FMH_COSEM_apu, options.tr_offsets, options.fmh_weights, options.weights, options.Nx, options.Ny, options.Nz, N, ...
            options.Ndx, options.Ndy, options.Ndz, options.epps, options.med_no_norm);
        if options.COSEM_OSL == 1
            [options.im_vectors.FMH_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.FMH_COSEM_apu, options.D, options.beta_fmh_cosem, med, options.epps, A, uu, ...
                options.C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        else
            [options.im_vectors.FMH_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.FMH_COSEM_apu, options.D, options.beta_fmh_cosem, med, options.epps, A, uu, ...
                options.C_osl, 0, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-OSL FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute OSL with weighted mean prior
    if options.weighted_mean && options.OSEM_im
        if options.verbose
            tStart = tic;
        end
        med = Weighted_mean(options.im_vectors.Weighted_OSL_apu, options.tr_offsets, options.weighted_weights, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, ...
            options.mean_type, options.epps, options.w_sum, options.med_no_norm);
        options.im_vectors.Weighted_OSL_apu = OSEM_im(options.im_vectors.Weighted_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_weighted_osem, med, options.epps), SinD, ...
            options.is_transposed, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.weighted_mean && options.MBSREM
        if options.verbose
            tStart = tic;
        end
        med = Weighted_mean(options.im_vectors.Weighted_MBSREM_apu, options.tr_offsets, options.weighted_weights, options.Nx, options.Ny, options.Nz, ...
            options.Ndx, options.Ndy, options.Ndz, options.mean_type, options.epps, options.w_sum, options.med_no_norm);
        options.im_vectors.Weighted_MBSREM_apu = MBSREM(options.im_vectors.Weighted_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, ...
            options.lam_mbsrem, iter, SinD, options.randoms_correction, options.is_transposed, options.beta_weighted_mbsrem, med, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['MBSREM weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.weighted_mean && options.BSREM
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.ramla
            options.im_vectors.Weighted_BSREM_apu = options.im_vectors.RAMLA_apu;
        else
            options.im_vectors.Weighted_BSREM_apu = BSREM_subiter(options.im_vectors.Weighted_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed, ...
                options, Nx, Ny, Nz, gaussK);
        end
        if any(options.im_vectors.Weighted_BSREM_apu < 0)
            error('Negative values in BSREM, lower options.lambda value')
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['BSREM weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.weighted_mean && options.ROSEM_MAP
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.rosem
            options.im_vectors.Weighted_ROSEM_apu = options.im_vectors.ROSEM_apu;
        else
            options.im_vectors.Weighted_ROSEM_apu = ROSEM_subiter(options.im_vectors.Weighted_ROSEM_apu, options.lam_rosem, A, uu, options.epps, iter, Summ, SinD, ...
                options.is_transposed, options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.weighted_mean && options.RBI_OSL
        if options.verbose
            tStart = tic;
        end
        med = Weighted_mean(options.im_vectors.Weighted_RBI_apu, options.tr_offsets, options.weighted_weights, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, ...
            options.mean_type, options.epps, options.w_sum, options.med_no_norm);
        options.im_vectors.Weighted_RBI_apu = RBI_subiter(options.im_vectors.Weighted_RBI_apu, A, uu, options.epps, Summ, options.beta_weighted_rbi, ...
            med, options.D, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-OSL weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.weighted_mean && any(options.COSEM_OSL)
        if options.verbose
            tStart = tic;
        end
        med = Weighted_mean(options.im_vectors.Weighted_COSEM_apu, options.tr_offsets, options.weighted_weights, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, ...
            options.mean_type, options.epps, options.w_sum, options.med_no_norm);
        if options.COSEM_OSL == 1
            [options.im_vectors.Weighted_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.Weighted_COSEM_apu, options.D, options.beta_weighted_cosem, ...
                med, options.epps, A, uu, options.C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        else
            [options.im_vectors.Weighted_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.Weighted_COSEM_apu, options.D, options.beta_weighted_cosem, ...
                med, options.epps, A, uu, options.C_osl, 0, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-OSL weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute OSL with TV prior
    if options.TV && options.OSEM_im
        if options.verbose
            tStart = tic;
        end
        grad = TVpriorFinal(options.im_vectors.TV_OSL_apu, TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, ...
            options.tr_offsets);
        options.im_vectors.TV_OSL_apu = OSEM_im(options.im_vectors.TV_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_TV_osem, grad, options.epps), SinD, ...
            options.is_transposed, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.TV && options.MBSREM
        if options.verbose
            tStart = tic;
        end
        grad = TVpriorFinal(options.im_vectors.TV_MBSREM_apu, TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, ...
            options.tr_offsets);
        options.im_vectors.TV_MBSREM_apu = MBSREM(options.im_vectors.TV_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
            iter, SinD, options.randoms_correction, options.is_transposed, options.beta_TV_mbsrem, grad, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['MBSREM TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.TV && options.BSREM
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.ramla
            options.im_vectors.TV_BSREM_apu = options.im_vectors.RAMLA_apu;
        else
            options.im_vectors.TV_BSREM_apu = BSREM_subiter(options.im_vectors.TV_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed, options, ...
                Nx, Ny, Nz, gaussK);
        end
        if any(options.im_vectors.TV_BSREM_apu < 0)
            error('Negative values in BSREM, lower options.lambda value')
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['BSREM TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.TV && options.ROSEM_MAP
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.rosem
            options.im_vectors.TV_ROSEM_apu = options.im_vectors.ROSEM_apu;
        else
            options.im_vectors.TV_ROSEM_apu = ROSEM_subiter(options.im_vectors.TV_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, options.is_transposed, ...
                options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.TV && options.RBI_OSL
        if options.verbose
            tStart = tic;
        end
        grad = TVpriorFinal(options.im_vectors.TV_RBI_apu, TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, ...
            options.tr_offsets);
        options.im_vectors.TV_RBI_apu = RBI_subiter(options.im_vectors.TV_RBI_apu, A, uu, options.epps, Summ, options.beta_TV_rbi, grad, options.D, SinD, options.is_transposed, ...
            options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-OSL TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.TV && any(options.COSEM_OSL)
        if options.verbose
            tStart = tic;
        end
        grad = TVpriorFinal(options.im_vectors.TV_COSEM_apu, TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, ...
            options.tr_offsets);
        if options.COSEM_OSL == 1
            [options.im_vectors.TV_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.TV_COSEM_apu, options.D, options.beta_TV_cosem, grad, options.epps, A, uu, ...
                options.C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        else
            [options.im_vectors.TV_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.TV_COSEM_apu, options.D, options.beta_TV_cosem, grad, options.epps, A, uu, ...
                options.C_osl, 0, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-OSL TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute OSL with MRP-AD prior
    if options.AD && options.OSEM_im
        if options.verbose
            tStart = tic;
        end
        if osa_iter > 1
            med = AD(options.im_vectors.AD_OSL_apu, options.FluxType, options.Nx, options.Ny, options.Nz, options);
            options.im_vectors.AD_OSL_apu = OSEM_im(options.im_vectors.AD_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_ad_osem, med, options.epps), SinD, ...
            options.is_transposed, options, Nx, Ny, Nz, gaussK);
        else
            options.im_vectors.AD_OSL_apu = OSEM_im(options.im_vectors.AD_OSL_apu, A, options.epps, uu, Summ, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.AD && options.MBSREM
        if options.verbose
            tStart = tic;
        end
        if osa_iter > 1
            med = AD(options.im_vectors.AD_MBSREM_apu, options.FluxType, options.Nx, options.Ny, options.Nz, options);
            options.im_vectors.AD_MBSREM_apu = MBSREM(options.im_vectors.AD_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                iter, SinD, options.randoms_correction, options.is_transposed, options.beta_ad_mbsrem, med, options, Nx, Ny, Nz, gaussK);
        else
            options.im_vectors.AD_MBSREM_apu = MBSREM(options.im_vectors.AD_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
                iter, SinD, options.randoms_correction, options.is_transposed, [], [], options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['MBSREM AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.AD && options.BSREM
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.ramla
            options.im_vectors.AD_BSREM_apu = options.im_vectors.RAMLA_apu;
        else
            options.im_vectors.AD_BSREM_apu = BSREM_subiter(options.im_vectors.AD_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed, ...
                options, Nx, Ny, Nz, gaussK);
        end
        if any(options.im_vectors.AD_BSREM_apu < 0)
            error('Negative values in BSREM, lower options.lambda value')
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['BSREM AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.AD && options.ROSEM_MAP
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.rosem
            options.im_vectors.AD_ROSEM_apu = options.im_vectors.ROSEM_apu;
        else
            options.im_vectors.AD_ROSEM_apu = ROSEM_subiter(options.im_vectors.AD_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, options.is_transposed, ...
                options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.AD && options.RBI_OSL
        if options.verbose
            tStart = tic;
        end
        if osa_iter > 1
            med = AD(options.im_vectors.AD_RBI_apu, options.FluxType, options.Nx, options.Ny, options.Nz, options);
        else
            med = ones(size(options.im_vectors.AD_RBI_apu,1),1);
        end
        options.im_vectors.AD_RBI_apu = RBI_subiter(options.im_vectors.AD_RBI_apu, A, uu, options.epps, Summ, options.beta_ad_rbi, med, options.D, SinD, options.is_transposed, ...
            options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-OSL AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.AD && any(options.COSEM_OSL)
        if options.verbose
            tStart = tic;
        end
        if osa_iter > 1
            med = AD(options.im_vectors.AD_COSEM_apu, options.FluxType, options.Nx, options.Ny, options.Nz, options);
        else
            med = ones(size(options.im_vectors.AD_COSEM_apu,1),1);
        end
        if options.COSEM_OSL == 1
            [options.im_vectors.AD_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.AD_COSEM_apu, options.D, options.beta_ad_cosem, med, options.epps, A, uu, ...
                options.C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        else
            [options.im_vectors.AD_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.AD_COSEM_apu, options.D, options.beta_ad_cosem, med, options.epps, A, uu, ...
                options.C_osl, 0, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-OSL AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute OSL with APLS prior
    if options.APLS && options.OSEM_im
        if options.verbose
            tStart = tic;
        end
        grad = TVpriorFinal(options.im_vectors.APLS_OSL_apu, [], options.Nx, options.Ny, options.Nz, true, options, 4);
        options.im_vectors.APLS_OSL_apu = OSEM_im(options.im_vectors.APLS_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_APLS_osem, grad, options.epps), SinD, ...
            options.is_transposed, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.APLS && options.MBSREM
        if options.verbose
            tStart = tic;
        end
        grad = TVpriorFinal(options.im_vectors.APLS_MBSREM_apu, [], options.Nx, options.Ny, options.Nz, true, options, 4);
        options.im_vectors.APLS_MBSREM_apu = MBSREM(options.im_vectors.APLS_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
            iter, SinD, options.randoms_correction, options.is_transposed, options.beta_APLS_mbsrem, grad, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['MBSREM APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.APLS && options.BSREM
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.ramla
            options.im_vectors.APLS_BSREM_apu = options.im_vectors.RAMLA_apu;
        else
            options.im_vectors.APLS_BSREM_apu = BSREM_subiter(options.im_vectors.APLS_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed, ...
                options, Nx, Ny, Nz, gaussK);
        end
        if any(options.im_vectors.APLS_BSREM_apu < 0)
            error('Negative values in BSREM, lower options.lambda value')
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['BSREM APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.APLS && options.ROSEM_MAP
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.rosem
            options.im_vectors.APLS_ROSEM_apu = options.im_vectors.ROSEM_apu;
        else
            options.im_vectors.APLS_ROSEM_apu = ROSEM_subiter(options.im_vectors.APLS_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, options.is_transposed, ...
                options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.APLS && options.RBI_OSL
        if options.verbose
            tStart = tic;
        end
        grad = TVpriorFinal(options.im_vectors.APLS_RBI_apu, [], options.Nx, options.Ny, options.Nz, true, options, 4);
        options.im_vectors.APLS_RBI_apu = RBI_subiter(options.im_vectors.APLS_RBI_apu, A, uu, options.epps, Summ, SinD, options.beta_APLS_rbi, grad, options.D, ...
            options.is_transposed, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-OSL APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.APLS && any(options.COSEM_OSL)
        if options.verbose
            tStart = tic;
        end
        grad = TVpriorFinal(options.im_vectors.APLS_COSEM_apu, [], options.Nx, options.Ny, options.Nz, true, options, 4);
        if options.COSEM_OSL == 1
            [options.im_vectors.APLS_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.APLS_COSEM_apu, options.D, options.beta_APLS_cosem, grad, A, uu, ...
                options.epps, options.C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        else
            [options.im_vectors.APLS_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.APLS_COSEM_apu, options.D, options.beta_APLS_cosem, grad, A, uu, ...
                options.epps, options.C_osl, 0, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-OSL APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute OSL with TGV prior
    if options.TGV && options.OSEM_im
        if options.verbose
            tStart = tic;
        end
        grad = TGV(options.im_vectors.TGV_OSL_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
        options.im_vectors.TGV_OSL_apu = OSEM_im(options.im_vectors.TGV_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_TGV_osem, grad, options.epps), SinD, ...
            options.is_transposed, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.TGV && options.MBSREM
        if options.verbose
            tStart = tic;
        end
        grad = TGV(options.im_vectors.TGV_MBSREM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
        options.im_vectors.TGV_MBSREM_apu = MBSREM(options.im_vectors.TGV_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
            iter, SinD, options.randoms_correction, options.is_transposed, options.beta_TGV_mbsrem, grad, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['MBSREM TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.TGV && options.BSREM
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.ramla
            options.im_vectors.TGV_BSREM_apu = options.im_vectors.RAMLA_apu;
        else
            options.im_vectors.TGV_BSREM_apu = BSREM_subiter(options.im_vectors.TGV_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed, ...
                options, Nx, Ny, Nz, gaussK);
        end
        if any(options.im_vectors.TGV_BSREM_apu < 0)
            error('Negative values in BSREM, lower options.lambda value')
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['BSREM TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.TGV && options.ROSEM_MAP
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.rosem
            options.im_vectors.TGV_ROSEM_apu = options.im_vectors.ROSEM_apu;
        else
            options.im_vectors.TGV_ROSEM_apu = ROSEM_subiter(options.im_vectors.TGV_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, options.is_transposed, ...
                options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.TGV && options.RBI_OSL
        if options.verbose
            tStart = tic;
        end
        grad = TGV(options.im_vectors.TGV_RBI_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
        options.im_vectors.TGV_RBI_apu = RBI_subiter(options.im_vectors.TGV_RBI_apu, A, uu, options.epps, Summ, options.beta_TGV_rbi, grad, options.D, SinD, options.is_transposed, ...
            options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-OSL TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.TGV && any(options.COSEM_OSL)
        if options.verbose
            tStart = tic;
        end
        grad = TGV(options.im_vectors.TGV_COSEM_apu,options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
        if options.COSEM_OSL == 1
            [options.im_vectors.TGV_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.TGV_COSEM_apu, options.D, options.beta_TGV_cosem, grad, A, uu, ...
                options.epps, options.C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        else
            [options.im_vectors.TGV_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.TGV_COSEM_apu, options.D, options.beta_TGV_cosem, grad, A, uu, ...
                options.epps, options.C_osl, 0, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-OSL TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    % Compute OSL with NLM prior
    if options.NLM && options.OSEM_im
        if options.verbose
            tStart = tic;
        end
        med = NLM(options.im_vectors.NLM_OSL_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
            options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
        options.im_vectors.NLM_OSL_apu = OSEM_im(options.im_vectors.NLM_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_NLM_osem, med, options.epps), SinD, ...
            options.is_transposed, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.NLM && options.MBSREM
        if options.verbose
            tStart = tic;
        end
        med = NLM(options.im_vectors.NLM_MBSREM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
            options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
        options.im_vectors.NLM_MBSREM_apu = MBSREM(options.im_vectors.NLM_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, ...
            iter, SinD, options.randoms_correction, options.is_transposed, options.beta_NLM_mbsrem, med, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['MBSREM NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.NLM && options.BSREM
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.ramla
            options.im_vectors.NLM_BSREM_apu = options.im_vectors.RAMLA_apu;
        else
            options.im_vectors.NLM_BSREM_apu = BSREM_subiter(options.im_vectors.NLM_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed, ...
                options, Nx, Ny, Nz, gaussK);
        end
        if any(options.im_vectors.NLM_BSREM_apu < 0)
            error('Negative values in BSREM, lower options.lambda value')
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['BSREM NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.NLM && options.ROSEM_MAP
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.rosem
            options.im_vectors.NLM_ROSEM_apu = options.im_vectors.ROSEM_apu;
        else
            options.im_vectors.NLM_ROSEM_apu = ROSEM_subiter(options.im_vectors.NLM_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, options.is_transposed, ...
                options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.NLM && options.RBI_OSL
        if options.verbose
            tStart = tic;
        end
        med = NLM(options.im_vectors.NLM_RBI_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
            options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
        options.im_vectors.NLM_RBI_apu = RBI_subiter(options.im_vectors.NLM_RBI_apu, A, uu, options.epps, Summ, options.beta_NLM_rbi, med, options.D, SinD, options.is_transposed, ...
            options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-OSL NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.NLM && any(options.COSEM_OSL)
        if options.verbose
            tStart = tic;
        end
        med = NLM(options.im_vectors.NLM_COSEM_apu, options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, ...
            options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
        if options.COSEM_OSL == 1
            [options.im_vectors.NLM_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.NLM_COSEM_apu, options.D, options.beta_NLM_cosem, med, A, uu, ...
                options.epps, options.C_osl, options.h, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        else
            [options.im_vectors.NLM_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.NLM_COSEM_apu, options.D, options.beta_NLM_cosem, med, A, uu, ...
                options.epps, options.C_osl, 0, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-OSL NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.OSEM_im && options.custom
        if options.verbose
            tStart = tic;
        end
        options.im_vectors.custom_OSL_apu = OSEM_im(options.im_vectors.custom_OSL_apu, A, options.epps, uu, OSL(Summ, options.beta_custom_osem, options.grad_OSEM, options.epps), SinD, ...
            options.is_transposed, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL custom prior sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.MBSREM && options.custom
        if options.verbose
            tStart = tic;
        end
        options.im_vectors.custom_MBSREM_apu = MBSREM(options.im_vectors.custom_MBSREM_apu, options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, ...
            options.lam_mbsrem, iter, SinD, options.randoms_correction, options.is_transposed, ptions.beta_custom_mbsrem, options.grad_MBSREM, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['MBSREM custom prior sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.BSREM && options.custom
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.ramla
            options.im_vectors.custom_BSREM_apu = options.im_vectors.RAMLA_apu;
        else
            options.im_vectors.custom_BSREM_apu = BSREM_subiter(options.im_vectors.custom_BSREM_apu, options.lam, options.epps, iter, A, uu, SinD, options.is_transposed, ...
                options, Nx, Ny, Nz, gaussK);
        end
        if any(options.im_vectors.custom_BSREM(:,iter) < 0)
            error('Negative values in BSREM, lower options.lambda value')
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['BSREM custom prior sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.ROSEM_MAP && options.custom
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.rosem
            options.im_vectors.custom_ROSEM_apu = options.im_vectors.ROSEM_apu;
        else
            options.im_vectors.custom_ROSEM_apu = ROSEM_subiter(options.im_vectors.custom_ROSEM_apu, options.lam_rosem, iter, Summ, options.epps, A, uu, SinD, ...
                options.is_transposed, options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP custom prior sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.RBI_OSL && options.custom
        if options.verbose
            tStart = tic;
        end
        options.im_vectors.custom_RBI_apu = RBI_subiter(options.im_vectors.custom_RBI_apu, A, uu, options.epps, Summ, options.D, SinD, options.is_transposed, ...
            options.beta_custom_rbi, options.grad_RBI, options, Nx, Ny, Nz, gaussK);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-OSL custom prior sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if any(options.COSEM_OSL) && options.custom
        if options.verbose
            tStart = tic;
        end
        if options.COSEM_OSL == 1
            [options.im_vectors.custom_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.custom_COSEM_apu, options.D, options.beta_custom_cosem, options.grad_COSEM, ...
                A, uu, options.epps, options.C_aco, options.h, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        else
            [options.im_vectors.custom_COSEM_apu, options.C_osl] = COSEM_OSL(options.im_vectors.custom_COSEM_apu, options.D, options.beta_custom_cosem, options.grad_COSEM, ...
                A, uu, options.epps, options.C_co, 0, options.COSEM_OSL, osa_iter, SinD, options.is_transposed, options, Nx, Ny, Nz, gaussK);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-OSL custom prior sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
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
    if isfield(options, 'grad_RBI') && options.RBI
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
        if issparse(options.SinDelayed{1})
            for kk = 1 : length(options.SinDelayed)
                options.SinDelayed{kk} = single(full(options.SinDelayed{kk}));
            end
        end
    end
    if options.use_raw_data
        options.xy_index = uint32(0);
        options.z_index = uint16(0);
    else
        if isempty(options.pseudot)
            options.pseudot = uint32(100000);
        end
        options.LL = uint16(0);
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
    header_directory = strcat('-I "', header_directory);
    header_directory = strcat(header_directory,'"');
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
    rekot = [options.rekot; false; false];
    tic
    if ~options.use_CUDA
        [pz] = OpenCL_matrixfree( kernel_path, options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, options.bx, options.bz, options.z_det, options.x, ...
            options.y, options.dy, options.yy(end), options.xx(end) , options.NSinos, single(options.NSlices), options.size_x, options.zmax, options.NSinos, ...
            options.verbose, options.LL, options.pseudot, options.det_per_ring, uint32(options.use_device), uint8(options.use_raw_data), filename, uint32(0), options.use_psf, ...
            header_directory, options.vaimennus, options.normalization, options.pituus, uint32(options.attenuation_correction), uint32(options.normalization_correction), ...
            uint32(options.Niter), uint32(options.subsets), uint8(rekot), single(options.epps), options.lor_a, options.xy_index, options.z_index, any(n_rekos), tube_width_xy, ...
            crystal_size_z, options.x_center, options.y_center, options.z_center, options.SinDelayed, randoms, uint32(options.projector_type), options.precompute_lor, ...
            int32(options.dec), n_rays, n_rays3D, dc_z, options, options.SinM, uint32(options.partitions), logical(options.use_64bit_atomics), n_rekos, n_rekos_mlem, reko_type, ...
            reko_type_mlem, options.global_correction_factor, options.bmin, options.bmax, options.Vmax, options.V, gaussK);
    else
        [pz] = CUDA_matrixfree( kernel_path, options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, options.bx,options. bz, options.z_det, options.x, ...
            options.y, options.dy, options.yy(end), options.xx(end) , options.NSinos, single(options.NSlices), options.size_x, options.zmax, options.NSinos, ...
            options.verbose, options.LL, options.pseudot, options.det_per_ring, uint32(options.use_device), uint8(options.use_raw_data), filename, uint32(0), options.use_psf, ...
            header_directory, options.vaimennus, options.normalization, options.pituus, uint32(options.attenuation_correction), uint32(options.normalization_correction), ...
            uint32(options.Niter), uint32(options.subsets), uint8(rekot), single(options.epps), options.lor_a, options.xy_index, options.z_index, any(n_rekos), tube_width_xy, ...
            crystal_size_z, options.x_center, options.y_center, options.z_center, options.SinDelayed, randoms, uint32(options.projector_type), options.precompute_lor, ...
            int32(options.dec), n_rays, n_rays3D, dc_z, options, options.SinM, uint32(options.partitions), logical(options.use_64bit_atomics), n_rekos, n_rekos_mlem, reko_type, ...
            reko_type_mlem, options.global_correction_factor, options.bmin, options.bmax, options.Vmax, options.V, gaussK);
    end
    toc
    
    options.im_vectors = transfer_im_vectors(options.im_vectors, pz, options, iter);
    pz{end} = 0;
    %     if (options.cosem || options.ecosem)
    %         options.C_co = pz{end-3};
    %     end
    %     if options.acosem
    %         options.C_aco = pz{end-2};
    %     end
    if any(options.COSEM_OSL)
        options.C_osl = pz{end-1};
    end
    if (options.mramla || options.MBSREM || options.RBI_OSL || options.rbi || options.cosem || options.ecosem...
            || options.acosem || any(options.COSEM_OSL)) && options.MBSREM_prepass && osa_iter == 1 && iter == 1 && t == 1
        options.D = pz{end};
    end
    
end
