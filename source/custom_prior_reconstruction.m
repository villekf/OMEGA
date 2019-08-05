function options = custom_prior_reconstruction(options, t, iter, osa_iter)
%% Main reconstruction file for the custom prior file
% This function is used to compute various reconstructions with the
% selected methods and the custom prior gradient.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi, Samuli Summala
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

%%
if options.reconstruction_method == 1
    
    
    if ~options.use_raw_data
        if isempty(options.pseudot)
            options.pseudot = int32(0);
        end
    end
    
    if options.precompute_lor == false
        if options.use_raw_data == false
            [ lor, indices, alkiot] = improved_Siddon_algorithm( options.verbose, options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, ...
                options.bx, options.bz, options.z_det, options.x, options.y, options.dy, options.yy, options.xx , options.NSinos, options.NSlices, ...
                options.size_x, options.zmax, options.NSinos, options.ind_size, options.vaimennus, options.index{osa_iter}, options.pituus(osa_iter), ...
                options.attenuation_correction, options.use_raw_data, uint16(0), options.pseudot, options.block1, options.blocks, options.det_per_ring);
        else
            L = options.LL(options.index{osa_iter},:);
            L = L';
            L = L(:);
            [ lor, indices, alkiot] = improved_Siddon_algorithm(options.verbose, options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, ...
                options.bx, options.bz, options.z_det, options.x, options.y, options.dy, options.yy, options.xx , options.NSinos, options.NSlices, ...
                options.size_x, options.zmax, options.NSinos, options.ind_size, options.vaimennus, uint32(0), int32(0), options.attenuation_correction, ...
                options.use_raw_data, L, options.pseudot, options.block1, options.blocks, options.det_per_ring);
        end
        lor = reshape(lor,[],2);
        lor=repelem(int32((lor(:,1))),lor(:,2));
        uu = double(Sino(options.index{osa_iter}));
        
        A_length = length(uu);
        indices=indices + 1;
        if options.verbose
            tStart = tic;
        end
        if options.use_fsparse == false
            A = sparse(double(lor),double(indices),double(alkiot), A_length, double(options.N));
        else
            A = fsparse(lor,indices,double(alkiot),[A_length double(options.N) length(alkiot)]);
        end
        clear indices alkiot lor
        if options.verbose
            tElapsed = toc(tStart);
            disp(['Sparse matrix formation took ' num2str(tElapsed) ' seconds'])
        end
        A = A';
    else
        lor2 = [0; cumsum(uint32(options.lor_a(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1))))];
        if options.use_raw_data == false
            [A,~] = improved_Siddon_algorithm_array(options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, options.bx, options.bz, options.z_det, ...
                options.x, options.y, options.dy, options.yy, options.xx , options.NSinos, options.NSlices, options.size_x, options.zmax, options.vaimennus, ...
                lor2, options.pituus(osa_iter + 1) - options.pituus(osa_iter), options.attenuation_correction, options.lor_a(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1)),...
                options.summa(osa_iter), options.xy_index(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1)), options.z_index(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1)), options.NSinos, ...
                uint16(0), options.pseudot, options.det_per_ring, options.verbose, options.use_raw_data, false);
        else
            [A,~] = improved_Siddon_algorithm_array(options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, options.bx, options.bz, options.z_det, ...
                options.x, options.y, options.dy, options.yy, options.xx , options.NSinos, options.NSlices, options.size_x, options.zmax, options.vaimennus, ...
                lor2, options.pituus(osa_iter + 1) - options.pituus(osa_iter), options.attenuation_correction, options.lor_a(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1)),...
                options.summa(osa_iter), uint32(0), uint32(0), options.NSinos, options.LL(options.pituus(osa_iter) * 2 + 1 : options.pituus(osa_iter + 1) * 2), ...
                options.pseudot, options.det_per_ring, options.verbose, options.use_raw_data, false);
        end
        uu = double(Sino(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1)));
        clear lor2
    end
    
    Summ = full(sum(A,2));
    if options.osem || options.ecosem
        if options.verbose
            tStart = tic;
        end
        options.im_vectors.OSEM(:,iter) = OSEM_im(options.im_vectors.OSEM(:,iter), A, options.epps, uu, Summ);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.mramla
        if options.verbose
            tStart = tic;
        end
        options.im_vectors.MRAMLA(:,iter) = MBSREM(options.im_vectors.MRAMLA(:,iter), options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, 0, 0);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['MRAMLA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.ramla
        if options.verbose
            tStart = tic;
        end
        options.im_vectors.RAMLA(:,iter) = BSREM_subiter(options.im_vectors.RAMLA(:,iter), options.lam, A, uu, options.epps, iter);
        if any(options.im_vectors.RAMLA(:,iter) < 0)
            error('Negative values in RAMLA, lower options.lambda value!')
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RAMLA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.rosem
        if options.verbose
            tStart = tic;
        end
        options.im_vectors.ROSEM(:,iter) = ROSEM_subiter(options.im_vectors.ROSEM(:,iter), options.lam_rosem, A, uu, options.epps, iter, Summ);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.rbi
        if options.verbose
            tStart = tic;
        end
        options.im_vectors.RBI(:,iter) = RBI_subiter(options.im_vectors.RBI(:,iter), A, uu, options.epps, Summ, 0, 0, options.D);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.drama
        if options.verbose
            tStart = tic;
        end
        options.im_vectors.DRAMA(:,iter) = DRAMA_subiter(options.im_vectors.DRAMA(:,iter), options.lam_drama, A, uu, options.epps, iter, Summ, osa_iter);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['DRAMA sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.cosem || options.ecosem
        if options.verbose
            tStart = tic;
        end
        [options.im_vectors.COSEM(:,iter), options.C_co] = COSEM_im(options.im_vectors.COSEM(:,iter), A, options.epps, uu, options.C_co, options.D, osa_iter);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.ecosem
        if options.verbose
            tStart = tic;
        end
        options.im_vectors.ECOSEM(:,iter) = ECOSEM_im(options.im_vectors.ECOSEM(:,iter), options.epps, options.D, options.im_vectors.COSEM(:,iter), options.im_vectors.OSEM(:,iter));
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ECOSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.acosem
        if options.verbose
            tStart = tic;
        end
        [options.im_vectors.ACOSEM(:,iter), options.C_aco] = ACOSEM_im(options.im_vectors.ACOSEM(:,iter), A, options.epps, uu, options.C_aco, options.D, options.h, osa_iter);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ACOSEM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.MRP && options.OSL_OSEM
        if options.verbose
            tStart = tic;
        end
        med = MRP(options.im_vectors.MRP_OSL(:,iter), options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
        options.im_vectors.MRP_OSL(:,iter) = OSL_OSEM(options.im_vectors.MRP_OSL(:,iter), Summ, options.beta_mrp_osem, med, A, uu, options.epps);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.MRP && options.MBSREM
        if options.verbose
            tStart = tic;
        end
        med = MRP(options.im_vectors.MRP_MBSREM(:,iter), options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
        options.im_vectors.MRP_MBSREM(:,iter) = MBSREM(options.im_vectors.MRP_MBSREM(:,iter), options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, options.beta_mrp_mbsrem, med);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['MBSREM MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.MRP && options.BSREM
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.ramla
            options.im_vectors.MRP_BSREM(:,iter) = options.im_vectors.RAMLA(:,iter);
        else
            options.im_vectors.MRP_BSREM(:,iter) = BSREM_subiter(options.im_vectors.MRP_BSREM(:,iter), options.lam, A, uu, options.epps, iter);
        end
        if any(options.im_vectors.MRP_BSREM(:,iter) < 0)
            error('Negative values in BSREM, lower options.lambda value!')
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['BSREM MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.MRP && options.ROSEM_MAP
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.rosem
            options.im_vectors.MRP_ROSEM(:,iter) = options.im_vectors.ROSEM(:,iter);
        else
            options.im_vectors.MRP_ROSEM(:,iter) = ROSEM_subiter(options.im_vectors.MRP_ROSEM(:,iter), options.lam_rosem, A, uu, options.epps, iter, Summ);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.MRP && options.RBI_MAP
        if options.verbose
            tStart = tic;
        end
        med = MRP(options.im_vectors.MRP_RBI(:,iter), options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
        options.im_vectors.MRP_RBI(:,iter) = RBI_subiter(options.im_vectors.MRP_RBI(:,iter), A, uu, options.epps, Summ, options.beta_mrp_rbi, med, options.D);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-MAP MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.MRP && any(options.COSEM_MAP)
        if options.verbose
            tStart = tic;
        end
        med = MRP(options.im_vectors.MRP_COSEM(:,iter), options.medx, options.medy, options.medz, options.Nx, options.Ny, options.Nz, options.epps, options.tr_offsets, options.med_no_norm);
        if options.COSEM_MAP == 1
            [options.im_vectors.MRP_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.MRP_COSEM(:,iter), options.D, options.beta_mrp_cosem, med, A, uu, options.epps, options.C_osl, options.h, options.COSEM_MAP, osa_iter);
        else
            [options.im_vectors.MRP_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.MRP_COSEM(:,iter), options.D, options.beta_mrp_cosem, med, A, uu, options.epps, options.C_osl, 0, options.COSEM_MAP, osa_iter);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-MAP MRP sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.quad && options.OSL_OSEM
        if options.verbose
            tStart = tic;
        end
        med = Quadratic_prior(options.im_vectors.Quad_OSL(:,iter), options.tr_offsets, options.weights, options.weights_quad, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz);
        options.im_vectors.Quad_OSL(:,iter) = OSL_OSEM(options.im_vectors.Quad_OSL(:,iter), Summ, options.beta_quad_osem, med, A, uu, options.epps);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.quad && options.MBSREM
        if options.verbose
            tStart = tic;
        end
        med = Quadratic_prior(options.im_vectors.Quad_MBSREM(:,iter), options.tr_offsets, options.weights, options.weights_quad, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz);
        options.im_vectors.Quad_MBSREM(:,iter) = MBSREM(options.im_vectors.Quad_MBSREM(:,iter), options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, options.beta_quad_mbsrem, med);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['MBSREM Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.quad && options.BSREM
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.ramla
            options.im_vectors.Quad_BSREM(:,iter) = options.im_vectors.RAMLA(:,iter);
        else
            options.im_vectors.Quad_BSREM(:,iter) = BSREM_subiter(options.im_vectors.Quad_BSREM(:,iter), options.lam, A, uu, options.epps, iter);
        end
        if any(options.im_vectors.Quad_BSREM(:,iter) < 0)
            error('Negative values in BSREM, lower options.lambda value')
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['BSREM Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.quad && options.ROSEM_MAP
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.rosem
            options.im_vectors.Quad_ROSEM(:,iter) = options.im_vectors.ROSEM(:,iter);
        else
            options.im_vectors.Quad_ROSEM(:,iter) = ROSEM_subiter(options.im_vectors.Quad_ROSEM(:,iter), options.lam_rosem, A, uu, options.epps, iter, Summ);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.quad && options.RBI_MAP
        if options.verbose
            tStart = tic;
        end
        med = Quadratic_prior(options.im_vectors.Quad_RBI(:,iter), options.tr_offsets, options.weights, options.weights_quad, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz);
        options.im_vectors.Quad_RBI(:,iter) = RBI_subiter(options.im_vectors.Quad_RBI(:,iter), A, uu, options.epps, Summ, options.beta_quad_rbi, med, options.D);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-MAP Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.quad && any(options.COSEM_MAP)
        if options.verbose
            tStart = tic;
        end
        med = Quadratic_prior(options.im_vectors.Quad_COSEM(:,iter), options.tr_offsets, options.weights, options.weights_quad, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz);
        if options.COSEM_MAP == 1
            [options.im_vectors.Quad_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.Quad_COSEM(:,iter), options.D, options.beta_quad_cosem, med, A, uu, options.epps, options.C_osl, options.h, options.COSEM_MAP, osa_iter);
        else
            [options.im_vectors.Quad_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.Quad_COSEM(:,iter), options.D, options.beta_quad_cosem, med, A, uu, options.epps, options.C_osl, 0, options.COSEM_MAP, osa_iter);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-MAP Quadratic sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.L && options.OSL_OSEM
        if options.verbose
            tStart = tic;
        end
        med = L_filter(options.im_vectors.L_OSL(:,iter), options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, options.epps, options.med_no_norm);
        options.im_vectors.L_OSL(:,iter) = OSL_OSEM(options.im_vectors.L_OSL(:,iter), Summ, options.beta_L_osem, med, A, uu, options.epps);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.L && options.MBSREM
        if options.verbose
            tStart = tic;
        end
        med = L_filter(options.im_vectors.L_MBSREM(:,iter), options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, options.epps, options.med_no_norm);
        options.im_vectors.L_MBSREM(:,iter) = MBSREM(options.im_vectors.L_MBSREM(:,iter), options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, options.beta_L_mbsrem, med);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['MBSREM L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.L && options.BSREM
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.ramla
            options.im_vectors.L_BSREM(:,iter) = options.im_vectors.RAMLA(:,iter);
        else
            options.im_vectors.L_BSREM(:,iter) = BSREM_subiter(options.im_vectors.L_BSREM(:,iter), options.lam, A, uu, options.epps, iter);
        end
        if any(options.im_vectors.L_BSREM(:,iter) < 0)
            error('Negative values in BSREM, lower options.lambda value')
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['BSREM L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.L && options.ROSEM_MAP
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.rosem
            options.im_vectors.L_ROSEM(:,iter) = options.im_vectors.ROSEM(:,iter);
        else
            options.im_vectors.L_ROSEM(:,iter) = ROSEM_subiter(options.im_vectors.L_ROSEM(:,iter), options.lam_rosem, A, uu, options.epps, iter, Summ);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.L && options.RBI_MAP
        if options.verbose
            tStart = tic;
        end
        med = L_filter(options.im_vectors.L_RBI(:,iter), options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, options.epps, options.med_no_norm);
        options.im_vectors.L_RBI(:,iter) = RBI_subiter(options.im_vectors.L_RBI(:,iter), A, uu, options.epps, Summ, options.beta_L_rbi, med, options.D);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-MAP L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.L && any(options.COSEM_MAP)
        if options.verbose
            tStart = tic;
        end
        med = L_filter(options.im_vectors.L_COSEM(:,iter), options.tr_offsets, options.a_L, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, options.epps, options.med_no_norm);
        if options.COSEM_MAP == 1
            [options.im_vectors.L_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.L_COSEM(:,iter), options.D, options.beta_L_cosem, med, A, uu, options.epps, options.C_osl, options.h, options.COSEM_MAP, osa_iter);
        else
            [options.im_vectors.L_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.L_COSEM(:,iter), options.D, options.beta_L_cosem, med, A, uu, options.epps, options.C_osl, 0, options.COSEM_MAP, osa_iter);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-MAP L-filter sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.FMH && options.OSL_OSEM
        if options.verbose
            tStart = tic;
        end
        med = FMH(options.im_vectors.FMH_OSL(:,iter), options.tr_offsets, options.fmh_weights, options.weights, options.Nx, options.Ny, options.Nz, options.N, options.Ndx, options.Ndy, options.Ndz, options.epps, options.med_no_norm);
        options.im_vectors.FMH_OSL(:,iter) = OSL_OSEM(options.im_vectors.FMH_OSL(:,iter), Summ, options.beta_fmh_osem, med, A, uu, options.epps);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.FMH && options.MBSREM
        if options.verbose
            tStart = tic;
        end
        med = FMH(options.im_vectors.FMH_MBSREM(:,iter), options.tr_offsets, options.fmh_weights, options.weights, options.Nx, options.Ny, options.Nz, options.N, options.Ndx, options.Ndy, options.Ndz, options.epps, options.med_no_norm);
        options.im_vectors.FMH_MBSREM(:,iter) = MBSREM(options.im_vectors.FMH_MBSREM(:,iter), options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, options.beta_fmh_mbsrem, med);
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
            options.im_vectors.FMH_BSREM(:,iter) = options.im_vectors.RAMLA(:,iter);
        else
            options.im_vectors.FMH_BSREM(:,iter) = BSREM_subiter(options.im_vectors.FMH_BSREM(:,iter), options.lam, A, uu, options.epps, iter);
        end
        if any(options.im_vectors.FMH_BSREM(:,iter) < 0)
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
            options.im_vectors.FMH_ROSEM(:,iter) = options.im_vectors.ROSEM(:,iter);
        else
            options.im_vectors.FMH_ROSEM(:,iter) = ROSEM_subiter(options.im_vectors.FMH_ROSEM(:,iter), options.lam_rosem, A, uu, options.epps, iter, Summ);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.FMH && options.RBI_MAP
        if options.verbose
            tStart = tic;
        end
        med = FMH(options.im_vectors.FMH_RBI(:,iter), options.tr_offsets, options.fmh_weights, options.weights, options.Nx, options.Ny, options.Nz, options.N, options.Ndx, options.Ndy, options.Ndz, options.epps, options.med_no_norm);
        options.im_vectors.FMH_RBI(:,iter) = RBI_subiter(options.im_vectors.FMH_RBI(:,iter), A, uu, options.epps, Summ, options.beta_fmh_rbi, med, options.D);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-MAP FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.FMH && any(options.COSEM_MAP)
        if options.verbose
            tStart = tic;
        end
        med = FMH(options.im_vectors.FMH_COSEM(:,iter), options.tr_offsets, options.fmh_weights, options.weights, options.Nx, options.Ny, options.Nz, options.N, options.Ndx, options.Ndy, options.Ndz, options.epps, options.med_no_norm);
        if options.COSEM_MAP == 1
            [options.im_vectors.FMH_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.FMH_COSEM(:,iter), options.D, options.beta_fmh_cosem, med, A, uu, options.epps, options.C_osl, options.h, options.COSEM_MAP, osa_iter);
        else
            [options.im_vectors.FMH_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.FMH_COSEM(:,iter), options.D, options.beta_fmh_cosem, med, A, uu, options.epps, options.C_osl, 0, options.COSEM_MAP, osa_iter);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-MAP FMH sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.weighted_mean && options.OSL_OSEM
        if options.verbose
            tStart = tic;
        end
        med = Weighted_mean(options.im_vectors.Weighted_OSL(:,iter), options.tr_offsets, options.weighted_weights, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, options.mean_type, options.epps, options.w_sum, options.med_no_norm);
        options.im_vectors.Weighted_OSL(:,iter) = OSL_OSEM(options.im_vectors.Weighted_OSL(:,iter), Summ, options.beta_weighted_osem, med, A, uu, options.epps);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.weighted_mean && options.MBSREM
        if options.verbose
            tStart = tic;
        end
        med = Weighted_mean(options.im_vectors.Weighted_MBSREM(:,iter), options.tr_offsets, options.weighted_weights, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, options.mean_type, options.epps, options.w_sum, options.med_no_norm);
        options.im_vectors.Weighted_MBSREM(:,iter) = MBSREM(options.im_vectors.Weighted_MBSREM(:,iter), options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, options.beta_weighted_mbsrem, med);
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
            options.im_vectors.Weighted_BSREM(:,iter) = options.im_vectors.RAMLA(:,iter);
        else
            options.im_vectors.Weighted_BSREM(:,iter) = BSREM_subiter(options.im_vectors.Weighted_BSREM(:,iter), options.lam, A, uu, options.epps, iter);
        end
        if any(options.im_vectors.Weighted_BSREM(:,iter) < 0)
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
            options.im_vectors.Weighted_ROSEM(:,iter) = options.im_vectors.ROSEM(:,iter);
        else
            options.im_vectors.Weighted_ROSEM(:,iter) = ROSEM_subiter(options.im_vectors.Weighted_ROSEM(:,iter), options.lam_rosem, A, uu, options.epps, iter, Summ);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.weighted_mean && options.RBI_MAP
        if options.verbose
            tStart = tic;
        end
        med = Weighted_mean(options.im_vectors.Weighted_RBI(:,iter), options.tr_offsets, options.weighted_weights, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, options.mean_type, options.epps, options.w_sum, options.med_no_norm);
        options.im_vectors.Weighted_RBI(:,iter) = RBI_subiter(options.im_vectors.Weighted_RBI(:,iter), A, uu, options.epps, Summ, options.beta_weighted_rbi, med, options.D);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-MAP weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.weighted_mean && any(options.COSEM_MAP)
        if options.verbose
            tStart = tic;
        end
        med = Weighted_mean(options.im_vectors.Weighted_COSEM(:,iter), options.tr_offsets, options.weighted_weights, options.Nx, options.Ny, options.Nz, options.Ndx, options.Ndy, options.Ndz, options.mean_type, options.epps, options.w_sum, options.med_no_norm);
        if options.COSEM_MAP == 1
            [options.im_vectors.Weighted_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.Weighted_COSEM(:,iter), options.D, options.beta_weighted_cosem, med, A, uu, options.epps, options.C_osl, options.h, options.COSEM_MAP, osa_iter);
        else
            [options.im_vectors.Weighted_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.Weighted_COSEM(:,iter), options.D, options.beta_weighted_cosem, med, A, uu, options.epps, options.C_osl, 0, options.COSEM_MAP, osa_iter);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-MAP weighted mean sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.TV && options.OSL_OSEM
        if options.verbose
            tStart = tic;
        end
        grad = TVpriorFinal(options.im_vectors.TV_OSL(:,iter), options.TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, options.tr_offsets);
        options.im_vectors.TV_OSL(:,iter) = OSL_OSEM(options.im_vectors.TV_OSL(:,iter), Summ, options.beta_TV_osem, grad, A, uu, options.epps);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.TV && options.MBSREM
        if options.verbose
            tStart = tic;
        end
        grad = TVpriorFinal(options.im_vectors.TV_MBSREM(:,iter), options.TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, options.tr_offsets);
        options.im_vectors.TV_MBSREM(:,iter) = MBSREM(options.im_vectors.TV_MBSREM(:,iter), options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, options.beta_TV_mbsrem, grad);
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
            options.im_vectors.TV_BSREM(:,iter) = options.im_vectors.RAMLA(:,iter);
        else
            options.im_vectors.TV_BSREM(:,iter) = BSREM_subiter(options.im_vectors.TV_BSREM(:,iter), options.lam, A, uu, options.epps, iter);
        end
        if any(options.im_vectors.TV_BSREM(:,iter) < 0)
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
            options.im_vectors.TV_ROSEM(:,iter) = options.im_vectors.ROSEM(:,iter);
        else
            options.im_vectors.TV_ROSEM(:,iter) = ROSEM_subiter(options.im_vectors.TV_ROSEM(:,iter), options.lam_rosem, A, uu, options.epps, iter, Summ);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.TV && options.RBI_MAP
        if options.verbose
            tStart = tic;
        end
        grad = TVpriorFinal(options.im_vectors.TV_RBI(:,iter), options.TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, options.tr_offsets);
        options.im_vectors.TV_RBI(:,iter) = RBI_subiter(options.im_vectors.TV_RBI(:,iter), A, uu, options.epps, Summ, options.beta_TV_rbi, grad, options.D);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-MAP TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.TV && any(options.COSEM_MAP)
        if options.verbose
            tStart = tic;
        end
        grad = TVpriorFinal(options.im_vectors.TV_COSEM(:,iter), options.TVdata, options.Nx, options.Ny, options.Nz, options.TV_use_anatomical, options, options.TVtype, options.tr_offsets);
        if options.COSEM_MAP == 1
            [options.im_vectors.TV_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.TV_COSEM(:,iter), options.D, options.beta_TV_cosem, grad, A, uu, options.epps, options.C_aco, options.h, options.COSEM_MAP, osa_iter);
        else
            [options.im_vectors.TV_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.TV_COSEM(:,iter), options.D, options.beta_TV_cosem, grad, A, uu, options.epps, options.C_co, 0, options.COSEM_MAP, osa_iter);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-MAP TV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.AD && options.OSL_OSEM
        if options.verbose
            tStart = tic;
        end
        if osa_iter > 1
            med = AD(options.im_vectors.AD_OSL(:,iter), options.FluxType, options.Nx, options.Ny, options.Nz, options);
            options.im_vectors.AD_OSL(:,iter) = OSL_OSEM(options.im_vectors.AD_OSL(:,iter), Summ, options.beta_ad_osem, med, A, uu, options.epps);
        else
            options.im_vectors.AD_OSL(:,iter) = OSEM_im(options.im_vectors.AD_OSL(:,iter), A, options.epps, uu, Summ);
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
        med = AD(options.im_vectors.AD_MBSREM(:,iter), options.FluxType, options.Nx, options.Ny, options.Nz, options);
        options.im_vectors.AD_MBSREM(:,iter) = MBSREM(options.im_vectors.AD_MBSREM(:,iter), options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, options.beta_ad_mbsrem, med);
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
            options.im_vectors.AD_BSREM(:,iter) = options.im_vectors.RAMLA(:,iter);
        else
            options.im_vectors.AD_BSREM(:,iter) = BSREM_subiter(options.im_vectors.AD_BSREM(:,iter), options.lam, A, uu, options.epps, iter);
        end
        if any(options.im_vectors.AD_BSREM(:,iter) < 0)
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
            options.im_vectors.AD_ROSEM(:,iter) = options.im_vectors.ROSEM(:,iter);
        else
            options.im_vectors.AD_ROSEM(:,iter) = ROSEM_subiter(options.im_vectors.AD_ROSEM(:,iter), options.lam_rosem, A, uu, options.epps, iter, Summ);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.AD && options.RBI_MAP
        if options.verbose
            tStart = tic;
        end
        med = AD(options.im_vectors.AD_RBI(:,iter), options.FluxType, options.Nx, options.Ny, options.Nz, options);
        options.im_vectors.AD_RBI(:,iter) = RBI_subiter(options.im_vectors.AD_RBI(:,iter), A, uu, options.epps, Summ, options.beta_ad_rbi, med, options.D);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-MAP AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.AD && any(options.COSEM_MAP)
        if options.verbose
            tStart = tic;
        end
        med = AD(options.im_vectors.AD_COSEM(:,iter), options.FluxType, options.Nx, options.Ny, options.Nz, options);
        if options.COSEM_MAP == 1
            [options.im_vectors.AD_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.AD_COSEM(:,iter), options.D, options.beta_ad_cosem, med, A, uu, options.epps, options.C_osl, options.h, options.COSEM_MAP, osa_iter);
        else
            [options.im_vectors.AD_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.AD_COSEM(:,iter), options.D, options.beta_ad_cosem, med, A, uu, options.epps, options.C_osl, 0, options.COSEM_MAP, osa_iter);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-MAP AD sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.APLS && options.OSL_OSEM
        if options.verbose
            tStart = tic;
        end
        grad = TVpriorFinal(options.im_vectors.APLS_OSL(:,iter), [], options.Nx, options.Ny, options.Nz, true, options, 4);
        options.im_vectors.APLS_OSL(:,iter) = OSL_OSEM(options.im_vectors.APLS_OSL(:,iter), Summ, options.beta_APLS_osem, grad, A, uu, options.epps);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.APLS && options.MBSREM
        if options.verbose
            tStart = tic;
        end
        grad = TVpriorFinal(options.im_vectors.APLS_MBSREM(:,iter), [], options.Nx, options.Ny, options.Nz, true, options, 4);
        options.im_vectors.APLS_MBSREM(:,iter) = MBSREM(options.im_vectors.APLS_MBSREM(:,iter), options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, options.beta_APLS_mbsrem, grad);
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
            options.im_vectors.APLS_BSREM(:,iter) = options.im_vectors.RAMLA(:,iter);
        else
            options.im_vectors.APLS_BSREM(:,iter) = BSREM_subiter(options.im_vectors.APLS_BSREM(:,iter), options.lam, A, uu, options.epps, iter);
        end
        if any(options.im_vectors.APLS_BSREM(:,iter) < 0)
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
            options.im_vectors.APLS_ROSEM(:,iter) = options.im_vectors.ROSEM(:,iter);
        else
            options.im_vectors.APLS_ROSEM(:,iter) = ROSEM_subiter(options.im_vectors.APLS_ROSEM(:,iter), options.lam_rosem, A, uu, options.epps, iter, Summ);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.APLS && options.RBI_MAP
        if options.verbose
            tStart = tic;
        end
        grad = TVpriorFinal(options.im_vectors.APLS_RBI(:,iter), [], options.Nx, options.Ny, options.Nz, true, options, 4);
        options.im_vectors.APLS_RBI(:,iter) = RBI_subiter(options.im_vectors.APLS_RBI(:,iter), A, uu, options.epps, Summ, options.beta_APLS_rbi, grad, options.D);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-MAP APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.APLS && any(options.COSEM_MAP)
        if options.verbose
            tStart = tic;
        end
        grad = TVpriorFinal(options.im_vectors.APLS_COSEM(:,iter), [], options.Nx, options.Ny, options.Nz, true, options, 4);
        if options.COSEM_MAP == 1
            [options.im_vectors.APLS_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.APLS_COSEM(:,iter), options.D, options.beta_APLS_cosem, grad, A, uu, options.epps, options.C_aco, options.h, options.COSEM_MAP, osa_iter);
        else
            [options.im_vectors.APLS_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.APLS_COSEM(:,iter), options.D, options.beta_APLS_cosem, grad, A, uu, options.epps, options.C_co, 0, options.COSEM_MAP, osa_iter);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-MAP APLS sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.TGV && options.OSL_OSEM
        if options.verbose
            tStart = tic;
        end
        grad = TGV(options.im_vectors.TGV_OSL(:,iter),options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
        options.im_vectors.TGV_OSL(:,iter) = OSL_OSEM(options.im_vectors.TGV_OSL(:,iter), Summ, options.beta_TGV_osem, grad, A, uu, options.epps);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.TGV && options.MBSREM
        if options.verbose
            tStart = tic;
        end
        grad = TGV(options.im_vectors.TGV_MBSREM(:,iter),options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
        options.im_vectors.TGV_MBSREM(:,iter) = MBSREM(options.im_vectors.TGV_MBSREM(:,iter), options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, options.beta_TGV_mbsrem, grad);
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
            options.im_vectors.TGV_BSREM(:,iter) = options.im_vectors.RAMLA(:,iter);
        else
            options.im_vectors.TGV_BSREM(:,iter) = BSREM_subiter(options.im_vectors.TGV_BSREM(:,iter), options.lam, A, uu, options.epps, iter);
        end
        if any(options.im_vectors.TGV_BSREM(:,iter) < 0)
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
            options.im_vectors.TGV_ROSEM(:,iter) = options.im_vectors.ROSEM(:,iter);
        else
            options.im_vectors.TGV_ROSEM(:,iter) = ROSEM_subiter(options.im_vectors.TGV_ROSEM(:,iter), options.lam_rosem, A, uu, options.epps, iter, Summ);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.TGV && options.RBI_MAP
        if options.verbose
            tStart = tic;
        end
        grad = TGV(options.im_vectors.TGV_RBI(:,iter),options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
        options.im_vectors.TGV_RBI(:,iter) = RBI_subiter(options.im_vectors.TGV_RBI(:,iter), A, uu, options.epps, Summ, options.beta_TGV_rbi, grad, options.D);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-MAP TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.TGV && any(options.COSEM_MAP)
        if options.verbose
            tStart = tic;
        end
        grad = TGV(options.im_vectors.TGV_COSEM(:,iter),options.NiterTGV,options.alphaTGV,options.betaTGV, options.Nx, options.Ny, options.Nz);
        if options.COSEM_MAP == 1
            [options.im_vectors.TGV_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.TGV_COSEM(:,iter), options.D, options.beta_TGV_cosem, grad, A, uu, options.epps, options.C_aco, options.h, options.COSEM_MAP, osa_iter);
        else
            [options.im_vectors.TGV_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.TGV_COSEM(:,iter), options.D, options.beta_TGV_cosem, grad, A, uu, options.epps, options.C_co, 0, options.COSEM_MAP, osa_iter);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-MAP TGV sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.NLM && options.OSL_OSEM
        if options.verbose
            tStart = tic;
        end
        med = NLM(options.im_vectors.NLM_OSL(:,iter), options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
        options.im_vectors.NLM_OSL(:,iter) = OSL_OSEM(options.im_vectors.NLM_OSL(:,iter), Summ, options.beta_NLM_osem, med, A, uu, options.epps);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.NLM && options.MBSREM
        if options.verbose
            tStart = tic;
        end
        med = NLM(options.im_vectors.NLM_MBSREM(:,iter), options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
        options.im_vectors.NLM_MBSREM(:,iter) = MBSREM(options.im_vectors.NLM_MBSREM(:,iter), options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, options.beta_NLM_mbsrem, med);
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
            options.im_vectors.NLM_BSREM(:,iter) = options.im_vectors.RAMLA(:,iter);
        else
            options.im_vectors.NLM_BSREM(:,iter) = BSREM_subiter(options.im_vectors.NLM_BSREM(:,iter), options.lam, A, uu, options.epps, iter);
        end
        if any(options.im_vectors.NLM_BSREM(:,iter) < 0)
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
            options.im_vectors.NLM_ROSEM(:,iter) = options.im_vectors.ROSEM(:,iter);
        else
            options.im_vectors.NLM_ROSEM(:,iter) = ROSEM_subiter(options.im_vectors.NLM_ROSEM(:,iter), options.lam_rosem, A, uu, options.epps, iter, Summ);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.NLM && options.RBI_MAP
        if options.verbose
            tStart = tic;
        end
        med = NLM(options.im_vectors.NLM_RBI(:,iter), options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
        options.im_vectors.NLM_RBI(:,iter) = RBI_subiter(options.im_vectors.NLM_RBI(:,iter), A, uu, options.epps, Summ, options.beta_NLM_rbi, med, options.D);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-MAP NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.NLM && any(options.COSEM_MAP)
        if options.verbose
            tStart = tic;
        end
        med = NLM(options.im_vectors.NLM_COSEM(:,iter), options.Ndx, options.Ndy, options.Ndz, options.Nlx, options.Nly, options.Nlz, options.sigma, options.epps, options.Nx, options.Ny, options.Nz, options);
        if options.COSEM_MAP == 1
            [options.im_vectors.NLM_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.NLM_COSEM(:,iter), options.D, options.beta_NLM_cosem, med, A, uu, options.epps, options.C_aco, options.h, options.COSEM_MAP, osa_iter);
        else
            [options.im_vectors.NLM_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.NLM_COSEM(:,iter), options.D, options.beta_NLM_cosem, med, A, uu, options.epps, options.C_co, 0, options.COSEM_MAP, osa_iter);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-MAP NLM sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.OSL_OSEM
        if options.verbose
            tStart = tic;
        end
        options.im_vectors.custom_OSEM(:,iter) = OSL_OSEM(options.im_vectors.custom_OSEM(:,iter), Summ, options.beta_custom_osem, options.grad_OSEM, A, uu, options.epps);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['OSL custom prior sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.MBSREM
        if options.verbose
            tStart = tic;
        end
        options.im_vectors.custom_MBSREM(:,iter) = MBSREM(options.im_vectors.custom_MBSREM(:,iter), options.U, options.pj3, A, options.epps, uu, options.epsilon_mramla, options.lam_mbsrem, iter, options.beta_custom_mbsrem, options.grad_MBSREM);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['MBSREM custom prior sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.BSREM
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.ramla
            options.im_vectors.custom_BSREM(:,iter) = options.im_vectors.RAMLA(:,iter);
        else
            options.im_vectors.custom_BSREM(:,iter) = BSREM_subiter(options.im_vectors.custom_BSREM(:,iter), options.lam, A, uu, options.epps, iter);
        end
        if any(options.im_vectors.custom_BSREM(:,iter) < 0)
            error('Negative values in BSREM, lower options.lambda value')
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['BSREM custom prior sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.ROSEM_MAP
        if options.verbose
            tStart = tic;
        end
        if iter == 1 && options.rosem
            options.im_vectors.custom_ROSEM(:,iter) = options.im_vectors.ROSEM(:,iter);
        else
            options.im_vectors.custom_ROSEM(:,iter) = ROSEM_subiter(options.im_vectors.custom_ROSEM(:,iter), options.lam_rosem, A, uu, options.epps, iter, Summ);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['ROSEM-MAP custom prior sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if options.RBI_MAP
        if options.verbose
            tStart = tic;
        end
        options.im_vectors.custom_RBI(:,iter) = RBI_subiter(options.im_vectors.custom_RBI(:,iter), A, uu, options.epps, Summ, options.beta_custom_rbi, options.grad_RBI, options.D);
        if options.verbose
            tElapsed = toc(tStart);
            disp(['RBI-MAP custom prior sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    if any(options.COSEM_MAP)
        if options.verbose
            tStart = tic;
        end
        if options.COSEM_MAP == 1
            [options.im_vectors.custom_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.custom_COSEM(:,iter), options.D, options.beta_custom_cosem, options.grad_COSEM, A, uu, options.epps, options.C_aco, options.h, options.COSEM_MAP, osa_iter);
        else
            [options.im_vectors.custom_COSEM(:,iter), options.C_osl] = COSEM_OSL(options.im_vectors.custom_COSEM(:,iter), options.D, options.beta_custom_cosem, options.grad_COSEM, A, uu, options.epps, options.C_co, 0, options.COSEM_MAP, osa_iter);
        end
        if options.verbose
            tElapsed = toc(tStart);
            disp(['COSEM-MAP custom prior sub-iteration ' num2str(osa_iter) ' took ' num2str(tElapsed) ' seconds'])
        end
    end
    clear A
else
    %%
    options = double_to_single(options);
        
        Sino = single(full(Sino));
        if options.use_raw_data
            options.xy_index = uint32(0);
            options.z_index = uint32(0);
        else
            if isempty(options.pseudot)
                options.pseudot = int32(100000);
            end
            options.LL = uint16(0);
        end
        if isfield(options, 'grad_OSEM')
            options.grad_OSEM = single(options.grad_OSEM);
            options.beta_custom_osem = single(options.beta_custom_osem);
        end
        if isfield(options, 'grad_MLEM')
            options.grad_MLEM = single(options.grad_MLEM);
            options.beta_custom_mlem = single(options.beta_custom_mlem);
        end
        if isfield(options, 'grad_BSREM')
            options.grad_BSREM = single(options.grad_BSREM);
            options.beta_custom_bsrem = single(options.beta_custom_bsrem);
        end
        if isfield(options, 'grad_MBSREM')
            options.grad_MBSREM = single(options.grad_MBSREM);
            options.beta_custom_mbsrem = single(options.beta_custom_mbsrem);
        end
        if isfield(options, 'grad_ROSEM')
            options.grad_ROSEM = single(options.grad_ROSEM);
            options.beta_custom_rosem = single(options.beta_custom_rosem);
        end
        if isfield(options, 'grad_RBI')
            options.grad_RBI = single(options.grad_RBI);
            options.beta_custom_rbi = single(options.beta_custom_rbi);
        end
        if isfield(options, 'grad_COSEM')
            options.grad_COSEM = single(options.grad_COSEM);
            options.beta_custom_cosem = single(options.beta_custom_cosem);
        end
        %         size_x = int32(numel(x));
        
        
        options.im_vectors = initialize_im_vectors(options.im_vectors, iter, options);
        
        if (options.cosem || options.ecosem) && t == 1 && osa_iter == 1 && iter == 1
            options.C_co = zeros(options.N, options.subsets, 'single');
        end
        if options.acosem && t == 1 && osa_iter == 1 && iter == 1
            options.C_aco = zeros(options.N, options.subsets, 'single');
        end
        if any(options.COSEM_MAP) && t == 1 && osa_iter == 1 && iter == 1
            options.C_osl = zeros(options.N, options.subsets, 'single');
        end
        
        
        %         if use_raw_data == false
        kernel_file = 'siddon_kernel_matrixfree_custom.cl';
%         kernel_file = 'siddon_kernel_matrixfree_GPU - Copy (2).cl';
        kernel_path = which(kernel_file);
        kernel_path = strrep(kernel_path, '\', '/');
        filename = 'OMEGA_matrix_free_OpenCL_custom_binary_device';
        filename = [kernel_path(1:end-length(kernel_file)), filename];
        mlem_bool = options.mlem || options.OSL_MLEM;
        if osa_iter == 1 && mlem_bool
            [zz] = openCL_matrixfree_custom( kernel_path, options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, options.bx, options.bz, ...
                options.z_det, options.x, options.y, options.dy, options.yy(end), options.xx(end), options.NSinos, ...
                single(options.NSlices), options.size_x, options.zmax, options.NSinos, options.vaimennus, options.pituus, int32(options.attenuation_correction), int32(iter - 1), int32(options.subsets), ...
                uint8(options.rekot), single(options.epps), Sino, options, options.lor_a, options.xy_index, options.z_index, options.verbose, options.LL, options.pseudot, options.det_per_ring, ...
                int32(options.use_device), int32(osa_iter - 1), uint8(options.use_raw_data), false, filename, options.force_build, int32(t - 1), ...
                int32(subsets), mlem_bool);
        end
        if any(sum(options.rekot(11:end)))
            [zz] = openCL_matrixfree_custom( kernel_path, options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, options.bx, options.bz, ...
                options.z_det, options.x, options.y, options.dy, options.yy(end), options.xx(end), options.NSinos, ...
                single(options.NSlices), options.size_x, options.zmax, options.NSinos, options.vaimennus, options.pituus, int32(options.attenuation_correction), int32(iter - 1), int32(options.subsets), ...
                uint8(options.rekot), single(options.epps), Sino, options, options.lor_a, options.xy_index, options.z_index, options.verbose, options.LL, options.pseudot, options.det_per_ring, ...
                int32(options.use_device), int32(osa_iter - 1), uint8(options.use_raw_data), any(sum(options.rekot(11:end))), filename, options.force_build, int32(t - 1), ...
                int32(1), false);
        end
        
        options.im_vectors = transfer_im_vectors(options.im_vectors, zz, iter, options);
        if (options.cosem || options.ecosem)
            options.C_co = zz{end-3};
        end
        if options.acosem
            options.C_aco = zz{end-2};
        end
        if any(options.COSEM_MAP)
            options.C_osl = zz{end-1};
        end
        if (options.mramla || options.MBSREM || options.rbi || options.RBI_MAP || options.cosem || options.ecosem...
                || options.acosem || any(options.COSEM_MAP)) && options.MBSREM_prepass && osa_iter == 1 && iter == 1 && t == 1
            options.D = zz{end};
        end
        
end
