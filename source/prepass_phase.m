function [options, D, C_co, C_aco, C_osl, Amin, E] = prepass_phase(options, pituus, index, SinM, pseudot, x, y, xx, yy, z_det, dz, dx, dy, bz, bx, by, NSlices, zmax, size_x, block1, blocks,...
    normalization_correction, randoms_correction, xy_index, z_index, lor_a, lor_orth, summa, LL, is_transposed, x_center, y_center, z_center, ind_size, gaussK, bmin, bmax, Vmax, V)
%PREPASS_PHASE Prepass step for various priors and algorithms
% Computes the necessary variables (e.g. weights) for certain
% algorithms/priors if they have been selected. Also converts various
% values to single precision if implementation 2 has been selected.

if options.sampling > 1 && ~options.use_raw_data && ~options.precompute_lor
    options.Ndist = options.Ndist * options.sampling;
end
Nx = options.Nx;
Ny = options.Ny;
Nz = options.Nz;
epps = options.epps;
attenuation_correction = options.attenuation_correction;
NSinos = options.NSinos;
det_per_ring = options.det_per_ring;
use_raw_data = options.use_raw_data;
verbose = options.verbose;
Ny = uint32(Ny);
Nx = uint32(Nx);
Nz = uint32(Nz);
N = (Nx)*(Ny)*(Nz);
det_per_ring = uint32(det_per_ring);
D = [];
C_co = [];
C_aco = [];
C_osl = [];
Amin = [];
E = [];

if (options.MRP || options.quad || options.Huber || options.TV ||options. FMH || options.L || options.weighted_mean || options.APLS || options.BSREM ...
        || options.ramla || options.MBSREM || options.mramla || options.rosem || options.drama || options.ROSEM_MAP || options.ecosem ...
        || options.cosem || options.acosem || options.AD || any(options.COSEM_OSL) || options.NLM || options.RBI_OSL || options.rbi)
    
    % Compute and/or load necessary variables for the TV regularization
    if options.TV && options.MAP
        options = TVPrepass(options);
    end
    
    if options.TV && options.MAP && options.implementation == 2
        options.alphaTGV = single(options.alphaTGV);
        options.betaTGV = single(options.betaTGV);
        options.NiterTGV = uint32(options.NiterTGV);
    end
    
    % Load necessary variables for the APLS regularization
    if options.APLS && options.MAP
        options = APLSPrepass(options);
    end
    
    % Compute the necessary variables for MRAMLA, RBI-OSL and/or various
    % COSEM algorithms
    % E.g. for COSEM compute the complete data matrix, for RBI-OSL compute
    % the sum of all the rows of the system matrix
    if ((options.mramla || options.MBSREM || options.RBI_OSL || options.rbi) && options.MBSREM_prepass || options.ecosem || options.cosem ...
            || options.acosem || any(options.COSEM_OSL))  && options.implementation == 1
        
        if options.acosem
            C_aco = zeros(double(N), options.subsets);
        end
        if options.cosem || options.ecosem
            C_co = zeros(double(N), options.subsets);
        end
        if any(options.COSEM_OSL)
            C_osl = zeros(double(N), options.subsets);
        end
        if options.acosem || options.cosem || options.ecosem || any(options.COSEM_OSL)
            if options.use_psf
                im_apu = computeConvolution(options.x0(:), options, Nx, Ny, Nz, gaussK);
            else
                im_apu = options.x0(:);
            end
            if options.acosem || options.COSEM_OSL == 1
                im = power(options.x0(:), 1/options.h);
            end
        end
        if options.mramla || options.MBSREM
            if options.precompute_lor == false
                Amin = zeros(options.Nang*options.Ndist*options.NSinos,1);
            else
                Amin = zeros(pituus(end),1);
            end
        end
        
        if ~use_raw_data
            if isempty(pseudot)
                pseudot = uint32(0);
            end
        end
        
        D = zeros(N,1);
        if options.precompute_lor
            if normalization_correction || options.attenuation_correction
                E = zeros(length(lor_a),1);
            else
                E = ones(length(lor_a),1);
            end
        else
            if normalization_correction || options.attenuation_correction
                E = zeros(options.Nang*options.Ndist*options.NSinos,1);
            else
                E = ones(options.Nang*options.Ndist*options.NSinos,1);
            end
        end
        
        if verbose
            disp('Prepass phase for MRAMLA, COSEM, ACOSEM and ECOSEM started')
        end
        if iscell(SinM)
            Sino = SinM{1};
        else
            Sino = SinM;
        end
        
        Sino = Sino(:);
        
        if issparse(Sino)
            Sino = (full(Sino));
        end
        for osa_iter = 1 : options.subsets
            if randoms_correction
                if iscell(options.SinDelayed)
                    SinD = double(options.SinDelayed{1}(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                else
                    SinD = double(options.SinDelayed(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                end
                if issparse(SinD)
                    SinD = (full(SinD));
                end
                SinD = SinD(:);
            else
                SinD = 0;
            end
            if normalization_correction
                norm_input = options.normalization(pituus(osa_iter)+1:pituus(osa_iter + 1));
            else
                norm_input = 0;
            end
            if options.scatter_correction && ~options.subtract_scatter
                scatter_input = double(options.ScatterC(pituus(osa_iter)+1:pituus(osa_iter + 1)));
            else
                scatter_input = 0;
            end
            if options.precompute_lor == false
                if use_raw_data == false
                    if options.projector_type == 1 || options.projector_type == 0
                        if exist('OCTAVE_VERSION','builtin') == 0 && exist('projector_mex','file') == 3
                            [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                zmax, options.vaimennus, options.normalization, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, ...
                                randoms_correction, options.scatter, scatter_input, uint16(0), uint32(0), uint32(0), NSinos, uint16(0), pseudot, det_per_ring, false, int64(0), 0, 0, int64(options.TOF_bins), uint32(0), options.verbose, ...
                                use_raw_data, uint32(2), ind_size, block1, blocks, index(pituus(osa_iter) + 1 : pituus(osa_iter + 1)), uint32(options.projector_type));
                        elseif exist('OCTAVE_VERSION','builtin') == 5 && exist('projector_oct','file') == 3
                            [ lor, indices, alkiot] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                zmax, options.vaimennus, options.normalization, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, ...
                                randoms_correction, options.scatter, scatter_input, uint16(0), uint32(0), uint32(0), NSinos, uint16(0), pseudot, det_per_ring, false, int64(0), 0, 0, int64(options.TOF_bins), uint32(0), options.verbose, ...
                                use_raw_data, uint32(2), ind_size, block1, blocks, index(pituus(osa_iter) + 1 : pituus(osa_iter + 1)), uint32(options.projector_type));
                        else
                            % The below lines allow for pure MATLAB
                            % implemention, i.e. no MEX-files will be
                            % used. Currently the below function
                            % uses parfor-loops (requires parallel
                            % computing toolbox).
                            % NOTE: The below method is not
                            % recommended since it is much slower
                            % method.
                            [ lor, indices, alkiot, discard] = improved_siddon_atten( int32(Ny), int32(Nx), int32(Nz), dx, dz, by, bx, bz, z_det, x, y, yy, xx, ...
                                NSinos, NSlices, options.vaimennus, index(pituus(osa_iter) + 1 : pituus(osa_iter + 1)), pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction);
                            alkiot = cell2mat(alkiot);
                            indices = indices(discard);
                            indices = cell2mat(indices) - 1;
                            lor = lor(:,2);
                        end
                    else
                        error('Unsupported projector type')
                    end
                else
                    L = LL(pituus(osa_iter) * 2 + 1 : pituus(osa_iter + 1) * 2);
                    if options.projector_type == 1
                        if exist('OCTAVE_VERSION','builtin') == 0
                            [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                zmax, options.vaimennus, options.normalization, SinD, uint32(0), attenuation_correction, normalization_correction, ...
                                randoms_correction, options.scatter, scatter_input, uint16(0), uint32(0), uint32(0), NSinos, L, pseudot, det_per_ring, false, int64(0), 0, 0, int64(options.TOF_bins), uint32(0), options.verbose, ...
                                use_raw_data, uint32(2), ind_size, block1, blocks, uint32(0), uint32(options.projector_type));
                        elseif exist('OCTAVE_VERSION','builtin') == 5
                            [ lor, indices, alkiot] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                zmax, options.vaimennus, options.normalization, SinD, uint32(0), attenuation_correction, normalization_correction, ...
                                randoms_correction, options.scatter, scatter_input, uint16(0), uint32(0), uint32(0), NSinos, L, pseudot, det_per_ring, false, int64(0), 0, 0, int64(options.TOF_bins), uint32(0), options.verbose, ...
                                use_raw_data, uint32(2), ind_size, block1, blocks, uint32(0), uint32(options.projector_type));
                        end
                    elseif options.projector_type == 2
                        if exist('OCTAVE_VERSION','builtin') == 0
                            [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                zmax, options.vaimennus, options.normalization, SinD, uint32(0), attenuation_correction, normalization_correction, ...
                                randoms_correction, options.scatter, scatter_input, uint16(0), uint32(0), uint32(0), NSinos, L, pseudot, det_per_ring, false, int64(0), 0, 0, int64(options.TOF_bins), uint32(0), options.verbose, ...
                                use_raw_data, uint32(2), ind_size, block1, blocks, uint32(0), uint32(options.projector_type), options.tube_width_xy, ...
                                x_center, y_center, z_center, options.tube_width_z, int32(0));
                        elseif exist('OCTAVE_VERSION','builtin') == 5
                            [ lor, indices, alkiot] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                zmax, options.vaimennus, options.normalization, SinD, uint32(0), attenuation_correction, normalization_correction, ...
                                randoms_correction, options.scatter, scatter_input, uint16(0), uint32(0), uint32(0), NSinos, L, pseudot, det_per_ring, false, int64(0), 0, 0, int64(options.TOF_bins), uint32(0), options.verbose, ...
                                use_raw_data, uint32(2), ind_size, block1, blocks, uint32(0), uint32(options.projector_type), options.tube_width_xy, ...
                                x_center, y_center, z_center, options.tube_width_z, int32(0));
                        end
                    else
                        error('Unsupported projector type')
                    end
                end
                if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                    lor = repeat_elem(int32(1:length(lor))',int32(lor));
                elseif exist('OCTAVE_VERSION','builtin') == 5
                    lor = repelem(int32(1:length(lor)),int32(lor));
                else
                    lor = repelem(int32(1:length(lor)),int32(lor))';
                end
                uu = double(Sino(index(pituus(osa_iter) + 1 : pituus(osa_iter + 1))));
                
                A_length = length(uu);
                indices=indices + 1;
                if verbose
                    tStart = tic;
                end
                if options.use_fsparse == false
                    A = sparse(double(lor),double(indices),double(alkiot), A_length, double(N));
                else
                    A = fsparse(lor,int32(indices),double(alkiot),[A_length double(N) length(alkiot)]);
                end
                clear indices alkiot lor
                if verbose
                    tElapsed = toc(tStart);
                    disp(['Sparse matrix formation took ' num2str(tElapsed) ' seconds'])
                end
            else
                if use_raw_data
                    L_input = LL(pituus(osa_iter) * 2 + 1 : pituus(osa_iter + 1) * 2);
                    xy_index_input = uint32(0);
                    z_index_input = uint16(0);
                else
                    L_input = uint16(0);
                    xy_index_input = xy_index(pituus(osa_iter)+1:pituus(osa_iter + 1));
                    z_index_input = z_index(pituus(osa_iter)+1:pituus(osa_iter + 1));
                end
                if options.projector_type == 2 || options.projector_type == 3
                    lor2 = [uint64(0); cumsum(uint64(lor_orth(pituus(osa_iter)+1:pituus(osa_iter + 1))))];
                else
                    lor2 = [uint64(0); cumsum(uint64(lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1))))];
                end
                if exist('OCTAVE_VERSION','builtin') == 0
                    [A, ~] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , uint32(NSinos), NSlices, size_x, zmax, options.vaimennus, ...
                        norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction,...
                        randoms_correction, options.scatter, scatter_input, options.global_correction_factor, lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1)), xy_index_input, ...
                        z_index_input, uint32(NSinos), L_input, pseudot, det_per_ring, false, int64(0), 0, 0, int64(options.TOF_bins), uint32(0), options.verbose, use_raw_data, uint32(0), lor2, summa(osa_iter), ...
                        options.attenuation_phase, uint32(options.projector_type), options.tube_width_xy, x_center, y_center, z_center, ...
                        options.tube_width_z, bmin, bmax, Vmax, V);
                elseif exist('OCTAVE_VERSION','builtin') == 5
                    [A, ~] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , uint32(NSinos), NSlices, size_x, zmax, options.vaimennus, ...
                        norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction,...
                        randoms_correction, options.scatter, scatter_input, options.global_correction_factor, lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1)), xy_index_input, ...
                        z_index_input, uint32(NSinos), L_input, pseudot, det_per_ring, false, int64(0), 0, 0, int64(options.TOF_bins), uint32(0), options.verbose, use_raw_data, uint32(0), lor2, summa(osa_iter), ...
                        options.attenuation_phase, uint32(options.projector_type), options.tube_width_xy, x_center, y_center, z_center, ...
                        options.tube_width_z, bmin, bmax, Vmax, V);
                end
%                 uu = double(Sino(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                clear lor2
            end
            if is_transposed
                % Sensitivity image
                D = D + A * ones(size(A,2),1,'double');
                % Required for MRAMLA/MBSREM epsilon value
                if (normalization_correction || options.attenuation_correction) && (options.mramla || options.MBSREM)
                    E(pituus(osa_iter)+1:pituus(osa_iter + 1)) = full(sum(A,1))';
                end
            else
                D = D + full(sum(A,1))';
                if normalization_correction || options.attenuation_correction && (options.mramla || options.MBSREM)
                    E(pituus(osa_iter)+1:pituus(osa_iter + 1)) = full(sum(A,2))';
                end
            end
            if options.ecosem || options.cosem || options.acosem || any(options.COSEM_OSL)
                uu = double(Sino(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                if is_transposed
                    FP = A' * im_apu + epps + SinD;
                    RHS = A * (uu ./ FP);
                else
                    FP = A * im_apu + epps + SinD;
                    RHS = A' * (uu ./ FP);
                end
                if options.use_psf
                    RHS = computeConvolution(RHS, options, Nx, Ny, Nz, gaussK);
                end
            end
            if options.cosem || options.ecosem
                if osa_iter > 1
                    if options.verbose
                        tic
                    end
                    C_co(:,osa_iter) = options.x0(:) .* RHS;
                    if options.verbose
                        disp(['COSEM complete data calculation took ' num2str(toc) ' seconds'])
                    end
                end
            end
            if options.acosem
                if osa_iter > 1
                    if options.verbose
                        tic
                    end
                    C_aco(:,osa_iter) = im .* RHS;
                    if options.verbose
                        disp(['ACOSEM complete data calculation took ' num2str(toc) ' seconds'])
                    end
                end
            end
            if any(options.COSEM_OSL)
                if options.COSEM_OSL == 2
                    if osa_iter > 1
                        C_osl(:,osa_iter) = options.x0(:) .* RHS;
                    end
                else
                    if osa_iter > 1
                        C_osl(:,osa_iter) = im .* RHS;
                    end
                end
            end
            % Required for upper bound of MRAMLA/MBSREM
            if options.MBSREM_prepass && options.U == 0 && (options.MBSREM || options.mramla)
                %%%% This particular piece of code was taken from:
                %%%% https://se.mathworks.com/matlabcentral/answers/35309-max-min-of-sparse-matrices
                if is_transposed
                    [~,m] = size(A);
                    rowMin = nan(m, 1);
                    [~,I,S] = find(A);
                else
                    [m,~] = size(A);
                    rowMin = nan(m, 1);
                    [I,~,S] = find(A);
                end
                I = I(S>1e-10);
                S = S(S>1e-10);
                [I,K] = sort(I);
                S = S(K);
                markers = [find([1; diff(I)]); numel(I)+1];
                iRows = uint32(I(markers(1:end-1)));
                for i = 1:numel(iRows)
                    s = S(markers(i):(markers(i+1)-1));
                    rowMin(iRows(i)) = min(s);
                end
                rowMin(isnan(rowMin)) = epps;
                if options.precompute_lor == false
                    Amin(index{osa_iter}) = rowMin;
                else
                    Amin(pituus(osa_iter)+1:pituus(osa_iter + 1)) = rowMin;
                end
                clear I K S markers rowMin s iRows
            end
            clear A
        end
        %             D = sum(pj,2);
        if verbose
            disp('Prepass phase for COSEM, ACOSEM and ECOSEM completed')
        end
    end
    
    % Lambda values (relaxation parameters)
    if (options.BSREM || options.ramla) && length(options.lambda0) == 1
        lam = zeros(options.Niter,1);
        lam(1) = options.lambda0;
        %             orig_lam = lam;
        %             if lam(1) > 1/max(max(pj))
        %                 lam(1) = min(min(pj));
        %             end
        for i=2:options.Niter
            %                 lam(i) = 0.5*lam(i-1);
            lam(i) = lam(1)/i;
            %                 lam(i) = lam(1)/1.01;
        end
        if options.implementation == 2
            options.lam = single(lam);
        else
            options.lam = lam;
        end
    elseif (options.BSREM || options.ramla) && options.implementation == 2
        options.lam = single(options.lam);
    end
    if (options.MBSREM || options.mramla) && length(options.lambda0_mbsrem) == 1
        lam_mbsrem = zeros(options.Niter,1);
        lam_mbsrem(1) = options.lambda0_mbsrem;
        for i=2:options.Niter
            lam_mbsrem(i) = lam_mbsrem(1)/(i);
        end
        if options.implementation == 2
            options.lam_mbsrem = single(lam_mbsrem);
        else
            options.lam_mbsrem = lam_mbsrem;
        end
    elseif (options.MBSREM || options.mramla) && options.implementation == 2
        options.lam_mbsrem = single(options.lam_mbsrem);
    end
    if (options.ROSEM_MAP || options.rosem) && length(options.lambda0_rosem) == 1
        lam_rosem = zeros(options.Niter,1);
        lam_rosem(1) = options.lambda0_rosem;
        for i=2:options.Niter
            lam_rosem(i) = lam_rosem(1)/i;
        end
        if options.implementation == 2
            options.lam_rosem = single(lam_rosem);
        else
            options.lam_rosem = lam_rosem;
        end
    elseif (options.ROSEM_MAP || options.rosem) && options.implementation == 2
        options.lambda0_rosem = single(options.lambda0_rosem);
    end
    if options.drama
        lam_drama = zeros(options.Niter,options.subsets);
        lam_drama(1,1) = options.beta_drama/(options.alpha_drama*options.beta0_drama);
        r = 1;
        for i=1:options.Niter
            for j = 1 : options.subsets
                lam_drama(i,j) = options.beta_drama/(options.alpha_drama*options.beta0_drama + r);
                r = r + 1;
            end
        end
        if options.implementation == 2
            options.lam_drama = single(lam_drama);
        else
            options.lam_drama = lam_drama;
        end
    end
    % Sensitivity image for MRAMLA/MBSREM
    if (options.MBSREM || options.mramla) && options.implementation == 1
        options.pj3 = D/options.subsets;
    end
    % Compute the weights
    if (options.quad || options.L || options.FMH || options.weighted_mean || options.MRP || (options.TV && options.TVtype == 3) || options.Huber) && options.MAP
        if options.quad || options.L || options.FMH || options.weighted_mean || (options.TV && options.TVtype == 3) || options.Huber
            options = computeWeights(options);
        end
        % These values are needed in order to vectorize the calculation of
        % certain priors
        % Specifies the indices of the center pixel and its neighborhood
%         if (options.MRP && (options.implementation == 2 || ~license('test', 'image_toolbox'))) || options.L || options.FMH || (options.TV && options.TVtype == 3)
            options = computeOffsets(options);
%         else
%             options.tr_offsets = uint32(0);
%         end
        if options.quad || (options.TV && options.TVtype == 3)
            options = quadWeights(options, options.empty_weight);
        end
        if options.Huber
            options = huberWeights(options);
        end
        if options.L
            if isempty(options.a_L)
                options.a_L = lfilter_weights(options.Ndx, options.Ndy, options.Ndz, dx, dy, dz, options.oneD_weights);
            end
            if options.implementation == 2
                options.a_L = single(options.a_L);
            end
        end
        if options.FMH
            options = fmhWeights(options);
        end
        if (options.FMH || options.quad || options.Huber) && options.implementation == 2
            options.weights = single(options.weights);
            options.inffi = uint32(find(isinf(options.weights)) - 1);
        end
        if options.weighted_mean
            options = weightedWeights(options);
        end
        if verbose
            disp('Prepass phase for MRP, quadratic prior, L-filter, FMH and weighted mean completed')
        end
    end
    if options.AD && options.MAP
        if options.implementation == 2
            options.NiterAD = uint32(options.NiterAD);
            options.KAD = single(options.KAD);
            options.TimeStepAD = single(options.TimeStepAD);
            options.FluxType = uint32(options.FluxType);
            options.DiffusionType = uint32(options.DiffusionType);
        else
            if options.FluxType == 1
                options.FluxType = 'exponential';
            elseif options.FluxType == 2
                options.FluxType = 'quadratic';
            end
        end
    end
    if options.NLM && options.MAP
        options = NLMPrepass(options);
    end
end
if options.sampling > 1 && ~options.use_raw_data && ~options.precompute_lor
    options.Ndist = options.Ndist / options.sampling;
end