function [options, D, C_co, C_aco, C_osl, Amin, E] = prepass_phase(options, pituus, index, SinM, pseudot, x, y, xx, yy, z_det, dz, dx, dy, bz, bx, by, NSlices, zmax, size_x, block1, blocks,...
    normalization_correction, randoms_correction, xy_index, z_index, lor_a, lor_orth, summa, LL, is_transposed, x_center, y_center, z_center)
%PREPASS_PHASE Prepass step for various priors and algorithms
% Computes the necessary variables (e.g. weights) for certain
% algorithms/priors if they have been selected

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

if (options.MRP || options.quad || options.TV ||options. FMH || options.L || options.weighted_mean || options.APLS || options.BSREM ...
        || options.ramla || options.MBSREM || options.mramla || options.rosem || options.drama || options.ROSEM_MAP || options.ecosem ...
        || options.cosem || options.acosem || options.AD || any(options.COSEM_MAP) || (options.NLM && options.NLM_use_anatomical))
    
    % Compute and/or load necessary variables for the TV regularization
    if options.TV && options.MAP
        if options.TV_use_anatomical
            apu = load(options.TV_reference_image);
            variables = fields(apu);
            alkuarvo = double(apu.(variables{1}));
            if size(alkuarvo,2) == 1
                koko_apu = sqrt(length(alkuarvo)/double(Nz));
                if floor(koko_apu) ~= koko_apu
                    error('Reference image has to be square')
                else
                    alkuarvo = reshape(alkuarvo, koko_apu,koko_apu,Nz);
                    if koko_apu ~= Nx
                        alkuarvo = imresize(alkuarvo, Nx, Ny, Nz);
                    end
                end
            else
                if size(alkuarvo,2) ~= Nx
                    alkuarvo = imresize(alkuarvo, Nx, Ny, Nz);
                end
            end
            alkuarvo = alkuarvo - min(min(min(alkuarvo)));
            alkuarvo = alkuarvo/max(max(max(alkuarvo)));
            if options.TVtype == 1
                S = assembleS(alkuarvo,options.T,Ny,Nx,Nz);
                if options.implementation == 2
                    TVdata.s1 = single(S(1:3:end,1));
                    TVdata.s2 = single(S(1:3:end,2));
                    TVdata.s3 = single(S(1:3:end,3));
                    TVdata.s4 = single(S(2:3:end,1));
                    TVdata.s5 = single(S(2:3:end,2));
                    TVdata.s6 = single(S(2:3:end,3));
                    TVdata.s7 = single(S(3:3:end,1));
                    TVdata.s8 = single(S(3:3:end,2));
                    TVdata.s9 = single(S(3:3:end,3));
                else
                    TVdata.s1 = S(1:3:end,1);
                    TVdata.s2 = S(1:3:end,2);
                    TVdata.s3 = S(1:3:end,3);
                    TVdata.s4 = S(2:3:end,1);
                    TVdata.s5 = S(2:3:end,2);
                    TVdata.s6 = S(2:3:end,3);
                    TVdata.s7 = S(3:3:end,1);
                    TVdata.s8 = S(3:3:end,2);
                    TVdata.s9 = S(3:3:end,3);
                end
            end
            if options.implementation == 2
                TVdata.reference_image = single(alkuarvo);
                TVdata.T = single(options.T);
                TVdata.C = single(options.C);
            else
                TVdata.reference_image = alkuarvo;
                TVdata.T = options.T;
                TVdata.C = options.C;
            end
            clear apu variables alkuarvo S
        end
        if options.implementation == 2
            options.tau = single(options.tau);
            TVdata.beta = single(options.TVsmoothing);
            TVdata.C = single(options.C);
            TVdata.T = single(options.T);
            options.TVdata = TVdata;
        else
            TVdata.beta = options.TVsmoothing;
            options.TVdata = TVdata;
        end
        clear TVdata;
    end
    
    if options.TV && options.MAP && options.implementation == 2
        options.alphaTGV = single(options.alphaTGV);
        options.betaTGV = single(options.betaTGV);
        options.NiterTGV = uint32(options.NiterTGV);
    end
    
    % Load necessary variables for the APLS regularization
    if options.APLS && options.MAP
        apu = load(options.APLS_reference_image);
        variables = fields(apu);
        alkuarvo = double(apu.(variables{1}));
        if size(alkuarvo,2) == 1
            koko_apu = sqrt(length(alkuarvo)/double(Nz));
            if floor(koko_apu) ~= koko_apu
                error('Reference image has to be square')
            else
                alkuarvo = reshape(alkuarvo, koko_apu,koko_apu,Nz);
                if koko_apu ~= Nx
                    alkuarvo = imresize(alkuarvo, Nx, Ny, Nz);
                end
            end
        else
            if size(alkuarvo,2) ~= Nx
                alkuarvo = imresize(alkuarvo, Nx, Ny, Nz);
            end
        end
        alkuarvo = alkuarvo - min(min(min(alkuarvo)));
        alkuarvo = alkuarvo/max(max(max(alkuarvo)));
        if options.implementation == 2
            options.APLS_ref_image = single(alkuarvo);
            options.eta = single(options.eta);
            options.APLSsmoothing = single(options.APLSsmoothing);
        else
            options.APLS_ref_image = double(alkuarvo);
        end
        clear alkuarvo apu variables
    end
    
    % Compute the necessary variables for MRAMLA, RBI and/or various
    % COSEM algorithms
    % E.g. for COSEM compute the complete data matrix, for RBI compute
    % the sum of all the rows of the system matrix
    if ((options.mramla || options.MBSREM || options.rbi || options.RBI_MAP) && options.MBSREM_prepass || options.ecosem || options.cosem ...
            || options.acosem || any(options.COSEM_MAP))  && options.implementation == 1
        
        if options.acosem
            C_aco = zeros(double(N), options.subsets);
        end
        if options.cosem || options.ecosem
            C_co = zeros(double(N), options.subsets);
        end
        if any(options.COSEM_MAP)
            C_osl = zeros(double(N), options.subsets);
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
        if normalization_correction || options.attenuation_correction
            E = zeros(options.Nang*options.Ndist*options.NSinos,1);
        else
            E = ones(options.Nang*options.Ndist*options.NSinos,1);
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
            if options.precompute_lor == false
                if use_raw_data == false
                    if options.projector_type == 1 || options.projector_type == 0
                        if exist('projector_mex','file') == 3
                            [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                                zmax, options.vaimennus, options.normalization, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, ...
                                randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, uint16(0), pseudot, det_per_ring, options.verbose, ...
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
                            [ lor, indices, alkiot, discard] = improved_siddon_atten( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, yy, xx, ...
                                NSinos, NSlices, options.vaimennus, index(pituus(osa_iter) + 1 : pituus(osa_iter + 1)), pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction);
                            alkiot = cell2mat(alkiot);
                            indices = indices(discard);
                            indices = cell2mat(indices) - 1;
                        end
                    elseif options.projector_type == 2
                        [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                            zmax, options.vaimennus, options.normalization, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, ...
                            randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, uint16(0), pseudot, det_per_ring, options.verbose, ...
                            use_raw_data, uint32(2), ind_size, block1, blocks, index(pituus(osa_iter) + 1 : pituus(osa_iter + 1)), uint32(options.projector_type), ...
                            options.tube_width_xy, x_center, y_center, z_center, options.tube_width_z, int32(options.accuracy_factor));
                    else
                        error('Unsupported projector type')
                    end
                else
                    L = LL(pituus(osa_iter) * 2 + 1 : pituus(osa_iter + 1) * 2);
                    if options.projector_type == 1
                        [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                            zmax, options.vaimennus, options.normalization, SinD, uint32(0), attenuation_correction, normalization_correction, ...
                            randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, L, pseudot, det_per_ring, options.verbose, ...
                            use_raw_data, uint32(2), ind_size, block1, blocks, uint32(0), uint32(options.projector_type));
                    elseif options.projector_type == 2
                        [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                            zmax, options.vaimennus, options.normalization, SinD, uint32(0), attenuation_correction, normalization_correction, ...
                            randoms_correction, uint16(0), uint32(0), uint32(0), NSinos, L, pseudot, det_per_ring, options.verbose, ...
                            use_raw_data, uint32(2), ind_size, block1, blocks, uint32(0), uint32(options.projector_type), options.tube_width_xy, ...
                            x_center, y_center, z_center, options.tube_width_z, int32(options.accuracy_factor));
                    else
                        error('Unsupported projector type')
                    end
                end
                %                     lor = reshape(lor,[],2);
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
                if options.projector_type == 2
                    lor2 = [uint64(0); cumsum(uint64(lor_orth(pituus(osa_iter)+1:pituus(osa_iter + 1))))];
                else
                    lor2 = [uint64(0); cumsum(uint64(lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1))))];
                end
                [A, ~] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
                    options.normalization, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, randoms_correction, ...
                    lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1)), xy_index_input, z_index_input, NSinos, L_input, pseudot, det_per_ring, options.verbose, ...
                    use_raw_data, uint32(0), lor2, summa(osa_iter), options.attenuation_phase, uint32(options.projector_type), options.tube_width_xy, x_center, ...
                    y_center, z_center, options.tube_width_z, int32(options.accuracy_factor));
                uu = double(Sino(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                clear lor2
            end
            if is_transposed
                D = D + A * ones(size(A,2),1,'double');
                if normalization_correction || options.attenuation_correction
                    if options.precompute_lor
                        E(index{osa_iter}) = full(sum(A,1))';
                    else
                        E(pituus(osa_iter)+1:pituus(osa_iter + 1)) = full(sum(A,1))';
                    end
                end
            else
                D = D + full(sum(A,1))';
                if normalization_correction || options.attenuation_correction
                    if options.precompute_lor
                        E(index{osa_iter}) = full(sum(A,2))';
                    else
                        E(pituus(osa_iter)+1:pituus(osa_iter + 1)) = full(sum(A,2))';
                    end
                end
            end
            if options.ecosem || options.cosem || options.acosem || any(options.COSEM_MAP)
                %                     if options.precompute_lor == false
                %                         uu = double(Sino(index{osa_iter}));
                %                     else
                uu = double(Sino(pituus(osa_iter)+1:pituus(osa_iter + 1)));
                %                     end
            end
            if options.cosem || options.ecosem
                if osa_iter > 1
                    if is_transposed
                        C_co(:,osa_iter) = full(sum(spdiags(options.x0(:),0,size(A,1),size(A,1)) * A * spdiags(uu ./ (A' * options.x0(:) + epps + SinD),...
                            0,size(A,2),size(A,2)),2));
                    else
                        C_co(:,osa_iter) = full(sum(spdiags(options.x0(:),0,size(A,2),size(A,2)) * A' * spdiags(uu ./ (A * options.x0(:) + epps + SinD),...
                            0,size(A,1),size(A,1)),2));
                    end
                end
            end
            if options.acosem
                if osa_iter > 1
                    if is_transposed
                        C_aco(:,osa_iter) = full(sum(spdiags(power(options.x0(:),1/options.h),0,size(A,1),size(A,1)) * A * ...
                            spdiags(uu ./ (A' * options.x0(:) + epps + SinD),0,size(A,2),size(A,2)),2));
                    else
                        C_aco(:,osa_iter) = full(sum(spdiags(power(options.x0(:),1/options.h),0,size(A,2),size(A,2)) * A' * ...
                            spdiags(uu ./ (A * options.x0(:) + epps + SinD),0,size(A,1),size(A,1)),2));
                    end
                end
            end
            if any(options.COSEM_MAP)
                if options.COSEM_MAP == 2
                    if osa_iter > 1
                        if is_transposed
                            C_osl(:,osa_iter) = full(sum(spdiags(options.x0(:),0,size(A,1),size(A,1)) * A * spdiags(uu ./ (A' * options.x0(:) + epps + SinD),...
                                0,size(A,2),size(A,2)),2));
                        else
                            C_osl(:,osa_iter) = full(sum(spdiags(options.x0(:),0,size(A,2),size(A,2)) * A' * spdiags(uu ./ (A * options.x0(:) + epps + SinD),...
                                0,size(A,1),size(A,1)),2));
                        end
                    end
                else
                    if osa_iter > 1
                        if is_transposed
                            C_osl(:,osa_iter) = full(sum(spdiags(power(options.x0(:),1/options.h),0,size(A,1),size(A,1)) * A * ...
                                spdiags(uu ./ (A' * options.x0(:) + epps + SinD),0,size(A,2),size(A,2)),2));
                        else
                            C_osl(:,osa_iter) = full(sum(spdiags(power(options.x0(:),1/options.h),0,size(A,2),size(A,2)) * A' * ...
                                spdiags(uu ./ (A * options.x0(:) + epps + SinD),0,size(A,1),size(A,1)),2));
                        end
                    end
                end
            end
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
                iRows = I(markers(1:end-1));
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
            options.lambda0 = lam;
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
    if (options.MBSREM || options.mramla) && options.implementation == 1
        options.pj3 = D/options.subsets;
    end
    % Compute the weights
    if (options.quad || options.L || options.FMH || options.weighted_mean || options.MRP || (options.TV && options.TVtype == 3)) && options.MAP
        if options.quad || options.L || options.FMH || options.weighted_mean || (options.TV && options.TVtype == 3)
            distX = options.FOVa_x/double(Nx);
            distY = options.FOVa_y/double(Ny);
            distZ = (double(options.axial_fov)/double(Nz));
            if isempty(options.weights)
                options.weights = zeros(((options.Ndx*2+1) * (options.Ndy*2+1) * (options.Ndz*2+1)),1);
                edist = zeros((options.Ndx*2+1),1);
                cc = zeros((options.Ndy*2+1)*(options.Ndx*2+1),1);
                lt = 0;
                for jj = options.Ndz : -1 : -options.Ndz
                    lt = lt + 1;
                    ll = 0;
                    for kk = options.Ndy : -1 : -options.Ndy
                        ll = ll + 1;
                        if options.Ndz == 0 || Nz == 1
                            if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                                apu = [((options.Ndx:-1:-options.Ndx) * distX)', (repeat_elem(kk,options.Ndy*2+1) * distY)];
                            elseif exist('OCTAVE_VERSION','builtin') == 5
                                apu = [((options.Ndx:-1:-options.Ndx) * distX)', (repelem(kk,options.Ndy*2+1) * distY)];
                            else
                                apu = [((options.Ndx:-1:-options.Ndx) * distX)', (repelem(kk,options.Ndy*2+1) * distY)'];
                            end
                        else
                            if options.Ndz ~= options.Ndx
                                if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                                    apu = [((options.Ndx:-1:-options.Ndx) * distX)', (repeat_elem(kk,options.Ndy*2+1) * distY), [zeros(options.Ndx-options.Ndz,1),(repeat_elem(jj,options.Ndz*2+1) * distZ),...
                                        zeros(options.Ndx-options.Ndz,1)]];
                                elseif exist('OCTAVE_VERSION','builtin') == 5
                                    apu = [((options.Ndx:-1:-options.Ndx) * distX)', (repelem(kk,options.Ndy*2+1) * distY), [zeros(options.Ndx-options.Ndz,1),(repelem(jj,options.Ndz*2+1) * distZ),...
                                        zeros(options.Ndx-options.Ndz,1)]];
                                else
                                    apu = [((options.Ndx:-1:-options.Ndx) * distX)', (repelem(kk,options.Ndy*2+1) * distY)', [zeros(options.Ndx-options.Ndz,1),(repelem(jj,options.Ndz*2+1) * distZ),...
                                        zeros(options.Ndx-options.Ndz,1)]'];
                                end
                            else
                                if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
                                    apu = [((options.Ndx:-1:-options.Ndx) * distX)', (repeat_elem(kk,options.Ndy*2+1) * distY), (repeat_elem(jj,options.Ndz*2+1) * distZ)];
                                elseif exist('OCTAVE_VERSION','builtin') == 5
                                    apu = [((options.Ndx:-1:-options.Ndx) * distX)', (repelem(kk,options.Ndy*2+1) * distY), (repelem(jj,options.Ndz*2+1) * distZ)];
                                else
                                    apu = [((options.Ndx:-1:-options.Ndx) * distX)', (repelem(kk,options.Ndy*2+1) * distY)', (repelem(jj,options.Ndz*2+1) * distZ)'];
                                end
                            end
                        end
                        for ii = 1 : length(apu)
                            edist(ii) = sqrt(apu(ii,:)*apu(ii,:)');
                        end
                        cc((options.Ndy*2+1)*(ll-1)+1:(options.Ndy*2+1)*ll) = edist;
                    end
                    options.weights((options.Ndx*2+1) * (options.Ndy*2+1) * (lt - 1) + 1: (options.Ndx*2+1) * (options.Ndy*2+1) * lt) = cc;
                end
                options.weights = 1./options.weights;
            end
        end
        %             pz_pad = padding(reshape(options.x0(:),Nx,Ny,Nz),[options.Ndx options.Ndy options.Ndz]);
        s = [Nx + options.Ndx*2 Ny + options.Ndy*2 Nz + options.Ndz*2];
        N_pad = min(3, options.Ndx + options.Ndy + options.Ndz);
        [c1{1:N_pad}]=ndgrid(1:(options.Ndx*2+1));
        c2(1:N_pad)={options.Ndy+1};
        if options.Ndz > options.Ndx && options.Ndz > 1
            c1{1} = cat(3, c1{1}, zeros(size(c1{1},1), size(c1{1},2), options.Ndz));
            c1{2} = cat(3, c1{2}, zeros(size(c1{2},1), size(c1{2},2), options.Ndz));
            c1{3} = cat(3, c1{3}, zeros(size(c1{3},1), size(c1{3},2), options.Ndz));
            %                 apu2 = c1{2};
            %                 apu3 = c1{3};
            for kk = options.Ndz - 1 : - 1 : 0
                %                     apu(:,:,end+1) = apu(:,:,end);
                %                     apu2(:,:,end+1) = apu2(:,:,end);
                %                     apu3(:,:,end+1) = apu3(:,:,end) + 1;
                c1{1}(:,:,end-kk) = c1{1}(:,:,end - kk - 1);
                c1{2}(:,:,end-kk) = c1{2}(:,:,end - kk - 1);
                c1{3}(:,:,end-kk) = c1{3}(:,:,end - kk - 1) + 1;
            end
            %                 c1{1} = apu;
            %                 c1{2} = apu2;
            %                 c1{3} = apu3;
            c2(end) = {options.Ndz+1};
        elseif options.Ndz < options.Ndx && options.Ndz > 1
            %                 apu = c1{1};
            %                 apu2 = c1{2};
            %                 apu3 = c1{3};
            c1{1}(:,:,end-2*(options.Ndx-options.Ndz) + 1) = [];
            c1{2}(:,:,end-2*(options.Ndx-options.Ndz) + 1) = [];
            c1{3}(:,:,end-2*(options.Ndx-options.Ndz) + 1) = [];
            %                 c1{1} = apu;
            %                 c1{2} = apu2;
            %                 c1{3} = apu3;
            c2(end) = {options.Ndz+1};
        end
        offsets=sub2ind(s,c1{:}) - sub2ind(s,c2{:});
        if Nz == 1
            tr_ind = sub2ind([Nx+options.Ndx*2 Ny + options.Ndy*2],mod((1:N)'-1,Nx)+(options.Ndx + 1),mod(floor(((1:double(N))'-1)/double(Nx)),double(Ny))+(options.Ndy + 1));
        else
            tr_ind = sub2ind([Nx+options.Ndx*2 Ny+options.Ndy*2 Nz+options.Ndz*2], mod((1:N)' - 1, Nx) + (options.Ndx + 1), mod(floor(((1:double(N))' - 1)/double(Nx)), ...
                double(Ny)) + (options.Ndy + 1), floor(((1:double(N))' - 1)/double(Nx * Ny)) + (options.Ndz + 1));
        end
        options.tr_offsets = uint32(bsxfun(@plus,tr_ind,offsets(:)'));
        if options.implementation == 2
            options.tr_offsets = options.tr_offsets - 1;
            options.options.Ndx = uint32(options.Ndx);
            options.options.Ndy = uint32(options.Ndy);
            options.options.Ndz = uint32(options.Ndz);
            clear tr_offsets
        end
        %             if (options.OSL_OSEM || options.OSL_MLEM) && options.quad
        %                 pz_pad_osl = pz_pad;
        %                 if options.implementation == 2
        %                     options.pz_pad_osl = single(pz_pad_osl);
        %                     clear pz_pad_osl
        %                 end
        %             end
        if options.quad || (options.TV && options.TVtype == 3)
            if options.empty_weight
                options.weights_quad = options.weights/sum(options.weights(~isinf(options.weights)));
                options.weights_quad = [options.weights_quad(1:floor(length(options.weights_quad) / 2)); ...
                    options.weights_quad(ceil(length(options.weights_quad)/2) + 1 : end)];
            else
                options.weights_quad = options.weights;
            end
            if options.implementation == 2
                options.weights_quad(isinf(options.weights_quad)) = [];
                options.weights_quad = single(options.weights_quad);
                clear weights_quad
            end
        end
        if options.L
            if isempty(options.a_L)
                options.a_L = lfilter_weights(options.Ndx, options.Ndy, options.Ndz, dx, dy, dz, options.oneD_weights);
            end
            if options.implementation == 2
                options.a_L = single(options.a_L);
                clear a_L
            end
            clear dd
        end
        if options.FMH
            if isempty(options.fmh_weights)
                kerroin = sqrt(2)*distX;
                if Nz == 1 || options.Ndz == 0
                    options.fmh_weights = zeros(options.Ndx*2+1, 4);
                    for jjj = 1 : 4
                        lll = lll + 1;
                        apu = zeros(options.Ndx*2+1,1);
                        hhh = 0;
                        if jjj == 1 || jjj == 3
                            for iii = options.Ndx : -1 : -options.Ndx
                                hhh = hhh + 1;
                                if iii == 0
                                    apu(hhh) = options.fmh_center_weight;
                                else
                                    apu(hhh) = kerroin/sqrt((distX*iii)^2+(distY*iii)^2);
                                end
                            end
                        elseif jjj == 2
                            for iii = options.Ndx : -1 : -options.Ndx
                                hhh = hhh + 1;
                                if iii == 0
                                    apu(hhh) = options.fmh_center_weight;
                                else
                                    apu(hhh) = kerroin/abs(distX*iii);
                                end
                            end
                        elseif jjj == 4
                            for iii = options.Ndx : -1 : -options.Ndx
                                hhh = hhh + 1;
                                if iii == 0
                                    apu(hhh) = options.fmh_center_weight;
                                else
                                    apu(hhh) = kerroin/abs(distY*iii);
                                end
                            end
                        end
                        options.fmh_weights(:, jjj) = apu;
                    end
                else
                    options.fmh_weights = zeros(max([options.Ndx*2+1,options.Ndz*2+1]), 13);
                    lll = 0;
                    for kkk = 1 : -1 : 0
                        for jjj = 1 : 9
                            lll = lll + 1;
                            if kkk == 1
                                apu = zeros(options.Ndz*2+1,1);
                                hhh = 0;
                                if jjj == 1 || jjj == 3 || jjj == 7 || jjj == 9
                                    for iii = options.Ndz : -1 : -options.Ndz
                                        hhh = hhh + 1;
                                        if iii == 0
                                            apu(hhh) = options.fmh_center_weight;
                                        else
                                            apu(hhh) = kerroin/sqrt(sqrt((distZ*iii)^2+(distX*iii)^2)^2+(distY*iii)^2);
                                        end
                                    end
                                elseif jjj == 2 || jjj == 8
                                    for iii = options.Ndz : -1 : -options.Ndz
                                        hhh = hhh + 1;
                                        if iii == 0
                                            apu(hhh) = options.fmh_center_weight;
                                        else
                                            apu(hhh) = kerroin/sqrt((distZ*iii)^2+(distX*iii)^2);
                                        end
                                    end
                                elseif jjj == 4 || jjj == 6
                                    for iii = options.Ndz : -1 : -options.Ndz
                                        hhh = hhh + 1;
                                        if iii == 0
                                            apu(hhh) = options.fmh_center_weight;
                                        else
                                            apu(hhh) = kerroin/sqrt((distZ*iii)^2+(distY*iii)^2);
                                        end
                                    end
                                elseif jjj == 5
                                    for iii = options.Ndz : -1 : -options.Ndz
                                        hhh = hhh + 1;
                                        if iii == 0
                                            apu(hhh) = options.fmh_center_weight;
                                        else
                                            apu(hhh) = kerroin/abs(distZ*iii);
                                        end
                                    end
                                end
                                options.fmh_weights(:, lll) = apu;
                            else
                                apu = zeros(options.Ndx*2+1,1);
                                hhh = 0;
                                if jjj == 1 || jjj == 3
                                    for iii = options.Ndx : -1 : -options.Ndx
                                        hhh = hhh + 1;
                                        if iii == 0
                                            apu(hhh) = options.fmh_center_weight;
                                        else
                                            apu(hhh) = kerroin/sqrt((distX*iii)^2+(distY*iii)^2);
                                        end
                                    end
                                elseif jjj == 2
                                    for iii = options.Ndx : -1 : -options.Ndx
                                        hhh = hhh + 1;
                                        if iii == 0
                                            apu(hhh) = options.fmh_center_weight;
                                        else
                                            apu(hhh) = kerroin/abs(distX*iii);
                                        end
                                    end
                                elseif jjj == 4
                                    for iii = options.Ndx : -1 : -options.Ndx
                                        hhh = hhh + 1;
                                        if iii == 0
                                            apu(hhh) = options.fmh_center_weight;
                                        else
                                            apu(hhh) = kerroin/abs(distY*iii);
                                        end
                                    end
                                else
                                    break
                                end
                                options.fmh_weights(:, lll) = apu;
                            end
                        end
                    end
                end
                options.fmh_weights = options.fmh_weights./sum(options.fmh_weights,1);
            end
            %                 if options.OSL_OSEM || options.OSL_MLEM
            %                     pz_pad_fmh = pz_pad;
            %                 end
            
            if options.implementation == 2
                options.fmh_weights = single(options.fmh_weights);
                clear fmh_weights pz_pad_fmh
            end
        end
        if (options.FMH || options.quad) && options.implementation == 2
            options.weights = single(options.weights);
            options.inffi = uint32(find(isinf(options.weights)) - 1);
        end
        if options.MRP
            options.medx = options.Ndx*2 + 1;
            options.medy = options.Ndy*2 + 1;
            options.medz = options.Ndz*2 + 1;
        end
        if options.weighted_mean
            if isempty(options.weighted_weights)
                kerroin = sqrt(2)*distX;
                options.weighted_weights = kerroin.*options.weights;
                options.weighted_weights(isinf(options.weighted_weights)) = options.weighted_center_weight;
                %                 options.weighted_weights = options.weighted_weights/sum(options.weighted_weights);
            end
            %                 if options.OSL_OSEM || options.OSL_MLEM
            %                     pz_pad_weighted = pz_pad;
            %                 end
            options.w_sum = sum(options.weighted_weights);
            if options.implementation == 2
                options.weighted_weights = single(options.weighted_weights);
                %                     if options.OSL_OSEM || options.OSL_MLEM
                %                         options.pz_pad_weighted = single(pz_pad_weighted);
                %                     end
                clear weighted_weights pz_pad_weighted
            end
        end
        clear tr_ind offsets c1 c2 apu apu2 apu3 N_pad cc pz_pad
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
    if (options.NLM && options.NLM_use_anatomical) && options.MAP
        apu = load(options.NLM_reference_image);
        variables = fields(apu);
        options.NLM_ref = double(apu.(variables{1}));
        options.NLM_ref = reshape(options.NLM_ref, Nx, Ny, Nz);
    end
end
end

