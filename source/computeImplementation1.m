function [A,ll, varargout] = computeImplementation1(options,use_raw_data,randoms_correction, pituus,osa_iter, normalization_correction,...
    Nx, Ny, Nz, dx, dy, dz, bx, by, bz, x, y, z_det, xx, yy, size_x, NSinos, NSlices, zmax, attenuation_correction, pseudot, det_per_ring, ...
    TOF, sigma_x, TOFCenter, dec, nCores, ind_size, block1, blocks, index, iij, jji, kkj, LL, N, summa, lor_a, xy_index, z_index, ...
    x_center, y_center, z_center, bmin, bmax, Vmax, V, lor_orth, gaussK,is_transposed, scatter_input, norm_input, SinD, koko)
%COMPUTEIMPLEMENTATION1 Computes the system matrix used by implementation
%1
%   Utility function
if options.precompute_lor == false && ~options.listmode
    if use_raw_data == false
        TOFSize = int64(pituus(osa_iter + 1) - pituus(osa_iter));
        % Siddon
        if options.projector_type == 1 || options.projector_type == 0
            if exist('OCTAVE_VERSION','builtin') == 0 && exist('projector_mex','file') == 3
                [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                    zmax, options.vaimennus, norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, ...
                    randoms_correction, options.scatter, scatter_input, options.global_correction_factor, uint16(0), uint32(0), uint16(0), NSinos, uint16(0), pseudot, det_per_ring, ...
                    TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, nCores, ...
                    use_raw_data, uint32(2), options.listmode, ind_size, block1, blocks, index(pituus(osa_iter)+1:pituus(osa_iter + 1)), ...
                    uint32(options.projector_type), iij, jji, kkj);
            elseif exist('OCTAVE_VERSION','builtin') == 5 && exist('projector_oct','file') == 3
                [ lor, indices, alkiot] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                    zmax, options.vaimennus, norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction, ...
                    randoms_correction, options.scatter, scatter_input, options.global_correction_factor, uint16(0), uint32(0), uint16(0), NSinos, uint16(0), pseudot, det_per_ring, ...
                    TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, nCores, ...
                    use_raw_data, uint32(2), options.listmode, ind_size, block1, blocks, index(pituus(osa_iter)+1:pituus(osa_iter + 1)), ...
                    uint32(options.projector_type), iij, jji, kkj);
            else
                % The below lines allow for pure MATLAB
                % implemention, i.e. no MEX-files will be
                % used. Currently the below function
                % uses parfor-loops (requires parallel
                % computing toolbox).
                % NOTE: The below method is not
                % recommended since it is much slower
                % method.
                [ lor, indices, alkiot, ~] = improved_siddon_atten( int32(Ny), int32(Nx), int32(Nz), dx, dz, by, bx, bz, z_det, x, y, yy, xx, ...
                    NSinos, NSlices, options.vaimennus, index(pituus(osa_iter)+1:pituus(osa_iter + 1)), pituus(osa_iter + 1) - pituus(osa_iter), ...
                    attenuation_correction);
                alkiot = cell2mat(alkiot);
                indices = cell2mat(indices) - 1;
                lor = lor(:,2);
            end
            % Orthogonal distance based
        elseif options.projector_type == 2
            error('Unsupported projector type')
        else
            error('Unsupported projector type')
        end
    else
        L = LL(pituus(osa_iter) * 2 + 1 : pituus(osa_iter + 1) * 2);
        TOFSize = int64(size(L,1));
        if options.projector_type == 1
            if exist('OCTAVE_VERSION','builtin') == 0
                [ lor, indices, alkiot] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                    zmax, options.vaimennus, norm_input, SinD, uint32(0), attenuation_correction, normalization_correction, ...
                    randoms_correction, options.scatter, scatter_input, options.global_correction_factor, uint16(0), uint32(0), uint16(0), NSinos, L, pseudot, det_per_ring, ...
                    TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, nCores, ...
                    use_raw_data, uint32(2), options.listmode, ind_size, block1, blocks, uint32(0), uint32(options.projector_type), iij, jji, kkj);
            elseif exist('OCTAVE_VERSION','builtin') == 5
                [ lor, indices, alkiot] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, ...
                    zmax, options.vaimennus, norm_input, SinD, uint32(0), attenuation_correction, normalization_correction, ...
                    randoms_correction, options.scatter, scatter_input, options.global_correction_factor, uint16(0), uint32(0), uint16(0), NSinos, L, pseudot, det_per_ring, ...
                    TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, nCores, ...
                    use_raw_data, uint32(2), options.listmode, ind_size, block1, blocks, uint32(0), uint32(options.projector_type), iij, jji, kkj);
            end
        elseif options.projector_type == 2
            error('Unsupported projector type')
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
    
    A_length = koko;
    if options.verbose
        tStart = tic;
    end
    % Form the sparse matrix
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
    ll = [];
    clear indices alkiot lor
    if options.verbose
        tElapsed = toc(tStart);
        disp(['Sparse matrix formation took ' num2str(tElapsed) ' seconds'])
    end
    % Precomputation performed
    % Parallel
    % Faster
    % Only C++ code (no pure MATLAB implementation)
else
    if use_raw_data
        L_input = LL(pituus(osa_iter) * 2 + 1 : pituus(osa_iter + 1) * 2);
        xy_index_input = uint32(0);
        z_index_input = uint16(0);
        TOFSize = int64(size(L_input,1));
    else
        L_input = uint16(0);
        xy_index_input = xy_index(pituus(osa_iter)+1:pituus(osa_iter + 1));
        z_index_input = z_index(pituus(osa_iter)+1:pituus(osa_iter + 1));
        TOFSize = int64(size(xy_index_input,1));
    end
    if options.projector_type == 2 || options.projector_type == 3
        if exist('OCTAVE_VERSION','builtin') == 0
            lor2 = [uint64(0);cumsum(uint64(lor_orth(pituus(osa_iter)+1:pituus(osa_iter + 1))))];
        else
            lor2 = [uint64(0);cumsum(uint64(lor_orth(pituus(osa_iter)+1:pituus(osa_iter + 1))),'native')];
        end
    else
        if exist('OCTAVE_VERSION','builtin') == 0
            lor2 = [uint64(0);cumsum(uint64(lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1))))];
        else
            lor2 = [uint64(0);cumsum(uint64(lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1))),'native')];
        end
    end
    if exist('OCTAVE_VERSION','builtin') == 0
        [A, ll] = projector_mex( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
            norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction,...
            randoms_correction, options.scatter, scatter_input, options.global_correction_factor, lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1)), xy_index_input, ...
            z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
            TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, nCores, use_raw_data, uint32(0), options.listmode, lor2, summa(osa_iter), ...
            options.attenuation_phase, uint32(options.projector_type), options.tube_width_xy, x_center, y_center, z_center, ...
            options.tube_width_z, bmin, bmax, Vmax, V);
    elseif exist('OCTAVE_VERSION','builtin') == 5
        [A, ll] = projector_oct( Ny, Nx, Nz, dx, dz, by, bx, bz, z_det, x, y, dy, yy, xx , NSinos, NSlices, size_x, zmax, options.vaimennus, ...
            norm_input, SinD, pituus(osa_iter + 1) - pituus(osa_iter), attenuation_correction, normalization_correction,...
            randoms_correction, options.scatter, scatter_input, options.global_correction_factor, lor_a(pituus(osa_iter)+1:pituus(osa_iter + 1)), xy_index_input, ...
            z_index_input, NSinos, L_input, pseudot, det_per_ring, ...
            TOF, TOFSize, sigma_x, TOFCenter, int64(options.TOF_bins), dec, options.verbose, nCores, use_raw_data, uint32(0), options.listmode, lor2, summa(osa_iter), ...
            options.attenuation_phase, uint32(options.projector_type), options.tube_width_xy, x_center, y_center, z_center, ...
            options.tube_width_z, bmin, bmax, Vmax, V);
    end
    clear lor2
end
if nargout >= 3
    % Sensitivity image
    if is_transposed
        varargout{1} = full(sum(A,2));
    else
        varargout{1} = full(sum(A,1))';
    end
    varargout{1}(varargout{1} < options.epps) = options.epps;
    if options.use_psf
        varargout{1} = computeConvolution(varargout{1}, options, Nx, Ny, Nz, gaussK);
    end
end
end

