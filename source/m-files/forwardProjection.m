function [outputFP, A] = forwardProjection(options, recApu, x, z, koko, nMeas, xy_index, z_index, norm_input, corr_input, L_input, TOF, lor2, lor1, summa, loopVar, subIter, varargin)
% FORWARDPROJECTION Computes the forward projection

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2023-2024 Ville-Veikko Wettenhovi
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

% Detector pitch for CT and SPECT, crystal size (transaxial, axial) for PET. The multi-ray shifts in the
% projectors use these as the transaxial and axial extents, whereas options.dPitchX and options.dPitchY hold (cr_pz, cr_p) for PET
if options.CT || options.SPECT
    dPitchXY = options.dPitchX;
    dPitchZ = options.dPitchY;
else
    dPitchXY = options.cr_p;
    dPitchZ = options.cr_pz;
end
if isempty(varargin)
    nCores = 0;
else
    nCores = varargin{1};
end
if ismac && options.use_64bit_atomics
    options.use_64bit_atomics = false;
end
if iscell(recApu)
    useCell = true;
else
    useCell = false;
end
projType = options.projector_type;
if projType == 1 || (projType >= 10 && projType < 20)
    projType = 1;
end
if options.additionalCorrection && isempty(corr_input)
    error('Additional correction selected, but no data inserted!')
end
if options.normalization_correction
    if isempty(norm_input)
        error('Normalization correction selected, but no normalization data inserted! Insert normalization coefficients to param.normalization of the class object')
    end
end
if options.attenuation_correction
    if ~isfield(options,'vaimennus') || numel(options.vaimennus) <= 1
        error('Attenuation correction selected, but no attenuation data inserted! Insert attenuation coefficients to param.atten of the class object')
    end
end

outputFP = [];
A = [];

function outputFP = forwardProjectionType6(options, recApu, loopVar, koko)
    % Matches the GPU forwardProjectionType6 (functions.hpp): GPU rotate kernel convention, integer (rounded) zero-fill
    % translations, the attenuation map is rotated and translated with the emission image
    outputFP = zeros(options.nRowsD, options.nColsD, koko, options.cType);
    for ii = loopVar
        u1 = options.uu;
        apuArr = reshape(recApu, options.Nx(ii), options.Ny(ii), options.Nz(ii));

        for kk = 1:koko % todo add support for non-square fov
            % Rounded integer shifts, as the int32 distInt/distInt2 of the GPU code
            blurPlane = round(double(options.blurPlanes(u1)));
            panelShift = round(double(options.blurPlanes2(u1)));

            % 1. Rotate the image
            kuvaRot = rotateType6(apuArr, -options.swivelAngles(u1));

            % 2. Translate the image
            if panelShift ~= 0
                kuvaRot = translateType6(kuvaRot, 0, -panelShift, 0);
            end

            % 3. Attenuation correction (rotation, translation, accumulation and exponentiation of the map)
            if options.attenuation_correction
                attenuationImage = attenuationType6(options.vaimennus, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.swivelAngles(u1), panelShift);
                % 3.1 pointwise multiply with image volume
                kuvaRot = kuvaRot .* attenuationImage;
            end

            % 4. Convolve with detector PSF
            kuvaRot = convolvePSFType6(kuvaRot, options.gFilter, blurPlane, options.Nx(ii), options.Ny(ii), options.Nz(ii));

            % 5. Sum
            kuvaRot = reshape(sum(kuvaRot, 1), options.Ny(ii), options.Nz(ii));
            kuvaRot = kuvaRot / double(options.Nx(ii));
            outputFP(:, :, kk) = outputFP(:, :, kk) + kuvaRot;
            u1 = u1 + 1;
        end
    end
end

if iscell(recApu)
    if ~isempty(recApu{1}) && numel(recApu{1}) ~= options.Nx(1) * options.Ny(1) * options.Nz(1)
        error('Input image has incorrect dimensions! Must equal to Nx*Ny*Nz!')
    end
else
    if ~isempty(recApu) && numel(recApu) ~= options.Nx(1) * options.Ny(1) * options.Nz(1)
        error('Input image has incorrect dimensions! Must equal to Nx*Ny*Nz!')
    end
end
if (~ismac && (options.implementation == 4 || options.implementation == 1)) || (ismac && options.projector_type == 6)
    if options.projector_type == 6
        outputFP = forwardProjectionType6(options, recApu, loopVar, koko);
        A = [];
        outputFP = outputFP(:);
    else
        % The CPU MEX files read the attenuation image without type conversion
        if isfield(options, 'vaimennus') && isnumeric(options.vaimennus)
            if options.implementation == 4 && options.useSingles
                if ~isa(options.vaimennus, 'single')
                    options.vaimennus = single(options.vaimennus);
                end
            elseif ~isa(options.vaimennus, 'double')
                options.vaimennus = double(options.vaimennus);
            end
        end
        for ii = loopVar
            if ii == 1 || isscalar(loopVar)
                if useCell
                    if options.use_psf && ~isempty(recApu{ii})
                        recApu{ii} = computeConvolution(recApu{ii}, options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.gaussK);
                    end
                    if options.implementation == 4
                        if options.useSingles
                            [outputFP, A] = projector_mexSingle( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ...
                                options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...
                                options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ...
                                options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, recApu{ii}, true, 1, dPitchXY, ...
                                dPitchZ, options.nProjections);
                        else
                            [outputFP, A] = projector_mex( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ...
                                options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...
                                options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ...
                                options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, recApu{ii}, true, 1, dPitchXY, ...
                                dPitchZ, options.nProjections);
                        end
                    else
                        error('Unsupported implementation! This is unintended behavior!')
                    end
                else
                    if options.use_psf && ~isempty(recApu)
                        recApu = computeConvolution(recApu, options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.gaussK);
                    end
                    if options.implementation == 1
                        [A] = projector_mex( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ... % 8
                            options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...% 18
                            options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ... % 28
                            options.TOF_bins, options.verbose, nCores, options.use_raw_data, 0, options.listmode, projType, subIter, nMeas, lor2, lor1, summa, dPitchZ, options.nProjections); % 42
                        if ~isempty(recApu)
                            outputFP = A' * recApu;
                        else
                            outputFP = [];
                        end
                    elseif options.implementation == 4
                        if options.useSingles
                            [outputFP, A] = projector_mexSingle( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ... % 8
                                options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ... % 18
                                options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ... % 28
                                options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, recApu, true, 1, dPitchXY, ... % 40
                                dPitchZ, options.nProjections); % 43
                        else
                            [outputFP, A] = projector_mex( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ... % 8
                                options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ... % 18
                                options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ... % 28
                                options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, recApu, true, 1, dPitchXY, ... % 40
                                dPitchZ, options.nProjections); % 43
                        end
                    else
                        error('Unsupported implementation! This is unintended behavior!')
                    end
                end
            else
                if useCell
                    if options.use_psf && ~isempty(recApu{ii})
                        recApu{ii} = computeConvolution(recApu{ii}, options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.gaussK);
                    end
                    if options.implementation == 4
                        if options.useSingles
                            [apu, A] = projector_mexSingle( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ...
                                options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...
                                options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ...
                                options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, recApu{ii}, true, 1, dPitchXY, ...
                                dPitchZ, options.nProjections);
                        else
                            [apu, A] = projector_mex( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ...
                                options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...
                                options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ...
                                options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, recApu{ii}, true, 1, dPitchXY, ...
                                dPitchZ, options.nProjections);
                        end
                    else
                        error('Unsupported implementation! This is unintended behavior!')
                    end
                else
                    if options.use_psf && ~isempty(recApu)
                        recApu = computeConvolution(recApu, options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.gaussK);
                    end
                    if options.implementation == 4
                        if options.useSingles
                            [apu, A] = projector_mexSingle( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ...
                                options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...
                                options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ...
                                options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, recApu, true, 1, dPitchXY, ...
                                dPitchZ, options.nProjections);
                        else
                            [apu, A] = projector_mex( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ...
                                options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...
                                options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ...
                                options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, recApu, true, 1, dPitchXY, ...
                                dPitchZ, options.nProjections);
                        end
                    else
                        error('Unsupported implementation! This is unintended behavior!')
                    end
                end
                outputFP = outputFP + apu;
            end
        end
    end
elseif options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || ismac
    A = [];
    % Determined here, before recApu is reshaped/converted below, so that
    % the gpuArray-ness of the original input is not lost
    useGPUArray = isa(recApu, 'gpuArray') || (iscell(recApu) && ~isempty(recApu) && isa(recApu{1}, 'gpuArray'));
    if useGPUArray
        if options.implementation ~= 5
            error('gpuArray input is only supported with implementation 5!')
        end
        if exist('CUDA_matrixfree_multi_gpu','file') ~= 3
            error(['gpuArray input was used, but the CUDA build of implementation 5 (CUDA_matrixfree_multi_gpu) was not found. ' ...
                'Run install_mex to build it. gpuArray input requires a CUDA-enabled installation and the CUDA toolkit.'])
        end
    end
    if ~isfield(options,'orthTransaxial') && (options.projector_type == 2 || options.projector_type == 3 || options.projector_type == 22 || options.projector_type == 33)
        if options.projector_type == 3 || options.projector_type == 33
            options.orthTransaxial = true;
        elseif (options.projector_type == 2 || options.projector_type == 22) && (isfield(options,'tube_width_xy') && options.tube_width_xy > 0 || options.SPECT)
            options.orthTransaxial = true;
        else
            options.orthTransaxial = false;
        end
    end
    if ~isfield(options,'orthAxial') && (options.projector_type == 2 || options.projector_type == 3 || options.projector_type == 22 || options.projector_type == 33)
        if options.projector_type == 3 || options.projector_type == 33
            options.orthAxial = true;
        elseif (options.projector_type == 2 || options.projector_type == 22) && (isfield(obj.param,'tube_width_z') && obj.param.tube_width_z > 0 || options.SPECT)
            options.orthAxial = true;
        else
            options.orthAxial = false;
        end
    end
    if options.use_32bit_atomics && options.use_64bit_atomics
        options.use_32bit_atomics = false;
    end
    if options.use_raw_data
        TOFSize = int64(size(options.LL,1));
    else
        TOFSize = int64(size(options.xy_index,1));
    end
    n_rays = uint16(options.n_rays_transaxial);
    n_rays3D = uint16(options.n_rays_axial);
    if ~isfield(options, 'vaimennus')
        options.vaimennus = single(0);
    end
    if iscell(recApu)
        if options.use_psf
            for kk = 1 : size(recApu,1)
                recApu{kk} = computeConvolution(recApu{kk}, options, options.Nx(kk), options.Ny(kk), options.Nz(kk), options.gaussK);
            end
        end
        if options.projector_type == 5 || options.projector_type == 55 || options.projector_type == 51 || options.projector_type == 54
            kopio = recApu;
            meanFPX = cell(size(recApu,1), 1);
            meanFPY = cell(size(recApu,1), 1);
            for kk = 1 : size(recApu,1)
                recApu{kk} = reshape(recApu{kk}, options.Nx(kk), options.Ny(kk), options.Nz(kk));
                kopio{kk} = reshape(kopio{kk}, options.Nx(kk), options.Ny(kk), options.Nz(kk));
                kopio{kk} = permute(kopio{kk}, [1 3 2]);
                % Matches the single-volume branch below (permute(recApu,[2 3 1]) before
                % computeIntegralImage): without this, recApu{kk} stays in (Nx,Ny,Nz)
                % orientation and its integral image/mean are computed along the wrong
                % axes instead of producing the YZ image (Ny+1,Nz+1,Nx) the C++ side
                % (multi_gpu_reconstruction.h, d_image_os_int) expects.
                recApu{kk} = permute(recApu{kk}, [2 3 1]);
                [recApu{kk}, meanFPX{kk}] = computeIntegralImage(recApu{kk}, options.meanFP);
                [kopio{kk}, meanFPY{kk}] = computeIntegralImage(kopio{kk}, options.meanFP);
            end
        end
        % cell2mat does not support gpuArray cells ("CELL2MAT does not support
        % cell arrays containing cell arrays or objects"), so vertcat(c{:}) is
        % used instead. This works identically to cell2mat here because every
        % cell element is already a column vector (recApu{kk}/kopio{kk} come
        % out of computeIntegralImage as cast(output(:), type), and
        % meanFPBlocks{kk} is built with (:) below), and it works for both
        % plain double/single arrays and gpuArray.
        options.x0 = vertcat(recApu{:});
        if options.projector_type == 5 || options.projector_type == 55 || options.projector_type == 51 || options.projector_type == 54
            options.x0 = [options.x0;vertcat(kopio{:})];
            if options.meanFP
                % Per volume: [Nx(kk) means of the YZ integral image; Ny(kk) means of the XZ integral image],
                % concatenated volume by volume (matches the C++ implementation-5 per-volume upload order)
                meanFPBlocks = cell(size(recApu,1), 1);
                for kk = 1 : size(recApu,1)
                    meanFPBlocks{kk} = [meanFPX{kk}(:); meanFPY{kk}(:)];
                end
                options.meanFPvec = single(vertcat(meanFPBlocks{:}));
                % meanFPvec is read by the MEX from the options struct via the
                % regular (host) mxArray API, so it must be a host array even
                % when x0 stays on the device: gatherHostFields below gathers
                % every gpuArray field except 'x0', which covers meanFPvec.
            end
        end
    else
        if options.use_psf
            recApu = computeConvolution(recApu, options, options.Nx(1), options.Ny(1), options.Nz(1), options.gaussK);
        end
        if options.projector_type == 5 || options.projector_type == 55 || options.projector_type == 51 || options.projector_type == 54
            recApu = reshape(recApu, options.Nx(1), options.Ny(1), options.Nz(1));
            kopio = permute(recApu, [1 3 2]);
            recApu = permute(recApu, [2 3 1]);
            [recApu, meanFPX] = computeIntegralImage(recApu, options.meanFP);
            [kopio, meanFPY] = computeIntegralImage(kopio, options.meanFP);
            recApu = [recApu(:);kopio(:)];
            if options.meanFP
                % [Nx means of the YZ integral image (d_image_os_int); Ny means of the XZ integral image (d_image_os)]
                options.meanFPvec = single([meanFPX(:); meanFPY(:)]);
            end
        end
        options.x0 = recApu;
    end
    if ~isa(options.x0,'single')
        options.x0 = single(options.x0);
    end
    options.use_device = uint32(options.use_device);

    if options.orthAxial
        crystal_size_z = (options.tube_width_z);
    else
        crystal_size_z = (options.tube_width_xy);
    end
    if options.projector_type == 1 || options.projector_type == 11 ...
            || options.projector_type == 2 || options.projector_type == 3 || options.projector_type == 22 || options.projector_type == 33
        kernel_file = 'projectorType123.cl';
        kernel_path = which(kernel_file);
        kernel_path = strrep(kernel_path, '\', '/');
        kernel_path = strrep(kernel_path, '.cl', '');
        header_directory = strrep(kernel_path,'projectorType123','');
    elseif options.projector_type == 4 || options.projector_type == 44 || options.projector_type == 41 || options.projector_type == 14 || options.projector_type == 45
        kernel_file = 'projectorType4.cl';
        kernel_path = which(kernel_file);
        kernel_path = strrep(kernel_path, '\', '/');
        kernel_path = strrep(kernel_path, '.cl', '');
        header_directory = strrep(kernel_path,'projectorType4','');
    elseif options.projector_type == 5 || options.projector_type == 55 || options.projector_type == 51 || options.projector_type == 15 || options.projector_type == 54
        kernel_file = 'projectorType5.cl';
        kernel_path = which(kernel_file);
        kernel_path = strrep(kernel_path, '\', '/');
        kernel_path = strrep(kernel_path, '.cl', '');
        header_directory = strrep(kernel_path,'projectorType5','');
    elseif options.projector_type == 7 || options.projector_type == 71 || options.projector_type == 17 || options.projector_type == 74 ...
            || options.projector_type == 47 || options.projector_type == 75 || options.projector_type == 57
        kernel_file = 'projectorType7.cl';
        kernel_path = which(kernel_file);
        kernel_path = strrep(kernel_path, '\', '/');
        kernel_path = strrep(kernel_path, '.cl', '');
        header_directory = strrep(kernel_path,'projectorType7','');
    elseif options.projector_type == 6
        kernel_file = 'auxKernels.cl';
        kernel_path = which(kernel_file);
        kernel_path = strrep(kernel_path, '\', '/');
        kernel_path = strrep(kernel_path, '.cl', '');
        header_directory = strrep(kernel_path,'auxKernels','');
    else
        error('Invalid projector for OpenCL')
    end
    if numel(nMeas) == 1
        nMeas = [0;nMeas];
    end
    if useGPUArray
        options = gatherHostFields(options, {'x0'});
        mexProjector = @CUDA_matrixfree_multi_gpu;
    else
        mexProjector = @OpenCL_matrixfree_multi_gpu;
    end
    [outputFP] = mexProjector( options.Nx, options.Ny, options.Nz, options.dx, options.dy, options.dz, options.bx, options.by, options.bz, ...
        z, x, options.nRowsD, options.verbose, options.LL, options.TOF, ... % 15
        TOFSize, options.sigma_x, options.TOFCenter, options.TOF_bins, options.platform, options.use_raw_data, options.use_psf, header_directory, options.vaimennus, ... % 24
        options.normalization, nMeas, options.attenuation_correction, options.normalization_correction, 1, options.subsets, options.epps, options.xy_index, ...
        options.z_index, crystal_size_z, ... % 34
        options.x_center, options.y_center, options.z_center, single(0), 0, options.projector_type, n_rays, n_rays3D, ... % 42
        options, single(0), options.partitions, options.use_64bit_atomics, options.bmin, options.bmax, options.Vmax, options.V, options.gaussK, 1, 1); % 51
    end
end

function out = rotateType6(im, angleDeg)
% Rotates every axial plane (dimensions 1 and 2) of the 3D array im by angleDeg degrees. This is a port of the
% GPU rotate kernel (auxKernels.cl, ROTATE) used by the projector type 6: rotation about the image center
% (floor(N/2) - 0.5 in zero-based pixel coordinates), bilinear interpolation, neighbors clamped to the edge and
% zero for destination pixels whose source lies outside the image. Note that this is NOT the convention of
% imrotate: a rotation by angleDeg here equals imrotate with -angleDeg
[Nx, Ny, Nz] = size(im);
% The GPU kernel (rotateCustomAF in functions.hpp) rotates by -angle
angleRad = -angleDeg * pi / 180;
cosa = cos(angleRad);
sina = sin(angleRad);
hx = floor(Nx / 2);
hy = floor(Ny / 2);
[ix, iy] = ndgrid(0 : Nx - 1, 0 : Ny - 1);
xA = ix - hx + 0.5;
yA = iy - hy + 0.5;
srcX = xA * cosa - yA * sina + hx - 0.5;
srcY = xA * sina + yA * cosa + hy - 0.5;
valid = srcX >= 0 & srcX < Nx & srcY >= 0 & srcY < Ny;
x0 = floor(srcX);
y0 = floor(srcY);
sx = srcX - x0;
sy = srcY - y0;
x1 = min(max(x0 + 1, 0), Nx - 1);
y1 = min(max(y0 + 1, 0), Ny - 1);
x0 = min(max(x0, 0), Nx - 1);
y0 = min(max(y0, 0), Ny - 1);
imFlat = reshape(im, Nx * Ny, Nz);
out = (1 - sx(:)) .* (1 - sy(:)) .* imFlat(x0(:) + y0(:) * Nx + 1, :) ...
    + sx(:) .* (1 - sy(:)) .* imFlat(x1(:) + y0(:) * Nx + 1, :) ...
    + (1 - sx(:)) .* sy(:) .* imFlat(x0(:) + y1(:) * Nx + 1, :) ...
    + sx(:) .* sy(:) .* imFlat(x1(:) + y1(:) * Nx + 1, :);
out(~valid(:), :) = 0;
out = cast(reshape(out, Nx, Ny, Nz), 'like', im);
end

function out = translateType6(im, shift1, shift2, shift3)
% Integer translation with zero fill (circular shift followed by zeroing of the wrapped elements), same as the
% translateHelperType6 of the GPU code. Positive shifts move the data toward higher indices
shifts = [shift1 shift2 shift3];
sz = [size(im, 1) size(im, 2) size(im, 3)];
if any(abs(shifts) >= sz & shifts ~= 0)
    out = zeros(size(im), 'like', im);
    return
end
out = circshift(im, shifts(1 : ndims(im)));
if shift1 > 0
    out(1 : shift1, :, :) = 0;
elseif shift1 < 0
    out(end + shift1 + 1 : end, :, :) = 0;
end
if shift2 > 0
    out(:, 1 : shift2, :) = 0;
elseif shift2 < 0
    out(:, end + shift2 + 1 : end, :) = 0;
end
if shift3 > 0
    out(:, :, 1 : shift3) = 0;
elseif shift3 < 0
    out(:, :, end + shift3 + 1 : end) = 0;
end
end

function attenuationImage = attenuationType6(vaimennus, Nx, Ny, Nz, dx, swivelAngle, shift2)
% Attenuation factors of the projector type 6: the map is rotated and translated with the emission image,
% accumulated along the projection direction (dimension 1) and exponentiated
attenuationImage = reshape(vaimennus, Nx, Ny, Nz);
attenuationImage = rotateType6(attenuationImage, -swivelAngle);
if shift2 ~= 0
    attenuationImage = translateType6(attenuationImage, 0, -shift2, 0);
end
attenuationImage = cumsum(attenuationImage, 1);
attenuationImage = exp(-dx * attenuationImage);
end

function kuvaRot = convolvePSFType6(kuvaRot, gFilter, shift3, Nx, Ny, Nz)
% Convolves each x-plane of the volume with the detector PSF plane of the same index. The PSF stack is shifted
% (zero fill) along its third dimension by the rounded blur plane value
PSF = translateType6(gFilter, 0, 0, shift3);
PSF = PSF(:, :, 1 : Nx);
PSF = permute(PSF, [3 2 1]);
for jj = 1 : Nx
    kuvaRot(jj, :, :) = reshape(conv2(reshape(kuvaRot(jj, :, :), Ny, Nz), reshape(PSF(jj, :, :), size(PSF, 2), size(PSF, 3)), 'same'), 1, Ny, Nz);
end
end
