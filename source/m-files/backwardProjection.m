function [output, sensIm] = backwardProjection(options, input, x, z, koko, nMeas, xy_index, z_index, norm_input, corr_input, L_input, TOF, noSensIm, subIter, varargin)
% BACKWARDPROJECTION Computes the backprojection

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
inputCell = false;
if options.nMultiVolumes > 0
    outputCell = true;
    output = cell(options.nMultiVolumes + 1, 1);
    sensIm = cell(options.nMultiVolumes + 1, 1);
    if iscell(input)
        inputCell = true;
    end
else
    outputCell = false;
end
projType = options.projector_type;
if projType == 1 || projType == 11
    projType = 1;
end
listmodeSPECTSensitivity = options.listmode > 0 && options.compute_sensitivity_image && options.SPECT;
if listmodeSPECTSensitivity
	if isfield(options, 'sensitivityViewWeights') && ~isempty(options.sensitivityViewWeights)
		if isfield(options, 'xSens') && isfield(options, 'zSens') && ...
				~isempty(options.xSens) && ~isempty(options.zSens)
			x = options.xSens;
			z = options.zSens;
		else
			[x, z] = get_coordinates_SPECT(options);
		end
	else
		[x, z, options.sensitivityViewWeights] = prepareSPECTListmodeSensitivity(options);
	end
elseif options.listmode > 0 && options.compute_sensitivity_image && isfield(options,'xSens') && isfield(options,'zSens') && ~isempty(options.xSens) && ~isempty(options.zSens)
    x = options.xSens;
    z_det = options.zSens;
	z = z_det;
	if ~options.SPECT
		options.det_per_ring = numel(x) / 2 / options.nLayers; % x and z contain the coordinates of all layers
		options.rings = numel(z_det) / options.nLayers;
	end
elseif options.listmode > 0 && options.compute_sensitivity_image && ~options.useIndexBasedReconstruction
	if options.SPECT
		[x, z] = get_coordinates_SPECT(options);
	else
		options.use_raw_data = true;
		[x, ~, z, ~] = get_coordinates(options);
		options.use_raw_data = false;
	end
end

function output = backwardProjectionType6(options, input, koko, norm_input, corr_input, nMeas, subIter)
    % Implementations that do not split subsets in projectorClass pass the full
    % measurement boundaries here. Restrict the input and view weights to the
    % requested subset; options.ub already points at its first view.
    nRowsProj = options.nRowsD * options.nColsD;
    if numel(nMeas) > 1 && subIter >= 0 && subIter + 2 <= numel(nMeas)
        firstView = double(nMeas(subIter + 1));
        lastView = double(nMeas(subIter + 2));
        subsetViews = lastView - firstView;
        totalViews = double(nMeas(end));
        if numel(input) == nRowsProj * totalViews
            firstElement = firstView * nRowsProj + 1;
            lastElement = lastView * nRowsProj;
            input = input(firstElement:lastElement);
            if numel(norm_input) == nRowsProj * totalViews
                norm_input = norm_input(firstElement:lastElement);
            end
            if numel(corr_input) == nRowsProj * totalViews
                corr_input = corr_input(firstElement:lastElement);
            end
        end
        if numel(input) == nRowsProj * subsetViews
            koko = subsetViews;
        end
    end
    fProj = reshape(input, options.nRowsD, options.nColsD, koko);
    viewIndices = options.ub:(options.ub + koko - 1);
    fProj = fProj .* type6CentralRayWeights(options, viewIndices);
    fProj = fProj .* type6MeasurementWeights(options, norm_input, corr_input, viewIndices);
    nVolumes = double(options.nMultiVolumes) + 1;
    if options.nMultiVolumes > 0
        output = cell(nVolumes, 1);
    else
        output = zeros(prod(double([options.Nx(1), options.Ny(1), options.Nz(1)])), 1, options.cType);
    end
    for volume = 1:nVolumes
        [geom, gridSize] = type6NativeOperator('volume-geometry', options, volume);
        gridNx = gridSize(1); gridNy = gridSize(2);
        [gridNz, padLo, padHi] = type6NativeOperator('axial-geometry', options, volume);
        transA = double(options.dy(1)) / double(options.dy(volume));
        transB = 0.5 * gridNy + (0.5 - 0.5 * options.nRowsD) * transA - 0.5;
        axialA = double(options.dz(1)) / double(options.dz(volume));
        axialB = -double(options.bz(volume)) / double(options.dz(volume)) + ...
            (0.5 - 0.5 * options.nColsD) * axialA - 0.5 + padLo;
        volumeOutput = zeros(prod(double([options.Nx(volume), options.Ny(volume), options.Nz(volume)])), 1, options.cType);
        u1 = options.ub;
        for kk = 1:koko
            panelShift = double(options.blurPlanes2F(u1)) * double(options.dx(1)) / double(options.dx(volume));
            detector = fProj(:, :, kk);
            detector = type6NativeOperator('detector-to-grid-adjoint', detector, gridNy, transA, transB);
            detector = type6NativeOperator('detector-to-grid-adjoint', detector.', gridNz, axialA, axialB).';
            kuvaRot = repmat(reshape(detector, 1, gridNy, gridNz), gridNx, 1, 1);
            kuvaRot = kuvaRot * (options.dx(volume) / options.dx(1)) / (options.dy(1) * options.dz(1));

            if isfield(options, 'CDRFOnTheFly') && options.CDRFOnTheFly
                PSF = type6NativeOperator('cdrf', options, volume, u1);
            elseif ~isempty(options.gFilter)
                PSF = imtranslate(options.gFilter, [0, 0, options.blurPlanes(u1)], FillValues=0);
                PSF = PSF(:, :, 1:gridNx);
                PSF = permute(PSF, [3 2 1]);
            else
                PSF = [];
            end
            if ~isempty(PSF)
                for depth = 1:gridNx
                    kuvaRot(depth, :, :) = type6NativeOperator('conv-adjoint', ...
                        squeeze(kuvaRot(depth, :, :)), squeeze(PSF(depth,:,:)));
                end
            end

            if options.attenuation_correction
                attenuationImage = type6NativeOperator('attenuation-factors', options, volume, u1, gridSize, gridNz, panelShift);
                kuvaRot = kuvaRot .* attenuationImage;
            end

            [offX, offY] = type6NativeOperator('rotation-parameters', options, volume, u1, ...
                gridSize, gridNz, padLo, padHi);
            kuvaRot = type6NativeOperator('type6-rotate-adjoint', kuvaRot, ...
                -options.swivelAngles(u1), offX, offY, geom);
            kuvaRot = kuvaRot(:, :, padLo + (1:double(options.Nz(volume))));
            volumeOutput = volumeOutput + kuvaRot(:);
            u1 = u1 + 1;
        end
        if options.nMultiVolumes > 0
            output{volume} = volumeOutput;
        else
            output = volumeOutput;
        end
    end
end

if options.projector_type == 6 || (options.implementation == 4 && ~ismac)
    if options.projector_type == 6
        if isa(input, 'gpuArray')
            error('Projector type 6 MATLAB implementation does not support gpuArray input.')
        end
        output = backwardProjectionType6(options, input, koko, norm_input, corr_input, nMeas, subIter);
        if noSensIm == 0
            sensIm = backwardProjectionType6(options, ones(size(input), 'like', input), koko, norm_input, corr_input, nMeas, subIter);
            if iscell(sensIm)
                for volume = 1:numel(sensIm)
                    sensIm{volume}(sensIm{volume} < options.epps) = 1;
                end
            else
                sensIm(sensIm < options.epps) = 1;
            end
        else
            sensIm = 0;
        end
        options.ub = options.ub + koko;
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
        for ii = 1 : options.nMultiVolumes + 1
            if outputCell
                if inputCell
                    if options.useSingles
                        [output{ii}, sensIm{ii}] = projector_mexSingle( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ...
                            options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...
                            options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ...
                            options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, input{ii}, noSensIm, 2, ...
                            dPitchXY, dPitchZ, options.nProjections);
                    else
                        [output{ii}, sensIm{ii}] = projector_mex( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ...
                            options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...
                            options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ...
                            options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, input{ii}, noSensIm, 2, ...
                            dPitchXY, dPitchZ, options.nProjections);
                    end
                else
                    if options.useSingles
                        [output{ii}, sensIm{ii}] = projector_mexSingle( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ...
                            options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...
                            options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ...
                            options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, input, noSensIm, 2, ...
                            dPitchXY, dPitchZ, options.nProjections);
                    else
                        [output{ii}, sensIm{ii}] = projector_mex( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ...
                            options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...
                            options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ...
                            options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, input, noSensIm, 2, ...
                            dPitchXY, dPitchZ, options.nProjections);
                    end
                end
                if options.use_psf
                    output{ii} = computeConvolution(output{ii}, options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.gaussK);
                    if numel(sensIm{ii}) > 1
                        sensIm{ii} = computeConvolution(sensIm{ii}, options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.gaussK);
                    end
                end
            else
                if options.useSingles
                    [output, sensIm] = projector_mexSingle( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ...
                        options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...
                        options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ...
                        options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, input, noSensIm, 2, ...
                        dPitchXY, dPitchZ, options.nProjections);
                else
                    [output, sensIm] = projector_mex( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ...
                        options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...
                        options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ...
                        options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, input, noSensIm, 2, ...
                        dPitchXY, dPitchZ, options.nProjections);
                end
                if options.use_psf
                    output = computeConvolution(output, options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.gaussK);
                    if numel(sensIm) > 1
                        sensIm = computeConvolution(sensIm, options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.gaussK);
                    end
                end
            end
        end
    end
elseif options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || ismac
    % Determined here, before input is reshaped/converted below, so that
    % the gpuArray-ness of the original input is not lost
    useGPUArray = isa(input, 'gpuArray');
    if useGPUArray
        if options.implementation ~= 5
            error('gpuArray input is only supported with implementation 5!')
        end
        if exist('CUDA_matrixfree_multi_gpu','file') ~= 3
            error(['gpuArray input was used, but the CUDA build of implementation 5 (CUDA_matrixfree_multi_gpu) was not found. ' ...
                'Run install_mex to build it. gpuArray input requires a CUDA-enabled installation and the CUDA toolkit.'])
        end
    end
    if ~isfield(options,'orthTransaxial') && (options.projector_type == 2 || options.projector_type == 3 || options.projector_type == 21 || options.projector_type == 22 || options.projector_type == 33)
        if options.projector_type == 3 || options.projector_type == 33
            options.orthTransaxial = true;
        elseif (options.projector_type == 2 || options.projector_type == 21 || options.projector_type == 22) && (isfield(options,'tube_width_xy') && options.tube_width_xy > 0 || options.SPECT)
            options.orthTransaxial = true;
        else
            options.orthTransaxial = false;
        end
    end
    if ~isfield(options,'orthAxial') && (options.projector_type == 2 || options.projector_type == 3 || options.projector_type == 21 || options.projector_type == 22 || options.projector_type == 33)
        if options.projector_type == 3 || options.projector_type == 33
            options.orthAxial = true;
        elseif (options.projector_type == 2 || options.projector_type == 21 || options.projector_type == 22) && (isfield(obj.param,'tube_width_z') && obj.param.tube_width_z > 0 || options.SPECT)
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
    options.x0 = zeros(options.Nx(1) * options.Ny(1) * options.Nz(1), 1, 'single');
    options.use_device = uint32(options.use_device);

    if options.orthAxial
        crystal_size_z = (options.tube_width_z);
    else
        crystal_size_z = (options.tube_width_xy);
    end
    if options.projector_type == 1 || options.projector_type == 11 ...
            || options.projector_type == 2 || options.projector_type == 3 || options.projector_type == 21 || options.projector_type == 22 || options.projector_type == 33
        kernel_file = 'projectorType123.cl';
        kernel_path = which(kernel_file);
        kernel_path = strrep(kernel_path, '\', '/');
        kernel_path = strrep(kernel_path, '.cl', '');
        header_directory = strrep(kernel_path,'projectorType123','');
    elseif options.projector_type == 4 || options.projector_type == 41 || options.projector_type == 14 || options.projector_type == 45
        kernel_file = 'projectorType4.cl';
        kernel_path = which(kernel_file);
        kernel_path = strrep(kernel_path, '\', '/');
        kernel_path = strrep(kernel_path, '.cl', '');
        header_directory = strrep(kernel_path,'projectorType4','');
    elseif options.projector_type == 5 || options.projector_type == 51 || options.projector_type == 15 || options.projector_type == 54 || options.projector_type == 45
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
    if listmodeSPECTSensitivity
        % The sensitivity kernel consumes one weight per full projection view,
        % with geometry stored independently in x/z. Flatten the per-frame
        % request to a one-timestep, one-subset MEX call.
        nSensitivityViews = size(options.sensitivityViewWeights, 1);
        if isfield(options, 'sensitivityTimestep') && ~isempty(options.sensitivityTimestep)
            sensitivityTimestep = double(options.sensitivityTimestep) + 1;
        elseif isfield(options, 'currentTimestep') && ~isempty(options.currentTimestep)
            sensitivityTimestep = double(options.currentTimestep) + 1;
        else
            sensitivityTimestep = 1;
        end
        if sensitivityTimestep < 1 || sensitivityTimestep > size(options.sensitivityViewWeights, 2) || ...
                sensitivityTimestep ~= floor(sensitivityTimestep)
            error('The requested SPECT sensitivity timestep is outside sensitivityViewWeights.')
        end
        input = options.sensitivityViewWeights(:, sensitivityTimestep);
        options.sensitivityViewWeights = input;
        options.Nt = 1;
        options.partitions = 1;
        options.subsets = 1;
        options.currentSubset = 0;
        options.currentTimestep = 0;
        options.nProjections = nSensitivityViews;
        options.totMeas = nSensitivityViews;
        nMeas = int64([0; nSensitivityViews]);
        subIter = 0;
        noSensIm = true;
    end
    if numel(nMeas) == 1
        nMeas = [0;nMeas];
    end
    if options.projector_type == 5 || options.projector_type == 15 || options.projector_type == 45
        input = reshape(input, options.nRowsD, options.nColsD, nMeas(subIter + 2) - nMeas(subIter + 1));
        [input, meanBP] = computeIntegralImage(input, options.meanBP);
        if options.meanBP
            % One mean per projection in this subset's measurement block
            options.meanBPvec = single(meanBP(:));
        end
    end
    if ~isa(input,'single')
        input = single(input);
    end
    if useGPUArray
        options = gatherHostFields(options, {});
        mexProjector = @CUDA_matrixfree_multi_gpu;
    else
        mexProjector = @OpenCL_matrixfree_multi_gpu;
    end
    [output, sensIm] = mexProjector( options.Nx, options.Ny, options.Nz, options.dx, options.dy, options.dz, options.bx, options.by, options.bz, ...
        z, x, options.nRowsD, options.verbose, options.LL, options.TOF, ... % 15
        TOFSize, options.sigma_x, options.TOFCenter, options.TOF_bins, options.platform, options.use_raw_data, options.use_psf, header_directory, options.vaimennus, ... % 24
        options.normalization, nMeas, options.attenuation_correction, options.normalization_correction, 1, options.subsets, options.epps, options.xy_index, ...
        options.z_index, crystal_size_z, ... % 34
        options.x_center, options.y_center, options.z_center, single(0), 0, options.projector_type, n_rays, n_rays3D, ... % 42
        options, input, options.partitions, options.use_64bit_atomics, options.bmin, options.bmax, options.Vmax, options.V, options.gaussK, 2, noSensIm); % 51
    if options.use_64bit_atomics
        output = single(output) / 99999997952;
        if ~isempty(sensIm) && numel(sensIm) > 1
            sensIm = single(sensIm) / 99999997952;
        end
    elseif options.use_32bit_atomics
        output = single(output) / 100000;
        if ~isempty(sensIm) && numel(sensIm) > 1
            sensIm = single(sensIm) / 100000;
        end
    end
    if outputCell
        temp = output;
        output = cell(options.nMultiVolumes + 1, 1);
        alku = 1;
        for kk = 1 : options.nMultiVolumes + 1
            output{kk} = temp(alku : alku - 1 + prod(options.N(kk)));
            % Was previously "alku = prod(options.N(kk)) + 1;", which discards the
            % offset accumulated over previous volumes instead of advancing past this
            % volume's block. That is a no-op bug for exactly 2 volumes (nMultiVolumes
            % == 1, since alku isn't read again after the loop's last iteration) but
            % reads the wrong slice of temp for 3 or more volumes.
            alku = alku + prod(options.N(kk));
            if options.use_psf
                output{kk} = computeConvolution(output{kk}, options, options.Nx(kk), options.Ny(kk), options.Nz(kk), options.gaussK);
            end
        end
        if numel(sensIm) > 1
            temp = sensIm;
            sensIm = cell(options.nMultiVolumes + 1, 1);
            alku = 1;
            for kk = 1 : options.nMultiVolumes + 1
                sensIm{kk} = temp(alku : alku - 1 + prod(options.N(kk)));
                % Same fix as the output{} split above (this variant was additionally
                % missing the "+ 1", so it was wrong even for exactly 2 volumes).
                alku = alku + prod(options.N(kk));
                if options.use_psf
                    sensIm{kk} = computeConvolution(sensIm{kk}, options, options.Nx(kk), options.Ny(kk), options.Nz(kk), options.gaussK);
                end
            end
        end
    else
        if options.use_psf
            for ii = 1 : options.nMultiVolumes + 1
                output = computeConvolution(output, options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.gaussK);
                if numel(sensIm) > 1
                    sensIm = computeConvolution(sensIm, options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.gaussK);
                end
            end
        end
    end
else
    error('Unsupported implementation! This is unintended behavior!')
end
end
