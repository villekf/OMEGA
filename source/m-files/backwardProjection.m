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
    % 1. Apply central-ray, normalization and additional correction weights.
    fProj = fProj .* type6CentralRayWeights(options, viewIndices);
    fProj = fProj .* type6MeasurementWeights(options, norm_input, corr_input, viewIndices);
    nVolumes = double(options.nMultiVolumes) + 1;
    if options.nMultiVolumes > 0
        output = cell(nVolumes, 1);
    else
        output = zeros(prod(double([options.Nx(1), options.Ny(1), options.Nz(1)])), 1, options.cType);
    end
    for volume = 1:nVolumes
        [geom, gridSize] = type6VolumeGeometry(options, volume);
        gridNx = gridSize(1); gridNy = gridSize(2);
        [gridNz, padLo, padHi] = type6AxialGeometry(options, volume, gridNx);
        transA = double(options.dy(1)) / double(options.dy(volume));
        transB = 0.5 * gridNy + (0.5 - 0.5 * options.nRowsD) * transA - 0.5;
        axialA = double(options.dz(1)) / double(options.dz(volume));
        axialB = -double(options.bz(volume)) / double(options.dz(volume)) + ...
            (0.5 - 0.5 * options.nColsD) * axialA - 0.5 + padLo;
        volumeOutput = zeros(prod(double([options.Nx(volume), options.Ny(volume), options.Nz(volume)])), 1, options.cType);
        u1 = options.ub;
        for kk = 1:koko
            panelShift = double(options.blurPlanes2F(u1)) * double(options.dx(1)) / double(options.dx(volume));

            % 2. Spread this detector view onto the transaxial and axial grids.
            detector = fProj(:, :, kk);
            detectorRows = size(detector, 1); rest = size(detector, 2);
            detectorValues = reshape(detector, detectorRows, rest);
            gridValues = zeros(gridNy, rest, 'like', detector);
            for detectorRow = 0:detectorRows-1
                position = transA * detectorRow + transB;
                lower = floor(position); fraction = position - lower;
                if lower >= 0 && lower < gridNy
                    gridValues(lower + 1, :) = gridValues(lower + 1, :) + ...
                        detectorValues(detectorRow + 1, :) .* cast(1 - fraction, 'like', detector);
                end
                upper = lower + 1;
                if upper >= 0 && upper < gridNy
                    gridValues(upper + 1, :) = gridValues(upper + 1, :) + ...
                        detectorValues(detectorRow + 1, :) .* cast(fraction, 'like', detector);
                end
            end
            detector = reshape(gridValues, gridNy, size(detector, 2));

            detectorRows = size(detector, 2); rest = size(detector, 1);
            detectorValues = reshape(detector.', detectorRows, rest);
            gridValues = zeros(gridNz, rest, 'like', detector);
            for detectorRow = 0:detectorRows-1
                position = axialA * detectorRow + axialB;
                lower = floor(position); fraction = position - lower;
                if lower >= 0 && lower < gridNz
                    gridValues(lower + 1, :) = gridValues(lower + 1, :) + ...
                        detectorValues(detectorRow + 1, :) .* cast(1 - fraction, 'like', detector);
                end
                upper = lower + 1;
                if upper >= 0 && upper < gridNz
                    gridValues(upper + 1, :) = gridValues(upper + 1, :) + ...
                        detectorValues(detectorRow + 1, :) .* cast(fraction, 'like', detector);
                end
            end
            detector = reshape(gridValues, gridNz, size(detector, 1)).';

            % 3. Extend across depth and apply voxel-to-ray scaling.
            kuvaRot = repmat(reshape(detector, 1, gridNy, gridNz), gridNx, 1, 1);
            kuvaRot = kuvaRot * (options.dx(volume) / options.dx(1)) / (options.dy(1) * options.dz(1));

            % 4. Apply the adjoint detector response on each depth plane.
            if isfield(options, 'CDRFOnTheFly') && options.CDRFOnTheFly
                PSF = type6CDRFKernel(options, volume, u1, gridSize, gridNz);
            elseif ~isempty(options.gFilter)
                PSF = imtranslate(options.gFilter, [0, 0, options.blurPlanes(u1)], FillValues=0);
                PSF = PSF(:, :, 1:gridNx);
                PSF = permute(PSF, [3 2 1]);
            else
                PSF = [];
            end
            if ~isempty(PSF)
                for depth = 1:gridNx
                    image = squeeze(kuvaRot(depth, :, :)); kernel = squeeze(PSF(depth,:,:));
                    [nRows, nCols] = size(image); [kRows, kCols] = size(kernel);
                    rowOrigin = floor(kRows / 2) + 1; colOrigin = floor(kCols / 2) + 1;
                    convolved = zeros(size(image), 'like', image);
                    for kernelRow = 1:kRows
                        rowOffset = kernelRow - rowOrigin;
                        if rowOffset >= 0
                            sourceRows = (1 + rowOffset):nRows;
                            targetRows = 1:(nRows - rowOffset);
                        else
                            sourceRows = 1:(nRows + rowOffset);
                            targetRows = (1 - rowOffset):nRows;
                        end
                        if isempty(sourceRows), continue, end
                        for kernelCol = 1:kCols
                            colOffset = kernelCol - colOrigin;
                            if colOffset >= 0
                                sourceCols = (1 + colOffset):nCols;
                                targetCols = 1:(nCols - colOffset);
                            else
                                sourceCols = 1:(nCols + colOffset);
                                targetCols = (1 - colOffset):nCols;
                            end
                            if ~isempty(sourceCols)
                                convolved(targetRows, targetCols) = convolved(targetRows, targetCols) + ...
                                    image(sourceRows, sourceCols) .* cast(kernel(kernelRow, kernelCol), 'like', image);
                            end
                        end
                    end
                    kuvaRot(depth, :, :) = convolved;
                end
            end

            % 5. Apply attenuation factors along each rotated ray.
            if options.attenuation_correction
                attenuationImage = type6AttenuationFactors(options, volume, u1, gridSize, gridNz, padLo, panelShift);
                kuvaRot = kuvaRot .* attenuationImage;
            end

            % 6. Rotate back with the exact transpose of the forward gather.
            theta = double(options.swivelAngles(u1));
            centerX = double(options.bx(volume)) + 0.5 * double(options.Nx(volume)) * double(options.dx(volume));
            centerY = double(options.by(volume)) + 0.5 * double(options.Ny(volume)) * double(options.dy(volume));
            mainCenter = double([options.bx(1), options.by(1), options.bz(1)]) + ...
                0.5 * double([options.Nx(1), options.Ny(1), options.Nz(1)]) .* ...
                double([options.dx(1), options.dy(1), options.dz(1)]);
            mainSpacing = double([options.dx(1), options.dy(1), options.dz(1)]);
            mainCentered = all(abs(mainCenter) < 1e-3 * mainSpacing);
            needsWindow = volume > 1 || ~mainCentered || padLo ~= padHi || ...
                double(options.nRowsD) > gridNy || double(options.nColsD) > gridNz;
            if ~needsWindow && mod(gridNy - double(options.nRowsD), 2) ~= 0
                panelShift = panelShift + 0.5;
            end
            offX = -sind(theta) * panelShift - centerX / double(options.dx(volume));
            offY = cosd(theta) * panelShift - centerY / double(options.dy(volume));
            kuvaRot = rotateType6Adjoint(kuvaRot, -theta, offX, offY, geom);
            kuvaRot = kuvaRot(:, :, padLo + (1:double(options.Nz(volume))));
            volumeOutput = volumeOutput + kuvaRot(:);
            u1 = u1 + 1;
        end
        % 7. Keep only voxels selected by the backprojection mask.
        volumeOutput = volumeOutput .* cast(backprojectionMaskType6(options, volume), 'like', volumeOutput);
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

function output = rotateType6Forward(input, angleDeg, offX, offY, geom)
% Gather the source volume on the rotated type-6 grid with zero outside it.
srcNx = size(input, 1); srcNy = size(input, 2); nPlanes = size(input, 3);
outNx = double(geom.outputSize(1)); outNy = double(geom.outputSize(2));
[rowOut, colOut] = ndgrid(1:outNx, 1:outNy);
qRow = rowOut - 0.5 - outNx / 2; qCol = colOut - 0.5 - outNy / 2;
c = cosd(-double(angleDeg)); s = sind(-double(angleDeg));
rowSrc = c .* qRow - s .* qCol + double(offX) + srcNx / 2 + 0.5;
colSrc = s .* qRow + c .* qCol + double(offY) + srcNy / 2 + 0.5;
inside = rowSrc > 0 & rowSrc < srcNx + 1 & colSrc > 0 & colSrc < srcNy + 1;
row0 = floor(rowSrc); col0 = floor(colSrc);
fr = rowSrc - row0; fc = colSrc - col0;
rows = {row0, row0 + 1, row0, row0 + 1};
cols = {col0, col0, col0 + 1, col0 + 1};
wr = {(1 - fr) .* (1 - fc), fr .* (1 - fc), (1 - fr) .* fc, fr .* fc};
output = zeros(outNx, outNy, nPlanes, 'like', input);
for plane = 1:nPlanes
    values = zeros(outNx, outNy, 'like', input);
    for neighbor = 1:4
        valid = inside & rows{neighbor} >= 1 & rows{neighbor} <= srcNx & ...
            cols{neighbor} >= 1 & cols{neighbor} <= srcNy;
        if any(valid(:))
            destination = find(valid(:));
            source = sub2ind([srcNx, srcNy], rows{neighbor}(valid), cols{neighbor}(valid));
            flat = values(:);
            flat(destination) = flat(destination) + input(source) .* cast(wr{neighbor}(valid), 'like', input);
            values = reshape(flat, outNx, outNy);
        end
    end
    output(:, :, plane) = values;
end
end

function output = rotateType6Adjoint(input, angleDeg, offX, offY, geom)
% Scatter the rotated image back to the source grid using the transpose weights.
srcNx = double(geom.inputSize(1)); srcNy = double(geom.inputSize(2));
outNx = double(geom.outputSize(1)); outNy = double(geom.outputSize(2));
nPlanes = size(input, 3);
[rowOut, colOut] = ndgrid(1:outNx, 1:outNy);
qRow = rowOut - 0.5 - outNx / 2; qCol = colOut - 0.5 - outNy / 2;
c = cosd(-double(angleDeg)); s = sind(-double(angleDeg));
rowSrc = c .* qRow - s .* qCol + double(offX) + srcNx / 2 + 0.5;
colSrc = s .* qRow + c .* qCol + double(offY) + srcNy / 2 + 0.5;
inside = rowSrc > 0 & rowSrc < srcNx + 1 & colSrc > 0 & colSrc < srcNy + 1;
row0 = floor(rowSrc); col0 = floor(colSrc);
fr = rowSrc - row0; fc = colSrc - col0;
rows = {row0, row0 + 1, row0, row0 + 1};
cols = {col0, col0, col0 + 1, col0 + 1};
wr = {(1 - fr) .* (1 - fc), fr .* (1 - fc), (1 - fr) .* fc, fr .* fc};
output = zeros(srcNx, srcNy, nPlanes, 'like', input);
for plane = 1:nPlanes
    values = zeros(srcNx, srcNy, 'like', input);
    for neighbor = 1:4
        valid = inside & rows{neighbor} >= 1 & rows{neighbor} <= srcNx & ...
            cols{neighbor} >= 1 & cols{neighbor} <= srcNy;
        if any(valid(:))
            destination = sub2ind([srcNx, srcNy], rows{neighbor}(valid), cols{neighbor}(valid));
            weight = cast(wr{neighbor}(valid), 'like', input);
            contribution = input(valid) .* weight;
            summed = accumarray(double(destination), double(contribution(:)), [srcNx * srcNy, 1], @sum, 0);
            values = values + reshape(cast(summed, 'like', input), srcNx, srcNy);
        end
    end
    output(:, :, plane) = values;
end
end

function mask = backprojectionMaskType6(options, volume)
% Select the type-6 BP mask for one volume and frame.
nx = double(options.Nx(volume)); ny = double(options.Ny(volume)); nz = double(options.Nz(volume));
mask = true(nx * ny * nz, 1);
if ~isfield(options, 'useMaskBP') || ~options.useMaskBP || ...
        ~isfield(options, 'maskBP') || isempty(options.maskBP) || isscalar(options.maskBP)
    return
end
source = options.maskBP;
if iscell(source)
    frame = 1;
    if isfield(options, 'currentTimestep'), frame = double(options.currentTimestep) + 1; end
    source = source{min(frame, numel(source))};
end
source = source(:);
counts = double(options.Nx(:)) .* double(options.Ny(:)) .* double(options.Nz(:));
nFrames = 1;
if isfield(options, 'Nt'), nFrames = double(options.Nt); end
if options.nMultiVolumes > 0 && any(numel(source) == sum(counts) * [1 nFrames]) && ...
        options.maskBPZ > 1 && options.maskBPZ == max(options.Nz)
    frame = 0;
    if numel(source) > sum(counts) && isfield(options, 'currentTimestep')
        frame = min(double(options.currentTimestep), nFrames - 1);
    end
    first = frame * sum(counts) + sum(counts(1:volume-1)) + 1;
    mask = source(first:first + counts(volume) - 1) ~= 0;
elseif options.nMultiVolumes == 0 && numel(source) == nx * ny
    mask = repmat(source ~= 0, nz, 1);
elseif options.nMultiVolumes == 0 && numel(source) == nx * ny * nz
    mask = source ~= 0;
elseif options.nMultiVolumes > 0
    origins = single([options.bx(:), options.by(:), options.bz(:)]);
    steps = single([options.dx(:), options.dy(:), options.dz(:)]);
    dims = single([options.Nx(:), options.Ny(:), options.Nz(:)]);
    low = min(origins, [], 1);
    grid = double(round((max(origins + dims .* steps, [], 1) - low) ./ steps(1, :)));
    is3D = options.maskBPZ > 1;
    if numel(source) ~= prod(grid(1:2)) * (1 + (grid(3) - 1) * is3D)
        error('backprojectionMaskType6:InvalidSize', 'The BP mask does not match the volume union.');
    end
    if ~is3D, grid(3) = 1; end
    coordinates = cell(1, 3);
    for axis = 1:3
        coordinates{axis} = floor((double(origins(volume, axis)) + ...
            ((0:double(dims(volume, axis))-1) + .5) * double(steps(volume, axis)) - double(low(axis))) ...
            / double(steps(1, axis)) + 1e-4) + 1;
    end
    if ~is3D, coordinates{3}(:) = 1; end
    [x, y, z] = ndgrid(coordinates{:});
    valid = x >= 1 & x <= grid(1) & y >= 1 & y <= grid(2) & z >= 1 & z <= grid(3);
    mask = false(nx, ny, nz);
    mask(valid) = source(sub2ind(grid, x(valid), y(valid), z(valid))) ~= 0;
    mask = mask(:);
else
    error('backprojectionMaskType6:InvalidSize', 'The BP mask does not match the image.');
end
end

function [geom, gridSize] = type6VolumeGeometry(options, volume)
% Physical XY geometry of the smallest axis-centered square grid that contains
% one reconstruction volume, following rotGridType6 in precomp.h.
volume = double(volume);
Nx = double(options.Nx(volume)); Ny = double(options.Ny(volume));
if volume == 1 && numel(options.Nx) == 1 && ...
        abs(double(options.bx(1)) + 0.5*Nx*double(options.dx(1))) < 1e-3*double(options.dx(1)) && ...
        abs(double(options.by(1)) + 0.5*Ny*double(options.dy(1))) < 1e-3*double(options.dy(1))
    geom=struct('inputSize',[Nx,Ny],'outputSize',[Nx,Ny], ...
        'inputSpacing',[double(options.dx(1)),double(options.dy(1))], ...
        'outputSpacing',[double(options.dx(1)),double(options.dy(1))], ...
        'inputOffset',[0,0],'outputOffset',[0,0]);
    gridSize=[Nx,Ny]; return
end
dx = double(options.dx(volume)); dy = double(options.dy(volume));
bx = double(options.bx(volume)); by = double(options.by(volume));
x0 = bx; x1 = bx + Nx * dx; y0 = by; y1 = by + Ny * dy;
R = sqrt(max(x0.^2, x1.^2) + max(y0.^2, y1.^2));
Rcol = R;
if 0.5 * double(options.nRowsD) * double(options.dy(1)) > R
    shifts = abs(double(options.radiusPerProj(:)) .* sind(double(options.swivelAngles(:)) - double(options.angles(:)) + 180));
    if ~isempty(shifts), Rcol = R + max(shifts); end
end
nxGrid = fitType6GridDimension(R, dx, double(options.Nx(1)), volume == 1);
nyGrid = fitType6GridDimension(Rcol, dy, double(options.Ny(1)), volume == 1);
centerX = bx + 0.5 * Nx * dx;
centerY = by + 0.5 * Ny * dy;
geom = struct('inputSize', [Nx, Ny], 'outputSize', [nxGrid, nyGrid], ...
    'inputSpacing', [dx, dy], 'outputSpacing', [dx, dy], ...
    'inputOffset', [centerY, centerX], 'outputOffset', [0, 0]);
gridSize = [nxGrid, nyGrid];
end

function n = fitType6GridDimension(radius, spacing, mainSize, isMain)
n = ceil(2 * radius / spacing);
if isMain, parity = mod(round(mainSize), 2); else, parity = 0; end
if mod(n, 2) ~= parity, n = n + 1; end
n = max(n, 1);
end

function [gridNz, padLo, padHi] = type6AxialGeometry(options, volume, gridNx)
% Axial blur spill padding from computePadZType6 in precomp.h.
volume = double(volume);
Nz = double(options.Nz(volume)); dz = double(options.dz(volume));
padLo = 0; padHi = 0; padZ = 0;
needsPad = numel(options.Nx) > 1 || ...
    abs(double(options.bx(1)) + 0.5 * double(options.Nx(1)) * double(options.dx(1))) > 1e-3 * double(options.dx(1)) || ...
    abs(double(options.by(1)) + 0.5 * double(options.Ny(1)) * double(options.dy(1))) > 1e-3 * double(options.dy(1)) || ...
    abs(double(options.bz(1)) + 0.5 * double(options.Nz(1)) * double(options.dz(1))) > 1e-3 * double(options.dz(1));
if needsPad && isfield(options, 'CDRFCoeffs6') && numel(options.CDRFCoeffs6) >= 5 && ...
        isfield(options, 'blurPlanesF') && ~isempty(options.blurPlanesF)
    coeff = double(options.CDRFCoeffs6(:));
    minDist = min(double(options.blurPlanesF(:)));
    gridNx = double(gridNx);
    collimator = (minDist - 0.5 * double(options.Nx(1))) * double(options.dx(1));
    dist = (collimator + 0.5 * gridNx * double(options.dx(volume))) / double(options.dx(volume));
    zMax = (gridNx - 0.5 - dist) * double(options.dx(volume));
    sigmaCoeff = 1 / (2 * sqrt(2 * log(2)));
    extentMax = 0;
    model = 0;
    if isfield(options, 'cdrfModel6'), model = double(options.cdrfModel6); end
    for z = [0, max(zMax, 0)]
        Rg = max(coeff(1) * z + coeff(2), 0);
        C = coeff(5);
        if model == 1
            extent = Rg / dz + 5.2565 * (C / dz * sigmaCoeff);
        else
            extent = 5.2565 * (sqrt(Rg^2 + C^2) / dz * sigmaCoeff);
        end
        extentMax = max(extentMax, extent);
    end
    radius = ceil(extentMax);
    padZ = ceil(radius / 8) * 8;
    zColLo = -0.5 * double(options.nColsD) * double(options.dz(1)) + 0.5 * double(options.dz(1));
    zColHi = 0.5 * double(options.nColsD) * double(options.dz(1)) - 0.5 * double(options.dz(1));
    posLo = (zColLo - double(options.bz(volume))) / dz - 0.5;
    posHi = (zColHi - double(options.bz(volume))) / dz - 0.5;
    needLo = -floor(posLo + 1e-4);
    needHi = ceil(posHi - 1e-4) - (Nz - 1);
    padLo = max(min(padZ, needLo), 0);
    padHi = max(min(padZ, needHi), 0);
end
gridNz = Nz + padLo + padHi;
end

function kernel = type6CDRFKernel(options, volume, view, gridXY, gridNz)
% Per-view Gaussian CDRF on the expanded rotation grid (Python reference).
volume=double(volume); view=double(view);
if isfield(options,'cdrfModel6') && double(options.cdrfModel6) ~= 0
    error('type6CDRFKernel:UnsupportedCDRF','The MATLAB type-6 path currently supports only the Gaussian on-the-fly CDRF model.');
end
if isfield(options,'CDRFOnTheFly') && ~options.CDRFOnTheFly
    error('type6CDRFKernel:UnsupportedCDRF','Manual type-6 PSF banks are not supported by the MATLAB type-6 path.');
end
gridXY=double(gridXY(:).'); gridNz=double(gridNz);
gridNx=gridXY(1); gridNy=gridXY(2); gridNz=double(gridNz);
dx=double(options.dx(volume)); dy=double(options.dy(volume)); dz=double(options.dz(volume)); dx0=double(options.dx(1));
distFloat=double(options.blurPlanesF(view)); collimator=(distFloat-.5*double(options.Nx(1)))*dx0;
distance=(collimator+.5*gridNx*dx)/dx;
firstPlane=max(0,min(gridNx-1,ceil(distance-.5))); firstPlane=floor(firstPlane/16)*16;
coeff=double(options.CDRFCoeffs6(:)); Az=coeff(1); Bz=coeff(2); Axy=coeff(3); Bxy=coeff(4); intrinsic=coeff(5);
sigmaCoeff=1/(2*sqrt(2*log(2)));
planes=(0:gridNx-1).'; depth=(planes+.5-distance)*dx;
rgZ=max(Az.*depth+Bz,0); rgXY=max(Axy.*depth+Bxy,0);
sigmaZ=sqrt(rgZ.^2+intrinsic^2)./dz.*sigmaCoeff;
fwhmXY=sqrt(rgXY.^2+intrinsic^2)./dy;
if ~isfield(options,'type6RotationComp') || options.type6RotationComp
    comp=type6RotationVariance(options,volume,view,gridXY);
else
    comp=sigmaCoeff^2;
end
sigmaXY=sqrt(max((fwhmXY.*sigmaCoeff).^2-comp,1e-16*sigmaCoeff^2));
extent=0;
for zmm=[(firstPlane+.5-distance)*dx,(gridNx-.5-distance)*dx]
    rz=max(Az*zmm+Bz,0); rxy=max(Axy*zmm+Bxy,0);
    extent=max(extent,5.2565*max(sqrt(rz^2+intrinsic^2)/dz,sqrt(rxy^2+intrinsic^2)/dy)*sigmaCoeff);
end
radius=min(ceil(extent/8)*8,max(gridNy,gridNz));
taps=-radius:radius;
kernel=zeros(gridNx,2*radius+1,2*radius+1,options.cType);
for ix=1:gridNx
    if ix-1<firstPlane || depth(ix)<=0, continue, end
    sy=max(sigmaXY(ix),1e-8); sz=max(sigmaZ(ix),1e-8);
    trans=exp(-(taps.^2)/(2*sy^2)); axial=exp(-(taps.^2)/(2*sz^2));
    trans(abs(taps)>ceil(5.2565*sy))=0;
    axial(abs(taps)>ceil(5.2565*sz))=0;
    if sum(trans)>0, trans=trans/sum(trans); end
    if sum(axial)>0, axial=axial/sum(axial); end
    kernel(ix,:,:)=cast(trans(:)*axial(:).','like',kernel);
end
end

function variance=type6RotationVariance(options,volume,view,grid)
% Port of Python _type6_rotation_variance; exact tent-weighted variance.
nx=double(options.Nx(volume)); ny=double(options.Ny(volume)); outNx=grid(1); outNy=grid(2);
theta=deg2rad(double(options.swivelAngles(view))); c=cos(theta); s=sin(theta);
panel=double(options.blurPlanes2F(view))*double(options.dx(1))/double(options.dx(volume));
axisX=-(double(options.bx(volume))+.5*nx*double(options.dx(volume)))/double(options.dx(volume));
axisY=-(double(options.by(volume))+.5*ny*double(options.dy(volume)))/double(options.dy(volume));
offX=-s*panel+axisX; offY=c*panel+axisY; hh=abs(c)+abs(s);
csx=.5*nx-.5; csy=.5*ny-.5; cox=.5*outNx-.5; coy=.5*outNy-.5;
if abs(s)<1e-6 || abs(c)<1e-6
    samples=[floor(nx/2),floor(ny/2)];
elseif nx*ny<=1024
    [ii,jj]=ndgrid(0:nx-1,0:ny-1); samples=[ii(:),jj(:)];
else
    % Deterministic Halton sampling, matching the Python fallback.
    samples=zeros(1024,2);
    for k=1:1024, samples(k,:)=[min(floor(type6Halton(k,2)*nx),nx-1),min(floor(type6Halton(k,3)*ny),ny-1)]; end
end
num=0; den=0;
for k=1:size(samples,1)
    ux=samples(k,1)-csx; uy=samples(k,2)-csy;
    qx=c*(ux-offX)+s*(uy-offY); qy=-s*(ux-offX)+c*(uy-offY);
    jx0=ceil(qx+cox-hh); jx1=floor(qx+cox+hh); jy0=ceil(qy+coy-hh); jy1=floor(qy+coy+hh);
    for jy=jy0:jy1
        ay=jy-coy;
        for jx=jx0:jx1
            ax=jx-cox; sx=c*ax-s*ay+offX; sy=s*ax+c*ay+offY;
            weight=max(1-abs(sx-ux),0)*max(1-abs(sy-uy),0);
            num=num+weight*(ay-qy)^2; den=den+weight;
        end
    end
end
if den>0, variance=num/den; else, variance=1/6; end
if volume>1
    nDet=double(options.nRowsD); a=double(options.dy(1))/double(options.dy(volume)); b=.5*outNy+(.5-.5*nDet)*a-.5;
    num=0; den=0;
    for row=0:nDet-1
        t=a*row+b; low=floor(t);
        for off=0:1, d=t-(low+off); w=max(1-abs(d),0); num=num+w*d^2; den=den+w; end
    end
    if den>0, variance=variance+num/den; end
end
end

function value=type6Halton(index,base)
value=0; factor=1/base;
while index>0, value=value+factor*mod(index,base); index=floor(index/base); factor=factor/base; end
end

function factors = type6AttenuationFactors(options, volume, view, gridXY, gridNz, padLo, panelShift)
% Native type-6 per-view attenuation factors from the union image-domain map.
volume=double(volume); view=double(view);
gridXY=double(gridXY(:).'); gridNz=double(gridNz); padLo=double(padLo);
if ~isfield(options,'vaimennus') || isempty(options.vaimennus)
    factors=ones([gridXY,gridNz],options.cType); return
end
nVol=double(options.nMultiVolumes)+1;
Nx=double(options.Nx(1:nVol)); Ny=double(options.Ny(1:nVol)); Nz=double(options.Nz(1:nVol));
dx=double(options.dx(1:nVol)); dy=double(options.dy(1:nVol)); dz=double(options.dz(1:nVol));
bx=double(options.bx(1:nVol)); by=double(options.by(1:nVol)); bz=double(options.bz(1:nVol));
lo=[min(bx),min(by),min(bz)]; hi=[max(bx+Nx.*dx),max(by+Ny.*dy),max(bz+Nz.*dz)];
map=options.vaimennus;
if iscell(map)
    if numel(map)~=1, error('type6AttenuationFactors:DynamicAttenuationUnsupported','Per-frame type-6 image attenuation is not supported in MATLAB.'); end
    map=map{1};
end
if isfield(options,'multiResUnionDims') && ~isempty(options.multiResUnionDims)
    mapDims=double(options.multiResUnionDims(:).');
else
    mapDims=double([size(map,1),size(map,2),size(map,3)]);
end
if numel(mapDims)~=3 || any(mapDims<1) || numel(map)~=prod(mapDims)
    error('type6AttenuationFactors:AttenuationDimensions','The type-6 attenuation map must cover its union grid.');
end
map=reshape(single(map),mapDims(1),mapDims(2),mapDims(3));
mapSpacing=[dx(1),dy(1),dz(1)]; mapOrigin=lo;
% In the native public map the first voxel is the low corner of the volume union.
if any(abs((hi-lo)./mapSpacing-mapDims)>1.1) && isfield(options,'multiResPadLow')
    mapOrigin=[bx(1),by(1),bz(1)]-double(options.multiResPadLow(:).').*mapSpacing;
end
runion=0;
for jj=1:nVol
    runion=max(runion,sqrt(max(bx(jj)^2,(bx(jj)+Nx(jj)*dx(jj))^2)+max(by(jj)^2,(by(jj)+Ny(jj)*dy(jj))^2)));
end
gridNx=gridXY(1); gridNy=gridXY(2);
minNx=ceil(2*runion/dx(volume)); minNy=ceil(2*runion/dy(volume));
attNx=max(1,minNx+mod(minNx-gridNx,2)); srcNy=max(1,minNy+mod(minNy-gridNy,2));
srcNz=max(1,ceil((hi(3)-lo(3))/dz(volume)-1e-10));
sourceDims=[attNx,srcNy,srcNz]; sourceSpacing=[dx(volume),dy(volume),dz(volume)];
sourceOrigin=[-.5*attNx*sourceSpacing(1),-.5*srcNy*sourceSpacing(2),lo(3)];
source=resamplePhysical3D(map,mapOrigin,mapSpacing,sourceDims,sourceOrigin,sourceSpacing);
% Rotate onto the emission grid's axial-centred row lattice. This source offset
% is the native offsetsFPType6 panel shift expressed in physical [column,row].
swivel=double(options.swivelAngles(view));
panel=double(panelShift);
offX=-sind(swivel)*panel; offY=cosd(swivel)*panel;
geom=struct('inputSize',[attNx,srcNy],'outputSize',[attNx,gridNy]);
rotated=rotateType6Forward(source,-swivel,offX,offY,geom);
cumulative=cumsum(rotated,1); previous=zeros(size(cumulative),'like',cumulative); previous(2:end,:,:)=cumulative(1:end-1,:,:);
mode=2; envMode=getenv('TYPE6_ATTCONV'); if ~isempty(envMode), parsed=str2double(envMode(1)); if isfinite(parsed)&&ismember(parsed,[0,1,2]), mode=parsed; end, end
if mode==0
    factorsGrid=exp(-dx(volume).*cumulative);
elseif mode==1
    factorsGrid=exp(-dx(volume).*(cumulative+previous).*0.5);
else
    interval=dx(volume).*(cumulative-previous); average=ones(size(interval),'like',interval);
    nonzero=abs(interval)>1e-6; average(nonzero)=-expm1(-interval(nonzero))./interval(nonzero);
    small=~nonzero; average(small)=1-interval(small).*0.5+interval(small).^2/6;
    factorsGrid=exp(-dx(volume).*previous).*average;
end
rowStart=floor((attNx-gridNx)/2)+1; rows=rowStart:(rowStart+gridNx-1);
factors=ones(gridNx,gridNy,gridNz,'single');
for iz=1:gridNz
    zPhysical=bz(volume)+(iz-padLo-0.5)*dz(volume);
    zMap=(zPhysical-sourceOrigin(3))/sourceSpacing(3)+0.5;
    plane=sampleVolumePlane(factorsGrid,zMap);
    if rows(1)>=1 && rows(end)<=attNx, factors(:,:,iz)=plane(rows,:); end
end
factors=cast(factors,options.cType);
end

function output = resamplePhysical3D(input,inputOrigin,inputSpacing,outputDims,outputOrigin,outputSpacing)
% Mean-preserving box resampling via a trilinearly sampled summed-area table.
[nx,ny,nz]=size(input);
ex= min(max((outputOrigin(1)+(0:outputDims(1))*outputSpacing(1)-inputOrigin(1))/inputSpacing(1),0),nx);
ey= min(max((outputOrigin(2)+(0:outputDims(2))*outputSpacing(2)-inputOrigin(2))/inputSpacing(2),0),ny);
ez= min(max((outputOrigin(3)+(0:outputDims(3))*outputSpacing(3)-inputOrigin(3))/inputSpacing(3),0),nz);
sat=zeros(nx+1,ny+1,nz+1,'double');
sat(2:end,2:end,2:end)=cumsum(cumsum(cumsum(double(input),1),2),3);
[X,Y,Z]=ndgrid(ex+1,ey+1,ez+1); corners=trilinearSample(sat,X,Y,Z);
boxSum=corners(2:end,2:end,2:end)-corners(1:end-1,2:end,2:end)-corners(2:end,1:end-1,2:end)-corners(2:end,2:end,1:end-1) ...
    +corners(1:end-1,1:end-1,2:end)+corners(1:end-1,2:end,1:end-1)+corners(2:end,1:end-1,1:end-1)-corners(1:end-1,1:end-1,1:end-1);
nominalVolume=prod(outputSpacing./inputSpacing);
output=cast(boxSum./max(nominalVolume,1e-6),'like',input);
end

function output = sampleVolumePlane(input, zIndex)
lo=floor(zIndex); f=zIndex-lo; output=zeros(size(input,1),size(input,2),'like',input);
if lo>=1 && lo<=size(input,3), output=output+input(:,:,lo).*(1-f); end
hi=lo+1; if hi>=1 && hi<=size(input,3), output=output+input(:,:,hi).*f; end
end

function output = trilinearSample(input, x, y, z)
[nx,ny,nz]=size(input); x0=floor(x); y0=floor(y); z0=floor(z); fx=x-x0; fy=y-y0; fz=z-z0;
output=zeros(size(x),'like',input);
for oz=0:1, for oy=0:1, for ox=0:1
    ix=x0+ox; iy=y0+oy; iz=z0+oz;
    valid=ix>=1 & ix<=nx & iy>=1 & iy<=ny & iz>=1 & iz<=nz;
    if any(valid(:))
        idx=sub2ind([nx,ny,nz],ix(valid),iy(valid),iz(valid));
        wx=fx; if ox==0, wx=1-fx; end
        wy=fy; if oy==0, wy=1-fy; end
        wz=fz; if oz==0, wz=1-fz; end
        output(valid)=output(valid)+input(idx).*cast(wx(valid).*wy(valid).*wz(valid),'like',input);
    end
end, end, end
end
