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

function outputFP = forwardProjectionType6(options, recApu, loopVar, koko, nMeas, subIter)
    % Some implementations keep the full measurement table and identify the
    % requested subset through nMeas/subIter instead of calling splitInput.
    if numel(nMeas) > 1 && subIter >= 0 && subIter + 2 <= numel(nMeas)
        firstView = double(nMeas(subIter + 1));
        lastView = double(nMeas(subIter + 2));
        subsetViews = lastView - firstView;
        if ~isscalar(koko) || koko ~= subsetViews
            koko = subsetViews;
        end
        if numel(norm_input) == options.nRowsD * options.nColsD * double(nMeas(end))
            norm_input = norm_input(firstView * options.nRowsD * options.nColsD + 1:lastView * options.nRowsD * options.nColsD);
        end
        if numel(corr_input) == options.nRowsD * options.nColsD * double(nMeas(end))
            corr_input = corr_input(firstView * options.nRowsD * options.nColsD + 1:lastView * options.nRowsD * options.nColsD);
        end
    end
    outputFP = zeros(options.nRowsD, options.nColsD, koko, options.cType);
    if isempty(recApu)
        return
    end
    viewIndices = options.uu:(options.uu + koko - 1);
    for ii = loopVar
        u1 = options.uu;
        if iscell(recApu)
            if ii > numel(recApu) || isempty(recApu{ii}), continue, end
            apuArr = reshape(recApu{ii}, options.Nx(ii), options.Ny(ii), options.Nz(ii));
        else
            apuArr = reshape(recApu, options.Nx(ii), options.Ny(ii), options.Nz(ii));
        end
        [geom, gridSize] = type6VolumeGeometry(options, ii);
        gridNx = gridSize(1);
        gridNy = gridSize(2);
        [gridNz, padLo, padHi] = type6AxialGeometry(options, ii, gridNx);
        transA = double(options.dy(1)) / double(options.dy(ii));
        transB = 0.5 * gridNy + (0.5 - 0.5 * options.nRowsD) * transA - 0.5;
        axialA = double(options.dz(1)) / double(options.dz(ii));
        axialB = -double(options.bz(ii)) / double(options.dz(ii)) + ...
            (0.5 - 0.5 * options.nColsD) * axialA - 0.5 + padLo;

        for kk = 1:koko
            % 1. Rotate activity into the type-6 XY grid.
            theta = double(options.swivelAngles(u1));
            panel = double(options.blurPlanes2F(u1)) * double(options.dx(1)) / double(options.dx(ii));
            centerX = double(options.bx(ii)) + 0.5 * double(options.Nx(ii)) * double(options.dx(ii));
            centerY = double(options.by(ii)) + 0.5 * double(options.Ny(ii)) * double(options.dy(ii));
            mainCenter = double([options.bx(1), options.by(1), options.bz(1)]) + ...
                0.5 * double([options.Nx(1), options.Ny(1), options.Nz(1)]) .* ...
                double([options.dx(1), options.dy(1), options.dz(1)]);
            mainSpacing = double([options.dx(1), options.dy(1), options.dz(1)]);
            mainCentered = all(abs(mainCenter) < 1e-3 * mainSpacing);
            needsWindow = ii > 1 || ~mainCentered || padLo ~= padHi || ...
                double(options.nRowsD) > gridNy || double(options.nColsD) > gridNz;
            if ~needsWindow && mod(gridNy - double(options.nRowsD), 2) ~= 0
                panel = panel + 0.5;
            end
            offX = -sind(theta) * panel - centerX / double(options.dx(ii));
            offY = cosd(theta) * panel - centerY / double(options.dy(ii));
            kuvaRot = rotateType6Forward(apuArr, -theta, offX, offY, geom);

            % 2. Add axial padding for the detector response.
            kuvaRot = cat(3, zeros(gridNx, gridNy, padLo, 'like', kuvaRot), kuvaRot, ...
                zeros(gridNx, gridNy, padHi, 'like', kuvaRot));

            % 3. Apply attenuation factors along each rotated ray.
            if options.attenuation_correction
                attenuationImage = type6AttenuationFactors(options, ii, u1, gridSize, gridNz, padLo, ...
                    double(options.blurPlanes2F(u1)) * double(options.dx(1)) / double(options.dx(ii)));
                kuvaRot = kuvaRot .* attenuationImage;
            end

            % 4. Apply the detector response on each depth plane.
            if isfield(options, 'CDRFOnTheFly') && options.CDRFOnTheFly
                PSF = type6CDRFKernel(options, ii, u1, gridSize, gridNz);
            elseif ~isempty(options.gFilter)
                PSF = imtranslate(options.gFilter, [0, 0, options.blurPlanes(u1)], FillValues=0);
                PSF = PSF(:, :, 1:gridNx);
                PSF = permute(PSF, [3 2 1]);
            else
                PSF = [];
            end
            if ~isempty(PSF)
                for jj = 1 : gridNx
                    kuvaRot(jj, :, :) = conv2(squeeze(kuvaRot(jj, :, :)), squeeze(PSF(jj,:,:)), 'same');
                end
            end
            
            % 5. Sum depth and resample onto detector rows and columns.
            kuvaRot = reshape(sum(kuvaRot, 1), gridNy, gridNz);
            gridRows = size(kuvaRot, 1);
            detectorRows = double(options.nRowsD);
            gridValues = reshape(kuvaRot, gridRows, []);
            rowValues = zeros(detectorRows, size(gridValues, 2), 'like', kuvaRot);
            for detectorRow = 0:detectorRows-1
                position = transA * detectorRow + transB;
                lower = floor(position); fraction = position - lower;
                if lower >= 0 && lower < gridRows
                    rowValues(detectorRow + 1, :) = rowValues(detectorRow + 1, :) + ...
                        gridValues(lower + 1, :) .* cast(1 - fraction, 'like', kuvaRot);
                end
                upper = lower + 1;
                if upper >= 0 && upper < gridRows
                    rowValues(detectorRow + 1, :) = rowValues(detectorRow + 1, :) + ...
                        gridValues(upper + 1, :) .* cast(fraction, 'like', kuvaRot);
                end
            end
            kuvaRot = reshape(rowValues, detectorRows, size(kuvaRot, 2));

            gridRows = size(kuvaRot, 2);
            detectorCols = double(options.nColsD);
            gridValues = reshape(kuvaRot.', gridRows, []);
            colValues = zeros(detectorCols, size(gridValues, 2), 'like', kuvaRot);
            for detectorCol = 0:detectorCols-1
                position = axialA * detectorCol + axialB;
                lower = floor(position); fraction = position - lower;
                if lower >= 0 && lower < gridRows
                    colValues(detectorCol + 1, :) = colValues(detectorCol + 1, :) + ...
                        gridValues(lower + 1, :) .* cast(1 - fraction, 'like', kuvaRot);
                end
                upper = lower + 1;
                if upper >= 0 && upper < gridRows
                    colValues(detectorCol + 1, :) = colValues(detectorCol + 1, :) + ...
                        gridValues(upper + 1, :) .* cast(fraction, 'like', kuvaRot);
                end
            end
            kuvaRot = reshape(colValues, detectorCols, size(kuvaRot, 1)).';
            kuvaRot = kuvaRot * (options.dx(ii) / options.dx(1)) / (options.dy(1) * options.dz(1));
            outputFP(:, :, kk) = outputFP(:, :, kk) + kuvaRot;
            u1 = u1 + 1;
        end
    end
    % 6. Apply central-ray, normalization and additional correction weights.
    outputFP = outputFP .* type6CentralRayWeights(options, viewIndices);
    outputFP = outputFP .* type6MeasurementWeights(options, norm_input, corr_input, viewIndices);
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
if options.projector_type == 6
    outputFP = forwardProjectionType6(options, recApu, loopVar, koko, nMeas, subIter);
    A = [];
    outputFP = outputFP(:);
    return
end
if (~ismac && (options.implementation == 4 || options.implementation == 1))
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
            || options.projector_type == 2 || options.projector_type == 3 || options.projector_type == 21 || options.projector_type == 22 || options.projector_type == 33
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
panel=double(options.blurPlanes2F(view))*dx(1)/dx(volume);
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
