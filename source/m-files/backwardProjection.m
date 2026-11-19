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
% The sensitivity image of list-mode data is the backprojection over all detector pairs of the scanner (or of
% options.xSens/options.zSens), not over the events: the projectors read the detector coordinates (x: x and y of each
% detector, z: z of each ring, layers concatenated) from x and z. The CPU MEX (implementation 4) enumerates the
% nLayers^2 * det_per_ring^2 * rings^2 detector pairs itself (det_per_ring and rings per layer) and ignores koko and the
% input
if options.listmode > 0 && options.compute_sensitivity_image && isfield(options,'xSens') && isfield(options,'zSens') && ~isempty(options.xSens) && ~isempty(options.zSens) && ~options.SPECT
    x = cast(options.xSens, options.cType);
    z = cast(options.zSens, options.cType);
	options.det_per_ring = numel(x) / 2 / options.nLayers; % x and z contain the coordinates of all layers
	options.rings = numel(z) / options.nLayers;
elseif options.listmode > 0 && options.compute_sensitivity_image && ~options.useIndexBasedReconstruction
    options.use_raw_data = true;
    [x, ~, z, ~] = get_coordinates(options);
    options.use_raw_data = false;
end

function output = backwardProjectionType6(options, input, koko)
    % Matches the GPU backProjectionType6Helper (functions.hpp): smear, PSF convolution, attenuation, translation back
    % (same axis, opposite direction to the forward projection), rotation back and the backprojection mask
    fProj = reshape(input, options.nRowsD, options.nColsD, koko);
    output = zeros(options.Nx(1) * options.Ny(1) * options.Nz(1), 1, options.cType);
    u1 = options.ub;
    for kk = 1 : koko
        % Rounded integer shifts, as the int32 distInt/distInt2 of the GPU code
        blurPlane = round(double(options.blurPlanes(u1)));
        panelShift = round(double(options.blurPlanes2(u1)));

        % 1. Smear the input FP across the image volume
        kuvaRot = fProj(:, :, kk);
        kuvaRot = permute(kuvaRot, [2 1 3]);
        kuvaRot = repmat(kuvaRot, 1, 1, options.Nx(1));
        kuvaRot = permute(kuvaRot, [3, 2, 1]);
        kuvaRot = kuvaRot / double(options.Nx(1));

        % 2. Convolve with detector PSF. The forward projection applies the attenuation before the convolution, so
        % the convolution has to be done first here
        kuvaRot = convolvePSFType6(kuvaRot, options.gFilter, blurPlane, options.Nx(1), options.Ny(1), options.Nz(1));

        % 3. Attenuation correction. The same attenuation factors as in the forward projection
        if options.attenuation_correction
            attenuationImage = attenuationType6(options.vaimennus, options.Nx(1), options.Ny(1), options.Nz(1), options.dx(1), options.swivelAngles(u1), panelShift);
            kuvaRot = kuvaRot .* attenuationImage;
        end

        % 4. Translate the image back
        if panelShift ~= 0
            kuvaRot = translateType6(kuvaRot, 0, panelShift, 0);
        end

        % 5. Rotate the image back
        kuvaRot = rotateType6(kuvaRot, options.swivelAngles(u1));

        output = output + kuvaRot(:);
        u1 = u1 + 1;
    end
    % 6. Backprojection mask
    if isfield(options, 'useMaskBP') && options.useMaskBP && isfield(options, 'maskBP') && numel(options.maskBP) > 1
        maskBP = double(reshape(options.maskBP, options.Nx(1), options.Ny(1), []));
        if size(maskBP, 3) == 1
            maskBP = repmat(maskBP, 1, 1, options.Nz(1));
        end
        output = output .* maskBP(:);
    end
end

if (options.implementation == 4 && ~ismac) || (ismac && options.projector_type == 6)
    if options.projector_type == 6
        output = backwardProjectionType6(options, input, koko);
        output(output < options.epps & output >= 0) = options.epps;
        if noSensIm == 0
            sensIm = backwardProjectionType6(options, ones(size(input)), koko);
            sensIm(sensIm < options.epps) = 1;
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
    options.x0 = zeros(options.Nx(1) * options.Ny(1) * options.Nz(1), 1, 'single');
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
    elseif options.projector_type == 5 || options.projector_type == 55 || options.projector_type == 51 || options.projector_type == 15 || options.projector_type == 54 || options.projector_type == 45
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
    if options.projector_type == 5 || options.projector_type == 55 || options.projector_type == 15 || options.projector_type == 45
        input = reshape(input, options.nRowsD, options.nColsD, nMeas(subIter + 2) - nMeas(subIter + 1));
        % The backprojection is the transpose of the forward projection, which multiplies its output by the normalization, so the
        % projections of the current subset are multiplied by the normalization before the integral image is formed, since the
        % box sum cannot weight the individual detector pixels
        if options.CT && options.normalization_correction && numel(options.normalization) > 1
            nDetPix = double(options.nRowsD) * double(options.nColsD);
            normBP5 = options.normalization(double(nMeas(subIter + 1)) * nDetPix + 1 : double(nMeas(subIter + 2)) * nDetPix);
            normBP5 = reshape(normBP5, options.nRowsD, options.nColsD, size(input, 3));
            if isa(input, 'gpuArray')
                normBP5 = gpuArray(cast(normBP5, classUnderlying(input)));
            else
                normBP5 = cast(normBP5, class(input));
            end
            input = input .* normBP5;
        end
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
