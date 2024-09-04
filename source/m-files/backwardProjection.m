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

if isempty(varargin)
    nCores = 0;
else
    nCores = varargin{1};
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
if projType == 1 || projType == 11 || projType == 21 || projType == 31
    projType = 1;
elseif projType == 2 || projType == 12 || projType == 22 || projType == 32
    projType = 2;
elseif projType == 3 || projType == 13 || projType == 33 || projType == 23
    projType = 3;
end
if options.compute_sensitivity_image
    options.use_raw_data = true;
    [x, ~, z, ~] = get_coordinates(options);
    options.use_raw_data = false;
end

if options.projector_type == 6
    fProj = reshape(input, options.nRowsD, options.nColsD, koko);
    apuBP2 = zeros(options.Nx(1) * options.Ny(1) * options.Nz(1), koko, options.cType);
    u1 = options.ub;
    for kk = 1 : koko
        apuBP = zeros(options.Nx(1), options.Ny(1), options.Nz(1), options.cType);
        kuvaRot = fProj(:, :, kk);
        kuvaRot = permute(kuvaRot, [2, 1, 3]);
        apu = kuvaRot;
        uu = 1;
        for ll = 1 : options.Ny
            kuvaRot(:,:,uu) = conv2(apu, options.gFilter(:, :, ll, u1),'same');
            uu = uu + 1;
        end
        kuvaRot = kuvaRot(:, :, options.blurPlanes(u1):end);
        kuvaRot = permute(kuvaRot, [3, 2, 1]);
        apuBP(options.blurPlanes(u1):end, :, :) = kuvaRot;
        apuBP = imrotate(apuBP, options.angles(u1), 'bilinear','crop');
        apuBP2(:, kk) = apuBP(:);
        u1 = u1 + 1;
    end
    output = sum(apuBP2, 2);
    output(output < options.epps & output >= 0) = options.epps;
    if noSensIm == 0
        apuBP2 = zeros(options.Nx(1) * options.Ny(1) * options.Nz(1), koko, options.cType);
        u1 = options.ub;
        for kk = 1 : koko
            apuSumm = zeros(options.Nx(1), options.Ny(1), options.Nz(1), options.cType);
            kuvaRot = ones(options.nColsD, options.nRowsD, options.Ny, options.cType);
            apu = kuvaRot(:,:,1);
            for ll = 1 : options.Ny
                kuvaRot(:,:,ll) = conv2(apu, options.gFilter(:, :, ll, u1),'same');
            end
            kuvaRot = kuvaRot(:, :, options.blurPlanes(u1):end);
            kuvaRot = permute(kuvaRot, [3, 2, 1]);
            apuSumm(options.blurPlanes(u1):end, :, :) = kuvaRot;
            apuSumm = imrotate(apuSumm, options.angles(u1), 'bilinear', 'crop');
            apuBP2(:, kk) = apuSumm(:);
            u1 = u1 + 1;
        end
        sensIm = sum(apuBP2, 2);
        sensIm(sensIm < options.epps) = 1;
    end
    options.ub = options.ub + koko;
elseif options.implementation == 4
    for ii = 1 : options.nMultiVolumes + 1
        if outputCell
            if inputCell
                if options.useSingles
                    [output{ii}, sensIm{ii}] = projector_mexSingle( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ...
                        options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...
                        options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ...
                        options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, input{ii}, noSensIm, 2, ...
                        options.dPitchX, options.dPitchY, options.nProjections);
                else
                    [output{ii}, sensIm{ii}] = projector_mex( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ...
                        options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...
                        options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ...
                        options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, input{ii}, noSensIm, 2, ...
                        options.dPitchX, options.dPitchY, options.nProjections);
                end
            else
                if options.useSingles
                    [output{ii}, sensIm{ii}] = projector_mexSingle( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ...
                        options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...
                        options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ...
                        options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, input, noSensIm, 2, ...
                        options.dPitchX, options.dPitchY, options.nProjections);
                else
                    [output{ii}, sensIm{ii}] = projector_mex( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ...
                        options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...
                        options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ...
                        options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, input, noSensIm, 2, ...
                        options.dPitchX, options.dPitchY, options.nProjections);
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
                    options.dPitchX, options.dPitchY, options.nProjections);
            else
                [output, sensIm] = projector_mex( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ...
                    options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...
                    options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ...
                    options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, input, noSensIm, 2, ...
                    options.dPitchX, options.dPitchY, options.nProjections);
            end
            if options.use_psf
                output = computeConvolution(output, options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.gaussK);
                if numel(sensIm) > 1
                    sensIm = computeConvolution(sensIm, options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.gaussK);
                end
            end
        end
    end
elseif options.implementation == 2 || options.implementation == 3 || options.implementation == 5
    if ~isfield(options,'orthTransaxial') && (options.projector_type == 2 || options.projector_type == 3 || options.projector_type == 22 || options.projector_type == 33)
        if options.projector_type == 3 || options.projector_type == 33
            options.orthTransaxial = true;
        elseif (options.projector_type == 2 || options.projector_type == 22) && isfield(options,'tube_width_xy') && options.tube_width_xy > 0
            options.orthTransaxial = true;
        else
            options.orthTransaxial = false;
        end
    end
    if ~isfield(options,'orthAxial') && (options.projector_type == 2 || options.projector_type == 3 || options.projector_type == 22 || options.projector_type == 33)
        if options.projector_type == 3 || options.projector_type == 33
            options.orthAxial = true;
        elseif (options.projector_type == 2 || options.projector_type == 22) && isfield(options,'tube_width_z') && options.tube_width_z > 0
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
    elseif options.projector_type == 6
        header_directory = '';
    else
        error('Invalid projector for OpenCL')
    end
    if numel(nMeas) == 1
        nMeas = [0;nMeas];
    end
    if options.projector_type == 5 || options.projector_type == 15 || options.projector_type == 45
        input = reshape(input, options.nRowsD, options.nColsD, nMeas(subIter + 2) - nMeas(subIter + 1));
        [input, meanBP] = computeIntegralImage(input, options.meanBP);
    end
    if ~isa(input,'single')
        input = single(input);
    end
    [output, sensIm] = OpenCL_matrixfree_multi_gpu( options.Nx, options.Ny, options.Nz, options.dx, options.dy, options.dz, options.bx, options.by, options.bz, ...
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
            alku = prod(options.N(kk)) + 1;
        end
        if numel(sensIm) > 1
            temp = sensIm;
            sensIm = cell(options.nMultiVolumes + 1, 1);
            alku = 1;
            for kk = 1 : options.nMultiVolumes + 1
                sensIm{kk} = temp(alku : alku - 1 + prod(options.N(kk)));
                alku = prod(options.N(kk));
            end
        end
    end
else
    error('Unsupported implementation! This is unintended behavior!')
end