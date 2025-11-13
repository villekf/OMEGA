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

if isempty(varargin)
    nCores = 0;
else
    nCores = varargin{1};
end
if iscell(recApu)
    useCell = true;
else
    useCell = false;
end
projType = options.projector_type;
if projType == 1 || (projType >= 10 && projType < 20)
    projType = 1;
elseif projType == 2 || (projType >= 20 && projType < 30)
    projType = 2;
elseif projType == 3 || (projType >= 30 && projType < 40)
    projType = 3;
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
    outputFP = zeros(options.nRowsD, options.nColsD, koko, options.cType);
    for ii = loopVar
        voxelSize = options.FOVa_x / single(options.Nx(ii));
        u1 = options.uu;
        apuArr = reshape(recApu, options.Nx(ii), options.Ny(ii), options.Nz(ii));

        for kk = 1:koko % todo add support for non-square fov
            kuvaRot = apuArr;
            panelTilt = options.swivelAngles(u1) - options.angles(u1) + 180;
            panelShift = options.radiusPerProj(u1) * sind(panelTilt) / voxelSize;
            PSFshift = (options.FOVa_y/2 - (options.radiusPerProj(u1) * cosd(panelTilt) - options.CORtoDetectorSurface)) / voxelSize;

            % 1. Rotate the image
            kuvaRot = imrotate(kuvaRot, -options.swivelAngles(u1), 'bilinear','crop');
            %volume3Dviewer(kuvaRot)
            
            % 2. Translate the image
            kuvaRot = imtranslate(kuvaRot, [-panelShift, 0, 0]);
            %volume3Dviewer(kuvaRot)

            % 3. Attenuation correction


            % 4. Convolve with detector PSF
            PSF = imtranslate(options.gFilter, [0, 0, options.blurPlanes(u1)], FillValues=0);
            PSF = PSF(:, :, 1:options.Nx(1));
            PSF = permute(PSF, [3 2 1]);
            for jj = 1 : options.Ny(1)
                kuvaRot(jj, :, :) = conv2(squeeze(kuvaRot(jj, :, :)), squeeze(PSF(jj,:,:)), 'same');
            end
            
            % 5. Sum
            kuvaRot = sum(kuvaRot, 1);
            kuvaRot = squeeze(kuvaRot);
            kuvaRot = kuvaRot / double(options.Nx(1));
            outputFP(:, :, kk) = outputFP(:, :, kk) + kuvaRot;
            u1 = u1 + 1;
        end
    end
end

if (~ismac && (options.implementation == 4 || options.implementation == 1)) || (ismac && options.projector_type == 6)
    if options.projector_type == 6
        outputFP = forwardProjectionType6(options, recApu, loopVar, koko);
        A = [];
        outputFP = outputFP(:);
        outputFP(outputFP < options.epps) = options.epps;
    else
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
                                options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, recApu{ii}, true, 1, options.dPitchX, ...
                                options.dPitchY, options.nProjections);
                        else
                            [outputFP, A] = projector_mex( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ...
                                options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...
                                options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ...
                                options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, recApu{ii}, true, 1, options.dPitchX, ...
                                options.dPitchY, options.nProjections);
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
                            options.TOF_bins, options.verbose, nCores, options.use_raw_data, 0, options.listmode, projType, subIter, nMeas, lor2, lor1, summa, options.dPitchY, options.nProjections); % 42
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
                                options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, recApu, true, 1, options.dPitchX, ... % 40
                                options.dPitchY, options.nProjections); % 43
                        else
                            [outputFP, A] = projector_mex( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ... % 8
                                options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ... % 18
                                options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ... % 28
                                options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, recApu, true, 1, options.dPitchX, ... % 40
                                options.dPitchY, options.nProjections); % 43
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
                                options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, recApu{ii}, true, 1, options.dPitchX, ...
                                options.dPitchY, options.nProjections);
                        else
                            [apu, A] = projector_mex( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ...
                                options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...
                                options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ...
                                options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, recApu{ii}, true, 1, options.dPitchX, ...
                                options.dPitchY, options.nProjections);
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
                                options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, recApu, true, 1, options.dPitchX, ...
                                options.dPitchY, options.nProjections);
                        else
                            [apu, A] = projector_mex( options, options.Nx(ii), options.Ny(ii), options.Nz(ii), options.dx(ii), options.dy(ii), options.dz(ii), options.bx(ii), ...
                                options.by(ii), options.bz(ii), z, x, options.size_x, options.vaimennus, norm_input, koko, options.attenuation_correction, options.normalization_correction, ...
                                options.scatter, corr_input, options.global_correction_factor, xy_index, z_index, L_input, options.det_per_ring, TOF, options.sigma_x, options.TOFCenter, ...
                                options.TOF_bins, options.verbose, nCores, options.use_raw_data, 1, options.listmode, projType, subIter, nMeas, options.epps, recApu, true, 1, options.dPitchX, ...
                                options.dPitchY, options.nProjections);
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
        if options.projector_type == 5 || options.projector_type == 51 || options.projector_type == 54
            kopio = recApu;
            for kk = 1 : size(recApu,1)
                recApu{kk} = reshape(recApu{kk}, options.Nx(kk), options.Ny(kk), options.Nz(kk));
                kopio{kk} = reshape(kopio{kk}, options.Nx(kk), options.Ny(kk), options.Nz(kk));
                kopio{kk} = permute(kopio{kk}, [1 3 2]);
                [recApu{kk}, meanFP] = computeIntegralImage(recApu{kk}, options.meanFP);
                [kopio{kk}, meanFP] = computeIntegralImage(kopio{kk}, options.meanFP);
            end
        end
        options.x0 = cell2mat(recApu);
        if options.projector_type == 5 || options.projector_type == 51 || options.projector_type == 54
            options.x0 = [options.x0;cell2mat(kopio)];
        end
    else
        if options.use_psf
            recApu = computeConvolution(recApu, options, options.Nx(1), options.Ny(1), options.Nz(1), options.gaussK);
        end
        if options.projector_type == 5 || options.projector_type == 51 || options.projector_type == 54
            recApu = reshape(recApu, options.Nx(1), options.Ny(1), options.Nz(1));
            kopio = permute(recApu, [1 3 2]);
            recApu = permute(recApu, [2 3 1]);
            [recApu, meanFP] = computeIntegralImage(recApu, options.meanFP);
            [kopio, meanFP] = computeIntegralImage(kopio, options.meanFP);
            recApu = [recApu(:);kopio(:)];
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
    elseif options.projector_type == 4 || options.projector_type == 41 || options.projector_type == 14 || options.projector_type == 45
        kernel_file = 'projectorType4.cl';
        kernel_path = which(kernel_file);
        kernel_path = strrep(kernel_path, '\', '/');
        kernel_path = strrep(kernel_path, '.cl', '');
        header_directory = strrep(kernel_path,'projectorType4','');
    elseif options.projector_type == 5 || options.projector_type == 51 || options.projector_type == 15 || options.projector_type == 54
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
    [outputFP] = OpenCL_matrixfree_multi_gpu( options.Nx, options.Ny, options.Nz, options.dx, options.dy, options.dz, options.bx, options.by, options.bz, ...
        z, x, options.nRowsD, options.verbose, options.LL, options.TOF, ... % 15
        TOFSize, options.sigma_x, options.TOFCenter, options.TOF_bins, options.platform, options.use_raw_data, options.use_psf, header_directory, options.vaimennus, ... % 24
        options.normalization, nMeas, options.attenuation_correction, options.normalization_correction, 1, options.subsets, options.epps, options.xy_index, ...
        options.z_index, crystal_size_z, ... % 34
        options.x_center, options.y_center, options.z_center, single(0), 0, options.projector_type, n_rays, n_rays3D, ... % 42
        options, single(0), options.partitions, options.use_64bit_atomics, options.bmin, options.bmax, options.Vmax, options.V, options.gaussK, 1, 1); % 51
    end
end