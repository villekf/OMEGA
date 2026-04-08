function options = setUpCorrections(options)
%SETUPCORRECTIONS Performs the necessary steps for multi-resolution
%reconstruction and offset correction
%   Computes the volume sizes for the multi-resolution volumes, adjust the
%   FOVs accordingly and resizes the initial value when using
%   multi-resolution reconstruction. For extended FOV adjusts the masks for
%   prior computations to ensure that the priors are only computed for the
%   original volume size. For offset correction, computes the indices
%   affected by the offset correction and computes the offset weights. Both
%   the input and output variables are the param struct from the
%   projectorClass class object.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2023-2025 Ville-Veikko Wettenhovi
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Below is the volume ordering for multi-resolution volume cases, assuming
% both extended transaxial and axial volumes
% 0 is the main volume located in the middle
% Volume 1 is in the front, in the negative region assuming the origin is
% is in (0,0,0)
% If only axial FOV is extended, then volumes 3-6 are omitted completely
% Likewise, if only transaxial FOV is extended volumes 1-2 are omitted
%           ________________________________________
%           |                                       |
%           |                                       |
%           |                  4                    |
%           |                                       |
%           |_______________________________________|
%           |           |               |           |
%           |           |               |           |
%           |           |               |           |
%           |     5     |      1/2      |     6     |
%           |           |               |           |
%           |           |               |           |
%           |___________|_______________|___________|
%           |                                       |
%           |                                       |
%           |                  3                    |
%           |                                       |
%           |_______________________________________|


if options.useEFOV
    nx = options.Nx;
    ny = options.Ny;
    nz = options.Nz;
    options.NxFull = nx;
    options.NyFull = ny;
    options.NzFull = nz;
    if options.useMultiResolutionVolumes
        if ~isfield(options, 'multiResolutionScale')
            warning('No scale value input for multi-resolution reconstruction. Using default value of 1/4 of the original voxel size.')
            options.multiResolutionScale = .25;
        end

        % x: up-down direction in figure above
        % y: left-right direction in figure above
        dxM = options.FOVxOrig / (options.NxOrig * options.multiResolutionScale); % Multiresolution voxel sizes
        dyM = options.FOVyOrig / (options.NyOrig * options.multiResolutionScale);
        dzM = options.axialFOVOrig / (options.NzOrig * options.multiResolutionScale); 
            
        FOVxM0 = options.FOVxOrig; % Size of FOV 0 (main volume)
        FOVyM0 = options.FOVyOrig;
        FOVzM0 = options.axialFOVOrig;

        NxM0 = options.NxOrig;
        NyM0 = options.NyOrig;
        NzM0 = options.NzOrig;

        if options.axialEFOV
            FOVxM1 = options.FOVxOrig;
            FOVxM2 = options.FOVxOrig;
            NxM1 = round(FOVxM1 / dxM);
            NxM2 = round(FOVxM2 / dxM);

            FOVyM1 = options.FOVyOrig;
            FOVyM2 = options.FOVyOrig;
            NyM1 = round(FOVyM1 / dyM);
            NyM2 = round(FOVyM2 / dyM);

            FOVzM1 = (options.axial_fov - options.axialFOVOrig) / 2 - options.eFOVShift(3);
            FOVzM2 = (options.axial_fov - options.axialFOVOrig) / 2 + options.eFOVShift(3);
            NzM1 = round(FOVzM1 / dzM);
            NzM2 = round(FOVzM2 / dzM);
        end

        if options.transaxialEFOV
            FOVxM3 = (options.FOVa_x - options.FOVxOrig) / 2 - options.eFOVShift(1); % Multiresolution FOV size x-direction (volume 3)
            FOVxM4 = (options.FOVa_x - options.FOVxOrig) / 2 + options.eFOVShift(1); % Multiresolution FOV size x-direction (volume 4)
            FOVxM5 = options.FOVxOrig;
            FOVxM6 = options.FOVxOrig;

            NxM3 = round(FOVxM3 / dxM); % Multiresolution amount of voxels x-direction (volume 3)
            NxM4 = round(FOVxM4 / dxM); % Multiresolution amount of voxels x-direction (volume 4)
            NxM5 = round(options.NxOrig * options.multiResolutionScale);
            NxM6 = round(options.NxOrig * options.multiResolutionScale);

            FOVyM3 = options.FOVa_y;
            FOVyM4 = options.FOVa_y;
            FOVyM5 = (options.FOVa_y - options.FOVyOrig) / 2 - options.eFOVShift(2);
            FOVyM6 = (options.FOVa_y - options.FOVyOrig) / 2 + options.eFOVShift(2);

            NyM3 = round(FOVyM3 / dyM);
            NyM4 = round(FOVyM4 / dyM);
            NyM5 = round(FOVyM5 / dyM);
            NyM6 = round(FOVyM6 / dyM);
        end

        if options.axialEFOV && options.transaxialEFOV % Both transaxial and axial FOVs are extended.
            options.nMultiVolumes = 6;

            options.FOVa_x = ([
                FOVxM0;
                FOVxM1;
                FOVxM2;
                FOVxM3;
                FOVxM4;
                FOVxM5;
                FOVxM6;
            ]');

            options.Nx = uint32([
                NxM0;
                NxM1;
                NxM2;
                NxM3;
                NxM4;
                NxM5;
                NxM6;
            ]');

            options.FOVa_y = [
                FOVyM0;
                FOVyM1;
                FOVyM2;
                FOVyM3;
                FOVyM4;
                FOVyM5;
                FOVyM6;
            ]';

            options.Ny = uint32([
                NyM0;
                NyM1;
                NyM2;
                NyM3;
                NyM4;
                NyM5;
                NyM6;
            ]');

            FOVzM3 = options.axial_fov;
            FOVzM4 = options.axial_fov;
            FOVzM5 = options.axial_fov;
            FOVzM6 = options.axial_fov;

            options.axial_fov = [
                FOVzM0;
                FOVzM1;
                FOVzM2;
                FOVzM3;
                FOVzM4;
                FOVzM5;
                FOVzM6;
            ]';

            NzM3 = round(FOVzM3 / dzM);
            NzM4 = round(FOVzM4 / dzM);
            NzM5 = round(FOVzM5 / dzM);
            NzM6 = round(FOVzM6 / dzM);

            options.Nz = uint32([
                NzM0;
                NzM1;
                NzM2;
                NzM3;
                NzM4;
                NzM5;
                NzM6;
            ]');

            disp(['Extended FOV is ' num2str(FOVyM3/options.FOVyOrig * 100) ' % of the original'])
            % If the initial value only covers the original, dense, volume,
            % it needs to be resized to match the extended FOV with
            % multi-resolution
            if size(options.x0, 1) == nx && min(options.x0(:)) ~= max(options.x0(:))
                if exist('imresize3','file') == 2
                    apu = single(imresize3(options.x0,options.multiResolutionScale));
                else
                    apu = single(imresize(options.x0,options.multiResolutionScale));
                    [XX, YY, ZZ] = meshgrid(1:size(apu,1), 1:size(apu,1), 1:round(1/options.multiResolutionScale):size(options.x0,3));
                    apu = interp3(apu, XX, YY, ZZ);
                end
                options.x1 = apu(options.Nx(4) + 1 : options.Nx(4) + options.Nx(2), options.Ny(6) + 1 : options.Ny(6) + options.Ny(2), ...
                    1 : options.Nz(2));
                options.x2 = apu(options.Nx(5) + 1 : options.Nx(5) + options.Nx(3), options.Ny(7) + 1 : options.Ny(7) + options.Ny(3), ...
                    end - options.Nz(2) + 1:end);
                options.x3 = apu(1 : options.Nx(4), :, :);
                if mod(size(apu,1),2) == 0
                    options.x4 = apu(options.Nx(5) + options.Nx(3) : end, :, :);
                else
                    options.x4 = apu(1 + options.Nx(5) + options.Nx(3) : end, :, :);
                end
                options.x5 = apu(options.Nx(4) + 1 : options.Nx(4) + options.Nx(2), 1:options.Ny(6), ...
                    1 : options.Nz(4));
                if mod(size(apu,2),2) == 0
                    options.x6 = apu(options.Nx(5) + 1 : options.Nx(5) + options.Nx(3), options.Ny(7) + options.Ny(3) : end, ...
                        1 : options.Nz(5));
                else
                    options.x6 = apu(options.Nx(5) + 1 : options.Nx(5) + options.Nx(3), 1 + options.Ny(7) + options.Ny(3) : end, ...
                        1 : options.Nz(5));
                end
                options.x0 = single(options.x0(1 + (size(options.x0,1) - options.NxOrig) / 2 : ...
                    (size(options.x0,1) - options.NxOrig) / 2 + options.NxOrig, ...
                    1 + (size(options.x0,2) - options.NyOrig) / 2 : ...
                    (size(options.x0,2) - options.NyOrig) / 2 + options.NyOrig,...
                    1 + (size(options.x0,3) - options.NzOrig) / 2 : ...
                    (size(options.x0,3) - options.NzOrig) / 2 + options.NzOrig));
                if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
                    options.x0 = [options.x0(:);options.x1(:); options.x2(:); options.x3(:); options.x4(:); options.x5(:); options.x6(:)];
                    options = rmfield(options, {'x1','x2','x3','x4','x5','x6'});
                end
            elseif size(options.x0, 1) == options.NxOrig || min(options.x0(:)) == max(options.x0(:))
                val = min(options.x0(:));
                options.x1 = ones(options.Nx(2), options.Ny(2), options.Nz(2), options.cType) * val;
                options.x2 = ones(options.Nx(3), options.Ny(3), options.Nz(3), options.cType) * val;
                options.x3 = ones(options.Nx(4), options.Ny(4), options.Nz(4), options.cType) * val;
                options.x4 = ones(options.Nx(5), options.Ny(5), options.Nz(5), options.cType) * val;
                options.x5 = ones(options.Nx(6), options.Ny(6), options.Nz(6), options.cType) * val;
                options.x6 = ones(options.Nx(7), options.Ny(7), options.Nz(7), options.cType) * val;
                if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
                    options.x0 = [options.x0(:);options.x1(:); options.x2(:); options.x3(:); options.x4(:); options.x5(:); options.x6(:)];
                    options = rmfield(options, {'x1','x2','x3','x4','x5','x6'});
                end
            end
        elseif options.transaxialEFOV && ~options.axialEFOV % Only the transaxial FOV is extended
            options.nMultiVolumes = 4;

            options.FOVa_x = [
                FOVxM0;
                FOVxM3;
                FOVxM4;
                FOVxM5;
                FOVxM6;
            ]';

            options.Nx = uint32([
                NxM0;
                NxM3;
                NxM4;
                NxM5;
                NxM6;
            ]');

            options.FOVa_y = [
                FOVyM0;
                FOVyM3;
                FOVyM4;
                FOVyM5;
                FOVyM6;
            ]';

            options.Ny = uint32([
                NyM0;
                NyM3;
                NyM4;
                NyM5;
                NyM6;
            ]');

            % z: not extended
            NzM = round(options.Nz * options.multiResolutionScale);

            options.axial_fov = [
                options.axialFOVOrig;
                options.axialFOVOrig;
                options.axialFOVOrig;
                options.axialFOVOrig;
                options.axialFOVOrig
            ]';

            options.Nz = uint32([
                NzM0;
                NzM;
                NzM;
                NzM;
                NzM;
            ]');

            disp(['Extended FOV is ' num2str(FOVyM3/options.FOVyOrig * 100) ' % of the original'])
            if size(options.x0, 1) == nx && min(options.x0(:)) ~= max(options.x0(:))
                if exist('imresize3','file') == 2
                    apu = single(imresize3(options.x0,options.multiResolutionScale));
                else
                    apu = single(imresize(options.x0,options.multiResolutionScale));
                    [XX, YY, ZZ] = meshgrid(1:size(apu,1), 1:size(apu,1), 1:round(1/options.multiResolutionScale):size(options.x0,3));
                    apu = interp3(apu, XX, YY, ZZ);
                end
                options.x1 = apu(1 : options.Nx(2), :, :);
                if mod(size(apu,1),2) == 0
                    options.x2 = apu(options.Nx(2) + options.Nx(4) : end, :, :);
                else
                    options.x2 = apu(options.Nx(2) + options.Nx(4) + 1 : end, :, :);
                end
                options.x3 = apu(options.Nx(2) + 1 : options.Nx(2) + options.Nx(4), 1:options.Ny(4), :);
                options.x4 = apu(options.Nx(2) + 1 : options.Nx(2) + options.Nx(5), 1 + options.Ny(2) - options.Ny(4) : end, :);
                options.x0 = single(options.x0(1 + (size(options.x0,1) - options.NxOrig) / 2 : ...
                    (size(options.x0,1) - options.NxOrig) / 2 + options.NxOrig, ...
                    1 + (size(options.x0,2) - options.NyOrig) / 2 : ...
                    (size(options.x0,2) - options.NyOrig) / 2 + options.NyOrig, :));
                if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
                    options.x0 = [options.x0(:);options.x1(:); options.x2(:); options.x3(:); options.x4(:)];
                    options = rmfield(options, {'x1','x2','x3','x4'});
                end
            elseif size(options.x0, 1) == options.NxOrig || min(options.x0(:)) == max(options.x0(:))
                val = min(options.x0(:));
                options.x1 = ones(options.Nx(2), options.Ny(2), options.Nz(2), options.cType) * val;
                options.x2 = ones(options.Nx(3), options.Ny(3), options.Nz(3), options.cType) * val;
                options.x3 = ones(options.Nx(4), options.Ny(4), options.Nz(4), options.cType) * val;
                options.x4 = ones(options.Nx(5), options.Ny(5), options.Nz(5), options.cType) * val;
                if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
                    options.x0 = [options.x0(:);options.x1(:); options.x2(:); options.x3(:); options.x4(:)];
                    options = rmfield(options, {'x1','x2','x3','x4'});
                end
            end
        else % Only axial FOV is extended
            options.nMultiVolumes = 2;

            options.FOVa_x = [
                FOVxM0;
                FOVxM1;
                FOVxM2;
            ]';

            options.Nx = uint32([
                NxM0;
                NxM1;
                NxM2;
            ]');

            options.FOVa_y = [
                FOVyM0;
                FOVyM1;
                FOVyM2;
            ]';

            options.Ny = uint32([
                NyM0;
                NyM1;
                NyM2;
            ]');

            options.axial_fov = [
                FOVzM0;
                FOVzM1;
                FOVzM2;
            ]';

            options.Nz = uint32([
                NzM0;
                NzM1;
                NzM2;
            ]');

            disp(['Extended FOV is ' num2str((FOVzM1 + FOVzM2 + options.axialFOVOrig)/options.axialFOVOrig * 100) ' % of the original'])
            if size(options.x0, 1) == nx && min(options.x0(:)) ~= max(options.x0(:))
                if exist('imresize3','file') == 2
                    apu = single(imresize3(options.x0,options.multiResolutionScale));
                else
                    apu = single(imresize(options.x0,options.multiResolutionScale));
                    [XX, YY, ZZ] = meshgrid(1:size(apu,1), 1:size(apu,1), 1:round(1/options.multiResolutionScale):size(options.x0,3));
                    apu = interp3(apu, XX, YY, ZZ);
                end
                options.x1 = apu(:, :, 1 : options.Nz(2));
                options.x2 = apu(:, :, end - options.Nz(3) + 1:end);
                options.x0 = single(options.x0(:, :, ...
                    1 + (size(options.x0,3) - options.NzOrig) / 2 : ...
                    (size(options.x0,3) - options.NzOrig) / 2 + options.NzOrig));
                if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
                    options.x0 = [options.x0(:);options.x1(:); options.x2(:)];
                    options = rmfield(options, {'x1','x2'});
                end
            elseif size(options.x0, 1) == options.NxOrig || min(options.x0(:)) == max(options.x0(:))
                val = min(options.x0(:));
                options.x1 = ones(options.Nx(2), options.Ny(2), options.Nz(2), options.cType) * val;
                options.x2 = ones(options.Nx(3), options.Ny(3), options.Nz(3), options.cType) * val;
                if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
                    options.x0 = [options.x0(:);options.x1(:); options.x2(:)];
                    options = rmfield(options, {'x1','x2'});
                end
            end
        end
        options.NxPrior = options.Nx(1);
        options.NyPrior = options.Ny(1);
        options.NzPrior = options.Nz(1);
    else
        % This is for non-multiresolution case
        % The idea is that priors/regularization is not computed in the
        % extended region
        if ~isfield(options, 'eFOVIndices') || numel(options.eFOVIndices) < 1
            options.eFOVIndices = zeros(options.Nz,1);
            options.eFOVIndices(((options.Nz - options.NzOrig)/2 + 1 + options.eFOVShift_Nz) : (end - (options.Nz - options.NzOrig)/2 + options.eFOVShift_Nz)) = 1;
        end
        options.eFOVIndices = uint8(options.eFOVIndices);
        options.NzPrior = sum(options.eFOVIndices);
        options.maskPrior = zeros(options.Nx, options.Ny, 'uint8');
        options.maskPrior(((options.Nx - options.NxOrig)/2 + 1 - options.eFOVShift_Nx) : (end - (options.Nx - options.NxOrig)/2 - options.eFOVShift_Nx), ((options.Ny - options.NyOrig)/2 + 1 - options.eFOVShift_Ny) : (end - (options.Ny - options.NyOrig)/2 - options.eFOVShift_Ny)) = uint8(1);
        options.NxPrior = sum(options.maskPrior(:,round(end/2)));
        options.NyPrior = sum(options.maskPrior(round(end/2),:));
        if options.useMaskBP
            options.maskPrior = options.maskPrior - uint8(~options.maskBP);
        end
    end
else
    options.eFOVShift = [0 0 0];
    options.Nx = uint32(options.Nx);
    options.Ny = uint32(options.Ny);
    options.Nz = uint32(options.Nz);
    options.NxFull = options.Nx;
    options.NyFull = options.Ny;
    options.NzFull = options.Nz;
end
if options.offsetCorrection
    options.OffsetLimit = zeros(options.nProjections, 1, options.cType);
    for kk = 1 : options.nProjections
        sx = options.x(1,kk);
        sy = options.x(2,kk);
        sz = options.x(3,kk);
        s = [sx;sy;sz];
        dx = options.x(4,kk);
        dy = options.x(5,kk);
        dz = options.x(6,kk);
        d = [dx;dy;dz];
        ii = (0 : 0.25 : options.nRowsD)' - options.nRowsD / 2;
        if options.pitch
            dx = dx + options.z(1,kk) * ii + options.z(4,kk) * ii;
            dy = dy + options.z(2,kk) * ii + options.z(5,kk) * ii;
        else
            dx = dx + options.z(1,kk) * ii;
            dy = dy + options.z(2,kk) * ii;
        end
        dist = abs((dx - sx) .* (sy) - ((sx) .* (dy - sy))) ./ sqrt((dx - sx).^2 + (dy - sy).^2);
        [~, ind] = min(dist);
        options.OffsetLimit(kk) = (ind) * options.dPitchY/4;
    end
end
end