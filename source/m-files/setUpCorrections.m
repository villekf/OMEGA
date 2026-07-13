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

            options.FOVa_x = ([FOVxM0, FOVxM1, FOVxM2, FOVxM3, FOVxM4, FOVxM5, FOVxM6]);
            options.Nx = uint32([NxM0, NxM1, NxM2, NxM3, NxM4, NxM5, NxM6]);

            options.FOVa_y = [FOVyM0, FOVyM1, FOVyM2, FOVyM3, FOVyM4, FOVyM5, FOVyM6];
            options.Ny = uint32([NyM0, NyM1, NyM2, NyM3, NyM4, NyM5, NyM6]);

            FOVzM3 = options.axial_fov;
            FOVzM4 = options.axial_fov;
            FOVzM5 = options.axial_fov;
            FOVzM6 = options.axial_fov;
            options.axial_fov = [FOVzM0, FOVzM1, FOVzM2, FOVzM3, FOVzM4, FOVzM5, FOVzM6];
            NzM3 = round(FOVzM3 / dzM);
            NzM4 = round(FOVzM4 / dzM);
            NzM5 = round(FOVzM5 / dzM);
            NzM6 = round(FOVzM6 / dzM);
            options.Nz = uint32([NzM0, NzM1, NzM2, NzM3, NzM4, NzM5, NzM6]);
        elseif options.transaxialEFOV && ~options.axialEFOV % Only the transaxial FOV is extended
            options.nMultiVolumes = 4;

            options.FOVa_x = [FOVxM0, FOVxM3, FOVxM4, FOVxM5, FOVxM6];
            options.Nx = uint32([NxM0, NxM3, NxM4, NxM5, NxM6]);

            options.FOVa_y = [FOVyM0, FOVyM3, FOVyM4, FOVyM5, FOVyM6];
            options.Ny = uint32([NyM0, NyM3, NyM4, NyM5, NyM6]);

            % z: not extended
            NzM = round(options.Nz * options.multiResolutionScale);
            options.axial_fov = [options.axialFOVOrig, options.axialFOVOrig, options.axialFOVOrig, options.axialFOVOrig, options.axialFOVOrig];
            options.Nz = uint32([NzM0, NzM, NzM, NzM, NzM]);
        else % Only axial FOV is extended
            options.nMultiVolumes = 2;

            options.FOVa_x = [FOVxM0, FOVxM1, FOVxM2];
            options.Nx = uint32([NxM0, NxM1, NxM2]);

            options.FOVa_y = [FOVyM0, FOVyM1, FOVyM2];
            options.Ny = uint32([NyM0, NyM1, NyM2]);

            options.axial_fov = [FOVzM0, FOVzM1, FOVzM2];
            options.Nz = uint32([NzM0, NzM1, NzM2]);
        end

        if options.axialEFOV
            axialExtension = sum(options.axial_fov(1:3)) / options.axialFOVOrig * 100;
        else
            axialExtension = 100;
        end
        if options.transaxialEFOV
            if options.nMultiVolumes == 6
                transaxialExtensionX = (options.FOVa_x(4) + options.FOVa_x(1) + options.FOVa_x(5)) / options.FOVxOrig * 100;
                transaxialExtensionY = (options.FOVa_y(6) + options.FOVa_y(1) + options.FOVa_y(7)) / options.FOVyOrig * 100;
            else
                transaxialExtensionX = (options.FOVa_x(2) + options.FOVa_x(1) + options.FOVa_x(3)) / options.FOVxOrig * 100;
                transaxialExtensionY = (options.FOVa_y(4) + options.FOVa_y(1) + options.FOVa_y(5)) / options.FOVyOrig * 100;
            end
        else
            transaxialExtensionX = 100;
            transaxialExtensionY = 100;
        end

        disp(['Axial FOV extension is ' num2str(axialExtension) ' % (z) of the original'])
        disp(['Transaxial FOV extension is ' num2str(transaxialExtensionX) ' % (x), ' num2str(transaxialExtensionY) ' % (y) of the original'])

        if options.useMaskBP && isfield(options, 'maskBP') && (iscell(options.maskBP) || numel(options.maskBP) > 1)
            if isfield(options, 'partitions')
                if numel(options.partitions) > 1
                    partitions = numel(options.partitions);
                elseif isempty(options.partitions)
                    partitions = 1;
                else
                    partitions = options.partitions;
                end
            else
                partitions = 1;
            end
            maskBP = options.maskBP;
            if iscell(maskBP)
                for tt = 1 : numel(maskBP)
                    maskVol = maskBP{tt};
                    if size(maskVol, 3) == 1 && options.NzFull > 1
                        maskVol = repmat(maskVol, 1, 1, double(options.NzFull));
                    end
                    options.maskBP{tt} = uint8(multiResReshape(uint8(maskVol), options) > 0);
                end
            else
                maskVol = maskBP;
                if size(maskVol, 3) == 1 && options.NzFull > 1
                    maskVol = repmat(maskVol, 1, 1, double(options.NzFull));
                end
                maskBP = uint8(multiResReshape(uint8(maskVol), options) > 0);
                if partitions > 1
                    options.maskBP = repmat({maskBP}, partitions, 1);
                else
                    options.maskBP = maskBP;
                end
            end
            options.maskBPZ = max(options.Nz);
        end

        options.x0 = multiResReshape(options.x0, options);

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

function outputVol = multiResReshape(inputVol, options)
    scale = double(options.multiResolutionScale);
    maskInput = islogical(inputVol) || isa(inputVol, 'uint8');
    resizeMethod = 'linear';
    if maskInput
        resizeMethod = 'nearest';
    end
    if options.nMultiVolumes == 6
        lowResSize = [
            double(options.Nx(4)) + round(double(options.Nx(1)) * scale) + double(options.Nx(5)), ...
            double(options.Ny(6)) + round(double(options.Ny(1)) * scale) + double(options.Ny(7)), ...
            double(options.Nz(2)) + round(double(options.Nz(1)) * scale) + double(options.Nz(3))
        ];
        lowResSize = max(lowResSize, double([max(options.Nx(2:end)), max(options.Ny(2:end)), max(options.Nz(2:end))]));
        highResSize = [
            round(double(options.Nx(4)) / scale) + double(options.Nx(1)) + round(double(options.Nx(5)) / scale), ...
            round(double(options.Ny(6)) / scale) + double(options.Ny(1)) + round(double(options.Ny(7)) / scale), ...
            round(double(options.Nz(2)) / scale) + double(options.Nz(1)) + round(double(options.Nz(3)) / scale)
        ];

        if exist('imresize3','file') == 2
            highResVol = single(imresize3(inputVol, highResSize, resizeMethod));
            lowResVol = single(imresize3(inputVol, lowResSize, resizeMethod));
        else
            [row, col, plane] = ndgrid( ...
                linspace(1, size(inputVol, 1), highResSize(1)), ...
                linspace(1, size(inputVol, 2), highResSize(2)), ...
                linspace(1, size(inputVol, 3), highResSize(3)));
            highResVol = single(interp3(single(inputVol), col, row, plane, resizeMethod));
            [row, col, plane] = ndgrid( ...
                linspace(1, size(inputVol, 1), lowResSize(1)), ...
                linspace(1, size(inputVol, 2), lowResSize(2)), ...
                linspace(1, size(inputVol, 3), lowResSize(3)));
            lowResVol = single(interp3(single(inputVol), col, row, plane, resizeMethod));
        end

        highStart = [
            round(double(options.Nx(4)) / scale) + 1, ...
            round(double(options.Ny(6)) / scale) + 1, ...
            round(double(options.Nz(2)) / scale) + 1
        ];
        lowStart = [
            double(options.Nx(4)) + 1, ...
            double(options.Ny(6)) + 1, ...
            double(options.Nz(2)) + 1
        ];
        lowMainSize = [
            round(double(options.Nx(1)) * scale), ...
            round(double(options.Ny(1)) * scale), ...
            round(double(options.Nz(1)) * scale)
        ];

        vol0 = highResVol(highStart(1) : highStart(1) + double(options.Nx(1)) - 1, highStart(2) : highStart(2) + double(options.Ny(1)) - 1, highStart(3) : highStart(3) + double(options.Nz(1)) - 1);
        vol1 = lowResVol(lowStart(1) : lowStart(1) + double(options.Nx(2)) - 1, lowStart(2) : lowStart(2) + double(options.Ny(2)) - 1, 1 : double(options.Nz(2)));
        vol2 = lowResVol(lowStart(1) : lowStart(1) + double(options.Nx(3)) - 1, lowStart(2) : lowStart(2) + double(options.Ny(3)) - 1, lowStart(3) + lowMainSize(3) : lowStart(3) + lowMainSize(3) + double(options.Nz(3)) - 1);
        vol3 = lowResVol(1 : double(options.Nx(4)), 1 : double(options.Ny(4)), 1 : double(options.Nz(4)));
        vol4 = lowResVol(lowStart(1) + lowMainSize(1) : lowStart(1) + lowMainSize(1) + double(options.Nx(5)) - 1, 1 : double(options.Ny(5)), 1 : double(options.Nz(5)));
        vol5 = lowResVol(lowStart(1) : lowStart(1) + double(options.Nx(6)) - 1, 1 : double(options.Ny(6)), 1 : double(options.Nz(6)));
        vol6 = lowResVol(lowStart(1) : lowStart(1) + double(options.Nx(7)) - 1, lowStart(2) + lowMainSize(2) : lowStart(2) + lowMainSize(2) + double(options.Ny(7)) - 1, 1 : double(options.Nz(7)));

        outputVol = [vol0(:); vol1(:); vol2(:); vol3(:); vol4(:); vol5(:); vol6(:)];
        if maskInput
            outputVol = uint8(outputVol > 0);
        else
            outputVol = single(outputVol);
        end
    elseif options.nMultiVolumes == 4
        lowResSize = [
            double(options.Nx(2)) + round(double(options.Nx(1)) * scale) + double(options.Nx(3)), ...
            double(options.Ny(4)) + round(double(options.Ny(1)) * scale) + double(options.Ny(5)), ...
            round(double(options.Nz(1)) * scale)
        ];
        lowResSize = max(lowResSize, double([max(options.Nx(2:end)), max(options.Ny(2:end)), max(options.Nz(2:end))]));
        highResSize = [
            round(double(options.Nx(2)) / scale) + double(options.Nx(1)) + round(double(options.Nx(3)) / scale), ...
            round(double(options.Ny(4)) / scale) + double(options.Ny(1)) + round(double(options.Ny(5)) / scale), ...
            double(options.Nz(1))
        ];

        if exist('imresize3','file') == 2
            highResVol = single(imresize3(inputVol, highResSize, resizeMethod));
            lowResVol = single(imresize3(inputVol, lowResSize, resizeMethod));
        else
            [row, col, plane] = ndgrid( ...
                linspace(1, size(inputVol, 1), highResSize(1)), ...
                linspace(1, size(inputVol, 2), highResSize(2)), ...
                linspace(1, size(inputVol, 3), highResSize(3)));
            highResVol = single(interp3(single(inputVol), col, row, plane, resizeMethod));
            [row, col, plane] = ndgrid( ...
                linspace(1, size(inputVol, 1), lowResSize(1)), ...
                linspace(1, size(inputVol, 2), lowResSize(2)), ...
                linspace(1, size(inputVol, 3), lowResSize(3)));
            lowResVol = single(interp3(single(inputVol), col, row, plane, resizeMethod));
        end

        highStart = [
            round(double(options.Nx(2)) / scale) + 1, ...
            round(double(options.Ny(4)) / scale) + 1, ...
            1
        ];
        lowStart = [
            double(options.Nx(2)) + 1, ...
            double(options.Ny(4)) + 1, ...
            1
        ];
        lowMainSize = [
            round(double(options.Nx(1)) * scale), ...
            round(double(options.Ny(1)) * scale), ...
            round(double(options.Nz(1)) * scale)
        ];

        vol0 = highResVol(highStart(1) : highStart(1) + double(options.Nx(1)) - 1, highStart(2) : highStart(2) + double(options.Ny(1)) - 1, 1 : double(options.Nz(1)));
        vol1 = lowResVol(1 : double(options.Nx(2)), 1 : double(options.Ny(2)), 1 : double(options.Nz(2)));
        vol2 = lowResVol(lowStart(1) + lowMainSize(1) : lowStart(1) + lowMainSize(1) + double(options.Nx(3)) - 1, 1 : double(options.Ny(3)), 1 : double(options.Nz(3)));
        vol3 = lowResVol(lowStart(1) : lowStart(1) + double(options.Nx(4)) - 1, 1 : double(options.Ny(4)), 1 : double(options.Nz(4)));
        vol4 = lowResVol(lowStart(1) : lowStart(1) + double(options.Nx(5)) - 1, lowStart(2) + lowMainSize(2) : lowStart(2) + lowMainSize(2) + double(options.Ny(5)) - 1, 1 : double(options.Nz(5)));

        outputVol = [vol0(:); vol1(:); vol2(:); vol3(:); vol4(:)];
        if maskInput
            outputVol = uint8(outputVol > 0);
        else
            outputVol = single(outputVol);
        end
    elseif options.nMultiVolumes == 2
        lowResSize = [
            round(double(options.Nx(1)) * scale), ...
            round(double(options.Ny(1)) * scale), ...
            double(options.Nz(2)) + round(double(options.Nz(1)) * scale) + double(options.Nz(3))
        ];
        lowResSize = max(lowResSize, double([max(options.Nx(2:end)), max(options.Ny(2:end)), max(options.Nz(2:end))]));
        highResSize = [
            double(options.Nx(1)), ...
            double(options.Ny(1)), ...
            round(double(options.Nz(2)) / scale) + double(options.Nz(1)) + round(double(options.Nz(3)) / scale)
        ];

        if exist('imresize3','file') == 2
            highResVol = single(imresize3(inputVol, highResSize, resizeMethod));
            lowResVol = single(imresize3(inputVol, lowResSize, resizeMethod));
        else
            [row, col, plane] = ndgrid( ...
                linspace(1, size(inputVol, 1), highResSize(1)), ...
                linspace(1, size(inputVol, 2), highResSize(2)), ...
                linspace(1, size(inputVol, 3), highResSize(3)));
            highResVol = single(interp3(single(inputVol), col, row, plane, resizeMethod));
            [row, col, plane] = ndgrid( ...
                linspace(1, size(inputVol, 1), lowResSize(1)), ...
                linspace(1, size(inputVol, 2), lowResSize(2)), ...
                linspace(1, size(inputVol, 3), lowResSize(3)));
            lowResVol = single(interp3(single(inputVol), col, row, plane, resizeMethod));
        end

        highStart = round(double(options.Nz(2)) / scale) + 1;
        lowStart = double(options.Nz(2)) + 1;
        lowMainSize = round(double(options.Nz(1)) * scale);

        vol0 = highResVol(1 : double(options.Nx(1)), 1 : double(options.Ny(1)), highStart : highStart + double(options.Nz(1)) - 1);
        vol1 = lowResVol(1 : double(options.Nx(2)), 1 : double(options.Ny(2)), 1 : double(options.Nz(2)));
        vol2 = lowResVol(1 : double(options.Nx(3)), 1 : double(options.Ny(3)), lowStart + lowMainSize : lowStart + lowMainSize + double(options.Nz(3)) - 1);

        outputVol = [vol0(:); vol1(:); vol2(:)];
        if maskInput
            outputVol = uint8(outputVol > 0);
        else
            outputVol = single(outputVol);
        end
    end
end
