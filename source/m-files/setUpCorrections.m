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
% Copyright (C) 2023-2024 Ville-Veikko Wettenhovi
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
        if options.axialEFOV && options.transaxialEFOV
            options.nMultiVolumes = 6;
            if ceil(options.Nx * options.multiResolutionScale) == floor(options.NxOrig * options.multiResolutionScale) + ceil((options.Nx - options.NxOrig) / 2 * options.multiResolutionScale)*2
                options.Nx = uint32([options.NxOrig, floor(options.NxOrig * options.multiResolutionScale), floor(options.NxOrig * options.multiResolutionScale), ...
                    ceil((options.Nx - options.NxOrig) / 2 * options.multiResolutionScale), ceil((options.Nx - options.NxOrig) / 2 * options.multiResolutionScale), ...
                    floor(options.NxOrig * options.multiResolutionScale), floor(options.NxOrig * options.multiResolutionScale)]);
            elseif ceil(options.Nx * options.multiResolutionScale) == ceil(options.NxOrig * options.multiResolutionScale) + floor((options.Nx - options.NxOrig) / 2 * options.multiResolutionScale)*2
                options.Nx = uint32([options.NxOrig, ceil(options.NxOrig * options.multiResolutionScale), ceil(options.NxOrig * options.multiResolutionScale), ...
                    floor((options.Nx - options.NxOrig) / 2 * options.multiResolutionScale), floor((options.Nx - options.NxOrig) / 2 * options.multiResolutionScale), ...
                    ceil(options.NxOrig * options.multiResolutionScale), ceil(options.NxOrig * options.multiResolutionScale)]);
            else
                options.Nx = uint32([options.NxOrig, ceil(options.NxOrig * options.multiResolutionScale), ceil(options.NxOrig * options.multiResolutionScale), ...
                    ceil((options.Nx - options.NxOrig) / 2 * options.multiResolutionScale), ceil((options.Nx - options.NxOrig) / 2 * options.multiResolutionScale), ...
                    ceil(options.NxOrig * options.multiResolutionScale), ceil(options.NxOrig * options.multiResolutionScale)]);
            end
            if ceil(options.Ny * options.multiResolutionScale) == floor(options.NyOrig * options.multiResolutionScale) + ceil((options.Ny - options.NyOrig) / 2 * options.multiResolutionScale)*2
                options.Ny = uint32([options.NyOrig, floor(options.NyOrig * options.multiResolutionScale), floor(options.NyOrig * options.multiResolutionScale), ...
                    ceil(options.Ny * options.multiResolutionScale), ceil(options.Ny * options.multiResolutionScale),...
                    ceil((options.Ny - options.NyOrig) / 2 * options.multiResolutionScale), ceil((options.Ny - options.NyOrig) / 2 * options.multiResolutionScale)]);
            elseif ceil(options.Ny * options.multiResolutionScale) == ceil(options.NyOrig * options.multiResolutionScale) + floor((options.Ny - options.NyOrig) / 2 * options.multiResolutionScale)*2
                options.Ny = uint32([options.NyOrig, ceil(options.NyOrig * options.multiResolutionScale), ceil(options.NyOrig * options.multiResolutionScale), ...
                    ceil(options.Ny * options.multiResolutionScale), ceil(options.Ny * options.multiResolutionScale),...
                    floor((options.Ny - options.NyOrig) / 2 * options.multiResolutionScale), floor((options.Ny - options.NyOrig) / 2 * options.multiResolutionScale)]);
            else
                options.Ny = uint32([options.NyOrig, ceil(options.NyOrig * options.multiResolutionScale), ceil(options.NyOrig * options.multiResolutionScale), ...
                    ceil(options.Ny * options.multiResolutionScale), ceil(options.Ny * options.multiResolutionScale),...
                    ceil((options.Ny - options.NyOrig) / 2 * options.multiResolutionScale), ceil((options.Ny - options.NyOrig) / 2 * options.multiResolutionScale)]);
            end
            options.Nz = uint32([options.NzOrig, ceil((options.Nz - options.NzOrig) / 2 * options.multiResolutionScale), ceil((options.Nz - options.NzOrig) / 2 * options.multiResolutionScale), ...
                ceil(options.Nz * options.multiResolutionScale), ceil(options.Nz * options.multiResolutionScale), ...
                ceil(options.Nz * options.multiResolutionScale), ceil(options.Nz * options.multiResolutionScale)]);
            options.FOVa_x = ([options.FOVxOrig, options.FOVxOrig, options.FOVxOrig, ...
                (options.FOVa_x - options.FOVxOrig) / 2, (options.FOVa_x - options.FOVxOrig) / 2, ...
                options.FOVxOrig, options.FOVxOrig]);
            options.FOVa_y = ([options.FOVyOrig, options.FOVyOrig, options.FOVyOrig, ...
                options.FOVa_y, options.FOVa_y, ...
                (options.FOVa_y - options.FOVyOrig) / 2, (options.FOVa_y - options.FOVyOrig) / 2]);
            options.axial_fov = ([options.axialFOVOrig, (options.axial_fov - options.axialFOVOrig) /2, (options.axial_fov - options.axialFOVOrig) /2, ...
                options.axial_fov, options.axial_fov,...
                options.axial_fov, options.axial_fov]);
            %             dimX = options.Nx(2) + options.Nx(4) * 2;
            %             dimY =
            if size(options.x0, 1) == nx
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
            end
        elseif options.transaxialEFOV && ~options.axialEFOV
            options.nMultiVolumes = 4;
            options.Nx = uint32([options.NxOrig, ceil((options.Nx - options.NxOrig) / 2 * options.multiResolutionScale), ceil((options.Nx - options.NxOrig) / 2 * options.multiResolutionScale), ...
                ceil(options.NxOrig * options.multiResolutionScale), ceil(options.NxOrig * options.multiResolutionScale)]);
            options.Ny = uint32([options.NyOrig, ceil(options.Ny * options.multiResolutionScale), ceil(options.Ny * options.multiResolutionScale),...
                ceil((options.Ny - options.NyOrig) / 2 * options.multiResolutionScale), ceil((options.Ny - options.NyOrig) / 2 * options.multiResolutionScale)]);
            options.Nz = uint32([options.Nz, ceil(options.Nz * options.multiResolutionScale), ceil(options.Nz * options.multiResolutionScale), ...
                ceil(options.Nz * options.multiResolutionScale), ceil(options.Nz * options.multiResolutionScale)]);
            options.FOVa_x = ([options.FOVxOrig, (options.FOVa_x - options.FOVxOrig) / 2, (options.FOVa_x - options.FOVxOrig) / 2, ...
                options.FOVxOrig, options.FOVxOrig]);
            options.FOVa_y = ([options.FOVyOrig, options.FOVa_y, options.FOVa_y, ...
                (options.FOVa_y - options.FOVyOrig) / 2, (options.FOVa_y - options.FOVyOrig) / 2]);
            options.axial_fov = ([options.axial_fov, options.axial_fov, options.axial_fov,...
                options.axial_fov, options.axial_fov]);
            if size(options.x0, 1) == nx
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
            end
        else
            options.nMultiVolumes = 2;
            options.Nx = uint32([options.Nx, ceil(options.Nx * options.multiResolutionScale), ceil(options.Nx * options.multiResolutionScale)]);
            options.Ny = uint32([options.Ny, ceil(options.Ny * options.multiResolutionScale), ceil(options.Ny * options.multiResolutionScale)]);
            options.Nz = uint32([options.NzOrig, ceil((options.Nz - options.NzOrig) / 2 * options.multiResolutionScale), ceil((options.Nz - options.NzOrig) / 2 * options.multiResolutionScale)]);
            options.FOVa_x = ([options.FOVa_x, options.FOVa_x, options.FOVa_x]);
            options.FOVa_y = ([options.FOVa_y, options.FOVa_y, options.FOVa_y]);
            options.axial_fov = ([options.axialFOVOrig, (options.axial_fov - options.axialFOVOrig) /2, (options.axial_fov - options.axialFOVOrig) /2]);
            if size(options.x0, 1) == nx
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
            end
        end
        options.NxPrior = options.Nx(1);
        options.NyPrior = options.Ny(1);
        options.NzPrior = options.Nz(1);
    else
        if ~isfield(options, 'eFOVIndices') || numel(options.eFOVIndices) < 1
            options.eFOVIndices = zeros(options.Nz,1);
            options.eFOVIndices((options.Nz - options.NzOrig)/2 + 1 : end - (options.Nz - options.NzOrig)/2) = 1;
        end
        options.eFOVIndices = uint8(options.eFOVIndices);
        options.NzPrior = sum(options.eFOVIndices);
        options.maskPrior = zeros(options.Nx, options.Ny, 'uint8');
        options.maskPrior((options.Nx - options.NxOrig)/2 + 1 : end - (options.Nx - options.NxOrig)/2, (options.Ny - options.NyOrig)/2 + 1 : end - (options.Ny - options.NyOrig)/2) = uint8(1);
        options.NxPrior = sum(options.maskPrior(:,round(end/2)));
        options.NyPrior = sum(options.maskPrior(round(end/2),:));
        if options.useMaskBP
            options.maskPrior = options.maskPrior - uint8(~options.maskBP);
        end
    end
else
    options.Nx = uint32(options.Nx);
    options.Ny = uint32(options.Ny);
    options.Nz = uint32(options.Nz);
    options.NxFull = options.Nx;
    options.NyFull = options.Ny;
    options.NzFull = options.Nz;
end
if options.offsetCorrection
    %     distance = zeros(options.nRowsD, options.nProjections);
    options.OffsetLimit = zeros(options.nProjections, 1, options.cType);
    for kk = 1 : options.nProjections
        % kk = 1;
        sx = options.x(1,kk);
        sy = options.x(2,kk);
%         sx = [sx;sx];
%         sy = [sy;sy];
        dx = options.x(4,kk);
        dy = options.x(5,kk);
        ii = (0 : 0.25 : options.nRowsD)' - options.nRowsD / 2;
%         jj = .5;
        dx = dx + options.z(1,kk) * ii + options.z(4,kk) * ii;
        dy = dy + options.z(2,kk) * ii + options.z(5,kk) * ii;
%         dx = [dx + options.z(1,kk) * (options.nRowsD / 2) + options.z(4,kk) * (options.nRowsD / 2); dx - options.z(1,kk) * (options.nRowsD / 2) - options.z(4,kk) * (options.nRowsD / 2)];
%         dy = [dy + options.z(2,kk) * (options.nRowsD / 2) + options.z(5,kk) * (options.nRowsD / 2); dy - options.z(2,kk) * (options.nRowsD / 2) - options.z(5,kk) * (options.nRowsD / 2)];
        dist = abs((dx - sx) .* (sy) - ((sx) .* (dy - sy))) ./ sqrt((dx - sx).^2 + (dy - sy).^2);
%         options.OffsetLimit(kk) = min(dist) * 1.85;
        [~, ind] = min(dist);
        %     distance(:,kk) = dist;
        %     indeksit(kk) = ind;
        %     end
        % if ind > options.nRowsD / 2
        options.OffsetLimit(kk) = (ind) * options.dPitchY/4;
        % else
        %     options.OffsetLimit = (options.nRowsD - ind) * options.dPitchY;
        % end
    end
end
end