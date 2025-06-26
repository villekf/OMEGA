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

%                       _________________
%                       |               |           
%                       |               |           
%                       |      4        |           
%                       |               |           
%            ___________|_______________|___________
%           |           |               |           |
%           |           |               |           |
%           |           |               |           |
%           |     5     |      1/2      |     6     |
%           |           |               |           |
%           |           |               |           |
%           |___________|_______________|___________|
%                       |               |           
%                       |               |           
%                       |      3        |           
%                       |               |           
%                       |_______________|           


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

                NxM = round(options.NxOrig * options.multiResolutionScale);
                dxM = options.FOVxOrig / NxM;
                NyM = round(options.NyOrig * options.multiResolutionScale);
                dyM = options.FOVyOrig / NyM;
                NzM = round(options.NzOrig * options.multiResolutionScale);
                dzM = options.axialFOVOrig / NzM;

                NxM2 = round((options.FOVa_x - options.FOVxOrig) / 2 / dxM) * 2;
                FOVxM = NxM2 * dxM;
                NyM2 = round((options.FOVa_y - options.FOVyOrig) / 2 / dyM) * 2;
                FOVyM = NyM2 * dyM;
                NzM2 = round((options.axial_fov - options.axialFOVOrig) / 2 / dzM) * 2;
                FOVzM = NzM2 * dzM;

                options.FOVa_x = ([options.FOVxOrig, options.FOVxOrig, options.FOVxOrig, ...
                    FOVxM / 2, FOVxM / 2, ...
                    options.FOVxOrig, options.FOVxOrig]);

                options.FOVa_y = ([options.FOVyOrig, options.FOVyOrig, options.FOVyOrig, ...
                    FOVyM + options.FOVyOrig, FOVyM + options.FOVyOrig, ...
                    FOVyM / 2, FOVyM / 2]);
                % options.FOVa_y = ([options.FOVyOrig, options.FOVyOrig, options.FOVyOrig, ...
                %     options.FOVyOrig, options.FOVyOrig, ...
                %     FOVyM / 2, FOVyM / 2]);

                options.axial_fov = ([options.axialFOVOrig, FOVzM / 2, FOVzM / 2, ...
                    FOVzM + options.axialFOVOrig, FOVzM + options.axialFOVOrig,...
                    FOVzM + options.axialFOVOrig, FOVzM + options.axialFOVOrig]);
                % options.axial_fov = ([options.axialFOVOrig, FOVzM / 2, FOVzM / 2, ...
                %     options.axialFOVOrig, options.axialFOVOrig,...
                %     options.axialFOVOrig, options.axialFOVOrig]);

                options.Nx = uint32([options.NxOrig, NxM, NxM, ...
                    NxM2 / 2, NxM2 / 2, ...
                    NxM, NxM]);

                options.Ny = uint32([options.NyOrig, NyM, NyM, ...
                    NyM + NyM2, NyM + NyM2, ...
                    NyM2 / 2, NyM2 / 2]);
                % options.Ny = uint32([options.NyOrig, NyM, NyM, ...
                %     NyM, NyM, ...
                %     NyM2 / 2, NyM2 / 2]);

                options.Nz = uint32([options.NzOrig, NzM2 / 2, NzM2 / 2, ...
                    NzM + NzM2, NzM + NzM2, ...
                    NzM + NzM2, NzM + NzM2]);
                % options.Nz = uint32([options.NzOrig, NzM2 / 2, NzM2 / 2, ...
                %     NzM, NzM, ...
                %     NzM, NzM]);

                disp(['Extended FOV is ' num2str((FOVxM + options.FOVxOrig)/options.FOVxOrig * 100) ' % of the original'])
            % If the initial value only covers the original, dense, volume,
            % it needs to be resized to match the extended FOV with
            % multi-resolution
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
            elseif size(options.x0, 1) == options.NxOrig
                options.x1 = ones(options.Nx(2), options.Ny(2), options.Nz(2)) * 1e-5;
                options.x2 = ones(options.Nx(3), options.Ny(3), options.Nz(3)) * 1e-5;
                options.x3 = ones(options.Nx(4), options.Ny(4), options.Nz(4)) * 1e-5;
                options.x4 = ones(options.Nx(5), options.Ny(5), options.Nz(5)) * 1e-5;
                options.x5 = ones(options.Nx(6), options.Ny(6), options.Nz(6)) * 1e-5;
                options.x6 = ones(options.Nx(7), options.Ny(7), options.Nz(7)) * 1e-5;
                if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
                    options.x0 = [options.x0(:);options.x1(:); options.x2(:); options.x3(:); options.x4(:); options.x5(:); options.x6(:)];
                    options = rmfield(options, {'x1','x2','x3','x4','x5','x6'});
                end
            end
        % The above case handles a situation where both the transaxial and
        % axial FOVs are extended. Here only the transaxial FOV is extended
        elseif options.transaxialEFOV && ~options.axialEFOV
            options.nMultiVolumes = 4;
                NxM = round(options.NxOrig * options.multiResolutionScale);
                dxM = options.FOVxOrig / NxM;
                NyM = round(options.NyOrig * options.multiResolutionScale);
                dyM = options.FOVyOrig / NyM;
                NzM = round(options.NzOrig * options.multiResolutionScale);

                NxM2 = round((options.FOVa_x - options.FOVxOrig) / 2 / dxM) * 2;
                FOVxM = NxM2 * dxM;
                NyM2 = round((options.FOVa_y - options.FOVyOrig) / 2 / dyM) * 2;
                FOVyM = NyM2 * dyM;

                options.FOVa_x = ([options.FOVxOrig, ...
                    FOVxM / 2, FOVxM / 2, ...
                    options.FOVxOrig, options.FOVxOrig]);

                options.FOVa_y = ([options.FOVyOrig, ...
                    FOVyM + options.FOVyOrig, FOVyM + options.FOVyOrig, ...
                    FOVyM / 2, FOVyM / 2]);

                options.axial_fov = ([options.axialFOVOrig, ...
                    options.axialFOVOrig, options.axialFOVOrig,...
                    options.axialFOVOrig, options.axialFOVOrig]);

                options.Nx = uint32([options.NxOrig, ...
                    NyM2 / 2, NyM2 / 2, ...
                    NxM, NxM]);

                options.Ny = uint32([options.NyOrig, ...
                    NxM + NyM2, NxM + NyM2, ...
                    NyM2 / 2, NyM2 / 2]);

                options.Nz = uint32([options.NzOrig, NzM, NzM, ...
                    NzM, NzM]);

                disp(['Extended FOV is ' num2str((FOVxM + options.FOVxOrig)/options.FOVxOrig * 100) ' % of the original'])
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
            elseif size(options.x0, 1) == options.NxOrig
                options.x1 = ones(options.Nx(2), options.Ny(2), options.Nz(2)) * 1e-5;
                options.x2 = ones(options.Nx(3), options.Ny(3), options.Nz(3)) * 1e-5;
                options.x3 = ones(options.Nx(4), options.Ny(4), options.Nz(4)) * 1e-5;
                options.x4 = ones(options.Nx(5), options.Ny(5), options.Nz(5)) * 1e-5;
                if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
                    options.x0 = [options.x0(:);options.x1(:); options.x2(:); options.x3(:); options.x4(:); options.x5(:); options.x6(:)];
                    options = rmfield(options, {'x1','x2','x3','x4','x5','x6'});
                end
            end
        % Only axial FOV is extended
        else
            options.nMultiVolumes = 2;
                NxM = round(options.NxOrig * options.multiResolutionScale);
                NyM = round(options.NyOrig * options.multiResolutionScale);
                NzM = round(options.NzOrig * options.multiResolutionScale);
                dzM = options.axialFOVOrig / NzM;

                NzM2 = round((options.axial_fov - options.axialFOVOrig) / 2 / dzM) * 2;
                FOVzM = NzM2 * dzM;

                options.FOVa_x = ([options.FOVxOrig, options.FOVxOrig, options.FOVxOrig]);

                options.FOVa_y = ([options.FOVyOrig, options.FOVyOrig, options.FOVyOrig]);

                options.axial_fov = ([options.axialFOVOrig, FOVzM / 2, FOVzM / 2]);

                options.Nx = uint32([options.NxOrig, NxM, NxM]);

                options.Ny = uint32([options.NyOrig, NyM, NyM]);

                options.Nz = uint32([options.NzOrig, NzM2 / 2, NzM2 / 2]);

                disp(['Extended FOV is ' num2str((FOVzM + options.axialFOVOrig)/options.axialFOVOrig * 100) ' % of the original'])
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
            elseif size(options.x0, 1) == options.NxOrig
                options.x1 = ones(options.Nx(2), options.Ny(2), options.Nz(2)) * 1e-5;
                options.x2 = ones(options.Nx(3), options.Ny(3), options.Nz(3)) * 1e-5;
                if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
                    options.x0 = [options.x0(:);options.x1(:); options.x2(:); options.x3(:); options.x4(:); options.x5(:); options.x6(:)];
                    options = rmfield(options, {'x1','x2','x3','x4','x5','x6'});
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
        sz = options.x(3,kk);
        s = [sx;sy;sz];
%         sx = [sx;sx];
%         sy = [sy;sy];
        dx = options.x(4,kk);
        dy = options.x(5,kk);
        dz = options.x(6,kk);
        % p = [options.oOffsetX;options.oOffsetY;options.oOffsetZ];
        d = [dx;dy;dz];
        % dist = norm(cross((d -s),(s - p))) / norm(d - s);
        ii = (0 : 0.25 : options.nRowsD)' - options.nRowsD / 2;
%         jj = .5;
        if options.pitch
        dx = dx + options.z(1,kk) * ii + options.z(4,kk) * ii;
        dy = dy + options.z(2,kk) * ii + options.z(5,kk) * ii;
        else
        dx = dx + options.z(1,kk) * ii;
        dy = dy + options.z(2,kk) * ii;
        end
        % dz = dz + options.z(3,kk) * ii + options.z(6,kk) * ii;
        % apuX = [options.z(1,kk);options.z(2,kk);options.z(3,kk)] * options.nRowsD / 2;
        % apuY = [options.z(4,kk);options.z(5,kk);options.z(6,kk)] * options.nColsD / 2;
        % apuX = [options.z(1,kk);options.z(2,kk);options.z(3,kk)];
        % apuY = [options.z(4,kk);options.z(5,kk);options.z(6,kk)];

        % d1 = [dx(10);dy(10);dz(10)];
        % d2 = apuX - apuY;
        % d3 = d - apuX - apuY;
        % n = cross(d2 - d, d3 - d);
        % n = cross(apuX, apuY);
        % n = n ./ norm(n);
        % v = s - d;
        % point = s - dot(v, n) .* n;
        % boundary = [dx(end);dy(end);dz];
        % dist = norm(point - boundary);
        % d = [dx;dy;dz];
        % n = abs(s - d) / norm(s-d);
        % p = [options.oOffsetX;options.oOffsetY;options.oOffsetZ];
%         dx = [dx + options.z(1,kk) * (options.nRowsD / 2) + options.z(4,kk) * (options.nRowsD / 2); dx - options.z(1,kk) * (options.nRowsD / 2) - options.z(4,kk) * (options.nRowsD / 2)];
%         dy = [dy + options.z(2,kk) * (options.nRowsD / 2) + options.z(5,kk) * (options.nRowsD / 2); dy - options.z(2,kk) * (options.nRowsD / 2) - options.z(5,kk) * (options.nRowsD / 2)];
        dist = abs((dx - sx) .* (sy) - ((sx) .* (dy - sy))) ./ sqrt((dx - sx).^2 + (dy - sy).^2);
%         % dist = norm((s-p) - (dot((s-p),n) .* n));
% %         options.OffsetLimit(kk) = min(dist) * 1.85;
        [~, ind] = min(dist);
%         %     distance(:,kk) = dist;
%         %     indeksit(kk) = ind;
%         %     end
%         if ind > options.nRowsD / 2 * 4
%             ind = options.nRowsD * 4 - ind;
%         end
        options.OffsetLimit(kk) = (ind + 0) * options.dPitchY/4;
        % options.OffsetLimit(kk) = dist;
        % options.OffsetLimit(kk) = (options.nRowsD * 4 - ind) * options.dPitchY/4;
        % else
        %     options.OffsetLimit(kk) = (options.nRowsD - ind) * options.dPitchY;
        % end
    end
end
end