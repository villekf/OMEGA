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
        % All the side volumes use exactly the coarse voxel size of the main volume (dxM, dyM, dzM), so that they all have the same
        % voxel size and lattice. The length of a side volume is therefore rounded UP to a whole number of coarse voxels: the
        % extended FOV never shrinks and grows by less than one coarse voxel on each side, the side volumes stay adjacent to the main
        % volume and only their outer faces move. The tolerance keeps a length that is a whole number of voxels (up to rounding error)
        % from getting one voxel more
        sideVoxels = @(len, d) ceil(double(len) / double(d) - 1e-6);

        % If the axial (or transaxial) extended FOV is smaller than one
        % multi-resolution voxel, the corresponding side volume(s) would
        % end up with zero thickness (NzM2/NxM2/NyM2 = 0), which crashes
        % later on. In that case disable the EFOV direction in question and
        % fall back to the other direction's multi-resolution branch (or to
        % no multi-resolution volumes at all if neither direction remains).
        if options.axialEFOV
            NzMchk = round(options.NzOrig * options.multiResolutionScale);
            dzMchk = options.axialFOVOrig / NzMchk;
            NzM2chk = sideVoxels((options.axial_fov - options.axialFOVOrig) / 2, dzMchk) * 2;
            if NzM2chk <= 0
                warning('Axial extended FOV is smaller than one multi-resolution voxel, disabling axial EFOV.')
                options.axialEFOV = false;
                options.Nz = options.NzOrig;
                options.axial_fov = options.axialFOVOrig;
            end
        end
        if options.transaxialEFOV
            NxMchk = round(options.NxOrig * options.multiResolutionScale);
            dxMchk = options.FOVxOrig / NxMchk;
            NxM2chk = sideVoxels((options.FOVa_x - options.FOVxOrig) / 2, dxMchk) * 2;
            NyMchk = round(options.NyOrig * options.multiResolutionScale);
            dyMchk = options.FOVyOrig / NyMchk;
            NyM2chk = sideVoxels((options.FOVa_y - options.FOVyOrig) / 2, dyMchk) * 2;
            if NxM2chk <= 0 || NyM2chk <= 0
                warning('Transaxial extended FOV is smaller than one multi-resolution voxel, disabling transaxial EFOV.')
                options.transaxialEFOV = false;
                options.Nx = options.NxOrig;
                options.Ny = options.NyOrig;
                options.FOVa_x = options.FOVxOrig;
                options.FOVa_y = options.FOVyOrig;
            end
        end
        if ~options.axialEFOV && ~options.transaxialEFOV
            % Neither direction has a usable extension left: skip
            % multi-resolution volume creation entirely.
            options.useMultiResolutionVolumes = false;
            options.nMultiVolumes = 0;
            options.NxPrior = options.Nx;
            options.NyPrior = options.Ny;
            options.NzPrior = options.Nz;
        else
        if options.axialEFOV && options.transaxialEFOV
            options.nMultiVolumes = 6;

            NxM = round(options.NxOrig * options.multiResolutionScale);
            dxM = options.FOVxOrig / NxM;
            NyM = round(options.NyOrig * options.multiResolutionScale);
            dyM = options.FOVyOrig / NyM;
            NzM = round(options.NzOrig * options.multiResolutionScale);
            dzM = options.axialFOVOrig / NzM;

            NxM2 = sideVoxels((options.FOVa_x - options.FOVxOrig) / 2, dxM) * 2;
            FOVxM = NxM2 * dxM;
            NyM2 = sideVoxels((options.FOVa_y - options.FOVyOrig) / 2, dyM) * 2;
            FOVyM = NyM2 * dyM;
            NzM2 = sideVoxels((options.axial_fov - options.axialFOVOrig) / 2, dzM) * 2;
            FOVzM = NzM2 * dzM;

            options = setMultiResPadding(options, [(options.FOVa_x - options.FOVxOrig) / 2, (options.FOVa_y - options.FOVyOrig) / 2, ...
                (options.axial_fov - options.axialFOVOrig) / 2], [FOVxM, FOVyM, FOVzM] / 2, [dxM, dyM, dzM]);

            options.FOVa_x = ([options.FOVxOrig, options.FOVxOrig, options.FOVxOrig, ...
                FOVxM / 2, FOVxM / 2, ...
                options.FOVxOrig, options.FOVxOrig]);

            options.FOVa_y = ([options.FOVyOrig, options.FOVyOrig, options.FOVyOrig, ...
                FOVyM + options.FOVyOrig, FOVyM + options.FOVyOrig, ...
                FOVyM / 2, FOVyM / 2]);

            options.axial_fov = ([options.axialFOVOrig, FOVzM / 2, FOVzM / 2, ...
                FOVzM + options.axialFOVOrig, FOVzM + options.axialFOVOrig,...
                FOVzM + options.axialFOVOrig, FOVzM + options.axialFOVOrig]);

            options.Nx = uint32([options.NxOrig, NxM, NxM, ...
                NxM2 / 2, NxM2 / 2, ...
                NxM, NxM]);

            options.Ny = uint32([options.NyOrig, NyM, NyM, ...
                NyM + NyM2, NyM + NyM2, ...
                NyM2 / 2, NyM2 / 2]);

            options.Nz = uint32([options.NzOrig, NzM2 / 2, NzM2 / 2, ...
                NzM + NzM2, NzM + NzM2, ...
                NzM + NzM2, NzM + NzM2]);

            disp(['Extended FOV is ' num2str((FOVxM + options.FOVxOrig)/options.FOVxOrig * 100) ' % of the original'])
            % If the initial value only covers the original, dense, volume,
            % it needs to be resized to match the extended FOV with
            % multi-resolution
            if size(options.x0, 1) == nx && min(options.x0(:)) ~= max(options.x0(:))
                % The initial value is for the requested extended FOV (padded in multiResReshape with the enlargement of it by the side
                % volumes, up to one coarse voxel). The main volume is its central crop at the full resolution and the side volumes are
                % the regions of it resized to their coarse grids
                options.x0 = multiResReshape(options.x0, options);
                options = unpackMultiResX0(options);
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
            % The above case handles a situation where both the transaxial and
            % axial FOVs are extended. Here only the transaxial FOV is extended
        elseif options.transaxialEFOV && ~options.axialEFOV
            options.nMultiVolumes = 4;
            NxM = round(options.NxOrig * options.multiResolutionScale);
            dxM = options.FOVxOrig / NxM;
            NyM = round(options.NyOrig * options.multiResolutionScale);
            dyM = options.FOVyOrig / NyM;
            NzM = round(options.NzOrig * options.multiResolutionScale);

            NxM2 = sideVoxels((options.FOVa_x - options.FOVxOrig) / 2, dxM) * 2;
            FOVxM = NxM2 * dxM;
            NyM2 = sideVoxels((options.FOVa_y - options.FOVyOrig) / 2, dyM) * 2;
            FOVyM = NyM2 * dyM;

            dzM = options.axialFOVOrig / NzM;
            options = setMultiResPadding(options, [(options.FOVa_x - options.FOVxOrig) / 2, (options.FOVa_y - options.FOVyOrig) / 2, 0], ...
                [FOVxM, FOVyM, 0] / 2, [dxM, dyM, dzM]);

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
                NxM2 / 2, NxM2 / 2, ...
                NxM, NxM]);

            options.Ny = uint32([options.NyOrig, ...
                NyM + NyM2, NyM + NyM2, ...
                NyM2 / 2, NyM2 / 2]);

            options.Nz = uint32([options.NzOrig, NzM, NzM, ...
                NzM, NzM]);

            disp(['Extended FOV is ' num2str((FOVxM + options.FOVxOrig)/options.FOVxOrig * 100) ' % of the original'])
            if size(options.x0, 1) == nx && min(options.x0(:)) ~= max(options.x0(:))
                % The initial value is for the requested extended FOV (padded in multiResReshape with the enlargement of it by the side
                % volumes, up to one coarse voxel). The main volume is its central crop at the full resolution and the side volumes are
                % the regions of it resized to their coarse grids
                options.x0 = multiResReshape(options.x0, options);
                options = unpackMultiResX0(options);
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
            % Only axial FOV is extended
        else
            options.nMultiVolumes = 2;
            NxM = round(options.NxOrig * options.multiResolutionScale);
            NyM = round(options.NyOrig * options.multiResolutionScale);
            NzM = round(options.NzOrig * options.multiResolutionScale);
            dzM = options.axialFOVOrig / NzM;

            NzM2 = sideVoxels((options.axial_fov - options.axialFOVOrig) / 2, dzM) * 2;
            FOVzM = NzM2 * dzM;

            dxM = options.FOVxOrig / NxM;
            dyM = options.FOVyOrig / NyM;
            options = setMultiResPadding(options, [0, 0, (options.axial_fov - options.axialFOVOrig) / 2], [0, 0, FOVzM] / 2, [dxM, dyM, dzM]);

            options.FOVa_x = ([options.FOVxOrig, options.FOVxOrig, options.FOVxOrig]);

            options.FOVa_y = ([options.FOVyOrig, options.FOVyOrig, options.FOVyOrig]);

            options.axial_fov = ([options.axialFOVOrig, FOVzM / 2, FOVzM / 2]);

            options.Nx = uint32([options.NxOrig, NxM, NxM]);

            options.Ny = uint32([options.NyOrig, NyM, NyM]);

            options.Nz = uint32([options.NzOrig, NzM2 / 2, NzM2 / 2]);

            disp(['Extended FOV is ' num2str((FOVzM + options.axialFOVOrig)/options.axialFOVOrig * 100) ' % of the original'])
            if size(options.x0, 1) == nx && min(options.x0(:)) ~= max(options.x0(:))
                % The initial value is for the requested extended FOV (padded in multiResReshape with the enlargement of it by the side
                % volumes, up to one coarse voxel). The main volume is its central crop at the full resolution and the side volumes are
                % the regions of it resized to their coarse grids
                options.x0 = multiResReshape(options.x0, options);
                options = unpackMultiResX0(options);
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
        end
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
    nP = options.nProjections;
    s = double(options.x([1 2], 1:nP));
    d = double(options.x([4 5], 1:nP));
    u = double(options.z([1 2], 1:nP));
    % Position t (in detector pixels, from the detector center along the row direction u) at which the line from the source s to
    % d + t * u passes through the rotation axis (0,0). The central ray is thus computed exactly instead of searched on a grid.
    t = (s(2,:) .* d(1,:) - s(1,:) .* d(2,:)) ./ (s(1,:) .* u(2,:) - s(2,:) .* u(1,:));
    koko = double(options.nRowsD) * double(options.dPitchX);
    % Distance (mm) of the central ray from the row-0 edge of the detector
    T = (t(:) + double(options.nRowsD) / 2) * double(options.dPitchX);
    if any(T <= 0 | T >= koko)
        warning('Offset correction: the central ray is outside the detector for some projections, clamping the offset limit.')
    end
    options.OffsetLimit = cast(min(max(T, 1e-3 * double(options.dPitchX)), koko - 1e-3 * double(options.dPitchX)), options.cType);
end
end

function options = setMultiResPadding(options, extRequested, extUsed, dM)
% The side volumes cover the requested extended FOV and, since their lengths are rounded up to a whole number of coarse voxels
% (dM, in mm), up to one coarse voxel more on each outer face. extRequested and extUsed are the length (mm) of the extension of the
% requested FOV on one side (x, y, z, the extension is symmetric) and the length that the side volumes have. The whole-grid inputs
% (initial image, type 6 attenuation map) are given for the requested extended FOV, options.NxFull x NyFull x NzFull voxels of the
% main volume, and are padded with the enlargement, in voxels of the main volume (multiResPadLow/High), see padMultiResVolume and
% padType6Attenuation. The padding is exact when the requested FOV and the enlargement are whole numbers of voxels
    dxMain = [options.FOVxOrig / double(options.NxOrig), options.FOVyOrig / double(options.NyOrig), options.axialFOVOrig / double(options.NzOrig)];
    enlarge = double(extUsed) - double(extRequested);
    gridFull = double([options.NxFull, options.NyFull, options.NzFull]);
    FOVRequested = double([options.FOVxOrig, options.FOVyOrig, options.axialFOVOrig]) + 2 * double(extRequested);
    options.multiResPadLow = enlarge ./ dxMain;
    options.multiResPadHigh = enlarge ./ dxMain;
    options.multiResPadExact = all(abs(options.multiResPadLow - round(options.multiResPadLow)) < 1e-3) && ...
        all(abs(gridFull .* dxMain - FOVRequested) < 1e-3 * dxMain);
    options.multiResUnionDims = round(gridFull + options.multiResPadLow + options.multiResPadHigh);
    if options.verbose > 0 && any(enlarge > 1e-4 * dM)
        disp(['The extended FOV was enlarged to a whole number of voxels (' num2str(dM(1)) ' x ' num2str(dM(2)) ' x ' num2str(dM(3)) ' mm) in the side volumes. The enlargement (mm) on both sides of x is ' ...
            num2str(enlarge(1)) ', y ' num2str(enlarge(2)) ' and z ' num2str(enlarge(3))])
    end
end


function options = unpackMultiResX0(options)
% options.x0 is the main volume and all the side volumes as one column vector (volume 0 first, each in the linear order of
% Nx(k) x Ny(k) x Nz(k)). With the implementations 2, 3 and 5 it is used as it is, otherwise the main volume is options.x0 and the
% side volumes options.x1, options.x2, ...
if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
    return
end
vec = options.x0;
nVol = double(options.nMultiVolumes) + 1;
offset = 0;
for kk = 1 : nVol
    sz = double([options.Nx(kk), options.Ny(kk), options.Nz(kk)]);
    vol = reshape(vec(offset + 1 : offset + prod(sz)), sz);
    offset = offset + prod(sz);
    if kk == 1
        options.x0 = vol;
    else
        options.(['x' num2str(kk - 1)]) = vol;
    end
end
end

function outputVol = multiResReshape(inputVol, options)
    scale = double(options.multiResolutionScale);
    inputVol = padMultiResVolume(inputVol, options, false);
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
