function [xx,yy,zz,dx,dy,dz,bx,by,bz] = computePixelSize(FOV, N, offset, useMultiResolution, multiResolutionShift, cType)
%COMPUTEPIXELSIZE Computes the pixel size and distance from origin
%   Utility function for OMEGA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021-2025 Ville-Veikko Wettenhovi
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

% Pixel boundaries
etaisyys = -(FOV) / 2;
dx = zeros(1,size(FOV,2));
dy = zeros(1,size(FOV,2));
dz = zeros(1,size(FOV,2));
bx = zeros(1,size(FOV,2));
by = zeros(1,size(FOV,2));
bz = zeros(1,size(FOV,2));
% Derive spacing directly from the requested field of view and voxel count.
% Taking the first difference of a single-precision linspace can give
% different spacings for two slabs that are intended to share the same
% lattice (for example, one 64-slice grid and two 32-slice grids).
voxelSpacing = double(FOV) ./ double(N);
if useMultiResolution && size(FOV, 2) > 1
    fineSpacing = voxelSpacing(:, 1);
    tolerance = max(1e-8, abs(fineSpacing) .* 1e-5);
    sameResolution = abs(voxelSpacing - fineSpacing) <= tolerance;
    matchingSpacing = repmat(fineSpacing, 1, size(FOV, 2));
    voxelSpacing(sameResolution) = matchingSpacing(sameResolution);
end
for kk = size(FOV,2) : - 1 : 1
    % We compute the pixel boundaries for each volume if a multi-resolution
    % volume, otherwise for the main volume only
    % xx(1), yy(1) and zz(1) is the boundary coordinate of the first voxel
    % These are the boundary coordinates, not center
    xx = double(linspace(etaisyys(1,kk) + offset(1), -etaisyys(1,kk) + offset(1), N(1,kk) + 1));
    yy = double(linspace(etaisyys(2,kk) + offset(2), -etaisyys(2,kk) + offset(2), N(2,kk) + 1));
    zz = double(linspace(etaisyys(3,kk) + offset(3), -etaisyys(3,kk) + offset(3), N(3,kk) + 1));

    % Distance of adjacent voxels, i.e. voxel size
    dx(kk) = voxelSpacing(1, kk);
    dy(kk) = voxelSpacing(2, kk);
    dz(kk) = voxelSpacing(3, kk);

    % Distance of image from the origin
    % Different cases for different volumes
    if kk == 1
        bx(kk) = xx(1,1);
        by(kk) = yy(1,1);
        bz(kk) = zz(1,1);
        if ~useMultiResolution % eFOV uses only main volume, shift
            bx(kk) = bx(kk) + multiResolutionShift(1);
            by(kk) = by(kk) + multiResolutionShift(2);
            bz(kk) = bz(kk) + multiResolutionShift(3);
        end
    else
        % Side volumes
        if kk > 5 || (size(FOV,2) == 5 && kk > 3)
            if mod(kk,2) == 1
                by(kk) = offset(2) + FOV(2,1) / 2;
            else
                by(kk) = offset(2) - FOV(2,1) / 2 - FOV(2,kk);
            end
            bx(kk) = xx(1,1);
            bz(kk) = zz(1,1) + multiResolutionShift(3);
        % Top and bottom volumes
        elseif (kk > 3 && kk < 6) || (size(FOV,2) == 5 && kk > 1)
            if mod(kk,2) == 1
                bx(kk) = etaisyys(1,1) + offset(1) + FOV(1,1);
            else
                bx(kk) = etaisyys(1,1) + offset(1) - FOV(1,kk);
            end
            by(kk) = yy(1,1) + multiResolutionShift(2); % Has to be shifted as side volumes are attached to main volume X and Y; top and bottom volumes are attached only by X
            bz(kk) = zz(1,1) + multiResolutionShift(3);
        % Front and back volumes
        elseif kk > 1 && kk < 4
            bx(kk) = xx(1,1);
            by(kk) = yy(1,1);
            if mod(kk,2) == 1
                bz(kk) = etaisyys(3,1) + offset(3) + FOV(3,1);
            else
                bz(kk) = etaisyys(3,1) + offset(3) - FOV(3,kk);
            end
        end
    end
end

xx = cast(xx, cType);
yy = cast(yy, cType);
zz = cast(zz, cType);
dx = cast(dx, cType);
dy = cast(dy, cType);
dz = cast(dz, cType);
bx = cast(bx, cType);
by = cast(by, cType);
bz = cast(bz, cType);

% Set neighboring multi-resolution slab origins from the central volume's
% actual voxel boundaries. This makes the shared face identical after the
% geometry is cast to the projector precision, even when FOV/N arithmetic
% has tiny rounding differences.
if useMultiResolution
    nVolumes = size(FOV, 2);
    if nVolumes == 3
        bz(2) = slabOrigin(bz(1), N(3, 2), dz(2), -1, cType);
        bz(3) = slabOrigin(bz(1), N(3, 1), dz(1), 1, cType);
    elseif nVolumes == 5
        bx(2) = slabOrigin(bx(1), N(1, 2), dx(2), -1, cType);
        bx(3) = slabOrigin(bx(1), N(1, 1), dx(1), 1, cType);
        by(4) = slabOrigin(by(1), N(2, 4), dy(4), -1, cType);
        by(5) = slabOrigin(by(1), N(2, 1), dy(1), 1, cType);
    elseif nVolumes == 7
        bz(2) = slabOrigin(bz(1), N(3, 2), dz(2), -1, cType);
        bz(3) = slabOrigin(bz(1), N(3, 1), dz(1), 1, cType);
        bx(4) = slabOrigin(bx(1), N(1, 4), dx(4), -1, cType);
        bx(5) = slabOrigin(bx(1), N(1, 1), dx(1), 1, cType);
        by(6) = slabOrigin(by(1), N(2, 6), dy(6), -1, cType);
        by(7) = slabOrigin(by(1), N(2, 1), dy(1), 1, cType);
    end
end
end

function origin = slabOrigin(mainOrigin, voxelCount, spacing, direction, cType)
% Match the single-precision arithmetic used to form each projector bmax.
extent = cast(cast(voxelCount, cType) .* spacing, cType);
origin = cast(mainOrigin + cast(direction, cType) .* extent, cType);
end
