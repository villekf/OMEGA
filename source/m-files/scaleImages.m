function pz = scaleImages(pz, varargin)
%SCALEIMAGES This function scales the input images to Bq/mL
%   This function can be used to scale the images reconstructed by OMEGA
%   to Bq/mL units (by default, OMEGA outputs images with number of counts
%   / voxel units). This function requires the total time of the
%   measurement and also requires it to be finite. If you have used
%   infinite time, you can input the correct total time in the second
%   input. If you used finite total time when reconstructing the data, no
%   other inputs other than the cell array pz is needed.
%
% EXAMPLES:
%   pz = scaleImages(pz)
%   pz = scaleImages(pz, options)
%   pz = scaleImages(pz, time, dx, dy, dz)
%
% INPUTS:
%   pz = Input image(s), can be the cell array pz created by the main files
%   or any single ND image.
%   options = The options struct, use this if you are inputting only a
%   single image (i.e. not the cell array). Total time, FOV info and image
%   sizes are taken from this.
%   time = Total time of the measurement, required only if using infinite
%   total time. For dynamic cases 
%   dx = Distance between adjacent voxels in x-direction (mm) (optional)
%   dy = Distance between adjacent voxels in y-direction (mm) (optional)
%   dz = Distance between adjacent voxels in z-direction (mm) (optional)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ville-Veikko Wettenhovi
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

if iscell(pz)
    totTime = pz{end}.total_time;
    FOVx = pz{end}.FOV_x;
    FOVy = pz{end}.FOV_y;
    FOVz = pz{end}.axial_FOV;
    Nx = pz{end}.Nx;
    Ny = pz{end}.Ny;
    Nz = pz{end}.Nz;
end
if nargin >= 2
    if isstruct(varargin{1})
        totTime = options.tot_time;
        FOVx = options.FOVa_x;
        FOVy = options.FOVa_y;
        FOVz = options.axial_fov;
        Nx = options.Nx;
        Ny = options.Ny;
        Nz = options.Nz;
    else
        totTime = varargin{1};
    end
end
if isinf(totTime)
    error('Total time is infinite, input finite total time with pz = scaleImages(pz, total_time)');
end
if nargin >= 3
    dx = varargin{2};
    dy = varargin{3};
    dz = varargin{4};
else
    dx = FOVx / Nx;
    dy = FOVy / Ny;
    dz = FOVz / Nz;
end
V = dx * dy * dz;
scale = V * totTime;
if iscell(pz)
    for kk = 1 : size(pz,1) - 1
        if isempty(pz{kk})
            continue
        else
            pz{kk} = pz{kk} / scale;
        end
    end
else
    pz = pz / scale;
end
end

