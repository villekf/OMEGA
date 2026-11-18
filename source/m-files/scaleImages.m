function pz = scaleImages(pz, varargin)
%SCALEIMAGES This function scales the input images to Bq/mL
%   This function can be used to scale the images reconstructed by OMEGA
%   to Bq/mL units (by default, OMEGA outputs images with number of counts
%   / voxel units). The image is divided by the voxel volume in mL (the
%   voxel size is in mm, i.e. 1 mL = 1000 mm^3) and by the duration of the
%   time window the counts were collected in (s).
%
%   Static data: the duration is the length of the time window,
%   min(options.end, options.tot_time) - options.start, which is the total
%   time of the measurement when the whole measurement is reconstructed.
%   The duration has to be finite. If you have used infinite total time,
%   input the correct total time in the second input.
%
%   Dynamic data: each time step is divided by its own duration. The
%   durations are taken, in order of priority, from the optional input
%   frameTimes, from options.partitions when it is a vector of time-step
%   durations (as it is after the data has been loaded) or, when
%   options.partitions is the number of time steps, by dividing the time
%   window equally. The time dimension of the images is the 5th (X x Y x Z
%   x iterations x time) or, when the iteration dimension has been
%   removed, the 4th dimension (X x Y x Z x time).
%
% EXAMPLES:
%   pz = scaleImages(pz)
%   pz = scaleImages(pz, options)
%   pz = scaleImages(pz, options, frameTimes)
%   pz = scaleImages(pz, time, dx, dy, dz)
%   pz = scaleImages(pz, time, dx, dy, dz, frameTimes)
%
% INPUTS:
%   pz = Input image(s), can be the cell array pz created by the main files
%   or any single ND image.
%   options = The options struct, use this if you are inputting only a
%   single image (i.e. not the cell array). Total time, start and end
%   times, time steps, FOV info and image sizes are taken from this.
%   time = Total time of the measurement (s), required only if using
%   infinite total time. For dynamic data the time window is divided
%   between the time steps (unless frameTimes or a vector of durations in
%   options.partitions is available).
%   dx = Distance between adjacent voxels in x-direction (mm) (optional)
%   dy = Distance between adjacent voxels in y-direction (mm) (optional)
%   dz = Distance between adjacent voxels in z-direction (mm) (optional)
%   frameTimes = Duration of each time step (s), a vector with one value
%   per time step (optional). Overrides the durations from the options.
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

totTime = [];
startTime = 0;
endTime = Inf;
nSteps = [];
stepTimes = [];
FOVx = [];
FOVy = [];
FOVz = [];
Nx = [];
Ny = [];
Nz = [];
dx = [];
dy = [];
dz = [];
frameTimes = [];
if iscell(pz)
    prop = pz{end};
    totTime = getValue(prop, {'totalTime'});
    startTime = getValue(prop, {'startTime'}, 0);
    endTime = getValue(prop, {'endTime'}, Inf);
    steps = getValue(prop, {'nTimeSteps'});
    if numel(steps) > 1
        stepTimes = steps;
    else
        nSteps = steps;
    end
    FOVx = prop.FOVx;
    FOVy = prop.FOVy;
    FOVz = prop.axialFOV;
    Nx = prop.Nx;
    Ny = prop.Ny;
    Nz = prop.Nz;
end
args = varargin;
if ~isempty(args)
    if isstruct(args{1})
        options = args{1};
        totTime = getValue(options, {'tot_time'});
        startTime = getValue(options, {'start'}, 0);
        endTime = getValue(options, {'end'}, Inf);
        steps = getValue(options, {'partitions'});
        stepTimes = [];
        nSteps = [];
        if numel(steps) > 1
            stepTimes = steps;
        else
            nSteps = steps;
        end
        FOVx = options.FOVa_x;
        FOVy = options.FOVa_y;
        FOVz = options.axial_fov;
        Nx = options.Nx;
        Ny = options.Ny;
        Nz = options.Nz;
    else
        totTime = args{1};
        % A total time given as input is the length of the time window
        startTime = 0;
        endTime = Inf;
    end
    if numel(args) == 2
        frameTimes = args{2};
    elseif numel(args) == 3
        error('Input all of dx, dy and dz: scaleImages(pz, time, dx, dy, dz)');
    elseif numel(args) >= 4
        dx = args{2};
        dy = args{3};
        dz = args{4};
        if numel(args) >= 5
            frameTimes = args{5};
        end
    end
end
if isempty(dx) && (isempty(FOVx) || isempty(Nx) || isempty(FOVy) || isempty(FOVz))
    error('Voxel size not found, input it with pz = scaleImages(pz, time, dx, dy, dz) or give the options');
end

% Duration(s) of the time step(s)
if ~isempty(frameTimes)
    durations = double(frameTimes(:));
elseif ~isempty(stepTimes)
    durations = double(stepTimes(:));
else
    if isempty(totTime)
        error('Total time not found, input it with pz = scaleImages(pz, total_time) or give the options');
    end
    if isempty(nSteps) || nSteps < 1
        nSteps = 1;
    end
    windowEnd = min(endTime, totTime);
    if isinf(windowEnd)
        error('Total time is infinite, input finite total time with pz = scaleImages(pz, total_time)');
    end
    durations = repmat((windowEnd - startTime) / double(nSteps), double(nSteps), 1);
end
if any(~isfinite(durations)) || any(durations <= 0)
    error('The durations of the time steps have to be positive and finite');
end
if isscalar(durations)
    durations = durations(1);
end

if iscell(pz)
    for kk = 1 : size(pz,1) - 1
        if isempty(pz{kk})
            continue
        else
            pz{kk} = scaleOne(pz{kk}, voxelVolume(kk, dx, dy, dz, FOVx, FOVy, FOVz, Nx, Ny, Nz), durations);
        end
    end
else
    pz = scaleOne(pz, voxelVolume(1, dx, dy, dz, FOVx, FOVy, FOVz, Nx, Ny, Nz), durations);
end
end

function val = getValue(s, names, default)
%GETVALUE Returns the value of the first existing non-empty field of s listed in names
if nargin < 3
    default = [];
end
val = default;
for kk = 1 : numel(names)
    if isfield(s, names{kk}) && ~isempty(s.(names{kk}))
        val = double(s.(names{kk}));
        return
    end
end
end

function V = voxelVolume(idx, dx, dy, dz, FOVx, FOVy, FOVz, Nx, Ny, Nz)
%VOXELVOLUME Voxel volume (mm^3) of the image idx (multi-resolution images can have different voxel sizes)
if isempty(dx)
    dx = pick(FOVx, idx) / pick(Nx, idx);
    dy = pick(FOVy, idx) / pick(Ny, idx);
    dz = pick(FOVz, idx) / pick(Nz, idx);
end
V = dx(1) * dy(1) * dz(1);
end

function v = pick(x, idx)
%PICK Returns the element idx of x or the first element if x has fewer elements
if numel(x) >= idx
    v = double(x(idx));
else
    v = double(x(1));
end
end

function img = scaleOne(img, V, durations)
%SCALEONE Scales the image to Bq/mL, V is the voxel volume in mm^3 and durations the duration(s) of the time step(s) in s
if isscalar(durations)
    img = img * (1000 / (V * durations));
    return
end
Nt = numel(durations);
nd = ndims(img);
if nd >= 5 && size(img, 5) == Nt
    timeDim = 5;
elseif nd == 4 && size(img, 4) == Nt
    timeDim = 4;
else
    error(['The number of time steps (' num2str(Nt) ' durations) does not match the image size (' num2str(size(img)) '). ' ...
        'The time dimension has to be the 5th (X x Y x Z x iterations x time) or the 4th (X x Y x Z x time) dimension.']);
end
shape = ones(1, timeDim);
shape(timeDim) = Nt;
img = img .* cast(reshape(1000 ./ (V * durations), shape), 'like', img);
end
