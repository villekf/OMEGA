function [varargout] = lor_pixel_count_prepass(options, z, x, detIndices)
%% Count the number of voxels each LOR traverses
% This function counts the number of voxels that each LOR traverses, i.e.
% the total number of voxels along each LOR. This is needed for 
% implementation 1.
% Separate codes for the sinogram and raw data.
%
% OUTPUTS:
%   lor = The number of voxels each LOR traverses (double precision)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019-2024 Ville-Veikko Wettenhovi
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

% if nargout > 4
%     error('Too many output arguments')
% end
listmode = false;
nMeas = options.size_x * options.size_y * options.nProjections;
folder = fileparts(which('lor_pixel_count_prepass.m'));
folder = [folder(1:end-(6 + 8)), 'mat-files/'];
folder = strrep(folder, '\','/');

% if exist('feature','builtin') == 5
%     nCores = uint32(feature('numcores'));
% else
%     nCores = uint32(1);
% end

%% Raw data
if options.use_raw_data
    if isempty(options.machine_name)
        file_string = [folder options.name '_detector_locations_' num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '_raw_' num2str(options.det_per_ring) 'x' num2str(options.rings) '.mat'];
    else
        file_string = [folder options.machine_name '_detector_locations_' num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '_raw_' num2str(options.det_per_ring) 'x' num2str(options.rings) '.mat'];
    end
else %% Sinogram data
    if isempty(options.machine_name)
        file_string = [folder options.name '_lor_pixel_count_' num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '_sino_' num2str(options.Ndist) 'x' ...
            num2str(options.Nang) 'x' num2str(options.TotSinos) '.mat'];
    else
        file_string = [folder options.machine_name '_lor_pixel_count_' num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '_sino_' num2str(options.Ndist) 'x' ...
            num2str(options.Nang) 'x' num2str(options.TotSinos) '.mat'];
    end
end
lor = projector_mex( options, options.Nx, options.Ny, options.Nz, options.dx, options.dy, options.dz, options.bx, options.by, options.bz, z, x, options.size_x, 0, 0, nMeas, false, false, false, 0, 0, uint32(0), uint16(0), detIndices, options.det_per_ring, ... % 25
    false, 0, 0, 0, options.verbose, 0, options.use_raw_data, uint32(3), listmode, 1, 0, nMeas, options.dPitchY, options.nProjections);

if exist('OCTAVE_VERSION', 'builtin') == 0
    save(file_string,'lor','-v7.3')
else
    save(file_string,'lor','-v7')
end
varargout{1} = lor;
end