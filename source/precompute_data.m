function precompute_data(options)
%% PRECOMPUTE PHASE FOR GATE DATA
% This function should be run when a specific machine is used for a first
% time and the precompute_lor is set to true. However, if the required
% files are not found, then the corresponding functions will be
% automatically run.
%
% This function should also be run if sinogram geometries or image sizes
% are later changed or if you later decide to use raw data, 
% precomputed LORs or OpenCL reconstruction and haven't run this function
% before. However, as mentioned above, this function is automatically run
% if the necessary data are not found.
%
% Input the machine, image, sinogram and reconstruction parameters
%
% All output data are saved in MAT-files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi
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


if options.precompute_lor || options.precompute_all
    if options.verbose
        disp('Beginning LOR precomputation phase. This may take several minutes.')
        tic
    end
    % Compute the number of voxels each LOR traverses.
    % This step will speed up the computations as it will allow to skip
    % beforehand all the LORs that do not intersect with the FOV. It will
    % also enhance the performance of matrix versions of the system matrix
    % formation due to preallocated memory.
    lor_pixel_count_prepass(options);
    if options.verbose
        disp('LOR precomputation phase complete')
        toc
    end
end