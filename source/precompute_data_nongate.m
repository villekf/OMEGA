function precompute_data_nongate(options)
%% PRECOMPUTE PHASE FOR PET SINOGRAM DATA
% This function needs to be run when a specific machine is used for a first
% time (non-GATE data)
%
% This function also needs to be run if sinogram geometries or image sizes
% are later changed
%
% You also need to run this function if you later decide to use precomputed
% LORs or OpenCL reconstruction and haven't run this function before
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


if options.verbose
    tic
end
detector_coordinates(options);
if options.verbose
    disp('Detector coordinates computed')
    toc
end
if options.use_raw_data || options.precompute_all
    if options.verbose
        tic
    end
    detector_pair_formation_raw(options);
    if options.verbose
        disp('Detector pairs formed for raw data')
        toc
    end
end


% Precompute sinogram coordinates for reconstruction
if options.use_raw_data == false || options.precompute_all
    if options.verbose
        tic
    end
    sinogram_coordinates_nongate(options);
    if options.verbose
        disp('Sinogram coordinates computed')
        toc
    end
    if options.precompute_lor || options.reconstruction_method == 2 || options.precompute_all || options.reconstruction_method == 4
        if options.verbose
            tic
        end
        lor_pixel_count_prepass(options);
        if options.verbose
            disp('LOR prepass phase complete')
            toc
        end
    end
end
