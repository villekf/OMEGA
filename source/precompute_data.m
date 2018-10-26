function precompute_data(options)
%% PRECOMPUTE PHASE FOR GATE DATA
% This function needs to be run when a specific machine is used for a first
% time
%
% This function also needs to be run if sinogram geometries or image sizes
% are later changed
%
% You also need to run this function if you later decide to use raw data, 
% precomputed LORs or OpenCL reconstruction and haven't run this function
% before
%
% Input the machine, image, sinogram and reconstruction parameters
%
% All output data are saved in MAT-files


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
    sinogram_coordinates(options);
    if options.verbose
        disp('Sinogram coordinates computed')
        toc
    end
    if options.precompute_lor || options.reconstruction_method == 2 || options.reconstruction_method == 4 || options.precompute_all
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