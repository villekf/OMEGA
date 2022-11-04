function options = setMissingValues(options)
%SETMISSINGVALUES Sets default values for variables that are missing
%   Utility function

if ~isfield(options, 'verticalOffset')
    options.verticalOffset = 0;
end
if ~isfield(options, 'horizontalOffset')
    options.horizontalOffset = 0;
end
if ~isfield(options, 'subsets')
    options.subsets = 1;
end
if ~isfield(options, 'subset_type')
    options.subset_type = 9;
end
if ~isfield(options, 'useMaskFP')
    options.useMaskFP = false;
end
if ~isfield(options, 'useMaskBP')
    options.useMaskBP = false;
end
if ~isfield(options, 'dPitchX')
    options.dPitchX = options.dPitch;
end
if ~isfield(options, 'dPitchY')
    options.dPitchY = options.dPitch;
end
if ~isfield(options,'bedOffset')
    options.bedOffset = [];
end
if ~isfield(options,'uCenter')
    options.uCenter = [];
end
if ~isfield(options,'vCenter')
    options.vCenter = [];
end
if ~isfield(options,'nBed')
    options.nBed = 1;
end
if ~isfield(options,'flip_image')
    options.flip_image = false;
end
if ~isfield(options,'offangle')
    options.offangle = 0;
end
if ~isfield(options, 'oOffsetX')
    options.oOffsetX = 0;
end
if ~isfield(options, 'oOffsetY')
    options.oOffsetY = 0;
end
if ~isfield(options, 'oOffsetZ')
    options.oOffsetZ = 0;
end
if ~isfield(options, 'tube_width_z')
    options.tube_width_z = 0;
end
if ~isfield(options, 'tube_width_xy')
    options.tube_width_xy = 0;
end
if ~isfield(options, 'use_psf')
    options.use_psf = false;
end
if ~isfield(options, 'save_iter')
    options.save_iter = false;
end
if ~isfield(options, 'apply_acceleration')
    options.apply_acceleration = false;
end
if ~isfield(options, 'deblurring')
    options.deblurring = false;
end
if ~isfield(options, 'use_64bit_atomics')
    options.use_64bit_atomics = false;
end
if ~isfield(options, 'use_CUDA')
    options.use_CUDA = false;
end
if ~isfield(options, 'n_rays_transaxial')
    options.n_rays_transaxial = 1;
end
if ~isfield(options, 'n_rays_axial')
    options.n_rays_axial = 1;
end
if ~isfield(options, 'cpu_to_gpu_factor')
    options.cpu_to_gpu_factor = 1;
end
if options.projector_type == 5
    if ~isfield(options,'meanFP')
        options.meanFP = false;
    end
    if ~isfield(options,'meanBP')
        options.meanBP = false;
    end
end
end