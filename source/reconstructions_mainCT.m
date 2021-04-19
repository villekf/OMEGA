function pz = reconstructions_mainCT(options)
%RECONSTRUCTIONS_MAINCT Utility function for CT reconstruction
%   This function simply adds PET specific variables that are needed for
%   error-free functionality and also converts some input variables to the
%   correct format (e.g. the projection angles are converted to radians if
%   they are input as degrees).


var = recNames(4);
ll = 0;
kk = 1;
while ll == 0 && kk <= numel(var)
    ll = ll + options.(var{kk});
    kk = kk +1;
end
OS_bool = ll > 0;
if options.implementation == 2 && options.MLEM && OS_bool
    warning('MLEM reconstruction does not work when non-MLEM method is selected as well!')
end
options.NSinos = options.nProjections;
options.TotSinos = options.nProjections;
options.span = 3;
options.segment_table = [];
options.Ndist = options.ySize;
options.Nang = options.xSize;
options.use_raw_data = false;
options.randoms_correction = false;
options.scatter_correction = false;
options.attenuation_correction = false;
options.normalization_correction = false;
if ~isfield(options, 'partitions')
    options.partitions = 1;
end
options.reconstruct_trues = false;
options.reconstruct_scatter = false;
options.machine_name = 'CT';
options.corrections_during_reconstruction = false;
options.precompute_obs_matrix = false;
options.precompute_all = false;
options.compute_normalization = false;
options.diameter = 0;
options.ring_difference = 0;
options.ndist_side = 0;
options.sampling_raw = 1;
options.pseudot = [];
if ~isfield(options, 'nBed')
    options.nBed = 1;
end
options.rings = options.xSize * options.nBed;
options.det_per_ring = options.ySize * options.nProjections;
options.precompute_lor = false;
options.sampling = 1;
options.CT = true;
options.start = 0;
options.end = 0;
options.arc_correction = false;
options.tot_time = 0;
options.fill_sinogram_gaps = false;
options.det_w_pseudo = options.det_per_ring;
options.blocks_per_ring = 1;
options.linear_multip = 0;
options.cryst_per_block = options.xSize * options.ySize;
options.cr_pz = options.dPitch;
if ~isfield(options, 'n_rays_transaxial')
    options.n_rays_transaxial = 1;
end
if ~isfield(options, 'n_rays_axial')
    options.n_rays_axial = 1;
end
if ~isfield(options, 'horizontalOffset')
    options.horizontalOffset = 0;
end
if ~isfield(options, 'verticalOffset')
    options.verticalOffset = 0;
end
if ~isfield(options, 'bedOffset')
    options.bedOffset = 0;
end
if ~isfield(options, 'flip_image')
    options.flip_image = false;
end
if ~isfield(options, 'offangle')
    options.offangle = 0;
end
if max(abs(options.angles(:))) > 2*pi
    options.angles = options.angles * (pi / 180);
end
if ~isfield(options, 'uCenter')
    options.uCenter = [];
end
if ~isfield(options, 'vCenter')
    options.vCenter = [];
end
if (max(options.SinM(:)) > 1 && (min(options.SinM(:)) < 1 || min(options.SinM(:)) > 10)) && options.implementation ~= 1
    options.SinM = single(options.SinM) ./ single(max(options.SinM(:)));
elseif (max(options.SinM(:)) > 1 && (min(options.SinM(:)) < 1 || min(options.SinM(:)) > 10)) && options.implementation == 1
    options.SinM = max(options.SinM(:)) ./ options.SinM;
elseif options.implementation ~= 1 && max(options.SinM(:)) > 1 && min(options.SinM(:)) >= 1
    options.SinM = 1 ./ options.SinM;
end
if options.flip_image
    options.angles = -options.angles;
end
if options.offangle > 0
    options.angles = options.angles + options.offangle;
end
options = OMEGA_error_check(options);

if ~isfield(options,'x') && ~isfield(options,'y') && ~isfield(options,'z') && ~isfield(options,'z_det')
    [options.x,options.y,options.z] = CTDetectorCoordinates(options.angles,options.sourceToDetector,options.sourceToCRot,options.dPitch,options.xSize,...
        options.ySize,options.horizontalOffset,options.verticalOffset,options.bedOffset, options.uCenter, options.vCenter);
    if numel(options.z)/2 > numel(options.angles)
        if size(options.angles,1) == 1
            options.angles = reshape(options.angles, [],1);
        end
        options.angles = repmat(options.angles,numel(options.z)/2/numel(options.angles),1);
    end
end
if isfield(options,'z_det')
    options.z = options.z_det;
    options = rmfield(options,'z_det');
end
% load CBCT_coord.mat x y
% apu = options.x;
% options.x = options.y;
% options.y = apu;
% options.x = fliplr(options.x);
% options.y = fliplr(options.y);
% options.z = fliplr(options.z);
% [options.x,options.y,options.z] = CTDetectorCoordinatesFull(-options.angles,options.sourceToDetector,options.sourceToCRot,options.dPitch,options.xSize,...
%         options.ySize,options.horizontalOffset,options.verticalOffset,options.bedOffset);
if options.implementation == 2 || options.implementation == 3 || options.implementation == 5
    options.x = single(options.x);
    options.y = single(options.y);
    options.z = single(options.z);
else
    options.x = double(options.x);
    options.y = double(options.y);
    options.z = double(options.z);
end
pz = reconstructions_main(options);
end

