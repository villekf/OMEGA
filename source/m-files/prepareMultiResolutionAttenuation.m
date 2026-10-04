function options = prepareMultiResolutionAttenuation(options)
%PREPAREMULTIRESOLUTIONATTENUATION Pack centered full-FOV attenuation maps.
%   The input remains an image-domain map. For multi-resolution geometry,
%   the packed layout is [fine; coarse] for each frame, with both grids
%   centered at the world origin. Shifted emission slabs can enlarge the
%   support symmetrically; padding outside the input map is air (mu = 0).

if ~isfield(options, 'useMultiResolutionVolumes') || ~options.useMultiResolutionVolumes || ...
        ~isfield(options, 'nMultiVolumes') || options.nMultiVolumes < 1 || ...
        ~isfield(options, 'attenuation_correction') || ~options.attenuation_correction || ...
        ~isfield(options, 'CT_attenuation') || ~options.CT_attenuation || ...
        ~isfield(options, 'vaimennus') || isempty(options.vaimennus)
    return
end

if ~isfield(options, 'SPECT') || ~options.SPECT
    error('Multi-resolution image-domain attenuation is currently supported only for SPECT.')
end
if ~isfield(options, 'projector_type') || ~ismember(options.projector_type, [1, 2, 11, 21, 22])
    error('Multi-resolution image-domain attenuation requires SPECT Siddon or orthogonal projectors (types 1, 2, 11, 21, or 22).')
end
if isfield(options, 'implementation') && (options.implementation == 1 || options.implementation == 4)
    error('Multi-resolution image-domain attenuation is not implemented by the CPU projector backends.')
end
if isfield(options, 'implementation') && options.implementation == 2 && ...
        isfield(options, 'use_CPU') && options.use_CPU
    error('Multi-resolution image-domain attenuation is not implemented by the native C++ CPU backend.')
end
if ismac && isfield(options, 'implementation') && ismember(options.implementation, [2, 3, 5])
    error('Multi-resolution image-domain attenuation is not implemented by the Metal/MPS projector backend.')
end

if isfield(options, 'imageAttenuationIsMultiResolution') && options.imageAttenuationIsMultiResolution && ...
        isfield(options, 'imageAttenuationGridDims') && isfield(options, 'imageAttenuationGridSpacing') && ...
        isfield(options, 'imageAttenuationGridOrigin')
    return
end

baseDims = double([options.NxFull, options.NyFull, options.NzFull]);
if isfield(options, 'imageAttenuationFineDims')
    baseDims = double(options.imageAttenuationFineDims(:).');
end
if any(~isfinite(baseDims)) || any(baseDims < 1) || any(baseDims ~= floor(baseDims))
    error('The full-FOV fine attenuation grid dimensions must be positive integers.')
end

if isfield(options, 'imageAttenuationFineSpacing')
    fineSpacing = double(options.imageAttenuationFineSpacing(:).');
elseif isfield(options, 'imageAttenuationFineFOV')
    fineFOV = double(options.imageAttenuationFineFOV(:).');
    fineSpacing = fineFOV ./ baseDims;
else
    fineFOV = double([options.FOVa_x(1), options.FOVa_y(1), options.axial_fov(1)]);
    fineSpacing = fineFOV ./ baseDims;
end
if any(~isfinite(fineSpacing)) || any(fineSpacing <= 0)
    error('The full-FOV fine attenuation grid spacing must be positive.')
end

nVolumes = double(options.nMultiVolumes) + 1;
origins = [reshape(double(options.bx(1:nVolumes)), [], 1), reshape(double(options.by(1:nVolumes)), [], 1), reshape(double(options.bz(1:nVolumes)), [], 1)];
volumeDims = [reshape(double(options.Nx(1:nVolumes)), [], 1), reshape(double(options.Ny(1:nVolumes)), [], 1), reshape(double(options.Nz(1:nVolumes)), [], 1)];
volumeSpacing = [reshape(double(options.dx(1:nVolumes)), [], 1), reshape(double(options.dy(1:nVolumes)), [], 1), reshape(double(options.dz(1:nVolumes)), [], 1)];
if size(origins, 1) ~= nVolumes || any(~isfinite(origins(:))) || any(~isfinite(volumeSpacing(:)))
    error('Multi-resolution attenuation requires geometry for each active volume.')
end
spacingTolerance = max(1e-5, abs(fineSpacing) .* 1e-5);
sameResolution = all(all(abs(volumeSpacing - repmat(fineSpacing, nVolumes, 1)) <= ...
    repmat(spacingTolerance, nVolumes, 1)));
tileMinimum = min(origins, [], 1);
tileMaximum = max(origins + volumeDims .* volumeSpacing, [], 1);

% Preserve the pre-split grid spacing and world center. Expand the centered
% grid only when an offset slab extends beyond that support; matching parity
% allows the original map to be padded equally on both sides.
fineDims = baseDims;
tileExtent = 2 .* max(abs(tileMinimum), abs(tileMaximum));
% The projector stores origins and spacings as single precision. Their
% multiply/add roundoff can put an otherwise exact FOV boundary a few ulps
% outside the input map. Discount only that coordinate representation error
% before ceil(), so it cannot trigger symmetric padding by two full voxels.
coordinateMagnitude = max(abs(tileMinimum), abs(tileMaximum));
coordinateRoundoff = 2 .* double(eps(single(coordinateMagnitude)));
requiredDims = ceil(tileExtent ./ fineSpacing - coordinateRoundoff ./ fineSpacing);
fineDims = max(fineDims, requiredDims);
oddPadding = mod(fineDims - baseDims, 2) ~= 0;
fineDims(oddPadding) = fineDims(oddPadding) + 1;
fineExtent = fineDims .* fineSpacing;
fineOrigin = -fineExtent ./ 2;

coarseNominalSpacing = [double(options.dx(2)), double(options.dy(2)), double(options.dz(2))];
if any(~isfinite(coarseNominalSpacing)) || any(coarseNominalSpacing <= 0)
    error('Coarse multi-resolution voxel spacings must be positive.')
end
if sameResolution
    % At scale 1, use one identical grid for both emission resolutions.
    % Tiny precision differences must not make ceil() add a voxel and shift
    % the coarse attenuation map by half a voxel.
    coarseDims = fineDims;
    coarseSpacing = fineSpacing;
    coarseOrigin = fineOrigin;
else
    coarseDims = max(1, ceil(fineExtent ./ coarseNominalSpacing - 1e-10));
    coarseSpacing = coarseNominalSpacing;
    coarseExtent = coarseDims .* coarseSpacing;
    coarseOrigin = -coarseExtent ./ 2;
end

% Normalize frame layout. 3D is static; a 4D map has time in its final axis.
attenuation = options.vaimennus;
if iscell(attenuation)
    if isempty(attenuation)
        error('The image-domain attenuation map is empty.')
    end
    if numel(attenuation) == 1
        attenuation = attenuation{1};
    else
        frameVolumes = attenuation(:);
        for frame = 1:numel(frameVolumes)
            volume = full(frameVolumes{frame});
            if numel(volume) == prod(baseDims)
                volume = reshape(volume, baseDims);
            elseif ~isequal(double(size3(volume)), baseDims)
                error('Each image-domain attenuation frame must match the full-FOV fine-grid dimensions.')
            end
            frameVolumes{frame} = volume;
        end
        attenuation = cat(4, frameVolumes{:});
    end
end
attenuation = full(attenuation);

if isvector(attenuation) && numel(attenuation) ~= prod(baseDims)
    nFramesFlat = numel(attenuation) / prod(baseDims);
    if nFramesFlat ~= floor(nFramesFlat)
        error('A flattened image-domain attenuation map must match the full-FOV grid dimensions.')
    end
    attenuation = reshape(attenuation, [baseDims, nFramesFlat]);
elseif isvector(attenuation) && numel(attenuation) == prod(baseDims)
    attenuation = reshape(attenuation, baseDims);
end

if ndims(attenuation) <= 3
    attenuation = reshape(attenuation, [size(attenuation, 1), size(attenuation, 2), size(attenuation, 3), 1]);
    dynamic = false;
else
    nFrames = size(attenuation, 4);
    if nFrames == 1
        dynamic = false;
    elseif nFrames == double(options.Nt)
        dynamic = true;
    else
        error('A 4D attenuation map must contain one static frame or Nt=%d frames; got %d.', double(options.Nt), nFrames)
    end
end

frames = zeros([fineDims, size(attenuation, 4)], 'like', attenuation);
for tt = 1:size(attenuation, 4)
    frame = attenuation(:, :, :, tt);
    if ~isequal(double(size3(frame)), baseDims)
        frame = resampleCentered(frame, baseDims, double(size3(frame)));
    end

    % Match existing MATLAB map orientation controls. SPECT map transforms
    % are deferred from projectorClass for multi-resolution maps so the same
    % transform is applied here to every static or dynamic frame.
    if options.SPECT
        if isfield(options, 'offangle') && options.offangle ~= 0
            frame = imrotate(frame, options.offangle, 'bilinear', 'crop');
        end
        if isfield(options, 'flipImageX') && options.flipImageX
            frame = flip(frame, 2);
        end
        if isfield(options, 'flipImageY') && options.flipImageY
            frame = flip(frame, 1);
        end
        if isfield(options, 'flipImageZ') && options.flipImageZ
            frame = flip(frame, 3);
        end
    else
        if isfield(options, 'rotateAttImage') && options.rotateAttImage ~= 0
            frame = rot90(frame, options.rotateAttImage);
        end
        if isfield(options, 'flipAttImageXY') && options.flipAttImageXY
            frame = fliplr(frame);
        end
        if isfield(options, 'flipAttImageZ') && options.flipAttImageZ
            frame = flip(frame, 3);
        end
    end
    if ~isequal(double(size3(frame)), baseDims)
        frame = resampleCentered(frame, baseDims, double(size3(frame)));
    end

    if any(fineDims > baseDims)
        padded = zeros(fineDims, 'like', frame);
        starts = floor((fineDims - baseDims) ./ 2) + 1;
        ix = starts(1) : starts(1) + baseDims(1) - 1;
        iy = starts(2) : starts(2) + baseDims(2) - 1;
        iz = starts(3) : starts(3) + baseDims(3) - 1;
        padded(ix, iy, iz) = frame;
        frame = padded;
    end
    if isfield(options, 'attIncm') && options.attIncm
        frame = frame ./ 10;
    end
    frames(:, :, :, tt) = frame;
end

packed = cell(2 * size(frames, 4), 1);
for tt = 1:size(frames, 4)
    fine = frames(:, :, :, tt);
    if sameResolution
        coarse = fine;
    else
        coarse = resampleCenteredSpacing(fine, coarseDims, fineSpacing, coarseSpacing, fineOrigin, coarseOrigin);
    end
    packed{2 * tt - 1} = fine(:);
    packed{2 * tt} = coarse(:);
end
if ~dynamic
    packed = packed(1:2);
end

if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
    options.vaimennus = single(vertcat(packed{:}));
else
    options.vaimennus = double(vertcat(packed{:}));
end
options.imageAttenuationGridDims = uint32([fineDims; coarseDims]);
options.imageAttenuationGridSpacing = single([fineSpacing; coarseSpacing]);
options.imageAttenuationGridOrigin = single([fineOrigin; coarseOrigin]);
options.imageAttenuationMapSizes = uint64([prod(fineDims); prod(coarseDims)]);
options.imageAttenuationFrames = uint32(size(frames, 4) * double(dynamic) + double(~dynamic));
options.imageAttenuationIsMultiResolution = true;
options.imageAttenuationIsDynamic = dynamic;
end

function dims = size3(volume)
dims = [size(volume, 1), size(volume, 2), size(volume, 3)];
end

function output = resampleCentered(input, outputDims, inputDims)
% Resample with the volume center fixed and nearest-edge values at the bounds.
if isequal(double(outputDims), double(inputDims))
    output = input;
    return
end
query = cell(1, 3);
for axis = 1:3
    query{axis} = ((1:outputDims(axis)) - 0.5) .* inputDims(axis) ./ outputDims(axis) + 0.5;
    query{axis} = min(max(query{axis}, 1), inputDims(axis));
end
[xq, yq, zq] = ndgrid(query{1}, query{2}, query{3});
output = interpn(input, xq, yq, zq, 'linear');
end

function output = resampleCenteredSpacing(input, outputDims, inputSpacing, outputSpacing, inputOrigin, outputOrigin)
% Sample each output voxel center in world coordinates. Outside the source
% map, use mu=0 (air), including the fractional boundary support introduced
% when the exact coarse voxel spacing needs one extra voxel.
query = cell(1, 3);
for axis = 1:3
    outputCenters = outputOrigin(axis) + ((0:(outputDims(axis) - 1)) + 0.5) .* outputSpacing(axis);
    query{axis} = (outputCenters - inputOrigin(axis)) ./ inputSpacing(axis) + 0.5;
end
[xq, yq, zq] = ndgrid(query{1}, query{2}, query{3});
output = interpn(input, xq, yq, zq, 'linear', 0);
end
