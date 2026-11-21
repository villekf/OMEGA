function [x, z, viewWeights] = prepareSPECTListmodeSensitivity(options)
%PREPARESPECTLISTMODESENSITIVITY Return full-view geometry and frame weights.

if isfield(options, 'xSens') && isfield(options, 'zSens') && ...
		~isempty(options.xSens) && ~isempty(options.zSens)
	x = options.xSens;
	z = options.zSens;
else
	[x, z] = get_coordinates_SPECT(options);
end

nViews = numel(x) / 6;
if nViews ~= floor(nViews) || numel(z) ~= 2 * nViews
	error('SPECT sensitivity geometry must contain six x coordinates and two z coordinates per projection frame.')
end

nTimesteps = double(options.Nt);
if nTimesteps < 1 || nTimesteps ~= floor(nTimesteps)
	error('Nt must be a positive integer when computing listmode SPECT sensitivity.')
end

if isfield(options, 'sensitivityViewWeights') && ~isempty(options.sensitivityViewWeights)
	providedWeights = options.sensitivityViewWeights;
	if ~isnumeric(providedWeights) || ~isreal(providedWeights) || ...
			size(providedWeights, 1) ~= nViews || size(providedWeights, 2) ~= nTimesteps || ...
			any(~isfinite(providedWeights(:))) || any(providedWeights(:) < 0)
		error('sensitivityViewWeights must be a finite, non-negative nViews-by-Nt matrix.')
	end
	viewWeights = single(providedWeights);
	return
end

frameIndex = [];
if isfield(options, 'temporalBinIndex') && ~isempty(options.temporalBinIndex)
	frameIndex = double(options.temporalBinIndex(:));
	if numel(frameIndex) ~= nViews || any(~isfinite(frameIndex)) || ...
			any(frameIndex ~= floor(frameIndex)) || any(frameIndex < 0) || any(frameIndex >= nTimesteps)
		error('temporalBinIndex must contain one zero-based timeframe index per SPECT sensitivity projection frame.')
	end
end

viewWeights = ones(nViews, nTimesteps, 'single');

% When acquisition intervals and timeframe boundaries are available, split
% each view's exposure by its temporal overlap. This handles both several
% views inside one frame and one view that spans multiple frames.
[captureStart, captureEnd, hasCaptureIntervals] = captureIntervals(options, nViews);
[frameStart, frameEnd, hasFrameIntervals] = timeframeIntervals(options, nTimesteps, captureStart, captureEnd);
canSplitByTime = hasCaptureIntervals && hasFrameIntervals && ...
	~(isfield(options, 'gated') && options.gated) && ...
	~(isfield(options, 'dynamicGated') && options.dynamicGated);

if canSplitByTime
	duration = captureEnd - captureStart;
	positiveDuration = duration > 0;
	viewWeights(:) = 0;
	for timestep = 1:nTimesteps
		overlap = max(0, min(captureEnd, frameEnd(timestep)) - max(captureStart, frameStart(timestep)));
		valid = positiveDuration & isfinite(overlap);
		viewWeights(valid, timestep) = single(overlap(valid) ./ duration(valid));
	end
	% Zero-duration capture intervals do not define a fractional exposure. Use
	% the loader's per-frame assignment for those rows.
	zeroDuration = ~positiveDuration;
	if any(zeroDuration)
		if isempty(frameIndex)
			viewWeights(zeroDuration, :) = 1;
		else
			for timestep = 1:nTimesteps
				viewWeights(zeroDuration & frameIndex == timestep - 1, timestep) = 1;
			end
		end
	end
elseif ~isempty(frameIndex)
	% Timing boundaries may be unavailable for gated data. In that case, use
	% the loader's explicit timeframe assignment and distribute repeated
	% copies of a view according to their recorded integration durations.
	projectionWeights = repeatedViewWeights(options, nViews, frameIndex);
	viewWeights(:) = 0;
	for timestep = 1:nTimesteps
		inTimestep = frameIndex == timestep - 1;
		viewWeights(inTimestep, timestep) = projectionWeights(inTimestep);
	end
elseif nTimesteps > 1
	error(['Dynamic listmode SPECT sensitivity requires per-view timing, temporalBinIndex, ' ...
		'or explicit sensitivityViewWeights for every frame.'])
end
end

function [captureStart, captureEnd, valid] = captureIntervals(options, nViews)
captureStart = [];
captureEnd = [];
valid = false;
startField = firstMatchingField(options, {'measurementStartMs', 'frameStartMs'}, nViews);
endField = firstMatchingField(options, {'measurementEndMs', 'frameEndMs'}, nViews);
if isempty(startField) || isempty(endField)
	return
end
captureStart = double(options.(startField)(:));
captureEnd = double(options.(endField)(:));
valid = all(isfinite(captureStart)) && all(isfinite(captureEnd)) && all(captureEnd >= captureStart);
end

function [frameStart, frameEnd, valid] = timeframeIntervals(options, nTimesteps, captureStart, captureEnd)
frameStart = [];
frameEnd = [];
valid = false;

if isfield(options, 'dynamicPartitionStartMs') && isfield(options, 'dynamicPartitionEndMs') && ...
		numel(options.dynamicPartitionStartMs) == nTimesteps && numel(options.dynamicPartitionEndMs) == nTimesteps
	frameStart = double(options.dynamicPartitionStartMs(:));
	frameEnd = double(options.dynamicPartitionEndMs(:));
	valid = validIntervals(frameStart, frameEnd, nTimesteps);
	if valid
		return
	end
end

if isfield(options, 'timeframes') && numel(options.timeframes) == nTimesteps + 1 && ~isempty(captureStart)
	boundaries = double(options.timeframes(:));
	% OMEGA timeframes are expressed in seconds relative to acquisition start.
	boundaries = min(captureStart) + boundaries .* 1000;
	frameStart = boundaries(1:end-1);
	frameEnd = boundaries(2:end);
	valid = all(diff(boundaries) > 0) && ~any(isnan(boundaries)) && ...
		~any(isinf(boundaries(1:end-1))) && ~any(isnan(frameStart)) && ~any(isnan(frameEnd));
	if valid
		return
	end
end

if nTimesteps == 1 && ~isempty(captureStart)
	frameStart = min(captureStart);
	frameEnd = max(captureEnd);
	valid = isfinite(frameStart) && isfinite(frameEnd) && frameEnd >= frameStart;
end
end

function valid = validIntervals(frameStart, frameEnd, nTimesteps)
valid = numel(frameStart) == nTimesteps && numel(frameEnd) == nTimesteps && ...
	all(isfinite(frameStart)) && all(frameEnd >= frameStart) && ...
	all(diff(frameStart) >= 0);
end

function field = firstMatchingField(options, candidates, expectedLength)
field = '';
for index = 1:numel(candidates)
	name = candidates{index};
	if isfield(options, name) && numel(options.(name)) == expectedLength
		field = name;
		return
	end
end
end

function weights = repeatedViewWeights(options, nViews, frameIndex)
weights = ones(nViews, 1, 'single');
if ~isfield(options, 'measurementDurationMs') || numel(options.measurementDurationMs) ~= nViews || ...
		~isfield(options, 'viewIndex') || numel(options.viewIndex) ~= nViews
	return
end

duration = double(options.measurementDurationMs(:));
viewIndex = double(options.viewIndex(:));
if any(~isfinite(duration)) || any(duration < 0) || any(~isfinite(viewIndex))
	return
end

if isfield(options, 'DetectorVector') && numel(options.DetectorVector) == nViews
	panelIndex = double(options.DetectorVector(:));
elseif isfield(options, 'panelIndex') && numel(options.panelIndex) == nViews
	panelIndex = double(options.panelIndex(:));
else
	panelIndex = zeros(nViews, 1);
end
if isempty(frameIndex)
	frameIndex = zeros(nViews, 1);
end
[~, ~, viewPanelGroup] = unique([viewIndex, panelIndex, frameIndex], 'rows');
groupDuration = accumarray(viewPanelGroup, duration, [], @sum);
validDuration = groupDuration(viewPanelGroup) > 0;
weights(:) = 0;
weights(validDuration) = single(duration(validDuration) ./ groupDuration(viewPanelGroup(validDuration)));
weights(~validDuration) = 1;
end
