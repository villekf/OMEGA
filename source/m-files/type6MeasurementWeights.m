function weights = type6MeasurementWeights(options, norm_input, corr_input, views)
%TYPE6MEASUREMENTWEIGHTS Build view-specific multiplicative measurement weights.
% This mirrors type6Weights in source/cpp/functions.hpp: normalization may be
% one image per detector head or one image per projection, while scatter and
% the extra correction vector are indexed in measurement order.

nRows = double(options.nRowsD);
nCols = double(options.nColsD);
nDet = nRows * nCols;
views = double(views(:).');
nViews = numel(views);
weights = ones(nRows, nCols, nViews, options.cType);

if options.normalization_correction && ~isempty(norm_input) && numel(norm_input) > 1
    norm_input = reshape(norm_input, nRows, nCols, []);
    nNorm = size(norm_input, 3);
    nHeads = 1;
    if isfield(options, 'nHeads') && ~isempty(options.nHeads)
        nHeads = double(options.nHeads);
    end
    normZ = nNorm;
    if isfield(options, 'normZ') && ~isempty(options.normZ)
        normZ = double(options.normZ);
    end
    if options.SPECT && normZ == nHeads && nNorm >= nHeads
        headIndices = zeros(1, nViews);
        detectorVector = double(options.DetectorVector(:));
        geomStart = 0;
        if isfield(options, 'currentTimestep') && isscalar(options.currentTimestep) && ...
                options.currentTimestep > 0 && isfield(options, 'nMeas') && numel(options.nMeas) > 1
            geomStart = double(options.nMeas(options.currentTimestep + 1));
        end
        for kk = 1:nViews
            geomIndex = geomStart + views(kk);
            if geomIndex <= numel(detectorVector)
                headIndices(kk) = detectorVector(geomIndex) + 1;
            else
                headIndices(kk) = 1;
            end
        end
        headIndices = max(1, min(nNorm, headIndices));
        for kk = 1:nViews
            weights(:, :, kk) = weights(:, :, kk) .* cast(norm_input(:, :, headIndices(kk)), 'like', weights);
        end
    else
        normCount = nNorm;
        if normCount == nViews
            normViews = 1:nViews;
        elseif max(views) <= normCount
            normViews = views;
        else
            error('type6MeasurementWeights:NormalizationSize', ...
                'Normalization data do not contain the requested type-6 views.');
        end
        weights = weights .* cast(norm_input(:, :, normViews), 'like', weights);
    end
end

useCorr = isfield(options, 'additionalCorrection') && options.additionalCorrection;
if isfield(options, 'scatter') && options.scatter == 1
    useCorr = true;
end
if useCorr && ~isempty(corr_input) && numel(corr_input) > 1
    corr_input = reshape(corr_input, nRows, nCols, []);
    nCorr = size(corr_input, 3);
    if nCorr == nViews
        corrViews = 1:nViews;
    elseif max(views) <= nCorr
        corrViews = views;
    else
        error('type6MeasurementWeights:CorrectionSize', ...
            'Additional-correction data do not contain the requested type-6 views.');
    end
    weights = weights .* cast(corr_input(:, :, corrViews), 'like', weights);
end

weights = weights .* cast(options.global_correction_factor, 'like', weights);
end
