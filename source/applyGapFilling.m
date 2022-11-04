function options = applyGapFilling(options)
%APPLYGAPFILLING Performs sinogram gap filling for input data
%   Utility function

[~, ~, xp, yp] = detector_coordinates(options);
for llo = 1 : options.partitions
    if llo == 1
        gaps = [];
    end
    if options.partitions > 1
        Sin = options.SinM{llo};
    else
        Sin = options.SinM;
    end
    [Sin, gaps] = gapFilling(options, Sin, xp, yp, llo, gaps);
    if options.partitions > 1
        options.SinM{llo} = Sin;
    else
        options.SinM = Sin;
    end
end
end