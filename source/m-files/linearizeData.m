function options = linearizeData(options)
%LINEARIZEDATA Linearizes the input CT data
%   The measurement is floored at flat * 1e-6 before the logarithm, such that
%   zero or negative (corrected) intensities do not produce Inf/NaN
flat = single(options.flat);
if iscell(options.SinM)
    for kk = 1 : numel(options.SinM)
        options.SinM{kk} = log(flat ./ max(single(options.SinM{kk}), flat * 1e-6));
    end
else
    options.SinM = log(flat ./ max(single(options.SinM), flat * 1e-6));
end
end
