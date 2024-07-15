function options = linearizeData(options)
%LINEARIZEDATA Linearizes the input CT data
%   Detailed explanation goes here
if iscell(options.SinM)
    for kk = 1 : numel(options.SinM)
        options.SinM{kk} = log(single(options.flat) ./ single(options.SinM{kk}));
    end
else
    options.SinM = log(single(options.flat) ./ single(options.SinM));
end