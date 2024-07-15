function options = quadWeights(options, varargin)
%QUADWEIGHTS Compute weights for the quadratic prior
%   options.weights are needed as the input data. They can be formed with
%   computeWeights.

if nargin >= 2
    empty_weight = varargin{1};
else
    empty_weight = true;
end
if empty_weight
    options.weights_quad = options.weights/sum(options.weights(~isinf(options.weights)));
    if ~options.GGMRF
        options.weights_quad = [options.weights_quad(1:floor(length(options.weights_quad) / 2)); ...
            options.weights_quad(ceil(length(options.weights_quad)/2) + 1 : end)];
    end
else
    options.weights_quad = options.weights;
end
if options.implementation == 2
    if ~options.GGMRF
        options.weights_quad(isinf(options.weights_quad)) = [];
    end
    options.weights_quad = single(options.weights_quad);
    clear weights_quad
else
    options.weights_quad = [options.weights_quad(1 : floor(end/2)); abs(sum(options.weights_quad(~isinf(options.weights_quad)))); options.weights_quad(ceil(end/2) + 1: end)];
    options.weights_quad = reshape(options.weights_quad, options.Ndx * 2 + 1, options.Ndy * 2 + 1, options.Ndz * 2 + 1);
end
end

