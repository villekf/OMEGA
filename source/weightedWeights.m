function options = weightedWeights(options)
%WEIGHTEDWEIGHTS Compute weights for the weighted mean prior
%   options.weights are needed as the input data. They can be formed with
%   computeWeights.

if isempty(options.weighted_weights)
    distX = options.FOVa_x/double(options.Nx);
    kerroin = sqrt(2)*distX;
    options.weighted_weights = kerroin.*options.weights;
    options.weighted_weights(isinf(options.weighted_weights)) = options.weighted_center_weight;
    options.weighted_weights = options.weighted_weights ./ sum(options.weighted_weights(:));
end
options.w_sum = sum(options.weighted_weights(:));
if options.implementation == 2
    options.weighted_weights = single(options.weighted_weights);
    options.mean_type = int32(options.mean_type);
    options.w_sum = single(options.w_sum);
end
options.weighted_weights = reshape(options.weighted_weights, options.Ndx*2 + 1, options.Ndy*2 + 1, options.Ndz*2 + 1);