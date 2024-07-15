function options = huberWeights(options)
%HUBERWEIGHTS Compute weights for the Huber prior
%   options.weights are needed as the input data. They can be formed with
%   computeWeights.


if isempty(options.weights_huber)
    options.weights_huber = options.weights/sum(options.weights(~isinf(options.weights)));
    options.weights_huber = [options.weights_huber(1:floor(length(options.weights_huber) / 2)); ...
        options.weights_huber(ceil(length(options.weights_huber)/2) + 1 : end)];
end
if options.implementation == 2
    options.weights_huber(isinf(options.weights_huber)) = [];
    options.weights_huber = single(options.weights_huber);
    options.huber_delta = single(options.huber_delta);
else
    options.weights_huber = options.weights_huber * -1;
    options.weights_huber = [options.weights_huber(1 : end/2); abs(sum(options.weights_huber)); options.weights_huber(end/2 + 1: end)];
    options.weights_huber = reshape(options.weights_huber, options.Ndx * 2 + 1, options.Ndy * 2 + 1, options.Ndz * 2 + 1);
end
end

