function options = RDPWeights(options)
%RDPWEIGHTS Compute weights for the relative difference prior

if ~isfield(options,'weights_RDP') || isempty(options.weights_RDP)
    if ~isfield(options,'weights') || isempty(options.weights)
        options = computeWeights(options);
    end
    options.weights_RDP = options.weights/sum(options.weights(~isinf(options.weights)));
    
end
if options.implementation == 2
    options.weights_RDP = single(options.weights_RDP);
    if any(isinf(options.weights_RDP))
        options.weights_RDP(isinf(options.weights_RDP)) = [];
    else
        options.weights_RDP(ceil(numel(options.weights_RDP) / 2)) = [];
    end
end