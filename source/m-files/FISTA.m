function [im_vectors, im, options] = FISTA(im, rhs, options, im_vectors, iter, ii, varargin)
% FISTA (Beck & Teboulle 2009), optionally with L1 regularization (soft thresholding with threshold varargin{1}):
% x_{k+1} = prox(y_k - tau * grad f(y_k)),  y_{k+1} = x_{k+1} + beta_k * (x_{k+1} - x_k)
% im is the extrapolated point y_k (initializationStep copies uFISTA to recApu), im_vectors.uFISTA holds the
% extrapolated point and im_vectors.xFISTA the previous iterate x_k
if nargin >= 7
    l1Threshold = varargin{1};
else
    l1Threshold = 0;
end
rhs = applyImagePreconditioning(options, rhs, im, iter, ii);
if iscell(im_vectors.uFISTA)
    im = im_vectors.uFISTA{ii} - options.tauCP(ii) * rhs;
    xPrev = im_vectors.xFISTA{ii};
else
    im = im_vectors.uFISTA - options.tauCP(ii) * rhs;
    xPrev = im_vectors.xFISTA;
end
% Proximal step (soft thresholding) before the momentum extrapolation
if l1Threshold > 0
    im = sign(im) .* max(abs(im) - l1Threshold, 0);
end
if options.enforcePositivity
    im(im < options.epps) = options.epps;
end
if (ii == 1)
    if ~isfield(options,'tNFista')
        options.tNFista = 1;
    end
    if isfield(options, 'FISTAType') && options.FISTAType == 1
        options.tFISTA = (1 + sqrt(1 + 4 * options.tNFista * options.tNFista)) / 2;
        options.betaFISTA = (options.tNFista - 1) / options.tFISTA;
        options.tNFista = options.tFISTA;
    else
        it = iter;
        options.betaFISTA = (it - 1) / (it + 2);
        if (options.betaFISTA <= 0)
            options.tFISTA = (1 + sqrt(1 + 4 * options.tNFista * options.tNFista)) / 2;
            options.betaFISTA = (options.tNFista - 1) / options.tFISTA;
            options.tNFista = options.tFISTA;
        end
    end
end
if iscell(im_vectors.uFISTA)
    im_vectors.uFISTA{ii} = im + options.betaFISTA * (im - xPrev);
    im_vectors.xFISTA{ii} = im;
else
    im_vectors.uFISTA = im + options.betaFISTA * (im - xPrev);
    im_vectors.xFISTA = im;
end
end
