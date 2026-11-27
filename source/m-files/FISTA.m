function [im_vectors, im, options] = FISTA(im, rhs, options, im_vectors, iter, osa_iter, ii, varargin)
% FISTA (Beck & Teboulle 2009), optionally with L1 regularization (soft thresholding with threshold varargin{1}):
% x_{k+1} = prox(y_k - tau * grad f(y_k)),  y_{k+1} = x_{k+1} + beta_k * (x_{k+1} - x_k)
% im_vectors.uFISTA holds the extrapolated point y_k and im_vectors.xFISTA the previous iterate x_k (FISTA in algorithms.h).
% Without subsets im is y_k (initializationStep copies uFISTA to recApu on every iteration). With subsets the first
% subsets - 1 subiterations are plain (proximal) gradient steps starting from y_k (copied at the first subset of each
% iteration) and the momentum step (outer iteration k) is taken only on the last subset
% iter is the outer iteration, osa_iter the subset (sub-iteration) number, both 1-based
if nargin >= 8
    l1Threshold = varargin{1};
else
    l1Threshold = 0;
end
kk = (iter - 1) * options.subsets + osa_iter;
[rhs, options] = applyImagePreconditioning(options, rhs, im, kk, ii);
if (options.subsets == 1 || osa_iter == options.subsets)
    if (options.subsets == 1)
        if iscell(im_vectors.uFISTA)
            im = im_vectors.uFISTA{ii} - options.tauCP(ii) * rhs;
        else
            im = im_vectors.uFISTA - options.tauCP(ii) * rhs;
        end
    else
        im = im - options.tauCP(ii) * rhs;
    end
    if iscell(im_vectors.xFISTA)
        xPrev = im_vectors.xFISTA{ii};
    else
        xPrev = im_vectors.xFISTA;
    end
    % Proximal step (soft thresholding) and positivity before the momentum extrapolation
    if l1Threshold > 0
        im = sign(im) .* max(abs(im) - l1Threshold, 0);
    end
    if options.enforcePositivity
        im = max(im, options.epps);
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
            % Outer iteration count (C++ iter + 1)
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
else
    % Plain gradient step on the subsets 1, ..., subsets - 1 (the positivity is enforced after the subiteration)
    im = im - options.tauCP(ii) * rhs;
    if l1Threshold > 0
        im = sign(im) .* max(abs(im) - l1Threshold, 0);
    end
end
end
