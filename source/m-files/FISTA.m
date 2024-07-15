function [im_vectors, im, options] = FISTA(im, rhs, options, im_vectors, iter, ii)
uPrev = im;
rhs = applyImagePreconditioning(options, rhs, im, iter, ii);
if iscell(im_vectors.uFISTA)
    im = im_vectors.uFISTA{ii} - options.tauCP(ii) * rhs;
else
    im = im_vectors.uFISTA - options.tauCP(ii) * rhs;
end
if (ii == 1)
    if ~isfield(options,'tNFista')
        options.tNFista = 1;
    end
    it = iter;
    options.betaFISTA = (it - 1) / (it + 2);
    if (options.betaFISTA <= 0)
        options.tFISTA = (1 + sqrt(1 + 4 * options.tNFista * options.tNFista)) / 2;
        options.betaFISTA = (options.tNFista - 1) / options.tFISTA;
        options.tNFista = options.tFISTA;
    end
end
if iscell(im_vectors.uFISTA)
    im_vectors.uFISTA{ii} = im + options.betaFISTA * (im - uPrev);
else
    im_vectors.uFISTA = im + options.betaFISTA * (im - uPrev);
end