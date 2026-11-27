function [im_vectors, im, options] = FISTAL1(im, rhs, options, im_vectors, beta, iter, osa_iter, ii)
% FISTA with L1 regularization: the soft thresholding (threshold tau * beta) is the proximal step of FISTA,
% applied before the momentum extrapolation (FISTAL1 in algorithms.h)
a = options.tauCP(ii) * beta;
[im_vectors, im, options] = FISTA(im, rhs, options, im_vectors, iter, osa_iter, ii, a);
end
