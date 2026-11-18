function [im_vectors, im, options] = FISTAL1(im, rhs, options, im_vectors, beta, iter, ii)
% FISTA with L1 regularization: the soft thresholding (threshold tau * beta) is the proximal step of FISTA,
% applied before the momentum extrapolation
a = options.tauCP(ii) * beta;
[im_vectors, im, options] = FISTA(im, rhs, options, im_vectors, iter, ii, a);
end
