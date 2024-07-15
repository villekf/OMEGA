function [im_vectors, im, options] = FISTAL1(im, rhs, options, im_vectors, beta, iter, ii)
[im_vectors, im, options] = FISTA(im, rhs, options, im_vectors, iter, ii);
a = w_vec.tauCP(ii) * beta;
im(abs(im) <= a) = 0;
apu = -sign(im);
apu(apu == 0) = 1;
im = apu * max(abs(im) - a, 0);