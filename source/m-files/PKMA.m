function output = PKMA(im, rhs, options, iter, osa_iter, ii)

kk = (iter - 1) * options.subsets + osa_iter;
rhs = applyImagePreconditioning(options, rhs, im, kk, ii);
im_apu = im - options.lambda(iter) * rhs;
if (options.enforcePositivity)
    im_apu(im_apu < options.epps) = options.epps;
end
output = (1 - options.alpha_PKMA(kk)) * im + options.alpha_PKMA(kk) * im_apu;