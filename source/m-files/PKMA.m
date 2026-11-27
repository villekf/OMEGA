function [output, options] = PKMA(im, rhs, options, iter, osa_iter, ii)
%PKMA Preconditioned Krasnoselskii-Mann algorithm update (PKMA in algorithms.h)
%   The relaxation parameter options.lambda(iter) can be computed
%   (computeRelaxationParameters) and/or scaled (relaxationScaling) here,
%   the modified values are returned in options

kk = (iter - 1) * options.subsets + osa_iter;
[rhs, options] = applyImagePreconditioning(options, rhs, im, kk, ii);
if isfield(options, 'computeRelaxationParameters') && options.computeRelaxationParameters
    if (kk == 1 && ii == 1)
        options.lambda(iter) = norm(im) / norm(rhs) * .25;
    elseif (iter > 1 && osa_iter == 1 && ii == 1)
        % The decay is applied to the initial lambda (C++ iter / 35 + 1 with the 0-based iteration)
        options.lambda(iter) = options.lambda(1) / ((iter - 1) / 35 + 1);
    end
end
if isfield(options, 'relaxationScaling') && options.relaxationScaling
    kerroin = norm(im) / norm(rhs * options.lambda(iter));
    if (kerroin < 1.5 && kerroin > 0)
        options.lambda(iter) = options.lambda(iter) * (kerroin / 1.5);
    end
end
im_apu = im - options.lambda(iter) * rhs;
if (options.enforcePositivity)
    im_apu(im_apu < options.epps) = options.epps;
end
output = (1 - options.alpha_PKMA(kk)) * im + options.alpha_PKMA(kk) * im_apu;
