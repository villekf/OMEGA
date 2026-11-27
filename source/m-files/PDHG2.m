function [output, options, im_vectors] = PDHG2(im, rhs, im_vectors, options, iter, subIter, ii)
%PDHG2 Final PDHG (and PDDY) primal update (PDHG2 in algorithms.h)
if (options.PDDY)
    % PDDY (Salim, Condat, Mishchenko and Richtarik, JOTA 2022), one subset. On entry im is x = p - tau * A^T y+ (p is
    % kept in uCP, see computeOSEstimates) and rhs = A^T y+ + grad R(x). The iteration is
    %   s+ = prox_+(x - tau * (A^T y+ + grad R(x))),  p+ = p + s+ - x,
    % and the next forward projection is taken at s+. prox_+ is the positivity clamp if enforcePositivity, otherwise the
    % identity. No image preconditioning or adaptive steps (rejected for PDDY)
    s = im - options.tauCP(ii) * rhs;
    if (options.enforcePositivity)
        s(s < options.epps) = options.epps;
    end
    if iscell(im_vectors.uCP)
        im_vectors.uCP{ii} = im_vectors.uCP{ii} + (s - im);
    else
        im_vectors.uCP = im_vectors.uCP + (s - im);
    end
    output = s;
    return
end
kk = (iter - 1) * options.subsets + subIter;
if (options.PDAdaptiveType == 1)
    im_old = im;
end
% The image-based preconditioners are applied to the main volume only
if (ii == 1)
    [rhs, options] = applyImagePreconditioning(options, rhs, im, kk, ii);
end
if (options.subsets > 1)
    if (options.verbose >= 3)
        disp("Using PDHG w/ subsets");
    end
    output = im - options.tauCP(ii) * rhs;
    if (options.enforcePositivity)
        output(output < options.epps) = options.epps;
    end
else
    if (options.verbose >= 3)
        disp("Using PDHG W/O subsets");
    end
    if iscell(im_vectors.uCP)
        uPrev = im_vectors.uCP{ii};
        im_vectors.uCP{ii} = im_vectors.uCP{ii} - options.tauCP(ii) * rhs;
        if (options.enforcePositivity)
            im_vectors.uCP{ii}(im_vectors.uCP{ii} < options.epps) = options.epps;
        end
        output = im_vectors.uCP{ii} + options.thetaCP(kk) * (im_vectors.uCP{ii} - uPrev);
    else
        uPrev = im_vectors.uCP;
        im_vectors.uCP = im_vectors.uCP - options.tauCP(ii) * rhs;
        if (options.enforcePositivity)
            im_vectors.uCP(im_vectors.uCP < options.epps) = options.epps;
        end
        output = im_vectors.uCP + options.thetaCP(kk) * (im_vectors.uCP - uPrev);
    end
end
% The step sizes are not adapted during the measurement-based filtering phase (C++: 0-based subiteration >= filterIter).
% computeForwardStep disables precondTypeMeas(2) after the filtering phase
if ((options.precondTypeMeas(2) && kk - 1 >= options.filteringIterations) || ~options.precondTypeMeas(2))
    if (ii == 1 && options.PDAdaptiveType == 1)
        if iscell(im_vectors.rhsCP)
            q = (im_old - output) / options.tauCP(ii) + options.subsets * im_vectors.rhsCP{ii};
        else
            q = (im_old - output) / options.tauCP(ii) + options.subsets * im_vectors.rhsCP;
        end
        w = dot((im_old - output), q) / (norm(im_old - output) * norm(q));
        if (w < 0)
            options.tauCP(ii) = options.tauCP(ii) / (1 + options.alphaCP(ii));
            options.sigmaCP(ii) = options.sigmaCP(ii) * (1 + options.alphaCP(ii));
            options.alphaCP(ii) = options.alphaCP(ii) * 0.99;
        elseif (w >= .999)
            options.sigmaCP(ii) = options.sigmaCP(ii) / (1 + options.alphaCP(ii));
            options.tauCP(ii) = options.tauCP(ii) * (1 + options.alphaCP(ii));
            options.alphaCP(ii) = options.alphaCP(ii) * 0.99;
        end
        options.sigma2CP(ii) = options.sigmaCP(ii);
    end
end
