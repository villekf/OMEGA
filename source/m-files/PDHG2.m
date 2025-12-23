function [output, options, im_vectors] = PDHG2(im, rhs, im_vectors, options, iter, subIter, timestep, ii)
kk = (iter - 1) * options.subsets + subIter;
if (options.PDAdaptiveType == 1)
    im_old = im;
end
rhs = applyImagePreconditioning(options, rhs, im, kk, ii);
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
    % if iscell(im_vectors.uCP)
        uPrev = im_vectors.uCP{timestep, ii};
        im_vectors.uCP{timestep, ii} = im_vectors.uCP{timestep, ii} - options.tauCP(ii) * rhs;
        if (options.enforcePositivity)
            im_vectors.uCP{timestep, ii}(im_vectors.uCP{timestep, ii} < options.epps) = options.epps;
        end
        output = im_vectors.uCP{timestep, ii} + options.thetaCP(kk) * (im_vectors.uCP{timestep, ii} - uPrev);
    % else
    %     uPrev = im_vectors.uCP;
    %     im_vectors.uCP = im_vectors.uCP - options.tauCP(ii) * rhs;
    %     if (options.enforcePositivity)
    %         im_vectors.uCP(im_vectors.uCP < options.epps) = options.epps;
    %     end
    %     output = im_vectors.uCP + options.thetaCP(kk) * (im_vectors.uCP - uPrev);
    % end
end
if (ii == 1 && options.PDAdaptiveType == 1)
    % if iscell(im_vectors.rhsCP)
        q = (im_old - im) / options.tauCP(ii) + options.subsets * im_vectors.rhsCP{timestep, ii};
    % else
    %     q = (im_old - im) / options.tauCP(ii) + options.subsets * im_vectors.rhsCP;
    % end
    w = dot((im_old - im), q) / (norm((im_old - im) * norm(q)));
    if (w < 0)
        options.tauCP(ii) = options.tauCP(ii) / (1 + options.alphaCP(ii));
        options.sigmaCP(ii) = options.sigmaCP(ii) * (1 + options.alphaCP(ii));
        options.alphaCP(ii) = options.alphaCP(ii) * 0.99;
    elseif (w >= .999)
        options.sigmaCP(ii) = options.sigmaCP(ii) / (1 + options.alphaCP(ii));
        options.tauCP(ii) = options.tauCP(ii) * (1 + options.alphaCP(ii));
        options.alphaCP(ii) = options.alphaCP(ii) * 0.99;
    end
end