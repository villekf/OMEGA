function [input, options] = applyImagePreconditioning(options, input, im, kk, ii)
%APPLYIMAGEPRECONDITIONING Applies the selected image-based preconditioners
%   kk is the 1-based global subiteration number, (iter - 1) * subsets +
%   osa_iter, i.e. the C++ (0-based) subiteration counter plus one. The
%   gradient-based preconditioner (type 4) stores its weights into
%   options.gradF, so the returned options has to be kept by the caller
VAL = 0.00001;
if (options.precondTypeImage(5) && kk >= options.gradInitIter + 1)
    if (options.verbose >= 3)
        disp("Applying gradient-based preconditioner, type 4");
    end
    if (kk <= options.gradLastIter + 1)
        if ~isfield(options, 'gradF') || ~iscell(options.gradF)
            options.gradF = cell(options.nMultiVolumes + 1, 1);
        end
        if isfield(options, 'derivativeType')
            derivType = options.derivativeType;
        else
            derivType = 0;
        end
        % Same as gradientPreconditioner in functions.hpp
        imR = reshape(im, options.Nx(ii), options.Ny(ii), options.Nz(ii));
        f = zeros(size(imR), 'like', imR);
        g = f;
        h = f;
        [f, g, h] = computeGradient(imR, options, f, g, h, derivType);
        f = max(VAL, sqrt(f .* f + g .* g + h .* h) ./ mean(imR(:)));
        f = mean(f(:)) ./ f;
        options.gradF{ii} = min(options.gradV2, max(options.gradV1, f(:)));
    end
    if iscell(options.gradF)
        input = input .* options.gradF{ii};
    else
        input = input .* options.gradF;
    end
end
if (options.precondTypeImage(4))
    if (options.verbose >= 3)
        disp("Applying momentum-like preconditioner, type 3");
    end
    input = input .* options.alphaPrecond(kk);
end
if (options.precondTypeImage(1) || options.precondTypeImage(2) || options.precondTypeImage(3))
    if (options.precondTypeImage(1))
        if (options.verbose >= 3)
            disp("Applying diagonal normalization preconditioner , type 0");
        end
        if iscell(options.D)
            input = input ./ options.D{ii};
        else
            input = input ./ options.D;
        end
    elseif (options.precondTypeImage(2))
        if (options.verbose >= 3)
            disp("Applying EM preconditioner, type 1");
        end
        if iscell(options.D)
            input = input .* (im ./ options.D{ii});
        else
            input = input .* (im ./ options.D);
        end
    elseif (options.precondTypeImage(3))
        if (options.verbose >= 3)
            disp("Applying IEM preconditioner, type 2");
        end
        if ii == 1
            fieldname = 'referenceImage';
        else
            fieldname = ['referenceImage' num2str(ii)];
        end
        if iscell(options.D)
            input = input .* (max(im, max(VAL, options.(fieldname))) ./ options.D{ii});
        else
            input = input .* (max(im, max(VAL, options.(fieldname))) ./ options.D);
        end
    end
end
if (options.precondTypeImage(7))
    if (options.verbose >= 3)
        disp("Applying curvature preconditioner , type 6");
    end
    if iscell(options.dP)
        input = input .* options.dP{ii};
    else
        input = input .* options.dP;
    end
end
% C++: 0-based kk <= filterIter. computeOSEstimates disables the filter (precondTypeImage(6) = false) when kk - 1 == filteringIterations,
% before the update of that subiteration, so the filter is used for the first filteringIterations subiterations
if (options.precondTypeImage(6) && kk - 1 <= options.filteringIterations)
    if (options.verbose >= 3)
        disp("Applying filtering-based preconditioner, type 5");
    end
    input = reshape(input, options.Nx(ii), options.Ny(ii), options.Nz(ii));
    input = filtering2D(options.filterIm, input, options.Nf);
end
if (options.verbose >= 3 && (options.precondTypeImage(1) || options.precondTypeImage(2) || options.precondTypeImage(3) || options.precondTypeImage(4) || ...
        (options.precondTypeImage(5) && kk >= options.gradInitIter + 1) || options.precondTypeImage(6) || options.precondTypeImage(7)))
    disp("Image-based preconditioning applied");
end
end
