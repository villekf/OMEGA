function [input] = applyImagePreconditioning(options, input, im, kk, ii)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
VAL = 0.00001;
if (options.precondTypeImage(5) && kk >= options.gradInitIter + 1)
    if (options.verbose >= 3)
        disp("Applying gradient-based preconditioner, type 4");
    end
    if (kk <= options.gradLastIter + 1)
        options = gradientPreconditioner(options, reshape(im, options.Nx(ii), options.Ny(ii), options.Nz(ii)), ii, VAL);
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
            input = input * (max(im, max(VAL, options.(fieldname))) ./ options.D{ii});
        else
            input = input * (max(im, max(VAL, options.(fieldname))) ./ options.D);
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
if (options.precondTypeImage(6) && kk <= options.filteringIterations)
    if (options.verbose >= 3)
        disp("Applying filtering-based preconditioner, type 5");
    end
    input = reshape(input, options.Nx(ii), options.Ny(ii), options.Nz(ii));
    input = filtering2D(options.filterIm, input, options.Nf);
end
if (options.verbose >= 3 && (options.precondTypeImage(1) || options.precondTypeImage(2) || options.precondTypeImage(3) || options.precondTypeImage(4) || ...
        (options.precondTypeImage(5) && kk >= options.gradInitIter) || options.precondTypeImage(6) || options.precondTypeImage(7)))
    disp("Image-based preconditioning applied");
end
end