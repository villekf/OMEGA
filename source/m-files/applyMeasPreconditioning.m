function [input] = applyMeasPreconditioning(options, input, length, subIter)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if (options.precondTypeMeas(1) || options.precondTypeMeas(2))
    if (options.verbose >= 3)
        disp("Applying measurement-based preconditioning");
    end
    if (options.precondTypeMeas(2))
        if (options.verbose >= 3)
            disp("Applying filtering-based preconditioner, type 1");
        end
        if options.subset_type == 4
            input = reshape(input, options.nRowsD, length / options.nRowsD);
        else
            input = reshape(input, options.nRowsD, options.nColsD, length / (options.nRowsD * options.nColsD));
        end
        input = filtering(options.filter, input, options.Nf);
    end
    if (options.precondTypeMeas(1))
        if (options.verbose >= 3)
            disp("Applying diagonal normalization preconditioner (1 / (A1)), type 0");
        end
        input = input ./ options.M{subIter};
    end
    if (options.verbose >= 3)
        disp("Measurement-based preconditioning applied");
    end
end