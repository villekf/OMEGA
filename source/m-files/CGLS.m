function [im_vectors,options] = CGLS(options, iter, im_vectors, ii)
if (ii == options.nMultiVolumes + 1)
    if iscell(im_vectors.rhs)
        gamma_ = 0;
        for ll = 1 : options.nMultiVolumes + 1
            gamma_ = gamma_ + im_vectors.rhs{ll}' * im_vectors.rhs{ll};
        end
    else
        gamma_ = im_vectors.rhs' * im_vectors.rhs;
    end
    beta = gamma_ / options.gammaCGLS;
    if iscell(im_vectors.rhs)
        for ll = 1 : options.nMultiVolumes + 1
            im_vectors.fCGLS{ll} = im_vectors.fCGLS{ll} + options.alphaCGLS * im_vectors.recApu{ll};
            if (iter == options.Niter)
                im_vectors.recApu{ll} = im_vectors.fCGLS{ll};
            else
                im_vectors.recApu{ll} = im_vectors.rhs{ll} + beta * im_vectors.recApu{ll};
            end
        end
    else
        im_vectors.fCGLS = im_vectors.fCGLS + options.alphaCGLS * im_vectors.recApu;
        if (iter == options.Niter)
            im_vectors.recApu = im_vectors.fCGLS;
        else
            im_vectors.recApu = im_vectors.rhs + beta * im_vectors.recApu;
        end
    end
    options.gammaCGLS = gamma_;
end