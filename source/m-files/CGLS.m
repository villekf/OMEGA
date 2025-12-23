function [im_vectors,options] = CGLS(options, iter, im_vectors, ii, timestep)
if (ii == options.nMultiVolumes + 1)
    gamma_ = 0;
    for ll = 1 : options.nMultiVolumes + 1
        gamma_ = gamma_ + im_vectors.rhs{timestep, ll}' * im_vectors.rhs{timestep, ll};
    end
    beta = gamma_ / options.gammaCGLS;
    for ll = 1 : options.nMultiVolumes + 1
        im_vectors.fCGLS{timestep, ll} = im_vectors.fCGLS{timestep, ll} + options.alphaCGLS * im_vectors.recApu{timestep, ll};
        if (iter == options.Niter)
            im_vectors.recApu{timestep, ll} = im_vectors.fCGLS{timestep, ll};
        else
            im_vectors.recApu{timestep, ll} = im_vectors.rhs{timestep, ll} + beta * im_vectors.recApu{timestep, ll};
        end
    end
    options.gammaCGLS = gamma_;
end