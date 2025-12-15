function im_vectors = PDHG1(options, im_vectors, subIter, ii, timestep)
    if (options.PDAdaptiveType == 1)
        im_vectors.rhsCP{timestep, ii} = im_vectors.rhs{timestep, ii};
    end
    if (options.subsets > 1)
        if (options.verbose >= 3)
            disp("Using PDHG w/ subsets");
        end
        im_vectors.uCP{timestep, ii} = im_vectors.uCP{timestep, ii} + im_vectors.rhs{timestep, ii};
        im_vectors.rhs{timestep, ii} = im_vectors.uCP{timestep, ii} + (options.subsets * options.thetaCP(subIter)) .* im_vectors.rhs{timestep, ii};
    end
end