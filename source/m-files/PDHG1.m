function im_vectors = PDHG1(options, im_vectors, subIter, ii)
if (options.PDAdaptiveType == 1)
    if iscell(im_vectors.rhsCP)
        im_vectors.rhsCP{ii} = im_vectors.rhs{ii};
    else
        im_vectors.rhsCP = im_vectors.rhs;
    end
end
if (options.subsets > 1)
    if (options.verbose >= 3)
        disp("Using PDHG w/ subsets");
    end
    if iscell(im_vectors.uCP)
        im_vectors.uCP{ii} = im_vectors.uCP{ii} + im_vectors.rhs{ii};
        im_vectors.rhs{ii} = im_vectors.uCP{ii} + (options.subsets * options.thetaCP(subIter)) .* im_vectors.rhs{ii};
    else
        im_vectors.uCP = im_vectors.uCP + im_vectors.rhs;
        im_vectors.rhs = im_vectors.uCP + (options.subsets * options.thetaCP(subIter)) .* im_vectors.rhs;
    end
end
end