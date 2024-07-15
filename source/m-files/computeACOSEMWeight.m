function ACOSEM_rhs = computeACOSEMWeight(options, recApu, mData, osa_iter, corrim_vectors, ii)
uu = sum(mData);
outputFP = forwardProject(options, recApu, osa_iter, corrim_vectors, ii);
if (options.param.CT)
    ACOSEM_rhs = sum(exp(-outputFP)) ./ uu;
else
    ACOSEM_rhs = uu ./ sum(outputFP);
end
end