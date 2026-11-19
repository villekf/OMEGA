function ACOSEM_rhs = computeACOSEMWeight(options, recApu, mData, osa_iter, corrim_vectors, ii, varargin)
%COMPUTEACOSEMWEIGHT Scaling weight of ACOSEM and OSL-COSEM type 1 (computeACOSEMWeight in functions.hpp)
%   The optional seventh input is the sum of the additive randoms/scatter
%   term of the current subset (already multiplied by the number of TOF
%   bins when the same randoms are used for every TOF bin). The forward
%   projection does not include the additive term, so it is removed from
%   the data sum, uu = sum(y) - sum(r) (PET/SPECT only)
uu = sum(mData);
if ~options.param.CT && ~isempty(varargin) && ~isempty(varargin{1})
    uu = uu - varargin{1};
end
outputFP = forwardProject(options, recApu, osa_iter, corrim_vectors, ii);
if (options.param.CT)
    ACOSEM_rhs = sum(exp(-outputFP)) ./ uu;
else
    ACOSEM_rhs = uu ./ sum(outputFP);
end
end
