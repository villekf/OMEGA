function offsetSqrt = resolveOffsetSqrt(options)
%RESOLVEOFFSETSQRT Resolves options.offsetSqrtWeights to a logical flag
%   Determines whether the square root of the offset redundancy weights is
%   applied on both the forward and the backward projection. This makes the
%   projection operator the linear sqrt(W) * A, which is needed by the
%   linearised least-squares algorithms (LSQR, CGLS, FISTA, FISTAL1, PDHG,
%   PDHGL1, PDDY, SART, ASD_POCS and BB). If options.offsetSqrtWeights is
%   empty (default), the weighting is enabled automatically when only
%   such algorithms are selected, if it is true the weighting is enabled
%   (also for custom algorithms) and if it is false the offset weights are
%   applied only in the backprojection. The weighting is available only
%   for CT data with offset correction, GPU implementations 2, 3 or 5,
%   backprojector type 4 or 5 and subsets formed from whole projections.
%   The native code (reconstructionAF.h) can still turn the weighting off
%   if the selected configuration is not supported.
%
% Input:
%   options = The options struct (obj.param of projectorClass)
% Output:
%   offsetSqrt = true if the sqrt(W) weighting is used

offsetSqrt = false;
if ~options.offsetCorrection || ~options.CT || options.FDK || options.largeDim || options.listmode
    return
end
if ~ismember(options.implementation, [2 3 5])
    return
end
% The backprojector (last digit of the projector type) has to be 4 or 5
if options.projector_type >= 10
    bpType = mod(options.projector_type, 10);
else
    bpType = options.projector_type;
end
if bpType ~= 4 && bpType ~= 5
    return
end
% The subsets must be formed from whole projections
if options.subsets > 1 && ~(options.subset_type == 0 || (options.subset_type >= 8 && options.subset_type <= 11))
    return
end
% Dynamic data are not supported
if numel(options.partitions) > 1 || options.partitions > 1
    return
end
if isfield(options, 'offsetSqrtWeights') && ~isempty(options.offsetSqrtWeights)
    offsetSqrt = logical(options.offsetSqrtWeights(1));
    return
end
% Automatic selection: only the linearised least-squares algorithms, no other algorithm
lsNames = setdiff(recNames(5), {'FDK'});
otherNames = setdiff(recNames(), recNames(5));
isLS = false;
for kk = 1 : numel(lsNames)
    isLS = isLS || (isfield(options, lsNames{kk}) && options.(lsNames{kk}) > 0);
end
isOther = false;
for kk = 1 : numel(otherNames)
    isOther = isOther || (isfield(options, otherNames{kk}) && options.(otherNames{kk}) > 0);
end
offsetSqrt = isLS && ~isOther;
end
