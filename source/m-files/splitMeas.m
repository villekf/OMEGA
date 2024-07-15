function [splitMeas, varargout] = splitMeas(options, meas, pituus, ll, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if (options.subset_type >= 8 || options.subsets == 1) && ~options.listmode && options.projector_type ~= 6
    kerroin = options.nRowsD * options.nColsD;
else
    kerroin = 1;
end
if options.TOF
    meas = reshape(meas, numel(meas) / options.TOF_bins, options.TOF_bins);
    splitMeas = meas(pituus(ll) * kerroin + 1 : pituus(ll + 1) * kerroin, :);
    splitMeas = splitMeas(:);
else
    if options.projector_type == 6
        meas = reshape(meas, options.nRowsD, options.nColsD, []);
        splitMeas = meas(:,:,pituus(ll) * kerroin + 1 : pituus(ll + 1) * kerroin);
        splitMeas = splitMeas(:);
    else
        splitMeas = meas(pituus(ll) * kerroin + 1 : pituus(ll + 1) * kerroin);
    end
end
if options.randoms_correction && nargout >= 2
    varargout{1} = varargin{1}(pituus(ll) * kerroin + 1 : pituus(ll + 1) * kerroin);
elseif nargout >= 2
    varargout{1} = cast([], options.cType);
end