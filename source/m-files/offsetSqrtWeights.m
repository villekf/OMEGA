function w = offsetSqrtWeights(T, nRowsD, dPitchX)
%OFFSETSQRTWEIGHTS Square root of the offset (redundancy) weights
%   Computes the square root of the offset weights at the detector pixel
%   centers. These match the weights that the projector type 4 and 5
%   backprojection kernels apply when OFFSETSQRT is defined, i.e.
%   sqrt(.5 * (1 + sin(pi * px / (2 * T)))).
%
% Inputs:
%   T = The distance (mm) of the central ray from the row-0 edge of the
%   detector for each projection (options.OffsetLimit)
%   nRowsD = The number of detector rows (pixels in the offset direction)
%   dPitchX = The detector pixel size (mm) in the offset direction
% Outputs:
%   w = The weights, size nRowsD x number of projections (single). In the
%   vectorized measurement data the detector row is the fastest changing
%   index, i.e. w(:,kk) corresponds to the rows of projection kk

T = double(T(:))';
nRowsD = double(nRowsD);
koko = nRowsD * double(dPitchX);
% Normalized detector coordinate of the pixel centers
px = repmat(((0 : nRowsD - 1)' + 0.5) / nRowsD, 1, numel(T));
% Same flipping as in the kernels when the central ray is on the far half of the detector
flip = T > koko / 2;
px(:, flip) = abs(px(:, flip) - 1);
TT = T;
TT(flip) = koko - T(flip);
TT = repmat(TT, nRowsD, 1);
p = px * koko - TT;
w = sin(pi * (p + TT) ./ (4 * TT));
w(p > TT) = 1;
w(p < -TT) = 0;
w = single(w);
end
