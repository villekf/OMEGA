function options = ParkerWeights(options)
% ParkerWeights
%
% Computes Parker weights as defined in DOI: 10.1118/1.1450132
% These are slightly different from the original Parker weights and allow
% for the optional "q" parameter, defined with options.ParkerWeight here.
%
% INPUT:
%   options.angles            [options.nProjections x 1] or [1 x options.nProjections], projection angles [rad]
%   options.sourceToDetector  scalar, source-to-detector distance [mm]
%   options.nRowsD            scalar, detector size in row direction [number of pixels]
%   options.dPitchX           scalar, detector pixel pitch in row direction [mm]
%
% OPTIONAL:
%   options.ParkerWeight      scalar in (0,1], transition fraction, default 0.25. This is the q parameter from the article
%   options.detOffsetRow      scalar, row detector center offset [same units as dPitchX], default 0 [mm]
%
% OUTPUT:
%   w    [nRowsD x nColsD x nProj] or [nRowsD x nProj] weights
%

% Required inputs
betaIn = options.angles(:);

DSD   = options.sourceToDetector;
nU    = options.nRowsD;
du    = options.dPitchX;

if isfield(options,'ParkerWeight') && options.ParkerWeight > 0
    q = options.ParkerWeight;
else
    q = 0.25;
end

if isfield(options,'detOffsetRow ') && sum(options.detOffsetRow) > 0
    detOffset = options.detOffsetRow;
else
    detOffset = 0;
end

if ~(isscalar(DSD) && isscalar(nU) && isscalar(du))
    error('sourceToDetector, nRowsD, and dPitchX must be scalars.');
end
if q <= 0 || q > 1
    error('options.ParkerWeight must satisfy 0 < q <= 1.');
end

% Flat-panel detector coordinate u [same units as DSD]
u = ((0:nU-1) - (nU-1)/2) * du + detOffset;

% Horizontal fan angle alpha for each detector pixel
alpha = atan2(u, DSD);

% Use actual maximum fan half-angle from detector edges
delta = max(abs(alpha));

if numel(betaIn) < 2
    error('options.angles must contain at least two projection angles.');
end

betaRel = betaIn - betaIn(1);
if betaRel(end) < 0
    betaRel = abs(betaRel);
    alpha = -alpha;
end

scanRange = betaRel(end);

% Short-scan validity
minShortScan = pi + 2*delta;
if scanRange + 1e-8 < minShortScan
    error(['Angular range is too short. Need at least pi + 2*delta = %.6f rad, ', ...
        'got %.6f rad.'], minShortScan, scanRange);
end
if scanRange > 2*pi + 1e-8
    warning('Angular range exceeds 2*pi. Weights are intended for short/over-scan up to 2*pi.');
end

epsilon = max(scanRange - (pi + 2*delta), 0);

% S-function
    function y = S(x)
        y = zeros(size(x));
        m1 = (x > -0.5) & (x < 0.5);
        m2 = (x >= 0.5);
        y(m1) = 0.5 * (1 + sin(pi * x(m1)));
        y(m2) = 1;
    end

% Compute weights

for iu = 1:nU
    g = alpha(iu);
    B = epsilon + 2*delta - 2*g;
    b = q * B;
    B2 = epsilon + 2*delta + 2*g;
    b2 = q * B2;

    if b <= 0 || b2 <= 0
        error('Encountered nonpositive transition length b. Check geometry.');
    end

    x1 =  betaRel / b - 0.5;
    x2 = (betaRel - B) / b + 0.5;
    x3 = (betaRel - pi + 2*g) / b2 - 0.5;
    x4 = (betaRel - pi - 2*delta - epsilon) / b2 + 0.5;

    w2 = 0.5 * ( S(x1) + S(x2) - S(x3) - S(x4) );
    w2 = max(0, min(1, w2));
    options.SinM(iu, :, :) = options.SinM(iu, :, :) .* permute(w2, [3 2 1]);
end

end