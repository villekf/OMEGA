function options = ParkerWeights(options)
% ParkerWeights
%
% Computes Parker weights as defined in DOI: 10.1118/1.1450132
% These are slightly different from the original Parker weights and allow
% for the optional "q" parameter, defined with options.ParkerWeight here.
%
% The weights are computed from the actual source and detector positions of
% every projection, i.e. the source polar angle, the rotation direction and
% the fan angle of each detector pixel (including detector offsets, mirrored
% detector axes and arbitrary start angles/arcs that cross 0 or 2*pi) are
% derived from options.x/options.z (set by setCTCoordinates). If the
% coordinates are not available, the built-in geometry is assumed and the
% projection angles are used.
%
% The angular quadrature is also included in the weights: the FDK
% backprojection integrates with the constant pi/nProjections (the full scan
% step 2*pi/nProjections times the 1/2 of the 360 degree redundancy), and the
% Parker weights of two conjugate rays sum to 1 instead of the 1/2 + 1/2 of a
% full scan. Each weight is therefore multiplied by dBeta_k * nProjections / pi,
% where dBeta_k is the actual angular step of projection k. For a full 2*pi
% scan dBeta_k = 2*pi/nProjections and the factor is 2.
%
% INPUT:
%   options.angles            [options.nProjections x 1] or [1 x options.nProjections], projection angles [rad]
%   options.sourceToDetector  scalar, source-to-detector distance [mm]
%   options.nRowsD            scalar, detector size in row direction [number of pixels]
%   options.dPitchX           scalar, detector pixel pitch in row direction [mm]
%
% OPTIONAL:
%   options.ParkerWeight      scalar in (0,1], transition fraction, default 0.25. This is the q parameter from the article
%   options.detOffsetRow      scalar or one value per projection, detector center offset [same units as dPitchX], default 0 [mm].
%                             Only used if the source/detector coordinates are not available
%
% OUTPUT:
%   options.SinM weighted with the Parker weights
%

% Required inputs
betaIn = options.angles(:);

DSD   = options.sourceToDetector;
nU    = options.nRowsD;
du    = options.dPitchX;
nProj = numel(betaIn);

if isfield(options,'ParkerWeight') && options.ParkerWeight > 0
    q = options.ParkerWeight;
else
    q = 0.25;
end

if ~(isscalar(DSD) && isscalar(nU) && isscalar(du))
    error('sourceToDetector, nRowsD, and dPitchX must be scalars.');
end
if q <= 0 || q > 1
    error('options.ParkerWeight must satisfy 0 < q <= 1.');
end
if nProj < 2
    error('options.angles must contain at least two projection angles.');
end

[hasGeom, sx, sy, dx, dy, ux, uy] = projectionGeometry(options, nProj);

if hasGeom
    % Source polar angle (unwrapped), rotation direction and fan angle of every pixel from the actual geometry
    beta = unwrap(atan2(sy, sx));
    ccw = beta(end) >= beta(1);
    sad = hypot(sx, sy);
    % Unit vector from the source towards the axis of rotation
    cx = -sx ./ sad;
    cy = -sy ./ sad;
    j = (0:nU-1) - (nU-1)/2;
    rx = dx + j * du .* ux - sx;
    ry = dy + j * du .* uy - sy;
    % Fan angle [nProj x nU], positive counter-clockwise from the central ray
    alpha = atan2(-rx .* cy + ry .* cx, rx .* cx + ry .* cy);
    if ~ccw
        % Clockwise scan: mirror both the source angle and the fan angle
        alpha = -alpha;
    end
else
    % Built-in geometry from the angles. Detector offset: any sign, scalar or one value per projection
    beta = unwrap(betaIn);
    ccw = beta(end) >= beta(1);
    detOffset = 0;
    if isfield(options,'detOffsetRow') && ~isempty(options.detOffsetRow) && any(options.detOffsetRow(:) ~= 0)
        off = double(options.detOffsetRow(:));
        if isscalar(off)
            detOffset = off;
        else
            detOffset = repmat(off, ceil(nProj / numel(off)), 1);
            detOffset = detOffset(1:nProj);
        end
    end
    % Flat-panel detector coordinate u [same units as DSD] (the built-in detector u axis points in the direction of the
    % counter-clockwise fan angle)
    u = ((0:nU-1) - (nU-1)/2) * du + detOffset;
    alpha = atan2(u .* ones(nProj, 1), DSD);
    if ~ccw
        alpha = -alpha;
    end
end

% Use actual maximum fan half-angle from detector edges
delta = max(abs(alpha(:)));

% Relative source angle along the scan
betaRel = abs(beta - beta(1));
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

% Angular quadrature, see above
dBeta = abs(gradient(beta));
angFactor = dBeta * nProj / pi;

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
    g = alpha(:, iu);
    B = epsilon + 2*delta - 2*g;
    b = q * B;
    B2 = epsilon + 2*delta + 2*g;
    b2 = q * B2;

    if any(b <= 0) || any(b2 <= 0)
        error('Encountered nonpositive transition length b. Check geometry.');
    end

    x1 =  betaRel ./ b - 0.5;
    x2 = (betaRel - B) ./ b + 0.5;
    x3 = (betaRel - pi + 2*g) ./ b2 - 0.5;
    x4 = (betaRel - pi - 2*delta - epsilon) ./ b2 + 0.5;

    w2 = 0.5 * ( S(x1) + S(x2) - S(x3) - S(x4) );
    w2 = max(0, min(1, w2)) .* angFactor;
    options.SinM(iu, :, :) = options.SinM(iu, :, :) .* permute(w2, [3 2 1]);
end

end

function [hasGeom, sx, sy, dx, dy, ux, uy] = projectionGeometry(options, nProj)
% In-plane (x/y) geometry of every projection as set by setCTCoordinates: the source position, the detector centre and the
% unit vector of the detector row (u) direction (column vectors). hasGeom is false if the coordinates are not available in the
% expected form (for example helical CT), in which case only the projection angles can be used.
hasGeom = false;
sx = []; sy = []; dx = []; dy = []; ux = []; uy = [];
if ~isfield(options,'x') || ~isfield(options,'z') || isempty(options.x) || isempty(options.z)
    return;
end
x = double(options.x);
z = double(options.z);
% setCTCoordinates stores the coordinates as [6 x nProj] (x/y/z of the source, x/y/z of the detector) and the direction vectors
% as [2 or 6 x nProj]
if size(x, 1) ~= 6 && size(x, 2) == 6
    x = x';
end
if size(z, 1) ~= 2 && size(z, 1) ~= 6 && (size(z, 2) == 2 || size(z, 2) == 6)
    z = z';
end
if size(x, 1) ~= 6 || size(x, 2) < nProj || size(z, 2) < nProj || size(z, 1) < 2
    return;
end
sx = x(1, 1:nProj)';
sy = x(2, 1:nProj)';
dx = x(4, 1:nProj)';
dy = x(5, 1:nProj)';
un = hypot(z(1, 1:nProj)', z(2, 1:nProj)');
if any(un <= 0) || any(hypot(sx, sy) <= 0)
    return;
end
ux = z(1, 1:nProj)' ./ un;
uy = z(2, 1:nProj)' ./ un;
hasGeom = true;
end
