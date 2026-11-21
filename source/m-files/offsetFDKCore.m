function [ypad, xNew, info] = offsetFDKCore(y, x, u, pitch, deltaMm, pad)
%OFFSETFDKCORE Offset-detector FDK preprocessing of the line integrals
%   Pure array version of the redundancy weighting with conjugate
%   completion used by offsetFDKPreprocess. The redundancy weights are
%   applied before the ramp filtering and the missing side of the detector
%   is completed with virtual detector rows: a virtual pixel gets the value
%   of the measured conjugate ray (theta + pi, -rho) of the same detector
%   column (linear interpolation in the (theta, rho) space of the measured
%   rays on a Delaunay triangulation). The physical and the virtual rows
%   are multiplied with the weight w(s) (1 for s < -delta, 0 for s > delta,
%   0.5 (1 - sin(pi s / (2 delta))) in between, s being the signed distance
%   from the central ray in pixels, positive towards the missing side) and
%   with 2 (the weights of the conjugate rays sum to one, the full-scan
%   FDK expects two). The detector is then padded with zeros after the
%   virtual rows and its centre is shifted so that the physical rows keep
%   their exact world positions.
%
% Inputs:
%   y = Line integrals, size nRowsD x nColsD x nProjections (the detector
%   rows are the transaxial direction)
%   x = Source and detector centre coordinates, 6 x nProjections
%   [sx;sy;sz;dx;dy;dz]
%   u = The detector row direction vector scaled with the pixel size, 3 x
%   nProjections (or 2 x nProjections)
%   pitch = Detector pixel size in the row direction (mm)
%   deltaMm = Half-width of the transition region (mm, on the detector)
%   pad = Number of padded rows, empty for automatic (internal, testing)
% Outputs:
%   ypad = The weighted, completed and padded data (single), size
%   (nRowsD + P) x nColsD x nProjections. If the high-index side is
%   missing, the physical rows are the first ones, otherwise they start at
%   the row P + 1.
%   xNew = The coordinates x with the detector centre shifted (double)
%   info = Details (P, nVirtual, delta, ...)

if nargin < 6
    pad = [];
end
x = double(x);
u = double(u);
if size(u, 1) < 3
    u = [u; zeros(3 - size(u, 1), size(u, 2))];
end
N = size(y, 1);
K = size(y, 3);
if size(x, 2) ~= K || size(u, 2) ~= K
    error('offsetFDK:input', 'The number of the source/detector coordinates does not match the number of projections')
end
T = centralRay(x, u, N);
high = mean(T) > N / 2;
if high
    [ypad, P, info] = padHighSide(y, x, u, pitch, deltaMm, pad);
    uShift = u;
else
    % The low-index side is missing: mirror the detector (reverse the rows and the row direction), process it as the high-index side and mirror
    % back. The physical rows then start at the row P + 1 of the padded detector.
    uShift = -u;
    [ypadF, P, info] = padHighSide(flip(y, 1), x, uShift, pitch, deltaMm, pad);
    ypad = flip(ypadF, 1);
    clear ypadF
    info.T = N - info.T;
end
xNew = x;
xNew(4:6, :) = x(4:6, :) + (P / 2) * uShift;
info.high = high;
info.N = N;
end

function [ypad, P, info] = padHighSide(y, x, u, pitch, deltaMm, pad)
% Weighting, conjugate completion and padding for the missing side being the high-index side of the detector
edgeMargin = 2;     % pixels kept between the transition region and the far edge of the detector
windowMargin = 16;  % extra rows in the interpolation window
zeroTail = 256;     % zero rows beyond the virtual rows (the zero tail of the ramp filter)
holdFirst = 8;      % the hold value is the mean of the rows N - 8, ..., N - 4 (zero based)
holdLast = 4;
N = size(y, 1);
M = size(y, 2);
K = size(y, 3);
T = centralRay(x, u, N);
if any(T <= 0) || any(T >= N)
    error('offsetFDK:geometry', 'The central ray is outside the detector for some projections')
end
TT = N - T;
delta = deltaMm / pitch;
raised = false;
clamped = false;
if delta < max(TT)
    delta = max(TT);
    raised = true;
end
deltaMax = min(T) - edgeMargin;
if delta > deltaMax
    if deltaMax < max(TT)
        error('offsetFDK:geometry', 'The detector is too close to being centred for the offset FDK preprocessing')
    end
    delta = deltaMax;
    clamped = true;
end
nVirtual = ceil(delta - min(TT)) + 2;
if isempty(pad)
    P = min(N, nVirtual + zeroTail);
else
    P = round(pad);
end
if P < 1
    error('offsetFDK:input', 'The padding must be at least one row')
end
nVirtual = min(nVirtual, P);
nWindow = min(N, ceil(2 * delta) + windowMargin);
[S, outside] = interpolationMatrix(x, u, N, nVirtual, nWindow, pitch, K);
w = 2 * weightProfile(N + nVirtual, T, delta);
ypad = zeros(N + P, M, K, 'single');
ypad(1:N, :, :) = single(y) .* reshape(single(w(1:N, :)), N, 1, K);
% The hold value of the queries outside the hull: the mean of the last reliable physical rows
h0 = max(0, N - holdFirst);
holdRows = (h0 + 1) : max(h0 + 1, N - holdLast + 1);
wv = reshape(single(w(N + 1 : N + nVirtual, :)), nVirtual, 1, K);
chunk = 64;
for j0 = 1 : chunk : M
    jc = j0 : min(M, j0 + chunk - 1);
    mc = numel(jc);
    % Measured values of the interpolation window, row index k + K * r (projection k, window row r)
    vals = double(reshape(permute(y(N - nWindow + 1 : N, jc, :), [3 1 2]), K * nWindow, mc));
    out = reshape(S * vals, K, nVirtual, mc);
    holdVal = reshape(mean(single(y(holdRows, jc, :)), 1), mc, K)';
    mask3 = repmat(outside, 1, 1, mc);
    holdR = repmat(reshape(holdVal, K, 1, mc), 1, nVirtual, 1);
    out(mask3) = holdR(mask3);
    ypad(N + 1 : N + nVirtual, jc, :) = single(permute(out, [2 3 1])) .* wv;
end
info = struct('P', P, 'nVirtual', nVirtual, 'delta', delta, 'deltaMm', delta * pitch, 'raised', raised, 'clamped', clamped, ...
    'nWindow', nWindow, 'outsideFraction', mean(outside(:)), 'T', T);
end

function T = centralRay(x, u, nRows)
% Position of the central ray (the ray through the rotation axis) of every projection in pixels from the row-0 edge of the detector. The same
% closed form as OffsetLimit (divided by the pixel size) in setUpCorrections.
t = (x(2, :) .* x(4, :) - x(1, :) .* x(5, :)) ./ (x(1, :) .* u(2, :) - x(2, :) .* u(1, :));
T = t(:) + nRows / 2;
end

function w = weightProfile(nRows, T, delta)
% The redundancy weight at the pixel centres of the first nRows rows (before the factor of two), nRows x K. s = (i + 0.5) - T is the signed
% distance from the central ray in pixels, positive towards the missing side.
s = ((0 : nRows - 1)' + 0.5) - T(:)';
w = 0.5 * (1 - sin(pi * s / (2 * delta)));
w(s < -delta) = 1;
w(s > delta) = 0;
end

function [th, rho] = rayParameters(x, u, rows, nRows)
% The transaxial ray parameters of the rays through the centres of the detector rows (zero based indices, may be beyond the physical detector):
% the direction angle theta and the signed distance rho = s_x e_y - s_y e_x of the ray from the rotation axis. Both are K x numel(rows).
t = double(rows(:))' + 0.5 - nRows / 2;
qx = x(4, :)' + t .* u(1, :)';
qy = x(5, :)' + t .* u(2, :)';
ex = qx - x(1, :)';
ey = qy - x(2, :)';
nn = hypot(ex, ey);
ex = ex ./ nn;
ey = ey ./ nn;
th = atan2(ey, ex);
rho = x(1, :)' .* ey - x(2, :)' .* ex;
end

function [S, outside] = interpolationMatrix(x, u, N, nVirtual, nWindow, pitch, K)
% Sparse linear interpolation from the measured rays to the conjugate rays of the virtual pixels. The measured rays are those of the last
% nWindow rows of all projections, the queries are the conjugate rays (theta + pi, -rho) of the nVirtual virtual rows. The interpolation is the
% piecewise linear one on the Delaunay triangulation of the measured rays in the (theta / (2 pi / K), rho / (pitch SOD / SDD)) plane, theta being
% periodic. S is (K * nVirtual) x (K * nWindow), the row index is k + K * e (projection k, virtual row e) and the column index k + K * r
% (projection k, measured row N - nWindow + r). outside (K x nVirtual) marks the queries outside the hull of the measured rays (empty rows of S).
thetaWrap = 0.5;  % the part of the periodic theta axis that is duplicated over the seam (radians)
sod = mean(hypot(x(1, :), x(2, :)));
sdd = mean(sqrt(sum((x(1:3, :) - x(4:6, :)) .^ 2, 1)));
scTheta = 2 * pi / K;
scRho = pitch * sod / sdd;
[th, rho] = rayParameters(x, u, (N - nWindow) : (N - 1), N);
thF = th(:);
rhoF = rho(:);
nPts = numel(thF);
lo = find(thF < -pi + thetaWrap);
hi = find(thF > pi - thetaWrap);
src = [(1 : nPts)'; lo; hi];
pts = [[thF; thF(lo) + 2 * pi; thF(hi) - 2 * pi] / scTheta, [rhoF; rhoF(lo); rhoF(hi)] / scRho];
% Coincident points are not allowed in the triangulation (the first occurrence is kept)
[~, first] = unique(pts, 'rows', 'stable');
pts = pts(first, :);
src = src(first);
DT = delaunayTriangulation(pts(:, 1), pts(:, 2));
[thq, rhoq] = rayParameters(x, u, N : (N + nVirtual - 1), N);
thc = mod(thq, 2 * pi) - pi;
q = [thc(:) / scTheta, -rhoq(:) / scRho];
ti = pointLocation(DT, q);
ids = find(~isnan(ti));
bary = cartesianToBarycentric(DT, ti(ids), q(ids, :));
good = all(isfinite(bary), 2);
ids = ids(good);
bary = bary(good, :);
verts = DT.ConnectivityList(ti(ids), :);
inside = false(size(q, 1), 1);
inside(ids) = true;
S = sparse(repmat(ids, 1, 3), src(verts), bary, size(q, 1), nPts);
outside = reshape(~inside, K, nVirtual);
end
