function W = arcInterpWeights(x_o, y_o, x, y, Ndist, Nang, method)
%ARCINTERPWEIGHTS Sparse interpolation weights between two LOR grids
%   W = ARCINTERPWEIGHTS(X_O, Y_O, X, Y, NDIST, NANG, METHOD) builds a
%   sparse matrix W such that NEW_VALUES = W * OLD_VALUES maps sinogram
%   values living on the original (non-arc-corrected) LOR grid, given by
%   the LOR endpoint coordinates X_O/Y_O (N x 2 each, N = NDIST*NANG), onto
%   the arc-corrected LOR grid given by the endpoint coordinates X/Y (N x 2
%   each). Both LOR sets are first converted to a canonical (angle, signed
%   distance) parametrisation so that the interpolation happens in an
%   absolute frame shared by both grids (no rotation offset between the
%   interpolated sinogram and the returned coordinates).
%
%   METHOD is 'linear' (barycentric interpolation over a Delaunay
%   triangulation, default) or 'nearest' (nearest-neighbor). For any other
%   method W is returned empty ([]) so that the caller can fall back to a
%   slower, generic interpolation (see arcCorrection).
%
%   If delaunayTriangulation does not exist (e.g. GNU Octave), delaunay,
%   tsearch and dsearchn are used instead. The fallback triangulation may
%   split cocircular point quadruples differently, which changes the linear
%   weights slightly inside such quadrilaterals.
%
% See also arcCorrection

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2026 Ville-Veikko Wettenhovi
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 7 || isempty(method)
    method = 'linear';
end
if ~(strcmpi(method, 'linear') || strcmpi(method, 'nearest'))
    W = [];
    return
end

N = Ndist * Nang;

% Canonical (angle, signed distance) parameters, same function for both
% the original and the arc-corrected LOR set
[th_o, s_o] = lorParams(x_o, y_o);
[th_n, s_n] = lorParams(x, y);

% Units: angular bins and distance bins (ds always from the original set)
ds = (max(s_o) - min(s_o)) / (Ndist - 1);
u_o = th_o / (pi / Nang);
v_o = s_o / ds;
u_n = th_n / (pi / Nang);
v_n = s_n / ds;

% Periodic extension of the original points: shifting the angle by +-Nang
% bins (i.e. by pi) represents the same physical line with the opposite
% signed distance
P = [u_o, v_o; u_o - Nang, -v_o; u_o + Nang, -v_o];
srcIdx = repmat((1 : N)', 3, 1);

% Merge duplicate points (within numerical tolerance). S (numUnique x N)
% maps the original values to the (averaged) unique-point values:
% unique_values = S * old_values
Pr = round(P * 1e6) / 1e6;
[uniqueP, ~, ic] = unique(Pr, 'rows', 'stable');
numUnique = size(uniqueP, 1);
S = sparse(ic, srcIdx, 1, numUnique, N);
rowSums = sum(S, 2);
rowSums(rowSums == 0) = 1;
S = spdiags(1 ./ rowSums, 0, numUnique, numUnique) * S;

Q = [u_n, v_n];
Nq = size(Q, 1);

useDT = exist('delaunayTriangulation') ~= 0;
if useDT
    DT = delaunayTriangulation(uniqueP(:,1), uniqueP(:,2));
    tri = DT.ConnectivityList;
else
    tri = delaunay(uniqueP(:,1), uniqueP(:,2));
end

if strcmpi(method, 'nearest')
    if useDT
        nearestIdx = nearestNeighbor(DT, Q);
    else
        nearestIdx = dsearchn(uniqueP, Q);
    end
    T = sparse(1 : Nq, nearestIdx, 1, Nq, numUnique);
else
    if useDT
        [ti, bc] = pointLocation(DT, Q);
    else
        ti = tsearch(uniqueP(:,1), uniqueP(:,2), tri, Q(:,1), Q(:,2));
        bc = zeros(Nq, 3);
        in = ~isnan(ti);
        if any(in)
            tv = tri(ti(in), :);
            Qi = Q(in, :);
            x1 = uniqueP(tv(:,1), 1); y1 = uniqueP(tv(:,1), 2);
            x2 = uniqueP(tv(:,2), 1); y2 = uniqueP(tv(:,2), 2);
            x3 = uniqueP(tv(:,3), 1); y3 = uniqueP(tv(:,3), 2);
            detT = (y2 - y3) .* (x1 - x3) + (x3 - x2) .* (y1 - y3);
            l1 = ((y2 - y3) .* (Qi(:,1) - x3) + (x3 - x2) .* (Qi(:,2) - y3)) ./ detT;
            l2 = ((y3 - y1) .* (Qi(:,1) - x3) + (x1 - x3) .* (Qi(:,2) - y3)) ./ detT;
            bc(in, :) = [l1, l2, 1 - l1 - l2];
        end
    end
    inHull = ~isnan(ti);
    rows = zeros(0, 1);
    cols = zeros(0, 1);
    vals = zeros(0, 1);
    if any(inHull)
        verts = tri(ti(inHull), :);
        w = bc(inHull, :);
        w(w < 0 & w > -1e-12) = 0;
        qIdx = find(inHull);
        rows = [rows; repmat(qIdx, 3, 1)];
        cols = [cols; verts(:)];
        vals = [vals; w(:)];
    end
    if any(~inHull)
        outIdx = find(~inHull);
        if useDT
            nearestIdx = nearestNeighbor(DT, Q(outIdx, :));
        else
            nearestIdx = dsearchn(uniqueP, Q(outIdx, :));
        end
        rows = [rows; outIdx];
        cols = [cols; nearestIdx];
        vals = [vals; ones(numel(outIdx), 1)];
    end
    T = sparse(rows, cols, vals, Nq, numUnique);
end

W = T * S;
end

function [th, s] = lorParams(x, y)
% Canonical (angle, signed distance) parametrisation of a set of LORs,
% each given by two endpoints (columns 1 and 2 of x and y). The angle is
% folded into [0, pi) with the sign of the distance flipped accordingly,
% so that a line and its reverse-order equivalent map to the same point.
x1 = x(:,1); x2 = x(:,2);
y1 = y(:,1); y2 = y(:,2);
dx = x2 - x1;
dy = y2 - y1;
th = atan2(dy, dx);
s = (x1 .* y2 - x2 .* y1) ./ hypot(dx, dy);
idx = th < 0;
th(idx) = th(idx) + pi;
s(idx) = -s(idx);
idx = th >= pi;
th(idx) = th(idx) - pi;
s(idx) = -s(idx);
end
