function [x, y, options] = arcCorrection(options, interpolateSinogram)
%ARCCORRECTION Performs arc correction on the detector coordinates
%   This function projects the (origin-centred) detector coordinates onto
%   an ideal circle of radius (diameter + 2*DOI)/2 to obtain an
%   equidistant, arc-corrected LOR geometry, and optionally interpolates
%   the measurement data (and, when applicable, the normalization,
%   randoms/scatter, and scatter correction data used during
%   reconstruction) from the original LOR grid onto this new grid. The
%   interpolation is performed with sparse barycentric (Delaunay) or
%   nearest-neighbor weights (see arcInterpWeights), computed in an
%   absolute angle/signed-distance frame shared by both LOR grids, so
%   that the interpolated data are not rotated relative to the returned
%   coordinates. Works only with sinogram data and without precomputation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020-2026 Ville-Veikko Wettenhovi
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

if isfield(options, 'DOI') && ~isempty(options.DOI)
    DOI = options.DOI;
else
    DOI = 0;
end
% The detectors reside at radius diameter/2 + DOI, not diameter/2
D = options.diameter + 2 * DOI;
cpb = options.cryst_per_block(1);

% The original (non-arc-corrected) LOR coordinates, used only for the
% interpolation, are computed with the TRUE options (flip as given)
[~, ~, orig_xp, orig_yp] = detector_coordinates(options);
[x_o, y_o] = sinogram_coordinates_2D(options, orig_xp, orig_yp);

% The arc-corrected geometry itself (new_xp/new_yp, J, eka/vika,
% rotations) is built with flip_image forced to false: flip_image breaks
% the quadrant mirror-fill below, and flipping is not needed here since
% the interpolation (below) works with absolute positions
geomOptions = options;
geomOptions.flip_image = false;
[~, ~, xp, yp] = detector_coordinates(geomOptions);

new_xp = zeros(size(xp));
new_yp = zeros(size(yp));
% The angles of the blocks/buckets
if mod(options.blocks_per_ring, 4) == 0
    l_angles = linspace(0,90, options.blocks_per_ring / 4 + 1);
else
    l_angles = linspace(0,180, options.blocks_per_ring / 2 + 1);
    l_angles = l_angles(1 : ceil(options.blocks_per_ring / 4));
end
ll = 1;

% Shift to the zero angle (with x-axis)
if options.det_w_pseudo > options.det_per_ring
    shift = round(options.offangle) + floor((cpb + 1) / 2);
else
    shift = round(options.offangle) + floor(cpb / 2);
end
xp = circshift(xp, -shift);
yp = circshift(yp, -shift);

% Determine the points along the circle that reside (approximately) on the
% same line as the original LOR
for kk = 1 : length(xp) / 4
    if options.det_w_pseudo > options.det_per_ring
        vali = kk + options.det_w_pseudo / 2 - (cpb + 1) : kk + options.det_w_pseudo / 2 + (cpb + 1);
    else
        vali = kk + options.det_w_pseudo / 2 - cpb : kk + options.det_w_pseudo / 2 + cpb;
    end
    kulma = l_angles(ll);
    angle1 = atand((yp(kk)-yp(vali))./(xp(kk)-xp(vali)));
    angle1 = round(angle1*10^4)/10^4;
    if mod(options.blocks_per_ring, 4) == 0
        ind1 = find(abs(angle1) == kulma);
        if isempty(ind1)
            [~, ind1] = min(abs(abs(angle1) - kulma));
        end
    else
        [~, ind1] = min(abs(abs(angle1) - kulma));
    end
    if options.det_w_pseudo > options.det_per_ring
        x2 = xp(kk + options.det_w_pseudo / 2 - (cpb + 1) + ind1 - 1);
        y2 = yp(kk + options.det_w_pseudo / 2 - (cpb + 1) + ind1 - 1);
    else
        x2 = xp(kk + options.det_w_pseudo / 2 - cpb + ind1 - 1);
        y2 = yp(kk + options.det_w_pseudo / 2 - cpb + ind1 - 1);
    end
    p = [xp(kk); yp(kk)];
    q = [x2; y2];
    d = q - p;
    l = ((-dot(2*p,d) - sqrt(dot(2*p,d)^2 - 4*dot(d,d)*(dot(p,p) - (D/2)^2)))/ (2*dot(d,d)));
    lx = p + l * d;
    new_xp(kk) = lx(1);
    new_yp(kk) = lx(2);
end

new_xp(1 : length(xp) / 4) = new_xp(1 : length(xp) / 4) + D/2;
new_yp(1 : length(yp) / 4) = new_yp(1 : length(yp) / 4) + D/2;

new_yp(kk + 1: kk * 2) = flip(new_yp(1: kk));
diffi = diff([flip(new_yp(1 : kk)) ; D/2]);
diffi = D/2 + cumsum(flip(diffi));
new_yp(kk * 2 + 1 : end) = [diffi ; flip(diffi)];

diffi = diff([new_xp(1 : kk) ; D/2]);
diffi = D/2 + cumsum(flip(diffi));
new_xp(kk + 1: kk * 2) = diffi;
new_xp(kk * 2 + 1 : end) = flip(new_xp(1 : kk * 2));
xp = new_xp;
yp = new_yp;

xp = circshift(xp, shift);
yp = circshift(yp, shift);

[x, y] = sinogram_coordinates_2D(geomOptions, xp, yp);


xx1 = reshape(x(:,1),options.Ndist,options.Nang);
xx2 = reshape(x(:,2),options.Ndist,options.Nang);
yy1 = reshape(y(:,1),options.Ndist,options.Nang);
yy2 = reshape(y(:,2),options.Ndist,options.Nang);

angle = atand((yy1-yy2)./(xx1-xx2)) + 90;
angle(angle == 180) = 0;
[~, J] = min(mean(angle));

% Create the arc corrected coordinates from the perpendicular LOR
% coordinates
eka = xx1(1,J);
vika = xx1(end,J);
alkux = linspace(eka, vika, options.Ndist);

% Determine the y-coordinates
alkuy = sqrt((D/2)^2 - (alkux - D/2).^2) + D/2;

alku = [alkux; alkuy];

alku2 = alku - D/2;

alku = [alkux; abs(alkuy - D)];

alku1 = alku - D/2;


angles = linspace(0, 180, options.Nang + 1);
angles = angles(1 : end - 1);
angles = circshift(angles, J);
angles = reshape(angles, 1, 1, []);

% Rotation matrix
rot_matrix = [cosd(angles) -sind(angles); sind(angles) cosd(angles)];

rot_matrix = squeeze(num2cell(rot_matrix, [1 2]))';

% Arc correction
% Rotate the arc corrected coordinates for the whole ring
new_xy1 = cell2mat(cellfun(@(x) x * alku1, rot_matrix, 'UniformOutput', false))' + D/2;
new_xy2 = cell2mat(cellfun(@(x) x * alku2, rot_matrix, 'UniformOutput', false))' + D/2;


x = [new_xy1(:,1), new_xy2(:,1)];
y = [new_xy1(:,2), new_xy2(:,2)];

xx1 = reshape(x(:,1),options.Ndist,options.Nang);
xx2 = reshape(x(:,2),options.Ndist,options.Nang);
apu = xx1(:, 1: J);
xx1(:, 1: J) = flipud(xx2(:, 1: J));
xx2(:, 1: J) = flipud(apu);
yy1 = reshape(y(:,1),options.Ndist,options.Nang);
yy2 = reshape(y(:,2),options.Ndist,options.Nang);
apu = yy1(:, 1: J);
yy1(:, 1: J) = flipud(yy2(:, 1: J));
yy2(:, 1: J) = flipud(apu);

x(:,1) = xx1(:);
x(:,2) = xx2(:);
y(:,1) = yy1(:);
y(:,2) = yy2(:);

x = x - D / 2;
y = y - D / 2;

% Equidistance check: the signed distance of each LOR from the origin
% must change by a constant step within each angular bin (column). This
% uses the UNFOLDED signed distance (not lorDistance, which folds the
% angle into [0,pi) and flips the sign of s along with it): a column
% whose LOR angle lies exactly at the 0/pi fold would otherwise have its
% sign flip mid-column and trigger a false equidistance warning.
% Endpoint order (columns 1/2 of x,y) is consistent within each column.
s_chk = (x(:,1) .* y(:,2) - x(:,2) .* y(:,1)) ./ hypot(x(:,2) - x(:,1), y(:,2) - y(:,1));
s_chk = reshape(s_chk, options.Ndist, options.Nang);
d_chk = diff(s_chk, 1, 1);
colSpread = max(d_chk, [], 1) - min(d_chk, [], 1);
if max(colSpread) > 1e-3 * median(abs(d_chk(:)))
    warning('Arc correction failed to make all the LORs equidistant.')
end

if interpolateSinogram
    tic
    W = arcInterpWeights(x_o, y_o, x, y, options.Ndist, options.Nang, options.arc_interpolation);

    N = options.Ndist * options.Nang;
    applyToField = @(A) applyArc(A, W, options.Ndist, options.Nang, options.arc_interpolation, x_o, y_o, x, y);

    % Interpolate the measured sinogram(s) (numeric array, or a cell
    % array of partitions)
    if iscell(options.SinM)
        for kk = 1 : options.partitions
            options.SinM{kk} = applyToField(options.SinM{kk});
        end
    else
        options.SinM = applyToField(options.SinM);
    end

    % When corrections are applied during reconstruction, the
    % normalization / randoms (SinDelayed) / scatter correction data are
    % still in the original (non-arc-corrected) geometry and must be
    % interpolated with the same weights. If normalization has already
    % been precorrected into SinM (not during reconstruction), it must
    % not be touched here.
    if options.normalization_correction && options.corrections_during_reconstruction ...
            && isfield(options, 'normalization') && numel(options.normalization) > 1 ...
            && mod(numel(options.normalization), N) == 0
        options.normalization = applyToField(options.normalization);
    end

    if isfield(options, 'SinDelayed')
        if iscell(options.SinDelayed)
            for kk = 1 : numel(options.SinDelayed)
                if numel(options.SinDelayed{kk}) > 1 && mod(numel(options.SinDelayed{kk}), N) == 0
                    options.SinDelayed{kk} = applyToField(options.SinDelayed{kk});
                end
            end
        else
            if numel(options.SinDelayed) > 1 && mod(numel(options.SinDelayed), N) == 0
                options.SinDelayed = applyToField(options.SinDelayed);
            end
        end
    end

    if isfield(options, 'ScatterC') && options.scatter_correction
        if iscell(options.ScatterC)
            for kk = 1 : numel(options.ScatterC)
                if numel(options.ScatterC{kk}) > 1 && mod(numel(options.ScatterC{kk}), N) == 0
                    options.ScatterC{kk} = applyToField(options.ScatterC{kk});
                end
            end
        else
            if numel(options.ScatterC) > 1 && mod(numel(options.ScatterC), N) == 0
                options.ScatterC = applyToField(options.ScatterC);
            end
        end
    end
    endTime = toc;
    if options.verbose
        disp(['Arc correction complete in ' num2str(endTime) ' seconds'])
    end
end
end

function out = applyArc(A, W, Ndist, Nang, method, x_o, y_o, x, y)
% Applies the arc-correction interpolation to A if (and only if) A is
% sinogram-shaped, i.e. numel(A) is a positive multiple of Ndist*Nang.
% Any trailing dimensions (TOF bins, multiple sinograms, dynamic frames,
% etc.) are treated as independent columns and interpolated identically.
N = Ndist * Nang;
if isempty(A) || numel(A) == 0 || mod(numel(A), N) ~= 0
    out = A;
    return
end
sz = size(A);
ncols = numel(A) / N;
Ar = reshape(A, N, ncols);
if ~isempty(W)
    B = W * double(Ar);
else
    B = arcGriddataFallback(x_o, y_o, x, y, Ndist, Nang, method, Ar);
end
if isa(A, 'double')
    B = double(B);
else
    B = single(B);
end
out = reshape(B, sz);
end

function out = arcGriddataFallback(x_o, y_o, x, y, Ndist, Nang, method, Ar)
% Fallback interpolation (used only for methods other than 'linear' and
% 'nearest', e.g. 'natural', 'cubic', 'v4'), based on the same periodic
% (angle, signed distance) points as arcInterpWeights, but evaluated with
% griddata. This may be considerably slower than the sparse-matrix path.
N = Ndist * Nang;
[th_o, s_o] = lorDistance(x_o, y_o);
[th_n, s_n] = lorDistance(x, y);
ds = (max(s_o) - min(s_o)) / (Ndist - 1);
u_o = th_o / (pi / Nang);
v_o = s_o / ds;
u_n = th_n / (pi / Nang);
v_n = s_n / ds;

Pu = [u_o - Nang; u_o; u_o + Nang];
Pv = [-v_o; v_o; -v_o];
srcIdx = repmat((1 : N)', 3, 1);

key = round([Pu, Pv] * 1e6);
[~, ~, ic] = unique(key, 'rows', 'stable');
numUnique = max(ic);
Pu_u = accumarray(ic, Pu, [numUnique, 1], @mean);
Pu_v = accumarray(ic, Pv, [numUnique, 1], @mean);

ncols = size(Ar, 2);
out = zeros(size(Ar));
for cc = 1 : ncols
    vals = double(Ar(:, cc));
    valsRep = vals(srcIdx);
    avgVals = accumarray(ic, valsRep, [numUnique, 1], @mean);
    outC = griddata(Pu_u, Pu_v, avgVals, u_n, v_n, method);
    outC(isnan(outC)) = 0;
    out(:, cc) = outC;
end
end

function [th, s] = lorDistance(x, y)
% Canonical (angle, signed distance) parametrisation of a set of LORs,
% each given by two endpoints (columns 1 and 2 of x and y), folded into
% the [0, pi) angular range (see arcInterpWeights).
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
