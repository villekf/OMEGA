function [index, pituus, varargout] = subset_angles(options, varargin)
%SUBSET_ANGLES Form the subset indices for the subset_type 6 and angles for subset_type 7

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019-2026 Ville-Veikko Wettenhovi
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if options.use_raw_data
    [x,y] = detector_coordinates(options);
    % Determine the rings of each detector
    apu = idivide(varargin{1}, uint16(options.det_per_ring),'ceil');
    % Same ring
    ind = apu(:,1) == apu(:,2);
    indi = find(ind);
    % Make all detector numbers constrained in [1, det_per_ring]
    varargin{1} = mod(varargin{1}, options.det_per_ring);
    varargin{1}(varargin{1} == 0) = options.det_per_ring;
    x2 = single(x);
    y2 = single(y);
    % Coordinates of detectors that are on the same rings (2D)
    x = single(x(varargin{1}(ind,:)));
    y = single(y(varargin{1}(ind,:)));
else
    if options.sampling > 1 && ~options.use_raw_data && ~options.precompute_lor
        options.Ndist = options.Ndist * options.sampling;
    end
    [~, ~, xp, yp] = detector_coordinates(options);
    [x,y] = sinogram_coordinates_2D(options, xp, yp);
    sino_length = options.Nang*options.Ndist;
end
% Lengths of the triangle sides
jokux = abs(x(:,1)-x(:,2));
jokuy = abs(y(:,1)-y(:,2));
% Indices on the other side of the circle (angle over 90 degrees)
ind = x(:,1) > x(:,2) & y(:,1) > y(:,2) | x(:,2) > x(:,1) & y(:,2) > y(:,1);
% LORs that are vertical (or horizontal) up to numerical precision are
% treated as exactly vertical (0 degrees) or horizontal (90 degrees), so that
% rounding errors (e.g. single vs. double precision) cannot move them to the
% other side of the circle
vertical = jokux <= 1e-6 * max(jokuy);
horizontal = jokuy <= 1e-6 * max(jokux);
ind(vertical | horizontal) = false;
% Angle
kulma = atand(jokux./jokuy);
kulma(horizontal) = 90;
kulma(vertical) = 0;
kulma(ind) = kulma(ind) + 90;
if nargout >= 3
    % Only the LOR angles are requested (e.g. goldenAngleSubsets), the
    % subset indices and options.n_angles are not needed
    varargout{1} = kulma;
    index = [];
    pituus = [];
    return
end
% Combine similar angles together. The angles are first snapped to 0.001
% degrees, so that LORs with an angle that is an exact half degree (up to
% numerical precision) are rounded identically in single and double precision
kulma = round(round(kulma * 1e3) / 1e3);
kulma(kulma==180) = 179;
kulmaOrig = kulma;
[kulma, I] = sort(kulma);
n_angles = options.n_angles;
if mod(180, n_angles) ~= 0
    error('options.n_angles has to divide 180 evenly!')
end
maara = 180/n_angles;

if options.use_raw_data
    maara2 = 360/n_angles;
    index1 = cell(maara + maara2 * (options.rings-1),1);
    pituus = zeros(maara + maara2 * (options.rings-1),1,'uint32');
    % 2D case
    for ii = 1 : maara
        % Take the indices of angles that are between the defined number of
        % angles
        index1{ii} = uint32(indi(I(kulma >= (ii-1)*n_angles & kulma < ii*n_angles)));
        pituus(ii) = length(index1{ii});
    end
else
    % Sinogram data. 2D data (only the direct planes are used): 180/n_angles
    % subsets, the subset j contains the bins with a (rounded) angle in
    % [j*n_angles, (j+1)*n_angles) of all (direct) planes.
    % 3D data: 360/n_angles subsets. The oblique segments of all segment
    % pairs belong to the subset that contains their full angle (the angle of
    % the first member of a pair, angle + 180 for the second member). The
    % direct planes are split by the plane parity: the bins with the angle
    % range j of the even planes (0, 2, ...) go to the subset j and those of
    % the odd planes to the subset j + 180/n_angles, i.e. every subset spans
    % the whole axial range. Order inside a subset: the direct planes first,
    % then the oblique segments in increasing segment order; every block is
    % plane-major (a plane after the other), the bins of a plane in
    % ascending sinogram order. The indices are 1-based.
    segment_table = options.segment_table(:);
    nPairs = floor(length(segment_table)/2);
    use3D = nPairs > 0 && ~(options.NSinos < options.TotSinos);
    if use3D
        nSubsets = 2 * maara;
        nDirect = segment_table(1);
        % First plane of the first member of each segment pair
        firstPlane = zeros(nPairs,1);
        nPlanes = zeros(nPairs,1);
        for rd = 1 : nPairs
            firstPlane(rd) = sum(segment_table(1:rd*2-1));
            nPlanes(rd) = segment_table(2+(rd-1)*2);
        end
    else
        nSubsets = maara;
        nDirect = options.NSinos;
    end
    index1 = cell(nSubsets,1);
    pituus = zeros(nSubsets,1,'uint32');
    for kk = 1 : nSubsets
        j = mod(kk-1, maara);
        par = floor((kk-1) / maara);
        bins = find(kulmaOrig >= j*n_angles & kulmaOrig < (j+1)*n_angles);
        if use3D
            planes = par : 2 : nDirect - 1;
            for rd = 1 : nPairs
                planes = [planes, firstPlane(rd) + nPlanes(rd)*par + (0 : nPlanes(rd) - 1)]; %#ok<AGROW>
            end
        else
            planes = 0 : nDirect - 1;
        end
        ind1 = double(bins(:)) + double(planes(:))' * sino_length;
        index1{kk} = uint32(ind1(:));
        pituus(kk) = length(index1{kk});
    end
    index = cell2mat(index1);
end

% 3D case
if options.use_raw_data
    ll = 1;
    ring_diff = zeros(options.rings,1);
    % Constant increment scheme
    for kk = 1 : options.rings
        ring_diff(kk+1) = mod(ring_diff(kk) + round(0.7*options.rings), options.rings + 1);
        while any(ismember(ring_diff(kk+1), ring_diff(1:kk)))
            ring_diff(kk+1) = ring_diff(kk+1) + 1;
        end
    end
    ring_diff(1) = [];
    for kk = ring_diff
        % Choose the detectors that have an angle of most 180 degrees, and
        % angles that are betweem 180 and 360.
        ind1 = apu(:,1) == (apu(:,2) + kk);
        ind2 = apu(:,2) == (apu(:,1) + kk);
        ind = logical(ind1 + ind2);
        indi = find(ind);
        xx = x2(varargin{1}(ind,:));
        yy = y2(varargin{1}(ind,:));
        jokux = abs(xx(:,1)-xx(:,2));
        jokuy = abs(yy(:,1)-yy(:,2));
        ind = xx(:,1) > xx(:,2) & yy(:,1) > yy(:,2) | xx(:,2) > xx(:,1) & yy(:,2) > yy(:,1);
        ind2 = ind2(ind);
        kulma = atand(jokux./jokuy);
        kulma(ind) = kulma(ind) + 90;
        kulma(ind2) = kulma(ind2) + 180;
        kulma = round(kulma);
        [kulma, I] = sort(kulma);
        for ii = 1 : maara2
            index1{maara + ll} = uint32(indi(I(kulma >= (ii-1)*n_angles & kulma < ii*n_angles)));
            pituus(ii) = length(index1{ii});
            ll = ll + 1;
        end
    end
    index = cell2mat(index1);
end
