function [xyLookup, zLookup] = sinogramPairLookup(det_per_ring, rings, Ndist, Nang, ring_difference, span, segment_table, NSinos, ndist_side)
%SINOGRAMPAIRLOOKUP Lookup tables from detector/ring pairs to the sinogram
%
% Computes the lookup tables from the detector pairs and the ring pairs to
% the sinogram bins and planes. The tables reproduce the sinogram indexing
% of saveSinogram (source/cpp/saveSinogram.h) for the raw detector order,
% i.e. the table entries give the sinogram element that the loaders would
% store an event with the given detector and ring indices to. The sinogram
% element of an event is
%   bin + plane * Ndist * Nang
% where the bin is xyLookup (bits 0 - 30) and the plane is zLookup. The
% bit 31 of xyLookup tells that the two rings of the event have to be
% exchanged before using zLookup.
%
% INPUTS:
%   det_per_ring = Number of detectors per ring
%   rings = Number of rings
%   Ndist = Number of radial bins of the sinogram
%   Nang = Number of angles of the sinogram
%   ring_difference = Maximum ring difference
%   span = Span (axial mashing)
%   segment_table = Number of planes in each segment (options.segment_table,
%   not cumulative). Not used with span <= 1.
%   NSinos = Number of sinogram planes
%   ndist_side = (optional) Which side the extra radial bin of an even
%   Ndist is on (1 or -1). The default is 1.
%
% OUTPUTS:
%   xyLookup = uint32 column vector with det_per_ring^2 elements. The element
%   of the detector pair (a, b), where a and b are zero-based, is
%   xyLookup(a + b * det_per_ring + 1). 4294967295 (0xFFFFFFFF) if the pair
%   is outside the sinogram.
%   zLookup = uint32 column vector with rings^2 elements. The element of the
%   ring pair (r1, r2), zero-based, is zLookup(r1 + r2 * rings + 1). The
%   sinogram plane (zero-based) or 4294967295 if the pair is outside the
%   sinogram (ring difference).
%
% See also saveSinogram.h, loadCorrections

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2026 Ville-Veikko Wettenhovi
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
if nargin < 9 || isempty(ndist_side)
    ndist_side = 1;
end
xyLookup = transaxialLookup(double(det_per_ring), double(Ndist), double(Nang), double(ndist_side));
zLookup = axialLookup(double(rings), double(ring_difference), double(span), double(segment_table(:)), double(NSinos));
end

function xyLookup = transaxialLookup(d, Ndist, Nang, ndist_side)
% The sinogram bin of every detector pair, see saveSinogram.h
half = floor(d / 2);
if floor(half / Nang) < 1
    error(['The number of sinogram angles (' num2str(Nang) ') can be at most the number of detectors per ring divided by two (' num2str(half) ')!'])
end
pos = (0 : d - 1)';
% The first dimension is the first detector, the second dimension the second detector
ringPos1 = pos;
ringPos2 = pos';
xa = max(ringPos1, ringPos2);
ya = min(ringPos1, ringPos2);
j = floor(mod(xa + ya + half + 1, d) / 2);
i = abs(xa - ya - half);
negative = ya < j | j + half < xa;
i(negative) = -i(negative);
% The "swap corners" of the sinogram
swap = (j * 2) < -i | i <= (j - half) * 2;
% The default value for an even Ndist is the same as in saveSinogram
if mod(Ndist, 2) == 0 && ndist_side == 0
    side = 1;
else
    side = ndist_side;
end
if mod(Ndist, 2) == 0
    accepted = i <= floor(Ndist / 2) + min(0, side) & i >= -floor(Ndist / 2) + max(0, side);
    i = i + floor(Ndist / 2) - max(0, side);
else
    accepted = i <= floor(Ndist / 2) & i >= -floor(Ndist / 2);
    i = i + floor(Ndist / 2);
end
j = floor(j / floor(half / Nang));
% The sinogram bin, with the bit 31 set when the rings of the event have to be exchanged
flip = xor(swap, ringPos2 > ringPos1);
xyLookup = i + j * Ndist + 2^31 * flip;
xyLookup(~accepted) = 2^32 - 1;
xyLookup = uint32(xyLookup(:));
end

function zLookup = axialLookup(rings, ring_difference, span, segment_table, NSinos)
% The sinogram plane of every ring pair, see saveSinogram.h
ring = (0 : rings - 1)';
% The first dimension is the first ring, the second dimension the second ring
ring1 = ring;
ring2 = ring';
if span <= 1
    plane = ring2 * rings + ring1;
else
    seg = cumsum(segment_table);
    erotus = ring1 - ring2;
    summa = ring1 + ring2;
    absErotus = abs(erotus);
    segment = floor((absErotus + floor(span / 2)) / span);
    outsideSegment0 = absErotus > floor(span / 2);
    needed = outsideSegment0 & absErotus <= ring_difference;
    if any(needed(:)) && numel(seg) < max(segment(needed)) * 2
        error(['The segment table has too few elements for the ring difference (' num2str(ring_difference) ')!'])
    end
    segmentIndex = min(max(segment - 1, 0) * 2 + (erotus >= 0) + 1, numel(seg));
    plane = summa;
    planeOutside = summa - (floor(span / 2) * (segment * 2 - 1) + segment) + seg(segmentIndex);
    plane(outsideSegment0) = planeOutside(outsideSegment0);
end
valid = abs(ring1 - ring2) <= ring_difference & plane >= 0 & plane < NSinos;
plane(~valid) = 2^32 - 1;
zLookup = uint32(plane(:));
end
