function [x,y] = computeCoordinates(blocks_per_ring, transaxial_multip,cryst_per_block, diameter, cr_p, usePseudo, transaxialGaps, firstBlockAngle)
%COMPUTECOORDINATES Computes the coordinates for each detector on a single
%ring
%   This function computes the transaxial detector coordinates (x,y) for
%   all crystals on one ring. Can be used also for cases with pseudo
%   detectors.
%
%   The ring is a regular polygon with N = blocks_per_ring *
%   transaxial_multip flat facets (blocks). The front face of each block is
%   at the perpendicular distance diameter/2 from the center, i.e. diameter
%   is the distance between two opposing (parallel) block faces. The normal
%   of facet 0 is at the angle firstBlockAngle (degrees, default 0, i.e.
%   facet 0 is centered on the +x axis) and the facets follow each other
%   counter-clockwise. The crystals of a block are placed at the pitch cr_p
%   along the facet, centered on the facet, the first crystal of block 0
%   having the smallest y coordinate.
%
%   Coordinate system: the origin is at the scanner center and z is the
%   scanner axis. +x is the direction of the image's first dimension (Nx, the
%   row direction) and +y the direction of its second dimension (Ny, the
%   column direction); angles are counter-clockwise, i.e. from +x towards +y.
%   With firstBlockAngle = 0 the center of block 0 is at (x, y) =
%   (diameter/2, 0) (the DOI is added to the radius) and its crystals run
%   towards +y. The angle is applied before offangle (detector renumbering)
%   and before flip_image.
%
%   A pseudo detector (usePseudo) is placed one cr_p after the last crystal
%   of each block.
%
%   transaxialGaps (optional) is a vector of cryst_per_block - 1 elements,
%   the extra spacing (mm) after crystal j of each block. The block stays
%   centered on its facet. Empty or zero = no gaps.
%
%   firstBlockAngle (optional) is the angle (degrees, counter-clockwise from
%   the +x axis, i.e. from the direction of the image's first dimension
%   towards its second) of the normal of the first block. Default 0.
%
% See also detector_coordinates

if nargin < 7 || isempty(transaxialGaps) || ~any(transaxialGaps(:))
    transaxialGaps = [];
elseif numel(transaxialGaps) ~= cryst_per_block - 1
    error(['transaxialGaps has to contain cryst_per_block - 1 = ' num2str(cryst_per_block - 1) ' elements, but it has ' num2str(numel(transaxialGaps)) ' elements!'])
end
if nargin < 8 || isempty(firstBlockAngle)
    firstBlockAngle = 0;
end
nBlocks = blocks_per_ring * transaxial_multip;

% Tangential position of each crystal along the facet, centered on the facet
t = ((0:cryst_per_block - 1) - (cryst_per_block - 1) / 2) * cr_p;
if ~isempty(transaxialGaps)
    g = double(transaxialGaps(:))';
    t = t + [0, cumsum(g)] - sum(g) / 2;
end
if usePseudo
    % One pseudo detector, one cr_p after the last crystal
    t = [t, t(end) + cr_p];
end

% Facet normal angles (facet 0 at firstBlockAngle)
theta = 2 * pi * (0:nBlocks - 1)' / nBlocks + firstBlockAngle * pi / 180;
R = diameter / 2;
x = reshape((R * cos(theta) - sin(theta) * t)', [], 1);
y = reshape((R * sin(theta) + cos(theta) * t)', [], 1);
end
