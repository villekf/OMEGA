function [x,y,z] = CTDetSource(angles, nProjections, sourceToDetector, sourceToCRot, varargin)
%CTDETSOURCE Computes the center coordinates for the X-ray source and
%detector panel for each projection
%   If not input by the user, this function will compute the center
%   coordinates for both the source and the detector panel. It is assumed
%   that the rotation is only in two dimensions. Offset values for the
%   source and/or detector can be input as well. Bed offset values can be
%   input as well for step-and-shoot measurement.
%   Examples:
%       [x,y,z] = CTDetSource(angles,nProjections,sourceToDetector,sourceToCRot)
%       [x,y,z] = CTDetSource(angles,nProjections,sourceToDetector,sourceToCRot, horizontalOffset, verticalOffset)
%       [x,y,z] = CTDetSource(angles,nProjections,sourceToDetector,sourceToCRot, horizontalOffset, verticalOffset, bedOffset)
%       [x,y,z] = CTDetSource(angles,nProjections,sourceToDetector,sourceToCRot, horizontalOffset, verticalOffset, bedOffset, uCenter, vCenter)
%   Inputs:
%       angles = The projection angles
%       nProjections = Total number of projections
%       sourceToDetector = Distance from the source to the detector panel
%       (mm)
%       sourceToCRot = Distance from the source to the center of rotation
%       (mm)
%       horizontalOffset = (optional) Possible horizontal offset of the
%       SOURCE (i.e. the center of rotation is not in the origin)
%       verticalOffset = (optional) same as above, but for vertical
%       direction
%       bedOffset = The bed offsets for each different bed positions in
%       step-and-shoot mode. I.e. this is the amount that the bed moves in
%       each bed position. First one should be 0. (mm)
%       uCenter = (optional) Possible horizontal offset of the DETECTOR
%       (i.e. the center of rotation is not in the origin)
%       vCenter = (optional) Possible vertical offset of the DETECTOR
%   Outputs:
%       x, y, z = The center coordinates for each three axes for the
%       detector panel and the source.
%
% See also get_coordinates, CTDetectorCoordinates, CTDetectorCoordinatesFull
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2022 Ville-Veikko Wettenhovi
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

if nargin >= 5 && ~isempty(varargin) && ~isempty(varargin{1})
    horizontalOffset = (varargin{1});
    if nargin >= 6 && ~isempty(varargin{2})
        verticalOffset = varargin{2};
    else
        verticalOffset = 0;
    end
    if numel(horizontalOffset) > numel(verticalOffset)
        horizontalOffset = horizontalOffset(1:numel(angles));
        verticalOffset = repmat(verticalOffset, numel(horizontalOffset)/numel(verticalOffset),1);
    elseif numel(horizontalOffset) < numel(verticalOffset)
        verticalOffset = verticalOffset(1:numel(angles));
        horizontalOffset = repmat(horizontalOffset, numel(verticalOffset)/numel(horizontalOffset),1);
    end
else
    horizontalOffset = 0;
    verticalOffset = 0;
end
if nargin >= 8 && ~isempty(varargin) && ~isempty(varargin{5}) && ~isempty(varargin{4})
    uCenter = (varargin{4});
    vCenter = (varargin{5});
else
    uCenter = [];
    vCenter = [];
end
if numel(angles) == nProjections
    angles = reshape(angles, 1, 1, []);
else
    if size(angles,1) == nProjections
        angles = angles(:,1);
    else
        angles = angles(1,:);
        angles = angles(:);
    end
end

R = [cos(angles) -sin(angles); sin(angles) cos(angles)];
% R = permute([cos(angles) sin(angles); -sin(angles) cos(angles)], [2 1 3]);

sourceCoordY = sourceToCRot;
sourceCoordX = horizontalOffset;
sourceCoordZ = verticalOffset;
if numel(sourceCoordX) == 1
    sourceCoordX = repmat(sourceCoordX, nProjections, 1);
end
if numel(sourceCoordZ) == 1
    sourceCoordZ = repmat(sourceCoordZ, nProjections, 1);
end
sourceCoordY = repmat(sourceCoordY, nProjections, 1);

testi = reshape([sourceCoordX,sourceCoordY]', 2, 1, []);
if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','9.1')
    sXY = squeeze(sum(bsxfun(@times, R, permute(testi, [2 1 3])),2))';
else
    sXY = squeeze(sum(R .* permute(testi, [2 1 3]),2))';
end

detCoordY = - (sourceToDetector - sourceToCRot);
detCoordX = 0;
detCoordZ = 0;
if ~isempty(vCenter)
    detCoordZ = vCenter;
end
if ~isempty(uCenter)
    detCoordX = uCenter;
end
if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','9.1')
    XY = squeeze(sum(bsxfun(@times, R, [detCoordX,detCoordY]),2))';
else
    XY = squeeze(sum(R .* [detCoordX,detCoordY],2))';
end
x = [sXY(:,1) XY(:,1)];
y = [sXY(:,2) XY(:,2)];
if numel(sourceCoordZ) == 1
    z = [repmat(sourceCoordZ,numel(detCoordZ),1) reshape(detCoordZ,[],1)];
elseif numel(detCoordZ) == 1
    z = [reshape(sourceCoordZ,[],1) repmat(detCoordZ,numel(sourceCoordZ),1)];
else
    z = [sourceCoordZ reshape(detCoordZ,[],1)];
end

if nargin >= 7 && ~isempty(varargin) && ~isempty(varargin{3}) && numel(varargin{3}) > 1
    x = repmat(x, numel(varargin{3}),1);
    y = repmat(y, numel(varargin{3}),1);
    z = bsxfun(@plus, repmat(z, numel(varargin{3}),1), repelem(varargin{3} - varargin{3}(end) / 2, numel(angles),1));
end