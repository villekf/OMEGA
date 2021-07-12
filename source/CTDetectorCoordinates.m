function [x,y,z] = CTDetectorCoordinates(angles,sourceToDetector,sourceToCRot,dPitch,xSize,ySize,varargin)
%CTDETECTORCOORDINATES Computes the detector coordinates of the CT
%projections
%   This function computes the detector coordinates for one detector panel
%   in each projection/bed position. x and y coordinates correspond to the
%   columns and rows of the projection image. z coordinates correspond to
%   the axial coordinates.
%
%   Examples:
%       [x,y,z] = CTDetectorCoordinates(angles,sourceToDetector,sourceToCRot,dPitch,xSize,ySize)
%       [x,y,z] = CTDetectorCoordinates(angles,sourceToDetector,sourceToCRot,dPitch,xSize,ySize, horizontalOffset, verticalOffset)
%       [x,y,z] = CTDetectorCoordinates(angles,sourceToDetector,sourceToCRot,dPitch,xSize,ySize, horizontalOffset, verticalOffset, bedOffset)
%   Inputs:
%       angles = The projection angles
%       sourceToDetector = Distance from the source to the detector panel
%       (mm)
%       sourceToCRot = Distance from the source to the center of rotation
%       (mm)
%       dPitch = Distance between adjacent detectors (detector pitch) (mm)
%       xSize = Number of columns in the projection image
%       ySize = Number of rows in the projection image
%       horizontalOffset = (optional) Possible horizontal offset of the
%       source (i.e. the center of rotation is not in the origin)
%       verticalOffset = (optional) same as above, but for vertical
%       direction
%       bedOffset = The bed offsets for each different bed positions in
%       step-and-shoot mode. I.e. this is the amount that the bed moves in
%       each bed position. First one should be 0. (mm)
%
% See also get_coordinates, CTDetectorCoordinatesFull
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021 Ville-Veikko Wettenhovi
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
detSizeTr = ySize * dPitch;
detSizeAx = xSize * dPitch;

if nargin >= 7 && ~isempty(varargin) && ~isempty(varargin{1})
    horizontalOffset = (varargin{1});
    if nargin >= 8 && ~isempty(varargin{2})
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
if nargin >= 10 && ~isempty(varargin) && ~isempty(varargin{5}) && ~isempty(varargin{4})
    uCenter = (varargin{4});
    vCenter = (varargin{5});
else
    uCenter = [];
    vCenter = [];
end

detCoordY = - (sourceToDetector - sourceToCRot);
detCoordX = detSizeTr / 2 - dPitch / 2;
if isempty(vCenter)
    detCoordZ = repmat(dPitch/2,numel(angles),1);
else
    detCoordZ = repmat(dPitch/2,numel(angles),1) - vCenter;
end
angles = reshape(angles, 1, 1, []);

R = permute([cos(angles) -sin(angles); sin(angles) cos(angles)], [2 1 3]);
if isempty(uCenter)
    if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','9.1')
        XY = squeeze(sum(bsxfun(@times, R, [detCoordX,detCoordY]),2))';
    else
        XY = squeeze(sum(R .* [detCoordX,detCoordY],2))';
    end
else
    R2 = permute([cos(pi/2 - angles) -sin(pi/2 - angles); sin(pi/2 - angles) cos(pi/2 - angles)], [2 1 3]);
    if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','9.1')
        XY = squeeze(sum(bsxfun(@times, R, [detCoordX,detCoordY]) + bsxfun(@times, R2, reshape([uCenter,uCenter], 1, 2, [])),2))';
    else
        XY = squeeze(sum(R .* [detCoordX,detCoordY] + R2 .* reshape([uCenter,uCenter], 1, 2, []),2))';
    end
end

sourceCoordY = sourceToCRot - verticalOffset * 1;
sourceCoordX = -horizontalOffset * 1;
sourceCoordZ = detSizeAx / 2;

testi = reshape([sourceCoordX,sourceCoordY]', 2, 1, []);
sXY = squeeze(sum(R .* permute(testi, [2 1 3]),2))';

x = [XY(:,1) sXY(:,1)];
y = -[XY(:,2) sXY(:,2)];
z = [reshape(detCoordZ,[],1) repmat(sourceCoordZ,numel(detCoordZ),1)];

if nargin >= 9 && ~isempty(varargin) && ~isempty(varargin{3}) && numel(varargin{3}) > 1
    x = repmat(x, numel(varargin{3}),1);
    y = repmat(y, numel(varargin{3}),1);
    z = bsxfun(@plus, repmat(z, numel(varargin{3}),1), repelem(varargin{3}, numel(angles),1));
end


end

