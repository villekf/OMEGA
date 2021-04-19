function [x,y,z] = CTDetectorCoordinatesFull(angles,sourceToDetector,sourceToCRot,dPitch,xSize,ySize,varargin)
%CTDETECTORCOORDINATESFULL Computes the detector coordinates of the CT
%projections
%   This function behaves identically to CTDetectorCoordinates except that
%   the detector coordinates are output for each detector crystal/pixel and
%   for each projection. The total size is thus xSize * ySize *
%   numel(angles).  Bed offset is not supported.
%
%   Examples:
%       [x,y,z] = CTDetectorCoordinatesFull(angles,sourceToDetector,sourceToCRot,dPitch,xSize,ySize)
%       [x,y,z] = CTDetectorCoordinatesFull(angles,sourceToDetector,sourceToCRot,dPitch,xSize,ySize, horizontalOffset, verticalOffset)
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
%
% See also get_coordinates, CTDetectorCoordinates
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

if nargin >= 8 && ~isempty(varargin) && ~isempty(varargin{1}) && ~isempty(varargin{2})
    horizontalOffset = (varargin{1});
    verticalOffset = varargin{2};
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

if max(abs(angles(:))) > 2*pi
    angles = angles * (pi / 180);
end

detCoordY = - (sourceToDetector - sourceToCRot);
detCoordX = -detSizeTr / 2 + dPitch / 2;

X = linspace(-detCoordX,detCoordX,ySize)';
Y = repmat(detCoordY,ySize,1);
X = repmat(X, xSize,1);
Y = repmat(Y, xSize,1);
detCoordZ = repmat(repelem(linspace(dPitch/2,detSizeAx-dPitch/2,xSize)',ySize),numel(angles),1);
angles = reshape(angles, 1, 1, []);

R = permute([cos(angles) -sin(angles); sin(angles) cos(angles)], [1 2 3]);
XY = zeros(numel(angles) * xSize * ySize,2);

sourceCoordY = sourceToCRot - verticalOffset * 1;
sourceCoordX = -horizontalOffset * 1;
sourceCoordZ = detSizeAx / 2;

for kk = 1 : numel(angles)
    XY(1 + (kk - 1) * xSize * ySize : kk * xSize * ySize,:) = (R(:,:,kk) * [X Y]')';
end

testi = reshape([sourceCoordX,sourceCoordY]', 2, 1, []);
sXY = repelem(squeeze(sum(R .* permute(testi, [2 1 3]),2))',xSize*ySize,1);

x = [XY(:,1) sXY(:,1)];
y = -[XY(:,2) sXY(:,2)];
z = [reshape(detCoordZ,[],1) repmat(sourceCoordZ,numel(detCoordZ),1)];
end