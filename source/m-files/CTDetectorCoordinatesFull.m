function [x,y,z] = CTDetectorCoordinatesFull(angles,nProjections,xSize,ySize,varargin)
%CTDETECTORCOORDINATESFULL Computes the detector coordinates of the CT
%projections
%   This function behaves identically to CTDetectorCoordinates except that
%   the detector coordinates are output for each detector crystal/pixel and
%   for each projection. The total size is thus xSize * ySize *
%   numel(angles).  Bed offset is not supported.
%
%   Examples:
%       [x,y,z] = CTDetectorCoordinatesFull(angles,nProjections,xSize,ySize,x,y,z,dPitch)
%       [x,y,z] = CTDetectorCoordinatesFull(angles,nProjections, xSize, ySize, sourceToDetector,sourceToCRot,dPitch,horizontalOffset, verticalOffset, bedOffset, uCenter, vCenter)
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

if nargin >= 7 && ~isempty(varargin{1}) && ~isempty(varargin{2}) && ~isempty(varargin{3}) && numel(varargin{1}) > 1
    x = varargin{1};
    y = varargin{2};
    z = varargin{3};
    if numel(varargin{4}) == 2
        dPitchX = varargin{4}(1);
        dPitchY = varargin{4}(2);
    else
        dPitchX = varargin{4};
        dPitchY = varargin{4};
    end
    if nargin >= 9 && ~isempty(varargin{5})
        pitchRoll = varargin{5};
    else
        pitchRoll = [];
    end
    if size(x,2) > 1 && size(x,1) > 1
        xx = x(:,2);
        yy = y(:,2);
        zz = z(:,2);
    elseif size(x,2) > 1 && size(x,1) == 1
        xx = x';
        yy = y';
        zz = z';
        clear x y z
    else
        xx = x;
        yy = y;
        zz = z;
        clear x y z
    end
else
    if nargin >= 5 && ~isempty(varargin{1})
        sourceToDetector = varargin{1};
    else
        error('Source to detector distance required!')
    end
    if nargin >= 6 && ~isempty(varargin{2})
        sourceToCRot = varargin{2};
    else
        error('Source to center-of-rotation distance required!')
    end
    if nargin >= 7 && ~isempty(varargin{3})
        if numel(varargin{3}) == 2
            dPitchX = varargin{3}(1);
            dPitchY = varargin{3}(2);
        else
            dPitchX = varargin{3};
            dPitchY = varargin{3};
        end
    end
    if nargin >= 8 && ~isempty(varargin{4})
        horizontalOffset = varargin{4};
    else
        horizontalOffset = [];
    end
    if nargin >= 9 && ~isempty(varargin{5})
        verticalOffset = varargin{5};
    else
        verticalOffset = [];
    end
    if nargin >= 10 && ~isempty(varargin{6})
        bedOffset = varargin{6};
    else
        bedOffset = [];
    end
    if nargin >= 11 && ~isempty(varargin{7})
        uCenter = varargin{7};
    else
        uCenter = [];
    end
    if nargin >= 12 && ~isempty(varargin{8})
        vCenter = varargin{8};
    else
        vCenter = [];
    end
    if nargin >= 13 && ~isempty(varargin{9})
        pitchRoll = varargin{9};
    else
        pitchRoll = [];
    end
    [x,y,z] = CTDetSource(angles,nProjections,sourceToDetector,sourceToCRot, horizontalOffset, verticalOffset, bedOffset, uCenter, vCenter);
    xx = x(:,2);
    yy = y(:,2);
    zz = z(:,2);
end
uV = CTDetectorCoordinates(angles,pitchRoll);

indeksiX = repmat((-ySize / 2 + .5 : ySize / 2 - .5)', nProjections, 1); 
indeksiZ = repmat((-xSize / 2 + .5 : xSize / 2 - .5)', nProjections, 1);
if size(uV,1) == 2
    uV = [repelem(uV(1,:)' * dPitchX, ySize, 1) .* indeksiX, repelem(uV(2,:)' * dPitchX, ySize, 1) .* indeksiX];
    xx = repelem(xx, ySize, 1) + uV(:,1);
    yy = repelem(yy, ySize, 1) + uV(:,2);
    xx = single(repmat(reshape(xx, ySize, []), xSize, 1));
    xx = xx(:);
    yy = single(repmat(reshape(yy, ySize, []), xSize, 1));
    yy = yy(:);
    zz = single(repelem(repelem(zz(:,1), xSize, 1) + indeksiZ * dPitchY, ySize, 1));
else
    uV = [repelem(uV(1,:)' * dPitchX, ySize, 1) .* indeksiX, repelem(uV(2,:)' * dPitchX, ySize, 1) .* indeksiX, repelem(uV(3,:)' * dPitchY, xSize, 1) .* indeksiZ,...
        repelem(uV(4,:)' * dPitchX, ySize, 1) .* indeksiX, repelem(uV(5,:)' * dPitchX, ySize, 1) .* indeksiX, repelem(uV(6,:)' * dPitchY, xSize, 1) .* indeksiZ];
    xx = repelem(xx, ySize, 1) + uV(:,1) + uV(:,4);
    yy = repelem(yy, ySize, 1) + uV(:,2) + uV(:,5);
    xx = single(repmat(reshape(xx, ySize, []), xSize, 1));
    xx = xx(:);
    yy = single(repmat(reshape(yy, ySize, []), xSize, 1));
    yy = yy(:);
    zz = single(repelem(repelem(zz(:,1), xSize, 1) + uV(:,3) + uV(:,6), ySize, 1));
end
if exist('x','var')
    x = [repelem(x(:,1), xSize * ySize, 1) xx];
    y = [repelem(y(:,1), xSize * ySize, 1) yy];
    z = [repelem(z(:,1), xSize * ySize, 1) zz];
else
    x = xx;
    y = yy;
    z = zz;
end