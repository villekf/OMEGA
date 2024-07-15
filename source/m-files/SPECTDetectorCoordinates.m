function [x,y,z] = SPECTDetectorCoordinates(angles,headToCRot,dPitchTr,dPitchAx,xSize,ySize,nHeads,varargin)
%SPECTDETECTORCOORDINATES Computes detector-source coordinates for SPECT
%imaging
%   Detailed explanation goes here
detSizeTr = ySize * dPitchTr;
detSizeAx = xSize * dPitchAx;
if nHeads > 1
    headAngles = varargin{1};
else
    headAngles = 0;
end

angles = reshape(angles, 1, 1, []) * pi / 180;
% R = [cosd(angles) -sind(angles); sind(angles) cosd(angles)];
z = repmat(repmat(-detSizeAx / 2 + dPitchAx/2, ySize, 1) + linspace(0, dPitchAx * (xSize - 1), xSize), 1, 1, numel(angles) * nHeads);
detCoordY = repmat(headToCRot, ySize, 1);
detCoordY2 = repmat(-headToCRot, ySize, 1);
detCoordX = repmat(detSizeTr / 2 - dPitchTr / 2, ySize, 1) - linspace(0, dPitchTr * (ySize - 1), ySize)';
x = [];
y = [];

for kk = 1 : nHeads
    rAngles = angles + headAngles(kk);
    if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','9.1')
        XY = bsxfun(@times,exp(rAngles*1i), ([1 1i]*[detCoordX';detCoordY']));
    else
        XY = exp(rAngles*1i) .* ([1 1i]*[detCoordX';detCoordY']);
        XY2 = exp(rAngles*1i) .* ([1 1i]*[detCoordX';detCoordY2']);
    end
    XY = [real(XY); imag(XY)];
    XY = repmat(XY, 1, xSize, 1);
    XY2 = [real(XY2); imag(XY2)];
    XY2 = repmat(XY2, 1, xSize, 1);
    
    xx = XY(1, :, :);
    xx2 = XY2(1, :, :);
    x = [x; xx(:) xx2(:)];
    yy = XY(2, :, :);
    yy2 = XY2(2, :, :);
    y = [y; yy(:) yy2(:)];
end
z = [z(:) z(:)];

end