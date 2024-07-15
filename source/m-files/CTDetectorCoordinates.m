function uV = CTDetectorCoordinates(angles,varargin)
%CTDETECTORCOORDINATES Computes the direction vectors for each direction in
%the detector array.
%   This function computes the direction vectors for each detector
%   direction in the detector array. If the detector panel has no pitch
%   and/or roll, then only two vectors are computed. If the pitch and/or
%   roll angles are input, six direction vectors are output. Even if either
%   pitch or roll angle is zero, it still has to be input if the other one
%   is input. If both are zero, they can (and should) be omitted. Note that
%   the direction vectors have to be multiplied by the detector size to get
%   the correct vectors, i.e. the vectors are output as unit vectors.
%
%   Examples:
%       uV = CTDetectorCoordinates(angles)
%       uV = CTDetectorCoordinates(angles,pitchRoll)
%   Inputs:
%       angles = The projection angles (in radians).
%       pitchRoll = (Optional) Pitch and roll angles for the detector. The
%       first column should include the pitch and the second the roll. As
%       above, this needs to be in radians.
%   Outputs:
%       uV = Direction unit vectors for the detector array.
%
% See also get_coordinates, CTDetSource, CTDetectorCoordinatesFull
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

angles = reshape(angles, [], 1);
if isempty(varargin) || isempty(varargin{1})
    uV = single([-sin(angles) cos(angles)]');
    %uV = single([-cos(angles) -sin(angles)]'); % Simuloidulle?
else
    pitchRoll = reshape(varargin{1}, [], 2);
    uV = single([-sin(angles).*cos(pitchRoll(:,1,:,:)) - cos(angles).*sin(pitchRoll(:,1,:,:)).*sin(pitchRoll(:,2,:,:)) ...
        cos(angles).*cos(pitchRoll(:,1,:)) - sin(angles).*sin(pitchRoll(:,1,:,:)).*sin(pitchRoll(:,2,:,:))...
        sin(pitchRoll(:,1,:,:)).*cos(pitchRoll(:,2,:,:)),...
        sin(angles).*sin(pitchRoll(:,1,:,:)) - cos(angles).*cos(pitchRoll(:,1,:,:)) .* sin(pitchRoll(:,2,:,:))...
        -(cos(angles).*sin(pitchRoll(:,1,:,:)) + sin(angles).*cos(pitchRoll(:,1,:,:)).* sin(pitchRoll(:,2,:,:)))...
        cos(pitchRoll(:,2,:,:)) .* cos(pitchRoll(:,1,:,:))]');
end

end

