function LL = form_detector_pairs_raw(rings, det_per_ring)
%form_detector_pairs_raw Forms the detector pairs for raw data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi
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

temp = uint16(rings*det_per_ring);

LL = (uint16(1):temp)';
LL = repelem(LL, flip(LL));

apu = zeros(size(LL),'uint16');

loc = uint32(1);
for kk = uint32(1) : uint32(temp)
    apu(loc: loc + (uint32(temp(end)) - kk)) = (uint16(kk) :  temp)';
    loc = loc + (uint32(temp(end)) - kk) + 1;
end

LL = [LL,apu];
end

