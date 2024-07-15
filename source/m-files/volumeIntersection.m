function V = volumeIntersection(R, r, b)
%VOLUMEINTERSECTION Computes the volume of intersection between a sphere
%and a cylinder. 
%   This is based on https://doi.org/10.1016/0010-4655(90)90184-3
%   Examples:
%       V = volumeIntersection(R, r, b)
%   Inputs:
%       R = The radius of the tube of response (TOR)
%       r = The radius of the spherical voxel
%       b = The distances from the center of the TOR to the center of the
%       spherical voxel
%   Outputs:
%       V = The volumes corresponding to the distances b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021 Ville-Veikko Wettenhovi
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = max(r^2, (b + R).^2);
B = min(r^2, (b + R).^2);
C = (b - R).^2;
k = abs((B - C) ./ (A - C));
alpha = (B - C) ./ C;
s = (b + R) .* (b - R);
Gamma_f = @(x) 1 ./ ((1 + alpha .* x.^2) .* sqrt(1 - x.^2) .* sqrt(1 - k .* x.^2));
if exist('integral','file') == 2
    Gamma = integral(Gamma_f, 0, 1, 'ArrayValued', true,'AbsTol',1e-6);
else
    Gamma = quadv(Gamma_f, 0, 1);
end
[K,E] = ellipke(k);
heaviside_f = zeros(size(b));
heaviside_f(R - b > 0) = 1;
heaviside_f(R - b == 0) = 1/2;
V = zeros(size(b));
A1 = A(r > b + R);
B1 = B(r > b + R);
C1 = C(r > b + R);
s1 = s(r > b + R);
K1 = K(r > b + R);
E1 = E(r > b + R);
Gamma1 = Gamma(r > b + R);
V(r > b + R) = (4 .* pi)/3 .* r^3 .* heaviside_f(r > b + R) + 4 ./ (3 .* sqrt(A1 - C1)) .* (Gamma1 .* (A1.^2 .* s1) ./ C1 - K1 .*(A1 .* s1 - ((A1 - B1) .* (A1 - C1)) ./ 3)...
    - E1 .* (A1 - C1) .* (s1 + (4*A1 - 2*B1 - 2*C1) ./ 3));

A2 = A(r < b + R);
B2 = B(r < b + R);
C2 = C(r < b + R);
s2 = s(r < b + R);
K2 = K(r < b + R);
E2 = E(r < b + R);
Gamma2 = Gamma(r < b + R);
V(r < b + R) = (4 .* pi)/3 .* r^3 .* heaviside_f(r < b + R) + ((4/3) ./ (sqrt(A2 - C2))) .* (Gamma2 .* ((B2.^2 .* s2) ./ C2) + K2 .* (s2 .* (A2 - 2*B2) + ...
    (A2 - B2) .* ((3*B2 - C2 - 2*A2) / 3)) + E2 .* (A2 - C2) .* (-s2 + (2*A2 - 4*B2 + 2*C2) ./ 3));
end