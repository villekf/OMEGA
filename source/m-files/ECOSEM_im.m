function im = ECOSEM_im(im, epps, D, COSEM_apu, OSEM_apu)
%ECOSEM_IM Computes the Enhanced COSEM (ECOSEM) estimates
%
% Example:
%   im = ECOSEM_im(im, epps, D, COSEM_apu, OSEM_apu)
% INPUTS:
%   im = The current estimate
%   epps = Small constant to prevent division by zero
%   D = Sum of the complete data system matrix (D = sum(B,2), where B =
%   [A_1,A_2,...,A_subsets]
%   COSEM_apu = Current COSEM estimate
%   OSEM_apu = Current OSEM estimate
%
% OUTPUTS:
%   im = The updated estimate
%
%   See also COSEM_im, ACOSEM_im, COSEM_OSL, OSEM_im

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi, Samuli Summala
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
pz_eco_apuw = im;
alpha_eco = 1;
im = alpha_eco*OSEM_apu + (1 - alpha_eco)*COSEM_apu;
apu = sum(D.*(-COSEM_apu.*log(pz_eco_apuw + epps) + pz_eco_apuw));
apu2 = sum(D.*(-COSEM_apu.*log(im + epps) + im));
while (alpha_eco > 0.0096 && apu < apu2)
    alpha_eco = alpha_eco*0.9;
    im = alpha_eco*OSEM_apu + (1 - alpha_eco)*COSEM_apu;
    apu2 = sum(D.*(-COSEM_apu.*log(im + epps) + im));
end
if alpha_eco <= 0.0096
    im = COSEM_apu;
end