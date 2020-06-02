function epsilon_mramla = MBSREM_epsilon(Sino, epps, randoms_correction, randoms, D)
%MBSREM_epsilon Computes the epsilon value for MBSREM/MRAMLA
%   INPUTS:
%   Sino = Measurements at current time-step
%   epps = Small constant to prevent division by zero
%   randoms_correction = Is randoms and/or scatter correction included
%   randoms = Number of randoms and/or scattered coincidences
%   D = Sum of the system matrix (measurements)

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
% I = find(Sino);
if randoms_correction == 1
    P_Sino = double(Sino(Sino > 0 & randoms == 0));
    apu = D + double(randoms);
    apu = sum(double(Sino) .* log(apu) - apu);
else
    P_Sino = double(Sino(Sino > 0));
    apu = sum(D .* log(D) - D);
end
% hk_summa = zeros(nnz(I),1);
% lhk = 1;
% hk_summa = zeros(size(Sino));
% hk_summa(Sino > 0) = double(Sino(Sino > 0)) .* log(double(Sino(Sino > 0))) - double(Sino(Sino > 0));
hk_summa = double(Sino) .* log(double(Sino)) - double(Sino);
hk_summa(isnan(hk_summa)) = 0;
% hk_summa = double(P_Sino) .* log(P_Sino) - double(P_Sino);
if randoms_correction == 1
    hk_summa = bsxfun(@minus, sum(hk_summa), hk_summa(Sino > 0 & randoms == 0));
%     hk_summa = bsxfun(@minus, sum(hk_summa), hk_summa);
else
    hk_summa = bsxfun(@minus, sum(hk_summa), hk_summa(Sino > 0));
%     hk_summa = bsxfun(@minus, sum(hk_summa), hk_summa);
end
% apu = sum(P_Sino .* log(P_Sino) - P_Sino);
% apu = length(P_Sino);
% apu = sum(P_Sino .* log(1) - 1);
% for hkl = I
%     apu = I ~= hkl;
%     hk_summa(lhk) = sum(double(Sino(I(apu))).*log(double(Sino(I(apu)))) - double(Sino(I(apu))));
%     lhk = lhk + 1;
% end
epsilon_mramla = min(min([P_Sino';exp((apu - hk_summa)./double(P_Sino))']));
if epsilon_mramla <= 0
    epsilon_mramla = epps;
end

