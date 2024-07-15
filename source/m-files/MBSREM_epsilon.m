function epsilon_mramla = MBSREM_epsilon(Sino, epps, randoms_correction, rand, D, TOF, nBins, varargin)
%MBSREM_epsilon Computes the epsilon value for MBSREM/MRAMLA
%   INPUTS:
%   Sino = Measurements at current time-step
%   epps = Small constant to prevent division by zero (unused)
%   randoms_correction = Is randoms and/or scatter correction included
%   rand = Number of randoms and/or scattered coincidences
%   D = Sum of the system matrix (measurements)
%   TOF = True if TOF data is used, false otherwise
%   nBins = Number of TOF bins
%   CT = Whether the CT likelihood is used

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020-2024 Ville-Veikko Wettenhovi
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
if ~isempty(varargin) && ~isempty(varargin{1})
    CT = varargin{1};
else
    CT = false;
end
if (CT)
    hk_summa = -exp(-Sino) ./ Sino - Sino;
    hk_summa(isnan(hk_summa)) = 0;
    if randoms_correction == 1
        Iind = (Sino > 0 & rand == 0);
        if (sum(Iind) == 0)
            epsilon_mramla = 1e8;
            return;
        end
        P_Sino = Sino(Iind);
        apu = D + rand;
        apu = sum(-exp(-apu) ./ Sino - apu);
    else
        Iind = (Sino > 0);
        P_Sino = Sino(Iind);
        apu = sum(-exp(-D) ./ Sino - D);
    end
    hk_summa = bsxfun(@minus, sum(hk_summa), hk_summa(Iind));
    epsilon = min(P_Sino, log(bsxfun(@minus, apu, hk_summa) ./ P_Sino));
    epsilon_mramla = min(epsilon);
else
    hk_summa = Sino .* log(Sino) - Sino;
    hk_summa(isnan(hk_summa)) = 0;
    if (TOF && randoms_correction)
        rInd = rand == 0;
        P_Sino = Sino(Sino > 0 & repmat(rInd, nBins,1));
        apu = D + repmat(rand, nBins, 1);
        apu = sum(Sino .* log(apu) - apu);
        hk_summa = bsxfun(@minus, sum(hk_summa), hk_summa(Sino > 0 & repmat(rInd, nBins, 1)));
    else
        if randoms_correction == 1
            Iind = (Sino > 0 & rand == 0);
            if (sum(Iind) == 0)
                epsilon_mramla = 1e8;
                return;
            end
            P_Sino = Sino(Iind);
            apu = D + rand;
            apu = sum(Sino .* log(apu) - apu);
        else
            Iind = (Sino > 0);
            P_Sino = Sino(Iind);
            apu = sum(Sino .* log(D) - D);
        end
        hk_summa = bsxfun(@minus, sum(hk_summa), hk_summa(Iind));
    end
    epsilon = min(P_Sino, exp(bsxfun(@minus, apu, hk_summa) ./ P_Sino));
    epsilon_mramla = min(epsilon);
end
if epsilon_mramla <= 0
    epsilon_mramla = epps;
end