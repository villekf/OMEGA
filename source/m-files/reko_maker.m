function rekot = reko_maker(options)
% Form the boolean vector that tells what algorithms and/or priors are used
% Indices: 
% 1 = MLEM, 2 = OSEM, 3 = MRAMLA, 4 = RAMLA, 5 = ROSEM, 6 = RBI, 7 = DRAMA,
% 8 = COSEM, 9 = ECOSEM, 10 = ACOSEM, 11 = OSL-OSEM MRP, 12 = OSL-MLEM MRP,
% 13 = BSREM MRP, 14 = MBSREM MRP, 15 = ROSEM-MAP MRP, 16 = RBI-OSL MRP, 17
% - 22 = Quadratic, 23 - 28 = L-filter, 29 - 34 = FMH, 35 - 40 = Weighted
% mean, 41 - 46 = TV, 47 - 52 = AD, 53 - 58 = APLS
% Total ML-methods + Total MAP-methods * Total number of priors + extra
% space for image properties
%
% 
% 
% Copyright (C) 2020 Ville-Veikko Wettenhovi
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
% ML-methods = 10, MAP-methods = 8, Priors = 11 + custom
if ~isfield(options,'custom')
    options.custom = false;
end

varML = recNames();
% var = variableNames(0);
varPrior = recNames(1);

prior = numel(varPrior);
MAP = numel(recNames(2));
ML = numel([recNames(6)]);

rekot = false(ML+MAP*prior + 1,1);
for gg = 1 : numel(varML) %ML
    if options.(varML{gg})
        rekot(gg) = true;
    end
end
% ll = 1;
% for gg = gg + 1 : MAP : numel(rekot) - 1 - prior
%     tt = ML + 1;
%     for uu = gg : gg + MAP - 1
%         if options.(varML{tt}) && options.(varPrior{ll})
%             rekot(uu) = true;
%         end
%         tt = tt + 1;
%     end
%     ll = ll + 1;
% end
% tt = ML + 1;
% for gg = uu + 1 : numel(rekot) - 1
%     if options.(varML{tt}) && options.(varPrior{end})
%         rekot(gg) = true;
%     end
%     tt = tt + 1;
% end