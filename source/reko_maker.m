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
% ML-methods = 10, MAP-methods = 7, Priors = 11 + custom
if ~isfield(options,'custom')
    options.custom = false;
end
rekot = false(10+7*12 + 1,1);
gg = 1;
if options.MLEM
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSEM
    rekot(gg) = true;
end
gg = gg + 1;
if options.MRAMLA
    rekot(gg) = true;
end
gg = gg + 1;
if options.RAMLA
    rekot(gg) = true;
end
gg = gg + 1;
if options.ROSEM
    rekot(gg) = true;
end
gg = gg + 1;
if options.RBI
    rekot(gg) = true;
end
gg = gg + 1;
if options.DRAMA
    rekot(gg) = true;
end
gg = gg + 1;
if options.COSEM
    rekot(gg) = true;
end
gg = gg + 1;
if options.ECOSEM
    rekot(gg) = true;
end
gg = gg + 1;
if options.ACOSEM
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_OSEM && options.MRP
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_MLEM && options.MRP
    rekot(gg) = true;
end
gg = gg + 1;
if options.BSREM && options.MRP
    rekot(gg) = true;
end
gg = gg + 1;
if options.MBSREM && options.MRP
    rekot(gg) = true;
end
gg = gg + 1;
if options.ROSEM_MAP && options.MRP
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_RBI && options.MRP
    rekot(gg) = true;
end
gg = gg + 1;
if any(options.OSL_COSEM) && options.MRP
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_OSEM && options.quad
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_MLEM && options.quad
    rekot(gg) = true;
end
gg = gg + 1;
if options.BSREM && options.quad
    rekot(gg) = true;
end
gg = gg + 1;
if options.MBSREM && options.quad
    rekot(gg) = true;
end
gg = gg + 1;
if options.ROSEM_MAP && options.quad
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_RBI && options.quad
    rekot(gg) = true;
end
gg = gg + 1;
if any(options.OSL_COSEM) && options.quad
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_OSEM && options.Huber
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_MLEM && options.Huber
    rekot(gg) = true;
end
gg = gg + 1;
if options.BSREM && options.Huber
    rekot(gg) = true;
end
gg = gg + 1;
if options.MBSREM && options.Huber
    rekot(gg) = true;
end
gg = gg + 1;
if options.ROSEM_MAP && options.Huber
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_RBI && options.Huber
    rekot(gg) = true;
end
gg = gg + 1;
if any(options.OSL_COSEM) && options.Huber
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_OSEM && options.L
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_MLEM && options.L
    rekot(gg) = true;
end
gg = gg + 1;
if options.BSREM && options.L
    rekot(gg) = true;
end
gg = gg + 1;
if options.MBSREM && options.L
    rekot(gg) = true;
end
gg = gg + 1;
if options.ROSEM_MAP && options.L
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_RBI && options.L
    rekot(gg) = true;
end
gg = gg + 1;
if any(options.OSL_COSEM) && options.L
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_OSEM && options.FMH
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_MLEM && options.FMH
    rekot(gg) = true;
end
gg = gg + 1;
if options.BSREM && options.FMH
    rekot(gg) = true;
end
gg = gg + 1;
if options.MBSREM && options.FMH
    rekot(gg) = true;
end
gg = gg + 1;
if options.ROSEM_MAP && options.FMH
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_RBI && options.FMH
    rekot(gg) = true;
end
gg = gg + 1;
if any(options.OSL_COSEM) && options.FMH
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_OSEM && options.weighted_mean
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_MLEM && options.weighted_mean
    rekot(gg) = true;
end
gg = gg + 1;
if options.BSREM && options.weighted_mean
    rekot(gg) = true;
end
gg = gg + 1;
if options.MBSREM && options.weighted_mean
    rekot(gg) = true;
end
gg = gg + 1;
if options.ROSEM_MAP && options.weighted_mean
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_RBI && options.weighted_mean
    rekot(gg) = true;
end
gg = gg + 1;
if any(options.OSL_COSEM) && options.weighted_mean
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_OSEM && options.TV
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_MLEM && options.TV
    rekot(gg) = true;
end
gg = gg + 1;
if options.BSREM && options.TV
    rekot(gg) = true;
end
gg = gg + 1;
if options.MBSREM && options.TV
    rekot(gg) = true;
end
gg = gg + 1;
if options.ROSEM_MAP && options.TV
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_RBI && options.TV
    rekot(gg) = true;
end
gg = gg + 1;
if any(options.OSL_COSEM) && options.TV
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_OSEM && options.AD
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_MLEM && options.AD
    rekot(gg) = true;
end
gg = gg + 1;
if options.BSREM && options.AD
    rekot(gg) = true;
end
gg = gg + 1;
if options.MBSREM && options.AD
    rekot(gg) = true;
end
gg = gg + 1;
if options.ROSEM_MAP && options.AD
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_RBI && options.AD
    rekot(gg) = true;
end
gg = gg + 1;
if any(options.OSL_COSEM) && options.AD
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_OSEM && options.APLS
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_MLEM && options.APLS
    rekot(gg) = true;
end
gg = gg + 1;
if options.BSREM && options.APLS
    rekot(gg) = true;
end
gg = gg + 1;
if options.MBSREM && options.APLS
    rekot(gg) = true;
end
gg = gg + 1;
if options.ROSEM_MAP && options.APLS
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_RBI && options.APLS
    rekot(gg) = true;
end
gg = gg + 1;
if any(options.OSL_COSEM) && options.APLS
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_OSEM && options.TGV
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_MLEM && options.TGV
    rekot(gg) = true;
end
gg = gg + 1;
if options.BSREM && options.TGV
    rekot(gg) = true;
end
gg = gg + 1;
if options.MBSREM && options.TGV
    rekot(gg) = true;
end
gg = gg + 1;
if options.ROSEM_MAP && options.TGV
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_RBI && options.TGV
    rekot(gg) = true;
end
gg = gg + 1;
if any(options.OSL_COSEM) && options.TGV
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_OSEM && options.NLM
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_MLEM && options.NLM
    rekot(gg) = true;
end
gg = gg + 1;
if options.BSREM && options.NLM
    rekot(gg) = true;
end
gg = gg + 1;
if options.MBSREM && options.NLM
    rekot(gg) = true;
end
gg = gg + 1;
if options.ROSEM_MAP && options.NLM
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_RBI && options.NLM
    rekot(gg) = true;
end
gg = gg + 1;
if any(options.OSL_COSEM) && options.NLM
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_OSEM && options.custom
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_MLEM && options.custom
    rekot(gg) = true;
end
gg = gg + 1;
if options.BSREM && options.custom
    rekot(gg) = true;
end
gg = gg + 1;
if options.MBSREM && options.custom
    rekot(gg) = true;
end
gg = gg + 1;
if options.ROSEM_MAP && options.custom
    rekot(gg) = true;
end
gg = gg + 1;
if options.OSL_RBI && options.custom
    rekot(gg) = true;
end
gg = gg + 1;
if any(options.OSL_COSEM) && options.custom
    rekot(gg) = true;
end