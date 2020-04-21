function [algo_char] = algorithms_char()
%ALGORITHMS_CHAR Returns the available reconstruction algorithms as a char
%array

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ville-Veikko Wettenhovi
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License; or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful;
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not; see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

algo_char = {'MLEM';'OSEM';'MRAMLA';'RAMLA';'ROSEM'; ... % 1 - 5
    'RBI'; 'DRAMA'; 'COSEM';'ECOSEM';'ACOSEM';... % 6 - 10
    'MRP-OSL-OSEM';'MRP-OSL-MLEM';'MRP-BSREM';'MRP-MBSREM';'MRP-ROSEM';... % 11 - 15
    'MRP-RBI';'MRP-OSL-COSEM';'QP (OSL-OSEM)';'QP (OSL-MLEM)';'QP (BSREM)';... % 16 - 20
    'QP (MBSREM)';'QP (ROSEM)';'QP (OSL-RBI)'; 'QP (OSL-COSEM)';'HP (OSL-OSEM)';... % 21 - 25
    'HP (OSL-MLEM)';'HP (BSREM)';'HP (MBSREM)';'HP (ROSEM)';'HP (OSL-RBI)'; ... % 26 - 30
    'HP (OSL-COSEM)';'L-filter (OSL-OSEM)'; 'L-filter (OSL-MLEM)'; 'L-filter (BSREM)'; 'L-filter (MBSREM)'; ... % 31 - 35
    'L-filter (ROSEM)';'L-filter (OSL-RBI)'; 'L-filter (OSL-COSEM)'; 'FMH (OSL-OSEM)'; 'FMH (OSL-MLEM)';... % 36 - 40
    'FMH (BSREM)'; 'FMH (MBSREM)'; 'FMH (ROSEM)'; 'FMH (OSL-RBI)';'FMH (OSL-COSEM)';... % 41 - 45
    'Weighted mean (OSL-OSEM)'; 'Weighted mean (OSL-MLEM)'; 'Weighted mean (BSREM)'; 'Weighted mean (MBSREM)'; 'Weighted mean (ROSEM)';... % 46 - 50
    'Weighted mean (OSL-RBI)'; 'Weighted mean (OSL-COSEM)'; 'Total variation (OSL-OSEM)';'Total variation (OSL-MLEM)';'Total variation (BSREM)';... % 51 - 55
    'Total variation (MBSREM)';'Total variation (ROSEM)';'Total variation (OSL-RBI)';'Total variation (OSL-COSEM)';'Anisotropic Diffusion (OSL-OSEM)';... % 56 - 60
    'Anisotropic Diffusion (OSL-MLEM)'; 'Anisotropic Diffusion (BSREM)';'Anisotropic Diffusion (MBSREM)';'Anisotropic Diffusion (ROSEM)';'Anisotropic Diffusion (OSL-RBI)';... % 61 - 65
    'Anisotropic Diffusion (OSL-COSEM)';'APLS (OSL-OSEM)'; 'APLS (OSL-MLEM)';'APLS (BSREM)';'APLS (MBSREM)';... % 66 - 70
    'APLS (ROSEM)';'APLS (OSL-RBI)';'APLS (OSL-COSEM)';'TGV (OSL-OSEM)'; 'TGV (OSL-MLEM)';... % 71 - 75
    'TGV (BSREM)';'TGV (MBSREM)';'TGV (ROSEM)';'TGV (OSL-RBI)';'TGV (OSL-COSEM)';... % 76 - 80
    'NLM (OSL-OSEM)'; 'NLM (OSL-MLEM)';'NLM (BSREM)';'NLM (MBSREM)';'NLM (ROSEM)';... % 81 - 85
    'NLM (OSL-RBI)';'NLM (OSL-COSEM)';'Custom prior (OSL-OSEM)'; 'Custom prior (OSL-MLEM)';'Custom prior (BSREM)';...% 86 - 90
    'Custom prior (MBSREM)';'Custom prior (ROSEM)';'Custom prior (OSL-RBI)';'Custom prior (OSL-COSEM)';'Image properties'}; % 91 - 95
end

