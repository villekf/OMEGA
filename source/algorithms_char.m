function [algo_char] = algorithms_char()
%ALGORITHMS_CHAR Returns the available reconstruction algorithms as a char
%array

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

algo_char = ["MLEM","OSEM","MRAMLA","RAMLA","ROSEM", ... % 1 - 5
    "RBI", "DRAMA", "COSEM","ECOSEM","ACOSEM",... % 6 - 10
    "MRP-OSL-OSEM","MRP-OSL-MLEM","MRP-MBSREM","MRP-BSREM","MRP-ROSEM",... % 11 - 15
    "MRP-RBI","MRP-OSL-COSEM","QP (OSL-OSEM)","QP (OSL-MLEM)","QP (MBSREM)",... % 16 - 20
    "QP (BSREM)","QP (ROSEM)","QP (RBI)", "QP (OSL-COSEM)","L-filter (OSL-OSEM)", ... % 21 - 25
    "L-filter (OSL-MLEM)", "L-filter (MBSREM)", "L-filter (BSREM)", "L-filter (ROSEM)","L-filter (RBI)", ... % 26 - 30
    "L-filter (OSL-COSEM)", "FMH (OSL-OSEM)", "FMH (OSL-MLEM)","FMH (MBSREM)", "FMH (BSREM)", ... % 31 - 35
    "FMH (ROSEM)", "FMH (RBI)","FMH (OSL-COSEM)","Weighted mean (OSL-OSEM)", "Weighted mean (OSL-MLEM)", ... % 36 - 40
    "Weighted mean (MBSREM)", "Weighted mean (BSREM)", "Weighted mean (ROSEM)","Weighted mean (RBI)", "Weighted mean (OSL-COSEM)", ... % 41 - 45
    "Total variation (OSL-OSEM)","Total variation (OSL-MLEM)","Total variation (MBSREM)","Total variation (BSREM)","Total variation (ROSEM)",... % 46 - 50
    "Total variation (RBI)","Total variation (OSL-COSEM)","Anisotropic Diffusion (OSL-OSEM)","Anisotropic Diffusion (OSL-MLEM)", "Anisotropic Diffusion (MBSREM)",... % 51 - 55
    "Anisotropic Diffusion (BSREM)","Anisotropic Diffusion (ROSEM)","Anisotropic Diffusion (RBI)","Anisotropic Diffusion (OSL-COSEM)","APLS (OSL-OSEM)", ... % 56 - 60
    "APLS (OSL-MLEM)","APLS (MBSREM)","APLS (BSREM)","APLS (ROSEM)","APLS (RBI)",... % 61 - 65
    "APLS (OSL-COSEM)","TGV (OSL-OSEM)", "TGV (OSL-MLEM)","TGV (MBSREM)","TGV (BSREM)",... % 66 - 70
    "TGV (ROSEM)","TGV (RBI)","TGV (OSL-COSEM)","NLM (OSL-OSEM)", "NLM (OSL-MLEM)",... % 71 - 75
    "NLM (MBSREM)","NLM (BSREM)","NLM (ROSEM)","NLM (RBI)","NLM (OSL-COSEM)",... % 76 - 80
    "Custom prior (OSL-MLEM)", "Custom prior (OSL-OSEM)","Custom prior (MBSREM)","Custom prior (BSREM)","Custom prior (ROSEM)",...% 81 - 85
    "Custom prior (RBI)","Custom prior (OSL-COSEM)"];
end

