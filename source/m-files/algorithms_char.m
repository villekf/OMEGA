function [algo_char] = algorithms_char(varargin)
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

varNonMAP = [recNames(6)];
varMAP = [recNames(2)];
varAll = recNames(11);
% varMAP = strrep(varMAP,'_','-');
varPriorNames = recNames(10);

algo_char = cell(numel(varNonMAP) + numel(varPriorNames) * numel(varMAP) + 1,1);
uu = 1;
ii = numel(varNonMAP) + 1;
for kk = 1 : numel(algo_char) - 1
    if kk <= numel(varNonMAP)
        algo_char{kk} = varAll{kk};
    else
        algo_char{kk} = [varPriorNames{uu} ' (' varAll{ii} ')'];
        if ii == numel(varNonMAP) + numel(varMAP)
            uu = uu + 1;
            ii = numel(varNonMAP) + 1;
        else
            ii = ii + 1;
        end
    end
end
algo_char{end} = 'FDK';
% algo_char{end} = 'Image properties';
dispi = [];
if nargin >= 1
    for kk = 1 : numel(algo_char) - 1
        dispi = [dispi, [num2str(kk) ' = ' algo_char{kk}]];
        if kk < numel(algo_char) - 1
            dispi = [dispi, ', '];
        end
    end
    disp(dispi)
end
% 
% algo_char = {'MLEM';'OSEM';'MRAMLA';'RAMLA';'ROSEM'; ... % 1 - 5
%     'RBI'; 'DRAMA'; 'COSEM';'ECOSEM';'ACOSEM';... % 6 - 10
%     'MRP-OSL-MLEM';'MRP-OSL-OSEM';'MRP-BSREM';'MRP-MBSREM';'MRP-ROSEM';... % 11 - 15
%     'MRP-RBI';'MRP-OSL-COSEM';'QP (OSL-MLEM)';'QP (OSL-OSEM)';'QP (BSREM)';... % 16 - 20
%     'QP (MBSREM)';'QP (ROSEM)';'QP (OSL-RBI)'; 'QP (OSL-COSEM)';'HP (OSL-MLEM)';'HP (OSL-OSEM)';... % 21 - 25
%     'HP (BSREM)';'HP (MBSREM)';'HP (ROSEM)';'HP (OSL-RBI)'; ... % 26 - 30
%     'HP (OSL-COSEM)'; 'L-filter (OSL-MLEM)';'L-filter (OSL-OSEM)'; 'L-filter (BSREM)'; 'L-filter (MBSREM)'; ... % 31 - 35
%     'L-filter (ROSEM)';'L-filter (OSL-RBI)'; 'L-filter (OSL-COSEM)'; 'FMH (OSL-MLEM)'; 'FMH (OSL-OSEM)';... % 36 - 40
%     'FMH (BSREM)'; 'FMH (MBSREM)'; 'FMH (ROSEM)'; 'FMH (OSL-RBI)';'FMH (OSL-COSEM)';... % 41 - 45
%     'Weighted mean (OSL-MLEM)'; 'Weighted mean (OSL-OSEM)'; 'Weighted mean (BSREM)'; 'Weighted mean (MBSREM)'; 'Weighted mean (ROSEM)';... % 46 - 50
%     'Weighted mean (OSL-RBI)'; 'Weighted mean (OSL-COSEM)';'Total variation (OSL-MLEM)'; 'Total variation (OSL-OSEM)';'Total variation (BSREM)';... % 51 - 55
%     'Total variation (MBSREM)';'Total variation (ROSEM)';'Total variation (OSL-RBI)';'Total variation (OSL-COSEM)';... % 56 - 60
%     'Anisotropic Diffusion (OSL-MLEM)'; 'Anisotropic Diffusion (OSL-OSEM)';'Anisotropic Diffusion (BSREM)';'Anisotropic Diffusion (MBSREM)';'Anisotropic Diffusion (ROSEM)';'Anisotropic Diffusion (OSL-RBI)';... % 61 - 65
%     'Anisotropic Diffusion (OSL-COSEM)'; 'APLS (OSL-MLEM)';'APLS (OSL-OSEM)';'APLS (BSREM)';'APLS (MBSREM)';... % 66 - 70
%     'APLS (ROSEM)';'APLS (OSL-RBI)';'APLS (OSL-COSEM)'; 'TGV (OSL-MLEM)';'TGV (OSL-OSEM)';... % 71 - 75
%     'TGV (BSREM)';'TGV (MBSREM)';'TGV (ROSEM)';'TGV (OSL-RBI)';'TGV (OSL-COSEM)';... % 76 - 80
%     'NLM (OSL-MLEM)';'NLM (OSL-OSEM)'; 'NLM (BSREM)';'NLM (MBSREM)';'NLM (ROSEM)';... % 81 - 85
%     'NLM (OSL-RBI)';'NLM (OSL-COSEM)'; 'Custom prior (OSL-MLEM)';'Custom prior (OSL-OSEM)';'Custom prior (BSREM)';...% 86 - 90
%     'Custom prior (MBSREM)';'Custom prior (ROSEM)';'Custom prior (OSL-RBI)';'Custom prior (OSL-COSEM)';'Image properties'}; % 91 - 95
% %'';'';'';'';'';'';'';'';'';'';'';'';
end

