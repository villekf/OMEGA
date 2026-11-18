function options = checkAlgorithmSelection(options)
%% Algorithm and prior selection checks
% Checks that the selected reconstruction algorithm(s) and prior(s) form an
% allowed combination (e.g. a MAP algorithm requires a prior). Some
% selections are corrected with a warning (the prior is disabled if no MAP
% algorithm has been selected and beta is set to zero if a MAP algorithm has
% no prior). The checks are skipped when only the data are loaded
% (options.only_sinos) or when the normalization coefficients are computed
% (options.compute_normalization).
%
% INPUT:
%   options = The options-struct (with the default values set, see
%   setMissingValues)
%
% OUTPUT:
%   options = The (possibly corrected) options-struct
%
% See also OMEGA_error_check, printReconstructionSettings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2026 Ville-Veikko Wettenhovi
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
if options.compute_normalization || options.only_sinos
    return
end
MAP = checkAlgorithmsPriors(options, 2);
PRIOR = checkAlgorithmsPriors(options, 1);
PRIOR_summa = checkAlgorithmsPriors(options, 1, 0);
OS = checkAlgorithmsPriors(options, 6);

reko = {};
varNonMAP = [recNames(6)];
for kk = 1 : numel(varNonMAP)
    if options.OSEM && options.subsets == 1
        reko = [reko;{'MLEM'}];
        break
    elseif options.(varNonMAP{kk})
        ch = strrep(varNonMAP{kk},'_','-');
        reko = [reko;ch];
        break
    end
end
varMAP = recNames(2);
for kk = 1 : numel(varMAP)
    if options.(varMAP{kk}) && strcmp(varMAP{kk},'OSL_OSEM') && PRIOR && options.subsets == 1
        reko = [reko;{'MLEM-OSL'}];
    elseif options.(varMAP{kk}) && strcmp(varMAP{kk},'OSL_COSEM')
        if options.OSL_COSEM == 1 && PRIOR
            reko = [reko;{'ACOSEM-OSL'}];
        elseif options.OSL_COSEM == 1 && ~PRIOR
            error('ACOSEM-OSL selected, but no prior has been selected.')
        elseif options.OSL_COSEM == 2 && PRIOR
            reko = [reko;{'COSEM-OSL'}];
        elseif options.OSL_COSEM == 2 && ~PRIOR
            error('COSEM-OSL selected, but no prior has been selected.')
        elseif options.OSL_COSEM > 2 || options.OSL_COSEM < 0
            error('Unsupported COSEM-OSL method selected!')
        end
        if PRIOR_summa > 1
            error('Only one prior can be used at a time with COSEM-OSL')
        end
    elseif options.(varMAP{kk}) && PRIOR
        ch = strrep(varMAP{kk},'_','-');
        reko = [reko;ch];
    elseif options.(varMAP{kk}) && ~PRIOR
        ch = strrep(varMAP{kk},'_','-');
        warning([ch ' selected, but no prior has been selected.'])
        reko = [reko;ch];
        options.beta = 0;
    end
end
if isempty(reko)
    error('No reconstruction method selected!')
end
priori = {};

varPrior = recNames(1);
varPriorName = recNames(10);
for kk = 1 : numel(varPrior)
    if options.(varPrior{kk}) && MAP
        if strcmp(varPrior{kk},'TV')
            if options.TV_use_anatomical
                if options.TVtype == 1
                    priori = [priori;{'Anatomically weighted TV'}];
                elseif options.TVtype == 2
                    priori = [priori;{'Joint TV'}];
                elseif options.TVtype == 3
                    priori = [priori;{'Weighted joint TV'}];
                else
                    error('Unsupported TV type selected.')
                end
            else
                if options.TVtype == 1 || options.TVtype == 2
                    priori = [priori;{'TV'}];
                elseif options.TVtype == 3
                    priori = [priori;{'Weighted TV'}];
                elseif options.TVtype == 4
                    priori = [priori;{'SATV'}];
                else
                    error('Unsupported TV type selected.')
                end
            end
        elseif strcmp(varPrior{kk},'weighted_mean')
            if options.mean_type == 1 || options.mean_type == 4
                priori = [priori;{'Weighted (arithmetic) mean'}];
            elseif options.mean_type == 2 || options.mean_type == 5
                priori = [priori;{'Weighted (harmonic) mean'}];
            elseif options.mean_type == 3 || options.mean_type == 6
                priori = [priori;{'Weighted (geometric) mean'}];
            else
                error('Unsupported mean type selected.')
            end
        elseif strcmp(varPrior{kk},'AD')
            if options.FluxType > 2 || options.FluxType < 1
                error('FluxType has to be either 1 or 2.')
            end
            if (options.DiffusionType > 2 || options.DiffusionType < 1) && options.implementation == 2
                error('DiffusionType has to be either 1 or 2.')
            end
            priori = [priori;{'AD-MRP'}];
        else
            ch = strrep(varPriorName{kk},'_','-');
            priori = [priori;ch];
        end
    elseif options.(varPrior{kk}) && ~MAP
        ch = strrep(varPriorName{kk},'_','-');
        warning([ch ' selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.'])
        options.(varPrior{kk}) = false;
    end
end
if ~OS && ~MAP && ~options.MLEM && ~(options.only_sinos || options.compute_normalization)
    error('No reconstruction algorithm selected.')
end
