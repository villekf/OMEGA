function varList = recNames(varargin)
%RECNAMES Obtain the OMEGA reconstruction variable names 
%   0 returns all variable names in uppercase
%   1 returns all prior names
%   2 returns all MAP names
%   6 returns all non-MAP names
%   9 returns all algorithms not supported by implementation 3
%   10 returns all full prior names, for display purposes only
%   otherwise the variable names are returned as in the options-struct
%   Add your own algorithm/prior name to the lists below

% Add Any MAP algorithm here (algorithms that can use gradient-based
% priors)
varMAP = {'OSL_OSEM';'BSREM';'MBSREM';'ROSEM_MAP';'OSL_RBI';'OSL_COSEM';'PKMA'};
% Add any non-MAP/prior-based algorithm to this list (or algorithms with
% fixed prior)
OS = {'OSEM';'MRAMLA';'RAMLA';'ROSEM';'RBI';'DRAMA';'COSEM';'ECOSEM';'ACOSEM';'LSQR';'CPLS';'CGLS';'CPTV'};
% Add your new prior to this list (before 'custom')
varPrior = {'MRP';'quad';'Huber';'L';'FMH';'weighted_mean';'TV';'AD';'APLS';'TGV';'NLM';'RDP';'custom'};
% Add your prior to this list also (this is used for display purposes only
% and can be unabbreviated and contain any characters):
varPriorFull = {'Median Root';'Quadratic';'Huber';'L-filter';'FIR Median Hybrid';'Weighted mean';'Total Variation';'Anisotropic Diffusion';...
    'Asymmetric Parallel Level Sets';'Total Generalized Variation';'Non-Local Means';'Relative Difference';'Custom'};
if ~isempty(varargin) && varargin{1} == 0
    varList = upper([OS(:);varMAP(:)]);
elseif ~isempty(varargin) && varargin{1} == 1
    varList = varPrior;
elseif ~isempty(varargin) && varargin{1} == 2
    varList = varMAP;
elseif ~isempty(varargin) && varargin{1} == 6
    varList = OS;
elseif ~isempty(varargin) && varargin{1} == 9
    varList = [OS(:);varMAP(:)];
    varList = varList(~strcmp(varList,'OSEM'));
elseif ~isempty(varargin) && varargin{1} == 10
    varList = varPriorFull;
else
    varList = [OS(:);varMAP(:)];
end
end