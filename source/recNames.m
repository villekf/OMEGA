function varList = recNames(varargin)
%RECNAMES Obtain the OMEGA reconstruction variable names 
%   0 returns all variable names in uppercase
%   1 returns all prior names
%   2 returns all MAP names
%   3 returns all MLEM-based names
%   4 returns all OS-based names
%   5 returns all MLEM-based non-MAP names
%   6 returns all OS-based non-MAP names
%   7 returns all MLEM-based MAP names
%   8 returns all OS-based MAP names
%   9 returns all algorithms not supported by implementation 3
%   10 returns all full prior names, for display purposes only
%   otherwise the variable names are returned as in the options-struct
%   Add your own algorithm/prior name to the lists below

% Add Any MAP algorithm here
% Add non-subset MAP algorithm before 'OSL_OSEM'
varMAP = {'OSL_MLEM';'OSL_OSEM';'BSREM';'MBSREM';'ROSEM_MAP';'OSL_RBI';'OSL_COSEM';'PKMA'};
% Increment this if the added algorithm is MLEM (no subsets) type MAP
% algorithm
nMAPMLEM = 1;
% Increment this if the added algorithm is MLEM (no subsets) type algorithm
% (non-MAP)
nMLEM = 1;
% Include the (non-MAP/prior-based) MLEM algorithm here
varMLEM = {'MLEM'};
% This  will be automatically filled
varML = {varMLEM{1:nMLEM};varMAP{1:nMAPMLEM}};
% This  will be automatically filled
varMAPML = {varMAP{1:nMAPMLEM}};
% Add an OS-based (non-MAP/prior-based) algorithm to this list)
OS = {'OSEM';'MRAMLA';'RAMLA';'ROSEM';'RBI';'DRAMA';'COSEM';'ECOSEM';'ACOSEM'};
% This  will be automatically filled
varOS = [OS;varMAP(nMAPMLEM + 1:end)];
% This  will be automatically filled
varMAPOS = [varMAP(nMAPMLEM + 1:end)];
% Add your new prior to this list (before 'custom')
varPrior = {'MRP';'quad';'Huber';'L';'FMH';'weighted_mean';'TV';'AD';'APLS';'TGV';'NLM';'custom'};
% Add your prior to this list also (this is used for display purposes only
% and can be unabbreviated and contain any characters):
varPriorFull = {'Median Root';'Quadratic';'Huber';'L-filter';'FIR Median Hybrid';'Weighted mean';'Total Variation';'Anisotropic Diffusion';...
    'Asymmetric Parallel Level Sets';'Total Generalized Variation';'Non-Local Means';'Custom'};
if ~isempty(varargin) && varargin{1} == 0
    varList = upper([upper(varMLEM(:));OS(:);varMAP(:)]);
elseif ~isempty(varargin) && varargin{1} == 1
    varList = varPrior;
elseif ~isempty(varargin) && varargin{1} == 2
    varList = varMAP;
elseif ~isempty(varargin) && varargin{1} == 3
    varList = varML;
elseif ~isempty(varargin) && varargin{1} == 4
    varList = varOS;
elseif ~isempty(varargin) && varargin{1} == 5
    varList = varMLEM;
elseif ~isempty(varargin) && varargin{1} == 6
    varList = OS;
elseif ~isempty(varargin) && varargin{1} == 7
    varList = varMAPML;
elseif ~isempty(varargin) && varargin{1} == 8
    varList = varMAPOS;
elseif ~isempty(varargin) && varargin{1} == 9
    varList = [varMLEM(:);OS(:);varMAP(:)];
    varList = varList(~strcmp(varList,'MLEM'));
    varList = varList(~strcmp(varList,'OSEM'));
elseif ~isempty(varargin) && varargin{1} == 10
    varList = varPriorFull;
else
    varList = [varMLEM(:);OS(:);varMAP(:)];
end
end