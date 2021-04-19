function var = recNames(varargin)
%RECNAMES Obtain the OMEGA reconstruction variable names 
%   0 returns all variable names in uppercase
%   1 returns all prior names
%   2 returns all MAP names
%   otherwise the variable names are returned as in the options-struct
%   Add your own algorithm/prior name to the lists below
varMAP = {'OSL_MLEM';'OSL_OSEM';'BSREM';'MBSREM';'ROSEM_MAP';'RBI_OSL';'COSEM_OSL';'PKMA'};
nMAPMLEM = 1;
nMLEM = 1;
varMLEM = {'mlem'};
varML = {varMLEM{1:nMLEM};varMAP{1:nMAPMLEM}};
OS = {'osem';'mramla';'ramla';'rosem';'rbi';'drama';'cosem';'ecosem';'acosem'};
varOS = [OS;varMAP(nMAPMLEM + 1:end)];
varPrior = {'MRP';'quad';'Huber';'L';'FMH';'weighted_mean';'TV';'AD';'APLS';'TGV';'NLM';'custom'};
if ~isempty(varargin) && varargin{1} == 0
    var = upper([{upper(varML{1})};OS;varMAP(:)]);
elseif ~isempty(varargin) && varargin{1} == 1
    var = varPrior;
elseif ~isempty(varargin) && varargin{1} == 2
    var = varMAP;
elseif ~isempty(varargin) && varargin{1} == 3
    var = varML;
elseif ~isempty(varargin) && varargin{1} == 4
    var = varOS;
elseif ~isempty(varargin) && varargin{1} == 5
    var = varMLEM;
elseif ~isempty(varargin) && varargin{1} == 6
    var = OS;
else
    var = [varML(1);OS(:);varMAP(:)];
end
end