function output = checkAlgorithmsPriors(options, varargin)
%UNTITLED Summary of this function goes here
%   0 returns all variable names in uppercase
%   1 returns all prior names
%   2 returns all MAP names
%   3 returns all MLEM-based names
%   4 returns all OS-based names
%   5 returns all MLEM-based non-MAP names
%   6 returns all OS-based non-MAP names
%   7 returns all MLEM-based MAP names
%   8 returns all OS-based MAP names
%   otherwise the variable names are returned as in the options-struct
if nargin >= 2 && ~isempty(varargin{1})
    type = varargin{1};
else
    type = 0;
end
var = recNames(type);
ll = 0;
for kk = 1 : numel(var)
    ll = ll + min(options.(var{kk}),1);
end
if nargin >= 3
    output = ll;
else
    output = ll > 0;
end
end

