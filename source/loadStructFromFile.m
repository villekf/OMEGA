function [result, varargout] = loadStructFromFile(fileName, environmentName, varargin)
tmp = load(fileName, environmentName);
result = tmp.(environmentName);
if nargin > 2
    try
        tmp = load(fileName, varargin{1});
        varargout{1} = tmp.(varargin{1});
    catch
        varargout{1} = [];
    end
end