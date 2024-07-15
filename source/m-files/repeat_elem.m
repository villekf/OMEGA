function output = repeat_elem(A, R, varargin)
%REPEAT_ELEM Repeats the elements of vector/matrix A by the values in vector R
%
% Example:
%   B = repeat_elem(A,R)
% INPUTS: 
%   A = Vector/matrix to be repeated
%   R = Scalar/vector containing the amounts that each of the elements will be
%   repeated
%   d = Dimension where repeating is performed (zero-based indexing, e.g. 0
%   is the first dimension)
% OUTPUTS:
%   B = Vector/matrix A with repeated elements
%
% See also repelem

% Source: 
% https://se.mathworks.com/matlabcentral/answers/346102-repelem-substitute-for-older-matlab-version

if nargin > 3
    error('Invalid number of input arguments')
end
numvarargs = length(varargin);
optargs = {0};
optargs(1:numvarargs) = varargin;
d = optargs{1};

if size(R,1) == 1
    R = repmat(R,size(A,1),size(A,2));
    if d == 0
        output = cell2mat(arrayfun(@(a,r)repmat(a,r,1),A,R,'uni',0));
    elseif d == 1
        output = cell2mat(arrayfun(@(a,r)repmat(a,1,r),A,R,'uni',0));
    else
        output = cell2mat(arrayfun(@(a,r)repmat(a,1,1,r),A,R,'uni',0));
    end
else
    if d == 0
        output = cell2mat(arrayfun(@(a,r)repmat(a,r,1),A,R,'uni',0));
    elseif d == 1
        output = cell2mat(arrayfun(@(a,r)repmat(a,1,r),A,R,'uni',0));
    else
        output = cell2mat(arrayfun(@(a,r)repmat(a,1,1,r),A,R,'uni',0));
    end
end
end