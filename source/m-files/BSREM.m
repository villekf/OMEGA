function output = BSREM(im, rhs, lam, iter, varargin)
%BSREM BSREM/RAMLA subset update (BSREM in algorithms.h)
%   output = BSREM(im, rhs, lam, iter)
%   output = BSREM(im, rhs, lam, iter, epps, enforcePositivity)
%   With enforcePositivity the result is limited to at least epps (as in
%   the GPU PoissonUpdate kernel)
output = im + lam(iter) .* im .* rhs;
if nargin >= 6 && varargin{2}
    output = max(varargin{1}, output);
end
end
