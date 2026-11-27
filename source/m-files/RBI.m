function output = RBI(im, Sens, rhs, D, varargin)
%RBI Rescaled block iterative EM (RBI) and one-step-late RBI update (RBI in algorithms.h)
%   output = RBI(im, Sens, rhs, D)
%   output = RBI(im, Sens, rhs, D, beta, dU, epps)
%   With beta ~= 0 the denominator (D + dU) is limited to at least 10 % of D
%   and at least epps, and the floored denominator replaces D + dU in all
%   places: rhs + D - denom equals rhs - dU and Sens - D + denom equals
%   Sens + dU when the floor is not active
if nargin <= 4
    beta = 0;
else
    beta = varargin{1};
    dU = varargin{2};
    if nargin >= 7
        epps = varargin{3};
    else
        epps = 1e-8;
    end
end
if (beta == 0)
    Summa = 1 ./ max(Sens ./ D);
    output = im + Summa * (im ./ D) .* (rhs);
else
    denom = max(D + dU, max(0.1 * D, epps));
    Summa = 1 ./ max((Sens - D + denom) ./ denom);
    output = im + Summa .* (im ./ denom) .* (rhs + D - denom);
end
end
