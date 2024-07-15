function output = RBI(im, Sens, rhs, D, varargin)
if nargin == 4
    beta = 0;
else
    beta = varargin{1};
    dU = varargin{2};
end
if (beta == 0)
    Summa = 1 ./ max(Sens ./ D);
    output = im + Summa * (im ./ D) .* (rhs);
else
    Summa = 1 ./ max((Sens + dU) ./ (D + dU));
    output = im + Summa .* (im ./ (D + dU)) .* (rhs - dU);
end
end