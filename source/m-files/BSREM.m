function output = BSREM(im, rhs, lam, iter)
output = im + lam(iter) .* im .* rhs;
end