function output = ROSEM(im, Sens, rhs, lam, iter)
output = im + lam(iter) .* im ./ Sens .* (rhs - Sens);
end