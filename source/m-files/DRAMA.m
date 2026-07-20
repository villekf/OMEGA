function output = DRAMA(im, Sens, rhs, lam, iter, sub_iter, subsets)
output = im + lam(iter, sub_iter) .* im ./ Sens .* rhs;
end