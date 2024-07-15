function output = DRAMA(im, Sens, rhs, lam, iter, sub_iter, subsets)
output = im + lam(iter * subsets + sub_iter) .* im ./ Sens .* rhs;
end