function [output] = EM(im, sens, rhs)
% min(sens)
% max(sens)
output = im ./ sens .* rhs;
end