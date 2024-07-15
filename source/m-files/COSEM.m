function [output, C_co] = COSEM(im, rhs, C_co, D, h, COSEM_TYPE, osa_iter)
if (COSEM_TYPE == 1)
    C_co(:, osa_iter) = rhs .* (im.^(1/h));
	output = (sum(C_co, 2) ./ D).^h;
else
    C_co(:, osa_iter) = rhs .* im;
	output = (sum(C_co, 2) ./ D);
end
end