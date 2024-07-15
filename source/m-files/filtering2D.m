function [input] = filtering2D(filter, input, dimmi)
temp = fft2(input, dimmi, dimmi);
temp = bsxfun(@times, temp, filter);
temp = ifft2(temp);
input = real(temp(1 : size(input,1), 1 : size(input,2), :));
input = input(:);