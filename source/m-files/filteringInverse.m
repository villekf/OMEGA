function [input] = filteringInverse(filter, input, Nf)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

temp = fft(input, Nf);
temp = bsxfun(@rdivide, temp, filter);
temp = ifft(temp);
input = real(temp(1 : size(input,1), :, :));
input = input(:);
end