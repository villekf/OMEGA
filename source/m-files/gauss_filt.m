function output = gauss_filt(hsize, sigma)
%GAUSSFILT Gaussian filter
%   Simple Gaussian filter. Prevents the need to have image processing
%   toolbox. Works essentially the same as fspecial with 'gaussian'.
%
% Example:
%   filter = gauss_filt([size1 size2], sigma)
%
% Inputs:
%   hsize = The 2D size of the filter. E.g. [2 2].
%   sigma = Standard deviation
%
% See also fspecial
ind1 = -floor(hsize(1)/2) : floor(hsize(1)/2);
ind2 = -floor(hsize(2)/2) : floor(hsize(2)/2);
[X, Y] = meshgrid(ind1, ind2);
output = exp(-(X.^2 + Y.^2) / (2*sigma*sigma));
output = output / sum(output(:));
end

