function [f, g, h] = computeGradient(im, options, f, g, h, varargin)
if (options.verbose >= 3)
    disp("Starting forward difference gradient")
end
f(1 : end - 1, :, :) = diff(im);
f(end, :, :) = -1 * im(end, :, :);
g(:, 1 : end - 1, :) = diff(im, 1, 1);
g(:, end, :) = -1 * im(:, end, :);
h(:, :, 1 : end - 1) = diff(im, 1, 2);
h(:, :, end) = -1 * im(:, :, end);