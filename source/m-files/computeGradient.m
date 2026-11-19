function [f, g, h] = computeGradient(im, options, f, g, h, varargin)
%COMPUTEGRADIENT Computes the first order differences of a 3D image along
%the three dimensions (computeGradient in functions.hpp)
%   The optional sixth input is the difference type: 0 = forward
%   differences (default), 1 = backward differences, 2 = central
%   differences (one-sided differences at the boundaries)
if nargin >= 6 && ~isempty(varargin{1})
    type = varargin{1};
else
    type = 0;
end
if type == 0
    % 1st order forward differences
    if (options.verbose >= 3)
        disp("Starting forward difference gradient")
    end
    f(1 : end - 1, :, :) = diff(im);
    f(end, :, :) = -1 * im(end, :, :);
    g(:, 1 : end - 1, :) = diff(im, 1, 2);
    g(:, end, :) = -1 * im(:, end, :);
    h(:, :, 1 : end - 1) = diff(im, 1, 3);
    h(:, :, end) = -1 * im(:, :, end);
elseif type == 1
    % 1st order backward differences
    if (options.verbose >= 3)
        disp("Starting backward difference gradient")
    end
    f(2 : end, :, :) = -diff(im);
    f(1, :, :) = -1 * im(1, :, :);
    g(:, 2 : end, :) = -diff(im, 1, 2);
    g(:, 1, :) = -1 * im(:, 1, :);
    h(:, :, 2 : end) = -diff(im, 1, 3);
    h(:, :, 1) = -1 * im(:, :, 1);
else
    % 1st order central differences, one-sided differences at the boundaries
    % (a dimension of length 1 has zero difference)
    if (options.verbose >= 3)
        disp("Starting central difference gradient")
    end
    f = (circshift(im, -1, 1) - circshift(im, 1, 1)) * .5;
    if size(im, 1) > 1
        f(1, :, :) = im(2, :, :) - im(1, :, :);
        f(end, :, :) = im(end, :, :) - im(end - 1, :, :);
    end
    g = (circshift(im, -1, 2) - circshift(im, 1, 2)) * .5;
    if size(im, 2) > 1
        g(:, 1, :) = im(:, 2, :) - im(:, 1, :);
        g(:, end, :) = im(:, end, :) - im(:, end - 1, :);
    end
    h = (circshift(im, -1, 3) - circshift(im, 1, 3)) * .5;
    if size(im, 3) > 1
        h(:, :, 1) = im(:, :, 2) - im(:, :, 1);
        h(:, :, end) = im(:, :, end) - im(:, :, end - 1);
    end
end
if (options.verbose >= 3)
    disp("Gradient computed")
end
