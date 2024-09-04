function filt = rampFilt(N,varargin)
%RAMPFILT Outputs the (optionally) windowed ramp filter
%   Hamming, Hann, Blackman, Nuttal, Gaussian, Shepp-Logan, cosine and
%   Parzen windowing methods available (second input). Third input is the
%   optional cut-off frequency value. Fourth is the optional STD for the
%   Gaussian filter. Fift input is the optional boolean for 2D filtering,
%   i.e. if true 2D filter is output instead.
if nargin > 2 && ~isempty(varargin) && ~isempty(varargin{2})
    c = varargin{2};
else
    c = 1;
end
if nargin > 3 && ~isempty(varargin) && ~isempty(varargin{3})
    sigma = varargin{3};
else
    sigma = 0.25;
end
if nargin > 1 && ~isempty(varargin) && ~isempty(varargin{1})
    ww = varargin{1};
else
    ww = 'rect';
end
if nargin > 4 && ~isempty(varargin) && ~isempty(varargin{4})
    use2D = varargin{4};
else
    use2D = false;
end

filt = linspace(0, N, N / 2 + 1)' ./ N;
w = 2 * pi * (0:size(filt,1) - 1)' / N;
switch ww
    case 'hamming'
        w = (.54 + .46 * cos(w / c));
        filt(2:end) = filt(2:end) .* w(2:end);
    case 'hann'
        w = (.5 + .5 * cos(w / c));
        filt(2:end) = filt(2:end) .* w(2:end);
    case 'blackman'
        w = (.42659 - .49656 * cos(w / c) + .076849 * cos(2*w / c));
        filt(2:end) = filt(2:end) .* w(2:end);
    case 'nuttal'
        w = (.355768 - .487396 * cos(w / c) + .144232 * cos(2*w / c) - .012604 * cos(3*w / c));
        filt(2:end) = filt(2:end) .* w(2:end);
    case 'gaussian'
        w = exp((-1/2) * (((size(filt,1)-1:-1:0)' - N/2)/(sigma*(N/2))) .^ 2);
        filt(2:end) = filt(2:end) .* w(2:end);
    case 'shepp-logan'
        w = sin(w / (2 * c)) ./ (w / (2 * c));
        filt(2:end) = filt(2:end) .* w(2:end);
    case 'cosine'
        filt(2:end) = filt(2:end) .* cos(w(2:end) / (2 * c));
    case 'parzen' % de la Vallée Poussin
        L = N + 1;
        w = (0:size(filt,1) - 1)';
        w(abs(w) <= L/4) = 1 - 6 * (w(abs(w) <= L/4) / (L/2)).^2 .* (1 - abs(w(abs(w) <= L/4)) / (L/2));
        w(abs(w) > L/4 & abs(w) <= L/2) = 2 * (1 - abs(w(abs(w) > L/4 & abs(w) <= L/2)) / (L/2)).^3;
        filt(2:end) = filt(2:end) .* w(2:end);
end
filt(w>pi*c) = 0;
if mod(N,2) > 0
    filt = [filt; filt(end-1:-1:1)];
else
    filt = [filt; filt(end-1:-1:2)];
end
if use2D
    filt = repmat(filt, 1, N);
    filt = filt .* filt';
    filt = filt ./ max(filt(:));
end
end
