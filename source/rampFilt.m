function filt = rampFilt(N,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if ~isempty(varargin) && ~isempty(varargin{1})
    c = varargin{1};
else
    c = 1;
end
nn = -(N/2):N/2-1;
filt = zeros(N,1);
filt(N/2 + 1) = 1/4;
filt(2:2:end) = -1 ./ (pi * nn(2:2:end)).^2;
filt = abs(fft(filt))*2;
filt = filt(1:N/2+1);
w = 2*pi*(0:size(filt,1)-1)'/N;
filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end) / c));
filt(w>pi*c) = 0;
filt = [filt; filt(end-1:-1:2)];
end

nn = -(N/2):N/2-1;
mm = -(N/2):N/2-1;
[NN, MM] = meshgrid(nn, mm);
% filt = zeros(N,N);
% filt = -1 ./ (2*pi.^2 * MM.^2 .* sinc(-pi * NN));
filt = -1 ./ (2*pi.^2 * (MM.^2 .* NN.^2));
filt(N/2 + 1,:) = 0;
filt(:,N/2 + 1) = 0;
filt(N/2 + 1,N/2 + 1) = 1/4;