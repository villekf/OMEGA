function visualizeRawData(data, detectors, varargin)
%VISUALIZERAWDATA Raw data visualization
%   Visualizes the raw list-mode data created by load_data. Input the raw
%   data and the total number of detectors on the machine. Optionally input
%   a sequence of rings from where to visualize along with the total number
%   of rings.
%
% Examples:
%   visualizeRawData(coincidences{1}, options.detectors)
%   visualizeRawData(coincidences{1}, options.detectors, 3:8, 80)

[K, ~, V] = find(data);
L = find(tril(true(detectors,detectors), 0));
L = L(K);

[I,J] = ind2sub([detectors detectors], L);
clear L
true_coincidences = sparse(I,J,V,detectors,detectors);
clear I J V

figure
if isempty(varargin)
    imagesc(true_coincidences)
else
    det_per_ring = detectors / varargin{2};
    alku = varargin{1}(1) - 1;
    if length(varargin{1}) == 1
        loppu = alku + 1;
    else
        loppu = varargin{1}(end) - 1;
    end
    imagesc(true_coincidences(1 + alku * det_per_ring : loppu * det_per_ring,1 + alku * det_per_ring : loppu * det_per_ring ))
end
end