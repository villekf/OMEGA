function visualizeRawData(data, detectors)
%VISUALIZERAWDATA Raw data visualization
%   Visualizes the raw list-mode data created by load_data. Input the raw
%   data and the total number of detectors on the machine.
%
% Example:
%   visualizeRawData(coincidences{1}, options.detectors)

[K, ~, V] = find(data);
L = find(tril(true(detectors,detectors), 0));
L = L(K);

[I,J] = ind2sub([detectors detectors], L);
clear L
true_coincidences = sparse(I,J,V,detectors,detectors);
clear I J V
    
imagesc(true_coincidences)
end

