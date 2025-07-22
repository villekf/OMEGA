function [HU_a, HU_b] = HU_rescale_coeff(fpathCT)
    files = dir(fullfile(fpathCT, '*.dcm'));
    n = numel(files);
    HU_a = zeros(n, 1);
    HU_b = zeros(n, 1);
    for ii = 1:n
        f = fullfile(fpathCT, files(ii).name);
        info = dicominfo(f);
        HU_a(ii) = double(info.RescaleSlope);
        HU_b(ii) = double(info.RescaleIntercept);
    end
    HU_a = reshape(HU_a, 1, 1, n);
    HU_b = reshape(HU_b, 1, 1, n);
end