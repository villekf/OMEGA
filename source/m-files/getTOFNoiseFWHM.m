function TOFNoiseFWHM = getTOFNoiseFWHM(options)
%GETTOFNOISEFWHM FWHM (s) of the temporal noise added to the TOF time differences
%   Uses options.TOF_noise_FWHM if it exists and is non-negative (zero means
%   that no noise is added). Otherwise (missing or negative, the default is
%   -1) the noise has the same FWHM as the TOF data, options.TOF_FWHM.

if isfield(options, 'TOF_noise_FWHM') && ~isempty(options.TOF_noise_FWHM) && options.TOF_noise_FWHM >= 0
    TOFNoiseFWHM = options.TOF_noise_FWHM;
elseif isfield(options, 'TOF_FWHM') && ~isempty(options.TOF_FWHM)
    TOFNoiseFWHM = options.TOF_FWHM;
else
    TOFNoiseFWHM = 0;
end
end
