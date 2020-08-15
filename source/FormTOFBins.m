function [bins, discard] = FormTOFBins(options, timeDiff)
%FORMTOFBINS Forms the TOF bins from the input time differences
%   Adds Gaussian noise with the specified FWHM.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Jonna Kangasniemi, Ville-Veikko Wettenhovi
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var_t = (options.TOF_FWHM / (2 * sqrt(2 * log(2))))^2;%options.TOF_FWHM is the wanted time accuracy
TOF_data = timeDiff + sqrt(var_t) * randn(size(timeDiff)); %TOF_data with added error
discard = abs(TOF_data) > (options.TOF_width / 2 * options.TOF_bins);
TOF_data(discard) = [];
bins = uint16(floor((abs(TOF_data) + options.TOF_width) / options.TOF_width));
tInd = TOF_data < 0;
bins(bins > 1 & ~tInd) = bins(bins > 1 & ~tInd) * 2 - 1;
bins(tInd) = bins(tInd) * 2;
end