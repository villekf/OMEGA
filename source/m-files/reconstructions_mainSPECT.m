function [pz,varargout] = reconstructions_mainSPECT(options)
%RECONSTRUCTIONS_MAINSPECT Utility function for SPECT reconstruction
%   This function simply adds PET specific variables that are needed for
%   error-free functionality and also converts some input variables to the
%   correct format (e.g. the projection angles are converted to radians if
%   they are input as degrees).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2023-2024 Ville-Veikko Wettenhovi, Niilo Saarlemo
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

options.SPECT = true;
options.errorChecking = true;
if nargout > 1
    if nargout == 3
        [pz,varargout{1},varargout{2}] = reconstructions_main(options);
    elseif nargout == 4
        [pz,varargout{1},varargout{2},varargout{3}] = reconstructions_main(options);
    elseif nargout == 5
        [pz,varargout{1},varargout{2},varargout{3},varargout{4}] = reconstructions_main(options);
    else
        [pz,varargout{1}] = reconstructions_main(options);
    end
else
    pz = reconstructions_main(options);
end
end

