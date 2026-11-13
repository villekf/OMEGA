function lm = isListModeCandidate(options)
%% Detect list-mode data with custom detector coordinates
% Returns true if the options-struct contains custom detector coordinates
% (options.x, options.y, options.z or options.z_det) that are sized to match
% the measurement data (options.SinM), i.e. list-mode data without a scanner
% geometry.
%
% INPUT:
%   options = The options-struct (with the default values set, see
%   setMissingValues)
%
% OUTPUT:
%   lm = True if the data look like list-mode data with custom coordinates
%
% See also OMEGA_error_check, checkSinogramOptions,
% printReconstructionSettings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2026 Ville-Veikko Wettenhovi
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
% List-mode data with custom detector coordinates (options.x/y/z(_det) sized
% to match options.SinM) or index-based reconstruction (options.useIndexBasedReconstruction)
% collapses Nang/Ndist to 1 with no matching scanner geometry (det_w_pseudo can be 0), so it
% must never be validated against the sinogram angle-count limit. This mirrors
% projectorClass.m's own list-mode/index-based detection (~line 561-582), replicated here
% (rather than relying on options.listmode) because the callers of this function run,
% both directly from main-files and from projectorClass.m (~line 421), BEFORE that detection runs -- at
% that point options.listmode still only reflects whatever the caller explicitly set (default
% false), not the auto-detected custom-coordinates case.
lm = (isfield(options,'x') || isfield(options,'y') || (isfield(options,'z') || isfield(options,'z_det'))) && ...
    (((numel(options.x) / 2 == numel(options.SinM) || numel(options.x) / 6 == numel(options.SinM))) || ...
    (numel(options.SinM) == 0 && numel(options.x) >= 6) || ...
    (iscell(options.x) && (numel(options.x{1}) / 2 == numel(options.SinM{1}) || numel(options.x{1}) / 6 == numel(options.SinM{1}))));
end
