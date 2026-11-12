function options = checkSinogramOptions(options)
%% Sinogram related error checking
% Checks the sinogram parameters (span, ring difference, number of angles,
% number of sinograms, ndist_side, sampling, arc correction and gap filling)
% of PET data. Some of the values are corrected if they are not allowed.
% This function is called by OMEGA_error_check and can also be called
% directly by the functions that form the sinograms.
%
% INPUT:
%   options = The options-struct (with the default values set, see
%   setMissingValues)
%
% OUTPUT:
%   options = The (possibly corrected) options-struct
%
% See also OMEGA_error_check, form_sinograms
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

if ~options.CT && ~options.SPECT && options.span > options.ring_difference && options.NSinos > 1 && ~options.use_raw_data
    error(['Span value cannot be larger than ring difference (' num2str(options.ring_difference) ')!'])
end
if ~options.CT && ~options.SPECT && (mod(options.span,2) == 0 || options.span <= 0) && ~options.use_raw_data
    error('Span value has to be odd and positive!')
end
if isfield(options, 'rings') && ~options.CT && ~options.SPECT && options.ring_difference >= options.rings(end) && ~options.use_raw_data
    warning(['Ring difference can be at most ' num2str(options.rings-1) '. Setting it to the maximum possible.'])
    options.ring_difference = options.rings - 1;
end
if ~options.CT && ~options.SPECT && options.ring_difference < 0 && ~options.use_raw_data
    error('Ring difference has to be at least 0!')
end
% List-mode data with custom detector coordinates (options.x/y/z(_det) sized
% to match options.SinM) or index-based reconstruction (options.useIndexBasedReconstruction)
% collapses Nang/Ndist to 1 with no matching scanner geometry (det_w_pseudo can be 0), so it
% must never be validated against the sinogram angle-count limit below. This mirrors
% projectorClass.m's own list-mode/index-based detection (~line 561-582), replicated here
% (rather than relying on options.listmode) because THIS function is called, both directly
% from main-files and from projectorClass.m (~line 421), BEFORE that detection runs -- at
% that point options.listmode still only reflects whatever the caller explicitly set (default
% false), not the auto-detected custom-coordinates case.
isListModeCandidate = (isfield(options,'x') || isfield(options,'y') || (isfield(options,'z') || isfield(options,'z_det'))) && ...
    (((numel(options.x) / 2 == numel(options.SinM) || numel(options.x) / 6 == numel(options.SinM))) || ...
    (numel(options.SinM) == 0 && numel(options.x) >= 6) || ...
    (iscell(options.x) && (numel(options.x{1}) / 2 == numel(options.SinM{1}) || numel(options.x{1}) / 6 == numel(options.SinM{1}))));
if ~options.CT && ~options.SPECT && options.Nang(end) > options.det_w_pseudo(end)/2 && ~options.use_raw_data ...
        && ~isListModeCandidate && ~options.listmode && ~options.useIndexBasedReconstruction
    error(['Number of sinogram angles can be at most the number of detectors per ring divided by two(' num2str(options.det_w_pseudo(1)/2) ')!'])
end
if ~options.CT && ~options.SPECT && options.TotSinos < options.NSinos && ~options.use_raw_data
    error(['The number of sinograms used (' num2str(options.NSinos) ') is larger than the total number of sinograms (' num2str(options.TotSinos) ')!'])
end
if ~options.CT && ~options.SPECT && (options.ndist_side > 1 && mod(options.Ndist,2) == 0 || options.ndist_side < -1 && mod(options.Ndist,2) == 0) && ~options.use_raw_data
    error('ndist_side can be either 1 or -1!')
end
if ~options.CT && ~options.SPECT && options.ndist_side == 0 && mod(options.Ndist,2) == 0 && ~options.use_raw_data
    error('ndist_side cannot be 0 when Ndist is even!')
end
if ~options.CT && ~options.SPECT && ((mod(options.sampling, 2) > 0 && options.sampling ~= 1) || options.sampling < 0) && ~options.use_raw_data
    error('Sampling rate has to be divisible by two and positive or one!')
end
if ~options.CT && ~options.SPECT && ((mod(options.sampling_raw, 2) > 0 && options.sampling_raw ~= 1) || options.sampling_raw < 0) && options.use_raw_data
    error('Sampling rate has to be divisible by two and positive or one!')
end
if ~options.CT && ~options.SPECT && ((options.sampling > 1 && ~options.use_raw_data) || (options.sampling_raw > 1 && options.use_raw_data)) && options.implementation == 1
    warning('Increased sampling rate is not supported with implementation 1. Using sample rate of 1.')
    options.sampling = 1;
end
if ~options.CT && ~options.SPECT && options.arc_correction && options.use_raw_data
    warning('Arc correction is not supported for raw data. Disabling arc correction.')
    options.arc_correction = false;
end
if ~options.CT && ~options.SPECT && options.arc_correction && options.implementation == 1
    warning('Arc correction is not supported with implementation 1. Disabling arc correction.')
    options.arc_correction = false;
end
if ~options.CT && ~options.SPECT && options.det_per_ring(1) == options.det_w_pseudo(1) && options.fill_sinogram_gaps
    error('Gap filling is only supported with pseudo detectors!')
end
