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


% var = recNames(4);
% ll = 0;
% kk = 1;
% while ll == 0 && kk <= numel(var)
%     ll = ll + options.(var{kk});
%     kk = kk +1;
% end

if ~isfield(options, 'CORtoDetectorSurface')
    options.CORtoDetectorSurface = 0;
end
if ~isfield(options, 'swivelAngles')
    options.swivelAngles = options.angles+180;
end
if ~isfield(options, 'rayShiftsDetector')
    %options.rayShiftsDetector = single(options.colR*(2*rand(2*options.nRays, 1)-1)/options.crXY);
    options.rayShiftsDetector = single(zeros(2*options.nRays, 1));
    options.rayShiftsDetector(1:2) = 0;
end
if ~isfield(options, 'rayShiftsSource')
    options.rayShiftsSource = single(0.2*randn(2*options.nRays, 1));
    options.rayShiftsSource(1:2) = 0;
    %options.rayShiftsSource = single(options.colR*(2*rand(2*options.nRays, 1)-1)/options.crXY);

end
if ~isfield(options, 'homeAngles')
    options.homeAngles = 0 * options.angles;
end
if ~isfield(options, 'coneOfResponseStdCoeff')
    options.coneOfResponseStdCoeff = 0.1;
end
if options.flipImageZ % Flip image by flipping sinograms' z-axis in image space
    options.SinM = flip(options.SinM, 2);
end

options.n_rays_transaxial = options.nRays;
options.NSinos = options.nProjections;
options.TotSinos = options.nProjections;
options.span = 3;
options.segment_table = [];
options.Ndist = options.nRowsD;
options.Nang = options.nColsD;
options.use_raw_data = false;
options.randoms_correction = false;
options.scatter_correction = false;
%options.attenuation_correction = false;
options.normalization_correction = false;
if ~isfield(options, 'partitions')
    options.partitions = 1;
end
options.use_psf = false;
options.reconstruct_trues = false;
options.reconstruct_scatter = false;
options.machine_name = 'CT';
options.corrections_during_reconstruction = false;
options.precompute_obs_matrix = false;
options.precompute_all = false;
options.compute_normalization = false;
options.diameter = 0;
options.ring_difference = 0;
options.ndist_side = 0;
options.sampling_raw = 1;
options.pseudot = [];
if ~isfield(options, 'nBed')
    options.nBed = 1;
end
options.rings = options.nColsD * options.nBed;
options.det_per_ring = options.nRowsD * options.nProjections;
options.precompute_lor = false;
options.sampling = 1;
options.CT = false;
options.PET = false;
options.SPECT = true;
options.start = 0;
options.end = 0;
options.arc_correction = false;
options.tot_time = 0;
options.fill_sinogram_gaps = false;
options.det_w_pseudo = options.det_per_ring;
options.blocks_per_ring = 1;
options.linear_multip = 0;
options.cryst_per_block = options.nColsD * options.nRowsD;
options.dPitch = options.cr_p;
%options.tube_width_xy = 0;
%options.tube_width_z = 0;
if options.projector_type == 2 && isfield(options, 'gFilter')
    options = rmfield(options, 'gFilter');
end
if ~isfield(options, 'cr_pz')
    % options.cr_pz = options.cr_p;
    options.cr_pz = options.crXY;
end
if ~isfield(options, 'n_rays_transaxial')
    options.n_rays_transaxial = 1;
end
if ~isfield(options, 'n_rays_axial')
    options.n_rays_axial = 1;
end
if ~isfield(options, 'horizontalOffset')
    options.horizontalOffset = 0;
end
if ~isfield(options, 'verticalOffset')
    options.verticalOffset = 0;
end
if ~isfield(options, 'bedOffset')
    options.bedOffset = 0;
end
if ~isfield(options, 'flip_image')
    options.flip_image = false;
end
if ~isfield(options, 'offangle')
    options.offangle = 0;
else
    if isfield(options, 'vaimennus')
        options.vaimennus = imrotate3(options.vaimennus, options.offangle, [0 0 1], 'crop');
    end
end
if ~isfield(options, 'uCenter')
    options.uCenter = [];
end
if ~isfield(options, 'vCenter')
    options.vCenter = [];
end
options = OMEGA_error_check(options);
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

