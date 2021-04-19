function A = observation_matrix_formation(options, varargin)
%% observation_matrix_formation
% Precomputes the system matrix for PET reconstruction to be used in
% equations of form y = Ax.
%
% Input arguments are obtained from custom_prior_prepass, see
% main_PET.m for an example on how to use this function.
%
% The output is the system matrix. The matrix is transposed if
% options.precompute_lor = true, otherwise it is in its original form.
%
% Examples:
%   A = observation_matrix_formation(options, current_subset)
%   A = observation_matrix_formation(options)
% INPUTS:
%   options = Contains the machine and (optionally) sinogram specific
%   information (e.g. detectors per ring, image (matrix) size, FOV size).
%   current_subset = The current subset number (optional), if omitted then
%   the entire system matrix is computed.
% OUTPUTS:
%   A = The system matrix with all the selected corrections applied.
%
% See also custom_prior_prepass

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ville-Veikko Wettenhovi
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


if nargin > 1 && options.subsets > 1
    osa_iter = varargin{1};
else
    osa_iter = 1;
end
if ~isfield(options,'TOF_bins')
    options.TOF_bins = 1;
end
if ~isfield(options,'listmode')
    options.listmode = false;
end

SinD = 0;

if ~options.use_raw_data
    if isempty(options.pseudot)
        options.pseudot = int32(0);
    end
end


iij = double(0:options.Nx);
jji = double(0:options.Ny);
kkj = double(0:options.Nz);

if exist('feature','builtin') == 5
    nCores = uint32(feature('numcores'));
else
    nCores = uint32(1);
end

if options.subsets > 1
    if options.normalization_correction
        norm_input = options.normalization(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1));
    else
        norm_input = 0;
    end
    if options.scatter_correction && ~options.subtract_scatter
        scatter_input = options.ScatterC(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1));
    else
        scatter_input = 0;
    end
else
    if options.normalization_correction
        norm_input = options.normalization;
    else
        norm_input = 0;
    end
    if options.scatter_correction && ~options.subtract_scatter
        scatter_input = options.ScatterC;
    else
        scatter_input = 0;
    end
end

if options.precompute_lor == false
    if ~isfield(options, 'index') || length(options.index) == 1 || (iscell(options.index) && length(options.index{osa_iter}) == 1)
        if options.use_raw_data
            options.index = {uint32(1 : options.detectors ^2/2 + options.detectors/2)'};
        else
            options.index = {uint32(1 : options.Nang * options.Ndist * options.NSinos)'};
        end
    end
    if iscell(options.index)
        options.index = cell2mat(options.index);
    end
end
koko = options.pituus(osa_iter + 1) - options.pituus(osa_iter);
[A,~] = computeImplementation1(options,options.use_raw_data,options.randoms_correction, options.pituus,osa_iter, options.normalization_correction,...
    options.Nx, options.Ny, options.Nz, options.dx, options.dy, options.dz, options.bx, options.by, options.bz, options.x, options.y, options.z_det, ...
    options.xx, options.yy, options.size_x, options.NSinos, options.NSlices, options.zmax, options.attenuation_correction, options.pseudot, options.det_per_ring, ...
    false, 0, 0, uint32(0), nCores, options.ind_size, options.block1, options.blocks, options.index, iij, jji, kkj, options.LL, options.N, options.summa, options.lor_a, options.xy_index, options.z_index, ...
    options.x_center, options.y_center, options.z_center, options.bmin, options.bmax, options.Vmax, options.V, options.lor_orth, options.gaussK,options.is_transposed, scatter_input, norm_input, SinD, koko);
end