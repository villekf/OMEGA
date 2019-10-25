function A = observation_matrix_formation(options, varargin)
%% observation_matrix_formation
% Precomputes the system matrix for PET reconstruction to be used in
% equations of form y = Ax.
%
% Input arguments are obtained from custom_prior_prepass, see
% main_nongate.m for an example on how to use this function. 
%
% The output is the system matrix. The matrix is transposed if
% options.precompute_lor = true, otherwise it is in its original form.
%
% See also custom_prior_prepass

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi
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


if nargin > 2 && options.subsets > 1
    osa_iter = varargin{1};
else
    osa_iter = 1;
end

SinD = 0;

if ~options.use_raw_data
    if isempty(options.pseudot)
        options.pseudot = int32(0);
    end
end

if options.precompute_lor == false
    if ~isfield(options, 'index') || length(options.index) == 1 || length(options.index{osa_iter}) == 1
        if options.use_raw_data
            options.index = {uint32(1 : options.detectors ^2/2 + options.detectors/2)'};
        else
            options.index = {uint32(1 : options.Nang * options.Ndist * options.NSinos)'};
        end
    end
    if options.use_raw_data == false
        if options.projector_type == 1 || options.projector_type == 0
            if exist('projector_mex','file') == 3
                [ lor, indices, alkiot] = projector_mex( options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, options.bx, options.bz, options.z_det, options.x, options.y, options.dy, options.yy, options.xx , options.NSinos, options.NSlices, options.size_x, ...
                    options.zmax, options.vaimennus, options.normalization, SinD, options.pituus(osa_iter), options.attenuation_correction, options.normalization_correction, ...
                    options.randoms_correction, uint16(0), uint32(0), uint32(0), options.NSinos, uint16(0), options.pseudot, options.det_per_ring, options.verbose, ...
                    options.use_raw_data, uint32(2), options.ind_size, options.block1, options.blocks, options.index{osa_iter}, uint32(options.projector_type));
            else
                % The below lines allow for pure MATLAB
                % implemention, i.e. no MEX-files will be
                % used. Currently the below function
                % uses parfor-loops (requires parallel
                % computing toolbox).
                % NOTE: The below method is not
                % recommended since it is much slower
                % method.
                [ lor, indices, alkiot, discard] = improved_siddon_atten( options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, options.bx, options.bz, options.z_det, options.x, options.y, options.yy, options.xx, ...
                    options.NSinos, options.NSlices, options.vaimennus, options.index{osa_iter}, options.pituus(osa_iter), options.attenuation_correction);
                alkiot = cell2mat(alkiot);
                indices = indices(discard);
                indices = cell2mat(indices) - 1;
            end
        elseif options.projector_type == 2
            [ lor, indices, alkiot] = projector_mex( options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, options.bx, options.bz, options.z_det, options.x, options.y, options.dy, options.yy, options.xx , options.NSinos, options.NSlices, options.size_x, ...
                options.zmax, options.vaimennus, options.normalization, SinD, options.pituus(osa_iter), options.attenuation_correction, options.normalization_correction, ...
                options.randoms_correction, uint16(0), uint32(0), uint32(0), options.NSinos, uint16(0), options.pseudot, options.det_per_ring, options.verbose, ...
                options.use_raw_data, uint32(2), options.ind_size, options.block1, options.blocks, options.index{osa_iter}, uint32(options.projector_type), ...
                options.tube_width_xy, options.x_center, options.y_center, options.z_center, options.tube_width_z);
        else
            error('Unsupported projector type')
        end
    else
        if options.subsets > 1
            L = options.LL(options.pituus(osa_iter) * 2 + 1 : options.pituus(osa_iter + 1) * 2);
        else
            L = options.LL;
        end
        if options.projector_type == 1
            [ lor, indices, alkiot] = projector_mex( options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, options.bx, options.bz, options.z_det, options.x, options.y, dy, options.yy, options.xx , options.NSinos, options.NSlices, options.size_x, ...
                options.zmax, options.vaimennus, options.normalization, SinD, uint32(0), options.attenuation_correction, options.normalization_correction, ...
                options.randoms_correction, uint16(0), uint32(0), uint32(0), options.NSinos, L, options.pseudot, options.det_per_ring, options.verbose, ...
                options.use_raw_data, uint32(2), options.ind_size, options.block1, options.blocks, uint32(0), uint32(options.projector_type));
        elseif options.projector_type == 2
            [ lor, indices, alkiot] = projector_mex( options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, options.bx, options.bz, options.z_det, options.x, options.y, dy, options.yy, options.xx , options.NSinos, options.NSlices, options.size_x, ...
                options.zmax, options.vaimennus, options.normalization, SinD, uint32(0), options.attenuation_correction, options.normalization_correction, ...
                options.randoms_correction, uint16(0), uint32(0), uint32(0), options.NSinos, L, options.pseudot, options.det_per_ring, options.verbose, ...
                options.use_raw_data, uint32(2), options.ind_size, options.block1, options.blocks, uint32(0), uint32(options.projector_type), options.tube_width_xy, ...
                options.x_center, options.y_center, options.z_center, options.tube_width_z);
        else
            error('Unsupported projector type')
        end
    end
%     lor = reshape(lor,[],2);
    A_length = length(lor);
    if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
        lor = repeat_elem(int32(1:length(lor))',int32(lor));
    elseif exist('OCTAVE_VERSION','builtin') == 5
        lor = repelem(int32(1:length(lor)),int32(lor));
    else
        lor = repelem(int32(1:length(lor)),int32(lor))';
    end
    
    indices=indices + 1;
    if options.verbose
        tStart = tic;
    end
    if options.use_fsparse && exist('fsparse','file') == 3
        A = fsparse(lor,int32(indices),double(alkiot),[A_length double(options.N) length(alkiot)]);
    elseif options.use_fsparse && exist('fsparse','file') == 0
        warning('options.fsparse set to true, but no FSparse mex-file found. Using regular sparse')
        A = sparse(double(lor),double(indices),double(alkiot), A_length, double(options.N));
    else
        A = sparse(double(lor),double(indices),double(alkiot), A_length, double(options.N));
    end
    clear indices alkiot lor
    if options.verbose
        tElapsed = toc(tStart);
        disp(['Sparse matrix formation took ' num2str(tElapsed) ' seconds'])
    end
else
    if options.use_raw_data && options.subsets > 1
        L_input = options.LL(options.pituus(osa_iter) * 2 + 1 : options.pituus(osa_iter + 1) * 2);
        xy_index_input = uint32(0);
        z_index_input = uint16(0);
    elseif options.subsets > 1
        L_input = uint16(0);
        xy_index_input = options.xy_index(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1));
        z_index_input = options.z_index(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1));
    end
    if options.projector_type == 2 && options.subsets > 1
        lor2 = [0; cumsum(uint64(options.lor_orth(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1))))];
    elseif options.projector_type == 2 && options.subsets == 1
        lor2 = [0; uint64(options.lor_orth)];
    elseif options.projector_type == 1 && options.subsets == 1
        lor2 = [0; uint64(options.lor_a)];
    else
        lor2 = [0; cumsum(uint64(options.lor_a(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1))))];
    end
    if options.subsets > 1
        [A, ~] = projector_mex( options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, options.bx, options.bz, options.z_det, options.x, options.y, options.dy, options.yy, options.xx , options.NSinos, options.NSlices, options.size_x, options.zmax, options.vaimennus, ...
            options.normalization, SinD, options.pituus(osa_iter + 1) - options.pituus(osa_iter), options.attenuation_correction, options.normalization_correction, options.randoms_correction, ...
            options.lor_a(options.pituus(osa_iter)+1:options.pituus(osa_iter + 1)), xy_index_input, z_index_input, options.NSinos, L_input, options.pseudot, options.det_per_ring, options.verbose, ...
            options.use_raw_data, uint32(0), lor2, options.summa(osa_iter), false, uint32(options.projector_type), options.tube_width_xy, options.x_center, ...
            options.y_center, options.z_center, options.tube_width_z, int32(options.accuracy_factor));
    else
        [A, ~] = projector_mex( options.Ny, options.Nx, options.Nz, options.dx, options.dz, options.by, options.bx, options.bz, options.z_det, options.x, options.y, options.dy, options.yy, options.xx , options.NSinos, options.NSlices, options.size_x, options.zmax, options.vaimennus, ...
            options.normalization, SinD, options.pituus(osa_iter + 1) - options.pituus(osa_iter), options.attenuation_correction, options.normalization_correction, options.randoms_correction, ...
            options.lor_a, options.xy_index, options.z_index, options.NSinos, options.LL, options.pseudot, options.det_per_ring, options.verbose, ...
            options.use_raw_data, uint32(0), lor2, options.summa, false, uint32(options.projector_type), options.tube_width_xy, options.x_center, ...
            options.y_center, options.z_center, options.tube_width_z, int32(options.accuracy_factor));
    end
    clear lor2
end
end