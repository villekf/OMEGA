function options = OMEGA_error_check(options)
%% Error checking file
% This function is used to check that all the input values are allowed
%

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

MAP = (options.OSL_MLEM || options.OSL_OSEM || options.BSREM || options.MBSREM || options.ROSEM_MAP || options.RBI_MAP || any(options.COSEM_MAP));
MAPOS = (options.OSL_OSEM || options.BSREM || options.MBSREM || options.ROSEM_MAP || options.RBI_MAP || any(options.COSEM_MAP));
PRIOR = (options.MRP || options.quad || options.L || options.FMH || options.weighted_mean || options.TV || options.AD || options.APLS || options.TGV || options.NLM);
PRIOR_summa = sum([options.MRP, options.quad, options.L, options.FMH, options.weighted_mean, options.TV, options.AD, options.APLS, options.TGV, options.NLM]);
OS = (options.osem || options.ramla || options.mramla || options.cosem || options.rosem || options.rbi || options.ecosem || options.acosem || options.drama || ...
    options.OSL_OSEM || options.BSREM || options.MBSREM || options.ROSEM_MAP || options.RBI_MAP || any(options.COSEM_MAP));
MLOS = (options.mlem || options.osem);
NMLOS = (options.ramla || options.mramla || options.cosem || options.rosem || options.rbi || options.ecosem || options.acosem || options.drama || ...
    options.OSL_OSEM || options.BSREM || options.MBSREM || options.ROSEM_MAP || options.RBI_MAP || any(options.COSEM_MAP) || options.OSL_MLEM);
OS_I4 = (options.mramla || options.cosem || options.rbi || options.ecosem || options.acosem || options.MBSREM || options.RBI_MAP || any(options.COSEM_MAP));
OS_I4_summa = sum([options.osem, options.ramla, options.rosem, options.drama, options.OSL_OSEM, options.BSREM, options.ROSEM_MAP]);
if options.FOVa_x >= options.diameter || options.FOVa_y >= options.diameter
    error(['Transaxial FOV is larger than the machine diameter (' num2str(options.diameter) ')'])
end
if (options.axial_fov) < (options.linear_multip * options.cryst_per_block * options.cr_pz - options.cr_pz)
    error('Axial FOV is too small, crystal ring(s) on the boundary have no slices')
end
if (options.axial_fov) > (options.linear_multip * options.cryst_per_block * options.cr_pz + options.axial_fov/options.Nz*2 + options.cr_pz*sum(options.pseudot))
    error('Axial FOV is too large, not all the slices have LORs going through them')
end
if options.use_LMF && options.data_bytes < 10
    error('Too little data bytes in LMF format, minimum allowed is 10 bytes (time + detector indices)')
end
if options.use_LMF && options.data_bytes > 21
    warning(['LMF format uses more bytes than the supported 21 bytes (time + detector indices + source coordinates + event indices + Compton scattering in phantom). '...
    'If these extra bytes are before the bytes that are used, output data will be incorrect.'])
end
if options.use_LMF && options.R_bits + options.C_bits + options.M_bits + options.S_bits + options.L_bits > 16
    error('Number of bits used in LMF is more than 16 bits. OMEGA supports only 16 bit detector indices')
end
if options.span > options.ring_difference
    error(['Span value cannot be larger than ring difference (' num2str(options.ring_difference) ')'])
end
if options.span < 1
    error('Span value has to be at least 1')
end
if mod(options.span,2) == 0
    error('Span value has to be odd')
end
if options.ring_difference >= options.rings
    error(['Ring difference can be at most ' num2str(options.rings-1)])
end
if options.ring_difference < 0
    error('Ring difference has to be at least 0')
end
if options.Nang > options.det_w_pseudo/2
    error(['Number of sinogram angles can be at most the number of detectors per ring divided by two(' num2str(options.det_w_pseudo/2) ')'])
end
if options.TotSinos < options.NSinos
    error(['The numnber of sinograms used (' num2str(options.NSinos) ') is larger than the total number of sinograms (' num2str(options.TotSinos) ')'])
end
if options.ndist_side > 1 && mod(options.Ndist,2) == 0 || options.ndist_side < -1 && mod(options.Ndist,2) == 0
    error('ndist_side can be either 1 or -1')
end
if options.ndist_side == 0 && mod(options.Ndist,2) == 0
    error('ndist_side cannot be 0 when Ndist is even')
end
if options.partitions < 1
    warning('Number of partitions is less than one. Using one partition.')
    options.partitions = 1;
end
if options.start > options.end
    error('Start time is later than end time')
end
if options.start > options.tot_time
    error('Start time is larger than the total time of the measurement')
end
if options.Niter < 1
    error('Number of iterations is less than one')
end
if options.subsets < 2 && OS
    warning('Number of subsets is less than two. Subset has to be at least 2 when using OS-methods. Using 2 subsets.')
    options.subsets = 2;
end
if size(options.x0,1)*size(options.x0,2)*size(options.x0,3) < options.Nx*options.Ny*options.Nz
    error(['Initial value has a matrix size smaller (' num2str(size(options.x0,1)*size(options.x0,2)*size(options.x0,3)) ') than the actual image size ('...
        num2str(options.Nx*options.Ny*options.Nz) ')'])
end
if size(options.x0,1)*size(options.x0,2)*size(options.x0,3) > options.Nx*options.Ny*options.Nz
    error(['Initial value has a matrix size larger (' num2str(size(options.x0,1)*size(options.x0,2)*size(options.x0,3)) ') than the actual image size (' ...
        num2str(options.Nx*options.Ny*options.Nz) ')'])
end
if isempty(options.x0)
    warning('Initial value is an empty array, using the default values (1)')
    options.x0 = ones(options.Nx, options.Ny, options.Nz);
end
if options.attenuation_correction && isempty(options.attenuation_datafile)
    error('Attenuation correction is selected, but no attenuation file has been specified')
end
if options.attenuation_correction && exist(options.attenuation_datafile, 'file') ~= 2
    error('Attenuation file not found on path')
end
if isunix
    if ~strcmp('/',options.fpath(end))
        options.fpath = [options.fpath '/'];
    end
elseif ispc
    if ~strcmp('\',options.fpath(end)) && ~strcmp('/',options.fpath(end))
        options.fpath = [options.fpath '\'];
    end
else
    if ~strcmp('/',options.fpath(end))
        options.fpath = [options.fpath '/'];
    end
end
if options.use_LMF && options.randoms_correction
    warning('Randoms correction is set to true although LMF input is selected. No randoms correction will be performed.')
    options.randoms_correction = false;
end
if options.TV_use_anatomical && options.TV_OSL && exist(options.TV_reference_image,'file') ~= 2 && MAP
    error('Anatomical reference image for TV was not found on path')
end
% if options.scatter_correction && ~options.only_reconstructions && ~options.use_raw_data
%     if iscell(options.ScatterC)
%         if size(options.ScatterC{1},1)*size(options.ScatterC{1},2)*size(options.ScatterC{1},3) ~= options.Nang*options.Ndist*options.NSlices
%             error('Size mismatch between scatter sinogram and the specified sinogram')
%         end
%     else
%         if size(options.ScatterC,1)*size(options.ScatterC,2)*size(options.ScatterC,3) ~= options.Nang*options.Ndist*options.NSlices
%             error('Size mismatch between scatter sinogram and the specified sinogram')
%         end
%     end
% end
if options.APLS && exist(options.APLS_reference_image,'file') ~= 2 && MAP
    error('APLS selected, but the anatomical reference image was not found on path')
end
if options.epps <= 0
    warning('Epsilon value is zero or less than zero; must be a positive value. Using the default value (1e-8).')
    options.epps = 1e-8;
end
if numel(options.epps) > 1
    warning('Epsilon has to be a scalar value! Using the default value (1e-8).')
    options.epps = 1e-8;
end
if options.store_scatter && sum(options.scatter_components) <= 0
    error('Store scatter selected, but no scatter components have been selected')
end
if options.reconstruction_method == 3
    if options.osem && options.mlem
        warning('Both OSEM and MLEM selected with implementation 3, using only MLEM');
        options.osem = false;
    end
    if options.mlem && options.subsets > 1
        options.subsets = 1;
    end
    if options.osem && options.subsets == 1
        warning('OSEM selected with only 1 subset. Switching to MLEM instead')
        options.osem = false;
        options.mlem = true;
    end
end
if ~isfield(options,'use_machine')
    options.use_machine = 0;
end
if options.use_machine == 2 && options.use_raw_data
    warning('Sinogram data cannot be used when raw data is set to true, using list-mode data instead')
    options.use_machine = 1;
end
% if (options.blank || options.transmission_scan) && options.attenuation_correction
%     error('Both transmission and image based attenuation correction selected. Use only one correction method.')
% end
if options.source && options.use_ASCII && (options.source_index1 == 0 || isempty(options.source_index1) || options.source_index2 == 0 || isempty(options.source_index2))
    error('Source image selected with ASCII data, but no source index column numbers are provided.')
end
if options.reconstruct_trues && options.reconstruct_scatter
    warning('Both reconstruct trues and scatter selected, reconstructing only trues.')
    options.reconstruct_scatter = false;
end
% if options.precompute_lor == false && options.reconstruction_method > 1
%     error('precompute_lor must be set to true if any other reconstruction method than 1 is used.')
% end
if options.reconstruction_method == 1 && exist('projector_mex','file') ~= 3 && options.precompute_lor
    error('MEX-file not found. Run install_mex first.')
end
if options.reconstruction_method == 4 && exist('projector_mex','file') ~= 3
    error('MEX-file not found. Run install_mex first.')
end
if options.use_root && exist('GATE_root_matlab','file') ~= 3
    warning('ROOT selected, but no MEX-file for ROOT data load found. Run install_mex to build ROOT MEX-file.')
end
if options.use_LMF && exist('gate_lmf_matlab','file') ~= 3
    error('LMF selected, but no MEX-file for LMF data load found. Run install_mex to build LMF MEX-file.')
end
if options.reconstruction_method == 2 && exist('OpenCL_matrixfree','file') ~= 3
    error('OpenCL reconstruction selected, but OpenCL MEX-files were not installed. Run install_mex to build OpenCL MEX-files.')
end
if options.reconstruction_method == 3 && exist('OpenCL_matrixfree_multi_gpu','file') ~= 3
    error('OpenCL reconstruction selected, but OpenCL MEX-files were not installed. Run install_mex to build OpenCL MEX-files.')
end
if options.reconstruction_method == 5 && exist('improved_Siddon_openCL','file') ~= 3
    error('OpenCL reconstruction selected, but OpenCL MEX-files were not installed. Run install_mex to build OpenCL MEX-files.')
end
if options.reconstruction_method == 3 && NMLOS
    warning(['Implementation ' num2str(options.reconstruction_method) ' selected with reconstruction algorithms other than MLEM or OSEM. '...
        'Only MLEM or OSEM reconstruction can be performed'])
end
if options.reconstruction_method == 3 && ~MLOS
    error(['Implementation ' num2str(options.reconstruction_method) ' selected, but neither MLEM nor OSEM algorithm has been selected.'])
end
if options.reconstruction_method == 4 && OS_I4
    error(['Implementation ' num2str(options.reconstruction_method) ' selected with unsupported algorithm. Only MLEM, OSEM, ROSEM,'...
       ' RAMLA (and their MAP-methods) and DRAMA are supported!'])
end
if options.reconstruction_method == 4 && PRIOR_summa > 1
    error(['Implementation ' num2str(options.reconstruction_method) ' supports only one prior at a time.'])
end
if options.reconstruction_method == 4 && (PRIOR_summa == 1 && ((options.osem && options.OSL_OSEM) || (options.mlem && options.OSL_MLEM)) ...
        || OS_I4_summa > 1)
    error(['Implementation ' num2str(options.reconstruction_method) ' supports only one OS and one MLEM algorithm at a time.'])
end
if options.reconstruction_method == 1 && ~options.precompute_lor
    warning(['Implementation 1 without precomputation is NOT recommended as it is extremely memory demanding and slow! Either set '...
        'precompute_lor to true or use another implementation.'])
    if options.projector_type == 2
        warning('Orthogonal distance based projector is NOT recommended when using implementation 1 without precomputation!')
    end
end
if options.verbose
    if options.mlem
        if (options.reconstruction_method == 1 && ~options.precompute_obs_matrix) || options.reconstruction_method == 5
            warning('MLEM is not supported with implementation 5 or with implementation 1 without precomputed observation matrix.')
            options.mlem = false;
            if ~OS && ~MAPOS
                error('No other reconstruction algorithms selected. Select an ordered subsets algorithm.')
            end
        else
            disp('MLEM selected.')
        end
    end
    if options.osem
        disp('OSEM selected.')
    end
    if options.ramla
        disp('RAMLA selected.')
    end
    if options.mramla
        disp('MRAMLA selected.')
    end
    if options.rosem
        disp('ROSEM selected.')
    end
    if options.rbi
        disp('RBI selected.')
    end
    if options.drama
        disp('DRAMA selected.')
    end
    if options.cosem
        disp('COSEM selected.')
    end
    if options.acosem
        disp('ACOSEM selected.')
    end
    if options.ecosem
        disp('ECOSEM selected.')
    end
    if options.OSL_MLEM && PRIOR
        if options.reconstruction_method ~= 2
            warning('MLEM-OSL is not supported with implementations 1, 3 or 4.')
            options.OSL_MLEM = false;
        else
            disp('MLEM-OSL selected.')
        end
    elseif options.OSL_MLEM && ~PRIOR
        warning('MLEM-OSL selected, but no prior has been selected. No MAP reconstruction will be performed.')
        options.OSL_MLEM = false;
    end
    if options.OSL_OSEM && PRIOR
        disp('OSEM-OSL selected.')
    elseif options.OSL_OSEM && ~PRIOR
        warning('OSEM-OSL selected, but no prior has been selected. No MAP reconstruction will be performed.')
        options.OSL_OSEM = false;
    end
    if options.BSREM && PRIOR
        disp('BSREM selected.')
    elseif options.BSREM && ~PRIOR
        warning('BSREM selected, but no prior has been selected. No MAP reconstruction will be performed.')
        options.BSREM = false;
    end
    if options.MBSREM && PRIOR
        disp('MBSREM selected.')
    elseif options.MBSREM && ~PRIOR
        warning('MBSREM selected, but no prior has been selected. No MAP reconstruction will be performed.')
        options.MBSREM = false;
    end
    if options.ROSEM_MAP && PRIOR
        disp('ROSEM-MAP selected.')
    elseif options.ROSEM_MAP && ~PRIOR
        warning('ROSEM_MAP selected, but no prior has been selected. No MAP reconstruction will be performed.')
        options.ROSEM_MAP = false;
    end
    if options.RBI_MAP && PRIOR
        disp('RBI-MAP selected.')
    elseif options.RBI_MAP && ~PRIOR
        warning('RBI_MAP selected, but no prior has been selected. No MAP reconstruction will be performed.')
        options.RBI_MAP = false;
    end
    if any(options.COSEM_MAP) && PRIOR
        if options.COSEM_MAP == 1 && PRIOR
            disp('ACOSEM-OSL selected.')
        elseif options.COSEM_MAP == 1 && ~PRIOR
            warning('ACOSEM-OSL selected, but no prior has been selected. No MAP reconstruction will be performed.')
            options.COSEM_MAP = 0;
        elseif options.COSEM_MAP == 2 && PRIOR
            disp('COSEM-OSL selected.')
        elseif options.COSEM_MAP == 2 && ~PRIOR
            warning('COSEM-OSL selected, but no prior has been selected. No MAP reconstruction will be performed.')
            options.COSEM_MAP = 0;
        elseif options.COSEM_MAP > 2 || options.COSEM_MAP < 0
            error('Unsupported COSEM-MAP method selected!')
        end
    end
    if options.MRP && MAP
        disp('MRP selected.')
    elseif options.MRP && ~MAP
        warning('MRP selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.')
        options.MRP = false;
    end
    if options.quad && MAP
        disp('Quadratic prior selected.')
    elseif options.quad && ~MAP
        warning('Quadratic prior selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.')
        options.quad = false;
    end
    if options.L && MAP
        disp('L-filter selected.')
    elseif options.L && ~MAP
        warning('L-filter selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.')
        options.L = false;
    end
    if options.FMH && MAP
        disp('FMH selected.')
    elseif options.FMH && ~MAP
        warning('FMH selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.')
        options.FMH = false;
    end
    if options.weighted_mean && MAP
        if options.mean_type == 1
            disp('Weighted (arithmetic) mean selected.')
        elseif options.mean_type == 2
            disp('Weighted (harmonic) mean selected.')
        elseif options.mean_type == 3
            disp('Weighted (geometric) mean selected.')
        else
            error('Unsupported mean type selected.')
        end
    elseif options.weighted_mean && ~MAP
        warning('Weighted mean selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.')
        options.weighted_mean = false;
    end
    if options.TV && MAP
        if options.TV_use_anatomical
            if options.TVtype == 1
                disp('Anatomically weighted TV selected.')
            elseif options.TVtype == 2
                disp('Joint TV selected.')
            elseif options.TVtype == 3
                disp('Weighted joint TV selected.')
            else
                error('Unsupported TV type selected.')
            end
        else
            if options.TVtype == 1 || options.TVtype == 2
                disp('TV selected.')
            elseif options.TVtype == 3
                disp('Weighted TV selected.')
            else
                error('Unsupported TV type selected.')
            end
        end
    elseif options.TV && ~MAP
        warning('TV selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.')
        options.TV = false;
    end
    if options.AD && MAP
        if options.FluxType > 2 || options.FluxType < 1
            error('FluxType has to be either 1 or 2')
        end
        if (options.DiffusionType > 2 || options.DiffusionType < 1) && options.reconstruction_method == 2
            error('DiffusionType has to be either 1 or 2')
        end
        disp('AD selected.')
    elseif options.AD && ~MAP
        warning('AD selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.')
        options.AD = false;
    end
    if options.APLS && MAP
        disp('APLS selected.')
    elseif options.APLS && ~MAP
        warning('APLS selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.')
        options.APLS = false;
    end
    if options.TGV && MAP
        disp('TGV selected.')
    elseif options.TGV && ~MAP
        warning('TGV selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.')
        options.TGV = false;
    end
    if options.NLM && MAP
        if options.reconstruction_method == 2
            warning('NLM is not supported with OpenCL reconstruction')
            options.NLM = false;
        else
            warning('NLM selected. NLM is an experimental feature.')
        end
    elseif options.NLM && ~MAP
        warning('NLM selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.')
        options.NLM = false;
    end
    if ~OS && ~MAP && ~options.mlem && ~options.only_sinos
        error('No reconstruction algorithm selected')
    end
end