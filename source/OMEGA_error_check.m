function options = OMEGA_error_check(options)
%% Error checking file
% This function is used to check that all the input values are allowed. It
% also prints several variables that were chosen to inform the user of the 
%

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

% First parts checks if certain elements are missing from the struct and
% assigns default values to them
if ~isfield(options, 'custom')
    options.custom = false;
end
if ~isfield(options, 'store_raw_data')
    options.store_raw_data = options.use_raw_data;
end
if ~isfield(options, 'Huber')
    options.Huber = false;
end
if ~isfield(options, 'NLM')
    options.NLM = false;
end
if ~isfield(options, 'no_data_load')
    options.no_data_load = false;
end
if ~isfield(options, 'TOF_bins')
    options.TOF_bins = 1;
end
if ~isfield(options, 'TOF_FWHM')
    options.TOF_FWHM = 0;
end
if ~isfield(options, 'cryst_per_block_axial')
    options.cryst_per_block_axial = options.cryst_per_block;
end
if ~isfield(options, 'transaxial_multip')
    options.transaxial_multip = 1;
end
if ~isfield(options,'use_machine')
    options.use_machine = 0;
end

% Determine whether various different reconstruction modes are used (e.g.
% MAP reconstruction and any prior)
MAP = (options.OSL_MLEM || options.OSL_OSEM || options.BSREM || options.MBSREM || options.ROSEM_MAP || options.RBI_OSL || any(options.COSEM_OSL));
MAPOS = (options.OSL_OSEM || options.BSREM || options.MBSREM || options.ROSEM_MAP || options.RBI_OSL || any(options.COSEM_OSL));
PRIOR = (options.MRP || options.quad || options.Huber || options.L || options.FMH || options.weighted_mean || options.TV || options.AD || options.APLS ...
    || options.TGV || options.NLM || options.custom);
PRIOR_summa = sum([options.MRP, options.quad, options.Huber, options.L, options.FMH, options.weighted_mean, options.TV, options.AD, options.APLS, ...
    options.TGV, options.NLM, options.custom]);
OS = (options.osem || options.ramla || options.mramla || options.cosem || options.rosem || options.rbi || options.ecosem || options.acosem || options.drama || ...
    options.OSL_OSEM || options.BSREM || options.MBSREM || options.ROSEM_MAP || options.RBI_OSL || any(options.COSEM_OSL));
MLOS = (options.mlem || options.osem);
OS_I3 = (options.ramla || options.mramla || options.cosem || options.rosem || options.rbi || options.ecosem || options.acosem || options.drama || ...
    options.OSL_OSEM || options.BSREM || options.MBSREM || options.ROSEM_MAP || options.RBI_OSL || any(options.COSEM_OSL));
NMLOS = (options.ramla || options.mramla || options.cosem || options.rosem || options.rbi || options.ecosem || options.acosem || options.drama || ...
    options.OSL_OSEM || options.BSREM || options.MBSREM || options.ROSEM_MAP || options.RBI_OSL || any(options.COSEM_OSL) || options.OSL_MLEM);
% OS_I4 = (options.mramla || options.cosem || options.rbi || options.ecosem || options.acosem || options.MBSREM || options.RBI_OSL || any(options.COSEM_OSL));
OS_I4_summa = sum([options.osem, options.ramla, options.rosem, options.drama, options.rbi, options.cosem, options.ecosem, options.acosem, ...
    options.OSL_OSEM, options.BSREM, options.ROSEM_MAP, options.RBI_OSL, any(options.COSEM_OSL)]);
N_PRIORS = (options.MRP + options.quad + options.Huber + options.L + options.FMH + options.weighted_mean + options.TV + options.AD + options.APLS ...
    + options.TGV + options.NLM + options.custom);

% Check for various illegal values
if options.FOVa_x >= options.diameter || options.FOVa_y >= options.diameter
    error(['Transaxial FOV is larger than the machine diameter (' num2str(options.diameter) ')'])
end
if (options.axial_fov) < (options.rings * options.cr_pz - options.cr_pz)
    error('Axial FOV is too small, crystal ring(s) on the boundary have no slices')
end
% if (options.axial_fov) > (options.linear_multip * options.cryst_per_block * options.cr_pz + options.axial_fov/options.Nz*2 + options.cr_pz*sum(options.pseudot))
%     error('Axial FOV is too large, not all the slices have LORs going through them')
% end
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
if options.span > options.ring_difference && options.NSinos > 1 && ~options.use_raw_data
    error(['Span value cannot be larger than ring difference (' num2str(options.ring_difference) ')'])
end
if options.span == 1 && ~options.use_raw_data
    warning('Span value of 1 is not recommended. Use raw data if you want uncompressed reconstruction!')
end
if (mod(options.span,2) == 0 || options.span <= 0) && ~options.use_raw_data
    error('Span value has to be odd and positive.')
end
if options.ring_difference >= options.rings && ~options.use_raw_data
    error(['Ring difference can be at most ' num2str(options.rings-1)])
end
if options.ring_difference < 0 && ~options.use_raw_data
    error('Ring difference has to be at least 0.')
end
if options.Nang > options.det_w_pseudo/2 && ~options.use_raw_data
    error(['Number of sinogram angles can be at most the number of detectors per ring divided by two(' num2str(options.det_w_pseudo/2) ')'])
end
if options.TotSinos < options.NSinos && ~options.use_raw_data
    error(['The numnber of sinograms used (' num2str(options.NSinos) ') is larger than the total number of sinograms (' num2str(options.TotSinos) ')'])
end
if (options.ndist_side > 1 && mod(options.Ndist,2) == 0 || options.ndist_side < -1 && mod(options.Ndist,2) == 0) && ~options.use_raw_data
    error('ndist_side can be either 1 or -1')
end
if options.ndist_side == 0 && mod(options.Ndist,2) == 0 && ~options.use_raw_data
    error('ndist_side cannot be 0 when Ndist is even')
end
if (mod(options.sampling, 2) > 0 && options.sampling ~= 1) || options.sampling < 0
    error('Sampling rate has to be divisible by two and positive or one')
end
if options.sampling > 1 && options.precompute_lor
    warning('Increased sampling rate is not supported for precomputed data')
end
if options.arc_correction && options.use_raw_data
    warning('Arc correction is not supported for raw data')
end
if options.arc_correction && options.precompute_lor
    warning('Arc correction is not supported with precomputed data')
    options.arc_correction = false;
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
% if options.subsets < 2 && OS
%     warning('Number of subsets is less than two. Subset has to be at least 2 when using OS-methods. Using 2 subsets.')
%     options.subsets = 2;
% end
if options.det_per_ring == options.det_w_pseudo && options.fill_sinogram_gaps
    error('Gap filling is only supported with pseudo detectors!')
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
if isunix
    if length(options.fpath) > 1 && ~strcmp('/',options.fpath(end))
        options.fpath = [options.fpath '/'];
    end
elseif ispc
    if length(options.fpath) > 1 && ~strcmp('\',options.fpath(end)) && ~strcmp('/',options.fpath(end))
        options.fpath = [options.fpath '\'];
    end
else
    if length(options.fpath) > 1 && ~strcmp('/',options.fpath(end))
        options.fpath = [options.fpath '/'];
    end
end
if options.use_LMF && options.randoms_correction
    warning('Randoms correction is set to true although LMF input is selected. No randoms correction will be performed.')
    options.randoms_correction = false;
end
if options.TV_use_anatomical && options.TV && exist(options.TV_reference_image,'file') ~= 2 && MAP
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
if options.implementation == 1 && ~options.precompute_lor && (options.n_rays_transaxial > 1 || options.n_rays_axial > 1)
    error('Multiray Siddon is not supported with implementation 1')
end
if options.implementation == 3
    if options.osem && options.mlem
        warning('Both OSEM and MLEM selected with implementation 3, using only MLEM');
        options.osem = false;
    end
    if options.mlem && options.subsets > 1 && (options.implementation == 3 || options.implementation == 1)
        options.subsets = 1;
    end
    if options.osem && options.subsets == 1
        warning('OSEM selected with only 1 subset. Switching to MLEM instead')
        options.osem = false;
        options.mlem = true;
    end
end
if options.use_machine == 2 && options.use_raw_data
    warning('Sinogram data cannot be used when raw data is set to true, using list-mode data instead')
    options.use_machine = 1;
end
% if (options.blank || options.transmission_scan) && options.attenuation_correction
%     error('Both transmission and image based attenuation correction selected. Use only one correction method.')
% end
% if options.source && options.use_ASCII && (options.source_index1 == 0 || isempty(options.source_index1) || options.source_index2 == 0 || isempty(options.source_index2))
%     error('Source image selected with ASCII data, but no source index column numbers are provided.')
% end
if options.reconstruct_trues && options.reconstruct_scatter
    warning('Both reconstruct trues and scatter selected, reconstructing only trues.')
    options.reconstruct_scatter = false;
end
if options.implementation == 1 && exist('projector_mex','file') ~= 3 && options.precompute_lor
    error('MEX-file for implementation 1 not found. Run install_mex first.')
end
if options.implementation == 4 && exist('projector_mex','file') ~= 3
    error('MEX-file for implementation 4 not found. Run install_mex first.')
end
if exist('OCTAVE_VERSION','builtin') == 0 && options.use_root && ((exist('GATE_root_matlab','file') ~= 3 && ~verLessThan('matlab', '9.6')) || exist('GATE_root_matlab_C','file') ~= 3) && options.use_machine == 0
    warning(['ROOT selected, but no MEX-file for ROOT data load found. Run install_mex to build ROOT MEX-file. Ignore this warning if you are ' ...
        'simply loading a mat-file containing measurement data from ROOT files.'])
end
if options.use_root && exist('GATE_root_matlab_oct','file') ~= 3 && options.use_machine == 0 && exist('OCTAVE_VERSION','builtin') == 5
    warning(['ROOT selected, but no OCT-file for ROOT data load found. Run install_mex to build ROOT OCT-file. Ignore this warning if you are ' ...
        'simply loading a mat-file containing measurement data from ROOT files.'])
end
if options.use_LMF && exist('gate_lmf_matlab','file') ~= 3 && options.use_machine == 0
    error('LMF selected, but no MEX-file for LMF data load found. Run install_mex to build LMF MEX-file.')
end
if options.implementation == 2 && exist('OpenCL_matrixfree','file') ~= 3
    error('OpenCL reconstruction selected, but OpenCL MEX-files were not installed. Run install_mex to build OpenCL MEX-files.')
end
if options.implementation == 3 && exist('OpenCL_matrixfree_multi_gpu','file') ~= 3
    error('OpenCL reconstruction selected, but OpenCL MEX-files were not installed. Run install_mex to build OpenCL MEX-files.')
end
if options.implementation == 5 && exist('improved_Siddon_openCL','file') ~= 3
    error('OpenCL reconstruction selected, but OpenCL MEX-files were not installed. Run install_mex to build OpenCL MEX-files.')
end
if options.implementation == 3 && NMLOS
    warning(['Implementation ' num2str(options.implementation) ' selected with reconstruction algorithms other than MLEM or OSEM. '...
        'Only MLEM or OSEM reconstruction can be performed'])
end
if options.implementation == 3 && ~MLOS
    error(['Implementation ' num2str(options.implementation) ' selected, but neither MLEM nor OSEM algorithm has been selected.'])
end
if options.implementation == 3 && OS_I3
    warning(['Implementation ' num2str(options.implementation) ' supports only MLEM and OSEM, any other algorithms will be ignored.'])
end
if options.implementation == 4 && (options.mramla || options.MBSREM)
    error(['Implementation ' num2str(options.implementation) ' selected with unsupported algorithm. MRAMLA or MBSREM are not supported!'])
end
if options.implementation == 4 && PRIOR_summa > 1 && MAP
    error(['Implementation ' num2str(options.implementation) ' supports only one prior at a time.'])
end
if options.implementation == 4 && (PRIOR_summa == 1 && ((options.mlem && options.OSL_MLEM)) || OS_I4_summa > 1)
    error(['Implementation ' num2str(options.implementation) ' supports only one OS and one MLEM algorithm at a time.'])
end
if options.implementation == 1 && ~options.precompute_lor
    if options.projector_type == 2 || options.projector_type == 3
        error('Orthogonal distance-based/volume-based projector is NOT supported when using implementation 1 without precomputation!')
    end
    warning(['Implementation 1 without precomputation is NOT recommended as it is extremely memory demanding and slow! It is highly recommended to either set '...
        'precompute_lor to true or use another implementation.'])
end
% Print various options that were selected if verbosity has been enabled
if options.verbose
    if options.use_ASCII && options.use_machine == 0
        disp('Using ASCII data.')
    elseif options.use_LMF && options.use_machine == 0
        disp('Using LMF data.')
    elseif options.use_root && options.use_machine == 0
        disp('Using ROOT data.')
    elseif options.use_machine == 1
        disp('Using data obtained from list-mode file.')
    elseif options.use_machine == 2
        disp('Using machine created sinogram data.')
    elseif options.use_machine == 3
        disp('Using 32-bit list-mode data.')
    end
    if options.only_sinos
        disp('Loading only data.')
    end
    if ~(options.compute_normalization || options.only_sinos)
        dispaus = (['Using implementation ' num2str(options.implementation)]);
        if options.implementation == 1 || options.implementation == 4
            dispaus = [dispaus, '.'];
        elseif options.implementation == 2
            inffo = ArrayFire_OpenCL_device_info();
            k2 = strfind(inffo, 'MB');
            if options.use_device == 0
                k1 = strfind(inffo, '[0]');
                dispaus = [dispaus, ' with ', inffo(k1 + 4:k2(1)+1)];
            else
                k1 = strfind(inffo, ['-' num2str(options.use_device) '-']);
                dispaus = [dispaus, ' with ', inffo(k1 + 4:k2(options.use_device + 1)+1)];
            end
        elseif options.implementation == 3
            inffo = OpenCL_device_info();
            k6 = strfind(inffo, 'No. platforms: ');
            k2 = strfind(inffo, 'MB');
            k1 = strfind(inffo, 'Platform') + 12;
            k3 = strfind(inffo, 'Device');
            k4 = strfind(inffo, 'No. devices:');
            n_platforms = str2double(inffo(k6+15:k1(1)-15));
            n_devices = [];
            hh = 1;
            for kk = 1 : n_platforms
                n_devices = [n_devices ; str2double(inffo(k4(kk)+14:k3(hh)-2))];
                hh = hh + n_devices(kk);
            end
            k5 = strfind(inffo, 'with ');
            MB = [];
            hh = 1;
            for kk = 1 : n_platforms
                if kk == options.use_device + 1
                    for ll = 1 : n_devices(kk)
                        MB = [MB; str2double(inffo(k5(hh + ll - 1)+5:k2(hh + ll - 1)-2))];
                    end
                end
                hh = hh + n_devices(kk);
            end
            kumulatiivinen = cumsum(n_devices);
            if options.cpu_to_gpu_factor == 0
                [~, loc] = max(MB);
                ll = kumulatiivinen(options.use_device + 1) - n_devices(options.use_device + 1) + loc;
                if n_devices(options.use_device + 1) > 1
                    dispaus = [dispaus, ' with ', inffo(k1(options.use_device + 1):k4(options.use_device + 1) - 2), ...
                        sprintf('\n'), inffo(k3(ll):k3(ll + 1) - 2)];
                else
                    if options.use_device + 1 == n_platforms
                        dispaus = [dispaus, ' with ', inffo(k1(options.use_device + 1):k4(options.use_device + 1) - 2), ...
                            sprintf('\n'), inffo(k3(kumulatiivinen(end)):end - 1)];
                    else
                        dispaus = [dispaus, ' with ', inffo(k1(options.use_device + 1):k4(options.use_device + 1) - 2), ...
                            sprintf('\n'), inffo(k3(kumulatiivinen(options.use_device + 1)):k1(options.use_device + 2) - 15)];
                    end
                end
            else
                loc = MB > 2048;
                if n_devices(options.use_device + 1) > 1
                    dispaus = [dispaus, ' with ', inffo(k1(options.use_device + 1):k4(options.use_device + 1) - 2)];
                    for hh = 1 : n_devices(options.use_device + 1)
                        if loc(hh) == 0
                            continue;
                        end
                        ll = kumulatiivinen(options.use_device + 1) - n_devices(options.use_device + 1) + hh;
                        dispaus = [dispaus, sprintf('\n'), inffo(k3(ll):k3(ll + 1) - 2)];
                    end
                else
                    if options.use_device + 1 == n_platforms
                        dispaus = [dispaus, ' with ', inffo(k1(options.use_device + 1):k4(options.use_device + 1) - 2), ...
                            sprintf('\n'), inffo(k3(kumulatiivinen(end)):end - 1)];
                    else
                        dispaus = [dispaus, ' with ', inffo(k1(options.use_device + 1):k4(options.use_device + 1) - 2), ...
                            sprintf('\n'), inffo(k3(kumulatiivinen(options.use_device + 1)):k1(options.use_device + 2) - 15)];
                    end
                end
            end
        end
        disp(dispaus);
        reko = {};
        if options.mlem
            if (options.implementation == 1 && ~options.precompute_obs_matrix) || options.implementation == 5
                warning('MLEM is not supported with implementation 1 or with implementation 5 without precomputed observation matrix.')
                options.mlem = false;
                if ~OS && ~MAPOS
                    error('No other reconstruction algorithms selected. Select an ordered subsets algorithm.')
                end
            else
                if ~OS || (options.implementation ~= 2 && options.implementation ~= 4)
                    options.subsets = 1;
                end
                reko = [reko;{'MLEM'}];
            end
        end
        if options.osem
            reko = [reko;{'OSEM'}];
        end
        if options.ramla
            reko = [reko;{'RAMLA'}];
        end
        if options.mramla
            reko = [reko;{'MRAMLA'}];
        end
        if options.rosem
            reko = [reko;{'ROSEM'}];
        end
        if options.rbi
            reko = [reko;{'RBI'}];
        end
        if options.drama
            reko = [reko;{'DRAMA'}];
        end
        if options.cosem
            reko = [reko;{'COSEM'}];
        end
        if options.acosem
            reko = [reko;{'ACOSEM'}];
        end
        if options.ecosem
            reko = [reko;{'ECOSEM'}];
        end
        if options.OSL_MLEM && PRIOR
            if options.implementation ~= 2
                warning('MLEM-OSL is not supported with implementations 1 and 3.')
                options.OSL_MLEM = false;
            else
                reko = [reko;{'MLEM-OSL'}];
            end
        elseif options.OSL_MLEM && ~PRIOR
            warning('MLEM-OSL selected, but no prior has been selected. No MAP reconstruction will be performed.')
            options.OSL_MLEM = false;
        end
        if options.OSL_OSEM && PRIOR
            reko = [reko;{'OSEM-OSL'}];
        elseif options.OSL_OSEM && ~PRIOR
            warning('OSEM-OSL selected, but no prior has been selected. No MAP reconstruction will be performed.')
            options.OSL_OSEM = false;
        end
        if options.BSREM && PRIOR
            reko = [reko;{'BSREM'}];
        elseif options.BSREM && ~PRIOR
            warning('BSREM selected, but no prior has been selected. No MAP reconstruction will be performed.')
            options.BSREM = false;
        end
        if options.MBSREM && PRIOR
            reko = [reko;{'MBSREM'}];
        elseif options.MBSREM && ~PRIOR
            warning('MBSREM selected, but no prior has been selected. No MAP reconstruction will be performed.')
            options.MBSREM = false;
        end
        if options.ROSEM_MAP && PRIOR
            reko = [reko;{'ROSEM-MAP'}];
        elseif options.ROSEM_MAP && ~PRIOR
            warning('ROSEM_MAP selected, but no prior has been selected. No MAP reconstruction will be performed.')
            options.ROSEM_MAP = false;
        end
        if options.RBI_OSL && PRIOR
            reko = [reko;{'RBI-OSL'}];
        elseif options.RBI_OSL && ~PRIOR
            warning('RBI_OSL selected, but no prior has been selected. No MAP reconstruction will be performed.')
            options.RBI_OSL = false;
        end
        if any(options.COSEM_OSL) && PRIOR
            if options.COSEM_OSL == 1 && PRIOR
                reko = [reko;{'ACOSEM-OSL'}];
            elseif options.COSEM_OSL == 1 && ~PRIOR
                warning('ACOSEM-OSL selected, but no prior has been selected. No MAP reconstruction will be performed.')
                options.COSEM_OSL = 0;
            elseif options.COSEM_OSL == 2 && PRIOR
                reko = [reko;{'COSEM-OSL'}];
            elseif options.COSEM_OSL == 2 && ~PRIOR
                warning('COSEM-OSL selected, but no prior has been selected. No MAP reconstruction will be performed.')
                options.COSEM_OSL = 0;
            elseif options.COSEM_OSL > 2 || options.COSEM_OSL < 0
                error('Unsupported COSEM-OSL method selected!')
            end
            if N_PRIORS > 1
                error('Only one prior can be used at a time with COSEM-OSL')
            end
        end
        for kk = 1 : length(reko)
            if kk == 1
                dispi = reko{kk};
            else
                dispi = [dispi, ', ' reko{kk}];
            end
            if kk == length(reko) && length(reko) > 1
                dispi = [dispi, ' reconstruction methods selected.'];
            elseif kk == length(reko) && length(reko) == 1
                dispi = [dispi, ' reconstruction method selected.'];
            end
        end
        if ~isempty(reko)
            disp(dispi)
        else
            error('No reconstruction method selected')
        end
        priori = {};
        if options.MRP && MAP
            priori = [priori;{'MRP'}];
        elseif options.MRP && ~MAP
            warning('MRP selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.')
            options.MRP = false;
        end
        if options.quad && MAP
            priori = [priori;{'Quadratic'}];
        elseif options.quad && ~MAP
            warning('Quadratic prior selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.')
            options.quad = false;
        end
        if options.Huber && MAP
            priori = [priori;{'Huber'}];
        elseif options.Huber && ~MAP
            warning('Huber prior selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.')
            options.Huber = false;
        end
        if options.L && MAP
            priori = [priori;{'L-filter'}];
        elseif options.L && ~MAP
            warning('L-filter selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.')
            options.L = false;
        end
        if options.FMH && MAP
            priori = [priori;{'FMH'}];
        elseif options.FMH && ~MAP
            warning('FMH selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.')
            options.FMH = false;
        end
        if options.weighted_mean && MAP
            if options.mean_type == 1 || options.mean_type == 4
                priori = [priori;{'Weighted (arithmetic) mean'}];
            elseif options.mean_type == 2 || options.mean_type == 5
                priori = [priori;{'Weighted (harmonic) mean'}];
            elseif options.mean_type == 3 || options.mean_type == 6
                priori = [priori;{'Weighted (geometric) mean'}];
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
                    priori = [priori;{'Anatomically weighted TV'}];
                elseif options.TVtype == 2
                    priori = [priori;{'Joint TV'}];
                elseif options.TVtype == 3
                    priori = [priori;{'Weighted joint TV'}];
                else
                    error('Unsupported TV type selected.')
                end
            else
                if options.TVtype == 1 || options.TVtype == 2
                    priori = [priori;{'TV'}];
                elseif options.TVtype == 3
                    priori = [priori;{'Weighted TV'}];
                elseif options.TVtype == 4
                    priori = [priori;{'SATV'}];
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
                error('FluxType has to be either 1 or 2.')
            end
            if (options.DiffusionType > 2 || options.DiffusionType < 1) && options.implementation == 2
                error('DiffusionType has to be either 1 or 2.')
            end
            priori = [priori;{'AD-MRP'}];
        elseif options.AD && ~MAP
            warning('AD-MRP selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.')
            options.AD = false;
        end
        if options.APLS && MAP
            priori = [priori;{'APLS'}];
        elseif options.APLS && ~MAP
            warning('APLS selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.')
            options.APLS = false;
        end
        if options.TGV && MAP
            priori = [priori;{'TGV'}];
        elseif options.TGV && ~MAP
            warning('TGV selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.')
            options.TGV = false;
        end
        if options.NLM && MAP
            priori = [priori;{'NLM'}];
        elseif options.NLM && ~MAP
            warning('NLM selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.')
            options.NLM = false;
        end
        if isfield(options, 'custom') && options.custom && MAP
            priori = [priori;{'Custom'}];
        elseif ~MAP && isfield(options, 'custom') && options.custom
            error('Custom prior selected, but no MAP-method selected!')
        end
        for kk = 1 : length(priori)
            if kk == 1
                dispi2 = priori{kk};
            else
                dispi2 = [dispi2, ', ' priori{kk}];
            end
            if kk == length(priori) && length(priori) > 1
                dispi2 = [dispi2, ' priors selected.'];
            elseif kk == length(priori) && length(priori) == 1
                dispi2 = [dispi2, ' prior selected.'];
            end
        end
        if ~isempty(priori)
            disp(dispi2)
        end
        if ~OS && ~MAP && ~options.mlem && ~(options.only_sinos || options.compute_normalization)
            error('No reconstruction algorithm selected.')
        end
        if options.precompute_lor
            disp('Precomputed LOR voxel counts used.')
        else
            disp('No precomputed data will be used.')
        end
        if options.projector_type == 1 && ~options.precompute_lor
            if options.implementation == 1
                disp('Improved Siddon''s algorithm used with 1 ray.')
            else
                if options.n_rays_transaxial > 1
                    ray = 'rays';
                else
                    ray = 'ray';
                end
                if options.n_rays_axial > 1
                    aray = 'rays';
                else
                    aray = 'ray';
                end
                disp(['Improved Siddon''s algorithm used with ' num2str(options.n_rays_transaxial) ' transaxial ' ray ' and ' ...
                    num2str(options.n_rays_axial) ' axial ' aray '.'])
            end
        elseif options.projector_type == 1
            disp('Improved Siddon''s algorithm used with 1 ray.')
        elseif options.projector_type == 2
            dispi = 'Orthogonal distance-based ray tracer used';
            if options.tube_width_z > 0
                dispi = [dispi, ' in 3D mode.'];
            else
                dispi = [dispi, ' in 2.5D mode.'];
            end
            disp(dispi)
        elseif options.projector_type == 3
            disp('Volume of intersection based ray tracer used.');
        end
        if options.use_psf
            if options.deblurring
                disp('PSF ON with deblurring phase.');
            else
                disp('PSF ON.');
            end
            
        end
        if options.attenuation_correction
            disp('Attenuation correction ON.')
        end
        if options.randoms_correction
            dispi = 'Randoms correction ON';
            if options.variance_reduction
                dispi = [dispi, ' with variance reduction'];
                if options.randoms_smoothing
                    dispi = [dispi, ' and smoothing'];
                end
            elseif options.randoms_smoothing
                dispi = [dispi, ' with smoothing'];
            end
            dispi = [dispi, '.'];
            disp(dispi)
        end
        if options.scatter_correction
            dispi = 'Scatter correction ON';
            if options.scatter_smoothing
                dispi = [dispi, ' with smoothing'];
            end
            dispi = [dispi, '.'];
            disp(dispi)
        end
        if options.fill_sinogram_gaps
            disp('Sinogram gap filling ON.')
        end
    end
    if options.normalization_correction && ~(options.compute_normalization || options.only_sinos)
        disp('Normalization correction ON.')
    elseif options.normalization_correction && options.compute_normalization
        warning('Normalization correction cannot be applied when computing normalization coefficients. Disabling normalization correction.')
        options.normalization_correction = false;
    elseif options.compute_normalization
        disp('Computing normalization coefficients.')
    end
    if options.compute_normalization && sum(options.normalization_options) == 0
        error('Normalization computation selected, but no normalization components selected.')
    end
    if ~(options.compute_normalization || options.only_sinos)
        if options.corrections_during_reconstruction && (options.normalization_correction || options.randoms_correction || options.scatter_correction)
            disp('Corrections applied during reconstruction (ordinary Poisson).')
        elseif ~options.corrections_during_reconstruction && (options.normalization_correction || options.randoms_correction || options.scatter_correction)
            disp('Corrections applied to the measurement data.')
        end
        if options.arc_correction && ~options.use_raw_data
            disp('Arc correction ON.')
        end
        if options.use_raw_data
            if options.partitions == 1
                dispi = 'Using STATIC raw list-mode data';
            else
                dispi = 'Using DYNAMIC raw list-mode data';
            end
            if options.reconstruct_trues
                dispi = strcat(dispi, ' (trues)');
            elseif options.reconstruct_scatter
                dispi = strcat(dispi, ' (scatter)');
            else
                dispi = strcat(dispi, ' (prompts)');
            end
            if options.sampling_raw > 1
                dispi = strcat(dispi, [' with ' num2str(options.sampling_raw) 'x sampling']);
            end
            if options.partitions > 1
                if options.sampling_raw > 1
                    dispi = strcat(dispi, [' and ' num2str(options.partitions) ' time steps']);
                else
                    dispi = strcat(dispi, [' with ' num2str(options.partitions) ' time steps']);
                end
            end
            dispi = strcat(dispi, '.');
            disp(dispi)
        else
            if options.partitions == 1
                dispi = 'Using STATIC sinogram data';
            else
                dispi = 'Using DYNAMIC sinogram data';
            end
            if options.reconstruct_trues
                dispi = strcat(dispi, ' (trues)');
            elseif options.reconstruct_scatter
                dispi = strcat(dispi, ' (scatter).');
            else
                dispi = strcat(dispi, ' (prompts)');
            end
            if options.sampling > 1
                dispi = strcat(dispi, [' with ' num2str(options.sampling) 'x sampling']);
            end
            if options.partitions > 1
                if options.sampling > 1
                    dispi = strcat(dispi, [' and ' num2str(options.partitions) ' time steps']);
                else
                    dispi = strcat(dispi, [' with ' num2str(options.partitions) ' time steps']);
                end
            end
            dispi = strcat(dispi, '.');
            disp(dispi)
        end
        disp(['Using an image (matrix) size of ' num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) ' with ' num2str(options.Niter) ...
            ' iterations and ' num2str(options.subsets) ' subsets.'])
        if options.use_CUDA && options.projector_type > 1
            warning('CUDA is not recommended with orthogonal or volume-based projectors')
        end
    end
end