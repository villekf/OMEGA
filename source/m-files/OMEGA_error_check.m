function options = OMEGA_error_check(options)
%% Error checking file
% This function is used to check that all the input values are allowed. It
% also prints several variables that were chosen to inform the user of the
% selected options.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019-2025 Ville-Veikko Wettenhovi
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
options = setMissingValues(options);
options = convertOptions(options);

% Determine whether various different reconstruction modes are used (e.g.
% MAP reconstruction and any prior)
MAP = checkAlgorithmsPriors(options, 2);
PRIOR = checkAlgorithmsPriors(options, 1);
PRIOR_summa = checkAlgorithmsPriors(options, 1, 0);
OS = checkAlgorithmsPriors(options, 6);
MLOS = options.OSEM;
NMLOS = checkAlgorithmsPriors(options, 9);
OS_I4_summa = checkAlgorithmsPriors(options, 2,0) + checkAlgorithmsPriors(options, 6, 0);
preCondImAlg = checkAlgorithmsPriors(options, 7);
preCondMeasAlg = checkAlgorithmsPriors(options, 8);

if numel(options.partitions) > 1
    partitions = numel(options.partitions);
else
    partitions = options.partitions;
end
if options.only_sinos && options.only_reconstructions
    error('options.only_sinos and options.only_reconstructions cannot be both set to true')
end
% Check for various illegal values
if ~options.CT && ~options.SPECT && (options.FOVa_x >= options.diameter || options.FOVa_y >= options.diameter)
    warning(['Transaxial FOV is larger than the scanner diameter (' num2str(options.diameter) ')!'])
end
if ~options.CT && ~options.SPECT && (options.axial_fov) < (options.rings(1) * options.cr_pz - options.cr_pz)
    warning('Axial FOV is too small, crystal ring(s) on the boundary have no slices!')
end
if ~preCondImAlg && any(options.precondTypeImage)
    warning('Image-based preconditioning selected, but the selected algorithm(s) do not support preconditioning. No preconditioning will be performed.')
    disp('Supported algorithms are:')
    disp(recNames(7))
    options.precondTypeImage = [false;false;false;false;false;false];
end
if sum(options.precondTypeImage(1:3)) > 1
    error('Only one of the first 3 image-based preconditioners can be selected at a time!')
end
if ~preCondMeasAlg && any(options.precondTypeMeas)
    warning('Measurement-based preconditioning selected, but the selected algorithm does not support preconditioning. No preconditioning will be performed.')
    disp('Supported algorithms are:')
    disp(recNames(8))
    options.precondTypeMeas = [false;false];
end
if ~options.CT && ~options.SPECT && options.use_LMF && options.data_bytes < 10
    error('Too little data bytes in LMF format, minimum allowed is 10 bytes (time + detector indices)!')
end
if ~options.CT && ~options.SPECT && options.use_LMF && options.data_bytes > 21
    warning(['LMF format uses more bytes than the supported 21 bytes (time + detector indices + source coordinates + event indices + Compton scattering in phantom). '...
        'If these extra bytes are before the bytes that are used, output data will be incorrect.'])
end
if ~options.CT && ~options.SPECT && options.use_LMF && options.R_bits + options.C_bits + options.M_bits + options.S_bits + options.L_bits > 16
    error('Number of bits used in LMF is more than 16 bits. OMEGA supports only 16 bit detector indices!')
end
if ~options.CT && ~options.SPECT && options.span > options.ring_difference && options.NSinos > 1 && ~options.use_raw_data
    error(['Span value cannot be larger than ring difference (' num2str(options.ring_difference) ')!'])
end
if ~options.CT && ~options.SPECT && (mod(options.span,2) == 0 || options.span <= 0) && ~options.use_raw_data
    error('Span value has to be odd and positive!')
end
if ~options.CT && ~options.SPECT && options.ring_difference >= options.rings(end) && ~options.use_raw_data
    warning(['Ring difference can be at most ' num2str(options.rings-1) '. Setting it to the maximum possible.'])
    options.ring_difference = options.rings - 1;
end
if ~options.CT && ~options.SPECT && options.ring_difference < 0 && ~options.use_raw_data
    error('Ring difference has to be at least 0!')
end
if ~options.CT && ~options.SPECT && options.Nang(end) > options.det_w_pseudo(end)/2 && ~options.use_raw_data
    error(['Number of sinogram angles can be at most the number of detectors per ring divided by two(' num2str(options.det_w_pseudo(1)/2) ')!'])
end
if ~options.CT && ~options.SPECT && options.TotSinos < options.NSinos && ~options.use_raw_data
    error(['The numnber of sinograms used (' num2str(options.NSinos) ') is larger than the total number of sinograms (' num2str(options.TotSinos) ')!'])
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
if options.SPECT && options.implementation == 1
    error('Implementation 1 is not supported with SPECT data.')
end
if options.SPECT && options.implementation == 3 && options.projector_type == 6
    error('Implementation 3 is not supported with projector type 6.')
end
if options.SPECT && options.implementation == 5 && options.projector_type == 6
    error('Implementation 5 is not supported with projector type 6.')
end
if ~options.CT && ~options.SPECT && options.arc_correction && options.use_raw_data
    warning('Arc correction is not supported for raw data. Disabling arc correction.')
    options.arc_correction = false;
end
if ~options.CT && ~options.SPECT && options.arc_correction && options.implementation == 1
    warning('Arc correction is not supported with implementation 1. Disabling arc correction.')
    options.arc_correction = false;
end
if partitions < 1
    warning('Number of time steps is less than one. Using one time step.')
    options.partitions = 1;
end
if options.start > options.end
    error('Start time is later than end time!')
end
if options.start > options.tot_time
    error('Start time is larger than the total time of the measurement!')
end
if options.Niter < 1
    warning('Number of iterations is less than one! Using one iteration.')
    options.Niter = 1;
end
if isfield(options, 'maskFP') && numel(options.maskFP) > 1 && (options.subset_type == 3 || options.subset_type == 6 || options.subset_type == 7)
    error('Forward projection mask is only available with subset types 1, 2, 4, 5, and >= 8')
end
if options.implementation == 1 && options.useMultiResolutionVolumes
    error('Multi-resolution reconstruction is not supported with implementation 1!')
end
if options.useMultiResolutionVolumes && ~options.useEFOV
    warning('Multi-resolution reconstruction selected, but extended FOV is not selected! Disabling multi-resolution volumes.')
    options.useMultiResolutionVolumes = false;
end
if (options.LSQR || options.CGLS) && options.useMultiResolutionVolumes
    error('Multi-resolution reconstruction is not supported with LSQR or CGLS!')
end
if options.FDK && options.storeFP
    warning('Forward projections cannot be stored with FDK/FBP!')
end
if ~options.CT && ~options.SPECT && options.det_per_ring(1) == options.det_w_pseudo(1) && options.fill_sinogram_gaps
    error('Gap filling is only supported with pseudo detectors!')
end
if ~options.largeDim && isempty(options.x0)
    if ~options.CT
        warning('Initial value is an empty array, using the default values (1)')
        options.x0 = ones(options.Nx, options.Ny, options.Nz);
    else
        warning('Initial value is an empty array, using the default values (1e-4)')
        options.x0 = ones(options.Nx, options.Ny, options.Nz) * 1e-4;
    end
end
if ~options.largeDim && size(options.x0,1)*size(options.x0,2)*size(options.x0,3) < options.Nx*options.Ny*options.Nz
    error(['Initial value has a matrix size smaller (' num2str(size(options.x0,1)*size(options.x0,2)*size(options.x0,3)) ') than the actual image size ('...
        num2str(options.Nx*options.Ny*options.Nz) ')!'])
end
if ~options.largeDim && size(options.x0,1)*size(options.x0,2)*size(options.x0,3) > options.Nx*options.Ny*options.Nz
    warning(['Initial value has a matrix size larger (' num2str(size(options.x0,1)*size(options.x0,2)*size(options.x0,3)) ') than the actual image size (' ...
        num2str(options.Nx*options.Ny*options.Nz) ')! Attempting automatic resize.'])
    apuVar = numel(options.x0) / (options.Nx * options.Ny * options.Nz);
    apuVar = (apuVar^(1/3));
    options.x0 = imresize3(reshape(options.x0, round(options.Nx * apuVar), round(options.Ny * apuVar), ceil(options.Nz * apuVar)), [options.Nx, options.Ny, options.Nz]);
end
if ~options.CT && ~options.SPECT && options.use_LMF && options.randoms_correction
    warning('Randoms correction is set to true although LMF input is selected. No randoms correction will be performed.')
    options.randoms_correction = false;
end
if options.TV_use_anatomical && options.TV && exist(options.TV_reference_image,'file') ~= 2 && MAP && exist(options.TV_reference_image,'var') ~= 1
    error('Anatomical reference image for TV was not found on the specified path!')
end
if options.NLM_use_anatomical && options.NLM && MAP && exist(options.NLM_reference_image,'file') ~= 2 && exist(options.NLM_reference_image,'var') ~= 1
    error('Anatomical reference image for NLM was not found on the specified path!')
end
if options.RDP_use_anatomical && options.RDP && options.RDPIncludeCorners && options.implementation == 2 && exist(options.RDP_reference_image,'file') ~= 2 && MAP && exist(options.RDP_reference_image,'var') ~= 1
    error('Reference image for RDP was not found on the specified path!')
end
if options.RDP_use_anatomical && options.RDP && ~options.RDPIncludeCorners
    warning('Reference image for RDP is only supported with options.RDPIncludeCorners = true')
end
if options.precondTypeImage(3) && exist(options.referenceImage,'file') ~= 2 && exist(options.referenceImage,'var') ~= 1
    error('Reference image for precondititiong was not found on the specified path!')
end
if options.TV && options.TVtype == 2 && ~options.TV_use_anatomical
    warning('Using TV type = 2, but no anatomical reference set. Using TV type 1 instead.')
    options.TVtype = 1;
end
if options.projector_type > 6 && options.projector_type ~= 11 && options.projector_type ~= 14 && options.projector_type ~= 12 && options.projector_type ~= 13 && ...
        options.projector_type ~= 21 && options.projector_type ~= 22 && options.projector_type ~= 31 && options.projector_type ~= 32 ...
        && options.projector_type ~= 33 && options.projector_type ~= 41 && options.projector_type ~= 51 && options.projector_type ~= 15 && options.projector_type ~= 45 ...
        && options.projector_type ~= 54 && options.projector_type ~= 55 && options.projector_type ~= 44
    error('The selected projector type is not supported!')
end
if options.use_CPU && (options.projector_type == 5 || options.projector_type == 4 || options.projector_type == 14 || options.projector_type == 41 || options.projector_type == 45 ...
        || options.projector_type == 54 || options.projector_type == 51 || options.projector_type == 15)
    error('Selected projector type is not supported with CPU implementation!')
end
if sum(options.precondTypeImage) == 0 && (options.PKMA || options.MRAMLA || options.MBSREM) && ~options.largeDim
    warning('No image-based preconditioner selected with PKMA/MRAMLA/MBSREM. EM preconditioner is highly recommended!')
end
if options.APLS && exist(options.APLS_reference_image,'file') ~= 2 && MAP
    error('APLS selected, but the anatomical reference image was not found on path!')
end
if options.epps <= 0
    warning('Epsilon value is zero or less than zero; must be a positive value. Using the default value (1e-8).')
    options.epps = 1e-8;
end
if numel(options.epps) > 1
    warning('Epsilon has to be a scalar value! Using the default value (1e-8).')
    options.epps = 1e-8;
end
if ~options.CT && ~options.SPECT && options.store_scatter && sum(options.scatter_components) <= 0
    error('Store scatter selected, but no scatter components have been selected!')
end
if options.implementation == 1 && ~options.precompute_lor && (options.n_rays_transaxial > 1 || options.n_rays_axial > 1)
    error('Multiray Siddon is not supported with implementation 1!')
end
if options.implementation == 3 && partitions > 1
    error('Implementation 3 does not support dynamic reconstruction!')
end
if options.precondTypeMeas(2) && (options.subset_type < 8 && options.subset_type ~= 4)
    error('Filtering-based preconditioner only works for subset types 4 and 8-11!')
end
if options.hyperbolic && options.implementation ~= 2
    error('Hyperbolic prior is only available when using implementation 2!')
elseif options.implementation == 2 && options.TV && options.TVtype == 3
    options.TV = false;
    options.hyperbolic = true;
    if options.TV_use_anatomical
        warning('Hyperbolic prior does not support anatomic weighting at the moment when using implementation 2')
        options.TV_use_anatomical = false;
    end
end
if ~options.CT && ~options.SPECT && options.use_machine == 2 && options.use_raw_data
    warning('Sinogram data cannot be used when raw data is set to true, using list-mode data instead')
    options.use_machine = 1;
end
if options.useIndexBasedReconstruction && options.projector_type > 3
    error('Index-based recpnstruction only supports projector types 1-3!')
end
if ~options.CT && ~options.SPECT && options.reconstruct_trues && options.reconstruct_scatter
    warning('Both reconstruct trues and scatter selected, reconstructing only trues.')
    options.reconstruct_scatter = false;
end
if options.implementation == 1 && exist('projector_mex','file') ~= 3
    warning('MEX-file for implementation 1 not found. It is recommended to run install_mex first.')
end
if options.implementation == 4 && exist('projector_mex','file') ~= 3
    error('MEX-file for implementation 4 not found. Run install_mex first.')
end
% if (options.CGLS || options.LSQR || options.FISTA || options.FISTAL1) && options.subsets > 1
if (options.CGLS || options.LSQR) && options.subsets > 1
    warning('CGLS or LSQR do not support subsets! Setting subsets to 1.')
    options.subsets = 1;
end
if options.subsets <= 0
    warning('Subsets set to 0 or less than 0. Setting subsets to 1.')
    options.subsets = 1;
end
if exist('OCTAVE_VERSION','builtin') == 0 && options.use_root && ((exist('GATE_root_matlab','file') ~= 3 && ~verLessThan('matlab', '9.6')) ...
        || (exist('GATE_root_matlab_C','file') ~= 3 && verLessThan('matlab', '9.6'))) && options.use_machine == 0 && ~options.CT
    warning(['ROOT selected, but no MEX-file for ROOT data load found. Run install_mex to build ROOT MEX-file. Ignore this warning if you are ' ...
        'simply loading a mat-file containing measurement data from ROOT files.'])
end
if options.use_root && exist('GATE_root_matlab_oct','file') ~= 3 && options.use_machine == 0 && exist('OCTAVE_VERSION','builtin') == 5 && ~options.CT
    warning(['ROOT selected, but no OCT-file for ROOT data load found. Run install_mex to build ROOT OCT-file. Ignore this warning if you are ' ...
        'simply loading a mat-file containing measurement data from ROOT files.'])
end
if options.use_LMF && exist('gate_lmf_matlab','file') ~= 3 && options.use_machine == 0 && ~options.CT
    error('LMF selected, but no MEX-file for LMF data load found. Run install_mex to build LMF MEX-file.')
end
if options.implementation == 2 && exist('OpenCL_matrixfree','file') ~= 3
    error(['OpenCL reconstruction (implementation 2) selected, but OpenCL MEX-files were not installed. Run install_mex to build OpenCL MEX-files.' sprintf('\n') ...
        'If you already ran install_mex, make sure you have installed both OpenCL and ArrayFire and they can be found on path.'])
end
if options.implementation == 3 && exist('OpenCL_matrixfree_multi_gpu','file') ~= 3
    error(['OpenCL reconstruction (implementation 3) selected, but OpenCL MEX-files were not installed. Run install_mex to build OpenCL MEX-files.' sprintf('\n') ...
        'If you already ran install_mex, make sure you have installed OpenCL and it can be found on path.'])
end
if options.implementation == 5 && exist('OpenCL_matrixfree_multi_gpu','file') ~= 3
    error(['OpenCL reconstruction (implementation 5) selected, but OpenCL MEX-files were not installed. Run install_mex to build OpenCL MEX-files.' sprintf('\n') ...
        'If you already ran install_mex, make sure you have installed OpenCL and it can be found on path.'])
end
if options.implementation == 3 && NMLOS
    warning(['Implementation ' num2str(options.implementation) ' selected with reconstruction algorithms other than OSEM. '...
        'Only OSEM reconstruction can be performed'])
end
if options.implementation == 3 && ~MLOS
    error(['Implementation ' num2str(options.implementation) ' selected, but OSEM algorithm has not been selected.'])
end
if options.implementation ~= 2 && options.RDP && options.RDPIncludeCorners
    error('RDP with include corners is supported only by implementation 2!')
end
if options.implementation == 2 && options.use_CPU && options.RDP && options.RDPIncludeCorners
    error('RDP with include corners is supported only on OpenCL and CUDA!')
end
if PRIOR_summa > 1 && MAP
    error('Only one prior at a time can be used!')
end
if OS_I4_summa > 1

    reko = recNames(0);
    for kk = 1 : length(reko)
        if kk == 1
            dispi = reko{kk};
        else
            dispi = [dispi, ', ' reko{kk}];
        end
        if kk == length(reko) && length(reko) > 1
            dispi = [dispi, ' reconstruction methods selected.'];
        elseif kk == length(reko) && isscalar(reko)
            dispi = [dispi, ' reconstruction method selected.'];
        end
    end
    error(['Implementation ' num2str(options.implementation) ' supports only one algorithm at a time. ' dispi])
end
if options.TOF_bins > 1 && options.TOF_bins_used == 1 && ~options.CT && ~options.SPECT
    disp('Summing TOF bins.')
end
if options.implementation == 1 && ~options.CT
    if options.projector_type == 2 || options.projector_type == 3 || options.projector_type == 31 || options.projector_type == 32 || options.projector_type == 21 || options.projector_type == 12
        error('Orthogonal distance-based/volume-based projector is NOT supported when using implementation 1!')
    end
    if options.TOF_bins_used > 1
        error('TOF is not supported with implementation 1!')
    end
end
if options.TOF_bins_used > 1 && options.implementation == 3
    error('TOF is not supported with implementation 3!')
end
if options.projector_type == 2 && options.implementation == 1 && options.precompute_lor
    warning('Orthogonal distance-based projector is not recommended with implementation 1!')
end
if options.projector_type == 2 && options.CT
    error('Orthogonal distance-based projector is NOT supported when using CT data!')
end
if (options.projector_type == 5 || options.projector_type == 15 || options.projector_type == 51 || options.projector_type == 45 || options.projector_type == 54) && ~options.CT
    error('Projector type 5 is only supported with CT data!')
end
if (options.projector_type == 6) && ~options.SPECT
    error('Projector type 6 is only supported with SPECT data!')
end
if (options.projector_type ~= 6 && options.projector_type ~= 1 && options.projector_type ~= 11 && options.projector_type ~= 2 && options.projector_type ~= 22 && options.projector_type ~= 7) && options.SPECT
    disp(options.projector_type)
    error('SPECT only supports projector types 1, 2 and 6!')
end
if (options.projector_type == 6)
    if options.Nx(1) ~= options.nRowsD
        error('options.Nx has to be the same as options.nRowsD when using projector type 6')
    end
    if options.Ny(1) ~= options.nRowsD
        error('options.Ny has to be the same as options.nRowsD when using projector type 6')
    end
    if options.Nz(1) ~= options.nColsD
        error('options.Nz has to be the same as options.nColsD when using projector type 6')
    end
    if options.subsets > 1 && options.subset_type < 8
        error('Subset types 0-7 are not supported with projector type 6!')
    end
end
if options.FDK && (options.Niter > 1 || options.subsets > 1)
    if options.largeDim
        options.Niter = 1;
    else
        warning('When using FDK/FBP, the number of iterations and subsets must be set as 1. Setting both to 1.')
        options.subsets = 1;
        options.Niter = 1;
    end
end
if options.use_CUDA && options.use_CPU && options.implementation == 2
    error('Both CUDA and CPU selected! Select only one!')
end
if options.TOF_bins_used > 1 && (options.projector_type ~= 1 && options.projector_type ~= 11 && options.projector_type ~= 3 && options.projector_type ~= 33 && options.projector_type ~= 31 ...
        && options.projector_type ~= 13 && options.projector_type ~= 4 && options.projector_type ~= 41 && options.projector_type ~= 14) && ~options.CT && ~options.SPECT
    error('TOF is currently only supported with improved Siddon (projector_type = 1), interpolation-based projector (projector_type = 4) and volume of intersection (projector_type = 3)')
end
if options.TOF_bins_used > 1 && options.TOF_width <= 0 && ~options.CT && ~options.SPECT
    error('TOF width (options.TOF_width) must be greater than zero.')
end
if options.TOF_bins_used > 1 && options.TOF_FWHM == 0 && ~options.CT && ~options.SPECT
    error('TOF enabled, but the TOF FWHM (options.TOF_FWHM) is zero. FWHM must be nonzero.')
end
if options.TOF_bins_used > 1 && options.use_raw_data && ~options.CT && ~options.SPECT
    error('TOF data is only available with sinogram data. Disable raw data (options.use_raw_data = false).')
end
if options.corrections_during_reconstruction && (options.scatter_correction || options.randoms_correction) && (options.PDHG || options.PDHGL1 || options.FISTA || options.LSQR || options.CGLS || options.FISTAL1)
    error('Randoms/scatter correction cannot be applied during the reconstruction with the selected algorithm!')
end
% Print various options that were selected if verbosity has been enabled
if options.verbose > 0
    if options.use_ASCII && options.use_machine == 0
        dispi = 'Using ASCII data';
    elseif options.use_LMF && options.use_machine == 0
        dispi = 'Using LMF data';
    elseif options.use_root && options.use_machine == 0
        dispi = 'Using ROOT data';
    elseif options.use_binary && options.use_machine == 0
        dispi = 'Using BINARY data';
    elseif options.use_machine == 1
        dispi = 'Using data obtained from list-mode file';
    elseif options.use_machine == 2
        dispi = 'Using scanner created sinogram data';
    elseif options.use_machine == 3
        dispi = 'Using 32-bit list-mode data';
    else
        dispi = [];
    end
    if options.TOF_bins_used > 1 && options.TOF_FWHM > 0
        dispi = [dispi ' with TOF (' num2str(options.TOF_bins_used) ' bins).'];
    else
        dispi = [dispi '.'];
    end
    if ~strcmp(dispi,'.')
        disp(dispi);
    end
    if options.only_sinos
        disp('Loading only data.')
    end
    if ~(options.compute_normalization || options.only_sinos)
        dispaus = (['Using implementation ' num2str(options.implementation)]);
        if options.implementation == 1 || options.implementation == 4
            dispaus = [dispaus, '.'];
        elseif options.implementation == 2
            if ~options.use_CUDA && ~options.use_CPU
                inffo = ArrayFire_OpenCL_device_info();
                k2 = strfind(inffo, 'MB');
                if options.use_device == 0
                    k1 = strfind(inffo, '[0]');
                    dispaus = [dispaus, ' with ', inffo(k1 + 4:k2(1)+1)];
                else
                    k1 = strfind(inffo, ['-' num2str(options.use_device) '-']);
                    dispaus = [dispaus, ' with ', inffo(k1 + 4:k2(options.use_device + 1)+1)];
                end
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
        varNonMAP = [recNames(6)];
        for kk = 1 : numel(varNonMAP)
            if options.OSEM && options.subsets == 1
                reko = [reko;{'MLEM'}];
                break
            elseif options.(varNonMAP{kk})
                ch = strrep(varNonMAP{kk},'_','-');
                reko = [reko;ch];
                break
            end
        end

        if PRIOR_summa > 1
            error('Only one prior can be used at a time!')
        end
        varMAP = recNames(2);
        for kk = 1 : numel(varMAP)
            if options.OSL_OSEM && PRIOR && options.subsets == 1
                reko = [reko;{'MLEM-OSL'}];
            elseif options.(varMAP{kk}) && strcmp(options.(varMAP{kk}),'OSL_COSEM')
                if options.OSL_COSEM == 1 && PRIOR
                    reko = [reko;{'ACOSEM-OSL'}];
                elseif options.OSL_COSEM == 1 && ~PRIOR
                    error('ACOSEM-OSL selected, but no prior has been selected.')
                elseif options.OSL_COSEM == 2 && PRIOR
                    reko = [reko;{'COSEM-OSL'}];
                elseif options.OSL_COSEM == 2 && ~PRIOR
                    error('COSEM-OSL selected, but no prior has been selected.')
                elseif options.OSL_COSEM > 2 || options.OSL_COSEM < 0
                    error('Unsupported COSEM-OSL method selected!')
                end
                if PRIOR_summa > 1
                    error('Only one prior can be used at a time with COSEM-OSL')
                end
            elseif options.(varMAP{kk}) && PRIOR
                ch = strrep(varMAP{kk},'_','-');
                reko = [reko;ch];
            elseif options.(varMAP{kk}) && ~PRIOR
                ch = strrep(varMAP{kk},'_','-');
                warning([ch ' selected, but no prior has been selected.'])
                reko = [reko;ch];
                options.beta = 0;
            end
        end
        for kk = 1 : length(reko)
            if kk == 1
                dispi = reko{kk};
            else
                dispi = [dispi, ', ' reko{kk}];
            end
            dispi = [dispi, ' reconstruction method selected.'];
        end
        if ~isempty(reko)
            disp(dispi)
        else
            error('No reconstruction method selected!')
        end
        priori = {};

        varPrior = recNames(1);
        varPriorName = recNames(10);
        for kk = 1 : numel(varPrior)
            if options.(varPrior{kk}) && MAP
                if strcmp(options.(varPrior{kk}),'TV')
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
                elseif strcmp(options.(varPrior{kk}),'weighted_mean')
                    if options.mean_type == 1 || options.mean_type == 4
                        priori = [priori;{'Weighted (arithmetic) mean'}];
                    elseif options.mean_type == 2 || options.mean_type == 5
                        priori = [priori;{'Weighted (harmonic) mean'}];
                    elseif options.mean_type == 3 || options.mean_type == 6
                        priori = [priori;{'Weighted (geometric) mean'}];
                    else
                        error('Unsupported mean type selected.')
                    end
                elseif strcmp(options.(varPrior{kk}),'AD')
                    if options.FluxType > 2 || options.FluxType < 1
                        error('FluxType has to be either 1 or 2.')
                    end
                    if (options.DiffusionType > 2 || options.DiffusionType < 1) && options.implementation == 2
                        error('DiffusionType has to be either 1 or 2.')
                    end
                    priori = [priori;{'AD-MRP'}];
                else
                    ch = strrep(varPriorName{kk},'_','-');
                    priori = [priori;ch];
                end
            elseif options.(varPrior{kk}) && ~MAP
                ch = strrep(varPriorName{kk},'_','-');
                warning([ch ' selected, but no MAP algorithm has been selected. No MAP reconstruction will be performed.'])
                options.(varPrior{kk}) = false;
            end
        end
        for kk = 1 : length(priori)
            if kk == 1
                dispi2 = priori{kk};
            else
                dispi2 = [dispi2, ', ' priori{kk}];
            end
            if kk == length(priori) && length(priori) > 1
                dispi2 = [dispi2, ' priors selected.'];
            elseif kk == length(priori) && isscalar(priori)
                if ~options.NLM
                    dispi2 = [dispi2, ' prior selected.'];
                else
                    if options.NLTV
                        dispi2 = [dispi2, ' prior selected with NLTV.'];
                    elseif options.NLRD 
                        dispi2 = [dispi2, ' prior selected with NLRD.'];
                    elseif options.NLLange
                        dispi2 = [dispi2, ' prior selected with NL Lange.'];
                    elseif options.NLGGMRF
                        dispi2 = [dispi2, ' prior selected with NLGGMRF.'];
                    elseif options.NLM_MRP
                        dispi2 = [dispi2, ' prior selected with filtering mode.'];
                    else
                        dispi2 = [dispi2, ' prior selected.'];
                    end
                    if options.NLAdaptive
                        dispi2 = [dispi2, ' Using adaptive weighting.'];
                    end
                end
            end
        end
        if ~isempty(priori)
            disp(dispi2)
            if options.NLM || options.MRP || options.quad || options.Huber || (options.RDP && options.RDPIncludeCorners) || options.GGMRF
                disp(['Using neighborhood/search window size of ' num2str(options.Ndx * 2 + 1) 'x' num2str(options.Ndy * 2 + 1) 'x' num2str(options.Ndz * 2 + 1) '.'])
            end
            if options.NLM
                disp(['Using patch size of ' num2str(options.Nlx * 2 + 1) 'x' num2str(options.Nly * 2 + 1) 'x' num2str(options.Nlz * 2 + 1) '.'])
            end
        end
        if ~OS && ~MAP && ~options.MLEM && ~(options.only_sinos || options.compute_normalization)
            error('No reconstruction algorithm selected.')
        end
        if options.projector_type == 1 && ~options.precompute_lor
            if options.implementation == 1
                disp('Improved Siddon''s algorithm selected with 1 ray.')
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
                if options.SPECT
                    disp(['Improved Siddon''s algorithm selected with ' num2str(options.n_rays_transaxial) ' ' ray '.'])
                else
                    disp(['Improved Siddon''s algorithm selected with ' num2str(options.n_rays_transaxial) ' transaxial ' ray ' and ' ...
                    num2str(options.n_rays_axial) ' axial ' aray '.'])
                end
            end
        elseif options.projector_type == 1 || options.projector_type == 11
            disp('Improved Siddon''s algorithm selected with 1 ray.')
        elseif options.projector_type == 2 || options.projector_type == 22
            dispi = 'Orthogonal distance-based ray tracer selected';
            if options.tube_width_z > 0 || options.implementation == 4
                dispi = [dispi, ' in 3D mode.'];
            elseif options.tube_width_z == 0 && options.implementation == 4
                warning('2.5D mode not available with implementation 4! Switching to 3D mode!')
                options.tube_width_z = options.tube_width_xy;
                dispi = [dispi, ' in 3D mode.'];
            else
                dispi = [dispi, ' in 2.5D mode.'];
            end
            disp(dispi)
        elseif options.projector_type == 21
            disp('Improved Siddon''s algorithm selected for forward projection, orthogonal for backprojection.')
        elseif options.projector_type == 12
            disp('Orthogonal selected for forward projection, improved Siddon''s algorithm for backprojection.')
        elseif options.projector_type == 3 || options.projector_type == 33
            disp('Volume of intersection based ray tracer selected.');
        elseif options.projector_type == 31
            disp('Improved Siddon''s algorithm selected for forward projection, Volume of intersection based ray tracer for backprojection.')
        elseif options.projector_type == 13
            disp('Volume of intersection based ray tracer selected for forward projection, improved Siddon''s algorithm for backprojection.')
        elseif options.projector_type == 4
            disp('Interpolation-based projector selected.')
        elseif options.projector_type == 5
            disp('Branchless distance-driven based projector selected.')
        elseif options.projector_type == 41
            disp('Interpolation-based projector selected for forward projection, improved Siddon for backprojection.')
        elseif options.projector_type == 14
            disp('Improved Siddon projector selected for forward projection, interpolation-based projector for backprojection.')
        elseif options.projector_type == 15
            disp('Improved Siddon projector selected for forward projection, branchless distance-driven projector for backprojection.')
        elseif options.projector_type == 45
            disp('Interpolation-based projector selected for forward projection, branchless distance-driven projector for backprojection.')
        elseif options.projector_type == 54
            disp('Branchless distance-driven projector selected for forward projection, interpolation-based projector for backprojection.')
        elseif options.projector_type == 51
            disp('Branchless distance-driven projector selected for forward projection, improved Siddon for backprojection.')
        elseif options.projector_type == 6
            disp('Rotation-based projector selected (SPECT).')
        end
        if options.use_psf
            if options.deblurring
                disp('PSF ON with deblurring phase.');
            else
                disp('PSF ON.');
            end

        end
        if options.attenuation_correction && ~options.CT
            disp('Attenuation correction ON.')
        end
        if options.randoms_correction && ~options.CT
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
        if options.fill_sinogram_gaps && ~options.CT
            disp('Sinogram gap filling ON.')
        end
        if any(options.precondTypeImage)
            if options.precondTypeImage(1)
                disp('Using image-based preconditioning with diagonal preconditioner.')
            elseif options.precondTypeImage(2)
                disp('Using image-based preconditioning with EM preconditioner.')
            elseif options.precondTypeImage(3)
                disp('Using image-based preconditioning with IEM preconditioner.')
            end
            if options.precondTypeImage(4)
                disp('Using image-based preconditioning with momentum preconditioner.')
            end
            if options.precondTypeImage(5)
                disp('Using image-based preconditioning with normalized gradient preconditioner.')
            end
            if options.precondTypeImage(6)
                dispP = ('Using image-based preconditioning with filtering preconditioner');
                if strcmp(options.filterWindow,'hamming')
                    dispP = [dispP, ' with Hamming window.'];
                elseif strcmp(options.filterWindow,'hann')
                    dispP = [dispP, ' with Hann window.'];
                elseif strcmp(options.filterWindow,'blackman')
                    dispP = [dispP, ' with Blackman window.'];
                elseif strcmp(options.filterWindow,'nuttal')
                    dispP = [dispP, ' with Nuttal window.'];
                elseif strcmp(options.filterWindow,'gaussian')
                    dispP = [dispP, ' with Gaussian window.'];
                elseif strcmp(options.filterWindow,'shepp-logan')
                    dispP = [dispP, ' with Shepp-Logan window.'];
                elseif strcmp(options.filterWindow,'cosine')
                    dispP = [dispP, ' with cosine window.'];
                elseif strcmp(options.filterWindow,'parzen')
                    dispP = [dispP, ' with Parzen window.'];
                else
                    dispP = [dispP, ' (no windowing).'];
                end
                disp(dispP)
            end
        end
        if any(options.precondTypeMeas)
            if options.precondTypeMeas(1)
                disp('Using measurement-based preconditioning with diagonal preconditioner.')
            end
            if options.precondTypeMeas(2)
                dispP = ('Using measurement-based preconditioning with filtering preconditioner');
                if strcmp(options.filterWindow,'hamming')
                    dispP = [dispP, ' with Hamming window.'];
                elseif strcmp(options.filterWindow,'hann')
                    dispP = [dispP, ' with Hann window.'];
                elseif strcmp(options.filterWindow,'blackman')
                    dispP = [dispP, ' with Blackman window.'];
                elseif strcmp(options.filterWindow,'nuttal')
                    dispP = [dispP, ' with Nuttal window.'];
                elseif strcmp(options.filterWindow,'gaussian')
                    dispP = [dispP, ' with Gaussian window.'];
                elseif strcmp(options.filterWindow,'shepp-logan')
                    dispP = [dispP, ' with Shepp-Logan window.'];
                elseif strcmp(options.filterWindow,'cosine')
                    dispP = [dispP, ' with cosine window.'];
                elseif strcmp(options.filterWindow,'parzen')
                    dispP = [dispP, ' with Parzen window.'];
                else
                    dispP = [dispP, ' (no windowing).'];
                end
                disp(dispP)
            end
        end
        if options.oOffsetZ ~= 0 || options.oOffsetX ~= 0 || options.oOffsetY ~= 0
            disp(['Object offset is [' num2str(options.oOffsetX) ', ' num2str(options.oOffsetY) ', ' num2str(options.oOffsetZ) '] (XYZ).'])
        end
    end
    if options.normalization_correction && ~(options.compute_normalization || options.only_sinos) && ~options.CT
        disp('Normalization correction ON.')
    elseif options.normalization_correction && options.compute_normalization && ~options.CT
        warning('Normalization correction cannot be applied when computing normalization coefficients. Disabling normalization correction.')
        options.normalization_correction = false;
    elseif options.compute_normalization && ~options.CT
        disp('Computing normalization coefficients.')
    end
    if options.compute_normalization && sum(options.normalization_options) == 0 && ~options.CT
        error('Normalization computation selected, but no normalization components selected.')
    end
    if options.useMultiResolutionVolumes
        if options.transaxialEFOV
            Nx = options.NxOrig + round((options.Nx - options.NxOrig) * options.multiResolutionScale);
            Ny = options.NyOrig + round((options.Ny - options.NyOrig) * options.multiResolutionScale);
        else
            Nx = options.Nx;
            Ny = options.Ny;
        end
        if options.axialEFOV
            Nz = options.NzOrig + round((options.Nz - options.NzOrig) * options.multiResolutionScale);
        else
            Nz = options.Nz;
        end
    else
        Nx = options.Nx;
        Ny = options.Ny;
        Nz = options.Nz;
    end
    if ~(options.compute_normalization || options.only_sinos)
        if options.corrections_during_reconstruction && (options.normalization_correction || options.randoms_correction || options.scatter_correction)
            disp('Corrections applied during reconstruction (ordinary Poisson).')
        elseif ~options.corrections_during_reconstruction && (options.normalization_correction || options.randoms_correction || options.scatter_correction)
            disp('Corrections applied to the measurement data.')
        end
        if options.arc_correction && ~options.use_raw_data && ~options.CT
            disp('Arc correction ON.')
        end
        if options.use_raw_data && ~options.CT
            if partitions == 1 || isempty(partitions)
                dispi = 'Using STATIC raw data';
            else
                dispi = 'Using DYNAMIC raw data';
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
            if partitions > 1
                if options.sampling_raw > 1
                    dispi = strcat(dispi, [' and ' num2str(partitions) ' time steps']);
                else
                    dispi = strcat(dispi, [' with ' num2str(partitions) ' time steps']);
                end
            end
            dispi = strcat(dispi, '.');
            disp(dispi)
        else
            if partitions == 1
                if options.CT
                    dispi = 'Using STATIC projection data';
                else
                    dispi = 'Using STATIC sinogram data';
                end
            else
                if options.CT
                    dispi = 'Using DYNAMIC projection data';
                else
                    dispi = 'Using DYNAMIC sinogram data';
                end
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
            if numel(options.partitions) > 1
                if options.sampling > 1
                    dispi = strcat(dispi, [' and ' num2str(numel(options.partitions)) ' time steps']);
                else
                    dispi = strcat(dispi, [' with ' num2str(numel(options.partitions)) ' time steps']);
                end
            elseif numel(options.partitions) == 1 && options.partitions > 1
                if options.sampling > 1
                    dispi = strcat(dispi, [' and ' num2str(options.partitions) ' time steps']);
                else
                    dispi = strcat(dispi, [' with ' num2str(options.partitions) ' time steps']);
                end
            end
            dispi = strcat(dispi, '.');
            disp(dispi)
        end
        if options.subsets > 1
            if options.subset_type == 1
                disp(['Every ' num2str(options.subsets) 'th column measurement is taken per subset.'])
            elseif options.subset_type == 2
                disp(['Every ' num2str(options.subsets) 'th row measurement is taken per subset.'])
            elseif options.subset_type == 3
                disp('Using random subset sampling.')
            elseif options.subset_type == 4
                disp(['Every ' num2str(options.subsets) 'th sinogram column is taken per subset.'])
            elseif options.subset_type == 5
                disp(['Every ' num2str(options.subsets) 'th sinogram row is taken per subset.'])
            elseif options.subset_type == 6
                disp(['Using angle-based subset sampling with ' num2str(options.n_angles) ' angles combined per subset.'])
            elseif options.subset_type == 7
                disp('Using golden angle-based subset sampling.')
            elseif options.subset_type == 8
                disp(['Using every ' num2str(options.subsets) 'th sinogram/projection image.'])
            elseif options.subset_type == 9
                disp('Using sinograms/projection images in random order.')
            elseif options.subset_type == 10
                disp('Using golde angle sampling with sinograms/projections.')
            elseif options.subset_type == 11
                disp('Using prime factor ordering of projections/sinograms into subsets.')
            elseif options.subset_type == 0
                disp(['Dividing data into ' num2str(options.subsets) ' segments.'])
            end
        end
        disp(['Using an image (matrix) size of ' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) ' with ' num2str(options.Niter) ...
            ' iterations and ' num2str(options.subsets) ' subsets.'])
    elseif options.CT
        disp(['Using an image (matrix) size of ' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) ' with ' num2str(options.Niter) ...
            ' iterations and ' num2str(options.subsets) ' subsets.'])
    end
end