function options = OMEGA_error_check(options)
%% Error checking file
% This function is used to check that all the input values are allowed.
% The selected settings are printed by printReconstructionSettings.

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


% SPECT collimator ray shifts are stored once per detector head and detector
% element.
if options.SPECT
    if options.normalization_correction && isfield(options, 'normalization') && numel(options.normalization) > 1
        options.normZ = numel(options.normalization) / (options.nRowsD * options.nColsD);
        if mod(options.normZ, 1) ~= 0 || ~(options.normZ == options.nProjections || options.normZ == options.nHeads)
            error(['Normalization array has an invalid size. Must contain one image per projection [' ...
                num2str(options.nRowsD) ' ' num2str(options.nColsD) ' ' num2str(options.nProjections) ...
                '] or one image per detector head [' num2str(options.nRowsD) ' ' num2str(options.nColsD) ' ' ...
                num2str(options.nHeads) '].'])
        end
    end
    if ~isfield(options, 'DetectorVector') || isempty(options.DetectorVector)
        options.DetectorVector = zeros(options.nProjections, 1, 'uint32');
    elseif numel(options.DetectorVector) ~= options.nProjections
        error(['DetectorVector must contain one detector-head index for each projection (' ...
            num2str(options.nProjections) ' values required).'])
    else
        options.DetectorVector = uint32(options.DetectorVector(:));
    end
    if any(options.DetectorVector >= options.nHeads)
        error('DetectorVector contains an index outside the available detector heads.')
    end
    if options.useMaskFP && options.maskFPZ == options.nHeads && ...
            numel(options.maskFP) ~= options.nRowsD * options.nColsD * options.nHeads
        error('Detector-indexed forward mask must contain one image for each detector head.')
    end

    % Implementation 4 (projector_mex.cpp) and implementation 2 with use_CPU (CPU_matrixfree, built
    % from OpenCL_matrixfree.cpp with -DCPU and no GPU backend macro) both route through
    % ProjectorClassCPU.h / projector_functions.h (see functions.hpp's #elif defined(CPU) branch),
    % a separate plain C++ projector with no DetectorVector/per-head indexing at all. Implementation
    % 1 already errors out for all SPECT data above.
    if (options.implementation == 4 || (options.implementation == 2 && options.use_CPU)) && ...
            ((options.normalization_correction && options.normZ == options.nHeads && options.normZ ~= options.nProjections) || ...
             (options.useMaskFP && options.maskFPZ > 1 && options.maskFPZ == options.nHeads && options.maskFPZ ~= options.nProjections))
        error(['Detector-head indexed normalization/forward projection mask (nRowsD x nColsD x nHeads) is not ' ...
            'supported by the CPU implementation. Use a GPU implementation or supply one image per projection.'])
    end

    if ismember(options.projector_type, [1, 11, 12, 2, 21, 22, 3, 13, 23, 31, 32, 33])
        if isfield(options, 'nHeads') && ~isempty(options.nHeads)
            nHeads = double(options.nHeads);
        else
            nHeads = 1;
        end
        compactRayShiftSize = 2 * double(options.n_rays_transaxial) * double(options.n_rays_axial) * double(options.nRowsD) * ...
            double(options.nColsD) * nHeads;

        detectorSize = numel(options.rayShiftsDetector);
        sourceSize = numel(options.rayShiftsSource);
        if detectorSize > 0 && detectorSize ~= compactRayShiftSize
            error(['rayShiftsDetector has an invalid size. Expected compact size ' ...
                num2str(compactRayShiftSize) ', got ' num2str(detectorSize) '.'])
        end
        if sourceSize > 0 && sourceSize ~= compactRayShiftSize
            error(['rayShiftsSource has an invalid size. Expected compact size ' ...
                num2str(compactRayShiftSize) ', got ' num2str(sourceSize) '.'])
        end
        if detectorSize > 0 && sourceSize > 0 && detectorSize ~= sourceSize
            error('rayShiftsDetector and rayShiftsSource must have the same compact size.')
        end

        if ~(options.ellipsePower == 2 || isinf(options.ellipsePower))
            error(['ellipsePower must be 2 (elliptic cylinder) or Inf (box); other superellipse ' ...
                'powers are not supported.'])
        end
        if options.ellipseRadiusX <= 0 || options.ellipseRadiusY <= 0 || options.ellipseRadiusZ <= 0
            error('ellipseRadiusX, ellipseRadiusY and ellipseRadiusZ must all be greater than zero.')
        end
    end
end

if ismac && options.use_64bit_atomics
    warning(['The current Metal compiler does not expose native atomic_ulong ' ...
        'addition. Disabling 64-bit atomics; use 32-bit atomics for older Metal hardware.'])
    options.use_64bit_atomics = false;
end

% Determine whether various different reconstruction modes are used (e.g.
% MAP reconstruction and any prior)
MAP = checkAlgorithmsPriors(options, 2);
PRIOR_summa = checkAlgorithmsPriors(options, 1, 0);
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
if partitions > 1 && options.implementation ~= 2 % TODO dynamic reconstruction for implementation 5 (macOS)
    error('Dynamic reconstruction is supported only with implementation 2.')
end
if options.SPECT && mod(sqrt(options.nRays), 1) ~= 0
    error('With SPECT, options.nRays has to be a square')
end
if options.SPECT && ismember(options.projector_type, [2, 12, 21, 22, 26, 62]) && options.n_rays_transaxial * options.n_rays_axial > 1
    error('Multiple rays are not supported with the orthogonal distance-based projector (projector type 2) in SPECT. Use n_rays_transaxial = n_rays_axial = 1.')
end
if options.only_sinos && options.only_reconstructions
    error('options.only_sinos and options.only_reconstructions cannot be both set to true')
end
% Check for various illegal values
if ~options.CT && ~options.SPECT && (options.FOVa_x >= options.diameter || options.FOVa_y >= options.diameter)
    warning(['Transaxial FOV is larger than the scanner diameter (' num2str(options.diameter) ')!'])
end
if ~options.CT && ~options.SPECT && (options.axial_fov) < (options.rings(1) * options.cr_pz - options.cr_pz) && options.rings > 0
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
options = checkSinogramOptions(options);
listModeCandidate = isListModeCandidate(options);
if options.SPECT && options.implementation == 1
    error('Implementation 1 is not supported with SPECT data.')
end
if options.SPECT && options.implementation == 3 && options.projector_type == 6
    error('Implementation 3 is not supported with projector type 6.')
end
if options.SPECT && options.implementation == 5 && options.projector_type == 6
    error('Implementation 5 is not supported with projector type 6.')
end
if partitions < 1
    warning('Number of time steps is less than one. Using one time step.')
    options.partitions = 1;
    partitions = 1;
end
if partitions == 1 && (options.temporal_smoothness || options.temporalTV)
    warning('Temporal regularization selected, but using static data. No temporal regularization will be done.')
    options.temporal_smoothness = false;
    options.temporalTV = false;
elseif partitions > 1 && options.temporal_smoothness && options.temporalTV
    error('Both temporal smoothness prior and temporal TV selected! Select only one!')
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
if options.implementation == 2 && options.use_CPU && options.useMultiResolutionVolumes
    error('Multi-resolution reconstruction is not supported with the CPU implementation!')
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
if ~options.largeDim && isempty(options.x0)
    if ~options.CT
        warning('Initial value is an empty array, using the default values (1)')
        options.x0 = ones(options.Nx, options.Ny, options.Nz);
    else
        warning('Initial value is an empty array, using the default values (1e-4)')
        options.x0 = ones(options.Nx, options.Ny, options.Nz) * 1e-4;
    end
end
if ~options.largeDim && size(options.x0,1)*size(options.x0,2)*size(options.x0,3) < options.Nx*options.Ny*options.Nz && ~options.useEFOV
    error(['Initial value has a matrix size smaller (' num2str(size(options.x0,1)*size(options.x0,2)*size(options.x0,3)) ') than the actual image size ('...
        num2str(options.Nx*options.Ny*options.Nz) ')!'])
elseif ~options.largeDim && options.useEFOV && (size(options.x0,1)*size(options.x0,2)*size(options.x0,3) < options.Nx*options.Ny*options.Nz && size(options.x0,1)*size(options.x0,2)*size(options.x0,3) < options.NxOrig*options.NyOrig*options.NzOrig)
    error(['Initial value has a matrix size smaller (' num2str(size(options.x0,1)*size(options.x0,2)*size(options.x0,3)) ') than the actual image size ('...
        num2str(options.NxOrig*options.NyOrig*options.NzOrig) ')!'])
end
if ~options.largeDim && size(options.x0,1)*size(options.x0,2)*size(options.x0,3) > options.Nx*options.Ny*options.Nz && ~options.useEFOV
    warning(['Initial value has a matrix size larger (' num2str(size(options.x0,1)*size(options.x0,2)*size(options.x0,3)) ') than the actual image size (' ...
        num2str(options.Nx*options.Ny*options.Nz) ')! Attempting automatic resize.'])
    apuVar = numel(options.x0) / (options.Nx * options.Ny * options.Nz);
    apuVar = (apuVar^(1/3));
    options.x0 = imresize3(reshape(options.x0, round(options.Nx * apuVar), round(options.Ny * apuVar), ceil(options.Nz * apuVar)), [options.Nx, options.Ny, options.Nz]);
end
if options.TV_use_anatomical && options.TV && MAP && (ischar(options.TV_reference_image) || isstring(options.TV_reference_image)) && exist(options.TV_reference_image,'file') ~= 2
    error('Anatomical reference image for TV was not found on the specified path!')
end
if options.NLM_use_anatomical && options.NLM && MAP && (ischar(options.NLM_reference_image) || isstring(options.NLM_reference_image)) && exist(options.NLM_reference_image,'file') ~= 2
    error('Anatomical reference image for NLM was not found on the specified path!')
end
if options.RDP_use_anatomical && options.RDP && options.RDPIncludeCorners && options.implementation == 2 && MAP && (ischar(options.RDP_reference_image) || isstring(options.RDP_reference_image)) && exist(options.RDP_reference_image,'file') ~= 2
    error('Reference image for RDP was not found on the specified path!')
end
if options.RDP_use_anatomical && options.RDP && ~options.RDPIncludeCorners
    error('Reference image for RDP is only supported with options.RDPIncludeCorners = true')
end
if options.precondTypeImage(3) && (ischar(options.referenceImage) || isstring(options.referenceImage)) && exist(options.referenceImage,'file') ~= 2
    error('Reference image for precondititiong was not found on the specified path!')
end
if options.TV && options.TVtype == 2 && ~options.TV_use_anatomical
    warning('Using TV type = 2, but no anatomical reference set. Using TV type 1 instead.')
    options.TVtype = 1;
end
if options.implementation == 1 && ((isfield(options,'maskBP') && ~isscalar(options.maskBP)) || (isfield(options,'maskFP') && ~isscalar(options.maskFP)))
    warning('Mask images are not supported with implementation 1!')
end
if options.projector_type > 6 && ~ismember(options.projector_type, ...
        [11 12 13 14 15 16 21 22 23 24 25 26 31 32 33 34 35 41 42 43 44 45 51 54 55 61 62 66])
    error('The selected projector type is not supported!')
end
if options.use_CPU && ~ismember(options.projector_type, [11 1]) && options.implementation == 2
    error('Selected projector type is not supported with CPU implementation!')
end
if options.use_device < 0
    error('Device number (options.use_device) has to be zero or positive!')
end
if options.implementation == 4 && ~ismember(options.projector_type, [1 11 6 16 61 66])
    error('Selected projector type is not supported with CPU implementation!')
end
if (isfield(options,'ProxRDP') && options.ProxRDP) || (isfield(options,'ProxNLM') && options.ProxNLM)
    error('ProxRDP and ProxNLM have been removed. Use RDP or NLM instead.')
end
if sum(options.precondTypeImage) == 0 && (options.PKMA || options.MRAMLA || options.MBSREM)
    warning('No image-based preconditioner selected with PKMA/MRAMLA/MBSREM. EM preconditioner is highly recommended!')
end
if options.APLS && MAP && (ischar(options.APLS_reference_image) || isstring(options.APLS_reference_image)) && ...
        ~isempty(char(options.APLS_reference_image)) && exist(char(options.APLS_reference_image),'file') ~= 2
    error('APLS selected, but the anatomical reference image was not found on path!')
end
if options.epps <= 0
    warning('Epsilon value is zero or less than zero; must be a positive value. Using the default value (1e-6).')
    options.epps = 1e-6;
end
if numel(options.epps) > 1
    warning('Epsilon has to be a scalar value! Using the default value (1e-6).')
    options.epps = 1e-6;
end
if options.implementation == 1 && ~options.precompute_lor && (options.n_rays_transaxial > 1 || options.n_rays_axial > 1)
    error('Multiray Siddon is not supported with implementation 1!')
end
if options.implementation == 3 && partitions > 1
    error('Implementation 3 does not support dynamic reconstruction!')
end
if options.precondTypeMeas(2) && (options.subset_type < 8 && options.subset_type ~= 4 && options.subset_type ~= 0)
    error('Filtering-based preconditioner only works for subset types 0, 4 and 8-11!')
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
% Forward and backward projector type of the (possibly combined) projector type
ptDigits = [floor(options.projector_type / 10), mod(options.projector_type, 10)];
if options.useIndexBasedReconstruction && any(ismember(ptDigits, [5 6]))
    error('Index-based reconstruction only supports projector types 1-4!')
end
% List-mode data (custom coordinates or index-based) can only be divided into
% subsets with subset types 0 (sequential), 1 (every n-th event) and 3
% (random). Other types (e.g. the default 8) are sinogram/projection-based
% and are rejected here, before any data is prepared, instead of being
% silently replaced. listModeCandidate only detects the custom coordinates
% case when options.x and options.SinM have already been set; projectorClass
% repeats this check once all the data are available.
if ~options.CT && ~options.SPECT && (options.useIndexBasedReconstruction || options.listmode || listModeCandidate) ...
        && options.subsets > 1 && ~ismember(options.subset_type, [0, 1, 3])
    error(['Subset type ' num2str(options.subset_type) ' is not supported with list-mode and index-based data! Only subset types 0 (sequential), 1 (every n-th event) and 3 (random) can be used when subsets > 1. Set options.subset_type to 0, 1 or 3 (or use options.subsets = 1).'])
end
% The sinogram plane index of subset types 3, 6 and 7 (z_index) is stored as 16-bit
if ~options.CT && ~options.SPECT && ~options.listmode && ~listModeCandidate && ~options.useIndexBasedReconstruction && ~options.use_raw_data ...
        && options.subsets > 1 && ismember(options.subset_type, [3, 6, 7]) && options.NSinos > 65536
    error('Subset types 3, 6 and 7 support at most 65536 sinogram planes (the plane index is stored as 16-bit); use another subset type, fewer planes (larger span) or a single subset.')
end
% Subset type 6 groups n_angles degrees per subset (180/n_angles subsets for 2D and 360/n_angles for 3D data)
if ~options.CT && ~options.SPECT && ~options.listmode && ~listModeCandidate && ~options.useIndexBasedReconstruction && ~options.use_raw_data ...
        && options.subsets > 1 && options.subset_type == 6 && (~isfield(options, 'n_angles') || ~isscalar(options.n_angles) ...
        || options.n_angles < 1 || mod(180, options.n_angles) ~= 0)
    error('Subset type 6 requires options.n_angles to be a positive number that divides 180 evenly!')
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
if (options.CGLS || options.LSQR) && options.subsets > 1 && ~options.largeDim
    warning('CGLS and LSQR do not support subsets! Setting subsets to 1.')
    options.subsets = 1;
end
if options.subsets <= 0
    warning('Subsets set to 0 or less than 0. Setting subsets to 1.')
    options.subsets = 1;
end
if options.implementation == 2 && exist('OpenCL_matrixfree','file') ~= 3
    error(['OpenCL reconstruction (implementation 2) selected, but OpenCL MEX-files were not installed. Run install_mex to build OpenCL MEX-files.' sprintf('\n') ...
        'If you already ran install_mex, make sure you have installed both OpenCL and ArrayFire and they can be found on path.'])
end
if options.implementation == 3 && exist('OpenCL_matrixfree_multi_gpu','file') ~= 3
    error(['OpenCL reconstruction (implementation 3) selected, but OpenCL MEX-files were not installed. Run install_mex to build OpenCL MEX-files.' sprintf('\n') ...
        'If you already ran install_mex, make sure you have installed OpenCL and it can be found on path.'])
end
if options.implementation == 5 && exist('OpenCL_matrixfree_multi_gpu','file') ~= 3 && exist('CUDA_matrixfree_multi_gpu','file') ~= 3
    error(['Implementation 5 selected, but neither the OpenCL nor the CUDA MEX-files were installed. Run install_mex to build them.' sprintf('\n') ...
        'If you already ran install_mex, make sure you have installed OpenCL (or, for gpuArray support, the CUDA toolkit) and it can be found on path.'])
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
if options.NLM_use_anatomical && options.use_CPU && options.implementation == 2 && options.NLM
    error('Reference image weighting for NLM is not supported with CPU!')
end
if options.useIndexBasedReconstruction && options.implementation == 2 && options.use_CPU
    error('Index-based reconstruction is not supported on CPU!')
end
% The sensitivity image of list-mode data has to include the same corrections as the measured events. The normalization is
% applied to all detector pairs either with the normalization sinogram that is looked up for each pair (index-based data without
% options.xSens and options.zSens, see loadCorrections) or with a pair table (options.normSens, any list-mode data). Dynamic
% data has the normalization disabled (loadCorrections). The list-mode data does not use options.additionalCorrection, the
% multiplicative correction is the non-subtracted scatter.
if ~options.CT && ~options.SPECT && (options.useIndexBasedReconstruction || options.listmode || listModeCandidate) ...
        && isfield(options, 'compute_sensitivity_image') && options.compute_sensitivity_image
    if numel(options.partitions) > 1
        nFramesEarly = numel(options.partitions);
    elseif isempty(options.partitions)
        nFramesEarly = 1;
    else
        nFramesEarly = options.partitions;
    end
    if options.normalization_correction && nFramesEarly <= 1
        if options.implementation == 4 || (options.implementation == 2 && options.use_CPU)
            error('Normalization correction in the sensitivity image of list-mode data is not supported with the CPU implementation! Disable the normalization correction or use GPU.')
        end
        % The counts of the pair table are taken from options.xSens and options.zSens, if they are given (they replace the scanner
        % coordinates: det_per_ring = numel(xSens) / 2 / nLayers and rings = numel(zSens) / nLayers), and otherwise from the scanner
        userCoord = isfield(options, 'xSens') && isfield(options, 'zSens') && ~isempty(options.xSens) && ~isempty(options.zSens);
        pairTable = isfield(options, 'normSens') && ~isempty(options.normSens) && ~(isfield(options, 'sensXY') && ~isempty(options.sensXY));
        if userCoord || ~options.useIndexBasedReconstruction || pairTable
            nLayersPair = 1;
            if isfield(options, 'nLayers')
                nLayersPair = double(options.nLayers);
            end
            if userCoord
                if mod(numel(options.xSens), 2 * nLayersPair) ~= 0 || mod(numel(options.zSens), nLayersPair) ~= 0
                    error(['options.xSens (' num2str(numel(options.xSens)) ' elements) has to contain two values (x and y) for each detector of all ' ...
                        num2str(nLayersPair) ' layer(s) and options.zSens (' num2str(numel(options.zSens)) ' elements) one value for each ring of all layers!'])
                end
                detPairs = numel(options.xSens) / 2 / nLayersPair;
                ringsPair = numel(options.zSens) / nLayersPair;
                originPair = 'from options.xSens and options.zSens';
            else
                detPairs = double(options.det_per_ring(1));
                ringsPair = double(options.rings(1));
                originPair = 'of the scanner (set options.xSens and options.zSens for a system that is not described by the scanner settings)';
            end
            expectedPairs = detPairs^2 * ringsPair^2 * nLayersPair^2;
            formulaPair = ['det_per_ring^2 * rings^2 * nLayers^2 = ' num2str(detPairs) '^2 * ' num2str(ringsPair) '^2 * ' num2str(nLayersPair) ...
                '^2 = ' num2str(expectedPairs) ' elements (' originPair ')'];
            layoutPair = ['with the index a + b * det_per_ring + (ringA * rings + ringB) * det_per_ring^2 for the detectors a, b and the rings ringA, ' ...
                'ringB (an array of size [det_per_ring, det_per_ring, rings, rings] indexed (a, b, ringB, ringA); with several layers the layer pairs ' ...
                '(layerA * nLayers + layerB) follow as the slowest dimension)'];
            if ~pairTable
                error(['Normalization correction in the sensitivity image of list-mode data with custom detector coordinates or user-input ' ...
                    'sensitivity image coordinates (options.xSens and options.zSens) needs options.normSens, the normalization of all ordered ' ...
                    'detector/ring pairs: ' formulaPair ' ' layoutPair '. Only the pairs of one orientation are read (detector a < b, or ' ...
                    'a == b and ringA < ringB). Alternatively disable the normalization correction.'])
            end
            if ~isnumeric(options.normSens)
                error('options.normSens must be a numeric array!')
            end
            if numel(options.normSens) ~= expectedPairs
                error(['options.normSens has ' num2str(numel(options.normSens)) ' elements, but the pair normalization table must have ' ...
                    formulaPair ' ' layoutPair '!'])
            end
            if ~all(isfinite(options.normSens(:)))
                error('options.normSens (the pair normalization table) contains NaN or infinite values!')
            end
            if any(options.normSens(:) < 0)
                error(['options.normSens (the pair normalization table) contains negative values (minimum ' num2str(min(options.normSens(:))) ')!'])
            end
            options.normSens = single(options.normSens(:));
        end
    end
    if options.scatter_correction && ~options.subtract_scatter
        error(['Multiplicative scatter correction (options.scatter_correction with options.subtract_scatter = false) is not supported in the ' ...
            'sensitivity image of list-mode data! Disable the correction or the sensitivity image.'])
    end
end
if OS_I4_summa > 1
    dispi = [];
    reko = recNames(0);
    for kk = 1 : length(reko)
        if isempty(dispi) && options.(reko{kk})
            dispi = reko{kk};
        elseif options.(reko{kk})
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
% Projector type 2 is rejected also as the forward or backward projector of a combined type (e.g. 12, 21, 23, 24, 42)
if options.CT && (options.projector_type == 2 || (options.projector_type > 10 && (floor(options.projector_type / 10) == 2 || mod(options.projector_type, 10) == 2)))
    error('Projector type 2 (orthogonal distance) is not supported with CT (also as the forward or backward projector of a combined projector type)!')
end
if any(ptDigits == 5) && ~options.CT
    error('Projector type 5 is only supported with CT data!')
end
if any(ptDigits == 6) && ~options.SPECT
    error('Projector type 6 is only supported with SPECT data!')
end
if ~ismember(options.projector_type, [1 2 6 11 16 21 22 26 61 62 66]) && options.SPECT
    error('SPECT only supports projector types 1, 2 and 6, plus supported hybrid variants!')
end
if ~ismember(options.projector_type, [4 5 44 45 54 55]) && options.useHelical
    error('Helical CT only supports projector types 4, 5, 44, 45, 54 and 55!')
end
if options.useParallelBeam && options.useHelical
    error('Parallel beam is not supported with helical CT!')
end
if options.CT && options.useParallelBeam && (mod(options.projector_type, 10) == 5 || (options.projector_type >= 10 && floor(options.projector_type / 10) == 5))
    error('Parallel beam is not supported with projector type 5!')
end
if options.CT && options.useParallelBeam && (options.implementation == 4 || options.use_CPU)
    error('Parallel beam is not supported with the CPU implementation!')
end
if options.subset_type == 12 && ~options.offsetCorrection
    error('Subset type 12 can only be used with offset correction (options.offsetCorrection = true).')
end
if options.offsetCorrection && options.useHelical
    error('Offset correction is not supported with helical CT!')
end
if options.offangle > 0 && options.useHelical
    warning('Rotation is not yet supported with helical CT data')
end
if options.useHelical && (options.implementation == 4 || options.use_CPU)
    error('Helical CT is not supported with the CPU implementation!')
end
if options.useHelical && options.FDK
    error('FDK/FBP is not supported with curved helical CT!')
end
if options.CT && options.useHelical && options.helicalRadius <= 0
    error('Helical reconstruction requires options.helicalRadius (the radius of the circle formed by the curved detector array, in mm) to be set!')
end
if (options.projector_type == 6)
    if options.subsets > 1 && options.subset_type < 8
        error('Subset types 0-7 are not supported with projector type 6!')
    end
end
if (~options.use_raw_data && options.SPECT) && options.subsets > 1 && options.subset_type < 8 && options.projector_type ~= 6
    error('Subset types 0-7 are not supported with SPECT data. Use subset types 8-11.')
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
if options.largeDim
    if ~options.PDHG && ~options.FDK && ~options.PKMA && ~options.PDDY && ~options.PDHGL1 && ~options.PDHGKL
        error('Large dimension support is only available for PDHG, PKMA, and FDK!')
    end
    if options.MRP || options.quad || options.Huber || options.weighted_mean || options.FMH || options.ProxTV || options.TGV || options.L || options.AD || options.APLS
        error('Large dimension support is only available for non-local methods, RDP, GGMRF, hyperbolic prior and TV!')
    end
    if ismember(options.projector_type,[15 51 45 5 54 55 2 3])
        error('Large dimension support currently supports only projector types 14 and 4!')
    end
    % Anatomical weighting is not supported by largeDim
    if options.RDP && options.RDP_use_anatomical
        error('Anatomical weighting for RDP (RDP_use_anatomical) is not supported with largeDim!')
    end
    if options.TV && options.TV_use_anatomical
        error('Anatomical weighting for TV (TV_use_anatomical) is not supported with largeDim!')
    end
    if options.NLM && options.NLM_use_anatomical
        error('Anatomical weighting for NLM (NLM_use_anatomical) is not supported with largeDim!')
    end
    if isfield(options, 'maskBP') && size(options.maskBP, 3) > 1
        error('A 3D backprojection/prior mask is not supported with largeDim! Use a 2D (transaxial) mask instead.')
    end
end
if options.use_CUDA && options.use_CPU && options.implementation == 2
    error('Both CUDA and CPU selected! Select only one!')
end
if options.TOF_bins_used > 1 && (options.projector_type ~= 1 && options.projector_type ~= 11 && options.projector_type ~= 2 && options.projector_type ~= 22 && options.projector_type ~= 12 && options.projector_type ~= 21 ...
        && options.projector_type ~= 23 && options.projector_type ~= 32 && options.projector_type ~= 24 && options.projector_type ~= 42 && options.projector_type ~= 3 && options.projector_type ~= 33 && options.projector_type ~= 31 ...
        && options.projector_type ~= 13 && options.projector_type ~= 4 && options.projector_type ~= 44 && options.projector_type ~= 41 && options.projector_type ~= 14 && options.projector_type ~= 34 && options.projector_type ~= 43) && ~options.CT && ~options.SPECT
    error('TOF is currently only supported with improved Siddon (projector_type = 1), orthogonal distance-based ray tracer (projector_type = 2), volume of intersection (projector_type = 3) and interpolation-based projector (projector_type = 4)')
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
% Randoms and/or scatter correction cannot be applied during the reconstruction with these algorithms
if options.corrections_during_reconstruction && (options.scatter_correction || options.randoms_correction) && (options.PDHG || options.PDHGL1 || options.FISTA || options.LSQR || options.CGLS || options.FISTAL1 || options.PDDY)
    warning('Randoms/scatter correction cannot be applied during the reconstruction with the selected algorithm! The measurement data will be precorrected instead.')
    options.ordinaryPoisson = false;
end
if options.useIndexBasedReconstruction && (options.randoms_correction || (options.scatter_correction && ~options.subtract_scatter)) && options.TOF_bins > 1
    % Subtracted scatter is the additive term of each event, which is allowed with TOF (the value of the TOF bin of the event)
    error('Randoms and/or multiplicative scatter correction cannot be used with index-based reconstruction with TOF data!')
end
if options.CT && (options.COSEM || options.ECOSEM || options.ACOSEM || options.OSL_COSEM > 0)
    error('COSEM, ECOSEM, ACOSEM and OSL-COSEM are not supported with CT data (no transmission tomography formulation).')
end
if options.implementation == 2 && options.use_CPU && options.ECOSEM
    error('ECOSEM is not supported with implementation 2 when using CPU!')
end
if options.implementation == 2 && options.use_CPU
    warning('CPU functionality is limited and might not work correctly in all cases, use at your own risk! It is recommended to install OpenCL libraries for CPU and disable options.use_CPU.')
end
% The normalization correction cannot be used when the normalization coefficients are computed
if options.normalization_correction && options.compute_normalization && ~options.CT
    warning('Normalization correction cannot be applied when computing normalization coefficients. Disabling normalization correction.')
    options.normalization_correction = false;
end
if options.compute_normalization && sum(options.normalization_options) == 0 && ~options.CT
    error('Normalization computation selected, but no normalization components selected.')
end
