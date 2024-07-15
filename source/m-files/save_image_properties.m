function image_properties = save_image_properties(options)
%SAVE_IMAGE_PROPERTIES Saves the various reconstruction related variables
%in a struct variable

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019-2022 Ville-Veikko Wettenhovi
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
if ~isfield(options,'CT')
    options.CT = false;
end

if options.implementation == 1
    image_properties.Implementation = 'MATLAB (1)';
elseif options.implementation == 3
    image_properties.Implementation = 'Matrix-free multi-device OpenCL (3)';
    image_properties.deviceNumber = options.use_device;
elseif options.implementation == 4
    image_properties.Implementation = 'CPU matrix-free (4)';
elseif options.implementation == 2
    image_properties.Implementation = 'Matrix-free ArrayFire OpenCL (2)';
    image_properties.deviceNumber = options.use_device;
end
if options.projector_type == 0
    image_properties.Projector = 'Original Siddon''s ray tracer';
elseif options.projector_type == 1 || options.projector_type == 11
    if options.precompute_lor
        image_properties.Projector = 'Improved Siddon''s ray tracer with 1 ray';
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
        image_properties.Projector = ['Improved Siddon''s algorithm used with ' num2str(options.n_rays_transaxial) ' transaxial ' ray ' and ' ...
            num2str(options.n_rays_axial) ' axial ' aray '.'];
    end
elseif options.projector_type == 2 || options.projector_type == 22
    if options.tube_width_z == 0
        image_properties.Projector = 'Orthogonal distance-based ray tracer in 2D';
    else
        image_properties.Projector = 'Orthogonal distance-based ray tracer in 3D';
        image_properties.orthogonalDistanceZ = options.tube_width_z;
    end
elseif options.projector_type == 3 || options.projector_type == 33
    image_properties.Projector = 'Volume-based ray tracer';
    image_properties.tubeRadius = options.tube_radius;
    image_properties.voxelRadius = options.voxel_radius;
elseif options.projector_type == 4
    image_properties.Projector = 'Interpolation-based ray tracer';
elseif options.projector_type == 14
    image_properties.Projector = 'Improved Siddon for forward projection, interpolation-based for backprojection';
elseif options.projector_type == 41
    image_properties.Projector = 'Improved Siddon for backprojection, interpolation-based for forward projection';
elseif options.projector_type == 15
    image_properties.Projector = 'Improved Siddon for forward projection, branchless distance-driven for backprojection';
elseif options.projector_type == 45
    image_properties.Projector = 'Interpolation-based for forward projection, branchless distance-driven for backprojection';
elseif options.projector_type == 54
    image_properties.Projector = 'Branchless distance-driven for forward projection, interpolation-based for backprojection';
elseif options.projector_type == 51
    image_properties.Projector = 'Branchless distance-driven for forward projection, improved Siddon for backprojection';
elseif options.projector_type == 5
    image_properties.Projector = 'Branchless distance-driven projector';
elseif options.projector_type == 6
    image_properties.Projector = 'SPECT rotation-based projector';
end
image_properties.Projector = [image_properties.Projector ' (' num2str(options.projector_type) ')'];
if options.use_psf
    if options.deblurring
        image_properties.PSF = 'PSF enabled with deblurring';
        image_properties.deblur_iterations = options.deblur_iterations;
    else
        image_properties.PSF = 'PSF enabled';
    end
    image_properties.PSF_FWHM = options.FWHM;
else
    image_properties.PSF = 'No PSF';
end
image_properties.nIter = options.Niter;
image_properties.nSubsets = options.subsets;
if isfield(options,'TOF_bins')
    if options.TOF_bins > 1
        image_properties.TOF = 'TOF enabled';
        image_properties.N_TOF_bins = options.TOF_bins;
        image_properties.TOF_bin_width = options.TOF_width;
        image_properties.TOF_FWHM = options.TOF_FWHM;
    end
end
if ~options.useMultiResolutionVolumes
    if isfield(options,'NxOrig') && options.NxOrig < options.Nx
        image_properties.Nx = options.NxOrig;
        image_properties.NxExtended = options.Nx;
    else
        image_properties.Nx = options.Nx;
    end
    if isfield(options,'NyOrig') && options.NyOrig < options.Ny
        image_properties.Ny = options.NyOrig;
        image_properties.NyExtended = options.Ny;
    else
        image_properties.Ny = options.Ny;
    end
    if isfield(options,'NzOrig') && options.NzOrig < options.Nz
        image_properties.Nz = options.NzOrig;
        image_properties.NzExtended = options.Nz;
    else
        image_properties.Nz = options.Nz;
    end
    if isfield(options,'FOVxOrig') && options.FOVxOrig < options.FOVa_x
        image_properties.FOVx = options.FOVxOrig;
        image_properties.FOVxExtended = options.FOVa_x;
    else
        image_properties.FOVx = options.FOVa_x;
    end
    if isfield(options,'FOVyOrig') && options.FOVyOrig < options.FOVa_y
        image_properties.FOVy = options.FOVyOrig;
        image_properties.FOVyExtended = options.FOVa_y;
    else
        image_properties.FOVy = options.FOVa_y;
    end
    if isfield(options,'axialFOVOrig') && options.axialFOVOrig < options.axial_fov
        image_properties.axialFOV = options.axialFOVOrig;
        image_properties.axialFOVExtended = options.axial_fov;
    else
        image_properties.axialFOV = options.axial_fov;
    end
    image_properties.voxelSizeX = image_properties.FOVx / image_properties.Nx;
    image_properties.voxelSizeY = image_properties.FOVy / image_properties.Ny;
    image_properties.voxelSizeZ = image_properties.axialFOV / image_properties.Nz;
else
    image_properties.Nx = options.Nx(1);
    image_properties.Ny = options.Ny(1);
    image_properties.Nz = options.Nz(1);
    image_properties.NxExtended = options.Nx(2:end);
    image_properties.NyExtended = options.Ny(2:end);
    image_properties.NzExtended = options.Nz(2:end);
    if isfield(options,'FOVxOrig') && options.FOVxOrig < options.FOVa_x(1)
        image_properties.FOVx = options.FOVxOrig;
    else
        image_properties.FOVx = options.FOVa_x;
    end
    image_properties.FOVxExtended = options.FOVa_x(2:end);
    if isfield(options,'FOVyOrig') && options.FOVyOrig < options.FOVa_y(1)
        image_properties.FOVy = options.FOVyOrig;
    else
        image_properties.FOVy = options.FOVa_y;
    end
    image_properties.FOVyExtended = options.FOVa_y(2:end);
    if isfield(options,'axialFOVOrig') && options.axialFOVOrig < options.axial_fov(1)
        image_properties.axialFOV = options.axialFOVOrig;
    else
        image_properties.axialFOV = options.axial_fov;
    end
    image_properties.axialFOVExtended = options.axial_fov(2:end);
    image_properties.voxelSizeX = image_properties.FOVx / double(image_properties.Nx(1));
    image_properties.voxelSizeY = image_properties.FOVy / double(image_properties.Ny(1));
    image_properties.voxelSizeZ = image_properties.axialFOV / double(image_properties.Nz(1));
    image_properties.multiResolutionUsed = true;
    image_properties.nMultiResolutionVolumes = options.nMultiVolumes;
    image_properties.multiResolutionScalingValue = options.multiResolutionScale;
    image_properties.voxelSizeXExtended = image_properties.FOVxExtended ./ double(image_properties.NxExtended);
    image_properties.voxelSizeYExtended = image_properties.FOVyExtended ./ double(image_properties.NyExtended);
    image_properties.voxelSizeZExtended = image_properties.axialFOVExtended ./ double(image_properties.NzExtended);
end
image_properties.name = options.name;
image_properties.scannerName = options.machine_name;
image_properties.nTimeSteps = options.partitions;
image_properties.offsetAngle = options.offangle;
image_properties.flippedImage = options.flip_image;
image_properties.scatterCorrection = options.scatter_correction;
if options.CT
    image_properties.Binning = options.binning;
    image_properties.detectorPitchX = options.dPitchX;
    image_properties.detectorPitchY = options.dPitchY;
    image_properties.nProjections = options.nProjections;
    if isfield(options,'nColsDOrig') && options.nColsDOrig < options.nColsD
        image_properties.nColumnsProjections = options.nColsDOrig;
        image_properties.nColumnsProjectionsExtrapolated = options.nColsD;
    else
        image_properties.nColumnsProjections = options.nColsD;
    end
    if isfield(options,'nRowsDOrig') && options.nRowsDOrig < options.nRowsD
        image_properties.nRowsProjections = options.nRowsDOrig;
        image_properties.nRowsProjectionsExtrapolated = options.nRowsD;
    else
        image_properties.nRowsProjections = options.nRowsD;
    end
    image_properties.sourceOffsetRow = options.sourceOffsetRow;
    image_properties.sourceOffsetCol = options.sourceOffsetCol;
    image_properties.sourceToDetectorDistance = options.sourceToDetector;
    image_properties.sourceToCenterRotation = options.sourceToCRot;
    image_properties.nBedPositions = options.nBed;
    image_properties.bedOffset = options.bedOffset;
    if isfield(options,'flatFieldScaling')
        image_properties.flatFieldScaling = options.flatFieldScaling;
    end
else
    image_properties.nSinogramAngles = options.Nang;
    image_properties.nSinogramRadial = options.Ndist;
    image_properties.nSinograms = options.NSinos;
    image_properties.rawDataUsed = options.use_raw_data;
    image_properties.startTime = options.start;
    image_properties.nRings = options.rings;
    if options.use_raw_data
        image_properties.ringDifference = options.ring_difference_raw;
        if options.sampling_raw > 1
            image_properties.samplingIncrease = options.sampling_raw;
        end
    else
        image_properties.ringDifference = options.ring_difference;
        if options.sampling > 1
            image_properties.samplingIncrease = options.sampling;
        end
    end
    image_properties.detectorsPerRing = options.det_w_pseudo;
    image_properties.boreDiameter = options.diameter;
    image_properties.blocksPerRing = options.blocks_per_ring;
    image_properties.axialBlocks = options.linear_multip;
    image_properties.transaxialBlocks = options.transaxial_multip;
    image_properties.crystalsPerBlock = options.cryst_per_block;
    image_properties.crystalsPerBlockAxial = options.cryst_per_block_axial;
    image_properties.attenuationCorrection = options.attenuation_correction;
    image_properties.arcCorrection = options.arc_correction;
    image_properties.totalTime = options.tot_time;
    image_properties.normalizationCorrection = options.normalization_correction;
    image_properties.randomsCorrection = options.randoms_correction;
    if options.randoms_correction
        image_properties.randomsVarianceReduction = options.variance_reduction;
        image_properties.randomsSmoothing = options.randoms_smoothing;
    end
    image_properties.gapFilling = options.fill_sinogram_gaps;
    if options.reconstruct_trues
        image_properties.reconstructTrues = options.reconstruct_trues;
    end
    if options.reconstruct_scatter
        image_properties.reconstructScatter = options.reconstruct_scatter;
    end
    if isfield(options,'DOI')
        image_properties.DOI = options.DOI;
    end
    if options.use_ASCII
        image_properties.dataType = 'Used ASCII GATE data';
    elseif options.use_LMF
        image_properties.dataType = 'Used LMF GATE data';
    elseif options.use_root
        image_properties.dataType = 'Used ROOT GATE data';
    elseif options.use_machine == 1
        image_properties.dataType = 'Used data obtained from a list-mode file';
    elseif options.use_machine == 2
        image_properties.dataType = 'Used scanner created sinogram data';
    elseif options.use_machine == 3
        image_properties.dataType = 'Used 32-bit list-mode data';
    end
end
if options.subsets > 1
    image_properties.subsetType = options.subset_type;
end
if options.corrections_during_reconstruction && (options.scatter_correction || options.randoms_correction)
    image_properties.ordinaryPoissonReconstruction = true;
else
    image_properties.ordinaryPoissonReconstruction = false;
end
% if isfield(options,'rekot')
%     image_properties.implementations = options.rekot;
% end
if options.ACOSEM || options.OSL_COSEM == 1
    image_properties.hACOSEM = options.h;
end
if options.RAMLA || options.BSREM
    image_properties.relaxationParameter = options.lambda;
end
if options.MRAMLA || options.MBSREM
    image_properties.relaxationParameter = options.lambda;
end
if options.ROSEM || options.ROSEM_MAP
    image_properties.relaxationParameter = options.lambda;
end
if options.DRAMA
    image_properties.relaxationParameter = options.lam_drama;
end
if options.MRAMLA || options.MBSREM
    image_properties.U = options.U;
end
if options.PKMA
    image_properties.rho_PKMA = options.rho_PKMA;
    image_properties.delta_PKMA = options.delta_PKMA;
    % image_properties.delta2_PKMA = options.delta2_PKMA;
    image_properties.relaxationParameter = options.lambda;
    image_properties.MomentumParameterPKMA = options.alpha_PKMA;
end
recA = recNames(0);
recP = recNames(1);
recN = recNames(11);
for kk = 1 : numel(recA)
    if options.(recA{kk})
        image_properties.algorithm = recN{kk};
    end
end
recN = recNames(10);
for kk = 1 : numel(recP)
    if options.(recP{kk})
        image_properties.prior = recN{kk};
    end
end
image_properties.regularizationParameter = options.beta;
if options.quad
    image_properties.weights = options.weights;
    image_properties.quadratic_prior_weights = options.weights_quad;
end
if options.Huber
    image_properties.weights_huber = options.weights_huber;
    image_properties.huber_delta = options.huber_delta;
end
if options.L
    image_properties.L_weights = options.a_L;
end
if options.L || options.quad || options.weighted_mean || options.FMH || options.MRP || options.NLM || options.GGMRF
    image_properties.NeighborhoodSizeX = options.Ndx;
    image_properties.NeighborhoodSizeY = options.Ndy;
    image_properties.NeighborhoodSizeZ = options.Ndz;
end
if options.FMH
    image_properties.fmh_center_weight = options.fmh_center_weight;
    image_properties.fmh_weights = options.fmh_weights;
end
if options.weighted_mean
    image_properties.weighted_mean_center_weight = options.weighted_center_weight;
    image_properties.weighted_mean_weights = options.weighted_weights;
    image_properties.mean_type = options.mean_type;
end
if options.TV
    image_properties.TVsmoothing_beta = options.TVsmoothing;
    image_properties.TV_use_anatomical = options.TV_use_anatomical;
    if options.TV_use_anatomical
%         image_properties.TVtype = options.TVtype;
        image_properties.TV_reference_image = options.TV_reference_image;
        image_properties.T_weighting_parameter = options.T;
        if options.TVtype == 3
            image_properties.C_original_image_weight = options.C;
        end
        if options.TVtype == 1
            image_properties.TVType = 'Anatomically weighted TV (TV type 1)';
        elseif options.TVtype == 2
            image_properties.TVType = 'Anatomically weighted joint TV (TV type 2)';
        elseif options.TVtype == 3
            image_properties.TVType = 'Anatomically weighted joint prior (TV type 3)';
        end
%         image_properties.tau = options.tau;
    else
        if options.TVtype == 1
            image_properties.TVType = 'Regular TV (TV type 1)';
        elseif options.TVtype == 2
            image_properties.TVType = 'Regular TV (TV type 2)';
        elseif options.TVtype == 3
            image_properties.TVType = 'Joint prior (hyperbolic function) (TV type 3)';
            image_properties.HyperbolicWeight = options.SATVPhi;
        elseif options.TVtype == 4
            image_properties.TVType = 'Smoothed anisotropic TV (Lange function) (TV type 4)';
            image_properties.LangePhi = options.SATVPhi;
        elseif options.TVtype == 6
            image_properties.TVType = 'Weighted TV (TV type 6)';
        end
    end

end
if options.AD
    image_properties.TimeStepAD = options.TimeStepAD;
    image_properties.KAD = options.KAD;
    image_properties.NiterAD = options.NiterAD;
    if options.FluxType == 1
        image_properties.FluxTypeAD = 'exponential';
    else
        image_properties.FluxTypeAD = 'quadratic';
    end
    image_properties.DiffusionTypeAD = options.DiffusionType;
end
if options.APLS
    image_properties.APLSsmoothing_beta = options.APLSsmoothing;
    image_properties.APLS_reference_image = options.APLS_reference_image;
    image_properties.eta = options.eta;
end
if options.TGV
    image_properties.alpha0TGV = options.alpha0TGV;
    image_properties.alpha1TGV = options.alpha1TGV;
%     if options.TGV
%         image_properties.NiterTGV = options.NiterTGV;
%     end
end
if options.NLM
    image_properties.NLPatchSizeX = options.Nlx;
    image_properties.NLPatchSizeY = options.Nly;
    image_properties.NLPatchSizeZ = options.Nlz;
    image_properties.NLFilterParameter = options.sigma;
    image_properties.NLM_gaussSTD = options.NLM_gauss;
    image_properties.NLM_use_anatomical = options.NLM_use_anatomical;
    if options.NLM_MRP
        image_properties.NLMType = 'NLM with NLM filtered image subtracted from the input image (NLM_MRP)';
    elseif options.NLTV
        image_properties.NLMType = 'Non-local TV (NLTV)';
    elseif options.NLRD
        image_properties.NLMType = 'Non-local relative difference prior (NLRDP)';
        image_properties.RDP_gamma = options.RDP_gamma;
    elseif options.NLLange
        image_properties.NLMType = 'Non-local Lange function prior (NLL)';
        image_properties.LangePhi = options.SATVPhi;
    else
        image_properties.NLMType = 'MRF-type NLM prior';
    end
end
if options.PDHG || options.PDHGKL || options.PDHGL1 || options.FISTA || options.FISTAL1
    image_properties.tauCP = options.tauCP;
    image_properties.sigmaCP = options.sigmaCP;
    image_properties.thetaCP = options.thetaCP;
    if options.precondTypeMeas(2) == 1
        image_properties.tauCPFilteredIterations = options.tauCPFilt;
    end
end
if options.RDP
    image_properties.RDP_gamma = options.RDP_gamma;
    image_properties.RDPIncludeCorners = options.RDPIncludeCorners;
end
if options.GGMRF
    image_properties.GGMRF_p = options.GGMRF_p;
    image_properties.GGMRF_q = options.GGMRF_q;
    image_properties.GGMRF_c = options.GGMRF_c;
    image_properties.GGMRFweights = options.weights_quad;
end
if (options.PDHG || options.PDHGKL || options.PDHGL1) && (options.ProxTV || options.TGV)
    image_properties.useL2Ball = options.useL2Ball;
end
precondIm = {'Diagonal';'EM';'IEM';'Momentum';'Normalized gradient';'Filtering';'Curvature'};
precondMeas = {'Diagonal';'Filtering'};
if any(options.precondTypeImage)
    image_properties.imageBasedPreconditioners = precondIm(logical(options.precondTypeImage));
end
if options.precondTypeMeas(2) == 1
    image_properties.filteringIterations = options.filteringIterations;
    image_properties.filterWindow = options.filterWindow;
    if options.PKMA || options.MBSREM || options.MRAMLA || options.SPS
        image_properties.relaxationParameterFiltering = options.lambdaFiltered;
    end
end
if options.precondTypeImage(4) == 1
    image_properties.Precond4gradV1 = options.gradV1;
    image_properties.Precond4gradV2 = options.gradV2;
    image_properties.Precond4gradInitIter = options.gradInitIter;
end
if any(options.precondTypeMeas)
    image_properties.measurementBasedPreconditioners = precondMeas(logical(options.precondTypeMeas));
end
if options.oOffsetZ ~= 0 || options.oOffsetX ~= 0 || options.oOffsetY ~= 0
    image_properties.objectOffsetsXYZ = [options.oOffsetX; options.oOffsetY; options.oOffsetZ];
end
image_properties.enforcePositivity = options.enforcePositivity;
if options.implementation == 2
    image_properties.savedIterationNumbers = options.saveNIter + 1;
    image_properties.ImagesWereUsed = options.useImages;
    image_properties.MADWasUsed = options.useMAD;
end
% pz{end,1} = image_properties;
