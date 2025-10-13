function [pz, t, recPar] = PET_TOF_example(type, varargin)
%% Function used to call each example individually
% If you want a standalone script for the examples, see the PET_exampleX.m
% files, where X=1 corresponds to the left figure, X=2 the center, etc.

clear mex

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SCANNER PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% R-sectors/modules/blocks/buckets in transaxial direction
options.blocks_per_ring = (38);

%%% R-sectors/modules/blocks/buckets in axial direction (i.e. number of physical
%%% scanner/crystal rings) 
% Multiplying this with the below cryst_per_block_axial should equal the
% total number of crystal rings. 
options.linear_multip = (8);

%%% Number of detectors on the side of R-sector/block/module (transaxial
%%% direction)
% (e.g. 13 if 13x13, 20 if 20x10)
options.cryst_per_block = (20);

%%% Number of detectors on the side of R-sector/block/module (axial
%%% direction)
% (e.g. 13 if 13x13, 10 if 20x10)
options.cryst_per_block_axial = 10;

%%% Crystal pitch/size in x- and y-directions (transaxial) (mm)
options.cr_p = 3.2;

%%% Crystal pitch/size in z-direction (axial) (mm)
options.cr_pz = 3.2;

%%% Ring diameter (distance between perpendicular detectors) (mm)
options.diameter = 780;

%%% Transaxial FOV size (mm), this is the length of the x (vertical/row) side
% of the FOV
options.FOVa_x = (704/sqrt(2));

%%% Transaxial FOV size (mm), this is the length of the y (horizontal/column) side
% of the FOV
options.FOVa_y = options.FOVa_x;

%%% Axial FOV (mm)
options.axial_fov = 261;

%%% Ring gaps (mm)
% The length of the gap between each ring, as defined by
% options.linear_multip
options.ringGaps = repmat(0.7142857,7,1);

%%% Number of detectors per ring (without pseudo detectors)
options.det_per_ring = options.blocks_per_ring*options.cryst_per_block;

%%% Number of detectors per ring (with pseudo detectors)
% NOTE: Vision normally has one pseudo detector per block, but it is
% omitted here by default
options.det_w_pseudo = options.blocks_per_ring*(options.cryst_per_block);

%%% Number of crystal rings
options.rings = options.linear_multip * options.cryst_per_block_axial;

%%% Total number of detectors
options.detectors = options.det_per_ring*options.rings;

%%% Scanner name
options.machine_name = 'Siemens_Vision';
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% IMAGE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Note that non-square transaxial image sizes can be unreliable just as the
% non-square transaxial FOV, but they should, generally, work
%%% Reconstructed image pixel count (X/row-direction)
options.Nx = 220;

%%% Y/column-direction
options.Ny = 220;

%%% Z-direction (number of slices) (axial)
options.Nz = 88;

%%% Flip the image (in column direction)?
options.flip_image = true;

%%% How much is the image rotated?
% NOTE: The rotation is done in the detector space (before reconstruction).
% This current setting is for scanner list-mode data.
% Positive values perform the rotation in clockwise direction
% The unit is the detector element, i.e. the below 390 means that the
% rotation is done by the size of 390 detector elements
options.offangle = 390;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SINOGRAM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% Span factor/axial compression
options.span = 1;

%%% Maximum ring difference
options.ring_difference = options.rings - 1;

%%% Number of radial positions (views) in sinogram
% You should primarily use the same number as the scanner uses.
% However, if that information is not available you can use ndist_max
% function to determine potential values (see help ndist_max for usage).
% This is the ROW dimension, i.e. the number of rows in the sinogram
options.Ndist = 520;

%%% Number of angles (tangential positions) in sinogram 
% This is the final amount after possible mashing, maximum allowed is the
% number of detectors per ring/2.
% This is the COLUMN dimension, i.e. the number of columns in the sinogram
% NOTE: On Vision, the tangential angles are mashed significantly with the use
% of TOF data. This is not possible in OMEGA.
% Using mashing of 4
options.Nang = options.det_w_pseudo/8;

%%% Specify the amount of sinograms contained on each segment
% (this should total the total number of sinograms).
% Currently this is computed automatically, but you can also manually
% specify the segment sizes.
options.segment_table = [options.rings*2-1, options.rings*2-1 - (options.span + 1):-options.span*2:max(options.Nz - options.ring_difference*2, options.rings - options.ring_difference)];
if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','8.5')
    options.segment_table = [options.segment_table(1), repeat_elem(options.segment_table(2:end),2,1)];
else
    options.segment_table = [options.segment_table(1), repelem(options.segment_table(2:end),2)];
end

%%% Total number of sinograms
options.TotSinos = sum(options.segment_table);

%%% Number of sinograms used in reconstruction
% The first NSinos sinograms will be used for the image reconstruction.
% This is the number of slices in the sinogram
options.NSinos = options.TotSinos;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CORRECTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%% Randoms correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If set to true, performs randoms correction during reconstruction or 
% performs precorrection, depending on the selections below
% If you are using your own data, the randoms
% data can be input either manually into options.SinDelayed or input when
% prompted (has to be stored in a mat-file beforehand!)
if type < 2
    options.randoms_correction = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%% Attenuation correction %%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Image-based attenuation correction
% Include attenuation correction from images (e.g. CT-images) (for this you
% need attenuation images of each slice correctly rotated and scaled for
% 511 keV) or from attenuation sinograms. Note that all the attenuation
% data has to be correctly scaled before reconstruction.
options.attenuation_correction = true;

%%% Attenuation image data file
% Specify the path (if not in MATLAB path) and filename.
% NOTE: the attenuation data must be the only variable in the file and
% have the dimensions of the final reconstructed image.
% Alternatively, just input the attenuation data into options.vaimennus
options.attenuation_datafile = '511keV_attenuation_coefficients_for_VisionDerenzoAttenuation.mat';
 
%%%%%%%%%%%%%%%%%%%% Corrections during reconstruction %%%%%%%%%%%%%%%%%%%%
% If set to true, all the corrections are performed during the
% reconstruction step, otherwise the corrections are performed to the
% sinogram/raw data before reconstruction (i.e. precorrected). I.e. this
% can be considered as e.g. normalization weighted reconstruction if
% normalization correction is applied.
% NOTE: Attenuation correction is always performed during reconstruction
% regardless of the choice here.
% If you have manually precorrected the data, do not put those corrections
% to true that have already been applied! Otherwise, the data will be 
% precorrected twice. This obviously only applies when this is set to false
options.corrections_during_reconstruction = true;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TOF PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Total number of TOF bins
options.TOF_bins = 33;
% options.TOF_bins = 1;

%%% Length of each TOF bin (s)
% The time length of each TOF bin in seconds
% This multiplied with the number of bins total the entire time frame that
% the TOF data contains. For example with 10 bins of size 400 ps all time
% differences of at most 4 ns will be included in the TOF data. The
% multiplied value should be, at most, the size of the coincidence window.
options.TOF_width = 143e-12;

%%% TOF offset (s)
% If your TOF bins are not centered on zero (center of FOV) you can specify
% the offset value here.
options.TOF_offset = 4.3333e-12;

%%% FWHM of the temporal noise/FWHM of the TOF data (s)
% This parameter has two properties. The first one applies to any TOF data
% that is saved by OMEGA (GATE, Inveon/Biograph list-mode), the second only
% to GATE data.
% Firstly this specifies the FWHM of TOF data used for file naming and
% loading purposes. This value is included in the filename when data is
% imported/saved and also used when that same data is later loaded. 
% Secondly, this is the FWHM of the ADDED temporal noise to the time
% differences. If you are using GATE data and have set a custom temporal
% blurring in GATE then you should set to this zero if you wish to use the
% same temporal resolution. If no custom temporal blurring was applied then
% use this value to control the accuracy of the TOF data. For example if
% you want to have TOF data with 500 ps FWHM then set this value to
% 500e-12. 
options.TOF_noise_FWHM = 214e-12;

%%% FWHM of the TOF data (s)
% Applies to ALL data.
% This value specifies the TOF accuracy during the reconstruction process
% and thus can be different from above. If you are using GATE data with
% temporal blurring, you need to multiply that FWHM with sqrt(2) here.
options.TOF_FWHM = 214e-12;

%%% Number of TOF bins used in reconstruction
% Number of TOF bins used during reconstruction phase.
% NOTE: Currently supports only either all bins specified by
% options.TOF_bins or 1 (or 0) which converts the TOF data into non-TOF
% data during reconstruction phase.
options.TOF_bins_used = options.TOF_bins;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MISC PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% Name of current datafile/examination
% Used in general for naming purposes, but not utilized here
options.name = 'Vision_derenzo';

%%% Show status messages
% These are e.g. time elapsed on various functions and what steps have been
% completed. It is recommended to keep this at 1 or 2. With value of 2, 
% you get more detailed timing information. Maximum is 3, minimum 0.
options.verbose = 2;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPLEMENTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reconstruction implementation used
% 1 = Reconstructions in MATLAB (projector in a MEX-file), uses matrices.
% (Slow and memory intensive)
% 2 = Matrix-free reconstruction with OpenCL/CUDA (Recommended)
% (Requires ArrayFire).
% 3 = Multi-GPU/device matrix-free OpenCL (OSEM & MLEM only).
% 4 = Matrix-free reconstruction with OpenMP (CPU, parallel), standard C++
% 5 = Matrix-free reconstruction with OpenCL (parallel)
% See the docs for more information: 
% https://omega-doc.readthedocs.io/en/latest/implementation.html
options.implementation = 2;

% Applies to implementations 2, 3 and 5 ONLY
%%% OpenCL/CUDA device used 
% NOTE: Use ArrayFire_OpenCL_device_info() to determine the device numbers
% with implementation 2.
% NOTE: Use OpenCL_device_info() to determine the platform numbers and
% their respective devices with implementations 3 or 5.
% NOTE: The device numbers might be different between implementation 2 and
% implementations 3 and 5
% NOTE: if you switch devices then you need to run the below line
% (uncommented) as well:
% clear mex
if ~isempty(varargin)
    options.use_device = varargin{1};
else
    options.use_device = 0;
end

% Implementation 2 ONLY
%%% Use CUDA
% Selecting this to true will use CUDA kernels/code instead of OpenCL. This
% only works if the CUDA code was successfully built. Recommended only for
% Siddon as the orthogonal/volume-based ray tracer are slower in CUDA.
options.use_CUDA = checkCUDA(options.use_device);

% Applies to implementations 2, 3 and 5 ONLY
%%% Use 32-bit integer atomic functions
% If true, then 32-bit integer atomic functions (atomic add) will be used
% if they are supported by the selected device.
% Setting this to true will make computations faster on GPUs that support
% the functions, but might make results slightly less reliable due to
% floating point rounding. Recommended for GPUs.
if ~options.use_CUDA
    options.use_32bit_atomics = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROJECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Type of projector to use for the geometric matrix
% 1 = Improved/accelerated Siddon's algorithm
% 2 = Orthogonal distance based ray tracer
% 3 = Volume of intersection based ray tracer
% 4 = Interpolation-based projector
% NOTE: You can mix and match most of the projectors. I.e. 41 will use
% interpolation-based projector for forward projection while improved
% Siddon is used for backprojection.
% NOTE 2: The below additional options apply also in hybrid cases as long
% as the other projector is the corresponding projector.
% See the doc for more information:
% https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
if type ~= 1
    options.projector_type = 1;
else
    options.projector_type = 3;
end

%%% Interpolation length (projector type = 4 only)
% This specifies the length after which the interpolation takes place. This
% value will be multiplied by the voxel size which means that 1 is
% the interpolation length corresponding to a single voxel (transaxial)
% length. Larger values lead to faster computation but at the cost of
% accuracy. Recommended values are between [0.5 1].
options.dL = 1;

%%% Use point spread function (PSF) blurring
% Applies PSF blurring through convolution to the image space. This is the
% same as multiplying the geometric matrix with an image blurring matrix.
if type ~= 1
    options.use_psf = true;
else
    options.use_psf = false;
end

% FWHM (mm) of the Gaussian used in PSF blurring in all three dimensions
options.FWHM = [options.cr_p*1. options.cr_p*1. options.cr_pz]*1.;

% Orthogonal ray tracer (projector_type = 2 only)
%%% The 2D (XY) width (mm) of the "strip/tube" where the orthogonal distances are
% included. If tube_width_z below is non-zero, then this value is ignored.
options.tube_width_xy = options.cr_p;

% Orthogonal ray tracer (projector_type = 2 only)
%%% The 3D (Z) width (mm) of the "tube" where the orthogonal distances are
% included. If set to 0, then the 2D orthogonal ray tracer is used. If this
% value is non-zero then the above value is IGNORED.
% If you want the projector to be a tube, use this, if you want it to be 
% strip, use the above
% This slows down the reconstruction, but makes it more accurate
options.tube_width_z = options.cr_pz;

% Volume ray tracer (projector_type = 3 only)
%%% Radius of the tube-of-response (cylinder)
% The radius (mm) of the cylinder that approximates the tube-of-response.
options.tube_radius = sqrt(2) * (options.cr_pz / 2);

% Volume ray tracer (projector_type = 3 only)
%%% Relative size of the voxel (sphere)
% In volume ray tracer, the voxels are modeled as spheres. This value
% specifies the relative radius of the sphere such that with 1 the sphere
% is just large enough to encompass an entire cubic voxel, i.e. the
% corners of the cubic voxel intersect with the sphere shell. Larger values
% create larger spheres, while smaller values create smaller spheres.
options.voxel_radius = 1;
 
%%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Number of iterations (all reconstruction methods)
options.Niter = 1;

%%% Save specific intermediate iterations
% You can specify the intermediate iterations you wish to save here. Note
% that this uses zero-based indexing, i.e. 0 is the first iteration (not
% the initial value). By default only the last iteration is saved.
options.saveNIter = [];
% Alternatively you can save ALL intermediate iterations by setting the
% below to true and uncommenting it
% options.save_iter = false;

%%% Number of subsets (excluding subset_type = 6)
options.subsets = 30;

%%% Subset type (n = subsets)
% 1 = Every nth (column) measurement is taken
% 2 = Every nth (row) measurement is taken (e.g. if subsets = 3, then
% first subset has measurements 1, 4, 7, etc., second 2, 5, 8, etc.) 
% 3 = Measurements are selected randomly
% 4 = (Sinogram only) Take every nth column in the sinogram
% 5 = (Sinogram only) Take every nth row in the sinogram
% 6 = Sort the LORs according to their angle with positive X-axis, combine
% n_angles together and have 180/n_angles subsets for 2D slices and
% 360/n_angles for 3D, see docs for more information:
% https://omega-doc.readthedocs.io/en/latest/algorithms.html#type-6
% 7 = Form the subsets by using golden angle sampling
% 8 = Use every nth sinogram
% 9 = Randomly select the full sinograms
% 10 = Use golden angle sampling to select the subsets (not recommended for
% PET)
% 11 = Use prime factor sampling to select the full sinograms
% Most of the time subset_type 1 or 4 is sufficient.
if type == 2
    options.subset_type = 0;
else
    options.subset_type = 1;
end

%%% Initial value for the reconstruction
options.x0 = ones(options.Nx, options.Ny, options.Nz);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION ALGORITHMS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction algorithms to use (choose only one algorithm and
% optionally one prior)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ML-METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are non-regularized versions

%%% Ordered Subsets Expectation Maximization (OSEM) OR Maximum-Likelihood
%%% Expectation Maximization (MLEM) (if subsets = 1)
% Supported by all implementations
if type < 2
    options.OSEM = true;
else
    options.OSEM = false;
end

%%% Modified Row-Action Maximum Likelihood Algorithm (MRAMLA)
% Supported by implementations 1, 2, 4, and 5
options.MRAMLA = false;

%%% Row-Action Maximum Likelihood Algorithm (RAMLA)
% Supported by implementations 1, 2, 4, and 5
options.RAMLA = false;

%%% Relaxed Ordered Subsets Expectation Maximization (ROSEM)
% Supported by implementations 1, 2, 4, and 5
options.ROSEM = false;

%%% Rescaled Block Iterative Expectation Maximization (RBI-EM)
% Supported by implementations 1, 2, 4, and 5
options.RBI = false;

%%% Dynamic RAMLA (DRAMA)
% Supported by implementations 1, 2, 4, and 5
options.DRAMA = false;

%%% Complete data OSEM (COSEM)
% Supported by implementations 1, 2, 4, and 5
options.COSEM = false;

%%% Enhanced COSEM (ECOSEM)
% Supported by implementations 1, 2, 4, and 5
options.ECOSEM = false;

%%% Accelerated COSEM (ACOSEM)
% Supported by implementations 1, 2, 4, and 5
options.ACOSEM = false;

%%% FISTA
% Supported by implementations 1, 2, 4, and 5
options.FISTA = false;

%%% FISTA with L1 regularization (FISTAL1)
% Supported by implementations 1, 2, 4, and 5
options.FISTAL1 = false;

%%% LSQR
% Supported by implementations 1, 2, 4, and 5
options.LSQR = false;

%%% CGLS
% Supported by implementations 1, 2, 4, and 5
options.CGLS = false;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAP-METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Any algorithm selected here will utilize any of the priors selected below
% this. Note that only one algorithm and prior combination is allowed! You
% can also use most of these algorithms without priors (such as PKMA or
% PDHG).
%%% One-Step Late OSEM (OSL-OSEM) or MLEM (if subsets = 1)
% Supported by implementations 1, 2, 4, and 5
options.OSL_OSEM = false;

%%% Modified BSREM (MBSREM)
% Supported by implementations 1, 2, 4, and 5
options.MBSREM = false;

%%% Block Sequential Regularized Expectation Maximization (BSREM)
% Supported by implementations 1, 2, 4, and 5
options.BSREM = false;

%%% ROSEM-MAP
% Supported by implementations 1, 2, 4, and 5
options.ROSEM_MAP = false;

%%% RBI-OSL
% Supported by implementations 1, 2, 4, and 5
options.OSL_RBI = false;

%%% (A)COSEM-OSL
% 0/false = No COSEM-OSL, 1/true = ACOSEM-OSL, 2 = COSEM-OSL
% Supported by implementations 1, 2, 4, and 5
options.OSL_COSEM = false;

%%% Preconditioned Krasnoselskii-Mann algorithm (PKMA)
% Supported by implementations 1, 2, 4, and 5
if type > 1
    options.PKMA = true;
else
    options.PKMA = false;
end

%%% Primal-dual hybrid gradient (PDHG)
% Supported by implementations 1, 2, 4, and 5
options.PDHG = false;

%%% Primal-dual hybrid gradient (PDHG) with L1 minimization
% Supported by implementations 1, 2, 4, and 5
options.PDHGL1 = false;

%%% Primal-dual hybrid gradient (PDHG) with Kullback-Leibler minimization
% Supported by implementations 1, 2, 4, and 5
options.PDHGKL = false;

%%% Primal-dual Davis-Yin (PDDY)
% Supported by implementation 2
options.PDDY = false;

options.SAGA = false;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Median Root Prior (MRP)
options.MRP = false;

%%% Quadratic Prior (QP)
options.quad = false;

%%% Huber Prior (QP)
options.Huber = false;

%%% L-filter prior
options.L = false;

%%% Finite impulse response (FIR) Median Hybrid (FMH) prior
options.FMH = false;

%%% Weighted mean prior
options.weighted_mean = false;

%%% Total Variation (TV) prior
options.TV = false;

%%% Anisotropic Diffusion Median Root Prior (ADMRP)
options.AD = false;

%%% Asymmetric Parallel Level Set (APLS) prior
options.APLS = false;

%%% Hyperbolic prior
options.hyperbolic = false;

%%% Total Generalized Variation (TGV) prior
options.TGV = false;

%%% Non-local Means (NLM) prior
options.NLM = false;

%%% Relative difference prior
if type > 1
    options.RDP = true;
else
    options.RDP = false;
end

%%% Generalized Gaussian Markov random field (GGMRF) prior
options.GGMRF = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%% ENFORCE POSITIVITY %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Applies to PDHG, PDHGL1, PDDY, FISTA, FISTAL1, MBSREM, MRAMLA, PKMA
% Enforces positivity in the estimate after each iteration
options.enforcePositivity = true;


%%%%%%%%%%%%%%%%%%%%%%%%%% RELAXATION PARAMETER %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Relaxation parameter for MRAMLA, RAMLA, ROSEM, BSREM, MBSREM and PKMA
% Use scalar if you want it to decrease as
% lambda / ((current_iteration - 1)/20 + 1). Use vector (length = Niter) if
% you want your own relaxation parameters. Use empty array or zero if you
% want to OMEGA to compute the relaxation parameter using the above formula
% with lambda = 1. Note that current_iteration is one-based, i.e. it starts
% at 1.
options.lambda = 1;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PKMA PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step size (alpha) parameter for PKMA
% If a scalar (or an empty) value is used, then the alpha parameter is
% computed automatically as alpha_PKMA(oo) = 1 + (options.rho_PKMA *((i -
% 1) * options.subsets + ll)) / ((i - 1) * options.subsets + ll +
% options.delta_PKMA), where i is the iteration number and l the subset
% number. The input number thus has no effect. options.rho_PKMA and
% options.delta_PKMA are defined below.
% If, on the other hand, a vector is input then the input alpha values are
% used as is without any modifications (the length has to be at least the
% number of iterations * number of subsets).
options.alpha_PKMA = 0;

%%% rho_PKMA
% This value is ignored if a vector input is used with alpha_PKMA
options.rho_PKMA = 0.95;

%%% delta_PKMA
% This value is ignored if a vector input is used with alpha_PKMA
options.delta_PKMA = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRECONDITIONERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Applies to PDHG, PDHGL1, PDHGKL, PKMA, MBSREM, MRAMLA, PDDY, FISTA and
%%% FISTAL1
% Measurement-based preconditioners
% precondTypeMeas(1) = Diagonal normalization preconditioner (1 / (A1))
% precondTypeMeas(2) = Filtering-based preconditioner
options.precondTypeMeas = [false;false];

% Image-based preconditioners
% Setting options.precondTypeImage(2) = true when using PKMA, MRAMLA or
% MBSREM is recommended
% precondTypeImage(1) = Diagonal normalization preconditioner (division with
% the sensitivity image 1 / (A^T1), A is the system matrix) 
% precondTypeImage(2) = EM preconditioner (f / (A^T1), where f is the current
% estimate) 
% precondTypeImage(3) = IEM preconditioner (max(n, fhat, f)/ (A^T1), where
% fhat is an estimate of the final image and n is a small positive number) 
% precondTypeImage(4) = Momentum-like preconditioner (basically a step size
% inclusion) 
% precondTypeImage(5) = Gradient-based preconditioner (Uses the normalized
% divergence (sum of the gradient) of the current estimate) 
% precondTypeImage(6) = Filtering-based preconditioner
% precondTypeImage(7) = Curvature-based preconditioner
if type > 1
    options.precondTypeImage = [false;true;false;false;false;false;false];
end

% Number of filtering iterations
% Applies to both precondTypeMeas(2) and precondTypeImage(6)
options.filteringIterations = 10;


%%%%%%%%%%%%%%%%%%%%%%%%% REGULARIZATION PARAMETER %%%%%%%%%%%%%%%%%%%%%%%%
%%% The regularization parameter for ALL regularization methods (priors)
if type == 2
    options.beta = 0.1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RDP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Edge weighting factor
options.RDP_gamma = 1;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if type < 2
    load SiemensVision_DerenzoPhantom_TOF214psFWHM_33bins_listmode.mat SinDelayed
    options.SinDelayed = SinDelayed;
    clear SinDelayed
    if isempty(varargin) || nargin == 2 || (nargin > 2 && isempty(varargin{2}))
        fid = fopen('SiemensVision_DerenzoPhantom_TOF214psFWHM_33bins_sinogram_520x95x6400_span1.bin');
    else
        fid = fopen(fullfile(varargin{2}, 'SiemensVision_DerenzoPhantom_TOF214psFWHM_33bins_sinogram_520x95x6400_span1.bin'));
    end
    sinogram = fread(fid, inf, 'uint8=>uint8');
    fclose(fid);
    sinogram = reshape(sinogram, [options.Ndist, options.Nang, options.rings^2, options.TOF_bins]);
    options.SinM = sinogram;
    clear sinogram
else
    load SiemensVision_DerenzoPhantom_TOF214psFWHM_33bins_listmode.mat xy z trIndex axIndex TOFIndices
    options.trIndex = trIndex;
    clear trIndex
    options.axIndex = axIndex;
    clear axIndex
    options.TOFIndices = TOFIndices;
    options.x = xy;
    options.z = z;
    clear TOFIndices xy z
end

options.epps = 1e-3;

if type == 2
    options.useIndexBasedReconstruction = true;
else
    options.useIndexBasedReconstruction = false;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ERROR CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic error checking is done here
options = OMEGA_error_check(options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if options.useIndexBasedReconstruction
    options.SinM = ones(size(options.trIndex, 2),1,'uint8');
    options.compute_sensitivity_image = true;
end

%% Reconstructions

% Only load the current subset into the device (GPU) memory
options.loadTOF = false;


tStart = tic;
[pz,recPar] = reconstructions_main(options);
t = toc(tStart);
disp(['Reconstruction process took ' num2str(t) ' seconds'])

[columnsInImage, rowsInImage] = meshgrid(1:options.Nx, 1:options.Ny);
centerX = options.Nx/2;
centerY = options.Ny/2;
radius = options.Nx/2;
mask = ((rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2);


pz(repmat(~mask,1,1,size(pz,3))) = 0;


volume3Dviewer(pz, [], [0 0 1])