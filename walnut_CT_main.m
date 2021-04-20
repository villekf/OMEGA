%% MATLAB codes for CT reconstruction using FIPS walnut projection images
% This example script uses the FIPS walnut projection images
% (https://zenodo.org/record/1254206). This is very similar to
% walnut2D_CT_main.m, but instead of being 2D this is a full 3D example.
% Same examples are included, i.e. this script contains both the structure
% for built-in reconstruction as well as examples on how to use the
% forward/backward projection class.
% NOTE: This script uses implementation 4 (parallel CPU computation with
% OpenMP) by default. To use GPU accelerated reconstruction use either
% implementation 2 or 3.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SCANNER PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Binning
% The level of binning used for the raw data. For example binning of 2
% reduces the size of the projections by two from both dimensions (e.g.
% 2048x3072 becomes 1024x1536).
options.binning = 4;

%%% Number of detector pixels (horizontal)
% The number of detector pixels in the detector panel (horizontal
% direction)
% NOTE: if you use binning, this value has to use the final binned
% dimensions
options.xSize = 2368/options.binning;

%%% Number of detector pixels (vertical)
% The number of detector pixels in the detector panel (vertical
% direction)
% NOTE: if you use binning, this value has to use the final binned
% dimensions
options.ySize = 2240/options.binning;

%%% Number of projections
% Total number of projections used
options.nProjections = 721;

%%% Projection angles (degree or radian)
% The angles corresponding to the projections
options.angles = -linspace(0, 360, options.nProjections);

%%% Detector pixel pitch/size (mm)
% The size of the detector/distance between adjacent detectors
% NOTE: if you use binning, this value has to use the final binned
% dimensions
options.dPitch = 0.05*options.binning;

%%% Source to detector distance (mm)
% The orthogonal distance from the source to the detector panel
options.sourceToDetector = 553.74;

%%% Source to center of rotation distance (mm)
% The distance from the source to the center of rotation/object/origin
options.sourceToCRot = 210.66;

%%% Name of current datafile/examination
% This is used for naming purposes only
options.name = 'Walnut3DCT_data';

%%% Compute only the reconstructions
% If this file is run with this set to true, then the data load will be
% skipped if the options.SinM variable exists
options.only_reconstructions = false;

%%% Show status messages
% These are e.g. time elapsed on various functions and what steps have been
% completed. It is recommended to keep this true.
options.verbose = true;

%%% Transaxial FOV size (mm), this is the length of the x (horizontal) side
% of the FOV
options.FOVa_x = 40.1;

%%% Transaxial FOV size (mm), this is the length of the y (vertical) side
% of the FOV
options.FOVa_y = options.FOVa_x;

%%% Axial FOV (mm)
options.axial_fov = 40;

%%% Source horizontal offset
% The center of rotation is not exactly in the origin. With this parameter
% the source location can be offset by the specifed amount (horizontally).
% This has a similar effect as circulary shifting the projection images.
% NOTE: The default value has been obtained experimentally and is not based
% on any known value.
options.horizontalOffset = -0.05*3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if ~options.only_reconstructions || ~isfield(options,'SinM')
    options.SinM = loadProjectionImages(options.nProjections,options.binning);
    options.SinM = single(options.SinM) ./ single(max(max(max(options.SinM(4:end-3,:,:)))));
    options.SinM(options.SinM > 1) = single(1);
    options.SinM = permute(options.SinM, [2 1 3]);
end
% NOTE: If you want to reduce the number of projections, you need to do
% this manually as outlined below:
% options.SinM = options.SinM(:,:,1:4:options.nProjections);
% options.angles = options.angles(1:4:numel(options.angles));
% options.nProjections = numel(options.angles);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% IMAGE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Reconstructed image pixel size (X-direction)
options.Nx = 280;

%%% Y-direction
options.Ny = 280;

%%% Z-direction (number of slices) (axial)
options.Nz = 140;

%%% Flip the image (in vertical direction)?
options.flip_image = false;

%%% How much is the image rotated (radians)?
% The angle (in radians) on how much the image is rotated BEFORE
% reconstruction, i.e. the rotation is performed in the detector space.
options.offangle = (3*pi)/2;

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
% 2 = Matrix-free reconstruction with OpenCL/ArrayFire (Recommended)
% (Requires ArrayFire. Compiles with MinGW ONLY when ArrayFire was compiled
% with MinGW as well (cannot use the prebuilt binaries)).
% 3 = Multi-GPU/device matrix-free OpenCL (OSEM & MLEM only).
% 4 = Matrix-free reconstruction with OpenMP (parallel), standard C++
% (Supports only one algorithm at a time)
% See the wiki for more information:
% https://github.com/villekf/OMEGA/wiki/Useful-information#selecting-the-correct-implementation
options.implementation = 4;

% Applies to implementations 2 and 3 ONLY
%%% Device used (this is applicable to implementation 2), or platform used
% (implementation 3)
% In implementation 2 this determines the device used for both system
% matrix formation and image reconstruction.
% NOTE: Use ArrayFire_OpenCL_device_info() to determine the device numbers.
% In implementation 3, this determines the platform from where the
% device(s) are taken.
% NOTE: Use OpenCL_device_info() to determine the platform numbers and
% their respective devices.
% NOTE: if you switch devices then you need to run the below line
% (uncommented) as well:
% clear mex
% NOTE: Using implementation 3 after using 2 might fail until MATLAB is
% restarted.
options.use_device = 0;

% Applies to implementations 2 and 3 ONLY
%%% Use 64-bit integer atomic functions
% If true, then 64-bit integer atomic functions (atomic add) will be used
% if they are supported by the selected device.
% Setting this to true will make computations faster on GPUs that support
% the functions, but might make results slightly less reliable due to
% floating point rounding. Recommended for GPUs.
options.use_64bit_atomics = true;

% Applies to implementations 2 and 3 ONLY
%%% Use 32-bit integer atomic functions
% If true, then 32-bit integer atomic functions (atomic add) will be used.
% This is even faster than the above 64-bit atomics version, but will also
% have significantly higher reduction in numerical/floating point accuracy.
% This should be about 20-30% faster than the above 64-bit version, but
% might lead to integer overflow if you have a high count measurement
% (thousands of coincidences per sinogram bin). Use this only if speed is
% of utmost importance. 64-bit atomics take precedence over 32-bit ones,
% i.e. if options.use_64bit_atomics = true then this will be always set as
% false.
options.use_32bit_atomics = false;

% Implementation 2 ONLY
%%% Use CUDA
% Selecting this to true will use CUDA kernels/code instead of OpenCL. This
% only works if the CUDA code was successfully built. Recommended only for
% Siddon as the orthogonal/volume-based ray tracer are slower in CUDA.
options.use_CUDA = false;

% Implementation 3 ONLY
%%% How many times more measurements/LORs are in the GPU part (applicable if
% heterogeneous computing (CPU + GPU) is used).
% Alternatively, set this to 0 to use only a single device on the specific
% platform in the multi-GPU or heterogenous case (the one with the highest
% memory count will be used). 
options.cpu_to_gpu_factor = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROJECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Type of projector to use for the geometric matrix
% 1 = Improved/accelerated Siddon's algorithm
% 3 = Volume of intersection based ray tracer
% See the wiki for more information:
% https://github.com/villekf/OMEGA/wiki/Useful-information#selecting-the-projector
options.projector_type = 1;

%%% Use point spread function (PSF) blurring
% Applies PSF blurring through convolution to the image space. This is the
% same as multiplying the geometric matrix with an image blurring matrix.
options.use_psf = false;

% FWHM of the Gaussian used in PSF blurring in all three dimensions
options.FWHM = [options.dPitch options.dPitch options.dPitch];

% Use deblurring phase
% If enabled, a deblurring phase is performed once the reconstruction has
% completed. This step is performed for all iterations (deblurred estimates
% are NOT used in the reconstruction phase). This is used ONLY when PSF
% blurring is used.
options.deblurring = false;
% Number of deblurring iterations
% How many iterations of the deblurring step is performed
options.deblur_iterations = 10;

% Orthogonal ray tracer (projector_type = 2) only
%%% The 2.5D (XY) width of the "strip/tube" where the orthogonal distances are
% included. If the tube_width_z is non-zero, then this value is ignored.
options.tube_width_xy = options.dPitch;

% Orthogonal ray tracer (projector_type = 2) only
%%% The 3D (Z) width of the "tube" where the orthogonal distances are
% included. If set to 0, then the 2.5D orthogonal ray tracer is used. If this
% value is non-zero then the above value is IGNORED.
options.tube_width_z = options.dPitch;

% Volume ray tracer (projector_type = 3) only
%%% Radius of the tube-of-response (cylinder)
% The radius of the cylinder that approximates the tube-of-response.
options.tube_radius = sqrt(2) * (options.dPitch / 2);

% Volume ray tracer (projector_type = 3) only
%%% Relative size of the voxel (sphere)
% In volume ray tracer, the voxels are modeled as spheres. This value
% specifies the relative radius of the sphere such that with 1 the sphere
% is just large enoough to encompass an entire cubic voxel, i.e. the
% corners of the cubic voxel intersect with the sphere shell. Larger values
% create larger spheres, while smaller values create smaller spheres.
options.voxel_radius = 1;

% Siddon (projector_type = 1) only
%%% Number of rays
% Number of rays used per detector if projector_type = 1 (i.e. Improved
% Siddon is used) and precompute_lor = false. I.e. when using precomputed
% LOR data, only 1 rays is always used.
% Number of rays in transaxial direction
options.n_rays_transaxial = 1;
% Number of rays in axial direction
options.n_rays_axial = 1;

% Orthogonal and volume ray tracers (projector_type = 2 and 3) only
% Implementations 2 and 3 ONLY
%%% Apply acceleration
% If true, then intermediate results are saved in memory. If you run out
% memory, you should set this to false. Does not apply to improved Siddon.
% Applies only if using either implementation 2 or 3. Can speed up
% computations by around 30% if set to true.
options.apply_acceleration = true;


%%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Number of iterations (all reconstruction methods)
options.Niter = 2;
% Save ALL iterations
% Set this to false if you do not want to save all the intermediate
% iterations, but only the very last one.
options.save_iter = false;

%%% Number of subsets (all excluding MLEM and subset_type = 5)
options.subsets = 16;

%%% Subset type (n = subsets)
% 1 = Every nth (column) measurement is taken
% 2 = Every nth (row) measurement is taken (e.g. if subsets = 3, then
% first subset has measurements 1, 4, 7, etc., second 2, 5, 8, etc.)
% 3 = Measurements are selected randomly
% 4 = (Sinogram only) Take every nth column in the sinogram
% 5 = (Sinogram only) Take every nth row in the sinogram
% 6 = Sort the LORs according to their angle with positive X-axis, combine
% n_angles together and have 180/n_angles subsets for 2D slices and
% 360/n_angles for 3D, see GitHub wiki for more information:
% https://github.com/villekf/OMEGA/wiki/Function-help#reconstruction-settings
% 7 = Form the subsets by using golden angle sampling
options.subset_type = 3;

%%% How many angles are combined in subset_type = 6
% E.g. there are 180 angles, in n_angles = 2, then angles 0 and 1 are
% combined to the same subset, 2 and 3, etc.
options.n_angles = 2;

%%% Initial value for the reconstruction
options.x0 = ones(options.Nx, options.Ny, options.Nz) * 1e-2;

%%% Epsilon value
% A small value to prevent division by zero and square root of zero. Should
% not be smaller than eps.
options.epps = 1e-8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MISC SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Use Shuffle (recommended)
% NOTE: Applies only when using subset_type = 3; accelerates the subset
% formation and uses less memory. Not included in OMEGA, needs to be
% manually downloaded and installed. Enabling this will automatically cause
% the function to be called if it is found on MATLAB path.
% Download from:
% https://se.mathworks.com/matlabcentral/fileexchange/27076-shuffle
options.use_Shuffle = false;

%%% Use fast sparse
% Not included in OMEGA, needs to be manually downloaded and installed.
% Download from: https://github.com/stefanengblom/stenglib
% NOTE: This applies only to implementation 1 when precompute_lor is false.
% Enabling this will automatically cause the function to be called if it is
% found on MATLAB path.
% NOTE: Suggested only for MATLAB 2019b and earlier.
options.use_fsparse = false;

%%% Skip the normalization phase in MRP, FMH, L-filter, ADMRP and
%%% weighted mean
% E.g. if set to true the MRP prior is (x - median(x))
% E.g. if set to false the MRP prior is (x - median(x)) / median(x)
% The published MRP uses the one that is obtained when this is set to
% false, however, you might get better results with true. I.e. this should
% be set to false if you wish to use the original prior implementation.
options.med_no_norm = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION ALGORITHMS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction algorithms to use (you can choose several)
% NOTE: MLEM requires precomputed observation matrix or a matrix-free
% method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ML-METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are non-regularized versions
%%% Maximum-Likelihood Expectation Maximization (MLEM)
% Supported by all implementations
% For implementation 1 requires the precomputed observation/system matrix
options.MLEM = false;

%%% Ordered Subsets Expectation Maximization (OSEM)
% Supported by all implementations
options.OSEM = true;

%%% Modified Row-Action Maximum Likelihood Algorithm (MRAMLA)
% Supported by implementations 1 and 2
options.MRAMLA = false;

%%% Row-Action Maximum Likelihood Algorithm (RAMLA)
% Supported by implementations 1, 2 and 4
options.RAMLA = false;

%%% Relaxed Ordered Subsets Expectation Maximization (ROSEM)
% Supported by implementations 1, 2 and 4
options.ROSEM = false;

%%% Rescaled Block Iterative Expectation Maximization (RBI-EM)
% Supported by implementations 1 and 2
options.RBI = false;

%%% Dynamic RAMLA (DRAMA)
% Supported by implementations 1, 2 and 4
options.DRAMA = false;

%%% Complete data OSEM (COSEM)
% Supported by implementations 1, 2 and 4
options.COSEM = false;

%%% Enhanced COSEM (ECOSEM)
% Supported by implementations 1, 2 and 4
options.ECOSEM = false;

%%% Accelerated COSEM (ACOSEM)
% Supported by implementations 1, 2 and 4
options.ACOSEM = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAP-METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Any algorithm selected here will utilize all the priors selected below.
% For example, if OSL-OSEM is set to true and MRP and Quad are set to true,
% then OSL-OSEM estimates will be computed for both MRP and Quadratic
% prior.
%%% One-Step Late MLEM (OSL-MLEM)
% Supported by implementations 2 and 4
options.OSL_MLEM = false;

%%% One-Step Late OSEM (OSL-OSEM)
% Supported by implementations 1, 2 and 4
options.OSL_OSEM = false;

%%% Modified BSREM (MBSREM)
% Supported by implementations 1 and 2
options.MBSREM = false;

%%% Block Sequential Regularized Expectation Maximation (BSREM)
% Supported by implementations 1, 2 and 4
options.BSREM = false;

%%% ROSEM-MAP
% Supported by implementations 1, 2 and 4
options.ROSEM_MAP = false;

%%% RBI-OSL
% Supported by implementations 1, 2 and 4
options.OSL_RBI = false;

%%% (A)COSEM-OSL
% 0/false = No COSEM-OSL, 1/true = ACOSEM-OSL, 2 = COSEM-OSL
% Supported by implementations 1, 2 and 4
options.OSL_COSEM = false;

%%% PKMA
% Supported by implementations 1, 2 and 4
options.PKMA = false;


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

%%% Total Generalized Variation (TGV) prior
options.TGV = false;

%%% Non-local Means (NLM) prior
options.NLM = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%% ACOSEM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Acceleration parameter for ACOSEM (1 equals COSEM)
options.h = 2;


%%%%%%%%%%%%%%%%%%%%%%%% MRAMLA & MBSREM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%
%%% Relaxation parameter for MRAMLA and MBSREM
% Use scalar if you want it to decrease as
% lambda0_mbsrem/current_iteration_number. Use vector (length = Niter) if
% you want your own relaxation parameters.
options.lambda0_MBSREM = 0.9;

%%% Upper bound for MRAMLA/MBSREM (use 0 for default (computed) value)
options.U = 0;


%%%%%%%%%%%%%%%%%%%%%%%%% RAMLA & BSREM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%
%%% Relaxation parameter for RAMLA and BSREM
% Use scalar if you want it to decrease as
% lambda0/current_iteration_number. Use vector (length = Niter) if you want
% your own relaxation parameters.
options.lambda0 = 0.2;


%%%%%%%%%%%%%%%%%%%%%%% ROSEM & ROSEM-MAP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%
%%% Relaxation parameter for ROSEM and ROSEM-MAP
% Use scalar if you want it to decrease as
% lambda0_rosem/current_iteration_number. Use vector (length = Niter) if
% you want your own relaxation parameters.
options.lambda0_ROSEM = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%% DRAMA PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Beta_0 value
options.beta0_drama = 0.1;
%%% Beta value
options.beta_drama = 1;
%%% Alpha value
options.alpha_drama = 0.1;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PKMA PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Relaxation parameter for PKMA
% If a scalar (or an empty) value is used, then the relaxation parameter is
% computed automatically as lambda(i) = 1 / ((i - 1)/12 + 1), where i is
% the iteration number. The input number thus has no effect.
% If, on the other hand, a vector is input then the input lambda values are
% used as is without any modifications (the length has to be at least the
% number of iterations).
options.lambda0_PKMA = 0;

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

%%% Additional step size (sigma) parameter for PKMA
% If a non-zero value is used, then the sigma parameter is
% computed automatically as sigma_PKMA = 1 - options.alpha_PKMA. If the
% input value is empty or zero, then this value is set as 1 for each subset
% and iteration, i.e. it has no effect on the reconstruction. 
% If, on the other hand, a vector is input then the input sigma values are
% used as is without any modifications (the length has to be at least the
% number of iterations * number of subsets).
% This value is based on the lambda value in
% https://doi.org/10.3846/1392-6292.2010.15.265-274
options.sigma_PKMA = 0;


%%%%%%%%%%%%%%%%%%%%%%%%% NEIGHBORHOOD PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%
%%% How many neighboring pixels are considered
% With MRP, QP, L, FMH, NLM and weighted mean
% E.g. if Ndx = 1, Ndy = 1, Ndz = 0, then you have 3x3 square area where
% the pixels are taken into account (I.e. (Ndx*2+1)x(Ndy*2+1)x(Ndz*2+1)
% area).
% NOTE: Currently Ndx and Ndy must be identical.
% For NLM this is often called the "search window".
options.Ndx = 1;
options.Ndy = 1;
options.Ndz = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MRP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for MRP with OSL-MLEM
options.beta_MRP_OSL_MLEM = 1.5;
%%% Regularization parameter for MRP with OSL-OSEM
options.beta_MRP_OSL_OSEM = 0.1;
%%% Regularization parameter for MRP with BSREM
options.beta_MRP_BSREM = 0.1;
%%% Regularization parameter for MRP with MBSREM
options.beta_MRP_MBSREM = 0.3;
%%% Regularization parameter for MRP with ROSEM
options.beta_MRP_ROSEM_MAP = 2;
%%% Regularization parameter for MRP with RBI
options.beta_MRP_OSL_RBI =  0.1;
%%% Regularization parameter for MRP with OSL-(A)COSEM
options.beta_MRP_OSL_COSEM = 1;
%%% Regularization parameter for MRP with PKMA
options.beta_MRP_PKMA =  0.1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for quadratic prior with OSL-OSEM
options.beta_quad_OSL_OSEM = 0.01;
%%% Regularization parameter for quadratic prior with OSL-MLEM
options.beta_quad_OSL_MLEM = 0.1;
%%% Regularization parameter for quadratic prior with MBSREM
options.beta_quad_MBSREM = 0.05;
%%% Regularization parameter for quadratic prior with BSREM
options.beta_quad_BSREM = 0.03;
%%% Regularization parameter for quadratic prior with ROSEM
options.beta_quad_ROSEM_MAP = 0.1;
%%% Regularization parameter for quadratic prior with RBI
options.beta_quad_OSL_RBI =  0.05;
%%% Regularization parameter for quadratic prior (OSL-(A)COSEM)
options.beta_quad_OSL_COSEM = 0.01;
%%% Regularization parameter for quadratic prior with PKMA
options.beta_quad_PKMA =  0.1;

%%% Pixel weights for quadratic prior
% The number of pixels need to be the amount of neighboring pixels,
% e.g. if the above Nd values are all 1, then 27 weights need to be
% included where the center pixel (if Nd values are 1, element 14) should
% be Inf. Size is (Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1). If left empty then
% they will be calculated by the algorithm and are based on the distance of
% the voxels from the center.
options.weights = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for Huber prior with OSL-OSEM
options.beta_Huber_OSL_OSEM = 0.01;
%%% Regularization parameter for Huber prior with OSL-MLEM
options.beta_Huber_OSL_MLEM = 0.1;
%%% Regularization parameter for Huber prior with MBSREM
options.beta_Huber_MBSREM = 0.05;
%%% Regularization parameter for Huber prior with BSREM
options.beta_Huber_BSREM = 0.03;
%%% Regularization parameter for Huber prior with ROSEM
options.beta_Huber_ROSEM_MAP = 0.1;
%%% Regularization parameter for Huber prior with RBI
options.beta_Huber_OSL_RBI =  0.05;
%%% Regularization parameter for Huber prior (OSL-(A)COSEM)
options.beta_Huber_OSL_COSEM = 0.01;
%%% Regularization parameter for Huber prior with PKMA
options.beta_Huber_PKMA =  0.1;

%%% Delta parameter for Huber prior
% Upper and lower bounds for the prior
options.huber_delta = 5;

%%% Pixel weights for Huber prior
% Same rules apply as with quadratic prior weights.
% If left empty then they will be calculated by the algorithm and are based
% on the distance of the voxels from the center.
options.weights_huber = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%% L-FILTER PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for L-filter with OSL-OSEM
options.beta_L_OSL_OSEM = 0.1;
%%% Regularization parameter for L-filter with OSL-MLEM
options.beta_L_OSL_MLEM = 0.1;
%%% Regularization parameter for L-filter with MBSREM
options.beta_L_MBSREM = 0.1;
%%% Regularization parameter for L-filter with BSREM
options.beta_L_BSREM = 0.03;
%%% Regularization parameter for L-filter with ROSEM
options.beta_L_ROSEM_MAP = 3;
%%% Regularization parameter for L-filter with RBI
options.beta_L_OSL_RBI =  0.09;
%%% Regularization parameter for L-filter (OSL-(A)COSEM)
options.beta_L_OSL_COSEM = 0.1;
%%% Regularization parameter for L-filter with PKMA
options.beta_L_PKMA =  0.1;

%%% Weighting factors for the L-filter pixels
% Otherwise the same as in quadratic prior, but center pixel is not Inf.
% If left empty then they will be calculated by the algorithm such that the
% weights resemble a Laplace distribution.
options.a_L = [];

%%% If the weighting factors are set empty, then this option will determine
% whether the computed weights follow a 1D weighting scheme (true) or 2D
% (false).
% See the wiki for more information:
% https://github.com/villekf/OMEGA/wiki/Function-help#reconstruction-algorithms
options.oneD_weights = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FMH PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for FMH with OSL-OSEM
options.beta_FMH_OSL_OSEM = 0.1;
%%% Regularization parameter for FMH with OSL-MLEM
options.beta_FMH_OSL_MLEM = 0.1;
%%% Regularization parameter for FMH with MBSREM
options.beta_FMH_MBSREM = 0.6;
%%% Regularization parameter for FMH with BSREM
options.beta_FMH_BSREM = 5;
%%% Regularization parameter for FMH with ROSEM
options.beta_FMH_ROSEM_MAP = 8;
%%% Regularization parameter for FMH with RBI
options.beta_FMH_OSL_RBI =  0.5;
%%% Regularization parameter for FMH (OSL-(A)COSEM)
options.beta_FMH_OSL_COSEM = 0.1;
%%% Regularization parameter for FMH with PKMA
options.beta_FMH_PKMA =  0.1;

%%% Pixel weights for FMH
% The matrix size needs to be [Ndx*2+1, 4] if Nz = 1 or Ndz = 0, or
% [Ndx*2+1, 13] otherwise.
% The center pixel weight should be in the middle of the weight matrix.
% If the sum of each column is > 1, then the weights will be normalized
% such that the sum = 1.
% If left empty then they will be calculated by the algorithm such that the
% weights follow the same pattern as in the original article.
options.fmh_weights = [];

%%% Weighting value for the center pixel
% Default value is 4, which was used in the original article.
% NOTE: This option is ignored if you provide your own weights.
options.fmh_center_weight = 4;


%%%%%%%%%%%%%%%%%%%%%%%%% WEIGHTED MEAN PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%
%%% Mean type
% 1 = Arithmetic mean, 2 = Harmonic mean, 3 = Geometric mean
options.mean_type = 1;

%%% Regularization parameter for weighted mean with OSL-OSEM
options.beta_weighted_mean_OSL_OSEM = 0.2;
%%% Regularization parameter for weighted mean with OSL-MLEM
options.beta_weighted_mean_OSL_MLEM = 0.1;
%%% Regularization parameter for weighted mean with MBSREM
options.beta_weighted_mean_MBSREM = 0.1;
%%% Regularization parameter for weighted mean with BSREM
options.beta_weighted_mean_BSREM = 5;
%%% Regularization parameter for weighted mean with ROSEM
options.beta_weighted_mean_ROSEM_MAP = 3;
%%% Regularization parameter for weighted mean with RBI
options.beta_weighted_mean_OSL_RBI =  0.04;
%%% Regularization parameter for weighted mean (OSL-(A)COSEM)
options.beta_weighted_mean_OSL_COSEM = 0.2;
%%% Regularization parameter for weighted mean with PKMA
options.beta_weighted_mean_PKMA =  0.1;

%%% Pixel weights for weighted mean
% The number of pixels needs to be the amount of neighboring pixels,
% e.g. if the above Ndx/y/z values are all 1, then 27 weights need to be
% included. Size is (Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1). If left empty then
% they will be calculated by the algorithm such that the weights are
% dependent on the distance from the center pixel to the neighboring
% pixels.
options.weighted_weights = [];

%%% Center pixel weight for weighted mean.
% NOTE: This option is ignored if you provide your own weights.
options.weighted_center_weight = 4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TV PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for TV OSL-OSEM
options.beta_TV_OSL_OSEM = 0.5;
%%% Regularization parameter for TV with OSL-MLEM
options.beta_TV_OSL_MLEM = 0.1;
%%% Regularization parameter for TV with MBSREM
options.beta_TV_MBSREM = 0.01;
%%% Regularization parameter for TV with BSREM
options.beta_TV_BSREM = 0.05;
%%% Regularization parameter for TV with ROSEM
options.beta_TV_ROSEM_MAP = 0.07;
%%% Regularization parameter for TV with RBI
options.beta_TV_OSL_RBI =  0.002;
%%% Regularization parameter for TV (OSL-(A)COSEM)
options.beta_TV_OSL_COSEM = 0.003;
%%% Regularization parameter for TV with PKMA
options.beta_TV_PKMA =  0.1;

%%% "Smoothing" parameter
% Also used to prevent zero values in square root.
options.TVsmoothing = 1e-2;

%%% Whether to use an anatomical reference/weighting image with the TV
options.TV_use_anatomical = false;

%%% If the TV_use_anatomical value is set to true, specify filename for the
% reference image here (same rules apply as with attenuation correction
% above).
options.TV_reference_image = 'reference_image.mat';

%%% Three different TV methods are available.
% Value can be 1, 2 or 3.
% Types 1 and 2 are the same if anatomical prior is not included
% Type 3 uses the same weights as quadratic prior
% See the wiki for more information:
% https://github.com/villekf/OMEGA/wiki/Function-help#reconstruction-algorithms
options.TVtype = 1;

%%% Weighting parameters for the TV prior.
% Applicable only if use_anatomical = true. T-value is specific to the used
% TVtype, e.g. for type 1 it is the edge threshold parameter. See the wiki
% for more details:
% https://github.com/villekf/OMEGA/wiki/Function-help#reconstruction-algorithms
options.T = 0.5;

%%% C is the weight for the original image in type 3 and is ignored with
% other types
options.C = 1;

%%% Tuning parameter for TV and APLS
options.tau = 1e-8;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADMRP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for AD with OSL-OSEM
options.beta_AD_OSL_OSEM = 0.1;
%%% Regularization parameter for AD with OSL-MLEM
options.beta_AD_OSL_MLEM = 0.1;
%%% Regularization parameter for AD with MBSREM
options.beta_AD_MBSREM = 0.3;
%%% Regularization parameter for AD with BSREM
options.beta_AD_BSREM = 0.2;
%%% Regularization parameter for AD with ROSEM
options.beta_AD_ROSEM_MAP = 0.0003;
%%% Regularization parameter for AD with RBI
options.beta_AD_OSL_RBI =  0.05;
%%% Regularization parameter for AD with (OSL-(A)COSEM)
options.beta_AD_OSL_COSEM = 0.1;
%%% Regularization parameter for AD with PKMA
options.beta_AD_PKMA =  0.1;

%%% Time step variable for AD (implementation 2 only)
options.TimeStepAD = 0.0625;

%%% Conductivity/connectivity for AD (edge threshold)
options.KAD = 2;

%%% Number of iterations for AD filter
% NOTE: This refers to the AD smoothing part, not the actual reconstruction
% phase.
options.NiterAD = 1;

%%% Flux/conduction type for AD filter
% 1 = Exponential
% 2 = Quadratic
options.FluxType = 1;

%%% Diffusion type for AD (implementation 2 only)
% 1 = Gradient
% 2 = Modified curvature
options.DiffusionType = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% APLS PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for APLS with OSL-OSEM
options.beta_APLS_OSL_OSEM = 0.01;
%%% Regularization parameter for APLS with OSL-MLEM
options.beta_APLS_OSL_MLEM = 0.1;
%%% Regularization parameter for APLS with MBSREM
options.beta_APLS_MBSREM = 0.1;
%%% Regularization parameter for APLS with BSREM
options.beta_APLS_BSREM = 0.005;
%%% Regularization parameter for APLS with ROSEM
options.beta_APLS_ROSEM_MAP = 0.1;
%%% Regularization parameter for APLS with RBI
options.beta_APLS_OSL_RBI =  0.1;
%%% Regularization parameter for APLS (OSL-(A)COSEM)
options.beta_APLS_OSL_COSEM = 0.01;
%%% Regularization parameter for APLS with PKMA
options.beta_APLS_PKMA =  0.1;

%%% Scaling parameter (eta)
% See the wiki for details:
% https://github.com/villekf/OMEGA/wiki/Function-help#reconstruction-algorithms
options.eta = 1e-5;

%%% "Smoothing" parameter (beta)
% Also used to prevent zero values in square root.
options.APLSsmoothing = 1e-5;

%%% Specify filename for the reference image here (same rules apply as with
% attenuation correction above)
% NOTE: For APSL, the prior is required.
options.APLS_reference_image = 'reference_image.mat';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TGV PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for TGV with OSL-OSEM
options.beta_TGV_OSL_OSEM = 0.1;
%%% Regularization parameter for TGV with OSL-MLEM
options.beta_TGV_OSL_MLEM = 0.1;
%%% Regularization parameter for TGV with MBSREM
options.beta_TGV_MBSREM = 0.1;
%%% Regularization parameter for TGV with BSREM
options.beta_TGV_BSREM = 1;
%%% Regularization parameter for TGV with ROSEM
options.beta_TGV_ROSEM_MAP = 0.25;
%%% Regularization parameter for TGV with RBI
options.beta_TGV_OSL_RBI =  0.1;
%%% Regularization parameter for TGV (OSL-(A)COSEM)
options.beta_TGV_OSL_COSEM = 0.05;
%%% Regularization parameter for TGV with PKMA
options.beta_TGV_PKMA =  0.1;

%%% TGV weights
% First part
options.betaTGV = 1;
% Second part (symmetrized derivative)
options.alphaTGV = 2;

%%% Number of TGV iterations
options.NiterTGV = 30;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NLM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for NLM with OSL-OSEM
options.beta_NLM_OSL_OSEM = 0.01;
%%% Regularization parameter for NLM with OSL-MLEM
options.beta_NLM_OSL_MLEM = 0.01;
%%% Regularization parameter for NLM with MBSREM
options.beta_NLM_MBSREM = 0.05;
%%% Regularization parameter for NLM with BSREM
options.beta_NLM_BSREM = 0.01;
%%% Regularization parameter for NLM with ROSEM
options.beta_NLM_ROSEM_MAP = 0.1;
%%% Regularization parameter for NLM with RBI
options.beta_NLM_OSL_RBI =  0.01;
%%% Regularization parameter for NLM (OSL-(A)COSEM)
options.beta_NLM_OSL_COSEM = 0.01;
%%% Regularization parameter for NLM with PKMA
options.beta_NLM_PKMA =  0.01;

%%% Filter parameter
options.sigma = 10;

%%% Patch radius
options.Nlx = 9;
options.Nly = 9;
options.Nlz = 0;

%%% Standard deviation of the Gaussian filter
options.NLM_gauss = 1;

% Search window radius is controlled by Ndx, Ndy and Ndz parameters
% Use anatomical reference image for the patches
options.NLM_use_anatomical = false;

%%% Specify filename for the reference image here (same rules apply as with
% attenuation correction above)
options.NLM_reference_image = 'reference_image.mat';

%%% Use Non-local total variation (NLTV)
% If selected, will overwrite regular NLM regularization as well as the
% below MRP version.
options.NLTV = false;

%%% Use MRP algorithm (without normalization)
% I.e. gradient = im - NLM_filtered(im)
options.NLM_MRP = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% OPENCL DEVICE INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Implementation 2
% Uncomment the below line and run it to determine the available device
% numbers
% ArrayFire_OpenCL_device_info()

%%% Implementation 3
% Uncomment the below line and run it to determine the available platforms,
% their respective numbers and device numbers
% OpenCL_device_info()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% 2D (sinogram) reconstruction can be enabled with the following changes:
% options.SinM = squeeze(sum(options.SinM,2));
% options.xSize = 1;
% options.axial_fov = options.dPitch;
% options.Nz = 1;
% options.x0 = ones(options.Nx, options.Ny, options.Nz) * 1e-2;

%%

pz = reconstructions_mainCT(options);

%% MLEM example and class constructor example

raw_SinM = options.SinM(:);
options.subsets = 1;

% There are multiple ways to construct the class object.
% The recommended way is using the options-struct created by the first
% section.
A = forwardBackwardProjectCT(options);
% However, it is also possible to create the class object by inputting the
% detector coordinates as well as image and FOV sizes. However, currently
% it is not possible to select the projector when using this method
% (improved Siddon is always used) and PSF blurring is not available.
% [x,y,z] = CTDetectorCoordinatesFull(-options.angles,options.sourceToDetector,options.sourceToCRot,options.dPitch,options.xSize,...
%     options.ySize);
% A = forwardBackwardProjectCT(x,y,z,options.Nx,options.Ny,options.Nz,options.FOVa_x,options.axial_fov, options.subsets, options.implementation, options.use_device);
% Lastly, it is also possible to create the class object by using the
% projection angles. Same restrictions apply as when using detector
% coordinates.
% A = forwardBackwardProjectCT(options.angles, options.ySize, options.xSize,options.Nx,options.Ny,options.Nz,options.FOVa_x,options.axial_fov, ...
%     options.dPitch, options.sourceToDetector, options.sourceToCRot, 0, 0, options.subsets, options.implementation, options.use_device);

f = options.x0(:);
for iter = 1 : 5
    % The result is stored in y
    y = exp(-(A*f));
    % The result is stored in x
    if iter == 1
        % Sensitivity image can be computed during the first iteration
        % Computed ONLY if the second output is present
        [x, A] = backwardProject(A, y ./ raw_SinM);
    else
        x = A' * (y ./ raw_SinM);
    end
    f = (f ./ (A.sens + options.epps)) .* (x + options.epps);
end
% PSF deblurring phase
if options.use_psf && options.deblurring
    f = deblur(f, options, gaussK, options.Nx, options.Ny, options.Nz);
end
ff = reshape(f, options.Nx,options.Ny,options.Nz);

% Any variable named f_osem will be automatically visualized (as OSEM) by
% visualize_pet.m as long as the pz-variable does not exist (and OSEM has
% been selected in visualize_pet.m)
f_osem = ff;

clear pz

%% OSEM (subset) example

raw_SinM = options.SinM(:);
options.subsets = 16;
A = forwardBackwardProjectCT(options);
% [x,y,z] = CTDetectorCoordinatesFull(-options.angles,options.sourceToDetector,options.sourceToCRot,options.dPitch,options.xSize,...
%     options.ySize);
% A = forwardBackwardProjectCT(x,y,z,options.Nx,options.Ny,options.Nz,options.FOVa_x,options.axial_fov, options.subsets, options.implementation, options.use_device);
% clear x y z
% A = forwardBackwardProjectCT(options.angles, options.ySize, options.xSize,options.Nx,options.Ny,options.Nz,options.FOVa_x,options.axial_fov, ...
%     options.dPitch, options.sourceToDetector, options.sourceToCRot, options.horizontalOffset, 0, options.subsets, options.implementation, options.use_device);

f = ones(options.Nx * options.Ny * options.Nz,1) * 1e-2;
% Default uses random subset sampling
raw_SinM = single(raw_SinM(A.index));
for iter = 1 : 2
    for osa_iter = 1 : options.subsets
        % You can input the subset to the class object either with
        A.subset = osa_iter;
        % Or you can simply call the functions forwardProject and/or
        % backwardProject with the subset index
        % y = exp(-forwardProject(A, f, osa_iter));
        % If the first method is used, you can simply use multiplication
        y = exp(-(A*f));
        % The result is stored in x
        if iter == 1
            % Sensitivity image can be computed during the first iteration
            % Computed ONLY if the second output is present
            [x, A] = backwardProject(A, y ./ raw_SinM(A.nn(osa_iter) + 1:A.nn(osa_iter+1)),osa_iter);
        else
            % x = backwardProject(A, y ./ raw_SinM(A.nn(osa_iter) + 1:A.nn(osa_iter+1)), osa_iter);
            x = A' * (y ./ raw_SinM(A.nn(osa_iter) + 1:A.nn(osa_iter+1)));
        end
        f = (f ./ (A.sens(:,osa_iter) + 1e-8)) .* (x + 1e-8);
    end
end
% PSF deblurring phase
if options.use_psf && options.deblurring
    f = deblur(f, options, gaussK, options.Nx, options.Ny, options.Nz);
end
ff = reshape(f, options.Nx,options.Ny,options.Nz);

f_osem = ff;

clear pz



%% Chambolle-Pock example with or without subsets

raw_SinM = log(1./options.SinM(:));
options.subsets = 1;
A = forwardBackwardProjectCT(options);

% This is computed with 
% tau = 1/powerMethod(A, 10, false);
tau = 0.0059534;
alpha = 0.01;
Niter = 50;
sigma = tau;
sigma2 = sigma;
theta = 1;

N = [options.Nx options.Ny options.Nz];

% Preconditioners
% sigma = 1 ./ (A * ones(prod(N),1));
% tau = 1 ./ (A' * ones(size(raw_SinM)));
% sigma2 = 1 ./ alpha;

if options.subsets > 1
    raw_SinM = raw_SinM(A.index);
end
M = size(raw_SinM,1);
Ndim = N;
sDims = length(Ndim); % spatial dims
N = prod(Ndim);

% Primal variables
u = zeros(N,1);
u2 = zeros(N,1);

% Data dual term
p = zeros(M,1);

% Spatial reg help variable
qHelp = zeros(sDims*N,1);
q = qHelp;
qHelp2 = zeros(N,1);

gradx = zeros(Ndim);
grady = zeros(Ndim);
gradz = zeros(Ndim);
divx = zeros(Ndim);
divy = zeros(Ndim);
divz = zeros(Ndim);

res = -raw_SinM;
%%
tic
for ii = 1:Niter
    for osa_iter = 1 : options.subsets
        uPrev = u;
        
        if options.subsets > 1
            A.subset = osa_iter;
        end
        temp2 = A * u2(:);
        res((A.nn(osa_iter) + 1:A.nn(osa_iter+1))) = (temp2)-raw_SinM((A.nn(osa_iter) + 1:A.nn(osa_iter+1)));
        % Dual variable updates
        p(A.nn(osa_iter) + 1:A.nn(osa_iter+1)) = (p(A.nn(osa_iter) + 1:A.nn(osa_iter+1)) + sigma*res(A.nn(osa_iter) + 1:A.nn(osa_iter+1)))/(1 + sigma);
        % KL version
%         p(A.nn(osa_iter) + 1:A.nn(osa_iter+1)) = 0.5 .* (ones(size(p(A.nn(osa_iter) + 1:A.nn(osa_iter+1)))) + p(A.nn(osa_iter) + 1:A.nn(osa_iter+1))...
%             + temp2 .* sigma - sqrt((p(A.nn(osa_iter) + 1:A.nn(osa_iter+1)) + temp2 .* sigma - ones(size(p(A.nn(osa_iter) + 1:A.nn(osa_iter+1))))).^2 ...
%             + 4 .* sigma .* raw_SinM));
        
        qHelp2 = abs(qHelp).^2;
        q = qHelp./max(1,sqrt(qHelp2 + circshift(qHelp2,N) + circshift(qHelp2,2*N))/alpha);
        
        % u update step
        uUpd1 = A'*(p(A.nn(osa_iter) + 1:A.nn(osa_iter+1)));
        qx = reshape(q(1:N),Ndim);
        qy = reshape(q(N+1:2*N),Ndim);
        divx(:,2:Ndim(2),:) = -diff(qx,1,2);
        divx(:,1,:) = -qx(:,1,:);
        divx(:,Ndim(2),:) = qx(:,Ndim(2)-1,:);
        divy(2:Ndim(1),:,:) = -diff(qy,1,1);
        divy(1,:,:) = -qy(1,:,:);
        divy(Ndim(1),:,:) = qy(Ndim(1)-1,:,:);
        if sDims >= 3
            qz = reshape(q(N*2+1:end),Ndim);
            divz(:,:,2:Ndim(3)) = -diff(qz,1,3);
            divz(:,:,1) = -qz(:,:,1);
            divz(:,:,Ndim(3)) = qz(:,:,Ndim(3)-1);
            uUpd2 = (divx + divy + divz);
        else
            uUpd2 = (divx + divy);
        end
        uUpd = uUpd1 + uUpd2(:);
        
        % Update u and u2
        
        u = u - tau .* uUpd;
        u(u < 0) = 0;
        u2 = u + theta * (u - uPrev);
        
        % Update help variables
        u2temp = reshape(u2(:),Ndim);
        gradx(:,1:Ndim(2)-1,:) = diff(u2temp,1,2);
        grady(1:Ndim(1)-1,:,:) = diff(u2temp,1,1);
        if sDims >= 3
            gradz(:,:,1:Ndim(3)-1) = diff(u2temp,1,3);
            tvPart = [gradx(:);grady(:);gradz(:)];
        else
            tvPart = [gradx(:);grady(:)];
        end
        qHelp = q(:) + sigma2 .* tvPart(:);
        
    end
    
end
toc
ff = reshape(u, options.Nx,options.Ny,options.Nz);

f_osem = ff;

clear pz

%% Matrix example using OSEM/MLEM

% For non-subset algorithm, use 1
options.subsets = 16;
% For matrix version, implementation 1 MUST be used
options.implementation = 1;
C = forwardBackwardProjectCT(options);
raw_SinM = options.SinM(:);

f = options.x0(:);
% If subsets are used
raw_SinM = double(raw_SinM(A.index));
% Form the full system matrix with this:
% B = formMatrix(C);
for iter = 1 : 5
    for osa_iter = 1 : options.subsets
        % Form a subset of the system matrix with this:
        B = formMatrix(C, osa_iter);
        
        if options.use_psf
            f = computeConvolution(f, options, options.Nx, options.Ny, options.Nz, gaussK);
        end
        
        y = exp(-(B' * f));
        Summ = full(sum(B,2));
        if options.use_psf
            Summ = computeConvolution(Summ, options, options.Nx, options.Ny, options.Nz, gaussK);
        end
        x = B * (y ./ raw_SinM(C.nn(osa_iter) + 1:C.nn(osa_iter+1)));
        if options.use_psf
            x = computeConvolution(x, options, options.Nx, options.Ny, options.Nz, gaussK);
        end
        f = (f ./ (Summ + options.epps)) .* (x + options.epps);
        % For subset version it is advised to clear the previous matrix
        clear B
    end
end
ff = reshape(f, options.Nx,options.Ny,options.Nz);

f_osem = ff;

clear pz