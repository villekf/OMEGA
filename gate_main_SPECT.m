%% MATLAB codes for PET reconstruction using sinogram input from any scanner

%clear
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SCANNER PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Number of detector heads
options.nHeads = 2;

%%% Number of projections (for all heads, divided by the number of heads)
options.nProjections = 120;

%%% Starting angle for each head
% Use multidimensional array for multiple heads
% Both radians and degrees are supported
options.startAngle = [0;180];

%%% Angle increment
options.angleIncrement = 6;

%%% Radius from the center of rotation to the detector face
% The radius for each projection
options.radiusPerProj = repmat(230, options.nProjections * options.nHeads, 1);

%%% Crystal pitch/size in x- and y-directions (transaxial) (mm)
options.cr_p = 0.904;

%%% Transaxial FOV size (mm), this is the length of the x (horizontal) side
% of the FOV
options.FOVa_x = 0.904*128;

%%% Transaxial FOV size (mm), this is the length of the y (vertical) side
% of the FOV
options.FOVa_y = options.FOVa_x;

%%% Axial FOV (mm)
options.axial_fov = 0.904*128;

%%% Scanner name
% Used for naming purposes (measurement data)
options.machine_name = 'Two_Heads_SPECT_example';
 
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
 
%%% Reconstructed image pixel count (X-direction)
% NOTE: Non-square image sizes (X- and Y-direction) may not work
options.Nx = 128;

%%% Y-direction
options.Ny = 128;

%%% Z-direction (number of slices) (axial)
options.Nz = 128;

%%% Flip the image (in vertical direction)?
options.flip_image = false;

%%% How much is the image rotated?
% You need to run the precompute phase again if you modify this
% NOTE: The rotation is done in the detector space (before reconstruction).
% This current setting is for systems whose detector blocks start from the
% right hand side when viewing the device from front.
% Positive values perform the rotation in clockwise direction
options.offangle = 0;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% COLLIMATOR PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Collimator-detector response function (CDRF)
% You can either input the (Gaussian) PSF filter, or the standard
% deviations for both transaxial and axial directions or simply the
% collimator parameters (see below) for an analytic solution for round (and
% hexagonal) holes (this may be unoptimal).

% If you want to compute the CDRF analytically, input the following values:
% Collimator hole length (mm)
options.colL = 30;
% Collimator hole radius
options.colR = 1.5;
% Distance from collimator to the detector
options.colD = 0;
% Intrinsic resolution
options.iR = 3;

% If you have the standard deviations for transaxial (XY) and axial (Z)
% directions, you can input them here instead of the above values (the
% dimensions need to be options.nProjections x options.Nx):
% Transaxial standard deviation
% options.sigmaXY = repmat(0, options.nProjection, options.Nx);
% Axial standard deviation
% options.sigmaZ = repmat(0, options.nProjection, options.Nx);

% Lastly, you can input the filter for the CDRF directly. This should be
% filterSizeXY x filterSizeZ x options.nProjections:
% options.gFilter = ones(10,10,options.nProjections);
 
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
% This is used to name the saved measurement data and also load it in
% future sessions.
options.name = 'spect_example';

%%% Show status messages
% These are e.g. time elapsed on various functions and what steps have been
% completed. It is recommended to keep this true.
options.verbose = true;
 
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
% 2 = Matrix-free reconstruction with OpenCL/ArrayFire (Recommended)
% (Requires ArrayFire. Compiles with MinGW ONLY when ArrayFire was compiled
% with MinGW as well (cannot use the prebuilt binaries)).
options.implementation = 2;

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

% Implementation 2 ONLY
%%% Use CUDA
% Selecting this to true will use CUDA kernels/code instead of OpenCL. This
% only works if the CUDA code was successfully built. Recommended only for
% Siddon as the orthogonal/volume-based ray tracer are slower in CUDA.
options.use_CUDA = false;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROJECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Type of projector to use for the geometric matrix
% 6 = Rotation-based projector
options.projector_type = 6;
 
%%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Number of iterations (all reconstruction methods)
options.Niter = 40;
% Save ALL iterations
% Set this to false if you do not want to save all the intermediate
% iterations, but only the very last one.
options.save_iter = false;

%%% Number of subsets (all excluding MLEM and subset_type = 6)
options.subsets = 8;

%%% Subset type (n = subsets)
% 8 = Every nth projection is taken
% 9 = The projection images are selected randomly
options.subset_type = 8;

%%% Initial value for the reconstruction
options.x0 = ones(options.Nx, options.Ny, options.Nz);

%%% Epsilon value 
% A small value to prevent division by zero and square root of zero. Should
% not be smaller than eps.
options.epps = 1e-8;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MISC SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Use Shuffle (recommended)
% NOTE: Applies only when using subset_type = 3. 
% Accelerates the subset formation and uses less memory. Not included in
% OMEGA, needs to be manually downloaded and installed.
% Download from: 
% https://se.mathworks.com/matlabcentral/fileexchange/27076-shuffle
options.use_Shuffle = false;

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
% Supported by implementations 1, 2 and 4
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
% Any algorithm selected here will utilize all the priors selected below
% this. For example, if OSL-OSEM is set to true and MRP and Quad are set to
% true, then OSL-OSEM estimates will be computed for both MRP and Quadratic
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

%%% Block Sequential Regularized Expectation Maximization (BSREM)
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

%%% Relative Difference Prior (RDP)
options.RDP = false;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% ACOSEM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Acceleration parameter for ACOSEM (1 equals COSEM)
options.h = 2;
 

%%%%%%%%%%%%%%%%%%%%%%%% MRAMLA & MBSREM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%
%%% Relaxation parameter for MRAMLA and MBSREM
% Use scalar if you want it to decrease as
% lambda0_mbsrem/current_iteration_number. Use vector (length = Niter) if
% you want your own relaxation parameters.
options.lambda0_MBSREM = 0.2;

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
% options.delta_PKMA), where i is the iteration number and ll the subset
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
options.Ndx = 2;
options.Ndy = 2;
options.Ndz = 2;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MRP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for MRP with OSL-OSEM
options.beta_MRP_OSL_OSEM = 0.005;
%%% Regularization parameter for MRP with OSL-MLEM
options.beta_MRP_OSL_MLEM = 1.5;
%%% Regularization parameter for MRP with MBSREM
options.beta_MRP_MBSREM = 0.3;
%%% Regularization parameter for MRP with BSREM
options.beta_MRP_BSREM = 0.1;
%%% Regularization parameter for MRP with ROSEM
options.beta_MRP_ROSEM_MAP = 2;
%%% Regularization parameter for MRP with RBI
options.beta_MRP_OSL_RBI = 0.1;
%%% Regularization parameter for MRP with OSL-(A)COSEM
options.beta_MRP_OSL_COSEM = 1;
%%% Regularization parameter for MRP with PKMA
options.beta_MRP_PKMA =  0.001;
 
 
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
options.beta_quad_OSL_RBI = 0.05;
%%% Regularization parameter for quadratic prior (OSL-(A)COSEM)
options.beta_quad_OSL_COSEM = 0.01;
%%% Regularization parameter for quadratic prior with PKMA
options.beta_quad_PKMA =  0.1;

%%% Pixel weights for quadratic prior
% The number of pixels need to be the amount of neighboring pixels,
% e.g. if the above Nd values are all 1, then 27 weights need to be
% included where the center pixel (if Nd values are 1, element 14) should
% be Inf. Size is (options.Ndx*2+1) * (options.Ndy*2+1) * (options.Ndz*2+1). If left empty then
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
options.beta_Huber_OSL_RBI = 0.05;
%%% Regularization parameter for Huber prior (OSL-(A)COSEM)
options.beta_Huber_OSL_COSEM = 0.01;
%%% Regularization parameter for Huber prior with PKMA
options.beta_Huber_PKMA = 0.1;

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
options.beta_L_OSL_RBI = 0.09;
%%% Regularization parameter for L-filter (OSL-(A)COSEM)
options.beta_L_OSL_COSEM = 0.1;
%%% Regularization parameter for L-filter with PKMA
options.beta_L_PKMA = 0.1;

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
options.beta_FMH_OSL_RBI = 0.5;
%%% Regularization parameter for FMH (OSL-(A)COSEM)
options.beta_FMH_OSL_COSEM = 0.1;
%%% Regularization parameter for FMH with PKMA
options.beta_FMH_PKMA = 0.1;

%%% Pixel weights for FMH
% The matrix size needs to be [options.Ndx*2+1, 4] if Nz = 1 or options.Ndz = 0, or
% [options.Ndx*2+1, 13] otherwise.
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
options.beta_weighted_mean_OSL_RBI = 0.04;
%%% Regularization parameter for weighted mean (OSL-(A)COSEM)
options.beta_weighted_mean_OSL_COSEM = 0.2;
%%% Regularization parameter for weighted mean with PKMA
options.beta_weighted_mean_PKMA = 0.1;

%%% Pixel weights for weighted mean
% The number of pixels needs to be the amount of neighboring pixels,
% e.g. if the above options.Ndx/y/z values are all 1, then 27 weights need to be
% included. Size is (options.Ndx*2+1) * (options.Ndy*2+1) * (options.Ndz*2+1). If left empty then
% they will be calculated by the algorithm such that the weights are
% dependent on the distance from the center pixel to the neighboring
% pixels.
options.weighted_weights = [];

%%% Center pixel weight for weighted mean.
% NOTE: This option is ignored if you provide your own weights.
options.weighted_center_weight = 4;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TV PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for TV OSL-OSEM
options.beta_TV_OSL_OSEM = 0.0005;
%%% Regularization parameter for TV with OSL-MLEM
options.beta_TV_OSL_MLEM = 0.1;
%%% Regularization parameter for TV with MBSREM
options.beta_TV_MBSREM = 0.01;
%%% Regularization parameter for TV with BSREM
options.beta_TV_BSREM = 0.05;
%%% Regularization parameter for TV with ROSEM
options.beta_TV_ROSEM_MAP = 0.07;
%%% Regularization parameter for TV with RBI
options.beta_TV_OSL_RBI = 0.002;
%%% Regularization parameter for TV (OSL-(A)COSEM)
options.beta_TV_OSL_COSEM = 0.003;
%%% Regularization parameter for TV with PKMA
options.beta_TV_PKMA = 0.002;

%%% "Smoothing" parameter
% Also used to prevent zero values in square root.
options.TVsmoothing = 1e-1;

%%% Whether to use an anatomical reference/weighting image with the TV
options.TV_use_anatomical = false;

%%% If the TV_use_anatomical value is set to true, specify filename for the
% reference image here (same rules apply as with attenuation correction
% above).
options.TV_reference_image = 'reference_image.mat';

%%% Three different TV methods are available.
% Value can be 1, 2, 3 or 4.
% Types 1 and 2 are the same if anatomical prior is not included
% Type 3 uses the same weights as quadratic prior
% Type 4 is smoothed anisotropic TV, does not support anatomic weighting,
% can also be used as a regular anisotropic TV
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

%%% Tuning parameter for Lange function in SATV
% Setting this to 0 gives regular anisotropic TV
options.SATVPhi = 0.2;
 
 
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
options.beta_AD_OSL_RBI = 0.05;
%%% Regularization parameter for AD with (OSL-(A)COSEM)
options.beta_AD_OSL_COSEM = 0.1;
%%% Regularization parameter for AD with PKMA
options.beta_AD_PKMA = 0.1;

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
options.beta_APLS_OSL_RBI = 0.1;
%%% Regularization parameter for APLS (OSL-(A)COSEM)
options.beta_APLS_OSL_COSEM = 0.01;
%%% Regularization parameter for APLS with PKMA
options.beta_APLS_PKMA = 0.1;

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
options.beta_TGV_OSL_RBI = 0.1;
%%% Regularization parameter for TGV (OSL-(A)COSEM)
options.beta_TGV_OSL_COSEM = 0.05;
%%% Regularization parameter for TGV with PKMA
options.beta_TGV_PKMA = 0.1;

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
%%% Regularization parameter for NLM with MBSREM
options.beta_NLM_MBSREM = 0.05;
%%% Regularization parameter for NLM with BSREM
options.beta_NLM_BSREM = 0.01;
%%% Regularization parameter for NLM with ROSEM
options.beta_NLM_ROSEM_MAP = 0.1;
%%% Regularization parameter for NLM with RBI
options.beta_NLM_OSL_RBI = 0.01;
%%% Regularization parameter for NLM (OSL-(A)COSEM)
options.beta_NLM_OSL_COSEM = 0.01;
%%% Regularization parameter for NLM with PKMA
options.beta_NLM_PKMA = 0.1;

%%% Filter parameter
options.sigma = 10;

%%% Patch radius
options.Nlx = 2;
options.Nly = 2;
options.Nlz = 2;

%%% Standard deviation of the Gaussian filter
options.NLM_gauss = 1;

% Search window radius is controlled by options.Ndx, options.Ndy and options.Ndz parameters
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
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RDP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Regularization parameter for quadratic prior with OSL-OSEM
options.beta_RDP_OSL_OSEM = 0.1;
%%% Regularization parameter for RDP with OSL-MLEM
options.beta_RDP_OSL_MLEM = 0.1;
%%% Regularization parameter for RDP with MBSREM
options.beta_RDP_MBSREM = 0.001;
%%% Regularization parameter for RDP with BSREM
options.beta_RDP_BSREM = 0.03;
%%% Regularization parameter for RDP with ROSEM
options.beta_RDP_ROSEM_MAP = 0.1;
%%% Regularization parameter for RDP with RBI
options.beta_RDP_OSL_RBI = 0.0001;
%%% Regularization parameter for RDP (OSL-(A)COSEM)
options.beta_RDP_OSL_COSEM = 0.01;
%%% Regularization parameter for RDP with PKMA
options.beta_RDP_PKMA = 0.0001;

%%% Edge preservation constant
options.RDP_gamma = 0.1;

%%% Pixel weights for RDP
% The number of pixels need to be the amount of neighboring pixels,
% e.g. if the above Nd(x/y/z) values are all 1, then 27 weights need to be
% included where the center pixel (if Nd values are 1, element 14) has no
% effect, i.e. it will always be one. Size is (Ndx*2+1) * (Ndy*2+1) *
% (Ndz*2+1). If left empty then  they will be calculated by OMEGA and are
% based on the distance of the voxels from the center voxel. If you do not
% wish to use weights, simply use ones, i.e.  
% options.weights_RDP = ones((options.Ndx*2+1) * (options.Ndy*2+1) * (options.Ndz*2+1),1);
%options.weights_RDP = [];


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
% ArrayFire_OpenCL_device_info();
 
%%% Implementation 3
% Uncomment the below line and run it to determine the available platforms,
% their respective numbers and device numbers
% OpenCL_device_info();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% CUSTOM DETECTOR COORDINATES %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load your custom detector coordinates and replace the below examples. In
% the example below x contains the detector coordinates for the x-direction
% and so on. For sinogram data, the dimensions need to be the following:
% size(x) = [(number of elements in a SINGLE sinogram) 2]
% size(y) = [(number of elements in a SINGLE sinogram) 2]
% size(z_det) = [(total number of sinograms) 2]
% E.g. each column represents the detector coordinates of one of the
% detectors in a LOR.
% For raw data:
% size(x) = [(detectors per ring) 1]
% size(y) = [(detectors per ring) 1]
% size(z_det) = [(total number of rings) 1]
% For raw data, as long as it is formatted as in OMEGA, you only need each
% detector coordinate once.

% options.x = x;
% options.y = y;
% options.z_det = z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD MEASUREMENT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if options.only_system_matrix == false
%     % Load the measurement data
%     options = loadMeasurementData(options);
% else
%     if options.use_raw_data == false
%         options.SinM = false(options.Ndist, options.Nang, options.NSinos);
%     else
%         options.coincidences = {false(options.detectors ^2/2 + options.detectors/2,1)};
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% DEPTH OF INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uncomment the below value to set a depth of interaction (mm)
% NOTE: By default this is set to 0, i.e. it is assumed that all the
% interactions occur at the surface of the detector crystals. What this
% value changes is the depth of where the interactions are assumed to
% occur, i.e. it only changes the detector coordinates such that the
% transaxial coordinates are "deeper" in the crystal.
% options.DOI = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data

% If you use SIMIND data or other data, you need to load the data manually.
% Simply input the projection images into options.SinM variable instead of
% loading the GATE data.
options = loadGATESPECTData(options);

%% Reconstructions

    
    tStart = tic;
    pz = reconstructions_mainSPECT(options);
    tElapsed = toc(tStart);
    disp(['Reconstruction process took ' num2str(tElapsed) ' seconds'])

	% save([options.name '_reconstruction_' num2str(options.subsets) 'subsets_' num2str(options.Niter) 'iterations_' ...
	%     num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '.mat'], 'pz');
    