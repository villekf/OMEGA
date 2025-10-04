%% MATLAB codes for CT reconstruction using Skyscan µCT projection images
% This example uses Skyscan µCT projection images. Furthermore, this example
% showcases the use of the "scalable" reconstruction for high-dimensional
% cases. This means that data of dozens of gigabytes can be successfully
% reconstructed on GPUs that cannot store the projection data and the final
% reconstructed image volume. The caveat is that functionality is limited.
% Only some of the algorithms are supported, such as PDHG, FDK and PKMA.
% Furthermore, only filtering-based preconditioner is supported.
% Multi-resolution reconstruction is not supported. The data and image are
% divided into options.subsets number of segments where only one segment is
% present at the GPU at a time. This means that the more subsets you use,
% the less memory will be used on the GPU side. The intermediate data will
% be stored in host (CPU) so high physical memory amount is recommended.
% Only subset types 1 and 8 are supported, though 8 should be used with CT
% data. Furthermore, only projector_type = 4 is supported. Only some of the
% regularization methods are supported, such as RDP, TV, NLM and GGMRF.
% You can use https://doi.org/10.5281/zenodo.12744181 as example data

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SCANNER PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This setting determines whether the high-dimensional scalable
% reconstruction is used (if set as true). Otherwise, the regular
% reconstruction is performed.
% NOTE: Currently the high-dimensional reconstructions are scaled
% differently than the regular ones
options.largeDim = true;

%%% Binning
% The level of binning used for the raw data. For example binning of 2
% reduces the size of the projections by two from both dimensions (e.g.
% 2048x3072 becomes 1024x1536).
options.binning = 1;

%%% Name of current datafile/examination
% This is used for naming purposes only
options.name = 'Jyvat_Skyscan';

%%% Show status messages
% These are e.g. time elapsed on various functions and what steps have been
% completed. It is recommended to keep this at 1 or 2. With value of 2, 
% you get more detailed timing information. Maximum is 3.
options.verbose = 1;

%%% Transaxial FOV size (mm), this is the length of the x (vertical/row) side
% of the FOV
options.FOVa_x = 13;

%%% Transaxial FOV size (mm), this is the length of the y (horizontal/column) side
% of the FOV
options.FOVa_y = options.FOVa_x;

%%% Axial FOV (mm)
options.axial_fov = 7.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if left blank, you will be prompted for the log-file. Alternatively,
% input the full path to the desired log file here
options.fpath = '/path/to/jyvat.log';


options = loadSkyscanData(options);
%%% Projection angles (degree or radian)
% The angles corresponding to the projections
options.angles = -(0:0.2:0.2*options.nProjections - 0.2)';

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

%%% Reconstructed image pixel count (X/row-direction)
options.Nx = 1000 * 2;

%%% Y/column-direction
options.Ny = 1000 * 2;

%%% Z-direction (number of slices) (axial)
options.Nz = 481 * 2;

%%% Flip the image (in column direction)?
options.flip_image = false;

%%% How much is the image rotated (radians)?
% The angle (in radians) on how much the image is rotated BEFORE
% reconstruction, i.e. the rotation is performed in the detector space.
% Positive values perform the rotation in counter-clockwise direction
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
% 2 = Matrix-free reconstruction with OpenCL/CUDA (Recommended)
% (Requires ArrayFire)
% See the doc for more information:
% https://omega-doc.readthedocs.io/en/latest/implementation.html
% No other implementation is supported when using largeDim!
options.implementation = 2;

%%% Device (GPU) used
% In implementation 2 this determines the device used image reconstruction.
% NOTE: Use ArrayFire_OpenCL_device_info() to determine the device numbers.
% NOTE: if you switch devices then you might need to run the below line
% (uncommented) as well:
% clear mex
options.use_device = 0;

%%% Use CUDA
% Selecting this to true will use CUDA kernels/code instead of OpenCL. This
% only works if the CUDA code was successfully built. This is recommended
% if you have CUDA-capable device.
options.use_CUDA = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROJECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Type of projector to use for the geometric matrix
% 4 = Interpolation-based projector (ray- and voxel-based)
% See the doc for more information:
% https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
options.projector_type = 4;

%%% Use mask
% The mask needs to be a binary mask (uint8 or logical) where 1 means that
% the pixel is included while 0 means it is skipped. Separate masks can be
% used for both forward and backward projection and either one or both can
% be utilized at the same time. E.g. if only backprojection mask is input,
% then only the voxels which have 1 in the mask are reconstructed.
% The mask can be either a 2D image that is applied identically to each slice
% or a 3D mask that is applied as-is
% Forward projection mask
% If nonempty, the mask will be applied. If empty, or completely omitted, no
% mask will be considered.
% options.maskFP = true(options.nRowsD,options.nColsD);
% Backprojection mask
% If nonempty, the mask will be applied. If empty, or completely omitted, no
% mask will be considered.
% The default mask uses the mask radius parameter from the Nikon xtekct-file
[columnsInImage, rowsInImage] = meshgrid(1:options.Nx, 1:options.Ny);
centerX = options.Nx/2;
centerY = options.Ny/2;
radius = options.Nx/2;
options.maskBP = uint8((rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2);


%%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Number of iterations (all reconstruction methods except FDK)
% Note that using FDK automatically sets this to 1
options.Niter = 5;

%%% Number of subsets
% Note that with high-dimensional data this is required for FDK as well.
% As mentioned above, for high-dimensional data this controls the amount of
% memory required by the GPU. More subsets, less memory, but using too many
% subsets can lead to reduced performance.
options.subsets = 20;

%%% Subset type (n = subsets)
% 8 = Use every nth projection image
% For high-dimensional data, do not change this!
options.subset_type = 8;

%%% Initial value for the reconstruction
% Should not be used with high-dimensional reconstruction
if ~options.largeDim
    options.x0 = ones(options.Nx, options.Ny, options.Nz,'single') * 1e-4;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction algorithms to use (you can choose one)

%%% Feldkamp-Davis-Kress (FDK)
options.FDK = true;


%%% Preconditioned Krasnoselskii-Mann algorithm (PKMA)
options.PKMA = false;

%%% Primal-dual hybrid gradient (PDHG)
options.PDHG = false;

%%% Primal-dual hybrid gradient (PDHG) with L1 minimization
options.PDHGL1 = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These priors do not work with FDK!
%%% Proximal TV
options.ProxTV = false;

%%% Non-local Means (NLM) prior
options.NLM = false;

%%% Relative difference prior
options.RDP = false;

%%% Generalized Gaussian Markov random field (GGMRF) prior
options.GGMRF = false;


%%%%%%%%%%%%%%%%%%%%%%%%% REGULARIZATION PARAMETER %%%%%%%%%%%%%%%%%%%%%%%%
%%% The regularization parameter for ALL regularization methods (priors)
% ~.1 is good starting region for RDP and NLRD
options.beta = .1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%% ENFORCE POSITIVITY %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Applies to PDHG, PDHGL1, PKMA
% Enforces positivity in the estimate after each iteration
options.enforcePositivity = true;


%%%%%%%%%%%%%%%%%%%% MEASUREMENT-DOMAIN PRECONDITIONERS %%%%%%%%%%%%%%%%%%%
% 2 = Filtering-based preconditioner
% Note: At the moment it is not recommended to use filtering with PKMA, use
% at your own risk!
options.precondTypeMeas = [false;true];

%%% Filtering-steps
% The number of filtering steps for image- and measurement-domain
% filtering. Includes subsets/sub-iterations as well
% E.g. if you have 10 filtering steps, 3 iterations, and 5 subsets, the
% first 2 iterations will use filtering
% There are no restrictions on when to stop the filtering
options.filteringIterations = 200;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PDHG PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Primal value
% If left zero, or empty, it will be automatically computed
% Note that if you change any of the model parameters, i.e. image volume
% size, number of projections or use binning, this needs to be recomputed
% or scaled accordingly!
% The computed largest eigenvalue is printed if verbose > 0. This can be 
% used as the below value as long as one is divided by it. For example, 
% if "Largest eigenvalue for volume 0 is 100" then options.tauCP should be 
% 1/100 (if you use filtering-based preconditioner this is the "without 
% filtering" value)
options.tauCP = 0;
% Primal value for filtered iterations, applicable only if
% options.precondTypeMeas[2] = true. As with above, automatically computed
% if left zero or empty. Same restrictions apply here as above.
% Use the "Largest eigenvalue for volume 0 with filtering" value here!
options.tauCPFilt = 0;
% Dual value. Recommended to set at 1.
options.sigmaCP = 1;
% Next estimate update variable, recommended to keep at 1.
options.thetaCP = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PKMA PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Relaxation parameter for PKMA
% If a scalar (or an empty) value is used, then the relaxation parameter is
% computed automatically as lambda(i) = (1 / ((i - 1)/20 + 1)) / 10000,
% where i is the iteration number. The input number thus has no effect.
% If, on the other hand, a vector is input then the input lambda values are
% used as is without any modifications (the length has to be at least the
% number of iterations).
options.lambda = 0;

% If the reconstruction doesn't work or there are holes or bright spots, 
% then this relaxation value is too high. Reduce the last value (default 
% is 8.1e-7) in such cases and try again.
options.lambda = zeros(options.Niter,1);
for i = 0 : options.Niter - 1
    if sum(options.precondTypeMeas) == 0
        options.lambda(i + 1) = 1 / ((1/3500) * i + 1) / 1 * 8.1e-7;
    else
        options.lambda(i + 1) = 1 / ((1/3500) * i + 1) / 1 * 5.1e-3;
    end
end

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
options.delta_PKMA = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NLM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filter parameter
% Higher values smooth the image, smaller values make it sharper
options.sigma = 6.00e-3;

%%% Patch radius
% Works exactly the same as the neighborhood size
options.Nlx = 1;
options.Nly = 1;
options.Nlz = 1;

%%% Standard deviation of the Gaussian-weighted Euclidean norm
options.NLM_gauss = 2;

% By default, the original NLM is used. You can, however, use another
% potential function by selecting ONE of the options below.
% Note that only one of the below options for NLM can be selected!
% If all the below ones are false, regular NLM is used!
%%% Use Non-local total variation (NLTV)
options.NLTV = false;

%%% Use Non-local relative difference (NLRD)
options.NLRD = false;

%%% Use Non-local Lange prior (NLLange)
options.NLLange = false;

% Tuning parameter for Lange
options.SATVPhi = 5;

%%% Use Non-local GGMRF (NLGGMRF)
options.NLGGMRF = false;

%%% Use MRP algorithm (without normalization)
% I.e. gradient = im - NLM_filtered(im)
options.NLM_MRP = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RDP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Edge weighting factor
% Higher values sharpen the image, smaller values make it smoother
% Note that this affects NLRD as well
options.RDP_gamma = 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GGMRF PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GGMRF parameters
% See the original article for details
% These affect the NLGGMRF as well
% https://omega-doc.readthedocs.io/en/latest/algorithms.html#ggmrf
options.GGMRF_p = 1.5;
options.GGMRF_q = 1;
options.GGMRF_c = 5;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

pz = reconstructions_mainCT(options);
volume3Dviewer(pz, [], [0 0 1])