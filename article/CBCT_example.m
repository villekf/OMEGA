function [z, t, recPar] = CBCT_example(type, varargin)

%% MATLAB/Octave codes for CBCT reconstruction
% This example file computes the results shown in the OMEGA v2 paper. The
% input variable "type" determines the example case. Type 0 corresponds to
% the top left, type 1 to top right, type 2 to bottom left and type 3 to
% bottom right. The input data
% Planmeca_VisoG7_100kV_80mAs_500proj_kneePhantom.mat has to be in
% MATLAB/Octave path!
% Implementation 2 is used by default, but the device number can be set
% after the type, e.g. CBCT_example(type, 1)
clear mex

if ~isempty(varargin) && nargin > 2 && ~isempty(varargin{2})
    largeDim = varargin{2};
else
    largeDim = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SCANNER PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Name of current datafile/examination
% This is used for (potential) naming purposes only
options.name = 'Planmeca_CT_data';

%%% Show status messages
% completed. It is recommended to keep this at 1 or 2. With value of 2, 
% you get more detailed timing information. Maximum is 3. Minimum is 0.
options.verbose = 2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Planmeca_VisoG7_100kV_80mAs_500proj_kneePhantom.mat

% Flat field corrected projections
options.SinM = proj;
clear proj

% Number of projections
options.nProjections = nProjections;

% Field of view
options.FOVa_x = FOV(1);
options.FOVa_y = FOV(2);
options.axial_fov = FOV(3);

% Number of rows and columns in a single projection image
options.nRowsD = nRowsD;
options.nColsD = nColsD;

% Object offset values from the origin, i.e. how much is the origin of the
% FOV shifted
% That is row, column and slice directions
options.oOffsetX = oOffset(1);
options.oOffsetY = oOffset(2);
options.oOffsetZ = oOffset(3);

% Flat value
options.flat = flatValue;

% Distance to center of rotation and detector
options.sourceToCRot = sourceToCRot;
options.sourceToDetector = sourceToDetector;

% Detector pixel size
options.dPitchX = dPitch(1);
options.dPitchY = dPitch(2);

% Projection angles
options.angles = projAngles;

% Rotation of the detector panel
options.pitchRoll = panelRot;

% Coordinates for the source and center of the detector
options.x = xCoord;
options.y = yCoord;
options.z = zCoord;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% IMAGE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The image size is taken from a conf file, but you can manually adjust
% these if desired
%%% Reconstructed image pixel count (X/row-direction)
options.Nx = 801;

%%% Y/column-direction
options.Ny = 801;

%%% Z-direction (number of slices) (axial)
options.Nz = 668;

% Use these two to rotate/flip the final image
%%% Flip the image (in horizontal direction)?
options.flip_image = true;

%%% How much is the image rotated (radians)?
% The angle (in radians) on how much the image is rotated BEFORE
% reconstruction, i.e. the rotation is performed in the detector space.
% Positive values perform the rotation in counter-clockwise direction
options.offangle = (3*pi)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Use projection extrapolation
% If true, extrapolates the projection data. You can select below whether
% this extrapolation is done only in the axial or transaxial directions, or
% both. Default extrapolation length is 20% of the original length, for
% both sides. For example if axial extrapolation is enabled, then the left
% and right regions of the projection get 20% increase in size. This value
% can be adjusted in CTEFOVCorrection. The values are scaled to air with
% the use of logarithmic scaling.
options.useExtrapolation = false;

%%% Use extended FOV
% Similar to above, but expands the FOV. The benefit of expanding the FOV
% this way is to enable to the use of multi-resolution reconstruction or
% computation of the priors/regularization only in the original FOV. The
% default extension is 40% per side (see below).
options.useEFOV = true;

% Use transaxial extended FOV (this is off by default)
options.transaxialEFOV = false;

% Use axial extended FOV (this is on by default. If both this and
% transaxialEFOV are false but useEFOV is true, the axial EFOV will be
% turned on)
options.axialEFOV = true;

% Same as above, but for extrapolation. Same default behavior exists.
options.transaxialExtrapolation = false;

% Same as above, but for extrapolation. Same default behavior exists.
options.axialExtrapolation = true;

% Setting this to true uses multi-resolution reconstruction when using
% extended FOV. Only applies to extended FOV!
options.useMultiResolutionVolumes = true;

% This is the scale value for the multi-resolution volumes. The original
% voxel size is divided by this value and then used as the voxel size for
% the multi-resolution volumes. Default is 4 times the original voxel size.
% This means that the multi-resolution regions have larger voxel sizes if
% this is < 1, i.e. 1/4 = 4 times the original voxel size.
options.multiResolutionScale = 1/4;

% Performs the extrapolation and adjusts the image size accordingly
options = CTEFOVCorrection(options);


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
% See the doc for more information:
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
if ~isempty(varargin) && ~isempty(varargin{1})
    options.use_device = varargin{1};
else
    options.use_device = 0;
end

% Implementation 2 ONLY
%%% Use CUDA
% Selecting this to true will use CUDA kernels/code instead of OpenCL. This
% only works if the CUDA code was successfully built. This is recommended
% if you have CUDA-capable device.
options.use_CUDA = checkCUDA(options.use_device);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROJECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Type of projector to use for the geometric matrix
% 1 = Improved/accelerated Siddon's algorithm
% 2 = Orthogonal distance based ray tracer
% 3 = Volume of intersection based ray tracer
% 4 = Interpolation-based projector (ray- and voxel-based)
% 5 = Branchless distance-driven projector
% NOTE: You can mix and match most of the projectors. I.e. 45 will use
% interpolation-based projector for forward projection while branchless
% distance-driven is used for backprojection
% NOTE 2: The below additional options apply also in hybrid cases as long
% as the other projector is the corresponding projector.
% See the doc for more information:
% https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
options.projector_type = 4;

%%% Interpolation length (projector type = 4 only)
% This specifies the length after which the interpolation takes place. This
% value will be multiplied by the voxel size which means that 1 is
% the interpolation length corresponding to a single voxel (transaxial)
% length. Larger values lead to faster computation but at the cost of
% accuracy. Recommended values are between [0.5 1].
options.dL = 1;


%%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Number of iterations (all reconstruction methods)
if type == 0 || type == 1
    options.Niter = 8;
else
    options.Niter = 50;
end
%%% Save specific intermediate iterations
% You can specify the intermediate iterations you wish to save here. Note
% that this uses zero-based indexing, i.e. 0 is the first iteration (not
% the initial value). By default only the last iteration is saved. Only
% full iterations (epochs) can be saved.
options.saveNIter = [];
% Alternatively you can save ALL intermediate iterations by setting the
% below to true and uncommenting it. As above, only full iterations
% (epochs) are saved.
% options.save_iter = false;

%%% Number of subsets (excluding subset_type = 6)
options.subsets = 20;

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
% 8 = Use every nth projection image
% 9 = Randomly select the projection images
% 10 = Use golden angle sampling to select the subsets (not recommended for
% PET)
% 11 = Use prime factor sampling to select the projection images
% Most of the time subset_type 8 is sufficient.
options.subset_type = 8;

%%% Initial value for the reconstruction
options.x0 = ones(options.Nx, options.Ny, options.Nz) * 1e-5;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NON-REGULARIZED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% LSQR
options.LSQR = false;

%%% Conjugate Gradient Least-squares (CGLS)
options.CGLS = false;

options.SART = false;

%%% Feldkamp-Davis-Kress (FDK)
options.FDK = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAP-METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These algorithms can utilize any of the selected priors, though only one
% prior can be used at a time

%%% Preconditioned Krasnoselskii-Mann algorithm (PKMA)
% Supported by implementations 1, 2, 4, and 5
if type == 2
    options.PKMA = true;
else
    options.PKMA = false;
end

%%% Primal-dual hybrid gradient (PDHG)
% Supported by implementations 1, 2, 4, and 5
if type == 0 || type == 1
    options.PDHG = true;
else
    options.PDHG = false;
end

%%% Primal-dual hybrid gradient (PDHG) with L1 minimization
% Supported by implementations 1, 2, 4, and 5
options.PDHGL1 = false;

%%% Primal-dual Davis-Yin (PDDY)
% Supported by implementation 2
options.PDDY = false;

options.ASD_POCS = false;

options.SART = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRIORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.TV = false;

%%% Total Generalized Variation (TGV) prior
options.TGV = false;

%%% Proximal TV
options.ProxTV = false;

%%% Non-local Means (NLM) prior
if type == 1
    options.NLM = true;
else
    options.NLM = false;
end

%%% Relative difference prior
if type == 2
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
% if you have a multi-resolution situation, you should input the values
% for each volume or use zero/empty
% options.tauCP = 0;
options.tauCP = [1/7540.685059;1/117180.726562;1/92894.101562];
% Primal value for filtered iterations, applicable only if
% options.precondTypeMeas[2] = true. As with above, automatically computed
% if left zero or empty. Same restrictions apply here as above.
% Use the "Largest eigenvalue for volume 0 with filtering" value here!
% if you have a multi-resolution situation, you should input the values
% for each volume or use zero/empty
% options.tauCPFilt = 0;
options.tauCPFilt = [1/127.564476;1/1606.458252;1/1479.921875];
% Dual value. Recommended to set at 1.
options.sigmaCP = 1;
% Next estimate update variable
options.thetaCP = 1;

% Use adaptive update of the primal and dual variables
% Currently two methods available
% Setting this to 1 or 2 uses an adaptive update for both the primal and 
% dual variables.
% Can lead to unstable behavior when using with multi-resolution
% Minimal to none use with filtering-based preconditioner
options.PDAdaptiveType = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRECONDITIONERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Applies to PDHG, PDHGL1, PDHGKL, PKMA, MBSREM, MRAMLA, PDDY, FISTA and
%%% FISTAL1
% Measurement-based preconditioners
% precondTypeMeas(1) = Diagonal normalization preconditioner (1 / (A1))
% precondTypeMeas(2) = Filtering-based preconditioner
if type == 0 || type == 1
    options.precondTypeMeas = [false;true];
else
    options.precondTypeMeas = [false;false];
end

% Number of filtering iterations
% Applies to both precondTypeMeas(2) and precondTypeImage(6)
options.filteringIterations = 80;

if type == 2
    options.precondTypeImage = [false;true;false;false;false;false;false];
else
    options.precondTypeImage = [false;false;false;false;false;false;false];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PKMA PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Relaxation parameter for PKMA
% If a scalar (or an empty) value is used, then the relaxation parameter is
% computed automatically as lambda(i) = (1 / ((i - 1)/20 + 1)) / 10000,
% where i is the iteration number. The input number thus has no effect.
% If, on the other hand, a vector is input then the input lambda values are
% used as is without any modifications (the length has to be at least the
% number of iterations).
options.lambda = 0;
if type == 2
    options.lambda = zeros(options.Niter,1);
    for i = 0 : options.Niter - 1
        options.lambda(i + 1) = 1 / ((1/300) * i + 1) / 1 * 1.0e-4;
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


%%%%%%%%%%%%%%%%%%%%%%%%% REGULARIZATION PARAMETER %%%%%%%%%%%%%%%%%%%%%%%%
%%% The regularization parameter for ALL regularization methods (priors)
if type == 1
    options.beta = 4; % NLRDP
elseif type == 2
    options.beta = 600; % RDP
end


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
options.Ndz = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NLM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filter parameter
% Higher values smooth the image, smaller values make it sharper
options.sigma = 1.5e-3;

%%% Patch radius
options.Nlx = 1;
options.Nly = 1;
options.Nlz = 1;

%%% Standard deviation of the Gaussian filter
options.NLM_gauss = 2;

% By default, the original NLM is used. You can, however, use another
% potential function by selecting ONE of the options below.
%%% Use Non-local total variation (NLTV)
% If selected, will overwrite regular NLM regularization as well as the
% below MRP version.
options.NLTV = false;

%%% Use Non-local relative difference (NLRD)
options.NLRD = true;

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
% Note that this affects NLRD as well
if type == 1
    options.RDP_gamma = 20; % NLRDP/CV
elseif type == 2
    options.RDP_gamma = 200; % RDP/PKMA
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Only load the current subset into the device (GPU) memory
options.loadTOF = false;

if largeDim
    options.largeDim = varargin{2};
end

%%

t = tic;
% pz is the reconstructed image volume
% recPar is a short struct of various reconstruction parameters used,
% useful when saving the reconstructed data with metadata
% classObject is the used class object in the reconstructions,
% essentially modified version of the input options-struct
% fp are the forward projections, if stored
% the primal-dual gap can be also be stored and is the variable after
% fp
[pz, recPar] = reconstructions_mainCT(options);
t = toc(t);

% Convert to HU-values
z = int16(pz(:,:,:,end) * 55000) - 1000;

volume3Dviewer(z, [-1000 2000], [0 0 1])
