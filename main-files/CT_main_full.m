%% MATLAB/Octave code for CT reconstruction
% This example file lists ALL adjustable parameters
% New parameters are in scanner properties and reconstruction parameters,
% and new section below reconstruction parameters (and above OpenCL device
% info)
% You can use the FIPS walnut data as an example data:
% https://zenodo.org/records/6986012
 
clear

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

%%% Number of detector pixels (vertical/row direction)
% The number of detector pixels in the detector panel (vertical
% direction/number of rows)
% NOTE: if you use binning, this value has to use the final binned
% dimensions
options.nRowsD = 2240/options.binning;

%%% Number of detector pixels (horizontal/column direction)
% The number of detector pixels in the detector panel (horizontal
% direction/number of columns)
% NOTE: if you use binning, this value has to use the final binned
% dimensions
options.nColsD = 2368/options.binning;

%%% Number of projections
% Total number of projections used
options.nProjections = 721;

%%% Projection angles (degree or radian)
% The angles corresponding to the projections
options.angles = -linspace(0, 360, options.nProjections);

%%% Detector pixel pitch/size (mm), row direction
% The size of the detector/distance between adjacent detectors
% NOTE: if you use binning, this value has to use the final binned
% dimensions
options.dPitchX = 0.05*options.binning;

%%% Crystal pitch/size in x- and y-directions (transaxial) (mm)
% Same as above, but different name for PET
options.cr_p = 0.05*options.binning;

%%% Detector pixel pitch/size (mm), column direction
% The size of the detector/distance between adjacent detectors
% NOTE: if you use binning, this value has to use the final binned
% dimensions
options.dPitchY = 0.05*options.binning;

%%% Crystal pitch/size in z-direction (axial) (mm)
% Same as above, but different name for PET
options.cr_pz = 0.05*options.binning;

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
% completed. It is recommended to keep this at 1 or 2. With value of 2, 
% you get more detailed timing information. Maximum is 3. Minimum is 0.
options.verbose = 1;

% Note that non-square transaxial FOV sizes should work, but might not work
% always. Square transaxial FOV is thus recommended.
%%% Transaxial FOV size (mm), this is the length of the x (vertical/row) side
% of the FOV
options.FOVa_x = 40.1;

%%% Transaxial FOV size (mm), this is the length of the y (horizontal/column) side
% of the FOV
options.FOVa_y = options.FOVa_x;

% The above recommendation doesn't apply to axial FOV, i.e. this can be
% different from the transaxial FOV size(s). 
%%% Axial FOV (mm)
options.axial_fov = 40;

%%% Source row offset (mm)
% The center of rotation is not exactly in the origin. With this parameter
% the source location can be offset by the specified amount (row direction).
% This has a similar effect as circularly shifting the projection images.
% Use vector values if these vary with each projection (this is untested at
% the moment).
% If inputing a vector, use a different value for each projection
% NOTE: The default value has been obtained experimentally and is not based
% on any known value.
options.sourceOffsetRow = -0.16;

%%% Source column offset (mm)
% Same as above, but for column direction.
options.sourceOffsetCol = 0;

%%% Detector panel row offset (mm)
% Same as above, but the offset value for the detector panel.
options.detOffsetRow = 0;

%%% Detector panel column offset (mm)
% Same as above, but for column direction.
options.detOffsetCol = 0;

%%% Bed offset (mm)
% The offset values for multi-bed step-and-shoot examinations. Each bed
% position should have its own offset value.
options.bedOffset = [];

%%% Pitch/roll angles for the detector panel (radian)
% Sometimes the detector panel is slightly rotated in all three directions.
% The main rotation is included in the above options.angles, but the other
% two directions can be included here. pitchRoll should be column vector,
% where the first column corresponds to the rotation in the XY-plane and
% the second to rotation in the ZY-plane. 
options.pitchRoll = [];

%%% Direction vectors (normalized)
% This one is optional, but you can also input straight the direction
% vectors for all dimensions. The number of required dimensions depends on
% the axes where rotation occurs. If pitchRoll would be empty, i.e.
% rotation is only in the XY-plane (angles) then only two dimensions are
% needed, one for X- and one for Y-direction (row direction). Z-direction
% (column direction) is handled automatically. If the other two rotations
% are included, then six dimensions are needed. The first three are the
% direction vectors for the placement in the row-direction. I.e. they are
% used to determine the current detector pixel location in the
% row-direction. The latter three are for the detector pixel location in
% the column direction. The dimensions should be such that the number of
% rows for uV should be either 2 or 6, while the number of columns should
% equal the number of projections. Note that these values have to be
% normalized values as they are later multiplied with the size of the
% detector pixel. If you input the above pitchRoll, these are computed
% automatically or if there is no rotation other than the general rotation
% of angles, then this is also generally not required.
options.uV = [];

% Note: Origin is assumed to be at the center. If this is not the case, you
% can shift it with options.oOffsetX, options.oOffsetY and options.oOffsetZ
% That is row, column and slice directions
% options.oOffsetX = 0;
% options.oOffsetY = 0;
% options.oOffsetZ = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Path to the first image, alternatively, you can leave this as is and 
% simply input the data when prompted
options.fpath = '/path/to/20201111_walnut_0001.tif';

if ~options.only_reconstructions || ~isfield(options,'SinM')
    options.SinM = loadProjectionImages(options.nProjections,options.binning,options.fpath);
    % options.SinM = single(options.SinM) ./ single(max(max(max(options.SinM(4:end-3,:,:)))));
    % options.SinM(options.SinM > 1) = single(1);
    options.SinM = permute(options.SinM, [2 1 3]);
end
% NOTE: If you want to reduce the number of projections, you need to do
% this manually as outlined below:
% options.SinM = options.SinM(:,:,1:4:options.nProjections);
% options.angles = options.angles(1:4:numel(options.angles));
% options.pitchRoll = pitchRoll(1:4:options.nProjections, :);
% options.uV = options.uV(:, 1:4:options.nProjections);
% options.nProjections = numel(options.angles);


% Flat value
% Needed for both linearized and Poisson-based data
% If omitted, the maximum value will be used automatically
% options.flat = [];


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
options.Nx = 280;

%%% Y/column-direction
options.Ny = 280;

% The above, again, doesn't apply to axial direction
% i.e. the number of slices can differ from Nx or Ny
%%% Z-direction (number of slices) (axial)
options.Nz = 280;

%%% Flip the image (in column direction)?
options.flip_image = false;

%%% How much is the image rotated (radians)?
% The angle (in radians) on how much the image is rotated BEFORE
% reconstruction, i.e. the rotation is performed in the detector space.
% Positive values perform the rotation in counter-clockwise direction
options.offangle = (2*pi)/2;

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

% The extrapolation length per side, i.e. the total size is this multiplied
% by two!
options.extrapLength = 0.2;

%%% Use extended FOV
% Similar to above, but expands the FOV. The benefit of expanding the FOV
% this way is to enable to the use of multi-resolution reconstruction or
% computation of the priors/regularization only in the original FOV. The
% default extension is 40% per side (see below).
options.useEFOV = false;

% The extended FOV length per side, i.e. the total size is this multiplied
% by two!
options.eFOVLength = 0.4;

% Use transaxial extended FOV (this is off by default)
options.transaxialEFOV = false;

% Use axial extended FOV (this is on by default. If both this and
% transaxialEFOV are false but useEFOV is true, the axial EFOV will be
% turned on)
options.axialEFOV = false;

% Same as above, but for extrapolation. Same default behavior exists.
options.transaxialExtrapolation = false;

% Same as above, but for extrapolation. Same default behavior exists.
options.axialExtrapolation = false;

% Setting this to true uses multi-resolution reconstruction when using
% extended FOV. Only applies to extended FOV!
options.useMultiResolutionVolumes = true;

% This is the scale value for the multi-resolution volumes. The original
% voxel size is divided by this value and then used as the voxel size for
% the multi-resolution volumes. Default is 1/4 of the original voxel size.
% This means that the multi-resolution regions have smaller voxel sizes if
% this is < 1.
options.multiResolutionScale = 1/4;

% Performs the extrapolation and adjusts the image size accordingly
options = CTEFOVCorrection(options);

% Use offset-correction
% If you use offset imaging, i.e. the center of rotation is not in the
% origin but rather a circle around the origin, you can enable automatic
% offset weighting by setting this to true.
options.offsetCorrection = false;


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

% Applies to implementations 3 and 5 ONLY
%%% OpenCL platform used 
% NOTE: Use OpenCL_device_info() to determine the platform numbers and
% their respective devices with implementations 3 or 5.
options.platform = 0;

% Applies to implementations 2, 3 and 5 ONLY
%%% OpenCL/CUDA device used 
% NOTE: Use ArrayFire_OpenCL_device_info() to determine the device numbers
% with implementation 2.
% NOTE: Use OpenCL_device_info() to determine the platform numbers and
% their respective devices with implementations 3 or 5.
% NOTE: The device numbers might be different between implementation 2 and
% implementations 3 and 5
% NOTE: if you switch devices then you might need to run the below line
% (uncommented) as well:
% clear mex
options.use_device = 0;

% Applies to implementations 2, 3 and 5 ONLY
% Applies to projector type 1-3 backprojection ONLY
%%% Use 64-bit integer atomic functions
% If true, then 64-bit integer atomic functions (atomic add) will be used
% if they are supported by the selected device.
% Setting this to true will make computations faster on GPUs that support
% the functions, but might make results slightly less reliable due to
% floating point rounding. Recommended for OpenCL GPUs.
options.use_64bit_atomics = false;

% Applies to implementations 2, 3 and 5 ONLY
% Applies to projector type 1-3 backprojection ONLY
%%% Use 32-bit integer atomic functions
% If true, then 32-bit integer atomic functions (atomic add) will be used.
% This is even faster than the above 64-bit atomics version, but will also
% have significantly higher reduction in numerical/floating point accuracy.
% This should be about 20-30% faster than the above 64-bit version, but
% might lead to integer overflow if you have a high count measurement
% (thousands of coincidences per sinogram bin). Use this only if speed is
% of utmost importance. 32-bit atomics take precedence over 64-bit ones,
% i.e. if options.use_32bit_atomics = true then the 64-bit version will be 
% always set as false.
options.use_32bit_atomics = false;

% Implementation 2 ONLY
%%% Use CUDA
% Selecting this to true will use CUDA kernels/code instead of OpenCL. This
% only works if the CUDA code was successfully built. This is recommended
% if you have CUDA-capable device.
options.use_CUDA = false;

% Implementation 2 ONLY
%%% Use CPU
% Selecting this to true will use CPU-based code instead of OpenCL or CUDA.
% Not recommended, even OpenCL with CPU should be used before this.
options.use_CPU = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROJECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Type of projector to use for the geometric matrix
% 1 = Improved/accelerated Siddon's algorithm
% 2 = Orthogonal distance based ray tracer (not recommended in CT)
% 3 = Volume of intersection based ray tracer (not recommended in CT)
% 4 = Interpolation-based projector (ray- and voxel-based)
% 5 = Branchless distance-driven projector
% NOTE: You can mix and match most of the projectors. I.e. 45 will use
% interpolation-based projector for forward projection while branchless
% distance-driven is used for backprojection
% NOTE 2: The below additional options apply also in hybrid cases as long
% as the other projector is the corresponding projector.
% See the documentation for more information:
% https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
% NOTE: Projector types 1-3 should not be used as backprojectors for CT!
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
% Create a circle that fills the FOV:
% [columnsInImage, rowsInImage] = meshgrid(1:options.Nx, 1:options.Ny);
% centerX = options.Nx/2;
% centerY = options.Ny/2;
% radius = options.Nx/2;
% options.maskBP = uint8((rowsInImage - centerY).^2 ...
%     + (columnsInImage - centerX).^2 <= radius.^2);

%%% Interpolation length (projector type = 4 only)
% This specifies the length after which the interpolation takes place. This
% value will be multiplied by the voxel size which means that 1 is
% the interpolation length corresponding to a single voxel (transaxial)
% length. Larger values lead to faster computation but at the cost of
% accuracy. Recommended values are between [0.5 1].
options.dL = 0.5;

%%% Use point spread function (PSF) blurring
% Applies PSF blurring through convolution to the image space. This is the
% same as multiplying the geometric matrix with an image blurring matrix.
options.use_psf = false;

% FWHM (mm) of the Gaussian used in PSF blurring in all three dimensions (X/Y/Z)
options.FWHM = [options.cr_p options.cr_p options.cr_pz];

% Use deblurring phase
% If enabled, a deblurring phase is performed once the reconstruction has
% completed. This step is performed for all iterations (deblurred estimates
% are NOT used in the reconstruction phase). This is used ONLY when PSF
% blurring is used.
options.deblurring = false;
% Number of deblurring iterations
% How many iterations of the deblurring step is performed
options.deblur_iterations = 10;

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

% projector_type = 1 and 4 only
%%% Number of rays
% Number of rays used per detector if projector_type = 1 (i.e. Improved
% Siddon is used) or projector_type = 4 (interpolation).
% The total number of rays per detector is the multiplication of the two
% below values!
% Number of rays in transaxial (row) direction
options.n_rays_transaxial = 1;
% Number of rays in axial (column) direction
options.n_rays_axial = 1;


%%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Number of iterations (all reconstruction methods)
options.Niter = 2;

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
% Note that with high-dimensional data this is required for FDK as well.
% For high-dimensional data this controls the amount of memory required by 
% the GPU. More subsets, less memory, but using too many subsets can lead 
% to reduced performance.
options.subsets = 1;

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
% 8 = Use every nth sinogram/projection
% 9 = Randomly select the full sinograms/projections
% 10 = Use golden angle sampling to select the subsets (not recommended for
% PET)
% 11 = Use prime factor sampling to select the full sinograms/projections
% Most of the time subset_type 8 is sufficient.
options.subset_type = 8;

%%% Stochastic subset selection
% If true, the subsets are selected stochastically
% This means that the subset numbers are selected randomly
% For example, if using subset_type = 8, the subsets are still grouped into
% groups with every nth projection but the group is selected randomly
% For example, if we have three subsets with 9 projections, the first group
% will have projection images 1, 4, and 7, second 2, 5, and 8, and the third
% 3, 6, and 9. During the reconstruction, the group is selected randomly, but
% the projections within the groups remain the same so first group always has
% projections 1, 4, and 7, but the first subiteration might use group three.
options.stochasticSubsetSelection = false;

%%% How many angles are combined in subset_type = 6
% E.g. there are 180 angles, in n_angles = 2, then angles 0 and 1 are
% combined to the same subset, 2 and 3, etc.
options.n_angles = 90;

%%% Initial value for the reconstruction
options.x0 = ones(options.Nx, options.Ny, options.Nz) * 1e-4;

%%% Epsilon value 
% A small value to prevent division by zero and square root of zero. Should
% not be smaller than machine epsilon (eps).
options.epps = 1e-5;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MISC SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Use Shuffle
% NOTE: Applies only when using subset_type = 3. 
% Accelerates the subset formation and uses less memory. Not included in
% OMEGA, needs to be manually downloaded and installed.
% Download from: 
% https://www.mathworks.com/matlabcentral/fileexchange/27076-shuffle
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
% Reconstruction algorithms to use (choose only one algorithm and
% optionally one prior)
% High-dimensional case ONLY supports FDK, PKMA, PDHG and PDHGL1!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ML-METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ordered Subsets Expectation Maximization (OSEM) OR Maximum-Likelihood
%%% Expectation Maximization (MLEM) (if subsets = 1)
% Supported by all implementations
options.OSEM = false;

%%% Modified Row-Action Maximum Likelihood Algorithm (MRAMLA)
% Supported by implementations 1, 2, 4, and 5
options.MRAMLA = false;

%%% Row-Action Maximum Likelihood Algorithm (RAMLA)
% Supported by implementations 1, 2, 4, and 5
options.RAMLA = false;

%%% Relaxed Ordered Subsets Expectation Maximization (ROSEM)
% Supported by implementations 1, 2, 4, and 5
options.ROSEM = false;

%%% LSQR
% Supported by implementations 1, 2, 4, and 5
options.LSQR = false;

%%% Conjugate Gradient Least-squares (CGLS)
% Supported by implementations 1, 2, 4, and 5
options.CGLS = false;

%%% Feldkamp-Davis-Kress (FDK)
% Supported by implementation 2
options.FDK = false;


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

%%% Preconditioned Krasnoselskii-Mann algorithm (PKMA)
% Supported by implementations 1, 2, 4, and 5
options.PKMA = false;

%%% Primal-dual hybrid gradient (PDHG)
% Supported by implementations 1, 2, 4, and 5
options.PDHG = true;

%%% Primal-dual hybrid gradient (PDHG) with L1 minimization
% Supported by implementations 1, 2, 4, and 5
options.PDHGL1 = false;

%%% Primal-dual Davis-Yin (PDDY)
% Supported by implementation 2
options.PDDY = false;

%%% Simultaneous ART
% Supported by implementation 2
options.SART = false;

%%% SAGA
% Supported by implementation 2
options.SAGA = false;

%%% Barzilai-Borwein
% Supported by implementation 2
options.BB = false;

%%% ASD-POCS
% Supported by implementation 2
options.ASD_POCS = false;



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

%%% Proximal TV
options.ProxTV = false;

%%% Total Generalized Variation (TGV) prior
options.TGV = false;

%%% Non-local Means (NLM) prior
options.NLM = false;

%%% Relative difference prior
options.RDP = false;

%%% Generalized Gaussian Markov random field (GGMRF) prior
options.GGMRF = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%% ENFORCE POSITIVITY %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Applies to PDHG, PDHGL1, PDDY, FISTA, FISTAL1, MBSREM, MRAMLA, PKMA
% Enforces positivity in the estimate after each iteration
options.enforcePositivity = false;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% ACOSEM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Acceleration parameter for ACOSEM (1 equals COSEM)
options.h = 2;


%%%%%%%%%%%%%%%%%%%%%%%%%% RELAXATION PARAMETER %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Relaxation parameter for MRAMLA, RAMLA, ROSEM, BSREM, MBSREM and PKMA
% If a scalar (or an empty) value is used, then the relaxation parameter is
% computed automatically as lambda(i) = (1 / ((i - 1)/20 + 1)) / 10000,
% where i is the iteration number. The input number thus has no effect.
% If, on the other hand, a vector is input then the input lambda values are
% used as is without any modifications (the length has to be at least the
% number of iterations).
options.lambda = 0;
 

%%%%%%%%%%%%%%%%%%%%%%%% MRAMLA & MBSREM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%
%%% Upper bound for MRAMLA/MBSREM (use 0 for default value)
options.U = 0;
 

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
options.tauCP = 0;
% Primal value for filtered iterations, applicable only if
% options.precondTypeMeas[2] = true. As with above, automatically computed
% if left zero or empty. Same restrictions apply here as above.
% Use the "Largest eigenvalue for volume 0 with filtering" value here!
% if you have a multi-resolution situation, you should input the values
% for each volume or use zero/empty
options.tauCPFilt = 0;
% Dual value. Recommended to set at 1.
options.sigmaCP = 1;
% Next estimate update variable
options.thetaCP = 1;
% Dual value for TV and/or TGV. For faster convergence, set this to higher
% than 1.
options.sigma2CP = 1;

% Use adaptive update of the primal and dual variables
% Currently two methods available
% Setting this to 1 or 2 uses an adaptive update for both the primal and 
% dual variables.
% Can lead to unstable behavior when using with multi-resolution
% Minimal to none use with filtering-based preconditioner
options.PDAdaptiveType = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FISTA PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FISTA step-size/acceleration type
% There are two different ways to compute the step-size/acceleration value
% for FISTA reconstructions. There should be slight convergence rate 
% differences between the two. Values 1 and 2 are supported.
options.FISTAType = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRECONDITIONERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Applies to PDHG, PDHGL1, PDHGKL, PKMA, MBSREM, MRAMLA, PDDY, FISTA, 
%%% FISTAL1, and SAGA
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
options.precondTypeImage = [false;false;false;false;false;false;false];

% Reference image for precondTypeImage(3). Can be either a mat-file or a
% variable
options.referenceImage = '';

% Momentum parameter for precondTypeImage(4)
% Set the desired momentum parameters to the following variable (note that
% the length should be options.Niter * options.subsets): 
% options.alphaPrecond = [];
% Otherwise set the following parameters:
options.rhoPrecond = options.rho_PKMA;
options.delta1Precond = options.delta_PKMA;

% Parameters for precondTypeImage(5)
% See the article for details:
% https://omega-doc.readthedocs.io/en/latest/algorithms.html#gradient-based-preconditioner
options.gradV1 = 1.5;
options.gradV2 = 2;
% Note that these include subiterations (options.Niter * options.subsets)
% The first iteration where to start the gradient computation
options.gradInitIter = 1;
% Last iteration of the gradient computation
options.gradLastIter = 100;

% Number of filtering iterations
% Applies to both precondTypeMeas(2) and precondTypeImage(6)
% The filtering is applies to this many (sub)iterations
% Note that this include subiterations (options.Niter * options.subsets)
options.filteringIterations = 100;


%%%%%%%%%%%%%%%%%%%%%%%%% REGULARIZATION PARAMETER %%%%%%%%%%%%%%%%%%%%%%%%
%%% The regularization parameter for ALL regularization methods (priors)
options.beta = 1;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%% NEIGHBORHOOD PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%
%%% How many neighboring pixels are considered 
% With MRP, QP, L, FMH, NLM, (RDP), GGMRF and weighted mean
% E.g. if Ndx = 1, Ndy = 1, Ndz = 0, then you have 3x3 square area where
% the pixels are taken into account (I.e. (Ndx*2+1)x(Ndy*2+1)x(Ndz*2+1)
% area).
% NOTE: Currently Ndx and Ndy must be identical.
% For NLM this is often called the "search window".
options.Ndx = 2;
options.Ndy = 2;
options.Ndz = 1;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pixel weights for quadratic prior
% The number of pixels need to be the amount of neighboring pixels,
% e.g. if the above Nd values are all 1, then 27 weights need to be
% included where the center pixel (if Nd values are 1, element 14) should
% be Inf. Size is (Ndx*2+1) * (Ndy*2+1) * (Ndz*2+1). If left empty then
% they will be calculated by the algorithm and are based on the distance of
% the voxels from the center.
options.weights = [];
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Delta parameter for Huber prior
% Upper and lower bounds for the prior
options.huber_delta = 5;

%%% Pixel weights for Huber prior
% Same rules apply as with quadratic prior weights.
% If left empty then they will be calculated by the algorithm and are based
% on the distance of the voxels from the center.
options.weights_huber = [];
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%% L-FILTER PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Weighting factors for the L-filter pixels
% Otherwise the same as in quadratic prior, but center pixel is not Inf.
% If left empty then they will be calculated by the algorithm such that the
% weights resemble a Laplace distribution.
options.a_L = [];

%%% If the weighting factors are set empty, then this option will determine
% whether the computed weights follow a 1D weighting scheme (true) or 2D 
% (false).
% See the docs for more information:
% https://omega-doc.readthedocs.io/en/latest/algorithms.html#l-filter
options.oneD_weights = false;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FMH PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Types 1-3 compute the weighted mean just as MRP is computed, but the
% median is replaced with the weighted mean.
% 1 = Arithmetic mean (MRP), 2 = Harmonic mean (MRP), 3 = Geometric mean
% (MRP)
% Types 4-6 compute the weighted mean around the neighborhood of the voxel
% and use joint estimation to compute the gradient where the other variable
% corresponds to the chosen mean value and the other is based on the chosen
% mean value. See the docs for more information.
% 4 = Arithmetic mean, 5 = Harmonic mean, 6 = Geometric mean
options.mean_type = 1;

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
% NOTE: These have no effect with proximal TV
%%% "Smoothing" parameter
% Also used to prevent zero values in square root.
options.TVsmoothing = 1e-5;

%%% Whether to use an anatomical reference/weighting image with the TV
options.TV_use_anatomical = false;

%%% If the TV_use_anatomical value is set to true, specify filename for the
% reference image here (same rules apply as with attenuation correction
% above). Alternatively you can specify the variable that holds the
% reference image.
options.TV_reference_image = 'reference_image.mat';

%%% Five different TV methods are available.
% Value can be 1, 2, 3, 4 or 6.
% Type 3 is not recommended!
% Types 1 and 2 are the same if anatomical prior is not included
% Type 3 uses the same weights as quadratic prior
% Type 4 is the Lange prior, does not support anatomic weighting.
% Type 6 is a weighted TV, does not support anatomic weighting.
% See the docs for more information:
% https://omega-doc.readthedocs.io/en/latest/algorithms.html#tv
options.TVtype = 1;

%%% Weighting parameters for the TV prior. 
% Applicable only if use_anatomical = true. T-value is specific to the used
% TVtype, e.g. for type 1 it is the edge threshold parameter. See the docs
% for more details:
% https://omega-doc.readthedocs.io/en/latest/algorithms.html#tv
options.T = 0.5;

%%% C is the weight for the original image in type 3 and is ignored with
% other types
options.C = 1;

%%% Tuning parameter for TV and APLS
options.tau = 1e-8;

%%% Tuning parameter for Lange function in SATV (type 4) or weight factor
%%% for weighted TV (type 6)
% Setting this to 0 gives regular anisotropic TV with type 4
% This affects also non-local Lange
options.SATVPhi = 0.2;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADMRP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Time step variable for AD (implementation 2 only)
options.TimeStepAD = 0.0625;

%%% Conductivity/connectivity for AD (edge threshold)
options.KAD = 2;

%%% Number of iterations for AD filter
% NOTE: This refers to the AD smoothing part, not the actual reconstruction
% phase.
options.NiterAD = 10;

%%% Flux/conduction type for AD filter
% 1 = Exponential
% 2 = Quadratic
options.FluxType = 1;

%%% Diffusion type for AD (implementation 2 only)
% 1 = Gradient
% 2 = Modified curvature
options.DiffusionType = 1;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% APLS PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Scaling parameter (eta)
% See the docs for details:
% https://omega-doc.readthedocs.io/en/latest/algorithms.html#tv
options.eta = 1e-5;

%%% "Smoothing" parameter (beta)
% Also used to prevent zero values in square root.
options.APLSsmoothing = 1e-5;

%%% Specify filename for the reference image here (same rules apply as with
% attenuation correction above). As before, this can also be a variable
% instead.
% NOTE: For APSL, the reference image is required!
options.APLS_reference_image = 'reference_image.mat';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HYPERBOLIC PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Edge weighting factor
options.hyperbolicDelta = 800;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TGV PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TGV weights
% First part
options.alpha0TGV = 1;
% Second part (symmetrized derivative)
options.alpha1TGV = 2;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NLM PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filter parameter
% Higher values smooth the image, smaller values make it sharper
options.sigma = 6e-3;

%%% Patch radius
% Works exactly the same as the neighborhood size
options.Nlx = 1;
options.Nly = 1;
options.Nlz = 1;

%%% Standard deviation of the Gaussian-weighted Euclidean norm
options.NLM_gauss = 2;

%%% Adaptive NL method
options.NLAdaptive = false;

%%% Summed constant for adaptive NL
options.NLAdaptiveConstant = 2.0e-7;

% Search window radius is controlled by Ndx, Ndy and Ndz parameters
% Use anatomical reference image for the patches
options.NLM_use_anatomical = false;

%%% Specify filename for the reference image here or the variable containing 
% the reference image or the variable containing the reference image
options.NLM_reference_image = 'reference_image.mat';

% Note that only one of the below options for NLM can be selected!
% If all the below ones are false, regular NLM is used!
%%% Use Non-local total variation (NLTV)
options.NLTV = false;

%%% Use Non-local Lange prior (NLLange)
options.NLLange = false;

%%% Use MRP algorithm (without normalization)
% I.e. gradient = im - NLM_filtered(im)
options.NLM_MRP = false;

%%% Use non-local relative difference prior (NLRD)
options.NLRD = false;

%%% Use non-local GGMRF (NLGGMRF)
options.NLGGMRF = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RDP PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Edge weighting factor
% Higher values sharpen the image, smaller values make it smoother
% Note that this affects NLRD as well
options.RDP_gamma = 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GGMRF PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GGMRF parameters
% These affect the NLGGMRF as well
% See the original article for details
% https://omega-doc.readthedocs.io/en/latest/algorithms.html#ggmrf
options.GGMRF_p = 1.5;
options.GGMRF_q = 1;
options.GGMRF_c = 5;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Store the intermediate forward projections. Unlike image estimates, this
% also stores subiteration results.
options.storeFP = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OTHER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If false, loads only the current subset of measurements to the selected
% device when using implementation 2
% Can be useful when dealing with large datasets, such as TOF data
% Unlike the below one (largeDim), this one has no restrictions on algorithms 
% or other features
% If you use listmode data or custom detector coordinates, this also
% affects the amount of coordinates transfered
% If false, will slow down computations but consume less memory, depending
% on the number of subsets
% Note: false will only have an effect when subsets are used!
% Default is true
options.loadTOF = true;

% This setting determines whether the high-dimensional scalable
% reconstruction is used (if set as true). Otherwise, the regular
% reconstruction is performed.
% NOTE: Currently the high-dimensional reconstructions are scaled
% differently than the regular ones
% Supports only some algorithms, such as FDK/FBP, PDHG and PKMA
% Default is false
options.largeDim = false;

% The number of power method iterations
% Applies only to PDHG and FISTA, and their variants as it is used to
% compute the Lipschitz values
% Default is 20, though 10 should be enough in most cases already
options.powerIterations = 20;

% If true, includes also the "diagonal" corners in the neighborhood in RDP
% By default, only the sides which the current voxel shares a side are
% included
% See https://omega-doc.readthedocs.io/en/latest/algorithms.html#rdp for
% details
% Default is false
options.RDPIncludeCorners = false;

% Whether to use L2 or L1 balls with proximal TV and TGV regularization
% Default is true
options.useL2Ball = true;

% If true, scales the relaxation parameters such that they are not too
% large
% Will slow down convergence
% Not recommended!
% Default is false
options.relaxationScaling = false;

% If true, will try to compute the relaxation parameters on-the-fly
% Not particularly reliable
% Not recommended!
% Default is false
options.computeRelaxationParameters = false;

% The window type for the filtering, both FDK and both the filtering-based
% preconditioners
% Default is Hamming window
% Available windows are: none, hamming, hann, blackman, nuttal, parzen, cosine, gaussian, and shepp-logan
% None gives the noisiest image
% Hamming and Hann should both give somewhat smoothed image
% Blackmann and Nuttal may, or may not, work
% Parzen gives more smoothed image than Hamming or Hann
% Cosine is slightly less smooth than Hamming or Hann
% Gaussian depends on the sigma value (see below), default 0.5 produces
% very smooth image
% Shepp-Logan smooths only very little, less than cosine, but slightly more
% than none
% The above are for FDK reconstructions, filtering behaves differently
% For filtering, there is much less differences
options.filterWindow = 'hamming';

% Related to above, when using Gaussian window this is the sigma value used
% in that filter
options.normalFilterSigma = 0.5;

% If true, uses images/textures in OpenCL/CUDA computations whenever
% possible
% If false, uses regular vectors/buffers
% On GPUs, it is recommended to keep this true (which is default value)
% On CPUs, you might get performance boost the other way around (false)
% Some OpenCL devices might not support images, in which case setting this
% to false should help
% Default is true
options.useImages = true;

% Use "fast math"
% If true, uses more inaccurate but faster math functions
% Default is true
options.useMAD = true;

% If true, TGV is only performed on each 2D slice even when using 3D inputs
% Default, the TGV computes everything in 3D
options.use2DTGV = false;

% Stores the primal-dual gaps during PDHG computations
% Default is false
options.storeResidual = false;

% Applies to implementation 4 ONLY
% Uses singe precision values if true, double precision otherwise
options.useSingles = true;

% Applies to implementation 2 ONLY
% Use FISTA acceleration with the selected algorithm
% If true, FISTA acceleration is used for the selected algorithm, even if
% the algorithm doesn't inherently use same acceleration scheme
options.FISTA_acceleration = false;


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

tStart = tic;
% pz is the reconstructed image volume
% recPar is a short struct of various reconstruction parameters used,
% useful when saving the reconstructed data with metadata
% classObject is the used class object in the reconstructions,
% essentially modified version of the input options-struct
% fp are the forward projections, if stored
% the primal-dual gap can be also be stored and is the variable after
% fp
[pz,recPar,classObject,fp] = reconstructions_mainCT(options);
tElapsed = toc(tStart);
disp(['Reconstruction process took ' num2str(tElapsed) ' seconds'])


volume3Dviewer(pz, [], [0 0 1])