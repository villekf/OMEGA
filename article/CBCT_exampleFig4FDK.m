%% MATLAB/Octave codes for CBCT reconstruction
% This example file computes the Figure 4 of the OMEGA V2 article with FDK.
% DOI will be added later.
% The input data Planmeca_VisoG7_100kV_80mAs_500proj_kneePhantom.mat 
% has to be in MATLAB/Octave path!
% Used data available from: https://doi.org/10.5281/zenodo.12722386
% Implementation 2 is used by default, but the device number can be set
% below (options.use_device)
clear

% Set the path (to the folder) of the above mat-file to here (if NOT in
% MATLAB/Octave path):
path = ''

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
options.verbose = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fullfile(path, 'Planmeca_VisoG7_100kV_80mAs_500proj_kneePhantom.mat'))

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
% If omitted, will use the maximum value from the input
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
% As source-detector pairs
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
options.useEFOV = false;

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
% 2 = Matrix-free reconstruction with OpenCL/CUDA (Recommended)
% (Requires ArrayFire).
% See the doc for more information:
% https://omega-doc.readthedocs.io/en/latest/implementation.html
options.implementation = 2;

% Applies to implementations 2, 3 and 5 ONLY
%%% OpenCL/CUDA device used 
% NOTE: Use ArrayFire_OpenCL_device_info() to determine the device numbers
% with implementation 2.
options.use_device = 0;

% Implementation 2 ONLY
%%% Use CUDA
% Selecting this to true will use CUDA kernels/code instead of OpenCL. This
% only works if the CUDA code was successfully built. This is recommended
% if you have CUDA-capable device.
options.use_CUDA = checkCUDA(options.use_device);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROJECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Type of projector to use for the geometric matrix
% 1 = Improved/accelerated Siddon's algorithm
% 4 = Interpolation-based projector (ray- and voxel-based)
% NOTE: You can mix and match most of the projectors. I.e. 45 will use
% interpolation-based projector for forward projection while branchless
% distance-driven is used for backprojection
% NOTE 2: The below additional options apply also in hybrid cases as long
% as the other projector is the corresponding projector.
% See the doc for more information:
% https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
% NOTE 3: FDK currently only supports projector type 4 as the backprojector
% and types 1 or 4 as the forward projectors. Thus only projector_types 14
% and 4 are valid for FDK.
options.projector_type = 4;

%%% Interpolation length (projector type = 4 only)
% This specifies the length after which the interpolation takes place. This
% value will be multiplied by the voxel size which means that 1 is
% the interpolation length corresponding to a single voxel (transaxial)
% length. Larger values lead to faster computation but at the cost of
% accuracy. Recommended values are between [0.5 1].
options.dL = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION ALGORITHMS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction algorithms to use (choose only one algorithm)

%%% Feldkamp-Davis-Kress (FDK)
% Supported by implementation 2
options.FDK = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FDK PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use Parker weights
% By default, FDK does NOT use Parker weights. However, if you are
% reconstructing a scan with less than 2pi coverage, it is HIGHLY
% recommended to set the Parker weights true as they are here. Note that
% you can also set an optional Parker weight value which affects the Parker
% weights as described in https://doi.org/10.1118/1.1450132
options.useParkerWeights = true;
% Specify the Parker weight below. The default value is 0.25 as in the
% above article, but 1 gives the "original" Parker weights.
options.ParkerWeight = 0.25;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
