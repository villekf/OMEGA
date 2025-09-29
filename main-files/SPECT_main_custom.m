%% MATLAB codes for SPECT reconstruction from custom sinogram data
% This example outlines the reconstruction of SPECT data.

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SCANNER PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Required for SPECT data
options.SPECT = true;

%%% Crystal thickness (mm)
options.cr_p = 9.525;

%%% Crystal width (mm)
options.dPitchX = 4.7952;
options.dPitchY = 4.7952;

%%% Scanner name
% Used for naming purposes (measurement data)
options.machine_name = 'Example';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% IMAGE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Reconstructed image pixel count
% NOTE: Non-square image sizes (X- and Y-direction) may not work
options.Nx = 128; % X-direction
options.Ny = 128; % Y-direction
options.Nz = 128; % Z-direction (number of axial slices)

%%% FOV size [mm]
% NOTE: Non-cubical voxels may not work
options.FOVa_x = options.dPitchX*64; % [mm], x-axis of FOV (transaxial)
options.FOVa_y = options.dPitchX*64; % [mm], y-axis of FOV (transaxial)
options.axial_fov = options.dPitchY*64; % [mm], z-axis of FOV (axial)

%%% Flip the image?
options.flipImageX = false;
options.flipImageY = false;
options.flipImageZ = false;

%%% Use back projection mask?
options.useMaskBP = false;

%%% Attenuation correction
options.attenuation_correction = false;

% Linear attenuation coefficients, size and dimensions should match FOV
options.vaimennus = []; 

%%% How much is the image rotated in degrees?
% NOTE: The rotation is done in the detector space (before reconstruction).
% Positive values perform the rotation in counterclockwise direction
options.offangle = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% PROJECTION DATA &&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gantry angles
options.angles = linspace(0, 90, 12)';

% Detector swivel angles
options.swivelAngles = options.angles+180;

% Distance between detector surface and FOV centre (origin)
options.radiusPerProj = 48*options.dPitchX*ones(size(options.angles)); 

% Initial value for the forward projection example
x0 = zeros(options.Nx, options.Ny, options.Nz);
x0(64, 64, 64) = 1;

% Projection images for backward projection example
y0 = zeros(128, 128, numel(options.angles), 'single'); 
y0(64, 64, :) = 1;

% Number of rows in a projection image
options.nRowsD = size(y0, 1);

% Number of columns in a projection image
options.nColsD = size(y0, 2);

% Number of projections
options.nProjections = size(y0, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% COLLIMATOR PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Collimator-detector response function (CDRF)
% For projector types 1, 2 and 6 you can either input either:
% 1. the collimator parameters (default) for an analytic solution for round (and hexagonal) holes (this may be unoptimal),
% 2. the standard deviations for both transaxial and axial directions or
% 3. the (Gaussian) PSF filter
% 4. the shifts of each ray traced 

% NOTE: For projector type 1 the CDRF is determined by
% options.rayShiftsDetector and options.rayShiftsSource defined in option 
% 4. These can also be calculated automatically when collimator parameters
% (1.) are input.

% NOTE: With projector_type == 2 (orthogonal distance projector), only the
% collimator parameters below are used for CDR calculation i.e. the
% collimator hole is assumed to be a circle. Thus only 1. below is
% supported with projector_type == 2

% 1. The collimator parameters (projector types 1, 2 and 6)
% Collimator hole length (mm)
options.colL = 24.05;
% Collimator hole radius (mm)
options.colR = 1.11/2;
% Distance from collimator to the detector (mm)
options.colD = 0;
% Intrinsic resolution (mm)
options.iR = 3.8;

% 2. If you have the standard deviations for transaxial (XY) and axial (Z)
% directions, you can input them here instead of the above values The
% dimensions need to be options.nProjections x options.Nx. Only for
% projector type 6.
% options.sigmaZ = repmat(1, options.nProjections, options.Nx); 
% options.sigmaXY = repmat(1, options.nProjections, options.Nx); % Transaxial standard deviation, (projector type 6)

% 3. You can input the filter for the CDRF directly. This should be of the
% size filterSizeXY x filterSizeZ x options.nProjections. Only for
% projector type 6.
% options.gFilter = ones(1, 1, options.Nx, options.nProjections);

% 4. For the Siddon ray tracer, the CDRF is defined by shifting the rays to
% the shape of the collimator hole. The below example is for random
% (uniform distribution) rays with one square collimator hole at the centre
% of each detector element. For 1 ray, the ray perpendicular to the
% detector element.
options.nRays = 1; % Number of rays traced per detector element
% options.rayShiftsDetector = options.colR*(2*rand(2*options.nRays, 1, 'single')-1)/options.dPitchX; % The relative shifts (dx1, dy1, dx2, dy2, ...) at the collimator-detector interface
% options.rayShiftsSource = options.colR*(2*rand(2*options.nRays, 1, 'single')-1)/options.dPitchY; % The relative shifts (dx1, dy1, dx2, dy2, ...) at the other end of the collimator

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
% completed. It is recommended to keep this 1.  Maximum value of 3 is
% supported.
options.verbose = 1;

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
% 4 = Matrix-free reconstruction with OpenMP (parallel), standard C++
% See the documentation for more information:
% https://omega-doc.readthedocs.io/en/latest/implementation.html
options.implementation = 2;

% Applies to implementation 2 ONLY
%%% OpenCL/CUDA device used
% NOTE: Use ArrayFire_OpenCL_device_info() to determine the device numbers
% with implementation 2.
options.use_device = 0;

% Implementation 2 ONLY
%%% Use CUDA
% Selecting this to true will use CUDA kernels/code instead of OpenCL. This
% only works if the CUDA code was successfully built.
options.use_CUDA = false;

% Implementation 2 ONLY
%%% Use CPU
% Selecting this to true will use CPU-based code instead of OpenCL or CUDA.
options.use_CPU = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROJECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Type of projector to use for the geometric matrix
% 1 = (Improved) Siddon ray-based projector
% 2 = Orthogonal distance ray tracing
% 6 = Rotation-based projector
% See the documentation on some details on the projectors:
% https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
options.projector_type = 2;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Projections
A = projectorClass(options);

x = A'*y0;
x = reshape(x, options.Nx, options.Ny, options.Nz);

y = A*x0;
y = reshape(y, options.nRowsD, options.nColsD, options.nProjections);


%% Plots
imshow(x(:, :, 64), []) % Backprojection
volume3Dviewer(y) % Forward projection