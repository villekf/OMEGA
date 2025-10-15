%% MATLAB codes for SPECT reconstruction from SIMIND output
% This example has been tested using SIMIND v8.0 and the respective
% tutorial sections 8(NEMA image quality phantom) and 9 (Brain CBF).
% Rather than built-in algorithms, this example highlights the projection 
% operators A*x and A'*y.

% Pinhole data can be generated with
% mpirun -np 5 simind_mpi nema nema_pinhole/fz:phantom/45:3/tr:11/tr:15/cc:ge-ph02/76:128/77:128/28:0.33/mp/55:1/53:1/42:-5/29:60/in:x22,5x/84:1

% Parallel-hole data can be generated with
% mpirun -np 5 simind_mpi nema nema_parallel/fz:phantom/45:3/tr:11/tr:15/cc:g8-luhr/76:128/77:128/28:0.33/mp/55:0/53:1/42:-5/29:60/in:x22,5x/84:1

clear
options.SPECT = true;

% All paths are input without file ending
options.fpath = '/path/to/data/nema_parallel_tot_w2'; % Main data (.a00 + .h00 files)
options.fpathCor = '/path/to/data/corfile'; % .cor file
options.fpathCT = '/path/to/data/nema_parallel'; % CT data (.hct + .ict files)
options.fpathScatterLower = '/path/to/data/nema_parallel_tot_w1'; % Lower scatter window data (.a00 + .h00 files)
options.fpathScatterUpper = '/path/to/data/nema_parallel_tot_w3'; % Upper scatter window data (.a00 + .h00 files)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SCANNER PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Crystal size (mm)
options.dPitchX = 3.3;
options.dPitchY = 3.3;

%%% Scanner name
% Used for naming purposes (measurement data)
options.machine_name = 'SIMIND';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% IMAGE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Reconstructed image pixel count
% NOTE: Square image sizes (X- and Y-direction) are recommended
options.Nx = 128; % X-direction
options.Ny = 128; % Y-direction
options.Nz = 128; % Z-direction (number of axial slices)

%%% FOV size [mm]
% NOTE: Cubical voxels are recommended
options.FOVa_x = options.dPitchX*128; % [mm], x-axis of FOV (transaxial)
options.FOVa_y = options.dPitchX*128; % [mm], y-axis of FOV (transaxial)
options.axial_fov = options.dPitchY*128; % [mm], z-axis of FOV (axial)

%%% Flip the image?
options.flipImageX = false;
options.flipImageY = false;
options.flipImageZ = false;

%%% Use back projection mask?
options.useMaskBP = false;

%%% How much is the image rotated in degrees?
% NOTE: The rotation is done in the detector space (before reconstruction).
% Positive values perform the rotation in counterclockwise direction
options.offangle = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CORRECTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% Attenuation correction %%%%%%%%%%%%%%%%%%%%%%%%%%
% Currently only a set of DICOM files is supported for the attenuation map.
options.attenuation_correction = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Scatter correction %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses linear interpolation between scatter windows. 
% See for example: 10.1371/journal.pone.0269542
options.scatter_correction = false;

%%%%%%%%%%%%%%%%%%%%%%%% Normalization correction %%%%%%%%%%%%%%%%%%%%%%%%%
% If set to true, normalization correction is applied to either the
% projection data or in the image reconstruction by using predefined
% normalization coefficients.
options.normalization_correction = false;
options.normalization = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%% Resolution recovery %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Collimator-detector response function (CDRF)
% For projector types 2 and 6 you can either input either:
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

% 1. The collimator parameters (projector types 1, 2 and 6).
% % 1.1 Pinhole collimator
%options.colD = 167.5; % Separation of collimator and detector
%options.colFxy = 0; % Focal distance XY, 0 for pinhole
%options.colFz = 0; % Focal distance Z, 0 for pinhole

% % 1.2 Parallel-hole collimator
options.colD = 0;
options.colFxy = Inf;
options.colFz = Inf;

% 2. If you have the standard deviations for transaxial (XY) and axial (Z)
% directions, you can input them here instead of the above values The
% dimensions need to be options.nProjections x options.Nx. Only for
% projector type 6.
% options.sigmaZ = repmat(1, options.nProjections, options.Nx); 
% options.sigmaXY = repmat(1, options.nProjections, options.Nx);

% 3. You can input the filter for the CDRF directly. This should be of the
% size filterSizeXY x filterSizeZ. Only for
% projector type 6.
% options.gFilter = ones(1, 1, options.Nx);

% 4. For the Siddon ray tracer, the CDRF is defined by shifting the rays to
% the shape of the collimator hole. The values of rayShiftsDetector and
% rayShiftsSource represent [shift1XY, shift1Z, shift2XY, ...] in mm. Size
% should be 2*nRays x nColsD x nRowsD x nProjections. If not input, values
% are calculated automatically.
options.nRays = 1; % Number of rays traced per detector element
% options.rayShiftsDetector = [];
% options.rayShiftsSource = [];


%%%%%%%%%%%%%%%%%%%% Corrections during reconstruction %%%%%%%%%%%%%%%%%%%%
% If set to true, all the corrections are performed during the
% reconstruction step, otherwise the corrections are performed to the
% sinogram/raw data before reconstruction (i.e. precorrected). I.e. this
% can be considered as e.g. normalization weighted reconstruction if
% normalization correction is applied.
% NOTE: Attenuation correction and resolution recovery are always performed
% during reconstruction regardless of the choice here.
options.corrections_during_reconstruction = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = loadSIMINDSPECTData(options);


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


%% Projections
A = projectorClass(options);

x = A'*y0;
x = reshape(x, options.Nx, options.Ny, options.Nz);

y = A*x0;
y = reshape(y, options.nRowsD, options.nColsD, options.nProjections);


%% Plots
volume3Dviewer(x) % Backprojection
volume3Dviewer(y) % Forward projection