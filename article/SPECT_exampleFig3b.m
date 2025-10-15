%% MATLAB codes for SPECT reconstruction from DICOM data
% This example outlines the reconstruction of SPECT data. In this case the
% data is Siemens Pro.specta projection data available at DOI
% 10.5281/zenodo.17315440

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('jaszczak_spectct_projection_data.mat')
options.SinM = projection_data;
options.angles = angular_position;
options.radiusPerProj = radial_position;
options.nRowsD = size(options.SinM, 1);
options.nColsD = size(options.SinM, 2);
options.nProjections = size(options.SinM, 3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% SCANNER PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Crystal thickness (mm)
options.cr_p = detector_thickness;

%%% Crystal width (mm)
options.dPitchX = pixel_spacing(1);
options.dPitchY = pixel_spacing(2);

%%% Scanner name
% Used for naming purposes (measurement data)
options.machine_name = 'Prospecta';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% IMAGE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Reconstructed image pixel count
% NOTE: Non-square image sizes (X- and Y-direction) may not work
options.Nx = 64; % X-direction
options.Ny = 64; % Y-direction
options.Nz = 128; % Z-direction (number of axial slices)

%%% FOV size [mm]
% NOTE: Non-cubical voxels may not work
options.FOVa_x = options.dPitchX*64; % [mm], x-axis of FOV (transaxial)
options.FOVa_y = options.dPitchX*64; % [mm], y-axis of FOV (transaxial)
options.axial_fov = options.dPitchY*128; % [mm], z-axis of FOV (axial)

%%% How much is the image rotated in degrees?
% NOTE: The rotation is done in the detector space (before reconstruction).
% Positive values perform the rotation in counterclockwise direction
options.offangle = 0;

% SPECT FOV spatial referencing
xLimits = options.FOVa_x * [-0.5, 0.5]; 
yLimits = options.FOVa_y * [-0.5, 0.5];
zLimits = options.axial_fov * [-0.5, 0.5];
refSPECT = imref3d([options.Nx, options.Ny, options.Nz], xLimits, yLimits, zLimits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CORRECTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% Attenuation correction %%%%%%%%%%%%%%%%%%%%%%%%%%
% Currently only a set of DICOM files is supported for the attenuation map.
options.attenuation_correction = true;

if options.attenuation_correction % Convert to LAC
    attenuation_map(attenuation_map < -1000) = -1000;
    MUvol = HU_to_mu(attenuation_map, 141);
    muAir = HU_to_mu(-1000, 141);
    tform = affinetform3d(eye(4)); % No scaling or rotation

    [MUvol, ~] = imwarp(MUvol, spatial_referencing_CT, tform, OutputView=refSPECT, FillValue=muAir, InterpolationMethod='linear');
    options.vaimennus = MUvol;
end

%%%%%%%%%%%%%%%%%%%%%%%% Normalization correction %%%%%%%%%%%%%%%%%%%%%%%%%
% If set to true, normalization correction is applied to either the
% projection data or in the image reconstruction by using predefined
% normalization coefficients.
options.normalization_correction = false;
options.normalization = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Scatter correction %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses linear interpolation between scatter windows. options.ScatterC{1} 
% contains the lower scatter window and options.ScatterC{2} contains the 
% upper scatter window (sizes equal options.SinM).
% See for example: 10.1371/journal.pone.0269542
options.scatter_correction = false;
options.ScatterC = {};
options.eWin = energy_window; % Main energy window: [lowerLimit upperLimit]
options.eWinL = []; % Lower energy window: [lowerLimit upperLimit]
options.eWinU = []; % Upper energy window: [lowerLimit upperLimit]

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

% 1. The collimator parameters (projector types 1, 2 and 6)
% Collimator hole length (mm)
options.colL = collimator_thickness;
% Collimator hole radius (mm)
options.colR = collimator_hole_radius;
% Distance from collimator to the detector (mm)
options.colD = 0;
% Intrinsic resolution (mm)
options.iR = detector_intrinsic_resolution;
% Focal distance (XY)
options.colFxy = Inf;
% Focal distance (Z)
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

%%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Number of iterations (all reconstruction methods)
options.Niter = 5;
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

%%% Number of subsets (all excluding MLEM and subset_type = 5)
options.subsets = 8;

%%% Subset type (n = subsets)
% 8 = Use every nth projection image (recommended for projector_type = 6)
% 9 = Randomly select the projection images
% 10 = Use golden angle sampling to select the subsets (not recommended for
% PET)
% 11 = Use prime factor sampling to select the projection images
options.subset_type = 8;

%%% Initial value for the reconstruction
options.x0 = ones(options.Nx, options.Ny, options.Nz);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION ALGORITHMS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% See examples in main-files folder for more algorithms.
options.OSEM = true;


%% Reconstructions
tStart = tic;
% pz contains the 3D or 4D image
% recPar contains certain reconstruction parameters in a struct that can be
% easily included with the reconstruction
% options is the modified options struct (some modifications might be
% applied before reconstruction)
% fp are the stored forward projections, if options.storeFP = true
[pz, recPar, ~, fp] = reconstructions_mainSPECT(options);
tElapsed = toc(tStart);
disp(['Reconstruction process took ' num2str(tElapsed) ' seconds'])

% save([options.name '_reconstruction_' num2str(options.subsets) 'subsets_' num2str(options.Niter) 'iterations_' ...
%     num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '.mat'], 'pz');

%% Plot
volume3Dviewer(pz, [], [0 0 1]) % Axial slices
volume3Dviewer(pz(:,:,17:80), [], [0 1 0]) % Sagittal slices