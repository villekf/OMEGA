%% MATLAB/Octave code for reconstruction of any ray-tracing compatible data
% This example showcases how to use data that is not in (standard) PET,
% SPECT or CT format. Essentially the only things that are needed are the
% FOV size, the number of voxels in each direction and the source/detector
% coordinates. The example data is a cylindrical PET data but in reality it
% could be anything. This is a simplified example.
% Documentation: https://omega-doc.readthedocs.io/en/latest/customcoordinates.html
% NOTE: This example has no error checking!

clear


%%% Transaxial FOV size (mm), this is the length of the x (vertical) side
% of the FOV
options.FOVa_x = 300;

%%% Transaxial FOV size (mm), this is the length of the y (horizontal) side
% of the FOV
options.FOVa_y = options.FOVa_x;

%%% Axial FOV (mm)
options.axial_fov = floor(76.8 - 2.4/10);

% Note: Origin is assumed to be at the center. If this is not the case, you
% can shift it with options.oOffsetX, options.oOffsetY and options.oOffsetZ
% options.oOffsetX = 0;
% options.oOffsetY = 0;
% options.oOffsetZ = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% IMAGE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Reconstructed image pixel count (X/row-direction)
% NOTE: Non-square image sizes (X- and Y-direction) may not work
options.Nx = 256;

%%% Y/column-direction
options.Ny = 256;

%%% Z-direction (number of slices) (axial)
options.Nz = 63;

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

%%% Show status messages
% These are e.g. time elapsed on various functions and what steps have been
% completed. It is recommended to keep this 1. This can be at most 3.
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
% 1 = Reconstructions in MATLAB (projector in a MEX-file), uses matrices.
% (Slow and memory intensive)
% 2 = Matrix-free reconstruction with OpenCL/CUDA ArrayFire (Recommended)
% (Requires ArrayFire. Compiles with MinGW ONLY when ArrayFire was compiled
% with MinGW as well (cannot use the prebuilt binaries)).
% 3 = Multi-GPU/device matrix-free OpenCL (OSEM & MLEM only).
% 4 = Matrix-free reconstruction with OpenMP (parallel), standard C++
% 5 = Matrix-free reconstruction with OpenCL (parallel)
% See the documentation for more information:
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
options.use_device = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROJECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Type of projector to use for the geometric matrix
% 1 = Improved/accelerated Siddon's algorithm
% 2 = Orthogonal distance based ray tracer
% 3 = Volume of intersection based ray tracer
% 4 = Interpolation-based projector
% NOTE: You can mix and match most of the projectors. I.e. 41 will use
% interpolation-based projector for forward projection while improved
% Siddon is used for backprojection.
% See the documentation for more information:
% https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
options.projector_type = 1;

%%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Number of iterations (all reconstruction methods)
options.Niter = 1;

%%% Number of subsets (excluding subset_type = 6 and algorithms that do not
%%% support subsets)
options.subsets = 8;

%%% Subset type (n = subsets)
% 0 = Measurements are divided into n segments
% 1 = Every nth measurement is taken
% 3 = Measurements are selected randomly
options.subset_type = 3;

%%% Initial value for the reconstruction
options.x0 = ones(options.Nx, options.Ny, options.Nz);


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
%%% Ordered Subsets Expectation Maximization (OSEM) OR Maximum-Likelihood
%%% Expectation Maximization (MLEM) (if subsets = 1)
% Supported by all implementations
% Note: OSEM, MRAMLA, RAMLA, ROSEM, RBI, COSEM, ACOSEM, ECOSEM are not
% recommended for non-Poisson data
options.OSEM = false;

%%% LSQR
% Supported by implementations 1, 2, 4, and 5
options.LSQR = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAP-METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Any algorithm selected here will utilize any of the priors selected below
% this. Note that only one algorithm and prior combination is allowed! You
% can also use most of these algorithms without priors (such as PKMA or
% PDHG).
%%% Preconditioner Krasnoselskii-Mann algorithm (PKMA)
% Supported by implementations 1, 2, 4, and 5
options.PKMA = false;

%%% Primal-dual hybrid gradient (PDHG)
% Supported by implementations 1, 2, 4, and 5
options.PDHG = true;

% You can input other reconstruction parameters and priors as in the other
% examples

% Input the measurement data
load Cylindrical_PET_example_cylpet_example_new_sinograms_combined_static_200x168x703_span3.mat

options.SinM = raw_SinM;

% Set the coordinates for EACH measurement
% The format for options.x variable is
% [sourceX1;sourceY1;sourceZ1;detectorX1;detectorY1;detectorZ1;sourceX2;sourceY2;sourceZ2;detectorX2;...]
% where the number 1 refers to the first measurement, 2 to the second
% measurement, etc. X/Y/Z refers to the X/Y/Z (Cartesian) coordinates. This
% means that there should be SIX (6) coordinates for EACH measurement. Note
% that MATLAB is COLUMN-MAJOR, which means that column values are read
% first. As such, if you wish to use a matrix (which is fine) use
% 6xNumberOfMeasurements matrix. Vector format is recommended though.
load cylpet_example_det_coord.mat

sizeX = size(x,2);
x = repmat(x, 1, size(z_det,2));

z = repelem(z_det, 1, sizeX);

options.x = [x(1,:);x(2,:);z(1,:);x(3,:);x(4,:);z(2,:)];

clear x z

% Under normal situations, you would want to enable the "CT" mode which
% uses the length of intersection rather than the probability
% options.CT = true;

%% Reconstructions


tStart = tic;
% [pz, recPar, options, fp] = reconstructions_main(options);
% pz contains the 3D or 4D image
% recPar contains certain reconstruction parameters in a struct that can be
% easily included with the reconstruction
% options is the modified options struct (some modifications might be
% applied before reconstruction)
% fp are the stored forward projections, if options.storeFP = true
[pz, recPar, ~, fp] = reconstructions_main(options);
tElapsed = toc(tStart);
disp(['Reconstruction process took ' num2str(tElapsed) ' seconds'])

%     save([options.name '_reconstruction_' num2str(options.subsets) 'subsets_' num2str(options.Niter) 'iterations_' ...
%         num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '_' num2str(options.implementation) '.mat'], 'pz');

volume3Dviewer(pz, [], [0 0 1])