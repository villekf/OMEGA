%% MATLAB/Octave code for custom algorithm reconstruction of any ray-tracing compatible data
% This example showcases how to use data that is not in (standard) PET,
% SPECT or CT format. Essentially the only things that are needed are the
% FOV size, the number of voxels in each direction and the source/detector
% coordinates. The example data is a cylindrical PET data but in reality it
% could be anything. This is a simplified example and computes the
% reconstructions by using the built-in class object rather than the
% built-in algorithms. This example can be modified to compute your own
% algorithm.

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
options.implementation = 5;

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
options.use_device = 0;

% Applies to implementations 2, 3 and 5 ONLY
%%% Use 64-bit integer atomic functions
% If true, then 64-bit integer atomic functions (atomic add) will be used
% if they are supported by the selected device.
% Setting this to true will make computations faster on GPUs that support
% the functions, but might make results slightly less reliable due to
% floating point rounding. Recommended for OpenCL GPUs. Not recommended for
% CUDA. Doesn't apply for CPU.
options.use_64bit_atomics = true;

% Implementation 2 ONLY
%%% Use CUDA
% Selecting this to true will use CUDA kernels/code instead of OpenCL. This
% only works if the CUDA code was successfully built.
options.use_CUDA = false;

% Implementation 2 ONLY
%%% Use CPU
% Selecting this to true will use CPU-based code instead of OpenCL or CUDA.
% Some features are not supported by CPU such as projector_type 4 and 5.
options.use_CPU = false;

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

%%% Use mask
% The mask needs to be a binary mask (uint8 or logical) where 1 means that
% the pixel is included while 0 means it is skipped. Separate masks can be
% used for both forward and backward projection and either one or both can
% be utilized at the same time. E.g. if only backprojection mask is input,
% then only the voxels which have 1 in the mask are reconstructed.
% Currently the masks need to be a 2D image that is applied identically at
% each slice/sinogram/projection.
% Forward projection mask
% If nonempty, the mask will be applied. If empty, or completely omitted, no
% mask will be considered.
% options.maskFP = true(options.Ndist,options.Nang);
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
% value will be multiplied by the voxel size which means that the
% interpolation length of 1 corresponds to a single voxel (transaxial) 
% length. Larger values lead to faster computation but at the cost of
% accuracy. Recommended values are between [0.5 1], though values up to 2
% should be fine.
options.dL = 0.5;

%%% Use point spread function (PSF) blurring
% Applies PSF blurring through convolution to the image space. This is the
% same as multiplying the geometric matrix with an image blurring matrix.
options.use_psf = false;

% FWHM of the Gaussian used in PSF blurring in all three dimensions (X/Y/Z)
options.FWHM = [2.4 2.4 2.4];

% Orthogonal ray tracer (projector_type = 2 only)
%%% The 2D (XY) width of the "strip/tube" where the orthogonal distances are
% included. If tube_width_z below is non-zero, then this value is ignored.
options.tube_width_xy = 2.4;

% Orthogonal ray tracer (projector_type = 2 only)
%%% The 3D (Z) width of the "tube" where the orthogonal distances are
% included. If set to 0, then the 2D orthogonal ray tracer is used. If this
% value is non-zero then the above value is IGNORED.
options.tube_width_z = 2.4;

% Volume ray tracer (projector_type = 3 only)
%%% Radius of the tube-of-response (cylinder)
% The radius of the cylinder that approximates the tube-of-response.
% Default uses circle size that is just large enough to fit one detector
% crystal
options.tube_radius = sqrt(2) * (2.4 / 2);

% Volume ray tracer (projector_type = 3 only)
%%% Relative size of the voxel (sphere)
% In volume ray tracer, the voxels are modeled as spheres. This value
% specifies the relative radius of the sphere such that with 1 the sphere
% is just large enoough to encompass an entire cubic voxel, i.e. the
% corners of the cubic voxel intersect with the sphere shell. Larger values
% create larger spheres, while smaller values create smaller spheres.
options.voxel_radius = 1;

% Siddon (projector_type = 1 only)
%%% Number of rays
% Number of rays used per detector if projector_type = 1 (i.e. Improved
% Siddon is used).
% Number of rays in transaxial direction
options.n_rays_transaxial = 1;
% Number of rays in axial direction
options.n_rays_axial = 1;

%%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Number of iterations (all reconstruction methods)
options.Niter = 1;

%%% Save specific intermediate iterations
% You can specify the intermediate iterations you wish to save here. Note
% that this uses zero-based indexing, i.e. 0 is the first iteration (not
% the initial value). By default only the last iteration is saved.
% Note: Subiterations cannot be saved!
options.saveNIter = [];
% Alternatively you can save ALL intermediate iterations by setting the
% below to true and uncommenting it
% Note: Only affects full iterations (epochs)
% options.save_iter = false;

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

%% Class example (OSEM)

% Here is an example of how to obtain the same results as above by using a
% specific MATLAB class. This is a bit more simplified from above and also
% allows more easily to use other properties files (such as
% Inveon_PET_main.m). PSF blurring will be performed automatically if it
% has been selected.

% Load data
load Cylindrical_PET_example_cylpet_example_new_sinograms_combined_static_200x168x703_span3.mat
% When using custom coordinates, it is important to store the measurements
% in options.SinM BEFORE constructing the class object
if options.implementation == 1
    options.SinM = double(raw_SinM(:));
else
    options.SinM = single(raw_SinM(:));
end


% Construct the forward and backward projections object (you need to rerun
% this if you make any changes to the system):
A = projectorClass(options);

% Important if you use subsets!
options.SinM = options.SinM(A.index);

% Compute the OSEM-estimate for the specified number of iterations and
% subsets (use 1 subset if you want to use a method that doesn't use
% subsets):
f = options.x0(:);
sens = cell(options.subsets, 1);
for iter = 1 : options.Niter
    for osa_iter = 1 : options.subsets
        A.subset = osa_iter;
        % If you use implementation 1, you can form the system matrix for
        % the current subset with H = formMatrix(A,osa_iter);
        % The matrix is a regular sparse matrix.
        % Note that the system matrix is the TRANSPOSE of the matrix
        % The forward projection is stored in y
        y = A * f;
        % The backprojection is stored in x
        if iter == 1
            % Sensitivity image can be computed during the first iteration
            % Computed ONLY if the second output is present
            [x, S] = backwardProject(A, options.SinM(A.nMeas(osa_iter) + 1:A.nMeas(osa_iter+1)) ./ (y + A.param.epps), osa_iter);
            sens{osa_iter} = S;
        else
            % NOTE: Implementations 3/5 and 4 include randoms/scatter
            % correction to y automatically. With implementation 1
            % randoms/scatter must be added manually.
            x = A' * (options.SinM(A.nMeas(osa_iter) + 1:A.nMeas(osa_iter+1)) ./ (y +  A.param.epps));
        end
        f = (f ./ (sens{osa_iter} +  A.param.epps)) .* (x +  A.param.epps);
    end
end
ff = reshape(f, options.Nx,options.Ny,options.Nz);
volume3Dviewer(ff, [], [0 0 1])

%% Class example (MLEM)

% Same as above but without subsets

options.subsets = 1;

% Load data
load Cylindrical_PET_example_cylpet_example_new_sinograms_combined_static_200x168x703_span3.mat
% When using custom coordinates, it is important to store the measurements
% in options.SinM BEFORE constructing the class object
if options.implementation == 1
    options.SinM = double(raw_SinM(:));
else
    options.SinM = single(raw_SinM(:));
end


% Construct the forward and backward projections object (you need to rerun
% this if you make any changes to the system):
A = projectorClass(options);

% Compute the OSEM-estimate for the specified number of iterations and
% subsets (use 1 subset if you want to use a method that doesn't use
% subsets):
f = options.x0(:);
for iter = 1 : options.Niter
    % If you use implementation 1, you can form the system matrix with H =
    % formMatrix(A);
    % The matrix is a regular sparse matrix.
    % Note that the system matrix is the TRANSPOSE of the matrix
    % The forward projection is stored in y
    y = A * f;
    % The backprojection is stored in x
    if iter == 1
        % Sensitivity image can be computed during the first iteration
        % Computed ONLY if the second output is present
        [x, S] = backwardProject(A, options.SinM ./ (y + A.param.epps));
    else
        % NOTE: Implementations 3/5 and 4 include randoms/scatter
        % correction to y automatically. With implementation 1
        % randoms/scatter must be added manually.
        x = A' * (options.SinM ./ (y +  A.param.epps));
    end
    f = (f ./ (S +  A.param.epps)) .* (x +  A.param.epps);
end
ff = reshape(f, options.Nx,options.Ny,options.Nz);
volume3Dviewer(ff, [], [0 0 1])