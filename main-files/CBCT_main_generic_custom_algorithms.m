%% MATLAB codes for CBCT custom algorithm reconstruction for Planmeca data
% This example contains a simplified example for custom algorithm
% reconstruction using Planmeca CBCT data. Currently the support for some
% of the additional features is limited. The default configuration uses
% CGLS without multi-resolution reconstruction, but PDHG with
% multi-resolution is included below CGLS, but has been commented. PDHG
% should also work without multi-resolution reconstruction.

clear
clear mex

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
options.binning = 1;

%%% Name of current datafile/examination
% This is used for naming purposes only
options.name = 'Planmeca_CT_data';

%%% Compute only the reconstructions
% If this file is run with this set to true, then the data load and
% sinogram formation steps are always skipped. Precomputation step is
% only performed if precompute_lor = true and precompute_all = true
% (below). Normalization coefficients are not computed even if selected.
options.only_reconstructions = false;

%%% Show status messages
% These are e.g. time elapsed on various functions and what steps have been
% completed. It is recommended to keep this true.
options.verbose = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This should be the full path to the metadata-file of the examination you
% wish to reconstruct.
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

% NOTE: If you want to reduce the number of projections, you need to do
% this manually as outlined below:
% lasku = 1;
% options.SinM = options.SinM(:,:,1:lasku:options.nProjections);
% options.angles = options.angles(1:lasku:options.nProjections);
% options.pitchRoll = options.pitchRoll(1:lasku:options.nProjections,:);
% options.x = options.x(1:lasku:options.nProjections,:);
% options.y = options.y(1:lasku:options.nProjections,:);
% options.z = options.z(1:lasku:options.nProjections,:);
% options.nProjections = numel(options.angles);

% NOTE: If you want to reduce the number of projections, you need to do
% this manually as outlined below:
lasku = 1;
options.SinM = options.SinM(:,:,1:lasku:options.nProjections);
options.angles = options.angles(1:lasku:options.nProjections);
options.pitchRoll = options.pitchRoll(1:lasku:options.nProjections,:);
options.x = options.x(1:lasku:options.nProjections,:);
options.y = options.y(1:lasku:options.nProjections,:);
options.z = options.z(1:lasku:options.nProjections,:);
options.nProjections = numel(options.angles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% IMAGE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The image size is taken from the conf file, but you can manually adjust
% these if desired
%%% Reconstructed image pixel size (X-direction)
options.Nx = 801;

%%% Y-direction
options.Ny = 801;

%%% Z-direction (number of slices) (axial)
options.Nz = 668;

% Use these two to rotate/flip the final image
%%% Flip the image (in vertical direction)?
options.flip_image = true;

%%% How much is the image rotated (radians)?
% The angle (in radians) on how much the image is rotated BEFORE
% reconstruction, i.e. the rotation is performed in the detector space.
options.offangle = (3*pi)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% When using extended FOV with custom reconstruction, only forward and
% backward projections automatically compute the multi-resolution volumes.
% The volume is assumed to be a cell-array, where the first cell element is
% the original FOV. The below, default, CGLS example does not use
% multi-resolution reconstruction. The commented PDHG example, however,
% does use multi-resolution imaging.

%%% Use projection extrapolation
options.useExtrapolation = false;

%%% Use extended FOV
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
% the multi-resolution volumes. Default is 1/4 of the original voxel size.
options.multiResolutionScale = 1/4;

% Performs the extrapolation and adjusts the image size accordingly
options = CTEFOVCorrection(options);

% Use offset-correction
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
% 2 = Matrix-free reconstruction with OpenCL/ArrayFire (Recommended)
% (Requires ArrayFire. Compiles with MinGW ONLY when ArrayFire was compiled
% with MinGW as well (cannot use the prebuilt binaries)).
% 3 = Multi-GPU/device matrix-free OpenCL (OSEM & MLEM only).
% 4 = Matrix-free reconstruction with OpenMP (parallel), standard C++
% 5 = Matrix-free reconstruction with OpenCL (parallel)
% See the doc for more information:
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
% NOTE: if you switch devices then you need to run the below line
% (uncommented) as well:
% clear mex
options.use_device = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROJECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Type of projector to use for the geometric matrix
% 1 = Improved/accelerated Siddon's algorithm
% 2 = Orthogonal distance based ray tracer
% 3 = Volume of intersection based ray tracer
% 4 = Interpolation-based projector
% 5 = Branchless distance-driven projector
% NOTE: You can mix and match most of the projectors. I.e. 41 will use
% interpolation-based projector for forward projection while improved
% Siddon is used for backprojection.
% See the doc for more information:
% https://omega-doc.readthedocs.io/en/latest/selectingprojector.html
options.projector_type = 4;

%%% Use mask
% The mask needs to be a binary mask (uint8 or logical) where 1 means that
% the pixel is included while 0 means it is skipped. Separate masks can be
% used for both forward and backward projection and either one or both can
% be utilized at the same time. E.g. if only backprojection mask is input,
% then only the voxels which have 1 in the mask are reconstructed.
% Currently the masks need to be a 2D image that is applied identically at
% each slice.
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
% value will be multiplied by the voxel size which means that 1 means that
% the interpolation length corresponds to a single voxel (transaxial)
% length. Larger values lead to faster computation but at the cost of
% accuracy. Recommended values are between [0.5 1].
options.dL = 0.5;


%%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Number of iterations (all reconstruction methods)
options.Niter = 10;

%%% Number of subsets (all excluding MLEM and subset_type = 5)
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
% 360/n_angles for 3D, see GitHub wiki for more information:
% https://github.com/villekf/OMEGA/wiki/Function-help#reconstruction-settings
% 7 = Form the subsets by using golden angle sampling
% 8 = Use every nth sinogram
% 9 = Randomly select the full sinograms
% 10 = Use golden angle sampling to select the subsets (not recommended for
% PET)
% 11 = Use prime factor sampling to select the full sinograms
% Most of the time subset_type 8 is sufficient.
options.subset_type = 8;

%%% Initial value for the reconstruction
options.x0 = ones(options.Nx, options.Ny, options.Nz,'single') * 1e-4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CGLS

% Required for CT data
options.CT = true;
% Linearize the input data
raw_SinM = log(single(options.flat)./single(options.SinM));
% CGLS doesn't support subsets
options.subsets = 1;
% Create the class object
A = projectorClass(options);
f = options.x0(:);
r = raw_SinM(:);
s = A' * r;
p = s;
gamma = norm(s)^2;
tic
for iter = 1 : options.Niter
    q = A * p;
    q = reshape(q, size(raw_SinM));
    q = q(:);

    alpha = gamma / norm(q)^2;
    f_old = reshape(f, options.Nx,options.Ny,options.Nz);
    f = f + alpha * p;
    r = r - alpha * q;
    s = A' * r;
    gamma_ = norm(s)^2;
    beta = gamma_ / gamma;
    p = s + beta * p;
    gamma = gamma_;
    ff = reshape(f, options.Nx,options.Ny,options.Nz);
    figure(1)
    clim = [0 max(max(ff(:,:,200)))];
    imagesc(flipud(ff(:,:,200)),clim)
    axis image
end
toc

% %% PDHG (multi-resolution support)
% 
% % Required for CT data
% options.CT = true;
% options.subsets = 10;
% options.powerIterations = 10;
% A = projectorClass(options);
% raw_SinM = log(single(options.flat)./single(options.SinM));
% if A.param.nMultiVolumes > 0
%     u = cell(A.param.nMultiVolumes + 1, 1);
%     for kk = 1 : A.param.nMultiVolumes + 1
%         N = prod(A.param.N(kk));
%         u{kk} = zeros(N,1,'single')*1e-4;
%     end
% else
%     N = prod(A.param.N(1));
%     u = zeros(N,1,'single')*1e-4;
% end
% u2 = u;
% uO = u;
% tau = 1./powerMethod(A);
% sigma = 1;
% theta = 1;
% L1 = false;
% 
% % Data dual term
% if options.subsets == 1
%     M = size(raw_SinM,1) * size(raw_SinM,2) * A.nMeas(2);
%     p = zeros(M,options.subsets,'single');
%     pO = p;
% else
%     p = cell(options.subsets,1);
%     for s = 1 : options.subsets
%         M = A.nMeasSubset(s);
%         p{s} = zeros(M,1,'single');
%     end
%     pO = p;
% end
% if options.subsets == 1
%     raw_SinM = raw_SinM(:);
% else
%     raw_SinM = raw_SinM(:,:,A.index);
% end
% 
% tic
% for ii = 1:options.Niter
%     for osa_iter = 1 : options.subsets
%         uPrev = u;
% 
%         if options.subsets > 1
%             A.subset = osa_iter;
%         end
%         temp2 = A * u2;
%         if options.subsets == 1
%             res = (temp2)-raw_SinM;
%         else
%             rr = raw_SinM(:,:,A.nMeas(osa_iter) + 1 : A.nMeas(osa_iter + 1));
%             rr = rr(:);
%             res = (temp2)-rr;
%         end
%         res = res(:);
%         % Dual variable updates
%         if options.subsets > 1
%             p{osa_iter} = (pO{osa_iter} + sigma.*res);
%         else
%             p = (p + sigma.*res);
%         end
%         if L1
%             if options.subsets == 1
%                 p = p ./ max(1, abs(p));
%             else
%                 p{osa_iter} = p{osa_iter} ./ max(1, abs(p{osa_iter}));
%             end
%         else
%             if options.subsets == 1
%                 p = p ./ (1 + sigma);
%             else
%                 p{osa_iter} = p{osa_iter}./(1 + sigma);
%             end
%         end
%         if osa_iter == 1
%             if options.subsets == 1
%                 uUpd = (A'*((p)));
%             else
%                 deltaZ = A'*(p{osa_iter} - pO{osa_iter});
%             end
%         else
%             deltaZ = A'*(p{osa_iter} - pO{osa_iter});
%         end
%         if options.subsets > 1
%             pO{osa_iter} = p{osa_iter};
%             if A.param.nMultiVolumes > 0
%                 for kk = 1 : A.param.nMultiVolumes + 1
%                     u{kk} = u{kk} + deltaZ{kk};
%                     uO{kk} = u{kk} + options.subsets * deltaZ{kk};
%                     u2{kk} = u2{kk} - tau(kk) .* uO{kk};
%                 end
%             else
%                 u = u + deltaZ;
%                 uO = u + options.subsets * deltaZ;
%                 u2 = u2 - tau .* uO;
%             end
%         end
%         if options.subsets == 1
%             if A.param.nMultiVolumes > 0
%                 for kk = 1 : A.param.nMultiVolumes + 1
%                     u{kk} = u{kk} - tau(kk) .* uUpd{kk};
%                     u2{kk} = u{kk} + theta * (u{kk} - uPrev{kk});
%                 end
%             else
%                 u = u - tau .* uUpd;
%                 u2 = u + theta * (u - uPrev);
%             end
%         else
%             %     u2 = u;
%         end
%         if A.param.nMultiVolumes > 0
%             if sum(u{1} - uPrev{1}) == 0
%                 break;
%             end
%             if sum(isnan(u{1}(:))) > 0
%                 break;
%             end
%         else
%             if sum(u - uPrev) == 0
%                 break;
%             end
%             if sum(isnan(u(:))) > 0
%                 break;
%             end
%         end
%     end
% 
%     % likelihood2(ii) = 1/2*norm((res(:)))^2;
%     if A.param.nMultiVolumes > 0
%         ff = reshape(u2{1}, A.param.Nx(1),A.param.Ny(1),A.param.Nz(1));
%     else
%         ff = reshape(u2, A.param.Nx,A.param.Ny,A.param.Nz);
%     end
%     figure(1)
%     clim = [0 max(max(ff(:,:,255)))];
%     imagesc(flipud(ff(:,:,255)),clim)
%     axis image
%     title([num2str(ii)])
% end
% toc