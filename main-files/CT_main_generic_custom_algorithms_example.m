%% MATLAB/Octave code for CBCT custom algorithm reconstruction
% This example contains a simplified example for custom algorithm
% reconstruction using TIFF projection CBCT data. Currently the support for
% some of the additional features is limited. The default configuration
% uses PDHG with multi-resolution, but CGLS is also included (use without
% multi-resolution!). PDHG should also work without multi-resolution 
% reconstruction.
% You can use the FIPS walnut data as an example data:
% https://zenodo.org/record/1254206

% Note that custom algorithm refers to your own algorithms and not the
% built-in algorithms. This example merely has the CGLS/PDHG algorithms
% shown as examples. The forward and/or backward projections of OMEGA are
% utilized for the computation of these algorithms. CGLS has been commented
% by default and this example uses multi-resolution reconstruction w/ PDHG.

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

%%% Detector pixel pitch/size (mm), column direction
% The size of the detector/distance between adjacent detectors
% NOTE: if you use binning, this value has to use the final binned
% dimensions
options.dPitchY = 0.05*options.binning;

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

%%% Transaxial FOV size (mm), this is the length of the x (vertical/row) side
% of the FOV
options.FOVa_x = 40.1;

%%% Transaxial FOV size (mm), this is the length of the y (horizontal/column) side
% of the FOV
options.FOVa_y = options.FOVa_x;

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
% options.nProjections = numel(options.angles);

% Flat value
% Needed for both linearized and Poisson-based data
% If omitted, the maximum value will be used automatically
options.flat = max(options.SinM(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% IMAGE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Reconstructed image pixel count (X/row-direction)
options.Nx = 280*2;

%%% Y/column-direction
options.Ny = 280*2;

%%% Z-direction (number of slices) (axial)
options.Nz = 280*2;

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

% When using extended FOV with custom reconstruction, only forward and
% backward projections automatically compute the multi-resolution volumes.
% The volume is assumed to be a cell-array, where the first cell element is
% the original FOV. The below, default, CGLS example does not use
% multi-resolution reconstruction. The commented PDHG example, however,
% does use multi-resolution imaging.

%%% Use projection extrapolation
options.useExtrapolation = false;

%%% Use extended FOV
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
options.implementation = 5;

% Applies to implementations 3 and 5 ONLY
%%% OpenCL platform used
% NOTE: Use OpenCL_device_info() to determine the platform numbers and
% their respective devices with implementations 3 or 5.
options.platform = 0;

% Applies to implementations 3 and 5 ONLY
%%% OpenCL device used
% NOTE: Use OpenCL_device_info() to determine the platform numbers and
% their respective devices with implementations 3 or 5.
% NOTE: if you switch devices then you need to run the below line
% (uncommented) as well:
% clear mex
options.use_device = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROJECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Type of projector to use for the geometric matrix
% 1 = Improved/accelerated Siddon's algorithm
% 4 = Interpolation-based projector
% 5 = Branchless distance-driven projector
% NOTE: You can mix and match most of the projectors. I.e. 41 will use
% interpolation-based projector for forward projection while improved
% Siddon is used for backprojection.
% NOTE 2: The below additional options apply also in hybrid cases as long
% as the other projector is the corresponding projector.
% See the docs for more information:
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


%%%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%
%%% Number of iterations (all reconstruction methods)
options.Niter = 2;

%%% Number of subsets
options.subsets = 20;

%%% Subset type (n = subsets)
% 8 = Use every nth sinogram
% 9 = Randomly select the full sinograms
% 10 = Use golden angle sampling to select the subsets (not recommended for
% PET)
% 11 = Use prime factor sampling to select the full sinograms
% Most of the time subset_type 8 is sufficient.
options.subset_type = 8;

%%% Initial value for the reconstruction
options.x0 = ones(options.Nx, options.Ny, options.Nz) * 1e-4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% CGLS

% % Required for CT data
% options.CT = true;
% % Linearize the input data
% raw_SinM = log(single(options.flat)./single(options.SinM));
% % CGLS doesn't support subsets
% options.subsets = 1;
% % Create the class object
% A = projectorClass(options);
% f = options.x0(:);
% r = raw_SinM(:);
% s = A' * r;
% p = s;
% gamma = norm(s)^2;
% tic
% for iter = 1 : options.Niter
    % q = A * p;
    % q = reshape(q, size(raw_SinM));
    % q = q(:);

    % alpha = gamma / norm(q)^2;
    % f_old = reshape(f, options.Nx,options.Ny,options.Nz);
    % f = f + alpha * p;
    % r = r - alpha * q;
    % s = A' * r;
    % gamma_ = norm(s)^2;
    % beta = gamma_ / gamma;
    % p = s + beta * p;
    % gamma = gamma_;
    % ff = reshape(f, options.Nx,options.Ny,options.Nz);
    % figure(1)
    % clim = [0 max(max(ff(:,:,200)))];
    % imagesc(flipud(ff(:,:,200)),clim)
    % axis image
% end
% toc

%% PDHG (multi-resolution support)

% Required for CT data
options.CT = true;
options.subsets = 10;
options.powerIterations = 10;
A = projectorClass(options);
raw_SinM = log(single(options.flat)./single(options.SinM));
if A.param.nMultiVolumes > 0
    u = cell(A.param.nMultiVolumes + 1, 1);
    for kk = 1 : A.param.nMultiVolumes + 1
        N = prod(A.param.N(kk));
        u{kk} = zeros(N,1,'single')*1e-4;
    end
else
    N = prod(A.param.N(1));
    u = zeros(N,1,'single')*1e-4;
end
u2 = u;
uO = u;
tau = 1./powerMethod(A);
sigma = 1;
theta = 1;
L1 = false;

% Data dual term
if options.subsets == 1
    M = size(raw_SinM,1) * size(raw_SinM,2) * A.nMeas(2);
    p = zeros(M,options.subsets,'single');
    pO = p;
else
    p = cell(options.subsets,1);
    for s = 1 : options.subsets
        M = A.nMeasSubset(s);
        p{s} = zeros(M,1,'single');
    end
    pO = p;
end
if options.subsets == 1
    raw_SinM = raw_SinM(:);
else
    raw_SinM = raw_SinM(:,:,A.index);
end

tic
for ii = 1:options.Niter
    for osa_iter = 1 : options.subsets
        uPrev = u;

        if options.subsets > 1
            A.subset = osa_iter;
        end
        temp2 = A * u2;
        if options.subsets == 1
            res = (temp2)-raw_SinM;
        else
            rr = raw_SinM(:,:,A.nMeas(osa_iter) + 1 : A.nMeas(osa_iter + 1));
            rr = rr(:);
            res = (temp2)-rr;
        end
        res = res(:);
        % Dual variable updates
        if options.subsets > 1
            p{osa_iter} = (pO{osa_iter} + sigma.*res);
        else
            p = (p + sigma.*res);
        end
        if L1
            if options.subsets == 1
                p = p ./ max(1, abs(p));
            else
                p{osa_iter} = p{osa_iter} ./ max(1, abs(p{osa_iter}));
            end
        else
            if options.subsets == 1
                p = p ./ (1 + sigma);
            else
                p{osa_iter} = p{osa_iter}./(1 + sigma);
            end
        end
        if osa_iter == 1
            if options.subsets == 1
                uUpd = (A'*((p)));
            else
                deltaZ = A'*(p{osa_iter} - pO{osa_iter});
            end
        else
            deltaZ = A'*(p{osa_iter} - pO{osa_iter});
        end
        if options.subsets > 1
            pO{osa_iter} = p{osa_iter};
            if A.param.nMultiVolumes > 0
                for kk = 1 : A.param.nMultiVolumes + 1
                    u{kk} = u{kk} + deltaZ{kk};
                    uO{kk} = u{kk} + options.subsets * deltaZ{kk};
                    u2{kk} = u2{kk} - tau(kk) .* uO{kk};
                end
            else
                u = u + deltaZ;
                uO = u + options.subsets * deltaZ;
                u2 = u2 - tau .* uO;
            end
        end
        if options.subsets == 1
            if A.param.nMultiVolumes > 0
                for kk = 1 : A.param.nMultiVolumes + 1
                    u{kk} = u{kk} - tau(kk) .* uUpd{kk};
                    u2{kk} = u{kk} + theta * (u{kk} - uPrev{kk});
                end
            else
                u = u - tau .* uUpd;
                u2 = u + theta * (u - uPrev);
            end
        else
            %     u2 = u;
        end
        if A.param.nMultiVolumes > 0
            if sum(u{1} - uPrev{1}) == 0
                break;
            end
            if sum(isnan(u{1}(:))) > 0
                break;
            end
        else
            if sum(u - uPrev) == 0
                break;
            end
            if sum(isnan(u(:))) > 0
                break;
            end
        end
    end

    % likelihood2(ii) = 1/2*norm((res(:)))^2;
    if A.param.nMultiVolumes > 0
        ff = reshape(u2{1}, A.param.Nx(1),A.param.Ny(1),A.param.Nz(1));
    else
        ff = reshape(u2, A.param.Nx,A.param.Ny,A.param.Nz);
    end
    figure(1)
    clim = [0 max(max(ff(:,:,255)))];
    imagesc(flipud(ff(:,:,255)),clim)
    axis image
    title([num2str(ii)])
end
toc
