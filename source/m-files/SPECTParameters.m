function options = SPECTParameters(options)
%SPECTParameters Computes the necessary variables for projector types 2 and 6
%(SPECT)
%   Computes the PSF standard deviations for the current SPECT collimator
%   if omitted. Also computes the interaction planes with the rotated
%   projector and the projection angles.
%
% Copyright (C) 2022-2025 Ville-Veikko Wettenhovi, Matti Kortelainen, Niilo Saarlemo
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ismember(options.projector_type, [1, 11, 12, 16, 2, 21, 22, 26, 61, 62]) % Ray tracing projectors, including the ray tracing side of the hybrid projectors
    % The central source point is placed at the transaxial midpoint, while axial ray shifts use the axial midpoint.
    if options.colD + 0.5 * options.colLxy <= 0 || options.colD + 0.5 * options.colLz <= 0
        error(['SPECT projector types 1 and 2 require a positive collimator length and/or collimator-to-detector distance: ' ...
            'colD + 0.5 * colLxy and/or colD + 0.5 * colLz would result in zero-length rays. ' ...
            'Set options.colLxy, options.colLz, and options.colD to positive collimator geometry values.'])
    end
end

if ismember(options.projector_type, [1, 11, 12, 2, 21, 22]) % Collimator modelling, ray tracing projectors
    nRays = options.n_rays_transaxial * options.n_rays_axial;
    if numel(options.rayShiftsDetector) == 0
        options.rayShiftsDetector = [0; 0];
        options.rayShiftsDetector = repmat(options.rayShiftsDetector, [nRays, options.nRowsD, options.nColsD, options.nHeads]);

        if options.colFxy == 0 && options.colFz == 0 % Pinhole collimator
            dx = linspace(-(options.nRowsD/2-0.5)*options.dPitchX, (options.nRowsD/2-0.5)*options.dPitchX, options.nRowsD);
            dy = linspace(-(options.nColsD/2-0.5)*options.dPitchY, (options.nColsD/2-0.5)*options.dPitchY, options.nColsD);
            
            for ii = 1:options.nRowsD
                for jj = 1:options.nColsD
                    for kk = 1:nRays
                        options.rayShiftsDetector(2*(kk-1)+1, ii, jj, :) = -dx(ii);
                        options.rayShiftsDetector(2*(kk-1)+2, ii, jj, :) = -dy(jj);
                    end
                end
            end
        end
    end
    if numel(options.rayShiftsSource) == 0
        options.rayShiftsSource = [0; 0];
        options.rayShiftsSource = repmat(options.rayShiftsSource, [nRays, options.nRowsD, options.nColsD, options.nHeads]);

        if nRays > 1 % Multiray shifts
            [tmp_x, tmp_y] = ndgrid(linspace(-0.5, 0.5, options.n_rays_transaxial), ...
                linspace(-0.5, 0.5, options.n_rays_axial));
            if options.colFxy == 0 && options.colFz == 0 % Pinhole collimator
                tmp_x = options.dPitchX * tmp_x;
                tmp_y = options.dPitchY * tmp_y;
            elseif ismember(options.colFxy, [-Inf, Inf]) && ismember(options.colFz, [-Inf, Inf]) % Parallel-hole collimator
                tmp_x = 2 * options.colR * tmp_x;
                tmp_y = 2 * options.colR * tmp_y;
            end

            tmp_shift = reshape([tmp_x(:), tmp_y(:)].', 1, [])';

            for kk = 1:nRays
                options.rayShiftsSource(2*(kk-1)+1,:,:,:) = tmp_shift(2*(kk-1)+1);
                options.rayShiftsSource(2*(kk-1)+2,:,:,:) = tmp_shift(2*(kk-1)+2);
            end
        end
    end
    if ismember(options.projector_type, [1, 11, 12, 21])
        % Use the XY source point as the shared reference; axial source shifts scale against colLz below.
        referenceLength = options.colD + 0.5 * options.colLxy;
        lengthXY = options.colD + 0.5 * options.colLxy;
        lengthZ = options.colD + 0.5 * options.colLz;
        % The first entry of every shift pair moves the ray along the detector panel (transaxial, XY) and the second one along the
        % rotation axis (axial, Z) (see getDetectorCoordinatesSPECT in general_opencl_functions.h), hence the first entries use
        % the transaxial septa length and the second ones the axial septa length.
        if referenceLength ~= lengthXY
            detectorShiftXY = options.rayShiftsDetector(1:2:end,:,:,:);
            options.rayShiftsSource(1:2:end,:,:,:) = detectorShiftXY + ...
                (options.rayShiftsSource(1:2:end,:,:,:) - detectorShiftXY) * (referenceLength / lengthXY);
        end
        if referenceLength ~= lengthZ
            detectorShiftZ = options.rayShiftsDetector(2:2:end,:,:,:);
            options.rayShiftsSource(2:2:end,:,:,:) = detectorShiftZ + ...
                (options.rayShiftsSource(2:2:end,:,:,:) - detectorShiftZ) * (referenceLength / lengthZ);
        end
    end
    if ismember(options.implementation, [2, 3, 5])
        options.rayShiftsDetector = single(options.rayShiftsDetector(:));
        options.rayShiftsSource = single(options.rayShiftsSource(:));
    end
end
if ismember(options.projector_type, [12, 2, 21, 22]) % Orthogonal distance ray tracer
    % Frey, E. C., & Tsui, B. M. W. (n.d.). Collimator-Detector Response Compensation in SPECT. Quantitative Analysis in Nuclear Medicine Imaging, 141–166. doi:10.1007/0-387-25444-7_5 
    % ODRT currently models an isotropic scalar response; use the transaxial length until the kernel supports both axes.
    if (~isfield(options,'coneOfResponseStdCoeffA') || options.coneOfResponseStdCoeffA < 0)
        options.coneOfResponseStdCoeffA = 2*options.colR/options.colLxy; % See equation (6) of book chapter
    end
    if (~isfield(options,'coneOfResponseStdCoeffB') || options.coneOfResponseStdCoeffB < 0)
        options.coneOfResponseStdCoeffB = 2*options.colR/options.colLxy*(options.colLxy+options.colD+options.cr_p/2);
    end
    if (~isfield(options,'coneOfResponseStdCoeffC') || options.coneOfResponseStdCoeffC < 0)
        options.coneOfResponseStdCoeffC = options.iR;
    end
    % Now the collimator response FWHM is sqrt((az+b)^2+c^2) where z is distance along detector element normal vector
end
if ismember(options.projector_type, [2, 12, 21, 22])
    % Optional precomputed SPECT ODRT lookup. The table axes are ray-local
    % (u, v, depth), stored in MATLAB/Fortran order (u is contiguous).
    % Lateral coordinate zero lies at sample index (N-1)/2, so even-sized
    % lateral axes are centered between their two middle samples. Depth index
    % zero is at the original ray start before ellipse clipping; the kernel
    % restores the per-ray distance removed by clipping. gFilterSpacing is
    % (du, dv, dd) in mm.
    % An empty filter preserves the analytic CoR kernel. Both type 6 and ODRT
    % model the detector response, but their arrays are not drop-in: type 6
    % uses a depth-shifted image-grid convolution kernel, while ODRT samples a
    % ray-local (u,v,depth) table with explicit mm spacing and direct weights.
    % The public gFilter input is shared for pure modes. Hybrids 26/62 need
    % both representations, so they retain type-6 gFilter semantics and use
    % the analytic ODRT response until both inputs can be supplied/converted.
    if ~isfield(options, 'gFilter') || isempty(options.gFilter)
        options.gFilter = single([]);
        options.gFilterSpacing = single([]);
        options.gFilterCustom = false;
        options.gFilterNu = uint32(0);
        options.gFilterNv = uint32(0);
        options.gFilterNd = uint32(0);
    else
        if ~isnumeric(options.gFilter) || ~isreal(options.gFilter) || ndims(options.gFilter) > 3
            error('ODRT gFilter must be a real numeric 2-D or 3-D array with axes (u, v, depth).');
        end
        psf = single(options.gFilter);
        psfSize = [size(psf, 1), size(psf, 2), size(psf, 3)];
        if any(psfSize < 1) || any(double(psfSize) > double(intmax('uint32')))
            error('ODRT gFilter dimensions must be positive and fit in uint32.');
        end
        if any(~isfinite(psf(:))) || any(psf(:) < 0) || ~any(psf(:) > 0)
            error('ODRT gFilter must contain finite, non-negative weights and at least one positive value.');
        end
        if ~isfield(options, 'gFilterSpacing') || ~isnumeric(options.gFilterSpacing) || ...
                ~isreal(options.gFilterSpacing) || numel(options.gFilterSpacing) ~= 3
            error('gFilterSpacing must contain three positive finite values (du, dv, dd) in mm.');
        end
        spacing = single(options.gFilterSpacing(:).');
        if any(~isfinite(spacing)) || any(spacing <= 0)
            error('gFilterSpacing must contain three positive finite values representable as single precision.');
        end
        options.gFilter = psf;
        options.gFilterSpacing = spacing;
        options.gFilterCustom = true;
        options.gFilterNu = uint32(psfSize(1));
        options.gFilterNv = uint32(psfSize(2));
        options.gFilterNd = uint32(psfSize(3));
    end
elseif ismember(options.projector_type, [6, 16, 26, 61, 62, 66])
    % Clear only ray-local ODRT metadata if a parameter struct is reused.
    % Keep the established type-6 gFilter contents and interpretation.
    options.gFilterSpacing = single([]);
    options.gFilterCustom = false;
    options.gFilterNu = uint32(0);
    options.gFilterNv = uint32(0);
    options.gFilterNd = uint32(0);
end
if options.projector_type == 6
    % Nx/Ny/Nz/dx/dy/dz may be vectors when multi-resolution volumes are in
    % use (options.useMultiResolutionVolumes). The C++/MEX side only
    % consumes a single options.gFilter/blurPlanes shared by all volumes
    % (see mfunctions.h/libHeader.h), so the CDRF and blur-plane geometry
    % below are always computed from the main (first) volume's geometry,
    % matching the previous scalar-only behavior bit-for-bit when Nx/dx
    % etc. are scalars.
    panelTilt = options.swivelAngles - options.angles + 180;

    % Effective septa lengths of the axial and transaxial collimator responses (resolved above)
    colLxy = options.colLxy;
    colLz = options.colLz;

    % Sub-voxel accurate geometry of the detector panel in the rotated coordinate system.
    % blurPlanesF is the distance from the collimator surface to the first image plane and
    % blurPlanes2F the lateral shift of the panel with respect to the rotation axis, both in
    % voxels. With detector heads that can swivel, panelTilt is different for every projection
    % and neither of these is, in general, an integer number of voxels. The geometry is given
    % with respect to the main volume, which the projector converts to the other volumes.
    options.blurPlanesF = (options.FOVa_x(1)/2 - (options.radiusPerProj .* cosd(panelTilt) - options.CORtoDetectorSurface) + options.colLxy + options.colD) / options.dx(1);
    options.blurPlanes2F = options.radiusPerProj .* sind(panelTilt) / options.dx(1);

    % Coefficients of the collimator-detector response used by the on-the-fly CDRF (and always
    % passed to the projector): FWHM_axis(z) = sqrt((A_axis*z + B_axis)^2 + C^2), where z is the
    % distance from the collimator surface in mm
    Az = 2*options.colR/colLz;
    Axy = 2*options.colR/colLxy;
    Bz = Az*(colLz+options.colD+options.cr_p/2);
    Bxy = Axy*(colLxy+options.colD+options.cr_p/2);
    options.CDRFCoeffs6 = single([Az Bz Axy Bxy options.iR]);

    if ~isfield(options, 'CDRFOnTheFly') || numel(options.CDRFOnTheFly) == 0
        % The on-the-fly response is the default with implementation 2 whenever the collimator
        % parameters are used, i.e. the CDRF is not input manually (gFilter/sigmaZ). It is
        % considerably faster than the precomputed filter bank (with the fused GPU pipeline several
        % times so) and more accurate: the bank can only be shifted by an integer number of voxels,
        % which is accurate only when the distance from the detector to the FOV is the same for every
        % projection, and it is sampled at the voxel size of the main volume, which the
        % multi-resolution volumes do not share, so those always need the on-the-fly response.
        options.CDRFOnTheFly = (~isfield(options, 'gFilter') || isempty(options.gFilter)) && ~isfield(options, 'sigmaZ') && ...
            options.implementation == 2;
    end
    options.CDRFOnTheFly = logical(options.CDRFOnTheFly);
    if options.CDRFOnTheFly && options.implementation ~= 2
        warning('On-the-fly computation of the collimator-detector response is only supported by implementation 2. Using a precomputed filter bank instead.')
        options.CDRFOnTheFly = false;
    end
    % The on-the-fly response uses the coefficients above (Az = 2 * colR / colLz, Axy = 2 * colR / colLxy), which need positive septa lengths
    if options.CDRFOnTheFly && (colLxy <= 0 || colLz <= 0)
        error(['SPECT projector type 6 computes the collimator-detector response from the collimator hole lengths: ' ...
            'options.colLxy and options.colLz have to be positive, or sigmaZ/sigmaXY/gFilter have to be input.'])
    end
    if ~options.CDRFOnTheFly && options.nMultiVolumes > 0
        error(['Multi-resolution reconstruction with projector type 6 requires the collimator-detector response to be computed on-the-fly, ' ...
            'since a precomputed filter bank only applies to the voxel size of the main volume. Remove the manually input filter bank (gFilter) ' ...
            'and sigmaZ, or disable multi-resolution volumes.'])
    end

    % The model of the collimator response: 'gaussian' (0, default) or 'square' (1, square holes). The number is what the projector
    % reads (cdrfModel6). Only the on-the-fly response has a model, a manual filter bank (gFilter) or sigmaZ is used as it is
    cdrfModelNames = {'gaussian', 'square'};
    if isempty(options.CDRFModel)
        options.CDRFModel = 'gaussian';
    end
    if ~(ischar(options.CDRFModel) || (isstring(options.CDRFModel) && isscalar(options.CDRFModel))) || ...
            ~any(strcmpi(strtrim(char(options.CDRFModel)), cdrfModelNames))
        error(['Unknown options.CDRFModel. The supported models are: ' strjoin(cellfun(@(x) ['''' x ''''], cdrfModelNames, 'UniformOutput', false), ', ') '.'])
    end
    cdrfModelInd = find(strcmpi(strtrim(char(options.CDRFModel)), cdrfModelNames), 1) - 1;
    if cdrfModelInd > 0 && ~options.CDRFOnTheFly
        error(['options.CDRFModel = ''' char(options.CDRFModel) ''' requires the on-the-fly collimator-detector response (options.CDRFOnTheFly, implementation 2, ' ...
            'no manually input gFilter or sigmaZ).'])
    end
    options.cdrfModel6 = uint32(cdrfModelInd);

    % How the backprojection rotates a view back into the image: 'adjoint' (0, default) or 'fast' (1). The number is what the projector
    % reads (bpRot6)
    bpRotationNames = {'adjoint', 'fast'};
    if isempty(options.type6BPRotation)
        options.type6BPRotation = 'adjoint';
    end
    if ~(ischar(options.type6BPRotation) || (isstring(options.type6BPRotation) && isscalar(options.type6BPRotation))) || ...
            ~any(strcmpi(strtrim(char(options.type6BPRotation)), bpRotationNames))
        error(['Unknown options.type6BPRotation. The supported values are: ' strjoin(cellfun(@(x) ['''' x ''''], bpRotationNames, 'UniformOutput', false), ', ') '.'])
    end
    options.bpRot6 = uint32(find(strcmpi(strtrim(char(options.type6BPRotation)), bpRotationNames), 1) - 1);

    DistanceToFirstRow = 0.5*options.dx(1);
    Distances = repmat(DistanceToFirstRow,1,options.Nx(1)*4)+repmat((0:double(options.Nx(1)*4)-1)*double(options.dx(1)),length(DistanceToFirstRow),1);
    if options.CDRFOnTheFly
        % The filters are computed separately for every projection inside the projector
        options.gFilter = single([]);
    elseif (~isfield(options,'gFilter'))
        if ~isfield(options, 'sigmaZ')
            % Axial and transaxial collimator responses use their own septa
            % lengths (colLz/colLxy), matching the Python projector type 6 CDRF computation.
            if colLxy <= 0 || colLz <= 0
                error(['SPECT projector type 6 computes the collimator-detector response from the collimator hole lengths: ' ...
                    'options.colLxy and options.colLz have to be positive, or sigmaZ/sigmaXY/gFilter have to be input.'])
            end
            Rg_z = 2*options.colR*(Distances+options.cr_p/2)/colLz; %Anger, "Scintillation Camera with Multichannel Collimators", J Nucl Med 5:515-531 (1964)
            Rg_xy = 2*options.colR*(Distances+options.cr_p/2)/colLxy;
            Rg_z(Rg_z<0) = 0;
            Rg_xy(Rg_xy<0) = 0;
            FWHMrot = 1;

            FWHM_z = sqrt(Rg_z.^2+options.iR^2);
            FWHM_xy = sqrt(Rg_xy.^2+options.iR^2);
            FWHM_z_pixel = FWHM_z/options.dz(1);
            FWHM_xy_pixel = FWHM_xy/options.dy(1);
            expr = FWHM_xy_pixel.^2-FWHMrot^2;
            expr(expr<=0) = 10^-16;
            FWHM_WithinPlane = sqrt(expr);

            %Parametrit CDR-mallinnukseen
            options.sigmaZ = FWHM_z_pixel./(2*sqrt(2*log(2)));
            options.sigmaXY = FWHM_WithinPlane./(2*sqrt(2*log(2)));
        end
        maxI = max([options.Nx(1), options.Ny(1), options.Nz(1)]);
        y = double(double(maxI) / 2 - 1:-1:-double(maxI) / 2 + 1);
        x = double(double(maxI) / 2 - 1:-1:-double(maxI) / 2 + 1);
        % xx = repmat(x', 1,options.Nz);
        xx = repmat(x', 1,size(x,2));
        yy = repmat(y, size(xx,1),1);
        if ~isfield(options, 'sigmaXY')
            s1 = double(repmat(permute(options.sigmaZ,[4 3 2 1]), size(xx,1), size(yy,2), 1));
            options.gFilter = exp(-(xx.^2 + yy.^2)./(2*s1.^2)); % .*  (1 / (2*pi*s1);
        else
            s1 = double(repmat(permute(options.sigmaZ,[4 3 2 1]), size(xx,1), size(yy,2), 1));
            s2 = double(repmat(permute(options.sigmaXY,[4 3 2 1]), size(xx,1), size(yy,2), 1));
            options.gFilter = exp(-(xx.^2./(2*s1.^2) + yy.^2./(2*s2.^2))); % .* (1 / (2*pi*s1.*s2)) ;
        end
        [rowE,colE] = find(options.gFilter(:,:,max(1,floor(end/4))) > 1e-6);
        [rowS,colS] = find(options.gFilter(:,:,max(1,floor(end/4))) > 1e-6);
        rowS = min(rowS);
        colS = min(colS);
        rowE = max(rowE);
        colE = max(colE);
        options.gFilter = options.gFilter(rowS:rowE, colS:colE, :);
        options.gFilter = options.gFilter ./ sum(options.gFilter, [1 2]); % Due to truncation the sum of every slice is not 1. Thus normalize by slice.
    end

    options.blurPlanes = (options.FOVa_x(1)/2 - (options.radiusPerProj .* cosd(panelTilt) - options.CORtoDetectorSurface)) / options.dx(1); % PSF shift
    options.blurPlanes2 = options.radiusPerProj .* sind(panelTilt) / options.dx(1); % Panel shift

    if options.implementation == 2
        options.blurPlanes = int32(options.blurPlanes);
        options.blurPlanes2 = int32(options.blurPlanes2);
    end
    if ~isfield(options,'angles') || numel(options.angles) == 0
        options.angles = (repelem(options.startAngle, options.nProjections / options.nHeads, 1) + repmat((0:options.angleIncrement:options.angleIncrement * (options.nProjections / options.nHeads - 1))', options.nHeads, 1) + options.offangle);
    end
    options.uu = 1;
    options.ub = 1;
    options.blurPlanesF = single(options.blurPlanesF);
    options.blurPlanes2F = single(options.blurPlanes2F);
    if options.useSingles
        options.gFilter = single(options.gFilter);
        options.angles = single(options.angles);
        options.swivelAngles = single(options.swivelAngles);
        options.radiusPerProj = single(options.radiusPerProj);
    end
end
end
