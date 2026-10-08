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
        % Keep the central detector-normal vector unchanged. Scale each source-detector shift difference so that the two angular response components use their corresponding septa lengths.
        referenceLength = options.colD + 0.5 * options.colL;
        lengthXY = options.colD + 0.5 * options.colLxy;
        lengthZ = options.colD + 0.5 * options.colLz;
        if referenceLength ~= lengthZ
            detectorShiftXY = options.rayShiftsDetector(1:2:end,:,:,:);
            options.rayShiftsSource(1:2:end,:,:,:) = detectorShiftXY + ...
                (options.rayShiftsSource(1:2:end,:,:,:) - detectorShiftXY) * (referenceLength / lengthZ);
        end
        if referenceLength ~= lengthXY
            detectorShiftZ = options.rayShiftsDetector(2:2:end,:,:,:);
            options.rayShiftsSource(2:2:end,:,:,:) = detectorShiftZ + ...
                (options.rayShiftsSource(2:2:end,:,:,:) - detectorShiftZ) * (referenceLength / lengthXY);
        end
    end
    if ismember(options.implementation, [2, 3, 5])
        options.rayShiftsDetector = single(options.rayShiftsDetector(:));
        options.rayShiftsSource = single(options.rayShiftsSource(:));
    end
end
if ismember(options.projector_type, [12, 2, 21, 22]) % Orthogonal distance ray tracer
    % Frey, E. C., & Tsui, B. M. W. (n.d.). Collimator-Detector Response Compensation in SPECT. Quantitative Analysis in Nuclear Medicine Imaging, 141–166. doi:10.1007/0-387-25444-7_5 
    if (~isfield(options,'coneOfResponseStdCoeffA') || options.coneOfResponseStdCoeffA < 0)
        options.coneOfResponseStdCoeffA = 2*options.colR/options.colL; % See equation (6) of book chapter
    end
    if (~isfield(options,'coneOfResponseStdCoeffB') || options.coneOfResponseStdCoeffB < 0)
        options.coneOfResponseStdCoeffB = 2*options.colR/options.colL*(options.colL+options.colD+options.cr_p/2);
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
    DistanceToFirstRow = 0.5*options.dx(1);
    Distances = repmat(DistanceToFirstRow,1,options.Nx(1)*4)+repmat((0:double(options.Nx(1)*4)-1)*double(options.dx(1)),length(DistanceToFirstRow),1);
    Distances = Distances-options.colL-options.colD; %these are distances to the actual detector surface

    if (~isfield(options,'gFilter'))
        if ~isfield(options, 'sigmaZ')
            % Axial and transaxial collimator responses use their own septa
            % lengths (colLz/colLxy), each falling back to the uniform-hole
            % length colL when not positive. Ported from the newer Python
            % behavior (detcoord.py, projector type 6 CDRF computation) so
            % that asymmetric collimator holes (colLxy ~= colLz) are handled
            % identically in MATLAB and Python.
            if options.colLxy > 0
                colLxy = options.colLxy;
            else
                colLxy = options.colL;
            end
            if options.colLz > 0
                colLz = options.colLz;
            else
                colLz = options.colL;
            end
            Rg_z = 2*options.colR*(options.colL+options.colD+(Distances)+options.cr_p/2)/colLz; %Anger, "Scintillation Camera with Multichannel Collimators", J Nucl Med 5:515-531 (1964)
            Rg_xy = 2*options.colR*(options.colL+options.colD+(Distances)+options.cr_p/2)/colLxy;
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

    panelTilt = options.swivelAngles - options.angles + 180;
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
    if options.useSingles
        options.gFilter = single(options.gFilter);
        options.angles = single(options.angles);
        options.swivelAngles = single(options.swivelAngles);
        options.radiusPerProj = single(options.radiusPerProj);
    end
end
end
