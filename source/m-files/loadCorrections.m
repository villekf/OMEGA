function [options] = loadCorrections(options, RandProp, ScatterProp)
%% SET UP CORRECTIONS
% This function sets up (loads) the necessary corrections that are applied.
% Included are attenuation correction, normalization correction and randoms
% and/or scatter correction. The variables used for correction are saved in
% the options-struct that also contains the input variables except for the
% number of rings in the current machine.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2022-2026 Ville-Veikko Wettenhovi
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

folder = fileparts(which('set_up_corrections.m'));
folder = [folder(1:end-(6 + 8)), 'mat-files/'];
folder = strrep(folder, '\','/');

x = [];
y = [];

if ~isfield(options,'reconstruct_trues')
    options.reconstruct_trues = false;
end
if ~isfield(options,'reconstruct_scatter')
    options.reconstruct_scatter = false;
end
if ~isfield(options,'TOF_bins') || options.TOF_bins == 0
    options.TOF_bins = 1;
end
if ~isfield(options,'use_machine')
    options.use_machine = 0;
end
if ~isfield(options,'Nt')
    options.Nt = 1;
end
% The scatter data is checked with TOF data below (checkTOFScatter), which sets options.randomsTOF if the scatter data is TOF-resolved
options.randomsTOF = false;

if options.listmode && options.normalization_correction && options.Nt > 1 && ~options.CT
    warning('Normalization correction is not (yet) supported for dynamic list-mode data! Disabling it.')
    options.normalization_correction = false;
end

% The sinogram dimensions of the normalization data. Index-based data has the dimensions collapsed (see projectorClass), so
% the values saved before the collapse are used for it
if options.useIndexBasedReconstruction && isfield(options, 'indexSinogramDims')
    normDims = options.indexSinogramDims;
else
    normDims = [options.Ndist, options.Nang, options.NSinos, options.TotSinos];
end

if ~options.use_raw_data && options.randoms_correction
    if options.partitions == 1
        load_string = [options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
            num2str(options.TotSinos) '_span' num2str(options.span)];
        if options.use_machine == 0
            sinoFile = [load_string '.mat'];
        elseif options.use_machine == 1
            sinoFile = [load_string '_listmode.mat'];
        elseif options.use_machine == 2
            sinoFile = [load_string '_machine_sinogram.mat'];
        elseif options.use_machine == 3
            sinoFile = [load_string '_listmode_sinogram.mat'];
        end
    else
        load_string = [options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_' ...
            num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
            num2str(options.span)];
        load_string2 = [options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
            num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
            num2str(options.span)];
        if options.use_machine == 0
            sinoFile = [load_string '.mat'];
            if exist(sinoFile, 'file') == 0
                sinoFile = [load_string2 '.mat'];
            end
        elseif options.use_machine == 1
            sinoFile = [load_string '_listmode.mat'];
            if exist(sinoFile, 'file') == 0
                sinoFile = [load_string2 '_listmode.mat'];
            end
        elseif options.use_machine == 2
            sinoFile = [load_string '_machine_sinogram.mat'];
            if exist(sinoFile, 'file') == 0
                sinoFile = [load_string2 '_machine_sinogram.mat'];
            end
        elseif options.use_machine == 3
            sinoFile = [load_string '_listmode_sinogram.mat'];
            if exist(sinoFile, 'file') == 0
                sinoFile = [load_string2 '_listmode_sinogram.mat'];
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% Attenuation correction %%%%%%%%%%%%%%%%%%%%%%%%%%
attenuationInmm = false; % True if the attenuation map was generated here (CT DICOM) and is already in 1/mm
if options.attenuation_correction && ~options.SPECT % PET attenuation
    if ~isfield(options,'vaimennus') || isempty(options.vaimennus)
        if ~isempty(options.attenuation_datafile) && strcmp(options.attenuation_datafile(end-2:end), 'mhd')
            [options.vaimennus, apuStruct] = loadMetaImage(options.attenuation_datafile);
            if options.CT_attenuation
                if round(apuStruct.EleSpacing(1)*100)/100 > round(options.FOVa_x(1) / double(options.Nx(1))*100)/100 ||  round(apuStruct.EleSpacing(1)*100)/100 < round(options.FOVa_x(1) / double(options.Nx(1))*100)/100
                    options.vaimennus = options.vaimennus .* (apuStruct.EleSpacing(1) / (options.FOVa_x(1) / double(options.Nx(1))));
                end
            end
        elseif ~isempty(options.attenuation_datafile)
            data = load(options.attenuation_datafile);
            variables = fieldnames(data);
            options.vaimennus = double(data.(variables{1}));
            clear data
        else
            [options.file, options.fpath] = uigetfile('*.*','Select attenuation file');
            if isequal(options.file, 0)
                error('No file was selected')
            end
            nimi = [options.fpath options.file];
            if strcmp(nimi(end-2:end), 'mhd')
                [options.vaimennus, ~] = loadMetaImage(nimi);
            elseif strcmp(nimi(end-2:end), 'h33')
                [options.vaimennus] = loadInterfile(nimi);
            else
                data = load(nimi);
                variables = fieldnames(data);
                options.vaimennus = double(data.(variables{1}));
                clear data
            end
        end
    end
    if options.CT_attenuation && options.useMultiResolutionVolumes && options.nMultiVolumes > 0
        options = prepareMultiResolutionAttenuation(options);
    elseif options.CT_attenuation
        if size(options.vaimennus,1) ~= options.Nx(1) || size(options.vaimennus,2) ~= options.Ny(1) || size(options.vaimennus,3) ~= options.Nz(1)
            if size(options.vaimennus,1) ~= options.N(1)
                warning('Error: Attenuation data is of different size than the reconstructed image. Attempting resize.')
                try
                    options.vaimennus = imresize3(options.vaimennus, [options.Nx(1), options.Ny(1), options.Nz(1)]);
                catch ME
                    error('Resize failed!')
                end
            end
        end
        options.vaimennus = options.vaimennus(:);
        if isfield(options,'rotateAttImage') && options.rotateAttImage ~= 0
            atn = reshape(options.vaimennus, options.Nx(1), options.Ny(1), options.Nz(1));
            atn = rot90(atn,options.rotateAttImage);
            options.vaimennus = atn(:);
        end
        if isfield(options,'flipAttImageXY') && options.flipAttImageXY
            atn = reshape(options.vaimennus, options.Nx(1), options.Ny(1), options.Nz(1));
            atn = fliplr(atn);
            options.vaimennus = atn(:);
        end
        if isfield(options,'flipAttImageZ') && options.flipAttImageZ
            atn = reshape(options.vaimennus, options.Nx(1), options.Ny(1), options.Nz(1));
            atn = flip(atn,3);
            options.vaimennus = atn(:);
        end
    end
elseif options.attenuation_correction && options.SPECT % SPECT attenuation
    if ~isfield(options, 'vaimennus') || (isempty(options.vaimennus) && ~iscell(options.vaimennus))
        if isempty(options.fpathCT) % Get folder if necessary
            options.fpathCT = uigetdir(matlabroot,'Select the folder containing attenuation DICOM files');
        end

        % Read volume and remove dimensions with size 1
        [CTvol, spatial, ~] = dicomreadVolume(options.fpathCT);
        CTvol = double(CTvol);
        CTvol = squeeze(CTvol);

        % DICOM PixelSpacing is [row spacing, column spacing] and ImageSize is [rows, columns, slices]. Columns run along x
        % and rows along y (as in imref3d)
        dx = spatial.PixelSpacings(1, 2); % Column spacing
        dy = spatial.PixelSpacings(1, 1); % Row spacing
        if size(spatial.PatientPositions, 1) < 2
            error('The CT DICOM series needs at least two slices!');
        end
        dz = abs(spatial.PatientPositions(2, 3) - spatial.PatientPositions(1, 3));
        nColsCT = spatial.ImageSize(2);
        nRowsCT = spatial.ImageSize(1);

        xLimits = spatial.PatientPositions(1,1) + [0, dx * nColsCT] - 0.5 * dx;
        yLimits = spatial.PatientPositions(1,2) + [0, dy * nRowsCT] - 0.5 * dy;
        zLimits = [min(spatial.PatientPositions(:,3)) max(spatial.PatientPositions(:,3))] + dz * [-0.5, 0.5];

        refCT = imref3d(size(CTvol), xLimits, yLimits, zLimits);

        % Apply DICOM Rescale slope and intercept (of each slice, in the slice order of dicomreadVolume)
        [HU_a, HU_b] = HU_rescale_coeff(options.fpathCT, spatial.PatientPositions);
        CTvol = HU_a .* CTvol + HU_b;
        CTvol(CTvol < -1000) = -1000;
        % Convert to linear attenuation coefficients
        MUvol = HU_to_mu(CTvol, options.keV);
        muAir = HU_to_mu(-1000, options.keV);

        % Use the pre-split full-FOV reference for multi-resolution maps.
        % Shared attenuation preparation aligns maps to the active volume grids.
        if options.useMultiResolutionVolumes && options.nMultiVolumes > 0
            fineDims = double([options.NxFull, options.NyFull, options.NzFull]);
            if isfield(options, 'imageAttenuationFineDims')
                fineDims = double(options.imageAttenuationFineDims(:).');
            end
            fineFOV = double([options.FOVa_x(1), options.FOVa_y(1), options.axial_fov(1)]);
            if isfield(options, 'imageAttenuationFineFOV')
                fineFOV = double(options.imageAttenuationFineFOV(:).');
            end
            fineCenter = double([options.oOffsetX, options.oOffsetY, options.oOffsetZ]);
            if options.useEFOV
                fineCenter = fineCenter + double(options.eFOVShift(:).');
            end
            attenuationRef = imref3d(fineDims, fineCenter(1) + fineFOV(1) .* [-0.5, 0.5], ...
                fineCenter(2) + fineFOV(2) .* [-0.5, 0.5], ...
                fineCenter(3) + fineFOV(3) .* [-0.5, 0.5]);
        else
            attenuationRef = options.refSPECT;
        end

        % Use imwarp to change the CT volume limits
        tform = affinetform3d(eye(4)); % No scaling or rotation
        [MUvol, ~] = imwarp(MUvol, refCT, tform, OutputView=attenuationRef, FillValue=muAir, InterpolationMethod='linear');
        options.vaimennus = 1*MUvol;
        % The CT-derived attenuation map (HU_to_mu) is already in 1/mm
        attenuationInmm = true;
    end
    if options.CT_attenuation && options.useMultiResolutionVolumes && options.nMultiVolumes > 0
        options = prepareMultiResolutionAttenuation(options);
    end
    if options.Nt > 1 && iscell(options.vaimennus) && numel(options.vaimennus) == 1 && ...
            ~options.CT_attenuation && iscell(options.SinM) && ~isempty(options.SinM)
        projectionCounts = cellfun(@(x) size(x, 3), options.SinM);
        if all(projectionCounts == projectionCounts(1)) && ...
                numel(options.vaimennus{1}) == options.nRowsD * options.nColsD * projectionCounts(1)
            % A single cell can also represent one shared measurement-based
            % attenuation frame. Expand it so downstream per-frame reordering
            % can use each timeframe's projection order.
            options.vaimennus = repmat(options.vaimennus(:), options.Nt, 1);
        end
    end
    if options.Nt > 1
        if ~iscell(options.vaimennus)
            isPreparedImageAttenuation = options.CT_attenuation && ...
                ((isfield(options, 'imageAttenuationIsMultiResolution') && options.imageAttenuationIsMultiResolution) || ...
                (isfield(options, 'projector_type') && options.projector_type == 6 && ...
                isfield(options, 'multiResAttenuationPadded') && options.multiResAttenuationPadded));
            isDynamicMeasurementAttenuation = false;
            isStaticMeasurementAttenuation = false;
            if ~options.CT_attenuation && iscell(options.SinM) && ~isempty(options.SinM)
                projectionCounts = cellfun(@(x) size(x, 3), options.SinM);
                detectorImageSize = options.nRowsD * options.nColsD;
                isDynamicMeasurementAttenuation = numel(options.vaimennus) == ...
                    detectorImageSize * sum(projectionCounts);
                isStaticMeasurementAttenuation = all(projectionCounts == projectionCounts(1)) && ...
                    numel(options.vaimennus) == detectorImageSize * projectionCounts(1);
            end
            if ~isDynamicMeasurementAttenuation && ~isStaticMeasurementAttenuation && ...
                    ~isPreparedImageAttenuation
                error("With dynamic reconstruction the attenuation map needs to be a cell type");
            end
        elseif numel(options.vaimennus) ~= options.Nt
            error("No attenuation map for each timestep")
        end
    end
else
    options.vaimennus = 0;
end
% attIncm (1/cm to 1/mm) applies only to attenuation images. A measurement-domain (sinogram) attenuation correction is
% unitless and is never scaled. SPECT attenuation is always an image. Data of the size of the reconstructed image is an image
% as well even if CT_attenuation is false (same test as in parseInputData, which switches CT_attenuation to true for it).
% The attenuation map generated here from the CT DICOM series is already in 1/mm and is never scaled
if options.attenuation_correction && options.attIncm && attenuationInmm
    disp('The attenuation map generated from the CT DICOM series is already in 1/mm. attIncm is ignored for it.')
end
if iscell(options.vaimennus) && ~isempty(options.vaimennus)
    attenuationFirst = options.vaimennus{1};
else
    attenuationFirst = options.vaimennus;
end
attenuationIsImage = options.attenuation_correction && (options.SPECT || options.CT_attenuation || (numel(attenuationFirst) == options.Nx(1) * options.Ny(1) * options.Nz(1) ...
    && ~(options.nRowsD == options.Nx(1) && options.nColsD == options.Ny(1))));
if options.attenuation_correction && options.attIncm && attenuationIsImage && ~attenuationInmm && ...
        ~(options.CT_attenuation && options.useMultiResolutionVolumes && options.nMultiVolumes > 0)
    if iscell(options.vaimennus)
        for kk = 1 : numel(options.vaimennus)
            options.vaimennus{kk} = options.vaimennus{kk} ./ 10;
        end
    else
        options.vaimennus = options.vaimennus ./ 10;
    end
end


if ~iscell(options.vaimennus) 
    if (options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles)
        options.vaimennus = single(options.vaimennus);
    else
        options.vaimennus = double(options.vaimennus);
    end
else
    for kk = 1:numel(options.vaimennus)
        if (options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles)
            options.vaimennus{kk} = single(options.vaimennus{kk}(:));
        else
            options.vaimennus{kk} = double(options.vaimennus{kk}(:));
        end
    end
end



if ~options.SPECT
    % Load (or compute) normalization correction coefficients. This is done before the scatter is normalized (normalize_scatter),
    % which needs the coefficients. A normalization that is already input (options.normalization) is used as it is.
    % normalizationInverted is true if and only if the normalization was loaded from a file and inverted in this call
    normalizationInverted = false;
    % CT normalization is always the user input array options.normalization (no file lookup or dialog), a multiplicative factor of
    % each measurement. ctNormMode selects how it is applied, see ctNormalizationSetup (0 = no CT normalization)
    ctNormMode = 0;
    if options.CT && options.normalization_correction
        [options, ctNormMode] = ctNormalizationSetup(options);
    end
    if (options.normalization_correction && options.corrections_during_reconstruction) && ~options.use_user_normalization && ~options.CT
        if ~isfield(options,'normalization') || isempty(options.normalization)
            if options.use_raw_data
                norm_file = [folder options.machine_name '_normalization_listmode.mat'];
                if exist(norm_file, 'file') == 2
                    options.normalization = loadStructFromFile(norm_file,'normalization');
                else
                    error('Normalization correction selected, but no normalization data found')
                end
            else
                norm_file = [folder options.machine_name '_normalization_' num2str(normDims(1)) 'x' num2str(normDims(2)) '_span' num2str(options.span) '.mat'];
                if exist(norm_file, 'file') == 2
                    options.normalization = loadStructFromFile(norm_file,'normalization');
                else
                    error('Normalization correction selected, but no normalization data found')
                end
            end
            % Coefficients are inverted to efficiencies; zero coefficients (bins without efficiency) stay zero
            options.normalization(options.normalization < 0) = 1;
            normCoef = options.normalization(:);
            options.normalization = zeros(size(normCoef), 'like', normCoef);
            normNz = normCoef ~= 0;
            options.normalization(normNz) = 1 ./ normCoef(normNz);
            clear normCoef normNz
            if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                options.normalization = single(options.normalization);
            end
            normalizationInverted = true;
        end
        if options.useIndexBasedReconstruction
            % One value for each event instead of a sinogram, see indexBasedNormalization
            options = indexBasedNormalization(options, true);
        end
    elseif options.normalization_correction && options.use_user_normalization && options.corrections_during_reconstruction && ~options.CT
        if ~isfield(options,'normalization') || isempty(options.normalization)
            [file, fpath] = selectNormalizationFile(options);
            [~, ~, normExt] = fileparts(file);
            if any(strfind(file, '.nrm'))
                if options.use_raw_data
                    error('Inveon normalization data cannot be used with raw list-mode data')
                end
                fid = fopen(fullfile(fpath, file));
                options.normalization = fread(fid, inf, 'single=>single',0,'l');
                fclose(fid);
                % The size of index-based data is checked when the normalization is converted
                if numel(options.normalization) ~= normDims(1) * normDims(2) * normDims(4) && ~options.use_raw_data && ~options.useIndexBasedReconstruction
                    error('Size mismatch between the current data and the normalization datafile')
                end
                options.InveonNorm = true;
            elseif any(strcmpi(normExt, {'.n', '.ptd'}))
                % Biograph (mCT/Vision) normalization component file, the coefficients multiply the data
                if options.use_raw_data
                    error('Biograph normalization data cannot be used with raw list-mode data')
                end
                options.normalization = loadNormalizationBiograph(fullfile(fpath, file), normalizationOptions(options, normDims));
                if numel(options.normalization) ~= normDims(1) * normDims(2) * normDims(4)
                    error('Size mismatch between the current data and the normalization datafile')
                end
            else
                data = load(fullfile(fpath, file));
                variables = fieldnames(data);
                if any(strcmp(variables,'normalization'))
                    options.normalization = data.(variables{strcmp(variables,'normalization')});
                else
                    options.normalization = data.(variables{1});
                end
                clear data
                % The size of index-based data is checked when the normalization is converted
                if numel(options.normalization) ~= normDims(1) * normDims(2) * normDims(4) && ~options.use_raw_data && ~options.useIndexBasedReconstruction
                    error('Size mismatch between the current data and the normalization datafile')
                elseif numel(options.normalization) ~= sum(1:options.detectors) && options.use_raw_data
                    error('Size mismatch between the current data and the normalization datafile')
                end
            end
            if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                options.normalization = single(options.normalization);
            else
                options.normalization = double(options.normalization);
            end
            % Coefficients are inverted to efficiencies; zero coefficients (bins without efficiency) stay zero
            normCoef = options.normalization(:);
            options.normalization = zeros(size(normCoef), 'like', normCoef);
            normNz = normCoef ~= 0;
            options.normalization(normNz) = 1 ./ normCoef(normNz);
            clear normCoef normNz
            normalizationInverted = true;
        end
        if options.useIndexBasedReconstruction
            % One value for each event instead of a sinogram, see indexBasedNormalization
            options = indexBasedNormalization(options, true);
        end
    end

    % Pair normalization table of the sensitivity image of list-mode data with custom coordinates (user input, see OMEGA_error_check and
    % indexBasedNormalization for index-based data): single precision vector, in the order of the sensitivity image kernel. The
    % normalization has to be given for the events
    if options.listmode && ~options.CT && options.normalization_correction && options.compute_sensitivity_image && ~options.useIndexBasedReconstruction ...
            && isfield(options, 'normSens') && ~isempty(options.normSens) && ~(isfield(options, 'sensXY') && ~isempty(options.sensXY))
        if iscell(options.SinM)
            nEventsPair = numel(options.SinM{1});
        else
            nEventsPair = numel(options.SinM);
        end
        if nEventsPair == 0
            nEventsPair = numel(options.x) / 6;
        end
        if numel(options.normalization) ~= nEventsPair
            error(['With a pair normalization table (options.normSens) options.normalization has to contain one value for each event (' ...
                num2str(nEventsPair) ' elements), but it has ' num2str(numel(options.normalization)) ' elements! The table is used only in ' ...
                'the sensitivity image, a sinogram normalization is not converted to the events.'])
        end
        options.normSens = single(options.normSens(:));
    end

    % List-mode and index-based PET data: options.ScatterC is the additive term of each event (validated in projectorClass, see
    % checkListModeScatter), used as it is: no division by the TOF bins (see scatterTOFApplies) and no sinogram operations. The
    % randoms are not used with list-mode data (the randoms are added into options.ScatterC by the user or input as minus one measurements)
    if options.listmode && ~options.CT && ~options.SPECT && options.scatter_correction && options.subtract_scatter ...
            && ~options.reconstruct_trues && ~options.reconstruct_scatter
        % A SinDelayed that is left from an earlier call with the returned options (randoms_correction is then set) is formed again below
        if isfield(options,'SinDelayed') && numel(options.SinDelayed) > 1
            if ~options.randoms_correction
                warning('Randoms (options.SinDelayed) are not used with list-mode and index-based data! Add the randoms into options.ScatterC. Ignoring options.SinDelayed.')
            end
            options.SinDelayed = [];
        end
        options.randoms_correction = false;
        if iscell(options.ScatterC) || numel(options.ScatterC) <= 1
            error('Scatter correction of list-mode and index-based data needs options.ScatterC as a single array with one value for each event!')
        end
        options.ScatterC = options.ScatterC(:);
    end

    % Apply randoms and scatter (if available) correction
    % Apply randoms smoothing/variance reduction
    if (options.randoms_correction || options.scatter_correction) && options.ordinaryPoisson && ~options.reconstruct_trues...
            && ~options.reconstruct_scatter
        randoms_correction = true;
        r_exist = isfield(options,'SinDelayed') && numel(options.SinDelayed) > 1;
        s_exist = isfield(options,'ScatterC') && numel(options.ScatterC) > 1;
        if r_exist && isfield(options,'SinDelayed') && ~iscell(options.SinDelayed) && isscalar(options.SinDelayed) && randoms_correction
            r_exist = false;
        elseif r_exist && isfield(options,'SinDelayed') && iscell(options.SinDelayed) && isscalar(options.SinDelayed{1}) && randoms_correction
            r_exist = false;
        end
        if exist('RandProp','var') == 0 || isempty(RandProp)
            RandProp.smoothing = false;
            RandProp.variance_reduction = false;
        end
        if exist('ScatterProp','var') == 0 || isempty(ScatterProp)
            ScatterProp.smoothing = false;
            ScatterProp.variance_reduction = false;
            ScatterProp.normalization = false;
        end
        % The scatter has been normalized by an earlier call with the same options (the normalization is then a prefilled one)
        scatterAlreadyNormalized = isfield(options,'scatterNormalized') && options.scatterNormalized;
        if options.scatter_correction && options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization
            if ~options.corrections_during_reconstruction
                % The normalization is loaded in the precorrection after the scatter has been normalized
                if ~isfield(options,'normalization') || isempty(options.normalization)
                    error('normalize_scatter with a normalization that is loaded from a file needs corrections during reconstruction (options.corrections_during_reconstruction = true)!')
                elseif options.useIndexBasedReconstruction
                    error('Index-based reconstruction with normalize_scatter needs corrections during reconstruction (options.corrections_during_reconstruction = true)!')
                end
            elseif options.useIndexBasedReconstruction && isfield(options,'ScatterC')
                % The normalization is one value for each event (indexBasedNormalization above)
                if iscell(options.ScatterC)
                    nScatter = numel(options.ScatterC{1});
                else
                    nScatter = numel(options.ScatterC);
                end
                if nScatter > 1 && nScatter ~= numel(options.normalization)
                    error(['Index-based reconstruction with normalize_scatter needs the scatter (options.ScatterC) as one value for each event (' ...
                        num2str(numel(options.normalization)) ' elements), but ' num2str(nScatter) ' elements were input!'])
                end
            end
        end
        if r_exist && randoms_correction
            if (~options.variance_reduction && RandProp.variance_reduction) || (~options.randoms_smoothing && RandProp.smoothing)
                options.SinDelayed = loadStructFromFile(sinoFile, 'raw_SinDelayed');
                RandProp.variance_reduction = false;
                RandProp.smoothing = false;
            end
        end
        if options.randoms_correction && ~r_exist
            options = loadDelayedData(options);
            r_exist = true;
        end
        if options.scatter_correction && ((~isfield(options,'ScatterC') || numel(options.ScatterC) == 0) || (ScatterProp.normalization && ~options.normalize_scatter) || ...
                (ScatterProp.smoothing && ~options.scatter_smoothing) || (ScatterProp.variance_reduction && ~options.scatter_variance_reduction))
            if ScatterProp.normalization && ~options.normalize_scatter
                warning('Scatter data is normalized, but no normalized scatter selected. Input non-normalized scatter data or put options.normalize_scatter = true!')
            end
            if ScatterProp.smoothing && ~options.scatter_smoothing
                warning('Scatter data is smoothed, but no scatter smoothing selected. Input non-smoothed scatter data or put options.scatter_smoothing = true!')
            end
            if ScatterProp.variance_reduction && ~options.scatter_variance_reduction
                warning('Scatter data is variance corrected, but no scatter variance correction selected. Input non-variance corrected scatter data or put options.scatter_variance_reduction = true!')
            end
            options = loadScatterData(options);
            s_exist = true;
        end
        % Check the size of the scatter data with TOF data (also prints the message that the scatter data is divided among the TOF bins)
        options = checkTOFScatter(options);
        if r_exist || s_exist
            if r_exist && iscell(options.SinDelayed)
                for kk = 1 : options.partitions
                    if r_exist
                        if options.variance_reduction && ~RandProp.variance_reduction
                            options.SinDelayed{kk} = Randoms_variance_reduction(single(options.SinDelayed{kk}), options);
                        end
                        if options.randoms_smoothing && ~RandProp.smoothing
                            options.SinDelayed{kk} = randoms_smoothing(single(options.SinDelayed{kk}), options);
                        end
                        options.SinDelayed{kk} = options.SinDelayed{kk}(:);
                        if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                            options.SinDelayed{kk} = options.SinDelayed{kk}(1:options.NSinos*options.Ndist*options.Nang);
                        end
                    end
                    if options.scatter_correction
                        % A numeric ScatterC is a single scatter sinogram shared by all the frames, it is processed only once
                        if sum(size(options.ScatterC)) > 1 && ~iscell(options.ScatterC) && kk == 1
                            if size(options.ScatterC,1) ~= options.Ndist && ~options.use_raw_data && ~options.listmode
                                options.ScatterC = permute(options.ScatterC,[2 1 3 4]);
                            end
                            if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization && ~scatterAlreadyNormalized
                                options.ScatterC = normalizeScatter(options.ScatterC, options, normalizationInverted);
                            end
                            if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                                options.ScatterC = Randoms_variance_reduction(single(options.ScatterC), options);
                            end
                            if options.scatter_smoothing && ~ScatterProp.smoothing
                                options.ScatterC = randoms_smoothing(single(options.ScatterC), options);
                            end
                            if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                                options.ScatterC = cropScatterSinograms(options.ScatterC, options);
                            end
                        elseif iscell(options.ScatterC)
                            if sum(size(options.ScatterC{1})) > 1
                                if size(options.ScatterC{kk},1) ~= options.Ndist && ~options.use_raw_data && ~options.listmode
                                    options.ScatterC{kk} = permute(options.ScatterC{kk},[2 1 3 4]);
                                end
                            end
                            if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization && ~scatterAlreadyNormalized
                                options.ScatterC{kk} = normalizeScatter(options.ScatterC{kk}, options, normalizationInverted);
                            end
                            if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                                options.ScatterC{kk} = Randoms_variance_reduction(single(options.ScatterC{kk}), options);
                            end
                            if options.scatter_smoothing && ~ScatterProp.smoothing
                                options.ScatterC{kk} = randoms_smoothing(single(options.ScatterC{kk}), options);
                            end
                            if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                                options.ScatterC{kk} = cropScatterSinograms(options.ScatterC{kk}, options);
                            end
                        end
                        if options.subtract_scatter
                            if iscell(options.ScatterC)
                                scatterFrame = options.ScatterC{kk};
                            else
                                scatterFrame = options.ScatterC;
                            end
                            if r_exist
                                if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.implementation == 4
                                    options.SinDelayed{kk} = addScatterPerTOFBin(single(options.SinDelayed{kk}) / options.TOF_bins, single(scatterPerTOFBin(scatterFrame(:), options)));
                                else
                                    options.SinDelayed{kk} = addScatterPerTOFBin(double(options.SinDelayed{kk}) / options.TOF_bins, double(scatterPerTOFBin(scatterFrame(:), options)));
                                end
                            else
                                if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.implementation == 4
                                    options.SinDelayed{kk} = single(scatterPerTOFBin(scatterFrame(:), options));
                                else
                                    options.SinDelayed{kk} = double(scatterPerTOFBin(scatterFrame(:), options));
                                end
                            end
                        end
                    end
                    if (~options.scatter_correction || ~options.subtract_scatter) && options.TOF_bins > 1 && r_exist
                        if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.implementation == 4
                            options.SinDelayed{kk} = single(options.SinDelayed{kk}) / options.TOF_bins;
                        else
                            options.SinDelayed{kk} = double(options.SinDelayed{kk}) / options.TOF_bins;
                        end
                    end
                end
            else
                if r_exist
                    if options.variance_reduction && ~RandProp.variance_reduction
                        options.SinDelayed = Randoms_variance_reduction(single(options.SinDelayed), options);
                    end
                    if options.randoms_smoothing && ~RandProp.smoothing
                        options.SinDelayed = randoms_smoothing(single(options.SinDelayed), options);
                    end
                    options.SinDelayed = options.SinDelayed(:);
                    if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                        options.SinDelayed = options.SinDelayed(1:options.NSinos*options.Ndist*options.Nang);
                    end
                end
                if options.scatter_correction && isfield(options,'ScatterC')
                    if isempty(ScatterProp)
                        ScatterProp.smoothing = false;
                        ScatterProp.variance_reduction = false;
                        ScatterProp.normalization = false;
                    end
                    if sum(size(options.ScatterC)) > 1 && ~iscell(options.ScatterC)
                        if size(options.ScatterC,1) ~= options.Ndist && ~options.use_raw_data && ~options.listmode
                            options.ScatterC = permute(options.ScatterC,[2 1 3 4]);
                        end
                        if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization && ~scatterAlreadyNormalized
                            options.ScatterC = normalizeScatter(options.ScatterC, options, normalizationInverted);
                        end
                        if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                            options.ScatterC = Randoms_variance_reduction(single(options.ScatterC), options);
                        end
                        if options.scatter_smoothing && ~ScatterProp.smoothing
                            options.ScatterC = randoms_smoothing(single(options.ScatterC), options);
                        end
                        if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                            options.ScatterC = cropScatterSinograms(options.ScatterC, options);
                        end
                        if options.subtract_scatter
                            if r_exist
                                if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                                    options.SinDelayed = addScatterPerTOFBin(single(options.SinDelayed) / options.TOF_bins, single(scatterPerTOFBin(options.ScatterC(:), options)));
                                else
                                    options.SinDelayed = addScatterPerTOFBin(double(options.SinDelayed) / options.TOF_bins, double(scatterPerTOFBin(options.ScatterC(:), options)));
                                end
                            else
                                if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                                    options.SinDelayed = single(scatterPerTOFBin(options.ScatterC(:), options));
                                else
                                    options.SinDelayed = double(scatterPerTOFBin(options.ScatterC(:), options));
                                end
                            end
                        end
                    elseif iscell(options.ScatterC)
                        if options.subtract_scatter
                            if r_exist
                                % The numeric randoms are a single sinogram shared by all the frames of the scatter data
                                randomsData = options.SinDelayed;
                            end
                            options.SinDelayed = cell(numel(options.ScatterC),1);
                        end
                        for kk = 1 : options.partitions
                            if sum(size(options.ScatterC{kk})) > 1
                                if size(options.ScatterC{kk},1) ~= options.Ndist && ~options.use_raw_data && ~options.listmode
                                    options.ScatterC{kk} = permute(options.ScatterC{kk},[2 1 3 4]);
                                end
                            end
                            if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization && ~scatterAlreadyNormalized
                                options.ScatterC{kk} = normalizeScatter(options.ScatterC{kk}, options, normalizationInverted);
                            end
                            if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                                options.ScatterC{kk} = Randoms_variance_reduction(single(options.ScatterC{kk}), options);
                            end
                            if options.scatter_smoothing && ~ScatterProp.smoothing
                                options.ScatterC{kk} = randoms_smoothing(single(options.ScatterC{kk}), options);
                            end
                            if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                                options.ScatterC{kk} = cropScatterSinograms(options.ScatterC{kk}, options);
                            end
                            if options.subtract_scatter
                                if r_exist
                                    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                                        options.SinDelayed{kk} = addScatterPerTOFBin(single(randomsData) / options.TOF_bins, single(scatterPerTOFBin(options.ScatterC{kk}(:), options)));
                                    else
                                        options.SinDelayed{kk} = addScatterPerTOFBin(double(randomsData) / options.TOF_bins, double(scatterPerTOFBin(options.ScatterC{kk}(:), options)));
                                    end
                                else
                                    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                                        options.SinDelayed{kk} = single(scatterPerTOFBin(options.ScatterC{kk}(:), options));
                                    else
                                        options.SinDelayed{kk} = double(scatterPerTOFBin(options.ScatterC{kk}(:), options));
                                    end
                                end
                            end
                        end
                        if options.subtract_scatter && r_exist && options.partitions == 1
                            options.SinDelayed = options.SinDelayed{1};
                        end
                    end
                end
                if (~options.scatter_correction || ~options.subtract_scatter) && options.TOF_bins > 1 && r_exist
                    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                        options.SinDelayed = single(options.SinDelayed) / options.TOF_bins;
                    else
                        options.SinDelayed = double(options.SinDelayed) / options.TOF_bins;
                    end
                end
            end
        else
            if options.randoms_correction && iscell(options.SinDelayed)
                for kk = 1 : options.partitions
                    if options.variance_reduction && ~RandProp.variance_reduction
                        options.SinDelayed{kk} = Randoms_variance_reduction(single(options.SinDelayed{kk}), options);
                    end
                    if options.randoms_smoothing && ~RandProp.smoothing
                        options.SinDelayed{kk} = randoms_smoothing(single(options.SinDelayed{kk}), options);
                    end
                    options.SinDelayed{kk} = options.SinDelayed{kk}(:);
                    if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                        options.SinDelayed{kk} = options.SinDelayed{kk}(1:options.NSinos*options.Ndist*options.Nang);
                    end
                    if options.scatter_correction
                        if sum(size(options.ScatterC)) > 1 && ~iscell(options.ScatterC)
                            if size(options.ScatterC,1) ~= options.Ndist && ~options.use_raw_data && ~options.listmode
                                options.ScatterC = permute(options.ScatterC,[2 1 3 4]);
                            end
                            if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization && ~scatterAlreadyNormalized
                                options.ScatterC = normalizeScatter(options.ScatterC, options, normalizationInverted);
                            end
                            if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                                options.ScatterC = Randoms_variance_reduction(single(options.ScatterC), options);
                            end
                            if options.scatter_smoothing && ~ScatterProp.smoothing
                                options.ScatterC = randoms_smoothing(single(options.ScatterC), options);
                            end
                        elseif iscell(options.ScatterC)
                            if sum(size(options.ScatterC{1})) > 1
                                if size(options.ScatterC{kk},1) ~= options.Ndist && ~options.use_raw_data && ~options.listmode
                                    options.ScatterC{kk} = permute(options.ScatterC{kk},[2 1 3 4]);
                                end
                            end
                            if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization && ~scatterAlreadyNormalized
                                options.ScatterC{kk} = normalizeScatter(options.ScatterC{kk}, options, normalizationInverted);
                            end
                            if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                                options.ScatterC{kk} = Randoms_variance_reduction(single(options.ScatterC{kk}), options);
                            end
                            if options.scatter_smoothing && ~ScatterProp.smoothing
                                options.ScatterC{kk} = randoms_smoothing(single(options.ScatterC{kk}), options);
                            end
                        end
                        if options.subtract_scatter
                            if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                                options.SinDelayed{kk} = addScatterPerTOFBin(single(options.SinDelayed{kk}) / options.TOF_bins, single(scatterPerTOFBin(options.ScatterC{kk}(:), options)));
                            else
                                options.SinDelayed{kk} = addScatterPerTOFBin(double(options.SinDelayed{kk}) / options.TOF_bins, double(scatterPerTOFBin(options.ScatterC{kk}(:), options)));
                            end
                        end
                    elseif ~options.scatter_correction && options.TOF_bins > 1
                        if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                            options.SinDelayed{kk} = single(options.SinDelayed{kk}) / single(options.TOF_bins);
                        else
                            options.SinDelayed{kk} = double(options.SinDelayed{kk}) / options.TOF_bins;
                        end
                    end
                end
            elseif options.randoms_correction
                if options.variance_reduction && ~RandProp.variance_reduction
                    options.SinDelayed = Randoms_variance_reduction(single(options.SinDelayed), options);
                end
                if options.randoms_smoothing && ~RandProp.smoothing
                    options.SinDelayed = randoms_smoothing(single(options.SinDelayed), options);
                end
                options.SinDelayed = options.SinDelayed(:);
                if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                    options.SinDelayed = options.SinDelayed(1:options.NSinos*options.Ndist*options.Nang);
                end
                if options.scatter_correction && isfield(options,'ScatterC')
                    if sum(size(options.ScatterC)) > 1 && ~iscell(options.ScatterC)
                        if size(options.ScatterC,1) ~= options.Ndist && ~options.use_raw_data && ~options.listmode
                            options.ScatterC = permute(options.ScatterC,[2 1 3 4]);
                        end
                        if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization && ~scatterAlreadyNormalized
                            options.ScatterC = normalizeScatter(options.ScatterC, options, normalizationInverted);
                        end
                        if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                            options.ScatterC = Randoms_variance_reduction(single(options.ScatterC), options);
                        end
                        if options.scatter_smoothing && ~ScatterProp.smoothing
                            options.ScatterC = randoms_smoothing(single(options.ScatterC), options);
                        end
                        if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                            options.ScatterC = cropScatterSinograms(options.ScatterC, options);
                        end
                        if options.subtract_scatter
                            if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                                options.SinDelayed = addScatterPerTOFBin(single(options.SinDelayed) / options.TOF_bins, single(scatterPerTOFBin(options.ScatterC(:), options)));
                            else
                                options.SinDelayed = addScatterPerTOFBin(double(options.SinDelayed) / options.TOF_bins, double(scatterPerTOFBin(options.ScatterC(:), options)));
                            end
                        end
                    elseif iscell(options.ScatterC)
                        for kk = 1 : options.partitions
                            if sum(size(options.ScatterC{kk})) > 1
                                if size(options.ScatterC{kk},1) ~= options.Ndist && ~options.use_raw_data && ~options.listmode
                                    options.ScatterC{kk} = permute(options.ScatterC{kk},[2 1 3 4]);
                                end
                            end
                            if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization && ~scatterAlreadyNormalized
                                options.ScatterC{kk} = normalizeScatter(options.ScatterC{kk}, options, normalizationInverted);
                            end
                            if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                                options.ScatterC{kk} = Randoms_variance_reduction(single(options.ScatterC{kk}), options);
                            end
                            if options.scatter_smoothing && ~ScatterProp.smoothing
                                options.ScatterC{kk} = randoms_smoothing(single(options.ScatterC{kk}), options);
                            end
                            if options.use_raw_data == false && options.NSinos ~= options.TotSinos
                                options.ScatterC{kk} = cropScatterSinograms(options.ScatterC{kk}, options);
                            end
                            if options.subtract_scatter
                                if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                                    options.SinDelayed = addScatterPerTOFBin(single(options.SinDelayed) / options.TOF_bins, single(scatterPerTOFBin(options.ScatterC{kk}(:), options)));
                                else
                                    options.SinDelayed = addScatterPerTOFBin(double(options.SinDelayed) / options.TOF_bins, double(scatterPerTOFBin(options.ScatterC{kk}(:), options)));
                                end
                            end
                        end
                    end
                elseif ~options.scatter_correction && options.TOF_bins > 1
                    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                        options.SinDelayed = single(options.SinDelayed) / single(options.TOF_bins);
                    else
                        options.SinDelayed = double(options.SinDelayed) / options.TOF_bins;
                    end
                end
            elseif options.scatter_correction
                if sum(size(options.ScatterC)) > 1 && ~iscell(options.ScatterC)
                    if size(options.ScatterC,1) ~= options.Ndist && ~options.use_raw_data && ~options.listmode
                        options.ScatterC = permute(options.ScatterC,[2 1 3 4]);
                    end
                    if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization && ~scatterAlreadyNormalized
                        options.ScatterC = normalizeScatter(options.ScatterC, options, normalizationInverted);
                    end
                    if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                        options.ScatterC = Randoms_variance_reduction(single(options.ScatterC), options);
                    end
                    if options.scatter_smoothing && ~ScatterProp.smoothing
                        options.ScatterC = randoms_smoothing(single(options.ScatterC), options);
                    end
                    if options.subtract_scatter
                        if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                            options.SinDelayed = single(scatterPerTOFBin(options.ScatterC(:), options));
                        else
                            options.SinDelayed = double(scatterPerTOFBin(options.ScatterC(:), options));
                        end
                    end
                elseif iscell(options.ScatterC)
                    options.SinDelayed = cell(length(options.ScatterC),1);
                    for kk = 1 : options.partitions
                        if sum(size(options.ScatterC{1})) > 1
                            if size(options.ScatterC{kk},1) ~= options.Ndist && ~options.use_raw_data && ~options.listmode
                                options.ScatterC{kk} = permute(options.ScatterC{kk},[2 1 3 4]);
                            end
                        end
                        if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization && ~scatterAlreadyNormalized
                            options.ScatterC{kk} = normalizeScatter(options.ScatterC{kk}, options, normalizationInverted);
                        end
                        if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                            options.ScatterC{kk} = Randoms_variance_reduction(single(options.ScatterC{kk}), options);
                        end
                        if options.scatter_smoothing && ~ScatterProp.smoothing
                            options.ScatterC{kk} = randoms_smoothing(single(options.ScatterC{kk}), options);
                        end
                        if options.subtract_scatter
                            if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                                options.SinDelayed{kk} = single(scatterPerTOFBin(options.ScatterC{kk}(:), options));
                            else
                                options.SinDelayed{kk} = double(scatterPerTOFBin(options.ScatterC{kk}(:), options));
                            end
                        end
                    end
                end
            end
        end
        if options.sampling > 1 && options.ordinaryPoisson
            options.SinDelayed = reshape(options.SinDelayed, options.Ndist, options.Nang, options.NSinos);
            options.SinDelayed = interpolateSinog(options.SinDelayed, options.sampling, options.Ndist, options.partitions, options.sampling_interpolation_method);
            options.SinDelayed = options.SinDelayed(:);
        end
        % Mark the scatter as normalized, so that another call with the returned options does not normalize it again
        if options.scatter_correction && options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization ...
                && ~scatterAlreadyNormalized && isfield(options,'ScatterC')
            if iscell(options.ScatterC)
                nScatter = numel(options.ScatterC{1});
            else
                nScatter = numel(options.ScatterC);
            end
            if nScatter > 1
                options.scatterNormalized = true;
            end
        end
    elseif (options.randoms_correction || options.scatter_correction) && ~options.reconstruct_trues && ~options.reconstruct_scatter
        r_exist = isfield(options,'SinDelayed');
        if r_exist && isfield(options,'SinDelayed') && ~iscell(options.SinDelayed) && isscalar(options.SinDelayed) && options.randoms_correction
            r_exist = false;
        elseif r_exist && isfield(options,'SinDelayed') && iscell(options.SinDelayed) && isscalar(options.SinDelayed{1}) && options.randoms_correction
            r_exist = false;
        end
        options.scatter = false;
        if exist('RandProp','var') == 0 || isempty(RandProp)
            RandProp.smoothing = false;
            RandProp.variance_reduction = false;
        end
        if exist('ScatterProp','var') == 0 || isempty(ScatterProp)
            ScatterProp.smoothing = false;
            ScatterProp.variance_reduction = false;
            ScatterProp.normalization = false;
        end
        if r_exist && options.randoms_correction
            if (~options.variance_reduction && RandProp.variance_reduction) || (~options.randoms_smoothing && RandProp.smoothing)
                %             if options.partitions == 1
                %                 options.SinDelayed = loadStructFromFile(sinoFile,'raw_SinDelayed');
                %             else
                options.SinDelayed = loadStructFromFile(sinoFile, 'raw_SinDelayed');
                %             end
                RandProp.variance_reduction = false;
                RandProp.smoothing = false;
            end
        end
        if options.randoms_correction && ~r_exist
            try
                options.SinDelayed = loadStructFromFile(sinoFile, 'SinDelayed');
            catch
                options = loadDelayedData(options);
            end
        end
        if options.randoms_correction
            [options,~, ~] = applyScatterRandoms(options,RandProp, options.SinDelayed, 0);
        end
        % randoms_correction = false;
        options.randoms_correction = false;
        if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
            options.SinDelayed = {single(0)};
        else
            options.SinDelayed = {0};
        end
        if options.scatter_correction
            [options,~, ~] = applyScatterRandoms(options,ScatterProp, options.ScatterC, 1);
            options = rmfield(options,'ScatterC');
            options.scatter_correction = false;
        end
        if options.verbose
            disp('Precorrections completed')
        end
    else
        % randoms_correction = false;
        options.randoms_correction = false;
        options.scatter_correction = false;
        options.scatter = false;
        if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
            options.SinDelayed = {single(0)};
        else
            options.SinDelayed = {0};
        end
    end
    % The normalization (options.normalization) is loaded and inverted before the randoms and scatter correction above, see
    % normalizationInverted. A file-loaded normalization of sinogram data is cut and interpolated here, after the scatter has been
    % normalized with it, as the scatter has not been cut or interpolated at that stage. Index-based data has been converted
    % to one value for each event already
    if normalizationInverted && ~options.useIndexBasedReconstruction
        if ~options.use_raw_data && options.NSinos ~= options.TotSinos
            options.normalization = options.normalization(1 : options.Ndist * options.Nang *options.NSinos);
        end
        if options.sampling > 1 && options.corrections_during_reconstruction
            options.normalization = reshape(options.normalization, options.Ndist, options.Nang, options.NSinos);
            options.normalization = interpolateSinog(options.normalization, options.sampling, options.Ndist, options.partitions, options.sampling_interpolation_method);
            options.normalization = options.normalization(:);
        end
    end
    if ctNormMode == 2
        % CT precorrection: the measurements are divided by the normalization (before the linearization in reconstructions_main).
        % The normalization correction is disabled afterwards, i.e. the normalization is not used by the projectors
        options = ctNormalizationPrecorrect(options);
    end
    if options.normalization_correction && ~options.corrections_during_reconstruction
        if ~isfield(options,'normalization') || isempty(options.normalization)
            if options.use_user_normalization
                [file, fpath] = selectNormalizationFile(options);
                [~, ~, normExt] = fileparts(file);
                if any(strfind(file, '.nrm'))
                    fid = fopen(fullfile(fpath, file));
                    options.normalization = fread(fid, inf, 'single=>single',0,'l');
                    fclose(fid);
                    options.InveonNorm = true;
                elseif any(strcmpi(normExt, {'.n', '.ptd'}))
                    % Biograph (mCT/Vision) normalization component file, the coefficients multiply the data
                    options.normalization = loadNormalizationBiograph(fullfile(fpath, file), normalizationOptions(options, normDims));
                    % The coefficients cover all sinograms (TotSinos) and have no TOF dimension. They are cropped
                    % to NSinos below and applied to every TOF bin
                    if numel(options.normalization) ~= normDims(1) * normDims(2) * normDims(4)
                        error('Size mismatch between the current data and the normalization datafile')
                    end
                else
                    data = load(fullfile(fpath, file));
                    variables = fieldnames(data);
                    if any(strcmp(variables,'normalization'))
                        options.normalization = data.(variables{strcmp(variables,'normalization')});
                    else
                        options.normalization = data.(variables{1});
                    end
                    clear data
                    % The size of index-based data is checked when the normalization is converted
                    if ~options.useIndexBasedReconstruction
                        if iscell(options.SinM)
                            if numel(options.normalization) ~= numel(options.SinM{1})
                                error('Size mismatch between the current data and the normalization datafile')
                            end
                        else
                            if numel(options.normalization) ~= numel(options.SinM)
                                error('Size mismatch between the current data and the normalization datafile')
                            end
                        end
                    end
                end
                options.normalization = options.normalization(:);
            else
                if options.use_raw_data
                    norm_file = [folder options.machine_name '_normalization_listmode.mat'];
                else
                    norm_file = [folder options.machine_name '_normalization_' num2str(normDims(1)) 'x' num2str(normDims(2)) '_span' num2str(options.span) '.mat'];
                end
                if exist(norm_file, 'file') == 2
                    options.normalization = loadStructFromFile(norm_file,'normalization');
                else
                    error('Normalization correction selected, but no normalization data found')
                end
                options.normalization = options.normalization(:);
            end
            if ~options.use_raw_data && options.NSinos ~= options.TotSinos
                options.normalization = options.normalization(1 : options.Ndist * options.Nang *options.NSinos);
            end
        end
        if options.useIndexBasedReconstruction
            % One value for each event instead of a sinogram, see indexBasedNormalization. The precorrected data does not
            % use the normalization in the sensitivity image. The normalization multiplies the data, i.e. it is a correction factor
            options = indexBasedNormalization(options, false, true);
        end
        NN = options.Ndist * options.Nang * options.NSinos;
        if iscell(options.SinM)
            for kk = 1 : options.partitions
                options.SinM{kk} = options.SinM{kk}(:);
                if options.TOF_bins > 1
                    % Convert to floating point first, otherwise the TOF bins would be rounded to the integer type of the raw data
                    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                        options.SinM{kk} = single(full(options.SinM{kk}));
                    else
                        options.SinM{kk} = double(full(options.SinM{kk}));
                    end
                    for uu = 1 : options.TOF_bins
                        if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                            options.SinM{kk}(NN * (uu - 1) + 1: NN * uu) = single(full(options.SinM{kk}(NN * (uu - 1) + 1: NN * uu))) .* single(options.normalization);
                        else
                            options.SinM{kk}(NN * (uu - 1) + 1: NN * uu) = full(options.SinM{kk}(NN * (uu - 1) + 1: NN * uu)) .* double(options.normalization);
                        end
                    end
                else
                    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                        options.SinM{kk} = single(full(options.SinM{kk})) .* single(options.normalization);
                    else
                        options.SinM{kk} = full(options.SinM{kk}) .* double(options.normalization);
                    end
                end
            end
        else
            if options.TOF_bins > 1
                % Convert to floating point first, otherwise the TOF bins would be rounded to the integer type of the raw data
                if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                    options.SinM = single(full(options.SinM(:)));
                else
                    options.SinM = double(full(options.SinM(:)));
                end
                for uu = 1 : options.TOF_bins
                    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                        options.SinM(NN * (uu - 1) + 1: NN * uu) = single(full(options.SinM(NN * (uu - 1) + 1: NN * uu))) .* single(options.normalization);
                    else
                        options.SinM(NN * (uu - 1) + 1: NN * uu) = full(options.SinM(NN * (uu - 1) + 1: NN * uu)) .* double(options.normalization);
                    end
                end
            else
                options.SinM = options.SinM(:);
                if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                    options.SinM = single(full(options.SinM)) .* single(options.normalization);
                else
                    options.SinM = full(options.SinM) .* double(options.normalization);
                end
            end
        end
        if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
            options.normalization = single(0);
        else
            options.normalization = 0;
        end
    elseif ~options.normalization_correction
        if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
            options.normalization = single(0);
        else
            options.normalization = 0;
        end
    end
    if options.normalization_correction && options.corrections_during_reconstruction && numel(options.normalization) > 1
        if iscell(options.SinM)
            SinM_numel = numel(options.SinM{1});
        else
            SinM_numel = numel(options.SinM);
        end
        if mod(SinM_numel, numel(options.normalization)) ~= 0
            warning('Normalization coefficient vector/matrix is of different size than the measurement data. Normalization might not work correctly and might cause a crash.')
        end
    end
    if options.arc_correction && ~options.precompute_lor
        [x, y, options] = arcCorrection(options, true);
    end
    if options.sampling > 1 && ~options.precompute_lor
        [~, ~, options] = increaseSampling(options, x, y, true);
    end
end

% Other SPECT corrections
if options.SPECT
    if options.scatter_correction && numel(options.SinDelayed) <= 1 && options.subtract_scatter % See 10.1088/0031-9155/56/14/R01
        if iscell(options.SinM) % SinM is cell (size = options.partitions)
            if options.corrections_during_reconstruction && ~iscell(options.SinDelayed)
                options.SinDelayed = cell(options.partitions,1);
            end
            for timestep = 1:options.partitions
                if numel(options.ScatterC) == 1 % DEW
                    k = 1;
                    if ~options.corrections_during_reconstruction
                        options.SinM{timestep} = options.SinM{timestep} - k * options.ScatterC{1}{timestep};
                        options.scatter_correction = false;
                    else
                        options.SinDelayed{timestep} = k * squeeze(options.ScatterC{1}{timestep});
                        if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                            options.SinDelayed{timestep} = single(options.SinDelayed{timestep} );
                        end
                    end
                elseif numel(options.ScatterC) == 2 % TEW
                    kLower = diff(options.eWin) / diff(options.eWinL);
                    kUpper = diff(options.eWin) / diff(options.eWinU);
                    if ~options.corrections_during_reconstruction
                        options.SinM{timestep} = options.SinM{timestep} - 0.5 * (kLower * squeeze(options.ScatterC{1}{timestep}) + kUpper * squeeze(options.ScatterC{2}{timestep}));
                        options.scatter_correction = false;
                    else
                        options.SinDelayed{timestep} = 0.5 * (kLower * squeeze(options.ScatterC{1}{timestep}) + kUpper * squeeze(options.ScatterC{2}{timestep}));
                        if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                            options.SinDelayed{timestep} = single(options.SinDelayed{timestep});
                        end
                    end
                end
            end
        else % SinM is not cell
            if numel(options.ScatterC) == 1 % DEW
                k = 1;
                if ~options.corrections_during_reconstruction
                    options.SinM = options.SinM - k * squeeze(options.ScatterC{1});
                    options.scatter_correction = false;
                else
                    options.SinDelayed = k * squeeze(options.ScatterC{1});
                    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                        options.SinDelayed = single(options.SinDelayed);
                    end
                    % randoms_correction = true;
                end
            elseif numel(options.ScatterC) == 2 % TEW
                kLower = diff(options.eWin) / diff(options.eWinL);
                kUpper = diff(options.eWin) / diff(options.eWinU);
                if ~options.corrections_during_reconstruction
                    options.SinM = options.SinM - 0.5 * (kLower * squeeze(options.ScatterC{1}) + kUpper * squeeze(options.ScatterC{2}));
                    options.scatter_correction = false;
                else
                    options.SinDelayed = 0.5 * (kLower * squeeze(options.ScatterC{1}) + kUpper * squeeze(options.ScatterC{2}));
                    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
                        options.SinDelayed = single(options.SinDelayed);
                    end
                    % randoms_correction = true;
                end
            end
        end
    end

    %error("break")
    if options.normalization_correction && options.corrections_during_reconstruction
        if numel(options.normalization) <= 1
            error('Normalization correction selected, but no normalization data input!')
        end
        % normalization_correction = true;
    elseif options.normalization_correction && ~options.corrections_during_reconstruction
        if numel(options.normalization) == numel(options.SinM)
            options.SinM = options.SinM ./ options.normalization;
        else
            if numel(options.normalization) <= 1
                error('Normalization correction selected, but no normalization data input!')
            end
            options.normalization = reshape(options.normalization, size(options.SinM));
            options.SinM = options.SinM ./ options.normalization;
        end
        % normalization_correction = false;
    end

    % randoms_correction = false;
end

if ~isfield(options,'global_correction_factor') || isempty(options.global_correction_factor) || options.global_correction_factor <= 0
    options.global_correction_factor = 1;
end
if ~isfield(options, 'ScatterC')
    if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
        options.ScatterC = {single(0)};
    else
        options.ScatterC = 0;
    end
end
if ~isfield(options, 'scatter')
    options.scatter = false;
end
if options.scatter_correction && ~options.subtract_scatter && ((iscell(options.ScatterC) && numel(options.ScatterC{1}) > 1) || (~iscell(options.ScatterC) && numel(options.ScatterC) > 1))
    options.scatter = true;
else
    options.scatter = false;
end

if options.corrections_during_reconstruction
    if numel(options.normalization) > 1
        options.normalization_correction = true;
    end
    if isfield(options,'SinDelayed') && ((iscell(options.SinDelayed) && numel(options.SinDelayed{1}) > 1) || (~iscell(options.SinDelayed) && numel(options.SinDelayed) > 1)) && options.ordinaryPoisson
        options.randoms_correction = true;
    end
else
    options.randoms_correction = false;
    options.scatter_correction = false;
    options.normalization_correction = false;
end
if ~options.randoms_correction
    options.randomsTOF = false;
end

function [options, ctNormMode] = ctNormalizationSetup(options)
% Checks the CT normalization (options.normalization_correction is true) and selects how it is applied. The CT normalization is a
% multiplicative factor of each measurement, input only as the array options.normalization: one map of the detector or one value for
% each measurement (the map is tiled to all projections here, before the subset permutation).
% The output ctNormMode is
%   1 (host gain): intensity data, a count-based algorithm (OSEM, etc.) and corrections during reconstruction. The normalization is a
%     detector gain of the blank scan, ybar = flat * n * exp(-l), and it is used by the native library (permuted to the subset order
%     like the measurements)
%   2 (precorrection): FDK, corrections applied to the data or a linear algorithm (LSQR, PDHG, etc.) with intensity data. The
%     measurements are divided by the normalization before the linearization, the normalization is not used by the projectors
%   3 (kernel): pre-linearized data with a linear algorithm (also a projection operator without an algorithm). The normalization
%     is used in the projector kernels, y = n * (A x)
% The randoms/scatter convention of CT is unchanged: the randoms and scatter estimates are in units of the flat value.
if ~isfield(options, 'normalization') || numel(options.normalization) <= 1
    error(['CT normalization correction is selected (options.normalization_correction = true), but options.normalization is empty. The CT normalization is ' ...
        'not loaded from a file, input it as an array of the detector size [' num2str(options.nRowsD) ' ' num2str(options.nColsD) '] or the size of all projections [' ...
        num2str(options.nRowsD) ' ' num2str(options.nColsD) ' ' num2str(options.nProjections) ']'])
end
options = tileCTNormalization(options);
if options.Nt > 1
    error('CT normalization is supported only with static data (one time frame)')
end
if options.largeDim
    error('CT normalization is not supported with largeDim')
end
if options.useEFOV
    error('CT normalization is not supported with extended FOV (options.useEFOV)')
end
if options.normalize_scatter
    error('CT normalization is not supported with options.normalize_scatter')
end
% Selected algorithms
algNames = recNames();
active = {};
for kk = 1 : numel(algNames)
    if isfield(options, algNames{kk}) && any(options.(algNames{kk})(:))
        active{end + 1} = algNames{kk};
    end
end
% Algorithms that use the exponential (intensity) model exp(-FP) * flat. All other algorithms are treated as linear
countNames = {'OSEM';'ROSEM';'OSL_OSEM';'ROSEM_MAP';'RAMLA';'BSREM';'RBI';'OSL_RBI';'DRAMA';'MBSREM';'MRAMLA';'SPS';'PKMA';'SAGA'};
isCount = ismember(active, countNames);
anyCount = any(isCount);
anyLinear = any(~isCount);
countData = ~options.usingLinearizedData;
if ~isempty(active)
    if options.implementation ~= 2 || options.use_CPU
        error('CT normalization is supported only with implementation 2 (OpenCL/CUDA/HIP), not with implementation 1, 3 or 4 or with options.use_CPU')
    end
    if anyCount && anyLinear
        error('CT normalization supports only one algorithm at a time, but both count-based and linear algorithms are selected')
    end
    if anyCount && ~countData
        error('CT normalization with count-based algorithms needs intensity data (options.usingLinearizedData = false), the linearized data cannot be used')
    end
end
useFDK = any(strcmp(active, 'FDK'));
if ~options.corrections_during_reconstruction || useFDK || (anyLinear && countData)
    ctNormMode = 2;
    normMeas = double(options.normalization(:));
    if any(~isfinite(normMeas)) || any(normMeas <= 0)
        error(['CT normalization has non-positive or non-finite values, but the data is divided by the normalization. Set the normalization of ' ...
            'dead or missing detector pixels to a positive value and exclude the pixels with the forward projection mask (options.maskFP)'])
    end
elseif anyCount
    ctNormMode = 1;
else
    ctNormMode = 3;
    if mod(options.projector_type, 10) == 5
        error('CT normalization with linearized data is not supported with backprojection type 5')
    end
end

function options = ctNormalizationPrecorrect(options)
% CT precorrection: divides the measurements by the normalization and disables the normalization correction (the normalization
% is not used by the projectors afterwards)
if options.implementation == 2 || options.implementation == 3 || options.implementation == 5 || options.useSingles
    normMeas = single(options.normalization(:));
else
    normMeas = double(options.normalization(:));
end
if iscell(options.SinM)
    for kk = 1 : numel(options.SinM)
        if numel(options.SinM{kk}) ~= numel(normMeas)
            error('Size mismatch between the CT measurements and the CT normalization')
        end
        options.SinM{kk} = reshape(cast(full(options.SinM{kk}(:)), 'like', normMeas) ./ normMeas, size(options.SinM{kk}));
    end
else
    if numel(options.SinM) ~= numel(normMeas)
        error('Size mismatch between the CT measurements and the CT normalization')
    end
    options.SinM = reshape(cast(full(options.SinM(:)), 'like', normMeas) ./ normMeas, size(options.SinM));
end
options.normalization_correction = false;
options.normalization = zeros(1, 1, 'like', normMeas);

function [file, fpath] = selectNormalizationFile(options)
% Returns the user normalization file: options.normalization_datafile if it is set, otherwise the user
% is prompted for the file
if isfield(options, 'normalization_datafile') && ~isempty(options.normalization_datafile)
    normFile = options.normalization_datafile;
    if exist(normFile, 'file') ~= 2
        error(['The normalization datafile ' normFile ' does not exist'])
    end
    [fpath, name, ext] = fileparts(normFile);
    file = [name ext];
    if isempty(fpath)
        % Relative file name, the file is in the current folder or on the MATLAB path
        fullName = which(normFile);
        if isempty(fullName)
            fullName = fullfile(pwd, normFile);
        end
        fpath = fileparts(fullName);
    end
else
    [file, fpath] = uigetfile({'*.nrm;*.n;*.ptd;*.mat'},'Select normalization datafile');
    if isequal(file, 0)
        error('No file was selected')
    end
end

function options = normalizationOptions(options, normDims)
% Returns the options with the sinogram dimensions of the normalization data (normDims = [Ndist, Nang, NSinos, TotSinos]),
% which index-based data has collapsed, for the functions that read the sinogram normalization
options.Ndist = normDims(1);
options.Nang = normDims(2);
options.NSinos = normDims(3);
options.TotSinos = normDims(4);

function tf = scatterTOFApplies(options)
% True if the size rules of the scatter data with TOF data apply, i.e. the data is a PET sinogram with TOF bins. Raw data,
% index-based reconstruction, list-mode data and SPECT keep their behavior
tf = options.TOF_bins > 1 && ~options.use_raw_data && ~options.useIndexBasedReconstruction && ~options.listmode && ~options.SPECT;

function options = checkTOFScatter(options)
% Checks the size of the scatter data with TOF data when the scatter is subtracted during the reconstruction (ordinary
% Poisson). A non-TOF scatter sinogram (Ndist * Nang * NSinos elements) is repeated for all the TOF bins and divided evenly
% among them (see scatterPerTOFBin), a message about this is printed here, once for all the time frames. TOF-resolved scatter data
% (Ndist * Nang * (NSinos or TotSinos) * TOF_bins elements, the TOF bin is the slow index as with the measurements) is used as it
% is in each TOF bin and options.randomsTOF is set (the randoms are divided evenly among the TOF bins and added to it, see
% addScatterPerTOFBin). It is supported only with implementation 2, one time frame, no increased sampling, and no scatter smoothing
% or variance reduction. The non-subtracted scatter (corrVector) is not affected.
options.randomsTOF = false;
if ~scatterTOFApplies(options) || ~options.scatter_correction || ~options.subtract_scatter || ~isfield(options, 'ScatterC')
    return
end
if iscell(options.ScatterC)
    scatterData = options.ScatterC;
else
    scatterData = {options.ScatterC};
end
sinogramSizes = options.Ndist * options.Nang * unique([options.NSinos, options.TotSinos]);
noTOF = false;
tofSized = false;
for kk = 1 : numel(scatterData)
    nScatter = numel(scatterData{kk});
    if nScatter <= 1
        continue
    end
    if any(nScatter == sinogramSizes)
        noTOF = true;
    elseif any(nScatter == sinogramSizes * options.TOF_bins)
        tofSized = true;
    else
        error(['The scatter data has ' num2str(nScatter) ' elements, which does not match the TOF data. With TOF the scatter data has to be a non-TOF sinogram of ' ...
            num2str(options.Ndist * options.Nang * options.NSinos) ' elements (Ndist * Nang * NSinos), or a TOF-resolved sinogram of ' ...
            num2str(options.Ndist * options.Nang * options.NSinos * options.TOF_bins) ' elements (Ndist * Nang * NSinos * TOF_bins)'])
    end
end
if tofSized && noTOF
    error('The scatter data must be either TOF-resolved or non-TOF in all the time frames')
end
if tofSized
    if numel(scatterData) > 1 || any(options.partitions(:) > 1) || options.Nt > 1
        error('TOF-resolved scatter data with corrections during reconstruction is supported only with static data (one time frame); use precorrection for dynamic data')
    end
    if options.implementation ~= 2
        error('TOF-resolved scatter data with corrections during reconstruction is supported only with implementation 2; use precorrection or a non-TOF scatter sinogram')
    end
    if options.sampling > 1
        error('TOF-resolved scatter data with corrections during reconstruction is not supported with increased sampling (options.sampling > 1)')
    end
    if options.scatter_smoothing || options.scatter_variance_reduction
        error('Scatter smoothing and scatter variance reduction are not supported with TOF-resolved scatter data; disable them or use a non-TOF scatter sinogram')
    end
    if options.largeDim
        error('TOF-resolved scatter data with corrections during reconstruction is not supported with largeDim')
    end
    options.randomsTOF = true;
    disp(['The scatter data is TOF-resolved and is used as it is in each of the ' num2str(options.TOF_bins) ' TOF bins (the randoms are divided evenly among the TOF bins)'])
elseif noTOF
    disp(['The scatter data has no TOF bins and is divided evenly among the ' num2str(options.TOF_bins) ' TOF bins'])
end

function scatterData = scatterPerTOFBin(scatterData, options)
% Returns the scatter data as it is added to the additive term of each TOF bin: the non-TOF scatter sinogram is divided evenly
% among the TOF bins (the size has been checked in checkTOFScatter), TOF-resolved scatter data (options.randomsTOF) is used as it is.
% options.ScatterC itself is not modified.
if scatterTOFApplies(options) && numel(scatterData) > 1 && ~options.randomsTOF
    scatterData = scatterData / options.TOF_bins;
end

function randoms = addScatterPerTOFBin(randoms, scatterData)
% Adds the scatter data of scatterPerTOFBin to the randoms that have already been divided among the TOF bins. With TOF-resolved
% scatter data the randoms are repeated for all the TOF bins first (TOF bin is the slow index, as with the TOF measurements)
nRandoms = numel(randoms);
nScatter = numel(scatterData);
if nScatter > nRandoms && nRandoms > 0 && mod(nScatter, nRandoms) == 0
    randoms = repmat(randoms(:), nScatter / nRandoms, 1);
end
randoms = randoms + scatterData;

function scatterData = normalizeScatter(scatterData, options, normalizationInverted)
% Multiplies the scatter data with the normalization factor (see normalizationScatterFactor). The normalization applies to the
% LORs, so TOF-resolved scatter data (options.randomsTOF) is multiplied with it in each TOF bin
nBins = 1;
if scatterTOFApplies(options) && options.randomsTOF
    nBins = options.TOF_bins;
end
factor = normalizationScatterFactor(options, normalizationInverted, numel(scatterData) / nBins);
if nBins > 1
    scatterData = reshape(reshape(single(scatterData), numel(factor), nBins) .* factor(:), size(scatterData));
else
    scatterData = single(scatterData) .* reshape(factor, size(scatterData, 1), size(scatterData, 2), size(scatterData, 3));
end

function scatterData = cropScatterSinograms(scatterData, options)
% Keeps the first NSinos sinograms of the scatter data that has TotSinos sinograms. TOF-resolved scatter data (options.randomsTOF)
% is cropped in each TOF bin, or left as it is if it has NSinos sinograms already
nCrop = options.Ndist * options.Nang * options.NSinos;
if scatterTOFApplies(options) && options.randomsTOF
    nTotal = options.Ndist * options.Nang * options.TotSinos;
    if numel(scatterData) == nTotal * options.TOF_bins
        scatterData = reshape(scatterData, nTotal, options.TOF_bins);
        scatterData = reshape(scatterData(1 : nCrop, :), [], 1);
    end
else
    scatterData = scatterData(1 : nCrop);
end

function factor = normalizationScatterFactor(options, normalizationInverted, nScatter)
% Returns the factor that the scatter is multiplied with when it is normalized. The normalization is a multiplier in the
% reconstruction (it has been inverted if it was loaded from a file in this call), so the scatter is multiplied with its
% reciprocal. Only a prefilled Inveon normalization (raw coefficients) with Inveon scatter is used as it is. A normalization
% loaded from a file has all the sinograms (TotSinos), as it is cut to NSinos only after the scatter has been normalized, so
% the first nScatter factors are used with a scatter that has NSinos sinograms.
if ~normalizationInverted && isfield(options,'InveonScatter') && isfield(options,'InveonNorm')
    factor = options.normalization;
else
    factor = reciprocalOrOne(options.normalization);
end
if normalizationInverted && nScatter < numel(factor)
    factor = factor(1 : nScatter);
end

function r = reciprocalOrOne(x)
% Returns 1 ./ x, except that the elements where x is zero are one. The scatter is multiplied with this reciprocal, so the
% scatter keeps its value where the normalization is zero (the zero normalization does not divide the scatter).
r = ones(size(x), 'like', x);
nz = x ~= 0;
r(nz) = 1 ./ x(nz);
