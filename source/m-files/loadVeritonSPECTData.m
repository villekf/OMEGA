function options = loadVeritonSPECTData(options)
    options = loadInfo(options); % Read imaging parameters
    options = loadSensitivityMaps(options); % Read sensitivity maps
    options = loadProjectionData(options); % Load sinograms
    options = rmfield(options, "maskFP"); % Remove FP mask
    options = loadBlockLocations(options); % Load BlockLocations.csv data
    if options.attenuation_correction
        options = readCT(options); % CT attenuation correction
    end
    
    function options = loadInfo(options)
        fprintf("Loading Veriton info\n")
        C_info = readcell(fullfile(options.fpath, options.infoFileName));
    
        % Variable parameters 
        options.numOrbits = C_info{12, 2};
        options.ellipsoid.x0 = C_info{62, 2}; % Center x (horizontal)
        options.ellipsoid.y0 = C_info{63, 2}; % Center y (horizontal)
        options.ellipsoid.a = C_info{64, 2}; % Minor axis
        options.ellipsoid.b = C_info{65, 2}; % Major axis
        options.ellipsoid.n = C_info{66, 2}; % Power
    
        % Fixed parameters
        options.nHeads = C_info{8, 2};
        options.nRowsD = 16;
        options.nColsD = 128;
    end
    
    function options = loadSensitivityMaps(options)
        fprintf("Loading Veriton sensitivity maps      ")
        options.maskFP = ones(options.nRowsD, options.nColsD, options.nHeads);
        SensitivityMapsTable = readtable(fullfile(options.fpath, options.sensitivityMapFileName));
        for ii = 1:options.nHeads
            fprintf("\b\b\b\b\b")
            fprintf("%2u/12", ii)
            sensMap = SensitivityMapsTable(ii, 2:end);
            sensMap = table2array(sensMap);
            sensMap = reshape(sensMap, 16, 128);
            options.maskFP(:,:,ii) = sensMap;
        end
        fprintf("\n")
    end
    
    function options = loadProjectionData(options)
        radiusPerProj = cell(options.nHeads, 1);
        SinM = cell(options.nHeads, 1);
        angles = cell(options.nHeads, 1);
        swivelAngles = cell(options.nHeads, 1);
    
        fprintf("Loading Veriton projection data      ")
        for ii = 0:options.nHeads-1
            fprintf("\b\b\b\b\b")
            fprintf("%2u/12", ii+1)
    
            fullpath = fullfile(options.fpath, strcat(num2str(ii), options.fname));
            T0 = readtable(fullpath);
    
            T0head = T0(:, 1:9);
            T0data = T0(:, 10:end-2);
            T0data = table2array(T0data);
            sinograms = zeros(16, 128, size(T0data, 1));
            for jj = 1:size(T0data, 1)
                sinogram = reshape(T0data(jj, :), options.nRowsD, options.nColsD); % Read and reshape sinogram
                if isfield(options, 'maskFP')
                    sinogram = sinogram ./ options.maskFP(:,:,ii+1); % Apply sensitivity map
                end
                sinograms(:, :, jj) = sinogram; % Add to sinogram array
            end
            SinM{ii+1} = sinograms;
            radiusPerProj{ii+1} = T0head.Var6(:);
            angles{ii+1} = T0head.Var7(:);
            swivelAngles{ii+1} = T0head.Var9(:);
        end
        fprintf("\n")
    
        options.radiusPerProj = cat(1, radiusPerProj{:});
        options.angles = cat(1, angles{:});
        options.swivelAngles = cat(1, swivelAngles{:});
        options.swivelAngles = options.swivelAngles + options.angles;
        options.swivelAngles = options.swivelAngles + 90;
        options.SinM = flipud(cat(3, SinM{:}));
        options.nProjections = size(options.SinM, 3);
    
        if isfield(options, 'maskFP')
            options.maskFP = uint8(logical(options.maskFP));
        end
    end
    
    function options = loadBlockLocations(options)
        fprintf("Loading Veriton detector locations\n")
        warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')
        T_blockLocations = readtable(fullfile(options.fpath, options.blockLocationsFileName));
        
        options.blockIndex = zeros(options.nProjections, 1);
        options.gantryIndex = zeros(options.nProjections, 1);
        options.homeAngles = zeros(options.nProjections, 1);
    
        idx = 1;
        nProjPerDetector = options.nProjections / options.nHeads;
        nProjPerOrbit = nProjPerDetector / options.numOrbits;
        for head = 0:options.nHeads-1
            for gantry = 0:options.numOrbits-1
                for proj = 0:nProjPerOrbit-1
                    options.blockIndex(idx) = head;
                    options.gantryIndex(idx) = gantry;
                    options.homeAngles(idx) = table2array(T_blockLocations(12*gantry+head+1, 5));
                    idx = idx + 1;
                end
            end
        end
        options.gantryHomeAngle = table2array(T_blockLocations(end, 2));
    end

    function options = readCT(options)
        fprintf("Loading Veriton CT images        ")
        listing = dir(options.fpathCT);
        
        CTvol = zeros(512, 512, numel(listing)-2);
        for ii = 1:numel(listing)
            fprintf("\b\b\b\b\b\b\b")
            fprintf("%3u/%3u", ii, numel(listing))
            f = listing(ii);
            if strcmp(f.name, ".") || strcmp(f.name, "..")
                continue
            end
            CTinfo = dicominfo(fullfile(f.folder, f.name));
            CTimg = dicomread(fullfile(f.folder, f.name));
            CTvol(:,:,ii) = CTimg;
        end
        fprintf("\n")
        
        % CT volume cropping
        scale = [options.FOVa_x options.FOVa_y options.axial_fov] / CTinfo.ReconstructionDiameter; % Scale of volume to be cropped
        win = centerCropWindow3d(size(CTvol), floor(scale.*size(CTvol))); % Crop window
        CTvol = imcrop3(CTvol, win);
        CTvol = flipud(CTvol);

        % Conversion from HU to linear attenuation coefficients
        options.vaimennus = 0.006 * (1 + CTvol./1000);
        options.vaimennus(options.vaimennus < 0) = 0;
        options.vaimennus = imresize3(options.vaimennus, [options.Nx, options.Ny, options.Nz], 'Method', 'linear');
        options.CT_attenuation = true;
        options.flipAttImageZ = 1;
    end
end