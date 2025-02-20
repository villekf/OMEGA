function options = loadProSpectaData(options)
% LOADPROSPECTADATA Loads Siemens Pro.Specta data. CT support requires MATLAB Medical Imaging toolbox.
    options = loadProjectionData(options);
    if options.attenuation_correction
        options = loadCTData(options);
    end

    function options = loadProjectionData(options)
        % Header file location
        info = dicominfo(options.fpath);
        options.SPECTinfo = info;
        % Load projection images
        options.SinM = squeeze(dicomread(options.fpath));

        options.SinM = permute(options.SinM, [2 1 3]);

        % options.SinM = options.SinM(:,64/4+1:128-64/4,:);

        % Number of rows in a projection image
        options.nRowsD = size(options.SinM, 1);

        % Number of columns in a projection image
        options.nColsD = size(options.SinM, 2);

        % Total number of projections
        options.nProjections = size(options.SinM,3);

        % Number of detector heads
        options.nHeads = 2;

        % Starting angles for both heads
        startAngle1 = info.DetectorInformationSequence.Item_1.StartAngle;
        startAngle2 = info.DetectorInformationSequence.Item_2.StartAngle;

        % Increment value for the projection angles
        angleIncrement = info.RotationInformationSequence.Item_1.AngularStep;

        % Projection angles
        options.angles = [(startAngle2 : angleIncrement : startAngle2 + angleIncrement * (options.nProjections / options.nHeads - 1))';(startAngle1 : angleIncrement : startAngle1 + angleIncrement * (options.nProjections / options.nHeads - 1))'];

        % Distances from the panel to the center of rotation
        options.radiusPerProj = [info.DetectorInformationSequence.Item_1.RadialPosition;info.DetectorInformationSequence.Item_2.RadialPosition];
    end

    function options = loadCTData(options)
        fprintf("Loading Pro.Specta CT images... ")

        m = dicominfo(options.fpath); % Load DICOM metadata
        sinogramZmax = m.DetectorInformationSequence.Item_1.ImagePositionPatient(3); % Maximum Z value (mm) of detector panel (pixel centre)

        % Load CT volume
        CTvol = medicalVolume(options.fpathCT); 

        % Calculate voxel sizes
        voxelSizeX = options.FOVa_x / options.Nx;
        voxelSizeY = options.FOVa_y / options.Ny;
        voxelSizeZ = options.axial_fov / options.Nz;
        
        refpos =  [ % CT volume correct positioning
            CTvol.VolumeGeometry.Position(1, 1),...
            CTvol.VolumeGeometry.Position(1, 2),...
            sinogramZmax
        ];
        refpos = repmat(refpos, options.Nz, 1);
        refpos(:,3) = refpos(:,3) - (0:options.Nz-1)' * options.crXY;
        %refpos = flipud(refpos);
        %disp(refpos)
        pixelSpacing = repmat([voxelSizeX, voxelSizeY], options.Nz, 1);
        cosines = repmat([1 0 0; 0 1 0], 1, 1, options.Nz);

        % Create ref object
        R = medicalref3d([options.Nx, options.Ny, options.Nz], refpos, pixelSpacing, cosines);
        R.PatientCoordinateSystem = "LPS+";

        % Resample (fill with HU=-1000)
        CTvol = resample(CTvol, R, Method='linear', FillValue=-1000);
        CTvol = CTvol.Voxels;

        % Assign to attenuation map variable
        options.vaimennus = single(CTvol);
        options.vaimennus = imrotate(options.vaimennus, options.offangle * 360 / (2*pi));

        % Map to linear attenuation coefficients
        options.vaimennus = options.vaimennus * options.HU.slope + options.HU.intercept;
        
        % Convert to single
        options.vaimennus = single(options.vaimennus);
        options.vaimennus(options.vaimennus < 0) = 0;

        fprintf("Ready\n")
    end
end