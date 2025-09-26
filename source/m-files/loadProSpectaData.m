function options = loadProSpectaData(options)
% LOADPROSPECTADATA Loads Siemens Pro.Specta data.
    fprintf("Loading Pro.Specta projection data... ")
    info = dicominfo(options.fpath);

    %% FOV alignment
    dz = info.PixelSpacing(1);
    Nz = double(info.Rows);
    sZ = info.DetectorInformationSequence.Item_1.ImagePositionPatient(3); % Z coordinate of the first sinogram pixel
    cZ = sZ - 0.5*dz*(Nz-1); % Center of sinogram Z in world coordinates

    tableHeight = info.RotationInformationSequence.Item_1.TableHeight;
    %tableTraverse = info.RotationInformationSequence.Item_1.TableTraverse;

    % SPECT world-limits along each axis
    xLimits = options.FOVa_x * [-0.5, 0.5];
    yLimits = options.FOVa_y * [-0.5, 0.5] - tableHeight;
    zLimits = options.axial_fov * [-0.5, 0.5] + cZ;

    % build the SPECT imref3d
    options.refSPECT = imref3d([options.Nx, options.Ny, options.Nz], xLimits, yLimits, zLimits);

    %% Projection images
    options.SinM = squeeze(dicomread(options.fpath));
    options.SinM = rot90(options.SinM, 3);

    % Number of rows in a projection image
    options.nRowsD = size(options.SinM, 1);

    % Number of columns in a projection image
    options.nColsD = size(options.SinM, 2);

    % Total number of projections
    options.nProjections = size(options.SinM,3);

    % Number of detector heads
    options.nHeads = 2;

    %% Projection angles
    % Starting angles for both heads
    startAngle1 = info.DetectorInformationSequence.Item_1.StartAngle;
    startAngle2 = info.DetectorInformationSequence.Item_2.StartAngle;

    % Increment value for the projection angles
    angleIncrement = info.RotationInformationSequence.Item_1.AngularStep;

    % Additional offangle for CT
    additionalOffAngleCT = -5;

    % Actual angles
    options.angles = -[(startAngle1 : angleIncrement : startAngle1 + angleIncrement * (options.nProjections / options.nHeads - 1))';(startAngle2 : angleIncrement : startAngle2 + angleIncrement * (options.nProjections / options.nHeads - 1))'] + additionalOffAngleCT;

    %% Projection distances
    % Distances from the panel to the center of rotation
    options.radiusPerProj = [info.DetectorInformationSequence.Item_1.RadialPosition;info.DetectorInformationSequence.Item_2.RadialPosition];


    fprintf("Ready\n")
end