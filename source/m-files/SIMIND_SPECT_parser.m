function options = SIMIND_SPECT_parser(fname)
    % Load header
    options.fpath = strcat(fname, '.h00');
    fid = fopen(options.fpath);
    hdr = textscan(fid,'%s','Delimiter','=');
    hdr = hdr{1};
    fclose(fid);

    % Number of projections
    ind = find(contains(hdr, '!number of projections')) + 1;
    options.nProjections = str2double(hdr{ind});
    
    % Start angle
    ind = find(contains(hdr, 'start angle')) + 1;
    options.startAngle = str2double(hdr{ind});

    % Number of detector heads
    ind = find(contains(hdr, 'number of detector heads')) + 1;
    options.nHeads = str2double(hdr{ind});

    % Extent of rotation
    ind = find(contains(hdr, 'extent of rotation')) + 1;
    extRot = str2double(hdr{ind});

    % Angle increment per view
    angleIncrement = extRot / options.nProjections;

    % Projection angles
    options.angles = (0:(options.nProjections - 1)) * angleIncrement + options.startAngle;

    % Detector pixel count X
    ind = find(contains(hdr, '!matrix size [1]')) + 1;
    matrixSizeX = str2double(hdr{ind});

    % Detector pixel count Y
    ind = find(contains(hdr, '!matrix size [2]')) + 1;
    matrixSizeY = str2double(hdr{ind});

    % Crystal size XY
    ind = find(contains(hdr, 'scaling factor (mm/pixel) [1]')) + 1;
    options.crXY = str2double(hdr{ind});

    % Crystal thickness (mm)
    ind = find(contains(hdr, ';# Crystal Thickness')) + 1;
    options.cr_p = str2double(hdr{ind});
    
    % Collimator hole length (mm)
    ind = find(contains(hdr, ';# Collimator thickness')) + 1;
    options.collimatorLength = str2double(hdr{ind});
    options.colL = str2double(hdr{ind});

    % Collimator hole radius (mm), larger inner radius
    ind = find(contains(hdr, ';# Collimator hole diameter')) + 1;
    options.colR = 0.5 * str2double(hdr{ind});

    % Septal thickness (mm)
    ind = find(contains(hdr, ';# Collimator hole septa')) + 1;  
    options.dSeptal = str2double(hdr{ind});

    % Intrinsic resolution
    ind = find(contains(hdr, ';# Intrinsic FWHM for the camera')) + 1;  
    options.iR = str2double(hdr{ind});

    % Load detector radius per projection
    options.radiusPerProj = load(strcat(fname, '.cor'));
    options.radiusPerProj = options.radiusPerProj(:, 1) .* 10;

    % Load data
    DataPath = strcat(fname, '.a00');
    fid = fopen(DataPath);
    data = fread(fid, inf, "float32");
    fclose(fid);
    options.SinM = reshape(data, matrixSizeX, matrixSizeY, options.nProjections);

    % Number of rows in a projection image
    options.nRowsD = size(options.SinM, 1);
    
    % Number of columns in a projection image
    options.nColsD = size(options.SinM, 2);

end