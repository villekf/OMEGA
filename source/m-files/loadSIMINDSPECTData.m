function options = loadSIMINDSPECTData(options)
%LOADSIMINDSPECTDATA Reads voxel-based SIMIND SPECT data to be reconstructed by OMEGA. 
    %   Utility function for OMEGA
    contains = @(str, pattern) ~cellfun('isempty', strfind(str, pattern));

    % Load header
    fid = fopen(strcat(options.fpath, '.h00'));
    hdr = textscan(fid,'%s','Delimiter','=');
    hdr = hdr{1};
    fclose(fid);

    % Energy window
    ind1 = find(contains(hdr, ';energy window lower level')) + 1;
    ind2 = find(contains(hdr, ';energy window upper level')) + 1;
    options.eWin = [str2double(hdr{ind1}) str2double(hdr{ind2})];

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
    options.dPitchX = str2double(hdr{ind});
    options.dPitchY = str2double(hdr{ind});

    % Crystal thickness (mm)
    ind = find(contains(hdr, ';# Crystal Thickness')) + 1;
    options.cr_p = str2double(hdr{ind});

    % Collimator hole length (mm)
    ind = find(contains(hdr, ';# Collimator thickness')) + 1;
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
    options.radiusPerProj = readmatrix(strcat(options.fpathCor, '.cor'), FileType="text");
    options.radiusPerProj = 10 * options.radiusPerProj(:, 1);

    % Load data
    fid = fopen(strcat(options.fpath, '.a00'));
    data = fread(fid, inf, "float32");
    fclose(fid);
    options.SinM = reshape(data, matrixSizeX, matrixSizeY, options.nProjections);

    % Number of rows in a projection image
    options.nRowsD = size(options.SinM, 1);

    % Number of columns in a projection image
    options.nColsD = size(options.SinM, 2);

    if isfile(strcat(options.fpathCT, '.hct')) && isfile(strcat(options.fpathCT, '.ict')) % Attenuation map
        fid = fopen(strcat(options.fpathCT, '.hct'));
        hdrCT = textscan(fid,'%s','Delimiter','=');
        hdrCT = hdrCT{1};
        fclose(fid);

        ind = find(contains(hdrCT, '!matrix size [1]')) + 1;
        Nx = str2double(hdrCT{ind});

        ind = find(contains(hdrCT, '!matrix size [2]')) + 1;
        Ny = str2double(hdrCT{ind});

        ind = find(contains(hdrCT, '!matrix size [3]')) + 1;
        Nz = str2double(hdrCT{ind});

        fid = fopen(strcat(options.fpathCT, '.ict'));
        I = fread(fid, 'float32');
        fclose(fid);
        I = reshape(I, [Nx, Ny, Nz]);
        I = rot90(I, 3);
        I = flip(I, 1);
        options.vaimennus = .1 * single(I);
    end

    if isfile(strcat(options.fpathScatterLower, '.h00')) % Lower energy window
        dataFilePath = strcat(options.fpathScatterLower, '.a00');
        fid = fopen(dataFilePath);
        data = fread(fid, inf, "float32");
        fclose(fid);
        options.ScatterC{1} = reshape(data, matrixSizeX, matrixSizeY, options.nProjections);

        fid = fopen(strcat(options.fpathScatterLower, '.h00'));
        hdrScatterL = textscan(fid,'%s','Delimiter','=');
        hdrScatterL = hdrScatterL{1};
        fclose(fid);

        % Energy window
        ind1 = find(contains(hdrScatterL, ';energy window lower level')) + 1;
        ind2 = find(contains(hdrScatterL, ';energy window upper level')) + 1;
        options.eWinL = [str2double(hdrScatterL{ind1}) str2double(hdrScatterL{ind2})];
    end
    if isfile(strcat(options.fpathScatterUpper, '.h00')) % Upper energy window
        dataFilePath = strcat(options.fpathScatterUpper, '.a00');
        fid = fopen(dataFilePath);
        data = fread(fid, inf, "float32");
        fclose(fid);
        options.ScatterC{2} = reshape(data, matrixSizeX, matrixSizeY, options.nProjections);
        
        fid = fopen(strcat(options.fpathScatterUpper, '.h00'));
        hdrScatterU = textscan(fid,'%s','Delimiter','=');
        hdrScatterU = hdrScatterU{1};
        fclose(fid);

        % Energy window
        ind1 = find(contains(hdrScatterU, ';energy window lower level')) + 1;
        ind2 = find(contains(hdrScatterU, ';energy window upper level')) + 1;
        options.eWinU = [str2double(hdrScatterU{ind1}) str2double(hdrScatterU{ind2})];
    end


    % SIMIND cannot simulate swiveling detector heads:
    options.swivelAngles = options.angles + 180;
    options.CORtoDetectorSurface = 0;

    % SIMIND measures radius to top of collimator
    options.radiusPerProj = options.radiusPerProj + options.colD + options.colL;
end
