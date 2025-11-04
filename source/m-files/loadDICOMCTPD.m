function [proj, vars] = loadDICOMCTPD(path)
%LOADDICOMCTPD Automatically loads DICOM CT projection data
%   This file will automatically load the DICOM CT projection data
%   available from:
%   https://aapm.app.box.com/s/eaw4jddb53keg1bptavvvd1sf4x3pe9h/folder/144226105715
%   Loads both projections and the necessary variables/coordinates
%   The only input is the path to the folder containing the DICOM images
%   Output is the projections and a struct containing the
%   variables/corodinates
arguments (Input)
    path (1,:) char {mustBeFolder}
end

arguments (Output)
    proj single
    vars struct
end

files = dir(fullfile(path, '*dcm'));
nFiles = numel(files);
xs = zeros(nFiles,1,'single');
xd = zeros(nFiles,1,'single');
ys = zeros(nFiles,1,'single');
yd = zeros(nFiles,1,'single');
zs = zeros(nFiles,1,'single');
zd = zeros(nFiles,1,'single');
angles = zeros(nFiles,1,'single');
for kk = 1 : nFiles
    info = dicominfo(fullfile(path, files(kk).name));
    nColsD = info.Width;
    nRowsD = info.Height;
    r = typecast(info.Private_7031_1003, 'single'); % DetectorFocalCenterRadialDistance
    deltar = typecast(info.Private_7033_100d, 'single'); % SourceRadialDistanceShift
    rho = r + deltar; 
    angles(kk) = typecast(info.Private_7031_1001, 'single') - pi/2; % DetectorFocalCenterAngularPosition
    deltaphi = typecast(info.Private_7033_100b, 'single'); % SourceAngularPositionShift
    phi = angles(kk) + deltaphi;
    deltaz = typecast(info.Private_7033_100c, 'single'); % SourceAxialPositionShift
    zs(kk) = typecast(info.Private_7031_1002, 'single') + deltaz; % DetectorFocalCenterAxialPosition
    ys(kk) = rho * sin(phi);
    xs(kk) = -rho * cos(phi);
    % nProjections = typecast(info.Private_7033_1013, 'uint16'); % NumberofSourceAngularSteps
    data = dicomread(fullfile(path, files(kk).name));
    if kk == 1
        proj = zeros(nRowsD, nColsD, nFiles, 'single');
    end
    proj(:,:,kk) = fliplr(single(data) * info.RescaleSlope + info.RescaleIntercept);
    dPitchX = typecast(info.Private_7029_1002, 'single'); % DetectorElementTransverseSpacing
    dPitchY = typecast(info.Private_7029_1006, 'single'); % DetectorElementAxialSpacing
    L = double(nRowsD) * dPitchX;
    theta = L / r;
    dtheta = theta / double(nRowsD);
    detElements = typecast(info.Private_7031_1033, 'single'); % DetectorCentralElement
    rowStart = detElements(1) - double(nRowsD)/2 - .5;
    colStart = detElements(2) - double(nColsD)/2 - .5;
    d = typecast(info.Private_7031_1031, 'single'); % ConstantRadialDistance
    yd(kk) = r*sin(angles(kk)) - d * sin(angles(kk) + dtheta * colStart) - sin(rowStart * dtheta) * rowStart * dPitchX;
    xd(kk) = -r*cos(angles(kk)) + d * cos(angles(kk) + dtheta * rowStart) - cos(rowStart * dtheta) * rowStart * dPitchX;
    zd(kk) = typecast(info.Private_7031_1002, 'single') + colStart * dPitchY;
end
nProjections = nFiles;
vars.xs = xs;
vars.ys = ys;
vars.zs = zs;
vars.xd = xd;
vars.yd = yd;
vars.zd = zd;
vars.angles = angles;
vars.nProjections = nProjections;
vars.dPitchX = dPitchX;
vars.dPitchY = dPitchY;
vars.r = d;
vars.sourceToCRot = r;

end