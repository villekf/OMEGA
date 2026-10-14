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
    r = double(typecast(info.Private_7031_1003, 'single')); % DetectorFocalCenterRadialDistance
    deltar = double(typecast(info.Private_7033_100d, 'single')); % SourceRadialDistanceShift
    rho = r + deltar; 
    a = double(typecast(info.Private_7031_1001, 'single')) - pi/2; % DetectorFocalCenterAngularPosition
    deltaphi = double(typecast(info.Private_7033_100b, 'single')); % SourceAngularPositionShift
    phi = a + deltaphi;
    deltaz = double(typecast(info.Private_7033_100c, 'single')); % SourceAxialPositionShift
    zf = double(typecast(info.Private_7031_1002, 'single')); % DetectorFocalCenterAxialPosition
    zs(kk) = zf + deltaz;
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
    d = double(typecast(info.Private_7031_1031, 'single')); % ConstantRadialDistance
    % The detector is an arc of radius d centered at the focal center
    % (-r cos(a), r sin(a)), the channels are spaced by dPitchX / d radians
    dtheta = double(dPitchX) / d;
    detElements = double(typecast(info.Private_7031_1033, 'single')); % DetectorCentralElement
    % Central channel offset from the geometric center of the channels
    deltaCh = detElements(1) - double(nRowsD)/2 - .5;
    % The central element is not the middle of the detector, so the
    % detector center is rotated along the arc
    aD = a - deltaCh * dtheta;
    % Axial offset of the central row (sign not verified, zero for the AAPM data)
    colStart = detElements(2) - double(nColsD)/2 - .5;
    xd(kk) = -r*cos(a) + d * cos(aD);
    yd(kk) = r*sin(a) - d * sin(aD);
    zd(kk) = zf + colStart * double(dPitchY);
    % The projector derives the arc center from the detector center and the
    % angle, hence the rotated angle has to be stored
    angles(kk) = aD;
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
vars.r = single(d);
vars.sourceToCRot = single(r);

end