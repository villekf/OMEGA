function options = loadPlanmecaData(options)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[file, fpath] = uigetfile('*.*','Select Planmeca metadata file');
if isequal(file, 0)
    error('No file was selected')
end
nimi = [fpath file];
fid = fopen(nimi);
hdr = textscan(fid,'%s','Delimiter','=');
hdr = hdr{1};
fclose(fid);
options.nProjections = str2double(hdr{find(strcmp(hdr,'Nframes')) + 1});
magnification = str2double(hdr{find(strcmp(hdr,'magnification')) + 1});
file = 'geometryFile';
nimi = [fpath file];
fid = fopen(nimi);
if fid == -1
    error('Planmeca geometryFile not found. Make sure the geometryFile is in the same folder as the metadata.')
end
geom = textscan(fid,'%f','Delimiter',' ');
geom = geom{1};
fclose(fid);
geom = reshape(geom, 9, options.nProjections)';
% file = 'geometryFile_calm';
% nimi = [fpath file];
% fid = fopen(nimi);
% geom2 = textscan(fid,'%f','Delimiter',' ');
% geom2 = geom2{1};
% fclose(fid);
% geom = reshape(geom2, 9, options.nProjections)';
% geom = [geom2(:,1:5),geom(:,6),geom2(:,7:end)];
options.angles = ((geom(:,end))) / 180 * pi;
% options.pitchRoll = ([(geom(:,end-2)) (geom(:,end-1))]);
% options.pitchRoll = zeros(size(options.pitchRoll),'single');
options.sourceToCRot = options.sourceToDetector / (magnification / 100);
% file = 'geometryFileVec';
% nimi = [fpath file];
% fid = fopen(nimi);
% u = textscan(fid,'%f','Delimiter',' ');
% u = u{1};
% fclose(fid);
% u = reshape(u, 12, options.nProjections)';
% options.uV = ([u(:,8), u(:,7), u(:,9), u(:,11), u(:,10), u(:,12)])';
options.x = ([geom(:,1), geom(:,4)]);
% options.x(:,1) = options.x(:,1) + options.dPitch * options.ySize / 2 * cosd(options.angles(:,1));
options.y = ([geom(:,2), geom(:,5)]);
% options.y(:,1) = options.y(:,1) + options.dPitch * options.ySize / 2 * sind(options.angles(:,1));
options.z = ([geom(:,3), geom(:,6)]);
% options.z(:,1) = options.z(:,1) - options.dPitch * options.xSize / 2;
options.pitchRoll = [geom(:,7) geom(:,8)] / 180 * pi;
fid = fopen([fpath 'flat.raw']);
flat = fread(fid, [options.ySize * options.binning, options.xSize * options.binning], 'uint16=>double') * options.flatFieldScaling;
fclose(fid);
options.SinM = zeros(options.ySize, options.xSize,options.nProjections,'single');
fpath = [fpath 'corrected/'];
for kk = 1 : options.nProjections
    fid = fopen([fpath 'frame-' num2str(kk + 99) '.raw']);
    apu = fread(fid, [options.ySize * options.binning, options.xSize * options.binning], 'uint16=>double');
    options.SinM(:,:,kk) = single(imresize(apu, 1/options.binning) ./ flat(1));
    fclose(fid);
end
% options.SinM = circshift(flip(options.SinM,3),0,3);
end