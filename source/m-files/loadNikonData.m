function options = loadNikonData(options)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if ~isfield(options, 'fpath') || isempty(options.fpath)
    [file, fpath] = uigetfile('*.xtekct','Select Nikon xtekct file');
    if isequal(file, 0)
        error('No file was selected')
    end
    nimi = [fpath file];
    [fpath,file,ext] = fileparts(nimi);
else
    if exist(options.fpath, 'file') == 0
        warning(['Input file ' options.fpath ' not found!'])
        [file, fpath] = uigetfile('*.xtekct','Select Nikon xtekct file');
        if isequal(file, 0)
            error('No file was selected')
        end
        nimi = [fpath file];
        [fpath,file,ext] = fileparts(nimi);
    else
        [fpath,file,ext] = fileparts(options.fpath);
        nimi = [fpath filesep file ext];
    end
    % nimi = options.fpath;
    % fpath = nimi(1:end-7);
end
fid = fopen(nimi);
hdr = textscan(fid,'%s','Delimiter','=');
hdr = hdr{1};
fclose(fid);
if isfield(options, 'binning') 
    binning = options.binning;
else
    binning = 1;
end
if isfield(options, 'only_reconstructions') 
    only_reconstructions = options.only_reconstructions;
else
    only_reconstructions = false;
end
if ~isfield(options,'Nx')
    ind = find(strcmp(hdr,'VoxelsX')) + 1;
    options.Nx = str2double(hdr{ind(1)});
end
if ~isfield(options,'Ny')
    ind = find(strcmp(hdr,'VoxelsY')) + 1;
    options.Ny = str2double(hdr{ind(1)});
end
if ~isfield(options,'Nz')
    ind = find(strcmp(hdr,'VoxelsZ')) + 1;
    options.Nz = str2double(hdr{ind(1)});
end
ind = find(strcmp(hdr,'DetectorPixelsY')) + 1;
options.nRowsD = str2double(hdr{ind(1)}) / binning;
ind = find(strcmp(hdr,'DetectorPixelsX')) + 1;
options.nColsD = str2double(hdr{ind(1)}) / binning;
ind = find(strcmp(hdr,'DetectorPixelSizeY')) + 1;
options.dPitchX = str2double(hdr{ind(1)}) * binning;
ind = find(strcmp(hdr,'DetectorPixelSizeX')) + 1;
options.dPitchY = str2double(hdr{ind(1)}) * binning;
ind = find(strcmp(hdr,'SrcToDetector')) + 1;
options.sourceToDetector = str2double(hdr{ind(1)});
ind = find(strcmp(hdr,'SrcToObject')) + 1;
options.sourceToCRot = str2double(hdr{ind(1)});
ind = find(strcmp(hdr,'ObjectOffsetX')) + 1;
options.oOffsetX = str2double(hdr{ind(1)});
ind = find(strcmp(hdr,'ObjectOffsetY')) + 1;
options.oOffsetY = str2double(hdr{ind(1)});
ind = find(strcmp(hdr,'DetectorOffsetY')) + 1;
options.detOffsetRow = str2double(hdr{ind(1)});
ind = find(strcmp(hdr,'DetectorOffsetX')) + 1;
options.detOffsetCol = str2double(hdr{ind(1)});
ind = find(strcmp(hdr,'MaskRadius')) + 1;
options.MaskRadius = str2double(hdr{ind(1)});
fileID = fopen([fpath filesep file '.ang']);
data = textscan(fileID,'%s');
options.nProjections = (numel(data{1}) - 2) / 2;
fclose(fileID);
options.angles = zeros(options.nProjections,1);
data = data{1}(3:end);
for kk = 1 : options.nProjections
    options.angles(kk) = -str2double(data{kk * 2});
end
if ~only_reconstructions || ~isfield(options,'SinM')
    fpath = [fpath filesep file '_0001.tif'];
    options.SinM = loadProjectionImages(options.nProjections,options.binning, fpath);
    options.SinM = permute(options.SinM, [2 1 3]);
end
end