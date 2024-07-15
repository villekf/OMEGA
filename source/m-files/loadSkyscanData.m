function options = loadSkyscanData(options)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if ~isfield(options, 'fpath') || isempty(options.fpath)
    [file, fpath] = uigetfile('*.log','Select Skyscan log file');
    if isequal(file, 0)
        error('No file was selected')
    end
    nimi = [fpath file];
    [fpath,file,ext] = fileparts(nimi);
else
    if exist(options.fpath, 'file') == 0
        warning(['Input file ' options.fpath ' not found!'])
        [file, fpath] = uigetfile('*.log','Select Skyscan log file');
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
if isfield(options, 'only_reconstructions') 
    only_reconstructions = options.only_reconstructions;
else
    only_reconstructions = false;
end
ind = find(strcmp(hdr,'Number of Columns')) + 1;
options.nRowsD = str2double(hdr{ind(1)}) / options.binning;
ind = find(strcmp(hdr,'Number of Rows')) + 1;
options.nColsD = str2double(hdr{ind(1)}) / options.binning;
ind = find(strcmp(hdr,'Camera Pixel Size (um)')) + 1;
options.dPitchX = str2double(hdr{ind(1)}) * options.binning / 1000;
ind = find(strcmp(hdr,'Camera Pixel Size (um)')) + 1;
options.dPitchY = str2double(hdr{ind(1)}) * options.binning / 1000;
ind = find(strcmp(hdr,'Camera to Source (mm)')) + 1;
options.sourceToDetector = str2double(hdr{ind(1)});
ind = find(strcmp(hdr,'Object to Source (mm)')) + 1;
options.sourceToCRot = str2double(hdr{ind(1)});
ind = find(strcmp(hdr,'Number of Files')) + 1;
options.nProjections = str2double(hdr{ind(1)}) - 1;
if ~only_reconstructions || ~isfield(options,'SinM')
    fpath = [fpath filesep file '0000.tif'];
    options.SinM = loadProjectionImages(options.nProjections,options.binning, fpath);
    options.SinM = permute(options.SinM, [2 1 3]);
end
end