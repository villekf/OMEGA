function saveInterfile(filename, img, varargin)
%SAVEINTERFILE Saves the input matrix as an interfile format image
%  This file saves the input 3D or 4D image into interfile format. By
%  default the header file is always saved in UNIX text format. You can
%  (optionally) change this into Windows format, regardless of your current
%  OS. This code is based on the code from MathWorks file exchange by Josh
%  Schaefferkoetter:
%  https://se.mathworks.com/matlabcentral/fileexchange/53745-medical-image-reader-and-viewer
%
% Examples:
%  saveInterfile(filename, img)
%  saveInterfile('image_name', img)
%  saveInterfile(filename, img, reko)
%  saveInterfile(filename, img, reko, type)
%  saveInterfile(filename, img, [], [], options)
%  saveInterfile(filename, img, [], [], pixelSize)
%  saveInterfile(filename, img, reko, type, options)
%  saveInterfile(filename, img, reko, type, options, windows_format)
%  saveInterfile(filename, img, [], [], pixelSize, windows_format)
%
% Input:
%  filename = Name of the image and header files (without file type).
%  Saves the header file as [filename '.h33'] and image file as [filename
%  '.i33']. Needs to be char input.
%
%  img = The 3D or 4D image. Note that vector format images are not
%  supported.
%
%  reko = Name of the current reconstruction (can be an empty array). For
%  naming purposes only. If omitted, will use an empty array (no name).
%  This is optional variable.
%
%  type = Data type, e.g. 'single', 'int8', 'uint32', etc. Default is
%  the same type as the input image. If omitted, will use the default
%  value.
%
%  options/image_properties = Either the options struct created by the
%  main files or the image_properties struct saved in the cell-matrix
%  (optional). Necessary if any of the optional values are to be saved
%  (voxel size, number of time steps, total duration, etc.).
%
%  pixelSize = If saving data where the above struct is not available, this
%  can be used to input the pixel sizes (in mm) for x, y and (possibly) z.
%  This should be input as a vector with the same number of dimensions as
%  the image. Use either the struct or the pixel size, not both. This
%  variable is optional, if omitted values of 1 mm are used by default.
%
%  windows_format = Text file format. Default is UNIX style (LF), but
%  setting this value to true uses Windows file format (CR LF). Omitting
%  will use UNIX style.
%
% See also saveImage, saveMetaImage

if nargin > 2 && ~isempty(varargin{1})
    reko = varargin{1};
else
    reko = [];
end
if nargin > 3 && ~isempty(varargin{2})
    type = varargin{2};
else
    type = class(img);
end
if nargin > 4 && ~isempty(varargin{3})
    prop = varargin{3};
else
    prop = [];
end
if nargin > 5 && ~isempty(varargin{4})
    windows_format = varargin{4};
else
    windows_format = false;
end
img = cast(img,type);
koko = size(img);
if length(koko) == 2
    koko = [koko, 1];
end
if ~isstruct(prop)
    if length(prop) == 2
        prop = [prop, 1];
    end
end
fid = fopen([filename '.i33'],'w');
fwrite(fid, img(:), type, 0, 'l');
fclose(fid);
hdrFile = [filename '.h33'];
if isfloat(img)
    imgFinite = img(isfinite(img));
else
    imgFinite = img(:);
end
if isempty(imgFinite)
    maxmin = [0 0];
else
    maxmin = [max(imgFinite) min(imgFinite)];
end

writeInterfileHeader(hdrFile, koko, maxmin, reko, type, windows_format, prop);
end

function writeInterfileHeader(hdrFile, koko, maxmin, rekot, type, windows_format, prop)
% The data file is referred to by its name only (no path)
[~, dataName] = fileparts(hdrFile);
dataFile = [dataName '.i33'];
if isempty(rekot)
    rekot = '';
end
fid = fopen(hdrFile,'w');

%write header file information
fprintf(fid,'!INTERFILE := \r\n');
fprintf(fid,'%%comment := \r\n');
fprintf(fid,['!originating system := OMEGA\r\n']);
fprintf(fid,'%%SMS-MI header name space := \r\n');
fprintf(fid,'%%SMS-MI version number := \r\n');
fprintf(fid,'\r\n');
fprintf(fid,'!GENERAL DATA := \r\n');
fprintf(fid,'%%sinogram header file := \r\n');
fprintf(fid,'%%sinogram data file := \r\n');
fprintf(fid,'!name of data file := %s\r\n', dataFile);
fprintf(fid,'\r\n');
fprintf(fid,'!GENERAL IMAGE DATA := \r\n');
% Number of time steps: nTimeSteps (image_properties), n_time_steps (old image_properties) or partitions (options)
n_time_steps = numTimeSteps(prop);
if ~isempty(n_time_steps) && n_time_steps > 1
    fprintf(fid,'!type of data := Dynamic\r\n');
    fprintf(fid,['!total number of images := ' num2str(n_time_steps * koko(3)) '\r\n']);
else
    fprintf(fid,'!type of data := Static\r\n');
    fprintf(fid,['!total number of images := ' num2str(koko(3)) '\r\n']);
end
fprintf(fid,['!number of slices := ' num2str(koko(3)) '\r\n']);
fprintf(fid,['%%study date (yyyy:mm:dd) := ' datestr(datetime('now', 'Format', 'yyyy:M:dd'),'yyyy:mm:dd') '\r\n']);
fprintf(fid,['%%study time (hh:mm:ss GMT+00:00) := ' datestr(datetime('now', 'Format', 'hh:mm:ss'),'HH:MM:SS') '\r\n']);
fprintf(fid,'isotope name := \r\n');
fprintf(fid,'isotope gamma halflife (sec) := \r\n');
fprintf(fid,'isotope branching factor := \r\n');
fprintf(fid,['radiopharmaceutical := \r\n']);
fprintf(fid,['%%tracer injection date (yyyy:mm:dd) := \r\n']);
fprintf(fid,['%%tracer injection time (hh:mm:ss GMT+00:00) := \r\n']);
fprintf(fid,['tracer activity at time of injection (Bq) := \r\n']);
fprintf(fid,'relative time of tracer injection (sec) := \r\n');
fprintf(fid,'injected volume (ml) := \r\n');
fprintf(fid,'!imagedata byte order := LITTLEENDIAN\r\n');
fprintf(fid,['%%patient orientation := \r\n']);
fprintf(fid,['%%image orientation := \r\n']);
fprintf(fid,'!PET data type := image\r\n');
if strcmp(type, 'single') || strcmp(type, 'double')
    if strcmp(type, 'single')
        fprintf(fid,'!number of bytes per pixel := 4\r\n');
        fprintf(fid,'!number format := float\r\n');
    else
        fprintf(fid,'!number of bytes per pixel := 8\r\n');
        fprintf(fid,'!number format := double\r\n');
    end
elseif strcmp(type, 'int8') || strcmp(type, 'int16') || strcmp(type, 'int32') || strcmp(type, 'int64') || strcmp(type, 'char')
    fprintf(fid,'!number format := signed integer\r\n');
    if strcmp(type, 'int8') || strcmp(type, 'char')
        fprintf(fid,'!number of bytes per pixel := 1\r\n');
    elseif strcmp(type, 'int16')
        fprintf(fid,'!number of bytes per pixel := 2\r\n');
    elseif strcmp(type, 'int32')
        fprintf(fid,'!number of bytes per pixel := 4\r\n');
    else
        fprintf(fid,'!number of bytes per pixel := 8\r\n');
    end
elseif strcmp(type, 'uint8') || strcmp(type, 'uint16') || strcmp(type, 'uint32') || strcmp(type, 'uint64') || strcmp(type, 'uchar')
    fprintf(fid,'!number format := unsigned integer\r\n');
    if strcmp(type, 'uint8') || strcmp(type, 'uchar')
        fprintf(fid,'!number of bytes per pixel := 1\r\n');
    elseif strcmp(type, 'uint16')
        fprintf(fid,'!number of bytes per pixel := 2\r\n');
    elseif strcmp(type, 'uint32')
        fprintf(fid,'!number of bytes per pixel := 4\r\n');
    else
        fprintf(fid,'!number of bytes per pixel := 8\r\n');
    end
end
if length(koko) > 3 && koko(4) > 1
    fprintf(fid,'number of dimensions := 4\r\n');
else
    fprintf(fid,'number of dimensions := 3\r\n');
end
fprintf(fid,'matrix axis label[1] := x\r\n');
fprintf(fid,'matrix axis label[2] := y\r\n');
fprintf(fid,'matrix axis label[3] := z\r\n');
if length(koko) > 3 && koko(4) > 1
    fprintf(fid,'matrix axis label[4] := t\r\n');
end
fprintf(fid,['!matrix size [1] := ' num2str(koko(1)) '\r\n']);
fprintf(fid,['!matrix size [2] := ' num2str(koko(2)) '\r\n']);
fprintf(fid,['!matrix size [3] := ' num2str(koko(3)) '\r\n']);
if length(koko) > 3 && koko(4) > 1
    fprintf(fid,['!matrix size [4] := ' num2str(koko(4)) '\r\n']);
end
% Voxel size (mm): either from the field of view of the options/image_properties or from the pixel size vector
voxel = [];
if isstruct(prop)
    % Field of view: FOVa_x/FOVa_y/axial_fov (options), FOVx/FOVy/axialFOV (image_properties) or FOV_x/FOV_y/axial_FOV (old image_properties)
    FOV_x = firstValue(getField(prop, {'FOVa_x', 'FOVx', 'FOV_x'}));
    FOV_y = firstValue(getField(prop, {'FOVa_y', 'FOVy', 'FOV_y'}));
    axial_FOV = firstValue(getField(prop, {'axial_fov', 'axialFOV', 'axial_FOV'}));
    if ~isempty(FOV_x) && ~isempty(FOV_y) && ~isempty(axial_FOV) && FOV_x > 0 && FOV_y > 0 && axial_FOV > 0
        voxel = [FOV_x / koko(1), FOV_y / koko(2), axial_FOV / koko(3)];
    end
elseif isnumeric(prop) && numel(prop) >= 3
    voxel = double(prop(1:3));
end
if ~isempty(voxel)
    fprintf(fid,['scaling factor (mm/pixel) [1] := +' num2str(voxel(1), '%e') '\r\n']);
    fprintf(fid,['scaling factor (mm/pixel) [2] := +' num2str(voxel(2), '%e') '\r\n']);
    fprintf(fid,['scaling factor (mm/pixel) [3] := +' num2str(voxel(3), '%e') '\r\n']);
    fprintf(fid,['slice thickness (pixels) := +' num2str(voxel(3), '%e') '\r\n']);
else
    fprintf(fid,'scaling factor (mm/pixel) [1] := 1\r\n');
    fprintf(fid,'scaling factor (mm/pixel) [2] := 1\r\n');
    fprintf(fid,'scaling factor (mm/pixel) [3] := 1\r\n');
    fprintf(fid,'slice thickness (pixels) := 1\r\n');
end
fprintf(fid,'horizontal bed translation := stepped\r\n');
fprintf(fid,['start horizontal bed position (mm) := \r\n']);
fprintf(fid,['end horizontal bed position (mm) := \r\n']);
fprintf(fid,['start vertical bed position (mm) := \r\n']);
fprintf(fid,'%%reconstruction diameter (mm) := \r\n');
fprintf(fid,['quantification units := \r\n']);
fprintf(fid,['%%scanner quantification factor (Bq*s/ECAT counts) := \r\n']);
fprintf(fid,'%%decay correction := \r\n');
fprintf(fid,['%%decay correction reference date (yyyy:mm:dd) := \r\n']);
fprintf(fid,['%%decay correction reference time (hh:mm:ss GMT+00:00) := \r\n']);
fprintf(fid,'slice orientation := \r\n');
if isfield(prop, 'simple')
    fprintf(fid,'method of reconstruction := OSEM\r\n');
else
    fprintf(fid,'method of reconstruction := %s\r\n', rekot);
end
fprintf(fid,'%%PSF axial sigma (mm) := \r\n');
fprintf(fid,'%%PSF axial cutoff := \r\n');
fprintf(fid,'%%number of radial PSF coefficients := \r\n');
fprintf(fid,'%%PSF radial left coefficient (bins) [1] := \r\n');
fprintf(fid,'%%PSF radial left coefficient (1/mm) [2] := \r\n');
fprintf(fid,'%%PSF radial left coefficient (1/mm^2) [3] := \r\n');
fprintf(fid,'%%PSF radial left coefficient (1/mm^3) [4] := \r\n');
fprintf(fid,'%%PSF radial left coefficient (1/mm^4) [5] := \r\n');
fprintf(fid,'%%PSF radial right coefficient (bins) [1] := \r\n');
fprintf(fid,'%%PSF radial right coefficient (1/mm) [2] := \r\n');
fprintf(fid,'%%PSF radial right coefficient (1/mm^2) [3] := \r\n');
fprintf(fid,'%%PSF radial right coefficient (1/mm^3) [4] := \r\n');
fprintf(fid,'%%PSF radial right coefficient (1/mm^4) [5] := \r\n');
fprintf(fid,'%%PSF radial cutoff := \r\n');
fprintf(fid,'%%gantry offset (mm) [1] := 0\r\n');
fprintf(fid,'%%gantry offset (mm) [2] := 0\r\n');
fprintf(fid,'%%gantry offset (mm) [3] := 0\r\n');
fprintf(fid,'%%gantry offset pitch (degrees) := \r\n');
fprintf(fid,'%%gantry offset yaw (degrees) := \r\n');
fprintf(fid,'%%gantry offset roll (degrees) := \r\n');
% Niter/subsets (options) or nIter/nSubsets (image_properties)
Niter = firstValue(getField(prop, {'Niter', 'nIter'}));
if ~isempty(Niter)
    fprintf(fid,['number of iterations := ' num2str(Niter) '\r\n']);
else
    fprintf(fid,'number of iterations := \r\n');
end
subsets = firstValue(getField(prop, {'subsets', 'nSubsets'}));
if ~isempty(subsets)
    fprintf(fid,['%%number of subsets := ' num2str(subsets) '\r\n']);
else
    fprintf(fid,'%%number of subsets := \r\n');
end
fprintf(fid,'filter name := \r\n');
fprintf(fid,'%%xy-filter (mm) := \r\n');
fprintf(fid,'%%z-filter (mm) := \r\n');
fprintf(fid,'%%filter order := \r\n');
fprintf(fid,'%%image zoom := \r\n');
fprintf(fid,'%%x-offset (mm) := \r\n');
fprintf(fid,'%%y-offset (mm) := \r\n');
corrections = [];
% The correction flags are normalizationCorrection etc. in image_properties and normalization_correction etc. in options
if isOne(getField(prop, {'normalization', 'normalization_correction', 'normalizationCorrection'}))
    corrections = [corrections 'normalization,'];
end
if isOne(getField(prop, {'randoms', 'randoms_correction', 'randomsCorrection'}))
    corrections = [corrections 'randoms,'];
end
if isOne(getField(prop, {'attenuation', 'attenuation_correction', 'attenuationCorrection'}))
    corrections = [corrections 'attenuation,'];
end
if isOne(getField(prop, {'scatter', 'scatter_correction', 'scatterCorrection'}))
    corrections = [corrections 'scatter'];
end
fprintf(fid,'applied corrections := %s\r\n', corrections);
fprintf(fid,'method of attenuation correction := image\r\n');
fprintf(fid,'%%CT coverage := \r\n');
fprintf(fid,'method of scatter correction := \r\n');
fprintf(fid,'%%method of random correction := delayed coincidences\r\n');
fprintf(fid,'%%TOF mashing factor := \r\n');
fprintf(fid,'number of energy windows := \r\n');
fprintf(fid,'%%energy window lower level (keV) [1] := \r\n');
fprintf(fid,'%%energy window upper level (keV) [1] := \r\n');
fprintf(fid,'%%coincidence window width (ns) := \r\n');
fprintf(fid,'\r\n');
fprintf(fid,'!IMAGE DATA DESCRIPTION := \r\n');
% Total time: totalTime (image_properties), total_time (old image_properties) or tot_time (options)
total_time = firstValue(getField(prop, {'totalTime', 'total_time', 'tot_time'}));
if ~isempty(n_time_steps)
    fprintf(fid,['!total number of data sets := ' num2str(n_time_steps) '\r\n']);
    if ~isempty(total_time) && isfinite(total_time)
        fprintf(fid,['!image duration (sec) := ' num2str(total_time/n_time_steps) '\r\n']);
    else
        fprintf(fid,'!image duration (sec) := 0\r\n');
    end
else
    fprintf(fid,'!total number of data sets := 1\r\n');
    fprintf(fid,'!image duration (sec) := 0\r\n');
end
% Start time: start (options) or startTime (image_properties)
start_time = firstValue(getField(prop, {'start', 'startTime'}));
if ~isempty(start_time)
    fprintf(fid,['!image relative start time (sec) := ' num2str(start_time) '\r\n']);
else
    fprintf(fid,'!image relative start time (sec) := 0\r\n');
end
fprintf(fid,['total prompts := 0\r\n']);
fprintf(fid,['%%total randoms := 0\r\n']);
fprintf(fid,['%%total net trues := 0\r\n']);
fprintf(fid,'%%GIM loss fraction := \r\n');
fprintf(fid,'%%PDR loss fraction := \r\n');
fprintf(fid,'%%total uncorrected singles rate := \r\n');
fprintf(fid,'%%image slope := \r\n');
fprintf(fid,'%%image intercept := \r\n');
if maxmin(1) >= 0
    fprintf(fid,['maximum pixel count := +' num2str(maxmin(1), '%e') '\r\n']);
else
    fprintf(fid,['maximum pixel count := ' num2str(maxmin(1), '%e') '\r\n']);
end
if maxmin(2) >= 0
    fprintf(fid,['minimum pixel count := +' num2str(maxmin(2), '%e') '\r\n']);
else
    fprintf(fid,['minimum pixel count := ' num2str(maxmin(2), '%e') '\r\n']);
end
fprintf(fid,'\r\n');
fprintf(fid,'%%SUPPLEMENTARY ATTRIBUTES := \r\n');
scanner = scannerParameters(prop);
for kk = 1 : numel(scanner)
    fprintf(fid,'%s\r\n', scanner{kk});
end
fclose(fid);

% The header was written with CR LF line endings. Convert to the requested format (LF by default)
fid = fopen(hdrFile,'r');
fcontent = fread(fid,'uint8');
fclose(fid);
fcontent(fcontent == 13) = [];
if windows_format
    fcontent = uint8(strrep(char(fcontent(:)'), newline, [char(13) newline]))';
end
fid = fopen(hdrFile,'w');
fwrite(fid,fcontent,'uint8');
fclose(fid);

end

function lines = scannerParameters(prop)
%SCANNERPARAMETERS Returns the lines of the scanner parameters block of the header (empty if the scanner is not described by
%the input struct). Accepts both the options struct and the image_properties struct.
lines = {};
if ~isstruct(prop) || isOne(getField(prop, {'CT'})) || isOne(getField(prop, {'SPECT'}))
    return
end
if isfield(prop, 'rings')
    % Options struct
    nRings = firstValue(prop.rings);
    perRing = firstValue(getField(prop, {'det_w_pseudo'}));
    diameter = firstValue(getField(prop, {'diameter'}));
    axialBlocks = firstValue(getField(prop, {'linear_multip'}));
    cryst = firstValue(getField(prop, {'cryst_per_block'}));
    crystAxial = firstValue(getField(prop, {'cryst_per_block_axial'}));
    angles = firstValue(getField(prop, {'Nang'}));
    blocks = firstValue(getField(prop, {'blocks_per_ring'}));
elseif isfield(prop, 'nRings') || isfield(prop, 'machine_rings')
    % image_properties struct (nRings, detectorsPerRing, ...) or old image_properties struct (machine_rings, detectors_per_ring, ...)
    nRings = firstValue(getField(prop, {'nRings', 'machine_rings'}));
    perRing = firstValue(getField(prop, {'detectorsPerRing', 'detectors_per_ring'}));
    diameter = firstValue(getField(prop, {'boreDiameter', 'bore_diameter'}));
    axialBlocks = firstValue(getField(prop, {'axialBlocks', 'axial_blocks'}));
    cryst = firstValue(getField(prop, {'crystalsPerBlock', 'crystals_per_block'}));
    crystAxial = firstValue(getField(prop, {'crystalsPerBlockAxial'}));
    angles = firstValue(getField(prop, {'nSinogramAngles', 'Nang'}));
    blocks = firstValue(getField(prop, {'blocksPerRing', 'blocks_per_ring'}));
else
    return
end
if isempty(crystAxial)
    crystAxial = cryst;
end
if isempty(nRings) || isempty(perRing) || isempty(diameter) || isempty(axialBlocks) || isempty(cryst) || isempty(angles) || isempty(blocks) ...
        || nRings <= 0 || perRing <= 0
    return
end
DOI = firstValue(getField(prop, {'DOI'}));
name = getField(prop, {'scannerName', 'machine_name'});
if isstring(name) && isscalar(name)
    name = char(name);
elseif ~ischar(name)
    name = '';
end
lines{end + 1} = 'Scanner parameters := ';
lines{end + 1} = ['Scanner type := ' name];
lines{end + 1} = ['Number of rings := ' num2str(nRings)];
lines{end + 1} = ['Number of detectors per ring := ' num2str(perRing)];
lines{end + 1} = ['Inner ring diameter (cm) := ' num2str(diameter / 10)];
if ~isempty(DOI)
    lines{end + 1} = ['Average depth of interaction (cm) := ' num2str(DOI / 10)];
else
    lines{end + 1} = 'Average depth of interaction (cm) := 0';
end
lines{end + 1} = 'View offset (degrees) := 0';
lines{end + 1} = ['Maximum number of non-arc-corrected bins := ' num2str(angles)];
lines{end + 1} = ['Default number of arc-corrected bins := ' num2str(angles)];
lines{end + 1} = ['Number of blocks per bucket in transaxial direction := ' num2str(blocks)];
lines{end + 1} = ['Number of blocks per bucket in axial direction := ' num2str(axialBlocks)];
lines{end + 1} = ['Number of crystals per block in axial direction := ' num2str(crystAxial)];
lines{end + 1} = ['Number of crystals per block in transaxial direction := ' num2str(cryst)];
lines{end + 1} = 'Number of detector layers := 1';
lines{end + 1} = 'end scanner parameters := ';
end

function val = getField(prop, names)
%GETFIELD Returns the value of the first of the fields in names that exists in the struct prop, [] if none exist
val = [];
if ~isstruct(prop)
    return
end
for kk = 1 : numel(names)
    if isfield(prop, names{kk})
        val = prop.(names{kk});
        return
    end
end
end

function val = firstValue(x)
%FIRSTVALUE Returns the first element of x as a double (empty if x is empty)
if isempty(x) || (~isnumeric(x) && ~islogical(x))
    val = [];
else
    val = double(x(1));
end
end

function n = numTimeSteps(prop)
%NUMTIMESTEPS Returns the number of time steps: nTimeSteps (image_properties), n_time_steps (old image_properties) or partitions
%(options). If the value is a vector, its number of elements is used.
n = getField(prop, {'nTimeSteps', 'n_time_steps', 'partitions'});
if isempty(n) || ~isnumeric(n)
    n = [];
elseif isscalar(n)
    n = double(n);
else
    n = numel(n);
end
end

function tf = isOne(x)
%ISONE True if x is non-empty and all its elements are equal to one
tf = ~isempty(x) && (isnumeric(x) || islogical(x)) && all(x(:) == 1);
end
