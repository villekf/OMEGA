function saveInterfile(filename, img, reko, varargin)
%SAVEINTERFILE Saves the input matrix as an interfile format image
%   This file saves the input 3D or 4D image into interfile format (32-bit
%   float). This code is based on the code from MathWorks file exchange by
%   Josh Schaefferkoetter: 
%   https://se.mathworks.com/matlabcentral/fileexchange/53745-medical-image-reader-and-viewer
%
%   If you wish to convert the image into a type other than 32-bit float,
%   change "number format:=float" and "!number of bytes per pixel:=4"
%   accordingly (e.g. number format:=signed integer, !number of bytes per
%   pixel:=4) and change img = single(img); to correspond to the correct 
%   type (e.g. img = int32(img);).
%
% Example:
%   saveInterfile(filename, img, reko)
%
% Input:
%   filename = Name of the image and header files (without file prefix)
%
%   img = The 3D or 4D image
%
%   reko = Name of the current reconstruction (can be an empty array)
%
%   options/image_properties = Either the options struct created by the
%   main files or the image_properties struct saved in the cell-matrix
%   (optional). Necessary if many of the optional values are to be saved
%   (voxel size, number of time steps, total duration, etc.).

if nargin > 3
    prop = varargin{1};
else
    prop = [];
end
img = single(img);
koko = size(img);
fid = fopen([filename '.img'],'w');
fwrite(fid, img, 'single');
fclose(fid);
hdrFile = [filename '.hdr'];
maxmin = [max(img(:)) min(img(:))];

if nargin > 3
    writeInterfileHeader(hdrFile, koko, maxmin, reko, prop);
else
    writeInterfileHeader(hdrFile, koko, maxmin, reko);
end
end

function writeInterfileHeader(hdrFile, koko, maxmin, rekot, varargin)
if nargin >= 5 && ~isempty(varargin{1})
    prop = varargin{1};
else
    prop = [];
end
fid = fopen(hdrFile,'w');

%write header file information
fprintf(fid,'!INTERFILE:=\r\n');
fprintf(fid,'%%comment:=\r\n');
fprintf(fid,['!originating system:=OMEGA\r\n']);
fprintf(fid,'%%SMS-MI header name space:=\r\n');
fprintf(fid,'%%SMS-MI version number:=\r\n');
fprintf(fid,'\r\n');
fprintf(fid,'!GENERAL DATA:=\r\n');
fprintf(fid,'%%sinogram header file:=\r\n');
fprintf(fid,'%%sinogram data file:=\r\n');
fprintf(fid,['!name of data file:=' [hdrFile(1:end-2) 'img'] '\r\n']);
fprintf(fid,'\r\n');
fprintf(fid,'!GENERAL IMAGE DATA:=\r\n');
fprintf(fid,['%%study date (yyyy:mm:dd):=' datestr(datetime('now', 'Format', 'yyyy:M:dd'),'yyyy:mm:dd') '\r\n']);
fprintf(fid,['%%study time (hh:mm:ss GMT+00:00):=' datestr(datetime('now', 'Format', 'hh:mm:ss'),'HH:MM:SS') '\r\n']);
fprintf(fid,'isotope name:=\r\n');
fprintf(fid,'isotope gamma halflife (sec):=\r\n');
fprintf(fid,'isotope branching factor:=\r\n');
fprintf(fid,['radiopharmaceutical:=\r\n']);
fprintf(fid,['%%tracer injection date (yyyy:mm:dd):=\r\n']);
fprintf(fid,['%%tracer injection time (hh:mm:ss GMT+00:00):=\r\n']);
fprintf(fid,['tracer activity at time of injection (Bq):=\r\n']);
fprintf(fid,'relative time of tracer injection (sec):=\r\n');
fprintf(fid,'injected volume (ml):=\r\n');
fprintf(fid,'image data byte order:=LITTLEENDIAN\r\n');
fprintf(fid,['%%patient orientation:=\r\n']);
fprintf(fid,['%%image orientation:=\r\n']);
fprintf(fid,'!PET data type:=image\r\n');
fprintf(fid,'number format:=float\r\n');
fprintf(fid,'!number of bytes per pixel:=4\r\n');
if length(koko) > 3 && koko(4) > 1
    fprintf(fid,'number of dimensions:=4\r\n');
else
    fprintf(fid,'number of dimensions:=3\r\n');
end
fprintf(fid,'matrix axis label[1]:=x\r\n');
fprintf(fid,'matrix axis label[2]:=y\r\n');
fprintf(fid,'matrix axis label[3]:=z\r\n');
if length(koko) > 3 && koko(4) > 1
    fprintf(fid,'matrix axis label[4]:=t\r\n');
end
fprintf(fid,['matrix size[1]:=' num2str(koko(1)) '\r\n']);
fprintf(fid,['matrix size[2]:=' num2str(koko(2)) '\r\n']);
fprintf(fid,['matrix size[3]:=' num2str(koko(3)) '\r\n']);
if length(koko) > 3 && koko(4) > 1
    fprintf(fid,['matrix size[4]:=' num2str(koko(4)) '\r\n']);
end
if nargin >= 5
    if isfield(prop,'FOV_x')
        fprintf(fid,['scale factor (mm/pixel) [1]:=' num2str(prop.FOV_x/koko(1)) '\r\n']);
        fprintf(fid,['scale factor (mm/pixel) [2]:=' num2str(prop.FOV_y/koko(2)) '\r\n']);
        fprintf(fid,['scale factor (mm/pixel) [3]:=' num2str(prop.axial_FOV/koko(3)) '\r\n']);
    elseif isfield(prop,'FOVa_x')
        fprintf(fid,['scale factor (mm/pixel) [1]:=' num2str(prop.FOVa_x/koko(1)) '\r\n']);
        fprintf(fid,['scale factor (mm/pixel) [2]:=' num2str(prop.FOVa_y/koko(2)) '\r\n']);
        fprintf(fid,['scale factor (mm/pixel) [3]:=' num2str(prop.axial_fov/koko(3)) '\r\n']);
    else
        fprintf(fid,['scale factor (mm/pixel) [1]:=0\r\n']);
        fprintf(fid,['scale factor (mm/pixel) [2]:=0\r\n']);
        fprintf(fid,['scale factor (mm/pixel) [3]:=0\r\n']);
    end
else
    fprintf(fid,['scale factor (mm/pixel) [1]:=0\r\n']);
    fprintf(fid,['scale factor (mm/pixel) [2]:=0\r\n']);
    fprintf(fid,['scale factor (mm/pixel) [3]:=0\r\n']);
end
fprintf(fid,'horizontal bed translation:=stepped\r\n');
fprintf(fid,['start horizontal bed position (mm):=\r\n']);
fprintf(fid,['end horizontal bed position (mm):=\r\n']);
fprintf(fid,['start vertical bed position (mm):=\r\n']);
fprintf(fid,'%%reconstruction diameter (mm):=\r\n');
fprintf(fid,['quantification units:=\r\n']);
fprintf(fid,['%%scanner quantification factor (Bq*s/ECAT counts):=\r\n']);
fprintf(fid,'%%decay correction:=\r\n');
fprintf(fid,['%%decay correction reference date (yyyy:mm:dd):=\r\n']);
fprintf(fid,['%%decay correction reference time (hh:mm:ss GMT+00:00):=\r\n']);
fprintf(fid,'slice orientation:=\r\n');
if isfield(prop, 'simple')
    fprintf(fid,'method of reconstruction:=OSEM\r\n');
else
    fprintf(fid,['method of reconstruction:=' rekot '\r\n']);
end
fprintf(fid,'%%PSF axial sigma (mm):=\r\n');
fprintf(fid,'%%PSF axial cutoff:=\r\n');
fprintf(fid,'%%number of radial PSF coefficients:=\r\n');
fprintf(fid,'%%PSF radial left coefficient (bins) [1]:=\r\n');
fprintf(fid,'%%PSF radial left coefficient (1/mm) [2]:=\r\n');
fprintf(fid,'%%PSF radial left coefficient (1/mm^2) [3]:=\r\n');
fprintf(fid,'%%PSF radial left coefficient (1/mm^3) [4]:=\r\n');
fprintf(fid,'%%PSF radial left coefficient (1/mm^4) [5]:=\r\n');
fprintf(fid,'%%PSF radial right coefficient (bins) [1]:=\r\n');
fprintf(fid,'%%PSF radial right coefficient (1/mm) [2]:=\r\n');
fprintf(fid,'%%PSF radial right coefficient (1/mm^2) [3]:=\r\n');
fprintf(fid,'%%PSF radial right coefficient (1/mm^3) [4]:=\r\n');
fprintf(fid,'%%PSF radial right coefficient (1/mm^4) [5]:=\r\n');
fprintf(fid,'%%PSF radial cutoff:=\r\n');
fprintf(fid,'%%gantry offset (mm) [1]:=0\r\n');
fprintf(fid,'%%gantry offset (mm) [2]:=0\r\n');
fprintf(fid,'%%gantry offset (mm) [3]:=0\r\n');
fprintf(fid,'%%gantry offset pitch (degrees):=\r\n');
fprintf(fid,'%%gantry offset yaw (degrees):=\r\n');
fprintf(fid,'%%gantry offset roll (degrees):=\r\n');
if isfield(prop,'Niter')
    fprintf(fid,['number of iterations:=' num2str(prop.Niter) '\r\n']);
else
    fprintf(fid,'number of iterations:=\r\n');
end
if isfield(prop,'subsets')
    fprintf(fid,['%%number of subsets:=' num2str(prop.subsets) '\r\n']);
else
    fprintf(fid,'%%number of subsets:=\r\n');
end
fprintf(fid,'filter name:=\r\n');
fprintf(fid,'%%xy-filter (mm):=\r\n');
fprintf(fid,'%%z-filter (mm):=\r\n');
fprintf(fid,'%%filter order:=\r\n');
fprintf(fid,'%%image zoom:=\r\n');
fprintf(fid,'%%x-offset (mm):=\r\n');
fprintf(fid,'%%y-offset (mm):=\r\n');
corrections = [];
if isfield(prop,'normalization') && prop.normalization == 1
    corrections = [corrections 'normalization,'];
elseif isfield(prop,'normalization_correction') && prop.normalization_correction == 1
    corrections = [corrections 'normalization,'];
end
if isfield(prop,'randoms') && prop.randoms == 1
    corrections = [corrections 'randoms,'];
elseif isfield(prop,'randoms_correction') && prop.randoms_correction == 1
    corrections = [corrections 'randoms,'];
end
if isfield(prop,'attenuation') && prop.attenuation == 1
    corrections = [corrections 'attenuation,'];
elseif isfield(prop,'attenuation_correction') && prop.attenuation_correction == 1
    corrections = [corrections 'attenuation,'];
end
if isfield(prop,'scatter') && prop.scatter == 1
    corrections = [corrections 'scatter'];
elseif isfield(prop,'scatter_correction') && prop.scatter_correction == 1
    corrections = [corrections 'scatter'];
end
fprintf(fid,['applied corrections:=' corrections '\r\n']);
fprintf(fid,'method of attenuation correction:=image\r\n');
fprintf(fid,'%%CT coverage:=\r\n');
fprintf(fid,'method of scatter correction:=\r\n');
fprintf(fid,'%%method of random correction:=delayed coincidences\r\n');
fprintf(fid,'%%TOF mashing factor:=\r\n');
fprintf(fid,'number of energy windows:=\r\n');
fprintf(fid,'%%energy window lower level (keV) [1]:=\r\n');
fprintf(fid,'%%energy window upper level (keV) [1]:=\r\n');
fprintf(fid,'%%coincidence window width (ns):=\r\n');
fprintf(fid,'\r\n');
fprintf(fid,'!IMAGE DATA DESCRIPTION:=\r\n');
if isfield(prop,'n_time_steps')
    fprintf(fid,['!total number of data sets:=' num2str(prop.n_time_steps) '\r\n']);
    fprintf(fid,['!image duration (sec):=' num2str(prop.total_time/prop.n_time_steps) '\r\n']);
elseif isfield(prop,'partitions')
    fprintf(fid,['!total number of data sets:=' num2str(prop.partitions) '\r\n']);
    fprintf(fid,['!image duration (sec):=' num2str(prop.tot_time/prop.partitions) '\r\n']);
else
    fprintf(fid,'!total number of data sets:=1\r\n');
    fprintf(fid,'!image duration (sec):=0\r\n');
end
if isfield(prop,'start')
    fprintf(fid,['!image relative start time (sec):=' num2str(prop.start) '\r\n']);
else
    fprintf(fid,'!image relative start time (sec):=0\r\n');
end
fprintf(fid,['total prompts:=0\r\n']);
fprintf(fid,['%%total randoms:=0\r\n']);
fprintf(fid,['%%total net trues:=0\r\n']);
fprintf(fid,'%%GIM loss fraction:=\r\n');
fprintf(fid,'%%PDR loss fraction:=\r\n');
fprintf(fid,'%%total uncorrected singles rate:=\r\n');
fprintf(fid,'%%image slope:=\r\n');
fprintf(fid,'%%image intercept:=\r\n');
fprintf(fid,['maximum pixel count:=' num2str(maxmin(1)) '\r\n']);
fprintf(fid,['minimum pixel count:=' num2str(maxmin(2)) '\r\n']);
fprintf(fid,'\r\n');
fprintf(fid,'%%SUPPLEMENTARY ATTRIBUTES:=\r\n');
fclose(fid);

end