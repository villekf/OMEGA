function options = loadElektaXVIData(options)
%LOADELEKTAXVIDATA Loads Elekta XVI CBCT projection data (.his + _Frames.xml)
%   Loads the raw projection images (Heimann HIS files, one per projection)
%   and the geometry information (_Frames.xml) of an Elekta XVI CBCT scan,
%   e.g. an XVI "projection export" folder. The projections are converted to
%   intensities (the XVI raw values are inverted, i.e. air has the highest
%   value after the conversion) and the OMEGA standard CT geometry fields are
%   set (angles, sourceToCRot, sourceToDetector, detector pitch, per
%   projection detector offsets), so no custom source/detector coordinates are
%   needed. Only frames with Exposed = True are used.
%
%   Inputs (options fields):
%       fpath = The folder that contains the .his files and the
%       *_Frames.xml file. If omitted, or not found, a folder selection
%       dialog is shown.
%       binning = (Optional, default 1) Binning factor. The projections are
%       binned in both dimensions by averaging (mean) binning.
%       sourceToCRot = (Optional) Source to center of rotation distance
%       (mm). Default is the nominal XVI value, 1000 mm. A value that has
%       already been set (> 1) is kept.
%       sourceToDetector = (Optional) Source to detector distance (mm).
%       Default is the nominal XVI value, 1536 mm. A value that has already
%       been set (> 1) is kept.
%       dPitchX, dPitchY = (Optional) Detector pixel pitch (mm) AFTER
%       binning, i.e. raw pitch * binning. Default is the nominal XVI raw
%       pitch 0.8 mm * binning. Values that have already been set (> 0) are
%       kept.
%       elektaUCentreScale, elektaVCentreScale = (Optional, default 1)
%       Scale factors that convert the UCentre/VCentre values (mm) of the
%       _Frames.xml file to the detector offsets. The nominal convention
%       (scale = 1) is used by default. The sign convention is fixed:
%       detOffsetRow = elektaUCentreScale * UCentre and
%       detOffsetCol = -elektaVCentreScale * VCentre (mm, in the OMEGA
%       detector frame). Other values (also fitted effective values) can be
%       used if the nominal ones do not fit the data well.
%       elektaFlatPercentile = (Optional, default 99.9) The percentile of
%       the binned intensity values (all projections) that is used as the
%       global air (flat) value, options.flat. If options.flat has already
%       been set (> 0) it is not modified.
%
%   Outputs:
%       options.SinM = The projection intensities, nRowsD x nColsD x
%       nProjections, single. Rows (first dimension) are the transaxial
%       detector direction, columns (second dimension) the axial direction.
%       options.nRowsD, nColsD, nProjections, angles (radians),
%       sourceToCRot, sourceToDetector, dPitchX, dPitchY, detOffsetRow and
%       detOffsetCol (nProjections x 1, mm) and flat (scalar).
%
%   Geometry notes: The XVI gantry angle a corresponds to the OMEGA
%   projection angle a + pi/2 (the source is at SAD*(sin a, -cos a, 0) in the
%   XVI frame, OMEGA places it at -SAD*(cos t, sin t, 0)). OMEGA's detector
%   row direction is opposite to the direction of increasing column index in
%   the .his files, hence the transaxial dimension of the images is
%   reversed. The last detector row of the XVI panel is dead and is replaced
%   by the previous row.
%
%   Example:
%       options.fpath = 'C:\data\XVI_scan'; options.binning = 2;
%       options = loadElektaXVIData(options);
%
% See also loadNikonData, loadProjectionImages
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2026 Ville-Veikko Wettenhovi
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(options, 'fpath') || isempty(options.fpath)
    fpath = uigetdir('', 'Select the Elekta XVI projection folder (.his files and _Frames.xml)');
    if isequal(fpath, 0)
        error('No folder was selected')
    end
else
    fpath = options.fpath;
    if exist(fpath, 'dir') == 0
        warning(['Input folder ' fpath ' not found!'])
        fpath = uigetdir('', 'Select the Elekta XVI projection folder (.his files and _Frames.xml)');
        if isequal(fpath, 0)
            error('No folder was selected')
        end
    end
end
if ~isfield(options, 'binning') || isempty(options.binning)
    binning = 1;
else
    binning = options.binning;
end
if isfield(options, 'only_reconstructions')
    only_reconstructions = options.only_reconstructions;
else
    only_reconstructions = false;
end
if ~isfield(options, 'elektaUCentreScale') || isempty(options.elektaUCentreScale)
    uScale = 1;
else
    uScale = options.elektaUCentreScale;
end
if ~isfield(options, 'elektaVCentreScale') || isempty(options.elektaVCentreScale)
    vScale = 1;
else
    vScale = options.elektaVCentreScale;
end
if ~isfield(options, 'elektaFlatPercentile') || isempty(options.elektaFlatPercentile)
    flatPercentile = 99.9;
else
    flatPercentile = options.elektaFlatPercentile;
end

%% Frames file
xmlFile = dir(fullfile(fpath, '*_Frames.xml'));
if isempty(xmlFile)
    error(['No *_Frames.xml file found in ' fpath])
end
txt = fileread(fullfile(fpath, xmlFile(1).name));
W = str2double(getTag(txt, 'Width'));
H = str2double(getTag(txt, 'Height'));
if isnan(W) || isnan(H)
    W = 512;
    H = 512;
end
if mod(W, binning) ~= 0 || mod(H, binning) ~= 0
    error('The binning factor has to divide the image size')
end
frames = regexp(txt, '<Frame>(.*?)</Frame>', 'tokens');
nFrames = numel(frames);
seq = zeros(nFrames, 1);
gantry = zeros(nFrames, 1);
uC = zeros(nFrames, 1);
vC = zeros(nFrames, 1);
exposed = false(nFrames, 1);
for kk = 1 : nFrames
    fr = frames{kk}{1};
    seq(kk) = str2double(getTag(fr, 'Seq'));
    gantry(kk) = str2double(getTag(fr, 'GantryAngle'));
    uC(kk) = str2double(getTag(fr, 'UCentre'));
    vC(kk) = str2double(getTag(fr, 'VCentre'));
    exposed(kk) = strcmpi(strtrim(getTag(fr, 'Exposed')), 'True');
end
[~, order] = sort(seq);
order = order(exposed(order));
seq = seq(order);
gantry = gantry(order);
uC = uC(order);
vC = vC(order);
nProj = numel(seq);

%% Geometry
if ~isfield(options, 'sourceToCRot') || isempty(options.sourceToCRot) || ~(options.sourceToCRot(1) > 1)
    options.sourceToCRot = 1000;
end
if ~isfield(options, 'sourceToDetector') || isempty(options.sourceToDetector) || ~(options.sourceToDetector(1) > 1)
    options.sourceToDetector = 1536;
end
if ~isfield(options, 'dPitchX') || isempty(options.dPitchX) || ~(options.dPitchX(1) > 0)
    options.dPitchX = 0.8 * binning;
end
if ~isfield(options, 'dPitchY') || isempty(options.dPitchY) || ~(options.dPitchY(1) > 0)
    options.dPitchY = 0.8 * binning;
end
options.nRowsD = W / binning;
options.nColsD = H / binning;
options.nProjections = nProj;
% The XVI gantry angle a is the OMEGA angle a + pi/2 (see above)
options.angles = deg2rad(gantry) + pi / 2;
options.detOffsetRow = uScale * uC;
options.detOffsetCol = -vScale * vC;

%% Projections
if ~only_reconstructions || ~isfield(options, 'SinM')
    files = dir(fullfile(fpath, '*.his'));
    fileSeq = nan(numel(files), 1);
    for kk = 1 : numel(files)
        tok = regexp(files(kk).name, '^(\d+)\.', 'tokens', 'once');
        if ~isempty(tok)
            fileSeq(kk) = str2double(tok{1});
        end
    end
    SinM = zeros(W / binning, H / binning, nProj, 'single');
    for kk = 1 : nProj
        ind = find(fileSeq == seq(kk), 1);
        if isempty(ind)
            error(['No .his file found for frame ' num2str(seq(kk))])
        end
        fname = fullfile(fpath, files(ind).name);
        % HIS file: 16-bit words at byte 0 (0x7000, file ID), byte 2 (header
        % size, 68) and byte 10 (image header size, 32); the image data
        % (uint16, little-endian, rows in the axial direction) follows
        % after both (100 bytes). If this does not match the file size, the
        % image is assumed to be the last part of the file.
        fid = fopen(fname, 'r', 'ieee-le');
        h = fread(fid, 6, 'uint16');
        hdrBytes = h(2) + h(6);
        if hdrBytes + W * H * 2 ~= files(ind).bytes
            hdrBytes = files(ind).bytes - W * H * 2;
        end
        fseek(fid, hdrBytes, 'bof');
        A = fread(fid, [W H], 'uint16=>single');
        fclose(fid);
        % A is (transaxial u, axial v). Raw value inversion: air has the
        % highest intensity
        A = single(65535) - A;
        % Dead last detector row
        A(:, H) = A(:, H - 1);
        % OMEGA's detector row direction is opposite to the file column order
        A = flip(A, 1);
        if binning > 1
            A = reshape(A, binning, W / binning, binning, H / binning);
            A = squeeze(sum(sum(A, 1), 3)) / single(binning^2);
        end
        SinM(:, :, kk) = A;
    end
    options.SinM = SinM;
    clear SinM
end

%% Global flat (air) value
if ~isfield(options, 'flat') || isempty(options.flat) || ~(options.flat(1) > 0)
    a = sort(options.SinM(:));
    k = ceil(flatPercentile / 100 * numel(a));
    k = min(max(k, 1), numel(a));
    options.flat = a(k);
    clear a
end
end

function val = getTag(txt, tag)
% Returns the text content of the first <tag>...</tag> in txt
tok = regexp(txt, ['<' tag '>([^<]*)</' tag '>'], 'tokens', 'once');
if isempty(tok)
    val = '';
else
    val = tok{1};
end
end
