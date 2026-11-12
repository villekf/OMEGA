function [factor, components] = quantificationFactorBiograph(normFile, acquisitionFile, options, varargin)
%QUANTIFICATIONFACTORBIOGRAPH Factor that converts a Siemens Biograph Vision
%image reconstructed with OMEGA to Bq/mL
%   Computes the scalar K such that
%       image_in_Bq_per_mL = K * x,
%   where x is an image that was reconstructed with OMEGA from measured
%   Vision data with the normalization of loadNormalizationBiograph (the
%   Siemens PET_CALIBRATION .ptd file, options.use_user_normalization = true)
%   and with the same attenuation, randoms and scatter corrections as the
%   Siemens reconstruction. The factor is
%
%       K = C * QF * xcal * DF * DT / (BF * T * mRef),
%
%   where
%       QF   = scanner quantification factor (Bq*s/ECAT counts) of the
%              normalization file header,
%       xcal = cross calibration factor of the normalization file header,
%       BF   = branching factor of the isotope (acquisition header),
%       T    = frame duration in seconds (acquisition header),
%       DF   = lambda * T / (1 - exp(-lambda * T)), the decay of the isotope
%              during the frame, lambda = log(2) / half-life. The reference
%              time is the start of the frame (the Siemens decay correction
%              to the scan start); decay from the time of the injection is
%              NOT applied,
%       DT   = 1 / L^2, the dead-time factor. L is the mean live-time
%              fraction of the crystals, from the dead-time table (component
%              5) of the normalization file and the block singles rates of
%              the acquisition header (the rate of a block divided by 210
%              is the table index of its 20 x 10 crystals; the divisor 210
%              is empirical and reproduces the Siemens dead-time factors),
%       mRef = the mean efficiency of the LORs of the normalization with the
%              LOR length factor (Dref / D)^n, i.e. the third output of
%              loadNormalizationBiograph, which depends on the sinogram size
%              (Ndist) and on options.visionLORLengthExponent (n) and
%              Dref = options.diameter + 2 * options.DOI,
%       C    = an empirical constant, see below.
%   The result multiplies the final image (e.g. pz of the reconstruction),
%   it is not applied inside the reconstruction.
%
%   The constant C is empirical (it has no physical interpretation). It was
%   determined by comparing OMEGA reconstructions with the Siemens (e-soft)
%   reconstructions in Bq/mL of three measured datasets (a Derenzo phantom,
%   a torso phantom and a pig, mean of 2.178e-3, 2.218e-3 and 2.183e-3).
%
%   Validity of C: it was determined with non-TOF reconstructions (projector
%   type 1, 3 iterations of 10 subsets) of sinograms with 520 radial bins and
%   span 19, with the attenuation, normalization, randoms and scatter
%   corrections. For the same data the factor was unchanged (within 0.2 %) by
%   the projector type (1, 2 or 3) of non-TOF reconstructions, and, through
%   mRef, within 1 % by Ndist (520 or 400), the angular mashing, the span (1 or
%   19) and options.visionLORLengthExponent (0 or 2). C has not been validated
%   with TOF reconstructions. The TOF projectors have the same scale as the
%   non-TOF ones, but TOF images need the randoms and scatter corrected per TOF
%   bin, and uncorrected randoms or scatter change the image level differently
%   with and without TOF. The factor has only been validated for non-TOF
%   reconstructions.
%
%   Siemens headers can be wrong about the isotope. The isotope name, the
%   half-life and the branching factor are read from the acquisition header
%   and can be replaced with the name-value pairs below, e.g. a Ga-68 scan
%   that was labelled F-18 needs 'halfLife', 4062.6, 'branchingFactor',
%   0.8891.
%
%   Only the Biograph Vision is supported (the dead-time table layout of the
%   mCT is different).
%
% Inputs:
%   normFile = path to the Vision normalization file, i.e. the
%   PET_CALIBRATION .ptd file (the header with QF and xcal is required)
%   acquisitionFile = path to the PET_EM_SINO .ptd file of the frame. The
%   header of a sinogram file contains the isotope, the frame duration and
%   the block singles rates. A PET_LISTMODE .ptd file contains the isotope
%   and the duration, but not the singles rates: a warning is then shown
%   and DT = 1 is used.
%   options = the options struct that is used for the reconstruction (the
%   sinogram size and the detector geometry as in loadNormalizationBiograph)
%
% Name-value pairs (optional):
%   'halfLife' = half-life of the isotope (s), replaces the header value
%   'branchingFactor' = branching factor, replaces the header value
%   'frameDuration' = frame duration (s), replaces the header value
%
% Outputs:
%   factor = K, double. Multiply the reconstructed image with it to get Bq/mL
%   components = struct with the fields C, QF, xcal, DF, DT, BF, T,
%   halfLife, mRef, livetime (L, NaN if the singles are not available),
%   isotope (the isotope name of the header), Dref and singlesFound
%
% Usage:
%   [K, comp] = quantificationFactorBiograph(options.normalization_datafile, ...
%       acquisitionFile, options);
%   pz = pz * K;
%   % Ga-68 scan that is labelled as F-18
%   K = quantificationFactorBiograph(normFile, acquisitionFile, options, ...
%       'halfLife', 4062.6, 'branchingFactor', 0.8891);

% Empirical constant (see above), unitless
C = 2.19e-3;
% Divisor of the block singles rate (empirical, see above)
rateDivisor = 210;
% Vision block layout
nAxialBlocks = 8;
nTransBlocks = 38;
cryst_ax = 10; % Axial crystals per block
cryst_tr = 20; % Transaxial crystals per block
nRings = nAxialBlocks * cryst_ax;
nCryst = nTransBlocks * cryst_tr;

% Optional overrides
halfLifeIn = [];
bfIn = [];
tIn = [];
if mod(numel(varargin), 2) ~= 0
    error('quantificationFactorBiograph:args', 'The optional inputs have to be name-value pairs')
end
for kk = 1 : 2 : numel(varargin)
    if ~ischar(varargin{kk})
        error('quantificationFactorBiograph:args', 'The optional input names have to be character vectors')
    end
    val = varargin{kk + 1};
    if ~isnumeric(val) || ~isscalar(val) || ~isfinite(val) || val <= 0
        error('quantificationFactorBiograph:args', 'The value of ''%s'' has to be a positive finite scalar', varargin{kk})
    end
    switch lower(varargin{kk})
        case 'halflife'
            halfLifeIn = double(val);
        case 'branchingfactor'
            bfIn = double(val);
        case 'frameduration'
            tIn = double(val);
        otherwise
            error('quantificationFactorBiograph:args', 'Unknown input ''%s''', varargin{kk})
    end
end

% Normalization file header
if exist(normFile, 'file') ~= 2
    error('quantificationFactorBiograph:file', 'The normalization file %s does not exist', normFile)
end
if exist(acquisitionFile, 'file') ~= 2
    error('quantificationFactorBiograph:file', 'The acquisition file %s does not exist', acquisitionFile)
end
txtN = readHeaderText(normFile, true);
if isempty(txtN)
    error('quantificationFactorBiograph:header', ['No Interfile header was found in %s. The normalization file has to be the ' ...
        'PET_CALIBRATION .ptd file of the Vision.'], normFile)
end
tokO = regexp(txtN, 'data offset in bytes\[(\d+)\]\s*:=\s*(\d+)', 'tokens');
tokS = regexp(txtN, 'matrix size\[(\d+)\]\s*:=\s*\{([^}]*)\}', 'tokens');
offs = zeros(1, 5);
dims = cell(1, 5);
for kk = 1 : numel(tokO)
    k = str2double(tokO{kk}{1});
    if k <= 5
        offs(k) = str2double(tokO{kk}{2});
    end
end
for kk = 1 : numel(tokS)
    k = str2double(tokS{kk}{1});
    if k <= 5
        dims{k} = sscanf(strrep(tokS{kk}{2}, ',', ' '), '%d')';
    end
end
if numel(tokO) < 5 || numel(tokS) < 5 || ~isequal(dims{3}, [798 nRings]) || ~isequal(dims{5}, [6000 42 nRings])
    error('quantificationFactorBiograph:scanner', ['%s is not a Biograph Vision normalization file with the dead-time table ' ...
        '(five components, crystal efficiencies 798 x 80, dead-time table 6000 x 42 x 80). Only the Vision is supported.'], normFile)
end
QF = headerNumber(txtN, 'scanner quantification factor \(Bq\*s/ECAT counts\)', 'the scanner quantification factor', normFile);
xcal = headerNumber(txtN, 'cross calibration factor', 'the cross calibration factor', normFile);

% Acquisition header
txtA = readHeaderText(acquisitionFile, false);
if isempty(txtA)
    error('quantificationFactorBiograph:header', ['No Interfile header was found in %s. The acquisition file has to be the ' ...
        'PET_EM_SINO (or PET_LISTMODE) .ptd file.'], acquisitionFile)
end
isotope = headerString(txtA, 'isotope name');
if isempty(halfLifeIn)
    halfLife = headerNumber(txtA, 'isotope gamma halflife \(sec\)', 'the half-life (give ''halfLife'')', acquisitionFile);
else
    halfLife = halfLifeIn;
end
if isempty(bfIn)
    BF = headerNumber(txtA, 'isotope branching factor', 'the branching factor (give ''branchingFactor'')', acquisitionFile);
else
    BF = bfIn;
end
if isempty(tIn)
    T = headerNumber(txtA, 'image duration \(sec\)', 'the frame duration (give ''frameDuration'')', acquisitionFile);
else
    T = tIn;
end
if ~(halfLife > 0 && BF > 0 && T > 0)
    error('quantificationFactorBiograph:header', 'The half-life, branching factor and frame duration have to be positive')
end

% Decay during the frame, referred to the start of the frame
lam = log(2) / halfLife;
DF = lam * T / (-expm1(-lam * T));

% Dead-time factor
singles = [];
tokR = regexp(txtA, 'singles rate\[(\d+)\]\s*:=\s*\{([^}]*)\}', 'tokens');
if ~isempty(tokR)
    nAx = headerNumber(txtA, 'number of axial blocks', 'the number of axial blocks', acquisitionFile);
    nTr = headerNumber(txtA, 'number of transaxial blocks', 'the number of transaxial blocks', acquisitionFile);
    if nAx ~= nAxialBlocks || nTr ~= nTransBlocks
        error('quantificationFactorBiograph:scanner', ['The acquisition header has %d x %d detector blocks, the Vision has %d x %d. ' ...
            'Only the Vision is supported.'], nAx, nTr, nAxialBlocks, nTransBlocks)
    end
    singles = zeros(nAxialBlocks * nTransBlocks, numel(tokR));
    for kk = 1 : numel(tokR)
        v = sscanf(strrep(tokR{kk}{2}, ',', ' '), '%f');
        % The first value is not a block (the blocks follow it)
        if numel(v) == nAxialBlocks * nTransBlocks + 1
            v = v(2:end);
        end
        if numel(v) ~= nAxialBlocks * nTransBlocks
            error('quantificationFactorBiograph:header', 'The singles rate of the acquisition header has %d values, %d expected', ...
                numel(v), nAxialBlocks * nTransBlocks)
        end
        singles(:, kk) = v;
    end
    % Several sampling points: the mean rate of the frame
    singles = mean(singles, 2);
end
if isempty(singles)
    warning('quantificationFactorBiograph:noSingles', ['The block singles rates were not found in the header of %s (a list-mode ' ...
        'file does not contain them), the dead-time factor is set to one. Use the PET_EM_SINO file for the dead-time factor.'], acquisitionFile)
    DT = 1;
    livetime = NaN;
else
    % Block rates [axial block, transaxial block]
    rates = reshape(singles, nTransBlocks, nAxialBlocks)';
    fid = fopen(normFile, 'r');
    fseek(fid, offs(5), 'bof');
    tab = fread(fid, 6000 * 42 * nRings, 'single=>double', 0, 'l');
    fclose(fid);
    if numel(tab) ~= 6000 * 42 * nRings
        error('quantificationFactorBiograph:file', 'The normalization file ended before the end of the dead-time table')
    end
    tab = reshape(tab, 6000, 42, nRings);
    % Live-time fraction of every crystal: the table index is the block rate / 210 (rounded half up) and the
    % transaxial position of the table skips the gap crystals 0 and 21; odd blocks use the positions 22-41
    ring = (0 : nRings - 1)';
    t = 0 : nCryst - 1;
    x = rates(floor(ring / cryst_ax) + 1, floor(t / cryst_tr) + 1) / rateDivisor;
    xi = min(max(floor(x + 0.5), 0), 5999);
    pos = mod(floor(t / cryst_tr), 2) * 21 + 1 + mod(t, cryst_tr);
    idx = sub2ind([6000 42 nRings], xi + 1, repmat(pos + 1, nRings, 1), repmat(ring + 1, 1, nCryst));
    livetime = mean(tab(idx(:)));
    if ~(livetime > 0) || ~isfinite(livetime)
        error('quantificationFactorBiograph:header', 'The live-time fraction from the dead-time table is not positive')
    end
    DT = 1 / livetime^2;
end

% Mean efficiency of the LORs with the length factor (Dref / D)^n
[~, ~, mRef] = loadNormalizationBiograph(normFile, options);
if isempty(mRef) || ~(mRef > 0)
    error('quantificationFactorBiograph:scanner', 'The mean efficiency of the normalization could not be determined')
end

factor = C * QF * xcal * DF * DT / (BF * T * mRef);

if isfield(options, 'DOI') && ~isempty(options.DOI)
    Dref = double(options.diameter(1)) + 2 * double(options.DOI(1));
else
    Dref = double(options.diameter(1));
end
components = struct('C', C, 'QF', QF, 'xcal', xcal, 'DF', DF, 'DT', DT, 'BF', BF, 'T', T, 'halfLife', halfLife, ...
    'mRef', mRef, 'livetime', livetime, 'isotope', isotope, 'Dref', Dref, 'singlesFound', ~isempty(singles));
end

function txt = readHeaderText(f_path, sidecar)
%READHEADERTEXT Reads the Interfile header that is embedded in a .ptd file
%(at the end of the file; the beginning of the file is checked as well) or,
%if sidecar is true, stored in the sidecar <file>.hdr. Returns an empty
%string if there is no header. Same logic as readNormHeader in
%loadNormalizationBiograph.
fid = fopen(f_path, 'r');
if fid < 0
    error('quantificationFactorBiograph:file', 'Could not open the file %s', f_path)
end
fseek(fid, 0, 'eof');
fsz = ftell(fid);
nBuf = min(fsz, 4 * 2^20);
fseek(fid, -nBuf, 'eof');
buf = fread(fid, nBuf, 'uint8=>uint8')';
pos = strfind(char(buf), '!INTERFILE');
if isempty(pos) && fsz > nBuf
    fseek(fid, 0, 'bof');
    buf = fread(fid, nBuf, 'uint8=>uint8')';
    pos = strfind(char(buf), '!INTERFILE');
end
fclose(fid);
txt = '';
if ~isempty(pos)
    txt = char(buf(pos(1):end));
    % The header ends where the binary (DICOM) data continues
    z = find(txt == 0, 1);
    if ~isempty(z)
        txt = txt(1:z-1);
    end
elseif sidecar
    [fp, fn] = fileparts(f_path);
    hdrs = {[f_path '.hdr'], fullfile(fp, [fn '.hdr'])};
    for kk = 1 : numel(hdrs)
        if exist(hdrs{kk}, 'file') == 2
            txt = fileread(hdrs{kk});
            break
        end
    end
end
end

function s = headerString(txt, key)
%HEADERSTRING Value of the Interfile key (regular expression) as a string,
%empty if the key is not in the header
tok = regexp(txt, ['^[!%]?' key '\s*:=[ \t]*([^\r\n]*)'], 'tokens', 'once', 'lineanchors');
if isempty(tok)
    s = '';
else
    s = strtrim(tok{1});
end
end

function v = headerNumber(txt, key, what, f_path)
%HEADERNUMBER Numeric value of the Interfile key (regular expression), error
%if the key is missing or not a number
s = headerString(txt, key);
v = str2double(s);
if isempty(s) || ~isfinite(v)
    error('quantificationFactorBiograph:header', 'The header of %s does not contain %s', f_path, what)
end
end
