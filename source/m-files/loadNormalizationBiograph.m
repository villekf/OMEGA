function [normalization_u, efficiency, meanEfficiency] = loadNormalizationBiograph(f_path, options)
%LOADNORMALIZATIONBIOGRAPH Loads the normalization coefficients created by
%Biograph machine, either the mCT or Vision
%   Reads the component-based normalization file of the scanner (.n or the
%   .ptd file of the Vision calibration) and converts it into sinogram
%   normalization coefficients with the sinogram dimensions of OMEGA.
%
%   The component layout is read from the Interfile header that follows the
%   data in a .ptd file or, for a plain .n file, from the sidecar header
%   <file>.hdr. Without a header the mCT layout is assumed for small files
%   and the standard Vision layout for large ones.
%
%   Vision: the efficiency of every LOR is
%       eff = E(s1,ring1) * E(s2,ring2) * I(p,r) / (G(r,P) * A(P)),
%   where E are the crystal efficiencies, I the crystal interference, G the
%   geometric effects and A the axial effects of the span-1 plane P of the
%   scanner. The efficiencies of all LORs that are binned into the same
%   sinogram bin by OMEGA (same binning as in saveSinogram, including the
%   mashing, span and the sinogram size) are summed, i.e. the multiplicity of
%   the bin is included, and the sum is normalized such that the mean
%   efficiency of the LORs is one. The component 5 of the file (dead-time
%   parameters) depends on the count rate and is ignored; a global dead-time
%   factor can be given with options.global_correction_factor.
%
%   The efficiency of every LOR is also multiplied by D^(-n), where D is the
%   transaxial (x-y) distance between the two detectors of the LOR and
%   n = options.visionLORLengthExponent (default 2, 0 = no length factor).
%   The factor is empirical: it was validated against the Siemens
%   reconstruction of measured Vision data, where the radial RMS error
%   decreased from 4.66 to 2.15 % with n = 2. The detector positions are
%   those of detector_coordinates, so the detector geometry of the options
%   (diameter, DOI, cr_p, blocks_per_ring, cryst_per_block, transaxialGaps,
%   firstBlockAngle, flip_image and offangle) has to be set.
%
%   Vision requirements: 80 rings, 760 detectors per ring without or 798 with
%   pseudo detectors (one at the end of each block of 20 crystals). The
%   Ndist, Nang, span, ring_difference, segment_table and ndist_side can be
%   freely selected.
%
%   LORs whose radial distance is outside of the 520 radial bins of the
%   scanner (possible with the 760 crystal sinogram and a large Ndist) are
%   given the crystal interference one and the geometric effect of the
%   nearest edge bin. A warning with the number of such LORs and their
%   share of the total efficiency is shown.
%
% Inputs:
%   f_path = path to the normalization file
%   options = the options struct with the scanner properties (Ndist, Nang,
%   rings, det_per_ring, det_w_pseudo, cryst_per_block, span,
%   ring_difference, segment_table, ndist_side and visionLORLengthExponent
%   are used for Vision)
%
% Outputs:
%   normalization_u = the normalization coefficients, single
%   [Ndist, Nang, TotSinos]. For Vision the coefficients MULTIPLY the
%   measured sinogram (normalized = measured .* normalization_u), and
%   1 ./ normalization_u is the efficiency (options.normalization during the
%   reconstruction). Bins without any efficiency (e.g. only pseudo detector
%   LORs) are zero and their efficiency should be taken as zero, not Inf.
%   efficiency = (optional) the summed efficiency of each sinogram bin,
%   i.e. the Vision normalization_u is 1 ./ efficiency (zero where the
%   efficiency is zero).
%   meanEfficiency = (optional, Vision only, empty otherwise) the scale that
%   was divided out when the mean efficiency of the LORs was normalized to
%   one: the mean over the LORs of E1 * E2 * I / (G * A) * (Dref / D)^n,
%   with Dref = options.diameter + 2 * options.DOI (the effective detector
%   ring diameter in mm) and n = options.visionLORLengthExponent. With n = 0
%   it is the mean of the plain efficiencies. It is used by
%   quantificationFactorBiograph to convert the images to Bq/mL.

efficiency = [];
meanEfficiency = [];
comp = readNormHeader(f_path);
if isVisionNorm(comp, f_path)
    [normalization_u, efficiency, meanEfficiency] = visionNormalization(f_path, comp, options);
    return
end
fid = fopen(f_path);
if fid < 0
    error('Could not open the normalization file %s', f_path)
end
if isempty(comp)
    nRead = inf;
else
    % Only the normalization components, not the header that follows them
    nRead = max([comp.offset] + [comp.nbytes]) / 4;
end
normal = fread(fid, nRead, 'single=>single',0,'l');
fclose(fid);
if size(normal,1) <= 1255709
    mashing = options.det_w_pseudo / options.Nang / 2;
    det_w_pseudo = options.det_w_pseudo;
    no_pseudo = 0;
    if det_w_pseudo == options.det_per_ring
        det_w_pseudo = det_w_pseudo + 1;
        no_pseudo = 1;
    end
    normalization = ones(options.Ndist, options.Nang*mashing, options.TotSinos,'single');
    geom = reshape(normal(1:options.Ndist*options.Nz), options.Ndist, options.Nz);
    alku = options.Ndist*options.Nz + 1;
    interf = reshape(normal(alku:alku + (det_w_pseudo / options.blocks_per_ring) * options.Ndist - 1), ...
        det_w_pseudo / options.blocks_per_ring, options.Ndist);
    if no_pseudo
        interf = interf(1:options.cryst_per_block,:);
    end
    alku = alku + (det_w_pseudo / options.blocks_per_ring) * options.Ndist;
    ceff = reshape(normal(alku:alku + det_w_pseudo * options.rings - 1), ...
        det_w_pseudo, options.rings);
    ceff = 1./ceff;
    if no_pseudo
        ceff(options.cryst_per_block + 1 : options.cryst_per_block + 1 : end,:) = [];
    end
%     ceff = ceff / mean(ceff(:));
    alku = alku + det_w_pseudo * options.rings;
    axial = normal(alku:alku + options.TotSinos-1);
    alku = alku + options.TotSinos;
    parringDT = normal(alku:alku + options.rings-1);
    alku = alku + options.rings;
    nonparringDT = normal(alku:alku + options.rings-1);
    alku = alku + options.rings;
    TXcrystDT = normal(alku:end);
    z = sinogram_coordinates_3D(options);
    zz = z(1:options.Nz,1);
    z = z(options.Nz + 1: end,:);
    [~,loc] = ismember(z,zz);
    geom = [geom, (geom(:,loc(:,1)) + geom(:,loc(:,2))) / 2];
    geom = permute(geom, [1 3 2]);
    normalization = bsxfun(@rdivide, normalization, geom);
%     if mashing > 1
%         interf = reshape(accumarray(ceil((1:numel(interf))'/2),interf(:),[],@mean), (options.det_w_pseudo / options.blocks_per_ring) / 2, options.Ndist);
%     end
    interf = repmat(1./interf, options.det_w_pseudo / size(interf,1) / 2, 1)';
    %     interf = repmat(interf, options.blocks_per_ring, 1);
    axial = permute(axial, [3 2 1]);
    normalization = bsxfun(@times, normalization, axial);
    
    if mashing > 1
        options.Nang = options.Nang * mashing;
    end
    L = zeros(sum(1:options.det_w_pseudo),2,'int32');
    jh = int32(1);
    for kk = int32(1) : (options.det_w_pseudo)
        if exist('OCTAVE_VERSION','builtin') == 0 && exist('repelem', 'builtin') == 0
            L(jh:(jh + (options.det_w_pseudo) - kk),:) = [repeat_elem((kk), options.det_w_pseudo-(kk-1)), ((kk):options.det_w_pseudo)'];
        elseif exist('OCTAVE_VERSION','builtin') == 5
            L(jh:(jh + (options.det_w_pseudo) - kk),:) = [repelem((kk), options.det_w_pseudo-(kk-1)), ((kk):options.det_w_pseudo)'];
        else
            L(jh:(jh + (options.det_w_pseudo) - kk),:) = [repelem((kk), options.det_w_pseudo-(kk-1))', ((kk):options.det_w_pseudo)'];
        end
        jh = jh + (options.det_w_pseudo) -kk + 1;
    end
    L(L(:,1) == 0,:) = [];
    
    L = L - 1;
    
    xa = max(L,[],2);
    ya = min(L,[],2);
    
    j = idivide(mod(xa+ya+options.det_w_pseudo/2+1,options.det_w_pseudo),2);
    
    b = j+options.det_w_pseudo/2;
    
    i = abs(xa-ya-options.det_w_pseudo/2);
    for kk = 1 : length(ya)
        if (ya(kk)<j(kk)) || (b(kk)<xa(kk))
            i(kk) = -i(kk);
        end
    end
    
    % Determine the accepted LORs (distances that are within the predefined
    % value)
    if mod(options.Ndist,2) == 0
        accepted_lors = (i <= (options.Ndist/2 + min(0,options.ndist_side)) & i >= (-options.Ndist/2 + max(0,options.ndist_side)));
    else
        accepted_lors = (i <= options.Ndist/2 & i >= (-options.Ndist/2));
    end
    
    j = idivide(j,options.det_w_pseudo/2/options.Nang);
    
    i = i(accepted_lors);
    j = j(accepted_lors);
    if min(i) <= 0
        i = i + abs(min(i)) + 1;
    end
    j = j + 1;
    
    L = L(accepted_lors,:);
    
    L = L + 1;
%     if mashing > 1
%         options.Nang = options.Nang / mashing;
%     end
    
    ceff_n1 = ceff(L(:,1),:);
    ceff_n2 = ceff(L(:,2),:);
    cell_ceff1 = arrayfun(@(x) accumarray([i j],ceff_n1(:,x)),1:size(ceff_n1,2),'un',0);
    cell_ceff2 = arrayfun(@(x) accumarray([i j],ceff_n2(:,x)),1:size(ceff_n2,2),'un',0);
    ceff_nn = cellfun(@times, repmat(cell_ceff1,size(cell_ceff1,2),1), repmat(cell_ceff2',1,size(cell_ceff2,2)), 'UniformOutput', false);
%     if mashing > 1
%         apu = cell2mat(ceff_nn);
%         ceff_nn = squeeze(mean(reshape(apu',[mashing,size(apu,2)/mashing,size(apu,1)]),1))';
%         [r, c, ~] = size(ceff_nn);
%         ceff_nn = mat2cell(ceff_nn, options.Ndist * ones(1,r/options.Ndist), options.Nang * ones(1,c/options.Nang));
%     end
    ceff_nn = cat(3,ceff_nn{:});
    ceff_nn = bsxfun(@times, ceff_nn, interf);
    ceff = ones(options.Ndist,options.Nang,options.TotSinos,'single');
    
    %     interf_n1 = interf(L(:,1),:);
    %     interf_n2 = interf(L(:,2),:);
    %     interfx = accumarray([i j], interf_n1(:,1), [options.Ndist options.Nang],@mean, NaN);
    %     interfy = accumarray([i j], interf_n2, [options.Ndist options.Nang],@mean, NaN);
    %     normalization = bsxfun(@rdivide, normalization, interf');
    
    kkj = zeros(floor((options.ring_difference-ceil(options.span/2))/options.span) + 1, 1);
    for kk = 1 : floor((options.ring_difference-ceil(options.span/2))/options.span) + 1
        kkj(kk) = ceil(options.span/2) + options.span*(kk - 1);
    end
    offset2 = cumsum(options.segment_table);
    % Create the sinograms
    % First the detectors on the same ring
    ceff(:,:,1:2:options.Nz) = ceff_nn(:,:,1:options.rings+1:options.rings^2);
    % Then the detectors on adjacent rings
    mean_jh = zeros(options.TotSinos,1);
    mean_jh(1:2:options.Nz) = 1;
    for jh=1:floor(options.span/2)
        apu = ceff_nn(:,:,jh*options.rings+1:options.rings+1:options.rings^2);
        apu2 = ceff_nn(:,:,jh+1:options.rings+1:(options.rings-jh)*options.rings);
        loc = ismember(1:options.Nz,jh+1:2:offset2(1)-jh);
        mean_jh(loc) = mean_jh(loc) + 2;
        ceff(:,:,jh+1:2:offset2(1)-jh) = ceff(:,:,jh+1:2:offset2(1)-jh) .* apu .* apu2;
    end
    %     mean_jh(mean_jh == 0) = 1;
    %     mean_jh = permute(mean_jh, [3 2 1]);
    %     ceff(:,:,1:options.Nz) = bsxfun(@rdivide, ceff(:,:,1:options.Nz), mean_jh);
    % Lastly the rest of the detectors with the amount of combined LORs
    % specified with the span value
    %     mean_jh = zeros(options.TotSinos,1);
    for ih=1:floor(length(options.segment_table)/2)
        for jh=1:options.span
            apu = ceff_nn(:,:,(kkj(ih)+jh-1)*options.rings+1:options.rings+1:end);
            loc = ismember(1:options.TotSinos,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1);
            mean_jh(loc) = mean_jh(loc) + 1;
            ceff(:,:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) = ceff(:,:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) .* (apu);
            apu2 = ceff_nn(:,:,kkj(ih)+jh:options.rings+1:(options.rings-kkj(ih)-jh+1)*options.rings);
            loc = ismember(1:options.TotSinos,offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1);
            mean_jh(loc) = mean_jh(loc) + 1;
            ceff(:,:,offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1) = ceff(:,:,offset2(2*ih)+jh:2:offset2(2*ih + 1)-jh+1) .* (apu2);
        end
    end
    mean_jh(mean_jh == 0) = 1;
    mean_jh = permute(mean_jh, [3 2 1]);
    ceff = min(ceff(:));
    ceff = bsxfun(@rdivide, ceff, mean_jh);
%     normalization = normalization ./ ceff;
    normalization = normalization .* ceff;
%     fad = fopen('norm3d_00.a');
%     A = fread(fad, inf, 'single=>single');
%     normalization_u = reshape(A, options.Ndist, options.Nang, options.TotSinos);
    normalization_u = cell2mat(arrayfun(@(k) mean(normalization(:,k:min(size(normalization,1),k+2-1),:),2), 1:2:size(normalization,2), 'un', 0));
    %     normalization = circshift(normalization,1);
    %     normalization = 1 / normalization;
else
    error('Unrecognized normalization file layout (%d floats and no usable header)', size(normal,1))
end
end

function comp = readNormHeader(f_path)
%READNORMHEADER Reads the component offsets and sizes from the Interfile
%header (embedded in the .ptd file or in the sidecar <file>.hdr). Returns an
%empty struct array if there is no header.
comp = struct('offset', {}, 'dims', {}, 'nbytes', {});
fid = fopen(f_path, 'r');
if fid < 0
    error('Could not open the normalization file %s', f_path)
end
fseek(fid, 0, 'eof');
fsz = ftell(fid);
nTail = min(fsz, 4 * 2^20);
fseek(fid, -nTail, 'eof');
buf = fread(fid, nTail, 'uint8=>uint8')';
pos = strfind(char(buf), '!INTERFILE');
if isempty(pos) && fsz > nTail
    fseek(fid, 0, 'bof');
    buf = fread(fid, fsz, 'uint8=>uint8')';
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
else
    [fp, fn] = fileparts(f_path);
    hdrs = {[f_path '.hdr'], fullfile(fp, [fn '.hdr'])};
    for kk = 1 : numel(hdrs)
        if exist(hdrs{kk}, 'file') == 2
            txt = fileread(hdrs{kk});
            break
        end
    end
end
if isempty(txt)
    return
end
tokO = regexp(txt, 'data offset in bytes\[(\d+)\]\s*:=\s*(\d+)', 'tokens');
tokS = regexp(txt, 'matrix size\[(\d+)\]\s*:=\s*\{([^}]*)\}', 'tokens');
if isempty(tokO) || numel(tokO) ~= numel(tokS)
    warning('The normalization header could not be parsed, using the default layout')
    return
end
for kk = 1 : numel(tokO)
    k = str2double(tokO{kk}{1});
    comp(k).offset = str2double(tokO{kk}{2});
end
for kk = 1 : numel(tokS)
    k = str2double(tokS{kk}{1});
    comp(k).dims = sscanf(strrep(tokS{kk}{2}, ',', ' '), '%d')';
    comp(k).nbytes = 4 * prod(comp(k).dims);
end
if any(cellfun(@isempty, {comp.offset})) || any(cellfun(@isempty, {comp.dims}))
    warning('The normalization header could not be parsed, using the default layout')
    comp = struct('offset', {}, 'dims', {}, 'nbytes', {});
end
end

function tf = isVisionNorm(comp, f_path)
%ISVISIONNORM True if the normalization file has the Vision layout
if ~isempty(comp)
    tf = numel(comp) >= 4 && isequal(comp(3).dims, [798 80]);
else
    d = dir(f_path);
    tf = d.bytes / 4 > 1255709;
end
end

function x = readComponent(fid, offset, dims)
%READCOMPONENT Reads one float32 component (first dimension fastest)
fseek(fid, offset, 'bof');
x = fread(fid, prod(dims), 'single=>double', 0, 'l');
if numel(x) ~= prod(dims)
    error('The normalization file ended before the end of the normalization component')
end
x = reshape(x, dims);
end

function [t, v] = siemensBin(s1, s2)
%SIEMENSBIN Siemens (CTI) sinogram bin of the crystal pair (s1, s2) of a
%ring of 798 crystals, oriented such that s1 = mod(v + floor(t/2), 798) and
%s2 = mod(v - floor((t+1)/2) + 399, 798). The pair is swapped if needed.
t = mod(s1 - s2 + 399, 798);
t(t >= 399) = t(t >= 399) - 798;
v = mod(s1 + s2 - 399 + mod(t, 2), 798) / 2;
bad = mod(v + floor(t / 2), 798) ~= s1;
if any(bad)
    a = s2(bad);
    b = s1(bad);
    t(bad) = mod(a - b + 399, 798);
    t(bad & t >= 399) = t(bad & t >= 399) - 798;
    v(bad) = mod(a + b - 399 + mod(t(bad), 2), 798) / 2;
end
end

function sp = sinogramPlane(rn1, rn2, span, rings, seg)
%SINOGRAMPLANE Zero-based sinogram plane of the ring numbers rn1 and rn2
%(the same rule as in saveSinogram)
if span <= 1
    sp = rn2 * rings + rn1;
else
    h = floor(span / 2);
    erotus = rn1 - rn2;
    summa = rn1 + rn2;
    sp = summa;
    far = abs(erotus) > h;
    k = floor((abs(erotus) + h) / span);
    if any(far(:)) && max(k(far)) * 2 > numel(seg)
        error('The segment table has too few segments for the ring difference')
    end
    off = zeros(size(sp));
    neg = far & erotus < 0;
    pos = far & erotus >= 0;
    off(neg) = seg((k(neg) - 1) * 2 + 1);
    off(pos) = seg((k(pos) - 1) * 2 + 2);
    sp(far) = summa(far) - (h * (2 * k(far) - 1) + k(far)) + off(far);
end
end

function h = lorLengthFactor(options, W, x1, x2, n)
%LORLENGTHFACTOR Transaxial LOR length factor D^(-n) of the detector pairs
%(x1, x2) (zero-based), D = distance between the two detectors, from the
%detector coordinates of OMEGA (detector_coordinates). Pseudo detectors have
%finite coordinates; non-finite factors are set to zero.
req = {'cr_p', 'diameter', 'blocks_per_ring', 'cryst_per_block'};
for kk = 1 : numel(req)
    if ~isfield(options, req{kk}) || isempty(options.(req{kk})) || any(options.(req{kk})(:) <= 0)
        error(['The Vision LOR length factor needs the detector geometry, but options.%s is missing or not positive. ' ...
            'Set the detector geometry or set options.visionLORLengthExponent = 0 to disable the factor.'], req{kk})
    end
end
if ~isfield(options, 'flip_image') || isempty(options.flip_image)
    options.flip_image = false;
end
if ~isfield(options, 'offangle') || isempty(options.offangle)
    options.offangle = 0;
end
% Double precision coordinates
options.listmode = false;
options.implementation = 4;
[~, ~, xp, yp] = detector_coordinates(options);
xp = double(xp(:));
yp = double(yp(:));
if numel(xp) ~= W
    error(['The number of transaxial detector coordinates (%d) differs from det_w_pseudo (%d). ' ...
        'Set the detector geometry correctly or set options.visionLORLengthExponent = 0 to disable the LOR length factor.'], numel(xp), W)
end
h = hypot(xp(x1 + 1) - xp(x2 + 1), yp(x1 + 1) - yp(x2 + 1)) .^ (-n);
h(~isfinite(h)) = 0;
end

function [coef, effSum, meanEff] = visionNormalization(f_path, comp, options)
%VISIONNORMALIZATION Sinogram normalization coefficients of the Vision from
%the normalization components

NR = 520; % Radial bins of the scanner
NS = 798; % Crystals per ring, including the gaps
R = 80; % Rings
% Required options
Nd = double(options.Ndist);
Nang = double(options.Nang);
rings = double(options.rings);
W = double(options.det_w_pseudo);
span = double(options.span);
rd = double(options.ring_difference);
if isfield(options, 'ndist_side') && ~isempty(options.ndist_side)
    nSide = double(options.ndist_side);
else
    nSide = 1;
end
if isfield(options, 'visionLORLengthExponent') && ~isempty(options.visionLORLengthExponent)
    nLen = double(options.visionLORLengthExponent);
else
    nLen = 2;
end
if ~isscalar(nLen) || ~isreal(nLen) || ~isfinite(nLen)
    error('options.visionLORLengthExponent has to be a finite real scalar')
end
if rings ~= R
    error('The Vision normalization requires 80 crystal rings')
end
if W == options.det_per_ring
    if W ~= 760 || options.cryst_per_block ~= 20
        error('The Vision normalization requires 760 detectors per ring (38 blocks of 20 crystals)')
    end
elseif W ~= NS
    error('The Vision normalization requires det_w_pseudo of 760 (no pseudo detectors) or 798 (one pseudo detector per block)')
end
if span <= 1
    totSinos = rings^2;
    seg = [];
else
    seg = cumsum(double(options.segment_table(:)));
    totSinos = seg(end);
end
if rd > rings - 1
    error('Ring difference cannot be larger than the number of rings - 1')
end

% Read the components
if isempty(comp)
    d = dir(f_path);
    if d.bytes < 13636640
        error('Unrecognized normalization file layout')
    end
    off = [0 13312000 13355680 13611040];
else
    dimsC = {[NR 6400], [21 NR], [NS R], 6400};
    for kk = 1 : 4
        if numel(comp(kk).dims) < 1 || prod(comp(kk).dims) ~= prod(dimsC{kk}) || (kk < 4 && ~isequal(comp(kk).dims, dimsC{kk}))
            error('Unexpected size of the Vision normalization component %d', kk)
        end
    end
    off = [comp(1:4).offset];
end
fid = fopen(f_path, 'r');
if fid < 0
    error('Could not open the normalization file %s', f_path)
end
G = readComponent(fid, off(1), [NR 6400]);
I = readComponent(fid, off(2), [21 NR]);
E = readComponent(fid, off(3), [NS R]);
A = readComponent(fid, off(4), [6400 1]);
fclose(fid);

% Siemens span-1 planes: blocks of [80, 79, 79, 78, 78, ..., 1, 1] planes
sizes = [R; reshape(repmat((R-1:-1:1), 2, 1), [], 1)];
starts = [0; cumsum(sizes(1:end-1))];
[RA, RB] = ndgrid(0:R-1, 0:R-1);
dd = abs(RA - RB);
blk = max(2 * dd - 1, 0);
Pmat = starts(blk + 1) + min(RA, RB);
% The blocks 2d-1 and 2d should be identical
dmax = 0;
for d = 1 : R - 1
    n = R - d;
    ia = starts(2*d) + (1:n);
    ib = starts(2*d + 1) + (1:n);
    dmax = max([dmax; abs(A(ia) - A(ib)) / max(abs(A(ia))); max(abs(G(:, ia) - G(:, ib)), [], 1)' / max(max(abs(G(:, ia))))]);
end
if dmax > 1e-6
    warning('The axial or geometric effects of the plane blocks 2d-1 and 2d differ by up to %g (relative), only the block 2d-1 is used', dmax)
end

% Transaxial LOR table: all detector pairs x1 < x2, binned as in saveSinogram
[x2, x1] = find(tril(true(W), -1));
x2 = x2 - 1;
x1 = x1 - 1;
half = W / 2;
j = floor(mod(x2 + x1 + half + 1, W) / 2);
i = abs(x2 - x1 - half);
i((x1 < j) | (j + half < x2)) = -i((x1 < j) | (j + half < x2));
swap = (j * 2) < -i | i <= (j - half) * 2;
if mod(Nd, 2) == 0
    if nSide == 0
        nSide = 1;
    end
    acc = i <= (Nd / 2 + min(0, nSide)) & i >= (-Nd / 2 + max(0, nSide));
    ib = i + Nd / 2 - max(0, nSide);
else
    acc = i <= floor(Nd / 2) & i >= -floor(Nd / 2);
    ib = i + floor(Nd / 2);
end
jb = floor(j / floor(half / Nang));
if any(acc & (jb >= Nang))
    error('Nang is not compatible with the number of detectors per ring')
end
% Siemens crystals of the OMEGA detectors
x = (0:W-1)';
if W == NS
    sx = mod(x + 1, NS);
else
    sx = x + floor(x / options.cryst_per_block) + 1;
end
s1 = sx(x1 + 1);
s2 = sx(x2 + 1);
% Transaxial LOR length factor D^(-n)
if nLen ~= 0
    hLen = lorLengthFactor(options, W, x1, x2, nLen);
else
    hLen = ones(size(x1));
end
% Pseudo detectors (gap crystals) have no efficiency
acc = acc & mod(s1, 21) ~= 0 & mod(s2, 21) ~= 0;
swap = swap(acc);
s1 = s1(acc);
s2 = s2(acc);
hLen = hLen(acc);
bin = ib(acc) + jb(acc) * Nd + 1;
% Orientation-corrected Siemens bin
[t, v] = siemensBin(s1, s2);
r = t + 260;
p = mod(v, 21);
oor = r < 0 | r > NR - 1;
rc = min(max(r, 0), NR - 1);
Ifac = ones(size(r));
Ifac(~oor) = I(sub2ind([21 NR], p(~oor) + 1, r(~oor) + 1));
% LORs with the swap flag first
[~, ord] = sort(double(~swap));
s1 = s1(ord);
s2 = s2(ord);
bin = bin(ord);
rc = rc(ord);
Ifac = Ifac(ord) .* hLen(ord);
oor = oor(ord);
nS = sum(swap);
nTr = numel(bin);
oorIdx = find(oor);
bS = bin(1:nS);
bN = bin(nS+1:end);
NdNa = Nd * Nang;
clear x1 x2 i j ib jb t v r p ord acc

E1 = E(s1 + 1, :);
E2 = E(s2 + 1, :);
Winv = 1 ./ (G .* reshape(A, 1, []));
Winv(~isfinite(Winv)) = 0;
ri = rc + 1;

% Sinogram planes of the ring pairs: the detector 1 is at x1 and the detector 2 at x2
% (ring numbers are swapped for the swapped LORs and then again since x2 > x1)
okPair = dd <= rd;
% (the ring pairs beyond the ring difference are not binned and are set to a valid pair)
planeS = sinogramPlane(RA .* okPair, RB .* okPair, span, rings, seg);
planeN = sinogramPlane(RB .* okPair, RA .* okPair, span, rings, seg);
if any(planeS(okPair) < 0 | planeS(okPair) >= totSinos | planeN(okPair) < 0 | planeN(okPair) >= totSinos)
    error('The segment table does not match the number of rings and the ring difference')
end
nRP = sum(okPair(:));

S = zeros(NdNa, totSinos, 'single');
oorSum = 0;
tic;
for ra = 0 : R - 1
    E1a = E1(:, ra + 1);
    for rb = 0 : R - 1
        if ~okPair(ra + 1, rb + 1)
            continue
        end
        w = Winv(:, Pmat(ra + 1, rb + 1) + 1);
        e = E1a .* E2(:, rb + 1) .* Ifac .* w(ri);
        if nS > 0
            pl = planeS(ra + 1, rb + 1) + 1;
            S(:, pl) = S(:, pl) + single(accumarray(bS, e(1:nS), [NdNa 1]));
        end
        if nS < nTr
            pl = planeN(ra + 1, rb + 1) + 1;
            S(:, pl) = S(:, pl) + single(accumarray(bN, e(nS+1:end), [NdNa 1]));
        end
        if ~isempty(oorIdx)
            oorSum = oorSum + sum(e(oorIdx));
        end
    end
end
if isfield(options, 'verbose') && options.verbose
    fprintf('Vision normalization summed in %.1f s\n', toc)
end

% Mean efficiency of the LORs is one
totEff = sum(S(:), 'double');
if totEff <= 0
    error('The Vision normalization has no non-zero efficiencies')
end
nLOR = nTr * nRP;
effSum = S ./ single(totEff / nLOR);
% Mean efficiency with the length factor referred to the effective ring
% diameter, (Dref / D)^n instead of D^(-n) (a missing DOI is zero as in
% detector_coordinates)
meanEff = totEff / nLOR;
if nLen ~= 0
    DOI = 0;
    if isfield(options, 'DOI') && ~isempty(options.DOI)
        DOI = double(options.DOI(1));
    end
    meanEff = meanEff * (double(options.diameter(1)) + 2 * DOI)^nLen;
end
if ~isempty(oorIdx)
    warning(['%d LORs (%.3g %% of all LORs) are outside of the %d radial bins of the Vision normalization; ' ...
        'their crystal interference is set to one and the geometric effect of the edge bin is used. They carry %.3g %% of the total efficiency.'], ...
        numel(oorIdx) * nRP, 100 * numel(oorIdx) / nTr, NR, 100 * oorSum / totEff)
end
coef = zeros(size(effSum), 'single');
nz = effSum > 0;
coef(nz) = 1 ./ effSum(nz);
coef = reshape(coef, Nd, Nang, totSinos);
effSum = reshape(effSum, Nd, Nang, totSinos);
end
