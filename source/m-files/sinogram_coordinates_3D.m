function z = sinogram_coordinates_3D(options, varargin)
%% Coordinates for the sinogram detectors
% This code is used to compute the z-coordinates for the detector
% coordinates in sinogram space.
%
% OUTPUT:
%   z = Z-coordinates for the sinograms
%
% See also sinogram_coordinates_2D, detector_coordinates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2022 Samuli Summala, Ville-Veikko Wettenhovi
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

if ~isempty(varargin) && ~isempty(varargin{1})
    layer1 = varargin{1}(1);
    layer2 = varargin{1}(2);
else
    layer1 = 1;
    layer2 = 1;
end
cr_pz = options.cr_pz;
Nz = options.rings*2-1;
% span = options.span;
% ring_difference = options.ring_difference;

if isfield(options, 'ringGaps') && numel(options.ringGaps) > 0 && sum(options.ringGaps) > 0
    % z_length = double(options.rings + 1+ numel(options.ringGaps)) * cr_pz;
    z_length = double(options.rings + 1) * cr_pz + sum(options.ringGaps);
else
    z_length = double(options.rings + 1) * cr_pz;
end
% if options.nLayers > 1 && options.cryst_per_block_axial(1) ~= options.cryst_per_block_axial(2)
%     maxZ = z_length + cr_pz * (options.linear_multip - 1);
% else
    maxZ = z_length;
% end

if options.nLayers > 1 && options.cryst_per_block_axial(1) ~= options.cryst_per_block_axial(2)
    % Layers with a different number of axial crystals: each layer has its own ring positions
    z = unequalLayerCoordinates3D(options, layer1, layer2, maxZ);
    return
end
%% Compute the 3D coordinates

if options.span > 1
    % Ring coordinates
    % The ring step is exactly cr_pz, the gaps between the axial blocks are added below
    z = linspace(cr_pz, double(options.rings + 1) * cr_pz + cr_pz, options.rings + 2)';
    z = z(1 : options.rings);
    if isfield(options, 'ringGaps') && numel(options.ringGaps) > 0 && sum(options.ringGaps) > 0
        gaps = cumsum(options.ringGaps);
        cpba = options.cryst_per_block_axial(1);
        for kk = 1 : options.rings / cpba - 1
            z(cpba*kk + 1: cpba * (kk + 1)) = z(cpba*kk + 1: cpba * (kk + 1)) + gaps(kk);
        end
    end
    z = z - maxZ / 2;
    ringsp = options.rings;
    z_ring = zeros(options.rings, options.rings, 2);
    % Create ring combinations
    z_ring(:,:,1) = repmat(z, 1, options.rings);
    z_ring(:,:,2) = repmat(z', options.rings, 1);
    z_ring = reshape(z_ring, options.rings * options.rings, 2);
    % One segment (a pair of the negative and positive ring differences) for every entry pair of the segment table,
    % i.e. the planes follow the table of the loader also when it extends beyond the ring difference
    kkj = ceil(options.span/2) + options.span * (0 : floor(length(options.segment_table)/2) - 1)';
    offset2 = cumsum(options.segment_table);

    % Perpendicular rings
    totSinos = planesPerBlock(options);
    z = zeros(totSinos,2);
    z(1:2:Nz,1) = z_ring(1:ringsp+1:ringsp^2,1);
    z(1:2:Nz,2) = z_ring(1:ringsp+1:ringsp^2,2);
    mean_jh = zeros(totSinos,1);
    mean_jh(1:2:Nz) = 1;
    % Then the detectors on adjacent rings
    for jh=1:floor(options.span/2)
        apu = z_ring(jh*ringsp+1:ringsp+1:ringsp^2, 1);
        apu2 = z_ring(jh+1:ringsp+1:(ringsp-jh)*ringsp, 1);
        z(jh+1:2:offset2(1)-jh,1) = z(jh+1:2:offset2(1)-jh, 1) + apu + apu2;
        apu = z_ring(jh*ringsp+1:ringsp+1:ringsp^2, 2);
        apu2 = z_ring(jh+1:ringsp+1:(ringsp-jh)*ringsp, 2);
        z(jh+1:2:offset2(1)-jh,2) = z(jh+1:2:offset2(1)-jh, 2) + apu + apu2;
        loc = ismember(1:Nz,jh+1:2:offset2(1)-jh);
        mean_jh(loc) = mean_jh(loc) + 2;
    end
    % Lastly the rest of the detectors with the amount of combined LORs
    % specified with the span value
    for ih=1:floor(length(options.segment_table)/2)
        for jh=1:options.span
            % Ring differences larger than ring_difference are not stored by the loader
            if kkj(ih) + jh - 1 > options.ring_difference
                continue
            end
            apu = z_ring((kkj(ih)+jh-1)*options.rings+1:options.rings+1:end, 1);
            z(offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1, 1) = z(offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1, 1) + (apu);
            apu2 = z_ring(kkj(ih)+jh:options.rings+1:(options.rings-kkj(ih)-jh+1)*options.rings, 1);
            z(offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1, 1) = z(offset2(2*ih)+jh:2:offset2(2*ih + 1)-jh+1, 1) + (apu2);

            apu = z_ring((kkj(ih)+jh-1)*options.rings+1:options.rings+1:end, 2);
            z(offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1, 2) = z(offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1, 2) + (apu);
            apu2 = z_ring(kkj(ih)+jh:options.rings+1:(options.rings-kkj(ih)-jh+1)*options.rings, 2);
            z(offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1, 2) = z(offset2(2*ih)+jh:2:offset2(2*ih + 1)-jh+1, 2) + (apu2);

            loc = ismember(1:totSinos,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1);
            mean_jh(loc) = mean_jh(loc) + 1;
            loc = ismember(1:totSinos,offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1);
            mean_jh(loc) = mean_jh(loc) + 1;
        end
    end
    % Planes without any ring pairs (segments beyond the ring difference)
    noPairs = mean_jh == 0;
    mean_jh(noPairs) = 1;
    % Take the mean value of coordinates
    z(:,1) = z(:,1) ./ mean_jh;
    z(:,2) = z(:,2) ./ mean_jh;
    ind1 = z(:,1)>z(:,2);
    z(z(:,1)<z(:,2),:) = fliplr(z(z(:,1)<z(:,2),:));
    z(ind1,:) = fliplr(z(ind1,:));
    z(noPairs,:) = Inf;

else

    dif = cr_pz;

    z = zeros(options.rings^2,2);
    loppu = options.rings;
    ind1 = ones(loppu,1);
    ind2 = (0 : cr_pz : cr_pz * (loppu - 1))';
    if isfield(options, 'ringGaps') && sum(options.ringGaps) > 0
        gaps = repelem([0;cumsum(options.ringGaps(:))], options.cryst_per_block_axial(1), 1);
    end
    for t=1:loppu
        z((t-1) * loppu + 1 : t * loppu,:) = [ dif * (t-1) * ind1 ind2];
        if isfield(options, 'ringGaps') && sum(options.ringGaps) > 0
            z((t-1) * loppu + 1 : t * loppu,2) = z((t-1) * loppu + 1 : t * loppu,2) + gaps(1:loppu);
            z((t-1) * loppu + 1 : t * loppu,1) = z((t-1) * loppu + 1 : t * loppu,1) + gaps(t);
        end
    end
    z(isnan(z)) = inf;
    z = z - (maxZ/2 - cr_pz);
end
end

function n = planesPerBlock(options)
% Number of sinogram planes of one layer combination (span > 1). With multiple layers options.TotSinos may include all
% nLayers^2 layer combinations, the number of planes of one combination is then the sum of the segment table.
if options.nLayers > 1
    n = double(sum(options.segment_table));
else
    n = double(options.TotSinos);
end
end

function zr = unequalLayerRings(options, layer, maxZ)
% Axial ring positions of one layer of a multi-layer scanner whose layers have a different number of axial crystals
% (cryst_per_block_axial). The ring slots of every axial block are as many as in the layer with the most crystals, the
% layers with fewer crystals are centered in the block and have a pseudo ring (Inf) at the end of the block. This is
% the same ring numbering as the loader uses (the module number is added to the ring number of the layers with fewer
% axial crystals).
cr_pz = options.cr_pz;
rings = double(options.rings);
cpba = double(options.cryst_per_block_axial(:));
M = max(cpba);
n = cpba(layer);
s = (0 : rings - 1)';
zr = (s + (M - n) / 2) * cr_pz;
if isfield(options, 'ringGaps') && numel(options.ringGaps) > 0 && sum(options.ringGaps) > 0
    % Gaps between the axial blocks
    gaps = repelem([0;cumsum(options.ringGaps(:))], M, 1);
    zr = zr + gaps(1 : rings);
end
zr(mod(s, M) >= n) = Inf;
zr = zr - (maxZ / 2 - cr_pz);
end

function z = unequalLayerCoordinates3D(options, layer1, layer2, maxZ)
% Axial sinogram coordinates of a multi-layer scanner whose layers have a different number of axial crystals. Every
% sinogram plane is the mean of the ring pairs (existing rings of the two layers only) that the loader places into the
% plane: the first ring (first column) is the ring of the first detector of the sinogram and the second ring that of
% the second detector. The planes that contain no ring pairs are Inf.
rings = double(options.rings);
zA = unequalLayerRings(options, layer1, maxZ);
zB = unequalLayerRings(options, layer2, maxZ);
[b, a] = meshgrid(0 : rings - 1, 0 : rings - 1);
a = a(:);
b = b(:);
if options.span > 1
    span = options.span;
    half = floor(span / 2);
    seg = cumsum(options.segment_table(:));
    nPlanes = planesPerBlock(options);
    % The same plane index as in the loader (saveSinogram.h)
    erotus = b - a;
    summa = a + b;
    ae = abs(erotus);
    k = floor((ae + half) / span);
    kk = max(k - 1, 0);
    valid = ae <= options.ring_difference & 2 * kk + 1 < numel(seg);
    ind1 = min(2 * kk + 1, numel(seg));
    ind2 = min(2 * kk + 2, numel(seg));
    off = seg(ind1);
    off(erotus >= 0) = seg(ind2(erotus >= 0));
    P = summa - (half * (2 * k - 1) + k) + off;
    P(ae <= half) = summa(ae <= half);
else
    nPlanes = rings^2;
    P = a * rings + b;
    valid = true(size(P));
end
valid = valid & isfinite(zA(a + 1)) & isfinite(zB(b + 1)) & P >= 0 & P < nPlanes;
P = P(valid) + 1;
count = accumarray(P, 1, [nPlanes 1]);
populated = count > 0;
z = Inf(nPlanes, 2);
sA = accumarray(P, zA(a(valid) + 1), [nPlanes 1]);
sB = accumarray(P, zB(b(valid) + 1), [nPlanes 1]);
z(populated, 1) = sA(populated) ./ count(populated);
z(populated, 2) = sB(populated) ./ count(populated);
end
