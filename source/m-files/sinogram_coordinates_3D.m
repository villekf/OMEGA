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
    z_length = double(options.rings + 1+ numel(options.ringGaps)) * cr_pz;
else
    z_length = double(options.rings + 1) * cr_pz;
end
% if options.nLayers > 1 && options.cryst_per_block_axial(1) ~= options.cryst_per_block_axial(2)
%     maxZ = z_length + cr_pz * (options.linear_multip - 1);
% else
    maxZ = z_length;
% end

if options.nLayers > 1 && options.cryst_per_block_axial(1) ~= options.cryst_per_block_axial(2)
    apu = zeros(options.linear_multip,1);
    for kk = 0 : options.linear_multip - 1
        apu(kk + 1) = kk * cr_pz;
    end
    apu = repelem(apu, sum(options.cryst_per_block_axial));
end
%% Compute the 3D coordinates

if options.span > 1
    % Ring coordinates
    if isfield(options, 'ringGaps') && numel(options.ringGaps) > 0 && sum(options.ringGaps) > 0
        z = linspace(cr_pz, z_length + cr_pz, options.rings + 2 + numel(options.ringGaps))';
    else
        z = linspace(cr_pz, z_length + cr_pz, options.rings + 2)';
    end
    %     if min(z(:)) == 0
    %         z = z + (options.axial_fov - options.rings * options.cr_pz)/2 + options.cr_pz/2;
    %     end
    if isfield(options, 'ringGaps') && numel(options.ringGaps) > 0 && sum(options.ringGaps) > 0
        z = z(1 : options.rings + numel(options.ringGaps));
        z(options.ringGaps + 1) = [];
    else
        z = z(1 : options.rings);
    end
    if options.nLayers > 1 && options.cryst_per_block_axial(1) ~= options.cryst_per_block_axial(2)
        z = z + apu;
    end
    z = z - maxZ / 2;
    ringsp = options.rings;
    z_ring = zeros(options.rings, options.rings, 2);
    % Create ring combinations
    z_ring(:,:,1) = repmat(z, 1, options.rings);
    z_ring(:,:,2) = repmat(z', options.rings, 1);
    z_ring = reshape(z_ring, options.rings * options.rings, 2);
    kkj = zeros(floor((options.ring_difference-ceil(options.span/2))/options.span) + 1, 1);
    for kk = 1 : floor((options.ring_difference-ceil(options.span/2))/options.span) + 1
        kkj(kk) = ceil(options.span/2) + options.span*(kk - 1);
    end
    offset2 = cumsum(options.segment_table);

    % Perpendicular rings
    z = zeros(options.TotSinos,2);
    z(1:2:Nz,1) = z_ring(1:ringsp+1:ringsp^2,1);
    z(1:2:Nz,2) = z_ring(1:ringsp+1:ringsp^2,2);
    mean_jh = zeros(options.TotSinos,1);
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
            apu = z_ring((kkj(ih)+jh-1)*options.rings+1:options.rings+1:end, 1);
            z(offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1, 1) = z(offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1, 1) + (apu);
            apu2 = z_ring(kkj(ih)+jh:options.rings+1:(options.rings-kkj(ih)-jh+1)*options.rings, 1);
            z(offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1, 1) = z(offset2(2*ih)+jh:2:offset2(2*ih + 1)-jh+1, 1) + (apu2);

            apu = z_ring((kkj(ih)+jh-1)*options.rings+1:options.rings+1:end, 2);
            z(offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1, 2) = z(offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1, 2) + (apu);
            apu2 = z_ring(kkj(ih)+jh:options.rings+1:(options.rings-kkj(ih)-jh+1)*options.rings, 2);
            z(offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1, 2) = z(offset2(2*ih)+jh:2:offset2(2*ih + 1)-jh+1, 2) + (apu2);

            loc = ismember(1:options.TotSinos,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1);
            mean_jh(loc) = mean_jh(loc) + 1;
            loc = ismember(1:options.TotSinos,offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1);
            mean_jh(loc) = mean_jh(loc) + 1;
        end
    end
    mean_jh(mean_jh == 0) = 1;
    % Take the mean value of coordinates
    z(:,1) = z(:,1) ./ mean_jh;
    z(:,2) = z(:,2) ./ mean_jh;
    ind1 = z(:,1)>z(:,2);
    z(z(:,1)<z(:,2),:) = fliplr(z(z(:,1)<z(:,2),:));
    z(ind1,:) = fliplr(z(ind1,:));

else

    dif = cr_pz;

    if options.nLayers > 1 && options.cryst_per_block_axial(1) ~= options.cryst_per_block_axial(2)
        gap = ((double(options.rings) * cr_pz) - (options.cryst_per_block_axial(1) * options.linear_multip)*cr_pz) / options.linear_multip / 2;
        gap = gap : gap * 2 : gap * options.linear_multip * 2;
    end
    % Michelogram row/column indices for each segment
    %     p = zeros(floor((ring_difference-ceil(span/2))/span) + 2,1);
    %     for kk = 0 : floor((ring_difference-ceil(span/2))/span)
    %         p(kk+2) = ceil(span/2)*2 + (span*2)*kk;
    %     end

    if isfield(options, 'ringGaps') && numel(options.ringGaps) > 0 && sum(options.ringGaps) > 0
        z = zeros((options.rings + numel(options.ringGaps))^2,2);
        z2 = zeros((options.rings + numel(options.ringGaps))^2,2);
        loppu = options.rings + numel(options.ringGaps);
    else
        z = zeros(options.rings^2,2);
        z2 = zeros(options.rings^2,2);
        loppu = options.rings;
    end
    % if options.cryst_per_block_axial(1) > options.cryst_per_block_axial(end) && layer1 == 2
    %     ind1(options.cryst_per_block_axial(1) : options.cryst_per_block_axial(1) : end) = inf;
    % elseif options.cryst_per_block_axial(end) > options.cryst_per_block_axial(1) && layer1 == 1
    %     ind1(options.cryst_per_block_axial(end) : options.cryst_per_block_axial(end) : end) = inf;
    % end
    ind1 = ones(loppu,1);
    if options.nLayers > 1
        % apu = ind2;
        % ind3 = ind2;
        if options.cryst_per_block_axial(1) > options.cryst_per_block_axial(end) && layer2 == 1
            r = options.cryst_per_block_axial(end) * options.linear_multip;
            apu = repelem(gap, options.cryst_per_block_axial(end))';
            ind2 = (0 : cr_pz : cr_pz * (r - 1))' + apu;
            insert_indices = setdiff(1:length(ind2), options.cryst_per_block(1):options.cryst_per_block(1):length(ind2));
            ind3 = Inf(loppu,1);
            ind3(insert_indices) = ind2;
            ind2 = ind3;
        elseif options.cryst_per_block_axial(end) > options.cryst_per_block_axial(1) && layer2 == 1
            r = options.cryst_per_block_axial(1) * options.linear_multip;
            apu = repelem(gap, options.cryst_per_block_axial(1))';
            ind2 = (0 : cr_pz : cr_pz * (r - 1))' + apu;
            ind3 = Inf(loppu,1);
            insert_indices = setdiff(1:length(ind3), options.cryst_per_block(end):options.cryst_per_block(end):length(ind3));
            ind3(insert_indices) = ind2;
            ind2 = ind3;
        else
            ind2 = (0 : cr_pz : cr_pz * (loppu - 1))';
        end
        if ((options.cryst_per_block_axial(1) > options.cryst_per_block_axial(end)) || (options.cryst_per_block_axial(end) > options.cryst_per_block_axial(1))) && layer1 == 1
            apu1 = Inf(loppu,1);
            if options.cryst_per_block_axial(end) > options.cryst_per_block_axial(1)
                apu = repelem(gap, options.cryst_per_block_axial(1))';
                insert_indices = setdiff(1:length(ind1), options.cryst_per_block(end):options.cryst_per_block(end):length(ind1));
            else
                apu = repelem(gap, options.cryst_per_block_axial(end))';
                insert_indices = setdiff(1:length(ind1), options.cryst_per_block(1):options.cryst_per_block(1):length(ind1));
            end
            apu1(insert_indices) = apu;
            % ind1 = ind1 + apu;
        end
        % if options.cryst_per_block_axial(1) > options.cryst_per_block_axial(end) && layer1 == 2
        %     ind3(options.cryst_per_block_axial(1) : options.cryst_per_block_axial(1) : end) = inf;
        % elseif options.cryst_per_block_axial(end) > options.cryst_per_block_axial(1) && layer1 == 1
        %     ind3(options.cryst_per_block_axial(end) : options.cryst_per_block_axial(end) : end) = inf;
        % end
        % apu = apu(isinf(ind3));
    else
        ind2 = (0 : cr_pz : cr_pz * (loppu - 1))';
    end
    % if options.nLayers > 1 && options.cryst_per_block_axial(1) ~= options.cryst_per_block_axial(2)
    %     ind = ind + apu;
    %     layers = repmat([0;1], options.cryst_per_block_axial(2), 1);
    %     layers = repmat([layers;0], options.linear_multip, 1);
    % end
    uu = 0;
    for t=1:loppu
        % if options.nLayers > 1 && options.cryst_per_block_axial(1) ~= options.cryst_per_block_axial(2)
        %     z((t-1) * loppu + 1 : t * loppu,:) = [ dif * (t-1) * ones(loppu,1) + apu(t) ind];
        %     z2((t-1) * loppu + 1 : t * loppu,:) = [ layers(t) * ones(loppu,1) layers];
        % else
        if options.nLayers > 1 && options.cryst_per_block_axial(1) ~= options.cryst_per_block_axial(2) && layer1 == 1
            z((t-1) * loppu + 1 : t * loppu,:) = [ dif * (uu) * ind1 + apu1(t) ind2];
            uu = uu + 1;
        else
            z((t-1) * loppu + 1 : t * loppu,:) = [ dif * (t-1) * ind1 ind2];
        end
            if options.nLayers > 1 && ((options.cryst_per_block_axial(1) > options.cryst_per_block_axial(end) && layer1 == 1) || (options.cryst_per_block_axial(end) > options.cryst_per_block_axial(1) && layer1 == 1))
                % KalPa = dif * (t-1);
                if isinf(apu1(t))
                    z((t-1) * loppu + 1 : t * loppu,:) = [ inf * ind1 ind2];
                    uu = uu - 1;
                end
            end
        % end
    end
    z(isnan(z)) = inf;
    if isfield(options, 'ringGaps') && numel(options.ringGaps) > 0 && sum(options.ringGaps) > 0
        ind = repmat(options.ringGaps,loppu,1) + repelem(0:loppu - 1,3)' * loppu + 1;
        ind = ind + repmat((0:numel(options.ringGaps)-1)',numel(ind)/numel(options.ringGaps),1);
        z(ind,:) = [];
        for kk = numel(options.ringGaps) : - 1 : 1
            z(options.rings * (options.ringGaps(kk) + kk - 1) + 1 : options.rings * (options.ringGaps(kk) + 1 + kk - 1),:) = [];
        end
    end
    z = z - (maxZ/2 - cr_pz);
    % if options.nLayers > 1 && options.cryst_per_block_axial(1) ~= options.cryst_per_block_axial(2)
    %     apuZ = zeros(options.rings^2,1);
    %     apuZ(ismember(z2,[1 1], 'rows')) = 3;
    %     apuZ(ismember(z2,[1 0], 'rows')) = 1;
    %     apuZ(ismember(z2,[0 1], 'rows')) = 2;
    %     z = [apuZ z];
    % end
end
