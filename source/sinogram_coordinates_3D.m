function z = sinogram_coordinates_3D(options)
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

cr_pz = options.cr_pz;
Nz = options.rings*2-1;
span = options.span;
ring_difference = options.ring_difference;

%% Compute the 3D coordinates

if options.span > 1
    % Ring coordinates
    z_length = double(options.rings + 1) * options.cr_pz;
    z = linspace(0, z_length, options.rings + 2)';
    if min(z(:)) == 0
        z = z + (options.axial_fov - options.rings * options.cr_pz)/2 + options.cr_pz/2;
    end
    z = z(1 : options.rings);
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
    mean_jh(1:2:options.Nz) = 1;
    % Then the detectors on adjacent rings
    for jh=1:floor(options.span/2)
        apu = z_ring(jh*ringsp+1:ringsp+1:ringsp^2, 1);
        apu2 = z_ring(jh+1:ringsp+1:(ringsp-jh)*ringsp, 1);
        z(jh+1:2:offset2(1)-jh,1) = z(jh+1:2:offset2(1)-jh, 1) + apu + apu2;
        apu = z_ring(jh*ringsp+1:ringsp+1:ringsp^2, 2);
        apu2 = z_ring(jh+1:ringsp+1:(ringsp-jh)*ringsp, 2);
        z(jh+1:2:offset2(1)-jh,2) = z(jh+1:2:offset2(1)-jh, 2) + apu + apu2;
        loc = ismember(1:options.Nz,jh+1:2:offset2(1)-jh);
        mean_jh(loc) = mean_jh(loc) + 2;
    end
    % Lastly the rest of the detectors with the amount of combined LORs
    % specified with the span value
    for ih=1:floor(length(options.segment_table)/2)
        for jh=1:options.span
            apu = z_ring((kkj(ih)+jh-1)*options.rings+1:options.rings+1:end, 1);
            loc = ismember(1:options.TotSinos,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1);
            mean_jh(loc) = mean_jh(loc) + 1;
            z(offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1, 1) = z(offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1, 1) + (apu);
            apu2 = z_ring(kkj(ih)+jh:options.rings+1:(options.rings-kkj(ih)-jh+1)*options.rings, 1);
            loc = ismember(1:options.TotSinos,offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1);
            mean_jh(loc) = mean_jh(loc) + 1;
            z(offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1, 1) = z(offset2(2*ih)+jh:2:offset2(2*ih + 1)-jh+1, 1) + (apu2);
            apu = z_ring((kkj(ih)+jh-1)*options.rings+1:options.rings+1:end, 2);
            z(offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1, 2) = z(offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1, 2) + (apu);
            apu2 = z_ring(kkj(ih)+jh:options.rings+1:(options.rings-kkj(ih)-jh+1)*options.rings, 2);
            z(offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1, 2) = z(offset2(2*ih)+jh:2:offset2(2*ih + 1)-jh+1, 2) + (apu2);
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
    % Michelogram row/column indices for each segment
    p = zeros(floor((ring_difference-ceil(span/2))/span) + 2,1);
    for kk = 0 : floor((ring_difference-ceil(span/2))/span)
        p(kk+2) = ceil(span/2)*2 + (span*2)*kk;
    end
    
    z = zeros(options.rings^2,2);
    
    for t=1:options.rings
        z((t-1) * options.rings + 1 : t * options.rings,:) = [ dif * (t-1) * ones(options.rings,1) (0 : dif : dif * (options.rings - 1))'];
    end
end
