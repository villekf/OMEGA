function [index, pituus, varargout] = subset_angles(options, varargin)
%SUBSET_ANGLES Form the subset indices for the subset_type 5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi
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

if options.use_raw_data
    [x,y] = detector_coordinates(options);
    % Determine the rings of each detector
    apu = idivide(varargin{1}, uint16(options.det_per_ring),'ceil');
    % Same ring
    ind = apu(:,1) == apu(:,2);
    indi = find(ind);
    % Make all detector numbers constrained in [1, det_per_ring]
    varargin{1} = mod(varargin{1}, options.det_per_ring);
    varargin{1}(varargin{1} == 0) = options.det_per_ring;
    x2 = single(x);
    y2 = single(y);
    % Coordinates of detectors that are on the same rings (2D)
    x = single(x(varargin{1}(ind,:)));
    y = single(y(varargin{1}(ind,:)));
else
    [x,y] = sinogram_coordinates_2D(options);
    sino_length = options.Nang*options.Ndist;
end
% Lengths of the triangle sides
jokux = abs(x(:,1)-x(:,2));
jokuy = abs(y(:,1)-y(:,2));
% Indices on the other side of the circle (angle over 90 degrees)
ind = x(:,1) > x(:,2) & y(:,1) > y(:,2) | x(:,2) > x(:,1) & y(:,2) > y(:,1);
% Angle
kulma = atand(jokux./jokuy);
kulma(ind) = kulma(ind) + 90;
if nargout >= 3
    varargout{1} = kulma;
end
% Combine similar angles together
kulma = round(kulma);
kulma(kulma==180) = 179;
[kulma, I] = sort(kulma);
n_angles = options.n_angles;
maara = 180/n_angles;

if options.use_raw_data
    maara2 = 360/n_angles;
    index1 = cell(maara + maara2 * (options.rings-1),1);
    pituus = zeros(maara + maara2 * (options.rings-1),1,'uint32');
    % 2D case
    for ii = 1 : maara
        % Take the indices of angles that are between the defined number of
        % angles
        index1{ii} = uint32(indi(I(kulma >= (ii-1)*n_angles & kulma < ii*n_angles)));
        pituus(ii) = length(index1{ii});
    end
else
    maara2 = maara * 2;
    index1 = cell(maara,1);
    pituus = zeros(maara*floor(length(options.segment_table)/2)*2 + maara,1,'uint32');
    index0 = cell(maara,1);
    % 2D case
    for ii = 1 : maara
        % Take the indices of angles that are between the defined number of
        % angles
        index0{ii} = uint32(I(kulma >= (ii-1)*n_angles & kulma < ii*n_angles));
        index1{ii} = repmat(index0{ii},options.segment_table(1),1) + uint32(repelem((0:options.segment_table(1)-1)*sino_length,length(index0{ii}))');
        pituus(ii) = length(index1{ii});
    end
    index1 = cell2mat(index1);
end

% 3D case
if options.use_raw_data
    ll = 1;
    ring_diff = zeros(options.rings,1);
    % Constant increment scheme
    for kk = 1 : options.rings
        ring_diff(kk+1) = mod(ring_diff(kk) + round(0.7*options.rings), options.rings + 1);
        while any(ismember(ring_diff(kk+1), ring_diff(1:kk)))
            ring_diff(kk+1) = ring_diff(kk+1) + 1;
        end
    end
    ring_diff(1) = [];
    for kk = ring_diff
        % Choose the detectors that have an angle of most 180 degrees, and
        % angles that are betweem 180 and 360.
        ind1 = apu(:,1) == (apu(:,2) + kk);
        ind2 = apu(:,2) == (apu(:,1) + kk);
        ind = logical(ind1 + ind2);
        indi = find(ind);
        xx = x2(varargin{1}(ind,:));
        yy = y2(varargin{1}(ind,:));
        jokux = abs(xx(:,1)-xx(:,2));
        jokuy = abs(yy(:,1)-yy(:,2));
        ind = xx(:,1) > xx(:,2) & yy(:,1) > yy(:,2) | xx(:,2) > xx(:,1) & yy(:,2) > yy(:,1);
        ind2 = ind2(ind);
        kulma = atand(jokux./jokuy);
        kulma(ind) = kulma(ind) + 90;
        kulma(ind2) = kulma(ind2) + 180;
        kulma = round(kulma);
        [kulma, I] = sort(kulma);
        for ii = 1 : maara2
            index1{maara + ll} = uint32(indi(I(kulma >= (ii-1)*n_angles & kulma < ii*n_angles)));
            pituus(ii) = length(index1{ii});
            ll = ll + 1;
        end
    end
    index = cell2mat(index1);
else
    ring_diff = zeros(ceil(length(options.segment_table)/2),1);
    index2 = cell(maara*2*floor(length(options.segment_table)/2),1);
    for kk = 1 : floor(length(options.segment_table)/2)
        ring_diff(kk+1) = mod(ring_diff(kk) + round(0.7*length(options.segment_table)/2), floor(length(options.segment_table)/2) + 1);
        while any(ismember(ring_diff(kk+1), ring_diff(1:kk)))
            ring_diff(kk+1) = ring_diff(kk+1) + 1;
        end
        ll = 1;
        alku = sum(options.segment_table(1:ring_diff(kk+1)*2-1)) * sino_length;
        kerroin = 0;
        for ii = 1 : maara2
            index2{maara2*(kk-1) + ii} = repmat(index0{ll},options.segment_table(2+(ring_diff(kk+1)-1)*2),1) + ...
                uint32(repelem((options.segment_table(2+(ring_diff(kk+1)-1)*2)*kerroin:options.segment_table(2+(ring_diff(kk+1)-1)*2)...
                *(kerroin+1)-1)*sino_length + alku,length(index0{ll}))');
            pituus(maara2*(kk-1) + ii + maara) = length(index2{maara2*(kk-1) + ii});
            % If the "other" segment (segment tables have two 3D components
            % of identical sizes) is reached, the angle is increased from
            % 180 to 360 
            if mod(ll,maara) == 0
                ll = 1;
                kerroin = 1;
            else
                ll = ll + 1;
            end
        end
    end
    index = [index1;cell2mat(index2)];
end