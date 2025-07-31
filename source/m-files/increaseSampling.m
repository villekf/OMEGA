function [x, y, options] = increaseSampling(options, x, y, interpolateSinogram)
%INCREASESAMPLING Increase the sampling rate
%   This function increases the sampling rate of the input detectors and
%   sinogram/raw data vector. Can be used to prevent aliasing artifacts
%   when using too high image resolution. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ville-Veikko Wettenhovi
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Raw list-mode data
if options.use_raw_data
    
    % Interpolate the coordinates
    x = interp1(1:options.sampling_raw:options.det_per_ring*options.sampling_raw,x,0:options.det_per_ring*options.sampling_raw-1,'linear','extrap')';
    y = interp1(1:options.sampling_raw:options.det_per_ring*options.sampling_raw,y,0:options.det_per_ring*options.sampling_raw-1,'linear','extrap')';
    if interpolateSinogram
        for hh = 1 : options.partitions
            % Create a (sparse) matrix and interpolate new detectors
            [K, ~, V] = find(options.SinM{hh});
            L = find(tril(true(options.detectors,options.detectors), 0));
            L = L(K);
            
            [I,J] = ind2sub([options.detectors options.detectors], L);
            clear L
            coincidences = sparse(I,J,V,options.detectors,options.detectors);
            clear I J V
            joku = single(full(coincidences));
            SinM_uus = zeros(size(coincidences,1)*options.sampling_raw, size(coincidences,2)*options.sampling_raw,'single');
            SinM_uus(1:end-1,1:end-1) = interp2(joku, options.sampling_raw / 2, options.sampling_interpolation_method_raw);
            SinM_uus = single(SinM_uus(tril(true(size(SinM_uus)), 0)));
            options.SinM{hh} = SinM_uus;
        end
        if options.verbose
            disp(['Raw data sampling increased by ' num2str(options.sampling_raw) 'x'])
        end
    end
else
    if isempty(x)
        [~, ~, x, y] = detector_coordinates(options);
        [x, y] = sinogram_coordinates_2D(options, x, y);
    end
    % Sinogram data
    
    xx1 = reshape(x(:,1),options.Ndist,options.Nang);
    xx2 = reshape(x(:,2),options.Ndist,options.Nang);
    yy1 = reshape(y(:,1),options.Ndist,options.Nang);
    yy2 = reshape(y(:,2),options.Ndist,options.Nang);

    [X,Y] = meshgrid(1:options.Nang, 1:options.sampling:options.Ndist*options.sampling+1);
    [Xq,Yq] = meshgrid(1:options.Nang, 1:options.Ndist*options.sampling+1);

    if strcmp(options.sampling_interpolation_method,'spline') || strcmp(options.sampling_interpolation_method, 'makima')
        x1 = interp2(X(1:end-1,:), Y(1:end-1,:) ,double(xx1), Xq(1:end-1,:), Yq(1:end-1,:), options.sampling_interpolation_method);
        x2 = interp2(X(1:end-1,:), Y(1:end-1,:) ,double(xx2), Xq(1:end-1,:), Yq(1:end-1,:), options.sampling_interpolation_method);
        y1 = interp2(X(1:end-1,:), Y(1:end-1,:) ,double(yy1), Xq(1:end-1,:), Yq(1:end-1,:), options.sampling_interpolation_method);
        y2 = interp2(X(1:end-1,:), Y(1:end-1,:) ,double(yy2), Xq(1:end-1,:), Yq(1:end-1,:), options.sampling_interpolation_method);
    else
        x1 = interp2(X, Y,double([xx1;xx1(end-1,:)]), Xq, Yq, options.sampling_interpolation_method);
        x1 = x1(1:end-1,:);
        x2 = interp2(X, Y,double([xx2;xx2(end-1,:)]), Xq, Yq, options.sampling_interpolation_method);
        x2 = x2(1:end-1,:);
        y1 = interp2(X, Y,double([yy1;yy1(end-1,:)]), Xq, Yq, options.sampling_interpolation_method);
        y1 = y1(1:end-1,:);
        y2 = interp2(X, Y,double([yy2;yy2(end-1,:)]), Xq, Yq, options.sampling_interpolation_method);
        y2 = y2(1:end-1,:);
    end
    x = [x1(:),x2(:)];
    y = [y1(:),y2(:)];
    
    % Interpolate the sinogram
    if interpolateSinogram
        
        if iscell(options.SinM)
            if numel(options.SinM{1}) == size(options.SinM{1},1)
                for kk = 1 : options.partitions
                    options.SinM{kk} = reshape(options.SinM{kk}, options.Ndist, options.Nang, numel(options.SinM)/(options.Ndist * options.Nang));
                end
            end
        else
            if numel(options.SinM) == size(options.SinM,1)
                options.SinM = reshape(options.SinM, options.Ndist, options.Nang, numel(options.SinM)/(options.Ndist * options.Nang));
            end
        end
        options.SinM = interpolateSinog(options.SinM, options.sampling, options.Ndist, options.partitions, options.sampling_interpolation_method);
        options.Ndist = options.Ndist * options.sampling;
        options.nRowsD = options.nRowsD * options.sampling;
        if options.verbose
            disp(['Sinogram sampling increased by ' num2str(options.sampling) 'x'])
        end
    end
end
end