function [x, y, options] = increaseSampling(options, x, y, interpolateSinogram)
%INCREASESAMPLING Increase the sampling rate
%   This function increases the sampling rate of the input detectors and
%   sinogram. Works only for sinogram data. Can be used to prevent aliasing
%   artifacts when using too high image resolution.

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
            [K, ~, V] = find(options.coincidences{hh});
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
            options.coincidences{hh} = SinM_uus;
        end
        if options.verbose
            disp(['Raw data sampling increased by ' num2str(options.sampling) 'x'])
        end
    end
else
    % Sinogram data
    
    xx1 = reshape(x(:,1),options.Ndist,options.Nang);
    xx2 = reshape(x(:,2),options.Ndist,options.Nang);
    yy1 = reshape(y(:,1),options.Ndist,options.Nang);
    yy2 = reshape(y(:,2),options.Ndist,options.Nang);
    
    
    % Form the new detector coordinates
    if exist('OCTAVE_VERSION','builtin') == 0 && exist('repelem', 'builtin') == 0
        joku = repeat_elem(xx1,options.sampling);
    else
        joku = repelem(xx1,options.sampling,1);
    end
%     % Gaps between detectors
    vali = [zeros(1, options.Nang);diff(xx1)];
    joku(1:options.sampling:end,:) = xx1 - (vali / options.sampling);
    xx = zeros(length(joku(:)), 2);
    xx(:,1) = joku(:);
    if exist('OCTAVE_VERSION','builtin') == 0 && exist('repelem', 'builtin') == 0
        joku = repeat_elem(xx2,options.sampling);
    else
        joku = repelem(xx2,options.sampling,1);
    end
    vali = [zeros(1, options.Nang);diff(xx2)];
    joku(1:options.sampling:end,:) = xx2 - (vali / options.sampling);
    xx(:,2) = joku(:);
    if exist('OCTAVE_VERSION','builtin') == 0 && exist('repelem', 'builtin') == 0
        joku = repeat_elem(yy1,options.sampling);
    else
        joku = repelem(yy1,options.sampling,1);
    end
    vali = [zeros(1, options.Nang);diff(yy1)];
    joku(1:options.sampling:end,:) = yy1 - (vali / options.sampling);
    yy = zeros(length(joku(:)), 2);
    yy(:,1) = joku(:);
    if exist('OCTAVE_VERSION','builtin') == 0 && exist('repelem', 'builtin') == 0
        joku = repeat_elem(yy2,options.sampling);
    else
        joku = repelem(yy2,options.sampling,1);
    end
    vali = [zeros(1, options.Nang);diff(yy2)];
    joku(1:options.sampling:end,:) = yy2 - (vali / options.sampling);
    yy(:,2) = joku(:);
    x = xx;
    y = yy;
    
    % Interpolate the sinogram
    if interpolateSinogram
        
        options.SinM = interpolateSinog(options.SinM, options.sampling, options.Ndist, options.partitions, options.sampling_interpolation_method);
        if options.verbose
            disp(['Sinogram sampling increased by ' num2str(options.sampling) 'x'])
        end
    end
end
end
