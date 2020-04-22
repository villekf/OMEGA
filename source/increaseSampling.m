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
    mashing = options.det_w_pseudo / options.Nang / 2;
    
    % If no arc correction is present
    if ~options.arc_correction
        % Swap corners, make the data continuous
        testi = abs(diff(xx1));
        [I,J] = find(testi > options.cr_p*2);
        testi2 = abs(diff(J));
        ind1 = find(testi2 > mean(testi2)*2, 1, 'first');
        if xx1(1,J(ind1)) <= xx1(1,J(ind1) + 1)
            indices2 = J(ind1) + 1 : - 1 : 1;
            indices2 = [indices2(1);repelem(indices2(2:end),mashing * 2)'];
        elseif xx1(1,J(ind1)) <= xx1(2,J(ind1))
            indices2 = J(ind1) : - 1 : 1;
            indices2 = [indices2(1);indices2(1);repelem(indices2(2:end),mashing * 2)'];
        else
            indices2 = J(ind1) : - 1 : 1;
            indices2 = [indices2(1);repelem(indices2(2:end),mashing * 2)'];
        end
        indices1 = false(size(xx1));
        for kk = 1 : I(1)
            indices1(kk,1:indices2(kk)) = true(indices2(kk),1);
        end
        testi = abs(diff(fliplr(xx1)));
        [I,J] = find(testi > options.cr_p*2);
        testi2 = abs(diff(J));
        ind1 = find(testi2 > mean(testi2)*2, 1, 'first');
        if xx1(1,options.Nang - J(ind1) + 1) >= xx1(2,options.Nang - J(ind1) + 1)
            indices2 = options.Nang - J(ind1) + 1 : options.Nang;
            indices2 = [indices2(1);indices2(1);repelem(indices2(2:end),mashing * 2)'];
        elseif xx1(1,options.Nang - J(ind1) + 1) >= xx1(1,options.Nang - J(ind1))
            indices2 = options.Nang - J(ind1) + 1 : options.Nang;
            indices2 = [indices2(1);repelem(indices2(2:end),mashing * 2)'];
        else
            indices2 = options.Nang - J(ind1) + 1 : options.Nang;
            indices2 = [indices2(1);repelem(indices2(2:end),mashing * 2)'];
        end
        for kk = 1 : I(1)
            indices1(kk,indices2(kk):end) = true(length(indices1(kk,indices2(kk):end)),1);
        end
        temp = xx1(indices1);
        xx1(indices1) = xx2(indices1);
        xx2(indices1) = temp;
        temp = yy1(indices1);
        yy1(indices1) = yy2(indices1);
        yy2(indices1) = temp;
        
        
        minimix = min(x(:));
        minimiy = min(y(:));
        maksimiy = max(y(:));
        maksimix = max(x(:));
        
        % Mark the new coordinates as NaN and interpolate them
        testi = [diff(xx1);ones(1, options.Nang)];
        xx1(testi == 0 & xx1 > minimix & xx1 < maksimix) = NaN;
        xx1 = fillmissing(xx1,'linear');
        xx1(1,:) = circshift(xx1(2,:),1);
        xx1(1,1) = max(xx1(:)) - xx1(1,1);
        xx1(1,:) = (xx1(1,:) + xx1(2,:)) / 2;
        testi = [zeros(1, options.Nang) ; diff(xx2)];
        xx2(testi == 0 & xx2 > minimix & xx2 < maksimix) = NaN;
        xx2 = fillmissing(xx2,'linear');
        testi = [diff(yy1); ones(1, options.Nang)];
        yy1(testi == 0 & yy1 > minimiy & yy1 < maksimiy) = NaN;
        yy1 = fillmissing(yy1,'linear');
        yy1(1,:) = circshift(yy1(2,:),1);
        yy1(1,1) = maksimiy - yy1(1,1);
        yy1(1,:) = (yy1(1,:) + yy1(2,:)) / 2;
        testi = [zeros(1, options.Nang) ; diff(yy2)];
        yy2(testi == 0 & yy2 > minimiy & yy2 < maksimiy) = NaN;
        yy2 = fillmissing(yy2,'linear');
        
        x(:,1) = xx1(:);
        x(:,2) = xx2(:);
        y(:,1) = yy1(:);
        y(:,2) = yy2(:);
    end
    
    % Form the new detector coordinates
    xx1 = reshape(x(:,1),options.Ndist,options.Nang);
    if exist('OCTAVE_VERSION','builtin') == 0 && exist('repelem', 'builtin') == 0
        joku = repeat_elem(xx1,options.sampling);
    else
        joku = repelem(xx1,options.sampling,1);
    end
    % Gaps between detectors
    vali = (xx1(1:2:end,:) - xx1(2:2:end,:));
    for kk = 1 : options.sampling - 1
        joku(kk + 1:options.sampling *2:end,:) = xx1(1:2:end,:) - (vali / options.sampling) * kk;
    end
    vali = (xx1(2:2:end,:) - [xx1(3:2:end,:);xx1(end,:) + abs(xx1(end-1,:) - xx1(end,:))]);
    for kk = 1 : options.sampling - 1
        joku(options.sampling + 2 + (kk - 1):options.sampling *2:end,:) = xx1(2:2:end,:) - (vali / options.sampling) * kk;
    end
    xx = zeros(length(joku(:)), 2);
    xx(:,1) = joku(:);
    xx1 = reshape(x(:,2),options.Ndist,options.Nang);
    if exist('OCTAVE_VERSION','builtin') == 0 && exist('repelem', 'builtin') == 0
        joku = repeat_elem(xx1,options.sampling);
    else
        joku = repelem(xx1,options.sampling,1);
    end
    vali = (xx1(1:2:end,:) - xx1(2:2:end,:));
    for kk = 1 : options.sampling - 1
        joku(kk + 1:options.sampling *2:end,:) = xx1(1:2:end,:) - (vali / options.sampling) * kk;
    end
    vali = (xx1(2:2:end,:) - [xx1(3:2:end,:);xx1(end,:) + abs(xx1(end-1,:) - xx1(end,:))]);
    for kk = 1 : options.sampling - 1
        joku(options.sampling + 2 + (kk - 1):options.sampling *2:end,:) = xx1(2:2:end,:) - (vali / options.sampling) * kk;
    end
    xx(:,2) = joku(:);
    xx1 = reshape(y(:,1),options.Ndist,options.Nang);
    if exist('OCTAVE_VERSION','builtin') == 0 && exist('repelem', 'builtin') == 0
        joku = repeat_elem(xx1,options.sampling);
    else
        joku = repelem(xx1,options.sampling,1);
    end
    vali = (xx1(1:2:end,:) - xx1(2:2:end,:));
    for kk = 1 : options.sampling - 1
        joku(kk + 1:options.sampling *2:end,:) = xx1(1:2:end,:) - (vali / options.sampling) * kk;
    end
    vali = (xx1(2:2:end,:) - [xx1(3:2:end,:);xx1(end,:) + abs(xx1(end-1,:) - xx1(end,:))]);
    for kk = 1 : options.sampling - 1
        joku(options.sampling + 2 + (kk - 1):options.sampling *2:end,:) = xx1(2:2:end,:) - (vali / options.sampling) * kk;
    end
    yy = zeros(length(joku(:)), 2);
    yy(:,1) = joku(:);
    xx1 = reshape(y(:,2),options.Ndist,options.Nang);
    if exist('OCTAVE_VERSION','builtin') == 0 && exist('repelem', 'builtin') == 0
        joku = repeat_elem(xx1,options.sampling);
    else
        joku = repelem(xx1,options.sampling,1);
    end
    vali = (xx1(1:2:end,:) - xx1(2:2:end,:));
    for kk = 1 : options.sampling - 1
        joku(kk + 1:options.sampling *2:end,:) = xx1(1:2:end,:) - (vali / options.sampling) * kk;
    end
    vali = (xx1(2:2:end,:) - [xx1(3:2:end,:);xx1(end,:) + abs(xx1(end-1,:) - xx1(end,:))]);
    for kk = 1 : options.sampling - 1
        joku(options.sampling + 2 + (kk - 1):options.sampling *2:end,:) = xx1(2:2:end,:) - (vali / options.sampling) * kk;
    end
    yy(:,2) = joku(:);
    x = xx;
    y = yy;
    
    % Interpolat the sinogram
    if interpolateSinogram
        
        if iscell(options.SinM)
            SinM_uus = cell(size(options.SinM));
            for hh = 1 : options.partitions
                SinM_uus{hh} = zeros(size(options.SinM{hh},1)*options.sampling, size(options.SinM{hh},2),size(options.SinM{hh},3),'single');
                for kk = 1 : size(options.SinM{hh},3)
                    for ll = 1 : size(options.SinM{hh},2)
                        SinM_uus{hh}(:,ll,kk) = single(interp1(1:options.sampling:options.Ndist*options.sampling+1, double([options.SinM{hh}(:,ll,kk);options.SinM{hh}(end - 1,ll,kk)]), ...
                            1:options.Ndist*options.sampling, options.sampling_interpolation_method));
                    end
                end
            end
            options.SinM = SinM_uus;
        else
            SinM_uus = zeros(size(options.SinM,1)*options.sampling, size(options.SinM,2),size(options.SinM,3));
            for kk = 1 : size(options.SinM,3)
                for ll = 1 : size(options.SinM,2)
                    SinM_uus(:,ll,kk) = interp1(1:options.sampling:options.Ndist*options.sampling+1, double([options.SinM(:,ll,kk);options.SinM(end - 1,ll,kk)]), ...
                        1:options.Ndist*options.sampling, options.sampling_interpolation_method);
                end
            end
            options.SinM = single(SinM_uus);
        end
        if options.verbose
            disp(['Sinogram sampling increased by ' num2str(options.sampling) 'x'])
        end
    end
end
end
