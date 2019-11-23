function [x, y, options] = arcCorrection(options, xp, yp, interpolateSinogram)
%ARCCORRECTION Performs arc correction on the detector coordinates
%   This function outputs the arc corrected detector coordinates and
%   sinogram. Works only with sinogram data.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi
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

mashing = options.det_w_pseudo / options.Nang / 2;
orig_xp = xp;
orig_yp = yp;
[x_o, y_o] = sinogram_coordinates_2D(options, orig_xp, orig_yp);

new_xp = zeros(size(xp));
new_yp = zeros(size(yp));
% The angles of the blocks/buckets
if mod(options.blocks_per_ring, 4) == 0
    l_angles = linspace(0,90, options.blocks_per_ring / 4 + 1);
else
    l_angles = linspace(0,180, options.blocks_per_ring / 2 + 1);
    l_angles = l_angles(1 : ceil(options.blocks_per_ring / 4));
end
ll = 1;

% Shift to the zero angle (with x-axis)
if options.det_w_pseudo > options.det_per_ring
    xp = circshift(xp, -options.offangle - (options.cryst_per_block + 1) / 2);
    yp = circshift(yp, -options.offangle - (options.cryst_per_block + 1) / 2);
else
    xp = circshift(xp, -options.offangle - options.cryst_per_block / 2);
    yp = circshift(yp, -options.offangle - options.cryst_per_block / 2);
end

xp = xp - options.diameter / 2;
yp = yp - options.diameter / 2;

% Determine the points along the circle that reside (approximately) on the
% same line as the original LOR
for kk = 1 : length(xp) / 4
    if options.det_w_pseudo > options.det_per_ring
        vali = kk + options.det_w_pseudo / 2 - (options.cryst_per_block + 1) : kk + options.det_w_pseudo / 2 + (options.cryst_per_block + 1);
    else
        vali = kk + options.det_w_pseudo / 2 - options.cryst_per_block : kk + options.det_w_pseudo / 2 + options.cryst_per_block;
    end
    kulma = l_angles(ll);
    angle1 = atand((yp(kk)-yp(vali))./(xp(kk)-xp(vali)));
    angle1 = round(angle1*10^4)/10^4;
    if mod(options.blocks_per_ring, 4) == 0
        ind1 = find(abs(angle1) == kulma);
        if isempty(ind1)
            [~, ind1] = min(abs(abs(angle1) - kulma));
        end
    else
        [~, ind1] = min(abs(abs(angle1) - kulma));
    end
    if options.det_w_pseudo > options.det_per_ring
        x2 = xp(kk + options.det_w_pseudo / 2 - (options.cryst_per_block + 1) + ind1 - 1);
        y2 = yp(kk + options.det_w_pseudo / 2 - (options.cryst_per_block + 1) + ind1 - 1);
    else
        x2 = xp(kk + options.det_w_pseudo / 2 - options.cryst_per_block + ind1 - 1);
        y2 = yp(kk + options.det_w_pseudo / 2 - options.cryst_per_block + ind1 - 1);
    end
    p = [xp(kk); yp(kk)];
    q = [x2; y2];
    d = q - p;
    l = ((-dot(2*p,d) - sqrt(dot(2*p,d)^2 - 4*dot(d,d)*(dot(p,p) - (options.diameter/2)^2)))/ (2*dot(d,d)));
    lx = p + l * d;
    new_xp(kk) = lx(1);
    new_yp(kk) = lx(2);
end

new_xp(1 : length(xp) / 4) = new_xp(1 : length(xp) / 4) + options.diameter/2;
new_yp(1 : length(yp) / 4) = new_yp(1 : length(yp) / 4) + options.diameter/2;

new_yp(kk + 1: kk * 2) = flip(new_yp(1: kk));
diffi = diff([flip(new_yp(1 : kk)) ; options.diameter/2]);
diffi = options.diameter/2 + cumsum(flip(diffi));
new_yp(kk * 2 + 1 : end) = [diffi ; flip(diffi)];

diffi = diff([new_xp(1 : kk) ; options.diameter/2]);
diffi = options.diameter/2 + cumsum(flip(diffi));
new_xp(kk + 1: kk * 2) = diffi;
new_xp(kk * 2 + 1 : end) = flip(new_xp(1 : kk * 2));
xp = new_xp;
yp = new_yp;

if options.det_w_pseudo > options.det_per_ring
    xp = circshift(xp, options.offangle + (options.cryst_per_block + 1) / 2);
    yp = circshift(yp, options.offangle + (options.cryst_per_block + 1) / 2);
else
    xp = circshift(xp, options.offangle + options.cryst_per_block / 2);
    yp = circshift(yp, options.offangle + options.cryst_per_block / 2);
end

[x, y] = sinogram_coordinates_2D(options, xp, yp);


xx1 = reshape(x(:,1),options.Ndist,options.Nang);
xx2 = reshape(x(:,2),options.Ndist,options.Nang);
yy1 = reshape(y(:,1),options.Ndist,options.Nang);
yy2 = reshape(y(:,2),options.Ndist,options.Nang);

% Exchange coordinates
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

x(:,1) = xx1(:);
x(:,2) = xx2(:);
y(:,1) = yy1(:);
y(:,2) = yy2(:);

% y_orig = y;
% x_orig = x;

angle = atand((yy1-yy2)./(xx1-xx2)) + 90;
angle(angle == 180) = 0;
[~, J] = min(mean(angle));

% Create the arc corrected coordinates from the perpendicular LOR
% coordinates
eka = xx1(1,J);
vika = xx1(end,J);
alkux = linspace(eka, vika, options.Ndist);

% Determine the y-coordinates
alkuy = sqrt((options.diameter/2)^2 - (alkux - options.diameter/2).^2) + options.diameter/2;

alku = [alkux; alkuy];

alku2 = alku - options.diameter/2;

alku = [alkux; abs(alkuy - options.diameter)];

alku1 = alku - options.diameter/2;


angles = linspace(0, 180, options.Nang + 1);
angles = angles(1 : end - 1);
% [~, JJ] = min(abs(angles));
% angles = circshift(angles, JJ-1);
angles = circshift(angles, J);
angles = reshape(angles, 1, 1, []);

% Rotation matrix
rot_matrix = [cosd(angles) -sind(angles); sind(angles) cosd(angles)];

rot_matrix = squeeze(num2cell(rot_matrix, [1 2]))';

% Arc correction
new_xy1 = cell2mat(cellfun(@(x) x * alku1, rot_matrix, 'UniformOutput', false))' + options.diameter/2;
new_xy2 = cell2mat(cellfun(@(x) x * alku2, rot_matrix, 'UniformOutput', false))' + options.diameter/2;


x = [new_xy1(:,1), new_xy2(:,1)];
y = [new_xy1(:,2), new_xy2(:,2)];

xx1 = reshape(x(:,1),options.Ndist,options.Nang);
xx2 = reshape(x(:,2),options.Ndist,options.Nang);
apu = xx1(:, 1: J);
xx1(:, 1: J) = flipud(xx2(:, 1: J));
xx2(:, 1: J) = flipud(apu);
yy1 = reshape(y(:,1),options.Ndist,options.Nang);
yy2 = reshape(y(:,2),options.Ndist,options.Nang);
apu = yy1(:, 1: J);
yy1(:, 1: J) = flipud(yy2(:, 1: J));
yy2(:, 1: J) = flipud(apu);

x(:,1) = xx1(:);
x(:,2) = xx2(:);
y(:,1) = yy1(:);
y(:,2) = yy2(:);

xx1_o = reshape(x_o(:,1),options.Ndist,options.Nang);
xx2_o = reshape(x_o(:,2),options.Ndist,options.Nang);
yy1_o = reshape(y_o(:,1),options.Ndist,options.Nang);
yy2_o = reshape(y_o(:,2),options.Ndist,options.Nang);

if interpolateSinogram
    
    
    angle_o = atand((yy1_o-yy2_o)./(xx1_o-xx2_o)) + 90;
    angle = atand((yy1-yy2)./(xx1-xx2)) + 90;
    
    
    x0 = options.diameter/2;
    y0 = options.diameter/2;
    distance_o = ((abs((y_o(:,1)-y_o(:,2))*x0 - (x_o(:,1) - x_o(:,2))*y0 + x_o(:,1).*y_o(:,2) - y_o(:,1).*x_o(:,2))./sqrt((y_o(:,1)-y_o(:,2)).^2 + (x_o(:,1)-x_o(:,2)).^2)));
    distance_o = reshape(distance_o, options.Ndist, options.Nang);
    distance_o(options.Ndist/2 + 1 : end,:) = -distance_o(options.Ndist/2 + 1 : end,:);
    distance = ((abs((y(:,1)-y(:,2))*x0 - (x(:,1) - x(:,2))*y0 + x(:,1).*y(:,2) - y(:,1).*x(:,2))./sqrt((y(:,1)-y(:,2)).^2 + (x(:,1)-x(:,2)).^2)));
    distance = reshape(distance, options.Ndist, options.Nang);
    distance(options.Ndist/2 + 1 : end,:) = -distance(options.Ndist/2 + 1 : end,:);
    
    % Interpolate the sinogram
    tic
    uus_SinM = zeros(size(options.SinM));
    if license('test','Distrib_Computing_Toolbox')
        temp = options.SinM;
        Ndist = options.Ndist;
        Nang = options.Nang;
        arc_interpolation = options.arc_interpolation;
        try
            parfor kk = 1 : size(options.SinM,3)
                if exist('scatteredInterpolant', 'file') == 2 && strcmp('cubic',arc_interpolation) && strcmp('v4',arc_interpolation)
                    F = scatteredInterpolant(angle_o(:), distance_o(:), reshape(double(temp(:,:,kk)), Ndist*Nang,1));
                    uus_SinM(:,:,kk) = F(angle, distance);
                else
                    uus_SinM(:,:,kk) = griddata(angle_o, distance_o, double(temp(:,:,kk)), angle, distance, arc_interpolation);
                end
            end
        catch
            for kk = 1 : size(options.SinM,3)
                if exist('scatteredInterpolant', 'file') == 2 && strcmp('cubic',options.arc_interpolation) && strcmp('v4',options.arc_interpolation)
                    F = scatteredInterpolant(angle_o(:), distance_o(:), reshape(double(options.SinM(:,:,kk)), options.Ndist*options.Nang,1));
                    uus_SinM(:,:,kk) = F(angle, distance);
                else
                    uus_SinM(:,:,kk) = griddata(angle_o, distance_o, double(options.SinM(:,:,kk)), angle, distance, options.arc_interpolation);
                end
            end
        end
    else
        for kk = 1 : size(options.SinM,3)
            if exist('scatteredInterpolant', 'file') == 2 && strcmp('cubic',options.arc_interpolation) && strcmp('v4',options.arc_interpolation)
                F = scatteredInterpolant(angle_o(:), distance_o(:), reshape(double(options.SinM(:,:,kk)), options.Ndist*options.Nang,1));
                uus_SinM(:,:,kk) = F(angle, distance);
            else
                uus_SinM(:,:,kk) = griddata(angle_o, distance_o, double(options.SinM(:,:,kk)), angle, distance, options.arc_interpolation);
            end
        end
    end
    endTime = toc;
    options.SinM = uus_SinM;
    if options.verbose
        disp(['Arc correction complete in ' num2str(endTime) ' seconds'])
    end
end
end