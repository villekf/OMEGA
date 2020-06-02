%% Visualization for the PET reconstructions
% Each section has a visualization code for a different purpose
% Only a specific section should be run at a time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Ville-Veikko Wettenhovi
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

% Load the saved reconstruction and machine specific variables
if exist('pz','var')
    image_properties = pz{end,1};
end



%% Visualize several reconstructions for one time step for all slices, last iterations

algo_char = algorithms_char();

% 1 = MLEM, 2 = OSEM, 3 = MRAMLA, 4 = RAMLA, 5 = ROSEM, 6 = RBI, 7 = DRAMA, 
% 8 = COSEM, 9 = ECOSEM, 10 = ACOSEM, 11 = MRP-OSL-OSEM, 12 = MRP-OSL-MLEM, 
% 13 = MRP-BSREM, 14 = MRP-MBSREM, 15 = MRP-ROSEM, 16 = MRP-RBI, 17 = MRP-OSL-COSEM, 
% 18 = QP (OSL-OSEM), 19 = QP (OSL-MLEM), 20 = QP (BSREM), 21 = QP (MBSREM), 
% 22 = QP (ROSEM), 23 = QP (RBI), 24 =  QP (OSL-COSEM), % 25 = HP (OSL-OSEM), 
% 26 = HP (OSL-MLEM), 27 = HP (BSREM), 28 = HP (MBSREM), 29 = HP (ROSEM), 
% 30 = HP (RBI), 31 =  HP (OSL-COSEM), 32 = L-filter (OSL-OSEM), 
% 33 = L-filter (OSL-MLEM), 34 =  L-filter (BSREM), 35 =  L-filter (MBSREM), 
% 36 = L-filter (ROSEM), 37 = L-filter (RBI), 38 = L-filter (OSL-COSEM), 
% 39 = FMH (OSL-OSEM), 40 = FMH (OSL-MLEM), 41 = FMH (BSREM), 42 = FMH (MBSREM), 
% 43 = FMH (ROSEM), 44 = FMH (RBI), 45 = FMH (OSL-COSEM), 46 = Weighted mean (OSL-OSEM), 
% 47 = Weighted mean (OSL-MLEM), 48 = Weighted mean (BSREM), 49 = Weighted mean (MBSREM), 
% 50 = Weighted mean (ROSEM), 51 = Weighted mean (RBI), 52 = Weighted mean (OSL-COSEM), 
% 53 = Total variation (OSL-OSEM), 54 = Total variation (OSL-MLEM), 55 = Total variation (BSREM), 
% 56 = Total variation (MBSREM), 57 = Total variation (ROSEM), 58 = Total variation (RBI), 
% 59 = Total variation (OSL-COSEM), 60 = Anisotropic Diffusion (OSL-OSEM), 
% 61 = Anisotropic Diffusion (OSL-MLEM), 62 = Anisotropic Diffusion (BSREM), 
% 63 = Anisotropic Diffusion (MBSREM), 64 = Anisotropic Diffusion (ROSEM), 
% 65 = Anisotropic Diffusion (RBI), 66 = Anisotropic Diffusion (OSL-COSEM), 
% 67 = APLS (OSL-OSEM), 68 = APLS (OSL-MLEM), 69 = APLS (BSREM), 70 = APLS (MBSREM), 
% 71 = APLS (ROSEM), 72 = APLS (RBI), 73 = APLS (OSL-COSEM), 74 = TGV (OSL-OSEM), 
% 75 = TGV (OSL-MLEM), 76 = TGV (BSREM), 77 = TGV (MBSREM), 78 = TGV (ROSEM), 
% 79 = TGV (RBI), 80 = TGV (OSL-COSEM), 81 = NLM (OSL-OSEM), 82 =  NLM (OSL-MLEM), 
% 83 = NLM (BSREM), 84 = NLM (MBSREM), 85 = NLM (ROSEM), 86 = NLM (RBI), 87 = NLM (OSL-COSEM), 
% 88 = Custom prior (OSL-OSEM), 89 =  Custom prior (OSL-MLEM), 90 = Custom prior (BSREM), 
% 92 = Custom prior (MBSREM), 92 = Custom prior (ROSEM), 93 = Custom prior (RBI), 
% 94 = Custom prior (OSL-COSEM)
% Inputing algorithm number that does not exist in the cell array shows all
% the available algorithms present in the cell array
algorithms = [2];
% Use this value to scale the color scale in the image (higher values make
% low count areas brighter)
color_scale = 1;
% From which reconstruction should the color scale be taken
% If zero, then each algorithm has its own color scala (from zero to their
% own maximum value, i.e. there is no global limit)
% NOTE: The numbering is according to the length of the above algorithms
% vector, e.g. if you have algorithms = [2, 4, 5] and color_from_algo = 2
% then the scale will be taken from RAMLA reconstruction (second element of
% algorithms) 
color_from_algo = 0;
% Visualization plane
% Choose the plane where the visualization takes place
% 1 = Transverse, 2 = Coronal/frontal, 3 = Sagittal
v_plane = 1;

if exist('f_osem','var') && ~exist('pz','var')
    pz = cell(95,1);
    pz{2} = f_osem;
end
img = check_algorithms(pz, algorithms, color_from_algo);
if isempty(img)
    return
end


if length(algorithms) > 3
    hh = 2;
else
    hh = 1;
end
if length(algorithms) < 4
    jj = min(3, length(algorithms));
else
    jj = 3;
end
set(0,'units','pixels')
gg = get(0,'ScreenSize');
if jj > 4
    im_size = gg(4)/(2.5 + (jj - 4)/2);
else
    im_size = gg(4)/2.5;
end
figure
set(gcf, 'Position', [min(gg(3)/2-im_size/2*jj,gg(3)), min(gg(4)/2-im_size/2*hh, im_size * jj), gg(4), im_size * hh]);

clim = [0 max(max(max(max(img(:,:,:,end)))))/color_scale];
if v_plane == 2
    koko = size(img,1);
elseif v_plane == 3
    koko = size(img,2);
else
    koko = size(img,3);
end
for kk = 1 : koko
    for ll = 1 : length(algorithms)
        
        img = pz{algorithms(ll)};
        if v_plane == 2
            img = rot90(permute(img, [3 2 1 4]),2);
        elseif v_plane == 3
            img = permute(img, [1 3 2 4]);
        end
        subplot(hh, jj, ll)
        if color_from_algo == 0
            clim = [0 max(max(max(max(img(:,:,:,end)))))/color_scale];
            imagesc(img(:,:,kk,end),clim)
        else
            imagesc(img(:,:,kk,end),clim)
        end
        axis image off
        title([char(algo_char(algorithms(ll))) ' slice = ' num2str(kk)])
    end
    pause%(0.1)
    drawnow
end


%% Visualize N iterations of a single reconstruction for one time step for all slices

algo_char = algorithms_char();

algorithms = 2;
% Use this value to scale the color scale in the image (higher values make
% low count areas brighter)
color_scale = 1;
% Visualization plane
% Choose the plane where the visualization takes place
% 1 = Transverse, 2 = Coronal/frontal, 3 = Sagittal
v_plane = 1;

% How many iterations in the image 
% Initial values can be included as iterations as they are saved as such
% N_iter LAST iterations will be used for visualization
N_iter = 4;


color_from_algo = 1;
if exist('f_osem','var') && ~exist('pz','var')
    pz = cell(95,1);
    pz{2} = f_osem;
end
img = check_algorithms(pz, algorithms, color_from_algo);
if isempty(img)
    return
end


if N_iter > 3
    hh = 2;
else
    hh = 1;
end
if N_iter < 4
    jj = min(3, N_iter);
else
    jj = 3;
end
set(0,'units','pixels')
gg = get(0,'ScreenSize');
if jj > 4
    im_size = gg(4)/(2.5 + (jj - 4)/2);
else
    im_size = gg(4)/2.5;
end
figure
set(gcf, 'Position', [min(gg(3)/2-im_size/2*jj,gg(3)), min(gg(4)/2-im_size/2*hh, im_size * jj), gg(4), im_size * hh]);
if v_plane == 2
    koko = size(img,1);
    img = rot90(permute(img, [3 2 1 4]),2);
elseif v_plane == 3
    koko = size(img,2);
    img = permute(img, [1 3 2 4]);
else
    koko = size(img,3);
end
for kk = 1 : koko
    for ll = 1 : N_iter
        clim = [0 max(max(max(max(img(:,:,:,end - ll + 1)))))/color_scale];
        subplot(hh, jj, ll)
        imagesc(img(:,:,kk,end - ll + 1),clim)
        axis image
        title([char(algo_char(algorithms)) ', iteration = ' num2str(size(img,4) - ll) ', slice = ' num2str(kk)])
    end
    pause(0.25)
    drawnow
end

%% Compare several reconstructions for one time step for all slices, last iteration, with the source image obtained from GATE data
% NOTE: This is valid only for GATE data

algo_char = algorithms_char();

% This can be used to compare the achieved reconstruction with the
% coordinates from which the actual radioactive decay happened
% I.e. these allow for the error checking of the reconstructed data

algorithms = [2];
% Use this value to scale the color scale in the image (higher values make
% low count areas brighter)
color_scale = 1;
% From according to which reconstruction should the color scale be taken
% If zero, then each algorithm has its own color scala (from zero to their
% own maximum value, i.e. there is no global limit)
% NOTE: The numbering is according to the length of algorithms vector, e.g.
% if you have algorithms = [2, 4, 5] and color_from_algo = 2 then the scale
% will be taken from RAMLA reconstruction (second element of algorithms)
color_from_algo = 1;
% How is the source image formed?
% 1 = Form the source image by using only coincidences that originate from
% the very same location (source coordinates are the same) 
% 2 = Form the source image by using only the first single
% 3 = Form the source image by using only the second single
% 4 = Form the source image by using both singles (singles mode)
% 5 = Form the source image by using both singles and then dividing the
% counts by two
% 6 = Form the source image by using the average coordinates from both
% singles
% 7 = Form the source image by using the true coincidences (requires
% obtain_trues = true)
source_coordinates = 1;
% Visualization plane
% Choose the plane where the visualization takes place
% 1 = Transverse, 2 = Coronal/frontal, 3 = Sagittal
v_plane = 1;

% The source data was obtained from
% 1 = ASCII, 2 = LMF, 3 = ROOT
source = 1;

if exist('f_osem','var') && ~exist('pz','var')
    pz = cell(95,1);
    pz{2} = f_osem;
end
img = check_algorithms(pz, algorithms, color_from_algo);
if isempty(img)
    return
end

if length(algorithms) + 1 > 3
    hh = 2;
else
    hh = 1;
end
if length(algorithms) + 1 < 4
    jj = min(3, length(algorithms) + 1);
else
    jj = 3;
end
set(0,'units','pixels')
gg = get(0,'ScreenSize');
if jj > 4
    im_size = gg(4)/(2.5 + (jj - 4)/2);
else
    im_size = gg(4)/2.5;
end
figure
set(gcf, 'Position', [min(gg(3)/2-im_size/2*jj,gg(3)), min(gg(4)/2-im_size/2*hh, im_size * jj), gg(4), im_size * hh]);

if source == 1
    load([image_properties.machine_name '_Ideal_image_coordinates_' image_properties.name '_ASCII.mat'])
elseif source == 2
    load([image_properties.machine_name '_Ideal_image_coordinates_' image_properties.name '_LMF.mat'])
elseif source == 3
    load([image_properties.machine_name '_Ideal_image_coordinates_' image_properties.name '_ROOT.mat'])
end

if source_coordinates == 1
    FOV = C{1};
elseif source_coordinates == 2
    FOV = C{2};
elseif source_coordinates == 3
    FOV = C{3};
elseif source_coordinates == 4
    FOV = C{4};
elseif source_coordinates == 5
    FOV = C{5};
elseif source_coordinates == 6
    FOV = C{6};
elseif source_coordinates == 7
    FOV = C{7};
end


clim = [0 max(max(max(max(img(:,:,:,end)))))/color_scale];
clim2 = [0 max(max(max(max(FOV(:,:,:)))))/color_scale];
if v_plane == 2
    koko = size(img,1);
    FOV = rot90(permute(FOV, [3 2 1]),2);
elseif v_plane == 3
    koko = size(img,2);
    FOV = permute(FOV, [1 3 2]);
else
    koko = size(img,3);
end
for kk = 1 : koko
    for ll = 1 : length(algorithms)
        img = pz{algorithms(ll)};
        if v_plane == 2
            img = rot90(permute(img, [3 2 1 4]),2);
        elseif v_plane == 3
            img = permute(img, [1 3 2 4]);
        end
        subplot(hh, jj, ll)
        if color_from_algo == 0
            clim = [0 max(max(max(max(img(:,:,:,end)))))/color_scale];
            imagesc(img(:,:,kk,end),clim)
        else
            imagesc(img(:,:,kk,end),clim)
        end
        axis image
        title([char(algo_char(algorithms(ll))) ' slice = ' num2str(kk)])
    end
    subplot(hh, jj, ll + 1)
    imagesc(FOV(:,:,kk),clim2)
    axis image
    title(['Original decay image, slice = ' num2str(kk)])
    pause(0.25)
    drawnow
end

%% Examine the entire volume for one reconstruction
% NOTE: Use of this section requires vol3D v2
% Download: 
% https://se.mathworks.com/matlabcentral/fileexchange/22940-vol3d-v2

algo_char = algorithms_char();

algorithms = 2;
% The scale value for the pixel alpha values. Higher values will make the
% pixels more transparent, allowing areas of higher activity to be seen
% through background noise
alpha_scale = 10;

if exist('f_osem','var') && ~exist('pz','var')
    pz = cell(95,1);
    pz{2} = f_osem;
end
img = check_algorithms(pz, algorithms, 1);
if isempty(img)
    return
end

alpha_scaling = max(max(max(img(:,:,:,end)))) * alpha_scale;
alpha = permute(img(:,:,:,end), [3 2 1 4]);
alpha = alpha./alpha_scaling;
alpha(alpha > 1) = 1;

figure;vol3d('CData', permute(img(:,:,:,end), [3 2 1 4]), 'Alpha', alpha);
set(gca, 'View', [45 30]);
set(gca, 'XLim', [60 200]);
set(gca, 'ZLim', [80 200]);
set(gcf, 'Color', 'w');

%% Visualize simultanously all the views for n algorithms
% NOTE: Due to likely different number of slices, the smallers views will
% not be updated once they reach maximum number of slices

algo_char = algorithms_char();

algorithms = [2];
% Use this value to scale the color scale in the image (higher values make
% low count areas brighter)
color_scale = 1;
% From which reconstruction should the color scale be taken
% If zero, then each algorithm has its own color scala (from zero to their
% own maximum value, i.e. there is no global limit)
% NOTE: The numbering is according to the length of the above algorithms
% vector, e.g. if you have algorithms = [2, 4, 5] and color_from_algo = 2
% then the scale will be taken from RAMLA reconstruction (second element of
% algorithms)
color_from_algo = 1;

if exist('f_osem','var') && ~exist('pz','var')
    pz = cell(95,1);
    pz{2} = f_osem;
end
img = check_algorithms(pz, algorithms, color_from_algo);
if isempty(img)
    return
end

hh = numel(algorithms);
jj = 3;
set(0,'units','pixels')
gg = get(0,'ScreenSize');
im_size = gg(4)/2.5;
figure
set(gcf, 'Position', [min(gg(3)/2-im_size/2*jj,gg(3)), min(gg(4)/2-im_size/2*hh, gg(4)), im_size * jj, im_size * hh]);
img = pz{algorithms(1)};
koko1 = size(img,1);
koko2 = size(img,2);
koko3 = size(img,3);

clim = [0 max(max(max(max(img(:,:,:,end)))))/color_scale];
for kk = 1 : max([koko1, koko2, koko3])
    for ll = 1 : numel(algorithms)
        img = pz{algorithms(ll)};
        koko1 = size(img,1);
        koko2 = size(img,2);
        koko3 = size(img,3);
        if kk <= koko3
            subplot(hh, jj, 1 + jj*(ll - 1))
            if color_from_algo == 0
                clim = [0 max(max(max(max(img(:,:,:,end)))))/color_scale];
                imagesc(img(:,:,kk,end),clim)
            else
                imagesc(img(:,:,kk,end),clim)
            end
            title([char(algo_char(algorithms(ll))) ' slice = ' num2str(kk) ', transverse'])
        else
            subplot(hh, jj, 1 + jj*(ll - 1))
            if color_from_algo == 0
                clim = [0 max(max(max(max(img(:,:,:,end)))))/color_scale];
                imagesc(img(:,:,koko3,end),clim)
            else
                imagesc(img(:,:,koko3,end),clim)
            end
            title([char(algo_char(algorithms(ll))) ' slice = ' num2str(koko1) ', transverse'])
        end
        axis image
        if kk <= koko1
            img = rot90(permute(img, [3 2 1 4]),2);
            subplot(hh, jj, 2 + jj*(ll - 1))
            if color_from_algo == 0
                clim = [0 max(max(max(max(img(:,:,:,end)))))/color_scale];
                imagesc(img(:,:,kk,end),clim)
            else
                imagesc(img(:,:,kk,end),clim)
            end
            title([char(algo_char(algorithms(ll))) ' slice = ' num2str(kk) ', frontal'])
        else
            img = rot90(permute(img, [3 2 1 4]),2);
            subplot(hh, jj, 2 + jj*(ll - 1))
            if color_from_algo == 0
                clim = [0 max(max(max(max(img(:,:,:,end)))))/color_scale];
                imagesc(img(:,:,koko1,end),clim)
            else
                imagesc(img(:,:,koko1,end),clim)
            end
            title([char(algo_char(algorithms(ll))) ' slice = ' num2str(koko1) ', frontal'])
        end
        axis image
        if kk <= koko2
            img = permute(img, [3 1 2 4]);
            subplot(hh, jj, 3 + jj*(ll - 1))
            if color_from_algo == 0
                clim = [0 max(max(max(max(img(:,:,:,end)))))/color_scale];
                imagesc(img(:,:,kk,end),clim)
            else
                imagesc(img(:,:,kk,end),clim)
            end
            title([char(algo_char(algorithms(ll))) ' slice = ' num2str(kk) ', sagittal'])
        else
            img = permute(img, [3 1 2 4]);
            subplot(hh, jj, 3 + jj*(ll - 1))
            if color_from_algo == 0
                clim = [0 max(max(max(max(img(:,:,:,end)))))/color_scale];
                imagesc(img(:,:,koko2,end),clim)
            else
                imagesc(img(:,:,koko2,end),clim)
            end
            title([char(algo_char(algorithms(ll))) ' slice = ' num2str(koko2) ', sagittal'])
        end
        axis image
    end
    pause(0.05)
    drawnow
end

%% Dynamic visualization
% Time series of images, n reconstructions, optionally also the
% "true" image

algo_char = algorithms_char();

algorithms = [2];
% Use this value to scale the color scale in the image (higher values make
% low count areas brighter)
color_scale = 1;
% From according to which reconstruction should the color scale be taken
% If zero, then each algorithm has its own color scala (from zero to their
% own maximum value, i.e. there is no global limit)
% NOTE: The numbering is according to the length of algorithms vector, e.g.
% if you have algorithms = [2, 4, 5] and color_from_algo = 2 then the scale
% will be taken from RAMLA reconstruction (second element of algorithms)
color_from_algo = 1;
% Visualization plane
% Choose the plane where the visualization takes place
% 1 = Transverse, 2 = Coronal/frontal, 3 = Sagittal
v_plane = 1;
% From which slice is the dynamic time series obtained?
slice = 40;
% The source data was obtained from
% 0 = No source image, 1 = ASCII, 2 = LMF, 3 = ROOT
source = 0;
% How is the source image formed?
% 1 = Form the source image by using only coincidences that originate from
% the very same location (source coordinates are the same) 
% 2 = Form the source image by using only the first single
% 3 = Form the source image by using only the second single
% 4 = Form the source image by using both singles (singles mode)
% 5 = Form the source image by using both singles and then dividing the
% counts by two
% 6 = Form the source image by using the average coordinates from both
% singles
% 7 = Form the source image by using the true coincidences (requires
% obtain_trues = true)
source_coordinates = 1;

if exist('f_osem','var') && ~exist('pz','var')
    pz = cell(95,1);
    pz{2} = f_osem;
end
img = check_algorithms(pz, algorithms, color_from_algo);
if isempty(img)
    return
end

if length(algorithms) + nnz(source) > 3
    hh = 3;
else
    hh = 1;
end
if length(algorithms) + nnz(source) < 4
    jj = min(3, length(algorithms) + nnz(source));
else
    jj = 3;
end
set(0,'units','pixels')
gg = get(0,'ScreenSize');
if jj > 4
    im_size = gg(4)/(2.5 + (jj - 4)/2);
else
    im_size = gg(4)/2.5;
end
figure
set(gcf, 'Position', [min(gg(3)/2-im_size/2*jj,gg(3)), min(gg(4)/2-im_size/2*hh, im_size * jj), gg(4), im_size * hh]);

if source == 1
    load([image_properties.machine_name '_Ideal_image_coordinates_' image_properties.name '_ASCII.mat'])
elseif source == 2
    load([image_properties.machine_name '_Ideal_image_coordinates_' image_properties.name '_LMF.mat'])
elseif source == 3
    load([image_properties.machine_name '_Ideal_image_coordinates_' image_properties.name '_ROOT.mat'])
end

if source > 0
    clim2 = [0 max(max(max(max(FOV(:,:,:)))))/color_scale];
end
clim = [0 max(max(max(max(img(:,:,:,end)))))/color_scale];
if v_plane == 2
    koko = size(img,1);
elseif v_plane == 3
    koko = size(img,2);
else
    koko = size(img,3);
end
if slice > koko
    error("Selected slice exceeds image size in the specified dimension/plane")
end
for kk = 1 : image_properties.n_time_steps
    for ll = 1 : length(algorithms)
        
        img = pz{algorithms(ll),kk};
        if v_plane == 2
            img = rot90(permute(img, [3 2 1 4]),2);
        elseif v_plane == 3
            img = permute(img, [1 3 2 4]);
        end
        subplot(hh, jj, ll)
            if color_from_algo == 0
                clim = [0 max(max(max(max(img(:,:,:,end)))))/color_scale];
                imagesc(img(:,:,slice,end),clim)
            else
                imagesc(img(:,:,slice,end),clim)
            end
        axis image
        title([char(algo_char(algorithms(ll))) ' time step = ' num2str(kk)])
    end
    if source > 0
        if source_coordinates == 1
            FOV = C{1,kk};
        elseif source_coordinates == 2
            FOV = C{2,kk};
        elseif source_coordinates == 3
            FOV = C{3,kk};
        elseif source_coordinates == 4
            FOV = C{4,kk};
        elseif source_coordinates == 5
            FOV = C{5,kk};
        elseif source_coordinates == 6
            FOV = C{6,kk};
        elseif source_coordinates == 7
            FOV = C{7};
        end
        if v_plane == 2
            FOV = permute(FOV, [1 3 2]);
        elseif v_plane == 3
            FOV = rot90(permute(FOV, [3 2 1]),2);
        end
        subplot(hh, jj, ll + 1)
        imagesc(FOV(:,:,slice),clim2)
        axis image
        title(['Original decay image, time step = ' num2str(kk)])
    end
    pause%(0.25)
    drawnow
end