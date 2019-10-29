%% Visualization for the PET reconstructions
% Each section has a visualization code for a different purpose
% Only a specific section should be run at a time

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

% Load the saved reconstruction and machine specific variables
image_properties = pz{end,1};



%% Visualize several reconstructions for one time step for all slices, last iterations

algo_char = algorithms_char();

% 1 = MLEM, 2 = OSEM, 3 = MRAMLA, 4 = RAMLA, 5 = ROSEM, 6 = RBI, 7 = DRAMA, 
% 8 = COSEM, 9 = ECOSEM, 10 = ACOSEM, 11 = MRP-OSL-OSEM, 12 = MRP-OSL-MLEM, 
% 13 = MRP-BSREM, 14 = MRP-MBSREM, 15 = MRP-ROSEM, 16 = MRP-RBI, 17 = MRP-OSL-COSEM, 
% 18 = QP (OSL-OSEM), 19 = QP (OSL-MLEM), 20 = QP (BSREM), 21 = QP (MBSREM), 
% 22 = QP (ROSEM), 23 = QP (RBI), 24 =  QP (OSL-COSEM), 25 = L-filter (OSL-OSEM), 
% 26 = L-filter (OSL-MLEM), 27 =  L-filter (BSREM), 28 =  L-filter (MBSREM), 
% 29 = L-filter (ROSEM), 30 = L-filter (RBI), 31 = L-filter (OSL-COSEM), 
% 32 = FMH (OSL-OSEM), 33 = FMH (OSL-MLEM), 34 = FMH (BSREM), 35 = FMH (MBSREM), 
% 36 = FMH (ROSEM), 37 = FMH (RBI), 38 = FMH (OSL-COSEM), 39 = Weighted mean (OSL-OSEM), 
% 40 = Weighted mean (OSL-MLEM), 41 = Weighted mean (BSREM), 42 = Weighted mean (MBSREM), 
% 43 = Weighted mean (ROSEM), 44 = Weighted mean (RBI), 45 = Weighted mean (OSL-COSEM), 
% 46 = Total variation (OSL-OSEM), 447 = Total variation (OSL-MLEM), 48 = Total variation (BSREM), 
% 49 = Total variation (MBSREM), 50 = Total variation (ROSEM), 51 = Total variation (RBI), 
% 52 = Total variation (OSL-COSEM), 53 = Anisotropic Diffusion (OSL-OSEM), 
% 54 = Anisotropic Diffusion (OSL-MLEM), 55 = Anisotropic Diffusion (BSREM), 
% 56 = Anisotropic Diffusion (MBSREM), 57 = Anisotropic Diffusion (ROSEM), 
% 58 = Anisotropic Diffusion (RBI), 59 = Anisotropic Diffusion (OSL-COSEM), 
% 60 = APLS (OSL-OSEM), 61 = APLS (OSL-MLEM), 62 = APLS (BSREM), 63 = APLS (MBSREM), 
% 64 = APLS (ROSEM), 65 = APLS (RBI), 66 = APLS (OSL-COSEM), 67 = TGV (OSL-OSEM), 
% 68 = TGV (OSL-MLEM), 69 = TGV (BSREM), 70 = TGV (MBSREM), 71 = TGV (ROSEM), 
% 72 = TGV (RBI), 73 = TGV (OSL-COSEM), 74 = NLM (OSL-OSEM), 75 =  NLM (OSL-MLEM), 
% 76 = NLM (BSREM), 77 = NLM (MBSREM), 78 = NLM (ROSEM), 79 = NLM (RBI), 80 = NLM (OSL-COSEM), 
% 81 = Custom prior (OSL-OSEM), 82 =  Custom prior (OSL-MLEM), 83 = Custom prior (BSREM), 
% 84 = Custom prior (MBSREM), 85 = Custom prior (ROSEM), 86 = Custom prior (RBI), 
% 87 = Custom prior (OSL-COSEM)
algorithms = [2];
% Use this value to scale the color scale in the image (higher values make
% low count areas brighter)
color_scale = 1;
% From which reconstruction should the color scale be taken
% NOTE: The numbering is according to the length of the above algorithms
% vector, e.g. if you have algorithms = [2, 4, 5] and color_from_algo = 2
% then the scale will be taken from RAMLA reconstruction (second element of
% algorithms) 
color_from_algo = 1;
% Visualization plane
% Choose the plane where the visualization takes place
% 1 = Transverse, 2 = Coronal/frontal, 3 = Sagittal
v_plane = 1;


img = check_algorithms(pz, algorithms, color_from_algo);
if isempty(img)
    return
end


if length(algorithms) > 2
    hh = 2;
else
    hh = 1;
end
if length(algorithms) > 1
    jj = max(2, round(length(algorithms) / 2));
else
    jj = 1;
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
        imagesc(img(:,:,kk,end),clim)
        axis image
        title([char(algo_char(algorithms(ll))) ' slice = ' num2str(kk)])
    end
    pause%(0.25)
    drawnow
end


%% Visualize N iterations of a single reconstruction for one time step for all slices

algo_char = algorithms_char();

algorithm = 2;
% Use this value to scale the color scale in the image (higher values make
% low count areas brighter)
color_scale = 1;
% Visualization plane
% Choose the plane where the visualization takes place
% 1 = Transverse, 2 = Coronal/frontal, 3 = Sagittal
v_plane = 1;

% How many iterations in the image 
% Initial values can be included as iterations as they are saved as such
N_iter = 2;

img = check_algorithms(pz, algorithms, color_from_algo);
if isempty(img)
    return
end

% clim = [0 max(max(max(max(img(:,:,:,end)))))/color_scale];

if N_iter > 2
    hh = 2;
else
    hh = 1;
end
if N_iter > 1
    jj = max(2, round(N_iter / 2));
else
    jj = 1;
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
        title([char(algo_char(algorithm)) ', iteration = ' num2str(size(img,4) - ll) ', slice = ' num2str(kk)])
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

% The source data was obtained from
% 1 = ASCII, 2 = LMF, 3 = ROOT
source = 1;

img = check_algorithms(pz, algorithms, color_from_algo);
if isempty(img)
    return
end

if length(algorithms) + 1 > 2
    hh = 2;
else
    hh = 1;
end
jj = max(2, round((length(algorithms)+1) / 2));
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
        imagesc(img(:,:,kk,end),clim)
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

algorithm = 2;
% The scale value for the pixel alpha values. Higher values will make the
% pixels more transparent, allowing areas of higher activity to be seen
% through background noise
alpha_scale = 10;

img = check_algorithms(pz, algorithms, color_from_algo);
if isempty(img)
    return
end

alpha_scaling = max(max(max(img(:,:,:,end)))) * alpha_scale;
alpha = img(:,:,:,end);
alpha = alpha./alpha_scaling;
alpha(alpha > 1) = 1;

vol3d('CData', img(:,:,:,end), 'Alpha', alpha);

%% Visualize simultanously all the views for n algorithms
% NOTE: Due to likely different number of slices, the smallers views will
% not be updated once they reach maximum number of slices

algo_char = algorithms_char();

algorithms = [1,2];
% algorithms = [2,3,4,5,6,7];
% algorithms = [8,9,10,11,12,13,14,15];
% Use this value to scale the color scale in the image (higher values make
% low count areas brighter)
color_scale = 1;
% From which reconstruction should the color scale be taken
% NOTE: The numbering is according to the length of the above algorithms
% vector, e.g. if you have algorithms = [2, 4, 5] and color_from_algo = 2
% then the scale will be taken from RAMLA reconstruction (second element of
% algorithms)
color_from_algo = 1;

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

clim = [0 max(max(max(max(img(:,:,:,end)))))/color_scale];
for kk = 1 : max([koko1, koko2, koko3])
    for ll = 1 : numel(algorithms)
        img = pz{algorithms(ll)};
        koko1 = size(img,1);
        koko2 = size(img,2);
        koko3 = size(img,3);
        if kk <= koko3
            subplot(hh, jj, 1 + jj*(ll - 1))
            imagesc(img(:,:,kk,end),clim)
            title([char(algo_char(algorithms(ll))) ' slice = ' num2str(kk) ', transverse'])
        else
            subplot(hh, jj, 1 + jj*(ll - 1))
            imagesc(img(:,:,koko3,end),clim)
            title([char(algo_char(algorithms(ll))) ' slice = ' num2str(koko1) ', transverse'])
        end
        axis image
        if kk <= koko1
            img = rot90(permute(img, [3 2 1 4]),2);
            subplot(hh, jj, 2 + jj*(ll - 1))
            imagesc(img(:,:,kk,end),clim)
            title([char(algo_char(algorithms(ll))) ' slice = ' num2str(kk) ', frontal'])
        else
            img = rot90(permute(img, [3 2 1 4]),2);
            subplot(hh, jj, 2 + jj*(ll - 1))
            imagesc(img(:,:,koko1,end),clim)
            title([char(algo_char(algorithms(ll))) ' slice = ' num2str(koko1) ', frontal'])
        end
        axis image
        if kk <= koko2
            img = permute(img, [1 3 2 4]);
            subplot(hh, jj, 3 + jj*(ll - 1))
            imagesc(img(:,:,kk,end),clim)
            title([char(algo_char(algorithms(ll))) ' slice = ' num2str(kk) ', sagittal'])
        else
            img = permute(img, [1 3 2 4]);
            subplot(hh, jj, 3 + jj*(ll - 1))
            imagesc(img(:,:,koko2,end),clim)
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

img = pz{algorithms(color_from_algo)};
for jj = 1:numel(algorithms)
    if isempty(pz{algorithms(jj)})
        warning('The current selected algorithm does not contain any estimes!')
        fprintf('The following are contained in the input array:\n')
        char_ar = algo_char(~cellfun(@isempty,pz(:,1)));
        loc = find(~cellfun(@isempty,pz(1:end-1,1)));
        for kk = 1 : nnz(~cellfun(@isempty,pz(:,1)))-1
            fprintf('Element %d: %s\n', loc(kk), char_ar{kk})
        end
        return
    end
end

if length(algorithms) + nnz(source) > 2
    hh = 2;
else
    hh = 1;
end
if length(algorithms) + nnz(source) > 1
    jj = max(2, round((length(algorithms) + nnz(source)) / 2));
else
    jj = 1;
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
        imagesc(img(:,:,slice,end),clim)
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