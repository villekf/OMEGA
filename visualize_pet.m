%% Visualization for the PET reconstructions
% Each section has a visualization code for a different purpose
% Only a specific section should be run at a time
% This visualization file is no longer maintained. It is recommended to
% use volume3Dviewer instead. See help volume3Dviewer

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

% The below list can be generated (without line endings) with: algorithms_char(0)
% 1 = MLEM, 2 = OSEM, 3 = MRAMLA, 4 = RAMLA, 5 = ROSEM, 6 = RBI, 7 = DRAMA, 
% 8 = COSEM, 9 = ECOSEM, 10 = ACOSEM, 11 = Median Root (OSL-MLEM), 
% 12 = Median Root (OSL-OSEM), 13 = Median Root (BSREM), 14 = Median Root (MBSREM), 
% 15 = Median Root (ROSEM-MAP), 16 = Median Root (OSL-RBI), 17 = Median Root (OSL-COSEM), 
% 18 = Median Root (PKMA), 19 = Quadratic (OSL-MLEM), 20 = Quadratic (OSL-OSEM), 
% 21 = Quadratic (BSREM), 22 = Quadratic (MBSREM), 23 = Quadratic (ROSEM-MAP), 
% 24 = Quadratic (OSL-RBI), 25 = Quadratic (OSL-COSEM), 26 = Quadratic (PKMA), 
% 27 = Huber (OSL-MLEM), 28 = Huber (OSL-OSEM), 29 = Huber (BSREM), 30 = Huber (MBSREM), 
% 31 = Huber (ROSEM-MAP), 32 = Huber (OSL-RBI), 33 = Huber (OSL-COSEM), 
% 34 = Huber (PKMA), 35 = L-filter (OSL-MLEM), 36 = L-filter (OSL-OSEM), 
% 37 = L-filter (BSREM), 38 = L-filter (MBSREM), 39 = L-filter (ROSEM-MAP), 
% 40 = L-filter (OSL-RBI), 41 = L-filter (OSL-COSEM), 42 = L-filter (PKMA), 
% 43 = FIR Median Hybrid (OSL-MLEM), 44 = FIR Median Hybrid (OSL-OSEM), 
% 45 = FIR Median Hybrid (BSREM), 46 = FIR Median Hybrid (MBSREM), 
% 47 = FIR Median Hybrid (ROSEM-MAP), 48 = FIR Median Hybrid (OSL-RBI), 
% 49 = FIR Median Hybrid (OSL-COSEM), 50 = FIR Median Hybrid (PKMA), 
% 51 = Weighted mean (OSL-MLEM), 52 = Weighted mean (OSL-OSEM), 53 = Weighted mean (BSREM), 
% 54 = Weighted mean (MBSREM), 55 = Weighted mean (ROSEM-MAP), 56 = Weighted mean (OSL-RBI), 
% 57 = Weighted mean (OSL-COSEM), 58 = Weighted mean (PKMA), 
% 59 = Total Variation (OSL-MLEM), 60 = Total Variation (OSL-OSEM), 
% 61 = Total Variation (BSREM), 62 = Total Variation (MBSREM), 
% 63 = Total Variation (ROSEM-MAP), 64 = Total Variation (OSL-RBI), 
% 65 = Total Variation (OSL-COSEM), 66 = Total Variation (PKMA), 
% 67 = Anisotropic Diffusion (OSL-MLEM), 68 = Anisotropic Diffusion (OSL-OSEM), 
% 69 = Anisotropic Diffusion (BSREM), 70 = Anisotropic Diffusion (MBSREM), 
% 71 = Anisotropic Diffusion (ROSEM-MAP), 72 = Anisotropic Diffusion (OSL-RBI), 
% 73 = Anisotropic Diffusion (OSL-COSEM), 74 = Anisotropic Diffusion (PKMA), 
% 75 = Asymmetric Parallel Level Sets (OSL-MLEM), 76 = Asymmetric Parallel Level Sets (OSL-OSEM), 
% 77 = Asymmetric Parallel Level Sets (BSREM), 78 = Asymmetric Parallel Level Sets (MBSREM), 
% 79 = Asymmetric Parallel Level Sets (ROSEM-MAP), 80 = Asymmetric Parallel Level Sets (OSL-RBI), 
% 81 = Asymmetric Parallel Level Sets (OSL-COSEM), 82 = Asymmetric Parallel Level Sets (PKMA), 
% 83 = Total Generalized Variation (OSL-MLEM), 84 = Total Generalized Variation (OSL-OSEM), 
% 85 = Total Generalized Variation (BSREM), 86 = Total Generalized Variation (MBSREM), 
% 87 = Total Generalized Variation (ROSEM-MAP), 88 = Total Generalized Variation (OSL-RBI), 
% 89 = Total Generalized Variation (OSL-COSEM), 90 = Total Generalized Variation (PKMA), 
% 91 = Non-Local Means (OSL-MLEM), 92 = Non-Local Means (OSL-OSEM), 93 = Non-Local Means (BSREM), 
% 94 = Non-Local Means (MBSREM), 95 = Non-Local Means (ROSEM-MAP), 96 = Non-Local Means (OSL-RBI), 
% 97 = Non-Local Means (OSL-COSEM), 98 = Non-Local Means (PKMA), 99 = Custom (OSL-MLEM), 
% 100 = Custom (OSL-OSEM), 101 = Custom (BSREM), 102 = Custom (MBSREM), 103 = Custom (ROSEM-MAP), 
% 104 = Custom (OSL-RBI), 105 = Custom (OSL-COSEM), 106 = Custom (PKMA)

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


if length(algorithms) >= 4
    hh = 2;
else
    hh = 1;
end
if length(algorithms) < 4
    jj = min(3, length(algorithms));
elseif length(algorithms) == 4
    jj = 2;
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

clim = [0 max(max(max(max(img(:,:,2:end-1,end)))))/color_scale];
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
            clim = [0 max(max(max(max(img(:,:,2:end-1,end)))))/color_scale];
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
        clim = [0 max(max(max(max(img(:,:,2:end-1,end - ll + 1)))))/color_scale];
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


clim = [0 max(max(max(max(img(:,:,2:end-1,end)))))/color_scale];
clim2 = [0 max(max(max(max(FOV(:,:,2:end-1)))))/color_scale];
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
            clim = [0 max(max(max(max(img(:,:,2:end-1,end)))))/color_scale];
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
alpha_scale = 1;

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
set(gca, 'XLim', [0 128]);
set(gca, 'ZLim', [0 128]);
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
                clim = [0 max(max(max(max(img(:,:,2:end-1,end)))))/color_scale];
                imagesc(img(:,:,kk,end),clim)
            else
                imagesc(img(:,:,kk,end),clim)
            end
            title([char(algo_char(algorithms(ll))) ' slice = ' num2str(kk) ', transverse'])
        else
            subplot(hh, jj, 1 + jj*(ll - 1))
            if color_from_algo == 0
                clim = [0 max(max(max(max(img(:,:,2:end-1,end)))))/color_scale];
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
                clim = [0 max(max(max(max(img(:,:,2:end-1,end)))))/color_scale];
                imagesc(img(:,:,kk,end),clim)
            else
                imagesc(img(:,:,kk,end),clim)
            end
            title([char(algo_char(algorithms(ll))) ' slice = ' num2str(kk) ', frontal'])
        else
            img = rot90(permute(img, [3 2 1 4]),2);
            subplot(hh, jj, 2 + jj*(ll - 1))
            if color_from_algo == 0
                clim = [0 max(max(max(max(img(:,:,2:end-1,end)))))/color_scale];
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
                clim = [0 max(max(max(max(img(:,:,2:end-1,end)))))/color_scale];
                imagesc(img(:,:,kk,end),clim)
            else
                imagesc(img(:,:,kk,end),clim)
            end
            title([char(algo_char(algorithms(ll))) ' slice = ' num2str(kk) ', sagittal'])
        else
            img = permute(img, [3 1 2 4]);
            subplot(hh, jj, 3 + jj*(ll - 1))
            if color_from_algo == 0
                clim = [0 max(max(max(max(img(:,:,2:end-1,end)))))/color_scale];
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
% "true" image. This section requires image_properties to be loaded.

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

if v_plane == 2
    koko = size(img,1);
	if source > 0
		clim2 = [0 max(max(max(max(FOV(slice,:,:)))))/color_scale];
	end
	clim = [0 max(max(max(max(img(slice,:,:,end)))))/color_scale];
elseif v_plane == 3
    koko = size(img,2);
	if source > 0
		clim2 = [0 max(max(max(max(FOV(:,slice,:)))))/color_scale];
	end
	clim = [0 max(max(max(max(img(:,slice,:,end)))))/color_scale];
else
    koko = size(img,3);
	if source > 0
		clim2 = [0 max(max(max(max(FOV(:,:,slice)))))/color_scale];
	end
	clim = [0 max(max(max(max(img(:,:,slice,end)))))/color_scale];
end
if slice > koko
    error("Selected slice exceeds image size in the specified dimension/plane")
end
for kk = 1 : size(pz,2)
    for ll = 1 : length(algorithms)
        
        img = pz{algorithms(ll),kk};
        if v_plane == 2
            img = rot90(permute(img, [3 2 1 4]),2);
        elseif v_plane == 3
            img = permute(img, [1 3 2 4]);
        end
        subplot(hh, jj, ll)
            if color_from_algo == 0
                clim = [0 max(max(max(max(img(:,:,2:end-1,end)))))/color_scale];
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