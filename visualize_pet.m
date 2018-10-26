%% Visualization for the PET reconstructions
% Each section has a visualization code for a different purpose
% Only a specific section should be run at a time

image_properties = pz{end,1};

    algo_char = ["MLEM","OSEM","MRAMLA","RAMLA","ECOSEM","COSEM","ACOSEM","MRP-OSL",...
        "MRP-BSREM","Quadratic prior (OSL)","Quadratic prior (BSREM)", "L-filter (OSL)",...
        "L-filter (BSREM)", "FIR FMH (OSL)", "Weighted mean (OSL)"];

%% Visualize a single reconstruction for one time step for all slices, last iteration

% 1 = MLEM, 2 = OSEM, 3 = MRAMLA, 4 = RAMLA, 5 = ECOSEM, 6 = COSEM, 7 =
% ACOSEM, 8 = MRP-OSL, 9 = MRP-BSREM, 10 = Quadratic prior OSL, 11 =
% Quadratic prior BSREM 12 = L-filter OSL, 13 = L-filter BSREM, 14 = FIR
% FMH OSL, 15 = Weighted mean OSL 
algorithm = 2;
% Use this value to scale the color scale in the image (higher values make
% low count areas brighter)
color_scale = 1;

img = pz{algorithm};

clim = [0 max(max(max(max(img(:,:,:,end)))))/color_scale];
for kk = 1 : size(img,3)
    imagesc(img(:,:,kk,end),clim)
    axis image
    title([char(algo_char(algorithm)) ', slice = ' num2str(kk)])
    pause(0.25)
    drawnow
end

%% Visualize a single reconstruction for one time step for all slices, up to 8 iterations

% 1 = MLEM, 2 = OSEM, 3 = MRAMLA, 4 = RAMLA, 5 = ECOSEM, 6 = COSEM, 7 =
% ACOSEM, 8 = MRP-OSL, 9 = MRP-BSREM, 10 = Quadratic prior OSL, 11 =
% Quadratic prior BSREM 12 = L-filter OSL, 13 = L-filter BSREM, 14 = FIR
% FMH OSL, 15 = Weighted mean OSL 
algorithm = 2;
% Use this value to scale the color scale in the image (higher values make
% low count areas brighter)
color_scale = 1;

% How many iterations in the image (if larger than 8, then the last 8 will
% be used)
% Initial values can be included as iterations as they are saved as such
N_iter = 2;

img = pz{algorithm};

clim = [0 max(max(max(max(img(:,:,:,end)))))/color_scale];

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
figure
set(gcf, 'Position', [200, 200, 400 * jj, 400 * hh]);
for kk = 1 : size(img,3)
    for ll = 1 : min(N_iter, 8)
        subplot(hh, jj, ll)
        imagesc(img(:,:,kk,end - ll + 1),clim)
        axis image
        title([char(algo_char(algorithm)) ', iteration = ' num2str(ll) ', slice = ' num2str(kk)])
    end
    pause(0.25)
    drawnow
end

%% Visualize several reconstructions for one time step for all slices, last iterations

% 1 = MLEM, 2 = OSEM, 3 = MRAMLA, 4 = RAMLA, 5 = ECOSEM, 6 = COSEM, 7 =
% ACOSEM, 8 = MRP-OSL, 9 = MRP-BSREM, 10 = Quadratic prior OSL, 11 =
% Quadratic prior BSREM 12 = L-filter OSL, 13 = L-filter BSREM, 14 = FIR
% FMH OSL, 15 = Weighted mean OSL 
algorithms = [2,3,4,5,6,7];
% algorithms = [8,9,10,11,12,13,14,15];
% Use this value to scale the color scale in the image (higher values make
% low count areas brighter)
color_scale = 1;
% From according to which reconstruction should the color scale be taken
% NOTE: The numbering is according to the length of algorithms vector, e.g.
% if you have algorithms = [2, 4, 5] and color_from_algo = 2 then the scale
% will be taken from RAMLA reconstruction (second element of algorithms)
color_from_algo = 1;

img = pz{algorithms(color_from_algo)};

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
figure
set(gcf, 'Position', [0, -100, 400 * jj, 400 * hh]);

clim = [0 max(max(max(max(img(:,:,:,end)))))/color_scale];
for kk = 1 : size(img,3)
    for ll = 1 : length(algorithms)
        img = pz{algorithms(ll)};
        subplot(hh, jj, ll)
        imagesc(img(:,:,kk,end),clim)
        axis image
        title([char(algo_char(algorithms(ll))) ' slice = ' num2str(kk)])
    end
    pause(0.25)
    drawnow
end

%% Compare a single reconstruction for one time step for all slices, last iteration, with the source image obtained from GATE data
% NOTE: This is valid only for GATE data

% This can be used to compare the achieved reconstruction with the
% coordinates from which the actual radioactive decay happened
% I.e. these allow for the error checking of the reconstructed data

% 1 = MLEM, 2 = OSEM, 3 = MRAMLA, 4 = RAMLA, 5 = ECOSEM, 6 = COSEM, 7 =
% ACOSEM, 8 = MRP-OSL, 9 = MRP-BSREM, 10 = Quadratic prior OSL, 11 =
% Quadratic prior BSREM 12 = L-filter OSL, 13 = L-filter BSREM, 14 = FIR
% FMH OSL, 15 = Weighted mean OSL 
algorithm = 1;
% Use this value to scale the color scale in the image (higher values make
% low count areas brighter)
color_scale = 1;
% How is the source image formed?
% 1 = Form the source image by using only singles (coincidences) that
% originate from the very same location (source coordinates are the same)
% 2 = Form the source image by using only the first single
% 3 = Form the source image by using only the second single
% 4 = Form the source image by using both singles and then dividing the
% counts by two
% 5 = Form the source image by using the average coordinates from both
% singles
source_coordinates = 5;

% The source data was obtained from
% 1 = ASCII, 2 = LMF
source = 1;

img = pz{algorithm};

if source == 1
    load([image_properties.machine_name '_Ideal_image_coordinates_' image_properties.name '_ASCII.mat'])
elseif source == 2
    load([image_properties.machine_name '_Ideal_image_coordinates_' image_properties.name '_LMF.mat'])
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
end

clim = [0 max(max(max(max(img(:,:,:,end)))))/color_scale];
clim2 = [0 max(max(max(max(FOV(:,:,:)))))/color_scale];
for kk = 1 : size(img,3)
    subplot 121
    imagesc(img(:,:,kk,end),clim)
    axis image
    title([char(algo_char(algorithm)) ' slice = ' num2str(kk)])
    subplot 122
    imagesc(FOV(:,:,kk),clim2)
    axis image
    title(['Original decay image, slice = ' num2str(kk)])
    pause(0.25)
    drawnow
end

%% Compare several reconstructions for one time step for all slices, last iteration, with the source image obtained from GATE data
% NOTE: This is valid only for GATE data

% This can be used to compare the achieved reconstruction with the
% coordinates from which the actual radioactive decay happened
% I.e. these allow for the error checking of the reconstructed data

% 1 = MLEM, 2 = OSEM, 3 = MRAMLA, 4 = RAMLA, 5 = ECOSEM, 6 = COSEM, 7 =
% ACOSEM, 8 = MRP-OSL, 9 = MRP-BSREM, 10 = Quadratic prior OSL, 11 =
% Quadratic prior BSREM 12 = L-filter OSL, 13 = L-filter BSREM, 14 = FIR
% FMH OSL, 15 = Weighted mean OSL 
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
% 1 = Form the source image by using only singles (coincidences) that
% originate from the very same location (source coordinates are the same)
% 2 = Form the source image by using only the first single
% 3 = Form the source image by using only the second single
% 4 = Form the source image by using both singles and then dividing the
% counts by two
% 5 = Form the source image by using the average coordinates from both
% singles
source_coordinates = 1;

% The source data was obtained from
% 1 = ASCII, 2 = LMF, 3 = Root
source = 1;

img = pz{algorithms(color_from_algo)};

if length(algorithms) + 1 > 2
    hh = 2;
else
    hh = 1;
end
jj = max(2, round((length(algorithms)+1) / 2));
figure
set(gcf, 'Position', [200, 200, 400 * jj, 400 * hh]);

if source == 1
    load([image_properties.machine_name '_Ideal_image_coordinates_' image_properties.name '_ASCII.mat'])
elseif source == 2
    load([image_properties.machine_name '_Ideal_image_coordinates_' image_properties.name '_LMF.mat'])
elseif source == 3
    load([image_properties.machine_name '_Ideal_image_coordinates_' image_properties.name '_Root.mat'])
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
end


clim = [0 max(max(max(max(img(:,:,:,end)))))/color_scale];
clim2 = [0 max(max(max(max(FOV(:,:,:)))))/color_scale];
for kk = 1 : size(img,3)
    for ll = 1 : length(algorithms)
        img = pz{algorithms(ll)};
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
% Download: https://se.mathworks.com/matlabcentral/fileexchange/22940-vol3d-v2

% 1 = MLEM, 2 = OSEM, 3 = MRAMLA, 4 = RAMLA, 5 = ECOSEM, 6 = COSEM, 7 =
% ACOSEM, 8 = MRP-OSL, 9 = MRP-BSREM, 10 = Quadratic prior OSL, 11 =
% Quadratic prior BSREM 12 = L-filter OSL, 13 = L-filter BSREM, 14 = FIR
% FMH OSL, 15 = Weighted mean OSL 
algorithm = 2;
% The scale value for the pixel alpha values. Higher values will make the
% pixels more transparent, allowing areas of higher activity to be seen
% through background noise
alpha_scale = 2;

img = pz{algorithm};

alpha_scaling = max(max(max(img(:,:,:,end)))) * alpha_scale;
alpha = img(:,:,:,end);
alpha = alpha./alpha_scaling;
alpha(alpha > 1) = 1;

vol3d('CData', img(:,:,:,end), 'Alpha', alpha);