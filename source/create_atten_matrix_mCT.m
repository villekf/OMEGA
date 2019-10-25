function attenuation_factors = create_atten_matrix_mCT(options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kansio, jossa vaimennuskuvat
% path = 'D:\Samulin kansio\University of Eastern Finland\Ville-Veikko Wettenhovi - Samulin kesatyo\vaimennuskuvat\CT_AC_SYDÄN_B19F_0002\CT_AC_SYDÄN_B19F_0002';
% path = 'C:\Users\samsu\University of Eastern Finland\Ville-Veikko Wettenhovi - Samulin kesatyo\vaimennuskuvat\CT_AC_SYDÄN_B19F_0002\CT_AC_SYDäN_B19F_0002';
% addpath(path);
fpath = uigetdir([], 'Select folder with Biograph mCT CT-images');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fnames = dir([fpath '\*.IMA']);
numfids = length(fnames);
info = dicominfo([fnames(1).folder '/' fnames(1).name]);
W = info.Width;    %width of image
H = info.Height;    %height of image
atten_image = zeros(H,W,length(numfids));
for k = 1:numfids
  atten_image(:,:,k) = dicomread([fnames(k).folder '/' fnames(k).name]);
end
ct_img_pos = info.ImagePositionPatient;
pet_img_pos = [0 141 2473.51];
z_difference = pet_img_pos(3)-ct_img_pos(3);
pix_space = info.PixelSpacing;
slice_thick = info.SliceThickness;
Pet_pixels = options.Nx;
Nslices = options.Nz;
pix_space_pet = options.FOVa_x/options.Nx;
slice_lkm_ct = floor(options.axial_fov/slice_thick) + 2*(ceil(((numfids*slice_thick - options.axial_fov)/2)/slice_thick) - floor(((numfids*slice_thick - options.axial_fov)/2)/slice_thick));
atten_image = atten_image(:,:,slice_lkm_ct+1:222-slice_lkm_ct)-1000;
alkukohta = -3*slice_lkm_ct/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%skaalataan CT-kuvat sopivan kokoiseksi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resized = zeros(W,H,Nslices);
for i=1:W
    resized(:,i,:) = reshape(imresize(reshape(atten_image(:,i,:),W,slice_lkm_ct),[W Nslices]),W,1,Nslices);
end
% diag = options.FOVa_x * sqrt(2);         %FOV:n diagonaalin pituus
% scaling_fact = diag/(pix_space(1)*W);
% pizel_width = diag/Pet_pixels;
pixel_lkm_ct = ceil(double(W) - ((double(W)*pix_space(1)) - options.FOVa_x) / pix_space(1));
dd = floor((double(W)-pixel_lkm_ct)/2);
resized = resized(dd+2:end-dd,dd+2:end-dd,:);
reresized = zeros(Pet_pixels,Pet_pixels,Nslices);
for i=1:Nslices
    reresized(:,:,i) = reshape(imresize(reshape(resized(:,:,i),pixel_lkm_ct,pixel_lkm_ct),[Pet_pixels Pet_pixels]),Pet_pixels,Pet_pixels,1);
end
% save(['CT_kuvat_resized_' num2str(Pet_pixels) '_' num2str(Pet_pixels) '_' num2str(Nslices) '.mat'],'reresized');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% vaimennuskertoimien skaalaaminen %%%
%%% trilinear interpolation          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Artikkeli "Accuracy of CT-Based Attenuation Correction...
% lineaarisen sovituksen kertoimet (kun KVP=100)!
if info.KVP == 100
    a = [9.3 4 0.5]*10^(-5);
    b = [0.093 0.093 0.128];
elseif info.KVP == 80
    a = [9.3 3.28 0.41]*10^(-5);
    b = [0.093 0.093 0.122];
elseif info.KVP == 120
    a = [9.3 4.71 0.589]*10^(-5);
    b = [0.093 0.093 0.134];
elseif info.KVP == 140
    a = [9.3 5.59 0.698]*10^(-5);
    b = [0.093 0.093 0.142];
else
    error('Unsupported kVp')
end
%risteyskohdat:
x(1,:) = [-1000,b(1)-1000*a(1)];
x(2,:) = [0,b(2)];
x(3,:) = [1000,b(2)+a(2)*1000];
x(4,:) = [3000,b(2)+a(2)*1000+a(3)*2000];
tarkkuus = 0.1;
inter = interp1(x(:,1),x(:,2),-1000:tarkkuus:3000,'linear');
vali = -1000:tarkkuus:3000;
plot(vali,inter)
%asetetaan sopivasti
attenuation_factors = zeros(size(reresized));
for ii=1:Nslices
%     ii
    temppi=reresized(:,:,ii);
    temppi=temppi(:);
    tic
    for ll=1:Pet_pixels^2
        apu=temppi(ll);
        [~,idx]=min(abs(vali-apu));
        apu=inter(idx);
        temppi(ll)=apu;
    end
    toc
    attenuation_factors(:,:,ii)=reshape(temppi,Pet_pixels,Pet_pixels);  
end
imagesc(attenuation_factors(:,:,55))
save([options.machine_name '_attenuation_coefficients_for_' options.name '.mat'], attenuation_factors)