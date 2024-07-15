function attenuation_datafile = attenuation_correction_factors(options)
%% ATTENUATION CORRECTION IMAGE
% This function creates an attenuation image from the Siemens Inveon
% .atn-file. First the data is reconstructed in a smaller size than the
% actual measurement, then scaled such that the majority of the values lie
% within the 122 keV tissue range. After that, the values are interpolated
% for 511 keV values. Lastly, the attenuation image is scaled to the actual
% image size. Outputs the name of the attenuation image datafile.

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


options.randoms_correction = false;
options.obtain_trues = false;
options.store_scatter = false;
options.store_randoms = false;
options.source = false;
options.partitions = 1;
options.root_source = false;
options.source_index1 = [];
if options.CT_attenuation
    [options.file, options.fpath] = uigetfile('*.img','Select Inveon CT UMAP file');
    if isequal(options.file, 0)
        error('No file was selected')
    end
    nimi = [options.fpath options.file];
    fid = fopen(nimi);
    attenuation_factors = fread(fid, inf, 'single',0,'l');
    if numel(attenuation_factors) < options.Nx*options.Ny*options.Nz
        val = sqrt(numel(attenuation_factors) / options.Nz);
        attenuation_factors = reshape(attenuation_factors, val, val, options.Nz);
        if license('test', 'image_toolbox') || exist('imresize3','file') == 2
            attenuation_factors = imresize3(attenuation_factors, [options.Nx, options.Ny, options.Nz]);
        else
            apu = zeros(options.Nx, options.Ny, options.Nz);
            for kk = 1 : options.Nz
                apu(:,:,kk) = imresize(attenuation_factors(:,:,kk), [options.Nx, options.Ny]);
            end
            attenuation_factors = apu;
        end
    end
    attenuation_factors = rot90(reshape(attenuation_factors, options.Nx, options.Ny, options.Nz),-1);
    fclose(fid);
else
    [options.file, options.fpath] = uigetfile('*.atn','Select Inveon attenuation file');
    if isequal(options.file, 0)
        error('No file was selected')
    end
    nimi = [options.fpath options.file];
    fid = fopen(nimi);
    attenuation_factors = fread(fid, inf, 'single',0,'l');
    aSize = numel(attenuation_factors) / (options.Ndist * options.Nang);
    attenuation_factors = log(reshape(attenuation_factors, options.Ndist, options.Nang, aSize));
    attenuation_factors = attenuation_factors(:,:,1:options.rings*2 - 1);
    fclose(fid);
    options.attenuation_phase = true;
    options.SinM = attenuation_factors;
    options.NSinos = size(attenuation_factors,3);
    rec = recNames(0);
    for kk = 1 : numel(rec)
        options.(rec{kk}) = false;
    end
    options.OSEM = true;
    options.verbose = false;
    options.implementation = 4;
    %     if options.subsets == 1
    options.subsets = 10;
    %     end
    options.Niter = 6;
    options.save_iter = false;
    options.use_raw_data = false;
    options.attenuation_correction = false;
    options.projector_type = 1;
    options.precompute_lor = false;
    options.x0 = ones(options.Nx, options.Ny, options.Nz);
    options.use_psf = true;
    options.deblurring = false;
    options.corrections_during_reconstruction = false;
    options.normalization_correction = false;
    options.randoms_correction = false;
    options.scatter_correction = false;
    options.arc_correction = false;
    pz = reconstructions_main(options);
    img = pz(:,:,:,end);
%     img = median_filter3d(img) / 10;
    attenuation_factors = attenuation122_to_511(img) / 2;
%     attenuation_factors = ones(options.Nx * 2, options.Ny * 2, options.Nz);
%     for kk = 1 : options.Nz
%         attenuation_factors(:,:,kk) = imresize(attenuation_factors_(:,:,kk), [options.Nx * 2 options.Ny * 2]);
%     end
%     options.Nx = options.Nx * 2;
%     options.Ny = options.Ny * 2;
end
attenuation_datafile = [options.machine_name '_attenuation_coefficients_for_' options.name '_' ...
    num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '.mat'];
save(attenuation_datafile, 'attenuation_factors')
disp('Attenuation correction image formed')
end