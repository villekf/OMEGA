function attenuation_datafile = attenuation_correction_factors(options)
%% ATTENUATION CORRECTION IMAGE
% This function creates an attenuation image from the Siemens Inveon
% .atn-file. First the data is reconstructed in a smaller size than the
% actual measurement, then scaled such that the majority of the values lie
% within the 122 keV tissue range. After that, the values are interpolated
% for 511 keV values. Lastly, the attenuation image is scaled to the actual
% image size. Outputs the name of the attenuation image datafile.

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


options.randoms_correction = false;
options.obtain_trues = false;
options.store_scatter = false;
options.store_randoms = false;
options.source = false;
options.partitions = 1;
options.root_source = false;
options.source_index1 = [];
% attenuation_factors = [];
% if options.blank
%     options.name = 'attenuation_blank_data';
%     if options.use_machine == 1
%         [options.file, options.fpath] = uigetfile('*.lst','Select Inveon list-mode blank datafile');
%         if isequal(options.file, 0)
%             error('No file was selected')
%         end
% %     elseif options.use_machine == 2
% %         [options.file, options.fpath] = uigetfile('*.scn','Select Inveon blank sinogram');
% %         if isequal(options.file, 0)
% %             error('No file was selected')
% %         end
% %         blank_sino = form_sinograms(options, false);
% %     elseif options.use_machine == 0
% %         options.fpath = uigetdir([], 'Select folder with GATE blank data');
% %         if isequal(options.fpath, 0)
% %             error('No folder was selected')
% %         end
%     end
%     if options.use_machine == 1
%         options.coincidences = load_data(options);
%         if ~options.use_raw_data
%             blank_sino = form_sinograms(options, false);
%             blank_sino(blank_sino == 0) = 1;
%         end
%     end
% end
% if options.transmission_scan
%     options.name = 'attenuation_transmission_data';
%     if options.use_machine == 1
%         [options.file, options.fpath] = uigetfile('*.lst','Select Inveon list-mode transmission datafile');
%         if isequal(options.file, 0)
%             error('No file was selected')
%         end
%     elseif options.use_machine == 2
        [options.file, options.fpath] = uigetfile('*.atn','Select Inveon attenuation file');
        if isequal(options.file, 0)
            error('No file was selected')
        end
        nimi = [options.fpath options.file];
        fid = fopen(nimi);
        attenuation_factors = fread(fid, options.Nang*options.Ndist*options.Nz, 'single',0,'l');
        attenuation_factors = log(reshape(attenuation_factors, options.Ndist, options.Nang, options.Nz));
%         attenuation_factors = permute(attenuation_factors, [2 1 3]);
%     else
%         options.fpath = uigetdir([], 'Select folder with GATE transmission data');
%         if isequal(options.fpath, 0)
%             error('No folder was selected')
%         end
%     end
%     if options.use_machine == 1
%         transmission_coincidences = load_data(options);
%     end
%     if (options.use_raw_data || options.precompute_all) && options.use_machine == 1
%         if options.precompute_all
%             if options.use_machine == 0
%                 load_string = [machine_name '_measurements_' name '_static_raw.mat'];
%             elseif options.use_machine == 1
%                 load_string = [machine_name '_measurements_' name '_static_raw_listmode.mat'];
%             end
%             load(load_string,'coincidences')
%             transmission_coincidences = coincidences;
%             clear coincidences
%         end
%         if options.use_machine == 0
%             load_string = [options.machine_name '_measurements_attenuation_blank_data_static_raw.mat'];
%         elseif options.use_machine == 1
%             load_string = [options.machine_name '_measurements_attenuation_blank_data_static_raw_listmode.mat'];
%         end
%         if exist(load_string,'file') == 2
%             load(load_string,'coincidences');
%             attenuation_factors = coincidences ./ transmission_coincidences;
%         else
%             error('Blank raw list-mode data file not found')
%         end
%     end
%     if ~options.use_raw_data || options.precompute_all || options.use_machine == 1
%         if options.precompute_all && options.use_machine < 2
%             if options.use_machine == 0
%                 load_string = [machine_name '_measurements_' name '_static.mat'];
%             elseif options.use_machine == 1
%                 load_string = [machine_name '_measurements_' name '_static_listmode.mat'];
%             end
%             load(load_string,'coincidences')
%             options.coincidences = coincidences;
%             clear coincidences
%         end
%         transmission_sino = form_sinograms(options, false);
%         transmission_sino(transmission_sino == 0) = 1;
%         if options.use_machine == 0
%             load_string = [options.machine_name '_attenuation_blank_data_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.NSlices) '_span' num2str(options.span) '.mat'];
%         elseif options.use_machine == 2
%             load_string = [options.machine_name '_attenuation_blank_data_sinograms_original_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.NSlices) '_span' num2str(options.span) '_machine_sinogram.mat'];
%         else
%             load_string = [options.machine_name '_attenuation_blank_data_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.NSlices) '_span' num2str(options.span) '_listmode.mat'];
%         end
%         if exist('blank_sino','var')
%             attenuation_factors = blank_sino ./ transmission_sino;
%         elseif exist(load_string,'file') == 2
%             load(load_string,'SinM');
%             SinM(SinM == 0) = 1;
%             attenuation_factors = SinM ./ transmission_sino;
%         else
%             error('Blank sinogram file not found')
%         end
%     end
    options.attenuation_phase = true;
    options.SinM = attenuation_factors;
    options.NSinos = size(attenuation_factors,3);
    options = set_implementations(options);
    options.osem = true;
    options.verbose = false;
    options.implementation = 4;
%     if options.subsets == 1
        options.subsets = 10;
%     end
    options.Niter = 6;
    options.use_raw_data = false;
    options.attenuation_correction = false;
    options.projector_type = 1;
    options.precompute_lor = false;
    options.Nx = options.Nx/2;
    options.Ny = options.Ny/2;
    options.x0 = ones(options.Nx, options.Ny, options.Nz);
    pz = reconstructions_main(options);
    for kk = 1 : length(pz)
        if isempty(pz{kk})
            continue
        else
            img = pz{kk}(:,:,:,end);
            break
        end
    end
    img = median_filter3d(img) / 10;
    attenuation_factors_ = attenuation122_to_511(img);
    attenuation_factors = ones(options.Nx * 2, options.Ny * 2, options.Nz);
    for kk = 1 : options.Nz
        attenuation_factors(:,:,kk) = imresize(attenuation_factors_(:,:,kk), [options.Nx * 2 options.Ny * 2]);
    end
    attenuation_datafile = [options.machine_name '_attenuation_coefficients_for_' options.name '_' ...
        num2str(options.Nx *2) 'x' num2str(options.Ny * 2) 'x' num2str(options.Nz) '.mat'];
    save(attenuation_datafile, 'attenuation_factors')
    disp('Attenuation correction image formed')
end