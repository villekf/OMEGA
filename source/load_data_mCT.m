function [varargout] = load_data_mCT(options)
%% Load Siemens Biograph mCT PET data
% Loads the 64-bit or 32-bit list-mode data into MATLAB/Octave.
%
% Outputs data formatted for pure, raw list-mode, data. This data can be
% later transformed into sinogram format or used as-is. Loads also the
% delayed coincidences if it is selected. Saves the corresponding sinograms
% unless raw list-mode data has been explicitly selected.
%
% The input struct options should include all the relevant information
% (machine properties, path, names, what are saved, etc.).
%
% The output data is saved in a mat-file in the current working directory.
%
% OUTPUT:
%   coincidences = A cell matrix containing the raw list-mode data for each
%   time step. The contents of the cell matrix are a vector containing the
%   coincidences of each unique LOR (no duplicates, e.g. coincidences
%   from detector 1 to 2 and 2 to 1 are summed up)
%   delayed_coincidences = Same as above, but for delayed coincidences
%
% See also form_sinograms, initial_michelogram, load_data, load_data_Vision

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

if nargout > 2
    error('Too many output arguments')
end

folder = fileparts(which('load_data_mCT.m'));
folder = strrep(folder, 'source','mat-files/');
folder = strrep(folder, '\','/');
% 
% Nx = uint32(options.Nx);
% Ny = uint32(options.Ny);
% Nz = uint32(options.Nz);

disp('Beginning data load')

if options.use_machine == 1 || options.use_machine == 3
    
    [file, fpath] = uigetfile({'*.ptd;*.l;*.lst'},'Select mCT list-mode datafile');
    if isequal(file, 0)
        error('No file was selected!')
    end
    
    if options.verbose
        tic
    end
    
    machine_name = options.machine_name;
    name = options.name;
    tot_time = options.tot_time;
%     detectors = options.detectors;
    
    if isinf(tot_time)
        tot_time = 1e9;
    end
    
    if options.use_machine == 1
        type = int32(0);
    else
        type = int32(1);
    end
    
    partitions = options.partitions;
    
    loppu = options.end * 1e3;
    alku = options.start * 1e3;
    if isinf(loppu)
        loppu = 1e9;
    end
    vali = (loppu - alku)/options.partitions;
    
    nimi = [fpath file];
    f = dir(nimi);
    v_size = uint64(f.bytes / 8);
    if v_size == 0
        error('File is empty!')
    end
    
    if options.use_machine == 1
        % 64-bit list-mode data
        
        if options.randoms_correction
            variableList = {'SinM','delayed_coincidences'};
        else
            variableList = {'SinM'};
        end
        if partitions > 1
            save_string_raw = [machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_' num2str(loppu - alku) 's_raw_listmode.mat'];
        else
            save_string_raw = [machine_name '_measurements_' name '_static_raw_listmode.mat'];
        end
        
        if options.verbose
            disp('Data load started')
        end
        
        % Load the data
        [LL1, LL2, DD1, DD2, tpoints] = list2sinogram(nimi, uint64(vali), uint64(alku), uint64(loppu), logical(options.randoms_correction), type, v_size);
        clear mex
        
        %     TOF(TOF == 0) = [];
        %     TOF = TOF - 1;
        
        
        if partitions == 1
            LL1(LL1 == 0) = [];
            LL2(LL2 == 0) = [];
            LL3 = LL2;
            LL2(LL2 > LL1) = LL1(LL2 > LL1);
            LL1(LL3 > LL1) = LL3(LL3 > LL1);
            clear LL3
            prompts = accumarray([LL1 LL2],1,[options.detectors options.detectors],[],[],true);
            if options.randoms_correction
                DD1(DD1 == 0) = [];
                DD2(DD2 == 0) = [];
                DD3 = DD2;
                DD2(DD2 > DD1) = DD1(DD2 > DD1);
                DD1(DD3 > DD1) = DD3(DD3 > DD1);
                clear DD3
                delays = accumarray([DD1 DD2],1,[options.detectors options.detectors],[],[],true);
            end
            clear LL1 LL2 DD1 DD2
        else
            prompts = cell(partitions,1);
            if options.randoms_correction
                delays = cell(partitions,1);
            end
            ll = 1;
            for kk = 1 : partitions
                apu1 = LL1(ll:tpoints(kk));
                apu1(apu1 == 0) = [];
                apu2 = LL2(ll:tpoints(kk));
                apu2(apu2 == 0) = [];
                apu3 = apu2;
                apu2(apu2 > apu1) = apu1(apu2 > apu1);
                apu1(apu3 > apu1) = apu3(apu3 > apu1);
                clear apu3
                prompts{kk} = accumarray([apu1 apu2],1,[options.detectors options.detectors],[],[],true);
                if options.randoms_correction
                    apu1 = DD1(ll:tpoints(kk));
                    apu1(apu1 == 0) = [];
                    apu2 = DD2(ll:tpoints(kk));
                    apu2(apu2 == 0) = [];
                    apu3 = apu2;
                    apu2(apu2 > apu1) = apu1(apu2 > apu1);
                    apu1(apu3 > apu1) = apu3(apu3 > apu1);
                    clear apu3
                    delays{kk} = accumarray([apu1 apu2],1,[options.detectors options.detectors],[],[],true);
                end
                ll = 1 + tpoints(kk);
            end
            clear LL1 LL2 DD1 DD2 apu1 apu2
        end
        %% Save the measurement data
        
        coincidences = cell(partitions);
        if options.randoms_correction
            delayed_coincidences = cell(partitions,1);
        end
        
        tot_counts = 0;
        if options.randoms_correction
            tot_delayed = 0;
        end
            
        if ~options.use_raw_data
            if options.randoms_correction
                form_sinograms(options, false, prompts, delays);
            else
                form_sinograms(options, false, prompts);
            end
        end
        for llo=1:partitions
            
            if iscell(prompts)
                prompts{llo} = (prompts{llo}(tril(true(size(prompts{llo})), 0)));
                if options.verbose
                    counts = full(sum(sum(prompts{llo})));
                end
                coincidences{llo, 1} = prompts{llo};
            else
                prompts = (prompts(tril(true(size(prompts)), 0)));
                if options.verbose
                    counts = full(sum(sum(prompts)));
                end
                coincidences{llo, 1} = prompts;
            end
            
            if options.verbose
                disp(['Total number of prompts at time point ' num2str(llo) ' is ' num2str(counts) '.'])
                tot_counts = tot_counts + counts;
            end
            
            if options.randoms_correction
                if ~iscell(prompts)
                    delays = (delays(tril(true(size(delays)), 0)));
                    delayed_coincidences{partitions, 1} = delays;
                    if options.verbose
                        disp(['Total number of delayed coincidences at time point ' num2str(llo) ' is ' num2str(full(sum(sum(delays)))) '.'])
                        tot_delayed = tot_delayed + full(sum(sum(delays)));
                    end
                else
                    delays{llo} = (delays{llo}(tril(true(size(delays{llo})), 0)));
                    delayed_coincidences{partitions, 1} = delays{llo};
                    if options.verbose
                        disp(['Total number of delayed coincidences at time point ' num2str(llo) ' is ' num2str(full(sum(sum(delays{llo})))) '.'])
                        tot_delayed = tot_delayed + full(sum(sum(delays{llo})));
                    end
                end
            end
        end
        if partitions == 1
            if exist('OCTAVE_VERSION', 'builtin') == 0
                save(save_string_raw, 'coincidences', '-v7.3')
            else
                save(save_string_raw, 'coincidences', '-v7')
            end
            for variableIndex = 1:length(variableList)
                if exist(variableList{variableIndex},'var')
                    if exist('OCTAVE_VERSION', 'builtin') == 0
                        save(save_string_raw,variableList{variableIndex},'-append')
                    else
                        save(save_string_raw,variableList{variableIndex},'-append', '-v7')
                    end
                end
            end
        else
            if options.verbose
                disp(['Total number of prompts at all time points is ' num2str(tot_counts) '.'])
                if options.randoms_correction
                    disp(['Total number of delayed coincidences at all time points is ' num2str(tot_delayed) '.'])
                end
            end
            if exist('OCTAVE_VERSION', 'builtin') == 0
                save(save_string_raw, 'coincidences', '-v7.3')
            else
                save(save_string_raw, 'coincidences', '-v7')
            end
            for variableIndex = 1:length(variableList)
                if exist(variableList{variableIndex},'var')
                    if exist('OCTAVE_VERSION', 'builtin') == 0
                        save(save_string_raw,variableList{variableIndex},'-append')
                    else
                        save(save_string_raw,variableList{variableIndex},'-append', '-v7')
                    end
                end
            end
        end
        %     end
        clear LL1 DD1
        if nargout >= 1
            varargout{1} = coincidences;
        end
        if nargout >= 2
            varargout{2} = delayed_coincidences;
        end
        if options.verbose
            disp('Measurements loaded and saved')
            toc
        end
    elseif options.use_machine == 3
        % 32-bit list-mode data
        
        if options.randoms_correction
            variableList = {'SinM','delayed_coincidences','appliedCorrections'};
        else
            variableList = {'SinM','appliedCorrections'};
        end
        appliedCorrections.gapFilling = false;
        appliedCorrections.normalization = false;
        appliedCorrections.randoms = false;
        appliedCorrections.scatter = false;
        appliedCorrections.globalFactor = options.global_correction_factor;
        RandProp.variance_reduction = false;
        RandProp.smoothing = false;
        if partitions == 1
            save_string = [options.machine_name '_' options.name '_sinogram_original_static_' num2str(options.Ndist) 'x' num2str(options.Nang)...
                'x' num2str(options.TotSinos) '_span' num2str(options.span) '_listmode_sinogram.mat'];
        else
            save_string = [machine_name '_' name '_sinograms_combined_' num2str(partitions) 'timepoints_for_total_of_' num2str(tot_time) 's_' ...
                num2str(Ndist) 'x' num2str(Nang) 'x' num2str(NSlices) '_span' num2str(span) '_listmode_sinogram.mat'];
        end
        
        [raw_SinM, delayed_coincidences] = list2sinogram(nimi, uint64(vali), uint64(alku), uint64(loppu), logical(options.randoms_correction), type);
        clear mex
        
        if partitions == 1
            raw_SinM = permute(reshape(raw_SinM, options.Nang, options.Ndist, options.TotSinos), [2 1 3 4]);
            if options.randoms_correction
                delayed_coincidences = permute(reshape(delayed_coincidences, options.Nang, options.Ndist, options.TotSinos), [2 1 3 4]);
            end
            SinM = single(raw_SinM);
            if (~options.corrections_during_reconstruction && (options.randoms_correction || options.normalization_correction)) || ...
                    options.fill_sinogram_gaps && (sum(options.pseudot) > 0 || options.det_per_ring < options.det_w_pseudo)
                if ~options.corrections_during_reconstruction && options.randoms_correction
                    SinM = SinM - single(delayed_coincidences);
                    appliedCorrections.randoms = true;
                end
                if options.normalization_correction && ~options.corrections_during_reconstruction && options.use_user_normalization
                    [n_file, fpath] = uigetfile({'*.n;*.ptd'},'Select normalization datafile');
                    if isequal(n_file, 0)
                        error('No file was selected')
                    end
                    f_path = [fpath n_file];
                    normalization = loadNormalizationBiograph(f_path, options);
%                     normalization = circshift(normalization,1);
                    SinM = single(SinM) .* normalization;
                    appliedCorrections.normalization = true;
                elseif options.normalization_correction && ~options.corrections_during_reconstruction && ~options.use_user_normalization
                    norm_file = [folder options.machine_name '_normalization_' num2str(options.Ndist) 'x' num2str(options.Nang) '_span' num2str(options.span) '.mat'];
                    if exist(norm_file, 'file') == 2
                        normalization = loadStructFromFile(norm_file,'normalization');
                    else
                        normalization = normalization_coefficients(options);
                    end
                    normalization = circshift(normalization,1);
                    SinM = single(SinM) .* normalization;
                    appliedCorrections.normalization = true;
                end
                if options.fill_sinogram_gaps && (sum(options.pseudot) > 0 || options.det_per_ring < options.det_w_pseudo)
                    [~, ~, xp, yp] = detector_coordinates(options);
                    gaps = [];
                    llo = 1;
                    [SinM, gaps] = gapFilling(options, SinM, xp, yp, llo, gaps);
                    appliedCorrections.gapFilling = true;
                end
            end
        else
            apu = reshape(raw_SinM, options.Nang, options.Ndist, options.TotSinos, options.partitions);
            raw_SinM = cell(options.partitions,1);
            if ~options.corrections_during_reconstruction && (options.randoms_correction || options.normalization_correction || ...
                    options.fill_sinogram_gaps && (sum(options.pseudot) > 0 || options.det_per_ring < options.det_w_pseudo))
                SinM = cell(options.partitions,1);
            end
            if options.randoms_correction
                apu_d = reshape(delayed_coincidences, options.Nang, options.Ndist, options.TotSinos, options.partitions);
                delayed_coincidences = cell(options.partitions,1);
            end
            for kk = 1 : options.partitions
                raw_SinM{kk} = permute(apu(:,:,:,kk), [2 1 3 4]);
                if options.randoms_correction
                    delayed_coincidences{kk} = permute(reshape(apu_d, options.Nang, options.Ndist, options.TotSinos), [2 1 3 4]);
                end
                if (~options.corrections_during_reconstruction && (options.randoms_correction || options.normalization_correction)) || ...
                        options.fill_sinogram_gaps && (sum(options.pseudot) > 0 || options.det_per_ring < options.det_w_pseudo)
                    Sino = single(raw_SinM);
                    if ~options.corrections_during_reconstruction && options.randoms_correction
                        Sino = Sino - single(delayed_coincidences);
                        appliedCorrections.randoms = true;
                    end
                    if options.normalization_correction && ~options.corrections_during_reconstruction
                        [n_file, fpath] = uigetfile({'*.n;*.ptd'},'Select normalization datafile');
                        if isequal(n_file, 0)
                            error('No file was selected')
                        end
                        f_path = [fpath n_file];
                        normalization = loadNormalizationBiograph(f_path, options);
                        Sino = single(Sino) .* normalization;
                        appliedCorrections.normalization = true;
                    end
                    if options.fill_sinogram_gaps && (sum(options.pseudot) > 0 || options.det_per_ring < options.det_w_pseudo)
                        [~, ~, xp, yp] = detector_coordinates(options);
                        gaps = [];
                        [Sino, gaps] = gapFilling(options, Sino, xp, yp, 1, gaps);
                        appliedCorrections.gapFilling = true;
                    end
                    SinM{kk} = Sino;
                end
            end
        end
        
        if partitions == 1
            if options.verbose
                disp(['Total number of prompts is ' num2str(sum(SinM(:))) '.'])
                if options.randoms_correction
                    disp(['Total number of delayed coincidences is ' num2str(sum(delayed_coincidences(:))) '.'])
                end
            end
        else
            if options.verbose
                disp(['Total number of prompts at all time points is ' num2str(sum(apu(:))) '.'])
                if options.randoms_correction
                    disp(['Total number of delayed coincidences at all time points is ' num2str(sum(apu_d(:))) '.'])
                end
            end
        end
        if exist('OCTAVE_VERSION','builtin') == 0
            save(save_string,'raw_SinM','-v7.3')
            for variableIndex = 1:length(variableList)
                if exist(variableList{variableIndex},'var')
                    save(save_string,variableList{variableIndex},'-append')
                end
            end
        else
            save(save_string,'raw_SinM','-v7')
            for variableIndex = 1:length(variableList)
                if exist(variableList{variableIndex},'var')
                    save(save_string,variableList{variableIndex},'-append','-v7')
                end
            end
        end
        if nargout >= 1
            varargout{1} = SinM;
        end
        if nargout >= 2
            varargout{2} = delayed_coincidences;
        end
    end
else
    % Machine created sinogram
    if options.randoms_correction
        variableList = {'delayed_coincidences','appliedCorrections'};
    else
        variableList = {'appliedCorrections'};
    end
    appliedCorrections.gapFilling = false;
    appliedCorrections.normalization = false;
    appliedCorrections.randoms = false;
    appliedCorrections.scatter = false;
    appliedCorrections.globalFactor = options.global_correction_factor;
    tot_time = options.tot_time;
    
    if isinf(tot_time)
        tot_time = 1e9;
    end
    if options.partitions == 1
        save_string = [options.machine_name '_' options.name '_sinogram_original_static_' num2str(options.Ndist) 'x' num2str(options.Nang)...
            'x' num2str(options.TotSinos) '_span' num2str(options.span) '_machine_sinogram.mat'];
    else
        save_string = [options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_' num2str(tot_time) 's_' ...
            num2str(Ndist) 'x' num2str(Nang) 'x' num2str(NSlices) '_span' num2str(span) '_machine_sinogram.mat'];
    end
    [file, fpath] = uigetfile({'*.ptd;*.s'},'Select mCT sinogram datafile');
    if isequal(file, 0)
        error('No file was selected!')
    end
    nimi = [fpath file];
    fid = fopen(nimi);
    if any(strfind(file, '.ptd'))
        Sino = fread(fid, inf, 'char=>char');
        Sinot = Sino(end-20000:end)';
        idx = strfind(Sinot,'data offset in bytes[1]:=');
        idx2 = strfind(Sinot,'scan data type description[2]:=prompts') - 3;
        characteri = 'data offset in bytes[1]:=';
        offset_alku = str2double(Sinot(idx + length(characteri) : idx2));
        characteri = '%number of TOF time bins:=';
        idx = strfind(Sinot,'%number of TOF time bins:=');
        idx2 = strfind(Sinot,'%TOF mashing factor:=') - 3;
        TOF_timebins = str2double(Sinot(idx + length(characteri) : idx2)) + 1;
        offset_loppu = length(Sino) - offset_alku - options.Nang*options.Ndist*options.TotSinos*options.partitions * TOF_timebins * 2;
        fseek(fid, offset_alku, 0);
        Sino = fread(fid, inf, 'int16=>int16');
        Sino(end - offset_loppu + 1 : end) = [];
        Sino = reshape(Sino, options.Ndist, options.Nang, options.TotSinos, TOF_timebins, options.partitions);
        if TOF_timebins > 1
            delayed_coincidences = squeeze(Sino(:,:,:,end,:));
            Sino = Sino(:,:,:,1:end-1,:);
        end
        Sino = squeeze(sum(Sino,4));
    else
        Sino = fread(fid, inf, 'int16=>int16');
        if any(Sino < 0)
            fid = fopen(nimi);
            Sino = fread(fid, inf, 'single=>single');
        end
        TOF_timebins = length(Sino) / (options.Ndist*options.Nang*options.TotSinos*options.partitions);
        Sino = reshape(Sino, options.Ndist, options.Nang, options.TotSinos, TOF_timebins, options.partitions);
        if TOF_timebins > 1
            delayed_coincidences = squeeze(Sino(:,:,:,end,:));
            Sino = Sino(:,:,:,1:end-1,:);
        end
        Sino = squeeze(sum(Sino,4));
        Sino(Sino < 0) = 0;
    end
    fclose(fid);
    if options.partitions > 1
        raw_SinM = cell(options.partitions,1);
        for kk = 1 : options.partitions
            raw_SinM{kk} = Sino(:,:,:,kk);
        end
    else
        raw_SinM = Sino;
    end
    if exist('OCTAVE_VERSION','builtin') == 0
        save(save_string,'raw_SinM','-v7.3')
        for variableIndex = 1:length(variableList)
            if exist(variableList{variableIndex},'var')
                save(save_string,variableList{variableIndex},'-append')
            end
        end
    else
        save(save_string,'raw_SinM','-v7')
        for variableIndex = 1:length(variableList)
            if exist(variableList{variableIndex},'var')
                save(save_string,variableList{variableIndex},'-append','-v7')
            end
        end
    end
    if nargout >= 1
        varargout{1} = raw_SinM;
    end
    if nargout >= 2
        varargout{2} = delayed_coincidences;
    end
    
end
disp('Data load complete')
