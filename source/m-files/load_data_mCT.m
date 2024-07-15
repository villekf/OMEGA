function [varargout] = load_data_mCT(options)
%% Load Siemens Biograph mCT PET data
% Loads the 64-bit or 32-bit list-mode data into MATLAB/Octave.
%
% Outputs the sinogram data and optionally the detector coordinates for
% list-mode reconstruction. Loads also the delayed coincidences if it is
% selected.
%
% The input struct options should include all the relevant information
% (scanner properties, path, names, what are saved, etc.).
%
% The output sinogram data is saved in a mat-file in the current working
% directory.
%
% OUTPUT:
%   Sino = A matrix containing the sinograms. If TOF is selected contains
%   the TOF bins as well. With dynamic data, the sinogram is a cell matrix
%   instead where each cell element is a single time step.
%   SinoD = Same as above, but for delayed coincidences. Note that the
%   structure is otherwise the same, but TOF bins are never included.
%   Coord = Detector coordinates for each event (optional), [xs, ys, zs,
%   xd, yd, zd] where s refers to source (or in PET detector 1) and d to
%   detector (or detector 2). x/y/z are the corresponding axes.
%   Rcoord = Delayed coincidence detector coordinates for each event
%   (optional), same structure as above.
%
% See also form_sinograms, initial_michelogram, load_data, load_data_Vision

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020-2024 Ville-Veikko Wettenhovi
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

if nargout > 4
    error('Too many output arguments')
end
if nargout >= 3
    store_coordinates = true;
else
    store_coordinates = false;
end

folder = fileparts(which('load_data_mCT.m'));
folder = [folder(1:end-14), 'mat-files/'];
folder = strrep(folder, '\','/');
%
% Nx = uint32(options.Nx);
% Ny = uint32(options.Ny);
% Nz = uint32(options.Nz);

if ~isfield(options,'TOF_bins') || options.TOF_bins == 0
    options.TOF_bins = 1;
end

disp('Beginning data load')

alku = double(options.start);
loppu = options.end;
if isinf(loppu)
    loppu = 1e9;
end
vali = double((loppu - alku)/options.partitions);

if numel(options.partitions) == 1 && options.partitions > 1
    if isinf(options.end)
        error('End time is infinity, but more than one time-step selected. Use either one time-step or input a finite end time.')
    end
    options.partitions = repmat(vali, options.partitions, 1);
end

Nt = numel(options.partitions);

if options.use_machine == 1 || options.use_machine == 3

    if exist(options.fpath,'file') ~= 2
        [options.file, options.fpath] = uigetfile({'*.ptd;*.l;*.lst'},'Select mCT list-mode datafile');
        if isequal(options.file, 0)
            error('No file was selected!')
        end
        nimi = [options.fpath options.file];
    else
        nimi = [options.fpath];
    end
    [~,~,ext] = fileparts(nimi);
    if strcmp(ext,'ptd')
        usePTD = true;
    else
        usePTD = false;
    end

    if options.verbose
        tic
    end

    machine_name = options.machine_name;
    name = options.name;
    tot_time = options.tot_time;
    %     detectors = options.detectors;
    totSinos = options.TotSinos;
    if options.span == 1
        totSinos = options.rings^2;
    end
    sinoSize = uint64(options.Ndist * options.Nang * totSinos);
    pseudoD = options.det_per_ring < options.det_w_pseudo;
    pseudoR = sum(options.pseudot) > 0;

    if isinf(tot_time)
        tot_time = 1e9;
    end

    if options.use_machine == 1
        type = int32(0);
    else
        type = int32(1);
    end

    partitions = options.partitions;


    f = dir(nimi);
    v_size = uint64(f.bytes / 8);
    if v_size == 0
        error('File is empty!')
    end

    if options.use_machine == 1
        % 64-bit list-mode data

        % if options.randoms_correction
        %     variableList = {'SinM','delayed_coincidences'};
        % else
        %     variableList = {'SinM'};
        % end
        % if partitions > 1
        %     save_string_raw = [machine_name '_measurements_' name '_' num2str(Nt) 'timepoints_for_total_of_' num2str(loppu - alku) 's_raw_listmode.mat'];
        % else
        %     save_string_raw = [machine_name '_measurements_' name '_static_raw_listmode.mat'];
        % end

        if options.verbose
            disp('Data load started')
        end

        % Load the data
        [raw_SinM, SinD, LL1, LL2, tpoints, DD1, DD2] = list2sinogram(nimi, partitions, (alku), (loppu), logical(options.randoms_correction), v_size, store_coordinates, ...
            uint32(options.det_w_pseudo), sinoSize, uint32(options.Ndist), uint32(options.Nang), uint32(options.ring_difference), uint32(options.span), uint32(cumsum(options.segment_table)), ...
            Nt, sinoSize * uint64(options.TOF_bins), int32(options.ndist_side), pseudoD, pseudoR, type);
        clear mex


        if store_coordinates
            if Nt > 1
                coordinate = cell(Nt,1);
            end
            LL1(LL1==0) = [];
            LL2(LL2==0) = [];
            ring_number1 = uint16(idivide(LL1-1, options.det_per_ring)) + 1;
            ring_number2 = uint16(idivide(LL2-1, options.det_per_ring)) + 1;
            ring_pos1 = uint16(mod(LL1-1, options.det_per_ring)) + 1;
            ring_pos2 = uint16(mod(LL2-1, options.det_per_ring)) + 1;
            clear LL1 LL2
            [x, y] = detector_coordinates(options);
            z_length = single((options.rings + 3) * options.cr_pz);
            z = linspace(single(0), z_length, options.rings + 1 + 3)';
            z = z(2:end) - z(2) / 2;
            z([14;28;42]) = [];
            z = z - options.axial_fov/2 + (options.axial_fov - (options.cr_pz*(options.rings + 3)))/2;
            % if min(z(:)) == 0
            %     z = z + (options.axial_fov - (options.rings) * options.cr_pz)/2 + options.cr_pz/2;
            % end
            x = single(x);
            y = single(y);
            if Nt > 1
                for uu = 1 : Nt
                    ind1 = find(tpoints == uu, 1, 'first');
                    ind2 = find(tpoints == uu, 1, 'last');
                    tempPos1 = ring_pos1(ind1:ind2);
                    tempPos2 = ring_pos2(ind1:ind2);
                    tempNum1 = ring_number1(ind1:ind2);
                    tempNum2 = ring_number2(ind1:ind2);
                    coordinate{uu + 1} = [x(tempPos1)'; y(tempPos1)'; z(tempNum1)'; x(tempPos2)'; y(tempPos2)'; z(tempNum2)'];
                end
            else
                coordinate = [x(ring_pos1)'; y(ring_pos1)'; z(ring_number1)'; x(ring_pos2)'; y(ring_pos2)'; z(ring_number2)'];
            end
            clear ring_pos1 ring_pos2 ring_number1 ring_number2
            if options.randoms_correction
                if Nt > 1
                    Rcoordinate = cell(Nt,1);
                end
                DD1(DD1==0) = [];
                DD2(DD2==0) = [];
                ring_number1 = uint16(idivide(DD1-1, options.det_per_ring)) + 1;
                ring_number2 = uint16(idivide(DD2-1, options.det_per_ring)) + 1;
                ring_pos1 = uint16(mod(DD1-1, options.det_per_ring)) + 1;
                ring_pos2 = uint16(mod(DD2-1, options.det_per_ring)) + 1;
                clear DD1 DD2
                if Nt > 1
                    tpoints(tpoints==0) = [];
                    for uu = 1 : Nt
                        ind1 = find(tpoints == uu, 1, 'first');
                        ind2 = find(tpoints == uu, 1, 'last');
                        tempPos1 = ring_pos1(ind1:ind2);
                        tempPos2 = ring_pos2(ind1:ind2);
                        tempNum1 = ring_number1(ind1:ind2);
                        tempNum2 = ring_number2(ind1:ind2);
                        Rcoordinate{uu + 1} = [x(tempPos1)'; y(tempPos1)'; z(tempNum1)'; x(tempPos2)'; y(tempPos2)'; z(tempNum2)'];
                    end
                else
                    Rcoordinate = [x(ring_pos1)'; y(ring_pos1)'; z(ring_number1)'; x(ring_pos2)'; y(ring_pos2)'; z(ring_number2)'];
                end
                clear ring_pos1 ring_pos2 ring_number1 ring_number2
            end
        end
        % if options.store_raw_data
        %     if partitions == 1
        %         LL1(LL1 == 0) = [];
        %         LL2(LL2 == 0) = [];
        %         LL3 = LL2;
        %         LL2(LL2 > LL1) = LL1(LL2 > LL1);
        %         LL1(LL3 > LL1) = LL3(LL3 > LL1);
        %         clear LL3
        %         prompts = accumarray([LL1 LL2],1,[options.detectors options.detectors],[],[],true);
        %         if options.randoms_correction
        %             DD1(DD1 == 0) = [];
        %             DD2(DD2 == 0) = [];
        %             DD3 = DD2;
        %             DD2(DD2 > DD1) = DD1(DD2 > DD1);
        %             DD1(DD3 > DD1) = DD3(DD3 > DD1);
        %             clear DD3
        %             delays = accumarray([DD1 DD2],1,[options.detectors options.detectors],[],[],true);
        %         end
        %         clear LL1 LL2 DD1 DD2
        %     else
        %         prompts = cell(partitions,1);
        %         if options.randoms_correction
        %             delays = cell(partitions,1);
        %         end
        %         ll = 1;
        %         for kk = 1 : partitions
        %             apu1 = LL1(ll:tpoints(kk));
        %             apu1(apu1 == 0) = [];
        %             apu2 = LL2(ll:tpoints(kk));
        %             apu2(apu2 == 0) = [];
        %             apu3 = apu2;
        %             apu2(apu2 > apu1) = apu1(apu2 > apu1);
        %             apu1(apu3 > apu1) = apu3(apu3 > apu1);
        %             clear apu3
        %             prompts{kk} = accumarray([apu1 apu2],1,[options.detectors options.detectors],[],[],true);
        %             if options.randoms_correction
        %                 apu1 = DD1(ll:tpoints(kk));
        %                 apu1(apu1 == 0) = [];
        %                 apu2 = DD2(ll:tpoints(kk));
        %                 apu2(apu2 == 0) = [];
        %                 apu3 = apu2;
        %                 apu2(apu2 > apu1) = apu1(apu2 > apu1);
        %                 apu1(apu3 > apu1) = apu3(apu3 > apu1);
        %                 clear apu3
        %                 delays{kk} = accumarray([apu1 apu2],1,[options.detectors options.detectors],[],[],true);
        %             end
        %             ll = 1 + tpoints(kk);
        %         end
        %         clear LL1 LL2 DD1 DD2 apu1 apu2
        %     end
        % end
        %% Save the measurement data

        % coincidences = cell(partitions);
        % if options.randoms_correction
        %     delayed_coincidences = cell(partitions,1);
        % end
        %
        % tot_counts = 0;
        % if options.randoms_correction
        %     tot_delayed = 0;
        % end

        if ~options.use_raw_data
            if size(raw_SinM,1) == numel(raw_SinM)
                raw_SinM = reshape(raw_SinM, options.Ndist, options.Nang, totSinos, options.TOF_bins, Nt);
            end
            if options.randoms_correction && size(SinD,1) == numel(SinD)
                SinD = reshape(SinD, options.Ndist, options.Nang, totSinos, Nt);
            end
            if Nt > 1
                raw_SinM = squeeze(num2cell(raw_SinM, (1 : ndims(raw_SinM) - 1)));
                if options.randoms_correction
                    SinD = squeeze(num2cell(SinD, [1 2 3]));
                end
            end
            if ~options.randoms_correction
                SinD = [];
            end
            form_sinograms(options, true, raw_SinM, [], [], [], SinD);
        end
        % if options.store_raw_data
        %     for llo=1:partitions
        %
        %         if iscell(prompts)
        %             prompts{llo} = (prompts{llo}(tril(true(size(prompts{llo})), 0)));
        %             if options.verbose
        %                 counts = full(sum(sum(prompts{llo})));
        %             end
        %             coincidences{llo, 1} = prompts{llo};
        %         else
        %             prompts = (prompts(tril(true(size(prompts)), 0)));
        %             if options.verbose
        %                 counts = full(sum(sum(prompts)));
        %             end
        %             coincidences{llo, 1} = prompts;
        %         end
        %
        %         if options.verbose
        %             disp(['Total number of prompts at time point ' num2str(llo) ' is ' num2str(counts) '.'])
        %             tot_counts = tot_counts + counts;
        %         end
        %
        %         if options.randoms_correction
        %             if ~iscell(prompts)
        %                 delays = (delays(tril(true(size(delays)), 0)));
        %                 delayed_coincidences{partitions, 1} = delays;
        %                 if options.verbose
        %                     disp(['Total number of delayed coincidences at time point ' num2str(llo) ' is ' num2str(full(sum(sum(delays)))) '.'])
        %                     tot_delayed = tot_delayed + full(sum(sum(delays)));
        %                 end
        %             else
        %                 delays{llo} = (delays{llo}(tril(true(size(delays{llo})), 0)));
        %                 delayed_coincidences{partitions, 1} = delays{llo};
        %                 if options.verbose
        %                     disp(['Total number of delayed coincidences at time point ' num2str(llo) ' is ' num2str(full(sum(sum(delays{llo})))) '.'])
        %                     tot_delayed = tot_delayed + full(sum(sum(delays{llo})));
        %                 end
        %             end
        %         end
        %     end
        %     if partitions == 1
        %         if exist('OCTAVE_VERSION', 'builtin') == 0
        %             save(save_string_raw, 'coincidences', '-v7.3')
        %         else
        %             save(save_string_raw, 'coincidences', '-v7')
        %         end
        %         for variableIndex = 1:length(variableList)
        %             if exist(variableList{variableIndex},'var')
        %                 if exist('OCTAVE_VERSION', 'builtin') == 0
        %                     save(save_string_raw,variableList{variableIndex},'-append')
        %                 else
        %                     save(save_string_raw,variableList{variableIndex},'-append', '-v7')
        %                 end
        %             end
        %         end
        %     else
        %         if options.verbose
        %             disp(['Total number of prompts at all time points is ' num2str(tot_counts) '.'])
        %             if options.randoms_correction
        %                 disp(['Total number of delayed coincidences at all time points is ' num2str(tot_delayed) '.'])
        %             end
        %         end
        %         if exist('OCTAVE_VERSION', 'builtin') == 0
        %             save(save_string_raw, 'coincidences', '-v7.3')
        %         else
        %             save(save_string_raw, 'coincidences', '-v7')
        %         end
        %         for variableIndex = 1:length(variableList)
        %             if exist(variableList{variableIndex},'var')
        %                 if exist('OCTAVE_VERSION', 'builtin') == 0
        %                     save(save_string_raw,variableList{variableIndex},'-append')
        %                 else
        %                     save(save_string_raw,variableList{variableIndex},'-append', '-v7')
        %                 end
        %             end
        %         end
        %     end
        %     %     end
        %     clear LL1 DD1
        %     if nargout >= 1
        %         varargout{1} = coincidences;
        %     end
        %     if nargout >= 2
        %         varargout{2} = delayed_coincidences;
        %     end
        % else
        if nargout >= 1
            varargout{1} = raw_SinM;
        end
        if nargout >= 2 && options.randoms_correction
            varargout{2} = SinD;
        elseif nargout >= 2
            varargout{2} = [];
        end
        if nargout >= 3
            varargout{3} = coordinate;
        end
        if nargout >= 4 && options.randoms_correction
            varargout{4} = Rcoordinate;
        elseif nargout >= 4
            varargout{4} = [];
        end
        % end
        if options.verbose
            disp('Measurements loaded and saved')
            toc
        end
    elseif options.use_machine == 3
        % 32-bit list-mode data

        if options.randoms_correction
            variableList = {'SinM','SinDelayed','appliedCorrections'};
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
        if Nt == 1
            if options.TOF_bins > 1
                save_string = [machine_name '_' name '_TOFsinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.NSinos) ...
                    '_span' num2str(options.span) '_' num2str(options.TOF_bins) 'bins_' num2str(options.TOF_width*1e12) 'psBinSize_' num2str(options.TOF_noise_FWHM*1e12) ...
                    'psFWHM_listmode_sinogram.mat'];
            else
                save_string = [machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang)...
                    'x' num2str(options.TotSinos) '_span' num2str(options.span) '_listmode_sinogram.mat'];
            end
        else
            if options.TOF_bins > 1
                save_string = [machine_name '_' name '_TOFsinograms_combined_' num2str(Nt) 'timepoints_for_total_of_' num2str(tot_time) 's_' ...
                    num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.NSinos) '_span' num2str(options.span) '_' num2str(options.TOF_bins) 'bins_' ...
                    num2str(options.TOF_width*1e12) 'psBinSize_' num2str(options.TOF_noise_FWHM*1e12) 'psFWHM_listmode_sinogram.mat'];
            else
                save_string = [machine_name '_' name '_sinograms_combined_' num2str(Nt) 'timepoints_for_total_of_' num2str(tot_time) 's_' ...
                    num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.NSinos) '_span' num2str(options.span) '_listmode_sinogram.mat'];
            end
        end

        [raw_SinM, SinDelayed] = list2sinogram(nimi, partitions, (alku), (loppu), logical(options.randoms_correction), v_size, options.store_raw_data, ...
            uint32(options.det_w_pseudo), sinoSize, uint32(options.Ndist), uint32(options.Nang), uint32(options.ring_difference), uint32(options.span), uint32(options.segment_table), ...
            Nt, sinoSize * uint64(options.TOF_bins), int32(options.ndist_side), pseudoD, pseudoR, type, usePTD);
        clear mex

        if Nt == 1
            raw_SinM = permute(reshape(raw_SinM, options.Nang, options.Ndist, options.TotSinos, options.TOF_bins), [2 1 3 4]);
            if options.NSinos < options.TotSinos
                raw_SinM = raw_SinM(:,:, 1: options.NSinos, :);
            end
            if options.randoms_correction
                SinDelayed = permute(reshape(SinDelayed, options.Nang, options.Ndist, options.TotSinos), [2 1 3]);
            end
            SinM = single(raw_SinM);
            form_sinograms(options, true, raw_SinM, [], [], [], SinDelayed);
            % if (~options.corrections_during_reconstruction && (options.randoms_correction || options.normalization_correction)) || ...
            %         options.fill_sinogram_gaps && (sum(options.pseudot) > 0 || options.det_per_ring < options.det_w_pseudo)
            %     if ~options.corrections_during_reconstruction && options.randoms_correction
            %         SinM = SinM - single(SinDelayed) / options.TOF_bins;
            %         appliedCorrections.randoms = true;
            %     end
            %     if options.normalization_correction && ~options.corrections_during_reconstruction && options.use_user_normalization
            %         [n_file, fpath] = uigetfile({'*.n;*.ptd'},'Select normalization datafile');
            %         if isequal(n_file, 0)
            %             error('No file was selected')
            %         end
            %         f_path = [fpath n_file];
            %         normalization = loadNormalizationBiograph(f_path, options);
            %         %                     normalization = circshift(normalization,1);
            %         SinM = bsxfun(@times, single(SinM), normalization);
            %         appliedCorrections.normalization = true;
            %     elseif options.normalization_correction && ~options.corrections_during_reconstruction && ~options.use_user_normalization
            %         norm_file = [folder options.machine_name '_normalization_' num2str(options.Ndist) 'x' num2str(options.Nang) '_span' num2str(options.span) '.mat'];
            %         if exist(norm_file, 'file') == 2
            %             normalization = loadStructFromFile(norm_file,'normalization');
            %         else
            %             normalization = normalization_coefficients(options);
            %         end
            %         normalization = circshift(normalization,1);
            %         SinM = bsxfun(@times, single(SinM), normalization);
            %         appliedCorrections.normalization = true;
            %     end
            %     if options.fill_sinogram_gaps && (sum(options.pseudot) > 0 || options.det_per_ring < options.det_w_pseudo)
            %         [~, ~, xp, yp] = detector_coordinates(options);
            %         gaps = [];
            %         llo = 1;
            %         [SinM, gaps] = gapFilling(options, SinM, xp, yp, llo, gaps);
            %         appliedCorrections.gapFilling = true;
            %     end
            % end
        else
            if options.TOF_bins > 1
                apu = reshape(raw_SinM, options.Nang, options.Ndist, options.TotSinos, options.TOF_bins, Nt);
            else
                apu = reshape(raw_SinM, options.Nang, options.Ndist, options.TotSinos, Nt);
            end
            if options.NSinos < options.TotSinos
                apu = apu(:,:, 1: options.NSinos, :, :);
            end
            raw_SinM = cell(Nt,1);
            if ~options.corrections_during_reconstruction && (options.randoms_correction || options.normalization_correction || ...
                    options.fill_sinogram_gaps && (sum(options.pseudot) > 0 || options.det_per_ring < options.det_w_pseudo))
                SinM = cell(Nt,1);
            end
            if options.randoms_correction
                apu_d = reshape(SinDelayed, options.Nang, options.Ndist, options.TotSinos, Nt);
                SinDelayed = cell(Nt,1);
            end
            for kk = 1 : Nt
                if options.TOF_bins > 1
                    raw_SinM{kk} = permute(apu(:,:,:,:,kk), [2 1 3 4 5]);
                else
                    raw_SinM{kk} = permute(apu(:,:,:,kk), [2 1 3 4]);
                end
                if options.randoms_correction
                    SinDelayed{kk} = permute(apu_d, [2 1 3 4]);
                end
                % if (~options.corrections_during_reconstruction && (options.randoms_correction || options.normalization_correction)) || ...
                %         options.fill_sinogram_gaps && (sum(options.pseudot) > 0 || options.det_per_ring < options.det_w_pseudo)
                %     Sino = single(raw_SinM);
                %     if ~options.corrections_during_reconstruction && options.randoms_correction
                %         Sino = Sino - single(SinDelayed) / options.TOF_bins;
                %         appliedCorrections.randoms = true;
                %     end
                %     if options.normalization_correction && ~options.corrections_during_reconstruction
                %         [n_file, fpath] = uigetfile({'*.n;*.ptd'},'Select normalization datafile');
                %         if isequal(n_file, 0)
                %             error('No file was selected')
                %         end
                %         f_path = [fpath n_file];
                %         normalization = loadNormalizationBiograph(f_path, options);
                %         Sino = bsxfun(@times, single(Sino), normalization);
                %         appliedCorrections.normalization = true;
                %     end
                %     if options.fill_sinogram_gaps && (sum(options.pseudot) > 0 || options.det_per_ring < options.det_w_pseudo)
                %         [~, ~, xp, yp] = detector_coordinates(options);
                %         gaps = [];
                %         [Sino, gaps] = gapFilling(options, Sino, xp, yp, 1, gaps);
                %         appliedCorrections.gapFilling = true;
                %     end
                %     SinM{kk} = Sino;
                % end
            end
            form_sinograms(options, true, raw_SinM, [], [], [], SinDelayed);
        end

        if Nt == 1
            if options.verbose
                disp(['Total number of prompts is ' num2str(sum(raw_SinM(:))) '.'])
                if options.randoms_correction
                    disp(['Total number of delayed coincidences is ' num2str(sum(SinDelayed(:))) '.'])
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
        % if exist('OCTAVE_VERSION','builtin') == 0
        %     save(save_string,'raw_SinM','-v7.3')
        %     for variableIndex = 1:length(variableList)
        %         if exist(variableList{variableIndex},'var')
        %             save(save_string,variableList{variableIndex},'-append')
        %         end
        %     end
        % else
        %     save(save_string,'raw_SinM','-v7')
        %     for variableIndex = 1:length(variableList)
        %         if exist(variableList{variableIndex},'var')
        %             save(save_string,variableList{variableIndex},'-append','-v7')
        %         end
        %     end
        % end
        if nargout >= 1
            varargout{1} = SinM;
        end
        if nargout >= 2
            varargout{2} = SinDelayed;
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
    if Nt == 1
        save_string = [options.machine_name '_' options.name '_sinogram_original_static_' num2str(options.Ndist) 'x' num2str(options.Nang)...
            'x' num2str(options.TotSinos) '_span' num2str(options.span) '_machine_sinogram.mat'];
    else
        save_string = [options.machine_name '_' options.name '_sinograms_combined_' num2str(Nt) 'timepoints_for_total_of_' num2str(tot_time) 's_' ...
            num2str(Ndist) 'x' num2str(Nang) 'x' num2str(NSlices) '_span' num2str(span) '_machine_sinogram.mat'];
    end
    if exist(options.fpath,'file') ~= 2
        [options.file, options.fpath] = uigetfile({'*.ptd;*.s'},'Select mCT sinogram datafile');
        if isequal(options.file, 0)
            error('No file was selected!')
        end
        nimi = [options.fpath options.file];
    else
        nimi = [options.fpath];
    end
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
        offset_loppu = length(Sino) - offset_alku - options.Nang*options.Ndist*options.TotSinos*Nt * TOF_timebins * 2;
        fseek(fid, offset_alku, 0);
        Sino = fread(fid, inf, 'int16=>int16');
        Sino(end - offset_loppu + 1 : end) = [];
        Sino = reshape(Sino, options.Ndist, options.Nang, options.TotSinos, TOF_timebins, Nt);
        if TOF_timebins > 1
            SinDelayed = squeeze(Sino(:,:,:,end,:));
            Sino = Sino(:,:,:,1:end-1,:);
        end
        if options.TOF_bins > 1
            Sino = squeeze(Sino);
        else
            Sino = squeeze(sum(Sino,4));
        end
    else
        Sino = fread(fid, inf, 'int16=>int16');
        if any(Sino < 0)
            fid = fopen(nimi);
            Sino = fread(fid, inf, 'single=>single');
        end
        TOF_timebins = length(Sino) / (options.Ndist*options.Nang*options.TotSinos*Nt);
        Sino = reshape(Sino, options.Ndist, options.Nang, options.TotSinos, TOF_timebins, Nt);
        if TOF_timebins > 1
            SinDelayed = squeeze(Sino(:,:,:,end,:));
            Sino = Sino(:,:,:,1:end-1,:);
        end
        if options.TOF_bins > 1
            Sino = squeeze(Sino);
        else
            Sino = squeeze(sum(Sino,4));
        end
        Sino(Sino < 0) = 0;
    end
    fclose(fid);
    if Nt > 1
        raw_SinM = cell(Nt,1);
        for kk = 1 : Nt
            if options.TOF_bins > 1
                raw_SinM{kk} = Sino(:,:,:,:,kk);
            else
                raw_SinM{kk} = Sino(:,:,:,kk);
            end
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
        varargout{2} = SinDelayed;
    end

end
disp('Data load complete')
