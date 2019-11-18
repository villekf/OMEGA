function [varargout] = load_data(options)
%% Load PET data
% Loads either ASCII, ROOT, or LMF GATE format data or the Inveon list-mode
% (.lst) data into MATLAB/Octave.
%
% Outputs data formatted for pure, raw list-mode, data. This data can be
% later transformed into sinogram format or used as-is. Loads also the
% delayed coincidences, true randoms, true scatter and true coincidences
% separately if they are selected.
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
%   true_coincidences = True coincidences (GATE only)
%   scattered_coincidences = Scattered coincidences (GATE only)
%   random_coincidences = True random coincidences (GATE only)
%
% See also form_sinograms, initial_michelogram

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

if nargout > 5
    error('Too many output arguments, there can be at most five')
end

% folder = fileparts(which('load_data.m'));
% folder = strrep(folder, 'source','mat-files/');
% folder = strrep(folder, '\','/');

Nx = uint32(options.Nx);
Ny = uint32(options.Ny);
Nz = uint32(options.Nz);

partitions = options.partitions;

alku = options.start;
loppu = options.end;
if isinf(loppu)
    loppu = 1e9;
end
vali = (loppu - alku)/options.partitions;
tot_time = options.tot_time;
    
machine_name = options.machine_name;
name = options.name;
detectors = options.detectors;

disp('Beginning data load')

%% Load Inveon list-mode data
if options.use_machine == 1
    
    % Select the list-mode file
    [options.file, options.fpath] = uigetfile('*.lst','Select Inveon list-mode datafile');
    if isequal(options.file, 0)
        error('No file was selected')
    end
    
    if options.verbose
        tic
    end
    
    if partitions > 1
        save_string_raw = [machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_listmode.mat'];
    else
        save_string_raw = [machine_name '_measurements_' name '_static_raw_listmode.mat'];
    end
    
    variableList = {'delayed_coincidences'};
    
    
    nimi = [options.fpath options.file];
    f = dir(nimi);
    pituus = uint64(f.bytes/6);
    
    if options.verbose
        disp('Data load started')
    end
    
    % Load the data from the binary file
    [LL1, LL2, tpoints, DD1, DD2, ~] = inveon_list2matlab(nimi,(vali),(alku),(loppu),pituus, uint32(detectors), options.randoms_correction);
%     clear mex
    
    
    % If no time steps, the output matirx LL1 is a options.detectors x
    % options.detectors matrix in unsigned short (u16) format
    % Contains counts for each detector pair
    if partitions == 1
        prompts = LL1;
        if options.randoms_correction
            delays = DD1;
        end
        clear LL1 LL2 DD1 DD2
        % Else the outputs are vectors LL1 and LL2 containing the detector
        % pairs for the counts
        % Form a similar (sparse) detectors x detectors matrix
    else
        prompts = cell(partitions,1);
        if options.randoms_correction
            delays = cell(partitions,1);
        end
%         ll = 1;
        % Obtain the counts that were within a specified time-window
        for kk = 1 : partitions
            apu1 = LL1(tpoints(kk) + 1:tpoints(kk+1));
            apu1(apu1 == 0) = [];
            apu2 = LL2(tpoints(kk) + 1:tpoints(kk+1));
            apu2(apu2 == 0) = [];
            prompts{kk} = accumarray([apu1 apu2],1,[options.detectors options.detectors],[],[],true);
            if options.randoms_correction
                apu1 = DD1(tpoints(kk) + 1:tpoints(kk+1));
                apu1(apu1 == 0) = [];
                apu2 = DD2(tpoints(kk) + 1:tpoints(kk+1));
                apu2(apu2 == 0) = [];
                delays{kk} = accumarray([apu1 apu2],1,[options.detectors options.detectors],[],[],true);
            end
%             ll = ll + tpoints(kk+1);
        end
        clear LL1 LL2 DD1 DD2 apu1 apu2
    end
    
    %% GATE data
elseif options.use_machine == 0
    if options.use_ASCII && options.use_LMF && options.use_root
        warning('ASCII, LMF and ROOT selected, using only ASCII')
        options.use_LMF = false;
        options.use_root = false;
    elseif options.use_ASCII && options.use_LMF
        warning('Both ASCII and LMF selected, using only ASCII')
        options.use_LMF = false;
    elseif options.use_ASCII && options.use_root
        warning('Both ASCII and ROOT selected, using only ASCII')
        options.use_root = false;
    elseif options.use_LMF && options.use_root
        warning('Both LMF and ROOT selected, using only LMF')
        options.use_root = false;
    elseif options.use_ASCII == false && options.use_LMF == false && options.use_root == false
        error('Error: Neither ASCII, LMF nor ROOT data selected')
    end
    if options.verbose
        tic
    end
    
    fpath = options.fpath;
    blocks_per_ring = options.blocks_per_ring;
    cryst_per_block = options.cryst_per_block;
    machine_name = options.machine_name;
    rings = options.rings;
    det_per_ring = options.det_per_ring;
    pseudot = options.pseudot;
    source = options.source;
    
    temp = pseudot;
    if ~isempty(temp) && sum(temp) > 0
        for kk = 1 : temp
            pseudot(kk) = (options.cryst_per_block + 1) * kk;
        end
    end
    variableList = {'true_coincidences','scattered_coincidences','random_coincidences','delayed_coincidences'};
    
    % Total number of detectors
    detectors = det_per_ring*rings;
    
    time_intervals = linspace(alku, loppu, options.partitions + 1);
    
    
    
    %%
    if options.use_LMF
        blocks_per_ring = uint32(blocks_per_ring);
        cryst_per_block = uint32(cryst_per_block);
        R_bits = int32(options.R_bits);
        M_bits = int32(options.M_bits);
        S_bits = int32(options.S_bits);
        C_bits = int32(options.C_bits);
        L_bits = int32(options.L_bits);
        coincidence_window = options.coincidence_window;
        header_bytes = int32(options.header_bytes);
        data_bytes = int32(options.data_bytes);
        det_per_ring = uint32(det_per_ring);
        
        coincidence_window = coincidence_window / options.clock_time_step;
        
        if partitions > 1
            save_string_raw = [machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_LMF.mat'];
        else
            save_string_raw = [machine_name '_measurements_' name '_static_raw_LMF.mat'];
        end
        
        % Take coincidence data
        fnames = dir([fpath '*.ccs']);
        if size(fnames,1) == 0
            fpath = pwd;
            fpath = [fpath '/'];
            fnames = dir([fpath '*.ccs']);
            if size(fnames,1) == 0
                error('No LMF (.ccs) files were found. Check your filepath (options.fpath) or current folder.')
            end
        end
        
%         bina = [];
        bina = char(zeros(1,R_bits));
        for kk = 1 : R_bits
%             bina = [bina, '1'];
            bina(kk) = '1';
        end
        R_length = int32(bin2dec(bina));
        bina = char(zeros(1,M_bits));
        for kk = 1 : M_bits
%             bina = [bina, '1'];
            bina(kk) = '1';
        end
        M_length = int32(bin2dec(bina));
        bina = char(zeros(1,C_bits));
%         bina = [];
        for kk = 1 : C_bits
%             bina = [bina, '1'];
            bina(kk) = '1';
        end
        C_length = int32(bin2dec(bina));
        
        
        if options.partitions > 1
            prompts = cell(partitions,1);
            if options.obtain_trues
                trues = cell(partitions,1);
            end
            if options.store_scatter
                scatter = cell(partitions,1);
            end
            if options.store_randoms
                randoms = cell(partitions,1);
            end
        else
            prompts = zeros(detectors, detectors, 'uint16');
            if options.obtain_trues
                trues = zeros(detectors, detectors, 'uint16');
            end
            if options.store_scatter
                scatter = zeros(detectors, detectors, 'uint16');
            end
            if options.store_randoms
                randoms = zeros(detectors, detectors, 'uint16');
            end
        end
        if source
            if options.obtain_trues
                C = cell(7,partitions);
            else
                C = cell(6,partitions);
            end
            C(:) = {zeros(options.Nx, options.Ny, options.Nz, 'uint16')};
        end
        
        if options.verbose
            disp('Data load started')
        end
        
        % Go through all the files
        for lk=1:length(fnames)
            
            pituus = int64((fnames(lk).bytes - header_bytes) / data_bytes);
            
            nimi = [fpath fnames(lk).name];
            
            [L1, L2, tpoints, A, int_loc, Ltrues, Lscatter, Lrandoms, trues_index] = gate_lmf_matlab(nimi,vali,alku,loppu,pituus, uint32(detectors), blocks_per_ring, cryst_per_block, ...
                det_per_ring, uint32(options.linear_multip), header_bytes, data_bytes, R_bits, M_bits, S_bits, C_bits, L_bits, R_length, M_length, C_length, source, ...
                coincidence_window, options.clock_time_step, time_intervals, options.obtain_trues, options.store_scatter, options.store_randoms);
            
            trues_index = logical(trues_index);
            
%             int_loc = int_loc + 1;
            
            % An ideal image can be formed with this, see the file
            % visualize_pet.m
            ll = 1;
            if source && int_loc(1) > 0
                
                for jj = int_loc(1) : min(int_loc(2),partitions)
                    if partitions > 1
                        S = A(tpoints(ll) + 1:tpoints(ll+1),:);
                        if obtain_trues
                            trues_index = logical(Ltrues(tpoints(ll) + 1:tpoints(ll+1),:));
                            inde = any(S,2);
                            trues_index = trues_index(inde);
                        end
                        S = single(S(inde,:));
                        ll = ll + 1;
                    else
                        inde = any(A,2);
                        trues_index = logical(trues_index(inde));
                        A = A(inde,:);
                        S = single(A);
                    end
                    pixel_width_x = options.FOVa_x/options.Nx;
                    pixel_width_y = options.FOVa_y/options.Ny;
                    z_width = options.axial_fov/options.Nz;
                    x = single((-pixel_width_x*(options.Nx-1)/2:pixel_width_x:pixel_width_x*(options.Nx-1)/2)');
                    y = single((-pixel_width_y*(options.Ny-1)/2:pixel_width_y:pixel_width_y*(options.Ny-1)/2)');
                    z = single((-z_width*(options.Nz-1)/2:z_width:z_width*(options.Nz-1)/2)');
                    max_x = max(x) + pixel_width_x/2;
                    max_y = max(y) + pixel_width_y/2;
                    max_z = max(z) + z_width/2;
                    
                    CC = ones(length(S),1,'int16');
                    t1 = all(S(:,1:3)==S(:,4:6),2);
                    t2 = all([abs(S(:,1))<=max_x abs(S(:,2))<=max_y abs(S(:,3))<=max_z],2);
                    t3 = all([abs(S(:,4))<=max_x abs(S(:,5))<=max_y abs(S(:,6))<=max_z],2);
                    t = (t1+t2+t3==3);
                    CC = CC(t);
                    [~,i] = min(abs(bsxfun(@minus,x.',S(t,1))),[],2);
                    [~,j] = min(abs(bsxfun(@minus,y.',S(t,2))),[],2);
                    [~,k] = min(abs(bsxfun(@minus,z.',S(t,3))),[],2);
                    FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                    %                 FOV = circshift(FOV, round(size(FOV,1)/6));
                    C{1,jj} = C{1,jj} + FOV;
                    
                    CC = ones(length(S),1,'int16');
                    t1 = true(size(S,1),1);
                    t2 = all([abs(S(:,1))<=max_x abs(S(:,2))<=max_y abs(S(:,3))<=max_z],2);
                    t = (t1+t2==2);
                    CC = CC(t);
                    [~,i] = min(abs(bsxfun(@minus,x.',S(t,1))),[],2);
                    [~,j] = min(abs(bsxfun(@minus,y.',S(t,2))),[],2);
                    [~,k] = min(abs(bsxfun(@minus,z.',S(t,3))),[],2);
                    FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                    %                 FOV = circshift(FOV, round(size(FOV,1)/6));
                    C{2,jj} = C{2,jj} + FOV;
                    
                    CC = ones(length(S),1,'int16');
                    t1 = true(size(S,4),1);
                    t2 = all([abs(S(:,4))<=max_x abs(S(:,5))<=max_y abs(S(:,6))<=max_z],2);
                    t = (t1+t2==2);
                    CC = CC(t);
                    [~,i] = min(abs(bsxfun(@minus,x.',S(t,4))),[],2);
                    [~,j] = min(abs(bsxfun(@minus,y.',S(t,5))),[],2);
                    [~,k] = min(abs(bsxfun(@minus,z.',S(t,6))),[],2);
                    FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                    %                 FOV = circshift(FOV, round(size(FOV,1)/6));
                    C{3,jj} = C{3,jj} + FOV;
                    
                    C{4,jj} = C{2,jj} + C{3,jj};
                    
                    CC = ones(length(S),1,'int16');
                    t1 = true(size(S,1),1);
                    t2 = all([abs(S(:,1))<=max_x abs(S(:,2))<=max_y abs(S(:,3))<=max_z],2);
                    t3 = all([abs(S(:,4))<=max_x abs(S(:,5))<=max_y abs(S(:,6))<=max_z],2);
                    t = (t1+t2+t3==3);
                    CC = CC(t);
                    [~,i] = min(abs(bsxfun(@minus,(x).',S(t,1))),[],2);
                    [~,j] = min(abs(bsxfun(@minus,(y).',S(t,2))),[],2);
                    [~,k] = min(abs(bsxfun(@minus,(z).',S(t,3))),[],2);
                    [~,i2] = min(abs(bsxfun(@minus,(x).',S(t,4))),[],2);
                    [~,j2] = min(abs(bsxfun(@minus,(y).',S(t,5))),[],2);
                    [~,k2] = min(abs(bsxfun(@minus,(z).',S(t,6))),[],2);
                    i = [i;i2];
                    j = [j;j2];
                    k = [k;k2];
                    FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],[CC;CC],[Nx Ny Nz])),-2));
                    %                 FOV = circshift(FOV, round(size(FOV,1)/6));
                    FOV = FOV/2;
                    C{5,jj} = C{5,jj} + FOV;
                    
                    CC = ones(length(S),1,'int16');
                    t1 = true(size(S,1),1);
                    t2 = all([abs(S(:,1))<=max_x abs(S(:,2))<=max_y abs(S(:,3))<=max_z],2);
                    t3 = all([abs(S(:,4))<=max_x abs(S(:,5))<=max_y abs(S(:,6))<=max_z],2);
                    t = (t1+t2+t3==3);
                    CC = CC(t);
                    C1 = mean([S(t,1) S(t,4)],2);
                    C2 = mean([S(t,2) S(t,5)],2);
                    C3 = mean([S(t,3) S(t,6)],2);
                    [~,i] = min(abs(bsxfun(@minus,(x).',C1)),[],2);
                    [~,j] = min(abs(bsxfun(@minus,(y).',C2)),[],2);
                    [~,k] = min(abs(bsxfun(@minus,(z).',C3)),[],2);
                    FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                    %                 FOV = circshift(FOV, round(size(FOV,1)/6));
                    C{6,jj} = C{6,jj} + FOV;
                    
                    if options.obtain_trues
                        CC = ones(length(S),1,'int16');
                        CC = CC(trues_index);
                        SS = S(trues_index,:);
                        t= all([abs(SS(:,1))<=max_x abs(SS(:,2))<=max_y abs(SS(:,3))<=max_z],2);
                        CC = CC(t);
                        [~,i] = min(abs(bsxfun(@minus,x.',SS(t,1))),[],2);
                        [~,j] = min(abs(bsxfun(@minus,y.',SS(t,2))),[],2);
                        [~,k] = min(abs(bsxfun(@minus,z.',SS(t,3))),[],2);
                        if isempty(i)
                            FOV = zeros([Nx Ny Nz],'uint16');
                        else
                            FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                        end
                        %                     FOV = circshift(FOV, round(size(FOV,1)/6));
                        C{7,jj} = C{7,jj} + FOV;
                    end
                end
            end
            
            if partitions == 1
                prompts = prompts + L1;
                if options.obtain_trues
                    trues = trues + Ltrues;
                end
                if options.store_scatter
                    scatter = scatter + Lscatter;
                end
                if options.store_randoms
                    randoms = randoms + Lrandoms;
                end
            else
                % forms a sparse matrix containing all the events at different LORs
                % e.g. a value at [325 25100] contains the number of events detected at
                % detectors 325 (detector 1) and 25100 (detector 2)
                % a value at [25100 325] on the other hand tells the number of events when
                % 25100 is detector 1 and 325 detector 2
                if  int_loc(1) > 0
                    ll = 1;
                    for kk = int_loc(1) : min(int_loc(2),partitions)
                        %                 index = find(M(:,time_index)<tpoints(kk),1,'last');
                        apu1 = L1(tpoints(ll) + 1:tpoints(ll+1));
                        apu2 = L2(tpoints(ll) + 1:tpoints(ll+1));
                        inde = any(apu1,2);
                        apu1 = apu1(inde,:);
                        apu2 = apu2(inde,:);
                        if isempty(prompts{kk})
                            prompts{kk} = accumarray([apu1 apu2],1,[detectors detectors],[],[],true);
                        else
                            prompts{kk} = prompts{kk} + accumarray([apu1 apu2],1,[detectors detectors],[],[],true);
                        end
                        if options.obtain_trues
                            trues_index = logical(Ltrues(tpoints(ll) + 1:tpoints(ll+1)));
                            trues_index = trues_index(inde);
                            if isempty(trues{kk})
                                trues{kk} = accumarray([apu1(trues_index) apu2(trues_index)],1,[detectors detectors],[],[],true);
                            else
                                trues{kk} = trues{kk} + accumarray([apu1(trues_index) apu2(trues_index)],1,[detectors detectors],[],[],true);
                            end
                        end
                        if options.store_scatter
                            trues_index = logical(Lscatter(tpoints(ll) + 1:tpoints(ll+1)));
                            trues_index = trues_index(inde);
                            if isempty(scatter{kk})
                                scatter{kk} = accumarray([apu1(trues_index) apu2(trues_index)],1,[detectors detectors],[],[],true);
                            else
                                scatter{kk} = scatter{kk} + accumarray([apu1(trues_index) apu2(trues_index)],1,[detectors detectors],[],[],true);
                            end
                        end
                        if options.store_randoms
                            trues_index = logical(Lrandoms(tpoints(ll) + 1:tpoints(ll+1)));
                            trues_index = trues_index(inde);
                            if isempty(randoms{kk})
                                randoms{kk} = accumarray([apu1(trues_index) apu2(trues_index)],1,[detectors detectors],[],[],true);
                            else
                                randoms{kk} = randoms{kk} + accumarray([apu1(trues_index) apu2(trues_index)],1,[detectors detectors],[],[],true);
                            end
                        end
                        ll = ll + 1;
                    end
                end
            end
            if options.verbose
                disp(['File ' fnames(lk).name ' loaded'])
            end
            
        end
        clear Lrandoms apu2 Lscatter Ltrues apu1 trues_index M FOV CC i j k t t1 t2 t3 A S inde
        if source
            save([machine_name '_Ideal_image_coordinates_' name '_LMF.mat'],'C')
        end
        
    elseif options.use_ASCII
        %%
        ascii_ind = get_ascii_indices(options.coincidence_mask);
        rsector_ind1 = ascii_ind.rsector_ind1;
        rsector_ind2 = ascii_ind.rsector_ind2;
        crs_ind1 = ascii_ind.crs_ind1;
        crs_ind2 = ascii_ind.crs_ind2;
        time_index = ascii_ind.time_index;
        source_index1 = ascii_ind.source_index1;
        source_index2 = ascii_ind.source_index2;
        det_per_ring = options.det_per_ring;
        
        if partitions > 1
            save_string_raw = [machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_ASCII.mat'];
        else
            save_string_raw = [machine_name '_measurements_' name '_static_raw_ASCII.mat'];
        end
        
        % Take coincidence data
        fnames = dir([fpath '*Coincidences*.dat']);
        
        if size(fnames,1) == 0
            fpath = pwd;
            fpath = [fpath '/'];
            fnames = dir([fpath '*Coincidences*.dat']);
            if size(fnames,1) == 0
                error('No ASCII (.dat) coincidence files were found. Check your filepath (options.fpath) or current folder.')
            end
        end
        
        if options.randoms_correction
            delay_names = dir([fpath '*delay*.dat']);
            if size(delay_names,1) == 0
                fpath = pwd;
                fpath = [fpath '/'];
                delay_names = dir([fpath '*delay*.dat']);
                if size(delay_names,1) == 0
                    error('No ASCII (.dat) coincidence files were found. Check your filepath (options.fpath) or current folder.')
                end
            end
        end
        
        if partitions > 1
            prompts = cell(partitions,1);
            if options.obtain_trues
                trues = cell(partitions,1);
            end
            if options.store_scatter
                scatter = cell(partitions,1);
            end
            if options.store_randoms
                randoms = cell(partitions,1);
            end
            if options.randoms_correction
                delays = cell(partitions,1);
            end
        else
            prompts = zeros(detectors,detectors,'uint16');
            if options.obtain_trues
                trues = zeros(detectors,detectors,'uint16');
            end
            if options.store_scatter
                scatter = zeros(detectors,detectors,'uint16');
            end
            if options.store_randoms
                randoms = zeros(detectors,detectors,'uint16');
            end
            if options.randoms_correction
                delays = zeros(detectors,detectors,'uint16');
            end
        end
        if source
            if options.obtain_trues
                C = cell(7,partitions);
            else
                C = cell(6,partitions);
            end
            C(:) = {zeros(options.Nx, options.Ny, options.Nz, 'uint16')};
            if options.store_scatter
                SC = zeros(options.Nx, options.Ny, options.Nz, 'uint16');
            end
            if options.store_randoms
                RA = zeros(options.Nx, options.Ny, options.Nz, 'uint16');
            end
        end
        li = 1;
        
        if options.verbose
            disp('Data load started')
        end
        
        % Go through all the files
        for lk=1:length(fnames)
            
            if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','9.6') || exist('OCTAVE_VERSION','builtin') == 5
                M = importdata([fpath fnames(lk).name]);
                if any(any(isnan(M))) > 0
                    header = size(M,1);
                    warning(['Line ' num2str(header) ' corrupted, skipping.'])
                    B = importdata([fpath fnames(lk).name],' ', header);
                    while any(any(isnan(B.data))) > 0
                        header = header + size(B.data,1);
                        M = [M; B.data(1:end-1,:)];
                        B = importdata([fpath fnames(lk).name],' ', header);
                        warning(['Line ' num2str(header) ' corrupted, skipping.'])
                    end
                    M = [M; B.data];
                    clear B
                end
            else
                M = readmatrix([fpath fnames(lk).name]);
                if sum(isnan(M(:,end))) > 0 && sum(isnan(M(:,end))) > round(size(M,1)/2)
                    [rows, ~] = find(~isnan(M(:,end)));
                    for kg = 1 : length(rows)
                        warning(['Line ' num2str(rows(kg)) ' corrupted, skipping.'])
                        M(rows(kg),:) = [];
                    end
                else
                    if any(any(isnan(M))) > 0
                        [rows, ~] = find(isnan(M));
                        rows = unique(rows);
                        for kg = 1 : length(rows)
                            warning(['Line ' num2str(rows(kg)) ' corrupted, skipping.'])
                            M(rows(kg),:) = [];
                        end
                    end
                end
            end
            if ascii_ind.module_ind1 > 0 && sum(M(:,ascii_ind.module_ind1)) == 0
                ascii_ind.module_ind1 = 0;
                ascii_ind.module_ind2 = 0;
            end
            
            if partitions > 1
                if isempty(find(M(:,time_index) > time_intervals(1),1,'last'))
                    int_loc(1) = 0;
                elseif isempty(find(M(1,time_index) > time_intervals,1,'last'))
                    int_loc(1) = 1;
                else
                    int_loc(1) = find(M(1,time_index) > time_intervals,1,'last');
                end
                if isempty(find(M(:,time_index) < time_intervals(end),1,'first'))
                    int_loc(2) = 0;
                elseif isempty(find(M(end,time_index) < time_intervals,1,'first'))
                    int_loc(2) = length(time_intervals) - 1;
                else
                    int_loc(2) = find(M(end,time_index) < time_intervals,1,'first') - 1;
                end
            else
                int_loc(1) = 1;
                int_loc(2) = 1;
            end
            
            if options.obtain_trues
                trues_index = true(size(M,1),1);
                ind = (M(:,ascii_ind.event_index1) == M(:,ascii_ind.event_index2)) & (M(:,ascii_ind.scatter_index_cp1) == 0 & M(:,ascii_ind.scatter_index_cp2) == 0);
                trues_index(~ind) = false;
                
                ind = (M(:,ascii_ind.event_index1) == M(:,ascii_ind.event_index2)) & (M(:,ascii_ind.scatter_index_cd1) == 0 & M(:,ascii_ind.scatter_index_cd1) == 0);
                trues_index(~ind) = false;

                ind = (M(:,ascii_ind.event_index1) == M(:,ascii_ind.event_index2)) & (M(:,ascii_ind.scatter_index_rp1) == 0 & M(:,ascii_ind.scatter_index_rp1) == 0);
                trues_index(~ind) = false;

                ind = (M(:,ascii_ind.event_index1) == M(:,ascii_ind.event_index2)) & (M(:,ascii_ind.scatter_index_rd1) == 0 & M(:,ascii_ind.scatter_index_rd1) == 0);
                trues_index(~ind) = false;
                %                 end
            end
            if options.store_scatter
                scatter_index = false(size(M,1),1);
                if options.scatter_components(1)
                    ind = (M(:,ascii_ind.scatter_index_cp1) > 0 | M(:,ascii_ind.scatter_index_cp2) > 0) & (M(:,ascii_ind.event_index1) == M(:,ascii_ind.event_index2));
                    scatter_index(ind) = true;
                end
                if options.scatter_components(2)
                    ind = (M(:,ascii_ind.scatter_index_cd1) > 0 | M(:,ascii_ind.scatter_index_cd2) > 0) & (M(:,ascii_ind.event_index1) == M(:,ascii_ind.event_index2));
                    scatter_index(ind) = true;
                end
                if options.scatter_components(3)
                    ind = (M(:,ascii_ind.scatter_index_rp1) > 0 | M(:,ascii_ind.scatter_index_rp2) > 0) & (M(:,ascii_ind.event_index1) == M(:,ascii_ind.event_index2));
                    scatter_index(ind) = true;
                end
                if options.scatter_components(4)
                    ind = (M(:,ascii_ind.scatter_index_rd1) > 0 | M(:,ascii_ind.scatter_index_rd2) > 0) & (M(:,ascii_ind.event_index1) == M(:,ascii_ind.event_index2));
                    scatter_index(ind) = true;
                end
            end
            if options.store_randoms
                randoms_index = M(:,ascii_ind.event_index1) ~= M(:,ascii_ind.event_index2);
            end
            clear ind
            
            % An ideal image can be formed with this, see the file
            % Compare_Phantom_vs_Reconstructions.m
            % Takes the columns containing the source locations for both singles
            if source_index1 ~= 0 && ~isempty(source_index1) && source_index2 ~= 0 && ~isempty(source_index2) && int_loc(1) > 0 && source
                
                A = single(M(:,[source_index1:source_index1+2 source_index2:source_index2+2]));
                ll = 1;
                for jj = int_loc(1) : int_loc(2)
                    index = find(M(:,time_index)<time_intervals(jj+1),1,'last');
                    S = A(ll:index,:);
                    pixel_width_x = options.FOVa_x/options.Nx;
                    pixel_width_y = options.FOVa_y/options.Ny;
                    z_width = options.axial_fov/options.Nz;
                    x = single((-pixel_width_x*(options.Nx-1)/2:pixel_width_x:pixel_width_x*(options.Nx-1)/2)');
                    y = single((-pixel_width_y*(options.Ny-1)/2:pixel_width_y:pixel_width_y*(options.Ny-1)/2)');
                    z = single((-z_width*(options.Nz-1)/2:z_width:z_width*(options.Nz-1)/2)');
                    max_x = max(x) + pixel_width_x/2;
                    max_y = max(y) + pixel_width_y/2;
                    max_z = max(z) + z_width/2;
                    
                    CC = ones(length(S),1,'int16');
                    t1 = all(S(:,1:3)==S(:,4:6),2);
                    t2 = all([abs(S(:,1))<=max_x abs(S(:,2))<=max_y abs(S(:,3))<=max_z],2);
                    t3 = all([abs(S(:,4))<=max_x abs(S(:,5))<=max_y abs(S(:,6))<=max_z],2);
                    t = (t1+t2+t3==3);
                    CC = CC(t);
                    [~,i] = min(abs(bsxfun(@minus,x.',S(t,1))),[],2);
                    [~,j] = min(abs(bsxfun(@minus,y.',S(t,2))),[],2);
                    [~,k] = min(abs(bsxfun(@minus,z.',S(t,3))),[],2);
                    FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                    %                 FOV = circshift(FOV, round(size(FOV,1)/6));
                    C{1,jj} = C{1,jj} + FOV;
                    
                    CC = ones(length(S),1,'int16');
                    t1 = true(size(S,1),1);
                    t2 = all([abs(S(:,1))<=max_x abs(S(:,2))<=max_y abs(S(:,3))<=max_z],2);
                    t = (t1+t2==2);
                    CC = CC(t);
                    [~,i] = min(abs(bsxfun(@minus,x.',S(t,1))),[],2);
                    [~,j] = min(abs(bsxfun(@minus,y.',S(t,2))),[],2);
                    [~,k] = min(abs(bsxfun(@minus,z.',S(t,3))),[],2);
                    FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                    %                 FOV = circshift(FOV, round(size(FOV,1)/6));
                    C{2,jj} = C{2,jj} + FOV;
                    
                    CC = ones(length(S),1,'int16');
                    t1 = true(size(S,4),1);
                    t2 = all([abs(S(:,4))<=max_x abs(S(:,5))<=max_y abs(S(:,6))<=max_z],2);
                    t = (t1+t2==2);
                    CC = CC(t);
                    [~,i] = min(abs(bsxfun(@minus,x.',S(t,4))),[],2);
                    [~,j] = min(abs(bsxfun(@minus,y.',S(t,5))),[],2);
                    [~,k] = min(abs(bsxfun(@minus,z.',S(t,6))),[],2);
                    FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                    %                 FOV = circshift(FOV, round(size(FOV,1)/6));
                    C{3,jj} = C{3,jj} + FOV;
                    
                    C{4,jj} = C{2,jj} + C{3,jj};
                    
                    CC = ones(length(S),1,'int16');
                    t1 = true(size(S,1),1);
                    t2 = all([abs(S(:,1))<=max_x abs(S(:,2))<=max_y abs(S(:,3))<=max_z],2);
                    t3 = all([abs(S(:,4))<=max_x abs(S(:,5))<=max_y abs(S(:,6))<=max_z],2);
                    t = (t1+t2+t3==3);
                    CC = CC(t);
                    [~,i] = min(abs(bsxfun(@minus,(x).',S(t,1))),[],2);
                    [~,j] = min(abs(bsxfun(@minus,(y).',S(t,2))),[],2);
                    [~,k] = min(abs(bsxfun(@minus,(z).',S(t,3))),[],2);
                    [~,i2] = min(abs(bsxfun(@minus,(x).',S(t,4))),[],2);
                    [~,j2] = min(abs(bsxfun(@minus,(y).',S(t,5))),[],2);
                    [~,k2] = min(abs(bsxfun(@minus,(z).',S(t,6))),[],2);
                    i = [i;i2];
                    j = [j;j2];
                    k = [k;k2];
                    CC = repmat(CC,uint32(2),uint32(1));
                    FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                    %                 FOV = circshift(FOV, round(size(FOV,1)/6));
                    FOV = FOV/2;
                    C{5,jj} = C{5,jj} + FOV;
                    
                    CC = ones(length(S),1,'int16');
                    t1 = true(size(S,1),1);
                    t2 = all([abs(S(:,1))<=max_x abs(S(:,2))<=max_y abs(S(:,3))<=max_z],2);
                    t3 = all([abs(S(:,4))<=max_x abs(S(:,5))<=max_y abs(S(:,6))<=max_z],2);
                    t = (t1+t2+t3==3);
                    CC = CC(t);
                    C1 = mean([S(t,1) S(t,4)],2);
                    C2 = mean([S(t,2) S(t,5)],2);
                    C3 = mean([S(t,3) S(t,6)],2);
                    [~,i] = min(abs(bsxfun(@minus,(x).',C1)),[],2);
                    [~,j] = min(abs(bsxfun(@minus,(y).',C2)),[],2);
                    [~,k] = min(abs(bsxfun(@minus,(z).',C3)),[],2);
                    FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                    %                 FOV = circshift(FOV, round(size(FOV,1)/6));
                    C{6,jj} = C{6,jj} + FOV;
                    
                    if options.obtain_trues
                        CC = ones(length(S),1,'int16');
                        CC = CC(trues_index);
                        SS = S(trues_index,:);
                        t= all([abs(SS(:,1))<=max_x abs(SS(:,2))<=max_y abs(SS(:,3))<=max_z],2);
                        CC = CC(t);
                        [~,i] = min(abs(bsxfun(@minus,x.',SS(t,1))),[],2);
                        [~,j] = min(abs(bsxfun(@minus,y.',SS(t,2))),[],2);
                        [~,k] = min(abs(bsxfun(@minus,z.',SS(t,3))),[],2);
                        if isempty(i)
                            FOV = zeros([Nx Ny Nz],'uint16');
                        else
                            FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                        end
                        % FOV = circshift(FOV, round(size(FOV,1)/6));
                        C{7,jj} = C{7,jj} + FOV;
                        clear SS
                    end
                    if options.store_scatter
                        CC = ones(length(S),1,'int16');
                        CC = CC(scatter_index);
                        SS = S(scatter_index,:);
                        t= all([abs(SS(:,1))<=max_x abs(SS(:,2))<=max_y abs(SS(:,3))<=max_z],2);
                        CC = CC(t);
                        [~,i] = min(abs(bsxfun(@minus,x.',SS(t,1))),[],2);
                        [~,j] = min(abs(bsxfun(@minus,y.',SS(t,2))),[],2);
                        [~,k] = min(abs(bsxfun(@minus,z.',SS(t,3))),[],2);
                        if isempty(i)
                            FOV = zeros([Nx Ny Nz],'uint16');
                        else
                            FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                        end
                        % FOV = circshift(FOV, round(size(FOV,1)/6));
                        SC = SC + FOV;
                        CC = ones(length(S),1,'int16');
                        CC = CC(scatter_index);
                        SS = S(scatter_index,:);
                        t= all([abs(SS(:,4))<=max_x abs(SS(:,5))<=max_y abs(SS(:,6))<=max_z],2);
                        CC = CC(t);
                        [~,i] = min(abs(bsxfun(@minus,x.',SS(t,4))),[],2);
                        [~,j] = min(abs(bsxfun(@minus,y.',SS(t,5))),[],2);
                        [~,k] = min(abs(bsxfun(@minus,z.',SS(t,6))),[],2);
                        if isempty(i)
                            FOV = zeros([Nx Ny Nz],'uint16');
                        else
                            FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                        end
                        % FOV = circshift(FOV, round(size(FOV,1)/6));
                        SC = SC + FOV;
                        clear SS
                    end
                    if options.store_randoms
                        CC = ones(length(S),1,'int16');
                        CC = CC(randoms_index);
                        SS = S(randoms_index,:);
                        t= all([abs(SS(:,1))<=max_x abs(SS(:,2))<=max_y abs(SS(:,3))<=max_z],2);
                        CC = CC(t);
                        [~,i] = min(abs(bsxfun(@minus,x.',SS(t,1))),[],2);
                        [~,j] = min(abs(bsxfun(@minus,y.',SS(t,2))),[],2);
                        [~,k] = min(abs(bsxfun(@minus,z.',SS(t,3))),[],2);
                        if isempty(i)
                            FOV = zeros([Nx Ny Nz],'uint16');
                        else
                            FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                        end
                        % FOV = circshift(FOV, round(size(FOV,1)/6));
                        RA = RA + FOV;
                        CC = ones(length(S),1,'int16');
                        CC = CC(randoms_index);
                        t= all([abs(SS(:,4))<=max_x abs(SS(:,5))<=max_y abs(SS(:,6))<=max_z],2);
                        CC = CC(t);
                        [~,i] = min(abs(bsxfun(@minus,x.',SS(t,4))),[],2);
                        [~,j] = min(abs(bsxfun(@minus,y.',SS(t,5))),[],2);
                        [~,k] = min(abs(bsxfun(@minus,z.',SS(t,6))),[],2);
                        if isempty(i)
                            FOV = zeros([Nx Ny Nz],'uint16');
                        else
                            FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                        end
                        % FOV = circshift(FOV, round(size(FOV,1)/6));
                        RA = RA + FOV;
                        clear SS
                    end
                    ll = ll + index;
                end
            end
            
            % The number of crystal ring of each coincidence event (e.g. the first single hits
            % a detector on crystal ring 5 and second a detector on crystal ring 10)
            % [0 rings - 1]
            if int_loc(1) > 0
                if (ascii_ind.module_ind1 == 0 || ascii_ind.module_ind2 == 0) && options.linear_multip > 1
                    ring_number1 = uint16(floor(M(:,rsector_ind1) / blocks_per_ring) * cryst_per_block + floor(M(:,crs_ind1)/cryst_per_block));
                    ring_number2 = uint16(floor(M(:,rsector_ind2) / blocks_per_ring) * cryst_per_block + floor(M(:,crs_ind2)/cryst_per_block));
                elseif (ascii_ind.module_ind1 == 0 || ascii_ind.module_ind2 == 0) && options.linear_multip == 1
                    ring_number1 = uint16(M(:,rsector_ind1));
                    ring_number2 = uint16(M(:,rsector_ind2));
                elseif options.linear_multip == 1
                    ring_number1 = uint16(M(:,ascii_ind.module_ind1));
                    ring_number2 = uint16(M(:,ascii_ind.module_ind2));
                else
                    ring_number1 = uint16(mod(M(:,ascii_ind.module_ind1),options.linear_multip) * cryst_per_block + floor(M(:,crs_ind1)/cryst_per_block));
                    ring_number2 = uint16(mod(M(:,ascii_ind.module_ind2),options.linear_multip) * cryst_per_block + floor(M(:,crs_ind2)/cryst_per_block));
                end
                
                % detector number of the single at the above ring (e.g. first single hits a
                % detector number 54 on ring 5)
                % [0 det_per_ring - 1]
                ring_pos1 = uint16(mod(M(:,rsector_ind1),blocks_per_ring)*cryst_per_block+mod(M(:,crs_ind1),cryst_per_block));
                ring_pos2 = uint16(mod(M(:,rsector_ind2),blocks_per_ring)*cryst_per_block+mod(M(:,crs_ind2),cryst_per_block));
                
                % detector number (MATLAB numbering) of each single
                % [1 det_per_ring * rings]
                if partitions == 1
                    LL1 = ring_number1*uint16(det_per_ring) + ring_pos1 + 1;
                    LL2 = ring_number2*uint16(det_per_ring) + ring_pos2 + 1;
                    prompts = prompts + accumarray([LL1 LL2],uint16(1),[detectors detectors],@(x) sum(x,'native'));
                    if options.obtain_trues
                        if any(trues_index)
                            trues = trues + accumarray([LL1(trues_index) LL2(trues_index)],uint16(1),[detectors detectors],@(x) sum(x,'native'));
                        end
                    end
                    if options.store_scatter
                        if any(scatter_index)
                            scatter = scatter + accumarray([LL1(scatter_index) LL2(scatter_index)],uint16(1),[detectors detectors],@(x) sum(x,'native'));
                        end
                    end
                    if options.store_randoms
                        if any(randoms_index)
                            randoms = randoms + accumarray([LL1(randoms_index) LL2(randoms_index)],uint16(1),[detectors detectors],@(x) sum(x,'native'));
                        end
                    end
                    li = li + length(LL1);
                else
                    LL1 = ring_number1*uint16(det_per_ring) + ring_pos1 + 1;
                    LL2 = ring_number2*uint16(det_per_ring) + ring_pos2 + 1;
                    ll = 1;
                    for kk = int_loc(1) : int_loc(2)
                        index = find(M(:,time_index)<time_intervals(kk+1),1,'last');
                        if isempty(prompts{kk})
                            prompts{kk} = accumarray([LL1(ll:index) LL2(ll:index)],1,[detectors detectors],[],[],true);
                        else
                            prompts{kk} = prompts{kk} + accumarray([LL1(ll:index) LL2(ll:index)],1,[detectors detectors],[],[],true);
                        end
                        if options.obtain_trues
                            apu_trues = trues_index(ll:index);
                            L1trues = LL1(ll:index);
                            L1trues = L1trues(apu_trues);
                            L2trues = LL2(ll:index);
                            L2trues = L2trues(apu_trues);
                            if isempty(trues{kk})
                                trues{kk} = accumarray([L1trues L2trues],1,[detectors detectors],[],[],true);
                            else
                                trues{kk} = trues{kk} + accumarray([L1trues L2trues],1,[detectors detectors],[],[],true);
                            end
                        end
                        if options.store_scatter
                            apu_scatter = scatter_index(ll:index);
                            L1scatter = LL1(ll:index);
                            L1scatter = L1scatter(apu_scatter);
                            L2scatter = LL2(ll:index);
                            L2scatter = L2scatter(apu_scatter);
                            if isempty(scatter{kk})
                                scatter{kk} = accumarray([L1scatter L2scatter],1,[detectors detectors],[],[],true);
                            else
                                scatter{kk} = scatter{kk} + accumarray([L1scatter L2scatter],1,[detectors detectors],[],[],true);
                            end
                        end
                        if options.store_randoms
                            apu_randoms = randoms_index(ll:index);
                            L1randoms = LL1(ll:index);
                            L1randoms = L1randoms(apu_randoms);
                            L2randoms = LL2(ll:index);
                            L2randoms = L2randoms(apu_randoms);
                            if isempty(randoms{kk})
                                randoms{kk} = accumarray([L1randoms L2randoms],1,[detectors detectors],[],[],true);
                            else
                                randoms{kk} = randoms{kk} + accumarray([L1randoms L2randoms],1,[detectors detectors],[],[],true);
                            end
                        end
                        ll = index + 1;
                    end
                end
            end
            
            if options.verbose
                disp(['File ' fnames(lk).name ' loaded'])
            end
            
            if options.randoms_correction && length(delay_names) >= lk
                if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','9.6') || exist('OCTAVE_VERSION','builtin') == 5
                    M = importdata([fpath delay_names(lk).name]);
                    % Check for corrupted data
                    if any(any(isnan(M))) > 0
                        header = size(M,1);
                        warning(['Line ' num2str(header) ' corrupted, skipping.'])
                        B = importdata([fpath delay_names(lk).name],' ', header);
                        while any(any(isnan(B.data))) > 0
                            header = header + size(B.data,1);
                            M = [M; B.data(1:end-1,:)];
                            B = importdata([fpath delay_names(lk).name],' ', header);
                            warning(['Line ' num2str(header) ' corrupted, skipping.'])
                        end
                        M = [M; B.data];
                        clear B
                    end
                else
                    M = readmatrix([fpath delay_names(lk).name]);
                    % Check for corrupted data
                    if sum(isnan(M(:,end))) > 0 && sum(isnan(M(:,end))) > round(size(M,1)/2)
                        [rows, ~] = find(~isnan(M(:,end)));
                        for kg = 1 : length(rows)
                            warning(['Line ' num2str(rows(kg)) ' corrupted, skipping.'])
                            M(rows(kg),:) = [];
                        end
                    else
                        if any(any(isnan(M))) > 0
                            [rows, ~] = find(isnan(M));
                            rows = unique(rows);
                            for kg = 1 : length(rows)
                                warning(['Line ' num2str(rows(kg)) ' corrupted, skipping.'])
                                M(rows(kg),:) = [];
                            end
                        end
                    end
                end
                
                if partitions > 1
                    if isempty(find(M(:,time_index) > time_intervals(1),1,'last'))
                        int_loc(1) = 0;
                    elseif isempty(find(M(1,time_index) > time_intervals,1,'last'))
                        int_loc(1) = 1;
                    else
                        int_loc(1) = find(M(1,time_index) > time_intervals,1,'last');
                    end
                    if isempty(find(M(:,time_index) < time_intervals(end),1,'first'))
                        int_loc(2) = 0;
                    elseif isempty(find(M(end,time_index) < time_intervals,1,'first'))
                        int_loc(2) = length(time_intervals) - 1;
                    else
                        int_loc(2) = find(M(end,time_index) < time_intervals,1,'first') - 1;
                    end
                else
                    int_loc(1) = 1;
                    int_loc(2) = 1;
                end
                
                if int_loc(1) > 0
                    if options.module_ind1 == 0 || options.module_ind2 == 0
                        ring_number1 = uint16(floor(M(:,rsector_ind1)/blocks_per_ring)*cryst_per_block+floor(M(:,crs_ind1)/cryst_per_block));
                        ring_number2 = uint16(floor(M(:,rsector_ind2)/blocks_per_ring)*cryst_per_block+floor(M(:,crs_ind2)/cryst_per_block));
                    else
                        ring_number1 = uint16(mod(M(:,options.module_ind1),options.linear_multip)*cryst_per_block+floor(M(:,crs_ind1)/cryst_per_block));
                        ring_number2 = uint16(mod(M(:,options.module_ind2),options.linear_multip)*cryst_per_block+floor(M(:,crs_ind2)/cryst_per_block));
                    end
                    
                    % detector number of the single at the above ring (e.g. first single hits a
                    % detector number 54 on ring 5)
                    % [0 det_per_ring - 1]
                    ring_pos1 = uint16(mod(M(:,rsector_ind1),blocks_per_ring)*cryst_per_block+mod(M(:,crs_ind1),cryst_per_block));
                    ring_pos2 = uint16(mod(M(:,rsector_ind2),blocks_per_ring)*cryst_per_block+mod(M(:,crs_ind2),cryst_per_block));
                    
                    % detector number (MATLAB numbering) of each single
                    % [1 det_per_ring * rings]
                    if partitions == 1
                        Ldelay1 = ring_number1*uint16(det_per_ring) + ring_pos1 + 1;
                        Ldelay2 = ring_number2*uint16(det_per_ring) + ring_pos2 + 1;
                        delays = delays + accumarray([Ldelay1 Ldelay2],uint16(1),[detectors detectors],@(x) sum(x,'native'));
                    else
                        LL1 = ring_number1*uint16(det_per_ring)+ring_pos1+1;
                        LL2 = ring_number2*uint16(det_per_ring)+ring_pos2+1;
                        ll = 1;
                        for kk = int_loc(1) : int_loc(2)
                            index = find(M(:,time_index)<time_intervals(kk+1),1,'last');
                            if isempty(delays{kk})
                                delays{kk} = accumarray([LL1(ll:index) LL2(ll:index)],1,[detectors detectors],[],[],true);
                            else
                                delays{kk} = delays{kk} + accumarray([LL1(ll:index) LL2(ll:index)],1,[detectors detectors],[],[],true);
                            end
                            ll = index + 1;
                        end
                    end
                end
                
                
                if options.verbose
                    disp(['File ' delay_names(lk).name ' loaded'])
                end
            end
            
        end
        clear randoms_index scatter_index apu_randoms L1randoms L2randoms apu_scatter L2scatter L1scatter L2trues L1trues apu_trues trues_index...
            M FOV CC i j k t t1 t2 t3 A S
        if source_index1 ~= 0 && ~isempty(source_index1) && source_index2 ~= 0 && ~isempty(source_index2) && options.source
            save([machine_name '_Ideal_image_coordinates_' name '_ASCII.mat'],'C')
            if options.store_randoms
                save([machine_name '_Ideal_image_coordinates_' name '_ASCII.mat'],'RA','-append')
            end
            if options.store_scatter
                save([machine_name '_Ideal_image_coordinates_' name '_ASCII.mat'],'SC','-append')
            end
        end
        
        % forms a sparse matrix containing all the events at different LORs
        % e.g. a value at [325 25100] contains the number of events detected at
        % detectors 325 (detector 1) and 25100 (detector 2)
        % a value at [25100 325] on the other hand tells the number of events when
        % 25100 is detector 1 and 325 detector 2
        %         if partitions == 1
        %             if options.obtain_trues
        %                 trues = accumarray([LL1(Ltrues) LL2(Ltrues)],uint16(1),[detectors detectors],@(x) sum(x,'native'));
        %             end
        %             if options.store_scatter
        %                 scatter = accumarray([LL1(Lscatter) LL2(Lscatter)],uint16(1),[detectors detectors],@(x) sum(x,'native'));
        %             end
        %             if options.store_randoms
        %                 randoms = accumarray([LL1(Lrandoms) LL2(Lrandoms)],uint16(1),[detectors detectors],@(x) sum(x,'native'));
        %             end
        %             if options.randoms_correction
        %                 delays = accumarray([Ldelay1 Ldelay2],uint16(1),[detectors detectors],@(x) sum(x,'native'));
        %             end
        %         end
        clear Lrandoms Lscatter Ltrues LL1 LL2 C Ldelay2 Ldelay1
        %     return
        
    elseif options.use_root
        %%
        blocks_per_ring = uint32(blocks_per_ring);
        cryst_per_block = uint32(cryst_per_block);
        det_per_ring = uint32(det_per_ring);
        scatter_components = logical(options.scatter_components);
        
        if partitions > 1
            save_string_raw = [machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_root.mat'];
        else
            save_string_raw = [machine_name '_measurements_' name '_static_raw_root.mat'];
        end
        
        
        % Take coincidence data
        fnames = dir([fpath '*.root']);
        if size(fnames,1) == 0
            fpath = pwd;
            fpath = [fpath '/'];
            fnames = dir([fpath '*.root']);
            if size(fnames,1) == 0
                error('No ROOT (.root) files were found. Check your filepath (options.fpath) or current folder.')
            end
        end
        
        if options.partitions > 1
            prompts = cell(partitions,1);
            if options.obtain_trues
                trues = cell(partitions,1);
            end
            if options.store_scatter
                scatter = cell(partitions,1);
            end
            if options.store_randoms
                randoms = cell(partitions,1);
            end
            if options.randoms_correction
                delays = cell(partitions,1);
            end
        else
            prompts = zeros(detectors, detectors, 'uint16');
            if options.obtain_trues
                trues = zeros(detectors, detectors, 'uint16');
            end
            if options.store_scatter
                scatter = zeros(detectors, detectors, 'uint16');
            end
            if options.store_randoms
                randoms = zeros(detectors, detectors, 'uint16');
            end
            if options.randoms_correction
                delays = zeros(detectors, detectors, 'uint16');
            end
        end
        if source
            if options.obtain_trues
                C = cell(7,partitions);
            else
                C = cell(6,partitions);
            end
            C(:) = {zeros(options.Nx, options.Ny, options.Nz, 'uint16')};
            if options.store_scatter
                SC = zeros(options.Nx, options.Ny, options.Nz, 'uint16');
            end
            if options.store_randoms
                RA = zeros(options.Nx, options.Ny, options.Nz, 'uint16');
            end
        end
        
        if options.verbose
            disp('Data load started')
        end
        
        % Go through all the files
        for lk=1:length(fnames)
            
            
            nimi = [fpath fnames(lk).name];
            
            % Use the new features of MATLAB if version is R2019a or newer
            if exist('OCTAVE_VERSION','builtin') == 0 && verLessThan('matlab','9.6') || exist('OCTAVE_VERSION','builtin') == 5
                [L1, L2, tpoints, A, int_loc, Ltrues, Lscatter, Lrandoms, trues_index, Ldelay1, Ldelay2, int_loc_delay, tpoints_delay, randoms_index, scatter_index] = GATE_root_matlab_C(nimi,vali,alku,loppu, ...
                    uint32(detectors), blocks_per_ring, cryst_per_block, det_per_ring, uint32(options.linear_multip), source, time_intervals, options.obtain_trues, ...
                    options.store_scatter, options.store_randoms, scatter_components, options.randoms_correction);
            else
                [L1, L2, tpoints, A, int_loc, Ltrues, Lscatter, Lrandoms, trues_index, Ldelay1, Ldelay2, int_loc_delay, tpoints_delay, randoms_index, scatter_index] ...
                    = GATE_root_matlab_MEX(nimi,vali,alku,loppu, detectors, blocks_per_ring, cryst_per_block, det_per_ring, source, time_intervals, options, scatter_components);
                
            end
            
%             int_loc = int_loc + 1;
%             
%             if options.randoms_correction
%                 int_loc_delay = int_loc_delay + 1;
%             end
            
            
            % An ideal image can be formed with this. Takes the variables
            % containing the source locations for both singles
            ll = 1;
            if source
                for jj = int_loc(1) : min(int_loc(2),partitions)
                    if partitions > 1
                        S = A(tpoints(ll) + 1:tpoints(ll+1),:);
                        if options.obtain_trues
                            trues_index = logical(Ltrues(tpoints(ll) + 1:tpoints(ll+1),:));
                            %                             inde = any(S,2);
                            %                             trues_index = trues_index(inde);
                        end
                        if options.store_scatter
                            scatter_index = logical(Lscatter(tpoints(ll) + 1:tpoints(ll+1),:));
                            %                             inde = any(S,2);
                            %                             scatter_index = scatter_index(inde);
                        end
                        if options.store_randoms
                            randoms_index = logical(Lrandoms(tpoints(ll) + 1:tpoints(ll+1),:));
                            %                             inde = any(S,2);
                            %                             randoms_index = randoms_index(inde);
                        end
                        ll = ll + 1;
                    else
                        S = single(A);
                    end
                    pixel_width_x = options.FOVa_x/options.Nx;
                    pixel_width_y = options.FOVa_y/options.Ny;
                    z_width = options.axial_fov/options.Nz;
                    x = single((-pixel_width_x*(options.Nx-1)/2:pixel_width_x:pixel_width_x*(options.Nx-1)/2)');
                    y = single((-pixel_width_y*(options.Ny-1)/2:pixel_width_y:pixel_width_y*(options.Ny-1)/2)');
                    z = single((-z_width*(options.Nz-1)/2:z_width:z_width*(options.Nz-1)/2)');
                    max_x = max(x) + pixel_width_x/2;
                    max_y = max(y) + pixel_width_y/2;
                    max_z = max(z) + z_width/2;
                    
                    CC = ones(length(S),1,'int16');
                    t1 = all(S(:,1:3)==S(:,4:6),2);
                    t2 = all([abs(S(:,1))<=max_x abs(S(:,2))<=max_y abs(S(:,3))<=max_z],2);
                    t3 = all([abs(S(:,4))<=max_x abs(S(:,5))<=max_y abs(S(:,6))<=max_z],2);
                    t = (t1+t2+t3==3);
                    CC = CC(t);
                    [~,i] = min(abs(bsxfun(@minus,x.',S(t,1))),[],2);
                    [~,j] = min(abs(bsxfun(@minus,y.',S(t,2))),[],2);
                    [~,k] = min(abs(bsxfun(@minus,z.',S(t,3))),[],2);
                    FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                    %             FOV = circshift(FOV, round(size(FOV,1)/6));
                    C{1,jj} = C{1,jj} + FOV;
                    
                    CC = ones(length(S),1,'int16');
                    t1 = true(size(S,1),1);
                    t2 = all([abs(S(:,1))<=max_x abs(S(:,2))<=max_y abs(S(:,3))<=max_z],2);
                    t = (t1+t2==2);
                    CC = CC(t);
                    [~,i] = min(abs(bsxfun(@minus,x.',S(t,1))),[],2);
                    [~,j] = min(abs(bsxfun(@minus,y.',S(t,2))),[],2);
                    [~,k] = min(abs(bsxfun(@minus,z.',S(t,3))),[],2);
                    FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                    %             FOV = circshift(FOV, round(size(FOV,1)/6));
                    C{2,jj} = C{2,jj} + FOV;
                    
                    CC = ones(length(S),1,'int16');
                    t1 = true(size(S,4),1);
                    t2 = all([abs(S(:,4))<=max_x abs(S(:,5))<=max_y abs(S(:,6))<=max_z],2);
                    t = (t1+t2==2);
                    CC = CC(t);
                    [~,i] = min(abs(bsxfun(@minus,x.',S(t,4))),[],2);
                    [~,j] = min(abs(bsxfun(@minus,y.',S(t,5))),[],2);
                    [~,k] = min(abs(bsxfun(@minus,z.',S(t,6))),[],2);
                    FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                    %             FOV = circshift(FOV, round(size(FOV,1)/6));
                    C{3,jj} = C{3,jj} + FOV;
                    
                    C{4,jj} = C{2,jj} + C{3,jj};
                    
                    CC = ones(length(S),1,'int16');
                    t1 = true(size(S,1),1);
                    t2 = all([abs(S(:,1))<=max_x abs(S(:,2))<=max_y abs(S(:,3))<=max_z],2);
                    t3 = all([abs(S(:,4))<=max_x abs(S(:,5))<=max_y abs(S(:,6))<=max_z],2);
                    t = (t1+t2+t3==3);
                    CC = CC(t);
                    [~,i] = min(abs(bsxfun(@minus,(x).',S(t,1))),[],2);
                    [~,j] = min(abs(bsxfun(@minus,(y).',S(t,2))),[],2);
                    [~,k] = min(abs(bsxfun(@minus,(z).',S(t,3))),[],2);
                    [~,i2] = min(abs(bsxfun(@minus,(x).',S(t,4))),[],2);
                    [~,j2] = min(abs(bsxfun(@minus,(y).',S(t,5))),[],2);
                    [~,k2] = min(abs(bsxfun(@minus,(z).',S(t,6))),[],2);
                    i = [i;i2];
                    j = [j;j2];
                    k = [k;k2];
                    CC = repmat(CC,uint32(2),uint32(1));
                    FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                    %             FOV = circshift(FOV, round(size(FOV,1)/6));
                    FOV = FOV/2;
                    C{5,jj} = C{5,jj} + FOV;
                    
                    CC = ones(length(S),1,'int16');
                    t1 = true(size(S,1),1);
                    t2 = all([abs(S(:,1))<=max_x abs(S(:,2))<=max_y abs(S(:,3))<=max_z],2);
                    t3 = all([abs(S(:,4))<=max_x abs(S(:,5))<=max_y abs(S(:,6))<=max_z],2);
                    t = (t1+t2+t3==3);
                    CC = CC(t);
                    C1 = mean([S(t,1) S(t,4)],2);
                    C2 = mean([S(t,2) S(t,5)],2);
                    C3 = mean([S(t,3) S(t,6)],2);
                    [~,i] = min(abs(bsxfun(@minus,(x).',C1)),[],2);
                    [~,j] = min(abs(bsxfun(@minus,(y).',C2)),[],2);
                    [~,k] = min(abs(bsxfun(@minus,(z).',C3)),[],2);
                    FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                    %             FOV = circshift(FOV, round(size(FOV,1)/6));
                    C{6,jj} = C{6,jj} + FOV;
                    
                    if options.obtain_trues
                        CC = ones(nnz(trues_index),1,'int16');
                        SS = S(trues_index,:);
                        t= all([abs(SS(:,1))<=max_x abs(SS(:,2))<=max_y abs(SS(:,3))<=max_z],2);
                        CC = CC(t);
                        [~,i] = min(abs(bsxfun(@minus,x.',SS(t,1))),[],2);
                        [~,j] = min(abs(bsxfun(@minus,y.',SS(t,2))),[],2);
                        [~,k] = min(abs(bsxfun(@minus,z.',SS(t,3))),[],2);
                        if isempty(i)
                            FOV = zeros([Nx Ny Nz],'uint16');
                        else
                            FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                        end
                        %                     FOV = circshift(FOV, round(size(FOV,1)/6));
                        C{7,jj} = C{7,jj} + FOV;
                    end
                    if options.store_scatter
                        CC = ones(nnz(scatter_index),1,'int16');
                        SS = S(scatter_index,:);
                        t= all([abs(SS(:,1))<=max_x abs(SS(:,2))<=max_y abs(SS(:,3))<=max_z],2);
                        CC = CC(t);
                        [~,i] = min(abs(bsxfun(@minus,x.',SS(t,1))),[],2);
                        [~,j] = min(abs(bsxfun(@minus,y.',SS(t,2))),[],2);
                        [~,k] = min(abs(bsxfun(@minus,z.',SS(t,3))),[],2);
                        if isempty(i)
                            FOV = zeros([Nx Ny Nz],'uint16');
                        else
                            FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                        end
                        % FOV = circshift(FOV, round(size(FOV,1)/6));
                        SC = SC + FOV;
                        CC = ones(length(S),1,'int16');
                        CC = CC(scatter_index);
                        SS = S(scatter_index,:);
                        t= all([abs(SS(:,4))<=max_x abs(SS(:,5))<=max_y abs(SS(:,6))<=max_z],2);
                        CC = CC(t);
                        [~,i] = min(abs(bsxfun(@minus,x.',SS(t,4))),[],2);
                        [~,j] = min(abs(bsxfun(@minus,y.',SS(t,5))),[],2);
                        [~,k] = min(abs(bsxfun(@minus,z.',SS(t,6))),[],2);
                        if isempty(i)
                            FOV = zeros([Nx Ny Nz],'uint16');
                        else
                            FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                        end
                        % FOV = circshift(FOV, round(size(FOV,1)/6));
                        SC = SC + FOV;
                        clear SS
                    end
                    if options.store_randoms
                        CC = ones(nnz(randoms_index),1,'int16');
                        SS = S(randoms_index,:);
                        t= all([abs(SS(:,1))<=max_x abs(SS(:,2))<=max_y abs(SS(:,3))<=max_z],2);
                        CC = CC(t);
                        [~,i] = min(abs(bsxfun(@minus,x.',SS(t,1))),[],2);
                        [~,j] = min(abs(bsxfun(@minus,y.',SS(t,2))),[],2);
                        [~,k] = min(abs(bsxfun(@minus,z.',SS(t,3))),[],2);
                        if isempty(i)
                            FOV = zeros([Nx Ny Nz],'uint16');
                        else
                            FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                        end
                        % FOV = circshift(FOV, round(size(FOV,1)/6));
                        RA = RA + FOV;
                        CC = ones(length(S),1,'int16');
                        CC = CC(randoms_index);
                        t= all([abs(SS(:,4))<=max_x abs(SS(:,5))<=max_y abs(SS(:,6))<=max_z],2);
                        CC = CC(t);
                        [~,i] = min(abs(bsxfun(@minus,x.',SS(t,1))),[],2);
                        [~,j] = min(abs(bsxfun(@minus,y.',SS(t,2))),[],2);
                        [~,k] = min(abs(bsxfun(@minus,z.',SS(t,3))),[],2);
                        if isempty(i)
                            FOV = zeros([Nx Ny Nz],'uint16');
                        else
                            FOV = (rot90(uint16(accumarray([uint32(i) uint32(j) uint32(k)],CC,[Nx Ny Nz])),1));
                        end
                        % FOV = circshift(FOV, round(size(FOV,1)/6));
                        RA = RA + FOV;
                        clear SS
                    end
                end
            end
            
            
            if partitions == 1
                prompts = prompts + L1;
                if options.obtain_trues
                    trues = trues + Ltrues;
                end
                if options.store_scatter
                    scatter = scatter + Lscatter;
                end
                if options.store_randoms
                    randoms = randoms + Lrandoms;
                end
                if options.randoms_correction
                    delays = delays + Ldelay1;
                end
            else
                if  int_loc(1) > 0
                    ll = 1;
                    for kk = int_loc(1) : min(int_loc(2),partitions)
                        apu1 = L1(tpoints(ll) + 1:tpoints(ll+1));
                        apu2 = L2(tpoints(ll) + 1:tpoints(ll+1));
                        inde = any(apu1,2);
                        apu1 = apu1(inde,:);
                        apu2 = apu2(inde,:);
                        if isempty(prompts{kk})
                            prompts{kk} = accumarray([apu1 apu2],1,[detectors detectors],[],[],true);
                        else
                            prompts{kk} = prompts{kk} + accumarray([apu1 apu2],1,[detectors detectors],[],[],true);
                        end
                        if options.obtain_trues
                            trues_index = logical(Ltrues(tpoints(ll) + 1:tpoints(ll+1)));
                            trues_index = trues_index(inde);
                            if isempty(trues{kk})
                                trues{kk} = accumarray([apu1(trues_index) apu2(trues_index)],1,[detectors detectors],[],[],true);
                            else
                                trues{kk} = trues{kk} + accumarray([apu1(trues_index) apu2(trues_index)],1,[detectors detectors],[],[],true);
                            end
                        end
                        if options.store_scatter
                            trues_index = logical(Lscatter(tpoints(ll) + 1:tpoints(ll+1)));
                            trues_index = trues_index(inde);
                            if isempty(scatter{kk})
                                scatter{kk} = accumarray([apu1(trues_index) apu2(trues_index)],1,[detectors detectors],[],[],true);
                            else
                                scatter{kk} = scatter{kk} + accumarray([apu1(trues_index) apu2(trues_index)],1,[detectors detectors],[],[],true);
                            end
                        end
                        if options.store_randoms
                            trues_index = logical(Lrandoms(tpoints(ll) + 1:tpoints(ll+1)));
                            trues_index = trues_index(inde);
                            if isempty(randoms{kk})
                                randoms{kk} = accumarray([apu1(trues_index) apu2(trues_index)],1,[detectors detectors],[],[],true);
                            else
                                randoms{kk} = randoms{kk} + accumarray([apu1(trues_index) apu2(trues_index)],1,[detectors detectors],[],[],true);
                            end
                        end
                        ll = ll + 1;
                    end
                end
                if options.randoms_correction
                    if  int_loc_delay(1) > 0
                        ll = 1;
                        for kk = int_loc_delay(1) : min(int_loc_delay(2),partitions)
                            apu1 = Ldelay1(tpoints_delay(ll) + 1:tpoints_delay(ll+1));
                            apu2 = Ldelay2(tpoints_delay(ll) + 1:tpoints_delay(ll+1));
                            inde = any(apu1,2);
                            apu1 = apu1(inde,:);
                            apu2 = apu2(inde,:);
                            if isempty(delays{kk})
                                delays{kk} = accumarray([apu1 apu2],1,[detectors detectors],[],[],true);
                            else
                                delays{kk} = delays{kk} + accumarray([apu1 apu2],1,[detectors detectors],[],[],true);
                            end
                            ll = ll + 1;
                        end
                    end
                end
            end
            if options.verbose
                disp(['File ' fnames(lk).name ' loaded'])
            end
            
        end
        clear Lrandoms apu2 Lscatter Ltrues apu1 trues_index FOV CC i j k t t1 t2 t3 A S inde
        if source
            save([machine_name '_Ideal_image_coordinates_' name '_ROOT.mat'],'C')
            if options.store_randoms
                save([machine_name '_Ideal_image_coordinates_' name '_ROOT.mat'],'RA','-append')
            end
            if options.store_scatter
                save([machine_name '_Ideal_image_coordinates_' name '_ROOT.mat'],'SC','-append')
            end
        end
        clear C GATE_root_matlab_MEX.m
        
        
    end
end
    %%
    
    coincidences = cell(partitions,1);
    if options.use_machine == 0 && options.obtain_trues
        true_coincidences = cell(partitions,1);
    end
    if options.use_machine == 0 && options.store_scatter
        scattered_coincidences = cell(partitions,1);
    end
    if options.use_machine == 0 && options.store_randoms
        random_coincidences = cell(partitions,1);
    end
    if options.randoms_correction && (~options.use_LMF || options.use_machine == 1)
        delayed_coincidences = cell(partitions,1);
    end
    
    tot_counts = 0;
    %     dis = 0;
    if options.use_machine == 0 && options.obtain_trues
        tot_trues = 0;
    end
    if options.use_machine == 0 && options.store_scatter
        tot_scatter = 0;
    end
    if options.use_machine == 0 && options.store_randoms
        tot_randoms = 0;
    end
    if options.randoms_correction && (~options.use_LMF || options.use_machine == 1)
        tot_delayed = 0;
    end
    
    for llo=1:partitions
        
        if partitions > 1
            P = uint16(full(prompts{llo}));
            if options.use_machine == 0 && options.obtain_trues
                T = uint16(full(trues{llo}));
            end
            if options.use_machine == 0 && options.store_scatter
                S = uint16(full(scatter{llo}));
            end
            if options.use_machine == 0 && options.store_randoms
                R = uint16(full(randoms{llo}));
            end
            if options.randoms_correction && (~options.use_LMF || options.use_machine == 1)
                D = uint16(full(delays{llo}));
            end
        else
            P = uint16(full(prompts));
            if options.use_machine == 0 && options.obtain_trues
                T = uint16(full(trues));
            end
            if options.use_machine == 0 && options.store_scatter
                S = uint16(full(scatter));
            end
            if options.use_machine == 0 && options.store_randoms
                R = uint16(full(randoms));
            end
            if options.randoms_correction && (~options.use_LMF || options.use_machine == 1)
                D = uint16(full(delays));
            end
        end
        
        % Sum the upper and lower parts together (LOR from detector 1 to 2
        % is the same as LOR from 2 to 1)
        P = tril(P,0) + triu(P,1)';
        % Take only the lower triangular part
        P = (P(tril(true(size(P)), 0)));
        
        if options.verbose
            counts = sum(P(:));
        end
        
        coincidences{llo, 1} = sparse(double(P));
        
        if options.verbose
            disp(['Total number of prompts at time point ' num2str(llo) ' is ' num2str(counts) '.'])
            tot_counts = tot_counts + counts;
        end
        
        if options.use_machine == 0 && options.obtain_trues
            T = tril(T,0) + triu(T,1)';
            T = (T(tril(true(size(T)), 0)));
            true_coincidences{llo, 1} = sparse(double(T));
            if options.verbose
                disp(['Total number of trues at time point ' num2str(llo) ' is ' num2str(sum(T(:))) '.'])
                tot_trues = tot_trues + sum(T(:));
            end
        end
        if options.use_machine == 0 && options.store_scatter
            S = tril(S,0) + triu(S,1)';
            S = (S(tril(true(size(S)), 0)));
            scattered_coincidences{llo, 1} = sparse(double(S));
            if options.verbose
                disp(['Total number of scattered coincidences at time point ' num2str(llo) ' is ' num2str(sum(S(:))) '.'])
                tot_scatter = tot_scatter + sum(S(:));
            end
        end
        if options.use_machine == 0 && options.store_randoms
            R = tril(R,0) + triu(R,1)';
            %                 R = tril(R,0);
            R = (R(tril(true(size(R)), 0)));
            %                 R = R(discard);
            random_coincidences{llo, 1} = sparse(double(R));
            if options.verbose
                disp(['Total number of randoms at time point ' num2str(llo) ' is ' num2str(sum(R(:))) '.'])
                tot_randoms = tot_randoms + sum(R(:));
            end
        end
        if options.randoms_correction && (~options.use_LMF || options.use_machine == 1)
            D = tril(D,0) + triu(D,1)';
            D = (D(tril(true(size(D)), 0)));
            delayed_coincidences{llo, 1} = sparse(double(D));
            if options.verbose
                disp(['Total number of delayed coincidences at time point ' num2str(llo) ' is ' num2str(sum(D(:))) '.'])
                tot_delayed = tot_delayed + sum(D(:));
            end
        end
    end
    if options.verbose
        disp('Saving data...')
    end
    if partitions == 1
      if exist('OCTAVE_VERSION', 'builtin') == 0
        save(save_string_raw, 'coincidences', '-v7.3')
      else
        save(save_string_raw, 'coincidences')
      end
        for variableIndex = 1:length(variableList)
            if exist(variableList{variableIndex},'var')
                save(save_string_raw,variableList{variableIndex},'-append')
            end
        end
    else
        if options.verbose
            disp(['Total number of prompts at all time points is ' num2str(tot_counts) '.'])
            if options.obtain_trues && options.use_machine == 0
                disp(['Total number of trues at all time points is ' num2str(tot_trues) '.'])
            end
            if options.store_scatter && options.use_machine == 0
                disp(['Total number of scattered coincidences at all time points is ' num2str(tot_scatter) '.'])
            end
            if options.store_randoms && options.use_machine == 0
                disp(['Total number of randoms at all time points is ' num2str(tot_randoms) '.'])
            end
            if options.randoms_correction && (~options.use_LMF || options.use_machine == 1)
                disp(['Total number of delayed coincidences at all time points is ' num2str(tot_delayed) '.'])
            end
        end
       if exist('OCTAVE_VERSION', 'builtin') == 0
          save(save_string_raw, 'coincidences', '-v7.3')
       else
          save(save_string_raw, 'coincidences')
        end
        for variableIndex = 1:length(variableList)
            if exist(variableList{variableIndex},'var')
                save(save_string_raw,variableList{variableIndex},'-append')
            end
        end
    end
    if nargout >= 1
        varargout{1} = coincidences;
    end
    if nargout >= 2
        varargout{2} = delayed_coincidences;
    end
    if nargout >= 3
        varargout{3} = true_coincidences;
    end
    if nargout >= 4
        varargout{4} = scattered_coincidences;
    end
    if nargout >= 5
        varargout{5} = random_coincidences;
    end
    if options.verbose
        disp('Measurements loaded and saved')
        toc
    end
end