function coincidences = load_GATE_data(options)
%% LOAD GATE OUTPUT DATA INTO MATLAB
% Loads either ASCII or LMF format data into MATLAB
%
% Can output data for sinogram formation or for pure, raw
% list-mode, data
%
% For raw data, only the accepted LORs are included
%
% For sinogram format, includes all possible LORs that are not duplicates
%
% Simply input the machine and file format specific parameters into the
% function
%
% For raw data, the prepass phase has to be done beforehand
%
% For all cases, the output data is saved in a mat-file
if options.use_ASCII && options.use_LMF && options.use_root
    disp('ASCII, LMF and root selected, using only ASCII')
    options.use_LMF = false;
    options.use_root = false;
elseif options.use_ASCII && options.use_LMF
    disp('Both ASCII and LMF selected, using only ASCII')
    options.use_LMF = false;
elseif options.use_ASCII && options.use_root
    disp('Both ASCII and root selected, using only ASCII')
    options.use_root = false;
elseif options.use_LMF && options.root
    disp('Both LMF and root selected, using only LMF')
    options.use_root = false;
elseif options.use_ASCII == false && options.use_LMF == false && options.use_root == false
    error('Error: Neither ASCII, LMF nor root data selected')
end
if options.verbose
    tic
end

tot_time = options.tot_time;
partitions = options.partitions;
name = options.name;
fpath = options.fpath;
blocks_per_ring = options.blocks_per_ring;
cryst_per_block = options.cryst_per_block;
machine_name = options.machine_name;
rings = options.rings;
time_index = options.time_index;
det_per_ring = options.det_per_ring;

% Total number of detectors
detectors = det_per_ring*rings;

vali = tot_time;
alku = 0;
loppu = tot_time;

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
    det_per_ring = uint32(options.det_per_ring);
    source = options.source;
    
    coincidence_window = coincidence_window * options.clock_time_step;
    
    
    % Take coincidence data
    fnames = dir([fpath '*.ccs']);
    C = [];
    LL1 = [];
    LL2 = [];
    
    
    bina = [];
    for kk = 1 : R_bits
        bina = [bina, '1'];
    end
    R_length = int32(bin2dec(bina));
    bina = [];
    for kk = 1 : M_bits
        bina = [bina, '1'];
    end
    M_length = int32(bin2dec(bina));
    bina = [];
    for kk = 1 : C_bits
        bina = [bina, '1'];
    end
    C_length = int32(bin2dec(bina));
    
    if partitions == 1
        LL1 = zeros(detectors, detectors, 'uint16');
    else
        LL1 = [];
        LL2 = [];
    end
    
    C = cell(5,1);
    C(:) = {zeros(options.Nx, options.Ny, options.Nz, 'uint16')};
    
    % Go through all the files
    for lk=1:length(fnames)
        
        pituus = int64((fnames(lk).bytes - header_bytes) / data_bytes);
        
        nimi = [fpath fnames(lk).name];
        
        if source
            [L1, L2, tpoints, S] = gate_lmf_matlab(nimi,vali,alku,loppu,pituus, uint32(detectors), blocks_per_ring, cryst_per_block, ...
                det_per_ring, uint32(options.linear_multip), header_bytes, data_bytes, R_bits, M_bits, S_bits, C_bits, L_bits, R_length, M_length, C_length, source, ...
                coincidence_window);
        else
            [L1, L2, tpoints, output] = gate_lmf_matlab(nimi,vali,alku,loppu,pituus, uint32(detectors), blocks_per_ring, cryst_per_block, ...
                det_per_ring, uint32(options.linear_multip), header_bytes, data_bytes, R_bits, M_bits, S_bits, C_bits, L_bits, R_length, M_length, C_length, source, ...
                coincidence_window);
        end
        clear mex
        
        % An ideal image can be formed with this, see the file
        % visualize_pet.m
        if source
            S = S(any(S,2),:);
            
            pixel_width = options.FOVa/options.Nx;
            z_width = options.axial_fov/options.Nz;
            x = single((-pixel_width*(options.Nx-1)/2:pixel_width:pixel_width*(options.Nx-1)/2)');
            y = x;
            z = single((-z_width*(options.Nz-1)/2:z_width:z_width*(options.Nz-1)/2)');
            CC = ones(length(S),1,'int16');
            t1 = all(S(:,1:3)==S(:,4:6),2);
            t2 = all([abs(S(:,1))<=max(x) abs(S(:,2))<=max(y) abs(S(:,3))<=max(z)],2);
            t3 = all([abs(S(:,4))<=max(x) abs(S(:,5))<=max(y) abs(S(:,6))<=max(z)],2);
            t = (t1+t2+t3==3);
            CC = CC(t);
            [~,i] = min(abs(bsxfun(@minus,x.',S(t,1))),[],2);
            [~,j] = min(abs(bsxfun(@minus,y.',S(t,2))),[],2);
            [~,k] = min(abs(bsxfun(@minus,z.',S(t,3))),[],2);
            FOV = (rot90(uint16(accumarray([i j k],CC,[options.Nx options.Ny options.Nz])),-2));
            FOV = circshift(FOV, round(size(FOV,1)/6));
            C{1} = C{1} + FOV;
            
            CC = ones(length(S),1,'int16');
            t1 = true(size(S,1),1);
            t2 = all([abs(S(:,1))<=max(x) abs(S(:,2))<=max(y) abs(S(:,3))<=max(z)],2);
            t = (t1+t2==2);
            CC = CC(t);
            [~,i] = min(abs(bsxfun(@minus,x.',S(t,1))),[],2);
            [~,j] = min(abs(bsxfun(@minus,y.',S(t,2))),[],2);
            [~,k] = min(abs(bsxfun(@minus,z.',S(t,3))),[],2);
            FOV = (rot90(uint16(accumarray([i j k],CC,[options.Nx options.Ny options.Nz])),-2));
            FOV = circshift(FOV, round(size(FOV,1)/6));
            C{2} = C{2} + FOV;
            
            CC = ones(length(S),1,'int16');
            t1 = true(size(S,4),1);
            t2 = all([abs(S(:,4))<=max(x) abs(S(:,5))<=max(y) abs(S(:,6))<=max(z)],2);
            t = (t1+t2==2);
            CC = CC(t);
            [~,i] = min(abs(bsxfun(@minus,x.',S(t,4))),[],2);
            [~,j] = min(abs(bsxfun(@minus,y.',S(t,5))),[],2);
            [~,k] = min(abs(bsxfun(@minus,z.',S(t,6))),[],2);
            FOV = (rot90(uint16(accumarray([i j k],CC,[options.Nx options.Ny options.Nz])),-2));
            FOV = circshift(FOV, round(size(FOV,1)/6));
            C{3} = C{3} + FOV;
            
            CC = ones(length(S),1,'int16');
            t1 = true(size(S,1),1);
            t2 = all([abs(S(:,1))<=max(x) abs(S(:,2))<=max(y) abs(S(:,3))<=max(z)],2);
            t3 = all([abs(S(:,4))<=max(x) abs(S(:,5))<=max(y) abs(S(:,6))<=max(z)],2);
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
            FOV = (rot90(uint16(accumarray([i j k],CC,[options.Nx options.Ny options.Nz])),-2));
            FOV = circshift(FOV, round(size(FOV,1)/6));
            FOV = FOV/2;
            C{4} = C{4} + FOV;
            
            CC = ones(length(S),1,'int16');
            t1 = true(size(S,1),1);
            t2 = all([abs(S(:,1))<=max(x) abs(S(:,2))<=max(y) abs(S(:,3))<=max(z)],2);
            t3 = all([abs(S(:,4))<=max(x) abs(S(:,5))<=max(y) abs(S(:,6))<=max(z)],2);
            t = (t1+t2+t3==3);
            CC = CC(t);
            C1 = mean([S(t,1) S(t,4)],2);
            C2 = mean([S(t,2) S(t,5)],2);
            C3 = mean([S(t,3) S(t,6)],2);
            [~,i] = min(abs(bsxfun(@minus,(x).',C1)),[],2);
            [~,j] = min(abs(bsxfun(@minus,(y).',C2)),[],2);
            [~,k] = min(abs(bsxfun(@minus,(z).',C3)),[],2);
            FOV = (rot90(uint16(accumarray([i j k],CC,[options.Nx options.Ny options.Nz])),-2));
            FOV = circshift(FOV, round(size(FOV,1)/6));
            C{5} = C{5} + FOV;
        end
        
        if partitions == 1
            LL1 = LL1 + L1;
        else
            LL1 = [LL1;L1];
            LL2 = [LL2;L2];
        end
        if options.verbose
            disp(['File ' fnames(lk).name ' loaded'])
        end
        
    end
    if source
        save([machine_name '_Ideal_image_coordinates_' name '_LMF.mat'],'C','-v7.3')
    end
    
    % forms a sparse matrix containing all the events at different LORs
    % e.g. a value at [325 25100] contains the number of events detected at
    % detectors 325 (detector 1) and 25100 (detector 2)
    % a value at [25100 325] on the other hand tells the number of events when
    % 25100 is detector 1 and 325 detector 2
    if partitions == 1
        prompts = LL1;
        clear LL1 LL2
    else
        prompts = cell(partitions,1);
        ll = 1;
        for kk = 1 : partitions
            index = find(M(:,time_index)<tpoints(kk),1,'last');
            prompts{kk} = accumarray([LL1(ll:index) LL2(ll:index)],1,[detectors detectors],[],[],true);
            ll = ll + index;
        end
    end
    
    coincidences = cell(partitions,1);
    
    if options.use_raw_data
        load([machine_name '_detector_locations_' num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '_raw.mat'],'discard')
        for llo=1:partitions
            
            if partitions > 1
                P = uint16(full(prompts{llo}));
            else
                P = uint16(full(prompts));
            end
            
            %     tic
            P = tril(P,0) + triu(P,1)';
            P = (P(tril(true(size(P)), 0)));
            P = P(discard);
            
            coincidences{partitions, 1} = sparse(double(P));
            toc
            
        end
        if partitions == 1
            save([machine_name '_measurements_' name '_static_raw_LMF.mat'], 'coincidences', '-v7.3')
        else
            save([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_LMF.mat'], 'coincidences', '-v7.3')
        end
    else
        load([num2str(machine_name) '_swap_corners_' num2str(options.Ndist) '.mat']);
        for llo=1:partitions
            if partitions > 1
                prompt = prompts{llo};
            else
                prompt = prompts;
            end
            
            eka=1;
            counts = 0;
            ll = 1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            P1 = cell(options.rings + length(options.pseudot),(options.rings + length(options.pseudot) - eka + 1));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            tic
            for ring=eka:options.rings % Ring at a time
                for luuppi=1:(options.rings-eka+1)-(ring-1) % All the LORs with other options.rings as well
                    if luuppi==1 % Both of the detectors are on the same ring
                        % If pseudo ring then place the measurements to the next
                        % ring
                        if ismember(ll,options.pseudot)
                            ll = ll + 1;
                        end
                        lk = ll + 1;
                        % Observations equaling the ones detected on the current
                        % ring
                        temppiP=prompt((1+options.det_per_ring*(ring-1)):(options.det_per_ring+options.det_per_ring*(ring-1)),1+options.det_per_ring*(ring-1):(options.det_per_ring+options.det_per_ring*(ring-1)));
                        % combine same LORs, but with different detector order
                        % (i.e. combine values at [325 25100] and [25100 325])
                        temppi=(temppiP'+temppiP);
                        temppi= temppi(tril(logical(true(options.det_per_ring)),0)); % Finally take only the other side
                        P1{ll,ll} = temppi(:);
                        if options.verbose
                            counts = counts + sum(temppi(:));
                        end
                        
                        
                    else % Detectors on different rings
                        if ismember(lk,options.pseudot)
                            lk = lk + 1;
                        end
                        temppiP=(prompt((1+options.det_per_ring*(ring-1)):(options.det_per_ring+options.det_per_ring*(ring-1)),1+options.det_per_ring*(ring-1)+options.det_per_ring*(luuppi-1):(options.det_per_ring+options.det_per_ring*(ring-1)+options.det_per_ring*(luuppi-1))))';
                        temppiP2=(prompt(1+options.det_per_ring*(ring-1)+options.det_per_ring*(luuppi-1):(options.det_per_ring+options.det_per_ring*(ring-1)+options.det_per_ring*(luuppi-1)),(1+options.det_per_ring*(ring-1)):(options.det_per_ring+options.det_per_ring*(ring-1))));
                        temppiP3 = temppiP;
                        
                        % Combine LORs
                        temppiP = triu(temppiP) + triu(temppiP2);
                        temppiP2 = tril(temppiP3) + tril(temppiP2);
                        
                        temppiP3 = temppiP;
                        
                        % Swap corners
                        temppiP4 = zeros(size(temppiP2),'uint16');
                        temppiP4(apu3) = temppiP2(apu3);
                        temppiP(apu1) = 0;
                        temppiP = temppiP + temppiP4';
                        temppiP4 = zeros(size(temppiP2),'uint16');
                        temppiP4(apu4) = temppiP2(apu4);
                        temppiP(apu2) = 0;
                        temppiP = temppiP + temppiP4';
                        temppiP = temppiP';
                        
                        temppiP4 = zeros(size(temppiP3),'uint16');
                        temppiP4(apu1) = temppiP3(apu1);
                        temppiP2(apu3) = 0;
                        temppiP2 = temppiP2 + temppiP4';
                        temppiP4 = zeros(size(temppiP3),'uint16');
                        temppiP4(apu2) = temppiP3(apu2);
                        temppiP2(apu4) = 0;
                        temppiP2 = temppiP2 + temppiP4';
                        
                        % Take only the other side
                        temppiP2 = temppiP2(tril(true(options.det_per_ring),0));
                        temppiP = temppiP(tril(true(options.det_per_ring),0));
                        
                        P1{ll,lk} = temppiP(:);
                        P1{lk,ll} = temppiP2(:);
                        if options.verbose
                            counts = counts + sum(temppiP(:)) + sum(temppiP2(:));
                        end
                        lk = lk + 1;
                    end
                end
                ll = ll + 1;
            end
            toc
            clear prompt
            P1(cellfun(@isempty,P1)) = {zeros(size(temppi),'uint16')};
            coincidences{partitions} = P1;
            
            if options.verbose
                disp(['Total number of counts at time point ' num2str(llo) ' is ' num2str(counts) '.'])
            end
            
        end
        if partitions == 1
            save([machine_name '_measurements_' name '_static_LMF.mat'], 'coincidences', '-v7.3')
        else
            save([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_LMF.mat'], 'coincidences', '-v7.3')
        end
    end
elseif options.use_ASCII
    rsector_ind1 = options.rsector_ind1;
    rsector_ind2 = options.rsector_ind2;
    crs_ind1 = options.crs_ind1;
    crs_ind2 = options.crs_ind2;
    det_per_ring = options.det_per_ring;
    time_index = options.time_index;
    source_index1 = options.source_index1;
    source_index2 = options.source_index2;
    
    % Take coincidence data
    fnames = dir([fpath '*Coincidences*.dat']);
%     C = [];
    LL1 = [];
    LL2 = [];
    
    C = cell(5,1);
    C(:) = {zeros(options.Nx, options.Ny, options.Nz, 'uint16')};
    
    % Go through all the files
    for lk=1:length(fnames)
        
        M = dlmread([fpath fnames(lk).name]);
        
        
        % An ideal image can be formed with this, see the file
        % Compare_Phantom_vs_Reconstructions.m
        % Takes the columns containing the source locations for both singles
        if source_index1 ~= 0 || ~isempty(source_index1) || source_index2 ~= 0 || ~isempty(source_index2)
            
            S = [single(M(:,[source_index1:source_index1+2 source_index2:source_index2+2]))];
            
            pixel_width = options.FOVa/options.Nx;
            z_width = options.axial_fov/options.Nz;
            x = single((-pixel_width*(options.Nx-1)/2:pixel_width:pixel_width*(options.Nx-1)/2)');
            y = x;
            z = single((-z_width*(options.Nz-1)/2:z_width:z_width*(options.Nz-1)/2)');
            
            CC = ones(length(S),1,'int16');
            t1 = all(S(:,1:3)==S(:,4:6),2);
            t2 = all([abs(S(:,1))<=max(x) abs(S(:,2))<=max(y) abs(S(:,3))<=max(z)],2);
            t3 = all([abs(S(:,4))<=max(x) abs(S(:,5))<=max(y) abs(S(:,6))<=max(z)],2);
            t = (t1+t2+t3==3);
            CC = CC(t);
            [~,i] = min(abs(bsxfun(@minus,x.',S(t,1))),[],2);
            [~,j] = min(abs(bsxfun(@minus,y.',S(t,2))),[],2);
            [~,k] = min(abs(bsxfun(@minus,z.',S(t,3))),[],2);
            FOV = (rot90(uint16(accumarray([i j k],CC,[options.Nx options.Ny options.Nz])),-2));
            FOV = circshift(FOV, round(size(FOV,1)/6));
            C{1} = C{1} + FOV;
            
            CC = ones(length(S),1,'int16');
            t1 = true(size(S,1),1);
            t2 = all([abs(S(:,1))<=max(x) abs(S(:,2))<=max(y) abs(S(:,3))<=max(z)],2);
            t = (t1+t2==2);
            CC = CC(t);
            [~,i] = min(abs(bsxfun(@minus,x.',S(t,1))),[],2);
            [~,j] = min(abs(bsxfun(@minus,y.',S(t,2))),[],2);
            [~,k] = min(abs(bsxfun(@minus,z.',S(t,3))),[],2);
            FOV = (rot90(uint16(accumarray([i j k],CC,[options.Nx options.Ny options.Nz])),-2));
            FOV = circshift(FOV, round(size(FOV,1)/6));
            C{2} = C{2} + FOV;
            
            CC = ones(length(S),1,'int16');
            t1 = true(size(S,4),1);
            t2 = all([abs(S(:,4))<=max(x) abs(S(:,5))<=max(y) abs(S(:,6))<=max(z)],2);
            t = (t1+t2==2);
            CC = CC(t);
            [~,i] = min(abs(bsxfun(@minus,x.',S(t,4))),[],2);
            [~,j] = min(abs(bsxfun(@minus,y.',S(t,5))),[],2);
            [~,k] = min(abs(bsxfun(@minus,z.',S(t,6))),[],2);
            FOV = (rot90(uint16(accumarray([i j k],CC,[options.Nx options.Ny options.Nz])),-2));
            FOV = circshift(FOV, round(size(FOV,1)/6));
            C{3} = C{3} + FOV;
            
            CC = ones(length(S),1,'int16');
            t1 = true(size(S,1),1);
            t2 = all([abs(S(:,1))<=max(x) abs(S(:,2))<=max(y) abs(S(:,3))<=max(z)],2);
            t3 = all([abs(S(:,4))<=max(x) abs(S(:,5))<=max(y) abs(S(:,6))<=max(z)],2);
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
            CC = repmat(CC,2,1);
            FOV = (rot90(uint16(accumarray([i j k],CC,[options.Nx options.Ny options.Nz])),-2));
            FOV = circshift(FOV, round(size(FOV,1)/6));
            FOV = FOV/2;
            C{4} = C{4} + FOV;
            
            CC = ones(length(S),1,'int16');
            t1 = true(size(S,1),1);
            t2 = all([abs(S(:,1))<=max(x) abs(S(:,2))<=max(y) abs(S(:,3))<=max(z)],2);
            t3 = all([abs(S(:,4))<=max(x) abs(S(:,5))<=max(y) abs(S(:,6))<=max(z)],2);
            t = (t1+t2+t3==3);
            CC = CC(t);
            C1 = mean([S(t,1) S(t,4)],2);
            C2 = mean([S(t,2) S(t,5)],2);
            C3 = mean([S(t,3) S(t,6)],2);
            [~,i] = min(abs(bsxfun(@minus,(x).',C1)),[],2);
            [~,j] = min(abs(bsxfun(@minus,(y).',C2)),[],2);
            [~,k] = min(abs(bsxfun(@minus,(z).',C3)),[],2);
            FOV = (rot90(uint16(accumarray([i j k],CC,[options.Nx options.Ny options.Nz])),-2));
            FOV = circshift(FOV, round(size(FOV,1)/6));
            C{5} = C{5} + FOV;
        end
        
        
        interval = tot_time/partitions; % Time per time point
        timepoints = 1/interval:1/interval:tot_time; % Maximum times of each time point
        
        % The number of rings of each coincidence event (e.g. the first single hits
        % a detector on ring 5 and second a detector on ring 10)
        % [0 rings - 1]
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
        LL1 = [LL1;ring_number1*uint16(det_per_ring)+ring_pos1+1];
        LL2 = [LL2;ring_number2*uint16(det_per_ring)+ring_pos2+1];
        
        if options.verbose
            disp(['File ' fnames(lk).name ' loaded'])
        end
        
    end
    if source_index1 ~= 0 || ~isempty(source_index1) || source_index2 ~= 0 || ~isempty(source_index2)
        save([machine_name '_Ideal_image_coordinates_' name '_ASCII.mat'],'C','-v7.3')
    end
    
    % forms a sparse matrix containing all the events at different LORs
    % e.g. a value at [325 25100] contains the number of events detected at
    % detectors 325 (detector 1) and 25100 (detector 2)
    % a value at [25100 325] on the other hand tells the number of events when
    % 25100 is detector 1 and 325 detector 2
    if partitions == 1
        prompts = accumarray([LL1 LL2],uint16(1),[detectors detectors],@(x) sum(x,'native'));
    else
        prompts = cell(partitions,1);
        ll = 1;
        for kk = 1 : partitions
            index = find(M(:,time_index)<timepoints(kk),1,'last');
            prompts{kk} = accumarray([LL1(ll:index) LL2(ll:index)],1,[detectors detectors],[],[],true);
            ll = ll + index;
        end
    end
    
    coincidences = cell(partitions,1);
    if options.use_raw_data
        load([machine_name '_detector_locations_' num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '_raw.mat'],'discard')
        
        for llo=1:partitions
            
            if partitions > 1
                P = uint16(full(prompts{llo}));
            else
                P = uint16(full(prompts));
            end
            
            P = tril(P,0) + triu(P,1)';
            P = (P(tril(true(size(P)), 0)));
            P = P(discard);
            
            coincidences{partitions, 1} = sparse(double(P));
            
        end
        if partitions == 1
            save([machine_name '_measurements_' name '_static_raw_ASCII.mat'], 'coincidences', '-v7.3')
        else
            save([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_ASCII.mat'], 'coincidences', '-v7.3')
        end
    else
        % Forms a raw Michelogram that is still in the detector space
        load([num2str(machine_name) '_swap_corners_' num2str(options.Ndist) '.mat']);
        for llo=1:partitions
            %     tic
            
            if partitions > 1
                prompt = prompts{llo};
            else
                prompt = prompts;
            end
            
            eka=1;
            counts = 0;
            ll = 1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            P1 = cell(options.rings + length(options.pseudot),(options.rings + length(options.pseudot) - eka + 1));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            for ring=eka:options.rings % Ring at a time
                for luuppi=1:(options.rings-eka+1)-(ring-1) % All the LORs with other options.rings as well
                    if luuppi==1 % Both of the detectors are on the same ring
                        % If pseudo ring then place the measurements to the next
                        % ring
                        if ismember(ll,options.pseudot)
                            ll = ll + 1;
                        end
                        lk = ll + 1;
                        % Observations equaling the ones detected on the current
                        % ring
                        temppiP=prompt((1+options.det_per_ring*(ring-1)):(options.det_per_ring+options.det_per_ring*(ring-1)),1+options.det_per_ring*(ring-1):(options.det_per_ring+options.det_per_ring*(ring-1)));
                        % combine same LORs, but with different detector order
                        % (i.e. combine values at [325 25100] and [25100 325])
                        temppi=(temppiP'+temppiP);
                        temppi= temppi(tril(logical(true(options.det_per_ring)),0)); % Finally take only the other side
                        P1{ll,ll} = temppi(:);
                        if options.verbose
                            counts = counts + sum(temppi(:));
                        end
                        
                        
                    else % Detectors on different rings
                        if ismember(lk,options.pseudot)
                            lk = lk + 1;
                        end
                        temppiP=(prompt((1+options.det_per_ring*(ring-1)):(options.det_per_ring+options.det_per_ring*(ring-1)),1+options.det_per_ring*(ring-1)+options.det_per_ring*(luuppi-1):(options.det_per_ring+options.det_per_ring*(ring-1)+options.det_per_ring*(luuppi-1))))';
                        temppiP2=(prompt(1+options.det_per_ring*(ring-1)+options.det_per_ring*(luuppi-1):(options.det_per_ring+options.det_per_ring*(ring-1)+options.det_per_ring*(luuppi-1)),(1+options.det_per_ring*(ring-1)):(options.det_per_ring+options.det_per_ring*(ring-1))));
                        temppiP3 = temppiP;
                        
                        % Combine LORs
                        temppiP = triu(temppiP) + triu(temppiP2);
                        temppiP2 = tril(temppiP3) + tril(temppiP2);
                        
                        temppiP3 = temppiP;
                        
                        % Swap corners
                        temppiP4 = zeros(size(temppiP2),'uint16');
                        temppiP4(apu3) = temppiP2(apu3);
                        temppiP(apu1) = 0;
                        temppiP = temppiP + temppiP4';
                        temppiP4 = zeros(size(temppiP2),'uint16');
                        temppiP4(apu4) = temppiP2(apu4);
                        temppiP(apu2) = 0;
                        temppiP = temppiP + temppiP4';
                        temppiP = temppiP';
                        
                        temppiP4 = zeros(size(temppiP3),'uint16');
                        temppiP4(apu1) = temppiP3(apu1);
                        temppiP2(apu3) = 0;
                        temppiP2 = temppiP2 + temppiP4';
                        temppiP4 = zeros(size(temppiP3),'uint16');
                        temppiP4(apu2) = temppiP3(apu2);
                        temppiP2(apu4) = 0;
                        temppiP2 = temppiP2 + temppiP4';
                        
                        % Take only the other side
                        temppiP2 = temppiP2(tril(true(options.det_per_ring),0));
                        temppiP = temppiP(tril(true(options.det_per_ring),0));
                        
                        P1{ll,lk} = temppiP(:);
                        P1{lk,ll} = temppiP2(:);
                        if options.verbose
                            counts = counts + sum(temppiP(:)) + sum(temppiP2(:));
                        end
                        lk = lk + 1;
                    end
                end
                ll = ll + 1;
            end
            clear prompt
            P1(cellfun(@isempty,P1)) = {zeros(size(temppi),'uint16')};
            coincidences{partitions} = P1;
            
            if options.verbose
                disp(['Total number of counts at time point ' num2str(llo) ' is ' num2str(counts) '.'])
            end
            
        end
        if partitions == 1
            save([machine_name '_measurements_' name '_static_ASCII.mat'], 'coincidences', '-v7.3')
        else
            save([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_ASCII.mat'], 'coincidences', '-v7.3')
        end
    end
elseif options.use_root
    blocks_per_ring = uint32(blocks_per_ring);
    cryst_per_block = uint32(cryst_per_block);
    det_per_ring = uint32(options.det_per_ring);
    source = options.source;
    
    
    % Take coincidence data
    fnames = dir([fpath '*.root']);
    
    
    if partitions == 1
        LL1 = zeros(detectors, detectors, 'uint16');
    else
        LL1 = [];
        LL2 = [];
    end
    
    C = cell(5,1);
    C(:) = {zeros(options.Nx, options.Ny, options.Nz, 'uint16')};
    
    % Go through all the files
    for lk=1:length(fnames)
        
        
        nimi = [fpath fnames(lk).name];
        
        if source
            [L1, L2, tpoints, S] = GATE_root_matlab(nimi,vali,alku,loppu, uint32(detectors), blocks_per_ring, cryst_per_block, ...
                det_per_ring, uint32(options.linear_multip), source);
        else
            [L1, L2, tpoints, output] = GATE_root_matlab(nimi,vali,alku,loppu, uint32(detectors), blocks_per_ring, cryst_per_block, ...
                det_per_ring, uint32(options.linear_multip), source);
        end
        clear mex
        
        
        % An ideal image can be formed with this. Takes the variables
        % containing the source locations for both singles 
        if source
            pixel_width = options.FOVa/options.Nx;
            z_width = options.axial_fov/options.Nz;
            x = single((-pixel_width*(options.Nx-1)/2:pixel_width:pixel_width*(options.Nx-1)/2)');
            y = x;
            z = single((-z_width*(options.Nz-1)/2:z_width:z_width*(options.Nz-1)/2)');
            
            CC = ones(length(S),1,'int16');
            t1 = all(S(:,1:3)==S(:,4:6),2);
            t2 = all([abs(S(:,1))<=max(x) abs(S(:,2))<=max(y) abs(S(:,3))<=max(z)],2);
            t3 = all([abs(S(:,4))<=max(x) abs(S(:,5))<=max(y) abs(S(:,6))<=max(z)],2);
            t = (t1+t2+t3==3);
            CC = CC(t);
            [~,i] = min(abs(bsxfun(@minus,x.',S(t,1))),[],2);
            [~,j] = min(abs(bsxfun(@minus,y.',S(t,2))),[],2);
            [~,k] = min(abs(bsxfun(@minus,z.',S(t,3))),[],2);
            FOV = (rot90(uint16(accumarray([i j k],CC,[options.Nx options.Ny options.Nz])),0));
%             FOV = circshift(FOV, round(size(FOV,1)/6));
            C{1} = C{1} + FOV;
            
            CC = ones(length(S),1,'int16');
            t1 = true(size(S,1),1);
            t2 = all([abs(S(:,1))<=max(x) abs(S(:,2))<=max(y) abs(S(:,3))<=max(z)],2);
            t = (t1+t2==2);
            CC = CC(t);
            [~,i] = min(abs(bsxfun(@minus,x.',S(t,1))),[],2);
            [~,j] = min(abs(bsxfun(@minus,y.',S(t,2))),[],2);
            [~,k] = min(abs(bsxfun(@minus,z.',S(t,3))),[],2);
            FOV = (rot90(uint16(accumarray([i j k],CC,[options.Nx options.Ny options.Nz])),0));
%             FOV = circshift(FOV, round(size(FOV,1)/6));
            C{2} = C{2} + FOV;
            
            CC = ones(length(S),1,'int16');
            t1 = true(size(S,4),1);
            t2 = all([abs(S(:,4))<=max(x) abs(S(:,5))<=max(y) abs(S(:,6))<=max(z)],2);
            t = (t1+t2==2);
            CC = CC(t);
            [~,i] = min(abs(bsxfun(@minus,x.',S(t,4))),[],2);
            [~,j] = min(abs(bsxfun(@minus,y.',S(t,5))),[],2);
            [~,k] = min(abs(bsxfun(@minus,z.',S(t,6))),[],2);
            FOV = (rot90(uint16(accumarray([i j k],CC,[options.Nx options.Ny options.Nz])),0));
%             FOV = circshift(FOV, round(size(FOV,1)/6));
            C{3} = C{3} + FOV;
            
            CC = ones(length(S),1,'int16');
            t1 = true(size(S,1),1);
            t2 = all([abs(S(:,1))<=max(x) abs(S(:,2))<=max(y) abs(S(:,3))<=max(z)],2);
            t3 = all([abs(S(:,4))<=max(x) abs(S(:,5))<=max(y) abs(S(:,6))<=max(z)],2);
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
            CC = repmat(CC,2,1);
            FOV = (rot90(uint16(accumarray([i j k],CC,[options.Nx options.Ny options.Nz])),0));
%             FOV = circshift(FOV, round(size(FOV,1)/6));
            FOV = FOV/2;
            C{4} = C{4} + FOV;
            
            CC = ones(length(S),1,'int16');
            t1 = true(size(S,1),1);
            t2 = all([abs(S(:,1))<=max(x) abs(S(:,2))<=max(y) abs(S(:,3))<=max(z)],2);
            t3 = all([abs(S(:,4))<=max(x) abs(S(:,5))<=max(y) abs(S(:,6))<=max(z)],2);
            t = (t1+t2+t3==3);
            CC = CC(t);
            C1 = mean([S(t,1) S(t,4)],2);
            C2 = mean([S(t,2) S(t,5)],2);
            C3 = mean([S(t,3) S(t,6)],2);
            [~,i] = min(abs(bsxfun(@minus,(x).',C1)),[],2);
            [~,j] = min(abs(bsxfun(@minus,(y).',C2)),[],2);
            [~,k] = min(abs(bsxfun(@minus,(z).',C3)),[],2);
            FOV = (rot90(uint16(accumarray([i j k],CC,[options.Nx options.Ny options.Nz])),0));
%             FOV = circshift(FOV, round(size(FOV,1)/6));
            C{5} = C{5} + FOV;
        end
        
        
        if partitions == 1
            LL1 = LL1 + L1;
        else
            LL1 = [LL1;L1];
            LL2 = [LL2;L2];
        end
%         if options.verbose
%             disp(['File ' fnames(lk).name ' loaded'])
%         end
        
    end
    if source
        save([machine_name '_Ideal_image_coordinates_' name '_Root.mat'],'C','-v7.3')
    end
    
    % forms a sparse matrix containing all the events at different LORs
    % e.g. a value at [325 25100] contains the number of events detected at
    % detectors 325 (detector 1) and 25100 (detector 2)
    % a value at [25100 325] on the other hand tells the number of events when
    % 25100 is detector 1 and 325 detector 2
    if partitions == 1
        prompts = LL1;
        clear LL1 LL2
    else
        prompts = cell(partitions,1);
        ll = 1;
        for kk = 1 : partitions
            index = find(M(:,time_index)<tpoints(kk),1,'last');
            prompts{kk} = accumarray([LL1(ll:index) LL2(ll:index)],1,[detectors detectors],[],[],true);
            ll = ll + index;
        end
    end
    
    coincidences = cell(partitions,1);
    
    if options.use_raw_data
        load([machine_name '_detector_locations_' num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '_raw.mat'],'discard')
        for llo=1:partitions
            
            if partitions > 1
                P = uint16(full(prompts{llo}));
            else
                P = uint16(full(prompts));
            end
            
            %     tic
            P = tril(P,0) + triu(P,1)';
            P = (P(tril(true(size(P)), 0)));
            P = P(discard);
            
            coincidences{partitions, 1} = sparse(double(P));
            toc
            
        end
        if partitions == 1
            save([machine_name '_measurements_' name '_static_raw_root.mat'], 'coincidences', '-v7.3')
        else
            save([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_root.mat'], 'coincidences', '-v7.3')
        end
    else
        load([num2str(machine_name) '_swap_corners_' num2str(options.Ndist) '.mat']);
        for llo=1:partitions
            if partitions > 1
                prompt = prompts{llo};
            else
                prompt = prompts;
            end
            
            eka=1;
            counts = 0;
            ll = 1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            P1 = cell(options.rings + length(options.pseudot),(options.rings + length(options.pseudot) - eka + 1));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            tic
            for ring=eka:options.rings % Ring at a time
                for luuppi=1:(options.rings-eka+1)-(ring-1) % All the LORs with other options.rings as well
                    if luuppi==1 % Both of the detectors are on the same ring
                        % If pseudo ring then place the measurements to the next
                        % ring
                        if ismember(ll,options.pseudot)
                            ll = ll + 1;
                        end
                        lk = ll + 1;
                        % Observations equaling the ones detected on the current
                        % ring
                        temppiP=prompt((1+options.det_per_ring*(ring-1)):(options.det_per_ring+options.det_per_ring*(ring-1)),1+options.det_per_ring*(ring-1):(options.det_per_ring+options.det_per_ring*(ring-1)));
                        % combine same LORs, but with different detector order
                        % (i.e. combine values at [325 25100] and [25100 325])
                        temppi=(temppiP'+temppiP);
                        temppi= temppi(tril(logical(true(options.det_per_ring)),0)); % Finally take only the other side
                        P1{ll,ll} = temppi(:);
                        if options.verbose
                            counts = counts + sum(temppi(:));
                        end
                        
                        
                    else % Detectors on different rings
                        if ismember(lk,options.pseudot)
                            lk = lk + 1;
                        end
                        temppiP=(prompt((1+options.det_per_ring*(ring-1)):(options.det_per_ring+options.det_per_ring*(ring-1)),1+options.det_per_ring*(ring-1)+options.det_per_ring*(luuppi-1):(options.det_per_ring+options.det_per_ring*(ring-1)+options.det_per_ring*(luuppi-1))))';
                        temppiP2=(prompt(1+options.det_per_ring*(ring-1)+options.det_per_ring*(luuppi-1):(options.det_per_ring+options.det_per_ring*(ring-1)+options.det_per_ring*(luuppi-1)),(1+options.det_per_ring*(ring-1)):(options.det_per_ring+options.det_per_ring*(ring-1))));
                        temppiP3 = temppiP;
                        
                        % Combine LORs
                        temppiP = triu(temppiP) + triu(temppiP2);
                        temppiP2 = tril(temppiP3) + tril(temppiP2);
                        
                        temppiP3 = temppiP;
                        
                        % Swap corners
                        temppiP4 = zeros(size(temppiP2),'uint16');
                        temppiP4(apu3) = temppiP2(apu3);
                        temppiP(apu1) = 0;
                        temppiP = temppiP + temppiP4';
                        temppiP4 = zeros(size(temppiP2),'uint16');
                        temppiP4(apu4) = temppiP2(apu4);
                        temppiP(apu2) = 0;
                        temppiP = temppiP + temppiP4';
                        temppiP = temppiP';
                        
                        temppiP4 = zeros(size(temppiP3),'uint16');
                        temppiP4(apu1) = temppiP3(apu1);
                        temppiP2(apu3) = 0;
                        temppiP2 = temppiP2 + temppiP4';
                        temppiP4 = zeros(size(temppiP3),'uint16');
                        temppiP4(apu2) = temppiP3(apu2);
                        temppiP2(apu4) = 0;
                        temppiP2 = temppiP2 + temppiP4';
                        
                        % Take only the other side
                        temppiP2 = temppiP2(tril(true(options.det_per_ring),0));
                        temppiP = temppiP(tril(true(options.det_per_ring),0));
                        
                        P1{ll,lk} = temppiP(:);
                        P1{lk,ll} = temppiP2(:);
                        if options.verbose
                            counts = counts + sum(temppiP(:)) + sum(temppiP2(:));
                        end
                        lk = lk + 1;
                    end
                end
                ll = ll + 1;
            end
            toc
            clear prompt
            P1(cellfun(@isempty,P1)) = {zeros(size(temppi),'uint16')};
            coincidences{partitions} = P1;
            
%             if options.verbose
%                 disp(['Total number of counts at time point ' num2str(llo) ' is ' num2str(counts) '.'])
%             end
            
        end
        if partitions == 1
            save([machine_name '_measurements_' name '_static_root.mat'], 'coincidences', '-v7.3')
        else
            save([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_root.mat'], 'coincidences', '-v7.3')
        end
    end
end

if options.verbose
    disp('Measurements loaded and saved')
    toc
end