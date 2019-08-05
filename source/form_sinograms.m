function [varargout] = form_sinograms(options)
%% FORM MICHELOGRAMS AND SINOGRAMS FROM RAW DETECTOR PAIR DATA
% This code forms the sinograms for the current machine from the raw
% list-mode data using the raw detector pairs. It first creates the
% Michelogram and later the sinograms with the dimensions specified by the
% user. 
% 
% Input the machine and sinogram properties and the raw detector pair
% data. 
% Output is the sinograms for each time point, saved in a mat-file in the
% current working directory.
% 
% Randoms, scatter and normalization corrections are applied if applicable.
%
% OUTPUT: 
%   raw_SinM = The raw sinogram without any corrections or, if Inveon data
%   is used, the Inveon sinogram file. In the former case this is the first
%   output ONLY if no corrections are applied OR if corrections are applied
%   during reconstruction
%   SinM = Sinogram with all the selected corrections applied. First input
%   if any corrections are applied and corrections are NOT performed during
%   image reconstruction, otherwise the second input
%   SinD = Delayed coincidence sinogram
%   SinTrues = True coincidences (GATE only)
%   SinScatter = Scattered coincidence sinogram (GATE only)
%   SinRandoms = Randoms coincidence sinogram (GATE only)
%
% See also michelogram_formation, load_data, sinogram_coordinates_2D,
% sinogram_coordinates_3D



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi, Samuli Summala
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
    error('Too many output arguments')
end
if options.use_machine < 2
    
    machine_name = options.machine_name;
    name = options.name;
    rings = options.rings;
    Nang = options.Nang;
    Ndist = options.Ndist;
    pseudot = options.pseudot;
    segment_table = options.segment_table;
    span = options.span;
    NSlices = options.TotSinos;
    Nz = options.Nz;
    ring_difference = options.ring_difference;
    partitions = options.partitions;
    tot_time = options.tot_time;
    
    ringsp = rings + sum(pseudot);
    
    variableList = {'SinM','SinTrues','SinScatter','SinRandoms','SinDelayed','SinSC'};
    
    folder = fileparts(which('form_sinograms.m'));
    folder = strrep(folder, 'source','mat-files/');
    folder = strrep(folder, '\','/');
    
    if options.verbose
        disp('Forming initial Michelogram')
        tic
    end
    if options.use_machine == 1
        if exist( 'options.coincidences' , 'var') == 0 || options.precompute_all
            if options.partitions == 1
                load([machine_name '_measurements_' name '_static_raw_listmode.mat'], 'coincidences')
            else
                load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_listmode.mat'], 'coincidences')
            end
            options.coincidences = coincidences;
            clear coincidences
        end
        if options.randoms_correction
            if options.partitions == 1
                load([machine_name '_measurements_' name '_static_raw_listmode.mat'], 'delayed_coincidences')
            else
                load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_listmode.mat'], 'delayed_coincidences')
            end
            delayed_coincidences = initial_michelogram(options, delayed_coincidences);
        end
    else
        if exist( 'options.coincidences' , 'var') == 0 || options.precompute_all
            if options.partitions == 1
                if options.use_ASCII
                    load([machine_name '_measurements_' name '_static_raw_ASCII.mat'], 'coincidences')
                elseif options.use_LMF
                    load([machine_name '_measurements_' name '_static_raw_LMF.mat'], 'coincidences')
                elseif options.use_root
                    load([machine_name '_measurements_' name '_static_raw_root.mat'], 'coincidences')
                else
                    load([machine_name '_measurements_' name '_static_raw_real.mat'], 'coincidences')
                end
            else
                if options.use_ASCII
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_ASCII.mat'], 'coincidences')
                elseif options.use_LMF
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_LMF.mat'], 'coincidences')
                elseif options.use_root
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_root.mat'], 'coincidences')
                else
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_real.mat'], 'coincidences')
                end
            end
            options.coincidences = coincidences;
            clear coincidences
        end
        if options.obtain_trues
            if options.partitions == 1
                if options.use_ASCII
                    load([machine_name '_measurements_' name '_static_raw_ASCII.mat'], 'true_coincidences')
                elseif options.use_LMF
                    load([machine_name '_measurements_' name '_static_raw_LMF.mat'], 'true_coincidences')
                elseif options.use_root
                    load([machine_name '_measurements_' name '_static_raw_root.mat'], 'true_coincidences')
                else
                    load([machine_name '_measurements_' name '_static_raw_real.mat'], 'true_coincidences')
                end
            else
                if options.use_ASCII
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_ASCII.mat'], 'true_coincidences')
                elseif options.use_LMF
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_LMF.mat'], 'true_coincidences')
                elseif options.use_root
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_root.mat'], 'true_coincidences')
                else
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_real.mat'], 'true_coincidences')
                end
            end
            true_coincidences = initial_michelogram(options, true_coincidences);
            options.true_coincidences = true_coincidences;
            clear true_coincidences
        end
        if options.store_scatter
            if options.partitions == 1
                if options.use_ASCII
                    load([machine_name '_measurements_' name '_static_raw_ASCII.mat'], 'scattered_coincidences')
                elseif options.use_LMF
                    load([machine_name '_measurements_' name '_static_raw_LMF.mat'], 'scattered_coincidences')
                elseif options.use_root
                    load([machine_name '_measurements_' name '_static_raw_root.mat'], 'scattered_coincidences')
                else
                    load([machine_name '_measurements_' name '_static_raw_real.mat'], 'scattered_coincidences')
                end
            else
                if options.use_ASCII
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_ASCII.mat'], 'scattered_coincidences')
                elseif options.use_LMF
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_LMF.mat'], 'scattered_coincidences')
                elseif options.use_root
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_root.mat'], 'scattered_coincidences')
                else
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_real.mat'], 'scattered_coincidences')
                end
            end
            scattered_coincidences = initial_michelogram(options, scattered_coincidences);
            options.scattered_coincidences = scattered_coincidences;
            clear scattered_coincidences
        end
        if options.store_randoms
            if options.partitions == 1
                if options.use_ASCII
                    load([machine_name '_measurements_' name '_static_raw_ASCII.mat'], 'random_coincidences')
                elseif options.use_LMF
                    load([machine_name '_measurements_' name '_static_raw_LMF.mat'], 'random_coincidences')
                elseif options.use_root
                    load([machine_name '_measurements_' name '_static_raw_root.mat'], 'random_coincidences')
                else
                    load([machine_name '_measurements_' name '_static_raw_real.mat'], 'random_coincidences')
                end
            else
                if options.use_ASCII
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_ASCII.mat'], 'random_coincidences')
                elseif options.use_LMF
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_LMF.mat'], 'random_coincidences')
                elseif options.use_root
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_root.mat'], 'random_coincidences')
                else
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_real.mat'], 'random_coincidences')
                end
            end
            random_coincidences = initial_michelogram(options, random_coincidences);
            options.random_coincidences = random_coincidences;
            clear random_coincidences
        end
        
        if options.randoms_correction && (options.use_ASCII || options.use_LMF || options.use_root)
            if options.partitions == 1
                if options.use_ASCII
                    load([machine_name '_measurements_' name '_static_raw_ASCII.mat'], 'delayed_coincidences')
                elseif options.use_LMF
                    load([machine_name '_measurements_' name '_static_raw_LMF.mat'], 'delayed_coincidences')
                elseif options.use_root
                    load([machine_name '_measurements_' name '_static_raw_root.mat'], 'delayed_coincidences')
                else
                    load([machine_name '_measurements_' name '_static_raw_real.mat'], 'delayed_coincidences')
                end
            else
                if options.use_ASCII
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_ASCII.mat'], 'delayed_coincidences')
                elseif options.use_LMF
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_LMF.mat'], 'delayed_coincidences')
                elseif options.use_root
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_root.mat'], 'delayed_coincidences')
                else
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_real.mat'], 'delayed_coincidences')
                end
            end
            delayed_coincidences = initial_michelogram(options, delayed_coincidences);
        end
    end
    options.coincidences = initial_michelogram(options, options.coincidences);
    
    if options.verbose
        disp('Initial Michelogram formed')
        toc
        tic
    end
    
if exist([folder machine_name '_app_coordinates_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'], 'file') == 2
    load([folder machine_name '_app_coordinates_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'i', 'j', 'accepted_lors');
else
    sinogram_coordinates_2D(options);
    load([folder machine_name '_app_coordinates_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'i', 'j', 'accepted_lors');
end
    
    
    if partitions > 1
        if options.randoms_correction || options.scatter_correction || options.normalization_correction
            SinM = cell(partitions,1);
        end
        if options.obtain_trues && options.use_machine == 0
            SinTrues = cell(partitions,1);
        end
        if options.store_scatter && options.use_machine == 0
            SinScatter = cell(partitions,1);
        end
        if options.store_randoms && options.use_machine == 0
            SinRandoms = cell(partitions,1);
        end
        if options.randoms_correction
            SinDelayed = cell(partitions,1);
        end
        raw_SinM = cell(partitions,1);
    end
    
    for llo = 1 : partitions
        
        
        if options.verbose
            disp(['Sinogram formation at timestep ' num2str(llo) ' began'])
        end
        
        P1 = options.coincidences{llo};
        
        Sinog = cell(ringsp,ringsp);
        
        ix=cellfun('isempty',Sinog);
        Sinog(ix)={zeros(Ndist,Nang,'uint16')};
        
        if options.obtain_trues && options.use_machine == 0
            T1 = options.true_coincidences{llo};
            SinogT = cell(ringsp,ringsp);
            SinogT(ix)={zeros(Ndist,Nang,'uint16')};
        end
        if options.store_scatter && options.use_machine == 0
            S1 = options.scattered_coincidences{llo};
            SinogS = cell(ringsp,ringsp);
            SinogS(ix)={zeros(Ndist,Nang,'uint16')};
        end
        if options.store_randoms && options.use_machine == 0
            R1 = options.random_coincidences{llo};
            SinogR = cell(ringsp,ringsp);
            SinogR(ix)={zeros(Ndist,Nang,'uint16')};
        end
        if options.randoms_correction && (options.use_ASCII || options.use_LMF || options.use_root || options.use_machine == 1)
            D1 = delayed_coincidences{llo};
            SinogD = cell(ringsp,ringsp);
            SinogD(ix)={zeros(Ndist,Nang,'uint16')};
        end
        
        % Create the Michelograms
        tic
        for ii=1:ringsp
            if any(ii==pseudot)
                continue
            end
            for jj=1:ringsp
                if any(jj==pseudot)
                    continue
                end
                if issparse(P1{ii,jj})
                    CC = uint16(full(P1{ii,jj}));
                else
                    CC = P1{ii,jj};
                end
                CC = CC(accepted_lors);
%                 Sinog{ii,jj}(ind) = CC;
                Sinog{ii,jj} = uint16(accumarray([i j],CC,[Ndist Nang]));
                if options.obtain_trues && options.use_machine == 0
                    if issparse(T1{ii,jj})
                        CC = uint16(full(T1{ii,jj}));
                    else
                        CC = T1{ii,jj};
                    end
                    CC = CC(accepted_lors);
                    SinogT{ii,jj} = uint16(accumarray([i j],CC,[Ndist Nang]));
                end
                if options.store_scatter && options.use_machine == 0
                    if issparse(S1{ii,jj})
                        CC = uint16(full(S1{ii,jj}));
                    else
                        CC = S1{ii,jj};
                    end
                    CC = CC(accepted_lors);
                    SinogS{ii,jj} = uint16(accumarray([i j],CC,[Ndist Nang]));
                end
                if options.store_randoms && options.use_machine == 0
                    if issparse(R1{ii,jj})
                        CC = uint16(full(R1{ii,jj}));
                    else
                        CC = R1{ii,jj};
                    end
                    CC = CC(accepted_lors);
                    SinogR{ii,jj} = uint16(accumarray([i j],CC,[Ndist Nang]));
                end
                if options.randoms_correction && (options.use_ASCII || options.use_LMF || options.use_root || options.use_machine == 1)
                    if issparse(D1{ii,jj})
                        CC = uint16(full(D1{ii,jj}));
                    else
                        CC = D1{ii,jj};
                    end
                    CC = CC(accepted_lors);
                    SinogD{ii,jj} = uint16(accumarray([i j],CC,[Ndist Nang]));
                end
            end
        end
        %     toc
        
        %%
        Sinog = cat(3,Sinog{:});
        Sin = zeros(Ndist,Nang,NSlices,'uint16');
        
        
        kkj = zeros(floor((ring_difference-ceil(span/2))/span) + 1, 1);
        for kk = 1 : floor((ring_difference-ceil(span/2))/span) + 1
            kkj(kk) = ceil(span/2) + span*(kk - 1);
        end
        offset2 = cumsum(segment_table);
        % Create the sinograms
        % First the detectors on the same ring
        Sin(:,:,1:2:Nz) = Sinog(:,:,1:ringsp+1:ringsp^2);
        % Then the detectors on adjacent rings
        for jh=1:floor(span/2)
            apu = Sinog(:,:,jh*ringsp+1:ringsp+1:ringsp^2);
            apu2 = Sinog(:,:,jh+1:ringsp+1:(ringsp-jh)*ringsp);
            Sin(:,:,jh+1:2:offset2(1)-jh) = Sin(:,:,jh+1:2:offset2(1)-jh) + apu + apu2;
        end
        % Lastly the rest of the detectors with the amount of combined LORs
        % specified with the span value
        for ih=1:floor(length(segment_table)/2)
            for jh=1:span
                apu = Sinog(:,:,(kkj(ih)+jh-1)*ringsp+1:ringsp+1:end);
                Sin(:,:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) = Sin(:,:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) + (apu);
                apu2 = Sinog(:,:,kkj(ih)+jh:ringsp+1:(ringsp-kkj(ih)-jh+1)*ringsp);
                Sin(:,:,offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1) = Sin(:,:,offset2(2*ih)+jh:2:offset2(2*ih + 1)-jh+1) + (apu2);
            end
        end
        % Fill sinogram gaps if pseudo detectors are present and the user
        % has seleced gap filling
        if options.fill_sinogram_gaps && sum(pseudot) > 0
            load([folder machine_name '_app_coordinates_' num2str(options.Ndist) 'x' num2str(options.Nang) '.mat'],'gaps')
            if strcmp('fillmissing',options.gap_filling_method)
                for kk = 1 : size(Sin,3)
                    apu = Sin(:,:,kk);
                    apu(gaps) = NaN;
                    Sin(:,:,kk) = fillmissing(apu, options.interpolation_method_fillmissing);
                end
            elseif strcmp('inpaint_nans',options.gap_filling_method)
                for kk = 1 : size(Sin,3)
                    apu = Sin(:,:,kk);
                    apu(gaps) = NaN;
                    Sin(:,:,kk) = inpaint_nans(apu, options.interpolation_method_inpaint);
                end
            else
                warning('Unsupported gap filling method! No gap filling was performed!')
            end
        end
        if partitions > 1
            raw_SinM{llo} = Sin;
        else
            raw_SinM = Sin;
        end
        % Apply normalization correction
        if options.normalization_correction
            if options.use_user_normalization
                [options.file, options.fpath] = uigetfile({'*.nrm;*.mat'},'Select normalization datafile');
                if isequal(options.file, 0)
                    error('No file was selected')
                end
                if any(strfind(options.file, '.nrm'))
                    fid = fopen(options.file);
                    normalization = fread(fid, inf, 'single=>single',0,'l');
                    fclose(fid);
                    normalization = reshape(normalization, options.Ndist, options.Nang, options.TotSinos);
                else
                    data = load(options.file);
                    variables = fields(data);
                    normalization = data.(variables{1});
                    clear data
                    if numel(normalization) ~= options.Ndist*options.Nang*options.TotSinos
                        error('Size mismatch between the current data and the normalization data file')
                    end
                    if size(normalization,3) < options.TotSinos
                        normalization = reshape(normalization, options.Ndist, options.Nang, options.TotSinos);
                    end
                end
            else
                
            end
            Sin = Sin .* normalization;
        end
        % Apply randoms correction
        if options.randoms_correction && (options.use_ASCII || options.use_LMF || options.use_root || options.use_machine == 1)
            SinogD = cat(3,SinogD{:});
            SinD = zeros(Ndist,Nang,NSlices,'uint16');
            SinD(:,:,1:2:Nz) = SinogD(:,:,1:ringsp+1:ringsp^2);
            for jh=1:floor(span/2)
                apu = SinogD(:,:,jh*ringsp+1:ringsp+1:ringsp^2);
                apu2 = SinogD(:,:,jh+1:ringsp+1:(ringsp-jh)*ringsp);
                SinD(:,:,jh+1:2:offset2(1)-jh) = SinD(:,:,jh+1:2:offset2(1)-jh) + apu + apu2;
            end
            for ih=1:floor(length(segment_table)/2)
                for jh=1:span
                    apu = SinogD(:,:,(kkj(ih)+jh-1)*ringsp+1:ringsp+1:end);
                    SinD(:,:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) = SinD(:,:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) + (apu);
                    apu2 = SinogD(:,:,kkj(ih)+jh:ringsp+1:(ringsp-kkj(ih)-jh+1)*ringsp);
                    SinD(:,:,offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1) = SinD(:,:,offset2(2*ih)+jh:2:offset2(2*ih + 1)-jh+1) + (apu2);
                end
            end
            Sin = Sin - SinD;
            if partitions > 1
                SinDelayed{llo} = SinD;
            else
                SinDelayed = SinD;
            end
            varargout{3} = SinD;
            clear SinogD SinD
        elseif options.randoms_correction && ~options.use_ASCII && ~options.use_LMF && ~options.use_root && ~options.use_machine
            if iscell(options.SinD)
                if size(options.SinD{llo},2) ~= Nang && size(options.SinD{llo},1) ~= Ndist
                    options.SinD{llo} = reshape(options.SinD{llo}, Ndist,Nang,NSlices);
                elseif size(options.SinD{llo},1)*size(options.SinD{llo},2)*size(options.SinD{llo},3) ~= Ndist*Nang*NSlices
                    error('Size mismatch between randoms sinogram and and specified sinogram')
                end
                Sin = Sin - uint16(options.SinD{llo});
            else
                if size(options.SinD,2) ~= Nang && size(options.SinD,1) ~= Ndist
                    options.SinD = reshape(options.SinD, Ndist,Nang,NSlices);
                elseif size(options.SinD,1)*size(options.SinD,2)*size(options.SinD,3) ~= Ndist*Nang*NSlices
                    error('Size mismatch between randoms sinogram and and specified sinogram')
                end
                Sin = Sin - uint16(options.SinD);
            end
        end
        
        % Apply scatter correction
        if options.scatter_correction
            if iscell(options.ScatterC)
                if size(options.ScatterC{llo},2) ~= Nang && size(options.ScatterC{llo},1) ~= Ndist && numel(options.ScatterC{llo}) == numel(Sin)
                    options.ScatterC{llo} = reshape(options.ScatterC{llo}, Ndist,Nang,NSlices);
                elseif size(options.ScatterC{llo},1)*size(options.ScatterC{llo},2)*size(options.ScatterC{llo},3) ~= Ndist*Nang*NSlices
                    error('Size mismatch between scatter sinogram and and specified sinogram')
                end
                Sin = Sin - uint16(options.ScatterC{llo});
%                 if partitions > 1
%                     SinSC{llo} = SinC;
%                 else
%                     SinSC = SinC;
%                 end
            else
                if size(options.ScatterC,2) ~= Nang && size(options.ScatterC,1) ~= Ndist && numel(options.ScatterC) == numel(Sin)
                    options.ScatterC = reshape(options.ScatterC, Ndist,Nang,NSlices);
                elseif size(options.ScatterC,1)*size(options.ScatterC,2)*size(options.ScatterC,3) ~= Ndist*Nang*NSlices
                    error('Size mismatch between scatter sinogram and and specified sinogram')
                end
                Sin = Sin - uint16(options.ScatterC);
%                 if partitions > 1
%                     SinSC{llo} = SinC;
%                 else
%                     SinSC = SinC;
%                 end
            end
        end
        if options.randoms_correction || options.scatter_correction || options.normalization_correction
            if partitions > 1
                SinM{llo} = Sin;
            else
                SinM = Sin;
            end
        end
        
        % Form the sinogram of true coincidences
        if options.obtain_trues && options.use_machine == 0
            SinogT = cat(3,SinogT{:});
            Sin = zeros(Ndist,Nang,NSlices,'uint16');
            Sin(:,:,1:2:Nz) = SinogT(:,:,1:ringsp+1:ringsp^2);
            for jh=1:floor(span/2)
                apu = SinogT(:,:,jh*ringsp+1:ringsp+1:ringsp^2);
                apu2 = SinogT(:,:,jh+1:ringsp+1:(ringsp-jh)*ringsp);
                Sin(:,:,jh+1:2:offset2(1)-jh) = Sin(:,:,jh+1:2:offset2(1)-jh) + apu + apu2;
            end
            for ih=1:floor(length(segment_table)/2)
                for jh=1:span
                    apu = SinogT(:,:,(kkj(ih)+jh-1)*ringsp+1:ringsp+1:end);
                    Sin(:,:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) = Sin(:,:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) + (apu);
                    apu2 = SinogT(:,:,kkj(ih)+jh:ringsp+1:(ringsp-kkj(ih)-jh+1)*ringsp);
                    Sin(:,:,offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1) = Sin(:,:,offset2(2*ih)+jh:2:offset2(2*ih + 1)-jh+1) + (apu2);
                end
            end
            if partitions > 1
                SinTrues{llo} = Sin;
            else
                SinTrues = Sin;
            end
            varargout{4} = Sin;
        end
        % Form the sinogram of scattered coincidences
        if options.store_scatter && options.use_machine == 0
            SinogS = cat(3,SinogS{:});
            Sin = zeros(Ndist,Nang,NSlices,'uint16');
            Sin(:,:,1:2:Nz) = SinogS(:,:,1:ringsp+1:ringsp^2);
            for jh=1:floor(span/2)
                apu = SinogS(:,:,jh*ringsp+1:ringsp+1:ringsp^2);
                apu2 = SinogS(:,:,jh+1:ringsp+1:(ringsp-jh)*ringsp);
                Sin(:,:,jh+1:2:offset2(1)-jh) = Sin(:,:,jh+1:2:offset2(1)-jh) + apu + apu2;
            end
            for ih=1:floor(length(segment_table)/2)
                for jh=1:span
                    apu = SinogS(:,:,(kkj(ih)+jh-1)*ringsp+1:ringsp+1:end);
                    Sin(:,:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) = Sin(:,:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) + (apu);
                    apu2 = SinogS(:,:,kkj(ih)+jh:ringsp+1:(ringsp-kkj(ih)-jh+1)*ringsp);
                    Sin(:,:,offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1) = Sin(:,:,offset2(2*ih)+jh:2:offset2(2*ih + 1)-jh+1) + (apu2);
                end
            end
            if partitions > 1
                SinScatter{llo} = Sin;
            else
                SinScatter = Sin;
            end
            varargout{5} = Sin;
        end
        % Form the sinogram of true random coincidences
        if options.store_randoms && options.use_machine == 0
            SinogR = cat(3,SinogR{:});
            Sin = zeros(Ndist,Nang,NSlices,'uint16');
            Sin(:,:,1:2:Nz) = SinogR(:,:,1:ringsp+1:ringsp^2);
            for jh=1:floor(span/2)
                apu = SinogR(:,:,jh*ringsp+1:ringsp+1:ringsp^2);
                apu2 = SinogR(:,:,jh+1:ringsp+1:(ringsp-jh)*ringsp);
                Sin(:,:,jh+1:2:offset2(1)-jh) = Sin(:,:,jh+1:2:offset2(1)-jh) + apu + apu2;
            end
            for ih=1:floor(length(segment_table)/2)
                for jh=1:span
                    apu = SinogR(:,:,(kkj(ih)+jh-1)*ringsp+1:ringsp+1:end);
                    Sin(:,:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) = Sin(:,:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) + (apu);
                    apu2 = SinogR(:,:,kkj(ih)+jh:ringsp+1:(ringsp-kkj(ih)-jh+1)*ringsp);
                    Sin(:,:,offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1) = Sin(:,:,offset2(2*ih)+jh:2:offset2(2*ih + 1)-jh+1) + (apu2);
                end
            end
            if partitions > 1
                SinRandoms{llo} = Sin;
            else
                SinRandoms = Sin;
            end
            varargout{6} = Sin;
        end
        if options.verbose
            disp(['Sinogram for timestep ' num2str(llo) ' formed'])
        end
    end
    %%
    if partitions == 1
        if options.use_machine == 0
            save_string = [machine_name '_' name '_sinograms_combined_static_' num2str(Ndist) 'x' num2str(Nang) 'x' num2str(NSlices) '_span' num2str(span) '.mat'];
        else
            save_string = [machine_name '_' name '_sinograms_combined_static_' num2str(Ndist) 'x' num2str(Nang) 'x' num2str(NSlices) '_span' num2str(span) '_listmode.mat'];
        end
    else
        if options.use_machine == 0
            save_string = [machine_name '_' name '_sinograms_combined_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_' num2str(Ndist) 'x' num2str(Nang) 'x' num2str(NSlices) '_span' num2str(span) '.mat'];
        else
            save_string = [machine_name '_' name '_sinograms_combined_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_' num2str(Ndist) 'x' num2str(Nang) 'x' num2str(NSlices) '_span' num2str(span) '_listmode.mat'];
        end
    end
    save(save_string,'raw_SinM')
    for variableIndex = 1:length(variableList)
        if exist(variableList{variableIndex},'var')
            save(save_string,variableList{variableIndex},'-append')
        end
    end
    if (options.randoms_correction || options.scatter_correction || options.normalization_correction) && ~options.corrections_during_reconstruction
        varargout{1} = SinM;
        varargout{2} = raw_SinM;
    else
        varargout{1} = raw_SinM;
    end
    if options.verbose
        disp('Sinograms saved')
        toc
    end
    
else
    % Load Inveon sinogram file
    if nargout > 1
        error('Inveon sinogram format can only have one output variable')
    end
    [options.file, options.fpath] = uigetfile('*.scn','Select Inveon sinogram datafile');
    if isequal(options.file, 0)
        error('No file was selected')
    end
    
    if options.verbose
        tic
    end
    
    nimi = [options.fpath options.file];
    fid = fopen(nimi);
%     SinM = fread(fid, inf, 'int32=>uint16',0,'l');
    SinM = fread(fid, inf, 'single=>single',0,'l');
    tof_length = length(SinM)/(options.Nang*options.Ndist*options.TotSinos);
    SinM = reshape(SinM,options.Ndist,options.Nang,options.TotSinos,tof_length);
    if options.partitions == 1
        SinM = sum(SinM,4,'native');
    end
%     Sin = SinM;
    varargout{1} = SinM;
    
    %%
    save_string = [options.machine_name '_' options.name '_sinogram_original_static_' num2str(options.Ndist) 'x' num2str(options.Nang)...
        'x' num2str(options.TotSinos) '_span' num2str(options.span) '_machine_sinogram.mat'];
    save(save_string,'SinM')
    if options.verbose
        disp('Sinograms formed')
        toc
    end
end