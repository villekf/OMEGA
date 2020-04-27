function [varargout] = form_sinograms(options, varargin)
%% FORM MICHELOGRAMS AND SINOGRAMS FROM RAW DETECTOR PAIR DATA
% This code forms the sinograms for the current machine from the raw
% list-mode data using the raw detector pairs. It first creates the
% Michelogram and later the sinograms with the dimensions specified by the
% user.
%
% Can also be used to simply perform the selected corrections to the
% specified uncorrected sinogram.
%
% Input the machine and sinogram properties and the raw detector pair
% data.
% Output is the sinograms for each time point, saved in a mat-file in the
% current working directory. Any previous corrections are overwritten.
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
%   SinDelayed = Delayed coincidence sinogram (with optional variance
%   reduction and/or smoothing)
%   raw_SinDelayed = Delayed coincidence sinogram (no modifications)
%   SinTrues = True coincidences (GATE only)
%   SinScatter = Scattered coincidence sinogram (GATE only)
%   SinRandoms = Randoms coincidence sinogram (GATE only)
%
% See also michelogram_formation, load_data, sinogram_coordinates_2D,
% sinogram_coordinates_3D, gapFilling, arcCorrection



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Ville-Veikko Wettenhovi
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout > 5
    error('Too many output arguments')
end
if nargin > 1
    corrections = varargin{1};
else
    corrections = false;
end
if nargin >= 3
    data_load = true;
else
    data_load = false;
end

folder = fileparts(which('form_sinograms.m'));
folder = strrep(folder, 'source','mat-files/');
folder = strrep(folder, '\','/');

% Variables that specify which corrections are applied
appliedCorrections.normalization = false;
appliedCorrections.randoms = '';
appliedCorrections.scatter = '';
appliedCorrections.gapFilling = false;
appliedCorrections.mashing = 1;
appliedCorrections.globalFactor = options.global_correction_factor;
ScatterProp.smoothing = false;

if options.use_machine ~= 2
    
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
    
    ringsp = rings;
    ScatterProp.smoothing = false;
    RandProp.variance_reduction = false;
    RandProp.smoothing = false;
    
    temp = pseudot;
    if ~isempty(temp) && temp > 0
        for kk = uint32(1) : temp
            pseudot(kk) = uint32(options.cryst_per_block + 1) * kk;
        end
    elseif temp == 0
        pseudot = [];
    end
    
    % Sinograms are reformed
    if ~corrections && ~data_load
        variableList = {'SinM','SinTrues','SinScatter','SinRandoms','SinDelayed','raw_SinDelayed','SinSC','ScatterProp','RandProp','appliedCorrections'};
        
        [GATE_vars, I] = sort([options.use_ASCII, options.use_root, options.use_LMF],'descend');
        GATE_char = {'ASCII';'root';'LMF'};
        GATE_char = GATE_char(I);
        
        % Load the data and form the initial Michelograms
        if options.verbose
            disp('Forming initial Michelogram')
            tic
        end
        if options.use_machine == 1
            if exist( 'options.coincidences' , 'var') == 0
                if options.partitions == 1
                    options.coincidences = loadStructFromFile([machine_name '_measurements_' name '_static_raw_listmode.mat'], 'coincidences');
                else
                    options.coincidences = loadStructFromFile([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_listmode.mat'], 'coincidences');
                end
            end
            if options.randoms_correction
                if options.partitions == 1
                    load([machine_name '_measurements_' name '_static_raw_listmode.mat'], 'delayed_coincidences')
                else
                    load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_listmode.mat'], 'delayed_coincidences')
                end
                delayed_coincidences = initial_michelogram(options, delayed_coincidences, 'delays');
            end
        else
            if exist( 'options.coincidences' , 'var') == 0
                if options.partitions == 1
                    if any(GATE_vars)
                        filename = [machine_name '_measurements_' name '_static_raw_'];
                        options.coincidences = loadGATEData(filename,'coincidences', GATE_char);
                    else
                        options.coincidences = loadStructFromFile([machine_name '_measurements_' name '_static_raw_real.mat'], 'coincidences');
                    end
                else
                    if any(GATE_vars)
                        filename = [machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_'];
                        options.coincidences = loadGATEData(filename,'coincidences', GATE_char);
                    else
                        options.coincidences = loadStructFromFile([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_real.mat'], 'coincidences');
                    end
                end
            end
            if options.obtain_trues
                if options.partitions == 1
                    if any(GATE_vars)
                        filename = [machine_name '_measurements_' name '_static_raw_'];
                        true_coincidences = loadGATEData(filename,'true_coincidences', GATE_char);
                    else
                        load([machine_name '_measurements_' name '_static_raw_real.mat'], 'true_coincidences')
                    end
                else
                    if any(GATE_vars)
                        filename = [machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_'];
                        true_coincidences = loadGATEData(filename,'true_coincidences', GATE_char);
                    else
                        load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_real.mat'], 'true_coincidences')
                    end
                end
                true_coincidences = initial_michelogram(options, true_coincidences, 'trues');
            end
            if options.store_scatter
                if options.partitions == 1
                    if any(GATE_vars)
                        filename = [machine_name '_measurements_' name '_static_raw_'];
                        scattered_coincidences = loadGATEData(filename,'scattered_coincidences', GATE_char);
                    else
                        load([machine_name '_measurements_' name '_static_raw_real.mat'], 'scattered_coincidences')
                    end
                else
                    if any(GATE_vars)
                        filename = [machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_'];
                        scattered_coincidences = loadGATEData(filename,'scattered_coincidences', GATE_char);
                    else
                        load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_real.mat'], 'scattered_coincidences')
                    end
                end
                scattered_coincidences = initial_michelogram(options, scattered_coincidences, 'scatter');
            end
            if options.store_randoms
                if options.partitions == 1
                    if any(GATE_vars)
                        filename = [machine_name '_measurements_' name '_static_raw_'];
                        random_coincidences = loadGATEData(filename,'random_coincidences', GATE_char);
                    else
                        load([machine_name '_measurements_' name '_static_raw_real.mat'], 'random_coincidences')
                    end
                else
                    if any(GATE_vars)
                        filename = [machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_'];
                        random_coincidences = loadGATEData(filename,'random_coincidences', GATE_char);
                    else
                        load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_real.mat'], 'random_coincidences')
                    end
                end
                random_coincidences = initial_michelogram(options, random_coincidences, 'randoms');
            end
            
            if options.randoms_correction && (options.use_ASCII || options.use_LMF || options.use_root)
                if options.partitions == 1
                    if any(GATE_vars)
                        filename = [machine_name '_measurements_' name '_static_raw_'];
                        delayed_coincidences = loadGATEData(filename,'delayed_coincidences', GATE_char);
                    else
                        load([machine_name '_measurements_' name '_static_raw_real.mat'], 'delayed_coincidences')
                    end
                else
                    if any(GATE_vars)
                        filename = [machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_'];
                        delayed_coincidences = loadGATEData(filename,'delayed_coincidences', GATE_char);
                    else
                        load([machine_name '_measurements_' name '_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_raw_real.mat'], 'delayed_coincidences')
                    end
                end
                delayed_coincidences = initial_michelogram(options, delayed_coincidences, 'delays');
            end
        end
        options.coincidences = initial_michelogram(options, options.coincidences, 'prompts');
        
        if options.verbose
            disp('Initial Michelogram formed')
            toc
            tic
        end
    elseif ~data_load
        % Corrections are applied
        if options.partitions == 1
            if options.use_machine == 0
                sinoFile = [options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                    num2str(options.TotSinos) '_span' num2str(options.span) '.mat'];
            elseif options.use_machine == 1
                sinoFile = [options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                    num2str(options.TotSinos) '_span' num2str(options.span) '_listmode.mat'];
            elseif options.use_machine == 2
                sinoFile = [options.machine_name '_' options.name '_sinogram_original_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                    num2str(options.TotSinos) '_span' num2str(options.span) '_machine_sinogram.mat'];
            elseif options.use_machine == 3
                sinoFile = [options.machine_name '_' options.name '_sinogram_original_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                    num2str(options.TotSinos) '_span' num2str(options.span) '_listmode_sinogram.mat'];
            end
        else
            if options.use_machine == 0
                sinoFile = [options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                    num2str(options.span) '.mat'];
            elseif options.use_machine == 1
                sinoFile = [options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                    num2str(options.span) '_listmode.mat'];
            elseif options.use_machine == 2
                sinoFile = [options.machine_name '_' options.name '_sinograms_original_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                    num2str(options.span) '_machine_sinogram.mat'];
            elseif options.use_machine == 3
                sinoFile = [options.machine_name '_' options.name '_sinograms_original_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                    num2str(options.span) '_listmode_sinogram.mat'];
            end
        end
        Sino = loadStructFromFile(sinoFile, 'raw_SinM');
        variableList = {'appliedCorrections'};
        if options.randoms_correction
            try
                RandProp = loadStructFromFile(sinoFile, 'RandProp');
                variableList = [variableList, {'RandProp'}];
            end
            if RandProp.variance_reduction || RandProp.smoothing
                try
                    SinDelayed = loadStructFromFile(sinoFile, 'SinDelayed');
                catch
                    RandProp.variance_reduction = false;
                    RandProp.smoothing = false;
                    try
                        SinDelayed = loadStructFromFile(sinoFile, 'raw_SinDelayed');
                    catch
                        options = loadDelayedData(options);
                    end
                end
            else
                try
                    SinDelayed = loadStructFromFile(sinoFile, 'raw_SinDelayed');
                catch
                    options = loadDelayedData(options);
                end
            end
        end
        if options.scatter_correction
            options.ScatterC = loadStructFromFile(sinoFile, 'SinScatter');
            try
                ScatterProp = loadStructFromFile(sinoFile, 'ScatterProp');
                variableList = [variableList, {'ScatterProp'}];
            end
        end
        SinM = cell(partitions,1);
    elseif data_load
        % Sinograms are formed when loading raw data
        variableList = {'SinM','SinTrues','SinScatter','SinRandoms','SinDelayed','raw_SinDelayed','SinSC','ScatterProp','RandProp','appliedCorrections'};
        
        if options.verbose
            disp('Forming initial Michelogram')
            tic
        end
        options.coincidences = initial_michelogram(options, varargin{2}, 'prompts', true);
        if options.obtain_trues && options.use_machine == 0 && ~isempty(varargin{4})
            true_coincidences = initial_michelogram(options, varargin{4}, 'trues', true);
        end
        if options.store_scatter && options.use_machine == 0 && ~isempty(varargin{5})
            scattered_coincidences = initial_michelogram(options, varargin{5}, 'scatter', true);
        end
        if options.store_randoms && options.use_machine == 0 && ~isempty(varargin{6})
            random_coincidences = initial_michelogram(options, varargin{6}, 'randoms', true);
        end
        if options.randoms_correction && (options.use_ASCII || options.use_LMF || options.use_root || options.use_machine == 1) && ~isempty(varargin{3})
            delayed_coincidences = initial_michelogram(options, varargin{3}, 'delays', true);
        end
        if options.verbose
            disp('Initial Michelogram formed')
            toc
            tic
        end
    end
    
    [~, ~, xp, yp] = detector_coordinates(options);
    if ~corrections
        [~, ~, i, j, accepted_lors] = sinogram_coordinates_2D(options, xp, yp);
        
        
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
                if options.randoms_smoothing || options.variance_reduction
                    raw_SinDelayed = cell(partitions,1);
                end
            end
            raw_SinM = cell(partitions,1);
        end
    end
    
    % Go through all the time steps
    for llo = 1 : partitions
        
        if ~corrections
            if options.verbose
                disp(['Sinogram formation at timestep ' num2str(llo) ' began'])
            end
            
            if partitions > 1 || ~data_load
                P1 = options.coincidences{llo};
            else
                P1 = options.coincidences;
            end
            
            Sinog = cell(ringsp,ringsp);
            
            ix=cellfun('isempty',Sinog);
            Sinog(ix)={zeros(Ndist,Nang,'uint16')};
            
            if options.obtain_trues && options.use_machine == 0
                if partitions > 1 || ~data_load
                    T1 = true_coincidences{llo};
                else
                    T1 = true_coincidences;
                end
                SinogT = cell(ringsp,ringsp);
                SinogT(ix)={zeros(Ndist,Nang,'uint16')};
            end
            if options.store_scatter && options.use_machine == 0
                if partitions > 1 || ~data_load
                    S1 = scattered_coincidences{llo};
                else
                    S1 = scattered_coincidences;
                end
                SinogS = cell(ringsp,ringsp);
                SinogS(ix)={zeros(Ndist,Nang,'uint16')};
            end
            if options.store_randoms && options.use_machine == 0
                if partitions > 1 || ~data_load
                    R1 = random_coincidences{llo};
                else
                    R1 = random_coincidences;
                end
                SinogR = cell(ringsp,ringsp);
                SinogR(ix)={zeros(Ndist,Nang,'uint16')};
            end
            if options.randoms_correction && (options.use_ASCII || options.use_LMF || options.use_root || options.use_machine == 1)
                if partitions > 1 || ~data_load
                    D1 = delayed_coincidences{llo};
                else
                    D1 = delayed_coincidences;
                end
                SinogD = cell(ringsp,ringsp);
                SinogD(ix)={zeros(Ndist,Nang,'uint16')};
            end
            
            % Create the Michelograms
            tic
            for ii=1:ringsp
                if any(ii == pseudot)
                    continue
                end
                for jj=1:ringsp
                    if any(jj == pseudot)
                        continue
                    end
                    if issparse(P1{ii,jj})
                        CC = uint16(full(P1{ii,jj}));
                    else
                        CC = P1{ii,jj};
                    end
                    CC = CC(accepted_lors);
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
            clear Sinog
            if partitions > 1
                raw_SinM{llo} = Sin;
            else
                raw_SinM = Sin;
            end
        else
            if iscell(Sino)
                Sin = Sino{llo};
            else
                Sin = Sino;
            end
        end
        % Apply randoms correction
        if options.randoms_correction && (options.use_ASCII || options.use_LMF || options.use_root || options.use_machine == 1) ...
                && (exist('SinDelayed','var') || exist('SinogD','var'))
            if ~corrections
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
                if partitions > 1
                    if options.variance_reduction || options.randoms_smoothing
                        raw_SinDelayed{llo} = SinD;
                    end
                    SinDelayed{llo} = SinD;
                else
                    if options.variance_reduction || options.randoms_smoothing
                        raw_SinDelayed = SinD;
                    end
                    SinDelayed = SinD;
                end
                varargout{4} = SinD;
            else
                if iscell(SinDelayed)
                    SinD = SinDelayed{llo};
                else
                    SinD = SinDelayed;
                end
            end
            if (options.variance_reduction || options.randoms_smoothing) && ~options.corrections_during_reconstruction
                if options.variance_reduction && ~RandProp.variance_reduction
                    SinD = single(Randoms_variance_reduction(double(SinD), options));
                    RandProp.variance_reduction = true;
                    appliedCorrections.randoms = 'randoms correction with variance reduction';
                end
                if options.randoms_smoothing && ~RandProp.smoothing
                    SinD = single(randoms_smoothing(double(SinD), options));
                    RandProp.smoothing = true;
                    if options.variance_reduction
                        appliedCorrections.randoms = 'and smoothing';
                    else
                        appliedCorrections.randoms = 'randoms correction with smoothing';
                    end
                end
                if ~options.corrections_during_reconstruction
                    Sin = single(Sin) - single(SinD);
                end
                if ~corrections
                    if partitions > 1
                        SinDelayed{llo} = SinD;
                    else
                        SinDelayed = SinD;
                    end
                end
            elseif ~options.corrections_during_reconstruction
                Sin = single(Sin) - single(SinD);
            end
            if ~corrections
                varargout{3} = SinD;
            end
            clear SinogD SinD
        elseif options.randoms_correction && ~options.corrections_during_reconstruction && isfield(options,'SinDelayed')
            Sin = single(Sin);
            if iscell(options.SinDelayed)
                if size(options.SinDelayed{llo},2) ~= Nang && size(options.SinDelayed{llo},1) ~= Ndist
                    options.SinDelayed{llo} = reshape(options.SinDelayed{llo}, Ndist,Nang,NSlices);
                elseif size(options.SinDelayed{llo},1)*size(options.SinDelayed{llo},2)*size(options.SinDelayed{llo},3) ~= Ndist*Nang*NSlices
                    error('Size mismatch between randoms sinogram and and specified sinogram')
                end
                if options.variance_reduction && ~RandProp.variance_reduction
                    options.SinDelayed{llo} = Randoms_variance_reduction(double(options.SinDelayed{llo}), options);
                    RandProp.variance_reduction = true;
                    appliedCorrections.randoms = 'randoms correction with variance reduction';
                end
                if options.randoms_smoothing && ~RandProp.smoothing
                    options.SinDelayed{llo} = randoms_smoothing(options.SinDelayed{llo}, options);
                    RandProp.smoothing = true;
                    if options.variance_reduction
                        appliedCorrections.randoms = 'and smoothing';
                    else
                        appliedCorrections.randoms = 'randoms correction with smoothing';
                    end
                end
                Sin = Sin - single(options.SinDelayed{llo});
            else
                if size(options.SinDelayed,2) ~= Nang && size(options.SinDelayed,1) ~= Ndist
                    options.SinDelayed = reshape(options.SinDelayed, Ndist,Nang,NSlices);
                elseif size(options.SinDelayed,1)*size(options.SinDelayed,2)*size(options.SinDelayed,3) ~= Ndist*Nang*NSlices
                    error('Size mismatch between randoms sinogram and and specified sinogram')
                end
                if options.variance_reduction && ~RandProp.variance_reduction
                    options.SinDelayed = Randoms_variance_reduction(double(options.SinDelayed), options);
                    RandProp.variance_reduction = true;
                    appliedCorrections.randoms = 'randoms correction with variance reduction';
                end
                if options.randoms_smoothing && ~RandProp.smoothing
                    options.SinDelayed = randoms_smoothing(options.SinDelayed, options);
                    RandProp.smoothing = true;
                    if options.variance_reduction
                        appliedCorrections.randoms = 'and smoothing';
                    else
                        appliedCorrections.randoms = 'randoms correction with smoothing';
                    end
                end
                Sin = Sin - single(options.SinDelayed);
            end
            if options.verbose
                disp('Randoms correction applied to sinogram')
            end
        end
        
        % Apply scatter correction
        if options.scatter_correction && ~options.corrections_during_reconstruction
            Sin = single(Sin);
            if options.scatter_correction && ~isfield(options,'ScatterC')
                options = loadScatterData(options);
            end
            if iscell(options.ScatterC)
                if size(options.ScatterC{llo},2) ~= Nang && size(options.ScatterC{llo},1) ~= Ndist && numel(options.ScatterC{llo}) == numel(Sin)
                    options.ScatterC{llo} = reshape(options.ScatterC{llo}, Ndist,Nang,NSlices);
                elseif size(options.ScatterC{llo},1)*size(options.ScatterC{llo},2)*size(options.ScatterC{llo},3) ~= Ndist*Nang*NSlices
                    error('Size mismatch between scatter sinogram and specified sinogram')
                end
                if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                    options.ScatterC{llo} = Randoms_variance_reduction(double(options.ScatterC{llo}), options);
                    ScatterProp.variance_reduction = true;
                    appliedCorrections.scatter = 'scatter correction with variance reduction';
                end
                if options.scatter_smoothing && ~ScatterProp.smoothing
                    options.ScatterC{llo} = randoms_smoothing(options.ScatterC{llo}, options);
                    ScatterProp.smoothing = true;
                    if options.scatter_variance_reduction
                        appliedCorrections.scatter = 'and smoothing';
                    else
                        appliedCorrections.scatter = 'scatter correction with smoothing';
                    end
                else
                    if ~options.scatter_variance_reduction
                        appliedCorrections.scatter = 'scatter correction';
                    end
                end
                Sin = Sin - single(options.ScatterC{llo});
            else
                if size(options.ScatterC,2) ~= Nang && size(options.ScatterC,1) ~= Ndist && numel(options.ScatterC) == numel(Sin)
                    options.ScatterC = reshape(options.ScatterC, Ndist,Nang,NSlices);
                elseif size(options.ScatterC,1)*size(options.ScatterC,2)*size(options.ScatterC,3) ~= Ndist*Nang*NSlices
                    error('Size mismatch between scatter sinogram and specified sinogram')
                end
                if options.variance_reduction && ~ScatterProp.variance_reduction
                    options.ScatterC = Randoms_variance_reduction(double(options.ScatterC), options);
                    ScatterProp.variance_reduction = true;
                    appliedCorrections.scatter = 'scatter correction with variance reduction';
                end
                if options.scatter_smoothing && ~ScatterProp.smoothing
                    options.ScatterC = randoms_smoothing(options.ScatterC, options);
                    ScatterProp.smoothing = true;
                    if options.scatter_variance_reduction
                        appliedCorrections.scatter = 'and smoothing';
                    else
                        appliedCorrections.scatter = 'scatter correction with smoothing';
                    end
                else
                    if ~options.scatter_variance_reduction
                        appliedCorrections.scatter = 'scatter correction';
                    end
                end
                Sin = Sin - single(options.ScatterC);
            end
            if options.verbose
                disp('Scatter correction applied to sinogram')
            end
        end
        % Apply normalization correction
        if options.normalization_correction && ~options.corrections_during_reconstruction
            Sin = single(Sin);
            if options.use_user_normalization
                [n_file, fpath] = uigetfile({'*.nrm;*.n;*.ptd;*.mat'},'Select normalization datafile');
                if isequal(n_file, 0)
                    error('No file was selected')
                end
                f_path = [fpath n_file];
                if any(strfind(n_file, '.nrm'))
                    fid = fopen(f_path);
                    normalization = fread(fid, inf, 'single=>single',0,'l');
                    fclose(fid);
                    normalization = reshape(normalization, options.Ndist, options.Nang, options.TotSinos);
                elseif any(strfind(n_file, '.n')) || any(strfind(n_file, '.ptd'))
                    normalization = loadNormalizationBiograph(f_path, options);
                else
                    data = load(f_path);
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
                norm_file = [folder options.machine_name '_normalization_' num2str(options.Ndist) 'x' num2str(options.Nang) '_span' num2str(options.span) '.mat'];
                if exist(norm_file, 'file') == 2
                    normalization = loadStructFromFile(norm_file,'normalization');
                else
                    normalization = normalization_coefficients(options);
                    %                         normalization = loadStructFromFile(norm_file,'normalization');
                end
                %                 options.normalization = options.normalization(:);
                %                 if (options.implementation == 1 || options.implementation == 4)
                %                     options.normalization = single(options.normalization);
                %                     options.normalization = 1 ./ options.normalization;
                %                 else
                %                     options.normalization = single(1) ./ options.normalization;
                %                 end
            end
            Sin = Sin .* normalization;
            if options.verbose
                disp('Normalization correction applied to sinogram')
            end
            appliedCorrections.normalization = true;
        end
        % Normalized scatter
        if options.scatter_correction && ~options.corrections_during_reconstruction && options.normalize_scatter
            Sin = single(Sin);
            if options.scatter_correction && ~isfield(options,'ScatterC')
                options = loadScatterData(options);
            end
            if iscell(options.ScatterC)
                if size(options.ScatterC{llo},2) ~= Nang && size(options.ScatterC{llo},1) ~= Ndist && numel(options.ScatterC{llo}) == numel(Sin)
                    options.ScatterC{llo} = reshape(options.ScatterC{llo}, Ndist,Nang,NSlices);
                elseif size(options.ScatterC{llo},1)*size(options.ScatterC{llo},2)*size(options.ScatterC{llo},3) ~= Ndist*Nang*NSlices
                    error('Size mismatch between scatter sinogram and specified sinogram')
                end
                if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization
                    options.ScatterC{kk} = single(options.ScatterC{kk}) .* reshape(normalization, size(options.ScatterC{kk}, 1), size(options.ScatterC{kk}, 2), size(options.ScatterC{kk}, 3));
                    ScatterProp.normalization = true;
                end
                if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                    options.ScatterC{llo} = Randoms_variance_reduction(double(options.ScatterC{llo}), options);
                    ScatterProp.variance_reduction = true;
                    appliedCorrections.scatter = 'normalized scatter correction with variance reduction';
                end
                if options.scatter_smoothing && ~ScatterProp.smoothing
                    options.ScatterC{llo} = randoms_smoothing(options.ScatterC{llo}, options);
                    ScatterProp.smoothing = true;
                    if options.scatter_variance_reduction
                        appliedCorrections.scatter = 'and smoothing';
                    else
                        appliedCorrections.scatter = 'normalized scatter correction with smoothing';
                    end
                else
                    if ~options.scatter_variance_reduction
                        appliedCorrections.scatter = 'normalized scatter correction';
                    end
                end
                Sin = Sin - single(options.ScatterC{llo});
            else
                if size(options.ScatterC,2) ~= Nang && size(options.ScatterC,1) ~= Ndist && numel(options.ScatterC) == numel(Sin)
                    options.ScatterC = reshape(options.ScatterC, Ndist,Nang,NSlices);
                elseif size(options.ScatterC,1)*size(options.ScatterC,2)*size(options.ScatterC,3) ~= Ndist*Nang*NSlices
                    error('Size mismatch between scatter sinogram and specified sinogram')
                end
                if options.normalization_correction && options.normalize_scatter && ~ScatterProp.normalization
                    options.ScatterC = single(options.ScatterC) .* reshape(normalization, size(options.ScatterC, 1), size(options.ScatterC, 2), size(options.ScatterC, 3));
                    ScatterProp.normalization = true;
                end
                if options.scatter_variance_reduction && ~ScatterProp.variance_reduction
                    options.ScatterC = Randoms_variance_reduction(double(options.ScatterC), options);
                    ScatterProp.variance_reduction = true;
                    appliedCorrections.scatter = 'normalized scatter correction with variance reduction';
                end
                if options.scatter_smoothing && ~ScatterProp.smoothing
                    options.ScatterC = randoms_smoothing(options.ScatterC, options);
                    ScatterProp.smoothing = true;
                    if options.scatter_variance_reduction
                        appliedCorrections.scatter = 'and smoothing';
                    else
                        appliedCorrections.scatter = 'normalized scatter correction with smoothing';
                    end
                else
                    if ~options.scatter_variance_reduction
                        appliedCorrections.scatter = 'normalized scatter correction';
                    end
                end
                Sin = Sin - single(options.ScatterC);
            end
            if options.verbose
                disp('Scatter correction applied to sinogram')
            end
        end
        if isfield(options,'global_correction_factor') && ~isempty(options.global_correction_factor) && options.global_correction_factor > 0 && ~options.corrections_during_reconstruction
            Sin = Sin * options.global_correction_factor;
        end
        % Fill sinogram gaps if pseudo detectors are present and the user
        % has seleced gap filling
        
        if options.fill_sinogram_gaps && sum(pseudot) > 0
            if llo == 1
                gaps = [];
            end
            [Sin, gaps] = gapFilling(options, Sin, xp, yp, llo, gaps);
            if partitions > 1
                SinM{llo} = Sin;
            else
                SinM = Sin;
            end
            appliedCorrections.gapFilling = true;
        end
        if (options.randoms_correction || options.scatter_correction || options.normalization_correction) && ~options.corrections_during_reconstruction
            Sin(Sin < 0) = 0;
            if partitions > 1
                SinM{llo} = Sin;
            else
                SinM = Sin;
            end
        end
        
        if ~corrections
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
                varargout{5} = Sin;
                if options.verbose
                    disp('Trues stored')
                end
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
                varargout{6} = Sin;
                if options.verbose
                    disp('Scattered coincidences stored')
                end
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
                varargout{7} = Sin;
                if options.verbose
                    disp('Randoms stored')
                end
            end
            if options.verbose
                disp(['Sinogram for timestep ' num2str(llo) ' formed'])
            end
        end
    end
    %%
    if partitions == 1
        if options.use_machine == 0
            save_string = [machine_name '_' name '_sinograms_combined_static_' num2str(Ndist) 'x' num2str(Nang) 'x' num2str(NSlices) '_span' num2str(span) '.mat'];
        elseif options.use_machine == 3
            save_string = [machine_name '_' name '_sinogram_original_static_' num2str(Ndist) 'x' num2str(Nang) 'x' num2str(NSlices) '_span' num2str(span) '_listmode_sinogram.mat'];
        else
            save_string = [machine_name '_' name '_sinograms_combined_static_' num2str(Ndist) 'x' num2str(Nang) 'x' num2str(NSlices) '_span' num2str(span) '_listmode.mat'];
        end
    else
        if options.use_machine == 0
            save_string = [machine_name '_' name '_sinograms_combined_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_' ...
                num2str(Ndist) 'x' num2str(Nang) 'x' num2str(NSlices) '_span' num2str(span) '.mat'];
        elseif options.use_machine == 3
            save_string = [machine_name '_' name '_sinogram_original_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_' ...
                num2str(Ndist) 'x' num2str(Nang) 'x' num2str(NSlices) '_span' num2str(span) '_listmode_sinogram.mat'];
        else
            save_string = [machine_name '_' name '_sinograms_combined_' num2str(partitions) 'timepoints_for_total_of_ ' num2str(tot_time) 's_' ...
                num2str(Ndist) 'x' num2str(Nang) 'x' num2str(NSlices) '_span' num2str(span) '_listmode.mat'];
        end
    end
    if options.verbose
        disp('Saving sinogram...')
    end
    if ~corrections
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
        if (options.randoms_correction || options.scatter_correction || options.normalization_correction) && ~options.corrections_during_reconstruction
            varargout{1} = SinM;
            varargout{2} = raw_SinM;
        else
            varargout{1} = raw_SinM;
        end
        options = rmfield(options,'coincidences');
    else
        if exist('OCTAVE_VERSION','builtin') == 0
            save(save_string,'SinM','-append')
            for variableIndex = 1:length(variableList)
                if exist(variableList{variableIndex},'var')
                    save(save_string,variableList{variableIndex},'-append')
                end
            end
        else
            save(save_string,'SinM','-append','-v7')
            for variableIndex = 1:length(variableList)
                if exist(variableList{variableIndex},'var')
                    save(save_string,variableList{variableIndex},'-append','-v7')
                end
            end
        end
        varargout{1} = SinM;
    end
    if options.verbose
        disp('Sinograms saved')
        if ~corrections
            toc
        end
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
    
    appliedCorrections.randoms = 'randoms correction';
    
    if options.verbose
        tic
    end
    
    nimi = [options.fpath options.file];
    fid = fopen(nimi);
    %     SinM = fread(fid, inf, 'int32=>uint16',0,'l');
    raw_SinM = fread(fid, inf, 'int32=>int32',0,'l');
    tof_length = length(raw_SinM)/(options.Nang*options.Ndist*options.TotSinos);
    raw_SinM = reshape(raw_SinM,options.Ndist,options.Nang,options.TotSinos,tof_length);
    if options.partitions == 1
        raw_SinM = sum(raw_SinM,4,'native');
        if (options.scatter_correction || options.normalization_correction) && ~options.corrections_during_reconstruction
            SinM = raw_SinM;
        end
    else
        if options.partitions > tof_length
            error('Number of dynamic frames specified in options.partitions is larger than the number of frames in the sinogram')
        end
        temp = cell(tof_length, 1);
        for kk = options.partitions : -1 : 1
            temp{kk} = raw_SinM(:,:,:,kk);
            raw_SinM = raw_SinM(:,:,:,1:kk-1);
        end
        raw_SinM = temp;
        clear temp
        if (options.scatter_correction || options.normalization_correction) && ~options.corrections_during_reconstruction
            SinM = raw_SinM;
        end
    end
    
    if options.scatter_correction && ~options.corrections_during_reconstruction
        if ~isfield(options,'ScatterC')
            [scatter_file, s_fpath] = uigetfile('*.mat; *.scn','Select scatter correction data');
            if isequal(scatter_file, 0)
                error('No file was selected')
            end
            if strcmp(scatter_file(end-3:end), '.scn')
                nimi = [s_fpath scatter_file];
                fid = fopen(nimi);
                options.ScatterC = fread(fid, inf, 'single=>single',0,'l');
                s_length = length(options.ScatterC)/(options.Nang*options.Ndist*options.TotSinos);
                options.ScatterC = reshape(raw_SinM,options.Ndist,options.Nang,options.TotSinos,s_length);
                if s_length > 1
                    temp = cell(s_length, 1);
                    for kk = options.partitions : -1 : 1
                        temp{kk} = options.ScatterC(:,:,:,kk);
                        options.ScatterC = options.ScatterC(:,:,:,1:end-1);
                    end
                    options.ScatterC = temp;
                end
            else
                FileName = fullfile(s_fpath, scatter_file);
                storedStructure = load(FileName);
                variables = fields(storedStructure);
                
                options.ScatterC = storedStructure.(variables{1});
            end
        end
        if ~options.normalize_scatter
            if iscell(options.ScatterC)
                for llo = 1 : length(options.ScatterC)
                    if size(options.ScatterC{llo},2) ~= options.Nang && size(options.ScatterC{llo},1) ~= options.Ndist && numel(options.ScatterC{llo}) == numel(SinM)
                        options.ScatterC{llo} = reshape(options.ScatterC{llo}, options.Ndist,options.Nang,options.NSinos);
                    elseif size(options.ScatterC{llo},1)*size(options.ScatterC{llo},2)*size(options.ScatterC{llo},3) ~= numel(SinM)
                        error('Size mismatch between scatter sinogram and and specified sinogram')
                    end
                    if options.scatter_variance_reduction
                        options.ScatterC{llo} = Randoms_variance_reduction(double(options.ScatterC{llo}), options);
                        ScatterProp.variance_reduction = true;
                        appliedCorrections.scatter = 'scatter correction with variance reduction';
                    end
                    if options.scatter_smoothing
                        options.ScatterC{llo} = randoms_smoothing(options.ScatterC{llo}, options);
                        ScatterProp.smoothing = true;
                        if options.scatter_variance_reduction
                            appliedCorrections.scatter = 'and smoothing';
                        else
                            appliedCorrections.scatter = 'scatter correction with smoothing';
                        end
                    else
                        if ~options.scatter_variance_reduction
                            appliedCorrections.scatter = 'scatter correction';
                        end
                    end
                    
                    if options.partitions > 1
                        SinM{llo} = single(SinM{llo}) - single(options.ScatterC{llo});
                    else
                        SinM = single(SinM) - single(options.ScatterC{1});
                    end
                end
            else
                if size(options.ScatterC,2) ~= options.Nang && size(options.ScatterC,1) ~= options.Ndist && numel(options.ScatterC) == numel(SinM)
                    options.ScatterC = reshape(options.ScatterC, options.Ndist,options.Nang,options.NSlices);
                elseif size(options.ScatterC,1)*size(options.ScatterC,2)*size(options.ScatterC,3) ~= numel(SinM)
                    error('Size mismatch between scatter sinogram and and specified sinogram')
                end
                if options.scatter_variance_reduction
                    options.ScatterC = Randoms_variance_reduction(double(options.ScatterC), options);
                    ScatterProp.variance_reduction = true;
                    appliedCorrections.scatter = 'scatter correction with variance reduction';
                end
                if options.scatter_smoothing
                    options.ScatterC = randoms_smoothing(options.ScatterC, options);
                    ScatterProp.smoothing = true;
                    if options.scatter_variance_reduction
                        appliedCorrections.scatter = 'and smoothing';
                    else
                        appliedCorrections.scatter = 'scatter correction with smoothing';
                    end
                else
                    if ~options.scatter_variance_reduction
                        appliedCorrections.scatter = 'scatter correction';
                    end
                end
                if options.partitions > 1
                    for llo = 1 : options.partitions
                        SinM{llo} = single(SinM{llo}) - single(options.ScatterC);
                    end
                else
                    SinM = single(SinM) - single(options.ScatterC);
                end
            end
            if options.verbose
                disp('Scatter correction applied to sinogram')
            end
        end
    end
    
    if options.normalization_correction && ~options.corrections_during_reconstruction
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
            norm_file = [folder options.machine_name '_normalization_' num2str(options.Ndist) 'x' num2str(options.Nang) '_span' num2str(options.span) '.mat'];
            if exist(norm_file, 'file') == 2
                normalization = loadStructFromFile(norm_file,'normalization');
            else
                normalization = normalization_coefficients(options);
                %                 normalization = loadStructFromFile(norm_file,'normalization');
            end
        end
        if options.partitions > 1
            for llo = 1 : options.partitions
                SinM{llo} = single(SinM{llo}) .* normalization;
            end
        else
            SinM = single(SinM) .* normalization;
        end
        if options.verbose
            disp('Normalization correction applied to sinogram')
        end
        appliedCorrections.normalization = true;
    end
    if options.scatter_correction && options.normalize_scatter
        if iscell(options.ScatterC)
            if size(options.ScatterC{llo},2) ~= options.Nang && size(options.ScatterC{llo},1) ~= options.Ndist && numel(options.ScatterC{llo}) == numel(SinM)
                options.ScatterC{llo} = reshape(options.ScatterC{llo}, options.Ndist,options.Nang,options.NSlices);
            elseif size(options.ScatterC{llo},1)*size(options.ScatterC{llo},2)*size(options.ScatterC{llo},3) ~= options.Ndist*options.Nang*options.NSlices
                error('Size mismatch between scatter sinogram and specified sinogram')
            end
            if options.normalization_correction && options.normalize_scatter
                options.ScatterC{kk} = single(options.ScatterC{kk}) .* reshape(normalization, size(options.ScatterC{kk}, 1), size(options.ScatterC{kk}, 2), size(options.ScatterC{kk}, 3));
                ScatterProp.normalization = true;
            end
            if options.scatter_variance_reduction
                options.ScatterC{llo} = Randoms_variance_reduction(double(options.ScatterC{llo}), options);
                ScatterProp.variance_reduction = true;
                appliedCorrections.scatter = 'normalized scatter correction with variance reduction';
            end
            if options.scatter_smoothing
                options.ScatterC{llo} = randoms_smoothing(options.ScatterC{llo}, options);
                ScatterProp.smoothing = true;
                if options.scatter_variance_reduction
                    appliedCorrections.scatter = 'and smoothing';
                else
                    appliedCorrections.scatter = 'normalized scatter correction with smoothing';
                end
            else
                if ~options.scatter_variance_reduction
                    appliedCorrections.scatter = 'normalized scatter correction';
                end
            end
            SinM = SinM - single(options.ScatterC{llo});
        else
            if size(options.ScatterC,2) ~= options.Nang && size(options.ScatterC,1) ~= options.Ndist && numel(options.ScatterC) == numel(SinM)
                options.ScatterC = reshape(options.ScatterC, options.Ndist,options.Nang,options.NSlices);
            elseif size(options.ScatterC,1)*size(options.ScatterC,2)*size(options.ScatterC,3) ~= options.Ndist*options.Nang*options.NSlices
                error('Size mismatch between scatter sinogram and specified sinogram')
            end
            if options.normalization_correction && options.normalize_scatter
                options.ScatterC = single(options.ScatterC) .* reshape(normalization, size(options.ScatterC, 1), size(options.ScatterC, 2), size(options.ScatterC, 3));
                ScatterProp.normalization = true;
            end
            if options.scatter_variance_reduction
                options.ScatterC = Randoms_variance_reduction(double(options.ScatterC), options);
                ScatterProp.variance_reduction = true;
                appliedCorrections.scatter = 'normalized scatter correction with variance reduction';
            end
            if options.scatter_smoothing
                options.ScatterC = randoms_smoothing(options.ScatterC, options);
                ScatterProp.smoothing = true;
                if options.scatter_variance_reduction
                    appliedCorrections.scatter = 'and smoothing';
                else
                    appliedCorrections.scatter = 'normalized scatter correction with smoothing';
                end
            else
                if ~options.scatter_variance_reduction
                    appliedCorrections.scatter = 'normalized scatter correction';
                end
            end
            SinM = SinM - single(options.ScatterC);
        end
        if options.verbose
            disp('Scatter correction applied to sinogram')
        end
    end
    if isfield(options,'global_correction_factor') && ~isempty(options.global_correction_factor) && options.global_correction_factor > 0 && ~options.corrections_during_reconstruction
        SinM = SinM * options.global_correction_factor;
    end
    
    %%
    if options.verbose
        disp('Saving sinogram...')
    end
    variableList = {'SinM','SinSC','ScatterProp'};
    save_string = [options.machine_name '_' options.name '_sinogram_original_static_' num2str(options.Ndist) 'x' num2str(options.Nang)...
        'x' num2str(options.TotSinos) '_span' num2str(options.span) '_machine_sinogram.mat'];
    save(save_string,'raw_SinM','-v7')
    for variableIndex = 1:length(variableList)
        if exist(variableList{variableIndex},'var')
            if exist('OCTAVE_VERSION','builtin') == 0
                save(save_string,variableList{variableIndex},'-append')
            else
                save(save_string,variableList{variableIndex},'-append','-v7')
            end
        end
    end
    if exist('OCTAVE_VERSION','builtin') == 0
        save(save_string,'appliedCorrections','-append')
    else
        save(save_string,'appliedCorrections','-append','-v7')
    end
    if (options.scatter_correction || options.normalization_correction) && ~options.corrections_during_reconstruction
        varargout{1} = SinM;
        varargout{2} = raw_SinM;
    else
        varargout{1} = raw_SinM;
    end
    if options.verbose
        disp('Sinograms formed and saved')
        toc
    end
end