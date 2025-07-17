function [varargout] = normalization_coefficients(options)
%% Args
% options: Sinogram data, detector pair coordinates, ring heights, apparatus radius
% normalization_attenuation_correction: apply attenuation correction
% yes = [inner_cylinder_radius (cm) cylinder_attenuation_constant (cm^2/m)] | no = empty)
% If inner_cylinder_radius=inf, cylinder is assumed to cover entire FOV
% If left empty uniform illumination for each LOR is assumed
% options.normalization_scatter_correction: fit gaussian to scatter tail from cylinder
% normalization data (Too narrow scatter tail may lead to unaccurate fit).
% Not supported for list-mode data.
% options.normalization_options(1): apply axial geometric correction (yes = 1 | no = 0)
% options.normalization_options(2): apply detector effiency correction. (Fansum = 1 | SPC = 2 | no = 0)
% Fan_sum uses 3-D fansum method for both data types or SPC "single-plane
% Casey" method for list mode-data (SPC computationally more expensive). SPC is
% supposed to be used with FOV covering source
% Fansum version includes block profile correction. Using
% options.normalization_options(2)=2 with fansum uses block profile correction
% before detector effiency correction
% options.normalization_options(4): apply transaxial geometric correction for plane source data (yes = 1 | no = 0)
% With cylinders transaxial correction produces appropriate coeffs for LORs
% passing cylinder (LORs passing near cylinder edges can also be inaccurate)
% Accuracy of radial sorting can be adjusted in transaxial correction section
% options.normalization_options(3): apply block profile correction. If using fansum
% correction is nested with effiency correction
% If using SPC block profile correction is done separately
% TRANAXIAL AND BLOCK PROFILE CORRECTION GROUP SORTING CAN BE ADJUSTED IN
% THEIR SECTIONS
%% Returns
% Normalization matrix containing all coeffs
% Corrected sinogram
% Individual correction coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2020 Anssi Manninen, Ville-Veikko Wettenhovi
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


folder = fileparts(which('normalization_coefficients.m'));
folder = [folder(1:end-(6 + 8)), 'mat-files/'];
folder = strrep(folder, '\','/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


uniform_source=0;

%if using uniform source for instance infinidesimal thickness cylinder use
%uniform_source=1;



varargout = cell(nargout,1);

disp('Computing normalization coefficients')


%% Coordinates

if options.span == 1
    error('Normalization correction doesn''t yet work with span of 1')
end

cryst_per_block = options.cryst_per_block;
if ~options.use_raw_data
    Nang = options.Nang;
    mashing = 1;
    if options.det_w_pseudo / options.Nang > 2
        mashing = (options.det_w_pseudo / options.Nang / 2);
        %         Nang = options.Nang * mashing;
        %         options.Nang = Nang;
    end
    
    z = sinogram_coordinates_3D(options);
    z_length = options.rings * options.cr_pz;
    z_true = linspace(0, z_length, options.rings + 1)./10;
    if min(z_true(:)) == 0
        z_true = z_true + ((options.axial_fov - options.rings * options.cr_pz)/2 + options.cr_pz/2) / 10;
    end
    %mm --> cm
    z = single(z./10);
    z = z - min(z(:));
    
    [~, ~, xp, yp] = detector_coordinates(options);
    [x, y] = sinogram_coordinates_2D(options, xp, yp);
    
    %mm --> cm
    x = single(x./10);
    y = single(y./10);
    x = x - min(x(:));
    y = y - min(y(:));
    
    if options.det_w_pseudo > options.det_per_ring
        cryst_per_block = cryst_per_block + 1;
    end
    %     if mashing > 1
    %         options.Nang = options.Nang / mashing;
    %     end
    
    rings = options.rings;
else
    
    temp = options.pseudot;
    if ~isempty(temp) && sum(temp) > 0
        for kk = uint32(1) : temp
            pseudot(kk) = uint32(options.cryst_per_block + 1) * kk;
        end
    elseif temp == 0
        pseudot = [];
    end
    
    z_length = double(options.rings + 1 + sum(options.pseudot)) * options.cr_pz;
    z = linspace(0, z_length, options.rings + 2 + sum(options.pseudot))';
    if sum(pseudot) > 0
        z(pseudot) = [];
    end
    if min(z(:)) == 0
        z = z + (options.axial_fov - (options.rings + sum(options.pseudot)) * options.cr_pz)/2 + options.cr_pz/2;
    end
    
    
    z = single(z./10);
    
end

%x and y detector coordinates
if options.use_raw_data
    [detectors_x, detectors_y] = detector_coordinates(options);
    detectors_ring = single(options.det_per_ring);
else
    [~, ~, detectors_x, detectors_y] = detector_coordinates(options);
    detectors_ring = single(options.det_w_pseudo);
end

%mm --> cm
detectors_x = single(detectors_x./10);
detectors_y = single(detectors_y./10);
detectors_x = detectors_x - min(detectors_x(:));
detectors_y = detectors_y - min(detectors_y(:));


sino_amount = single(sum(options.segment_table));
segment_amount = single(length(options.segment_table));

%Ring radius

R = single(options.diameter/10/2); %cm

normalization_attenuation_correction = [options.normalization_phantom_radius options.normalization_attenuation];

%Inner ring radius


if ~isempty(normalization_attenuation_correction)
    
    r = normalization_attenuation_correction(1); %cm
    if length(normalization_attenuation_correction) == 2
        attenuation_coeff = normalization_attenuation_correction(2);  %mass attenuation times density
    end
    
else
    r = inf;
    if length(normalization_attenuation_correction) == 2
        normalization_attenuation_correction(2) = [];
    end
end


%% Scale stacked data (when using sinograms)

[GATE_vars, I] = sort([options.use_ASCII, options.use_root, options.use_LMF],'descend');
GATE_char = {'ASCII';'root';'LMF'};
GATE_char = GATE_char(I);

if options.use_raw_data
    if options.partitions == 1
        load_string = [options.machine_name '_measurements_' options.name '_static_raw'];
        if options.use_ASCII && options.use_machine == 0
            load_string =  [load_string '_ASCII.mat'];
        elseif options.use_LMF && options.use_machine == 0
            load_string =  [load_string '_LMF.mat'];
        elseif options.use_root && options.use_machine == 0
            load_string =  [load_string '_root.mat'];
        else
            load_string =  [load_string '_listmode.mat'];
        end
    else
        load_string = [options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_' ...
            num2str(options.tot_time) 's_raw'];
        load_string2 = [options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
            num2str(options.tot_time) 's_raw'];
        if options.use_ASCII && options.use_machine == 0
            if exist([load_string '_ASCII.mat'], 'file') == 0
                load_string = load_string2;
            end
            load_string =  [load_string '_ASCII.mat'];
        elseif options.use_LMF && options.use_machine == 0
            if exist([load_string '_LMF.mat'], 'file') == 0
                load_string = load_string2;
            end
            load_string =  [load_string '_LMF.mat'];
        elseif options.use_root && options.use_machine == 0
            if exist([load_string '_root.mat'], 'file') == 0
                load_string = load_string2;
            end
            load_string =  [load_string '_root.mat'];
        else
            if exist([load_string '_listmode.mat'], 'file') == 0
                load_string = load_string2;
            end
            load_string =  [load_string '_listmode.mat'];
        end
    end
    if (isfield(options, 'coincidences') == 0 && ~exist('coincidences','var')) && options.use_machine < 2
        if options.partitions == 1
            if any(GATE_vars) && options.use_machine == 0
                options.coincidences = loadGATEData(load_string,'coincidences', GATE_char);
            else
                load(load_string, 'coincidences')
                options.coincidences = coincidences;
                clear coincidences
            end
        else
            if any(GATE_vars) && options.use_machine == 0
                options.coincidences = loadGATEData({load_string;load_string2},'coincidences', GATE_char);
            else
                if exist(load_string, 'file') == 0
                    load_string = load_string2;
                end
                load(load_string, 'coincidences')
                options.coincidences = coincidences;
                clear coincidences
            end
        end
    end
    % Sinogram data
else
    if (~options.reconstruct_trues && ~options.reconstruct_scatter) || options.use_machine > 0
        if options.partitions == 1
            if options.span > 1
                load_string = [options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                    num2str(options.TotSinos) '_span' num2str(options.span)];
            else
                load_string = [options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                    num2str(options.rings^2) '_span' num2str(options.span)];
            end
            if options.use_machine == 0
                sinoFile = [load_string '.mat'];
            elseif options.use_machine == 1
                sinoFile = [load_string '_listmode.mat'];
            elseif options.use_machine == 2
                sinoFile = [load_string '_machine_sinogram.mat'];
            elseif options.use_machine == 3
                sinoFile = [load_string '_listmode_sinogram.mat'];
            end
        else
            load_string = [options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_' ...
                num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                num2str(options.span)];
            load_string2 = [options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                num2str(options.span)];
            if options.use_machine == 0
                sinoFile = [load_string '.mat'];
                if exist(sinoFile, 'file') == 0
                    sinoFile = [load_string2 '.mat'];
                end
            elseif options.use_machine == 1
                sinoFile = [load_string '_listmode.mat'];
                if exist(sinoFile, 'file') == 0
                    sinoFile = [load_string2 '_listmode.mat'];
                end
            elseif options.use_machine == 2
                sinoFile = [load_string '_machine_sinogram.mat'];
                if exist(sinoFile, 'file') == 0
                    sinoFile = [load_string2 '_machine_sinogram.mat'];
                end
            elseif options.use_machine == 3
                sinoFile = [load_string '_listmode_sinogram.mat'];
                if exist(sinoFile, 'file') == 0
                    sinoFile = [load_string2 '_listmode_sinogram.mat'];
                end
            end
        end
        if options.partitions == 1 && (isfield(options, 'SinM') == 0 || (isfield(options, 'SinM') == 1 && numel(options.SinM) <= 1))
            if options.use_machine < 2
                if options.use_machine == 0
                    try
                        load(sinoFile,'raw_SinM')
                    catch ME
                        if mashing > 1
                            error('When computing normalization coefficients for mashed data, the input sinogram has to be formed with mashing of 1!')
                        else
                            error(ME)
                        end
                    end
                elseif  options.use_machine == 1
                    try
                        load(sinoFile,'raw_SinM')
                    catch ME
                        if mashing > 1
                            error('When computing normalization coefficients for mashed data, the input sinogram has to be formed with mashing of 1!')
                        else
                            error(ME)
                        end
                    end
                end
                options.SinM = raw_SinM;
                clear raw_SinM
            else
                try
                    load(sinoFile,'SinM')
                catch ME
                    if mashing > 1
                        error('When computing normalization coefficients for mashed data, the input sinogram has to be formed with mashing of 1!')
                    else
                        error(ME)
                    end
                end
                options.SinM = SinM;
                clear SinM
            end
        elseif isfield(options, 'SinM') == 0 || (isfield(options, 'SinM') == 1 && numel(options.SinM) <= 1)
            if options.use_machine < 2
                if options.use_machine == 0
                    try
                        load(sinoFile, 'raw_SinM')
                    catch ME
                        if mashing > 1
                            error('When computing normalization coefficients for mashed data, the input sinogram has to be formed with mashing of 1!')
                        else
                            error(ME)
                        end
                    end
                elseif  options.use_machine == 1
                    try
                        load(sinoFile, 'raw_SinM')
                    catch ME
                        if mashing > 1
                            error('When computing normalization coefficients for mashed data, the input sinogram has to be formed with mashing of 1!')
                        else
                            error(ME)
                        end
                    end
                end
                options.SinM = raw_SinM;
                clear raw_SinM
            else
                try
                    load(sinoFile, 'SinM')
                catch ME
                    if mashing > 1
                        error('When computing normalization coefficients for mashed data, the input sinogram has to be formed with mashing of 1!')
                    else
                        error(ME)
                    end
                end
                options.SinM = SinM;
                clear SinM
            end
        end
    end
end


if ~options.use_raw_data
    
    tic
    
    %Convert coincidence data
    
    if isstruct(options.SinM)
        options.SinM = cell2mat(struct2cell(options.SinM));
    end
    SinDouble = single(options.SinM);
    
    %Initialize normalization matrix
    
    normalization = ones(size(options.SinM), 'single');
    
    if mashing > 1
        %         SinDouble(SinDouble == 0) = NaN;
        %         apu = cell2mat(arrayfun(@(i) mean(SinDouble(:,i:i+mashing-1,:),2),1:mashing:size(SinDouble,2)-mashing+1,'UniformOutput',false));
        %         gaps = apu == 0;
    end
    
    %Amount of stacked planes in even and odd numbered plane sinograms
    
    %     if rem(floor(options.span/2),2)==0
    %
    %         odds=ceil(options.span/2);
    %         evens=floor(options.span/2);
    %
    %     else
    %
    %         odds=floor(options.span/2);
    %         evens=ceil(options.span/2);
    %
    %     end
    %     last_index=0;
    
    %     for t=1:segment_amount
    %
    %         if t==1
    %             for i=2:options.segment_table(1)-1
    %
    %                 %Michelogram corners
    %
    %                 if i<ceil(options.span/2) || options.segment_table(1)-i<ceil(options.span/2)
    %
    %                     SinDouble(:,:,i)=SinDouble(:,:,i)/min(i,options.segment_table(1)-i+1);
    %                     normalization(:,:,i)=normalization(:,:,i)/min(i,options.segment_table(1)-i+1);
    %
    %                     %Michelogram mid-section
    %                 else
    %
    %                     if rem(i,2)==0
    %                         SinDouble(:,:,i)=SinDouble(:,:,i)/evens;
    %                         normalization(:,:,i)=normalization(:,:,i)/evens;
    %                     else
    %                         SinDouble(:,:,i)=SinDouble(:,:,i)/odds;
    %                         normalization(:,:,i)=normalization(:,:,i)/odds;
    %                     end
    %                 end
    %             end
    %
    %             last_index=options.segment_table(1);
    %
    %         else
    %
    %             for i=last_index+3:last_index+options.segment_table(t)-2
    %
    %                 %Michelogram corners
    %                 if i-last_index>min(evens,odds)*2 && (last_index+options.segment_table(t))-i+1>min(evens,odds)*2
    %
    %
    %                     if rem(i-last_index,2)==0
    %                         SinDouble(:,:,i)=SinDouble(:,:,i)/min(odds,evens);
    %                         normalization(:,:,i)=normalization(:,:,i)/min(odds,evens);
    %
    %                     else
    %                         SinDouble(:,:,i)=SinDouble(:,:,i)/max(odds,evens);
    %                         normalization(:,:,i)=normalization(:,:,i)/max(odds,evens);
    %                     end
    %
    %                 else
    %                     %Michelogram mid-section
    %
    %                     SinDouble(:,:,i)=SinDouble(:,:,i)/min(ceil((i-last_index)/2),(ceil((last_index+options.segment_table(t)-i+1)/2)));
    %                     normalization(:,:,i)=normalization(:,:,i)/min(ceil(i/2),(ceil(last_index+options.segment_table(t)-i+1/2)));
    %                 end
    %
    %             end
    %
    %             last_index=last_index+options.segment_table(t);
    %
    %         end
    %     end
    
    %     time=toc;
    %     if options.verbose
    %         disp(['Spanned data scaled (', num2str(time),'s)'])
    %     end
    
end


%% Extended sparse list-mode matrix (when using list-mode data)

if options.use_raw_data
    
    tic
    
    %Convert coincidences to matrix from
    if sum(options.pseudot) > 0 || options.det_w_pseudo > options.det_per_ring
        koko = rings * options.det_per_ring;
    else
        koko = options.detectors;
    end
    
    true_coincidences = zeros(koko, koko,'single');
    if isa(options.coincidences,'cell')
        true_coincidences(tril(true(size(true_coincidences)), 0)) = single(full(options.coincidences{1}));
    else
        true_coincidences(tril(true(size(true_coincidences)), 0)) = single(full(options.coincidences));
    end
    
    
    x0 = max(detectors_x)/2;
    y0 = max(detectors_y)/2;
    [X,Y] = meshgrid(detectors_x, detectors_y);
    Y = Y';
    
    distance = ((abs((Y-detectors_y)*x0 - (X - detectors_x)*y0 + X.*detectors_y - Y.*detectors_x)./sqrt((Y-detectors_y).^2 + (X-detectors_x).^2)));
    
    if r == inf
        start_ind = find(distance(:,1) <= (options.FOVa_x/2/10),1,'first');
        end_ind = find(distance(:,1) <= (options.FOVa_x/2/10),1,'last');
    else
        start_ind = find(distance(:,1) <= r,1,'first');
        end_ind = find(distance(:,1) <= r,1,'last');
    end
    
    true_det_indices = false(detectors_ring);
    for kk = start_ind : end_ind
        true_det_indices(kk : detectors_ring + 1: detectors_ring * (detectors_ring - kk + 1)) = true;
    end
    true_det_apu = true_det_indices' + true_det_indices;
    true_det_indices = [true_det_indices; repmat(true_det_apu, rings-1,1)];
    true_det_indices2 = false(size(true_coincidences));
    for kk = 1 : rings
        true_det_indices2(1 + detectors_ring * (kk - 1) : end, 1 + detectors_ring * (kk - 1) : detectors_ring * kk) = ...
            true_det_indices(1 + detectors_ring * (kk - 1) : end, :);
    end
    true_det_indices=false(size(true_det_indices2));
    true_det_indices(tril(true(size(true_det_indices)), 0)) = true_det_indices2(tril(true(size(true_det_indices2)), 0));
    clear true_det_indices2 true_det_apu
    
    true_coincidences(~true_det_indices) = 0;
    
    %Initialize normalization matrix
    
    normalization = ones(size(true_coincidences),'single');
    
    
    %Get starting indices of FOV detector pairs for each column (for speed up)
    %first one of mirror projections (ring diff=0 have only 1 projection)
    
    start_ind_low_col=zeros(1,detectors_ring, 'single');
    end_ind_low_col=start_ind_low_col;
    for j=1:detectors_ring
        
        i=1;
        p=detectors_ring;
        out_of=1;
        out_of2=1;
        
        while true_det_indices(i,j)==0 && out_of
            
            i=i+1;
            
            if i>detectors_ring
                start_ind_low_col(j)=0;
                i=1;
                out_of=0;
            end
        end
        
        if out_of==1
            
            start_ind_low_col(j)=i;
            
        end
        
        while true_det_indices(p,j)==0 && out_of2
            
            p=p-1;
            
            if p==0
                end_ind_low_col(j)=0;
                p=1;
                out_of2=0;
            end
        end
        
        if out_of2==1
            
            end_ind_low_col(j)=p;
            
        end
        
    end
    
    start_ind_up_col=zeros(1,detectors_ring, 'single');
    end_ind_up_col=start_ind_up_col;
    
    for j=1:detectors_ring
        
        i=1;
        if start_ind_low_col(j)~=0
            p=start_ind_low_col(j)-1;
        else
            p=detectors_ring;
        end
        
        out_of=1;
        out_of2=1;
        
        while true_det_indices(i+detectors_ring,j)==0 && out_of
            
            i=i+1;
            
            if (i>=start_ind_low_col(j) && start_ind_low_col(j)~=0) || i==detectors_ring
                start_ind_up_col(j)=0;
                i=1;
                out_of=0;
            end
            
        end
        
        if out_of==1
            
            start_ind_up_col(j)=i;
            
        end
        
        while true_det_indices(p+detectors_ring,j)==0 && out_of2
            
            p=p-1;
            
            if p==0
                end_ind_up_col(j)=0;
                p=1;
                out_of2=0;
            end
        end
        
        if out_of2==1
            
            end_ind_up_col(j)=p;
            
        end
        
    end
    
    
    if options.normalization_options(2)==2
        
        %Get starting indices of FOV detector pairs for each row
        %first one of mirror projections (ring diff=0 have only 1 projection)
        
        start_ind_left_row=zeros(1,detectors_ring, 'single');
        end_ind_left_row=start_ind_left_row;
        
        for i=1:detectors_ring
            
            j=1;
            p=detectors_ring;
            out_of=1;
            out_of2=1;
            
            while true_det_indices(i,j)==0 && out_of
                
                j=j+1;
                
                if j>detectors_ring
                    start_ind_left_row(i)=0;
                    j=1;
                    out_of=0;
                end
            end
            
            if out_of==1
                
                start_ind_left_row(i)=j;
                
            end
            
            while true_det_indices(i,p)==0 && out_of2
                
                p=p-1;
                
                if p==0
                    end_ind_left_row(i)=0;
                    p=1;
                    out_of2=0;
                end
            end
            
            if out_of2==1
                
                end_ind_left_row(i)=p;
                
            end
            
        end
        
        start_ind_right_row=zeros(1,detectors_ring, 'single');
        end_ind_right_row=start_ind_right_row;
        
        for j=1:detectors_ring
            
            
            if end_ind_left_row(j)~=0
                i=end_ind_left_row(j)+1;
            else
                i=1;
            end
            
            p=detectors_ring;
            
            out_of=1;
            out_of2=1;
            
            while true_det_indices(j+detectors_ring,i)==0 && out_of
                
                i=i+1;
                
                if  i>detectors_ring
                    start_ind_right_row(j)=0;
                    i=1;
                    out_of=0;
                end
                
            end
            
            if out_of==1
                
                start_ind_right_row(j)=i;
                
            end
            
            if start_ind_right_row(j)~=0
                
                while true_det_indices(j+detectors_ring,p)==0 && out_of2
                    
                    p=p-1;
                    
                    if p==start_ind_right_row(j)
                        end_ind_right_row(j)=0;
                        p=1;
                        out_of2=0;
                    end
                end
                
                if out_of2==1
                    
                    end_ind_right_row(j)=p;
                    
                end
                
                
            end
            
        end
        
    end
    
    time=toc;
    if options.verbose
        disp(['Sparse matrices extended and false detector pairs disposed (', num2str(time),'s)'])
    end
else
    koko = options.detectors;
    x0 = max(x(:))/2;
    y0 = max(y(:))/2;
    
    distance = ((abs((y(:,1)-y(:,2))*x0 - (x(:,1) - x(:,2))*y0 + x(:,1).*y(:,2) - y(:,1).*x(:,2))./sqrt((y(:,1)-y(:,2)).^2 + (x(:,1)-x(:,2)).^2)));
    distance = reshape(distance, options.Ndist, Nang);
    
    if r == inf
        start_inds = find(distance(:,1) <= (options.FOVa_x/1.75/10),1,'first');
        end_inds = find(distance(:,1) <= (options.FOVa_x/1.75/10),1,'last');
    else
        start_inds = find(distance(:,1) <= r,1,'first');
        end_inds = find(distance(:,1) <= r,1,'last');
    end
end






%% Attenuation correction & scatter correction


if options.use_raw_data==0
    
    true_radius = sqrt((x(:,1)-x(:,2)).^2+(y(:,1)-y(:,2)).^2);
    true_radius = reshape(true_radius, options.Ndist, Nang);
    
    radius_in_cylinder=true_radius;
    
else
    
    true_radius=ones(detectors_ring,'single');
    
    radius_in_cylinder=true_radius;
    
end

if ~isempty(normalization_attenuation_correction)
    
    if r~=inf && ~options.use_raw_data
        
        %blur specs
        
        sigm = 2.0;
        %Window size
        sz = 4;
        x_gauss=-sz:sz;
        
        Exp_comp = -(x_gauss.^2)/(2*sigm^2);
        Kernel= exp(Exp_comp)/(sqrt(2*pi*sigm^2));
        
        if r~=inf
            
            radius_in_cylinder=zeros(1,Nang*options.Ndist, 'single');
            a=zeros(1,Nang*options.Ndist, 'single');
            c=a;
            intersects=zeros(Nang*options.Ndist,2, 'single');
            
            
            for i=1:Nang*options.Ndist
                
                %each LOR detector coordinates
                x_start=x(i,1);
                x_end=x(i,2);
                y_start=y(i,1);
                y_end=y(i,2);
                
                if (x_end-x_start)~=0
                    a(i)=(y_end-y_start)/(x_end-x_start);
                    if x_start<x_end
                        c(i)=y_start-a(i)*x_start;
                    else
                        c(i)=y_end-a(i)*x_end;
                    end
                else
                    a(i)=NaN;
                    c(i)=NaN;
                end
                
                %intersection points with inner cylinder
                if (x_end-x_start)~=0
                    intersects(i,:)=roots([(a(i)^2+1) (2*a(i)*c(i)-2*R-2*R*a(i)) (c(i)^2+2*R^2-2*R*c(i)-r^2)]);
                    
                    if imag(intersects(i,1))<0.0001 && imag(intersects(i,2))<0.0001
                        radius_in_cylinder(i)=sqrt((intersects(i,1)-intersects(i,2))^2+((intersects(i,1)-intersects(i,2))*a(i))^2);
                    else
                        radius_in_cylinder(i)=NaN;
                    end
                    
                else
                    intersects(i,:)=roots([1 -2*R (x_start-R)^2-r^2+R^2]);
                    if imag(intersects(i,:))<0.0001
                        radius_in_cylinder(i)=(intersects(i,1)-intersects(i,2));
                    else
                        radius_in_cylinder(i)=NaN;
                    end
                    
                end
                
            end
            true_radius=reshape(radius_in_cylinder,options.Ndist,Nang);
            
        end
        % Out of cylinder indices
        
        cut_off_start=zeros(1,Nang, 'single');
        cut_off_end=cut_off_start;
        
        for j=1:Nang
            i=1;
            while i<=options.Ndist/2 && isnan(true_radius(i,j))
                
                i=i+1;
                
            end
            
            cut_off_start(j)=i;
            
        end
        cut_off_start = cut_off_start(j);
        
        
        for j=1:Nang
            i=options.Ndist/2;
            while i<=options.Ndist && ~isnan(true_radius(i,j))
                
                i=i+1;
                
            end
            cut_off_end(j)=i-1;
        end
        cut_off_end = cut_off_end(j);
        
    elseif r~=inf && options.use_raw_data
        intersects=zeros(detectors_ring,detectors_ring,2, 'single');
        a=zeros(detectors_ring, 'single');
        c=a;
        radius_in_cylinder=c;
        
        for j=1:detectors_ring
            
            if (start_ind_low_col(j)~=0 && end_ind_low_col(j)~=0) || (start_ind_up_col(j)~=0 && end_ind_up_col(j)~=0)
                
                for i=[start_ind_up_col(j):end_ind_up_col(j) start_ind_low_col(j):end_ind_low_col(j)]
                    
                    if i~=0
                        
                        %each LOR detector coordinates
                        x_start=detectors_x(j);
                        x_end=detectors_x(i);
                        y_start=detectors_y(j);
                        y_end=detectors_y(i);
                        
                        if (x_end-x_start)~=0
                            a(i,j)=(y_end-y_start)/(x_end-x_start);
                            if x_start<x_end
                                c(i,j)=y_start-a(i,j)*x_start;
                            else
                                c(i,j)=y_end-a(i,j)*x_end;
                            end
                        else
                            a(i,j)=NaN;
                            c(i,j)=NaN;
                        end
                        
                        %intersection points with inner cylinder
                        if (x_end-x_start)~=0
                            intersects(i,j,:)=roots([(a(i,j)^2+1) (2*a(i,j)*c(i,j)-2*R-2*R*a(i,j)) (c(i,j)^2+2*R^2-2*R*c(i,j)-r^2)]);
                            
                            if imag(intersects(i,j,1))<0.0001 && imag(intersects(i,j,2))<0.0001
                                radius_in_cylinder(i,j)=sqrt((intersects(i,j,1)-intersects(i,j,2))^2+((intersects(i,j,1)-intersects(i,j,2))*a(i,j))^2);
                            else
                                radius_in_cylinder(i,j)=NaN;
                            end
                            
                        else
                            intersects(i,j,:)=roots([1 -2*R (x_start-R)^2-r^2+R^2]);
                            if imag(intersects(i,j,1))<0.0001 && imag(intersects(i,j,2))<0.0001
                                radius_in_cylinder(i,j)=(intersects(i,j,1)-intersects(i,j,2));
                            else
                                radius_in_cylinder(i,j)=NaN;
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
    end
    
    if (options.normalization_options(4)==1 || options.normalization_options(2)==2) && options.use_raw_data
        
        %Min length in cylinder for each ring difference
        
        rad_min=min(min(nonzeros(radius_in_cylinder)));
        
        true_radius=radius_in_cylinder;
        
        rad_profile_axial=zeros(detectors_ring,detectors_ring,options.ring_difference+1, 'single');
        h=zeros(1,options.ring_difference+1, 'single');
        activity_coeffs=zeros(size(true_coincidences), 'single');
        
        for u=0:rings - 1
            
            %Height from first ring of segment u
            
            h(u+1)=abs(z(1)-z(1+u));
            
            if h(u+1)~=0
                
                rad_profile_axial(:,:,u+1)=sqrt(radius_in_cylinder.^2+((radius_in_cylinder~=0).*(~isnan(radius_in_cylinder)).*h(u+1)).^2);
                
            else
                
                rad_profile_axial(:,:,u+1)=radius_in_cylinder;
                
            end
            
            rad_profile_axial(isnan(rad_profile_axial))=rad_min;
            rad_profile_axial(rad_profile_axial==0)=rad_min;
            
            if r~=inf && (options.normalization_options(4)==1 || options.normalization_options(2)==2) && uniform_source==0
                
                for q=1:rings-u
                    
                    activity_coeffs((u+q-1)*detectors_ring+1:(u+q)*detectors_ring,(q-1)...
                        *detectors_ring+1:q*detectors_ring)=max(max(rad_profile_axial(:,:,u+1)))./rad_profile_axial(:,:,u+1);
                    
                end
            end
        end
    end
    
end

if (options.normalization_scatter_correction || ~isempty(normalization_attenuation_correction)) && ~options.use_raw_data
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Scatter correction for sinogram%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if options.normalization_scatter_correction
        
        if r~=inf
            
            tic
            
            % Transaxial scatter counts
            
            Scatters=zeros(1,options.Ndist, 'single');
            for j=1:options.Ndist
                
                if sum(isnan(true_radius(j,:)))==Nang
                    
                    Scatters(j)=sum(sum(SinDouble(j,:,:)));
                    
                end
                
            end
            
            Sino_counts=zeros(1,sino_amount, 'single');
            
            cut_pix=ceil(options.Ndist/100)+1;
            Scatters=[Scatters(1:cut_off_start-cut_pix) Scatters(cut_off_end+cut_pix:end)]'./sino_amount;
            
            %fit gaussian to scatter tail (closest transaxial sums to cylinder are cut off)
            
            sigm=zeros(1,1000, 'single');
            const=sigm;
            [~ , mu]=max(nanmax(true_radius));
            xx=[1:cut_off_start-cut_pix cut_off_end+cut_pix:options.Ndist]';
            [const(1), max_ind]=max(Scatters);
            sigm(1)=sqrt((-(xx(max_ind)-mu)^2)/(2*log(1/3)));
            yy=1:options.Ndist;
            i=1;
            old_coeffs=[0 0];
            new_coeffs=old_coeffs;
            
            while (norm(new_coeffs-old_coeffs)>0.0001 || i==1) && i<1000
                
                old_coeffs=[sigm(i) ; const(i)];
                
                J_h(:,1)=const(i).*(xx-mu*ones(size(xx),'single')).^2/(sigm(i)^3).*exp(-(xx-mu*ones(size(xx),'single')).^2./(2*sigm(i)^2));
                
                J_h(:,2)=exp(-(xx-mu*ones(size(xx),'single')).^2./(2*sigm(i)^2));
                
                h_i=const(i).*exp(-(xx-mu*ones(size(xx),'single')).^2./(2*sigm(i)^2));
                
                new_coeffs=Gauss_Newton(Scatters,h_i,J_h,[sigm(i) ; const(i)],1);
                sigm(i+1)=new_coeffs(1);
                const(i+1)=new_coeffs(2);
                
                i=i+1;
                
            end
            
            values=const(i).*exp(-(yy'-mu*ones(size(yy),'single')').^2./(2*sigm(i)^2));
            
            if i==1000
                
                error('Gaussian fit for scatter not converging')
                
            end
            
            for u=1:sino_amount
                
                
                if sum(isnan(true_radius(j,:)))==Nang
                    
                    Sino_counts(u)=sum(sum(SinDouble(1:cut_off_start-cut_pix,:,u)))+sum(sum(SinDouble(cut_off_end+cut_pix:end,:,u)));
                    
                end
                
            end
            
            Scaling_factors=Sino_counts./mean(Sino_counts);
            
            
            %export scattered counts
            Scatter_matrix=zeros(size(SinDouble), 'single');
            for u=1:sino_amount
                
                Sino_gauss=Scaling_factors(u).*values'./Nang;
                
                Scatter_matrix(:,:,u)=repmat(Sino_gauss', 1, Nang);
                
            end
            SinDouble = SinDouble-Scatter_matrix;
            
            SinDouble(SinDouble<0)=0;
            
            time=toc;
            if options.verbose
                disp(['Scatter correction done (', num2str(time),'s)'])
            end
            
        else
            
            warning('Scatter correction is not supported for non-cylinder source')
            
        end
        
    else
        
        disp('Normalization coefficients are calculated without scatter correction')
        
    end
    
    
    if ~isempty(normalization_attenuation_correction) && length(normalization_attenuation_correction)==2
        
        tic
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  Attenuation for sinogram data  %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        rad_coeff_matrix = zeros(size(SinDouble), 'single');
        
        C = [x(:,1) - options.diameter/2/10 y(:,1) - options.diameter/2/10 zeros(size(x,1),1, 'single')]';
        D = [x(:,2) - options.diameter/2/10 y(:,2) - options.diameter/2/10 zeros(size(x,1),1, 'single')]';
        t = (dot(C, (C-D)))./(dot((C - D), (C - D)));
        d = dot(C,C) - t.^2.*dot((C - D),(C-D));
        d(d < 0) = 0;
        d = sqrt(d);
        d(d > r) = NaN;
        k = sqrt((r^2 - d.^2)./(dot(C - D, C - D)));
        
        for u=1:segment_amount
            
            
            %height from first ring of segment u
            if u==1
                ring=1;
            else
                ring=sum(options.segment_table(1:u-1))+1;
            end
            
            
            
            if u~=1
                start_ind=sum(options.segment_table(1:(u-1)))+1;
            else
                start_ind=1;
            end
            end_ind=sum(options.segment_table(1:u));
            
            
            C2 = [x(:,1) - options.diameter/2/10 y(:,1) - options.diameter/2/10 repmat(z(ring,1), 1, size(x,1))'];
            D2 = [x(:,2) - options.diameter/2/10 y(:,2) - options.diameter/2/10 repmat(z(ring,2), 1, size(x,1))'];
            P1 = C2 + (t + k)' .* (D2 - C2);
            P2 = C2 + (t - k)' .* (D2 - C2);
            rad_profile_axial = sqrt(sum(abs((P1' - P2')).^2,1));
            rad_profile_axial(isnan(rad_profile_axial)) = 0;
            
            rad_coeffs=exp(-rad_profile_axial * attenuation_coeff);
            
            %Blur attenuation coefficient sharp altering
            rad_coeffs=reshape(rad_coeffs,options.Ndist,Nang);
            I=[rad_coeffs; (ones(sz+1,Nang,'single').*rad_coeffs(1:sz+1,:))];
            I=I(:);
            blurred=filter(Kernel,1,I');
            %                 blurred = I;
            
            new_coeffs=reshape(blurred',options.Ndist+sz+1,Nang);
            new_coeffs=new_coeffs(sz+2:end,:);
            
            
            rad_coeff_matrix(:,:,start_ind:end_ind)=repmat(new_coeffs,1,1,options.segment_table(u));
        end
        
        
        time=toc;
        if options.verbose
            disp(['Attenuation correction done (', num2str(time),'s)'])
        end
    end
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Attenuation correction for  list-mode data %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(normalization_attenuation_correction) && options.use_raw_data && length(normalization_attenuation_correction)==2
    
    tic
    attenuation_coeff=normalization_attenuation_correction(2);  %mass attenuation times density
    %blur specs
    
    sigm = 2.0;
    %Window size
    sz = 4;
    x_gauss=-sz:sz;
    
    Exp_comp = -(x_gauss.^2)/(2*sigm^2);
    Kernel= exp(Exp_comp)/(sqrt(2*pi*sigm^2));
    
    [X,Y] = meshgrid(detectors_x, detectors_y);
    XX = X';
    YY = Y';
    
    C = [X(:) - options.diameter/2/10 YY(:) - options.diameter/2/10 zeros(size(X(:),1),1, 'single')]';
    D = [XX(:) - options.diameter/2/10 Y(:) - options.diameter/2/10 zeros(size(XX(:),1),1, 'single')]';
    t = (dot(C, (C-D)))./(dot((C - D), (C - D)));
    d = dot(C,C) - t.^2.*dot((C - D),(C-D));
    d(d < 0) = 0;
    d = sqrt(d);
    d(d > r) = NaN;
    k = sqrt((r^2 - d.^2)./(dot(C - D, C - D)));
    
    if r~=inf
        C2 = [X(:) - options.diameter/2/10 YY(:) - options.diameter/2/10 repmat(z(1), 1, size(X(:),1))'];
        D2 = [XX(:) - options.diameter/2/10 Y(:) - options.diameter/2/10 repmat(z(1), 1, size(X(:),1))'];
        P1 = C2 + (t + k)' .* (D2 - C2);
        P2 = C2 + (t - k)' .* (D2 - C2);
        radius_in_cylinder = sqrt(sum(abs((P1' - P2')).^2,1));
    end
    
    true_radius = reshape(radius_in_cylinder, detectors_ring, detectors_ring);
    
    activity_coeffs=ones(size(true_coincidences),'single');
    
    for u=0:rings-1
        
        C2 = [X(:) - options.diameter/2/10 YY(:) - options.diameter/2/10 repmat(z(1), 1, size(X(:),1))'];
        D2 = [XX(:) - options.diameter/2/10 Y(:) - options.diameter/2/10 repmat(z(u+1), 1, size(X(:),1))'];
        P1 = C2 + (t + k)' .* (D2 - C2);
        P2 = C2 + (t - k)' .* (D2 - C2);
        rad_profile_axial = sqrt(sum(abs((P1' - P2')).^2,1));
        
        rad_profile_axial(isnan(rad_profile_axial))=0;
        
        rad_coeffs=reshape(exp(-rad_profile_axial.*attenuation_coeff), detectors_ring, detectors_ring);
        
        
        %blur coeffs with gaussian window
        temp_vec=rad_coeffs;
        [sorted,indices]=sort(temp_vec(:));
        I=[sorted(1).*ones(sz*2+1,1, 'single') ; sorted];
        
        blurred=filter(Kernel,1,I');
        
        new_coeffs=blurred(sz*2+2:end);
        
        rad_coeffs(indices)=new_coeffs;
        
        %Blur radial coeffs with
        
        for q=1:rings-u
            
            activity_coeffs((u+q-1)*detectors_ring+1:(u+q)*detectors_ring,(q-1)*detectors_ring+1:q*detectors_ring)=rad_coeffs;
            
        end
        
        
    end
    
    time=toc;
    if options.verbose
        disp(['Attenuation correction done (', num2str(time),'s)'])
    end
end

if options.use_raw_data
    
    if options.normalization_options(2)==2
        
        true_radius(isnan(true_radius))=0;
        lower_radius=tril(true_radius);
        
        start_ind_low_col_cylinder=zeros(1,detectors_ring, 'single');
        end_ind_low_col_cylinder=start_ind_low_col_cylinder;
        
        for j=1:detectors_ring
            
            i=1;
            p=detectors_ring;
            out_of=1;
            out_of2=1;
            
            while lower_radius(i,j)==0 && out_of
                
                i=i+1;
                
                if i>detectors_ring
                    start_ind_low_col_cylinder(j)=0;
                    i=1;
                    out_of=0;
                end
            end
            
            if out_of==1
                
                start_ind_low_col_cylinder(j)=i;
                
            end
            
            while lower_radius(p,j)==0 && out_of2
                
                p=p-1;
                
                if p==0
                    end_ind_low_col_cylinder(j)=0;
                    p=1;
                    out_of2=0;
                end
            end
            
            if out_of2==1
                
                end_ind_low_col_cylinder(j)=p;
                
            end
            
        end
        
        start_ind_up_col_cylinder=zeros(1,detectors_ring, 'single');
        end_ind_up_col_cylinder=start_ind_up_col_cylinder;
        
        for j=1:detectors_ring
            
            i=1;
            if start_ind_low_col_cylinder(j)~=0
                p=start_ind_low_col_cylinder(j)-1;
            else
                p=detectors_ring;
            end
            
            out_of=1;
            out_of2=1;
            
            while true_radius(i,j)==0 && out_of
                
                i=i+1;
                
                if (i>=start_ind_low_col_cylinder(j) && start_ind_low_col_cylinder(j)~=0) || i==detectors_ring
                    start_ind_up_col_cylinder(j)=0;
                    i=1;
                    out_of=0;
                end
                
            end
            
            if out_of==1
                
                start_ind_up_col_cylinder(j)=i;
                
            end
            
            while true_radius(p,j)==0 && out_of2
                
                p=p-1;
                
                if p==0
                    end_ind_up_col_cylinder(j)=0;
                    p=1;
                    out_of2=0;
                end
                
            end
            
            if out_of2==1
                
                end_ind_up_col_cylinder(j)=p;
                
            end
            
        end
        
        
        %Get starting indices of FOV detector pairs for each row (for speed up)
        %first one of mirror projections (ring diff=0 have only 1 projection)
        start_ind_left_row_cylinder=zeros(1,detectors_ring, 'single');
        end_ind_left_row_cylinder=start_ind_left_row_cylinder;
        
        for i=1:detectors_ring
            
            j=1;
            p=detectors_ring;
            out_of=1;
            out_of2=1;
            
            while lower_radius(i,j)==0 && out_of
                
                j=j+1;
                
                if j>detectors_ring
                    start_ind_left_row_cylinder(i)=0;
                    j=1;
                    out_of=0;
                end
            end
            
            if out_of==1
                
                start_ind_left_row_cylinder(i)=j;
                
            end
            
            while lower_radius(i,p)==0 && out_of2
                
                p=p-1;
                
                if p==0
                    end_ind_left_row_cylinder(i)=0;
                    p=1;
                    out_of2=0;
                end
            end
            
            if out_of2==1
                
                end_ind_left_row_cylinder(i)=p;
                
            end
            
        end
        
        start_ind_right_row_cylinder=zeros(1,detectors_ring, 'single');
        end_ind_right_row_cylinder=start_ind_right_row_cylinder;
        
        for j=1:detectors_ring
            
            if end_ind_left_row_cylinder(j)~=0
                i=end_ind_left_row_cylinder(j)+1;
            else
                i=1;
            end
            
            p=detectors_ring;
            
            out_of=1;
            out_of2=1;
            
            while true_radius(j,i)==0 && out_of
                
                i=i+1;
                
                if  i>detectors_ring
                    start_ind_right_row_cylinder(j)=0;
                    i=1;
                    out_of=0;
                end
                
            end
            
            if out_of==1
                
                start_ind_right_row_cylinder(j)=i;
                
            end
            
            if start_ind_right_row_cylinder(j)~=0
                
                while true_radius(j,p)==0 && out_of2
                    
                    p=p-1;
                    
                    if p==start_ind_right_row_cylinder(j)
                        end_ind_right_row_cylinder(j)=0;
                        p=1;
                        out_of2=0;
                    end
                end
                
                if out_of2==1
                    
                    end_ind_right_row_cylinder(j)=p;
                    
                end
                
                
            end
            
        end
        
    end
    
    
end


%% Axial block profiles and geom. factors for each ring

if options.normalization_options(1)==1
    tic
    
    if ~options.use_raw_data
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% Sinogram %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if r ~= inf && length(normalization_attenuation_correction)==2
            Sincounts = sum(sum(SinDouble ./ rad_coeff_matrix,2),1);
        else
            if exist('OCTAVE_VERSION','builtin') == 5
                if any(isnan(SinDouble))
                    temp = SinDouble;
                    temp(isnan(temp)) = 0;
                    Sincounts = sum(sum(temp,2),1);
                    clear temp
                else
                    Sincounts = sum(sum(SinDouble,2),1);
                end
            else
                Sincounts = sum(sum(SinDouble,2,'omitnan'),1,'omitnan');
            end
        end
        axial_block_profile = sqrt(mean(Sincounts(1:options.segment_table(1)))./Sincounts(1:options.segment_table(1)));
        index1 = round((z(:,1) + (z(2,1) - z(1,1))) / (z(2,1) - z(1,1)));
        index2 = round((z(:,2) + (z(2,2) - z(1,2))) / (z(2,2) - z(1,2)));
        if min(index1(:) > 1)
            index1 = index1 - min(index1(:)) + 1;
        end
        if min(index2(:) > 1)
            index2 = index2 - min(index2(:)) + 1;
        end
        SinDouble = bsxfun(@times, SinDouble, axial_block_profile(index1) .* axial_block_profile(index2));
        normalization = bsxfun(@times, normalization, axial_block_profile(index1) .* axial_block_profile(index2));
        
        if r ~= inf && length(normalization_attenuation_correction)==2
            Sincounts = sum(sum(SinDouble ./ rad_coeff_matrix,2),1);
        else
            if exist('OCTAVE_VERSION','builtin') == 5
                if any(isnan(SinDouble))
                    temp = SinDouble;
                    temp(isnan(temp)) = 0;
                    Sincounts = sum(sum(temp,2),1);
                    clear temp
                else
                    Sincounts = sum(sum(SinDouble,2),1);
                end
            else
                Sincounts = sum(sum(SinDouble,2,'omitnan'),1,'omitnan');
            end
        end
        
        axial_geom_coeffs = mean(Sincounts)./Sincounts;
        
        SinDouble = bsxfun(@times, SinDouble, axial_geom_coeffs);
        normalization = bsxfun(@times, normalization, axial_geom_coeffs);
        
        time=toc;
        if options.verbose
            disp(['Axial correction done (', num2str(time),'s)'])
        end
        
        if nargout >= 3
            varargout{3} = axial_geom_coeffs;
        end
        if nargout >= 4
            varargout{4} = axial_block_profile;
        end
        
    else
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% List Mode %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tic
        
        plane_counts = zeros(rings,1, 'single');
        
        for u=1:rings
            if r ~= inf && length(normalization_attenuation_correction)==2
                plane_counts(u) = sum(sum(true_coincidences(detectors_ring*(u-1)+1 : detectors_ring*(u), detectors_ring*(u-1)+1 : detectors_ring*u)...
                    ./ activity_coeffs(detectors_ring*(u-1)+1 : detectors_ring*(u), detectors_ring*(u-1)+1 : detectors_ring*u)));
            else
                plane_counts(u) = sum(sum(true_coincidences(detectors_ring*(u-1)+1 : detectors_ring*(u), detectors_ring*(u-1)+1 : detectors_ring*u)));
            end
        end
        
        axial_block_profile = sqrt(mean(plane_counts)./plane_counts);
        
        if exist('OCTAVE_VERSION','builtin') == 0 && exist('repelem', 'builtin') == 0
            normalization = bsxfun(@times, normalization, repeat_elem(axial_block_profile, detectors_ring));
            true_coincidences = bsxfun(@times, true_coincidences, repeat_elem(axial_block_profile, detectors_ring));
            normalization = bsxfun(@times, normalization, repeat_elem(axial_block_profile, detectors_ring)');
            true_coincidences = bsxfun(@times, true_coincidences, repeat_elem(axial_block_profile, detectors_ring)');
        else
            normalization = bsxfun(@times, normalization, repelem(axial_block_profile, detectors_ring, 1));
            true_coincidences = bsxfun(@times, true_coincidences, repelem(axial_block_profile, detectors_ring, 1));
            normalization = bsxfun(@times, normalization, repelem(axial_block_profile, detectors_ring, 1)');
            true_coincidences = bsxfun(@times, true_coincidences, repelem(axial_block_profile, detectors_ring, 1)');
        end
        if nargout >= 4
            varargout{4} = axial_block_profile;
        end
        
        axial_geom_coeffs = zeros(rings, rings, 'single');
        for u = 1 : rings
            for v = u : rings
                if u == v
                    doubling = 1;
                else
                    doubling = 2;
                end
                if r ~= inf && length(normalization_attenuation_correction)==2
                    axial_geom_coeffs(v,u) = sum(true_coincidences(1 + detectors_ring*(v-1) : detectors_ring*v, 1 + detectors_ring*(u-1) : detectors_ring*u)...
                        ./ activity_coeffs(1 + detectors_ring*(v-1) : detectors_ring*v, 1 + detectors_ring*(u-1) : detectors_ring*u),'all') / doubling;
                else
                    axial_geom_coeffs(v,u) = sum(true_coincidences(1 + detectors_ring*(v-1) : detectors_ring*v, 1 + detectors_ring*(u-1) : detectors_ring*u),'all') / doubling;
                end
            end
        end
        axial_geom_coeffs(axial_geom_coeffs == 0) = NaN;
        axial_mean = nanmean(axial_geom_coeffs(:));
        axial_geom_coeffs = axial_mean ./ axial_geom_coeffs;
        
        for u = 1 : rings
            for v = u : rings
                true_coincidences(1 + detectors_ring*(v-1) : detectors_ring*v, 1 + detectors_ring*(u-1) : detectors_ring*u) = ...
                    true_coincidences(1 + detectors_ring*(v-1) : detectors_ring*v, 1 + detectors_ring*(u-1) : detectors_ring*u) * axial_geom_coeffs(v,u);
                normalization(1 + detectors_ring*(v-1) : detectors_ring*v, 1 + detectors_ring*(u-1) : detectors_ring*u) = ...
                    normalization(1 + detectors_ring*(v-1) : detectors_ring*v, 1 + detectors_ring*(u-1) : detectors_ring*u) * axial_geom_coeffs(v,u);
            end
        end
        
        
        if nargout >= 3
            varargout{3} = axial_geom_coeffs;
        end
        
        time=toc;
        if options.verbose
            disp(['Axial correction done (', num2str(time),'s)'])
        end
        
    end
    
end




%% Crystal interference correction
if options.normalization_options(4)==1 || options.normalization_options(2)==1 || options.normalization_options(2)==2 || options.normalization_options(3)==1
    
    tic
    
    %Detector block profiles
    
    %Find first detectors position in block
    
    for i=1:koko
        
        if detectors_x(i+1)-detectors_x(i)~=0 && detectors_x(i+2)-detectors_x(i+1)~=0
            
            k1=(detectors_y(i+1)-detectors_y(i))/(detectors_x(i+1)-detectors_x(i));
            k2=(detectors_y(i+2)-detectors_y(i+1))/(detectors_x(i+2)-detectors_x(i+1));
            
        else
            
            k1=(detectors_x(i+1)-detectors_x(i));
            k2=(detectors_x(i+2)-detectors_x(i+1));
            
        end
        
        if norm(k1-k2)>abs(k1)/10
            break
        end
        
    end
    
    %First detectors position is i+1 from blocks edge
    
    starting_block=cryst_per_block-i;
    current_block=starting_block;
    det_position=zeros(1,detectors_ring, 'single');
    
    %Determine detectors positions
    
    for i=1:length(detectors_x)
        
        if floor(cryst_per_block/2)<current_block
            
            det_position(i)=floor(norm(floor(cryst_per_block/2)-current_block));
            
        else
            
            det_position(i)=floor(norm(floor(cryst_per_block/2)-current_block))+1;
            
        end
        
        current_block=current_block+1;
        
        if current_block==cryst_per_block+1
            
            current_block=1;
            
        end
        
    end
    
    %Averaged data for each block/radial profile
    
    %Different block position combinations = (cryst_per_block/2)^2/2
    
    if ~options.use_raw_data
        
        %         avg_profile=zeros(cryst_per_block/2,cryst_per_block/2,options.Ndist);
        avg_profile_LOR_amount = zeros(ceil(cryst_per_block/2),ceil(cryst_per_block/2),options.Ndist, 'single');
        %         not_found=0;
        detector_index_start=zeros(1,Nang*options.Ndist, 'single');
        detector_index_end=detector_index_start;
        index=zeros(options.Ndist*Nang,2, 'single');
        
        x = permute(reshape(x, options.Ndist, Nang, 2), [2 1 3]);
        x = reshape(x, options.Ndist*Nang, 2);
        y = permute(reshape(y, options.Ndist, Nang, 2), [2 1 3]);
        y = reshape(y, options.Ndist*Nang, 2);
        
        for j=1:options.Ndist*Nang
            
            p=1;
            t=0;
            alku = true;
            loppu = true;
            
            
            if mashing > 1
                locx = abs(detectors_x-x(j,1)) - min(abs(detectors_x-x(j,1))) < 1e-5;
                locy = abs(detectors_y-y(j,1)) - min(abs(detectors_y-y(j,1))) < 1e-5;
                locx = find(locx);
                locy = find(locy);
                detector_index_start(j) = det_position(min(locx(ismember(locx,locy))));
                locx = abs(detectors_x-x(j,2)) - min(abs(detectors_x-x(j,2))) < 1e-5;
                locy = abs(detectors_y-y(j,2)) - min(abs(detectors_y-y(j,2))) < 1e-5;
                locx = find(locx);
                locy = find(locy);
                detector_index_end(j) = det_position(max(locx(ismember(locx,locy))));
            else
                
                while t~=2 &&  p<=length(det_position)
                    
                    if x(j,1)==detectors_x(p) && y(j,1)==detectors_y(p) && alku
                        
                        detector_index_start(j)=det_position(p);
                        t=t+1;
                        alku = false;
                        
                    end
                    
                    
                    if x(j,2)==detectors_x(p) && y(j,2)==detectors_y(p) && loppu
                        
                        detector_index_end(j)=det_position(p);
                        t=t+1;
                        loppu = false;
                        
                    end
                    p=p+1;
                end
                
                if t~=2
                    %                 not_found=not_found+1;
                    error('No corresponding coordinates found')
                end
            end
            
            index(j,:)=[max(detector_index_start(j),detector_index_end(j))...
                min(detector_index_start(j),detector_index_end(j))];
            
            
            avg_profile_LOR_amount(index(j,1),index(j,2),ceil(j/Nang))=...
                avg_profile_LOR_amount(index(j,1),index(j,2),ceil(j/Nang)) + 1;
            
        end
        
        block_profile_coeffs=zeros(size(avg_profile_LOR_amount), 'single');
        
        
        
        %Calculate block profile coeffs by averaging over all sinograms%
        
        
        if options.normalization_options(3)==1 && options.normalization_options(2)==2
            
            block_profile_matrix=zeros(Nang,options.Ndist, 'single');
            
            avg_profile=zeros(size(avg_profile_LOR_amount), 'single');
            
            if r ~= inf && length(normalization_attenuation_correction)==2
                SinDouble = permute(SinDouble ./ rad_coeff_matrix, [2 1 3]);
            else
                SinDouble = permute(SinDouble, [2 1 3]);
            end
            
            for j=1:options.Ndist*Nang
                
                avg_profile(index(j,1),index(j,2),ceil(j/Nang))=avg_profile...
                    (index(j,1),index(j,2),ceil(j/Nang))+sum(SinDouble(j-floor((j-1)/Nang)...
                    *Nang,ceil(j/Nang),:));
                
            end
            
            
            
            for j=1:options.Ndist*Nang
                if avg_profile(index(j,1),index(j,2),ceil(j/Nang))==0
                    
                    avg_profile(index(j,1),index(j,2),ceil(j/Nang))=...
                        min(min(avg_profile(:,:,ceil(j/Nang))));
                    
                end
            end
            
            avg_profile=avg_profile./(avg_profile_LOR_amount.*sino_amount);
            
            for j=1:options.Ndist
                temp_mat=avg_profile(:,:,j);
                
                if nansum(nansum(avg_profile(:,:,j)))~=0
                    block_profile_coeffs(:,:,j)=nanmean(temp_mat(:))./avg_profile(:,:,j);
                    
                end
            end
            %             index1=reshape(index(:,1),options.Ndist,Nang);
            %             index2=reshape(index(:,2),options.Ndist,Nang);
            index1=reshape(index(:,1),Nang,options.Ndist);
            index2=reshape(index(:,2),Nang,options.Ndist);
            
            %coeff matrix for sinograms
            for p=1:Nang
                
                for i=1:options.Ndist
                    
                    if block_profile_coeffs(index1(p,i),index2(p,i),i)==inf
                        temp_coeffs=block_profile_coeffs(:,:,i);
                        block_profile_coeffs(index1(p,i),index2(p,i),i)=max(max(temp_coeffs(temp_coeffs~=inf)));
                    end
                    if isnan(block_profile_coeffs(index1(p,i),index2(p,i),i))
                        temp_coeffs=block_profile_coeffs(:,:,i);
                        block_profile_coeffs(index1(p,i),index2(p,i),i)=max(max(temp_coeffs(~isnan(temp_coeffs))));
                    end
                    block_profile_matrix(p,i)=block_profile_coeffs(index1(p,i),index2(p,i),i);
                    
                end
            end
            
            x = permute(reshape(x, Nang, options.Ndist, 2), [2 1 3]);
            x = reshape(x, options.Ndist*Nang, 2);
            y = permute(reshape(y, Nang, options.Ndist, 2), [2 1 3]);
            y = reshape(y, options.Ndist*Nang, 2);
            
            block_profile_matrix=repmat(block_profile_matrix,1,1,sino_amount);
            
            SinDouble=SinDouble.*block_profile_matrix;
            
            SinDouble = permute(SinDouble, [2 1 3]);
            block_profile_matrix = permute(block_profile_matrix, [2 1 3]);
            normalization=normalization.*block_profile_matrix;
            
            if nargout >= 5
                varargout{5} = block_profile_matrix;
            end
            
            time=toc;
            if options.verbose
                disp(['Transaxial block profile correction done (', num2str(time),'s)'])
            end
            
        end
        
    else
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% List-mode form correction %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if options.normalization_options(2)==2 || options.normalization_options(4)==1 || options.normalization_options(3)==1
            
            radius_raw_data=zeros(detectors_ring, 'single');
            for j=1:detectors_ring
                
                if (start_ind_low_col(j)~=0 && end_ind_low_col(j)~=0) || (start_ind_up_col(j)~=0 && end_ind_up_col(j)~=0)
                    
                    for i=[start_ind_up_col(j):end_ind_up_col(j) start_ind_low_col(j):end_ind_low_col(j)]
                        
                        if i~=0
                            
                            radius_raw_data(i,j)=norm([detectors_x(j)-detectors_x(i) detectors_y(j)-detectors_y(i)]);
                            
                        end
                        
                    end
                    
                end
                
            end
            
            sorting_factor=cryst_per_block;
            
            max_rad=max(max(radius_raw_data));
            min_rad=min(min(nonzeros(radius_raw_data)));
            environment=(max_rad-min_rad)/(ceil(detectors_ring/sorting_factor)-1);
            
            %radial sorting
            radial_index=zeros(1,ceil(detectors_ring/sorting_factor)-1, 'single');
            
            for i=1:length(radial_index)
                
                radial_index(i)=min_rad+environment/2+(i-1)*environment;
                
            end
            
            radial_index=nonzeros(radial_index)';
            det_index=zeros(size(radius_raw_data), 'single');
            
            
            %Sort LORs by radial profile
            
            for j=1:detectors_ring
                
                if (start_ind_low_col(j)~=0 && end_ind_low_col(j)~=0) || (start_ind_up_col(j)~=0 && end_ind_up_col(j)~=0)
                    
                    for i=[start_ind_up_col(j):end_ind_up_col(j) start_ind_low_col(j):end_ind_low_col(j)]
                        
                        if i~=0
                            
                            p=1;
                            over_index=0;
                            while norm(radius_raw_data(i,j)-radial_index(p))>environment/2+10^-8 && over_index==0
                                
                                p=p+1;
                                
                                if  p>length(radial_index)
                                    over_index=1;
                                    p=1;
                                end
                                
                            end
                            
                            if over_index==0
                                
                                det_index(i,j)=p;
                                
                            else
                                
                                det_index(i,j)=0;
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
            
            if options.normalization_options(2)==2 || options.normalization_options(2)==0
                
                profile_hits=zeros(max(det_position),max(det_position),max(max((det_index))), 'single');
                profile_avg=profile_hits;
                
                %                 start=1;
                
                for j=1:detectors_ring
                    
                    if (start_ind_low_col(j)~=0 && end_ind_low_col(j)~=0) || (start_ind_up_col(j)~=0 && end_ind_up_col(j)~=0)
                        
                        
                        for i=[start_ind_up_col(j):end_ind_up_col(j) start_ind_low_col(j):end_ind_low_col(j)]
                            
                            if i~=0
                                
                                if det_index(i,j)~=0
                                    
                                    y_ind=max([det_position(i),det_position(j)]);
                                    x_ind=min([det_position(i),det_position(j)]);
                                    z_ind=det_index(i,j);
                                    
                                    if start_ind_low_col(j)<=i && i<=end_ind_low_col(j)
                                        
                                        start=2;
                                        
                                    else
                                        
                                        start=1;
                                        
                                    end
                                    
                                    for p=start:rings
                                        
                                        for t=p:rings
                                            
                                            sec_ind=j+(p-1)*detectors_ring;
                                            
                                            profile_avg(y_ind,x_ind,z_ind)=...
                                                profile_avg(y_ind,x_ind,z_ind)+true_coincidences(i+(t-1)*detectors_ring,sec_ind);
                                            
                                        end
                                    end
                                    
                                    if start==2
                                        profile_hits(y_ind,x_ind,z_ind)...
                                            =profile_hits(y_ind,x_ind,z_ind)+(rings^2/2-rings/2);
                                    else
                                        profile_hits(y_ind,x_ind,z_ind)...
                                            =profile_hits(y_ind,x_ind,z_ind)+(rings^2/2+rings/2);
                                        
                                    end
                                    
                                end
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
                profile_avg=profile_avg./profile_hits;
                profile_avg(isnan(profile_avg))=0;
                %Mean counts for different radial groups
                mean_profile=zeros(1,max(max(det_index)), 'single');
                block_profile_coeffs=zeros(size(profile_avg),'single');
                for u=1:max(max(det_index))
                    
                    mean_profile(u)=sum(sum(profile_avg(:,:,u)))./numel(nonzeros(profile_avg(:,:,u)));
                    
                    %coeffs
                    
                    block_profile_coeffs(:,:,u)=mean_profile(u)./profile_avg(:,:,u);
                    
                end
                block_profile_coeffs(block_profile_coeffs==inf)=0;
                block_profile_coeffs(isnan(block_profile_coeffs))=0;
                
                %correct list mode data
                
                block_matrix=zeros(size(true_coincidences),'single');
                
                for j=1:detectors_ring
                    
                    for i=[start_ind_low_col(j):end_ind_low_col(j) start_ind_up_col(j):end_ind_up_col(j)]
                        
                        if i~=0
                            
                            y_ind=max([det_position(i),det_position(j)]);
                            x_ind=min([det_position(i),det_position(j)]);
                            
                            
                            if det_index(i,j)~=0
                                
                                current_coeff=block_profile_coeffs(y_ind,x_ind,det_index(i,j));
                                if current_coeff==0
                                    current_coeff=1;
                                end
                                
                                for p=1:rings
                                    
                                    sec_ind=j+(p-1)*detectors_ring;
                                    
                                    for t=p:rings
                                        
                                        true_coincidences(i+(t-1)*detectors_ring,sec_ind)=true_coincidences(i+(t-1)*detectors_ring,sec_ind)...
                                            *current_coeff;
                                        
                                        normalization(i+(t-1)*detectors_ring,sec_ind)=normalization(i+(t-1)*detectors_ring,sec_ind)...
                                            *current_coeff;
                                        
                                        block_matrix(i+(t-1)*detectors_ring,sec_ind)=current_coeff;
                                        
                                    end
                                    
                                end
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
                if nargout >= 5
                    varargout{5} = block_matrix;
                end
            end
            
            
            time=toc;
            if options.verbose
                disp(['Block profile correction done (', num2str(time),'s)'])
            end
        end
        
        
    end
    
end

%% Detector effiency factors

if options.normalization_options(2)~=0
    
    tic
    
    if ~options.use_raw_data
        
        
        %  with 3-D fan-sum algorithm
        
        
        %effiency for LOR(i,j) = detector_effiency(i) x detector_effiency(j)
        
        z_rings=z(1:options.segment_table(1),1);
        det_num=zeros(options.Ndist,Nang,2,'single');
        counts_det=zeros(detectors_ring,options.segment_table(1),'single');
        %         coeff_matrix=zeros(size(SinDouble));
        ring = round((z - z(1,1)) ./ (z(2,1) - z(1,1))) + 1;
        
        %Determine each LORS detector numbers
        
        for u=1:options.Ndist*Nang
            
            found_it=0;
            t=1;
            while found_it==0
                
                if mashing > 1
                    if norm([x(u,1)-detectors_x(t) y(u,1)-detectors_y(t)])<0.3
                        %                 if norm([x(u,1)-detectors_x(t) y(u,1)-detectors_y(t)]) == min(vecnorm([x(u,1)-detectors_x y(u,1)-detectors_y]'))
                        
                        if norm([x(u,1)-detectors_x(t+1) y(u,1)-detectors_y(t+1)])<0.3
                            det_num(floor((u-1)/Nang)+1,u-floor((u-1)/Nang)*Nang,1)=t + 1;
                        else
                            det_num(floor((u-1)/Nang)+1,u-floor((u-1)/Nang)*Nang,1)=t;
                        end
                        found_it=1;
                        
                    end
                else
                    if norm([x(u,1)-detectors_x(t) y(u,1)-detectors_y(t)])<0.001
                        
                        det_num(floor((u-1)/Nang)+1,u-floor((u-1)/Nang)*Nang,1)=t;
                        found_it=1;
                        
                    end
                end
                t=t+1;
            end
            found_it=0;
            t=1;
            while found_it==0
                
                if mashing > 1
                    if norm([x(u,2)-detectors_x(t) y(u,2)-detectors_y(t)])<0.3
                        % if norm([x(u,1)-detectors_x(t) y(u,1)-detectors_y(t)]) == min(vecnorm([x(u,1)-detectors_x y(u,1)-detectors_y]'))
                        
                        if norm([x(u,2)-detectors_x(t+1) y(u,2)-detectors_y(t+1)])<0.3
                            det_num(floor((u-1)/Nang)+1,u-floor((u-1)/Nang)*Nang,2)=t + 1;
                        else
                            det_num(floor((u-1)/Nang)+1,u-floor((u-1)/Nang)*Nang,2)=t;
                        end
                        found_it=1;
                        
                    end
                else
                    
                    if norm([x(u,2)-detectors_x(t) y(u,2)-detectors_y(t)])<0.001
                        
                        det_num(floor((u-1)/Nang)+1,u-floor((u-1)/Nang)*Nang,2)=t;
                        found_it=1;
                        
                    end
                end
                t=t+1;
            end
            
        end
        
        hits_det=zeros(detectors_ring,options.segment_table(1),'single');
        testi1 = det_num(:,:,1);
        testi2 = det_num(:,:,2);
        testi1 = testi1(:);
        testi2 = testi2(:);
        testi1 = bsxfun(@plus, testi1, (ring(:,1)' - 1) * detectors_ring);
        testi2 = bsxfun(@plus, testi2, (ring(:,2)' - 1) * detectors_ring);
        testi1 = testi1(:);
        testi2 = testi2(:);
        maksimi1 = max(testi1(:));
        maksimi2 = max(testi2(:));
        
        if r ~= inf && length(normalization_attenuation_correction)==2
            counts_det(1:maksimi1) = accumarray(testi1, SinDouble(:) ./ rad_coeff_matrix(:));
            counts_det(1:maksimi2) = counts_det(1:maksimi2)' + accumarray(testi2, SinDouble(:) ./ rad_coeff_matrix(:));
        else
            counts_det(1:maksimi1) = accumarray(testi1, SinDouble(:));
            counts_det(1:maksimi2) = counts_det(1:maksimi2)' + accumarray(testi2, SinDouble(:));
        end
        hits_det(1:maksimi1) = accumarray(testi1(:),1);
        hits_det(1:maksimi2) = hits_det(1:maksimi2)' + accumarray(testi2(:),1);
        
        %Calculate coeffs (mean inverse)
        
        counts_det = counts_det./hits_det;
        if detectors_ring > options.det_per_ring
            keskiarvo = mean(counts_det(any(counts_det,2),:),1);
            coeffs_detectors = bsxfun(@rdivide, keskiarvo, counts_det);
            coeffs_detectors(isinf(coeffs_detectors)) = 0;
        else
            coeffs_detectors = bsxfun(@rdivide, mean(counts_det,1), counts_det);
        end
        
        
        if options.normalization_options(3) == 1 && options.normalization_options(2)~=2
            
            detector_index_start=reshape(detector_index_start,options.Ndist,Nang);
            detector_index_end=reshape(detector_index_end,options.Ndist,Nang);
            det_number_block=zeros(1,detectors_ring,'single');
            
            block_profile_counts = accumarray(repmat(detector_index_start(:), sino_amount, 1), coeffs_detectors(testi1));
            block_profile_counts = block_profile_counts' + accumarray(repmat(detector_index_end(:), sino_amount, 1), coeffs_detectors(testi2))';
            block_profile_hits = accumarray(detector_index_start(:),1);
            block_profile_hits = block_profile_hits + accumarray(detector_index_end(:),1);
            block_profile_hits = block_profile_hits' * sino_amount;
        end
        coeff_matrix = reshape(coeffs_detectors(det_num(:,:,1),ring(:,1)) .* coeffs_detectors(det_num(:,:,2),ring(:,2)), options.Ndist, Nang, sino_amount);
        if options.normalization_options(3) == 1 && options.normalization_options(2)~=2
            
            p=1;
            
            for k=1:detectors_ring
                
                for i=1:options.Ndist
                    
                    for j=1:Nang
                        
                        if det_num(i,j,1)==p
                            
                            det_number_block(p)=detector_index_start(i,j);
                            p=p+1;
                            
                        elseif det_num(i,j,2)==p
                            
                            det_number_block(p)=detector_index_end(i,j);
                            p=p+1;
                            
                        end
                        
                        if p>detectors_ring
                            
                            break
                            
                        end
                        
                    end
                    
                    if p>detectors_ring
                        
                        break
                        
                    end
                    
                end
                
                if p>detectors_ring
                    
                    break
                    
                end
                
            end
            
        end
        %Correct data
        SinDouble=SinDouble.*coeff_matrix(1:options.Ndist,:,:);
        
        %Coeffs to normalization matrix
        normalization=normalization.*coeff_matrix(1:options.Ndist,:,:);
        
        if nargout >= 6
            varargout{6} = coeff_matrix;
        end
        
        if options.normalization_options(3) == 1 && options.normalization_options(2)~=2 && nargout >= 5
            block_profile_counts=block_profile_counts./block_profile_hits;
            block_profile_coeffs=mean(block_profile_counts)./block_profile_counts;
            %Separate block profile coeffs from effiency coeffs
            block_matrix=zeros(detectors_ring,rings,'single');
            compared=(z_rings(2)-z_rings(1));
            
            for t=1:rings
                k=1;
                while abs(z_true(t)-z_rings(k))>compared
                    k=k+1;
                end
                
                block_matrix(:,t) = coeffs_detectors(:,k)./block_profile_coeffs(det_number_block)';
                
            end
            
            if nargout >= 5
                varargout{5} = block_matrix;
            end
        end
        
    else
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% List mode correction %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if isempty(normalization_attenuation_correction)
            
            true_radius=ones(detectors_ring,'single');
            
        else
            
            if r~=inf
                
                if r==inf
                    true_radius=ones(detectors_ring,'single');
                end
                
            end
            
        end
        
        if options.normalization_options(2)==2
            
            %Detector effiency correction with SP-C method (https://doi.org/10.1088/0031-9155/43/1/012, page: 193)
            
            elements_x_ax=zeros(1,detectors_ring,'single');
            elements_y_ax=elements_x_ax;
            
            
            %elements numbers for detectors and projections between groups
            
            if r~=inf
                
                end_ind_low_col=end_ind_low_col_cylinder;
                start_ind_low_col=start_ind_low_col_cylinder;
                end_ind_up_col=end_ind_up_col_cylinder;
                start_ind_up_col=start_ind_up_col_cylinder;
                end_ind_left_row=end_ind_left_row_cylinder;
                start_ind_left_row=start_ind_left_row_cylinder;
                end_ind_right_row=end_ind_right_row_cylinder;
                start_ind_right_row=start_ind_right_row_cylinder;
                
            end
            %Element amounts of each detectors projection group
            
            for j=1:detectors_ring
                
                if end_ind_low_col(j)~=0
                    elements_y_ax(j)=end_ind_low_col(j)-start_ind_low_col(j)+1;
                end
                if end_ind_up_col(j)~=0
                    elements_y_ax(j)=elements_y_ax(j)+end_ind_up_col(j)-start_ind_up_col(j)+1;
                end
                
                if end_ind_left_row(j)~=0
                    elements_x_ax(j)=end_ind_left_row(j)-start_ind_left_row(j)+1;
                end
                
                if end_ind_right_row(j)~=0
                    elements_x_ax(j)=elements_x_ax(j)+end_ind_right_row(j)-start_ind_right_row(j);
                end
                
            end
            
            current_plane=true_det_indices(detectors_ring+1:detectors_ring*2,1:detectors_ring);
            
            elements=zeros(detectors_ring,'single');
            elements_plane=elements;
            
            %elements of projections between groups A and B
            
            for j=1:detectors_ring
                
                for i=[start_ind_low_col(j):end_ind_low_col(j) start_ind_up_col(j):end_ind_up_col(j)]
                    
                    if i~=0
                        
                        if end_ind_low_col(j)~=0 && start_ind_left_row(i)~=0
                            
                            elements_plane(i,j)=elements_plane(i,j)+numel(nonzeros(current_plane...
                                (start_ind_low_col(j):end_ind_low_col(j),start_ind_left_row(i):end_ind_left_row(i))));
                            
                            elements(i,j)=elements(i,j)+numel(nonzeros(current_plane...
                                (start_ind_low_col(j):end_ind_low_col(j),start_ind_left_row(i):end_ind_left_row(i))));
                            
                        end
                        
                        if start_ind_right_row(i)~=0 && end_ind_low_col(j)~=0
                            
                            elements(i,j)=elements(i,j)+numel(nonzeros(current_plane...
                                (start_ind_low_col(j):end_ind_low_col(j),start_ind_right_row(i):end_ind_right_row(i))));
                            
                        end
                        
                        if end_ind_up_col(j)~=0 && start_ind_left_row(i)~=0
                            
                            elements(i,j)=elements(i,j)+numel(nonzeros(current_plane...
                                (start_ind_up_col(j):end_ind_up_col(j),start_ind_left_row(i):end_ind_left_row(i))));
                            
                        end
                        
                        if start_ind_right_row(i)~=0 && end_ind_up_col(j)~=0
                            
                            elements(i,j)=elements(i,j)+numel(nonzeros(current_plane...
                                (start_ind_up_col(j):end_ind_up_col(j),start_ind_right_row(i):end_ind_right_row(i))));
                            
                        end
                        
                    end
                    
                end
                
            end
            
            %Get FOV projections between Groups A and B
            
            mean_detector_pairs=zeros(size(current_plane),'single');
            
            effiencies_y_mean=zeros(1,detectors_ring,'single');
            effiencies_x_mean=effiencies_y_mean;
            
            
            total_elements=0;
            
            %Total true detector pairs
            for j=1:detectors_ring
                
                total_elements=total_elements+(end_ind_up_col(j)-start_ind_up_col(j))+(end_ind_low_col(j)-start_ind_low_col(j));
                
                if end_ind_up_col(j)~=0
                    
                    total_elements=total_elements+1;
                    
                end
                
                if end_ind_low_col(j)~=0
                    
                    total_elements=total_elements+1;
                    
                end
                
            end
            
            %Mean counts (detector sensitivity)
            
            cylinder_indices=true_radius~=0;
            
            det_coeffs = zeros(size(normalization),'single');
            
            for u=1:options.ring_difference+1
                
                for k=1:rings-u+1
                    
                    Sum_AB=zeros(detectors_ring,'single');
                    
                    if r~=inf && uniform_source==0
                        
                        current_plane=true_coincidences(detectors_ring*(u+k-2)+1:detectors_ring*(u+k-1),detectors_ring*(k-1)+1:detectors_ring*k)...
                            .*activity_coeffs(detectors_ring*(u+k-2)+1:detectors_ring*(u+k-1),detectors_ring*(k-1)+1:detectors_ring*k).*cylinder_indices;
                        
                    else
                        
                        current_plane=true_coincidences(detectors_ring*(u+k-2)+1:detectors_ring*(u+k-1),detectors_ring*(k-1)+1:detectors_ring*k);
                        
                    end
                    
                    if u==1
                        
                        for j=1:detectors_ring
                            
                            for i=j:detectors_ring
                                
                                current_plane(j,i)=current_plane(i,j);
                                
                            end
                        end
                    end
                    
                    mean_counts=sum(sum(current_plane(cylinder_indices)))/total_elements;
                    
                    effiencies=mean_counts./current_plane;
                    effiencies=effiencies.*cylinder_indices;
                    effiencies(effiencies==inf)=1;
                    effiencies(isnan(effiencies))=0;
                    
                    for j=1:detectors_ring
                        
                        effiencies_y_mean(j)=sum(effiencies(j,:))/elements_x_ax(j);
                        effiencies_x_mean(j)=sum(effiencies(:,j))/elements_y_ax(j);
                        
                    end
                    
                    for j=1:detectors_ring
                        
                        for i=[start_ind_low_col(j):end_ind_low_col(j) start_ind_up_col(j):end_ind_up_col(j)]
                            
                            if i~=0
                                mean_detector_pairs(i,j)=effiencies_y_mean(i)*effiencies_x_mean(j);
                            end
                            
                        end
                        
                    end
                    
                    
                    first_sum_low=1;
                    first_sum_up=1;
                    
                    first_sum_box_low=0;
                    first_sum_box_up=0;
                    
                    first_sum_low2=1;
                    first_sum_up2=1;
                    
                    first_sum_box_low2=0;
                    first_sum_box_up2=0;
                    
                    
                    for j=1:max(end_ind_left_row)
                        
                        Sum_box_low=first_sum_box_low;
                        Sum_box_up=first_sum_box_up;
                        
                        if start_ind_low_col(j)~=0
                            
                            new_i=start_ind_low_col(j);
                            
                            %%% MOVE SUMBOX TO NEW COLUMN COORDINATES %%%
                            
                            if j>1
                                
                                last_i=start_ind_low_col(j-1);
                                
                                column_diff=new_i-last_i;
                                
                                if start_ind_right_row(new_i)~=0 && start_ind_left_row(new_i)
                                    
                                    if column_diff~=0
                                        
                                        Sum_box_low=Sum_box_low-sum(sum(current_plane(last_i:...
                                            new_i-1,start_ind_right_row(last_i):end_ind_right_row(new_i))))...
                                            -sum(sum(current_plane(last_i:...
                                            new_i-1,start_ind_left_row(last_i):end_ind_left_row(new_i))));
                                        
                                    end
                                    
                                else
                                    
                                    if start_ind_left_row(last_i)~=0
                                        
                                        if column_diff~=0
                                            
                                            Sum_box_low=Sum_box_low-sum(sum(current_plane(last_i:...
                                                new_i-1,start_ind_left_row(last_i):end_ind_left_row(new_i))));
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            for i=new_i:end_ind_low_col(j)
                                
                                if i~=0
                                    
                                    if first_sum_low==1
                                        
                                        if end_ind_right_row(i)~=0
                                            
                                            Sum_box_low=sum(sum(current_plane(new_i:end_ind_low_col(j),start_ind_left_row(i):end_ind_left_row(i))));
                                            Sum_box_low=Sum_box_low+sum(sum(current_plane(new_i:end_ind_low_col(j),start_ind_right_row(i):end_ind_right_row(i))));
                                            
                                        else
                                            
                                            Sum_box_low=sum(sum(current_plane(new_i:end_ind_low_col(j),start_ind_left_row(i):end_ind_left_row(i))));
                                            
                                        end
                                        
                                        
                                        
                                        first_sum_low=0;
                                        first_sum_box_low=Sum_box_low;
                                        
                                        Sum_AB(i,j)=Sum_box_low;
                                        
                                    else
                                        
                                        
                                        %%% MOVE SUMBOX TO NEW ROW COORDINATES %%%
                                        
                                        if i==new_i
                                            
                                            row_diff=[start_ind_left_row(i)-start_ind_left_row(last_i) ...
                                                end_ind_left_row(i)-end_ind_left_row(last_i) start_ind_right_row(i)-...
                                                start_ind_right_row(last_i) end_ind_right_row(i)-end_ind_right_row(last_i)];
                                            
                                            last_row_start_left=start_ind_left_row(last_i);
                                            last_row_end_left=end_ind_left_row(last_i);
                                            
                                            last_row_start_right=start_ind_right_row(last_i);
                                            last_row_end_right=end_ind_right_row(last_i);
                                            
                                            temp_ind=last_i;
                                            if i~=1
                                                
                                                if  end_ind_right_row(i)==0 && end_ind_right_row(temp_ind)~=0
                                                    
                                                    row_diff(4)=0;
                                                    row_diff(3)=0;
                                                    
                                                    Sum_box_low=Sum_box_low-sum(sum(current_plane(start_ind_low_col(j-1):...
                                                        end_ind_low_col(j-1),start_ind_right_row(temp_ind):end_ind_right_row(temp_ind))));
                                                    
                                                end
                                                
                                            end
                                            
                                        else
                                            
                                            temp_col=j;
                                            
                                            row_diff=[start_ind_left_row(i)-start_ind_left_row(i-1) end_ind_left_row(i)-end_ind_left_row(i-1)...
                                                start_ind_right_row(i)-start_ind_right_row(i-1) end_ind_right_row(i)-end_ind_right_row(i-1)];
                                            
                                            last_row_start_left=start_ind_left_row(i-1);
                                            last_row_end_left=end_ind_left_row(i-1);
                                            
                                            last_row_start_right=start_ind_right_row(i-1);
                                            last_row_end_right=end_ind_right_row(i-1);
                                            
                                            temp_ind=i-1;
                                            
                                            if i~=1
                                                
                                                if  end_ind_right_row(i)==0 && end_ind_right_row(temp_ind)~=0
                                                    
                                                    row_diff(4)=0;
                                                    row_diff(3)=0;
                                                    
                                                    Sum_box_low=Sum_box_low-sum(sum(current_plane(start_ind_low_col(temp_col):...
                                                        end_ind_low_col(temp_col),start_ind_right_row(temp_ind):end_ind_right_row(temp_ind))));
                                                    
                                                end
                                                
                                            end
                                            
                                        end
                                        
                                        %Start left box
                                        
                                        if row_diff(1)~=0
                                            
                                            Sum_box_low=Sum_box_low-sum(sum(current_plane(start_ind_low_col(j):...
                                                end_ind_low_col(j),last_row_start_left:start_ind_left_row(i)-1)));
                                            
                                        end
                                        
                                        %End left box
                                        
                                        if row_diff(2)~=0
                                            
                                            Sum_box_low=Sum_box_low+sum(sum(current_plane(start_ind_low_col(j):...
                                                end_ind_low_col(j),last_row_end_left+1:end_ind_left_row(i))));
                                            
                                        end
                                        
                                        
                                        if last_row_start_right~=0 && start_ind_right_row(i)~=0
                                            
                                            % Start right box
                                            
                                            if row_diff(3)~=0
                                                
                                                Sum_box_low=Sum_box_low-sum(sum(current_plane(start_ind_low_col(j):...
                                                    end_ind_low_col(j),last_row_start_right:start_ind_right_row(i)-1)));
                                                
                                            end
                                            
                                            % End right box
                                            
                                            if row_diff(4)~=0
                                                
                                                Sum_box_low=Sum_box_low+sum(sum(current_plane(start_ind_low_col(j):end_ind_low_col...
                                                    (j),last_row_end_right+1:end_ind_right_row(i))));
                                                
                                            end
                                            
                                        end
                                        
                                        Sum_AB(i,j)=Sum_box_low;
                                        
                                        
                                        if i==start_ind_low_col(j)
                                            
                                            first_sum_box_low=Sum_box_low;
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            
                        end
                        
                        %%% UPPER SUMBOX
                        
                        if start_ind_up_col(j)~=0
                            
                            new_i_up=start_ind_up_col(j);
                            new_i_low=start_ind_low_col(j);
                            %%% MOVE SUMBOX TO NEW COLUMN COORDINATES %%%
                            
                            
                            if start_ind_up_col(j-1)~=0
                                
                                last_i_up=start_ind_up_col(j-1);
                                last_i_low=start_ind_low_col(j-1);
                                
                                column_diff=start_ind_up_col(j)-start_ind_up_col(j-1);
                                
                                if start_ind_right_row(last_i_low)~=0 && start_ind_left_row(last_i_low)~=0
                                    
                                    if column_diff~=0
                                        
                                        Sum_box_up=Sum_box_up-sum(sum(current_plane(last_i_up:...
                                            new_i_up-1,start_ind_right_row(last_i_low):end_ind_right_row(new_i_low))))...
                                            -sum(sum(current_plane(last_i_up:...
                                            new_i_up-1,start_ind_left_row(last_i_low):end_ind_left_row(new_i_low))));
                                        
                                    end
                                    
                                else
                                    
                                    if start_ind_left_row(last_i_low)~=0
                                        
                                        if column_diff~=0
                                            
                                            Sum_box_up=Sum_box_up-sum(sum(current_plane(last_i_up:...
                                                new_i_up-1,start_ind_left_row(last_i_low):end_ind_left_row(new_i_low))));
                                            
                                        end
                                        
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            for i=start_ind_low_col(j):end_ind_low_col(j)
                                
                                if i~=0
                                    
                                    if first_sum_up==1
                                        
                                        if start_ind_right_row(i)~=0
                                            
                                            Sum_box_up=sum(sum(current_plane(new_i_up:end_ind_up_col(j),...
                                                start_ind_left_row(i):end_ind_left_row(i))));
                                            Sum_box_up=Sum_box_up+sum(sum(current_plane(new_i_up:end_ind_up_col(j),...
                                                start_ind_right_row(i):end_ind_right_row(i))));
                                            
                                        else
                                            
                                            Sum_box_up=sum(sum(current_plane(new_i_up:end_ind_up_col(j),...
                                                start_ind_left_row(i):end_ind_left_row(i))));
                                            
                                        end
                                        
                                        first_sum_up=0;
                                        first_sum_box_up=Sum_box_up;
                                        
                                        Sum_AB(i,j)=Sum_AB(i,j)+Sum_box_up;
                                        
                                    else
                                        
                                        
                                        %%% MOVE SUMBOX TO NEW ROW COORDINATES %%%
                                        
                                        if i==start_ind_low_col(j)
                                            
                                            row_diff=[start_ind_left_row(i)-start_ind_left_row(start_ind_low_col(j-1)) ...
                                                end_ind_left_row(i)-end_ind_left_row(start_ind_low_col(j-1)) start_ind_right_row(i)-...
                                                start_ind_right_row(start_ind_low_col(j-1)) end_ind_right_row(i)-end_ind_right_row(start_ind_low_col(j-1))];
                                            
                                            last_row_start_left=start_ind_left_row(start_ind_low_col(j-1));
                                            last_row_end_left=end_ind_left_row(start_ind_low_col(j-1));
                                            
                                            last_row_start_right=start_ind_right_row(start_ind_low_col(j-1));
                                            last_row_end_right=end_ind_right_row(start_ind_low_col(j-1));
                                            
                                            temp_ind=start_ind_low_col(j-1);
                                            
                                            if  end_ind_right_row(i)==0 && end_ind_right_row(temp_ind)~=0
                                                
                                                row_diff(4)=0;
                                                row_diff(3)=0;
                                                
                                                Sum_box_up=Sum_box_up-sum(sum(current_plane(start_ind_low_col(j-1):...
                                                    end_ind_low_col(j-1),start_ind_right_row(temp_ind):end_ind_right_row(temp_ind))));
                                                
                                            end
                                            
                                            
                                        else
                                            
                                            row_diff=[start_ind_left_row(i)-start_ind_left_row(i-1) end_ind_left_row(i)-end_ind_left_row(i-1)...
                                                start_ind_right_row(i)-start_ind_right_row(i-1) end_ind_right_row(i)-end_ind_right_row(i-1)];
                                            
                                            
                                            last_row_start_left=start_ind_left_row(i-1);
                                            last_row_end_left=end_ind_left_row(i-1);
                                            
                                            last_row_start_right=start_ind_right_row(i-1);
                                            last_row_end_right=end_ind_right_row(i-1);
                                            
                                            temp_ind=i-1;
                                            
                                            if  end_ind_right_row(i)==0 && end_ind_right_row(temp_ind)~=0
                                                
                                                row_diff(4)=0;
                                                row_diff(3)=0;
                                                
                                                Sum_box_up=Sum_box_up-sum(sum(current_plane(start_ind_up_col(j):...
                                                    end_ind_up_col(j),start_ind_right_row(temp_ind):end_ind_right_row(temp_ind))));
                                                
                                            end
                                            
                                        end
                                        
                                        
                                        %Start left box
                                        
                                        %When both indices exist
                                        
                                        if row_diff(1)~=0
                                            
                                            Sum_box_up=Sum_box_up-sum(sum(current_plane(start_ind_up_col(j):...
                                                end_ind_up_col(j),last_row_start_left:start_ind_left_row(i)-1)));
                                            
                                        end
                                        
                                        %End left box
                                        
                                        if row_diff(2)~=0
                                            
                                            Sum_box_up=Sum_box_up+sum(sum(current_plane(start_ind_up_col(j):...
                                                end_ind_up_col(j),last_row_end_left+1:end_ind_left_row(i))));
                                            
                                        end
                                        
                                        if last_row_start_right~=0 && start_ind_right_row(i)~=0
                                            
                                            % Start right box
                                            
                                            if row_diff(3)~=0
                                                
                                                Sum_box_up=Sum_box_up-sum(sum(current_plane(start_ind_up_col(j):...
                                                    end_ind_up_col(j),last_row_start_right:start_ind_right_row(i)-1)));
                                                
                                            end
                                            
                                            % End right box
                                            
                                            if row_diff(4)~=0
                                                
                                                Sum_box_up=Sum_box_up+sum(sum(current_plane(start_ind_up_col(j):end_ind_up_col...
                                                    (j),last_row_end_right+1:end_ind_right_row(i))));
                                                
                                            end
                                            
                                        end
                                        
                                        Sum_AB(i,j)=Sum_AB(i,j)+Sum_box_up;
                                        
                                        
                                        if i==start_ind_low_col(j)
                                            
                                            first_sum_box_up=Sum_box_up;
                                            
                                        end
                                        
                                    end
                                    
                                    
                                    
                                end
                                
                                
                            end
                            
                            
                        end
                        
                    end
                    
                    if u~=1
                        
                        %Upper "rod" counts
                        
                        for j=start_ind_right_row(1):max(end_ind_right_row)
                            
                            
                            Sum_box_low2=first_sum_box_low2;
                            Sum_box_up2=first_sum_box_up2;
                            
                            
                            if start_ind_low_col(j)==0 && start_ind_low_col(j-1)~=0
                                
                                
                                Sum_box_low2=Sum_box_low2-sum(sum(current_plane(start_ind_low_col(j-1):...
                                    end_ind_low_col(j-1),start_ind_right_row(start_ind_up_col(j-1)):...
                                    end_ind_right_row(start_ind_up_col(j-1)))));
                                
                                if start_ind_left_row(start_ind_up_col(j-1))~=0
                                    
                                    Sum_box_low2=Sum_box_low2-sum(sum(current_plane(start_ind_low_col(j-1):...
                                        end_ind_low_col(j-1),start_ind_left_row(start_ind_up_col(j-1)):...
                                        end_ind_left_row(start_ind_up_col(j-1)))));
                                    
                                end
                                
                            end
                            
                            
                            if  start_ind_low_col(j)~=0
                                
                                new_i=start_ind_low_col(j);
                                
                                %%% MOVE SUMBOX TO NEW COLUMN COORDINATES %%%
                                
                                
                                if j~=start_ind_right_row(1)
                                    
                                    if start_ind_up_col(j-1)~=0
                                        
                                        last_i=start_ind_low_col(j-1);
                                        
                                        column_diff=[new_i-last_i end_ind_low_col(j)-end_ind_low_col(j-1)];
                                        
                                        if start_ind_right_row(start_ind_up_col(j-1))~=0 && start_ind_left_row(start_ind_up_col(j-1))
                                            
                                            if column_diff(1)~=0
                                                
                                                Sum_box_low2=Sum_box_low2-sum(sum(current_plane(last_i:...
                                                    new_i-1,start_ind_right_row(start_ind_up_col(j-1)):end_ind_right_row(start_ind_up_col(j-1)))))...
                                                    -sum(sum(current_plane(last_i:...
                                                    new_i-1,start_ind_left_row(start_ind_up_col(j-1)):end_ind_left_row(start_ind_up_col(j-1)))));
                                                
                                            end
                                            
                                            
                                            if column_diff(2)~=0
                                                
                                                Sum_box_low2=Sum_box_low2+sum(sum(current_plane(end_ind_low_col(j-1)+1:...
                                                    end_ind_low_col(j),start_ind_right_row(start_ind_up_col(j-1)):end_ind_right_row(start_ind_up_col(j-1)))))...
                                                    +sum(sum(current_plane(end_ind_low_col(j-1)+1:...
                                                    end_ind_low_col(j),start_ind_left_row(start_ind_up_col(j-1)):end_ind_left_row(start_ind_up_col(j-1)))));
                                                
                                            end
                                            
                                        else
                                            
                                            
                                            if start_ind_right_row(start_ind_up_col(j-1))~=0
                                                
                                                if column_diff(1)~=0
                                                    
                                                    Sum_box_low2=Sum_box_low2-sum(sum(current_plane(last_i:...
                                                        new_i-1,start_ind_right_row(start_ind_up_col(j-1)):end_ind_right_row(start_ind_up_col(j-1)))));
                                                    
                                                    
                                                    
                                                end
                                                
                                                if column_diff(2)~=0
                                                    
                                                    Sum_box_low2=Sum_box_low2+sum(sum(current_plane(end_ind_low_col(j-1)+1:...
                                                        end_ind_low_col(j),start_ind_right_row(start_ind_up_col(j-1)):end_ind_right_row(start_ind_up_col(j-1)))));
                                                    
                                                end
                                                
                                            elseif start_ind_left_row(start_ind_up_col(j-1))~=0
                                                
                                                if column_diff(1)~=0
                                                    
                                                    
                                                    Sum_box_low2=Sum_box_low2-sum(sum(current_plane(last_i:...
                                                        new_i-1,start_ind_left_row(start_ind_up_col(j-1)):end_ind_left_row(start_ind_up_col(j-1)))));
                                                    
                                                    
                                                end
                                                
                                                if column_diff(2)~=0
                                                    
                                                    Sum_box_low2=Sum_box_low2+sum(sum(current_plane(end_ind_low_col(j-1)+1:...
                                                        end_ind_low_col(j),start_ind_left_row(start_ind_up_col(j-1)):end_ind_left_row(last_i))));
                                                    
                                                    
                                                end
                                                
                                            end
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                                for i=start_ind_up_col(j):end_ind_up_col(j)
                                    
                                    if first_sum_low2==1
                                        
                                        if end_ind_left_row(i)~=0
                                            
                                            Sum_box_low2=sum(sum(current_plane(new_i:end_ind_low_col(j),start_ind_left_row(i)...
                                                :end_ind_left_row(i))))+sum(sum(current_plane(new_i:end_ind_low_col(j)...
                                                ,start_ind_right_row(i):end_ind_right_row(i))));
                                            
                                            
                                        else
                                            
                                            Sum_box_low2=Sum_box_low2+sum(sum(current_plane(new_i:end_ind_low_col(j)...
                                                ,start_ind_right_row(i):end_ind_right_row(i))));
                                            
                                        end
                                        
                                        
                                        first_sum_low2=0;
                                        first_sum_box_low2=Sum_box_low2;
                                        
                                        Sum_AB(i,j)=Sum_box_low2;
                                        
                                    else
                                        
                                        
                                        %%% MOVE SUMBOX TO NEW ROW COORDINATES %%%
                                        
                                        if i==start_ind_up_col(j)
                                            
                                            row_diff=[start_ind_left_row(i)-start_ind_left_row(start_ind_up_col(j-1)) ...
                                                end_ind_left_row(i)-end_ind_left_row(start_ind_up_col(j-1)) start_ind_right_row(i)-...
                                                start_ind_right_row(start_ind_up_col(j-1))...
                                                end_ind_right_row(i)-end_ind_right_row(start_ind_up_col(j-1))];
                                            
                                            last_row_start_left=start_ind_left_row(start_ind_up_col(j-1));
                                            last_row_end_left=end_ind_left_row(start_ind_up_col(j-1));
                                            
                                            last_row_start_right=start_ind_right_row(start_ind_up_col(j-1));
                                            last_row_end_right=end_ind_right_row(start_ind_up_col(j-1));
                                            
                                            temp_ind=start_ind_up_col(j-1);
                                            
                                        else
                                            
                                            row_diff=[start_ind_left_row(i)-start_ind_left_row(i-1) end_ind_left_row(i)-end_ind_left_row(i-1)...
                                                start_ind_right_row(i)-start_ind_right_row(i-1) end_ind_right_row(i)-end_ind_right_row(i-1)];
                                            
                                            
                                            last_row_start_left=start_ind_left_row(i-1);
                                            last_row_end_left=end_ind_left_row(i-1);
                                            
                                            last_row_start_right=start_ind_right_row(i-1);
                                            last_row_end_right=end_ind_right_row(i-1);
                                            
                                            temp_ind=i-1;
                                            
                                        end
                                        
                                        
                                        
                                        if  start_ind_left_row(i)~=0 && start_ind_left_row(temp_ind)==0
                                            
                                            row_diff(1)=0;
                                            row_diff(2)=0;
                                            
                                            Sum_box_low2=Sum_box_low2+sum(sum(current_plane(start_ind_low_col(j):...
                                                end_ind_low_col(j),start_ind_left_row(i):end_ind_left_row(i))));
                                            
                                        end
                                        
                                        
                                        %Start left box
                                        
                                        
                                        %When both indices exist
                                        
                                        if last_row_start_left~=0 && start_ind_left_row(i)~=0
                                            
                                            if row_diff(1)~=0
                                                
                                                Sum_box_low2=Sum_box_low2-sum(sum(current_plane(start_ind_low_col(j):...
                                                    end_ind_low_col(j),last_row_start_left:start_ind_left_row(i)-1)));
                                                
                                            end
                                            
                                            %End left box
                                            
                                            if row_diff(2)~=0
                                                
                                                if row_diff(2)>0
                                                    
                                                    Sum_box_low2=Sum_box_low2+sum(sum(current_plane(start_ind_low_col(j):...
                                                        end_ind_low_col(j),last_row_end_left+1:end_ind_left_row(i))));
                                                    
                                                else
                                                    
                                                    Sum_box_low2=Sum_box_low2-sum(sum(current_plane(start_ind_low_col(j):...
                                                        end_ind_low_col(j),end_ind_left_row(i)+1:last_row_end_left)));
                                                    
                                                end
                                                
                                            end
                                            
                                            % only 1 index exists
                                            %asd check
                                            
                                        elseif last_row_start_left~=0 && start_ind_left_row(i)==0
                                            
                                            Sum_box_low2=Sum_box_low2-sum(sum(current_plane(start_ind_low_col(j):...
                                                end_ind_low_col(j),last_row_start_left:last_row_end_left)));
                                            
                                        elseif start_ind_left_row(i)~=0 && last_row_start_left==0
                                            
                                            Sum_box_low2=Sum_box_low2+sum(sum(current_plane(start_ind_low_col(j):...
                                                end_ind_low_col(j),start_ind_left_row(i):end_ind_left_row(i))));
                                            
                                            
                                        end
                                        
                                        if last_row_start_right~=0 && start_ind_right_row(i)~=0
                                            
                                            % Start right box
                                            
                                            if row_diff(3)~=0
                                                
                                                if row_diff(3)>0
                                                    
                                                    Sum_box_low2=Sum_box_low2-sum(sum(current_plane(start_ind_low_col(j):...
                                                        end_ind_low_col(j),last_row_start_right:start_ind_right_row(i)-1)));
                                                    
                                                else
                                                    
                                                    Sum_box_low2=Sum_box_low2+sum(sum(current_plane(start_ind_low_col(j):...
                                                        end_ind_low_col(j),start_ind_right_row(i):last_row_start_right-1)));
                                                    
                                                end
                                                
                                            end
                                            
                                            % End right box
                                            
                                            if row_diff(4)~=0
                                                
                                                if row_diff(4)>0
                                                    
                                                    Sum_box_low2=Sum_box_low2+sum(sum(current_plane(start_ind_low_col(j):end_ind_low_col...
                                                        (j),last_row_end_right:end_ind_right_row(i))));
                                                    
                                                else
                                                    
                                                    Sum_box_low2=Sum_box_low2-sum(sum(current_plane(start_ind_low_col(j):end_ind_low_col...
                                                        (j),end_ind_right_row(i):last_row_end_right)));
                                                    
                                                end
                                                
                                            end
                                            
                                        elseif last_row_start_right~=0 && start_ind_right_row(i)==0
                                            
                                            Sum_box_low2=Sum_box_low2-sum(sum(current_plane(start_ind_low_col(j):end_ind_low_col...
                                                (j),last_row_start_right:last_row_end_right)));
                                            
                                        elseif start_ind_right_row(i)~=0 && last_row_start_right==0
                                            
                                            Sum_box_low2=Sum_box_low2+sum(sum(current_plane(start_ind_low_col(j):end_ind_low_col...
                                                (j),start_ind_right_row(i):end_ind_right_row(i))));
                                            
                                        end
                                        
                                        
                                        Sum_AB(i,j)=Sum_box_low2;
                                        
                                        
                                        if i==start_ind_up_col(j)
                                            
                                            first_sum_box_low2=Sum_box_low2;
                                            
                                        end
                                        
                                        
                                    end
                                    
                                    
                                end
                                
                                
                            end
                            
                            
                            %%% UPPER SUMBOX
                            
                            new_i_up=start_ind_up_col(j);
                            
                            %%% MOVE SUMBOX TO NEW COLUMN COORDINATES %%%
                            
                            if j~=start_ind_right_row(1) && start_ind_up_col(j-1)~=0
                                
                                last_i_up=start_ind_up_col(j-1);
                                
                                column_diff=[start_ind_up_col(j)-start_ind_up_col(j-1) end_ind_up_col(j)-end_ind_up_col(j-1)];
                                
                                if start_ind_right_row(last_i_up)~=0 && start_ind_left_row(last_i_up)
                                    
                                    if column_diff(1)~=0
                                        
                                        Sum_box_up2=Sum_box_up2-sum(sum(current_plane(last_i_up:...
                                            new_i_up-1,start_ind_right_row(last_i_up):end_ind_right_row(new_i_up))))...
                                            -sum(sum(current_plane(last_i_up:...
                                            new_i_up-1,start_ind_left_row(last_i_up):end_ind_left_row(new_i_up))));
                                        
                                    end
                                    
                                    if column_diff(2)~=0
                                        
                                        Sum_box_up2=Sum_box_up2+sum(sum(current_plane(end_ind_up_col(j-1)+1:...
                                            end_ind_up_col(j),start_ind_right_row(last_i_up):end_ind_right_row(last_i_up))))...
                                            +sum(sum(current_plane(end_ind_up_col(j-1)+1:...
                                            end_ind_up_col(j),start_ind_left_row(last_i_up):end_ind_left_row(last_i_up))));
                                        
                                    end
                                    
                                else
                                    
                                    
                                    if start_ind_right_row(last_i_up)~=0
                                        
                                        if column_diff(1)~=0
                                            
                                            Sum_box_up2=Sum_box_up2-sum(sum(current_plane(last_i_up:...
                                                new_i_up-1,start_ind_right_row(last_i_up):end_ind_right_row(last_i_up))));
                                            
                                        end
                                        
                                        if column_diff(2)~=0
                                            
                                            Sum_box_up2=Sum_box_up2+sum(sum(current_plane(end_ind_up_col(j-1)+1:...
                                                end_ind_up_col(j),start_ind_right_row(last_i_up):end_ind_right_row(last_i_up))));
                                            
                                        end
                                        
                                    elseif start_ind_left_row(last_i_up)~=0
                                        
                                        if column_diff(1)~=0
                                            
                                            
                                            Sum_box_up2=Sum_box_up2-sum(sum(current_plane(last_i_up:...
                                                new_i_up-1,start_ind_left_row(last_i_up):end_ind_left_row(last_i_up))));
                                            
                                        end
                                        
                                        if column_diff(2)~=0
                                            
                                            Sum_box_up2=Sum_box_up2+sum(sum(current_plane(end_ind_up_col(j-1)+1:...
                                                end_ind_up_col(j),start_ind_left_row(last_i_up):end_ind_left_row(last_i_up))));
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            for i=new_i_up:end_ind_up_col(j)
                                
                                if i~=0
                                    
                                    if first_sum_up2==1
                                        
                                        if end_ind_left_row(i)~=0
                                            
                                            Sum_box_up2=Sum_box_up2+sum(sum(current_plane(new_i_up:end_ind_up_col(j)...
                                                ,start_ind_left_row(i):end_ind_left_row(i))));
                                            Sum_box_up2=Sum_box_up2+sum(sum(current_plane(new_i_up:end_ind_up_col(j)...
                                                ,start_ind_right_row(i):end_ind_right_row(i))));
                                            
                                        else
                                            
                                            Sum_box_up2=Sum_box_up2+sum(sum(current_plane(new_i_up:end_ind_up_col(j)...
                                                ,start_ind_right_row(i):end_ind_right_row(i))));
                                            
                                        end
                                        
                                        first_sum_up2=0;
                                        first_sum_box_up2=Sum_box_up2;
                                        
                                        Sum_AB(i,j)=Sum_AB(i,j)+Sum_box_up2;
                                        
                                    else
                                        
                                        
                                        %%% MOVE SUMBOX TO NEW ROW COORDINATES %%%
                                        
                                        if i==new_i_up
                                            
                                            row_diff=[start_ind_left_row(start_ind_up_col(j))-start_ind_left_row(last_i_up) ...
                                                end_ind_left_row(start_ind_up_col(j))-end_ind_left_row(last_i_up) start_ind_right_row(start_ind_up_col(j))-...
                                                start_ind_right_row(last_i_up) end_ind_right_row(start_ind_up_col(j))-end_ind_right_row(last_i_up)];
                                            
                                            last_row_start_left=start_ind_left_row(last_i_up);
                                            last_row_end_left=end_ind_left_row(last_i_up);
                                            
                                            last_row_start_right=start_ind_right_row(last_i_up);
                                            last_row_end_right=end_ind_right_row(last_i_up);
                                            
                                            temp_ind=last_i_up;
                                            
                                            
                                        else
                                            
                                            row_diff=[start_ind_left_row(i)-start_ind_left_row(i-1) end_ind_left_row(i)-end_ind_left_row(i-1)...
                                                start_ind_right_row(i)-start_ind_right_row(i-1) end_ind_right_row(i)-end_ind_right_row(i-1)];
                                            
                                            
                                            last_row_start_left=start_ind_left_row(i-1);
                                            last_row_end_left=end_ind_left_row(i-1);
                                            
                                            last_row_start_right=start_ind_right_row(i-1);
                                            last_row_end_right=end_ind_right_row(i-1);
                                            
                                            temp_ind=i-1;
                                            
                                        end
                                        
                                        if  start_ind_left_row(start_ind_up_col(j))~=0 && start_ind_left_row(temp_ind)==0
                                            
                                            row_diff(1)=0;
                                            row_diff(2)=0;
                                            
                                            Sum_box_up2=Sum_box_up2+sum(sum(current_plane(start_ind_up_col(j):...
                                                end_ind_up_col(j),start_ind_left_row(start_ind_up_col(j)):end_ind_left_row(start_ind_up_col(j)))));
                                            
                                        end
                                        
                                        %Start left box
                                        
                                        
                                        %When both indices exist
                                        
                                        if last_row_start_left~=0 && start_ind_left_row(i)~=0
                                            
                                            if row_diff(1)~=0
                                                
                                                Sum_box_up2=Sum_box_up2-sum(sum(current_plane(start_ind_up_col(j):...
                                                    end_ind_up_col(j),last_row_start_left:start_ind_left_row(i)-1)));
                                                
                                            end
                                            
                                            %End left box
                                            
                                            if row_diff(2)~=0
                                                
                                                
                                                Sum_box_up2=Sum_box_up2+sum(sum(current_plane(start_ind_up_col(j):...
                                                    end_ind_up_col(j),last_row_end_left+1:end_ind_left_row(i))));
                                                
                                            end
                                            
                                            % only 1 index exists
                                            
                                        elseif last_row_start_left~=0 && start_ind_left_row(i)==0
                                            
                                            Sum_box_up2=Sum_box_up2-sum(sum(current_plane(start_ind_up_col(j):...
                                                end_ind_up_col(j),last_row_start_left:last_row_end_left)));
                                            
                                        elseif start_ind_left_row(i)~=0 && last_row_start_left==0
                                            
                                            Sum_box_up2=Sum_box_up2+sum(sum(current_plane(start_ind_up_col(j):...
                                                end_ind_up_col(j),start_ind_left_row(i):end_ind_left_row(i))));
                                            
                                            
                                        end
                                        
                                        if last_row_start_right~=0 && start_ind_right_row(i)~=0
                                            
                                            % Start right box
                                            
                                            if row_diff(3)~=0
                                                
                                                Sum_box_up2=Sum_box_up2-sum(sum(current_plane(start_ind_up_col(j):...
                                                    end_ind_up_col(j),last_row_start_right:start_ind_right_row(i)-1)));
                                                
                                            end
                                            
                                            % End right box
                                            
                                            if row_diff(4)~=0
                                                
                                                Sum_box_up2=Sum_box_up2+sum(sum(current_plane(start_ind_up_col(j):end_ind_up_col...
                                                    (j),last_row_end_right+1:end_ind_right_row(i))));
                                                
                                            end
                                            
                                        end
                                        
                                        
                                        Sum_AB(i,j)=Sum_AB(i,j)+Sum_box_up2;
                                        
                                        
                                        if i==new_i_up
                                            
                                            first_sum_box_up2=Sum_box_up2;
                                            
                                        end
                                        
                                    end
                                    
                                    
                                    
                                end
                                
                            end
                            
                        end
                        
                        
                        
                    end
                    
                    coeffs_AB=mean_detector_pairs./(mean_counts./(Sum_AB./elements));
                    coeffs_AB(isnan(coeffs_AB))=1;
                    coeffs_AB(coeffs_AB==inf)=1;
                    coeffs_AB(coeffs_AB==0)=1;
                    
                    det_coeffs(detectors_ring*(u+k-2)+1:detectors_ring*(u+k-1),detectors_ring*(k-1)+1:detectors_ring*k) = coeffs_AB;
                    
                end
                normalization = normalization .* det_coeffs;
                true_coincidences = true_coincidences .* det_coeffs;
                
            end
            if nargout >= 6
                varargout{6} = det_coeffs;
            end
            
            
        elseif options.normalization_options(2)==1
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% apply 3-D fansum to calculate effiency factors %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            elements_detectors = sum(true_det_indices(:,1) > 0);
            
            if r ~= inf && length(normalization_attenuation_correction)==2
                detector_counts = sum(true_coincidences ./ activity_coeffs) + sum(true_coincidences ./ activity_coeffs,2)';
            else
                detector_counts = sum(true_coincidences) + sum(true_coincidences,2)';
            end
            
            det_coeffs = mean(detector_counts./elements_detectors) ./ (detector_counts./elements_detectors);
            
            %Correct list-mode data and form coeff matrix
            
            true_coincidences = bsxfun(@times, true_coincidences, det_coeffs);
            true_coincidences = bsxfun(@times, true_coincidences, det_coeffs');
            normalization = bsxfun(@times, normalization, det_coeffs);
            normalization = bsxfun(@times, normalization, det_coeffs');
            if nargout >= 6
                varargout{6} = det_coeffs;
            end
            
            
            %Separate block profile coeffs from effiency coeffs
            if options.normalization_options(3)==1
                
                position_sum=zeros(1,max(det_position),'single');
                
                for u=1:rings
                    for i=1:detectors_ring
                        
                        current_detector=(u-1)*detectors_ring+i;
                        
                        position_sum(det_position(i))=position_sum(det_position(i))+det_coeffs(current_detector);
                        
                    end
                end
                
                position_coeffs=mean(position_sum)./position_sum;
                
                
                block_matrix = det_coeffs./repmat(position_coeffs(det_position), 1, rings);
                
                if nargout >= 5
                    varargout{5} = block_matrix;
                end
                
            end
            
        end
        
        
    end
    
    time=toc;
    if options.verbose
        disp(['Detector efficiency corrected (', num2str(time),'s)'])
    end
    
end

%% Transaxial geom. factors

if options.normalization_options(4)==1
    
    
    tic
    
    if ~options.use_raw_data
        
        Ns = Nang/cryst_per_block;
        
        tr_geom_matrix = zeros(options.Ndist, cryst_per_block, options.TotSinos,'single');
        for i = 1 : cryst_per_block
            k = mod(i,cryst_per_block);
            if k == 0
                k = cryst_per_block;
            end
            for phi = 1 : Ns
                if r ~= inf && length(normalization_attenuation_correction)==2
                    tr_geom_matrix(start_inds:end_inds,k,:) = tr_geom_matrix(start_inds:end_inds,k,:) + ...
                        SinDouble(start_inds:end_inds, k + (phi - 1)*cryst_per_block, :) ./ rad_coeff_matrix(start_inds:end_inds, k + (phi - 1)*cryst_per_block, :);
                else
                    tr_geom_matrix(start_inds:end_inds,k,:) = tr_geom_matrix(start_inds:end_inds,k,:) + SinDouble(start_inds:end_inds, k + (phi - 1)*cryst_per_block, :);
                end
            end
            tr_geom_matrix(:,k,:)= tr_geom_matrix(:,k,:) / Ns;
        end
        tr_geom_matrix(tr_geom_matrix == 0) = NaN;
        if r ~= inf
            tr_geom_matrix2 = tr_geom_matrix;
            tr_geom_matrix2(tr_geom_matrix2 < mean(tr_geom_matrix2,1)/2) = NaN;
            keskiarvo = nanmean(tr_geom_matrix2,1);
            tr_geom_matrix = bsxfun(@rdivide, keskiarvo, tr_geom_matrix);
            tr_geom_matrix(isinf(tr_geom_matrix(:))) = 1;
        else
            keskiarvo = nanmean(tr_geom_matrix,1);
            tr_geom_matrix = bsxfun(@rdivide, keskiarvo, tr_geom_matrix);
            tr_geom_matrix(isinf(tr_geom_matrix(:))) = 1;
        end
        if nanmean(tr_geom_matrix(:)) > 1.1
            tr_geom_matrix(tr_geom_matrix > nanmean(tr_geom_matrix(:))) = nanmean(tr_geom_matrix(:));
        end
        if r ~= inf
            tr_geom_matrix(isnan(tr_geom_matrix2(:))) = 1;
        else
            tr_geom_matrix(isnan(tr_geom_matrix(:))) = 1;
        end
        
        
        if nargout >= 7
            varargout{7} = tr_geom_matrix;
        end
        
        tr_geom_matrix = repmat(tr_geom_matrix, 1, Ns, 1);
        
        
        
        SinDouble=SinDouble.*tr_geom_matrix;
        
        normalization=normalization.*tr_geom_matrix;
        
    else
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% List-mode correction %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tr_geom_matrix = ones(koko, koko, 'single');
        
        for j = 1 : rings
            KalPa = zeros(detectors_ring, cryst_per_block,'single');
            ind = start_ind;
            for k = start_ind + (j - 1)*detectors_ring : end_ind + (j - 1)*detectors_ring
                if r ~= inf && length(normalization_attenuation_correction)==2
                    temp = diag(true_coincidences,-(k-1)) ./ diag(activity_coeffs,-(k-1));
                else
                    temp = diag(true_coincidences,-(k-1));
                end
                l = 0;
                for ll = 1 : cryst_per_block
                    KalPa(ind + l,ll) = sum(temp(ll:cryst_per_block:end));
                    l = l + 1;
                end
                ind = ind + 1;
            end
            KalPa(KalPa == 0) = NaN;
            if r ~= inf
                KalPa2 = KalPa;
                KalPa2(KalPa < nanmean(KalPa,1)/1.5) = NaN;
                keskiarvo = nanmean(KalPa2,1);
            else
                keskiarvo = nanmean(KalPa,1);
            end
            keskiarvo(isnan(keskiarvo)) = 1;
            apu = keskiarvo ./ KalPa;
            if nanmean(apu(:)) > 1.1
                apu(apu > nanmean(apu(:))) = nanmean(apu(:));
            end
            apu(isnan(apu)) = 1;
            l = start_ind;
            for k = start_ind + (j - 1)*detectors_ring : end_ind + (j - 1)*detectors_ring
                temp = diag(apu,-(l-1));
                temp = repmat(temp, options.blocks_per_ring, 1);
                temp = repmat(temp, rings, 1);
                tr_geom_matrix(k:detectors_ring*rings + 1 : koko * (koko - (detectors_ring * (j-1)) - l)) = ...
                    temp(1:length(k:detectors_ring*rings + 1 : koko * (koko - (detectors_ring * (j-1)) - l)));
                l = l + 1;
            end
        end
        
        true_coincidences = true_coincidences .* tr_geom_matrix;
        normalization = normalization .* tr_geom_matrix;
        
        if nargout >= 7
            varargout{7} = tr_geom_matrix;
        end
        
    end
    
    time=toc;
    if options.verbose
        disp(['Transaxial geometric correction done (', num2str(time),'s)'])
    end
    
    
end

if ~options.use_raw_data
    if mashing > 1
        %         normalization = cell2mat(arrayfun(@(i) mean(normalization(:,i:i+mashing-1,:),2),1:mashing:size(normalization,2)-mashing+1,'UniformOutput',false));
        %         if detectors_ring > options.det_per_ring
        %             normalization(gaps) = 0;
        %         end
        %         Nang = Nang / mashing;
    end
end


if options.verbose
    disp('Saving normalization data')
end

norm_components = options.normalization_options;

if options.use_raw_data
    norm_file = [folder options.machine_name '_normalization_listmode.mat'];
    normalization = ((normalization(tril(true(size(normalization)), 0))));
else
    norm_file = [folder options.machine_name '_normalization_' num2str(options.Ndist) 'x' num2str(Nang) '_span' num2str(options.span) '.mat'];
end
if exist('OCTAVE_VERSION','builtin') == 0
    save(norm_file, 'normalization','norm_components','-v7.3')
else
    save(norm_file, 'normalization','norm_components','-v7')
end


if options.verbose
    disp('Normalization matrix saved')
end


%% Return corrected data (convert to sparse if list-mode)

if nargout >= 1
    varargout{1} = normalization;
end

if nargout >= 2
    if ~options.use_raw_data
        varargout{2} = SinDouble;
    else
        %form sparse matrix from
        true_coincidences = sparse(double(true_coincidences(tril(true(size(true_coincidences)), 0))));
        
        varargout{2} = true_coincidences;
        
    end
    
end

disp('Normalization complete')