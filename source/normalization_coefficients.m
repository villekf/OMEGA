function [varargout] = normalization_coefficients(options)
%% Args

%options: Sinogram data, detector pair coordinates, ring heights, apparatus radius

%normalization_attenuation_correction: apply attenuation correction
%yes = [inner_cylinder_radius (cm) cylinder_attenuation_constant (cm^2/m)] | no = empty)
%If inner_cylinder_radius=inf, cylinder is assumed to cover entire FOV
%If left empty uniform illumination for each LOR is assumed

%options.normalization_scatter_correction: fit gaussian to scatter tail from cylinder
%normalization data (Too narrow scatter tail may lead to unaccurate fit).
%Not supported for list-mode data.

%options.normalization_options(1): apply axial geometric correction (yes = 1 | no = 0)

%options.normalization_options(2): apply detector effiency correction. (Fansum = 1 | SPC = 2 | no = 0)
%Fan_sum uses 3-D fansum method for both data types or SPC "single-plane
%Casey" method for list mode-data (SPC computationally more expensive). SPC is
%supposed to be used with FOV covering source
%Fansum version includes block profile correction. Using
%options.normalization_options(2)=2 with fansum uses block profile correction
%before detector effiency correction

%options.normalization_options(4): apply transaxial geometric correction for plane source data (yes = 1 | no = 0)
%With cylinders transaxial correction produces appropriate coeffs for LORs
%passing cylinder (LORs passing near cylinder edges can also be inaccurate)
%Accuracy of radial sorting can be adjusted in transaxial correction section

%options.normalization_options(3): apply block profile correction. If using fansum
%correction is nested with effiency correction
%If using SPC block profile correction is done separately

%TRANAXIAL AND BLOCK PROFILE CORRECTION GROUP SORTING CAN BE ADJUSTED IN
%THEIR SECTIONS

%% Returns
% Normalization matrix containing all coeffs
% Corrected sinogram
% Individual correction coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Anssi Manninen, Ville-Veikko Wettenhovi
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
folder = strrep(folder, 'source','mat-files/');
folder = strrep(folder, '\','/');

if ~isfield(options,'lor') && options.use_raw_data
    
    
    
    lor_file = [folder options.machine_name '_detector_locations_' num2str(options.Nx) 'x' num2str(options.Ny) 'x' num2str(options.Nz) '_raw.mat'];
    if exist(lor_file, 'file') == 2
        variableInfo = who('-file', lor_file);
        if options.implementation == 1 || options.implementation == 4
            if ismember('lor', variableInfo)
                true_detectors = loadStructFromFile(lor_file,'lor');
            else
                lor_pixel_count_prepass(options);
                true_detectors = loadStructFromFile(lor_file,'lor');
            end
        else
            if ismember('lor_opencl', variableInfo)
                load(lor_file,'lor_opencl')
            else
                lor_pixel_count_prepass(options);
                load(lor_file,'lor_opencl')
            end
            true_detectors = lor_opencl;
            clear lor_opencl
        end
    else
        lor_pixel_count_prepass(options);
        if options.implementation == 1 || options.implementation == 4
            true_detectors = loadStructFromFile(lor_file,'lor');
        else
            load(lor_file,'lor_opencl')
            true_detectors = lor_opencl;
            clear lor_opencl
        end
    end
    
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


uniform_source=0;

%if using uniform source for instance infinidesimal thickness cylinder use
%uniform_source=1;



varargout = cell(nargout,1);


%% Coordinates

if ~options.use_raw_data
    
    z=sinogram_coordinates_3D(options);
    z_length = options.rings * options.cr_pz;
    z_true = linspace(0, z_length, options.rings + 1)./10;
    %mm --> cm
    z=z./10;
    
    [~, ~, xp, yp] = detector_coordinates(options);
    [x,y]=sinogram_coordinates_2D(options, xp, yp);
    
    %mm --> cm
    x=x./10;
    y=y./10;
    
else
    
    %z
    z_length = options.rings * options.cr_pz;
    z = linspace(0, z_length, options.rings + 1);
    
    z=z./10;
    
end

%x and y detector coordinates
[detectors_x,detectors_y]=detector_coordinates(options);

%mm --> cm
detectors_x=detectors_x./10;
detectors_y=detectors_y./10;


detectors_ring=options.detectors/options.rings;
sino_amount=sum(options.segment_table);
segment_amount=length(options.segment_table);

%Ring radius

R=options.diameter/10/2; %cm

normalization_attenuation_correction = options.normalization_phantom_radius;

%Inner ring radius


if ~isempty(normalization_attenuation_correction)
    
    r=normalization_attenuation_correction(1); %cm
    
else
    r = inf;
    
end


%% Scale stacked data (when using sinograms)

if options.use_raw_data
    if (isfield(options, 'coincidences') == 0 && ~exist('coincidences','var')) && options.use_machine < 2
        if options.partitions == 1
            if options.use_ASCII && options.use_machine == 0
                load([options.machine_name '_measurements_' options.name '_static_raw_ASCII.mat'], 'coincidences')
            elseif options.use_LMF && options.use_machine == 0
                load([options.machine_name '_measurements_' options.name '_static_raw_LMF.mat'], 'coincidences')
            elseif options.use_root && options.use_machine == 0
                load([options.machine_name '_measurements_' options.name '_static_raw_root.mat'], 'coincidences')
            else
                load([options.machine_name '_measurements_' options.name '_static_raw_listmode.mat'], 'coincidences')
            end
        else
            if options.use_ASCII && options.use_machine == 0
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_ASCII.mat'], 'coincidences')
            elseif options.use_LMF && options.use_machine == 0
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_LMF.mat'], 'coincidences')
            elseif options.use_root && options.use_machine == 0
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_root.mat'], 'coincidences')
            else
                load([options.machine_name '_measurements_' options.name '_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_raw_listmode.mat'], 'coincidences')
            end
        end
        options.coincidences = coincidences;
    end
    if iscell(options.coincidences)
        I = find(options.coincidences{1});
    else
        I = find(options.coincidences);
    end
    clear coincidences
% Sinogram data
else
    if (~options.reconstruct_trues && ~options.reconstruct_scatter) || options.use_machine > 0
        if options.partitions == 1 && isfield(options, 'SinM') == 0
            if options.use_machine < 2
                if options.use_machine == 0
                    load([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                        num2str(options.TotSinos) '_span' num2str(options.span) '.mat'],'raw_SinM')
                elseif  options.use_machine == 1
                    load([options.machine_name '_' options.name '_sinograms_combined_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                        num2str(options.TotSinos) '_span' num2str(options.span) '_listmode.mat'],'raw_SinM')
                end
                options.SinM = raw_SinM;
                clear raw_SinM
            else
                load([options.machine_name '_' options.name '_sinogram_original_static_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' ...
                    num2str(options.TotSinos) '_span' num2str(options.span) '_machine_sinogram.mat'],'SinM')
                options.SinM = SinM;
                clear SinM
            end
        elseif isfield(options, 'SinM') == 0
            if options.use_machine < 2
                if options.use_machine == 0
                    load([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ '
                        num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                        num2str(options.span) '.mat'], 'raw_SinM')
                elseif  options.use_machine == 1
                    load([options.machine_name '_' options.name '_sinograms_combined_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                        num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                        num2str(options.span) '_listmode.mat'], 'raw_SinM')
                end
                options.SinM = raw_SinM;
                clear raw_SinM
            else
                load([options.machine_name '_' options.name '_sinograms_original_' num2str(options.partitions) 'timepoints_for_total_of_ ' ...
                    num2str(options.tot_time) 's_' num2str(options.Ndist) 'x' num2str(options.Nang) 'x' num2str(options.TotSinos) '_span' ...
                    num2str(options.span) '_machine_sinogram.mat'], 'SinM')
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
        options.SinM=cell2mat(struct2cell(options.SinM));
    end
    SinDouble=double(options.SinM);
    
    %Initialize normalization matrix
    
    norm_matrix=ones(size(options.SinM));
    
    %Amount of stacked planes in even and odd numbered plane sinograms
    
    if rem(floor(options.span/2),2)==0
        
        odds=ceil(options.span/2);
        evens=floor(options.span/2);
        
    else
        
        odds=floor(options.span/2);
        evens=ceil(options.span/2);
        
    end
    last_index=0;
    
    for t=1:segment_amount
        
        if t==1
            for i=2:options.segment_table(1)-1
                
                %Michelogram corners
                
                if i<ceil(options.span/2) || options.segment_table(1)-i<ceil(options.span/2)
                    
                    SinDouble(:,:,i)=SinDouble(:,:,i)/min(i,options.segment_table(1)-i+1);
                    norm_matrix(:,:,i)=norm_matrix(:,:,i)/min(i,options.segment_table(1)-i+1);
                    
                    %Michelogram mid-section
                else
                    
                    if rem(i,2)==0
                        SinDouble(:,:,i)=SinDouble(:,:,i)/evens;
                        norm_matrix(:,:,i)=norm_matrix(:,:,i)/evens;
                    else
                        SinDouble(:,:,i)=SinDouble(:,:,i)/odds;
                        norm_matrix(:,:,i)=norm_matrix(:,:,i)/odds;
                    end
                end
            end
            
            last_index=options.segment_table(1);
            
        else
            
            for i=last_index+3:last_index+options.segment_table(t)-2
                
                %Michelogram corners
                if i-last_index>min(evens,odds)*2 && (last_index+options.segment_table(t))-i+1>min(evens,odds)*2
                    
                    
                    if rem(i-last_index,2)==0
                        SinDouble(:,:,i)=SinDouble(:,:,i)/min(odds,evens);
                        norm_matrix(:,:,i)=norm_matrix(:,:,i)/min(odds,evens);
                        
                    else
                        SinDouble(:,:,i)=SinDouble(:,:,i)/max(odds,evens);
                        norm_matrix(:,:,i)=norm_matrix(:,:,i)/max(odds,evens);
                    end
                    
                else
                    %Michelogram mid-section
                    
                    SinDouble(:,:,i)=SinDouble(:,:,i)/min(ceil((i-last_index)/2),(ceil((last_index+options.segment_table(t)-i+1)/2)));
                    norm_matrix(:,:,i)=norm_matrix(:,:,i)/min(ceil(i/2),(ceil(last_index+options.segment_table(t)-i+1/2)));
                end
                
            end
            
            last_index=last_index+options.segment_table(t);
            
        end
    end
    
    time=toc;
    if options.verbose
        disp(['Spanned data scaled (', num2str(time),'s)'])
    end
    
end


%% Extended sparse list-mode matrix (when using list-mode data)

if options.use_raw_data
    
    tic
    
    %Convert coincidences to matrix from
    
    if isa(options.coincidences,'cell')
        Fov_coincidences=options.coincidences{1};
    else
        Fov_coincidences = options.coincidences;
    end
    
    
        koko = options.detectors;
        true_coincidences = zeros(koko, koko,'single');
        true_coincidences(tril(true(size(true_coincidences)), 0)) = full(Fov_coincidences);
%     cell_data{1}=[];
    
    %Extend sparse matrix
    
%     true_coincidences=zeros(options.detectors,options.detectors);
    true_det_indices=true_coincidences;
    index_list_mode=zeros(1,options.detectors+1);
    index_list_mode(1)=0;
    
    %Extend sparse matrix
    
%     Fov_coincidences=full(Fov_coincidences)';
    Fov_detectors=double(true_detectors);
    
    %Reshape coincidences
    
    for i=1:options.detectors
        
        index_list_mode(i+1)=index_list_mode(i)+options.detectors-(i-1);
        
%         true_coincidences(i:end,i)=Fov_coincidences(index_list_mode(i)+1:index_list_mode(i+1));
        
        true_det_indices(i:end,i)=Fov_detectors(index_list_mode(i)+1:index_list_mode(i+1));
        
    end
    
    %dispose of out of FOV counts
    
    true_det_indices(true_det_indices~=0)=1;
    true_coincidences=true_coincidences.*true_det_indices;
    
    %Initialize normalization matrix
    
    norm_matrix=true_det_indices;
    
    
    %Get starting indices of FOV detector pairs for each column (for speed up)
    %first one of mirror projections (ring diff=0 have only 1 projection)
    
    start_ind_low_col=zeros(1,detectors_ring);
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
    
    start_ind_up_col=zeros(1,detectors_ring);
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
    
    %First non-zero elements
    
    start_upper=find(start_ind_up_col,1);
    end_upper=detectors_ring+1-find(fliplr(start_ind_up_col),1);
    start_lower=find(start_ind_low_col,1);
    end_lower=detectors_ring+1-find(fliplr(start_ind_low_col),1);
    
    
    if options.normalization_options(2)==2
        
        %Get starting indices of FOV detector pairs for each row
        %first one of mirror projections (ring diff=0 have only 1 projection)
        
        start_ind_left_row=zeros(1,detectors_ring);
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
        
        start_ind_right_row=zeros(1,detectors_ring);
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
%     figure(5135)
%     subplot(1,2,1)
%     imagesc(true_coincidences(51*detectors_ring+1:52*detectors_ring,50*detectors_ring+1:51*detectors_ring))
    
    time=toc;
    if options.verbose
        disp(['Sparse matrices extended and false detector pairs disposed (', num2str(time),'s)'])
    end
end






%% Attenuation correction & scatter correction

if isempty(normalization_attenuation_correction)
    
    if options.use_raw_data==0
        
        true_radius = sqrt((x(:,1)-x(:,2)).^2+(y(:,1)-y(:,2)).^2);
        true_radius = reshape(true_radius, options.Ndist, options.Nang);
        
        radius_in_cylinder=true_radius;
        
    else
        
        true_radius=ones(detectors_ring);
        
    end
end

if ~isempty(normalization_attenuation_correction)
    
    if normalization_attenuation_correction(1)==inf && options.use_raw_data==0
        
        true_radius = sqrt((x(:,1)-x(:,2)).^2+(y(:,1)-y(:,2)).^2);
        true_radius = reshape(true_radius, options.Ndist, options.Nang);
        
        radius_in_cylinder=true_radius;
        
    elseif normalization_attenuation_correction(1)==inf && options.use_raw_data
        
        
        true_radius=ones(detectors_ring);
        
    end
    
    if r~=inf && ~options.use_raw_data
        
        %blur specs
        
        sigm = 2.0;
        %Window size
        sz = 4;
        x_gauss=-sz:sz;
        
        Exp_comp = -(x_gauss.^2)/(2*sigm^2);
        Kernel= exp(Exp_comp)/(sqrt(2*pi*sigm^2));
        
        if r~=inf
            
            radius_in_cylinder=zeros(1,options.Nang*options.Ndist);
            a=zeros(1,options.Nang*options.Ndist);
            c=a;
            intersects=zeros(options.Nang*options.Ndist,2);
            
            
            for i=1:options.Nang*options.Ndist
                
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
            true_radius=reshape(radius_in_cylinder,options.Ndist,options.Nang);
            
        end
        %min length in cylinder for each ring difference
        
        rad_min=min(radius_in_cylinder);
        
        
        
        % Out of cylinder indices
        
        cut_off_start=zeros(1,options.Nang);
        cut_off_end=cut_off_start;
        
        for j=1:options.Nang
            i=1;
            while i<=options.Ndist/2 && isnan(true_radius(i,j))
                
                i=i+1;
                
            end
            
            cut_off_start(j)=i;
            
        end
        
        
        for j=1:options.Nang
            i=options.Ndist/2;
            while i<=options.Ndist && ~isnan(true_radius(i,j))
                
                i=i+1;
                
            end
            cut_off_end(j)=i-1;
        end
        
    end
    
end

if (options.normalization_scatter_correction || ~isempty(normalization_attenuation_correction)) && ~options.use_raw_data
    
    if length(normalization_attenuation_correction)==2
        attenuation_coeff=normalization_attenuation_correction(2);  %mass attenuation times density
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Scatter correction for sinogram%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if options.normalization_scatter_correction
        
        if r~=inf
            
            tic
            
            % Transaxial scatter counts
            
            Scatters=zeros(1,options.Ndist);
            for j=1:options.Ndist
                
                if sum(isnan(true_radius(j,:)))==options.Nang
                    
                    Scatters(j)=sum(sum(SinDouble(j,:,:)));
                    
                end
                
            end
            
            Sino_counts=zeros(1,sino_amount);
            
            cut_pix=ceil(options.Ndist/100)+1;
            Scatters=[Scatters(1:cut_off_start-cut_pix) Scatters(cut_off_end+cut_pix:end)]'./sino_amount;
            
            %fit gaussian to scatter tail (closest transaxial sums to cylinder are cut off)
            
            sigm=zeros(1,1000);
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
                
                J_h(:,1)=const(i).*(xx-mu*ones(size(xx))).^2/(sigm(i)^3).*exp(-(xx-mu*ones(size(xx))).^2./(2*sigm(i)^2));
                
                J_h(:,2)=exp(-(xx-mu*ones(size(xx))).^2./(2*sigm(i)^2));
                
                h_i=const(i).*exp(-(xx-mu*ones(size(xx))).^2./(2*sigm(i)^2));
                
                new_coeffs=Gauss_Newton(Scatters,h_i,J_h,[sigm(i) ; const(i)],1);
                sigm(i+1)=new_coeffs(1);
                const(i+1)=new_coeffs(2);
                
                i=i+1;
                
            end
            
            values=const(i).*exp(-(yy'-mu*ones(size(yy))').^2./(2*sigm(i)^2));
            
            if i==1000
                
                error('Gaussian fit for scatter not converging')
                
            end
            
            for u=1:sino_amount
                
                
                if sum(isnan(true_radius(j,:)))==options.Nang
                    
                    Sino_counts(u)=sum(sum(SinDouble(1:cut_off_start-cut_pix,:,u)))+sum(sum(SinDouble(cut_off_end+cut_pix:end,:,u)));
                    
                end
                
            end
            
            Scaling_factors=Sino_counts./mean(Sino_counts);
            
            
            %export scattered counts
            Scatter_matrix=zeros(size(SinDouble));
            for u=1:sino_amount
                
                Sino_gauss=Scaling_factors(u).*values'./options.Nang;
                
                Scatter_matrix(:,:,u)=repmat(Sino_gauss', 1, options.Nang);
                SinDouble(:,:,u)=SinDouble(:,:,u)-Scatter_matrix(:,:,u);
                
            end
            
            SinDouble(SinDouble<0)=0;
            
            time=toc;
            if options.verbose
            disp(['Scatter correction done (', num2str(time),'s)'])
            end
            
        else
            
            warning("Scatter correction is not supported for non-cylinder source")
            
        end
        
    else
        
        disp("Normalization coefficients are calculated without scatter correction")
        
    end
    
    
%     if ~isempty(normalization_attenuation_correction)
%         
%         tic
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%  Attenuation for sinogram data (UNUSED)  %%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         activity_coeffs=zeros(size(SinDouble));
%         rad_coeff_matrix = zeros(size(SinDouble));
%         
%         h=zeros(1,sino_amount);
%         for u=1:segment_amount
%             
%             
%             %height from first ring of segment u
%             if u==1
%                 ring=1;
%             else
%                 ring=sum(options.segment_table(1:u-1))+1;
%             end
%             
%             
%             
%             if u~=1
%                 start_ind=sum(options.segment_table(1:(u-1)))+1;
%             else
%                 start_ind=1;
%             end
%             end_ind=sum(options.segment_table(1:u));
%             
%             
%             
%             h(u)=abs(z(ring,1)-z(ring,2))/10;
%             if h(u)~=0
%                 rad_profile_axial=sqrt(radius_in_cylinder.^2+(isnan(radius_in_cylinder).*h(u)).^2);
%             else
%                 rad_profile_axial=radius_in_cylinder;
%             end
%             
%             rad_profile_axial(isnan(rad_profile_axial))=rad_min;
%             rad_profile_axial(rad_profile_axial==0)=rad_min;
%             
%             if length(normalization_attenuation_correction)==2
%                 
%                 rad_coeffs=exp(-rad_min.*attenuation_coeff)./exp(-rad_profile_axial.*attenuation_coeff);
%                 
%                 %Blur attenuation coefficient sharp altering
%                 rad_coeffs=reshape(rad_coeffs,options.Ndist,options.Nang);
%                 I=[rad_coeffs; (ones(sz+1,options.Nang).*rad_coeffs(1:sz+1,:))];
%                 I=I(:);
%                 blurred=filter(Kernel,1,I');
% %                 blurred = I;
%                 
%                 new_coeffs=reshape(blurred',options.Ndist+sz+1,options.Nang);
%                 new_coeffs=new_coeffs(sz+2:end,:);
%                 
%                 
%                 rad_coeff_matrix(:,:,start_ind:end_ind)=repmat(new_coeffs,1,1,options.segment_table(u));
%                 
%                 
%                 %Correct Sinograms of segment u
%                 
%                 
%                 
%                 SinDouble(:,:,start_ind:end_ind)=SinDouble(:,:,start_ind:end_ind).*rad_coeff_matrix(:,:,start_ind:end_ind);
%             end
%             
%             %scale substance depth differences if applying transaxial correction (activity correction)
%             
%             radius_axial_mat=reshape(rad_profile_axial,options.Ndist,options.Nang);
%             
%             if r~=inf && uniform_source==0
%                 
%                 activity_coeffs(:,:,start_ind:end_ind)=max(max(radius_axial_mat))./repmat(radius_axial_mat,1,1,options.segment_table(u));
%                 activity_coeffs(isnan(true_radius))=1;
%                 
%                 
%                 
%                 cylinder_rad_cut=ceil((cut_off_end(1)-cut_off_start(1))/50);
%                 
%                 activity_coeffs(1:cut_off_start-cylinder_rad_cut,:,:)=0;
%                 
%                 activity_coeffs(cut_off_end+cylinder_rad_cut:end,:,:)=0;
%                 
%             end
%         end
%         
%     end
%     
%     
%     if length(normalization_attenuation_correction)==2
%         time=toc;
%         if options.verbose
%             disp(['Attenuation correction done (', num2str(time),'s)'])
%         end
%     end
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Attenuation correction for  list-mode data (UNUSED) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if ~isempty(normalization_attenuation_correction) && options.use_raw_data
%     
%     tic
%     if length(normalization_attenuation_correction)==2
%         attenuation_coeff=normalization_attenuation_correction(2);  %mass attenuation times density
%     end
%     %blur specs
%     
%     sigm = 2.0;
%     %Window size
%     sz = 4;
%     x_gauss=-sz:sz;
%     
%     Exp_comp = -(x_gauss.^2)/(2*sigm^2);
%     Kernel= exp(Exp_comp)/(sqrt(2*pi*sigm^2));
%     
%     if r~=inf
%         intersects=zeros(detectors_ring,detectors_ring,2);
%         a=zeros(detectors_ring);
%         c=a;
%         radius_in_cylinder=c;
%         
%         for j=1:detectors_ring
%             
%             if (start_ind_low_col(j)~=0 && end_ind_low_col(j)~=0) || (start_ind_up_col(j)~=0 && end_ind_up_col(j)~=0)
%                 
%                 for i=[start_ind_up_col(j):end_ind_up_col(j) start_ind_low_col(j):end_ind_low_col(j)]
%                     
%                     if i~=0
%                         
%                         %each LOR detector coordinates
%                         x_start=detectors_x(j);
%                         x_end=detectors_x(i);
%                         y_start=detectors_y(j);
%                         y_end=detectors_y(i);
%                         
%                         if (x_end-x_start)~=0
%                             a(i,j)=(y_end-y_start)/(x_end-x_start);
%                             if x_start<x_end
%                                 c(i,j)=y_start-a(i,j)*x_start;
%                             else
%                                 c(i,j)=y_end-a(i,j)*x_end;
%                             end
%                         else
%                             a(i,j)=NaN;
%                             c(i,j)=NaN;
%                         end
%                         
%                         %intersection points with inner cylinder
%                         if (x_end-x_start)~=0
%                             intersects(i,j,:)=roots([(a(i,j)^2+1) (2*a(i,j)*c(i,j)-2*R-2*R*a(i,j)) (c(i,j)^2+2*R^2-2*R*c(i,j)-r^2)]);
%                             
%                             if imag(intersects(i,j,1))<0.0001 && imag(intersects(i,j,2))<0.0001
%                                 radius_in_cylinder(i,j)=sqrt((intersects(i,j,1)-intersects(i,j,2))^2+((intersects(i,j,1)-intersects(i,j,2))*a(i,j))^2);
%                             else
%                                 radius_in_cylinder(i,j)=NaN;
%                             end
%                             
%                         else
%                             intersects(i,j,:)=roots([1 -2*R (x_start-R)^2-r^2+R^2]);
%                             if imag(intersects(i,j,1))<0.0001 && imag(intersects(i,j,2))<0.0001
%                                 radius_in_cylinder(i,j)=(intersects(i,j,1)-intersects(i,j,2));
%                             else
%                                 radius_in_cylinder(i,j)=NaN;
%                             end
%                             
%                         end
%                         
%                     end
%                     
%                 end
%                 
%             end
%             
%         end
%     end
%     
%     %Min length in cylinder for each ring difference
%     
%     rad_min=min(min(nonzeros(radius_in_cylinder)));
%     
%     true_radius=radius_in_cylinder;
%     
%     %
%     
%     if ~isempty(normalization_attenuation_correction) && (options.normalization_options(4)==1 || options.normalization_options(2)==2)
%         
%         rad_profile_axial=zeros(detectors_ring,detectors_ring,options.ring_difference+1);
%         h=zeros(1,options.ring_difference+1);
%         activity_coeffs=zeros(size(true_coincidences));
%         
%         for u=0:options.ring_difference
%             
%             %Height from first ring of segment u
%             
%             h(u+1)=abs(z(1)-z(1+u));
%             
%             if h(u+1)~=0
%                 
%                 rad_profile_axial(:,:,u+1)=sqrt(radius_in_cylinder.^2+((radius_in_cylinder~=0).*(~isnan(radius_in_cylinder)).*h(u+1)).^2);
%                 
%             else
%                 
%                 rad_profile_axial(:,:,u+1)=radius_in_cylinder;
%                 
%             end
%             
%             rad_profile_axial(isnan(rad_profile_axial))=rad_min;
%             rad_profile_axial(rad_profile_axial==0)=rad_min;
%             if length(normalization_attenuation_correction)==2
%                 
%                 rad_coeffs(:,:,u+1)=exp(-rad_min.*attenuation_coeff)./exp(-rad_profile_axial(:,:,u+1).*attenuation_coeff);
%                 
%                 
%                 %blur coeffs with gaussian window
%                 temp_vec=rad_coeffs(:,:,u+1);
%                 [sorted,indices]=sort(temp_vec(:));
%                 I=[sorted(1).*ones(sz*2+1,1) ; sorted];
%                 
%                 blurred=filter(Kernel,1,I');
%                 
%                 new_coeffs=blurred(sz*2+2:end);
%                 temp_vec2=zeros(1,length(new_coeffs));
%                 
%                 for i=1:length(new_coeffs)
%                     
%                     temp_vec2(indices(i))=new_coeffs(i);
%                     
%                 end
%                 
%                 new_coeffs_rad=reshape(temp_vec2,detectors_ring,detectors_ring);
%                 
%                 %Blur radial coeffs with
%                 
%                 for i=1:options.rings-u
%                     
%                     true_coincidences((u+i-1)*detectors_ring+1:(u+i)*detectors_ring,(i-1)*detectors_ring+1:i*detectors_ring)=new_coeffs_rad.*...
%                         true_coincidences((u+i-1)*detectors_ring+1:(u+i)*detectors_ring,(i-1)*detectors_ring+1:i*detectors_ring);
%                     
%                 end
%                 
%             end
%             if r~=inf && (options.normalization_options(4)==1 || options.normalization_options(2)==2) && uniform_source==0
%                 
%                 for q=1:options.rings-u
%                     
%                     activity_coeffs((u+q-1)*detectors_ring+1:(u+q)*detectors_ring,(q-1)...
%                         *detectors_ring+1:q*detectors_ring)=max(max(rad_profile_axial(:,:,u+1)))./rad_profile_axial(:,:,u+1);
%                     
%                 end
%                 
% %                 cylinder_rad_cut=ceil((cut_off_end(1)-cut_off_start(1))/50);
% %                 
% %                 activity_coeffs(1:cut_off_start-cylinder_rad_cut,:,:)=0;
% %                 
% %                 activity_coeffs(cut_off_end+cylinder_rad_cut:end,:,:)=0;
%                 
%             end
%             
%         end
%         
%         if length(normalization_attenuation_correction)==2
%             time=toc;
%             if options.verbose
%                 disp(['Attenuation correction done (', num2str(time),'s)'])
%             end
%         end
%         
%     end
    
    
% end

if options.use_raw_data
    
    if options.normalization_options(2)==2
        
        true_radius(isnan(true_radius))=0;
        lower_radius=tril(true_radius);
        
        start_ind_low_col_cylinder=zeros(1,detectors_ring);
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
        
        start_ind_up_col_cylinder=zeros(1,detectors_ring);
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
        start_ind_left_row_cylinder=zeros(1,detectors_ring);
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
        
        start_ind_right_row_cylinder=zeros(1,detectors_ring);
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


% if isempty(normalization_attenuation_correction)
%     
%     disp("Normalization coefficients are calculated without attenuation correction")
%     
% end


%% Axial block profiles and geom. factors for each ring

if options.normalization_options(1)==1
    tic
    
    if ~options.use_raw_data
        
        %%%Sinogram%%%
        
%         Sincounts=zeros(1,length(SinDouble(1,1,:)));
%         
%         for u=1:length(SinDouble(1,1,:))
%             
%             Sincounts(u)=sum(sum(SinDouble(:,:,u)));
%             
%         end
        
        Sincounts=sum(sum(SinDouble,2),1);
        
%         axial_geom_coeffs=zeros(1,1,sino_amount);
%         axial_block_profile=zeros(1,1,options.segment_table(1));
        
%         for u = 1 : options.segment_table(1)
%             axial_block_profile(u) = sqrt(mean(Sincounts)./Sincounts(u));
%         end
            axial_block_profile = sqrt(mean(Sincounts(1:options.segment_table(1)))./Sincounts(1:options.segment_table(1)));
        %         tic
%         for u = 1 : sino_amount
%             index1 = round(z(u,1), 8);
%             index2 = round(z(u,2), 8);
%             ring1 = (index1 == round(z(1:options.segment_table(1),1), 8));
%             ring2 = (index2 == round(z(1:options.segment_table(1),1), 8));
%             SinDouble(:,:,u) = SinDouble(:,:,u) * axial_block_profile(ring1) * axial_block_profile(ring2);
%             norm_matrix(:,:,u) = norm_matrix(:,:,u) * axial_block_profile(ring1) * axial_block_profile(ring2);
%         end
        %         toc
        %         tic
    if exist('OCTAVE_VERSION', 'builtin') == 0
        index1 = round((z(:,1) + (z(2,1) - z(1,1))) / (z(2,1) - z(1,1)), 8);
        index2 = round((z(:,2) + (z(2,2) - z(1,2))) / (z(2,2) - z(1,2)), 8);
    else
        index1 = round((z(:,1) + (z(2,1) - z(1,1))) / (z(2,1) - z(1,1)).*10e8)./10e8;
        index2 = round((z(:,2) + (z(2,2) - z(1,2))) / (z(2,2) - z(1,2)).*10e8)./10e8;
    end
        SinDouble = bsxfun(@times, SinDouble, axial_block_profile(index1) .* axial_block_profile(index2));
        norm_matrix = bsxfun(@times, norm_matrix, axial_block_profile(index1) .* axial_block_profile(index2));
        %             toc
        
%         Sincounts=zeros(1,length(SinDouble(1,1,:)));
%         
%         for u=1:length(SinDouble(1,1,:))
%             
%             Sincounts(u)=sum(sum(SinDouble(:,:,u)));
%             
%         end
%         Sincounts=sum(sum(SinDouble,2),1);
        
        
        %calculate coeffs (inverse ratio to mean counts)
            
            axial_geom_coeffs=mean(Sincounts)./Sincounts;
        
%         for u=1:sino_amount
%             
%             axial_geom_coeffs(u)=mean(Sincounts)/Sincounts(u);
%             
%             %correct sinogram
%             SinDouble(:,:,u)=axial_geom_coeffs(u).*SinDouble(:,:,u);
%             
%             %add coeffs to normalization matrix
%             norm_matrix(:,:,u)=axial_geom_coeffs(u).*norm_matrix(:,:,u);
%             
%         end
        SinDouble = bsxfun(@times, SinDouble, axial_geom_coeffs);
        norm_matrix = bsxfun(@times, norm_matrix, axial_geom_coeffs);
        
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% List Mode %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        tic
        
        plane_counts = zeros(options.rings,1);
        
        for u=1:options.rings
                if u==1
                    doubling=1;
                else
                    doubling=1;
                end
            
            plane_counts(u) = sum(sum(true_coincidences(detectors_ring*(u-1)+...
                    1:detectors_ring*(u),detectors_ring*(u-1)+1:detectors_ring*u)))*doubling;
            
        end
        
        axial_block_profile = sqrt(mean(plane_counts)./plane_counts);
        for u=1:options.rings
            norm_matrix(:,1 + detectors_ring*(u-1) : detectors_ring*u) = norm_matrix(:,1 + detectors_ring*(u-1) : detectors_ring*u) * axial_block_profile(u);
            true_coincidences(:,1 + detectors_ring*(u-1) : detectors_ring*u) = true_coincidences(:,1 + detectors_ring*(u-1) : detectors_ring*u) * axial_block_profile(u);
        end
        if nargout >= 4
            varargout{4} = axial_block_profile;
        end
        
        
        j=1:options.rings;
        %total for all planes
        plane_counts = sum(sum(true_coincidences)) + sum(sum(true_coincidences(detectors_ring*(j-1)+...
            1:detectors_ring*(j),detectors_ring*(j-1)+1:detectors_ring*j)));
        
%         for u=1:options.rings
%             
%             for j=1:options.rings-u+1
%                 if u==1
%                     doubling=2;
%                 else
%                     doubling=1;
%                 end
%                 
%                 plane_counts=plane_counts+sum(sum(true_coincidences(detectors_ring*(j+u-2)+...
%                     1:detectors_ring*(j+u-1),detectors_ring*(j-1)+1:detectors_ring*j)))*doubling;
%                 
%             end
%             
%         end
        
        plane_counts=plane_counts./(options.rings^2/2+options.rings/2);
        axial_coeffs=zeros(1,options.rings-1);
        
        for u=1:options.rings
            
            if u~=1
                doubling=1;
            else
                doubling=2;
            end
            
            for j=1:options.rings-u+1
                
                
                %Ring indices
                
                y_ring=detectors_ring*(j+u-2);
                y_ring2=detectors_ring*(j+u-1);
                x_ring=detectors_ring*(j-1);
                x_ring2=detectors_ring*j;
                
                %Coeffs for each plane
                
                axial_coeffs(u,j)=plane_counts./(sum(sum(true_coincidences(y_ring+1:y_ring2,x_ring+1:x_ring2)))*doubling);
                
                %Correct list-mode data axially
                
                true_coincidences(y_ring+1:y_ring2,x_ring+1:x_ring2)=axial_coeffs(u,j)*true_coincidences(y_ring+1:y_ring2,x_ring+1:x_ring2);
                
                %Coeffs to normalization matrix
                
                norm_matrix(y_ring+1:y_ring2,x_ring+1:x_ring2)=axial_coeffs(u,j)*norm_matrix(y_ring+1:y_ring2,x_ring+1:x_ring2);
                
            end
            
        end
        
        if nargout >= 3
            varargout{3} = axial_coeffs;
        end
        
        %Coeffs between cross projection planes
        
        
%         if nargout >= 4
            cross_coefs_save = zeros(options.rings-2, options.rings,2);
%         end
        
        for u=2:options.rings
            
            
            for j=1:options.rings-u+1
                
                Lower_sum=0;
                Upper_sum=0;
                
                y_ring=detectors_ring*(j+u-2);
                x_ring=detectors_ring*(j-1);
                
                %Counts for first cross plane
                
                for p=start_lower:end_lower
% p=start_lower:end_lower;

% Lower_sum = sum(sum(true_coincidences(y_ring+start_ind_low_col(p):y_ring+end_ind_low_col(p),x_ring+p)));
                    
                    Lower_sum=Lower_sum+sum(true_coincidences(y_ring+start_ind_low_col(p):y_ring+end_ind_low_col(p),x_ring+p));
                    
                end
                
                %Counts for second cross plane
                
                for p=start_upper:end_upper
% p=start_upper:end_upper;
                    
                    Upper_sum=Upper_sum+sum(true_coincidences(y_ring+start_ind_up_col(p):y_ring+end_ind_up_col(p),x_ring+p));
                    
                end
                %Calculate coeffs
                
                cross_coeffs=mean([Upper_sum Lower_sum])./[Upper_sum Lower_sum];
                
                cross_coefs_save(u-1,j,:) = cross_coeffs;
                
                %Correct cross planes and add coeffs to normalization matrix
                
                for p=start_lower:end_lower
                    
                    true_coincidences(y_ring+start_ind_low_col(p):y_ring+end_ind_low_col(p),x_ring+p)=...
                        true_coincidences(y_ring+start_ind_low_col(p):y_ring+end_ind_low_col(p),x_ring+p)*cross_coeffs(2);
                    
                    norm_matrix(y_ring+start_ind_low_col(p):y_ring+end_ind_low_col(p),x_ring+p)=...
                        norm_matrix(y_ring+start_ind_low_col(p):y_ring+end_ind_low_col(p),x_ring+p)*cross_coeffs(2);
                    
                end
                
                for p=start_upper:end_upper
                    
                    true_coincidences(y_ring+start_ind_up_col(p):y_ring+end_ind_up_col(p),x_ring+p)=...
                        true_coincidences(y_ring+start_ind_up_col(p):y_ring+end_ind_up_col(p),x_ring+p)*cross_coeffs(1);
                    
                    norm_matrix(y_ring+start_ind_up_col(p):y_ring+end_ind_up_col(p),x_ring+p)=...
                        norm_matrix(y_ring+start_ind_up_col(p):y_ring+end_ind_up_col(p),x_ring+p)*cross_coeffs(1);
                    
                end
                
            end
            
        end
%         if nargout >= 4
%             varargout{4} = cross_coefs_save;
%         end
        
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
    
    for i=1:options.detectors
        
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
    
    starting_block=options.cryst_per_block-i;
    current_block=starting_block;
    det_position=zeros(1,options.det_per_ring);
    
    %Determine detectors positions
    
    for i=1:length(detectors_x)
        
        if options.cryst_per_block/2<current_block
            
            det_position(i)=floor(norm(options.cryst_per_block/2-current_block));
            
        else
            
            det_position(i)=floor(norm(options.cryst_per_block/2-current_block))+1;
            
        end
        
        current_block=current_block+1;
        
        if current_block==options.cryst_per_block+1
            
            current_block=1;
            
        end
        
    end
    
    %Averaged data for each block/radial profile
    
    %Different block position combinations = (options.cryst_per_block/2)^2/2
    
    if ~options.use_raw_data
        
%         avg_profile=zeros(options.cryst_per_block/2,options.cryst_per_block/2,options.Ndist);
        avg_profile_LOR_amount = zeros(options.cryst_per_block/2,options.cryst_per_block/2,options.Ndist);
        not_found=0;
        detector_index_start=zeros(1,options.Nang*options.Ndist);
        detector_index_end=detector_index_start;
        index=zeros(options.Ndist*options.Nang,2);
        
        x = permute(reshape(x, options.Ndist, options.Nang, 2), [2 1 3]);
        x = reshape(x, options.Ndist*options.Nang, 2);
        y = permute(reshape(y, options.Ndist, options.Nang, 2), [2 1 3]);
        y = reshape(y, options.Ndist*options.Nang, 2);
        
        for j=1:options.Ndist*options.Nang
            
            p=1;
            t=0;
            
            while t~=2 &&  p<=length(det_position)
                
                if x(j,1)==detectors_x(p) && y(j,1)==detectors_y(p)
                    
                    detector_index_start(j)=det_position(p);
                    t=t+1;
                    
                end
                
                
                if x(j,2)==detectors_x(p) && y(j,2)==detectors_y(p)
                    
                    detector_index_end(j)=det_position(p);
                    t=t+1;
                    
                end
                p=p+1;
            end
            
            if t~=2
                not_found=not_found+1;
                warning('No corresponding coordinates found')
            end
            
            index(j,:)=[max(detector_index_start(j),detector_index_end(j))...
                min(detector_index_start(j),detector_index_end(j))];
            
            
            avg_profile_LOR_amount(index(j,1),index(j,2),ceil(j/options.Nang))=...
                avg_profile_LOR_amount(index(j,1),index(j,2),ceil(j/options.Nang)) + 1;
            
        end
        
        block_profile_coeffs=zeros(size(avg_profile_LOR_amount));
        
        
        
        %Calculate block profile coeffs by averaging over all sinograms%
        
        
        if options.normalization_options(3)==1 && options.normalization_options(2)==2
            
            block_profile_matrix=zeros(options.Nang,options.Ndist);
            
            avg_profile=zeros(size(avg_profile_LOR_amount));
            
            SinDouble = permute(SinDouble, [2 1 3]);
            
            for j=1:options.Ndist*options.Nang
                
                avg_profile(index(j,1),index(j,2),ceil(j/options.Nang))=avg_profile...
                    (index(j,1),index(j,2),ceil(j/options.Nang))+sum(SinDouble(j-floor((j-1)/options.Nang)...
                    *options.Nang,ceil(j/options.Nang),:));
                
            end
            
            
            
            for j=1:options.Ndist*options.Nang
                if avg_profile(index(j,1),index(j,2),ceil(j/options.Nang))==0
                    
                    avg_profile(index(j,1),index(j,2),ceil(j/options.Nang))=...
                        min(min(avg_profile(:,:,ceil(j/options.Nang))));
                    
                end
            end
            
            avg_profile=avg_profile./(avg_profile_LOR_amount.*sino_amount);
            
            for j=1:options.Ndist
                temp_mat=avg_profile(:,:,j);
                
                if nansum(nansum(avg_profile(:,:,j)))~=0
                    block_profile_coeffs(:,:,j)=nanmean(temp_mat(:))./avg_profile(:,:,j);
                    
                end
            end
%             index1=reshape(index(:,1),options.Ndist,options.Nang);
%             index2=reshape(index(:,2),options.Ndist,options.Nang);
            index1=reshape(index(:,1),options.Nang,options.Ndist);
            index2=reshape(index(:,2),options.Nang,options.Ndist);
            
            %coeff matrix for sinograms
            for p=1:options.Nang
                
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
            
            x = permute(reshape(x, options.Nang, options.Ndist, 2), [2 1 3]);
            x = reshape(x, options.Ndist*options.Nang, 2);
            y = permute(reshape(y, options.Nang, options.Ndist, 2), [2 1 3]);
            y = reshape(y, options.Ndist*options.Nang, 2);
            
            block_profile_matrix=repmat(block_profile_matrix,1,1,sino_amount);
            
            SinDouble=SinDouble.*block_profile_matrix;
            
            SinDouble = permute(SinDouble, [2 1 3]);
            block_profile_matrix = permute(block_profile_matrix, [2 1 3]);
            norm_matrix=norm_matrix.*block_profile_matrix;
            
            if nargout >= 5
                varargout{5} = block_profile_matrix;
            end
            
            time=toc;
            if options.verbose
                disp(['Transaxial block profile correction done (', num2str(time),'s)'])
            end
            
        end
%         
%         Ns = options.Nang / options.cryst_per_block;
%         
%         crystal_interference = zeros(options.Ndist, options.cryst_per_block, options.TotSinos);
%         for i = 1 : options.cryst_per_block
%             k = mod(i,options.cryst_per_block);
%             if k == 0
%                 k = options.cryst_per_block;
%             end
%             for phi = 1 : Ns
%                 crystal_interference(:,k,:) = crystal_interference(:,k,:) + SinDouble(:, k + (phi - 1)*options.cryst_per_block, :) * (1/Ns);
%             end
%         end
        
    else
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% List-mode form correction %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if options.normalization_options(2)==2 || options.normalization_options(4)==1 || options.normalization_options(3)==1
            
            radius_raw_data=zeros(detectors_ring);
            for j=1:detectors_ring
                
                if (start_ind_low_col(j)~=0 && end_ind_low_col(j)~=0) || (start_ind_up_col(j)~=0 && end_ind_up_col(j)~=0)
                    
                    for i=[start_ind_up_col(j):end_ind_up_col(j) start_ind_low_col(j):end_ind_low_col(j)]
                        
                        if i~=0
                            
                            radius_raw_data(i,j)=norm([detectors_x(j)-detectors_x(i) detectors_y(j)-detectors_y(i)]);
                            
                        end
                        
                    end
                    
                end
                
            end
            
            sorting_factor=6;
            
            max_rad=max(max(radius_raw_data));
            min_rad=min(min(nonzeros(radius_raw_data)));
            environment=(max_rad-min_rad)/(ceil(detectors_ring/sorting_factor)-1);
            
            %radial sorting
            radial_index=zeros(1,ceil(detectors_ring/sorting_factor)-1);
            
            for i=1:length(radial_index)
                
                radial_index(i)=min_rad+environment/2+(i-1)*environment;
                
            end
            
            radial_index=nonzeros(radial_index)';
            det_index=zeros(size(radius_raw_data));
            
            
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
            
            
            if options.normalization_options(2)==2
                
                profile_hits=zeros(max(det_position),max(det_position),max(max((det_index))));
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
                                    
                                    for p=start:options.rings
                                        
                                        for t=p:options.rings
                                            
                                            sec_ind=j+(p-1)*detectors_ring;
                                            
                                            profile_avg(y_ind,x_ind,z_ind)=...
                                                profile_avg(y_ind,x_ind,z_ind)+true_coincidences(i+(t-1)*detectors_ring,sec_ind);
                                            
                                        end
                                    end
                                    
                                    if start==2
                                        profile_hits(y_ind,x_ind,z_ind)...
                                            =profile_hits(y_ind,x_ind,z_ind)+(options.rings^2/2-options.rings/2);
                                    else
                                        profile_hits(y_ind,x_ind,z_ind)...
                                            =profile_hits(y_ind,x_ind,z_ind)+(options.rings^2/2+options.rings/2);
                                        
                                    end
                                    
                                end
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
                profile_avg=profile_avg./profile_hits;
                profile_avg(isnan(profile_avg))=0;
                %Mean counts for different radial groups
                mean_profile=zeros(1,max(max(det_index)));
                block_profile_coeffs=zeros(size(profile_avg));
                for u=1:max(max(det_index))
                    
                    mean_profile(u)=sum(sum(profile_avg(:,:,u)))./numel(nonzeros(profile_avg(:,:,u)));
                    
                    %coeffs
                    
                    block_profile_coeffs(:,:,u)=mean_profile(u)./profile_avg(:,:,u);
                    
                end
                block_profile_coeffs(block_profile_coeffs==inf)=0;
                block_profile_coeffs(isnan(block_profile_coeffs))=0;
                
                %correct list mode data
                
                block_matrix=zeros(size(true_coincidences));
                
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
                                
                                for p=1:options.rings
                                    
                                    sec_ind=j+(p-1)*detectors_ring;
                                    
                                    for t=p:options.rings
                                        
                                        true_coincidences(i+(t-1)*detectors_ring,sec_ind)=true_coincidences(i+(t-1)*detectors_ring,sec_ind)...
                                            *current_coeff;
                                        
                                        norm_matrix(i+(t-1)*detectors_ring,sec_ind)=norm_matrix(i+(t-1)*detectors_ring,sec_ind)...
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
        
%         if options.normalization_options(2)==2
            
%         end
        
    end
    
end

%% Detector effiency factors

if options.normalization_options(2)~=0
    
    tic
    
    if ~options.use_raw_data
        
        
        %  with 3-D fan-sum algorithm
        
        
        %effiency for LOR(i,j) = detector_effiency(i) x detector_effiency(j)
        
        z_rings=z(1:options.segment_table(1),1);
        det_num=zeros(options.Ndist,options.Nang,2);
        counts_det=zeros(detectors_ring,options.segment_table(1));
        coeff_matrix=zeros(size(SinDouble));
        
        %Determine each LORS detector numbers
        
        for u=1:options.Ndist*options.Nang
            
            found_it=0;
            t=1;
            while found_it==0
                
                if norm([x(u,1)-detectors_x(t) y(u,1)-detectors_y(t)])<0.001
                    
                    det_num(floor((u-1)/options.Nang)+1,u-floor((u-1)/options.Nang)*options.Nang,1)=t;
                    found_it=1;
                    
                end
                t=t+1;
            end
            found_it=0;
            t=1;
            while found_it==0
                
                if norm([x(u,2)-detectors_x(t) y(u,2)-detectors_y(t)])<0.001
                    
                    det_num(floor((u-1)/options.Nang)+1,u-floor((u-1)/options.Nang)*options.Nang,2)=t;
                    found_it=1;
                    
                end
                t=t+1;
            end
            
        end
        
        coeffs_detectors=zeros(detectors_ring,options.segment_table(1));
        hits_det=zeros(detectors_ring,options.segment_table(1));
        
        for k=1:sino_amount
            
            found1=0;
            found2=0;
            t=1;
            while found1==0
                
                if norm(z_rings(t)-z(k,1))<0.001
                    
                    ring_1=t;
                    found1=1;
                    
                end
                
                t=t+1;
                
            end
            
            t=1;
            
            while found2==0
                
                if norm(z_rings(t)-z(k,2))<0.001
                    
                    ring_2=t;
                    found2=1;
                    
                end
                
                t=t+1;
            end
            
            
            
            for i=1:options.Ndist
                
                for u=1:options.Nang
                    
                    counts_det(det_num(i,u,1),ring_1)=counts_det(det_num(i,u,1),ring_1)+SinDouble(i,u,k);
                    counts_det(det_num(i,u,2),ring_2)=counts_det(det_num(i,u,2),ring_2)+SinDouble(i,u,k);
                    
                    hits_det(det_num(i,u,1),ring_1)=hits_det(det_num(i,u,1),ring_1)+1;
                    hits_det(det_num(i,u,2),ring_2)=hits_det(det_num(i,u,2),ring_2)+1;
                    
                end
                
                
            end
            
        end
        
        %Calculate coeffs (mean inverse)
        
        for k=1:options.segment_table(1)
            
            counts_det(:,k)=counts_det(:,k)./hits_det(:,k);
            coeffs_detectors(:,k)=mean(counts_det(:,k))./counts_det(:,k);
            
        end
        
        
        if options.normalization_options(2)~=2
            
            block_profile_counts=zeros(1,max(detector_index_start));
            block_profile_hits=block_profile_counts;
            detector_index_start=reshape(detector_index_start,options.Ndist,options.Nang);
            detector_index_end=reshape(detector_index_end,options.Ndist,options.Nang);
            det_number_block=zeros(1,detectors_ring);
        end
        for k=1:sino_amount
            
            for i=1:options.Ndist
                
                for j=1:options.Nang
                    
                    coeff_matrix(i,j,k)=coeffs_detectors(det_num(i,j,1),ring_1)*coeffs_detectors(det_num(i,j,2),ring_2);
                    
                    if options.normalization_options(2)~=2
                        
                        block_profile_counts(detector_index_start(i,j))=block_profile_counts(detector_index_start(i,j))+coeffs_detectors(det_num(i,j,1),ring_1);
                        block_profile_counts(detector_index_end(i,j))=block_profile_counts(detector_index_end(i,j))+coeffs_detectors(det_num(i,j,2),ring_2);
                        block_profile_hits(detector_index_start(i,j))=block_profile_hits(detector_index_start(i,j))+1;
                        block_profile_hits(detector_index_end(i,j))=block_profile_hits(detector_index_end(i,j))+1;
                        
                    end
                end
                
            end
            
        end
        if options.normalization_options(2)~=2
            
            p=1;
            
            for k=1:detectors_ring
                
                for i=1:options.Ndist
                    
                    for j=1:options.Nang
                        
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
        norm_matrix=norm_matrix.*coeff_matrix(1:options.Ndist,:,:);
        
        if nargout >= 6
            varargout{6} = coeff_matrix;
        end
        
        if options.normalization_options(3) == 1 && options.normalization_options(2)~=2
            block_profile_counts=block_profile_counts./block_profile_hits;
            block_profile_coeffs=mean(block_profile_counts)./block_profile_counts;
            %Separate block profile coeffs from effiency coeffs
            detector_coeffs=zeros(detectors_ring,options.rings);
            compared=(z_rings(2)-z_rings(1))/10;
            
            for t=1:options.rings
                k=1;
                while abs(z_true(t)-z_rings(k))>compared
                    k=k+1;
                end
                
                for p=1:detectors_ring
                    
                    detector_coeffs(p,t)=coeffs_detectors(p,k)/block_profile_coeffs(det_number_block(p));
                    
                end
                
            end
            
            if nargout >= 5
                varargout{5} = detector_coeffs;
            end
        end
        
    else
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% List mode correction %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if isempty(normalization_attenuation_correction)
            
            true_radius=ones(detectors_ring);
            
        else
            
            if r~=inf
                
                if r==inf
                    true_radius=ones(detectors_ring);
                end
                
            end
            
        end
        
        if options.normalization_options(2)==2
            
            %Detector effiency correction with SP-C method (https://doi.org/10.1088/0031-9155/43/1/012, page: 193)
            
            elements_x_ax=zeros(1,detectors_ring);
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
            
            elements=zeros(detectors_ring);
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
            
            mean_detector_pairs=zeros(size(current_plane));
            
            effiencies_y_mean=zeros(1,detectors_ring);
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
            
            det_coeffs = zeros(size(norm_matrix));
            
            for u=1:options.ring_difference+1
                
                for k=1:options.rings-u+1
                    
                    Sum_AB=zeros(detectors_ring);
                    
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
                    
                    true_coincidences(detectors_ring*(u+k-2)+1:detectors_ring*(u+k-1),detectors_ring*(k-1)+1:detectors_ring*k)=...
                        true_coincidences(detectors_ring*(u+k-2)+1:detectors_ring*(u+k-1),detectors_ring*(k-1)+1:detectors_ring*k).*coeffs_AB;
                    norm_matrix(detectors_ring*(u+k-2)+1:detectors_ring*(u+k-1),detectors_ring*(k-1)+1:detectors_ring*k)=...
                        norm_matrix(detectors_ring*(u+k-2)+1:detectors_ring*(u+k-1),detectors_ring*(k-1)+1:detectors_ring*k).*coeffs_AB;
                    
                    det_coeffs(detectors_ring*(u+k-2)+1:detectors_ring*(u+k-1),detectors_ring*(k-1)+1:detectors_ring*k) = coeffs_AB;
                    
                end
                
            end
            if nargout >= 6
                varargout{6} = det_coeffs;
            end
            
%             figure(5135)
%             subplot(1,2,2)
%             imagesc(true_coincidences(51*detectors_ring+1:52*detectors_ring,50*detectors_ring+1:51*detectors_ring))
            
            
        end
        
        if options.normalization_options(2)==1
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%apply 3-D fansum to calculate effiency factors
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            true_radius(isnan(true_radius))=0;
            true_radius(true_radius~=0)=1;
            
            
            for u=1:options.ring_difference+1
                
                for k=1:options.rings-u+1
                    
                    cylinder_counts(detectors_ring*(u+k-2)+1:detectors_ring*(u+k-1),detectors_ring*(k-1)+1:detectors_ring*k)=true_coincidences(detectors_ring*(u+k-2)+1:detectors_ring*(u+k-1),detectors_ring*(k-1)+1:detectors_ring*k).*true_radius;
                    
                end
                
            end
            
            
            detector_counts=zeros(1,detectors_ring);
            det_coeffs=detector_counts;
            elements_detectors=numel(nonzeros(true_radius(:,1)))*options.rings;
            
            for u=1:options.rings
                
                det_indices=(u-1)*detectors_ring+1:u*detectors_ring;
                
                for i=1:detectors_ring
                    
                    current_detector=(u-1)*detectors_ring+i;
                    
                    detector_counts(current_detector)=sum(cylinder_counts(:,current_detector))...
                        +sum(cylinder_counts(current_detector,:));
                    
                end
                
                det_coeffs(det_indices)=mean(detector_counts(det_indices)./elements_detectors)...
                    ./(detector_counts(det_indices)./elements_detectors);
                
                
            end
            
            %Correct list-mode data and form coeff matrix
            
            for u=1:options.rings
                
                start_ind_col=detectors_ring*(u-1)+1;
                end_ind_row=u*detectors_ring;
                
                for i=1:detectors_ring
                    
                    
                    current_detector=i+(u-1)*detectors_ring;
                    
                    true_coincidences(start_ind_col:end,current_detector)=true_coincidences...
                        (start_ind_col:end,current_detector)*det_coeffs(current_detector);
                    
                    
                    norm_matrix(start_ind_col:end,current_detector)=norm_matrix(start_ind_col:...
                        end,current_detector)*det_coeffs(current_detector);
                    
                    
                    true_coincidences(current_detector,1:end_ind_row)=...
                        true_coincidences(current_detector,1:end_ind_row)*det_coeffs(current_detector);
                    
                    norm_matrix(current_detector,1:end_ind_row)=...
                        norm_matrix(current_detector,1:end_ind_row)*det_coeffs(current_detector);
                    
                    
                end
                
            end
            if nargout >= 6
                varargout{6} = det_coeffs;
            end
            
            
            %Separate block profile coeffs from effiency coeffs
            if options.normalization_options(3)==1
                
                position_sum=zeros(1,max(det_position));
                
                for u=1:options.rings
                    
                    for i=1:detectors_ring
                        
                        current_detector=(u-1)*detectors_ring+i;
                        
                        position_sum(det_position(i))=position_sum(det_position(i))+det_coeffs(current_detector);
                        
                    end
                    
                end
                
                position_coeffs=mean(position_sum)./position_sum;
                
                t=1;
                for i=1:options.detectors
                    
                    
                    det_coeffs(i)=det_coeffs(i)/position_coeffs(det_position(t));
                    
                    t=t+1;
                    
                    
                    if t>detectors_ring
                        t=1;
                    end
                    
                end
                
                if nargout >= 5
                    varargout{5} = det_coeffs;
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
        
        
        %scale activities to same for each radial profile (if using cylinder data)
%         if r~=inf && uniform_source==0
%             
%             Temp_sino=SinDouble.*activity_coeffs;
%             
%         else
            
%             Temp_sino=SinDouble;
%             
% %         end
%         % Radial profiles
%         
%         radius=zeros(options.Nang*options.Ndist,1);
%         
%         
%         %All radial lenghts
%         for i=1:options.Ndist*options.Nang
%             if ~isnan(true_radius(i))
%                 
%                 radius(i)=sqrt((x(i,1)-x(i,2))^2+(y(i,1)-y(i,2))^2);
%                 
%             end
%         end
%         
%         %Form radial groups
%         
%         %%%%%%%%%%% Group fining can be adjusted here if needed %%%%%%%%%%%%%%%%%%%
%         
%         %increasing factor increases groups
%         fining_factor=options.cryst_per_block;
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         
%         maximum=max(radius);
%         minimum=min(nonzeros(radius));
%         environment=norm(maximum-minimum)/(fining_factor*options.Ndist-1);
%         radial_index=zeros(fining_factor*options.Ndist-1,1);
%         
%         %Group radiusses
%         
%         for i=1:length(radial_index)
%             
%             radial_index(i)=minimum+environment/2+(i-1)*environment;
%             
%         end
%         
%         radial_index=nonzeros(radial_index);
%         radial_times=zeros(length(radial_index),1);
%         
%         
%         %Check if cylinder covers entire FOV
%         
%         if ~r~=inf
%             
%             true_radius=ones(options.Ndist*options.Nang,1);
%             
%         end
%         
%         %Sort LORs by radial profile
%         LOR_index=zeros(options.Nang*options.Ndist,1);
%         
%         for i=1:options.Ndist*options.Nang
%             
%             if ~isnan(true_radius(i))
%                 
%                 p=1;
%                 while norm(radius(i)-radial_index(p))>environment/2+10^-8 && radial_index(p)~=0
%                     p=p+1;
%                 end
%                 LOR_index(i)=p;
%                 
%                 radial_times(p)=radial_times(p)+1;
%                 
%             end
%             
%         end
%         
%         avg_counts_transaxial=zeros(1,max(LOR_index));
%         tr_geom_matrix=zeros(options.Ndist,options.Nang);
%         
%         %Calculate average counts for each radial group
%         
% %         LOR_index = reshape(LOR_index, options.Ndist, options.Nang)';
% %         LOR_index = LOR_index(:);
%         
%         for j=1:options.Ndist*options.Nang
%             
%             if LOR_index(j)~=0
%                 
%                 avg_counts_transaxial(LOR_index(j))=avg_counts_transaxial(LOR_index(j))+...
%                     (Temp_sino(ceil(j/options.Nang),j-floor((j-1)/options.Nang)*options.Nang,:));
%                 
% %                 if ceil(j/options.Nang) == 99 && j-floor((j-1)/options.Nang)*options.Nang == 9
% %                     return
% %                 end
%                 
%             end
%             
%         end
%         
%         avg_counts_transaxial=avg_counts_transaxial./(radial_times'*sino_amount);
% %         avg_counts_transaxial=avg_counts_transaxial./(sino_amount);
%         
%         %Get rid of undetermined coeffs
%         
%         avg_counts_transaxial(isnan(avg_counts_transaxial))=0;
%         avg_counts_transaxial(avg_counts_transaxial==inf)=0;
%         non_zeros=nonzeros(avg_counts_transaxial);
%         
%         %transaxial coeffs
%         
%         tr_geom_coeffs=(sum(non_zeros)/length(non_zeros))./(avg_counts_transaxial);
%         tr_geom_coeffs(tr_geom_coeffs==inf)=0;
%         tr_geom_coeffs(isnan(tr_geom_coeffs))=0;
%         
%         %Form matrix from coeffs
%         
%         for i=1:options.Ndist*options.Nang
%             
%             if LOR_index(i)~=0
%                 
%                 tr_geom_matrix(ceil(i/options.Nang),i-floor((i-1)/options.Nang)*options.Nang)=tr_geom_coeffs(LOR_index(i));
%                 
%             end
%         end
        
%         testi1 = mean(SinDouble,2);
%         testi2 = mean(testi1,1);
%         testi3 = 1./bsxfun(@rdivide, testi1, testi2);
%         testi4 = reshape(testi3, options.Ndist, options.TotSinos);
        
        %Cut columns near edges of cylinder when cylinder not covering entire FOV
        
%         if ~isempty(normalization_attenuation_correction) && normalization_attenuation_correction(1)~=inf
%             
% %             cylinder_rad_cut=ceil((cut_off_end(1)-cut_off_start(1))/50);
% %             cylinder_rad_cut = 3;
%             
%             min_coeff=min(min(tr_geom_matrix(cut_off_start(1):cut_off_end,:)));
%             
%             tr_geom_matrix(tr_geom_matrix == 0) = min_coeff;
%             
% %             tr_geom_matrix(1:cut_off_start-cylinder_rad_cut,:)=min_coeff;
% %             
% %             tr_geom_matrix(cut_off_end+cylinder_rad_cut:end,:)=min_coeff;
%             
%             
% %             min_coeff=min(min(testi3(cut_off_start(1):cut_off_end,:)));
% %             
% %             testi3(1:cut_off_start-cylinder_rad_cut,:,:)=min_coeff;
% %             
% %             testi3(cut_off_end+cylinder_rad_cut:end,:,:)=min_coeff;
%             
%             tr_geom_matrix(tr_geom_matrix <= 0) = min_coeff;
%             
%             tr_geom_matrix(tr_geom_matrix > mean(abs(nonzeros(tr_geom_matrix(:))))) = mean(abs(nonzeros(tr_geom_matrix(:))));
%             
%         end
%         keskiarvo = mean(mean(mean(tr_geom_matrix,3)));
%         tr_geom_matrix(tr_geom_matrix > keskiarvo * 2) = keskiarvo*2;
%         testi = diff(testi3);
% %         tr_geom_matrix3 = tr_geom_matrix;
%         testi3(testi3 > mean(abs((testi(:)))) / 10) = mean(abs((testi(:))) /10);
%         testi3(testi3 > mean(abs(nonzeros(testi3(:))))) = mean(abs(nonzeros(testi3(:))));
%         tr_geom_matrix = testi3;
%         tr_geom_matrix=repmat(testi3,1,options.Nang,1);
        
%         testi = diff(tr_geom_matrix);
%         testii = diff(flipud(tr_geom_matrix));
%         testi = [zeros(1, size(testi,2)) ; testi];
%         testi(end/2:end,:) = testii(round(end/2 - 1):end,:);
        
%         tr_geom_matrix(testi > mean(abs(nonzeros(testi(:))))) = mean(abs(nonzeros(tr_geom_matrix(:))));
        
        Ns = options.Nang/options.cryst_per_block;
        
        tr_geom_matrix = zeros(options.Ndist, options.cryst_per_block, options.TotSinos);
        for i = 1 : options.cryst_per_block
            k = mod(i,options.cryst_per_block);
            if k == 0
                k = options.cryst_per_block;
            end
            for phi = 1 : Ns
                tr_geom_matrix(:,k,:) = tr_geom_matrix(:,k,:) + SinDouble(:, k + (phi - 1)*options.cryst_per_block, :);
            end
            tr_geom_matrix(:,k,:)= tr_geom_matrix(:,k,:) / Ns;
        end
        keskiarvo = mean(tr_geom_matrix,1);
%         ind = tr_geom_matrix < (mean(keskiarvo(:)) * .1);
%         minimi = min(tr_geom_matrix(~ind));
        tr_geom_matrix = bsxfun(@rdivide, keskiarvo, tr_geom_matrix);
%         tr_geom_matrix(ind) = 1 / minimi;
        
%         if ~isempty(normalization_attenuation_correction) && normalization_attenuation_correction(1)~=inf
            
%             cylinder_rad_cut=ceil((cut_off_end(1)-cut_off_start(1))/50);
%             cylinder_rad_cut = 3;
            
%             min_coeff=min(min(min(tr_geom_matrix(cut_off_start(1):cut_off_end,:,:))));
%             
%             tr_geom_matrix(1:cut_off_start(1)-1,:,:) = min_coeff;
%             tr_geom_matrix(cut_off_end+1:end,:,:) = min_coeff;
            
%             tr_geom_matrix(tr_geom_matrix == 0) = min_coeff;
            
%             tr_geom_matrix(1:cut_off_start-cylinder_rad_cut,:)=min_coeff;
            
%             tr_geom_matrix(cut_off_end+cylinder_rad_cut:end,:)=min_coeff;
            
            
%             min_coeff=min(min(testi3(cut_off_start(1):cut_off_end,:)));
%             
%             testi3(1:cut_off_start-cylinder_rad_cut,:,:)=min_coeff;
%             
%             testi3(cut_off_end+cylinder_rad_cut:end,:,:)=min_coeff;
            
%             tr_geom_matrix(tr_geom_matrix <= 0) = min_coeff;
            
%             tr_geom_matrix(tr_geom_matrix > mean(abs(nonzeros(tr_geom_matrix(:))))) = mean(abs(nonzeros(tr_geom_matrix(:))));
            
%         end
        
%         keskiarvo = mean(mean(mean(tr_geom_matrix,3)));
%         tr_geom_matrix(tr_geom_matrix > keskiarvo * 2) = keskiarvo*2;
        tr_geom_matrix(isinf(tr_geom_matrix(:))) = 1;
        if mean(tr_geom_matrix(:)) > 1
            tr_geom_matrix(tr_geom_matrix > 1) = 1;
        end
        
%         testi = mean(SinDouble,2);
%         keskiarvo = mean(testi,1);
%         testi = bsxfun(@rdivide, keskiarvo, testi);
        
        if nargout >= 7
            varargout{7} = tr_geom_matrix;
        end
        
        tr_geom_matrix = repmat(tr_geom_matrix, 1, Ns, 1);
%         tr_geom_matrix=repmat(tr_geom_matrix,1,1,sino_amount);

        
        
        SinDouble=SinDouble.*tr_geom_matrix;
        
        norm_matrix=norm_matrix.*tr_geom_matrix;
        
    else
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% List-mode correction %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        x0 = max(detectors_x)/2;
        y0 = max(detectors_y)/2;
        [X,Y] = meshgrid(detectors_x, detectors_y);
        Y = Y';
        
        distance = ((abs((Y-detectors_y)*x0 - (X - detectors_x)*y0 + X.*detectors_y - Y.*detectors_x)./sqrt((Y-detectors_y).^2 + (X-detectors_x).^2)));
        
        start_ind = find(distance(:,1) <= options.FOVa_x/2/10*sqrt(2),1,'first');
        end_ind = find(distance(:,1) <= options.FOVa_x/2/10*sqrt(2),1,'last');

        tr_geom_matrix = zeros(options.detectors, options.detectors);

        for j = 1 : options.rings
            for i = j : options.rings
                if i == j
                    apu = true_coincidences(1 + (j-1)*options.det_per_ring:j*options.det_per_ring,1 + (i-1)*options.det_per_ring:i*options.det_per_ring);
                    KalPa = zeros(options.det_per_ring, options.cryst_per_block);
%                     hh = 1;
                    for k = start_ind : end_ind
                        temp = diag(apu,-(k-1));
                        for ll = 1 : options.cryst_per_block
                            KalPa(k,ll) = sum(temp(ll:options.cryst_per_block:end));
                        end
%                         hh = hh + 1;
                    end
                    KalPa(KalPa == 0) = NaN;
                    keskiarvo = nanmean(KalPa(start_ind:end_ind,:),1);
                    apu = keskiarvo ./ KalPa;
                    apu(isinf(apu)) = 1;
                    if nanmean(apu(:)) > 1
                        apu(apu > 1) = 1;
                    end
                    apu(isnan(apu)) = 1;
                    apu = repmat(apu,1, options.det_per_ring/options.cryst_per_block);
                    for kk = 1 : options.det_per_ring
                        apu(:,kk) = circshift(apu(:,kk),kk);
                    end
                    tr_geom_matrix(1 + (j-1)*options.det_per_ring:j*options.det_per_ring,1 + (i-1)*options.det_per_ring:i*options.det_per_ring) = apu;
                else
                    apu = true_coincidences(1 + (i-1)*options.det_per_ring:i*options.det_per_ring,1 + (j-1)*options.det_per_ring:j*options.det_per_ring);
                    KalPa = zeros(options.det_per_ring, options.cryst_per_block);
%                     hh = 1;
                    for k = start_ind : end_ind
                        temp1 = [diag(apu,-(k-1)) ; diag(apu,(k-1))];
                        for ll = 1 : options.cryst_per_block
                            KalPa(k,ll) = sum(temp1(ll:options.cryst_per_block:end));
                        end
%                         hh = hh + 1;
                    end
                    KalPa(KalPa == 0) = NaN;
                    keskiarvo = nanmean(KalPa(start_ind:end_ind,:),1);
                    apu = keskiarvo ./ KalPa;
                    apu(isinf(apu)) = 1;
                    if nanmean(apu(:)) > 1
                        apu(apu > 1) = 1;
                    end
                    apu(isnan(apu)) = 1;
                    apu = repmat(apu,1, options.det_per_ring/options.cryst_per_block);
                    for kk = 1 : options.det_per_ring
                        apu(:,kk) = circshift(apu(:,kk),kk);
                    end
                    tr_geom_matrix(1 + (i-1)*options.det_per_ring:i*options.det_per_ring,1 + (j-1)*options.det_per_ring:j*options.det_per_ring) = apu;
                end
            end
        end
        
        true_coincidences = true_coincidences .* tr_geom_matrix;
        norm_matrix = norm_matrix .* tr_geom_matrix;
        
        if nargout >= 7
            varargout{7} = tr_geom_matrix;
        end
        
    end
    
    time=toc;
    if options.verbose
        disp(['Transaxial geometric correction done (', num2str(time),'s)'])
    end
    
    
end


if options.verbose
    disp('Saving normalization data')
end

if options.use_raw_data
    norm_file = [folder options.machine_name '_normalization_listmode.mat'];
    norm_matrix = (double(norm_matrix(tril(true(size(norm_matrix)), 0))));
else
    norm_file = [folder options.machine_name '_normalization_' num2str(options.Ndist) 'x' num2str(options.Nang) '_span' num2str(options.span) '.mat'];
end
normalization = norm_matrix;
save(norm_file, 'normalization')


if options.verbose
    disp('Normalization matrix saved')
end


%% Return corrected data (convert to sparse if list-mode)

if nargout >= 1
    varargout{1} = norm_matrix;
end

if nargout >= 2
    if ~options.use_raw_data
        varargout{2} =SinDouble;
    else
        %         tic
        %form sparse matrix from
        
        %         corrected_data=zeros(1,numel(Fov_coincidences));
        %         for i=1:options.detectors
        %
        %             corrected_data(index_list_mode(i)+1:index_list_mode(i+1))=true_coincidences(i:end,i);
        %
        %             %%Norm matrix to vector form
        %
        %             %norm_coeffs(index_list_mode(i)+1:index_list_mode(i+1))=norm_matrix(i:end,i);
        %
        %         end
        %         corrected_data =
        true_coincidences = sparse(double(true_coincidences(tril(true(size(true_coincidences)), 0))));
        
        varargout{2} = true_coincidences;

%         varargout{2} = sparse(true_coincidences);
        
        %%Form sparse coeff_matrix%%
        % norm_coeffs=sparse(norm_coeffs);
        
        
        %Sparse corrected data
%         varargout{2} =sparse(corrected_data);
        
%         time=toc;
%         disp(['Data packed and returned (', num2str(time),'s)'])
        
    end
    
% end

end

disp('Normalization complete')