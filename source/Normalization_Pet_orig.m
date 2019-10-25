%% Function form

function [corrected_data, norm_matrix, detector_coeffs]=Normalization(options, attenuation_correction, scatter_correction, axial_geom_correction,...
    detector_effiency_correction, transaxial_geom_correction, return_data)


%%%%%%%%%%%%%%%%%%%%%%%% Vaihda tähän %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(options,'lor') && options.use_raw_data==1
    
    true_detectors=load('Inveon_AIVI_detector_locations_options.Ndistxoptions.Ndistx159_raw.mat', 'lor');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


uniform_source=0;

%if using uniform source for instance infinidesimal thickness cylinder use
%uniform_source=1;

%% Args

%options: Sinogram data, detector pair coordinates, ring heights, apparatus radius

%attenuation_correction: apply attenuation correction
%yes = [inner_cylinder_radius (cm) cylinder_attenuation_constant (cm^2/m)] | no = empty)
%If inner_cylinder_radius=inf, cylinder is assumed to cover entire FOV
%If left empty uniform illumination for each LOR is assumed

%scatter_correction: fit gaussian to scatter tail from cylinder
%normalization data (Too narrow scatter tail may lead to unaccurate fit).
%Not supported for list-mode data.

%axial_geom_correction: apply axial geometric correction (yes = 1 | no = 0)

%detector_effiency_correction: apply detector effiency correction. (Fansum = 1 | SPC = 2 | no = 0)
%Fan_sum uses 3-D fansum method for both data types or SPC "single-plane
%Casey" method for list mode-data (SPC computationally more expensive). SPC is
%supposed to be used with FOV covering source
%Fansum version includes block profile correction. Using
%detector_effiency_correction=2 with fansum uses block profile correction
%before detector effiency correction

%transaxial_geom_correction: apply transaxial geometric correction for plane source data (yes = 1 | no = 0)
%With cylinders transaxial correction produces appropriate coeffs for LORs
%passing cylinder (LORs passing near cylinder edges can also be inaccurate)
%Accuracy of radial sorting can be adjusted in transaxial correction section

%block_profile_correction: apply block profile correction. If using fansum
%correction is nested with effiency correction
%If using SPC block profile correction is done separately

%TRANAXIAL AND BLOOCK PROFILE CORRECTION GROUP SORTING CAN BE ADJUSTED IN
%THEIR SECTIONS


%return_data: return normalization corrected sinogram/list-mode data (yes = 1 | no = 0)

%% Returns
%Corrected sinogram
%Normalization matrix containing all coeffs


%% Coordinates

if ~options.use_raw_data
    
    z=sinogram_coordinates_3D(options);
    z_length = options.rings * options.cr_pz;
    z_true = linspace(0, z_length, options.rings + 1)./10;
    %mm --> cm
    z=z./10;
    
    [x,y]=sinogram_coordinates_2D(options);
    
    %mm --> cm
    x=x./10;
    y=y./10;
    
    x = permute(reshape(x, options.Ndist, options.Nang, 2), [2 1 3]);
    x = reshape(x, options.Ndist*options.Nang, 2);
    y = permute(reshape(y, options.Ndist, options.Nang, 2), [2 1 3]);
    y = reshape(y, options.Ndist*options.Nang, 2);
    
else
    
    %z
    z_length = options.rings * options.cr_pz;
    z = linspace(0, z_length, options.rings + 1);
    
    z=z./10;
    
end

%x and y detector coordinates
[detectors_x,detectors_y,~,~]=detector_coordinates(options);

%mm --> cm
detectors_x=detectors_x./10;
detectors_y=detectors_y./10;


detectors_ring=options.detectors/options.rings;
sino_amount=sum(options.segment_table);
segment_amount=length(options.segment_table);

%Ring radius

R=options.diameter/10/2; %cm

%Inner ring radius


if ~isempty(attenuation_correction)
    
    r=attenuation_correction(1); %cm
    
end


%% Scale stacked data (when using sinograms)


if ~options.use_raw_data
    
    tic
    
    %Convert coincidence data
    
    if isstruct(options.SinM)
        options.SinM=cell2mat(struct2cell(options.SinM));
    end
    SinDouble=permute(double(options.SinM), [2 1 3]);
    
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
    disp(['Spanned data scaled (', num2str(time),'s)'])
    
end


%% Extended sparse list-mode matrix (when using list-mode data)

if options.use_raw_data
    
    tic
    
    %Convert coincidences to matrix from
    
    if isa(options.coincidences,'cell')
        Fov_coincidences=cell2mat(options.coincidences);
    elseif isa(options.coincidences,'struct')
        cell_data=struct2cell(options.coincidences);
        Fov_coincidences=cell2mat(cell_data{1});
    end
    
    cell_data{1}=[];
    
    %Extend sparse matrix
    
    true_coincidences=zeros(options.rings*options.Nang*2);
    true_det_indices=true_coincidences;
    index_list_mode=zeros(1,options.detectors+1);
    index_list_mode(1)=0;
    
    %Extend sparse matrix
    
    Fov_coincidences=full(Fov_coincidences)';
    Fov_detectors=double(true_detectors.lor);
    
    %Reshape coincidences
    
    for i=1:options.detectors
        
        index_list_mode(i+1)=index_list_mode(i)+options.detectors-(i-1);
        
        true_coincidences(i:end,i)=Fov_coincidences(index_list_mode(i)+1:index_list_mode(i+1));
        
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
    
    
    if detector_effiency_correction==2
        
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
    figure(5135)
    subplot(1,2,1)
    imagesc(true_coincidences(51*detectors_ring+1:52*detectors_ring,50*detectors_ring+1:51*detectors_ring))
    
    time=toc;
    disp(['Sparse matrices extended and false detector pairs disposed (', num2str(time),'s)'])
    
end






%% Attenuation correction & scatter correction

if isempty(attenuation_correction)
    
    if options.use_raw_data==0
        true_radius=zeros(options.Nang,options.Ndist);
        
        for i=1:options.Ndist*options.Nang
            if ~isnan(true_radius(i))
                
                true_radius(i)=sqrt((x(i,1)-x(i,2))^2+(y(i,1)-y(i,2))^2);
                
            end
        end
        
        radius_in_cylinder=true_radius;
        
    else
        
        true_radius=ones(detectors_ring);
        
    end
end

if ~isempty(attenuation_correction)
    
    if attenuation_correction(1)==inf && options.use_raw_data==0
        true_radius=zeros(options.Nang,options.Ndist);
        
        for i=1:options.Ndist*options.Nang
            if ~isnan(true_radius(i))
                
                true_radius(i)=sqrt((x(i,1)-x(i,2))^2+(y(i,1)-y(i,2))^2);
                
            end
        end
        
        radius_in_cylinder=true_radius;
        
    elseif attenuation_correction(1)==inf && options.use_raw_data==1
        
        
        true_radius=ones(detectors_ring);
        
    end
    
    if exist('r','var') && ~options.use_raw_data
        
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
            true_radius=reshape(radius_in_cylinder,options.Nang,options.Ndist);
            
        end
        %min length in cylinder for each ring difference
        
        rad_min=min(radius_in_cylinder);
        
        
        
        % Out of cylinder indices
        
        cut_off_start=zeros(1,options.Nang);
        cut_off_end=cut_off_start;
        
        for j=1:options.Nang
            i=1;
            while i<=options.Ndist/2 && isnan(true_radius(j,i))
                
                i=i+1;
                
            end
            
            cut_off_start(j)=i;
            
        end
        
        
        for j=1:options.Nang
            i=options.Ndist/2;
            while i<=options.Ndist && ~isnan(true_radius(j,i))
                
                i=i+1;
                
            end
            cut_off_end(j)=i-1;
        end
        
    end
    
end

if (scatter_correction==1 || ~isempty(attenuation_correction)) && ~options.use_raw_data
    
    if length(attenuation_correction)==2
        attenuation_coeff=attenuation_correction(2);  %mass attenuation times density
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Scatter correction for sinogram%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if scatter_correction==1
        
        if r~=inf
            
            tic
            
            % Transaxial scatter counts
            
            Scatters=zeros(1,options.Ndist);
            for j=1:options.Ndist
                
                if sum(isnan(true_radius(:,j)))==options.Nang
                    
                    Scatters(j)=sum(sum(SinDouble(:,j,:)));
                    
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
                
                
                if sum(isnan(true_radius(:,j)))==options.Nang
                    
                    Sino_counts(u)=sum(sum(SinDouble(:,1:cut_off_start-cut_pix,u)))+sum(sum(SinDouble(:,cut_off_end+cut_pix:end,u)));
                    
                end
                
            end
            
            Scaling_factors=Sino_counts./mean(Sino_counts);
            
            
            %export scattered counts
            Scatter_matrix=zeros(size(SinDouble));
            for u=1:sino_amount
                
                Sino_gauss=Scaling_factors(u).*values'./options.Nang;
                
                Scatter_matrix(:,:,u)=repmat(Sino_gauss,options.Nang,1);
                SinDouble(:,:,u)=SinDouble(:,:,u)-Scatter_matrix(:,:,u);
                
            end
            
            SinDouble(SinDouble<0)=0;
            
            time=toc;
            disp(['Scatter correction done (', num2str(time),'s)'])
            
        else
            
            warning("Scatter correction is not supported for FOV covering source")
            
        end
        
    else
        
        warning("Normalization coefficients are calculated without scatter correction")
        
    end
    
    
    if ~isempty(attenuation_correction)
        
        tic
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Attenuation corrected sinogram form data %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        activity_coeffs=zeros(size(SinDouble));
        rad_coeff_matrix = zeros(size(SinDouble));
        
        h=zeros(1,sino_amount);
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
            
            
            
            h(u)=abs(z(ring,1)-z(ring,2))/10;
            if h(u)~=0
                rad_profile_axial=sqrt(radius_in_cylinder.^2+(isnan(radius_in_cylinder).*h(u)).^2);
            else
                rad_profile_axial=radius_in_cylinder;
            end
            
            rad_profile_axial(isnan(rad_profile_axial))=rad_min;
            rad_profile_axial(rad_profile_axial==0)=rad_min;
            
            if length(attenuation_correction)==2
                
                rad_coeffs=exp(-rad_min.*attenuation_coeff)./exp(-rad_profile_axial.*attenuation_coeff);
                
                %Blur attenuation coefficient sharp altering
                rad_coeffs=reshape(rad_coeffs,options.Nang,options.Ndist);
                I=[(ones(options.Nang,sz+1).*rad_coeffs(:,1:sz+1)) rad_coeffs];
                I=I(:);
                blurred=filter(Kernel,1,I');
                
                new_coeffs=reshape(blurred',options.Nang,options.Ndist+sz+1);
                new_coeffs=new_coeffs(:,sz+2:end);
                
                
                rad_coeff_matrix(:,:,start_ind:end_ind)=repmat(new_coeffs,1,1,options.segment_table(u));
                
                
                %Correct Sinograms of segment u
                
                
                
                SinDouble(:,:,start_ind:end_ind)=SinDouble(:,:,start_ind:end_ind).*rad_coeff_matrix(:,:,start_ind:end_ind);
            end
            
            %scale substance depth differences if applying transaxial correction (activity correction)
            
            radius_axial_mat=reshape(rad_profile_axial,options.Nang,options.Ndist);
            
            if exist('r','var') && uniform_source==0
                
                activity_coeffs(:,:,start_ind:end_ind)=max(max(radius_axial_mat))./repmat(radius_axial_mat,1,1,options.segment_table(u));
                activity_coeffs(isnan(true_radius))=1;
                
            end
        end
        
    end
    
    
    time=toc;
    disp(['Attenuation correction done (', num2str(time),'s)'])
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Attenuation correction for  list-mode data %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(attenuation_correction) && options.use_raw_data
    
    tic
    if length(attenuation_correction)==2
        attenuation_coeff=attenuation_correction(2);  %mass attenuation times density
    end
    %blur specs
    
    sigm = 2.0;
    %Window size
    sz = 4;
    x_gauss=-sz:sz;
    
    Exp_comp = -(x_gauss.^2)/(2*sigm^2);
    Kernel= exp(Exp_comp)/(sqrt(2*pi*sigm^2));
    
    if r~=inf
        intersects=zeros(detectors_ring,detectors_ring,2);
        a=zeros(detectors_ring);
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
    
    %Min length in cylinder for each ring difference
    
    rad_min=min(min(nonzeros(radius_in_cylinder)));
    
    true_radius=radius_in_cylinder;
    
    %
    
    if ~isempty(attenuation_correction) && (transaxial_geom_correction==1 || detector_effiency_correction==2)
        
        rad_profile_axial=zeros(detectors_ring,detectors_ring,options.ring_difference+1);
        h=zeros(1,options.ring_difference+1);
        activity_coeffs=zeros(size(true_coincidences));
        
        for u=0:options.ring_difference
            
            %Height from first ring of segment u
            
            h(u+1)=abs(z(1)-z(1+u));
            
            if h(u+1)~=0
                
                rad_profile_axial(:,:,u+1)=sqrt(radius_in_cylinder.^2+((radius_in_cylinder~=0).*(~isnan(radius_in_cylinder)).*h(u+1)).^2);
                
            else
                
                rad_profile_axial(:,:,u+1)=radius_in_cylinder;
                
            end
            
            rad_profile_axial(isnan(rad_profile_axial))=rad_min;
            rad_profile_axial(rad_profile_axial==0)=rad_min;
            if length(attenuation_correction)==2
                
                rad_coeffs(:,:,u+1)=exp(-rad_min.*attenuation_coeff)./exp(-rad_profile_axial(:,:,u+1).*attenuation_coeff);
                
                
                %blur coeffs with gaussian window
                temp_vec=rad_coeffs(:,:,u+1);
                [sorted,indices]=sort(temp_vec(:));
                I=[sorted(1).*ones(sz*2+1,1) ; sorted];
                
                blurred=filter(Kernel,1,I');
                
                new_coeffs=blurred(sz*2+2:end);
                temp_vec2=zeros(1,length(new_coeffs));
                
                for i=1:length(new_coeffs)
                    
                    temp_vec2(indices(i))=new_coeffs(i);
                    
                end
                
                new_coeffs_rad=reshape(temp_vec2,detectors_ring,detectors_ring);
                
                %Blur radial coeffs with
                
                for i=1:options.rings-u
                    
                    true_coincidences((u+i-1)*detectors_ring+1:(u+i)*detectors_ring,(i-1)*detectors_ring+1:i*detectors_ring)=new_coeffs_rad.*...
                        true_coincidences((u+i-1)*detectors_ring+1:(u+i)*detectors_ring,(i-1)*detectors_ring+1:i*detectors_ring);
                    
                end
                
            end
            if exist('r','var') && (transaxial_geom_correction==1 || detector_effiency_correction==2) && uniform_source==0
                
                for q=1:options.rings-u
                    
                    activity_coeffs((u+q-1)*detectors_ring+1:(u+q)*detectors_ring,(q-1)...
                        *detectors_ring+1:q*detectors_ring)=max(max(rad_profile_axial(:,:,u+1)))./rad_profile_axial(:,:,u+1);
                    
                end
                
            end
            
        end
        
        if length(attenuation_correction)==2
            time=toc;
            disp(['Attenuation correction done (', num2str(time),'s)'])
        end
        
    end
    
    
end

if options.use_raw_data==1
    
    if detector_effiency_correction==2
        
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


if isempty(attenuation_correction)
    
    warning("Normalization coefficients are calculated without attenuation correction")
    
end


%% Axial geom. factors for each ring

if axial_geom_correction==1
    tic
    
    if ~options.use_raw_data
        
        %%%Sinogram%%%
        
        Sincounts=zeros(1,length(SinDouble(1,1,:)));
        
        for u=1:length(SinDouble(1,1,:))
            
            Sincounts(u)=sum(sum(SinDouble(:,:,u)));
            
        end
        
        axial_geom_coeffs=zeros(1,sino_amount);
        
        %calculate coeffs (inverse ratio to mean counts)
        
        for u=1:sino_amount
            
            axial_geom_coeffs(u)=mean(Sincounts)/Sincounts(u);
            
            %correct sinogram
            SinDouble(:,:,u)=axial_geom_coeffs(u).*SinDouble(:,:,u);
            
            %add coeffs to normalization matrix
            norm_matrix(:,:,u)=axial_geom_coeffs(u).*norm_matrix(:,:,u);
            
        end
        
        time=toc;
        disp(['Axial correction done (', num2str(time),'s)'])
        
    else
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% List Mode %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        tic
        
        
        %total for all planes
        plane_counts=0;
        
        for u=1:options.ring_difference+1
            
            for j=1:options.rings-u+1
                if u==1
                    doubling=2;
                else
                    doubling=1;
                end
                
                plane_counts=plane_counts+sum(sum(true_coincidences(detectors_ring*(j+u-2)+...
                    1:detectors_ring*(j+u-1),detectors_ring*(j-1)+1:detectors_ring*j)))*doubling;
                
            end
            
        end
        
        plane_counts=plane_counts./(options.rings^2/2+options.rings/2);
        axial_coeffs=zeros(1,options.ring_difference);
        
        for u=1:options.ring_difference+1
            
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
        
        %Coeffs between cross projection planes
        
        for u=2:options.ring_difference+1
            
            
            for j=1:options.rings-u+1
                
                Lower_sum=0;
                Upper_sum=0;
                
                y_ring=detectors_ring*(j+u-2);
                x_ring=detectors_ring*(j-1);
                
                %Counts for first cross plane
                
                for p=start_lower:end_lower
                    
                    Lower_sum=Lower_sum+sum(true_coincidences(y_ring+start_ind_low_col(p):y_ring+end_ind_low_col(p),x_ring+p));
                    
                end
                
                %Counts for second cross plane
                
                for p=start_upper:end_upper
                    
                    Upper_sum=Upper_sum+sum(true_coincidences(y_ring+start_ind_up_col(p):y_ring+end_ind_up_col(p),x_ring+p));
                    
                end
                %Calculate coeffs
                
                cross_coeffs=mean([Upper_sum Lower_sum])./[Upper_sum Lower_sum];
                
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
        
        time=toc;
        disp(['Axial correction done (', num2str(time),'s)'])
        
    end
    
end




%% Crystal interference correction
if block_profile_correction==1 || transaxial_geom_correction==1 || detector_effiency_correction==1
    
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
    det_position=zeros(1,options.detectors/options.rings);
    
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
        
        avg_profile=zeros(options.cryst_per_block/2,options.cryst_per_block/2,options.Ndist);
        avg_profile_LOR_amount=avg_profile;
        not_found=0;
        detector_index_start=zeros(1,options.Nang*options.Ndist);
        detector_index_end=detector_index_start;
        index=zeros(options.Ndist*options.Nang,2);
        
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
                disp('error no corresponding coordinates found')
            end
            
            index(j,:)=[max(detector_index_start(j),detector_index_end(j))...
                min(detector_index_start(j),detector_index_end(j))];
            
            
            avg_profile_LOR_amount(index(j,1),index(j,2),ceil(j/options.Nang))=...
                avg_profile_LOR_amount(index(j,1),index(j,2),ceil(j/options.Nang))+1;
            
        end
        
        block_profile_coeffs=zeros(size(avg_profile_LOR_amount));
        
        
        
        %Calculate block profile coeffs by averaging over all sinograms%
        
        
        if block_profile_correction==1 && detector_effiency_correction==2
            
            block_profile_matrix=zeros(options.Nang,options.Ndist);
            
            avg_profile=zeros(size(avg_profile_LOR_amount));
            
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
            index1=reshape(index(:,1),options.Nang,options.Ndist);
            index2=reshape(index(:,2),options.Nang,options.Ndist);
            
            %coeff matrix for sinograms
            
            for p=1:options.Nang
                
                for i=1:options.Ndist
                    
                    if block_profile_coeffs(index1(p,i),index2(p,i),i)==inf
                        temp_coeffs=block_profile_coeffs(:,:,i);
                        block_profile_coeffs(index1(p,i),index2(p,i),i)=max(max(temp_coeffs(temp_coeffs~=inf)));
                    end
                    block_profile_matrix(p,i)=block_profile_coeffs(index1(p,i),index2(p,i),i);
                    
                end
                
                
            end
            
            block_profile_matrix=repmat(block_profile_matrix,1,1,sino_amount);
            SinDouble=SinDouble.*block_profile_matrix;
            norm_matrix=norm_matrix.*block_profile_matrix;
            
            
            time=toc;
            disp(['Block profile correction done (', num2str(time),'s)'])
            
        end
        
    else
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% List-mode form correction %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if detector_effiency_correction==2 || transaxial_geom_correction==1
            
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
            enviroment=(max_rad-min_rad)/(ceil(detectors_ring/sorting_factor)-1);
            
            %radial sorting
            radial_index=zeros(1,ceil(detectors_ring/sorting_factor)-1);
            
            for i=1:length(radial_index)
                
                radial_index(i)=min_rad+enviroment/2+(i-1)*enviroment;
                
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
                            while norm(radius_raw_data(i,j)-radial_index(p))>enviroment/2+10^-8 && over_index==0
                                
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
            
            
            if detector_effiency_correction==2
                
                profile_hits=zeros(max(det_position),max(det_position),max(max((det_index))));
                profile_avg=profile_hits;
                
                start=1;
                
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
                
            end
            
        end
        
        if detector_effiency_correction==2
            
            time=toc;
            disp(['Block profile correction done (', num2str(time),'s)'])
            
        end
        
    end
    
end

%% Transaxial geom. factors

if transaxial_geom_correction==1
    
    
    tic
    
    if ~options.use_raw_data
        
        
        %scale activities to same for each radial profile (if using cylinder data)
%         if exist('r','var') && uniform_source==0
%             
%             Temp_sino=SinDouble.*activity_coeffs;
%             
%         else
%             
            Temp_sino=SinDouble;
            
%         end
        % Radial profiles
        
        radial_index=zeros(1,options.Nang*options.Ndist);
        LOR_index=zeros(1,options.Nang*options.Ndist);
        radius=zeros(1,options.Nang*options.Ndist);
        
        
        %All radial lenghts
        for i=1:options.Ndist*options.Nang
            if ~isnan(true_radius(i))
                
                radius(i)=sqrt((x(i,1)-x(i,2))^2+(y(i,1)-y(i,2))^2);
                
            end
        end
        
        %Form radial groups
        
        %%%%%%%%%%% Group fining can be adjusted here if needed %%%%%%%%%%%%%%%%%%%
        
        %increasing factor increases groups
        fining_factor=2;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        maximum=max(radius);
        minimum=min(nonzeros(radius));
        enviroment=norm(maximum-minimum)/(fining_factor*options.Ndist-1);
        
        %Group radiusses
        
        for i=1:fining_factor*options.Ndist-1
            
            radial_index(i)=minimum+enviroment/2+(i-1)*enviroment;
            
        end
        
        radial_index=nonzeros(radial_index)';
        radial_times=zeros(1,length(radial_index));
        
        
        %Check if cylinder covers entire FOV
        
        if ~exist('r','var')
            
            true_radius=ones(1,options.Ndist*options.Nang);
            
        end
        
        %Sort LORs by radial profile
        
        for i=1:options.Ndist*options.Nang
            
            if ~isnan(true_radius(i))
                
                p=1;
                while norm(radius(i)-radial_index(p))>enviroment/2+10^-8 && radial_index(p)~=0
                    p=p+1;
                end
                LOR_index(i)=p;
                
                radial_times(p)=radial_times(p)+1;
                
            end
            
        end
        
        avg_counts_transaxial=zeros(max(LOR_index),1);
        tr_geom_matrix=zeros(options.Nang,options.Ndist);
        
        %Calculate average counts for each radial group
        
        for j=1:options.Ndist*options.Nang
            
            if LOR_index(j)~=0
                
                avg_counts_transaxial(LOR_index(j))=avg_counts_transaxial(LOR_index(j))+...
                    sum(Temp_sino(j-floor((j-1)/options.Nang)*options.Nang,ceil(j/options.Nang),:));
                
            end
            
        end
        
        avg_counts_transaxial=avg_counts_transaxial./(radial_times'*sino_amount);
        
        %Get rid of undetermined coeffs
        
        avg_counts_transaxial(isnan(avg_counts_transaxial))=0;
        avg_counts_transaxial(avg_counts_transaxial==inf)=0;
        non_zeros=nonzeros(avg_counts_transaxial);
        
        %transaxial coeffs
        
        tr_geom_coeffs=(sum(non_zeros)/length(non_zeros))./(avg_counts_transaxial);
        tr_geom_coeffs(tr_geom_coeffs==inf)=0;
        tr_geom_coeffs(isnan(tr_geom_coeffs))=0;
        
        %Form matrix from coeffs
        
        for i=1:options.Ndist*options.Nang
            
            if LOR_index(i)~=0
                
                tr_geom_matrix(i-floor((i-1)/options.Nang)*options.Nang,ceil(i/options.Nang))=tr_geom_coeffs(LOR_index(i));
                
            end
        end
        
        
        %Cut columns near edges of cylinder when cylinder not covering entire FOV
        
        if attenuation_correction(1)~=inf
            
            cylinder_rad_cut=ceil((cut_off_end(1)-cut_off_start(1))/50);
            
            min_coeff=min(min(tr_geom_matrix(:,cut_off_start(1):cut_off_end)));
            
            tr_geom_matrix(:,1:cut_off_start-cylinder_rad_cut)=min_coeff;
            
            tr_geom_matrix(:,cut_off_end+cylinder_rad_cut:end)=min_coeff;
            
        end
        
        tr_geom_matrix=repmat(tr_geom_matrix,1,1,sino_amount);
        
        SinDouble=SinDouble.*tr_geom_matrix;
        
        norm_matrix=norm_matrix.*tr_geom_matrix;
        
    else
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% List-mode correction %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Calculate radiusses of LORs passing cylinder
        
        radius_raw_data_true=zeros(size(true_detectors));
        
        for j=1:detectors_ring
            
            for i=[start_ind_up_col(j):end_ind_up_col(j) start_ind_low_col(j):end_ind_low_col(j)]
                
                if i~=0
                    
                    if true_radius(i,j)~=0 && ~isnan(true_radius(i,j))~=0
                        
                        radius_raw_data_true(i,j)=norm([detectors_x(j)-detectors_x(i) detectors_y(j)-detectors_y(i)]);
                        
                    end
                    
                end
                
            end
            
        end
        
        %Form radial groups
        
        max_rad=max(max(radius_raw_data_true));
        min_rad_true=min(min(nonzeros(radius_raw_data_true)));
        enviroment_true=(max_rad-min_rad_true)/(detectors_ring/5-1);
        
        %Radial sorting
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Fining factor can be adjusted to alter amount of radial groups
        
        fining_factor=5;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        radial_index_true=zeros(1,detectors_ring/fining_factor-1);
        
        for i=1:length(radial_index_true)
            
            radial_index_true(i)=min_rad_true+enviroment_true/2+(i-1)*enviroment_true;
            
        end
        
        radial_index_true=nonzeros(radial_index_true)';
        det_index_true=zeros(size(radius_raw_data_true));
        
        %Sort LORs by radial profile
        
        for j=1:detectors_ring
            
            for i=[start_ind_up_col(j):end_ind_up_col(j) start_ind_low_col(j):end_ind_low_col(j)]
                if i~=0
                    
                    if true_radius(i,j)~=0 && ~isnan(true_radius(i,j))~=0
                        
                        p=1;
                        over_index=0;
                        while norm(radius_raw_data_true(i,j)-radial_index_true(p))>enviroment_true/2+10^-8 && over_index==0
                            
                            p=p+1;
                            
                            if  p>length(radial_index_true)
                                over_index=1;
                                p=1;
                            end
                            
                        end
                        
                        if over_index==0
                            
                            det_index_true(i,j)=p;
                            
                        else
                            
                            det_index_true(i,j)=0;
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
        %Scale counts temporarily if using cylinder source
        
        if exist('r','var') && uniform_source==0
            
            temp_coincidences=activity_coeffs.*true_coincidences;
            
        else
            
            temp_coincidences=true_coincidences;
            
        end
        
        
        radial_sum=zeros(1,max(max(det_index_true)));
        radial_hits_true=radial_sum;
        
        for j=1:detectors_ring
            
            for i=[start_ind_low_col(j):end_ind_low_col(j) start_ind_up_col(j):end_ind_up_col(j)]
                
                if i~=0
                    
                    if det_index_true(i,j)~=0
                        if start_ind_up_col(j)<=i && i<=end_ind_up_col(j)
                            start=2;
                        else
                            start=1;
                        end
                        
                        for u=start:options.rings
                            
                            for p=1:options.rings-u+1
                                
                                radial_sum(det_index_true(i,j))=radial_sum(det_index_true(i,j))+...
                                    temp_coincidences((p+u-2)*detectors_ring+i,(p-1)*detectors_ring+j);
                                
                            end
                            
                        end
                        if start==2
                            radial_hits_true(det_index_true(i,j))=radial_hits_true(det_index_true(i,j))+options.rings^2/2-options.rings/2;
                        else
                            radial_hits_true(det_index_true(i,j))=radial_hits_true(det_index_true(i,j))+options.rings^2/2+options.rings/2;
                        end
                        
                    end
                end
                
            end
            
        end
        
        radial_sum=radial_sum./radial_hits_true;
        
        radial_sum(radial_sum==inf)=0;
        radial_sum(isnan(radial_sum))=0;
        
        radial_coeffs=mean(nonzeros(radial_sum))./radial_sum;
        
        radial_coeffs(radial_coeffs==inf)=0;
        
        
        for j=1:detectors_ring
            
            for i=[start_ind_low_col(j):end_ind_low_col(j) start_ind_up_col(j):end_ind_up_col(j)]
                
                if i~=0 && true_radius(i,j)~=0
                    
                    for u=1:options.rings
                        
                        for p=1:options.rings-u+1
                            
                            true_coincidences((p+u-2)*detectors_ring+i,(p-1)*detectors_ring+j)...
                                =true_coincidences((p+u-2)*detectors_ring+i,(p-1)*detectors_ring+j)*radial_coeffs(det_index(i,j));
                            
                            norm_matrix((p+u-2)*detectors_ring+i,(p-1)*detectors_ring+j)...
                                =norm_matrix((p+u-2)*detectors_ring+i,(p-1)*detectors_ring+j)*radial_coeffs(det_index(i,j));
                            
                        end
                        
                    end
                    
                end
            end
            
        end
        
        
    end
    
    time=toc;
    disp(['Transaxial geometric correction done (', num2str(time),'s)'])
    
    
end


%% Detector effiency factors

if detector_effiency_correction~=0
    
    tic
    
    if ~options.use_raw_data
        
        
        %  with 3-D fan-sum algorithm
        
        
        %effiency for LOR(i,j) = detector_effiency(i) x detector_effiency(j)
        
        z_rings=z(1:options.segment_table(1),1);
        det_num=zeros(options.Nang,options.Ndist,2);
        counts_det=zeros(detectors_ring,options.segment_table(1));
        coeff_matrix=zeros(size(SinDouble));
        
        %Determine each LORS detector numbers
        
        for u=1:options.Ndist*options.Nang
            
            found_it=0;
            t=1;
            while found_it==0
                
                if norm([x(u,1)-detectors_x(t) y(u,1)-detectors_y(t)])<0.001
                    
                    det_num(u-floor((u-1)/options.Nang)*options.Nang,floor((u-1)/options.Nang)+1,1)=t;
                    found_it=1;
                    
                end
                t=t+1;
            end
            found_it=0;
            t=1;
            while found_it==0
                
                if norm([x(u,2)-detectors_x(t) y(u,2)-detectors_y(t)])<0.001
                    
                    det_num(u-floor((u-1)/options.Nang)*options.Nang,floor((u-1)/options.Nang)+1,2)=t;
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
            
            
            
            for i=1:options.Nang
                
                for u=1:options.Ndist
                    
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
        
        
        if detector_effiency_correction~=2
            
            block_profile_counts=zeros(1,max(detector_index_start));
            block_profile_hits=block_profile_counts;
            detector_index_start=reshape(detector_index_start,options.Nang,options.Ndist);
            detector_index_end=reshape(detector_index_end,options.Nang,options.Ndist);
            det_number_block=zeros(1,detectors_ring);
        end
        for k=1:sino_amount
            
            for i=1:options.Nang
                
                for j=1:options.Ndist
                    
                    coeff_matrix(i,j,k)=coeffs_detectors(det_num(i,j,1),ring_1)*coeffs_detectors(det_num(i,j,2),ring_2);
                    
                    if detector_effiency_correction~=2
                        
                        block_profile_counts(detector_index_start(i,j))=block_profile_counts(detector_index_start(i,j))+coeffs_detectors(det_num(i,j,1),ring_1);
                        block_profile_counts(detector_index_end(i,j))=block_profile_counts(detector_index_end(i,j))+coeffs_detectors(det_num(i,j,2),ring_2);
                        block_profile_hits(detector_index_start(i,j))=block_profile_hits(detector_index_start(i,j))+1;
                        block_profile_hits(detector_index_end(i,j))=block_profile_hits(detector_index_end(i,j))+1;
                        
                    end
                end
                
            end
            
        end
        if detector_effiency_correction~=2
            block_profile_counts=block_profile_counts./block_profile_hits;
            block_profile_coeffs=mean(block_profile_counts)./block_profile_counts;
            
            p=1;
            
            for k=1:detectors_ring
                
                for i=1:options.Nang
                    
                    for j=1:options.Ndist
                        
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
        SinDouble=SinDouble.*coeff_matrix(1:options.Nang,:,:);
        
        %Coeffs to normalization matrix
        norm_matrix=norm_matrix.*coeff_matrix(1:options.Nang,:,:);
        
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
        
    else
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% List mode correction %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if isempty(attenuation_correction)
            
            true_radius=ones(detectors_ring);
            
        else
            
            if exist('r','var')
                
                if r==inf
                    true_radius=ones(detectors_ring);
                end
                
            end
            
        end
        
        if detector_effiency_correction==2
            
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
            
            for u=1:options.ring_difference+1
                
                for k=1:options.rings-u+1
                    
                    Sum_AB=zeros(detectors_ring);
                    
                    if exist('r','var') && uniform_source==0
                        
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
                    
                end
                
            end
            
            figure(5135)
            subplot(1,2,2)
            imagesc(true_coincidences(51*detectors_ring+1:52*detectors_ring,50*detectors_ring+1:51*detectors_ring))
            
            
        end
        
        if detector_effiency_correction==1
            
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
            
            
            %Separate block profile coeffs from effiency coeffs
            if block_profile_correction==1
                
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
                
                detector_coeffs=det_coeffs;
                
            end
            
        end
        
        
    end
    
    time=toc;
    disp(['Detector effiency corrected (', num2str(time),'s)'])
    
end
%% Return corrected data (convert to sparse if list-mode)

if return_data==1
    
    if ~options.use_raw_data
        corrected_data=SinDouble;
    else
        tic
        %form sparse matrix from
        
        corrected_data=zeros(1,numel(Fov_coincidences));
        for i=1:options.detectors
            
            corrected_data(index_list_mode(i)+1:index_list_mode(i+1))=true_coincidences(i:end,i);
            
            %%Norm matrix to vector form
            
            %norm_coeffs(index_list_mode(i)+1:index_list_mode(i+1))=norm_matrix(i:end,i);
            
        end
        
        %%Form sparse coeff_matrix%%
        % norm_coeffs=sparse(norm_coeffs);
        
        
        %Sparse corrected data
        corrected_data=sparse(corrected_data);
        
        time=toc;
        disp(['Data packed and returned (', num2str(time),'s)'])
        
    end
    
else
    corrected_data=[];
end

end