function [New_randoms] = Randoms_variance_reduction(Randoms,options)
%% RANDOMS VARIANCE REDUCTION
% This function applies the randoms variance reduction on the input randoms
% or scatter measurements. The variance reduction is computed by using the
% 3D fan sum algorithm.
% Example:
%   new_randoms = Randoms_variance_reduction(Randoms, options)
% INPUTS:
%   Randoms = The randoms measurement sinogram/raw list-mode data. The
%   input data format needs to be the same as used elsewhere in OMEGA.
%   options = Contains the machine and (optionally) sinogram specific
%   information (e.g. detectors per ring, image (matrix) size, FOV size).
% OUTPUT:
%   New_randoms = Randoms with variance reduction applied.
%
% See also normalization_coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019  Anssi Manninen
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

%Sinogram data


%coordinates
if ~options.use_raw_data
    
    z=sinogram_coordinates_3D(options);
    
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
% R=options.diameter/10/2;

if ~options.use_raw_data
    
    %  with 3-D fan-sum algorithm
    
    if size(Randoms,3) == 1
        Randoms = reshape(Randoms, options.Ndist, options.Nang, options.TotSinos);
    end
    
    
    %effiency for LOR(i,j) = detector_effiency(i) x detector_effiency(j)
    
    det_num=zeros(options.Ndist,options.Nang,2);
    randoms_det=zeros(detectors_ring,options.segment_table(1));
    coeff_matrix=zeros(size(Randoms));
    sino_amount=sum(options.segment_table);
    %Determine each LORS detector numbers
    z = z - min(z(:));
    ring = round(z ./ z(2,1)) + 1;
    
    for u=1:options.Ndist*options.Nang
        
        found_it=0;
        t=1;
        while found_it==0
            
            if norm([x(u,1)-detectors_x(t) y(u,1)-detectors_y(t)])<0.001
                
                det_num(u-floor((u-1)/options.Ndist)*options.Ndist,floor((u-1)/options.Ndist)+1,1)=t;
                found_it=1;
                
            end
            t=t+1;
        end
        found_it=0;
        t=1;
        while found_it==0
            
            if norm([x(u,2)-detectors_x(t) y(u,2)-detectors_y(t)])<0.001
                
                det_num(u-floor((u-1)/options.Ndist)*options.Ndist,floor((u-1)/options.Ndist)+1,2)=t;
                found_it=1;
                
            end
            t=t+1;
        end
        
    end
    
    hits_det=zeros(detectors_ring,options.segment_table(1));
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
    
    randoms_det(1:maksimi1) = accumarray(testi1, Randoms(:));
    randoms_det(1:maksimi2) = randoms_det(1:maksimi2)' + accumarray(testi2,  Randoms(:));
    hits_det(1:maksimi1) = accumarray(testi1(:),1);
    hits_det(1:maksimi2) = hits_det(1:maksimi2)' + accumarray(testi2(:),1);
    
    %Calculate coeffs (mean inverse)
    
    randoms_det = randoms_det./hits_det;
    coeffs_detectors = bsxfun(@rdivide, mean(randoms_det,1), randoms_det);
    
    for k=1:sino_amount
        for i=1:options.Ndist
            for j=1:options.Nang
                
                coeff_matrix(i,j,k) = coeffs_detectors(det_num(i,j,1),ring(k,1))*coeffs_detectors(det_num(i,j,2),ring(k,2));
                
            end
        end
    end
    
    New_randoms=Randoms.*coeff_matrix;
    
    %Scale noise reduced sinograms total counts to original amount
    New_randoms = New_randoms .* (sum(sum(Randoms)) ./ sum(sum(New_randoms)));
    
else
    %% List-mode data
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Apply 3-D fansum to reduce randoms variance %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~isa(Randoms,'sparse')
        
        tic
        
        %Convert delayed_coincidences to matrix from
        
        if isa(Randoms,'cell')
            Fov_coincidences=cell2mat(Randoms);
        else
            Fov_coincidences = Randoms;
%             cell_data=struct2cell(Randoms);
%             Fov_coincidences=cell2mat(cell_data{1});
        end
        
        %Extend sparse matrix
        
        Randoms=zeros(options.rings*options.det_per_ring);
        index_list_mode=zeros(1,options.detectors+1);
        index_list_mode(1)=0;
        
        %Extend sparse matrix
        
        Fov_coincidences=full(Fov_coincidences)';
        
        %Reshape delayed_coincidences
        
        for i=1:options.detectors
            
            index_list_mode(i+1)=index_list_mode(i)+options.detectors-(i-1);
            
            Randoms(i:end,i)=Fov_coincidences(index_list_mode(i)+1:index_list_mode(i+1));
            
        end
        
    end
    
    
%     if r==inf
        
        true_radius=ones(detectors_ring);
        
%     else
%         
%         
%         intersects=zeros(detectors_ring,detectors_ring,2);
%         a=zeros(detectors_ring);
%         c=a;
%         radius_in_cylinder=c;
%         
%         for j=1:detectors_ring
%             
%             
%             for i=1:detectors_ring
%                 
%                 if i~=j
%                     
%                     %each LOR detector coordinates
%                     x_start=detectors_x(j);
%                     x_end=detectors_x(i);
%                     y_start=detectors_y(j);
%                     y_end=detectors_y(i);
%                     
%                     if (x_end-x_start)~=0
%                         a(i,j)=(y_end-y_start)/(x_end-x_start);
%                         if x_start<x_end
%                             c(i,j)=y_start-a(i,j)*x_start;
%                         else
%                             c(i,j)=y_end-a(i,j)*x_end;
%                         end
%                     else
%                         a(i,j)=NaN;
%                         c(i,j)=NaN;
%                     end
%                     
%                     %intersection points with inner cylinder
%                     if (x_end-x_start)~=0
%                         intersects(i,j,:)=roots([(a(i,j)^2+1) (2*a(i,j)*c(i,j)-2*R-2*R*a(i,j)) (c(i,j)^2+2*R^2-2*R*c(i,j)-r^2)]);
%                         
%                         if imag(intersects(i,j,1))<0.0001 && imag(intersects(i,j,2))<0.0001
%                             radius_in_cylinder(i,j)=sqrt((intersects(i,j,1)-intersects(i,j,2))^2+((intersects(i,j,1)-intersects(i,j,2))*a(i,j))^2);
%                         else
%                             radius_in_cylinder(i,j)=NaN;
%                         end
%                         
%                     else
%                         intersects(i,j,:)=roots([1 -2*R (x_start-R)^2-r^2+R^2]);
%                         if imag(intersects(i,j,1))<0.0001 && imag(intersects(i,j,2))<0.0001
%                             radius_in_cylinder(i,j)=(intersects(i,j,1)-intersects(i,j,2));
%                         else
%                             radius_in_cylinder(i,j)=NaN;
%                         end
%                         
%                     end
%                     
%                 end
%                 
%             end
%             
%             
%         end
%         
%         true_radius=radius_in_cylinder;
%         true_radius(isnan(true_radius))=0;
%         true_radius(true_radius~=0)=1;
%         
%         
%     end
    
    New_randoms=zeros(size(Randoms));
    cylinder_counts=zeros(size(Randoms));
    
    for u=1:options.rings
        
        for k=1:options.rings-u+1
            
            cylinder_counts(detectors_ring*(u+k-2)+1:detectors_ring*(u+k-1),detectors_ring*(k-1)+1:detectors_ring*k)=Randoms(detectors_ring*(u+k-2)+1:detectors_ring*(u+k-1),detectors_ring*(k-1)+1:detectors_ring*k).*true_radius;
            
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
            
            New_randoms(start_ind_col:end,current_detector)=Randoms...
                (start_ind_col:end,current_detector)*det_coeffs(current_detector);
            
            
            New_randoms(current_detector,1:end_ind_row)=...
                Randoms(current_detector,1:end_ind_row)*det_coeffs(current_detector);
            
            
        end
        
    end
    
    %Scale new randoms
    
    for u=1:options.rings
        
        for k=1:options.rings-u+1
            
            start=detectors_ring*(u+k-2)+1;
            end_col=detectors_ring*(u+k-1);
            start_row=detectors_ring*(k-1)+1;
            end_row=detectors_ring*k;
            
            
            New_randoms(start:end_col,start_row:end_row)=...
                New_randoms(start:end_col,start_row:end_row)*...
                sum(sum(cylinder_counts(start:end_col,start_row:end_row)))/sum(sum(New_randoms(start:end_col,start_row:end_row)));
            
        end
        
    end
    
    New_randoms = New_randoms(tril(true(options.detectors,options.detectors),0));
    
    New_randoms = sparse(New_randoms);
    
end
if options.verbose
    disp('Randoms variance reduction completed')
end
