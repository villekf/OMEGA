function [ lor, indices, alkiot, discard] = siddon_atten_for_KYS( pikselikoko, Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, iij, jji, kkj, yy, xx , NSinos, NSlices, vaimennus)
%
%SIDDON Function to calculate the Siddon's algorithm
%   This function calculates the line intersections in a pixel space by
%   using the Siddon's algorithm and converts them to probabilities. This
%   function also determines the LORs that are not in the desired slices.
%   Only those LORs that are inside the slice interval are taken into
%   account. Returns the LORs, the pixel coordinates, and the
%   probabilities, in addition to the detector pairs that are not included
%   (i.e. they are discarded). This code doesn't take into account the LORs that
%   do not originate from the specified block-interval.
%
% INPUTS: 
%
% LL contains the detector numbers of each LOR in the measurement vector
%
% LOOP_VAR_PAR is the length of the subset used for the current
% subiteration
%
% PIKSELIKOKO defines the size of the image matrix in x and y directions
% 
% NY is the length of the y direction
% 
% NX same for x
%
% NZ same for z
% 
% D defines the square (pixel) length, that is the distance between pixel
% boundaries in x and y directions
% 
% DZ same for z direction
%
% BY the distance from the origin to the pixel space in y direction
%
% BX same for x
% 
% BZ same for z
%
% OSAJOUKOT the subset intervals used for OSEM reconstruction
%
% TEST current subiteration (subset)
%
% Z_DET detector ring locations in z dimension (axial)
%
% X detector ring locations in x direction (radial)
% 
% Y same for y
% 
% IIJ pixel boundaries in x directon (iij=0:Nx)
% 
% JJI same for y
%
% KKJ same for z
% 
% YY 
z_det = z_det(1:NSinos,:);
loop_var_par = NSinos*size(x,1);
indices=cell((loop_var_par),1);
alkiot=cell((loop_var_par),1);
% indices=zeros(100,(loop_var_par)/100,'int32');
% alkiot=zeros(100,(loop_var_par)/100,'int32');
lor=zeros((loop_var_par),2,'int32');
% lor=zeros((loop_var_par),2,'single');
discard=true((loop_var_par),1);
% apu_var=osajoukot(test);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%vaimennusasiat:
vaimennuskuva = vaimennus(:,:,1:NSlices);
vaimennuskuva=vaimennuskuva(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 436142
%1587419
ll=int32(1:size(x,1));
ll = repmat(ll,1,NSinos);
lz = int16(1:NSinos);
lz = repelem(lz,size(x,1));
zmax = max(max(z_det));
if zmax==0
    zmax = 1;
end
Nyy = (Ny-1);
Nxx = (Nx-1);
for lo=int32(1:loop_var_par)
%     if lo > 1 && ~isempty((indices{lo-1})) && min(indices{lo-1})<1
%         lo
%         return
%     elseif lo > 1 && ~isempty((indices{lo-1})) && max(indices{lo-1})>N
%         lo
%         return
%     end
%       lo
    tx_n=0;
    ty_n=0;
    tz_n=0;
    imax=int32(0);
    imin=int32(0);
    jmax=int32(0);
    jmin=int32(0);

    % Jos detektorit vastakkaiset niin laskeminen ei onnistu Siddonin
    % algoritmilla
    
%     if random_order==1
%         detektorit=LL(lo+apu_var,:);
%     else
%          detektorit=int32(LL(lo,:));
%     end

%     if detektorit(1)==detektorit(2)
%         discard(lo)=false;
%         continue
%     end
%     
%     %break
%     loop=idivide(detektorit,320,'ceil');
%     
%     if loop(1)==loop(2)
%         if  loop(1)>blocks || loop(1)<block1 || loop(2)>blocks || loop(2)<block1
%             %disp('skip')
%             discard(lo)=false;
%             continue
%         end
% %         luuppi=1;
%         zs=loop(1);
%         zd=zs;
% %         break
%     else
%         if loop(1)>blocks || loop(1)<block1 || loop(2)>blocks || loop(2)<block1
%             %disp('skip')
%             discard(lo)=false;
%             continue
%         end
% %         luuppi=2;
%         zs=z_det(loop(1));
%         zd=z_det(loop(2));
%     end
    
    xs=x(ll(lo),1);
    xd=x(ll(lo),2);
    ys=y(ll(lo),1);
    yd=y(ll(lo),2);
    zs = z_det(lz(lo),1);
    zd = z_det(lz(lo),2);
    z_loop=int32((zs/zmax)*(NSlices-1)+1);
%     if lo == 1008289
%         ((zs/zmax)*(NSlices-1)+1)
%         zmax
%         z_loop
%         break
%     end
    
    x_diff = xd - xs;
    y_diff = yd - ys;
    z_diff = zd - zs;
    
    if zs==zd
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
%         discard(lo)=false;
%         continue
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if abs(y_diff)<1e-8
            
            if yd<= max(yy) && yd>=min(yy)
                tempi=(1:pikselikoko)';
                [~,apu]=min(abs(yy(1:pikselikoko)-yd));
                tempj=ones(size(tempi),'int32')*(int32(apu));
                templ_ijk=ones(size(tempi),'single')*d;
                tempk=(z_loop)*ones(size(tempi),'int32');
                temp=sum(templ_ijk);
                [tempk, tempi] = sort(((tempj-1).*Nyy+tempi+(Nxx.*Nyy.*(tempk-1))));
                indices{lo}=tempk;                
                alkiot{lo}=single(templ_ijk(tempi).*(1/temp))*exp(-vaimennuskuva(indices{lo})'*templ_ijk(tempi));
                lor(lo,:)=[lo,int32(int32(length(alkiot{lo})))];
                continue
            end
            % Jos aiemmat ehdot ei t�yty (LOR ei mene
            % pikseliavaruuden l�pi), siirry seuraavaan
            % detektoriin
            discard(lo)=false;
            continue
        elseif abs(x_diff)<1e-8
            if xd<= max(xx) && xd>=min(xx)
                tempj=(1:pikselikoko)';
                [~,apu]=min(abs(xx(1:pikselikoko)-xd));
                tempi=ones(size(tempj),'int32')*(int32(apu));
                templ_ijk=ones(size(tempj),'single')*d;
                tempk=(z_loop)*ones(size(tempj),'int32');
                temp=sum(templ_ijk);
                [tempk, tempi] = sort(((tempj-1).*Nyy+tempi+(Nxx.*Nyy.*(tempk-1))));
                indices{lo}=tempk;
                alkiot{lo}=single(templ_ijk(tempi).*(1/temp))*exp(-vaimennuskuva(indices{lo})'*templ_ijk(tempi));
                lor(lo,:)=[lo,int32(length(alkiot{lo}))];
                continue
            end
            discard(lo)=false;
            continue
        end
        % Parametriesitykseen tarvittavat t-arvot
        
        tx=((bx+iij.*d)-xs)./x_diff;
        
        % Sama y:lle
        
        ty=((by+jji.*d)-ys)./y_diff;
        txmin=min(tx(1),tx(end));
        txmax=max(tx(1),tx(end));
        tymin=min(ty(1),ty(end));
        tymax=max(ty(1),ty(end));
        
        % Yhteiset arvot
        tmin=max(txmin,tymin);
        tmax=min(txmax,tymax);
        
%         tmin
        
%         lor(lo,:)=[single(lo),tmin];
        
        if tmin >= tmax
           %disp('skip')
            discard(lo)=false;
            continue
        end

        
        % PET:ss� l�hdekin on detektori
        % i kuvastaa x-akselin pikseliviivova, j y-akselin
        if xs<xd % Ensiksi jos l�hde on alempana
            if tmin==txmin % Jos yhteinen minimi t on sama kuin x-akselin parametrisoidun t niin tulopiste on i=1
                imin=int32(1);
            else % Muulloin lasketaan pisteen arvo pienimm�ss� t:n arvossa (esim. leikkaus y-akselilla), v�hennet��n
                % et�isyys origosta kuvaan ja jaetaan pikselien v�lisell�
                % et�isyydell�
                pxt=xs+tmin*x_diff;
                imin=int32(ceil((pxt-bx)/d));
            end
            if tmax==txmax % Sama poistumistielle
                imax=int32(Nx - 1);
            else
                pxt=xs+tmax*x_diff;
                imax=int32(floor((pxt-bx)/d));
            end
            apu=(imin:1:imax)+1;
            tx_n=tx((apu));
        else % Jos l�hde on ylemp�n�
            if tmin==txmin
                imax=int32(Nx-2);
            else
                pxt=xs+tmin*x_diff;
                imax=int32(floor((pxt-bx)/d));
            end
            if tmax==txmax
                imin=int32(0);
            else
                pxt=xs+tmax*x_diff;
                imin=int32(ceil((pxt-bx)/d));
            end
            apu=(imax:-1:imin)+1;
            tx_n=tx((apu));
        end
        
        % Sama y-akselille
        if ys<yd
            if tmin==tymin
                jmin=int32(1);
            else
                pyt=ys+tmin*y_diff;
                jmin=int32(ceil((pyt-by)/d));
            end
            if tmax==tymax
                jmax=int32(Ny - 1);
            else
                pyt=ys+tmax*y_diff;
                jmax=int32(floor((pyt-by)/d));
            end
            apu=(jmin:1:jmax)+1;
            ty_n=ty((apu));
        else
            if tmin==tymin
                jmax=int32(Ny-2);
            else
                pyt=ys+tmin*y_diff;
                jmax=int32(floor((pyt-by)/d));
            end
            
            if tmax==tymax
                jmin=int32(0);
            else
                pyt=ys+tmax*y_diff;
                jmin=int32(ceil((pyt-by)/d));
            end
            apu=(jmax:-1:jmin)+1;
            ty_n=ty((apu));
        end
        
%         if lo == 67489
%             jmax
%             jmin
%             imax
%             imin
%             break
%         end

%         if lo == 1075
%             jmax
%             jmin
%             break
%         end
        
%         if lo == 1075
%             ys
%             yd
%             jmax
%             jmin
%             break
%         end
        
        % Otetaan vain ne t:n arvot, jotka ovat s�teen kulkureitill�
        

        % Koko s�teen pituus
        L=sqrt(x_diff^2+y_diff^2);
        
%         if lo == 14575
%             ty_n
%             break
%         end
        
        % Pistet��n parametri t:n arvot pienemm�st� suurempaan
%         tt=sort([tmin;tx_n';ty_n']);
        % Poistetaan tuplat
        tt=unique([tmin;tx_n';ty_n']);
        
%         if lo == 14575
%             tt
%             length(tt)
%             break
%         end
        
        % Pisteet joiden l�pi s�de kulkee (riste�� akselien kanssa)
        
        pxt=xs+((tt(2:end)+tt(1:end-1))./2).*x_diff;
        tempi=int32(floor((pxt-bx)./d)+1);
        pyt=ys+((tt(2:end)+tt(1:end-1))./2).*y_diff;
        tempj=int32(floor((pyt-by)./d));
%         tempj(tempj>63)=63;
        
        tempk=(z_loop - 1)*ones(size(tempj),'int32');
        
%         if lo == 13103601
%             tempk(1)
%             z_loop
%             break
%         end
        
        % Lasketaan kuljettu matka jokaisessa pikseliss�
        
        templ_ijk=(diff(tt)).*L;
        
        temp=sum(templ_ijk);
        [tempk, tempi] = sort(((tempj).*Nyy+tempi+(Nxx.*Nyy.*(tempk))));
%         if lo == 33605
%             tempk
%         end
%         if tempk(1) == 65345 && length(tempk) > 1 && tempk(2) == 65409
%             lo
%             break
%         end
        indices{lo}=tempk;
%         if max(tempk) == 4031
%             lo
%         end
        alkiot{lo}=single(templ_ijk(tempi)).*((1/temp)*exp(-vaimennuskuva(indices{lo})'*templ_ijk(tempi)));
        lor(lo,:)=[lo,int32(length(alkiot{lo}))];
%         break
%         lo
        
%         if lo > 291
%             break
%         end
%         if sum(tempk>4096)>0
%             lo
%             break
%         end
        
    else

        
        
        imax=int32(0);
        imin=int32(0);
        jmax=int32(0);
        jmin=int32(0);
        kmax=int32(0);
        kmin=int32(0);

        
        
        
        % Jos detektorit vastakkaiset niin laskeminen ei onnistu Siddonin
        % algoritmilla
        
        if abs(y_diff)<1e-8
            if yd<= max(yy) && yd>=min(yy)

                tx=((bx+iij.*d)-xs)./x_diff;
                tz=((bz+kkj.*dz)-zs)./z_diff;
                
                
                txmin=min(tx(1),tx(end));
                txmax=max(tx(1),tx(end));
                tzmin=min(tz(1),tz(end));
                tzmax=max(tz(1),tz(end));
                
                % Yhteiset arvot
                tmin=max(txmin,tzmin);
                tmax=min(txmax,tzmax);
                
                if tmin >= tmax
                   %disp('skip')
                    discard(lo)=false;
                    continue
                end

                
                if xs<xd % Ensiksi jos l�hde on alempana
                    if tmin==txmin % Jos yhteinen minimi t on sama kuin x-akselin parametrisoidun t niin tulopiste on i=1
                        imin=int32(1);
                    else
                        pxt=xs+tmin*x_diff;
                        imin=int32(ceil((pxt-bx)/d));
                    end
                    if tmax==txmax % Sama poistumistielle
                        imax=int32(Nx - 1);
                    else
                        pxt=xs+tmax*x_diff;
                        imax=int32(floor((pxt-bx)/d));
                    end
                elseif xs>xd % Jos l�hde on ylemp�n�
                    if tmin==txmin
                        imax=int32(Nx-2);
                    else
                        pxt=xs+tmin*x_diff;
                        imax=int32(floor((pxt-bx)/d));
                    end
                    if tmax==txmax
                        imin=int32(0);
                    else
                        pxt=xs+tmax*x_diff;
                        imin=int32(ceil((pxt-bx)/d));
                    end
                end
                
                if zs<zd
                    if tmin==tzmin
                        kmin=int32(1);
                    else
                        pzt=zs+tmin*z_diff;
                        kmin=int32(ceil((pzt-bz)/dz));
                    end
                    if tmax==tzmax
                        kmax=int32(Nz - 1);
                    else
                        pzt=zs+tmax*z_diff;
                        kmax=int32(floor((pzt-bz)/dz));
                    end
                elseif zs>zd
                    if tmin==tzmin
                        kmax=int32(Nz-2);
                    else
                        pzt=zs+tmin*z_diff;
                        kmax=int32(floor((pzt-bz)/dz));
                    end
                    if tmax==tzmax
                        kmin=int32(0);
                    else
                        pzt=zs+tmax*z_diff;
                        kmin=int32(ceil((pzt-bz)/dz));
                    end
                end
                
                if xs<xd
                    apu=(imin:1:imax)+1;
                    tx_n=tx((apu));
                elseif xs>xd
                    apu=(imax:-1:imin)+1;
                    tx_n=tx((apu));
                end
                
                
                if zs<zd
                    apu=(kmin:1:kmax)+1;
                    tz_n=tz((apu));
                elseif zs>zd
                    apu=(kmax:-1:kmin)+1;
                    tz_n=tz((apu));
                end
                
                L=sqrt(x_diff^2+z_diff^2);
                
%                 tt=sort([tmin;tx_n';tz_n']);
                tt=unique([tmin;tx_n';tz_n']);
                
                pxt=xs+((tt(2:end)+tt(1:end-1))./2).*x_diff;
                tempi=int32(floor((pxt-bx)./d) + 1);
                pzt=zs+((tt(2:end)+tt(1:end-1))./2).*z_diff;
                tempk=int32(floor((pzt-bz)./dz));
                
                [~,apu]=min(abs(yy(1:pikselikoko)-yd));
%                 break
                tempj=ones(size(tempi),'int32')*((apu));
                
                
                templ_ijk=(diff(tt)).*L;
                temp=sum(templ_ijk);
                [tempk, tempi] = sort(((tempj).*Nyy+tempi+(Nxx.*Nyy.*(tempk))));
                indices{lo}=tempk;
                alkiot{lo}=single(templ_ijk(tempi).*(1/temp))*exp(-vaimennuskuva(indices{lo})'*templ_ijk(tempi));
                lor(lo,:)=[lo,int32(length(alkiot{lo}))];
                
                continue
            
            end
            discard(lo)=false;
            continue
        elseif abs(x_diff)<1e-8
            if xd<= max(xx) && xd>=min(xx)
                
                ty=((by+iij.*d)-ys)./y_diff;
                
                tz=((bz+kkj.*dz)-zs)./z_diff;
                
                tymin=min(ty(1),ty(end));
                tymax=max(ty(1),ty(end));
                tzmin=min(tz(1),tz(end));
                tzmax=max(tz(1),tz(end));
                
                % Yhteiset arvot
                tmin=max(tymin,tzmin);
                tmax=min(tymax,tzmax);
                
                if tmin >= tmax
                   %disp('skip')
                    discard(lo)=false;
                    continue
                end
                
                if ys<yd % Ensiksi jos l�hde on alempana
                    if tmin==tymin % Jos yhteinen minimi t on sama kuin y-akselin parametrisoidun t niin tulopiste on i=1
                        jmin=int32(1);
                    else
                        pyt=ys+tmin*y_diff;
                        jmin=int32(ceil((pyt-by)/d));
                    end
                    if tmax==tymax % Sama poistumistielle
                        jmax=int32(Ny - 1);
                    else
                        pyt=ys+tmax*y_diff;
                        jmax=int32(floor((pyt-by)/d));
                    end
                elseif ys>yd % Jos l�hde on ylemp�n�
                    if tmin==tymin
                        jmax=int32(Ny-2);
                    else
                        pyt=ys+tmin*y_diff;
                        jmax=int32(floor((pyt-by)/d));
                    end
                    if tmax==tymax
                        jmin=int32(0);
                    else
                        pyt=ys+tmax*y_diff;
                        jmin=int32(ceil((pyt-by)/d));
                    end
                end
                
                if zs<zd
                    if tmin==tzmin
                        kmin=int32(1);
                    else
                        pzt=zs+tmin*z_diff;
                        kmin=int32(ceil((pzt-bz)/dz));
                    end
                    if tmax==tzmax
                        kmax=int32(Nz - 1);
                    else
                        pzt=zs+tmax*z_diff;
                        kmax=int32(floor((pzt-bz)/dz));
                    end
                elseif zs>zd
                    if tmin==tzmin
                        kmax=int32(Nz-2);
                    else
                        pzt=zs+tmin*z_diff;
                        kmax=int32(floor((pzt-bz)/dz));
                    end
                    if tmax==tzmax
                        kmin=int32(0);
                    else
                        pzt=zs+tmax*z_diff;
                        kmin=int32(ceil((pzt-bz)/dz));
                    end
                end
                
                if ys<yd
                    apu=(jmin:1:jmax)+1;
                    ty_n=ty((apu));
                elseif ys>yd
                    apu=(jmax:-1:jmin)+1;
                    ty_n=ty((apu));
                end
                
                if zs<zd
                    apu=(kmin:1:kmax)+1;
                    tz_n=tz((apu));
                elseif zs>zd
                    apu=(kmax:-1:kmin)+1;
                    tz_n=tz((apu));
                end
                
                L=sqrt(y_diff^2+z_diff^2);
                
                %tt=sort([tmin;ty_n';tz_n']);
                tt=unique([tmin;ty_n';tz_n']);
                
                pyt=ys+((tt(2:end)+tt(1:end-1))./2).*y_diff;
                tempj=int32(floor((pyt-by)./d));
                pzt=zs+((tt(2:end)+tt(1:end-1))./2).*z_diff;
                tempk=int32(floor((pzt-bz)./dz));
                
                [~,apu]=min(abs(xx(1:pikselikoko)-xd));
                tempi=ones(size(tempj),'int32')*(int32(apu));
                
                
                templ_ijk=(diff(tt)).*L;
                
                temp=sum(templ_ijk);
                [tempk, tempi] = sort(((tempj).*Nyy+tempi+(Nxx.*Nyy.*(tempk))));
                indices{lo}=tempk;
                alkiot{lo}=single(templ_ijk(tempi).*(1/temp))*exp(-vaimennuskuva(indices{lo})'*templ_ijk(tempi));
                lor(lo,:)=[lo,int32(length(alkiot{lo}))];
                
                continue
            
            end
            discard(lo)=false;
            continue
        end
        
        % Parametriesitykseen tarvittavat t-arvot
        
        tx=((bx+iij.*d)-xs)./x_diff;
        
        ty=((by+jji.*d)-ys)./y_diff;
        
        tz=((bz+kkj.*dz)-zs)./z_diff;
        
        
        % Etsit��n minimit ja maksimit, eli ts. t:n arvot joilla s�de tulee ja
        % poistuu kuvasta x-akselilla ja y-akselilla
        txmin=min(tx(1),tx(end));
        txmax=max(tx(1),tx(end));
        tymin=min(ty(1),ty(end));
        tymax=max(ty(1),ty(end));
        tzmin=min(tz(1),tz(end));
        tzmax=max(tz(1),tz(end));
        
        % Yhteiset arvot
        tmin=max([txmin;tymin;tzmin]);
        tmax=min([txmax;tymax;tzmax]);
        
%         if lo == 13103601
%             tmin
%             tmax
%         end
        
        if tmin >= tmax
           %disp('skip')
            discard(lo)=false;
            continue
        end
        
        % Ensiksi tutkitaan onko l�hde "alempana" vai "ylemp�n�" kuin detektori,
        % PET:ss� l�hdekin on detektori
        % i kuvastaa x-akselin pikseliviivova, j y-akselin
        if xs<xd % Ensiksi jos l�hde on alempana
            if tmin==txmin % Jos yhteinen minimi t on sama kuin x-akselin parametrisoidun t niin tulopiste on i=1
                imin=int32(1);
            else % Muulloin lasketaan pisteen arvo pienimm�ss� t:n arvossa (esim. leikkaus y-akselilla), v�hennet��n
                % et�isyys origosta kuvaan ja jaetaan pikselien v�lisell�
                % et�isyydell�
                pxt=xs+tmin*x_diff;
                imin=int32(ceil((pxt-bx)/d));
            end
            if tmax==txmax % Sama poistumistielle
                imax=int32(Nx - 1);
            else
                pxt=xs+tmax*x_diff;
                imax=int32(floor((pxt-bx)/d));
            end
        elseif xs>xd % Jos l�hde on ylemp�n�
            if tmin==txmin
                imax=int32(Nx-2);
            else
                pxt=xs+tmin*x_diff;
                imax=int32(floor((pxt-bx)/d));
            end
            if tmax==txmax
                imin=int32(0);
            else
                pxt=xs+tmax*x_diff;
                imin=int32(ceil((pxt-bx)/d));
            end
        end
        
        % Sama y-akselille
        if ys<yd
            if tmin==tymin
                jmin=int32(1);
            else
                pyt=ys+tmin*y_diff;
                jmin=int32(ceil((pyt-by)/d));
            end
            if tmax==tymax
                jmax=int32(Ny - 1);
            else
                pyt=ys+tmax*y_diff;
                jmax=int32(floor((pyt-by)/d));
            end
        elseif ys>yd
            if tmin==tymin
                jmax=int32(Ny-2);
            else
                pyt=ys+tmin*y_diff;
                jmax=int32(floor((pyt-by)/d));
            end
            if tmax==tymax
                jmin=int32(0);
            else
                pyt=ys+tmax*y_diff;
                jmin=int32(ceil((pyt-by)/d));
            end
        end
        
        % Sama z-akselille
        if zs<zd
            if tmin==tzmin
                kmin=int32(1);
            else
                pzt=zs+tmin*z_diff;
                kmin=int32(ceil((pzt-bz)/dz));
            end
            if tmax==tzmax
                kmax=int32(Nz - 1);
            else
                pzt=zs+tmax*z_diff;
                kmax=int32(floor((pzt-bz)/dz));
            end
        elseif zs>zd
            if tmin==tzmin
                kmax=int32(Nz-2);
            else
                pzt=zs+tmin*z_diff;
                kmax=int32(floor((pzt-bz)/dz));
            end
            if tmax==tzmax
                kmin=int32(0);
            else
                pzt=zs+tmax*z_diff;
                kmin=int32(ceil((pzt-bz)/dz));
            end
        end
        
        % Pikselien numerot, joissa s�de kulkee (+1 lis�ttty, koska Matlabissa ei
        % voi olla 0 alkavia)
        
        
        % Otetaan vain ne t:n arvot, jotka ovat s�teen kulkureitill�
        if xs<xd
            apu=(imin:1:imax)+1;
        elseif xs>xd
            apu=(imax:-1:imin)+1;
        end
        tx_n=tx((apu));
        
        if ys<yd
            apu=(jmin:1:jmax)+1;
        elseif ys>yd
            apu=(jmax:-1:jmin)+1;
        end
        ty_n=ty((apu));
        
        if zs<zd
            apu=(kmin:1:kmax)+1;
        elseif zs>zd
            apu=(kmax:-1:kmin)+1;
        end
        tz_n=tz((apu));
        
        % Koko s�teen pituus
        L=sqrt(x_diff^2+y_diff^2+z_diff^2);
        
        % Pistet��n parametri t:n arvot pienemm�st� suurempaan
%        tt=sort([tmin;tx_n';ty_n';tz_n']);
        % Poistetaan tuplat
        tt=unique([tmin;tx_n';ty_n';tz_n']);
        
%         if lo == 13103601
%             tt
%         end
        
        % Pisteet joiden l�pi s�de kulkee (riste�� akselien kanssa)
        
        pxt=xs+((tt(2:end)+tt(1:end-1))./2).*x_diff;
        tempi=int32(floor((pxt-bx)./d)+1);
        pyt=ys+((tt(2:end)+tt(1:end-1))./2).*y_diff;
        tempj=int32(floor((pyt-by)./d));
        pzt=zs+((tt(2:end)+tt(1:end-1))./2).*z_diff;
        tempk=int32(floor((pzt-bz)./dz));
        
        
        % Lasketaan kuljettu matka jokaisessa pikseliss�
        
        templ_ijk=(diff(tt)).*L;
        
        

        temp=sum(templ_ijk);
        [tempk, tempi] = sort(((tempj).*Nyy+tempi+(Nxx.*Nyy.*(tempk))));
%         if lo == 13103601
%             tempk
%         end
        
        indices{lo}=tempk;
%         apu = indices{lo}
%         joku=templ_ijk(tempi)
%         atten= (-vaimennuskuva(apu(1))'*joku(1))
%         atten_vali = (-vaimennuskuva(apu(2))'*joku(2))
%         atten= atten + (-vaimennuskuva(apu(2))'*joku(2))
%         jelppi1=(-vaimennuskuva(indices{lo})'*templ_ijk)
%         jelppi2=(1/temp)*exp(-vaimennuskuva(indices{lo})'*templ_ijk)
%         templ_ijk(tempi).*(1/temp)*exp(-vaimennuskuva(indices{lo})'*templ_ijk(tempi))
        alkiot{lo}=single(templ_ijk(tempi).*(1/temp))*exp(-vaimennuskuva(indices{lo})'*templ_ijk(tempi));
        
%         alkiot{lo}=single(templ_ijk(tempi).*(1/temp));
        lor(lo,:)=[lo,int32(length(alkiot{lo}))];
%         if length(indices{lo})>3
%         return
%         end
    end
end


end

