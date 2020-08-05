function [ lor, indices, alkiot, discard] = improved_siddon_atten( Ny, Nx, Nz, d, dz, by, bx, bz, z_det, x, y, yy, xx , NSinos, NSlices, vaimennus, index, pituus, attenuation_correction)
%IMPROVED_SIDDON_ATTEN Function to calculate the Siddon's algorithm
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
% LOOP_VAR_PAR is the length of the subset used for the current
% subiteration
% PIKSELIKOKO defines the size of the image matrix in x and y directions
% NY is the length of the y direction
% NX same for x
% NZ same for z
% D defines the square (pixel) length, that is the distance between pixel
% boundaries in x and y directions
% DZ same for z direction
% BY the distance from the origin to the pixel space in y direction
% BX same for x
% BZ same for z
% Z_DET detector ring locations in z dimension (axial)
% X detector ring locations in x direction (radial)
% Y same for y
% YY coordinates of the pixel boundaries in y-direction
% XX same for x
% NSINOS Number of sinograms included
% NSLICES Number of image slices included
% VAIMENNUS Attenuation images
% INDEX LOR indices to be used
% PITUUS Number of LORs
bxb = bx + double(Nx) * d;
byb = by + double(Ny) * d;
bzb = bz + double(Nz) * dz;
z_det = z_det(1:NSinos,:);
size_x = size(x,1);
% loop_var_par = NSinos*size(x,1);
indices=cell((pituus),1);
alkiot=cell((pituus),1);
% indices=zeros(100,(loop_var_par)/100,'int32');
% alkiot=zeros(100,(loop_var_par)/100,'int32');
lor=zeros((pituus),2,'int32');
discard=true((pituus),1);
% apu_var=osajoukot(test);
vaimennuskuva = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%vaimennusasiat:
if attenuation_correction
    vaimennuskuva = vaimennus(:,:,1:NSlices);
    vaimennuskuva = vaimennuskuva(:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 436142
%1587419
% ll=int32(1:size(x,1));
% ll = repmat(ll,1,NSinos);
% lz = int16(1:NSinos);
% lz = repelem(lz,size(x,1));
zmax = max(max(z_det));
if zmax==0
    zmax = 1;
end
parfor lo=int32(1:pituus)
    imax=int32(0);
    imin=int32(0);
    jmax=int32(0);
    jmin=int32(0);
    
    ju = 0;
    iu = 0;
    ku = 0;

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
    
    ll = mod(index(lo), uint32(size_x));
    if ll == 0
        ll = size_x;
    end
    lz = idivide(index(lo), uint32(size_x), 'ceil');
    xs=x(ll,1);
    xd=x(ll,2);
    ys=y(ll,1);
    yd=y(ll,2);
    zs = z_det(lz,1);
    zd = z_det(lz,2);
    z_loop=int32((zs/zmax)*(NSlices-1)+1);
    x_diff = xd-xs;
    y_diff = yd-ys;
    z_diff = zd-zs;
    
    if zs==zd
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
%         discard(lo)=false;
%         continue
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if abs(y_diff)<1e-8
            
            if yd<= max(yy) && yd>=min(yy)
                tempi=(1:Ny)';
                [~,apu]=min(abs(yy(1:Ny)-yd));
                tempj = ones(Ny,1,'int32')*(int32(apu));
                templ_ijk = ones(Ny,1,'double')*d;
                tempk = (z_loop)*ones(size(tempi),'int32');
                temp=sum(templ_ijk);
                [tempk, tempi] = sort(((tempj-1).*Ny+tempi+(Nx.*Ny.*(tempk-1))));
                
                indices{lo}=tempk;
                if attenuation_correction
                    alkiot{lo}=double(templ_ijk(tempi).*(1/temp))*exp(-vaimennuskuva(indices{lo})'*templ_ijk(tempi));
                else
                    alkiot{lo}=double(templ_ijk(tempi).*(1/temp));
                end
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
                tempj=(1:Nx)';
                [~,apu]=min(abs(xx(1:Nx)-xd));
                tempi=ones(Nx,1,'int32')*(int32(apu));
                templ_ijk = ones(Nx,1,'double')*d;
                tempk=(z_loop)*ones(size(tempj),'int32');
                temp=sum(templ_ijk);
                [tempk, tempi] = sort(((tempj-1).*Ny+tempi+(Nx.*Ny.*(tempk-1))));
                
                indices{lo}=tempk;
                if attenuation_correction
                    alkiot{lo}=double(templ_ijk(tempi).*(1/temp))*exp(-vaimennuskuva(indices{lo})'*templ_ijk(tempi));
                else
                    alkiot{lo}=double(templ_ijk(tempi).*(1/temp));
                end
                lor(lo,:)=[lo,int32(length(alkiot{lo}))];
                continue
            end
            discard(lo)=false;
            continue
        end
        % Parametriesitykseen tarvittavat t-arvot
        
        tx0=(bx-xs)/x_diff;
        txback=(bxb-xs)/x_diff;
        
        % Sama y:lle
        
        ty0=(by-ys)/y_diff;
        tyback=(byb-ys)/y_diff;
        txmin=min(tx0,txback);
        txmax=max(tx0,txback);
        tymin=min(ty0,tyback);
        tymax=max(ty0,tyback);
        
        % Yhteiset arvot
        tmin=max(txmin,tymin);
        tmax=min(txmax,tymax);
        
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
                imax=int32(Nx);
            else
                pxt=xs+tmax*x_diff;
                imax=int32(floor((pxt-bx)/d));
            end
            iu = 1;
            tx0=((bx+double(imin)*d)-xs)/x_diff;
%             apu=(imin:1:imax)+1;
%             tx_n=tx((apu));
        elseif xs>xd % Jos l�hde on ylemp�n�
            if tmin==txmin
                imax=int32(Nx-1);
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
            iu = -1;
            tx0=((bx+double(imax)*d)-xs)/x_diff;
%             apu=(imax:-1:imin)+1;
%             tx_n=tx((apu));
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
                jmax=int32(Ny);
            else
                pyt=ys+tmax*y_diff;
                jmax=int32(floor((pyt-by)/d));
            end
            ju = 1;
            ty0=((by+double(jmin)*d)-ys)/y_diff;
        elseif ys>yd
            if tmin==tymin
                jmax=int32(Ny-1);
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
            ju = -1;
            ty0=((by+double(jmax)*d)-ys)/y_diff;
        end
        
        Np = (imax - imin + 1) + (jmax - jmin + 1);
        
        % Otetaan vain ne t:n arvot, jotka ovat s�teen kulkureitill� (ei käytössä enää)
        

        % Koko s�teen pituus
        L=sqrt(x_diff^2+y_diff^2);
        
        % Pistet��n parametri t:n arvot pienemm�st� suurempaan
%         tt=sort([tmin;tx_n';ty_n']);
        % Poistetaan tuplat
        pt = (min(tx0, ty0) + tmin) /2;
        
        % Pisteet joiden l�pi s�de kulkee (riste�� akselien kanssa)
        
        tempi = int32(floor(((xs + pt * x_diff) - bx) / d) + 1);
        tempj = int32(floor(((ys + pt * y_diff) - by) / d));
        
        txu = d / abs(x_diff);
        tyu = d / abs(y_diff);
        
        tempk=(z_loop - 1) * Nx * Ny;
        
        tc = tmin;
        
        % Lasketaan kuljettu matka jokaisessa pikseliss�
        
        templ_ijk=zeros(Np,1,'double');
        idx=zeros(Np,1,'int32');
        hpk = 0;
       
        
        for ii = 1 : Np
            
            if tx0 < ty0
                apu = (tx0 - tc);
                idx(ii) = (tempj)*Nx+tempi+tempk;
                templ_ijk(ii) = apu * L;
                tempi = tempi + iu;
                tc = tx0;
                tx0 = tx0 + txu;
            elseif ty0 <= tx0
                apu = (ty0 - tc);
                idx(ii) = (tempj)*Nx+tempi+tempk;
                templ_ijk(ii) = apu * L;
                tempj = tempj + ju;
                tc = ty0;
                ty0 = ty0 + tyu;
            end
            hpk = hpk + 1;
            if tempj < 0 || tempi < 1 || tempj >= Ny || tempi > Nx
                break
            end
            
        end
            
        templ_ijk = templ_ijk(1:hpk);
        idx = idx(1:hpk);
        
        [idx, tempi] = sort(idx);
        temp=sum(templ_ijk);
        indices{lo}=idx;
        if attenuation_correction
            alkiot{lo}=double(templ_ijk(tempi).*(1/temp))*exp(-vaimennuskuva(idx)'*templ_ijk(tempi));
        else
            alkiot{lo}=double(templ_ijk(tempi).*(1/temp));
        end
        lor(lo,:)=[lo,int32(length(idx))];
        
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

                tx0=(bx-xs)/x_diff;
                txback=(bxb-xs)/x_diff;
        
                % Sama y:lle
        
                tz0=(bz-zs)/z_diff;
                tzback=(bzb-zs)/z_diff;
                txmin=min(tx0,txback);
                txmax=max(tx0,txback);
                tzmin=min(tz0,tzback);
                tzmax=max(tz0,tzback);
                
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
                        imax=int32(Nx);
                    else
                        pxt=xs+tmax*x_diff;
                        imax=int32(floor((pxt-bx)/d));
                    end
                    tx0=((bx+double(imin)*d)-xs)/x_diff;
                    iu = 1;
                elseif xs>xd % Jos l�hde on ylemp�n�
                    if tmin==txmin
                        imax=int32(Nx-1);
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
                    tx0=((bx+double(imax)*d)-xs)/x_diff;
                    iu = -1;
                end
                
                if zs<zd
                    if tmin==tzmin
                        kmin=int32(1);
                    else
                        pzt=zs+tmin*z_diff;
                        kmin=int32(ceil((pzt-bz)/dz));
                    end
                    if tmax==tzmax
                        kmax=int32(Nz);
                    else
                        pzt=zs+tmax*z_diff;
                        kmax=int32(floor((pzt-bz)/dz));
                    end
                    tz0=((bz+double(kmin)*dz)-zs)/z_diff;
                    ku = 1;
                elseif zs>zd
                    if tmin==tzmin
                        kmax=int32(Nz-1);
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
                    tz0=((bz+double(kmax)*dz)-zs)/z_diff;
                    ku = -1;
                end
                
                Np = (imax - imin + 1) + (kmax - kmin + 1);
                
                L=sqrt(x_diff^2+z_diff^2);
                
                pt = (min(tx0, tz0) + tmin) /2;
                
                % Pisteet joiden l�pi s�de kulkee (riste�� akselien kanssa)
                
                tempi = int32(floor(((xs + pt * x_diff) - bx) / d) + 1);
                tempk = int32(floor(((zs + pt * z_diff) - bz) / dz));
                
                txu = d / abs(x_diff);
                tzu = dz / abs(z_diff);

                tc = tmin;
                
                [~,apu]=min(abs(yy(1:Ny)-yd));
                tempj=int32((apu - 1)*Ny);
                
                templ_ijk=zeros(Np,1,'double');
                idx=zeros(Np,1,'int32');
                hpk = 0;
                
                for ii = 1 : Np
                    if tx0 < tz0
                        idx(ii) = (tempj)+tempi+(Nx*Ny*(tempk));
                        templ_ijk(ii) = (tx0 - tc) * L;
                        tempi = tempi + iu;
                        tc = tx0;
                        tx0 = tx0 + txu;
                    else
                        idx(ii) = (tempj)+tempi+(Nx*Ny*(tempk));
                        templ_ijk(ii) = (tz0 - tc) * L;
                        tempk = tempk + ku;
                        tc = tz0;
                        tz0 = tz0 + tzu;
                    end
                    hpk = hpk + 1;
                    if tempk < 0 || tempi < 1 || tempk >= Nz || tempi > Nx
                        break
                    end
                end
                
                templ_ijk = templ_ijk(1:hpk);
                idx = idx(1:hpk);
                
                [idx, tempi] = sort(idx);
                
                temp=sum(templ_ijk);
                indices{lo}=idx;
                if attenuation_correction
                    alkiot{lo}=double(templ_ijk(tempi).*(1/temp))*exp(-vaimennuskuva(idx)'*templ_ijk(tempi));
                else
                    alkiot{lo}=double(templ_ijk(tempi).*(1/temp));
                end
                lor(lo,:)=[lo,int32(length(idx))];
                
                continue
            
            end
            discard(lo)=false;
            continue
        elseif abs(x_diff)<1e-8
            if xd<= max(xx) && xd>=min(xx)
                
                ty0=(by-ys)/y_diff;
                tyback=(byb-ys)/y_diff;
        
                % Sama y:lle
        
                tz0=(bz-zs)/z_diff;
                tzback=(bzb-zs)/z_diff;
                tymin=min(ty0,tyback);
                tymax=max(ty0,tyback);
                tzmin=min(tz0,tzback);
                tzmax=max(tz0,tzback);
                
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
                        jmax=int32(Ny);
                    else
                        pyt=ys+tmax*y_diff;
                        jmax=int32(floor((pyt-by)/d));
                    end
                    ty0=((by+double(jmin)*d)-ys)/y_diff;
                    ju = 1;
                elseif ys>yd % Jos l�hde on ylemp�n�
                    if tmin==tymin
                        jmax=int32(Ny-1);
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
                    ty0=((by+double(jmax)*d)-ys)/y_diff;
                    ju = -1;
                end
                
                if zs<zd
                    if tmin==tzmin
                        kmin=int32(1);
                    else
                        pzt=zs+tmin*z_diff;
                        kmin=int32(ceil((pzt-bz)/dz));
                    end
                    if tmax==tzmax
                        kmax=int32(Nz);
                    else
                        pzt=zs+tmax*z_diff;
                        kmax=int32(floor((pzt-bz)/dz));
                    end
                    tz0=((bz+double(kmin)*dz)-zs)/z_diff;
                    ku = 1;
                elseif zs>zd
                    if tmin==tzmin
                        kmax=int32(Nz-1);
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
                    tz0=((bz+double(kmax)*dz)-zs)/z_diff;
                    ku = -1;
                end
                
                Np = (kmax - kmin + 1) + (jmax - jmin + 1);
                
                L=sqrt(y_diff^2+z_diff^2);
                
                pt = (min(ty0, tz0) + tmin) /2;
                
                % Pisteet joiden l�pi s�de kulkee (riste�� akselien kanssa)
                
                tempj = int32(floor(((ys + pt * y_diff) - by) / d));
                tempk = int32(floor(((zs + pt * z_diff) - bz) / dz));
                
                tyu = d / abs(y_diff);
                tzu = dz / abs(z_diff);

                tc = tmin;
                
                [~,apu]=min(abs(xx(1:Nx)-xd));
                tempi=int32(apu);
                
                templ_ijk=zeros(Np,1,'double');
                idx=zeros(Np,1,'int32');
                hpk = 0;
                
                for ii = 1 : Np
                    if ty0 < tz0
                        idx(ii) = (tempj)*Ny+tempi+(Nx*Ny*(tempk));
                        templ_ijk(ii) = (ty0 - tc) * L;
                        tempj = tempj + ju;
                        tc = ty0;
                        ty0 = ty0 + tyu;
                    else
                        idx(ii) = (tempj)*Ny+tempi+(Nx*Ny*(tempk));
                        templ_ijk(ii) = (tz0 - tc) * L;
                        tempk = tempk + ku;
                        tc = tz0;
                        tz0 = tz0 + tzu;
                    end
                    hpk = hpk + 1;
                    if tempi < 1 || tempj < 0 || tempk < 0 || tempj >= Ny || tempi > Nx || tempz >= Nz
                        break
                    end
                end
                
                templ_ijk = templ_ijk(1:hpk);
                idx = idx(1:hpk);
                
                [idx, tempi] = sort(idx);
                
                temp=sum(templ_ijk);
                indices{lo}=idx;
                if attenuation_correction
                    alkiot{lo}=double(templ_ijk(tempi).*(1/temp))*exp(-vaimennuskuva(idx)'*templ_ijk(tempi));
                else
                    alkiot{lo}=double(templ_ijk(tempi).*(1/temp));
                end
                lor(lo,:)=[lo,int32(length(idx))];
                
                continue
            
            end
            discard(lo)=false;
            continue
        end
        
        % Parametriesitykseen tarvittavat t-arvot
        
        ty0=(by-ys)/y_diff;
        tyback=(byb-ys)/y_diff;
        tx0=(bx-xs)/x_diff;
        txback=(bxb-xs)/x_diff;
        
        % Sama y:lle
        
        tz0=(bz-zs)/z_diff;
        tzback=(bzb-zs)/z_diff;
        
        txmin=min(tx0,txback);
        txmax=max(tx0,txback);
        tymin=min(ty0,tyback);
        tymax=max(ty0,tyback);
        tzmin=min(tz0,tzback);
        tzmax=max(tz0,tzback);
        
        % Yhteiset arvot
        tmin=max([txmin;tymin;tzmin]);
        tmax=min([txmax;tymax;tzmax]);
        
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
                imax=int32(Nx);
            else
                pxt=xs+tmax*x_diff;
                imax=int32(floor((pxt-bx)/d));
            end
            tx0=((bx+double(imin)*d)-xs)/x_diff;
            iu = 1;
        elseif xs>xd % Jos l�hde on ylemp�n�
            if tmin==txmin
                imax=int32(Nx-1);
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
            tx0=((bx+double(imax)*d)-xs)/x_diff;
            iu = -1;
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
                jmax=int32(Ny);
            else
                pyt=ys+tmax*y_diff;
                jmax=int32(floor((pyt-by)/d));
            end
            ty0=((by+double(jmin)*d)-ys)/y_diff;
            ju = 1;
        elseif ys>yd
            if tmin==tymin
                jmax=int32(Ny-1);
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
            ty0=((by+double(jmax)*d)-ys)/y_diff;
            ju = -1;
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
                kmax=int32(Nz);
            else
                pzt=zs+tmax*z_diff;
                kmax=int32(floor((pzt-bz)/dz));
            end
            tz0=((bz+double(kmin)*dz)-zs)/z_diff;
            ku = 1;
        elseif zs>zd
            if tmin==tzmin
                kmax=int32(Nz-1);
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
            tz0=((bz+double(kmax)*dz)-zs)/z_diff;
            ku = -1;
        end
        
        Np = (imax - imin + 1) + (jmax - jmin + 1)+ (kmax - kmin + 1);
        
        % Koko s�teen pituus
        L=sqrt(x_diff^2+y_diff^2+z_diff^2);
        
        pt = (min(min(ty0, tz0), tx0) + tmin) /2;
        
        % Pisteet joiden l�pi s�de kulkee (riste�� akselien kanssa)
        
        tempi = int32(floor(((xs + pt * x_diff) - bx) / d) + 1);
        tempj = int32(floor(((ys + pt * y_diff) - by) / d));
        tempk = int32(floor(((zs + pt * z_diff) - bz) / dz));
        
        txu = d / abs(x_diff);
        tyu = d / abs(y_diff);
        tzu = dz / abs(z_diff);
        
        tc = tmin;
        
        templ_ijk=zeros(Np,1,'double');
        idx=zeros(Np,1,'int32');
        hpk = 0;
        
        for ii = 1 : Np
            if ty0 < tz0 && ty0 < tx0
                idx(ii) = (tempj)*Ny+tempi+(Nx*Ny*(tempk));
                templ_ijk(ii) = (ty0 - tc) * L;
                tempj = tempj + ju;
                tc = ty0;
                ty0 = ty0 + tyu;
            elseif tz0 < tx0
                idx(ii) = (tempj)*Ny+tempi+(Nx*Ny*(tempk));
                templ_ijk(ii) = (tz0 - tc) * L;
                tempk = tempk + ku;
                tc = tz0;
                tz0 = tz0 + tzu;
            else
                idx(ii) = (tempj)*Ny+tempi+(Nx*Ny*(tempk));
                templ_ijk(ii) = (tx0 - tc) * L;
                tempi = tempi + iu;
                tc = tx0;
                tx0 = tx0 + txu;
            end
            hpk = hpk + 1;
            if tempi < 1 || tempj < 0 || tempk < 0 || tempj >= Ny || tempi > Nx || tempz >= Nz
                break
            end
        end
        
        
        templ_ijk = templ_ijk(1:hpk);
        idx = idx(1:hpk);

        [idx, tempi] = sort(idx);
        
        temp=sum(templ_ijk);
        indices{lo}=idx;
        if attenuation_correction
            alkiot{lo}=double(templ_ijk(tempi).*(1/temp))*exp(-vaimennuskuva(idx)'*templ_ijk(tempi));
        else
            alkiot{lo}=double(templ_ijk(tempi).*(1/temp));
        end
        lor(lo,:)=[lo,int32(length(idx))];
    end
end


end

