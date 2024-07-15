function normalization_u = loadNormalizationBiograph(f_path, options)
%LOADNORMALIZATIONBIOGRAPH Loads the normalization coefficients created by
%Biograph machine, either the mCT or Vision
%   Detailed explanation goes here

fid = fopen(f_path);
normal = fread(fid, inf, 'single=>single',0,'l');
fclose(fid);
fid = fopen('normalization_test.n','w');
fwrite(fid, normal,'single');
fclose(fid);
if size(normal,1) <= 1255709
    mashing = options.det_w_pseudo / options.Nang / 2;
    det_w_pseudo = options.det_w_pseudo;
    no_pseudo = 0;
    if det_w_pseudo == options.det_per_ring
        det_w_pseudo = det_w_pseudo + 1;
        no_pseudo = 1;
    end
    normalization = ones(options.Ndist, options.Nang*mashing, options.TotSinos,'single');
    geom = reshape(normal(1:options.Ndist*options.Nz), options.Ndist, options.Nz);
    alku = options.Ndist*options.Nz + 1;
    interf = reshape(normal(alku:alku + (det_w_pseudo / options.blocks_per_ring) * options.Ndist - 1), ...
        det_w_pseudo / options.blocks_per_ring, options.Ndist);
    if no_pseudo
        interf = interf(1:options.cryst_per_block,:);
    end
    alku = alku + (det_w_pseudo / options.blocks_per_ring) * options.Ndist;
    ceff = reshape(normal(alku:alku + det_w_pseudo * options.rings - 1), ...
        det_w_pseudo, options.rings);
    ceff = 1./ceff;
    if no_pseudo
        ceff(options.cryst_per_block + 1 : options.cryst_per_block + 1 : end,:) = [];
    end
%     ceff = ceff / mean(ceff(:));
    alku = alku + det_w_pseudo * options.rings;
    axial = normal(alku:alku + options.TotSinos-1);
    alku = alku + options.TotSinos;
    parringDT = normal(alku:alku + options.rings-1);
    alku = alku + options.rings;
    nonparringDT = normal(alku:alku + options.rings-1);
    alku = alku + options.rings;
    TXcrystDT = normal(alku:end);
    z = sinogram_coordinates_3D(options);
    zz = z(1:options.Nz,1);
    z = z(options.Nz + 1: end,:);
    [~,loc] = ismember(z,zz);
    geom = [geom, (geom(:,loc(:,1)) + geom(:,loc(:,2))) / 2];
    geom = permute(geom, [1 3 2]);
    normalization = bsxfun(@rdivide, normalization, geom);
%     if mashing > 1
%         interf = reshape(accumarray(ceil((1:numel(interf))'/2),interf(:),[],@mean), (options.det_w_pseudo / options.blocks_per_ring) / 2, options.Ndist);
%     end
    interf = repmat(1./interf, options.det_w_pseudo / size(interf,1) / 2, 1)';
    %     interf = repmat(interf, options.blocks_per_ring, 1);
    axial = permute(axial, [3 2 1]);
    normalization = bsxfun(@times, normalization, axial);
    
    if mashing > 1
        options.Nang = options.Nang * mashing;
    end
    L = zeros(sum(1:options.det_w_pseudo),2,'int32');
    jh = int32(1);
    for kk = int32(1) : (options.det_w_pseudo)
        if exist('OCTAVE_VERSION','builtin') == 0 && exist('repelem', 'builtin') == 0
            L(jh:(jh + (options.det_w_pseudo) - kk),:) = [repeat_elem((kk), options.det_w_pseudo-(kk-1)), ((kk):options.det_w_pseudo)'];
        elseif exist('OCTAVE_VERSION','builtin') == 5
            L(jh:(jh + (options.det_w_pseudo) - kk),:) = [repelem((kk), options.det_w_pseudo-(kk-1)), ((kk):options.det_w_pseudo)'];
        else
            L(jh:(jh + (options.det_w_pseudo) - kk),:) = [repelem((kk), options.det_w_pseudo-(kk-1))', ((kk):options.det_w_pseudo)'];
        end
        jh = jh + (options.det_w_pseudo) -kk + 1;
    end
    L(L(:,1) == 0,:) = [];
    
    L = L - 1;
    
    xa = max(L,[],2);
    ya = min(L,[],2);
    
    j = idivide(mod(xa+ya+options.det_w_pseudo/2+1,options.det_w_pseudo),2);
    
    b = j+options.det_w_pseudo/2;
    
    i = abs(xa-ya-options.det_w_pseudo/2);
    for kk = 1 : length(ya)
        if (ya(kk)<j(kk)) || (b(kk)<xa(kk))
            i(kk) = -i(kk);
        end
    end
    
    % Determine the accepted LORs (distances that are within the predefined
    % value)
    if mod(options.Ndist,2) == 0
        accepted_lors = (i <= (options.Ndist/2 + min(0,options.ndist_side)) & i >= (-options.Ndist/2 + max(0,options.ndist_side)));
    else
        accepted_lors = (i <= options.Ndist/2 & i >= (-options.Ndist/2));
    end
    
    j = idivide(j,options.det_w_pseudo/2/options.Nang);
    
    i = i(accepted_lors);
    j = j(accepted_lors);
    if min(i) <= 0
        i = i + abs(min(i)) + 1;
    end
    j = j + 1;
    
    L = L(accepted_lors,:);
    
    L = L + 1;
%     if mashing > 1
%         options.Nang = options.Nang / mashing;
%     end
    
    ceff_n1 = ceff(L(:,1),:);
    ceff_n2 = ceff(L(:,2),:);
    cell_ceff1 = arrayfun(@(x) accumarray([i j],ceff_n1(:,x)),1:size(ceff_n1,2),'un',0);
    cell_ceff2 = arrayfun(@(x) accumarray([i j],ceff_n2(:,x)),1:size(ceff_n2,2),'un',0);
    ceff_nn = cellfun(@times, repmat(cell_ceff1,size(cell_ceff1,2),1), repmat(cell_ceff2',1,size(cell_ceff2,2)), 'UniformOutput', false);
%     if mashing > 1
%         apu = cell2mat(ceff_nn);
%         ceff_nn = squeeze(mean(reshape(apu',[mashing,size(apu,2)/mashing,size(apu,1)]),1))';
%         [r, c, ~] = size(ceff_nn);
%         ceff_nn = mat2cell(ceff_nn, options.Ndist * ones(1,r/options.Ndist), options.Nang * ones(1,c/options.Nang));
%     end
    ceff_nn = cat(3,ceff_nn{:});
    ceff_nn = bsxfun(@times, ceff_nn, interf);
    ceff = ones(options.Ndist,options.Nang,options.TotSinos,'single');
    
    %     interf_n1 = interf(L(:,1),:);
    %     interf_n2 = interf(L(:,2),:);
    %     interfx = accumarray([i j], interf_n1(:,1), [options.Ndist options.Nang],@mean, NaN);
    %     interfy = accumarray([i j], interf_n2, [options.Ndist options.Nang],@mean, NaN);
    %     normalization = bsxfun(@rdivide, normalization, interf');
    
    kkj = zeros(floor((options.ring_difference-ceil(options.span/2))/options.span) + 1, 1);
    for kk = 1 : floor((options.ring_difference-ceil(options.span/2))/options.span) + 1
        kkj(kk) = ceil(options.span/2) + options.span*(kk - 1);
    end
    offset2 = cumsum(options.segment_table);
    % Create the sinograms
    % First the detectors on the same ring
    ceff(:,:,1:2:options.Nz) = ceff_nn(:,:,1:options.rings+1:options.rings^2);
    % Then the detectors on adjacent rings
    mean_jh = zeros(options.TotSinos,1);
    mean_jh(1:2:options.Nz) = 1;
    for jh=1:floor(options.span/2)
        apu = ceff_nn(:,:,jh*options.rings+1:options.rings+1:options.rings^2);
        apu2 = ceff_nn(:,:,jh+1:options.rings+1:(options.rings-jh)*options.rings);
        loc = ismember(1:options.Nz,jh+1:2:offset2(1)-jh);
        mean_jh(loc) = mean_jh(loc) + 2;
        ceff(:,:,jh+1:2:offset2(1)-jh) = ceff(:,:,jh+1:2:offset2(1)-jh) .* apu .* apu2;
    end
    %     mean_jh(mean_jh == 0) = 1;
    %     mean_jh = permute(mean_jh, [3 2 1]);
    %     ceff(:,:,1:options.Nz) = bsxfun(@rdivide, ceff(:,:,1:options.Nz), mean_jh);
    % Lastly the rest of the detectors with the amount of combined LORs
    % specified with the span value
    %     mean_jh = zeros(options.TotSinos,1);
    for ih=1:floor(length(options.segment_table)/2)
        for jh=1:options.span
            apu = ceff_nn(:,:,(kkj(ih)+jh-1)*options.rings+1:options.rings+1:end);
            loc = ismember(1:options.TotSinos,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1);
            mean_jh(loc) = mean_jh(loc) + 1;
            ceff(:,:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) = ceff(:,:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) .* (apu);
            apu2 = ceff_nn(:,:,kkj(ih)+jh:options.rings+1:(options.rings-kkj(ih)-jh+1)*options.rings);
            loc = ismember(1:options.TotSinos,offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1);
            mean_jh(loc) = mean_jh(loc) + 1;
            ceff(:,:,offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1) = ceff(:,:,offset2(2*ih)+jh:2:offset2(2*ih + 1)-jh+1) .* (apu2);
        end
    end
    mean_jh(mean_jh == 0) = 1;
    mean_jh = permute(mean_jh, [3 2 1]);
    ceff = min(ceff(:));
    ceff = bsxfun(@rdivide, ceff, mean_jh);
%     normalization = normalization ./ ceff;
    normalization = normalization .* ceff;
%     fad = fopen('norm3d_00.a');
%     A = fread(fad, inf, 'single=>single');
%     normalization_u = reshape(A, options.Ndist, options.Nang, options.TotSinos);
    normalization_u = cell2mat(arrayfun(@(k) mean(normalization(:,k:min(size(normalization,1),k+2-1),:),2), 1:2:size(normalization,2), 'un', 0));
    %     normalization = circshift(normalization,1);
    %     normalization = 1 / normalization;
else
    det_w_pseudo = options.det_w_pseudo;
    no_pseudo = 0;
    if det_w_pseudo == options.det_per_ring
        det_w_pseudo = options.blocks_per_ring * (options.cryst_per_block + 1);
        no_pseudo = 1;
    end
    geom = reshape(normal(1:options.Ndist*(options.rings * options.rings)), options.Ndist, (options.rings * options.rings));
    alku = options.Ndist*(options.rings * options.rings) + 1;
    interf = reshape(normal(alku:alku + ceil(det_w_pseudo / options.blocks_per_ring) * options.Ndist - 1), ...
        ceil(det_w_pseudo / options.blocks_per_ring), options.Ndist);
    %     if no_pseudo
    %         interf = interf(1:options.cryst_per_block,:);
    %     end
    interf = repmat(interf, options.Nang / size(interf,1),1);
    alku = alku + ceil(det_w_pseudo / options.blocks_per_ring) * options.Ndist;
    ceff = reshape(normal(alku:alku + 798 * options.rings - 1), ...
        798, options.rings);
    if no_pseudo
        ceff(1 : options.cryst_per_block + 1 : end,:) = [];
    end
    alku = alku + 798 * options.rings;
    axial = normal(alku:alku + (options.rings * options.rings)-1);
    alku = alku + (options.rings * options.rings);
    parringDT = normal(alku:alku + 6000*42*options.rings-1);
    parringDT = reshape(parringDT, 6000, 42, options.rings);
    
    if options.span > 1
        normalization = ones(options.Ndist, options.Nang, size(axial,1),'single');
        L = zeros(sum(1:options.det_w_pseudo),2,'int32');
        jh = int32(1);
        for kk = int32(1) : (options.det_w_pseudo)
            if exist('OCTAVE_VERSION','builtin') == 0 && exist('repelem', 'builtin') == 0
                L(jh:(jh + (options.det_w_pseudo) - kk),:) = [repeat_elem((kk), options.det_w_pseudo-(kk-1)), ((kk):options.det_w_pseudo)'];
            elseif exist('OCTAVE_VERSION','builtin') == 5
                L(jh:(jh + (options.det_w_pseudo) - kk),:) = [repelem((kk), options.det_w_pseudo-(kk-1)), ((kk):options.det_w_pseudo)'];
            else
                L(jh:(jh + (options.det_w_pseudo) - kk),:) = [repelem((kk), options.det_w_pseudo-(kk-1))', ((kk):options.det_w_pseudo)'];
            end
            jh = jh + (options.det_w_pseudo) -kk + 1;
        end
        L(L(:,1) == 0,:) = [];
        
        L = L - 1;
        
        xa = max(L,[],2);
        ya = min(L,[],2);
        
        j = idivide(mod(xa+ya+options.det_w_pseudo/2+1,options.det_w_pseudo),2);
        
        b = j+options.det_w_pseudo/2;
        
        i = abs(xa-ya-options.det_w_pseudo/2);
        for kk = 1 : length(ya)
            if (ya(kk)<j(kk)) || (b(kk)<xa(kk))
                i(kk) = -i(kk);
            end
        end
        
        % Determine the accepted LORs (distances that are within the predefined
        % value)
        if mod(options.Ndist,2) == 0
            accepted_lors = (i <= (options.Ndist/2 + min(0,options.ndist_side)) & i >= (-options.Ndist/2 + max(0,options.ndist_side)));
        else
            accepted_lors = (i <= options.Ndist/2 & i >= (-options.Ndist/2));
        end
        
        j = idivide(j,options.det_w_pseudo/2/options.Nang);
        
        i = i(accepted_lors);
        j = j(accepted_lors);
        if min(i) <= 0
            i = i + abs(min(i)) + 1;
        end
        j = j + 1;
        
        L = L(accepted_lors,:);
        
        L = L + 1;
        %     if mashing > 1
        %         options.Nang = options.Nang / mashing;
        %     end
        
        ceff_n1 = ceff(L(:,1),:);
        ceff_n2 = ceff(L(:,2),:);
        cell_ceff1 = arrayfun(@(x) accumarray([i j],ceff_n1(:,x)),1:size(ceff_n1,2),'un',0);
        cell_ceff2 = arrayfun(@(x) accumarray([i j],ceff_n2(:,x)),1:size(ceff_n2,2),'un',0);
        ceff_nn = cellfun(@times, repmat(cell_ceff1,size(cell_ceff1,2),1), repmat(cell_ceff2',1,size(cell_ceff2,2)), 'UniformOutput', false);
        %     if mashing > 1
        %         apu = cell2mat(ceff_nn);
        %         ceff_nn = squeeze(mean(reshape(apu',[mashing,size(apu,2)/mashing,size(apu,1)]),1))';
        %         [r, c, ~] = size(ceff_nn);
        %         ceff_nn = mat2cell(ceff_nn, options.Ndist * ones(1,r/options.Ndist), options.Nang * ones(1,c/options.Nang));
        %     end
        ceff_nn = cat(3,ceff_nn{:});
        interf = circshift(interf, -1, 2);
        normalization = normalization .* permute(axial, [3 2 1]) .* ceff_nn .* permute(geom, [1 3 2]) .* interf';
        
        kkj = zeros(floor((options.ring_difference-ceil(options.span/2))/options.span) + 1, 1);
        for kk = 1 : floor((options.ring_difference-ceil(options.span/2))/options.span) + 1
            kkj(kk) = ceil(options.span/2) + options.span*(kk - 1);
        end
        offset2 = cumsum(options.segment_table);
        normalization_u = zeros(options.Ndist,options.Nang, options.TotSinos,'single');
        % Create the sinograms
        % First the detectors on the same ring
        normalization_u(:,:,1:2:options.Nz) = normalization(:,:,1:options.rings+1:options.rings^2);
        % Then the detectors on adjacent rings
        mean_jh = zeros(options.TotSinos,1);
        mean_jh(1:2:options.Nz) = 1;
        for jh=1:floor(options.span/2)
            apu = normalization(:,:,jh*options.rings+1:options.rings+1:options.rings^2);
            apu2 = normalization(:,:,jh+1:options.rings+1:(options.rings-jh)*options.rings);
            loc = ismember(1:options.Nz,jh+1:2:offset2(1)-jh);
            mean_jh(loc) = mean_jh(loc) + 2;
            normalization_u(:,:,jh+1:2:offset2(1)-jh) = normalization_u(:,:,jh+1:2:offset2(1)-jh) + apu + apu2;
        end
        %     mean_jh(mean_jh == 0) = 1;
        %     mean_jh = permute(mean_jh, [3 2 1]);
        %     ceff(:,:,1:options.Nz) = bsxfun(@rdivide, ceff(:,:,1:options.Nz), mean_jh);
        % Lastly the rest of the detectors with the amount of combined LORs
        % specified with the span value
        %     mean_jh = zeros(options.TotSinos,1);
        for ih=1:floor(length(options.segment_table)/2)
            for jh=1:options.span
                apu = normalization(:,:,(kkj(ih)+jh-1)*options.rings+1:options.rings+1:end);
                loc = ismember(1:options.TotSinos,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1);
                mean_jh(loc) = mean_jh(loc) + 1;
                normalization_u(:,:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) = normalization_u(:,:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) + (apu);
                apu2 = normalization(:,:,kkj(ih)+jh:options.rings+1:(options.rings-kkj(ih)-jh+1)*options.rings);
                loc = ismember(1:options.TotSinos,offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1);
                mean_jh(loc) = mean_jh(loc) + 1;
                normalization_u(:,:,offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1) = normalization_u(:,:,offset2(2*ih)+jh:2:offset2(2*ih + 1)-jh+1) + (apu2);
            end
        end
        mean_jh(mean_jh == 0) = 1;
        mean_jh = permute(mean_jh, [3 2 1]);
        normalization_u = bsxfun(@rdivide, normalization_u, mean_jh);
        normalization_u = circshift(normalization_u, [0 -1]);
        
        
%         kkj = zeros(floor((options.ring_difference-ceil(options.span/2))/options.span) + 1, 1);
%         for kk = 1 : floor((options.ring_difference-ceil(options.span/2))/options.span) + 1
%             kkj(kk) = ceil(options.span/2) + options.span*(kk - 1);
%         end
%         offset2 = cumsum(options.segment_table);
%         geom_u = zeros(options.Ndist,options.TotSinos,'single');
%         % Create the sinograms
%         % First the detectors on the same ring
%         geom_u(:,1:2:options.Nz) = geom(:,1:options.rings+1:options.rings^2);
%         % Then the detectors on adjacent rings
%         mean_jh = zeros(options.TotSinos,1);
%         mean_jh(1:2:options.Nz) = 1;
%         for jh=1:floor(options.span/2)
%             apu = geom(:,jh*options.rings+1:options.rings+1:options.rings^2);
%             apu2 = geom(:,jh+1:options.rings+1:(options.rings-jh)*options.rings);
%             loc = ismember(1:options.Nz,jh+1:2:offset2(1)-jh);
%             mean_jh(loc) = mean_jh(loc) + 2;
%             geom_u(:,jh+1:2:offset2(1)-jh) = geom_u(:,jh+1:2:offset2(1)-jh) + apu + apu2;
%         end
%         %     mean_jh(mean_jh == 0) = 1;
%         %     mean_jh = permute(mean_jh, [3 2 1]);
%         %     ceff(:,:,1:options.Nz) = bsxfun(@rdivide, ceff(:,:,1:options.Nz), mean_jh);
%         % Lastly the rest of the detectors with the amount of combined LORs
%         % specified with the span value
%         %     mean_jh = zeros(options.TotSinos,1);
%         for ih=1:floor(length(options.segment_table)/2)
%             for jh=1:options.span
%                 apu = geom(:,(kkj(ih)+jh-1)*options.rings+1:options.rings+1:end);
%                 loc = ismember(1:options.TotSinos,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1);
%                 mean_jh(loc) = mean_jh(loc) + 1;
%                 geom_u(:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) = geom_u(:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) + (apu);
%                 apu2 = geom(:,kkj(ih)+jh:options.rings+1:(options.rings-kkj(ih)-jh+1)*options.rings);
%                 loc = ismember(1:options.TotSinos,offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1);
%                 mean_jh(loc) = mean_jh(loc) + 1;
%                 geom_u(:,offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1) = geom_u(:,offset2(2*ih)+jh:2:offset2(2*ih + 1)-jh+1) + (apu2);
%             end
%         end
%         mean_jh(mean_jh == 0) = 1;
%         mean_jh = permute(mean_jh, [2 1]);
%         geom_u = bsxfun(@rdivide, geom_u, mean_jh);
%         
%         %     if mashing > 1
%         %         options.Nang = options.Nang * mashing;
%         %     end
%         L = zeros(sum(1:options.det_w_pseudo),2,'int32');
%         jh = int32(1);
%         for kk = int32(1) : (options.det_w_pseudo)
%             if exist('OCTAVE_VERSION','builtin') == 0 && exist('repelem', 'builtin') == 0
%                 L(jh:(jh + (options.det_w_pseudo) - kk),:) = [repeat_elem((kk), options.det_w_pseudo-(kk-1)), ((kk):options.det_w_pseudo)'];
%             elseif exist('OCTAVE_VERSION','builtin') == 5
%                 L(jh:(jh + (options.det_w_pseudo) - kk),:) = [repelem((kk), options.det_w_pseudo-(kk-1)), ((kk):options.det_w_pseudo)'];
%             else
%                 L(jh:(jh + (options.det_w_pseudo) - kk),:) = [repelem((kk), options.det_w_pseudo-(kk-1))', ((kk):options.det_w_pseudo)'];
%             end
%             jh = jh + (options.det_w_pseudo) -kk + 1;
%         end
%         L(L(:,1) == 0,:) = [];
%         
%         L = L - 1;
%         
%         xa = max(L,[],2);
%         ya = min(L,[],2);
%         
%         j = idivide(mod(xa+ya+options.det_w_pseudo/2+1,options.det_w_pseudo),2);
%         
%         b = j+options.det_w_pseudo/2;
%         
%         i = abs(xa-ya-options.det_w_pseudo/2);
%         for kk = 1 : length(ya)
%             if (ya(kk)<j(kk)) || (b(kk)<xa(kk))
%                 i(kk) = -i(kk);
%             end
%         end
%         
%         % Determine the accepted LORs (distances that are within the predefined
%         % value)
%         if mod(options.Ndist,2) == 0
%             accepted_lors = (i <= (options.Ndist/2 + min(0,options.ndist_side)) & i >= (-options.Ndist/2 + max(0,options.ndist_side)));
%         else
%             accepted_lors = (i <= options.Ndist/2 & i >= (-options.Ndist/2));
%         end
%         
%         j = idivide(j,options.det_w_pseudo/2/options.Nang);
%         
%         i = i(accepted_lors);
%         j = j(accepted_lors);
%         if min(i) <= 0
%             i = i + abs(min(i)) + 1;
%         end
%         j = j + 1;
%         
%         L = L(accepted_lors,:);
%         
%         L = L + 1;
%         %     if mashing > 1
%         %         options.Nang = options.Nang / mashing;
%         %     end
%         
%         ceff_n1 = ceff(L(:,1),:);
%         ceff_n2 = ceff(L(:,2),:);
%         cell_ceff1 = arrayfun(@(x) accumarray([i j],ceff_n1(:,x)),1:size(ceff_n1,2),'un',0);
%         cell_ceff2 = arrayfun(@(x) accumarray([i j],ceff_n2(:,x)),1:size(ceff_n2,2),'un',0);
%         ceff_nn = cellfun(@times, repmat(cell_ceff1,size(cell_ceff1,2),1), repmat(cell_ceff2',1,size(cell_ceff2,2)), 'UniformOutput', false);
%         %     if mashing > 1
%         %         apu = cell2mat(ceff_nn);
%         %         ceff_nn = squeeze(mean(reshape(apu',[mashing,size(apu,2)/mashing,size(apu,1)]),1))';
%         %         [r, c, ~] = size(ceff_nn);
%         %         ceff_nn = mat2cell(ceff_nn, options.Ndist * ones(1,r/options.Ndist), options.Nang * ones(1,c/options.Nang));
%         %     end
%         ceff_nn = cat(3,ceff_nn{:});
%         ceff = zeros(options.Ndist,options.Nang,options.TotSinos,'single');
%         
%         %     interf_n1 = interf(L(:,1),:);
%         %     interf_n2 = interf(L(:,2),:);
%         %     interfx = accumarray([i j], interf_n1(:,1), [options.Ndist options.Nang],@mean, NaN);
%         %     interfy = accumarray([i j], interf_n2, [options.Ndist options.Nang],@mean, NaN);
%         %     normalization = bsxfun(@times, normalization, interf');
%         
%         kkj = zeros(floor((options.ring_difference-ceil(options.span/2))/options.span) + 1, 1);
%         for kk = 1 : floor((options.ring_difference-ceil(options.span/2))/options.span) + 1
%             kkj(kk) = ceil(options.span/2) + options.span*(kk - 1);
%         end
%         offset2 = cumsum(options.segment_table);
%         % Create the sinograms
%         % First the detectors on the same ring
%         ceff(:,:,1:2:options.Nz) = ceff_nn(:,:,1:options.rings+1:options.rings^2);
%         % Then the detectors on adjacent rings
%         mean_jh = zeros(options.TotSinos,1);
%         mean_jh(1:2:options.Nz) = 1;
%         for jh=1:floor(options.span/2)
%             apu = ceff_nn(:,:,jh*options.rings+1:options.rings+1:options.rings^2);
%             apu2 = ceff_nn(:,:,jh+1:options.rings+1:(options.rings-jh)*options.rings);
%             loc = ismember(1:options.Nz,jh+1:2:offset2(1)-jh);
%             mean_jh(loc) = mean_jh(loc) + 2;
%             ceff(:,:,jh+1:2:offset2(1)-jh) = ceff(:,:,jh+1:2:offset2(1)-jh) + apu + apu2;
%         end
%         %     mean_jh(mean_jh == 0) = 1;
%         %     mean_jh = permute(mean_jh, [3 2 1]);
%         %     ceff(:,:,1:options.Nz) = bsxfun(@rdivide, ceff(:,:,1:options.Nz), mean_jh);
%         % Lastly the rest of the detectors with the amount of combined LORs
%         % specified with the span value
%         %     mean_jh = zeros(options.TotSinos,1);
%         for ih=1:floor(length(options.segment_table)/2)
%             for jh=1:options.span
%                 apu = ceff_nn(:,:,(kkj(ih)+jh-1)*options.rings+1:options.rings+1:end);
%                 loc = ismember(1:options.TotSinos,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1);
%                 mean_jh(loc) = mean_jh(loc) + 1;
%                 ceff(:,:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) = ceff(:,:,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) + (apu);
%                 apu2 = ceff_nn(:,:,kkj(ih)+jh:options.rings+1:(options.rings-kkj(ih)-jh+1)*options.rings);
%                 loc = ismember(1:options.TotSinos,offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1);
%                 mean_jh(loc) = mean_jh(loc) + 1;
%                 ceff(:,:,offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1) = ceff(:,:,offset2(2*ih)+jh:2:offset2(2*ih + 1)-jh+1) + (apu2);
%             end
%         end
%         mean_jh(mean_jh == 0) = 1;
%         mean_jh = permute(mean_jh, [3 2 1]);
%         ceff = bsxfun(@rdivide, ceff, mean_jh);
%         
%         kkj = zeros(floor((options.ring_difference-ceil(options.span/2))/options.span) + 1, 1);
%         for kk = 1 : floor((options.ring_difference-ceil(options.span/2))/options.span) + 1
%             kkj(kk) = ceil(options.span/2) + options.span*(kk - 1);
%         end
%         offset2 = cumsum(options.segment_table);
%         axial_u = zeros(options.TotSinos,1,'single');
%         % Create the sinograms
%         % First the detectors on the same ring
%         axial_u(1:2:options.Nz) = axial(1:options.rings+1:options.rings^2);
%         % Then the detectors on adjacent rings
%         mean_jh = zeros(options.TotSinos,1);
%         mean_jh(1:2:options.Nz) = 1;
%         for jh=1:floor(options.span/2)
%             apu = axial(jh*options.rings+1:options.rings+1:options.rings^2);
%             apu2 = axial(jh+1:options.rings+1:(options.rings-jh)*options.rings);
%             loc = ismember(1:options.Nz,jh+1:2:offset2(1)-jh);
%             mean_jh(loc) = mean_jh(loc) + 2;
%             axial_u(jh+1:2:offset2(1)-jh) = axial_u(jh+1:2:offset2(1)-jh) + apu + apu2;
%         end
%         %     mean_jh(mean_jh == 0) = 1;
%         %     mean_jh = permute(mean_jh, [3 2 1]);
%         %     ceff(:,:,1:options.Nz) = bsxfun(@rdivide, ceff(:,:,1:options.Nz), mean_jh);
%         % Lastly the rest of the detectors with the amount of combined LORs
%         % specified with the span value
%         %     mean_jh = zeros(options.TotSinos,1);
%         for ih=1:floor(length(options.segment_table)/2)
%             for jh=1:options.span
%                 apu = axial((kkj(ih)+jh-1)*options.rings+1:options.rings+1:end);
%                 loc = ismember(1:options.TotSinos,offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1);
%                 mean_jh(loc) = mean_jh(loc) + 1;
%                 axial_u(offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) = axial_u(offset2(2*ih-1)+jh:2:offset2(2*ih)-jh+1) + (apu);
%                 apu2 = axial(kkj(ih)+jh:options.rings+1:(options.rings-kkj(ih)-jh+1)*options.rings);
%                 loc = ismember(1:options.TotSinos,offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1);
%                 mean_jh(loc) = mean_jh(loc) + 1;
%                 axial_u(offset2(2*ih)+jh:2:offset2(2*ih+1)-jh+1) = axial_u(offset2(2*ih)+jh:2:offset2(2*ih + 1)-jh+1) + (apu2);
%             end
%         end
%         mean_jh(mean_jh == 0) = 1;
%         axial_u = bsxfun(@rdivide, axial_u, mean_jh);
%         
%         normalization = normalization .* permute(axial_u, [3 2 1]) .* ceff .* permute(geom_u, [1 3 2]) .* interf';
%         normalization = normalization .* ceff;
    else
        L = zeros(sum(1:options.det_w_pseudo),2,'int32');
        jh = int32(1);
        for kk = int32(1) : (options.det_w_pseudo)
            if exist('OCTAVE_VERSION','builtin') == 0 && exist('repelem', 'builtin') == 0
                L(jh:(jh + (options.det_w_pseudo) - kk),:) = [repeat_elem((kk), options.det_w_pseudo-(kk-1)), ((kk):options.det_w_pseudo)'];
            elseif exist('OCTAVE_VERSION','builtin') == 5
                L(jh:(jh + (options.det_w_pseudo) - kk),:) = [repelem((kk), options.det_w_pseudo-(kk-1)), ((kk):options.det_w_pseudo)'];
            else
                L(jh:(jh + (options.det_w_pseudo) - kk),:) = [repelem((kk), options.det_w_pseudo-(kk-1))', ((kk):options.det_w_pseudo)'];
            end
            jh = jh + (options.det_w_pseudo) -kk + 1;
        end
        L(L(:,1) == 0,:) = [];
        
        L = L - 1;
        
        xa = max(L,[],2);
        ya = min(L,[],2);
        
        j = idivide(mod(xa+ya+options.det_w_pseudo/2+1,options.det_w_pseudo),2);
        
        b = j+options.det_w_pseudo/2;
        
        i = abs(xa-ya-options.det_w_pseudo/2);
        for kk = 1 : length(ya)
            if (ya(kk)<j(kk)) || (b(kk)<xa(kk))
                i(kk) = -i(kk);
            end
        end
        
        % Determine the accepted LORs (distances that are within the predefined
        % value)
        if mod(options.Ndist,2) == 0
            accepted_lors = (i <= (options.Ndist/2 + min(0,options.ndist_side)) & i >= (-options.Ndist/2 + max(0,options.ndist_side)));
        else
            accepted_lors = (i <= options.Ndist/2 & i >= (-options.Ndist/2));
        end
        
        j = idivide(j,options.det_w_pseudo/2/options.Nang);
        
        i = i(accepted_lors);
        j = j(accepted_lors);
        if min(i) <= 0
            i = i + abs(min(i)) + 1;
        end
        j = j + 1;
        
        L = L(accepted_lors,:);
        
        L = L + 1;
        %     if mashing > 1
        %         options.Nang = options.Nang / mashing;
        %     end
        
        ceff_n1 = ceff(L(:,1),:);
        ceff_n2 = ceff(L(:,2),:);
        cell_ceff1 = arrayfun(@(x) accumarray([i j],ceff_n1(:,x)),1:size(ceff_n1,2),'un',0);
        cell_ceff2 = arrayfun(@(x) accumarray([i j],ceff_n2(:,x)),1:size(ceff_n2,2),'un',0);
        ceff_nn = cellfun(@times, repmat(cell_ceff1,size(cell_ceff1,2),1), repmat(cell_ceff2',1,size(cell_ceff2,2)), 'UniformOutput', false);
        %     if mashing > 1
        %         apu = cell2mat(ceff_nn);
        %         ceff_nn = squeeze(mean(reshape(apu',[mashing,size(apu,2)/mashing,size(apu,1)]),1))';
        %         [r, c, ~] = size(ceff_nn);
        %         ceff_nn = mat2cell(ceff_nn, options.Ndist * ones(1,r/options.Ndist), options.Nang * ones(1,c/options.Nang));
        %     end
        ceff_nn = cat(3,ceff_nn{:});
        normalization = ones(options.Ndist, options.Nang, options.rings^2,'single');
        normalization_u = normalization .* permute(geom, [1 3 2]) .* ceff_nn .* permute(axial, [3 2 1]);
    end
%     normalization = circshift(normalization,1);
%     normalization_u = 1 / normalization_u;
end
end

