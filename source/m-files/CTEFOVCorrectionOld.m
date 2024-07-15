function options = CTEFOVCorrection(options)
%CTEFOVCORRECTION Extended FOV correction for (CB)CT data
%   Detailed explanation goes here
if ~isfield(options, 'useExtrapolation')
    options.useExtrapolation = false;
end
if ~isfield(options, 'useEFOV')
    options.useEFOV = false;
end
if ~isfield(options, 'useInpaint')
    options.useInpaint = false;
end
if ~isfield(options, 'extrapLength')
    extrapLength = .1;
else
    extrapLength = options.extrapLength;
end
if options.useEFOV
    if ~isfield(options, 'transaxialEFOV')
        if ~isfield(options, 'axialEFOV')
            warning('Neither transaxial nor axial extended FOV selected! Defaulting to axial EFOV!')
            options.transaxialEFOV = false;
            options.axialEFOV = true;
        else
            options.transaxialEFOV = false;
        end
    elseif ~isfield(options, 'axialEFOV')
        options.axialEFOV = false;
    end
end
if options.useExtrapolation
    if ~isfield(options, 'transaxialExtrapolation')
        if ~isfield(options, 'axialExtrapolation')
            warning('Neither transaxial nor axial extrapolation selected! Defaulting to axial extrapolation!')
            options.transaxialExtrapolation = false;
            options.axialExtrapolation = true;
        else
            options.transaxialExtrapolation = false;
        end
    elseif ~isfield(options, 'axialExtrapolation')
        options.axialExtrapolation = false;
    end
end
if options.useExtrapolation
    disp('Extrapolating the projections')
    Pn = floor(size(options.SinM, 1) * extrapLength) * 2;
    if options.useInpaint
        Vq = zeros(size(options.SinM,1) + Pn*2, size(options.SinM,2) + Pn*2, options.nProjections);
        for kk = 1 : options.nProjections
            apu = zeros(size(options.SinM,1) + Pn*2, size(options.SinM,2) + Pn*2);
            apu(Pn + 1 : size(options.SinM,1) + Pn, Pn + 1 : size(options.SinM,2) + Pn) = log(options.flat ./ options.SinM(:,:,kk));
            testi = apu(21:end-20,21:end-20);
            testi(testi == 0) = NaN;
            apu(21:end-20,21:end-20) = testi;
%             apuG = gpuArray(apu);
%             tic
            testi2 = inpaint_nans(apu,1);
%             toc
            testi2(testi2 < 0) = 0;
            testi2(Pn + 1 : size(options.SinM,1) + Pn, Pn + 1 : size(options.SinM,2) + Pn) = log(options.flat ./ options.SinM(:,:,kk));
            Vq(:,:,kk) = testi2;
            kk
        end
        Vq = options.flat ./ exp(Vq);
        options.SinM = Vq;
        % testi = ones(size(options.SinM, 1) + Pn, size(options.SinM,2), size(options.SinM,3));
        % testi(1+Pn:end,:,:) = options.SinM;
    else
    if exist('griddedInterpolant','file') == 2
        if options.transaxialExtrapolation
            % [X,~] = ndgrid(1:size(options.SinM,1),1:size(options.SinM,2));
            X = (1:size(options.SinM,1))';
%             Y = meshgrid(1:size(options.SinM,2));
            F = griddedInterpolant(X,options.SinM,'previous');
%             F = griddedInterpolant(X, Y,options.SinM,'spline');
            % [XX,~] = ndgrid(1:size(options.SinM,1)+Pn,1:size(options.SinM,2)+Pn);
            XX = (1:size(options.SinM,1)+Pn)';
            Vq = F(XX);
            Vq = log(single(options.flat) ./ single(Vq));
            pituus = round(Pn * (8 / 8));
            kerroin = sin(linspace(pi/2, 0, pituus));
%             kerroin = log(linspace(exp(1),1, pituus));
            Vq(size(options.SinM,1) + 1 : size(options.SinM,1) + pituus,:,:) = Vq(size(options.SinM,1) + 1 : size(options.SinM,1) + pituus,:,:);
            Vq(size(options.SinM,1) + pituus + 1 : end,:,:) = 0;
            if options.scatter_correction && options.corrections_during_reconstruction
                FS = griddedInterpolant(X,options.ScatterC,'previous');
                Sq = FS(XX);
            end
            % [X,~] = ndgrid(1:size(Vq,1),1:size(Vq,2));
            X = (1:size(Vq,1))';
            F = griddedInterpolant(X(:,1),Vq,'next');
            % [XX,~] = ndgrid(-Pn+1:size(Vq,1),-Pn+1:size(Vq));
            XX = (-Pn+1:size(Vq,1))';
            Vq = F(XX);
%             kerroin = sin(linspace(pi/2, 0, Pn/2));
            Vq(Pn - pituus + 1 : Pn,:,:) = Vq(Pn - pituus + 1 : Pn,:,:) .* flipud(kerroin');
            Vq(1 : Pn - pituus,:,:) = 0;
            Vq = single(single(options.flat) ./ exp(Vq));
%             x = linspace(1,0,Pn-ceil(Pn*.1)).^2';
%             x = (linspace(1,0,Pn-ceil(Pn*.00))).^(0)';
            % x = flipud(-log(linspace(1,1e-1,Pn-ceil(Pn*.1))'));
            % x = sqrt(linspace(1,0,Pn-ceil(Pn*.1)) + linspace(1,0,Pn-ceil(Pn*.1)).^2)';
            % x = x ./ max(x);
%             x = exp(-linspace(0,4,Pn-ceil(Pn*.1)))';
%             x = [zeros(ceil(Pn*.00),1);flipud(x);ones(size(options.SinM,1),1);x;zeros(ceil(Pn*.00),1)];
%             Vq = 1 - (1 - Vq) .* x;
            options.SinM = Vq;
            if options.scatter_correction && options.corrections_during_reconstruction
                FS = griddedInterpolant(X(:,1),Sq,'next');
                Sq = FS(XX);
                options.ScatterC = Sq;
            end
        end
        if options.axialExtrapolation
            Pn = floor(size(options.SinM, 2) * extrapLength) * 2;
            % testi = ones(size(options.SinM, 1), size(options.SinM,2), size(options.SinM,3));
            % testi(options.ySize + Pn + 1 :end, : ,:) = repmat(linspace(1,min(Vq(end,:,100)), Pn)'.^2, 1, options.xSize, options.nProjections);
            % options.SinM = options.SinM ./ testi;
            % options.SinM = options.SinM ./ flipud(testi);
            % options.SinM(options.SinM > 1) = 1;
            % testi = ones(size(options.SinM, 1), size(options.SinM,2) + Pn, size(options.SinM,3));
            % testi(:,1+Pn:end,:) = options.SinM;
            % [X,Y] = ndgrid(1:size(options.SinM,1),1:size(options.SinM,2));
            X = (1:size(options.SinM,2))';
            XX = (1:size(options.SinM,2)+Pn)';
            F = griddedInterpolant(X,permute(options.SinM, [2 1 3]),'previous');
            % [XX,YY] = ndgrid(1:size(options.SinM,1)+Pn,1:size(options.SinM,2)+Pn);
            Vq = F(XX);
            Vq = log(single(options.flat) ./ single(Vq));
            pituus = round(Pn * (8 / 8));
%             kerroin = sin(linspace(pi/2, 0, pituus));
            kerroin = log(linspace(exp(1),1, pituus));
            Vq(size(options.SinM,2) + 1 : size(options.SinM,2) + pituus,:,:) = Vq(size(options.SinM,2) + 1 : size(options.SinM,2) + pituus,:,:) .* kerroin';
            Vq(size(options.SinM,2) + pituus + 1 : end,:,:) = 0;
            if options.scatter_correction && options.corrections_during_reconstruction
                FS = griddedInterpolant(X,permute(options.ScatterC, [2 1 3]),'previous');
                Sq = FS(XX);
                Sq(size(options.SinM,2) + 1 : size(options.SinM,2) + pituus,:,:) = Sq(size(options.SinM,2) + 1 : size(options.SinM,2) + pituus,:,:) .* kerroin';
                Sq(size(options.SinM,2) + pituus + 1 : end,:,:) = 0;
            end
            % [X,Y] = ndgrid(1:size(Vq,1),1:size(Vq,2));
            X = (1:size(Vq,1))';
            F = griddedInterpolant(X,Vq,'next');
            % [XX,YY] = ndgrid(-Pn+1:size(Vq,1),-Pn+1:size(Vq));
            XX = (-Pn+1:size(Vq,1))';
            Vq = F(XX);
            Vq(Pn - pituus + 1 : Pn,:,:) = Vq(Pn - pituus + 1 : Pn,:,:) .* flipud(kerroin');
            Vq(1 : Pn - pituus,:,:) = 0;
            Vq = single(single(options.flat) ./ exp(Vq));
%             x = linspace(1,0,Pn-ceil(Pn*.1)).^2';
            % x = flipud(-log(linspace(1,1e-1,Pn-ceil(Pn*.1))'));
%             x = (linspace(1,0,Pn-ceil(Pn*.00))).^(0)';
%             x = sqrt(linspace(1,0,Pn-ceil(Pn*.1)) + linspace(1,0,Pn-ceil(Pn*.1)).^2)';
            % x = x ./ max(x);
%             x = exp(-linspace(0,4,Pn-ceil(Pn*.1)))';
%             x = [zeros(ceil(Pn*.00),1);flipud(x);ones(size(options.SinM,2),1);x;zeros(ceil(Pn*.00),1)];
%             Vq = 1 - (1 - Vq) .* x;
            options.SinM = permute(Vq, [2 1 3]);
            if options.scatter_correction && options.corrections_during_reconstruction
                FS = griddedInterpolant(X,Sq,'next');
                Sq = FS(XX);
                Sq(Pn - pituus + 1 : Pn,:,:) = Sq(Pn - pituus + 1 : Pn,:,:) .* flipud(kerroin');
                Sq(1 : Pn - pituus,:,:) = 0;
%                 Sq = 1 - (1 - Sq) .* x;
                options.ScatterC = permute(Sq, [2 1 3]);
            end
            % testi = ones(size(options.SinM, 1), size(options.SinM,2), size(options.SinM,3));
            % testi(: ,options.xSize + Pn + 1 :end, :) = repmat(permute(linspace(1,min(Vq(end,:,100)), Pn)'.^2, [2 1]), size(testi,1), 1, options.nProjections);
            % options.SinM = options.SinM ./ testi;
            % options.SinM = options.SinM ./ fliplr(testi);
            % options.SinM(options.SinM > 1) = 1;
        end
    else
        if options.transaxialExtrapolation
            X = (1:size(options.SinM,1))';
            XX = (1:size(options.SinM,1)+Pn)';
            Vq = interp1(X,options.SinM,XX,'previous','extrap');
            X = (1:size(Vq,1))';
            XX = (-Pn+1:size(Vq,1))';
            Vq = interp1(X,Vq,XX,'next','extrap');
%             x = exp(-linspace(0,4,Pn-ceil(Pn*.1)))';
            x = (linspace(1,0,Pn-ceil(Pn*.00))).^(0)';
            x = [zeros(ceil(Pn*.1),1);flipud(x);ones(size(options.SinM,1),1);x;zeros(ceil(Pn*.1),1)];
            Vq = 1 - (1 - Vq) .* x;
            options.SinM = Vq;
        end
        if options.axialExtrapolation
            Pn = floor(size(options.SinM, 2) * extrapLength) * 2;
            X = (1:size(options.SinM,2))';
            XX = (1:size(options.SinM,2)+Pn)';
            Vq = interp1(X,permute(options.SinM, [2 1 3]),XX,'previous','extrap');
            X = (1:size(Vq,1))';
            XX = (-Pn+1:size(Vq,1))';
            Vq = interp1(X,Vq,XX,'next','extrap');
%             x = exp(-linspace(0,4,Pn-ceil(Pn*.1)))';
            x = (linspace(1,0,Pn-ceil(Pn*.00))).^(0)';
            x = [zeros(ceil(Pn*.1),1);flipud(x);ones(size(options.SinM,2),1);x;zeros(ceil(Pn*.1),1)];
            Vq = 1 - (1 - Vq) .* x;
            options.SinM = permute(Vq, [2 1 3]);
        end
    end
    end
    options.nRowsDOrig = options.nRowsD;
    options.nColsDOrig = options.nColsD;
    options.nRowsD = size(options.SinM,1);
    options.nColsD = size(options.SinM,2);
end
if options.useEFOV
    disp('Extending the FOV')
    eFOVLength = .4;
    if options.transaxialEFOV
        nTransaxial = floor(options.Nx * eFOVLength) * 2;
        options.NxOrig = options.Nx;
        options.NyOrig = options.Ny;
        options.Nx = options.Nx + nTransaxial;
        options.Ny = options.Ny + nTransaxial;
        options.FOVxOrig = options.FOVa_x;
        options.FOVyOrig = options.FOVa_y;
        options.FOVa_x = options.FOVa_x + options.FOVa_x/options.NxOrig * nTransaxial;
        options.FOVa_y = options.FOVa_y + options.FOVa_y/options.NyOrig * nTransaxial;
    end
    if options.axialEFOV
        nAxial = floor(options.Nz * eFOVLength) * 2;
        options.NzOrig = options.Nz;
        options.Nz = options.Nz + nAxial;
        options.axialFOVOrig = options.axial_fov;
%         options.oOffsetZ = options.oOffsetZ - (options.axial_fov/options.NzOrig * nAxial);
        options.axial_fov = options.axial_fov + options.axial_fov/options.NzOrig * nAxial;
    end
end
end