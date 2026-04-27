function options = CTEFOVCorrection(options)
%CTEFOVCORRECTION Extended FOV correction for (CB)CT data
%   Detailed explanation goes here
if ~isfield(options, 'useExtrapolation')
    options.useExtrapolation = false;
end
if ~isfield(options, 'scatter_correction')
    options.scatter_correction = false;
end
if ~isfield(options, 'useEFOV')
    options.useEFOV = false;
end
if ~isfield(options, 'useInpaint')
    options.useInpaint = false;
end
if ~isfield(options, 'extrapLengthAxial')
    if ~isfield(options, 'extrapLength')
        extrapLengthAxial = 0.25;
    else
        extrapLengthAxial = options.extrapLength;
    end
else
    extrapLengthAxial = options.extrapLengthAxial;
end
if ~isfield(options, 'extrapLength')
    extrapLength = .25;
else
    extrapLength = options.extrapLength;
end
if ~isfield(options, 'extrapLengthTransaxial')
    extrapLengthTransaxial = extrapLength;
else
    extrapLengthTransaxial = options.extrapLengthTransaxial;
end
if ~isfield(options, 'eFOVLengthAxial')
    if ~isfield(options, 'eFOVLength')
        eFOVLengthAxial = 0.3;
    else
        eFOVLengthAxial = options.eFOVLength;
    end
else
    eFOVLengthAxial = options.eFOVLengthAxial;
end
if ~isfield(options, 'eFOVLength')
    eFOVLength = .4;
else
    eFOVLength = options.eFOVLength;
end
if ~isfield(options, 'eFOVLengthTransaxial')
    eFOVLengthTransaxial = eFOVLength;
else
    eFOVLengthTransaxial = options.eFOVLengthTransaxial;
end
if ~isfield(options, 'useExtrapolationWeighting')
    weighting = false;
else
    weighting = options.useExtrapolationWeighting;
end
if ~isfield(options, 'offsetCorrection')
    offset = false;
else
    offset = options.offsetCorrection;
end
if options.useExtrapolation
    if ~isfield(options, 'transaxialExtrapolation') || ~options.transaxialExtrapolation
        if ~isfield(options, 'axialExtrapolation') || ~options.axialExtrapolation
            warning('Neither transaxial nor axial extrapolation selected! Defaulting to axial extrapolation!')
            options.transaxialExtrapolation = false;
            options.axialExtrapolation = true;
        else
            options.transaxialExtrapolation = false;
        end
    elseif ~isfield(options, 'axialExtrapolation')
        options.axialExtrapolation = false;
    end    
    disp('Extrapolating the projections')
    PnTr = floor(size(options.SinM, 1) * extrapLengthTransaxial);
    PnAx = floor(size(options.SinM, 2) * extrapLengthAxial);
    if options.useInpaint
        Vq = zeros(size(options.SinM,1) + PnTr*2, size(options.SinM,2) + PnAx*2, options.nProjections);
        for kk = 1 : options.nProjections
            apu = zeros(size(options.SinM,1) + PnTr*2, size(options.SinM,2) + PnAx*2);
            apu(PnTr + 1 : size(options.SinM,1) + PnTr, PnAx + 1 : size(options.SinM,2) + PnAx) = log(options.flat ./ options.SinM(:,:,kk));
            testi = apu(21:end-20,21:end-20);
            testi(testi == 0) = NaN;
            apu(21:end-20,21:end-20) = testi;
            %             apuG = gpuArray(apu);
            %             tic
            testi2 = inpaint_nans(apu,1);
            %             toc
            testi2(testi2 < 0) = 0;
            testi2(PnTr + 1 : size(options.SinM,1) + PnTr, PnAx + 1 : size(options.SinM,2) + PnAx) = log(options.flat ./ options.SinM(:,:,kk));
            Vq(:,:,kk) = testi2;
            kk
        end
        Vq = options.flat ./ exp(Vq);
        options.SinM = Vq;
        % testi = ones(size(options.SinM, 1) + Pn, size(options.SinM,2), size(options.SinM,3));
        % testi(1+Pn:end,:,:) = options.SinM;
    else
        if options.transaxialExtrapolation
            if offset
                size1 = size(options.SinM,1) + PnTr;
            else
                size1 = size(options.SinM,1) + PnTr * 2;
            end
        else
            size1 = size(options.SinM,1);
        end
        if options.axialExtrapolation
            size2 = size(options.SinM,2) + PnAx * 2;
        else
            size2 = size(options.SinM,2);
        end
        erotus1 = size1 - size(options.SinM,1);
        erotus2 = size2 - size(options.SinM,2);
        newProj = zeros(size1, size2, size(options.SinM,3), class(options.SinM));
        if offset
            newProj(1 + erotus1 : size(options.SinM,1) + erotus1, 1 + erotus2 / 2 : size(options.SinM,2) + erotus2 / 2,:) = options.SinM;
        else
            newProj(1 + erotus1 / 2 : size(options.SinM,1) + erotus1 / 2, 1 + erotus2 / 2 : size(options.SinM,2) + erotus2 / 2,:) = options.SinM;
        end
        if options.transaxialExtrapolation
            if offset
                apu = repmat(options.SinM(1,:,:), erotus1, 1, 1);
            else
                apu = repmat(options.SinM(1,:,:), erotus1 / 2, 1, 1);
            end
            if weighting
                apu = log(single(options.flat) ./ apu);
                pituus = round(size(apu,1) / (6/6));
                pituus2 = size(apu,1) - pituus;
                % apu = apu .* sin(linspace(0, pi/2, size(apu,1)))';
                apu = apu .* [zeros(pituus2, 1, class(apu));log(linspace(1, exp(1), pituus))'];
                apu = single(options.flat) ./ exp(apu);
                % apu = flipud(options.SinM(1:erotus1/2,:,:));
            end

            if ~offset
                newProj(1: erotus1 / 2, 1 + erotus2 / 2 : size(options.SinM,2) + erotus2 / 2, :) = apu;
                apu = repmat(options.SinM(end,:,:), erotus1 / 2, 1, 1);
                if weighting
                    apu = log(single(options.flat) ./ apu);
                    % apu = apu .* sin(linspace(pi/2, 0, size(apu,1)))';
                    apu = apu .* [log(linspace(exp(1), 1, pituus))';zeros(pituus2, 1, class(apu))];
                    apu = single(options.flat) ./ exp(apu);
                    % apu = flipud(options.SinM(size(options.SinM,1) - erotus1/2 + 1:end,:,:));
                end
                newProj(size(options.SinM,1) + erotus1 / 2 + 1 : end, 1 + erotus2 / 2 : size(options.SinM,2) + erotus2 / 2, :) = apu;
            else
                newProj(1: erotus1, 1 + erotus2 / 2 : size(options.SinM,2) + erotus2 / 2, :) = apu;
            end
        end
        if options.axialExtrapolation
            apu = repmat(newProj(:,erotus2 / 2 + 1,:), 1, erotus2 / 2, 1);
            if weighting
                apu = log(single(options.flat) ./ apu);
                % apu = apu .* sin(linspace(0, pi/2, size(apu,2)));
                apu = apu .* log(linspace(1, exp(1), size(apu,2)));
                apu = single(options.flat) ./ exp(apu);
                % apu = fliplr(newProj(:,erotus2/2 + 1 : erotus2,:));
            end
            newProj(:, 1: erotus2 / 2, :) = apu;
            apu = repmat(newProj(:,size(options.SinM,2) + erotus2 / 2,:), 1, erotus2 / 2, 1);
            if weighting
                apu = log(single(options.flat) ./ apu);
                % apu = apu .* sin(linspace(pi/2, 0, size(apu,2)));
                apu = apu .* log(linspace(exp(1), 1, size(apu,2)));
                apu = single(options.flat) ./ exp(apu);
                % apu = fliplr(newProj(:,size(newProj,2) - erotus2 + 1:size(newProj,2) - erotus2/2,:));
            end
            newProj(:, size(options.SinM,2) + erotus2 / 2 + 1 : end, :) = apu;
        end
    end
    options.SinM = newProj;
    clear newProj
    options.nRowsDOrig = options.nRowsD;
    options.nColsDOrig = options.nColsD;
    options.nRowsD = size(options.SinM,1);
    options.nColsD = size(options.SinM,2);
    if options.scatter_correction && options.corrections_during_reconstruction
        newProj = zeros(size1, size2, size(options.ScatterC,3), class(options.ScatterC));
        newProj(1 + erotus1 / 2 : size(options.ScatterC,1) + erotus1 / 2, 1 + erotus2 / 2 : size(options.ScatterC,2) + erotus2 / 2,:) = options.ScatterC;
        if options.transaxialExtrapolation
            apu = repmat(options.ScatterC(1,:,:), erotus1 / 2, 1, 1);
            if weighting
                apu = log(single(options.flat) ./ apu);
                pituus = round(size(apu,1) / (6/6));
                pituus2 = size(apu,1) - pituus;
                % apu = apu .* sin(linspace(0, pi/2, size(apu,1)))';
                apu = apu .* [zeros(pituus2, 1, class(apu));log(linspace(1, exp(1), pituus))'];
                apu = single(options.flat) ./ exp(apu);
                % apu = flipud(options.SinM(1:erotus1/2,:,:));
            end
            newProj(1: erotus1 / 2, 1 + erotus2 / 2 : size(options.ScatterC,2) + erotus2 / 2, :) = apu;
            apu = repmat(options.ScatterC(end,:,:), erotus1 / 2, 1, 1);
            if weighting
                apu = log(single(options.flat) ./ apu);
                % apu = apu .* sin(linspace(pi/2, 0, size(apu,1)))';
                apu = apu .* [log(linspace(exp(1), 1, pituus))';zeros(pituus2, 1, class(apu))];
                apu = single(options.flat) ./ exp(apu);
                % apu = flipud(options.SinM(size(options.SinM,1) - erotus1/2 + 1:end,:,:));
            end
            newProj(size(options.ScatterC,1) + erotus1 / 2 + 1 : end, 1 + erotus2 / 2 : size(options.ScatterC,2) + erotus2 / 2, :) = apu;
        end
        if options.axialExtrapolation
            apu = repmat(newProj(:,erotus2 / 2 + 1,:), 1, erotus2 / 2, 1);
            if weighting
                apu = log(single(options.flat) ./ apu);
                % apu = apu .* sin(linspace(0, pi/2, size(apu,2)));
                apu = apu .* log(linspace(1, exp(1), size(apu,2)));
                apu = single(options.flat) ./ exp(apu);
                % apu = fliplr(newProj(:,erotus2/2 + 1 : erotus2,:));
            end
            newProj(:, 1: erotus2 / 2, :) = apu;
            apu = repmat(newProj(:,size(options.ScatterC,2) + erotus2 / 2,:), 1, erotus2 / 2, 1);
            if weighting
                apu = log(single(options.flat) ./ apu);
                % apu = apu .* sin(linspace(pi/2, 0, size(apu,2)));
                apu = apu .* log(linspace(exp(1), 1, size(apu,2)));
                apu = single(options.flat) ./ exp(apu);
                % apu = fliplr(newProj(:,size(newProj,2) - erotus2 + 1:size(newProj,2) - erotus2/2,:));
            end
            newProj(:, size(options.ScatterC,2) + erotus2 / 2 + 1 : end, :) = apu;
        end
        options.ScatterC = newProj;
    end
end
if options.useEFOV
    options.axialEFOV = false; % Check if axial EFOV is inside FOV (after shift). If is outside (in both directions), set axialEFOV to true
    FOVmin_z = -options.axial_fov / 2;
    FOVmax_z = options.axial_fov / 2;
    eFOVmin_z = -options.eFOVSize(3) / 2 + options.eFOVShift(3);
    eFOVmax_z = options.eFOVSize(3) / 2 + options.eFOVShift(3);
    if FOVmin_z < eFOVmin_z || FOVmax_z > eFOVmax_z % FOV not entirely inside eFOV
        warning('The high-resolution FOV is not entirely inside the extended FOV in z-direction. No extension will be performed in the axial direction.')
        options.eFOVShift(3) = 0; % No extension; no shift required
        options.eFOVSize(3) = options.axial_fov;
    else % FOV entirely inside eFOV; extension required
        options.axialEFOV = true;
    end

    options.transaxialEFOV = false; % Check if transaxial EFOV is inside FOV (after shift). If is outside (in both directions), set transaxialEFOV to true
    FOVmin_x = -options.FOVa_x / 2;
    FOVmax_x = options.FOVa_x / 2;
    eFOVmin_x = -options.eFOVSize(1) / 2 + options.eFOVShift(1);
    eFOVmax_x = options.eFOVSize(1) / 2 + options.eFOVShift(1);
    FOVmin_y = -options.FOVa_y / 2;
    FOVmax_y = options.FOVa_y / 2;
    eFOVmin_y = -options.eFOVSize(2) / 2 + options.eFOVShift(2);
    eFOVmax_y = options.eFOVSize(2) / 2 + options.eFOVShift(2);
    if FOVmin_x < eFOVmin_x || FOVmax_x > eFOVmax_x || FOVmin_y < eFOVmin_y || FOVmax_y > eFOVmax_y % FOV not entirely inside eFOV 
        warning('The high-resolution FOV is not entirely inside the extended FOV in xy-direction. No extension will be performed in the transaxial direction.')
        options.eFOVShift(1) = 0; % No extension; no shift required
        options.eFOVSize(1) = options.FOVa_x;
        options.eFOVShift(2) = 0; % No extension; no shift required
        options.eFOVSize(2) = options.FOVa_y;
    else % FOV entirely inside eFOV; extension required
        options.transaxialEFOV = true;
    end

    if ~(options.axialEFOV || options.transaxialEFOV)
        options.useEFOV = false;
        warning('FOV extension is not performed; turning off options.useEFOV')
    end
end

if options.useEFOV
    disp('Extending the FOV')
    options.FOVxOrig = options.FOVa_x;
    options.FOVyOrig = options.FOVa_y;
    options.axialFOVOrig = options.axial_fov;    
    options.NxOrig = options.Nx;
    options.NyOrig = options.Ny;
    options.NzOrig = options.Nz;

    if options.transaxialEFOV
        options.FOVa_x = options.eFOVSize(1);
        options.FOVa_y = options.eFOVSize(2);
        options.Nx = ceil(options.Nx * options.FOVa_x / options.FOVxOrig);
        options.Ny = ceil(options.Ny * options.FOVa_y / options.FOVyOrig);
    end
    if options.axialEFOV
        options.axial_fov = options.eFOVSize(3);
        options.Nz = ceil(options.Nz * options.axial_fov / options.axialFOVOrig);
    end

    dx = options.FOVa_x / options.Nx;
    dy = options.FOVa_y / options.Ny;
    dz = options.axial_fov / options.Nz;
    options.eFOVShift_Nx = round(options.eFOVShift(1) / dx, TieBreaker="fromzero"); % Add half this amount to one side and subtract from another
    options.eFOVShift_Ny = round(options.eFOVShift(2) / dy, TieBreaker="fromzero");
    options.eFOVShift_Nz = round(options.eFOVShift(3) / dz, TieBreaker="fromzero");
end
end