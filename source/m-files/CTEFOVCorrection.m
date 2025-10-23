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
if options.useEFOV
    if ~isfield(options, 'transaxialEFOV') || ~options.transaxialEFOV
        if ~isfield(options, 'axialEFOV') || ~options.axialEFOV
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
end
if options.useExtrapolation
    disp('Extrapolating the projections')
    PnAx = floor(size(options.SinM, 1) * extrapLengthAxial);
    PnTr = floor(size(options.SinM, 2) * extrapLengthTransaxial);
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
    disp('Extending the FOV')
    if ~(options.transaxialEFOV) && ~(options.axialEFOV)
        warning('Neither transaxial nor axial EFOV selected, but EFOV itself is selected! Setting axial EFOV to True!')
        options.axialEFOV = true;
    end
    if options.transaxialEFOV
        options.FOVxOrig = options.FOVa_x;
        options.FOVyOrig = options.FOVa_y;
        options.NxOrig = options.Nx;
        options.NyOrig = options.Ny;
        nTransaxial = floor(options.Nx * eFOVLengthTransaxial) * 2;
        options.Nx = options.Nx + nTransaxial;
        options.Ny = options.Ny + nTransaxial;
        options.FOVa_x = options.FOVa_x + options.FOVa_x/options.NxOrig * nTransaxial;
        options.FOVa_y = options.FOVa_y + options.FOVa_y/options.NyOrig * nTransaxial;
    else
        options.FOVxOrig = options.FOVa_x;
        options.FOVyOrig = options.FOVa_y;
        options.NxOrig = options.Nx;
        options.NyOrig = options.Ny;
    end
    if options.axialEFOV
        nAxial = floor(options.Nz * eFOVLengthAxial) * 2;
        options.NzOrig = options.Nz;
        options.Nz = options.Nz + nAxial;
        options.axialFOVOrig = options.axial_fov;
        options.axial_fov = options.axial_fov + options.axial_fov/options.NzOrig * nAxial;
    else
        options.axialFOVOrig = options.axial_fov;
        options.NzOrig = options.Nz;
    end
end
end