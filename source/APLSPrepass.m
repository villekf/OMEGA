function options = APLSPrepass(options)
%APLSPREPASS Precomputation phase for APLS
%   Precomputes necessary (anatomical prior) variables
apu = load(options.APLS_reference_image);
variables = fields(apu);
alkuarvo = double(apu.(variables{1}));
if size(alkuarvo,2) == 1
    koko_apu = sqrt(length(alkuarvo)/double(options.Nz));
    if floor(koko_apu) ~= koko_apu
        error('Reference image has to be square')
    else
        alkuarvo = reshape(alkuarvo, koko_apu,koko_apu,options.Nz);
        if koko_apu ~= options.Nx
            alkuarvo = imresize(alkuarvo, options.Nx, options.Ny, options.Nz);
        end
    end
else
    if size(alkuarvo,2) ~= options.Nx
        alkuarvo = imresize(alkuarvo, options.Nx, options.Ny, options.Nz);
    end
end
alkuarvo = alkuarvo - min(min(min(alkuarvo)));
alkuarvo = alkuarvo/max(max(max(alkuarvo)));
if options.implementation == 2
    options.APLS_ref_image = single(alkuarvo);
    options.eta = single(options.eta);
    options.APLSsmoothing = single(options.APLSsmoothing);
else
    options.APLS_ref_image = double(alkuarvo);
end
options.TVtype = 4;
end

