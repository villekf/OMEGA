function options = APLSPrepass(options)
%APLSPREPASS Precomputation phase for APLS
%   Precomputes necessary (anatomical prior) variables
apu = load(options.APLS_reference_image);
variables = fieldnames(apu);
alkuarvo = double(apu.(variables{1}));
if size(alkuarvo,2) == 1
    koko_apu = sqrt(length(alkuarvo)/double(options.Nz));
    if floor(koko_apu) ~= koko_apu
        error('Reference image has to be square')
    else
        alkuarvo = reshape(alkuarvo, koko_apu,koko_apu,options.Nz);
        if koko_apu ~= options.Nx || size(alkuarvo,3) ~= options.Nz
            if (license('test', 'image_toolbox') || exist('imresize3','file') == 2) && options.Nz > 1
                alkuarvo = imresize3(alkuarvo, [options.Nx, options.Ny, options.Nz]);
            elseif options.Nz > 1 && options.Nz ~= size(alkuarvo,3) && ~license('test', 'image_toolbox') && exist('imresize3','file') ~= 2
                error(['The reference image has different number of slices and no imresize3 was found. Resize the reference image manually ' ...
                    'to the same size as the reconstructed image.'])
            elseif options.Nz == 1 || options.Nz == size(alkuarvo,3) && ~license('test', 'image_toolbox') && exist('imresize3','file') ~= 2
                alkuarvo_new = zeros([options.Nx, options.Ny]);
                for kk = 1 : options.Nz
                    alkuarvo_new(:,:,kk) = imresize(alkuarvo(:,:,kk), [options.Nx, options.Ny]);
                end
                alkuarvo = alkuarvo_new;
            end
        end
    end
else
    if size(alkuarvo,2) ~= options.Ny || size(alkuarvo,3) ~= options.Nz
        if (license('test', 'image_toolbox') || exist('imresize3','file') == 2) && options.Nz > 1
            alkuarvo = imresize3(alkuarvo, [options.Nx, options.Ny, options.Nz]);
        elseif options.Nz > 1 && options.Nz ~= size(alkuarvo,3) && ~license('test', 'image_toolbox') && exist('imresize3','file') ~= 2
            error(['The reference image has different number of slices and no imresize3 was found. Resize the reference image manually ' ...
                'to the same size as the reconstructed image.'])
        elseif options.Nz == 1 || options.Nz == size(alkuarvo,3) && ~license('test', 'image_toolbox') && exist('imresize3','file') ~= 2
            alkuarvo_new = zeros([options.Nx, options.Ny]);
            for kk = 1 : options.Nz
                alkuarvo_new(:,:,kk) = imresize(alkuarvo(:,:,kk), [options.Nx, options.Ny]);
            end
            alkuarvo = alkuarvo_new;
        end
    end
end
% alkuarvo = alkuarvo - min(min(min(alkuarvo)));
% alkuarvo = alkuarvo/max(max(max(alkuarvo)));
if options.implementation == 2
    options.APLS_ref_image = single(alkuarvo);
    options.eta = single(options.eta);
    options.APLSsmoothing = single(options.APLSsmoothing);
else
    options.APLS_ref_image = double(alkuarvo);
end
end

