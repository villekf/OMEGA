function options = APLSPrepass(options)
%APLSPREPASS Precomputation phase for APLS
%   Precomputes necessary (anatomical prior) variables
% Size of the (first) reconstructed volume
Nx = options.Nx(1);
Ny = options.Ny(1);
Nz = options.Nz(1);
% The reference image can be input either as a file (mat-file) or as a variable
if ischar(options.APLS_reference_image) || isstring(options.APLS_reference_image)
    if strlength(string(options.APLS_reference_image)) == 0
        error('APLS selected, but no reference image provided!')
    end
    apu = load(options.APLS_reference_image);
    variables = fieldnames(apu);
    alkuarvo = double(apu.(variables{1}));
else
    alkuarvo = double(options.APLS_reference_image);
end
if isempty(alkuarvo)
    error('APLS selected, but no reference image provided!')
end
if size(alkuarvo,2) == 1
    koko_apu = sqrt(length(alkuarvo)/double(Nz));
    if floor(koko_apu) ~= koko_apu
        error('Reference image has to be square')
    else
        alkuarvo = reshape(alkuarvo, koko_apu,koko_apu,Nz);
        if koko_apu ~= Nx || size(alkuarvo,3) ~= Nz
            if (license('test', 'image_toolbox') || exist('imresize3','file') == 2) && Nz > 1
                alkuarvo = imresize3(alkuarvo, [Nx, Ny, Nz]);
            elseif Nz > 1 && Nz ~= size(alkuarvo,3) && ~license('test', 'image_toolbox') && exist('imresize3','file') ~= 2
                error(['The reference image has different number of slices and no imresize3 was found. Resize the reference image manually ' ...
                    'to the same size as the reconstructed image.'])
            elseif Nz == 1 || Nz == size(alkuarvo,3) && ~license('test', 'image_toolbox') && exist('imresize3','file') ~= 2
                alkuarvo_new = zeros([Nx, Ny]);
                for kk = 1 : Nz
                    alkuarvo_new(:,:,kk) = imresize(alkuarvo(:,:,kk), [Nx, Ny]);
                end
                alkuarvo = alkuarvo_new;
            end
        end
    end
else
    if size(alkuarvo,1) ~= Nx || size(alkuarvo,2) ~= Ny || size(alkuarvo,3) ~= Nz
        if (license('test', 'image_toolbox') || exist('imresize3','file') == 2) && Nz > 1
            alkuarvo = imresize3(alkuarvo, [Nx, Ny, Nz]);
        elseif Nz > 1 && Nz ~= size(alkuarvo,3) && ~license('test', 'image_toolbox') && exist('imresize3','file') ~= 2
            error(['The reference image has different number of slices and no imresize3 was found. Resize the reference image manually ' ...
                'to the same size as the reconstructed image.'])
        elseif Nz == 1 || Nz == size(alkuarvo,3) && ~license('test', 'image_toolbox') && exist('imresize3','file') ~= 2
            alkuarvo_new = zeros([Nx, Ny]);
            for kk = 1 : Nz
                alkuarvo_new(:,:,kk) = imresize(alkuarvo(:,:,kk), [Nx, Ny]);
            end
            alkuarvo = alkuarvo_new;
        end
    end
end
if numel(alkuarvo) ~= Nx * Ny * Nz
    error(['The APLS reference image has ' num2str(numel(alkuarvo)) ' elements after resizing, but the reconstructed image requires ' ...
        num2str(Nx * Ny * Nz) ' (' num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) '). Resize the reference image manually to the same size as the reconstructed image.'])
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

