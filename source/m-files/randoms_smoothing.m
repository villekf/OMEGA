function randoms = randoms_smoothing(randoms,options)
%RANDOMS_SMOOTHING Performs a moving mean smoothing on randoms or scatter
%data
%   Input randoms/scatter sinogram/raw list-mode data is smoothed by using
%   a 7x7 moving mean smoothing. The parameters Ndx and Ndy can be used to
%   adjust the size of the mean window. Smoothing is done symmetrically by
%   taking the mirror images of the measurement data on all sides. Outputs
%   the smoothed sinogram/raw list-mode data.
if options.verbose
    disp('Beginning randoms/scatter smoothing')
end
Ndx = 7;
Ndy = 7;
Ndz = 0;
if Ndz == 0
    g = repmat(1 / (Ndx * Ndy), Ndx, Ndy);
else
    g = repmat(1 / (Ndx * Ndy * Ndz), Ndx, Ndy, Ndz);
end
if options.use_raw_data
    koko = options.det_per_ring*options.rings;
    L = uint32(find(tril(true(koko,koko), 0)));
    V = single(nonzeros(randoms));
    K = find(randoms);
    prompt = zeros(koko, koko, 'single');
    prompt(L(K)) = V;
    prompt = conv2(prompt, single(g));
    randoms = double(prompt(tril(true(options.detectors,options.detectors),0)));
else
    padd = padding(randoms,[floor(Ndx/2) floor(Ndy/2) floor(Ndz/2)]);
    randoms = convn(padd, g, 'valid');
    
end
if options.verbose
    disp('Smoothing complete')
end