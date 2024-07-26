function options = NLMPrepass(options)
%NLMPREPASS Prepass phase for the NLM prior
%   Computes the NLM Gaussian weights

g_x = linspace(-options.Nlx, options.Nlx, 2*options.Nlx + 1)';
g_y = linspace(-options.Nly, options.Nly, 2*options.Nly + 1);
g_z = zeros(1,1,options.Nlz*2+1);
g_z(1,1,:) = linspace(-options.Nlz, options.Nlz, 2*options.Nlz + 1);
gaussian = gaussianKernel(g_x, g_y, g_z, options.NLM_gauss, options.NLM_gauss, options.NLM_gauss);
options.gaussianNLM = gaussian(:);
if options.NLM_use_anatomical
    if ischar(options.NLM_reference_image)
        apu = load(options.NLM_reference_image);
        variables = fieldnames(apu);
        options.NLM_ref = double(apu.(variables{1}));
    else
        options.NLM_ref = options.NLM_reference_image;
    end
    options.NLM_ref = reshape(options.NLM_ref, options.Nx(1), options.Ny(1), options.Nz(1));
    if options.implementation == 2 || options.implementation == 3
        options.NLM_ref = single(options.NLM_ref);
    end
end
if options.implementation == 2 || options.implementation == 3
    options.gaussianNLM = single(options.gaussianNLM);
%     options.sigma = single(options.sigma);
%     options.Nlx = uint32(options.Nlx);
%     options.Nly = uint32(options.Nly);
%     options.Nlz = uint32(options.Nlz);
%     options.Ndx = uint32(options.Ndx);
%     options.Ndy = uint32(options.Ndy);
%     options.Ndz = uint32(options.Ndz);
end
end