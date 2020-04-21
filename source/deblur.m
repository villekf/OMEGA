function vec = deblur(vec, options, iter, subsets, gaussK, Nx, Ny, Nz)
%DEBLUR Computes the deblur phase for the PSF reconstruction for the input
%image/vector. Performs symmetric padding.

if size(vec,2) == 1
    jelppi = reshape(vec, Nx, Ny, Nz);
    apu = reshape(vec, Nx, Ny, Nz);
else
    jelppi = vec;
    apu = vec;
end
apu = padding(apu, [options.g_dim_x options.g_dim_y options.g_dim_z]);
for kk = 1 : iter * subsets
    apu2 = convn(padding(jelppi, [options.g_dim_x options.g_dim_y options.g_dim_z]), gaussK, 'valid');
    apu2 = padding(apu2, [options.g_dim_x options.g_dim_y options.g_dim_z]);
    jelppi = jelppi .* convn(apu ./ apu2, gaussK, 'valid');
end
vec = jelppi(:);