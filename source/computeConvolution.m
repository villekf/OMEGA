function vec = computeConvolution(vec, options, Nx, Ny, Nz, gaussK)
%COMPUTECONVOLUTION Computes symmetrically padded 3D convolution for the
%input vector (image)
%   
vec = reshape(vec, Nx, Ny, Nz);
gaussK = reshape(gaussK, options.g_dim_x * 2 + 1, options.g_dim_y * 2 + 1, options.g_dim_z * 2 + 1);
vec = padding(vec, [options.g_dim_x options.g_dim_y options.g_dim_z]);
vec = convn(vec, gaussK, 'valid');
vec = vec(:);
end

