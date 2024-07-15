function im_vectors = ACOSEM_coefficient(im_vectors, options, Nx, Ny, Nz, dx, dz, dy, bx, by, bz, z_det, x, y, yy, xx, NSinos, NSlices, size_x, ...
    zmax, norm_input, SinD, pituus, osa_iter, attenuation_correction, normalization_correction, randoms_correction, scatter, scatter_input, lor_a_input, ...
    xy_index_input, z_index_input, L_input, pseudot, det_per_ring, use_raw_data, epps, uu, dc_z, x_center, y_center, z_center, dec, bmin, ...
    bmax, Vmax, V, gaussK, TOF, TOFSize, sigma_x, TOFCenter, list_mode_format)

%ACOSEM_COEFFICIENT Computes the scaling coefficient for ACOSEM
%   This function computes the scaling coefficient for ACOSEM when using
%   implementation 4

if exist('feature','builtin') == 5
    nCores = uint32(feature('numcores'));
else
    nCores = uint32(1);
end

if options.use_psf
    OSEM_apu = computeConvolution(im_vectors.OSEM_apu, options, Nx, Ny, Nz, gaussK);
else
    OSEM_apu = im_vectors.OSEM_apu;
end
if list_mode_format
    x = x(pituus(osa_iter) * 2 + 1 : pituus(osa_iter + 1) * 2);
    y = y(pituus(osa_iter) * 2 + 1 : pituus(osa_iter + 1) * 2);
    z_det = z_det(pituus(osa_iter) * 2 + 1 : pituus(osa_iter + 1) * 2);
    det_per_ring = uint32(numel(x)/2);
end
no_norm_acosem = false;
[~,rhs_acosem] = computeImplementation4(options,use_raw_data,randoms_correction, pituus,osa_iter, normalization_correction,...
    Nx, Ny, Nz, dx, dy, dz, bx, by, bz, x, y, z_det, xx, yy, size_x, NSinos, NSlices, zmax, attenuation_correction, pseudot, det_per_ring, ...
    TOF, TOFSize, sigma_x, TOFCenter, dec, nCores, L_input, lor_a_input, xy_index_input, z_index_input, epps, uu, OSEM_apu, no_norm_acosem, ...
    x_center, y_center, z_center, bmin, bmax, Vmax, V, scatter_input, norm_input, SinD, dc_z);
[im_vectors.OSEM_apu, ~] = ACOSEM_im(im_vectors.OSEM_apu, rhs_acosem, uu, options, [], [], [], [], [], []);
end

