function im_vectors = computeDeblurMLEM(im_vectors, options, iter, gaussK, Nx, Ny, Nz)
%COMPUTEDEBLURMLEM Computes the PSF deblurring phase for all selected MLEM
%algorithms

if options.save_iter
    iter_n = iter;
else
    iter_n = 0;
end
if options.MLEM
    im_vectors.MLEM(:, iter_n + 1) = deblur(im_vectors.MLEM(:, iter_n + 1), options, gaussK, Nx, Ny, Nz);
elseif options.MRP && options.OSL_MLEM
    im_vectors.MRP_MLEM(:, iter_n + 1) = deblur(im_vectors.MRP_MLEM(:, iter_n + 1), options, gaussK, Nx, Ny, Nz);
elseif options.quad && options.OSL_MLEM
    im_vectors.Quad_MLEM(:, iter_n + 1) = deblur(im_vectors.Quad_MLEM(:, iter_n + 1), options, gaussK, Nx, Ny, Nz);
elseif options.Huber && options.OSL_MLEM
    im_vectors.Huber_MLEM(:, iter_n + 1) = deblur(im_vectors.Huber_MLEM(:, iter_n + 1), options, gaussK, Nx, Ny, Nz);
elseif options.L && options.OSL_MLEM
    im_vectors.L_MLEM(:, iter_n + 1) = deblur(im_vectors.L_MLEM(:, iter_n + 1), options, gaussK, Nx, Ny, Nz);
elseif options.FMH && options.OSL_MLEM
    im_vectors.FMH_MLEM(:, iter_n + 1) = deblur(im_vectors.FMH_MLEM(:, iter_n + 1), options, gaussK, Nx, Ny, Nz);
elseif options.weighted_mean && options.OSL_MLEM
    im_vectors.Weighted_MLEM(:, iter_n + 1) = deblur(im_vectors.Weighted_MLEM(:, iter_n + 1), options, gaussK, Nx, Ny, Nz);
elseif options.TV && options.OSL_MLEM
    im_vectors.TV_MLEM(:, iter_n + 1) = deblur(im_vectors.TV_MLEM(:, iter_n + 1), options, gaussK, Nx, Ny, Nz);
elseif options.AD && options.OSL_MLEM
    im_vectors.AD_MLEM(:, iter_n + 1) = deblur(im_vectors.AD_MLEM(:, iter_n + 1), options, gaussK, Nx, Ny, Nz);
elseif options.APLS && options.OSL_MLEM
    im_vectors.APLS_MLEM(:, iter_n + 1) = deblur(im_vectors.APLS_MLEM(:, iter_n + 1), options, gaussK, Nx, Ny, Nz);
elseif options.TGV && options.OSL_MLEM
    im_vectors.TGV_MLEM(:, iter_n + 1) = deblur(im_vectors.TGV_MLEM(:, iter_n + 1), options, gaussK, Nx, Ny, Nz);
elseif options.NLM && options.OSL_MLEM
    im_vectors.NLM_MLEM(:, iter_n + 1) = deblur(im_vectors.NLM_MLEM(:, iter_n + 1), options, gaussK, Nx, Ny, Nz);
end