function im_vectors = computeDeblur(im_vectors, options, iter, subsets, gaussK, Nx, Ny, Nz)
%COMPUTEDEBLUR Computes the PSF deblurring phase for all selected OS
%algorithms

if options.implementation == 4
    if options.osem
        im_vectors.OSEM(:, iter + 1) = deblur(im_vectors.OSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.mramla
        im_vectors.MRAMLA(:, iter + 1) = deblur(im_vectors.MRAMLA(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.ramla
        im_vectors.RAMLA(:, iter + 1) = deblur(im_vectors.RAMLA(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.rosem
        im_vectors.ROSEM(:, iter + 1) = deblur(im_vectors.ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.rbi
        im_vectors.RBI(:, iter + 1) = deblur(im_vectors.RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.drama
        im_vectors.DRAMA(:, iter + 1) = deblur(im_vectors.DRAMA(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.cosem
        im_vectors.COSEM(:, iter + 1) = deblur(im_vectors.COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.ecosem
        im_vectors.ECOSEM(:, iter + 1) = deblur(im_vectors.ECOSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.acosem
        im_vectors.ACOSEM(:, iter + 1) = deblur(im_vectors.ACOSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.MRP && options.OSL_OSEM
        im_vectors.MRP_OSL(:, iter + 1) = deblur(im_vectors.MRP_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.MRP && options.MBSREM
        im_vectors.MRP_MBSREM(:, iter + 1) = deblur(im_vectors.MRP_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.MRP && options.BSREM
        im_vectors.MRP_BSREM(:,iter+1) = deblur(im_vectors.MRP_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.MRP && options.ROSEM_MAP
        im_vectors.MRP_ROSEM(:,iter+1) = deblur(im_vectors.MRP_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.MRP && options.RBI_OSL
        im_vectors.MRP_RBI(:, iter + 1) = deblur(im_vectors.MRP_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.MRP && any(options.COSEM_OSL)
        im_vectors.MRP_COSEM(:, iter + 1) = deblur(im_vectors.MRP_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.quad && options.OSL_OSEM
        im_vectors.Quad_OSL(:, iter + 1) = deblur(im_vectors.Quad_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.quad && options.MBSREM
        im_vectors.Quad_MBSREM(:, iter + 1) = deblur(im_vectors.Quad_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.quad && options.BSREM
        im_vectors.Quad_BSREM(:,iter+1) = deblur(im_vectors.Quad_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.quad && options.ROSEM_MAP
        im_vectors.Quad_ROSEM(:,iter+1) = deblur(im_vectors.Quad_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.quad && options.RBI_OSL
        im_vectors.Quad_RBI(:, iter + 1) = deblur(im_vectors.Quad_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.quad && any(options.COSEM_OSL)
        im_vectors.Quad_COSEM(:, iter + 1) = deblur(im_vectors.Quad_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.Huber && options.OSL_OSEM
        im_vectors.Huber_OSL(:, iter + 1) = deblur(im_vectors.Huber_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.Huber && options.MBSREM
        im_vectors.Huber_MBSREM(:, iter + 1) = deblur(im_vectors.Huber_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.Huber && options.BSREM
        im_vectors.Huber_BSREM(:,iter+1) = deblur(im_vectors.Huber_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.Huber && options.ROSEM_MAP
        im_vectors.Huber_ROSEM(:,iter+1) = deblur(im_vectors.Huber_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.Huber && options.RBI_OSL
        im_vectors.Huber_RBI(:, iter + 1) = deblur(im_vectors.Huber_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.Huber && any(options.COSEM_OSL)
        im_vectors.Huber_COSEM(:, iter + 1) = deblur(im_vectors.Huber_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.L && options.OSL_OSEM
        im_vectors.L_OSL(:, iter + 1) = deblur(im_vectors.L_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.L && options.MBSREM
        im_vectors.L_MBSREM(:, iter + 1) = deblur(im_vectors.L_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.L && options.BSREM
        im_vectors.L_BSREM(:,iter+1) = deblur(im_vectors.L_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.L && options.ROSEM_MAP
        im_vectors.L_ROSEM(:,iter+1) = deblur(im_vectors.L_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.L && options.RBI_OSL
        im_vectors.L_RBI(:, iter + 1) = deblur(im_vectors.L_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.L && any(options.COSEM_OSL)
        im_vectors.L_COSEM(:, iter + 1) = deblur(im_vectors.L_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.FMH && options.OSL_OSEM
        im_vectors.FMH_OSL(:, iter + 1) = deblur(im_vectors.FMH_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.FMH && options.MBSREM
        im_vectors.FMH_MBSREM(:, iter + 1) = deblur(im_vectors.FMH_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.FMH && options.BSREM
        im_vectors.FMH_BSREM(:,iter+1) = deblur(im_vectors.FMH_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.FMH && options.ROSEM_MAP
        im_vectors.FMH_ROSEM(:,iter+1) = deblur(im_vectors.FMH_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.FMH && options.RBI_OSL
        im_vectors.FMH_RBI(:, iter + 1) = deblur(im_vectors.FMH_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.FMH && any(options.COSEM_OSL)
        im_vectors.FMH_COSEM(:, iter + 1) = deblur(im_vectors.FMH_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.weighted_mean && options.OSL_OSEM
        im_vectors.Weighted_OSL(:, iter + 1) = deblur(im_vectors.Weighted_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.weighted_mean && options.MBSREM
        im_vectors.Weighted_MBSREM(:, iter + 1) = deblur(im_vectors.Weighted_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.weighted_mean && options.BSREM
        im_vectors.Weighted_BSREM(:,iter+1) = deblur(im_vectors.Weighted_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.weighted_mean && options.ROSEM_MAP
        im_vectors.Weighted_ROSEM(:,iter+1) = deblur(im_vectors.Weighted_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.weighted_mean && options.RBI_OSL
        im_vectors.Weighted_RBI(:, iter + 1) = deblur(im_vectors.Weighted_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.weighted_mean && any(options.COSEM_OSL)
        im_vectors.Weighted_COSEM(:, iter + 1) = deblur(im_vectors.Weighted_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.TV && options.OSL_OSEM
        im_vectors.TV_OSL(:, iter + 1) = deblur(im_vectors.TV_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.TV && options.MBSREM
        im_vectors.TV_MBSREM(:, iter + 1) = deblur(im_vectors.TV_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.TV && options.BSREM
        im_vectors.TV_BSREM(:, iter + 1) = deblur(im_vectors.TV_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.TV && options.ROSEM_MAP
        im_vectors.TV_ROSEM(:, iter + 1) = deblur(im_vectors.TV_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.TV && options.RBI_OSL
        im_vectors.TV_RBI(:, iter + 1) = deblur(im_vectors.TV_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.TV && any(options.COSEM_OSL)
        im_vectors.TV_COSEM(:, iter + 1) = deblur(im_vectors.TV_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.AD && options.OSL_OSEM
        im_vectors.AD_OSL(:, iter + 1) = deblur(im_vectors.AD_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.AD && options.MBSREM
        im_vectors.AD_MBSREM(:, iter + 1) = deblur(im_vectors.AD_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.AD && options.BSREM
        im_vectors.AD_BSREM(:,iter+1) = deblur(im_vectors.AD_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.AD && options.ROSEM_MAP
        im_vectors.AD_ROSEM(:, iter + 1) = deblur(im_vectors.AD_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.AD && options.RBI_OSL
        im_vectors.AD_RBI(:, iter + 1) = deblur(im_vectors.AD_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.AD && any(options.COSEM_OSL)
        im_vectors.AD_COSEM(:, iter + 1) = deblur(im_vectors.AD_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.APLS && options.OSL_OSEM
        im_vectors.APLS_OSL(:, iter + 1) = deblur(im_vectors.APLS_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.APLS && options.MBSREM
        im_vectors.APLS_MBSREM(:, iter + 1) = deblur(im_vectors.APLS_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.APLS && options.BSREM
        im_vectors.APLS_BSREM(:, iter + 1) = deblur(im_vectors.APLS_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.APLS && options.ROSEM_MAP
        im_vectors.APLS_ROSEM(:, iter + 1) = deblur(im_vectors.APLS_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.APLS && options.RBI_OSL
        im_vectors.APLS_RBI(:, iter + 1) = deblur(im_vectors.APLS_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.APLS && any(options.COSEM_OSL)
        im_vectors.APLS_COSEM(:, iter + 1) = deblur(im_vectors.APLS_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.TGV && options.OSL_OSEM
        im_vectors.TGV_OSL(:, iter + 1) = deblur(im_vectors.TGV_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.TGV && options.MBSREM
        im_vectors.TGV_MBSREM(:, iter + 1) = deblur(im_vectors.TGV_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.TGV && options.BSREM
        im_vectors.TGV_BSREM(:, iter + 1) = deblur(im_vectors.TGV_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.TGV && options.ROSEM_MAP
        im_vectors.TGV_ROSEM(:, iter + 1) = deblur(im_vectors.TGV_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.TGV && options.RBI_OSL
        im_vectors.TGV_RBI(:, iter + 1) = deblur(im_vectors.TGV_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.TGV && any(options.COSEM_OSL)
        im_vectors.TGV_COSEM(:, iter + 1) = deblur(im_vectors.TGV_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.NLM && options.OSL_OSEM
        im_vectors.NLM_OSL(:, iter + 1) = deblur(im_vectors.NLM_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.NLM && options.MBSREM
        im_vectors.NLM_MBSREM(:, iter + 1) = deblur(im_vectors.NLM_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.NLM && options.BSREM
        im_vectors.NLM_BSREM(:, iter + 1) = deblur(im_vectors.NLM_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.NLM && options.ROSEM_MAP
        im_vectors.NLM_ROSEM(:, iter + 1) = deblur(im_vectors.NLM_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.NLM && options.RBI_OSL
        im_vectors.NLM_RBI(:, iter + 1) = deblur(im_vectors.NLM_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    elseif options.NLM && any(options.COSEM_OSL)
        im_vectors.NLM_COSEM(:, iter + 1) = deblur(im_vectors.NLM_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
elseif options.implementation == 1
    
    if options.osem
        im_vectors.OSEM(:, iter + 1) = deblur(im_vectors.OSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.mramla
        im_vectors.MRAMLA(:, iter + 1) = deblur(im_vectors.MRAMLA(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.ramla
        im_vectors.RAMLA(:, iter + 1) = deblur(im_vectors.RAMLA(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.rosem
        im_vectors.ROSEM(:, iter + 1) = deblur(im_vectors.ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.rbi
        im_vectors.RBI(:, iter + 1) = deblur(im_vectors.RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.drama
        im_vectors.DRAMA(:, iter + 1) = deblur(im_vectors.DRAMA(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.cosem
        im_vectors.COSEM(:, iter + 1) = deblur(im_vectors.COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.ecosem
        im_vectors.ECOSEM(:, iter + 1) = deblur(im_vectors.ECOSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.acosem
        im_vectors.ACOSEM(:, iter + 1) = deblur(im_vectors.ACOSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.MRP && options.OSL_OSEM
        im_vectors.MRP_OSL(:, iter + 1) = deblur(im_vectors.MRP_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.MRP && options.MBSREM
        im_vectors.MRP_MBSREM(:, iter + 1) = deblur(im_vectors.MRP_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.MRP && options.BSREM
        im_vectors.MRP_BSREM(:,iter+1) = deblur(im_vectors.MRP_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.MRP && options.ROSEM_MAP
        im_vectors.MRP_ROSEM(:,iter+1) = deblur(im_vectors.MRP_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.MRP && options.RBI_OSL
        im_vectors.MRP_RBI(:, iter + 1) = deblur(im_vectors.MRP_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.MRP && any(options.COSEM_OSL)
        im_vectors.MRP_COSEM(:, iter + 1) = deblur(im_vectors.MRP_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.quad && options.OSL_OSEM
        im_vectors.Quad_OSL(:, iter + 1) = deblur(im_vectors.Quad_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.quad && options.MBSREM
        im_vectors.Quad_MBSREM(:, iter + 1) = deblur(im_vectors.Quad_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.quad && options.BSREM
        im_vectors.Quad_BSREM(:,iter+1) = deblur(im_vectors.Quad_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.quad && options.ROSEM_MAP
        im_vectors.Quad_ROSEM(:,iter+1) = deblur(im_vectors.Quad_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.quad && options.RBI_OSL
        im_vectors.Quad_RBI(:, iter + 1) = deblur(im_vectors.Quad_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.quad && any(options.COSEM_OSL)
        im_vectors.Quad_COSEM(:, iter + 1) = deblur(im_vectors.Quad_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.Huber && options.OSL_OSEM
        im_vectors.Huber_OSL(:, iter + 1) = deblur(im_vectors.Huber_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.Huber && options.MBSREM
        im_vectors.Huber_MBSREM(:, iter + 1) = deblur(im_vectors.Huber_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.Huber && options.BSREM
        im_vectors.Huber_BSREM(:,iter+1) = deblur(im_vectors.Huber_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.Huber && options.ROSEM_MAP
        im_vectors.Huber_ROSEM(:,iter+1) = deblur(im_vectors.Huber_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.Huber && options.RBI_OSL
        im_vectors.Huber_RBI(:, iter + 1) = deblur(im_vectors.Huber_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.Huber && any(options.COSEM_OSL)
        im_vectors.Huber_COSEM(:, iter + 1) = deblur(im_vectors.Huber_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.L && options.OSL_OSEM
        im_vectors.L_OSL(:, iter + 1) = deblur(im_vectors.L_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.L && options.MBSREM
        im_vectors.L_MBSREM(:, iter + 1) = deblur(im_vectors.L_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.L && options.BSREM
        im_vectors.L_BSREM(:,iter+1) = deblur(im_vectors.L_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.L && options.ROSEM_MAP
        im_vectors.L_ROSEM(:,iter+1) = deblur(im_vectors.L_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.L && options.RBI_OSL
        im_vectors.L_RBI(:, iter + 1) = deblur(im_vectors.L_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.L && any(options.COSEM_OSL)
        im_vectors.L_COSEM(:, iter + 1) = deblur(im_vectors.L_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.FMH && options.OSL_OSEM
        im_vectors.FMH_OSL(:, iter + 1) = deblur(im_vectors.FMH_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.FMH && options.MBSREM
        im_vectors.FMH_MBSREM(:, iter + 1) = deblur(im_vectors.FMH_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.FMH && options.BSREM
        im_vectors.FMH_BSREM(:,iter+1) = deblur(im_vectors.FMH_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.FMH && options.ROSEM_MAP
        im_vectors.FMH_ROSEM(:,iter+1) = deblur(im_vectors.FMH_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.FMH && options.RBI_OSL
        im_vectors.FMH_RBI(:, iter + 1) = deblur(im_vectors.FMH_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.FMH && any(options.COSEM_OSL)
        im_vectors.FMH_COSEM(:, iter + 1) = deblur(im_vectors.FMH_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.weighted_mean && options.OSL_OSEM
        im_vectors.Weighted_OSL(:, iter + 1) = deblur(im_vectors.Weighted_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.weighted_mean && options.MBSREM
        im_vectors.Weighted_MBSREM(:, iter + 1) = deblur(im_vectors.Weighted_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.weighted_mean && options.BSREM
        im_vectors.Weighted_BSREM(:,iter+1) = deblur(im_vectors.Weighted_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.weighted_mean && options.ROSEM_MAP
        im_vectors.Weighted_ROSEM(:,iter+1) = deblur(im_vectors.Weighted_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.weighted_mean && options.RBI_OSL
        im_vectors.Weighted_RBI(:, iter + 1) = deblur(im_vectors.Weighted_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.weighted_mean && any(options.COSEM_OSL)
        im_vectors.Weighted_COSEM(:, iter + 1) = deblur(im_vectors.Weighted_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.TV && options.OSL_OSEM
        im_vectors.TV_OSL(:, iter + 1) = deblur(im_vectors.TV_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.TV && options.MBSREM
        im_vectors.TV_MBSREM(:, iter + 1) = deblur(im_vectors.TV_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.TV && options.BSREM
        im_vectors.TV_BSREM(:, iter + 1) = deblur(im_vectors.TV_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.TV && options.ROSEM_MAP
        im_vectors.TV_ROSEM(:, iter + 1) = deblur(im_vectors.TV_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.TV && options.RBI_OSL
        im_vectors.TV_RBI(:, iter + 1) = deblur(im_vectors.TV_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.TV && any(options.COSEM_OSL)
        im_vectors.TV_COSEM(:, iter + 1) = deblur(im_vectors.TV_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.AD && options.OSL_OSEM
        im_vectors.AD_OSL(:, iter + 1) = deblur(im_vectors.AD_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.AD && options.MBSREM
        im_vectors.AD_MBSREM(:, iter + 1) = deblur(im_vectors.AD_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.AD && options.BSREM
        im_vectors.AD_BSREM(:,iter+1) = deblur(im_vectors.AD_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.AD && options.ROSEM_MAP
        im_vectors.AD_ROSEM(:, iter + 1) = deblur(im_vectors.AD_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.AD && options.RBI_OSL
        im_vectors.AD_RBI(:, iter + 1) = deblur(im_vectors.AD_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.AD && any(options.COSEM_OSL)
        im_vectors.AD_COSEM(:, iter + 1) = deblur(im_vectors.AD_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.APLS && options.OSL_OSEM
        im_vectors.APLS_OSL(:, iter + 1) = deblur(im_vectors.APLS_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.APLS && options.MBSREM
        im_vectors.APLS_MBSREM(:, iter + 1) = deblur(im_vectors.APLS_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.APLS && options.BSREM
        im_vectors.APLS_BSREM(:, iter + 1) = deblur(im_vectors.APLS_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.APLS && options.ROSEM_MAP
        im_vectors.APLS_ROSEM(:, iter + 1) = deblur(im_vectors.APLS_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.APLS && options.RBI_OSL
        im_vectors.APLS_RBI(:, iter + 1) = deblur(im_vectors.APLS_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.APLS && any(options.COSEM_OSL)
        im_vectors.APLS_COSEM(:, iter + 1) = deblur(im_vectors.APLS_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.TGV && options.OSL_OSEM
        im_vectors.TGV_OSL(:, iter + 1) = deblur(im_vectors.TGV_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.TGV && options.MBSREM
        im_vectors.TGV_MBSREM(:, iter + 1) = deblur(im_vectors.TGV_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.TGV && options.BSREM
        im_vectors.TGV_BSREM(:, iter + 1) = deblur(im_vectors.TGV_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.TGV && options.ROSEM_MAP
        im_vectors.TGV_ROSEM(:, iter + 1) = deblur(im_vectors.TGV_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.TGV && options.RBI_OSL
        im_vectors.TGV_RBI(:, iter + 1) = deblur(im_vectors.TGV_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.TGV && any(options.COSEM_OSL)
        im_vectors.TGV_COSEM(:, iter + 1) = deblur(im_vectors.TGV_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.NLM && options.OSL_OSEM
        im_vectors.NLM_OSL(:, iter + 1) = deblur(im_vectors.NLM_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.NLM && options.MBSREM
        im_vectors.NLM_MBSREM(:, iter + 1) = deblur(im_vectors.NLM_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.NLM && options.BSREM
        im_vectors.NLM_BSREM(:, iter + 1) = deblur(im_vectors.NLM_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.NLM && options.ROSEM_MAP
        im_vectors.NLM_ROSEM(:, iter + 1) = deblur(im_vectors.NLM_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.NLM && options.RBI_OSL
        im_vectors.NLM_RBI(:, iter + 1) = deblur(im_vectors.NLM_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.NLM && any(options.COSEM_OSL)
        im_vectors.NLM_COSEM(:, iter + 1) = deblur(im_vectors.NLM_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.custom && options.OSL_OSEM
        im_vectors.custom_OSL(:, iter + 1) = deblur(im_vectors.custom_OSL(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.custom && options.MBSREM
        im_vectors.custom_MBSREM(:, iter + 1) = deblur(im_vectors.custom_MBSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.custom && options.BSREM
        im_vectors.custom_BSREM(:, iter + 1) = deblur(im_vectors.custom_BSREM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.custom && options.ROSEM_MAP
        im_vectors.custom_ROSEM(:, iter + 1) = deblur(im_vectors.custom_ROSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.custom && options.RBI_OSL
        im_vectors.custom_RBI(:, iter + 1) = deblur(im_vectors.custom_RBI(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
    if options.custom && any(options.COSEM_OSL)
        im_vectors.custom_COSEM(:, iter + 1) = deblur(im_vectors.custom_COSEM(:, iter + 1), options, iter, subsets, gaussK, Nx, Ny, Nz);
    end
end