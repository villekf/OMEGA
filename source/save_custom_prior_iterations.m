function [options, pz] = save_custom_prior_iterations(options, llo, pz)
% Save the current time step and all the iterations in a cell

options.im_vectors = reshape_vectors(options.im_vectors, options);

[pz, gg] = images_to_cell(options.im_vectors, llo, pz, options);
if options.save_iter
    Niter = options.Niter + 1;
else
    Niter = 1;
end
gg = gg + 1;
if options.OSL_OSEM && options.custom
    pz{gg, llo} = reshape(options.im_vectors.custom_OSL,options.Nx,options.Ny,options.Nz,Niter);
end
gg = gg + 1;
if options.implementation == 2
    if options.OSL_MLEM && options.custom
        pz{gg, llo} = reshape(options.im_vectors.custom_MLEM,options.Nx,options.Ny,options.Nz,Niter);
    end
end
gg = gg + 1;
if options.BSREM && options.custom
    pz{gg, llo} = reshape(options.im_vectors.custom_BSREM,options.Nx,options.Ny,options.Nz,Niter);
end
gg = gg + 1;
if options.MBSREM && options.custom
    pz{gg, llo} = reshape(options.im_vectors.custom_MBSREM,options.Nx,options.Ny,options.Nz,Niter);
end
gg = gg + 1;
if options.ROSEM_MAP && options.custom
    pz{gg, llo} = reshape(options.im_vectors.custom_ROSEM,options.Nx,options.Ny,options.Nz,Niter);
end
gg = gg + 1;
if options.OSL_RBI && options.custom
    pz{gg, llo} = reshape(options.im_vectors.custom_RBI,options.Nx,options.Ny,options.Nz,Niter);
end
gg = gg + 1;
if any(options.OSL_COSEM) && options.custom
    pz{gg, llo} = reshape(options.im_vectors.custom_COSEM,options.Nx,options.Ny,options.Nz,Niter);
end
pz{end,1} = 0;

if options.partitions > 1 && llo < options.partitions
    
    if options.OSEM
        options.im_vectors.OSEM(:, 1) = options.im_vectors.OSEM(:,end);
    end
    
    if options.MLEM
        options.im_vectors.MLEM(:, 1) = options.im_vectors.MLEM(:,end);
    end
    
    if options.MRAMLA
        options.im_vectors.MRAMLA(:, 1) = options.im_vectors.MRAMLA(:,end);
    end
    
    if options.RAMLA
        options.im_vectors.RAMLA(:, 1) = options.im_vectors.RAMLA(:,end);
    end
    
    if options.ROSEM
        options.im_vectors.ROSEM(:, 1) = options.im_vectors.ROSEM(:,end);
    end
    
    if options.RBI
        options.im_vectors.RBI(:, 1) = options.im_vectors.RBI(:,end);
    end
    
    if options.DRAMA
        options.im_vectors.DRAMA(:, 1) = options.im_vectors.DRAMA(:,end);
    end
    
    if options.COSEM
        options.im_vectors.COSEM(:, 1) = options.im_vectors.COSEM(:,end);
    end
    
    if options.ECOSEM
        options.im_vectors.ECOSEM(:, 1) = options.im_vectors.ECOSEM(:,end);
    end
    
    if options.ACOSEM
        options.im_vectors.ACOSEM(:, 1) = options.im_vectors.ACOSEM(:,end);
    end
    
    if options.MRP && options.OSL_OSEM
        options.im_vectors.MRP_OSL(:, 1) = options.im_vectors.MRP_OSL(:,end);
    end
    if options.MRP && options.OSL_MLEM
        options.im_vectors.MRP_MLEM(:, 1) = options.im_vectors.MRP_MLEM(:,end);
    end
    if options.MRP && options.MBSREM
        options.im_vectors.MRP_MBSREM(:, 1) = options.im_vectors.MRP_MBSREM(:,end);
    end
    
    if options.MRP && options.BSREM
        options.im_vectors.MRP_BSREM(:, 1) = options.im_vectors.MRP_BSREM(:, end);
    end
    
    if options.MRP && options.ROSEM_MAP
        options.im_vectors.MRP_ROSEM(:, 1) = options.im_vectors.MRP_ROSEM(:, end);
    end
    
    if options.MRP && options.OSL_RBI
        options.im_vectors.MRP_RBI(:, 1) = options.im_vectors.MRP_RBI(:,end);
    end
    
    if options.MRP && any(options.OSL_COSEM)
        options.im_vectors.MRP_COSEM(:, 1) = options.im_vectors.MRP_COSEM(:,end);
    end
    
    if options.quad && options.OSL_OSEM
        options.im_vectors.Quad_OSL(:, 1) = options.im_vectors.Quad_OSL(:,end);
    end
    if options.quad && options.OSL_MLEM
        options.im_vectors.Quad_MLEM(:, 1) = options.im_vectors.Quad_MLEM(:,end);
    end
    if options.quad && options.MBSREM
        options.im_vectors.Quad_MBSREM(:, 1) = options.im_vectors.Quad_MBSREM(:,end);
    end
    
    if options.quad && options.BSREM
        options.im_vectors.Quad_BSREM(:, 1) = options.im_vectors.Quad_BSREM(:, end);
    end
    
    if options.quad && options.ROSEM_MAP
        options.im_vectors.Quad_ROSEM(:, 1) = options.im_vectors.Quad_ROSEM(:, end);
    end
    
    if options.quad && options.OSL_RBI
        options.im_vectors.Quad_RBI(:, 1) = options.im_vectors.Quad_RBI(:,end);
    end
    
    if options.quad && any(options.OSL_COSEM)
        options.im_vectors.Quad_COSEM(:, 1) = options.im_vectors.Quad_COSEM(:,end);
    end
    
    if options.L && options.OSL_OSEM
        options.im_vectors.L_OSL(:, 1) = options.im_vectors.L_OSL(:,end);
    end
    if options.L && options.OSL_MLEM
        options.im_vectors.L_MLEM(:, 1) = options.im_vectors.L_MLEM(:,end);
    end
    if options.L && options.MBSREM
        options.im_vectors.L_MBSREM(:, 1) = options.im_vectors.L_MBSREM(:,end);
    end
    
    if options.L && options.BSREM
        options.im_vectors.L_BSREM(:, 1) = options.im_vectors.L_BSREM(:, end);
    end
    
    if options.L && options.ROSEM_MAP
        options.im_vectors.L_ROSEM(:, 1) = options.im_vectors.L_ROSEM(:, end);
    end
    
    if options.L && options.OSL_RBI
        options.im_vectors.L_RBI(:, 1) = options.im_vectors.L_RBI(:,end);
    end
    
    if options.L && any(options.OSL_COSEM)
        options.im_vectors.L_COSEM(:, 1) = options.im_vectors.L_COSEM(:,end);
    end
    
    if options.FMH && options.OSL_OSEM
        options.im_vectors.FMH_OSL(:, 1) = options.im_vectors.FMH_OSL(:,end);
    end
    if options.FMH && options.OSL_MLEM
        options.im_vectors.FMH_MLEM(:, 1) = options.im_vectors.FMH_MLEM(:,end);
    end
    if options.FMH && options.MBSREM
        options.im_vectors.FMH_MBSREM(:, 1) = options.im_vectors.FMH_MBSREM(:,end);
    end
    
    if options.FMH && options.BSREM
        options.im_vectors.FMH_BSREM(:, 1) = options.im_vectors.FMH_BSREM(:, end);
    end
    
    if options.FMH && options.ROSEM_MAP
        options.im_vectors.FMH_ROSEM(:, 1) = options.im_vectors.FMH_ROSEM(:, end);
    end
    
    if options.FMH && options.OSL_RBI
        options.im_vectors.FMH_RBI(:, 1) = options.im_vectors.FMH_RBI(:,end);
    end
    
    if options.FMH && any(options.OSL_COSEM)
        options.im_vectors.FMH_COSEM(:, 1) = options.im_vectors.FMH_COSEM(:,end);
    end
    
    if options.weighted_mean && options.OSL_OSEM
        options.im_vectors.Weighted_OSL(:, 1) = options.im_vectors.Weighted_OSL(:,end);
    end
    if options.weighted_mean && options.OSL_MLEM
        options.im_vectors.Weighted_MLEM(:, 1) = options.im_vectors.Weighted_MLEM(:,end);
    end
    if options.weighted_mean && options.MBSREM
        options.im_vectors.Weighted_MBSREM(:, 1) = options.im_vectors.Weighted_MBSREM(:,end);
    end
    
    if options.weighted_mean && options.BSREM
        options.im_vectors.Weighted_BSREM(:, 1) = options.im_vectors.Weighted_BSREM(:, end);
    end
    
    if options.weighted_mean && options.ROSEM_MAP
        options.im_vectors.Weighted_ROSEM(:, 1) = options.im_vectors.Weighted_ROSEM(:, end);
    end
    
    if options.weighted_mean && options.OSL_RBI
        options.im_vectors.Weighted_RBI(:, 1) = options.im_vectors.Weighted_RBI(:,end);
    end
    
    if options.weighted_mean && any(options.OSL_COSEM)
        options.im_vectors.Weighted_COSEM(:, 1) = options.im_vectors.Weighted_COSEM(:,end);
    end
    
    if options.TV && options.OSL_OSEM
        options.im_vectors.TV_OSL(:, 1) = options.im_vectors.TV_OSL(:,end);
    end
    if options.TV && options.OSL_MLEM
        options.im_vectors.TV_MLEM(:, 1) = options.im_vectors.TV_MLEM(:,end);
    end
    if options.TV && options.MBSREM
        options.im_vectors.TV_MBSREM(:, 1) = options.im_vectors.TV_MBSREM(:,end);
    end
    
    if options.TV && options.BSREM
        options.im_vectors.TV_BSREM(:, 1) = options.im_vectors.TV_BSREM(:, end);
    end
    
    if options.TV && options.ROSEM_MAP
        options.im_vectors.TV_ROSEM(:, 1) = options.im_vectors.TV_ROSEM(:, end);
    end
    
    if options.TV && options.OSL_RBI
        options.im_vectors.TV_RBI(:, 1) = options.im_vectors.TV_RBI(:,end);
    end
    
    if options.TV && any(options.OSL_COSEM)
        options.im_vectors.TV_COSEM(:, 1) = options.im_vectors.TV_COSEM(:,end);
    end
    
    if options.AD && options.OSL_OSEM
        options.im_vectors.AD_OSL(:, 1) = options.im_vectors.AD_OSL(:,end);
    end
    if options.AD && options.OSL_MLEM
        options.im_vectors.AD_MLEM(:, 1) = options.im_vectors.AD_MLEM(:,end);
    end
    if options.AD && options.MBSREM
        options.im_vectors.AD_MBSREM(:, 1) = options.im_vectors.AD_MBSREM(:,end);
    end
    
    if options.AD && options.BSREM
        options.im_vectors.AD_BSREM(:, 1) = options.im_vectors.AD_BSREM(:, end);
    end
    if options.AD && options.ROSEM_MAP
        options.im_vectors.AD_ROSEM(:, 1) = options.im_vectors.AD_ROSEM(:, end);
    end
    if options.AD && options.OSL_RBI
        options.im_vectors.AD_RBI(:, 1) = options.im_vectors.AD_RBI(:,end);
    end
    if options.AD && any(options.OSL_COSEM)
        options.im_vectors.AD_COSEM(:, 1) = options.im_vectors.AD_COSEM(:,end);
    end
    
    if options.APLS && options.OSL_OSEM
        options.im_vectors.APLS_OSL(:, 1) = options.im_vectors.APLS_OSL(:,end);
    end
    if options.APLS && options.OSL_MLEM
        options.im_vectors.APLS_MLEM(:, 1) = options.im_vectors.APLS_MLEM(:,end);
    end
    if options.APLS && options.MBSREM
        options.im_vectors.APLS_MBSREM(:, 1) = options.im_vectors.APLS_MBSREM(:,end);
    end
    if options.APLS && options.BSREM
        options.im_vectors.APLS_BSREM(:, 1) = options.im_vectors.APLS_BSREM(:, end);
    end
    if options.APLS && options.ROSEM_MAP
        options.im_vectors.APLS_ROSEM(:, 1) = options.im_vectors.APLS_ROSEM(:, end);
    end
    if options.APLS && options.OSL_RBI
        options.im_vectors.APLS_RBI(:, 1) = options.im_vectors.APLS_RBI(:,end);
    end
    if options.APLS && any(options.OSL_COSEM)
        options.im_vectors.APLS_COSEM(:, 1) = options.im_vectors.APLS_COSEM(:,end);
    end
    
    if options.TGV && options.OSL_OSEM
        options.im_vectors.TGV_OSL(:, 1) = options.im_vectors.TGV_OSL(:,end);
    end
    if options.TGV && options.OSL_MLEM
        options.im_vectors.TGV_MLEM(:, 1) = options.im_vectors.TGV_MLEM(:,end);
    end
    if options.TGV && options.MBSREM
        options.im_vectors.TGV_MBSREM(:, 1) = options.im_vectors.TGV_MBSREM(:,end);
    end
    if options.TGV && options.BSREM
        options.im_vectors.TGV_BSREM(:, 1) = options.im_vectors.TGV_BSREM(:, end);
    end
    if options.TGV && options.ROSEM_MAP
        options.im_vectors.TGV_ROSEM(:, 1) = options.im_vectors.TGV_ROSEM(:, end);
    end
    if options.TGV && options.OSL_RBI
        options.im_vectors.TGV_RBI(:, 1) = options.im_vectors.TGV_RBI(:,end);
    end
    if options.TGV && any(options.OSL_COSEM)
        options.im_vectors.TGV_COSEM(:, 1) = options.im_vectors.TGV_COSEM(:,end);
    end
    
    if options.NLM && options.OSL_OSEM
        options.im_vectors.NLM_OSL(:, 1) = options.im_vectors.NLM_OSL(:,end);
    end
    if options.NLM && options.MBSREM
        options.im_vectors.NLM_MBSREM(:, 1) = options.im_vectors.NLM_MBSREM(:,end);
    end
    if options.NLM && options.BSREM
        options.im_vectors.NLM_BSREM(:, 1) = options.im_vectors.NLM_BSREM(:, end);
    end
    if options.NLM && options.ROSEM_MAP
        options.im_vectors.NLM_ROSEM(:, 1) = options.im_vectors.NLM_ROSEM(:, end);
    end
    if options.NLM && options.OSL_RBI
        options.im_vectors.NLM_RBI(:, 1) = options.im_vectors.NLM_RBI(:,end);
    end
    if options.NLM && any(options.OSL_COSEM)
        options.im_vectors.NLM_COSEM(:, 1) = options.im_vectors.NLM_COSEM(:,end);
    end
    
    if options.OSL_OSEM && options.custom
        options.im_vectors.custom_OSL(:, 1) = options.im_vectors.custom_OSL(:,end);
    end
    if options.OSL_MLEM && options.custom
        options.im_vectors.custom_MLEM(:, 1) = options.im_vectors.custom_MLEM(:,end);
    end
    if options.MBSREM && options.custom
        options.im_vectors.custom_MBSREM(:, 1) = options.im_vectors.custom_MBSREM(:,end);
    end
    if options.BSREM && options.custom
        options.im_vectors.custom_BSREM(:, 1) = options.im_vectors.custom_BSREM(:, end);
    end
    if options.ROSEM_MAP && options.custom
        options.im_vectors.custom_ROSEM(:, 1) = options.im_vectors.custom_ROSEM(:, end);
    end
    if options.OSL_RBI && options.custom
        options.im_vectors.custom_RBI(:, 1) = options.im_vectors.custom_RBI(:,end);
    end
    if any(options.OSL_COSEM) && options.custom
        options.im_vectors.custom_COSEM(:, 1) = options.im_vectors.custom_COSEM(:,end);
    end
end