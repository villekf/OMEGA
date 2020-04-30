function im_vectors = transfer_im_vectors(im_vectors, pz, options, iter)

oo = 1;
if options.osem
    im_vectors.OSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
%     im_vectors.OSEM_apu = im_vectors.OSEM_apu;
end
oo = oo + 1;
if options.mlem
    im_vectors.MLEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.mramla
    im_vectors.MRAMLA_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.ramla
    im_vectors.RAMLA_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.rosem
    im_vectors.ROSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.rbi
    im_vectors.RBI_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.drama
    im_vectors.DRAMA_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.cosem
    im_vectors.COSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.ecosem
    im_vectors.ECOSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
    if ~options.osem
        im_vectors.OSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
    end
    if ~options.cosem
        im_vectors.COSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
    end
end
oo = oo + 1;
if options.acosem
    im_vectors.ACOSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;

if options.MRP && options.OSL_OSEM
    im_vectors.MRP_OSL_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.MRP && options.OSL_MLEM
    im_vectors.MRP_MLEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.MRP && options.MBSREM
    im_vectors.MRP_MBSREM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.MRP && options.BSREM
    im_vectors.MRP_BSREM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.MRP && options.ROSEM_MAP
    im_vectors.MRP_ROSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.MRP && options.RBI_OSL
    im_vectors.MRP_RBI_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.MRP && any(options.COSEM_OSL)
    im_vectors.MRP_COSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;

if options.quad && options.OSL_OSEM
    im_vectors.Quad_OSL_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.quad && options.OSL_MLEM
    im_vectors.Quad_MLEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.quad && options.MBSREM
    im_vectors.Quad_MBSREM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.quad && options.BSREM
    im_vectors.Quad_BSREM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.quad && options.ROSEM_MAP
    im_vectors.Quad_ROSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.quad && options.RBI_OSL
    im_vectors.Quad_RBI_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.quad && any(options.COSEM_OSL)
    im_vectors.Quad_COSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;

if options.L && options.OSL_OSEM
    im_vectors.L_OSL_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.L && options.OSL_MLEM
    im_vectors.L_MLEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.L && options.MBSREM
    im_vectors.L_MBSREM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.L && options.BSREM
    im_vectors.L_BSREM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.L && options.ROSEM_MAP
    im_vectors.L_ROSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.L && options.RBI_OSL
    im_vectors.L_RBI_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.L && any(options.COSEM_OSL)
    im_vectors.L_COSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;

if options.FMH && options.OSL_OSEM
    im_vectors.FMH_OSL_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.FMH && options.OSL_MLEM
    im_vectors.FMH_MLEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.FMH && options.MBSREM
    im_vectors.FMH_MBSREM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.FMH && options.BSREM
    im_vectors.FMH_BSREM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.FMH && options.ROSEM_MAP
    im_vectors.FMH_ROSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.FMH && options.RBI_OSL
    im_vectors.FMH_RBI_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.FMH && any(options.COSEM_OSL)
    im_vectors.FMH_COSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;

if options.weighted_mean && options.OSL_OSEM
    im_vectors.Weighted_OSL_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.weighted_mean && options.OSL_MLEM
    im_vectors.Weighted_MLEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.weighted_mean && options.MBSREM
    im_vectors.Weighted_MBSREM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.weighted_mean && options.BSREM
    im_vectors.Weighted_BSREM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.weighted_mean && options.ROSEM_MAP
    im_vectors.Weighted_ROSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.weighted_mean && options.RBI_OSL
    im_vectors.Weighted_RBI_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.weighted_mean && any(options.COSEM_OSL)
    im_vectors.Weighted_COSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;

if options.TV && options.OSL_OSEM
    im_vectors.TV_OSL_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.TV && options.OSL_MLEM
    im_vectors.TV_MLEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.TV && options.MBSREM
    im_vectors.TV_MBSREM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.TV && options.BSREM
    im_vectors.TV_BSREM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.TV && options.ROSEM_MAP
    im_vectors.TV_ROSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.TV && options.RBI_OSL
    im_vectors.TV_RBI_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.TV && any(options.COSEM_OSL)
    im_vectors.TV_COSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;

if options.AD && options.OSL_OSEM
    im_vectors.AD_OSL_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.AD && options.OSL_MLEM
    im_vectors.AD_MLEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.AD && options.MBSREM
    im_vectors.AD_MBSREM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.AD && options.BSREM
    im_vectors.AD_BSREM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.AD && options.ROSEM_MAP
    im_vectors.AD_ROSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.AD && options.RBI_OSL
    im_vectors.AD_RBI_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.AD && any(options.COSEM_OSL)
    im_vectors.AD_COSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;

if options.APLS && options.OSL_OSEM
    im_vectors.APLS_OSL_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.APLS && options.OSL_MLEM
    im_vectors.APLS_MLEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.APLS && options.MBSREM
    im_vectors.APLS_MBSREM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.APLS && options.BSREM
    im_vectors.APLS_BSREM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.APLS && options.ROSEM_MAP
    im_vectors.APLS_ROSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.APLS && options.RBI_OSL
    im_vectors.APLS_RBI_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.APLS && any(options.COSEM_OSL)
    im_vectors.APLS_COSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;

if options.TGV && options.OSL_OSEM
    im_vectors.TGV_OSL_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.TGV && options.OSL_MLEM
    im_vectors.TGV_MLEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.TGV && options.MBSREM
    im_vectors.TGV_MBSREM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.TGV && options.BSREM
    im_vectors.TGV_BSREM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.TGV && options.ROSEM_MAP
    im_vectors.TGV_ROSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.TGV && options.RBI_OSL
    im_vectors.TGV_RBI_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.TGV && any(options.COSEM_OSL)
    im_vectors.TGV_COSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;

if options.NLM && options.OSL_OSEM
    im_vectors.NLM_OSL_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.NLM && options.OSL_MLEM
    im_vectors.NLM_MLEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.NLM && options.MBSREM
    im_vectors.NLM_MBSREM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.NLM && options.BSREM
    im_vectors.NLM_BSREM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.NLM && options.ROSEM_MAP
    im_vectors.NLM_ROSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.NLM && options.RBI_OSL
    im_vectors.NLM_RBI_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.NLM && any(options.COSEM_OSL)
    im_vectors.NLM_COSEM_apu = reshape(pz{oo}(:,:,:,iter + 1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;


if options.OSL_OSEM && options.custom
    im_vectors.custom_OSL_apu = reshape(pz{oo}(:,:,:,1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.OSL_MLEM && options.custom
    im_vectors.custom_MLEM_apu = reshape(pz{oo}(:,:,:,1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.MBSREM && options.custom
    im_vectors.custom_MBSREM_apu = reshape(pz{oo}(:,:,:,1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.BSREM && options.custom
    im_vectors.custom_BSREM_apu = reshape(pz{oo}(:,:,:,1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.ROSEM_MAP && options.custom
    im_vectors.custom_ROSEM_apu = reshape(pz{oo}(:,:,:,1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if options.RBI_OSL && options.custom
    im_vectors.custom_RBI_apu = reshape(pz{oo}(:,:,:,1), options.Nx*options.Ny*options.Nz, 1);
end
oo = oo + 1;
if any(options.COSEM_OSL) && options.custom
    im_vectors.custom_COSEM_apu = reshape(pz{oo}(:,:,:,1), options.Nx*options.Ny*options.Nz, 1);
end