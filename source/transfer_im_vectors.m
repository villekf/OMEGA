function im_vectors = transfer_im_vectors(im_vectors, pz, iter, options)

oo = 1;
if options.osem
    im_vectors.OSEM(:,iter) = pz{oo};
%     im_vectors.OSEM_apu = im_vectors.OSEM(:,iter);
end
oo = oo + 1;
if options.mlem
    im_vectors.MLEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.mramla
    im_vectors.MRAMLA(:,iter) = pz{oo};
end
oo = oo + 1;
if options.ramla
    im_vectors.RAMLA(:,iter) = pz{oo};
end
oo = oo + 1;
if options.rosem
    im_vectors.ROSEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.rbi
    im_vectors.RBI(:,iter) = pz{oo};
end
oo = oo + 1;
if options.drama
    im_vectors.DRAMA(:,iter) = pz{oo};
end
oo = oo + 1;
if options.cosem
    im_vectors.COSEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.ecosem
    im_vectors.ECOSEM(:,iter) = pz{oo};
    if ~options.osem
        im_vectors.OSEM(:,iter) = pz{oo};
    end
    if ~options.cosem
        im_vectors.COSEM(:,iter) = pz{oo};
    end
end
oo = oo + 1;
if options.acosem
    im_vectors.ACOSEM(:,iter) = pz{oo};
end
oo = oo + 1;

if options.MRP && options.OSL_OSEM
    im_vectors.MRP_OSL(:,iter) = pz{oo};
end
oo = oo + 1;
if options.MRP && options.OSL_MLEM
    im_vectors.MRP_MLEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.MRP && options.MBSREM
    im_vectors.MRP_MBSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.MRP && options.BSREM
    im_vectors.MRP_BSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.MRP && options.ROSEM_MAP
    im_vectors.MRP_ROSEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.MRP && options.RBI_MAP
    im_vectors.MRP_RBI(:,iter) = pz{oo};
end
oo = oo + 1;
if options.MRP && any(options.COSEM_MAP)
    im_vectors.MRP_COSEM(:,iter) = pz{oo};
end
oo = oo + 1;

if options.quad && options.OSL_OSEM
    im_vectors.Quad_OSL(:,iter) = pz{oo};
end
oo = oo + 1;
if options.quad && options.OSL_MLEM
    im_vectors.Quad_MLEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.quad && options.MBSREM
    im_vectors.Quad_MBSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.quad && options.BSREM
    im_vectors.Quad_BSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.quad && options.ROSEM_MAP
    im_vectors.Quad_ROSEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.quad && options.RBI_MAP
    im_vectors.Quad_RBI(:,iter) = pz{oo};
end
oo = oo + 1;
if options.quad && any(options.COSEM_MAP)
    im_vectors.Quad_COSEM(:,iter) = pz{oo};
end
oo = oo + 1;

if options.L && options.OSL_OSEM
    im_vectors.L_OSL(:,iter) = pz{oo};
end
oo = oo + 1;
if options.L && options.OSL_MLEM
    im_vectors.L_MLEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.L && options.MBSREM
    im_vectors.L_MBSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.L && options.BSREM
    im_vectors.L_BSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.L && options.ROSEM_MAP
    im_vectors.L_ROSEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.L && options.RBI_MAP
    im_vectors.L_RBI(:,iter) = pz{oo};
end
oo = oo + 1;
if options.L && any(options.COSEM_MAP)
    im_vectors.L_COSEM(:,iter) = pz{oo};
end
oo = oo + 1;

if options.FMH && options.OSL_OSEM
    im_vectors.FMH_OSL(:,iter) = pz{oo};
end
oo = oo + 1;
if options.FMH && options.OSL_MLEM
    im_vectors.FMH_MLEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.FMH && options.MBSREM
    im_vectors.FMH_MBSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.FMH && options.BSREM
    im_vectors.FMH_BSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.FMH && options.ROSEM_MAP
    im_vectors.FMH_ROSEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.FMH && options.RBI_MAP
    im_vectors.FMH_RBI(:,iter) = pz{oo};
end
oo = oo + 1;
if options.FMH && any(options.COSEM_MAP)
    im_vectors.FMH_COSEM(:,iter) = pz{oo};
end
oo = oo + 1;

if options.weighted_mean && options.OSL_OSEM
    im_vectors.Weighted_OSL(:,iter) = pz{oo};
end
oo = oo + 1;
if options.weighted_mean && options.OSL_MLEM
    im_vectors.Weighted_MLEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.weighted_mean && options.MBSREM
    im_vectors.Weighted_MBSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.weighted_mean && options.BSREM
    im_vectors.Weighted_BSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.weighted_mean && options.ROSEM_MAP
    im_vectors.Weighted_ROSEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.weighted_mean && options.RBI_MAP
    im_vectors.Weighted_RBI(:,iter) = pz{oo};
end
oo = oo + 1;
if options.weighted_mean && any(options.COSEM_MAP)
    im_vectors.Weighted_COSEM(:,iter) = pz{oo};
end
oo = oo + 1;

if options.TV && options.OSL_OSEM
    im_vectors.TV_OSL(:,iter) = pz{oo};
end
oo = oo + 1;
if options.TV && options.OSL_MLEM
    im_vectors.TV_MLEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.TV && options.MBSREM
    im_vectors.TV_MBSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.TV && options.BSREM
    im_vectors.TV_BSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.TV && options.ROSEM_MAP
    im_vectors.TV_ROSEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.TV && options.RBI_MAP
    im_vectors.TV_RBI(:,iter) = pz{oo};
end
oo = oo + 1;
if options.TV && any(options.COSEM_MAP)
    im_vectors.TV_COSEM(:,iter) = pz{oo};
end
oo = oo + 1;

if options.AD && options.OSL_OSEM
    im_vectors.AD_OSL(:,iter) = pz{oo};
end
oo = oo + 1;
if options.AD && options.OSL_MLEM
    im_vectors.AD_MLEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.AD && options.MBSREM
    im_vectors.AD_MBSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.AD && options.BSREM
    im_vectors.AD_BSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.AD && options.ROSEM_MAP
    im_vectors.AD_ROSEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.AD && options.RBI_MAP
    im_vectors.AD_RBI(:,iter) = pz{oo};
end
oo = oo + 1;
if options.AD && any(options.COSEM_MAP)
    im_vectors.AD_COSEM(:,iter) = pz{oo};
end
oo = oo + 1;

if options.APLS && options.OSL_OSEM
    im_vectors.APLS_OSL(:,iter) = pz{oo};
end
oo = oo + 1;
if options.APLS && options.OSL_MLEM
    im_vectors.APLS_MLEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.APLS && options.MBSREM
    im_vectors.APLS_MBSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.APLS && options.BSREM
    im_vectors.APLS_BSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.APLS && options.ROSEM_MAP
    im_vectors.APLS_ROSEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.APLS && options.RBI_MAP
    im_vectors.APLS_RBI(:,iter) = pz{oo};
end
oo = oo + 1;
if options.APLS && any(options.COSEM_MAP)
    im_vectors.APLS_COSEM(:,iter) = pz{oo};
end
oo = oo + 1;

if options.TGV && options.OSL_OSEM
    im_vectors.TGV_OSL(:,iter) = pz{oo};
end
oo = oo + 1;
if options.TGV && options.OSL_MLEM
    im_vectors.TGV_MLEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.TGV && options.MBSREM
    im_vectors.TGV_MBSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.TGV && options.BSREM
    im_vectors.TGV_BSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.TGV && options.ROSEM_MAP
    im_vectors.TGV_ROSEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.TGV && options.RBI_MAP
    im_vectors.TGV_RBI(:,iter) = pz{oo};
end
oo = oo + 1;
if options.TGV && any(options.COSEM_MAP)
    im_vectors.TGV_COSEM(:,iter) = pz{oo};
end
oo = oo + 1;

if options.NLM && options.OSL_OSEM
    im_vectors.NLM_OSL(:,iter) = pz{oo};
end
oo = oo + 1;
if options.NLM && options.OSL_MLEM
    im_vectors.NLM_MLEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.NLM && options.MBSREM
    im_vectors.NLM_MBSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.NLM && options.BSREM
    im_vectors.NLM_BSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.NLM && options.ROSEM_MAP
    im_vectors.NLM_ROSEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.NLM && options.RBI_MAP
    im_vectors.NLM_RBI(:,iter) = pz{oo};
end
oo = oo + 1;
if options.NLM && any(options.COSEM_MAP)
    im_vectors.NLM_COSEM(:,iter) = pz{oo};
end
oo = oo + 1;


if options.OSL_MLEM
    im_vectors.custom_MLEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.OSL_OSEM
    im_vectors.custom_OSEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.MBSREM
    im_vectors.custom_MBSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.BSREM
    im_vectors.custom_BSREM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.ROSEM_MAP
    im_vectors.custom_ROSEM(:,iter) = pz{oo};
end
oo = oo + 1;
if options.RBI_MAP
    im_vectors.custom_RBI(:,iter) = pz{oo};
end
oo = oo + 1;
if any(options.COSEM_MAP)
    im_vectors.custom_COSEM(:,iter) = pz{oo};
end