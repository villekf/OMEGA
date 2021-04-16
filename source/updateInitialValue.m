function x0 = updateInitialValue(im_vectors, options)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if options.OSEM
    x0 = im_vectors.OSEM_apu;
end

if options.MLEM
    x0 = im_vectors.MLEM_apu;
end

if options.MRAMLA
    x0 = im_vectors.MRAMLA_apu;
end

if options.RAMLA
    x0 = im_vectors.RAMLA_apu;
end

if options.ROSEM
    x0 = im_vectors.ROSEM_apu;
end

if options.RBI
    x0 = im_vectors.RBI_apu;
end

if options.DRAMA
    x0 = im_vectors.DRAMA_apu;
end

if options.COSEM
    x0 = im_vectors.COSEM_apu;
end

if options.ECOSEM
    x0 = im_vectors.ECOSEM_apu;
end

if options.ACOSEM
    x0 = im_vectors.ACOSEM_apu;
end

if options.MRP && options.OSL_OSEM
    x0 = im_vectors.MRP_OSL_apu;
end
if options.MRP && options.OSL_MLEM
    x0 = im_vectors.MRP_MLEM_apu;
end
if options.MRP && options.MBSREM
    x0 = im_vectors.MRP_MBSREM_apu;
end

if options.MRP && options.BSREM
    x0 = im_vectors.MRP_BSREM_apu;
end

if options.MRP && options.ROSEM_MAP
    x0 = im_vectors.MRP_ROSEM_apu;
end

if options.MRP && options.OSL_RBI
    x0 = im_vectors.MRP_RBI_apu;
end

if options.MRP && any(options.OSL_COSEM)
    x0 = im_vectors.MRP_COSEM_apu;
end

if options.quad && options.OSL_OSEM
    x0 = im_vectors.Quad_OSL_apu;
end
if options.quad && options.OSL_MLEM
    x0 = im_vectors.Quad_MLEM_apu;
end
if options.quad && options.MBSREM
    x0 = im_vectors.Quad_MBSREM_apu;
end

if options.quad && options.BSREM
    x0 = im_vectors.Quad_BSREM_apu;
end

if options.quad && options.ROSEM_MAP
    x0 = im_vectors.Quad_ROSEM_apu;
end

if options.quad && options.OSL_RBI
    x0 = im_vectors.Quad_RBI_apu;
end

if options.quad && any(options.OSL_COSEM)
    x0 = im_vectors.Quad_COSEM_apu;
end

if options.L && options.OSL_OSEM
    x0 = im_vectors.L_OSL_apu;
end
if options.L && options.OSL_MLEM
    x0 = im_vectors.L_MLEM_apu;
end
if options.L && options.MBSREM
    x0 = im_vectors.L_MBSREM_apu;
end

if options.L && options.BSREM
    x0 = im_vectors.L_BSREM_apu;
end

if options.L && options.ROSEM_MAP
    x0 = im_vectors.L_ROSEM_apu;
end

if options.L && options.OSL_RBI
    x0 = im_vectors.L_RBI_apu;
end

if options.L && any(options.OSL_COSEM)
    x0 = im_vectors.L_COSEM_apu;
end

if options.FMH && options.OSL_OSEM
    x0 = im_vectors.FMH_OSL_apu;
end
if options.FMH && options.OSL_MLEM
    x0 = im_vectors.FMH_MLEM_apu;
end
if options.FMH && options.MBSREM
    x0 = im_vectors.FMH_MBSREM_apu;
end

if options.FMH && options.BSREM
    x0 = im_vectors.FMH_BSREM_apu;
end

if options.FMH && options.ROSEM_MAP
    x0 = im_vectors.FMH_ROSEM_apu;
end

if options.FMH && options.OSL_RBI
    x0 = im_vectors.FMH_RBI_apu;
end

if options.FMH && any(options.OSL_COSEM)
    x0 = im_vectors.FMH_COSEM_apu;
end

if options.weighted_mean && options.OSL_OSEM
    x0 = im_vectors.Weighted_OSL_apu;
end
if options.weighted_mean && options.OSL_MLEM
    x0 = im_vectors.Weighted_MLEM_apu;
end
if options.weighted_mean && options.MBSREM
    x0 = im_vectors.Weighted_MBSREM_apu;
end

if options.weighted_mean && options.BSREM
    x0 = im_vectors.Weighted_BSREM_apu;
end

if options.weighted_mean && options.ROSEM_MAP
    x0 = im_vectors.Weighted_ROSEM_apu;
end

if options.weighted_mean && options.OSL_RBI
    x0 = im_vectors.Weighted_RBI_apu;
end

if options.weighted_mean && any(options.OSL_COSEM)
    x0 = im_vectors.Weighted_COSEM_apu;
end

if options.TV && options.OSL_OSEM
    x0 = im_vectors.TV_OSL_apu;
end
if options.TV && options.OSL_MLEM
    x0 = im_vectors.TV_MLEM_apu;
end
if options.TV && options.MBSREM
    x0 = im_vectors.TV_MBSREM_apu;
end

if options.TV && options.BSREM
    x0 = im_vectors.TV_BSREM_apu;
end

if options.TV && options.ROSEM_MAP
    x0 = im_vectors.TV_ROSEM_apu;
end

if options.TV && options.OSL_RBI
    x0 = im_vectors.TV_RBI_apu;
end

if options.TV && any(options.OSL_COSEM)
    x0 = im_vectors.TV_COSEM_apu;
end

if options.AD && options.OSL_OSEM
    x0 = im_vectors.AD_OSL_apu;
end
if options.AD && options.OSL_MLEM
    x0 = im_vectors.AD_MLEM_apu;
end
if options.AD && options.MBSREM
    x0 = im_vectors.AD_MBSREM_apu;
end

if options.AD && options.BSREM
    x0 = im_vectors.AD_BSREM_apu;
end
if options.AD && options.ROSEM_MAP
    x0 = im_vectors.AD_ROSEM_apu;
end
if options.AD && options.OSL_RBI
    x0 = im_vectors.AD_RBI_apu;
end
if options.AD && any(options.OSL_COSEM)
    x0 = im_vectors.AD_COSEM_apu;
end

if options.APLS && options.OSL_OSEM
    x0 = im_vectors.APLS_OSL_apu;
end
if options.APLS && options.OSL_MLEM
    x0 = im_vectors.APLS_MLEM_apu;
end
if options.APLS && options.MBSREM
    x0 = im_vectors.APLS_MBSREM_apu;
end
if options.APLS && options.BSREM
    x0 = im_vectors.APLS_BSREM_apu;
end
if options.APLS && options.ROSEM_MAP
    x0 = im_vectors.APLS_ROSEM_apu;
end
if options.APLS && options.OSL_RBI
    x0 = im_vectors.APLS_RBI_apu;
end
if options.APLS && any(options.OSL_COSEM)
    x0 = im_vectors.APLS_COSEM_apu;
end

if options.TGV && options.OSL_OSEM
    x0 = im_vectors.TGV_OSL_apu;
end
if options.TGV && options.OSL_MLEM
    x0 = im_vectors.TGV_MLEM_apu;
end
if options.TGV && options.MBSREM
    x0 = im_vectors.TGV_MBSREM_apu;
end
if options.TGV && options.BSREM
    x0 = im_vectors.TGV_BSREM_apu;
end
if options.TGV && options.ROSEM_MAP
    x0 = im_vectors.TGV_ROSEM_apu;
end
if options.TGV && options.OSL_RBI
    x0 = im_vectors.TGV_RBI_apu;
end
if options.TGV && any(options.OSL_COSEM)
    x0 = im_vectors.TGV_COSEM_apu;
end
if options.NLM && options.OSL_OSEM
    x0 = im_vectors.NLM_OSL_apu;
end
if options.NLM && options.MBSREM
    x0 = im_vectors.NLM_MBSREM_apu;
end
if options.NLM && options.BSREM
    x0 = im_vectors.NLM_BSREM_apu;
end
if options.NLM && options.ROSEM_MAP
    x0 = im_vectors.NLM_ROSEM_apu;
end
if options.NLM && options.OSL_RBI
    x0 = im_vectors.NLM_RBI_apu;
end
if options.NLM && any(options.OSL_COSEM)
    x0 = im_vectors.NLM_COSEM_apu;
end
end

