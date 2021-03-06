function options = double_to_single(options)
%DOUBLE_TO_SINGLE Convert all necessary values in options to single
%precision
options.x0 = single(options.x0(:));
options.Ndx = uint32(options.Ndx);
options.Ndy = uint32(options.Ndy);
options.Ndz = uint32(options.Ndz);
options.h = single(options.h);
options.U = single(options.U);
options.MLEM = logical(options.MLEM);
options.OSEM = logical(options.OSEM);
options.MRAMLA = logical(options.MRAMLA);
options.RAMLA = logical(options.RAMLA);
options.ROSEM = logical(options.ROSEM);
options.RBI = logical(options.RBI);
options.DRAMA = logical(options.DRAMA);
options.COSEM = logical(options.COSEM);
options.ECOSEM = logical(options.ECOSEM);
options.ACOSEM = logical(options.ACOSEM);
options.OSL_MLEM = logical(options.OSL_MLEM);
options.OSL_OSEM = logical(options.OSL_OSEM);
options.MBSREM = logical(options.MBSREM);
options.BSREM = logical(options.BSREM);
options.ROSEM_MAP = logical(options.ROSEM_MAP);
options.OSL_RBI = logical(options.OSL_RBI);
options.MRP = logical(options.MRP);
options.quad = logical(options.quad);
options.L = logical(options.L);
options.FMH = logical(options.FMH);
options.weighted_mean = logical(options.weighted_mean);
options.TV = logical(options.TV);
options.AD = logical(options.AD);
options.APLS = logical(options.APLS);
options.TGV = logical(options.TGV);
options.OSL_COSEM = uint32(options.OSL_COSEM);
options.mean_type = int32(options.mean_type);

varPrior = recNames(1);
varMAP = recNames(2);

apu = false(numel(varPrior),1);
apu2 = false(numel(varMAP),1);
if numel(varMAP) > 0
    for kk = 1 : numel(varPrior)
        for ll = 1 : numel(varMAP)
            if options.(varPrior{kk}) && options.(varMAP{ll})
                apu(kk) = true;
                apu2(ll) = true;
            end
        end
    end
end
varPrior = varPrior(apu);
varMAP = varMAP(apu2);


for kk = 1 : numel(varPrior)
    for ll = 1 : numel(varMAP)
        options.(['beta_' [varPrior{kk} '_' varMAP{ll}]]) = single(options.(['beta_' [varPrior{kk} '_' varMAP{ll}]]));
    end
end
% options.beta_Huber_OSL_OSEM = single(options.beta_Huber_OSL_OSEM);
% options.beta_Huber_OSL_MLEM = single(options.beta_Huber_OSL_MLEM);
% options.beta_Huber_BSREM = single(options.beta_Huber_BSREM);
% options.beta_Huber_MBSREM = single(options.beta_Huber_MBSREM);
% options.beta_Huber_ROSEM_MAP = single(options.beta_Huber_ROSEM_MAP);
% options.beta_Huber_OSL_RBI = single(options.beta_Huber_OSL_RBI);
% options.beta_Huber_OSL_COSEM = single(options.beta_Huber_OSL_COSEM);
% options.beta_MRP_OSL_OSEM = single(options.beta_MRP_OSL_OSEM);
% options.beta_MRP_BSREM = single(options.beta_MRP_BSREM);
% options.beta_L_OSL_OSEM = single(options.beta_L_OSL_OSEM);
% options.beta_L_BSREM = single(options.beta_L_BSREM);
% options.beta_quad_OSL_OSEM = single(options.beta_quad_OSL_OSEM);
% options.beta_quad_BSREM = single(options.beta_quad_BSREM);
% options.beta_FMH_OSL_OSEM = single(options.beta_FMH_OSL_OSEM);
% options.beta_weighted_mean_OSL_OSEM = single(options.beta_weighted_mean_OSL_OSEM);
% options.beta_TV_OSL_OSEM = single(options.beta_TV_OSL_OSEM);
% options.beta_TGV_OSL_OSEM = single(options.beta_TGV_OSL_OSEM);
% options.beta_FMH_BSREM = single(options.beta_FMH_BSREM);
% options.beta_weighted_mean_BSREM = single(options.beta_weighted_mean_BSREM);
% options.beta_TV_BSREM = single(options.beta_TV_BSREM);
% options.beta_TGV_BSREM = single(options.beta_TGV_BSREM);
% options.beta_AD_OSL_OSEM = single(options.beta_AD_OSL_OSEM);
% options.beta_AD_BSREM = single(options.beta_AD_BSREM);
% options.beta_APLS_OSL_OSEM = single(options.beta_APLS_OSL_OSEM);
% options.beta_APLS_BSREM = single(options.beta_APLS_BSREM);
% options.beta_MRP_MBSREM = single(options.beta_MRP_MBSREM);
% options.beta_L_MBSREM = single(options.beta_L_MBSREM);
% options.beta_quad_MBSREM = single(options.beta_quad_MBSREM);
% options.beta_FMH_MBSREM = single(options.beta_FMH_MBSREM);
% options.beta_weighted_mean_MBSREM = single(options.beta_weighted_mean_MBSREM);
% options.beta_TV_MBSREM = single(options.beta_TV_MBSREM);
% options.beta_AD_MBSREM = single(options.beta_AD_MBSREM);
% options.beta_APLS_MBSREM = single(options.beta_APLS_MBSREM);
% options.beta_TGV_MBSREM = single(options.beta_TGV_MBSREM);
% options.beta_MRP_ROSEM_MAP = single(options.beta_MRP_ROSEM_MAP);
% options.beta_L_ROSEM_MAP = single(options.beta_L_ROSEM_MAP);
% options.beta_quad_ROSEM_MAP = single(options.beta_quad_ROSEM_MAP);
% options.beta_FMH_ROSEM_MAP = single(options.beta_FMH_ROSEM_MAP);
% options.beta_weighted_mean_ROSEM_MAP = single(options.beta_weighted_mean_ROSEM_MAP);
% options.beta_TV_ROSEM_MAP = single(options.beta_TV_ROSEM_MAP);
% options.beta_AD_ROSEM_MAP = single(options.beta_AD_ROSEM_MAP);
% options.beta_APLS_ROSEM_MAP = single(options.beta_APLS_ROSEM_MAP);
% options.beta_TGV_ROSEM_MAP = single(options.beta_TGV_ROSEM_MAP);
% options.beta_MRP_OSL_RBI = single(options.beta_MRP_OSL_RBI);
% options.beta_L_OSL_RBI = single(options.beta_L_OSL_RBI);
% options.beta_quad_OSL_RBI = single(options.beta_quad_OSL_RBI);
% options.beta_FMH_OSL_RBI = single(options.beta_FMH_OSL_RBI);
% options.beta_weighted_mean_OSL_RBI = single(options.beta_weighted_mean_OSL_RBI);
% options.beta_TV_OSL_RBI = single(options.beta_TV_OSL_RBI);
% options.beta_AD_OSL_RBI = single(options.beta_AD_OSL_RBI);
% options.beta_APLS_OSL_RBI = single(options.beta_APLS_OSL_RBI);
% options.beta_TGV_OSL_RBI = single(options.beta_TGV_OSL_RBI);
% options.beta_MRP_OSL_COSEM = single(options.beta_MRP_OSL_COSEM);
% options.beta_L_OSL_COSEM = single(options.beta_L_OSL_COSEM);
% options.beta_quad_OSL_COSEM = single(options.beta_quad_OSL_COSEM);
% options.beta_FMH_OSL_COSEM = single(options.beta_FMH_OSL_COSEM);
% options.beta_weighted_mean_OSL_COSEM = single(options.beta_weighted_mean_OSL_COSEM);
% options.beta_TV_OSL_COSEM = single(options.beta_TV_OSL_COSEM);
% options.beta_AD_OSL_COSEM = single(options.beta_AD_OSL_COSEM);
% options.beta_APLS_OSL_COSEM = single(options.beta_APLS_OSL_COSEM);
% options.beta_TGV_OSL_COSEM = single(options.beta_TGV_OSL_COSEM);
% options.beta_MRP_OSL_MLEM = single(options.beta_MRP_OSL_MLEM);
% options.beta_L_OSL_MLEM = single(options.beta_L_OSL_MLEM);
% options.beta_quad_OSL_MLEM = single(options.beta_quad_OSL_MLEM);
% options.beta_FMH_OSL_MLEM = single(options.beta_FMH_OSL_MLEM);
% options.beta_weighted_mean_OSL_MLEM = single(options.beta_weighted_mean_OSL_MLEM);
% options.beta_TV_OSL_MLEM = single(options.beta_TV_OSL_MLEM);
% options.beta_AD_OSL_MLEM = single(options.beta_AD_OSL_MLEM);
% options.beta_APLS_OSL_MLEM = single(options.beta_APLS_OSL_MLEM);
% options.beta_TGV_OSL_MLEM = single(options.beta_TGV_OSL_MLEM);
% options.beta_NLM_OSL_OSEM = single(options.beta_NLM_OSL_OSEM);
% options.beta_NLM_OSL_MLEM = single(options.beta_NLM_OSL_MLEM);
% options.beta_NLM_BSREM = single(options.beta_NLM_BSREM);
% options.beta_NLM_MBSREM = single(options.beta_NLM_MBSREM);
% options.beta_NLM_ROSEM_MAP = single(options.beta_NLM_ROSEM_MAP);
% options.beta_NLM_OSL_RBI = single(options.beta_NLM_OSL_RBI);
% options.beta_NLM_OSL_COSEM = single(options.beta_NLM_OSL_COSEM);
end

