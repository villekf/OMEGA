#include "functions.hpp"
// Use ArrayFire namespace for convenience
using namespace af;

void computeMLEstimates(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const uint32_t im_dim, const float epps, const uint32_t iter, const uint32_t subsets, 
	const Beta& beta, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const TVdata& data, const af::array& Summ_mlem, bool& break_iter, const kernelStruct& OpenCLStruct, 
	const bool saveIter)
{
	uint32_t ee = 0u;
	uint32_t it = 0U;
	if (saveIter)
		it = iter + 1U;
	// Compute the new estimates
	// MLEM
	if (MethodList.MLEM) {
		vec.im_mlem(seq(ee, ee + im_dim - 1u)) = EM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), Summ_mlem, vec.rhs_mlem(seq(ee, ee + im_dim - 1u)));
		ee += im_dim;
	}
	// OSL-MLEM with Median Root Prior
	if (MethodList.OSLMLEM && MethodList.MRP) {
		array dU = MRP(vec.im_mlem(seq(ee, ee + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.med_no_norm,
			im_dim, OpenCLStruct);
		vec.im_mlem(seq(ee, ee + im_dim - 1u)) = EM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), OSL(Summ_mlem, dU, beta.MRP_MLEM, epps), vec.rhs_mlem(seq(ee, ee + im_dim - 1u)));
		ee += im_dim;
	}
	// OSL-MLEM with Quadratic prior
	if (MethodList.OSLMLEM && MethodList.Quad) {
		array dU = Quadratic_prior(vec.im_mlem(seq(ee, ee + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, w_vec.tr_offsets,
			w_vec.weights_quad, im_dim);
		vec.im_mlem(seq(ee, ee + im_dim - 1u)) = EM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), OSL(Summ_mlem, dU, beta.Quad_MLEM, epps), vec.rhs_mlem(seq(ee, ee + im_dim - 1u)));
		ee += im_dim;
	}
	// OSL-MLEM with Huber prior
	if (MethodList.OSLMLEM && MethodList.Huber) {
		array dU = Huber_prior(vec.im_mlem(seq(ee, ee + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, w_vec.tr_offsets,
			w_vec.weights_huber, im_dim, w_vec.huber_delta);
		vec.im_mlem(seq(ee, ee + im_dim - 1u)) = EM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), OSL(Summ_mlem, dU, beta.Huber_MLEM, epps), vec.rhs_mlem(seq(ee, ee + im_dim - 1u)));
		ee += im_dim;
	}
	// OSL-MLEM with L-filter prior
	if (MethodList.OSLMLEM && MethodList.L) {
		array dU = L_filter(vec.im_mlem(seq(ee, ee + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L,
			w_vec.med_no_norm, im_dim);
		vec.im_mlem(seq(ee, ee + im_dim - 1u)) = EM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), OSL(Summ_mlem, dU, beta.L_MLEM, epps), vec.rhs_mlem(seq(ee, ee + im_dim - 1u)));
		ee += im_dim;
	}
	// OSL-MLEM with FIR Median Hybrid prior
	if (MethodList.OSLMLEM && MethodList.FMH) {
		array dU = FMH(vec.im_mlem(seq(ee, ee + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets,
			w_vec.fmh_weights, w_vec.med_no_norm, w_vec.alku_fmh, im_dim);
		vec.im_mlem(seq(ee, ee + im_dim - 1u)) = EM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), OSL(Summ_mlem, dU, beta.FMH_MLEM, epps), vec.rhs_mlem(seq(ee, ee + im_dim - 1u)));
		ee += im_dim;
	}
	// OSL-MLEM with Weighted Mean prior
	if (MethodList.OSLMLEM && MethodList.WeightedMean) {
		array dU = Weighted_mean(vec.im_mlem(seq(ee, ee + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.weighted_weights, 
			w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
		vec.im_mlem(seq(ee, ee + im_dim - 1u)) = EM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), OSL(Summ_mlem, dU, beta.Weighted_MLEM, epps), vec.rhs_mlem(seq(ee, ee + im_dim - 1u)));
		ee += im_dim;
	}
	// OSL-MLEM with Total Variation prior
	if (MethodList.OSLMLEM && MethodList.TV) {
		array dU = TVprior(Nx, Ny, Nz, data, vec.im_mlem(seq(ee, ee + im_dim - 1u)), epps, data.TVtype, w_vec, w_vec.tr_offsets);
		vec.im_mlem(seq(ee, ee + im_dim - 1u)) = EM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), OSL(Summ_mlem, dU, beta.TV_MLEM, epps), vec.rhs_mlem(seq(ee, ee + im_dim - 1u)));
	}
	// OSL-MLEM with Anisotropic Diffusion smoothing prior
	if (MethodList.OSLMLEM && MethodList.AD) {
		array dU = AD(vec.im_mlem(seq(ee, ee + im_dim - 1u)), Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType,
			w_vec.DiffusionType, w_vec.med_no_norm);
		vec.im_mlem(seq(ee, ee + im_dim - 1u)) = EM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), OSL(Summ_mlem, dU, beta.AD_MLEM, epps), vec.rhs_mlem(seq(ee, ee + im_dim - 1u)));
		ee += im_dim;
	}
	// OSL-MLEM with Asymmetric Parallel Level Sets prior
	if (MethodList.OSLMLEM && MethodList.APLS) {
		array dU = TVprior(Nx, Ny, Nz, data, vec.im_mlem(seq(ee, ee + im_dim - 1u)), epps, 4, w_vec, w_vec.tr_offsets);
		vec.im_mlem(seq(ee, ee + im_dim - 1u)) = EM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), OSL(Summ_mlem, dU, beta.APLS_MLEM, epps), vec.rhs_mlem(seq(ee, ee + im_dim - 1u)));
		ee += im_dim;
	}
	// OSL-MLEM with Total Generalized Variation prior
	if (MethodList.OSLMLEM && MethodList.TGV) {
		array dU = TGV(vec.im_mlem(seq(ee, ee + im_dim - 1u)), Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
		vec.im_mlem(seq(ee, ee + im_dim - 1u)) = EM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), OSL(Summ_mlem, dU, beta.TGV_MLEM, epps), vec.rhs_mlem(seq(ee, ee + im_dim - 1u)));
		ee += im_dim;
	}
	// OSL-MLEM with Non-Local Means prior
	if (MethodList.OSLMLEM && MethodList.NLM) {
		array dU = NLM(vec.im_os(seq(ee, ee + im_dim - 1u)), w_vec, epps, Nx, Ny, Nz, OpenCLStruct);
		vec.im_mlem(seq(ee, ee + im_dim - 1u)) = EM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), OSL(Summ_mlem, dU, beta.NLM_MLEM, epps), vec.rhs_mlem(seq(ee, ee + im_dim - 1u)));
		ee += im_dim;
	}
	// OSL-MLEM with custom prior
	if (MethodList.OSLMLEM && MethodList.CUSTOM) {
		vec.im_mlem(seq(ee, ee + im_dim - 1u)) = EM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), OSL(Summ_mlem, w_vec.dU_MLEM, beta.custom_MLEM, epps), vec.rhs_mlem(seq(ee, ee + im_dim - 1u)));
		ee += im_dim;
	}

	ee = 0u;
	//Save the current iteration
	if (MethodList.MLEM) {
		vec.MLEM(span, it) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
		ee += im_dim;
	}
	if (MethodList.OSLMLEM) {
		if (MethodList.MRP) {
			vec.MRP_MLEM(span, it) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
			ee += im_dim;
		}
		if (MethodList.Quad) {
			vec.Quad_MLEM(span, it) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
			ee += im_dim;
		}
		if (MethodList.Huber) {
			vec.Huber_MLEM(span, it) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
			ee += im_dim;
		}
		if (MethodList.L) {
			vec.L_MLEM(span, it) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
			ee += im_dim;
		}
		if (MethodList.FMH) {
			vec.FMH_MLEM(span, it) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
			ee += im_dim;
		}
		if (MethodList.WeightedMean) {
			vec.Weighted_MLEM(span, it) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
			ee += im_dim;
		}
		if (MethodList.TV) {
			vec.TV_MLEM(span, it) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
			ee += im_dim;
		}
		if (MethodList.AD) {
			vec.AD_MLEM(span, it) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
			ee += im_dim;
		}
		if (MethodList.APLS) {
			vec.APLS_MLEM(span, it) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
			ee += im_dim;
		}
		if (MethodList.TGV) {
			vec.TGV_MLEM(span, it) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
			ee += im_dim;
		}
		if (MethodList.NLM) {
			vec.NLM_MLEM(span, it) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
			ee += im_dim;
		}
	}
	if (MethodList.CUSTOM) {
		if (MethodList.OSLMLEM) {
			vec.custom_MLEM = vec.im_os(seq(ee, ee + im_dim - 1u)).copy();
			ee += im_dim;
		}
		break_iter = true;
	}
}
