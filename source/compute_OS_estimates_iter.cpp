#include "functions.hpp"
 //Use ArrayFire namespace for convenience
using namespace af;

void computeOSEstimatesIter(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const uint32_t im_dim, const float epps,
	const uint32_t iter, const uint32_t osa_iter0, const uint32_t subsets, const Beta& beta, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz,
	const TVdata& data, const uint32_t n_rekos2, const kernelStruct& OpenCLStruct, const bool saveIter) {
	uint32_t yy = 0u;
	// Compute BSREM and ROSEMMAP updates if applicable
	// Otherwise simply save the current iterate
	uint32_t it = 0U;
	if (saveIter)
		it = iter + 1U;
	if (MethodList.OSEM || MethodList.ECOSEM) {
		if (MethodList.OSEM)
			vec.OSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
		yy += im_dim;
	}

	if (MethodList.MRAMLA) {
		vec.MRAMLA(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
		yy += im_dim;
	}

	if (MethodList.RAMLA) {
		vec.RAMLA(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
		yy += im_dim;
	}

	if (MethodList.ROSEM) {
		vec.ROSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
		yy += im_dim;
	}

	if (MethodList.RBI) {
		vec.RBI(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
		yy += im_dim;
	}

	if (MethodList.DRAMA) {
		vec.DRAMA(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
		yy += im_dim;
	}

	if (MethodList.COSEM || MethodList.ECOSEM) {
		if (MethodList.COSEM)
			vec.COSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
		yy += im_dim;
	}

	if (MethodList.ECOSEM) {
		vec.ECOSEM(span, it) = vec.im_os(seq(im_dim * (n_rekos2 - 1u), im_dim * n_rekos2 - 1u)).copy();
		//yy += im_dim;
	}

	if (MethodList.ACOSEM) {
		vec.ACOSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
		yy += im_dim;
	}

	if (MethodList.MRP) {
		if (MethodList.OSLOSEM) {
			vec.MRP_OSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = MRP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
					w_vec.med_no_norm, im_dim, OpenCLStruct);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM[iter], beta.MRP_BSREM, dU, epps);
			}
			vec.MRP_BSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			vec.MRP_MBSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = MRP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
					w_vec.med_no_norm, im_dim, OpenCLStruct);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM[iter], beta.MRP_ROSEM, dU, epps);
			}
			vec.MRP_ROSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			vec.MRP_RBI(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.MRP_COSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
	}

	if (MethodList.Quad) {
		if (MethodList.OSLOSEM) {
			vec.Quad_OSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = Quadratic_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi,
					w_vec.tr_offsets, w_vec.weights_quad, im_dim);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM[iter], beta.Quad_BSREM, dU, epps);
			}
			vec.Quad_BSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			vec.Quad_MBSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = Quadratic_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi,
					w_vec.tr_offsets, w_vec.weights_quad, im_dim);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM[iter], beta.Quad_ROSEM, dU, epps);
			}
			vec.Quad_ROSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			vec.Quad_RBI(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.Quad_COSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
	}

	if (MethodList.Huber) {
		if (MethodList.OSLOSEM) {
			vec.Huber_OSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = Huber_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi,
					w_vec.tr_offsets, w_vec.weights_huber, im_dim, w_vec.huber_delta);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM[iter], beta.Huber_BSREM, dU, epps);
			}
			vec.Huber_BSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			vec.Huber_MBSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = Huber_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi,
					w_vec.tr_offsets, w_vec.weights_huber, im_dim, w_vec.huber_delta);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM[iter], beta.Huber_ROSEM, dU, epps);
			}
			vec.Huber_ROSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			vec.Huber_RBI(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.Huber_COSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
	}

	if (MethodList.L) {
		if (MethodList.OSLOSEM) {
			vec.L_OSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = L_filter(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L,
					w_vec.med_no_norm, im_dim);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM[iter], beta.L_BSREM, dU, epps);
			}
			vec.L_BSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			vec.L_MBSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = L_filter(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L,
					w_vec.med_no_norm, im_dim);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM[iter], beta.L_ROSEM, dU, epps);
			}
			vec.L_ROSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			vec.L_RBI(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.L_COSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
	}
	if (MethodList.FMH) {
		if (MethodList.OSLOSEM) {
			vec.FMH_OSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = FMH(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets,
					w_vec.fmh_weights, w_vec.med_no_norm, w_vec.alku_fmh, im_dim);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM[iter], beta.FMH_BSREM, dU, epps);
			}
			vec.FMH_BSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			vec.FMH_MBSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = FMH(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets,
					w_vec.fmh_weights, w_vec.med_no_norm, w_vec.alku_fmh, im_dim);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM[iter], beta.FMH_ROSEM, dU, epps);
			}
			vec.FMH_ROSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			vec.FMH_RBI(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.FMH_COSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
	}
	if (MethodList.WeightedMean) {
		if (MethodList.OSLOSEM) {
			vec.Weighted_OSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = Weighted_mean(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.weighted_weights,
					w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM[iter], beta.Weighted_BSREM, dU,
					epps);
			}
			vec.Weighted_BSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			vec.Weighted_MBSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = Weighted_mean(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.weighted_weights,
					w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM[iter], beta.Weighted_ROSEM, dU,
					epps);
			}
			vec.Weighted_ROSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			vec.Weighted_RBI(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.Weighted_COSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
	}
	if (MethodList.TV) {
		if (MethodList.OSLOSEM) {
			vec.TV_OSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, data.TVtype, w_vec, w_vec.tr_offsets);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM[iter], beta.TV_BSREM, dU, epps);
			}
			vec.TV_BSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			vec.TV_MBSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, data.TVtype, w_vec, w_vec.tr_offsets);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM[iter], beta.TV_ROSEM, dU, epps);
			}
			vec.TV_ROSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			vec.TV_RBI(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.TV_COSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
	}
	if (MethodList.AD) {
		if (MethodList.OSLOSEM) {
			vec.AD_OSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = AD(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType,
					w_vec.DiffusionType, w_vec.med_no_norm);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM[iter], beta.AD_BSREM, dU, epps);
			}
			vec.AD_BSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			vec.AD_MBSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = AD(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType,
					w_vec.DiffusionType, w_vec.med_no_norm);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM[iter], beta.AD_ROSEM, dU, epps);
			}
			vec.AD_ROSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			vec.AD_RBI(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.AD_COSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
	}
	if (MethodList.APLS) {
		if (MethodList.OSLOSEM) {
			vec.APLS_OSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, 5U, w_vec, w_vec.tr_offsets);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM[iter], beta.APLS_BSREM, dU, epps);
			}
			vec.APLS_BSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			vec.APLS_MBSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, 5U, w_vec, w_vec.tr_offsets);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM[iter], beta.APLS_ROSEM, dU, epps);
			}
			vec.APLS_ROSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			vec.APLS_RBI(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.APLS_COSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
	}
	if (MethodList.TGV) {
		if (MethodList.OSLOSEM) {
			vec.TGV_OSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = TGV(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM[iter], beta.TGV_BSREM, dU, epps);
			}
			vec.TGV_BSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			vec.TGV_MBSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = TGV(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM[iter], beta.TGV_ROSEM, dU, epps);
			}
			vec.TGV_ROSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			vec.TGV_RBI(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.TGV_COSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
	}
	if (MethodList.NLM) {
		if (MethodList.OSLOSEM) {
			vec.NLM_OSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = NLM(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec, epps, Nx, Ny, Nz, OpenCLStruct);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM[iter], beta.NLM_BSREM, dU, epps);
			}
			vec.NLM_BSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			vec.NLM_MBSREM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			if (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM)) {
				const array dU = NLM(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec, epps, Nx, Ny, Nz, OpenCLStruct);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM[iter], beta.NLM_ROSEM, dU, epps);
			}
			vec.NLM_ROSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			vec.NLM_RBI(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.NLM_COSEM(span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
	}
	if (MethodList.CUSTOM) {
		if (MethodList.OSLOSEM) {
			vec.custom_OSEM = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.BSREM) {
		//	if (osa_iter0 == subsets)
		//		vec.custom_BSREM = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM[iter], beta.custom_BSREM, w_vec.dU_BSREM, epps);
		//	else
				vec.custom_BSREM = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			vec.custom_MBSREM = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			//if (osa_iter0 == subsets)
			//	vec.custom_ROSEM = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM[iter], beta.custom_ROSEM, w_vec.dU_ROSEM, epps);
			//else
				vec.custom_ROSEM = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			vec.custom_RBI = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0) {
			vec.custom_COSEM = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			yy += im_dim;
		}
	}
}