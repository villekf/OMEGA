#include "functions.hpp"
 //Use ArrayFire namespace for convenience
using namespace af;

void computeOSEstimatesIter(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const uint32_t im_dim, const float epps,
	const uint32_t iter, const uint32_t osa_iter0, const uint32_t subsets, const std::vector<float>& beta, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz,
	const TVdata& data, const uint32_t n_rekos2, const kernelStruct& OpenCLStruct, const bool saveIter) {
	uint32_t yy = 0u;
	// Compute BSREM and ROSEMMAP updates if applicable
	// Otherwise simply save the current iterate
	uint32_t it = 0U;
	uint32_t dd = w_vec.nMLEM;
	if (saveIter)
		it = iter + 1U;
	if (MethodList.OSEM || MethodList.ECOSEM) {
		if (MethodList.OSEM) {
			vec.imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			dd++;
		}
		yy += im_dim;
	}

	if (MethodList.MRAMLA) {
		vec.imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
		yy += im_dim;
		dd++;
	}

	if (MethodList.RAMLA) {
		vec.imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
		yy += im_dim;
		dd++;
	}

	if (MethodList.ROSEM) {
		vec.imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
		yy += im_dim;
		dd++;
	}

	if (MethodList.RBI) {
		vec.imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
		yy += im_dim;
		dd++;
	}

	if (MethodList.DRAMA) {
		vec.imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
		yy += im_dim;
		dd++;
	}

	if (MethodList.COSEM || MethodList.ECOSEM) {
		if (MethodList.COSEM) {
			vec.imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			dd++;
		}
		yy += im_dim;
	}

	if (MethodList.ECOSEM) {
		vec.imEstimates[dd](span, it) = vec.im_os(seq(im_dim * (n_rekos2 - 1u), im_dim * n_rekos2 - 1u)).copy();
		//yy += im_dim;
		dd++;
	}

	if (MethodList.ACOSEM) {
		vec.imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
		yy += im_dim;
		dd++;
	}

	dd += w_vec.nMAPML;
	uint32_t ss = w_vec.nMAPML;
	uint32_t oo = 0U;

	RecMethods MethodListPrior = MethodList;
	array dU;
	for (uint32_t kk = 0; kk < w_vec.nPriorsTot; kk++) {
		RecMethods MethodListMAP = MethodList;
		for (uint32_t ll = 0; ll < w_vec.nMAPOS; ll++) {
			// PRIORS
			if (MethodListPrior.MRP && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM))) {
				dU = MRP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
					w_vec.med_no_norm, im_dim, OpenCLStruct);
			}
			else if (MethodListPrior.Quad && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM))) {
				dU = Quadratic_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi,
					w_vec.tr_offsets, w_vec.weights_quad, im_dim);
			}
			else if (MethodListPrior.Huber && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM))) {
				dU = Huber_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi,
					w_vec.tr_offsets, w_vec.weights_huber, im_dim, w_vec.huber_delta);
			}
			else if (MethodListPrior.L && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM))) {
				dU = L_filter(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
					w_vec.a_L, w_vec.med_no_norm, im_dim);
			}
			else if (MethodListPrior.FMH && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM))) {
				dU = FMH(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets,
					w_vec.fmh_weights, w_vec.med_no_norm, w_vec.alku_fmh, im_dim);
			}
			else if (MethodListPrior.WeightedMean && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM))) {
				dU = Weighted_mean(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.weighted_weights, w_vec.med_no_norm,
					im_dim, w_vec.mean_type, w_vec.w_sum);
			}
			else if (MethodListPrior.TV && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM))) {
				dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, data.TVtype, w_vec, w_vec.tr_offsets);
			}
			else if (MethodListPrior.AD && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM))) {
				dU = AD(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType,
					w_vec.DiffusionType, w_vec.med_no_norm);
			}
			else if (MethodListPrior.APLS && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM))) {
				dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, 5U, w_vec, w_vec.tr_offsets);
			}
			else if (MethodListPrior.TGV && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM))) {
				dU = TGV(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
			}
			else if (MethodListPrior.NLM && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (osa_iter0 == subsets && MethodList.CUSTOM))) {
				dU = NLM(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec, epps, Nx, Ny, Nz, OpenCLStruct);
			}
			else if (MethodListPrior.CUSTOM) {
				if ((ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && osa_iter0 == subsets)
					dU = w_vec.dU[oo];
				oo++;
			}

			// MAP/Prior-algorithms
			// Special case for BSREM and ROSEM-MAP
			if (MethodListMAP.BSREM && ll == w_vec.mIt[0]) {
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM[iter], beta[ss], dU, epps);
				vec.imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
				MethodListMAP.BSREM = false;
			}
			else if (MethodListMAP.ROSEMMAP && ll == w_vec.mIt[1]) {
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM[iter], beta[ss], dU, epps);
				vec.imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
				MethodListMAP.ROSEMMAP = false;
			}
			else {
				vec.imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
			}
			if (DEBUG) {
				mexPrintf("[dd] = %u\n", dd);
				mexPrintf("vec.imEstimates[dd](span, it) = %f\n", af::sum<float>(vec.imEstimates[dd](span, it)));
				//mexPrintf("beta[ss] = %f\n", beta[ss]);
				//mexPrintf("[ss] = %u\n", ss);
				mexEvalString("pause(.0001);");
			}
			dd++;
			ss++;
			yy += im_dim;
		}

		if (MethodListPrior.MRP) {
			MethodListPrior.MRP = false;
		}
		else if (MethodListPrior.Quad) {
			MethodListPrior.Quad = false;
		}
		else if (MethodListPrior.Huber) {
			MethodListPrior.Huber = false;
		}
		else if (MethodListPrior.L) {
			MethodListPrior.L = false;
		}
		else if (MethodListPrior.FMH) {
			MethodListPrior.FMH = false;
		}
		else if (MethodListPrior.WeightedMean) {
			MethodListPrior.WeightedMean = false;
		}
		else if (MethodListPrior.TV) {
			MethodListPrior.TV = false;
		}
		else if (MethodListPrior.AD) {
			MethodListPrior.AD = false;
		}
		else if (MethodListPrior.APLS) {
			MethodListPrior.APLS = false;
		}
		else if (MethodListPrior.TGV) {
			MethodListPrior.TGV = false;
		}
		else if (MethodListPrior.NLM) {
			MethodListPrior.NLM = false;
		}
		else if (MethodListPrior.CUSTOM) {
			MethodListPrior.CUSTOM = false;
		}
		dd += w_vec.nMAPML;
		ss += w_vec.nMAPML;
	}
}