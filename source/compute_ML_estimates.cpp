#include "functions.hpp"
// Use ArrayFire namespace for convenience
using namespace af;

void computeMLEstimates(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const uint32_t im_dim, const float epps, const uint32_t iter, const uint32_t subsets, 
	const std::vector<float>& beta, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const TVdata& data, const af::array& Summ_mlem, bool& break_iter, const kernelStruct& OpenCLStruct, 
	const bool saveIter)
{
	uint32_t ee = 0u;
	uint32_t it = 0U;
	uint32_t dd = 0U;
	if (saveIter)
		it = iter + 1U;
	// Compute the new estimates
	// MLEM
	if (MethodList.MLEM) {
		vec.im_mlem(seq(ee, ee + im_dim - 1u)) = EM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), Summ_mlem, vec.rhs_mlem(seq(ee, ee + im_dim - 1u)));
		vec.imEstimates[dd](span, it) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
		ee += im_dim;
		dd++;
	}
	dd += w_vec.nOS;

	uint32_t ss = 0U;
	RecMethods MethodListPrior = MethodList;
	array dU;
	for (uint32_t kk = 0; kk < w_vec.nPriorsTot; kk++) {
		RecMethods MethodListMAP = MethodList;
		for (uint32_t ll = 0; ll < w_vec.nMAPML; ll++) {
			// PRIORS
			if (MethodListPrior.MRP) {
				dU = MRP(vec.im_mlem(seq(ee, ee + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
					w_vec.med_no_norm, im_dim, OpenCLStruct);
			}
			else if (MethodListPrior.Quad) {
				dU = Quadratic_prior(vec.im_mlem(seq(ee, ee + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi,
					w_vec.tr_offsets, w_vec.weights_quad, im_dim);
			}
			else if (MethodListPrior.Huber) {
				dU = Huber_prior(vec.im_mlem(seq(ee, ee + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi,
					w_vec.tr_offsets, w_vec.weights_huber, im_dim, w_vec.huber_delta);
			}
			else if (MethodListPrior.L) {
				dU = L_filter(vec.im_mlem(seq(ee, ee + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
					w_vec.a_L, w_vec.med_no_norm, im_dim);
			}
			else if (MethodListPrior.FMH) {
				dU = FMH(vec.im_mlem(seq(ee, ee + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets,
					w_vec.fmh_weights, w_vec.med_no_norm, w_vec.alku_fmh, im_dim);
			}
			else if (MethodListPrior.WeightedMean) {
				dU = Weighted_mean(vec.im_mlem(seq(ee, ee + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.weighted_weights, w_vec.med_no_norm,
					im_dim, w_vec.mean_type, w_vec.w_sum);
			}
			else if (MethodListPrior.TV) {
				dU = TVprior(Nx, Ny, Nz, data, vec.im_mlem(seq(ee, ee + im_dim - 1u)), epps, data.TVtype, w_vec, w_vec.tr_offsets);
			}
			else if (MethodListPrior.AD) {
				if (iter == 0u) {
					dU = af::constant(0.f, im_dim, 1);
				}
				else {
					dU = AD(vec.im_mlem(seq(ee, ee + im_dim - 1u)), Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType,
						w_vec.DiffusionType, w_vec.med_no_norm);
				}
			}
			else if (MethodListPrior.APLS) {
				dU = TVprior(Nx, Ny, Nz, data, vec.im_mlem(seq(ee, ee + im_dim - 1u)), epps, 5U, w_vec, w_vec.tr_offsets);
			}
			else if (MethodListPrior.TGV) {
				dU = TGV(vec.im_mlem(seq(ee, ee + im_dim - 1u)), Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
			}
			else if (MethodListPrior.NLM) {
				dU = NLM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), w_vec, epps, Nx, Ny, Nz, OpenCLStruct);
			}
			else if (MethodListPrior.RDP) {
				dU = RDP(vec.im_mlem(seq(ee, ee + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.weights_RDP, im_dim, w_vec.RDP_gamma, w_vec.tr_offsets, w_vec.inffi);
			}
			else if (MethodListPrior.CUSTOM) {
				dU = w_vec.dU[ll];
			}
			//if (MethodListMAP.OSLMLEM) {
				vec.im_mlem(seq(ee, ee + im_dim - 1u)) = EM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), OSL(Summ_mlem, dU, beta[ss], epps), vec.rhs_mlem(seq(ee, ee + im_dim - 1u)));
				vec.imEstimates[dd](span, it) = vec.im_mlem(seq(ee, ee + im_dim - 1u)).copy();
				//MethodListMAP.OSLMLEM = false;
			//}
				if (DEBUG) {
					mexPrintf("beta[ss] = %f\n", beta[ss]);
					mexPrintf("[dd] = %u\n", dd);
					mexPrintf("[ss] = %u\n", ss);
					mexEvalString("pause(.0001);");
				}
			dd++;
			ss++;
			ee += im_dim;
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
		else if (MethodListPrior.RDP) {
			MethodListPrior.RDP = false;
		}
		else if (MethodListPrior.CUSTOM) {
			MethodListPrior.CUSTOM = false;
		}
		dd += w_vec.nMAPOS;
		ss += w_vec.nMAPOS;
	}
	if (MethodList.CUSTOM) {
		break_iter = true;
	}
}
