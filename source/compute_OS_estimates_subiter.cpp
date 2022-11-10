#include "functions.hpp"
// Use ArrayFire namespace for convenience
using namespace af;

void computeOSEstimates(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, af::array* testi, const uint32_t iter,
	const uint32_t osa_iter, const scalarStruct& inputScalars, const std::vector<float>& beta, const TVdata& data, std::vector<int64_t>& length, 
	bool& break_iter, const int64_t* pituus, std::vector<af::array>& Summ, af::array& E, const af::array& g, const af::array& D, const mxArray* Sin, 
	ProjectorClass& proj, const af::array& mData, const uint64_t m_size, const uint32_t subSum) {

	uint64_t yy = 0u;
	int status = 0;
	//float uu;
	//if (MethodList.ACOSEM || MethodList.OSLCOSEM == 1u) {

	//	//const af::array u1 = afcl::array(length[osa_iter], d_Sino[osa_iter](), f32, true);
	//	//clRetainMemObject(d_Sino[osa_iter]());
	//	//uu = af::sum<float>(u1);
	//	//if (DEBUG) {
	//	//	mexPrintf("uu = %f\n", uu);
	//	//	mexEvalString("pause(.0001);");
	//	//}
	//}

	// Compute the (matrix free) algorithms
	// Ordered Subsets Expectation Maximization (OSEM)
	if (MethodList.OSEM || MethodList.ECOSEM) {
		vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = EM(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), *testi, vec.rhs_os[0](seq(yy, yy + inputScalars.im_dim - 1u)));
		yy += inputScalars.im_dim;
	}

	// Modfied Row-action Maximum Likelihood (MRAMLA)
	if (MethodList.MRAMLA) {
		vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), vec.rhs_os[0](seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.U,
			D, w_vec.lambda_MBSREM, iter, 0.f, af::constant(0.f, 1, 1), *testi, inputScalars);
		yy += inputScalars.im_dim;
	}

	// Row-action Maximum Likelihood (RAMLA)
	if (MethodList.RAMLA) {
		vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), vec.rhs_os[0](seq(yy, yy + inputScalars.im_dim - 1u)),
			w_vec.lambda_BSREM, iter, *testi);
		yy += inputScalars.im_dim;
	}

	// Relaxed OSEM (ROSEM)
	if (MethodList.ROSEM) {
		vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), *testi, vec.rhs_os[0](seq(yy, yy + inputScalars.im_dim - 1u)),
			w_vec.lambda_ROSEM, iter);
		yy += inputScalars.im_dim;
	}

	// Rescaled Block Iterative EM (RBI)
	if (MethodList.RBI) {
		vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), *testi, vec.rhs_os[0](seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.D);
		yy += inputScalars.im_dim;
	}

	// Dynamic RAMLA
	if (MethodList.DRAMA) {
		vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = DRAMA(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), *testi, vec.rhs_os[0](seq(yy, yy + inputScalars.im_dim - 1u)),
			w_vec.lambda_DRAMA, iter, osa_iter, inputScalars.subsets);
		yy += inputScalars.im_dim;
	}

	// Complete data OSEM
	if (MethodList.COSEM || MethodList.ECOSEM) {
		vec.C_co(span, osa_iter) = vec.rhs_os[0](seq(yy, yy + inputScalars.im_dim - 1u)) * vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u));
		vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), vec.C_co, D, w_vec.h_ACOSEM, 2u);
		yy += inputScalars.im_dim;
	}

	// Enhanced COSEM
	if (MethodList.ECOSEM) {
		vec.im_os(seq(inputScalars.im_dim * (inputScalars.nRekos2 - 1u), inputScalars.im_dim * inputScalars.nRekos2 - 1u)) = 
			ECOSEM(vec.im_os(seq(inputScalars.im_dim * (inputScalars.nRekos2 - 1u), inputScalars.im_dim * inputScalars.nRekos2 - 1u)),
			w_vec.D, vec.im_os(seq(0, inputScalars.im_dim - 1u)), vec.im_os(seq(yy - inputScalars.im_dim, yy - 1u)), inputScalars.epps);
		//yy += inputScalars.im_dim;
	}

	// Accelerated COSEM
	if (MethodList.ACOSEM) {
		vec.C_aco(span, osa_iter) = vec.rhs_os[0](seq(yy, yy + inputScalars.im_dim - 1u)) * pow(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.h_ACOSEM_2);
		if (DEBUG) {
			mexPrintf("D = %f\n", af::sum<float>(D));
			mexPrintf("C_aco = %f\n", af::sum<float>(af::sum(vec.C_aco,1) / D));
			mexPrintf("C_aco = %f\n", af::sum<float>(vec.C_aco,1));
			mexPrintf("im_os = %f\n", af::sum<float>(vec.im_os));
			mexPrintf("min(D) = %f\n", af::min<float>(D));
			mexPrintf("h = %f\n", w_vec.h_ACOSEM);
		}
		float uu;
		vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), vec.C_aco, D, w_vec.h_ACOSEM, 1u); 
		if (inputScalars.use_psf)
			vec.im_os_blurred = computeConvolution(vec.im_os, g, inputScalars, w_vec, inputScalars.nRekos2);
		status = computeACOSEMWeight(inputScalars, length, uu, osa_iter, mData, m_size, w_vec, vec, proj, subSum);
		if (inputScalars.CT)
			vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) * (w_vec.ACOSEM_rhs / uu);
		else
			vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) * (uu / w_vec.ACOSEM_rhs);
		yy += inputScalars.im_dim;
	}

	RecMethods MethodListPrior = MethodList;
	uint32_t dd = w_vec.nMAPML;
	uint32_t oo = 0U;
	array dU;
	for (uint32_t kk = 0; kk < w_vec.nPriorsTot; kk++) {
		RecMethods MethodListMAP = MethodList;
		for (uint32_t ll = 0; ll < w_vec.nMAPOS; ll++) {
			//if (DEBUG) {
			//	mexPrintf("yy = %u\n", yy);
			//}
			// PRIORS
			if (MethodListPrior.MRP && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				dU = MRP(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.tr_offsets, w_vec.med_no_norm, proj);
			}
			else if (MethodListPrior.Quad && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				dU = Quadratic_prior(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.inffi,
					w_vec.tr_offsets, w_vec.weights_quad);
			}
			else if (MethodListPrior.Huber && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				dU = Huber_prior(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.inffi,
					w_vec.tr_offsets, w_vec.weights_huber, w_vec.huber_delta);
			}
			else if (MethodListPrior.L && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				dU = L_filter(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.tr_offsets,
					w_vec.a_L, w_vec.med_no_norm);
			}
			else if (MethodListPrior.FMH && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				dU = FMH(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.inffi, w_vec.tr_offsets,
					w_vec.fmh_weights, w_vec.med_no_norm, w_vec.alku_fmh);
			}
			else if (MethodListPrior.WeightedMean && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				dU = Weighted_mean(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.weighted_weights, w_vec.med_no_norm,
					w_vec.mean_type, w_vec.w_sum);
			}
			else if (MethodListPrior.TV && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				dU = TVprior(inputScalars, data, vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), data.TVtype, w_vec, w_vec.tr_offsets);
			}
			else if (MethodListPrior.AD && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				if (osa_iter == 0u) {
					dU = af::constant(0.f, inputScalars.im_dim, 1);
				}
				else {
					dU = AD(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), inputScalars, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType,
						w_vec.DiffusionType, w_vec.med_no_norm);
				}
			}
			else if (MethodListPrior.APLS && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				dU = TVprior(inputScalars, data, vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), 5U, w_vec, w_vec.tr_offsets);
			}
			else if (MethodListPrior.TGV && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				dU = TGV(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), inputScalars, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
			}
			else if (MethodListPrior.NLM && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				dU = NLM(proj, vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), w_vec, inputScalars);
			}
			else if (MethodListPrior.RDP && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				dU = RDP(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.weights_RDP, 
					w_vec.RDP_gamma, w_vec.tr_offsets, w_vec.inffi);
			}
			else if (MethodListPrior.CUSTOM) {
				if (ll != w_vec.mIt[0] && ll != w_vec.mIt[1])
					dU = w_vec.dU[oo];
				oo++;
			}

			if (DEBUG) {
				mexPrintf("vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = %f\n", af::sum<float>(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u))));
				mexPrintf("dU = %f\n", af::sum<float>(dU));
				mexPrintf("beta.size() = %d\n", beta.size());
				mexPrintf("vec.im_os.dims(0) = %d\n", vec.im_os.dims(0));
				mexPrintf("vec.rhs_os.dims(0) = %d\n", vec.rhs_os[0].dims(0));
				mexPrintf("dU.dims(0) = %d\n", dU.dims(0));
				mexPrintf("[dd] = %u\n", dd);
				mexPrintf("beta[dd] = %f\n", beta[dd]);
				mexEvalString("pause(.0001);");
			}
			// MAP/Prior-algorithms
			if (MethodListMAP.OSLOSEM) {
				vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = EM(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), OSL(*testi, dU, beta[dd], inputScalars.epps), vec.rhs_os[0](seq(yy, yy + inputScalars.im_dim - 1u)));
				//vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = dU;
				MethodListMAP.OSLOSEM = false;
			}
			else if (MethodListMAP.BSREM) {
				vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), vec.rhs_os[0](seq(yy, yy + inputScalars.im_dim - 1u)),
					w_vec.lambda_BSREM, iter, *testi);
				MethodListMAP.BSREM = false;
			}
			else if (MethodListMAP.MBSREM) {
				vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), vec.rhs_os[0](seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.U,
					D, w_vec.lambda_MBSREM, iter, beta[dd], dU, *testi, inputScalars);
				MethodListMAP.MBSREM = false;
			}
			else if (MethodListMAP.ROSEMMAP) {
				vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), *testi, vec.rhs_os[0](seq(yy, yy + inputScalars.im_dim - 1u)),
					w_vec.lambda_ROSEM, iter);
				MethodListMAP.ROSEMMAP = false;
			}
			else if (MethodListMAP.RBIOSL) {
				vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), *testi, vec.rhs_os[0](seq(yy, yy + inputScalars.im_dim - 1u)),
					w_vec.D, beta[dd], dU);
				MethodListMAP.RBIOSL = false;
			}
			else if (MethodListMAP.OSLCOSEM > 0u) {
				if (MethodList.OSLCOSEM == 1u)
					vec.C_osl(span, osa_iter) = vec.rhs_os[0](seq(yy, yy + inputScalars.im_dim - 1u)) * pow(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.h_ACOSEM_2);
				else
					vec.C_osl(span, osa_iter) = vec.rhs_os[0](seq(yy, yy + inputScalars.im_dim - 1u)) * vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u));
				vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), vec.C_osl, OSL(w_vec.D, dU, beta[dd], inputScalars.epps), w_vec.h_ACOSEM,
					MethodList.OSLCOSEM);
				if (MethodList.OSLCOSEM == 1u) {
					float uu = 0.f;
					//array apu = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u));
					if (inputScalars.use_psf)
						vec.im_os_blurred = computeConvolution(vec.im_os, g, inputScalars, w_vec, inputScalars.nRekos2);
					status = computeACOSEMWeight(inputScalars, length, uu, osa_iter, mData, m_size, w_vec, vec, proj, subSum);
					if (DEBUG) {
						mexPrintf("w_vec.ACOSEM_rhs1 = %f\n", w_vec.ACOSEM_rhs);
						mexEvalString("pause(.0001);");
					}
					//vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = apu;
					if (DEBUG) {
						mexPrintf("uu = %f\n", uu);
						mexEvalString("pause(.0001);");
					}
					if (inputScalars.CT)
						vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) * (w_vec.ACOSEM_rhs / uu);
					else
						vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) * (uu / w_vec.ACOSEM_rhs);
					if (DEBUG) {
						mexPrintf("w_vec.ACOSEM_rhs3 = %f\n", w_vec.ACOSEM_rhs);
						mexEvalString("pause(.0001);");
					}
				}
				MethodListMAP.OSLCOSEM = false;
			}
			else if (MethodListMAP.PKMA) {
				vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = PKMA(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), *testi, 
					vec.rhs_os[0](seq(yy, yy + inputScalars.im_dim - 1u)), w_vec, D, iter, osa_iter, inputScalars.subsets, inputScalars.epps, beta[dd], dU);
			}
			if (DEBUG) {
				mexPrintf("vec.rhs_os[0](seq(yy, yy + inputScalars.im_dim - 1u)) = %f\n", af::sum<float>(vec.rhs_os[0](seq(yy, yy + inputScalars.im_dim - 1u))));
				mexEvalString("pause(.0001);");
			}
			dd++;
			yy += inputScalars.im_dim;
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
			break_iter = true;
		}
		dd += w_vec.nMAPML;
	}
}