#include "functions.hpp"
 //Use ArrayFire namespace for convenience
using namespace af;

void computeOSEstimatesIter(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const scalarStruct& inputScalars, const uint32_t iter,
	const std::vector<float>& beta, const TVdata& data, std::vector<std::vector<float*>>& imEstimates, ProjectorClass& proj, const af::array& g) {
	uint64_t yy = 0u;

	if (MethodList.LSQR) {
		if (DEBUG) {
			mexPrintf("Starting LSQR\n");
			mexEvalString("pause(.0001);");
		}
		LSQR(vec.im_os, vec.rhs_os[0], inputScalars, w_vec, iter, vec);
	}
	if (MethodList.CGLS) {
		CGLS(vec.im_os, vec.rhs_os[0], inputScalars, w_vec, iter, vec);
	}
	if (MethodList.CPLS || MethodList.CPLSKL) {
		CPLS(vec.im_os, vec.rhs_os[0], inputScalars, w_vec, vec);
	}
	if (MethodList.CPTV || MethodList.CPTVKL) {
		CPTV(vec.im_os, vec.rhs_os[0], inputScalars, w_vec, vec, proj);
	}
	// Compute BSREM and ROSEMMAP updates if applicable
	// Otherwise simply save the current iterate
	uint32_t it = 0U;
	uint32_t dd = w_vec.nMLEM;
	if (inputScalars.saveIter) {
		it = iter + 1U;
		if (MethodList.OSEM || MethodList.ECOSEM) {
			if (MethodList.OSEM) {
				if (inputScalars.use_psf && w_vec.deconvolution) {
					af::array apu = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
					deblur(apu, g, inputScalars, w_vec);
					apu.host(imEstimates[dd][it]);
				}
				else
					vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).host(imEstimates[dd][it]);
				dd++;
			}
			yy += inputScalars.im_dim;
		}

		if (MethodList.MRAMLA) {
			if (inputScalars.use_psf && w_vec.deconvolution) {
				af::array apu = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
				deblur(apu, g, inputScalars, w_vec);
				apu.host(imEstimates[dd][it]);
			}
			else
				vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).host(imEstimates[dd][it]);
			//imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
			yy += inputScalars.im_dim;
			dd++;
		}

		if (MethodList.RAMLA) {
			if (inputScalars.use_psf && w_vec.deconvolution) {
				af::array apu = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
				deblur(apu, g, inputScalars, w_vec);
				apu.host(imEstimates[dd][it]);
			}
			else
				vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).host(imEstimates[dd][it]);
			//imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
			yy += inputScalars.im_dim;
			dd++;
		}

		if (MethodList.ROSEM) {
			if (inputScalars.use_psf && w_vec.deconvolution) {
				af::array apu = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
				deblur(apu, g, inputScalars, w_vec);
				apu.host(imEstimates[dd][it]);
			}
			else
				vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).host(imEstimates[dd][it]);
			//imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
			yy += inputScalars.im_dim;
			dd++;
		}

		if (MethodList.RBI) {
			if (inputScalars.use_psf && w_vec.deconvolution) {
				af::array apu = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
				deblur(apu, g, inputScalars, w_vec);
				apu.host(imEstimates[dd][it]);
			}
			else
				vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).host(imEstimates[dd][it]);
			//imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
			yy += inputScalars.im_dim;
			dd++;
		}

		if (MethodList.DRAMA) {
			if (inputScalars.use_psf && w_vec.deconvolution) {
				af::array apu = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
				deblur(apu, g, inputScalars, w_vec);
				apu.host(imEstimates[dd][it]);
			}
			else
				vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).host(imEstimates[dd][it]);
			//imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
			yy += inputScalars.im_dim;
			dd++;
		}

		if (MethodList.COSEM || MethodList.ECOSEM) {
			if (MethodList.COSEM) {
				if (inputScalars.use_psf && w_vec.deconvolution) {
					af::array apu = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
					deblur(apu, g, inputScalars, w_vec);
					apu.host(imEstimates[dd][it]);
				}
				else
					vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).host(imEstimates[dd][it]);
				//imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
				dd++;
			}
			yy += inputScalars.im_dim;
		}

		if (MethodList.ECOSEM) {
			if (inputScalars.use_psf && w_vec.deconvolution) {
				af::array apu = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
				deblur(apu, g, inputScalars, w_vec);
				apu.host(imEstimates[dd][it]);
			}
			else
				vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).host(imEstimates[dd][it]);
			//imEstimates[dd](span, it) = vec.im_os(seq(inputScalars.im_dim * (inputScalars.nRekos2 - 1u), inputScalars.im_dim * inputScalars.nRekos2 - 1u)).copy();
			//yy += inputScalars.im_dim;
			dd++;
		}

		if (MethodList.ACOSEM) {
			if (inputScalars.use_psf && w_vec.deconvolution) {
				af::array apu = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
				deblur(apu, g, inputScalars, w_vec);
				apu.host(imEstimates[dd][it]);
			}
			else
				vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).host(imEstimates[dd][it]);
			//imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
			yy += inputScalars.im_dim;
			dd++;
		}
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
			if (MethodListPrior.MRP && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (inputScalars.osa_iter0 == inputScalars.subsets && MethodList.CUSTOM))) {
				dU = MRP(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.tr_offsets, w_vec.med_no_norm, proj);
			}
			else if (MethodListPrior.Quad && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (inputScalars.osa_iter0 == inputScalars.subsets && MethodList.CUSTOM))) {
				dU = Quadratic_prior(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.inffi, w_vec.tr_offsets, w_vec.weights_quad);
			}
			else if (MethodListPrior.Huber && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (inputScalars.osa_iter0 == inputScalars.subsets && MethodList.CUSTOM))) {
				dU = Huber_prior(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.inffi,
					w_vec.tr_offsets, w_vec.weights_huber, w_vec.huber_delta);
			}
			else if (MethodListPrior.L && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (inputScalars.osa_iter0 == inputScalars.subsets && MethodList.CUSTOM))) {
				dU = L_filter(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.tr_offsets,
					w_vec.a_L, w_vec.med_no_norm);
			}
			else if (MethodListPrior.FMH && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (inputScalars.osa_iter0 == inputScalars.subsets && MethodList.CUSTOM))) {
				dU = FMH(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.inffi, w_vec.tr_offsets,
					w_vec.fmh_weights, w_vec.med_no_norm, w_vec.alku_fmh);
			}
			else if (MethodListPrior.WeightedMean && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (inputScalars.osa_iter0 == inputScalars.subsets && MethodList.CUSTOM))) {
				dU = Weighted_mean(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.weighted_weights, w_vec.med_no_norm,
					w_vec.mean_type, w_vec.w_sum);
			}
			else if (MethodListPrior.TV && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (inputScalars.osa_iter0 == inputScalars.subsets && MethodList.CUSTOM))) {
				dU = TVprior(inputScalars, data, vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), data.TVtype, w_vec, w_vec.tr_offsets);
			}
			else if (MethodListPrior.AD && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (inputScalars.osa_iter0 == inputScalars.subsets && MethodList.CUSTOM))) {
				dU = AD(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), inputScalars, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType,
					w_vec.DiffusionType, w_vec.med_no_norm);
			}
			else if (MethodListPrior.APLS && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (inputScalars.osa_iter0 == inputScalars.subsets && MethodList.CUSTOM))) {
				dU = TVprior(inputScalars, data, vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), 5U, w_vec, w_vec.tr_offsets);
			}
			else if (MethodListPrior.TGV && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (inputScalars.osa_iter0 == inputScalars.subsets && MethodList.CUSTOM))) {
				dU = TGV(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), inputScalars, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
			}
			else if (MethodListPrior.NLM && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (inputScalars.osa_iter0 == inputScalars.subsets && MethodList.CUSTOM))) {
				dU = NLM(proj, vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), w_vec, inputScalars);
			}
			else if (MethodListPrior.RDP && (ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && (!MethodList.CUSTOM || (inputScalars.osa_iter0 == inputScalars.subsets && MethodList.CUSTOM))) {
				dU = RDP(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.weights_RDP,
					w_vec.RDP_gamma, w_vec.tr_offsets, w_vec.inffi);
			}
			else if (MethodListPrior.CUSTOM) {
				if ((ll == w_vec.mIt[0] || ll == w_vec.mIt[1]) && inputScalars.osa_iter0 == inputScalars.subsets)
					dU = w_vec.dU[oo];
				oo++;
			}

			// MAP/Prior-algorithms
			// Special case for BSREM and ROSEM-MAP
			if (MethodListMAP.BSREM && ll == w_vec.mIt[0]) {
				vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.lambda_BSREM[iter], beta[ss], dU, inputScalars.epps);
				if (inputScalars.saveIter) {
					if (inputScalars.use_psf && w_vec.deconvolution) {
						af::array apu = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
						deblur(apu, g, inputScalars, w_vec);
						apu.host(imEstimates[dd][it]);
					}
					else
						vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).host(imEstimates[dd][it]);
				}
				//imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
				MethodListMAP.BSREM = false;
			}
			else if (MethodListMAP.ROSEMMAP && ll == w_vec.mIt[1]) {
				vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)) = MAP(vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)), w_vec.lambda_ROSEM[iter], beta[ss], dU, inputScalars.epps);
				if (inputScalars.saveIter) {
					if (inputScalars.use_psf && w_vec.deconvolution) {
						af::array apu = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
						deblur(apu, g, inputScalars, w_vec);
						apu.host(imEstimates[dd][it]);
					}
					else
						vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).host(imEstimates[dd][it]);
				}
				//imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
				MethodListMAP.ROSEMMAP = false;
			}
			else {
				if (inputScalars.saveIter) {
					if (inputScalars.use_psf && w_vec.deconvolution) {
						af::array apu = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
						deblur(apu, g, inputScalars, w_vec);
						apu.host(imEstimates[dd][it]);
					}
					else
						vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).host(imEstimates[dd][it]);
				}
				//imEstimates[dd](span, it) = vec.im_os(seq(yy, yy + inputScalars.im_dim - 1u)).copy();
			}
			if (DEBUG) {
				mexPrintf("[dd] = %u\n", dd);
				//mexPrintf("imEstimates[dd](span, it) = %f\n", af::sum<float>(imEstimates[dd](span, it)));
				//mexPrintf("beta[ss] = %f\n", beta[ss]);
				//mexPrintf("[ss] = %u\n", ss);
				mexEvalString("pause(.0001);");
			}
			dd++;
			ss++;
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
		}
		dd += w_vec.nMAPML;
		ss += w_vec.nMAPML;
	}
}