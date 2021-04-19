#include "AF_opencl_functions.hpp"
// Use ArrayFire namespace for convenience
using namespace af;

void computeOSEstimates(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const uint32_t im_dim, array* testi, const float epps, 
	const uint32_t iter, const uint32_t osa_iter, const uint32_t subsets, const std::vector<float>& beta, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, 
	const TVdata& data, std::vector<size_t>& length, std::vector<cl::Buffer>& d_Sino, bool& break_iter, array& pj3, const uint32_t n_rekos2, const int64_t* pituus, 
	const std::vector<cl::Buffer>& d_lor, const std::vector<cl::Buffer>& d_zindex, const std::vector<cl::Buffer>& d_xyindex, cl::Program& program_mbsrem, const cl::CommandQueue& af_queue,
	const cl::Context& af_context, std::vector<af::array>& Summ, cl::Kernel& kernel_mramla, const std::vector<cl::Buffer>& d_L, const uint8_t raw,
	const RecMethodsOpenCL& MethodListOpenCL, const size_t koko, const bool atomic_64bit, const bool atomic_32bit, const cl_uchar compute_norm_matrix, cl::Kernel& kernelNLM,
	const std::vector<cl::Buffer>& d_sc_ra, cl_uint kernelInd_MRAMLA, af::array& E, const std::vector<cl::Buffer>& d_norm, const std::vector<cl::Buffer>& d_scat, const bool use_psf,
	const af::array& g, const kernelStruct& OpenCLStruct, const bool TOF, const bool loadTOF, const mxArray* Sin, const int64_t nBins, const bool randoms_correction, const size_t local_size,
	const bool CT) {

	uint64_t yy = 0u;
	float uu;
	if (MethodList.ACOSEM || MethodList.OSLCOSEM == 1u) {
		const af::array u1 = afcl::array(length[osa_iter], d_Sino[osa_iter](), f32, true);
		clRetainMemObject(d_Sino[osa_iter]());
		uu = af::sum<float>(u1);
		//if (DEBUG) {
		//	mexPrintf("uu = %f\n", uu);
		//	mexEvalString("pause(.0001);");
		//}
	}

	// Compute the (matrix free) algorithms
	// Ordered Subsets Expectation Maximization (OSEM)
	if (MethodList.OSEM || MethodList.ECOSEM) {
		vec.im_os(seq(yy, yy + im_dim - 1u)) = EM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)));
		yy += im_dim;
	}

	// Modfied Row-action Maximum Likelihood (MRAMLA)
	if (MethodList.MRAMLA) {
		vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U,
			pj3, w_vec.lambda_MBSREM, iter, im_dim, 0.f, af::constant(0.f, 1, 1), *testi, epps);
		yy += im_dim;
	}

	// Row-action Maximum Likelihood (RAMLA)
	if (MethodList.RAMLA) {
		vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)),
			w_vec.lambda_BSREM, iter, *testi);
		yy += im_dim;
	}

	// Relaxed OSEM (ROSEM)
	if (MethodList.ROSEM) {
		vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
			w_vec.lambda_ROSEM, iter);
		yy += im_dim;
	}

	// Rescaled Block Iterative EM (RBI)
	if (MethodList.RBI) {
		vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.D);
		yy += im_dim;
	}

	// Dynamic RAMLA
	if (MethodList.DRAMA) {
		vec.im_os(seq(yy, yy + im_dim - 1u)) = DRAMA(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
			w_vec.lambda_DRAMA, iter, osa_iter, subsets);
		yy += im_dim;
	}

	// Complete data OSEM
	if (MethodList.COSEM || MethodList.ECOSEM) {
		vec.C_co(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * vec.im_os(seq(yy, yy + im_dim - 1u));
		vec.im_os(seq(yy, yy + im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_co, w_vec.D, w_vec.h_ACOSEM, 2u);
		yy += im_dim;
	}

	// Enhanced COSEM
	if (MethodList.ECOSEM) {
		vec.im_os(seq(im_dim * (n_rekos2 - 1u), im_dim * n_rekos2 - 1u)) = ECOSEM(vec.im_os(seq(im_dim * (n_rekos2 - 1u), im_dim * n_rekos2 - 1u)),
			w_vec.D, vec.im_os(seq(0, im_dim - 1u)), vec.im_os(seq(yy - im_dim, yy - 1u)), epps);
		//yy += im_dim;
	}

	// Accelerated COSEM
	if (MethodList.ACOSEM) {
		vec.C_aco(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * pow(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.h_ACOSEM_2);
		if (DEBUG) {
			mexPrintf("D = %f\n", af::sum<float>(w_vec.D));
			mexPrintf("C_aco = %f\n", af::sum<float>(af::sum(vec.C_aco,1) / w_vec.D));
			mexPrintf("C_aco = %f\n", af::sum<float>(vec.C_aco,1));
			mexPrintf("im_os = %f\n", af::sum<float>(vec.im_os));
			mexPrintf("min(D) = %f\n", af::min<float>(w_vec.D));
			mexPrintf("h = %f\n", w_vec.h_ACOSEM);
		}
		vec.im_os(seq(yy, yy + im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_aco, w_vec.D, w_vec.h_ACOSEM, 1u);
		array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
		MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, w_vec, Summ, d_Sino,
			koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length,
			atomic_64bit, atomic_32bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E, d_norm, d_scat, use_psf, g, Nx, Ny, Nz, epps, TOF, loadTOF, Sin, nBins, 
			koko, randoms_correction, local_size);
		w_vec.ACOSEM_rhs = w_vec.ACOSEM_rhs < epps ? epps : w_vec.ACOSEM_rhs;
		if (CT)
			vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (w_vec.ACOSEM_rhs / uu);
		else
			vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (uu / w_vec.ACOSEM_rhs);
		yy += im_dim;
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
				dU = MRP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
					w_vec.med_no_norm, im_dim, OpenCLStruct);
			}
			else if (MethodListPrior.Quad && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				dU = Quadratic_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi,
					w_vec.tr_offsets, w_vec.weights_quad, im_dim);
			}
			else if (MethodListPrior.Huber && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				dU = Huber_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi,
					w_vec.tr_offsets, w_vec.weights_huber, im_dim, w_vec.huber_delta);
			}
			else if (MethodListPrior.L && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				dU = L_filter(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
					w_vec.a_L, w_vec.med_no_norm, im_dim);
			}
			else if (MethodListPrior.FMH && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				dU = FMH(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets,
					w_vec.fmh_weights, w_vec.med_no_norm, w_vec.alku_fmh, im_dim);
			}
			else if (MethodListPrior.WeightedMean && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				dU = Weighted_mean(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.weighted_weights, w_vec.med_no_norm,
					im_dim, w_vec.mean_type, w_vec.w_sum);
			}
			else if (MethodListPrior.TV && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, data.TVtype, w_vec, w_vec.tr_offsets);
			}
			else if (MethodListPrior.AD && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				if (osa_iter == 0u) {
					dU = af::constant(0.f, im_dim, 1);
				}
				else {
					dU = AD(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType,
						w_vec.DiffusionType, w_vec.med_no_norm);
				}
			}
			else if (MethodListPrior.APLS && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, 5U, w_vec, w_vec.tr_offsets);
			}
			else if (MethodListPrior.TGV && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				dU = TGV(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
			}
			else if (MethodListPrior.NLM && ll != w_vec.mIt[0] && ll != w_vec.mIt[1]) {
				dU = NLM(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec, epps, Nx, Ny, Nz, OpenCLStruct);
			}
			else if (MethodListPrior.CUSTOM) {
				if (ll != w_vec.mIt[0] && ll != w_vec.mIt[1])
					dU = w_vec.dU[oo];
				oo++;
			}

			if (DEBUG) {
				mexPrintf("vec.im_os(seq(yy, yy + im_dim - 1u)) = %f\n", af::sum<float>(vec.im_os(seq(yy, yy + im_dim - 1u))));
				mexPrintf("beta.size() = %d\n", beta.size());
				mexPrintf("vec.im_os.dims(0) = %d\n", vec.im_os.dims(0));
				mexPrintf("vec.rhs_os.dims(0) = %d\n", vec.rhs_os.dims(0));
				mexPrintf("dU.dims(0) = %d\n", dU.dims(0));
				mexPrintf("[dd] = %u\n", dd);
				mexPrintf("beta[dd] = %f\n", beta[dd]);
				mexEvalString("pause(.0001);");
			}
			// MAP/Prior-algorithms
			if (MethodListMAP.OSLOSEM) {
				vec.im_os(seq(yy, yy + im_dim - 1u)) = EM(vec.im_os(seq(yy, yy + im_dim - 1u)), OSL(*testi, dU, beta[dd], epps), vec.rhs_os(seq(yy, yy + im_dim - 1u)));
				MethodListMAP.OSLOSEM = false;
				//vec.im_os(seq(yy + im_dim, yy + im_dim * 2 - 1u)) = dU;
			}
			else if (MethodListMAP.BSREM) {
				vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)),
					w_vec.lambda_BSREM, iter, *testi);
				MethodListMAP.BSREM = false;
			}
			else if (MethodListMAP.MBSREM) {
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U,
					pj3, w_vec.lambda_MBSREM, iter, im_dim, beta[dd], dU, *testi, epps);
				MethodListMAP.MBSREM = false;
			}
			else if (MethodListMAP.ROSEMMAP) {
				vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
					w_vec.lambda_ROSEM, iter);
				MethodListMAP.ROSEMMAP = false;
			}
			else if (MethodListMAP.RBIOSL) {
				vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
					w_vec.D, beta[dd], dU);
				MethodListMAP.RBIOSL = false;
			}
			else if (MethodListMAP.OSLCOSEM > 0u) {
				if (MethodList.OSLCOSEM == 1u)
					vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * pow(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.h_ACOSEM_2);
				else
					vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * vec.im_os(seq(yy, yy + im_dim - 1u));
				vec.im_os(seq(yy, yy + im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, OSL(w_vec.D, dU, beta[dd], epps), w_vec.h_ACOSEM,
					MethodList.OSLCOSEM);
				if (MethodList.OSLCOSEM == 1u) {
					array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
					if (DEBUG) {
						mexPrintf("w_vec.ACOSEM_rhs1 = %f\n", w_vec.ACOSEM_rhs);
						mexEvalString("pause(.0001);");
					}
					MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, w_vec, Summ,
						d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
						length, atomic_64bit, atomic_32bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E, d_norm, d_scat, use_psf, g, Nx, Ny, Nz, epps, TOF, loadTOF, Sin, nBins,
						koko, randoms_correction, local_size);
					//vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
					if (DEBUG) {
						mexPrintf("uu = %f\n", uu);
						mexEvalString("pause(.0001);");
					}
					if (CT)
						vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (w_vec.ACOSEM_rhs / uu);
					else
						vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (uu / w_vec.ACOSEM_rhs);
					if (DEBUG) {
						mexPrintf("w_vec.ACOSEM_rhs3 = %f\n", w_vec.ACOSEM_rhs);
						mexEvalString("pause(.0001);");
					}
				}
				MethodListMAP.OSLCOSEM = false;
			}
			else if (MethodListMAP.PKMA) {
				vec.im_os(seq(yy, yy + im_dim - 1u)) = PKMA(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_PKMA, w_vec.alpha_PKMA, 
					w_vec.sigma_PKMA, pj3, iter, osa_iter, subsets, epps, beta[dd], dU);
			}
			if (DEBUG) {
				mexPrintf("vec.rhs_os(seq(yy, yy + im_dim - 1u)) = %f\n", af::sum<float>(vec.rhs_os(seq(yy, yy + im_dim - 1u))));
				mexEvalString("pause(.0001);");
			}
			dd++;
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
			break_iter = true;
		}
		dd += w_vec.nMAPML;
	}
}