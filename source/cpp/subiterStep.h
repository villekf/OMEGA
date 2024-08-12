#pragma once
#include "functions.hpp"
#include "algorithms.h"
#include "priors.h"

inline int computeOSEstimates(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const uint32_t iter,
	uint32_t osa_iter, scalarStruct& inputScalars, std::vector<int64_t>& length, bool& break_iter, 
	const int64_t* pituus, const af::array& g, ProjectorClass& proj, const af::array& mData, uint64_t m_size, 
	const int64_t subSum, const uint8_t compute_norm_matrix) {

	int status = 0;

	af::array OSEMApu, COSEMApu, PDDYApu, FISTAApu;

	if (DEBUG) {
		mexPrint("Algo start\n");
		mexEval();
	}
	bool MAP = false;
	if (MethodList.MRP || MethodList.Quad || MethodList.Huber || MethodList.L || MethodList.FMH || MethodList.TV
		|| MethodList.WeightedMean || MethodList.AD || MethodList.APLS || MethodList.TGV || MethodList.NLM || MethodList.RDP || MethodList.ProxTGV
		|| MethodList.ProxTV || MethodList.ProxRDP || MethodList.ProxNLM || MethodList.GGMRF)
		MAP = true;

	if (DEBUG) {
		mexPrintBase("vec.rhs_os = %f\n", af::sum<float>(vec.rhs_os[0]));
		mexEval();
	}
	for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
		if (inputScalars.FISTAAcceleration)
			FISTAApu = vec.im_os[ii].copy();

		//af::array* Sens;
		//if (compute_norm_matrix == 1u) {
		//	Sens = &vec.Summ[ii][0];
		//}
		//else if (compute_norm_matrix == 2u) {
		//	Sens = &vec.Summ[ii][osa_iter];
		//}

		if (w_vec.precondTypeIm[5] && w_vec.filterIter > 0 && osa_iter + inputScalars.subsets * iter  == w_vec.filterIter) {
			if (inputScalars.verbose >= 3)
				mexPrint("Image-based filter iterations complete.");
			if (MethodList.CPType || MethodList.FISTA || MethodList.FISTAL1 || MethodList.ProxTGV || MethodList.ProxTV) {
				if (w_vec.sigmaCP[ii] == 1.f)
					w_vec.tauCP[ii] = w_vec.tauCP2[ii];
				else if (w_vec.sigmaCP[ii] == w_vec.tauCP[ii]) {
					w_vec.tauCP[ii] = w_vec.tauCP2[ii];
					w_vec.sigmaCP[ii] = w_vec.tauCP2[ii];
				}
				else
					w_vec.sigmaCP[ii] = w_vec.tauCP2[ii];
				if (MethodList.CPType)
					w_vec.sigma2CP = w_vec.sigmaCP;
			}
			if (MethodList.MRAMLA || MethodList.MBSREM || MethodList.SPS || MethodList.RAMLA || MethodList.BSREM || MethodList.ROSEM || MethodList.ROSEMMAP || MethodList.PKMA)
				w_vec.lambda = w_vec.lambdaFiltered;
			w_vec.precondTypeIm[5] = false;
		}
		if (MethodList.PDHG || MethodList.PDHGKL || MethodList.PDHGL1 || MethodList.CV || MethodList.PDDY)
			PDHG1(vec.rhs_os[ii], inputScalars, w_vec, vec, osa_iter + inputScalars.subsets * iter, ii);
		if (MethodList.PDDY && ii == 0 && MAP) {
			PDDYApu = vec.im_os[0].copy();
			vec.im_os[0] -= w_vec.tauCP[0] * vec.uCP[0];
			if (inputScalars.verbose == 3) {
				mexPrint("Computing PDDY step\n");
				mexEval();
			}
		}
		if (status != 0)
			return -1;
		vec.im_os[ii].eval();
		vec.rhs_os[ii].eval();
	}

	if (DEBUG) {
		mexPrint("Priori\n");
		mexEval();
	}
	if (!MethodList.BSREM && !MethodList.ROSEMMAP && !MethodList.POCS) {
		status = applyPrior(vec, w_vec, MethodList, inputScalars, proj, w_vec.beta, osa_iter + inputScalars.subsets * iter, compute_norm_matrix);
		if (status != 0)
			return -1;
	}
	if (MethodList.PDDY && MAP)
		vec.im_os[0] = PDDYApu.copy();
	//vec.im_os[0].eval();
	//if (DEBUG) {
	//	//mexPrintBase("dU.dims(0) = %d\n", dU.dims(0));
	//	//af::deviceGC();
	//	af_print_mem_info("mem info", -1);
	//	//mexEval();
	//}

	for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
		af::array* Sens = nullptr;
		if (compute_norm_matrix == 1u) {
			Sens = &vec.Summ[ii][0];
		}
		else if (compute_norm_matrix == 2u) {
			Sens = &vec.Summ[ii][osa_iter];
		}

		// Compute the (matrix free) algorithms
		// Ordered Subsets Expectation Maximization (OSEM)
		if (MethodList.OSEM || MethodList.ECOSEM) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing OSEM/ECOSEM");
			if (MethodList.ECOSEM)
				OSEMApu = EM(vec.im_os[ii], *Sens, vec.rhs_os[ii] + inputScalars.epps);
			else
				if (inputScalars.CT)
					vec.im_os[ii] = EM(vec.im_os[ii], *Sens, inputScalars.flat * vec.rhs_os[ii]);
				else
					vec.im_os[ii] = EM(vec.im_os[ii], *Sens, vec.rhs_os[ii] + inputScalars.epps);
		}

		// Modfied Row-action Maximum Likelihood (MRAMLA)
		if (MethodList.MRAMLA) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing MRAMLA");
			//vec.im_os[ii] = MBSREM(vec.im_os[ii], vec.rhs_os[ii], w_vec.U, w_vec.lambda, iter, osa_iter, inputScalars, w_vec, proj);
			status = MBSREM(vec.im_os[ii], vec.rhs_os[ii], w_vec.U, w_vec.lambda, iter, osa_iter, inputScalars, w_vec, proj, ii);
		}

		// Row-action Maximum Likelihood (RAMLA)
		if (MethodList.RAMLA || MethodList.BSREM) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing RAMLA");
			//vec.im_os[ii] = BSREM(vec.im_os[ii], vec.rhs_os[ii], w_vec.lambda, iter, *Sens);
			status = BSREM(vec.im_os[ii], vec.rhs_os[ii], w_vec.lambda, iter, inputScalars, proj, ii);
		}

		// Relaxed OSEM (ROSEM)
		if (MethodList.ROSEM) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing ROSEM");
			if (inputScalars.CT)
				vec.im_os[ii] = ROSEM(vec.im_os[ii], *Sens, inputScalars.flat * vec.rhs_os[ii], w_vec.lambda, iter);
			else
				vec.im_os[ii] = ROSEM(vec.im_os[ii], *Sens, vec.rhs_os[ii], w_vec.lambda, iter);
		}

		// Rescaled Block Iterative EM (RBI)
		if (MethodList.RBI) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing RBI");
			vec.im_os[ii] = RBI(vec.im_os[ii], *Sens, vec.rhs_os[ii], w_vec.D[ii]);
		}

		// Dynamic RAMLA
		if (MethodList.DRAMA) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing DRAMA");
			vec.im_os[ii] = DRAMA(vec.im_os[ii], *Sens, vec.rhs_os[ii], w_vec.lambda, iter, osa_iter, inputScalars.subsets);
		}

		// Complete data OSEM
		if (MethodList.COSEM || MethodList.ECOSEM) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing COSEM/ECOSEM");
			vec.C_co(af::span, osa_iter) = vec.rhs_os[ii] * vec.im_os[ii];
			if (MethodList.ECOSEM)
				COSEMApu = COSEM(vec.im_os[ii], vec.C_co, w_vec.D[ii], w_vec.h_ACOSEM, 2u);
			else
				vec.im_os[ii] = COSEM(vec.im_os[ii], vec.C_co, w_vec.D[ii], w_vec.h_ACOSEM, 2u);
		}

		// Enhanced COSEM
		if (MethodList.ECOSEM) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing ECOSEM");
			vec.im_os[ii] = ECOSEM(vec.im_os[ii], w_vec.D[ii], OSEMApu, COSEMApu, inputScalars.epps);
		}

		// Accelerated COSEM
		if (MethodList.ACOSEM) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing ACOSEM");
			if (DEBUG) {
				mexPrintBase("vec.rhs_os[ii].dims(0) = %u\n", vec.rhs_os[ii].dims(0));
				mexPrintBase("vec.im_os[ii].dims(0) = %u\n", vec.im_os[ii].dims(0));
				mexPrintBase("vec.C_co.dims(0) = %u\n", vec.C_co.dims(0));
				mexPrintBase("vec.C_co.dims(1) = %u\n", vec.C_co.dims(1));
				mexPrintBase("osa_iter = %u\n", osa_iter);
				mexEval();
			}
			vec.C_co(af::span, osa_iter) = vec.rhs_os[ii] * af::pow(vec.im_os[ii], w_vec.h_ACOSEM_2);
			if (DEBUG) {
				mexPrintBase("im_os = %f\n", af::sum<float>(vec.im_os[ii]));
				mexPrintBase("rhs_os = %f\n", af::sum<float>(vec.rhs_os[ii]));
				mexPrintBase("D = %f\n", af::sum<float>(w_vec.D[ii]));
				mexPrintBase("C_co / D = %f\n", af::sum<float>(af::sum(vec.C_co, 1) / w_vec.D[ii]));
				mexPrintBase("C_co = %f\n", af::sum<float>(vec.C_co, 1));
				mexPrintBase("C_co(:,osa_iter) = %f\n", af::sum<float>(vec.C_co(af::span, osa_iter), 1));
				mexPrintBase("min(D) = %f\n", af::min<float>(w_vec.D[ii]));
				mexPrintBase("h = %f\n", w_vec.h_ACOSEM);
				mexPrintBase("h2 = %f\n", w_vec.h_ACOSEM_2);
				mexEval();
			}
			float uu;
			vec.im_os[ii] = COSEM(vec.im_os[ii], vec.C_co, w_vec.D[ii], w_vec.h_ACOSEM, 1u);
			//if (inputScalars.use_psf)
			//	vec.im_os_blurred[ii] = computeConvolution(vec.im_os[ii], g, inputScalars, w_vec, inputScalars.nRekos2, ii);
			status = computeACOSEMWeight(inputScalars, length, uu, osa_iter, mData, m_size, w_vec, vec, proj, subSum, g);
			if (inputScalars.CT)
				vec.im_os[ii] *= (w_vec.ACOSEM_rhs / uu);
			else
				vec.im_os[ii] *= (uu / w_vec.ACOSEM_rhs);
			vec.im_os[ii].eval();
			if (DEBUG) {
				mexPrintBase("w_vec.ACOSEM_rhs = %f\n", w_vec.ACOSEM_rhs);
				mexPrintBase("uu = %f\n", uu);
				mexPrintBase("uu / w_vec.ACOSEM_rhs = %f\n", uu / w_vec.ACOSEM_rhs);
				mexEval();
			}
		}
		if (MethodList.OSLOSEM) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing OSL-OSEM");
			if (ii == 0) {
				if (inputScalars.CT)
					vec.im_os[ii] = EM(vec.im_os[ii], *Sens + vec.dU, inputScalars.flat * vec.rhs_os[ii]);
				else
					vec.im_os[ii] = EM(vec.im_os[ii], *Sens + vec.dU, vec.rhs_os[ii]);
			}
			else {

			}
		}
		else if (MethodList.BSREM) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing BSREM");
			//vec.im_os[ii] = BSREM(vec.im_os[ii], vec.rhs_os[ii],	w_vec.lambda, iter, *Sens);
			status = BSREM(vec.im_os[ii], vec.rhs_os[ii], w_vec.lambda, iter, inputScalars, proj, ii);
		}
		else if (MethodList.MBSREM) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing MBSREM");
			//vec.im_os[ii] = MBSREM(vec.im_os[ii], vec.rhs_os[ii], w_vec.U, w_vec.lambda, iter, osa_iter, inputScalars, w_vec, proj);
			status = MBSREM(vec.im_os[ii], vec.rhs_os[ii], w_vec.U, w_vec.lambda, iter, osa_iter, inputScalars, w_vec, proj, ii);
		}
		else if (MethodList.ROSEMMAP) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing ROSEMMAP");
			if (inputScalars.CT)
				vec.im_os[ii] = ROSEM(vec.im_os[ii], *Sens, inputScalars.flat * vec.rhs_os[ii], w_vec.lambda, iter);
			else
				vec.im_os[ii] = ROSEM(vec.im_os[ii], *Sens, vec.rhs_os[ii], w_vec.lambda, iter);
		}
		else if (MethodList.RBIOSL) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing RBIOSL");
			if (ii == 0)
				vec.im_os[ii] = RBI(vec.im_os[ii], *Sens, vec.rhs_os[ii], w_vec.D[ii], w_vec.beta, vec.dU);
			else
				vec.im_os[ii] = RBI(vec.im_os[ii], *Sens, vec.rhs_os[ii], w_vec.D[ii], 0.f);
		}
		else if (MethodList.OSLCOSEM > 0u) {
			if (inputScalars.verbose >= 3)
				mexPrintVar("Computing OSLCOSEM ", MethodList.OSLCOSEM);
			if (MethodList.OSLCOSEM == 1u)
				vec.C_co(af::span, osa_iter) = vec.rhs_os[ii] * pow(vec.im_os[ii], w_vec.h_ACOSEM_2);
			else
				vec.C_co(af::span, osa_iter) = vec.rhs_os[ii] * vec.im_os[ii];
			if (ii == 0)
				vec.im_os[ii] = COSEM(vec.im_os[ii], vec.C_co, w_vec.D[ii] + vec.dU, w_vec.h_ACOSEM, MethodList.OSLCOSEM);
			else
				vec.im_os[ii] = COSEM(vec.im_os[ii], vec.C_co, w_vec.D[ii], w_vec.h_ACOSEM, MethodList.OSLCOSEM);
			if (MethodList.OSLCOSEM == 1u) {
				float uu = 0.f;
				//if (inputScalars.use_psf)
				//	vec.im_os_blurred[ii] = computeConvolution(vec.im_os[ii], g, inputScalars, w_vec, inputScalars.nRekos2, ii);
				status = computeACOSEMWeight(inputScalars, length, uu, osa_iter, mData, m_size, w_vec, vec, proj, subSum, g);
				if (status != 0)
					return -1;
				if (DEBUG) {
					mexPrintBase("w_vec.ACOSEM_rhs1 = %f\n", w_vec.ACOSEM_rhs);
					mexPrintBase("uu = %f\n", uu);
					mexEval();
				}
				if (inputScalars.CT)
					vec.im_os[ii] = vec.im_os[ii] * (w_vec.ACOSEM_rhs / uu);
				else
					vec.im_os[ii] = vec.im_os[ii] * (uu / w_vec.ACOSEM_rhs);
				if (DEBUG) {
					mexPrintBase("w_vec.ACOSEM_rhs3 = %f\n", w_vec.ACOSEM_rhs);
					mexEval();
				}
			}
		}
		else if (MethodList.PKMA) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing PKMA");
			//vec.im_os[ii] = PKMA(vec.im_os[ii], vec.rhs_os[ii], w_vec, inputScalars, iter, osa_iter, proj);
			status = PKMA(vec.im_os[ii], vec.rhs_os[ii], w_vec, inputScalars, iter, osa_iter, proj, ii);
		}
		else if (MethodList.SPS) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing SPS");
			status = SPS(vec.im_os[ii], vec.rhs_os[ii], w_vec.U, w_vec.lambda, iter, osa_iter, inputScalars, w_vec, proj, ii);
		}
		else if (MethodList.LSQR) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing LSQR");
			LSQR(inputScalars, w_vec, iter, vec);
		}
		else if (MethodList.CGLS) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing CGLS");
			CGLS(inputScalars, w_vec, iter, vec);
		}
		else if (MethodList.SART || MethodList.POCS) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing SART");
			if (DEBUG) {
				mexPrintBase("w_vec.lambda[iter] = %f\n", w_vec.lambda[iter]);
				mexEval();
			}
			vec.im_os[ii] = SART(vec.im_os[ii], *Sens, vec.rhs_os[ii], w_vec.lambda[iter]);
			if (MethodList.POCS)
				POCS(vec.im_os[ii], inputScalars, w_vec, MethodList, vec, proj, mData, g, length, pituus, osa_iter, iter, ii);
		}
		//else if (MethodList.PDHG || MethodList.PDHGKL || MethodList.PDHGL1) {
		//	if (inputScalars.verbose >= 3)
		//		mexPrint("Computing PDHG/PDHGKL/PDHGL1");
		//	status = CPLS(vec.im_os[ii], vec.rhs_os[ii], inputScalars, w_vec, vec, proj, iter, osa_iter);
		//}	
		else if (MethodList.PDHG || MethodList.PDHGKL || MethodList.PDHGL1 || MethodList.CV || MethodList.PDDY) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing PDHG/PDHGKL/PDHGL1/PDDY");
			status = PDHG2(vec.im_os[ii], vec.rhs_os[ii], inputScalars, w_vec, vec, proj, iter, osa_iter, ii, pituus, g, m_size, length);
			//vec.im_os[ii] = vec.rhs_os[ii];
		}
		else if (MethodList.FISTA) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing FISTA");
			status = FISTA(vec.im_os[ii], vec.rhs_os[ii], inputScalars, w_vec, vec, proj, osa_iter + inputScalars.subsets * iter, ii);
		}
		else if (MethodList.FISTAL1) {
			if (inputScalars.verbose >= 3)
				mexPrint("Computing FISTAL1");
			status = FISTAL1(vec.im_os[ii], vec.rhs_os[ii], inputScalars, w_vec, vec, w_vec.beta, proj, osa_iter + inputScalars.subsets * iter, ii);
		}
		if (inputScalars.FISTAAcceleration) {
			//if (osa_iter == 0) {
				//const uint32_t it = iter + 1;
				//w_vec.betaFISTA = static_cast<float>(it - 1) / static_cast<float>(it + 2);
				//if (w_vec.betaFISTA <= 0.f) {
				//	w_vec.tFISTA = (1.f + std::sqrt(1.f + 4.f * w_vec.tNFista * w_vec.tNFista)) / 2.f;
				//	w_vec.betaFISTA = (w_vec.tNFista - 1.f) / w_vec.tFISTA;
				//	w_vec.tNFista = w_vec.tFISTA;
				//}
				//vec.im_os[ii] = vec.im_os[ii] + w_vec.betaFISTA * (vec.im_os[ii] - FISTAApu);
			const float t = w_vec.tFISTA;
			if (osa_iter == 0) {
				w_vec.tFISTA = (1.f + std::sqrt(1.f + 4.f * w_vec.tFISTA * w_vec.tFISTA)) / 2.f;
			}
			vec.im_os[ii] = vec.im_os[ii] + (t - 1.f) / w_vec.tFISTA * (vec.im_os[ii] - FISTAApu);
			//vec.im_os[ii](vec.im_os[ii] < inputScalars.epps) = inputScalars.epps;
			//}
			af::eval(vec.im_os[ii]);
		}
	}
	if (inputScalars.verbose >= 3)
		mexPrint("Iterative algorithm computed");
	return status;
}