#include "AF_cuda_functions.hpp"
// Use ArrayFire namespace for convenience
using namespace af;

void computeOSEstimatesCUDA(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, RecMethodsOpenCL& MethodListOpenCL, const uint32_t im_dim,
	af::array* testi, const float epps, const uint32_t iter, const uint32_t osa_iter, const uint32_t subsets, const std::vector<float>& beta, const uint32_t Nx, const uint32_t Ny,
	const uint32_t Nz, const TVdata& data, std::vector<size_t>& length, std::vector<CUdeviceptr>& d_Sino, bool& break_iter, af::array& pj3, const uint32_t n_rekos2,
	const int64_t* pituus, std::vector<CUdeviceptr>& d_lor, std::vector<CUdeviceptr>& d_zindex, std::vector<CUdeviceptr>& d_xyindex, const CUstream& af_cuda_stream,
	std::vector<af::array>& Summ, CUfunction& kernel_mramla, std::vector<CUdeviceptr>& d_L, const uint8_t raw, const size_t koko, const bool atomic_64bit,
	const uint8_t compute_norm_matrix, std::vector<CUdeviceptr>& d_sc_ra, af::array& E, std::vector<CUdeviceptr>& d_norm, std::vector<CUdeviceptr>& d_scat,
	const uint32_t det_per_ring, CUdeviceptr& d_pseudos, const uint32_t prows, const float dz, const float dx, const float dy, const float bz, const float bx,
	const float by, const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices, CUdeviceptr& d_x, CUdeviceptr& d_y, CUdeviceptr& d_z,
	const uint32_t size_x, const uint32_t TotSinos, CUdeviceptr& d_atten, const uint32_t Nxy, const float tube_width, const float crystal_size_z, const float bmin,
	const float bmax, const float Vmax, CUdeviceptr& d_xcenter, CUdeviceptr& d_ycenter, CUdeviceptr& d_zcenter, CUdeviceptr& d_V, const float dc_z, const uint16_t n_rays,
	const uint16_t n_rays3D, const bool precompute, const uint32_t projector_type, const float global_factor, CUdeviceptr& d_reko_type, CUfunction& kernel_mbsrem,
	const bool use_psf, const af::array& g, const kernelStruct& CUDAStruct, const bool TOF, const bool loadTOF, const mxArray* Sin, const int64_t nBins,
	const bool randoms_correction, const float sigma_x, CUdeviceptr& d_TOFCenter, CUdeviceptr& d_angles, const bool CT) {

	uint64_t yy = 0u;
	float a_Summa = 1.f;
	//float* hOut = reinterpret_cast<float*>(d_Sino[osa_iter]);
	if (MethodList.ACOSEM || MethodList.OSLCOSEM == 1u) {
		//hOut = reinterpret_cast<float*>(d_Sino[osa_iter]);
		float* hOut = new float[length[osa_iter]];
		cuMemcpyDtoH(hOut, d_Sino[osa_iter], length[osa_iter] * sizeof(float));
		cuCtxSynchronize();
		//array uu(length[osa_iter], d_Sino[osa_iter], afDevice);
		array uu(length[osa_iter], hOut, afHost);
		//array uu(length[osa_iter], hOut, afDevice);
		af::eval(uu);
		a_Summa = af::sum<float>(uu);
		delete[] hOut;
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
		vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.D,
			0.f, af::constant(0.f, 1, 1));
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
			mexPrintf("a_Summa = %f\n", a_Summa);
			mexPrintf("C_aco = %f\n", af::sum<float>(vec.C_aco, 1));
			mexPrintf("im_os = %f\n", af::sum<float>(vec.im_os));
		}
		vec.im_os(seq(yy, yy + im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_aco, w_vec.D, w_vec.h_ACOSEM, 1u);
		array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
		MRAMLA_prepass_CUDA(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, w_vec, Summ, d_Sino, koko, apu, vec.C_co,
			vec.C_aco, vec.C_osl, osa_iter + 1u, d_L, raw, MethodListOpenCL, length, compute_norm_matrix, d_sc_ra, E, det_per_ring, d_pseudos,
			prows, Nx, Ny, Nz, dz, dx, dy, bz, bx, by, bzb, maxxx, maxyy, zmax, NSlices, d_x, d_y, d_z, size_x, TotSinos,
			d_atten, d_norm, d_scat, epps, Nxy, tube_width, crystal_size_z, bmin, bmax, Vmax,
			d_xcenter, d_ycenter, d_zcenter, d_V, dc_z, n_rays, n_rays3D, precompute, projector_type, af_cuda_stream, global_factor, d_reko_type,
			kernel_mbsrem, atomic_64bit, use_psf, g, TOF, loadTOF, Sin, nBins, randoms_correction, sigma_x, d_TOFCenter, d_angles, 1U, CT);
		w_vec.ACOSEM_rhs = w_vec.ACOSEM_rhs < epps ? epps : w_vec.ACOSEM_rhs;
		if (CT)
			vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (w_vec.ACOSEM_rhs / a_Summa);
		else
			if (CT)
				vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (w_vec.ACOSEM_rhs / a_Summa);
			else
				vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (a_Summa / w_vec.ACOSEM_rhs);
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
					w_vec.med_no_norm, im_dim, CUDAStruct);
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
				dU = NLM(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec, epps, Nx, Ny, Nz, CUDAStruct);
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
					MRAMLA_prepass_CUDA(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, w_vec, Summ, d_Sino, koko, apu, vec.C_co,
						vec.C_aco, vec.C_osl, osa_iter + 1u, d_L, raw, MethodListOpenCL, length, compute_norm_matrix, d_sc_ra, E, det_per_ring, d_pseudos,
						prows, Nx, Ny, Nz, dz, dx, dy, bz, bx, by, bzb, maxxx, maxyy, zmax, NSlices, d_x, d_y, d_z, size_x, TotSinos,
						d_atten, d_norm, d_scat, epps, Nxy, tube_width, crystal_size_z, bmin, bmax, Vmax,
						d_xcenter, d_ycenter, d_zcenter, d_V, dc_z, n_rays, n_rays3D, precompute, projector_type, af_cuda_stream, global_factor, d_reko_type,
						kernel_mbsrem, atomic_64bit, use_psf, g, TOF, loadTOF, Sin, nBins, randoms_correction, sigma_x, d_TOFCenter, d_angles, 1U, CT);
					//vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
					if (CT)
						vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (w_vec.ACOSEM_rhs / a_Summa);
					else
						vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (a_Summa / w_vec.ACOSEM_rhs);
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