#include "AF_opencl_functions.hpp"
// Use ArrayFire namespace for convenience
using namespace af;

void computeOSEstimates(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const uint32_t im_dim, array* testi, const float epps, 
	const uint32_t iter, const uint32_t osa_iter, const uint32_t subsets, const Beta& beta, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, 
	const TVdata& data, std::vector<size_t>& length, std::vector<cl::Buffer>& d_Sino, bool& break_iter, array& pj3, const uint32_t n_rekos2, const int64_t* pituus, 
	const std::vector<cl::Buffer>& d_lor, const std::vector<cl::Buffer>& d_zindex, const std::vector<cl::Buffer>& d_xyindex, cl::Program& program_mbsrem, const cl::CommandQueue& af_queue,
	const cl::Context& af_context, std::vector<af::array>& Summ, cl::Kernel& kernel_mramla, const std::vector<cl::Buffer>& d_L, const uint8_t raw,
	const RecMethodsOpenCL& MethodListOpenCL, const size_t koko, const bool atomic_64bit, const cl_uchar compute_norm_matrix, cl::Kernel& kernelNLM, 
	const std::vector<cl::Buffer>& d_sc_ra, cl_uint kernelInd_MRAMLA, af::array& E, const std::vector<cl::Buffer>& d_norm, const std::vector<cl::Buffer>& d_scat, const bool use_psf,
	const af::array& g, const kernelStruct& OpenCLStruct, const bool TOF, const bool loadTOF, const mxArray* Sin, const int64_t nBins, const bool randoms_correction) {

	uint64_t yy = 0u;
	float uu;
	if (MethodList.ACOSEM || MethodList.OSLCOSEM == 1u) {
		array u1 = afcl::array(length[osa_iter], d_Sino[osa_iter].get(), f32, true);
		uu = af::sum<float>(u1);
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
			w_vec.lambda_BSREM, iter);
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
		vec.im_os(seq(yy, yy + im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_aco, w_vec.D, w_vec.h_ACOSEM, 1u);
		array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
		MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, w_vec, Summ, d_Sino,
			koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length,
			atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E, d_norm, d_scat, use_psf, g, Nx, Ny, Nz, epps, TOF, loadTOF, Sin, nBins, 
			koko, randoms_correction);
		w_vec.ACOSEM_rhs = w_vec.ACOSEM_rhs < epps ? epps : w_vec.ACOSEM_rhs;
		vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (uu / w_vec.ACOSEM_rhs);
		yy += im_dim;
	}
	// Median Root Prior
	if (MethodList.MRP) {
		// OSL-OSEM
		if (MethodList.OSLOSEM) {
			array dU = MRP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
				w_vec.med_no_norm, im_dim);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = EM(vec.im_os(seq(yy, yy + im_dim - 1u)), OSL(*testi, dU, beta.MRP_OSEM, epps), vec.rhs_os(seq(yy, yy + im_dim - 1u)));
			yy += im_dim;
		}
		if (MethodList.MAPOSEM) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = EM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)));
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_BSREM, iter);
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			array dU = MRP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
				w_vec.med_no_norm, im_dim);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U,
				pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.MRP_MBSREM, dU, *testi, epps);
			yy += im_dim;
		}
		if (MethodList.ROSEMOSL) {
			array dU = MRP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
				w_vec.med_no_norm, im_dim);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), OSL(*testi, dU, beta.MRP_ROSEMOSL, epps), vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_ROSEM, iter);
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_ROSEM, iter);
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			array dU = MRP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
				w_vec.med_no_norm, im_dim);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.D, beta.MRP_RBI, dU);
			yy += im_dim;
		}
		if (MethodList.RBIMAP) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)));
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0u) {
			array dU = MRP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
				w_vec.med_no_norm, im_dim);
			if (MethodList.OSLCOSEM == 1u)
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * pow(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.h_ACOSEM_2);
			else
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * vec.im_os(seq(yy, yy + im_dim - 1u));
			vec.im_os(seq(yy, yy + im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, OSL(w_vec.D, dU, beta.MRP_COSEM, epps), w_vec.h_ACOSEM,
				MethodList.OSLCOSEM);
			if (MethodList.OSLCOSEM == 1u) {
				array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
				MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, w_vec, Summ,
					d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
					length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E, d_norm, d_scat, use_psf, g, Nx, Ny, Nz, epps, TOF, loadTOF, Sin, nBins, 
					koko, randoms_correction);
				//vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
				vec.im_os(seq(yy, yy + im_dim - 1u)) = batchFunc(vec.im_os(seq(yy, yy + im_dim - 1u)), uu / w_vec.ACOSEM_rhs, batchMul);
			}
			yy += im_dim;
		}
		if (MethodList.MAPCOSEM > 0u) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, w_vec.D, w_vec.h_ACOSEM,
				MethodList.MAPCOSEM);
			if (MethodList.MAPCOSEM == 1u) {
				array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
				MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, w_vec, Summ,
					d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
					length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E, d_norm, d_scat, use_psf, g, Nx, Ny, Nz, epps, TOF, loadTOF, Sin, nBins, 
					koko, randoms_correction);
				//vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
				vec.im_os(seq(yy, yy + im_dim - 1u)) = batchFunc(vec.im_os(seq(yy, yy + im_dim - 1u)), uu / w_vec.ACOSEM_rhs, batchMul);
			}
			yy += im_dim;
		}
	}
	// Quadratic Prior
	if (MethodList.Quad) {
		if (MethodList.OSLOSEM) {
			array dU = Quadratic_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi,
				w_vec.tr_offsets, w_vec.weights_quad, im_dim);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = EM(vec.im_os(seq(yy, yy + im_dim - 1u)), OSL(*testi, dU, beta.Quad_OSEM, epps), vec.rhs_os(seq(yy, yy + im_dim - 1u)));
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_BSREM, iter);
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			array dU = Quadratic_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi,
				w_vec.tr_offsets, w_vec.weights_quad, im_dim);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U,
				pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.Quad_MBSREM, dU, *testi, epps);
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_ROSEM, iter);
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			array dU = Quadratic_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi,
				w_vec.tr_offsets, w_vec.weights_quad, im_dim);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.D, beta.Quad_RBI, dU);
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0u) {
			array dU = Quadratic_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi,
				w_vec.tr_offsets, w_vec.weights_quad, im_dim);
			if (MethodList.OSLCOSEM == 1u)
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * pow(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.h_ACOSEM_2);
			else
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * vec.im_os(seq(yy, yy + im_dim - 1u));
			vec.im_os(seq(yy, yy + im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, OSL(w_vec.D, dU, beta.Quad_COSEM, epps), w_vec.h_ACOSEM,
				MethodList.OSLCOSEM);
			if (MethodList.OSLCOSEM == 1u) {
				array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
				MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, w_vec, Summ,
					d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
					length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E, d_norm, d_scat, use_psf, g, Nx, Ny, Nz, epps, TOF, loadTOF, Sin, nBins, 
					koko, randoms_correction);
				//vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
				w_vec.ACOSEM_rhs = w_vec.ACOSEM_rhs < epps ? epps : w_vec.ACOSEM_rhs;
				vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (uu / w_vec.ACOSEM_rhs);
			}
			yy += im_dim;
		}
	}
	// Huber Prior
	if (MethodList.Huber) {
		if (MethodList.OSLOSEM) {
			array dU = Huber_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi,
				w_vec.tr_offsets, w_vec.weights_huber, im_dim, w_vec.huber_delta);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = EM(vec.im_os(seq(yy, yy + im_dim - 1u)), OSL(*testi, dU, beta.Huber_OSEM, epps), vec.rhs_os(seq(yy, yy + im_dim - 1u)));
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_BSREM, iter);
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			array dU = Huber_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi,
				w_vec.tr_offsets, w_vec.weights_huber, im_dim, w_vec.huber_delta);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U,
				pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.Huber_MBSREM, dU, *testi, epps);
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_ROSEM, iter);
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			array dU = Huber_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi,
				w_vec.tr_offsets, w_vec.weights_huber, im_dim, w_vec.huber_delta);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.D, beta.Huber_RBI, dU);
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0u) {
			array dU = Huber_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi,
				w_vec.tr_offsets, w_vec.weights_huber, im_dim, w_vec.huber_delta);
			if (MethodList.OSLCOSEM == 1u)
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * pow(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.h_ACOSEM_2);
			else
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * vec.im_os(seq(yy, yy + im_dim - 1u));
			vec.im_os(seq(yy, yy + im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, OSL(w_vec.D, dU, beta.Huber_COSEM, epps), w_vec.h_ACOSEM,
				MethodList.OSLCOSEM);
			if (MethodList.OSLCOSEM == 1u) {
				array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
				MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, w_vec, Summ,
					d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
					length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E, d_norm, d_scat, use_psf, g, Nx, Ny, Nz, epps, TOF, loadTOF, Sin, nBins, 
					koko, randoms_correction);
				//vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
				w_vec.ACOSEM_rhs = w_vec.ACOSEM_rhs < epps ? epps : w_vec.ACOSEM_rhs;
				vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (uu / w_vec.ACOSEM_rhs);
			}
			yy += im_dim;
		}
	}
	// L-filter prior
	if (MethodList.L) {
		if (MethodList.OSLOSEM) {
			array dU = L_filter(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
				w_vec.a_L, w_vec.med_no_norm, im_dim);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = EM(vec.im_os(seq(yy, yy + im_dim - 1u)), OSL(*testi, dU, beta.L_OSEM, epps), vec.rhs_os(seq(yy, yy + im_dim - 1u)));
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_BSREM, iter);
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			array dU = L_filter(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
				w_vec.a_L, w_vec.med_no_norm, im_dim);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U,
				pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.L_MBSREM, dU, *testi, epps);
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_ROSEM, iter);
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			array dU = L_filter(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
				w_vec.a_L, w_vec.med_no_norm, im_dim);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.D, beta.L_RBI, dU);
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0u) {
			array dU = L_filter(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
				w_vec.a_L, w_vec.med_no_norm, im_dim);
			if (MethodList.OSLCOSEM == 1u)
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * pow(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.h_ACOSEM_2);
			else
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * vec.im_os(seq(yy, yy + im_dim - 1u));
			vec.im_os(seq(yy, yy + im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, OSL(w_vec.D, dU, beta.L_COSEM, epps), w_vec.h_ACOSEM,
				MethodList.OSLCOSEM);
			if (MethodList.OSLCOSEM == 1u) {
				array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
				MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, w_vec, Summ,
					d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
					length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E, d_norm, d_scat, use_psf, g, Nx, Ny, Nz, epps, TOF, loadTOF, Sin, nBins, 
					koko, randoms_correction);
				//vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
				w_vec.ACOSEM_rhs = w_vec.ACOSEM_rhs < epps ? epps : w_vec.ACOSEM_rhs;
				vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (uu / w_vec.ACOSEM_rhs);
			}
			yy += im_dim;
		}
	}
	// FIR Median Hybrid prior
	if (MethodList.FMH) {
		if (MethodList.OSLOSEM) {
			array dU = FMH(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets,
				w_vec.fmh_weights, w_vec.med_no_norm, w_vec.alku_fmh, im_dim);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = EM(vec.im_os(seq(yy, yy + im_dim - 1u)), OSL(*testi, dU, beta.FMH_OSEM, epps), vec.rhs_os(seq(yy, yy + im_dim - 1u)));
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_BSREM, iter);
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			array dU = FMH(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets,
				w_vec.fmh_weights, w_vec.med_no_norm, w_vec.alku_fmh, im_dim);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U,
				pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.FMH_MBSREM, dU, *testi, epps);
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_ROSEM, iter);
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			array dU = FMH(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets,
				w_vec.fmh_weights, w_vec.med_no_norm, w_vec.alku_fmh, im_dim);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.D, beta.FMH_RBI, dU);
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0u) {
			array dU = FMH(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets,
				w_vec.fmh_weights, w_vec.med_no_norm, w_vec.alku_fmh, im_dim);
			if (MethodList.OSLCOSEM == 1u)
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * pow(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.h_ACOSEM_2);
			else
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * vec.im_os(seq(yy, yy + im_dim - 1u));
			vec.im_os(seq(yy, yy + im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, OSL(w_vec.D, dU, beta.FMH_COSEM, epps), w_vec.h_ACOSEM,
				MethodList.OSLCOSEM);
			if (MethodList.OSLCOSEM == 1u) {
				array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
				MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, w_vec, Summ,
					d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
					length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E, d_norm, d_scat, use_psf, g, Nx, Ny, Nz, epps, TOF, loadTOF, Sin, nBins, 
					koko, randoms_correction);
				//vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
				w_vec.ACOSEM_rhs = w_vec.ACOSEM_rhs < epps ? epps : w_vec.ACOSEM_rhs;
				vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (uu / w_vec.ACOSEM_rhs);
			}
			yy += im_dim;
		}
	}
	// Weighted Mean prior
	if (MethodList.WeightedMean) {
		if (MethodList.OSLOSEM) {
			array dU = Weighted_mean(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.weighted_weights, w_vec.med_no_norm, 
				im_dim, w_vec.mean_type, w_vec.w_sum);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = EM(vec.im_os(seq(yy, yy + im_dim - 1u)), OSL(*testi, dU, beta.Weighted_OSEM, epps), vec.rhs_os(seq(yy, yy + im_dim - 1u)));
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_BSREM, iter);
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			array dU = Weighted_mean(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.weighted_weights, w_vec.med_no_norm, 
				im_dim, w_vec.mean_type, w_vec.w_sum);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U,
				pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.Weighted_MBSREM, dU, *testi, epps);
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_ROSEM, iter);
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			array dU = Weighted_mean(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.weighted_weights, w_vec.med_no_norm, 
				im_dim, w_vec.mean_type, w_vec.w_sum);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.D, beta.Weighted_RBI, dU);
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0u) {
			array dU = Weighted_mean(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.weighted_weights, w_vec.med_no_norm, 
				im_dim, w_vec.mean_type, w_vec.w_sum);
			if (MethodList.OSLCOSEM == 1u)
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * pow(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.h_ACOSEM_2);
			else
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * vec.im_os(seq(yy, yy + im_dim - 1u));
			vec.im_os(seq(yy, yy + im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, OSL(w_vec.D, dU, beta.Weighted_COSEM, epps), w_vec.h_ACOSEM,
				MethodList.OSLCOSEM);
			if (MethodList.OSLCOSEM == 1u) {
				array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
				MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, w_vec, Summ,
					d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw,
					MethodListOpenCL, length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E, d_norm, d_scat, use_psf, g, Nx, Ny, Nz, epps, TOF, loadTOF, Sin, nBins, 
					koko, randoms_correction);
				//vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
				w_vec.ACOSEM_rhs = w_vec.ACOSEM_rhs < epps ? epps : w_vec.ACOSEM_rhs;
				vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (uu / w_vec.ACOSEM_rhs);
			}
			yy += im_dim;
		}
	}
	// Total Variation prior
	if (MethodList.TV) {
		if (MethodList.OSLOSEM) {
			array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, data.TVtype, w_vec, w_vec.tr_offsets);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = EM(vec.im_os(seq(yy, yy + im_dim - 1u)), OSL(*testi, dU, beta.TV_OSEM, epps), vec.rhs_os(seq(yy, yy + im_dim - 1u)));
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_BSREM, iter);
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, data.TVtype, w_vec, w_vec.tr_offsets);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U,
				pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.TV_MBSREM, dU, *testi, epps);
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_ROSEM, iter);
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, data.TVtype, w_vec, w_vec.tr_offsets);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.D, beta.TV_RBI, dU);
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0u) {
			array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, data.TVtype, w_vec, w_vec.tr_offsets);
			if (MethodList.OSLCOSEM == 1u)
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * pow(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.h_ACOSEM_2);
			else
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * vec.im_os(seq(yy, yy + im_dim - 1u));
			vec.im_os(seq(yy, yy + im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, OSL(w_vec.D, dU, beta.TV_COSEM, epps), w_vec.h_ACOSEM,
				MethodList.OSLCOSEM);
			if (MethodList.OSLCOSEM == 1u) {
				array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
				MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, w_vec, Summ,
					d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
					length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E, d_norm, d_scat, use_psf, g, Nx, Ny, Nz, epps, TOF, loadTOF, Sin, nBins, 
					koko, randoms_correction);
				//vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
				w_vec.ACOSEM_rhs = w_vec.ACOSEM_rhs < epps ? epps : w_vec.ACOSEM_rhs;
				vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (uu / w_vec.ACOSEM_rhs);
			}
			yy += im_dim;
		}
	}
	// Anisotropic Diffusion smoothing prior
	if (MethodList.AD) {
		if (MethodList.OSLOSEM) {
			if (osa_iter == 0u)
				if (MethodList.OSEM)
					vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(0, im_dim - 1u));
				else
					vec.im_os(seq(yy, yy + im_dim - 1u)) = EM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)));
			else {
				array dU = AD(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType,
					w_vec.DiffusionType, w_vec.med_no_norm);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = EM(vec.im_os(seq(yy, yy + im_dim - 1u)), OSL(*testi, dU, beta.AD_OSEM, epps), vec.rhs_os(seq(yy, yy + im_dim - 1u)));
			}
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_BSREM, iter);
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			if (osa_iter == 0u) {
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)),
					w_vec.U, pj3, w_vec.lambda_MBSREM, iter, im_dim, 0.f, constant(0.f, 1, 1), *testi, epps);
			}
			else {
				array dU = AD(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType,
					w_vec.DiffusionType, w_vec.med_no_norm);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)),
					w_vec.U, pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.AD_MBSREM, dU, *testi, epps);
			}
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_ROSEM, iter);
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			if (osa_iter == 0u) {
				vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
					w_vec.D, 0.f, constant(0.f, 1, 1));
			}
			else {
				array dU = AD(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType,
					w_vec.DiffusionType, w_vec.med_no_norm);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
					w_vec.D, beta.AD_RBI, dU);
			}
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0u) {
			if (MethodList.OSLCOSEM == 1u)
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * pow(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.h_ACOSEM_2);
			else
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * vec.im_os(seq(yy, yy + im_dim - 1u));
			if (osa_iter == 0u) {
				vec.im_os(seq(yy, yy + im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, w_vec.D, w_vec.h_ACOSEM,
					MethodList.OSLCOSEM);
				if (MethodList.OSLCOSEM == 1u) {
					array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
					MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, w_vec, Summ,
						d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
						length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E, d_norm, d_scat, use_psf, g, Nx, Ny, Nz, epps, TOF, loadTOF, Sin, nBins, 
						koko, randoms_correction);
					//vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
					w_vec.ACOSEM_rhs = w_vec.ACOSEM_rhs < epps ? epps : w_vec.ACOSEM_rhs;
					vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (uu / w_vec.ACOSEM_rhs);
				}
			}
			else {
				array dU = AD(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType,
					w_vec.DiffusionType, w_vec.med_no_norm);
				vec.im_os(seq(yy, yy + im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, OSL(w_vec.D, dU, beta.AD_COSEM, epps), w_vec.h_ACOSEM,
					MethodList.OSLCOSEM);
				if (MethodList.OSLCOSEM == 1u) {
					array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
					MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, w_vec, Summ,
						d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
						length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E, d_norm, d_scat, use_psf, g, Nx, Ny, Nz, epps, TOF, loadTOF, Sin, nBins, 
						koko, randoms_correction);
					//vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
					w_vec.ACOSEM_rhs = w_vec.ACOSEM_rhs < epps ? epps : w_vec.ACOSEM_rhs;
					vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (uu / w_vec.ACOSEM_rhs);
				}
			}
			yy += im_dim;
		}
	}
	// Asymmetric Parallel Level Sets prior
	if (MethodList.APLS) {
		if (MethodList.OSLOSEM) {
			array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, 4, w_vec, w_vec.tr_offsets);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = EM(vec.im_os(seq(yy, yy + im_dim - 1u)), OSL(*testi, dU, beta.APLS_OSEM, epps), vec.rhs_os(seq(yy, yy + im_dim - 1u)));
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_BSREM, iter);
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, 4, w_vec, w_vec.tr_offsets);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U,
				pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.APLS_MBSREM, dU, *testi, epps);
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_ROSEM, iter);
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, 4, w_vec, w_vec.tr_offsets);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.D, beta.APLS_RBI, dU);
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0u) {
			array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, 4, w_vec, w_vec.tr_offsets);
			if (MethodList.OSLCOSEM == 1u)
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * pow(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.h_ACOSEM_2);
			else
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * vec.im_os(seq(yy, yy + im_dim - 1u));
			vec.im_os(seq(yy, yy + im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, OSL(w_vec.D, dU, beta.APLS_COSEM, epps), w_vec.h_ACOSEM,
				MethodList.OSLCOSEM);
			if (MethodList.OSLCOSEM == 1u) {
				array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
				MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, w_vec, Summ, d_Sino,
					koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length, atomic_64bit,
					compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E, d_norm, d_scat, use_psf, g, Nx, Ny, Nz, epps, TOF, loadTOF, Sin, nBins, koko, randoms_correction);
				//vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
				w_vec.ACOSEM_rhs = w_vec.ACOSEM_rhs < epps ? epps : w_vec.ACOSEM_rhs;
				vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (uu / w_vec.ACOSEM_rhs);
			}
			yy += im_dim;
		}
	}
	// Total Generalized Variation prior
	if (MethodList.TGV) {
		if (MethodList.OSLOSEM) {
			array dU = TGV(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = EM(vec.im_os(seq(yy, yy + im_dim - 1u)), OSL(*testi, dU, beta.TGV_OSEM, epps), vec.rhs_os(seq(yy, yy + im_dim - 1u)));
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_BSREM, iter);
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			array dU = TGV(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U,
				pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.TGV_MBSREM, dU, *testi, epps);
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_ROSEM, iter);
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			array dU = TGV(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.D, beta.TGV_RBI, dU);
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0u) {
			array dU = TGV(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
			if (MethodList.OSLCOSEM == 1u)
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * pow(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.h_ACOSEM_2);
			else
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * vec.im_os(seq(yy, yy + im_dim - 1u));
			vec.im_os(seq(yy, yy + im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, OSL(w_vec.D, dU, beta.TGV_COSEM, epps), w_vec.h_ACOSEM,
				MethodList.OSLCOSEM);
			if (MethodList.OSLCOSEM == 1u) {
				array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
				MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, w_vec, Summ,
					d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
					length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E, d_norm, d_scat, use_psf, g, Nx, Ny, Nz, epps, TOF, loadTOF, Sin, nBins, 
					koko, randoms_correction);
				//vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
				w_vec.ACOSEM_rhs = w_vec.ACOSEM_rhs < epps ? epps : w_vec.ACOSEM_rhs;
				vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (uu / w_vec.ACOSEM_rhs);
			}
			yy += im_dim;
		}
	}
	// Non-local means prior
	if (MethodList.NLM) {
		if (MethodList.OSLOSEM) {
			array dU = NLM(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec, epps, Nx, Ny, Nz, OpenCLStruct);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = EM(vec.im_os(seq(yy, yy + im_dim - 1u)), OSL(*testi, dU, beta.NLM_OSEM, epps), vec.rhs_os(seq(yy, yy + im_dim - 1u)));
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_BSREM, iter);
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			array dU = NLM(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec, epps, Nx, Ny, Nz, OpenCLStruct);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U,
				pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.NLM_MBSREM, dU, *testi, epps);
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_ROSEM, iter);
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			array dU = NLM(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec, epps, Nx, Ny, Nz, OpenCLStruct);
			vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.D, beta.NLM_RBI, dU);
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0u) {
			array dU = NLM(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec, epps, Nx, Ny, Nz, OpenCLStruct);
			if (MethodList.OSLCOSEM == 1u)
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * pow(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.h_ACOSEM_2);
			else
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * vec.im_os(seq(yy, yy + im_dim - 1u));
			vec.im_os(seq(yy, yy + im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, OSL(w_vec.D, dU, beta.NLM_COSEM, epps), w_vec.h_ACOSEM,
				MethodList.OSLCOSEM);
			if (MethodList.OSLCOSEM == 1u) {
				array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
				MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, w_vec, Summ,
					d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
					length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E, d_norm, d_scat, use_psf, g, Nx, Ny, Nz, epps, TOF, loadTOF, Sin, nBins, 
					koko, randoms_correction);
				//vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
				w_vec.ACOSEM_rhs = w_vec.ACOSEM_rhs < epps ? epps : w_vec.ACOSEM_rhs;
				vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (uu / w_vec.ACOSEM_rhs);
			}
			yy += im_dim;
		}
	}
	// Custom prior
	if (MethodList.CUSTOM) {
		if (MethodList.OSLOSEM) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = EM(vec.im_os(seq(yy, yy + im_dim - 1u)), OSL(*testi, w_vec.dU_OSEM, beta.custom_OSEM, epps), vec.rhs_os(seq(yy, yy + im_dim - 1u)));
			yy += im_dim;
		}
		if (MethodList.BSREM) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_BSREM, iter);
			yy += im_dim;
		}
		if (MethodList.MBSREM) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U,
				pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.custom_MBSREM, w_vec.dU_MBSREM, *testi, epps);
			yy += im_dim;
		}
		if (MethodList.ROSEMMAP) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.lambda_ROSEM, iter);
			yy += im_dim;
		}
		if (MethodList.RBIOSL) {
			vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
				w_vec.D, beta.custom_RBI, w_vec.dU_RBI);
			yy += im_dim;
		}
		if (MethodList.OSLCOSEM > 0u) {
			if (MethodList.OSLCOSEM == 1u)
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * pow(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.h_ACOSEM_2);
			else
				vec.C_osl(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u)) * vec.im_os(seq(yy, yy + im_dim - 1u));
			vec.im_os(seq(yy, yy + im_dim - 1u)) = COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, OSL(w_vec.D, w_vec.dU_COSEM, beta.custom_COSEM, epps), w_vec.h_ACOSEM,
				MethodList.OSLCOSEM);
			if (MethodList.OSLCOSEM == 1u) {
				array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
				MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, w_vec, Summ,
					d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
					length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E, d_norm, d_scat, use_psf, g, Nx, Ny, Nz, epps, TOF, loadTOF, Sin, nBins, 
					koko, randoms_correction);
				//vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
				w_vec.ACOSEM_rhs = w_vec.ACOSEM_rhs < epps ? epps : w_vec.ACOSEM_rhs;
				vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(yy, yy + im_dim - 1u)) * (uu / w_vec.ACOSEM_rhs);
			}
			yy += im_dim;
		}
		break_iter = true;
	}
}