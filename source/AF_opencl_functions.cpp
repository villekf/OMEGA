/**************************************************************************
* All the functions needed for the matrix-free OpenCL image reconstruction
*
* Copyright(C) 2019  Ville - Veikko Wettenhovi
*
* This program is free software : you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <https://www.gnu.org/licenses/>.
***************************************************************************/
#include "AF_opencl_functions.hpp"


// Update the OpenCL kernel inputs for the current iteration/subset
// If a method is not used, do nothing
// Otherwise create an OpenCL buffer pointer from the ArrayFire image estimate and initialize the right-hand side vector


//// Reconsruction methods as cl_chars
//void OpenCLRecMethods(const RecMethods &MethodList, RecMethodsOpenCL &MethodListOpenCL)
//{
//	// Non-MAP/prior algorithms
//	MethodListOpenCL.MLEM = static_cast<cl_char>(MethodList.MLEM);
//	MethodListOpenCL.OSEM = static_cast<cl_char>(MethodList.OSEM);
//	MethodListOpenCL.RAMLA = static_cast<cl_char>(MethodList.RAMLA);
//	MethodListOpenCL.MRAMLA = static_cast<cl_char>(MethodList.MRAMLA);
//	MethodListOpenCL.ROSEM = static_cast<cl_char>(MethodList.ROSEM);
//	MethodListOpenCL.RBI = static_cast<cl_char>(MethodList.RBI);
//	MethodListOpenCL.DRAMA = static_cast<cl_char>(MethodList.DRAMA);
//	MethodListOpenCL.COSEM = static_cast<cl_char>(MethodList.COSEM);
//	MethodListOpenCL.ECOSEM = static_cast<cl_char>(MethodList.ECOSEM);
//	MethodListOpenCL.ACOSEM = static_cast<cl_char>(MethodList.ACOSEM);
//	MethodListOpenCL.LSQR = static_cast<cl_char>(MethodList.LSQR);
//
//	// Priors
//	MethodListOpenCL.MRP = static_cast<cl_char>(MethodList.MRP);
//	MethodListOpenCL.Quad = static_cast<cl_char>(MethodList.Quad);
//	MethodListOpenCL.Huber = static_cast<cl_char>(MethodList.Huber);
//	MethodListOpenCL.L = static_cast<cl_char>(MethodList.L);
//	MethodListOpenCL.FMH = static_cast<cl_char>(MethodList.FMH);
//	MethodListOpenCL.WeightedMean = static_cast<cl_char>(MethodList.WeightedMean);
//	MethodListOpenCL.TV = static_cast<cl_char>(MethodList.TV);
//	MethodListOpenCL.AD = static_cast<cl_char>(MethodList.AD);
//	MethodListOpenCL.APLS = static_cast<cl_char>(MethodList.APLS);
//	MethodListOpenCL.TGV = static_cast<cl_char>(MethodList.TGV);
//	MethodListOpenCL.NLM = static_cast<cl_char>(MethodList.NLM);
//
//	// MAP/prior-based algorithms
//	MethodListOpenCL.OSLMLEM = static_cast<cl_char>(MethodList.OSLMLEM);
//	MethodListOpenCL.OSLOSEM = static_cast<cl_char>(MethodList.OSLOSEM);
//	MethodListOpenCL.BSREM = static_cast<cl_char>(MethodList.BSREM);
//	MethodListOpenCL.MBSREM = static_cast<cl_char>(MethodList.MBSREM);
//	MethodListOpenCL.ROSEMMAP = static_cast<cl_char>(MethodList.ROSEMMAP);
//	MethodListOpenCL.RBIOSL = static_cast<cl_char>(MethodList.RBIOSL);
//	MethodListOpenCL.OSLCOSEM = static_cast<cl_char>(MethodList.OSLCOSEM);
//	MethodListOpenCL.PKMA = static_cast<cl_char>(MethodList.PKMA);
//	MethodListOpenCL.CP = static_cast<cl_char>(MethodList.CP);
//}



// Prepass phase for MRAMLA, COSEM, ACOSEM, ECOSEM
//void MRAMLA_prepass(const uint32_t subsets, const uint32_t im_dim, const int64_t* pituus, const std::vector<cl::Buffer> &lor, const std::vector<cl::Buffer> &zindex,
//	const std::vector<cl::Buffer> &xindex, cl::Program program, const cl::CommandQueue &af_queue, const cl::Context af_context, Weighting& w_vec,
//	std::vector<af::array>& Summ, std::vector<cl::Buffer> &d_Sino, const size_t koko_l, af::array& cosem, af::array& C_co,
//	af::array& C_aco, af::array& C_osl, const uint32_t alku, cl::Kernel &kernel_mramla, const std::vector<cl::Buffer> &L, const uint8_t raw,
//	const RecMethodsOpenCL MethodListOpenCL, const std::vector<int64_t> length, const bool atomic_64bit, const bool atomic_32bit, const cl_uchar compute_norm_matrix,
//	const std::vector<cl::Buffer>& d_sc_ra, cl_uint kernelInd_MRAMLA, af::array& E, const std::vector<cl::Buffer>& d_norm, const std::vector<cl::Buffer>& d_scat, const bool use_psf,
//	const af::array& g, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float epps, const bool TOF, const bool loadTOF, const mxArray* Sin, const int64_t nBins, 
//	const size_t koko, const bool randoms_correction, const size_t local_size[], const uint64_t* randSize, const uint32_t Nt, const bool CT) {
//
//	cl_int status = CL_SUCCESS;
//
//	cl_uchar MBSREM_prepass = static_cast<cl_uchar>(w_vec.MBSREM_prepass);
//
//	af::array apu_co, apu_aco, uu, sub_index_array, apu_summa, apu_Amin, apu_summa_m;
//
//	const bool U_skip = w_vec.U > 0.f ? false : true;
//
//	af::array cosem_psf;
//	if (use_psf) {
//		cosem_psf = computeConvolution(cosem, g, Nx, Ny, Nz, w_vec, 1u);
//		af::sync();
//	}
//	else
//		cosem_psf = cosem;
//
//	cl::Buffer d_cosem = cl::Buffer(*cosem_psf.device<cl_mem>(), true);
//
//	af::array apu;
//
//	uint32_t eka = alku;
//
//	if (alku > 0u) {
//		eka = alku - 1u;
//		apu = af::constant(0.f, length[eka] * nBins);
//	}
//	else {
//		apu = af::constant(0.f, 1);
//	}
//	cl::Buffer d_ACOSEM_lhs = cl::Buffer(*apu.device<cl_mem>(), true);
//
//	cl_ulong st = 0ULL;
//
////	for (uint32_t osa_iter = eka; osa_iter < subsets; osa_iter++) {
////
////		if (DEBUG) {
////			mexPrintf("prepass kernel iteration started\n");
////			mexEvalString("pause(.0001);");
////		}
////
////		cl_uint kernelInd_MRAMLA_sub = kernelInd_MRAMLA;
////
////		size_t global_size = length[osa_iter];
////		if ((CT || inputScalars.SPECT || inputScalars.PET) && listmode == 0) {
////			erotus[0] = w_vec.size_x % local_size[0];
////			erotus[1] = w_vec.size_y % local_size[1];
////			if (erotus[1] > 0)
////				erotus[1] = (local_size[1] - erotus[1]);
////		}
////		else
////			erotus[0] = global_size % local_size[0];
////
////		if (erotus[0] > 0)
////			erotus[0] = (local_size[0] - erotus[0]);
////
////		cl::NDRange local(local_size[0]);
////		cl::NDRange global(global_size);
////
////		if ((CT || inputScalars.SPECT || inputScalars.PET) && listmode == 0) {
////			global = { w_vec.size_x + erotus[0], w_vec.size_y + erotus[1], length[osa_iter] };
////			local = { local_size[0] , local_size[1] };
////			m_size = static_cast<uint64_t>(w_vec.size_x) * static_cast<uint64_t>(w_vec.size_y) * length[osa_iter];
////		}
////		else {
////			global = { length[osa_iter] + erotus[0], 1, 1 };
////			m_size = static_cast<uint64_t>(length[osa_iter]);
////		}
////
////		size_t erotus = length[osa_iter] % local_size;
////
////		if (erotus > 0)
////			erotus = (local_size - erotus);
////
////		const size_t global_size = length[osa_iter] + erotus;
////
////		const uint64_t m_size = static_cast<uint64_t>(length[osa_iter]);
////
////		if (subsets > 1) {
////
////
////			if ((MethodListOpenCL.MRAMLA || MethodListOpenCL.MBSREM) && w_vec.MBSREM_prepass && alku == 0u) {
////				apu_Amin = af::constant(0.f, length[osa_iter] * nBins, 1);
////				apu_summa_m = af::constant(0.f, length[osa_iter] * nBins, 1);
////			}
////			else {
////				apu_Amin = af::constant(0.f, 1, 1);
////				apu_summa_m = af::constant(0.f, 1, 1);
////			}
////
////			if (w_vec.MBSREM_prepass && alku == 0u) {
////				if (atomic_64bit) {
////					apu_summa = af::constant(std::int64_t{0}, im_dim, 1, s64);
////				}
////				else if (atomic_32bit) {
////					apu_summa = af::constant(0, im_dim, 1, s32);
////				}
////				else {
////					apu_summa = af::constant(0.f, im_dim, 1);
////				}
////			}
////			else {
////				if (atomic_64bit) {
////					apu_summa = af::constant(std::int64_t{0}, 1, 1, s64);
////				}
////				else if (atomic_32bit) {
////					apu_summa = af::constant(0, 1, 1, s32);
////				}
////				else {
////					apu_summa = af::constant(0.f, 1, 1);
////				}
////			}
////
////			if ((MethodListOpenCL.COSEM || MethodListOpenCL.ECOSEM || MethodListOpenCL.OSLCOSEM == 2) && alku == 0u) {
////				if (atomic_64bit)
////					apu_co = af::constant(std::int64_t{0}, im_dim, 1, s64);
////				else if (atomic_32bit)
////					apu_co = af::constant(0, im_dim, 1, s32);
////				else
////					apu_co = af::constant(0.f, im_dim, 1);
////			}
////			else {
////				if (atomic_64bit)
////					apu_co = af::constant(std::int64_t{0}, 1, 1, s64);
////				else if (atomic_32bit)
////					apu_co = af::constant(0, 1, 1, s32);
////				else
////					apu_co = af::constant(0.f, 1, 1);
////			}
////
////			if ((MethodListOpenCL.ACOSEM || MethodListOpenCL.OSLCOSEM == 1) && alku == 0u) {
////				if (atomic_64bit)
////					apu_aco = af::constant(std::int64_t{0}, im_dim, 1, s64);
////				else if (atomic_32bit)
////					apu_aco = af::constant(0, im_dim, 1, s32);
////				else
////					apu_aco = af::constant(0.f, im_dim, 1);
////			}
////			else {
////				if (atomic_64bit)
////					apu_aco = af::constant(std::int64_t{0}, 1, 1, s64);
////				else if (atomic_32bit)
////					apu_aco = af::constant(0, 1, 1, s32);
////				else
////					apu_aco = af::constant(0.f, 1, 1);
////			}
////			if (TOF && !loadTOF && osa_iter > 0) {
////#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
////				float* apuS = (float*)mxGetSingles(mxGetCell(Sin, 0));
////#else
////				float* apuS = (float*)mxGetData(mxGetCell(Sin, 0));
////#endif
////				d_Sino[0] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[osa_iter] * nBins, NULL, &status);
////				for (int64_t to = std::int64_t{0}; to < nBins; to++)
////					status = af_queue.enqueueWriteBuffer(d_Sino[0], CL_FALSE, sizeof(float) * length[osa_iter] * to, sizeof(float) * length[osa_iter], &apuS[pituus[osa_iter] + koko * to]);
////			}
////		}
////		else {
////
////			if (w_vec.MBSREM_prepass && alku == 0u) {
////				if (atomic_64bit) {
////					apu_summa = af::constant(std::int64_t{0}, im_dim, 1, s64);
////				}
////				else if (atomic_32bit) {
////					apu_summa = af::constant(0, im_dim, 1, s32);
////				}
////				else {
////					apu_summa = Summ[0];
////				}
////			}
////			else {
////				if (atomic_64bit) {
////					apu_summa = af::constant(std::int64_t{0}, 1, 1, s64);
////				}
////				else if (atomic_32bit) {
////					apu_summa = af::constant(0, 1, 1, s32);
////				}
////				else {
////					apu_summa = af::constant(0.f, 1, 1);
////				}
////			}
////			apu_summa_m = E;
////
////			if ((MethodListOpenCL.COSEM || MethodListOpenCL.ECOSEM || MethodListOpenCL.OSLCOSEM == 2) && alku == 0u) {
////				if (atomic_64bit)
////					apu_co = (C_co * TH).as(s64);
////				else if (atomic_32bit)
////					apu_co = (C_co * TH32).as(s32);
////				else
////					apu_co = C_co;
////			}
////			else {
////				if (atomic_64bit)
////					apu_co = af::constant(std::int64_t{0}, 1, 1, s64);
////				else
////					apu_co = af::constant(0.f, 1, 1);
////			}
////
////			if ((MethodListOpenCL.ACOSEM || MethodListOpenCL.OSLCOSEM == 1) && alku == 0) {
////				if (atomic_64bit)
////					apu_aco = (C_aco * TH).as(s64);
////				else if (atomic_32bit)
////					apu_aco = (C_aco * TH32).as(s32);
////				else
////					apu_aco = C_aco;
////			}
////			else {
////				if (atomic_64bit)
////					apu_aco = af::constant(std::int64_t{0}, 1, 1, s64);
////				else if (atomic_32bit)
////					apu_aco = af::constant(0, 1, 1, s32);
////				else
////					apu_aco = af::constant(0.f, 1, 1);
////			}
////
////			if ((MethodListOpenCL.MRAMLA || MethodListOpenCL.MBSREM) && w_vec.MBSREM_prepass && alku == 0u)
////				apu_Amin = w_vec.Amin;
////			else
////				apu_Amin = af::constant(0.f, 1, 1);
////		}
////
////		cl::Buffer d_apu_co = cl::Buffer(*apu_co.device<cl_mem>(), true);
////		cl::Buffer d_apu_aco = cl::Buffer(*apu_aco.device<cl_mem>(), true);
////		cl::Buffer d_E = cl::Buffer(*apu_summa_m.device<cl_mem>(), true);
////		cl::Buffer d_Amin = cl::Buffer(*apu_Amin.device<cl_mem>(), true);
////		cl::Buffer d_Summ = cl::Buffer(*apu_summa.device<cl_mem>(), true);
////		af::sync();
////
////		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_norm[osa_iter]);
////		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_scat[osa_iter]);
////		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_Summ);
////		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, lor[osa_iter]);
////		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, xindex[osa_iter]);
////		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, zindex[osa_iter]);
////		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, L[osa_iter]);
////		if (TOF && !loadTOF)
////			kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_Sino[0]);
////		else
////			kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_Sino[osa_iter]);
////		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_sc_ra[osa_iter]);
////		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_cosem);
////		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, alku);
////		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, MBSREM_prepass);
////		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_ACOSEM_lhs);
////		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_Amin);
////		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_apu_co);
////		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_apu_aco);
////		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_E);
////		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, m_size);
////		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, MethodListOpenCL);
////		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, st);
////		cl::NDRange local(local_size);
////		cl::NDRange global(global_size);
////		// Compute the kernel
////		status = af_queue.enqueueNDRangeKernel(kernel_mramla, cl::NullRange, global, local);
////		if (status != CL_SUCCESS) {
////			getErrorString(status);
////			mexPrintf("Failed to launch the prepass kernel\n");
////			return;
////		}
////		else if (DEBUG) {
////			mexPrintf("prepass kernel launched successfully\n");
////			mexEvalString("pause(.0001);");
////		}
////
////		status = af_queue.finish();
////		if (status != CL_SUCCESS) {
////			getErrorString(status);
////			mexPrintf("Queue finish failed\n");
////			return;
////		}
////		apu_co.unlock();
////		apu_aco.unlock();
////		apu_summa.unlock();
////		apu_summa_m.unlock();
////		apu.unlock();
////		apu_Amin.unlock();
////		cosem_psf.unlock();
////		af::sync();
////		st += length[osa_iter];
////
////		if (alku == 0u) {
////			if ((MethodListOpenCL.COSEM || MethodListOpenCL.ECOSEM)) {
////				if (atomic_64bit)
////					C_co(af::span, osa_iter) = apu_co.as(f32) / TH;
////				else if (atomic_32bit)
////					C_co(af::span, osa_iter) = apu_co.as(f32) / TH32;
////				else
////					C_co(af::span, osa_iter) = apu_co;
////				if (use_psf)
////					C_co(af::span, osa_iter) = computeConvolution(C_co(af::span, osa_iter), g, Nx, Ny, Nz, w_vec, 1u) * cosem;
////				else
////					C_co(af::span, osa_iter) = C_co(af::span, osa_iter) * cosem;
////			}
////			if (MethodListOpenCL.ACOSEM) {
////				if (atomic_64bit)
////					C_aco(af::span, osa_iter) = apu_aco.as(f32) / TH;
////				else if (atomic_32bit)
////					C_aco(af::span, osa_iter) = apu_aco.as(f32) / TH32;
////				else
////					C_aco(af::span, osa_iter) = apu_aco;
////				if (use_psf)
////					C_aco(af::span, osa_iter) = computeConvolution(C_aco(af::span, osa_iter), g, Nx, Ny, Nz, w_vec, 1u) * af::pow(cosem, w_vec.h_ACOSEM_2);
////				else
////					C_aco(af::span, osa_iter) = C_aco(af::span, osa_iter) * af::pow(cosem, w_vec.h_ACOSEM_2);
////			}
////			if (MethodListOpenCL.OSLCOSEM == 2u) {
////				if (atomic_64bit)
////					C_osl(af::span, osa_iter) = apu_co.as(f32) / TH;
////				else if (atomic_32bit)
////					C_osl(af::span, osa_iter) = apu_co.as(f32) / TH32;
////				else
////					C_osl(af::span, osa_iter) = apu_co;
////				if (use_psf)
////					C_osl(af::span, osa_iter) = computeConvolution(C_osl(af::span, osa_iter), g, Nx, Ny, Nz, w_vec, 1u) * cosem;
////				else
////					C_osl(af::span, osa_iter) = C_osl(af::span, osa_iter) * cosem;
////			}
////			else if (MethodListOpenCL.OSLCOSEM == 1) {
////				if (atomic_64bit)
////					C_osl(af::span, osa_iter) = apu_aco.as(f32) / TH;
////				else if (atomic_32bit)
////					C_osl(af::span, osa_iter) = apu_aco.as(f32) / TH32;
////				else
////					C_osl(af::span, osa_iter) = apu_aco;
////				if (use_psf)
////					C_osl(af::span, osa_iter) = computeConvolution(C_osl(af::span, osa_iter), g, Nx, Ny, Nz, w_vec, 1u) * af::pow(cosem, w_vec.h_ACOSEM_2);
////				else
////					C_osl(af::span, osa_iter) = C_osl(af::span, osa_iter) * af::pow(cosem, w_vec.h_ACOSEM_2);
////			}
////			//if (DEBUG) {
////			//	mexPrintf("co = %f\n", af::sum<float>(C_aco(af::span, osa_iter)));
////			//	mexPrintf("dim0 = %u\n", C_aco(af::span, osa_iter).dims(0));
////			//	mexPrintf("dim1 = %u\n", C_aco(af::span, osa_iter).dims(1));
////			//}
////
////			if ((MethodListOpenCL.MRAMLA || MethodListOpenCL.MBSREM) && w_vec.MBSREM_prepass && Nt > 1U) {
////				sub_index_array = af::range(af::dim4(length[osa_iter] * nBins), 0, s64) + pituus[osa_iter] * nBins;
////				if (subsets > 1)
////					w_vec.Amin(sub_index_array) = apu_Amin;
////				else
////					w_vec.Amin = apu_Amin;
////			}
////
////			if (w_vec.MBSREM_prepass) {
////				if (DEBUG) {
////					mexPrintf("D = %f\n", af::sum<float>(w_vec.D));
////					mexEvalString("pause(.0001);");
////				}
////				if (compute_norm_matrix == 0u) {
////					if (atomic_64bit)
////						Summ[osa_iter] = (apu_summa).as(f32) / TH;
////					else if (atomic_32bit)
////						Summ[osa_iter] = (apu_summa).as(f32) / TH32;
////					else
////						Summ[osa_iter] = apu_summa;
////					af::sync();
////					w_vec.D += Summ[osa_iter];
////					Summ[osa_iter](Summ[osa_iter] < epps) = epps;
////					if (use_psf)
////						Summ[osa_iter] = computeConvolution(Summ[osa_iter], g, Nx, Ny, Nz, w_vec, 1u);
////					if (DEBUG) {
////						mexPrintf("Summ[osa_iter] = %f\n", af::sum<float>(Summ[osa_iter]));
////						mexEvalString("pause(.0001);");
////					}
////				}
////				else {
////					if (atomic_64bit)
////						Summ[0] = (apu_summa).as(f32) / TH;
////					else if (atomic_32bit)
////						Summ[0] = (apu_summa).as(f32) / TH32;
////					else
////						Summ[0] = apu_summa;
////					w_vec.D += Summ[0];
////					if (DEBUG) {
////						mexPrintf("Summ[0] = %f\n", af::sum<float>(Summ[0]));
////						mexPrintf("atomic_64bit = %d\n", atomic_64bit);
////						mexEvalString("pause(.0001);");
////					}
////				}
////
////				if ((MethodListOpenCL.MRAMLA || MethodListOpenCL.MBSREM) && Nt > 1U) {
////					if (subsets > 1)
////						E(sub_index_array) = apu_summa_m;
////					else
////						E = apu_summa_m;
////				}
////			}
////		}
////		else {
////			w_vec.ACOSEM_rhs = af::sum<float>(apu);
////		}
////		if (alku == 0u && (MethodListOpenCL.MBSREM || MethodListOpenCL.MRAMLA) && w_vec.MBSREM_prepass && Nt == 1U) {
////			uint32_t H = osa_iter;
////			uint32_t L = 0U;
////			if (TOF && !loadTOF)
////				H = 0;
////			if (randoms_correction)
////				L = osa_iter;
////			const af::array Sino = afcl::array(length[osa_iter] * nBins, d_Sino[H](), f32, true);
////			clRetainMemObject(d_Sino[H]());
////			const af::array rand = afcl::array(randSize[L], d_sc_ra[L](), f32, true);
////			clRetainMemObject(d_sc_ra[L]());
////			if (U_skip) {
////				float UU = w_vec.U;
////				const af::array Aind = apu_Amin > 0.f;
////				w_vec.U = af::max<float>(Sino(Aind) / apu_Amin(Aind));
////				if (UU > w_vec.U || std::isinf(w_vec.U))
////					w_vec.U = UU;
////			}
////			float eps_mramla = w_vec.epsilon_mramla;
////			w_vec.epsilon_mramla = MBSREM_epsilon(Sino, epps, randoms_correction, rand, apu_summa_m, TOF, nBins, CT);
////			if (eps_mramla < w_vec.epsilon_mramla)
////				w_vec.epsilon_mramla = eps_mramla;
////			if (DEBUG) {
////				mexPrintf("w_vec.epsilon_mramla = %f\n", w_vec.epsilon_mramla);
////				mexEvalString("pause(.0001);");
////			}
////		}
////		af::deviceGC();
////	}
//	if (TOF && !loadTOF) {
//#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
//		float* apu = (float*)mxGetSingles(mxGetCell(Sin, 0));
//#else
//		float* apu = (float*)mxGetData(mxGetCell(Sin, 0));
//#endif
//		d_Sino[0] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[0] * nBins, NULL, &status);
//		for (int64_t to = std::int64_t{0}; to < nBins; to++)
//			status = af_queue.enqueueWriteBuffer(d_Sino[0], CL_FALSE, sizeof(float) * length[0] * to, sizeof(float) * length[0], &apu[pituus[0] + koko * to]);
//	}
//	if (use_psf && alku == 0 && w_vec.MBSREM_prepass) {
//		w_vec.D = computeConvolution(w_vec.D, g, Nx, Ny, Nz, w_vec, 1u);
//	}
//	if (w_vec.MBSREM_prepass)
//		w_vec.D(w_vec.D <= 0.f) = 1.f;
//	return;
//}

//void precomp_siddon(const cl_context &context, const cl_command_queue &commandQueues,
//	uint16_t* lor1, const float* z_det, const float* x, const float* y, const uint32_t Nx, const uint32_t Ny,
//	const uint32_t Nz, const float dx, const float dy, const float dz, const float bx, const float by, const float bz,
//	const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices,
//	const uint32_t size_x, const uint16_t TotSinos, const bool verbose, const uint32_t loop_var_par,
//	const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L, const uint8_t raw, const size_t size_z, const uint32_t im_dim,
//	const cl_kernel &kernel, const size_t numel_x) {
//
//	cl_int status = CL_SUCCESS;
//	const uint32_t Nxy = Nx * Ny;
//
//	size_t osa_length = loop_var_par;
//
//	cl_mem d_z, d_x, d_y, d_pseudos, d_L, d_lor;
//
//
//	// Create the necessary buffers
//	d_z = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_z, NULL, &status);
//	if (status != CL_SUCCESS) {
//		getErrorString(status);
//		return;
//	}
//	d_x = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * numel_x, NULL, &status);
//	if (status != CL_SUCCESS) {
//		getErrorString(status);
//		return;
//	}
//	d_y = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * numel_x, NULL, &status);
//	if (status != CL_SUCCESS) {
//		getErrorString(status);
//		return;
//	}
//	d_pseudos = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t) * prows, NULL, &status);
//	if (status != CL_SUCCESS) {
//		getErrorString(status);
//		return;
//	}
//	d_lor = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(uint16_t) * osa_length, NULL, &status);
//	if (status != CL_SUCCESS) {
//		getErrorString(status);
//		return;
//	}
//	if (raw) {
//		d_L = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * osa_length * 2, NULL, &status);
//		if (status != CL_SUCCESS) {
//			getErrorString(status);
//			return;
//		}
//	}
//	else {
//		d_L = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
//		if (status != CL_SUCCESS) {
//			getErrorString(status);
//			return;
//		}
//	}
//
//	status = clEnqueueWriteBuffer(commandQueues, d_x, CL_FALSE, 0, sizeof(float) * numel_x, x, 0, NULL, NULL);
//	if (status != CL_SUCCESS) {
//		getErrorString(status);
//		return;
//	}
//	status = clEnqueueWriteBuffer(commandQueues, d_y, CL_FALSE, 0, sizeof(float) * numel_x, y, 0, NULL, NULL);
//	if (status != CL_SUCCESS) {
//		getErrorString(status);
//		return;
//	}
//	status = clEnqueueWriteBuffer(commandQueues, d_z, CL_FALSE, 0, sizeof(float) * size_z, z_det, 0, NULL, NULL);
//	if (status != CL_SUCCESS) {
//		getErrorString(status);
//		return;
//	}
//	status = clEnqueueWriteBuffer(commandQueues, d_pseudos, CL_FALSE, 0, sizeof(uint32_t) * prows, pseudos, 0, NULL, NULL);
//	if (status != CL_SUCCESS) {
//		getErrorString(status);
//		return;
//	}
//
//	if (raw) {
//		status = clEnqueueWriteBuffer(commandQueues, d_L, CL_FALSE, 0, sizeof(uint16_t) * osa_length * 2, L, 0, NULL, NULL);
//		if (status != CL_SUCCESS) {
//			getErrorString(status);
//			return;
//		}
//	}
//	else {
//		status = clEnqueueWriteBuffer(commandQueues, d_L, CL_FALSE, 0, sizeof(uint16_t), L, 0, NULL, NULL);
//		if (status != CL_SUCCESS) {
//			getErrorString(status);
//			return;
//		}
//	}
//
//
//	clSetKernelArg(kernel, 0, sizeof(uint32_t), &Nxy);
//	clSetKernelArg(kernel, 1, sizeof(uint32_t), &im_dim);
//	clSetKernelArg(kernel, 2, sizeof(uint32_t), &Nx);
//	clSetKernelArg(kernel, 3, sizeof(uint32_t), &Ny);
//	clSetKernelArg(kernel, 4, sizeof(uint32_t), &Nz);
//	clSetKernelArg(kernel, 5, sizeof(float), &dz);
//	clSetKernelArg(kernel, 6, sizeof(float), &dx);
//	clSetKernelArg(kernel, 7, sizeof(float), &dy);
//	clSetKernelArg(kernel, 8, sizeof(float), &bz);
//	clSetKernelArg(kernel, 9, sizeof(float), &bx);
//	clSetKernelArg(kernel, 10, sizeof(float), &by);
//	clSetKernelArg(kernel, 11, sizeof(float), &bzb);
//	clSetKernelArg(kernel, 12, sizeof(float), &maxxx);
//	clSetKernelArg(kernel, 13, sizeof(float), &maxyy);
//	clSetKernelArg(kernel, 14, sizeof(float), &zmax);
//	clSetKernelArg(kernel, 15, sizeof(float), &NSlices);
//	clSetKernelArg(kernel, 16, sizeof(uint32_t), &size_x);
//	clSetKernelArg(kernel, 17, sizeof(uint16_t), &TotSinos);
//	clSetKernelArg(kernel, 18, sizeof(uint32_t), &det_per_ring);
//	clSetKernelArg(kernel, 19, sizeof(uint8_t), &raw);
//	clSetKernelArg(kernel, 20, sizeof(uint32_t), &prows);
//
//	cl_event event1;
//
//	clSetKernelArg(kernel, 21, sizeof(cl_mem), &d_pseudos);
//	clSetKernelArg(kernel, 22, sizeof(cl_mem), &d_x);
//	clSetKernelArg(kernel, 23, sizeof(cl_mem), &d_y);
//	clSetKernelArg(kernel, 24, sizeof(cl_mem), &d_z);
//	clSetKernelArg(kernel, 25, sizeof(cl_mem), &d_lor);
//	clSetKernelArg(kernel, 26, sizeof(cl_mem), &d_L);
//	status = clEnqueueNDRangeKernel(commandQueues, kernel, 1, NULL, &osa_length, NULL, 0, NULL, &event1);
//	if (status != CL_SUCCESS) {
//		getErrorString(status);
//		return;
//	}
//
//	status = clEnqueueReadBuffer(commandQueues, d_lor, CL_TRUE, 0, sizeof(uint16_t) * osa_length, lor1, 1, &event1, NULL);
//	if (status != CL_SUCCESS) {
//		getErrorString(status);
//		return;
//	}
//
//	af::sync();
//	return;
//}
