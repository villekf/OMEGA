/**************************************************************************
* All the functions needed for the matrix-free CUDA image reconstruction
*
* Copyright(C) 2020 Ville - Veikko Wettenhovi
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
#include "AF_cuda_functions.hpp"

// Update the OpenCL kernel inputs for the current iteration/subset
// If a method is not used, do nothing
// Otherwise create an OpenCL buffer pointer from the ArrayFire image estimate and initialize the right-hand side vector
void update_cuda_inputs(AF_im_vectors & vec, CUDA_im_vectors &vec_cuda, const bool mlem, const uint32_t im_dim, const uint32_t n_rekos, 
	const uint32_t n_rekos_mlem, const RecMethods MethodList, const bool atomic_64bit, const bool use_psf)
{
	//if (MethodList.CUSTOM) {
	//	if (mlem) {
	//		vec.im_mlem(af::seq(n_rekos * im_dim - im_dim, n_rekos * im_dim - 1u)) = vec.custom_MLEM;
	//	}
	//	else {
	//		uint32_t yy = n_rekos * im_dim;
	//		if (MethodList.OSLCOSEM > 0u) {
	//			vec.im_os(af::seq(yy - im_dim, yy - 1u)) = vec.custom_COSEM;
	//			yy -= im_dim;
	//		}
	//		if (MethodList.RBIOSL) {
	//			vec.im_os(af::seq(yy - im_dim, yy - 1u)) = vec.custom_RBI;
	//			yy -= im_dim;
	//		}
	//		if (MethodList.ROSEMMAP) {
	//			vec.im_os(af::seq(yy - im_dim, yy - 1u)) = vec.custom_ROSEM;
	//			yy -= im_dim;
	//		}
	//		if (MethodList.MBSREM) {
	//			vec.im_os(af::seq(yy - im_dim, yy - 1u)) = vec.custom_MBSREM;
	//			yy -= im_dim;
	//		}
	//		if (MethodList.BSREM) {
	//			vec.im_os(af::seq(yy - im_dim, yy - 1u)) = vec.custom_BSREM;
	//			yy -= im_dim;
	//		}
	//		if (MethodList.OSLOSEM) {
	//			vec.im_os(af::seq(yy - im_dim, yy - 1u)) = vec.custom_OSEM;
	//			yy -= im_dim;
	//		}
	//	}
	//}
	if (mlem) {
		if (use_psf)
			vec_cuda.d_im_mlem = vec.im_mlem_blurred.device<CUdeviceptr>();
		else
			vec_cuda.d_im_mlem = vec.im_mlem.device<CUdeviceptr>();
		if (atomic_64bit)
			vec.rhs_mlem = af::constant(0ULL, im_dim * n_rekos_mlem, 1, u64);
		else
			vec.rhs_mlem = af::constant(0.f, im_dim * n_rekos_mlem, 1);
		vec_cuda.d_rhs_mlem = vec.rhs_mlem.device<CUdeviceptr>();
	}
	else {
		if (use_psf)
			vec_cuda.d_im_os = vec.im_os_blurred.device<CUdeviceptr>();
		else
			vec_cuda.d_im_os = vec.im_os.device<CUdeviceptr>();
		if (atomic_64bit)
			vec.rhs_os = af::constant(0ULL, im_dim * n_rekos, 1, u64);
		else
			vec.rhs_os = af::constant(0.f, im_dim * n_rekos, 1);
		vec_cuda.d_rhs_os = vec.rhs_os.device<CUdeviceptr>();
	}
}

// Reconsruction methods as cl_chars
void OpenCLRecMethods(const RecMethods &MethodList, RecMethodsOpenCL &MethodListOpenCL)
{
	MethodListOpenCL.MLEM = static_cast<char>(MethodList.MLEM);
	MethodListOpenCL.OSEM = static_cast<char>(MethodList.OSEM);
	MethodListOpenCL.RAMLA = static_cast<char>(MethodList.RAMLA);
	MethodListOpenCL.MRAMLA = static_cast<char>(MethodList.MRAMLA);
	MethodListOpenCL.ROSEM = static_cast<char>(MethodList.ROSEM);
	MethodListOpenCL.RBI = static_cast<char>(MethodList.RBI);
	MethodListOpenCL.DRAMA = static_cast<char>(MethodList.DRAMA);
	MethodListOpenCL.COSEM = static_cast<char>(MethodList.COSEM);
	MethodListOpenCL.ECOSEM = static_cast<char>(MethodList.ECOSEM);
	MethodListOpenCL.ACOSEM = static_cast<char>(MethodList.ACOSEM);

	MethodListOpenCL.MRP = static_cast<char>(MethodList.MRP);
	MethodListOpenCL.Quad = static_cast<char>(MethodList.Quad);
	MethodListOpenCL.L = static_cast<char>(MethodList.L);
	MethodListOpenCL.FMH = static_cast<char>(MethodList.FMH);
	MethodListOpenCL.WeightedMean = static_cast<char>(MethodList.WeightedMean);
	MethodListOpenCL.TV = static_cast<char>(MethodList.TV);
	MethodListOpenCL.AD = static_cast<char>(MethodList.AD);
	MethodListOpenCL.APLS = static_cast<char>(MethodList.APLS);
	MethodListOpenCL.TGV = static_cast<char>(MethodList.TGV);
	MethodListOpenCL.NLM = static_cast<char>(MethodList.NLM);

	MethodListOpenCL.OSLMLEM = static_cast<char>(MethodList.OSLMLEM);
	MethodListOpenCL.OSLOSEM = static_cast<char>(MethodList.OSLOSEM);
	MethodListOpenCL.BSREM = static_cast<char>(MethodList.BSREM);
	MethodListOpenCL.MBSREM = static_cast<char>(MethodList.MBSREM);
	MethodListOpenCL.ROSEMMAP = static_cast<char>(MethodList.ROSEMMAP);
	MethodListOpenCL.RBIOSL = static_cast<char>(MethodList.RBIOSL);
	MethodListOpenCL.OSLCOSEM = static_cast<char>(MethodList.OSLCOSEM);
}


CUresult createAndWriteBuffers(CUdeviceptr& d_x, CUdeviceptr& d_y, CUdeviceptr& d_z, CUdeviceptr& d_angles, std::vector<CUdeviceptr>& d_lor, std::vector<CUdeviceptr>& d_L,
	std::vector<CUdeviceptr>& d_zindex, std::vector<CUdeviceptr>& d_xyindex, std::vector<CUdeviceptr>& d_Sino, std::vector<CUdeviceptr>& d_sc_ra, 
	const uint32_t size_x, const size_t size_z, const uint32_t TotSinos, const size_t size_atten, const size_t size_norm, const size_t size_scat, const uint32_t prows,
	std::vector<size_t>& length, const float* x, const float* y, const float* z_det, const uint32_t* xy_index, const uint16_t* z_index, const uint16_t* lor1, 
	const uint16_t* L, const float* Sin, const uint8_t raw, const uint32_t subsets, const int64_t* pituus, const float* atten, const float* norm, const float* scat, 
	const uint32_t* pseudos, const float* V, CUdeviceptr& d_atten, std::vector<CUdeviceptr>& d_norm, std::vector<CUdeviceptr>& d_scat, CUdeviceptr& d_pseudos, CUdeviceptr& d_V,
	CUdeviceptr& d_xcenter, CUdeviceptr& d_ycenter, CUdeviceptr& d_zcenter, const float* x_center, const float* y_center, const float* z_center, const size_t size_center_x, 
	const size_t size_center_y, const size_t size_center_z, const size_t size_of_x, const size_t size_V, const bool randoms_correction, const mxArray* sc_ra, const bool precompute, 
	CUdeviceptr& d_lor_mlem, CUdeviceptr& d_L_mlem, CUdeviceptr& d_zindex_mlem, CUdeviceptr& d_xyindex_mlem, CUdeviceptr& d_Sino_mlem, CUdeviceptr& d_sc_ra_mlem, 
	CUdeviceptr& d_reko_type, CUdeviceptr& d_reko_type_mlem, const bool osem_bool, const bool mlem_bool, const size_t koko, const uint8_t* reko_type, const uint8_t* reko_type_mlem, 
	const uint32_t n_rekos, const uint32_t n_rekos_mlem, CUdeviceptr& d_norm_mlem, CUdeviceptr& d_scat_mlem, const float* angles, const bool TOF, const int64_t nBins, const bool loadTOF,
	CUdeviceptr& d_TOFCenter, const float* TOFCenter, const uint32_t subsetsUsed, const uint32_t osa_iter0, const uint8_t listmode, const bool CT)
{
	// Create the necessary buffers
	// Detector coordinates
	CUresult status = CUDA_SUCCESS;
	status = cuMemAlloc(&d_V, sizeof(float) * size_V);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = cuMemAlloc(&d_z, sizeof(float) * size_z);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = cuMemAlloc(&d_x, sizeof(float) * size_of_x);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = cuMemAlloc(&d_y, sizeof(float) * size_of_x);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	if (CT) {
		status = cuMemAlloc(&d_angles, sizeof(float) * TotSinos);
		if (status != CUDA_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
	}
	status = cuMemAlloc(&d_xcenter, sizeof(float) * size_center_x);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = cuMemAlloc(&d_ycenter, sizeof(float) * size_center_y);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = cuMemAlloc(&d_zcenter, sizeof(float) * size_center_z);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	// Attenuation data for image-based attenuation
	status = cuMemAlloc(&d_atten, sizeof(float) * size_atten);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	// TOF bin centers
	status = cuMemAlloc(&d_TOFCenter, sizeof(float) * nBins);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	// Pseudo rings
	status = cuMemAlloc(&d_pseudos, sizeof(uint32_t) * prows);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	if (osem_bool) {
		status = cuMemAlloc(&d_reko_type, sizeof(uint8_t) * n_rekos);
		if (status != CUDA_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		for (uint32_t kk = osa_iter0; kk < subsets; kk++) {
			// How many voxels does each LOR traverse
			if (precompute)
				status = cuMemAlloc(&d_lor[kk], sizeof(uint16_t) * length[kk]);
			else
				status = cuMemAlloc(&d_lor[kk], sizeof(uint16_t));
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			if (size_norm > 1)
				status = cuMemAlloc(&d_norm[kk], sizeof(float) * length[kk]);
			else
				status = cuMemAlloc(&d_norm[kk], sizeof(float) * size_norm);
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			if (size_scat > 1)
				status = cuMemAlloc(&d_scat[kk], sizeof(float) * length[kk]);
			else
				status = cuMemAlloc(&d_scat[kk], sizeof(float) * size_scat);
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			// Measurement data
			if (TOF) {
				if (loadTOF) {
					status = cuMemAlloc(&d_Sino[kk], sizeof(float) * length[kk] * nBins);
				}
				else {
					if (kk == 0)
						status = cuMemAlloc(&d_Sino[kk], sizeof(float) * length[kk] * nBins);
				}
			}
			else
				status = cuMemAlloc(&d_Sino[kk], sizeof(float) * length[kk]);
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			if (randoms_correction == 1u)
				status = cuMemAlloc(&d_sc_ra[kk], sizeof(float) * length[kk]);
			else
				status = cuMemAlloc(&d_sc_ra[kk], sizeof(float));
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			// Indices corresponding to the detector index (Sinogram data) or the detector number (raw data) at each measurement
			if (raw && listmode != 1) {
				status = cuMemAlloc(&d_xyindex[kk], sizeof(uint32_t));
				if (status != CUDA_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				status = cuMemAlloc(&d_zindex[kk], sizeof(uint16_t));
				if (status != CUDA_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				status = cuMemAlloc(&d_L[kk], sizeof(uint16_t) * length[kk] * 2);
				if (status != CUDA_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
			}
			else if (listmode != 1 && (!CT || subsets > 1)) {
				status = cuMemAlloc(&d_xyindex[kk], sizeof(uint32_t) * length[kk]);
				if (status != CUDA_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				status = cuMemAlloc(&d_zindex[kk], sizeof(uint16_t) * length[kk]);
				if (status != CUDA_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				status = cuMemAlloc(&d_L[kk], sizeof(uint16_t));
				if (status != CUDA_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
			}
			else {
				status = cuMemAlloc(&d_xyindex[kk], sizeof(uint32_t));
				if (status != CUDA_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				status = cuMemAlloc(&d_zindex[kk], sizeof(uint16_t));
				if (status != CUDA_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				status = cuMemAlloc(&d_L[kk], sizeof(uint16_t));
				if (status != CUDA_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
			}
		}
	}
	if (mlem_bool) {
		status = cuMemAlloc(&d_reko_type_mlem, sizeof(uint8_t) * n_rekos_mlem);
		if (status != CUDA_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		if (precompute)
			status = cuMemAlloc(&d_lor_mlem, sizeof(uint16_t) * koko);
		else
			status = cuMemAlloc(&d_lor_mlem, sizeof(uint16_t));
		if (status != CUDA_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		// Measurement data
		if (TOF && listmode != 2)
			status = cuMemAlloc(&d_Sino_mlem, sizeof(float) * koko * nBins);
		else if (listmode != 2)
			status = cuMemAlloc(&d_Sino_mlem, sizeof(float) * koko);
		else
			status = cuMemAlloc(&d_Sino_mlem, sizeof(float));
		if (status != CUDA_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		status = cuMemAlloc(&d_norm_mlem, sizeof(float) * size_norm);
		if (status != CUDA_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		if (randoms_correction == 1u)
			status = cuMemAlloc(&d_sc_ra_mlem, sizeof(float) * koko);
		else
			status = cuMemAlloc(&d_sc_ra_mlem, sizeof(float));
		if (status != CUDA_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		status = cuMemAlloc(&d_scat_mlem, sizeof(float) * size_scat);
		if (status != CUDA_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		// Indices corresponding to the detector index (Sinogram data) or the detector number (raw data) at each measurement
		if (raw && listmode != 1) {
			status = cuMemAlloc(&d_xyindex_mlem, sizeof(uint32_t));
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			status = cuMemAlloc(&d_zindex_mlem, sizeof(uint16_t));
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			status = cuMemAlloc(&d_L_mlem, sizeof(uint16_t) * koko * 2);
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
		}
		else if (listmode != 1 && !CT) {
			status = cuMemAlloc(&d_xyindex_mlem, sizeof(uint32_t) * koko);
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			status = cuMemAlloc(&d_zindex_mlem, sizeof(uint16_t) * koko);
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			status = cuMemAlloc(&d_L_mlem, sizeof(uint16_t));
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
		}
		else {
			status = cuMemAlloc(&d_xyindex_mlem, sizeof(uint32_t));
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			status = cuMemAlloc(&d_zindex_mlem, sizeof(uint16_t));
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			status = cuMemAlloc(&d_L_mlem, sizeof(uint16_t));
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
		}
	}

	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		mexPrintf("Buffer creation failed\n");
		return status;
	}
	else if (DEBUG) {
		mexPrintf("Buffer creation succeeded\n");
		mexEvalString("pause(.0001);");
	}


	// assign values to the buffers
	status = cuMemcpyHtoD(d_V, V, sizeof(float) * size_V);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = cuMemcpyHtoD(d_x, x, sizeof(float) * size_of_x);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = cuMemcpyHtoD(d_y, y, sizeof(float) * size_of_x);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = cuMemcpyHtoD(d_z, z_det, sizeof(float) * size_z);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	if (CT) {
		status = cuMemcpyHtoD(d_angles, angles, sizeof(float) * TotSinos);
		if (status != CUDA_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
	}
	status = cuMemcpyHtoD(d_xcenter, x_center, sizeof(float) * size_center_x);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = cuMemcpyHtoD(d_ycenter, y_center, sizeof(float) * size_center_y);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = cuMemcpyHtoD(d_zcenter, z_center, sizeof(float) * size_center_z);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = cuMemcpyHtoD(d_atten, atten, sizeof(float) * size_atten);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = cuMemcpyHtoD(d_TOFCenter, TOFCenter, sizeof(float) * nBins);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = cuMemcpyHtoD(d_pseudos, pseudos, sizeof(uint32_t) * prows);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	if (osem_bool) {
		status = cuMemcpyHtoD(d_reko_type, reko_type, sizeof(uint8_t) * n_rekos);
		if (status != CUDA_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		for (uint32_t kk = osa_iter0; kk < subsets; kk++) {
			if (raw && listmode != 1) {
				status = cuMemcpyHtoD(d_xyindex[kk], xy_index, sizeof(uint32_t));
				if (status != CUDA_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				status = cuMemcpyHtoD(d_zindex[kk], z_index, sizeof(uint16_t));
				if (status != CUDA_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				status = cuMemcpyHtoD(d_L[kk], &L[pituus[kk] * 2], sizeof(uint16_t) * length[kk] * 2);
				if (status != CUDA_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
			}
			else if (listmode != 1 && (!CT || subsets > 1)) {
				status = cuMemcpyHtoD(d_xyindex[kk], &xy_index[pituus[kk]], sizeof(uint32_t) * length[kk]);
				if (status != CUDA_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				status = cuMemcpyHtoD(d_zindex[kk], &z_index[pituus[kk]], sizeof(uint16_t) * length[kk]);
				if (status != CUDA_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				status = cuMemcpyHtoD(d_L[kk], L, sizeof(uint16_t));
				if (status != CUDA_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
			}
			else {
				status = cuMemcpyHtoD(d_xyindex[kk], xy_index, sizeof(uint32_t));
				if (status != CUDA_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				status = cuMemcpyHtoD(d_zindex[kk], z_index, sizeof(uint16_t));
				if (status != CUDA_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				status = cuMemcpyHtoD(d_L[kk], L, sizeof(uint16_t));
				if (status != CUDA_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
			}
			if (size_norm > 1)
				status = cuMemcpyHtoD(d_norm[kk], &norm[pituus[kk]], sizeof(float) * length[kk]);
			else
				status = cuMemcpyHtoD(d_norm[kk], norm, sizeof(float) * size_norm);
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			if (size_scat > 1)
				status = cuMemcpyHtoD(d_scat[kk], &scat[pituus[kk]], sizeof(float) * length[kk]);
			else
				status = cuMemcpyHtoD(d_scat[kk], scat, sizeof(float) * size_scat);
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			if (precompute)
				status = cuMemcpyHtoD(d_lor[kk], &lor1[pituus[kk]], sizeof(uint16_t) * length[kk]);
			else
				status = cuMemcpyHtoD(d_lor[kk], lor1, sizeof(uint16_t));
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}

			if (TOF) {
				if (loadTOF) {
					for (int64_t to = 0LL; to < nBins; to++) {
						CUdeviceptr dst = reinterpret_cast<CUdeviceptr> (reinterpret_cast<char*>(d_Sino[kk]) + sizeof(float) * length[kk] * to);
						char* src = (char*)Sin + sizeof(float) * (pituus[kk] + koko * to);
						size_t bytes = sizeof(float) * length[kk];
						status = cuMemcpyHtoD(dst, src, bytes);
					}
				}
				else {
					if (kk == osa_iter0) {
						for (int64_t to = 0LL; to < nBins; to++) {
							CUdeviceptr dst = reinterpret_cast<CUdeviceptr> (reinterpret_cast<char*>(d_Sino[kk]) + sizeof(float) * length[kk] * to);
							char* src = (char*)Sin + sizeof(float) * (pituus[kk] + koko * to);
							size_t bytes = sizeof(float) * length[kk];
							status = cuMemcpyHtoD(dst, src, bytes);
						}
					}
				}
			}
			else if (listmode != 2)
				status = cuMemcpyHtoD(d_Sino[kk], &Sin[pituus[kk]], sizeof(float) * length[kk]);
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			float* apu = (float*)mxGetSingles(mxGetCell(sc_ra, 0));
#else
			float* apu = (float*)mxGetData(mxGetCell(sc_ra, 0));
#endif
			if (randoms_correction)
				status = cuMemcpyHtoD(d_sc_ra[kk], &apu[pituus[kk]], sizeof(float) * length[kk]);
			else
				status = cuMemcpyHtoD(d_sc_ra[kk], apu, sizeof(float));
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
		}
	}
	if (mlem_bool) {
		status = cuCtxSynchronize();
		status = cuMemcpyHtoD(d_reko_type_mlem, reko_type_mlem, sizeof(uint8_t) * n_rekos_mlem);
		if (status != CUDA_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		if (precompute)
			status = cuMemcpyHtoD(d_lor_mlem, lor1, sizeof(uint16_t) * koko);
		else
			status = cuMemcpyHtoD(d_lor_mlem, lor1, sizeof(uint16_t));
		if (status != CUDA_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}

		if (TOF)
			status = cuMemcpyHtoD(d_Sino_mlem, Sin, sizeof(float) * koko * nBins);
		else if (listmode != 2)
			status = cuMemcpyHtoD(d_Sino_mlem, Sin, sizeof(float) * koko);
		if (status != CUDA_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		status = cuMemcpyHtoD(d_norm_mlem, norm, sizeof(float) * size_norm);
		if (status != CUDA_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		status = cuMemcpyHtoD(d_scat_mlem, norm, sizeof(float) * size_scat);
		if (status != CUDA_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		float* apu = (float*)mxGetSingles(mxGetCell(sc_ra, 0));
#else
		float* apu = (float*)mxGetData(mxGetCell(sc_ra, 0));
#endif
		if (randoms_correction)
			status = cuMemcpyHtoD(d_sc_ra_mlem, apu, sizeof(float) * koko);
		else
			status = cuMemcpyHtoD(d_sc_ra_mlem, apu, sizeof(float));
		if (status != CUDA_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		if (raw && listmode != 1) {
			status = cuMemcpyHtoD(d_xyindex_mlem, xy_index, sizeof(uint32_t));
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			status = cuMemcpyHtoD(d_zindex_mlem, z_index, sizeof(uint16_t));
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			status = cuMemcpyHtoD(d_L_mlem, L, sizeof(uint16_t) * koko * 2ULL);
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
		}
		else if (listmode != 1 && !CT) {
			status = cuMemcpyHtoD(d_zindex_mlem, z_index, sizeof(uint16_t) * koko);
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			status = cuMemcpyHtoD(d_xyindex_mlem, xy_index, sizeof(uint32_t) * koko);
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			status = cuMemcpyHtoD(d_L_mlem, L, sizeof(uint16_t));
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
		}
		else {
			status = cuMemcpyHtoD(d_xyindex_mlem, xy_index, sizeof(uint32_t));
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			status = cuMemcpyHtoD(d_zindex_mlem, z_index, sizeof(uint16_t));
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			status = cuMemcpyHtoD(d_L_mlem, L, sizeof(uint16_t));
			if (status != CUDA_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
		}
	}
	status = cuCtxSynchronize();
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		mexPrintf("Buffer write failed\n");
		return status;
	}
	else if (DEBUG) {
		mexPrintf("Buffer write succeeded\n");
		mexEvalString("pause(.0001);");
	}
	return status;
}

 //Prepass phase for MRAMLA, COSEM, ACOSEM, ECOSEM
void MRAMLA_prepass_CUDA(const uint32_t subsets, const uint32_t im_dim, const int64_t* pituus, std::vector<CUdeviceptr>& d_lor, std::vector<CUdeviceptr>& d_zindex,
	std::vector<CUdeviceptr>& d_xyindex, Weighting& w_vec, std::vector<af::array>& Summ, std::vector<CUdeviceptr>& d_Sino, size_t koko_l, af::array& cosem,
	af::array& C_co, af::array& C_aco, af::array& C_osl, uint32_t alku, std::vector<CUdeviceptr>& d_L, uint8_t raw, RecMethodsOpenCL& MethodList,
	std::vector<size_t> length, uint8_t compute_norm_matrix, std::vector<CUdeviceptr>& d_sc_ra, af::array& E, const uint32_t det_per_ring, CUdeviceptr& d_pseudos,
	const uint32_t prows, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float dz, const float dx, const float dy, const float bz, const float bx,
	const float by, const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices, CUdeviceptr& d_x, CUdeviceptr& d_y, CUdeviceptr& d_z, 
	const uint32_t size_x, const uint32_t TotSinos, CUdeviceptr& d_atten, std::vector<CUdeviceptr>& d_norm, std::vector<CUdeviceptr>& d_scat, const float epps, const uint32_t Nxy, 
	const float tube_width, const float crystal_size_z, const float bmin, const float bmax, const float Vmax, CUdeviceptr& d_xcenter, CUdeviceptr& d_ycenter, 
	CUdeviceptr& d_zcenter, CUdeviceptr& d_V, const float dc_z, const uint16_t n_rays, const uint16_t n_rays3D, const bool precompute, const uint32_t projector_type, 
	const CUstream& af_cuda_stream, const float global_factor, CUdeviceptr& d_reko_type, CUfunction& kernel_mbsrem, const bool atomic_64bit, const bool use_psf, 
	const af::array& g, const bool TOF, const bool loadTOF, const mxArray* Sin, const int64_t nBins, const bool randoms_correction, const float sigma_x, CUdeviceptr& d_TOFCenter, 
	CUdeviceptr& d_angles, const uint32_t Nt, const bool CT) {

	CUresult status = CUDA_SUCCESS;

	unsigned char MBSREM_prepass = static_cast<unsigned char>(w_vec.MBSREM_prepass);

	af::array apu_co, apu_aco, uu, sub_index_array, apu_summa, apu_Amin, apu_summa_m;

	const bool U_skip = w_vec.U > 0.f ? false : true;

	af::array cosem_psf;
	if (use_psf) {
		cosem_psf = computeConvolution(cosem, g, Nx, Ny, Nz, w_vec, 1u);
		af::sync();
	}
	else
		cosem_psf = cosem;

	CUdeviceptr* d_cosem = cosem_psf.device<CUdeviceptr>();

	unsigned char fp = 0;

	af::array apu;

	uint32_t eka = alku;

	if (alku > 0u) {
		eka = alku - 1u;
		apu = af::constant(0.f, length[eka] * nBins);
	}
	else {
		apu = af::constant(0.f, 1);
	}
	CUdeviceptr* d_ACOSEM_lhs = apu.device<CUdeviceptr>();

	const size_t local_size = 64ULL;
	uint64_t st = 0ULL;

	for (uint32_t osa_iter = eka; osa_iter < subsets; osa_iter++) {

		size_t erotus = length[osa_iter] % local_size;

		if (erotus > 0)
			erotus = (local_size - erotus);

		size_t global_size = length[osa_iter] + erotus;

		global_size = global_size / local_size;

		const uint64_t m_size = static_cast<uint64_t>(length[osa_iter]);

		if (subsets > 1) {


			if ((MethodList.MRAMLA || MethodList.MBSREM) && w_vec.MBSREM_prepass && alku == 0u) {
				apu_Amin = af::constant(0.f, length[osa_iter] * nBins, 1);
				apu_summa_m = af::constant(0.f, length[osa_iter] * nBins, 1);
			}
			else {
				apu_Amin = af::constant(0.f, 1, 1);
				apu_summa_m = af::constant(0.f, 1, 1);
			}

			if (w_vec.MBSREM_prepass && alku == 0u) {
				if (atomic_64bit) {
					apu_summa = af::constant(0ULL, im_dim, 1, u64);
				}
				else {
					apu_summa = af::constant(0.f, im_dim, 1);
				}
			}
			else {
				if (atomic_64bit) {
					apu_summa = af::constant(0ULL, 1, 1, u64);
				}
				else {
					apu_summa = af::constant(0.f, 1, 1);
				}
			}

			if ((MethodList.COSEM || MethodList.ECOSEM || MethodList.OSLCOSEM == 2) && alku == 0u) {
				if (atomic_64bit)
					apu_co = af::constant(0ULL, im_dim, 1, u64);
				else
					apu_co = af::constant(0.f, im_dim, 1);
			}
			else {
				if (atomic_64bit)
					apu_co = af::constant(0ULL, 1, 1, u64);
				else
					apu_co = af::constant(0.f, 1, 1);
			}

			if ((MethodList.ACOSEM || MethodList.OSLCOSEM == 1) && alku == 0u) {
				if (atomic_64bit)
					apu_aco = af::constant(0ULL, 1, 1, u64);
				else
					apu_aco = af::constant(0.f, im_dim, 1);
			}
			else {
				if (atomic_64bit)
					apu_aco = af::constant(0ULL, 1, 1, u64);
				else
					apu_aco = af::constant(0.f, 1, 1);
			}
			if (TOF && !loadTOF && osa_iter > 0) {
				gpuErrchk(cuMemFree(d_Sino[0]));
#ifdef MX_HAS_INTERLEAVED_COMPLEX
				float* apu = (float*)mxGetSingles(mxGetCell(Sin, 0));
#else
				float* apu = (float*)mxGetData(mxGetCell(Sin, 0));
#endif
				status = cuCtxSynchronize();
				status = cuMemAlloc(&d_Sino[0], sizeof(float) * length[osa_iter] * nBins);
				for (int64_t to = 0LL; to < nBins; to++) {
					CUdeviceptr dst = reinterpret_cast<CUdeviceptr> (reinterpret_cast<char*>(d_Sino[0]) + sizeof(float) * length[osa_iter] * to);
					char* src = (char*)Sin + sizeof(float) * (pituus[osa_iter] + koko_l * to);
					size_t bytes = sizeof(float) * length[osa_iter];
					status = cuMemcpyHtoD(dst, src, bytes);
				}
			}
		}
		else {

			if (w_vec.MBSREM_prepass && alku == 0u) {
				if (atomic_64bit) {
					apu_summa = af::constant(0ULL, im_dim, 1, u64);
				}
				else {
					apu_summa = Summ[0];
				}
			}
			else {
				if (atomic_64bit) {
					apu_summa = af::constant(0ULL, 1, 1, u64);
				}
				else {
					apu_summa = af::constant(0.f, 1, 1);
				}
			}
			apu_summa_m = E;

			if ((MethodList.COSEM || MethodList.ECOSEM || MethodList.OSLCOSEM == 2) && alku == 0u) {
				if (atomic_64bit)
					apu_co = (C_co * TH).as(u64);
				else
					apu_co = C_co;
			}
			else {
				if (atomic_64bit)
					apu_co = af::constant(0ULL, 1, 1, u64);
				else
					apu_co = af::constant(0.f, 1, 1);
			}

			if ((MethodList.ACOSEM || MethodList.OSLCOSEM == 1) && alku == 0) {
				if (atomic_64bit)
					apu_aco = (C_aco * TH).as(u64);
				else
					apu_aco = C_aco;
			}
			else {
				if (atomic_64bit)
					apu_aco = af::constant(0ULL, 1, 1, u64);
				else
					apu_aco = af::constant(0.f, 1, 1);
			}

			if ((MethodList.MRAMLA || MethodList.MBSREM) && w_vec.MBSREM_prepass && alku == 0u)
				apu_Amin = w_vec.Amin;
			else
				apu_Amin = af::constant(0.f, 1, 1);
		}

		uint32_t H = osa_iter;
		if (TOF && !loadTOF)
			H = 0U;

		CUdeviceptr* d_apu_co = apu_co.device<CUdeviceptr>();
		CUdeviceptr* d_apu_aco = apu_aco.device<CUdeviceptr>();
		CUdeviceptr* d_E = apu_summa_m.device<CUdeviceptr>();
		CUdeviceptr* d_Amin = apu_Amin.device<CUdeviceptr>();
		CUdeviceptr* d_Summ = apu_summa.device<CUdeviceptr>();

		af::sync();

		if (CT) {
			void* args[] = { (void*)&global_factor, (void*)&epps, (void*)&im_dim, (void*)&Nx, (void*)&Ny, (void*)&Nz, (void*)&dz, (void*)&dx, (void*)&dy, (void*)&bz,
				(void*)&bx, (void*)&by,	(void*)&bzb, (void*)&maxxx, (void*)&maxyy, (void*)&zmax, (void*)&NSlices, (void*)&size_x, (void*)&TotSinos,
				(void*)&det_per_ring, (void*)&prows, (void*)&Nxy, (void*)&fp, (void*)&sigma_x, (void*)&tube_width, (void*)&crystal_size_z, (void*)&bmin, (void*)&bmax,
				(void*)&Vmax, &w_vec.epsilon_mramla,&d_TOFCenter, &d_atten, &d_pseudos, &d_x, &d_y, &d_z, &d_xcenter, &d_ycenter, &d_zcenter, &d_V, &d_reko_type,
				(void*)&subsets, &d_angles, (void*)&w_vec.size_y, (void*)&w_vec.dPitch, (void*)&w_vec.nProjections, &d_norm[osa_iter], &d_scat[osa_iter],
				reinterpret_cast<void*>(&d_Summ),&d_lor[osa_iter], &d_xyindex[osa_iter], &d_zindex[osa_iter], &d_L[osa_iter], &d_Sino[H], &d_sc_ra[osa_iter],
				reinterpret_cast<void*>(&d_cosem), &alku, &MBSREM_prepass, reinterpret_cast<void*>(&d_ACOSEM_lhs), reinterpret_cast<void*>(&d_Amin),
				reinterpret_cast<void*>(&d_apu_co), reinterpret_cast<void*>(&d_apu_aco), reinterpret_cast<void*>(&d_E), (void*)&m_size, &MethodList, (void*)&st };
			status = cuLaunchKernel(kernel_mbsrem, global_size, 1, 1, local_size, 1, 1, 0, af_cuda_stream, &args[0], 0);
		}
		else {
			if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1))) {
				void* args[] = { (void*)&global_factor, (void*)&epps, (void*)&im_dim, (void*)&Nx, (void*)&Ny, (void*)&Nz, (void*)&dz, (void*)&dx, (void*)&dy, (void*)&bz, 
					(void*)&bx, (void*)&by, (void*)&bzb, (void*)&maxxx, (void*)&maxyy, (void*)&zmax, (void*)&NSlices, (void*)&size_x, (void*)&TotSinos, 
					(void*)&det_per_ring, (void*)&prows, (void*)&Nxy, (void*)&fp, (void*)&sigma_x, (void*)&tube_width, (void*)&crystal_size_z, (void*)&bmin, (void*)&bmax,
					(void*)&Vmax, &w_vec.epsilon_mramla,&d_TOFCenter, &d_atten, &d_pseudos, &d_x, &d_y, &d_z, &d_xcenter, &d_ycenter, &d_zcenter, &d_V, &d_reko_type, 
					&d_norm[osa_iter], &d_scat[osa_iter], reinterpret_cast<void*>(&d_Summ),&d_lor[osa_iter], &d_xyindex[osa_iter], &d_zindex[osa_iter], &d_L[osa_iter], 
					&d_Sino[H], &d_sc_ra[osa_iter], reinterpret_cast<void*>(&d_cosem), &alku, &MBSREM_prepass, reinterpret_cast<void*>(&d_ACOSEM_lhs), 
					reinterpret_cast<void*>(&d_Amin), reinterpret_cast<void*>(&d_apu_co), reinterpret_cast<void*>(&d_apu_aco), reinterpret_cast<void*>(&d_E), (void*)&m_size, 
					&MethodList, (void*)&st };
				status = cuLaunchKernel(kernel_mbsrem, global_size, 1, 1, local_size, 1, 1, 0, af_cuda_stream, &args[0], 0);
			}
			else {
				void* args[] = { (void*)&global_factor, (void*)&epps, (void*)&im_dim, (void*)&Nx, (void*)&Ny, (void*)&Nz, (void*)&dz, (void*)&dx, (void*)&dy, (void*)&bz, 
					(void*)&bx, (void*)&by, (void*)&bzb, (void*)&maxxx, (void*)&maxyy, (void*)&zmax, (void*)&NSlices, (void*)&size_x, (void*)&TotSinos,
					(void*)&det_per_ring, (void*)&prows, (void*)&Nxy, (void*)&fp, (void*)&sigma_x, (void*)&dc_z, (void*)&n_rays, &w_vec.epsilon_mramla,&d_TOFCenter, &d_atten, 
					&d_pseudos, &d_x, &d_y, &d_z, &d_reko_type, &d_norm[osa_iter], &d_scat[osa_iter], reinterpret_cast<void*>(&d_Summ), &d_lor[osa_iter], &d_xyindex[osa_iter], 
					&d_zindex[osa_iter],&d_L[osa_iter],&d_Sino[H],&d_sc_ra[osa_iter], reinterpret_cast<void*>(&d_cosem),&alku,&MBSREM_prepass, reinterpret_cast<void*>(&d_ACOSEM_lhs), 
					reinterpret_cast<void*>(&d_Amin), reinterpret_cast<void*>(&d_apu_co), reinterpret_cast<void*>(&d_apu_aco), reinterpret_cast<void*>(&d_E), (void*)&m_size, 
					&MethodList, (void*)&st };
				status = cuLaunchKernel(kernel_mbsrem, global_size, 1, 1, local_size, 1, 1, 0, af_cuda_stream, &args[0], 0);
			}
		}
		if ((status != CUDA_SUCCESS)) {
			std::cerr << getErrorString(status) << std::endl;
			mexPrintf("Failed to launch the prepass kernel\n");
			return;
		}
		else if (DEBUG) {
			mexPrintf("Prepass kernel launched\n");
			mexEvalString("pause(.0001);");
		}
		status = cuCtxSynchronize();
		if ((status != CUDA_SUCCESS)) {
			std::cerr << getErrorString(status) << std::endl;
			mexPrintf("Failed to synchronize\n");
			return;
		}

		//clFinish(af_queue);
		apu_co.unlock();
		apu_aco.unlock();
		apu_summa.unlock();
		apu_summa_m.unlock();
		apu.unlock();
		apu_Amin.unlock();
		cosem_psf.unlock();
		st += length[osa_iter];

		if (alku == 0u) {
			if ((MethodList.COSEM || MethodList.ECOSEM)) {
				if (atomic_64bit)
					C_co(af::span, osa_iter) = apu_co.as(f32) / TH;
				else
					C_co(af::span, osa_iter) = apu_co;
				if (use_psf)
					C_co(af::span, osa_iter) = computeConvolution(C_co(af::span, osa_iter), g, Nx, Ny, Nz, w_vec, 1u) * cosem;
				else
					C_co(af::span, osa_iter) = C_co(af::span, osa_iter) * cosem;
			}
			if (MethodList.ACOSEM) {
				if (atomic_64bit)
					C_aco(af::span, osa_iter) = apu_aco.as(f32) / TH;
				else
					C_aco(af::span, osa_iter) = apu_aco;
				if (use_psf)
					C_aco(af::span, osa_iter) = computeConvolution(C_aco(af::span, osa_iter), g, Nx, Ny, Nz, w_vec, 1u) * af::pow(cosem, w_vec.h_ACOSEM_2);
				else
					C_aco(af::span, osa_iter) = C_aco(af::span, osa_iter) * af::pow(cosem, w_vec.h_ACOSEM_2);
			}
			if (MethodList.OSLCOSEM == 2u) {
				if (atomic_64bit)
					C_osl(af::span, osa_iter) = apu_co.as(f32) / TH;
				else
					C_osl(af::span, osa_iter) = apu_co;
				if (use_psf)
					C_osl(af::span, osa_iter) = computeConvolution(C_osl(af::span, osa_iter), g, Nx, Ny, Nz, w_vec, 1u) * cosem;
				else
					C_osl(af::span, osa_iter) = C_osl(af::span, osa_iter) * cosem;
			}
			else if (MethodList.OSLCOSEM == 1) {
				if (atomic_64bit)
					C_osl(af::span, osa_iter) = apu_aco.as(f32) / TH;
				else
					C_osl(af::span, osa_iter) = apu_aco;
				if (use_psf)
					C_osl(af::span, osa_iter) = computeConvolution(C_osl(af::span, osa_iter), g, Nx, Ny, Nz, w_vec, 1u) * af::pow(cosem, w_vec.h_ACOSEM_2);
				else
					C_osl(af::span, osa_iter) = C_osl(af::span, osa_iter) * af::pow(cosem, w_vec.h_ACOSEM_2);
			}


			if ((MethodList.MRAMLA || MethodList.MBSREM) && w_vec.MBSREM_prepass && Nt > 1U) {
				sub_index_array = af::range(af::dim4(length[osa_iter] * nBins), 0, u32) + pituus[osa_iter] * nBins;
				if (subsets > 1)
					w_vec.Amin(sub_index_array) = apu_Amin;
				else
					w_vec.Amin = apu_Amin;
			}

			if (w_vec.MBSREM_prepass) {
				if (compute_norm_matrix == 0u) {
					if (atomic_64bit)
						Summ[osa_iter] = (apu_summa).as(f32) / TH;
					else
						Summ[osa_iter] = apu_summa;
					w_vec.D += Summ[osa_iter];
					Summ[osa_iter](Summ[osa_iter] < epps) = epps;
					if (use_psf)
						Summ[osa_iter] = computeConvolution(Summ[osa_iter], g, Nx, Ny, Nz, w_vec, 1u);
				}
				else {
					Summ[0] = apu_summa;
					w_vec.D += apu_summa.as(f32);
				}

				if ((MethodList.MRAMLA || MethodList.MBSREM) && Nt > 1U) {
					if (subsets > 1)
						E(sub_index_array) = apu_summa_m;
					else
						E = apu_summa_m;
				}
			}
			if (alku == 0U && (MethodList.MBSREM || MethodList.MRAMLA) && w_vec.MBSREM_prepass && Nt == 1U) {
				uint32_t H = osa_iter;
				uint32_t L = 0U;
				if (TOF && !loadTOF)
					H = 0;
				if (randoms_correction)
					L = osa_iter;
				float* hOut = new float[length[osa_iter] * nBins];
				cuMemcpyDtoH(hOut, d_Sino[H], length[osa_iter] * sizeof(float) * nBins);
				cuCtxSynchronize();
				af::array Sino(length[osa_iter] * nBins, hOut, afHost);
				af::eval(Sino);
				af::array rand;
				if (randoms_correction) {
					float* hOutR = new float[length[osa_iter]];
					cuMemcpyDtoH(hOutR, d_sc_ra[L], length[osa_iter] * sizeof(float));
					cuCtxSynchronize();
					rand = af::array(length[osa_iter], hOutR, afHost);
					af::eval(rand);
					af::sync();
					delete[] hOutR;
				}
				if (U_skip) {
					float UU = w_vec.U;
					const af::array Aind = apu_Amin > 0.f;
					w_vec.U = af::max<float>(Sino(Aind) / apu_Amin(Aind));
					if (UU > w_vec.U || std::isinf(w_vec.U))
						w_vec.U = UU;
				}
				float eps_mramla = w_vec.epsilon_mramla;
				w_vec.epsilon_mramla = MBSREM_epsilon(Sino, epps, randoms_correction, rand, apu_summa_m, TOF, nBins, CT);
				if (eps_mramla < w_vec.epsilon_mramla)
					w_vec.epsilon_mramla = eps_mramla;
				delete[] hOut;
			}
			af::deviceGC();
		}
		else {
			w_vec.ACOSEM_rhs = af::sum<float>(apu);
		}
	}
	if (use_psf && alku == 0 && w_vec.MBSREM_prepass) {
		w_vec.D = computeConvolution(w_vec.D, g, Nx, Ny, Nz, w_vec, 1u);
	}
	if (w_vec.MBSREM_prepass)
		w_vec.D(w_vec.D <= 0.f) = 1.f;
	return;
}

nvrtcResult createProgramCUDA(const bool verbose, const char* k_path, const char* fileName,
	nvrtcProgram& program_os, nvrtcProgram& program_ml, nvrtcProgram& program_mbsrem, bool& atomic_64bit, const char* header_directory,
	const uint32_t projector_type, const float crystal_size_z, const bool precompute, const uint8_t raw, const uint32_t attenuation_correction,
	const uint32_t normalization_correction, const int32_t dec, const size_t local_size, const uint16_t n_rays, const uint16_t n_rays3D,
	const RecMethods MethodList, const bool osem_bool, const bool mlem_bool, const uint32_t n_rekos, const uint32_t n_rekos_mlem,
	const Weighting& w_vec, const uint32_t osa_iter0, const float cr_pz, const float dx, const bool use_psf, const uint32_t scatter,
	const uint32_t randoms_correction, const bool TOF, const int64_t nBins, const uint8_t listmode, const bool CT) {

	nvrtcResult status = NVRTC_SUCCESS;

	const char* options[50];
	int uu = 0;

	options[uu] = header_directory;
	uu++;


	char buffer1[30];
	char buffer2[30];
	char buffer3[30];
	char buffer4[30];
	char buffer5[30];
	char buffer6[30];
	char buffer7[30];
	char buffer8[30];
	char buffer9[30];
	if (crystal_size_z == 0.f && projector_type == 2u) {
		options[uu] = "-DCRYST";
		uu++;
	}
	if ((crystal_size_z > 0.f && projector_type == 2u) || projector_type == 3u) {
		options[uu] = "-DCRYSTZ";
		uu++;
	}
	if (projector_type == 3u) {
		options[uu] = "-DVOL";
		uu++;
	}
	if (precompute) {
		options[uu] = "-DPRECOMPUTE";
		uu++;
	}
	if (raw == 1) {
		options[uu] = "-DRAW";
		uu++;
	}
	if (projector_type == 1u) {
		options[uu] = "-DSIDDON";
		uu++;
	}
	else if (projector_type == 2u || projector_type == 3u) {
		options[uu] = "-DORTH";
		uu++;
	}
	if (attenuation_correction == 1u) {
		options[uu] = "-DATN";
		uu++;
	}
	if (normalization_correction == 1u) {
		options[uu] = "-DNORM";
		uu++;
	}
	if (randoms_correction == 1u) {
		options[uu] = "-DRANDOMS";
		uu++;
	}
	if (TOF && projector_type == 1u) {
		options[uu] = "-DTOF";
		uu++;
	}
	if (CT) {
		options[uu] = "-DCT";
		uu++;
	}
	if (listmode == 1) {
		options[uu] = "-DLISTMODE";
		uu++;
	}
	else if (listmode == 2) {
		options[uu] = "-DLISTMODE2";
		uu++;
	}
	options[uu] = "-DFP";
	uu++;

//#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__) || defined(_WIN64)) && defined(_MSC_VER)
//	sprintf_s(buffer2, 30, "-DNBINS=%ll", nBins);
//#else
	std::snprintf(buffer2, 30,"-DNBINS=%d", static_cast<int32_t>(nBins));
//#endif
	options[uu] = buffer2;
	uu++;

	if (scatter == 1u) {
		options[uu] = "-DSCATTER";
		uu++;
	}
	if (MethodList.NLM) {
		options[uu] = "-DNLM_";
		uu++;
	}
	if (projector_type == 1u && !precompute && (n_rays * n_rays3D) > 1) {
//#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__) || defined(_WIN64)) && defined(_MSC_VER)
//		sprintf_s(buffer3, 30, "-DN_RAYS=%u", n_rays* n_rays3D);
//#else
		std::snprintf(buffer3, 30, "-DN_RAYS=%u", n_rays* n_rays3D);
//#endif
		options[uu] = buffer3;
		uu++;
		if (n_rays > 1) {
			options[uu] = "-DN_RAYS2D";
			uu++;
		}
//#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__) || defined(_WIN64)) && defined(_MSC_VER)
//		sprintf_s(buffer1, 30, "-DN_RAYS3D=%u", n_rays3D);
//#else
		std::snprintf(buffer1, 30, "-DN_RAYS3D=%u", n_rays3D);
//#endif
		options[uu] = buffer1;
		uu++;
	}
	if (MethodList.MRP) {
		options[uu] = "-DMEDIAN";
		uu++;

//#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__) || defined(_WIN64)) && defined(_MSC_VER)
//		sprintf_s(buffer4, 30, "-DSEARCH_WINDOW_X=%u", w_vec.Ndx);
//#else
		std::snprintf(buffer4, 30, "-DSEARCH_WINDOW_X=%u", w_vec.Ndx);
//#endif
		options[uu] = buffer4;
		uu++;

//#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__) || defined(_WIN64)) && defined(_MSC_VER)
//		sprintf_s(buffer5, 30, "-DSEARCH_WINDOW_Y=%u", w_vec.Ndy);
//#else
		std::snprintf(buffer5, 30, "-DSEARCH_WINDOW_Y=%u", w_vec.Ndy);
//#endif
		options[uu] = buffer5;
		uu++;

//#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__) || defined(_WIN64)) && defined(_MSC_VER)
//		sprintf_s(buffer6, 30, "-DSEARCH_WINDOW_Z=%u", w_vec.Ndz);
//#else
		std::snprintf(buffer6, 30, "-DSEARCH_WINDOW_Z=%u", w_vec.Ndz);
//#endif
		options[uu] = buffer6;
		uu++;
	}
	if ((projector_type == 2U || projector_type == 3U || TOF) && dec > 0) {
//#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__) || defined(_WIN64)) && defined(_MSC_VER)
//		sprintf_s(buffer8, 30, "-DDEC=%u", dec);
//#else
		std::snprintf(buffer8, 30, "-DDEC=%u", dec);
//#endif
		options[uu] = buffer8;
		uu++;
	}
	if (DEBUG) {
		for (int gg = 0; gg < uu; gg++)
			mexPrintf("%s", options[gg]);
		mexPrintf("\n");
		mexPrintf("uu = %d\n", uu);
	}
	if (osem_bool) {
		int ll = uu;
		const char* os_options[50];
		
		for (int kk = 0; kk <= ll; kk++) {
			if (kk < ll) {
				os_options[kk] = options[kk];
			}
			else {
//#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__) || defined(_WIN64)) && defined(_MSC_VER)
//				sprintf_s(buffer7, 30, "-DN_REKOS=%u", n_rekos);
//#else
				std::snprintf(buffer7, 30, "-DN_REKOS=%u", n_rekos);
//#endif
				os_options[uu] = buffer7;
			}
		}
		ll++;
		os_options[ll] = "-DAF";
		ll++;
		if (n_rekos == 1) {
			os_options[ll] = "-DNREKOS1";
			ll++;
		}
		else if (n_rekos == 2) {
			os_options[ll] = "-DNREKOS2";
			ll++;
		}
		if (MethodList.MRAMLA || MethodList.MBSREM) {
			os_options[ll] = "-DMRAMLA";
			ll++;
		}
		if (MethodList.COSEM || MethodList.ACOSEM || MethodList.OSLCOSEM > 0u || MethodList.ECOSEM) {
			os_options[ll] = "-DCOSEM";
			ll++;
		}
		if (DEBUG) {
			for (int gg = 0; gg < ll; gg++)
				mexPrintf("%s\n", os_options[gg]);
			mexPrintf("ll = %d\n", ll);
		}

		status = buildProgramCUDA(verbose, k_path, program_os, atomic_64bit, os_options, ll);
	}
	if (mlem_bool) {
		int rr = uu;
		const char* ml_options[50];
		for (int kk = 0; kk <= rr; kk++) {
			if (kk < rr)
				ml_options[kk] = options[kk];
			else {
//#if (defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__) || defined(_WIN64)) && defined(_MSC_VER)
//				sprintf_s(buffer9, 30, "-DN_REKOS=%d", n_rekos_mlem);
//#else
				std::snprintf(buffer9, 30, "-DN_REKOS=%d", n_rekos_mlem);
//#endif
				ml_options[kk] = buffer9;
			}
		}
		rr++;
		ml_options[rr] = "-DAF";
		rr++;
		if (n_rekos_mlem == 1) {
			ml_options[rr] = "-DNREKOS1";
			rr++;
		}
		else if (n_rekos_mlem == 2) {
			ml_options[rr] = "-DNREKOS2";
			rr++;
		}
		status = buildProgramCUDA(verbose, k_path, program_ml, atomic_64bit, ml_options, rr);
	}
	if ((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL || MethodList.PKMA) && w_vec.MBSREM_prepass ||
		MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0) {
		options[uu] = "-DMBSREM";
		uu++;
		options[uu] = "-DAF";
		uu++;
		if (n_rekos == 1) {
			options[uu] = "-DNREKOS1";
			uu++;
		}
		else if (n_rekos == 2) {
			options[uu] = "-DNREKOS2";
			uu++;
		}
		if (MethodList.MRAMLA || MethodList.MBSREM) {
			options[uu] = "-DMRAMLA";
			uu++;
		}
		status = buildProgramCUDA(verbose, k_path, program_mbsrem, atomic_64bit, options, uu);
	}
	//delete[] buffer;
	return status;
}

nvrtcResult buildProgramCUDA(const bool verbose, const char* k_path, nvrtcProgram& program, bool& atomic_64bit, const char* options[], int uu) {
	nvrtcResult status = NVRTC_SUCCESS;
	size_t pituus;
	if (atomic_64bit) {
		options[uu] = "-DCAST=unsigned long long int";
		uu++;
		options[uu] = "-DATOMIC";
	}
	else {
		options[uu] = "-DCAST=float";
	}
	//for (int ll = uu + 1; ll < 50; ll++)
	//	options[ll] = "";
	//options.erase(options.end() - (options.size() - uu) + 1, options.end());
	if (DEBUG) {
		for (int ll = 0; ll <= uu; ll++)
			mexPrintf("%s", options[ll]);
		mexPrintf("\n");
		mexPrintf("uu = %d\n", uu);
		//mexPrintf("options.size() = %d\n", options.size());
	}
	if (atomic_64bit) {
		std::string kernel_path_atom;

		kernel_path_atom = k_path;
		kernel_path_atom += ".cu";
		// Load the source text file
		std::fstream sourceFile_atom(kernel_path_atom.c_str());
		std::string content_atom((std::istreambuf_iterator<char>(sourceFile_atom)), std::istreambuf_iterator<char>());
		const char* sourceCode_atom = new char[content_atom.size()];
		sourceCode_atom = content_atom.c_str();
		// Create the program from the source
		status = nvrtcCreateProgram(&program, sourceCode_atom, "64bit_atom", 0, NULL, NULL);

		//const char* testi = options.c_str();
		//const char* const* testi2 = reinterpret_cast<const char* const*>(testi);

		// Build the program
		//status = nvrtcCompileProgram(program, uu, testi2);
		status = nvrtcCompileProgram(program, uu + 1, options);
		if (status != NVRTC_SUCCESS) {
			mexPrintf("Failed to build 64-bit atomics program.\n");
			if (DEBUG) {
				std::cerr << nvrtcGetErrorString(status) << std::endl;
				mexPrintf("Failed to build CUDA program. Build log: \n");
				size_t len;
				char* buffer;
				nvrtcGetProgramLogSize(program, &len);
				buffer = (char*)calloc(len, sizeof(size_t));
				nvrtcGetProgramLog(program, buffer);
				mexPrintf("%s\n", buffer);
				free(buffer);
				nvrtcDestroyProgram(&program);
				return status;
			}
			options[uu] = "";
			uu--;
			options[uu] = "-DCAST=float";
			//options.erase(options.end() - 2, options.end());
			//options.emplace_back("-DCAST=float");
		}
		else if (verbose)
			mexPrintf("CUDA program (64-bit atomics) built\n");
	}
	else
		status = NVRTC_ERROR_PROGRAM_CREATION_FAILURE;
	// If not, use 32-bit atomic add (float)
	if (status != NVRTC_SUCCESS) {
		status = NVRTC_SUCCESS;
		atomic_64bit = false;

		std::string kernel_path;

		kernel_path = k_path;
		kernel_path += ".cu";
		std::fstream sourceFile(kernel_path.c_str());
		std::string content((std::istreambuf_iterator<char>(sourceFile)), std::istreambuf_iterator<char>());
		const char* sourceCode = new char[content.size()];
		sourceCode = content.c_str();
		status = nvrtcCreateProgram(&program, sourceCode, "32bit", 0, NULL, NULL);
		if (status != NVRTC_SUCCESS) {
			std::cerr << nvrtcGetErrorString(status) << std::endl;
			return status;
		}
		//const char* testi = options.c_str();
		//const char* const* testi2 = reinterpret_cast<const char* const*>(testi);
		//status = nvrtcCompileProgram(program, uu, testi2);

		// Build the program
		//status = nvrtcCompileProgram(program, options.size(), options.data());
		status = nvrtcCompileProgram(program, uu + 1, options);
		// Build log in case of failure
		if (status != NVRTC_SUCCESS) {
			std::cerr << nvrtcGetErrorString(status) << std::endl;
			mexPrintf("Failed to build CUDA program. Build log: \n");
			size_t len;
			char* buffer;
			nvrtcGetProgramLogSize(program, &len);
			buffer = (char*)calloc(len, sizeof(size_t));
			nvrtcGetProgramLog(program, buffer);
			mexPrintf("%s\n", buffer);
			free(buffer);
			nvrtcDestroyProgram(&program);
			return status;
		}
		else if (verbose)
			mexPrintf("CUDA program built\n");
	}
	return status;
}

nvrtcResult createKernelsCUDA(const bool verbose, nvrtcProgram& program_os, nvrtcProgram& program_ml, nvrtcProgram& program_mbsrem, CUfunction& kernel_os, CUfunction& kernel_ml,
	CUfunction& kernel_mbsrem, CUfunction& kernelNLM, CUfunction& kernelMed, const bool osem_bool, const bool mlem_bool, const RecMethods& MethodList, const Weighting& w_vec, const bool precompute, const uint32_t projector_type,
	const uint16_t n_rays, const uint16_t n_rays3D, CUmodule& moduleOS, CUmodule& moduleML, CUmodule& moduleMB) {

	nvrtcResult status = NVRTC_SUCCESS;
	CUresult status2 = CUDA_SUCCESS;

	if (osem_bool) {
		size_t ptxSize;
		status = nvrtcGetPTXSize(program_os, &ptxSize);
		if (status != NVRTC_SUCCESS) {
			std::cerr << nvrtcGetErrorString(status) << std::endl;
			mexPrintf("Unable to get the PTX size\n");
			return status;
		}
		char* ptx = new char[ptxSize];
		status = nvrtcGetPTX(program_os, ptx);
		if (status != NVRTC_SUCCESS) {
			std::cerr << nvrtcGetErrorString(status) << std::endl;
			mexPrintf("Unable to get the PTX\n");
			return status;
		}
		// Destroy the program.
		status = nvrtcDestroyProgram(&program_os);
		if (status != NVRTC_SUCCESS) {
			std::cerr << nvrtcGetErrorString(status) << std::endl;
			mexPrintf("Unable to destroy the program\n");
			return status;
		}
		status2 = cuModuleLoadData(&moduleOS, ptx);
		if (status2 != CUDA_SUCCESS) {
			std::cerr << getErrorString(status2) << std::endl;
			mexPrintf("Unable to load the module data\n");
			return NVRTC_ERROR_PROGRAM_CREATION_FAILURE;
		}
		if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && ((precompute || (n_rays * n_rays3D) == 1))))
			status2 = cuModuleGetFunction(&kernel_os, moduleOS, "kernel_multi");
		else
			status2 = cuModuleGetFunction(&kernel_os, moduleOS, "siddon_multi");
		if (status2 != CUDA_SUCCESS) {
			std::cerr << getErrorString(status2) << std::endl;
			mexPrintf("Unable to find the kernel function\n");
			return NVRTC_ERROR_PROGRAM_CREATION_FAILURE;
		}
		if (MethodList.NLM) {
			status2 = cuModuleGetFunction(&kernelNLM, moduleOS, "NLM");
			if (status2 != CUDA_SUCCESS) {
				std::cerr << getErrorString(status2) << std::endl;
				mexPrintf("Unable to find the NLM kernel function\n");
				return NVRTC_ERROR_PROGRAM_CREATION_FAILURE;
			}
		}
		if (MethodList.MRP) {
			status2 = cuModuleGetFunction(&kernelMed, moduleOS, "medianFilter3D");
			if (status2 != CUDA_SUCCESS) {
				std::cerr << getErrorString(status2) << std::endl;
				mexPrintf("Unable to find the median kernel function\n");
				return NVRTC_ERROR_PROGRAM_CREATION_FAILURE;
			}
		}
		delete[] ptx;
	}

	if (mlem_bool) {
		size_t ptxSize;
		status = nvrtcGetPTXSize(program_ml, &ptxSize);
		if (status != NVRTC_SUCCESS) {
			std::cerr << nvrtcGetErrorString(status) << std::endl;
			return status;
		}
		char* ptx = new char[ptxSize];
		status = nvrtcGetPTX(program_ml, ptx);
		if (status != NVRTC_SUCCESS) {
			std::cerr << nvrtcGetErrorString(status) << std::endl;
			return status;
		}
		// Destroy the program.
		status = nvrtcDestroyProgram(&program_ml);
		if (status != NVRTC_SUCCESS) {
			std::cerr << nvrtcGetErrorString(status) << std::endl;
			return status;
		}
		status2 = cuModuleLoadData(&moduleML, ptx);
		if (status2 != CUDA_SUCCESS) {
			std::cerr << getErrorString(status2) << std::endl;
			return NVRTC_ERROR_PROGRAM_CREATION_FAILURE;
		}
		if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && ((precompute || (n_rays * n_rays3D) == 1))))
			status2 = cuModuleGetFunction(&kernel_ml, moduleML, "kernel_multi");
		else
			status2 = cuModuleGetFunction(&kernel_ml, moduleML, "siddon_multi");
		if (status2 != CUDA_SUCCESS) {
			std::cerr << getErrorString(status2) << std::endl;
			return NVRTC_ERROR_PROGRAM_CREATION_FAILURE;
		}
		if (MethodList.NLM && !osem_bool) {
			status2 = cuModuleGetFunction(&kernelNLM, moduleML, "NLM");
			if (status2 != CUDA_SUCCESS) {
				std::cerr << getErrorString(status2) << std::endl;
				mexPrintf("Unable to find the NLM kernel function\n");
				return NVRTC_ERROR_PROGRAM_CREATION_FAILURE;
			}
		}
		if (MethodList.MRP && !osem_bool) {
			status2 = cuModuleGetFunction(&kernelMed, moduleML, "medianFilter3D");
			if (status2 != CUDA_SUCCESS) {
				std::cerr << getErrorString(status2) << std::endl;
				mexPrintf("Unable to find the median kernel function\n");
				return NVRTC_ERROR_PROGRAM_CREATION_FAILURE;
			}
		}
		delete[] ptx;
	}
	if ((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL) && w_vec.MBSREM_prepass ||
		MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.PKMA || MethodList.OSLCOSEM > 0) {

		size_t ptxSize;
		status = nvrtcGetPTXSize(program_mbsrem, &ptxSize);
		if (status != NVRTC_SUCCESS) {
			std::cerr << nvrtcGetErrorString(status) << std::endl;
			return status;
		}
		char* ptx = new char[ptxSize];
		status = nvrtcGetPTX(program_mbsrem, ptx);
		if (status != NVRTC_SUCCESS) {
			std::cerr << nvrtcGetErrorString(status) << std::endl;
			return status;
		}
		// Destroy the program.
		status = nvrtcDestroyProgram(&program_mbsrem);
		if (status != NVRTC_SUCCESS) {
			std::cerr << nvrtcGetErrorString(status) << std::endl;
			return status;
		}
		status2 = cuModuleLoadData(&moduleMB, ptx);
		if (status2 != CUDA_SUCCESS) {
			std::cerr << getErrorString(status2) << std::endl;
			return NVRTC_ERROR_PROGRAM_CREATION_FAILURE;
		}
		if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && ((precompute || (n_rays * n_rays3D) == 1))))
			status2 = cuModuleGetFunction(&kernel_mbsrem, moduleMB, "kernel_multi");
		else
			status2 = cuModuleGetFunction(&kernel_mbsrem, moduleMB, "siddon_multi");
		if (status2 != CUDA_SUCCESS) {
			std::cerr << getErrorString(status2) << std::endl;
			return NVRTC_ERROR_PROGRAM_CREATION_FAILURE;
		}
		delete[] ptx;
	}

	return status;
}


af::array NLM(const af::array& im, Weighting& w_vec, const float epps, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const kernelStruct& CUDAStruct)
{
	CUresult status = CUDA_SUCCESS;
	const int32_t ndx = static_cast<int32_t>(w_vec.Ndx);
	const int32_t ndy = static_cast<int32_t>(w_vec.Ndy);
	const int32_t ndz = static_cast<int32_t>(w_vec.Ndz);
	const int32_t nlx = static_cast<int32_t>(w_vec.Nlx);
	const int32_t nly = static_cast<int32_t>(w_vec.Nly);
	const int32_t nlz = static_cast<int32_t>(w_vec.Nlz);
	af::array grad;
	//af::array input = af::moddims(im, Nx, Ny, Nz);
	af::array padInput;
	const uint32_t padx = w_vec.Ndx + w_vec.Nlx;
	const uint32_t pady = w_vec.Ndy + w_vec.Nly;
	const uint32_t padz = w_vec.Ndz + w_vec.Nlz;
	if (w_vec.NLM_anatomical) {
		padInput = padding(w_vec.NLM_ref, Nx, Ny, Nz, padx, pady, padz);
	}
	else {
		padInput = padding(im, Nx, Ny, Nz, padx, pady, padz);
	}
	af::array input = padding(im, Nx, Ny, Nz, padx, pady, padz);
	padInput = af::flat(padInput);
	const uint32_t N = input.dims(0);
	const uint32_t M = input.dims(1);
	const uint32_t K = input.dims(2);
	const int32_t Nxy = N * M;
	int32_t type = 0;
	if (w_vec.NLTV)
		type = 1;
	else if (w_vec.NLM_MRP)
		type = 2;
	else
		type = 0;
	input = af::flat(input); 
	size_t local_size = 64ULL;
	size_t erotus = (N * M * K) % local_size;

	if (erotus > 0)
		erotus = (local_size - erotus);

	size_t global_size = (N * M * K) + erotus;

	global_size = global_size / local_size;
	const int window_x = w_vec.Ndx + w_vec.Nlx;
	const int window_y = w_vec.Ndy + w_vec.Nly;
	const int window_z = w_vec.Ndz + w_vec.Nlz;
	const int min_x = window_x;
	const int max_x = static_cast<int32_t>(N) - window_x;
	const int min_y = window_y;
	const int max_y = static_cast<int32_t>(M) - window_y;
	const int min_z = window_z;
	const int max_z = static_cast<int32_t>(K) - window_z;
	af::array W = af::constant(0.f, N, M, K);
	W = af::flat(W);
	CUdeviceptr* d_W = W.device<CUdeviceptr>();
	CUdeviceptr* d_input = input.device<CUdeviceptr>();
	CUdeviceptr* d_padInput = padInput.device<CUdeviceptr>();
	CUdeviceptr* d_gaussianNLM = w_vec.gaussianNLM.device<CUdeviceptr>();
	af::sync();
	// Compute the kernel
	void* args[] = { reinterpret_cast<void*>(&d_W), reinterpret_cast<void*>(&d_input) , reinterpret_cast<void*>(&d_padInput),
		reinterpret_cast<void*>(&d_gaussianNLM), (void*)&ndx, (void*)&ndy , (void*)&ndz, (void*)&nlx, (void*)&nly, (void*)&nlz,
		(void*)&N , (void*)&M , (void*)&K, &w_vec.h2, (void*)&epps, (void*)&Nxy, (void*)&min_x, (void*)&max_x, (void*)&min_y,
		(void*)&max_y , (void*)&min_z, (void*)&max_z, &type};
	status = cuLaunchKernel(CUDAStruct.kernelNLM, global_size, 1, 1, local_size, 1, 1, 0, *CUDAStruct.af_cuda_stream, &args[0], 0);
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		mexPrintf("Failed to launch the NLM kernel\n");
		mexEvalString("pause(.0001);");
	}
	status = cuCtxSynchronize();
	if (status != CUDA_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		mexPrintf("Queue finish failed after kernel\n");
		mexEvalString("pause(.0001);");
	}

	W.unlock();
	input.unlock();
	padInput.unlock();
	w_vec.gaussianNLM.unlock();
	W = af::moddims(W, N, M, K);
	W = W(af::seq(padx, Nx + padx - 1), af::seq(pady, Ny + pady - 1), af::seq(padz, Nz + padz - 1));
	W = af::flat(W);

	if (w_vec.NLM_MRP) {
		grad = im - W;
	}
	else
		grad = W;
	af::sync();
	return grad;
}