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
void update_opencl_inputs(AF_im_vectors & vec, OpenCL_im_vectors &vec_opencl, const bool mlem, const uint32_t im_dim, const uint32_t n_rekos, 
	const uint32_t n_rekos_mlem, const RecMethods MethodList, const bool atomic_64bit, const bool use_psf)
{
	if (MethodList.CUSTOM) {
		if (mlem) {
			vec.im_mlem(af::seq(n_rekos * im_dim - im_dim, n_rekos * im_dim - 1u)) = vec.custom_MLEM;
		}
		else {
			uint32_t yy = n_rekos * im_dim;
			if (MethodList.OSLCOSEM > 0u) {
				vec.im_os(af::seq(yy - im_dim, yy - 1u)) = vec.custom_COSEM;
				yy -= im_dim;
			}
			if (MethodList.RBIMAP) {
				vec.im_os(af::seq(yy - im_dim, yy - 1u)) = vec.custom_RBI;
				yy -= im_dim;
			}
			if (MethodList.ROSEMMAP) {
				vec.im_os(af::seq(yy - im_dim, yy - 1u)) = vec.custom_ROSEM;
				yy -= im_dim;
			}
			if (MethodList.MBSREM) {
				vec.im_os(af::seq(yy - im_dim, yy - 1u)) = vec.custom_MBSREM;
				yy -= im_dim;
			}
			if (MethodList.BSREM) {
				vec.im_os(af::seq(yy - im_dim, yy - 1u)) = vec.custom_BSREM;
				yy -= im_dim;
			}
			if (MethodList.OSLOSEM) {
				vec.im_os(af::seq(yy - im_dim, yy - 1u)) = vec.custom_OSEM;
				yy -= im_dim;
			}
		}
	}
	if (mlem) {
		if (use_psf)
			vec_opencl.d_im_mlem = vec.im_mlem_blurred.device<cl_mem>();
		else
			vec_opencl.d_im_mlem = vec.im_mlem.device<cl_mem>();
		if (atomic_64bit)
			vec.rhs_mlem = af::constant(0ULL, static_cast<size_t>(im_dim) * n_rekos_mlem, 1, u64);
		else
			vec.rhs_mlem = af::constant(0.f, static_cast<size_t>(im_dim) * n_rekos_mlem, 1);
		vec_opencl.d_rhs_mlem = vec.rhs_mlem.device<cl_mem>();
	}
	else {
		if (use_psf)
			vec_opencl.d_im_os = vec.im_os_blurred.device<cl_mem>();
		else
			vec_opencl.d_im_os = vec.im_os.device<cl_mem>();
		if (atomic_64bit)
			vec.rhs_os = af::constant(0ULL, static_cast<size_t>(im_dim) * static_cast<size_t>(n_rekos), 1, u64);
		else
			vec.rhs_os = af::constant(0.f, static_cast<size_t>(im_dim) * static_cast<size_t>(n_rekos), 1);
		vec_opencl.d_rhs_os = vec.rhs_os.device<cl_mem>();
	}
}

// Reconsruction methods as cl_chars
void OpenCLRecMethods(const RecMethods &MethodList, RecMethodsOpenCL &MethodListOpenCL)
{
	MethodListOpenCL.MLEM = static_cast<cl_char>(MethodList.MLEM);
	MethodListOpenCL.OSEM = static_cast<cl_char>(MethodList.OSEM);
	MethodListOpenCL.RAMLA = static_cast<cl_char>(MethodList.RAMLA);
	MethodListOpenCL.MRAMLA = static_cast<cl_char>(MethodList.MRAMLA);
	MethodListOpenCL.ROSEM = static_cast<cl_char>(MethodList.ROSEM);
	MethodListOpenCL.RBI = static_cast<cl_char>(MethodList.RBI);
	MethodListOpenCL.DRAMA = static_cast<cl_char>(MethodList.DRAMA);
	MethodListOpenCL.COSEM = static_cast<cl_char>(MethodList.COSEM);
	MethodListOpenCL.ECOSEM = static_cast<cl_char>(MethodList.ECOSEM);
	MethodListOpenCL.ACOSEM = static_cast<cl_char>(MethodList.ACOSEM);

	MethodListOpenCL.MRP = static_cast<cl_char>(MethodList.MRP);
	MethodListOpenCL.Quad = static_cast<cl_char>(MethodList.Quad);
	MethodListOpenCL.Huber = static_cast<cl_char>(MethodList.Huber);
	MethodListOpenCL.L = static_cast<cl_char>(MethodList.L);
	MethodListOpenCL.FMH = static_cast<cl_char>(MethodList.FMH);
	MethodListOpenCL.WeightedMean = static_cast<cl_char>(MethodList.WeightedMean);
	MethodListOpenCL.TV = static_cast<cl_char>(MethodList.TV);
	MethodListOpenCL.AD = static_cast<cl_char>(MethodList.AD);
	MethodListOpenCL.APLS = static_cast<cl_char>(MethodList.APLS);
	MethodListOpenCL.TGV = static_cast<cl_char>(MethodList.TGV);
	MethodListOpenCL.NLM = static_cast<cl_char>(MethodList.NLM);

	MethodListOpenCL.OSLMLEM = static_cast<cl_char>(MethodList.OSLMLEM);
	MethodListOpenCL.OSLOSEM = static_cast<cl_char>(MethodList.OSLOSEM);
	MethodListOpenCL.BSREM = static_cast<cl_char>(MethodList.BSREM);
	MethodListOpenCL.MBSREM = static_cast<cl_char>(MethodList.MBSREM);
	MethodListOpenCL.ROSEMMAP = static_cast<cl_char>(MethodList.ROSEMMAP);
	MethodListOpenCL.RBIOSL = static_cast<cl_char>(MethodList.RBIOSL);
	MethodListOpenCL.OSLCOSEM = static_cast<cl_char>(MethodList.OSLCOSEM);
}

cl_int createKernels(cl_kernel & kernel_ml, cl_kernel & kernel, cl_kernel & kernel_mramla, cl_kernel& kernelNLM, const bool osem_bool, const cl_program &program_os, const cl_program& program_ml,
	const cl_program& program_mbsrem, const RecMethods MethodList, const Weighting w_vec, const uint32_t projector_type, const bool mlem_bool, const bool precompute,
	const uint16_t n_rays, const uint16_t n_rays3D)
{
	cl_int status = CL_SUCCESS;
	// Kernel for the OS-methods (OSEM, RAMLA, RBI, BSREM, etc.)
	if (osem_bool) {


		if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && ((precompute || (n_rays * n_rays3D) == 1))))
			kernel = clCreateKernel(program_os, "kernel_multi", &status);
		else
			kernel = clCreateKernel(program_os, "siddon_multi", &status);

		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			mexPrintf("Failed to create OS-methods kernel\n");
			return status;
		}
		//else if (verbose) {
			//mexPrintf("OpenCL kernel successfully created\n");
			//mexEvalString("pause(.0001);");
		//}
	}

	// Kernel for the ML-methods (MLEM, OSL-ML)
	if (mlem_bool) {

		if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1)))
			kernel_ml = clCreateKernel(program_ml, "kernel_multi", &status);
		else
			kernel_ml = clCreateKernel(program_ml, "siddon_multi", &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			mexPrintf("Failed to create MLEM kernel\n");
			return status;
		}
		//else if (verbose) {
		//	mexPrintf("OpenCL kernel successfully created\n");
		//	mexEvalString("pause(.0001);");
		//}
	}

	// Kernel for the prepass phase needed for MRAMLA, MBSREM, RBI, COSEM, ACOSEM and ECOSEM
	if ((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL) && w_vec.MBSREM_prepass ||
		MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0) {

		// Create the prepass kernel
		if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1)))
			kernel_mramla = clCreateKernel(program_mbsrem, "kernel_multi", &status);
		else
			kernel_mramla = clCreateKernel(program_mbsrem, "siddon_multi", &status);

		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			mexPrintf("Failed to create prepass kernel\n");
			return status;
		}
		else {
			mexPrintf("Prepass kernel successfully created\n");
			mexEvalString("pause(.0001);");
		}
	}
	if (MethodList.NLM) {
		if (osem_bool)
			kernelNLM = clCreateKernel(program_os, "NLM", &status);
		else if (mlem_bool)
			kernelNLM = clCreateKernel(program_ml, "NLM", &status);

		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			mexPrintf("Failed to create NLM kernel\n");
			return status;
		}
		else {
			mexPrintf("NLM kernel successfully created\n");
			mexEvalString("pause(.0001);");
		}
	}
	return status;
}

cl_int createAndWriteBuffers(cl_mem& d_x, cl_mem& d_y, cl_mem& d_z, std::vector<cl_mem>& d_lor, std::vector<cl_mem>& d_L, std::vector<cl_mem>& d_zindex,
	std::vector<cl_mem>& d_xyindex, std::vector<cl_mem>& d_Sino, std::vector<cl_mem>& d_sc_ra, const uint32_t size_x, const size_t size_z,
	const uint32_t TotSinos, const size_t size_atten, const size_t size_norm, const size_t size_scat, const uint32_t prows, std::vector<size_t>& length, const float* x, const float* y,
	const float* z_det, const uint32_t* xy_index, const uint16_t* z_index, const uint16_t* lor1, const uint16_t* L, const float* Sin, const uint8_t raw,
	cl_context& af_context, const uint32_t subsets, const uint32_t* pituus, const float* atten, const float* norm, const float* scat, const uint32_t* pseudos, const float* V,
	cl_command_queue& af_queue, cl_mem& d_atten, std::vector<cl_mem>& d_norm, std::vector<cl_mem>& d_scat, cl_mem& d_pseudos, cl_mem& d_V, cl_mem& d_xcenter, cl_mem& d_ycenter, cl_mem& d_zcenter,
	const float* x_center, const float* y_center, const float* z_center, const size_t size_center_x, const size_t size_center_y, const size_t size_center_z,
	const size_t size_of_x, const size_t size_V, const bool atomic_64bit, const bool randoms_correction, const mxArray* sc_ra, const bool precompute, cl_mem& d_lor_mlem,
	cl_mem& d_L_mlem, cl_mem& d_zindex_mlem, cl_mem& d_xyindex_mlem, cl_mem& d_Sino_mlem, cl_mem& d_sc_ra_mlem, cl_mem& d_reko_type, cl_mem& d_reko_type_mlem, const bool osem_bool,
	const bool mlem_bool, const size_t koko, const uint8_t* reko_type, const uint8_t* reko_type_mlem, const uint32_t n_rekos, const uint32_t n_rekos_mlem, cl_mem& d_norm_mlem, cl_mem& d_scat_mlem)
{
	cl_int status = CL_SUCCESS;
	// Create the necessary buffers
	// Detector coordinates
	d_V = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_V, NULL, &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	d_z = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_z, NULL, &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	d_x = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_of_x, NULL, &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	d_y = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_of_x, NULL, &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	d_xcenter = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_center_x, NULL, &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	d_ycenter = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_center_y, NULL, &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	d_zcenter = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_center_z, NULL, &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	// Attenuation data for image-based attenuation
	d_atten = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_atten, NULL, &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	// Pseudo rings
	d_pseudos = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint32_t) * prows, NULL, &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	if (osem_bool) {
		d_reko_type = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint8_t) * n_rekos, NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		for (uint32_t kk = 0; kk < subsets; kk++) {
			// How many voxels does each LOR traverse
			if (precompute)
				d_lor[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk], NULL, &status);
			else
				d_lor[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			if (size_norm > 1) {
				d_norm[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk], NULL, &status);
			}
			else {
				d_norm[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_norm, NULL, &status);
			}
			if (size_scat > 1) {
				d_scat[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk], NULL, &status);
			}
			else {
				d_scat[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_scat, NULL, &status);
			}
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			// Measurement data
			d_Sino[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk], NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			if (randoms_correction == 1u)
				d_sc_ra[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk], NULL, &status);
			else
				d_sc_ra[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			// Indices corresponding to the detector index (Sinogram data) or the detector number (raw data) at each measurement
			if (raw) {
				d_xyindex[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint32_t), NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				d_zindex[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				d_L[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk] * 2, NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
			}
			else {
				d_xyindex[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint32_t) * length[kk], NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				d_zindex[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk], NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				d_L[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
			}
		}
	}
	if (mlem_bool) {
		d_reko_type_mlem = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint8_t) * n_rekos_mlem, NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		if (precompute)
			d_lor_mlem = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t) * koko, NULL, &status);
		else
			d_lor_mlem = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		// Measurement data
		d_Sino_mlem = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * koko, NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		d_norm_mlem = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_norm, NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		d_scat_mlem = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_scat, NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		if (randoms_correction == 1u)
			d_sc_ra_mlem = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * koko, NULL, &status);
		else
			d_sc_ra_mlem = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		// Indices corresponding to the detector index (Sinogram data) or the detector number (raw data) at each measurement
		if (raw) {
			d_xyindex_mlem = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint32_t), NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			d_zindex_mlem = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			d_L_mlem = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t) * koko * 2, NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
		}
		else {
			d_xyindex_mlem = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint32_t) * koko, NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			d_zindex_mlem = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t) * koko, NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			d_L_mlem = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
		}
	}

	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		mexPrintf("Buffer creation failed\n");
		return status;
	}
	//else {
	//	mexPrintf("Buffer creation succeeded\n");
	//	mexEvalString("pause(.0001);");
	//}


	// assign values to the buffers
	status = clEnqueueWriteBuffer(af_queue, d_V, CL_FALSE, 0, sizeof(float) * size_V, V, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = clEnqueueWriteBuffer(af_queue, d_x, CL_FALSE, 0, sizeof(float) * size_of_x, x, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = clEnqueueWriteBuffer(af_queue, d_y, CL_FALSE, 0, sizeof(float) * size_of_x, y, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = clEnqueueWriteBuffer(af_queue, d_z, CL_FALSE, 0, sizeof(float) * size_z, z_det, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = clEnqueueWriteBuffer(af_queue, d_xcenter, CL_FALSE, 0, sizeof(float) * size_center_x, x_center, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = clEnqueueWriteBuffer(af_queue, d_ycenter, CL_FALSE, 0, sizeof(float) * size_center_y, y_center, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = clEnqueueWriteBuffer(af_queue, d_zcenter, CL_FALSE, 0, sizeof(float) * size_center_z, z_center, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = clEnqueueWriteBuffer(af_queue, d_atten, CL_FALSE, 0, sizeof(float) * size_atten, atten, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	status = clEnqueueWriteBuffer(af_queue, d_pseudos, CL_FALSE, 0, sizeof(uint32_t) * prows, pseudos, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	if (osem_bool) {
		status = clEnqueueWriteBuffer(af_queue, d_reko_type, CL_FALSE, 0, sizeof(uint8_t) * n_rekos, reko_type, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		for (uint32_t kk = 0; kk < subsets; kk++) {
			if (raw) {
				status = clEnqueueWriteBuffer(af_queue, d_xyindex[kk], CL_FALSE, 0, sizeof(uint32_t), xy_index, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				status = clEnqueueWriteBuffer(af_queue, d_zindex[kk], CL_FALSE, 0, sizeof(uint16_t), z_index, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				status = clEnqueueWriteBuffer(af_queue, d_L[kk], CL_FALSE, 0, sizeof(uint16_t) * length[kk] * 2, &L[pituus[kk] * 2], 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
			}
			else {
				status = clEnqueueWriteBuffer(af_queue, d_xyindex[kk], CL_FALSE, 0, sizeof(uint32_t) * length[kk], &xy_index[pituus[kk]], 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				status = clEnqueueWriteBuffer(af_queue, d_zindex[kk], CL_FALSE, 0, sizeof(uint16_t) * length[kk], &z_index[pituus[kk]], 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				status = clEnqueueWriteBuffer(af_queue, d_L[kk], CL_FALSE, 0, sizeof(uint16_t), L, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
			}
			if (precompute)
				status = clEnqueueWriteBuffer(af_queue, d_lor[kk], CL_FALSE, 0, sizeof(uint16_t) * length[kk], &lor1[pituus[kk]], 0, NULL, NULL);
			else
				status = clEnqueueWriteBuffer(af_queue, d_lor[kk], CL_FALSE, 0, sizeof(uint16_t), lor1, 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			if (size_norm > 1) {
				status = clEnqueueWriteBuffer(af_queue, d_norm[kk], CL_FALSE, 0, sizeof(float) * length[kk], &norm[pituus[kk]], 0, NULL, NULL);
			}
			else {
				status = clEnqueueWriteBuffer(af_queue, d_norm[kk], CL_FALSE, 0, sizeof(float) * size_norm, norm, 0, NULL, NULL);
			}
			if (size_scat > 1) {
				status = clEnqueueWriteBuffer(af_queue, d_scat[kk], CL_FALSE, 0, sizeof(float) * length[kk], &scat[pituus[kk]], 0, NULL, NULL);
			}
			else {
				status = clEnqueueWriteBuffer(af_queue, d_scat[kk], CL_FALSE, 0, sizeof(float) * size_scat, scat, 0, NULL, NULL);
			}
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}

			status = clEnqueueWriteBuffer(af_queue, d_Sino[kk], CL_FALSE, 0, sizeof(float) * length[kk], &Sin[pituus[kk]], 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			float* apu = (float*)mxGetData(mxGetCell(sc_ra, 0));
			if (randoms_correction)
				status = clEnqueueWriteBuffer(af_queue, d_sc_ra[kk], CL_FALSE, 0, sizeof(float) * length[kk], &apu[pituus[kk]], 0, NULL, NULL);
			else
				status = clEnqueueWriteBuffer(af_queue, d_sc_ra[kk], CL_FALSE, 0, sizeof(float), apu, 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
		}
	}
	if (mlem_bool) {
		status = clFinish(af_queue);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		status = clEnqueueWriteBuffer(af_queue, d_reko_type_mlem, CL_FALSE, 0, sizeof(uint8_t) * n_rekos_mlem, reko_type_mlem, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		if (precompute)
			status = clEnqueueWriteBuffer(af_queue, d_lor_mlem, CL_FALSE, 0, sizeof(uint16_t) * koko, lor1, 0, NULL, NULL);
		else
			status = clEnqueueWriteBuffer(af_queue, d_lor_mlem, CL_FALSE, 0, sizeof(uint16_t), lor1, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}

		status = clEnqueueWriteBuffer(af_queue, d_Sino_mlem, CL_FALSE, 0, sizeof(float) * koko, Sin, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		status = clEnqueueWriteBuffer(af_queue, d_norm_mlem, CL_FALSE, 0, sizeof(float) * size_norm, norm, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		status = clEnqueueWriteBuffer(af_queue, d_scat_mlem, CL_FALSE, 0, sizeof(float) * size_scat, scat, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		float* apu = (float*)mxGetData(mxGetCell(sc_ra, 0));
		if (randoms_correction)
			status = clEnqueueWriteBuffer(af_queue, d_sc_ra_mlem, CL_FALSE, 0, sizeof(float) * koko, apu, 0, NULL, NULL);
		else
			status = clEnqueueWriteBuffer(af_queue, d_sc_ra_mlem, CL_FALSE, 0, sizeof(float), apu, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		if (raw) {
			status = clEnqueueWriteBuffer(af_queue, d_xyindex_mlem, CL_FALSE, 0, sizeof(uint32_t), xy_index, 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			status = clEnqueueWriteBuffer(af_queue, d_zindex_mlem, CL_FALSE, 0, sizeof(uint16_t), z_index, 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			status = clEnqueueWriteBuffer(af_queue, d_L_mlem, CL_FALSE, 0, sizeof(uint16_t) * koko * 2ULL, L, 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
		}
		else {
			status = clEnqueueWriteBuffer(af_queue, d_xyindex_mlem, CL_FALSE, 0, sizeof(uint32_t) * koko, xy_index, 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			status = clEnqueueWriteBuffer(af_queue, d_zindex_mlem, CL_FALSE, 0, sizeof(uint16_t) * koko, z_index, 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
			status = clEnqueueWriteBuffer(af_queue, d_L_mlem, CL_FALSE, 0, sizeof(uint16_t), L, 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return status;
			}
		}
	}

	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		mexPrintf("Buffer write failed\n");
		return status;
	}
	//else {
	//	mexPrintf("Buffer write succeeded\n");
	//	mexEvalString("pause(.0001);");
	//}
	return status;
}

// Prepass phase for MRAMLA, COSEM, ACOSEM, ECOSEM
void MRAMLA_prepass(const uint32_t subsets, const uint32_t im_dim, const uint32_t* pituus, const std::vector<cl_mem> &lor, const std::vector<cl_mem> &zindex,
	const std::vector<cl_mem> &xindex, cl_program program, const cl_command_queue &af_queue, const cl_context af_context, Weighting& w_vec, 
	std::vector<af::array>& Summ, const std::vector<cl_mem> &d_Sino, const size_t koko_l, af::array& cosem, af::array& C_co, 
	af::array& C_aco, af::array& C_osl, const uint32_t alku, cl_kernel &kernel_mramla, const std::vector<cl_mem> &L, const uint8_t raw, 
	const RecMethodsOpenCL MethodListOpenCL, const std::vector<size_t> length, const bool atomic_64bit, const cl_uchar compute_norm_matrix, 
	const std::vector<cl_mem>& d_sc_ra, cl_uint kernelInd_MRAMLA, af::array& E, const std::vector<cl_mem>& d_norm, const std::vector<cl_mem>& d_scat, const bool use_psf,
	const af::array& g, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float epps) {

	cl_int status = CL_SUCCESS;

	cl_uchar MBSREM_prepass = static_cast<cl_uchar>(w_vec.MBSREM_prepass);

	af::array apu_co, apu_aco, uu, sub_index_array, apu_summa, apu_Amin, apu_summa_m;

	af::array cosem_psf;
	if (use_psf) {
		cosem_psf = computeConvolution(cosem, g, Nx, Ny, Nz, w_vec, 1u);
		af::sync();
	}
	else
		cosem_psf = cosem;

	cl_mem * d_cosem = cosem_psf.device<cl_mem>();

	af::array apu;

	uint32_t eka = alku;

	if (alku > 0u) {
		eka = alku - 1u;
		apu = af::constant(0.f, length[eka]);
	}
	else {
		apu = af::constant(0.f, 1);
	}
	cl_mem * d_ACOSEM_lhs = apu.device<cl_mem>();

	const size_t local_size = 64ULL;

	for (uint32_t osa_iter = eka; osa_iter < subsets; osa_iter++) {

		cl_uint kernelInd_MRAMLA_sub = kernelInd_MRAMLA;

		size_t erotus = length[osa_iter] % local_size;

		if (erotus > 0)
			erotus = (local_size - erotus);

		const size_t global_size = length[osa_iter] + erotus;

		const uint64_t m_size = static_cast<uint64_t>(length[osa_iter]);

		if (subsets > 1) {


			if ((MethodListOpenCL.MRAMLA || MethodListOpenCL.MBSREM) && w_vec.MBSREM_prepass && alku == 0u) {
				apu_Amin = af::constant(0.f, length[osa_iter], 1);
				apu_summa_m = af::constant(0.f, length[osa_iter], 1);
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

			if ((MethodListOpenCL.COSEM || MethodListOpenCL.ECOSEM || MethodListOpenCL.OSLCOSEM == 2) && alku == 0u) {
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

			if ((MethodListOpenCL.ACOSEM || MethodListOpenCL.OSLCOSEM == 1) && alku == 0u) {
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

			if ((MethodListOpenCL.COSEM || MethodListOpenCL.ECOSEM || MethodListOpenCL.OSLCOSEM == 2) && alku == 0u) {
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

			if ((MethodListOpenCL.ACOSEM || MethodListOpenCL.OSLCOSEM == 1) && alku == 0) {
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

			if ((MethodListOpenCL.MRAMLA || MethodListOpenCL.MBSREM) && w_vec.MBSREM_prepass && alku == 0u)
				apu_Amin = w_vec.Amin;
			else
				apu_Amin = af::constant(0.f, 1, 1);
		}

		//af::sync();

		cl_mem * d_apu_co = apu_co.device<cl_mem>();
		cl_mem * d_apu_aco = apu_aco.device<cl_mem>();
		cl_mem * d_E = apu_summa_m.device<cl_mem>();
		cl_mem * d_Amin = apu_Amin.device<cl_mem>();
		cl_mem* d_Summ = apu_summa.device<cl_mem>();

		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA_sub++, sizeof(cl_mem), &d_norm[osa_iter]);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA_sub++, sizeof(cl_mem), &d_scat[osa_iter]);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA_sub++, sizeof(cl_mem), d_Summ);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA_sub++, sizeof(cl_mem), &lor[osa_iter]);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA_sub++, sizeof(cl_mem), &xindex[osa_iter]);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA_sub++, sizeof(cl_mem), &zindex[osa_iter]);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA_sub++, sizeof(cl_mem), &L[osa_iter]);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA_sub++, sizeof(cl_mem), &d_Sino[osa_iter]);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA_sub++, sizeof(cl_mem), &d_sc_ra[osa_iter]);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA_sub++, sizeof(cl_mem), d_cosem);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA_sub++, sizeof(uint32_t), &alku);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA_sub++, sizeof(cl_uchar), &MBSREM_prepass);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA_sub++, sizeof(cl_mem), d_ACOSEM_lhs);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA_sub++, sizeof(cl_mem), d_Amin);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA_sub++, sizeof(cl_mem), d_apu_co);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA_sub++, sizeof(cl_mem), d_apu_aco);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA_sub++, sizeof(cl_mem), d_E);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA_sub++, sizeof(uint64_t), &m_size);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA_sub++, sizeof(MethodListOpenCL), &MethodListOpenCL);
		status = clEnqueueNDRangeKernel(af_queue, kernel_mramla, 1u, NULL, &global_size, &local_size, 0u, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			mexPrintf("Failed to launch the prepass kernel\n");
			return;
		}

		clFinish(af_queue);
		apu_co.unlock();
		apu_aco.unlock();
		apu_summa.unlock();
		apu_summa_m.unlock();
		apu.unlock();
		apu_Amin.unlock();
		cosem_psf.unlock();

		if (alku == 0u) {
			if ((MethodListOpenCL.COSEM || MethodListOpenCL.ECOSEM)) {
				if (atomic_64bit)
					C_co(af::span, osa_iter) = apu_co.as(f32) / TH;
				else
					C_co(af::span, osa_iter) = apu_co;
				if (use_psf)
					C_co(af::span, osa_iter) = computeConvolution(C_co(af::span, osa_iter), g, Nx, Ny, Nz, w_vec, 1u) * cosem;
				else
					C_co(af::span, osa_iter) = C_co(af::span, osa_iter) * cosem;
				//mexPrintf("co = %f\n", af::sum<float>(C_co(af::span, osa_iter)));
				//mexPrintf("dim0 = %u\n", C_co(af::span, osa_iter).dims(0));
				//mexPrintf("dim1 = %u\n", C_co(af::span, osa_iter).dims(1));
			}
			if (MethodListOpenCL.ACOSEM) {
				if (atomic_64bit)
					C_aco(af::span, osa_iter) = apu_aco.as(f32) / TH;
				else
					C_aco(af::span, osa_iter) = apu_aco;
				if (use_psf)
					C_aco(af::span, osa_iter) = computeConvolution(C_aco(af::span, osa_iter), g, Nx, Ny, Nz, w_vec, 1u) * af::pow(cosem, w_vec.h_ACOSEM_2);
				else
					C_aco(af::span, osa_iter) = C_aco(af::span, osa_iter) * af::pow(cosem, w_vec.h_ACOSEM_2);
			}
			if (MethodListOpenCL.OSLCOSEM == 2u) {
				if (atomic_64bit)
					C_osl(af::span, osa_iter) = apu_co.as(f32) / TH;
				else
					C_osl(af::span, osa_iter) = apu_co;
				if (use_psf)
					C_osl(af::span, osa_iter) = computeConvolution(C_osl(af::span, osa_iter), g, Nx, Ny, Nz, w_vec, 1u) * cosem;
				else
					C_osl(af::span, osa_iter) = C_osl(af::span, osa_iter) * cosem;
			}
			else if (MethodListOpenCL.OSLCOSEM == 1) {
				if (atomic_64bit)
					C_osl(af::span, osa_iter) = apu_aco.as(f32) / TH;
				else
					C_osl(af::span, osa_iter) = apu_aco;
				if (use_psf)
					C_osl(af::span, osa_iter) = computeConvolution(C_osl(af::span, osa_iter), g, Nx, Ny, Nz, w_vec, 1u) * af::pow(cosem, w_vec.h_ACOSEM_2);
				else
					C_osl(af::span, osa_iter) = C_osl(af::span, osa_iter) * af::pow(cosem, w_vec.h_ACOSEM_2);
			}


			if ((MethodListOpenCL.MRAMLA || MethodListOpenCL.MBSREM) && w_vec.MBSREM_prepass) {
				sub_index_array = af::range(af::dim4(length[osa_iter]), 0, u32) + pituus[osa_iter];
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

				if ((MethodListOpenCL.MRAMLA || MethodListOpenCL.MBSREM) && w_vec.MBSREM_prepass) {
					if (subsets > 1)
						E(sub_index_array) = apu_summa_m;
					else
						E = apu_summa_m;
				}
			}
		}
		else {
			w_vec.ACOSEM_rhs = af::sum<float>(apu);
		}
	}
	if (use_psf && alku == 0 && w_vec.MBSREM_prepass) {
		w_vec.D = computeConvolution(w_vec.D, g, Nx, Ny, Nz, w_vec, 1u);
	}
	return;
}

void precomp_siddon(const cl_context &context, const cl_command_queue &commandQueues,
	uint16_t* lor1, const float* z_det, const float* x, const float* y, const uint32_t Nx, const uint32_t Ny,
	const uint32_t Nz, const float dx, const float dy, const float dz, const float bx, const float by, const float bz,
	const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices,
	const uint32_t size_x, const uint16_t TotSinos, const bool verbose, const uint32_t loop_var_par,
	const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L, const uint8_t raw, const size_t size_z, const uint32_t im_dim,
	const cl_kernel &kernel, const size_t numel_x) {

	cl_int status = CL_SUCCESS;
	const uint32_t Nxy = Nx * Ny;

	size_t osa_length = loop_var_par;

	cl_mem d_z, d_x, d_y, d_pseudos, d_L, d_lor;


	// Create the necessary buffers
	d_z = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_z, NULL, &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return;
	}
	d_x = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * numel_x, NULL, &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return;
	}
	d_y = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * numel_x, NULL, &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return;
	}
	d_pseudos = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t) * prows, NULL, &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return;
	}
	d_lor = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(uint16_t) * osa_length, NULL, &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return;
	}
	if (raw) {
		d_L = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * osa_length * 2, NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
	}
	else {
		d_L = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
	}

	status = clEnqueueWriteBuffer(commandQueues, d_x, CL_FALSE, 0, sizeof(float) * numel_x, x, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return;
	}
	status = clEnqueueWriteBuffer(commandQueues, d_y, CL_FALSE, 0, sizeof(float) * numel_x, y, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return;
	}
	status = clEnqueueWriteBuffer(commandQueues, d_z, CL_FALSE, 0, sizeof(float) * size_z, z_det, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return;
	}
	status = clEnqueueWriteBuffer(commandQueues, d_pseudos, CL_FALSE, 0, sizeof(uint32_t) * prows, pseudos, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return;
	}

	if (raw) {
		status = clEnqueueWriteBuffer(commandQueues, d_L, CL_FALSE, 0, sizeof(uint16_t) * osa_length * 2, L, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
	}
	else {
		status = clEnqueueWriteBuffer(commandQueues, d_L, CL_FALSE, 0, sizeof(uint16_t), L, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
	}


	clSetKernelArg(kernel, 0, sizeof(uint32_t), &Nxy);
	clSetKernelArg(kernel, 1, sizeof(uint32_t), &im_dim);
	clSetKernelArg(kernel, 2, sizeof(uint32_t), &Nx);
	clSetKernelArg(kernel, 3, sizeof(uint32_t), &Ny);
	clSetKernelArg(kernel, 4, sizeof(uint32_t), &Nz);
	clSetKernelArg(kernel, 5, sizeof(float), &dz);
	clSetKernelArg(kernel, 6, sizeof(float), &dx);
	clSetKernelArg(kernel, 7, sizeof(float), &dy);
	clSetKernelArg(kernel, 8, sizeof(float), &bz);
	clSetKernelArg(kernel, 9, sizeof(float), &bx);
	clSetKernelArg(kernel, 10, sizeof(float), &by);
	clSetKernelArg(kernel, 11, sizeof(float), &bzb);
	clSetKernelArg(kernel, 12, sizeof(float), &maxxx);
	clSetKernelArg(kernel, 13, sizeof(float), &maxyy);
	clSetKernelArg(kernel, 14, sizeof(float), &zmax);
	clSetKernelArg(kernel, 15, sizeof(float), &NSlices);
	clSetKernelArg(kernel, 16, sizeof(uint32_t), &size_x);
	clSetKernelArg(kernel, 17, sizeof(uint16_t), &TotSinos);
	clSetKernelArg(kernel, 18, sizeof(uint32_t), &det_per_ring);
	clSetKernelArg(kernel, 19, sizeof(uint8_t), &raw);
	clSetKernelArg(kernel, 20, sizeof(uint32_t), &prows);

	cl_event event1;

	clSetKernelArg(kernel, 21, sizeof(cl_mem), &d_pseudos);
	clSetKernelArg(kernel, 22, sizeof(cl_mem), &d_x);
	clSetKernelArg(kernel, 23, sizeof(cl_mem), &d_y);
	clSetKernelArg(kernel, 24, sizeof(cl_mem), &d_z);
	clSetKernelArg(kernel, 25, sizeof(cl_mem), &d_lor);
	clSetKernelArg(kernel, 26, sizeof(cl_mem), &d_L);
	status = clEnqueueNDRangeKernel(commandQueues, kernel, 1, NULL, &osa_length, NULL, 0, NULL, &event1);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;

		clReleaseEvent(event1);

		clReleaseMemObject(d_z);
		clReleaseMemObject(d_x);
		clReleaseMemObject(d_y);
		clReleaseMemObject(d_pseudos);
		clReleaseMemObject(d_L);
		clReleaseMemObject(d_lor);
		return;
	}

	status = clEnqueueReadBuffer(commandQueues, d_lor, CL_TRUE, 0, sizeof(uint16_t) * osa_length, lor1, 1, &event1, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;

		clReleaseEvent(event1);

		clReleaseMemObject(d_z);
		clReleaseMemObject(d_x);
		clReleaseMemObject(d_y);
		clReleaseMemObject(d_pseudos);
		clReleaseMemObject(d_L);
		clReleaseMemObject(d_lor);
		return;
	}

	clReleaseEvent(event1);

	clReleaseMemObject(d_z);
	clReleaseMemObject(d_x);
	clReleaseMemObject(d_y);
	clReleaseMemObject(d_pseudos);
	clReleaseMemObject(d_L);
	clReleaseMemObject(d_lor);

	af::sync();
	return;
}

cl_int createProgram(const bool verbose, const char* k_path, cl_context& af_context, cl_device_id& af_device_id, const char* fileName, 
	cl_program& program_os, cl_program& program_ml, cl_program& program_mbsrem, bool& atomic_64bit, const uint32_t device, const char* header_directory,
	const uint32_t projector_type, const float crystal_size_z, const bool precompute, const uint8_t raw, const uint32_t attenuation_correction, 
	const uint32_t normalization_correction, const int32_t dec, const size_t local_size, const uint16_t n_rays, const uint16_t n_rays3D, 
	const bool find_lors, const RecMethods MethodList, const bool osem_bool, const bool mlem_bool, const uint32_t n_rekos, const uint32_t n_rekos_mlem, 
	const Weighting& w_vec, const uint32_t osa_iter0, const float cr_pz, const float dx, const bool use_psf, const uint32_t scatter, const uint32_t randoms_correction) {

	cl_int status = CL_SUCCESS;

	std::string options = header_directory;
	options += " -cl-single-precision-constant";

	if (crystal_size_z == 0.f && projector_type == 2u)
		options += " -DCRYST";
	if ((crystal_size_z > 0.f && projector_type == 2u) || projector_type == 3u)
		options += " -DCRYSTZ";
	if (projector_type == 3u)
		options += " -DVOL";
	if (precompute)
		options += " -DPRECOMPUTE";
	if (raw == 1)
		options += " -DRAW";
	if (projector_type == 1u)
		options += " -DSIDDON";
	else if (projector_type == 2u || projector_type == 3u)
		options += " -DORTH";
	if (attenuation_correction == 1u)
		options += " -DATN";
	if (normalization_correction == 1u)
		options += " -DNORM";
	if (scatter == 1u)
		options += " -DSCATTER";
	if (randoms_correction == 1u)
		options += " -DRANDOMS";
	options += " -DFP";
	if (projector_type == 1u && !precompute && (n_rays * n_rays3D) > 1) {
		options += (" -DN_RAYS=" + std::to_string(n_rays * n_rays3D));
		options += (" -DN_RAYS2D=" + std::to_string(n_rays));
		options += (" -DN_RAYS3D=" + std::to_string(n_rays3D));
	}
	if (find_lors)
		options += " -DFIND_LORS";
	//if (projector_type == 1u && use_psf && (precompute || (n_rays * n_rays3D) == 1)) {
	//	options += " -DORTH";
	//	options += " -DCRYSTZ";
	//	uint32_t limit = static_cast<uint32_t>(std::floor(cr_pz / dx));
	//	options += (" -DPSF_LIMIT=" + std::to_string(limit));
	//	options += (" -DX=" + std::to_string(dx));
	//	options += (" -DSIGMA=" + std::to_string(cr_pz / 2.355f));
	//}
	options += (" -DLOCAL_SIZE=" + std::to_string(local_size));
	if (MethodList.NLM) {
		options += " -DNLM_";
	}
	if (osem_bool) {
		std::string os_options = options;
		if ((projector_type == 2 || projector_type == 3u) && dec > 0)
			os_options += (" -DDEC=" + std::to_string(dec));
		os_options += (" -DN_REKOS=" + std::to_string(n_rekos));
		os_options += " -DAF";
		if (n_rekos == 1)
			os_options += " -DNREKOS1";
		else if (n_rekos == 2)
			os_options += " -DNREKOS2";
		if (MethodList.MRAMLA || MethodList.MBSREM)
			os_options += " -DMRAMLA";
		if (MethodList.COSEM || MethodList.ACOSEM || MethodList.OSLCOSEM > 0u || MethodList.ECOSEM)
			os_options += " -DCOSEM";

		status = buildProgram(verbose, k_path, af_context, af_device_id, program_os, atomic_64bit, os_options);
	}
	if (mlem_bool) {
		std::string ml_options = options;
		if ((projector_type == 2 || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) > 2))) && dec > 0)
			ml_options += (" -DDEC=" + std::to_string(dec));
		ml_options += (" -DN_REKOS=" + std::to_string(n_rekos_mlem));
		ml_options += " -DAF";
		if (n_rekos_mlem == 1)
			ml_options += " -DNREKOS1";
		else if (n_rekos_mlem == 2)
			ml_options += " -DNREKOS2";

		status = buildProgram(verbose, k_path, af_context, af_device_id, program_ml, atomic_64bit, ml_options);
	}
	if ((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL) && w_vec.MBSREM_prepass ||
		MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0) {
		options += " -DAF";
		options += " -DMBSREM";
		if (n_rekos == 1)
			options += " -DNREKOS1";
		else if (n_rekos == 2)
			options += " -DNREKOS2";
		if (MethodList.MRAMLA || MethodList.MBSREM)
			options += " -DMRAMLA";
		status = buildProgram(verbose, k_path, af_context, af_device_id, program_mbsrem, atomic_64bit, options);
	}
	return status;
}

cl_int buildProgram(const bool verbose, const char* k_path, cl_context& af_context, cl_device_id& af_device_id, cl_program& program, 
	bool& atomic_64bit, std::string options) {
	cl_int status = CL_SUCCESS;
	size_t pituus;
	if (atomic_64bit) {
		pituus = options.length();
		options += " -DCAST=ulong";
		options += " -DATOMIC";
		options += (" -DTH=" + std::to_string(TH));
	}
	else
		options += " -DCAST=float";
	//mexPrintf("%s\n", options.c_str());
	if (atomic_64bit) {
		std::string kernel_path_atom;

		kernel_path_atom = k_path;
		kernel_path_atom += ".cl";
		// Load the source text file
		std::fstream sourceFile_atom(kernel_path_atom.c_str());
		std::string content_atom((std::istreambuf_iterator<char>(sourceFile_atom)), std::istreambuf_iterator<char>());
		const char* sourceCode_atom = new char[content_atom.size()];
		sourceCode_atom = content_atom.c_str();
		// Create the program from the source
		program = clCreateProgramWithSource(af_context, 1, (const char**)&sourceCode_atom, NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}

		// Build the program
		status = clBuildProgram(program, 1, &af_device_id, options.c_str(), NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			mexPrintf("Failed to build OpenCL program. Build log: \n");
			size_t len;
			char* buffer;
			clGetProgramBuildInfo(program, af_device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
			buffer = (char*)calloc(len, sizeof(size_t));
			clGetProgramBuildInfo(program, af_device_id, CL_PROGRAM_BUILD_LOG, len, buffer, NULL);
			mexPrintf("%s\n", buffer);
			////return status;
			options.erase(pituus, options.size() + 1);
			options += " -DCAST=float";
			//mexPrintf("%s\n", options.c_str());
			mexPrintf("Failed to build 64-bit atomics program.\n");
		}
		else if (verbose)
			mexPrintf("OpenCL program (64-bit atomics) built\n");
	}
	else
		status = -1;
	// If not, use 32-bit atomic add (float)
	if (status != CL_SUCCESS) {
		status = CL_SUCCESS;
		atomic_64bit = false;

		std::string kernel_path;

		kernel_path = k_path;
		kernel_path += ".cl";
		std::fstream sourceFile(kernel_path.c_str());
		std::string content((std::istreambuf_iterator<char>(sourceFile)), std::istreambuf_iterator<char>());
		const char* sourceCode = new char[content.size()];
		sourceCode = content.c_str();
		program = clCreateProgramWithSource(af_context, 1, (const char**)&sourceCode, NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		//mexPrintf("%s\n", options.c_str());
		status = clBuildProgram(program, 1, &af_device_id, options.c_str(), NULL, NULL);
		// Build log in case of failure
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			mexPrintf("Failed to build OpenCL program. Build log: \n");
			size_t len;
			char* buffer;
			clGetProgramBuildInfo(program, af_device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
			buffer = (char*)calloc(len, sizeof(size_t));
			clGetProgramBuildInfo(program, af_device_id, CL_PROGRAM_BUILD_LOG, len, buffer, NULL);
			mexPrintf("%s\n", buffer);
			free(buffer);
			return status;
		}
		else if (verbose)
			mexPrintf("OpenCL program built\n");
	}
	return status;
}

af::array NLM(const af::array& im, Weighting& w_vec, const float epps, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const kernelStruct& OpenCLStruct)
{
	cl_int status = CL_SUCCESS;
	const int32_t ndx = static_cast<int32_t>(w_vec.Ndx);
	const int32_t ndy = static_cast<int32_t>(w_vec.Ndy);
	const int32_t ndz = static_cast<int32_t>(w_vec.Ndz);
	const int32_t nlx = static_cast<int32_t>(w_vec.Nlx);
	const int32_t nly = static_cast<int32_t>(w_vec.Nly);
	const int32_t nlz = static_cast<int32_t>(w_vec.Nlz);
	size_t kernelIndNLM = 0ULL;
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
	//std::string joku = af::toString("gaussianNLM", w_vec.gaussianNLM);
	//mexPrintf("%s\n", joku.c_str());
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
	cl_mem* d_W = W.device<cl_mem>();
	cl_mem* d_input = input.device<cl_mem>();
	cl_mem* d_padInput = padInput.device<cl_mem>();
	cl_mem* d_gaussianNLM = w_vec.gaussianNLM.device<cl_mem>();
	//size_t global_size[] = { N, M, K };
	size_t global_size = N * M * K;
	//size_t global_offset[] = { window_z , 0, 0 };
	af::sync();
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(cl_mem), d_W);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(cl_mem), d_input);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(cl_mem), d_padInput);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(cl_mem), d_gaussianNLM);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(int32_t), &ndx);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(int32_t), &ndy);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(int32_t), &ndz);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(int32_t), &nlx);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(int32_t), &nly);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(int32_t), &nlz);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(uint32_t), &N);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(uint32_t), &M);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(uint32_t), &K);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(float), &w_vec.h2);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(float), &epps);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(int32_t), &Nxy);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(int32_t), &min_x);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(int32_t), &max_x);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(int32_t), &min_y);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(int32_t), &max_y);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(int32_t), &min_z);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(int32_t), &max_z);
	clSetKernelArg(OpenCLStruct.kernelNLM, kernelIndNLM++, sizeof(int32_t), &type);
	// Compute the kernel
	status = clEnqueueNDRangeKernel(*OpenCLStruct.af_queue, OpenCLStruct.kernelNLM, 1u, NULL, &global_size, NULL, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		mexPrintf("Failed to launch the NLM kernel\n");
		mexEvalString("pause(.0001);");
	}
	status = clFinish(*OpenCLStruct.af_queue);
	if (status != CL_SUCCESS) {
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