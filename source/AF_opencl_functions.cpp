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
void update_opencl_inputs(AF_im_vectors& vec, OpenCL_im_vectors& vec_opencl, const bool mlem, const uint32_t im_dim, const uint32_t n_rekos,
	const uint32_t n_rekos_mlem, const RecMethods MethodList, const bool atomic_64bit, const bool atomic_32bit, const bool use_psf, const uint32_t nMAPOS)
{

	if (mlem) {
		if (use_psf)
			vec_opencl.d_im_mlem = cl::Buffer(*vec.im_mlem_blurred.device<cl_mem>(), true);
		else
			vec_opencl.d_im_mlem = cl::Buffer(*vec.im_mlem.device<cl_mem>(), true);
		if (atomic_64bit)
			vec.rhs_mlem = af::constant(0LL, static_cast<size_t>(im_dim) * n_rekos_mlem, 1, s64);
		else if (atomic_32bit)
			vec.rhs_mlem = af::constant(0, static_cast<size_t>(im_dim) * n_rekos_mlem, 1, s32);
		else
			vec.rhs_mlem = af::constant(0.f, static_cast<size_t>(im_dim) * n_rekos_mlem, 1);
		vec_opencl.d_rhs_mlem = cl::Buffer(*vec.rhs_mlem.device<cl_mem>(), true);
	}
	else {
		if (use_psf)
			vec_opencl.d_im_os = cl::Buffer(*vec.im_os_blurred.device<cl_mem>(), true);
		else
			vec_opencl.d_im_os = cl::Buffer(*vec.im_os.device<cl_mem>(), true);
		if (atomic_64bit)
			vec.rhs_os = af::constant(0LL, static_cast<size_t>(im_dim) * static_cast<size_t>(n_rekos), 1, s64);
		else if (atomic_32bit)
			vec.rhs_os = af::constant(0, static_cast<size_t>(im_dim) * static_cast<size_t>(n_rekos), 1, s32);
		else
			vec.rhs_os = af::constant(0.f, static_cast<size_t>(im_dim) * static_cast<size_t>(n_rekos), 1);
		vec_opencl.d_rhs_os = cl::Buffer(*vec.rhs_os.device<cl_mem>(), true);
	}
}

// Reconsruction methods as cl_chars
void OpenCLRecMethods(const RecMethods &MethodList, RecMethodsOpenCL &MethodListOpenCL)
{
	// Non-MAP/prior algorithms
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

	// Priors
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

	// MAP/prior-based algorithms
	MethodListOpenCL.OSLMLEM = static_cast<cl_char>(MethodList.OSLMLEM);
	MethodListOpenCL.OSLOSEM = static_cast<cl_char>(MethodList.OSLOSEM);
	MethodListOpenCL.BSREM = static_cast<cl_char>(MethodList.BSREM);
	MethodListOpenCL.MBSREM = static_cast<cl_char>(MethodList.MBSREM);
	MethodListOpenCL.ROSEMMAP = static_cast<cl_char>(MethodList.ROSEMMAP);
	MethodListOpenCL.RBIOSL = static_cast<cl_char>(MethodList.RBIOSL);
	MethodListOpenCL.OSLCOSEM = static_cast<cl_char>(MethodList.OSLCOSEM);
	MethodListOpenCL.PKMA = static_cast<cl_char>(MethodList.PKMA);
}

cl_int createKernels(cl::Kernel& kernel_ml, cl::Kernel & kernel, cl::Kernel& kernel_mramla, cl::Kernel& kernelNLM, cl::Kernel& kernelMed, const bool osem_bool, const cl::Program &program_os, const cl::Program& program_ml,
	const cl::Program& program_mbsrem, const RecMethods MethodList, const Weighting w_vec, const uint32_t projector_type, const bool mlem_bool, const bool precompute,
	const uint16_t n_rays, const uint16_t n_rays3D)
{
	cl_int status = CL_SUCCESS;
	// Kernel for the OS-methods (OSEM, RAMLA, RBI, BSREM, etc.)
	if (osem_bool) {


		if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && ((precompute || (n_rays * n_rays3D) == 1)))) {
			kernel = cl::Kernel(program_os, "kernel_multi", &status);
		}
			//kernel = clCreateKernel(program_os, "kernel_multi", &status);
		else
			kernel = cl::Kernel(program_os, "siddon_multi", &status);

		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to create OS-methods kernel\n");
			return status;
		}
		else if (DEBUG) {
			mexPrintf("OpenCL kernel successfully created\n");
			mexEvalString("pause(.0001);");
		}
	}

	// Kernel for the ML-methods (MLEM, OSL-ML)
	if (mlem_bool) {

		if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1)))
			kernel_ml = cl::Kernel(program_ml, "kernel_multi", &status);
		else
			kernel_ml = cl::Kernel(program_ml, "siddon_multi", &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to create MLEM kernel\n");
			return status;
		}
		else if (DEBUG) {
			mexPrintf("OpenCL kernel successfully created\n");
			mexEvalString("pause(.0001);");
		}
	}

	// Kernel for the prepass phase needed for MRAMLA, MBSREM, RBI, COSEM, ACOSEM and ECOSEM
	if ((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI || MethodList.PKMA) && w_vec.MBSREM_prepass ||
		MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0) {

		// Create the prepass kernel
		if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1)))
			kernel_mramla = cl::Kernel(program_mbsrem, "kernel_multi", &status);
		else
			kernel_mramla = cl::Kernel(program_mbsrem, "siddon_multi", &status);

		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to create prepass kernel\n");
			return status;
		}
		else if (DEBUG) {
			mexPrintf("Prepass kernel successfully created\n");
			mexEvalString("pause(.0001);");
		}
	}
	if (MethodList.NLM) {
		if (osem_bool)
			kernelNLM = cl::Kernel(program_os, "NLM", &status);
		else if (mlem_bool)
			kernelNLM = cl::Kernel(program_ml, "NLM", &status);

		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to create NLM kernel\n");
			return status;
		}
		else if (DEBUG) {
			mexPrintf("NLM kernel successfully created\n");
			mexEvalString("pause(.0001);");
		}
	}
	if (MethodList.MRP) {
		if (osem_bool)
			kernelMed = cl::Kernel(program_os, "medianFilter3D", &status);
		else if (mlem_bool)
			kernelMed = cl::Kernel(program_ml, "medianFilter3D", &status);

		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to create Median kernel\n");
			return status;
		}
		else if (DEBUG) {
			mexPrintf("Median kernel successfully created\n");
			mexEvalString("pause(.0001);");
		}
	}
	return status;
}

cl_int createAndWriteBuffers(cl::Buffer& d_x, cl::Buffer& d_y, cl::Buffer& d_z, cl::Buffer& d_angles, std::vector<cl::Buffer>& d_lor, std::vector<cl::Buffer>& d_L, std::vector<cl::Buffer>& d_zindex,
	std::vector<cl::Buffer>& d_xyindex, std::vector<cl::Buffer>& d_Sino, std::vector<cl::Buffer>& d_sc_ra, const uint32_t size_x, const size_t size_z,
	const uint32_t TotSinos, const size_t size_atten, const size_t size_norm, const size_t size_scat, const uint32_t prows, std::vector<size_t>& length, const float* x, const float* y,
	const float* z_det, const uint32_t* xy_index, const uint16_t* z_index, const uint16_t* lor1, const uint16_t* L, const float* Sin, const uint8_t raw,
	cl::Context& af_context, const uint32_t subsets, const int64_t* pituus, const float* atten, const float* norm, const float* scat, const uint32_t* pseudos, const float* V,
	cl::CommandQueue& af_queue, cl::Buffer& d_atten, std::vector<cl::Buffer>& d_norm, std::vector<cl::Buffer>& d_scat, cl::Buffer& d_pseudos, cl::Buffer& d_V, cl::Buffer& d_xcenter, 
	cl::Buffer& d_ycenter, cl::Buffer& d_zcenter, const float* x_center, const float* y_center, const float* z_center, const size_t size_center_x, const size_t size_center_y, 
	const size_t size_center_z, const size_t size_of_x, const size_t size_V, const bool atomic_64bit, const bool atomic_32bit, const bool randoms_correction, const mxArray* sc_ra, const bool precompute,
	cl::Buffer& d_lor_mlem, cl::Buffer& d_L_mlem, cl::Buffer& d_zindex_mlem, cl::Buffer& d_xyindex_mlem, cl::Buffer& d_Sino_mlem, cl::Buffer& d_sc_ra_mlem, cl::Buffer& d_reko_type, 
	cl::Buffer& d_reko_type_mlem, const bool osem_bool,	const bool mlem_bool, const size_t koko, const uint8_t* reko_type, const uint8_t* reko_type_mlem, const uint32_t n_rekos, 
	const uint32_t n_rekos_mlem, cl::Buffer& d_norm_mlem, cl::Buffer& d_scat_mlem, const float* angles, const bool TOF, const int64_t nBins, const bool loadTOF, cl::Buffer& d_TOFCenter, 
	const float* TOFCenter, const uint32_t subsetsUsed, const uint32_t osa_iter0, const uint8_t listmode, const bool CT)
{
	cl_int status = CL_SUCCESS;
	// Create the necessary buffers
	// Detector coordinates
	d_V = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_V, NULL, &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	d_z = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_z, NULL, &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	d_x = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_of_x, NULL, &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	d_y = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_of_x, NULL, &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	if (CT) {
		d_angles = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * TotSinos, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
	}
	d_xcenter = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_center_x, NULL, &status);;
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	d_ycenter = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_center_y, NULL, &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	d_zcenter = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_center_z, NULL, &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	// Attenuation data for image-based attenuation
	d_atten = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_atten, NULL, &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	// TOF bin centers
	d_TOFCenter = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * nBins, NULL, &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	// Pseudo rings
	d_pseudos = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint32_t) * prows, NULL, &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	if (osem_bool) {
		d_reko_type = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint8_t) * n_rekos, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
		for (uint32_t kk = osa_iter0; kk < subsetsUsed; kk++) {
			// How many voxels does each LOR traverse
			if (precompute)
				d_lor[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk], NULL, &status);
			else
				d_lor[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			if (size_norm > 1) {
				d_norm[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk], NULL, &status);
			}
			else {
				d_norm[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
			}
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			if (size_scat > 1) {
				d_scat[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk], NULL, &status);
			}
			else {
				d_scat[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
			}
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			// Measurement data
			if (TOF) {
				if (loadTOF) {
					d_Sino[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * nBins, NULL, &status);
				}
				else {
					if (kk == 0)
						d_Sino[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * nBins, NULL, &status);
				}
			}
			else
				d_Sino[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk], NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			if (randoms_correction == 1u)
				d_sc_ra[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk], NULL, &status);
			else
				d_sc_ra[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			// Indices corresponding to the detector index (Sinogram data) or the detector number (raw data) at each measurement
			if (raw && listmode != 1) {
				d_xyindex[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint32_t), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
				d_zindex[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
				d_L[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk] * 2, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
			}
			else if (listmode != 1 && (!CT || subsets > 1)) {
				d_xyindex[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint32_t) * length[kk], NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
				d_zindex[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk], NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
				d_L[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
			}
			else {
				d_xyindex[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint32_t), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
				d_zindex[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
				d_L[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
			}
		}
	}
	if (mlem_bool) {
		d_reko_type_mlem = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint8_t) * n_rekos_mlem, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
		if (precompute)
			d_lor_mlem = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t) * koko, NULL, &status);
		else
			d_lor_mlem = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
		// Measurement data
		if (TOF && listmode != 2)
			d_Sino_mlem = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * koko * nBins, NULL, &status);
		else if (listmode != 2)
			d_Sino_mlem = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * koko, NULL, &status);
		else
			d_Sino_mlem = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
		d_norm_mlem = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_norm, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
		d_scat_mlem = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_scat, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
		if (randoms_correction == 1u)
			d_sc_ra_mlem = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * koko, NULL, &status);
		else
			d_sc_ra_mlem = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
		// Indices corresponding to the detector index (Sinogram data) or the detector number (raw data) at each measurement
		if (raw && listmode != 1) {
			d_xyindex_mlem = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint32_t), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			d_zindex_mlem = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			d_L_mlem = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t) * koko * 2, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
		}
		else if (listmode != 1 && !CT) {
			d_xyindex_mlem = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint32_t) * koko, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			d_zindex_mlem = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t) * koko, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			d_L_mlem = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
		}
		else {
			d_xyindex_mlem = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint32_t), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			d_zindex_mlem = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			d_L_mlem = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
		}
	}

	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Buffer creation failed\n");
		mexEvalString("pause(.0001);");
		return status;
	}
	else if (DEBUG) {
		mexPrintf("Buffer creation succeeded\n");
		mexEvalString("pause(.0001);");
	}


	// assign values to the buffers
	status = af_queue.enqueueWriteBuffer(d_V, CL_FALSE, 0, sizeof(float) * size_V, V);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	status = af_queue.enqueueWriteBuffer(d_x, CL_FALSE, 0, sizeof(float) * size_of_x, x);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	status = af_queue.enqueueWriteBuffer(d_y, CL_FALSE, 0, sizeof(float) * size_of_x, y);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	status = af_queue.enqueueWriteBuffer(d_z, CL_FALSE, 0, sizeof(float) * size_z, z_det);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	if (CT) {
		status = af_queue.enqueueWriteBuffer(d_angles, CL_FALSE, 0, sizeof(float) * TotSinos, angles);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
	}
	status = af_queue.enqueueWriteBuffer(d_xcenter, CL_FALSE, 0, sizeof(float) * size_center_x, x_center);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	status = af_queue.enqueueWriteBuffer(d_ycenter, CL_FALSE, 0, sizeof(float) * size_center_y, y_center);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	status = af_queue.enqueueWriteBuffer(d_zcenter, CL_FALSE, 0, sizeof(float) * size_center_z, z_center);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	status = af_queue.enqueueWriteBuffer(d_atten, CL_FALSE, 0, sizeof(float) * size_atten, atten);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	status = af_queue.enqueueWriteBuffer(d_TOFCenter, CL_FALSE, 0, sizeof(float) * nBins, TOFCenter);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	status = af_queue.enqueueWriteBuffer(d_pseudos, CL_FALSE, 0, sizeof(uint32_t) * prows, pseudos);
	status = af_queue.finish();
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return status;
	}
	if (osem_bool) {
		status = af_queue.enqueueWriteBuffer(d_reko_type, CL_FALSE, 0, sizeof(uint8_t) * n_rekos, reko_type);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
		for (uint32_t kk = osa_iter0; kk < subsetsUsed; kk++) {
			if (raw && listmode != 1) {
				status = af_queue.enqueueWriteBuffer(d_xyindex[kk], CL_FALSE, 0, sizeof(uint32_t), xy_index);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
				status = af_queue.enqueueWriteBuffer(d_zindex[kk], CL_FALSE, 0, sizeof(uint16_t), z_index);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
				status = af_queue.enqueueWriteBuffer(d_L[kk], CL_FALSE, 0, sizeof(uint16_t) * length[kk] * 2, &L[pituus[kk] * 2]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
			}
			else if (listmode != 1 && (!CT || subsets > 1)) {
				status = af_queue.enqueueWriteBuffer(d_zindex[kk], CL_FALSE, 0, sizeof(uint16_t) * length[kk], &z_index[pituus[kk]]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
				status = af_queue.enqueueWriteBuffer(d_xyindex[kk], CL_FALSE, 0, sizeof(uint32_t) * length[kk], &xy_index[pituus[kk]]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
				status = af_queue.enqueueWriteBuffer(d_L[kk], CL_FALSE, 0, sizeof(uint16_t), L);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
			}
			else {
				status = af_queue.enqueueWriteBuffer(d_xyindex[kk], CL_FALSE, 0, sizeof(uint32_t), xy_index);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
				status = af_queue.enqueueWriteBuffer(d_zindex[kk], CL_FALSE, 0, sizeof(uint16_t), z_index);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
				status = af_queue.enqueueWriteBuffer(d_L[kk], CL_FALSE, 0, sizeof(uint16_t), L);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return status;
				}
			}
			if (precompute)
				status = af_queue.enqueueWriteBuffer(d_lor[kk], CL_FALSE, 0, sizeof(uint16_t) * length[kk], &lor1[pituus[kk]]);
			else
				status = af_queue.enqueueWriteBuffer(d_lor[kk], CL_FALSE, 0, sizeof(uint16_t), lor1);
			status = af_queue.finish();
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			if (size_norm > 1ULL) {
				status = af_queue.enqueueWriteBuffer(d_norm[kk], CL_FALSE, 0, sizeof(float) * length[kk], &norm[pituus[kk]]);
			}
			else {
				status = af_queue.enqueueWriteBuffer(d_norm[kk], CL_FALSE, 0, sizeof(float) * size_norm, norm);
			}
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			if (size_scat > 1ULL) {
				status = af_queue.enqueueWriteBuffer(d_scat[kk], CL_FALSE, 0, sizeof(float) * length[kk], &scat[pituus[kk]]);
			}
			else {
				status = af_queue.enqueueWriteBuffer(d_scat[kk], CL_FALSE, 0, sizeof(float) * size_scat, scat);
			}
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			if (DEBUG) {
				mexPrintf("length[kk] = %d\n", length[kk]);
				mexPrintf("pituus[kk] = %d\n", pituus[kk]);
				mexPrintf("erotus = %d\n", pituus[kk + 1] - pituus[kk]);
				mexPrintf("kk = %d\n", kk);
				mexEvalString("pause(.0001);");
			}
			if (TOF) {
				if (loadTOF) {
					for (int64_t to = 0LL; to < nBins; to++)
						status = af_queue.enqueueWriteBuffer(d_Sino[kk], CL_FALSE, sizeof(float) * length[kk] * to, sizeof(float) * length[kk], &Sin[pituus[kk] + koko * to]);
				}
				else {
					if (kk == osa_iter0) {
						for (int64_t to = 0LL; to < nBins; to++)
							status = af_queue.enqueueWriteBuffer(d_Sino[kk], CL_FALSE, sizeof(float) * length[kk] * to, sizeof(float) * length[kk], &Sin[pituus[kk] + koko * to]);
					}
				}
			}
			else if (listmode != 2)
				status = af_queue.enqueueWriteBuffer(d_Sino[kk], CL_FALSE, 0, sizeof(float) * length[kk], &Sin[pituus[kk]]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
#if MX_HAS_INTERLEAVED_COMPLEX
			float* apu = (float*)mxGetSingles(mxGetCell(sc_ra, 0));
#else
			float* apu = (float*)mxGetData(mxGetCell(sc_ra, 0));
#endif
			if (randoms_correction)
				status = af_queue.enqueueWriteBuffer(d_sc_ra[kk], CL_FALSE, 0, sizeof(float) * length[kk], &apu[pituus[kk]]);
			else
				status = af_queue.enqueueWriteBuffer(d_sc_ra[kk], CL_FALSE, 0, sizeof(float), apu);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			status = af_queue.finish();
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
		}
	}
	if (mlem_bool) {
		status = af_queue.finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
		status = af_queue.enqueueWriteBuffer(d_reko_type_mlem, CL_FALSE, 0, sizeof(uint8_t) * n_rekos_mlem, reko_type_mlem);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
		if (precompute)
			status = af_queue.enqueueWriteBuffer(d_lor_mlem, CL_FALSE, 0, sizeof(uint16_t) * koko, lor1);
		else
			status = af_queue.enqueueWriteBuffer(d_lor_mlem, CL_FALSE, 0, sizeof(uint16_t), lor1);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
		if (TOF)
			status = af_queue.enqueueWriteBuffer(d_Sino_mlem, CL_FALSE, 0, sizeof(float) * koko * nBins, Sin);
		else if (listmode != 2)
			status = af_queue.enqueueWriteBuffer(d_Sino_mlem, CL_FALSE, 0, sizeof(float) * koko, Sin);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
		status = af_queue.enqueueWriteBuffer(d_norm_mlem, CL_FALSE, 0, sizeof(float) * size_norm, norm);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
		status = af_queue.enqueueWriteBuffer(d_scat_mlem, CL_FALSE, 0, sizeof(float) * size_scat, scat);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
#if MX_HAS_INTERLEAVED_COMPLEX
		float* apu = (float*)mxGetSingles(mxGetCell(sc_ra, 0));
#else
		float* apu = (float*)mxGetData(mxGetCell(sc_ra, 0));
#endif
		if (randoms_correction)
			status = af_queue.enqueueWriteBuffer(d_sc_ra_mlem, CL_FALSE, 0, sizeof(float) * koko, apu);
		else
			status = af_queue.enqueueWriteBuffer(d_sc_ra_mlem, CL_FALSE, 0, sizeof(float), apu);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return status;
		}
		if (raw && listmode != 1) {
			status = af_queue.enqueueWriteBuffer(d_xyindex_mlem, CL_FALSE, 0, sizeof(uint32_t), xy_index);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			status = af_queue.enqueueWriteBuffer(d_zindex_mlem, CL_FALSE, 0, sizeof(uint16_t), z_index);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			status = af_queue.enqueueWriteBuffer(d_L_mlem, CL_FALSE, 0, sizeof(uint16_t) * koko * 2ULL, L);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
		}
		else if (listmode != 1 && !CT) {
			status = af_queue.enqueueWriteBuffer(d_xyindex_mlem, CL_FALSE, 0, sizeof(uint32_t) * koko, xy_index);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			status = af_queue.enqueueWriteBuffer(d_zindex_mlem, CL_FALSE, 0, sizeof(uint16_t) * koko, z_index);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			status = af_queue.enqueueWriteBuffer(d_L_mlem, CL_FALSE, 0, sizeof(uint16_t), L);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
		}
		else {
			status = af_queue.enqueueWriteBuffer(d_xyindex_mlem, CL_FALSE, 0, sizeof(uint32_t), xy_index);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			status = af_queue.enqueueWriteBuffer(d_zindex_mlem, CL_FALSE, 0, sizeof(uint16_t), z_index);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
			status = af_queue.enqueueWriteBuffer(d_L_mlem, CL_FALSE, 0, sizeof(uint16_t), L);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return status;
			}
		}
	}

	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Buffer write failed\n");
		return status;
	}
	else if (DEBUG) {
		mexPrintf("Buffer write succeeded\n");
		mexEvalString("pause(.0001);");
	}
	return status;
}

// Prepass phase for MRAMLA, COSEM, ACOSEM, ECOSEM
void MRAMLA_prepass(const uint32_t subsets, const uint32_t im_dim, const int64_t* pituus, const std::vector<cl::Buffer> &lor, const std::vector<cl::Buffer> &zindex,
	const std::vector<cl::Buffer> &xindex, cl::Program program, const cl::CommandQueue &af_queue, const cl::Context af_context, Weighting& w_vec,
	std::vector<af::array>& Summ, std::vector<cl::Buffer> &d_Sino, const size_t koko_l, af::array& cosem, af::array& C_co,
	af::array& C_aco, af::array& C_osl, const uint32_t alku, cl::Kernel &kernel_mramla, const std::vector<cl::Buffer> &L, const uint8_t raw,
	const RecMethodsOpenCL MethodListOpenCL, const std::vector<size_t> length, const bool atomic_64bit, const bool atomic_32bit, const cl_uchar compute_norm_matrix, 
	const std::vector<cl::Buffer>& d_sc_ra, cl_uint kernelInd_MRAMLA, af::array& E, const std::vector<cl::Buffer>& d_norm, const std::vector<cl::Buffer>& d_scat, const bool use_psf,
	const af::array& g, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float epps, const bool TOF, const bool loadTOF, const mxArray* Sin, const int64_t nBins, 
	const size_t koko, const bool randoms_correction, const size_t local_size, const uint64_t* randSize, const uint32_t Nt, const bool CT) {

	cl_int status = CL_SUCCESS;

	cl_uchar MBSREM_prepass = static_cast<cl_uchar>(w_vec.MBSREM_prepass);

	af::array apu_co, apu_aco, uu, sub_index_array, apu_summa, apu_Amin, apu_summa_m;

	const bool U_skip = w_vec.U > 0.f ? false : true;

	af::array cosem_psf;
	if (use_psf) {
		cosem_psf = computeConvolution(cosem, g, Nx, Ny, Nz, w_vec, 1u);
		af::sync();
	}
	else
		cosem_psf = cosem;

	cl::Buffer d_cosem = cl::Buffer(*cosem_psf.device<cl_mem>(), true);

	af::array apu;

	uint32_t eka = alku;

	if (alku > 0u) {
		eka = alku - 1u;
		apu = af::constant(0.f, length[eka] * nBins);
	}
	else {
		apu = af::constant(0.f, 1);
	}
	cl::Buffer d_ACOSEM_lhs = cl::Buffer(*apu.device<cl_mem>(), true);

	cl_ulong st = 0ULL;

	for (uint32_t osa_iter = eka; osa_iter < subsets; osa_iter++) {

		if (DEBUG) {
			mexPrintf("prepass kernel iteration started\n");
			mexEvalString("pause(.0001);");
		}

		cl_uint kernelInd_MRAMLA_sub = kernelInd_MRAMLA;

		size_t erotus = length[osa_iter] % local_size;

		if (erotus > 0)
			erotus = (local_size - erotus);

		const size_t global_size = length[osa_iter] + erotus;

		const uint64_t m_size = static_cast<uint64_t>(length[osa_iter]);

		if (subsets > 1) {


			if ((MethodListOpenCL.MRAMLA || MethodListOpenCL.MBSREM) && w_vec.MBSREM_prepass && alku == 0u) {
				apu_Amin = af::constant(0.f, length[osa_iter] * nBins, 1);
				apu_summa_m = af::constant(0.f, length[osa_iter] * nBins, 1);
			}
			else {
				apu_Amin = af::constant(0.f, 1, 1);
				apu_summa_m = af::constant(0.f, 1, 1);
			}

			if (w_vec.MBSREM_prepass && alku == 0u) {
				if (atomic_64bit) {
					apu_summa = af::constant(0LL, im_dim, 1, s64);
				}
				else if (atomic_32bit) {
					apu_summa = af::constant(0, im_dim, 1, s32);
				}
				else {
					apu_summa = af::constant(0.f, im_dim, 1);
				}
			}
			else {
				if (atomic_64bit) {
					apu_summa = af::constant(0LL, 1, 1, s64);
				}
				else if (atomic_32bit) {
					apu_summa = af::constant(0, 1, 1, s32);
				}
				else {
					apu_summa = af::constant(0.f, 1, 1);
				}
			}

			if ((MethodListOpenCL.COSEM || MethodListOpenCL.ECOSEM || MethodListOpenCL.OSLCOSEM == 2) && alku == 0u) {
				if (atomic_64bit)
					apu_co = af::constant(0LL, im_dim, 1, s64);
				else if (atomic_32bit)
					apu_co = af::constant(0, im_dim, 1, s32);
				else
					apu_co = af::constant(0.f, im_dim, 1);
			}
			else {
				if (atomic_64bit)
					apu_co = af::constant(0LL, 1, 1, s64);
				else if (atomic_32bit)
					apu_co = af::constant(0, 1, 1, s32);
				else
					apu_co = af::constant(0.f, 1, 1);
			}

			if ((MethodListOpenCL.ACOSEM || MethodListOpenCL.OSLCOSEM == 1) && alku == 0u) {
				if (atomic_64bit)
					apu_aco = af::constant(0LL, im_dim, 1, s64);
				else if (atomic_32bit)
					apu_aco = af::constant(0, im_dim, 1, s32);
				else
					apu_aco = af::constant(0.f, im_dim, 1);
			}
			else {
				if (atomic_64bit)
					apu_aco = af::constant(0LL, 1, 1, s64);
				else if (atomic_32bit)
					apu_aco = af::constant(0, 1, 1, s32);
				else
					apu_aco = af::constant(0.f, 1, 1);
			}
			if (TOF && !loadTOF && osa_iter > 0) {
#if MX_HAS_INTERLEAVED_COMPLEX
				float* apuS = (float*)mxGetSingles(mxGetCell(Sin, 0));
#else
				float* apuS = (float*)mxGetData(mxGetCell(Sin, 0));
#endif
				d_Sino[0] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[osa_iter] * nBins, NULL, &status);
				for (int64_t to = 0LL; to < nBins; to++)
					status = af_queue.enqueueWriteBuffer(d_Sino[0], CL_FALSE, sizeof(float) * length[osa_iter] * to, sizeof(float) * length[osa_iter], &apuS[pituus[osa_iter] + koko * to]);
			}
		}
		else {

			if (w_vec.MBSREM_prepass && alku == 0u) {
				if (atomic_64bit) {
					apu_summa = af::constant(0LL, im_dim, 1, s64);
				}
				else if (atomic_32bit) {
					apu_summa = af::constant(0, im_dim, 1, s32);
				}
				else {
					apu_summa = Summ[0];
				}
			}
			else {
				if (atomic_64bit) {
					apu_summa = af::constant(0LL, 1, 1, s64);
				}
				else if (atomic_32bit) {
					apu_summa = af::constant(0, 1, 1, s32);
				}
				else {
					apu_summa = af::constant(0.f, 1, 1);
				}
			}
			apu_summa_m = E;

			if ((MethodListOpenCL.COSEM || MethodListOpenCL.ECOSEM || MethodListOpenCL.OSLCOSEM == 2) && alku == 0u) {
				if (atomic_64bit)
					apu_co = (C_co * TH).as(s64);
				else if (atomic_32bit)
					apu_co = (C_co * TH32).as(s32);
				else
					apu_co = C_co;
			}
			else {
				if (atomic_64bit)
					apu_co = af::constant(0LL, 1, 1, s64);
				else
					apu_co = af::constant(0.f, 1, 1);
			}

			if ((MethodListOpenCL.ACOSEM || MethodListOpenCL.OSLCOSEM == 1) && alku == 0) {
				if (atomic_64bit)
					apu_aco = (C_aco * TH).as(s64);
				else if (atomic_32bit)
					apu_aco = (C_aco * TH32).as(s32);
				else
					apu_aco = C_aco;
			}
			else {
				if (atomic_64bit)
					apu_aco = af::constant(0LL, 1, 1, s64);
				else if (atomic_32bit)
					apu_aco = af::constant(0, 1, 1, s32);
				else
					apu_aco = af::constant(0.f, 1, 1);
			}

			if ((MethodListOpenCL.MRAMLA || MethodListOpenCL.MBSREM) && w_vec.MBSREM_prepass && alku == 0u)
				apu_Amin = w_vec.Amin;
			else
				apu_Amin = af::constant(0.f, 1, 1);
		}

		cl::Buffer d_apu_co = cl::Buffer(*apu_co.device<cl_mem>(), true);
		cl::Buffer d_apu_aco = cl::Buffer(*apu_aco.device<cl_mem>(), true);
		cl::Buffer d_E = cl::Buffer(*apu_summa_m.device<cl_mem>(), true);
		cl::Buffer d_Amin = cl::Buffer(*apu_Amin.device<cl_mem>(), true);
		cl::Buffer d_Summ = cl::Buffer(*apu_summa.device<cl_mem>(), true);
		af::sync();

		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_norm[osa_iter]);
		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_scat[osa_iter]);
		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_Summ);
		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, lor[osa_iter]);
		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, xindex[osa_iter]);
		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, zindex[osa_iter]);
		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, L[osa_iter]);
		if (TOF && !loadTOF)
			kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_Sino[0]);
		else
			kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_Sino[osa_iter]);
		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_sc_ra[osa_iter]);
		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_cosem);
		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, alku);
		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, MBSREM_prepass);
		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_ACOSEM_lhs);
		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_Amin);
		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_apu_co);
		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_apu_aco);
		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, d_E);
		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, m_size);
		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, MethodListOpenCL);
		kernel_mramla.setArg(kernelInd_MRAMLA_sub++, st);
		cl::NDRange local(local_size);
		cl::NDRange global(global_size);
		// Compute the kernel
		status = af_queue.enqueueNDRangeKernel(kernel_mramla, cl::NullRange, global, local);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to launch the prepass kernel\n");
			return;
		}
		else if (DEBUG) {
			mexPrintf("prepass kernel launched successfully\n");
			mexEvalString("pause(.0001);");
		}

		status = af_queue.finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Queue finish failed\n");
			return;
		}
		apu_co.unlock();
		apu_aco.unlock();
		apu_summa.unlock();
		apu_summa_m.unlock();
		apu.unlock();
		apu_Amin.unlock();
		cosem_psf.unlock();
		af::sync();
		st += length[osa_iter];

		if (alku == 0u) {
			if ((MethodListOpenCL.COSEM || MethodListOpenCL.ECOSEM)) {
				if (atomic_64bit)
					C_co(af::span, osa_iter) = apu_co.as(f32) / TH;
				else if (atomic_32bit)
					C_co(af::span, osa_iter) = apu_co.as(f32) / TH32;
				else
					C_co(af::span, osa_iter) = apu_co;
				if (use_psf)
					C_co(af::span, osa_iter) = computeConvolution(C_co(af::span, osa_iter), g, Nx, Ny, Nz, w_vec, 1u) * cosem;
				else
					C_co(af::span, osa_iter) = C_co(af::span, osa_iter) * cosem;
			}
			if (MethodListOpenCL.ACOSEM) {
				if (atomic_64bit)
					C_aco(af::span, osa_iter) = apu_aco.as(f32) / TH;
				else if (atomic_32bit)
					C_aco(af::span, osa_iter) = apu_aco.as(f32) / TH32;
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
				else if (atomic_32bit)
					C_osl(af::span, osa_iter) = apu_co.as(f32) / TH32;
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
				else if (atomic_32bit)
					C_osl(af::span, osa_iter) = apu_aco.as(f32) / TH32;
				else
					C_osl(af::span, osa_iter) = apu_aco;
				if (use_psf)
					C_osl(af::span, osa_iter) = computeConvolution(C_osl(af::span, osa_iter), g, Nx, Ny, Nz, w_vec, 1u) * af::pow(cosem, w_vec.h_ACOSEM_2);
				else
					C_osl(af::span, osa_iter) = C_osl(af::span, osa_iter) * af::pow(cosem, w_vec.h_ACOSEM_2);
			}
			//if (DEBUG) {
			//	mexPrintf("co = %f\n", af::sum<float>(C_aco(af::span, osa_iter)));
			//	mexPrintf("dim0 = %u\n", C_aco(af::span, osa_iter).dims(0));
			//	mexPrintf("dim1 = %u\n", C_aco(af::span, osa_iter).dims(1));
			//}

			if ((MethodListOpenCL.MRAMLA || MethodListOpenCL.MBSREM) && w_vec.MBSREM_prepass && Nt > 1U) {
				sub_index_array = af::range(af::dim4(length[osa_iter] * nBins), 0, s64) + pituus[osa_iter] * nBins;
				if (subsets > 1)
					w_vec.Amin(sub_index_array) = apu_Amin;
				else
					w_vec.Amin = apu_Amin;
			}

			if (w_vec.MBSREM_prepass) {
				if (DEBUG) {
					mexPrintf("D = %f\n", af::sum<float>(w_vec.D));
					mexEvalString("pause(.0001);");
				}
				if (compute_norm_matrix == 0u) {
					if (atomic_64bit)
						Summ[osa_iter] = (apu_summa).as(f32) / TH;
					else if (atomic_32bit)
						Summ[osa_iter] = (apu_summa).as(f32) / TH32;
					else
						Summ[osa_iter] = apu_summa;
					af::sync();
					w_vec.D += Summ[osa_iter];
					Summ[osa_iter](Summ[osa_iter] < epps) = epps;
					if (use_psf)
						Summ[osa_iter] = computeConvolution(Summ[osa_iter], g, Nx, Ny, Nz, w_vec, 1u);
					if (DEBUG) {
						mexPrintf("Summ[osa_iter] = %f\n", af::sum<float>(Summ[osa_iter]));
						mexEvalString("pause(.0001);");
					}
				}
				else {
					if (atomic_64bit)
						Summ[0] = (apu_summa).as(f32) / TH;
					else if (atomic_32bit)
						Summ[0] = (apu_summa).as(f32) / TH32;
					else
						Summ[0] = apu_summa;
					w_vec.D += Summ[0];
					if (DEBUG) {
						mexPrintf("Summ[0] = %f\n", af::sum<float>(Summ[0]));
						mexPrintf("atomic_64bit = %d\n", atomic_64bit);
						mexEvalString("pause(.0001);");
					}
				}

				if ((MethodListOpenCL.MRAMLA || MethodListOpenCL.MBSREM) && Nt > 1U) {
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
		if (alku == 0u && (MethodListOpenCL.MBSREM || MethodListOpenCL.MRAMLA) && w_vec.MBSREM_prepass && Nt == 1U) {
			uint32_t H = osa_iter;
			uint32_t L = 0U;
			if (TOF && !loadTOF)
				H = 0;
			if (randoms_correction)
				L = osa_iter;
			const af::array Sino = afcl::array(length[osa_iter] * nBins, d_Sino[H](), f32, true);
			clRetainMemObject(d_Sino[H]());
			const af::array rand = afcl::array(randSize[L], d_sc_ra[L](), f32, true);
			clRetainMemObject(d_sc_ra[L]());
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
			if (DEBUG) {
				mexPrintf("w_vec.epsilon_mramla = %f\n", w_vec.epsilon_mramla);
				mexEvalString("pause(.0001);");
			}
		}
		af::deviceGC();
	}
	if (TOF && !loadTOF) {
#if MX_HAS_INTERLEAVED_COMPLEX
		float* apu = (float*)mxGetSingles(mxGetCell(Sin, 0));
#else
		float* apu = (float*)mxGetData(mxGetCell(Sin, 0));
#endif
		d_Sino[0] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[0] * nBins, NULL, &status);
		for (int64_t to = 0LL; to < nBins; to++)
			status = af_queue.enqueueWriteBuffer(d_Sino[0], CL_FALSE, sizeof(float) * length[0] * to, sizeof(float) * length[0], &apu[pituus[0] + koko * to]);
	}
	if (use_psf && alku == 0 && w_vec.MBSREM_prepass) {
		w_vec.D = computeConvolution(w_vec.D, g, Nx, Ny, Nz, w_vec, 1u);
	}
	if (w_vec.MBSREM_prepass)
		w_vec.D(w_vec.D <= 0.f) = 1.f;
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
		getErrorString(status);
		return;
	}
	d_x = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * numel_x, NULL, &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}
	d_y = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * numel_x, NULL, &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}
	d_pseudos = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t) * prows, NULL, &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}
	d_lor = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(uint16_t) * osa_length, NULL, &status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}
	if (raw) {
		d_L = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * osa_length * 2, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}
	else {
		d_L = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}

	status = clEnqueueWriteBuffer(commandQueues, d_x, CL_FALSE, 0, sizeof(float) * numel_x, x, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}
	status = clEnqueueWriteBuffer(commandQueues, d_y, CL_FALSE, 0, sizeof(float) * numel_x, y, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}
	status = clEnqueueWriteBuffer(commandQueues, d_z, CL_FALSE, 0, sizeof(float) * size_z, z_det, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}
	status = clEnqueueWriteBuffer(commandQueues, d_pseudos, CL_FALSE, 0, sizeof(uint32_t) * prows, pseudos, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}

	if (raw) {
		status = clEnqueueWriteBuffer(commandQueues, d_L, CL_FALSE, 0, sizeof(uint16_t) * osa_length * 2, L, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}
	else {
		status = clEnqueueWriteBuffer(commandQueues, d_L, CL_FALSE, 0, sizeof(uint16_t), L, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			getErrorString(status);
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
		getErrorString(status);
		return;
	}

	status = clEnqueueReadBuffer(commandQueues, d_lor, CL_TRUE, 0, sizeof(uint16_t) * osa_length, lor1, 1, &event1, NULL);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}

	af::sync();
	return;
}

// This function creates the OpenCL problems for the MBSREM/COSEM/etc. prepass, for OSEM program and for MLEM program
cl_int createProgram(const bool verbose, const char* k_path, cl::Context& af_context, cl::Device& af_device_id, const char* fileName, 
	cl::Program& program_os, cl::Program& program_ml, cl::Program& program_mbsrem, bool& atomic_64bit, const bool atomic_32bit, const uint32_t device, const char* header_directory,
	const uint32_t projector_type, const float crystal_size_z, const bool precompute, const uint8_t raw, const uint32_t attenuation_correction, 
	const uint32_t normalization_correction, const int32_t dec, const size_t local_size, const uint16_t n_rays, const uint16_t n_rays3D, 
	const bool find_lors, const RecMethods MethodList, const bool osem_bool, const bool mlem_bool, const uint32_t n_rekos, const uint32_t n_rekos_mlem, 
	const Weighting& w_vec, const uint32_t osa_iter0, const float cr_pz, const float dx, const bool use_psf, const uint32_t scatter, const uint32_t randoms_correction, 
	const bool TOF, const int64_t nBins, const uint8_t listmode, const bool CT) {

	cl_int status = CL_SUCCESS;

	//std::string options = header_directory;
	//options += " -cl-single-precision-constant";

	std::string kernelFile = header_directory;
	std::string kernel_path;

	kernel_path = k_path;
	kernel_path += ".cl";
	// Load the source text file
	std::ifstream sourceFile(kernel_path.c_str());
	std::string contentF((std::istreambuf_iterator<char>(sourceFile)), std::istreambuf_iterator<char>());
	// Load the header text file
	std::ifstream sourceHeader(kernelFile + "general_opencl_functions.h");
	std::string contentHeader((std::istreambuf_iterator<char>(sourceHeader)), std::istreambuf_iterator<char>());
	std::string content;
	//content = content + contentHeader;
	std::string options = "-cl-single-precision-constant";
	//options += " -cl-fast-relaxed-math";
	// Load correct header files
	if (crystal_size_z == 0.f && projector_type == 2u) {
		options += " -DCRYST";
		std::ifstream sourceHeader1(kernelFile + "general_orth_opencl_functions.h");
		std::string contentHeader1((std::istreambuf_iterator<char>(sourceHeader1)), std::istreambuf_iterator<char>());
		std::ifstream sourceHeader2(kernelFile + "opencl_functions_orth25D.h");
		std::string contentHeader2((std::istreambuf_iterator<char>(sourceHeader2)), std::istreambuf_iterator<char>());
		//content = content + contentHeader2;
		std::ifstream sourceHeader3(kernelFile + "opencl_functions_orth3D.h");
		std::string contentHeader3((std::istreambuf_iterator<char>(sourceHeader3)), std::istreambuf_iterator<char>());
		content = contentHeader + contentHeader1 + contentHeader2 + contentHeader3 + contentF;
	}
	else if ((crystal_size_z > 0.f && projector_type == 2u) || projector_type == 3u) {
		options += " -DCRYSTZ";
		std::ifstream sourceHeader1(kernelFile + "general_orth_opencl_functions.h");
		std::string contentHeader1((std::istreambuf_iterator<char>(sourceHeader1)), std::istreambuf_iterator<char>());
		std::ifstream sourceHeader3(kernelFile + "opencl_functions_orth3D.h");
		std::string contentHeader3((std::istreambuf_iterator<char>(sourceHeader3)), std::istreambuf_iterator<char>());
		//content = content + contentHeader3;
		content = contentHeader + contentHeader1 + contentHeader3 + contentF;
	}
	else
		content = contentHeader + contentF;
	// Set all preprocessor definitions
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
	if (TOF && projector_type == 1u) {
		options += " -DTOF";
	}
	if (CT)
		options += " -DCT";
	options += (" -DNBINS=" + std::to_string(nBins));
	if (listmode == 1)
		options += " -DLISTMODE";
	else if (listmode == 2)
		options += " -DLISTMODE2";
	options += " -DFP";
	if (projector_type == 1u && !precompute && (n_rays * n_rays3D) > 1) {
		options += (" -DN_RAYS=" + std::to_string(n_rays * n_rays3D));
		options += (" -DN_RAYS2D=" + std::to_string(n_rays));
		options += (" -DN_RAYS3D=" + std::to_string(n_rays3D));
	}
	if (find_lors)
		options += " -DFIND_LORS";
	if (MethodList.MRP) {
		options += " -DMEDIAN";
		options += (" -DSEARCH_WINDOW_X=" + std::to_string(w_vec.Ndx));
		options += (" -DSEARCH_WINDOW_Y=" + std::to_string(w_vec.Ndy));
		options += (" -DSEARCH_WINDOW_Z=" + std::to_string(w_vec.Ndz));
	}
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
	// Build subset-based program
	if (osem_bool) {
		std::string os_options = options;
		if ((projector_type == 2 || projector_type == 3u || TOF) && dec > 0)
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

		//status = buildProgram(verbose, k_path, af_context, af_device_id, program_os, atomic_64bit, os_options);
		status = buildProgram(verbose, content, af_context, af_device_id, program_os, atomic_64bit, atomic_32bit, os_options);
	}
	// Build MLEM (non subset) program
	if (mlem_bool) {
		std::string ml_options = options;
		if ((projector_type == 2 || projector_type == 3u || TOF) && dec > 0)
			ml_options += (" -DDEC=" + std::to_string(dec));
		ml_options += (" -DN_REKOS=" + std::to_string(n_rekos_mlem));
		ml_options += " -DAF";
		if (n_rekos_mlem == 1)
			ml_options += " -DNREKOS1";
		else if (n_rekos_mlem == 2)
			ml_options += " -DNREKOS2";

		//status = buildProgram(verbose, k_path, af_context, af_device_id, program_ml, atomic_64bit, ml_options);
		status = buildProgram(verbose, content, af_context, af_device_id, program_ml, atomic_64bit, atomic_32bit, ml_options);
	}
	// Build the prepass phase program
	if ((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI || MethodList.PKMA) && w_vec.MBSREM_prepass ||
		MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0) {
		options += " -DAF";
		options += " -DMBSREM";
		if (n_rekos == 1)
			options += " -DNREKOS1";
		else if (n_rekos == 2)
			options += " -DNREKOS2";
		if (MethodList.MRAMLA || MethodList.MBSREM)
			options += " -DMRAMLA";
		//status = buildProgram(verbose, k_path, af_context, af_device_id, program_mbsrem, atomic_64bit, options);
		status = buildProgram(verbose, content, af_context, af_device_id, program_mbsrem, atomic_64bit, atomic_32bit, options);
	}
	return status;
}

cl_int buildProgram(const bool verbose, std::string content, cl::Context& af_context, cl::Device& af_device_id, cl::Program& program,
	bool& atomic_64bit, const bool atomic_32bit, std::string options) {
	cl_int status = CL_SUCCESS;
	size_t pituus;
	if (atomic_64bit) {
		pituus = options.length();
		options += " -DCAST=long";
		options += " -DATOMIC";
		options += (" -DTH=" + std::to_string(TH));
	}
	else if (atomic_32bit) {
		options += " -DCAST=int";
		options += " -DATOMIC32";
		options += (" -DTH=" + std::to_string(TH32));
	}
	else
		options += " -DCAST=float";
	if (DEBUG)
		mexPrintf("%s\n", options.c_str());
	if (atomic_64bit) {
		cl::string apu = af_device_id.getInfo<CL_DEVICE_EXTENSIONS>();
		cl::string apu2 = "cl_khr_int64_base_atomics";
		int32_t var = apu.find(apu2);
		if (var < 0) {
			options.erase(pituus, options.size() + 1);
			options += " -DCAST=float";
			status = -1;
		}
		else {
			//std::string kernel_path_atom;

			//kernel_path_atom = k_path;
			//kernel_path_atom += ".cl";
			//// Load the source text file
			//std::ifstream sourceFile_atom(kernel_path_atom.c_str());
			//std::string content_atom((std::istreambuf_iterator<char>(sourceFile_atom)), std::istreambuf_iterator<char>());
			std::vector<std::string> testi;
			testi.push_back(content);
			cl::Program::Sources source(testi);
			program = cl::Program(af_context, source);
			//try {
			status = program.build(options.c_str());
			if (status == CL_SUCCESS) {
				mexPrintf("OpenCL program (64-bit atomics) built\n");
			}
			//}
			//catch (cl::Error& e) {
				//mexPrintf("%s\n", e.what());
			else {
				mexPrintf("Failed to build 64-bit atomics program.\n");
				if (DEBUG) {
					std::vector<cl::Device> dev;
					af_context.getInfo(CL_CONTEXT_DEVICES, &dev);
					for (int ll = 0; ll < dev.size(); ll++) {
						cl_build_status status = program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(dev[ll]);
						if (status != CL_BUILD_ERROR)
							continue;
						std::string name = dev[ll].getInfo<CL_DEVICE_NAME>();
						std::string buildlog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(dev[ll]);
						mexPrintf("Build log for %s:\n %s", name.c_str(), buildlog.c_str());
					}
					//return status;
				}
				options.erase(pituus, options.size() + 1);
				options += " -DCAST=float";
				//status = -1;
			}
			//}
		}
	}
	else
		status = -1;
	// If not, use 32-bit atomic add (float)
	if (status != CL_SUCCESS) {
		status = CL_SUCCESS;
		atomic_64bit = false;

		//std::string kernel_path;

		//kernel_path = k_path;
		//kernel_path += ".cl";
		//std::fstream sourceFile(kernel_path.c_str());
		//std::string content((std::istreambuf_iterator<char>(sourceFile)), std::istreambuf_iterator<char>());
		std::vector<std::string> testi;
		testi.push_back(content);
		cl::Program::Sources source(testi);
		program = cl::Program(af_context, source);
		//try {
		status = program.build(options.c_str());
		if (status == CL_SUCCESS) {
			mexPrintf("OpenCL program built\n");
		}
		//}
		//catch (cl::Error& e) {
		else {
			mexPrintf("Failed to build OpenCL program.\n");
			std::vector<cl::Device> dev;
			af_context.getInfo(CL_CONTEXT_DEVICES, &dev);
			for (int ll = 0; ll < dev.size(); ll++) {
				cl_build_status status = program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(dev[ll]);
				if (status != CL_BUILD_ERROR)
					continue;
				std::string name = dev[ll].getInfo<CL_DEVICE_NAME>();
				std::string buildlog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(dev[ll]);
				mexPrintf("Build log for %s:\n %s", name.c_str(), buildlog.c_str());
			}
		}
			//mexPrintf("%s\n", e.what());
			//status = -1;
		//}
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
	cl::Buffer d_W = cl::Buffer(*W.device<cl_mem>(), true);
	cl::Buffer d_input = cl::Buffer(*input.device<cl_mem>(), true);
	cl::Buffer d_padInput = cl::Buffer(*padInput.device<cl_mem>(), true);
	cl::Buffer d_gaussianNLM = cl::Buffer(*w_vec.gaussianNLM.device<cl_mem>(), true);
	size_t global_size = N * M * K;
	cl::NDRange local(1);
	cl::NDRange global(global_size);
	af::sync();
	cl::Kernel kernelNLM(OpenCLStruct.kernelNLM);
	kernelNLM.setArg(kernelIndNLM++, d_W);
	kernelNLM.setArg(kernelIndNLM++, d_input);
	kernelNLM.setArg(kernelIndNLM++, d_padInput);
	kernelNLM.setArg(kernelIndNLM++, d_gaussianNLM);
	kernelNLM.setArg(kernelIndNLM++, ndx);
	kernelNLM.setArg(kernelIndNLM++, ndy);
	kernelNLM.setArg(kernelIndNLM++, ndz);
	kernelNLM.setArg(kernelIndNLM++, nlx);
	kernelNLM.setArg(kernelIndNLM++, nly);
	kernelNLM.setArg(kernelIndNLM++, nlz);
	kernelNLM.setArg(kernelIndNLM++, N);
	kernelNLM.setArg(kernelIndNLM++, M);
	kernelNLM.setArg(kernelIndNLM++, K);
	kernelNLM.setArg(kernelIndNLM++, w_vec.h2);
	kernelNLM.setArg(kernelIndNLM++, epps);
	kernelNLM.setArg(kernelIndNLM++, Nxy);
	kernelNLM.setArg(kernelIndNLM++, min_x);
	kernelNLM.setArg(kernelIndNLM++, max_x);
	kernelNLM.setArg(kernelIndNLM++, min_y);
	kernelNLM.setArg(kernelIndNLM++, max_y);
	kernelNLM.setArg(kernelIndNLM++, min_z);
	kernelNLM.setArg(kernelIndNLM++, max_z);
	kernelNLM.setArg(kernelIndNLM++, type);
	// Compute the kernel
	status = (*OpenCLStruct.af_queue).enqueueNDRangeKernel(OpenCLStruct.kernelNLM, cl::NullRange, global, local);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Failed to launch the NLM kernel\n");
		mexEvalString("pause(.0001);");
	}

	status = (*OpenCLStruct.af_queue).finish();
	if (status != CL_SUCCESS) {
		getErrorString(status);
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