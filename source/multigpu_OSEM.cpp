/**************************************************************************
* Implementation 3 (OSEM/MLEM) matrix-free function.
*
* Copyright(C) 2020 Ville-Veikko Wettenhovi
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
#include "functions_multigpu.hpp"

// Implementation 3
void OSEM_MLEM(const cl_uint& num_devices_context, const float kerroin, const int cpu_device, const cl::Context& context, const std::vector<cl::CommandQueue>& commandQueues,
	const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, const mxArray* Sin, const mxArray* sc_ra, const uint32_t Nx,
	const uint32_t Ny, const uint32_t Nz, const uint32_t Niter, const mxArray* options, const float dx, const float dy, const float dz, const float bx,
	const float by, const float bz, const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices, const int64_t* pituus,
	const size_t koko_l, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t size_x, const uint32_t TotSinos, const bool verbose,
	const uint32_t randoms_correction, const uint32_t attenuation_correction, const uint32_t normalization, const float* atten, const size_t size_atten,
	const float* norm, const size_t size_norm, const uint32_t subsets, const float epps, const uint32_t Nt, const uint32_t* pseudos, const uint32_t det_per_ring,
	const uint32_t prows, const uint16_t* L, const uint8_t raw, const size_t size_z, const uint32_t im_dim, const cl::Kernel& kernel, const cl::Kernel& kernel_sum,
	const cl::Kernel& kernel_mlem, const cl::Kernel& kernel_convolution, const cl::Kernel& kernel_convolution_f, const cl::Kernel& kernel_vectorMult,
	const cl::Kernel& kernel_vectorDiv, const size_t numel_x, const float tube_width, const float crystal_size_z, const float* x_center, const float* y_center,
	const float* z_center, const size_t size_center_x, const size_t size_center_y, const size_t size_center_z, const bool atomic_64bit, const bool atomic_32bit, 
	const cl_uchar compute_norm_matrix, const bool precompute, const int32_t dec, const uint32_t projector_type, const uint16_t n_rays, const uint16_t n_rays3D,
	const float cr_pz, mxArray* cell, const bool osem_bool, const float global_factor, const float bmin, const float bmax, const float Vmax, const float* V,
	const size_t size_V, const size_t local_size, const bool use_psf, const float* gaussian, const size_t size_gauss, const uint32_t scatter, const bool TOF, 
	const int64_t TOFSize, const float sigma_x, const float* TOFCenter, const int64_t nBins, const cl::vector<cl::Device> devices) {

	cl_int status = CL_SUCCESS;
	cl_float zero = 0.f;
	cl_short zeroL = 0;
	cl_ulong zerou = 0ULL;
	cl_uint zeroU = 0u;
	cl_int zero32 = 0;
	const cl_ulong st = 0ULL;
	const cl_uchar fp = 0u;
	// Number of voxels in a single image
	const uint32_t Nxy = Nx * Ny;
	const bool deblur = (bool)mxGetScalar(mxGetField(options, 0, "deblurring"));
	const bool saveIter = (bool)mxGetScalar(mxGetField(options, 0, "save_iter"));
	const bool CT = (bool)mxGetScalar(mxGetField(options, 0, "CT"));

	size_t vSize = 0ULL;
	if (atomic_64bit)
		vSize = sizeof(cl_ulong);
	else if (atomic_32bit)
		vSize = sizeof(cl_uint);
	else
		vSize = sizeof(cl_float);

	bool loadTOF = true;

	size_t Ni = 0ULL;
	if (saveIter)
		Ni = static_cast<size_t>(Niter);
	const mwSize dimmi[2] = { static_cast<mwSize>(Nx) * static_cast<mwSize>(Ny) * static_cast<mwSize>(Nz), static_cast<mwSize>(Ni + 1ULL) };

	// Output matrix
	mxArray* mlem = mxCreateNumericArray(2, dimmi, mxSINGLE_CLASS, mxREAL);
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	float* ele_ml = (float*)mxGetSingles(mlem);
#else
	float* ele_ml = (float*)mxGetData(mlem);
#endif

	int64_t nProjections = 0LL;
	float dPitch = 0.f;
	float* angles = nullptr;
	uint32_t size_y = 0U;
	if (CT) {
		dPitch = (float)mxGetScalar(mxGetField(options, 0, "dPitch"));
		nProjections = (int64_t)mxGetScalar(mxGetField(options, 0, "nProjections"));
		size_y = (uint32_t)mxGetScalar(mxGetField(options, 0, "xSize"));
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		angles = (float*)mxGetSingles(mxGetField(options, 0, "angles"));
#else
		angles = (float*)mxGetData(mxGetField(options, 0, "angles"));
#endif
	}

	// Distance between rays in multi-ray Siddon
	const float dc_z = cr_pz / static_cast<float>(n_rays3D + 1);

	// How many LORs at each subset with each device
	std::vector<size_t> meas_per_gpu(subsets);
	std::vector<size_t> meas_per_cpu(subsets);
	// How many measurements should be in a single GPU and/or CPU part
	// Each subset is divided among the devices
	for (uint32_t kk = 0u; kk < subsets; kk++) {
		size_t osa_length = (pituus[kk + 1u] - pituus[kk]);
		// CPU is used if cpu_device >= 0
		if (cpu_device >= 0) {
			meas_per_cpu[kk] = static_cast<size_t>(static_cast<float>(osa_length) / (kerroin + 1.f));
			meas_per_gpu[kk] = (osa_length - meas_per_cpu[kk]) / static_cast<size_t>(num_devices_context - 1u);
		}
		else
			meas_per_gpu[kk] = osa_length / static_cast<size_t>(num_devices_context);
	}

	std::vector<cl_ulong> cumsum((num_devices_context + 1u) * subsets, 0);
	std::vector<size_t> length(num_devices_context * subsets, 0);
	// Number of measurements/LORs in the current subset in the current device
	// Heterogeneous cases are handled differently
	for (cl_uint kk = 0u; kk < subsets; kk++) {
		size_t osa_length = pituus[kk + 1] - pituus[kk];
		for (cl_uint i = 0u; i < num_devices_context; i++) {
			// Last device
			if (i == num_devices_context - 1u) {
				// If CPU is last, the GPU needs some additional LORs, based on the given options.cpu_to_gpu_factor
				if (i == cpu_device) {
					length[kk * num_devices_context + i] = meas_per_cpu[kk];
					length[kk * num_devices_context + i - 1u] = meas_per_gpu[kk] + (osa_length - meas_per_gpu[kk] * static_cast<size_t>(num_devices_context - 1u)
						- meas_per_cpu[kk]);
					cumsum[kk * num_devices_context + i] = cumsum[kk * num_devices_context + i - 1] + length[kk * num_devices_context + i - 1u];
				}
				else {
					// No CPU at all
					if (cpu_device < 0)
						length[kk * num_devices_context + i] = meas_per_gpu[kk];
					else
						length[kk * num_devices_context + i] = meas_per_gpu[kk] + (osa_length - meas_per_gpu[kk] * static_cast<size_t>(num_devices_context - 1u)
							- meas_per_cpu[kk]);
				}
			}
			// All the other devices
			else {
				// Device is a CPU
				if (i == cpu_device) {
					length[kk * num_devices_context + i] = meas_per_cpu[kk];
				}
				else
					length[kk * num_devices_context + i] = meas_per_gpu[kk];
			}
			// Cumulative sum of all LOR counts used at each subset by each device
			cumsum[kk * num_devices_context + i + 1u] = cumsum[kk * num_devices_context + i] + length[kk * num_devices_context + i];
		}
	}


	const uint8_t listmode = (uint8_t)mxGetScalar(mxGetField(options, 0, "listmode"));
	const bool computeSensImag = (bool)mxGetScalar(mxGetField(options, 0, "compute_sensitivity_image"));

	size_t size_scat = 1ULL;
	if (scatter == 1U) {
		size_scat = mxGetNumberOfElements(mxGetCell(mxGetField(options, 0, "ScatterC"), 0));
	}

	uint32_t deblur_iterations = 0U;
	if (use_psf && deblur) {
		deblur_iterations = (uint32_t)mxGetScalar(mxGetField(options, 0, "deblur_iterations"));
	}

	if (num_devices_context == 1 && TOF) {
		cl_ulong mem = devices[0].getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>(&status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}

		if (static_cast<double>(mem) * 0.75 < static_cast<double>(koko * nBins * sizeof(float)) && TOF)
			loadTOF = false;
	}

	uint32_t TOFsubsets = subsets;
	if (!loadTOF)
		TOFsubsets = 1U;

	// Memory allocation
	cl::Buffer d_gauss;
	std::vector<cl::Buffer> apu_sum(num_devices_context);
	std::vector<cl::Buffer> d0_rhs, d0_Summ, d_Summ;
	std::vector<cl::Buffer> d_z(num_devices_context);
	std::vector<cl::Buffer> d_x(num_devices_context);
	std::vector<cl::Buffer> d_y(num_devices_context);
	std::vector<cl::Buffer> d_angles(num_devices_context);
	std::vector<cl::Buffer> d_xcenter(num_devices_context);
	std::vector<cl::Buffer> d_ycenter(num_devices_context);
	std::vector<cl::Buffer> d_zcenter(num_devices_context);
	std::vector<cl::Buffer> d_V(num_devices_context);
	std::vector<cl::Buffer> d_atten(num_devices_context);
	std::vector<cl::Buffer> d_TOFCenter(num_devices_context);
	std::vector<cl::Buffer> d_pseudos(num_devices_context);
	std::vector<cl::Buffer> d_rhs(num_devices_context);
	std::vector<cl::Buffer> d_mlem(num_devices_context);
	std::vector<cl::Buffer> d_mlem_blurred(num_devices_context);
	std::vector<cl::Buffer> d_Sino(TOFsubsets * num_devices_context);
	std::vector<cl::Buffer> d_sc_ra(subsets * num_devices_context);
	std::vector<cl::Buffer> d_norm(subsets * num_devices_context);
	std::vector<cl::Buffer> d_scat(subsets * num_devices_context);
	std::vector<cl::Buffer> d_lor(subsets * num_devices_context);
	std::vector<cl::Buffer> d_xyindex(subsets * num_devices_context);
	std::vector<cl::Buffer> d_zindex(subsets * num_devices_context);
	std::vector<cl::Buffer> d_L(subsets * num_devices_context);
	std::vector<cl::Buffer> d_reko_type(num_devices_context);
	if (compute_norm_matrix == 0u)
		d_Summ.resize(subsets * num_devices_context);
	else
		d_Summ.resize(num_devices_context);
	if (num_devices_context > 1u) {
		d0_rhs.resize(num_devices_context - 1u);
		d0_Summ.resize(num_devices_context - 1u);
	}

	// Create the necessary buffers
	for (cl_uint kk = 0u; kk < subsets; kk++) {
		for (cl_uint i = 0u; i < num_devices_context; i++) {
			if (kk == 0u) {
				if (i == 0 && use_psf) {
					d_gauss = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_gauss, NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				apu_sum[i] = cl::Buffer(context, CL_MEM_READ_ONLY, vSize, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_reko_type[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint8_t), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_z[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_z, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_x[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * numel_x, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_y[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * numel_x, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				if (CT) {
					d_angles[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * nProjections, NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				d_xcenter[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_center_x, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_ycenter[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_center_y, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_zcenter[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_center_z, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_V[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_V, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_atten[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_atten, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				// TOF bin centers
				d_TOFCenter[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * nBins, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_pseudos[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t) * prows, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_mlem[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * im_dim, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_mlem_blurred[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * im_dim, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_rhs[i] = cl::Buffer(context, CL_MEM_READ_WRITE, vSize * im_dim, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				if (compute_norm_matrix == 1u) {
					d_Summ[i] = cl::Buffer(context, CL_MEM_READ_WRITE, vSize * im_dim, NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				if (i < num_devices_context - 1) {
					d0_rhs[i] = cl::Buffer(context, CL_MEM_READ_WRITE, vSize * im_dim, NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
					d0_Summ[i] = cl::Buffer(context, CL_MEM_READ_WRITE, vSize * im_dim, NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				const std::vector<cl::Memory> siirrettavat = { d_z[i] , d_x[i], d_y[i], d_xcenter[i], d_ycenter[i], d_zcenter[i], d_atten[i], d_pseudos[i], d_mlem[i],
					d_rhs[i] };
				status = commandQueues[i].enqueueMigrateMemObjects(siirrettavat, 0);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			if (TOF && num_devices_context == 1 && !loadTOF && listmode != 2) {
				if (kk == 0)
					d_Sino[kk] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * nBins, NULL, &status);
			}
			else if (listmode != 2)
				d_Sino[kk * num_devices_context + i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * length[kk * num_devices_context + i] * nBins, NULL, &status);
			else
				d_Sino[kk * num_devices_context + i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (normalization == 1u) {
				d_norm[kk * num_devices_context + i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * length[kk * num_devices_context + i], NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			else {
				d_norm[kk * num_devices_context + i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			if (randoms_correction == 1u) {
				d_sc_ra[kk * num_devices_context + i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * length[kk * num_devices_context + i], NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			else {
				d_sc_ra[kk * num_devices_context + i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			if (scatter == 1u) {
				d_scat[kk * num_devices_context + i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * length[kk * num_devices_context + i], NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			else {
				d_scat[kk * num_devices_context + i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			if (precompute)
				d_lor[kk * num_devices_context + i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk * num_devices_context + i], NULL, &status);
			else
				d_lor[kk * num_devices_context + i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (compute_norm_matrix == 0u) {
				d_Summ[kk * num_devices_context + i] = cl::Buffer(context, CL_MEM_READ_WRITE, vSize * im_dim, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			if (raw && listmode != 1) {
				d_xyindex[kk * num_devices_context + i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_zindex[kk * num_devices_context + i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_L[kk * num_devices_context + i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk * num_devices_context + i] * 2, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			else if (listmode != 1 && (!CT || subsets > 1)) {
				d_xyindex[kk * num_devices_context + i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t) * length[kk * num_devices_context + i], NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_zindex[kk * num_devices_context + i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk * num_devices_context + i], NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_L[kk * num_devices_context + i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			else {
				d_xyindex[kk * num_devices_context + i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_zindex[kk * num_devices_context + i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_L[kk * num_devices_context + i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			//const std::vector<cl::Memory> siirrettavat = { d_norm[kk * num_devices_context + i], d_sc_ra[kk * num_devices_context + i],
			//	d_lor[kk * num_devices_context + i], d_xyindex[kk * num_devices_context + i], d_zindex[kk * num_devices_context + i],
			//	d_L[kk * num_devices_context + i] };
			//status = commandQueues[i].enqueueMigrateMemObjects(siirrettavat, 0);
			//if (status != CL_SUCCESS) {
			//	getErrorString(status);
			//	return;
			//}
		}
	}

	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Buffer creation failed\n");
		mexEvalString("pause(.0001);");
		return;
	}
	else if (DEBUG) {
		mexPrintf("Buffer creation succeeded\n");
		mexEvalString("pause(.0001);");
	}

	// Fill the buffers with host data
	for (cl_uint kk = 0u; kk < subsets; kk++) {
		for (cl_uint i = 0u; i < num_devices_context; i++) {
			if (kk == 0u) {
				if (i == 0 && use_psf) {
					status = commandQueues[i].enqueueWriteBuffer(d_gauss, CL_FALSE, 0, sizeof(float) * size_gauss, gaussian);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				if (atomic_64bit)
					status = commandQueues[i].enqueueFillBuffer(apu_sum[i], zerou, 0, sizeof(cl_ulong));
				else if (atomic_32bit)
					status = commandQueues[i].enqueueFillBuffer(apu_sum[i], zero32, 0, sizeof(cl_int));
				else
					status = commandQueues[i].enqueueFillBuffer(apu_sum[i], zero, 0, sizeof(cl_float));
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d_reko_type[i], CL_FALSE, 0, sizeof(uint8_t), &fp);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d_x[i], CL_FALSE, 0, sizeof(float) * numel_x, x);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d_y[i], CL_FALSE, 0, sizeof(float) * numel_x, y);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				if (CT) {
					status = commandQueues[i].enqueueWriteBuffer(d_angles[i], CL_FALSE, 0, sizeof(float) * nProjections, angles);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				status = commandQueues[i].enqueueWriteBuffer(d_xcenter[i], CL_FALSE, 0, sizeof(float) * size_center_x, x_center);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d_ycenter[i], CL_FALSE, 0, sizeof(float) * size_center_y, y_center);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d_zcenter[i], CL_FALSE, 0, sizeof(float) * size_center_z, z_center);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d_V[i], CL_FALSE, 0, sizeof(float) * size_V, V);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d_z[i], CL_FALSE, 0, sizeof(float) * size_z, z_det);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d_atten[i], CL_FALSE, 0, sizeof(float) * size_atten, atten);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d_TOFCenter[i], CL_FALSE, 0, sizeof(float) * nBins, TOFCenter);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d_pseudos[i], CL_FALSE, 0, sizeof(uint32_t) * prows, pseudos);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
#ifdef MX_HAS_INTERLEAVED_COMPLEX
				status = commandQueues[i].enqueueWriteBuffer(d_mlem[i], CL_FALSE, 0, sizeof(float) * im_dim, (float*)mxGetSingles(mxGetField(options, 0, "x0")));
#else
				status = commandQueues[i].enqueueWriteBuffer(d_mlem[i], CL_FALSE, 0, sizeof(float) * im_dim, (float*)mxGetData(mxGetField(options, 0, "x0")));
#endif
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			if (precompute)
				status = commandQueues[i].enqueueWriteBuffer(d_lor[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint16_t) * length[kk * num_devices_context + i], 
					&lor1[cumsum[kk * num_devices_context + i]]);
			else
				status = commandQueues[i].enqueueWriteBuffer(d_lor[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint16_t), lor1);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (normalization == 1u) {
				status = commandQueues[i].enqueueWriteBuffer(d_norm[kk * num_devices_context + i], CL_FALSE, 0, sizeof(cl_float) * length[kk * num_devices_context + i], 
					&norm[cumsum[kk * num_devices_context + i]]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			else {
				status = commandQueues[i].enqueueWriteBuffer(d_norm[kk * num_devices_context + i], CL_FALSE, 0, sizeof(cl_float) * size_norm, norm);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}

			if (raw && listmode != 1) {
				status = commandQueues[i].enqueueWriteBuffer(d_xyindex[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint32_t), xy_index);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d_zindex[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint16_t), z_index);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d_L[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint16_t) * length[kk * num_devices_context + i] * 2, 
					&L[cumsum[kk * num_devices_context + i] * 2]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			else if (listmode != 1 && (!CT || subsets > 1)) {
				status = commandQueues[i].enqueueWriteBuffer(d_xyindex[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint32_t) * length[kk * num_devices_context + i], 
					&xy_index[cumsum[kk * num_devices_context + i]]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d_zindex[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint16_t) * length[kk * num_devices_context + i], 
					&z_index[cumsum[kk * num_devices_context + i]]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d_L[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint16_t), L);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			else {
				status = commandQueues[i].enqueueWriteBuffer(d_xyindex[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint32_t), xy_index);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d_zindex[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint16_t), z_index);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d_L[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint16_t), L);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}

			status = commandQueues[0].flush();
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
	}

	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Buffer write failed\n");
		return;
	}
	else if (DEBUG) {
		mexPrintf("Buffer write succeeded\n");
		mexEvalString("pause(.0001);");
	}

	size_t sum_dim = static_cast<size_t>(im_dim);

	// Temporary vectors for multi-device case (output data are stored here, then transferred to the primary device)
	std::vector<std::vector<float>> testi_summ;
	std::vector<std::vector<float>> testi_rhs;
	std::vector<std::vector<int64_t>> testi_summ_u;
	std::vector<std::vector<int64_t>> testi_rhs_u;
	std::vector<std::vector<int32_t>> testi_summ_32;
	std::vector<std::vector<int32_t>> testi_rhs_32;
	if (atomic_64bit) {
		testi_summ_u.resize(num_devices_context - 1u, std::vector<int64_t>(im_dim));
		testi_rhs_u.resize(num_devices_context - 1u, std::vector<int64_t>(im_dim));
	}
	else if (atomic_32bit) {
		testi_summ_32.resize(num_devices_context - 1u, std::vector<int32_t>(im_dim));
		testi_rhs_32.resize(num_devices_context - 1u, std::vector<int32_t>(im_dim));
	}
	else {
		testi_summ.resize(num_devices_context - 1u, std::vector<float>(im_dim));
		testi_rhs.resize(num_devices_context - 1u, std::vector<float>(im_dim));
	}

	for (cl_uint i = 0u; i < num_devices_context; i++) {
		commandQueues[i].finish();
	}

	cl_uint kernelInd = 0U;

	// Set the constant kernel arguments
	cl::Kernel kernel_ = kernel;
	cl::Kernel kernel_mlem_ = kernel_mlem;
	kernel_.setArg(kernelInd++, global_factor);
	kernel_.setArg(kernelInd++, epps);
	kernel_.setArg(kernelInd++, im_dim);
	kernel_.setArg(kernelInd++, Nx);
	kernel_.setArg(kernelInd++, Ny);
	kernel_.setArg(kernelInd++, Nz);
	kernel_.setArg(kernelInd++, dz);
	kernel_.setArg(kernelInd++, dx);
	kernel_.setArg(kernelInd++, dy);
	kernel_.setArg(kernelInd++, bz);
	kernel_.setArg(kernelInd++, bx);
	kernel_.setArg(kernelInd++, by);
	kernel_.setArg(kernelInd++, bzb);
	kernel_.setArg(kernelInd++, maxxx);
	kernel_.setArg(kernelInd++, maxyy);
	kernel_.setArg(kernelInd++, zmax);
	kernel_.setArg(kernelInd++, NSlices);
	kernel_.setArg(kernelInd++, size_x);
	kernel_.setArg(kernelInd++, TotSinos);
	kernel_.setArg(kernelInd++, det_per_ring);
	kernel_.setArg(kernelInd++, prows);
	kernel_.setArg(kernelInd++, Nxy);
	kernel_.setArg(kernelInd++, fp);
	kernel_.setArg(kernelInd++, sigma_x);
	if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1))) {
		kernel_.setArg(kernelInd++, tube_width);
		kernel_.setArg(kernelInd++, crystal_size_z);
		kernel_.setArg(kernelInd++, bmin);
		kernel_.setArg(kernelInd++, bmax);
		kernel_.setArg(kernelInd++, Vmax);
	}
	else if (projector_type == 1u && !precompute) {
		kernel_.setArg(kernelInd++, dc_z);
		kernel_.setArg(kernelInd++, n_rays);
	}
	kernel_.setArg(kernelInd++, zero);


	kernel_mlem_.setArg(3, im_dim);
	kernel_mlem_.setArg(4, epps);

	cl_uint event_count = (3u);

	int32_t w_size_x, w_size_y, w_size_z;
	if (use_psf) {
		w_size_x = (int32_t)mxGetScalar(mxGetField(options, 0, "g_dim_x"));
		w_size_y = (int32_t)mxGetScalar(mxGetField(options, 0, "g_dim_y"));
		w_size_z = (int32_t)mxGetScalar(mxGetField(options, 0, "g_dim_z"));
	}

	// If the normalization constant is only computed during the first iteration, create a sufficiently large buffer and fill it with zeros
	if (compute_norm_matrix == 0u || (listmode == 1 && computeSensImag)) {
		if (listmode == 1 && computeSensImag) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			float* apu = (float*)mxGetSingles(mxGetField(options, 0, "Summ"));
#else
			float* apu = (float*)mxGetData(mxGetField(options, 0, "Summ"));
#endif
			int64_t* summ_apu;
			int32_t* summ_apu32;
			if (atomic_64bit) {
				summ_apu = new int64_t[im_dim];
				for (int64_t aa = 0; aa < im_dim; aa++) {
					summ_apu[aa] = static_cast<int64_t>(apu[aa] * TH);
				}
			}
			else if (atomic_32bit) {
				summ_apu32 = new int32_t[im_dim];
				for (int64_t aa = 0; aa < im_dim; aa++) {
					summ_apu32[aa] = static_cast<int32_t>(apu[aa] * TH32);
				}
			}
			for (uint32_t osa_iter = 0u; osa_iter < subsets; osa_iter++) {
				for (cl_uint i = 0u; i < num_devices_context; i++) {
					if (atomic_64bit)
						status = commandQueues[i].enqueueWriteBuffer(d_Summ[osa_iter * num_devices_context + i], CL_TRUE, 0, sizeof(cl_long) * im_dim, summ_apu);
					else if (atomic_32bit)
						status = commandQueues[i].enqueueWriteBuffer(d_Summ[osa_iter * num_devices_context + i], CL_TRUE, 0, sizeof(cl_int) * im_dim, summ_apu32);
					else
						status = commandQueues[i].enqueueWriteBuffer(d_Summ[osa_iter * num_devices_context + i], CL_TRUE, 0, sizeof(cl_float) * im_dim, apu);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
			}

			for (cl_uint i = 0u; i < num_devices_context; i++) {
				status = commandQueues[i].finish();
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			if (atomic_64bit)
				delete[] summ_apu; 
			else if (atomic_32bit)
				delete[] summ_apu32;
			if (DEBUG) {
				mexPrintf("Summ write succeeded\n");
				mexEvalString("pause(.0001);");
			}
		}
		else {
			for (uint32_t osa_iter = 0u; osa_iter < subsets; osa_iter++) {
				for (cl_uint i = 0u; i < num_devices_context; i++) {
					if (atomic_64bit)
						status = commandQueues[i].enqueueFillBuffer(d_Summ[osa_iter * num_devices_context + i], zerou, 0, sizeof(cl_ulong) * im_dim);
					else if (atomic_32bit)
						status = commandQueues[i].enqueueFillBuffer(d_Summ[osa_iter * num_devices_context + i], zero32, 0, sizeof(cl_uint) * im_dim);
					else
						status = commandQueues[i].enqueueFillBuffer(d_Summ[osa_iter * num_devices_context + i], zero, 0, sizeof(cl_float) * im_dim);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
			}
		}
	}

	cl_uchar no_norm = 0u;

	// Time steps
	for (uint32_t tt = 0u; tt < Nt; tt++) {

		if ((tt > 0u && compute_norm_matrix == 0u) || (listmode == 1 && computeSensImag))
			no_norm = 1u;

		// Measurement data
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		float* Sino = (float*)mxGetSingles(mxGetCell(Sin, static_cast<mwIndex>(tt)));
#else
		float* Sino = (float*)mxGetData(mxGetCell(Sin, static_cast<mwIndex>(tt)));
#endif

		for (cl_uint i = 0u; i < num_devices_context; i++) {
			status = commandQueues[i].finish();
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}

		// Transfer data from host to device
		for (cl_uint kk = 0u; kk < subsets; kk++) {
			for (cl_uint i = 0u; i < num_devices_context; i++) {
				if (TOF && num_devices_context == 1 && !loadTOF && listmode != 2) {
					if (kk == 0) {
						for (int64_t to = 0LL; to < nBins; to++) {
							status = commandQueues[i].enqueueWriteBuffer(d_Sino[0], CL_FALSE, sizeof(float) * length[0] * to,
								sizeof(float) * length[0], &Sino[cumsum[0] + koko * to]);
						}
					}
				}
				else if (listmode != 2) {
					for (int64_t to = 0LL; to < nBins; to++) {
						status = commandQueues[i].enqueueWriteBuffer(d_Sino[kk * num_devices_context + i], CL_FALSE, sizeof(float) * length[kk * num_devices_context + i] * to,
							sizeof(float) * length[kk * num_devices_context + i], &Sino[cumsum[kk * num_devices_context + i] + koko * to]);
					}
				}
				else {
					status = commandQueues[i].enqueueFillBuffer(d_Sino[kk * num_devices_context + i], zero, 0, sizeof(cl_float));
				}
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				if (DEBUG) {
					mexPrintf("Sino write succeeded\n");
					mexEvalString("pause(.0001);");
				}
				// Randoms
				if (randoms_correction == 1u) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
					float* S_R = (float*)mxGetSingles(mxGetCell(sc_ra, static_cast<mwIndex>(tt)));
#else
					float* S_R = (float*)mxGetData(mxGetCell(sc_ra, tt));
#endif
					status = commandQueues[i].enqueueWriteBuffer(d_sc_ra[kk * num_devices_context + i], CL_FALSE, 0, sizeof(float) * length[kk * num_devices_context + i],
						&S_R[cumsum[kk * num_devices_context + i]]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				else {
					status = commandQueues[i].enqueueFillBuffer(d_sc_ra[kk * num_devices_context + i], zero, 0, sizeof(cl_float));
					//status = clEnqueueWriteBuffer(commandQueues[i], d_sc_ra[kk * num_devices_context + i], CL_FALSE, 0, sizeof(float), S_R, 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				if (scatter == 1u) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
					float* scat = (float*)mxGetSingles(mxGetCell(mxGetField(options, 0, "ScatterC"), tt));
#else
					float* scat = (float*)mxGetData(mxGetCell(mxGetField(options, 0, "ScatterC"), tt));
#endif
					status = commandQueues[i].enqueueWriteBuffer(d_scat[kk * num_devices_context + i], CL_FALSE, 0, sizeof(float) * length[kk * num_devices_context + i],
						&scat[cumsum[kk * num_devices_context + i]]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				else {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
					float* scat = (float*)mxGetSingles(mxGetCell(mxGetField(options, 0, "ScatterC"), static_cast<mwIndex>(0)));
#else
					float* scat = (float*)mxGetData(mxGetCell(mxGetField(options, 0, "ScatterC"), 0));
#endif
					//status = commandQueues[i].enqueueWriteBuffer(d_scat[kk * num_devices_context + i], CL_FALSE, 0, sizeof(float), scat);
					status = commandQueues[i].enqueueFillBuffer(d_scat[kk * num_devices_context + i], zero, 0, sizeof(cl_float));
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				status = commandQueues[i].flush();
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
		}


		// Reset the estimate with the initial value in the next time step
		if (tt > 0u) {
			for (cl_uint i = 0u; i < num_devices_context; i++) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
				status = commandQueues[i].enqueueWriteBuffer(d_mlem[i], CL_FALSE, 0, sizeof(float) * im_dim, (float*)mxGetSingles(mxGetField(options, 0, "x0")));
#else
				status = commandQueues[i].enqueueWriteBuffer(d_mlem[i], CL_FALSE, 0, sizeof(float) * im_dim, (float*)mxGetData(mxGetField(options, 0, "x0")));
#endif
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
		}


		for (cl_uint i = 0u; i < num_devices_context; i++) {
			status = commandQueues[i].finish();
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		uint32_t it = 0U;

		// Iterations
		for (uint32_t iter = 0u; iter < Niter; iter++) {
			if (saveIter)
				it = iter + 1U;

			// Do not recompute the normalization constant after the first iteration
			if (iter >= 1u && compute_norm_matrix == 0u)
				no_norm = 1u;

			// Subsets
			for (uint32_t osa_iter = 0u; osa_iter < subsets; osa_iter++) {


				if (num_devices_context == 1 && osa_iter > 0u && TOF && !loadTOF) {
					d_Sino[0] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * length[osa_iter] * nBins, NULL, &status);
					for (int64_t to = 0LL; to < nBins; to++)
						status = commandQueues[0].enqueueWriteBuffer(d_Sino[0], CL_FALSE, sizeof(float) * length[osa_iter] * to, sizeof(float) * length[osa_iter], &Sino[cumsum[osa_iter] + koko * to]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}

				std::vector<cl::Event> summ_event(1);
				std::vector<std::vector<cl::Event>> events(num_devices_context, std::vector<cl::Event>(1));
				std::vector<cl::Event> events2(1);
				std::vector<cl::Event> events3(1);

				// Reset RHS and normalization constants (fill with zeros)
				for (cl_uint i = 0u; i < num_devices_context; i++) {
					if (atomic_64bit) {
						//if (compute_norm_matrix == 1u) {
						//	status = commandQueues[i].enqueueFillBuffer(d_Summ[i], zerou, 0, sizeof(cl_ulong) * im_dim);
						//}
						//if (status != CL_SUCCESS) {
						//	getErrorString(status);
						//	return;
						//}
						status = commandQueues[i].enqueueFillBuffer(d_rhs[i], zerou, 0, sizeof(cl_ulong) * im_dim);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
					}
					else if (atomic_32bit) {
						status = commandQueues[i].enqueueFillBuffer(d_rhs[i], zero32, 0, sizeof(cl_uint) * im_dim);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
					}
					else {
						//if (compute_norm_matrix == 1u) {
						//	status = commandQueues[i].enqueueFillBuffer(d_Summ[i], zero, 0, sizeof(cl_float) * im_dim);
						//}
						//if (status != CL_SUCCESS) {
						//	getErrorString(status);
						//	return;
						//}
						status = commandQueues[i].enqueueFillBuffer(d_rhs[i], zero, 0, sizeof(cl_float) * im_dim);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
					}

					status = commandQueues[i].flush();
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}

				cl::Kernel kernel_convolution_f_ = kernel_convolution_f;

				for (cl_uint i = 0u; i < num_devices_context; i++) {
					status = commandQueues[i].finish();
				}

				// Loop through the devices
				for (cl_uint i = 0u; i < num_devices_context; i++) {


					if (use_psf) {
						status = commandQueues[i].finish();
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
						cl::NDRange convn_size(Nx, Ny, Nz);
						kernel_convolution_f_.setArg(0, d_mlem[i]);
						kernel_convolution_f_.setArg(1, d_mlem_blurred[i]);
						kernel_convolution_f_.setArg(2, d_gauss);
						kernel_convolution_f_.setArg(3, w_size_x);
						kernel_convolution_f_.setArg(4, w_size_y);
						kernel_convolution_f_.setArg(5, w_size_z);
						status = commandQueues[i].enqueueNDRangeKernel(kernel_convolution_f_, cl::NullRange, convn_size);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							mexPrintf("Failed to launch the convolution kernel\n");
							return;
						}
						else if (DEBUG) {
							mexPrintf("Convolution kernel launched successfully\n");
							mexEvalString("pause(.0001);");
						}
					}
					else {
						status = commandQueues[i].enqueueCopyBuffer(d_mlem[i], d_mlem_blurred[i], 0, 0, im_dim * sizeof(float));
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
					}
					status = commandQueues[i].finish();
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}

					cl_uint kernelIndSubIter = kernelInd;

					// Make sure the global size is divisible with 64
					size_t erotus = length[osa_iter * num_devices_context + i] % local_size;

					if (erotus > 0)
						erotus = (local_size - erotus);

					const size_t global_size = length[osa_iter * num_devices_context + i] + erotus;

					// The "true" global size
					const uint64_t m_size = static_cast<uint64_t>(length[osa_iter * num_devices_context + i]);

					cl::NDRange global(global_size);
					cl::NDRange local(local_size);

					if (DEBUG) {
						mexPrintf("global_size = %u\n", global_size);
						mexPrintf("local_size = %u\n", local_size);
						mexPrintf("st = %u\n", st);
						mexEvalString("pause(.0001);");
					}

					// Set dynamic kernel arguments
					kernel_.setArg(kernelIndSubIter++, d_TOFCenter[i]);
					kernel_.setArg(kernelIndSubIter++, d_atten[i]);
					kernel_.setArg(kernelIndSubIter++, d_pseudos[i]);
					kernel_.setArg(kernelIndSubIter++, d_x[i]);
					kernel_.setArg(kernelIndSubIter++, d_y[i]);
					kernel_.setArg(kernelIndSubIter++, d_z[i]);
					if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1))) {
						kernel_.setArg(kernelIndSubIter++, d_xcenter[i]);
						kernel_.setArg(kernelIndSubIter++, d_ycenter[i]);
						kernel_.setArg(kernelIndSubIter++, d_zcenter[i]);
						kernel_.setArg(kernelIndSubIter++, d_V[i]);
					}
					kernel_.setArg(kernelIndSubIter++, d_reko_type[i]);
					if (CT) {
						kernel_.setArg(kernelIndSubIter++, subsets);
						kernel_.setArg(kernelIndSubIter++, d_angles[i]);
						kernel_.setArg(kernelIndSubIter++, size_y);
						kernel_.setArg(kernelIndSubIter++, dPitch);
						kernel_.setArg(kernelIndSubIter++, nProjections);
					}
					kernel_.setArg(kernelIndSubIter++, d_norm[osa_iter * num_devices_context + i]);
					kernel_.setArg(kernelIndSubIter++, d_scat[osa_iter * num_devices_context + i]);
					if (compute_norm_matrix == 0u && no_norm == 0)
						kernel_.setArg(kernelIndSubIter++, d_Summ[osa_iter * num_devices_context + i]);
					else if (compute_norm_matrix == 0u && no_norm == 1)
						kernel_.setArg(kernelIndSubIter++, apu_sum[i]);
					else
						kernel_.setArg(kernelIndSubIter++, d_Summ[i]);
					kernel_.setArg(kernelIndSubIter++, d_lor[osa_iter * num_devices_context + i]);
					kernel_.setArg(kernelIndSubIter++, d_xyindex[osa_iter * num_devices_context + i]);
					kernel_.setArg(kernelIndSubIter++, d_zindex[osa_iter * num_devices_context + i]);
					kernel_.setArg(kernelIndSubIter++, d_L[osa_iter * num_devices_context + i]);
					if (TOF && !loadTOF && num_devices_context == 1U)
						kernel_.setArg(kernelIndSubIter++, d_Sino[0]);
					else
						kernel_.setArg(kernelIndSubIter++, d_Sino[osa_iter * num_devices_context + i]);
					kernel_.setArg(kernelIndSubIter++, d_sc_ra[osa_iter * num_devices_context + i]);
					kernel_.setArg(kernelIndSubIter++, d_mlem_blurred[i]);
					kernel_.setArg(kernelIndSubIter++, d_rhs[i]);
					kernel_.setArg(kernelIndSubIter++, no_norm);
					kernel_.setArg(kernelIndSubIter++, m_size);
					kernel_.setArg(kernelIndSubIter++, cumsum[osa_iter * num_devices_context + i]);
					// Compute the RHS and normalization constant
					status = commandQueues[i].enqueueNDRangeKernel(kernel_, cl::NullRange, global, local, NULL, &events[i][0]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						mexPrintf("Failed to launch the OS kernel\n");
						return;
					}
					else if (DEBUG) {
						mexPrintf("OS kernel launched successfully\n");
						mexEvalString("pause(.0001);");
					}
				}

				// Transfer kernel output from secondary devices to primary device
				if (num_devices_context > 1u) {
					for (cl_uint i = 1u; i < num_devices_context; i++) {
						if (atomic_64bit) {
							if (compute_norm_matrix == 1u)
								status = commandQueues[i].enqueueReadBuffer(d_Summ[i], CL_TRUE, 0, sizeof(cl_ulong) * im_dim, testi_summ_u[i - 1u].data(), &events[i]);
							else if (compute_norm_matrix == 0u && no_norm == 0u)
								status = commandQueues[i].enqueueReadBuffer(d_Summ[osa_iter * num_devices_context + i], CL_TRUE, 0, sizeof(cl_ulong) * im_dim, 
									testi_summ_u[i - 1u].data(), &events[i]);
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
							if ((compute_norm_matrix == 0u && no_norm == 0u) || compute_norm_matrix == 1u) {
								status = commandQueues[i].enqueueWriteBuffer(d0_Summ[i - 1u], CL_FALSE, 0, sizeof(cl_ulong) * im_dim, testi_summ_u[i - 1u].data());
								if (status != CL_SUCCESS) {
									getErrorString(status);
									return;
								}
							}
							status = commandQueues[i].enqueueReadBuffer(d_rhs[i], CL_TRUE, 0, sizeof(cl_ulong) * im_dim, testi_rhs_u[i - 1u].data(), &events[i]);
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
							status = commandQueues[i].enqueueWriteBuffer(d0_rhs[i - 1u], CL_FALSE, 0, sizeof(cl_ulong) * im_dim, testi_rhs_u[i - 1u].data());
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
						}
						else if (atomic_32bit) {
							if (compute_norm_matrix == 1u)
								status = commandQueues[i].enqueueReadBuffer(d_Summ[i], CL_TRUE, 0, sizeof(cl_uint) * im_dim, testi_summ_32[i - 1u].data(), &events[i]);
							else if (compute_norm_matrix == 0u && no_norm == 0u)
								status = commandQueues[i].enqueueReadBuffer(d_Summ[osa_iter * num_devices_context + i], CL_TRUE, 0, sizeof(cl_uint) * im_dim, 
									testi_summ_32[i - 1u].data(), &events[i]);
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
							if ((compute_norm_matrix == 0u && no_norm == 0u) || compute_norm_matrix == 1u) {
								status = commandQueues[i].enqueueWriteBuffer(d0_Summ[i - 1u], CL_FALSE, 0, sizeof(cl_uint) * im_dim, testi_summ_32[i - 1u].data());
								if (status != CL_SUCCESS) {
									getErrorString(status);
									return;
								}
							}
							status = commandQueues[i].enqueueReadBuffer(d_rhs[i], CL_TRUE, 0, sizeof(cl_uint) * im_dim, testi_rhs_32[i - 1u].data(), &events[i]);
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
							status = commandQueues[i].enqueueWriteBuffer(d0_rhs[i - 1u], CL_FALSE, 0, sizeof(cl_uint) * im_dim, testi_rhs_32[i - 1u].data());
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
						}
						else {
							if (compute_norm_matrix == 1u)
								status = commandQueues[i].enqueueReadBuffer(d_Summ[i], CL_TRUE, 0, sizeof(float) * im_dim, testi_summ[i - 1u].data(), &events[i]);
							else if (compute_norm_matrix == 0u && no_norm == 0u)
								status = commandQueues[i].enqueueReadBuffer(d_Summ[osa_iter * num_devices_context + i], CL_TRUE, 0, sizeof(float) * im_dim,
									testi_summ[i - 1u].data(), &events[i]);
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
							if ((compute_norm_matrix == 0u && no_norm == 0u) || compute_norm_matrix == 1u) {
								status = commandQueues[i].enqueueWriteBuffer(d0_Summ[i - 1u], CL_FALSE, 0, sizeof(float) * im_dim, testi_summ[i - 1u].data());
								if (status != CL_SUCCESS) {
									getErrorString(status);
									return;
								}
							}
							status = commandQueues[i].enqueueReadBuffer(d_rhs[i], CL_TRUE, 0, sizeof(float) * im_dim, testi_rhs_u[i - 1u].data(), &events[i]);
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
							status = commandQueues[i].enqueueWriteBuffer(d0_rhs[i - 1u], CL_FALSE, 0, sizeof(float) * im_dim, testi_rhs[i - 1u].data());
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
						}
					}
					for (cl_uint i = 0ULL; i < num_devices_context; i++) {
						status = commandQueues[i].finish();
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
					}
				}

				events[0][0].waitForEvents(events[0]);

				// Combine the data
				if (num_devices_context > 1u) {
					cl::Kernel kernel_sum_ = kernel_sum;
					for (cl_uint i = 1u; i < num_devices_context; i++) {
						kernel_sum_.setArg(0, d0_Summ[i - 1u]);
						if (compute_norm_matrix == 1u)
							kernel_sum_.setArg(1, d_Summ[0]);
						else
							kernel_sum_.setArg(1, d_Summ[osa_iter * num_devices_context]);
						kernel_sum_.setArg(2, d0_rhs[i - 1u]);
						kernel_sum_.setArg(3, d_rhs[0]);
						kernel_sum_.setArg(4, im_dim);
						kernel_sum_.setArg(5, no_norm);
						cl::NDRange global = sum_dim;
						status = commandQueues[i].enqueueNDRangeKernel(kernel_sum_, cl::NullRange, global);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							mexPrintf("Failed to launch the merge kernel\n");
							return;
						}
						else if (DEBUG) {
							mexPrintf("Merge kernel launched successfully\n");
							mexEvalString("pause(.0001);");
						}
						status = commandQueues[i].finish();
					}
				}

				//for (cl_uint i = 1u; i < num_devices_context; i++)
				//	clWaitForEvents(1, &events[i][3]);

				if (use_psf) {
					cl::Kernel kernel_convolution_ = kernel_convolution;
					cl::Buffer d_rhs_apu;
					status = commandQueues[0].finish();
						d_rhs_apu = cl::Buffer(context, CL_MEM_READ_ONLY, vSize * im_dim, NULL, &status);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
						status = commandQueues[0].enqueueCopyBuffer(d_rhs[0], d_rhs_apu, 0, 0, im_dim * vSize);
					status = commandQueues[0].finish();
					cl::NDRange convn_size(Nx, Ny, Nz);
					kernel_convolution_.setArg(0, d_rhs_apu);
					kernel_convolution_.setArg(1, d_rhs[0]);
					kernel_convolution_.setArg(2, d_gauss);
					kernel_convolution_.setArg(3, w_size_x);
					kernel_convolution_.setArg(4, w_size_y);
					kernel_convolution_.setArg(5, w_size_z);
					status = commandQueues[0].enqueueNDRangeKernel(kernel_convolution_, cl::NullRange, convn_size);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						mexPrintf("Failed to launch the convolution kernel\n");
						return;
					}
					else if (DEBUG) {
						mexPrintf("Convolution kernel launched successfully\n");
						mexEvalString("pause(.0001);");
					}
					if (no_norm == 0) {
						cl::Buffer d_Summ_apu;
							d_Summ_apu = cl::Buffer(context, CL_MEM_READ_ONLY, vSize * im_dim, NULL, &status);
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
							if (compute_norm_matrix == 1u)
								status = commandQueues[0].enqueueCopyBuffer(d_Summ[0], d_Summ_apu, 0, 0, im_dim * vSize);
							else
								status = commandQueues[0].enqueueCopyBuffer(d_Summ[osa_iter * num_devices_context], d_Summ_apu, 0, 0, im_dim * vSize);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
						status = commandQueues[0].finish();
						kernel_convolution_.setArg(0, d_Summ_apu);
						if (compute_norm_matrix == 1u)
							kernel_convolution_.setArg(1, d_Summ[0]);
						else
							kernel_convolution_.setArg(1, d_Summ[osa_iter * num_devices_context]);
						status = commandQueues[0].enqueueNDRangeKernel(kernel_convolution_, cl::NullRange, convn_size);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							mexPrintf("Failed to launch the convolution kernel\n");
							return;
						}
						else if (DEBUG) {
							mexPrintf("Convolution kernel launched successfully\n");
							mexEvalString("pause(.0001);");
						}
					}
					status = commandQueues[0].finish();
				}

				cl::Kernel kernel_mlem_ = kernel_mlem;
				// Compute MLEM/OSEM
				if (compute_norm_matrix == 1u)
					kernel_mlem_.setArg(0, d_Summ[0]);
				else
					kernel_mlem_.setArg(0, d_Summ[osa_iter * num_devices_context]);
				kernel_mlem_.setArg(1, d_rhs[0]);
				kernel_mlem_.setArg(2, d_mlem[0]);
				status = commandQueues[0].enqueueNDRangeKernel(kernel_mlem_, cl::NullRange, sum_dim, cl::NullRange, NULL, &events2[0]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}

				// Transfer primary device normalization constant to secondary devices
				if (compute_norm_matrix == 0u && iter == 0u && num_devices_context > 1u) {
					if (atomic_64bit) {
						status = commandQueues[0].enqueueReadBuffer(d_Summ[osa_iter * num_devices_context], CL_FALSE, 0, sizeof(cl_ulong) * im_dim, testi_summ_u[0].data(), &events2, &summ_event[0]);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
					}
					else if (atomic_32bit) {
						status = commandQueues[0].enqueueReadBuffer(d_Summ[osa_iter * num_devices_context], CL_FALSE, 0, sizeof(cl_uint) * im_dim, testi_summ_32[0].data(), &events2, &summ_event[0]);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
					}
					else {
						status = commandQueues[0].enqueueReadBuffer(d_Summ[osa_iter * num_devices_context], CL_FALSE, 0, sizeof(float) * im_dim, testi_summ[0].data(), &events2, &summ_event[0]);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
					}
					for (cl_uint i = 1u; i < num_devices_context; i++) {
						if (atomic_64bit) {
							status = commandQueues[0].enqueueWriteBuffer(d_Summ[osa_iter * num_devices_context + i], CL_FALSE, 0, sizeof(cl_ulong) * im_dim, testi_summ_u[0].data(), &summ_event);
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
						}
						else if (atomic_32bit) {
							status = commandQueues[0].enqueueWriteBuffer(d_Summ[osa_iter * num_devices_context + i], CL_FALSE, 0, sizeof(cl_ulong) * im_dim, testi_summ_32[0].data(), &summ_event);
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
						}
						else {
							status = commandQueues[0].enqueueWriteBuffer(d_Summ[osa_iter * num_devices_context + i], CL_FALSE, 0, sizeof(float) * im_dim, testi_summ[0].data(), &summ_event);
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
						}
					}
				}

				// Transfer estimate data back to host
				//status = clEnqueueReadBuffer(commandQueues[0], d_mlem[0], CL_TRUE, 0, sizeof(float) * im_dim, &ele_ml[im_dim * (iter + 1u)], 1, &events[0][1], 
				//	&events[0][2]);
				status = commandQueues[0].enqueueReadBuffer(d_mlem[0], CL_FALSE, 0, sizeof(cl_float) * im_dim, &ele_ml[im_dim * it], &events2, &events3[0]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}

				// Transfer estimate data to secondary devices
				for (cl_uint i = 1u; i < num_devices_context; i++) {
					status = commandQueues[0].enqueueWriteBuffer(d_mlem[0], CL_FALSE, 0, sizeof(cl_float) * im_dim, &ele_ml[im_dim * it], &events3);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				for (cl_uint i = 0u; i < num_devices_context; i++) {
					status = commandQueues[i].finish();
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}

				if (verbose && subsets > 1u) {
					mexPrintf("Sub-iteration %u complete\n", osa_iter + 1u);
					mexEvalString("pause(.0001);");
				}
			}
			if (use_psf && deblur && listmode != 2) {
				cl::Kernel kernel_convolution_f_ = kernel_convolution_f;
				cl::Kernel kernel_vectorDiv_ = kernel_vectorDiv;
				cl::Kernel kernel_vectorMult_ = kernel_vectorMult;
				cl::Buffer d_mlem_apu = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * im_dim, NULL, &status);
				cl::Buffer d_mlem_apu_neljas = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * im_dim, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				//cl::Buffer d_mlem_apu_toka = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * im_dim, NULL, &status);
				cl::Buffer d_mlem_apu_kolmas = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * im_dim, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				//status = clEnqueueCopyBuffer(commandQueues[0], d_mlem[0], d_mlem_apu_toka, 0, 0, im_dim * sizeof(float), NULL, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[0].enqueueCopyBuffer(d_mlem[0], d_mlem_apu_neljas, 0, 0, sizeof(cl_float) * im_dim);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				cl::NDRange convn_size(Nx, Ny, Nz);
				cl::NDRange vector_size(im_dim);
				kernel_convolution_f_.setArg(2, d_gauss);
				kernel_convolution_f_.setArg(3, w_size_x);
				kernel_convolution_f_.setArg(4, w_size_y);
				kernel_convolution_f_.setArg(5, w_size_z);
				for (int ss = 0; ss < deblur_iterations; ss++) {
					status = commandQueues[0].enqueueCopyBuffer(d_mlem[0], d_mlem_apu_kolmas, 0, 0, im_dim * sizeof(float));
					status = commandQueues[0].finish();
					kernel_convolution_f_.setArg(0, d_mlem_apu_neljas);
					kernel_convolution_f_.setArg(1, d_mlem_apu);
					status = commandQueues[0].enqueueNDRangeKernel(kernel_convolution_f_, cl::NullRange, convn_size);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						mexPrintf("Failed to launch the convolution kernel\n");
						return;
					}
					else if (DEBUG) {
						mexPrintf("Convolution kernel launched successfully\n");
						mexEvalString("pause(.0001);");
					}
					status = commandQueues[0].finish();
					kernel_vectorDiv_.setArg(0, d_mlem_apu);
					kernel_vectorDiv_.setArg(1, d_mlem_apu_kolmas);
					kernel_vectorDiv_.setArg(2, epps);
					status = commandQueues[0].enqueueNDRangeKernel(kernel_vectorDiv_, cl::NullRange, vector_size);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						mexPrintf("Failed to launch the division kernel\n");
						return;
					}
					else if (DEBUG) {
						mexPrintf("Division kernel launched successfully\n");
						mexEvalString("pause(.0001);");
					}
					status = commandQueues[0].finish();
					kernel_convolution_f_.setArg(0, d_mlem_apu_kolmas);
					kernel_convolution_f_.setArg(1, d_mlem_apu);
					status = commandQueues[0].enqueueNDRangeKernel(kernel_convolution_f_, cl::NullRange, convn_size);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						mexPrintf("Failed to launch the convolution kernel\n");
						return;
					}
					else if (DEBUG) {
						mexPrintf("Convolution kernel launched successfully\n");
						mexEvalString("pause(.0001);");
					}
					status = commandQueues[0].finish();
					kernel_vectorMult_.setArg(0, d_mlem_apu);
					kernel_vectorMult_.setArg(1, d_mlem_apu_neljas);
					status = commandQueues[0].enqueueNDRangeKernel(kernel_vectorMult_, cl::NullRange, vector_size);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						mexPrintf("Failed to launch the multiplication kernel\n");
						return;
					}
					else if (DEBUG) {
						mexPrintf("Multiplication kernel launched successfully\n");
						mexEvalString("pause(.0001);");
					}
					status = commandQueues[0].finish();
					//status = clEnqueueCopyBuffer(commandQueues[0], d_mlem_apu_toka, d_mlem_blurred[0], 0, 0, im_dim * sizeof(float), NULL, NULL, NULL);
				}
				status = commandQueues[0].finish();
				status = commandQueues[0].enqueueReadBuffer(d_mlem_apu_neljas, CL_TRUE, 0, sizeof(cl_float) * im_dim, &ele_ml[im_dim * it]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			for (cl_uint i = 0u; i < num_devices_context; i++) {
				status = commandQueues[i].finish();
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			if (verbose) {
				mexPrintf("Iteration %u complete\n", iter + 1u);
				mexEvalString("pause(.0001);");
			}
		}
		for (cl_uint i = 0u; i < num_devices_context; i++) {
			status = commandQueues[i].finish();
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		if (verbose) {
			mexPrintf("Time step %u complete\n", tt + 1u);
			mexEvalString("pause(.0001);");
		}
		status = commandQueues[0].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		mxSetCell(cell, static_cast<mwIndex>(tt), mxDuplicateArray(mlem));
	}

	for (cl_uint i = 0u; i < num_devices_context; i++) {
		status = commandQueues[i].flush();
		if (status != CL_SUCCESS) {
			getErrorString(status);
		}
		status = commandQueues[i].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
		}
	}
	mxDestroyArray(mlem);
	return;
}