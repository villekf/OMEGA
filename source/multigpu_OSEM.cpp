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
void OSEM_MLEM(const cl_uint& num_devices_context, const float kerroin, const int cpu_device, const cl_context& context, const cl_command_queue* commandQueues,
	const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, const mxArray* Sin, const mxArray* sc_ra, const uint32_t Nx,
	const uint32_t Ny, const uint32_t Nz, const uint32_t Niter, const mxArray* options, const float dx, const float dy, const float dz, const float bx,
	const float by, const float bz, const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices, const uint32_t* pituus,
	const size_t koko_l, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t size_x, const uint16_t TotSinos, const bool verbose,
	const uint32_t randoms_correction, const uint32_t attenuation_correction, const uint32_t normalization, const float* atten, const size_t size_atten,
	const float* norm, const size_t size_norm, const uint32_t subsets, const float epps, const uint32_t Nt, const uint32_t* pseudos, const uint32_t det_per_ring,
	const uint32_t prows, const uint16_t* L, const uint8_t raw, const size_t size_z, const uint32_t im_dim, const cl_kernel kernel, const cl_kernel& kernel_sum,
	const cl_kernel& kernel_mlem, const cl_kernel& kernel_convolution, const cl_kernel& kernel_convolution_f, const cl_kernel& kernel_vectorMult,
	const cl_kernel& kernel_vectorDiv, const size_t numel_x, const float tube_width, const float crystal_size_z, const float* x_center, const float* y_center,
	const float* z_center, const size_t size_center_x, const size_t size_center_y, const size_t size_center_z, const bool atomic_64bit,
	const cl_uchar compute_norm_matrix, const bool precompute, const int32_t dec, const uint32_t projector_type, const uint16_t n_rays, const uint16_t n_rays3D,
	const float cr_pz, mxArray* cell, const bool osem_bool, const float global_factor, const float bmin, const float bmax, const float Vmax, const float* V,
	const size_t size_V, const size_t local_size, const bool use_psf, const float* gaussian, const size_t size_gauss, const uint32_t scatter) {

	cl_int status = CL_SUCCESS;
	cl_float zero = 0.f;
	cl_short zeroL = 0;
	cl_ulong zerou = 0ULL;
	cl_uint zeroU = 0u;
	const cl_ulong st = 0ULL;
	const cl_uchar fp = 0u;
	// Number of voxels in a single image
	const uint32_t Nxy = Nx * Ny;
	const bool deblur = (bool)mxGetScalar(mxGetField(options, 0, "deblurring"));

	const mwSize dimmi[2] = { static_cast<mwSize>(Nx * Ny * Nz), static_cast<mwSize>(Niter + 1) };

	// Output matrix
	mxArray* mlem = mxCreateNumericArray(2, dimmi, mxSINGLE_CLASS, mxREAL);
	float* ele_ml = (float*)mxGetData(mlem);

	// Distance between rays in multi-ray Siddon
	const float dc_z = cr_pz / static_cast<float>(n_rays3D + 1);

	// How many LORs at each subset with each device
	std::vector<size_t> meas_per_gpu(subsets);
	std::vector<size_t> meas_per_cpu(subsets);
	// How many measurements should be in a single GPU and/or CPU part
	// Each subset is divided among the devices
	for (uint32_t kk = 0u; kk < subsets; kk++) {
		size_t osa_length = pituus[kk + 1u] - pituus[kk];
		// CPU is used if cpu_device >= 0
		if (cpu_device >= 0) {
			meas_per_cpu[kk] = static_cast<size_t>(static_cast<float>(osa_length) / (kerroin + 1.f));
			meas_per_gpu[kk] = (osa_length - meas_per_cpu[kk]) / static_cast<size_t>(num_devices_context - 1u);
		}
		else
			meas_per_gpu[kk] = osa_length / static_cast<size_t>(num_devices_context);
	}
	std::vector<size_t> cumsum((num_devices_context + 1u) * subsets, 0);
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

	size_t size_scat = 1ULL;
	if (scatter == 1U) {
		size_scat = mxGetNumberOfElements(mxGetCell(mxGetField(options, 0, "ScatterC"), 0));
	}

	uint32_t deblur_iterations = 0U;
	if (use_psf && deblur) {
		deblur_iterations = (uint32_t)mxGetScalar(mxGetField(options, 0, "deblur_iterations"));
	}

	// Memory allocation
	cl_mem d_gauss;
	std::vector<cl_mem> d0_rhs, d0_Summ, d_Summ;
	std::vector<cl_mem> d_z(num_devices_context, 0);
	std::vector<cl_mem> d_x(num_devices_context, 0);
	std::vector<cl_mem> d_y(num_devices_context, 0);
	std::vector<cl_mem> d_xcenter(num_devices_context, 0);
	std::vector<cl_mem> d_ycenter(num_devices_context, 0);
	std::vector<cl_mem> d_zcenter(num_devices_context, 0);
	std::vector<cl_mem> d_V(num_devices_context, 0);
	std::vector<cl_mem> d_atten(num_devices_context, 0);
	std::vector<cl_mem> d_pseudos(num_devices_context, 0);
	std::vector<cl_mem> d_rhs(num_devices_context, 0);
	std::vector<cl_mem> d_mlem(num_devices_context, 0);
	std::vector<cl_mem> d_mlem_blurred(num_devices_context, 0);
	std::vector<cl_mem> d_Sino(subsets * num_devices_context, 0);
	std::vector<cl_mem> d_sc_ra(subsets * num_devices_context, 0);
	std::vector<cl_mem> d_norm(subsets * num_devices_context, 0);
	std::vector<cl_mem> d_scat(subsets * num_devices_context, 0);
	std::vector<cl_mem> d_lor(subsets * num_devices_context, 0);
	std::vector<cl_mem> d_xyindex(subsets * num_devices_context, 0);
	std::vector<cl_mem> d_zindex(subsets * num_devices_context, 0);
	std::vector<cl_mem> d_L(subsets * num_devices_context, 0);
	std::vector<cl_mem> d_reko_type(num_devices_context, 0);
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
					d_gauss = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_gauss, NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				d_reko_type[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint8_t), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_z[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_z, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_x[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * numel_x, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_y[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * numel_x, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_xcenter[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_center_x, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_ycenter[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_center_y, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_zcenter[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_center_z, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_V[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_V, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_atten[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_atten, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_pseudos[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t) * prows, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_mlem[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * im_dim, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_mlem_blurred[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * im_dim, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				if (atomic_64bit) {
					d_rhs[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_ulong) * im_dim, NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
					if (compute_norm_matrix == 1u) {
						d_Summ[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_ulong) * im_dim, NULL, &status);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
					}
					if (i < num_devices_context - 1) {
						d0_rhs[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_ulong) * im_dim, NULL, &status);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
						d0_Summ[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_ulong) * im_dim, NULL, &status);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
					}
				}
				else {
					d_rhs[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * im_dim, NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
					if (compute_norm_matrix == 1u) {
						d_Summ[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * im_dim, NULL, &status);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
					}
					if (i < num_devices_context - 1) {
						d0_rhs[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * im_dim, NULL, &status);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
						d0_Summ[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * im_dim, NULL, &status);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
					}
				}
				const cl_mem siirrettavat[10] = { d_z[i] , d_x[i], d_y[i], d_xcenter[i], d_ycenter[i], d_zcenter[i], d_atten[i], d_pseudos[i], d_mlem[i],
					d_rhs[i] };
				status = clEnqueueMigrateMemObjects(commandQueues[i], 10, siirrettavat, 0, 0, nullptr, nullptr);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			d_Sino[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * length[kk * num_devices_context + i], NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (normalization == 1u) {
				d_norm[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * length[kk * num_devices_context + i], NULL,
					&status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			else {
				d_norm[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			if (randoms_correction == 1u) {
				d_sc_ra[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * length[kk * num_devices_context + i], NULL,
					&status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			else {
				d_sc_ra[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			if (scatter == 1u) {
				d_scat[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * length[kk * num_devices_context + i], NULL,
					&status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			else {
				d_scat[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			if (precompute)
				d_lor[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk * num_devices_context + i], NULL,
					&status);
			else
				d_lor[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (compute_norm_matrix == 0u) {
				if (atomic_64bit) {
					d_Summ[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_ulong) * im_dim, NULL, &status);
				}
				else {
					d_Summ[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_float) * im_dim, NULL, &status);
				}
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			if (raw) {
				d_xyindex[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_zindex[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_L[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk * num_devices_context + i] * 2, NULL,
					&status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			else {
				d_xyindex[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t) * length[kk * num_devices_context + i],
					NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_zindex[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk * num_devices_context + i],
					NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d_L[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			const cl_mem siirrettavat[7] = { d_Sino[kk * num_devices_context + i], d_norm[kk * num_devices_context + i], d_sc_ra[kk * num_devices_context + i],
				d_lor[kk * num_devices_context + i], d_xyindex[kk * num_devices_context + i], d_zindex[kk * num_devices_context + i],
				d_L[kk * num_devices_context + i] };
			status = clEnqueueMigrateMemObjects(commandQueues[i], 7, siirrettavat, 0, 0, nullptr, nullptr);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
	}

	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}

	// Fill the buffers with host data
	for (cl_uint kk = 0u; kk < subsets; kk++) {
		for (cl_uint i = 0u; i < num_devices_context; i++) {
			if (kk == 0u) {
				if (i == 0 && use_psf) {
					status = clEnqueueWriteBuffer(commandQueues[i], d_gauss, CL_FALSE, 0, sizeof(float) * size_gauss, gaussian, 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_reko_type[i], CL_FALSE, 0, sizeof(uint8_t), &fp, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_x[i], CL_FALSE, 0, sizeof(float) * numel_x, x, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_y[i], CL_FALSE, 0, sizeof(float) * numel_x, y, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_xcenter[i], CL_FALSE, 0, sizeof(float) * size_center_x, x_center, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_ycenter[i], CL_FALSE, 0, sizeof(float) * size_center_y, y_center, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_zcenter[i], CL_FALSE, 0, sizeof(float) * size_center_z, z_center, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_V[i], CL_FALSE, 0, sizeof(float) * size_V, V, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_z[i], CL_FALSE, 0, sizeof(float) * size_z, z_det, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_atten[i], CL_FALSE, 0, sizeof(float) * size_atten, atten, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_pseudos[i], CL_FALSE, 0, sizeof(uint32_t) * prows, pseudos, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_mlem[i], CL_FALSE, 0, sizeof(float) * im_dim, (float*)mxGetData(mxGetField(options, 0, "x0")),
					0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			if (precompute)
				status = clEnqueueWriteBuffer(commandQueues[i], d_lor[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint16_t) *
					length[kk * num_devices_context + i], &lor1[cumsum[kk * num_devices_context + i]], 0, NULL, NULL);
			else
				status = clEnqueueWriteBuffer(commandQueues[i], d_lor[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint16_t), lor1, 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (normalization == 1u) {
				status = clEnqueueWriteBuffer(commandQueues[i], d_norm[kk * num_devices_context + i], CL_FALSE, 0, sizeof(cl_float) *
					length[kk * num_devices_context + i], &norm[cumsum[kk * num_devices_context + i]], 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			else {
				status = clEnqueueWriteBuffer(commandQueues[i], d_norm[kk * num_devices_context + i], CL_FALSE, 0, sizeof(cl_float) * size_norm, norm, 0,
					NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}

			if (raw) {
				status = clEnqueueWriteBuffer(commandQueues[i], d_xyindex[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint32_t), xy_index, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_zindex[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint16_t), z_index, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_L[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint16_t) *
					length[kk * num_devices_context + i] * 2, &L[cumsum[kk * num_devices_context + i] * 2], 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			else {
				status = clEnqueueWriteBuffer(commandQueues[i], d_xyindex[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint32_t) *
					length[kk * num_devices_context + i], &xy_index[cumsum[kk * num_devices_context + i]], 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_zindex[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint16_t) *
					length[kk * num_devices_context + i], &z_index[cumsum[kk * num_devices_context + i]], 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_L[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint16_t), L, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}

			status = clFlush(commandQueues[i]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
	}

	size_t sum_dim = static_cast<size_t>(im_dim);

	// Temporary vectors for multi-device case (output data are stored here, then transferred to the primary device)
	std::vector<std::vector<float>> testi_summ;
	std::vector<std::vector<float>> testi_rhs;
	std::vector<std::vector<uint64_t>> testi_summ_u;
	std::vector<std::vector<uint64_t>> testi_rhs_u;
	if (atomic_64bit) {
		testi_summ_u.resize(num_devices_context - 1u, std::vector<uint64_t>(im_dim));
		testi_rhs_u.resize(num_devices_context - 1u, std::vector<uint64_t>(im_dim));
	}
	else {
		testi_summ.resize(num_devices_context - 1u, std::vector<float>(im_dim));
		testi_rhs.resize(num_devices_context - 1u, std::vector<float>(im_dim));
	}

	for (cl_uint i = 0u; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

	cl_uint kernelInd = 0U;

	// Set the constant kernel arguments
	clSetKernelArg(kernel, kernelInd++, sizeof(float), &global_factor);
	clSetKernelArg(kernel, kernelInd++, sizeof(float), &epps);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &im_dim);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &Nx);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &Ny);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &Nz);
	clSetKernelArg(kernel, kernelInd++, sizeof(float), &dz);
	clSetKernelArg(kernel, kernelInd++, sizeof(float), &dx);
	clSetKernelArg(kernel, kernelInd++, sizeof(float), &dy);
	clSetKernelArg(kernel, kernelInd++, sizeof(float), &bz);
	clSetKernelArg(kernel, kernelInd++, sizeof(float), &bx);
	clSetKernelArg(kernel, kernelInd++, sizeof(float), &by);
	clSetKernelArg(kernel, kernelInd++, sizeof(float), &bzb);
	clSetKernelArg(kernel, kernelInd++, sizeof(float), &maxxx);
	clSetKernelArg(kernel, kernelInd++, sizeof(float), &maxyy);
	clSetKernelArg(kernel, kernelInd++, sizeof(float), &zmax);
	clSetKernelArg(kernel, kernelInd++, sizeof(float), &NSlices);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &size_x);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint16_t), &TotSinos);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &det_per_ring);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &prows);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &Nxy);
	clSetKernelArg(kernel, kernelInd++, sizeof(cl_uchar), &fp);
	if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1))) {
		clSetKernelArg(kernel, kernelInd++, sizeof(float), &tube_width);
		clSetKernelArg(kernel, kernelInd++, sizeof(float), &crystal_size_z);
		clSetKernelArg(kernel, kernelInd++, sizeof(float), &bmin);
		clSetKernelArg(kernel, kernelInd++, sizeof(float), &bmax);
		clSetKernelArg(kernel, kernelInd++, sizeof(float), &Vmax);
	}
	else if (projector_type == 1u && !precompute) {
		clSetKernelArg(kernel, kernelInd++, sizeof(float), &dc_z);
		clSetKernelArg(kernel, kernelInd++, sizeof(uint16_t), &n_rays);
	}
	clSetKernelArg(kernel, kernelInd++, sizeof(float), &zero);


	clSetKernelArg(kernel_mlem, 3, sizeof(uint32_t), &im_dim);
	clSetKernelArg(kernel_mlem, 4, sizeof(float), &epps);

	cl_uint event_count = (3u);

	int32_t w_size_x, w_size_y, w_size_z;
	if (use_psf) {
		w_size_x = (int32_t)mxGetScalar(mxGetField(options, 0, "g_dim_x"));
		w_size_y = (int32_t)mxGetScalar(mxGetField(options, 0, "g_dim_y"));
		w_size_z = (int32_t)mxGetScalar(mxGetField(options, 0, "g_dim_z"));
	}

	// If the normalization constant is only computed during the first iteration, create a sufficiently large buffer and fill it with zeros
	if (compute_norm_matrix == 0u) {
		for (uint32_t osa_iter = 0u; osa_iter < subsets; osa_iter++) {
			for (cl_uint i = 0u; i < num_devices_context; i++) {
				if (atomic_64bit)
					status = clEnqueueFillBuffer(commandQueues[i], d_Summ[osa_iter * num_devices_context + i], &zerou, sizeof(cl_ulong), 0,
						sizeof(cl_ulong) * im_dim, 0, NULL, NULL);
				else
					status = clEnqueueFillBuffer(commandQueues[i], d_Summ[osa_iter * num_devices_context + i], &zero, sizeof(cl_float), 0,
						sizeof(cl_float) * im_dim, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
		}
	}

	cl_uchar no_norm = 0u;

	// Time steps
	for (uint32_t tt = 0u; tt < Nt; tt++) {

		if (tt > 0u && compute_norm_matrix == 0u)
			no_norm = 1u;

		// Measurement data
		float* Sino = (float*)mxGetData(mxGetCell(Sin, static_cast<mwIndex>(tt)));

		for (cl_uint i = 0u; i < num_devices_context; i++) {
			clFinish(commandQueues[i]);
		}

		// Transfer data from host to device
		for (cl_uint kk = 0u; kk < subsets; kk++) {
			for (cl_uint i = 0u; i < num_devices_context; i++) {
				status = clEnqueueWriteBuffer(commandQueues[i], d_Sino[kk * num_devices_context + i], CL_FALSE, 0, sizeof(float) *
					length[kk * num_devices_context + i], &Sino[cumsum[kk * num_devices_context + i]], 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				// Randoms
				if (randoms_correction == 1u) {
					float* S_R = (float*)mxGetData(mxGetCell(sc_ra, tt));
					status = clEnqueueWriteBuffer(commandQueues[i], d_sc_ra[kk * num_devices_context + i], CL_FALSE, 0, sizeof(float) *
						length[kk * num_devices_context + i], &S_R[cumsum[kk * num_devices_context + i]], 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				else {
					status = clEnqueueFillBuffer(commandQueues[i], d_sc_ra[kk * num_devices_context + i], &zero, sizeof(cl_float), 0, sizeof(cl_float), 0, NULL, NULL);
					//status = clEnqueueWriteBuffer(commandQueues[i], d_sc_ra[kk * num_devices_context + i], CL_FALSE, 0, sizeof(float), S_R, 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				if (scatter == 1u) {
					float* scat = (float*)mxGetData(mxGetCell(mxGetField(options, 0, "ScatterC"), tt));
					status = clEnqueueWriteBuffer(commandQueues[i], d_scat[kk * num_devices_context + i], CL_FALSE, 0, sizeof(float) *
						length[kk * num_devices_context + i], &scat[cumsum[kk * num_devices_context + i]], 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				else {
					status = clEnqueueFillBuffer(commandQueues[i], d_scat[kk * num_devices_context + i], &zero, sizeof(cl_float), 0, sizeof(cl_float), 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				status = clFlush(commandQueues[i]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
		}


		// Reset the estimate with the initial value in the next time step
		if (tt > 0u) {
			for (cl_uint i = 0u; i < num_devices_context; i++) {
				status = clEnqueueWriteBuffer(commandQueues[i], d_mlem[i], CL_TRUE, 0, sizeof(float) * im_dim, (float*)mxGetData(mxGetField(options, 0, "x0")),
					0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
		}


		for (cl_uint i = 0u; i < num_devices_context; i++) {
			clFinish(commandQueues[i]);
		}

		// Iterations
		for (uint32_t iter = 0u; iter < Niter; iter++) {

			// Do not recompute the normalization constant after the first iteration
			if (iter == 1u && compute_norm_matrix == 0u)
				no_norm = 1u;

			// Subsets
			for (uint32_t osa_iter = 0u; osa_iter < subsets; osa_iter++) {



				cl_event summ_event;
				cl_event** events = (cl_event**)malloc(num_devices_context * sizeof(cl_event**));
				for (cl_uint i = 0u; i < num_devices_context; i++) {
					events[i] = (cl_event*)malloc(event_count * sizeof(cl_event**));
				}

				// Reset RHS and normalization constants (fill with zeros)
				for (cl_uint i = 0u; i < num_devices_context; i++) {
					if (atomic_64bit) {
						if (compute_norm_matrix == 1u) {
							status = clEnqueueFillBuffer(commandQueues[i], d_Summ[i], &zerou, sizeof(cl_ulong), 0, sizeof(cl_ulong) * im_dim, 0, NULL, NULL);
						}
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
						status = clEnqueueFillBuffer(commandQueues[i], d_rhs[i], &zerou, sizeof(cl_ulong), 0, sizeof(cl_ulong) * im_dim, 0, NULL, NULL);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
					}
					else {
						if (compute_norm_matrix == 1u) {
							status = clEnqueueFillBuffer(commandQueues[i], d_Summ[i], &zero, sizeof(cl_float), 0, sizeof(cl_float) * im_dim, 0, NULL, NULL);
						}
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
						status = clEnqueueFillBuffer(commandQueues[i], d_rhs[i], &zero, sizeof(cl_float), 0, sizeof(cl_float) * im_dim, 0, NULL, NULL);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
					}

					status = clFlush(commandQueues[i]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}

				for (cl_uint i = 0u; i < num_devices_context; i++) {
					clFinish(commandQueues[i]);
				}

				// Loop through the devices
				for (cl_uint i = 0u; i < num_devices_context; i++) {


					if (use_psf) {
						status = clFinish(commandQueues[i]);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
						size_t convn_size[] = { Nx, Ny, Nz };
						clSetKernelArg(kernel_convolution_f, 0, sizeof(cl_mem), &d_mlem[i]);
						clSetKernelArg(kernel_convolution_f, 1, sizeof(cl_mem), &d_mlem_blurred[i]);
						clSetKernelArg(kernel_convolution_f, 2, sizeof(cl_mem), &d_gauss);
						clSetKernelArg(kernel_convolution_f, 3, sizeof(int32_t), &w_size_x);
						clSetKernelArg(kernel_convolution_f, 4, sizeof(int32_t), &w_size_y);
						clSetKernelArg(kernel_convolution_f, 5, sizeof(int32_t), &w_size_z);
						status = clEnqueueNDRangeKernel(commandQueues[i], kernel_convolution_f, 3ULL, NULL, convn_size, NULL, 0, NULL, NULL);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							mexPrintf("Failed to launch the convolution kernel\n");
							return;
						}
					}
					else {
						status = clEnqueueCopyBuffer(commandQueues[i], d_mlem[i], d_mlem_blurred[i], 0, 0, im_dim * sizeof(float), 0, NULL, NULL);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
					}
					status = clFinish(commandQueues[i]);
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

					// Set dynamic kernel arguments
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_atten[i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_pseudos[i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_x[i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_y[i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_z[i]);
					if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1))) {
						clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_xcenter[i]);
						clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_ycenter[i]);
						clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_zcenter[i]);
						clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_V[i]);
					}
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_reko_type[i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_norm[osa_iter * num_devices_context + i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_scat[osa_iter * num_devices_context + i]);
					if (compute_norm_matrix == 0u)
						clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_Summ[osa_iter * num_devices_context + i]);
					else
						clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_Summ[i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_lor[osa_iter * num_devices_context + i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_xyindex[osa_iter * num_devices_context + i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_zindex[osa_iter * num_devices_context + i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_L[osa_iter * num_devices_context + i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_Sino[osa_iter * num_devices_context + i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_sc_ra[osa_iter * num_devices_context + i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_mlem_blurred[i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_rhs[i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_uchar), &no_norm);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(uint64_t), &m_size);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_ulong), &st);
					// Compute the RHS and normalization constant
					status = clEnqueueNDRangeKernel(commandQueues[i], kernel, 1u, NULL, &global_size, &local_size, 0u, NULL, &events[i][0]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						mexPrintf("Failed to launch the OS kernel\n");
						for (cl_uint kk = 0u; kk < subsets; kk++) {
							for (cl_uint i = 0u; i < num_devices_context; i++) {
								if (kk == 0) {
									clReleaseMemObject(d_reko_type[i]);
									clReleaseMemObject(d_z[i]);
									clReleaseMemObject(d_x[i]);
									clReleaseMemObject(d_y[i]);
									clReleaseMemObject(d_xcenter[i]);
									clReleaseMemObject(d_ycenter[i]);
									clReleaseMemObject(d_zcenter[i]);
									clReleaseMemObject(d_V[i]);
									clReleaseMemObject(d_atten[i]);
									clReleaseMemObject(d_pseudos[i]);
									if (compute_norm_matrix == 1u)
										clReleaseMemObject(d_Summ[i]);
									clReleaseMemObject(d_rhs[i]);
									clReleaseMemObject(d_mlem[i]);
									clReleaseMemObject(d_mlem_blurred[i]);
									if (i < num_devices_context - 1) {
										clReleaseMemObject(d0_Summ[i]);
										clReleaseMemObject(d0_rhs[i]);
									}
								}
								if (compute_norm_matrix == 0u)
									clReleaseMemObject(d_Summ[kk * num_devices_context + i]);
								clReleaseMemObject(d_xyindex[kk * num_devices_context + i]);
								clReleaseMemObject(d_zindex[kk * num_devices_context + i]);
								clReleaseMemObject(d_L[kk * num_devices_context + i]);
								clReleaseMemObject(d_lor[kk * num_devices_context + i]);
								clReleaseMemObject(d_Sino[kk * num_devices_context + i]);
								clReleaseMemObject(d_sc_ra[kk * num_devices_context + i]);
								clReleaseMemObject(d_norm[kk * num_devices_context + i]);
								clReleaseMemObject(d_scat[kk * num_devices_context + i]);
							}
						}
						return;
					}
				}

				// Transfer kernel output from secondary devices to primary device
				if (num_devices_context > 1u) {
					for (cl_uint i = 1u; i < num_devices_context; i++) {
						if (atomic_64bit) {
							if (compute_norm_matrix == 1u)
								status = clEnqueueReadBuffer(commandQueues[i], d_Summ[i], CL_TRUE, 0, sizeof(cl_ulong) * im_dim, testi_summ_u[i - 1u].data(), 1,
									&events[i][0], NULL);
							else if (compute_norm_matrix == 0u && iter == 0u)
								status = clEnqueueReadBuffer(commandQueues[i], d_Summ[osa_iter * num_devices_context + i], CL_TRUE, 0, sizeof(cl_ulong) * im_dim,
									testi_summ_u[i - 1u].data(), 1, &events[i][0], NULL);
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
							if ((compute_norm_matrix == 0u && iter == 0u) || compute_norm_matrix == 1u) {
								status = clEnqueueWriteBuffer(commandQueues[0], d0_Summ[i - 1u], CL_FALSE, 0, sizeof(cl_ulong) * im_dim, testi_summ_u[i - 1u].data(),
									0, NULL, NULL);
								if (status != CL_SUCCESS) {
									getErrorString(status);
									return;
								}
							}
							status = clEnqueueReadBuffer(commandQueues[i], d_rhs[i], CL_TRUE, 0, sizeof(cl_ulong) * im_dim, testi_rhs_u[i - 1u].data(), 1,
								&events[i][0], NULL);
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
							status = clEnqueueWriteBuffer(commandQueues[0], d0_rhs[i - 1u], CL_FALSE, 0, sizeof(cl_ulong) * im_dim, testi_rhs_u[i - 1u].data(), 0,
								NULL, NULL);
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
						}
						else {
							if (compute_norm_matrix == 1u)
								status = clEnqueueReadBuffer(commandQueues[i], d_Summ[i], CL_TRUE, 0, sizeof(float) * im_dim, testi_summ[i - 1u].data(), 1,
									&events[i][0], NULL);
							else if (compute_norm_matrix == 0u && iter == 0u)
								status = clEnqueueReadBuffer(commandQueues[i], d_Summ[osa_iter * num_devices_context + i], CL_TRUE, 0, sizeof(float) * im_dim,
									testi_summ[i - 1u].data(), 1, &events[i][0], NULL);
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
							if ((compute_norm_matrix == 0u && iter == 0u) || compute_norm_matrix == 1u) {
								status = clEnqueueWriteBuffer(commandQueues[0], d0_Summ[i - 1u], CL_FALSE, 0, sizeof(float) * im_dim, testi_summ[i - 1u].data(),
									0, NULL, NULL);
								if (status != CL_SUCCESS) {
									getErrorString(status);
									return;
								}
							}
							status = clEnqueueReadBuffer(commandQueues[i], d_rhs[i], CL_TRUE, 0, sizeof(float) * im_dim, testi_rhs[i - 1u].data(), 1,
								&events[i][0], NULL);
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
							status = clEnqueueWriteBuffer(commandQueues[0], d0_rhs[i - 1u], CL_FALSE, 0, sizeof(float) * im_dim, testi_rhs[i - 1u].data(), 0,
								NULL, NULL);
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
						}
					}
					for (cl_uint i = 0ULL; i < num_devices_context; i++) {
						clFinish(commandQueues[i]);
					}
				}

				clWaitForEvents(1, &events[0][0]);

				// Combine the data
				if (num_devices_context > 1u) {
					for (cl_uint i = 1u; i < num_devices_context; i++) {
						clSetKernelArg(kernel_sum, 0, sizeof(cl_mem), &d0_Summ[i - 1u]);
						if (compute_norm_matrix == 1u)
							clSetKernelArg(kernel_sum, 1, sizeof(cl_mem), &d_Summ[0]);
						else
							clSetKernelArg(kernel_sum, 1, sizeof(cl_mem), &d_Summ[osa_iter * num_devices_context]);
						clSetKernelArg(kernel_sum, 2, sizeof(cl_mem), &d0_rhs[i - 1u]);
						clSetKernelArg(kernel_sum, 3, sizeof(cl_mem), &d_rhs[0]);
						clSetKernelArg(kernel_sum, 4, sizeof(uint32_t), &im_dim);
						clSetKernelArg(kernel_sum, 5, sizeof(cl_uchar), &no_norm);
						status = clEnqueueNDRangeKernel(commandQueues[0], kernel_sum, 1, NULL, &sum_dim, NULL, 0, NULL, NULL);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							mexPrintf("Failed to launch the merge kernel\n");
							return;
						}
						clFinish(commandQueues[0]);
					}
				}

				//for (cl_uint i = 1u; i < num_devices_context; i++)
				//	clWaitForEvents(1, &events[i][3]);

				if (use_psf) {
					cl_mem d_rhs_apu;
					if (atomic_64bit) {
						d_rhs_apu = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_ulong) * im_dim, NULL, &status);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
						status = clEnqueueCopyBuffer(commandQueues[0], d_rhs[0], d_rhs_apu, 0, 0, im_dim * sizeof(cl_ulong), 0, NULL, NULL);
					}
					else {
						d_rhs_apu = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * im_dim, NULL, &status);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
						status = clEnqueueCopyBuffer(commandQueues[0], d_rhs[0], d_rhs_apu, 0, 0, im_dim * sizeof(float), 0, NULL, NULL);
					}
					clFinish(commandQueues[0]);
					size_t convn_size[] = { Nx, Ny, Nz };
					clSetKernelArg(kernel_convolution, 0, sizeof(cl_mem), &d_rhs_apu);
					clSetKernelArg(kernel_convolution, 1, sizeof(cl_mem), &d_rhs[0]);
					clSetKernelArg(kernel_convolution, 2, sizeof(cl_mem), &d_gauss);
					clSetKernelArg(kernel_convolution, 3, sizeof(int32_t), &w_size_x);
					clSetKernelArg(kernel_convolution, 4, sizeof(int32_t), &w_size_y);
					clSetKernelArg(kernel_convolution, 5, sizeof(int32_t), &w_size_z);
					status = clEnqueueNDRangeKernel(commandQueues[0], kernel_convolution, 3, NULL, convn_size, NULL, 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						mexPrintf("Failed to launch the convolution kernel\n");
						return;
					}
					if (no_norm == 0) {
						cl_mem d_Summ_apu;
						if (atomic_64bit) {
							d_Summ_apu = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_ulong) * im_dim, NULL, &status);
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
							if (compute_norm_matrix == 1u)
								status = clEnqueueCopyBuffer(commandQueues[0], d_Summ[0], d_Summ_apu, 0, 0, im_dim * sizeof(cl_ulong), 0, NULL, NULL);
							else
								status = clEnqueueCopyBuffer(commandQueues[0], d_Summ[osa_iter * num_devices_context], d_Summ_apu, 0, 0, im_dim * sizeof(cl_ulong), 0, NULL, NULL);
						}
						else {
							d_Summ_apu = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * im_dim, NULL, &status);
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
							if (compute_norm_matrix == 1u)
								status = clEnqueueCopyBuffer(commandQueues[0], d_Summ[0], d_Summ_apu, 0, 0, im_dim * sizeof(float), 0, NULL, NULL);
							else
								status = clEnqueueCopyBuffer(commandQueues[0], d_Summ[osa_iter * num_devices_context], d_Summ_apu, 0, 0, im_dim * sizeof(float), 0, NULL, NULL);
						}
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
						clFinish(commandQueues[0]);
						clSetKernelArg(kernel_convolution, 0, sizeof(cl_mem), &d_Summ_apu);
						if (compute_norm_matrix == 1u)
							clSetKernelArg(kernel_convolution, 1, sizeof(cl_mem), &d_Summ[0]);
						else
							clSetKernelArg(kernel_convolution, 1, sizeof(cl_mem), &d_Summ[osa_iter * num_devices_context]);
						status = clEnqueueNDRangeKernel(commandQueues[0], kernel_convolution, 3, NULL, convn_size, NULL, 0, NULL, NULL);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							mexPrintf("Failed to launch the convolution kernel\n");
							return;
						}
						clFinish(commandQueues[0]);
						clReleaseMemObject(d_Summ_apu);
					}
					clFinish(commandQueues[0]);
					clReleaseMemObject(d_rhs_apu);
				}

				// Compute MLEM/OSEM
				if (compute_norm_matrix == 1u)
					clSetKernelArg(kernel_mlem, 0, sizeof(cl_mem), &d_Summ[0]);
				else
					clSetKernelArg(kernel_mlem, 0, sizeof(cl_mem), &d_Summ[osa_iter * num_devices_context]);
				clSetKernelArg(kernel_mlem, 1, sizeof(cl_mem), &d_rhs[0]);
				clSetKernelArg(kernel_mlem, 2, sizeof(cl_mem), &d_mlem[0]);
				status = clEnqueueNDRangeKernel(commandQueues[0], kernel_mlem, 1, NULL, &sum_dim, NULL, 0, NULL, &events[0][1]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					for (cl_uint kk = 0u; kk < subsets; kk++) {
						for (cl_uint i = 0u; i < num_devices_context; i++) {
							if (kk == 0) {
								clReleaseMemObject(d_reko_type[i]);
								clReleaseMemObject(d_z[i]);
								clReleaseMemObject(d_x[i]);
								clReleaseMemObject(d_y[i]);
								clReleaseMemObject(d_xcenter[i]);
								clReleaseMemObject(d_ycenter[i]);
								clReleaseMemObject(d_zcenter[i]);
								clReleaseMemObject(d_V[i]);
								clReleaseMemObject(d_atten[i]);
								clReleaseMemObject(d_pseudos[i]);
								if (compute_norm_matrix == 1u)
									clReleaseMemObject(d_Summ[i]);
								clReleaseMemObject(d_rhs[i]);
								clReleaseMemObject(d_mlem[i]);
								clReleaseMemObject(d_mlem_blurred[i]);
								if (i < num_devices_context - 1) {
									clReleaseMemObject(d0_Summ[i]);
									clReleaseMemObject(d0_rhs[i]);
								}
							}
							if (compute_norm_matrix == 0u)
								clReleaseMemObject(d_Summ[kk * num_devices_context + i]);
							clReleaseMemObject(d_xyindex[kk * num_devices_context + i]);
							clReleaseMemObject(d_zindex[kk * num_devices_context + i]);
							clReleaseMemObject(d_L[kk * num_devices_context + i]);
							clReleaseMemObject(d_lor[kk * num_devices_context + i]);
							clReleaseMemObject(d_Sino[kk * num_devices_context + i]);
							clReleaseMemObject(d_sc_ra[kk * num_devices_context + i]);
							clReleaseMemObject(d_norm[kk * num_devices_context + i]);
							clReleaseMemObject(d_scat[kk * num_devices_context + i]);
						}
					}
					return;
				}

				// Transfer primary device normalization constant to secondary devices
				if (compute_norm_matrix == 0u && iter == 0u && num_devices_context > 1u) {
					if (atomic_64bit) {
						status = clEnqueueReadBuffer(commandQueues[0], d_Summ[osa_iter * num_devices_context], CL_FALSE, 0, sizeof(cl_ulong) * im_dim,
							testi_summ[0].data(), 1, &events[0][1], &summ_event);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
					}
					else {
						status = clEnqueueReadBuffer(commandQueues[0], d_Summ[osa_iter * num_devices_context], CL_FALSE, 0, sizeof(float) * im_dim,
							testi_summ[0].data(), 1, &events[0][1], &summ_event);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
					}
					for (cl_uint i = 1u; i < num_devices_context; i++) {
						if (atomic_64bit) {
							status = clEnqueueWriteBuffer(commandQueues[0], d_Summ[osa_iter * num_devices_context + i], CL_FALSE, 0, sizeof(cl_ulong) * im_dim,
								testi_summ[0].data(), 1, &summ_event, NULL);
							if (status != CL_SUCCESS) {
								getErrorString(status);
								return;
							}
						}
						else {
							status = clEnqueueWriteBuffer(commandQueues[0], d_Summ[osa_iter * num_devices_context + i], CL_FALSE, 0, sizeof(float) * im_dim,
								testi_summ[0].data(), 1, &summ_event, NULL);
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
				status = clEnqueueReadBuffer(commandQueues[0], d_mlem[0], CL_TRUE, 0, sizeof(float) * im_dim, &ele_ml[im_dim * (iter + 1u)], 1, &events[0][1],
					&events[0][2]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}

				// Transfer estimate data to secondary devices
				for (cl_uint i = 1u; i < num_devices_context; i++) {
					status = clEnqueueWriteBuffer(commandQueues[i], d_mlem[i], CL_FALSE, 0, sizeof(cl_float) * im_dim, &ele_ml[im_dim * (iter + 1u)], 1,
						&events[0][2], NULL);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				for (cl_uint i = 0u; i < num_devices_context; i++) {
					clFinish(commandQueues[i]);
					if (i == 0u) {
						for (cl_uint k = 0u; k <= 2; k++)
							clReleaseEvent(events[i][k]);
					}
					else {
						//for (cl_uint k = 0u; k <= 3u; k++)
						clReleaseEvent(events[i][0]);
					}
					clFinish(commandQueues[i]);
				}

				if (compute_norm_matrix == 0u && iter == 0u && num_devices_context > 1u)
					clReleaseEvent(summ_event);
				free(events);
				if (verbose && subsets > 1u) {
					mexPrintf("Sub-iteration %u complete\n", osa_iter + 1u);
					mexEvalString("pause(.0001);");
				}
			}
			if (use_psf && deblur) {
				cl_mem d_mlem_apu = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * im_dim, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				//cl_mem d_mlem_apu_toka = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * im_dim, NULL, &status);
				cl_mem d_mlem_apu_kolmas = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * im_dim, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				//status = clEnqueueCopyBuffer(commandQueues[0], d_mlem[0], d_mlem_apu_toka, 0, 0, im_dim * sizeof(float), NULL, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueCopyBuffer(commandQueues[0], d_mlem[0], d_mlem_blurred[0], 0, 0, im_dim * sizeof(float), 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				size_t convn_size[] = { Nx, Ny, Nz };
				size_t vector_size = static_cast<size_t>(im_dim);
				clSetKernelArg(kernel_convolution_f, 2, sizeof(cl_mem), &d_gauss);
				clSetKernelArg(kernel_convolution_f, 3, sizeof(int32_t), &w_size_x);
				clSetKernelArg(kernel_convolution_f, 4, sizeof(int32_t), &w_size_y);
				clSetKernelArg(kernel_convolution_f, 5, sizeof(int32_t), &w_size_z);
				for (int ss = 0; ss < deblur_iterations; ss++) {
					status = clEnqueueCopyBuffer(commandQueues[0], d_mlem[0], d_mlem_apu_kolmas, 0, 0, im_dim * sizeof(float), 0, NULL, NULL);
					clFinish(commandQueues[0]);
					clSetKernelArg(kernel_convolution_f, 0, sizeof(cl_mem), &d_mlem_blurred[0]);
					clSetKernelArg(kernel_convolution_f, 1, sizeof(cl_mem), &d_mlem_apu);
					status = clEnqueueNDRangeKernel(commandQueues[0], kernel_convolution_f, 3ULL, NULL, convn_size, NULL, 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						mexPrintf("Failed to launch the convolution kernel\n");
						return;
					}
					clFinish(commandQueues[0]);
					clSetKernelArg(kernel_vectorDiv, 0, sizeof(cl_mem), &d_mlem_apu);
					clSetKernelArg(kernel_vectorDiv, 1, sizeof(cl_mem), &d_mlem_apu_kolmas);
					status = clEnqueueNDRangeKernel(commandQueues[0], kernel_vectorDiv, 1ULL, NULL, &vector_size, NULL, 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						mexPrintf("Failed to launch the division kernel\n");
						return;
					}
					clFinish(commandQueues[0]);
					clSetKernelArg(kernel_convolution_f, 0, sizeof(cl_mem), &d_mlem_apu_kolmas);
					clSetKernelArg(kernel_convolution_f, 1, sizeof(cl_mem), &d_mlem_apu);
					status = clEnqueueNDRangeKernel(commandQueues[0], kernel_convolution_f, 3ULL, NULL, convn_size, NULL, 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						mexPrintf("Failed to launch the convolution kernel\n");
						return;
					}
					clFinish(commandQueues[0]);
					clSetKernelArg(kernel_vectorMult, 0, sizeof(cl_mem), &d_mlem_apu);
					clSetKernelArg(kernel_vectorMult, 1, sizeof(cl_mem), &d_mlem_blurred[0]);
					status = clEnqueueNDRangeKernel(commandQueues[0], kernel_vectorMult, 1ULL, NULL, &vector_size, NULL, 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						mexPrintf("Failed to launch the multiplication kernel\n");
						return;
					}
					clFinish(commandQueues[0]);
					//status = clEnqueueCopyBuffer(commandQueues[0], d_mlem_apu_toka, d_mlem_blurred[0], 0, 0, im_dim * sizeof(float), NULL, NULL, NULL);
				}
				clFinish(commandQueues[0]);
				status = clEnqueueReadBuffer(commandQueues[0], d_mlem_blurred[0], CL_TRUE, 0, sizeof(float) * im_dim, &ele_ml[im_dim * (iter + 1u)], 0, NULL, NULL);
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
		if (verbose) {
			mexPrintf("Time step %u complete\n", tt + 1u);
			mexEvalString("pause(.0001);");
		}
		clFinish(commandQueues[0]);
		mxSetCell(cell, static_cast<mwIndex>(tt), mxDuplicateArray(mlem));
	}

	// Release memory
	for (cl_uint kk = 0u; kk < subsets; kk++) {
		for (cl_uint i = 0u; i < num_devices_context; i++) {
			if (kk == 0u) {
				if (i == 0 && use_psf) {
					clReleaseMemObject(d_gauss);
				}
				clReleaseMemObject(d_reko_type[i]);
				clReleaseMemObject(d_z[i]);
				clReleaseMemObject(d_x[i]);
				clReleaseMemObject(d_y[i]);
				clReleaseMemObject(d_xcenter[i]);
				clReleaseMemObject(d_ycenter[i]);
				clReleaseMemObject(d_zcenter[i]);
				clReleaseMemObject(d_V[i]);
				clReleaseMemObject(d_atten[i]);
				clReleaseMemObject(d_pseudos[i]);
				if (compute_norm_matrix == 1u)
					clReleaseMemObject(d_Summ[i]);
				clReleaseMemObject(d_rhs[i]);
				clReleaseMemObject(d_mlem[i]);
				clReleaseMemObject(d_mlem_blurred[i]);
				if (i < num_devices_context - 1) {
					clReleaseMemObject(d0_Summ[i]);
					clReleaseMemObject(d0_rhs[i]);
				}
			}
			if (compute_norm_matrix == 0u)
				clReleaseMemObject(d_Summ[kk * num_devices_context + i]);
			clReleaseMemObject(d_xyindex[kk * num_devices_context + i]);
			clReleaseMemObject(d_zindex[kk * num_devices_context + i]);
			clReleaseMemObject(d_L[kk * num_devices_context + i]);
			clReleaseMemObject(d_lor[kk * num_devices_context + i]);
			clReleaseMemObject(d_Sino[kk * num_devices_context + i]);
			clReleaseMemObject(d_sc_ra[kk * num_devices_context + i]);
			clReleaseMemObject(d_norm[kk * num_devices_context + i]);
			clReleaseMemObject(d_scat[kk * num_devices_context + i]);
		}
	}
	for (cl_uint i = 0u; i < num_devices_context; i++) {
		status = clFlush(commandQueues[i]);
		if (status != CL_SUCCESS) {
			getErrorString(status);
		}
		status = clFinish(commandQueues[i]);
		if (status != CL_SUCCESS) {
			getErrorString(status);
		}
	}
	mxDestroyArray(mlem);
	return;
}