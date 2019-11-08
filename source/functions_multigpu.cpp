/**************************************************************************
* All the functions needed for the matrix-free OpenCL image reconstruction
* in the multi-GPU/device case. First is the implementation 3 (OSEM_MLEM),
* second is the forward/backward projections (reconstruction_f_b_proj). 
* Lastly there are functions for context creation for both multi-device
* and single device case (clGetPlatformsContextSingle) and a function for
* program building and command queue creation (ClBuildProgramGetQueues).
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
#include "functions_multigpu.hpp"

// Implementation 3
void OSEM_MLEM(const cl_uint &num_devices_context, const float kerroin, const int cpu_device, const cl_context &context, const cl_command_queue *commandQueues,
	const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, const mxArray* Sin, const mxArray* sc_ra, const uint32_t Nx, 
	const uint32_t Ny, const uint32_t Nz, const uint32_t Niter, const mxArray* options, const float dx, const float dy, const float dz, const float bx, 
	const float by, const float bz, const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices, const uint32_t* pituus, 
	const size_t koko_l, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t size_x, const uint16_t TotSinos, const bool verbose, 
	const uint32_t randoms_correction, const uint32_t attenuation_correction, const uint32_t normalization, const float* atten, const size_t size_atten, 
	const float* norm, const size_t size_norm, const uint32_t subsets, const float epps, const uint32_t Nt,	const uint32_t* pseudos, const uint32_t det_per_ring, 
	const uint32_t prows, const uint16_t* L, const uint8_t raw, const size_t size_z, const uint32_t im_dim,	const cl_kernel &kernel, const cl_kernel &kernel_sum, 
	const cl_kernel &kernel_mlem, const size_t numel_x, const float tube_width, const float crystal_size_z, const float* x_center,
	const float* y_center, const float* z_center, const size_t size_center_x, const size_t size_center_y, const size_t size_center_z, const bool atomic_64bit, 
	const cl_uchar compute_norm_matrix, const bool precompute, const int32_t dec, const uint32_t projector_type, const uint16_t n_rays, const float cr_pz, 
	mxArray* cell, const bool osem_bool) {

	cl_int status = CL_SUCCESS;
	cl_float zero = 0.f;
	cl_short zeroL = 0;
	cl_ulong zerou = 0u;
	// Number of voxels in a single image
	const uint32_t Nxy = Nx * Ny;

	const mwSize dimmi[2] = { static_cast<mwSize>(Nx * Ny * Nz), static_cast<mwSize>(Niter + 1) };

	// Output matrix
	mxArray* mlem = mxCreateNumericArray(2, dimmi, mxSINGLE_CLASS, mxREAL);
	float* ele_ml = (float*)mxGetData(mlem);

	// Distance between rays in multi-ray Siddon
	const float dc_z = cr_pz / 3.f;

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

	// Memory allocation
	cl_mem* d0_rhs, * d0_Summ, * d_Summ;
	cl_mem* d_z = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_x = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_y = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_xcenter = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_ycenter = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_zcenter = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_atten = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_pseudos = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_rhs = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_mlem = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	if (compute_norm_matrix == 0u)
		d_Summ = (cl_mem*)malloc(subsets * num_devices_context * sizeof(cl_mem));
	else
		d_Summ = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_Sino = (cl_mem*)malloc((subsets * num_devices_context) * sizeof(cl_mem));
	cl_mem* d_sc_ra = (cl_mem*)malloc((subsets * num_devices_context) * sizeof(cl_mem));
	cl_mem* d_norm = (cl_mem*)malloc((subsets * num_devices_context) * sizeof(cl_mem));
	cl_mem* d_lor = (cl_mem*)malloc((subsets * num_devices_context) * sizeof(cl_mem));
	cl_mem* d_xyindex = (cl_mem*)malloc((subsets * num_devices_context) * sizeof(cl_mem));
	cl_mem* d_zindex = (cl_mem*)malloc((subsets * num_devices_context) * sizeof(cl_mem));
	cl_mem* d_L = (cl_mem*)malloc((subsets * num_devices_context) * sizeof(cl_mem));
	if (num_devices_context > 1u) {
		d0_rhs = (cl_mem*)malloc((num_devices_context - 1u) * sizeof(cl_mem));
		d0_Summ = (cl_mem*)malloc((num_devices_context - 1u) * sizeof(cl_mem));
	}

	// Create the necessary buffers
	for (cl_uint kk = 0u; kk < subsets; kk++) {
		for (cl_uint i = 0u; i < num_devices_context; i++) {
			if (kk == 0u) {
				d_z[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_z, NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				d_x[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * numel_x, NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				d_y[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * numel_x, NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				d_xcenter[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_center_x, NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				d_ycenter[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_center_y, NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				d_zcenter[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_center_z, NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				d_atten[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_atten, NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				d_pseudos[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t) * prows, NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				d_mlem[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * im_dim, NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				if (atomic_64bit) {
					d_rhs[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_ulong) * im_dim, NULL, &status);
					if (status != CL_SUCCESS) {
						std::cerr << getErrorString(status) << std::endl;
						return;
					}
					if (compute_norm_matrix == 1u) {
						d_Summ[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_ulong) * im_dim, NULL, &status);
						if (status != CL_SUCCESS) {
							std::cerr << getErrorString(status) << std::endl;
							return;
						}
					}
					if (i < num_devices_context - 1) {
						d0_rhs[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_ulong) * im_dim, NULL, &status);
						if (status != CL_SUCCESS) {
							std::cerr << getErrorString(status) << std::endl;
							return;
						}
						d0_Summ[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_ulong) * im_dim, NULL, &status);
						if (status != CL_SUCCESS) {
							std::cerr << getErrorString(status) << std::endl;
							return;
						}
					}
				}
				else {
					d_rhs[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * im_dim, NULL, &status);
					if (status != CL_SUCCESS) {
						std::cerr << getErrorString(status) << std::endl;
						return;
					}
					if (compute_norm_matrix == 1u) {
						d_Summ[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * im_dim, NULL, &status);
						if (status != CL_SUCCESS) {
							std::cerr << getErrorString(status) << std::endl;
							return;
						}
					}
					if (i < num_devices_context - 1) {
						d0_rhs[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * im_dim, NULL, &status);
						if (status != CL_SUCCESS) {
							std::cerr << getErrorString(status) << std::endl;
							return;
						}
						d0_Summ[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * im_dim, NULL, &status);
						if (status != CL_SUCCESS) {
							std::cerr << getErrorString(status) << std::endl;
							return;
						}
					}
				}
				const cl_mem siirrettavat[10] = { d_z[i] , d_x[i], d_y[i], d_xcenter[i], d_ycenter[i], d_zcenter[i], d_atten[i], d_pseudos[i], d_mlem[i], 
					d_rhs[i] };
				status = clEnqueueMigrateMemObjects(commandQueues[i], 10, siirrettavat, 0, 0, nullptr, nullptr);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
			}
			d_Sino[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * length[kk * num_devices_context + i], NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
			if (normalization == 1u) {
				d_norm[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * length[kk * num_devices_context + i], NULL, 
					&status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
			}
			else {
				d_norm[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
			}
			if (randoms_correction == 1u) {
				d_sc_ra[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * length[kk * num_devices_context + i], NULL, 
					&status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
			}
			else {
				d_sc_ra[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
			}
			if (precompute)
				d_lor[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk * num_devices_context + i], NULL, 
					&status);
			else
				d_lor[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
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
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
			}
			if (raw) {
				d_xyindex[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t), NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				d_zindex[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				d_L[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk * num_devices_context + i] * 2, NULL, 
					&status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
			}
			else {
				d_xyindex[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t) * length[kk * num_devices_context + i], 
					NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				d_zindex[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk * num_devices_context + i], 
					NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				d_L[kk * num_devices_context + i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
			}
			const cl_mem siirrettavat[7] = { d_Sino[kk * num_devices_context + i], d_norm[kk * num_devices_context + i], d_sc_ra[kk * num_devices_context + i], 
				d_lor[kk * num_devices_context + i], d_xyindex[kk * num_devices_context + i], d_zindex[kk * num_devices_context + i], 
				d_L[kk * num_devices_context + i]};
			status = clEnqueueMigrateMemObjects(commandQueues[i], 7, siirrettavat, 0, 0, nullptr, nullptr);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
		}
	}

	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return;
	}

	// Fill the buffers with host data
	for (cl_uint kk = 0u; kk < subsets; kk++) {
		for (cl_uint i = 0u; i < num_devices_context; i++) {
			if (kk == 0u) {
				status = clEnqueueWriteBuffer(commandQueues[i], d_x[i], CL_FALSE, 0, sizeof(float) * numel_x, x, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_y[i], CL_FALSE, 0, sizeof(float) * numel_x, y, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_xcenter[i], CL_FALSE, 0, sizeof(float) * size_center_x, x_center, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_ycenter[i], CL_FALSE, 0, sizeof(float) * size_center_y, y_center, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_zcenter[i], CL_FALSE, 0, sizeof(float) * size_center_z, z_center, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_z[i], CL_FALSE, 0, sizeof(float) * size_z, z_det, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_atten[i], CL_FALSE, 0, sizeof(float) * size_atten, atten, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_pseudos[i], CL_FALSE, 0, sizeof(uint32_t) * prows, pseudos, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_mlem[i], CL_FALSE, 0, sizeof(float) * im_dim, (float*)mxGetData(mxGetField(options, 0, "x0")), 
					0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}

			}
			if (precompute)
				status = clEnqueueWriteBuffer(commandQueues[i], d_lor[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint16_t) * 
					length[kk * num_devices_context + i], &lor1[cumsum[kk * num_devices_context + i]], 0, NULL, NULL);
			else
				status = clEnqueueWriteBuffer(commandQueues[i], d_lor[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint16_t), lor1, 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
			if (normalization == 1u) {
				status = clEnqueueWriteBuffer(commandQueues[i], d_norm[kk * num_devices_context + i], CL_FALSE, 0, sizeof(cl_float) * 
					length[kk * num_devices_context + i], &norm[cumsum[kk * num_devices_context + i]], 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
			}
			else {
				status = clEnqueueWriteBuffer(commandQueues[i], d_norm[kk * num_devices_context + i], CL_FALSE, 0, sizeof(cl_float) * size_norm, norm, 0, 
					NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
			}

			if (raw) {
				status = clEnqueueWriteBuffer(commandQueues[i], d_xyindex[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint32_t), xy_index, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_zindex[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint16_t), z_index, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_L[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint16_t) * 
					length[kk * num_devices_context + i] * 2, &L[cumsum[kk * num_devices_context + i] * 2], 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
			}
			else {
				status = clEnqueueWriteBuffer(commandQueues[i], d_xyindex[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint32_t) * 
					length[kk * num_devices_context + i], &xy_index[cumsum[kk * num_devices_context + i]], 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_zindex[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint16_t) * 
					length[kk * num_devices_context + i], &z_index[cumsum[kk * num_devices_context + i]], 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[i], d_L[kk * num_devices_context + i], CL_FALSE, 0, sizeof(uint16_t), L, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
			}

			status = clFlush(commandQueues[i]);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
		}
	}

	size_t sum_dim = static_cast<size_t>(im_dim);

	// Temporary vectors for multi-device case (output data are stored here, then transferred to the primary device)
	std::vector<std::vector<float>> testi_summ(num_devices_context - 1u, std::vector<float>(im_dim));
	std::vector<std::vector<float>> testi_rhs(num_devices_context - 1u, std::vector<float>(im_dim));

	for (cl_uint i = 0u; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

	cl_uint kernelInd = 0U;

	// Set the constant kernel arguments
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
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &attenuation_correction);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &normalization);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &randoms_correction);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &det_per_ring);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint8_t), &raw);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &prows);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &Nxy);
	if (projector_type == 2u) {
		clSetKernelArg(kernel, kernelInd++, sizeof(float), &tube_width);
		clSetKernelArg(kernel, kernelInd++, sizeof(float), &crystal_size_z);
		clSetKernelArg(kernel, kernelInd++, sizeof(int32_t), &dec);
	}
	else if (projector_type == 1u && !precompute) {
		clSetKernelArg(kernel, kernelInd++, sizeof(float), &dc_z);
		clSetKernelArg(kernel, kernelInd++, sizeof(uint16_t), &n_rays);
	}


	clSetKernelArg(kernel_mlem, 3, sizeof(uint32_t), &im_dim);
	clSetKernelArg(kernel_mlem, 4, sizeof(float), &epps);

	cl_uint event_count = (3u);

	// Fixed local size
	const size_t local_size = 64ULL;

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
					std::cerr << getErrorString(status) << std::endl;
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
		float * Sino = (float*)mxGetData(mxGetCell(Sin, static_cast<mwIndex>(tt)));

		for (cl_uint i = 0u; i < num_devices_context; i++) {
			clFinish(commandQueues[i]);
		}

		// Transfer data from host to device
		for (cl_uint kk = 0u; kk < subsets; kk++) {
			for (cl_uint i = 0u; i < num_devices_context; i++) {
				status = clEnqueueWriteBuffer(commandQueues[i], d_Sino[kk * num_devices_context + i], CL_FALSE, 0, sizeof(float) * 
					length[kk * num_devices_context + i], &Sino[cumsum[kk * num_devices_context + i]], 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}
				// Randoms
				if (randoms_correction == 1u) {
					float* S_R = (float*)mxGetData(mxGetCell(sc_ra, tt));
					status = clEnqueueWriteBuffer(commandQueues[i], d_sc_ra[kk * num_devices_context + i], CL_FALSE, 0, sizeof(float) * 
						length[kk * num_devices_context + i], &S_R[cumsum[kk * num_devices_context + i]], 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						std::cerr << getErrorString(status) << std::endl;
						return;
					}
				}
				else {
					status = clEnqueueFillBuffer(commandQueues[i], d_sc_ra[kk * num_devices_context + i], &zero, sizeof(cl_float), 0, sizeof(cl_float), 0, NULL, NULL);
					//status = clEnqueueWriteBuffer(commandQueues[i], d_sc_ra[kk * num_devices_context + i], CL_FALSE, 0, sizeof(float), S_R, 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						std::cerr << getErrorString(status) << std::endl;
						return;
					}
				}
				status = clFlush(commandQueues[i]);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
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
					std::cerr << getErrorString(status) << std::endl;
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
							std::cerr << getErrorString(status) << std::endl;
							return;
						}
						status = clEnqueueFillBuffer(commandQueues[i], d_rhs[i], &zerou, sizeof(cl_ulong), 0, sizeof(cl_ulong) * im_dim, 0, NULL, NULL);
						if (status != CL_SUCCESS) {
							std::cerr << getErrorString(status) << std::endl;
							return;
						}
					}
					else {
						if (compute_norm_matrix == 1u) {
							status = clEnqueueFillBuffer(commandQueues[i], d_Summ[i], &zero, sizeof(cl_float), 0, sizeof(cl_float) * im_dim, 0, NULL, NULL);
						}
						if (status != CL_SUCCESS) {
							std::cerr << getErrorString(status) << std::endl;
							return;
						}
						status = clEnqueueFillBuffer(commandQueues[i], d_rhs[i], &zero, sizeof(cl_float), 0, sizeof(cl_float) * im_dim, 0, NULL, NULL);
						if (status != CL_SUCCESS) {
							std::cerr << getErrorString(status) << std::endl;
							return;
						}
					}

					status = clFlush(commandQueues[i]);
					if (status != CL_SUCCESS) {
						std::cerr << getErrorString(status) << std::endl;
						return;
					}
				}

				for (cl_uint i = 0u; i < num_devices_context; i++) {
					clFinish(commandQueues[i]);
				}

				// Loop through the devices
				for (cl_uint i = 0u; i < num_devices_context; i++) {

					cl_uint kernelIndSubIter = kernelInd;

					// Make sure the global size is divisible with 64
					const size_t global_size = length[osa_iter * num_devices_context + i] + (local_size - length[osa_iter * num_devices_context + i] % local_size);

					// The "true" global size
					const uint64_t m_size = static_cast<uint64_t>(length[osa_iter * num_devices_context + i]);

					// Set dynamic kernel arguments
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_atten[i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_norm[osa_iter * num_devices_context + i]);
					if (compute_norm_matrix == 0u)
						clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_Summ[osa_iter * num_devices_context + i]);
					else
						clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_Summ[i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_lor[osa_iter * num_devices_context + i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_pseudos[i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_x[i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_y[i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_z[i]);
					if (projector_type == 2u) {
						clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_xcenter[i]);
						clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_ycenter[i]);
						clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_zcenter[i]);
					}
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_xyindex[osa_iter * num_devices_context + i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_zindex[osa_iter * num_devices_context + i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_L[osa_iter * num_devices_context + i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_Sino[osa_iter * num_devices_context + i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_sc_ra[osa_iter * num_devices_context + i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_mlem[i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_rhs[i]);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_uchar), &no_norm);
					clSetKernelArg(kernel, kernelIndSubIter++, sizeof(uint64_t), &m_size);
					// Compute the RHS and normalization constant
					status = clEnqueueNDRangeKernel(commandQueues[i], kernel, 1u, NULL, &global_size, &local_size, 0u, NULL, &events[i][0]);
					if (status != CL_SUCCESS) {
						std::cerr << getErrorString(status) << std::endl;
						for (cl_uint kk = 0u; kk < subsets; kk++) {
							for (cl_uint i = 0u; i < num_devices_context; i++) {
								if (kk == 0) {
									clReleaseMemObject(d_z[i]);
									clReleaseMemObject(d_x[i]);
									clReleaseMemObject(d_y[i]);
									clReleaseMemObject(d_xcenter[i]);
									clReleaseMemObject(d_ycenter[i]);
									clReleaseMemObject(d_zcenter[i]);
									clReleaseMemObject(d_atten[i]);
									clReleaseMemObject(d_pseudos[i]);
									if (compute_norm_matrix == 1u)
										clReleaseMemObject(d_Summ[i]);
									clReleaseMemObject(d_rhs[i]);
									clReleaseMemObject(d_mlem[i]);
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
								status = clEnqueueReadBuffer(commandQueues[i], d_Summ[i], CL_TRUE, 0, sizeof(cl_ulong) * im_dim, testi_summ[i - 1u].data(), 1, 
									&events[i][0], NULL);
							else if (compute_norm_matrix == 0u && iter == 0u)
								status = clEnqueueReadBuffer(commandQueues[i], d_Summ[osa_iter * num_devices_context + i], CL_TRUE, 0, sizeof(cl_ulong) * im_dim, 
									testi_summ[i - 1u].data(), 1, &events[i][0], NULL);
							if (status != CL_SUCCESS) {
								std::cerr << getErrorString(status) << std::endl;
								return;
							}
							if ((compute_norm_matrix == 0u && iter == 0u) || compute_norm_matrix == 1u) {
								status = clEnqueueWriteBuffer(commandQueues[0], d0_Summ[i - 1u], CL_FALSE, 0, sizeof(cl_ulong) * im_dim, testi_summ[i - 1u].data(), 
									0, NULL, NULL);
								if (status != CL_SUCCESS) {
									std::cerr << getErrorString(status) << std::endl;
									return;
								}
							}
							status = clEnqueueReadBuffer(commandQueues[i], d_rhs[i], CL_TRUE, 0, sizeof(cl_ulong) * im_dim, testi_rhs[i - 1u].data(), 1, 
								&events[i][0], NULL);
							if (status != CL_SUCCESS) {
								std::cerr << getErrorString(status) << std::endl;
								return;
							}
							status = clEnqueueWriteBuffer(commandQueues[0], d0_rhs[i - 1u], CL_FALSE, 0, sizeof(cl_ulong) * im_dim, testi_rhs[i - 1u].data(), 0, 
								NULL, NULL);
							if (status != CL_SUCCESS) {
								std::cerr << getErrorString(status) << std::endl;
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
								std::cerr << getErrorString(status) << std::endl;
								return;
							}
							if ((compute_norm_matrix == 0u && iter == 0u) || compute_norm_matrix == 1u) {
								status = clEnqueueWriteBuffer(commandQueues[0], d0_Summ[i - 1u], CL_FALSE, 0, sizeof(float) * im_dim, testi_summ[i - 1u].data(), 
									0, NULL, NULL);
								if (status != CL_SUCCESS) {
									std::cerr << getErrorString(status) << std::endl;
									return;
								}
							}
							status = clEnqueueReadBuffer(commandQueues[i], d_rhs[i], CL_TRUE, 0, sizeof(float) * im_dim, testi_rhs[i - 1u].data(), 1, 
								&events[i][0], NULL);
							if (status != CL_SUCCESS) {
								std::cerr << getErrorString(status) << std::endl;
								return;
							}
							status = clEnqueueWriteBuffer(commandQueues[0], d0_rhs[i - 1u], CL_FALSE, 0, sizeof(float) * im_dim, testi_rhs[i - 1u].data(), 0, 
								NULL, NULL);
							if (status != CL_SUCCESS) {
								std::cerr << getErrorString(status) << std::endl;
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
							std::cerr << getErrorString(status) << std::endl;
							return;
						}
						clFinish(commandQueues[0]);
					}
				}

				//for (cl_uint i = 1u; i < num_devices_context; i++)
				//	clWaitForEvents(1, &events[i][3]);

				// Compute MLEM/OSEM
				if (compute_norm_matrix == 1u)
					clSetKernelArg(kernel_mlem, 0, sizeof(cl_mem), &d_Summ[0]);
				else
					clSetKernelArg(kernel_mlem, 0, sizeof(cl_mem), &d_Summ[osa_iter * num_devices_context]);
				clSetKernelArg(kernel_mlem, 1, sizeof(cl_mem), &d_rhs[0]);
				clSetKernelArg(kernel_mlem, 2, sizeof(cl_mem), &d_mlem[0]);
				status = clEnqueueNDRangeKernel(commandQueues[0], kernel_mlem, 1, NULL, &sum_dim, NULL, 0, NULL, &events[0][1]);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					for (cl_uint kk = 0u; kk < subsets; kk++) {
						for (cl_uint i = 0u; i < num_devices_context; i++) {
							if (kk == 0) {
								clReleaseMemObject(d_z[i]);
								clReleaseMemObject(d_x[i]);
								clReleaseMemObject(d_y[i]);
								clReleaseMemObject(d_xcenter[i]);
								clReleaseMemObject(d_ycenter[i]);
								clReleaseMemObject(d_zcenter[i]);
								clReleaseMemObject(d_atten[i]);
								clReleaseMemObject(d_pseudos[i]);
								if (compute_norm_matrix == 1u)
									clReleaseMemObject(d_Summ[i]);
								clReleaseMemObject(d_rhs[i]);
								clReleaseMemObject(d_mlem[i]);
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
							std::cerr << getErrorString(status) << std::endl;
							return;
						}
					}
					else {
						status = clEnqueueReadBuffer(commandQueues[0], d_Summ[osa_iter * num_devices_context], CL_FALSE, 0, sizeof(float) * im_dim, 
							testi_summ[0].data(), 1, &events[0][1], &summ_event);
						if (status != CL_SUCCESS) {
							std::cerr << getErrorString(status) << std::endl;
							return;
						}
					}
					for (cl_uint i = 1u; i < num_devices_context; i++) {
						if (atomic_64bit) {
							status = clEnqueueWriteBuffer(commandQueues[0], d_Summ[osa_iter * num_devices_context + i], CL_FALSE, 0, sizeof(cl_ulong) * im_dim, 
								testi_summ[0].data(), 1, &summ_event, NULL);
							if (status != CL_SUCCESS) {
								std::cerr << getErrorString(status) << std::endl;
								return;
							}
						}
						else {
							status = clEnqueueWriteBuffer(commandQueues[0], d_Summ[osa_iter * num_devices_context + i], CL_FALSE, 0, sizeof(float) * im_dim, 
								testi_summ[0].data(), 1, &summ_event, NULL);
							if (status != CL_SUCCESS) {
								std::cerr << getErrorString(status) << std::endl;
								return;
							}
						}
					}
				}

				// Transfer estimate data back to host
				status = clEnqueueReadBuffer(commandQueues[0], d_mlem[0], CL_TRUE, 0, sizeof(float) * im_dim, &ele_ml[im_dim * (iter + 1u)], 1, &events[0][1], 
					&events[0][2]);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return;
				}

				// Transfer estimate data to secondary devices
				for (cl_uint i = 1u; i < num_devices_context; i++) {
					status = clEnqueueWriteBuffer(commandQueues[i], d_mlem[i], CL_FALSE, 0, sizeof(cl_float) * im_dim, &ele_ml[im_dim * (iter + 1u)], 1, 
						&events[0][2], NULL);
					if (status != CL_SUCCESS) {
						std::cerr << getErrorString(status) << std::endl;
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
				}

				if (compute_norm_matrix == 0u && iter == 0u && num_devices_context > 1u)
					clReleaseEvent(summ_event);
				if (verbose && subsets > 1u) {
					mexPrintf("Sub-iteration %u complete\n", osa_iter + 1u);
					mexEvalString("pause(.0001);");
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
				clReleaseMemObject(d_z[i]);
				clReleaseMemObject(d_x[i]);
				clReleaseMemObject(d_y[i]);
				clReleaseMemObject(d_xcenter[i]);
				clReleaseMemObject(d_ycenter[i]);
				clReleaseMemObject(d_zcenter[i]);
				clReleaseMemObject(d_atten[i]);
				clReleaseMemObject(d_pseudos[i]);
				if (compute_norm_matrix == 1u)
					clReleaseMemObject(d_Summ[i]);
				clReleaseMemObject(d_rhs[i]);
				clReleaseMemObject(d_mlem[i]);
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
		}
	}
	mxDestroyArray(mlem);
}

// Forward/backward projections
void f_b_project(const cl_uint &num_devices_context, const float kerroin, const int cpu_device, const cl_context &context, const cl_command_queue *commandQueues,
	const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, const float* rhs, const mxArray* sc_ra, const uint32_t Nx, 
	const uint32_t Ny, const uint32_t Nz, const float dx, const float dy, const float dz, const float bx, const float by, const float bz, const float bzb, 
	const float maxxx, const float maxyy, const float zmax, const float NSlices, const uint32_t* pituus, const size_t koko_l, const uint32_t* xy_index, 
	const uint16_t* z_index, const uint32_t size_x, const uint16_t TotSinos, const bool verbose, const uint32_t randoms_correction, 
	const uint32_t attenuation_correction, const uint32_t normalization, const float* atten, const size_t size_atten, const float* norm, const size_t size_norm,
	const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L, const uint8_t raw, const size_t size_z, const uint32_t im_dim,
	const cl_kernel &kernel_sum, const cl_kernel &kernel, float *output, float* normalizer, const size_t size_rhs, const bool no_norm, const size_t numel_x, 
	const float tube_width, const float crystal_size_z, const float* x_center, const float* y_center, const float* z_center, const size_t size_center_x, 
	const size_t size_center_y, const size_t size_center_z, const bool precompute, const int32_t dec, const uint32_t projector_type, const uint16_t n_rays, 
	const float cr_pz) {

	const uint32_t Nxy = Nx * Ny;
	cl_int status = CL_SUCCESS;
	cl_float zero = 0.f;
	cl_short zeroL = 0;

	// If 1, then the forward projection is computed
	uint8_t fp = size_rhs == im_dim;

	// Distance between rays in multi-ray Siddon
	const float dc_z = cr_pz / 3.f;

	size_t size_output;

	if (fp)
		size_output = pituus[0];
	else
		size_output = static_cast<size_t>(im_dim);

	// Essentially the same as above, but without subsets
	std::vector<size_t> cumsum((num_devices_context + 1U), 0);
	std::vector<size_t> length(num_devices_context, 0);

	size_t meas_per_gpu;
	size_t meas_per_cpu = 0ULL;
	if (cpu_device >= 0) {
		meas_per_cpu = static_cast<size_t>(static_cast<float>(pituus[0]) / (kerroin + 1.f));
		meas_per_gpu = (pituus[0] - meas_per_cpu) / static_cast<size_t>(num_devices_context - 1U);
	}
	else
		meas_per_gpu = pituus[0] / static_cast<size_t>(num_devices_context);
	cumsum[0] = 0ULL;

	for (cl_uint i = 0U; i < num_devices_context; i++) {
		if (i == num_devices_context - 1U) {
			if (i == cpu_device) {
				length[i] = meas_per_cpu;
				length[i - 1] = meas_per_gpu + (pituus[0] - meas_per_gpu * static_cast<size_t>(num_devices_context - 1U) - meas_per_cpu);
				cumsum[i] = cumsum[i - 1U] + length[i - 1U];
			}
			else {
				if (cpu_device < 0)
					length[i] = meas_per_gpu + (pituus[0] - meas_per_gpu * static_cast<size_t>(num_devices_context));
				else
					length[i] = meas_per_gpu + (pituus[0] - meas_per_gpu * static_cast<size_t>(num_devices_context - 1U) - meas_per_cpu);
			}
		}
		else {
			if (i == cpu_device) {
				length[i] = meas_per_cpu;
			}
			else
				length[i] = meas_per_gpu;
		}
		cumsum[i + 1U] = cumsum[i] + length[i];
	}

	cl_mem* d0_output, * d0_Summ;
	cl_mem* d_z = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_x = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_y = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_xcenter = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_ycenter = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_zcenter = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_atten = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_pseudos = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_rhs = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_output = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_Summ = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_lor = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	if (num_devices_context > 1u) {
		d0_output = (cl_mem*)malloc((num_devices_context - 1) * sizeof(cl_mem));
		d0_Summ = (cl_mem*)malloc((num_devices_context - 1) * sizeof(cl_mem));
	}
	cl_mem* d_xyindex = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_zindex = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_L = (cl_mem*)malloc(num_devices_context * sizeof(cl_mem));
	cl_mem* d_sc_ra = (cl_mem*)malloc((num_devices_context) * sizeof(cl_mem));
	cl_mem* d_norm = (cl_mem*)malloc((num_devices_context) * sizeof(cl_mem));

	// Create the necessary buffers
	for (cl_uint i = 0U; i < num_devices_context; i++) {
		d_z[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_z, NULL, &status);
		d_x[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * numel_x, NULL, &status);
		d_y[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * numel_x, NULL, &status);
		d_xcenter[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_center_x, NULL, &status);
		d_ycenter[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_center_y, NULL, &status);
		d_zcenter[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_center_z, NULL, &status);
		d_atten[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_atten, NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		d_pseudos[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t) * prows, NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		d_output[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * size_output, NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		d_rhs[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * size_rhs, NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		d_Summ[i] = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * im_dim, NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		if (randoms_correction == 1u) {
			d_sc_ra[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * length[i], NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
		}
		else {
			d_sc_ra[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
		}
		if (normalization == 1u) {
			d_norm[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_float) * length[i], NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
		}
		else {
			d_norm[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_float), NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
		}
		if (precompute)
			d_lor[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[i], NULL, &status);
		else
			d_lor[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		if (raw) {
			d_xyindex[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t), NULL, &status);
			d_zindex[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			d_L[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[i], NULL, &status);
		}
		else {
			d_xyindex[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t) * length[i], NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
			d_zindex[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[i], NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
			d_L[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
		}
		if (i < num_devices_context - 1) {
			d0_output[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * size_output, NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
			d0_Summ[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * im_dim, NULL, &status);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
		}
	}


	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return;
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		status = clEnqueueWriteBuffer(commandQueues[i], d_x[i], CL_FALSE, 0, sizeof(float) * numel_x, x, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		status = clEnqueueWriteBuffer(commandQueues[i], d_y[i], CL_FALSE, 0, sizeof(float) * numel_x, y, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		status = clEnqueueWriteBuffer(commandQueues[i], d_xcenter[i], CL_FALSE, 0, sizeof(float) * size_center_x, x_center, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		status = clEnqueueWriteBuffer(commandQueues[i], d_ycenter[i], CL_FALSE, 0, sizeof(float) * size_center_y, y_center, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		status = clEnqueueWriteBuffer(commandQueues[i], d_zcenter[i], CL_FALSE, 0, sizeof(float) * size_center_z, z_center, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		status = clEnqueueWriteBuffer(commandQueues[i], d_z[i], CL_FALSE, 0, sizeof(float) * size_z, z_det, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		status = clEnqueueWriteBuffer(commandQueues[i], d_atten[i], CL_FALSE, 0, sizeof(float) * size_atten, atten, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		status = clEnqueueWriteBuffer(commandQueues[i], d_pseudos[i], CL_FALSE, 0, sizeof(uint32_t) * prows, pseudos, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		status = clEnqueueWriteBuffer(commandQueues[i], d_rhs[i], CL_FALSE, 0, sizeof(float) * size_rhs, rhs, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}


		if (raw) {
			status = clEnqueueWriteBuffer(commandQueues[i], d_xyindex[i], CL_FALSE, 0, sizeof(uint32_t), xy_index, 0, NULL, NULL);
			status = clEnqueueWriteBuffer(commandQueues[i], d_zindex[i], CL_FALSE, 0, sizeof(uint16_t), z_index, 0, NULL, NULL);
			status = clEnqueueWriteBuffer(commandQueues[i], d_L[i], CL_FALSE, 0, sizeof(uint16_t) * length[i], &L[cumsum[i]], 0, NULL, NULL);
		}
		else {
			status = clEnqueueWriteBuffer(commandQueues[i], d_xyindex[i], CL_FALSE, 0, sizeof(uint32_t) * length[i], &xy_index[cumsum[i]], 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
			status = clEnqueueWriteBuffer(commandQueues[i], d_zindex[i], CL_FALSE, 0, sizeof(uint16_t) * length[i], &z_index[cumsum[i]], 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
			status = clEnqueueWriteBuffer(commandQueues[i], d_L[i], CL_FALSE, 0, sizeof(uint16_t), L, 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
		}
		if (precompute)
			status = clEnqueueWriteBuffer(commandQueues[i], d_lor[i], CL_FALSE, 0, sizeof(uint16_t) * length[i], &lor1[cumsum[i]], 0, NULL, NULL);
		else
			status = clEnqueueWriteBuffer(commandQueues[i], d_lor[i], CL_FALSE, 0, sizeof(uint16_t), lor1, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		if (normalization == 1u) {
			status = clEnqueueWriteBuffer(commandQueues[i], d_norm[i], CL_FALSE, 0, sizeof(cl_float) * length[i], &norm[cumsum[i]], 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
		}
		else {
			status = clEnqueueWriteBuffer(commandQueues[i], d_norm[i], CL_FALSE, 0, sizeof(cl_float), norm, 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
		}

		status = clFlush(commandQueues[i]);
	}

	uint64_t size_rhs_i = static_cast<uint64_t>(size_rhs);
	const size_t local_size = 64ULL;


	std::vector<std::vector<float>> testi_summ(num_devices_context - 1, std::vector<float>(im_dim));
	std::vector<std::vector<float>> testi_rhs(num_devices_context - 1, std::vector<float>(size_output));


	for (cl_uint i = 0; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

	cl_uint kernelInd = 0U;

	status = clSetKernelArg(kernel, kernelInd++, sizeof(uint8_t), &fp);
	status = clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &im_dim);
	status = clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &Nx);
	status = clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &Ny);
	status = clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &Nz);
	status = clSetKernelArg(kernel, kernelInd++, sizeof(float), &dz);
	status = clSetKernelArg(kernel, kernelInd++, sizeof(float), &dx);
	status = clSetKernelArg(kernel, kernelInd++, sizeof(float), &dy);
	status = clSetKernelArg(kernel, kernelInd++, sizeof(float), &bz);
	status = clSetKernelArg(kernel, kernelInd++, sizeof(float), &bx);
	status = clSetKernelArg(kernel, kernelInd++, sizeof(float), &by);
	clSetKernelArg(kernel, kernelInd++, sizeof(float), &bzb);
	clSetKernelArg(kernel, kernelInd++, sizeof(float), &maxxx);
	clSetKernelArg(kernel, kernelInd++, sizeof(float), &maxyy);
	clSetKernelArg(kernel, kernelInd++, sizeof(float), &zmax);
	clSetKernelArg(kernel, kernelInd++, sizeof(float), &NSlices);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &size_x);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint16_t), &TotSinos);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &attenuation_correction);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &normalization);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &randoms_correction);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &det_per_ring);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint8_t), &raw);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &prows);
	clSetKernelArg(kernel, kernelInd++, sizeof(cl_uchar), &no_norm);
	clSetKernelArg(kernel, kernelInd++, sizeof(uint32_t), &Nxy);
	if (projector_type == 2u) {
		clSetKernelArg(kernel, kernelInd++, sizeof(float), &tube_width);
		clSetKernelArg(kernel, kernelInd++, sizeof(float), &crystal_size_z);
		clSetKernelArg(kernel, kernelInd++, sizeof(int32_t), &dec);
	}
	else if (projector_type == 1u && !precompute) {
		clSetKernelArg(kernel, kernelInd++, sizeof(float), &dc_z);
		clSetKernelArg(kernel, kernelInd++, sizeof(uint16_t), &n_rays);
	}
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return;
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

	cl_uint event_count = (1u);


	cl_event** events = (cl_event * *)malloc(num_devices_context * sizeof(cl_event * *));
	for (cl_uint i = 0u; i < num_devices_context; i++) {
		events[i] = (cl_event*)malloc(event_count * sizeof(cl_event * *));
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		status = clEnqueueFillBuffer(commandQueues[i], d_Summ[i], &zero, sizeof(cl_float), 0, sizeof(cl_float) * im_dim, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		status = clEnqueueFillBuffer(commandQueues[i], d_output[i], &zero, sizeof(cl_float), 0, sizeof(cl_float) * size_output, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		if (randoms_correction == 1u) {
			float* S_R = (float*)mxGetData(mxGetCell(sc_ra, 0));
			status = clEnqueueWriteBuffer(commandQueues[i], d_sc_ra[i], CL_FALSE, 0, sizeof(float) * length[i], &S_R[cumsum[i]], 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
		}
		else {
			status = clEnqueueFillBuffer(commandQueues[i], d_sc_ra[i], &zero, sizeof(cl_float), 0, sizeof(cl_float), 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
		}

		//status = clFlush(commandQueues[i]);
		//if (status != CL_SUCCESS) {
		//	std::cerr << getErrorString(status) << std::endl;
		//	return;
		//}
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

	std::vector<cl_ulong> globals(num_devices_context, 0ULL);

	for (cl_uint i = 0; i < num_devices_context; i++) {

		cl_uint kernelIndSubIter = kernelInd;

		const size_t global_size = length[i] + (local_size - length[i] % local_size);
		const uint64_t m_size = static_cast<uint64_t>(length[i]);

		// These are needed to make sure the forward projection data is correctly combined
		if (i > 0)
			globals[i] = globals[i - 1ULL] + m_size;
		else
			globals[i] = m_size;

		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_atten[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_norm[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_Summ[i]);
		if (precompute)
			clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_lor[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_pseudos[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_x[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_y[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_z[i]);
		if (projector_type == 2u) {
			clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_xcenter[i]);
			clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_ycenter[i]);
			clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_zcenter[i]);
		}
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_xyindex[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_zindex[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_L[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_sc_ra[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_rhs[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_output[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(uint64_t), &m_size);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(uint64_t), &static_cast<uint64_t>(cumsum[i]));
		status = clEnqueueNDRangeKernel(commandQueues[i], kernel, 1, NULL, &global_size, &local_size, 0, NULL, &events[i][0]);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			for (cl_uint i = 0; i < num_devices_context; i++) {
				clReleaseMemObject(d_z[i]);
				clReleaseMemObject(d_x[i]);
				clReleaseMemObject(d_y[i]);
				clReleaseMemObject(d_xcenter[i]);
				clReleaseMemObject(d_ycenter[i]);
				clReleaseMemObject(d_zcenter[i]);
				clReleaseMemObject(d_atten[i]);
				clReleaseMemObject(d_norm[i]);
				clReleaseMemObject(d_pseudos[i]);
				clReleaseMemObject(d_Summ[i]);
				clReleaseMemObject(d_rhs[i]);
				clReleaseMemObject(d_output[i]);
				clReleaseMemObject(d_xyindex[i]);
				clReleaseMemObject(d_zindex[i]);
				clReleaseMemObject(d_L[i]);
				clReleaseMemObject(d_sc_ra[i]);
				clReleaseMemObject(d_lor[i]);
				if (i < num_devices_context - 1) {
					clReleaseMemObject(d0_Summ[i]);
					clReleaseMemObject(d0_output[i]);
				}
			}
			return;
		}
	}

	if (num_devices_context > 1u) {
		for (cl_uint i = 1; i < num_devices_context; i++) {
			status = clEnqueueReadBuffer(commandQueues[i], d_Summ[i], CL_TRUE, 0, sizeof(float) * im_dim, testi_summ[i - 1].data(), 1, &events[i][0], NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
			status = clEnqueueReadBuffer(commandQueues[i], d_output[i], CL_TRUE, 0, sizeof(float) * size_output, testi_rhs[i - 1].data(), 1, &events[i][0], NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
			status = clEnqueueWriteBuffer(commandQueues[0], d0_Summ[i - 1], CL_FALSE, 0, sizeof(float) * im_dim, testi_summ[i - 1].data(), 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
			status = clEnqueueWriteBuffer(commandQueues[0], d0_output[i - 1], CL_FALSE, 0, sizeof(float) * size_output, testi_rhs[i - 1].data(), 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				std::cerr << getErrorString(status) << std::endl;
				return;
			}
			for (cl_uint i = 0ULL; i < num_devices_context; i++) {
				clFinish(commandQueues[i]);
			}
		}
	}

	clWaitForEvents(1, &events[0][0]);

	for (cl_uint i = 1; i < num_devices_context; i++) {
		size_t summa_pituus;
		if (fp == 1)
			summa_pituus = length[i];
		else
			summa_pituus = size_rhs;
		clSetKernelArg(kernel_sum, 0, sizeof(cl_mem), &d0_Summ[i - 1]);
		clSetKernelArg(kernel_sum, 1, sizeof(cl_mem), &d_Summ[0]);
		clSetKernelArg(kernel_sum, 2, sizeof(cl_mem), &d0_output[i - 1]);
		clSetKernelArg(kernel_sum, 3, sizeof(cl_mem), &d_output[0]);
		clSetKernelArg(kernel_sum, 4, sizeof(uint64_t), &size_rhs_i);
		clSetKernelArg(kernel_sum, 5, sizeof(uint32_t), &im_dim);
		clSetKernelArg(kernel_sum, 6, sizeof(bool), &no_norm);
		clSetKernelArg(kernel_sum, 7, sizeof(uint64_t), &globals[i - 1ULL]);
		clSetKernelArg(kernel_sum, 8, sizeof(uint8_t), &fp);
		status = clEnqueueNDRangeKernel(commandQueues[0], kernel_sum, 1, NULL, &summa_pituus, NULL, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return;
		}
		clFinish(commandQueues[0]);
	}

	status = clEnqueueReadBuffer(commandQueues[0], d_output[0], CL_FALSE, 0, sizeof(float) * size_output, output, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return;
	}
	status = clEnqueueReadBuffer(commandQueues[0], d_Summ[0], CL_FALSE, 0, sizeof(float) * im_dim, normalizer, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return;
	}

	for (cl_uint i = 0u; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
		clReleaseEvent(events[i][0]);
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		clReleaseMemObject(d_z[i]);
		clReleaseMemObject(d_x[i]);
		clReleaseMemObject(d_y[i]);
		clReleaseMemObject(d_xcenter[i]);
		clReleaseMemObject(d_ycenter[i]);
		clReleaseMemObject(d_zcenter[i]);
		clReleaseMemObject(d_atten[i]);
		clReleaseMemObject(d_norm[i]);
		clReleaseMemObject(d_pseudos[i]);
		clReleaseMemObject(d_Summ[i]);
		clReleaseMemObject(d_rhs[i]);
		clReleaseMemObject(d_output[i]);
		clReleaseMemObject(d_xyindex[i]);
		clReleaseMemObject(d_zindex[i]);
		clReleaseMemObject(d_L[i]);
		clReleaseMemObject(d_sc_ra[i]);
		clReleaseMemObject(d_lor[i]);
		if (i < num_devices_context - 1) {
			clReleaseMemObject(d0_Summ[i]);
			clReleaseMemObject(d0_output[i]);
		}
	}
}

// Get the OpenCL context for the current platform
cl_int clGetPlatformsContext(const uint32_t device, const float kerroin, cl_context& context, size_t& size, int& cpu_device, 
	cl_uint& num_devices_context, cl_device_id* devices, bool& atomic_64bit, cl_uchar& compute_norm_matrix, const uint32_t Nxyz, 
	const uint32_t subsets, const uint8_t raw) {
	cl_int status = CL_SUCCESS;
	cl_uint num_platforms;
	cl_float mem_portions;
	// These variables are used to determine if the normalization constant is memory-wise possible to store
	if (raw == 1u)
		mem_portions = 0.1f;
	else
		mem_portions = 0.2f;
	cl_float image_bytes = static_cast<cl_float>(Nxyz) * 8.f;

	// Get the number of platforms
	status = clGetPlatformIDs(0, NULL, &num_platforms);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}

	cl_platform_id *platforms = new cl_platform_id[num_platforms];

	// Get the platform IDs
	status = clGetPlatformIDs(num_platforms, platforms, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		delete[] platforms;
		return status;
	}

	if (device > num_platforms) {
		std::cerr << "The specified platform number is greater than the available platform numbers!" << std::endl;
		delete[] platforms;
		return status;
	}

	// Get context properties from the chosen platform
	cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platforms[device], 0 };

	// Create context from the chosen platform
	if (kerroin == 0.f) {
		// If a single device was selected (options.cpu_to_gpu_factor = 0), use GPU if possible
		context = clCreateContextFromType(properties, CL_DEVICE_TYPE_GPU, NULL, NULL, &status);
		if (status != CL_SUCCESS) {
			// Otherwise CPU
			context = clCreateContextFromType(properties, CL_DEVICE_TYPE_CPU, NULL, NULL, &status);
		}
	}
	// Use all devices if options.cpu_to_gpu_factor > 0
	else
		context = clCreateContextFromType(properties, CL_DEVICE_TYPE_ALL, NULL, NULL, &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	// Get the size of the device ID variable
	status = clGetContextInfo(context, CL_CONTEXT_DEVICES, 0, NULL, &size);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}
	// Get device IDs
	cl_device_id * devices2 = (cl_device_id*)alloca(size);
	status = clGetContextInfo(context, CL_CONTEXT_DEVICES, size, devices2, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}

	// Ignore devices with less than 2GB of memory
	// Also check global memory size and the feasibility of storing the normalization constants
	cl_int n_ignores = 0;
	cl_int n_gpus = 0;
	std::vector<cl_int> ignores(size / sizeof(cl_device_id));
	if (size / sizeof(cl_device_id) > 1ULL) {
		cl_ulong mem_max = 0ULL;
		for (size_t i = 0ULL; i < size / sizeof(cl_device_id); i++)
		{
			cl_device_type type;
			cl_ulong mem;
			clGetDeviceInfo(devices2[i], CL_DEVICE_TYPE, sizeof(cl_device_type), &type, NULL);
			switch (type)
			{
			case CL_DEVICE_TYPE_GPU:
				status = clGetDeviceInfo(devices2[i], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &mem, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				if ((mem * 1024ULL * 1024ULL) < 2000ULL || (kerroin == 0.f && mem < mem_max)) {
					ignores[i] = 1;
					n_ignores++;
				}
				else {
					ignores[i] = 0;
					mem_max = mem;
					if ((static_cast<cl_float>(mem) * mem_portions) > image_bytes)
						compute_norm_matrix = 0u;
					else
						compute_norm_matrix = 1u;
				}
				n_gpus++;
				break;
			case CL_DEVICE_TYPE_CPU:
				status = clGetDeviceInfo(devices2[i], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &mem, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					return status;
				}
				cpu_device = static_cast<int>(i);
				if ((static_cast<cl_float>(mem) * mem_portions) > image_bytes)
					compute_norm_matrix = 0u;
				else
					compute_norm_matrix = 1u;
				break;
			case CL_DEVICE_TYPE_ACCELERATOR:
				break;
			}
		}
	}
	else {
		cl_ulong mem;
		status = clGetDeviceInfo(devices2[0], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &mem, NULL);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		if ((static_cast<cl_float>(mem) * mem_portions) > image_bytes)
			compute_norm_matrix = 0u;
		else
			compute_norm_matrix = 1u;
	}

	delete[] platforms;

	// Get the number of devices
	status = clGetContextInfo(context, CL_CONTEXT_NUM_DEVICES, sizeof(cl_uint), &num_devices_context, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}

	// Remove devices that were ignored above
	uint32_t ll = 0u;
	for (size_t i = 0; i < num_devices_context; i++) {
		if (ignores[i] == 1) {
			num_devices_context--;
		}
		else {
			devices[ll] = devices2[i];
		}
		ll++;
	}

	return status;
}

// Get context for a single device
cl_int clGetPlatformsContextSingle(const uint32_t device, cl_context& context, cl_uint& num_devices_context, cl_device_id* devices) {
	cl_int status = CL_SUCCESS;
	cl_uint num_platforms = 0;

	cl_uint num_devices = 0;

	status = clGetPlatformIDs(0, NULL, &num_platforms);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		return status;
	}

	if (device >= num_platforms || device < 0) {
		std::cerr << "The specified platform number is greater or smaller than the available platform numbers!" << std::endl;
		return -1;
	}

	cl_platform_id *platforms = new cl_platform_id[num_platforms];

	status = clGetPlatformIDs(num_platforms, platforms, NULL);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		delete[] platforms;
		return status;
	}

	clGetDeviceIDs(platforms[device], CL_DEVICE_TYPE_ALL, 0, NULL, &num_devices);

	cl_device_id * devices2 = (cl_device_id*)alloca(num_devices * sizeof(cl_device_id));

	clGetDeviceIDs(platforms[device], CL_DEVICE_TYPE_ALL, num_devices, devices2, NULL);

	// Choose the GPU with highest amount of memory
	// If no GPU, use CPU
	if (num_devices > 1) {
		cl_ulong mem_prev = 0;
		for (size_t i = 0; i < num_devices; i++) {
			cl_device_type type;
			clGetDeviceInfo(devices2[i], CL_DEVICE_TYPE, sizeof(cl_device_type), &type, NULL);
			switch (type) {
			case CL_DEVICE_TYPE_GPU:
				cl_ulong mem;
				status = clGetDeviceInfo(devices2[i], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &mem, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					delete[] platforms;
					return status;
				}
				if (mem_prev < mem) {
					devices[0] = devices2[i];
					mem_prev = mem;
				}
				break;
			case CL_DEVICE_TYPE_CPU:
				if (mem_prev == 0)
					devices[0] = devices2[i];
				break;
			case CL_DEVICE_TYPE_ACCELERATOR:
				break;
			}
		}
	}
	else
		devices[0] = devices2[0];


	cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platforms[device], 0 };

	context = clCreateContext(properties, 1, devices, NULL, NULL, &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		delete[] platforms;
		return status;
	}

	delete[] platforms;

	num_devices_context = 1;

	return CL_SUCCESS;
}

// Build the programs and get the command queues
cl_int ClBuildProgramGetQueues(cl_program& program, const char* k_path, const cl_context context, const cl_uint num_devices_context, const cl_device_id* devices,
	const bool verbose, cl_command_queue* commandQueues, bool& atomic_64bit, const uint32_t projector_type, const char* header_directory) {
	cl_int status = CL_SUCCESS;


	cl_device_id* devices2 = (cl_device_id*)alloca(num_devices_context * sizeof(cl_device_id));
	for (size_t i = 0; i < num_devices_context; i++) {
		devices2[i] = devices[i];
	}

	// If integer atomic 64-bit operations are enabled, check if they are supported by the device(s)
	if (atomic_64bit) {
		std::string kernel_path_atom;

		kernel_path_atom = k_path;
		kernel_path_atom += "_64atom.cl";
		// Load the source text file
		std::fstream sourceFile_atom(kernel_path_atom.c_str());
		std::string content_atom((std::istreambuf_iterator<char>(sourceFile_atom)), std::istreambuf_iterator<char>());
		const char* sourceCode_atom = new char[content_atom.size()];
		sourceCode_atom = content_atom.c_str();
		// Create the program from the source
		program = clCreateProgramWithSource(context, 1, (const char**)& sourceCode_atom, NULL, &status);

		// Build the program
		status = clBuildProgram(program, num_devices_context, devices2, header_directory, NULL, NULL);
		//if (status != CL_SUCCESS) {
			//std::cerr << getErrorString(status) << std::endl;
			//mexPrintf("Failed to build OpenCL program. Build log: \n");
			//size_t len;
			//char* buffer;
			//clGetProgramBuildInfo(program, devices2[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
			//buffer = (char*)calloc(len, sizeof(size_t));
			//clGetProgramBuildInfo(program, devices2[0], CL_PROGRAM_BUILD_LOG, len, buffer, NULL);
			//mexPrintf("%s\n", buffer);
			//return status;
		//}
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
		program = clCreateProgramWithSource(context, 1, (const char**)& sourceCode, NULL, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
		status = clBuildProgram(program, num_devices_context, devices2, header_directory, NULL, NULL);
		// Build log in case of failure
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			mexPrintf("Failed to build OpenCL program. Build log: \n");
			size_t len;
			char* buffer;
			clGetProgramBuildInfo(program, devices2[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
			buffer = (char*)calloc(len, sizeof(size_t));
			clGetProgramBuildInfo(program, devices2[0], CL_PROGRAM_BUILD_LOG, len, buffer, NULL);
			mexPrintf("%s\n", buffer);
			return status;
		}
	}
	if (verbose)
		mexPrintf("OpenCL program built\n");

	// Create the command queues
	// Enable out of order execution (devices can compute kernels at the same time)
	for (size_t i = 0; i < num_devices_context; i++) {
		commandQueues[i] = clCreateCommandQueue(context, devices[i], CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &status);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			return status;
		}
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

	return status;
}