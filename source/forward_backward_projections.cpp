/**************************************************************************
* Forward/backward projections (reconstruction_f_b_proj).
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

// Forward/backward projections
void f_b_project(const cl_uint& num_devices_context, const float kerroin, const int cpu_device, const cl_context& context, const cl_command_queue* commandQueues,
	const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, const float* rhs, const mxArray* sc_ra, const uint32_t Nx,
	const uint32_t Ny, const uint32_t Nz, const float dx, const float dy, const float dz, const float bx, const float by, const float bz, const float bzb,
	const float maxxx, const float maxyy, const float zmax, const float NSlices, const uint32_t* pituus, const size_t koko_l, const uint32_t* xy_index,
	const uint16_t* z_index, const uint32_t size_x, const uint16_t TotSinos, const bool verbose, const uint32_t randoms_correction,
	const uint32_t attenuation_correction, const uint32_t normalization, const float* atten, const size_t size_atten, const float* norm, const size_t size_norm,
	const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L, const uint8_t raw, const size_t size_z, const uint32_t im_dim,
	const cl_kernel& kernel_sum, const cl_kernel& kernel, mxArray* output, const size_t size_rhs, const bool no_norm, const size_t numel_x,
	const float tube_width, const float crystal_size_z, const float* x_center, const float* y_center, const float* z_center, const size_t size_center_x,
	const size_t size_center_y, const size_t size_center_z, const bool precompute, const int32_t dec, const uint32_t projector_type, const uint16_t n_rays, 
	const uint16_t n_rays3D, const float cr_pz, const mxArray* Sin, const bool atomic_64bit, const float global_factor, const float bmin, const float bmax, 
	const float Vmax, const float* V, const size_t size_V, const uint8_t fp, const size_t local_size, const mxArray* options, const uint32_t scatter) {

	const uint32_t Nxy = Nx * Ny;
	cl_int status = CL_SUCCESS;
	cl_float zero = 0.f;
	cl_short zeroL = 0;
	cl_ulong zeroULL = 0ULL;
	const float epps = 1e-8f;

	// Distance between rays in multi-ray Siddon
	const float dc_z = cr_pz / static_cast<float>(n_rays3D);

	size_t size_output;

	if (fp == 1)
		size_output = pituus[0];
	else
		size_output = static_cast<size_t>(im_dim);

	const mwSize dimmi[1] = { static_cast<mwSize>(size_output) };

	mxArray* output_m, * normalizer_m;
	float* output_f, * normalizer_f;
	uint64_t* output_u, * normalizer_u;
	if (atomic_64bit) {
		// Output matrix
		output_m = mxCreateNumericMatrix(size_output, 1, mxUINT64_CLASS, mxREAL);
		output_u = (uint64_t*)mxGetData(output_m);
		normalizer_m = mxCreateNumericMatrix(im_dim, 1, mxUINT64_CLASS, mxREAL);
		normalizer_u = (uint64_t*)mxGetData(normalizer_m);
	}
	else {
		output_m = mxCreateNumericMatrix(size_output, 1, mxSINGLE_CLASS, mxREAL);
		output_f = (float*)mxGetData(output_m);
		normalizer_m = mxCreateNumericMatrix(im_dim, 1, mxSINGLE_CLASS, mxREAL);
		normalizer_f = (float*)mxGetData(normalizer_m);
	}

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
	size_t size_scat = 1ULL;
	if (scatter == 1U) {
		size_scat = mxGetNumberOfElements(mxGetCell(mxGetField(options, 0, "ScatterFB"), 0));
	}


	std::vector<cl_mem> d0_output, d0_Summ;
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
	std::vector<cl_mem> d_output(num_devices_context, 0);
	std::vector<cl_mem> d_Summ(num_devices_context, 0);
	std::vector<cl_mem> d_Sino(num_devices_context, 0);
	std::vector<cl_mem> d_sc_ra(num_devices_context, 0);
	std::vector<cl_mem> d_norm(num_devices_context, 0);
	std::vector<cl_mem> d_scat(num_devices_context, 0);
	std::vector<cl_mem> d_lor(num_devices_context, 0);
	std::vector<cl_mem> d_xyindex(num_devices_context, 0);
	std::vector<cl_mem> d_zindex(num_devices_context, 0);
	std::vector<cl_mem> d_L(num_devices_context, 0);
	std::vector<cl_mem> d_reko_type(num_devices_context, 0);
	if (num_devices_context > 1u) {
		d0_output.resize(num_devices_context - 1u);
		d0_Summ.resize(num_devices_context - 1u);
	}

	// Create the necessary buffers
	for (cl_uint i = 0U; i < num_devices_context; i++) {
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
		d_rhs[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * size_rhs, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		if (atomic_64bit) {
			d_output[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_ulong) * size_output, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			d_Summ[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_ulong) * im_dim, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (i < num_devices_context - 1) {
				d0_output[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_ulong) * size_output, NULL, &status);
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
			d_output[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * size_output, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			d_Summ[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * im_dim, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (i < num_devices_context - 1) {
				d0_output[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * size_output, NULL, &status);
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
		if (randoms_correction == 1u) {
			d_sc_ra[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * length[i], NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else {
			d_sc_ra[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		if (normalization == 1u) {
			d_norm[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_float) * length[i], NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else {
			d_norm[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_float), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		if (scatter == 1u) {
			d_scat[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_float) * length[i], NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else {
			d_scat[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_float), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		d_Sino[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * length[i], NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		if (precompute)
			d_lor[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[i], NULL, &status);
		else
			d_lor[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		if (raw) {
			d_xyindex[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t), NULL, &status);
			d_zindex[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			d_L[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[i] * 2, NULL, &status);
		}
		else {
			d_xyindex[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t) * length[i], NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			d_zindex[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[i], NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			d_L[i] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
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

	const float* Sino = (float*)mxGetData(mxGetCell(Sin, static_cast<mwIndex>(0)));

	for (cl_uint i = 0; i < num_devices_context; i++) {
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
		status = clEnqueueWriteBuffer(commandQueues[i], d_rhs[i], CL_FALSE, 0, sizeof(float) * size_rhs, rhs, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = clEnqueueWriteBuffer(commandQueues[i], d_Sino[i], CL_FALSE, 0, sizeof(float) * length[i], &Sino[cumsum[i]], 0, NULL, NULL);


		if (raw) {
			status = clEnqueueWriteBuffer(commandQueues[i], d_xyindex[i], CL_FALSE, 0, sizeof(uint32_t), xy_index, 0, NULL, NULL);
			status = clEnqueueWriteBuffer(commandQueues[i], d_zindex[i], CL_FALSE, 0, sizeof(uint16_t), z_index, 0, NULL, NULL);
			status = clEnqueueWriteBuffer(commandQueues[i], d_L[i], CL_FALSE, 0, sizeof(uint16_t) * length[i] * 2, &L[cumsum[i] * 2], 0, NULL, NULL);
		}
		else {
			status = clEnqueueWriteBuffer(commandQueues[i], d_xyindex[i], CL_FALSE, 0, sizeof(uint32_t) * length[i], &xy_index[cumsum[i]], 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			status = clEnqueueWriteBuffer(commandQueues[i], d_zindex[i], CL_FALSE, 0, sizeof(uint16_t) * length[i], &z_index[cumsum[i]], 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			status = clEnqueueWriteBuffer(commandQueues[i], d_L[i], CL_FALSE, 0, sizeof(uint16_t), L, 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		if (precompute)
			status = clEnqueueWriteBuffer(commandQueues[i], d_lor[i], CL_FALSE, 0, sizeof(uint16_t) * length[i], &lor1[cumsum[i]], 0, NULL, NULL);
		else
			status = clEnqueueWriteBuffer(commandQueues[i], d_lor[i], CL_FALSE, 0, sizeof(uint16_t), lor1, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		if (normalization == 1u) {
			status = clEnqueueWriteBuffer(commandQueues[i], d_norm[i], CL_FALSE, 0, sizeof(cl_float) * length[i], &norm[cumsum[i]], 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else {
			status = clEnqueueWriteBuffer(commandQueues[i], d_norm[i], CL_FALSE, 0, sizeof(cl_float), norm, 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}

		status = clFlush(commandQueues[i]);
	}

	uint64_t size_rhs_i = static_cast<uint64_t>(size_rhs);


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


	for (cl_uint i = 0; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

	cl_uint kernelInd = 0U;

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
	clSetKernelArg(kernel, kernelInd++, sizeof(uint8_t), &fp);
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
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

	cl_uint event_count = (1u);


	//cl_event** events = (cl_event * *)malloc(num_devices_context * sizeof(cl_event * *));
	//for (cl_uint i = 0u; i < num_devices_context; i++) {
	//	events[i] = (cl_event*)malloc(event_count * sizeof(cl_event * *));
	//}
	std::vector<cl_event> events(num_devices_context, NULL);

	for (cl_uint i = 0; i < num_devices_context; i++) {
		status = clEnqueueFillBuffer(commandQueues[i], d_Summ[i], &zero, sizeof(cl_float), 0, sizeof(cl_float) * im_dim, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		if (atomic_64bit) {
			status = clEnqueueFillBuffer(commandQueues[i], d_output[i], &zeroULL, sizeof(cl_ulong), 0, sizeof(cl_ulong) * size_output, 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else {
			status = clEnqueueFillBuffer(commandQueues[i], d_output[i], &zero, sizeof(cl_float), 0, sizeof(cl_float) * size_output, 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		if (randoms_correction == 1u) {
			float* S_R = (float*)mxGetData(mxGetCell(sc_ra, 0));
			status = clEnqueueWriteBuffer(commandQueues[i], d_sc_ra[i], CL_FALSE, 0, sizeof(float) * length[i], &S_R[cumsum[i]], 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else {
			status = clEnqueueFillBuffer(commandQueues[i], d_sc_ra[i], &zero, sizeof(cl_float), 0, sizeof(cl_float), 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		if (scatter == 1u) {
			float* scat = (float*)mxGetData(mxGetCell(mxGetField(options, 0, "ScatterFB"), 0));
			status = clEnqueueWriteBuffer(commandQueues[i], d_scat[i], CL_FALSE, 0, sizeof(float) * length[i], &scat[cumsum[i]], 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else {
			status = clEnqueueFillBuffer(commandQueues[i], d_scat[i], &zero, sizeof(cl_float), 0, sizeof(cl_float), 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
	}

	std::vector<cl_ulong> globals(num_devices_context, 0ULL);

	for (cl_uint i = 0; i < num_devices_context; i++) {

		const uint64_t st = static_cast<uint64_t>(cumsum[i]);

		cl_uint kernelIndSubIter = kernelInd;

		size_t erotus = length[i] % local_size;

		if (erotus > 0)
			erotus = (local_size - erotus);

		const size_t global_size = length[i] + erotus;
		const uint64_t m_size = static_cast<uint64_t>(length[i]);

		// These are needed to make sure the forward projection data is correctly combined
		if (i > 0)
			globals[i] = globals[i - 1ULL] + m_size;
		else
			globals[i] = m_size;

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
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_norm[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_scat[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_Summ[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_lor[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_xyindex[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_zindex[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_L[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_Sino[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_sc_ra[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_rhs[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_mem), &d_output[i]);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(cl_uchar), &no_norm);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(uint64_t), &m_size);
		clSetKernelArg(kernel, kernelIndSubIter++, sizeof(uint64_t), &st);
		status = clEnqueueNDRangeKernel(commandQueues[i], kernel, 1, NULL, &global_size, &local_size, 0, NULL, &events[i]);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			for (cl_uint i = 0; i < num_devices_context; i++) {
				clReleaseMemObject(d_reko_type[i]);
				clReleaseMemObject(d_z[i]);
				clReleaseMemObject(d_x[i]);
				clReleaseMemObject(d_y[i]);
				clReleaseMemObject(d_xcenter[i]);
				clReleaseMemObject(d_ycenter[i]);
				clReleaseMemObject(d_zcenter[i]);
				clReleaseMemObject(d_V[i]);
				clReleaseMemObject(d_atten[i]);
				clReleaseMemObject(d_norm[i]);
				clReleaseMemObject(d_scat[i]);
				clReleaseMemObject(d_pseudos[i]);
				clReleaseMemObject(d_Summ[i]);
				clReleaseMemObject(d_rhs[i]);
				clReleaseMemObject(d_output[i]);
				clReleaseMemObject(d_xyindex[i]);
				clReleaseMemObject(d_zindex[i]);
				clReleaseMemObject(d_L[i]);
				clReleaseMemObject(d_sc_ra[i]);
				clReleaseMemObject(d_lor[i]);
				clReleaseMemObject(d_Sino[i]);
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
			if (atomic_64bit) {
				status = clEnqueueReadBuffer(commandQueues[i], d_Summ[i], CL_TRUE, 0, sizeof(uint64_t) * im_dim, testi_summ_u[i - 1].data(), 1, &events[i], NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueReadBuffer(commandQueues[i], d_output[i], CL_TRUE, 0, sizeof(uint64_t) * size_output, testi_rhs_u[i - 1].data(), 1, &events[i], NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[0], d0_Summ[i - 1], CL_FALSE, 0, sizeof(uint64_t) * im_dim, testi_summ_u[i - 1].data(), 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[0], d0_output[i - 1], CL_FALSE, 0, sizeof(uint64_t) * size_output, testi_rhs_u[i - 1].data(), 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			else {
				status = clEnqueueReadBuffer(commandQueues[i], d_Summ[i], CL_TRUE, 0, sizeof(float) * im_dim, testi_summ[i - 1].data(), 1, &events[i], NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueReadBuffer(commandQueues[i], d_output[i], CL_TRUE, 0, sizeof(float) * size_output, testi_rhs[i - 1].data(), 1, &events[i], NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[0], d0_Summ[i - 1], CL_FALSE, 0, sizeof(float) * im_dim, testi_summ[i - 1].data(), 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = clEnqueueWriteBuffer(commandQueues[0], d0_output[i - 1], CL_FALSE, 0, sizeof(float) * size_output, testi_rhs[i - 1].data(), 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			for (cl_uint i = 0ULL; i < num_devices_context; i++) {
				clFinish(commandQueues[i]);
			}
		}
	}

	clWaitForEvents(1, &events[0]);

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
		clSetKernelArg(kernel_sum, 6, sizeof(cl_uchar), &no_norm);
		clSetKernelArg(kernel_sum, 7, sizeof(uint64_t), &globals[i - 1ULL]);
		clSetKernelArg(kernel_sum, 8, sizeof(uint8_t), &fp);
		status = clEnqueueNDRangeKernel(commandQueues[0], kernel_sum, 1, NULL, &summa_pituus, NULL, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		clFinish(commandQueues[0]);
	}

	if (atomic_64bit) {
		status = clEnqueueReadBuffer(commandQueues[0], d_output[0], CL_FALSE, 0, sizeof(uint64_t) * size_output, output_u, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = clEnqueueReadBuffer(commandQueues[0], d_Summ[0], CL_FALSE, 0, sizeof(uint64_t) * im_dim, normalizer_u, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}
	else {
		status = clEnqueueReadBuffer(commandQueues[0], d_output[0], CL_FALSE, 0, sizeof(float) * size_output, output_f, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = clEnqueueReadBuffer(commandQueues[0], d_Summ[0], CL_FALSE, 0, sizeof(float) * im_dim, normalizer_f, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}

	for (cl_uint i = 0u; i < num_devices_context; i++) {
		clFinish(commandQueues[i]);
		clReleaseEvent(events[i]);
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		clReleaseMemObject(d_reko_type[i]);
		clReleaseMemObject(d_z[i]);
		clReleaseMemObject(d_x[i]);
		clReleaseMemObject(d_y[i]);
		clReleaseMemObject(d_xcenter[i]);
		clReleaseMemObject(d_ycenter[i]);
		clReleaseMemObject(d_zcenter[i]);
		clReleaseMemObject(d_V[i]);
		clReleaseMemObject(d_atten[i]);
		clReleaseMemObject(d_norm[i]);
		clReleaseMemObject(d_scat[i]);
		clReleaseMemObject(d_pseudos[i]);
		clReleaseMemObject(d_Summ[i]);
		clReleaseMemObject(d_rhs[i]);
		clReleaseMemObject(d_output[i]);
		clReleaseMemObject(d_xyindex[i]);
		clReleaseMemObject(d_zindex[i]);
		clReleaseMemObject(d_L[i]);
		clReleaseMemObject(d_sc_ra[i]);
		clReleaseMemObject(d_lor[i]);
		clReleaseMemObject(d_Sino[i]);
		if (i < num_devices_context - 1) {
			clReleaseMemObject(d0_Summ[i]);
			clReleaseMemObject(d0_output[i]);
		}
	}
	mxSetCell(output, 0, output_m);
	mxSetCell(output, 1, normalizer_m);
}