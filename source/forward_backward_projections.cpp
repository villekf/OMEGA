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
void f_b_project(const cl_uint& num_devices_context, const float kerroin, const int cpu_device, const cl::Context& context, const std::vector<cl::CommandQueue>& commandQueues,
	const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, const float* rhs, const mxArray* sc_ra, const uint32_t Nx,
	const uint32_t Ny, const uint32_t Nz, const float dx, const float dy, const float dz, const float bx, const float by, const float bz, const float bzb,
	const float maxxx, const float maxyy, const float zmax, const float NSlices, const uint32_t* pituus, const size_t koko_l, const uint32_t* xy_index,
	const uint16_t* z_index, const uint32_t size_x, const uint32_t TotSinos, const bool verbose, const uint32_t randoms_correction,
	const uint32_t attenuation_correction, const uint32_t normalization, const float* atten, const size_t size_atten, const float* norm, const size_t size_norm,
	const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L, const uint8_t raw, const size_t size_z, const uint32_t im_dim,
	const cl::Kernel& kernel_sum, const cl::Kernel& kernel, mxArray* output, const size_t size_rhs, const bool no_norm, const size_t numel_x,
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
	std::vector<size_t> length(num_devices_context);

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


	std::vector<cl::Buffer> d0_output, d0_Summ;
	std::vector<cl::Buffer> d_z(num_devices_context);
	std::vector<cl::Buffer> d_x(num_devices_context);
	std::vector<cl::Buffer> d_y(num_devices_context);
	std::vector<cl::Buffer> d_xcenter(num_devices_context);
	std::vector<cl::Buffer> d_ycenter(num_devices_context);
	std::vector<cl::Buffer> d_zcenter(num_devices_context);
	std::vector<cl::Buffer> d_V(num_devices_context);
	std::vector<cl::Buffer> d_atten(num_devices_context);
	std::vector<cl::Buffer> d_pseudos(num_devices_context);
	std::vector<cl::Buffer> d_rhs(num_devices_context);
	std::vector<cl::Buffer> d_output(num_devices_context);
	std::vector<cl::Buffer> d_Summ(num_devices_context);
	std::vector<cl::Buffer> d_Sino(num_devices_context);
	std::vector<cl::Buffer> d_sc_ra(num_devices_context);
	std::vector<cl::Buffer> d_norm(num_devices_context);
	std::vector<cl::Buffer> d_scat(num_devices_context);
	std::vector<cl::Buffer> d_lor(num_devices_context);
	std::vector<cl::Buffer> d_xyindex(num_devices_context);
	std::vector<cl::Buffer> d_zindex(num_devices_context);
	std::vector<cl::Buffer> d_L(num_devices_context);
	std::vector<cl::Buffer> d_reko_type(num_devices_context);
	if (num_devices_context > 1u) {
		d0_output.resize(num_devices_context - 1u);
		d0_Summ.resize(num_devices_context - 1u);
	}

	// Create the necessary buffers
	for (cl_uint i = 0U; i < num_devices_context; i++) {
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
		d_pseudos[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t) * prows, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		d_rhs[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * size_rhs, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		if (atomic_64bit) {
			d_output[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(cl_ulong) * size_output, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			d_Summ[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(cl_ulong) * im_dim, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (i < num_devices_context - 1) {
				d0_output[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(cl_ulong) * size_output, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d0_Summ[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(cl_ulong) * im_dim, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
		}
		else {
			d_output[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(cl_float) * size_output, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			d_Summ[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(cl_float) * im_dim, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (i < num_devices_context - 1) {
				d0_output[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(cl_float) * size_output, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d0_Summ[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(cl_float) * im_dim, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
		}
		if (randoms_correction == 1u) {
			d_sc_ra[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * length[i], NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else {
			d_sc_ra[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		if (normalization == 1u) {
			d_norm[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(cl_float) * length[i], NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else {
			d_norm[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(cl_float), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		if (scatter == 1u) {
			d_scat[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(cl_float) * length[i], NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else {
			d_scat[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(cl_float), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		d_Sino[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * length[i], NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		if (precompute)
			d_lor[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[i], NULL, &status);
		else
			d_lor[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[i], NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		if (raw) {
			d_xyindex[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t), NULL, &status);
			d_zindex[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			d_L[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[i] * 2, NULL, &status);
		}
		else {
			d_xyindex[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t) * length[i], NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			d_zindex[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[i], NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			d_L[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
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

	cl::Kernel kernel_ = kernel;
	cl::Kernel kernel_sum_ = kernel_sum;
	const float* Sino = (float*)mxGetData(mxGetCell(Sin, static_cast<mwIndex>(0)));

	for (cl_uint i = 0; i < num_devices_context; i++) {
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
		status = commandQueues[i].enqueueWriteBuffer(d_pseudos[i], CL_FALSE, 0, sizeof(uint32_t) * prows, pseudos);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = commandQueues[i].enqueueWriteBuffer(d_rhs[i], CL_FALSE, 0, sizeof(float) * size_rhs, rhs);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = commandQueues[i].enqueueWriteBuffer(d_Sino[i], CL_FALSE, 0, sizeof(float) * length[i], &Sino[cumsum[i]]);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}


		if (raw) {
			status = commandQueues[i].enqueueWriteBuffer(d_xyindex[i], CL_FALSE, 0, sizeof(uint32_t), xy_index);
			status = commandQueues[i].enqueueWriteBuffer(d_zindex[i], CL_FALSE, 0, sizeof(uint16_t), z_index);
			status = commandQueues[i].enqueueWriteBuffer(d_L[i], CL_FALSE, 0, sizeof(uint16_t) * length[i] * 2, &L[cumsum[i] * 2]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else {
			status = commandQueues[i].enqueueWriteBuffer(d_xyindex[i], CL_FALSE, 0, sizeof(uint32_t) * length[i], &xy_index[cumsum[i]]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			status = commandQueues[i].enqueueWriteBuffer(d_zindex[i], CL_FALSE, 0, sizeof(uint16_t) * length[i], &z_index[cumsum[i]]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			status = commandQueues[i].enqueueWriteBuffer(d_L[i], CL_FALSE, 0, sizeof(uint16_t), L);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		if (precompute)
			status = commandQueues[i].enqueueWriteBuffer(d_lor[i], CL_FALSE, 0, sizeof(uint16_t) * length[i], &lor1[cumsum[i]]);
		else
			status = commandQueues[i].enqueueWriteBuffer(d_lor[i], CL_FALSE, 0, sizeof(uint16_t), lor1);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		if (normalization == 1u) {
			status = commandQueues[i].enqueueWriteBuffer(d_norm[i], CL_FALSE, 0, sizeof(cl_float) * length[i], norm);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else {
			status = commandQueues[i].enqueueWriteBuffer(d_norm[i], CL_FALSE, 0, sizeof(cl_float), norm);
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
	}

	uint64_t size_rhs_i = static_cast<uint64_t>(size_rhs);


	std::vector<std::vector<float>> testi_summ;
	std::vector<std::vector<float>> testi_rhs;
	std::vector<std::vector<uint64_t>> testi_summ_u;
	std::vector<std::vector<uint64_t>> testi_rhs_u;
	if (num_devices_context > 1u) {
		if (atomic_64bit) {
			testi_summ_u.resize(num_devices_context - 1u, std::vector<uint64_t>(im_dim));
			testi_rhs_u.resize(num_devices_context - 1u, std::vector<uint64_t>(im_dim));
		}
		else {
			testi_summ.resize(num_devices_context - 1u, std::vector<float>(im_dim));
			testi_rhs.resize(num_devices_context - 1u, std::vector<float>(im_dim));
		}
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		status = commandQueues[i].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}

	cl_uint kernelInd = 0U;

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
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		status = commandQueues[i].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}

	cl_uint event_count = (1u);


	//cl_event** events = (cl_event * *)malloc(num_devices_context * sizeof(cl_event * *));
	//for (cl_uint i = 0u; i < num_devices_context; i++) {
	//	events[i] = (cl_event*)malloc(event_count * sizeof(cl_event * *));
	//}
	std::vector<std::vector<cl::Event>> events(num_devices_context, std::vector<cl::Event>(1));

	for (cl_uint i = 0; i < num_devices_context; i++) {
		if (atomic_64bit) {
			status = commandQueues[i].enqueueFillBuffer(d_output[i], zeroULL, 0, sizeof(cl_ulong) * size_output);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			status = commandQueues[i].enqueueFillBuffer(d_Summ[i], zeroULL, 0, sizeof(cl_ulong) * im_dim);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else {
			status = commandQueues[i].enqueueFillBuffer(d_Summ[i], zero, 0, sizeof(cl_float) * im_dim);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			status = commandQueues[i].enqueueFillBuffer(d_output[i], zero, 0, sizeof(cl_float) * size_output);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		if (randoms_correction == 1u) {
			float* S_R = (float*)mxGetData(mxGetCell(sc_ra, 0));
			status = commandQueues[i].enqueueWriteBuffer(d_sc_ra[i], CL_FALSE, 0, sizeof(float) * length[i], S_R);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else {
			status = commandQueues[i].enqueueFillBuffer(d_sc_ra[i], zero, 0, sizeof(cl_float));
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		if (scatter == 1u) {
			float* scat = (float*)mxGetData(mxGetCell(mxGetField(options, 0, "ScatterFB"), 0));
			status = commandQueues[i].enqueueWriteBuffer(d_scat[i], CL_FALSE, 0, sizeof(float) * length[i], scat);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else {
			status = commandQueues[i].enqueueFillBuffer(d_scat[i], zero, 0, sizeof(cl_float));
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
	}

	for (cl_uint i = 0; i < num_devices_context; i++) {
		commandQueues[i].finish();
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
		cl::NDRange global(global_size);
		cl::NDRange local(local_size);

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
		kernel_.setArg(kernelIndSubIter++, d_norm[i]);
		kernel_.setArg(kernelIndSubIter++, d_scat[i]);
		kernel_.setArg(kernelIndSubIter++, d_Summ[i]);
		kernel_.setArg(kernelIndSubIter++, d_lor[i]);
		kernel_.setArg(kernelIndSubIter++, d_xyindex[i]);
		kernel_.setArg(kernelIndSubIter++, d_zindex[i]);
		kernel_.setArg(kernelIndSubIter++, d_L[i]);
		kernel_.setArg(kernelIndSubIter++, d_Sino[i]);
		kernel_.setArg(kernelIndSubIter++, d_sc_ra[i]);
		kernel_.setArg(kernelIndSubIter++, d_rhs[i]);
		kernel_.setArg(kernelIndSubIter++, d_output[i]);
		kernel_.setArg(kernelIndSubIter++, no_norm);
		kernel_.setArg(kernelIndSubIter++, m_size);
		kernel_.setArg(kernelIndSubIter++, st);
		status = commandQueues[i].enqueueNDRangeKernel(kernel_, cl::NullRange, global, local, NULL, &events[i][0]);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}

	if (num_devices_context > 1u) {
		for (cl_uint i = 1; i < num_devices_context; i++) {
			if (atomic_64bit) {
				status = commandQueues[i].enqueueReadBuffer(d_Summ[i], CL_TRUE, 0, sizeof(uint64_t) * im_dim, testi_summ_u[i - 1u].data(), &events[i]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueReadBuffer(d_output[i], CL_TRUE, 0, sizeof(uint64_t) * size_output, testi_rhs_u[i - 1u].data(), &events[i]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d0_Summ[i - 1u], CL_FALSE, 0, sizeof(uint64_t) * im_dim, testi_summ_u[i - 1u].data());
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d0_output[i - 1u], CL_FALSE, 0, sizeof(uint64_t) * size_output, testi_rhs_u[i - 1u].data());
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			else {
				status = commandQueues[i].enqueueReadBuffer(d_Summ[i], CL_TRUE, 0, sizeof(float) * im_dim, testi_summ[i - 1u].data(), &events[i]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueReadBuffer(d_output[i], CL_TRUE, 0, sizeof(float) * size_output, testi_rhs[i - 1u].data(), &events[i]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d0_Summ[i - 1u], CL_FALSE, 0, sizeof(float) * im_dim, testi_summ_u[i - 1u].data());
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d0_output[i - 1u], CL_FALSE, 0, sizeof(float) * size_output, testi_rhs_u[i - 1u].data());
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			for (cl_uint i = 0ULL; i < num_devices_context; i++) {
				commandQueues[i].finish();
			}
		}
	}

	events[0][0].waitForEvents(events[0]);

	for (cl_uint i = 1; i < num_devices_context; i++) {
		size_t summa_pituus;
		if (fp == 1)
			summa_pituus = length[i];
		else
			summa_pituus = size_rhs;
		cl::NDRange global(summa_pituus);
		kernel_sum_.setArg(0, d0_Summ[i - 1]);
		kernel_sum_.setArg(1, d_Summ[0]);
		kernel_sum_.setArg(2, d0_output[i - 1]);
		kernel_sum_.setArg(3, d_output[0]);
		kernel_sum_.setArg(4, size_rhs_i);
		kernel_sum_.setArg(5, im_dim);
		kernel_sum_.setArg(6, no_norm);
		kernel_sum_.setArg(7, globals[i - 1ULL]);
		kernel_sum_.setArg(8, fp);
		status = commandQueues[0].enqueueNDRangeKernel(kernel_, cl::NullRange, global);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		commandQueues[i].finish();
	}

	if (atomic_64bit) {
		status = commandQueues[0].enqueueReadBuffer(d_output[0], CL_FALSE, 0, sizeof(cl_ulong) * size_output, output_u);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = commandQueues[0].enqueueReadBuffer(d_Summ[0], CL_FALSE, 0, sizeof(cl_ulong) * im_dim, normalizer_u);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}
	else {
		status = commandQueues[0].enqueueReadBuffer(d_output[0], CL_FALSE, 0, sizeof(float) * size_output, output_f);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = commandQueues[0].enqueueReadBuffer(d_Summ[0], CL_FALSE, 0, sizeof(float) * im_dim, normalizer_f);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}

	for (cl_uint i = 0u; i < num_devices_context; i++) {
		commandQueues[i].finish();
	}

	mxSetCell(output, 0, output_m);
	mxSetCell(output, 1, normalizer_m);
	return;
}