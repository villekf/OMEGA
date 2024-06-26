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
	const float maxxx, const float maxyy, const float zmax, const float NSlices, const int64_t* pituus, const size_t koko_l, const uint32_t* xy_index,
	const uint16_t* z_index, const uint32_t size_x, const uint32_t TotSinos, const bool verbose, const uint32_t randoms_correction,
	const uint32_t attenuation_correction, const uint32_t normalization, const float* atten, const size_t size_atten, const float* norm, const size_t size_norm,
	const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L, const uint8_t raw, const size_t size_z, const uint32_t im_dim,
	const cl::Kernel& kernel_sum, const cl::Kernel& kernel, mxArray* output, const size_t size_rhs, const cl_uchar no_norm, const size_t numel_x,
	const float tube_width, const float crystal_size_z, const float* x_center, const float* y_center, const float* z_center, const size_t size_center_x,
	const size_t size_center_y, const size_t size_center_z, const bool precompute, const int32_t dec, const uint32_t projector_type, const uint16_t n_rays, 
	const uint16_t n_rays3D, const float cr_pz, const mxArray* Sin, const bool atomic_64bit, const bool atomic_32bit, const float global_factor, const float bmin, const float bmax,
	const float Vmax, const float* V, const size_t size_V, const uint8_t fp, const size_t local_size, const mxArray* options, const uint32_t scatter, const bool TOF,
	const int64_t TOFSize, const float sigma_x, const float* TOFCenter, const int64_t nBins) {

	const uint32_t Nxy = Nx * Ny;
	cl_int status = CL_SUCCESS;
	cl_float zero = 0.f;
	cl_short zeroL = 0;
	cl_ulong zeroULL = 0ULL;
	cl_uint zero32 = 0;
	const float epps = 1e-8f;

	// Distance between rays in multi-ray Siddon
	const float dc_z = cr_pz / static_cast<float>(n_rays3D);
	const bool CT = (bool)mxGetScalar(mxGetField(options, 0, "CT"));

	int64_t nProjections = 0LL;
	float dPitch = 0.f;
	float* angles = nullptr;
	uint32_t size_y = 0U;
	uint32_t subsets = 0U;
	if (CT) {
		dPitch = (float)mxGetScalar(mxGetField(options, 0, "dPitch"));
		nProjections = (int64_t)mxGetScalar(mxGetField(options, 0, "nProjections"));
		size_y = (uint32_t)mxGetScalar(mxGetField(options, 0, "xSize"));
#if MX_HAS_INTERLEAVED_COMPLEX
		angles = (float*)mxGetSingles(mxGetField(options, 0, "angles"));
#else
		angles = (float*)mxGetData(mxGetField(options, 0, "angles"));
#endif
		if ((bool)mxGetScalar(mxGetField(options, 0, "useSubsets")))
			subsets = 2U;
	}

	size_t size_output;

	if (fp == 1)
		size_output = pituus[0] * nBins;
	else
		size_output = static_cast<size_t>(im_dim);

	const mwSize dimmi[1] = { static_cast<mwSize>(size_output) };

	mxArray* output_m, * normalizer_m;
	float* output_f, * normalizer_f;
	int64_t* output_64, * normalizer_64;
	int32_t* output_32, * normalizer_32;
	if (atomic_64bit) {
		// Output matrix
		output_m = mxCreateNumericMatrix(size_output, 1, mxINT64_CLASS, mxREAL);
#if MX_HAS_INTERLEAVED_COMPLEX
		output_64 = (int64_t*)mxGetInt64s(output_m);
#else
		output_64 = (int64_t*)mxGetData(output_m);
#endif
		normalizer_m = mxCreateNumericMatrix(im_dim, 1, mxINT64_CLASS, mxREAL);
#if MX_HAS_INTERLEAVED_COMPLEX
		normalizer_64 = (int64_t*)mxGetInt64s(normalizer_m);
#else
		normalizer_64 = (int64_t*)mxGetData(normalizer_m);
#endif
	}
	else if (atomic_32bit) {
		// Output matrix
		output_m = mxCreateNumericMatrix(size_output, 1, mxINT32_CLASS, mxREAL);
#if MX_HAS_INTERLEAVED_COMPLEX
		output_32 = (int32_t*)mxGetInt64s(output_m);
#else
		output_32 = (int32_t*)mxGetData(output_m);
#endif
		normalizer_m = mxCreateNumericMatrix(im_dim, 1, mxINT32_CLASS, mxREAL);
#if MX_HAS_INTERLEAVED_COMPLEX
		normalizer_32 = (int32_t*)mxGetInt64s(normalizer_m);
#else
		normalizer_32 = (int32_t*)mxGetData(normalizer_m);
#endif
	}
	else {
		output_m = mxCreateNumericMatrix(size_output, 1, mxSINGLE_CLASS, mxREAL);
#if MX_HAS_INTERLEAVED_COMPLEX
		output_f = (float*)mxGetSingles(output_m);
#else
		output_f = (float*)mxGetData(output_m);
#endif
		normalizer_m = mxCreateNumericMatrix(im_dim, 1, mxSINGLE_CLASS, mxREAL);
#if MX_HAS_INTERLEAVED_COMPLEX
		normalizer_f = (float*)mxGetSingles(normalizer_m);
#else
		normalizer_f = (float*)mxGetData(normalizer_m);
#endif
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

	const bool listmode = (bool)mxGetScalar(mxGetField(options, 0, "listmode"));
	size_t size_scat = 1ULL;
	if (scatter == 1U) {
		size_scat = mxGetNumberOfElements(mxGetCell(mxGetField(options, 0, "ScatterFB"), 0));
	}


	std::vector<cl::Buffer> d0_output, d0_Summ;
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
		d_reko_type[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(cl_uchar), NULL, &status);
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
		d_TOFCenter[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * nBins, NULL, &status);
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
		else if (atomic_32bit) {
			d_output[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(cl_int) * size_output, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			d_Summ[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(cl_uint) * im_dim, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (i < num_devices_context - 1) {
				d0_output[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(cl_int) * size_output, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				d0_Summ[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(cl_int) * im_dim, NULL, &status);
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
		if (listmode != 2)
			d_Sino[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * length[i] * nBins, NULL, &status);
		else
			d_Sino[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
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
		if (raw && listmode != 1) {
			d_xyindex[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t), NULL, &status);
			d_zindex[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			d_L[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[i] * 2, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else if (listmode != 1 && (!CT || subsets > 1)) {
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
		else {
			d_xyindex[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t), NULL, &status);
			d_zindex[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			d_L[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
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


	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}

	cl::Kernel kernel_ = kernel;
	cl::Kernel kernel_sum_ = kernel_sum;
#if MX_HAS_INTERLEAVED_COMPLEX
	const float* Sino = (float*)mxGetSingles(mxGetCell(Sin, static_cast<mwIndex>(0)));
#else
	const float* Sino = (float*)mxGetData(mxGetCell(Sin, static_cast<mwIndex>(0)));
#endif

	for (cl_uint i = 0; i < num_devices_context; i++) {
		status = commandQueues[i].enqueueWriteBuffer(d_reko_type[i], CL_FALSE, 0, sizeof(cl_uchar), &fp);
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
		status = commandQueues[i].enqueueWriteBuffer(d_TOFCenter[i], CL_FALSE, 0, sizeof(float) * nBins, TOFCenter);
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
		if (TOF && listmode != 2) {
			for (int64_t to = 0LL; to < nBins; to++) {
				status = commandQueues[i].enqueueWriteBuffer(d_Sino[i], CL_FALSE, sizeof(float) * length[i] * to, sizeof(float) * length[i], &Sino[cumsum[i] + koko * to]);
			}
		}
		else if (listmode != 2)
			status = commandQueues[i].enqueueWriteBuffer(d_Sino[i], CL_FALSE, 0, sizeof(float) * length[i], &Sino[cumsum[i]]);
		else
			status = commandQueues[i].enqueueFillBuffer(d_Sino[i], zero, 0, sizeof(cl_float));
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}


		if (raw && listmode != 1) {
			status = commandQueues[i].enqueueWriteBuffer(d_xyindex[i], CL_FALSE, 0, sizeof(uint32_t), xy_index);
			status = commandQueues[i].enqueueWriteBuffer(d_zindex[i], CL_FALSE, 0, sizeof(uint16_t), z_index);
			status = commandQueues[i].enqueueWriteBuffer(d_L[i], CL_FALSE, 0, sizeof(uint16_t) * length[i] * 2, &L[cumsum[i] * 2]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else if (listmode != 1 && (!CT || subsets > 1)){
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
		else {
			status = commandQueues[i].enqueueWriteBuffer(d_xyindex[i], CL_FALSE, 0, sizeof(uint32_t), xy_index);
			status = commandQueues[i].enqueueWriteBuffer(d_zindex[i], CL_FALSE, 0, sizeof(uint16_t), z_index);
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

	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Buffer write failed\n");
		return;
	}
	else if (DEBUG) {
		mexPrintf("Buffer write succeeded\n");
		mexEvalString("pause(.0001);");
	}


	std::vector<std::vector<float>> testi_summ;
	std::vector<std::vector<float>> testi_rhs;
	std::vector<std::vector<int64_t>> testi_summ_u;
	std::vector<std::vector<int64_t>> testi_rhs_u;
	std::vector<std::vector<int32_t>> testi_summ_32;
	std::vector<std::vector<int32_t>> testi_rhs_32;
	if (num_devices_context > 1u) {
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
		else if (atomic_32bit) {
			status = commandQueues[i].enqueueFillBuffer(d_output[i], zero32, 0, sizeof(cl_uint) * size_output);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			status = commandQueues[i].enqueueFillBuffer(d_Summ[i], zero32, 0, sizeof(cl_uint) * im_dim);
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
#if MX_HAS_INTERLEAVED_COMPLEX
			float* S_R = (float*)mxGetSingles(mxGetCell(sc_ra, static_cast<mwIndex>(0)));
#else
			float* S_R = (float*)mxGetData(mxGetCell(sc_ra, 0));
#endif
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
#if MX_HAS_INTERLEAVED_COMPLEX
			float* scat = (float*)mxGetSingles(mxGetCell(mxGetField(options, 0, "ScatterFB"), static_cast<mwIndex>(0)));
#else
			float* scat = (float*)mxGetData(mxGetCell(mxGetField(options, 0, "ScatterFB"), 0));
#endif
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
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
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

		if (DEBUG) {
			mexPrintf("global_size = %u\n", global_size);
			mexPrintf("m_size = %u\n", m_size);
			mexPrintf("size_output = %u\n", size_output);
			mexPrintf("size_rhs = %u\n", size_rhs);
			mexPrintf("num_devices_context = %u\n", num_devices_context);
			mexPrintf("st = %u\n", st);
			mexPrintf("length[0] * nBins = %u\n", length[0] * nBins);
			mexEvalString("pause(.0001);");
		}

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
		else if (DEBUG) {
			mexPrintf("Kernel launched successfully\n");
			mexEvalString("pause(.0001);");
		}
	}

	if (num_devices_context > 1u) {
		for (cl_uint i = 1; i < num_devices_context; i++) {
			if (atomic_64bit) {
				status = commandQueues[i].enqueueReadBuffer(d_Summ[i], CL_TRUE, 0, sizeof(int64_t) * im_dim, testi_summ_u[i - 1u].data(), &events[i]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueReadBuffer(d_output[i], CL_TRUE, 0, sizeof(int64_t) * size_output, testi_rhs_u[i - 1u].data(), &events[i]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d0_Summ[i - 1u], CL_FALSE, 0, sizeof(int64_t) * im_dim, testi_summ_u[i - 1u].data());
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d0_output[i - 1u], CL_FALSE, 0, sizeof(int64_t) * size_output, testi_rhs_u[i - 1u].data());
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			else if (atomic_32bit) {
				status = commandQueues[i].enqueueReadBuffer(d_Summ[i], CL_TRUE, 0, sizeof(int32_t) * im_dim, testi_summ_32[i - 1u].data(), &events[i]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueReadBuffer(d_output[i], CL_TRUE, 0, sizeof(int32_t) * size_output, testi_rhs_32[i - 1u].data(), &events[i]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d0_Summ[i - 1u], CL_FALSE, 0, sizeof(int32_t) * im_dim, testi_summ_32[i - 1u].data());
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = commandQueues[i].enqueueWriteBuffer(d0_output[i - 1u], CL_FALSE, 0, sizeof(int32_t) * size_output, testi_rhs_32[i - 1u].data());
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
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
		}
	}

	events[0][0].waitForEvents(events[0]);

	for (cl_uint i = 1; i < num_devices_context; i++) {
		size_t summa_pituus = size_output;
		//if (fp == 1)
			//summa_pituus = size_output;
		//else
			//summa_pituus = im_dim;
		cl::NDRange global(summa_pituus);
		kernel_sum_.setArg(0, d0_Summ[i - 1]);
		kernel_sum_.setArg(1, d_Summ[0]);
		kernel_sum_.setArg(2, d0_output[i - 1]);
		kernel_sum_.setArg(3, d_output[0]);
		kernel_sum_.setArg(4, summa_pituus);
		kernel_sum_.setArg(5, no_norm);
		//kernel_sum_.setArg(6, im_dim);
		//kernel_sum_.setArg(7, globals[i - 1ULL]);
		//kernel_sum_.setArg(8, fp);
		status = commandQueues[0].enqueueNDRangeKernel(kernel_sum_, cl::NullRange, global);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		else if (DEBUG) {
			mexPrintf("Merge kernel launched successfully\n");
			mexEvalString("pause(.0001);");
		}
		commandQueues[i].finish();
	}

	if (atomic_64bit) {
		status = commandQueues[0].enqueueReadBuffer(d_output[0], CL_FALSE, 0, sizeof(cl_ulong) * size_output, output_64);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = commandQueues[0].enqueueReadBuffer(d_Summ[0], CL_FALSE, 0, sizeof(cl_ulong) * im_dim, normalizer_64);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}
	else if (atomic_32bit) {
		status = commandQueues[0].enqueueReadBuffer(d_output[0], CL_FALSE, 0, sizeof(cl_uint) * size_output, output_32);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = commandQueues[0].enqueueReadBuffer(d_Summ[0], CL_FALSE, 0, sizeof(cl_uint) * im_dim, normalizer_32);
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
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}

	mxSetCell(output, 0, output_m);
	mxSetCell(output, 1, normalizer_m);
	return;
}