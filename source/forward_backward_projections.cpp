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
void f_b_project(const cl_uint& num_devices_context, const float kerroin, const int cpu_device, const cl::Context& context, 
	const std::vector<cl::CommandQueue>& commandQueues, const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, 
	const float* rhs, const mxArray* sc_ra, scalarStruct inputScalars, const int64_t* pituus, const size_t koko_l, const uint32_t* xy_index,
	const uint16_t* z_index, const uint32_t TotSinos, const bool verbose, const float* atten, const size_t size_atten, const float* norm, 
	const size_t size_norm,	const uint32_t* pseudos, const uint16_t* L, const uint32_t im_dim, const cl::Kernel& kernel_sum, const cl::Kernel& kernel, 
	mxArray* output, const size_t size_rhs, const cl_uchar no_norm, const size_t numel_x, const float* x_center, const float* y_center, 
	const float* z_center, const size_t size_center_x, const size_t size_center_y, const size_t size_center_z, const mxArray* Sin, 
	const bool atomic_64bit, const bool atomic_32bit, const float* V, const size_t size_V, const size_t local_size[], const mxArray* options, 
	const int64_t TOFSize, const float* TOFCenter) {

	// Number of pixels in one slice
	const uint32_t Nxy = inputScalars.Nx * inputScalars.Ny;
	// Initialize variables
	cl_int status = CL_SUCCESS;
	cl_float zero = 0.f;
	cl_short zeroL = 0;
	cl_ulong zeroULL = 0ULL;
	cl_uint zero32 = 0;
	const float epps = 1e-8f;
	int field_num;

	// Distance between rays in multi-ray Siddon
	inputScalars.dc_z = inputScalars.cr_pz / static_cast<float>(inputScalars.n_rays3D);
	const bool CT = (bool)mxGetScalar(getField(options, 0, "CT"));
	//const bool SPECT = (bool)mxGetScalar(mxGetField(options, 0, "SPECT"));
	// Is listmode data used
	const uint8_t listmode = (uint8_t)mxGetScalar(getField(options, 0, "listmode"));

	//int64_t nProjections = 0LL;
	//float dPitch = 0.f;
	float* angles = nullptr, * mask1 = nullptr, * mask2 = nullptr, * integralX = nullptr, * integralY = nullptr, * uv = nullptr, * meanVals = nullptr;
	uint8_t *mask = nullptr;
	uint64_t nAngles = 0;
	cl_float2 dPitch = { 0.f, 0.f };
	//size_t numelMean;
	if (DEBUG) {
		//mexPrintf("numelM = %u\n", numelA);
		mexPrintf("CT = %d\n", CT);
		mexEvalString("pause(.0001);");
	}
	//uint32_t size_y = 0U;
	//uint32_t subsets = 0U;
	if (inputScalars.maskFP && inputScalars.fp == 1)
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		mask = (uint8_t*)mxGetUint8s(getField(options, 0, "maskFP"));
#else
		mask = (uint8_t*)mxGetData(getField(options, 0, "maskFP"));
#endif
	if (CT || inputScalars.SPECT) {
		// Detector pitch
		dPitch.s[0] = (float)mxGetScalar(getField(options, 0, "dPitchX"));
		dPitch.s[1] = (float)mxGetScalar(getField(options, 0, "dPitchY"));
		//dPitch = { (float)mxGetScalar(mxGetField(options, 0, "dPitchX")), (float)mxGetScalar(mxGetField(options, 0, "dPitchY")) };
		// Number of projections
		inputScalars.nProjections = (int64_t)mxGetScalar(getField(options, 0, "nProjectionsS"));
		// Number of detector crystals
		inputScalars.size_y = (uint32_t)mxGetScalar(getField(options, 0, "xSize"));
//#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
//		angles = (float*)mxGetSingles(mxGetField(options, 0, "angles"));
//#else
//		angles = (float*)mxGetData(mxGetField(options, 0, "angles"));
//#endif
//		nAngles = mxGetNumberOfElements(mxGetField(options, 0, "angles"));
			//inputScalars.size_z = inputScalars.nProjections * 2LL;
			// Direction vectors of the detector panel
			// Different versions depending on the MEX-API version
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			uv = (float*)mxGetSingles(getField(options, 0, "uVS"));
#else
			uv = (float*)mxGetData(getField(options, 0, "uVS"));
#endif
		if (inputScalars.projector_type == 5 || inputScalars.projector_type == 4) {
			// (Optional) masks for both forward and backward projections
			if (inputScalars.maskBP && inputScalars.fp == 2) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
				mask = (uint8_t*)mxGetUint8s(getField(options, 0, "maskBP"));
#else
				mask = (uint8_t*)mxGetData(getField(options, 0, "maskBP"));
#endif
				const size_t nMask = mxGetNumberOfElements(getField(options, 0, "maskBP"));
				if (DEBUG) {
					mexPrintf("nMask = %d\n", nMask);
					mexEvalString("pause(.0001);");
				}
			}
		}
		if (inputScalars.projector_type == 5) {
			// Integral image for branchless distance-driven method
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			integralX = (float*)mxGetSingles(getField(options, 0, "integralX"));
			if (inputScalars.fp == 1)
				integralY = (float*)mxGetSingles(getField(options, 0, "integralY"));
			// Mean values to reduce the dynamic range
			if ((inputScalars.meanFP && inputScalars.fp == 1) || (inputScalars.meanBP && inputScalars.fp == 2))
				meanVals = (float*)mxGetSingles(getField(options, 0, "meanV"));
#else
			integralX = (float*)mxGetData(getField(options, 0, "integralX"));
			if (inputScalars.fp == 1)
				integralY = (float*)mxGetData(getField(options, 0, "integralY"));
			if ((inputScalars.meanFP && inputScalars.fp == 1) || (inputScalars.meanBP && inputScalars.fp == 2))
				meanVals = (float*)mxGetData(getField(options, 0, "meanV"));
#endif
			//numelMean = mxGetNumberOfElements(mxGetField(options, 0, "meanV")); 
		}
		//const size_t numelA = mxGetNumberOfElements(mxGetField(options, 0, "angles"));
		if (DEBUG) {
			//mexPrintf("numelM = %u\n", numelA);
			mexPrintf("nProjections = %d\n", inputScalars.nProjections);
			mexEvalString("pause(.0001);");
		}
//		if (inputScalars.SPECT) {
//#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
//			mask2 = (float*)mxGetSingles(mxGetField(options, 0, "mask1"));
//#else
//			mask2 = (float*)mxGetData(mxGetField(options, 0, "mask1"));
//#endif
//#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
//			mask1 = (float*)mxGetSingles(mxGetField(options, 0, "mask2"));
//#else
//			mask1 = (float*)mxGetData(mxGetField(options, 0, "mask2"));
//#endif
//			const size_t numelM = mxGetNumberOfElements(mxGetField(options, 0, "mask2")); 
//			if (DEBUG) {
//				mexPrintf("numelM = %u\n", numelM);
//				mexPrintf("d_Scale0 = %u\n", inputScalars.d_Scale.s[0]);
//				mexPrintf("d_Scale1 = %u\n", inputScalars.d_Scale.s[1]);
//				mexPrintf("d_Scale2 = %u\n", inputScalars.d_Scale.s[2]);
//				mexEvalString("pause(.0001);");
//			}
//		}
	}
	else if (inputScalars.PET) {
		inputScalars.nProjections = (int64_t)mxGetScalar(getField(options, 0, "nProjectionsS"));
		//inputScalars.size_z = inputScalars.nProjections * 2LL;
		inputScalars.size_y = (uint32_t)mxGetScalar(getField(options, 0, "Nang"));
		inputScalars.size_x = (uint32_t)mxGetScalar(getField(options, 0, "Ndist"));
		dPitch.s[0] = (float)mxGetScalar(getField(options, 0, "cr_p"));
		dPitch.s[1] = (float)mxGetScalar(getField(options, 0, "cr_pz"));
	}
	else {
		// Detector pitch
		dPitch.s[0] = (float)mxGetScalar(getField(options, 0, "cr_p"));
		dPitch.s[1] = (float)mxGetScalar(getField(options, 0, "cr_pz"));
	}

	size_t size_output;

	// Size of the output vector
	if (inputScalars.fp == 1)
		if ((CT || inputScalars.PET) && listmode == 0)
			size_output = inputScalars.size_x * inputScalars.size_y * inputScalars.nProjections;
		else
			size_output = pituus[0] * inputScalars.nBins;
	else
		size_output = static_cast<size_t>(im_dim);

	const mwSize dimmi[1] = { static_cast<mwSize>(size_output) };

	mxArray* output_m, * normalizer_m;
	float* output_f, * normalizer_f;
	int64_t* output_64, * normalizer_64;
	int32_t* output_32, * normalizer_32;
	if (atomic_64bit) {
		// Output matrix if 64-bit atomics are used
		output_m = mxCreateNumericMatrix(size_output, 1, mxINT64_CLASS, mxREAL);
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		output_64 = (int64_t*)mxGetInt64s(output_m);
#else
		output_64 = (int64_t*)mxGetData(output_m);
#endif
		normalizer_m = mxCreateNumericMatrix(im_dim, 1, mxINT64_CLASS, mxREAL);
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		normalizer_64 = (int64_t*)mxGetInt64s(normalizer_m);
#else
		normalizer_64 = (int64_t*)mxGetData(normalizer_m);
#endif
	}
	else if (atomic_32bit) {
		// Output matrix if 32-bit atomics are used
		output_m = mxCreateNumericMatrix(size_output, 1, mxINT32_CLASS, mxREAL);
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		output_32 = (int32_t*)mxGetInt64s(output_m);
#else
		output_32 = (int32_t*)mxGetData(output_m);
#endif
		normalizer_m = mxCreateNumericMatrix(im_dim, 1, mxINT32_CLASS, mxREAL);
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		normalizer_32 = (int32_t*)mxGetInt64s(normalizer_m);
#else
		normalizer_32 = (int32_t*)mxGetData(normalizer_m);
#endif
	}
	else {
		// Output matrix if float atomics are used
		output_m = mxCreateNumericMatrix(size_output, 1, mxSINGLE_CLASS, mxREAL);
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
		output_f = (float*)mxGetSingles(output_m);
#else
		output_f = (float*)mxGetData(output_m);
#endif
		if (no_norm == 0)
			normalizer_m = mxCreateNumericMatrix(im_dim, 1, mxSINGLE_CLASS, mxREAL);
		else
			normalizer_m = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
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

	// Divide measurements between the (optional) devices
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
	float kerroin4 = 0.f;
	// Weighting coefficient for projector type 4
	if (inputScalars.projector_type == 4U)
		kerroin4 = (float)mxGetScalar(getField(options, 0, "kerroin"));
	size_t size_scat = 1ULL;
	if (inputScalars.scatter == 1U) {
		size_scat = mxGetNumberOfElements(mxGetCell(getField(options, 0, "ScatterFB"), 0));
	}


	std::vector<cl::Buffer> d0_output, d0_Summ;
	//std::vector<cl::Buffer> d_offsetV(num_devices_context);
	std::vector<cl::Buffer> d_z(num_devices_context);
	std::vector<cl::Buffer> d_xy(num_devices_context);
	//std::vector<cl::Buffer> d_y(num_devices_context);
	//std::vector<cl::Buffer> d_angles(num_devices_context);
	std::vector<cl::Buffer> d_xcenter(num_devices_context);
	std::vector<cl::Buffer> d_ycenter(num_devices_context);
	std::vector<cl::Buffer> d_zcenter(num_devices_context);
	//std::vector<cl::Buffer> d_offsetH(num_devices_context);
	//std::vector<cl::Buffer> d_offsetB(num_devices_context);
	std::vector<cl::Buffer> d_V(num_devices_context);
	//std::vector<cl::Buffer> d_atten(num_devices_context);
	std::vector<cl::Image3D> d_atten(num_devices_context);
	std::vector<cl::Buffer> d_TOFCenter(num_devices_context);
	//std::vector<cl::Buffer> d_pseudos(num_devices_context);
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
	std::vector<cl::Buffer> d_meanV(num_devices_context);
	std::vector<cl::Image3D> d_inputImage(num_devices_context);
	std::vector<cl::Image2D> d_mask1(num_devices_context);
	std::vector<cl::Image2D> d_mask2(num_devices_context);
	std::vector<cl::Image3D> d_ImageX(num_devices_context);
	std::vector<cl::Image3D> d_ImageY(num_devices_context);
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
		//if (CT || inputScalars.SPECT) {
		//	d_angles[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * nAngles, NULL, &status);
		//	if (status != CL_SUCCESS) {
		//		getErrorString(status);
		//		return;
		//	}
		//}
		// Coordinates of the centers of the voxels
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
		//if (inputScalars.projector_type == 4) {
		//	d_offsetH[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.nProjections, NULL, &status);
		//	if (status != CL_SUCCESS) {
		//		getErrorString(status);
		//		return;
		//	}
		//	d_offsetB[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.nProjections, NULL, &status);
		//	if (status != CL_SUCCESS) {
		//		getErrorString(status);
		//		return;
		//	}
		//}
		d_V[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_V, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		d_TOFCenter[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.nBins, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		// Attenuation image
		//d_atten[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_atten, NULL, &status);
		cl::ImageFormat format;
		format.image_channel_order = CL_R;
		format.image_channel_data_type = CL_FLOAT;
		cl::size_type imX = inputScalars.Nx;
		cl::size_type imY = inputScalars.Ny;
		cl::size_type imZ = inputScalars.Nz;
		d_atten[i] = cl::Image3D(context, CL_MEM_READ_ONLY, format, imX, imY, imZ, 0, 0, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		//d_pseudos[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t) * inputScalars.pRows, NULL, &status);
		//if (status != CL_SUCCESS) {
		//	getErrorString(status);
		//	return;
		//}
		// Transaxial source/detector coordinates
		d_xy[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * numel_x, NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		// Axial source/detector coordinates
		if (CT && listmode != 1) {
			if (inputScalars.PITCH)
				d_z[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.nProjections * 6, NULL, &status);
			else
				d_z[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.nProjections * 2, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else {
			d_z[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_z, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		// Mask images for projector types 4 and 5
		if ((inputScalars.projector_type == 5 || inputScalars.projector_type == 4)) {
			cl::ImageFormat format;
			format.image_channel_order = CL_R;
			format.image_channel_data_type = CL_UNSIGNED_INT8;
			if (inputScalars.maskFP && inputScalars.fp == 1) {
				cl::size_type imX = inputScalars.size_x;
				cl::size_type imY = inputScalars.size_y;
				d_mask1[i] = cl::Image2D(context, CL_MEM_READ_ONLY, format, imX, imY, 0, NULL, &status);
			}
			else if (inputScalars.maskBP && inputScalars.fp == 2) {
				cl::size_type imX = inputScalars.Nx;
				cl::size_type imY = inputScalars.Ny;
				if (DEBUG) {
					mexPrintf("imX = %u\n", imX);
					mexPrintf("imY = %u\n", imY);
					mexEvalString("pause(.0001);");
				}
				d_mask1[i] = cl::Image2D(context, CL_MEM_READ_ONLY, format, imX, imY, 0, NULL, &status);
			}
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		// Integral image and mean values for projector type 5
		if (inputScalars.projector_type == 5) {
			cl::ImageFormat format;
			format.image_channel_order = CL_R;
			format.image_channel_data_type = CL_FLOAT;
			if (inputScalars.fp == 2) {
				cl::size_type imX = inputScalars.size_x;
				cl::size_type imY = inputScalars.size_y;
				cl::size_type imZ = inputScalars.nProjections;
				d_inputImage[i] = cl::Image3D(context, CL_MEM_READ_ONLY, format, imX + 1, imY + 1, imZ, 0, 0, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				if (inputScalars.meanBP) {
					d_meanV[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.nProjections, NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
			}
			else {
				cl::size_type imX = inputScalars.Nx;
				cl::size_type imY = inputScalars.Ny;
				cl::size_type imZ = inputScalars.Nz;
				//if (inputScalars.fp == 2) {
				//	imX = inputScalars.size_x;
				//	imY = inputScalars.size_y;
				//	imZ = inputScalars.nProjections;
				//}
				d_ImageY[i] = cl::Image3D(context, CL_MEM_READ_ONLY, format, imX + 1, imZ + 1, imY, 0, 0, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				//d_ImageX[i] = cl::Image3D(context, CL_MEM_READ_ONLY, format, imY + 1, imZ + 1, imX + imY, 0, 0, NULL, &status);
				d_ImageX[i] = cl::Image3D(context, CL_MEM_READ_ONLY, format, imY + 1, imZ + 1, imX, 0, 0, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				if (inputScalars.meanFP) {
					d_meanV[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * (inputScalars.Nx + inputScalars.Ny), NULL, &status);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
			}
		}
		// Current estimate/image when using forward projection or projection image for projector type 4
		else if (inputScalars.projector_type == 4 || inputScalars.fp == 1) {
			cl::ImageFormat format;
			format.image_channel_order = CL_R;
			format.image_channel_data_type = CL_FLOAT;
			cl::size_type imX = inputScalars.Nx;
			cl::size_type imY = inputScalars.Ny;
			cl::size_type imZ = inputScalars.Nz;
			if (inputScalars.fp == 2) {
				imX = inputScalars.size_x;
				imY = inputScalars.size_y;
				imZ = inputScalars.nProjections;
			}
			d_inputImage[i] = cl::Image3D(context, CL_MEM_READ_ONLY, format, imX, imY, imZ, 0, 0, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			//if (inputScalars.SPECT) {
			//	imX = inputScalars.cSizeX;
			//	imY = inputScalars.cSizeY;
			//	d_mask1[i] = cl::Image2D(context, CL_MEM_READ_ONLY, format, imX, imY);
			//	d_mask2[i] = cl::Image2D(context, CL_MEM_READ_ONLY, format, imX, imY);
			//}
		}
		else if (inputScalars.fp == 2) {
			d_rhs[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * size_rhs, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		// Output vector and sensitivity image
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
			if (no_norm == 0)
				d_Summ[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(cl_float) * im_dim, NULL, &status);
			else
				d_Summ[i] = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(cl_float), NULL, &status);
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
		if (inputScalars.randoms_correction == 1u) {
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
		if (inputScalars.normalization_correction == 1u) {
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
		if (inputScalars.scatter == 1u) {
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
			d_Sino[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float) * length[i] * inputScalars.nBins, NULL, &status);
		else
			d_Sino[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		if (inputScalars.precompute)
			d_lor[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[i], NULL, &status);
		else
			d_lor[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[i], NULL, &status);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		if (inputScalars.raw && listmode != 1) {
			d_xyindex[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint32_t), NULL, &status);
			d_zindex[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
			d_L[i] = cl::Buffer(context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[i] * 2, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else if (listmode != 1 && ((!CT && !inputScalars.SPECT && !inputScalars.PET) && inputScalars.subsets > 1)) {
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
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
	const float* Sino = (float*)mxGetSingles(mxGetCell(Sin, static_cast<mwIndex>(0)));
#else
	const float* Sino = (float*)mxGetData(mxGetCell(Sin, static_cast<mwIndex>(0)));
#endif

	for (cl_uint i = 0; i < num_devices_context; i++) {
		status = commandQueues[i].enqueueWriteBuffer(d_reko_type[i], CL_FALSE, 0, sizeof(cl_uchar), &inputScalars.fp);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		//if (CT || inputScalars.SPECT) {
		//	status = commandQueues[i].enqueueWriteBuffer(d_angles[i], CL_FALSE, 0, sizeof(float) * nAngles, angles);
		//	if (status != CL_SUCCESS) {
		//		getErrorString(status);
		//		return;
		//	}
		//}
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
		status = commandQueues[i].enqueueWriteBuffer(d_TOFCenter[i], CL_FALSE, 0, sizeof(float) * inputScalars.nBins, TOFCenter);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		if (inputScalars.attenuation_correction) {
			const cl::detail::size_t_array origin = { 0, 0, 0 };
			cl::detail::size_t_array region = { 0, 0, 0 };
			region[0] = inputScalars.Nx;
			region[1] = inputScalars.Ny;
			region[2] = inputScalars.Nz;
			status = commandQueues[i].enqueueWriteImage(d_atten[i], CL_FALSE, origin, region, 0, 0, atten);
			//status = commandQueues[i].enqueueWriteBuffer(d_atten[i], CL_FALSE, 0, sizeof(float) * size_atten, atten);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		//status = commandQueues[i].enqueueWriteBuffer(d_pseudos[i], CL_FALSE, 0, sizeof(uint32_t) * inputScalars.pRows, pseudos);
		//if (status != CL_SUCCESS) {
		//	getErrorString(status);
		//	return;
		//}
		status = commandQueues[i].enqueueWriteBuffer(d_xy[i], CL_FALSE, 0, sizeof(float) * numel_x, x);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		//return;
		if (CT && listmode != 1) {
			if (inputScalars.PITCH)
				status = commandQueues[i].enqueueWriteBuffer(d_z[i], CL_FALSE, 0, sizeof(float) * inputScalars.nProjections * 6, uv);
			else
				status = commandQueues[i].enqueueWriteBuffer(d_z[i], CL_FALSE, 0, sizeof(float) * inputScalars.nProjections * 2, uv);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else {
			status = commandQueues[i].enqueueWriteBuffer(d_z[i], CL_FALSE, 0, sizeof(float) * inputScalars.size_z, z_det);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		if ((inputScalars.projector_type == 5 || inputScalars.projector_type == 4)) {
			const cl::detail::size_t_array origin = { 0, 0, 0 };
			cl::detail::size_t_array region = { 1, 1, 1 };
			if (inputScalars.maskFP && inputScalars.fp == 1) {
				region[0] = inputScalars.size_x;
				region[1] = inputScalars.size_y;
				status = commandQueues[i].enqueueWriteImage(d_mask1[i], CL_FALSE, origin, region, 0, 0, mask);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			else if (inputScalars.maskBP && inputScalars.fp == 2) {
				region[0] = inputScalars.Nx;
				region[1] = inputScalars.Ny;
				if (DEBUG) {
					mexPrintf("size_rhs = %u\n", size_rhs);
					mexPrintf("region[0] = %u\n", region[0]);
					mexPrintf("region[1] = %u\n", region[1]);
					mexPrintf("region[2] = %u\n", region[2]);
					mexEvalString("pause(.0001);");
				}
				status = commandQueues[i].enqueueWriteImage(d_mask1[i], CL_FALSE, origin, region, 0, 0, mask);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
		}
		if (inputScalars.projector_type == 5) {
			const cl::detail::size_t_array origin = { 0, 0, 0 };
			cl::detail::size_t_array region = { 0, 0, 0 };
			if (inputScalars.fp == 2) {
				region[0] = inputScalars.size_x + 1;
				region[1] = inputScalars.size_y + 1;
				region[2] = inputScalars.nProjections;
				if (DEBUG) {
					mexPrintf("size_rhs = %u\n", size_rhs);
					mexPrintf("region[0] = %u\n", region[0]);
					mexPrintf("region[1] = %u\n", region[1]);
					mexPrintf("region[2] = %u\n", region[2]);
					mexEvalString("pause(.0001);");
				}
				status = commandQueues[i].enqueueWriteImage(d_inputImage[i], CL_FALSE, origin, region, 0, 0, integralX);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				if (inputScalars.meanBP) {
					status = commandQueues[i].enqueueWriteBuffer(d_meanV[i], CL_FALSE, 0, sizeof(float) * inputScalars.nProjections, meanVals);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
			}
			else {
				region[0] = inputScalars.Nx + 1;
				region[1] = inputScalars.Nz + 1;
				//region[2] = inputScalars.Ny + inputScalars.Nx;
				region[2] = inputScalars.Ny;
				if (DEBUG) {
					mexPrintf("size_rhs = %u\n", size_rhs);
					mexPrintf("region[0] = %u\n", region[0]);
					mexPrintf("region[1] = %u\n", region[1]);
					mexPrintf("region[2] = %u\n", region[2]);
					mexEvalString("pause(.0001);");
				}
				status = commandQueues[i].enqueueWriteImage(d_ImageX[i], CL_FALSE, origin, region, 0, 0, integralX);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				if (inputScalars.meanFP) {
					status = commandQueues[i].enqueueWriteBuffer(d_meanV[i], CL_FALSE, 0, sizeof(float) * (inputScalars.Nx + inputScalars.Ny), meanVals);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				region[0] = inputScalars.Ny + 1;
				region[1] = inputScalars.Nz + 1;
				region[2] = inputScalars.Nx;
				status = commandQueues[i].enqueueWriteImage(d_ImageY[i], CL_FALSE, origin, region, 0, 0, integralY);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
		}
		else if (inputScalars.projector_type == 4 || inputScalars.fp == 1) {
			const cl::detail::size_t_array origin = { 0, 0, 0 };
			cl::detail::size_t_array region = { 0, 0, 0 };
			if (inputScalars.fp == 2) {
				region[0] = inputScalars.size_x;
				region[1] = inputScalars.size_y;
				region[2] = inputScalars.nProjections;
			}
			else {
				region[0] = inputScalars.Nx;
				region[1] = inputScalars.Ny;
				region[2] = inputScalars.Nz;
			}
			status = commandQueues[i].enqueueWriteImage(d_inputImage[i], CL_FALSE, origin, region, 0, 0, rhs);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			//if (inputScalars.SPECT) {
			//	const cl::detail::size_t_array origin2 = { 0, 0, 0};
			//	const cl::detail::size_t_array region2 = { inputScalars.cSizeX, inputScalars.cSizeY, 1};
			//	status = commandQueues[i].enqueueWriteImage(d_mask1[i], CL_FALSE, origin2, region2, 0, 0, mask1);
			//	if (status != CL_SUCCESS) {
			//		getErrorString(status);
			//		return;
			//	}
			//	status = commandQueues[i].enqueueWriteImage(d_mask2[i], CL_FALSE, origin2, region2, 0, 0, mask2);
			//	if (status != CL_SUCCESS) {
			//		getErrorString(status);
			//		return;
			//	}
			//}
		}
		else if (inputScalars.fp == 2) {
			status = commandQueues[i].enqueueWriteBuffer(d_rhs[i], CL_FALSE, 0, sizeof(float) * size_rhs, rhs);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		if (DEBUG) {
			mexPrintf("Step 4 done\n");
			mexEvalString("pause(.0001);");
		}
		if (inputScalars.TOF && listmode != 2) {
			for (int64_t to = 0LL; to < inputScalars.nBins; to++) {
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
		if (DEBUG) {
			mexPrintf("Step 44 done\n");
			mexEvalString("pause(.0001);");
		}


		if (inputScalars.raw && listmode != 1) {
			if (DEBUG) {
				mexPrintf("Step 4420 done\n");
				mexPrintf("inputScalars.raw = %u\n", inputScalars.raw);
				mexPrintf("listmode = %u\n", listmode);
				mexEvalString("pause(.0001);");
			}
			return;
			status = commandQueues[i].enqueueWriteBuffer(d_xyindex[i], CL_FALSE, 0, sizeof(uint32_t), xy_index);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			status = commandQueues[i].enqueueWriteBuffer(d_zindex[i], CL_FALSE, 0, sizeof(uint16_t), z_index);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			status = commandQueues[i].enqueueWriteBuffer(d_L[i], CL_FALSE, 0, sizeof(uint16_t) * length[i] * 2, &L[cumsum[i] * 2]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		else if (listmode != 1 && ((!CT && !inputScalars.SPECT && !inputScalars.PET) && inputScalars.subsets > 1)){
			if (DEBUG) {
				mexPrintf("Step 4430 done\n");
				mexEvalString("pause(.0001);");
			}
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
			if (DEBUG) {
				mexPrintf("Step 4440 done\n");
				mexEvalString("pause(.0001);");
			}
			status = commandQueues[i].enqueueWriteBuffer(d_xyindex[i], CL_FALSE, 0, sizeof(uint32_t), xy_index);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (DEBUG) {
				mexPrintf("Step 4441 done\n");
				mexEvalString("pause(.0001);");
			}
			status = commandQueues[i].enqueueWriteBuffer(d_zindex[i], CL_FALSE, 0, sizeof(uint16_t), z_index);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (DEBUG) {
				mexPrintf("Step 4442 done\n");
				mexEvalString("pause(.0001);");
			}
			status = commandQueues[i].enqueueWriteBuffer(d_L[i], CL_FALSE, 0, sizeof(uint16_t), L);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (DEBUG) {
				mexPrintf("Step 4443 done\n");
				mexEvalString("pause(.0001);");
			}
		}
		if (DEBUG) {
			mexPrintf("Step 444 done\n");
			mexEvalString("pause(.0001);");
		}
		if (inputScalars.precompute)
			status = commandQueues[i].enqueueWriteBuffer(d_lor[i], CL_FALSE, 0, sizeof(uint16_t) * length[i], &lor1[cumsum[i]]);
		else
			status = commandQueues[i].enqueueWriteBuffer(d_lor[i], CL_FALSE, 0, sizeof(uint16_t), lor1);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		if (DEBUG) {
			mexPrintf("Step 4444 done\n");
			mexEvalString("pause(.0001);");
		}
		if (inputScalars.normalization_correction == 1u) {
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
		if (DEBUG) {
			mexPrintf("Step 44444 done\n");
			mexEvalString("pause(.0001);");
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
	const cl_float3 b = { inputScalars.bx, inputScalars.by, inputScalars.bz };
	//const cl_float2 crystalSize = { inputScalars.crystal_size_xy, inputScalars.crystal_size_z };
	cl_float3 d = { inputScalars.dx, inputScalars.dy, inputScalars.dz };
	const cl_uint3 d_N = { inputScalars.Nx, inputScalars.Ny, inputScalars.Nz };
	const cl_float3 bmax = { static_cast<float>(inputScalars.Nx) * inputScalars.dx + inputScalars.bx, 
		static_cast<float>(inputScalars.Ny) * inputScalars.dy + inputScalars.by, 
		static_cast<float>(inputScalars.Nz) * inputScalars.dz + inputScalars.bz };

	if (DEBUG) {
		mexPrintf("bmaxX = %f\n", bmax.s[0]);
		mexPrintf("bmaxY = %f\n", bmax.s[1]);
		mexPrintf("bmaxZ = %f\n", bmax.s[2]);
		mexPrintf("by = %f\n", inputScalars.by);
		mexPrintf("bz = %f\n", inputScalars.bz);
		mexPrintf("dx = %f\n", inputScalars.dx);
		mexPrintf("dy = %f\n", inputScalars.dy);
		mexPrintf("dz = %f\n", inputScalars.dz);
		mexEvalString("pause(.0001);");
	}

	if (inputScalars.projector_type == 4 || inputScalars.projector_type == 5) {
		kernel_.setArg(kernelInd++, d_N);
		kernel_.setArg(kernelInd++, b);
		kernel_.setArg(kernelInd++, inputScalars.size_x);
		kernel_.setArg(kernelInd++, inputScalars.size_y);
		kernel_.setArg(kernelInd++, dPitch);
	}
	if (inputScalars.projector_type == 4) {
		//if (!inputScalars.SPECT) {
		if (inputScalars.fp == 1 && !inputScalars.SPECT) {
			kernel_.setArg(kernelInd++, bmax);
			kernel_.setArg(kernelInd++, inputScalars.dL);
			kernel_.setArg(kernelInd++, inputScalars.d_Scale);
		}
		//}
		//else {
		//	kernel_.setArg(kernelInd++, inputScalars.Nx);
		//	kernel_.setArg(kernelInd++, inputScalars.Ny);
		//	kernel_.setArg(kernelInd++, inputScalars.Nz);
		//	kernel_.setArg(kernelInd++, inputScalars.dx);
		//	kernel_.setArg(kernelInd++, inputScalars.dy);
		//	kernel_.setArg(kernelInd++, inputScalars.dz);
		//}
		if ((inputScalars.fp == 2 || (inputScalars.fp == 1 && inputScalars.SPECT))) {
			kernel_.setArg(kernelInd++, d);
			kernel_.setArg(kernelInd++, kerroin4);
			//kernel_.setArg(kernelInd++, b);
			//if (inputScalars.fp == 2)
			//	kernel_.setArg(kernelInd++, inputScalars.detY);
			//if (inputScalars.SPECT) {
			//	status = kernel_.setArg(kernelInd++, inputScalars.cThickness);
			//	if (status != CL_SUCCESS) {
			//		getErrorString(status);
			//		return;
			//	}
			//}
		}
		//else {
		//	kernel_.setArg(kernelInd++, inputScalars.d_Scale);
		//	kernel_.setArg(kernelInd++, Nxy);
		//}
	}
	else if (inputScalars.projector_type == 5) {
		//kernel_.setArg(kernelInd++, inputScalars.Nx);
		//kernel_.setArg(kernelInd++, inputScalars.Ny);
		//kernel_.setArg(kernelInd++, inputScalars.Nz);
		//kernel_.setArg(kernelInd++, inputScalars.dL);
		//kernel_.setArg(kernelInd++, inputScalars.dx);
		//kernel_.setArg(kernelInd++, inputScalars.dy);
		//kernel_.setArg(kernelInd++, inputScalars.dz);
		kernel_.setArg(kernelInd++, d);
		kernel_.setArg(kernelInd++, inputScalars.d_Scale);
		if (inputScalars.fp == 2)
			kernel_.setArg(kernelInd++, inputScalars.dSizeBP);
		else
			kernel_.setArg(kernelInd++, inputScalars.dSize);
		//if (inputScalars.fp == 2) {
		//	//kernel_.setArg(kernelInd++, inputScalars.detY);
		//	//kernel_.setArg(kernelInd++, inputScalars.meanV);
		//	kernel_.setArg(kernelInd++, no_norm);
		//}
		//kernel_.setArg(kernelInd++, Nxy);
	}
	else {
		kernel_.setArg(kernelInd++, inputScalars.global_factor);
		kernel_.setArg(kernelInd++, epps);
		kernel_.setArg(kernelInd++, im_dim);
		kernel_.setArg(kernelInd++, d_N);
		kernel_.setArg(kernelInd++, d);
		kernel_.setArg(kernelInd++, b);
		kernel_.setArg(kernelInd++, bmax);
		kernel_.setArg(kernelInd++, inputScalars.size_x);
		kernel_.setArg(kernelInd++, inputScalars.det_per_ring);
		kernel_.setArg(kernelInd++, Nxy);
		kernel_.setArg(kernelInd++, inputScalars.fp);
		kernel_.setArg(kernelInd++, inputScalars.sigma_x);
		kernel_.setArg(kernelInd++, dPitch);
		if (inputScalars.projector_type == 2u || inputScalars.projector_type == 3u || 
			(inputScalars.projector_type == 1u && (inputScalars.precompute || (inputScalars.n_rays * inputScalars.n_rays3D) == 1))) {
			kernel_.setArg(kernelInd++, inputScalars.tube_width);
			kernel_.setArg(kernelInd++, inputScalars.bmin);
			kernel_.setArg(kernelInd++, inputScalars.bmax);
			kernel_.setArg(kernelInd++, inputScalars.Vmax);
		}
		//else if (inputScalars.projector_type == 1u && !inputScalars.precompute) {
		//	//kernel_.setArg(kernelInd++, inputScalars.dc_z);
		//	kernel_.setArg(kernelInd++, inputScalars.n_rays);
		//}
	}
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
			if (no_norm == 0)
				status = commandQueues[i].enqueueFillBuffer(d_Summ[i], zero, 0, sizeof(cl_float) * im_dim);
			else
				status = commandQueues[i].enqueueFillBuffer(d_Summ[i], zero, 0, sizeof(cl_float));
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
		if (inputScalars.randoms_correction == 1u) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
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
		if (inputScalars.scatter == 1u) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			float* scat = (float*)mxGetSingles(mxGetCell(getField(options, 0, "ScatterFB"), static_cast<mwIndex>(0)));
#else
			float* scat = (float*)mxGetData(mxGetCell(getField(options, 0, "ScatterFB"), 0));
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

	//std::vector<cl_ulong> globals(num_devices_context, 0ULL);

	for (cl_uint i = 0; i < num_devices_context; i++) {

		const uint64_t st = static_cast<uint64_t>(cumsum[i]);

		cl_uint kernelIndSubIter = kernelInd;

		size_t global_size;
		//const size_t global_size = length[i] + erotus;
		if (inputScalars.projector_type == 4 && (inputScalars.fp == 2 || (inputScalars.fp == 1 && inputScalars.SPECT)))
			global_size = inputScalars.Nx * inputScalars.Ny * inputScalars.Nz;
		else
			global_size = length[i];
		const uint64_t m_size = static_cast<uint64_t>(length[i]);
		size_t erotus[3];
		if (inputScalars.projector_type == 4 || inputScalars.projector_type == 5) {
			if (inputScalars.fp == 1 && !inputScalars.SPECT) {
				erotus[0] = inputScalars.size_x % local_size[0];
				erotus[1] = inputScalars.size_y % local_size[1];
			}
			else {
				erotus[0] = inputScalars.Nx % local_size[0];
				erotus[1] = inputScalars.Ny % local_size[1];
				//erotus[2] = inputScalars.Nz % NVOXELS;
			}
			if (erotus[0] > 0)
				erotus[0] = (local_size[0] - erotus[0]);
			if (erotus[1] > 0)
				erotus[1] = (local_size[1] - erotus[1]);
			//if (erotus[2] > 0)
			//	erotus[2] = (NVOXELS - erotus[2]);
		}
		else {
			//erotus[0] = length[i] % local_size[0];
			if ((CT || inputScalars.SPECT || inputScalars.PET) && listmode == 0) {
				erotus[0] = inputScalars.size_x % local_size[0];
				erotus[1] = inputScalars.size_y % local_size[1];
				if (erotus[1] > 0)
					erotus[1] = (local_size[1] - erotus[1]);
			}
			else
				erotus[0] = global_size % local_size[0];

			if (erotus[0] > 0)
				erotus[0] = (local_size[0] - erotus[0]);
		}

		cl::NDRange local(local_size[0]);
		cl::NDRange global(global_size);
		if (inputScalars.projector_type == 4 || inputScalars.projector_type == 5) {
			if (inputScalars.fp == 1 && !inputScalars.SPECT)
				global = { inputScalars.size_x + erotus[0], inputScalars.size_y + erotus[1], static_cast<cl::size_type>(inputScalars.nProjections) };
			else
				if (inputScalars.projector_type == 4)
					global = { inputScalars.Nx + erotus[0], inputScalars.Ny + erotus[1], (inputScalars.Nz + NVOXELS - 1) / NVOXELS };
				else
					global = { inputScalars.Nx + erotus[0], inputScalars.Ny + erotus[1], inputScalars.Nz };
				//global = { inputScalars.Nx + erotus[0], inputScalars.Ny + erotus[1], (inputScalars.Nz + erotus[2]) / NVOXELS };
			local = { local_size[0] , local_size[1] };
		}
		else if ((CT || inputScalars.SPECT || inputScalars.PET) && listmode == 0) {
			global = { inputScalars.size_x + erotus[0], inputScalars.size_y + erotus[1], static_cast<cl::size_type>(inputScalars.nProjections) };
			local = { local_size[0] , local_size[1] };
		}
		else
			global = { length[i] + erotus[0], 1, 1 };
		//else if (listmode == 0)
		//	global = { inputScalars.size_x + erotus[0], inputScalars.size_y, static_cast<cl::size_type>(inputScalars.nProjections) };
			//global = { inputScalars.Nx, inputScalars.Ny, inputScalars.Nz };
		//	global = Nx * Ny * Nz;

		// These are needed to make sure the forward projection data is correctly combined
		//if (i > 0)
		//	globals[i] = globals[i - 1ULL] + m_size;
		//else
		//	globals[i] = m_size;

		if (DEBUG) {
			mexPrintf("global[0] = %u\n", global[0]);
			mexPrintf("local[0] = %u\n", local[0]);
			mexPrintf("local[1] = %u\n", local[1]);
			//if (inputScalars.projector_type == 4 || inputScalars.projector_type == 5) {
				mexPrintf("global[1] = %u\n", global[1]);
				mexPrintf("global[2] = %u\n", global[2]);
			//}
			mexPrintf("erotus[0] = %u\n", erotus[0]);
			mexPrintf("erotus[1] = %u\n", erotus[1]);
			mexPrintf("global.dimensions() = %u\n", global.dimensions());
			mexPrintf("local.dimensions() = %u\n", local.dimensions());
			mexPrintf("m_size = %u\n", m_size);
			mexPrintf("inputScalars.Nz = %u\n", inputScalars.Nz);
			mexPrintf("size_output = %u\n", size_output);
			mexPrintf("size_rhs = %u\n", size_rhs);
			mexPrintf("size_x = %u\n", inputScalars.size_x);
			mexPrintf("size_y = %u\n", inputScalars.size_y);
			mexPrintf("num_devices_context = %u\n", num_devices_context);
			mexPrintf("inputScalars.fp = %u\n", inputScalars.fp);
			mexPrintf("listmode = %u\n", listmode);
			mexPrintf("maskBP = %u\n", inputScalars.maskBP);
			mexPrintf("st = %u\n", st);
			mexPrintf("im_dim = %u\n", im_dim);
			mexPrintf("no_norm = %u\n", no_norm);
			mexPrintf("NVOXELS = %u\n", NVOXELS);
			mexPrintf("length[0] * nBins = %u\n", length[0] * inputScalars.nBins);
			mexEvalString("pause(.0001);");
		}

		//if (inputScalars.projector_type == 4 || inputScalars.projector_type == 5) {
		if ((inputScalars.maskFP && inputScalars.fp == 1) || (inputScalars.maskBP && inputScalars.fp == 2)) {
			status = kernel_.setArg(kernelIndSubIter++, d_mask1[i]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
		}
		//}
		if (inputScalars.projector_type == 4) {
			kernel_.setArg(kernelIndSubIter++, d_inputImage[i]);
			//kernel_.setArg(kernelIndSubIter++, m_size);
			status = kernel_.setArg(kernelIndSubIter++, d_output[i]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if ((inputScalars.fp == 2 || (inputScalars.fp == 1 && inputScalars.SPECT)))
				kernel_.setArg(kernelIndSubIter++, d_Summ[i]);
			kernel_.setArg(kernelIndSubIter++, d_xy[i]);
			status = kernel_.setArg(kernelIndSubIter++, d_z[i]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if ((inputScalars.fp == 2 || (inputScalars.fp == 1 && inputScalars.SPECT)))
				kernel_.setArg(kernelIndSubIter++, no_norm);
			kernel_.setArg(kernelIndSubIter++, inputScalars.nProjections);
		}
		else if (inputScalars.projector_type == 5) {
			//if ((inputScalars.fp == 2 || (inputScalars.fp == 1 && inputScalars.SPECT)))
			//	kernel_.setArg(kernelIndSubIter++, no_norm);
			status = kernel_.setArg(kernelIndSubIter++, d_xy[i]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (DEBUG) {
				mexPrintf("step 1\n");
				mexEvalString("pause(.0001);");
			}
			//if (inputScalars.fp == 2)
			//	status = kernel_.setArg(kernelIndSubIter++, d_offsetV[i]);
			//else
			//	kernel_.setArg(kernelIndSubIter++, d_z[i]);
			status = kernel_.setArg(kernelIndSubIter++, d_z[i]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (DEBUG) {
				mexPrintf("step 2\n");
				mexEvalString("pause(.0001);");
			}
			if (inputScalars.fp == 2) {
				status = kernel_.setArg(kernelIndSubIter++, d_inputImage[i]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			else {
				status = kernel_.setArg(kernelIndSubIter++, d_ImageY[i]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				status = kernel_.setArg(kernelIndSubIter++, d_ImageX[i]);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
			}
			if (DEBUG) {
				mexPrintf("step 3\n");
				mexEvalString("pause(.0001);");
			}
			//kernel_.setArg(kernelIndSubIter++, inputScalars.subsets);
			//status = kernel_.setArg(kernelIndSubIter++, d_angles[i]);
			//if (status != CL_SUCCESS) {
			//	getErrorString(status);
			//	return;
			//}
			//kernel_.setArg(kernelIndSubIter++, inputScalars.size_y);
			//kernel_.setArg(kernelIndSubIter++, inputScalars.dPitch);
			//kernel_.setArg(kernelIndSubIter++, d_xyindex[i]);
			//kernel_.setArg(kernelIndSubIter++, d_zindex[i]);
			//kernel_.setArg(kernelIndSubIter++, d_rhs[i]);
			//kernel_.setArg(kernelIndSubIter++, m_size);
			status = kernel_.setArg(kernelIndSubIter++, d_output[i]);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				return;
			}
			if (DEBUG) {
				mexPrintf("step 4\n");
				mexEvalString("pause(.0001);");
			}
			if ((inputScalars.fp == 2 || (inputScalars.fp == 1 && inputScalars.SPECT)))
				kernel_.setArg(kernelIndSubIter++, d_Summ[i]);
			if ((inputScalars.fp == 1 && inputScalars.meanFP) || (inputScalars.meanBP && inputScalars.fp == 2))
				kernel_.setArg(kernelIndSubIter++, d_meanV[i]);
			if ((inputScalars.fp == 2 || (inputScalars.fp == 1 && inputScalars.SPECT)))
				kernel_.setArg(kernelInd++, no_norm);
			kernel_.setArg(kernelIndSubIter++, inputScalars.nProjections);
			if (DEBUG) {
				mexPrintf("step 5\n");
				mexEvalString("pause(.0001);");
			}
		}
		else {
			if (inputScalars.TOF)
				kernel_.setArg(kernelIndSubIter++, d_TOFCenter[i]);
			if (inputScalars.projector_type == 2u || inputScalars.projector_type == 3u) {
				kernel_.setArg(kernelIndSubIter++, d_xcenter[i]);
				kernel_.setArg(kernelIndSubIter++, d_ycenter[i]);
				kernel_.setArg(kernelIndSubIter++, d_zcenter[i]);
				kernel_.setArg(kernelIndSubIter++, d_V[i]);
			}
			kernel_.setArg(kernelIndSubIter++, d_reko_type[i]);
			kernel_.setArg(kernelIndSubIter++, zero);
			if (!CT && inputScalars.attenuation_correction)
				kernel_.setArg(kernelIndSubIter++, d_atten[i]);
			if (CT || inputScalars.SPECT || inputScalars.PET) {
				kernel_.setArg(kernelIndSubIter++, inputScalars.size_y);
				kernel_.setArg(kernelIndSubIter++, inputScalars.nProjections);
			}
			kernel_.setArg(kernelIndSubIter++, d_xy[i]);
			kernel_.setArg(kernelIndSubIter++, d_z[i]);
			if (DEBUG) {
				mexPrintf("step 1\n");
				mexEvalString("pause(.0001);");
			}
			if (DEBUG) {
				mexPrintf("step 2\n");
				mexEvalString("pause(.0001);");
			}
			kernel_.setArg(kernelIndSubIter++, d_norm[i]);
			kernel_.setArg(kernelIndSubIter++, d_scat[i]);
			kernel_.setArg(kernelIndSubIter++, d_Summ[i]);
			if (inputScalars.precompute)
				kernel_.setArg(kernelIndSubIter++, d_lor[i]);
			if (inputScalars.subsets > 1 && inputScalars.subsetType < 8) {
				kernel_.setArg(kernelIndSubIter++, d_xyindex[i]);
				kernel_.setArg(kernelIndSubIter++, d_zindex[i]);
			}
			if (inputScalars.raw)
				kernel_.setArg(kernelIndSubIter++, d_L[i]);
			kernel_.setArg(kernelIndSubIter++, d_Sino[i]);
			kernel_.setArg(kernelIndSubIter++, d_sc_ra[i]);
			if (DEBUG) {
				mexPrintf("step 3\n");
				mexEvalString("pause(.0001);");
			}
			if (inputScalars.fp == 1)
				kernel_.setArg(kernelIndSubIter++, d_inputImage[i]);
			else
				kernel_.setArg(kernelIndSubIter++, d_rhs[i]);
			if (DEBUG) {
				mexPrintf("step 4\n");
				mexEvalString("pause(.0001);");
			}
			kernel_.setArg(kernelIndSubIter++, d_output[i]);
			kernel_.setArg(kernelIndSubIter++, no_norm);
			kernel_.setArg(kernelIndSubIter++, m_size);
			kernel_.setArg(kernelIndSubIter++, st);
			if (DEBUG) {
				mexPrintf("step 5\n");
				mexEvalString("pause(.0001);");
			}
		}
		status = commandQueues[i].enqueueNDRangeKernel(kernel_, cl::NDRange(), global, local, NULL, &events[i][0]);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		else if (DEBUG) {
			if (inputScalars.fp == 1)
				mexPrintf("Forward projection kernel launched successfully\n");
			else if (inputScalars.fp == 2)
				mexPrintf("Backprojection kernel launched successfully\n");
			else
				mexPrintf("Kernel launched successfully\n");
			mexEvalString("pause(.0001);");
		}
	}
	for (cl_uint i = 0ULL; i < num_devices_context; i++) {
		commandQueues[i].finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
	}
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
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
				if (no_norm == 0)
					status = commandQueues[i].enqueueReadBuffer(d_Summ[i], CL_TRUE, 0, sizeof(float) * im_dim, testi_summ[i - 1u].data(), &events[i]);
				else
					status = commandQueues[i].enqueueReadBuffer(d_Summ[i], CL_TRUE, 0, sizeof(float), testi_summ[i - 1u].data(), &events[i]);
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
		if (no_norm == 0)
			status = commandQueues[0].enqueueReadBuffer(d_Summ[0], CL_FALSE, 0, sizeof(float) * im_dim, normalizer_f);
		else
			status = commandQueues[0].enqueueReadBuffer(d_Summ[0], CL_FALSE, 0, sizeof(float), normalizer_f);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = commandQueues[0].enqueueReadBuffer(d_output[0], CL_FALSE, 0, sizeof(float) * size_output, output_f);
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