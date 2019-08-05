/**************************************************************************
* Matrix free computations for OMEGA.
* In this file the OpenCL buffers are created, calls to other necessary 
* functions are made and the OpenCL kernels are launched. This file 
* contains the code for the matrix-free reconstructions in OMEGA using the
* implementation 2.
* Unlike the non-OpenCL versions, this one uses (32-bit) floats and thus
* can be slightly more inaccurate.
*
* Copyright (C) 2019  Ville-Veikko Wettenhovi
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <https://www.gnu.org/licenses/>.
***************************************************************************/
#include "functions.hpp"

// Use ArrayFire namespace for convenience
using namespace af;

// Main reconstruction function
void reconstruction_AF_matrixfree(const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, const mxArray* Sin, const uint32_t Nx, const uint32_t Ny,
	const uint32_t Nz, const uint32_t Niter, const mxArray* options, const float dx, const float dy, const float dz, const float bx, const float by, const float bz,
	const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices, const uint32_t* pituus, const size_t koko_l,
	const uint32_t* xy_index, const uint16_t* z_index, const uint32_t size_x, const uint32_t TotSinos, mxArray* cell, const mwSize* dimmi, const bool verbose,
	const uint32_t attenuation_correction, const uint32_t normalization, const float* atten, const size_t size_atten, const float* norm, const size_t size_norm, const uint32_t subsets, const float epps, const uint8_t* rekot,
	const char* k_path, const size_t size_rekot, const uint32_t Nt, const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L,
	const uint8_t raw, const size_t size_z, const bool osem_bool, const char* fileName, const bool force_build, const float tube_width, const float crystal_size_z,
	const float* x_center, const float* y_center, const float* z_center, const size_t size_center_x, const size_t size_center_y, const size_t size_of_x,
	const size_t size_center_z, const uint32_t projector_type, const uint32_t device) {



	// Number of voxels
	const uint32_t Nxy = Nx * Ny;
	const uint32_t im_dim = Nxy * Nz;
	bool atomic_64bit = true;
	cl_uchar compute_norm_matrix = 1u;
	cl_float mem_portions;
	if (raw == 1u)
		mem_portions = 0.1f;
	else
		mem_portions = 0.2f;
	cl_float image_bytes = static_cast<cl_float>(Nx * Ny * Nz) * 8.f;

	uint32_t oo = 0u;
	size_t ll = 0ULL;

	// Create the output structs
	matlabArrays ArrayList;
	// Create the struct containing the reconstruction methods used
	RecMethods MethodList;
	// Same as above, but as cl_ variables
	RecMethodsOpenCL MethodListOpenCL;

	// Obtain the reconstruction methods used
	get_rec_methods(options, MethodList);
	OpenCLRecMethods(MethodList, MethodListOpenCL);

	// Create the MATLAB output arrays
	create_matlab_output(ArrayList, dimmi, MethodList, 4);

	// Initial value
	array x00(Nx*Ny*Nz, (float*)mxGetData(mxGetField(options, 0, "x0")), afHost);

	array pj3, Sino, apu_sum;

	// Are ML-methods used?
	bool mlem_bool = (MethodList.MLEM || MethodList.OSLMLEM) ? true : false;

	// Number of measurements at each subset
	std::vector<size_t> length(subsets);

	for (uint32_t kk = 0u; kk < subsets; kk++)
		length[kk] = pituus[kk + 1u] - pituus[kk];

	// Struct containing ArrayFire arrays containing the image estimates and other necessary vectors/matrices
	AF_im_vectors vec;
	// Struct containing beta-values for the MAP-methods
	Beta beta;
	// Struct containing the necessary variables for the priors
	Weighting w_vec;
	// Struct for TV data
	TVdata data;
	// Struct containing the OpenCL image estimates
	OpenCL_im_vectors vec_opencl;

	// Load the necessary data from the MATLAB input and form the necessary variables
	form_data_variables(vec, beta, w_vec, options, Nx, Ny, Nz, Niter, x00, im_dim, koko_l, MethodList, data, subsets);

	// Power factor for ACOSEM
	const float hh = 1.f / w_vec.h_ACOSEM;

	// Create the OpenCL context and command queue and assign the device
	cl_context af_context = afcl::getContext();
	cl_device_id af_device_id = afcl::getDeviceId();
	cl_command_queue af_queue = afcl::getQueue();

	cl_program program = NULL;

	cl_int status = CL_SUCCESS;

	cl_ulong mem;
	status = clGetDeviceInfo(af_device_id, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &mem, NULL);
	if ((static_cast<cl_float>(mem) * mem_portions) > image_bytes)
		compute_norm_matrix = 0u;

	std::vector<array> Summ;

	if (compute_norm_matrix == 0u)
		Summ.assign(subsets, constant(0.f, im_dim, 1));
	else
		Summ.assign(1ULL, constant(0.f, im_dim, 1));

	std::string new_filename;

	// Build the OpenCL program and save the binaries
	// Force the build regardless of whether the binaries already exist
	if (force_build) {
		status = SaveProgramBinary(verbose, k_path, af_context, af_device_id, fileName, program, atomic_64bit, device);
		if (status != CL_SUCCESS) {
			std::cerr << "Error while saving binary" << std::endl;
			return;
		}
	}
	else {
		FILE *fp = NULL;// = fopen(fileName, "rb");
		new_filename = fileName;
		new_filename += (std::to_string(device) + "_64atom.bin");
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
		fopen_s(&fp, new_filename.c_str(), "rb");
#else
		fp = fopen(new_filename.c_str(), "rb");
#endif
		if (fp == NULL) {
			new_filename = fileName;
			new_filename += (std::to_string(device) + ".bin");
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
			fopen_s(&fp, new_filename.c_str(), "rb");
#else
			fp = fopen(new_filename.c_str(), "rb");
#endif
			// If the binaries do not yet exist, create them
			if (fp == NULL) {
				status = SaveProgramBinary(verbose, k_path, af_context, af_device_id, fileName, program, atomic_64bit, device);
				if (status != CL_SUCCESS)
					return;
			}
			else {
				// If the binaries exist, load them
				status = CreateProgramFromBinary(af_context, af_device_id, fp, program);
				fclose(fp);
				if (status != CL_SUCCESS) {
					mexPrintf("Failed to load OpenCL binaries\n");
					clReleaseProgram(program);
					return;
				}
				atomic_64bit = false;
			}
		}
		else {
			// If the binaries exist, load them
			status = CreateProgramFromBinary(af_context, af_device_id, fp, program);
			fclose(fp);
			if (status != CL_SUCCESS) {
				mexPrintf("Failed to load OpenCL binaries\n");
				clReleaseProgram(program);
				return;
			}
			//else {
			//	mexPrintf("OpenCL binaries successfully loaded\n");
			//	mexEvalString("pause(.0001);");
			//}
		}
	}

	// Create the kernels
	cl_kernel kernel_ml = NULL, kernel = NULL, kernel_mramla = NULL;

	status = createKernels(kernel_ml, kernel, kernel_mramla, osem_bool, program, MethodList, w_vec, projector_type); 
	if (status != CL_SUCCESS) {
		mexPrintf("Failed to create kernels\n");
		clReleaseProgram(program);
		clReleaseKernel(kernel);
		clReleaseKernel(kernel_ml);
		clReleaseKernel(kernel_mramla);
		return;
	}

	cl_mem d_x, d_y, d_z, d_pseudos, d_atten, d_xcenter, d_ycenter, d_zcenter, d_norm;
	cl_mem* d_Summ;

	std::vector<cl_mem> d_lor(subsets);
	std::vector<cl_mem> d_L(subsets);
	std::vector<cl_mem> d_zindex(subsets);
	std::vector<cl_mem> d_xyindex(subsets);
	std::vector<cl_mem> d_Sino(subsets);

	float *apu = (float*)mxGetData(mxGetCell(Sin, 0));

	status = createAndWriteBuffers(d_x, d_y, d_z, d_lor, d_L, d_zindex, d_xyindex, d_Sino, size_x, size_z, TotSinos, size_atten, size_norm, prows, length, x, y, z_det, xy_index,
		z_index, lor1, L, apu, raw, af_context, subsets, pituus, atten, norm, pseudos, af_queue, d_atten, d_norm, d_pseudos, d_xcenter, d_ycenter, d_zcenter, x_center, y_center, 
		z_center, size_center_x, size_center_y, size_center_z, size_of_x, atomic_64bit);
	if (status != CL_SUCCESS) {
		clReleaseProgram(program);
		clReleaseKernel(kernel);
		clReleaseKernel(kernel_ml);
		clReleaseKernel(kernel_mramla);
		return;
	}

	// Initialize the OpenCL variables
	initialize_opencl_inputs(vec, vec_opencl, MethodList, mlem_bool, osem_bool, af_context, af_queue, im_dim);

	// Compute the prepass phase for MRAMLA, MBSREM, RBI, COSEM, ACOSEM or ECOSEM if applicable
	if ((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBI || MethodList.RBIMAP) && w_vec.MBSREM_prepass ||
		MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0) {

		// Set the kernel parameters that do not change
		clSetKernelArg(kernel_mramla, 0, sizeof(uint32_t), &det_per_ring);
		clSetKernelArg(kernel_mramla, 1, sizeof(uint8_t), &raw);
		clSetKernelArg(kernel_mramla, 2, sizeof(cl_mem), &d_pseudos);
		clSetKernelArg(kernel_mramla, 3, sizeof(uint32_t), &prows);
		clSetKernelArg(kernel_mramla, 4, sizeof(float), &hh);
		clSetKernelArg(kernel_mramla, 5, sizeof(uint32_t), &Nx);
		clSetKernelArg(kernel_mramla, 6, sizeof(uint32_t), &Ny);
		clSetKernelArg(kernel_mramla, 7, sizeof(uint32_t), &Nz);
		clSetKernelArg(kernel_mramla, 8, sizeof(float), &dz);
		clSetKernelArg(kernel_mramla, 9, sizeof(float), &dx);
		clSetKernelArg(kernel_mramla, 10, sizeof(float), &dy);
		clSetKernelArg(kernel_mramla, 11, sizeof(float), &bz);
		clSetKernelArg(kernel_mramla, 12, sizeof(float), &bx);
		clSetKernelArg(kernel_mramla, 13, sizeof(float), &by);
		clSetKernelArg(kernel_mramla, 14, sizeof(float), &bzb);
		clSetKernelArg(kernel_mramla, 15, sizeof(float), &maxxx);
		clSetKernelArg(kernel_mramla, 16, sizeof(float), &maxyy);
		clSetKernelArg(kernel_mramla, 17, sizeof(float), &zmax);
		clSetKernelArg(kernel_mramla, 18, sizeof(float), &NSlices);
		clSetKernelArg(kernel_mramla, 19, sizeof(cl_mem), &d_x);
		clSetKernelArg(kernel_mramla, 20, sizeof(cl_mem), &d_y);
		clSetKernelArg(kernel_mramla, 21, sizeof(cl_mem), &d_z);
		clSetKernelArg(kernel_mramla, 22, sizeof(uint32_t), &size_x);
		clSetKernelArg(kernel_mramla, 23, sizeof(uint32_t), &TotSinos);
		clSetKernelArg(kernel_mramla, 24, sizeof(uint32_t), &attenuation_correction);
		clSetKernelArg(kernel_mramla, 25, sizeof(uint32_t), &normalization);
		clSetKernelArg(kernel_mramla, 26, sizeof(cl_mem), &d_atten);
		clSetKernelArg(kernel_mramla, 27, sizeof(cl_mem), &d_norm);
		clSetKernelArg(kernel_mramla, 28, sizeof(float), &epps);
		clSetKernelArg(kernel_mramla, 29, sizeof(uint32_t), &im_dim);
		clSetKernelArg(kernel_mramla, 30, sizeof(MethodListOpenCL), &MethodListOpenCL);
		clSetKernelArg(kernel_mramla, 31, sizeof(uint32_t), &Nxy);
		clSetKernelArg(kernel_mramla, 32, sizeof(float), &tube_width);
		clSetKernelArg(kernel_mramla, 33, sizeof(float), &crystal_size_z);
		clSetKernelArg(kernel_mramla, 34, sizeof(cl_mem), &d_xcenter);
		clSetKernelArg(kernel_mramla, 35, sizeof(cl_mem), &d_ycenter);
		clSetKernelArg(kernel_mramla, 36, sizeof(cl_mem), &d_zcenter);

		uint32_t alku = 0u;

		//if (w_vec.MBSREM_prepass || compute_norm_matrix == 0u)
		//	Summ = constant(0.f, im_dim, subsets);
		//if (w_vec.MBSREM_prepass || compute_norm_matrix == 0u)

		//else
		//	Summ = constant(0.f, 1, 1);

		if ((MethodList.MRAMLA || MethodList.MBSREM) && w_vec.MBSREM_prepass)
			w_vec.Amin = constant(0.f, koko_l, 1);
		// Run the prepass phase
		MRAMLA_prepass(subsets, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, x00, vec.C_co,
			vec.C_aco, vec.C_osl, alku, kernel_mramla, d_L, raw, MethodListOpenCL, length, atomic_64bit, compute_norm_matrix);

		//if (w_vec.MBSREM_prepass)
		//	w_vec.D = sum(Summ, 1);


		if (verbose) {
			mexPrintf("MRAMLA & COSEM prepass completed\n");
			mexEvalString("pause(.0001);");
		}
	}

	if (MethodList.MRAMLA || MethodList.MBSREM) {
		pj3 = w_vec.D / static_cast<float>(subsets);
	}

	//if (!w_vec.MBSREM_prepass && compute_norm_matrix == 0u)
	//	Summ = constant(0.f, im_dim, subsets);

	if (osem_bool) {
		// Set the kernel parameters that do not change
		clSetKernelArg(kernel, 0, sizeof(MethodListOpenCL), &MethodListOpenCL);
		clSetKernelArg(kernel, 1, sizeof(uint8_t), &raw);
		clSetKernelArg(kernel, 2, sizeof(float), &hh);
		clSetKernelArg(kernel, 3, sizeof(uint32_t), &Nx);
		clSetKernelArg(kernel, 4, sizeof(uint32_t), &Ny);
		clSetKernelArg(kernel, 5, sizeof(uint32_t), &Nz);
		clSetKernelArg(kernel, 6, sizeof(float), &dz);
		clSetKernelArg(kernel, 7, sizeof(float), &dx);
		clSetKernelArg(kernel, 8, sizeof(float), &dy);
		clSetKernelArg(kernel, 9, sizeof(float), &bz);
		clSetKernelArg(kernel, 10, sizeof(float), &bx);
		clSetKernelArg(kernel, 11, sizeof(float), &by);
		clSetKernelArg(kernel, 12, sizeof(float), &bzb);
		clSetKernelArg(kernel, 13, sizeof(float), &maxxx);
		clSetKernelArg(kernel, 14, sizeof(float), &maxyy);
		clSetKernelArg(kernel, 15, sizeof(float), &zmax);
		clSetKernelArg(kernel, 16, sizeof(float), &NSlices);
		clSetKernelArg(kernel, 17, sizeof(cl_mem), &d_x);
		clSetKernelArg(kernel, 18, sizeof(cl_mem), &d_y);
		clSetKernelArg(kernel, 19, sizeof(cl_mem), &d_z);
		clSetKernelArg(kernel, 20, sizeof(uint32_t), &size_x);
		clSetKernelArg(kernel, 21, sizeof(uint32_t), &TotSinos);
		clSetKernelArg(kernel, 22, sizeof(uint32_t), &attenuation_correction);
		clSetKernelArg(kernel, 23, sizeof(uint32_t), &normalization);
		clSetKernelArg(kernel, 24, sizeof(cl_mem), &d_atten);
		clSetKernelArg(kernel, 25, sizeof(cl_mem), &d_norm);
		clSetKernelArg(kernel, 26, sizeof(float), &epps);
		clSetKernelArg(kernel, 27, sizeof(uint32_t), &im_dim);
		clSetKernelArg(kernel, 28, sizeof(cl_mem), &d_pseudos);
		clSetKernelArg(kernel, 29, sizeof(uint32_t), &prows);
		clSetKernelArg(kernel, 30, sizeof(uint32_t), &Nxy);
		clSetKernelArg(kernel, 31, sizeof(uint32_t), &det_per_ring);
		clSetKernelArg(kernel, 32, sizeof(float), &tube_width);
		clSetKernelArg(kernel, 33, sizeof(float), &crystal_size_z);
		clSetKernelArg(kernel, 34, sizeof(cl_mem), &d_xcenter);
		clSetKernelArg(kernel, 35, sizeof(cl_mem), &d_ycenter);
		clSetKernelArg(kernel, 36, sizeof(cl_mem), &d_zcenter);
	}

	if (mlem_bool) {
		// Set the kernel parameters that do not change
		clSetKernelArg(kernel_ml, 0, sizeof(MethodListOpenCL), &MethodListOpenCL);
		clSetKernelArg(kernel_ml, 1, sizeof(uint32_t), &det_per_ring);
		clSetKernelArg(kernel_ml, 2, sizeof(uint8_t), &raw);
		clSetKernelArg(kernel_ml, 3, sizeof(uint32_t), &Nx);
		clSetKernelArg(kernel_ml, 4, sizeof(uint32_t), &Ny);
		clSetKernelArg(kernel_ml, 5, sizeof(uint32_t), &Nz);
		clSetKernelArg(kernel_ml, 6, sizeof(float), &dz);
		clSetKernelArg(kernel_ml, 7, sizeof(float), &dx);
		clSetKernelArg(kernel_ml, 8, sizeof(float), &dy);
		clSetKernelArg(kernel_ml, 9, sizeof(float), &bz);
		clSetKernelArg(kernel_ml, 10, sizeof(float), &bx);
		clSetKernelArg(kernel_ml, 11, sizeof(float), &by);
		clSetKernelArg(kernel_ml, 12, sizeof(float), &bzb);
		clSetKernelArg(kernel_ml, 13, sizeof(float), &maxxx);
		clSetKernelArg(kernel_ml, 14, sizeof(float), &maxyy);
		clSetKernelArg(kernel_ml, 15, sizeof(float), &zmax);
		clSetKernelArg(kernel_ml, 16, sizeof(float), &NSlices);
		clSetKernelArg(kernel_ml, 17, sizeof(cl_mem), &d_x);
		clSetKernelArg(kernel_ml, 18, sizeof(cl_mem), &d_y);
		clSetKernelArg(kernel_ml, 19, sizeof(cl_mem), &d_z);
		clSetKernelArg(kernel_ml, 20, sizeof(uint32_t), &size_x);
		clSetKernelArg(kernel_ml, 21, sizeof(uint32_t), &TotSinos);
		clSetKernelArg(kernel_ml, 22, sizeof(uint32_t), &attenuation_correction);
		clSetKernelArg(kernel_ml, 23, sizeof(uint32_t), &normalization);
		clSetKernelArg(kernel_ml, 24, sizeof(cl_mem), &d_atten);
		clSetKernelArg(kernel_ml, 25, sizeof(cl_mem), &d_norm);
		clSetKernelArg(kernel_ml, 26, sizeof(float), &epps);
		clSetKernelArg(kernel_ml, 27, sizeof(uint32_t), &im_dim);
		clSetKernelArg(kernel_ml, 28, sizeof(cl_mem), &d_pseudos);
		clSetKernelArg(kernel_ml, 29, sizeof(uint32_t), &prows);
		clSetKernelArg(kernel_ml, 30, sizeof(uint32_t), &Nxy);
		clSetKernelArg(kernel_ml, 31, sizeof(float), &tube_width);
		clSetKernelArg(kernel_ml, 32, sizeof(float), &crystal_size_z);
		clSetKernelArg(kernel_ml, 33, sizeof(cl_mem), &d_xcenter);
		clSetKernelArg(kernel_ml, 34, sizeof(cl_mem), &d_ycenter);
		clSetKernelArg(kernel_ml, 35, sizeof(cl_mem), &d_zcenter);
	}

	// Loop through each time-step
	for (uint32_t tt = 0u; tt < Nt; tt++) {

		cl_uchar no_norm = 0u;

		// Load the measurement data from the cell array
		if (tt > 0u) {
			float* apu = (float*)mxGetData(mxGetCell(Sin, tt));
			for (uint32_t kk = 0u; kk < subsets; kk++) {
				status = clEnqueueWriteBuffer(af_queue, d_Sino[kk], CL_TRUE, 0, sizeof(float) * length[kk], &apu[pituus[kk]], 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
				}
			}
		}

		// Compute values needed for MBSREM and MRAMLA
		if (MethodList.MBSREM || MethodList.MRAMLA) {
			Sino = array(koko_l, (float*)mxGetData(mxGetCell(Sin, tt)), afHost);
			if (*w_vec.U == 0.f) {
				array temppi = (max)(Sino / w_vec.Amin);
				temppi.host(w_vec.U);
			}
			w_vec.epsilon_mramla = MBSREM_epsilon(Sino, epps);
		}

		// Loop through each iteration
		for (uint32_t iter = 0u; iter < Niter; iter++) {

			// Compute MLEM separately
			//if (mlem_bool) {

			//	if (compute_norm_matrix == 1u) {
			//		if (atomic_64bit)
			//			Summ = constant(0ULL, im_dim, 1, u64);
			//		else
			//			Summ = constant(0.f, im_dim, 1);
			//		d_Summ = Summ.device<cl_mem>();
			//	}
			//	else if (compute_norm_matrix == 0u && osem_bool) {
			//		if (atomic_64bit) {
			//			apu_sum = sum(Summ, 1);
			//			apu_sum = apu_sum.as(u64);
			//			d_Summ = apu_sum.device<cl_mem>();
			//		}
			//		else {
			//			apu_sum = sum(Summ, 1);
			//			d_Summ = apu_sum.device<cl_mem>();
			//		}
			//	}
			//	else {
			//		if (atomic_64bit) {
			//			apu_sum = Summ.as(u64);
			//			d_Summ = apu_sum.device<cl_mem>();
			//		}
			//		else
			//			d_Summ = Summ.device<cl_mem>();
			//	}

			//	// Update the OpenCL inputs for this iteration (image estimates)
			//	update_opencl_inputs(vec, vec_opencl, im_dim, MethodList, true, 0, af_context, 0, af_queue);

			//	for (uint32_t kk = 0u; kk < subsets; kk++) {
			//		clSetKernelArg(kernel_ml, 36, sizeof(cl_mem), &d_Sino[kk]);
			//		clSetKernelArg(kernel_ml, 37, sizeof(cl_mem), &d_L[kk]);
			//		clSetKernelArg(kernel_ml, 38, sizeof(cl_mem), &d_xyindex[kk]);
			//		clSetKernelArg(kernel_ml, 39, sizeof(cl_mem), &d_zindex[kk]);
			//		clSetKernelArg(kernel_ml, 40, sizeof(cl_mem), d_Summ);
			//		clSetKernelArg(kernel_ml, 41, sizeof(cl_mem), &d_lor[kk]);
			//		clSetKernelArg(kernel_ml, 42, sizeof(cl_mem), vec_opencl.d_MLEM);
			//		clSetKernelArg(kernel_ml, 43, sizeof(cl_mem), vec_opencl.d_MRP_MLEM);
			//		clSetKernelArg(kernel_ml, 44, sizeof(cl_mem), vec_opencl.d_Quad_MLEM);
			//		clSetKernelArg(kernel_ml, 45, sizeof(cl_mem), vec_opencl.d_L_MLEM);
			//		clSetKernelArg(kernel_ml, 46, sizeof(cl_mem), vec_opencl.d_FMH_MLEM);
			//		clSetKernelArg(kernel_ml, 47, sizeof(cl_mem), vec_opencl.d_Weighted_MLEM);
			//		clSetKernelArg(kernel_ml, 48, sizeof(cl_mem), vec_opencl.d_TV_MLEM);
			//		clSetKernelArg(kernel_ml, 49, sizeof(cl_mem), vec_opencl.d_AD_MLEM);
			//		clSetKernelArg(kernel_ml, 50, sizeof(cl_mem), vec_opencl.d_APLS_MLEM);
			//		clSetKernelArg(kernel_ml, 51, sizeof(cl_mem), vec_opencl.d_TGV_MLEM);
			//		clSetKernelArg(kernel_ml, 52, sizeof(cl_mem), vec_opencl.d_MLEM_rhs);
			//		clSetKernelArg(kernel_ml, 53, sizeof(cl_mem), vec_opencl.d_MRP_MLEM_rhs);
			//		clSetKernelArg(kernel_ml, 54, sizeof(cl_mem), vec_opencl.d_Quad_MLEM_rhs);
			//		clSetKernelArg(kernel_ml, 55, sizeof(cl_mem), vec_opencl.d_L_MLEM_rhs);
			//		clSetKernelArg(kernel_ml, 56, sizeof(cl_mem), vec_opencl.d_FMH_MLEM_rhs);
			//		clSetKernelArg(kernel_ml, 57, sizeof(cl_mem), vec_opencl.d_Weighted_MLEM_rhs);
			//		clSetKernelArg(kernel_ml, 58, sizeof(cl_mem), vec_opencl.d_TV_MLEM_rhs);
			//		clSetKernelArg(kernel_ml, 59, sizeof(cl_mem), vec_opencl.d_AD_MLEM_rhs);
			//		clSetKernelArg(kernel_ml, 60, sizeof(cl_mem), vec_opencl.d_APLS_MLEM_rhs);
			//		clSetKernelArg(kernel_ml, 61, sizeof(cl_mem), vec_opencl.d_TGV_MLEM_rhs);
			//		clSetKernelArg(kernel_ml, 62, sizeof(cl_uchar), &no_norm);
			//		status = clEnqueueNDRangeKernel(af_queue, kernel_ml, 1, NULL, &length[kk], NULL, 0, NULL, NULL);
			//		clFinish(af_queue);
			//	}

			//	if (status != CL_SUCCESS) {
			//		std::cerr << getErrorString(status) << std::endl;
			//		mexPrintf("Failed to launch the MLEM kernel\n");
			//		break;
			//	}
			//	//else if (verbose) {
			//	//	mexPrintf("OpenCL kernel executed successfully\n");
			//	//	mexEvalString("pause(.0001);");
			//	//}
			//	clFinish(af_queue);
			//	// Transfer memory control back to ArrayFire (ML-methods)
			//	unlock_AF_im_vectors(vec, MethodList, false, true, false, 0);
			//	if (compute_norm_matrix == 1u) {
			//		Summ.unlock();
			//		if (atomic_64bit)
			//			Summ = Summ.as(f32);
			//	}
			//	else if (compute_norm_matrix == 0u && osem_bool) {
			//		apu_sum.unlock();
			//		if (atomic_64bit) {
			//			apu_sum = apu_sum.as(f32);
			//		}
			//	}
			//	else {
			//		if (atomic_64bit) {
			//			apu_sum.unlock();
			//			Summ = apu_sum.as(f32);
			//		}
			//		else
			//			Summ.unlock();
			//	}

			//	// Compute the new estimates
			//	// MLEM
			//	if (compute_norm_matrix == 0u && osem_bool) {
			//		// Prevents division by zero
			//		apu_sum(apu_sum == 0.f) = epps;

			//		if (MethodList.MLEM)
			//			vec.MLEM_apu = MLEM(vec.MLEM_apu, apu_sum, vec.MLEM_rhs);
			//		// OSL-MLEM with Median Root Prior
			//		if (MethodList.OSLMLEM && MethodList.MRP) {
			//			array dU = MRP(vec.MRP_MLEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.med_no_norm, im_dim);
			//			vec.MRP_MLEM_apu = OSL_MLEM(vec.MRP_MLEM_apu, apu_sum, vec.MRP_MLEM_rhs, dU, beta.MRP_MLEM);
			//		}
			//		// OSL-MLEM with Quadratic prior
			//		if (MethodList.OSLMLEM && MethodList.Quad) {
			//			array dU = Quadratic_prior(vec.Quad_MLEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, w_vec.tr_offsets, w_vec.weights_quad, im_dim);
			//			vec.Quad_MLEM_apu = OSL_MLEM(vec.Quad_MLEM_apu, apu_sum, vec.Quad_MLEM_rhs, dU, beta.Quad_MLEM);
			//		}
			//		// OSL-MLEM with L-filter prior
			//		if (MethodList.OSLMLEM && MethodList.L) {
			//			array dU = L_filter(vec.L_MLEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L, w_vec.med_no_norm, im_dim);
			//			vec.L_MLEM_apu = OSL_MLEM(vec.L_MLEM_apu, apu_sum, vec.L_MLEM_rhs, dU, beta.L_MLEM);
			//		}
			//		// OSL-MLEM with FIR Median Hybrid prior
			//		if (MethodList.OSLMLEM && MethodList.FMH) {
			//			array dU = FMH(vec.FMH_MLEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets, w_vec.fmh_weights, w_vec.med_no_norm,
			//				w_vec.alku_fmh, im_dim);
			//			vec.FMH_MLEM_apu = OSL_MLEM(vec.FMH_MLEM_apu, apu_sum, vec.FMH_MLEM_rhs, dU, beta.FMH_MLEM);
			//		}
			//		// OSL-MLEM with Weighted Mean prior
			//		if (MethodList.OSLMLEM && MethodList.WeightedMean) {
			//			array dU = Weighted_mean(vec.Weighted_MLEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.weighted_weights,
			//				w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
			//			vec.Weighted_MLEM_apu = OSL_MLEM(vec.Weighted_MLEM_apu, apu_sum, vec.Weighted_MLEM_rhs, dU, beta.Weighted_MLEM);
			//		}
			//		// OSL-MLEM with Total Variation prior
			//		if (MethodList.OSLMLEM && MethodList.TV) {
			//			array dU = TVprior(Nx, Ny, Nz, data, vec.TV_MLEM_apu, epps, data.TVtype, w_vec);
			//			vec.TV_MLEM_apu = OSL_MLEM(vec.TV_MLEM_apu, apu_sum, vec.TV_MLEM_rhs, dU, beta.TV_MLEM);
			//		}
			//		// OSL-MLEM with Anisotropic Diffusion smoothing prior
			//		if (MethodList.OSLMLEM && MethodList.AD) {
			//			array dU = AD(vec.AD_MLEM_apu, Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, w_vec.DiffusionType, w_vec.med_no_norm);
			//			vec.AD_MLEM_apu = OSL_MLEM(vec.AD_MLEM_apu, apu_sum, vec.AD_MLEM_rhs, dU, beta.AD_MLEM);
			//		}
			//		// OSL-MLEM with Asymmetric Parallel Level Sets prior
			//		if (MethodList.OSLMLEM && MethodList.APLS) {
			//			array dU = TVprior(Nx, Ny, Nz, data, vec.TV_MLEM_apu, epps, 4, w_vec);
			//			vec.APLS_MLEM_apu = OSL_MLEM(vec.APLS_MLEM_apu, apu_sum, vec.APLS_MLEM_rhs, dU, beta.APLS_MLEM);
			//		}
			//		// OSL-MLEM with Total Generalized Variation prior
			//		if (MethodList.OSLMLEM && MethodList.TGV) {
			//			array dU = TGV(vec.TGV_MLEM_apu, Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
			//			vec.TGV_MLEM_apu = OSL_MLEM(vec.TGV_MLEM_apu, apu_sum, vec.TGV_MLEM_rhs, dU, beta.TGV_MLEM);
			//		}
			//	}
			//	else {
			//		// Prevents division by zero
			//		Summ(Summ == 0.f) = epps;

			//		if (MethodList.MLEM)
			//			vec.MLEM_apu = MLEM(vec.MLEM_apu, Summ, vec.MLEM_rhs);
			//		// OSL-MLEM with Median Root Prior
			//		if (MethodList.OSLMLEM && MethodList.MRP) {
			//			array dU = MRP(vec.MRP_MLEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.med_no_norm, im_dim);
			//			vec.MRP_MLEM_apu = OSL_MLEM(vec.MRP_MLEM_apu, Summ, vec.MRP_MLEM_rhs, dU, beta.MRP_MLEM);
			//		}
			//		// OSL-MLEM with Quadratic prior
			//		if (MethodList.OSLMLEM && MethodList.Quad) {
			//			array dU = Quadratic_prior(vec.Quad_MLEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, w_vec.tr_offsets, w_vec.weights_quad, im_dim);
			//			vec.Quad_MLEM_apu = OSL_MLEM(vec.Quad_MLEM_apu, Summ, vec.Quad_MLEM_rhs, dU, beta.Quad_MLEM);
			//		}
			//		// OSL-MLEM with L-filter prior
			//		if (MethodList.OSLMLEM && MethodList.L) {
			//			array dU = L_filter(vec.L_MLEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L, w_vec.med_no_norm, im_dim);
			//			vec.L_MLEM_apu = OSL_MLEM(vec.L_MLEM_apu, Summ, vec.L_MLEM_rhs, dU, beta.L_MLEM);
			//		}
			//		// OSL-MLEM with FIR Median Hybrid prior
			//		if (MethodList.OSLMLEM && MethodList.FMH) {
			//			array dU = FMH(vec.FMH_MLEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets, w_vec.fmh_weights, w_vec.med_no_norm,
			//				w_vec.alku_fmh, im_dim);
			//			vec.FMH_MLEM_apu = OSL_MLEM(vec.FMH_MLEM_apu, Summ, vec.FMH_MLEM_rhs, dU, beta.FMH_MLEM);
			//		}
			//		// OSL-MLEM with Weighted Mean prior
			//		if (MethodList.OSLMLEM && MethodList.WeightedMean) {
			//			array dU = Weighted_mean(vec.Weighted_MLEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.weighted_weights,
			//				w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
			//			vec.Weighted_MLEM_apu = OSL_MLEM(vec.Weighted_MLEM_apu, Summ, vec.Weighted_MLEM_rhs, dU, beta.Weighted_MLEM);
			//		}
			//		// OSL-MLEM with Total Variation prior
			//		if (MethodList.OSLMLEM && MethodList.TV) {
			//			array dU = TVprior(Nx, Ny, Nz, data, vec.TV_MLEM_apu, epps, data.TVtype, w_vec);
			//			vec.TV_MLEM_apu = OSL_MLEM(vec.TV_MLEM_apu, Summ, vec.TV_MLEM_rhs, dU, beta.TV_MLEM);
			//		}
			//		// OSL-MLEM with Anisotropic Diffusion smoothing prior
			//		if (MethodList.OSLMLEM && MethodList.AD) {
			//			array dU = AD(vec.AD_MLEM_apu, Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, w_vec.DiffusionType, w_vec.med_no_norm);
			//			vec.AD_MLEM_apu = OSL_MLEM(vec.AD_MLEM_apu, Summ, vec.AD_MLEM_rhs, dU, beta.AD_MLEM);
			//		}
			//		// OSL-MLEM with Asymmetric Parallel Level Sets prior
			//		if (MethodList.OSLMLEM && MethodList.APLS) {
			//			array dU = TVprior(Nx, Ny, Nz, data, vec.TV_MLEM_apu, epps, 4, w_vec);
			//			vec.APLS_MLEM_apu = OSL_MLEM(vec.APLS_MLEM_apu, Summ, vec.APLS_MLEM_rhs, dU, beta.APLS_MLEM);
			//		}
			//		// OSL-MLEM with Total Generalized Variation prior
			//		if (MethodList.OSLMLEM && MethodList.TGV) {
			//			array dU = TGV(vec.TGV_MLEM_apu, Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
			//			vec.TGV_MLEM_apu = OSL_MLEM(vec.TGV_MLEM_apu, Summ, vec.TGV_MLEM_rhs, dU, beta.TGV_MLEM);
			//		}
			//	}

			//	// Save the current iteration
			//	if (MethodList.MLEM)
			//		vec.MLEM(span, iter + 1u) = vec.MLEM_apu;
			//	if (MethodList.OSLMLEM) {
			//		if (MethodList.MRP)
			//			vec.MRP_MLEM(span, iter + 1u) = vec.MRP_MLEM_apu;
			//		if (MethodList.Quad)
			//			vec.Quad_MLEM(span, iter + 1u) = vec.Quad_MLEM_apu;
			//		if (MethodList.L)
			//			vec.L_MLEM(span, iter + 1u) = vec.L_MLEM_apu;
			//		if (MethodList.FMH)
			//			vec.FMH_MLEM(span, iter + 1u) = vec.FMH_MLEM_apu;
			//		if (MethodList.WeightedMean)
			//			vec.Weighted_MLEM(span, iter + 1u) = vec.Weighted_MLEM_apu;
			//		if (MethodList.TV)
			//			vec.TV_MLEM(span, iter + 1u) = vec.TV_MLEM_apu;
			//		if (MethodList.AD)
			//			vec.AD_MLEM(span, iter + 1u) = vec.AD_MLEM_apu;
			//		if (MethodList.APLS)
			//			vec.APLS_MLEM(span, iter + 1u) = vec.APLS_MLEM_apu;
			//		if (MethodList.TGV)
			//			vec.TGV_MLEM(span, iter + 1u) = vec.TGV_MLEM_apu;
			//	}

			//	if (verbose) {
			//		mexPrintf("MLEM iteration %d complete\n", iter + 1u);
			//		mexEvalString("pause(.0001);");
			//	}
			//}

			// Compute any of the other algorithms, if applicable
			if (osem_bool) {


				cl_mem* d_epsilon_mramla = w_vec.epsilon_mramla.device<cl_mem>();

				// Loop through the subsets
				for (uint32_t osa_iter = 0u; osa_iter < subsets; osa_iter++) {

					array uu;

					if (compute_norm_matrix == 1u) {
						if (atomic_64bit) {
							Summ[0] = constant(0ULL, im_dim, 1, u64);
						}
						else
							Summ[0] = constant(0.f, im_dim, 1);
						d_Summ = Summ[0].device<cl_mem>();
					}
					else {
						ll = osa_iter;
						if (atomic_64bit && iter == 0u) {
							Summ[ll] = Summ[ll].as(u64);
							d_Summ = Summ[ll].device<cl_mem>();
						}
						else {
							d_Summ = Summ[ll].device<cl_mem>();
						}
					}

					update_opencl_inputs(vec, vec_opencl, im_dim, MethodList, false, osa_iter, af_context, length[osa_iter], af_queue);


					// Set kernel arguments
					clSetKernelArg(kernel, 37, sizeof(cl_mem), d_Summ);
					clSetKernelArg(kernel, 38, sizeof(cl_mem), &d_lor[osa_iter]);
					clSetKernelArg(kernel, 39, sizeof(cl_mem), &d_xyindex[osa_iter]);
					clSetKernelArg(kernel, 49, sizeof(cl_mem), &d_zindex[osa_iter]);
					clSetKernelArg(kernel, 41, sizeof(cl_mem), &d_L[osa_iter]);
					clSetKernelArg(kernel, 42, sizeof(cl_mem), d_epsilon_mramla);
					clSetKernelArg(kernel, 43, sizeof(cl_mem), &d_Sino[osa_iter]);
					clSetKernelArg(kernel, 44, sizeof(cl_mem), vec_opencl.d_OSEM);
					clSetKernelArg(kernel, 45, sizeof(cl_mem), vec_opencl.d_RAMLA);
					clSetKernelArg(kernel, 46, sizeof(cl_mem), vec_opencl.d_MRAMLA);
					clSetKernelArg(kernel, 47, sizeof(cl_mem), vec_opencl.d_ROSEM);
					clSetKernelArg(kernel, 48, sizeof(cl_mem), vec_opencl.d_RBI);
					clSetKernelArg(kernel, 49, sizeof(cl_mem), vec_opencl.d_DRAMA);
					clSetKernelArg(kernel, 50, sizeof(cl_mem), vec_opencl.d_COSEM);
					clSetKernelArg(kernel, 51, sizeof(cl_mem), vec_opencl.d_ACOSEM);
					clSetKernelArg(kernel, 52, sizeof(cl_mem), vec_opencl.d_MRP_OSEM);
					clSetKernelArg(kernel, 53, sizeof(cl_mem), vec_opencl.d_Quad_OSEM);
					clSetKernelArg(kernel, 54, sizeof(cl_mem), vec_opencl.d_L_OSEM);
					clSetKernelArg(kernel, 55, sizeof(cl_mem), vec_opencl.d_FMH_OSEM);
					clSetKernelArg(kernel, 56, sizeof(cl_mem), vec_opencl.d_Weighted_OSEM);
					clSetKernelArg(kernel, 57, sizeof(cl_mem), vec_opencl.d_TV_OSEM);
					clSetKernelArg(kernel, 58, sizeof(cl_mem), vec_opencl.d_AD_OSEM);
					clSetKernelArg(kernel, 59, sizeof(cl_mem), vec_opencl.d_APLS_OSEM);
					clSetKernelArg(kernel, 60, sizeof(cl_mem), vec_opencl.d_TGV_OSEM);
					clSetKernelArg(kernel, 61, sizeof(cl_mem), vec_opencl.d_MRP_BSREM);
					clSetKernelArg(kernel, 62, sizeof(cl_mem), vec_opencl.d_Quad_BSREM);
					clSetKernelArg(kernel, 63, sizeof(cl_mem), vec_opencl.d_L_BSREM);
					clSetKernelArg(kernel, 64, sizeof(cl_mem), vec_opencl.d_FMH_BSREM);
					clSetKernelArg(kernel, 65, sizeof(cl_mem), vec_opencl.d_Weighted_BSREM);
					clSetKernelArg(kernel, 66, sizeof(cl_mem), vec_opencl.d_TV_BSREM);
					clSetKernelArg(kernel, 67, sizeof(cl_mem), vec_opencl.d_AD_BSREM);
					clSetKernelArg(kernel, 68, sizeof(cl_mem), vec_opencl.d_APLS_BSREM);
					clSetKernelArg(kernel, 69, sizeof(cl_mem), vec_opencl.d_TGV_BSREM);
					clSetKernelArg(kernel, 70, sizeof(cl_mem), vec_opencl.d_MRP_MBSREM);
					clSetKernelArg(kernel, 71, sizeof(cl_mem), vec_opencl.d_Quad_MBSREM);
					clSetKernelArg(kernel, 72, sizeof(cl_mem), vec_opencl.d_L_MBSREM);
					clSetKernelArg(kernel, 73, sizeof(cl_mem), vec_opencl.d_FMH_MBSREM);
					clSetKernelArg(kernel, 74, sizeof(cl_mem), vec_opencl.d_Weighted_MBSREM);
					clSetKernelArg(kernel, 75, sizeof(cl_mem), vec_opencl.d_TV_MBSREM);
					clSetKernelArg(kernel, 76, sizeof(cl_mem), vec_opencl.d_AD_MBSREM);
					clSetKernelArg(kernel, 77, sizeof(cl_mem), vec_opencl.d_APLS_MBSREM);
					clSetKernelArg(kernel, 78, sizeof(cl_mem), vec_opencl.d_TGV_MBSREM);
					clSetKernelArg(kernel, 79, sizeof(cl_mem), vec_opencl.d_MRP_ROSEM);
					clSetKernelArg(kernel, 80, sizeof(cl_mem), vec_opencl.d_Quad_ROSEM);
					clSetKernelArg(kernel, 81, sizeof(cl_mem), vec_opencl.d_L_ROSEM);
					clSetKernelArg(kernel, 82, sizeof(cl_mem), vec_opencl.d_FMH_ROSEM);
					clSetKernelArg(kernel, 83, sizeof(cl_mem), vec_opencl.d_Weighted_ROSEM);
					clSetKernelArg(kernel, 84, sizeof(cl_mem), vec_opencl.d_TV_ROSEM);
					clSetKernelArg(kernel, 85, sizeof(cl_mem), vec_opencl.d_AD_ROSEM);
					clSetKernelArg(kernel, 86, sizeof(cl_mem), vec_opencl.d_APLS_ROSEM);
					clSetKernelArg(kernel, 87, sizeof(cl_mem), vec_opencl.d_TGV_ROSEM);
					clSetKernelArg(kernel, 88, sizeof(cl_mem), vec_opencl.d_MRP_RBI);
					clSetKernelArg(kernel, 89, sizeof(cl_mem), vec_opencl.d_Quad_RBI);
					clSetKernelArg(kernel, 90, sizeof(cl_mem), vec_opencl.d_L_RBI);
					clSetKernelArg(kernel, 91, sizeof(cl_mem), vec_opencl.d_FMH_RBI);
					clSetKernelArg(kernel, 92, sizeof(cl_mem), vec_opencl.d_Weighted_RBI);
					clSetKernelArg(kernel, 93, sizeof(cl_mem), vec_opencl.d_TV_RBI);
					clSetKernelArg(kernel, 94, sizeof(cl_mem), vec_opencl.d_AD_RBI);
					clSetKernelArg(kernel, 95, sizeof(cl_mem), vec_opencl.d_APLS_RBI);
					clSetKernelArg(kernel, 96, sizeof(cl_mem), vec_opencl.d_TGV_RBI);
					clSetKernelArg(kernel, 97, sizeof(cl_mem), vec_opencl.d_MRP_COSEM);
					clSetKernelArg(kernel, 98, sizeof(cl_mem), vec_opencl.d_Quad_COSEM);
					clSetKernelArg(kernel, 99, sizeof(cl_mem), vec_opencl.d_L_COSEM);
					clSetKernelArg(kernel, 100, sizeof(cl_mem), vec_opencl.d_FMH_COSEM);
					clSetKernelArg(kernel, 101, sizeof(cl_mem), vec_opencl.d_Weighted_COSEM);
					clSetKernelArg(kernel, 102, sizeof(cl_mem), vec_opencl.d_TV_COSEM);
					clSetKernelArg(kernel, 103, sizeof(cl_mem), vec_opencl.d_AD_COSEM);
					clSetKernelArg(kernel, 104, sizeof(cl_mem), vec_opencl.d_APLS_COSEM);
					clSetKernelArg(kernel, 105, sizeof(cl_mem), vec_opencl.d_TGV_COSEM);
					clSetKernelArg(kernel, 106, sizeof(cl_mem), vec_opencl.d_OSEM_rhs);
					clSetKernelArg(kernel, 107, sizeof(cl_mem), vec_opencl.d_RAMLA_rhs);
					clSetKernelArg(kernel, 108, sizeof(cl_mem), vec_opencl.d_MRAMLA_rhs);
					clSetKernelArg(kernel, 109, sizeof(cl_mem), vec_opencl.d_ROSEM_rhs);
					clSetKernelArg(kernel, 110, sizeof(cl_mem), vec_opencl.d_RBI_rhs);
					clSetKernelArg(kernel, 111, sizeof(cl_mem), vec_opencl.d_DRAMA_rhs);
					clSetKernelArg(kernel, 112, sizeof(cl_mem), vec_opencl.d_COSEM_rhs);
					clSetKernelArg(kernel, 113, sizeof(cl_mem), vec_opencl.d_ACOSEM_rhs);
					clSetKernelArg(kernel, 114, sizeof(cl_mem), vec_opencl.d_MRP_OSEM_rhs);
					clSetKernelArg(kernel, 115, sizeof(cl_mem), vec_opencl.d_Quad_OSEM_rhs);
					clSetKernelArg(kernel, 116, sizeof(cl_mem), vec_opencl.d_L_OSEM_rhs);
					clSetKernelArg(kernel, 117, sizeof(cl_mem), vec_opencl.d_FMH_OSEM_rhs);
					clSetKernelArg(kernel, 118, sizeof(cl_mem), vec_opencl.d_Weighted_OSEM_rhs);
					clSetKernelArg(kernel, 119, sizeof(cl_mem), vec_opencl.d_TV_OSEM_rhs);
					clSetKernelArg(kernel, 120, sizeof(cl_mem), vec_opencl.d_AD_OSEM_rhs);
					clSetKernelArg(kernel, 121, sizeof(cl_mem), vec_opencl.d_APLS_OSEM_rhs);
					clSetKernelArg(kernel, 122, sizeof(cl_mem), vec_opencl.d_TGV_OSEM_rhs);
					clSetKernelArg(kernel, 123, sizeof(cl_mem), vec_opencl.d_MRP_BSREM_rhs);
					clSetKernelArg(kernel, 124, sizeof(cl_mem), vec_opencl.d_Quad_BSREM_rhs);
					clSetKernelArg(kernel, 125, sizeof(cl_mem), vec_opencl.d_L_BSREM_rhs);
					clSetKernelArg(kernel, 126, sizeof(cl_mem), vec_opencl.d_FMH_BSREM_rhs);
					clSetKernelArg(kernel, 127, sizeof(cl_mem), vec_opencl.d_Weighted_BSREM_rhs);
					clSetKernelArg(kernel, 128, sizeof(cl_mem), vec_opencl.d_TV_BSREM_rhs);
					clSetKernelArg(kernel, 129, sizeof(cl_mem), vec_opencl.d_AD_BSREM_rhs);
					clSetKernelArg(kernel, 130, sizeof(cl_mem), vec_opencl.d_APLS_BSREM_rhs);
					clSetKernelArg(kernel, 131, sizeof(cl_mem), vec_opencl.d_TGV_BSREM_rhs);
					clSetKernelArg(kernel, 132, sizeof(cl_mem), vec_opencl.d_MRP_MBSREM_rhs);
					clSetKernelArg(kernel, 133, sizeof(cl_mem), vec_opencl.d_Quad_MBSREM_rhs);
					clSetKernelArg(kernel, 134, sizeof(cl_mem), vec_opencl.d_L_MBSREM_rhs);
					clSetKernelArg(kernel, 135, sizeof(cl_mem), vec_opencl.d_FMH_MBSREM_rhs);
					clSetKernelArg(kernel, 136, sizeof(cl_mem), vec_opencl.d_Weighted_MBSREM_rhs);
					clSetKernelArg(kernel, 137, sizeof(cl_mem), vec_opencl.d_TV_MBSREM_rhs);
					clSetKernelArg(kernel, 138, sizeof(cl_mem), vec_opencl.d_AD_MBSREM_rhs);
					clSetKernelArg(kernel, 139, sizeof(cl_mem), vec_opencl.d_APLS_MBSREM_rhs);
					clSetKernelArg(kernel, 140, sizeof(cl_mem), vec_opencl.d_TGV_MBSREM_rhs);
					clSetKernelArg(kernel, 141, sizeof(cl_mem), vec_opencl.d_MRP_ROSEM_rhs);
					clSetKernelArg(kernel, 142, sizeof(cl_mem), vec_opencl.d_Quad_ROSEM_rhs);
					clSetKernelArg(kernel, 143, sizeof(cl_mem), vec_opencl.d_L_ROSEM_rhs);
					clSetKernelArg(kernel, 144, sizeof(cl_mem), vec_opencl.d_FMH_ROSEM_rhs);
					clSetKernelArg(kernel, 145, sizeof(cl_mem), vec_opencl.d_Weighted_ROSEM_rhs);
					clSetKernelArg(kernel, 146, sizeof(cl_mem), vec_opencl.d_TV_ROSEM_rhs);
					clSetKernelArg(kernel, 147, sizeof(cl_mem), vec_opencl.d_AD_ROSEM_rhs);
					clSetKernelArg(kernel, 148, sizeof(cl_mem), vec_opencl.d_APLS_ROSEM_rhs);
					clSetKernelArg(kernel, 149, sizeof(cl_mem), vec_opencl.d_TGV_ROSEM_rhs);
					clSetKernelArg(kernel, 150, sizeof(cl_mem), vec_opencl.d_MRP_RBI_rhs);
					clSetKernelArg(kernel, 151, sizeof(cl_mem), vec_opencl.d_Quad_RBI_rhs);
					clSetKernelArg(kernel, 152, sizeof(cl_mem), vec_opencl.d_L_RBI_rhs);
					clSetKernelArg(kernel, 153, sizeof(cl_mem), vec_opencl.d_FMH_RBI_rhs);
					clSetKernelArg(kernel, 154, sizeof(cl_mem), vec_opencl.d_Weighted_RBI_rhs);
					clSetKernelArg(kernel, 155, sizeof(cl_mem), vec_opencl.d_TV_RBI_rhs);
					clSetKernelArg(kernel, 156, sizeof(cl_mem), vec_opencl.d_AD_RBI_rhs);
					clSetKernelArg(kernel, 157, sizeof(cl_mem), vec_opencl.d_APLS_RBI_rhs);
					clSetKernelArg(kernel, 158, sizeof(cl_mem), vec_opencl.d_TGV_RBI_rhs);
					clSetKernelArg(kernel, 159, sizeof(cl_mem), vec_opencl.d_MRP_COSEM_rhs);
					clSetKernelArg(kernel, 160, sizeof(cl_mem), vec_opencl.d_Quad_COSEM_rhs);
					clSetKernelArg(kernel, 161, sizeof(cl_mem), vec_opencl.d_L_COSEM_rhs);
					clSetKernelArg(kernel, 162, sizeof(cl_mem), vec_opencl.d_FMH_COSEM_rhs);
					clSetKernelArg(kernel, 163, sizeof(cl_mem), vec_opencl.d_Weighted_COSEM_rhs);
					clSetKernelArg(kernel, 164, sizeof(cl_mem), vec_opencl.d_TV_COSEM_rhs);
					clSetKernelArg(kernel, 165, sizeof(cl_mem), vec_opencl.d_AD_COSEM_rhs);
					clSetKernelArg(kernel, 166, sizeof(cl_mem), vec_opencl.d_APLS_COSEM_rhs);
					clSetKernelArg(kernel, 167, sizeof(cl_mem), vec_opencl.d_TGV_COSEM_rhs);
					clSetKernelArg(kernel, 168, sizeof(cl_uchar), &no_norm);
					//clSetKernelArg(kernel, 161, sizeof(float) * 62, NULL);
					//clSetKernelArg(kernel, 161, sizeof(cl_int)*max_lor, NULL);
					//clSetKernelArg(kernel, 162, sizeof(cl_float)*max_lor, NULL);
					// Compute the kernel
					status = clEnqueueNDRangeKernel(af_queue, kernel, 1u, NULL, &length[osa_iter], NULL, 0, NULL, NULL);

					if (status != CL_SUCCESS) {
						std::cerr << getErrorString(status) << std::endl;
						mexPrintf("Failed to launch the OS kernel\n");
						break;
					}
					//else if (verbose) {
					//	mexPrintf("OS kernel launched successfully\n");
					//	mexEvalString("pause(.0001);");
					//}

					if (MethodList.ACOSEM || MethodList.OSLCOSEM == 1u) {
						uu = afcl::array(length[osa_iter], d_Sino[osa_iter], f32, true);
						uu = sum(uu);
						//uu = sum(uu.as(f32));
					}
					clFinish(af_queue);
					// Transfer memory control back to ArrayFire (OS-methods)

					if (compute_norm_matrix == 1u) {
						Summ[0].unlock();
						if (atomic_64bit)
							Summ[0] = Summ[0].as(f32) / TH;
					}
					else {
						Summ[ll].unlock();
						if (atomic_64bit && iter == 0u) {
							Summ[ll] = Summ[ll].as(f32) / TH;
						}
					}
					unlock_AF_im_vectors(vec, MethodList, false, false, true, osa_iter);

					// Compute the (matrix free) algorithms
						// Prevent division by zero
					Summ[ll](Summ[ll] == 0.f) = epps;
					// Ordered Subsets Expectation Maximization (OSEM)
					if (MethodList.OSEM || MethodList.ECOSEM) {
						//vec.OSEM_rhs(vec.OSEM_rhs == 0.f) = epps;
						vec.OSEM_apu = OSEM(vec.OSEM_apu, Summ[ll], vec.OSEM_rhs);
					}
					// OSL-OSEM
					if (MethodList.OSLOSEM) {
						// Median Root Prior
						if (MethodList.MRP) {
							array dU = MRP(vec.MRP_OSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.med_no_norm, im_dim);
							vec.MRP_OSEM_apu = OSL_OSEM(vec.MRP_OSEM_apu, Summ[ll], vec.MRP_OSEM_rhs, dU, beta.MRP_OSEM);
						}
						// Quadratic Prior
						if (MethodList.Quad) {
							array dU = Quadratic_prior(vec.Quad_OSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, w_vec.tr_offsets, w_vec.weights_quad, im_dim);
							vec.Quad_OSEM_apu = OSL_OSEM(vec.Quad_OSEM_apu, Summ[ll], vec.Quad_OSEM_rhs, dU, beta.Quad_OSEM);
						}
						// L-filter prior
						if (MethodList.L) {
							array dU = L_filter(vec.L_OSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L, w_vec.med_no_norm, im_dim);
							vec.L_OSEM_apu = OSL_OSEM(vec.L_OSEM_apu, Summ[ll], vec.L_OSEM_rhs, dU, beta.L_OSEM);
						}
						// FIR Median Hybrid prior
						if (MethodList.FMH) {
							array dU = FMH(vec.FMH_OSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets, w_vec.fmh_weights, w_vec.med_no_norm,
								w_vec.alku_fmh, im_dim);
							vec.FMH_OSEM_apu = OSL_OSEM(vec.FMH_OSEM_apu, Summ[ll], vec.FMH_OSEM_rhs, dU, beta.FMH_OSEM);
						}
						// Weighted Mean prior
						if (MethodList.WeightedMean) {
							array dU = Weighted_mean(vec.Weighted_OSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.weighted_weights,
								w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
							vec.Weighted_OSEM_apu = OSL_OSEM(vec.Weighted_OSEM_apu, Summ[ll], vec.Weighted_OSEM_rhs, dU, beta.Weighted_OSEM);
						}
						// Total Variation prior
						if (MethodList.TV) {
							array dU = TVprior(Nx, Ny, Nz, data, vec.TV_OSEM_apu, epps, data.TVtype, w_vec);
							vec.TV_OSEM_apu = OSL_OSEM(vec.TV_OSEM_apu, Summ[ll], vec.TV_OSEM_rhs, dU, beta.TV_OSEM);
						}
						// Anisotropic Diffusion smoothing prior
						if (MethodList.AD) {
							if (osa_iter == 0u)
								if (MethodList.OSEM)
									vec.AD_OSEM_apu = vec.OSEM_apu;
								else
									vec.AD_OSEM_apu = OSEM(vec.AD_OSEM_apu, Summ[ll], vec.AD_OSEM_rhs);
							else {
								array dU = AD(vec.AD_OSEM_apu, Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, w_vec.DiffusionType, w_vec.med_no_norm);
								vec.AD_OSEM_apu = OSL_OSEM(vec.AD_OSEM_apu, Summ[ll], vec.AD_OSEM_rhs, dU, beta.AD_OSEM);
							}
						}
						// Asymmetric Parallel Level Sets prior
						if (MethodList.APLS) {
							array dU = TVprior(Nx, Ny, Nz, data, vec.TV_OSEM_apu, epps, 4, w_vec);
							vec.APLS_OSEM_apu = OSL_OSEM(vec.APLS_OSEM_apu, Summ[ll], vec.APLS_OSEM_rhs, dU, beta.APLS_OSEM);
						}
						// Total Generalized Variation prior
						if (MethodList.TGV) {
							array dU = TGV(vec.TGV_OSEM_apu, Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
							vec.TGV_OSEM_apu = OSL_OSEM(vec.TGV_OSEM_apu, Summ[ll], vec.TGV_OSEM_rhs, dU, beta.TGV_OSEM);
						}
					}

					// Modfied Row-action Maximum Likelihood (MRAMLA)
					if (MethodList.MRAMLA)
						vec.MRAMLA_apu = MBSREM(vec.MRAMLA_apu, vec.MRAMLA_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM, iter, im_dim, 0.f, af::constant(0.f, 1, 1),
							Summ[ll], epps);
					// Modfied Block Sequential Regularized Expectation Maximization (MBSREM)
					if (MethodList.MBSREM) {
						if (MethodList.MRP) {
							array dU = MRP(vec.MRP_MBSREM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.med_no_norm, im_dim);
							vec.MRP_MBSREM_apu = MBSREM(vec.MRP_MBSREM_apu, vec.MRP_MBSREM_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.MRP_MBSREM,
								dU, Summ[ll], epps);
						}
						if (MethodList.Quad) {
							array dU = Quadratic_prior(vec.Quad_MBSREM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, w_vec.tr_offsets, w_vec.weights_quad, im_dim);
							vec.Quad_MBSREM_apu = MBSREM(vec.Quad_MBSREM_apu, vec.MRAMLA_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM,
								iter, im_dim, beta.Quad_MBSREM, dU, Summ[ll], epps);
						}
						if (MethodList.L) {
							array dU = L_filter(vec.L_MBSREM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L, w_vec.med_no_norm, im_dim);
							vec.L_MBSREM_apu = MBSREM(vec.L_MBSREM_apu, vec.L_MBSREM_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM,
								iter, im_dim, beta.L_MBSREM, dU, Summ[ll], epps);
						}
						if (MethodList.FMH) {
							array dU = FMH(vec.FMH_MBSREM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets, w_vec.fmh_weights, w_vec.med_no_norm,
								w_vec.alku_fmh, im_dim);
							vec.FMH_MBSREM_apu = MBSREM(vec.FMH_MBSREM_apu, vec.FMH_MBSREM_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM,
								iter, im_dim, beta.FMH_MBSREM, dU, Summ[ll], epps);
						}
						if (MethodList.WeightedMean) {
							array dU = Weighted_mean(vec.Weighted_MBSREM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.weighted_weights,
								w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
							vec.Weighted_MBSREM_apu = MBSREM(vec.Weighted_MBSREM_apu, vec.Weighted_MBSREM_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM,
								iter, im_dim, beta.Weighted_MBSREM, dU, Summ[ll], epps);
						}
						if (MethodList.TV) {
							array dU = TVprior(Nx, Ny, Nz, data, vec.TV_MBSREM_apu, epps, data.TVtype, w_vec);
							vec.TV_MBSREM_apu = MBSREM(vec.TV_MBSREM_apu, vec.TV_MBSREM_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM,
								iter, im_dim, beta.TV_MBSREM, dU, Summ[ll], epps);
						}
						if (MethodList.AD) {
							array dU = AD(vec.AD_MBSREM_apu, Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, w_vec.DiffusionType, w_vec.med_no_norm);
							vec.AD_MBSREM_apu = MBSREM(vec.AD_MBSREM_apu, vec.AD_MBSREM_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM,
								iter, im_dim, beta.AD_MBSREM, dU, Summ[ll], epps);
						}
						if (MethodList.APLS) {
							array dU = TVprior(Nx, Ny, Nz, data, vec.TV_MBSREM_apu, epps, 4, w_vec);
							vec.TV_MBSREM_apu = MBSREM(vec.TV_MBSREM_apu, vec.TV_MBSREM_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM,
								iter, im_dim, beta.APLS_MBSREM, dU, Summ[ll], epps);
						}
						if (MethodList.TGV) {
							array dU = TGV(vec.TGV_MBSREM_apu, Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
							vec.TGV_MBSREM_apu = MBSREM(vec.TGV_MBSREM_apu, vec.TGV_MBSREM_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM,
								iter, im_dim, beta.TGV_MBSREM, dU, Summ[ll], epps);
						}
					}

					// Row-action Maximum Likelihood (RAMLA)
					if (MethodList.RAMLA)
						vec.RAMLA_apu = BSREM(vec.RAMLA_apu, vec.RAMLA_rhs, w_vec.lambda_BSREM, iter);
					// Block Sequential Regularized Expectation Maximization (BSREM)
					if (MethodList.BSREM) {
						if (MethodList.MRP) {
							vec.MRP_BSREM_apu = BSREM(vec.MRP_BSREM_apu, vec.MRP_BSREM_rhs, w_vec.lambda_BSREM, iter);
						}
						if (MethodList.Quad) {
							vec.Quad_BSREM_apu = BSREM(vec.Quad_BSREM_apu, vec.MRAMLA_rhs, w_vec.lambda_BSREM, iter);
						}
						if (MethodList.L) {
							vec.L_BSREM_apu = BSREM(vec.L_BSREM_apu, vec.L_BSREM_rhs, w_vec.lambda_BSREM, iter);
						}
						if (MethodList.FMH) {
							vec.FMH_BSREM_apu = BSREM(vec.FMH_BSREM_apu, vec.FMH_BSREM_rhs, w_vec.lambda_BSREM, iter);
						}
						if (MethodList.WeightedMean) {
							vec.Weighted_BSREM_apu = BSREM(vec.Weighted_BSREM_apu, vec.Weighted_BSREM_rhs, w_vec.lambda_BSREM, iter);
						}
						if (MethodList.TV) {
							vec.TV_BSREM_apu = BSREM(vec.TV_BSREM_apu, vec.TV_BSREM_rhs, w_vec.lambda_BSREM, iter);
						}
						if (MethodList.AD) {
							vec.AD_BSREM_apu = BSREM(vec.AD_BSREM_apu, vec.AD_BSREM_rhs, w_vec.lambda_BSREM, iter);
						}
						if (MethodList.APLS) {
							vec.TV_BSREM_apu = BSREM(vec.TV_BSREM_apu, vec.TV_BSREM_rhs, w_vec.lambda_BSREM, iter);
						}
						if (MethodList.TGV) {
							vec.TGV_BSREM_apu = BSREM(vec.TGV_BSREM_apu, vec.TGV_BSREM_rhs, w_vec.lambda_BSREM, iter);
						}
					}

					// Relaxed OSEM (ROSEM)
					if (MethodList.ROSEM)
						vec.ROSEM_apu = ROSEM(vec.ROSEM_apu, Summ[ll], vec.ROSEM_rhs, w_vec.lambda_ROSEM, iter);
					if (MethodList.ROSEMMAP) {
						if (MethodList.MRP) {
							array dU = MRP(vec.MRP_ROSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.med_no_norm, im_dim);
							vec.MRP_ROSEM_apu = ROSEM(vec.MRP_ROSEM_apu, Summ[ll], vec.MRP_ROSEM_rhs, w_vec.lambda_ROSEM, iter);
						}
						if (MethodList.Quad) {
							array dU = Quadratic_prior(vec.Quad_ROSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, w_vec.tr_offsets, w_vec.weights_quad, im_dim);
							vec.Quad_ROSEM_apu = ROSEM(vec.Quad_ROSEM_apu, Summ[ll], vec.ROSEM_rhs, w_vec.lambda_ROSEM, iter);
						}
						if (MethodList.L) {
							array dU = L_filter(vec.L_ROSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L, w_vec.med_no_norm, im_dim);
							vec.L_ROSEM_apu = ROSEM(vec.L_ROSEM_apu, Summ[ll], vec.L_ROSEM_rhs, w_vec.lambda_ROSEM, iter);
						}
						if (MethodList.FMH) {
							array dU = FMH(vec.FMH_ROSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets, w_vec.fmh_weights, w_vec.med_no_norm,
								w_vec.alku_fmh, im_dim);
							vec.FMH_ROSEM_apu = ROSEM(vec.FMH_ROSEM_apu, Summ[ll], vec.FMH_ROSEM_rhs, w_vec.lambda_ROSEM, iter);
						}
						if (MethodList.WeightedMean) {
							array dU = Weighted_mean(vec.Weighted_ROSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.weighted_weights,
								w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
							vec.Weighted_ROSEM_apu = ROSEM(vec.Weighted_ROSEM_apu, Summ[ll], vec.Weighted_ROSEM_rhs, w_vec.lambda_ROSEM, iter);
						}
						if (MethodList.TV) {
							array dU = TVprior(Nx, Ny, Nz, data, vec.TV_ROSEM_apu, epps, data.TVtype, w_vec);
							vec.TV_ROSEM_apu = ROSEM(vec.TV_ROSEM_apu, Summ[ll], vec.TV_ROSEM_rhs, w_vec.lambda_ROSEM, iter);
						}
						if (MethodList.AD) {
							array dU = AD(vec.AD_ROSEM_apu, Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, w_vec.DiffusionType, w_vec.med_no_norm);
							vec.AD_ROSEM_apu = ROSEM(vec.AD_ROSEM_apu, Summ[ll], vec.AD_ROSEM_rhs, w_vec.lambda_ROSEM, iter);
						}
						if (MethodList.APLS) {
							array dU = TVprior(Nx, Ny, Nz, data, vec.TV_ROSEM_apu, epps, 4, w_vec);
							vec.APLS_ROSEM_apu = ROSEM(vec.TV_ROSEM_apu, Summ[ll], vec.TV_ROSEM_rhs, w_vec.lambda_ROSEM, iter);
						}
						if (MethodList.TGV) {
							array dU = TGV(vec.TGV_ROSEM_apu, Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
							vec.TGV_ROSEM_apu = ROSEM(vec.TGV_ROSEM_apu, Summ[ll], vec.TGV_ROSEM_rhs, w_vec.lambda_ROSEM, iter);
						}
					}

					// Rescaled Block Iterative EM (RBI)
					if (MethodList.RBI)
						vec.RBI_apu = RBI(vec.RBI_apu, Summ[ll], vec.RBI_rhs, w_vec.D, 0.f, af::constant(0.f, 1, 1));
					if (MethodList.RBIMAP) {
						if (MethodList.MRP) {
							array dU = MRP(vec.MRP_RBI_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.med_no_norm, im_dim);
							vec.MRP_RBI_apu = RBI(vec.MRP_RBI_apu, Summ[ll], vec.MRP_RBI_rhs, w_vec.D, beta.MRP_RBI, dU);
						}
						if (MethodList.Quad) {
							array dU = Quadratic_prior(vec.Quad_RBI_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, w_vec.tr_offsets, w_vec.weights_quad, im_dim);
							vec.Quad_RBI_apu = RBI(vec.Quad_RBI_apu, Summ[ll], vec.RBI_rhs, w_vec.D, beta.Quad_RBI, dU);
						}
						if (MethodList.L) {
							array dU = L_filter(vec.L_RBI_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L, w_vec.med_no_norm, im_dim);
							vec.L_RBI_apu = RBI(vec.L_RBI_apu, Summ[ll], vec.L_RBI_rhs, w_vec.D, beta.L_RBI, dU);
						}
						if (MethodList.FMH) {
							array dU = FMH(vec.FMH_RBI_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets, w_vec.fmh_weights, w_vec.med_no_norm,
								w_vec.alku_fmh, im_dim);
							vec.FMH_RBI_apu = RBI(vec.FMH_RBI_apu, Summ[ll], vec.FMH_RBI_rhs, w_vec.D, beta.FMH_RBI, dU);
						}
						if (MethodList.WeightedMean) {
							array dU = Weighted_mean(vec.Weighted_RBI_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.weighted_weights,
								w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
							vec.Weighted_RBI_apu = RBI(vec.Weighted_RBI_apu, Summ[ll], vec.Weighted_RBI_rhs, w_vec.D, beta.Weighted_RBI, dU);
						}
						if (MethodList.TV) {
							array dU = TVprior(Nx, Ny, Nz, data, vec.TV_RBI_apu, epps, data.TVtype, w_vec);
							vec.TV_RBI_apu = RBI(vec.TV_RBI_apu, Summ[ll], vec.TV_RBI_rhs, w_vec.D, beta.TV_RBI, dU);
						}
						if (MethodList.AD) {
							array dU = AD(vec.AD_RBI_apu, Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, w_vec.DiffusionType, w_vec.med_no_norm);
							vec.AD_RBI_apu = RBI(vec.AD_RBI_apu, Summ[ll], vec.AD_RBI_rhs, w_vec.D, beta.AD_RBI, dU);
						}
						if (MethodList.APLS) {
							array dU = TVprior(Nx, Ny, Nz, data, vec.TV_RBI_apu, epps, 4, w_vec);
							vec.APLS_RBI_apu = RBI(vec.TV_RBI_apu, Summ[ll], vec.TV_RBI_rhs, w_vec.D, beta.APLS_RBI, dU);
						}
						if (MethodList.TGV) {
							array dU = TGV(vec.TGV_RBI_apu, Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
							vec.TGV_RBI_apu = RBI(vec.TGV_RBI_apu, Summ[ll], vec.TGV_RBI_rhs, w_vec.D, beta.TGV_RBI, dU);
						}
					}

					// Complete data OSEM
					if (MethodList.COSEM || MethodList.ECOSEM)
						vec.COSEM_apu = COSEM(vec.COSEM_apu, vec.C_co, w_vec.D);

					// Enhanced COSEM
					if (MethodList.ECOSEM)
						vec.ECOSEM_apu = ECOSEM(vec.ECOSEM_apu, w_vec.D, vec.OSEM_apu, vec.COSEM_apu, epps);

					// Accelerated COSEM
					if (MethodList.ACOSEM) {
						vec.ACOSEM_apu = ACOSEM(vec.ACOSEM_apu, vec.C_aco, w_vec.D, w_vec.h_ACOSEM);
						MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.ACOSEM_apu, vec.C_co,
							vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length, atomic_64bit, compute_norm_matrix);
						w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
						vec.ACOSEM_apu = batchFunc(vec.ACOSEM_apu, uu / w_vec.ACOSEM_rhs, batchMul);
						vec.ACOSEM_apu(vec.ACOSEM_apu < 0.f) = epps;
					}

					if (MethodList.OSLCOSEM > 0u) {
						if (MethodList.MRP) {
							array dU = MRP(vec.MRP_COSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.med_no_norm, im_dim);
							vec.MRP_COSEM_apu = OSL_COSEM(vec.MRP_COSEM_apu, vec.C_osl, w_vec.D, w_vec.h_ACOSEM, MethodList.OSLCOSEM, dU, beta.MRP_COSEM);
							if (MethodList.OSLCOSEM == 1u) {
								MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.MRP_COSEM_apu, vec.C_co,
									vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length, atomic_64bit, compute_norm_matrix);
								w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
								vec.MRP_COSEM_apu = batchFunc(vec.MRP_COSEM_apu, uu / w_vec.ACOSEM_rhs, batchMul);
								vec.MRP_COSEM_apu(vec.MRP_COSEM_apu < 0.f) = epps;
							}
						}
						if (MethodList.Quad) {
							array dU = Quadratic_prior(vec.Quad_COSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, w_vec.tr_offsets, w_vec.weights_quad, im_dim);
							vec.Quad_COSEM_apu = OSL_COSEM(vec.Quad_COSEM_apu, vec.C_osl, w_vec.D, w_vec.h_ACOSEM, MethodList.OSLCOSEM, dU, beta.Quad_COSEM);
							if (MethodList.OSLCOSEM == 1u) {
								MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.Quad_COSEM_apu, vec.C_co,
									vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length, atomic_64bit, compute_norm_matrix);
								w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
								vec.Quad_COSEM_apu = batchFunc(vec.Quad_COSEM_apu, uu / w_vec.ACOSEM_rhs, batchMul);
								vec.Quad_COSEM_apu(vec.Quad_COSEM_apu < 0.f) = epps;
							}
						}
						if (MethodList.L) {
							array dU = L_filter(vec.L_COSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L, w_vec.med_no_norm, im_dim);
							vec.L_COSEM_apu = OSL_COSEM(vec.L_COSEM_apu, vec.C_osl, w_vec.D, w_vec.h_ACOSEM, MethodList.OSLCOSEM, dU, beta.L_COSEM);
							if (MethodList.OSLCOSEM == 1u) {
								MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.L_COSEM_apu, vec.C_co,
									vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length, atomic_64bit, compute_norm_matrix);
								w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
								vec.L_COSEM_apu = batchFunc(vec.L_COSEM_apu, uu / w_vec.ACOSEM_rhs, batchMul);
								vec.L_COSEM_apu(vec.L_COSEM_apu < 0.f) = epps;
							}
						}
						if (MethodList.FMH) {
							array dU = FMH(vec.FMH_COSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets, w_vec.fmh_weights, w_vec.med_no_norm,
								w_vec.alku_fmh, im_dim);
							vec.FMH_COSEM_apu = OSL_COSEM(vec.FMH_COSEM_apu, vec.C_osl, w_vec.D, w_vec.h_ACOSEM, MethodList.OSLCOSEM, dU, beta.FMH_COSEM);
							if (MethodList.OSLCOSEM == 1u) {
								MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.FMH_COSEM_apu, vec.C_co,
									vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length, atomic_64bit, compute_norm_matrix);
								w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
								vec.FMH_COSEM_apu = batchFunc(vec.FMH_COSEM_apu, uu / w_vec.ACOSEM_rhs, batchMul);
								vec.FMH_COSEM_apu(vec.FMH_COSEM_apu < 0.f) = epps;
							}
						}
						if (MethodList.WeightedMean) {
							array dU = Weighted_mean(vec.Weighted_COSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.weighted_weights,
								w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
							vec.Weighted_COSEM_apu = OSL_COSEM(vec.Weighted_COSEM_apu, vec.C_osl, w_vec.D, w_vec.h_ACOSEM, MethodList.OSLCOSEM, dU, beta.Weighted_COSEM);
							if (MethodList.OSLCOSEM == 1u) {
								MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.Weighted_COSEM_apu, vec.C_co,
									vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length, atomic_64bit, compute_norm_matrix);
								w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
								vec.Weighted_COSEM_apu = batchFunc(vec.Weighted_COSEM_apu, uu / w_vec.ACOSEM_rhs, batchMul);
								vec.Weighted_COSEM_apu(vec.Weighted_COSEM_apu < 0.f) = epps;
							}
						}
						if (MethodList.TV) {
							array dU = TVprior(Nx, Ny, Nz, data, vec.TV_COSEM_apu, epps, data.TVtype, w_vec);
							vec.TV_COSEM_apu = OSL_COSEM(vec.TV_COSEM_apu, vec.C_osl, w_vec.D, w_vec.h_ACOSEM, MethodList.OSLCOSEM, dU, beta.TV_COSEM);
							if (MethodList.OSLCOSEM == 1u) {
								MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.TV_COSEM_apu, vec.C_co,
									vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length, atomic_64bit, compute_norm_matrix);
								w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
								vec.TV_COSEM_apu = batchFunc(vec.TV_COSEM_apu, uu / w_vec.ACOSEM_rhs, batchMul);
								vec.TV_COSEM_apu(vec.TV_COSEM_apu < 0.f) = epps;
							}
						}
						if (MethodList.AD) {
							array dU = AD(vec.AD_COSEM_apu, Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, w_vec.DiffusionType, w_vec.med_no_norm);
							vec.AD_COSEM_apu = OSL_COSEM(vec.AD_COSEM_apu, vec.C_osl, w_vec.D, w_vec.h_ACOSEM, MethodList.OSLCOSEM, dU, beta.AD_COSEM);
							if (MethodList.OSLCOSEM == 1u) {
								MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.AD_COSEM_apu, vec.C_co,
									vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length, atomic_64bit, compute_norm_matrix);
								w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
								vec.AD_COSEM_apu = batchFunc(vec.AD_COSEM_apu, uu / w_vec.ACOSEM_rhs, batchMul);
								vec.AD_COSEM_apu(vec.AD_COSEM_apu < 0.f) = epps;
							}
						}
						if (MethodList.APLS) {
							array dU = TVprior(Nx, Ny, Nz, data, vec.TV_COSEM_apu, epps, 4, w_vec);
							vec.APLS_COSEM_apu = OSL_COSEM(vec.APLS_COSEM_apu, vec.C_osl, w_vec.D, w_vec.h_ACOSEM, MethodList.OSLCOSEM, dU, beta.APLS_COSEM);
							if (MethodList.OSLCOSEM == 1u) {
								MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.APLS_COSEM_apu, vec.C_co,
									vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length, atomic_64bit, compute_norm_matrix);
								w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
								vec.APLS_COSEM_apu = batchFunc(vec.APLS_COSEM_apu, uu / w_vec.ACOSEM_rhs, batchMul);
								vec.APLS_COSEM_apu(vec.APLS_COSEM_apu < 0.f) = epps;
							}
						}
						if (MethodList.TGV) {
							array dU = TGV(vec.TGV_COSEM_apu, Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
							vec.TGV_COSEM_apu = OSL_COSEM(vec.TGV_COSEM_apu, vec.C_osl, w_vec.D, w_vec.h_ACOSEM, MethodList.OSLCOSEM, dU, beta.TGV_COSEM);
							if (MethodList.OSLCOSEM == 1u) {
								MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.TGV_COSEM_apu, vec.C_co,
									vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length, atomic_64bit, compute_norm_matrix);
								w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
								vec.TGV_COSEM_apu = batchFunc(vec.TGV_COSEM_apu, uu / w_vec.ACOSEM_rhs, batchMul);
								vec.TGV_COSEM_apu(vec.TGV_COSEM_apu < 0.f) = epps;
							}
						}
					}

					// Dynamic RAMLA
					if (MethodList.DRAMA)
						vec.DRAMA_apu = DRAMA(vec.DRAMA_apu, Summ[ll], vec.DRAMA_rhs, w_vec.lambda_DRAMA, iter, osa_iter, subsets);

					if (verbose) {
						mexPrintf("Sub-iteration %d complete\n", osa_iter + 1u);
						mexEvalString("pause(.0001);");
					}

					clFinish(af_queue);

				}

				// Compute BSREM and ROSEMMAP updates if applicable
				// Otherwise simply save the current iterate
				if (MethodList.OSEM)
					vec.OSEM(span, iter + 1u) = vec.OSEM_apu;
				if (MethodList.OSLOSEM) {
					if (MethodList.MRP)
						vec.MRP_OSEM(span, iter + 1u) = vec.MRP_OSEM_apu;
					if (MethodList.Quad)
						vec.Quad_OSEM(span, iter + 1u) = vec.Quad_OSEM_apu;
					if (MethodList.L)
						vec.L_OSEM(span, iter + 1u) = vec.L_OSEM_apu;
					if (MethodList.FMH)
						vec.FMH_OSEM(span, iter + 1u) = vec.FMH_OSEM_apu;
					if (MethodList.WeightedMean)
						vec.Weighted_OSEM(span, iter + 1u) = vec.Weighted_OSEM_apu;
					if (MethodList.TV)
						vec.TV_OSEM(span, iter + 1u) = vec.TV_OSEM_apu;
					if (MethodList.AD)
						vec.AD_OSEM(span, iter + 1u) = vec.AD_OSEM_apu;
					if (MethodList.APLS)
						vec.APLS_OSEM(span, iter + 1u) = vec.APLS_OSEM_apu;
					if (MethodList.TGV)
						vec.TGV_OSEM(span, iter + 1u) = vec.TGV_OSEM_apu;
				}

				if (MethodList.MRAMLA)
					vec.MRAMLA(span, iter + 1u) = vec.MRAMLA_apu;
				if (MethodList.MBSREM) {
					if (MethodList.MRP)
						vec.MRP_MBSREM(span, iter + 1u) = vec.MRP_MBSREM_apu;
					if (MethodList.Quad)
						vec.Quad_MBSREM(span, iter + 1u) = vec.Quad_MBSREM_apu;
					if (MethodList.L)
						vec.L_MBSREM(span, iter + 1u) = vec.L_MBSREM_apu;
					if (MethodList.FMH)
						vec.FMH_MBSREM(span, iter + 1u) = vec.FMH_MBSREM_apu;
					if (MethodList.WeightedMean)
						vec.Weighted_MBSREM(span, iter + 1u) = vec.Weighted_MBSREM_apu;
					if (MethodList.TV)
						vec.TV_MBSREM(span, iter + 1u) = vec.TV_MBSREM_apu;
					if (MethodList.AD)
						vec.AD_MBSREM(span, iter + 1u) = vec.AD_MBSREM_apu;
					if (MethodList.APLS)
						vec.APLS_MBSREM(span, iter + 1u) = vec.APLS_MBSREM_apu;
					if (MethodList.TGV)
						vec.TGV_MBSREM(span, iter + 1u) = vec.TGV_MBSREM_apu;
				}

				if (MethodList.RAMLA)
					vec.RAMLA(span, iter + 1u) = vec.RAMLA_apu;
				if (MethodList.BSREM) {
					if (MethodList.MRP) {
						array dU = MRP(vec.MRP_BSREM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.med_no_norm, im_dim);
						vec.MRP_BSREM(span, iter + 1u) = BSREM_MAP(vec.MRP_BSREM_apu, w_vec.lambda_BSREM, iter, beta.MRP_BSREM, dU, epps);
					}
					if (MethodList.Quad) {
						array dU = Quadratic_prior(vec.Quad_BSREM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, w_vec.tr_offsets, w_vec.weights_quad, im_dim);
						vec.Quad_BSREM(span, iter + 1u) = BSREM_MAP(vec.Quad_BSREM_apu, w_vec.lambda_BSREM, iter, beta.Quad_BSREM, dU, epps);
					}
					if (MethodList.L) {
						array dU = L_filter(vec.L_BSREM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L, w_vec.med_no_norm, im_dim);
						vec.L_BSREM(span, iter + 1u) = BSREM_MAP(vec.L_BSREM_apu, w_vec.lambda_BSREM, iter, beta.L_BSREM, dU, epps);
					}
					if (MethodList.FMH) {
						array dU = FMH(vec.FMH_BSREM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets, w_vec.fmh_weights, w_vec.med_no_norm,
							w_vec.alku_fmh, im_dim);
						vec.FMH_BSREM(span, iter + 1u) = BSREM_MAP(vec.FMH_BSREM_apu, w_vec.lambda_BSREM, iter, beta.FMH_BSREM, dU, epps);
					}
					if (MethodList.WeightedMean) {
						array dU = Weighted_mean(vec.Weighted_BSREM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.weighted_weights,
							w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
						vec.Weighted_BSREM(span, iter + 1u) = BSREM_MAP(vec.Weighted_BSREM_apu, w_vec.lambda_BSREM, iter, beta.Weighted_BSREM, dU, epps);
					}
					if (MethodList.TV) {
						array dU = TVprior(Nx, Ny, Nz, data, vec.TV_BSREM_apu, epps, data.TVtype, w_vec);
						vec.TV_BSREM(span, iter + 1u) = BSREM_MAP(vec.TV_BSREM_apu, w_vec.lambda_BSREM, iter, beta.TV_BSREM, dU, epps);
					}
					if (MethodList.AD) {
						array dU = AD(vec.AD_BSREM_apu, Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, w_vec.DiffusionType, w_vec.med_no_norm);
						vec.AD_BSREM(span, iter + 1u) = BSREM_MAP(vec.AD_BSREM_apu, w_vec.lambda_BSREM, iter, beta.AD_BSREM, dU, epps);
					}
					if (MethodList.APLS) {
						array dU = TVprior(Nx, Ny, Nz, data, vec.TV_BSREM_apu, epps, 4, w_vec);
						vec.APLS_BSREM(span, iter + 1u) = BSREM_MAP(vec.APLS_BSREM_apu, w_vec.lambda_BSREM, iter, beta.APLS_BSREM, dU, epps);
					}
					if (MethodList.TGV) {
						array dU = TGV(vec.TGV_BSREM_apu, Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
						vec.TGV_BSREM(span, iter + 1u) = BSREM_MAP(vec.TGV_BSREM_apu, w_vec.lambda_BSREM, iter, beta.TGV_BSREM, dU, epps);
					}
				}

				if (MethodList.ROSEM)
					vec.ROSEM(span, iter + 1u) = vec.ROSEM_apu;
				if (MethodList.ROSEMMAP) {
					if (MethodList.MRP) {
						array dU = MRP(vec.MRP_ROSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.med_no_norm, im_dim);
						vec.MRP_ROSEM(span, iter + 1u) = BSREM_MAP(vec.MRP_ROSEM_apu, w_vec.lambda_ROSEM, iter, beta.MRP_ROSEM, dU, epps);
					}
					if (MethodList.Quad) {
						array dU = Quadratic_prior(vec.Quad_ROSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, w_vec.tr_offsets, w_vec.weights_quad, im_dim);
						vec.Quad_ROSEM(span, iter + 1u) = BSREM_MAP(vec.Quad_ROSEM_apu, w_vec.lambda_ROSEM, iter, beta.Quad_ROSEM, dU, epps);
					}
					if (MethodList.L) {
						array dU = L_filter(vec.L_ROSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L, w_vec.med_no_norm, im_dim);
						vec.L_ROSEM(span, iter + 1u) = BSREM_MAP(vec.L_ROSEM_apu, w_vec.lambda_ROSEM, iter, beta.L_ROSEM, dU, epps);
					}
					if (MethodList.FMH) {
						array dU = FMH(vec.FMH_ROSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets, w_vec.fmh_weights, w_vec.med_no_norm,
							w_vec.alku_fmh, im_dim);
						vec.FMH_ROSEM(span, iter + 1u) = BSREM_MAP(vec.FMH_ROSEM_apu, w_vec.lambda_ROSEM, iter, beta.FMH_ROSEM, dU, epps);
					}
					if (MethodList.WeightedMean) {
						array dU = Weighted_mean(vec.Weighted_ROSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.weighted_weights,
							w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
						vec.Weighted_ROSEM(span, iter + 1u) = BSREM_MAP(vec.Weighted_ROSEM_apu, w_vec.lambda_ROSEM, iter, beta.Weighted_ROSEM, dU, epps);
					}
					if (MethodList.TV) {
						array dU = TVprior(Nx, Ny, Nz, data, vec.TV_ROSEM_apu, epps, data.TVtype, w_vec);
						vec.TV_ROSEM(span, iter + 1u) = BSREM_MAP(vec.TV_ROSEM_apu, w_vec.lambda_ROSEM, iter, beta.TV_ROSEM, dU, epps);
					}
					if (MethodList.AD) {
						array dU = AD(vec.AD_ROSEM_apu, Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, w_vec.DiffusionType, w_vec.med_no_norm);
						vec.AD_ROSEM(span, iter + 1u) = BSREM_MAP(vec.AD_ROSEM_apu, w_vec.lambda_ROSEM, iter, beta.AD_ROSEM, dU, epps);
					}
					if (MethodList.APLS) {
						array dU = TVprior(Nx, Ny, Nz, data, vec.TV_ROSEM_apu, epps, 4, w_vec);
						vec.APLS_ROSEM(span, iter + 1u) = BSREM_MAP(vec.APLS_ROSEM_apu, w_vec.lambda_ROSEM, iter, beta.APLS_ROSEM, dU, epps);
					}
					if (MethodList.TGV) {
						array dU = TGV(vec.TGV_ROSEM_apu, Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
						vec.TGV_ROSEM(span, iter + 1u) = BSREM_MAP(vec.TGV_ROSEM_apu, w_vec.lambda_ROSEM, iter, beta.TGV_ROSEM, dU, epps);
					}
				}

				if (MethodList.RBI)
					vec.RBI(span, iter + 1u) = vec.RBI_apu;
				if (MethodList.RBIMAP) {
					if (MethodList.MRP)
						vec.MRP_RBI(span, iter + 1u) = vec.MRP_RBI_apu;
					if (MethodList.Quad)
						vec.Quad_RBI(span, iter + 1u) = vec.Quad_RBI_apu;
					if (MethodList.L)
						vec.L_RBI(span, iter + 1u) = vec.L_RBI_apu;
					if (MethodList.FMH)
						vec.FMH_RBI(span, iter + 1u) = vec.FMH_RBI_apu;
					if (MethodList.WeightedMean)
						vec.Weighted_RBI(span, iter + 1u) = vec.Weighted_RBI_apu;
					if (MethodList.TV)
						vec.TV_RBI(span, iter + 1u) = vec.TV_RBI_apu;
					if (MethodList.AD)
						vec.AD_RBI(span, iter + 1u) = vec.AD_RBI_apu;
					if (MethodList.APLS)
						vec.APLS_RBI(span, iter + 1u) = vec.APLS_RBI_apu;
					if (MethodList.TGV)
						vec.TGV_RBI(span, iter + 1u) = vec.TGV_RBI_apu;
				}

				if (MethodList.COSEM)
					vec.COSEM(span, iter + 1u) = vec.COSEM_apu;
				if (MethodList.OSLCOSEM > 0) {
					if (MethodList.MRP)
						vec.MRP_COSEM(span, iter + 1u) = vec.MRP_COSEM_apu;
					if (MethodList.Quad)
						vec.Quad_COSEM(span, iter + 1u) = vec.Quad_COSEM_apu;
					if (MethodList.L)
						vec.L_COSEM(span, iter + 1u) = vec.L_COSEM_apu;
					if (MethodList.FMH)
						vec.FMH_COSEM(span, iter + 1u) = vec.FMH_COSEM_apu;
					if (MethodList.WeightedMean)
						vec.Weighted_COSEM(span, iter + 1u) = vec.Weighted_COSEM_apu;
					if (MethodList.TV)
						vec.TV_COSEM(span, iter + 1u) = vec.TV_COSEM_apu;
					if (MethodList.AD)
						vec.AD_COSEM(span, iter + 1u) = vec.AD_COSEM_apu;
					if (MethodList.APLS)
						vec.APLS_COSEM(span, iter + 1u) = vec.APLS_COSEM_apu;
					if (MethodList.TGV)
						vec.TGV_COSEM(span, iter + 1u) = vec.TGV_COSEM_apu;
				}

				if (MethodList.ECOSEM)
					vec.ECOSEM(span, iter + 1u) = vec.ECOSEM_apu;

				if (MethodList.ACOSEM)
					vec.ACOSEM(span, iter + 1u) = vec.ACOSEM_apu;

				if (MethodList.DRAMA)
					vec.DRAMA(span, iter + 1u) = vec.DRAMA_apu;

				w_vec.epsilon_mramla.unlock();

				if (osem_bool && compute_norm_matrix == 0u)
					no_norm = 1u;

				if (verbose) {
					mexPrintf("Iteration %d complete\n", iter + 1u);
					mexEvalString("pause(.0001);");
				}
			}
		}

		// Transfer the device data to host MATLAB cell array
		device_to_host_cell(ArrayList, MethodList, vec, oo, cell);

		if (verbose) {
			mexPrintf("Time step %d complete\n", tt + 1u);
			mexEvalString("pause(.0001);");
		}
	}

	// Transfer memory control of all variables that weren't used
	unlock_AF_im_vectors(vec, MethodList, true, mlem_bool, osem_bool, 0u);

	// Clear OpenCL buffers
	clReleaseMemObject(d_z);
	clReleaseMemObject(d_x);
	clReleaseMemObject(d_y);
	clReleaseMemObject(d_atten);
	clReleaseMemObject(d_norm);
	clReleaseMemObject(d_pseudos);
	clReleaseMemObject(d_xcenter);
	clReleaseMemObject(d_ycenter);
	clReleaseMemObject(d_zcenter);
	for (uint32_t kk = 0u; kk < subsets; kk++) {
		clReleaseMemObject(d_lor[kk]);
		clReleaseMemObject(d_xyindex[kk]);
		clReleaseMemObject(d_zindex[kk]);
		clReleaseMemObject(d_L[kk]);
		clReleaseMemObject(d_Sino[kk]);
	}

	// Release program and kernels
	status = clReleaseProgram(program);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		mexPrintf("Failed to release program\n");
	}
	if (osem_bool) {
		status = clReleaseKernel(kernel);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			mexPrintf("Failed to release kernel\n");
		}
	}
	if (MethodList.MLEM || MethodList.OSLMLEM) {
		status = clReleaseKernel(kernel_ml);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			mexPrintf("Failed to release kernel\n");
		}
	}
	if ((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBI || MethodList.RBIMAP) && w_vec.MBSREM_prepass ||
		MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0) {
		status = clReleaseKernel(kernel_mramla);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			mexPrintf("Failed to release kernel\n");
		}
	}
	return;
}

void find_LORs(uint16_t* lor, const float* z_det, const float* x, const float* y, const uint32_t Nx, const uint32_t Ny,
	const uint32_t Nz, const float dx, const float dy, const float dz, const float bx, const float by, const float bz,
	const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices,
	const uint32_t size_x, const uint16_t TotSinos, const bool verbose, const uint32_t loop_var_par,
	const char* k_path, const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L,
	const uint8_t raw, const size_t size_z, const char* fileName, const uint32_t device, const size_t numel_x, 
	const bool force_build) {


	const uint32_t im_dim = Nx * Ny * Nz;
	bool atomic_64bit = false;

	cl_int status = CL_SUCCESS;
	cl_kernel kernel;
	cl_context af_context = afcl::getContext();
	cl_device_id af_device_id = afcl::getDeviceId();
	cl_command_queue af_queue = afcl::getQueue();
	cl_program program = NULL;

	if (force_build) {
		status = SaveProgramBinary(verbose, k_path, af_context, af_device_id, fileName, program, atomic_64bit, device);
		if (status != CL_SUCCESS) {
			std::cerr << "Error while saving binary" << std::endl;
			return;
		}
	}
	else {
		FILE *fp = NULL;// = fopen(fileName, "rb");
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
		fopen_s(&fp, fileName, "rb");
#else
		fp = fopen(fileName, "rb");
#endif
		// If the binaries do not yet exist, create them
		if (fp == NULL) {
			status = SaveProgramBinary(verbose, k_path, af_context, af_device_id, fileName, program, atomic_64bit, device);
			if (status != CL_SUCCESS)
				return;
		}
		else {
			// If the binaries exist, load them
			status = CreateProgramFromBinary(af_context, af_device_id, fp, program);
			fclose(fp);
			if (status != CL_SUCCESS) {
				mexPrintf("Failed to load OpenCL binaries\n");
				clReleaseProgram(program);
				return;
			}
			//else {
			//	mexPrintf("OpenCL binaries successfully loaded\n");
			//	mexEvalString("pause(.0001);");
			//}
		}
	}


	kernel = clCreateKernel(program, "siddon_precomp", &status);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
		mexPrintf("Failed to create OpenCL kernel\n");
		return;
	}

	precomp_siddon(af_context, af_queue, lor, z_det, x, y, Nx, Ny, Nz, dx, dy, dz, bx, by, bz, bzb, maxxx, maxyy, zmax, 
		NSlices, size_x, TotSinos, verbose, loop_var_par, pseudos, det_per_ring, prows, L, raw, size_z, im_dim,
		kernel, numel_x);



	status = clReleaseProgram(program);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
	}

	status = clReleaseKernel(kernel);
	if (status != CL_SUCCESS) {
		std::cerr << getErrorString(status) << std::endl;
	}
	return;
}
