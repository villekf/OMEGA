/**************************************************************************
* Matrix free computations for OMEGA.
* In this file the OpenCL buffers are created, calls to other necessary
* functions are made and the OpenCL kernels are launched. This file 
* contains the code for the matrix-free custom prior reconstructions in 
* OMEGA using the implementation 2.
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

using namespace af;

// Main reconstruction function
void reconstruction_custom_matrixfree(const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, const float* Sin, const uint32_t Nx, const uint32_t Ny,
	const uint32_t Nz, const uint32_t iter, const mxArray* options, const float dx, const float dy, const float dz, const float bx, const float by, const float bz,
	const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices, uint32_t* pituus, const size_t koko_l,
	const uint32_t* xy_index, const uint16_t* z_index, const uint32_t size_x, const uint32_t TotSinos, mxArray* cell, const mwSize dimmi, const bool verbose,
	const uint32_t attenuation_correction, const uint32_t normalization, const float* atten, const size_t size_atten, const float* norm, const size_t size_norm, const uint32_t subsets, const float epps, const uint8_t* rekot,
	const char* k_path, const size_t size_rekot, const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L,
	const uint8_t raw, const size_t size_z, const bool osem_bool, const char* fileName, const bool force_build, const uint32_t osa_iter, const uint32_t tt,
	const uint32_t n_subsets, const bool mlem_bool, const float tube_width, const float crystal_size_z, const float* x_center, const float* y_center, const float* z_center, const size_t size_center_x,
	const size_t size_center_y, const size_t size_center_z, const size_t size_of_x, const uint32_t projector_type, const uint32_t device) {


	const uint32_t Niter = 0u;
	const uint32_t Nxy = Nx * Ny;
	// Number of voxels
	const uint32_t im_dim = Nx * Ny * Nz;
	uint32_t oo = 0u;
	bool atomic_64bit = true;

	matlabArrays ArrayList;
	RecMethods MethodList;
	RecMethodsOpenCL MethodListOpenCL;

	get_rec_methods(options, MethodList);

	OpenCLRecMethods(MethodList, MethodListOpenCL);

	//create_matlab_output(ArrayList, dimmi, MethodList);
	create_matlab_output(ArrayList, &dimmi, MethodList, 1);

	array x00(Nx*Ny*Nz, (float*)mxGetData(mxGetField(options, 0, "x0")), afHost);

	std::vector<size_t> length(subsets);

	for (uint32_t kk = 0u; kk < subsets; kk++)
		length[kk] = pituus[kk + 1u] - pituus[kk];

	mxArray *custom_osem, *custom_mlem, *custom_bsrem, *custom_mbsrem, *custom_rosem, *custom_rbi, *custom_cosem;
	float *ele_custom_osem, *ele_custom_mlem, *ele_custom_bsrem, *ele_custom_mbsrem, *ele_custom_rosem, *ele_custom_rbi, *ele_custom_cosem;

	if (MethodList.OSLOSEM)
		custom_osem = mxCreateNumericArray(1, &dimmi, mxSINGLE_CLASS, mxREAL);
	else
		custom_osem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.OSLMLEM)
		custom_mlem = mxCreateNumericArray(1, &dimmi, mxSINGLE_CLASS, mxREAL);
	else
		custom_mlem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.BSREM)
		custom_bsrem = mxCreateNumericArray(1, &dimmi, mxSINGLE_CLASS, mxREAL);
	else
		custom_bsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.MBSREM)
		custom_mbsrem = mxCreateNumericArray(1, &dimmi, mxSINGLE_CLASS, mxREAL);
	else
		custom_mbsrem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.ROSEMMAP)
		custom_rosem = mxCreateNumericArray(1, &dimmi, mxSINGLE_CLASS, mxREAL);
	else
		custom_rosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.RBIMAP)
		custom_rbi = mxCreateNumericArray(1, &dimmi, mxSINGLE_CLASS, mxREAL);
	else
		custom_rbi = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	if (MethodList.OSLCOSEM > 0u)
		custom_cosem = mxCreateNumericArray(1, &dimmi, mxSINGLE_CLASS, mxREAL);
	else
		custom_cosem = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);

	if (MethodList.OSLOSEM)
		ele_custom_osem = (float*)mxGetData(custom_osem);
	if (MethodList.OSLMLEM)
		ele_custom_mlem = (float*)mxGetData(custom_mlem);
	if (MethodList.BSREM)
		ele_custom_bsrem = (float*)mxGetData(custom_bsrem);
	if (MethodList.MBSREM)
		ele_custom_mbsrem = (float*)mxGetData(custom_mbsrem);
	if (MethodList.ROSEMMAP)
		ele_custom_rosem = (float*)mxGetData(custom_rosem);
	if (MethodList.RBIMAP)
		ele_custom_rbi = (float*)mxGetData(custom_rbi);
	if (MethodList.OSLCOSEM > 0u)
		ele_custom_cosem = (float*)mxGetData(custom_cosem);


	const array lor(koko_l, lor1, afHost);

	array zindex, xindex, LL, pj3;

	if (raw) {
		zindex = constant(0, 1, 1, u32);
		xindex = constant(0, 1, 1, u32);
		LL = array(koko, L, afHost);
	}
	else {
		zindex = array(koko, z_index, afHost);
		xindex = array(koko, xy_index, afHost);
		LL = constant(0, 1, 1, u16);
	}

	AF_im_vectors vec;
	Beta beta;
	Weighting w_vec;
	TVdata data;
	OpenCL_im_vectors vec_opencl;

	//form_data_variables(vec, beta, w_vec, options, Nx, Ny, Nz, Niter, x00, im_dim, koko_l, MethodList, data, subsets);
	form_data_variables_custom(vec, beta, w_vec, options, Nx, Ny, Nz, Niter, im_dim, koko_l, MethodList, data, subsets, iter);

	float beta_custom_MLEM, beta_custom_OSEM, beta_custom_MBSREM, beta_custom_BSREM, beta_custom_ROSEM, beta_custom_RBI, beta_custom_COSEM;
	array custom_MLEM = constant(0.f, 1, 1), custom_OSEM = constant(0.f, 1, 1), custom_ROSEM = constant(0.f, 1, 1), custom_COSEM = constant(0.f, 1, 1),
		custom_BSREM = constant(0.f, 1, 1), custom_MBSREM = constant(0.f, 1, 1), custom_RBI = constant(0.f, 1, 1);
	array custom_MLEM_rhs = constant(0.f, 1, 1), custom_OSEM_rhs = constant(0.f, 1, 1), custom_ROSEM_rhs = constant(0.f, 1, 1), custom_COSEM_rhs = constant(0.f, 1, 1),
		custom_BSREM_rhs = constant(0.f, 1, 1), custom_MBSREM_rhs = constant(0.f, 1, 1), custom_RBI_rhs = constant(0.f, 1, 1);
	array grad_OSEM, grad_MLEM, grad_BSREM, grad_MBSREM, grad_ROSEM, grad_RBI, grad_COSEM;

	mxArray* im_vectors = mxGetField(options, 0, "im_vectors");
	if (MethodList.OSLMLEM) {
		beta_custom_MLEM = (float)mxGetScalar(mxGetField(options, 0, "beta_custom_mlem"));
		custom_MLEM = array(im_dim, 1, (float*)mxGetData(mxGetField(im_vectors, 0, "custom_MLEM_apu")), afHost);
		custom_MLEM_rhs = constant(0.f, im_dim, 1);
		grad_MLEM = array(im_dim, 1, (float*)mxGetData(mxGetField(options, 0, "grad_MLEM")), afHost);
	}
	if (MethodList.OSLOSEM) {
		beta_custom_OSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_custom_osem"));
		custom_OSEM = array(im_dim, 1, (float*)mxGetData(mxGetField(im_vectors, 0, "custom_OSEM_apu")), afHost);
		grad_OSEM = array(im_dim, 1, (float*)mxGetData(mxGetField(options, 0, "grad_OSEM")), afHost);
		custom_OSEM_rhs = constant(0.f, im_dim, 1);
	}
	if (MethodList.MBSREM) {
		beta_custom_MBSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_custom_mbsrem"));
		custom_MBSREM = array(im_dim, 1, (float*)mxGetData(mxGetField(im_vectors, 0, "custom_MBSREM_apu")), afHost);
		custom_MBSREM_rhs = constant(0.f, im_dim, 1);
		grad_MBSREM = array(im_dim, 1, (float*)mxGetData(mxGetField(options, 0, "grad_MBSREM")), afHost);
	}
	if (MethodList.BSREM) {
		beta_custom_BSREM = (float)mxGetScalar(mxGetField(options, 0, "beta_custom_bsrem"));
		custom_BSREM = array(im_dim, 1, (float*)mxGetData(mxGetField(im_vectors, 0, "custom_BSREM_apu")), afHost);
		custom_BSREM_rhs = constant(0.f, im_dim, 1);
		grad_BSREM = array(im_dim, 1, (float*)mxGetData(mxGetField(options, 0, "grad_BSREM")), afHost);
	}
	if (MethodList.ROSEMMAP) {
		beta_custom_ROSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_custom_rosem"));
		custom_ROSEM = array(im_dim, 1, (float*)mxGetData(mxGetField(im_vectors, 0, "custom_ROSEM_apu")), afHost);
		custom_ROSEM_rhs = constant(0.f, im_dim, 1);
		grad_ROSEM = array(im_dim, 1, (float*)mxGetData(mxGetField(options, 0, "grad_ROSEM")), afHost);
	}
	if (MethodList.RBIMAP) {
		beta_custom_RBI = (float)mxGetScalar(mxGetField(options, 0, "beta_custom_rbi"));
		custom_RBI = array(im_dim, 1, (float*)mxGetData(mxGetField(im_vectors, 0, "custom_RBI_apu")), afHost);
		custom_RBI_rhs = constant(0.f, im_dim, 1);
		grad_RBI = array(im_dim, 1, (float*)mxGetData(mxGetField(options, 0, "grad_RBI")), afHost);
	}
	if (MethodList.OSLCOSEM > 0) {
		beta_custom_COSEM = (float)mxGetScalar(mxGetField(options, 0, "beta_custom_cosem"));
		custom_COSEM = array(im_dim, 1, (float*)mxGetData(mxGetField(im_vectors, 0, "custom_COSEM_apu")), afHost);
		custom_COSEM_rhs = constant(0.f, im_dim, 1);
		grad_COSEM = array(im_dim, 1, (float*)mxGetData(mxGetField(options, 0, "grad_COSEM")), afHost);
	}

	float hh = 1.f / w_vec.h_ACOSEM;

	// Create the OpenCL context and command queue and assign the device
	cl_context af_context = afcl::getContext();
	cl_device_id af_device_id = afcl::getDeviceId();
	cl_command_queue af_queue = afcl::getQueue();

	cl_program program = NULL;

	cl_int status = CL_SUCCESS;

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
		if (fp == NULL) {
			status = SaveProgramBinary(verbose, k_path, af_context, af_device_id, fileName, program, atomic_64bit, device);
			if (status != CL_SUCCESS) {
				return;
			}
		}
		else {
			status = CreateProgramFromBinary(af_context, af_device_id, fp, program);
			fclose(fp);
			if (status != CL_SUCCESS) {
				clReleaseProgram(program);
				mexPrintf("Failed to load OpenCL binaries\n");
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

	cl_mem d_x, d_y, d_z, d_pseudos, d_atten, d_xcenter, d_ycenter;

	std::vector<cl_mem> d_lor(n_subsets);
	std::vector<cl_mem> d_L(n_subsets);
	std::vector<cl_mem> d_zindex(n_subsets);
	std::vector<cl_mem> d_xyindex(n_subsets);
	std::vector<cl_mem> d_Sino(n_subsets);

	if (n_subsets == 1u) {
		length[0] = length[osa_iter];
		pituus[0] = pituus[osa_iter];
	}

	status = createAndWriteBuffers(d_x, d_y, d_z, d_lor, d_L, d_zindex, d_xyindex, d_Sino, size_x, size_z, TotSinos, size_atten, prows, length, x, y, z_det, xy_index,
		z_index, lor1, L, Sin, raw, af_context, n_subsets, pituus, atten, pseudos, af_queue, d_atten, d_pseudos, d_xcenter, d_ycenter, x_center, y_center,
		size_center_x, size_center_y, size_of_x);
	if (status != CL_SUCCESS) {
		clReleaseProgram(program);
		clReleaseKernel(kernel);
		clReleaseKernel(kernel_ml);
		clReleaseKernel(kernel_mramla);
		return;
	}

	// Create the necessary buffers
	//cl_mem d_Nx = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	//cl_mem d_Ny = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	//cl_mem d_Nz = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	//cl_mem d_dz = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	//cl_mem d_dx = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	//cl_mem d_dy = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	//cl_mem d_bz = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	//cl_mem d_bx = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	//cl_mem d_by = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	//cl_mem d_bzb = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	//cl_mem d_maxxx = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	//cl_mem d_maxyy = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	//cl_mem d_zmax = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	//cl_mem d_NSlices = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	//if (raw) {
	//	d_z = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_z, NULL, &status);
	//	d_x = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_x, NULL, &status);
	//	d_y = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_x, NULL, &status);
	//}
	//else {
	//	d_z = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * TotSinos * 2, NULL, &status);
	//	d_x = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_x * 2, NULL, &status);
	//	d_y = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_x * 2, NULL, &status);
	//}
	//cl_mem d_atten = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * size_atten, NULL, &status);
	//cl_mem d_epps = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	//cl_mem d_size_x = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	//cl_mem d_TotSinos = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	//cl_mem d_attenuation_correction = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	//cl_mem d_N = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	//cl_mem d_det_per_ring = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	//cl_mem d_pseudos = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int) * prows, NULL, &status);
	//cl_mem d_pRows = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(int), NULL, &status);
	//cl_mem d_raw = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint8_t), NULL, &status);
	//cl_mem d_h = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
	//cl_mem d_epsilon_mramla = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(float), NULL, &status);
		// How many voxels does each LOR traverse
	//for (int kk = 0; kk < n_subsets; kk++) {
	//	int ii;
	//	if (n_subsets == 1 && osem_bool)
	//		ii = osa_iter;
	//	else
	//		ii = kk;
	//	// How many voxels does each LOR traverse
	//	d_lor[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[ii], NULL, &status);
	//	if (status != CL_SUCCESS) {
	//		std::cerr << getErrorString(status) << std::endl;
	//	}
	//	// Measurement data
	//	d_Sino[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[ii], NULL, &status);
	//	// Indices corresponding to the detector index (Sinogram data) or the detector number (raw data) at each measurement
	//	if (raw) {
	//		d_xyindex[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint32_t), NULL, &status);
	//		d_zindex[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint32_t), NULL, &status);
	//		d_L[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t) * length[ii] * 2, NULL, &status);
	//	}
	//	else {
	//		d_xyindex[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint32_t) * length[ii], NULL, &status);
	//		if (status != CL_SUCCESS) {
	//			std::cerr << getErrorString(status) << std::endl;
	//		}
	//		d_zindex[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint32_t) * length[ii], NULL, &status);
	//		if (status != CL_SUCCESS) {
	//			std::cerr << getErrorString(status) << std::endl;
	//		}
	//		d_L[kk] = clCreateBuffer(af_context, CL_MEM_READ_ONLY, sizeof(uint16_t), NULL, &status);
	//		if (status != CL_SUCCESS) {
	//			std::cerr << getErrorString(status) << std::endl;
	//		}
	//	}
	//}

	//if (status != CL_SUCCESS) {
	//	std::cerr << getErrorString(status) << std::endl;
	//	mexPrintf("Buffer creation failed\n");
	//}
	//else if (verbose){
	//	mexPrintf("Buffers created\n");
	//	mexEvalString("pause(.0001);");
	//}

	// assign values to the buffers
	//status = clEnqueueWriteBuffer(af_queue, d_Nx, CL_FALSE, 0, sizeof(int), &Nx, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_Ny, CL_FALSE, 0, sizeof(int), &Ny, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_Nz, CL_FALSE, 0, sizeof(int), &Nz, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_dz, CL_FALSE, 0, sizeof(float), &dz, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_dx, CL_FALSE, 0, sizeof(float), &dx, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_dy, CL_FALSE, 0, sizeof(float), &dy, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_bx, CL_FALSE, 0, sizeof(float), &bx, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_by, CL_FALSE, 0, sizeof(float), &by, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_bz, CL_FALSE, 0, sizeof(float), &bz, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_bzb, CL_FALSE, 0, sizeof(float), &bzb, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_maxxx, CL_FALSE, 0, sizeof(float), &maxxx, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_maxyy, CL_FALSE, 0, sizeof(float), &maxyy, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_zmax, CL_FALSE, 0, sizeof(float), &zmax, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_NSlices, CL_FALSE, 0, sizeof(float), &NSlices, 0, NULL, NULL);
	//if (raw) {
	//	status = clEnqueueWriteBuffer(af_queue, d_x, CL_FALSE, 0, sizeof(float) * size_x, x, 0, NULL, NULL);
	//	status = clEnqueueWriteBuffer(af_queue, d_y, CL_FALSE, 0, sizeof(float) * size_x, y, 0, NULL, NULL);
	//	status = clEnqueueWriteBuffer(af_queue, d_z, CL_FALSE, 0, sizeof(float) * size_z, z_det, 0, NULL, NULL);
	//}
	//else {
	//	status = clEnqueueWriteBuffer(af_queue, d_x, CL_FALSE, 0, sizeof(float) * size_x * 2, x, 0, NULL, NULL);
	//	status = clEnqueueWriteBuffer(af_queue, d_y, CL_FALSE, 0, sizeof(float) * size_x * 2, y, 0, NULL, NULL);
	//	status = clEnqueueWriteBuffer(af_queue, d_z, CL_FALSE, 0, sizeof(float) * TotSinos * 2, z_det, 0, NULL, NULL);
	//}
	//status = clEnqueueWriteBuffer(af_queue, d_atten, CL_FALSE, 0, sizeof(float) * size_atten, atten, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_epps, CL_FALSE, 0, sizeof(float), &epps, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_size_x, CL_FALSE, 0, sizeof(int), &size_x, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_TotSinos, CL_FALSE, 0, sizeof(int), &TotSinos, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_attenuation_correction, CL_FALSE, 0, sizeof(int), &attenuation_correction, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_N, CL_FALSE, 0, sizeof(int), &im_dim, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_det_per_ring, CL_FALSE, 0, sizeof(int), &det_per_ring, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_pseudos, CL_FALSE, 0, sizeof(int) * prows, pseudos, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_pRows, CL_FALSE, 0, sizeof(int), &prows, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_raw, CL_FALSE, 0, sizeof(uint8_t), &raw, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_h, CL_FALSE, 0, sizeof(float), &hh, 0, NULL, NULL);
	//status = clEnqueueWriteBuffer(af_queue, d_epsilon_mramla, CL_FALSE, 0, sizeof(float), &w_vec.epsilon_mramla, 0, NULL, NULL);
	//for (int kk = 0; kk < n_subsets; kk++) {
	//	int ii;
	//	if (n_subsets == 1 && osem_bool)
	//		ii = osa_iter;
	//	else
	//		ii = kk;
	//	if (raw) {
	//		status = clEnqueueWriteBuffer(af_queue, d_xyindex[kk], CL_FALSE, 0, sizeof(uint32_t), xy_index, 0, NULL, NULL);
	//		status = clEnqueueWriteBuffer(af_queue, d_zindex[kk], CL_FALSE, 0, sizeof(uint32_t), z_index, 0, NULL, NULL);
	//		status = clEnqueueWriteBuffer(af_queue, d_L[kk], CL_FALSE, 0, sizeof(uint16_t) * length[ii] * 2, &L[pituus[ii] * 2], 0, NULL, NULL);
	//	}
	//	else {
	//		status = clEnqueueWriteBuffer(af_queue, d_xyindex[kk], CL_FALSE, 0, sizeof(uint32_t) * length[ii], &xy_index[pituus[ii]], 0, NULL, NULL);
	//		if (status != CL_SUCCESS) {
	//			std::cerr << getErrorString(status) << std::endl;
	//		}
	//		status = clEnqueueWriteBuffer(af_queue, d_zindex[kk], CL_FALSE, 0, sizeof(uint32_t) * length[ii], &z_index[pituus[ii]], 0, NULL, NULL);
	//		if (status != CL_SUCCESS) {
	//			std::cerr << getErrorString(status) << std::endl;
	//		}
	//		status = clEnqueueWriteBuffer(af_queue, d_L[kk], CL_FALSE, 0, sizeof(uint16_t), L, 0, NULL, NULL);
	//		if (status != CL_SUCCESS) {
	//			std::cerr << getErrorString(status) << std::endl;
	//		}
	//	}
	//	status = clEnqueueWriteBuffer(af_queue, d_lor[kk], CL_FALSE, 0, sizeof(uint16_t) * length[ii], &lor1[pituus[ii]], 0, NULL, NULL);
	//	if (status != CL_SUCCESS) {
	//		std::cerr << getErrorString(status) << std::endl;
	//	}
	//	status = clEnqueueWriteBuffer(af_queue, d_Sino[kk], CL_TRUE, 0, sizeof(float) * length[ii], &Sin[pituus[ii]], 0, NULL, NULL);
	//	if (status != CL_SUCCESS) {
	//		std::cerr << getErrorString(status) << std::endl;
	//	}
	//}

	//if (status != CL_SUCCESS) {
	//	std::cerr << getErrorString(status) << std::endl;
	//	mexPrintf("Buffer write failed\n");
	//}
	//else if (verbose) {
		//mexPrintf("Buffers written\n");
		//mexEvalString("pause(.0001);");
	//}

	array Summ, Sino, indeksi1, indeksi2;
	//size_t length;

	//bool mlem_bool = (MethodList.MLEM || MethodList.OSLMLEM) ? true : false;

	initialize_opencl_inputs(vec, vec_opencl, MethodList, mlem_bool, osem_bool, af_context, af_queue, im_dim);


	//Sino = array(koko_l, Sin, afHost);

	// Compute the prepass phase for MRAMLA, MBSREM, RBI, COSEM, ACOSEM or ECOSEM if applicable
	if (((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBI || MethodList.RBIMAP) && w_vec.MBSREM_prepass ||
		MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0) && osa_iter == 0u && iter == 0u && tt == 0u) {


		// Initial value
		//array x00(Nx*Ny*Nz, (float*)mxGetData(mxGetField(options, 0, "x0")), afHost);

		// Create the prepass kernel
		//kernel_mramla = clCreateKernel(program, "MRAMLA_prepass", &status);

		//if (status != CL_SUCCESS) {
		//	std::cerr << getErrorString(status) << std::endl;
		//	mexPrintf("Failed to create OpenCL kernel\n");
		//	return;
		//}
		//else if (verbose) {
		//	mexPrintf("OpenCL kernel (MRAMLA & COSEM prepass) successfully created\n");
		//	mexEvalString("pause(.0001);");
		//}

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
		clSetKernelArg(kernel_mramla, 25, sizeof(cl_mem), &d_atten);
		clSetKernelArg(kernel_mramla, 26, sizeof(float), &epps);
		clSetKernelArg(kernel_mramla, 27, sizeof(uint32_t), &im_dim);
		clSetKernelArg(kernel_mramla, 28, sizeof(MethodListOpenCL), &MethodListOpenCL);
		clSetKernelArg(kernel_mramla, 29, sizeof(uint32_t), &Nxy);
		clSetKernelArg(kernel_mramla, 30, sizeof(float), &tube_width);
		clSetKernelArg(kernel_mramla, 31, sizeof(cl_mem), &d_xcenter);
		clSetKernelArg(kernel_mramla, 32, sizeof(cl_mem), &d_ycenter);

		uint32_t alku = 0u;

		if (w_vec.MBSREM_prepass)
			Summ = constant(0.f, im_dim, 1);
		else
			Summ = constant(0.f, 1, 1);

		if ((MethodList.MRAMLA || MethodList.MBSREM) && w_vec.MBSREM_prepass)
			w_vec.Amin = constant(0.f, koko_l, 1);
		// Run the prepass phase
		//MRAMLA_prepass(subsets, im_dim, pituus, lor, zindex, xindex, program, &d_Nx, &d_Ny, &d_Nz, &d_dx, &d_dy, &d_dz, &d_bz, &d_bx, &d_by, &d_bzb, 
		//	&d_maxxx, &d_maxyy, &d_zmax, &d_NSlices, &d_x, &d_y, &d_z, &d_size_x, &d_TotSinos, &d_N, &d_epps, &d_atten, MethodListOpenCL, &d_attenuation_correction, 
		//	af_queue, af_context, w_vec, Summ, rekot, Sino, koko_l, x00, vec.C_co, vec.C_aco, vec.C_osl, &d_h, alku, kernel_mramla, &d_pseudos, &d_pRows, &d_det_per_ring, LL,
		//	raw, &d_raw);
		MRAMLA_prepass(subsets, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, x00, vec.C_co,
			vec.C_aco, vec.C_osl, alku, kernel_mramla, d_L, raw, MethodListOpenCL, length);

		if (w_vec.MBSREM_prepass)
			w_vec.D = sum(Summ, 1);


		if (verbose) {
			mexPrintf("MRAMLA & COSEM prepass completed\n");
			mexEvalString("pause(.0001);");
		}
	}

	if (MethodList.MRAMLA || MethodList.MBSREM) {
		pj3 = w_vec.D / static_cast<float>(subsets);
	}

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
		clSetKernelArg(kernel, 23, sizeof(cl_mem), &d_atten);
		clSetKernelArg(kernel, 24, sizeof(float), &epps);
		clSetKernelArg(kernel, 25, sizeof(uint32_t), &im_dim);
		clSetKernelArg(kernel, 26, sizeof(cl_mem), &d_pseudos);
		clSetKernelArg(kernel, 27, sizeof(uint32_t), &prows);
		clSetKernelArg(kernel, 28, sizeof(uint32_t), &Nxy);
		clSetKernelArg(kernel, 29, sizeof(uint32_t), &det_per_ring);
		clSetKernelArg(kernel, 30, sizeof(float), &tube_width);
		clSetKernelArg(kernel, 31, sizeof(cl_mem), &d_xcenter);
		clSetKernelArg(kernel, 32, sizeof(cl_mem), &d_ycenter);
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
		clSetKernelArg(kernel_ml, 23, sizeof(cl_mem), &d_atten);
		clSetKernelArg(kernel_ml, 24, sizeof(float), &epps);
		clSetKernelArg(kernel_ml, 25, sizeof(uint32_t), &im_dim);
		clSetKernelArg(kernel_ml, 26, sizeof(cl_mem), &d_pseudos);
		clSetKernelArg(kernel_ml, 27, sizeof(uint32_t), &prows);
		clSetKernelArg(kernel_ml, 28, sizeof(uint32_t), &Nxy);
		clSetKernelArg(kernel_ml, 29, sizeof(float), &tube_width);
		clSetKernelArg(kernel_ml, 30, sizeof(cl_mem), &d_xcenter);
		clSetKernelArg(kernel_ml, 31, sizeof(cl_mem), &d_ycenter);
	}

	// Load the measurement data from the cell array

	if (MethodList.MBSREM || MethodList.MRAMLA) {
		Sino = array(koko_l, Sin, afHost);
		if (*w_vec.U == 0.f) {
			array temppi = (max)(Sino / w_vec.Amin);
			temppi.host(w_vec.U);
		}
		w_vec.epsilon_mramla = MBSREM_epsilon(Sino, epps);
	}

	// Compute MLEM separately
	if (mlem_bool) {

		Summ = constant(0.f, im_dim, 1);

		update_opencl_inputs(vec, vec_opencl, im_dim, MethodList, true, 0, af_context, 0, af_queue);

		cl_mem * d_custom_MLEM = custom_MLEM.device<cl_mem>();
		cl_mem * d_custom_MLEM_rhs = custom_MLEM_rhs.device<cl_mem>();

		cl_mem * d_Summ = Summ.device<cl_mem>();

		for (uint32_t kk = 0u; kk < n_subsets; kk++) {
			clSetKernelArg(kernel_ml, 32, sizeof(cl_mem), &d_Sino[kk]);
			clSetKernelArg(kernel_ml, 33, sizeof(cl_mem), &d_L[kk]);
			clSetKernelArg(kernel_ml, 34, sizeof(cl_mem), &d_xyindex[kk]);
			clSetKernelArg(kernel_ml, 35, sizeof(cl_mem), &d_zindex[kk]);
			clSetKernelArg(kernel_ml, 36, sizeof(cl_mem), d_Summ);
			clSetKernelArg(kernel_ml, 37, sizeof(cl_mem), &d_lor[kk]);
			clSetKernelArg(kernel_ml, 38, sizeof(cl_mem), vec_opencl.d_MLEM);
			clSetKernelArg(kernel_ml, 39, sizeof(cl_mem), vec_opencl.d_MRP_MLEM);
			clSetKernelArg(kernel_ml, 40, sizeof(cl_mem), vec_opencl.d_Quad_MLEM);
			clSetKernelArg(kernel_ml, 41, sizeof(cl_mem), vec_opencl.d_L_MLEM);
			clSetKernelArg(kernel_ml, 42, sizeof(cl_mem), vec_opencl.d_FMH_MLEM);
			clSetKernelArg(kernel_ml, 43, sizeof(cl_mem), vec_opencl.d_Weighted_MLEM);
			clSetKernelArg(kernel_ml, 44, sizeof(cl_mem), vec_opencl.d_TV_MLEM);
			clSetKernelArg(kernel_ml, 45, sizeof(cl_mem), vec_opencl.d_AD_MLEM);
			clSetKernelArg(kernel_ml, 46, sizeof(cl_mem), vec_opencl.d_APLS_MLEM);
			clSetKernelArg(kernel_ml, 47, sizeof(cl_mem), vec_opencl.d_TGV_MLEM);
			clSetKernelArg(kernel_ml, 48, sizeof(cl_mem), vec_opencl.d_MLEM_rhs);
			clSetKernelArg(kernel_ml, 49, sizeof(cl_mem), vec_opencl.d_MRP_MLEM_rhs);
			clSetKernelArg(kernel_ml, 50, sizeof(cl_mem), vec_opencl.d_Quad_MLEM_rhs);
			clSetKernelArg(kernel_ml, 51, sizeof(cl_mem), vec_opencl.d_L_MLEM_rhs);
			clSetKernelArg(kernel_ml, 52, sizeof(cl_mem), vec_opencl.d_FMH_MLEM_rhs);
			clSetKernelArg(kernel_ml, 53, sizeof(cl_mem), vec_opencl.d_Weighted_MLEM_rhs);
			clSetKernelArg(kernel_ml, 54, sizeof(cl_mem), vec_opencl.d_TV_MLEM_rhs);
			clSetKernelArg(kernel_ml, 55, sizeof(cl_mem), vec_opencl.d_AD_MLEM_rhs);
			clSetKernelArg(kernel_ml, 56, sizeof(cl_mem), vec_opencl.d_APLS_MLEM_rhs);
			clSetKernelArg(kernel_ml, 57, sizeof(cl_mem), vec_opencl.d_TGV_MLEM_rhs);
			clSetKernelArg(kernel_ml, 58, sizeof(cl_mem), d_custom_MLEM);
			clSetKernelArg(kernel_ml, 59, sizeof(cl_mem), d_custom_MLEM_rhs);
			status = clEnqueueNDRangeKernel(af_queue, kernel_ml, 1u, NULL, &length[kk], NULL, 0, NULL, NULL);
			clFinish(af_queue);
		}

		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			mexPrintf("Failed to execute the OpenCL kernel\n");
		}
		//else if (verbose) {
		//	mexPrintf("OpenCL kernel executed successfully\n");
		//	mexEvalString("pause(.0001);");
		//}

		clFinish(af_queue);
		// Transfer memory control back to ArrayFire (ML-methods)
		unlock_AF_im_vectors(vec, MethodList, false, true, false, 0);
		Summ.unlock();

		// Prevents division by zero
		Summ(Summ == 0.f) = epps;
		if (MethodList.MLEM)
			vec.MLEM_apu = MLEM(vec.MLEM_apu, Summ, vec.MLEM_rhs);
		if (MethodList.OSLMLEM) {
			if (MethodList.MRP) {
				array dU = MRP(vec.MRP_MLEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.med_no_norm, im_dim);
				vec.MRP_MLEM_apu = OSL_MLEM(vec.MRP_MLEM_apu, Summ, vec.MRP_MLEM_rhs, dU, beta.MRP_MLEM);
			}
			if (MethodList.Quad) {
				array dU = Quadratic_prior(vec.Quad_MLEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, w_vec.tr_offsets, w_vec.weights_quad, im_dim);
				vec.Quad_MLEM_apu = OSL_MLEM(vec.Quad_MLEM_apu, Summ, vec.Quad_MLEM_rhs, dU, beta.Quad_MLEM);
			}
			if (MethodList.L) {
				array dU = L_filter(vec.L_MLEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L, w_vec.med_no_norm, im_dim);
				vec.L_MLEM_apu = OSL_MLEM(vec.L_MLEM_apu, Summ, vec.L_MLEM_rhs, dU, beta.L_MLEM);
			}
			if (MethodList.FMH) {
				array dU = FMH(vec.FMH_MLEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets, w_vec.fmh_weights, w_vec.med_no_norm,
					w_vec.alku_fmh, im_dim);
				vec.FMH_MLEM_apu = OSL_MLEM(vec.FMH_MLEM_apu, Summ, vec.FMH_MLEM_rhs, dU, beta.FMH_MLEM);
			}
			if (MethodList.WeightedMean) {
				array dU = Weighted_mean(vec.Weighted_MLEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.weighted_weights,
					w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
				vec.Weighted_MLEM_apu = OSL_MLEM(vec.Weighted_MLEM_apu, Summ, vec.Weighted_MLEM_rhs, dU, beta.Weighted_MLEM);
			}
			if (MethodList.TV) {
				array dU = TVprior(Nx, Ny, Nz, data, vec.TV_MLEM_apu, epps, data.TVtype, w_vec);
				vec.TV_MLEM_apu = OSL_MLEM(vec.TV_MLEM_apu, Summ, vec.TV_MLEM_rhs, dU, beta.TV_MLEM);
			}
			if (MethodList.AD) {
				array dU = AD(vec.AD_MLEM_apu, Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, w_vec.DiffusionType, w_vec.med_no_norm);
				vec.AD_MLEM_apu = OSL_MLEM(vec.AD_MLEM_apu, Summ, vec.AD_MLEM_rhs, dU, beta.AD_MLEM);
			}
			if (MethodList.APLS) {
				array dU = TVprior(Nx, Ny, Nz, data, vec.TV_MLEM_apu, epps, 4, w_vec);
				vec.APLS_MLEM_apu = OSL_MLEM(vec.APLS_MLEM_apu, Summ, vec.APLS_MLEM_rhs, dU, beta.APLS_MLEM);
			}
			if (MethodList.TGV) {
				array dU = TGV(vec.TGV_MLEM_apu, Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
				vec.TGV_MLEM_apu = OSL_MLEM(vec.TGV_MLEM_apu, Summ, vec.TGV_MLEM_rhs, dU, beta.TGV_MLEM);
			}
			custom_MLEM = OSL_MLEM(custom_MLEM, Summ, custom_MLEM_rhs, grad_MLEM, beta_custom_MLEM);
		}
	}

	//		// Compute any of the other algorithms, if applicable
	if (osem_bool) {

		//array du_apu;

		//array OSEM_rhs = constant(1.f, length, 1);

		cl_mem *d_epsilon_mramla = w_vec.epsilon_mramla.device<cl_mem>();

		//array Lo, zo, xo, sub_index_array, uu, apu_lor;

		Summ = constant(0.f, im_dim, 1);

		array uu;

		Summ = constant(0.f, im_dim, 1);

		update_opencl_inputs(vec, vec_opencl, im_dim, MethodList, false, osa_iter, af_context, length[osa_iter], af_queue);


		cl_mem *d_Summ = Summ.device<cl_mem>();


		cl_mem * d_custom_OSEM = custom_OSEM.device<cl_mem>();
		cl_mem * d_custom_OSEM_rhs = custom_OSEM_rhs.device<cl_mem>();
		cl_mem * d_custom_BSREM = custom_BSREM.device<cl_mem>();
		cl_mem * d_custom_BSREM_rhs = custom_BSREM_rhs.device<cl_mem>();
		cl_mem * d_custom_MBSREM = custom_MBSREM.device<cl_mem>();
		cl_mem * d_custom_MBSREM_rhs = custom_MBSREM_rhs.device<cl_mem>();
		cl_mem * d_custom_ROSEM = custom_ROSEM.device<cl_mem>();
		cl_mem * d_custom_ROSEM_rhs = custom_ROSEM_rhs.device<cl_mem>();
		cl_mem * d_custom_RBI = custom_RBI.device<cl_mem>();
		cl_mem * d_custom_RBI_rhs = custom_RBI_rhs.device<cl_mem>();
		cl_mem * d_custom_COSEM = custom_COSEM.device<cl_mem>();
		cl_mem * d_custom_COSEM_rhs = custom_COSEM_rhs.device<cl_mem>();


		//cl_mem *d_Summ = Summ.device<cl_mem>();
		//cl_mem *d_Sino = uu.device<cl_mem>();



		// Set kernel arugments
		//clSetKernelArg(kernel, 0, sizeof(MethodListOpenCL), &MethodListOpenCL);
		clSetKernelArg(kernel, 33, sizeof(cl_mem), d_Summ);
		clSetKernelArg(kernel, 34, sizeof(cl_mem), &d_lor[0]);
		clSetKernelArg(kernel, 35, sizeof(cl_mem), &d_xyindex[0]);
		clSetKernelArg(kernel, 36, sizeof(cl_mem), &d_zindex[0]);
		clSetKernelArg(kernel, 37, sizeof(cl_mem), &d_L[0]);
		clSetKernelArg(kernel, 38, sizeof(cl_mem), d_epsilon_mramla);
		clSetKernelArg(kernel, 39, sizeof(cl_mem), &d_Sino[0]);
		clSetKernelArg(kernel, 40, sizeof(cl_mem), vec_opencl.d_OSEM);
		clSetKernelArg(kernel, 41, sizeof(cl_mem), vec_opencl.d_RAMLA);
		clSetKernelArg(kernel, 42, sizeof(cl_mem), vec_opencl.d_MRAMLA);
		clSetKernelArg(kernel, 43, sizeof(cl_mem), vec_opencl.d_ROSEM);
		clSetKernelArg(kernel, 44, sizeof(cl_mem), vec_opencl.d_RBI);
		clSetKernelArg(kernel, 45, sizeof(cl_mem), vec_opencl.d_DRAMA);
		clSetKernelArg(kernel, 46, sizeof(cl_mem), vec_opencl.d_COSEM);
		clSetKernelArg(kernel, 47, sizeof(cl_mem), vec_opencl.d_ACOSEM);
		clSetKernelArg(kernel, 48, sizeof(cl_mem), vec_opencl.d_MRP_OSEM);
		clSetKernelArg(kernel, 49, sizeof(cl_mem), vec_opencl.d_Quad_OSEM);
		clSetKernelArg(kernel, 50, sizeof(cl_mem), vec_opencl.d_L_OSEM);
		clSetKernelArg(kernel, 51, sizeof(cl_mem), vec_opencl.d_FMH_OSEM);
		clSetKernelArg(kernel, 52, sizeof(cl_mem), vec_opencl.d_Weighted_OSEM);
		clSetKernelArg(kernel, 53, sizeof(cl_mem), vec_opencl.d_TV_OSEM);
		clSetKernelArg(kernel, 54, sizeof(cl_mem), vec_opencl.d_AD_OSEM);
		clSetKernelArg(kernel, 55, sizeof(cl_mem), vec_opencl.d_APLS_OSEM);
		clSetKernelArg(kernel, 56, sizeof(cl_mem), vec_opencl.d_TGV_OSEM);
		clSetKernelArg(kernel, 57, sizeof(cl_mem), vec_opencl.d_MRP_BSREM);
		clSetKernelArg(kernel, 58, sizeof(cl_mem), vec_opencl.d_Quad_BSREM);
		clSetKernelArg(kernel, 59, sizeof(cl_mem), vec_opencl.d_L_BSREM);
		clSetKernelArg(kernel, 60, sizeof(cl_mem), vec_opencl.d_FMH_BSREM);
		clSetKernelArg(kernel, 61, sizeof(cl_mem), vec_opencl.d_Weighted_BSREM);
		clSetKernelArg(kernel, 62, sizeof(cl_mem), vec_opencl.d_TV_BSREM);
		clSetKernelArg(kernel, 63, sizeof(cl_mem), vec_opencl.d_AD_BSREM);
		clSetKernelArg(kernel, 64, sizeof(cl_mem), vec_opencl.d_APLS_BSREM);
		clSetKernelArg(kernel, 65, sizeof(cl_mem), vec_opencl.d_TGV_BSREM);
		clSetKernelArg(kernel, 66, sizeof(cl_mem), vec_opencl.d_MRP_MBSREM);
		clSetKernelArg(kernel, 67, sizeof(cl_mem), vec_opencl.d_Quad_MBSREM);
		clSetKernelArg(kernel, 68, sizeof(cl_mem), vec_opencl.d_L_MBSREM);
		clSetKernelArg(kernel, 69, sizeof(cl_mem), vec_opencl.d_FMH_MBSREM);
		clSetKernelArg(kernel, 70, sizeof(cl_mem), vec_opencl.d_Weighted_MBSREM);
		clSetKernelArg(kernel, 71, sizeof(cl_mem), vec_opencl.d_TV_MBSREM);
		clSetKernelArg(kernel, 72, sizeof(cl_mem), vec_opencl.d_AD_MBSREM);
		clSetKernelArg(kernel, 73, sizeof(cl_mem), vec_opencl.d_APLS_MBSREM);
		clSetKernelArg(kernel, 74, sizeof(cl_mem), vec_opencl.d_TGV_MBSREM);
		clSetKernelArg(kernel, 75, sizeof(cl_mem), vec_opencl.d_MRP_ROSEM);
		clSetKernelArg(kernel, 76, sizeof(cl_mem), vec_opencl.d_Quad_ROSEM);
		clSetKernelArg(kernel, 77, sizeof(cl_mem), vec_opencl.d_L_ROSEM);
		clSetKernelArg(kernel, 78, sizeof(cl_mem), vec_opencl.d_FMH_ROSEM);
		clSetKernelArg(kernel, 79, sizeof(cl_mem), vec_opencl.d_Weighted_ROSEM);
		clSetKernelArg(kernel, 80, sizeof(cl_mem), vec_opencl.d_TV_ROSEM);
		clSetKernelArg(kernel, 81, sizeof(cl_mem), vec_opencl.d_AD_ROSEM);
		clSetKernelArg(kernel, 82, sizeof(cl_mem), vec_opencl.d_APLS_ROSEM);
		clSetKernelArg(kernel, 83, sizeof(cl_mem), vec_opencl.d_TGV_ROSEM);
		clSetKernelArg(kernel, 84, sizeof(cl_mem), vec_opencl.d_MRP_RBI);
		clSetKernelArg(kernel, 85, sizeof(cl_mem), vec_opencl.d_Quad_RBI);
		clSetKernelArg(kernel, 86, sizeof(cl_mem), vec_opencl.d_L_RBI);
		clSetKernelArg(kernel, 87, sizeof(cl_mem), vec_opencl.d_FMH_RBI);
		clSetKernelArg(kernel, 88, sizeof(cl_mem), vec_opencl.d_Weighted_RBI);
		clSetKernelArg(kernel, 89, sizeof(cl_mem), vec_opencl.d_TV_RBI);
		clSetKernelArg(kernel, 90, sizeof(cl_mem), vec_opencl.d_AD_RBI);
		clSetKernelArg(kernel, 91, sizeof(cl_mem), vec_opencl.d_APLS_RBI);
		clSetKernelArg(kernel, 92, sizeof(cl_mem), vec_opencl.d_TGV_RBI);
		clSetKernelArg(kernel, 93, sizeof(cl_mem), vec_opencl.d_MRP_COSEM);
		clSetKernelArg(kernel, 94, sizeof(cl_mem), vec_opencl.d_Quad_COSEM);
		clSetKernelArg(kernel, 95, sizeof(cl_mem), vec_opencl.d_L_COSEM);
		clSetKernelArg(kernel, 96, sizeof(cl_mem), vec_opencl.d_FMH_COSEM);
		clSetKernelArg(kernel, 97, sizeof(cl_mem), vec_opencl.d_Weighted_COSEM);
		clSetKernelArg(kernel, 98, sizeof(cl_mem), vec_opencl.d_TV_COSEM);
		clSetKernelArg(kernel, 99, sizeof(cl_mem), vec_opencl.d_AD_COSEM);
		clSetKernelArg(kernel, 100, sizeof(cl_mem), vec_opencl.d_APLS_COSEM);
		clSetKernelArg(kernel, 101, sizeof(cl_mem), vec_opencl.d_TGV_COSEM);
		clSetKernelArg(kernel, 102, sizeof(cl_mem), vec_opencl.d_OSEM_rhs);
		clSetKernelArg(kernel, 103, sizeof(cl_mem), vec_opencl.d_RAMLA_rhs);
		clSetKernelArg(kernel, 104, sizeof(cl_mem), vec_opencl.d_MRAMLA_rhs);
		clSetKernelArg(kernel, 105, sizeof(cl_mem), vec_opencl.d_ROSEM_rhs);
		clSetKernelArg(kernel, 106, sizeof(cl_mem), vec_opencl.d_RBI_rhs);
		clSetKernelArg(kernel, 107, sizeof(cl_mem), vec_opencl.d_DRAMA_rhs);
		clSetKernelArg(kernel, 108, sizeof(cl_mem), vec_opencl.d_COSEM_rhs);
		clSetKernelArg(kernel, 109, sizeof(cl_mem), vec_opencl.d_ACOSEM_rhs);
		clSetKernelArg(kernel, 110, sizeof(cl_mem), vec_opencl.d_MRP_OSEM_rhs);
		clSetKernelArg(kernel, 111, sizeof(cl_mem), vec_opencl.d_Quad_OSEM_rhs);
		clSetKernelArg(kernel, 112, sizeof(cl_mem), vec_opencl.d_L_OSEM_rhs);
		clSetKernelArg(kernel, 113, sizeof(cl_mem), vec_opencl.d_FMH_OSEM_rhs);
		clSetKernelArg(kernel, 114, sizeof(cl_mem), vec_opencl.d_Weighted_OSEM_rhs);
		clSetKernelArg(kernel, 115, sizeof(cl_mem), vec_opencl.d_TV_OSEM_rhs);
		clSetKernelArg(kernel, 116, sizeof(cl_mem), vec_opencl.d_AD_OSEM_rhs);
		clSetKernelArg(kernel, 117, sizeof(cl_mem), vec_opencl.d_APLS_OSEM_rhs);
		clSetKernelArg(kernel, 118, sizeof(cl_mem), vec_opencl.d_TGV_OSEM_rhs);
		clSetKernelArg(kernel, 119, sizeof(cl_mem), vec_opencl.d_MRP_BSREM_rhs);
		clSetKernelArg(kernel, 120, sizeof(cl_mem), vec_opencl.d_Quad_BSREM_rhs);
		clSetKernelArg(kernel, 121, sizeof(cl_mem), vec_opencl.d_L_BSREM_rhs);
		clSetKernelArg(kernel, 122, sizeof(cl_mem), vec_opencl.d_FMH_BSREM_rhs);
		clSetKernelArg(kernel, 123, sizeof(cl_mem), vec_opencl.d_Weighted_BSREM_rhs);
		clSetKernelArg(kernel, 124, sizeof(cl_mem), vec_opencl.d_TV_BSREM_rhs);
		clSetKernelArg(kernel, 125, sizeof(cl_mem), vec_opencl.d_AD_BSREM_rhs);
		clSetKernelArg(kernel, 126, sizeof(cl_mem), vec_opencl.d_APLS_BSREM_rhs);
		clSetKernelArg(kernel, 127, sizeof(cl_mem), vec_opencl.d_TGV_BSREM_rhs);
		clSetKernelArg(kernel, 128, sizeof(cl_mem), vec_opencl.d_MRP_MBSREM_rhs);
		clSetKernelArg(kernel, 129, sizeof(cl_mem), vec_opencl.d_Quad_MBSREM_rhs);
		clSetKernelArg(kernel, 130, sizeof(cl_mem), vec_opencl.d_L_MBSREM_rhs);
		clSetKernelArg(kernel, 131, sizeof(cl_mem), vec_opencl.d_FMH_MBSREM_rhs);
		clSetKernelArg(kernel, 132, sizeof(cl_mem), vec_opencl.d_Weighted_MBSREM_rhs);
		clSetKernelArg(kernel, 133, sizeof(cl_mem), vec_opencl.d_TV_MBSREM_rhs);
		clSetKernelArg(kernel, 134, sizeof(cl_mem), vec_opencl.d_AD_MBSREM_rhs);
		clSetKernelArg(kernel, 135, sizeof(cl_mem), vec_opencl.d_APLS_MBSREM_rhs);
		clSetKernelArg(kernel, 136, sizeof(cl_mem), vec_opencl.d_TGV_MBSREM_rhs);
		clSetKernelArg(kernel, 137, sizeof(cl_mem), vec_opencl.d_MRP_ROSEM_rhs);
		clSetKernelArg(kernel, 138, sizeof(cl_mem), vec_opencl.d_Quad_ROSEM_rhs);
		clSetKernelArg(kernel, 139, sizeof(cl_mem), vec_opencl.d_L_ROSEM_rhs);
		clSetKernelArg(kernel, 140, sizeof(cl_mem), vec_opencl.d_FMH_ROSEM_rhs);
		clSetKernelArg(kernel, 141, sizeof(cl_mem), vec_opencl.d_Weighted_ROSEM_rhs);
		clSetKernelArg(kernel, 142, sizeof(cl_mem), vec_opencl.d_TV_ROSEM_rhs);
		clSetKernelArg(kernel, 143, sizeof(cl_mem), vec_opencl.d_AD_ROSEM_rhs);
		clSetKernelArg(kernel, 144, sizeof(cl_mem), vec_opencl.d_APLS_ROSEM_rhs);
		clSetKernelArg(kernel, 145, sizeof(cl_mem), vec_opencl.d_TGV_ROSEM_rhs);
		clSetKernelArg(kernel, 146, sizeof(cl_mem), vec_opencl.d_MRP_RBI_rhs);
		clSetKernelArg(kernel, 147, sizeof(cl_mem), vec_opencl.d_Quad_RBI_rhs);
		clSetKernelArg(kernel, 148, sizeof(cl_mem), vec_opencl.d_L_RBI_rhs);
		clSetKernelArg(kernel, 149, sizeof(cl_mem), vec_opencl.d_FMH_RBI_rhs);
		clSetKernelArg(kernel, 150, sizeof(cl_mem), vec_opencl.d_Weighted_RBI_rhs);
		clSetKernelArg(kernel, 151, sizeof(cl_mem), vec_opencl.d_TV_RBI_rhs);
		clSetKernelArg(kernel, 152, sizeof(cl_mem), vec_opencl.d_AD_RBI_rhs);
		clSetKernelArg(kernel, 153, sizeof(cl_mem), vec_opencl.d_APLS_RBI_rhs);
		clSetKernelArg(kernel, 154, sizeof(cl_mem), vec_opencl.d_TGV_RBI_rhs);
		clSetKernelArg(kernel, 155, sizeof(cl_mem), vec_opencl.d_MRP_COSEM_rhs);
		clSetKernelArg(kernel, 156, sizeof(cl_mem), vec_opencl.d_Quad_COSEM_rhs);
		clSetKernelArg(kernel, 157, sizeof(cl_mem), vec_opencl.d_L_COSEM_rhs);
		clSetKernelArg(kernel, 158, sizeof(cl_mem), vec_opencl.d_FMH_COSEM_rhs);
		clSetKernelArg(kernel, 159, sizeof(cl_mem), vec_opencl.d_Weighted_COSEM_rhs);
		clSetKernelArg(kernel, 160, sizeof(cl_mem), vec_opencl.d_TV_COSEM_rhs);
		clSetKernelArg(kernel, 161, sizeof(cl_mem), vec_opencl.d_AD_COSEM_rhs);
		clSetKernelArg(kernel, 162, sizeof(cl_mem), vec_opencl.d_APLS_COSEM_rhs);
		clSetKernelArg(kernel, 163, sizeof(cl_mem), vec_opencl.d_TGV_COSEM_rhs);
		clSetKernelArg(kernel, 164, sizeof(cl_mem), d_custom_OSEM);
		clSetKernelArg(kernel, 165, sizeof(cl_mem), d_custom_OSEM_rhs);
		clSetKernelArg(kernel, 166, sizeof(cl_mem), d_custom_BSREM);
		clSetKernelArg(kernel, 167, sizeof(cl_mem), d_custom_BSREM_rhs);
		clSetKernelArg(kernel, 168, sizeof(cl_mem), d_custom_MBSREM);
		clSetKernelArg(kernel, 169, sizeof(cl_mem), d_custom_MBSREM_rhs);
		clSetKernelArg(kernel, 170, sizeof(cl_mem), d_custom_ROSEM);
		clSetKernelArg(kernel, 171, sizeof(cl_mem), d_custom_ROSEM_rhs);
		clSetKernelArg(kernel, 172, sizeof(cl_mem), d_custom_RBI);
		clSetKernelArg(kernel, 173, sizeof(cl_mem), d_custom_RBI_rhs);
		clSetKernelArg(kernel, 174, sizeof(cl_mem), d_custom_COSEM);
		clSetKernelArg(kernel, 175, sizeof(cl_mem), d_custom_COSEM_rhs);
		// Compute the kernel
		status = clEnqueueNDRangeKernel(af_queue, kernel, 1u, NULL, &length[0], NULL, 0, NULL, NULL);

		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			mexPrintf("Failed to execute the OpenCL kernel\n");
		}
		else if (verbose) {
			mexPrintf("OpenCL kernel executed successfully\n");
			mexEvalString("pause(.0001);");
		}

		if (MethodList.ACOSEM || MethodList.OSLCOSEM == 1) {
			uu = afcl::array(length[osa_iter], d_Sino[0], f32, true);
		}
		clFinish(af_queue);
		// Transfer memory control back to ArrayFire (OS-methods)
		Summ.unlock();
		unlock_AF_im_vectors(vec, MethodList, false, false, true, osa_iter);
		// Prevent division by zero
		Summ(Summ == 0.f) = epps;
		//Summ += epps;
//				// Compute the (matrix free) algorithms
		if (MethodList.OSEM || MethodList.ECOSEM)
			vec.OSEM_apu = OSEM(vec.OSEM_apu, Summ, vec.OSEM_rhs);
		if (MethodList.OSLOSEM) {
			if (MethodList.MRP) {
				array dU = MRP(vec.MRP_OSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.med_no_norm, im_dim);
				vec.MRP_OSEM_apu = OSL_OSEM(vec.MRP_OSEM_apu, Summ, vec.MRP_OSEM_rhs, dU, beta.MRP_OSEM);
			}
			if (MethodList.Quad) {
				array dU = Quadratic_prior(vec.Quad_OSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, w_vec.tr_offsets, w_vec.weights_quad, im_dim);
				vec.Quad_OSEM_apu = OSL_OSEM(vec.Quad_OSEM_apu, Summ, vec.Quad_OSEM_rhs, dU, beta.Quad_OSEM);
			}
			if (MethodList.L) {
				array dU = L_filter(vec.L_OSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L, w_vec.med_no_norm, im_dim);
				vec.L_OSEM_apu = OSL_OSEM(vec.L_OSEM_apu, Summ, vec.L_OSEM_rhs, dU, beta.L_OSEM);
			}
			if (MethodList.FMH) {
				array dU = FMH(vec.FMH_OSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets, w_vec.fmh_weights, w_vec.med_no_norm,
					w_vec.alku_fmh, im_dim);
				vec.FMH_OSEM_apu = OSL_OSEM(vec.FMH_OSEM_apu, Summ, vec.FMH_OSEM_rhs, dU, beta.FMH_OSEM);
			}
			if (MethodList.WeightedMean) {
				array dU = Weighted_mean(vec.Weighted_OSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.weighted_weights,
					w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
				vec.Weighted_OSEM_apu = OSL_OSEM(vec.Weighted_OSEM_apu, Summ, vec.Weighted_OSEM_rhs, dU, beta.Weighted_OSEM);
			}
			if (MethodList.TV) {
				array dU = TVprior(Nx, Ny, Nz, data, vec.TV_OSEM_apu, epps, data.TVtype, w_vec);
				vec.TV_OSEM_apu = OSL_OSEM(vec.TV_OSEM_apu, Summ, vec.TV_OSEM_rhs, dU, beta.TV_OSEM);
			}
			if (MethodList.AD) {
				if (osa_iter == 0)
					if (MethodList.OSEM)
						vec.AD_OSEM_apu = vec.OSEM_apu;
					else
						vec.AD_OSEM_apu = OSEM(vec.AD_OSEM_apu, Summ, vec.AD_OSEM_rhs);
				else {
					array dU = AD(vec.AD_OSEM_apu, Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, w_vec.DiffusionType, w_vec.med_no_norm);
					vec.AD_OSEM_apu = OSL_OSEM(vec.AD_OSEM_apu, Summ, vec.AD_OSEM_rhs, dU, beta.AD_OSEM);
				}
			}
			if (MethodList.APLS) {
				array dU = TVprior(Nx, Ny, Nz, data, vec.TV_OSEM_apu, epps, 4, w_vec);
				vec.APLS_OSEM_apu = OSL_OSEM(vec.APLS_OSEM_apu, Summ, vec.APLS_OSEM_rhs, dU, beta.APLS_OSEM);
			}
			if (MethodList.TGV) {
				array dU = TGV(vec.TGV_OSEM_apu, Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
				vec.TGV_OSEM_apu = OSL_OSEM(vec.TGV_OSEM_apu, Summ, vec.TGV_OSEM_rhs, dU, beta.TGV_OSEM);
			}
			custom_OSEM = OSL_OSEM(custom_OSEM, Summ, custom_OSEM_rhs, grad_OSEM, beta_custom_OSEM);
		}

		if (MethodList.MRAMLA)
			vec.MRAMLA_apu = MBSREM(vec.MRAMLA_apu, vec.MRAMLA_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM, iter, im_dim, 0.f, af::constant(0.f, 1, 1),
				Summ, epps);
		if (MethodList.MBSREM) {
			if (MethodList.MRP) {
				array dU = MRP(vec.MRP_MBSREM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.med_no_norm, im_dim);
				vec.MRP_MBSREM_apu = MBSREM(vec.MRP_MBSREM_apu, vec.MRP_MBSREM_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.MRP_MBSREM,
					dU, Summ, epps);
			}
			if (MethodList.Quad) {
				array dU = Quadratic_prior(vec.Quad_MBSREM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, w_vec.tr_offsets, w_vec.weights_quad, im_dim);
				vec.Quad_MBSREM_apu = MBSREM(vec.Quad_MBSREM_apu, vec.MRAMLA_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM,
					iter, im_dim, beta.Quad_MBSREM, dU, Summ, epps);
			}
			if (MethodList.L) {
				array dU = L_filter(vec.L_MBSREM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L, w_vec.med_no_norm, im_dim);
				vec.L_MBSREM_apu = MBSREM(vec.L_MBSREM_apu, vec.L_MBSREM_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM,
					iter, im_dim, beta.L_MBSREM, dU, Summ, epps);
			}
			if (MethodList.FMH) {
				array dU = FMH(vec.FMH_MBSREM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets, w_vec.fmh_weights, w_vec.med_no_norm,
					w_vec.alku_fmh, im_dim);
				vec.FMH_MBSREM_apu = MBSREM(vec.FMH_MBSREM_apu, vec.FMH_MBSREM_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM,
					iter, im_dim, beta.FMH_MBSREM, dU, Summ, epps);
			}
			if (MethodList.WeightedMean) {
				array dU = Weighted_mean(vec.Weighted_MBSREM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.weighted_weights,
					w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
				vec.Weighted_MBSREM_apu = MBSREM(vec.Weighted_MBSREM_apu, vec.Weighted_MBSREM_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM,
					iter, im_dim, beta.Weighted_MBSREM, dU, Summ, epps);
			}
			if (MethodList.TV) {
				array dU = TVprior(Nx, Ny, Nz, data, vec.TV_MBSREM_apu, epps, data.TVtype, w_vec);
				vec.TV_MBSREM_apu = MBSREM(vec.TV_MBSREM_apu, vec.TV_MBSREM_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM,
					iter, im_dim, beta.TV_MBSREM, dU, Summ, epps);
			}
			if (MethodList.AD) {
				array dU = AD(vec.AD_MBSREM_apu, Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, w_vec.DiffusionType, w_vec.med_no_norm);
				vec.AD_MBSREM_apu = MBSREM(vec.AD_MBSREM_apu, vec.AD_MBSREM_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM,
					iter, im_dim, beta.AD_MBSREM, dU, Summ, epps);
			}
			if (MethodList.APLS) {
				array dU = TVprior(Nx, Ny, Nz, data, vec.TV_MBSREM_apu, epps, 4, w_vec);
				vec.TV_MBSREM_apu = MBSREM(vec.TV_MBSREM_apu, vec.TV_MBSREM_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM,
					iter, im_dim, beta.APLS_MBSREM, dU, Summ, epps);
			}
			if (MethodList.TGV) {
				array dU = TGV(vec.TGV_MBSREM_apu, Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
				vec.TGV_MBSREM_apu = MBSREM(vec.TGV_MBSREM_apu, vec.TGV_MBSREM_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM,
					iter, im_dim, beta.TGV_MBSREM, dU, Summ, epps);
			}
			custom_MBSREM = MBSREM(custom_MBSREM, custom_MBSREM_rhs, w_vec.U, pj3, w_vec.lambda_MBSREM,
				iter, im_dim, beta_custom_MBSREM, grad_MBSREM, Summ, epps);
		}

		if (MethodList.RAMLA)
			vec.RAMLA_apu = BSREM(vec.RAMLA_apu, vec.RAMLA_rhs, w_vec.lambda_BSREM, iter);
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
			custom_BSREM = BSREM(custom_BSREM, custom_BSREM_rhs, w_vec.lambda_BSREM, iter);
		}

		if (MethodList.ROSEM)
			vec.ROSEM_apu = ROSEM(vec.ROSEM_apu, Summ, vec.ROSEM_rhs, w_vec.lambda_ROSEM, iter);
		if (MethodList.ROSEMMAP) {
			if (MethodList.MRP) {
				array dU = MRP(vec.MRP_ROSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.med_no_norm, im_dim);
				vec.MRP_ROSEM_apu = ROSEM(vec.MRP_ROSEM_apu, Summ, vec.MRP_ROSEM_rhs, w_vec.lambda_ROSEM, iter);
			}
			if (MethodList.Quad) {
				array dU = Quadratic_prior(vec.Quad_ROSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, w_vec.tr_offsets, w_vec.weights_quad, im_dim);
				vec.Quad_ROSEM_apu = ROSEM(vec.Quad_ROSEM_apu, Summ, vec.ROSEM_rhs, w_vec.lambda_ROSEM, iter);
			}
			if (MethodList.L) {
				array dU = L_filter(vec.L_ROSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L, w_vec.med_no_norm, im_dim);
				vec.L_ROSEM_apu = ROSEM(vec.L_ROSEM_apu, Summ, vec.L_ROSEM_rhs, w_vec.lambda_ROSEM, iter);
			}
			if (MethodList.FMH) {
				array dU = FMH(vec.FMH_ROSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets, w_vec.fmh_weights, w_vec.med_no_norm,
					w_vec.alku_fmh, im_dim);
				vec.FMH_ROSEM_apu = ROSEM(vec.FMH_ROSEM_apu, Summ, vec.FMH_ROSEM_rhs, w_vec.lambda_ROSEM, iter);
			}
			if (MethodList.WeightedMean) {
				array dU = Weighted_mean(vec.Weighted_ROSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.weighted_weights,
					w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
				vec.Weighted_ROSEM_apu = ROSEM(vec.Weighted_ROSEM_apu, Summ, vec.Weighted_ROSEM_rhs, w_vec.lambda_ROSEM, iter);
			}
			if (MethodList.TV) {
				array dU = TVprior(Nx, Ny, Nz, data, vec.TV_ROSEM_apu, epps, data.TVtype, w_vec);
				vec.TV_ROSEM_apu = ROSEM(vec.TV_ROSEM_apu, Summ, vec.TV_ROSEM_rhs, w_vec.lambda_ROSEM, iter);
			}
			if (MethodList.AD) {
				array dU = AD(vec.AD_ROSEM_apu, Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, w_vec.DiffusionType, w_vec.med_no_norm);
				vec.AD_ROSEM_apu = ROSEM(vec.AD_ROSEM_apu, Summ, vec.AD_ROSEM_rhs, w_vec.lambda_ROSEM, iter);
			}
			if (MethodList.APLS) {
				array dU = TVprior(Nx, Ny, Nz, data, vec.TV_ROSEM_apu, epps, 4, w_vec);
				vec.APLS_ROSEM_apu = ROSEM(vec.TV_ROSEM_apu, Summ, vec.TV_ROSEM_rhs, w_vec.lambda_ROSEM, iter);
			}
			if (MethodList.TGV) {
				array dU = TGV(vec.TGV_ROSEM_apu, Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
				vec.TGV_ROSEM_apu = ROSEM(vec.TGV_ROSEM_apu, Summ, vec.TGV_ROSEM_rhs, w_vec.lambda_ROSEM, iter);
			}
			custom_ROSEM = ROSEM(custom_ROSEM, Summ, custom_ROSEM_rhs, w_vec.lambda_ROSEM, iter);
		}

		if (MethodList.RBI)
			vec.RBI_apu = RBI(vec.RBI_apu, Summ, vec.RBI_rhs, w_vec.D, 0.f, af::constant(0.f, 1, 1));
		if (MethodList.RBIMAP) {
			if (MethodList.MRP) {
				array dU = MRP(vec.MRP_RBI_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.med_no_norm, im_dim);
				vec.MRP_RBI_apu = RBI(vec.MRP_RBI_apu, Summ, vec.MRP_RBI_rhs, w_vec.D, beta.MRP_RBI, dU);
			}
			if (MethodList.Quad) {
				array dU = Quadratic_prior(vec.Quad_RBI_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, w_vec.tr_offsets, w_vec.weights_quad, im_dim);
				vec.Quad_RBI_apu = RBI(vec.Quad_RBI_apu, Summ, vec.RBI_rhs, w_vec.D, beta.Quad_RBI, dU);
			}
			if (MethodList.L) {
				array dU = L_filter(vec.L_RBI_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L, w_vec.med_no_norm, im_dim);
				vec.L_RBI_apu = RBI(vec.L_RBI_apu, Summ, vec.L_RBI_rhs, w_vec.D, beta.L_RBI, dU);
			}
			if (MethodList.FMH) {
				array dU = FMH(vec.FMH_RBI_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets, w_vec.fmh_weights, w_vec.med_no_norm,
					w_vec.alku_fmh, im_dim);
				vec.FMH_RBI_apu = RBI(vec.FMH_RBI_apu, Summ, vec.FMH_RBI_rhs, w_vec.D, beta.FMH_RBI, dU);
			}
			if (MethodList.WeightedMean) {
				array dU = Weighted_mean(vec.Weighted_RBI_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.weighted_weights,
					w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
				vec.Weighted_RBI_apu = RBI(vec.Weighted_RBI_apu, Summ, vec.Weighted_RBI_rhs, w_vec.D, beta.Weighted_RBI, dU);
			}
			if (MethodList.TV) {
				array dU = TVprior(Nx, Ny, Nz, data, vec.TV_RBI_apu, epps, data.TVtype, w_vec);
				vec.TV_RBI_apu = RBI(vec.TV_RBI_apu, Summ, vec.TV_RBI_rhs, w_vec.D, beta.TV_RBI, dU);
			}
			if (MethodList.AD) {
				array dU = AD(vec.AD_RBI_apu, Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, w_vec.DiffusionType, w_vec.med_no_norm);
				vec.AD_RBI_apu = RBI(vec.AD_RBI_apu, Summ, vec.AD_RBI_rhs, w_vec.D, beta.AD_RBI, dU);
			}
			if (MethodList.APLS) {
				array dU = TVprior(Nx, Ny, Nz, data, vec.TV_RBI_apu, epps, 4, w_vec);
				vec.APLS_RBI_apu = RBI(vec.TV_RBI_apu, Summ, vec.TV_RBI_rhs, w_vec.D, beta.APLS_RBI, dU);
			}
			if (MethodList.TGV) {
				array dU = TGV(vec.TGV_RBI_apu, Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
				vec.TGV_RBI_apu = RBI(vec.TGV_RBI_apu, Summ, vec.TGV_RBI_rhs, w_vec.D, beta.TGV_RBI, dU);
			}
			custom_RBI = RBI(custom_RBI, Summ, custom_RBI_rhs, w_vec.D, beta_custom_RBI, grad_RBI);
		}

		if (MethodList.COSEM || MethodList.ECOSEM)
			vec.COSEM_apu = COSEM(vec.COSEM_apu, vec.C_co, w_vec.D);
		if (MethodList.OSLCOSEM > 0u) {
			if (MethodList.MRP) {
				mexPrintf("OSLCOSEM = %d\n", MethodList.OSLCOSEM);
				array dU = MRP(vec.MRP_COSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.med_no_norm, im_dim);
				vec.MRP_COSEM_apu = OSL_COSEM(vec.MRP_COSEM_apu, vec.C_osl, w_vec.D, w_vec.h_ACOSEM, MethodList.OSLCOSEM, dU, beta.MRP_COSEM);
				if (MethodList.OSLCOSEM == 1) {
					MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.MRP_COSEM_apu, vec.C_co,
						vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length);
					w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
					vec.MRP_COSEM_apu = batchFunc(vec.MRP_COSEM_apu, uu / w_vec.ACOSEM_rhs, batchMul);
					vec.MRP_COSEM_apu(vec.MRP_COSEM_apu < 0.f) = epps;
				}
			}
			if (MethodList.Quad) {
				array dU = Quadratic_prior(vec.Quad_COSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, w_vec.tr_offsets, w_vec.weights_quad, im_dim);
				vec.Quad_COSEM_apu = OSL_COSEM(vec.Quad_COSEM_apu, vec.C_osl, w_vec.D, w_vec.h_ACOSEM, MethodList.OSLCOSEM, dU, beta.Quad_COSEM);
				if (MethodList.OSLCOSEM == 1u) {
					MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.MRP_COSEM_apu, vec.C_co,
						vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length);
					w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
					vec.Quad_COSEM_apu = batchFunc(vec.Quad_COSEM_apu, uu / w_vec.ACOSEM_rhs, batchMul);
					vec.Quad_COSEM_apu(vec.Quad_COSEM_apu < 0.f) = epps;
				}
			}
			if (MethodList.L) {
				array dU = L_filter(vec.L_COSEM_apu, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L, w_vec.med_no_norm, im_dim);
				vec.L_COSEM_apu = OSL_COSEM(vec.L_COSEM_apu, vec.C_osl, w_vec.D, w_vec.h_ACOSEM, MethodList.OSLCOSEM, dU, beta.L_COSEM);
				if (MethodList.OSLCOSEM == 1u) {
					MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.MRP_COSEM_apu, vec.C_co,
						vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length);
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
					MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.MRP_COSEM_apu, vec.C_co,
						vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length);
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
					MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.MRP_COSEM_apu, vec.C_co,
						vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length);
					w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
					vec.Weighted_COSEM_apu = batchFunc(vec.Weighted_COSEM_apu, uu / w_vec.ACOSEM_rhs, batchMul);
					vec.Weighted_COSEM_apu(vec.Weighted_COSEM_apu < 0.f) = epps;
				}
			}
			if (MethodList.TV) {
				array dU = TVprior(Nx, Ny, Nz, data, vec.TV_COSEM_apu, epps, data.TVtype, w_vec);
				vec.TV_COSEM_apu = OSL_COSEM(vec.TV_COSEM_apu, vec.C_osl, w_vec.D, w_vec.h_ACOSEM, MethodList.OSLCOSEM, dU, beta.TV_COSEM);
				if (MethodList.OSLCOSEM == 1u) {
					MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.MRP_COSEM_apu, vec.C_co,
						vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length);
					w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
					vec.TV_COSEM_apu = batchFunc(vec.TV_COSEM_apu, uu / w_vec.ACOSEM_rhs, batchMul);
					vec.TV_COSEM_apu(vec.TV_COSEM_apu < 0.f) = epps;
				}
			}
			if (MethodList.AD) {
				array dU = AD(vec.AD_COSEM_apu, Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, w_vec.DiffusionType, w_vec.med_no_norm);
				vec.AD_COSEM_apu = OSL_COSEM(vec.AD_COSEM_apu, vec.C_osl, w_vec.D, w_vec.h_ACOSEM, MethodList.OSLCOSEM, dU, beta.AD_COSEM);
				if (MethodList.OSLCOSEM == 1u) {
					MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.MRP_COSEM_apu, vec.C_co,
						vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length);
					w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
					vec.AD_COSEM_apu = batchFunc(vec.AD_COSEM_apu, uu / w_vec.ACOSEM_rhs, batchMul);
					vec.AD_COSEM_apu(vec.AD_COSEM_apu < 0.f) = epps;
				}
			}
			if (MethodList.APLS) {
				array dU = TVprior(Nx, Ny, Nz, data, vec.TV_COSEM_apu, epps, 4, w_vec);
				vec.APLS_COSEM_apu = OSL_COSEM(vec.APLS_COSEM_apu, vec.C_osl, w_vec.D, w_vec.h_ACOSEM, MethodList.OSLCOSEM, dU, beta.APLS_COSEM);
				if (MethodList.OSLCOSEM == 1u) {
					MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.MRP_COSEM_apu, vec.C_co,
						vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length);
					w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
					vec.APLS_COSEM_apu = batchFunc(vec.APLS_COSEM_apu, uu / w_vec.ACOSEM_rhs, batchMul);
					vec.APLS_COSEM_apu(vec.APLS_COSEM_apu < 0.f) = epps;
				}
			}
			if (MethodList.TGV) {
				array dU = TGV(vec.TGV_COSEM_apu, Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
				vec.TGV_COSEM_apu = OSL_COSEM(vec.TGV_COSEM_apu, vec.C_osl, w_vec.D, w_vec.h_ACOSEM, MethodList.OSLCOSEM, dU, beta.TGV_COSEM);
				if (MethodList.OSLCOSEM == 1u) {
					MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.MRP_COSEM_apu, vec.C_co,
						vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length);
					w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
					vec.TGV_COSEM_apu = batchFunc(vec.TGV_COSEM_apu, uu / w_vec.ACOSEM_rhs, batchMul);
					vec.TGV_COSEM_apu(vec.TGV_COSEM_apu < 0.f) = epps;
				}
			}
			custom_COSEM = OSL_COSEM(custom_COSEM, vec.C_osl, w_vec.D, w_vec.h_ACOSEM, MethodList.OSLCOSEM, grad_COSEM, beta_custom_COSEM);
			if (MethodList.OSLCOSEM == 1u) {
				MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.MRP_COSEM_apu, vec.C_co,
					vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length);
				w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
				custom_COSEM = batchFunc(custom_COSEM, uu / w_vec.ACOSEM_rhs, batchMul);
				custom_COSEM(custom_COSEM < 0.f) = epps;
			}
		}

		if (MethodList.ECOSEM)
			vec.ECOSEM_apu = ECOSEM(vec.ECOSEM_apu, w_vec.D, vec.OSEM_apu, vec.COSEM_apu, epps);

		if (MethodList.ACOSEM) {
			vec.ACOSEM_apu = ACOSEM(vec.ACOSEM_apu, vec.C_aco, w_vec.D, w_vec.h_ACOSEM);
			MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, rekot, d_Sino, koko_l, vec.ACOSEM_apu, vec.C_co,
				vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length);
			w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
			vec.ACOSEM_apu = batchFunc(vec.ACOSEM_apu, uu / w_vec.ACOSEM_rhs, batchMul);
			vec.ACOSEM_apu(vec.ACOSEM_apu < 0.f) = epps;
		}

		if (MethodList.DRAMA)
			vec.DRAMA_apu = DRAMA(vec.DRAMA_apu, Summ, vec.DRAMA_rhs, w_vec.lambda_DRAMA, iter, osa_iter, subsets);

		clFinish(af_queue);

		w_vec.epsilon_mramla.unlock();
	}

	//device_to_host_cell(ArrayList, MethodList, vec, oo, cell);

	device_to_host_cell_custom(ArrayList, MethodList, vec, oo, cell);


	if (MethodList.OSLOSEM)
		custom_OSEM.host(ele_custom_osem);
	if (MethodList.OSLMLEM)
		custom_MLEM.host(ele_custom_mlem);
	if (MethodList.BSREM)
		custom_BSREM.host(ele_custom_bsrem);
	if (MethodList.MBSREM)
		custom_MBSREM.host(ele_custom_mbsrem);
	if (MethodList.ROSEMMAP)
		custom_ROSEM.host(ele_custom_rosem);
	if (MethodList.RBIMAP)
		custom_RBI.host(ele_custom_rbi);
	if (MethodList.OSLCOSEM > 0u)
		custom_COSEM.host(ele_custom_cosem);

	oo = size_rekot;

	if (MethodList.OSLMLEM)
		mxSetCell(cell, oo, custom_mlem);
	oo++;
	if (MethodList.OSLOSEM)
		mxSetCell(cell, oo, custom_osem);
	oo++;
	if (MethodList.BSREM)
		mxSetCell(cell, oo, custom_bsrem);
	oo++;
	if (MethodList.MBSREM)
		mxSetCell(cell, oo, custom_mbsrem);
	oo++;
	if (MethodList.ROSEMMAP)
		mxSetCell(cell, oo, custom_rosem);
	oo++;
	if (MethodList.RBIMAP)
		mxSetCell(cell, oo, custom_rbi);
	oo++;
	if (MethodList.OSLCOSEM > 0u)
		mxSetCell(cell, oo, custom_cosem);
	oo++;

	if (MethodList.COSEM || MethodList.ECOSEM) {
		mwSize dime[2] = { im_dim, subsets };
		mxArray* apu2 = mxCreateNumericArray(2, dime, mxSINGLE_CLASS, mxREAL);
		float* apu = (float*)mxGetData(apu2);
		apu = vec.C_co.host<float>();
		mxSetCell(cell, oo, apu2);
	}
	oo++;
	if (MethodList.ACOSEM) {
		mwSize dime[2] = { im_dim, subsets };
		mxArray* apu2 = mxCreateNumericArray(2, dime, mxSINGLE_CLASS, mxREAL);
		float* apu = (float*)mxGetData(apu2);
		apu = vec.C_aco.host<float>();
		mxSetCell(cell, oo, apu2);
	}
	oo++;
	if (MethodList.OSLCOSEM > 0) {
		mwSize dime[2] = { im_dim, subsets };
		mxArray* apu2 = mxCreateNumericArray(2, dime, mxSINGLE_CLASS, mxREAL);
		float* apu = (float*)mxGetData(apu2);
		apu = vec.C_osl.host<float>();
		mxSetCell(cell, oo, apu2);
	}
	oo++;

	if (((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBI || MethodList.RBIMAP) && w_vec.MBSREM_prepass ||
		MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0) && osa_iter == 0 && iter == 0 && tt == 0) {
		mwSize dime = im_dim;
		mxArray* apu2 = mxCreateNumericArray(1, &dime, mxSINGLE_CLASS, mxREAL);
		float* apu = (float*)mxGetData(apu2);
		apu = w_vec.D.host<float>();
		mxSetCell(cell, oo, apu2);
	}


	unlock_AF_im_vectors(vec, MethodList, true, mlem_bool, osem_bool, 0);

	// Clear OpenCL buffers
	//clReleaseMemObject(d_Nx);
	//clReleaseMemObject(d_Ny);
	//clReleaseMemObject(d_Nz);
	//clReleaseMemObject(d_dz);
	//clReleaseMemObject(d_dx);
	//clReleaseMemObject(d_dy);
	//clReleaseMemObject(d_bz);
	//clReleaseMemObject(d_bx);
	//clReleaseMemObject(d_by);
	//clReleaseMemObject(d_bzb);
	//clReleaseMemObject(d_maxxx);
	//clReleaseMemObject(d_maxyy);
	//clReleaseMemObject(d_zmax);
	//clReleaseMemObject(d_NSlices);
	clReleaseMemObject(d_z);
	clReleaseMemObject(d_x);
	clReleaseMemObject(d_y);
	clReleaseMemObject(d_atten);
	clReleaseMemObject(d_xcenter);
	clReleaseMemObject(d_ycenter);
	//clReleaseMemObject(d_epps);
	//clReleaseMemObject(d_size_x);
	//clReleaseMemObject(d_TotSinos);
	//clReleaseMemObject(d_attenuation_correction);
	//clReleaseMemObject(d_N);
	//clReleaseMemObject(d_raw);
	//clReleaseMemObject(d_det_per_ring);
	clReleaseMemObject(d_pseudos);
	//clReleaseMemObject(d_pRows);
	//clReleaseMemObject(d_h);
	//clReleaseMemObject(d_epsilon_mramla);


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
	if ((MethodList.MLEM || MethodList.OSLMLEM) && osa_iter == 0u) {
		status = clReleaseKernel(kernel_ml);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			mexPrintf("Failed to release kernel\n");
		}
	}
	if ((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBI || MethodList.RBIMAP) && w_vec.MBSREM_prepass ||
		MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0u) {
		status = clReleaseKernel(kernel_mramla);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			mexPrintf("Failed to release kernel\n");
		}
	}
}
