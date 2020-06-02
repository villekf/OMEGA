/**************************************************************************
* Matrix free computations for OMEGA.
* In this file the OpenCL buffers are created, calls to other necessary 
* functions are made and the OpenCL kernels are launched. This file 
* contains the code for the matrix-free reconstructions in OMEGA using the
* implementation 2.
* Unlike the non-OpenCL versions, this one uses (32-bit) floats and thus
* can be slightly more inaccurate.
*
* Copyright (C) 2020 Ville-Veikko Wettenhovi
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
#include "AF_opencl_functions.hpp"

// Use ArrayFire namespace for convenience
using namespace af;

// Main reconstruction function
void reconstruction_AF_matrixfree(const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, const mxArray* Sin,
	const mxArray* sc_ra, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const uint32_t Niter, const mxArray* options, const float dx,
	const float dy, const float dz, const float bx, const float by, const float bz, const float bzb, const float maxxx, const float maxyy, const float zmax,
	const float NSlices, const uint32_t* pituus, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t size_x, const uint32_t TotSinos, 
	mxArray* cell, const mwSize* dimmi, const bool verbose, const uint32_t randoms_correction, const uint32_t attenuation_correction,
	const uint32_t normalization, const float* atten, const size_t size_atten, const float* norm, const size_t size_norm, const uint32_t subsets,
	const float epps, const char* k_path, const uint32_t Nt, const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L, 
	const uint8_t raw, const size_t size_z, const bool osem_bool, const char* fileName, const bool use_psf, const float tube_width, 
	const float crystal_size_z, const float* x_center, const float* y_center, const float* z_center, const size_t size_center_x, const size_t size_center_y, 
	const size_t size_of_x, const size_t size_center_z, const uint32_t projector_type, const char* header_directory, const bool precompute, 
	const uint32_t device, const int32_t dec, const uint16_t n_rays, const uint16_t n_rays3D, const float cr_pz, const bool use_64bit_atomics, uint32_t n_rekos,
	const uint32_t n_rekos_mlem, const uint8_t* reko_type, const uint8_t* reko_type_mlem, const float global_factor, const float bmin, const float bmax, 
	const float Vmax, const float* V, const size_t size_V, const float* gaussian, const size_t size_gauss) {

	// Number of voxels
	const uint32_t Nxy = Nx * Ny;
	const uint32_t im_dim = Nxy * Nz;
	const cl_ulong st = 0ULL;
	const cl_uchar fp = 0;

	// Distance between rays in multi-ray Siddon
	const float dc_z = cr_pz / static_cast<float>(n_rays3D + 1);

	bool break_iter = false;

	uint32_t t0 = 0u;
	uint32_t iter0 = 0u;
	uint32_t osa_iter0 = 0u;

	// Hard-coded local size
	const uint64_t local_size = 64ULL;

	const cl_float zerof = 0.f;
	cl_uint kernelInd_MRAMLA = 0U;

	bool atomic_64bit = use_64bit_atomics;
	cl_uchar compute_norm_matrix = 1u;
	cl_float mem_portions;
	if (raw == 1u)
		mem_portions = 0.1f;
	else
		mem_portions = 0.2f;
	cl_float image_bytes = static_cast<cl_float>(Nx * Ny * Nz) * 8.f;

	uint32_t oo = 0u;
	size_t ll = 0ULL;

	// Create output structs
	matlabArrays ArrayList;
	// Create a struct containing the reconstruction methods used
	RecMethods MethodList;
	// Same as above, but as cl_ variables
	RecMethodsOpenCL MethodListOpenCL;

	kernelStruct OpenCLStruct;

	// Obtain the reconstruction methods used
	get_rec_methods(options, MethodList);
	OpenCLRecMethods(MethodList, MethodListOpenCL);

	// Create the OpenCL context and command queue and assign the device
	cl_context af_context = afcl::getContext();
	cl_device_id af_device_id = afcl::getDeviceId();
	cl_command_queue af_queue = afcl::getQueue();

	cl_program program_os = NULL;
	cl_program program_ml = NULL;
	cl_program program_mbsrem = NULL;

	OpenCLStruct.af_queue = &af_queue;

	cl_int status = CL_SUCCESS;

	cl_ulong mem;
	cl_ulong mem_loc;
	status = clGetDeviceInfo(af_device_id, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &mem, NULL);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}
	if ((static_cast<cl_float>(mem) * mem_portions) > image_bytes && !MethodList.CUSTOM)
		compute_norm_matrix = 0u;

	status = clGetDeviceInfo(af_device_id, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &mem_loc, NULL);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}

	// Create the MATLAB output arrays
	create_matlab_output(ArrayList, dimmi, MethodList, 4);

	// Initial value
	array x00(Nx * Ny * Nz, (float*)mxGetData(mxGetField(options, 0, "x0")), afHost);

	// For custom prior
	if (MethodList.CUSTOM) {
		osa_iter0 = (uint32_t)mxGetScalar(mxGetField(options, 0, "osa_iter"));
		iter0 = (uint32_t)mxGetScalar(mxGetField(options, 0, "iter"));
		t0 = (uint32_t)mxGetScalar(mxGetField(options, 0, "tt"));
	}

	array pj3, apu_sum, E, apu_sum_mlem;

	if (MethodList.MRAMLA || MethodList.MBSREM)
		E = constant(1.f, koko, 1);
	else
		E = constant(0.f, 1, 1);

	// Are ML-methods used?
	bool mlem_bool = (MethodList.MLEM || MethodList.OSLMLEM) ? true : false;

	// Number of measurements at each subset
	std::vector<size_t> length(subsets);

	for (uint32_t kk = 0; kk < subsets; kk++)
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
	form_data_variables(vec, beta, w_vec, options, Nx, Ny, Nz, Niter, x00, im_dim, koko, MethodList, data, subsets, osa_iter0, use_psf);

	// Power factor for ACOSEM
	w_vec.h_ACOSEM_2 = 1.f / w_vec.h_ACOSEM;


	// Adjust the number of reconstruction methods and create the output vector containing all the estimates
	// E.g. using ECOSEM, the OSEM and COSEM estimates are needed even if they were not selected
	uint32_t n_rekos2 = n_rekos;
	if (osem_bool) {
		if (MethodList.ECOSEM) {
			if (!MethodList.OSEM && !MethodList.COSEM) {
				n_rekos++;
				n_rekos2 += 2u;
			}
			else if (!MethodList.OSEM || !MethodList.COSEM)
				n_rekos2++;
			else if (MethodList.OSEM && MethodList.COSEM)
				n_rekos--;
		}
		vec.im_os = constant(0.f, im_dim * n_rekos2, 1);
		for (int kk = 0; kk < n_rekos2; kk++) {
			vec.im_os(seq(kk* im_dim, (kk + 1)* im_dim - 1)) = x00;
		}
	}
	
	if (mlem_bool) {
		vec.im_mlem = constant(0.f, im_dim * n_rekos_mlem, 1);
		for (int kk = 0; kk < n_rekos_mlem; kk++) {
			vec.im_mlem(seq(kk * im_dim, (kk + 1) * im_dim - 1)) = x00;
		}
	}

	float* scat = nullptr;
	const uint32_t scatter = static_cast<uint32_t>((bool)mxGetScalar(mxGetField(options, 0, "scatter")));
	size_t size_scat = 1ULL;
	scat = (float*)mxGetData(mxGetCell(mxGetField(options, 0, "ScatterC"), 0));
	if (scatter == 1U) {
		size_scat = mxGetNumberOfElements(mxGetCell(mxGetField(options, 0, "ScatterC"), 0));
	}

	status = createProgram(verbose, k_path, af_context, af_device_id, fileName, program_os, program_ml, program_mbsrem, atomic_64bit, device, header_directory,
		projector_type, crystal_size_z, precompute, raw, attenuation_correction, normalization, dec, local_size, n_rays, n_rays3D, false, MethodList, osem_bool, 
		mlem_bool, n_rekos2, n_rekos_mlem, w_vec, osa_iter0, cr_pz, dx, use_psf, scatter, randoms_correction);
	if (status != CL_SUCCESS) {
		std::cerr << "Error while creating program" << std::endl;
		return;
	}

	std::vector<array> Summ;
	array Summ_mlem;

	// Normalization constant
	// Save the constants if there was enough memory
	// 64-bit atomic operations require unsigned 64-bit integers
	if (compute_norm_matrix == 0u) {
		if (atomic_64bit && !w_vec.MBSREM_prepass)
			Summ.assign(subsets, constant(0ULL, im_dim, 1, u64));
		else
			Summ.assign(subsets, constant(0.f, im_dim, 1));
	}
	else {
		if (atomic_64bit)
			Summ.assign(1ULL, constant(0ULL, im_dim, 1, u64));
		else
			Summ.assign(1ULL, constant(0.f, im_dim, 1));
	}

	// Create the kernels
	cl_kernel kernel_ml = NULL, kernel = NULL, kernel_mramla = NULL;

	status = createKernels(kernel_ml, kernel, kernel_mramla, OpenCLStruct.kernelNLM, osem_bool, program_os, program_ml, program_mbsrem, MethodList, w_vec, projector_type, mlem_bool, precompute,
		n_rays, n_rays3D);
	if (status != CL_SUCCESS) {
		mexPrintf("Failed to create kernels\n");
		clReleaseProgram(program_os);
		clReleaseProgram(program_ml);
		clReleaseProgram(program_mbsrem);
		clReleaseKernel(kernel);
		clReleaseKernel(kernel_ml);
		clReleaseKernel(kernel_mramla);
		clReleaseKernel(OpenCLStruct.kernelNLM);
		return;
	}
	else {
		mexPrintf("OpenCL kernels successfully created\n");
		mexEvalString("pause(.0001);");
	}

	// Create and write buffers
	cl_mem d_x, d_y, d_z, d_pseudos, d_atten, d_xcenter, d_ycenter, d_zcenter, d_norm_mlem, d_reko_type, d_reko_type_mlem, d_V, d_scat_mlem;
	cl_mem* d_Summ, *d_Summ_mlem;
	cl_mem d_lor_mlem, d_L_mlem, d_zindex_mlem, d_xyindex_mlem, d_Sino_mlem, d_sc_ra_mlem;

	std::vector<cl_mem> d_lor(subsets);
	std::vector<cl_mem> d_L(subsets);
	std::vector<cl_mem> d_zindex(subsets);
	std::vector<cl_mem> d_xyindex(subsets);
	std::vector<cl_mem> d_Sino(subsets);
	std::vector<cl_mem> d_sc_ra(subsets);
	std::vector<cl_mem> d_norm(subsets);
	std::vector<cl_mem> d_scat(subsets);

	float *apu = (float*)mxGetData(mxGetCell(Sin, 0));

	status = createAndWriteBuffers(d_x, d_y, d_z, d_lor, d_L, d_zindex, d_xyindex, d_Sino, d_sc_ra, size_x, size_z, TotSinos, size_atten, size_norm, size_scat, prows, 
		length, x, y, z_det, xy_index, z_index, lor1, L, apu, raw, af_context, subsets, pituus, atten, norm, scat, pseudos, V, af_queue, d_atten, d_norm, d_scat, d_pseudos, d_V, 
		d_xcenter, d_ycenter, d_zcenter, x_center, y_center, z_center, size_center_x, size_center_y, size_center_z, size_of_x, size_V, atomic_64bit, randoms_correction, 
		sc_ra, precompute, d_lor_mlem, d_L_mlem, d_zindex_mlem, d_xyindex_mlem, d_Sino_mlem, d_sc_ra_mlem, d_reko_type, d_reko_type_mlem, osem_bool, mlem_bool, koko,
		reko_type, reko_type_mlem, n_rekos, n_rekos_mlem, d_norm_mlem, d_scat_mlem);
	if (status != CL_SUCCESS) {
		clReleaseProgram(program_os);
		clReleaseProgram(program_ml);
		clReleaseProgram(program_mbsrem);
		clReleaseKernel(kernel);
		clReleaseKernel(kernel_ml);
		clReleaseKernel(kernel_mramla);
		clReleaseKernel(OpenCLStruct.kernelNLM);
		mexPrintf("Buffer creation failed\n");
		return;
	}

	array g(size_gauss, gaussian, afHost);
	if (use_psf) {
		g = moddims(g, w_vec.g_dim_x * 2u + 1u, w_vec.g_dim_y * 2u + 1u, w_vec.g_dim_z * 2u + 1u);
	}
	uint32_t deblur_iterations = 0U;
	if (use_psf && w_vec.deconvolution) {
		deblur_iterations = (uint32_t)mxGetScalar(mxGetField(options, 0, "deblur_iterations"));
	}

	cl_uint kernelInd_OSEM = 0U;
	cl_uint kernelInd_MLEM = 0U;

	if (osem_bool) {
		// Set the kernel parameters that do not change
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &global_factor);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &epps);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint32_t), &im_dim);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint32_t), &Nx);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint32_t), &Ny);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint32_t), &Nz);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &dz);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &dx);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &dy);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &bz);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &bx);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &by);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &bzb);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &maxxx);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &maxyy);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &zmax);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &NSlices);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint32_t), &size_x);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint16_t), &TotSinos);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint32_t), &det_per_ring);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint32_t), &prows);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint32_t), &Nxy);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(cl_uchar), &fp);
		if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1))) {
			clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &tube_width);
			clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &crystal_size_z);
			clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &bmin);
			clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &bmax);
			clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &Vmax);
		}
		else if (projector_type == 1u && !precompute) {
			clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &dc_z);
			clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint16_t), &n_rays);
		}
	}

	if (mlem_bool) {
		// Set the kernel parameters that do not change
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &global_factor);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &epps);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &im_dim);
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &Nx);
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &Ny);
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &Nz);
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &dz);
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &dx);
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &dy);
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &bz);
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &bx);
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &by);
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &bzb);
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &maxxx);
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &maxyy);
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &zmax);
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &NSlices);
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &size_x);
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint16_t), &TotSinos);
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &det_per_ring);
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &prows);
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &Nxy);
		status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(cl_uchar), &fp);
		if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1))) {
			status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &tube_width);
			status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &crystal_size_z);
			status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &bmin);
			status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &bmax);
			status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &Vmax);
		}
		else if (projector_type == 1u && !precompute) {
			status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &dc_z);
			status = clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint16_t), &n_rays);
		}
	}

	if (mlem_bool) {
		if (atomic_64bit)
			Summ_mlem = constant(0ULL, im_dim, 1, u64);
		else
			Summ_mlem = constant(0.f, im_dim, 1);
	}

	// Loop through each time-step
	for (uint32_t tt = t0; tt < Nt; tt++) {
	// Compute the prepass phase for MRAMLA, MBSREM, RBI, COSEM, ACOSEM or ECOSEM if applicable
	if (((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI) && w_vec.MBSREM_prepass ||
		MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0) && (!MethodList.CUSTOM || osa_iter0 == 0u)) {

		// Set the kernel parameters that do not change
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &global_factor);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &epps);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint32_t), &im_dim);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint32_t), &Nx);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint32_t), &Ny);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint32_t), &Nz);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &dz);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &dx);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &dy);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &bz);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &bx);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &by);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &bzb);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &maxxx);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &maxyy);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &zmax);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &NSlices);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint32_t), &size_x);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint16_t), &TotSinos);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint32_t), &det_per_ring);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint32_t), &prows);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint32_t), &Nxy);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(cl_uchar), &fp);
		if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1))) {
			clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &tube_width);
			clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &crystal_size_z);
			clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &bmin);
			clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &bmax);
			clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &Vmax);
		}
		else if (projector_type == 1u && !precompute) {
			clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &dc_z);
			clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint16_t), &n_rays);
		}
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &w_vec.epsilon_mramla);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(cl_mem), &d_atten);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(cl_mem), &d_pseudos);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(cl_mem), &d_x);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(cl_mem), &d_y);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(cl_mem), &d_z);
		if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1))) {
			clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(cl_mem), &d_xcenter);
			clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(cl_mem), &d_ycenter);
			clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(cl_mem), &d_zcenter);
			clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(cl_mem), &d_V);
		}
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(cl_mem), &d_reko_type);

		uint32_t alku = 0u;

		if ((MethodList.MRAMLA || MethodList.MBSREM) && w_vec.MBSREM_prepass)
			w_vec.Amin = constant(0.f, koko, 1);

		mexPrintf("n_rekos2 = %d\n", n_rekos2);
		mexEvalString("pause(.0001);");

		// Run the prepass phase
		MRAMLA_prepass(subsets, im_dim, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, w_vec, Summ, d_Sino, koko, x00, vec.C_co,
			vec.C_aco, vec.C_osl, alku, kernel_mramla, d_L, raw, MethodListOpenCL, length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E, 
			d_norm, d_scat, use_psf, g, Nx, Ny, Nz, epps);


		if ((MethodList.MRAMLA || MethodList.MBSREM) && tt == 0) {
			pj3 = w_vec.D / static_cast<float>(subsets);
		}

		if (verbose) {
			mexPrintf("MRAMLA & COSEM prepass completed\n");
			mexEvalString("pause(.0001);");
		}
	}


		cl_uint kernelInd_OSEMTIter = kernelInd_OSEM;
		cl_uint kernelInd_MLEMT = kernelInd_MLEM;

		cl_uchar no_norm = 0u;
		cl_uchar no_norm_mlem = 0u;
		if (tt == 0u  && compute_norm_matrix == 0u) {
			if (w_vec.MBSREM_prepass)
				no_norm = 1u;
			if (osem_bool)
				no_norm_mlem = 1u;
		}
		else if (tt > 0u) {
			if (osem_bool && compute_norm_matrix == 0u)
				no_norm = 1u;
			if (mlem_bool)
				no_norm_mlem = 1u;
		}

		// Load the measurement and randoms data from the cell arrays
		if (tt > 0u) {
			float* apu = (float*)mxGetData(mxGetCell(Sin, tt));
			if (osem_bool) {
				for (uint32_t kk = 0u; kk < subsets; kk++) {
					status = clEnqueueWriteBuffer(af_queue, d_Sino[kk], CL_TRUE, 0, sizeof(float) * length[kk], &apu[pituus[kk]], 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
					if (randoms_correction) {
						float* ra_apu = (float*)mxGetData(mxGetCell(sc_ra, tt));
						status = clEnqueueWriteBuffer(af_queue, d_sc_ra[kk], CL_TRUE, 0, sizeof(float) * length[kk], &ra_apu[pituus[kk]], 0, NULL, NULL);
					}
					else
						status = clEnqueueFillBuffer(af_queue, d_sc_ra[kk], &zerof, sizeof(cl_float), 0, sizeof(cl_float), 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
					if (scatter == 1u) {
						scat = (float*)mxGetData(mxGetCell(mxGetField(options, 0, "ScatterC"), tt));
						status = clEnqueueWriteBuffer(af_queue, d_scat[kk], CL_TRUE, 0, sizeof(float) * length[kk], &scat[pituus[kk]], 0, NULL, NULL);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
					}
				}
				vec.im_os = constant(0.f, im_dim * n_rekos2, 1);
				for (int kk = 0; kk < n_rekos2; kk++) {
					vec.im_os(seq(kk * im_dim, (kk + 1) * im_dim - 1)) = x00;
				}
			}
			if (mlem_bool) {
				status = clEnqueueWriteBuffer(af_queue, d_Sino_mlem, CL_TRUE, 0, sizeof(float) * koko, apu, 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				if (randoms_correction) {
					float* ra_apu = (float*)mxGetData(mxGetCell(sc_ra, tt));
					status = clEnqueueWriteBuffer(af_queue, d_sc_ra_mlem, CL_TRUE, 0, sizeof(float) * koko, ra_apu, 0, NULL, NULL);
				}
				else
					status = clEnqueueFillBuffer(af_queue, d_sc_ra_mlem, &zerof, sizeof(cl_float), 0, sizeof(cl_float), 0, NULL, NULL);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				if (scatter == 1u) {
					scat = (float*)mxGetData(mxGetCell(mxGetField(options, 0, "ScatterC"), tt));
					status = clEnqueueWriteBuffer(af_queue, d_scat_mlem, CL_TRUE, 0, sizeof(float) * koko, scat, 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				vec.im_mlem = constant(0.f, im_dim * n_rekos_mlem, 1);
				for (int kk = 0; kk < n_rekos_mlem; kk++) {
					vec.im_mlem(seq(kk * im_dim, (kk + 1) * im_dim - 1)) = x00;
				}
			}
		}

		// Compute values needed for MBSREM and MRAMLA
		if (MethodList.MBSREM || MethodList.MRAMLA) {
			array Sino = array(koko, (float*)mxGetData(mxGetCell(Sin, tt)), afHost);
			array rand;
			if (randoms_correction)
				rand = array(koko, (float*)mxGetData(mxGetCell(sc_ra, tt)), afHost);
			if (w_vec.U == 0.f) {
				w_vec.U = max<float>(Sino / w_vec.Amin);
			}
			w_vec.epsilon_mramla = MBSREM_epsilon(Sino, epps, randoms_correction, rand, E);
		}
		if (osem_bool) {
			clSetKernelArg(kernel, kernelInd_OSEMTIter++, sizeof(float), &w_vec.epsilon_mramla);
			clSetKernelArg(kernel, kernelInd_OSEMTIter++, sizeof(cl_mem), &d_atten);
			clSetKernelArg(kernel, kernelInd_OSEMTIter++, sizeof(cl_mem), &d_pseudos);
			clSetKernelArg(kernel, kernelInd_OSEMTIter++, sizeof(cl_mem), &d_x);
			clSetKernelArg(kernel, kernelInd_OSEMTIter++, sizeof(cl_mem), &d_y);
			clSetKernelArg(kernel, kernelInd_OSEMTIter++, sizeof(cl_mem), &d_z);
			if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1))) {
				clSetKernelArg(kernel, kernelInd_OSEMTIter++, sizeof(cl_mem), &d_xcenter);
				clSetKernelArg(kernel, kernelInd_OSEMTIter++, sizeof(cl_mem), &d_ycenter);
				clSetKernelArg(kernel, kernelInd_OSEMTIter++, sizeof(cl_mem), &d_zcenter);
				clSetKernelArg(kernel, kernelInd_OSEMTIter++, sizeof(cl_mem), &d_V);
			}
			clSetKernelArg(kernel, kernelInd_OSEMTIter++, sizeof(cl_mem), &d_reko_type);
		}
		if (mlem_bool) {
			status = clSetKernelArg(kernel_ml, kernelInd_MLEMT++, sizeof(float), &w_vec.epsilon_mramla);
			status = clSetKernelArg(kernel_ml, kernelInd_MLEMT++, sizeof(cl_mem), &d_atten);
			status = clSetKernelArg(kernel_ml, kernelInd_MLEMT++, sizeof(cl_mem), &d_pseudos);
			status = clSetKernelArg(kernel_ml, kernelInd_MLEMT++, sizeof(cl_mem), &d_x);
			status = clSetKernelArg(kernel_ml, kernelInd_MLEMT++, sizeof(cl_mem), &d_y);
			status = clSetKernelArg(kernel_ml, kernelInd_MLEMT++, sizeof(cl_mem), &d_z);
			if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1))) {
				status = clSetKernelArg(kernel_ml, kernelInd_MLEMT++, sizeof(cl_mem), &d_xcenter);
				status = clSetKernelArg(kernel_ml, kernelInd_MLEMT++, sizeof(cl_mem), &d_ycenter);
				status = clSetKernelArg(kernel_ml, kernelInd_MLEMT++, sizeof(cl_mem), &d_zcenter);
				status = clSetKernelArg(kernel_ml, kernelInd_MLEMT++, sizeof(cl_mem), &d_V);
			}
			status = clSetKernelArg(kernel_ml, kernelInd_MLEMT++, sizeof(cl_mem), &d_reko_type_mlem);
		}


		status = clFinish(af_queue);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Queue finish failed\n");
			mexEvalString("pause(.0001);");
			return;
		}
		//else if (status == CL_SUCCESS) {
		//	//getErrorString(status);
		//	mexPrintf("Queue finish succeeded\n");
		//	mexEvalString("pause(.0001);");
		//}

		// Loop through each iteration
		for (uint32_t iter = iter0; iter < Niter; iter++) {


			// Compute any of the other algorithms, if applicable
			if (osem_bool) {


				// Loop through the subsets
				for (uint32_t osa_iter = osa_iter0; osa_iter < subsets; osa_iter++) {


					if (compute_norm_matrix == 1u) {
						if (atomic_64bit) {
							Summ[0] = constant(0ULL, im_dim, 1, u64);
						}
						else
							Summ[0] = constant(0.f, im_dim, 1);
						d_Summ = Summ[0].device<cl_mem>();
					}
					else {
						if (no_norm == 0u && tt > 0u) {
							if (atomic_64bit)
								Summ[osa_iter] = constant(0ULL, im_dim, 1, u64);
							else
								Summ[osa_iter] = constant(0.f, im_dim, 1);
							d_Summ = Summ[osa_iter].device<cl_mem>();
						}
						else if (no_norm == 1u) {
							if (atomic_64bit)
								apu_sum = constant(0ULL, 1, 1, u64);
							else
								apu_sum = constant(0.f, 1, 1, f32);
							d_Summ = apu_sum.device<cl_mem>();
						}
						else
							d_Summ = Summ[osa_iter].device<cl_mem>();
					}

					if (use_psf) {
						vec.im_os_blurred = computeConvolution(vec.im_os, g, Nx, Ny, Nz, w_vec, n_rekos2);
						af::sync();
					}


					update_opencl_inputs(vec, vec_opencl, false, im_dim, n_rekos2, n_rekos_mlem, MethodList, atomic_64bit, use_psf);

					cl_uint kernelInd_OSEMSubIter = kernelInd_OSEMTIter;

					size_t erotus = length[osa_iter] % local_size;

					if (erotus > 0)
						erotus = (local_size - erotus);

					const size_t global_size = length[osa_iter] + erotus;

					const uint64_t m_size = static_cast<uint64_t>(length[osa_iter]);

					af::sync();

					// Set kernel arguments
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_mem), &d_norm[osa_iter]);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_mem), &d_scat[osa_iter]);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_mem), d_Summ);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_mem), &d_lor[osa_iter]);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_mem), &d_xyindex[osa_iter]);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_mem), &d_zindex[osa_iter]);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_mem), &d_L[osa_iter]);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_mem), &d_Sino[osa_iter]);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_mem), &d_sc_ra[osa_iter]);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_mem), vec_opencl.d_im_os);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_mem), vec_opencl.d_rhs_os);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_uchar), &no_norm);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(uint64_t), &m_size);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_ulong), &st);
					// Compute the kernel
					status = clEnqueueNDRangeKernel(af_queue, kernel, 1u, NULL, &global_size, &local_size, 0, NULL, NULL);

					if (status != CL_SUCCESS) {
						getErrorString(status);
						mexPrintf("Failed to launch the OS kernel\n");
						mexEvalString("pause(.0001);");
						break;
					}
					//else if (verbose) {
					//	mexPrintf("OS kernel launched successfully\n");
					//	mexEvalString("pause(.0001);");
					//}
					status = clFinish(af_queue);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						mexPrintf("Queue finish failed after kernel\n");
						mexEvalString("pause(.0001);");
						break;
					}
					array* testi;

					// Transfer memory control back to ArrayFire (OS-methods)
					if (compute_norm_matrix == 1u) {
						Summ[0].unlock();
						if (atomic_64bit)
							Summ[0] = Summ[0].as(f32) / TH;
						// Prevent division by zero
						Summ[0](Summ[0] < epps) = epps;
						if (use_psf) {
							Summ[0] = computeConvolution(Summ[0], g, Nx, Ny, Nz, w_vec, 1u);
							af::sync();
						}
						testi = &Summ[0];
						eval(*testi);
					}
					else {
						if (no_norm == 0u) {
							Summ[osa_iter].unlock();
							if (atomic_64bit) {
								Summ[osa_iter] = Summ[osa_iter].as(f32) / TH;
							}
							Summ[osa_iter](Summ[osa_iter] < epps) = epps;
							if (use_psf) {
								Summ[osa_iter] = computeConvolution(Summ[osa_iter], g, Nx, Ny, Nz, w_vec, 1u);
								af::sync();
							}
						}
						else
							apu_sum.unlock();
						// Prevent division by zero
						testi = &Summ[osa_iter];
						eval(*testi);
					}
					if (use_psf)
						vec.im_os_blurred.unlock();
					else
						vec.im_os.unlock();
					vec.rhs_os.unlock();
					if (atomic_64bit)
						vec.rhs_os = vec.rhs_os.as(f32) / TH;
					vec.rhs_os(vec.rhs_os < epps) = epps;
					if (use_psf) {
						vec.rhs_os = computeConvolution(vec.rhs_os, g, Nx, Ny, Nz, w_vec, n_rekos2);
						af::sync();
					}

					computeOSEstimates(vec, w_vec, MethodList, im_dim, testi, epps, iter, osa_iter, subsets, beta, Nx, Ny, Nz, data, length, d_Sino, break_iter, pj3,
						n_rekos2, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, Summ, kernel_mramla, d_L, raw, MethodListOpenCL, koko, atomic_64bit,
						compute_norm_matrix, OpenCLStruct.kernelNLM, d_sc_ra, kernelInd_MRAMLA, E, d_norm, d_scat, use_psf, g, OpenCLStruct);

					vec.im_os(vec.im_os < epps) = epps;

					if (verbose) {
						mexPrintf("Sub-iteration %d complete\n", osa_iter + 1u);
						mexEvalString("pause(.0001);");
					}

					clFinish(af_queue);

					if (break_iter)
						break;

				}
				//vec.im_os = vec.rhs_os;

				computeOSEstimatesIter(vec, w_vec, MethodList, im_dim, epps, iter, osa_iter0, subsets, beta, Nx, Ny, Nz, data, n_rekos2, OpenCLStruct);
				
				if (use_psf && w_vec.deconvolution && osem_bool) {
					computeDeblur(vec, g, Nx, Ny, Nz, w_vec, MethodList, iter, deblur_iterations);
				}

				if (osem_bool && compute_norm_matrix == 0u)
					no_norm = 1u;

				if (verbose) {
					mexPrintf("Iteration %d complete\n", iter + 1u);
					mexEvalString("pause(.0001);");
				}
			}
			// Compute MLEM separately
			if (mlem_bool) {

				// Use previously computed normalization factor if available
				if (compute_norm_matrix == 0u) {
					if (osem_bool && no_norm == 1u && iter == 0u) {
						for (int kk = 0; kk < subsets; kk++)
							Summ_mlem += Summ[kk];
					}
					if (no_norm_mlem == 1u) {
						if (atomic_64bit)
							apu_sum_mlem = constant(0ULL, 1, 1, u64);
						else
							apu_sum_mlem = constant(0.f, 1, 1, f32);
						d_Summ_mlem = apu_sum_mlem.device<cl_mem>();
					}
					else
						d_Summ_mlem = Summ_mlem.device<cl_mem>();
				}
				else {
					if (no_norm_mlem == 1u) {
						if (atomic_64bit)
							apu_sum_mlem = constant(0ULL, 1, 1, u64);
						else
							apu_sum_mlem = constant(0.f, 1, 1, f32);
						d_Summ_mlem = apu_sum_mlem.device<cl_mem>();
					}
					else
						d_Summ_mlem = Summ_mlem.device<cl_mem>();
				}

				cl_uint kernelInd_MLEMSubIter = kernelInd_MLEMT;
				size_t erotus = koko % local_size;

				if (erotus > 0)
					erotus = (local_size - erotus);

				const size_t global_size = koko + erotus;

				const uint64_t m_size = static_cast<uint64_t>(koko);

				if (use_psf) {
					vec.im_mlem_blurred = computeConvolution(vec.im_mlem, g, Nx, Ny, Nz, w_vec, n_rekos_mlem);
					af::sync();
				}

				// Update the OpenCL inputs for this iteration (image estimates)
				update_opencl_inputs(vec, vec_opencl, true, im_dim, n_rekos, n_rekos_mlem, MethodList, atomic_64bit, use_psf);

				clSetKernelArg(kernel_ml, kernelInd_MLEMSubIter++, sizeof(cl_mem), &d_norm_mlem);
				clSetKernelArg(kernel_ml, kernelInd_MLEMSubIter++, sizeof(cl_mem), &d_scat_mlem);
				clSetKernelArg(kernel_ml, kernelInd_MLEMSubIter++, sizeof(cl_mem), d_Summ_mlem);
				clSetKernelArg(kernel_ml, kernelInd_MLEMSubIter++, sizeof(cl_mem), &d_lor_mlem);
				clSetKernelArg(kernel_ml, kernelInd_MLEMSubIter++, sizeof(cl_mem), &d_xyindex_mlem);
				clSetKernelArg(kernel_ml, kernelInd_MLEMSubIter++, sizeof(cl_mem), &d_zindex_mlem);
				clSetKernelArg(kernel_ml, kernelInd_MLEMSubIter++, sizeof(cl_mem), &d_L_mlem);
				clSetKernelArg(kernel_ml, kernelInd_MLEMSubIter++, sizeof(cl_mem), &d_Sino_mlem);
				clSetKernelArg(kernel_ml, kernelInd_MLEMSubIter++, sizeof(cl_mem), &d_sc_ra_mlem);
				clSetKernelArg(kernel_ml, kernelInd_MLEMSubIter++, sizeof(cl_mem), vec_opencl.d_im_mlem);
				clSetKernelArg(kernel_ml, kernelInd_MLEMSubIter++, sizeof(cl_mem), vec_opencl.d_rhs_mlem);
				clSetKernelArg(kernel_ml, kernelInd_MLEMSubIter++, sizeof(cl_uchar), &no_norm_mlem);
				clSetKernelArg(kernel_ml, kernelInd_MLEMSubIter++, sizeof(uint64_t), &m_size);
				clSetKernelArg(kernel_ml, kernelInd_MLEMSubIter++, sizeof(cl_ulong), &st);
				status = clEnqueueNDRangeKernel(af_queue, kernel_ml, 1, NULL, &global_size, &local_size, 0, NULL, NULL);

				if (status != CL_SUCCESS) {
					getErrorString(status);
					mexPrintf("Failed to launch the MLEM kernel\n");
					mexEvalString("pause(.0001);");
					break;
				}
				//else if (verbose) {
				//	mexPrintf("OpenCL kernel executed successfully\n");
				//	mexEvalString("pause(.0001);");
				//}
				status = clFinish(af_queue);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					mexPrintf("Queue finish failed after kernel\n");
					mexEvalString("pause(.0001);");
					break;
				}
				// Transfer memory control back to ArrayFire (ML-methods)
				if (no_norm_mlem == 0u) {
					Summ_mlem.unlock();
					if (atomic_64bit)
						Summ_mlem = Summ_mlem.as(f32) / TH;
					// Prevents division by zero
					Summ_mlem(Summ_mlem < epps) = epps;
					if (use_psf) {
						Summ_mlem = computeConvolution(Summ_mlem, g, Nx, Ny, Nz, w_vec, 1u);
						af::sync();
					}
				}
				else
					apu_sum_mlem.unlock();

				if (use_psf)
					vec.im_mlem_blurred.unlock();
				else
					vec.im_mlem.unlock();
				vec.rhs_mlem.unlock();
				if (atomic_64bit)
					vec.rhs_mlem = vec.rhs_mlem.as(f32) / TH;
				vec.rhs_mlem(vec.rhs_mlem < epps) = epps;

				if (use_psf) {
					vec.rhs_mlem = computeConvolution(vec.rhs_mlem, g, Nx, Ny, Nz, w_vec, n_rekos_mlem);
					af::sync();
				}

				computeMLEstimates(vec, w_vec, MethodList, im_dim, epps, iter, subsets, beta, Nx, Ny, Nz, data, Summ_mlem, break_iter, OpenCLStruct);

				if (no_norm_mlem == 0u)
					no_norm_mlem = 1u;

				if (use_psf && w_vec.deconvolution) {
					computeDeblurMLEM(vec, g, Nx, Ny, Nz, w_vec, MethodList, iter, deblur_iterations);
				}
				if (verbose) {
					mexPrintf("MLEM iteration %d complete\n", iter + 1u);
					mexEvalString("pause(.0001);");
				}
			}
			if (break_iter)
				break;
		}

		// Transfer the device data to host MATLAB cell array
		device_to_host_cell(ArrayList, MethodList, vec, oo, cell, w_vec);

		if (verbose) {
			mexPrintf("Time step %d complete\n", tt + 1u);
			mexEvalString("pause(.0001);");
		}
		w_vec.MBSREM_prepass = false;
	}

	// Transfer memory control of all variables that weren't used
	//unlock_AF_im_vectors(vec, MethodList, true, mlem_bool, osem_bool, 0u);

	// Clear OpenCL buffers
	status = clReleaseMemObject(d_x);
	if (status != CL_SUCCESS)
		getErrorString(status);
	status = clReleaseMemObject(d_y);
	if (status != CL_SUCCESS)
		getErrorString(status);
	status = clReleaseMemObject(d_z);
	if (status != CL_SUCCESS)
		getErrorString(status);
	status = clReleaseMemObject(d_atten);
	if (status != CL_SUCCESS)
		getErrorString(status);
	status = clReleaseMemObject(d_pseudos);
	if (status != CL_SUCCESS)
		getErrorString(status);
	status = clReleaseMemObject(d_xcenter);
	if (status != CL_SUCCESS)
		getErrorString(status);
	status = clReleaseMemObject(d_ycenter);
	if (status != CL_SUCCESS)
		getErrorString(status);
	status = clReleaseMemObject(d_zcenter);
	if (status != CL_SUCCESS)
		getErrorString(status);
	status = clReleaseMemObject(d_V);
	if (status != CL_SUCCESS)
		getErrorString(status);
	if (osem_bool) {
		status = clReleaseMemObject(d_reko_type);
		if (status != CL_SUCCESS)
			getErrorString(status);
		for (uint32_t kk = 0u; kk < subsets; kk++) {
			status = clReleaseMemObject(d_lor[kk]);
			if (status != CL_SUCCESS)
				getErrorString(status);
			status = clReleaseMemObject(d_norm[kk]);
			if (status != CL_SUCCESS)
				getErrorString(status);
			status = clReleaseMemObject(d_scat[kk]);
			if (status != CL_SUCCESS)
				getErrorString(status);
			status = clReleaseMemObject(d_xyindex[kk]);
			if (status != CL_SUCCESS)
				getErrorString(status);
			status = clReleaseMemObject(d_zindex[kk]);
			if (status != CL_SUCCESS)
				getErrorString(status);
			status = clReleaseMemObject(d_L[kk]);
			if (status != CL_SUCCESS)
				getErrorString(status);
			status = clReleaseMemObject(d_Sino[kk]);
			if (status != CL_SUCCESS)
				getErrorString(status);
			//if (randoms_correction)
			status = clReleaseMemObject(d_sc_ra[kk]);
			if (status != CL_SUCCESS)
				getErrorString(status);
		}
	}
	if (mlem_bool) {
		status = clReleaseMemObject(d_reko_type_mlem);
		if (status != CL_SUCCESS)
			getErrorString(status);
		status = clReleaseMemObject(d_lor_mlem);
		if (status != CL_SUCCESS)
			getErrorString(status);
		status = clReleaseMemObject(d_xyindex_mlem);
		if (status != CL_SUCCESS)
			getErrorString(status);
		status = clReleaseMemObject(d_zindex_mlem);
		if (status != CL_SUCCESS)
			getErrorString(status);
		status = clReleaseMemObject(d_L_mlem);
		if (status != CL_SUCCESS)
			getErrorString(status);
		status = clReleaseMemObject(d_norm_mlem);
		if (status != CL_SUCCESS)
			getErrorString(status);
		status = clReleaseMemObject(d_scat_mlem);
		if (status != CL_SUCCESS)
			getErrorString(status);
		status = clReleaseMemObject(d_Sino_mlem);
		if (status != CL_SUCCESS)
			getErrorString(status);
		//if (randoms_correction)
		status = clReleaseMemObject(d_sc_ra_mlem);
		if (status != CL_SUCCESS)
			getErrorString(status);
	}

	// Release program and kernels
	clReleaseProgram(program_os);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Failed to release OS program\n");
	}
	clReleaseProgram(program_ml);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Failed to release ML program\n");
	}
	clReleaseProgram(program_mbsrem);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		mexPrintf("Failed to release prepass program\n");
	}
	if (osem_bool) {
		status = clReleaseKernel(kernel);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to release OS kernel\n");
		}
	}
	if (mlem_bool) {
		status = clReleaseKernel(kernel_ml);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to release MLEM kernel\n");
		}
	}
	if ((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI) ||
		MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0) {
		status = clReleaseKernel(kernel_mramla);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to release prepass kernel\n");
		}
	}
	if (MethodList.NLM) {
		status = clReleaseKernel(OpenCLStruct.kernelNLM);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Failed to release NLM kernel\n");
		}
	}
	clFinish(af_queue);
	af::sync();

	return;
}