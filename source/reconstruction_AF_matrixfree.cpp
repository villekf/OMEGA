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
void reconstruction_AF_matrixfree(const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, const mxArray* Sin,
	const mxArray* sc_ra, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const uint32_t Niter, const mxArray* options, const float dx,
	const float dy, const float dz, const float bx, const float by, const float bz, const float bzb, const float maxxx, const float maxyy, const float zmax,
	const float NSlices, const uint32_t* pituus, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t size_x, const uint32_t TotSinos, 
	mxArray* cell, const mwSize* dimmi, const bool verbose, const uint32_t randoms_correction, const uint32_t attenuation_correction,
	const uint32_t normalization, const float* atten, const size_t size_atten, const float* norm, const size_t size_norm, const uint32_t subsets,
	const float epps, const char* k_path, const uint32_t Nt, const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L, 
	const uint8_t raw, const size_t size_z, const bool osem_bool, const char* fileName, const bool force_build, const float tube_width, 
	const float crystal_size_z, const float* x_center, const float* y_center, const float* z_center, const size_t size_center_x, const size_t size_center_y, 
	const size_t size_of_x, const size_t size_center_z, const uint32_t projector_type, const char* header_directory, const bool precompute, 
	const uint32_t device, const int32_t dec, const uint16_t n_rays, const float cr_pz, const bool use_64bit_atomics, uint32_t n_rekos, 
	const uint32_t n_rekos_mlem, const uint8_t* reko_type) {



	// Number of voxels
	const uint32_t Nxy = Nx * Ny;
	const uint32_t im_dim = Nxy * Nz;

	// Distance between rays in multi-ray Siddon
	const float dc_z = cr_pz / 3.f;

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

	// Obtain the reconstruction methods used
	get_rec_methods(options, MethodList);
	OpenCLRecMethods(MethodList, MethodListOpenCL);

	// Create the OpenCL context and command queue and assign the device
	cl_context af_context = afcl::getContext();
	cl_device_id af_device_id = afcl::getDeviceId();
	cl_command_queue af_queue = afcl::getQueue();

	cl_program program = NULL;

	cl_int status = CL_SUCCESS;

	cl_ulong mem;
	cl_ulong mem_loc;
	status = clGetDeviceInfo(af_device_id, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &mem, NULL);
	if ((static_cast<cl_float>(mem) * mem_portions) > image_bytes && !MethodList.CUSTOM)
		compute_norm_matrix = 0u;

	status = clGetDeviceInfo(af_device_id, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &mem_loc, NULL);

	// Check if there is local enough memory
	if (mem_loc <= static_cast<cl_ulong>(n_rekos) * sizeof(cl_float) * local_size) {
		mexErrMsgTxt("Too many algorithms selected.");
		mexPrintf("Reduce the number of algorithms to %u\n", mem_loc / (sizeof(cl_float) * local_size));
	}

	if (static_cast<uint64_t>(n_rekos) * static_cast<uint64_t>(im_dim) > 4294967295ULL) {
		mexErrMsgTxt("Too many algorithms selected.");
		mexPrintf("Reduce the number of algorithms to %u\n", 4294967295ULL / static_cast<uint64_t>(im_dim));
	}

	std::string new_filename;

	// Build the OpenCL program and save the binaries
	// Force the build regardless of whether the binaries already exist
	if (force_build) {
		status = SaveProgramBinary(verbose, k_path, af_context, af_device_id, fileName, program, atomic_64bit, device, header_directory, force_build);
		if (status != CL_SUCCESS) {
			std::cerr << "Error while saving binary" << std::endl;
			return;
		}
	}
	else {
		FILE* fp = NULL;// = fopen(fileName, "rb");
		new_filename = fileName;
		if (atomic_64bit)
			new_filename += (std::to_string(device) + "_64atom.bin");
		else
			new_filename += (std::to_string(device) + ".bin");
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
		fopen_s(&fp, new_filename.c_str(), "rb");
#else
		fp = fopen(new_filename.c_str(), "rb");
#endif
		if (fp == NULL) {
			// If the binaries do not yet exist, create them
			status = SaveProgramBinary(verbose, k_path, af_context, af_device_id, fileName, program, atomic_64bit, device, header_directory, false);
			if (status != CL_SUCCESS && atomic_64bit) {
				atomic_64bit = false;
				new_filename = fileName;
				new_filename += (std::to_string(device) + ".bin");
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
				fopen_s(&fp, new_filename.c_str(), "rb");
#else
				fp = fopen(new_filename.c_str(), "rb");
#endif
				if (fp == NULL) {
					status = SaveProgramBinary(verbose, k_path, af_context, af_device_id, fileName, program, atomic_64bit, device, header_directory, false);
					if (status != CL_SUCCESS) {
						std::cerr << "Error while saving binary" << std::endl;
						return;
					}
				}
				else {
					status = CreateProgramFromBinary(af_context, af_device_id, fp, program);
					fclose(fp);
					if (status != CL_SUCCESS) {
						mexPrintf("Failed to load OpenCL binaries\n");
						clReleaseProgram(program);
						return;
					}
					else {
						mexPrintf("OpenCL binaries successfully loaded\n");
						mexEvalString("pause(.0001);");
					}
				}
			}
			else if (status != CL_SUCCESS) {
				std::cerr << "Error while saving binary" << std::endl;
				return;
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
			else {
				mexPrintf("OpenCL binaries successfully loaded\n");
				mexEvalString("pause(.0001);");
			}
		}
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
	form_data_variables(vec, beta, w_vec, options, Nx, Ny, Nz, Niter, x00, im_dim, koko, MethodList, data, subsets, osa_iter0);

	// Power factor for ACOSEM
	const float hh = 1.f / w_vec.h_ACOSEM;


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

	status = createKernels(kernel_ml, kernel, kernel_mramla, osem_bool, program, MethodList, w_vec, projector_type, mlem_bool, precompute); 
	if (status != CL_SUCCESS) {
		mexPrintf("Failed to create kernels\n");
		clReleaseProgram(program);
		clReleaseKernel(kernel);
		clReleaseKernel(kernel_ml);
		clReleaseKernel(kernel_mramla);
		return;
	}
	else {
		mexPrintf("OpenCL kernels successfully created\n");
		mexEvalString("pause(.0001);");
	}

	// Create and write buffers
	cl_mem d_x, d_y, d_z, d_pseudos, d_atten, d_xcenter, d_ycenter, d_zcenter, d_norm, d_reko_type;
	cl_mem* d_Summ, *d_Summ_mlem;
	cl_mem d_lor_mlem, d_L_mlem, d_zindex_mlem, d_xyindex_mlem, d_Sino_mlem, d_sc_ra_mlem;

	std::vector<cl_mem> d_lor(subsets);
	std::vector<cl_mem> d_L(subsets);
	std::vector<cl_mem> d_zindex(subsets);
	std::vector<cl_mem> d_xyindex(subsets);
	std::vector<cl_mem> d_Sino(subsets);
	std::vector<cl_mem> d_sc_ra(subsets);

	float *apu = (float*)mxGetData(mxGetCell(Sin, 0));

	status = createAndWriteBuffers(d_x, d_y, d_z, d_lor, d_L, d_zindex, d_xyindex, d_Sino, d_sc_ra, size_x, size_z, TotSinos, size_atten, size_norm, prows, 
		length, x, y, z_det, xy_index, z_index, lor1, L, apu, raw, af_context, subsets, pituus, atten, norm, pseudos, af_queue, d_atten, d_norm, d_pseudos, 
		d_xcenter, d_ycenter, d_zcenter, x_center, y_center, z_center, size_center_x, size_center_y, size_center_z, size_of_x, atomic_64bit, randoms_correction, 
		sc_ra, precompute, d_lor_mlem, d_L_mlem, d_zindex_mlem, d_xyindex_mlem, d_Sino_mlem, d_sc_ra_mlem, d_reko_type, osem_bool, mlem_bool, koko, 
		reko_type, n_rekos);
	if (status != CL_SUCCESS) {
		clReleaseProgram(program);
		clReleaseKernel(kernel);
		clReleaseKernel(kernel_ml);
		clReleaseKernel(kernel_mramla);
		return;
	}

	// Compute the prepass phase for MRAMLA, MBSREM, RBI, COSEM, ACOSEM or ECOSEM if applicable
	if (((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBI || MethodList.RBIMAP) && w_vec.MBSREM_prepass ||
		MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0) && (!MethodList.CUSTOM || osa_iter0 == 0u)) {

		// Set the kernel parameters that do not change
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint32_t), &det_per_ring);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint8_t), &raw);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(cl_mem), &d_pseudos);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint32_t), &prows);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &hh);
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
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(cl_mem), &d_x);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(cl_mem), &d_y);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(cl_mem), &d_z);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint32_t), &size_x);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint32_t), &TotSinos);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint32_t), &attenuation_correction);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint32_t), &normalization);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint32_t), &randoms_correction);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(cl_mem), &d_atten);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(cl_mem), &d_norm);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &epps);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint32_t), &im_dim);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(MethodListOpenCL), &MethodListOpenCL);
		clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint32_t), &Nxy);
		if (projector_type == 2u) {
			clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &tube_width);
			clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &crystal_size_z);
			clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(int32_t), &dec);
			clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(cl_mem), &d_xcenter);
			clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(cl_mem), &d_ycenter);
			clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(cl_mem), &d_zcenter);
		}
		else if (projector_type == 1u && !precompute) {
			clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(float), &dc_z);
			clSetKernelArg(kernel_mramla, kernelInd_MRAMLA++, sizeof(uint16_t), &n_rays);
		}

		uint32_t alku = 0u;

		if ((MethodList.MRAMLA || MethodList.MBSREM) && w_vec.MBSREM_prepass)
			w_vec.Amin = constant(0.f, koko, 1);

		// Run the prepass phase
		MRAMLA_prepass(subsets, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, d_Sino, koko, x00, vec.C_co,
			vec.C_aco, vec.C_osl, alku, kernel_mramla, d_L, raw, MethodListOpenCL, length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E);


		if (verbose) {
			mexPrintf("MRAMLA & COSEM prepass completed\n");
			mexEvalString("pause(.0001);");
		}
	}

	if (MethodList.MRAMLA || MethodList.MBSREM) {
		pj3 = w_vec.D / static_cast<float>(subsets);
	}

	cl_uint kernelInd_OSEM = 0U;
	cl_uint kernelInd_MLEM = 0U;

	if (osem_bool) {
		// Set the kernel parameters that do not change
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(cl_mem), &d_reko_type);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint8_t), &raw);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &hh);
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
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(cl_mem), &d_x);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(cl_mem), &d_y);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(cl_mem), &d_z);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint32_t), &size_x);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint32_t), &TotSinos);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint32_t), &attenuation_correction);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint32_t), &normalization);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint32_t), &randoms_correction);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(cl_mem), &d_atten);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(cl_mem), &d_norm);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &epps);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint32_t), &im_dim);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(cl_mem), &d_pseudos);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint32_t), &prows);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint32_t), &Nxy);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint32_t), &det_per_ring);
		clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint32_t), &n_rekos);
		if (projector_type == 2u) {
			clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &tube_width);
			clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &crystal_size_z);
			clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(int32_t), &dec);
			clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(cl_mem), &d_xcenter);
			clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(cl_mem), &d_ycenter);
			clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(cl_mem), &d_zcenter);
		}
		else if (projector_type == 1u && !precompute) {
			clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(float), &dc_z);
			clSetKernelArg(kernel, kernelInd_OSEM++, sizeof(uint16_t), &n_rays);
		}
	}

	if (mlem_bool) {
		// Set the kernel parameters that do not change
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint8_t), &raw);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &Nx);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &Ny);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &Nz);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &dz);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &dx);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &dy);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &bz);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &bx);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &by);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &bzb);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &maxxx);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &maxyy);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &zmax);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &NSlices);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(cl_mem), &d_x);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(cl_mem), &d_y);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(cl_mem), &d_z);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &size_x);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &TotSinos);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &attenuation_correction);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &normalization);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &randoms_correction);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(cl_mem), &d_atten);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(cl_mem), &d_norm);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &epps);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &im_dim);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(cl_mem), &d_pseudos);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &prows);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &Nxy);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &det_per_ring);
		clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint32_t), &n_rekos_mlem);
		if (projector_type == 2u) {
			clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &tube_width);
			clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &crystal_size_z);
			clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(int32_t), &dec);
			clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(cl_mem), &d_xcenter);
			clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(cl_mem), &d_ycenter);
			clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(cl_mem), &d_zcenter);
		}
		else if (projector_type == 1u && !precompute) {
			clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(float), &dc_z);
			clSetKernelArg(kernel_ml, kernelInd_MLEM++, sizeof(uint16_t), &n_rays);
		}
	}

	if (mlem_bool) {
		if (!atomic_64bit && compute_norm_matrix == 0u && osem_bool)
			Summ_mlem = constant(0.f, im_dim, 1);
		else
			Summ_mlem = constant(0ULL, im_dim, 1, u64);
	}

	// Loop through each time-step
	for (uint32_t tt = t0; tt < Nt; tt++) {

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
						std::cerr << getErrorString(status) << std::endl;
					}
					if (randoms_correction) {
						float* ra_apu = (float*)mxGetData(mxGetCell(sc_ra, tt));
						status = clEnqueueWriteBuffer(af_queue, d_sc_ra[kk], CL_TRUE, 0, sizeof(float) * length[kk], &ra_apu[pituus[kk]], 0, NULL, NULL);
					}
					else
						status = clEnqueueFillBuffer(af_queue, d_sc_ra[kk], &zerof, sizeof(cl_float), 0, sizeof(cl_float), 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						std::cerr << getErrorString(status) << std::endl;
					}
				}
				vec.im_os = constant(0.f, im_dim * n_rekos2, 1);
				for (int kk = 0; kk < n_rekos2; kk++) {
					vec.im_os(seq(kk * im_dim, (kk + 1) * im_dim - 1)) = x00;
				}
			}
			if (mlem_bool) {
				status = clEnqueueWriteBuffer(af_queue, d_Sino_mlem, CL_TRUE, 0, sizeof(float) * koko, apu, 0, NULL, NULL);
				if (randoms_correction) {
					float* ra_apu = (float*)mxGetData(mxGetCell(sc_ra, tt));
					status = clEnqueueWriteBuffer(af_queue, d_sc_ra_mlem, CL_TRUE, 0, sizeof(float) * koko, ra_apu, 0, NULL, NULL);
				}
				else
					status = clEnqueueFillBuffer(af_queue, d_sc_ra_mlem, &zerof, sizeof(cl_float), 0, sizeof(cl_float), 0, NULL, NULL);
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


		status = clFinish(af_queue);
		if (status != CL_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			mexPrintf("Queue finish failed\n");
			mexEvalString("pause(.0001);");
			return;
		}

		// Loop through each iteration
		for (uint32_t iter = iter0; iter < Niter; iter++) {


			// Compute any of the other algorithms, if applicable
			if (osem_bool) {


				// Loop through the subsets
				for (uint32_t osa_iter = osa_iter0; osa_iter < subsets; osa_iter++) {

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

					update_opencl_inputs(vec, vec_opencl, false, im_dim, n_rekos2, n_rekos_mlem, MethodList, atomic_64bit);

					cl_uint kernelInd_OSEMSubIter = kernelInd_OSEM;

					const size_t global_size = length[osa_iter] + (local_size - length[osa_iter] % local_size);

					const uint64_t m_size = static_cast<uint64_t>(length[osa_iter]);


					// Set kernel arguments
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_mem), d_Summ);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_mem), &d_lor[osa_iter]);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_mem), &d_xyindex[osa_iter]);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_mem), &d_zindex[osa_iter]);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_mem), &d_L[osa_iter]);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(float), &w_vec.epsilon_mramla);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_mem), &d_Sino[osa_iter]);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_mem), &d_sc_ra[osa_iter]);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_mem), vec_opencl.d_im_os);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_mem), vec_opencl.d_rhs_os);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_uchar), &no_norm);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(uint64_t), &m_size);
					clSetKernelArg(kernel, kernelInd_OSEMSubIter++, sizeof(cl_float) * static_cast<size_t>(n_rekos) * local_size, NULL);
					// Compute the kernel
					status = clEnqueueNDRangeKernel(af_queue, kernel, 1u, NULL, &global_size, &local_size, 0, NULL, NULL);

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
					array* testi;

					// Transfer memory control back to ArrayFire (OS-methods)
					if (compute_norm_matrix == 1u) {
						Summ[0].unlock();
						if (atomic_64bit)
							Summ[0] = Summ[0].as(f32) / TH;
						// Prevent division by zero
						Summ[0](Summ[0] == 0.f) = epps;
						testi = &Summ[0];
						eval(*testi);
					}
					else {
						if (no_norm == 0u) {
							Summ[osa_iter].unlock();
							if (atomic_64bit) {
								Summ[osa_iter] = Summ[osa_iter].as(f32) / TH;
							}
						}
						else
							apu_sum.unlock();
						// Prevent division by zero
						Summ[osa_iter](Summ[osa_iter] == 0.f) = epps;
						testi = &Summ[osa_iter];
						eval(*testi);
					}
					vec.im_os.unlock();
					vec.rhs_os.unlock();
					uint32_t yy = 0u;
					if (atomic_64bit)
						vec.rhs_os = vec.rhs_os.as(f32) / TH;;
					vec.rhs_os(vec.rhs_os == 0.f) = epps;

					// Compute the (matrix free) algorithms
					// Ordered Subsets Expectation Maximization (OSEM)
					if (MethodList.OSEM || MethodList.ECOSEM) {
						vec.im_os(seq(yy, yy + im_dim - 1u)) = OSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)));
						yy += im_dim;
					}

					// Modfied Row-action Maximum Likelihood (MRAMLA)
					if (MethodList.MRAMLA) {
						vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U, 
							pj3, w_vec.lambda_MBSREM, iter, im_dim, 0.f, af::constant(0.f, 1, 1), *testi, epps);
						yy += im_dim;
					}

					// Row-action Maximum Likelihood (RAMLA)
					if (MethodList.RAMLA) {
						vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
							w_vec.lambda_BSREM, iter);
						yy += im_dim;
					}

					// Relaxed OSEM (ROSEM)
					if (MethodList.ROSEM) {
						vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
							w_vec.lambda_ROSEM, iter);
						yy += im_dim;
					}

					// Rescaled Block Iterative EM (RBI)
					if (MethodList.RBI) {
						vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.D, 
							0.f, af::constant(0.f, 1, 1));
						yy += im_dim;
					}

					// Dynamic RAMLA
					if (MethodList.DRAMA) {
						vec.im_os(seq(yy, yy + im_dim - 1u)) = DRAMA(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
							w_vec.lambda_DRAMA, iter, osa_iter, subsets);
						yy += im_dim;
					}

					// Complete data OSEM
					if (MethodList.COSEM || MethodList.ECOSEM) {
						vec.C_co(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u));
						vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_co, w_vec.D, w_vec.h_ACOSEM, 2u, 
							constant(0.f, 1,1), 0.f);
						yy += im_dim;
					}

					// Enhanced COSEM
					if (MethodList.ECOSEM) {
						vec.im_os(seq(im_dim * (n_rekos2 - 1u), im_dim * n_rekos2 - 1u)) = ECOSEM(vec.im_os(seq(im_dim * (n_rekos2 - 1u), im_dim * n_rekos2 - 1u)), 
							w_vec.D, vec.im_os(seq(0, im_dim - 1u)), vec.im_os(seq(yy - im_dim, yy - 1u)), epps);
						//yy += im_dim;
					}

					// Accelerated COSEM
					if (MethodList.ACOSEM) {
						vec.C_aco(span, osa_iter) = vec.rhs_os(seq(yy, yy + im_dim - 1u));
						vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_aco, w_vec.D, w_vec.h_ACOSEM, 1u, 
							constant(0.f, 1, 1), 0.f);
						array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
						MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, d_Sino,
							koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length,
							atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E);
						vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
						w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
						vec.im_os(seq(yy, yy + im_dim - 1u)) = batchFunc(vec.im_os(seq(yy, yy + im_dim - 1u)), uu / w_vec.ACOSEM_rhs, batchMul);
						yy += im_dim;
					}
					// Median Root Prior
					if (MethodList.MRP) {
						// OSL-OSEM
						if (MethodList.OSLOSEM) {
							array dU = MRP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, 
								w_vec.med_no_norm, im_dim);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_OSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								dU, beta.MRP_OSEM);
							yy += im_dim;
						}
						if (MethodList.BSREM) {
							vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.lambda_BSREM, iter);
							yy += im_dim;
						}
						if (MethodList.MBSREM) {
							array dU = MRP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, 
								w_vec.med_no_norm, im_dim);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U, 
								pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.MRP_MBSREM, dU, *testi, epps);
							yy += im_dim;
						}
						if (MethodList.ROSEMMAP) {
							array dU = MRP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, 
								w_vec.med_no_norm, im_dim);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.lambda_ROSEM, iter);
							yy += im_dim;
						}
						if (MethodList.RBIMAP) {
							array dU = MRP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, 
								w_vec.med_no_norm, im_dim);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.D, beta.MRP_RBI, dU);
							yy += im_dim;
						}
						if (MethodList.OSLCOSEM > 0u) {
							array dU = MRP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, 
								w_vec.med_no_norm, im_dim);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, w_vec.D, w_vec.h_ACOSEM, 
								MethodList.OSLCOSEM, dU, beta.MRP_COSEM);
							if (MethodList.OSLCOSEM == 1u) {
								array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
								MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, 
									d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
									length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E);
								vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
								vec.im_os(seq(yy, yy + im_dim - 1u)) = batchFunc(vec.im_os(seq(yy, yy + im_dim - 1u)), uu / w_vec.ACOSEM_rhs, batchMul);
							}
							yy += im_dim;
						}
					}
					// Quadratic Prior
					if (MethodList.Quad) {
						if (MethodList.OSLOSEM) {
							array dU = Quadratic_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, 
								w_vec.tr_offsets, w_vec.weights_quad, im_dim);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_OSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								dU, beta.Quad_OSEM);
							yy += im_dim;
						}
						if (MethodList.BSREM) {
							vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.lambda_BSREM, iter);
							yy += im_dim;
						}
						if (MethodList.MBSREM) {
							array dU = Quadratic_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, 
								w_vec.tr_offsets, w_vec.weights_quad, im_dim);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U, 
								pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.Quad_MBSREM, dU, *testi, epps);
							yy += im_dim;
						}
						if (MethodList.ROSEMMAP) {
							array dU = Quadratic_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, 
								w_vec.tr_offsets, w_vec.weights_quad, im_dim);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.lambda_ROSEM, iter);
							yy += im_dim;
						}
						if (MethodList.RBIMAP) {
							array dU = Quadratic_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, 
								w_vec.tr_offsets, w_vec.weights_quad, im_dim);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.D, beta.Quad_RBI, dU);
							yy += im_dim;
						}
						if (MethodList.OSLCOSEM > 0u) {
							array dU = Quadratic_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, 
								w_vec.tr_offsets, w_vec.weights_quad, im_dim);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, w_vec.D, w_vec.h_ACOSEM, 
								MethodList.OSLCOSEM, dU, beta.Quad_COSEM);
							if (MethodList.OSLCOSEM == 1u) {
								array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
								MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, 
									d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
									length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E);
								vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
								w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
								vec.im_os(seq(yy, yy + im_dim - 1u)) = batchFunc(vec.im_os(seq(yy, yy + im_dim - 1u)), uu / w_vec.ACOSEM_rhs, batchMul);
							}
							yy += im_dim;
						}
					}
					// L-filter prior
					if (MethodList.L) {
						if (MethodList.OSLOSEM) {
							array dU = L_filter(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, 
								w_vec.a_L, w_vec.med_no_norm, im_dim);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_OSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								dU, beta.L_OSEM);
							yy += im_dim;
						}
						if (MethodList.BSREM) {
							vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.lambda_BSREM, iter);
							yy += im_dim;
						}
						if (MethodList.MBSREM) {
							array dU = L_filter(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, 
								w_vec.a_L, w_vec.med_no_norm, im_dim);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U, 
								pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.L_MBSREM, dU, *testi, epps);
							yy += im_dim;
						}
						if (MethodList.ROSEMMAP) {
							array dU = L_filter(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, 
								w_vec.a_L, w_vec.med_no_norm, im_dim);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.lambda_ROSEM, iter);
							yy += im_dim;
						}
						if (MethodList.RBIMAP) {
							array dU = L_filter(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, 
								w_vec.a_L, w_vec.med_no_norm, im_dim);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.D, beta.L_RBI, dU);
							yy += im_dim;
						}
						if (MethodList.OSLCOSEM > 0u) {
							array dU = L_filter(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, 
								w_vec.a_L, w_vec.med_no_norm, im_dim);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, w_vec.D, w_vec.h_ACOSEM, 
								MethodList.OSLCOSEM, dU, beta.L_COSEM);
							if (MethodList.OSLCOSEM == 1u) {
								array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
								MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, 
									d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
									length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E);
								vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
								w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
								vec.im_os(seq(yy, yy + im_dim - 1u)) = batchFunc(vec.im_os(seq(yy, yy + im_dim - 1u)), uu / w_vec.ACOSEM_rhs, batchMul);
							}
							yy += im_dim;
						}
					}
					// FIR Median Hybrid prior
					if (MethodList.FMH) {
						if (MethodList.OSLOSEM) {
							array dU = FMH(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets,
								w_vec.fmh_weights, w_vec.med_no_norm, w_vec.alku_fmh, im_dim);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_OSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								dU, beta.FMH_OSEM);
							yy += im_dim;
						}
						if (MethodList.BSREM) {
							vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.lambda_BSREM, iter);
							yy += im_dim;
						}
						if (MethodList.MBSREM) {
							array dU = FMH(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets,
								w_vec.fmh_weights, w_vec.med_no_norm, w_vec.alku_fmh, im_dim);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U, 
								pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.FMH_MBSREM, dU, *testi, epps);
							yy += im_dim;
						}
						if (MethodList.ROSEMMAP) {
							array dU = FMH(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets,
								w_vec.fmh_weights, w_vec.med_no_norm, w_vec.alku_fmh, im_dim);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.lambda_ROSEM, iter);
							yy += im_dim;
						}
						if (MethodList.RBIMAP) {
							array dU = FMH(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets,
								w_vec.fmh_weights, w_vec.med_no_norm, w_vec.alku_fmh, im_dim);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.D, beta.FMH_RBI, dU);
							yy += im_dim;
						}
						if (MethodList.OSLCOSEM > 0u) {
							array dU = FMH(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets,
								w_vec.fmh_weights, w_vec.med_no_norm, w_vec.alku_fmh, im_dim);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, w_vec.D, w_vec.h_ACOSEM, 
								MethodList.OSLCOSEM, dU, beta.FMH_COSEM);
							if (MethodList.OSLCOSEM == 1u) {
								array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
								MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, 
									d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
									length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E);
								vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
								w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
								vec.im_os(seq(yy, yy + im_dim - 1u)) = batchFunc(vec.im_os(seq(yy, yy + im_dim - 1u)), uu / w_vec.ACOSEM_rhs, batchMul);
							}
							yy += im_dim;
						}
					}
					// Weighted Mean prior
					if (MethodList.WeightedMean) {
						if (MethodList.OSLOSEM) {
							array dU = Weighted_mean(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
								w_vec.weighted_weights, w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_OSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								dU, beta.Weighted_OSEM);
							yy += im_dim;
						}
						if (MethodList.BSREM) {
							vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.lambda_BSREM, iter);
							yy += im_dim;
						}
						if (MethodList.MBSREM) {
							array dU = Weighted_mean(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
								w_vec.weighted_weights, w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U, 
								pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.Weighted_MBSREM, dU, *testi, epps);
							yy += im_dim;
						}
						if (MethodList.ROSEMMAP) {
							array dU = Weighted_mean(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
								w_vec.weighted_weights, w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.lambda_ROSEM, iter);
							yy += im_dim;
						}
						if (MethodList.RBIMAP) {
							array dU = Weighted_mean(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
								w_vec.weighted_weights, w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.D, beta.Weighted_RBI, dU);
							yy += im_dim;
						}
						if (MethodList.OSLCOSEM > 0u) {
							array dU = Weighted_mean(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
								w_vec.weighted_weights, w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, w_vec.D, w_vec.h_ACOSEM, 
								MethodList.OSLCOSEM, dU, beta.Weighted_COSEM);
							if (MethodList.OSLCOSEM == 1u) {
								array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
								MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, 
									d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw,
									MethodListOpenCL, length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E);
								vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
								w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
								vec.im_os(seq(yy, yy + im_dim - 1u)) = batchFunc(vec.im_os(seq(yy, yy + im_dim - 1u)), uu / w_vec.ACOSEM_rhs, batchMul);
							}
							yy += im_dim;
						}
					}
					// Total Variation prior
					if (MethodList.TV) {
						if (MethodList.OSLOSEM) {
							array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, data.TVtype, w_vec, w_vec.tr_offsets);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_OSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								dU, beta.TV_OSEM);
							yy += im_dim;
						}
						if (MethodList.BSREM) {
							vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.lambda_BSREM, iter);
							yy += im_dim;
						}
						if (MethodList.MBSREM) {
							array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, data.TVtype, w_vec, w_vec.tr_offsets);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U, 
								pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.TV_MBSREM, dU, *testi, epps);
							yy += im_dim;
						}
						if (MethodList.ROSEMMAP) {
							vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.lambda_ROSEM, iter);
							yy += im_dim;
						}
						if (MethodList.RBIMAP) {
							array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, data.TVtype, w_vec, w_vec.tr_offsets);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.D, beta.TV_RBI, dU);
							yy += im_dim;
						}
						if (MethodList.OSLCOSEM > 0u) {
							array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, data.TVtype, w_vec, w_vec.tr_offsets);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, w_vec.D, w_vec.h_ACOSEM, 
								MethodList.OSLCOSEM, dU, beta.TV_COSEM);
							if (MethodList.OSLCOSEM == 1u) {
								array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
								MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, 
									d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
									length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E);
								vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
								w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
								vec.im_os(seq(yy, yy + im_dim - 1u)) = batchFunc(vec.im_os(seq(yy, yy + im_dim - 1u)), uu / w_vec.ACOSEM_rhs, batchMul);
							}
							yy += im_dim;
						}
					}
					// Anisotropic Diffusion smoothing prior
					if (MethodList.AD) {
						if (MethodList.OSLOSEM) {
							if (osa_iter == 0u)
								if (MethodList.OSEM)
									vec.im_os(seq(yy, yy + im_dim - 1u)) = vec.im_os(seq(0, im_dim - 1u));
								else
									vec.im_os(seq(yy, yy + im_dim - 1u)) = OSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)));
							else {
								array dU = AD(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType,
									w_vec.DiffusionType, w_vec.med_no_norm);
								vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_OSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
									dU, beta.AD_OSEM);
							}
							yy += im_dim;
						}
						if (MethodList.BSREM) {
							vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.lambda_BSREM, iter);
							yy += im_dim;
						}
						if (MethodList.MBSREM) {
							if (osa_iter == 0u) {
								vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
									w_vec.U, pj3, w_vec.lambda_MBSREM, iter, im_dim, 0.f, constant(0.f, 1, 1), *testi, epps);
							}
							else {
								array dU = AD(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, 
									w_vec.DiffusionType, w_vec.med_no_norm);
								vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
									w_vec.U, pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.AD_MBSREM, dU, *testi, epps);
							}
							yy += im_dim;
						}
						if (MethodList.ROSEMMAP) {
							vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.lambda_ROSEM, iter);
							yy += im_dim;
						}
						if (MethodList.RBIMAP) {
							if (osa_iter == 0u) {
								vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
									w_vec.D, 0.f, constant(0.f, 1, 1));
							}
							else {
								array dU = AD(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, 
									w_vec.DiffusionType, w_vec.med_no_norm);
								vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
									w_vec.D, beta.AD_RBI, dU);
							}
							yy += im_dim;
						}
						if (MethodList.OSLCOSEM > 0u) {
							if (osa_iter == 0u) {
								vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, w_vec.D, w_vec.h_ACOSEM, 
									MethodList.OSLCOSEM, constant(0.f, 1, 1), 0.f);
								if (MethodList.OSLCOSEM == 1u) {
									array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
									MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ,
										d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
										length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E);
									vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
									w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
									vec.im_os(seq(yy, yy + im_dim - 1u)) = batchFunc(vec.im_os(seq(yy, yy + im_dim - 1u)), uu / w_vec.ACOSEM_rhs, batchMul);
								}
							}
							else {
								array dU = AD(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, 
									w_vec.DiffusionType, w_vec.med_no_norm);
								vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, w_vec.D, w_vec.h_ACOSEM, 
									MethodList.OSLCOSEM, dU, beta.AD_COSEM);
								if (MethodList.OSLCOSEM == 1u) {
									array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
									MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ,
										d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
										length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E);
									vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
									w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
									vec.im_os(seq(yy, yy + im_dim - 1u)) = batchFunc(vec.im_os(seq(yy, yy + im_dim - 1u)), uu / w_vec.ACOSEM_rhs, batchMul);
								}
							}
							yy += im_dim;
						}
					}
					// Asymmetric Parallel Level Sets prior
					if (MethodList.APLS) {
						if (MethodList.OSLOSEM) {
							array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, 4, w_vec, w_vec.tr_offsets);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_OSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								dU, beta.APLS_OSEM);
							yy += im_dim;
						}
						if (MethodList.BSREM) {
							vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.lambda_BSREM, iter);
							yy += im_dim;
						}
						if (MethodList.MBSREM) {
							array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, 4, w_vec, w_vec.tr_offsets);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U, 
								pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.APLS_MBSREM, dU, *testi, epps);
							yy += im_dim;
						}
						if (MethodList.ROSEMMAP) {
							array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, 4, w_vec, w_vec.tr_offsets);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.lambda_ROSEM, iter);
							yy += im_dim;
						}
						if (MethodList.RBIMAP) {
							array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, 4, w_vec, w_vec.tr_offsets);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.D, beta.APLS_RBI, dU);
							yy += im_dim;
						}
						if (MethodList.OSLCOSEM > 0u) {
							array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, 4, w_vec, w_vec.tr_offsets);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, w_vec.D, w_vec.h_ACOSEM, 
								MethodList.OSLCOSEM, dU, beta.APLS_COSEM);
							if (MethodList.OSLCOSEM == 1u) {
								array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
								MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, d_Sino, 
									koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL, length, atomic_64bit, 
									compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E);
								vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
								w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
								vec.im_os(seq(yy, yy + im_dim - 1u)) = batchFunc(vec.im_os(seq(yy, yy + im_dim - 1u)), uu / w_vec.ACOSEM_rhs, batchMul);
							}
							yy += im_dim;
						}
					}
					// Total Generalized Variation prior
					if (MethodList.TGV) {
						if (MethodList.OSLOSEM) {
							array dU = TGV(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_OSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								dU, beta.TGV_OSEM);
							yy += im_dim;
						}
						if (MethodList.BSREM) {
							vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.lambda_BSREM, iter);
							yy += im_dim;
						}
						if (MethodList.MBSREM) {
							array dU = TGV(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U, 
								pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.TGV_MBSREM, dU, *testi, epps);
							yy += im_dim;
						}
						if (MethodList.ROSEMMAP) {
							array dU = TGV(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.lambda_ROSEM, iter);
							yy += im_dim;
						}
						if (MethodList.RBIMAP) {
							array dU = TGV(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.D, beta.TGV_RBI, dU);
							yy += im_dim;
						}
						if (MethodList.OSLCOSEM > 0u) {
							array dU = TGV(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, w_vec.D, w_vec.h_ACOSEM, 
								MethodList.OSLCOSEM, dU, beta.TGV_COSEM);
							if (MethodList.OSLCOSEM == 1u) {
								array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
								MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ, 
									d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
									length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E);
								vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
								w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
								vec.im_os(seq(yy, yy + im_dim - 1u)) = batchFunc(vec.im_os(seq(yy, yy + im_dim - 1u)), uu / w_vec.ACOSEM_rhs, batchMul);
							}
							yy += im_dim;
						}
					}
					// Non-local means prior
					if (MethodList.NLM) {
						if (MethodList.OSLOSEM) {
							array dU = NLM(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, w_vec.Nlx, w_vec.Nly, w_vec.Nlz, w_vec.h2, 
								epps, Nx, Ny, Nz, w_vec.NLM_anatomical, w_vec.NLM_gauss, w_vec.NLTV, w_vec.NLM_MRP, w_vec.NLM_ref);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_OSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
								dU, beta.NLM_OSEM);
							yy += im_dim;
						}
						if (MethodList.BSREM) {
							vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)),
								w_vec.lambda_BSREM, iter);
							yy += im_dim;
						}
						if (MethodList.MBSREM) {
							array dU = NLM(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, w_vec.Nlx, w_vec.Nly, w_vec.Nlz, w_vec.h2,
								epps, Nx, Ny, Nz, w_vec.NLM_anatomical, w_vec.NLM_gauss, w_vec.NLTV, w_vec.NLM_MRP, w_vec.NLM_ref);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U,
								pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.NLM_MBSREM, dU, *testi, epps);
							yy += im_dim;
						}
						if (MethodList.ROSEMMAP) {
							array dU = NLM(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, w_vec.Nlx, w_vec.Nly, w_vec.Nlz, w_vec.h2,
								epps, Nx, Ny, Nz, w_vec.NLM_anatomical, w_vec.NLM_gauss, w_vec.NLTV, w_vec.NLM_MRP, w_vec.NLM_ref);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
								w_vec.lambda_ROSEM, iter);
							yy += im_dim;
						}
						if (MethodList.RBIMAP) {
							array dU = NLM(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, w_vec.Nlx, w_vec.Nly, w_vec.Nlz, w_vec.h2,
								epps, Nx, Ny, Nz, w_vec.NLM_anatomical, w_vec.NLM_gauss, w_vec.NLTV, w_vec.NLM_MRP, w_vec.NLM_ref);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)),
								w_vec.D, beta.NLM_RBI, dU);
							yy += im_dim;
						}
						if (MethodList.OSLCOSEM > 0u) {
							array dU = NLM(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, w_vec.Nlx, w_vec.Nly, w_vec.Nlz, w_vec.h2,
								epps, Nx, Ny, Nz, w_vec.NLM_anatomical, w_vec.NLM_gauss, w_vec.NLTV, w_vec.NLM_MRP, w_vec.NLM_ref);
							vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, w_vec.D, w_vec.h_ACOSEM,
								MethodList.OSLCOSEM, dU, beta.NLM_COSEM);
							if (MethodList.OSLCOSEM == 1u) {
								array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
								MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ,
									d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
									length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E);
								vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
								w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
								vec.im_os(seq(yy, yy + im_dim - 1u)) = batchFunc(vec.im_os(seq(yy, yy + im_dim - 1u)), uu / w_vec.ACOSEM_rhs, batchMul);
							}
							yy += im_dim;
						}
					}
					// Custom prior
					if (MethodList.CUSTOM) {
						if (MethodList.OSLOSEM) {
							vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_OSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.dU_OSEM, beta.custom_OSEM);
							yy += im_dim;
						}
						if (MethodList.BSREM) {
							vec.im_os(seq(yy, yy + im_dim - 1u)) = BSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.lambda_BSREM, iter);
							yy += im_dim;
						}
						if (MethodList.MBSREM) {
							vec.im_os(seq(yy, yy + im_dim - 1u)) = MBSREM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.rhs_os(seq(yy, yy + im_dim - 1u)), w_vec.U, 
								pj3, w_vec.lambda_MBSREM, iter, im_dim, beta.custom_MBSREM, w_vec.dU_MBSREM, *testi, epps);
							yy += im_dim;
						}
						if (MethodList.ROSEMMAP) {
							vec.im_os(seq(yy, yy + im_dim - 1u)) = ROSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.lambda_ROSEM, iter);
							yy += im_dim;
						}
						if (MethodList.RBIMAP) {
							vec.im_os(seq(yy, yy + im_dim - 1u)) = RBI(vec.im_os(seq(yy, yy + im_dim - 1u)), *testi, vec.rhs_os(seq(yy, yy + im_dim - 1u)), 
								w_vec.D, beta.custom_RBI, w_vec.dU_RBI);
							yy += im_dim;
						}
						if (MethodList.OSLCOSEM > 0u) {
							vec.im_os(seq(yy, yy + im_dim - 1u)) = OSL_COSEM(vec.im_os(seq(yy, yy + im_dim - 1u)), vec.C_osl, w_vec.D, w_vec.h_ACOSEM, 
								MethodList.OSLCOSEM, w_vec.dU_COSEM, beta.custom_COSEM);
							if (MethodList.OSLCOSEM == 1u) {
								array apu = vec.im_os(seq(yy, yy + im_dim - 1u));
								MRAMLA_prepass(osa_iter + 1u, im_dim, pituus, d_lor, d_zindex, d_xyindex, program, af_queue, af_context, w_vec, Summ,
									d_Sino, koko, apu, vec.C_co, vec.C_aco, vec.C_osl, osa_iter + 1u, kernel_mramla, d_L, raw, MethodListOpenCL,
									length, atomic_64bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E);
								vec.im_os(seq(yy, yy + im_dim - 1u)) = apu;
								w_vec.ACOSEM_rhs(w_vec.ACOSEM_rhs <= 0.f) = epps;
								vec.im_os(seq(yy, yy + im_dim - 1u)) = batchFunc(vec.im_os(seq(yy, yy + im_dim - 1u)), uu / w_vec.ACOSEM_rhs, batchMul);
							}
							yy += im_dim;
						}
						break_iter = true;
					}


					vec.im_os(vec.im_os < 0.f) = epps;

					if (verbose) {
						mexPrintf("Sub-iteration %d complete\n", osa_iter + 1u);
						mexEvalString("pause(.0001);");
					}

					clFinish(af_queue);

					if (break_iter)
						break;

				}


				uint32_t yy = 0u;
				// Compute BSREM and ROSEMMAP updates if applicable
				// Otherwise simply save the current iterate
				if (MethodList.OSEM || MethodList.ECOSEM) {
					if (MethodList.OSEM)
						vec.OSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
					yy += im_dim;
				}

				if (MethodList.MRAMLA) {
					vec.MRAMLA(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
					yy += im_dim;
				}

				if (MethodList.RAMLA) {
					vec.RAMLA(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
					yy += im_dim;
				}

				if (MethodList.ROSEM) {
					vec.ROSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
					yy += im_dim;
				}

				if (MethodList.RBI) {
					vec.RBI(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
					yy += im_dim;
				}

				if (MethodList.DRAMA) {
					vec.DRAMA(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
					yy += im_dim;
				}

				if (MethodList.COSEM || MethodList.ECOSEM) {
					if (MethodList.COSEM)
						vec.COSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
					yy += im_dim;
				}

				if (MethodList.ECOSEM) {
					vec.ECOSEM(span, iter + 1u) = vec.im_os(seq(im_dim * (n_rekos2 - 1u), im_dim * n_rekos2 - 1u)).copy();
					//yy += im_dim;
				}

				if (MethodList.ACOSEM) {
					vec.ACOSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
					yy += im_dim;
				}

				if (MethodList.MRP) {
					if (MethodList.OSLOSEM) {
						vec.MRP_OSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.BSREM) {
						array dU = MRP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, 
							w_vec.med_no_norm, im_dim);
						vec.MRP_BSREM(span, iter + 1u) = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM, iter, beta.MRP_BSREM, dU, epps);
						yy += im_dim;
					}
					if (MethodList.MBSREM) {
						vec.MRP_MBSREM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.ROSEMMAP) {
						array dU = MRP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, 
							w_vec.med_no_norm, im_dim);
						vec.MRP_ROSEM(span, iter + 1u) = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM, iter, beta.MRP_ROSEM, dU, epps);
						yy += im_dim;
					}
					if (MethodList.RBIMAP) {
						vec.MRP_RBI(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.OSLCOSEM > 0) {
						vec.MRP_COSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
				}

				if (MethodList.Quad) {
					if (MethodList.OSLOSEM) {
						vec.Quad_OSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.BSREM) {
						array dU = Quadratic_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, 
							w_vec.tr_offsets, w_vec.weights_quad, im_dim);
						vec.Quad_BSREM(span, iter + 1u) = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM, iter, beta.Quad_BSREM, dU, epps);
						yy += im_dim;
					}
					if (MethodList.MBSREM) {
						vec.Quad_MBSREM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.ROSEMMAP) {
						array dU = Quadratic_prior(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, 
							w_vec.tr_offsets, w_vec.weights_quad, im_dim);
						vec.Quad_ROSEM(span, iter + 1u) = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM, iter, beta.Quad_ROSEM, dU, epps);
						yy += im_dim;
					}
					if (MethodList.RBIMAP) {
						vec.Quad_RBI(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.OSLCOSEM > 0) {
						vec.Quad_COSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
				}
				if (MethodList.L) {
					if (MethodList.OSLOSEM) {
						vec.L_OSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.BSREM) {
						array dU = L_filter(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L, 
							w_vec.med_no_norm, im_dim);
						vec.L_BSREM(span, iter + 1u) = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM, iter, beta.L_BSREM, dU, epps);
						yy += im_dim;
					}
					if (MethodList.MBSREM) {
						vec.L_MBSREM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.ROSEMMAP) {
						array dU = L_filter(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L, 
							w_vec.med_no_norm, im_dim);
						vec.L_ROSEM(span, iter + 1u) = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM, iter, beta.L_ROSEM, dU, epps);
						yy += im_dim;
					}
					if (MethodList.RBIMAP) {
						vec.L_RBI(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.OSLCOSEM > 0) {
						vec.L_COSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
				}
				if (MethodList.FMH) {
					if (MethodList.OSLOSEM) {
						vec.FMH_OSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.BSREM) {
						array dU = FMH(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets, 
							w_vec.fmh_weights, w_vec.med_no_norm, w_vec.alku_fmh, im_dim);
						vec.FMH_BSREM(span, iter + 1u) = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM, iter, beta.FMH_BSREM, dU, epps);
						yy += im_dim;
					}
					if (MethodList.MBSREM) {
						vec.FMH_MBSREM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.ROSEMMAP) {
						array dU = FMH(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets, 
							w_vec.fmh_weights, w_vec.med_no_norm, w_vec.alku_fmh, im_dim);
						vec.FMH_ROSEM(span, iter + 1u) = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM, iter, beta.FMH_ROSEM, dU, epps);
						yy += im_dim;
					}
					if (MethodList.RBIMAP) {
						vec.FMH_RBI(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.OSLCOSEM > 0) {
						vec.FMH_COSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
				}
				if (MethodList.WeightedMean) {
					if (MethodList.OSLOSEM) {
						vec.Weighted_OSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.BSREM) {
						array dU = Weighted_mean(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
							w_vec.weighted_weights, w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
						vec.Weighted_BSREM(span, iter + 1u) = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM, iter, beta.Weighted_BSREM, dU, 
							epps);
						yy += im_dim;
					}
					if (MethodList.MBSREM) {
						vec.Weighted_MBSREM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.ROSEMMAP) {
						array dU = Weighted_mean(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
							w_vec.weighted_weights, w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
						vec.Weighted_ROSEM(span, iter + 1u) = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM, iter, beta.Weighted_ROSEM, dU, 
							epps);
						yy += im_dim;
					}
					if (MethodList.RBIMAP) {
						vec.Weighted_RBI(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.OSLCOSEM > 0) {
						vec.Weighted_COSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
				}
				if (MethodList.TV) {
					if (MethodList.OSLOSEM) {
						vec.TV_OSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.BSREM) {
						array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, data.TVtype, w_vec, w_vec.tr_offsets);
						vec.TV_BSREM(span, iter + 1u) = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM, iter, beta.TV_BSREM, dU, epps);
						yy += im_dim;
					}
					if (MethodList.MBSREM) {
						vec.TV_MBSREM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.ROSEMMAP) {
						array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, data.TVtype, w_vec, w_vec.tr_offsets);
						vec.TV_ROSEM(span, iter + 1u) = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM, iter, beta.TV_ROSEM, dU, epps);
						yy += im_dim;
					}
					if (MethodList.RBIMAP) {
						vec.TV_RBI(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.OSLCOSEM > 0) {
						vec.TV_COSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
				}
				if (MethodList.AD) {
					if (MethodList.OSLOSEM) {
						vec.AD_OSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.BSREM) {
						array dU = AD(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, 
							w_vec.DiffusionType, w_vec.med_no_norm);
						vec.AD_BSREM(span, iter + 1u) = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM, iter, beta.AD_BSREM, dU, epps);
						yy += im_dim;
					}
					if (MethodList.MBSREM) {
						vec.AD_MBSREM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.ROSEMMAP) {
						array dU = AD(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, 
							w_vec.DiffusionType, w_vec.med_no_norm);
						vec.AD_ROSEM(span, iter + 1u) = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM, iter, beta.AD_ROSEM, dU, epps);
						yy += im_dim;
					}
					if (MethodList.RBIMAP) {
						vec.AD_RBI(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.OSLCOSEM > 0) {
						vec.AD_COSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
				}
				if (MethodList.APLS) {
					if (MethodList.OSLOSEM) {
						vec.APLS_OSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.BSREM) {
						array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, 4, w_vec, w_vec.tr_offsets);
						vec.APLS_BSREM(span, iter + 1u) = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM, iter, beta.APLS_BSREM, dU, epps);
						yy += im_dim;
					}
					if (MethodList.MBSREM) {
						vec.APLS_MBSREM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.ROSEMMAP) {
						array dU = TVprior(Nx, Ny, Nz, data, vec.im_os(seq(yy, yy + im_dim - 1u)), epps, 4, w_vec, w_vec.tr_offsets);
						vec.APLS_ROSEM(span, iter + 1u) = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM, iter, beta.APLS_ROSEM, dU, epps);
						yy += im_dim;
					}
					if (MethodList.RBIMAP) {
						vec.APLS_RBI(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.OSLCOSEM > 0) {
						vec.APLS_COSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
				}
				if (MethodList.TGV) {
					if (MethodList.OSLOSEM) {
						vec.TGV_OSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.BSREM) {
						array dU = TGV(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
						vec.TGV_BSREM(span, iter + 1u) = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM, iter, beta.TGV_BSREM, dU, epps);
						yy += im_dim;
					}
					if (MethodList.MBSREM) {
						vec.TGV_MBSREM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.ROSEMMAP) {
						array dU = TGV(vec.im_os(seq(yy, yy + im_dim - 1u)), Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
						vec.TGV_ROSEM(span, iter + 1u) = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM, iter, beta.TGV_ROSEM, dU, epps);
						yy += im_dim;
					}
					if (MethodList.RBIMAP) {
						vec.TGV_RBI(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.OSLCOSEM > 0) {
						vec.TGV_COSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
				}
				if (MethodList.NLM) {
					if (MethodList.OSLOSEM) {
						vec.NLM_OSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.BSREM) {
						array dU = NLM(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, w_vec.Nlx, w_vec.Nly, w_vec.Nlz, w_vec.h2,
							epps, Nx, Ny, Nz, w_vec.NLM_anatomical, w_vec.NLM_gauss, w_vec.NLTV, w_vec.NLM_MRP, w_vec.NLM_ref);
						vec.NLM_BSREM(span, iter + 1u) = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM, iter, beta.NLM_BSREM, dU, epps);
						yy += im_dim;
					}
					if (MethodList.MBSREM) {
						vec.NLM_MBSREM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.ROSEMMAP) {
						array dU = NLM(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, w_vec.Nlx, w_vec.Nly, w_vec.Nlz, w_vec.h2,
							epps, Nx, Ny, Nz, w_vec.NLM_anatomical, w_vec.NLM_gauss, w_vec.NLTV, w_vec.NLM_MRP, w_vec.NLM_ref);
						vec.NLM_ROSEM(span, iter + 1u) = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM, iter, beta.NLM_ROSEM, dU, epps);
						yy += im_dim;
					}
					if (MethodList.RBIMAP) {
						vec.NLM_RBI(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.OSLCOSEM > 0) {
						vec.NLM_COSEM(span, iter + 1u) = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
				}
				if (MethodList.CUSTOM) {
					if (MethodList.OSLOSEM) {
						vec.custom_OSEM = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.BSREM) {
						if (osa_iter0 == subsets)
							vec.custom_BSREM = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_BSREM, iter, beta.custom_BSREM, w_vec.dU_BSREM, epps);
						else
							vec.custom_BSREM = vec.im_os(seq(yy, yy + im_dim - 1u));
						yy += im_dim;
					}
					if (MethodList.MBSREM) {
						vec.custom_MBSREM = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.ROSEMMAP) {
						if (osa_iter0 == subsets)
							vec.custom_ROSEM = BSREM_MAP(vec.im_os(seq(yy, yy + im_dim - 1u)), w_vec.lambda_ROSEM, iter, beta.custom_ROSEM, w_vec.dU_ROSEM, epps);
						else
							vec.custom_BSREM = vec.im_os(seq(yy, yy + im_dim - 1u));
						yy += im_dim;
					}
					if (MethodList.RBIMAP) {
						vec.custom_RBI = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
					if (MethodList.OSLCOSEM > 0) {
						vec.custom_COSEM = vec.im_os(seq(yy, yy + im_dim - 1u)).copy();
						yy += im_dim;
					}
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

				cl_uint kernelInd_MLEMSubIter = kernelInd_MLEM;

				const size_t global_size = koko + (local_size - koko % local_size);

				const uint64_t m_size = static_cast<uint64_t>(koko);

				// Update the OpenCL inputs for this iteration (image estimates)
				update_opencl_inputs(vec, vec_opencl, true, im_dim, n_rekos, n_rekos_mlem, MethodList, atomic_64bit);

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
				clSetKernelArg(kernel_ml, kernelInd_MLEMSubIter++, sizeof(cl_float) * static_cast<size_t>(n_rekos_mlem) * local_size, NULL);
				status = clEnqueueNDRangeKernel(af_queue, kernel_ml, 1, NULL, &global_size, &local_size, 0, NULL, NULL);
				clFinish(af_queue);

				if (status != CL_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
					mexPrintf("Failed to launch the MLEM kernel\n");
					break;
				}
				else if (verbose) {
					mexPrintf("OpenCL kernel executed successfully\n");
					mexEvalString("pause(.0001);");
				}
				clFinish(af_queue);
				// Transfer memory control back to ArrayFire (ML-methods)
				if (no_norm_mlem == 0u) {
					Summ_mlem.unlock();
					if (atomic_64bit)
						Summ_mlem = Summ_mlem.as(f32) / TH;
				}
				else
					apu_sum_mlem.unlock();

				// Prevents division by zero
				Summ_mlem(Summ_mlem == 0.f) = epps;

				vec.im_mlem.unlock();
				vec.rhs_mlem.unlock();
				uint32_t ee = 0u;
				if (atomic_64bit)
					vec.rhs_mlem = vec.rhs_mlem.as(f32) / TH;
				vec.rhs_mlem(vec.rhs_mlem == 0.f) = epps;

				// Compute the new estimates
				// MLEM
				if (MethodList.MLEM) {
					vec.im_mlem(seq(ee, ee + im_dim - 1u)) = MLEM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), Summ_mlem, vec.rhs_mlem(seq(ee, ee + im_dim - 1u)));
					ee += im_dim;
				}
				// OSL-MLEM with Median Root Prior
				if (MethodList.OSLMLEM && MethodList.MRP) {
					array dU = MRP(vec.im_mlem(seq(ee, ee + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.med_no_norm, 
						im_dim);
					vec.im_mlem(seq(ee, ee + im_dim - 1u)) = OSL_MLEM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), Summ_mlem, vec.rhs_mlem(seq(ee, ee + im_dim - 1u)), 
						dU, beta.MRP_MLEM);
					ee += im_dim;
				}
				// OSL-MLEM with Quadratic prior
				if (MethodList.OSLMLEM && MethodList.Quad) {
					array dU = Quadratic_prior(vec.im_mlem(seq(ee, ee + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, w_vec.inffi, w_vec.tr_offsets,
						w_vec.weights_quad, im_dim);
					vec.im_mlem(seq(ee, ee + im_dim - 1u)) = OSL_MLEM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), Summ_mlem, vec.rhs_mlem(seq(ee, ee + im_dim - 1u)), 
						dU, beta.Quad_MLEM);
					ee += im_dim;
				}
				// OSL-MLEM with L-filter prior
				if (MethodList.OSLMLEM && MethodList.L) {
					array dU = L_filter(vec.im_mlem(seq(ee, ee + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets, w_vec.a_L, 
						w_vec.med_no_norm, im_dim);
					vec.im_mlem(seq(ee, ee + im_dim - 1u)) = OSL_MLEM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), Summ_mlem, vec.rhs_mlem(seq(ee, ee + im_dim - 1u)), 
						dU, beta.L_MLEM);
					ee += im_dim;
				}
				// OSL-MLEM with FIR Median Hybrid prior
				if (MethodList.OSLMLEM && MethodList.FMH) {
					array dU = FMH(vec.im_mlem(seq(ee, ee + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.inffi, w_vec.tr_offsets, 
						w_vec.fmh_weights, w_vec.med_no_norm, w_vec.alku_fmh, im_dim);
					vec.im_mlem(seq(ee, ee + im_dim - 1u)) = OSL_MLEM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), Summ_mlem, vec.rhs_mlem(seq(ee, ee + im_dim - 1u)), 
						dU, beta.FMH_MLEM);
					ee += im_dim;
				}
				// OSL-MLEM with Weighted Mean prior
				if (MethodList.OSLMLEM && MethodList.WeightedMean) {
					array dU = Weighted_mean(vec.im_mlem(seq(ee, ee + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, Nx, Ny, Nz, epps, w_vec.tr_offsets,
						w_vec.weighted_weights, w_vec.med_no_norm, im_dim, w_vec.mean_type, w_vec.w_sum);
					vec.im_mlem(seq(ee, ee + im_dim - 1u)) = OSL_MLEM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), Summ_mlem, vec.rhs_mlem(seq(ee, ee + im_dim - 1u)), 
						dU, beta.Weighted_MLEM);
					ee += im_dim;
				}
				// OSL-MLEM with Total Variation prior
				if (MethodList.OSLMLEM && MethodList.TV) {
					array dU = TVprior(Nx, Ny, Nz, data, vec.im_mlem(seq(ee, ee + im_dim - 1u)), epps, data.TVtype, w_vec, w_vec.tr_offsets);
					vec.im_mlem(seq(ee, ee + im_dim - 1u)) = OSL_MLEM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), Summ_mlem, vec.rhs_mlem(seq(ee, ee + im_dim - 1u)), 
						dU, beta.TV_MLEM);
				}
				// OSL-MLEM with Anisotropic Diffusion smoothing prior
				if (MethodList.OSLMLEM && MethodList.AD) {
					array dU = AD(vec.im_mlem(seq(ee, ee + im_dim - 1u)), Nx, Ny, Nz, epps, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType, 
						w_vec.DiffusionType, w_vec.med_no_norm);
					vec.im_mlem(seq(ee, ee + im_dim - 1u)) = OSL_MLEM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), Summ_mlem, vec.rhs_mlem(seq(ee, ee + im_dim - 1u)), 
						dU, beta.AD_MLEM);
					ee += im_dim;
				}
				// OSL-MLEM with Asymmetric Parallel Level Sets prior
				if (MethodList.OSLMLEM && MethodList.APLS) {
					array dU = TVprior(Nx, Ny, Nz, data, vec.im_mlem(seq(ee, ee + im_dim - 1u)), epps, 4, w_vec, w_vec.tr_offsets);
					vec.im_mlem(seq(ee, ee + im_dim - 1u)) = OSL_MLEM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), Summ_mlem, vec.rhs_mlem(seq(ee, ee + im_dim - 1u)), 
						dU, beta.APLS_MLEM);
					ee += im_dim;
				}
				// OSL-MLEM with Total Generalized Variation prior
				if (MethodList.OSLMLEM && MethodList.TGV) {
					array dU = TGV(vec.im_mlem(seq(ee, ee + im_dim - 1u)), Nx, Ny, Nz, data.NiterTGV, data.TGVAlpha, data.TGVBeta);
					vec.im_mlem(seq(ee, ee + im_dim - 1u)) = OSL_MLEM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), Summ_mlem, vec.rhs_mlem(seq(ee, ee + im_dim - 1u)), 
						dU, beta.TGV_MLEM);
					ee += im_dim;
				}
				// OSL-MLEM with Non-Local Means prior
				if (MethodList.OSLMLEM && MethodList.NLM) {
					array dU = NLM(vec.im_os(seq(ee, ee + im_dim - 1u)), w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, w_vec.Nlx, w_vec.Nly, w_vec.Nlz, w_vec.h2,
						epps, Nx, Ny, Nz, w_vec.NLM_anatomical, w_vec.NLM_gauss, w_vec.NLTV, w_vec.NLM_MRP, w_vec.NLM_ref);
					vec.im_mlem(seq(ee, ee + im_dim - 1u)) = OSL_MLEM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), Summ_mlem, vec.rhs_mlem(seq(ee, ee + im_dim - 1u)),
						dU, beta.NLM_MLEM);
					ee += im_dim;
				}
				// OSL-MLEM with custom prior
				if (MethodList.OSLMLEM && MethodList.CUSTOM) {
					vec.im_mlem(seq(ee, ee + im_dim - 1u)) = OSL_MLEM(vec.im_mlem(seq(ee, ee + im_dim - 1u)), Summ_mlem, vec.rhs_mlem(seq(ee, ee + im_dim - 1u)), 
						w_vec.dU_COSEM, beta.custom_MLEM);
					ee += im_dim;
				}

				ee = 0u;
				 //Save the current iteration
				if (MethodList.MLEM) {
					vec.MLEM(span, iter + 1u) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
					ee += im_dim;
				}
				if (MethodList.OSLMLEM) {
					if (MethodList.MRP) {
						vec.MRP_MLEM(span, iter + 1u) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
						ee += im_dim;
					}
					if (MethodList.Quad) {
						vec.Quad_MLEM(span, iter + 1u) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
						ee += im_dim;
					}
					if (MethodList.L) {
						vec.L_MLEM(span, iter + 1u) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
						ee += im_dim;
					}
					if (MethodList.FMH) {
						vec.FMH_MLEM(span, iter + 1u) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
						ee += im_dim;
					}
					if (MethodList.WeightedMean) {
						vec.Weighted_MLEM(span, iter + 1u) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
						ee += im_dim;
					}
					if (MethodList.TV) {
						vec.TV_MLEM(span, iter + 1u) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
						ee += im_dim;
					}
					if (MethodList.AD) {
						vec.AD_MLEM(span, iter + 1u) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
						ee += im_dim;
					}
					if (MethodList.APLS) {
						vec.APLS_MLEM(span, iter + 1u) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
						ee += im_dim;
					}
					if (MethodList.TGV) {
						vec.TGV_MLEM(span, iter + 1u) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
						ee += im_dim;
					}
					if (MethodList.NLM) {
						vec.NLM_MLEM(span, iter + 1u) = vec.im_mlem(seq(ee, ee + im_dim - 1u));
						ee += im_dim;
					}
				}
				if (MethodList.CUSTOM) {
					if (MethodList.OSLMLEM) {
						vec.custom_MLEM = vec.im_os(seq(ee, ee + im_dim - 1u)).copy();
						ee += im_dim;
					}
					break_iter = true;
				}

				if (no_norm_mlem == 0u)
					no_norm_mlem = 1u;

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
	}

	// Transfer memory control of all variables that weren't used
	//unlock_AF_im_vectors(vec, MethodList, true, mlem_bool, osem_bool, 0u);

	// Clear OpenCL buffers
	status = clReleaseMemObject(d_z);
	if (status != CL_SUCCESS)
		std::cerr << getErrorString(status) << std::endl;
	status = clReleaseMemObject(d_x);
	if (status != CL_SUCCESS)
		std::cerr << getErrorString(status) << std::endl;
	status = clReleaseMemObject(d_y);
	if (status != CL_SUCCESS)
		std::cerr << getErrorString(status) << std::endl;
	status = clReleaseMemObject(d_atten);
	if (status != CL_SUCCESS)
		std::cerr << getErrorString(status) << std::endl;
	status = clReleaseMemObject(d_norm);
	if (status != CL_SUCCESS)
		std::cerr << getErrorString(status) << std::endl;
	status = clReleaseMemObject(d_pseudos);
	if (status != CL_SUCCESS)
		std::cerr << getErrorString(status) << std::endl;
	status = clReleaseMemObject(d_xcenter);
	if (status != CL_SUCCESS)
		std::cerr << getErrorString(status) << std::endl;
	status = clReleaseMemObject(d_ycenter);
	if (status != CL_SUCCESS)
		std::cerr << getErrorString(status) << std::endl;
	status = clReleaseMemObject(d_zcenter);
	if (status != CL_SUCCESS)
		std::cerr << getErrorString(status) << std::endl;
	if (osem_bool) {
		status = clReleaseMemObject(d_reko_type);
		if (status != CL_SUCCESS)
			std::cerr << getErrorString(status) << std::endl;
		for (uint32_t kk = 0u; kk < subsets; kk++) {
			status = clReleaseMemObject(d_lor[kk]);
			if (status != CL_SUCCESS)
				std::cerr << getErrorString(status) << std::endl;
			status = clReleaseMemObject(d_xyindex[kk]);
			if (status != CL_SUCCESS)
				std::cerr << getErrorString(status) << std::endl;
			status = clReleaseMemObject(d_zindex[kk]);
			if (status != CL_SUCCESS)
				std::cerr << getErrorString(status) << std::endl;
			status = clReleaseMemObject(d_L[kk]);
			if (status != CL_SUCCESS)
				std::cerr << getErrorString(status) << std::endl;
			status = clReleaseMemObject(d_Sino[kk]);
			if (status != CL_SUCCESS)
				std::cerr << getErrorString(status) << std::endl;
			//if (randoms_correction)
			status = clReleaseMemObject(d_sc_ra[kk]);
			if (status != CL_SUCCESS)
				std::cerr << getErrorString(status) << std::endl;
		}
	}
	if (mlem_bool) {
		status = clReleaseMemObject(d_lor_mlem);
		if (status != CL_SUCCESS)
			std::cerr << getErrorString(status) << std::endl;
		status = clReleaseMemObject(d_xyindex_mlem);
		if (status != CL_SUCCESS)
			std::cerr << getErrorString(status) << std::endl;
		status = clReleaseMemObject(d_zindex_mlem);
		if (status != CL_SUCCESS)
			std::cerr << getErrorString(status) << std::endl;
		status = clReleaseMemObject(d_L_mlem);
		if (status != CL_SUCCESS)
			std::cerr << getErrorString(status) << std::endl;
		status = clReleaseMemObject(d_Sino_mlem);
		if (status != CL_SUCCESS)
			std::cerr << getErrorString(status) << std::endl;
		//if (randoms_correction)
		status = clReleaseMemObject(d_sc_ra_mlem);
		if (status != CL_SUCCESS)
			std::cerr << getErrorString(status) << std::endl;
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