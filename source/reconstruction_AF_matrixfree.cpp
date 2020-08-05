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

	af::setDevice(device);

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
	cl_int status = CL_SUCCESS;

	// Create the OpenCL context and command queue and assign the device
	cl::Context af_context(afcl::getContext(true));
	std::vector<cl::Device> devices = af_context.getInfo<CL_CONTEXT_DEVICES>(&status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}
	cl::Device af_device_id = devices[0];
	cl::CommandQueue af_queue(afcl::getQueue(true));

	cl::Program program_os;
	cl::Program program_ml;
	cl::Program program_mbsrem;

	OpenCLStruct.af_queue = &af_queue;


	cl_ulong mem;
	cl_ulong mem_loc;
	mem = af_device_id.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>(&status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}
	if ((static_cast<cl_float>(mem) * mem_portions) > image_bytes && !MethodList.CUSTOM)
		compute_norm_matrix = 0u;

	mem_loc = af_device_id.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>(&status);
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
			Summ.assign(subsets, constant(0LL, im_dim, 1, s64));
		else
			Summ.assign(subsets, constant(0.f, im_dim, 1));
	}
	else {
		if (atomic_64bit)
			Summ.assign(1ULL, constant(0LL, im_dim, 1, s64));
		else
			Summ.assign(1ULL, constant(0.f, im_dim, 1));
	}

	// Create the kernels
	cl::Kernel kernel_ml, kernel, kernel_mramla;

	status = createKernels(kernel_ml, kernel, kernel_mramla, OpenCLStruct.kernelNLM, osem_bool, program_os, program_ml, program_mbsrem, MethodList, w_vec, projector_type, mlem_bool, precompute,
		n_rays, n_rays3D);
	if (status != CL_SUCCESS) {
		mexPrintf("Failed to create kernels\n");
		return;
	}
	//else if (verbose) {
	//	mexPrintf("OpenCL kernels successfully created\n");
	//	mexEvalString("pause(.0001);");
	//}

	// Create and write buffers
	cl::Buffer d_x, d_y, d_z, d_pseudos, d_atten, d_xcenter, d_ycenter, d_zcenter, d_norm_mlem, d_reko_type, d_reko_type_mlem, d_V, d_scat_mlem;
	cl::Buffer d_Summ, d_Summ_mlem;
	cl::Buffer d_lor_mlem, d_L_mlem, d_zindex_mlem, d_xyindex_mlem, d_Sino_mlem, d_sc_ra_mlem;

	std::vector<cl::Buffer> d_lor(subsets);
	std::vector<cl::Buffer> d_L(subsets);
	std::vector<cl::Buffer> d_zindex(subsets);
	std::vector<cl::Buffer> d_xyindex(subsets);
	std::vector<cl::Buffer> d_Sino(subsets);
	std::vector<cl::Buffer> d_sc_ra(subsets);
	std::vector<cl::Buffer> d_norm(subsets);
	std::vector<cl::Buffer> d_scat(subsets);

	float *apu = (float*)mxGetData(mxGetCell(Sin, 0));

	status = createAndWriteBuffers(d_x, d_y, d_z, d_lor, d_L, d_zindex, d_xyindex, d_Sino, d_sc_ra, size_x, size_z, TotSinos, size_atten, size_norm, size_scat, prows, 
		length, x, y, z_det, xy_index, z_index, lor1, L, apu, raw, af_context, subsets, pituus, atten, norm, scat, pseudos, V, af_queue, d_atten, d_norm, d_scat, d_pseudos, d_V, 
		d_xcenter, d_ycenter, d_zcenter, x_center, y_center, z_center, size_center_x, size_center_y, size_center_z, size_of_x, size_V, atomic_64bit, randoms_correction, 
		sc_ra, precompute, d_lor_mlem, d_L_mlem, d_zindex_mlem, d_xyindex_mlem, d_Sino_mlem, d_sc_ra_mlem, d_reko_type, d_reko_type_mlem, osem_bool, mlem_bool, koko,
		reko_type, reko_type_mlem, n_rekos, n_rekos_mlem, d_norm_mlem, d_scat_mlem);
	if (status != CL_SUCCESS) {
		mexPrintf("Buffer creation failed\n");
		return;
	}
	//else {
	//	mexPrintf("Buffer creation succeeded\n");
	//}

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
		status = kernel.setArg(kernelInd_OSEM++, global_factor);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		kernel.setArg(kernelInd_OSEM++, epps);
		kernel.setArg(kernelInd_OSEM++, im_dim);
		kernel.setArg(kernelInd_OSEM++, Nx);
		kernel.setArg(kernelInd_OSEM++, Ny);
		kernel.setArg(kernelInd_OSEM++, Nz);
		kernel.setArg(kernelInd_OSEM++, dz);
		kernel.setArg(kernelInd_OSEM++, dx);
		kernel.setArg(kernelInd_OSEM++, dy);
		kernel.setArg(kernelInd_OSEM++, bz);
		kernel.setArg(kernelInd_OSEM++, bx);
		kernel.setArg(kernelInd_OSEM++, by);
		kernel.setArg(kernelInd_OSEM++, bzb);
		kernel.setArg(kernelInd_OSEM++, maxxx);
		kernel.setArg(kernelInd_OSEM++, maxyy);
		kernel.setArg(kernelInd_OSEM++, zmax);
		kernel.setArg(kernelInd_OSEM++, NSlices);
		kernel.setArg(kernelInd_OSEM++, size_x);
		kernel.setArg(kernelInd_OSEM++, TotSinos);
		kernel.setArg(kernelInd_OSEM++, det_per_ring);
		kernel.setArg(kernelInd_OSEM++, prows);
		kernel.setArg(kernelInd_OSEM++, Nxy);
		kernel.setArg(kernelInd_OSEM++, fp);
		if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1))) {
			kernel.setArg(kernelInd_OSEM++, tube_width);
			kernel.setArg(kernelInd_OSEM++, crystal_size_z);
			kernel.setArg(kernelInd_OSEM++, bmin);
			kernel.setArg(kernelInd_OSEM++, bmax);
			kernel.setArg(kernelInd_OSEM++, Vmax);
		}
		else if (projector_type == 1u && !precompute) {
			kernel.setArg(kernelInd_OSEM++, dc_z);
			kernel.setArg(kernelInd_OSEM++, n_rays);
		}
	}

	if (mlem_bool) {
		// Set the kernel parameters that do not change
		status = kernel_ml.setArg(kernelInd_MLEM++, global_factor);
		if (status != CL_SUCCESS) {
			getErrorString(status);
			return;
		}
		kernel_ml.setArg(kernelInd_MLEM++, epps);
		kernel_ml.setArg(kernelInd_MLEM++, im_dim);
		kernel_ml.setArg(kernelInd_MLEM++, Nx);
		kernel_ml.setArg(kernelInd_MLEM++, Ny);
		kernel_ml.setArg(kernelInd_MLEM++, Nz);
		kernel_ml.setArg(kernelInd_MLEM++, dz);
		kernel_ml.setArg(kernelInd_MLEM++, dx);
		kernel_ml.setArg(kernelInd_MLEM++, dy);
		kernel_ml.setArg(kernelInd_MLEM++, bz);
		kernel_ml.setArg(kernelInd_MLEM++, bx);
		kernel_ml.setArg(kernelInd_MLEM++, by);
		kernel_ml.setArg(kernelInd_MLEM++, bzb);
		kernel_ml.setArg(kernelInd_MLEM++, maxxx);
		kernel_ml.setArg(kernelInd_MLEM++, maxyy);
		kernel_ml.setArg(kernelInd_MLEM++, zmax);
		kernel_ml.setArg(kernelInd_MLEM++, NSlices);
		kernel_ml.setArg(kernelInd_MLEM++, size_x);
		kernel_ml.setArg(kernelInd_MLEM++, TotSinos);
		kernel_ml.setArg(kernelInd_MLEM++, det_per_ring);
		kernel_ml.setArg(kernelInd_MLEM++, prows);
		kernel_ml.setArg(kernelInd_MLEM++, Nxy);
		kernel_ml.setArg(kernelInd_MLEM++, fp);
		if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1))) {
			kernel_ml.setArg(kernelInd_MLEM++, tube_width);
			kernel_ml.setArg(kernelInd_MLEM++, crystal_size_z);
			kernel_ml.setArg(kernelInd_MLEM++, bmin);
			kernel_ml.setArg(kernelInd_MLEM++, bmax);
			kernel_ml.setArg(kernelInd_MLEM++, Vmax);
		}
		else if (projector_type == 1u && !precompute) {
			kernel_ml.setArg(kernelInd_MLEM++, dc_z);
			kernel_ml.setArg(kernelInd_MLEM++, n_rays);
		}
	}

	if (mlem_bool) {
		if (atomic_64bit)
			Summ_mlem = constant(0LL, im_dim, 1, s64);
		else
			Summ_mlem = constant(0.f, im_dim, 1);
	}

	// Loop through each time-step
	for (uint32_t tt = t0; tt < Nt; tt++) {
	// Compute the prepass phase for MRAMLA, MBSREM, RBI, COSEM, ACOSEM or ECOSEM if applicable
	if (((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI) && w_vec.MBSREM_prepass ||
		MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0) && (!MethodList.CUSTOM || osa_iter0 == 0u)) {

		// Set the kernel parameters that do not change
		kernel_mramla.setArg(kernelInd_MRAMLA++, global_factor);
		kernel_mramla.setArg(kernelInd_MRAMLA++, epps);
		kernel_mramla.setArg(kernelInd_MRAMLA++, im_dim);
		kernel_mramla.setArg(kernelInd_MRAMLA++, Nx);
		kernel_mramla.setArg(kernelInd_MRAMLA++, Ny);
		kernel_mramla.setArg(kernelInd_MRAMLA++, Nz);
		kernel_mramla.setArg(kernelInd_MRAMLA++, dz);
		kernel_mramla.setArg(kernelInd_MRAMLA++, dx);
		kernel_mramla.setArg(kernelInd_MRAMLA++, dy);
		kernel_mramla.setArg(kernelInd_MRAMLA++, bz);
		kernel_mramla.setArg(kernelInd_MRAMLA++, bx);
		kernel_mramla.setArg(kernelInd_MRAMLA++, by);
		kernel_mramla.setArg(kernelInd_MRAMLA++, bzb);
		kernel_mramla.setArg(kernelInd_MRAMLA++, maxxx);
		kernel_mramla.setArg(kernelInd_MRAMLA++, maxyy);
		kernel_mramla.setArg(kernelInd_MRAMLA++, zmax);
		kernel_mramla.setArg(kernelInd_MRAMLA++, NSlices);
		kernel_mramla.setArg(kernelInd_MRAMLA++, size_x);
		kernel_mramla.setArg(kernelInd_MRAMLA++, TotSinos);
		kernel_mramla.setArg(kernelInd_MRAMLA++, det_per_ring);
		kernel_mramla.setArg(kernelInd_MRAMLA++, prows);
		kernel_mramla.setArg(kernelInd_MRAMLA++, Nxy);
		kernel_mramla.setArg(kernelInd_MRAMLA++, fp);
		if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1))) {
			kernel_mramla.setArg(kernelInd_MRAMLA++, tube_width);
			kernel_mramla.setArg(kernelInd_MRAMLA++, crystal_size_z);
			kernel_mramla.setArg(kernelInd_MRAMLA++, bmin);
			kernel_mramla.setArg(kernelInd_MRAMLA++, bmax);
			kernel_mramla.setArg(kernelInd_MRAMLA++, Vmax);
		}
		else if (projector_type == 1u && !precompute) {
			kernel_mramla.setArg(kernelInd_MRAMLA++, dc_z);
			kernel_mramla.setArg(kernelInd_MRAMLA++, n_rays);
		}
		kernel_mramla.setArg(kernelInd_MRAMLA++, w_vec.epsilon_mramla);
		kernel_mramla.setArg(kernelInd_MRAMLA++, d_atten);
		kernel_mramla.setArg(kernelInd_MRAMLA++, d_pseudos);
		kernel_mramla.setArg(kernelInd_MRAMLA++, d_x);
		kernel_mramla.setArg(kernelInd_MRAMLA++, d_y);
		kernel_mramla.setArg(kernelInd_MRAMLA++, d_z);
		if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1))) {
			kernel_mramla.setArg(kernelInd_MRAMLA++, d_xcenter);
			kernel_mramla.setArg(kernelInd_MRAMLA++, d_ycenter);
			kernel_mramla.setArg(kernelInd_MRAMLA++, d_zcenter);
			kernel_mramla.setArg(kernelInd_MRAMLA++, d_V);
		}
		kernel_mramla.setArg(kernelInd_MRAMLA++, d_reko_type);

		uint32_t alku = 0u;

		if ((MethodList.MRAMLA || MethodList.MBSREM) && w_vec.MBSREM_prepass)
			w_vec.Amin = constant(0.f, koko, 1);

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
					af_queue.enqueueWriteBuffer(d_Sino[kk], CL_TRUE, 0, sizeof(float) * length[kk], &apu[pituus[kk]]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
					if (randoms_correction) {
						float* ra_apu = (float*)mxGetData(mxGetCell(sc_ra, tt));
						af_queue.enqueueWriteBuffer(d_sc_ra[kk], CL_TRUE, 0, sizeof(float)* length[kk], &ra_apu[pituus[kk]]);
					}
					else
						af_queue.enqueueWriteBuffer(d_sc_ra[kk], CL_TRUE, 0, sizeof(cl_float), &zerof);
						//status = clEnqueueFillBuffer(af_queue, d_sc_ra[kk], &zerof, sizeof(cl_float), 0, sizeof(cl_float), 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
					if (scatter == 1u) {
						scat = (float*)mxGetData(mxGetCell(mxGetField(options, 0, "ScatterC"), tt));
						af_queue.enqueueWriteBuffer(d_scat[kk], CL_TRUE, 0, sizeof(float)* length[kk], &scat[pituus[kk]]);
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
				af_queue.enqueueWriteBuffer(d_Sino_mlem, CL_TRUE, 0, sizeof(float) * koko, apu);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				if (randoms_correction) {
					float* ra_apu = (float*)mxGetData(mxGetCell(sc_ra, tt));
					af_queue.enqueueWriteBuffer(d_sc_ra_mlem, CL_TRUE, 0, sizeof(float)* koko, ra_apu);
				}
				else
					af_queue.enqueueWriteBuffer(d_sc_ra_mlem, CL_TRUE, 0, sizeof(cl_float), &zerof);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				if (scatter == 1u) {
					scat = (float*)mxGetData(mxGetCell(mxGetField(options, 0, "ScatterC"), tt));
					af_queue.enqueueWriteBuffer(d_scat_mlem, CL_TRUE, 0, sizeof(float)* koko, scat);
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
			kernel.setArg(kernelInd_OSEMTIter++, w_vec.epsilon_mramla);
			kernel.setArg(kernelInd_OSEMTIter++, d_atten);
			kernel.setArg(kernelInd_OSEMTIter++, d_pseudos);
			kernel.setArg(kernelInd_OSEMTIter++, d_x);
			kernel.setArg(kernelInd_OSEMTIter++, d_y);
			kernel.setArg(kernelInd_OSEMTIter++, d_z);
			if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1))) {
				kernel.setArg(kernelInd_OSEMTIter++, d_xcenter);
				kernel.setArg(kernelInd_OSEMTIter++, d_ycenter);
				kernel.setArg(kernelInd_OSEMTIter++, d_zcenter);
				kernel.setArg(kernelInd_OSEMTIter++, d_V);
			}
			kernel.setArg(kernelInd_OSEMTIter++, d_reko_type);
		}
		if (mlem_bool) {
			kernel_ml.setArg(kernelInd_MLEMT++, w_vec.epsilon_mramla);
			kernel_ml.setArg(kernelInd_MLEMT++, d_atten);
			kernel_ml.setArg(kernelInd_MLEMT++, d_pseudos);
			kernel_ml.setArg(kernelInd_MLEMT++, d_x);
			kernel_ml.setArg(kernelInd_MLEMT++, d_y);
			kernel_ml.setArg(kernelInd_MLEMT++, d_z);
			if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1))) {
				kernel_ml.setArg(kernelInd_MLEMT++, d_xcenter);
				kernel_ml.setArg(kernelInd_MLEMT++, d_ycenter);
				kernel_ml.setArg(kernelInd_MLEMT++, d_zcenter);
				kernel_ml.setArg(kernelInd_MLEMT++, d_V);
			}
			kernel_ml.setArg(kernelInd_MLEMT++, d_reko_type_mlem);
		}


		status = af_queue.finish();
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
							Summ[0] = constant(0LL, im_dim, 1, s64);
						}
						else
							Summ[0] = constant(0.f, im_dim, 1);
						d_Summ = cl::Buffer(*Summ[0].device<cl_mem>(), true);
						//d_Summ->operator=(*Summ[0].device<cl_mem>());
					}
					else {
						if (no_norm == 0u && tt > 0u) {
							if (atomic_64bit)
								Summ[osa_iter] = constant(0LL, im_dim, 1, s64);
							else
								Summ[osa_iter] = constant(0.f, im_dim, 1);
							//d_Summ->operator=(*Summ[osa_iter].device<cl_mem>());
							d_Summ = cl::Buffer(*Summ[osa_iter].device<cl_mem>(), true);
						}
						else if (no_norm == 1u) {
							if (atomic_64bit)
								apu_sum = constant(0LL, 1, 1, s64);
							else
								apu_sum = constant(0.f, 1, 1, f32);
							//d_Summ->operator=(*apu_sum.device<cl_mem>());
							d_Summ = cl::Buffer(*apu_sum.device<cl_mem>(), true);
						}
						else
							//d_Summ->operator=(*Summ[osa_iter].device<cl_mem>());
							d_Summ = cl::Buffer(*Summ[osa_iter].device<cl_mem>(), true);
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
					kernel.setArg(kernelInd_OSEMSubIter++, d_norm[osa_iter]);
					kernel.setArg(kernelInd_OSEMSubIter++, d_scat[osa_iter]);
					kernel.setArg(kernelInd_OSEMSubIter++, d_Summ);
					kernel.setArg(kernelInd_OSEMSubIter++, d_lor[osa_iter]);
					kernel.setArg(kernelInd_OSEMSubIter++, d_xyindex[osa_iter]);
					kernel.setArg(kernelInd_OSEMSubIter++, d_zindex[osa_iter]);
					kernel.setArg(kernelInd_OSEMSubIter++, d_L[osa_iter]);
					kernel.setArg(kernelInd_OSEMSubIter++, d_Sino[osa_iter]);
					kernel.setArg(kernelInd_OSEMSubIter++, d_sc_ra[osa_iter]);
					kernel.setArg(kernelInd_OSEMSubIter++, vec_opencl.d_im_os);
					kernel.setArg(kernelInd_OSEMSubIter++, vec_opencl.d_rhs_os);
					kernel.setArg(kernelInd_OSEMSubIter++, no_norm);
					kernel.setArg(kernelInd_OSEMSubIter++, m_size);
					kernel.setArg(kernelInd_OSEMSubIter++, st);
					cl::NDRange local(local_size);
					cl::NDRange global(global_size);
					// Compute the kernel
					status = af_queue.enqueueNDRangeKernel(kernel, cl::NullRange, global, local);

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
					status = af_queue.finish();
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

					status = af_queue.finish();

					if (break_iter)
						break;

				}
				//vec.im_os = vec.rhs_os;

				computeOSEstimatesIter(vec, w_vec, MethodList, im_dim, epps, iter, osa_iter0, subsets, beta, Nx, Ny, Nz, data, n_rekos2, OpenCLStruct);
				
				if (use_psf && w_vec.deconvolution && osem_bool) {
					computeDeblur(vec, g, Nx, Ny, Nz, w_vec, MethodList, iter, deblur_iterations, epps);
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
							apu_sum_mlem = constant(0LL, 1, 1, s64);
						else
							apu_sum_mlem = constant(0.f, 1, 1, f32);
						d_Summ_mlem = cl::Buffer(*apu_sum_mlem.device<cl_mem>(), true);
					}
					else
						d_Summ_mlem = cl::Buffer(*Summ_mlem.device<cl_mem>(), true);
				}
				else {
					if (no_norm_mlem == 1u) {
						if (atomic_64bit)
							apu_sum_mlem = constant(0LL, 1, 1, s64);
						else
							apu_sum_mlem = constant(0.f, 1, 1, f32);
						d_Summ_mlem = cl::Buffer(*apu_sum_mlem.device<cl_mem>(), true);
					}
					else
						d_Summ_mlem = cl::Buffer(*Summ_mlem.device<cl_mem>(), true);
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

				cl::NDRange local(local_size);
				cl::NDRange global(global_size);

				kernel_ml.setArg(kernelInd_MLEMSubIter++, d_norm_mlem);
				kernel_ml.setArg(kernelInd_MLEMSubIter++, d_scat_mlem);
				kernel_ml.setArg(kernelInd_MLEMSubIter++, d_Summ_mlem);
				kernel_ml.setArg(kernelInd_MLEMSubIter++, d_lor_mlem);
				kernel_ml.setArg(kernelInd_MLEMSubIter++, d_xyindex_mlem);
				kernel_ml.setArg(kernelInd_MLEMSubIter++, d_zindex_mlem);
				kernel_ml.setArg(kernelInd_MLEMSubIter++, d_L_mlem);
				kernel_ml.setArg(kernelInd_MLEMSubIter++, d_Sino_mlem);
				kernel_ml.setArg(kernelInd_MLEMSubIter++, d_sc_ra_mlem);
				kernel_ml.setArg(kernelInd_MLEMSubIter++, vec_opencl.d_im_mlem);
				kernel_ml.setArg(kernelInd_MLEMSubIter++, vec_opencl.d_rhs_mlem);
				kernel_ml.setArg(kernelInd_MLEMSubIter++, no_norm_mlem);
				kernel_ml.setArg(kernelInd_MLEMSubIter++, m_size);
				kernel_ml.setArg(kernelInd_MLEMSubIter++, st);
				status = af_queue.enqueueNDRangeKernel(kernel_ml, cl::NullRange, global, local);

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
				status = af_queue.finish();
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
					computeDeblurMLEM(vec, g, Nx, Ny, Nz, w_vec, MethodList, iter, deblur_iterations, epps);
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

	status = af_queue.finish();
	af::sync();
	af::deviceGC();

	return;
}