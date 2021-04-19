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
	const float NSlices, const int64_t* pituus, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t size_x, const uint32_t TotSinos, 
	mxArray* cell, const mwSize* dimmi, const bool verbose, const uint32_t randoms_correction, const uint32_t attenuation_correction,
	const uint32_t normalization, const float* atten, const size_t size_atten, const float* norm, const size_t size_norm, const uint32_t subsets,
	const float epps, const char* k_path, const uint32_t Nt, const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L, 
	const uint8_t raw, const size_t size_z, const bool osem_bool, const char* fileName, const bool use_psf, const float tube_width, 
	const float crystal_size_z, const float* x_center, const float* y_center, const float* z_center, const size_t size_center_x, const size_t size_center_y, 
	const size_t size_of_x, const size_t size_center_z, const uint32_t projector_type, const char* header_directory, const bool precompute, 
	const uint32_t device, const int32_t dec, const uint16_t n_rays, const uint16_t n_rays3D, const float cr_pz, const bool use_64bit_atomics, uint32_t n_rekos,
	const uint32_t n_rekos_mlem, const uint8_t* reko_type, const uint8_t* reko_type_mlem, const float global_factor, const float bmin, const float bmax, 
	const float Vmax, const float* V, const size_t size_V, const float* gaussian, const size_t size_gauss, const bool saveIter, const bool TOF, const int64_t TOFSize, 
	const float sigma_x, const float* TOFCenter, const int64_t nBins) {

	af::setDevice(device);

	// Number of voxels
	const uint32_t Nxy = Nx * Ny;
	const uint32_t im_dim = Nxy * Nz;
	const cl_uchar fp = 0;

	// Distance between rays in multi-ray Siddon
	const float dc_z = cr_pz / static_cast<float>(n_rays3D + 1);

	bool break_iter = false;

	bool loadTOF = true;

	uint32_t t0 = 0u;
	uint32_t iter0 = 0u;
	uint32_t osa_iter0 = 0u;

	// Hard-coded local size
	uint64_t local_size = 64ULL;

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
	//matlabArrays ArrayList;
	// Create a struct containing the reconstruction methods used
	RecMethods MethodList;
	// Same as above, but as cl_ variables
	RecMethodsOpenCL MethodListOpenCL;

	kernelStruct OpenCLStruct;

	// Obtain the reconstruction methods used
	get_rec_methods(options, MethodList);

	const uint8_t listmode = (uint8_t)mxGetScalar(mxGetField(options, 0, "listmode"));
	const bool computeSensImag = (bool)mxGetScalar(mxGetField(options, 0, "compute_sensitivity_image"));
	const bool CT = (bool)mxGetScalar(mxGetField(options, 0, "CT"));
	const bool atomic_32bit = (bool)mxGetScalar(mxGetField(options, 0, "use_32bit_atomics"));

	if (listmode == 2)
		MethodList.MLEM = true;
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
	std::string deviceName = af_device_id.getInfo<CL_DEVICE_VENDOR>(&status);
	std::string NV("NVIDIA Corporation");
	if (NV.compare(deviceName) == 0 && projector_type == 1)
		local_size = 32ULL;

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
	if (((static_cast<cl_float>(mem) * mem_portions) > image_bytes && !MethodList.CUSTOM) || (listmode == 1 && computeSensImag))
		compute_norm_matrix = 0u;

	mem_loc = af_device_id.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>(&status);
	if (status != CL_SUCCESS) {
		getErrorString(status);
		return;
	}

	// Create the MATLAB output arrays
	//create_matlab_output(ArrayList, dimmi, MethodList, 4);

	// Initial value
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	array x00(Nx * Ny * Nz, (float*)mxGetSingles(mxGetField(options, 0, "x0")), afHost);
#else
	array x00(Nx * Ny * Nz, (float*)mxGetData(mxGetField(options, 0, "x0")), afHost);
#endif


	array pj3, apu_sum, E, apu_sum_mlem;

	if ((MethodList.MRAMLA || MethodList.MBSREM) && Nt > 1U)
		E = constant(1.f, koko * nBins, 1);
	else
		E = constant(0.f, 1, 1);

	// Are ML-methods used?
	bool mlem_bool = (MethodList.MLEM || MethodList.OSLMLEM) ? true : false;

	// Number of measurements at each subset
	std::vector<size_t> length(subsets);

	for (uint32_t kk = 0; kk < subsets; kk++)
		length[kk] = pituus[kk + 1u] - pituus[kk];

	uint32_t subsetsUsed = subsets;
	// For custom prior
	if (MethodList.CUSTOM) {
		osa_iter0 = (uint32_t)mxGetScalar(mxGetField(options, 0, "osa_iter"));
		iter0 = (uint32_t)mxGetScalar(mxGetField(options, 0, "iter"));
		t0 = (uint32_t)mxGetScalar(mxGetField(options, 0, "tt"));
	}

	// Struct containing ArrayFire arrays containing the image estimates and other necessary vectors/matrices
	AF_im_vectors vec;
	// vector containing beta-values for the MAP-methods
	//Beta beta;
	std::vector<float> beta;
	// Struct containing the necessary variables for the priors
	Weighting w_vec;
	// Struct for TV data
	TVdata data;
	// Struct containing the OpenCL image estimates
	OpenCL_im_vectors vec_opencl;

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
		for (uint32_t kk = 0U; kk < n_rekos2; kk++) {
			vec.im_os(seq(kk * im_dim, (kk + 1) * im_dim - 1)) = x00;
		}
	}

	if (mlem_bool) {
		vec.im_mlem = constant(0.f, im_dim * n_rekos_mlem, 1);
		for (uint32_t kk = 0U; kk < n_rekos_mlem; kk++) {
			vec.im_mlem(seq(kk * im_dim, (kk + 1) * im_dim - 1)) = x00;
		}
	}
	if (DEBUG) {
		mexPrintf("n_rekos_mlem = %u\n", n_rekos_mlem);
		mexPrintf("n_rekos2 = %u\n", n_rekos2);
		mexPrintf("n_rekos = %u\n", n_rekos);
		mexPrintf("koko = %u\n", koko);
		mexPrintf("nBins = %u\n", nBins);
		mexEvalString("pause(.0001);");
	}

	// Load the necessary data from the MATLAB input and form the necessary variables
	form_data_variables(vec, beta, w_vec, options, Nx, Ny, Nz, Niter, x00, im_dim, koko, MethodList, data, subsets, osa_iter0, use_psf, saveIter, Nt, iter0, CT);

	// Power factor for ACOSEM
	w_vec.h_ACOSEM_2 = 1.f / w_vec.h_ACOSEM;

	if (MethodList.CUSTOM) {
		if (w_vec.MBSREM_prepass && osa_iter0 == 0U) {
			subsetsUsed = subsets;
		}
		else
			subsetsUsed = osa_iter0 + 1U;
	}
	if (DEBUG) {
		mexPrintf("nProjections = %d\n", w_vec.nProjections);
		mexPrintf("size_y = %u\n", w_vec.size_y);
		mexPrintf("subsetsUsed = %u\n", subsetsUsed);
		mexPrintf("CT = %u\n", CT);
		mexPrintf("listmode = %u\n", listmode);
		mexEvalString("pause(.0001);");
	}

	float* scat = nullptr;
	const uint32_t scatter = static_cast<uint32_t>((bool)mxGetScalar(mxGetField(options, 0, "scatter")));
	size_t size_scat = 1ULL;
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	scat = (float*)mxGetSingles(mxGetCell(mxGetField(options, 0, "ScatterC"), 0));
#else
	scat = (float*)mxGetData(mxGetCell(mxGetField(options, 0, "ScatterC"), 0));
#endif
	if (scatter == 1U) {
		size_scat = mxGetNumberOfElements(mxGetCell(mxGetField(options, 0, "ScatterC"), 0));
	}

	status = createProgram(verbose, k_path, af_context, af_device_id, fileName, program_os, program_ml, program_mbsrem, atomic_64bit, atomic_32bit, device, header_directory,
		projector_type, crystal_size_z, precompute, raw, attenuation_correction, normalization, dec, local_size, n_rays, n_rays3D, false, MethodList, osem_bool, 
		mlem_bool, n_rekos2, n_rekos_mlem, w_vec, osa_iter0, cr_pz, dx, use_psf, scatter, randoms_correction, TOF, nBins, listmode, CT);
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
		else if (atomic_32bit && !w_vec.MBSREM_prepass)
			Summ.assign(subsets, constant(0, im_dim, 1, s32));
		else
			Summ.assign(subsets, constant(0.f, im_dim, 1));
	}
	else {
		if (atomic_64bit)
			Summ.assign(1ULL, constant(0LL, im_dim, 1, s64));
		else if (atomic_32bit)
			Summ.assign(1ULL, constant(0, im_dim, 1, s32));
		else
			Summ.assign(1ULL, constant(0.f, im_dim, 1));
	}

	if (DEBUG) {
		mexPrintf("Summ[0] = %f\n", af::sum<float>(Summ[0]));
		mexPrintf("w_vec.MBSREM_prepass = %u\n", w_vec.MBSREM_prepass);
		mexEvalString("pause(.0001);");
	}

	// Create the kernels
	cl::Kernel kernel_ml, kernel, kernel_mramla;

	status = createKernels(kernel_ml, kernel, kernel_mramla, OpenCLStruct.kernelNLM, OpenCLStruct.kernelMed, osem_bool, program_os, program_ml, program_mbsrem, MethodList, w_vec, projector_type,
		mlem_bool, precompute, n_rays, n_rays3D);
	if (status != CL_SUCCESS) {
		mexPrintf("Failed to create kernels\n");
		return;
	}
	else if (DEBUG) {
		mexPrintf("OpenCL kernels successfully created\n");
		mexEvalString("pause(.0001);");
	}

	if (static_cast<double>(mem) * 0.75 < static_cast<double>(koko * nBins * sizeof(float)) && TOF)
		loadTOF = false;

	uint32_t TOFsubsets = subsets;
	if (!loadTOF)
		TOFsubsets = 1U;

	// Create and write buffers
	cl::Buffer d_x, d_y, d_z, d_pseudos, d_atten, d_xcenter, d_ycenter, d_zcenter, d_norm_mlem, d_reko_type, d_reko_type_mlem, d_V, d_scat_mlem, d_TOFCenter;
	cl::Buffer d_Summ, d_Summ_mlem;
	cl::Buffer d_lor_mlem, d_L_mlem, d_zindex_mlem, d_xyindex_mlem, d_Sino_mlem, d_sc_ra_mlem;
	cl::Buffer d_angles;

	std::vector<cl::Buffer> d_lor(subsets);
	std::vector<cl::Buffer> d_L(subsets);
	std::vector<cl::Buffer> d_zindex(subsets);
	std::vector<cl::Buffer> d_xyindex(subsets);
	std::vector<cl::Buffer> d_Sino(TOFsubsets);
	std::vector<cl::Buffer> d_sc_ra(subsets);
	std::vector<cl::Buffer> d_norm(subsets);
	std::vector<cl::Buffer> d_scat(subsets);

#ifdef MX_HAS_INTERLEAVED_COMPLEX
	float* apu = (float*)mxGetSingles(mxGetCell(Sin, 0));
#else
	float* apu = (float*)mxGetData(mxGetCell(Sin, 0));
#endif
	size_t ind = mxGetNumberOfElements(mxGetCell(Sin, 0));
	float* angles = nullptr;
	if (CT)
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		angles = (float*)mxGetSingles(mxGetField(options, 0, "angles"));
#else
		angles = (float*)mxGetData(mxGetField(options, 0, "angles"));
#endif


	if (DEBUG && !CT) {
		mexPrintf("ind = %d\n", ind);
		mexPrintf("TOFsubsets = %d\n", TOFsubsets);
		mexPrintf("TOF = %d\n", TOF);
		mexEvalString("pause(.0001);");
	}


	if (DEBUG && CT) {
		//for (int uu = 0; uu < w_vec.nProjections; uu++)
		//	mexPrintf("angles = %f\n", angles[uu]);
		mexPrintf("nProjections = %d\n", w_vec.nProjections);
		mexPrintf("size_y = %d\n", w_vec.size_y);
		mexPrintf("dPitch = %f\n", w_vec.dPitch);
		mexEvalString("pause(.0001);");
	}

	status = createAndWriteBuffers(d_x, d_y, d_z, d_angles, d_lor, d_L, d_zindex, d_xyindex, d_Sino, d_sc_ra, size_x, size_z, TotSinos, size_atten, size_norm, size_scat, prows,
		length, x, y, z_det, xy_index, z_index, lor1, L, apu, raw, af_context, subsets, pituus, atten, norm, scat, pseudos, V, af_queue, d_atten, d_norm, d_scat, d_pseudos, d_V, 
		d_xcenter, d_ycenter, d_zcenter, x_center, y_center, z_center, size_center_x, size_center_y, size_center_z, size_of_x, size_V, atomic_64bit, atomic_32bit, randoms_correction,
		sc_ra, precompute, d_lor_mlem, d_L_mlem, d_zindex_mlem, d_xyindex_mlem, d_Sino_mlem, d_sc_ra_mlem, d_reko_type, d_reko_type_mlem, osem_bool, mlem_bool, koko,
		reko_type, reko_type_mlem, n_rekos, n_rekos_mlem, d_norm_mlem, d_scat_mlem, angles, TOF, nBins, loadTOF, d_TOFCenter, TOFCenter, subsetsUsed, osa_iter0, listmode, CT);
	if (status != CL_SUCCESS) {
		mexPrintf("Buffer creation failed\n");
		return;
	}
	else if (DEBUG) {
		mexPrintf("Buffer creation succeeded\n");
		mexEvalString("pause(.0001);");
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
		if (DEBUG) {
			mexPrintf("im_dim = %d\n", im_dim);
			mexEvalString("pause(.0001);");
		}
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
		kernel.setArg(kernelInd_OSEM++, sigma_x);
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
		kernel_ml.setArg(kernelInd_MLEM++, sigma_x);
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
		else if (atomic_32bit)
			Summ_mlem = constant(0, im_dim, 1, s32);
		else
			Summ_mlem = constant(0.f, im_dim, 1);
	}

	if (!osem_bool && w_vec.MBSREM_prepass)
		w_vec.MBSREM_prepass = false;

	// Loop through each time-step
	for (uint32_t tt = t0; tt < Nt; tt++) {
	// Compute the prepass phase for MRAMLA, MBSREM, RBI, PKMA, COSEM, ACOSEM or ECOSEM if applicable
	if (((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI || MethodList.PKMA) && w_vec.MBSREM_prepass ||
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
		kernel_mramla.setArg(kernelInd_MRAMLA++, sigma_x);
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
		kernel_mramla.setArg(kernelInd_MRAMLA++, d_TOFCenter);
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
		if (CT) {
			kernel_mramla.setArg(kernelInd_MRAMLA++, subsets);
			kernel_mramla.setArg(kernelInd_MRAMLA++, d_angles);
			kernel_mramla.setArg(kernelInd_MRAMLA++, w_vec.size_y);
			kernel_mramla.setArg(kernelInd_MRAMLA++, w_vec.dPitch);
			kernel_mramla.setArg(kernelInd_MRAMLA++, w_vec.nProjections);
		}

		uint32_t alku = 0u;

#ifdef MX_HAS_INTERLEAVED_COMPLEX
		const uint64_t* randSize = (uint64_t*)mxGetUint64s(mxGetField(options, 0, "randSize"));
#else
		const uint64_t* randSize = (uint64_t*)mxGetData(mxGetField(options, 0, "randSize"));
#endif

		if ((MethodList.MRAMLA || MethodList.MBSREM) && w_vec.MBSREM_prepass && Nt > 1U)
			w_vec.Amin = constant(0.f, koko * nBins, 1);

		// Run the prepass phase
		MRAMLA_prepass(subsets, im_dim, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, w_vec, Summ, d_Sino, koko, x00, vec.C_co,
			vec.C_aco, vec.C_osl, alku, kernel_mramla, d_L, raw, MethodListOpenCL, length, atomic_64bit, atomic_32bit, compute_norm_matrix, d_sc_ra, kernelInd_MRAMLA, E,
			d_norm, d_scat, use_psf, g, Nx, Ny, Nz, epps, TOF, loadTOF, Sin, nBins, koko, randoms_correction, local_size, randSize, Nt, CT);


		if ((MethodList.MRAMLA || MethodList.MBSREM || MethodList.PKMA) && tt == 0) {
			pj3 = w_vec.D / static_cast<float>(subsets);
		}

		if (verbose) {
			mexPrintf("MRAMLA & COSEM prepass completed\n");
			mexEvalString("pause(.0001);");
		}
	}
	// For custom prior only
	if ((MethodList.MRAMLA || MethodList.MBSREM || MethodList.PKMA) && MethodList.CUSTOM && (iter0 > 0 || osa_iter0 > 0)) {
		pj3 = w_vec.D / static_cast<float>(subsets);
		if (DEBUG) {
			mexPrintf("pj3 = %f\n", af::sum<float>(pj3));
			mexPrintf("w_vec.D = %f\n", af::sum<float>(w_vec.D));
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

		// Load precomputed sensitivity image
		if (computeSensImag && listmode == 1) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			float* apu = (float*)mxGetSingles(mxGetField(options, 0, "Summ"));
#else
			float* apu = (float*)mxGetData(mxGetField(options, 0, "Summ"));
#endif
			if (mlem_bool) {
				Summ_mlem = array(im_dim, 1, apu);
				if (use_psf) {
					Summ_mlem = computeConvolution(Summ_mlem, g, Nx, Ny, Nz, w_vec, 1u);
					af::sync();
				}
				no_norm_mlem = 1u;
			}
			if (osem_bool) {
				for (uint32_t osa_iter = 0; osa_iter < subsets; osa_iter++) {
					Summ[osa_iter] = array(im_dim, 1, apu) / static_cast<float>(subsets);
					if (use_psf) {
						Summ[osa_iter] = computeConvolution(Summ[osa_iter], g, Nx, Ny, Nz, w_vec, 1u);
						af::sync();
					}
				}
				no_norm = 1u;
			}
		}

		if (listmode == 2)
			af_queue.enqueueFillBuffer(d_Sino_mlem, zerof, 0, sizeof(cl_float));

		// Load the measurement and randoms data from the cell arrays
		if (tt > 0u) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			float* apu = (float*)mxGetSingles(mxGetCell(Sin, tt));
#else
			float* apu = (float*)mxGetData(mxGetCell(Sin, tt));
#endif
			if (osem_bool) {
				for (uint32_t kk = osa_iter0; kk < subsetsUsed; kk++) {
					if (TOF) {
						if (!loadTOF && kk == 0) {
							d_Sino[kk] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[kk] * nBins, NULL, &status);
							for (int64_t to = 0LL; to < nBins; to++)
								status = af_queue.enqueueWriteBuffer(d_Sino[kk], CL_FALSE, sizeof(float) * length[kk] * to, sizeof(float) * length[kk], &apu[pituus[kk] + koko * to]);
						}
						else if (loadTOF) {
							for (int64_t to = 0LL; to < nBins; to++)
								status = af_queue.enqueueWriteBuffer(d_Sino[kk], CL_FALSE, sizeof(float) * length[kk] * to, sizeof(float) * length[kk], &apu[pituus[kk] + koko * to]);
						}
					}
					else
						af_queue.enqueueWriteBuffer(d_Sino[kk], CL_TRUE, 0, sizeof(float) * length[kk], &apu[pituus[kk]]);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
					if (randoms_correction) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
						float* ra_apu = (float*)mxGetSingles(mxGetCell(sc_ra, tt));
#else
						float* ra_apu = (float*)mxGetData(mxGetCell(sc_ra, tt));
#endif
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
#ifdef MX_HAS_INTERLEAVED_COMPLEX
						scat = (float*)mxGetSingles(mxGetCell(mxGetField(options, 0, "ScatterC"), tt));
#else
						scat = (float*)mxGetData(mxGetCell(mxGetField(options, 0, "ScatterC"), tt));
#endif
						af_queue.enqueueWriteBuffer(d_scat[kk], CL_TRUE, 0, sizeof(float)* length[kk], &scat[pituus[kk]]);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
					}
				}
				vec.im_os = constant(0.f, im_dim * n_rekos2, 1);
				for (uint32_t kk = 0U; kk < n_rekos2; kk++) {
					vec.im_os(seq(kk * im_dim, (kk + 1) * im_dim - 1)) = x00;
				}
			}
			if (mlem_bool) {
				af_queue.enqueueWriteBuffer(d_Sino_mlem, CL_TRUE, 0, sizeof(float) * koko * nBins, apu);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				if (randoms_correction) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
					float* ra_apu = (float*)mxGetSingles(mxGetCell(sc_ra, tt));
#else
					float* ra_apu = (float*)mxGetData(mxGetCell(sc_ra, tt));
#endif
					af_queue.enqueueWriteBuffer(d_sc_ra_mlem, CL_TRUE, 0, sizeof(float)* koko, ra_apu);
				}
				else
					af_queue.enqueueWriteBuffer(d_sc_ra_mlem, CL_TRUE, 0, sizeof(cl_float), &zerof);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					return;
				}
				if (scatter == 1u) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
					scat = (float*)mxGetSingles(mxGetCell(mxGetField(options, 0, "ScatterC"), tt));
#else
					scat = (float*)mxGetData(mxGetCell(mxGetField(options, 0, "ScatterC"), tt));
#endif
					af_queue.enqueueWriteBuffer(d_scat_mlem, CL_TRUE, 0, sizeof(float)* koko, scat);
					if (status != CL_SUCCESS) {
						getErrorString(status);
						return;
					}
				}
				vec.im_mlem = constant(0.f, im_dim * n_rekos_mlem, 1);
				for (uint32_t kk = 0U; kk < n_rekos_mlem; kk++) {
					vec.im_mlem(seq(kk * im_dim, (kk + 1) * im_dim - 1)) = x00;
				}
			}
		}

		// Compute values needed for MBSREM and MRAMLA
		if ((MethodList.MBSREM || MethodList.MRAMLA) && Nt > 1U) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			array Sino = array(koko, (float*)mxGetSingles(mxGetCell(Sin, tt)), afHost);
#else
			array Sino = array(koko, (float*)mxGetData(mxGetCell(Sin, tt)), afHost);
#endif
			array rand;
			if (randoms_correction)
#ifdef MX_HAS_INTERLEAVED_COMPLEX
				rand = array(koko, (float*)mxGetSingles(mxGetCell(sc_ra, tt)), afHost);
#else
				rand = array(koko, (float*)mxGetData(mxGetCell(sc_ra, tt)), afHost);
#endif
			if (w_vec.U == 0.f) {
				const array Aind = w_vec.Amin > 0.f;
				if (CT)
					w_vec.U = max<float>(-af::log(Sino(Aind)) / w_vec.Amin(Aind));
				else
					w_vec.U = max<float>(Sino(Aind) / w_vec.Amin(Aind));
			}
			w_vec.epsilon_mramla = MBSREM_epsilon(Sino, epps, randoms_correction, rand, E, TOF, nBins, CT);
		}
		if (DEBUG && (MethodList.MRAMLA || MethodList.MBSREM)) {
			mexPrintf("w_vec.epsilon_mramla = %f\n", w_vec.epsilon_mramla);
			mexEvalString("pause(.0001);");
		}
		if (DEBUG) {
			mexPrintf("osem_bool = %u\n", osem_bool);
			mexEvalString("pause(.0001);");
		}
		// Set kernel parameters
		if (osem_bool) {
			kernel.setArg(kernelInd_OSEMTIter++, w_vec.epsilon_mramla);
			kernel.setArg(kernelInd_OSEMTIter++, d_TOFCenter);
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
			if (CT) {
				kernel.setArg(kernelInd_OSEMTIter++, subsets);
				kernel.setArg(kernelInd_OSEMTIter++, d_angles);
				kernel.setArg(kernelInd_OSEMTIter++, w_vec.size_y);
				kernel.setArg(kernelInd_OSEMTIter++, w_vec.dPitch);
				kernel.setArg(kernelInd_OSEMTIter++, w_vec.nProjections);
			}
		}
		if (mlem_bool) {
			kernel_ml.setArg(kernelInd_MLEMT++, w_vec.epsilon_mramla);
			kernel_ml.setArg(kernelInd_MLEMT++, d_TOFCenter);
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
			if (CT) {
				kernel_ml.setArg(kernelInd_MLEMT++, 0U);
				kernel_ml.setArg(kernelInd_MLEMT++, d_angles);
				kernel_ml.setArg(kernelInd_MLEMT++, w_vec.size_y);
				kernel_ml.setArg(kernelInd_MLEMT++, w_vec.dPitch);
				kernel_ml.setArg(kernelInd_MLEMT++, w_vec.nProjections);
			}
		}


		status = af_queue.finish();
		if (status != CL_SUCCESS) {
			getErrorString(status);
			mexPrintf("Queue finish failed\n");
			mexEvalString("pause(.0001);");
			return;
		}
		else if (DEBUG) {
			//getErrorString(status);
			mexPrintf("Queue finish succeeded\n");
			mexEvalString("pause(.0001);");
		}

		// Loop through each iteration
		for (uint32_t iter = iter0; iter < Niter; iter++) {


			cl_ulong st = 0ULL;
			// Compute any of the other algorithms, if applicable
			if (osem_bool) {


				// Loop through the subsets
				for (uint32_t osa_iter = osa_iter0; osa_iter < subsets; osa_iter++) {

					if (osa_iter > osa_iter0 && TOF && !loadTOF) {
						d_Sino[0] = cl::Buffer(af_context, CL_MEM_READ_ONLY, sizeof(float) * length[osa_iter] * nBins, NULL, &status);
#ifdef MX_HAS_INTERLEAVED_COMPLEX
						float* apu = (float*)mxGetSingles(mxGetCell(Sin, tt));
#else
						float* apu = (float*)mxGetData(mxGetCell(Sin, tt));
#endif
						for (int64_t to = 0LL; to < nBins; to++)
							status = af_queue.enqueueWriteBuffer(d_Sino[0], CL_FALSE, sizeof(float) * length[osa_iter] * to, sizeof(float) * length[osa_iter], &apu[pituus[osa_iter] + koko * to]);
						if (status != CL_SUCCESS) {
							getErrorString(status);
							return;
						}
					}

					if (compute_norm_matrix == 1u) {
						if (atomic_64bit) {
							Summ[0] = constant(0LL, im_dim, 1, s64);
						}
						else if (atomic_32bit) {
							Summ[0] = constant(0, im_dim, 1, s32);
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
							else if (atomic_32bit)
								Summ[osa_iter] = constant(0, im_dim, 1, s32);
							else
								Summ[osa_iter] = constant(0.f, im_dim, 1);
							//d_Summ->operator=(*Summ[osa_iter].device<cl_mem>());
							d_Summ = cl::Buffer(*Summ[osa_iter].device<cl_mem>(), true);
						}
						else if (no_norm == 1u) {
							if (atomic_64bit)
								apu_sum = constant(0LL, 1, 1, s64);
							else if (atomic_32bit)
								apu_sum = constant(0, 1, 1, s32);
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


					update_opencl_inputs(vec, vec_opencl, false, im_dim, n_rekos2, n_rekos_mlem, MethodList, atomic_64bit, atomic_32bit, use_psf, w_vec.nMAPOS);

					cl_uint kernelInd_OSEMSubIter = kernelInd_OSEMTIter;

					size_t erotus = length[osa_iter] % local_size;

					if (erotus > 0)
						erotus = (local_size - erotus);

					const size_t global_size = length[osa_iter] + erotus;

					const uint64_t m_size = static_cast<uint64_t>(length[osa_iter]);

					if (DEBUG) {
						mexPrintf("global_size = %u\n", global_size);
						mexPrintf("local_size = %u\n", local_size);
						mexPrintf("size_x = %u\n", size_x);
						mexPrintf("st = %u\n", st);
						mexEvalString("pause(.0001);");
					}

					af::sync();

					// Set kernel arguments
					kernel.setArg(kernelInd_OSEMSubIter++, d_norm[osa_iter]);
					kernel.setArg(kernelInd_OSEMSubIter++, d_scat[osa_iter]);
					kernel.setArg(kernelInd_OSEMSubIter++, d_Summ);
					kernel.setArg(kernelInd_OSEMSubIter++, d_lor[osa_iter]);
					kernel.setArg(kernelInd_OSEMSubIter++, d_xyindex[osa_iter]);
					kernel.setArg(kernelInd_OSEMSubIter++, d_zindex[osa_iter]);
					kernel.setArg(kernelInd_OSEMSubIter++, d_L[osa_iter]);
					if (TOF && !loadTOF)
						kernel.setArg(kernelInd_OSEMSubIter++, d_Sino[0]);
					else
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
						if (compute_norm_matrix == 1u) {
							Summ[0].unlock();
						}
						else {
							if (no_norm == 0u) {
								Summ[osa_iter].unlock();
							}
							else
								apu_sum.unlock();
						}
						if (use_psf)
							vec.im_os_blurred.unlock();
						else
							vec.im_os.unlock();
						vec.rhs_os.unlock();
						break;
					}
					else if (DEBUG) {
						mexPrintf("OS kernel launched successfully\n");
						mexEvalString("pause(.0001);");
					}
					status = af_queue.finish();
					if (status != CL_SUCCESS) {
						getErrorString(status);
						mexPrintf("Queue finish failed after kernel\n");
						mexEvalString("pause(.0001);");
						if (compute_norm_matrix == 1u) {
							Summ[0].unlock();
						}
						else {
							if (no_norm == 0u) {
								Summ[osa_iter].unlock();
							}
							else
								apu_sum.unlock();
						}
						if (use_psf)
							vec.im_os_blurred.unlock();
						else
							vec.im_os.unlock();
						vec.rhs_os.unlock();
						break;
					}
					array* testi;

					// Transfer memory control back to ArrayFire (OS-methods)
					if (compute_norm_matrix == 1u) {
						Summ[0].unlock();
						if (atomic_64bit)
							Summ[0] = Summ[0].as(f32) / TH;
						else if (atomic_32bit)
							Summ[0] = Summ[0].as(f32) / TH32;
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
							else if (atomic_32bit) {
								Summ[osa_iter] = Summ[osa_iter].as(f32) / TH32;
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
					else if (atomic_32bit)
						vec.rhs_os = vec.rhs_os.as(f32) / TH32;
					vec.rhs_os(vec.rhs_os < epps && vec.rhs_os >= 0.f) = epps;
					if (use_psf) {
						vec.rhs_os = computeConvolution(vec.rhs_os, g, Nx, Ny, Nz, w_vec, n_rekos2);
						af::sync();
					}

					if (DEBUG) {
						mexPrintf("Summ = %f\n", af::sum<float>(*testi));
						if (MethodList.MRAMLA)
							mexPrintf("pj3 = %f\n", af::sum<float>(pj3));
						mexPrintf("vec.rhs_os = %f\n", af::sum<float>(vec.rhs_os));
						//af::array apu1 = (vec.im_os / *testi * vec.rhs_os);
						//mexPrintf("apu1 = %f\n", af::sum<float>(apu1));
						//mexPrintf("erotus = %f\n", af::sum<float>(af::abs(*testi - vec.rhs_os)));
						mexEvalString("pause(.0001);");
						//vec.im_os = vec.rhs_os;
					}


					computeOSEstimates(vec, w_vec, MethodList, im_dim, testi, epps, iter, osa_iter, subsets, beta, Nx, Ny, Nz, data, length, d_Sino, break_iter, pj3,
						n_rekos2, pituus, d_lor, d_zindex, d_xyindex, program_mbsrem, af_queue, af_context, Summ, kernel_mramla, d_L, raw, MethodListOpenCL, koko, atomic_64bit, atomic_32bit,
						compute_norm_matrix, OpenCLStruct.kernelNLM, d_sc_ra, kernelInd_MRAMLA, E, d_norm, d_scat, use_psf, g, OpenCLStruct, TOF, loadTOF, Sin, nBins, 
						randoms_correction, local_size, CT);


					if (DEBUG) {
						mexPrintf("vec.im_os = %f\n", af::sum<float>(vec.im_os));
					}

					vec.im_os(vec.im_os < epps) = epps;

					if (verbose) {
						mexPrintf("Sub-iteration %d complete\n", osa_iter + 1u);
						mexEvalString("pause(.0001);");
					}

					st += length[osa_iter];

					status = af_queue.finish();

					if (break_iter)
						break;

				}
				//vec.im_os = vec.rhs_os;

				computeOSEstimatesIter(vec, w_vec, MethodList, im_dim, epps, iter, osa_iter0, subsets, beta, Nx, Ny, Nz, data, n_rekos2, OpenCLStruct, saveIter);
				
				//if (use_psf && w_vec.deconvolution && osem_bool && (saveIter || (!saveIter && iter == Niter - 1))) {
				//	computeDeblur(vec, g, Nx, Ny, Nz, w_vec, MethodList, iter, deblur_iterations, epps, saveIter);
				//}

				if (osem_bool && compute_norm_matrix == 0u)
					no_norm = 1u;

				if (verbose) {
					mexPrintf("Iteration %d complete\n", iter + 1u);
					mexEvalString("pause(.0001);");
				}
				af::deviceGC();
			}
			// Compute MLEM separately
			if (mlem_bool) {

				// Use previously computed normalization factor if available
				if (compute_norm_matrix == 0u) {
					if (osem_bool && no_norm == 1u && iter == 0u) {
						for (uint32_t kk = 0U; kk < subsets; kk++)
							Summ_mlem += Summ[kk];
					}
					if (no_norm_mlem == 1u) {
						if (atomic_64bit)
							apu_sum_mlem = constant(0LL, 1, 1, s64);
						else if (atomic_32bit)
							apu_sum_mlem = constant(0, 1, 1, s32);
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
						else if (atomic_32bit)
							apu_sum_mlem = constant(0, 1, 1, s32);
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
				if (DEBUG) {
					mexPrintf("global_size = %u\n", global_size);
					mexPrintf("size_x = %u\n", size_x);
					mexPrintf("koko = %u\n", koko);
					mexPrintf("st = %u\n", st);
					mexPrintf("no_norm_mlem = %u\n", no_norm_mlem);
					mexEvalString("pause(.0001);");
				}
				// Update the OpenCL inputs for this iteration (image estimates)
				update_opencl_inputs(vec, vec_opencl, true, im_dim, n_rekos, n_rekos_mlem, MethodList, atomic_64bit, atomic_32bit, use_psf, w_vec.nMAPOS);

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
				else if (DEBUG) {
					mexPrintf("OpenCL MLEM kernel executed successfully\n");
					mexEvalString("pause(.0001);");
				}
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
					else if (atomic_32bit)
						Summ_mlem = Summ_mlem.as(f32) / TH32;
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
				else if (atomic_32bit)
					vec.rhs_mlem = vec.rhs_mlem.as(f32) / TH32;
				vec.rhs_mlem(vec.rhs_mlem < epps) = epps;

				if (use_psf) {
					vec.rhs_mlem = computeConvolution(vec.rhs_mlem, g, Nx, Ny, Nz, w_vec, n_rekos_mlem);
					af::sync();
				}

				if (listmode == 2) {
					vec.imEstimates[0] = Summ_mlem;
					break_iter = true;
					break;
				}

				if (DEBUG) {
					mexPrintf("Summ_mlem = %f\n", af::sum<float>(Summ_mlem));
					mexPrintf("vec.rhs_mlem = %f\n", af::sum<float>(vec.rhs_mlem));
					//af::array apu1 = (vec.im_os / *testi * vec.rhs_os);
					//mexPrintf("apu1 = %f\n", af::sum<float>(apu1));
					//mexPrintf("erotus = %f\n", af::sum<float>(af::abs(*testi - vec.rhs_os)));
					mexEvalString("pause(.0001);");
					//vec.im_mlem = vec.rhs_mlem;
					vec.imEstimates[0] = Summ_mlem;
				}

				computeMLEstimates(vec, w_vec, MethodList, im_dim, epps, iter, subsets, beta, Nx, Ny, Nz, data, Summ_mlem, break_iter, OpenCLStruct, saveIter);

				if (no_norm_mlem == 0u)
					no_norm_mlem = 1u;

				if (DEBUG) {
					mexPrintf("vec.im_mlem = %f\n", af::sum<float>(vec.im_mlem));
				}

				//if (use_psf && w_vec.deconvolution && (saveIter || (!saveIter && iter == Niter - 1))) {
				//	computeDeblurMLEM(vec, g, Nx, Ny, Nz, w_vec, MethodList, iter, deblur_iterations, epps, saveIter);
				//}
				if (verbose) {
					mexPrintf("MLEM iteration %d complete\n", iter + 1u);
					mexEvalString("pause(.0001);");
				}
			}
			if (use_psf && w_vec.deconvolution && (saveIter || (!saveIter && iter == Niter - 1))) {
				computeDeblur(vec, g, Nx, Ny, Nz, w_vec, MethodList, iter, deblur_iterations, epps, saveIter);
			}
			if (break_iter)
				break;
		}

		if (MethodList.CUSTOM && (MethodList.MBSREM || MethodList.RBI || MethodList.RBIOSL || MethodList.PKMA))
			af::eval(w_vec.D);

		// Transfer the device data to host MATLAB cell array
		device_to_host_cell(MethodList, vec, oo, cell, w_vec, dimmi, 4);

		if (verbose && listmode != 2) {
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