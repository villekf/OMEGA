/**************************************************************************
* Matrix free computations for OMEGA.
* In this file the CUDA buffers are created, calls to other necessary 
* functions are made and the CUDA kernels are launched. This file 
* contains the code for the matrix-free reconstructions in OMEGA using the
* implementation 2.
* Unlike the non-CUDA/OpenCL versions, this one uses (32-bit) floats and 
* thus can be slightly more inaccurate.
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
#include "AF_cuda_functions.hpp"

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
	uint32_t im_dim = Nxy * Nz;
	const unsigned long long st = 0ULL;
	const unsigned char fp = 0;

	// Distance between rays in multi-ray Siddon
	const float dc_z = cr_pz / static_cast<float>(n_rays3D + 1);

	bool break_iter = false;

	bool loadTOF = true;

	uint32_t t0 = 0u;
	uint32_t iter0 = 0u;
	uint32_t osa_iter0 = 0u;

	// Hard-coded local size
	uint64_t local_size = 32ULL;
	if (projector_type > 1)
		local_size = 128ULL;

	bool atomic_64bit = use_64bit_atomics;
	uint8_t compute_norm_matrix = 1u;
	float mem_portions;
	if (raw == 1u)
		mem_portions = 0.1f;
	else
		mem_portions = 0.2f;
	float image_bytes = static_cast<float>(Nx * Ny * Nz) * static_cast<float>(subsets);

	uint32_t oo = 0u;
	size_t ll = 0ULL;
	float zerof = 0.f;

	// Create output structs
	//matlabArrays ArrayList;
	// Create a struct containing the reconstruction methods used
	RecMethods MethodList;
	RecMethodsOpenCL MethodListOpenCL;

	kernelStruct CUDAStruct;

	// Obtain the reconstruction methods used
	get_rec_methods(options, MethodList);

	const uint8_t listmode = (uint8_t)mxGetScalar(mxGetField(options, 0, "listmode"));
	const bool computeSensImag = (bool)mxGetScalar(mxGetField(options, 0, "compute_sensitivity_image"));
	const bool CT = (bool)mxGetScalar(mxGetField(options, 0, "CT"));

	if (listmode == 2)
		MethodList.MLEM = true;
	OpenCLRecMethods(MethodList, MethodListOpenCL);

	int af_id = af::getDevice();
	CUdevice cuda_id = afcu::getNativeId(af_id);

	CUstream af_cuda_stream = afcu::getStream(cuda_id);


	nvrtcProgram program_os = NULL;
	nvrtcProgram program_ml = NULL;
	nvrtcProgram program_mbsrem = NULL;

	// Create the MATLAB output arrays
	//create_matlab_output(ArrayList, dimmi, MethodList, 4);

	CUDAStruct.af_cuda_stream = &af_cuda_stream;

	// Initial value
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	array x00(Nx * Ny * Nz, (float*)mxGetSingles(mxGetField(options, 0, "x0")), afHost);
#else
	array x00(Nx * Ny * Nz, (float*)mxGetData(mxGetField(options, 0, "x0")), afHost);
#endif

	uint32_t subsetsUsed = subsets;
	// For custom prior
	if (MethodList.CUSTOM) {
		osa_iter0 = (uint32_t)mxGetScalar(mxGetField(options, 0, "osa_iter"));
		iter0 = (uint32_t)mxGetScalar(mxGetField(options, 0, "iter"));
		t0 = (uint32_t)mxGetScalar(mxGetField(options, 0, "tt"));
	}

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

	// Struct containing ArrayFire arrays containing the image estimates and other necessary vectors/matrices
	AF_im_vectors vec;
	// Struct containing beta-values for the MAP-methods
	//Beta beta;
	std::vector<float> beta;
	// Struct containing the necessary variables for the priors
	Weighting w_vec;
	// Struct for TV data
	TVdata data;
	// Struct containing the CUDA image estimates
	CUDA_im_vectors vec_cuda;

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
			vec.im_os(seq(kk * im_dim, (kk + 1) * im_dim - 1)) = x00;
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
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	scat = (float*)mxGetSingles(mxGetCell(mxGetField(options, 0, "ScatterC"), 0));
#else
	scat = (float*)mxGetData(mxGetCell(mxGetField(options, 0, "ScatterC"), 0));
#endif
	if (scatter == 1U) {
		size_scat = mxGetNumberOfElements(mxGetCell(mxGetField(options, 0, "ScatterC"), 0));
	}

	std::vector<array> Summ;
	array Summ_mlem;

	size_t free, total;

	cuMemGetInfo(&free, &total);
	if (((static_cast<float>(free)) > image_bytes * 1.25f && !MethodList.CUSTOM) || (listmode == 1 && computeSensImag))
		compute_norm_matrix = 0u;

	nvrtcResult status1 = NVRTC_SUCCESS;

	status1 = createProgramCUDA(verbose, k_path, fileName, program_os, program_ml, program_mbsrem, atomic_64bit, header_directory,
		projector_type, crystal_size_z, precompute, raw, attenuation_correction, normalization, dec, local_size, n_rays, n_rays3D, MethodList, osem_bool,
		mlem_bool, n_rekos, n_rekos_mlem, w_vec, osa_iter0, cr_pz, dx, use_psf, scatter, randoms_correction, TOF, nBins, listmode, CT);
	if (status1 != NVRTC_SUCCESS) {
		std::cerr << "Error while creating program" << std::endl;
		return;
	}
	else if (DEBUG) {
		mexPrintf("Program created\n");
		mexEvalString("pause(.0001);");
	}

	status1 = NVRTC_SUCCESS;
	CUfunction kernel_os = NULL;
	CUfunction kernel_ml = NULL;
	CUfunction kernel_mbsrem = NULL;
	CUmodule moduleOS, moduleML, moduleMB;

	status1 = createKernelsCUDA(verbose, program_os, program_ml, program_mbsrem, kernel_os, kernel_ml, kernel_mbsrem, CUDAStruct.kernelNLM, CUDAStruct.kernelMed, osem_bool, mlem_bool, MethodList, w_vec,
		precompute, projector_type, n_rays, n_rays3D, moduleOS, moduleML, moduleMB);
	if (status1 != NVRTC_SUCCESS) {
		mexPrintf("Failed to create the kernels\n");
		mexEvalString("pause(.0001);");
		return;
	}
	else if (DEBUG) {
		mexPrintf("Kernels loaded\n");
		mexEvalString("pause(.0001);");
	}

	if (static_cast<double>(total) * 0.75 < static_cast<double>(koko * nBins * sizeof(float)) && TOF)
		loadTOF = false;

	uint32_t TOFsubsets = subsets;
	if (!loadTOF)
		TOFsubsets = 1U;

	// Normalization constant
	// Save the constants if there was enough memory
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

	CUresult status = CUDA_SUCCESS;

	// Create and write buffers
	CUdeviceptr d_x, d_y, d_z, d_atten, d_xcenter, d_ycenter, d_zcenter, d_norm_mlem, d_Sino_mlem, d_sc_ra_mlem, d_V, d_scat_mlem, d_TOFCenter;
	CUdeviceptr* d_Summ, * d_Summ_mlem;
	//CUdeviceptr d_Summ2;
	CUdeviceptr d_lor_mlem, d_L_mlem, d_zindex_mlem;
	CUdeviceptr d_xyindex_mlem, d_pseudos;
	CUdeviceptr d_reko_type, d_reko_type_mlem;
	CUdeviceptr d_angles;

	std::vector<CUdeviceptr> d_lor(subsets);
	std::vector<CUdeviceptr> d_L(subsets);
	std::vector<CUdeviceptr> d_zindex(subsets);
	std::vector<CUdeviceptr> d_xyindex(subsets);
	std::vector<CUdeviceptr> d_Sino(TOFsubsets);
	std::vector<CUdeviceptr> d_sc_ra(subsets);
	std::vector<CUdeviceptr> d_norm(subsets);
	std::vector<CUdeviceptr> d_scat(subsets);
	
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	float* apu = (float*)mxGetSingles(mxGetCell(Sin, 0));
#else
	float* apu = (float*)mxGetData(mxGetCell(Sin, 0));
#endif
	float* angles = nullptr;
	if (CT)
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		angles = (float*)mxGetSingles(mxGetField(options, 0, "angles"));
#else
		angles = (float*)mxGetData(mxGetField(options, 0, "angles"));
#endif

	status = createAndWriteBuffers(d_x, d_y, d_z, d_angles, d_lor, d_L, d_zindex, d_xyindex, d_Sino, d_sc_ra, size_x, size_z, TotSinos, size_atten, size_norm, size_scat, 
		prows, length, x, y, z_det, xy_index, z_index, lor1, L, apu, raw, subsets, pituus, atten, norm, scat, pseudos, V, d_atten, d_norm, d_scat, d_pseudos, d_V,
		d_xcenter, d_ycenter, d_zcenter, x_center, y_center, z_center, size_center_x, size_center_y, size_center_z, size_of_x, size_V, randoms_correction, sc_ra, 
		precompute, d_lor_mlem, d_L_mlem, d_zindex_mlem, d_xyindex_mlem, d_Sino_mlem, d_sc_ra_mlem, d_reko_type, d_reko_type_mlem, osem_bool, mlem_bool, koko, reko_type, 
		reko_type_mlem, n_rekos, n_rekos_mlem, d_norm_mlem, d_scat_mlem, angles, TOF, nBins, loadTOF, d_TOFCenter, TOFCenter, subsetsUsed, osa_iter0, listmode, CT);
	if (status != CUDA_SUCCESS) {
		mexPrintf("Failed to load buffers\n");
		mexEvalString("pause(.0001);");
		return;
	}
	else if (DEBUG) {
		mexPrintf("Buffers loaded\n");
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

	if (mlem_bool) {
		if (atomic_64bit)
			Summ_mlem = constant(0ULL, im_dim, 1, u64);
		else
			Summ_mlem = constant(0.f, im_dim, 1);
	}

	for (uint32_t tt = t0; tt < Nt; tt++) {

		// Compute the prepass phase for MRAMLA, MBSREM, RBI, COSEM, ACOSEM or ECOSEM if applicable
		if (((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI || MethodList.PKMA) && w_vec.MBSREM_prepass ||
			MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0) && (!MethodList.CUSTOM || osa_iter0 == 0u)) {

			uint32_t alku = 0u;

			if ((MethodList.MRAMLA || MethodList.MBSREM) && w_vec.MBSREM_prepass && Nt > 1U)
				w_vec.Amin = constant(0.f, koko * nBins, 1);

			// Run the prepass phase
			MRAMLA_prepass_CUDA(subsets, im_dim, pituus, d_lor, d_zindex, d_xyindex, w_vec, Summ, d_Sino, koko, x00, vec.C_co,
				vec.C_aco, vec.C_osl, alku, d_L, raw, MethodListOpenCL, length, compute_norm_matrix, d_sc_ra, E, det_per_ring, d_pseudos,
				prows, Nx, Ny, Nz, dz, dx, dy, bz, bx, by, bzb, maxxx, maxyy, zmax, NSlices, d_x, d_y, d_z, size_x, TotSinos,
				d_atten, d_norm, d_scat, epps, Nxy, tube_width, crystal_size_z, bmin, bmax, Vmax,
				d_xcenter, d_ycenter, d_zcenter, d_V, dc_z, n_rays, n_rays3D, precompute, projector_type, af_cuda_stream, global_factor, d_reko_type, kernel_mbsrem,
				atomic_64bit, use_psf, g, TOF, loadTOF, Sin, nBins, randoms_correction, sigma_x, d_TOFCenter, d_angles, Nt, CT);


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
		}

		// Loop through each time-step

		uint8_t no_norm = 0u;
		uint8_t no_norm_mlem = 0u;
		if (tt == 0u && compute_norm_matrix == 0u) {
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
			status = cuMemsetD32(d_Sino_mlem, zerof, sizeof(float) * koko);

		// Load the measurement and randoms data from the cell arrays
		if (tt > 0u) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			float* apu2 = (float*)mxGetSingles(mxGetCell(Sin, tt));
#else
			float* apu2 = (float*)mxGetData(mxGetCell(Sin, tt));
#endif
			if (osem_bool) {
				for (uint32_t kk = osa_iter0; kk < subsetsUsed; kk++) {
					if (TOF) {
						if (!loadTOF && kk == 0) {
							status = cuMemAlloc(&d_Sino[kk], sizeof(float) * length[kk] * nBins);
							for (int64_t to = 0LL; to < nBins; to++) {
								CUdeviceptr dst = reinterpret_cast<CUdeviceptr> (reinterpret_cast<char*>(d_Sino[kk]) + sizeof(float) * length[kk] * to);
								char* src = (char*)apu2 + sizeof(float) * (pituus[kk] + koko * to);
								size_t bytes = sizeof(float) * length[kk];
								status = cuMemcpyHtoD(dst, src, bytes);
							}
						}
						else if (loadTOF) {
							for (int64_t to = 0LL; to < nBins; to++) {
								CUdeviceptr dst = reinterpret_cast<CUdeviceptr> (reinterpret_cast<char*>(d_Sino[kk]) + sizeof(float) * length[kk] * to);
								char* src = (char*)apu2 + sizeof(float) * (pituus[kk] + koko * to);
								size_t bytes = sizeof(float) * length[kk];
								status = cuMemcpyHtoD(dst, src, bytes);
							}
						}
					}
					else
						status = cuMemcpyHtoD(d_Sino[kk], &apu2[pituus[kk]], sizeof(float) * length[kk]);
					if (status != CUDA_SUCCESS) {
						std::cerr << getErrorString(status) << std::endl;
					}
					if (randoms_correction) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
						float* ra_apu = (float*)mxGetSingles(mxGetCell(sc_ra, tt));
#else
						float* ra_apu = (float*)mxGetData(mxGetCell(sc_ra, tt));
#endif
						status = cuMemcpyHtoD(d_sc_ra[kk], &ra_apu[pituus[kk]], sizeof(float) * length[kk]);
					}
					else
						status = cuMemcpyHtoD(d_sc_ra[kk], &zerof, sizeof(float));
					if (status != CUDA_SUCCESS) {
						std::cerr << getErrorString(status) << std::endl;
					}
					if (scatter == 1u) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
						scat = (float*)mxGetSingles(mxGetCell(mxGetField(options, 0, "ScatterC"), tt));
#else
						scat = (float*)mxGetData(mxGetCell(mxGetField(options, 0, "ScatterC"), tt));
#endif
						status = cuMemcpyHtoD(d_scat[kk], &scat[pituus[kk]], sizeof(float) * length[kk]);
						if (status != CUDA_SUCCESS) {
							std::cerr << getErrorString(status) << std::endl;
						}
					}
				}
				vec.im_os = constant(0.f, im_dim * n_rekos2, 1);
				for (int kk = 0; kk < n_rekos2; kk++) {
					vec.im_os(seq(kk * im_dim, (kk + 1) * im_dim - 1)) = x00;
				}
			}
			if (mlem_bool) {
				status = cuMemcpyHtoD(d_Sino_mlem, apu2, sizeof(float) * koko * nBins);
				if (status != CUDA_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
				}
				if (randoms_correction) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
					float* ra_apu = (float*)mxGetSingles(mxGetCell(sc_ra, tt));
#else
					float* ra_apu = (float*)mxGetData(mxGetCell(sc_ra, tt));
#endif
					status = cuMemcpyHtoD(d_sc_ra_mlem, ra_apu, sizeof(float) * koko);
				}
				else
					status = cuMemcpyHtoD(d_sc_ra_mlem, &zerof, sizeof(float));
				if (status != CUDA_SUCCESS) {
					std::cerr << getErrorString(status) << std::endl;
				}
				if (scatter == 1u) {
#ifdef MX_HAS_INTERLEAVED_COMPLEX
					scat = (float*)mxGetSingles(mxGetCell(mxGetField(options, 0, "ScatterC"), tt));
#else
					scat = (float*)mxGetData(mxGetCell(mxGetField(options, 0, "ScatterC"), tt));
#endif
					status = cuMemcpyHtoD(d_scat_mlem, scat, sizeof(float) * koko);
					if (status != CUDA_SUCCESS) {
						std::cerr << getErrorString(status) << std::endl;
					}
				}
				vec.im_mlem = constant(0.f, im_dim * n_rekos_mlem, 1);
				for (int kk = 0; kk < n_rekos_mlem; kk++) {
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
			w_vec.epsilon_mramla = MBSREM_epsilon(Sino, epps, randoms_correction, rand, E, TOF, nBins);
		}

		status = cuCtxSynchronize();
		if (status != CUDA_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
			gpuErrchk(cuMemFree(d_z));
			gpuErrchk(cuMemFree(d_x));
			gpuErrchk(cuMemFree(d_y));
			gpuErrchk(cuMemFree(d_angles));
			gpuErrchk(cuMemFree(d_atten));
			gpuErrchk(cuMemFree(d_pseudos));
			gpuErrchk(cuMemFree(d_xcenter));
			gpuErrchk(cuMemFree(d_ycenter));
			gpuErrchk(cuMemFree(d_zcenter));
			gpuErrchk(cuMemFree(d_V));
			if (osem_bool) {
				gpuErrchk(cuModuleUnload(moduleOS));
				gpuErrchk(cuMemFree(d_reko_type));
				for (uint32_t kk = 0u; kk < subsets; kk++) {
					gpuErrchk(cuMemFree(d_lor[kk]));
					gpuErrchk(cuMemFree(d_xyindex[kk]));
					gpuErrchk(cuMemFree(d_zindex[kk]));
					gpuErrchk(cuMemFree(d_L[kk]));
					gpuErrchk(cuMemFree(d_norm[kk]));
					gpuErrchk(cuMemFree(d_scat[kk]));
					//if (randoms_correction)
					gpuErrchk(cuMemFree(d_sc_ra[kk]));
				}
				for (uint32_t kk = 0u; kk < TOFsubsets; kk++)
					gpuErrchk(cuMemFree(d_Sino[kk]));
			}
			if (mlem_bool) {
				gpuErrchk(cuMemFree(d_reko_type_mlem));
				gpuErrchk(cuMemFree(d_lor_mlem));
				gpuErrchk(cuMemFree(d_xyindex_mlem));
				gpuErrchk(cuMemFree(d_zindex_mlem));
				gpuErrchk(cuMemFree(d_L_mlem));
				gpuErrchk(cuMemFree(d_norm_mlem));
				gpuErrchk(cuMemFree(d_scat_mlem));
				gpuErrchk(cuMemFree(d_Sino_mlem));
				gpuErrchk(cuMemFree(d_sc_ra_mlem));
				gpuErrchk(cuModuleUnload(moduleML));
			}
			if ((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI || MethodList.PKMA) && w_vec.MBSREM_prepass ||
				MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0)
				gpuErrchk(cuModuleUnload(moduleMB));
			return;
		}

		//CUdeviceptr d_Summ;
		//status = cuMemAlloc(&d_Summ2, sizeof(float) * im_dim);
		if (status != CUDA_SUCCESS) {
			std::cerr << getErrorString(status) << std::endl;
		}

		// Loop through each iteration
		for (uint32_t iter = iter0; iter < Niter; iter++) {


			uint64_t st = 0ULL;
			// Compute any of the other algorithms, if applicable
			if (osem_bool) {


				// Loop through the subsets
				for (uint32_t osa_iter = osa_iter0; osa_iter < subsets; osa_iter++) {

					array a_Summa;

					if (osa_iter > osa_iter0 && TOF && !loadTOF) {
						status = cuMemAlloc(&d_Sino[0], sizeof(float) * length[osa_iter] * nBins);
#ifdef MX_HAS_INTERLEAVED_COMPLEX
						float* apu = (float*)mxGetSingles(mxGetCell(Sin, tt));
#else
						float* apu = (float*)mxGetData(mxGetCell(Sin, tt));
#endif
						for (int64_t to = 0LL; to < nBins; to++) {
							CUdeviceptr dst = reinterpret_cast<CUdeviceptr> (reinterpret_cast<char*>(d_Sino[0]) + sizeof(float) * length[osa_iter] * to);
							char* src = (char*)apu + sizeof(float) * (pituus[osa_iter] + koko * to);
							size_t bytes = sizeof(float) * length[osa_iter];
							status = cuMemcpyHtoD(dst, src, bytes);
						}
						if (status != CUDA_SUCCESS) {
							std::cerr << getErrorString(status) << std::endl;
							return;
						}
					}
					uint32_t H = osa_iter;
					if (TOF && !loadTOF)
						H = 0U;

					if (compute_norm_matrix == 1u) {
						if (atomic_64bit) {
							Summ[0] = constant(0ULL, im_dim, 1, u64);
						}
						else
							Summ[0] = constant(0.f, im_dim, 1);
						af::eval(Summ[0]);
						af::sync();
						d_Summ = Summ[0].device<CUdeviceptr>();
					}
					else {
						if (no_norm == 0u && tt > 0u) {
							if (atomic_64bit)
								Summ[osa_iter] = constant(0ULL, im_dim, 1, u64);
							else
								Summ[osa_iter] = constant(0.f, im_dim, 1);
							af::eval(Summ[osa_iter]);
							af::sync();
							d_Summ = Summ[osa_iter].device<CUdeviceptr>();
						}
						else if (no_norm == 1u) {
							if (atomic_64bit)
								apu_sum = constant(0ULL, 1, 1, u64);
							else
								apu_sum = constant(0.f, 1, 1, f32);
							af::eval(apu_sum);
							af::sync();
							d_Summ = apu_sum.device<CUdeviceptr>();
						}
						else {
							af::eval(Summ[osa_iter]);
							af::sync();
							d_Summ = Summ[osa_iter].device<CUdeviceptr>();
							//cuMemcpyDtoD(d_Summ2, *Summ[osa_iter].device<CUdeviceptr>(), sizeof(float) * im_dim);
						}
					}

					if (use_psf) {
						vec.im_os_blurred = computeConvolution(vec.im_os, g, Nx, Ny, Nz, w_vec, n_rekos2);
						af::sync();
					}

					update_cuda_inputs(vec, vec_cuda, false, im_dim, n_rekos2, n_rekos_mlem, MethodList, atomic_64bit, use_psf);

					size_t erotus = length[osa_iter] % local_size;

					if (erotus > 0)
						erotus = (local_size - erotus);

					size_t global_size = length[osa_iter] + erotus;

					global_size = global_size / local_size;

					const uint64_t m_size = static_cast<uint64_t>(length[osa_iter]);

					if (DEBUG) {
						mexPrintf("Summ[osa_iter].dims(0) = %u\n", Summ[osa_iter].dims(0));
						mexPrintf("Summ[osa_iter] = %f\n", af::sum<float>(Summ[osa_iter]));
						mexPrintf("global_size = %u\n", global_size);
						mexPrintf("length[osa_iter] = %u\n", length[osa_iter]);
						mexPrintf("erotus = %u\n", erotus);
						mexPrintf("size_x = %u\n", size_x);
						mexPrintf("st = %u\n", st);
						mexEvalString("pause(.0001);");
					}

					// Set kernel arguments
					// Compute the kernel
					if (CT) {
						void* args[] = { (void*)&global_factor, (void*)&epps, &im_dim, (void*)&Nx, (void*)&Ny, (void*)&Nz, (void*)&dz, (void*)&dx, (void*)&dy, (void*)&bz,
							(void*)&bx, (void*)&by, (void*)&bzb, (void*)&maxxx, (void*)&maxyy, (void*)&zmax, (void*)&NSlices, (void*)&size_x, (void*)&TotSinos,
							(void*)&det_per_ring, (void*)&prows, (void*)&Nxy, (void*)&fp, (void*)&sigma_x, (void*)&tube_width, (void*)&crystal_size_z, (void*)&bmin, (void*)&bmax,
							(void*)&Vmax, &w_vec.epsilon_mramla,&d_TOFCenter, &d_atten, &d_pseudos, &d_x, &d_y, &d_z, &d_xcenter, &d_ycenter, &d_zcenter, &d_V, &d_reko_type,
							(void*)&subsets, &d_angles, (void*)&w_vec.size_y, (void*)&w_vec.dPitch, (void*)&w_vec.nProjections, &d_norm[osa_iter], &d_scat[osa_iter],
							reinterpret_cast<void*>(&d_Summ), &d_lor[osa_iter], &d_xyindex[osa_iter], &d_zindex[osa_iter], &d_L[osa_iter], &d_Sino[H], &d_sc_ra[osa_iter],
							reinterpret_cast<void*>(&vec_cuda.d_im_os), reinterpret_cast<void*>(&vec_cuda.d_rhs_os), &no_norm, (void*)&m_size, (void*)&st };
						af::sync();
						status = cuCtxSynchronize();
						if ((status != CUDA_SUCCESS)) {
							std::cerr << getErrorString(status) << std::endl;
						}
						status = cuLaunchKernel(kernel_os, global_size, 1, 1, local_size, 1, 1, 0, af_cuda_stream, &args[0], 0);
					}
					else {
						if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1))) {
							void* args[] = { (void*)&global_factor, (void*)&epps, &im_dim, (void*)&Nx, (void*)&Ny, (void*)&Nz, (void*)&dz, (void*)&dx, (void*)&dy, (void*)&bz,
								(void*)&bx, (void*)&by, (void*)&bzb, (void*)&maxxx, (void*)&maxyy, (void*)&zmax, (void*)&NSlices, (void*)&size_x, (void*)&TotSinos,
								(void*)&det_per_ring, (void*)&prows, (void*)&Nxy, (void*)&fp, (void*)&sigma_x, (void*)&tube_width, (void*)&crystal_size_z, (void*)&bmin, (void*)&bmax, 
								(void*)&Vmax, &w_vec.epsilon_mramla,&d_TOFCenter, &d_atten, &d_pseudos, &d_x, &d_y, &d_z, &d_xcenter, &d_ycenter, &d_zcenter, &d_V, &d_reko_type, 
								&d_norm[osa_iter], &d_scat[osa_iter], reinterpret_cast<void*>(&d_Summ), &d_lor[osa_iter], &d_xyindex[osa_iter], &d_zindex[osa_iter],
								&d_L[osa_iter], &d_Sino[H], &d_sc_ra[osa_iter], reinterpret_cast<void*>(&vec_cuda.d_im_os), reinterpret_cast<void*>(&vec_cuda.d_rhs_os), &no_norm, 
								(void*)&m_size, (void*)&st };
							af::sync();
							status = cuCtxSynchronize();
							if ((status != CUDA_SUCCESS)) {
								std::cerr << getErrorString(status) << std::endl;
							}
							status = cuLaunchKernel(kernel_os, global_size, 1, 1, local_size, 1, 1, 0, af_cuda_stream, &args[0], 0);
						}
						else {
							void* args[] = { (void*)&global_factor, (void*)&epps, (void*)&im_dim, (void*)&Nx, (void*)&Ny, (void*)&Nz, (void*)&dz, (void*)&dx, (void*)&dy, (void*)&bz,
								(void*)&bx, (void*)&by, (void*)&bzb, (void*)&maxxx, (void*)&maxyy, (void*)&zmax, (void*)&NSlices, (void*)&size_x, (void*)&TotSinos,
								(void*)&det_per_ring, (void*)&prows, (void*)&Nxy, (void*)&fp, (void*)&sigma_x, (void*)&dc_z, (void*)&n_rays, &w_vec.epsilon_mramla,&d_TOFCenter, &d_atten, 
								&d_pseudos, &d_x, &d_y, &d_z, &d_reko_type, &d_norm[osa_iter], &d_scat[osa_iter], reinterpret_cast<void*>(&d_Summ), &d_lor[osa_iter],
								&d_xyindex[osa_iter], &d_zindex[osa_iter], &d_L[osa_iter], &d_Sino[H], &d_sc_ra[osa_iter], reinterpret_cast<void*>(&vec_cuda.d_im_os),
								reinterpret_cast<void*>(&vec_cuda.d_rhs_os), &no_norm, (void*)&m_size, (void*)&st };
							af::sync();
							status = cuCtxSynchronize();
							if ((status != CUDA_SUCCESS)) {
								std::cerr << getErrorString(status) << std::endl;
							}
							status = cuLaunchKernel(kernel_os, global_size, 1, 1, local_size, 1, 1, 0, af_cuda_stream, &args[0], 0);
						}
					}


					if ((status != CUDA_SUCCESS)) {
						std::cerr << getErrorString(status) << std::endl;
						gpuErrchk(cuMemFree(d_z));
						gpuErrchk(cuMemFree(d_x));
						gpuErrchk(cuMemFree(d_y));
						if (CT)
							gpuErrchk(cuMemFree(d_angles));
						gpuErrchk(cuMemFree(d_atten));
						gpuErrchk(cuMemFree(d_TOFCenter));
						gpuErrchk(cuMemFree(d_pseudos));
						gpuErrchk(cuMemFree(d_xcenter));
						gpuErrchk(cuMemFree(d_ycenter));
						gpuErrchk(cuMemFree(d_zcenter));
						gpuErrchk(cuMemFree(d_V));
						if (osem_bool) {
							gpuErrchk(cuModuleUnload(moduleOS));
							gpuErrchk(cuMemFree(d_reko_type));
							for (uint32_t kk = 0u; kk < subsets; kk++) {
								gpuErrchk(cuMemFree(d_lor[kk]));
								gpuErrchk(cuMemFree(d_xyindex[kk]));
								gpuErrchk(cuMemFree(d_zindex[kk]));
								gpuErrchk(cuMemFree(d_L[kk]));
								gpuErrchk(cuMemFree(d_norm[kk]));
								gpuErrchk(cuMemFree(d_scat[kk]));
								gpuErrchk(cuMemFree(d_Sino[kk]));
								//if (randoms_correction)
								gpuErrchk(cuMemFree(d_sc_ra[kk]));
							}
						}
						if (mlem_bool) {
							gpuErrchk(cuMemFree(d_reko_type_mlem));
							gpuErrchk(cuMemFree(d_lor_mlem));
							gpuErrchk(cuMemFree(d_xyindex_mlem));
							gpuErrchk(cuMemFree(d_zindex_mlem));
							gpuErrchk(cuMemFree(d_L_mlem));
							gpuErrchk(cuMemFree(d_norm_mlem));
							gpuErrchk(cuMemFree(d_scat_mlem));
							gpuErrchk(cuMemFree(d_Sino_mlem));
							gpuErrchk(cuMemFree(d_sc_ra_mlem));
							gpuErrchk(cuModuleUnload(moduleML));
						}
						mexPrintf("Failed to launch the OS kernel\n");
						return;
					}
					status = cuCtxSynchronize();
					if ((status != CUDA_SUCCESS)) {
						std::cerr << getErrorString(status) << std::endl;
						gpuErrchk(cuMemFree(d_z));
						gpuErrchk(cuMemFree(d_x));
						gpuErrchk(cuMemFree(d_y));
						if (CT)
							gpuErrchk(cuMemFree(d_angles));
						gpuErrchk(cuMemFree(d_atten));
						gpuErrchk(cuMemFree(d_TOFCenter));
						gpuErrchk(cuMemFree(d_pseudos));
						gpuErrchk(cuMemFree(d_xcenter));
						gpuErrchk(cuMemFree(d_ycenter));
						gpuErrchk(cuMemFree(d_zcenter));
						gpuErrchk(cuMemFree(d_V));
						if (osem_bool) {
							gpuErrchk(cuModuleUnload(moduleOS));
							gpuErrchk(cuMemFree(d_reko_type));
							for (uint32_t kk = 0u; kk < subsets; kk++) {
								gpuErrchk(cuMemFree(d_lor[kk]));
								gpuErrchk(cuMemFree(d_xyindex[kk]));
								gpuErrchk(cuMemFree(d_zindex[kk]));
								gpuErrchk(cuMemFree(d_L[kk]));
								gpuErrchk(cuMemFree(d_norm[kk]));
								gpuErrchk(cuMemFree(d_scat[kk]));
								//if (randoms_correction)
								gpuErrchk(cuMemFree(d_sc_ra[kk]));
							}
							for (uint32_t kk = 0u; kk < TOFsubsets; kk++)
								gpuErrchk(cuMemFree(d_Sino[kk]));
						}
						if (mlem_bool) {
							gpuErrchk(cuMemFree(d_reko_type_mlem));
							gpuErrchk(cuMemFree(d_lor_mlem));
							gpuErrchk(cuMemFree(d_xyindex_mlem));
							gpuErrchk(cuMemFree(d_zindex_mlem));
							gpuErrchk(cuMemFree(d_L_mlem));
							gpuErrchk(cuMemFree(d_norm_mlem));
							gpuErrchk(cuMemFree(d_scat_mlem));
							gpuErrchk(cuMemFree(d_Sino_mlem));
							gpuErrchk(cuMemFree(d_sc_ra_mlem));
							gpuErrchk(cuModuleUnload(moduleML));
						}
						if ((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI || MethodList.PKMA) && w_vec.MBSREM_prepass ||
							MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0)
							gpuErrchk(cuModuleUnload(moduleMB));
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
						mexPrintf("Failed to synchronize\n");
						return;
					}
					else if (DEBUG) {
						mexPrintf("OS kernel launched successfully\n");
						mexEvalString("pause(.0001);");
					}
					//float* hOut = new float[im_dim];
					//cuMemcpyDtoH(hOut, d_Summ2, im_dim * sizeof(float));
					//cuCtxSynchronize();
					//Summ[osa_iter] = array(im_dim, hOut, afHost);
					//af::eval(Summ[osa_iter]);
					//af::sync();
					//delete[] hOut;
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


					if (DEBUG) {
						mexPrintf("Summ = %f\n", af::sum<float>(*testi));
						if (MethodList.MRAMLA)
							mexPrintf("pj3 = %f\n", af::sum<float>(pj3));
						mexPrintf("vec.rhs_os = %f\n", af::sum<float>(vec.rhs_os));
						mexPrintf("vec.im_os = %f\n", af::sum<float>(vec.im_os));
						mexEvalString("pause(.0001);");
						//vec.im_os = vec.rhs_os;
					}

					computeOSEstimatesCUDA(vec, w_vec, MethodList, MethodListOpenCL, im_dim, testi, epps, iter, osa_iter, subsets, beta, Nx, Ny, Nz, data, length, 
						d_Sino, break_iter, pj3, n_rekos2, pituus, d_lor, d_zindex, d_xyindex, af_cuda_stream, Summ, kernel_mbsrem, d_L, raw, koko, atomic_64bit,
						compute_norm_matrix, d_sc_ra, E, d_norm, d_scat, det_per_ring, d_pseudos, prows, dz, dx, dy, bz, bx, by, bzb, maxxx, maxyy, zmax, NSlices, 
						d_x, d_y, d_z, size_x, TotSinos, d_atten, Nxy, tube_width, crystal_size_z, bmin, bmax, Vmax, d_xcenter, d_ycenter, d_zcenter, d_V, dc_z, n_rays, 
						n_rays3D, precompute, projector_type, global_factor, d_reko_type, kernel_mbsrem, use_psf, g, CUDAStruct, TOF, loadTOF, Sin, nBins, randoms_correction, 
						sigma_x, d_TOFCenter, d_angles, CT);

					vec.im_os(vec.im_os < epps) = epps;

					if (verbose) {
						mexPrintf("Sub-iteration %d complete\n", osa_iter + 1u);
						mexEvalString("pause(.0001);");
					}

					status = cuCtxSynchronize();
					if ((status != CUDA_SUCCESS)) {
						std::cerr << getErrorString(status) << std::endl;
					}
					st += length[osa_iter];

					if (break_iter)
						break;

				}

				computeOSEstimatesIter(vec, w_vec, MethodList, im_dim, epps, iter, osa_iter0, subsets, beta, Nx, Ny, Nz, data, n_rekos2, CUDAStruct, saveIter);


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
						for (uint32_t kk = 0; kk < subsets; kk++)
							Summ_mlem += Summ[kk];
					}
					if (no_norm_mlem == 1u) {
						if (atomic_64bit)
							apu_sum_mlem = constant(0LL, 1, 1, s64);
						else
							apu_sum_mlem = constant(0.f, 1, 1, f32);
						d_Summ_mlem = apu_sum_mlem.device<CUdeviceptr>();
					}
					else
						d_Summ_mlem = Summ_mlem.device<CUdeviceptr>();
				}
				else {
					if (no_norm_mlem == 1u) {
						if (atomic_64bit)
							apu_sum_mlem = constant(0LL, 1, 1, s64);
						else
							apu_sum_mlem = constant(0.f, 1, 1, f32);
						d_Summ_mlem = apu_sum_mlem.device<CUdeviceptr>();
					}
					else
						d_Summ_mlem = Summ_mlem.device<CUdeviceptr>();
				}
				
				size_t erotus = koko % local_size;

				if (erotus > 0)
					erotus = (local_size - erotus);

				const size_t global_size = koko + erotus;

				const uint64_t m_size = static_cast<uint64_t>(koko);

				if (DEBUG) {
					mexPrintf("global_size = %u\n", global_size);
					mexPrintf("size_x = %u\n", size_x);
					mexPrintf("n_rekos_mlem = %u\n", n_rekos_mlem);
					for (int rr = 0; rr < n_rekos_mlem; rr++)
						mexPrintf("reko_type_mlem = %u\n", reko_type_mlem[rr]);
					mexEvalString("pause(.0001);");
				}

				if (use_psf) {
					vec.im_mlem_blurred = computeConvolution(vec.im_mlem, g, Nx, Ny, Nz, w_vec, n_rekos_mlem);
					af::sync();
				}

				// Update the CUDA inputs for this iteration (image estimates)
				update_cuda_inputs(vec, vec_cuda, true, im_dim, n_rekos, n_rekos_mlem, MethodList, atomic_64bit, use_psf);

				if (CT) {
					void* args[] = { (void*)&global_factor, (void*)&epps, (void*)&im_dim, (void*)&Nx, (void*)&Ny, (void*)&Nz, (void*)&dz, (void*)&dx, (void*)&dy, (void*)&bz,
						(void*)&bx, (void*)&by, (void*)&bzb, (void*)&maxxx, (void*)&maxyy, (void*)&zmax, (void*)&NSlices, (void*)&size_x, (void*)&TotSinos,
						(void*)&det_per_ring, (void*)&prows, (void*)&Nxy, (void*)&fp, (void*)&sigma_x,&d_TOFCenter, (void*)&tube_width, (void*)&crystal_size_z, (void*)&bmin,
						(void*)&bmax, (void*)&Vmax, &w_vec.epsilon_mramla, &d_atten, &d_pseudos, &d_x, &d_y, &d_z, &d_xcenter, &d_ycenter, &d_zcenter, &d_V, &d_reko_type_mlem,
						(void*)&subsets,& d_angles, (void*)&w_vec.size_y, (void*)&w_vec.dPitch, (void*)&w_vec.nProjections, &d_norm_mlem, &d_scat_mlem, 
						reinterpret_cast<void*>(&d_Summ_mlem), &d_lor_mlem, &d_xyindex_mlem, &d_zindex_mlem, &d_L_mlem, &d_Sino_mlem, &d_sc_ra_mlem,
						reinterpret_cast<void*>(&vec_cuda.d_im_mlem), reinterpret_cast<void*>(&vec_cuda.d_rhs_mlem), &no_norm_mlem, (void*)&m_size };
					status = cuLaunchKernel(kernel_ml, global_size, 1, 1, local_size, 1, 1, 0, af_cuda_stream, &args[0], 0);
				}
				else {
					if (projector_type == 2u || projector_type == 3u || (projector_type == 1u && (precompute || (n_rays * n_rays3D) == 1))) {
						void* args[] = { (void*)&global_factor, (void*)&epps, (void*)&im_dim, (void*)&Nx, (void*)&Ny, (void*)&Nz, (void*)&dz, (void*)&dx, (void*)&dy, (void*)&bz, 
							(void*)&bx, (void*)&by, (void*)&bzb, (void*)&maxxx, (void*)&maxyy, (void*)&zmax, (void*)&NSlices, (void*)&size_x, (void*)&TotSinos,
							(void*)&det_per_ring, (void*)&prows, (void*)&Nxy, (void*)&fp, (void*)&sigma_x,&d_TOFCenter, (void*)&tube_width, (void*)&crystal_size_z, (void*)&bmin, 
							(void*)&bmax, (void*)&Vmax, &w_vec.epsilon_mramla, &d_atten, &d_pseudos, &d_x, &d_y, &d_z, &d_xcenter, &d_ycenter, &d_zcenter, &d_V, &d_reko_type_mlem, 
							&d_norm_mlem, &d_scat_mlem, reinterpret_cast<void*>(&d_Summ_mlem), &d_lor_mlem, &d_xyindex_mlem, &d_zindex_mlem, &d_L_mlem, &d_Sino_mlem, &d_sc_ra_mlem, 
							reinterpret_cast<void*>(&vec_cuda.d_im_mlem), reinterpret_cast<void*>(&vec_cuda.d_rhs_mlem), &no_norm_mlem, (void*)&m_size };
						status = cuLaunchKernel(kernel_ml, global_size, 1, 1, local_size, 1, 1, 0, af_cuda_stream, &args[0], 0);
					}
					else {
						void* args[] = { (void*)&global_factor, (void*)&epps, (void*)&im_dim, (void*)&Nx, (void*)&Ny, (void*)&Nz, (void*)&dz, (void*)&dx, (void*)&dy, (void*)&bz, 
							(void*)&bx, (void*)&by, (void*)&bzb, (void*)&maxxx, (void*)&maxyy, (void*)&zmax, (void*)&NSlices, (void*)&size_x, (void*)&TotSinos,
							(void*)&det_per_ring, (void*)&prows, (void*)&Nxy, (void*)&fp, (void*)&sigma_x,&d_TOFCenter, (void*)&dc_z, (void*)&n_rays, &w_vec.epsilon_mramla, &d_atten, 
							&d_pseudos, &d_x, &d_y, &d_z, &d_reko_type_mlem, &d_norm_mlem, &d_scat_mlem, reinterpret_cast<void*>(&d_Summ_mlem), &d_lor_mlem, &d_xyindex_mlem, 
							&d_zindex_mlem, &d_L_mlem, &d_Sino_mlem, &d_sc_ra_mlem, reinterpret_cast<void*>(&vec_cuda.d_im_mlem), reinterpret_cast<void*>(&vec_cuda.d_rhs_mlem), 
							&no_norm_mlem, (void*)&m_size };
						status = cuLaunchKernel(kernel_ml, global_size, 1, 1, local_size, 1, 1, 0, af_cuda_stream, &args[0], 0);
					}
				}

				status = cuCtxSynchronize();
				if ((status != CUDA_SUCCESS)) {
					std::cerr << getErrorString(status) << std::endl;
					gpuErrchk(cuMemFree(d_z));
					gpuErrchk(cuMemFree(d_x));
					gpuErrchk(cuMemFree(d_y));
					if (CT)
						gpuErrchk(cuMemFree(d_angles));
					gpuErrchk(cuMemFree(d_atten));
					gpuErrchk(cuMemFree(d_TOFCenter));
					gpuErrchk(cuMemFree(d_pseudos));
					gpuErrchk(cuMemFree(d_xcenter));
					gpuErrchk(cuMemFree(d_ycenter));
					gpuErrchk(cuMemFree(d_zcenter));
					gpuErrchk(cuMemFree(d_V));
					if (osem_bool) {
						gpuErrchk(cuModuleUnload(moduleOS));
						gpuErrchk(cuMemFree(d_reko_type));
						for (uint32_t kk = 0u; kk < subsets; kk++) {
							gpuErrchk(cuMemFree(d_lor[kk]));
							gpuErrchk(cuMemFree(d_xyindex[kk]));
							gpuErrchk(cuMemFree(d_zindex[kk]));
							gpuErrchk(cuMemFree(d_L[kk]));
							gpuErrchk(cuMemFree(d_norm[kk]));
							gpuErrchk(cuMemFree(d_scat[kk]));
							//if (randoms_correction)
							gpuErrchk(cuMemFree(d_sc_ra[kk]));
						}
						for (uint32_t kk = 0u; kk < TOFsubsets; kk++)
							gpuErrchk(cuMemFree(d_Sino[kk]));
					}
					if (mlem_bool) {
						gpuErrchk(cuMemFree(d_reko_type_mlem));
						gpuErrchk(cuMemFree(d_lor_mlem));
						gpuErrchk(cuMemFree(d_xyindex_mlem));
						gpuErrchk(cuMemFree(d_zindex_mlem));
						gpuErrchk(cuMemFree(d_L_mlem));
						gpuErrchk(cuMemFree(d_norm_mlem));
						gpuErrchk(cuMemFree(d_scat_mlem));
						gpuErrchk(cuMemFree(d_Sino_mlem));
						gpuErrchk(cuMemFree(d_sc_ra_mlem));
						gpuErrchk(cuModuleUnload(moduleML));
					}
					if ((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI || MethodList.PKMA) && w_vec.MBSREM_prepass ||
						MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0)
						gpuErrchk(cuModuleUnload(moduleMB));
					mexPrintf("Failed to launch the MLEM kernel\n");
					break;
				}
				//else if (verbose) {
				//	mexPrintf("MLEM kernel executed successfully\n");
				//	mexEvalString("pause(.0001);");
				//}
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

				if (listmode == 2) {
					vec.imEstimates[0] = Summ_mlem;
					break_iter = true;
					break;
				}

				computeMLEstimates(vec, w_vec, MethodList, im_dim, epps, iter, subsets, beta, Nx, Ny, Nz, data, Summ_mlem, break_iter, CUDAStruct, saveIter);

				if (no_norm_mlem == 0u)
					no_norm_mlem = 1u;

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

		// Transfer the device data to host MATLAB cell array
		device_to_host_cell(MethodList, vec, oo, cell, w_vec, dimmi, 4);

		if (verbose && listmode != 2) {
			mexPrintf("Time step %d complete\n", tt + 1u);
			mexEvalString("pause(.0001);");
		}
		w_vec.MBSREM_prepass = false;
	}

	gpuErrchk(cuCtxSynchronize());

	// Transfer memory control of all variables that weren't used
	//unlock_AF_im_vectors(vec, MethodList, true, mlem_bool, osem_bool, 0u);

	// Clear CUDA buffers
	gpuErrchk(cuMemFree(d_z));
	gpuErrchk(cuMemFree(d_x));
	gpuErrchk(cuMemFree(d_y));
	if (CT)
		gpuErrchk(cuMemFree(d_angles));
	gpuErrchk(cuMemFree(d_atten));
	gpuErrchk(cuMemFree(d_TOFCenter));
	gpuErrchk(cuMemFree(d_pseudos));
	gpuErrchk(cuMemFree(d_xcenter));
	gpuErrchk(cuMemFree(d_ycenter));
	gpuErrchk(cuMemFree(d_zcenter));
	gpuErrchk(cuMemFree(d_V));
	if (osem_bool) {
		gpuErrchk(cuModuleUnload(moduleOS));
		gpuErrchk(cuMemFree(d_reko_type));
		for (uint32_t kk = 0u; kk < subsets; kk++) {
			gpuErrchk(cuMemFree(d_lor[kk]));
			gpuErrchk(cuMemFree(d_xyindex[kk]));
			gpuErrchk(cuMemFree(d_zindex[kk]));
			gpuErrchk(cuMemFree(d_L[kk]));
			gpuErrchk(cuMemFree(d_norm[kk]));
			gpuErrchk(cuMemFree(d_scat[kk]));
			//if (randoms_correction)
			gpuErrchk(cuMemFree(d_sc_ra[kk]));
		}
		for (uint32_t kk = 0u; kk < TOFsubsets; kk++)
			gpuErrchk(cuMemFree(d_Sino[kk]));
	}
	if (mlem_bool) {
		gpuErrchk(cuMemFree(d_reko_type_mlem));
		gpuErrchk(cuMemFree(d_lor_mlem));
		gpuErrchk(cuMemFree(d_xyindex_mlem));
		gpuErrchk(cuMemFree(d_zindex_mlem));
		gpuErrchk(cuMemFree(d_L_mlem));
		gpuErrchk(cuMemFree(d_norm_mlem));
		gpuErrchk(cuMemFree(d_scat_mlem));
		gpuErrchk(cuMemFree(d_Sino_mlem));
		//if (randoms_correction)
		gpuErrchk(cuMemFree(d_sc_ra_mlem));
		gpuErrchk(cuModuleUnload(moduleML));
	}
	if ((MethodList.MRAMLA || MethodList.MBSREM || MethodList.RBIOSL || MethodList.RBI || MethodList.PKMA) && w_vec.MBSREM_prepass ||
		MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || MethodList.OSLCOSEM > 0)
		gpuErrchk(cuModuleUnload(moduleMB));


	return;
}