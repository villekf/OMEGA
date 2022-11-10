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
//#ifdef OPENCL
//#include "AF_opencl_functions.hpp"
//#else
//
//#endif
#include "functions.hpp"

// Use ArrayFire namespace for convenience
using namespace af;

// Main reconstruction function
void reconstruction_AF_matrixfree(const uint16_t* lor1, const float* z_det, const float* x, const float* y, const mxArray* Sin,
	const mxArray* sc_ra, scalarStruct inputScalars, const mxArray* options, const int64_t* pituus, const uint32_t* xy_index,
	const uint16_t* z_index, mxArray* cell, const mwSize* dimmi, const float* atten, const float* norm, const char* k_path, 
	const uint32_t Nt, const uint32_t* pseudos, const uint32_t prows, const uint16_t* L, const char* fileName, const float* x_center,
	const float* y_center, const float* z_center, const char* header_directory, const uint32_t device, uint32_t n_rekos,
	const uint8_t* reko_type, const float* V, const float* gaussian, const size_t size_gauss, const float* TOFCenter) {

	af::setDevice(device);

	// Number of voxels
	inputScalars.Nxy = inputScalars.Nx * inputScalars.Ny;
	inputScalars.im_dim = inputScalars.Nxy * inputScalars.Nz;
	int status = 0;

	// Distance between rays in multi-ray Siddon
	//const float dc_z = cr_pz / static_cast<float>(inputScalars.n_rays3D + 1);

	bool break_iter = false;

	uint32_t t0 = 0u;
	uint32_t iter0 = 0u;

	inputScalars.listmode = (uint8_t)mxGetScalar(getField(options, 0, "listmode"));
	const bool computeSensImag = (bool)mxGetScalar(getField(options, 0, "compute_sensitivity_image"));
	inputScalars.CT = (bool)mxGetScalar(getField(options, 0, "CT"));
	inputScalars.atomic_32bit = (bool)mxGetScalar(getField(options, 0, "use_32bit_atomics"));

	uint8_t compute_norm_matrix = 0u;
	float mem_portions;
	if (inputScalars.raw == 1u)
		mem_portions = 0.1f;
	else
		mem_portions = 0.2f;
	float image_bytes = static_cast<float>(inputScalars.im_dim) * 8.f;

	uint32_t oo = 0u;
	size_t ll = 0ULL;

	// Create a struct containing the reconstruction methods used
	RecMethods MethodList;
	// Same as above, but as cl_ variables
	//RecMethodsOpenCL MethodListOpenCL;

	std::vector<std::vector<float*>> imEstimates;

	//kernelStruct OpenCLStruct;

	// Obtain the reconstruction methods used
	get_rec_methods(options, MethodList);

	//if (inputScalars.listmode == 2)
	//	MethodList.MLEM = true;
	//OpenCLRecMethods(MethodList, MethodListOpenCL);

	// Number of measurements at each subset
	std::vector<int64_t> length(inputScalars.subsets);
	std::vector<int64_t> totLength(1);

	for (uint32_t kk = 0; kk < inputScalars.subsets; kk++)
		length[kk] = pituus[kk + 1u] - pituus[kk];
	totLength[0] = pituus[inputScalars.subsets];

	array D, apu_sum, E, indices, rowInd, values, meanBP, meanFP, outputFP;
	std::vector<array> Summ;

	if ((MethodList.MRAMLA || MethodList.MBSREM) && Nt > 1U)
		E = constant(1.f, inputScalars.koko * inputScalars.nBins, 1);
	//else
	//	E = constant(0.f, 1, 1);

	inputScalars.subsetsUsed = inputScalars.subsets;
	// For custom prior
	if (MethodList.CUSTOM) {
		inputScalars.osa_iter0 = (uint32_t)mxGetScalar(getField(options, 0, "osa_iter"));
		iter0 = (uint32_t)mxGetScalar(getField(options, 0, "iter"));
		t0 = (uint32_t)mxGetScalar(getField(options, 0, "tt"));
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

	uint32_t nIter = 1;
	if (inputScalars.saveIter) {
		imEstimates.resize(n_rekos);
		for (int kk = 0; kk < n_rekos; kk++) {
			imEstimates[kk].resize(inputScalars.Niter + 1);
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			imEstimates[kk][0] = (float*)mxGetSingles(getField(options, 0, "x0"));
#else
			imEstimates[kk][0] = (float*)mxGetData(getField(options, 0, "x0"));
#endif
		}
		//nIter = inputScalars.Niter;
	}
	// Adjust the number of reconstruction methods and create the output vector containing all the estimates
	// E.g. using ECOSEM, the OSEM and COSEM estimates are needed even if they were not selected
	inputScalars.nRekos = n_rekos;
	uint32_t n_rekos2 = n_rekos;
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
		vec.im_os = constant(0.f, inputScalars.im_dim * n_rekos2);
		vec.rhs_os.resize(1);
		//vec.rhs_os = constant(0.f, inputScalars.im_dim * n_rekos2, 1);
		//if (Nt > 1) {
		//	//for (uint32_t tt = 1U; tt < Nt; tt++) {
		//		for (uint32_t kk = 0U; kk < n_rekos2; kk++) {
		//			vec.im_os(seq(kk * inputScalars.im_dim, (kk + 1) * inputScalars.im_dim - 1)) = vec.im_os(seq(0, inputScalars.im_dim - 1), 1);
		//		}
		//	//}
		//}
	inputScalars.nRekos2 = n_rekos2;

	if (DEBUG) {
		mexPrintf("n_rekos2 = %u\n", n_rekos2);
		mexPrintf("n_rekos = %u\n", n_rekos);
		mexPrintf("koko = %u\n", inputScalars.koko);
		mexPrintf("nBins = %u\n", inputScalars.nBins);
		mexEvalString("pause(.0001);");
	}

	// Load the necessary data from the MATLAB input and form the necessary variables
	form_data_variables(vec, beta, w_vec, options, inputScalars, MethodList, data, Nt, iter0, imEstimates);

	uint64_t fullMSize = pituus[inputScalars.subsets];
	if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
		fullMSize *= (static_cast<uint64_t>(w_vec.size_x) * static_cast<uint64_t>(w_vec.size_y));

	ProjectorClass proj;
	if (inputScalars.projector_type != 6)
		status = proj.addProjector(inputScalars, w_vec, MethodList, k_path, header_directory, false);
	if (status != 0)
		return;

	const int64_t mem = proj.getGlobalMem();
	if (((static_cast<float>(mem) * mem_portions) < image_bytes && !MethodList.CUSTOM) || (inputScalars.listmode == 1 && computeSensImag)) {
		if (MethodList.OSEM || MethodList.ECOSEM || MethodList.ROSEM || MethodList.RBI || MethodList.RBIOSL || MethodList.DRAMA || 
			(inputScalars.listmode == 1 && computeSensImag))
			compute_norm_matrix = 1u;
	}
	else if (MethodList.OSEM || MethodList.ECOSEM || MethodList.ROSEM || MethodList.RBI || MethodList.RBIOSL || MethodList.DRAMA)
		compute_norm_matrix = 2u;

	if (compute_norm_matrix == 0)
		proj.no_norm = 1;

	if (MethodList.CUSTOM) {
		if (w_vec.MBSREM_prepass && inputScalars.osa_iter0 == 0U) {
			inputScalars.subsetsUsed = inputScalars.subsets;
		}
		else
			inputScalars.subsetsUsed = inputScalars.osa_iter0 + 1U;
	}
	if (DEBUG) {
		mexPrintf("nProjections = %d\n", w_vec.nProjections);
		mexPrintf("size_y = %u\n", w_vec.size_y);
		mexPrintf("subsetsUsed = %u\n", inputScalars.subsetsUsed);
		mexPrintf("CT = %u\n", inputScalars.CT);
		mexPrintf("listmode = %u\n", inputScalars.listmode);
		mexPrintf("im_dim = %u\n", inputScalars.im_dim);
		mexPrintf("compute_norm_matrix = %u\n", compute_norm_matrix);
		mexEvalString("pause(.0001);");
	}

	float* scat = nullptr;
	inputScalars.scatter = static_cast<uint32_t>((bool)mxGetScalar(getField(options, 0, "scatter")));
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
	scat = (float*)mxGetSingles(mxGetCell(getField(options, 0, "ScatterC"), 0));
#else
	scat = (float*)mxGetData(mxGetCell(getField(options, 0, "ScatterC"), 0));
#endif
	if (inputScalars.scatter == 1U) {
		inputScalars.size_scat = mxGetNumberOfElements(mxGetCell(getField(options, 0, "ScatterC"), 0));
	}

#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
	float* randomsData = (float*)mxGetSingles(mxGetCell(sc_ra, 0));
#else
	float* randomsData = (float*)mxGetData(mxGetCell(sc_ra, 0));
#endif

	// Normalization constant
	// Save the constants if there was enough memory
	// 64-bit atomic operations require unsigned 64-bit integers
	if (compute_norm_matrix == 2u) {
		//if (atomic_64bit && !w_vec.MBSREM_prepass)
		//	Summ.assign(inputScalars.subsets, constant(0LL, im_dim, 1, s64));
		//else if (atomic_32bit && !w_vec.MBSREM_prepass)
		//	Summ.assign(inputScalars.subsets, constant(0, im_dim, 1, s32));
		//else
		//	Summ.assign(inputScalars.subsets, constant(0.f, im_dim, 1));
		for (uint32_t osa_iter = inputScalars.osa_iter0; osa_iter < inputScalars.subsetsUsed; osa_iter++) {
			if (inputScalars.atomic_64bit)
				Summ.emplace_back(constant(0LL, inputScalars.im_dim, 1, s64));
			else if (inputScalars.atomic_32bit)
				Summ.emplace_back(constant(0, inputScalars.im_dim, 1, s32));
			else
				Summ.emplace_back(constant(0.f, inputScalars.im_dim, 1));
			//eval(Summ[osa_iter]);
			//if (DEBUG) {
			//	mexEvalString("pause(5);");
			//}
		}
	}
	//else {
	//	if (inputScalars.atomic_64bit)
	//		Summ.assign(1ULL, constant(0LL, inputScalars.im_dim, 1, s64));
	//	else if (inputScalars.atomic_32bit)
	//		Summ.assign(1ULL, constant(0, inputScalars.im_dim, 1, s32));
	//	else
	//		Summ.assign(1ULL, constant(0.f, inputScalars.im_dim, 1));
	//	//for (uint32_t osa_iter = osa_iter0; osa_iter < inputScalars.subsets; osa_iter++) {
	//	//	if (atomic_64bit)
	//	//		Summ.emplace_back(constant(0LL, im_dim, 1, s64));
	//	//	else if (atomic_32bit)
	//	//		Summ.emplace_back(constant(0, im_dim, 1, s32));
	//	//	else
	//	//		Summ.emplace_back(constant(0.f, im_dim, 1));
	//	//	eval(Summ[osa_iter]);
	//	//	if (DEBUG) {
	//	//		mexEvalString("pause(5);");
	//	//	}
	//	//}
	//}
	//return;

	if (DEBUG) {
		//mexPrintf("Summ[0] = %f\n", af::sum<float>(Summ[0]));
		mexPrintf("w_vec.MBSREM_prepass = %u\n", w_vec.MBSREM_prepass);
		mexPrintf("inputScalars.subsets = %u\n", inputScalars.subsets);
		mexEvalString("pause(.0001);");
	}

	if (static_cast<double>(mem) * 0.75 < static_cast<double>(inputScalars.koko * inputScalars.nBins * sizeof(float)) && inputScalars.TOF)
		inputScalars.loadTOF = false;

	inputScalars.TOFsubsets = inputScalars.subsets;
	if (!inputScalars.loadTOF)
		inputScalars.TOFsubsets = 1U;

	std::vector<af::array> mData(inputScalars.TOFsubsets);
	std::vector<af::array> aRand(inputScalars.TOFsubsets);
	std::vector<af::array> aScat(inputScalars.TOFsubsets);

#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
	float* apu = (float*)mxGetSingles(mxGetCell(Sin, 0));
#else
	float* apu = (float*)mxGetData(mxGetCell(Sin, 0));
#endif
	size_t ind = mxGetNumberOfElements(mxGetCell(Sin, 0));
	//float* angles = nullptr;
	//if (CT)
//#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
//		angles = (float*)mxGetSingles(getField(options, 0, "angles"));
//#else
//		angles = (float*)mxGetData(getField(options, 0, "angles"));
//#endif
	if (inputScalars.projector_type > 3) {
		for (int kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++)
			mData[kk] = array(w_vec.size_x * w_vec.size_y * length[kk], &apu[pituus[kk] * w_vec.size_x * w_vec.size_y], afHost);
		if (inputScalars.randoms_correction) {
			for (int kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++)
				aRand[kk] = array(w_vec.size_x * w_vec.size_y * length[kk], &randomsData[pituus[kk] * w_vec.size_x * w_vec.size_y], afHost);
		}
		if (inputScalars.scatter) {
			for (int kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++)
				aScat[kk] = array(w_vec.size_x * w_vec.size_y * length[kk], &scat[pituus[kk] * w_vec.size_x * w_vec.size_y], afHost);
		}
	}

	if (DEBUG && !inputScalars.CT) {
		mexPrintf("ind = %d\n", ind);
		mexPrintf("TOFsubsets = %d\n", inputScalars.TOFsubsets);
		mexPrintf("TOF = %d\n", inputScalars.TOF);
		mexEvalString("pause(.0001);");
	}

	if (DEBUG && inputScalars.CT) {
		//for (int uu = 0; uu < w_vec.nProjections; uu++)
		//	mexPrintf("angles = %f\n", angles[uu]);
		mexPrintf("nProjections = %d\n", w_vec.nProjections);
		mexPrintf("inputScalars.PITCH = %d\n", inputScalars.PITCH);
		mexPrintf("size_y = %d\n", w_vec.size_y);
		mexPrintf("size_x = %d\n", w_vec.size_x);
		mexPrintf("dPitchX = %f\n", w_vec.dPitchX);
		mexPrintf("dPitchY = %f\n", w_vec.dPitchY);
		mexEvalString("pause(.0001);");
	}

	status = proj.createBuffers(inputScalars, w_vec, x, z_det, xy_index, z_index, lor1, L, pituus, atten, norm, scat, V, x_center, y_center,
		z_center, randomsData, TOFCenter, length, apu, reko_type, MethodList);
	if (status != 0)
		return;

	array g(size_gauss, gaussian, afHost);
	if (inputScalars.use_psf) {
		g = moddims(g, w_vec.g_dim_x * 2u + 1u, w_vec.g_dim_y * 2u + 1u, w_vec.g_dim_z * 2u + 1u);
	}

	if (w_vec.MBSREM_prepass)
		w_vec.MBSREM_prepass = false;

	status = proj.initializeKernel(inputScalars, w_vec);
	if (status != 0)
		return;
	af::sync();
	proj.setDynamicKernelData(inputScalars, w_vec);

	//af_print_mem_info("mem info alussa", -1);

	// Loop through each time-step
	for (uint32_t tt = t0; tt < Nt; tt++) {

		if (tt == 0u && compute_norm_matrix == 2u) {
			if (w_vec.MBSREM_prepass)
				proj.no_norm = 1u;
		}
		else if (tt > 0u) {
			if (compute_norm_matrix == 2u)
				proj.no_norm = 1u;
		}

		// Load precomputed sensitivity image
		if (computeSensImag && inputScalars.listmode == 1) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			float* apu = (float*)mxGetSingles(getField(options, 0, "Summ"));
#else
			float* apu = (float*)mxGetData(getField(options, 0, "Summ"));
#endif
				for (uint32_t osa_iter = 0; osa_iter < inputScalars.subsets; osa_iter++) {
					Summ[osa_iter] = array(inputScalars.im_dim, 1, apu) / static_cast<float>(inputScalars.subsets);
					if (inputScalars.use_psf) {
						Summ[osa_iter] = computeConvolution(Summ[osa_iter], g, inputScalars, w_vec);
						af::sync();
					}
				}
				proj.no_norm = 1u;
		}
		for (uint32_t kk = 0U; kk < inputScalars.nRekos2; kk++) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			vec.im_os(seq(kk * inputScalars.im_dim, (kk + 1) * inputScalars.im_dim - 1), 1, 1) = array(inputScalars.im_dim, (float*)mxGetSingles(getField(options, 0, "x0")), afHost);
#else
			vec.im_os(seq(kk * inputScalars.im_dim, (kk + 1) * inputScalars.im_dim - 1), 1, 1) = array(inputScalars.im_dim, (float*)mxGetData(getField(options, 0, "x0")), afHost);
#endif
		}

		// Load the measurement and randoms data from the cell arrays
		if (tt > 0u) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			float* apu = (float*)mxGetSingles(mxGetCell(Sin, tt));
#else
			float* apu = (float*)mxGetData(mxGetCell(Sin, tt));
#endif
			if (inputScalars.randoms_correction) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
				randomsData = (float*)mxGetSingles(mxGetCell(sc_ra, tt));
#else
				randomsData = (float*)mxGetData(mxGetCell(sc_ra, tt));
#endif
			}
			if (inputScalars.scatter == 1u) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
				scat = (float*)mxGetSingles(mxGetCell(getField(options, 0, "ScatterC"), tt));
#else
				scat = (float*)mxGetData(mxGetCell(getField(options, 0, "ScatterC"), tt));
#endif
			}
			if (inputScalars.projector_type < 4) {
				status = proj.loadDynamicData(inputScalars, length, apu, randomsData, scat, pituus);
				if (status != 0)
					return;
			}
			else {
				for (int kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++)
					mData[kk] = array(w_vec.size_x * w_vec.size_y * length[kk], &apu[pituus[kk] * w_vec.size_x * w_vec.size_y], afHost);
				if (inputScalars.randoms_correction) {
					for (int kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++)
						aRand[kk] = array(w_vec.size_x * w_vec.size_y * length[kk], &randomsData[pituus[kk] * w_vec.size_x * w_vec.size_y], afHost);
				}
				if (inputScalars.scatter) {
					for (int kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++)
						aScat[kk] = array(w_vec.size_x * w_vec.size_y * length[kk], &scat[pituus[kk] * w_vec.size_x * w_vec.size_y], afHost);
				}
			}
		}

		// Compute values needed for MBSREM and MRAMLA
//		if ((MethodList.MBSREM || MethodList.MRAMLA) && Nt > 1U) {
//#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
//			array Sino = array(koko, (float*)mxGetSingles(mxGetCell(Sin, tt)), afHost);
//#else
//			array Sino = array(koko, (float*)mxGetData(mxGetCell(Sin, tt)), afHost);
//#endif
//			array rand;
//			if (inputScalars.randoms_correction)
//#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
//				rand = array(koko, (float*)mxGetSingles(mxGetCell(sc_ra, tt)), afHost);
//#else
//				rand = array(koko, (float*)mxGetData(mxGetCell(sc_ra, tt)), afHost);
//#endif
//			if (w_vec.U == 0.f) {
//				const array Aind = w_vec.Amin > 0.f;
//				if (inputScalars.CT)
//					w_vec.U = max<float>(-af::log(Sino(Aind)) / w_vec.Amin(Aind));
//				else
//					w_vec.U = max<float>(Sino(Aind) / w_vec.Amin(Aind));
//			}
//			w_vec.epsilon_mramla = MBSREM_epsilon(Sino, epps, inputScalars.randoms_correction, rand, E, inputScalars.TOF, inputScalars.nBins, inputScalars.CT);
//		}
		//if (DEBUG && (MethodList.MRAMLA || MethodList.MBSREM)) {
		if (DEBUG) {
			mexPrintf("w_vec.epsilon_mramla = %f\n", w_vec.epsilon_mramla);
			mexPrintf("proj.computeD = %d\n", proj.computeD);
			mexEvalString("pause(.0001);");
		}
		// Set kernel parameters
		//array testi2 = constant(0.f, inputScalars.im_dim);
		//array fLSQR = array(inputScalars.im_dim, (float*)mxGetSingles(getField(options, 0, "x0")), afHost);

		//if (proj.no_norm == 0u && tt > 0u && compute_norm_matrix == 0) {
		//	for (uint32_t osa_iter = inputScalars.osa_iter0; osa_iter < inputScalars.subsets; osa_iter++) {
		//		if (inputScalars.atomic_64bit)
		//			Summ[osa_iter] = constant(0LL, inputScalars.im_dim, 1, s64);
		//		else if (inputScalars.atomic_32bit)
		//			Summ[osa_iter] = constant(0, inputScalars.im_dim, 1, s32);
		//		else
		//			Summ[osa_iter] = constant(0.f, inputScalars.im_dim, 1);
		//	}
		//}
		//if (DEBUG) {
		//	mexPrintf("Sens image allocated\n");
		//}
		if (tt == 0) {
			if (proj.computeD) {
				//initializeRHS(vec, inputScalars);
				if (inputScalars.atomic_64bit)
					D = constant(0LL, inputScalars.im_dim, 1, s64);
				else if (inputScalars.atomic_32bit)
					D = constant(0, inputScalars.im_dim, 1, s32);
				else
					D = constant(0.f, inputScalars.im_dim, 1);
				af::array oneInput = constant(1.f, fullMSize);
				if (inputScalars.projector_type == 6)
					backprojectionSPECT(oneInput, Summ, w_vec, vec, inputScalars, totLength[0], 0, 0, 0, 0, 0);
				else
					status = proj.backwardProjection(vec, inputScalars, w_vec, oneInput, 0, totLength, 0, fullMSize, meanBP, true);
				if (status != 0) {
					vec.rhs_os[0].unlock();
					oneInput.unlock();
					return;
				}
				af::sync();
				vec.rhs_os[0].unlock();
				oneInput.unlock();
				af::sync();
				af::array Dy = vec.rhs_os[0];
				D = Dy.copy();
				if (inputScalars.use_psf) {
					D = computeConvolution(D, g, inputScalars, w_vec);
					//af::sync();
				}
				proj.releaseBuffer(inputScalars);
			}
			if ((MethodList.CPLS || MethodList.CPTV) && w_vec.tauCP == 0.f)
				status = powerMethod(inputScalars, w_vec, length, fullMSize, proj, vec, g, MethodList);
			if (status != 0) {
				return;
			}
		}

		// Loop through each iteration
		for (uint32_t iter = iter0; iter < inputScalars.Niter; iter++) {


			uint32_t curIter = 0;
			uint64_t st = 0ULL;
			if (inputScalars.saveIter)
				curIter = iter;
			// Compute any of the other algorithms, if applicable

				uint32_t uu = 0;

				// Loop through the subsets
				for (uint32_t osa_iter = inputScalars.osa_iter0; osa_iter < inputScalars.subsets; osa_iter++) {
				//for (uint32_t osa_iter = osa_iter0; osa_iter < 1; osa_iter++) {

					uint64_t m_size = length[osa_iter];
					if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
						m_size = static_cast<uint64_t>(w_vec.size_x) * static_cast<uint64_t>(w_vec.size_y) * length[osa_iter];

					if (osa_iter > inputScalars.osa_iter0 && inputScalars.TOF && !inputScalars.loadTOF) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
						float* apu = (float*)mxGetSingles(mxGetCell(Sin, tt));
#else
						float* apu = (float*)mxGetData(mxGetCell(Sin, tt));
#endif
						status = proj.loadTOFData(inputScalars, apu, length[osa_iter], pituus[osa_iter]);
						if (status != 0)
							return;
					}

					if (compute_norm_matrix == 1u) {
						if (inputScalars.atomic_64bit) {
							Summ[0] = constant(0LL, inputScalars.im_dim, 1, s64);
						}
						else if (inputScalars.atomic_32bit) {
							Summ[0] = constant(0, inputScalars.im_dim, 1, s32);
						}
						else
							Summ[0] = constant(0.f, inputScalars.im_dim, 1);
					}
					if (inputScalars.projector_type != 6) {
						if (compute_norm_matrix == 1u)
							proj.transferSensitivityImage(Summ[0]);
						else if (compute_norm_matrix == 2u && proj.no_norm == 0)
							proj.transferSensitivityImage(Summ[osa_iter]);
					}

					status = initializationStep(w_vec, mData[osa_iter], vec, proj, inputScalars, length, m_size, st, MethodList, iter, g, meanBP, Summ);
					if (status != 0)
						return;
					if (DEBUG) {
						mexPrintf("vec.im_os = %f\n", af::sum<float>(vec.im_os));
						mexEvalString("pause(.0001);");
					}
					//else {
					//	if (proj.no_norm == 0u && tt > 0u) {
					//		if (atomic_64bit)
					//			Summ[osa_iter] = constant(0LL, inputScalars.im_dim, 1, s64);
					//		else if (atomic_32bit)
					//			Summ[osa_iter] = constant(0, inputScalars.im_dim, 1, s32);
					//		else
					//			Summ[osa_iter] = constant(0.f, inputScalars.im_dim, 1);
					//		//d_Summ->operator=(*Summ[osa_iter].device<cl_mem>());
					//		if (inputScalars.projector_type < 6) {
					//			d_Summ = cl::Buffer(*Summ[osa_iter].device<cl_mem>(), true);
					//		}
					//	}
					//	else if (proj.no_norm == 1u) {
					//		if (atomic_64bit)
					//			apu_sum = constant(0LL, 1, 1, s64);
					//		else if (atomic_32bit)
					//			apu_sum = constant(0, 1, 1, s32);
					//		else
					//			apu_sum = constant(0.f, 1, 1, f32);
					//		//d_Summ->operator=(*apu_sum.device<cl_mem>());
					//		if (inputScalars.projector_type < 6) {
					//			d_Summ = cl::Buffer(*apu_sum.device<cl_mem>(), true);
					//		}
					//	}
					//	else {
					//		//d_Summ->operator=(*Summ[osa_iter].device<cl_mem>());
					//		if (inputScalars.projector_type < 6)
					//			d_Summ = cl::Buffer(*Summ[osa_iter].device<cl_mem>(), true);
					//	}
					//}

					//if (inputScalars.projector_type == 7) {
					//	indices = af::constant(-1, inputScalars.dec * length[osa_iter], s32);
					//	values = af::constant(0.f, inputScalars.dec * length[osa_iter], f32);
					//	rowInd = af::constant(0, length[osa_iter] + 1, s32);

					//	d_indices = cl::Buffer(*indices.device<cl_mem>(), true);
					//	d_values = cl::Buffer(*values.device<cl_mem>(), true);
					//	d_rowInd = cl::Buffer(*rowInd.device<cl_mem>(), true);
					//}
					//if (inputScalars.projector_type != 6) {
					//	af::sync();
					//	status = proj.update_opencl_inputs(vec, inputScalars);
					//	proj.transferRHS(vec);
					//	if (status != 0)
					//		return;
					//}

					af::sync();

					if (inputScalars.projector_type == 4 || inputScalars.projector_type == 5 || inputScalars.projector_type > 6) {
						//array outputFP = constant(0.f, w_vec.size_x * w_vec.size_y * length[osa_iter]);

						outputFP = constant(0.f, m_size);
						af::sync();
						status = proj.forwardProjection(vec, inputScalars, w_vec, outputFP, osa_iter, length, st, m_size);
						if (status != 0) {
							//if (compute_norm_matrix == 1u) {
							//	Summ[0].unlock();
							//}
							//else if (compute_norm_matrix == 2) {
							//	if (proj.no_norm == 0u) {
							//		Summ[osa_iter].unlock();
							//	}
							//	//else
							//	//	apu_sum.unlock();
							//}
							if (inputScalars.use_psf)
								vec.im_os_blurred.unlock();
							else
								vec.im_os.unlock();
							//vec.rhs_os.unlock();
							outputFP.unlock();
							return;
						}
						af::sync();

						outputFP.unlock();
						if (inputScalars.use_psf)
							vec.im_os_blurred.unlock();
						else
							vec.im_os.unlock();
						//af::array outputFP = constant(0.f, m_size);
						//af::array outputFP = afcl::array(m_size, d_output(), f32, true);
						//clRetainMemObject(d_output());
						//af::array y = afcl::array(m_size, d_Sino[osa_iter](), f32, true);
						//clRetainMemObject(d_Sino[osa_iter]());
						//af::array y = constant(1.f, m_size, 1);
						//vec.im_os(seq(0, (w_vec.size_x)* (w_vec.size_y) - 1)) = exp(-1.f * outputFP(seq(0, (w_vec.size_x) * (w_vec.size_y) - 1)));
						if (DEBUG) {
							mexPrintf("outputFP = %f\n", af::sum<float>(outputFP));
							mexPrintf("min(outputFP) = %f\n", af::min<float>(outputFP));
							mexEvalString("pause(.0001);");
						}
						//vec.im_os(seq((w_vec.size_x)* (w_vec.size_y) * 2, (w_vec.size_x)* (w_vec.size_y) * 3 - 1)) = mData[osa_iter](seq(0, (w_vec.size_x) * (w_vec.size_y) - 1));
						//af::sync();
						computeForwardStep(MethodList, mData[osa_iter], outputFP, m_size, inputScalars, w_vec, aRand[osa_iter], vec);
						////af::eval(outputFP);
						////vec.im_os(seq((w_vec.size_x) * (w_vec.size_y), (w_vec.size_x) * (w_vec.size_y) * 2 - 1)) = outputFP(seq(0, (w_vec.size_x) * (w_vec.size_y) - 1));
						computeIntegralImage(inputScalars, w_vec, length[osa_iter], outputFP, meanBP);
						af::sync();

						status = proj.backwardProjection(vec, inputScalars, w_vec, outputFP, osa_iter, length, st, m_size, meanBP);
						if (status != 0) {
							//getErrorString(status);
							if (compute_norm_matrix == 1u) {
								Summ[0].unlock();
							}
							else if (compute_norm_matrix == 2) {
								if (proj.no_norm == 0u) {
									Summ[osa_iter].unlock();
								}
								//else
								//	apu_sum.unlock();
							}
							//if (inputScalars.use_psf)
							//	vec.im_os_blurred.unlock();
							//else
							//	vec.im_os.unlock();
							vec.rhs_os[0].unlock();
							outputFP.unlock();
							if (inputScalars.meanBP)
								meanBP.unlock();
							//if (inputScalars.projector_type == 5)
							return;
						}
						af::sync();
						outputFP.unlock();
						if (inputScalars.meanBP)
							meanBP.unlock();
						vec.rhs_os[0].unlock();
						//vec.im_os = vec.rhs_os;
						//vec.im_os(seq(0, (w_vec.size_x + 1) * (w_vec.size_y + 1))) = outputFP(seq(0, (w_vec.size_x + 1) * (w_vec.size_y + 1)));
					}
					else if (inputScalars.projector_type == 6) {
						
						array fProj = constant(0.f, inputScalars.size_x, w_vec.size_y, length[osa_iter]);
						forwardProjectionSPECT(fProj, w_vec, vec, inputScalars, length[osa_iter], uu);
						computeForwardStep(MethodList, mData[osa_iter], fProj, m_size, inputScalars, w_vec, aRand[osa_iter], vec); 
						backprojectionSPECT(fProj, Summ, w_vec, vec, inputScalars, length[osa_iter], uu, osa_iter, iter, compute_norm_matrix, iter0);
						uu += length[osa_iter];
					}
					if (inputScalars.projector_type < 4) {
						//af::array outputFP;
						status = proj.backwardProjection(vec, inputScalars, w_vec, outputFP, osa_iter, length, st, m_size, meanBP);

						if (status != 0) {
							//getErrorString(status);
							mexPrintf("Failed to launch the OS kernel\n");
							mexEvalString("pause(.0001);");
							if (compute_norm_matrix == 1u) {
								Summ[0].unlock();
							}
							else if (compute_norm_matrix == 2) {
								if (proj.no_norm == 0u) {
									Summ[osa_iter].unlock();
								}
								//else
								//	apu_sum.unlock();
							}
							//if (inputScalars.use_psf)
							//	vec.im_os_blurred.unlock();
							//else
							//	vec.im_os.unlock();
							vec.rhs_os[0].unlock();
							return;
						}
					}
					af::sync();
					array* testi;
					if (inputScalars.projector_type != 6) {
						//if (status != CL_SUCCESS) {
						//	getErrorString(status);
						//	mexPrintf("Queue finish failed after kernel\n");
						//	mexEvalString("pause(.0001);");
						//	if (compute_norm_matrix == 1u) {
						//		Summ[0].unlock();
						//	}
						//	else if (compute_norm_matrix == 2) {
						//		if (proj.no_norm == 0u) {
						//			Summ[osa_iter].unlock();
						//		}
						//		//else
						//			//apu_sum.unlock();
						//	}
						//	if (inputScalars.use_psf)
						//		vec.im_os_blurred.unlock();
						//	else
						//		vec.im_os.unlock();
						//	vec.rhs_os.unlock();
						//	break;
						//}

						// Transfer memory control back to ArrayFire (OS-methods)
						if (compute_norm_matrix == 1u) {
							Summ[0].unlock();
							if (inputScalars.atomic_64bit)
								Summ[0] = Summ[0].as(f32) / TH;
							else if (inputScalars.atomic_32bit)
								Summ[0] = Summ[0].as(f32) / TH32;
							if (inputScalars.use_psf) {
								Summ[0] = computeConvolution(Summ[0], g, inputScalars, w_vec);
								//af::sync();
							}
							// Prevent division by zero
							Summ[0](Summ[0] < inputScalars.epps) = inputScalars.epps;
							testi = &Summ[0];
							eval(*testi);
							if (DEBUG) {
								mexPrintf("Sens image steps 1 done\n");
								mexEvalString("pause(.0001);");
							}
						}
						else if (compute_norm_matrix == 2) {
							if (proj.no_norm == 0u) {
								Summ[osa_iter].unlock();
								if (inputScalars.atomic_64bit) {
									Summ[osa_iter] = Summ[osa_iter].as(f32) / TH;
								}
								else if (inputScalars.atomic_32bit) {
									Summ[osa_iter] = Summ[osa_iter].as(f32) / TH32;
								}
								if (inputScalars.use_psf) {
									Summ[osa_iter] = computeConvolution(Summ[osa_iter], g, inputScalars, w_vec);
									af::sync();
								}
								Summ[osa_iter](Summ[osa_iter] < inputScalars.epps) = inputScalars.epps;
								if (DEBUG) {
									mexPrintf("Sens image steps 2 done\n");
									mexEvalString("pause(.0001);");
								}
							}
							if (DEBUG) {
								mexPrintf("inputScalars.epps = %f\n", inputScalars.epps);
								mexPrintf("min(Summ) = %f\n", af::min<float>(Summ[osa_iter]));
								mexEvalString("pause(.0001);");
							}
							//else
							//	apu_sum.unlock();
							// Prevent division by zero
							testi = &Summ[osa_iter];
							eval(*testi);
						}
						//if (inputScalars.use_psf)
						//	vec.im_os_blurred.unlock();
						//else
						//	vec.im_os.unlock();
						//vec.rhs_os.unlock();
						af::sync();
						if (DEBUG && inputScalars.atomic_64bit) {
							mexPrintf("min(rhs_os) = %d\n", af::min<int64_t>(vec.rhs_os[0]));
							mexPrintf("inputScalars.atomic_64bit = %d\n", inputScalars.atomic_64bit);
							mexEvalString("pause(.0001);");
						}
						if (inputScalars.atomic_64bit)
							vec.rhs_os[0] = vec.rhs_os[0].as(f32) / TH;
						else if (inputScalars.atomic_32bit)
							vec.rhs_os[0] = vec.rhs_os[0].as(f32) / TH32;
						if (inputScalars.use_psf) {
							//mexPrintf("min(rhs_os) = %f\n", af::min<float>(vec.rhs_os));
							vec.rhs_os[0] = computeConvolution(vec.rhs_os[0], g, inputScalars, w_vec);
							//af::sync();
						}
						if (!MethodList.LSQR && !MethodList.CGLS && !MethodList.CPLS && !MethodList.CPTV)
							vec.rhs_os[0](vec.rhs_os[0] < inputScalars.epps && vec.rhs_os[0] >= 0.f) = inputScalars.epps;

					}
					else {
						testi = &Summ[osa_iter];
						eval(*testi);
					}
					if (DEBUG) {
						if (compute_norm_matrix > 1u) {
							mexPrintf("Summ = %f\n", af::sum<float>(*testi));
							mexPrintf("min(Summ) = %f\n", af::min<float>(*testi));
						}
						mexPrintf("vec.im_os = %f\n", af::sum<float>(vec.im_os));
						mexPrintf("vec.rhs_os = %f\n", af::sum<float>(vec.rhs_os[0]));
						mexPrintf("min(rhs_os) = %f\n", af::min<float>(vec.rhs_os[0]));
						//af::array apu1 = (vec.im_os / *testi * vec.rhs_os);
						//mexPrintf("apu1 = %f\n", af::sum<float>(apu1));
						//mexPrintf("erotus = %f\n", af::sum<float>(af::abs(*testi - vec.rhs_os)));
						mexEvalString("pause(.0001);");
						//vec.im_os = vec.rhs_os;
					}

					////if (iter == 0)
					//if (osa_iter < inputScalars.subsets - 1 || iter < Niter - 1)
					//if (osa_iter < inputScalars.subsets - 1)
					computeOSEstimates(vec, w_vec, MethodList, testi, iter, osa_iter, inputScalars, beta, data, length, break_iter,  
						pituus, Summ, E, g, D, Sin, proj, mData[osa_iter], m_size, uu);


					if (DEBUG) {
						mexPrintf("vec.im_os = %f\n", af::sum<float>(vec.im_os));
						mexEvalString("pause(.0001);");
					}


					if (!MethodList.LSQR && !MethodList.CGLS && !MethodList.CPLS)
						vec.im_os(vec.im_os < inputScalars.epps) = inputScalars.epps;

					if (inputScalars.verbose && inputScalars.subsets > 1) {
						mexPrintf("Sub-iteration %d complete\n", osa_iter + 1u);
						mexEvalString("pause(.0001);");
					}

					//st += length[osa_iter];
					if (inputScalars.projector_type == 6)
						uu += length[osa_iter];

					//status = af_queue.finish();

					if (break_iter)
						break;

					//if (inputScalars.subsets > 1)
					//	vec.rhs_os.clear();
					//af::deviceGC();
				}
				//vec.im_os.lock();
				//vec.im_os.unlock();
				//vec.im_os = vec.rhs_os;
				//vec.im_os = Summ[0];
				//if (MethodList.LSQR) {
				//	//if (iter == 0)
				//	//	testi = vec.im_os;
				//	vec.im_os(seq(0, 0 + inputScalars.im_dim - 1u)) = vec.rhs_os(seq(0, 0 + inputScalars.im_dim - 1u)) - w_vec.betaLSQR * vec.im_os(seq(0, 0 + inputScalars.im_dim - 1u));
				//	w_vec.alphaLSQR = af::norm(vec.im_os(seq(0, 0 + inputScalars.im_dim - 1u)));
				//	vec.im_os(seq(0, 0 + inputScalars.im_dim - 1u)) = vec.im_os(seq(0, 0 + inputScalars.im_dim - 1u)) / w_vec.alphaLSQR;
				//	const float rho_ = sqrt(w_vec.rhoLSQR * w_vec.rhoLSQR + w_vec.betaLSQR * w_vec.betaLSQR);
				//	const float c = w_vec.rhoLSQR / rho_;
				//	const float s = w_vec.betaLSQR / rho_;
				//	w_vec.thetaLSQR = s * w_vec.alphaLSQR;
				//	w_vec.rhoLSQR = -c * w_vec.alphaLSQR;
				//	const float phi_ = c * w_vec.phiLSQR;
				//	w_vec.phiLSQR = s * w_vec.phiLSQR;
				//	//vec.wLSQR.unlock();
				//	//if (DEBUG) {
				//	//	mexPrintf("Step 2\n");
				//	//	mexEvalString("pause(.0001);");
				//	//}
				//	
				//	af::array Df = (phi_ / rho_) * vec.im_os + fLSQR;
				//	fLSQR = Df.copy();
				//	//fLSQR += (phi_ / rho_) * vec.im_os;
				//	//const af::array apuA = vec.im_os.copy();
				//	//eval(apuA);
				//	//const af::array* apuB = &apuA;
				//	//af::sync();
				//	//const af::array apuA = vec.im_os - (w_vec.thetaLSQR / rho_) * vec.wLSQR;
				//	//eval(apuA);
				//	//testi = vec.im_os(seq(0, 0 + inputScalars.im_dim - 1u)) - (w_vec.thetaLSQR / rho_) * testi;
				//	if (iter == inputScalars.Niter - 1)
				//		vec.im_os(seq(0, 0 + inputScalars.im_dim - 1u)) = fLSQR;
				//}
				//fLSQR.lock();
				//fLSQR.unlock();
				af_print_mem_info("mem info", -1);
				//vec.wLSQR = vec.im_os - vec.wLSQR;

				//if (inputScalars.saveIter || iter == inputScalars.Niter - 1)
					computeOSEstimatesIter(vec, w_vec, MethodList, inputScalars, iter, beta, data, imEstimates, proj, g);
					//if (DEBUG) {
					//	mexPrintf("vec.im_os.dims(0) = %d\n", vec.im_os.dims(0));
					//	mexPrintf("vec.im_os.dims(1) = %d\n", vec.im_os.dims(1));
					//	mexEvalString("pause(.0001);");
					//}
				
				//if (use_psf && w_vec.deconvolution && osem_bool && (saveIter || (!saveIter && iter == Niter - 1))) {
				//	computeDeblur(vec, g, Nx, Ny, Nz, w_vec, MethodList, iter, deblur_iterations, epps, saveIter);
				//}

				if (compute_norm_matrix == 2u)
					proj.no_norm = 1u;

				if (inputScalars.verbose) {
					mexPrintf("Iteration %d complete\n", iter + 1u);
					mexEvalString("pause(.0001);");
				}
				//if (inputScalars.subsets == 1)
				//	vec.rhs_os.clear();
				af::deviceGC();
			if (break_iter)
				break;
		}

		if (MethodList.CUSTOM && (MethodList.MBSREM || MethodList.RBI || MethodList.RBIOSL || MethodList.PKMA))
			af::eval(w_vec.D);

		// Transfer the device data to host MATLAB cell array
		device_to_host_cell(MethodList, vec, oo, cell, w_vec, dimmi, 4, imEstimates, inputScalars);

		if (inputScalars.verbose && inputScalars.listmode != 2) {
			mexPrintf("Time step %d complete\n", tt + 1u);
			mexEvalString("pause(.0001);");
		}
		w_vec.MBSREM_prepass = false;
	
	}
	// Transfer memory control of all variables that weren't used
	//unlock_AF_im_vectors(vec, MethodList, true, mlem_bool, osem_bool, 0u);

	//status = af_queue.finish();
	af::sync();
	af::deviceGC();

	return;
}