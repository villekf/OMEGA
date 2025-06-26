/*************************************************************************************************************************************************
* Matrix free computations for OMEGA.
* This is the main reconstruction file for ArrayFire reconstruction. Can utilize OpenCL, CUDA or OneAPI (in the future) as the backend. All the 
* main operations are performed in this file.
* Unlike the non-OpenCL versions, this one uses (32-bit) floats and thus can be slightly more inaccurate.
*
* Copyright (C) 2019-2024 Ville-Veikko Wettenhovi
*
* This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
*************************************************************************************************************************************************/
#pragma once
#define AFTYPE afHost
#include "functions.hpp"
#ifdef MATLAB
#include "mfunctions.h"
#endif
#include "subiterStep.h"
#include "computeForwardProject.h"
#include "iterStep.h"
#include <random>

/// <summary>
/// The main reconstruction function. This function performs all the image reconstruction duties. First it transfers the
/// necessary variables and builds the necessary kernels. After this, it performs any necessary precomputations and then
/// performs the image reconstruction with all the selected features. All the image reconstruction procedures are 
/// performed on the device (GPU). When the computations are done, the results are transfered back to host.
/// </summary>
/// <param name="z_det the z coordinates for the detectors (PET and SPECT) or the directional vectors for the detector panel pixels (CT)"></param>
/// <param name="x the x/y/z coordinates for the detectors (PET and SPECT) or source and detector (CT). z-coordinate applies only for CT"></param>
/// <param name="Sin measurement data (sinogram or projection images)"></param>
/// <param name="sc_ra randoms and/or scatter data. Input nullptr if no randoms/scatter correction."></param>
/// <param name="inputScalars various scalar parameters defining the build and projection parameters"></param>
/// <param name="device selected device number"></param>
/// <param name="pituus cumulative sum of the number of measurements/sinograms/projection images per subset. Has to start at 0"></param>
/// <param name="w_vec specifies algorithm and prior parameters"></param>
/// <param name="MethodList specifies the algorithms and priors used"></param>
/// <param name="header_directory the path to the kernel and header files"></param>
/// <param name="x0 Initial value(s), in case of multi-resolution reconstruction should contain ALL volumes"></param>
/// <param name="cell mxArray where the output is stored (only used in MATLAB/Octave reconstruction), otherwise the output (float*) pointer array for the reconstructions"></param>
/// <param name="FPcell mxArray where the forward projections are stored (only used in MATLAB/Octave reconstruction), otherwise the (float*) pointer array for the forward projections (optional)"></param>
/// <param name="atten attenuation correction images (optional)"></param>
/// <param name="norm normalization correction data (optional)"></param>
/// <param name="extraCorr additional correction vector (optional)"></param>
/// <param name="size_gauss For PSF, the length of the Gaussian kernel vector (optional, required for PSF)"></param>
/// <param name="xy_index subset indices for subset types 3, 6 and 7, x/y dimensions (optional, depends on the selected subset type)"></param>
/// <param name="z_index same as above but for z dimension (optional, depends on the selected subset type)"></param>
/// <param name="residual float* pointer to optional primal-dual gap values, valid only for PDHG type algorithms"></param>
/// <param name="L detector numbers for each ray/LOR. Used by raw data (optional, required for raw data)"></param>
/// <param name="n_rekos Number of selected reconstruction algorithms. NOT USED AT THE MOMENT! (optional)"></param>
/// <returns></returns>
template <typename T, typename F, typename R>
int reconstructionAF(const float* z_det, const float* x, const F* Sin, const R* sc_ra, scalarStruct inputScalars,
	const uint32_t device, const int64_t* pituus, Weighting& w_vec, RecMethods& MethodList, const char* header_directory,
	const float* x0, T* cell, T* FPcell = nullptr, const float* atten = nullptr, const float* norm = nullptr, const float* extraCorr = nullptr, const size_t size_gauss = 0,
	const uint32_t* xy_index = nullptr, const uint16_t* z_index = nullptr, float* residual = nullptr, const uint16_t* L = nullptr, uint32_t n_rekos = 1) {

	if (DEBUG) {
		mexPrint("start\n");
	}
	int curDev = af::getDevice();
	if (curDev != device)
		af::setDevice(device);

	if (inputScalars.verbose >= 3 || DEBUG) {
		mexPrintVar("Compute device set using device number ", af::getDevice());
	}

	// Struct containing ArrayFire arrays containing the image estimates and other necessary vectors/matrices
	AF_im_vectors vec;
	float* apuF, *apuM, *apuU;


	if (MethodList.COSEM || MethodList.ACOSEM || MethodList.OSLCOSEM > 0 || MethodList.ECOSEM) {
		// Complete data
		vec.C_co = af::constant(0.f, inputScalars.im_dim[0], inputScalars.subsets);
	}
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<int> distribution(0, inputScalars.subsets - 1);


	int status = 0;

	bool break_iter = false;

	uint32_t t0 = 0u;
	uint32_t iter0 = 0u;

	uint8_t compute_norm_matrix = 0u;
	float mem_portions;
	if (inputScalars.raw == 1u)
		mem_portions = 0.1f;
	else
		mem_portions = 0.2f;
	float image_bytes = static_cast<float>(inputScalars.im_dim[0] * inputScalars.subsets);

	int64_t oo = 0u;
	size_t ll = 0ULL;

	double totTime = 0., iterTime = 0., cumulativeTime = 0.;
	
	
	if (DEBUG) {
		mexPrintBase("dx = %f\n", inputScalars.dx[0]);
		mexPrintBase("dy = %f\n", inputScalars.dy[0]);
		mexPrintBase("dz = %f\n", inputScalars.dz[0]);
		mexPrintBase("bx = %f\n", inputScalars.bx[0]);
		mexPrintBase("by = %f\n", inputScalars.by[0]);
		mexPrintBase("bz = %f\n", inputScalars.bz[0]);
		mexPrintBase("Nx = %u\n", inputScalars.Nx[0]);
		mexPrintBase("Ny = %u\n", inputScalars.Ny[0]);
		mexPrintBase("Nz = %u\n", inputScalars.Nz[0]);
		mexPrintBase("nProjections = %u\n", inputScalars.nProjections);
		mexPrintBase("nRowsD = %u\n", inputScalars.nRowsD);
		mexPrintBase("nColsD = %u\n", inputScalars.nColsD);
		mexPrintBase("im_dim[0] = %u\n", inputScalars.im_dim[0]);
		mexPrintBase("subsetType = %u\n", inputScalars.subsetType);
		mexPrintBase("useImages = %u\n", inputScalars.useImages);
		mexPrintBase("useMAD = %u\n", inputScalars.useMAD);
		mexPrintBase("flat = %f\n", inputScalars.flat);
		mexPrintBase("dL = %f\n", inputScalars.dL);
		mexPrintBase("epps = %f\n", inputScalars.epps);
		if (inputScalars.projector_type == 4) {
			mexPrintBase("d_Scale4x = %f\n", inputScalars.d_Scale4[0].x);
			mexPrintBase("d_Scale4y = %f\n", inputScalars.d_Scale4[0].y);
			mexPrintBase("d_Scale4z = %f\n", inputScalars.d_Scale4[0].z);
		}
		mexEval();
	}

	if (DEBUG) {
		mexPrintBase("iter0 = %u\n", iter0);
		mexEval();
	}

	uint64_t totDim = std::accumulate(inputScalars.im_dim.begin(), inputScalars.im_dim.end(), (uint64_t)0);

	int64_t nBins = inputScalars.nBins;
	if (inputScalars.listmode && inputScalars.TOF)
		nBins = 1;

	// Number of measurements in each subset
	std::vector<int64_t> length(inputScalars.subsets);
	std::vector<int64_t> lengthFull(inputScalars.subsets);
	std::vector<int64_t> lengthTOF(inputScalars.subsets);
	std::vector<int64_t> totLength(1);

	for (uint32_t kk = 0; kk < inputScalars.subsets; kk++) {
		length[kk] = pituus[kk + 1u] - pituus[kk];
		lengthTOF[kk] = length[kk] * nBins;
		lengthFull[kk] = length[kk];
		if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
			lengthFull[kk] = static_cast<uint64_t>(inputScalars.nRowsD) * static_cast<uint64_t>(inputScalars.nColsD) * lengthFull[kk];
	}
	totLength[0] = pituus[inputScalars.subsets];

	af::array apu_sum, E, indices, rowInd, values, meanBP, meanFP, OSEMapu;

	if ((MethodList.MRAMLA || MethodList.MBSREM) && inputScalars.Nt > 1U)
		E = af::constant(1.f, inputScalars.kokoTOF, 1);

	inputScalars.subsetsUsed = inputScalars.subsets;
	// For custom prior
	// TODO: Re-add support for custom priors
	//if (MethodList.CUSTOM) {
	//	inputScalars.osa_iter0 = (uint32_t)mxGetScalar(getField(options, 0, "osa_iter"));
	//	iter0 = (uint32_t)mxGetScalar(getField(options, 0, "iter"));
	//	t0 = (uint32_t)mxGetScalar(getField(options, 0, "tt"));
	//}

	if (DEBUG) {
		mexPrintBase("subsetsUsed = %u\n", inputScalars.subsetsUsed);
		mexPrintBase("totLength[0] = %d\n", totLength[0]);
		mexEval();
	}

	uint32_t nIter = 1;
	size_t eInd = 0ULL;
	// Create a vector holding the forward projections, if applicable
	std::vector<std::vector<std::vector<float>>> FPEstimates;
	if (inputScalars.storeFP) {
		FPEstimates.resize(inputScalars.Niter);
		for (uint32_t ii = 0; ii < inputScalars.Niter; ii++) {
			FPEstimates[ii].resize(inputScalars.subsets);
			for (uint32_t kk = 0; kk < inputScalars.subsets; kk++)
				FPEstimates[ii][kk].resize(length[kk]);
		}
	}
	// Adjust the number of reconstruction methods and create the output vector containing all the estimates
	// E.g. using ECOSEM, the OSEM and COSEM estimates are needed even if they were not selected
	// Note that in OMEGA v2, the support for multiple algorithms at the same time has been dropped
	// This means that this step is currently unnecessary
	// TODO: Re-add support for multiple algorithms
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
	vec.rhs_os.resize(inputScalars.nMultiVolumes + 1);
	vec.im_os.resize(inputScalars.nMultiVolumes + 1);
	if (inputScalars.use_psf)
		vec.im_os_blurred.resize(inputScalars.nMultiVolumes + 1);
	inputScalars.nRekos2 = n_rekos2;

	if (DEBUG) {
		mexPrintBase("n_rekos2 = %u\n", n_rekos2);
		mexPrintBase("n_rekos = %u\n", n_rekos);
		mexPrintBase("kokoTOF = %u\n", inputScalars.kokoTOF);
		mexPrintBase("nBins = %u\n", nBins);
		mexEval();
	}

	// For the sensitivity image of the whole measurement domain, i.e. no subsets
	uint64_t fullMSize = pituus[inputScalars.subsetsUsed];
	if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
		fullMSize *= (static_cast<uint64_t>(inputScalars.nRowsD) * static_cast<uint64_t>(inputScalars.nColsD));

	// TODO: Add proper support for 64-bit indices
	if (MethodList.NLM || MethodList.MRP || MethodList.ProxTV || MethodList.RDP || (MethodList.TV && !w_vec.data.TV_use_anatomical) || MethodList.ProxTGV || MethodList.ProxRDP || MethodList.GGMRF)
		inputScalars.use64BitIndices = inputScalars.im_dim[0] > 2147483647LL;

	// Create the projector class that handles the projectors and all OpenCL/CUDA/OneAPI related functions
	ProjectorClass proj;
	status = proj.addProjector(inputScalars, w_vec, MethodList, header_directory);
	if (status != 0)
		return -1;

	if (inputScalars.verbose >= 3 || DEBUG)
		mexPrint("Projector class object created");

#ifndef CPU
	// Global memory
	const int64_t mem = proj.getGlobalMem();
	// Is the sensitivity image needed?
	if (((static_cast<float>(mem) * mem_portions) < image_bytes && !MethodList.CUSTOM) && inputScalars.listmode == 0) {
		if ((MethodList.OSEM || MethodList.ECOSEM || MethodList.ROSEM || MethodList.RBI || MethodList.RBIOSL || MethodList.DRAMA || MethodList.OSLOSEM || MethodList.ROSEMMAP || MethodList.SART || MethodList.POCS))
			compute_norm_matrix = 1u;
	}
	else if ((MethodList.OSEM || MethodList.ECOSEM || MethodList.ROSEM || MethodList.RBI || MethodList.RBIOSL || MethodList.DRAMA || MethodList.OSLOSEM || MethodList.ROSEMMAP || MethodList.SART || MethodList.POCS ||
		(inputScalars.listmode == 1 && inputScalars.computeSensImag)))
		compute_norm_matrix = 2u;
#else
	if ((MethodList.OSEM || MethodList.ECOSEM || MethodList.ROSEM || MethodList.RBI || MethodList.RBIOSL || MethodList.DRAMA || MethodList.OSLOSEM || MethodList.ROSEMMAP || MethodList.SART || MethodList.POCS ||
		(inputScalars.listmode == 1 && inputScalars.computeSensImag)))
		compute_norm_matrix = 2u;

#endif

	if (compute_norm_matrix == 0)
		proj.no_norm = 1;
	else if (compute_norm_matrix > 0 && (inputScalars.CT && (MethodList.OSEM || MethodList.ROSEM || MethodList.OSLOSEM || MethodList.ROSEMMAP))) {
		proj.no_norm = 1;
		if (inputScalars.randoms_correction)
			compute_norm_matrix = 1U;
		else
			compute_norm_matrix = 2U;
	}

	if (MethodList.NLM || MethodList.MRP || MethodList.ProxTV || MethodList.RDP || MethodList.TV || MethodList.ProxTGV || MethodList.ProxRDP || MethodList.GGMRF || MethodList.Quad || MethodList.APLS ||
		MethodList.FMH || MethodList.Huber || MethodList.hyperbolic || MethodList.L || MethodList.WeightedMean || MethodList.TGV)
		MethodList.prior = true;

	if (DEBUG) {
		mexPrintBase("nProjections = %d\n", w_vec.nProjections);
		mexPrintBase("nColsD = %u\n", inputScalars.nColsD);
		mexPrintBase("nRowsD = %u\n", inputScalars.nRowsD);
		mexPrintBase("subsetsUsed = %u\n", inputScalars.subsetsUsed);
		mexPrintBase("CT = %u\n", inputScalars.CT);
		mexPrintBase("PET = %u\n", inputScalars.PET);
		mexPrintBase("nLayers = %u\n", inputScalars.nLayers);
		mexPrintBase("listmode = %u\n", inputScalars.listmode);
		mexPrintBase("computeSensImag = %u\n", inputScalars.computeSensImag);
		mexPrintBase("im_dim = %u\n", inputScalars.im_dim[0]);
		mexPrintBase("compute_norm_matrix = %u\n", compute_norm_matrix);
		mexPrintBase("w_vec.computeD = %d\n", w_vec.computeD);
		mexPrintBase("w_vec.computeM = %d\n", w_vec.computeM);
		mexPrintBase("MethodList.prior = %d\n", MethodList.prior);
		mexEval();
	}

	// Sensitivity image
	// Save the images if there was enough memory
	// 64-bit atomic operations require signed 64-bit integers
	// 32-bit atomic operations require signed 32-bit integers
	if (compute_norm_matrix == 2u) {
		vec.Summ.resize(inputScalars.nMultiVolumes + 1);
		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
			vec.Summ[ii].resize(inputScalars.subsetsUsed);
			if (inputScalars.atomic_64bit)
				std::fill(vec.Summ[ii].begin(), vec.Summ[ii].end(), af::constant(0LL, inputScalars.im_dim[ii], 1, s64));
			else if (inputScalars.atomic_32bit)
				std::fill(vec.Summ[ii].begin(), vec.Summ[ii].end(), af::constant(0, inputScalars.im_dim[ii], 1, s32));
			else
				std::fill(vec.Summ[ii].begin(), vec.Summ[ii].end(), af::constant(0.f, inputScalars.im_dim[ii], 1));
			proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]) / 1048576ULL;
		}
	}
	else {
		if (inputScalars.computeSensImag && inputScalars.listmode == 1) {
			vec.Summ.resize(1);
			vec.Summ[0].resize(inputScalars.subsetsUsed);
		}
		else if (compute_norm_matrix == 1u) {
			vec.Summ.resize(inputScalars.nMultiVolumes + 1);
			for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++)
				vec.Summ[ii].resize(1);
		}
	}

	if (DEBUG) {
		mexPrintBase("w_vec.MBSREM_prepass = %u\n", w_vec.MBSREM_prepass);
		mexPrintBase("inputScalars.subsets = %u\n", inputScalars.subsets);
		mexPrintBase("inputScalars.osa_iter0 = %u\n", inputScalars.osa_iter0);
		mexPrintBase("inputScalars.size_of_x = %u\n", inputScalars.size_of_x);
		mexPrintBase("nProjections = %d\n", w_vec.nProjections);
		mexPrintBase("inputScalars.pitch = %d\n", inputScalars.pitch);
		mexPrintBase("inputScalars.randoms_correction = %d\n", inputScalars.randoms_correction);
		mexPrintBase("inputScalars.scatter = %d\n", inputScalars.scatter);
		mexPrintBase("nColsD = %d\n", inputScalars.nColsD);
		mexPrintBase("nRowsD = %d\n", inputScalars.nRowsD);
		mexPrintBase("dPitchX = %f\n", w_vec.dPitchX);
		mexPrintBase("dPitchY = %f\n", w_vec.dPitchY);
		mexEval();
	}

#ifndef CPU
	// Is there sufficient memory to load TOF data to the GPU memory?
	if (inputScalars.loadTOF && inputScalars.listmode == 0) {
		if (static_cast<double>(mem) * 0.7 < static_cast<double>(inputScalars.kokoTOF * sizeof(float)) && inputScalars.TOF)
			inputScalars.loadTOF = false;
		else if (MethodList.CPType && static_cast<double>(mem) * 0.4 < static_cast<double>(inputScalars.kokoTOF * sizeof(float)) && inputScalars.TOF)
			inputScalars.loadTOF = false;
	}
#else
	inputScalars.loadTOF = true;
#endif

	inputScalars.TOFsubsets = inputScalars.subsetsUsed;
	if (!inputScalars.loadTOF || inputScalars.largeDim)
		inputScalars.TOFsubsets = 1U;

	if (DEBUG && !inputScalars.CT) {
		mexPrintBase("TOFsubsets = %d\n", inputScalars.TOFsubsets);
		mexPrintBase("nBins = %d\n", nBins);
		mexPrintBase("TOF = %d\n", inputScalars.TOF);
		mexPrintBase("loadTOF = %d\n", inputScalars.loadTOF);
		mexEval();
	}

	// Measurement data
	std::vector<af::array> mData(inputScalars.TOFsubsets);
	if (!inputScalars.largeDim && !MethodList.FDK) {
		if (inputScalars.TOF && inputScalars.listmode == 0) {
			if (inputScalars.subsetsUsed > 1 && (inputScalars.subsetType < 8 && inputScalars.subsetType > 0))
				for (uint32_t to = 0; to < inputScalars.TOFsubsets; to++) {
					mData[to] = af::constant(0.f, lengthTOF[to]);
					proj.memSize += (sizeof(float) * lengthTOF[to]) / 1048576ULL;
				}
			else
				for (uint32_t to = 0; to < inputScalars.TOFsubsets; to++) {
					mData[to] = af::constant(0.f, inputScalars.nRowsD * inputScalars.nColsD * lengthTOF[to]);
					proj.memSize += (sizeof(float) * inputScalars.nRowsD * inputScalars.nColsD * lengthTOF[to]) / 1048576ULL;
				}
		}
	}
	// Randoms + Scatter
	std::vector<af::array> aRand(inputScalars.TOFsubsets);

	// Create OpenCL buffers, CUDA arrays or OneAPI buffers (in the future)
	status = proj.createBuffers(inputScalars, w_vec, x, z_det, xy_index, z_index, L, pituus, atten, norm, extraCorr, length, MethodList);
	if (status != 0)
		return -1;

	// Special case for large dimensional reconstructions, i.e. cases where the image volume/measurements don't fit into the device memory
	if (inputScalars.largeDim)
		largeDimCreate(inputScalars);

	// Load measurement data into device memory
	if (!inputScalars.largeDim) {
		for (int kk = inputScalars.osa_iter0; kk < inputScalars.TOFsubsets; kk++) {
			if ((inputScalars.subsetsUsed > 1 && inputScalars.subsetType < 8 && inputScalars.subsetType > 0) || inputScalars.listmode > 0)
				if (inputScalars.TOF && inputScalars.listmode == 0) {
					for (int64_t to = 0; to < inputScalars.nBins; to++) {
						mData[kk](af::seq(length[kk] * to, length[kk] * (to + 1) - 1)) = af::array(length[kk], &Sin[pituus[kk] + inputScalars.kokoNonTOF * to], AFTYPE);
					}
				}
				else {
					mData[kk] = af::array(length[kk], &Sin[pituus[kk]], AFTYPE);
					proj.memSize += (sizeof(float) * length[kk]) / 1048576ULL;
				}
			else
				if (inputScalars.TOF && inputScalars.listmode == 0) {
					for (int64_t to = 0; to < inputScalars.nBins; to++)
						mData[kk](af::seq(inputScalars.nRowsD* inputScalars.nColsD* length[kk] * to, inputScalars.nRowsD* inputScalars.nColsD* length[kk] * (to + 1) - 1)) =
						af::array(length[kk] * inputScalars.nRowsD * inputScalars.nColsD, &Sin[pituus[kk] * inputScalars.nRowsD * inputScalars.nColsD + inputScalars.kokoNonTOF * to], AFTYPE);
				}
				else {
					mData[kk] = af::array(inputScalars.nRowsD * inputScalars.nColsD * length[kk], &Sin[pituus[kk] * inputScalars.nRowsD * inputScalars.nColsD], AFTYPE);
					proj.memSize += (sizeof(float) * inputScalars.nRowsD * inputScalars.nColsD * length[kk]) / 1048576ULL;
				}
			if (DEBUG) {
				mexPrintBase("pituus[kk] = %u\n", pituus[kk]);
				mexPrintBase("mData[kk].elements() = %d\n", mData[kk].elements());
				mexPrintBase("inputScalars.nRowsD = %d\n", inputScalars.nRowsD);
				mexPrintBase("inputScalars.nColsD = %d\n", inputScalars.nColsD);
				mexPrintBase("length[kk] = %d\n", length[kk]);
				mexEval();
			}
		}
		if (inputScalars.randoms_correction) {
			for (int kk = inputScalars.osa_iter0; kk < inputScalars.TOFsubsets; kk++)
				if ((inputScalars.subsetsUsed > 1 && (inputScalars.subsetType < 8 && inputScalars.subsetType > 0)) || inputScalars.listmode > 0) {
					aRand[kk] = af::array(length[kk], &sc_ra[pituus[kk]], AFTYPE);
					proj.memSize += (sizeof(float) * length[kk]) / 1048576ULL;
				}
				else {
					aRand[kk] = af::array(inputScalars.nRowsD * inputScalars.nColsD * length[kk], &sc_ra[pituus[kk] * inputScalars.nRowsD * inputScalars.nColsD], AFTYPE);
					proj.memSize += (sizeof(float) * inputScalars.nRowsD * inputScalars.nColsD * length[kk]) / 1048576ULL;
				}
		}
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Measurement/randoms/scatter data loaded");
	}

	// Gaussian kernel for PSF
	af::array g;
	if (inputScalars.use_psf) {
		g = af::array(size_gauss, inputScalars.gaussian, AFTYPE);
		g = af::moddims(g, inputScalars.g_dim_x * 2u + 1u, inputScalars.g_dim_y * 2u + 1u, inputScalars.g_dim_z * 2u + 1u);
	}
	if (DEBUG) {
		mexPrintBase("g_elements = %u\n", g.elements());
		mexEval();
	}

	if (w_vec.MBSREM_prepass)
		w_vec.MBSREM_prepass = false;

#ifndef CPU
	if (inputScalars.projector_type != 6) {
		// Input constant data to the kernels
		status = proj.initializeKernel(inputScalars, w_vec);
		if (status != 0) {
			mexPrint("Failed to initialize kernel");
			return -1;
		}
		af::sync();
		status = proj.setDynamicKernelData(inputScalars, w_vec);
		if (status != 0)
			return -1;
	}
#endif

	initializeProxPriors(MethodList, inputScalars, vec);

	if (inputScalars.verbose >= 3 || DEBUG)
		mexPrint("Constant kernel variables set");
	

	fflush(stdout);
	// Loop through each time-step
	for (uint32_t tt = t0; tt < inputScalars.Nt; tt++) {
		//Compute values needed for MBSREM, MRAMLA and SPS
		if ((MethodList.MBSREM || MethodList.MRAMLA || MethodList.SPS || w_vec.precondTypeIm[6])) {
			if (inputScalars.verbose >= 3)
				mexPrint("Starting computation of E for MBSREM/MRAMLA/SPS");
			af::array Sino = af::array(inputScalars.kokoTOF, &Sin[inputScalars.kokoTOF * tt], AFTYPE);
			af::array rand;
			proj.memSize += (sizeof(float) * inputScalars.kokoTOF) / 1048576ULL;
			if (inputScalars.randoms_correction)
				rand = af::array(inputScalars.kokoNonTOF, &sc_ra[inputScalars.kokoNonTOF * tt], AFTYPE);
			if ((w_vec.precondTypeIm[6] || MethodList.SPS) || (!inputScalars.CT && ((inputScalars.randoms_correction && !af::allTrue<bool>(rand > 0)) || inputScalars.randoms_correction == 0))) {
				int64_t uu = 0;
				if (inputScalars.projector_type == 6)
					E = af::constant(0.f, inputScalars.nRowsD, inputScalars.nColsD, inputScalars.nProjections);
				else
					E = af::constant(0.f, inputScalars.kokoTOF);
				int64_t alku = 0;
				for (uint32_t ll = 0; ll < inputScalars.subsetsUsed; ll++) {
					for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
						af::array oneInput;
						if (inputScalars.use_psf) {
							vec.im_os_blurred[ii] = af::constant(1.f, inputScalars.im_dim[ii]);
						}
						vec.im_os[ii] = af::constant(1.f, inputScalars.im_dim[ii]);
						if (inputScalars.projector_type == 6) {
							oneInput = af::constant(1.f, inputScalars.nRowsD, inputScalars.nColsD, length[ll]);
							forwardProjectionType6(oneInput, w_vec, vec, inputScalars, length[ll], uu, proj, ii, atten);
							uu += length[ll];
						}
						else {
							oneInput = af::constant(1.f, lengthFull[ll] * nBins, 1);
							status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, oneInput, ll, length, g, lengthFull[ll], proj, ii, pituus);
							if (status != 0) {
								return -1;
							}
							af::sync();
						}
						E(af::seq(alku, alku + lengthFull[ll] * nBins - 1)) += oneInput;
						if (DEBUG) {
							mexPrintBase("E = %f\n", af::sum<float>(E));
							mexPrintBase("E.dims(0) = %d\n", E.dims(0));
							mexEval();
						}
					}
					alku += lengthFull[ll] * nBins;
					E.eval();
				}
				E = af::flat(E);
				if (!inputScalars.CT && ((inputScalars.randoms_correction && !af::allTrue<bool>(rand > 0)) || inputScalars.randoms_correction == 0) && (MethodList.MBSREM || MethodList.MRAMLA || MethodList.SPS)) {
					if (inputScalars.verbose >= 3)
						mexPrint("Computing epsilon value for MBSREM/MRAMLA");
					w_vec.epsilon_mramla = MBSREM_epsilon(Sino, E, inputScalars.epps, inputScalars.randoms_correction, rand, inputScalars.TOF, nBins, inputScalars.CT);
					if (inputScalars.verbose >= 3)
						mexPrint("Epsilon value for MBSREM/MRAMLA computed");
				}
				if (w_vec.precondTypeIm[6] || MethodList.SPS) {
					w_vec.dP.resize(inputScalars.nMultiVolumes + 1);
					if (inputScalars.verbose >= 3)
						mexPrint("Starting computation of dP for preconditioner type 6");
					af::array EE;
					if (inputScalars.CT) {
						const af::array apu6 = Sino / inputScalars.flat;
						if (inputScalars.randoms_correction)
							EE = (inputScalars.flat * apu6 - (Sino * rand * apu6) / ((apu6 + rand) * (apu6 + rand)));
						else
							EE = (inputScalars.flat * apu6);
						E *= EE;
					}
					else {
						EE = af::constant(0.f, Sino.elements());
						EE(Sino > 0.f) = E(Sino > 0.f) / Sino(Sino > 0.f);
						E = EE;
					}
					bool noNorm = false;
					if (proj.no_norm == 0) {
						proj.no_norm = 1;
						noNorm = true;
					}
					E.eval();
					computeIntegralImage(inputScalars, w_vec, inputScalars.nProjections, E, meanBP);
					af::sync();
					for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
						int64_t alku = 0;
						for (uint32_t subIter = 0; subIter < inputScalars.subsetsUsed; subIter++) {
							af::array inputM = E(af::seq(alku, alku + lengthFull[subIter] * nBins - 1));
							computeIntegralImage(inputScalars, w_vec, length[subIter], inputM, meanBP);
							af::sync();
							if (inputScalars.projector_type == 6)
								backprojectionType6(inputM, w_vec, vec, inputScalars, length[subIter], uu, proj, subIter, 0, 0, 0, ii, atten);
							else {
								status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, inputM, subIter, length, lengthFull[subIter], meanBP, g, proj, false, ii, pituus);
								if (status != 0) {
									return -1;
								}
								af::sync();
							}
							if (subIter == 0)
								w_vec.dP[ii] = static_cast<float>(inputScalars.subsetsUsed) / vec.rhs_os[ii];
							else {
								w_vec.dP[ii] += static_cast<float>(inputScalars.subsetsUsed) / vec.rhs_os[ii];
							}
							w_vec.dP[ii].eval();
							if (inputScalars.projector_type == 6)
								uu += length[subIter];
							alku += lengthFull[subIter] * nBins;
						}
						w_vec.dP[ii](af::isNaN(w_vec.dP[ii])) = 1.f;
						w_vec.dP[ii](af::isInf(w_vec.dP[ii])) = 1.f;
					}
					if (noNorm)
						proj.no_norm = 0;
					if (DEBUG || inputScalars.verbose >= 3)
						mexPrint("dP computation finished");
					if (DEBUG) {
						const float med = af::median<float>(w_vec.dP[0]);
						for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
							mexPrintBase("normi = %f\n", af::norm(w_vec.dP[ii]));
							mexPrintBase("median = %f\n", af::median<float>(w_vec.dP[ii]));
							mexPrintBase("mediansuhde = %f\n", af::median<float>(w_vec.dP[ii]) / med);
							mexPrintBase("ii = %d\n", ii);
							mexEval();
						}
					}
				}
			}
			if (inputScalars.verbose >= 3 || DEBUG)
				mexPrint("E computation finished");
			proj.memSize -= (sizeof(float) * inputScalars.kokoTOF) / 1048576ULL;
		}

		// Sensitivity image already computed, no need to recompute
		if (tt == 0u && compute_norm_matrix == 2u) {
			if (w_vec.MBSREM_prepass || (inputScalars.listmode > 0 && inputScalars.computeSensImag))
				proj.no_norm = 1u;
		}
		else if (tt > 0u) {
			if (compute_norm_matrix == 2u && !inputScalars.CT) {
				proj.no_norm = 1u;
			}
		}
		uint32_t ee = 0;
		if (inputScalars.largeDim) {
#ifdef MATLAB
			apuF = getSingles(cell, "solu");
#else
			apuF = cell;
#endif
			if (MethodList.PDHG || MethodList.PDHGKL || MethodList.PDHGL1 || MethodList.CV || MethodList.PDDY) {
				apuM = new float[inputScalars.kokoTOF];
				apuU = new float[inputScalars.im_dim[0]];
			}
		}

		// Load the measurement and randoms data for the next time step in case of dynamic data
		if (tt > 0u) {
			if (inputScalars.verbose >= 3 || DEBUG)
				mexPrint("Loading dynamic measurement/randoms/scatter data");
			for (int kk = inputScalars.osa_iter0; kk < inputScalars.TOFsubsets; kk++) {
				if ((inputScalars.subsetsUsed > 1 && inputScalars.subsetType < 8 && inputScalars.subsetType > 0) || inputScalars.listmode > 0) {
					if (inputScalars.TOF && inputScalars.listmode == 0) {
						for (int64_t to = 0; to < inputScalars.nBins; to++) {
							mData[kk](af::seq(length[kk] * to, length[kk] * (to + 1) - 1)) = af::array(length[kk], &Sin[pituus[kk] + inputScalars.kokoNonTOF * to + inputScalars.kokoTOF * tt], AFTYPE);
						}
					}
					else {
						mData[kk] = af::array(length[kk], &Sin[pituus[kk] + inputScalars.kokoNonTOF * tt], AFTYPE);
					}
				}
				else {
					if (inputScalars.TOF && inputScalars.listmode == 0) {
						for (int64_t to = 0; to < inputScalars.nBins; to++)
							mData[kk](af::seq(inputScalars.nRowsD* inputScalars.nColsD* length[kk] * to, inputScalars.nRowsD* inputScalars.nColsD* length[kk] * (to + 1) - 1)) =
							af::array(length[kk] * inputScalars.nRowsD * inputScalars.nColsD, &Sin[pituus[kk] * inputScalars.nRowsD * inputScalars.nColsD + inputScalars.kokoNonTOF * to + inputScalars.kokoTOF * tt], AFTYPE);
					}
					else {
						mData[kk] = af::array(inputScalars.nRowsD * inputScalars.nColsD * length[kk], &Sin[pituus[kk] * inputScalars.nRowsD * inputScalars.nColsD + inputScalars.kokoNonTOF * tt], AFTYPE);
					}
				}
			}
			if (inputScalars.randoms_correction) {
				if ((inputScalars.subsetsUsed > 1 && inputScalars.subsetType < 8 && inputScalars.subsetType > 0) || inputScalars.listmode > 0)
					for (int kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++)
						aRand[kk] = af::array(length[kk], &sc_ra[pituus[kk] + inputScalars.kokoNonTOF * tt], AFTYPE);
				else {
					for (int kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++)
						aRand[kk] = af::array(inputScalars.nRowsD * inputScalars.nColsD * length[kk], &sc_ra[pituus[kk] * inputScalars.nRowsD * inputScalars.nColsD + inputScalars.kokoNonTOF * tt], AFTYPE);
				}
			}
#ifndef CPU
			if (inputScalars.scatter) {
				status = proj.loadDynamicData(inputScalars, length, extraCorr, pituus, tt);
				if (status != 0)
					return -1;
			}
			if (DEBUG || inputScalars.verbose >= 3)
				mexPrint("Dynamic data loaded");
#endif
		}

		if (DEBUG) {
			mexPrintBase("w_vec.epsilon_mramla = %f\n", w_vec.epsilon_mramla);
			mexPrintBase("inputScalars.nMultiVolumes = %d\n", inputScalars.nMultiVolumes);
			mexEval();
		}
		if (tt == 0) {
			// Compute sensitivity image for the whole measurement domain
			if (w_vec.computeD || (inputScalars.listmode > 0 && inputScalars.computeSensImag)) {
				if (DEBUG || inputScalars.verbose >= 3)
					mexPrint("Starting computation of sensitivity image (D)");
				w_vec.D.resize(inputScalars.nMultiVolumes + 1);
				uint8_t apuN = proj.no_norm;
				proj.no_norm = 1;
				int64_t uu = 0;
				for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
					if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
						uint64_t m_size = inputScalars.det_per_ring * inputScalars.det_per_ring * inputScalars.rings * inputScalars.rings;
						if (inputScalars.nLayers > 1)
							m_size *= inputScalars.nLayers;
						if (DEBUG) {
							mexPrintBase("m_size = %u\n", m_size);
							mexPrintBase("inputScalars.det_per_ring = %u\n", inputScalars.det_per_ring);
							mexPrintBase("inputScalars.rings = %u\n", inputScalars.rings);
							mexEval();
						}
						af::array oneInput = af::constant(0.f, 1, 1);
						status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, oneInput, 0, length, m_size, meanBP, g, proj, true, ii, pituus);
						if (status != 0) {
							return -1;
						}
						af::sync();
						vec.rhs_os[ii](vec.rhs_os[ii] < inputScalars.epps) = inputScalars.epps;
						if (w_vec.computeD) {
							w_vec.D[ii] = vec.rhs_os[ii].as(f32);
							w_vec.D[ii].eval();
						}
						else {
							for (int kk = 0; kk < inputScalars.subsets; kk++) {
								vec.Summ[ii][kk] = vec.rhs_os[ii].as(f32) / static_cast<float>(inputScalars.subsets) + inputScalars.epps;
								vec.Summ[ii][kk](vec.Summ[ii][kk] == 0.f) = 1.f;
								vec.Summ[ii][kk].eval();
								if (DEBUG)
									mexPrintBase("sum(vec.Summ[ii][kk]) = %f\n", af::sum<float>(vec.Summ[ii][kk]));
							}
						}
					}
					else {
						for (uint32_t subIter = 0; subIter < inputScalars.subsetsUsed; subIter++) {
							af::array oneInput;
							if (SENSSCALE && inputScalars.CT)
								oneInput = af::flat(mData[subIter]);
							else
								oneInput = af::constant(1.f, lengthFull[subIter] * nBins);
							computeIntegralImage(inputScalars, w_vec, length[subIter], oneInput, meanBP);
							af::sync();
							oneInput.eval();
							if (inputScalars.projector_type == 6)
								backprojectionType6(oneInput, w_vec, vec, inputScalars, length[subIter], uu, proj, subIter, 0, 0, 0, ii, atten);
							else {
								status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, oneInput, subIter, length, lengthFull[subIter], meanBP, g, proj, false, ii, pituus);
								if (status != 0) {
									return -1;
								}
								af::sync();
							}
							if (subIter == 0)
								w_vec.D[ii] = vec.rhs_os[ii].as(f32);
							else {
								w_vec.D[ii] += vec.rhs_os[ii].as(f32);
							}
							w_vec.D[ii].eval();
							if (inputScalars.projector_type == 6)
								uu += length[subIter];
						}
					}
				}
				if (w_vec.computeD) {
					for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
						if (inputScalars.subsetsUsed > 1)
							w_vec.D[ii] /= static_cast<float>(inputScalars.subsetsUsed);
						if (DEBUG)
							mexPrintBase("sum(w_vec.D[ii]) = %f\n", af::sum<float>(w_vec.D[ii]));
						w_vec.D[ii](w_vec.D[ii] == 0.f) = 1.f;
						w_vec.D[ii].eval();
						proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]) / 1048576ULL;
					}
				}
				proj.no_norm = apuN;
				if (inputScalars.listmode > 0 && inputScalars.computeSensImag)
					proj.no_norm = 1;
				if (DEBUG || inputScalars.verbose >= 3)
					mexPrint("D computation finished");
			}
			// Compute the measurement sensitivity image for each subset
			if (w_vec.computeM) {
				if (inputScalars.verbose >= 3)
					mexPrint("Starting computation of measurement sensitivity image (M)");
				int64_t uu = 0;
				for (uint32_t ll = 0; ll < inputScalars.subsetsUsed; ll++) {
					w_vec.M.push_back(af::constant(0.f, lengthFull[ll] * nBins, 1));
					for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
						af::array oneInput;
						if (inputScalars.use_psf) {
							vec.im_os_blurred[ii] = af::constant(1.f, inputScalars.im_dim[ii]);
						}
						vec.im_os[ii] = af::constant(1.f, inputScalars.im_dim[ii]);
						if (inputScalars.projector_type == 6) {
							oneInput = af::constant(1.f, inputScalars.nRowsD, inputScalars.nColsD, length[ll]);
							forwardProjectionType6(oneInput, w_vec, vec, inputScalars, length[ll], uu, proj, ii, atten);
							uu += length[ll];
						}
						else {
							oneInput = af::constant(1.f, lengthFull[ll] * nBins, 1);
							status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, oneInput, ll, length, g, lengthFull[ll], proj, ii, pituus);
							if (status != 0) {
								return - 1;
							}
							af::sync();
						}
						w_vec.M[ll] += oneInput;
						if (DEBUG) {
							mexPrintBase("w_vec.M[ll] = %f\n", af::sum<float>(w_vec.M[ll]));
							mexPrintBase("w_vec.M[ll].dims(0) = %d\n", w_vec.M[ll].dims(0));
							mexEval();
						}
					}
					w_vec.M[ll].eval();
					w_vec.M[ll](w_vec.M[ll] == 0.f) = 1.f;
				}
				if (inputScalars.verbose >= 3)
					mexPrint("M computation finished");
			}
			// Use power method to compute the tau/primal value, if necessary
			if ((MethodList.CPType || MethodList.FISTA || MethodList.FISTAL1) && w_vec.tauCP[0] == 0.f)
				status = powerMethod(inputScalars, w_vec, length, proj, vec, MethodList, g, apuF);
			if (status != 0) {
				return -1;
			}
			if (inputScalars.verbose >= 3 && (MethodList.CPType || MethodList.FISTA || MethodList.FISTAL1)) {
				mexPrintVarf("Primal step size (tau) is ", w_vec.tauCP[0]);
				mexPrintVarf("Dual step size (sigma) is ", w_vec.sigmaCP[0]);
			}
		}
		// Initial values
		int64_t dSum = 0LL;
		if (!inputScalars.largeDim && !MethodList.FDK) {
			for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
				vec.im_os[ii] = af::array(inputScalars.im_dim[ii], &x0[dSum], AFTYPE);
				dSum += inputScalars.im_dim[ii];
				vec.im_os[ii].eval();
				if (DEBUG) {
					mexPrintBase("vec.im_os[ii = %f\n", af::sum<float>(vec.im_os[ii]));
					mexPrintBase("max(vec.im_os[ii) = %f\n", af::max<float>(vec.im_os[ii]));
					mexPrintBase("min(vec.im_os[ii) = %f\n", af::min<float>(vec.im_os[ii]));
					mexPrintBase("ii = %d\n", ii);
					mexPrintBase("inputScalars.im_dim[ii] = %d\n", inputScalars.im_dim[ii]);
					mexEval();
				}
			}
		}

		// Compute the necessary precomputations for the RDP reference image method: https://doi.org/10.1109/TMI.2019.2913889
		if (MethodList.RDP && w_vec.RDPLargeNeighbor && w_vec.RDP_anatomical) {
			int64_t uu = 0;
			w_vec.RDPref = af::constant(0.f, inputScalars.im_dim[0]);
			for (uint32_t ll = 0; ll < inputScalars.subsetsUsed; ll++) {
				af::array oneInput1, oneInput2;
				af::array temp = vec.im_os[0].copy();
				vec.im_os[0] = af::constant(1.f, inputScalars.im_dim[0]);
				if (inputScalars.projector_type == 6) {
					oneInput1 = af::constant(1.f, inputScalars.nRowsD, inputScalars.nColsD, length[ll]);
					forwardProjectionType6(oneInput1, w_vec, vec, inputScalars, length[ll], uu, proj, 0, atten);
					uu += length[ll];
				}
				else {
					oneInput1 = af::constant(1.f, lengthFull[ll] * nBins, 1);
					status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, oneInput1, ll, length, g, lengthFull[ll], proj, 0, pituus);
					if (status != 0) {
						return -1;
					}
					af::sync();
				}
				vec.im_os[0] = af::array(inputScalars.im_dim[0], w_vec.RDP_ref);
				if (inputScalars.projector_type == 6) {
					oneInput2 = af::constant(1.f, inputScalars.nRowsD, inputScalars.nColsD, length[ll]);
					forwardProjectionType6(oneInput2, w_vec, vec, inputScalars, length[ll], uu, proj, 0, atten);
					uu += length[ll];
				}
				else {
					oneInput2 = af::constant(1.f, lengthFull[ll] * nBins, 1);
					status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, oneInput2, ll, length, g, lengthFull[ll], proj, 0, pituus);
					if (status != 0) {
						return -1;
					}
					af::sync();
				}
				vec.im_os[0] = temp.copy();
				af::array apuData = af::constant(0.f, length[ll] * nBins);
				if (inputScalars.subsetsUsed > 1 && inputScalars.subsetType < 8)
					for (int64_t to = 0; to < nBins; to++)
						apuData = af::array(length[ll], &Sin[pituus[ll] + inputScalars.kokoNonTOF * to + inputScalars.kokoTOF * tt], AFTYPE);
				else
					for (int64_t to = 0; to < nBins; to++)
						apuData = af::array(length[ll] * inputScalars.nRowsD * inputScalars.nColsD, &Sin[pituus[ll] * inputScalars.nRowsD * inputScalars.nColsD + inputScalars.kokoNonTOF * to + inputScalars.kokoTOF * tt], AFTYPE);
				oneInput1 = (apuData / (oneInput2 * oneInput2)) * oneInput1;
				af::eval(oneInput1);
				if (inputScalars.projector_type == 6)
					backprojectionType6(oneInput1, w_vec, vec, inputScalars, length[ll], uu, proj, ll, 0, 0, 0, 0, atten);
				else {
					status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, oneInput1, ll, length, lengthFull[ll], meanBP, g, proj, false, 0, pituus);
					if (status != 0) {
						return -1;
					}
					af::sync();
				}
				w_vec.RDPref += af::sqrt(vec.rhs_os[0] + inputScalars.epps);
				af::eval(w_vec.RDPref);
			}
			if (DEBUG) {
				mexPrintBase("w_vec.RDPref = %f\n", af::sum<float>(w_vec.RDPref));
				mexEval();
				if (af::anyTrue<bool>(af::isNaN(w_vec.RDPref)))
					return -1;
				if (af::anyTrue<bool>(af::isInf(w_vec.RDPref)))
					return -1;
			}
		}

		if (inputScalars.stochastic) {
			inputScalars.subsets = 1;
			//inputScalars.subsetsUsed = 1;
		}
		// FDK only computations
		if (MethodList.FDK) {
			inputScalars.Niter = 0;
			// Regular FDK
			if (!inputScalars.largeDim) {
				status = applyMeasPreconditioning(w_vec, inputScalars, mData[0], proj, 0);
				if (status != 0)
					return -1;
				computeIntegralImage(inputScalars, w_vec, length[0], mData[0], meanBP);
				status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, mData[0], 0, length, fullMSize, meanBP, g, proj, false, 0, pituus);
				if (status != 0) {
					return -1;
				}
			}
			else {
				// Large-dimensional case
				bool FDK = false;
				if (inputScalars.CT)
					FDK = true;
				// loop through "subsets"
				// Both the measurements and the image volume are divided into inputScalars.subsets parts
				// Loop through the measurement "subsets"
				for (int kk = 0; kk < inputScalars.subsets; kk++) {
					fullMSize = length[kk] * (static_cast<uint64_t>(inputScalars.nRowsD) * static_cast<uint64_t>(inputScalars.nColsD));
					mData[0] = af::log(inputScalars.flat / af::array(inputScalars.nRowsD * inputScalars.nColsD * length[kk], &Sin[pituus[kk] * inputScalars.nRowsD * inputScalars.nColsD], AFTYPE).as(f32));
					if (DEBUG) {
						mexPrintBase("mData[0] = %f\n", af::sum<float>(mData[0]));
						mexPrintBase("mData[0] = %f\n", af::max<float>(mData[0]));
						mexPrintBase("mData[0] = %f\n", af::min<float>(mData[0]));
						mexPrintBase("kk = %d\n", kk);
						mexPrintBase("length[kk] = %d\n", length[kk]);
						mexPrintBase("pituus[kk] = %d\n", pituus[kk]);
						mexEval();
					}
					status = applyMeasPreconditioning(w_vec, inputScalars, mData[0], proj, 0);
					if (status != 0)
						return -1;
					computeIntegralImage(inputScalars, w_vec, length[kk], mData[0], meanBP);
					// Loop through the image volume for the current measurement "subset"
					for (int ii = 0; ii < inputScalars.subsets; ii++) {
						largeDimFirst(inputScalars, proj, ii);
						if (FDK)
							vec.rhs_os[0] = af::array(inputScalars.lDimStruct.imDim[ii], &apuF[inputScalars.lDimStruct.cumDim[ii]], afHost);
						status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, mData[0], kk, length, fullMSize, meanBP, g, proj, false, 0, pituus, FDK);
						if (status != 0) {
							return -1;
						}
						if (!FDK)
							vec.rhs_os[0] += af::array(inputScalars.lDimStruct.imDim[ii], &apuF[inputScalars.lDimStruct.cumDim[ii]], afHost);
						vec.rhs_os[0].eval();					
						if (DEBUG) {
							mexPrintBase("vec.rhs_os[0] = %f\n", af::sum<float>(vec.rhs_os[0]));
							mexPrintBase("vec.rhs_os[0] = %f\n", af::max<float>(vec.rhs_os[0]));
							mexPrintBase("vec.rhs_os[0] = %f\n", af::min<float>(vec.rhs_os[0]));
							mexPrintBase("ii = %d\n", ii);
							mexPrintBase("inputScalars.lDimStruct.imDim[ii] = %d\n", inputScalars.lDimStruct.imDim[ii]);
							mexPrintBase("inputScalars.lDimStruct.cumDim[ii]] = %d\n", inputScalars.lDimStruct.cumDim[ii]);
							mexEval();
						}
						vec.rhs_os[0].host(&apuF[inputScalars.lDimStruct.cumDim[ii]]);
					}
					largeDimLast(inputScalars, proj);
				}
			}
		}
		// Loop through each iteration
		for (uint32_t iter = iter0; iter < inputScalars.Niter; iter++) {

			if (inputScalars.verbose >= 3 || DEBUG)
				mexPrintVar("Starting iteration ", iter + 1);

			uint32_t curIter = 0;
			// Save the current iteration number if the current iteration is saved
			if (inputScalars.saveIter || inputScalars.saveIterationsMiddle > 0)
				curIter = iter;
			if (DEBUG) {
				mexPrintBase("curIter = %d\n", curIter);
				mexEval();
			}

			int64_t uu = 0;

			// Loop through the subsets
			for (uint32_t osa_iter = inputScalars.osa_iter0; osa_iter < inputScalars.subsetsUsed; osa_iter++) {
				inputScalars.currentSubset = osa_iter;
				if (inputScalars.stochastic) {
					osa_iter = distribution(rng);
					inputScalars.currentSubset = 0;
				}
				if (DEBUG || inputScalars.verbose >= 2) {
					af::sync();
					proj.tStartGlobal = std::chrono::steady_clock::now();
					if (DEBUG || inputScalars.verbose >= 3)
						mexPrintVar("Starting sub-iteration ", osa_iter + 1);
				}

				uint64_t m_size = length[osa_iter];
				uint32_t subIter = osa_iter;
				if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
					m_size = static_cast<uint64_t>(inputScalars.nRowsD) * static_cast<uint64_t>(inputScalars.nColsD) * length[osa_iter];

				// Load TOF/measurement data if it wasn't preloaded
				if (((osa_iter > inputScalars.osa_iter0 || iter > iter0) && !inputScalars.loadTOF) || inputScalars.largeDim) {
					if (inputScalars.subsetsUsed > 1 && (inputScalars.subsetType < 8))
						mData[0] = af::constant(0.f, lengthTOF[osa_iter]);
					else
						mData[0] = af::constant(0.f, inputScalars.nRowsD * inputScalars.nColsD * lengthTOF[osa_iter]);
					if (inputScalars.largeDim) {
						if (inputScalars.subsetsUsed > 1 && inputScalars.subsetType < 8) {
							const F** fptr = new const F* [length[osa_iter]];
							for (int ii = 0; ii < length[osa_iter]; ii++)
								fptr[ii] = const_cast<const F*>(&Sin[osa_iter + ii * inputScalars.subsetsUsed]);
							mData[0] = af::array(length[osa_iter], *fptr, AFTYPE).as(f32);
							af::sync();
							delete[] fptr;
						}
						else
							for (int ii = 0; ii < length[osa_iter]; ii++)
								mData[0](af::seq(inputScalars.nRowsD* inputScalars.nColsD* ii, inputScalars.nRowsD* inputScalars.nColsD* (ii + 1) - 1)) = af::array(inputScalars.nRowsD * inputScalars.nColsD, 
									&Sin[osa_iter * inputScalars.nRowsD * inputScalars.nColsD + inputScalars.nRowsD * inputScalars.nColsD * ii * inputScalars.subsetsUsed], afHost).as(f32);
						if (inputScalars.CT && MethodList.CPType)
							mData[0] = af::log(inputScalars.flat / mData[0]);
					}
#ifndef CPU
					else if (inputScalars.listmode > 0) {
						mData[0] = af::array(length[osa_iter], &Sin[pituus[osa_iter] + inputScalars.kokoNonTOF * tt], AFTYPE);
						if (inputScalars.indexBased)
							if (inputScalars.TOF)
								proj.loadCoord(inputScalars, length[osa_iter], &w_vec.trIndex[pituus[osa_iter] * 2], &w_vec.axIndex[pituus[osa_iter] * 2], &w_vec.TOFIndices[pituus[osa_iter]]);
							else
								proj.loadCoord(inputScalars, length[osa_iter], &w_vec.trIndex[pituus[osa_iter] * 2], &w_vec.axIndex[pituus[osa_iter] * 2]);
						else
							if (inputScalars.TOF)
								proj.loadCoord(inputScalars, length[osa_iter], &w_vec.listCoord[pituus[osa_iter] * 6], &w_vec.listCoord[pituus[osa_iter] * 6], &w_vec.TOFIndices[pituus[osa_iter]]);
							else
								proj.loadCoord(inputScalars, length[osa_iter], &w_vec.listCoord[pituus[osa_iter] * 6]);
					}
#endif
					else {
						if (inputScalars.subsetsUsed > 1 && inputScalars.subsetType < 8 && inputScalars.subsetType > 0)
							for (int64_t to = 0; to < inputScalars.nBins; to++)
								mData[0](af::seq(length[osa_iter] * to, length[osa_iter] * (to + 1) - 1)) = af::array(length[osa_iter], &Sin[pituus[osa_iter] + inputScalars.kokoNonTOF * to + inputScalars.kokoTOF * tt], AFTYPE);
						else
							for (int64_t to = 0; to < inputScalars.nBins; to++)
								mData[0](af::seq(inputScalars.nRowsD* inputScalars.nColsD* length[osa_iter] * to, inputScalars.nRowsD* inputScalars.nColsD* length[osa_iter] * (to + 1) - 1)) =
								af::array(length[osa_iter] * inputScalars.nRowsD * inputScalars.nColsD, &Sin[pituus[osa_iter] * inputScalars.nRowsD * inputScalars.nColsD + inputScalars.kokoNonTOF * to + inputScalars.kokoTOF * tt], AFTYPE);
					}
					if (inputScalars.randoms_correction) {
						if (inputScalars.subsetsUsed > 1 && (inputScalars.subsetType < 8 && inputScalars.subsetType > 0))
							aRand[0] = af::array(length[osa_iter], &sc_ra[pituus[osa_iter] + inputScalars.kokoNonTOF * tt], AFTYPE).as(f32);
						else
							aRand[0] = af::array(inputScalars.nRowsD * inputScalars.nColsD * length[osa_iter], &sc_ra[pituus[osa_iter] * inputScalars.nRowsD * inputScalars.nColsD + inputScalars.kokoNonTOF * tt], AFTYPE).as(f32);
					}
					subIter = 0U;
				}

				// Fill the sensitivity images with zeros, if necessary
				// Different versions depending on whether 64/32/float atomics are used
				if (compute_norm_matrix == 1u && proj.no_norm == 0) {
					for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
						if (inputScalars.atomic_64bit) {
							vec.Summ[ii][0] = af::constant(0LL, inputScalars.im_dim[ii], 1, s64);
						}
						else if (inputScalars.atomic_32bit) {
							vec.Summ[ii][0] = af::constant(0, inputScalars.im_dim[ii], 1, s32);
						}
						else
							vec.Summ[ii][0] = af::constant(0.f, inputScalars.im_dim[ii], 1);
					}
				}

				// Initialize some algorithms, such as the initial steps of LSQR
				for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
					status = initializationStep(w_vec, mData[subIter], vec, proj, inputScalars, length, m_size, MethodList, iter, meanBP, g, osa_iter, ii);
					if (status != 0)
						return -1;
				}

				af::sync();

				// Forward projections
				if (inputScalars.projector_type != 6) {

					float* FPapu = nullptr;
					if (!inputScalars.largeDim) {
						af::array outputFP;
						outputFP = af::constant(0.f, m_size * nBins);
						for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
							status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, osa_iter, length, g, m_size, proj, ii, pituus);
							if (status != 0) {
								return -1;
							}
						}
						af::sync();
						if (DEBUG) {
							mexPrintBase("outputFP.elements() = %d\n", outputFP.elements());
							mexPrintBase("mData[subIter].elements() = %d\n", mData[subIter].elements());
							mexPrintBase("outputFP = %f\n", af::sum<float>(outputFP));
							mexPrintBase("min(outputFP) = %f\n", af::min<float>(outputFP));
							mexPrintBase("max(outputFP) = %f\n", af::max<float>(outputFP));
							mexEval();
						}
						if (inputScalars.storeFP) {
							outputFP.eval();
							FPEstimates[iter][osa_iter].resize(outputFP.elements());
							FPapu = FPEstimates[iter][osa_iter].data();
							if (DEBUG) {
								mexPrintBase("FPEstimates[iter][osa_iter].size() = %d\n", FPEstimates[iter][osa_iter].size());
								mexEval();
							}
							outputFP.host(FPapu);
							outputFP.eval();
							af::sync();
						}
						//if (inputScalars.storeResidual) {
						//	residual[iter * inputScalars.subsets + osa_iter] = af::norm(outputFP - mData[subIter]);
						//	residual[iter * inputScalars.subsets + osa_iter] = residual[iter * inputScalars.subsets + osa_iter] * residual[iter * inputScalars.subsets + osa_iter] * .5f;
						//}
						status = computeForwardStep(MethodList, mData[subIter], outputFP, m_size, inputScalars, w_vec, aRand[subIter], vec, proj, iter, osa_iter, 0, residual);
						if (status != 0)
							return -1;
						if (DEBUG) {
							mexPrintBase("outputFP = %f\n", af::sum<float>(outputFP));
							mexPrintBase("max(outputFP) = %f\n", af::max<float>(outputFP));
							mexEval();
						}
						computeIntegralImage(inputScalars, w_vec, length[osa_iter], outputFP, meanBP);
						if (inputScalars.CT && (MethodList.ACOSEM || MethodList.OSLCOSEM > 0 || MethodList.OSEM || MethodList.COSEM || MethodList.ECOSEM ||
							MethodList.ROSEM || MethodList.OSLOSEM || MethodList.ROSEMMAP)) {
							if (inputScalars.randoms_correction) {
								OSEMapu = (mData[subIter] * outputFP) / (outputFP + aRand[subIter]);
								computeIntegralImage(inputScalars, w_vec, length[osa_iter], OSEMapu, meanBP);
								af::sync();
							}
							else if (iter == 0) {
								OSEMapu = mData[subIter];
								computeIntegralImage(inputScalars, w_vec, length[osa_iter], OSEMapu, meanBP);
								af::sync();
							}
						}
						for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
							if (ii == 0 && inputScalars.adaptiveType == 2 && MethodList.CPType)
								vec.adapTypeA = outputFP.copy();
							if (compute_norm_matrix == 1u)
								transferSensitivityImage(vec.Summ[ii][0], proj);
							else if (compute_norm_matrix == 2u)
								transferSensitivityImage(vec.Summ[ii][osa_iter], proj);

							if (inputScalars.CT && (MethodList.ACOSEM || MethodList.OSLCOSEM > 0 || MethodList.OSEM || MethodList.COSEM || MethodList.ECOSEM ||
								MethodList.ROSEM || MethodList.OSLOSEM || MethodList.ROSEMMAP)) {
								if (inputScalars.randoms_correction || iter == 0 || (compute_norm_matrix == 1u && proj.no_norm == 0)) {

									status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, OSEMapu, osa_iter, length, m_size, meanBP, g, proj, false, ii, pituus);
									if (compute_norm_matrix == 1u) {
										vec.Summ[ii][0] = vec.rhs_os[ii];
										vec.Summ[ii][0](vec.Summ[ii][0] < inputScalars.epps) = inputScalars.epps;
										vec.Summ[ii][0].eval();
									}
									else if (compute_norm_matrix == 2u) {
										vec.Summ[ii][osa_iter] = vec.rhs_os[ii];
										vec.Summ[ii][osa_iter](vec.Summ[ii][osa_iter] < inputScalars.epps) = inputScalars.epps;
										vec.Summ[ii][osa_iter].eval();
									}
								}
							}
							af::sync();

							status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, osa_iter, length, m_size, meanBP, g, proj, false, ii, pituus);
							if (status != 0) {
								if (compute_norm_matrix == 1u) {
									vec.Summ[ii][0].unlock();
								}
								else if (compute_norm_matrix == 2) {
									vec.Summ[ii][osa_iter].unlock();
								}
								return -1;
							}
							transferControl(vec, inputScalars, g, w_vec, compute_norm_matrix, proj.no_norm, osa_iter, ii);
						}
					}
					else {
						if (DEBUG)
							mexPrintVar("Starting largeDim ", inputScalars.largeDim);
						af::array outputFP = af::constant(0.f, m_size * nBins);
						for (int ii = 0; ii < inputScalars.subsetsUsed; ii++) {
							largeDimFirst(inputScalars, proj, ii);
							if (iter == 0 && osa_iter == 0)
								vec.im_os[0] = af::constant(1e-4f, inputScalars.lDimStruct.imDim[ii]);
							else
								vec.im_os[0] = af::array(inputScalars.lDimStruct.imDim[ii], &apuF[inputScalars.lDimStruct.cumDim[ii]], afHost);
							af::sync();
							status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, osa_iter, length, g, m_size, proj, 0, pituus);
							if (status != 0) {
								return -1;
							}
						}
						largeDimLast(inputScalars, proj);
						if (MethodList.PDHG || MethodList.PDHGKL || MethodList.PDHGL1 || MethodList.CV || MethodList.PDDY) {
							if (iter > 0) {
								if (inputScalars.subsetsUsed > 1 && (inputScalars.subsetType < 8))
									vec.pCP[0] = af::array(m_size * nBins, &apuM[pituus[osa_iter] * nBins], afHost);
								else
									vec.pCP[0] = af::array(m_size * nBins, &apuM[inputScalars.nRowsD * inputScalars.nColsD * pituus[osa_iter] * nBins], afHost);
							}
							else
								vec.pCP[0] = af::constant(0.f, m_size * nBins);
						}
						if (DEBUG || inputScalars.verbose >= 3) {
							proj.tStartLocal = std::chrono::steady_clock::now();
						}
						af::sync();
						status = computeForwardStep(MethodList, mData[0], outputFP, m_size, inputScalars, w_vec, aRand[0], vec, proj, iter, 0);
						if (status != 0)
							return -1;
						if (MethodList.PDHG || MethodList.PDHGKL || MethodList.PDHGL1 || MethodList.CV || MethodList.PDDY) {
							if (inputScalars.subsetsUsed > 1 && (inputScalars.subsetType < 8))
								vec.pCP[0].host(&apuM[pituus[osa_iter] * nBins]);
							else
								vec.pCP[0].host(&apuM[inputScalars.nRowsD * inputScalars.nColsD * pituus[osa_iter] * nBins]);
						}

						af::sync();

						for (int ii = 0; ii < inputScalars.subsetsUsed; ii++) {
							largeDimFirst(inputScalars, proj, ii);
							if (iter == 0 && osa_iter == 0) {
								vec.im_os[0] = af::constant(1e-4f, inputScalars.lDimStruct.imDim[ii]);
							}
							else
								vec.im_os[0] = af::array(inputScalars.lDimStruct.imDim[ii], &apuF[inputScalars.lDimStruct.cumDim[ii]], afHost);
							if (iter == 0 && osa_iter == 0 && (MethodList.PDHG || MethodList.PDHGKL || MethodList.PDHGL1 || MethodList.CV || MethodList.PDDY))
								vec.uCP[0] = vec.im_os[0].copy();
							else if (MethodList.PDHG || MethodList.PDHGKL || MethodList.PDHGL1 || MethodList.CV || MethodList.PDDY)
								vec.uCP[0] = af::array(inputScalars.lDimStruct.imDim[ii], &apuU[inputScalars.lDimStruct.cumDim[ii]], afHost);
							af::sync();
							status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, osa_iter, length, m_size, meanBP, g, proj, false, 0, pituus);
							status = computeOSEstimates(vec, w_vec, MethodList, iter, osa_iter, inputScalars, length, break_iter,
								pituus, g, proj, mData[0], m_size, uu, compute_norm_matrix);
							if (status != 0)
								return -1;
							if (inputScalars.enforcePositivity && inputScalars.subsetsUsed > 1 && !MethodList.CPType) {
								for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++)
									vec.im_os[ii](vec.im_os[ii] < inputScalars.epps) = inputScalars.epps;
							}
							af::sync();
							vec.im_os[0].host(&apuF[inputScalars.lDimStruct.cumDim[ii]]);
							if (MethodList.PDHG || MethodList.PDHGKL || MethodList.PDHGL1 || MethodList.CV || MethodList.PDDY)
								vec.uCP[0].host(&apuU[inputScalars.lDimStruct.cumDim[ii]]);
						}
						largeDimLast(inputScalars, proj);
					}
					af::sync();
				}
				// SPECT rotation-based projector
				else if (inputScalars.projector_type == 6) {

					float* FPapu = nullptr;
					af::sync();
					af::array fProj = af::constant(0.f, inputScalars.nRowsD, inputScalars.nColsD, length[osa_iter]);
					if (DEBUG) {
						mexPrintBase("length[osa_iter] = %d\n", length[osa_iter]);
						mexPrintBase("uu = %d\n", uu);
						mexEval();
					}
					for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++)
						forwardProjectionType6(fProj, w_vec, vec, inputScalars, length[osa_iter], uu, proj, ii, atten);
					fProj.eval();
					fProj = af::flat(fProj);
					fProj(fProj < inputScalars.epps) = inputScalars.epps;
					af::sync();
					if (inputScalars.storeFP) {
						fProj.eval();
						FPEstimates[iter][osa_iter].resize(fProj.elements());
						FPapu = FPEstimates[iter][osa_iter].data();
						if (DEBUG) {
							mexPrintBase("FPEstimates[iter][osa_iter].size() = %d\n", FPEstimates[iter][osa_iter].size());
							mexPrintBase("fProj.elements() = %d\n", fProj.elements());
							mexEval();
						}
						fProj.host(FPapu);
					}
					af::sync();
					status = computeForwardStep(MethodList, mData[subIter], fProj, m_size, inputScalars, w_vec, aRand[subIter], vec, proj, iter, osa_iter);
					if (status != 0)
						return -1;
					for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++)
						backprojectionType6(fProj, w_vec, vec, inputScalars, length[osa_iter], uu, proj, osa_iter, iter, compute_norm_matrix, iter0, ii, atten);
				}
				af::sync();
				if (DEBUG) {
					if (compute_norm_matrix > 1u) {
						mexPrintBase("Summ = %f\n", af::sum<float>(vec.Summ[0][osa_iter]));
						mexPrintBase("min(Summ) = %f\n", af::min<float>(vec.Summ[0][osa_iter]));
					}
					mexPrintBase("vec.im_os = %f\n", af::sum<float>(vec.im_os[0]));
					mexPrintBase("vec.rhs_os = %f\n", af::sum<float>(vec.rhs_os[0]));
					mexPrintBase("min(rhs_os) = %f\n", af::min<float>(vec.rhs_os[0]));
					mexEval();
				}
				if (DEBUG) {
					mexPrint("Algo\n");
				}
				if (!inputScalars.largeDim) {
					status = computeOSEstimates(vec, w_vec, MethodList, iter, osa_iter, inputScalars, length, break_iter,
						pituus, g, proj, mData[subIter], m_size, uu, compute_norm_matrix);
					if (status != 0)
						return -1;
					if (DEBUG) {
						mexPrintBase("vec.im_os = %f\n", af::sum<float>(vec.im_os[0]));
						mexEval();
						if (af::anyTrue<bool>(af::isNaN(vec.im_os[0])))
							return -1;
						if (af::anyTrue<bool>(af::isInf(vec.im_os[0])))
							return -1;
					}

					// Enforce positivity if applicable
					if (inputScalars.enforcePositivity && inputScalars.subsetsUsed > 1 && !MethodList.CPType && !MethodList.POCS) {
						for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++)
							vec.im_os[ii](vec.im_os[ii] < inputScalars.epps) = inputScalars.epps;
					}


				}

				if (inputScalars.verbose > 0 && inputScalars.subsetsUsed > 1 && inputScalars.stochastic == false) {
					
					if (DEBUG || inputScalars.verbose >= 2) {
						af::sync();
						proj.tEndGlobal = std::chrono::steady_clock::now();
						const std::chrono::duration<double> tDiff = proj.tEndGlobal - proj.tStartGlobal;
						mexPrintBase("Sub-iteration %d complete in %f seconds\n", osa_iter + 1u, tDiff);
						totTime += tDiff.count();
					}
					else
						mexPrintBase("Sub-iteration %d complete\n", osa_iter + 1u);
					mexEval();
				}

				if (inputScalars.projector_type == 6)
					uu += length[osa_iter];

				proj.memSize -= (sizeof(float) * totDim) / 1048576ULL;

				if (break_iter || inputScalars.stochastic)
					break;

				af::deviceGC();
			}
			// Compute some subset-based algorithms that require special operations after each iteration such as BSREM or ROSEM
			// Also copy the current iteration if needed
			if (!inputScalars.largeDim) {
				status = computeOSEstimatesIter(vec, w_vec, MethodList, inputScalars, iter, proj, g, cell, ee, eInd, x0);
				if (status != 0)
					return -1;
				if (inputScalars.enforcePositivity && inputScalars.subsetsUsed == 1 && !MethodList.CPType) {
					for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++)
						vec.im_os[ii](vec.im_os[ii] < inputScalars.epps) = inputScalars.epps;
				}
			}

			if (compute_norm_matrix == 2u)
				proj.no_norm = 1u;

			if (inputScalars.verbose > 0) {
				if (DEBUG || inputScalars.verbose >= 2) {
					af::sync();
					iterTime = totTime - cumulativeTime;
					cumulativeTime += iterTime;
					mexPrintBase("Iteration %d complete in %f seconds\n", iter + 1u, iterTime);
					mexPrintBase("Estimated time left: %f seconds\n", iterTime * (inputScalars.Niter - 1 - iter));
				}
				else
					mexPrintBase("Iteration %d complete\n", iter + 1u);
				mexEval();
			}
			af::deviceGC();
			if (break_iter)
				break;
		}

		if (!inputScalars.largeDim) {
#ifdef MATLAB
			// Transfer the device data to host MATLAB mxarray
			device_to_host(MethodList, vec, oo, cell, FPcell, w_vec, 4, inputScalars, FPEstimates);
			oo += inputScalars.im_dim[0] * static_cast<int64_t>(nIter);
#else
			device_to_host(MethodList, vec, oo, cell, FPcell, inputScalars, FPEstimates);
#endif
		}

		if (w_vec.filteringOrig)
			w_vec.precondTypeMeas[1] = w_vec.filteringOrig;

		if (inputScalars.verbose > 0 && inputScalars.listmode != 2 && inputScalars.Nt > 1) {
			mexPrintBase("Time step %d complete\n", tt + 1u);
			mexEval();
		}
		if (inputScalars.largeDim) {
			if (MethodList.PDHG || MethodList.PDHGKL || MethodList.PDHGL1 || MethodList.CV || MethodList.PDDY) {
				delete[] apuM;
				delete[] apuU;
			}
		}
	}
	if (w_vec.computeD)
		w_vec.D.clear();
	if (w_vec.computeM)
		w_vec.M.clear();
	af::sync();
	af::deviceGC();

	if (inputScalars.verbose >= 2 || DEBUG) {
		mexPrintBase("Reconstruction complete in %f seconds\n", totTime);
	}

	return 0;
}