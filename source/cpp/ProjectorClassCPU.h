/*******************************************************************************************************************************************
* Class object for forward and backward projections.
*
* Copyright (C) 2022-2024 Ville-Veikko Wettenhovi, Niilo Saarlemo
*
* This program is free software: you can redistribute it and/or modify  it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or  (at your option) any later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License  along with this program. If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************************************************************************/
#pragma once
#include "projector_functions.h"
#include "structs.h"
#include "RDP.h"
#include "NLM.h"
#include "TV.h"
#include "MRP.h"
#include "GGMRF.h"
#include "hyperbolicPrior.h"

// #define DEBUG false

class ProjectorClass {
	paramStruct<float> param;

public:
	uint8_t no_norm = 0;
	// d_xFull is the detector (crystal) coordinate table and is never overwritten (GPU d_x[0][0] / d_xFull). The list-mode event
	// coordinates of the current subset are a local pointer in forwardProjection/backwardProjection, so that the list-mode sensitivity
	// image (compSens) always uses the detector table, like the GPU (d_xFull)
	float* d_xFull = nullptr, * d_z, * input, * output, * SensIm = nullptr, *d_norm = nullptr, *d_atten = nullptr, * extraCor = nullptr;
	uint16_t* detIndices = nullptr, *zIndex = nullptr;
	uint32_t* xyIndex = nullptr;
	CPUVectors vec_opencl;
	std::vector<float*> d_Summ = { nullptr };
	float* d_output, * d_W, * d_gaussianNLM, * d_inputB, * d_refIm, * weights;
	size_t memSize = 0ULL;
	int32_t Ndx, Ndy, Ndz;
	// Prior mask / extended FOV slice indices, mirroring the mask buffers the GPU projector class holds
	const uint8_t* d_maskPrior = nullptr, * d_eFOVIndices = nullptr;
	uint32_t maskPriorZ = 1U;
	std::chrono::steady_clock::time_point tStartLocal, tStartGlobal, tStartAll;
	std::chrono::steady_clock::time_point tEndLocal, tEndGlobal, tEndAll;

	~ProjectorClass() {	}

	inline int addProjector(scalarStruct& inputScalars, Weighting& w_vec, const RecMethods& MethodList, const char* header_directory = "", const int type = -1) {
		// Number of detector indices
		param.size_x = inputScalars.nRowsD;
		param.attenuationCorrection = inputScalars.attenuation_correction;
		// Image-based (true) or measurement-based (false, GPU -DATNM) attenuation correction
		param.CTAttenuation = inputScalars.CTAttenuation;
		// Total length (true, GPU -DTOTLENGTH) or the in-FOV length (false) of the ray is used in the normalization of the PET rays
		param.useTotLength = inputScalars.useTotLength;
		param.normalizationCorrection = inputScalars.normalization_correction;
		param.scatterCorrectionMult = inputScalars.scatter;
		param.globalFactor = inputScalars.global_factor;
		param.det_per_ring = inputScalars.det_per_ring;
		param.TOF = inputScalars.TOF;
		param.sigma_x = inputScalars.sigma_x;
		param.TOFCenters = inputScalars.TOFCenter;
		param.nBins = inputScalars.nBins;
		param.raw = inputScalars.raw;
		param.listMode = inputScalars.listmode;
		param.projType = inputScalars.projector_type;
		param.pitch = inputScalars.pitch;
		param.subsetType = inputScalars.subsetType;
		param.subsets = inputScalars.subsetsUsed;
		param.size_y = inputScalars.nColsD;
		param.dPitchXY = w_vec.dPitchX;
		param.nRays2D = inputScalars.n_rays;
		param.nRays3D = inputScalars.n_rays3D;
		if (inputScalars.nLayers > 1)
			param.nLayers = true;
		// Number of detector layers, the list-mode sensitivity image enumerates all nLayers^2 layer pairs (GPU -DNLAYERS)
		param.nLayerCount = inputScalars.nLayers > 1 ? inputScalars.nLayers : 1U;
		param.computeSensIm = inputScalars.computeSensImag;
		param.rings = inputScalars.rings;
		param.epps = inputScalars.epps;
		param.dPitchZ = w_vec.dPitchY;
		param.nProjections = inputScalars.nProjections;
		param.useMaskBP = inputScalars.maskBP;
		param.useMaskFP = inputScalars.maskFP;
		if (param.useMaskFP)
			param.maskFP = w_vec.maskFP;
		if (param.useMaskBP)
			param.maskBP = w_vec.maskBP;
		if (inputScalars.priorMask) {
			d_maskPrior = w_vec.maskPrior;
			maskPriorZ = inputScalars.maskBPZ;
		}
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			d_eFOVIndices = w_vec.eFOVIndices;
		// SPECT EDIT
		param.rayShiftsDetector = w_vec.rayShiftsDetector;
		param.rayShiftsSource = w_vec.rayShiftsSource;
        param.coneOfResponseStdCoeffA = inputScalars.coneOfResponseStdCoeffA;
        param.coneOfResponseStdCoeffB = inputScalars.coneOfResponseStdCoeffB;
        param.coneOfResponseStdCoeffC = inputScalars.coneOfResponseStdCoeffC;
		// END SPECT EDIT
		if (inputScalars.SPECT) {
			// Support the SPECT rays are clipped to (GPU: extendRayToEllipse), ellipsePower is FLT_MAX for the box support
			param.ellipseCenterX = inputScalars.ellipseCenterX;
			param.ellipseCenterY = inputScalars.ellipseCenterY;
			param.ellipseCenterZ = inputScalars.ellipseCenterZ;
			param.ellipseRadiusX = inputScalars.ellipseRadiusX;
			param.ellipseRadiusY = inputScalars.ellipseRadiusY;
			param.ellipseRadiusZ = inputScalars.ellipseRadiusZ;
			param.ellipsePower = inputScalars.ellipsePower;
		}
		return 0;
	}

	// SPECT: the detector head vector of the current subset/time step (selects the ray shifts of the head). The detector head vector is
	// concatenated over the time steps and the current subset of the time step starts at pituus[osa_iter + timestep * subsets] (GPU:
	// d_detectorVector). The projection index of the projector already includes the subset offset pituus[osa_iter] (the same index as the
	// detector coordinates)
	inline void setSPECTData(const scalarStruct& inputScalars, const Weighting& w_vec, const uint32_t osa_iter, const uint32_t timestep, const int64_t* pituus) {
		if (!inputScalars.SPECT)
			return;
		if (w_vec.detectorVector != nullptr)
			param.detectorVector = w_vec.detectorVector + (pituus[osa_iter + static_cast<size_t>(timestep) * inputScalars.subsets] - pituus[osa_iter]);
		else
			param.detectorVector = nullptr;
	}

	// The attenuation data of the current subset/time step, the same data as the GPU projector class uses:
	// measurement-based attenuation (CTAttenuation == false) has one value per measurement, the current subset starts at the same
	// offset as the normalization data. Attenuation images are per time step if more than one image is input
	inline void setAttenuation(const scalarStruct& inputScalars, const uint32_t osa_iter, const uint32_t timestep, const int64_t* pituus, const size_t vecSize) {
		if (!inputScalars.attenuation_correction)
			return;
		if (!inputScalars.CTAttenuation)
			param.atten = &d_atten[pituus[osa_iter] * vecSize];
		else if (inputScalars.size_atten > static_cast<size_t>(inputScalars.im_dim[0]))
			param.atten = &d_atten[static_cast<size_t>(inputScalars.im_dim[0]) * static_cast<size_t>(timestep)];
		else
			param.atten = d_atten;
	}

	inline int createBuffers(scalarStruct& inputScalars, Weighting& w_vec, const float* x, const float* z_det, const uint32_t* xy_index,
		const uint16_t* z_index, const uint16_t* L, const int64_t* pituus, const float* atten, const float* norm, const float* extraCorr,
		const std::vector<int64_t>& length, const RecMethods& MethodList, const int type = 0) {
		d_xFull = const_cast<float*>(x);
		d_z = const_cast<float*>(z_det);
		if (inputScalars.size_norm > 1ULL && inputScalars.normalization_correction)
			d_norm = const_cast<float*>(norm);
		if (inputScalars.attenuation_correction)
			d_atten = const_cast<float*>(atten);
		if (inputScalars.size_scat > 1ULL && inputScalars.scatter == 1U)
			extraCor = const_cast<float*>(extraCorr);
		if (inputScalars.listmode != 1 && ((!inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET) && (inputScalars.subsets > 1 && (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7)))) {
			xyIndex = const_cast<uint32_t*>(xy_index);
			zIndex = const_cast<uint16_t*>(z_index);
		}
		if (MethodList.GGMRF || MethodList.hyperbolic) {
			weights = w_vec.weights;
			Ndx = w_vec.Ndx;
			Ndy = w_vec.Ndy;
			Ndz = w_vec.Ndz;
		}
		return 0;
	}

	inline int forwardProjection(const scalarStruct& inputScalars, Weighting& w_vec, const uint32_t osa_iter, uint32_t timestep, const std::vector<int64_t>& length, const int64_t* pituus, const int ii = 0, const int uu = 0) {
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrintVar("Starting forward projection for projector type = ", inputScalars.FPType);
		size_t vecSize = 1;
		if ((inputScalars.PET || inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0)
			vecSize = static_cast<size_t>(inputScalars.nRowsD) * static_cast<size_t>(inputScalars.nColsD);
		if (DEBUG) {
			mexPrintBase("FPType = %u\n", inputScalars.FPType);
			mexPrintBase("nRowsD = %u\n", inputScalars.nRowsD);
			mexPrintBase("nColsD = %u\n", inputScalars.nColsD);
			mexPrintBase("length[osa_iter] = %u\n", length[osa_iter]);
			mexPrintBase("pituus[osa_iter] = %u\n", pituus[osa_iter]);
			mexPrintBase("listmode = %u\n", inputScalars.listmode);
			mexPrintBase("no_norm = %u\n", no_norm);
			mexPrintBase("vecSize = %u\n", vecSize);
			mexPrintBase("ii = %u\n", ii);
			mexPrintBase("osa_iter = %u\n", osa_iter);
			mexPrintBase("subsetType = %u\n", inputScalars.subsetType);
			mexEval();
		}
		param.Nx = inputScalars.Nx[ii];
		param.Ny = inputScalars.Ny[ii];
		param.Nz = inputScalars.Nz[ii];
		param.dx = inputScalars.dx[ii];
		param.dy = inputScalars.dy[ii];
		param.dz = inputScalars.dz[ii];
		param.bx = inputScalars.bx[ii];
		param.by = inputScalars.by[ii];
		param.bz = inputScalars.bz[ii];
		setAttenuation(inputScalars, osa_iter, timestep, pituus, vecSize);
		if (inputScalars.size_norm > 1ULL && inputScalars.normalization_correction)
			param.normCoefs = &d_norm[pituus[osa_iter] * vecSize];
		if (inputScalars.size_scat > 1ULL && inputScalars.scatter == 1U)
			param.scatterCoefs = &extraCor[pituus[osa_iter] * vecSize];
		if (inputScalars.listmode != 1 && ((!inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET) && (inputScalars.subsets > 1 && (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7)))) {
			param.xy_index = &xyIndex[pituus[osa_iter] * vecSize];
			param.z_index = &zIndex[pituus[osa_iter] * vecSize];
		}
		// List-mode: the event coordinates of the current subset, otherwise the detector coordinate table
		const float* xCoord = d_xFull;
		if (param.listMode > 0)
			xCoord = &w_vec.listCoord[pituus[osa_iter] * 6];
		param.currentSubset = osa_iter;
		// Offset (in projections/sinograms) to the first measurement of the current subset, required by subset types >= 8
		param.nMeas = pituus[osa_iter];
		param.computeSensIm = false;
		param.projType = inputScalars.FPType;
		setSPECTData(inputScalars, w_vec, osa_iter, timestep, pituus);

		projectorType123Implementation4<float>(param, length[osa_iter] * vecSize, d_output, xCoord, d_z, vec_opencl.d_im_os, inputScalars.CT, inputScalars.SPECT, 1, SensIm, detIndices);

		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Forward projection completed");
		return 0;
	}

	inline int backwardProjection(const scalarStruct& inputScalars, Weighting& w_vec, const uint32_t osa_iter, uint32_t timestep, const std::vector<int64_t>& length, const int64_t* pituus, const bool compSens = false, const int ii = 0, const int uu = 0) {
		if (inputScalars.verbose >= 3)
			mexPrintVar("Starting backprojection for projector type = ", inputScalars.BPType);
		size_t vecSize = 1;
		if ((inputScalars.PET || inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0)
			vecSize = static_cast<size_t>(inputScalars.nRowsD) * static_cast<size_t>(inputScalars.nColsD);
		int64_t nMeas = length[osa_iter] * vecSize;
		// List-mode sensitivity image: all ordered detector pairs of all ring pairs and layer pairs (64-bit count, the same m_size as the GPU,
		// reconstructionAF.h), the pair index is enumerated in projectorType123Implementation4 as in projectorType123.cl (SENS)
		if (compSens && param.listMode > 0) {
			uint64_t nPairs = static_cast<uint64_t>(inputScalars.det_per_ring) * static_cast<uint64_t>(inputScalars.det_per_ring) *
				static_cast<uint64_t>(inputScalars.rings) * static_cast<uint64_t>(inputScalars.rings);
			if (inputScalars.nLayers > 1)
				nPairs *= static_cast<uint64_t>(inputScalars.nLayers) * static_cast<uint64_t>(inputScalars.nLayers);
			nMeas = static_cast<int64_t>(nPairs);
		}
		// The detector coordinate table, or the list-mode event coordinates of the current subset (set below)
		const float* xCoord = d_xFull;
		if (DEBUG) {
			mexPrintBase("nMeas = %u\n", nMeas);
			mexPrintBase("no_norm = %u\n", no_norm);
			mexPrintBase("compSens = %u\n", compSens);
			mexPrintBase("vecSize = %u\n", vecSize);
			mexPrintBase("osa_iter = %u\n", osa_iter);
			mexPrintBase("length[osa_iter] = %u\n", length[osa_iter]);
			mexEval();
		}
		param.Nx = inputScalars.Nx[ii];
		param.Ny = inputScalars.Ny[ii];
		param.Nz = inputScalars.Nz[ii];
		param.dx = inputScalars.dx[ii];
		param.dy = inputScalars.dy[ii];
		param.dz = inputScalars.dz[ii];
		param.bx = inputScalars.bx[ii];
		param.by = inputScalars.by[ii];
		param.bz = inputScalars.bz[ii];
		setAttenuation(inputScalars, osa_iter, timestep, pituus, vecSize);
		if (compSens) {
			if (inputScalars.size_norm > 1ULL && inputScalars.normalization_correction)
				param.normCoefs = d_norm;
			if (inputScalars.size_scat > 1ULL && inputScalars.scatter == 1U)
				param.scatterCoefs = extraCor;
			param.computeSensIm = compSens;
			param.noSensImage = compSens;
		}
		else {
			if (inputScalars.size_norm > 1ULL && inputScalars.normalization_correction)
				param.normCoefs = &d_norm[pituus[osa_iter] * vecSize];
			if (inputScalars.size_scat > 1ULL && inputScalars.scatter == 1U)
				param.scatterCoefs = &extraCor[pituus[osa_iter] * vecSize];
			if (inputScalars.listmode != 1 && ((!inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET) && (inputScalars.subsets > 1 && (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7)))) {
				param.xy_index = &xyIndex[pituus[osa_iter] * vecSize];
				param.z_index = &zIndex[pituus[osa_iter] * vecSize];
			}
			if (param.listMode > 0)
				xCoord = &w_vec.listCoord[pituus[osa_iter] * 6];
			param.noSensImage = no_norm;
			param.computeSensIm = false;
		}
		param.projType = inputScalars.BPType;
		param.currentSubset = osa_iter;
		// Offset (in projections/sinograms) to the first measurement of the current subset, required by subset types >= 8
		param.nMeas = pituus[osa_iter];
		if (!compSens)
			setSPECTData(inputScalars, w_vec, osa_iter, timestep, pituus);

		projectorType123Implementation4<float>(param, nMeas, vec_opencl.d_rhs_os[ii], xCoord, d_z, d_output, inputScalars.CT, inputScalars.SPECT, 2, d_Summ[uu], detIndices);
		return 0;
	}

	inline int computeMRP(const scalarStruct& inputScalars, const int32_t Ndx, const int32_t Ndy, const int32_t Ndz) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting OpenMP median root prior computation");
		medianFilter3D(d_inputB, d_W, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], inputScalars.NxOrig, inputScalars.NyOrig, inputScalars.NzOrig, Ndx, Ndy, Ndz,
			d_maskPrior, d_eFOVIndices, maskPriorZ);
		if (inputScalars.verbose >= 3)
			mexPrint("OpenMP median root prior computed");
		return 0;
	}

	inline int computeNLM(const scalarStruct& inputScalars, Weighting& w_vec, const float beta) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting OpenMP NLM gradient computation");
		int type = 0;
		if (w_vec.NLM_MRP)
			type = 2;
		else if (w_vec.NLTV)
			type = 1;
		else if (w_vec.NLRD)
			type = 3;
		else if (w_vec.NLLange)
			type = 4;
		else if (w_vec.NLLangeFiltered)
			type = 5;
		else if (w_vec.NLGGMRF)
			type = 6;
		else if (w_vec.NLGM)
			type = 7;
		else
			type = 0;

		// Input padded sizes
		const uint32_t pNx = inputScalars.Nx[0] + 2u * (w_vec.Ndx + w_vec.Nlx);
		const uint32_t pNy = inputScalars.Ny[0] + 2u * (w_vec.Ndy + w_vec.Nly);
		const uint32_t pNz = inputScalars.Nz[0] + 2u * (w_vec.Ndz + w_vec.Nlz);
		NLMFunc(d_W, w_vec.NLM_ref, d_inputB, d_gaussianNLM, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, w_vec.Nlx, w_vec.Nly, w_vec.Nlz, pNx,
			pNy, pNz, static_cast<int64_t>(pNx) * static_cast<int64_t>(pNy), w_vec.h2, type, w_vec.RDP_gamma, inputScalars.epps, w_vec.GGMRF_p, w_vec.GGMRF_q, w_vec.GGMRF_c,
			w_vec.NLM_anatomical, d_maskPrior, d_eFOVIndices, maskPriorZ);
		if (inputScalars.verbose >= 3)
			mexPrint("OpenMP NLM gradient computed");
		return 0;
	}

	inline int computeRDP(const scalarStruct& inputScalars, const float gamma, const float beta, const bool RDPLargeNeighbor = false, const bool useRDPRef = false) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting OpenMP RDP gradient computation");
		RDPKernel(d_W, d_inputB, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], inputScalars.NxOrig, inputScalars.NyOrig, inputScalars.NzOrig, gamma, inputScalars.epps, beta, d_maskPrior, d_eFOVIndices, maskPriorZ);
		if (inputScalars.verbose >= 3)
			mexPrint("OpenMP RDP gradient computed");
		return 0;
	}
	inline int TVGradient(const scalarStruct& inputScalars, const TVdata& data, const float sigma, const float smooth, const float beta, const float C = 0.f, const int type = 0) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting OpenMP TV gradient computation");
		if (DEBUG) {
			mexPrintBase("type = %u\n", type);
			mexPrintBase("data.TV_use_anatomical = %u\n", data.TV_use_anatomical);
			mexEval();
		}
		// With an anatomical reference image, the GPU only defines the reference image variants (TV types 1, 2 and 5). TV types 4 (SATV) and 6 (TVW1)
		// are therefore computed as the standard TV, same as the GPU
		int tvType = type;
		if (data.TV_use_anatomical && (type == 4 || type == 6))
			tvType = 0;
		TVKernel(d_W, d_inputB, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], inputScalars.NxOrig, inputScalars.NyOrig, inputScalars.NzOrig, sigma, smooth, beta, tvType, data.TV_use_anatomical, C, d_refIm,
			d_maskPrior, d_eFOVIndices, maskPriorZ);
		if (inputScalars.verbose >= 3)
			mexPrint("OpenMP TV gradient computed");
		return 0;
	}

	inline int computeGGMRF(const scalarStruct& inputScalars, const float p, const float q, const float c, const float pqc, const float beta) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting OpenMP GGMRF gradient computation");
		GGMRFKernel(d_W, d_inputB, weights, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], p, q, c, pqc, beta, Ndx, Ndy, Ndz, d_maskPrior, d_eFOVIndices, maskPriorZ);
		if (inputScalars.verbose >= 3)
			mexPrint("OpenMP GGMRF gradient computed");
		return 0;
	}

	inline int hyperGradient(const scalarStruct& inputScalars, const float sigma, const float beta) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting OpenMP hyperbolic gradient computation");
		hyperbolicKernel(d_W, d_inputB, weights, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], sigma, beta, Ndx, Ndy, Ndz, d_maskPrior, d_eFOVIndices, maskPriorZ);
		if (inputScalars.verbose >= 3)
			mexPrint("OpenMP hyperbolic gradient computed");
		return 0;
	}

};