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
	float* d_x, * d_z, * input, * output, * SensIm = nullptr, *d_norm = nullptr, *d_atten = nullptr, * extraCor = nullptr;
	uint16_t* detIndices = nullptr, *zIndex = nullptr;
	uint32_t* xyIndex = nullptr;
	CPUVectors vec_opencl;
	std::vector<float*> d_Summ = { nullptr };
	float* d_output, * d_W, * d_gaussianNLM, * d_inputB, * d_refIm, * weights;
	size_t memSize = 0ULL;
	int32_t Ndx, Ndy, Ndz;

	~ProjectorClass() {	}

	inline int addProjector(scalarStruct& inputScalars, Weighting& w_vec, const RecMethods& MethodList, const char* header_directory = "", const int type = -1) {
		// Number of detector indices
		param.size_x = inputScalars.nRowsD;
		param.attenuationCorrection = inputScalars.attenuation_correction;
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
		param.computeSensIm = inputScalars.computeSensImag;
		param.rings = inputScalars.rings;
		param.x_center = inputScalars.x_center;
		param.y_center = inputScalars.y_center;
		param.z_center = inputScalars.z_center;
		param.bmin = inputScalars.bmin;
		param.bmax = inputScalars.bmax;
		param.Vmax = inputScalars.Vmax;
		if (param.projType == 2 || param.projType == 21 || param.projType == 12)
			param.orthWidth = inputScalars.orthXY;
		else if (param.projType == 3 || param.projType == 31 || param.projType == 13) {
			param.orthWidth = inputScalars.cylRadiusProj3;
			param.V = inputScalars.V;
		}
		param.epps = inputScalars.epps;
		param.dPitchZ = w_vec.dPitchY;
		param.nProjections = inputScalars.nProjections;
		param.useMaskBP = inputScalars.maskBP;
		param.useMaskFP = inputScalars.maskFP;
		if (param.useMaskFP)
			param.maskFP = w_vec.maskFP;
		if (param.useMaskBP)
			param.maskBP = w_vec.maskBP;
		// SPECT EDIT
		param.colL = inputScalars.colL;
		param.colD = inputScalars.colD;
		param.dSeptal = inputScalars.dSeptal;
		param.nRaySPECT = inputScalars.nRaySPECT;
		param.hexOrientation = inputScalars.hexOrientation;
		param.coneMethod = inputScalars.coneMethod;
		// END SPECT EDIT
		return 0;
	}

	inline int createBuffers(scalarStruct& inputScalars, Weighting& w_vec, const float* x, const float* z_det, const uint32_t* xy_index,
		const uint16_t* z_index, const uint16_t* L, const int64_t* pituus, const float* atten, const float* norm, const float* extraCorr,
		const std::vector<int64_t>& length, const RecMethods& MethodList, const int type = 0) {
		d_x = const_cast<float*>(x);
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

	inline int forwardProjection(const scalarStruct& inputScalars, Weighting& w_vec, const uint32_t osa_iter, const std::vector<int64_t>& length, const int64_t* pituus, const int ii = 0, const int uu = 0) {
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrintVar("Starting forward projection for projector type = ", inputScalars.FPType);
		size_t vecSize = 1;
		if ((inputScalars.PET || inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0)
			vecSize = static_cast<size_t>(inputScalars.nRowsD) * static_cast<size_t>(inputScalars.nColsD);
		if (DEBUG) {
			mexPrintBase("nRaySPECT = %f\n", inputScalars.nRaySPECT);
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
		if (inputScalars.attenuation_correction)
			param.atten = d_atten;
		if (inputScalars.size_norm > 1ULL && inputScalars.normalization_correction)
			param.normCoefs = &d_norm[pituus[osa_iter] * vecSize];
		if (inputScalars.size_scat > 1ULL && inputScalars.scatter == 1U)
			param.scatterCoefs = &extraCor[pituus[osa_iter] * vecSize];
		if (inputScalars.listmode != 1 && ((!inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET) && (inputScalars.subsets > 1 && (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7)))) {
			param.xy_index = &xyIndex[pituus[osa_iter] * vecSize];
			param.z_index = &zIndex[pituus[osa_iter] * vecSize];
		}
		if (param.listMode > 0)
			d_x = &w_vec.listCoord[pituus[osa_iter] * 6];
		param.currentSubset = osa_iter;
		param.nMeas = length[osa_iter];
		param.computeSensIm = false;
		param.projType = inputScalars.FPType;

		projectorType123Implementation4<float>(param, length[osa_iter] * vecSize, d_output, d_x, d_z, vec_opencl.d_im_os, inputScalars.CT, inputScalars.SPECT, 1, SensIm, detIndices);

		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Forward projection completed");
		return 0;
	}

	inline int backwardProjection(const scalarStruct& inputScalars, Weighting& w_vec, const uint32_t osa_iter, const std::vector<int64_t>& length, const int64_t* pituus, const bool compSens = false, const int ii = 0, const int uu = 0) {
		if (inputScalars.verbose >= 3)
			mexPrintVar("Starting backprojection for projector type = ", inputScalars.BPType);
		size_t vecSize = 1;
		if ((inputScalars.PET || inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0)
			vecSize = static_cast<size_t>(inputScalars.nRowsD) * static_cast<size_t>(inputScalars.nColsD);
		int64_t nMeas = length[osa_iter] * vecSize;
		if (compSens)
			nMeas = inputScalars.det_per_ring * inputScalars.det_per_ring * inputScalars.rings;
		if (DEBUG) {
			mexPrintBase("nMeas = %u\n", nMeas);
			mexPrintBase("no_norm = %u\n", no_norm);
			mexPrintBase("compSens = %u\n", compSens);
			mexPrintBase("vecSize = %u\n", vecSize);
			mexPrintBase("osa_iter = %u\n", osa_iter);
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
		if (inputScalars.attenuation_correction)
			param.atten = d_atten;
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
				d_x = &w_vec.listCoord[pituus[osa_iter] * 6];
			param.noSensImage = no_norm;
			param.computeSensIm = false;
		}
		param.projType = inputScalars.BPType;
		param.currentSubset = osa_iter;
		param.nMeas = length[osa_iter];

		projectorType123Implementation4<float>(param, nMeas, vec_opencl.d_rhs_os[ii], d_x, d_z, d_output, inputScalars.CT, inputScalars.SPECT, 2, d_Summ[uu], detIndices);
		return 0;
	}

	inline int computeMRP(const scalarStruct& inputScalars, const int32_t Ndx, const int32_t Ndy, const int32_t Ndz) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting OpenMP median root prior computation");
		medianFilter3D(d_inputB, d_W, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], inputScalars.NxOrig, inputScalars.NyOrig, inputScalars.NzOrig, Ndx, Ndy, Ndz);
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
		else
			type = 0;

		NLMFunc(d_W, w_vec.NLM_ref, d_inputB, d_gaussianNLM, w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, w_vec.Nlx, w_vec.Nly, w_vec.Nlz, inputScalars.Nx[0],
			inputScalars.Ny[0], inputScalars.Nz[0], inputScalars.Nxy, w_vec.h2, type, w_vec.RDP_gamma, inputScalars.epps, w_vec.GGMRF_p, w_vec.GGMRF_q, w_vec.GGMRF_c);
		if (inputScalars.verbose >= 3)
			mexPrint("OpenMP NLM gradient computed");
		return 0;
	}

	inline int computeRDP(const scalarStruct& inputScalars, const float gamma, const float beta) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting OpenMP RDP gradient computation");
		RDPKernel(d_W, d_inputB, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], inputScalars.NxOrig, inputScalars.NyOrig, inputScalars.NzOrig, gamma, inputScalars.epps, beta);
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
		TVKernel(d_W, d_inputB, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], inputScalars.NxOrig, inputScalars.NyOrig, inputScalars.NzOrig, sigma, smooth, beta, type, data.TV_use_anatomical, C, d_refIm);
		if (inputScalars.verbose >= 3)
			mexPrint("OpenMP TV gradient computed");
		return 0;
	}

	inline int computeGGMRF(const scalarStruct& inputScalars, const float p, const float q, const float c, const float pqc, const float beta) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting OpenMP GGMRF gradient computation");
		GGMRFKernel(d_W, d_inputB, weights, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], p, q, c, pqc, beta, Ndx, Ndy, Ndz);
		if (inputScalars.verbose >= 3)
			mexPrint("OpenMP GGMRF gradient computed");
		return 0;
	}

	inline int hyperGradient(const scalarStruct& inputScalars, const float sigma, const float beta) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting OpenMP GGMRF gradient computation");
		hyperbolicKernel(d_W, d_inputB, weights, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], sigma, beta, Ndx, Ndy, Ndz);
		if (inputScalars.verbose >= 3)
			mexPrint("OpenMP GGMRF gradient computed");
		return 0;
	}

};