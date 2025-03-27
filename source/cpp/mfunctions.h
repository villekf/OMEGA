/*************************************************************************************************************************************************
* MATLAB/Octave related functions.
* This file contains functions to load data from MATLAB/Octave structs.
*
* Copyright (C) 2023-2024 Ville-Veikko Wettenhovi, Niilo Saarlemo
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
#include "structs.h"

inline void loadInput(scalarStruct& inputScalars, const mxArray* options, const int type = -1) {

	if (inputScalars.projector_type == 1 || inputScalars.projector_type == 11 || inputScalars.projector_type == 14 || inputScalars.projector_type == 15 || inputScalars.projector_type == 12 || inputScalars.projector_type == 13)
		inputScalars.FPType = 1;
	else if (inputScalars.projector_type == 2 || inputScalars.projector_type == 21 || inputScalars.projector_type == 24 || inputScalars.projector_type == 22 || inputScalars.projector_type == 23)
		inputScalars.FPType = 2;
	else if (inputScalars.projector_type == 3 || inputScalars.projector_type == 31 || inputScalars.projector_type == 34 || inputScalars.projector_type == 32 || inputScalars.projector_type == 33)
		inputScalars.FPType = 3;
	else if (inputScalars.projector_type == 41 || inputScalars.projector_type == 4 || inputScalars.projector_type == 42 || inputScalars.projector_type == 43 || inputScalars.projector_type == 45)
		inputScalars.FPType = 4;
	else if (inputScalars.projector_type == 51 || inputScalars.projector_type == 5 || inputScalars.projector_type == 54)
		inputScalars.FPType = 5;
	else if (inputScalars.projector_type == 6)
		inputScalars.FPType = 6;

	if (inputScalars.projector_type == 11 || inputScalars.projector_type == 41 || inputScalars.projector_type == 51 || inputScalars.projector_type == 21 || inputScalars.projector_type == 31 || inputScalars.projector_type == 1)
		inputScalars.BPType = 1;
	else if (inputScalars.projector_type == 12 || inputScalars.projector_type == 42 || inputScalars.projector_type == 22 || inputScalars.projector_type == 32 || inputScalars.projector_type == 2)
		inputScalars.BPType = 2;
	else if (inputScalars.projector_type == 13 || inputScalars.projector_type == 43 || inputScalars.projector_type == 23 || inputScalars.projector_type == 33 || inputScalars.projector_type == 3)
		inputScalars.BPType = 3;
	else if (inputScalars.projector_type == 14 || inputScalars.projector_type == 4 || inputScalars.projector_type == 24 || inputScalars.projector_type == 34 || inputScalars.projector_type == 54)
		inputScalars.BPType = 4;
	else if (inputScalars.projector_type == 15 || inputScalars.projector_type == 5 || inputScalars.projector_type == 45)
		inputScalars.BPType = 5;
	else if (inputScalars.projector_type == 6)
		inputScalars.BPType = 6;

	if (type == 1)
		inputScalars.BPType = 0;
	else if (type == 2)
		inputScalars.FPType = 0;

	inputScalars.det_per_ring = getScalarUInt32(options, 0, "det_per_ring");
	inputScalars.rings = getScalarUInt32(options, 0, "rings");
	inputScalars.global_factor = getScalarFloat(options, 0, "global_correction_factor");
	inputScalars.flat = getScalarFloat(options, 0, "flat");
	inputScalars.useMAD = getScalarBool(options, 0, "useMAD");
	inputScalars.useImages = getScalarBool(options, 0, "useImages");
	inputScalars.listmode = getScalarUInt8(options, 0, "listmode");
	inputScalars.indexBased = getScalarBool(options, 0, "useIndexBasedReconstruction");
	inputScalars.relaxScaling = getScalarBool(options, 0, "relaxationScaling");
	inputScalars.computeRelaxation = getScalarBool(options, 0, "computeRelaxationParameters");
	inputScalars.computeSensImag = getScalarBool(options, 0, "compute_sensitivity_image");
	inputScalars.CT = getScalarBool(options, 0, "CT");
	inputScalars.atomic_32bit = getScalarBool(options, 0, "use_32bit_atomics");
	inputScalars.scatter = static_cast<uint32_t>(getScalarBool(options, 0, "additionalCorrection"));
	inputScalars.CTAttenuation = getScalarBool(options, 0, "CT_attenuation");
	inputScalars.largeDim = getScalarBool(options, 0, "largeDim");
	inputScalars.loadTOF = getScalarBool(options, 0, "loadTOF");
	inputScalars.storeResidual = getScalarBool(options, 0, "storeResidual");
	inputScalars.FISTAAcceleration = getScalarBool(options, 0, "FISTA_acceleration");
	inputScalars.stochastic = getScalarBool(options, 0, "stochasticSubsetSelection");
	if (inputScalars.scatter == 1U) {
		inputScalars.size_scat = mxGetNumberOfElements(mxGetCell(getField(options, 0, "ScatterC"), 0));
	}
	inputScalars.PET = getScalarBool(options, 0, "PET");
	inputScalars.CT = getScalarBool(options, 0, "CT");
	inputScalars.SPECT = getScalarBool(options, 0, "SPECT");
	inputScalars.pitch = getScalarBool(options, 0, "pitch");
	inputScalars.enforcePositivity = getScalarBool(options, 0, "enforcePositivity");
	inputScalars.multiResolution = getScalarBool(options, 0, "useMultiResolutionVolumes");
	inputScalars.nMultiVolumes = getScalarUInt32(options, 0, "nMultiVolumes");
	if (inputScalars.FPType == 5 || inputScalars.BPType == 5) {
		inputScalars.meanFP = getScalarBool(options, 0, "meanFP");
		inputScalars.meanBP = getScalarBool(options, 0, "meanBP");
	}
	//if (inputScalars.BPType == 5 || inputScalars.BPType == 4) {
	inputScalars.maskBP = getScalarBool(options, 0, "useMaskBP");
	//}
	inputScalars.maskFP = getScalarBool(options, 0, "useMaskFP");
	if (type == 1)
		inputScalars.maskBP = false;
	else if (type == 2 && inputScalars.FPType > 3)
		inputScalars.maskFP = false;
	inputScalars.offset = getScalarBool(options, 0, "offsetCorrection");
	if (inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
		inputScalars.orthXY = getScalarBool(options, 0, "orthTransaxial");
		inputScalars.orthZ = getScalarBool(options, 0, "orthAxial");
		if (inputScalars.BPType == 3 || inputScalars.FPType == 3)
			inputScalars.cylRadiusProj3 = getScalarFloat(options, 0, "tube_radius");
	}
	if (inputScalars.offset)
		inputScalars.T = getSingles(options, "OffsetLimit");
	//inputScalars.T = getScalarFloat(options, 0, "OffsetLimit");
	inputScalars.nProjections = getScalarInt64(options, 0, "nProjections");
	inputScalars.subsetType = getScalarUInt32(options, 0, "subset_type");
	if (inputScalars.FPType == 4 || inputScalars.FPType == 5 || inputScalars.BPType == 4 || inputScalars.BPType == 5) {
		inputScalars.dL = getScalarFloat(options, 0, "dL");
		inputScalars.d_Scale4.resize(inputScalars.nMultiVolumes + 1);
		inputScalars.dSize.resize(inputScalars.nMultiVolumes + 1);
		inputScalars.d_Scale.resize(inputScalars.nMultiVolumes + 1);
		float* dScaleX4 = getSingles(options, "dScaleX4");
		float* dScaleY4 = getSingles(options, "dScaleY4");
		float* dScaleZ4 = getSingles(options, "dScaleZ4");
		float* dSizeX = nullptr, * dSizeY = nullptr, * dScaleX = nullptr, * dScaleY = nullptr, * dScaleZ = nullptr;
		if (inputScalars.FPType == 5 || inputScalars.BPType == 5) {
			dSizeX = getSingles(options, "dSizeX");
			dSizeY = getSingles(options, "dSizeY");
			//float* dSizeZ = getSingles(options, "dSizeZ");
			dScaleX = getSingles(options, "dScaleX");
			dScaleY = getSingles(options, "dScaleY");
			dScaleZ = getSingles(options, "dScaleZ");
		}
		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
			inputScalars.d_Scale4[ii].x = dScaleX4[ii];
			inputScalars.d_Scale4[ii].y = dScaleY4[ii];
			inputScalars.d_Scale4[ii].z = dScaleZ4[ii];
			//inputScalars.d_Scale.s[3] = 0.f;
			if (inputScalars.FPType == 5 || inputScalars.BPType == 5) {
				inputScalars.dSize[ii].x = dSizeX[ii];
				inputScalars.dSize[ii].y = dSizeY[ii];
				//inputScalars.dSize[ii].z = dSizeZ[ii];
				inputScalars.d_Scale[ii].x = dScaleX[ii];
				inputScalars.d_Scale[ii].y = dScaleY[ii];
				inputScalars.d_Scale[ii].z = dScaleZ[ii];
				if (ii == 0) {
					inputScalars.dSizeBP.x = getScalarFloat(options, 0, "dSizeXBP");
					inputScalars.dSizeBP.y = getScalarFloat(options, 0, "dSizeZBP");
				}
			}
		}
	}
	if (!inputScalars.CT && !inputScalars.SPECT)
		inputScalars.nLayers = getScalarUInt32(options, 0, "nLayers");
	inputScalars.Nf = getScalarUInt32(options, 0, "Nf");
	inputScalars.useExtendedFOV = getScalarBool(options, 0, "useEFOV");
	if (inputScalars.useExtendedFOV)
		inputScalars.eFOV = mxGetNumberOfElements(getField(options, 0, "eFOVIndices")) > 1;
	inputScalars.NxOrig = getScalarUInt32(options, 0, "NxOrig");
	inputScalars.NyOrig = getScalarUInt32(options, 0, "NyOrig");
	inputScalars.NzOrig = getScalarUInt32(options, 0, "NzOrig");
	inputScalars.NxPrior = getScalarUInt32(options, 0, "NxPrior");
	inputScalars.NyPrior = getScalarUInt32(options, 0, "NyPrior");
	inputScalars.NzPrior = getScalarUInt32(options, 0, "NzPrior");
	inputScalars.TGV2D = getScalarBool(options, 0, "use2DTGV");
	inputScalars.adaptiveType = getScalarUInt32(options, 0, "PDAdaptiveType");
	inputScalars.storeFP = getScalarBool(options, 0, "storeFP");
	inputScalars.useTotLength = getScalarBool(options, 0, "useTotLength");
	const uint32_t* devPointer = getUint32s(options, "use_device");
	size_t devLength = mxGetNumberOfElements(mxGetField(options, 0, "use_device"));
	inputScalars.usedDevices = std::vector<uint32_t>(devPointer, devPointer + devLength);
	if (inputScalars.CT) {
		inputScalars.nColsD = getScalarUInt32(getField(options, 0, "nColsD"), -10);
		inputScalars.nRowsD = getScalarUInt32(getField(options, 0, "nRowsD"), -10);
	} else if (inputScalars.SPECT && (inputScalars.projector_type == 1 || inputScalars.projector_type == 2)) {
		inputScalars.nColsD = getScalarUInt32(getField(options, 0, "nColsD"));
		inputScalars.nRowsD = getScalarUInt32(getField(options, 0, "nRowsD"));
        inputScalars.coneOfResponseStdCoeffA = getScalarFloat(getField(options, 0, "coneOfResponseStdCoeffA"));
        inputScalars.coneOfResponseStdCoeffB = getScalarFloat(getField(options, 0, "coneOfResponseStdCoeffB"));
        inputScalars.coneOfResponseStdCoeffC = getScalarFloat(getField(options, 0, "coneOfResponseStdCoeffC"));
	} else {
		inputScalars.nColsD = getScalarUInt32(getField(options, 0, "Nang"));
		inputScalars.nRowsD = getScalarUInt32(getField(options, 0, "Ndist"));
	}
	if (inputScalars.use_psf) {
		inputScalars.g_dim_x = getScalarUInt32(getField(options, 0, "g_dim_x"), -57);
		inputScalars.g_dim_y = getScalarUInt32(getField(options, 0, "g_dim_y"), -58);
		inputScalars.g_dim_z = getScalarUInt32(getField(options, 0, "g_dim_z"), -59);
		inputScalars.deconvolution = getScalarBool(getField(options, 0, "deblurring"), -60);
	}
	if (inputScalars.use_psf && inputScalars.deconvolution) {
		inputScalars.deblur_iterations = getScalarUInt32(getField(options, 0, "deblur_iterations"));
	}

	inputScalars.Nxy = inputScalars.Nx[0] * inputScalars.Ny[0];
	inputScalars.im_dim[0] = static_cast<int64_t>(inputScalars.Nxy) * static_cast<int64_t>(inputScalars.Nz[0]);
	if (inputScalars.multiResolution) {
		for (int ii = 1; ii <= inputScalars.nMultiVolumes; ii++)
			inputScalars.im_dim[ii] = static_cast<int64_t>(inputScalars.Nx[ii]) * static_cast<int64_t>(inputScalars.Ny[ii]) * static_cast<int64_t>(inputScalars.Nz[ii]);
	}
}

// Loads the input data and forms device data variables
inline void form_data_variables(Weighting& w_vec, const mxArray* options, scalarStruct& inputScalars, const RecMethods& MethodList) {
	// Load the number of priors, all MAP-algorithms, non-OS MAP algorithms, non-OS non-MAP algorithms, non-MAP algorithms and total number of algorithms
	uint32_t Ni = 1U;
	// Number of voxels
	if (inputScalars.saveIter)
		Ni = inputScalars.Niter + 1U;
	// Load the necessary variables if the corresponding reconstruction method is used
	int yy = 0;

	if (MethodList.DRAMA) {
		// Relaxation parameter
		w_vec.lambda = getSingles(options, "lam_drama");
	}

	// Load regularization parameter
	w_vec.beta = getScalarFloat(getField(options, 0, "beta"), -9);
	w_vec.betaReg = w_vec.beta;

	// Masks
	if (inputScalars.maskFP) {
		w_vec.maskFP = getUint8s(options, "maskFP");
		inputScalars.maskFPZ = getScalarUInt32(getField(options, 0, "maskFPZ"));
	}
	if (inputScalars.maskBP) {
		w_vec.maskBP = getUint8s(options, "maskBP");
		inputScalars.maskBPZ = getScalarUInt32(getField(options, 0, "maskBPZ"));
		if (DEBUG) {
			const size_t nMask = mxGetNumberOfElements(getField(options, 0, "maskBP"));
			mexPrintBase("nMask = %d\n", nMask);
			mexEval();
		}
	}
	//if (inputScalars.offset)
	//	w_vec.maskOffset = getUint8s(options, "maskOffset");
	if (inputScalars.eFOV && inputScalars.useExtendedFOV && !inputScalars.multiResolution)
		w_vec.eFOVIndices = getUint8s(options, "eFOVIndices");
	if (inputScalars.useExtendedFOV && !inputScalars.multiResolution)
		w_vec.maskPrior = getUint8s(options, "maskPrior");
	else if (inputScalars.maskBP)
		w_vec.maskPrior = getUint8s(options, "maskBP");
	//if (MethodList.NLM || MethodList.MRP || MethodList.CPTV || MethodList.CPTVL1 || MethodList.CPTVKL || MethodList.RDP || (MethodList.TV && !w_vec.data.TV_use_anatomical)
	//	|| MethodList.CPTGV || MethodList.CPTGVKL || MethodList.CPTGVL1)
	// CT-related variables such as number of projection images
	if (inputScalars.CT) {
		w_vec.nProjections = getScalarInt64(getField(options, 0, "nProjections"), -11);
		w_vec.dPitchX = getScalarFloat(getField(options, 0, "dPitchX"));
		w_vec.dPitchY = getScalarFloat(getField(options, 0, "dPitchY"));
	} else if (inputScalars.SPECT) {
		w_vec.nProjections = getScalarInt64(getField(options, 0, "nProjections"));
		w_vec.dPitchX = getScalarFloat(getField(options, 0, "crXY"));
		w_vec.dPitchY = getScalarFloat(getField(options, 0, "crXY"));
		if (inputScalars.projector_type == 1 || inputScalars.projector_type == 2) {
			w_vec.rayShiftsDetector = getSingles(options, "rayShiftsDetector");
			w_vec.rayShiftsSource = getSingles(options, "rayShiftsSource");
		}
	} else {
		w_vec.nProjections = getScalarInt64(getField(options, 0, "nProjections"));
		// Detector pitch
		w_vec.dPitchX = getScalarFloat(getField(options, 0, "cr_p"));
		w_vec.dPitchY = getScalarFloat(getField(options, 0, "cr_pz"));
	}
	if (inputScalars.FPType == 4 || inputScalars.BPType == 4)
		w_vec.kerroin4 = getSingles(options, "kerroin");

#ifdef AF
	if (inputScalars.projector_type == 6U) {
		const mwSize* ind = mxGetDimensions(getField(options, 0, "gFilter"));
		if (DEBUG) {
			mexPrintBase("indX = %d\n", ind[0]);
			mexPrintBase("indY = %d\n", ind[1]);
			mexEval();
		}
		w_vec.angles = getSingles(options, "angles");
		w_vec.gFilter = af::array(ind[0], ind[1], ind[2], ind[3], getSingles(options, "gFilter"));
		w_vec.distInt = getUint32s(options, "blurPlanes");
		if (DEBUG) {
			mexPrintBase("w_vec.gFilter.dims(0) = %d\n", w_vec.gFilter.dims(0));
			mexPrintBase("w_vec.gFilter.dims(1) = %d\n", w_vec.gFilter.dims(1));
			mexPrintBase("w_vec.gFilter.dims(2) = %d\n", w_vec.gFilter.dims(2));
			mexPrintBase("w_vec.gFilter.dims(3) = %d\n", w_vec.gFilter.dims(3));
			mexPrintBase("w_vec.distInt[0] = %d\n", w_vec.distInt[0]);
			mexEval();
		}
		if (DEBUG) {
			mexPrint("SPECT vars loaded");
		}
	}
#endif

	// True value means the preconditioner is included
	// precondTypeIm[0] = Diagonal normalization preconditioner (division with the sensitivity image 1 / (A^T1), A is the system matrix)
	// precondTypeIm[1] = EM preconditioner (f / (A^T1), where f is the current estimate)
	// precondTypeIm[2] = IEM preconditioner (max(n, fhat, f)/ (A^T1), where fhat is an estimate of the final image and n is a small positive number)
	// precondTypeIm[3] = Momentum-like preconditioner (basically a step size inclusion)
	// precondTypeIm[4] = Gradient-based preconditioner (Uses the normalized divergence (sum of the gradient) of the current estimate)
	// precondTypeIm[5] = Filtering-based preconditioner
	// precondTYpeIm[6] = Curvature-based preconditioner
	// The first three are mutually exclusive, but otherwise they can be mixed and matched
	const bool* precondIm = getBools(options, "precondTypeImage");
	w_vec.precondTypeIm.assign(precondIm, precondIm + w_vec.precondTypeIm.size());
	const bool* precondMe = getBools(options, "precondTypeMeas");
	w_vec.precondTypeMeas.assign(precondMe, precondMe + w_vec.precondTypeMeas.size());
	w_vec.filteringOrig = w_vec.precondTypeMeas[1];
	// precondTypeMeas[0] = Diagonal normalization preconditioner (1 / (A1))
	// precondTypeMeas[1] = Filtering-based preconditioner
#ifdef AF
	if (w_vec.precondTypeIm[2]) {
		if (DEBUG) {
			mexPrintBase("w_vec.precondTypeIm[2] = %u\n", w_vec.precondTypeIm[2]);
			mexEval();
		}
		w_vec.preRef.resize(inputScalars.nMultiVolumes + 1);
		const float* ref = getSingles(options, "referenceImage");
		w_vec.preRef[0] = af::array(inputScalars.im_dim[0], ref);
		if (inputScalars.multiResolution) {
			for (int kk = 1; kk <= inputScalars.nMultiVolumes; kk++) {
				size_t dimApu = inputScalars.im_dim[kk - 1];
				//std::string imText = "referenceImage" + std::to_string(kk);
				//w_vec.preRef[kk] = af::array(inputScalars.im_dim[kk], getSingles(options, imText.c_str()));
				w_vec.preRef[kk] = af::array(inputScalars.im_dim[kk], &ref[dimApu]);
			}
		}
		if (DEBUG) {
			mexPrint("Precond im [2] loaded");
		}
	}

	if (w_vec.precondTypeIm[4]) {
		if (DEBUG) {
			mexPrintBase("w_vec.precondTypeIm[4] = %u\n", w_vec.precondTypeIm[4]);
			mexEval();
		}
		w_vec.gradV1 = getScalarFloat(getField(options, 0, "gradV1"));
		w_vec.gradV2 = getScalarFloat(getField(options, 0, "gradV2"));
		w_vec.gradInitIter = getScalarUInt32(getField(options, 0, "gradInitIter"));
		w_vec.gradFinalIter = getScalarUInt32(getField(options, 0, "gradLastIter"));
		w_vec.gradF.resize(inputScalars.nMultiVolumes + 1);
		if (DEBUG) {
			mexPrint("Precond im [4] loaded");
		}
	}

	if (w_vec.precondTypeIm[5]) {
		if (DEBUG) {
			mexPrintBase("w_vec.precondTypeIm[5] = %u\n", w_vec.precondTypeIm[5]);
			mexEval();
		}
		w_vec.filterIm = af::array(inputScalars.Nf, inputScalars.Nf, getSingles(options, "filterIm"));
		w_vec.filterIter = getScalarUInt32(getField(options, 0, "filteringIterations"));
		if (DEBUG) {
			mexPrint("Precond im [5] loaded");
		}
	}

	if (w_vec.precondTypeMeas[1]) {
		if (DEBUG) {
			mexPrintBase("w_vec.precondTypeMeas[1] = %u\n", w_vec.precondTypeMeas[1]);
			mexEval();
		}
		if (FINVERSE)
			w_vec.filter = af::array(inputScalars.Nf, getSingles(options, "filter"));
		else
			if (inputScalars.subsetType == 5)
				w_vec.filter = af::array(inputScalars.nColsD, getSingles(options, "filter2"));
			else
				w_vec.filter = af::array(inputScalars.nRowsD, getSingles(options, "filter2"));
		w_vec.Ffilter = af::array(inputScalars.Nf, getSingles(options, "Ffilter"));
		w_vec.filterIter = getScalarUInt32(getField(options, 0, "filteringIterations"));
		if (DEBUG) {
			mexPrint("Precond meas [1] loaded");
		}
	}

	// The complete sensitivity image is computed
	if (MethodList.RBI || MethodList.RBIOSL || MethodList.COSEM || MethodList.ACOSEM || MethodList.ECOSEM || (w_vec.precondTypeIm[0]
		|| w_vec.precondTypeIm[1] || w_vec.precondTypeIm[2])) {
		if (DEBUG) {
			mexPrintBase("w_vec.precondTypeIm[1] = %u\n", w_vec.precondTypeIm[1]);
			mexEval();
		}
		w_vec.computeD = true;
	}

	if (w_vec.precondTypeMeas[0] || MethodList.SART || MethodList.POCS)
		w_vec.computeM = true;
#endif

	// Load TV related input data
	if (MethodList.TV && MethodList.MAP) {
		// Is anatomical reference image used
		w_vec.data.TV_use_anatomical = getScalarBool(getField(options, 0, "TV_use_anatomical"), -13);
		// Tau-value
		w_vec.data.tau = getScalarFloat(getField(options, 0, "tau"), -14);
		// "Smoothing" parameter, prevents zero values in the square root
		w_vec.data.TVsmoothing = getScalarFloat(getField(options, 0, "TVsmoothing"), -15);
		// The type of TV prior used
		w_vec.data.TVtype = getScalarUInt32(getField(options, 0, "TVtype"), -16);
		// If anatomical prior is used, load the necessary coefficients
#ifdef AF
		if (w_vec.data.TV_use_anatomical) {
			mxArray* TVdata_init = getField(options, 0, "TVdata");
			if (w_vec.data.TVtype == 1) {
				w_vec.data.refIm = af::array(inputScalars.im_dim[0], getSingles(TVdata_init, "s", 0), afHost);
				//w_vec.data.s1 = af::array(inputScalars.im_dim[0], getSingles(TVdata_init, "s1", 0), afHost);
				//w_vec.data.s2 = af::array(inputScalars.im_dim[0], getSingles(TVdata_init, "s2", 0), afHost);
				//w_vec.data.s3 = af::array(inputScalars.im_dim[0], getSingles(TVdata_init, "s3", 0), afHost);
				//w_vec.data.s4 = af::array(inputScalars.im_dim[0], getSingles(TVdata_init, "s4", 0), afHost);
				//w_vec.data.s5 = af::array(inputScalars.im_dim[0], getSingles(TVdata_init, "s5", 0), afHost);
				//w_vec.data.s6 = af::array(inputScalars.im_dim[0], getSingles(TVdata_init, "s6", 0), afHost);
				//w_vec.data.s7 = af::array(inputScalars.im_dim[0], getSingles(TVdata_init, "s7", 0), afHost);
				//w_vec.data.s8 = af::array(inputScalars.im_dim[0], getSingles(TVdata_init, "s8", 0), afHost);
				//w_vec.data.s9 = af::array(inputScalars.im_dim[0], getSingles(TVdata_init, "s9", 0), afHost); 
			}
			else {
				w_vec.data.refIm = af::array(inputScalars.im_dim[0], getSingles(TVdata_init, "reference_image", 0), afHost);
			}
			w_vec.data.T = getScalarFloat(getField(options, 0, "T"), -17);
			w_vec.data.C = getScalarFloat(getField(options, 0, "C"), -18);
		}
		// Additional weights for the TV type 3
		//if (w_vec.data.TVtype == 3 && !MethodList.Quad && w_vec.data.TV_use_anatomical) {
		//	w_vec.Ndx = getScalarUInt32(getField(options, 0, "Ndx"), -19);
		//	w_vec.Ndy = getScalarUInt32(getField(options, 0, "Ndy"), -20);
		//	w_vec.Ndz = getScalarUInt32(getField(options, 0, "Ndz"), -21);
		//	w_vec.dimmu = (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1);
		//	w_vec.weights_TV = af::array(w_vec.dimmu - 1, getSingles(options, "weights_quad", 0), afHost);
		//}
#endif
		if (w_vec.data.TVtype == 4 || w_vec.data.TVtype == 6)
			w_vec.data.SATVPhi = getScalarFloat(getField(options, 0, "SATVPhi"), -23);
		if (DEBUG) {
			mexPrintBase("w_vec.data.TVtype = %u\n", w_vec.data.TVtype);
			mexPrintBase("w_vec.data.SATVPhi = %f\n", w_vec.data.SATVPhi);
			mexEval();
		}
	}
	// General variables for neighborhood-based methods
	if ((MethodList.L || MethodList.FMH || MethodList.WeightedMean || MethodList.Quad || MethodList.Huber || MethodList.MRP || MethodList.NLM || MethodList.ProxNLM || MethodList.hyperbolic || MethodList.RDP) 
		&& MethodList.MAP) {
		// Neighborhood size
		w_vec.Ndx = getScalarUInt32(getField(options, 0, "Ndx"), -24);
		w_vec.Ndy = getScalarUInt32(getField(options, 0, "Ndy"), -25);
		w_vec.Ndz = getScalarUInt32(getField(options, 0, "Ndz"), -26);
		// Is normalization used in MRP, FMH, L, weighted mean or AD
		w_vec.med_no_norm = getScalarBool(getField(options, 0, "med_no_norm"), -27);
		w_vec.dimmu = (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1);
		if (DEBUG) {
			mexPrint("Neighborhood loaded");
		}
	}
#ifdef AF
	if ((MethodList.L || MethodList.FMH) && MethodList.MAP) {
		// Index values for the neighborhood
		w_vec.tr_offsets = af::array(inputScalars.im_dim[0], w_vec.dimmu, getUint32s(options, "tr_offsets", 0), afHost);
	}
#endif
	if (MethodList.FMH || MethodList.Quad || MethodList.Huber)
		w_vec.inffi = getScalarUInt32(getField(options, 0, "inffi"), -28);
	// Weights for the quadratic prior
#ifdef AF
	if (MethodList.Quad && MethodList.MAP) {
		w_vec.weights_quad = af::array(w_vec.dimmu - 1, getSingles(options, "weights_quad", 0), afHost);
		int pituus = (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1);
		af::array weight = w_vec.weights_quad * -1.f;
		af::array w_quad = af::constant(0.f, pituus);
		w_quad(af::seq(0, pituus / 2 - 1)) = weight(af::seq(0, pituus / 2 - 1));
		w_quad(af::seq(pituus / 2 + 1, af::end)) = weight(af::seq(pituus / 2, af::end));
		w_quad(pituus / 2) = af::abs(af::sum(weight));
		w_vec.weights_quad = af::moddims(w_quad, w_vec.Ndx * 2U + 1, w_vec.Ndy * 2U + 1, w_vec.Ndz * 2U + 1);
		if (DEBUG) {
			mexPrint("Quad loaded");
		}
	}
	if (MethodList.hyperbolic && MethodList.MAP) {
		w_vec.weights = getSingles(options, "weights_quad", 0);
		w_vec.data.SATVPhi = getScalarFloat(getField(options, 0, "hyperbolicDelta"), -23);
	}
#endif
	if ((MethodList.RDP && MethodList.MAP) || MethodList.ProxRDP) {
		w_vec.RDP_gamma = getScalarFloat(getField(options, 0, "RDP_gamma"), -29);
		w_vec.RDPLargeNeighbor = getScalarBool(getField(options, 0, "RDPIncludeCorners"), -29);
		w_vec.RDP_anatomical = getScalarBool(getField(options, 0, "RDP_use_anatomical"), -29);
		if (w_vec.RDP_anatomical)
			w_vec.RDP_ref = getSingles(options, "RDP_ref", 0);
		if (w_vec.RDPLargeNeighbor)
			w_vec.weights = getSingles(options, "weights_quad");
		if (DEBUG) {
			mexPrint("RDP loaded");
			mexPrintBase("w_vec.RDPLargeNeighbor = %d\n", w_vec.RDPLargeNeighbor);
		}
	}
	if (MethodList.GGMRF) {
		w_vec.GGMRF_p = getScalarFloat(getField(options, 0, "GGMRF_p"), -29);
		w_vec.GGMRF_q = getScalarFloat(getField(options, 0, "GGMRF_q"), -29);
		w_vec.GGMRF_c = getScalarFloat(getField(options, 0, "GGMRF_c"), -29);
		w_vec.GGMRF_pqc = (w_vec.GGMRF_p - w_vec.GGMRF_q) / std::pow(w_vec.GGMRF_c, w_vec.GGMRF_p - w_vec.GGMRF_q);
		w_vec.weights = getSingles(options, "weights_quad");
		if (DEBUG) {
			mexPrint("GGMRF loaded");
		}
	}
#ifdef AF
	if (MethodList.Huber && MethodList.MAP) {
		w_vec.weights_huber = af::array(w_vec.dimmu - 1, getSingles(options, "weights_huber", 0), afHost);
		int pituus = (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1);
		af::array weight = w_vec.weights_huber * -1.f;
		af::array w_quad = af::constant(0.f, pituus);
		w_quad(af::seq(0, pituus / 2 - 1)) = weight(af::seq(0, pituus / 2 - 1));
		w_quad(af::seq(pituus / 2 + 1, af::end)) = weight(af::seq(pituus / 2, af::end));
		w_quad(pituus / 2) = af::abs(af::sum(weight));
		w_vec.weights_huber = af::moddims(w_quad, w_vec.Ndx * 2U + 1, w_vec.Ndy * 2U + 1, w_vec.Ndz * 2U + 1);
		w_vec.huber_delta = getScalarFloat(getField(options, 0, "huber_delta"), -29);
		if (DEBUG) {
			mexPrint("Huber loaded");
		}
	}
	if (MethodList.L && MethodList.MAP)
		w_vec.a_L = af::array(w_vec.dimmu, getSingles(options, "a_L", 0), afHost);
	if (MethodList.FMH && MethodList.MAP) {
		if (inputScalars.Nz[0] == 1 || w_vec.Ndz == 0)
			w_vec.fmh_weights = af::array(w_vec.Ndx * 2 + 1, 4, getSingles(options, "fmh_weights", 0), afHost);
		else
			w_vec.fmh_weights = af::array((std::max)(w_vec.Ndz * 2 + 1, w_vec.Ndx * 2 + 1), 13, getSingles(options, "fmh_weights", 0), afHost);
		w_vec.alku_fmh = getScalarUInt32(getField(options, 0, "inffi"), -30);
	}
	if (MethodList.WeightedMean && MethodList.MAP) {
		w_vec.weighted_weights = af::moddims(af::array(w_vec.dimmu, getSingles(options, "weighted_weights", 0), afHost), w_vec.Ndx * 2U + 1U, w_vec.Ndy * 2U + 1U, w_vec.Ndz * 2U + 1U);
		// Type of mean used (arithmetic, harmonic or geometric)
		w_vec.mean_type = getScalarInt32(getField(options, 0, "mean_type"), -31);
		// Sum of the weights
		w_vec.w_sum = getScalarFloat(getField(options, 0, "w_sum"), -31);
		if (DEBUG) {
			mexPrint("WMean loaded");
		}
	}
	if (MethodList.AD && MethodList.MAP) {
		// Time-step value
		w_vec.TimeStepAD = getScalarFloat(getField(options, 0, "TimeStepAD"), -32);
		// Conductance (edge value)
		w_vec.KAD = getScalarFloat(getField(options, 0, "KAD"), -33);
		// Number of AD iterations
		w_vec.NiterAD = getScalarUInt32(getField(options, 0, "NiterAD"), -34);
		// Flux type
		uint32_t Flux = getScalarUInt32(getField(options, 0, "FluxType"), -35);
		// Diffusion type
		uint32_t Diffusion = getScalarUInt32(getField(options, 0, "DiffusionType"), -36);
		if (Flux == 2U)
			w_vec.FluxType = AF_FLUX_QUADRATIC;
		else
			w_vec.FluxType = AF_FLUX_EXPONENTIAL;
		if (Diffusion == 2U)
			w_vec.DiffusionType = AF_DIFFUSION_MCDE;
		else
			w_vec.DiffusionType = AF_DIFFUSION_GRAD;
		if (!MethodList.L && !MethodList.FMH && !MethodList.WeightedMean && !MethodList.Quad && !MethodList.MRP)
			w_vec.med_no_norm = getScalarBool(getField(options, 0, "med_no_norm"), -37);
	}
	// Asymmetric parallel level sets
	if (MethodList.APLS && MethodList.MAP) {
		// Eta value
		w_vec.data.eta = getScalarFloat(getField(options, 0, "eta"), -38);
		// Tau-value
		if (!MethodList.TV)
			w_vec.data.tau = getScalarFloat(getField(options, 0, "tau"), -39);
		// Smoothing value
		w_vec.data.TVsmoothing = getScalarFloat(getField(options, 0, "APLSsmoothing"), -40);
		// Anatomical reference image
		w_vec.data.refIm = af::array(inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], getSingles(options, "APLS_ref_image", 0), afHost);
		w_vec.data.TV_use_anatomical = true;
		w_vec.data.TVtype = 5;
	}
#endif
	if (MethodList.MRAMLA || MethodList.MBSREM || MethodList.SPS) {
		// Relaxation parameter
		w_vec.lambda = getSingles(options, "lambda");
		// Upper bound
		w_vec.U = getScalarFloat(getField(options, 0, "U"), -42);
	}
	// Relaxation parameters
	if (MethodList.RAMLA || MethodList.BSREM || MethodList.ROSEM || MethodList.ROSEMMAP || MethodList.SART || MethodList.POCS || MethodList.SAGA)
		w_vec.lambda = getSingles(options, "lambda");
	if (MethodList.PKMA) {
		w_vec.alphaM = getSingles(options, "alpha_PKMA");
		w_vec.lambda = getSingles(options, "lambda");
	}
	if ((w_vec.precondTypeIm[5] || w_vec.precondTypeMeas[1]) && (MethodList.MRAMLA || MethodList.MBSREM || MethodList.SPS || MethodList.RAMLA || MethodList.BSREM || MethodList.ROSEM || MethodList.ROSEMMAP || MethodList.PKMA || MethodList.SAGA)) {
		w_vec.lambdaFiltered = w_vec.lambda;
		w_vec.lambda = getSingles(options, "lambdaFiltered");
	}
	if (DEBUG && (MethodList.MRAMLA || MethodList.MBSREM || MethodList.SPS || MethodList.RAMLA || MethodList.BSREM || MethodList.ROSEM || MethodList.ROSEMMAP || MethodList.PKMA || MethodList.SAGA)) {
		mexPrintBase("w_vec.lambda[0] = %f\n", w_vec.lambda[0]);
		mexEval();
	}
	if (w_vec.precondTypeIm[3])
		w_vec.alphaPrecond = getSingles(options, "alphaPrecond");
	// Power factor for ACOSEM
	if (MethodList.ACOSEM || MethodList.OSLCOSEM == 1) {
		w_vec.h_ACOSEM = getScalarFloat(getField(options, 0, "h"), -43);
		w_vec.h_ACOSEM_2 = 1.f / w_vec.h_ACOSEM;
	}
	//if (MethodList.TGV && MethodList.MAP) {
	//	w_vec.data.TGVAlpha = getScalarFloat(getField(options, 0, "alpha1TGV"), -44);
	//	w_vec.data.TGVBeta = getScalarFloat(getField(options, 0, "beta"), -45);
	//	w_vec.data.NiterTGV = getScalarUInt32(getField(options, 0, "NiterTGV"), -46);
	//}
	if ((MethodList.NLM && MethodList.MAP) || MethodList.ProxNLM) {
		w_vec.NLM_anatomical = getScalarBool(getField(options, 0, "NLM_use_anatomical"), -47);
		w_vec.NLTV = getScalarBool(getField(options, 0, "NLTV"), -48);
		w_vec.NLRD = getScalarBool(getField(options, 0, "NLRD"), -48);
		w_vec.NLM_MRP = getScalarBool(getField(options, 0, "NLM_MRP"), -49);
		w_vec.NLLange = getScalarBool(getField(options, 0, "NLLange"), -49);
		w_vec.NLGGMRF = getScalarBool(getField(options, 0, "NLGGMRF"), -49);
		w_vec.NLAdaptive = getScalarBool(getField(options, 0, "NLAdaptive"), -49);
		if (w_vec.NLRD)
			w_vec.RDP_gamma = getScalarFloat(getField(options, 0, "RDP_gamma"), -29);
		else if (w_vec.NLLange)
			w_vec.RDP_gamma = getScalarFloat(getField(options, 0, "SATVPhi"), -23);
		else if (w_vec.NLGGMRF) {
			w_vec.GGMRF_p = getScalarFloat(getField(options, 0, "GGMRF_p"), -29);
			w_vec.GGMRF_q = getScalarFloat(getField(options, 0, "GGMRF_q"), -29);
			w_vec.GGMRF_c = getScalarFloat(getField(options, 0, "GGMRF_c"), -29);
			w_vec.RDP_gamma = (w_vec.GGMRF_p - w_vec.GGMRF_q) / std::pow(w_vec.GGMRF_c, w_vec.GGMRF_p - w_vec.GGMRF_q);
		}
		if (w_vec.NLM_anatomical)
			w_vec.NLM_ref = getSingles(options, "NLM_ref", 0);
		if (w_vec.NLAdaptive)
			w_vec.NLAdaptiveConstant = getScalarFloat(getField(options, 0, "NLAdaptiveConstant"), -29);
		w_vec.h2 = getScalarFloat(getField(options, 0, "sigma"), -50);
		w_vec.h2 = w_vec.h2 * w_vec.h2;
		w_vec.Nlx = getScalarUInt32(getField(options, 0, "Nlx"), -51);
		w_vec.Nly = getScalarUInt32(getField(options, 0, "Nly"), -52);
		w_vec.Nlz = getScalarUInt32(getField(options, 0, "Nlz"), -53);
#ifdef AF
		w_vec.gaussianNLM = af::array((2 * w_vec.Nlx + 1) * (2 * w_vec.Nly + 1) * (2 * w_vec.Nlz + 1), getSingles(options, "gaussianNLM", 0), afHost);
#endif
		if (DEBUG) {
			mexPrint("NLM loaded");
		}
	}
	if (MethodList.CPType || MethodList.FISTA || MethodList.FISTAL1) {
		w_vec.tauCP = getSingles(options, "tauCPFilt");
		w_vec.sigmaCP = getSingles(options, "sigmaCP");
		w_vec.powerIterations = getScalarUInt32(getField(options, 0, "powerIterations"), -63);
		if (DEBUG) {
			mexPrint("PIter loaded");
		}
	}
	if (MethodList.CPType || MethodList.FISTA || MethodList.FISTAL1 || MethodList.ProxTGV || MethodList.ProxTV) {
		if (w_vec.precondTypeMeas[1])
			w_vec.tauCP2 = getSingles(options, "tauCP");
		else
			w_vec.tauCP = getSingles(options, "tauCP");
		w_vec.sigma2CP = getSingles(options, "sigma2CP");
		w_vec.betaReg = getScalarFloat(getField(options, 0, "beta"), -63);
		w_vec.thetaCP = getSingles(options, "thetaCP", 0);
		w_vec.alpha0CPTGV = getScalarFloat(getField(options, 0, "alpha0TGV"), -63);
		w_vec.alpha1CPTGV = getScalarFloat(getField(options, 0, "alpha1TGV"), -63);
		w_vec.UseL2Ball = getScalarBool(getField(options, 0, "useL2Ball"), -63);
		inputScalars.FISTAType = getScalarUInt32(getField(options, 0, "FISTAType"), -63);
		if (inputScalars.adaptiveType == 1) {
			for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++)
				//w_vec.alphaCP.emplace_back(.3f);
				w_vec.alphaCP.emplace_back(1.f);
		}
		else if (inputScalars.adaptiveType == 2)
			for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++)
				w_vec.alphaCP.emplace_back(.95f);
		if (DEBUG) {
			mexPrint("CPType loaded");
		}
	}
	inputScalars.useFDKWeights = getScalarBool(getField(options, 0, "useFDKWeights"), -49);
	if (MethodList.FDK && inputScalars.CT && inputScalars.useFDKWeights) {
		w_vec.angles = getSingles(options, "angles");
		inputScalars.DSC = getScalarFloat(getField(options, 0, "sourceToCRot"), -63);
	}
	w_vec.derivType = getScalarUInt32(getField(options, 0, "derivativeType"), -63);
	if (MethodList.POCS) {
		w_vec.ng = getScalarUInt32(getField(options, 0, "POCS_NgradIter"), -63);
		w_vec.alphaPOCS = getScalarFloat(getField(options, 0, "POCS_alpha"), -63);
		w_vec.rMaxPOCS = getScalarFloat(getField(options, 0, "POCS_rMax"), -63);
		w_vec.POCSalphaRed = getScalarFloat(getField(options, 0, "POCS_alphaRed"), -63);
		w_vec.POCSepps = getScalarFloat(getField(options, 0, "POCSepps"), -63);
	}
}

// Obtain the reconstruction methods used
inline void get_rec_methods(const mxArray* options, RecMethods& MethodList) {
	// Non-MAP/prior or single prior algorithms
	MethodList.OSEM = getScalarBool(getField(options, 0, "OSEM"), -61);
	MethodList.RAMLA = getScalarBool(getField(options, 0, "RAMLA"), -61);
	MethodList.MRAMLA = getScalarBool(getField(options, 0, "MRAMLA"), -61);
	MethodList.ROSEM = getScalarBool(getField(options, 0, "ROSEM"), -61);
	MethodList.RBI = getScalarBool(getField(options, 0, "RBI"), -61);
	MethodList.DRAMA = getScalarBool(getField(options, 0, "DRAMA"), -61);
	MethodList.COSEM = getScalarBool(getField(options, 0, "COSEM"), -61);
	MethodList.ECOSEM = getScalarBool(getField(options, 0, "ECOSEM"), -61);
	MethodList.ACOSEM = getScalarBool(getField(options, 0, "ACOSEM"), -61);
	MethodList.LSQR = getScalarBool(getField(options, 0, "LSQR"), -61);
	MethodList.CGLS = getScalarBool(getField(options, 0, "CGLS"), -61);
	MethodList.SART = getScalarBool(getField(options, 0, "SART"), -61);
	MethodList.FISTA = getScalarBool(getField(options, 0, "FISTA"), -61);
	MethodList.FISTAL1 = getScalarBool(getField(options, 0, "FISTAL1"), -61);
	if (MethodList.LSQR || MethodList.CGLS)
		MethodList.initAlg = true;

	// Priors
	MethodList.MRP = getScalarBool(getField(options, 0, "MRP"), -61);
	MethodList.Quad = getScalarBool(getField(options, 0, "quad"), -61);
	MethodList.Huber = getScalarBool(getField(options, 0, "Huber"), -61);
	MethodList.L = getScalarBool(getField(options, 0, "L"), -61);
	MethodList.FMH = getScalarBool(getField(options, 0, "FMH"), -61);
	MethodList.WeightedMean = getScalarBool(getField(options, 0, "weighted_mean"), -61);
	MethodList.TV = getScalarBool(getField(options, 0, "TV"), -61);
	MethodList.AD = getScalarBool(getField(options, 0, "AD"), -61);
	MethodList.APLS = getScalarBool(getField(options, 0, "APLS"), -61);
	MethodList.TGV = getScalarBool(getField(options, 0, "TGV"), -61);
	MethodList.NLM = getScalarBool(getField(options, 0, "NLM"), -61);
	MethodList.hyperbolic = getScalarBool(getField(options, 0, "hyperbolic"), -61);
	MethodList.RDP = getScalarBool(getField(options, 0, "RDP"), -61);
	MethodList.GGMRF = getScalarBool(getField(options, 0, "GGMRF"), -61);
	MethodList.ProxTV = getScalarBool(getField(options, 0, "ProxTV"), -61);
	MethodList.ProxTGV = getScalarBool(getField(options, 0, "TGV"), -61);
	MethodList.ProxRDP = getScalarBool(getField(options, 0, "ProxRDP"), -61);
	MethodList.ProxNLM = getScalarBool(getField(options, 0, "ProxNLM"), -61);

	// MAP/prior-based algorithms
	MethodList.OSLOSEM = getScalarBool(getField(options, 0, "OSL_OSEM"), -61);
	MethodList.BSREM = getScalarBool(getField(options, 0, "BSREM"), -61);
	MethodList.MBSREM = getScalarBool(getField(options, 0, "MBSREM"), -61);
	MethodList.ROSEMMAP = getScalarBool(getField(options, 0, "ROSEM_MAP"), -61);
	MethodList.RBIOSL = getScalarBool(getField(options, 0, "OSL_RBI"), -61);
	MethodList.OSLCOSEM = getScalarUInt32(getField(options, 0, "OSL_COSEM"), -61);
	MethodList.PKMA = getScalarBool(getField(options, 0, "PKMA"), -61);
	MethodList.SPS = getScalarBool(getField(options, 0, "SPS"), -61);
	MethodList.PDHG = getScalarBool(getField(options, 0, "PDHG"), -61);
	MethodList.PDHGKL = getScalarBool(getField(options, 0, "PDHGKL"), -61);
	MethodList.PDHGL1 = getScalarBool(getField(options, 0, "PDHGL1"), -61);
	MethodList.CV = getScalarBool(getField(options, 0, "CV"), -61);
	MethodList.PDDY = getScalarBool(getField(options, 0, "PDDY"), -61);
	MethodList.POCS = getScalarBool(getField(options, 0, "ASD_POCS"), -61);
	MethodList.SAGA = getScalarBool(getField(options, 0, "SAGA"), -61);

	// Whether MAP/prior-based algorithms are used
	MethodList.MAP = getScalarBool(getField(options, 0, "MAP"), -61);

	// Custom prior
	MethodList.CUSTOM = getScalarBool(getField(options, 0, "custom"), -61);

	MethodList.FDK = getScalarBool(getField(options, 0, "FDK"), -61);

	// Primal-dual algorithms
	if (MethodList.PDHG || MethodList.PDHGKL || MethodList.PDHGL1 || MethodList.CV || MethodList.PDDY)
		MethodList.CPType = true;
}

#ifdef AF
// Transfers the device data to host
// First transfer the ArrayFire arrays from the device to the host pointers pointing to the mxArrays
// Transfer the mxArrays to the cell
inline void device_to_host(const RecMethods& MethodList, AF_im_vectors& vec, int64_t& oo, mxArray* cell, mxArray* FPcell, Weighting& w_vec,
	const uint32_t dim_n, const scalarStruct& inputScalars, std::vector<std::vector<std::vector<float>>>& FPEstimates) {
	if (DEBUG) {
		mexPrintBase("vec.im_os.dims(0) = %d\n", vec.im_os[0].dims(0));
		mexPrintBase("vec.im_os.dims(1) = %d\n", vec.im_os[0].dims(1));
		mexEval();
	}
	if (inputScalars.storeFP) {
		for (uint32_t ii = 0; ii < inputScalars.subsets * inputScalars.Niter; ii++) {
			const uint32_t jj = ii % inputScalars.subsets;
			const uint32_t kk = ii / inputScalars.subsets;
			const mwSize dim[1] = { static_cast<mwSize>(FPEstimates[kk][jj].size()) };
			mxArray* apu = mxCreateNumericArray(1, dim, mxSINGLE_CLASS, mxREAL);
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			float* apuF = (float*)mxGetSingles(apu);
#else
			float* apuF = (float*)mxGetData(apu);
#endif
			std::copy(FPEstimates[kk][jj].begin(), FPEstimates[kk][jj].end(), apuF);
			mxSetCell(FPcell, static_cast<mwIndex>(ii), mxDuplicateArray(apu));
		}
	}
	// Transfer data back to host
	if (CELL) {
		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
			const mwSize dim[3] = { static_cast<mwSize>(inputScalars.Nx[ii]), static_cast<mwSize>(inputScalars.Ny[ii]), static_cast<mwSize>(inputScalars.Nz[ii]) };
			if (DEBUG) {
				mexPrintBase("inputScalars.Nx[ii] = %d\n", inputScalars.Nx[ii]);
				mexPrintBase("inputScalars.Ny[ii] = %d\n", inputScalars.Ny[ii]);
				mexPrintBase("inputScalars.Nz[ii] = %d\n", inputScalars.Nz[ii]);
				mexEval();
			}
			mxArray* apu = mxCreateNumericArray(3, dim, mxSINGLE_CLASS, mxREAL);
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && TARGET_API_VERSION > 700
			float* apuF = (float*)mxGetSingles(apu);
#else
			float* apuF = (float*)mxGetData(apu);
#endif
			if (inputScalars.saveIter || inputScalars.saveIterationsMiddle > 0) {
			}
			else {
				if (MethodList.FDK)
					vec.rhs_os[ii].host(&apuF[oo]);
				else
					vec.im_os[ii].host(&apuF[oo]);
				if (inputScalars.verbose >= 3)
					mexPrint("Data transfered to host");
			}
			mxSetCell(cell, static_cast<mwIndex>(ii), mxDuplicateArray(apu));
		}
	}
	else {
		float* apuF = getSingles(cell, "solu");
		if (inputScalars.saveIter || inputScalars.saveIterationsMiddle > 0) {
		}
		else {
			if (MethodList.FDK && inputScalars.largeDim) {
				//vec.rhs_os[0].host(&output[oo]);
			}
			else if (MethodList.FDK && !inputScalars.largeDim) {
				vec.rhs_os[0].host(&apuF[oo]);
			}
			else
				vec.im_os[0].host(&apuF[oo]);
			if (inputScalars.verbose >= 3)
				mexPrint("Data transfered to host");
		}
	}
	af::sync();
}
#endif