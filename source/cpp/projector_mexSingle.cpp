/*******************************************************************************************************************************************
* Implements either the improved Siddon's algorithm, orthogonal distance based projector or volume-based projector for OMEGA.
*
* Implementation 1 supports only Siddon's algorithm and outputs the sparse system matrix. Implementation 4 supports projectors 1-3 and is 
* computed matrix-free.
* 
* This file contains the mex-functions for both the implementation type 1 and type 4.
* 
* A precomputed vector that has the number of voxels a LOR/ray passes, as well as the number of voxels that have been traversed in previous 
* LORs/rays is required if options.precompute_lor = true, but is optional otherwise. If the precompute_lor options is set to false and 
* implementation type 4 has been selected then all the measurements are investigated, but only those intercepting the FOV will be included.
* 
* Both implementations use OpenMP for parallellization (if available, otherwise the code will be sequential with no parallelization).
* 
* Copyright (C) 2022 Ville-Veikko Wettenhovi
*
* This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************************************************************************/
#define DEBUG false
#define CAST float
#include "projector_functions.h"
#include "mex.h"
#include "mexFunktio.h"
 
using namespace std;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	// Check for the number of input and output arguments
	if (nrhs < 39)
		mexErrMsgTxt("Too few input arguments. There must be at least 39.");
	else if (nrhs > 44)
		mexErrMsgTxt("Too many input arguments. There can be at most 44.");

	if (nlhs > 3 || nlhs < 1)
		mexErrMsgTxt("Invalid number of output arguments. There can be at most three.");

	int ind = 0;
	const mxArray* options = prhs[ind];
	ind++;
	paramStruct<CAST> param;
	// Load the input arguments
	// Image size in x-direction
	param.Nx = getScalarUInt32(prhs[ind], ind);
	ind++;
	
	if (DEBUG)
		mexPrintf("param.Nx = %u\n", param.Nx);

	// Image size in y-direction
	param.Ny = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Image size in z-direction
	param.Nz = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Distance between adjacent pixels in x-direction
	param.dx = getScalarFloat(prhs[ind], ind);
	ind++;

	// Distance between adjacent pixels in y-direction
	param.dy = getScalarFloat(prhs[ind], ind);
	ind++;

	// Distance between adjacent pixels in z-direction
	param.dz = getScalarFloat(prhs[ind], ind);
	ind++;

	// Distance of the pixel grid from the origin in x-direction
	param.bx = getScalarFloat(prhs[ind], ind);
	ind++;

	// Distance of the pixel grid from the origin in y-direction
	param.by = getScalarFloat(prhs[ind], ind);
	ind++;

	// Distance of the pixel grid from the origin in z-direction
	param.bz = getScalarFloat(prhs[ind], ind);
	ind++;

	// Coordinates of the detectors in z-direction
	const float* z = getSingles(prhs[ind], "solu");
	ind++;

	// Coordinates of the detectors in x/y-direction
	float* x = getSingles(prhs[ind], "solu");
	ind++;

	// Number of detector indices
	param.size_x = getScalarUInt32(prhs[ind], ind);
	ind++;

	// attenuation values (attenuation images)
	param.atten = getSingles(prhs[ind], "solu");
	ind++;

	// Normalization coefficients
	param.normCoefs = getSingles(prhs[ind], "solu");
	ind++;

	// Number of measurements/LORs
	const int64_t pituus = getScalarInt64(prhs[ind], ind);
	ind++;

	// Is the attenuation correction included
	param.attenuationCorrection = getScalarBool(prhs[ind], ind);
	ind++;

	// Is the normalization correction included
	param.normalizationCorrection = getScalarBool(prhs[ind], ind);
	ind++;

	// Is scatter correction included as multiplication (system matrix)
	param.scatterCorrectionMult = getScalarBool(prhs[ind], ind);
	ind++;

	// Scatter data
	param.scatterCoefs = getSingles(prhs[ind], "solu");
	ind++;

	// Global correction factor
	param.globalFactor = getScalarFloat(prhs[ind], ind);
	ind++;

	// For sinogram data, the indices of the detectors corresponding to the current sinogram bin
	param.xy_index = getUint32s(prhs[ind], "solu");
	ind++;

	// Same as above, but for z-direction
	param.z_index = getUint16s(prhs[ind], "solu");
	ind++;

	// Detector pair numbers, for raw list-mode data
	const uint16_t* detIndices = getUint16s(prhs[ind], "solu");
	const size_t numRows = mxGetNumberOfElements(prhs[ind]);
	ind++;

	// Number of detectors per ring
	param.det_per_ring = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Is TOF data used?
	param.TOF = getScalarBool(prhs[ind], ind);
	ind++;

	// Variance of the Gaussian TOF
	param.sigma_x = getScalarFloat(prhs[ind], ind);
	ind++;

	// Centers of the TOF-bins
	param.TOFCenters = getSingles(prhs[ind], "solu");
	ind++;

	// Index offset for TOF subsets
	param.nBins = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Are status messages displayed
	const uint8_t verbose = getScalarBool(prhs[ind], ind);
	ind++;

	// Number of physical CPU cores
	const uint32_t nCores = getScalarUInt32(prhs[ind], ind);
	ind++;
	//if (DEBUG)
	//	mexPrintf("nCores = %u\n", nCores);

	// Is raw list-mode data used
	param.raw = getScalarBool(prhs[ind], ind);
	ind++;

	const uint32_t type = getScalarUInt32(prhs[ind], ind);
	ind++;

	param.listMode = getScalarUInt8(prhs[ind], ind);
	ind++;

	param.projType = getScalarUInt32(prhs[ind], ind);
	ind++;

	param.currentSubset = getScalarUInt32(prhs[ind], ind);
	ind++;

	param.nMeas = getScalarUInt64(prhs[ind], ind);
	ind++;

	//46

	const bool CT = getScalarBool(getField(options, 0, "CT"), ind);
	const bool SPECT = getScalarBool(getField(options, 0, "SPECT"), ind);

	param.pitch = getScalarBool(getField(options, 0, "pitch"), ind);
	param.subsetType = getScalarInt32(getField(options, 0, "subset_type"), ind);
	param.subsets = getScalarUInt32(getField(options, 0, "subsets"), ind);
	if (CT) {
		param.size_y = getScalarUInt32(getField(options, 0, "nColsD"), ind);
		param.dPitchXY = getScalarFloat(getField(options, 0, "dPitchY"), ind);
	} else if (SPECT) {
		param.size_y = getScalarUInt32(getField(options, 0, "nColsD"), ind);
		param.dPitchXY = getScalarFloat(getField(options, 0, "crXY"), ind);
		param.colL = getScalarFloat(getField(options, 0, "colL"), ind);
		param.colD = 2 * getScalarFloat(getField(options, 0, "colR"), ind);
		param.dSeptal = getScalarFloat(getField(options, 0, "dSeptal"), ind);
		param.nRaySPECT = getScalarFloat(getField(options, 0, "nRaySPECT"), ind);
		param.hexOrientation = getScalarFloat(getField(options, 0, "hexOrientation"), ind);
		param.coneMethod = getScalarFloat(getField(options, 0, "coneMethod"), ind);
	} else {
		param.size_y = getScalarUInt32(getField(options, 0, "Nang"), ind);
		param.dPitchXY = getScalarFloat(getField(options, 0, "cr_p"), ind);
	}
	
	if (DEBUG)
		mexPrintf("param.colD = %f\n", param.colD);
	param.nRays2D = getScalarUInt16(options, 0, "n_rays_transaxial");
	param.nRays3D = getScalarUInt16(options, 0, "n_rays_axial");
	param.useMaskFP = getScalarBool(options, 0, "useMaskFP");
	param.useMaskBP = getScalarBool(options, 0, "useMaskBP");
	if (DEBUG)
		mexPrintf("param.nRays2D = %d\n", param.nRays2D);
	if (param.useMaskFP)
		param.maskFP = getUint8s(options, "maskFP");
	if (param.useMaskBP)
		param.maskBP = getUint8s(options, "maskBP");
	const uint32_t nLayers = getScalarUInt32(options, 0, "nLayers");
	if (nLayers > 1)
		param.nLayers = true;
	if (DEBUG)
		mexPrintf("param.nLayers = %d\n", param.nLayers);
	param.computeSensIm = getScalarBool(getField(options, 0, "compute_sensitivity_image"), ind);
	param.rings = getScalarInt32(getField(options, 0, "rings"), ind);
	if (DEBUG)
		mexPrintf("param.rings = %d\n", param.rings);

	if (param.listMode > 0 && !param.computeSensIm) {
		if (DEBUG)
			mexPrintf("param.listMode = %d\n", param.listMode);
		const uint64_t* index = getUint64s(options, "summa");
		if (DEBUG)
			mexPrintf("index = %d\n", index[0]);
		x = getSingles(options, "x", 0);
		if (DEBUG)
			mexPrintf("x = %f\n", x[0]);
		x = &x[index[param.currentSubset] * 6LL];
	}
	// Number of measurements/LORs
	int64_t loop_var_par = pituus;

	if (DEBUG)
		mexPrintf("loop_var_par = %u\n", loop_var_par);

	// Implementation 1
	if (type == 0u) {

		mexErrMsgTxt("Unsupported implementation!");

	}
	// Implementation 4
	else if (type == 1u) {

		if (nrhs < 44)
			mexErrMsgTxt("Too few input arguments.  There must be at least 44.");
		else if (nrhs > 44)
			mexErrMsgTxt("Too many input arguments.  There can be at most 44.");

		if (nlhs != 2)
			mexErrMsgTxt("Invalid number of output arguments. There has to be two.");

		param.x_center = getSingles(options, "x_center");
		param.y_center = getSingles(options, "y_center");
		param.z_center = getSingles(options, "z_center");
		param.bmin = getScalarFloat(options, 0, "bmin");
		param.bmax = getScalarFloat(options, 0, "bmax");
		param.Vmax = getScalarFloat(options, 0, "Vmax");
		if (param.projType == 2)
			param.orthWidth = getScalarFloat(options, 0, "tube_width_z");
		else if (param.projType == 3) {
			param.orthWidth = getScalarFloat(options, 0, "tube_radius");
			param.V = getSingles(options, "V");
		}

		// Small constant to prevent division by zero
		param.epps = getScalarFloat(prhs[ind], ind);
		ind++;

		// Current estimates
		CAST* input = (CAST*)mxGetData(prhs[ind]);
		ind++;

		const size_t N = static_cast<size_t>(param.Nx) * static_cast<size_t>(param.Ny) * static_cast<size_t>(param.Nz);

		// If 1, do not compute the normalization constant anymore
		param.noSensImage = getScalarBool(prhs[ind], ind);
		ind++;

		const uint8_t fp = getScalarUInt8(prhs[ind], ind);
		ind++;

		if (fp == 1)
			param.useMaskBP = false;

		// 53

		size_t imDim;

		if (param.noSensImage)
			imDim = 1ULL;
		else
			imDim = N;

		plhs[1] = mxCreateNumericMatrix(imDim, 1, mxSINGLE_CLASS, mxREAL);

		// Normalization constants
		float* SensIm = getSingles(plhs[1], "solu");

		if (fp == 1)
			plhs[0] = mxCreateNumericMatrix(pituus * param.nBins, 1, mxSINGLE_CLASS, mxREAL);
		else
			plhs[0] = mxCreateNumericMatrix(N, 1, mxSINGLE_CLASS, mxREAL);

		// Right-hand side
		float* output = getSingles(plhs[0], "solu");

		//clock_t time = clock();

		//if (nrhs < 56)
		//	mexErrMsgTxt("Incorrect number of input arguments. There has to be 56.");

		// Crystal pitch in z-direction (for multi-ray)
		param.dPitchXY = getScalarFloat(prhs[ind], ind);
		ind++;

		param.dPitchZ = getScalarFloat(prhs[ind], ind);
		ind++;

		param.nProjections = getScalarUInt32(prhs[ind], ind);
		ind++;

		std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

		projectorType123Implementation4(param, loop_var_par, output, x, z, input, CT, SPECT, fp, SensIm, detIndices);

		std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

		std::chrono::duration<float> time_span = std::chrono::duration_cast<std::chrono::duration<float>>(t2 - t1);

		if (verbose > 1) {
			mexPrintf("Improved Siddon took %f seconds\n", ((float)time_span.count()));
			mexEvalString("pause(.001);");
		}
	}
	// Precomputation phase for implementation 1
	else if (type == 3u) {
		mexErrMsgTxt("Unsupported implementation!");

	}
}