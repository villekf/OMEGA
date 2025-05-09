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
* Copyright (C) 2022-2024 Ville-Veikko Wettenhovi, Niilo Saarlemo
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
#define CAST double
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

	// Image size in y-direction
	param.Ny = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Image size in z-direction
	param.Nz = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Distance between adjacent pixels in x-direction
	param.dx = getScalarDouble(prhs[ind], ind);
	ind++;

	// Distance between adjacent pixels in y-direction
	param.dy = getScalarDouble(prhs[ind], ind);
	ind++;

	// Distance between adjacent pixels in z-direction
	param.dz = getScalarDouble(prhs[ind], ind);
	ind++;

	// Distance of the pixel grid from the origin in x-direction
	param.bx = getScalarDouble(prhs[ind], ind);
	ind++;

	// Distance of the pixel grid from the origin in y-direction
	param.by = getScalarDouble(prhs[ind], ind);
	ind++;

	// Distance of the pixel grid from the origin in z-direction
	param.bz = getScalarDouble(prhs[ind], ind);
	ind++;

	// Coordinates of the detectors in z-direction
	const double* z = getDoubles(prhs[ind], "solu");
	ind++;

	// Coordinates of the detectors in x/y-direction
	double* x = getDoubles(prhs[ind], "solu");
	ind++;

	// Number of detector indices
	param.size_x = getScalarUInt32(prhs[ind], ind);
	ind++;

	// attenuation values (attenuation images)
	param.atten = getDoubles(prhs[ind], "solu");
	ind++;

	// Normalization coefficients
	param.normCoefs = getDoubles(prhs[ind], "solu");
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
	param.scatterCoefs = getDoubles(prhs[ind], "solu");
	ind++;

	// Global correction factor
	param.globalFactor = getScalarDouble(prhs[ind], ind);
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
	param.sigma_x = getScalarDouble(prhs[ind], ind);
	ind++;

	// Centers of the TOF-bins
	param.TOFCenters = getDoubles(prhs[ind], "solu");
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
		param.dPitchXY = getScalarDouble(getField(options, 0, "dPitchY"), ind);
	} else if (SPECT) {
		param.size_y = getScalarUInt32(getField(options, 0, "nColsD"), ind);
		param.dPitchXY = getScalarDouble(getField(options, 0, "crXY"), ind);
		param.rayShiftsDetector = getDoubles(options, "rayShiftsDetector");
		param.rayShiftsSource = getDoubles(options, "rayShiftsSource");
        param.coneOfResponseStdCoeffA = getScalarDouble(getField(options, 0, "coneOfResponseStdCoeffA"), ind);
        param.coneOfResponseStdCoeffB = getScalarDouble(getField(options, 0, "coneOfResponseStdCoeffB"), ind);
        param.coneOfResponseStdCoeffC = getScalarDouble(getField(options, 0, "coneOfResponseStdCoeffC"), ind);
	} else {
		param.size_y = getScalarUInt32(getField(options, 0, "Nang"), ind);
		param.dPitchXY = getScalarDouble(getField(options, 0, "cr_p"), ind);
	}
	param.nRays2D = getScalarDouble(options, 0, "n_rays_transaxial");
	param.nRays3D = getScalarDouble(options, 0, "n_rays_axial");
	param.useMaskFP = getScalarBool(options, 0, "useMaskFP");
	if (param.useMaskFP)
		param.maskFP = getUint8s(options, "maskFP");
	if (param.useMaskBP)
		param.maskBP = getUint8s(options, "maskBP");
	const uint32_t nLayers = getScalarUInt32(options, 0, "nLayers");
	if (nLayers > 1)
		param.nLayers = true;
	param.computeSensIm = getScalarBool(getField(options, 0, "compute_sensitivity_image"), ind);
	param.rings = getScalarInt32(getField(options, 0, "rings"), ind);

	if (param.listMode > 0 && !param.computeSensIm) {
		const uint64_t* index = getUint64s(options, "summa");
		x = getDoubles(options, "x", 0);
		x = &x[index[param.currentSubset] * 6LL];
	}
	// Number of measurements/LORs
	int64_t loop_var_par = pituus;

	if (DEBUG)
		mexPrintf("loop_var_par = %u\n", loop_var_par);

	// Implementation 1
	if (type == 0u) {

		if (nrhs < 42)
			mexErrMsgTxt("Too few input arguments. There must be at least 42.");
		else if (nrhs > 42)
			mexErrMsgTxt("Too many input arguments. There can be at most 42.");

		if (nlhs != 1)
			mexErrMsgTxt("Invalid number of output arguments. There has to be one.");

		// Total number of voxels traversed by the LORs at the specific LOR
		// e.g. if the first LOR traverses through 10 voxels then at the second lor the value is 10
		const uint64_t* lor2 = getUint64s(prhs[ind], "solu");
		ind++;

		// Number of voxels the current LOR/ray traverses, precomputed data ONLY
		const uint16_t* lor1 = getUint16s(prhs[ind], "solu");
		ind++;

		// Total number of non-zero values
		const uint64_t summa = getScalarUInt64(prhs[ind], ind);
		ind++;

		//50

		// output sizes
		const mwSize N = static_cast<mwSize>(param.Nx) * static_cast<mwSize>(param.Ny) * static_cast<mwSize>(param.Nz);
		const mwSize nzmax = summa;
		const mwSize rows = pituus;

		// Create the MATLAB sparse matrix
		plhs[0] = mxCreateSparse(N, rows, nzmax, mxREAL);

		// Non-zero elements of the matrix
		double* elements = getDoubles(plhs[0], "solu");

		// Row indices
		size_t* indices = reinterpret_cast<size_t*>(mxGetIr(plhs[0]));

		// Column indices
		size_t* lor = reinterpret_cast<size_t*>(mxGetJc(plhs[0]));

		for (int64_t kk = 0; kk <= loop_var_par; kk++)
			lor[kk] = lor2[kk];

		// Timing
		std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
		if (param.projType == 1u || param.projType == 11) {

			param.dPitchZ = getScalarDouble(prhs[ind], ind);
			ind++;

			param.nProjections = getScalarUInt32(prhs[ind], ind);
			ind++;

			// run the Improved Siddon's algorithm, precomputed_lor = true
			improved_siddon_precomputed(param, loop_var_par, x, z, elements, indices, lor1, lor2, CT, detIndices);

			std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

			if (verbose > 0) {
				mexPrintf("Improved Siddon took %f seconds\n", ((float)time_span.count()));
				mexEvalString("pause(.001);");
			}
		}
		else
			mexErrMsgTxt("Unsupported projector!");

	}
	// Implementation 4
	else if (type == 1u) {

		if (nrhs < 44)
			mexErrMsgTxt("Too few input arguments.  There must be at least 44.");
		else if (nrhs > 44)
			mexErrMsgTxt("Too many input arguments.  There can be at most 44.");

		if (nlhs != 2)
			mexErrMsgTxt("Invalid number of output arguments. There has to be two.");

		param.x_center = getDoubles(options, "x_center");
		param.y_center = getDoubles(options, "y_center");
		param.z_center = getDoubles(options, "z_center");
		param.bmin = getScalarDouble(options, 0, "bmin");
		param.bmax = getScalarDouble(options, 0, "bmax");
		param.Vmax = getScalarDouble(options, 0, "Vmax");
		if (param.projType == 2)
			param.orthWidth = getScalarDouble(options, 0, "tube_width_z");
		else if (param.projType == 3) {
			param.orthWidth = getScalarDouble(options, 0, "tube_radius");
			param.V = getDoubles(options, "V");
		}

		// Small constant to prevent division by zero
		param.epps = getScalarDouble(prhs[ind], ind);
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

		// 53

		size_t imDim;

		if (param.noSensImage)
			imDim = 1ULL;
		else
			imDim = N;

		plhs[1] = mxCreateNumericMatrix(imDim, 1, mxDOUBLE_CLASS, mxREAL);

		// Normalization constants
		double* SensIm = getDoubles(plhs[1], "solu");

		if (fp == 1)
			plhs[0] = mxCreateNumericMatrix(pituus * param.nBins, 1, mxDOUBLE_CLASS, mxREAL);
		else
			plhs[0] = mxCreateNumericMatrix(N, 1, mxDOUBLE_CLASS, mxREAL);

		// Right-hand side
		double* output = getDoubles(plhs[0], "solu");

		// Crystal pitch in z-direction (for multi-ray)
		param.dPitchXY = getScalarDouble(prhs[ind], ind);
		ind++;

		param.dPitchZ = getScalarDouble(prhs[ind], ind);
		ind++;

		param.nProjections = getScalarUInt32(prhs[ind], ind);
		ind++;

		std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

		projectorType123Implementation4(param, loop_var_par, output, x, z, input, CT, SPECT, fp, SensIm, detIndices);

		std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

		if (verbose > 1) {
			mexPrintf("Improved Siddon took %f seconds\n", ((float)time_span.count()));
			mexEvalString("pause(.001);");
		}
	}
	// Precomputation phase for implementation 1
	else if (type == 3u) {
		if (nrhs < 39)
			mexErrMsgTxt("Too few input arguments.  There must be at least 39.");
		else if (nrhs > 39)
			mexErrMsgTxt("Too many input arguments.  There can be at most 39.");

		if (nlhs != 1)
			mexErrMsgTxt("Invalid number of output arguments. There has to be one.");

		if (DEBUG)
			mexPrintf("type = %u\n", type);

		param.dPitchZ = getScalarDouble(prhs[ind], ind);
		ind++;

		if (DEBUG)
			mexPrintf("ind = %u\n", ind);

		param.nProjections = getScalarUInt32(prhs[ind], ind);
		ind++;

		if (DEBUG)
			mexPrintf("ind = %u\n", ind);

		if (param.raw)
			loop_var_par = numRows / 2ULL;
		else
			loop_var_par = static_cast<size_t>(param.nProjections) * static_cast<size_t>(param.size_x) * static_cast<size_t>(param.size_y);

		if (DEBUG)
			mexPrintf("loop_var_par2 = %u\n", loop_var_par);
		uint16_t* lor;

		plhs[0] = mxCreateNumericMatrix(loop_var_par, 1, mxUINT16_CLASS, mxREAL);

		lor = getUint16s(plhs[0], "solu");

		improved_siddon_precomputation_phase(param, loop_var_par, x, z, lor, CT, detIndices);

	}
}