/**************************************************************************
* Implements either the improved Siddon's algorithm, orthogonal 
* distance based projector or volume-based projector for OMEGA.
*
* Implementation 1 supports only Siddon's algorithm and outputs the sparse
* system matrix. Implementation 4 supports all projectors and is computed
* matrix-free.
* 
* This file contains the mex-functions for both the implementation type 1
* and type 4.
* 
* A precomputed vector that has the number of voxels a LOR/ray passes, as 
* well as the number of voxels that have been traversed in previous 
* LORs/rays is required if options.precompute_lor = true, but is optional 
* otherwise. If the precompute_lor options is set to false and implementation 
* type 4 has been selected then all the measurements are investigated, but 
* only those intercepting the FOV will be included.
* 
* Both implementations use OpenMP for parallellization (if available, 
* otherwise the code will be sequential with no parallelization).
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
#include "projector_functions.h"
 
using namespace std;


void mexFunction(int nlhs, mxArray* plhs[],
	int nrhs, const mxArray* prhs[])

{
	// Check for the number of input and output arguments
	if (nrhs < 36)
		mexErrMsgTxt("Too few input arguments. There must be at least 46.");
	else if (nrhs > 65)
		mexErrMsgTxt("Too many input arguments. There can be at most 65.");

	if (nlhs > 3 || nlhs < 1)
		mexErrMsgTxt("Invalid number of output arguments. There can be at most three.");

	int ind = 0;
	// Load the input arguments
	// Image size in y-direction
	const uint32_t Ny = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Image size in x-direction
	const uint32_t Nx = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Image size in z-direction
	const uint32_t Nz = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Distance between adjacent pixels in x-direction
	const double d = getScalarDouble(prhs[ind], ind);
	ind++;

	// Distance between adjacent pixels in z-direction
	const double dz = getScalarDouble(prhs[ind], ind);
	ind++;

	// Distance of the pixel grid from the origin in y-direction
	const double by = getScalarDouble(prhs[ind], ind);
	ind++;

	// Distance of the pixel grid from the origin in x-direction
	const double bx = getScalarDouble(prhs[ind], ind);
	ind++;

	// Distance of the pixel grid from the origin in z-direction
	const double bz = getScalarDouble(prhs[ind], ind);
	ind++;

	// Coordinates of the detectors in z-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const double* z_det = (double*)mxGetDoubles(prhs[ind]);
#else
	const double* z_det = (double*)mxGetData(prhs[ind]);
#endif

	const vector<double> z_det_vec(z_det, z_det + mxGetNumberOfElements(prhs[ind]));
	ind++;

	// Coordinates of the detectors in x-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const double* x = (double*)mxGetDoubles(prhs[ind]);
#else
	const double* x = (double*)mxGetData(prhs[ind]);
#endif
	ind++;

	// Coordinates of the detectors in y-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const double* y = (double*)mxGetDoubles(prhs[ind]);
#else
	const double* y = (double*)mxGetData(prhs[ind]);
#endif
	ind++;

	// Distance between adjacent pixels in y-direction
	const double dy = getScalarDouble(prhs[ind], ind);
	ind++;

	// Coordinates of the pixel planes (boundaries of the pixels) in y-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const double* yy = (double*)mxGetDoubles(prhs[ind]);
#else
	const double* yy = (double*)mxGetData(prhs[ind]);
#endif

	// from array to std::vector
	const vector<double> yy_vec(yy, yy + mxGetNumberOfElements(prhs[ind]));
	ind++;

	// Coordinates of the pixel planes (boundaries of the pixels) in x-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const double* xx = (double*)mxGetDoubles(prhs[ind]);
#else
	const double* xx = (double*)mxGetData(prhs[ind]);
#endif

	const vector<double> xx_vec(xx, xx + mxGetNumberOfElements(prhs[ind]));
	ind++;

	// Number of sinograms used
	const uint32_t NSinos = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Number of slices included
	const uint32_t NSlices = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Number of detector indices
	const uint32_t size_x = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Maximum value of the z-direction detector coordinates
	const double zmax = getScalarDouble(prhs[ind], ind);
	ind++;

	// attenuation values (attenuation images)
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const double* atten = (double*)mxGetDoubles(prhs[ind]);
#else
	const double* atten = (double*)mxGetData(prhs[ind]);
#endif
	ind++;

	// Normalization coefficients
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const float* norm_coef = (float*)mxGetSingles(prhs[ind]);
#else
	const float* norm_coef = (float*)mxGetData(prhs[ind]);
#endif
	ind++;

	// Randoms
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const float* randoms = (float*)mxGetSingles(prhs[ind]);
#else
	const float* randoms = (float*)mxGetData(prhs[ind]);
#endif
	ind++;

	// Number of measurements/LORs
	const int64_t pituus = getScalarInt64(prhs[ind], ind);
	ind++;

	// Is the attenuation correction included
	const bool attenuation_correction = getScalarBool(prhs[ind], ind);
	ind++;

	// Is the normalization correction included
	const bool normalization = getScalarBool(prhs[ind], ind);
	ind++;

	// Is the randoms/scatter correction included
	const bool randoms_correction = getScalarBool(prhs[ind], ind);
	ind++;

	// Is scatter correction included as multiplication (system matrix)
	const bool scatter = getScalarBool(prhs[ind], ind);
	ind++;

	// Scatter data
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const double* scatter_coef = (double*)mxGetDoubles(prhs[ind]);
#else
	const double* scatter_coef = (double*)mxGetData(prhs[ind]);
#endif
	ind++;

	// Global correction factor
	const double global_factor = getScalarDouble(prhs[ind], ind);
	ind++;

	// Number of voxels the current LOR/ray traverses, precomputed data ONLY
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const uint16_t* lor1 = (uint16_t*)mxGetUint16s(prhs[ind]);
#else
	const uint16_t* lor1 = (uint16_t*)mxGetData(prhs[ind]);
#endif
	ind++;

	// For sinogram data, the indices of the detectors corresponding to the current sinogram bin
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const uint32_t* xy_index = (uint32_t*)mxGetUint32s(prhs[ind]);
#else
	const uint32_t* xy_index = (uint32_t*)mxGetData(prhs[ind]);
#endif
	ind++;

	// Same as above, but for z-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const uint16_t* z_index = (uint16_t*)mxGetUint16s(prhs[ind]);
#else
	const uint16_t* z_index = (uint16_t*)mxGetData(prhs[ind]);
#endif
	ind++;

	// Total number of sinograms
	const uint32_t TotSinos = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Detector pair numbers, for raw list-mode data
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const uint16_t* L = (uint16_t*)mxGetUint16s(prhs[ind]);
#else
	const uint16_t* L = (uint16_t*)mxGetData(prhs[ind]);
#endif
	const size_t numRows = mxGetNumberOfElements(prhs[ind]);
	ind++;

	// Location (ring numbers) of pseudo rings, if present
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const uint32_t* pseudos = (uint32_t*)mxGetUint32s(prhs[ind]);
#else
	const uint32_t* pseudos = (uint32_t*)mxGetData(prhs[ind]);
#endif
	const uint32_t pRows = (uint32_t)mxGetM(prhs[ind]);
	ind++;

	// Number of detectors per ring
	const uint32_t det_per_ring = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Is TOF data used?
	const bool TOF = getScalarBool(prhs[ind], ind);
	ind++;

	// Size of single TOF-subset
	const int64_t TOFSize = getScalarInt64(prhs[ind], ind);
	ind++;

	// Variance of the Gaussian TOF
	const double sigma_x = getScalarDouble(prhs[ind], ind);
	ind++;

	// Centers of the TOF-bins
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const double* TOFCenter = (double*)mxGetDoubles(prhs[ind]);
#else
	const double* TOFCenter = (double*)mxGetData(prhs[ind]);
#endif
	ind++;

	// Index offset for TOF subsets
	const int64_t nBins = getScalarInt64(prhs[ind], ind);
	ind++;


	const uint32_t dec_v = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Are status messages displayed
	const bool verbose = getScalarBool(prhs[ind], ind);
	ind++;

	// Number of physical CPU cores
	const uint32_t nCores = getScalarUInt32(prhs[ind], ind);
	ind++;

	// Is raw list-mode data used
	const bool raw = getScalarBool(prhs[ind], ind);
	ind++;

	const uint32_t type = getScalarUInt32(prhs[ind], ind);
	ind++;

	const uint8_t list_mode_format = getScalarUInt8(prhs[ind], ind);
	ind++;

	//46

	// Number of measurements/LORs
	int64_t loop_var_par = pituus;

	// The maximum elements of the pixel space in both x- and y-directions
	const double maxyy = yy_vec.back();
	const double maxxx = xx_vec.back();

	// Implementation 1, with precomputed_lor = true
	if (type == 0u) {

		if (nrhs < 50)
			mexErrMsgTxt("Too few input arguments. There must be at least 50.");
		else if (nrhs > 65)
			mexErrMsgTxt("Too many input arguments. There can be at most 65.");

		if (nlhs != 2)
			mexErrMsgTxt("Invalid number of output arguments. There has to be two.");

		// Total number of voxels traversed by the LORs at the specific LOR
		// e.g. if the first LOR traverses through 10 voxels then at the second lor the value is 10
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		const uint64_t* lor2 = (uint64_t*)mxGetUint64s(prhs[ind]);
#else
		const uint64_t* lor2 = (uint64_t*)mxGetData(prhs[ind]);
#endif
		ind++;

		// Total number of non-zero values
		const uint64_t summa = getScalarUInt64(prhs[ind], ind);
		ind++;

		// Is this a transmission image reconstruction (true) or emission image reconstruction (false)
		const bool attenuation_phase = getScalarBool(prhs[ind], ind);
		ind++;

		const uint32_t projector_type = getScalarUInt32(prhs[ind], ind);
		ind++;

		//50

		// output sizes
		const mwSize N = static_cast<mwSize>(Nx) * static_cast<mwSize>(Ny) * static_cast<mwSize>(Nz);
		const mwSize nzmax = summa;
		const mwSize rows = pituus;

		// Create the MATLAB sparse matrix
		plhs[0] = mxCreateSparse(N, rows, nzmax, mxREAL);

		// Non-zero elements of the matrix
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		double* elements = (double*)mxGetDoubles(plhs[0]);
#else
		double* elements = (double*)mxGetData(plhs[0]);
#endif

		// Row indices
		size_t* indices = reinterpret_cast<size_t*>(mxGetIr(plhs[0]));

		// Column indices
		size_t* lor = reinterpret_cast<size_t*>(mxGetJc(plhs[0]));

		for (size_t kk = 0; kk <= loop_var_par; kk++)
			lor[kk] = lor2[kk];

		// If doing Inveon attenuation
		if (attenuation_phase)
			plhs[1] = mxCreateNumericMatrix(pituus, 1, mxDOUBLE_CLASS, mxREAL);
		else
			plhs[1] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);

#ifdef MX_HAS_INTERLEAVED_COMPLEX
		double* ll = (double*)mxGetDoubles(plhs[1]);
#else
		double* ll = (double*)mxGetData(plhs[1]);
#endif

		// Timing
		std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

		if (projector_type == 2u) {

			// Width of the TOR
			const double crystal_size = getScalarDouble(prhs[ind], ind);
			ind++;

			// Coordinates of the pixel centers in y-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			double* x_center = (double*)mxGetDoubles(prhs[ind]);
#else
			double* x_center = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			// Coordinates of the pixel centers in x-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			double* y_center = (double*)mxGetDoubles(prhs[ind]);
#else
			double* y_center = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			// Coordinates of the pixel centers in z-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			const double* z_center = (double*)mxGetDoubles(prhs[ind]);
#else
			const double* z_center = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			const double crystal_size_z = getScalarDouble(prhs[ind], ind);
			ind++;

#ifndef CT
			// run the Orthogonal distance based ray tracer algorithm, precomputed_lor = true
			orth_siddon_precomputed(loop_var_par, size_x, zmax, indices, elements, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, x, y, z_det,
				NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, lor1, lor2, xy_index, z_index, TotSinos, L, pseudos,
				pRows, det_per_ring, raw, attenuation_phase, ll, crystal_size, crystal_size_z, y_center, x_center, z_center, global_factor, scatter, 
				scatter_coef, nCores, list_mode_format);
#endif

			std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

			if (verbose) {
				mexPrintf("Orthogonal distance based ray tracer took %f seconds\n", ((float)time_span.count()));
				mexEvalString("pause(.001);");
			}
		}
		else if (projector_type == 1u) {

#ifdef CT
			// Width of the TOR
			const double crystal_size = getScalarDouble(prhs[ind], ind);
			ind++;

			// Coordinates of the pixel centers in y-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			double* x_center = (double*)mxGetDoubles(prhs[ind]);
#else
			double* x_center = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			// Coordinates of the pixel centers in x-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			double* y_center = (double*)mxGetDoubles(prhs[ind]);
#else
			double* y_center = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			// Coordinates of the pixel centers in z-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			const double* z_center = (double*)mxGetDoubles(prhs[ind]);
#else
			const double* z_center = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			const double crystal_size_z = getScalarDouble(prhs[ind], ind);
			ind++;

			const double bmin = getScalarDouble(prhs[ind], ind);
			ind++;

			const double bmax = getScalarDouble(prhs[ind], ind);
			ind++;

			const double Vmax = getScalarDouble(prhs[ind], ind);
			ind++;

#ifdef MX_HAS_INTERLEAVED_COMPLEX
			const double* V = (double*)mxGetDoubles(prhs[ind]);
#else
			const double* V = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			const uint32_t subsets = getScalarUInt32(prhs[ind], ind);
			ind++;

#ifdef MX_HAS_INTERLEAVED_COMPLEX
			const double* angles = (double*)mxGetDoubles(prhs[ind]);
#else
			const double* angles = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			const uint32_t size_y = getScalarUInt32(prhs[ind], ind);
			ind++;

			const double dPitch = getScalarDouble(prhs[ind], ind);
			ind++;

			const int64_t nProjections = getScalarInt64(prhs[ind], ind);
			ind++;
#else
			const uint32_t subsets = 1U;
			const double* angles = nullptr;
			const uint32_t size_y = 1U;
			const double dPitch = 0.;
			const int64_t nProjections = 0LL;
#endif

			// run the Improved Siddon's algorithm, precomputed_lor = true
			improved_siddon_precomputed(loop_var_par, size_x, zmax, indices, elements, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, x, y, z_det,
				NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, lor1, lor2, xy_index, z_index, TotSinos, L, pseudos,
				pRows, det_per_ring, raw, attenuation_phase, ll, global_factor, scatter, scatter_coef, subsets, angles, size_y, dPitch, nProjections, 
				nCores, list_mode_format);

			//for (size_t gg = 0; gg < loop_var_par; gg++) {
			//	if (indices[gg] >= Nx * Ny * Nz) {
			//		mexPrintf("indices[gg] = %u\n", indices[gg]);
			//		mexEvalString("pause(.001);");
			//	}
			//}

			std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

			if (verbose) {
				mexPrintf("Improved Siddon took %f seconds\n", ((float)time_span.count()));
				mexEvalString("pause(.001);");
			}
		}
		else if ((projector_type == 3u)) {

			// Width of the TOR
			const double crystal_size = getScalarDouble(prhs[ind], ind);
			ind++;

			// Coordinates of the pixel centers in y-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			double* x_center = (double*)mxGetDoubles(prhs[ind]);
#else
			double* x_center = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			// Coordinates of the pixel centers in x-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			double* y_center = (double*)mxGetDoubles(prhs[ind]);
#else
			double* y_center = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			// Coordinates of the pixel centers in z-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			const double* z_center = (double*)mxGetDoubles(prhs[ind]);
#else
			const double* z_center = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			const double crystal_size_z = getScalarDouble(prhs[ind], ind);
			ind++;

			const double bmin = getScalarDouble(prhs[ind], ind);
			ind++;

			const double bmax = getScalarDouble(prhs[ind], ind);
			ind++;

			const double Vmax = getScalarDouble(prhs[ind], ind);
			ind++;

#ifdef MX_HAS_INTERLEAVED_COMPLEX
			const double* V = (double*)mxGetDoubles(prhs[ind]);
#else
			const double* V = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

#ifdef CT
			const uint32_t subsets = getScalarUInt32(prhs[ind], ind);
			ind++;

#ifdef MX_HAS_INTERLEAVED_COMPLEX
			const double* angles = (double*)mxGetDoubles(prhs[ind]);
#else
			const double* angles = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			const uint32_t size_y = getScalarUInt32(prhs[ind], ind);
			ind++;

			const double dPitch = getScalarDouble(prhs[ind], ind);
			ind++;

			const int64_t nProjections = getScalarInt64(prhs[ind], ind);
			ind++;
#else
			const uint32_t subsets = 1U;
			const double* angles = nullptr;
			const uint32_t size_y = 1U;
			const double dPitch = 0.;
			const int64_t nProjections = 0LL;
#endif

			// run the Orthogonal distance based ray tracer algorithm, precomputed_lor = true
			vol_siddon_precomputed(loop_var_par, size_x, zmax, indices, elements, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, x, y, z_det,
				NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, lor1, lor2, xy_index, z_index, TotSinos, L, pseudos,
				pRows, det_per_ring, raw, attenuation_phase, ll, crystal_size, crystal_size_z, y_center, x_center, z_center, global_factor, bmin, 
				bmax, Vmax, V, scatter, scatter_coef, subsets, angles, size_y, dPitch, nProjections, nCores, list_mode_format);

			std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

			if (verbose) {
				mexPrintf("Volume-based ray tracer took %f seconds\n", ((float)time_span.count()));
				mexEvalString("pause(.001);");
			}
		}
		else
			mexErrMsgTxt("Unsupported projector");

	}
	// Implementation 4
	else if (type == 1u) {

		if (nrhs < 53)
			mexErrMsgTxt("Too few input arguments.  There must be at least 53.");
		else if (nrhs > 65)
			mexErrMsgTxt("Too many input arguments.  There can be at most 65.");

		if (nlhs != 2)
			mexErrMsgTxt("Invalid number of output arguments. There has to be two.");

		// Small constant to prevent division by zero
		const double epps = getScalarDouble(prhs[ind], ind);
		ind++;

		// Measurement data
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		const float* Sino = (float*)mxGetSingles(prhs[ind]);
#else
		const float* Sino = (float*)mxGetData(prhs[ind]);
#endif
		ind++;

		// Current estimates
		double* osem_apu = (double*)mxGetData(prhs[ind]);
		ind++;

		// Projector used
		const uint32_t projector_type = getScalarUInt32(prhs[ind], ind);
		ind++;

		const size_t N = static_cast<size_t>(Nx) * static_cast<size_t>(Ny) * static_cast<size_t>(Nz);

		// If 1, do not compute the normalization constant anymore
		const bool no_norm = getScalarBool(prhs[ind], ind);
		ind++;

		// precomputed_lor = false
		const bool precompute = getScalarBool(prhs[ind], ind);
		ind++;

		const uint8_t fp = getScalarUInt8(prhs[ind], ind);
		ind++;

		// 53

		size_t imDim;

		if (no_norm)
			imDim = 1ULL;
		else
			imDim = N;

		plhs[0] = mxCreateNumericMatrix(imDim, 1, mxDOUBLE_CLASS, mxREAL);

		// Normalization constants
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		double* Summ = (double*)mxGetDoubles(plhs[0]);
#else
		double* Summ = (double*)mxGetData(plhs[0]);
#endif

		if (fp == 1)
			plhs[1] = mxCreateNumericMatrix(pituus, 1, mxDOUBLE_CLASS, mxREAL);
		else
			plhs[1] = mxCreateNumericMatrix(N, 1, mxDOUBLE_CLASS, mxREAL);

		// Right-hand side
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		double* rhs = (double*)mxGetDoubles(plhs[1]);
#else
		double* rhs = (double*)mxGetData(plhs[1]);
#endif

		//clock_t time = clock();
		//std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

		// Orthogonal
		if (projector_type == 2u) {

			if (nrhs < 58)
				mexErrMsgTxt("Incorrect number of input arguments. There has to be at least 58.");

			// Width of the strip in 2D case
			const double crystal_size = getScalarDouble(prhs[ind], ind);
			ind++;

			// Coordinates of the pixel centers in y-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			double* x_center = (double*)mxGetDoubles(prhs[ind]);
#else
			double* x_center = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			// Coordinates of the pixel centers in x-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			double* y_center = (double*)mxGetDoubles(prhs[ind]);
#else
			double* y_center = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			// Coordinates of the pixel centers in z-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			const double* z_center = (double*)mxGetDoubles(prhs[ind]);
#else
			const double* z_center = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			// Width of the TOR in 3D case
			const double crystal_size_z = getScalarDouble(prhs[ind], ind);
			ind++;

#ifndef CT
			if (precompute) {
				sequential_orth_siddon(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, randoms, x, y, z_det, 
					NSlices, Nx, Ny, Nz, d, dz,	bx, by, bz, attenuation_correction, normalization, randoms_correction, lor1, xy_index, z_index, 
					TotSinos, epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, crystal_size, x_center, y_center, z_center, crystal_size_z, 
					no_norm, dec_v, global_factor, fp, scatter, scatter_coef, nCores);
			}
			else {
				sequential_orth_siddon_no_precomp(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, randoms,
					x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, randoms_correction, xy_index, z_index, 
					TotSinos, epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, crystal_size, x_center, y_center, z_center, crystal_size_z, 
					no_norm, dec_v, global_factor, fp, list_mode_format, scatter, scatter_coef, nCores);
			}
#endif
		}
		// Improved Siddon
		else if (projector_type == 1u) {

			if (nrhs < 56)
				mexErrMsgTxt("Incorrect number of input arguments. There has to be 56.");

			// Number of rays in Siddon (transaxial)
			uint16_t n_rays = getScalarUInt16(prhs[ind], ind);
			ind++;

			// Number of rays in Siddon (axial)
			uint16_t n_rays3D = getScalarUInt16(prhs[ind], ind);
			ind++;

			// Crystal pitch in z-direction (for multi-ray)
			const double cr_pz = getScalarDouble(prhs[ind], ind);
			ind++;

#ifdef CT
			const uint32_t subsets = getScalarUInt32(prhs[ind], ind);
			ind++;

#ifdef MX_HAS_INTERLEAVED_COMPLEX
			const double* angles = (double*)mxGetDoubles(prhs[ind]);
#else
			const double* angles = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			const uint32_t size_y = getScalarUInt32(prhs[ind], ind);
			ind++;

			const double dPitch = getScalarDouble(prhs[ind], ind);
			ind++;

			const int64_t nProjections = getScalarInt64(prhs[ind], ind);
			ind++;
#else
			const uint32_t subsets = 1U;
			const double* angles = nullptr;
			const uint32_t size_y = 1U;
			const double dPitch = 0.;
			const int64_t nProjections = 0LL;
#endif

			if (precompute) {
				sequential_improved_siddon(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, randoms, x, y,
					z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, randoms_correction, lor1, xy_index, z_index, 
					TotSinos, epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, no_norm, global_factor, fp, scatter, scatter_coef, TOF, TOFSize, 
					sigma_x, TOFCenter, nBins, dec_v, subsets, angles, size_y, dPitch, nProjections, nCores);
			}
			else {
				sequential_improved_siddon_no_precompute(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, randoms, x, y,
					z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, randoms_correction, xy_index, z_index, TotSinos,
					epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, cr_pz, no_norm, n_rays, n_rays3D, global_factor, fp, list_mode_format, 
					scatter, scatter_coef, TOF, TOFSize, sigma_x, TOFCenter, nBins, dec_v, subsets, angles, size_y, dPitch, nProjections, nCores);
			}
		}
		else if ((projector_type == 3u)) {
			if (nrhs < 60)
				mexErrMsgTxt("Incorrect number of input arguments. There has to be 60.");

			// Coordinates of the pixel centers in y-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			double* x_center = (double*)mxGetDoubles(prhs[ind]);
#else
			double* x_center = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			// Coordinates of the pixel centers in x-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			double* y_center = (double*)mxGetDoubles(prhs[ind]);
#else
			double* y_center = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			// Coordinates of the pixel centers in z-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			const double* z_center = (double*)mxGetDoubles(prhs[ind]);
#else
			const double* z_center = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			// Width of the TOR in 3D case
			const double bmin = getScalarDouble(prhs[ind], ind);
			ind++;

			// Width of the TOR in 3D case
			const double bmax = getScalarDouble(prhs[ind], ind);
			ind++;

			// Width of the TOR in 3D case
			const double Vmax = getScalarDouble(prhs[ind], ind);
			ind++;

			// Coordinates of the pixel centers in y-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			const double* V = (double*)mxGetDoubles(prhs[ind]);
#else
			const double* V = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

#ifdef CT
			const uint32_t subsets = getScalarUInt32(prhs[ind], ind);
			ind++;

#ifdef MX_HAS_INTERLEAVED_COMPLEX
			const double* angles = (double*)mxGetDoubles(prhs[ind]);
#else
			const double* angles = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			const uint32_t size_y = getScalarUInt32(prhs[ind], ind);
			ind++;

			const double dPitch = getScalarDouble(prhs[ind], ind);
			ind++;

			const int64_t nProjections = getScalarInt64(prhs[ind], ind);
			ind++;
#else
			const uint32_t subsets = 1U;
			const double* angles = nullptr;
			const uint32_t size_y = 1U;
			const double dPitch = 0.;
			const int64_t nProjections = 0LL;
#endif

			if (precompute) {
				sequential_volume_siddon(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, randoms, x, y, z_det,
					NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, randoms_correction, lor1, xy_index, z_index,
					TotSinos, epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, Vmax, x_center, y_center, z_center, bmin, bmax, V,
					no_norm, dec_v, global_factor, fp, scatter, scatter_coef, subsets, angles, size_y, dPitch, nProjections, nCores);
			}
			else {
				sequential_volume_siddon_no_precomp(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, randoms,
					x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, randoms_correction, xy_index, z_index,
					TotSinos, epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, Vmax, x_center, y_center, z_center, bmin, bmax, V,
					no_norm, dec_v, global_factor, fp, list_mode_format, scatter, scatter_coef, subsets, angles, size_y, dPitch, nProjections, nCores);
			}
		}
	}
	// Implementation 1, precomputed_lor = false
	else if (type == 2u) {

		if (nrhs < 49)
			mexErrMsgTxt("Too few input arguments.  There must be at least 49.");
		else if (nrhs > 61)
			mexErrMsgTxt("Too many input arguments.  There can be at most 61.");

		if (nlhs != 3)
			mexErrMsgTxt("Invalid number of output arguments. There has to be three.");


		// How many elements are preallocated in memory
		const uint32_t ind_size = getScalarUInt32(prhs[ind], ind);
		ind++;

		// Starting ring
		const uint32_t block1 = getScalarUInt32(prhs[ind], ind);
		ind++;

		// End ring
		const uint32_t blocks = getScalarUInt32(prhs[ind], ind);
		ind++;

		// Subset indices
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		const uint32_t* index = (uint32_t*)mxGetUint32s(prhs[ind]);
#else
		const uint32_t* index = (uint32_t*)mxGetData(prhs[ind]);
#endif
		const size_t index_size = mxGetNumberOfElements(prhs[ind]);
		ind++;

		// Projector
		const uint32_t projector_type = getScalarUInt32(prhs[ind], ind);
		ind++;

		// 51

		// Voxel numbers in x-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		const double* iij = (double*)mxGetDoubles(prhs[ind]);
#else
		const double* iij = (double*)mxGetData(prhs[ind]);
#endif
		const vector<double> iij_vec(iij, iij + mxGetNumberOfElements(prhs[ind]));
		ind++;

		// Voxel numbers in y-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		const double* jji = (double*)mxGetDoubles(prhs[ind]);
#else
		const double* jji = (double*)mxGetData(prhs[ind]);
#endif
		const vector<double> jjk_vec(jji, jji + mxGetNumberOfElements(prhs[ind]));
		ind++;

		// Voxel numbers in z-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		const double* kkj = (double*)mxGetDoubles(prhs[ind]);
#else
		const double* kkj = (double*)mxGetData(prhs[ind]);
#endif
		const vector<double> kkj_vec(kkj, kkj + mxGetNumberOfElements(prhs[ind]));
		ind++;

#ifdef CT
		const uint32_t subsets = getScalarUInt32(prhs[ind], ind);
		ind++;

#ifdef MX_HAS_INTERLEAVED_COMPLEX
		const double* angles = (double*)mxGetDoubles(prhs[ind]);
#else
		const double* angles = (double*)mxGetData(prhs[ind]);
#endif
		ind++;

		const uint32_t size_y = getScalarUInt32(prhs[ind], ind);
		ind++;

		const double dPitch = getScalarDouble(prhs[ind], ind);
		ind++;

		const int64_t nProjections = getScalarInt64(prhs[ind], ind);
		ind++;
#else
		const uint32_t subsets = 1U;
		const double* angles = nullptr;
		const uint32_t size_y = 1U;
		const double dPitch = 0.;
		const int64_t nProjections = 0LL;
#endif

		// Number of LORs
		if (index_size > 1ULL && !raw && !list_mode_format) {
			loop_var_par = index_size;
		}
		else if (!raw || list_mode_format) {
			loop_var_par = static_cast<size_t>(NSinos) * static_cast<size_t>(size_x);
		}
		else {
			loop_var_par = numRows / 2ULL;
		}

		//mexPrintf("loop_var_par = %u\n", loop_var_par);

		plhs[0] = mxCreateNumericMatrix(loop_var_par, 1, mxUINT16_CLASS, mxREAL);

#ifdef MX_HAS_INTERLEAVED_COMPLEX
		uint16_t* lor = (uint16_t*)mxGetUint16s(plhs[0]);
#else
		uint16_t* lor = (uint16_t*)mxGetData(plhs[0]);
#endif

		vector<uint32_t> indices;

		vector<double> elements;

		// Reserve some memory
		indices.reserve(ind_size);
		elements.reserve(ind_size);

		// The maximum elements of the pixel space in both x- and y-directions
		const double maxyy = yy_vec.back();
		const double maxxx = xx_vec.back();

		//clock_t time = clock();
		std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

		uint32_t lj = 0U;

		if (projector_type == 1u) {

			// run the Improved Siddon's algorithm
			lj = improved_siddon_no_precompute(loop_var_par, size_x, zmax, TotSinos, indices, elements, lor, maxyy, maxxx, xx_vec, dy,
				yy_vec, atten, norm_coef, x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, index, attenuation_correction, normalization, raw, 
				det_per_ring, blocks, block1, L, pseudos, pRows, global_factor, scatter, scatter_coef, subsets, angles, xy_index, z_index, size_y, 
				dPitch, nProjections, list_mode_format);
		}
		else if (projector_type == 2u) {

			if (nrhs != 56)
				mexErrMsgTxt("Incorrect number of input arguments. There has to be 56.");

			const double crystal_size = getScalarDouble(prhs[ind], ind);
			ind++;
			
			// Coordinates of the pixel centers in y-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			double* x_center = (double*)mxGetDoubles(prhs[ind]);
#else
			double* x_center = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			// Coordinates of the pixel centers in x-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			double* y_center = (double*)mxGetDoubles(prhs[ind]);
#else
			double* y_center = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			// Coordinates of the pixel centers in z-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
			const double* z_center = (double*)mxGetDoubles(prhs[ind]);
#else
			const double* z_center = (double*)mxGetData(prhs[ind]);
#endif
			ind++;

			const double crystal_size_z = getScalarDouble(prhs[ind], ind);
			ind++;

			// run the Orthogonal Siddon algorithm
			//lj = orth_siddon_no_precompute(loop_var_par, size_x, zmax, TotSinos, indices, elements, lor, maxyy, maxxx, xx_vec, dy,
			//	yy_vec, atten, norm_coef, x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, index, attenuation_correction, normalization, raw, 
			//	det_per_ring, blocks, block1, L, pseudos, pRows, crystal_size, crystal_size_z, y_center, x_center, z_center, dec_v);
		}
		// Original Siddon's ray tracer
		else if (projector_type == 0u) {

			if (nrhs < 54)
				mexErrMsgTxt("Too few input arguments.  There must be at least 54.");


			// run the original Siddon's algorithm
			lj = original_siddon_no_precompute(loop_var_par, size_x, zmax, TotSinos, indices, elements, lor, maxyy, maxxx, xx_vec, dy,
				yy_vec, atten, norm_coef, x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, index, attenuation_correction, normalization, raw,
				det_per_ring, blocks, block1, L, pseudos, pRows, iij_vec, jjk_vec, kkj_vec, global_factor, scatter, scatter_coef, subsets, angles, xy_index, z_index, size_y,
				dPitch, nProjections, list_mode_format);
		}

		const size_t outSize1 = static_cast<size_t>(lj) * 2ULL;
		const size_t outSize2 = 1ULL;

		// Create the MATLAB output vectors (row and column indices, elements)

		plhs[1] = mxCreateNumericMatrix(indices.size(), outSize2, mxUINT32_CLASS, mxREAL);

#ifdef MX_HAS_INTERLEAVED_COMPLEX
		uint32_t* outputMatrix2 = (uint32_t*)mxGetUint32s(plhs[1]);
#else
		uint32_t* outputMatrix2 = (uint32_t*)mxGetData(plhs[1]);
#endif

		std::copy(indices.begin(), indices.end(), outputMatrix2);
		indices.erase(indices.begin(), indices.end());
		indices.shrink_to_fit();

		plhs[2] = mxCreateNumericMatrix(elements.size(), outSize2, mxDOUBLE_CLASS, mxREAL);

#ifdef MX_HAS_INTERLEAVED_COMPLEX
		double* outputMatrix3 = (double*)mxGetDoubles(plhs[2]);
#else
		double* outputMatrix3 = (double*)mxGetData(plhs[2]);
#endif

		std::copy(elements.begin(), elements.end(), outputMatrix3);

		std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

		if (verbose) {
			mexPrintf("Function elapsed time is %f seconds\n", ((float)time_span.count()));
			mexEvalString("pause(.001);");
		}

	}
	// Precomputation phase
	else if (type == 3u) {
		if (nrhs < 59)
			mexErrMsgTxt("Too few input arguments.  There must be at least 59.");
		else if (nrhs > 64)
			mexErrMsgTxt("Too many input arguments.  There can be at most 64.");

		if (nlhs != 3)
			mexErrMsgTxt("Invalid number of output arguments. There has to be three.");


		// Starting ring
		const uint32_t block1 = getScalarUInt32(prhs[ind], ind);
		ind++;

		// End ring
		const uint32_t blocks = getScalarUInt32(prhs[ind], ind);
		ind++;

		// Projector
		const uint32_t projector_type = getScalarUInt32(prhs[ind], ind);
		ind++;

		const double crystal_size = getScalarDouble(prhs[ind], ind);
		ind++;

		// Coordinates of the pixel centers in y-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		double* x_center = (double*)mxGetDoubles(prhs[ind]);
#else
		double* x_center = (double*)mxGetData(prhs[ind]);
#endif
		ind++;

		// Coordinates of the pixel centers in x-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		double* y_center = (double*)mxGetDoubles(prhs[ind]);
#else
		double* y_center = (double*)mxGetData(prhs[ind]);
#endif
		ind++;

		// Coordinates of the pixel centers in z-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		const double* z_center = (double*)mxGetDoubles(prhs[ind]);
#else
		const double* z_center = (double*)mxGetData(prhs[ind]);
#endif
		ind++;

		const double crystal_size_z = getScalarDouble(prhs[ind], ind);
		ind++;

		// Width of the TOR in 3D case
		const double bmin = getScalarDouble(prhs[ind], ind);
		ind++;

		// Width of the TOR in 3D case
		const double bmax = getScalarDouble(prhs[ind], ind);
		ind++;

		// Width of the TOR in 3D case
		const double Vmax = getScalarDouble(prhs[ind], ind);
		ind++;

		// Coordinates of the pixel centers in y-direction
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		const double* V = (double*)mxGetDoubles(prhs[ind]);
#else
		const double* V = (double*)mxGetData(prhs[ind]);
#endif
		ind++;

		const uint32_t tyyppi = getScalarUInt32(prhs[ind], ind);
		ind++;

#ifdef CT
		const uint32_t subsets = getScalarUInt32(prhs[ind], ind);
		ind++;

#ifdef MX_HAS_INTERLEAVED_COMPLEX
		const double* angles = (double*)mxGetDoubles(prhs[ind]);
#else
		const double* angles = (double*)mxGetData(prhs[ind]);
#endif
		ind++;

		const uint32_t size_y = getScalarUInt32(prhs[ind], ind);
		ind++;

		const double dPitch = getScalarDouble(prhs[ind], ind);
		ind++;

		const int64_t nProjections = getScalarInt64(prhs[ind], ind);
		ind++;
#else
		const uint32_t subsets = 1U;
		const double* angles = nullptr;
		const uint32_t size_y = 1U;
		const double dPitch = 0.;
		const int64_t nProjections = 0LL;
#endif

		if (raw)
			loop_var_par = numRows / 2ULL;
		else
			loop_var_par = static_cast<size_t>(NSinos) * static_cast<size_t>(size_x);
		uint16_t* lor;
		uint16_t* lor_orth;
		uint16_t* lor_vol;

		plhs[0] = mxCreateNumericMatrix(loop_var_par, 1, mxUINT16_CLASS, mxREAL);

		if (tyyppi == 0) {
			plhs[1] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
			plhs[2] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		}
		else if (tyyppi == 1 || tyyppi == 2) {

			if (tyyppi == 2u)
				plhs[1] = mxCreateNumericMatrix(loop_var_par * 2LL, 1, mxUINT16_CLASS, mxREAL);
			else
				plhs[1] = mxCreateNumericMatrix(loop_var_par, 1, mxUINT16_CLASS, mxREAL);
			plhs[2] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);

		}
		else {
			plhs[1] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
			plhs[2] = mxCreateNumericMatrix(loop_var_par, 1, mxUINT16_CLASS, mxREAL);
		}

#ifdef MX_HAS_INTERLEAVED_COMPLEX
		lor = (uint16_t*)mxGetUint16s(plhs[0]);
		lor_orth = (uint16_t*)mxGetUint16s(plhs[1]);
		lor_vol = (uint16_t*)mxGetUint16s(plhs[2]);
#else
		lor = (uint16_t*)mxGetData(plhs[0]);
		lor_orth = (uint16_t*)mxGetData(plhs[1]);
		lor_vol = (uint16_t*)mxGetData(plhs[2]);
#endif

		improved_siddon_precomputation_phase(loop_var_par, size_x, zmax, TotSinos, lor, maxyy, maxxx, xx_vec, z_det_vec, dy, yy_vec, x, y, z_det, NSlices, Nx, Ny, Nz,
			d, dz, bx, by, bz, block1, blocks, L, pseudos, raw, pRows, det_per_ring, tyyppi, lor_orth, lor_vol, crystal_size, crystal_size_z, x_center, y_center, z_center, 
			bmin, bmax, Vmax, V, angles, size_y, dPitch, nProjections, nCores, list_mode_format);

	}
}