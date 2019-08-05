/**************************************************************************
* Implements the either the improved Siddon's algorithm or orthogonal 
* distance based projector for OMEGA.
* 
* This file contains the mex-functions for both the implementation type 1
* and type 4. The former in the case of precomputed_lor = true.
* 
* A precomputed vector that has the number of voxels a LOR/ray passes, as 
* well as the number of voxels that have been traversed in previous 
* LORs/rays is required for the implemention type 1, but is optional for 
* type 4. If the precompute_lor options is set to false and implementation 
* type 4 has been selected then all the measurements are investigated, but 
* only those intercepting the FOV will be included.
* 
* Implementation type 1 uses C++11 threads for multithreading (parallel 
* for loop), type 4 uses OpenMP (if available, otherwise the code will be
* sequential with no parallelization).
* 
* Copyright (C) 2019  Ville-Veikko Wettenhovi
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
	if (nrhs < 33)
		mexErrMsgTxt("Too few input arguments. There must be at least 33.");
	else if (nrhs > 47)
		mexErrMsgTxt("Too many input arguments. There can be at most 47.");

	if (nlhs > 3 || nlhs < 1)
		mexErrMsgTxt("Invalid number of output arguments. There can be at most three.");

	// Load the input arguments
	// Image size in y-direction
	const uint32_t Ny = (uint32_t)mxGetScalar(prhs[0]);

	// Image size in x-direction
	const uint32_t Nx = (uint32_t)mxGetScalar(prhs[1]);

	// Image size in z-direction
	const uint32_t Nz = (uint32_t)mxGetScalar(prhs[2]);

	// Distance between adjacent pixels in x-direction
	const double d = (double)mxGetScalar(prhs[3]);

	// Distance between adjacent pixels in z-direction
	const double dz = (double)mxGetScalar(prhs[4]);

	// Distance of the pixel grid from the origin in y-direction
	const double by = (double)mxGetScalar(prhs[5]);

	// Distance of the pixel grid from the origin in x-direction
	const double bx = (double)mxGetScalar(prhs[6]);

	// Distance of the pixel grid from the origin in z-direction
	const double bz = (double)mxGetScalar(prhs[7]);

	// Coordinates of the detectors in z-direction
	const double* z_det = (double*)mxGetData(prhs[8]);

	// Coordinates of the detectors in x-direction
	const double* x = (double*)mxGetData(prhs[9]);

	// Coordinates of the detectors in y-direction
	const double* y = (double*)mxGetData(prhs[10]);

	// Distance between adjacent pixels in y-direction
	const double dy = (double)mxGetScalar(prhs[11]);

	// Coordinates of the pixel planes (boundaries of the pixels) in y-direction
	const double* yy = (double*)mxGetData(prhs[12]);

	// Coordinates of the pixel planes (boundaries of the pixels) in x-direction
	const double* xx = (double*)mxGetData(prhs[13]);

	// Number of sinograms used
	const uint32_t NSinos = (uint32_t)mxGetScalar(prhs[14]);

	// Number of slices included
	const uint32_t NSlices = (uint32_t)mxGetScalar(prhs[15]);

	// from array to std::vector
	const vector<double> yy_vec(yy, yy + mxGetNumberOfElements(prhs[12]));

	const vector<double> xx_vec(xx, xx + mxGetNumberOfElements(prhs[13]));

	// Number of detector indices
	const uint32_t size_x = (uint32_t)mxGetScalar(prhs[16]);

	// Maximum value of the z-direction detector coordinates
	const double zmax = (double)mxGetScalar(prhs[17]);

	// attenuation values
	const double* atten = (double*)mxGetData(prhs[18]);

	// Normalization coefficients
	const double* norm_coef = (double*)mxGetData(prhs[19]);

	// Randoms
	const double* randoms = (double*)mxGetData(prhs[20]);

	// Number of measurements/LORs
	const uint32_t pituus = (uint32_t)mxGetScalar(prhs[21]);

	// Is the attenuation correction included
	const bool attenuation_correction = (bool)mxGetScalar(prhs[22]);

	// Is the attenuation correction included
	const bool normalization = (bool)mxGetScalar(prhs[23]);

	// Is the attenuation correction included
	const bool randoms_correction = (bool)mxGetScalar(prhs[24]);

	// Number of voxels the current LOR/ray traverses
	const uint16_t* lor1 = (uint16_t*)mxGetData(prhs[25]);

	// For sinogram data, the indices of the detectors corresponding to the current sinogram bin
	const uint32_t* xy_index = (uint32_t*)mxGetData(prhs[26]);

	// Same as above, but for z-direction
	const uint16_t* z_index = (uint16_t*)mxGetData(prhs[27]);

	// Total number of sinograms
	const uint32_t TotSinos = (uint32_t)mxGetScalar(prhs[28]);

	// Detector pair numbers, for raw list-mode data
	const uint16_t* L = (uint16_t*)mxGetData(prhs[29]);
	const size_t numRows = mxGetM(prhs[29]);

	// Location (ring numbers) of pseudo rings, if present
	const uint32_t* pseudos = (uint32_t*)mxGetData(prhs[30]);
	const uint32_t pRows = (uint32_t)mxGetM(prhs[30]);

	// Number of detectors per ring
	const uint32_t det_per_ring = (uint32_t)mxGetScalar(prhs[31]);

	// Are status messages displayed
	const bool verbose = (bool)mxGetScalar(prhs[32]);

	// Is raw list-mode data used
	const bool raw = (bool)mxGetScalar(prhs[33]);

	const uint32_t type = (uint32_t)mxGetScalar(prhs[34]);

	// Number of measurements/LORs
	const size_t loop_var_par = pituus;

	// The maximum elements of the pixel space in both x- and y-directions
	const double maxyy = yy_vec.back();
	const double maxxx = xx_vec.back();

	if (type == 0u) {

		if (nrhs < 39)
			mexErrMsgTxt("Too few input arguments.  There must be at least 39.");
		else if (nrhs > 45)
			mexErrMsgTxt("Too many input arguments.  There can be at most 45.");

		if (nlhs != 2)
			mexErrMsgTxt("Invalid number of output arguments. There has to be two.");

		// Total number of voxels traversed by the LORs at the specific LOR
		// e.g. if the first LOR traverses through 10 voxels then at the second lor the value is 10
		const uint64_t* lor2 = (uint64_t*)mxGetData(prhs[35]);
		//mwIndex *lor2 = (mwIndex*)mxGetData(prhs[19]);

		// Total number of non-zero values
		const uint64_t summa = (uint64_t)mxGetScalar(prhs[36]);

		// Is this a transmission image reconstruction (true) or emission image reconstruction (false)
		const bool attenuation_phase = (bool)mxGetScalar(prhs[37]);

		const uint32_t projector_type = (uint32_t)mxGetScalar(prhs[38]);

		// output sizes
		const mwSize N = static_cast<mwSize>(Nx) * static_cast<mwSize>(Ny) * static_cast<mwSize>(Nz);
		const mwSize nzmax = summa;
		const mwSize rows = pituus;

		// Create the MATLAB sparse matrix
		//mwSize joku = pituus;
		plhs[0] = mxCreateSparse(N, rows, nzmax, mxREAL);

		//plhs[2] = mxCreateNumericMatrix(nzmax, 1, mxDOUBLE_CLASS, mxREAL);
		//plhs[0] = mxCreateNumericMatrix(nzmax, 1, mxUINT64_CLASS, mxREAL);
		//plhs[3] = mxCreateNumericMatrix(loop_var_par + 1, 1, mxUINT64_CLASS, mxREAL);

		// Non-zero elements of the matrix
		double* elements = (double*)mxGetData(plhs[0]);

		// Row indices
		//mwIndex* indices = (mwIndex*)mxGetData(plhs[0]);
		mwIndex* indices = mxGetIr(plhs[0]);

		// Column indices
		mwIndex* lor = mxGetJc(plhs[0]);
		//mxSetJc(plhs[0], lor2);
		//mwIndex* lor = (mwIndex*)mxGetData(plhs[3]);

		for (uint32_t kk = 0; kk <= loop_var_par; kk++)
			lor[kk] = lor2[kk];

		if (attenuation_phase)
			plhs[1] = mxCreateNumericMatrix(pituus, 1, mxDOUBLE_CLASS, mxREAL);
		else
			plhs[1] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);

		double* ll = (double*)mxGetData(plhs[1]);

		//clock_t time = clock();
		std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

		if (projector_type == 2u) {

			if (nrhs != 45)
				mexErrMsgTxt("Incorrect number of input arguments. There has to be 45.");

			// Width of the TOR
			const double crystal_size = (double)mxGetScalar(prhs[39]);

			// Coordinates of the pixel centers in y-direction
			const double* x_center = (double*)mxGetData(prhs[40]);

			// Coordinates of the pixel centers in x-direction
			const double* y_center = (double*)mxGetData(prhs[41]);

			// Coordinates of the pixel centers in z-direction
			const double* z_center = (double*)mxGetData(prhs[42]);

			const double crystal_size_z = (double)mxGetScalar(prhs[43]);

			const int32_t dec_v = (int32_t)mxGetScalar(prhs[44]);

			// run the Orthogonal distance based ray tracer algorithm
			orth_siddon_precomputed(loop_var_par, size_x, zmax, indices, elements, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, x, y, z_det, 
				NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, lor1, lor2, xy_index, z_index, TotSinos, L, pseudos, 
				pRows, det_per_ring, raw, attenuation_phase, ll, crystal_size, crystal_size_z, y_center, x_center, z_center, dec_v);

			//time = clock() - time;
			
			std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

			if (verbose) {
				//mexPrintf("Orthogonal distance based ray tracer took %f seconds\n", ((float)time) / CLOCKS_PER_SEC);
				mexPrintf("Orthogonal distance based ray tracer took %f seconds\n", ((float)time_span.count()));
				mexEvalString("pause(.001);");
			}
		}
		else {

			// run the Improved Siddon's algorithm
			improved_siddon_precomputed(loop_var_par, size_x, zmax, indices, elements, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, x, y, z_det, 
				NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, lor1, lor2, xy_index, z_index, TotSinos, L, pseudos, 
				pRows, det_per_ring, raw, attenuation_phase, ll);

			//time = clock() - time;
			
			std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

			if (verbose) {
				//mexPrintf("Improved Siddon took %f seconds\n", ((float)time) / CLOCKS_PER_SEC);
				mexPrintf("Improved Siddon took %f seconds\n", ((float)time_span.count()));
				mexEvalString("pause(.001);");
			}
		}

	}
	else if (type == 1u) {

		if (nrhs < 41)
			mexErrMsgTxt("Too few input arguments.  There must be at least 41.");
		else if (nrhs > 47)
			mexErrMsgTxt("Too many input arguments.  There can be at most 47.");

		if (nlhs != 2)
			mexErrMsgTxt("Invalid number of output arguments. There has to be two.");

		const double epps = (double)mxGetScalar(prhs[35]);

		const double* Sino = (double*)mxGetData(prhs[36]);

		double* osem_apu = (double*)mxGetData(prhs[37]);

		const uint32_t projector_type = (uint32_t)mxGetScalar(prhs[38]);

		const size_t N = static_cast<size_t>(Nx) * static_cast<size_t>(Ny) * static_cast<size_t>(Nz);

		const bool no_norm = (bool)mxGetScalar(prhs[39]);

		const bool precompute = (bool)mxGetScalar(prhs[40]);

		plhs[0] = mxCreateNumericMatrix(N, 1, mxDOUBLE_CLASS, mxREAL);

		double* Summ = (double*)mxGetData(plhs[0]);

		plhs[1] = mxCreateNumericMatrix(N, 1, mxDOUBLE_CLASS, mxREAL);

		double* rhs = (double*)mxGetData(plhs[1]);

		//clock_t time = clock();
		//std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

		if (projector_type == 2u) {

			if (nrhs != 47)
				mexErrMsgTxt("Incorrect number of input arguments. There has to be 47.");

			// Width of the TOR
			const double crystal_size = (double)mxGetScalar(prhs[41]);

			// Coordinates of the pixel centers in y-direction
			const double* x_center = (double*)mxGetData(prhs[42]);

			// Coordinates of the pixel centers in x-direction
			const double* y_center = (double*)mxGetData(prhs[43]);

			// Coordinates of the pixel centers in z-direction
			const double* z_center = (double*)mxGetData(prhs[44]);

			const double crystal_size_z = (double)mxGetScalar(prhs[45]);

			const int32_t dec_v = (int32_t)mxGetScalar(prhs[46]);

			if (precompute) {
				sequential_orth_siddon(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, randoms, x, y, z_det, 
					NSlices, Nx, Ny, Nz, d, dz,	bx, by, bz, attenuation_correction, normalization, randoms_correction, lor1, xy_index, z_index, 
					TotSinos, epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, crystal_size, x_center, y_center, z_center, crystal_size_z, 
					no_norm, dec_v);
			}
			else {
				sequential_orth_siddon_no_precomp(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, randoms, norm_coef, 
					x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, randoms_correction, xy_index, z_index, 
					TotSinos, epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, crystal_size, x_center, y_center, z_center, crystal_size_z, 
					no_norm, dec_v);
			}

			//time = clock() - time;
		}
		else {

			//if (nrhs != 41)
			//	mexErrMsgTxt("Incorrect number of input arguments. There has to be 41.");

			if (precompute) {
				sequential_improved_siddon(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, randoms, norm_coef, x, y, 
					z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, randoms_correction, lor1, xy_index, z_index, 
					TotSinos, epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, no_norm);
			}
			else {
				sequential_improved_siddon_no_precompute(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, randoms, 
					norm_coef, x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, randoms_correction, xy_index, 
					z_index, TotSinos, epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, no_norm);
			}

			//time = clock() - time;

		}

		//if (verbose) {
		//	mexPrintf("Matrix free reconstruction took %f seconds\n", ((float)time) / CLOCKS_PER_SEC);
		//	mexEvalString("pause(.001);");
		//}
	}
	else if (type == 2u) {

		if (nrhs < 40)
			mexErrMsgTxt("Too few input arguments.  There must be at least 40.");
		else if (nrhs > 46)
			mexErrMsgTxt("Too many input arguments.  There can be at most 46.");

		if (nlhs != 3)
			mexErrMsgTxt("Invalid number of output arguments. There has to be three.");


		// How many elements are preallocated in memory
		const uint32_t ind_size = (uint32_t)mxGetScalar(prhs[35]);

		// Starting ring
		const uint32_t block1 = (uint32_t)mxGetScalar(prhs[36]);

		// End ring
		const uint32_t blocks = (uint32_t)mxGetScalar(prhs[37]);

		// Subset indices
		const uint32_t* index = (uint32_t*)mxGetData(prhs[38]);
		const uint32_t index_size = static_cast<uint32_t>(mxGetNumberOfElements(prhs[38]));

		const uint32_t projector_type = (uint32_t)mxGetScalar(prhs[39]);

		uint32_t loop_var_par = 1u;

		if (index_size > 1ULL && !raw) {
			loop_var_par = pituus;
		}
		else if (!raw) {
			loop_var_par = NSinos * size_x;
		}
		else {
			loop_var_par = static_cast<uint32_t>(numRows / 2ULL);
		}

		plhs[0] = mxCreateNumericMatrix(loop_var_par, 1, mxUINT16_CLASS, mxREAL);

		uint16_t* lor = (uint16_t*)mxGetData(plhs[0]);

		vector<uint32_t> indices;

		vector<double> elements;

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
				det_per_ring, blocks, block1, L, pseudos, pRows);
		}
		else if (projector_type == 2u) {

			if (nrhs != 46)
				mexErrMsgTxt("Incorrect number of input arguments. There has to be 46.");

			const double crystal_size = (double)mxGetScalar(prhs[40]);

			const double* x_center = (double*)mxGetData(prhs[41]);

			const double* y_center = (double*)mxGetData(prhs[42]);

			// Coordinates of the pixel centers in z-direction
			const double* z_center = (double*)mxGetData(prhs[43]);

			const double crystal_size_z = (double)mxGetScalar(prhs[44]);

			const int32_t dec_v = (int32_t)mxGetScalar(prhs[45]);

			// run the Orthogonal Siddon algorithm
			lj = orth_siddon_no_precompute(loop_var_par, size_x, zmax, TotSinos, indices, elements, lor, maxyy, maxxx, xx_vec, dy,
				yy_vec, atten, norm_coef, x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, index, attenuation_correction, normalization, raw, 
				det_per_ring, blocks, block1, L, pseudos, pRows, crystal_size, crystal_size_z, y_center, x_center, z_center, dec_v);
		}
		else if (projector_type == 0u) {

			if (nrhs < 43)
				mexErrMsgTxt("Too few input arguments.  There must be at least 43.");

			const double* iij = (double*)mxGetData(prhs[40]);

			const double* jji = (double*)mxGetData(prhs[41]);

			const double* kkj = (double*)mxGetData(prhs[42]);

			const vector<double> iij_vec(iij, iij + mxGetNumberOfElements(prhs[40]));

			const vector<double> jjk_vec(jji, jji + mxGetNumberOfElements(prhs[41]));

			const vector<double> kkj_vec(kkj, kkj + mxGetNumberOfElements(prhs[42]));

			// run the original Siddon's algorithm
			lj = original_siddon_no_precompute(loop_var_par, size_x, zmax, TotSinos, indices, elements, lor, maxyy, maxxx, xx_vec, dy,
				yy_vec, atten, norm_coef, x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, index, attenuation_correction, normalization, raw,
				det_per_ring, blocks, block1, L, pseudos, pRows, iij_vec, jjk_vec, kkj_vec);
		}

		const size_t outSize1 = static_cast<size_t>(lj) * 2ULL;
		const size_t outSize2 = 1ULL;

		// Create the MATLAB output vectors (row and column indices, elements)

		plhs[1] = mxCreateNumericMatrix(indices.size(), outSize2, mxUINT32_CLASS, mxREAL);

		uint32_t* outputMatrix2 = (uint32_t*)mxGetData(plhs[1]);

		copy(indices.begin(), indices.end(), outputMatrix2);
		indices.erase(indices.begin(), indices.end());
		indices.shrink_to_fit();

		plhs[2] = mxCreateNumericMatrix(elements.size(), outSize2, mxDOUBLE_CLASS, mxREAL);

		double* outputMatrix3 = (double*)mxGetData(plhs[2]);

		copy(elements.begin(), elements.end(), outputMatrix3);

		//time = clock() - time;

		std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

		if (verbose) {
			//mexPrintf("Function elapsed time is %f seconds\n", ((float)time) / CLOCKS_PER_SEC);
			mexPrintf("Function elapsed time is %f seconds\n", ((float)time_span.count()));
			mexEvalString("pause(.001);");
		}

	}
}