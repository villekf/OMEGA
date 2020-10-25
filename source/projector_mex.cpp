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
	if (nrhs < 34)
		mexErrMsgTxt("Too few input arguments. There must be at least 34.");
	else if (nrhs > 61)
		mexErrMsgTxt("Too many input arguments. There can be at most 61.");

	if (nlhs > 3 || nlhs < 1)
		mexErrMsgTxt("Invalid number of output arguments. There can be at most three.");

	int ind = 0;
	// Load the input arguments
	// Image size in y-direction
	const uint32_t Ny = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	// Image size in x-direction
	const uint32_t Nx = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	// Image size in z-direction
	const uint32_t Nz = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	// Distance between adjacent pixels in x-direction
	const double d = (double)mxGetScalar(prhs[ind]);
	ind++;

	// Distance between adjacent pixels in z-direction
	const double dz = (double)mxGetScalar(prhs[ind]);
	ind++;

	// Distance of the pixel grid from the origin in y-direction
	const double by = (double)mxGetScalar(prhs[ind]);
	ind++;

	// Distance of the pixel grid from the origin in x-direction
	const double bx = (double)mxGetScalar(prhs[ind]);
	ind++;

	// Distance of the pixel grid from the origin in z-direction
	const double bz = (double)mxGetScalar(prhs[ind]);
	ind++;

	// Coordinates of the detectors in z-direction
	const double* z_det = (double*)mxGetData(prhs[ind]);

	const vector<double> z_det_vec(z_det, z_det + mxGetNumberOfElements(prhs[ind]));
	ind++;

	// Coordinates of the detectors in x-direction
	const double* x = (double*)mxGetData(prhs[ind]);
	ind++;

	// Coordinates of the detectors in y-direction
	const double* y = (double*)mxGetData(prhs[ind]);
	ind++;

	// Distance between adjacent pixels in y-direction
	const double dy = (double)mxGetScalar(prhs[ind]);
	ind++;

	// Coordinates of the pixel planes (boundaries of the pixels) in y-direction
	const double* yy = (double*)mxGetData(prhs[ind]);

	// from array to std::vector
	const vector<double> yy_vec(yy, yy + mxGetNumberOfElements(prhs[ind]));
	ind++;

	// Coordinates of the pixel planes (boundaries of the pixels) in x-direction
	const double* xx = (double*)mxGetData(prhs[ind]);

	const vector<double> xx_vec(xx, xx + mxGetNumberOfElements(prhs[ind]));
	ind++;

	// Number of sinograms used
	const uint32_t NSinos = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	// Number of slices included
	const uint32_t NSlices = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	// Number of detector indices
	const uint32_t size_x = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	// Maximum value of the z-direction detector coordinates
	const double zmax = (double)mxGetScalar(prhs[ind]);
	ind++;

	// attenuation values (attenuation images)
	const double* atten = (double*)mxGetData(prhs[ind]);
	ind++;

	// Normalization coefficients
	const double* norm_coef = (double*)mxGetData(prhs[ind]);
	ind++;

	// Randoms
	const double* randoms = (double*)mxGetData(prhs[ind]);
	ind++;

	// Number of measurements/LORs
	const int64_t pituus = (int64_t)mxGetScalar(prhs[ind]);
	ind++;

	// Is the attenuation correction included
	const bool attenuation_correction = (bool)mxGetScalar(prhs[ind]);
	ind++;

	// Is the normalization correction included
	const bool normalization = (bool)mxGetScalar(prhs[ind]);
	ind++;

	// Is the randoms/scatter correction included
	const bool randoms_correction = (bool)mxGetScalar(prhs[ind]);
	ind++;

	// Is scatter correction included as multiplication (system matrix)
	const bool scatter = (bool)mxGetScalar(prhs[ind]);
	ind++;

	// Scatter data
	const double* scatter_coef = (double*)mxGetData(prhs[ind]);
	ind++;

	// Global correction factor
	const double global_factor = (double)mxGetScalar(prhs[ind]);
	ind++;

	// Number of voxels the current LOR/ray traverses, precomputed data ONLY
	const uint16_t* lor1 = (uint16_t*)mxGetData(prhs[ind]);
	ind++;

	// For sinogram data, the indices of the detectors corresponding to the current sinogram bin
	const uint32_t* xy_index = (uint32_t*)mxGetData(prhs[ind]);
	ind++;

	// Same as above, but for z-direction
	const uint16_t* z_index = (uint16_t*)mxGetData(prhs[ind]);
	ind++;

	// Total number of sinograms
	const uint32_t TotSinos = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	// Detector pair numbers, for raw list-mode data
	const uint16_t* L = (uint16_t*)mxGetData(prhs[ind]);
	const size_t numRows = mxGetNumberOfElements(prhs[ind]);
	ind++;

	// Location (ring numbers) of pseudo rings, if present
	const uint32_t* pseudos = (uint32_t*)mxGetData(prhs[ind]);
	const uint32_t pRows = (uint32_t)mxGetM(prhs[ind]);
	ind++;

	// Number of detectors per ring
	const uint32_t det_per_ring = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	// Is TOF data used?
	const bool TOF = (bool)mxGetScalar(prhs[ind]);
	ind++;

	// Size of single TOF-subset
	const int64_t TOFSize = (int64_t)mxGetScalar(prhs[ind]);
	ind++;

	// Variance of the Gaussian TOF
	const double sigma_x = (double)mxGetScalar(prhs[ind]);
	ind++;

	// Centers of the TOF-bins
	const double* TOFCenter = (double*)mxGetData(prhs[ind]);
	ind++;

	// Index offset for TOF subsets
	const int64_t nBins = (int64_t)mxGetScalar(prhs[ind]);
	ind++;


	const uint32_t dec_v = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	// Are status messages displayed
	const bool verbose = (bool)mxGetScalar(prhs[ind]);
	ind++;

	// Is raw list-mode data used
	const bool raw = (bool)mxGetScalar(prhs[ind]);
	ind++;

	const uint32_t type = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	//44

	// Number of measurements/LORs
	int64_t loop_var_par = pituus;

	// The maximum elements of the pixel space in both x- and y-directions
	const double maxyy = yy_vec.back();
	const double maxxx = xx_vec.back();

	// Implementation 1, with precomputed_lor = true
	if (type == 0u) {

		if (nrhs < 48)
			mexErrMsgTxt("Too few input arguments.  There must be at least 48.");
		else if (nrhs > 57)
			mexErrMsgTxt("Too many input arguments.  There can be at most 57.");

		if (nlhs != 2)
			mexErrMsgTxt("Invalid number of output arguments. There has to be two.");

		// Total number of voxels traversed by the LORs at the specific LOR
		// e.g. if the first LOR traverses through 10 voxels then at the second lor the value is 10
		const uint64_t* lor2 = (uint64_t*)mxGetData(prhs[ind]);
		ind++;

		// Total number of non-zero values
		const uint64_t summa = (uint64_t)mxGetScalar(prhs[ind]);
		ind++;

		// Is this a transmission image reconstruction (true) or emission image reconstruction (false)
		const bool attenuation_phase = (bool)mxGetScalar(prhs[ind]);
		ind++;

		const uint32_t projector_type = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		// output sizes
		const mwSize N = static_cast<mwSize>(Nx) * static_cast<mwSize>(Ny) * static_cast<mwSize>(Nz);
		const mwSize nzmax = summa;
		const mwSize rows = pituus;

		// Create the MATLAB sparse matrix
		plhs[0] = mxCreateSparse(N, rows, nzmax, mxREAL);

		// Non-zero elements of the matrix
		double* elements = (double*)mxGetData(plhs[0]);

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

		double* ll = (double*)mxGetData(plhs[1]);

		// Timing
		std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

		if (projector_type == 2u) {

			if (nrhs != 57)
				mexErrMsgTxt("Incorrect number of input arguments. There has to be 57.");

			// Width of the TOR
			const double crystal_size = (double)mxGetScalar(prhs[ind]);
			ind++;

			// Coordinates of the pixel centers in y-direction
			double* x_center = (double*)mxGetData(prhs[ind]);
			ind++;

			// Coordinates of the pixel centers in x-direction
			double* y_center = (double*)mxGetData(prhs[ind]);
			ind++;

			// Coordinates of the pixel centers in z-direction
			const double* z_center = (double*)mxGetData(prhs[ind]);
			ind++;

			const double crystal_size_z = (double)mxGetScalar(prhs[ind]);
			ind++;

			// run the Orthogonal distance based ray tracer algorithm, precomputed_lor = true
			orth_siddon_precomputed(loop_var_par, size_x, zmax, indices, elements, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, x, y, z_det,
				NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, lor1, lor2, xy_index, z_index, TotSinos, L, pseudos,
				pRows, det_per_ring, raw, attenuation_phase, ll, crystal_size, crystal_size_z, y_center, x_center, z_center, global_factor, scatter, scatter_coef);

			std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

			if (verbose) {
				mexPrintf("Orthogonal distance based ray tracer took %f seconds\n", ((float)time_span.count()));
				mexEvalString("pause(.001);");
			}
		}
		else if (projector_type == 1u) {

			// run the Improved Siddon's algorithm, precomputed_lor = true
			improved_siddon_precomputed(loop_var_par, size_x, zmax, indices, elements, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, x, y, z_det,
				NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, lor1, lor2, xy_index, z_index, TotSinos, L, pseudos,
				pRows, det_per_ring, raw, attenuation_phase, ll, global_factor, scatter, scatter_coef);

			std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

			if (verbose) {
				mexPrintf("Improved Siddon took %f seconds\n", ((float)time_span.count()));
				mexEvalString("pause(.001);");
			}
		}
		else if ((projector_type == 3u)) {

			if (nrhs != 57)
				mexErrMsgTxt("Incorrect number of input arguments. There has to be 57.");

			// Width of the TOR
			const double crystal_size = (double)mxGetScalar(prhs[ind]);
			ind++;

			// Coordinates of the pixel centers in y-direction
			double* x_center = (double*)mxGetData(prhs[ind]);
			ind++;

			// Coordinates of the pixel centers in x-direction
			double* y_center = (double*)mxGetData(prhs[ind]);
			ind++;

			// Coordinates of the pixel centers in z-direction
			const double* z_center = (double*)mxGetData(prhs[ind]);
			ind++;

			const double crystal_size_z = (double)mxGetScalar(prhs[ind]);
			ind++;

			const double bmin = (double)mxGetScalar(prhs[ind]);
			ind++;

			const double bmax = (double)mxGetScalar(prhs[ind]);
			ind++;

			const double Vmax = (double)mxGetScalar(prhs[ind]);
			ind++;

			const double* V = (double*)mxGetData(prhs[ind]);
			ind++;

			// run the Orthogonal distance based ray tracer algorithm, precomputed_lor = true
			vol_siddon_precomputed(loop_var_par, size_x, zmax, indices, elements, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, x, y, z_det,
				NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, lor1, lor2, xy_index, z_index, TotSinos, L, pseudos,
				pRows, det_per_ring, raw, attenuation_phase, ll, crystal_size, crystal_size_z, y_center, x_center, z_center, global_factor, bmin, 
				bmax, Vmax, V, scatter, scatter_coef);

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
		else if (nrhs > 61)
			mexErrMsgTxt("Too many input arguments.  There can be at most 61.");

		if (nlhs != 2)
			mexErrMsgTxt("Invalid number of output arguments. There has to be two.");

		// Small constant to prevent division by zero
		const double epps = (double)mxGetScalar(prhs[ind]);
		ind++;

		// Measurement data
		const double* Sino = (double*)mxGetData(prhs[ind]);
		ind++;

		// Current estimates
		double* osem_apu = (double*)mxGetData(prhs[ind]);
		ind++;

		// Projector used
		const uint32_t projector_type = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		const size_t N = static_cast<size_t>(Nx) * static_cast<size_t>(Ny) * static_cast<size_t>(Nz);

		// If 1, do not compute the normalization constant anymore
		const bool no_norm = (bool)mxGetScalar(prhs[ind]);
		ind++;

		// precomputed_lor = false
		const bool precompute = (bool)mxGetScalar(prhs[ind]);
		ind++;

		const uint8_t fp = (uint8_t)mxGetScalar(prhs[ind]);
		ind++;

		const bool list_mode_format = (bool)mxGetScalar(prhs[ind]);
		ind++;

		plhs[0] = mxCreateNumericMatrix(N, 1, mxDOUBLE_CLASS, mxREAL);

		// Normalization constants
		double* Summ = (double*)mxGetData(plhs[0]);

		if (fp == 1)
			plhs[1] = mxCreateNumericMatrix(pituus, 1, mxDOUBLE_CLASS, mxREAL);
		else
			plhs[1] = mxCreateNumericMatrix(N, 1, mxDOUBLE_CLASS, mxREAL);

		// Right-hand side
		double* rhs = (double*)mxGetData(plhs[1]);

		//clock_t time = clock();
		//std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

		// Orthogonal
		if (projector_type == 2u) {

			if (nrhs < 56)
				mexErrMsgTxt("Incorrect number of input arguments. There has to be 56.");

			// Width of the strip in 2D case
			const double crystal_size = (double)mxGetScalar(prhs[ind]);
			ind++;

			// Coordinates of the pixel centers in y-direction
			double* x_center = (double*)mxGetData(prhs[ind]);
			ind++;

			// Coordinates of the pixel centers in x-direction
			double* y_center = (double*)mxGetData(prhs[ind]);
			ind++;

			// Coordinates of the pixel centers in z-direction
			const double* z_center = (double*)mxGetData(prhs[ind]);
			ind++;

			// Width of the TOR in 3D case
			const double crystal_size_z = (double)mxGetScalar(prhs[ind]);
			ind++;

			if (precompute) {
				sequential_orth_siddon(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, randoms, x, y, z_det, 
					NSlices, Nx, Ny, Nz, d, dz,	bx, by, bz, attenuation_correction, normalization, randoms_correction, lor1, xy_index, z_index, 
					TotSinos, epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, crystal_size, x_center, y_center, z_center, crystal_size_z, 
					no_norm, dec_v, global_factor, fp, scatter, scatter_coef);
			}
			else {
				sequential_orth_siddon_no_precomp(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, randoms,
					x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, randoms_correction, xy_index, z_index, 
					TotSinos, epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, crystal_size, x_center, y_center, z_center, crystal_size_z, 
					no_norm, dec_v, global_factor, fp, list_mode_format, scatter, scatter_coef);
			}
		}
		// Improved Siddon
		else if (projector_type == 1u) {

			if (nrhs < 54)
				mexErrMsgTxt("Incorrect number of input arguments. There has to be 54.");

			// Number of rays in Siddon (transaxial)
			uint16_t n_rays = (uint16_t)mxGetScalar(prhs[ind]);
			ind++;

			// Number of rays in Siddon (axial)
			uint16_t n_rays3D = (uint16_t)mxGetScalar(prhs[ind]);
			ind++;

			// Crystal pitch in z-direction (for multi-ray)
			const double cr_pz = (double)mxGetScalar(prhs[ind]);
			ind++;

			if (precompute) {
				sequential_improved_siddon(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, randoms, x, y,
					z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, randoms_correction, lor1, xy_index, z_index, 
					TotSinos, epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, no_norm, global_factor, fp, scatter, scatter_coef, TOF, TOFSize, 
					sigma_x, TOFCenter, nBins, dec_v);
			}
			else {
				sequential_improved_siddon_no_precompute(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, randoms, x, y,
					z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, randoms_correction, xy_index, z_index, TotSinos,
					epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, cr_pz, no_norm, n_rays, n_rays3D, global_factor, fp, list_mode_format, 
					scatter, scatter_coef, TOF, TOFSize, sigma_x, TOFCenter, nBins, dec_v);
			}


		}
		else if ((projector_type == 3u)) {
			if (nrhs < 59)
				mexErrMsgTxt("Incorrect number of input arguments. There has to be 59.");

			// Coordinates of the pixel centers in y-direction
			double* x_center = (double*)mxGetData(prhs[ind]);
			ind++;

			// Coordinates of the pixel centers in x-direction
			double* y_center = (double*)mxGetData(prhs[ind]);
			ind++;

			// Coordinates of the pixel centers in z-direction
			const double* z_center = (double*)mxGetData(prhs[ind]);
			ind++;

			// Width of the TOR in 3D case
			const double bmin = (double)mxGetScalar(prhs[ind]);
			ind++;

			// Width of the TOR in 3D case
			const double bmax = (double)mxGetScalar(prhs[ind]);
			ind++;

			// Width of the TOR in 3D case
			const double Vmax = (double)mxGetScalar(prhs[ind]);
			ind++;

			// Coordinates of the pixel centers in y-direction
			const double* V = (double*)mxGetData(prhs[ind]);
			ind++;

			if (precompute) {
				sequential_volume_siddon(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, randoms, x, y, z_det,
					NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, randoms_correction, lor1, xy_index, z_index,
					TotSinos, epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, Vmax, x_center, y_center, z_center, bmin, bmax, V,
					no_norm, dec_v, global_factor, fp, scatter, scatter_coef);
			}
			else {
				sequential_volume_siddon_no_precomp(loop_var_par, size_x, zmax, Summ, rhs, maxyy, maxxx, xx_vec, dy, yy_vec, atten, norm_coef, randoms,
					x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, normalization, randoms_correction, xy_index, z_index,
					TotSinos, epps, Sino, osem_apu, L, pseudos, pRows, det_per_ring, raw, Vmax, x_center, y_center, z_center, bmin, bmax, V,
					no_norm, dec_v, global_factor, fp, list_mode_format, scatter, scatter_coef);
			}
		}
	}
	// Implementation 1, precomputed_lor = false
	else if (type == 2u) {

		if (nrhs < 48)
			mexErrMsgTxt("Too few input arguments.  There must be at least 48.");
		else if (nrhs > 54)
			mexErrMsgTxt("Too many input arguments.  There can be at most 54.");

		if (nlhs != 3)
			mexErrMsgTxt("Invalid number of output arguments. There has to be three.");


		// How many elements are preallocated in memory
		const uint32_t ind_size = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		// Starting ring
		const uint32_t block1 = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		// End ring
		const uint32_t blocks = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		// Subset indices
		const uint32_t* index = (uint32_t*)mxGetData(prhs[ind]);
		const uint32_t index_size = static_cast<uint32_t>(mxGetNumberOfElements(prhs[ind]));
		ind++;

		// Projector
		const uint32_t projector_type = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		// Number of LORs
		if (index_size > 1ULL && !raw) {
			loop_var_par = index_size;
		}
		else if (!raw) {
			loop_var_par = static_cast<size_t>(NSinos) * static_cast<size_t>(size_x);
		}
		else {
			loop_var_par = numRows / 2ULL;
		}

		plhs[0] = mxCreateNumericMatrix(loop_var_par, 1, mxUINT16_CLASS, mxREAL);

		uint16_t* lor = (uint16_t*)mxGetData(plhs[0]);

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
				det_per_ring, blocks, block1, L, pseudos, pRows, global_factor, scatter, scatter_coef);
		}
		else if (projector_type == 2u) {

			if (nrhs != 54)
				mexErrMsgTxt("Incorrect number of input arguments. There has to be 54.");

			const double crystal_size = (double)mxGetScalar(prhs[ind]);
			ind++;

			double* x_center = (double*)mxGetData(prhs[ind]);
			ind++;

			double* y_center = (double*)mxGetData(prhs[ind]);
			ind++;

			// Coordinates of the pixel centers in z-direction
			const double* z_center = (double*)mxGetData(prhs[ind]);
			ind++;

			const double crystal_size_z = (double)mxGetScalar(prhs[ind]);
			ind++;

			// run the Orthogonal Siddon algorithm
			//lj = orth_siddon_no_precompute(loop_var_par, size_x, zmax, TotSinos, indices, elements, lor, maxyy, maxxx, xx_vec, dy,
			//	yy_vec, atten, norm_coef, x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, index, attenuation_correction, normalization, raw, 
			//	det_per_ring, blocks, block1, L, pseudos, pRows, crystal_size, crystal_size_z, y_center, x_center, z_center, dec_v);
		}
		// Original Siddon's ray tracer
		else if (projector_type == 0u) {

			if (nrhs < 51)
				mexErrMsgTxt("Too few input arguments.  There must be at least 51.");

			// Voxel numbers in x-direction
			const double* iij = (double*)mxGetData(prhs[ind]);
			const vector<double> iij_vec(iij, iij + mxGetNumberOfElements(prhs[ind]));
			ind++;

			// Voxel numbers in y-direction
			const double* jji = (double*)mxGetData(prhs[ind]);
			const vector<double> jjk_vec(jji, jji + mxGetNumberOfElements(prhs[ind]));
			ind++;

			// Voxel numbers in z-direction
			const double* kkj = (double*)mxGetData(prhs[ind]);
			const vector<double> kkj_vec(kkj, kkj + mxGetNumberOfElements(prhs[ind]));
			ind++;


			// run the original Siddon's algorithm
			lj = original_siddon_no_precompute(loop_var_par, size_x, zmax, TotSinos, indices, elements, lor, maxyy, maxxx, xx_vec, dy,
				yy_vec, atten, norm_coef, x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, index, attenuation_correction, normalization, raw,
				det_per_ring, blocks, block1, L, pseudos, pRows, iij_vec, jjk_vec, kkj_vec, global_factor, scatter, scatter_coef);
		}

		const size_t outSize1 = static_cast<size_t>(lj) * 2ULL;
		const size_t outSize2 = 1ULL;

		// Create the MATLAB output vectors (row and column indices, elements)

		plhs[1] = mxCreateNumericMatrix(indices.size(), outSize2, mxUINT32_CLASS, mxREAL);

		uint32_t* outputMatrix2 = (uint32_t*)mxGetData(plhs[1]);

		std::copy(indices.begin(), indices.end(), outputMatrix2);
		indices.erase(indices.begin(), indices.end());
		indices.shrink_to_fit();

		plhs[2] = mxCreateNumericMatrix(elements.size(), outSize2, mxDOUBLE_CLASS, mxREAL);

		double* outputMatrix3 = (double*)mxGetData(plhs[2]);

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
		if (nrhs < 57)
			mexErrMsgTxt("Too few input arguments.  There must be at least 57.");
		else if (nrhs > 57)
			mexErrMsgTxt("Too many input arguments.  There can be at most 57.");

		if (nlhs != 3)
			mexErrMsgTxt("Invalid number of output arguments. There has to be three.");


		// Starting ring
		const uint32_t block1 = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		// End ring
		const uint32_t blocks = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		// Projector
		const uint32_t projector_type = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		const double crystal_size = (double)mxGetScalar(prhs[ind]);
		ind++;

		double* x_center = (double*)mxGetData(prhs[ind]);
		ind++;

		double* y_center = (double*)mxGetData(prhs[ind]);
		ind++;

		// Coordinates of the pixel centers in z-direction
		const double* z_center = (double*)mxGetData(prhs[ind]);
		ind++;

		const double crystal_size_z = (double)mxGetScalar(prhs[ind]);
		ind++;

		// Width of the TOR in 3D case
		const double bmin = (double)mxGetScalar(prhs[ind]);
		ind++;

		// Width of the TOR in 3D case
		const double bmax = (double)mxGetScalar(prhs[ind]);
		ind++;

		// Width of the TOR in 3D case
		const double Vmax = (double)mxGetScalar(prhs[ind]);
		ind++;

		// Coordinates of the pixel centers in y-direction
		const double* V = (double*)mxGetData(prhs[ind]);
		ind++;

		const uint32_t tyyppi = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

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

		lor = (uint16_t*)mxGetData(plhs[0]);
		lor_orth = (uint16_t*)mxGetData(plhs[1]);
		lor_vol = (uint16_t*)mxGetData(plhs[2]);

		improved_siddon_precomputation_phase(loop_var_par, size_x, zmax, TotSinos, lor, maxyy, maxxx, xx_vec, z_det_vec, dy, yy_vec, x, y, z_det, NSlices, Nx, Ny, Nz,
			d, dz, bx, by, bz, block1, blocks, L, pseudos, raw, pRows, det_per_ring, tyyppi, lor_orth, lor_vol, crystal_size, crystal_size_z, x_center, y_center, z_center, 
			bmin, bmax, Vmax, V);

	}
}