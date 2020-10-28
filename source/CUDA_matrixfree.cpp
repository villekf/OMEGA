/**************************************************************************
* Matrix free computations for OMEGA (Implementation 2, CUDA).
* This is the main file for the matrix-free computations in OMEGA. In this
* file the MATLAB variables are loaded and either the matrix-free 
* reconstructions are computed or the custom prior reconstruction.
* Unlike the non-CUDA/OpenCL versions, this one uses (32-bit) floats and 
* thus can be slightly inaccurate.
* Both functions use ArrayFire functions and thus require the installation
* of the ArrayFire library.
*
* Uses NVRTC (real-time compilation).
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
#include "AF_cuda_functions.hpp"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[]) {
	// Check for the number of input and output arguments
	if (nrhs < 40)
		mexErrMsgTxt("Too few input arguments.  There must be at least 40.");
	else if (nrhs > 75)
		mexErrMsgTxt("Too many input arguments.  There can be at most 75.");

	if (nlhs != 1)
		mexErrMsgTxt("Invalid number of output arguments.  There must be exactly one.");

	int ind = 0;
	// Load the input arguments
	// These are mainly same as with implementation 3 (check the comments in OpenCL_matrixfree_multi_gpu.cpp)
	const char* k_path = mxArrayToString(prhs[ind]);
	ind++;

	const uint32_t Ny = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	const uint32_t Nx = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	const uint32_t Nz = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	const float dx = (float)mxGetScalar(prhs[ind]);
	ind++;

	const float dz = (float)mxGetScalar(prhs[ind]);
	ind++;

	const float by = (float)mxGetScalar(prhs[ind]);
	ind++;

	const float bx = (float)mxGetScalar(prhs[ind]);
	ind++;

	const float bz = (float)mxGetScalar(prhs[ind]);
	ind++;

	const float* z_det = (float*)mxGetData(prhs[ind]);
	const size_t size_z = mxGetNumberOfElements(prhs[ind]);
	ind++;

	const float* x = (float*)mxGetData(prhs[ind]);
	const size_t size_of_x = mxGetNumberOfElements(prhs[ind]);
	ind++;

	const float* y = (float*)mxGetData(prhs[ind]);
	ind++;

	const float dy = (float)mxGetScalar(prhs[ind]);
	ind++;

	const float maxyy = (float)mxGetScalar(prhs[ind]);
	ind++;

	const float maxxx = (float)mxGetScalar(prhs[ind]);
	ind++;

	const uint32_t NSinos = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	const float NSlices = (float)mxGetScalar(prhs[ind]);
	ind++;

	const uint32_t size_x = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	const float zmax = (float)mxGetScalar(prhs[ind]);
	ind++;

	const uint32_t TotSinos = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	const bool verbose = (bool)mxGetScalar(prhs[ind]);
	ind++;

	const uint16_t* L = (uint16_t*)mxGetData(prhs[ind]);
	const size_t numRows = mxGetM(prhs[ind]);
	ind++;

	const uint32_t* pseudos = (uint32_t*)mxGetData(prhs[ind]);
	const uint32_t pRows = (uint32_t)mxGetNumberOfElements(prhs[ind]);
	ind++;

	const uint32_t det_per_ring = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	// Is TOF data used?
	const bool TOF = (bool)mxGetScalar(prhs[ind]);
	ind++;

	// Size of single TOF-subset
	const int64_t TOFSize = (int64_t)mxGetScalar(prhs[ind]);
	ind++;

	// Variance of the Gaussian TOF
	const float sigma_x = (float)mxGetScalar(prhs[ind]);
	ind++;

	// Centers of the TOF-bins
	const float* TOFCenter = (float*)mxGetData(prhs[ind]);
	ind++;

	// Index offset for TOF subsets
	const int64_t nBins = (int64_t)mxGetScalar(prhs[ind]);
	ind++;


	const int32_t dec = (int32_t)mxGetScalar(prhs[ind]);
	ind++;

	// The device used
	const uint32_t device = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	const uint8_t raw = (uint8_t)mxGetScalar(prhs[ind]);
	ind++;

	//af::setDevice(device);

	const char* fileName = mxArrayToString(prhs[ind]);
	ind++;

	const uint32_t type = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	const bool use_psf = (bool)mxGetScalar(prhs[ind]);
	ind++;

	// Directory to look for OpenCL headers
	const char* header_directory = mxArrayToString(prhs[ind]);
	ind++;

	const float bzb = bz + static_cast<float>(Nz) * dz;

	if (type < 2) {

		if (nrhs != 75)
			mexErrMsgTxt("Invalid number of input arguments. There must be 75.");

		const float* atten = (float*)mxGetData(prhs[ind]);
		const size_t size_atten = mxGetNumberOfElements(prhs[ind]);
		ind++;

		const float* norm = (float*)mxGetData(prhs[ind]);
		const size_t size_norm = mxGetNumberOfElements(prhs[ind]);
		ind++;

		int64_t* pituus = (int64_t*)mxGetData(prhs[ind]);
		ind++;

		const uint32_t attenuation_correction = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		const uint32_t normalization = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		const uint32_t Niter = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		const uint32_t subsets = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		//const uint8_t* rekot = (uint8_t*)mxGetData(prhs[37]);
		const size_t size_reko = mxGetNumberOfElements(prhs[ind]);
		ind++;

		const float epps = (float)mxGetScalar(prhs[ind]);
		ind++;

		const uint16_t* lor1 = (uint16_t*)mxGetData(prhs[ind]);
		ind++;

		const uint32_t* xy_index = (uint32_t*)mxGetData(prhs[ind]);

		size_t koko;
		if (raw)
			koko = numRows / 2;
		else
			koko = mxGetNumberOfElements(prhs[ind]);
		ind++;

		const uint16_t* z_index = (uint16_t*)mxGetData(prhs[ind]);
		ind++;

		// Is any OS-method used
		const bool osem_bool = (bool)mxGetScalar(prhs[ind]);
		ind++;

		const float tube_width = (float)mxGetScalar(prhs[ind]);
		ind++;

		const float crystal_size_z = (float)mxGetScalar(prhs[ind]);
		ind++;

		const float* x_center = (float*)mxGetData(prhs[ind]);
		const size_t size_center_x = mxGetNumberOfElements(prhs[ind]);
		ind++;

		const float* y_center = (float*)mxGetData(prhs[ind]);
		const size_t size_center_y = mxGetNumberOfElements(prhs[ind]);
		ind++;

		const float* z_center = (float*)mxGetData(prhs[ind]);
		const size_t size_center_z = mxGetNumberOfElements(prhs[ind]);
		ind++;

		// Randoms
		const mxArray* sc_ra = prhs[ind];
		ind++;

		// Randoms corrections
		const uint32_t randoms_correction = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		// The type of projector used (Siddon or orthogonal)
		const uint32_t projector_type = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		// If true, then the precomputed LOR voxels counts are used
		const bool precompute_var = (bool)mxGetScalar(prhs[ind]);
		ind++;

		// Number of rays in Siddon
		uint16_t n_rays = (uint16_t)mxGetScalar(prhs[ind]);
		ind++;

		// Number of rays in Siddon (axial)
		uint16_t n_rays3D = (uint16_t)mxGetScalar(prhs[ind]);
		ind++;

		// Crystal pitch in z-direction
		const float cr_pz = (float)mxGetScalar(prhs[ind]);
		ind++;

		const mxArray* options = prhs[ind];
		ind++;

		const bool saveIter = (bool)mxGetScalar(mxGetField(options, 0, "save_iter"));
		size_t Ni = 0ULL;
		if (saveIter)
			Ni = static_cast<size_t>(Niter);
		const size_t outSize = static_cast<size_t>(Nx) * static_cast<size_t>(Ny) * static_cast<size_t>(Nz);
		const size_t outSize2 = Ni + 1ULL;

		//// Implementation 2
		//if (type == 0) {

			// Cell array containing the measurements
		const mxArray* Sin = prhs[ind];
		ind++;

		// Number of time steps
		const uint32_t Nt = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		// Use 64-bit integer atomic functions if possible
		const bool use_64bit_atomics = (bool)mxGetScalar(prhs[ind]);
		ind++;

		// Number of OS-reconstruction algorithms (including priors)
		const uint32_t n_rekos = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		// Number of MLEM algorithms (including priors)
		const uint32_t n_rekos_mlem = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		// What type of reconstruction, needed for the OpenCL kernel
		// E.g. 2 means COSEM (different computations)
		const uint8_t* reko_type = (uint8_t*)mxGetData(prhs[ind]);
		ind++;


		const uint8_t* reko_type_mlem = (uint8_t*)mxGetData(prhs[ind]);
		ind++;

		// Global correction factor
		const float global_factor = (float)mxGetScalar(prhs[ind]);
		ind++;

		const float bmin = (float)mxGetScalar(prhs[ind]);
		ind++;

		const float bmax = (float)mxGetScalar(prhs[ind]);
		ind++;

		const float Vmax = (float)mxGetScalar(prhs[ind]);
		ind++;

		const float* V = (float*)mxGetData(prhs[ind]);
		const size_t size_V = mxGetNumberOfElements(prhs[ind]);
		ind++;

		const float* gaussian = (float*)mxGetData(prhs[ind]);
		const size_t size_gauss = mxGetNumberOfElements(prhs[ind]);
		ind++;

		// Create the output cell array
		mxArray* cell_array_ptr = mxCreateCellMatrix(size_reko + 1ULL, Nt);

		// Output dimensions
		const mwSize dim[4] = { static_cast<mwSize>(Nx), static_cast<mwSize>(Ny), static_cast<mwSize>(Nz), static_cast<mwSize>(outSize2) };

		reconstruction_AF_matrixfree(koko, lor1, z_det, x, y, Sin, sc_ra, Nx, Ny, Nz, Niter, options, dx, dy, dz, bx, by, bz, bzb, maxxx, maxyy, zmax,
			NSlices, pituus, xy_index, z_index, size_x, TotSinos, cell_array_ptr, dim, verbose, randoms_correction, attenuation_correction,
			normalization, atten, size_atten, norm, size_norm, subsets, epps, k_path, Nt, pseudos, det_per_ring, pRows, L, raw, size_z, osem_bool,
			fileName, use_psf, tube_width, crystal_size_z, x_center, y_center, z_center, size_center_x, size_center_y, size_of_x, size_center_z,
			projector_type, header_directory, precompute_var, device, dec, n_rays, n_rays3D, cr_pz, use_64bit_atomics, n_rekos, n_rekos_mlem, reko_type, reko_type_mlem,
			global_factor, bmin, bmax, Vmax, V, size_V, gaussian, size_gauss, saveIter, TOF, TOFSize, sigma_x, TOFCenter, nBins);
		//}


		plhs[0] = cell_array_ptr;
	}
	// Compute the number of voxels each LOR traverses (using AF device)
	else if (type == 2) {

		if (nrhs != 40)
			mexErrMsgTxt("Invalid number of input arguments. There must be 40.");

		// Starting block
		const uint32_t block1 = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		// Ending block
		const uint32_t blocks = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		// Number of sinograms
		const uint32_t NSinos = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		// Total number of sinograms
		const uint16_t TotSinos = (uint16_t)mxGetScalar(prhs[ind]);
		ind++;

		size_t loop_var_par = 1ULL;

		if (raw)
			loop_var_par = numRows / 2ULL;
		else
			loop_var_par = NSinos * size_x;

		plhs[0] = mxCreateNumericMatrix(loop_var_par, 1, mxUINT16_CLASS, mxREAL);

		uint16_t* lor = (uint16_t*)mxGetData(plhs[0]);

		//find_LORs(lor, z_det, x, y, Nx, Ny, Nz, dx, dy, dz, bx, by, bz, bzb, maxxx, maxyy, zmax, NSlices, size_x, TotSinos, verbose, loop_var_par, 
		//	k_path, pseudos, det_per_ring, pRows, L, raw, size_z, fileName, device, size_of_x, use_psf, header_directory);

	}

	// Clear ArrayFire memory
	af::deviceGC();

	return;
}