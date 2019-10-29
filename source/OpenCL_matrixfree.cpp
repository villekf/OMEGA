/**************************************************************************
* Matrix free computations for OMEGA (Implementation 2).
* This is the main file for the matrix-free computations in OMEGA. In this
* file the MATLAB variables are loaded and either the matrix-free 
* reconstructions are computed or the custom prior reconstruction.
* Unlike the non-OpenCL versions, this one uses (32-bit) floats and thus
* can be slightly inaccurate.
* Both functions use ArrayFire functions and thus require the installation
* of the ArrayFire library.
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
#include "functions.hpp"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[]) {
	// Check for the number of input and output arguments
	if (nrhs < 34)
		mexErrMsgTxt("Too few input arguments.  There must be at least 34.");
	else if (nrhs > 62)
		mexErrMsgTxt("Too many input arguments.  There can be at most 62.");

	if (nlhs != 1)
		mexErrMsgTxt("Invalid number of output arguments.  There must be exactly one.");

	// Load the input arguments
	// These are mainly same as with implementation 3 (check the comments in OpenCL_matrixfree_multi_gpu.cpp)
	const char *k_path = mxArrayToString(prhs[0]);

	const uint32_t Ny = (uint32_t)mxGetScalar(prhs[1]);

	const uint32_t Nx = (uint32_t)mxGetScalar(prhs[2]);

	const uint32_t Nz = (uint32_t)mxGetScalar(prhs[3]);

	const float dx = (float)mxGetScalar(prhs[4]);

	const float dz = (float)mxGetScalar(prhs[5]);

	const float by = (float)mxGetScalar(prhs[6]);

	const float bx = (float)mxGetScalar(prhs[7]);

	const float bz = (float)mxGetScalar(prhs[8]);

	const float *z_det = (float*)mxGetData(prhs[9]);
	const size_t size_z = mxGetNumberOfElements(prhs[9]);

	const float *x = (float*)mxGetData(prhs[10]);
	const size_t size_of_x = mxGetNumberOfElements(prhs[10]);

	const float *y = (float*)mxGetData(prhs[11]);

	const float dy = (float)mxGetScalar(prhs[12]);

	const float maxyy = (float)mxGetScalar(prhs[13]);

	const float maxxx = (float)mxGetScalar(prhs[14]);

	const uint32_t NSinos = (uint32_t)mxGetScalar(prhs[15]);

	const float NSlices = (float)mxGetScalar(prhs[16]);

	const uint32_t size_x = (uint32_t)mxGetScalar(prhs[17]);

	const float zmax = (float)mxGetScalar(prhs[18]);

	const uint32_t TotSinos = (uint32_t)mxGetScalar(prhs[19]);

	const bool verbose = (bool)mxGetScalar(prhs[20]);

	const uint16_t *L = (uint16_t*)mxGetData(prhs[21]);
	const size_t numRows = mxGetM(prhs[21]);

	const uint32_t *pseudos = (uint32_t*)mxGetData(prhs[22]);
	const uint32_t pRows = (uint32_t)mxGetNumberOfElements(prhs[22]);

	const uint32_t det_per_ring = (uint32_t)mxGetScalar(prhs[23]);

	// The device used
	const uint32_t device = (uint32_t)mxGetScalar(prhs[24]);

	const uint8_t raw = (uint8_t)mxGetScalar(prhs[25]);

	af::setDevice(device);

	const char *fileName = mxArrayToString(prhs[26]);

	const uint32_t type = (uint32_t)mxGetScalar(prhs[27]);

	const bool force_build = (bool)mxGetScalar(prhs[28]);

	// Directory to look for OpenCL headers
	const char* header_directory = mxArrayToString(prhs[29]);

	const float bzb = bz + static_cast<float>(Nz) * dz;

	if (type < 2) {

		const float *atten = (float*)mxGetData(prhs[30]);
		const size_t size_atten = mxGetNumberOfElements(prhs[30]);

		const float* norm = (float*)mxGetData(prhs[31]);
		const size_t size_norm = mxGetNumberOfElements(prhs[31]);

		uint32_t *pituus = (uint32_t*)mxGetData(prhs[32]);

		const uint32_t attenuation_correction = (uint32_t)mxGetScalar(prhs[33]);

		const uint32_t normalization = (uint32_t)mxGetScalar(prhs[34]);

		const uint32_t Niter = (uint32_t)mxGetScalar(prhs[35]);

		const uint32_t subsets = (uint32_t)mxGetScalar(prhs[36]);

		//const uint8_t* rekot = (uint8_t*)mxGetData(prhs[37]);
		const size_t size_reko = mxGetNumberOfElements(prhs[37]);

		const float epps = (float)mxGetScalar(prhs[38]);

		const uint16_t* lor1 = (uint16_t*)mxGetData(prhs[39]);

		const uint32_t* xy_index = (uint32_t*)mxGetData(prhs[40]);

		const uint16_t* z_index = (uint16_t*)mxGetData(prhs[41]);

		// Is any OS-method used
		const bool osem_bool = (bool)mxGetScalar(prhs[42]);

		const float tube_width = (float)mxGetScalar(prhs[43]);

		const float crystal_size_z = (float)mxGetScalar(prhs[44]);

		const float *x_center = (float*)mxGetData(prhs[45]);
		const size_t size_center_x = mxGetNumberOfElements(prhs[45]);

		const float *y_center = (float*)mxGetData(prhs[46]);
		const size_t size_center_y = mxGetNumberOfElements(prhs[46]);

		const float* z_center = (float*)mxGetData(prhs[47]);
		const size_t size_center_z = mxGetNumberOfElements(prhs[47]);

		// Randoms
		const mxArray* sc_ra = prhs[48];

		// Randoms corrections
		const uint32_t randoms_correction = (uint32_t)mxGetScalar(prhs[49]);

		// The type of projector used (Siddon or orthogonal)
		const uint32_t projector_type = (uint32_t)mxGetScalar(prhs[50]);

		// If true, then the precomputed LOR voxels counts are used
		const bool precompute_var = (bool)mxGetScalar(prhs[51]);

		// Accuracy factor in orthogonal distance
		const int32_t dec = (int32_t)mxGetScalar(prhs[52]);

		// Number of rays in Siddon
		uint16_t n_rays = (uint16_t)mxGetScalar(prhs[53]);

		// Crystal pitch in z-direction
		const float cr_pz = (double)mxGetScalar(prhs[54]);

		size_t koko;
		if (raw)
			koko = numRows;
		else
			koko = mxGetM(prhs[40]);

		const size_t outSize = Nx * Ny * Nz;
		const size_t outSize2 = Niter + 1u;
		mxArray *cell_array_ptr;

		// Implementation 2
		if (type == 0) {

			if (nrhs != 62)
				mexErrMsgTxt("Invalid number of input arguments. There must be 62.");

			// Cell array containing the measurements
			const mxArray* Sin = prhs[56];

			// Number of time steps
			const uint32_t Nt = (uint32_t)mxGetScalar(prhs[57]);

			// Use 64-bit integer atomic functions if possible
			const bool use_64bit_atomics = (bool)mxGetScalar(prhs[58]);

			// Number of OS-reconstruction algorithms (including priors)
			const uint32_t n_rekos = (uint32_t)mxGetScalar(prhs[59]);

			// Number of MLEM algorithms (including priors)
			const uint32_t n_rekos_mlem = (uint32_t)mxGetScalar(prhs[60]);

			// What type of reconstruction, needed for the OpenCL kernel
			// E.g. 2 means COSEM (different computations)
			const uint8_t* reko_type = (uint8_t*)mxGetData(prhs[61]);

			// Create the output cell array
			cell_array_ptr = mxCreateCellMatrix(size_reko + 1ULL, Nt);

			// Output dimensions
			const mwSize dim[4] = { static_cast<mwSize>(Nx), static_cast<mwSize>(Ny), static_cast<mwSize>(Nz), static_cast<mwSize>(Niter + 1u) };

			reconstruction_AF_matrixfree(koko, lor1, z_det, x, y, Sin, sc_ra, Nx, Ny, Nz, Niter, prhs[55], dx, dy, dz, bx, by, bz, bzb, maxxx, maxyy, zmax, 
				NSlices, pituus, xy_index, z_index, size_x, TotSinos, cell_array_ptr, dim, verbose, randoms_correction, attenuation_correction, 
				normalization, atten, size_atten, norm, size_norm, subsets, epps, k_path, Nt, pseudos, det_per_ring, pRows, L, raw, size_z, osem_bool, 
				fileName, force_build, tube_width, crystal_size_z, x_center, y_center, z_center, size_center_x, size_center_y, size_of_x, size_center_z, 
				projector_type, header_directory, precompute_var, device, dec, n_rays, cr_pz, use_64bit_atomics, n_rekos, n_rekos_mlem, reko_type);
		}


		plhs[0] = cell_array_ptr;
	}
	// Compute the number of voxels each LOR traverses (using AF device)
	else if (type == 2) {

		if (nrhs != 34)
			mexErrMsgTxt("Invalid number of input arguments. There must be 34.");

		// Starting block
		const uint32_t block1 = (uint32_t)mxGetScalar(prhs[30]);

		// Ending block
		const uint32_t blocks = (uint32_t)mxGetScalar(prhs[31]);

		// Number of sinograms
		const uint32_t NSinos = (uint32_t)mxGetScalar(prhs[32]);

		// Total number of sinograms
		const uint16_t TotSinos = (uint16_t)mxGetScalar(prhs[33]);

		size_t loop_var_par = 1ULL;

		if (raw)
			loop_var_par = numRows / 2ULL;
		else
			loop_var_par = NSinos * size_x;

		plhs[0] = mxCreateNumericMatrix(loop_var_par, 1, mxUINT16_CLASS, mxREAL);

		uint16_t* lor = (uint16_t*)mxGetData(plhs[0]);

		find_LORs(lor, z_det, x, y, Nx, Ny, Nz, dx, dy, dz, bx, by, bz, bzb, maxxx, maxyy, zmax, NSlices, size_x, TotSinos, verbose, loop_var_par, 
			k_path, pseudos, det_per_ring, pRows, L, raw, size_z, fileName, device, size_of_x, force_build, header_directory);

	}

	// Clear ArrayFire memory
	af::deviceGC();

	return;
}