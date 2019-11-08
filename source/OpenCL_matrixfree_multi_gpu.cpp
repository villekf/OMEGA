/**************************************************************************
* Matrix free computations for OMEGA for the multi-GPU/device case.
* Supports heterogeneous computing, but only on the same platform (i.e. if 
* you have a CPU and an integrated GPU from the same vendor and OpenCL
* runtime for both you can utilize both of them at the same time). Mixing
* CPU and GPU requires the optimization of the GPU to CPU value.
* This code can also be run for single device.
* This code is very similar to the other matrix-free code, but this one
* can also be run without installing ArrayFire, i.e. this uses pure OpenCL.
* For purely OSEM or MLEM reconstructions this code should be the fastest,
* regardless of the amount of devices used.
* Unlike the non-OpenCL versions, this one uses (32-bit) floats and thus
* can be slightly more inaccurate.
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
#include "functions_multigpu.hpp"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[]) {
	// Check for the number of input and output arguments
	if (nrhs < 32)
		mexErrMsgTxt("Too few input arguments.  There must be at least 32.");
	else if (nrhs > 60)
		mexErrMsgTxt("Too many input arguments.  There can be at most 60.");

	if (nlhs < 1)
		mexErrMsgTxt("Invalid number of output arguments.  There must be at least one.");
	else if (nlhs > 2)
		mexErrMsgTxt("Too many output arguments.  There can be at most two.");

	// Load the input arguments

	// Path to the kernel (.cl) files
	const char *k_path = mxArrayToString(prhs[0]);

	// Image size in y-direction
	const uint32_t Ny = (uint32_t)mxGetScalar(prhs[1]);

	// Image size in x-direction
	const uint32_t Nx = (uint32_t)mxGetScalar(prhs[2]);

	// Image size in z-direction
	const uint32_t Nz = (uint32_t)mxGetScalar(prhs[3]);

	// Distance between adjacent pixels in x-direction
	const float dx = (float)mxGetScalar(prhs[4]);

	// Distance between adjacent pixels in z-direction
	const float dz = (float)mxGetScalar(prhs[5]);

	// Distance of the pixel grid from the origin in y-direction
	const float by = (float)mxGetScalar(prhs[6]);

	// Distance of the pixel grid from the origin in x-direction
	const float bx = (float)mxGetScalar(prhs[7]);

	// Distance of the pixel grid from the origin in z-direction
	const float bz = (float)mxGetScalar(prhs[8]);
	
	// Coordinates of the detectors in z-direction
	const float *z_det = (float*)mxGetData(prhs[9]);
	const size_t size_z = mxGetNumberOfElements(prhs[9]);

	// Coordinates of the detectors in x-direction
	const float *x = (float*)mxGetData(prhs[10]);
	const size_t numel_x = mxGetNumberOfElements(prhs[10]);

	// Coordinates of the detectors in y-direction
	const float *y = (float*)mxGetData(prhs[11]);

	// Distance between adjacent pixels in y-direction
	const float dy = (float)mxGetScalar(prhs[12]);

	// The maximum elements of the pixel space in both x- and y-directions
	const float maxyy = (float)mxGetScalar(prhs[13]);

	const float maxxx = (float)mxGetScalar(prhs[14]);

	// Number of slices included
	const float NSlices = (float)mxGetScalar(prhs[15]);

	// Number of detector indices
	const uint32_t size_x = (uint32_t)mxGetScalar(prhs[16]);

	// Maximum value of the z-direction detector coordinates
	const float zmax = (float)mxGetScalar(prhs[17]);

	// Are status messages displayed
	const bool verbose = (bool)mxGetScalar(prhs[18]);

	// Detector pair numbers, for raw list-mode data
	const uint16_t *L = (uint16_t*)mxGetData(prhs[19]);
	const size_t numRows = mxGetM(prhs[19]);

	// Location (ring numbers) of pseudo rings, if present
	const uint32_t *pseudos = (uint32_t*)mxGetData(prhs[20]);
	const uint32_t pRows = (uint32_t)mxGetNumberOfElements(prhs[20]);

	// Number of detectors per ring
	const uint32_t det_per_ring = (uint32_t)mxGetScalar(prhs[21]);

	// Platform used
	const uint32_t device = (uint32_t)mxGetScalar(prhs[22]);

	// Filename of the current kernel (.cl) file
	const char *fileName = mxArrayToString(prhs[23]);

	// Is raw list-mode data used
	const uint8_t raw = (uint8_t)mxGetScalar(prhs[24]);

	// Coefficient that determines how much more measurements the other device (usually GPU) compared to the other device
	const float kerroin = (float)mxGetScalar(prhs[25]);

	// Which of the below sections are used (i.e. is this a precomputation phase, implementation 3 or forward-backwards projection)
	const uint32_t type = (uint32_t)mxGetScalar(prhs[26]);
	
	// Directory to look for OpenCL headers
	const char* header_directory = mxArrayToString(prhs[27]);

	const float bzb = bz + static_cast<float>(Nz) * dz;

	std::string s(fileName);

	s += (std::to_string(device) + ".bin");

	fileName = s.c_str();

	if (type == 2) {

		if (nrhs != 32)
			mexErrMsgTxt("Incorrect number of input arguments. There has to be 32.");

		// Starting ring
		const uint32_t block1 = (uint32_t)mxGetScalar(prhs[28]);

		// Ending ring
		const uint32_t blocks = (uint32_t)mxGetScalar(prhs[29]);

		// Number of sinograms used
		const uint32_t NSinos = (uint32_t)mxGetScalar(prhs[30]);

		// The total number of sinograms
		const uint16_t TotSinos = (uint16_t)mxGetScalar(prhs[31]);
		
		uint32_t loop_var_par = 1u;

		if (raw)
			loop_var_par = numRows / 2ULL;
		else
			loop_var_par = NSinos * size_x;

		plhs[0] = mxCreateNumericMatrix(loop_var_par, 1, mxUINT16_CLASS, mxREAL);

		uint16_t* lor = (uint16_t*)mxGetData(plhs[0]);

		// Find the number of voxels each LOR traverses
		find_LORs(lor, z_det, x, y, Nx, Ny, Nz, dx, dy, dz, bx, by, bz, bzb, maxxx, maxyy, zmax, NSlices, size_x, TotSinos, verbose, loop_var_par, 
			k_path, pseudos, det_per_ring, pRows, L, raw, size_z, fileName, device, kerroin, numel_x, header_directory);

	}
	else if (type < 2) {

		// attenuation values
		const float *atten = (float*)mxGetData(prhs[28]);
		const size_t size_atten = mxGetNumberOfElements(prhs[28]);

		// Normalization coefficients
		const float* norm = (float*)mxGetData(prhs[29]);
		const size_t size_norm = mxGetNumberOfElements(prhs[29]);

		// Number of measurements/LORs
		const uint32_t *pituus = (uint32_t*)mxGetData(prhs[30]);

		// Is the attenuation correction included
		const uint32_t attenuation_correction = (uint32_t)mxGetScalar(prhs[31]);

		// Is the attenuation correction included
		const uint32_t normalization = (uint32_t)mxGetScalar(prhs[32]);

		// Number of voxels the current LOR/ray traverses
		const uint16_t* lor1 = (uint16_t*)mxGetData(prhs[33]);
		const size_t koko_l = mxGetNumberOfElements(prhs[33]);

		// XY-indices of the detector coordinates of each LOR
		const uint32_t* xy_index = (uint32_t*)mxGetData(prhs[34]);

		// Z-indices of the detector coordinates of each LOR
		const uint16_t* z_index = (uint16_t*)mxGetData(prhs[35]);

		// For "D orthogonal, the "strip width"
		const float tube_width = (float)mxGetScalar(prhs[36]);

		// For 3D orthogonal, the "tube width"
		const float crystal_size_z = (float)mxGetScalar(prhs[37]);

		// Center coordinates of voxels in the X-dimension
		const float *x_center = (float*)mxGetData(prhs[38]);
		const size_t size_center_x = mxGetNumberOfElements(prhs[38]);

		// Center coordinates of voxels in the Y-dimension
		const float *y_center = (float*)mxGetData(prhs[39]);
		const size_t size_center_y = mxGetNumberOfElements(prhs[39]);

		// Center coordinates of voxels in the Z-dimension
		const float* z_center = (float*)mxGetData(prhs[40]);
		const size_t size_center_z = mxGetNumberOfElements(prhs[40]);

		// Randoms
		const mxArray* sc_ra = prhs[41];

		// Randoms corrections
		const uint32_t randoms_correction = (uint32_t)mxGetScalar(prhs[42]);

		// The type of projector used (Siddon or orthogonal)
		const uint32_t projector_type = (uint32_t)mxGetScalar(prhs[43]);

		// If true, then the precomputed LOR voxels counts are used
		const bool precompute_var = (bool)mxGetScalar(prhs[44]);

		// Accuracy factor in orthogonal distance
		const int32_t dec = (int32_t)mxGetScalar(prhs[45]);

		// Number of rays in Siddon
		uint16_t n_rays = (uint16_t)mxGetScalar(prhs[46]);

		// Crystal pitch in z-direction
		const double cr_pz = (double)mxGetScalar(prhs[47]);

		size_t koko;
		if (raw)
			koko = mxGetNumberOfElements(prhs[19]);
		else
			koko = mxGetNumberOfElements(prhs[33]);

		if (type == 0) {

			if (nrhs != 50)
				mexErrMsgTxt("Incorrect number of input arguments. There has to be 50.");

			// Right hand side for forward or backprojections
			const float* rhs = (float*)mxGetData(prhs[48]);
			const size_t size_rhs = mxGetNumberOfElements(prhs[48]);

			// Is the normalization constant computed (sum(A))
			const bool no_norm = (bool)mxGetScalar(prhs[49]);

			size_t outSize;
			if (size_rhs == Nx * Ny*Nz)
				outSize = pituus[0];
			else
				outSize = Nx * Ny*Nz;
			const size_t outSize2 = 1;
			const size_t outSize3 = Nx * Ny*Nz;

			plhs[0] = mxCreateNumericMatrix(outSize, outSize2, mxSINGLE_CLASS, mxREAL);
			plhs[1] = mxCreateNumericMatrix(outSize3, outSize2, mxSINGLE_CLASS, mxREAL);

			float* output = (float*)mxGetData(plhs[0]);
			float* normalizer = (float*)mxGetData(plhs[1]);

			const uint16_t TotSinos = size_z / 2;

			// Forward/backward projection
			reconstruction_f_b_proj(koko, lor1, z_det, x, y, rhs, sc_ra, Nx, Ny, Nz, dx, dy, dz, bx, by, bz, bzb, maxxx, maxyy, zmax, NSlices, pituus, 
				koko_l,	xy_index, z_index, size_x, TotSinos, verbose, randoms_correction, attenuation_correction, normalization, atten, size_atten, norm, 
				size_norm, k_path, pseudos, det_per_ring, pRows, L, raw, size_z, fileName, device, kerroin, output, normalizer, size_rhs, no_norm, numel_x, 
				tube_width, crystal_size_z, x_center, y_center, z_center, size_center_x, size_center_y, size_center_z, projector_type, header_directory, 
				precompute_var, dec, n_rays, cr_pz);

		}
		else if (type == 1) {

			if (nrhs != 60)
				mexErrMsgTxt("Incorrect number of input arguments. There has to be 60.");

			// Number of sinograms used
			const uint32_t NSinos = (uint32_t)mxGetScalar(prhs[48]);

			// Total number of sinograms
			const uint16_t TotSinos = (uint16_t)mxGetScalar(prhs[49]);

			// Number of iterations
			const uint32_t Niter = (uint32_t)mxGetScalar(prhs[50]);

			// Number of subsets
			const uint32_t subsets = (uint32_t)mxGetScalar(prhs[51]);

			// Which reconstruction methods are used
			const uint8_t* rekot = (uint8_t*)mxGetData(prhs[52]);
			const size_t size_reko = mxGetNumberOfElements(prhs[52]);

			// Epsilon value
			const float epps = (float)mxGetScalar(prhs[53]);

			// Measurement data
			const mxArray* Sin = prhs[54];

			// Number of time steps
			const uint32_t Nt = (uint32_t)mxGetScalar(prhs[55]);

			// Is OSEM used
			const bool osem_bool = (bool)mxGetScalar(prhs[56]);

			// Force re-building of OpenCL kernels (unused)
			const bool force_build = (bool)mxGetScalar(prhs[57]);

			// Use 64-bit integer atomic functions if possible
			const bool use_64bit_atomics = (bool)mxGetScalar(prhs[59]);

			const size_t outSize = Nx * Ny * Nz;
			const size_t outSize2 = Niter + 1;

			// Create the output cell array
			mxArray *cell_array_ptr;
			cell_array_ptr = mxCreateCellMatrix(1, Nt);

			// Implementation 3 (multi-device OpenCL)
			reconstruction_multigpu(koko, lor1, z_det, x, y, Sin, sc_ra, Nx, Ny, Nz, Niter, prhs[58], dx, dy, dz, bx, by, bz, bzb, maxxx, maxyy, zmax, 
				NSlices, pituus, koko_l, xy_index, z_index, size_x, TotSinos, cell_array_ptr, verbose, randoms_correction, attenuation_correction, 
				normalization, atten, size_atten, norm, size_norm, subsets, epps, rekot, k_path, size_reko, Nt, pseudos, det_per_ring, pRows, L, raw, 
				size_z, osem_bool, fileName, force_build, device, kerroin, numel_x, tube_width, crystal_size_z, x_center, y_center, z_center, size_center_x, 
				size_center_y, size_center_z, projector_type, header_directory, precompute_var, dec, n_rays, cr_pz, use_64bit_atomics);

			plhs[0] = cell_array_ptr;
		}
	}

	return;
}