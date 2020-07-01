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
* Copyright(C) 2020 Ville-Veikko Wettenhovi
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
	else if (nrhs > 67)
		mexErrMsgTxt("Too many input arguments.  There can be at most 67.");

	if (nlhs < 1)
		mexErrMsgTxt("Invalid number of output arguments.  There must be at least one.");
	else if (nlhs > 1)
		mexErrMsgTxt("Too many output arguments.  There can be at most one.");

	// Load the input arguments
	int ind = 0;

	// Path to the kernel (.cl) files
	const char *k_path = mxArrayToString(prhs[ind]);
	ind++;

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
	const float dx = (float)mxGetScalar(prhs[ind]);
	ind++;

	// Distance between adjacent pixels in z-direction
	const float dz = (float)mxGetScalar(prhs[ind]);
	ind++;

	// Distance of the pixel grid from the origin in y-direction
	const float by = (float)mxGetScalar(prhs[ind]);
	ind++;

	// Distance of the pixel grid from the origin in x-direction
	const float bx = (float)mxGetScalar(prhs[ind]);
	ind++;

	// Distance of the pixel grid from the origin in z-direction
	const float bz = (float)mxGetScalar(prhs[ind]);
	ind++;
	
	// Coordinates of the detectors in z-direction
	const float *z_det = (float*)mxGetData(prhs[ind]);
	const size_t size_z = mxGetNumberOfElements(prhs[ind]);
	ind++;

	// Coordinates of the detectors in x-direction
	const float *x = (float*)mxGetData(prhs[ind]);
	const size_t numel_x = mxGetNumberOfElements(prhs[ind]);
	ind++;

	// Coordinates of the detectors in y-direction
	const float *y = (float*)mxGetData(prhs[ind]);
	ind++;

	// Distance between adjacent pixels in y-direction
	const float dy = (float)mxGetScalar(prhs[ind]);
	ind++;

	// The maximum elements of the pixel space in both x- and y-directions
	const float maxyy = (float)mxGetScalar(prhs[ind]);
	ind++;

	const float maxxx = (float)mxGetScalar(prhs[ind]);
	ind++;

	// Number of slices included
	const float NSlices = (float)mxGetScalar(prhs[ind]);
	ind++;

	// Number of detector indices
	const uint32_t size_x = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	// Maximum value of the z-direction detector coordinates
	const float zmax = (float)mxGetScalar(prhs[ind]);
	ind++;

	// Are status messages displayed
	const bool verbose = (bool)mxGetScalar(prhs[ind]);
	ind++;

	// Detector pair numbers, for raw list-mode data
	const uint16_t *L = (uint16_t*)mxGetData(prhs[ind]);
	const size_t numRows = mxGetM(prhs[ind]);
	ind++;

	// Location (ring numbers) of pseudo rings, if present
	const uint32_t *pseudos = (uint32_t*)mxGetData(prhs[ind]);
	const uint32_t pRows = (uint32_t)mxGetNumberOfElements(prhs[ind]);
	ind++;

	// Number of detectors per ring
	const uint32_t det_per_ring = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	// Platform used
	const uint32_t device = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;

	// Filename of the current kernel (.cl) file
	const char *fileName = mxArrayToString(prhs[ind]);
	ind++;

	// Is raw list-mode data used
	const uint8_t raw = (uint8_t)mxGetScalar(prhs[ind]);
	ind++;

	// Coefficient that determines how much more measurements the other device (usually GPU) compared to the other device
	const float kerroin = (float)mxGetScalar(prhs[ind]);
	ind++;

	// Which of the below sections are used (i.e. is this a precomputation phase, implementation 3 or forward-backwards projection)
	const uint32_t type = (uint32_t)mxGetScalar(prhs[ind]);
	ind++;
	
	// Directory to look for OpenCL headers
	const char* header_directory = mxArrayToString(prhs[ind]);
	ind++;

	const float bzb = bz + static_cast<float>(Nz) * dz;

	std::string s(fileName);

	s += (std::to_string(device) + ".bin");

	fileName = s.c_str();

	// Fixed local size
	const size_t local_size = 64ULL;

	if (type == 2) {

		if (nrhs != 32)
			mexErrMsgTxt("Incorrect number of input arguments. There has to be 32.");

		// Starting ring
		const uint32_t block1 = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		// Ending ring
		const uint32_t blocks = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		// Number of sinograms used
		const uint32_t NSinos = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		// The total number of sinograms
		const uint32_t TotSinos = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;
		
		uint32_t loop_var_par = 1u;

		if (raw)
			loop_var_par = numRows / 2ULL;
		else
			loop_var_par = NSinos * size_x;

		plhs[0] = mxCreateNumericMatrix(loop_var_par, 1, mxUINT16_CLASS, mxREAL);

		uint16_t* lor = (uint16_t*)mxGetData(plhs[0]);

		// Find the number of voxels each LOR traverses
		find_LORs(lor, z_det, x, y, Nx, Ny, Nz, dx, dy, dz, bx, by, bz, bzb, maxxx, maxyy, zmax, NSlices, size_x, TotSinos, verbose, loop_var_par, 
			k_path, pseudos, det_per_ring, pRows, L, raw, size_z, fileName, device, numel_x, header_directory, local_size);

	}
	else if (type < 2) {

		// attenuation values
		const float *atten = (float*)mxGetData(prhs[ind]);
		const size_t size_atten = mxGetNumberOfElements(prhs[ind]);
		ind++;

		// Normalization coefficients
		const float* norm = (float*)mxGetData(prhs[ind]);
		const size_t size_norm = mxGetNumberOfElements(prhs[ind]);
		ind++;

		// Number of measurements/LORs
		const uint32_t *pituus = (uint32_t*)mxGetData(prhs[ind]);
		ind++;

		// Is the attenuation correction included
		const uint32_t attenuation_correction = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		// Is the attenuation correction included
		const uint32_t normalization = (uint32_t)mxGetScalar(prhs[ind]);
		ind++;

		// Number of voxels the current LOR/ray traverses
		const uint16_t* lor1 = (uint16_t*)mxGetData(prhs[ind]);
		const size_t koko_l = mxGetNumberOfElements(prhs[ind]);
		ind++;

		// XY-indices of the detector coordinates of each LOR
		const uint32_t* xy_index = (uint32_t*)mxGetData(prhs[ind]);
		ind++;

		// Z-indices of the detector coordinates of each LOR
		const uint16_t* z_index = (uint16_t*)mxGetData(prhs[ind]);
		ind++;

		// For "D orthogonal, the "strip width"
		const float tube_width = (float)mxGetScalar(prhs[ind]);
		ind++;

		// For 3D orthogonal, the "tube width"
		const float crystal_size_z = (float)mxGetScalar(prhs[ind]);
		ind++;

		// Center coordinates of voxels in the X-dimension
		const float *x_center = (float*)mxGetData(prhs[ind]);
		const size_t size_center_x = mxGetNumberOfElements(prhs[ind]);
		ind++;

		// Center coordinates of voxels in the Y-dimension
		const float *y_center = (float*)mxGetData(prhs[ind]);
		const size_t size_center_y = mxGetNumberOfElements(prhs[ind]);
		ind++;

		// Center coordinates of voxels in the Z-dimension
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

		// Accuracy factor in orthogonal distance
		const int32_t dec = (int32_t)mxGetScalar(prhs[ind]);
		ind++;

		// Number of rays in Siddon (transaxial)
		uint16_t n_rays = (uint16_t)mxGetScalar(prhs[ind]);
		ind++;

		// Number of rays in Siddon (axial)
		uint16_t n_rays3D = (uint16_t)mxGetScalar(prhs[ind]);
		ind++;

		// Crystal pitch in z-direction
		const float cr_pz = (float)mxGetScalar(prhs[ind]);
		ind++;

		// Measurement data
		const mxArray* Sin = prhs[ind];
		ind++;

		// Use 64-bit integer atomic functions if possible
		const bool use_64bit_atomics = (bool)mxGetScalar(prhs[ind]);
		ind++;

		size_t koko;
		if (raw)
			koko = numRows;
		else
			koko = koko_l;

		if (type == 0) {

			if (nrhs != 60)
				mexErrMsgTxt("Incorrect number of input arguments. There has to be 60.");

			// Right hand side for forward or backprojections
			const float* rhs = (float*)mxGetData(prhs[ind]);
			const size_t size_rhs = mxGetNumberOfElements(prhs[ind]);
			ind++;

			// Is the normalization constant computed (sum(A))
			const bool no_norm = (bool)mxGetScalar(prhs[ind]);
			ind++;

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

			// Use PSF in Siddon
			const bool use_psf = (bool)mxGetScalar(prhs[ind]);
			ind++;

			size_t outSize;
			if (size_rhs == Nx * Ny * Nz)
				outSize = pituus[0];
			else
				outSize = Nx * Ny * Nz;
			const size_t outSize2 = 1;
			const size_t outSize3 = Nx * Ny * Nz;

			const uint16_t TotSinos = size_z / 2ULL;

			mxArray* output;
			output = mxCreateCellMatrix(2, 1);

			// Forward/backward projection
			reconstruction_f_b_proj(koko, lor1, z_det, x, y, rhs, sc_ra, Nx, Ny, Nz, dx, dy, dz, bx, by, bz, bzb, maxxx, maxyy, zmax, NSlices, pituus,
				koko_l, xy_index, z_index, size_x, TotSinos, verbose, randoms_correction, attenuation_correction, normalization, atten, size_atten, norm,
				size_norm, k_path, pseudos, det_per_ring, pRows, L, raw, size_z, fileName, device, kerroin, output, size_rhs, no_norm, numel_x,
				tube_width, crystal_size_z, x_center, y_center, z_center, size_center_x, size_center_y, size_center_z, projector_type, header_directory,
				precompute_var, dec, n_rays, n_rays3D, cr_pz, Sin, use_64bit_atomics, global_factor, bmin, bmax, Vmax, V, size_V, local_size, false, prhs[ind]);
			plhs[0] = output;

		}
		else if (type == 1) {

			if (nrhs != 67)
				mexErrMsgTxt("Incorrect number of input arguments. There has to be 67.");

			// Number of sinograms used
			const uint32_t NSinos = (uint32_t)mxGetScalar(prhs[ind]);
			ind++;

			// Total number of sinograms
			const uint16_t TotSinos = (uint16_t)mxGetScalar(prhs[ind]);
			ind++;

			// Number of iterations
			const uint32_t Niter = (uint32_t)mxGetScalar(prhs[ind]);
			ind++;

			// Number of subsets
			const uint32_t subsets = (uint32_t)mxGetScalar(prhs[ind]);
			ind++;

			// Which reconstruction methods are used
			const uint8_t* rekot = (uint8_t*)mxGetData(prhs[ind]);
			const size_t size_reko = mxGetNumberOfElements(prhs[ind]);
			ind++;

			// Epsilon value
			const float epps = (float)mxGetScalar(prhs[ind]);
			ind++;

			// Number of time steps
			const uint32_t Nt = (uint32_t)mxGetScalar(prhs[ind]);
			ind++;

			// Is OSEM used
			const bool osem_bool = (bool)mxGetScalar(prhs[ind]);
			ind++;

			// Use PSF in Siddon
			const bool use_psf = (bool)mxGetScalar(prhs[ind]);
			ind++;

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

			const size_t outSize = Nx * Ny * Nz;
			const size_t outSize2 = Niter + 1u;

			// Create the output cell array
			mxArray *cell_array_ptr;
			cell_array_ptr = mxCreateCellMatrix(1, Nt);

			// Implementation 3 (multi-device OpenCL)
			reconstruction_multigpu(koko, lor1, z_det, x, y, Sin, sc_ra, Nx, Ny, Nz, Niter, prhs[ind], dx, dy, dz, bx, by, bz, bzb, maxxx, maxyy, zmax, 
				NSlices, pituus, koko_l, xy_index, z_index, size_x, TotSinos, cell_array_ptr, verbose, randoms_correction, attenuation_correction, 
				normalization, atten, size_atten, norm, size_norm, subsets, epps, rekot, k_path, size_reko, Nt, pseudos, det_per_ring, pRows, L, raw, 
				size_z, osem_bool, fileName, use_psf, device, kerroin, numel_x, tube_width, crystal_size_z, x_center, y_center, z_center, size_center_x, 
				size_center_y, size_center_z, projector_type, header_directory, precompute_var, dec, n_rays, n_rays3D, cr_pz, use_64bit_atomics, global_factor,
				bmin, bmax, Vmax, V, size_V, local_size, gaussian, size_gauss);

			plhs[0] = cell_array_ptr;
		}
	}

	return;
}