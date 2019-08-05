/**************************************************************************
* Matrix free computations for OMEGA.
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
	if (nrhs < 45)
		mexErrMsgTxt("Too few input arguments.  There must be at least 45.");
	else if (nrhs > 49)
		mexErrMsgTxt("Too many input arguments.  There can be at most 49.");

	if (nlhs != 1)
		mexErrMsgTxt("Invalid number of output arguments.  There must be exactly one.");

	// Load the input arguments
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

	const float bzb = bz + static_cast<float>(Nz) * dz;

	//std::string s(fileName);

	//s += (std::to_string(device) + ".bin");

	//// Filename for the OpenCL binaries
	//fileName = s.c_str();

	if (type < 2) {

		const float *atten = (float*)mxGetData(prhs[29]);
		const size_t size_atten = mxGetNumberOfElements(prhs[29]);

		const float* norm = (float*)mxGetData(prhs[30]);
		const size_t size_norm = mxGetNumberOfElements(prhs[30]);

		uint32_t *pituus = (uint32_t*)mxGetData(prhs[31]);

		const uint32_t attenuation_correction = (uint32_t)mxGetScalar(prhs[32]);

		const uint32_t normalization = (uint32_t)mxGetScalar(prhs[33]);

		const uint32_t Niter = (uint32_t)mxGetScalar(prhs[34]);

		const uint32_t subsets = (uint32_t)mxGetScalar(prhs[35]);

		const uint8_t* rekot = (uint8_t*)mxGetData(prhs[36]);
		const size_t size_reko = mxGetNumberOfElements(prhs[36]);

		const float epps = (float)mxGetScalar(prhs[37]);

		const uint16_t* lor1 = (uint16_t*)mxGetData(prhs[38]);
		const size_t koko_l = mxGetNumberOfElements(prhs[38]);

		const uint32_t* xy_index = (uint32_t*)mxGetData(prhs[39]);

		const uint16_t* z_index = (uint16_t*)mxGetData(prhs[40]);

		const bool osem_bool = (bool)mxGetScalar(prhs[41]);

		const float tube_width = (float)mxGetScalar(prhs[42]);

		const float crystal_size_z = (float)mxGetScalar(prhs[43]);

		const float *x_center = (float*)mxGetData(prhs[44]);
		const size_t size_center_x = mxGetNumberOfElements(prhs[44]);

		const float *y_center = (float*)mxGetData(prhs[45]);
		const size_t size_center_y = mxGetNumberOfElements(prhs[45]);

		const float* z_center = (float*)mxGetData(prhs[46]);
		const size_t size_center_z = mxGetNumberOfElements(prhs[46]);

		size_t koko;
		if (raw)
			koko = mxGetNumberOfElements(prhs[21]);
		else
			koko = mxGetNumberOfElements(prhs[38]);

		const size_t outSize = Nx * Ny * Nz;
		const size_t outSize2 = Niter + 1u;
		mxArray *cell_array_ptr;

		if (type == 0) {

			if (nrhs != 51)
				mexErrMsgTxt("Invalid number of input arguments. There must be 51.");

			const mxArray* Sin = prhs[48];

			const uint32_t Nt = (uint32_t)mxGetScalar(prhs[49]);

			const uint32_t projector_type = (uint32_t)mxGetScalar(prhs[50]);

			// Create the output cell array
			cell_array_ptr = mxCreateCellMatrix(size_reko + 1u, Nt);

			// Output dimensions
			const mwSize dim[4] = { static_cast<mwSize>(Nx), static_cast<mwSize>(Ny), static_cast<mwSize>(Nz), static_cast<mwSize>(Niter + 1u) };

			reconstruction_AF_matrixfree(koko, lor1, z_det, x, y, Sin, Nx, Ny, Nz, Niter, prhs[47], dx, dy, dz, bx, by, bz, bzb, maxxx, maxyy, zmax, NSlices, pituus, koko_l,
				xy_index, z_index, size_x, TotSinos, cell_array_ptr, dim, verbose, attenuation_correction, normalization, atten, size_atten, norm, size_norm, subsets, epps, rekot,
				k_path, size_reko, Nt, pseudos, det_per_ring, pRows, L, raw, size_z, osem_bool, fileName, force_build, tube_width, crystal_size_z, x_center, 
				y_center, z_center, size_center_x, size_center_y, size_of_x, size_center_z, projector_type, device);
		}
		else if (type == 1) {

			if (nrhs != 54)
				mexErrMsgTxt("Invalid number of input arguments. There must be 54.");

			const float * Sin = (float*)mxGetData(prhs[48]);

			const uint32_t osa_iter = (uint32_t)mxGetScalar(prhs[49]);

			const uint32_t tt = (uint32_t)mxGetScalar(prhs[50]);

			const uint32_t n_subsets = (uint32_t)mxGetScalar(prhs[51]);

			const bool mlem_bool = (bool)mxGetScalar(prhs[52]);

			const uint32_t projector_type = (uint32_t)mxGetScalar(prhs[53]);

			// Create the output cell array
			cell_array_ptr = mxCreateCellMatrix(size_reko + 3 + 7, 1);

			const mwSize dim = static_cast<mwSize>(outSize);

			//reconstruction_custom_matrixfree(koko, lor1, z_det, x, y, Sin, Nx, Ny, Nz, Niter, prhs[47], dx, dy, dz, bx, by, bz, bzb, maxxx, maxyy, zmax, NSlices, pituus, koko_l,
			//	xy_index, z_index, size_x, TotSinos, cell_array_ptr, dim, verbose, attenuation_correction, normalization, atten, size_atten, norm, size_norm, subsets, epps, rekot,
			//	k_path, size_reko, pseudos, det_per_ring, pRows, L, raw, size_z, osem_bool, fileName, force_build, osa_iter, tt, n_subsets, mlem_bool, tube_width, crystal_size_z, x_center, y_center, z_center, size_center_x,
			//	size_center_y, size_of_x, size_center_z, projector_type, device);
		}


		plhs[0] = cell_array_ptr;
	}
	else if (type == 2) {

		if (nrhs != 33)
			mexErrMsgTxt("Invalid number of input arguments. There must be 33.");

		const uint32_t block1 = (uint32_t)mxGetScalar(prhs[29]);

		const uint32_t blocks = (uint32_t)mxGetScalar(prhs[30]);

		const uint32_t NSinos = (uint32_t)mxGetScalar(prhs[31]);

		const uint16_t TotSinos = (uint16_t)mxGetScalar(prhs[32]);

		int loop_var_par = 1;

		if (raw)
			loop_var_par = numRows / 2;
		else
			loop_var_par = NSinos * size_x;

		plhs[0] = mxCreateNumericMatrix(loop_var_par, 1, mxUINT16_CLASS, mxREAL);

		uint16_t* lor = (uint16_t*)mxGetData(plhs[0]);

		find_LORs(lor, z_det, x, y, Nx, Ny, Nz, dx, dy, dz, bx, by, bz, bzb, maxxx, maxyy, zmax, NSlices, size_x, TotSinos, verbose, loop_var_par, 
			k_path, pseudos, det_per_ring, pRows, L, raw, size_z, fileName, device, size_of_x, force_build);

	}

	// Clear ArrayFire memory
	af::deviceGC();

	return;
}