/**************************************************************************
* Computes the system matrix for OMEGA.
* This version also checks whether the LOR/ray intersects the pixel space,
* i.e. it doesn't require any precomputation.
* This version outputs the row and column indices and values that can be
* used to create a sparse matrix. Orthogonal, improved Siddon and original 
* are Siddon are available projectors.
* This is slower than the precomputed versions.
* 
* Copyright (C) 2019 Ville-Veikko Wettenhovi
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

void mexFunction(int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[])

{
	// Check for the number of input and output arguments
	if (nrhs < 32)
		mexErrMsgTxt("Too few input arguments.  There must be at least 32.");
	if (nrhs > 35)
		mexErrMsgTxt("Too many input arguments.  There can be at most 35.");

	if (nlhs != 3)
		mexErrMsgTxt("Invalid number of output arguments.  There must be exactly three.");

	// Load the input arguments
	const bool verbose = (bool)mxGetScalar(prhs[0]);

	const uint32_t Ny = (uint32_t)mxGetScalar(prhs[1]);

	const uint32_t Nx = (uint32_t)mxGetScalar(prhs[2]);

	const uint32_t Nz = (uint32_t)mxGetScalar(prhs[3]);

	const double d = (double)mxGetScalar(prhs[4]);

	const double dz = (double)mxGetScalar(prhs[5]);

	const double by = (double)mxGetScalar(prhs[6]);

	const double bx = (double)mxGetScalar(prhs[7]);

	const double bz = (double)mxGetScalar(prhs[8]);

	const double *z_det = (double*)mxGetData(prhs[9]);

	const double *x = (double*)mxGetData(prhs[10]);

	const double *y = (double*)mxGetData(prhs[11]);

	const double dy = (double)mxGetScalar(prhs[12]);

	const double *yy = (double*)mxGetData(prhs[13]);

	const double *xx = (double*)mxGetData(prhs[14]);

	const uint32_t NSinos = (uint32_t)mxGetScalar(prhs[15]);

	const uint32_t NSlices = (uint32_t)mxGetScalar(prhs[16]);

	// from array to std::vector

	vector<double> yy_vec(yy, yy + mxGetNumberOfElements(prhs[13]) - 1);

	vector<double> xx_vec(xx, xx + mxGetNumberOfElements(prhs[14]) - 1);

	const uint32_t size_x = (uint32_t)mxGetScalar(prhs[17]);

	const double zmax = (double)mxGetScalar(prhs[18]);

	const uint32_t TotSinos = (uint32_t)mxGetScalar(prhs[19]);

	const uint32_t ind_size = (uint32_t)mxGetScalar(prhs[20]);

	const double *atten = (double*)mxGetData(prhs[21]);

	const uint32_t *index = (uint32_t*)mxGetData(prhs[22]);
	const size_t index_size = mxGetNumberOfElements(prhs[22]);

	const uint32_t pituus = (uint32_t)mxGetScalar(prhs[23]);

	const bool attenuation_correction = (bool)mxGetScalar(prhs[24]);

	const bool raw = (bool)mxGetScalar(prhs[25]);

	const uint16_t *L = (uint16_t*)mxGetData(prhs[26]);
	const size_t numRows = mxGetM(prhs[26]);

	const uint32_t *pseudos = (uint32_t*)mxGetData(prhs[27]);
	const size_t pRows = mxGetM(prhs[27]);

	const uint32_t block1 = (uint32_t)mxGetScalar(prhs[28]);

	const uint32_t blocks = (uint32_t)mxGetScalar(prhs[29]);

	const uint32_t det_per_ring = (uint32_t)mxGetScalar(prhs[30]);

	const uint32_t proj_type = (uint32_t)mxGetScalar(prhs[31]);



	uint32_t loop_var_par = 1u;

	//vector<uint16_t> lor;

	if (index_size > 1u && !raw) {
		//loop_var_par = index[pituus - 1];
		loop_var_par = pituus;
		//lor.assign(pituus, static_cast<uint16_t>(0));
	}
	else if (!raw) {
		loop_var_par = NSinos * size_x;
		//lor.assign(loop_var_par, static_cast<uint16_t>(0));
	}
	else {
		loop_var_par = numRows / 2u;
	}

	plhs[0] = mxCreateNumericMatrix(loop_var_par, 1, mxUINT16_CLASS, mxREAL);

	uint16_t* lor = (uint16_t *)mxGetData(plhs[0]);

	vector<uint32_t> indices;

	vector<double> elements;

	indices.reserve(ind_size);
	elements.reserve(ind_size);

	// The maximum elements of the pixel space in both x- and y-directions
	const double maxyy = yy_vec.back();
	const double maxxx = xx_vec.back();

	clock_t time = clock();

	uint32_t lj;

	if (proj_type == 1u) {

		// run the Improved Siddon's algorithm
		lj = improved_siddon_no_precompute(loop_var_par, size_x, zmax, TotSinos, indices, elements, lor, maxyy, maxxx, xx_vec, dy,
			yy_vec, atten, x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, index_size, index, attenuation_correction, raw, det_per_ring, blocks, block1,
			L, pseudos, pRows);
	}
	else if (proj_type == 2u) {

		const double crystal_size = (double)mxGetScalar(prhs[32]);

		const double *y_center = (double*)mxGetData(prhs[33]);

		const double *x_center = (double*)mxGetData(prhs[34]);

		// run the Orthogonal Siddon algorithm
		lj = orth_siddon_no_precompute(loop_var_par, size_x, zmax, TotSinos, indices, elements, lor, maxyy, maxxx, xx_vec, dy,
			yy_vec, atten, x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, index_size, index, attenuation_correction, raw, det_per_ring, blocks, block1,
			L, pseudos, pRows, crystal_size, y_center, x_center);
	}
	else if(proj_type == 0u) {

		const double *iij = (double*)mxGetData(prhs[32]);

		const double *jji = (double*)mxGetData(prhs[33]);

		const double *kkj = (double*)mxGetData(prhs[34]);

		const vector<double> iij_vec(iij, iij + mxGetNumberOfElements(prhs[32]));

		const vector<double> jjk_vec(jji, jji + mxGetNumberOfElements(prhs[33]));

		const vector<double> kkj_vec(kkj, kkj + mxGetNumberOfElements(prhs[34]));

		// run the original Siddon's algorithm
		lj = original_siddon_no_precompute(loop_var_par, size_x, zmax, TotSinos, indices, elements, lor, maxyy, maxxx, xx_vec, dy,
			yy_vec, atten, x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, index_size, index, attenuation_correction, raw, det_per_ring, blocks, block1,
			L, pseudos, pRows, iij_vec, jjk_vec, kkj_vec);
	}
    
	time = clock() - time;

	if (verbose) {
		mexPrintf("Function elapsed time is %f seconds\n", ((float)time) / CLOCKS_PER_SEC);
		mexEvalString("pause(.001);");
	}

	const size_t outSize1 = lj * 2u;
	const size_t outSize2 = 1;

	// Create the MATLAB output vectors (row and column indices, elements)

	plhs[1] = mxCreateNumericMatrix(indices.size(), outSize2, mxUINT32_CLASS, mxREAL);

	uint32_t* outputMatrix2 = (uint32_t *)mxGetData(plhs[1]);

	time = clock();

	copy(indices.begin(), indices.end(), outputMatrix2);
	indices.erase(indices.begin(), indices.end());
	indices.shrink_to_fit();

	time = clock() - time;
	if (verbose)
		mexPrintf("indices copy and erase elapsed time is %f seconds\n", ((float)time) / CLOCKS_PER_SEC);
	time = clock();

	plhs[2] = mxCreateNumericMatrix(elements.size(), outSize2, mxDOUBLE_CLASS, mxREAL);

	double* outputMatrix3 = (double *)mxGetData(plhs[2]);

	copy(elements.begin(), elements.end(), outputMatrix3);

	time = clock() - time;
	if (verbose)
		mexPrintf("elements copy elapsed time is %f seconds\n", ((float)time) / CLOCKS_PER_SEC);
}