/**************************************************************************
* Loads the Siemens Inveon 48-bit list-mode data and stores the detector
* indices of prompts and delays at specific time steps.
*
* Copyright (C) 2020-2024 Ville-Veikko Wettenhovi
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
#define MATLAB
#include "inveon.h"

void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])

{


	/* Check for proper number of arguments */

	if (nrhs != 17) {
		mexErrMsgIdAndTxt("MATLAB:list2matlab_aivi:invalidNumInputs",
			"17 input arguments required.");
	}
	else if (nlhs > 7) {
		mexErrMsgIdAndTxt("MATLAB:list2matlab_aivi:maxlhs",
			"Too many output arguments.");
	}

	const double* vali = (double*)mxGetData(prhs[1]);
	const double alku = (double)mxGetScalar(prhs[2]);
	const double loppu = (double)mxGetScalar(prhs[3]);
	const size_t pituus = (size_t)mxGetScalar(prhs[4]);
	const uint32_t detectors = (uint32_t)mxGetScalar(prhs[5]);
	const bool randoms_correction = (bool)mxGetScalar(prhs[6]);
	const uint64_t sinoSize = (uint64_t)mxGetScalar(prhs[7]);
	const bool saveRawData = (bool)mxGetScalar(prhs[8]);
	const uint32_t Ndist = (uint32_t)mxGetScalar(prhs[9]);
	const uint32_t Nang = (uint32_t)mxGetScalar(prhs[10]);
	const uint32_t ringDifference = (uint32_t)mxGetScalar(prhs[11]);
	const uint32_t span = (uint32_t)mxGetScalar(prhs[12]);
	const uint32_t* seg = (uint32_t*)mxGetData(prhs[13]);
	const uint64_t NT = (uint64_t)mxGetScalar(prhs[14]);
	const int32_t nDistSide = (int32_t)mxGetScalar(prhs[15]);
	const bool storeCoordinates = (bool)mxGetScalar(prhs[16]);

	if (storeCoordinates) {
		plhs[0] = mxCreateNumericMatrix(pituus, 1, mxUINT16_CLASS, mxREAL);
		plhs[1] = mxCreateNumericMatrix(pituus, 1, mxUINT16_CLASS, mxREAL);
		if (randoms_correction) {
			plhs[3] = mxCreateNumericMatrix(pituus, 1, mxUINT16_CLASS, mxREAL);
			plhs[4] = mxCreateNumericMatrix(pituus, 1, mxUINT16_CLASS, mxREAL);
		}
		else {
			plhs[3] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
			plhs[4] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		}
	}
	else {
		plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		plhs[1] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		plhs[3] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		plhs[4] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
	}
	// tpoints is only used when the coordinates are stored
	if (NT > 1 && storeCoordinates)
		plhs[2] = mxCreateNumericMatrix(pituus, 1, mxUINT16_CLASS, mxREAL);
	else
		plhs[2] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
	if (!storeCoordinates)
		plhs[5] = mxCreateNumericMatrix(sinoSize * NT, 1, mxUINT16_CLASS, mxREAL);
	else
		plhs[5] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
	if (randoms_correction && !storeCoordinates)
		plhs[6] = mxCreateNumericMatrix(sinoSize * NT, 1, mxUINT16_CLASS, mxREAL);
	else
		plhs[6] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);


	/* Assign pointers to the various parameters */
	uint16_t* LL1 = (uint16_t*)mxGetData(plhs[0]);
	uint16_t* LL2 = (uint16_t*)mxGetData(plhs[1]);
	uint16_t* tpoints = (uint16_t*)mxGetData(plhs[2]);
	uint16_t* DD1 = (uint16_t*)mxGetData(plhs[3]);
	uint16_t* DD2 = (uint16_t*)mxGetData(plhs[4]);
	uint16_t* Sino = (uint16_t*)mxGetData(plhs[5]);
	uint16_t* SinoD = (uint16_t*)mxGetData(plhs[6]);

	// Check for char type
	if (!mxIsChar(prhs[0]))
		mexErrMsgTxt("First input argument is not char");

	/* Pointer to character array */
	const char* argv = mxArrayToString(prhs[0]);

	uint64_t nPrompts = 0, nDelays = 0;
	histogram(LL1, LL2, tpoints, argv, vali, alku, loppu, detectors, pituus, randoms_correction, DD1, DD2, Sino, SinoD, saveRawData,
		Ndist, Nang, ringDifference, span, sinoSize, seg, nDistSide, storeCoordinates, NT, &nPrompts, &nDelays);

	// Shrink the coordinate arrays to the actual number of stored events
	if (storeCoordinates) {
		const size_t nP = std::max(static_cast<size_t>(nPrompts), static_cast<size_t>(1));
		mxSetData(plhs[0], mxRealloc(LL1, nP * sizeof(uint16_t)));
		mxSetM(plhs[0], static_cast<size_t>(nPrompts));
		mxSetData(plhs[1], mxRealloc(LL2, nP * sizeof(uint16_t)));
		mxSetM(plhs[1], static_cast<size_t>(nPrompts));
		if (NT > 1) {
			mxSetData(plhs[2], mxRealloc(tpoints, nP * sizeof(uint16_t)));
			mxSetM(plhs[2], static_cast<size_t>(nPrompts));
		}
		if (randoms_correction) {
			const size_t nD = std::max(static_cast<size_t>(nDelays), static_cast<size_t>(1));
			mxSetData(plhs[3], mxRealloc(DD1, nD * sizeof(uint16_t)));
			mxSetM(plhs[3], static_cast<size_t>(nDelays));
			mxSetData(plhs[4], mxRealloc(DD2, nD * sizeof(uint16_t)));
			mxSetM(plhs[4], static_cast<size_t>(nDelays));
		}
	}
	return;

}
