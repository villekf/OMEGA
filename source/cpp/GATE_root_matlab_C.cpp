/**************************************************************************
* ROOT file import into MATLAB. This file contains the C implementation.
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
#include "rootImport.h"
#include "TChain.h"
extern "C" mxArray * mxCreateSharedDataCopy(const mxArray * pr);



void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])

{


	/* Check for proper number of arguments */

	if (nrhs != 50) {
		mexErrMsgIdAndTxt("MATLAB:GATE_root_matlab:invalidNumInputs",
			"50 input arguments required.");
	}
	else if (nlhs > 15) {
		mexErrMsgIdAndTxt("MATLAB:GATE_root_matlab:maxlhs",
			"Too many output arguments.");
	}

	const double* tPoints = (double*)mxGetData(prhs[1]);
	const double alku = (double)mxGetScalar(prhs[2]);
	const double loppu = (double)mxGetScalar(prhs[3]);
	const uint32_t* detectors = (uint32_t*)mxGetData(prhs[4]);
	const uint32_t blocks_per_ring = (uint32_t)mxGetScalar(prhs[5]);
	const uint32_t* cryst_per_block = (uint32_t*)mxGetData(prhs[6]);
	const uint32_t* det_per_ring = (uint32_t*)mxGetData(prhs[7]);
	const uint32_t linear_multp = (uint32_t)mxGetScalar(prhs[8]);
	bool source = (bool)mxGetScalar(prhs[9]);
	const int64_t Nt = (int64_t)mxGetScalar(prhs[10]);
	bool obtain_trues = (bool)mxGetScalar(prhs[11]);
	bool store_scatter = (bool)mxGetScalar(prhs[12]);
	bool store_randoms = (bool)mxGetScalar(prhs[13]);
	uint8_t* scatter_components = (uint8_t*)mxGetData(prhs[14]);
	bool randoms_correction = (bool)mxGetScalar(prhs[15]);
	bool store_coordinates = (bool)mxGetScalar(prhs[16]);
	const uint32_t transaxial_multip = (uint32_t)mxGetScalar(prhs[17]);
	const uint32_t* cryst_per_block_z = (uint32_t*)mxGetData(prhs[18]);
	const uint32_t* rings = (uint32_t*)mxGetData(prhs[19]);
	const uint64_t* sinoSize = (uint64_t*)mxGetData(prhs[20]);
	const bool TOF = (bool)mxGetScalar(prhs[21]);
	const bool verbose = (bool)mxGetScalar(prhs[22]);
	const int32_t nLayers = (int32_t)mxGetScalar(prhs[23]);
	const uint32_t Ndist = (uint32_t)mxGetScalar(prhs[24]);
	const uint32_t* Nang = (uint32_t*)mxGetData(prhs[25]);
	const uint32_t ringDifference = (uint32_t)mxGetScalar(prhs[26]);
	const uint32_t span = (uint32_t)mxGetScalar(prhs[27]);
	const uint32_t* seg = (uint32_t*)mxGetData(prhs[28]);
	const uint64_t TOFSize = (uint64_t)mxGetScalar(prhs[29]);
	const int32_t nDistSide = (int32_t)mxGetScalar(prhs[30]);
	const uint32_t* detWPseudo = (uint32_t*)mxGetData(prhs[36]);
	const int32_t nPseudos = (int32_t)mxGetScalar(prhs[37]);
	const double binSize = (double)mxGetScalar(prhs[38]);
	const double FWHM = (double)mxGetScalar(prhs[39]);
	const float dx = (float)mxGetScalar(prhs[40]);
	const float dy = (float)mxGetScalar(prhs[41]);
	const float dz = (float)mxGetScalar(prhs[42]);
	const float bx = (float)mxGetScalar(prhs[43]);
	const float by = (float)mxGetScalar(prhs[44]);
	const float bz = (float)mxGetScalar(prhs[45]);
	const int64_t Nx = (int64_t)mxGetScalar(prhs[46]);
	const int64_t Ny = (int64_t)mxGetScalar(prhs[47]);
	const int64_t Nz = (int64_t)mxGetScalar(prhs[48]);
	const bool dualLayerSubmodule = (bool)mxGetScalar(prhs[49]);
	const bool indexBased = (bool)mxGetScalar(prhs[50]);

	const int64_t imDim = Nx * Ny * Nz;

	const bool dynamic = Nt > 1;

	// Count inputs and check for char type
	if (!mxIsChar(prhs[0]))
		mexErrMsgTxt("Input argument is not char");

	/* Pointer to character array */
	const char* argv = mxArrayToString(prhs[0]);

	TChain *Coincidences = new TChain("Coincidences");
	Coincidences->Add(argv);

	TChain *delay = nullptr;
	int64_t Ndelays = 0LL;
	//bool dynamic = false;

	if (randoms_correction) {
		delay = new TChain("delay");
		delay->Add(argv);
		Ndelays = delay->GetEntries();
	}


	int64_t Nentries = Coincidences->GetEntries();
	delete Coincidences;
	if (randoms_correction)
		delete delay;
	if (source) {
		plhs[5] = mxCreateNumericMatrix(imDim * Nt, 1, mxUINT16_CLASS, mxREAL);
		plhs[6] = mxCreateNumericMatrix(imDim * Nt, 1, mxUINT16_CLASS, mxREAL);
		plhs[7] = mxCreateNumericMatrix(imDim * Nt, 1, mxUINT16_CLASS, mxREAL);
	}
	else {
		plhs[5] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		plhs[6] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		plhs[7] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
	}

	if (store_coordinates) {
		plhs[8] = mxCreateNumericMatrix(Nentries, 1, mxUINT16_CLASS, mxREAL);
		plhs[9] = mxCreateNumericMatrix(6, Nentries, mxSINGLE_CLASS, mxREAL);
		if (randoms_correction)
			plhs[10] = mxCreateNumericMatrix(6, Ndelays, mxSINGLE_CLASS, mxREAL);
		else
			plhs[10] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
		//plhs[10] = mxCreateNumericMatrix(Nentries, 1, mxSINGLE_CLASS, mxREAL);
		//plhs[11] = mxCreateNumericMatrix(Nentries, 1, mxSINGLE_CLASS, mxREAL);
		//plhs[12] = mxCreateNumericMatrix(Nentries, 1, mxSINGLE_CLASS, mxREAL);
		//plhs[13] = mxCreateNumericMatrix(Nentries, 1, mxSINGLE_CLASS, mxREAL);
		//plhs[14] = mxCreateNumericMatrix(Nentries, 1, mxSINGLE_CLASS, mxREAL);
	}
	else {
		plhs[8] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		plhs[9] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
		plhs[10] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
		//plhs[10] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
		//plhs[11] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
		//plhs[12] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
		//plhs[13] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
		//plhs[14] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	}
	if (indexBased) {
		plhs[11] = mxCreateNumericMatrix(2, Nentries, mxUINT16_CLASS, mxREAL);
		plhs[12] = mxCreateNumericMatrix(2, Nentries, mxUINT16_CLASS, mxREAL);
		if (randoms_correction) {
			plhs[13] = mxCreateNumericMatrix(2, Ndelays, mxUINT16_CLASS, mxREAL);
			plhs[14] = mxCreateNumericMatrix(2, Ndelays, mxUINT16_CLASS, mxREAL);
		}
		else {
			plhs[13] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
			plhs[14] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		}
	}
	else {
		plhs[11] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		plhs[12] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		plhs[13] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
		plhs[14] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);
	}

	/* Assign pointers to the various parameters */
	uint16_t * tIndex = (uint16_t*)mxGetData(plhs[8]);
	uint16_t* S = nullptr;
	uint16_t* SC = nullptr;
	uint16_t* RA = nullptr;
	uint16_t* trIndex = nullptr;
	uint16_t* axIndex = nullptr;
	uint16_t* DtrIndex = nullptr;
	uint16_t* DaxIndex = nullptr;
	//float* x1 = nullptr, * x2 = nullptr, * y1 = nullptr, * y2 = nullptr, * z1 = nullptr, * z2 = nullptr;
	float* coord = nullptr, * Dcoord = nullptr;
	if (source) {
		S = (uint16_t*)mxGetData(plhs[5]);
		SC = (uint16_t*)mxGetData(plhs[6]);
		RA = (uint16_t*)mxGetData(plhs[7]);
	}
	if (store_coordinates) {
		coord = (float*)mxGetData(plhs[9]);
		if (randoms_correction)
			Dcoord = (float*)mxGetData(plhs[10]);
		//x2 = (float*)mxGetData(plhs[10]);
		//y1 = (float*)mxGetData(plhs[11]);
		//y2 = (float*)mxGetData(plhs[12]);
		//z1 = (float*)mxGetData(plhs[13]);
		//z2 = (float*)mxGetData(plhs[14]);
	}
	if (indexBased) {
		trIndex = (uint16_t*)mxGetData(plhs[11]);
		axIndex = (uint16_t*)mxGetData(plhs[12]);
		DtrIndex = (uint16_t*)mxGetData(plhs[13]);
		DaxIndex = (uint16_t*)mxGetData(plhs[14]);
	}
	plhs[0] = mxCreateSharedDataCopy(prhs[31]);
	plhs[1] = mxCreateSharedDataCopy(prhs[32]);
	plhs[2] = mxCreateSharedDataCopy(prhs[33]);
	plhs[3] = mxCreateSharedDataCopy(prhs[34]);
	plhs[4] = mxCreateSharedDataCopy(prhs[35]);
	uint16_t* Sino = (uint16_t*)mxGetData(plhs[0]);
	uint16_t* SinoT = (uint16_t*)mxGetData(plhs[1]);
	uint16_t* SinoC = (uint16_t*)mxGetData(plhs[2]);
	uint16_t* SinoR = (uint16_t*)mxGetData(plhs[3]);
	uint16_t* SinoD = (uint16_t*)mxGetData(plhs[4]);

	const float matlabPtr = NULL;

	histogram(argv, tPoints, alku, loppu, source, linear_multp, cryst_per_block, blocks_per_ring, det_per_ring, S, SC, RA, trIndex, axIndex, DtrIndex, DaxIndex, obtain_trues, store_scatter, store_randoms,
		scatter_components, randoms_correction, coord, Dcoord, store_coordinates, dynamic, cryst_per_block_z, transaxial_multip, rings, sinoSize, Ndist, Nang, ringDifference, span,
		seg, Nt, TOFSize, nDistSide, Sino, SinoT, SinoC, SinoR, SinoD, detWPseudo, nPseudos, binSize, FWHM, verbose, nLayers, dx, dy, dz, bx, by, bz, Nx, Ny, Nz, dualLayerSubmodule, imDim, indexBased, tIndex, matlabPtr);


	mexEvalString("pause(.001);");
	gROOT->Reset();

	return;
}
