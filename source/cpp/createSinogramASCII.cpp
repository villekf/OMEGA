/**************************************************************************
* Constructs a 3D/4D/5D sinogram from the input ring number and positions.
* This code uses the old C MEX API. For MATLAB 2017b and older.
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

#include "mexFunktio.h"
#include "saveSinogram.h"
extern "C" mxArray * mxCreateSharedDataCopy(const mxArray * pr);



void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])

{


	/* Check for proper number of arguments */

	if (nrhs != 30) {
		mexErrMsgIdAndTxt("MATLAB:list2matlab_aivi:invalidNumInputs",
			"30 input arguments required.");
	}
	else if (nlhs > 4) {
		mexErrMsgIdAndTxt("MATLAB:list2matlab_aivi:maxlhs",
			"Too many output arguments.");
	}

	int ind = 0;
	const uint16_t* ringPos1 = getUint16s(prhs[ind], "solu");
	ind++;
	const uint16_t* ringPos2 = getUint16s(prhs[ind], "solu");
	ind++;
	const uint16_t* ringNumber1 = getUint16s(prhs[ind], "solu");
	ind++;
	const uint16_t* ringNumber2 = getUint16s(prhs[ind], "solu");
	ind++;
	const bool* truesIndex = getBools(prhs[ind], "solu");
	ind++;
	const bool* scatterIndex = getBools(prhs[ind], "solu");
	ind++;
	const bool* randomsIndex = getBools(prhs[ind], "solu");
	ind++;
	const uint64_t sinoSize = getScalarUInt64(prhs[ind], ind);
	ind++;
	const uint32_t Ndist = getScalarUInt32(prhs[ind], ind);
	ind++;
	const uint32_t Nang = getScalarUInt32(prhs[ind], ind);
	ind++;
	const uint32_t ringDifference = getScalarUInt32(prhs[ind], ind);
	ind++;
	const uint32_t span = getScalarUInt32(prhs[ind], ind);
	ind++;
	const uint32_t* seg = getUint32s(prhs[ind], "solu");
	const size_t pituus = mxGetNumberOfElements(prhs[ind]);
	ind++;
	const uint64_t TOFSize = getScalarUInt64(prhs[ind], ind);
	ind++;
	const uint16_t* time = getUint16s(prhs[ind], "solu");
	ind++;
	const uint64_t Nt = getScalarUInt64(prhs[ind], ind);
	ind++;
	const int32_t detPerRing = getScalarInt32(prhs[ind], ind);
	ind++;
	const int32_t rings = getScalarInt32(prhs[ind], ind);
	ind++;
	const uint16_t* bins = getUint16s(prhs[ind], "solu");
	ind += 5;
	const int32_t nDistSide = getScalarInt32(prhs[ind], ind);
	ind++;
	const int32_t detWPseudo = getScalarInt32(prhs[ind], ind);
	ind++;
	const int32_t nPseudos = getScalarInt32(prhs[ind], ind);
	ind++;
	const int32_t crystPerBlock = getScalarInt32(prhs[ind], ind);
	ind++;
	const int32_t nLayers = getScalarInt32(prhs[ind], ind);
	ind++;
	const uint8_t* layer1 = getUint8s(prhs[ind], "solu");
	ind++;
	const uint8_t* layer2 = getUint8s(prhs[ind], "solu");
	ind++;

	bool storeTrues = false;
	bool storeScatter = false;
	bool storeRandoms = false;
	const int64_t koko = mxGetNumberOfElements(prhs[1]);
	if (mxGetNumberOfElements(prhs[4]) > 1) {
		storeTrues = true;
	}
	if (mxGetNumberOfElements(prhs[5]) > 1) {
		storeScatter = true;
	}
	if (mxGetNumberOfElements(prhs[6]) > 1) {
		storeRandoms = true;
	}


	/* Assign pointers to the various parameters */
	plhs[0] = mxCreateSharedDataCopy(prhs[19]);
	plhs[1] = mxCreateSharedDataCopy(prhs[20]);
	plhs[2] = mxCreateSharedDataCopy(prhs[21]);
	plhs[3] = mxCreateSharedDataCopy(prhs[22]);
	//getUint16s(plhs[0], "solu");
	uint16_t* Sino = getUint16s(plhs[0], "solu");
	uint16_t* SinoT = getUint16s(plhs[1], "solu");
	uint16_t* SinoC = getUint16s(plhs[2], "solu");
	uint16_t* SinoR = getUint16s(plhs[3], "solu");

	openMPSino(ringPos1, ringPos2, ringNumber1, ringNumber2, truesIndex, scatterIndex, randomsIndex, sinoSize, Ndist, Nang, ringDifference,
		span, seg, time, Nt, TOFSize, Sino, SinoT, SinoC, SinoR, storeTrues, storeScatter, storeRandoms, detPerRing, rings, koko, 
		bins, nDistSide, pituus, detWPseudo, nPseudos, crystPerBlock, layer1, layer2, nLayers);

	return;

}
