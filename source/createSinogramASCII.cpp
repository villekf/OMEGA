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
#include <thread>
extern "C" mxArray * mxCreateSharedDataCopy(const mxArray * pr);


void openMPSino(const uint16_t* ringPos1, const uint16_t* ringPos2, const uint16_t* ringNumber1, const uint16_t* ringNumber2,
	const bool* trues_index, const bool* scatter_index, const bool* randoms_index, const uint64_t sinoSize, const uint32_t Ndist, 
	const uint32_t Nang, const uint32_t ring_difference, const uint32_t span, const uint32_t* seg, const double* time, const uint64_t NT, 
	const uint64_t TOFSize, const double vali, const double alku, uint16_t* Sino, uint16_t* SinoT, uint16_t* SinoC, uint16_t* SinoR, 
	const bool store_trues, const bool store_scatter, const bool store_randoms, const int32_t det_per_ring, const int32_t rings, 
	const int64_t koko, const uint16_t* bins, const int32_t nDistSide, const size_t pituus, const int32_t detWPseudo, const int32_t nPseudos, 
	const int32_t cryst_per_block) {

#ifdef _OPENMP
	if (omp_get_max_threads() == 1) {
		int n_threads = std::thread::hardware_concurrency();
		omp_set_num_threads(n_threads);
	}
#endif
	const bool pseudoD = detWPseudo > det_per_ring;
	const bool pseudoR = nPseudos > 0;
	int32_t gapSize = 0;
	if (pseudoR) {
		gapSize = rings / (nPseudos + 1);
	}
#ifdef _OPENMP
#if _OPENMP >= 201511
#pragma omp parallel for schedule(monotonic:dynamic, nChunks)
#else
#pragma omp parallel for schedule(dynamic, nChunks)
#endif
#endif
	for (int64_t kk = 0; kk < koko; kk++) {
		double aika = 0.;
		if (NT > 1)
			aika = time[kk];
		uint64_t binN = 0ULL;
		if (TOFSize > sinoSize)
			binN = static_cast<uint64_t>(bins[kk]);
		int32_t ring_pos1 = static_cast<int32_t>(ringPos1[kk]);
		int32_t ring_pos2 = static_cast<int32_t>(ringPos2[kk]);
		int32_t ring_number1 = static_cast<int32_t>(ringNumber1[kk]);
		int32_t ring_number2 = static_cast<int32_t>(ringNumber2[kk]);
		// Pseudo detectors
		if (pseudoD) {
			ring_pos1 += ring_pos1 / cryst_per_block;
			ring_pos2 += ring_pos2 / cryst_per_block;
		}
		// Pseudo rings
		if (pseudoR) {
			ring_number1 += ring_number1 / gapSize;
			ring_number2 += ring_number2 / gapSize;
		}
		bool swap = false;
		const int64_t indeksi = saveSinogram(ring_pos1, ring_pos2, ring_number1, ring_number2, sinoSize, Ndist, Nang, ring_difference, span, seg, aika, NT, TOFSize,
			vali, alku, detWPseudo, rings, binN, nDistSide, swap);
		if (indeksi >= 0) {
			// Trues
			if (store_trues && trues_index[kk]) {
#pragma omp atomic
				SinoT[indeksi]++;
			}
			// Scatter
			else if (store_scatter && scatter_index[kk]) {
#pragma omp atomic
				SinoC[indeksi]++;
			}
			// Randoms
			else if (store_randoms && randoms_index[kk]) {
#pragma omp atomic
				SinoR[indeksi]++;
			}
			// Prompts
#pragma omp atomic
			Sino[indeksi]++;
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])

{


	/* Check for proper number of arguments */

	if (nrhs != 29) {
		mexErrMsgIdAndTxt("MATLAB:list2matlab_aivi:invalidNumInputs",
			"29 input arguments required.");
	}
	else if (nlhs > 4) {
		mexErrMsgIdAndTxt("MATLAB:list2matlab_aivi:maxlhs",
			"Too many output arguments.");
	}

	int ind = 0;
	const double vali = getScalarDouble(prhs[ind], ind);
	ind++;
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const uint16_t* ringPos1 = (uint16_t*)mxGetUint16s(prhs[ind]);
#else
	const uint16_t* ringPos1 = (uint16_t*)mxGetData(prhs[ind]);
#endif
	ind++;
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const uint16_t* ringPos2 = (uint16_t*)mxGetUint16s(prhs[ind]);
#else
	const uint16_t* ringPos2 = (uint16_t*)mxGetData(prhs[ind]);
#endif
	ind++;
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const uint16_t* ringNumber1 = (uint16_t*)mxGetUint16s(prhs[ind]);
#else
	const uint16_t* ringNumber1 = (uint16_t*)mxGetData(prhs[ind]);
#endif
	ind++;
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const uint16_t* ringNumber2 = (uint16_t*)mxGetUint16s(prhs[ind]);
#else
	const uint16_t* ringNumber2 = (uint16_t*)mxGetData(prhs[ind]);
#endif
	ind++;
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const bool* truesIndex = (bool*)mxGetLogicals(prhs[ind]);
#else
	const bool* truesIndex = (bool*)mxGetData(prhs[ind]);
#endif
	ind++;
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const bool* scatterIndex = (bool*)mxGetLogicals(prhs[ind]);
#else
	const bool* scatterIndex = (bool*)mxGetData(prhs[ind]);
#endif
	ind++;
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const bool* randomsIndex = (bool*)mxGetLogicals(prhs[ind]);
#else
	const bool* randomsIndex = (bool*)mxGetData(prhs[ind]);
#endif
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
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const uint32_t* seg = (uint32_t*)mxGetUint32s(prhs[ind]);
#else
	const uint32_t* seg = (uint32_t*)mxGetData(prhs[ind]);
#endif
	const size_t pituus = mxGetNumberOfElements(prhs[ind]);
	ind++;
	const uint64_t TOFSize = getScalarUInt64(prhs[ind], ind);
	ind++;
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const double* time = (double*)mxGetDoubles(prhs[ind]);
#else
	const double* time = (double*)mxGetData(prhs[ind]);
#endif
	ind++;
	const uint64_t NT = getScalarUInt64(prhs[ind], ind);
	ind++;
	const double alku = getScalarDouble(prhs[ind], ind);
	ind++;
	const int32_t detPerRing = getScalarInt32(prhs[ind], ind);
	ind++;
	const int32_t rings = getScalarInt32(prhs[ind], ind);
	ind++;
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	const uint16_t* bins = (uint16_t*)mxGetUint16s(prhs[ind]);
#else
	const uint16_t* bins = (uint16_t*)mxGetData(prhs[ind]);
#endif
	ind++;
	const int32_t nDistSide = getScalarInt32(prhs[ind], ind);
	ind++;
	const int32_t detWPseudo = getScalarInt32(prhs[ind], ind);
	ind++;
	const int32_t nPseudos = getScalarInt32(prhs[ind], ind);
	ind++;
	const int32_t crystPerBlock = getScalarInt32(prhs[ind], ind);
	ind++;

	bool storeTrues = false;
	bool storeScatter = false;
	bool storeRandoms = false;
	const int64_t koko = mxGetNumberOfElements(prhs[1]);
	if (mxGetNumberOfElements(prhs[5]) > 1) {
		storeTrues = true;
	}
	if (mxGetNumberOfElements(prhs[6]) > 1) {
		storeScatter = true;
	}
	if (mxGetNumberOfElements(prhs[7]) > 1) {
		storeRandoms = true;
	}


	/* Assign pointers to the various parameters */
	plhs[0] = mxCreateSharedDataCopy(prhs[21]);
	plhs[1] = mxCreateSharedDataCopy(prhs[22]);
	plhs[2] = mxCreateSharedDataCopy(prhs[23]);
	plhs[3] = mxCreateSharedDataCopy(prhs[24]);
#ifdef MX_HAS_INTERLEAVED_COMPLEX
	uint16_t* Sino = (uint16_t*)mxGetUint16s(plhs[0]);
	uint16_t* SinoT = (uint16_t*)mxGetUint16s(plhs[1]);
	uint16_t* SinoC = (uint16_t*)mxGetUint16s(plhs[2]);
	uint16_t* SinoR = (uint16_t*)mxGetUint16s(plhs[3]);
#else
	uint16_t* Sino = (uint16_t*)mxGetData(plhs[0]);
	uint16_t* SinoT = (uint16_t*)mxGetData(plhs[1]);
	uint16_t* SinoC = (uint16_t*)mxGetData(plhs[2]);
	uint16_t* SinoR = (uint16_t*)mxGetData(plhs[3]);
#endif

	openMPSino(ringPos1, ringPos2, ringNumber1, ringNumber2, truesIndex, scatterIndex, randomsIndex, sinoSize, Ndist, Nang, ringDifference,
		span, seg, time, NT, TOFSize, vali, alku, Sino, SinoT, SinoC, SinoR, storeTrues, storeScatter, storeRandoms, detPerRing, rings, koko, 
		bins, nDistSide, pituus, detWPseudo, nPseudos, crystPerBlock);

	return;

}
