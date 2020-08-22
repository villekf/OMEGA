/**************************************************************************
* Constructs a 3D/4D/5D sinogram from the input ring number and positions.
* This code uses the Octave API.
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
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif
#include <octave/oct.h>
#include "saveSinogram.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <thread>


void openMPSino(const octave_uint16* ringPos1, const octave_uint16* ringPos2, const octave_uint16* ringNumber1, const octave_uint16* ringNumber2,
	const bool* trues_index, const bool* scatter_index, const bool* randoms_index, const uint64_t sinoSize, const uint32_t Ndist, 
	const uint32_t Nang, const uint32_t ring_difference, const uint32_t span, const octave_uint32* seg, const double* time, const uint64_t NT,
	const uint64_t TOFSize, const double vali, const double alku, uint16_t* Sino, uint16_t* SinoT, uint16_t* SinoC, uint16_t* SinoR,
	const bool store_trues, const bool store_scatter, const bool store_randoms, const int32_t det_per_ring, const int32_t rings, 
	const int64_t koko, const octave_uint16* bins, const int32_t nDistSide, const size_t pituus, const int32_t detWPseudo, const int32_t nPseudos,
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

#pragma omp parallel for
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
		if (pseudoD) {
			ring_pos1 += ring_pos1 / cryst_per_block;
			ring_pos2 += ring_pos2 / cryst_per_block;
		}
		if (pseudoR) {
			ring_number1 += ring_number1 / gapSize;
			ring_number2 += ring_number2 / gapSize;
		}
		const int64_t indeksi = saveSinogram(ring_pos1, ring_pos2, ring_number1, ring_number2, sinoSize, Ndist, Nang, ring_difference, span, seg, aika, NT, TOFSize,
			vali, alku, detWPseudo, rings, binN, nDistSide);
		if (indeksi >= 0) {
			if (store_trues && trues_index[kk]) {
#pragma omp atomic
				SinoT[indeksi]++;
			}
			else if (store_scatter && scatter_index[kk]) {
#pragma omp atomic
				SinoC[indeksi]++;
			}
			else if (store_randoms && randoms_index[kk]) {
#pragma omp atomic
				SinoR[indeksi]++;
			}
#pragma omp atomic
			Sino[indeksi]++;
		}
	}
}

DEFUN_DLD(createSinogramASCIIOct, prhs, nargout, "ASCII to sinogram help") {


	/* Check for proper number of arguments */

	const double vali = prhs(0).scalar_value();
	const uint16NDArray ringPos1 = prhs(1).uint16_array_value();
	const uint16NDArray ringPos2 = prhs(2).uint16_array_value();
	const uint16NDArray ringNumber1 = prhs(3).uint16_array_value();
	const uint16NDArray ringNumber2 = prhs(4).uint16_array_value();
	const boolNDArray truesIndex = prhs(5).bool_array_value();
	const boolNDArray scatterIndex = prhs(6).bool_array_value();
	const boolNDArray randomsIndex = prhs(7).bool_array_value();
	const uint64_t sinoSize = prhs(8).uint64_scalar_value();
	const uint32_t Ndist = prhs(9).uint32_scalar_value();
	const uint32_t Nang = prhs(10).uint32_scalar_value();
	const uint32_t ringDifference = prhs(11).uint32_scalar_value();
	const uint32_t span = prhs(12).uint32_scalar_value();
	const uint32NDArray seg = prhs(13).uint32_array_value();
	const size_t pituus = seg.numel();
	const uint64_t TOFSize = prhs(14).uint64_scalar_value();
	const NDArray time = prhs(15).array_value();
	const uint64_t NT = prhs(16).uint64_scalar_value();
	const double alku = prhs(17).scalar_value();
	const int32_t detPerRing = prhs(18).int32_scalar_value();
	const int32_t rings = prhs(19).int32_scalar_value();
	const uint16NDArray bins = prhs(20).uint16_array_value();
	const int32_t nDistSide = prhs(25).int32_scalar_value();
	const int32_t detWPseudo = prhs(26).int32_scalar_value();
	const int32_t nPseudos = prhs(27).int32_scalar_value();
	const int32_t crystPerBlock = prhs(28).int32_scalar_value();

	bool storeTrues = false;
	bool storeScatter = false;
	bool storeRandoms = false;
	const int64_t koko = ringPos1.numel();
	if (truesIndex.numel() > 1) {
		storeTrues = true;
	}
	if (scatterIndex.numel() > 1) {
		storeScatter = true;
	}
	if (randomsIndex.numel() > 1) {
		storeRandoms = true;
	}

	const octave_uint16* ring_pos1 = ringPos1.fortran_vec();
	const octave_uint16* ring_pos2 = ringPos2.fortran_vec();
	const octave_uint16* ring_number1 = ringNumber1.fortran_vec();
	const octave_uint16* ring_number2 = ringNumber2.fortran_vec();

	const octave_uint32* seg_p = seg.fortran_vec();
	const octave_uint16* bins_p = bins.fortran_vec();
	const double* time_p = time.fortran_vec();
	const bool* tIndex = truesIndex.fortran_vec();
	const bool* cIndex = scatterIndex.fortran_vec();
	const bool* rIndex = randomsIndex.fortran_vec();


	/* Assign pointers to the various parameters */
	uint16NDArray SinoO = prhs(21).uint16_array_value();
	uint16NDArray SinoTO = prhs(22).uint16_array_value();
	uint16NDArray SinoCO = prhs(23).uint16_array_value();
	uint16NDArray SinoRO = prhs(24).uint16_array_value();
	uint16_t* Sino = reinterpret_cast<uint16_t*>(SinoO.fortran_vec());
	uint16_t* SinoT = reinterpret_cast<uint16_t*>(SinoTO.fortran_vec());
	uint16_t* SinoC = reinterpret_cast<uint16_t*>(SinoCO.fortran_vec());
	uint16_t* SinoR = reinterpret_cast<uint16_t*>(SinoRO.fortran_vec());

	openMPSino(ring_pos1, ring_pos2, ring_number1, ring_number2, tIndex, cIndex, rIndex, sinoSize, Ndist, Nang, ringDifference,
		span, seg_p, time_p, NT, TOFSize, vali, alku, Sino, SinoT, SinoC, SinoR, storeTrues, storeScatter, storeRandoms, detPerRing, rings, koko, 
		bins_p, nDistSide, pituus, detWPseudo, nPseudos, crystPerBlock);

	octave_value_list retval(nargout);

	retval(0) = octave_value(SinoO);
	retval(1) = octave_value(SinoTO);
	retval(2) = octave_value(SinoCO);
	retval(3) = octave_value(SinoRO);

	return retval;

}
