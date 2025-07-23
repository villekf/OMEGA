/**************************************************************************
* Sinogram creation function for OMEGA. This function uses the provided
* ring number and position values to compute the corresponding sinogram
* index value. This index value is the output of this function.
*
* Copyright (C) 2020-2025 Ville-Veikko Wettenhovi
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
#pragma once
#include "dIndices.h"
#include <thread>

inline void mexPrint(const char* str) {
//#ifdef MATLAB
//	mexPrintf("%s\n", str);
//	mexEvalString("pause(.0001);");
//#else
	fprintf(stdout, "%s\n", str);
	fflush(stdout);
//#endif
}


template<typename T>
int64_t saveSinogram(const int32_t ring_pos1, const int32_t ring_pos2, int32_t ring_number1, int32_t ring_number2, const uint64_t sinoSize, const uint32_t Ndist, 
	const uint32_t Nang, const uint32_t ring_difference, const uint32_t span, const T* seg, const uint64_t TOFSize,
	const int32_t det_per_ring, const int32_t rings, const uint64_t bins, const int32_t nDistSide, bool& swap, uint16_t tPoint = 0, int32_t layer = 0,
	int64_t nLayers = 1) {
	// Initial TOF bin
	uint64_t binN = 0ULL;
	// Initial index value, negative value means the LOR does not belong to the sinogram FOV
	int64_t indeksi = -1;
	const int32_t xa = std::max(ring_pos1, ring_pos2);
	const int32_t ya = std::min(ring_pos1, ring_pos2);
	int32_t j = ((xa + ya + det_per_ring / 2 + 1) % det_per_ring) / 2;
	const int32_t b = j + det_per_ring / 2;
	int32_t i = std::abs(xa - ya - det_per_ring / 2);
	const bool ind = ya < j || b < xa;
	if (ind)
		i = -i;
	// Determine if the index belongs to the "swap corners"
	swap = (j * 2) < -i || i <= ((j - det_per_ring / 2) * 2);
	// Determine if the index is within the sinogram FOV
	bool accepted_lors;
	if (Ndist % 2U == 0)
		accepted_lors = (i <= (static_cast<int32_t>(Ndist) / 2 + std::min(0, nDistSide)) && i >= (-static_cast<int32_t>(Ndist) / 2 + std::max(0, nDistSide)));
	else
		accepted_lors = (i <= static_cast<int32_t>(Ndist) / 2 && i >= (-static_cast<int32_t>(Ndist) / 2));
	accepted_lors = accepted_lors && (std::abs(ring_number1 - ring_number2) <= ring_difference);
	// LOR is accepted
	if (accepted_lors) {
		int32_t sinoIndex = 0;
		// Determine the TOF bin
		if (TOFSize > sinoSize) {
			binN = bins;
			if (swap && binN != 0ULL) {
				if (binN % 2 == 0)
					binN--;
				else
					binN++;
			}
		}
		j = j / (det_per_ring / 2 / Nang);
		if (swap) {
			const int32_t ring_number3 = ring_number1;
			ring_number1 = ring_number2;
			ring_number2 = ring_number3;
		}
		i = i + Ndist / 2 - std::max(0, nDistSide);
		const bool swappi = ring_pos2 > ring_pos1;
		if (swappi) {
			const int32_t ring_number3 = ring_number1;
			ring_number1 = ring_number2;
			ring_number2 = ring_number3;
		}
		//// Simple index when span = 1
		if (span <= 1) {
			sinoIndex = ring_number2 * rings + ring_number1;
		}
		else {
			const int32_t erotus = ring_number1 - ring_number2;
			const int32_t summa = ring_number1 + ring_number2;
			if (std::abs(erotus) <= span / 2) {
				sinoIndex = summa;
			}
			else {
				sinoIndex = ((std::abs(erotus) + (span / 2)) / span);
				if (erotus < 0) {
					sinoIndex = (summa - ((span / 2) * (sinoIndex * 2 - 1) + sinoIndex)) + static_cast<uint32_t>(seg[(sinoIndex - 1) * 2]);
				}
				else {
					sinoIndex = (summa - ((span / 2) * (sinoIndex * 2 - 1) + sinoIndex)) + static_cast<uint32_t>(seg[(sinoIndex - 1) * 2 + 1]);
				}
			}
		}
		// Final index number
		indeksi = static_cast<int64_t>(i) + static_cast<int64_t>(j) * static_cast<int64_t>(Ndist) +
			static_cast<int64_t>(sinoIndex) * static_cast<int64_t>(Ndist) * static_cast<int64_t>(Nang) + sinoSize * static_cast<int64_t>(layer)
			+ sinoSize * binN * nLayers * nLayers + TOFSize * static_cast<int64_t>(tPoint) * nLayers * nLayers;
	}
	return indeksi;
}


template<typename T, typename K, typename C, typename B>
void openMPSino(const T* ringPos1, const T* ringPos2, const T* ringNumber1, const T* ringNumber2,
	const B* trues_index, const B* scatter_index, const B* randoms_index, const uint64_t sinoSize, const uint32_t Ndist,
	const uint32_t Nang, const uint32_t ring_difference, const uint32_t span, const K* seg, const T* time, const uint64_t Nt,
	const uint64_t TOFSize, T* Sino, T* SinoT, T* SinoC, T* SinoR, const bool store_trues, const bool store_scatter, 
	const bool store_randoms, const int32_t det_per_ring, const int32_t rings, const int64_t koko, const T* bins, 
	const int32_t nDistSide, const size_t pituus, const int32_t detWPseudo, const int32_t nPseudos,
	const int32_t cryst_per_block, const C* layer1, const C* layer2, const int32_t nLayers) {

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
#pragma omp parallel for schedule(monotonic:dynamic)
#else
#pragma omp parallel for schedule(dynamic)
#endif
#endif
	for (int64_t kk = 0; kk < koko; kk++) {
		uint16_t tPoint = 0.;
		if (Nt > 1)
			tPoint = time[kk];
		uint64_t binN = 0ULL;
		if (TOFSize > sinoSize)
			binN = static_cast<uint64_t>(bins[kk]);
		int32_t ring_pos1 = static_cast<int32_t>(ringPos1[kk]);
		int32_t ring_pos2 = static_cast<int32_t>(ringPos2[kk]);
		int32_t ring_number1 = static_cast<int32_t>(ringNumber1[kk]);
		int32_t ring_number2 = static_cast<int32_t>(ringNumber2[kk]);
		int32_t layer = 0;
		if (nLayers > 1) {
			if (layer1[kk] == 1 && layer2[kk] == 1)
				layer = 3;
			else if (layer1[kk] == 1 && layer2[kk] == 0)
				layer = 1;
			else if (layer1[kk] == 0 && layer2[kk] == 1)
				layer = 2;
			if (nLayers > 2) {
				if (layer1[kk] == 2 && layer2[kk] == 2)
					layer = 8;
				else if (layer1[kk] == 2 && layer2[kk] == 0)
					layer = 4;
				else if (layer1[kk] == 0 && layer2[kk] == 2)
					layer = 5;
				else if (layer1[kk] == 2 && layer2[kk] == 1)
					layer = 6;
				else if (layer1[kk] == 1 && layer2[kk] == 2)
					layer = 7;
			}
		}
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
		const int64_t indeksi = saveSinogram(ring_pos1, ring_pos2, ring_number1, ring_number2, sinoSize, Ndist, Nang, ring_difference, span, seg, TOFSize,
			detWPseudo, rings, binN, nDistSide, swap, tPoint, layer, nLayers);
		if (indeksi >= TOFSize * nLayers * Nt) {
			mexPrint("Sinogram index is larger than the maximum possible! Aborting!");
			continue;
		}
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