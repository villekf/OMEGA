/**************************************************************************
* Sinogram creation function for OMEGA. This function uses the provided
* ring number and position values to compute the corresponding sinogram
* index value. This index value is the output of this function.
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
#pragma once
#include "dIndices.h"


template<typename T>
int64_t saveSinogram(const int32_t ring_pos1, const int32_t ring_pos2, int32_t ring_number1, int32_t ring_number2, const uint64_t sinoSize, const uint32_t Ndist, 
	const uint32_t Nang, const uint32_t ring_difference, const uint32_t span, const T* seg, const double time, const uint64_t NT, const uint64_t TOFSize,
	const double vali, const double alku, const int32_t det_per_ring, const int32_t rings, const uint64_t bins, const int32_t nDistSide, bool& swap) {
	// Initial time step
	uint64_t tPoint = 0ULL;
	// Initial TOF bin
	uint64_t binN = 0ULL;
	// Initial index value, negative value means the LOR does not belong to the sinogram FOV
	int64_t indeksi = -1;
	// Compute the real time step if using dynamic data
	if (NT > 1) {
		tPoint = static_cast<uint64_t>(std::floor((time - alku) / vali));
	}
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
		// Simple index when span = 1
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
			static_cast<int64_t>(sinoIndex) * static_cast<int64_t>(Ndist) * static_cast<int64_t>(Nang) + sinoSize * binN + TOFSize * tPoint;
	}
	return indeksi;
}