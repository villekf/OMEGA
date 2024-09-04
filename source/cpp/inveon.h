#pragma once

#include <cstdint>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#ifdef MATLAB
#include "mex.h"
#endif

FILE* streami;
#define DET_PER_RING 320
#define RINGS 80

void saveSinogram(uint16_t L1, uint16_t L2, uint16_t* Sino, const uint32_t Ndist, const uint32_t Nang, const uint32_t ring_difference, const uint32_t span,
	const uint64_t sinoSize, const uint32_t* seg, const int32_t nDistSide, const int tPoint = 0) {
	int32_t ring_pos1 = L1 % DET_PER_RING;
	int32_t ring_pos2 = L2 % DET_PER_RING;
	int32_t ring_number1 = L1 / DET_PER_RING;
	int32_t ring_number2 = L2 / DET_PER_RING;
	const int32_t xa = std::max(ring_pos1, ring_pos2);
	const int32_t ya = std::min(ring_pos1, ring_pos2);
	int32_t j = ((xa + ya + DET_PER_RING / 2 + 1) % DET_PER_RING) / 2;
	const int32_t b = j + DET_PER_RING / 2;
	int32_t i = std::abs(xa - ya - DET_PER_RING / 2);
	const bool ind = ya < j || b < xa;
	if (ind)
		i = -i;
	const bool swap = (j * 2) < -i || i <= ((j - DET_PER_RING / 2) * 2);
	bool accepted_lors;
	if (Ndist % 2U == 0)
		accepted_lors = (i <= (static_cast<int32_t>(Ndist) / 2 + std::min(0, nDistSide)) && i >= (-static_cast<int32_t>(Ndist) / 2 + std::max(0, nDistSide)));
	else
		accepted_lors = (i <= static_cast<int32_t>(Ndist) / 2 && i >= (-static_cast<int32_t>(Ndist) / 2));
	accepted_lors = accepted_lors && (std::abs(ring_number1 - ring_number2) <= ring_difference);
	int32_t sinoIndex = 0;
	if (accepted_lors) {
		j = j / (DET_PER_RING / 2 / Nang);
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
		if (span <= 1) {
			sinoIndex = ring_number2 * RINGS + ring_number1;
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
		const uint64_t indeksi = static_cast<uint64_t>(i) + static_cast<uint64_t>(j) * static_cast<uint64_t>(Ndist) +
			static_cast<uint64_t>(sinoIndex) * static_cast<uint64_t>(Ndist) * static_cast<uint64_t>(Nang) + sinoSize * static_cast<uint64_t>(tPoint);
		Sino[indeksi]++;
	}
	else
		return;
}

int histogram(uint16_t* LL1, uint16_t* LL2, uint16_t* tpoints, const char* argv, const double* vali, const double alku, const double loppu, 
	const uint32_t detectors, const size_t pituus, const bool randoms_correction, uint16_t* DD1, uint16_t* DD2, uint16_t* Sino, uint16_t* SinoD, const bool saveRawData,
	const uint32_t Ndist, const uint32_t Nang, const uint32_t ringDifference, const uint32_t span, const uint64_t sinoSize, const uint32_t* seg, const int32_t nDistSide,
	const bool storeCoordinates, const uint64_t Nt)
{

	static uint64_t i;

	static char* in_file;

	static int16_t qb;
	static int prompt;
	static uint64_t ew1;
	static int tag;

	double ms = 0;		// seconds

#if (defined(WIN32) || defined(_WIN32) || (defined(__WIN32) && !defined(__CYGWIN__)) || defined(_WIN64)) && defined(_MSC_VER)
	errno_t err;
	err = fopen_s(&streami, argv, "rb");
	if (err != 0) {
#ifdef MATLAB
		mexErrMsgIdAndTxt("MATLAB:inveon_list2matlab:invalidFile",
			"Error opening file or no file opened");
#else
		fprintf(stdout, "Error opening file or no file opened\n");
		fflush(stdout);
#endif
		return 0;
	}
#else
	streami = fopen(argv, "rb");
	if (streami == NULL) {
#ifdef MATLAB
		mexErrMsgIdAndTxt("MATLAB:inveon_list2matlab:invalidFile",
			"Error opening file or no file opened");
#else
		fprintf(stdout, "Error opening file or no file opened\n");
		fflush(stdout);
#endif
		return 0;
	}
#endif

	uint32_t L1;    //detector 1
	uint32_t L2;    //detector 2
	int mscount = 0;
	int tPoint = 0;
	int ll = 0;
	int dd = 0;
	double aika = alku;

#ifdef MATLAB
	mexPrintf("File opened \n");
#else
	fprintf(stdout, "File opened \n");
	fflush(stdout);
#endif
	while ((i = fread(&ew1, sizeof(qb) * 3, 1, streami)) != 0 && ms <= loppu) {

		tag = ((ew1 >> (43)) & 1);

		prompt = 0;

		if (!tag) {
			if (ms >= alku) {
				prompt = (ew1 >> (42)) & 1;
				if (prompt) {
					L1 = (ew1 >> 19) & 0x1ffff;
					L2 = ew1 & 0x1ffff;
					if (L1 >= detectors || L2 >= detectors)
						continue;
					if (!storeCoordinates)
						saveSinogram(L1, L2, Sino, Ndist, Nang, ringDifference, span, sinoSize, seg, nDistSide, tPoint);
					if (saveRawData || storeCoordinates) {
						if (L2 > L1) {
							const uint32_t L3 = L1;
							L1 = L2;
							L2 = L3;
						}
						LL1[ll] = static_cast<uint16_t>(L1 + 1);
						LL2[ll] = static_cast<uint16_t>(L2 + 1);
						if (Nt > 1)
							tpoints[ll] = tPoint + 1;
						ll++;

					}
				}
				else if (randoms_correction && prompt == 0) {
					L1 = (ew1 >> 19) & 0x1ffff;
					L2 = ew1 & 0x1ffff;
					if (L1 >= detectors || L2 >= detectors)
						continue;
					if (!storeCoordinates)
						saveSinogram(L1, L2, SinoD, Ndist, Nang, ringDifference, span, sinoSize, seg, nDistSide, tPoint);
					if (storeCoordinates) {
						if (L2 > L1) {
							const uint32_t L3 = L1;
							L1 = L2;
							L2 = L3;
						}
						DD1[dd] = static_cast<uint16_t>(L1 + 1);
						DD2[dd] = static_cast<uint16_t>(L2 + 1);
						dd++;
					}
				}
			}
		}
		if (tag) {
			if (((ew1 >> (36)) & 0xff) == 160) { // Elapsed Time Tag Packet
				ms += 200e-6; // 200 microsecond increments

				if (Nt > 1 && tPoint < Nt && ms >= aika + vali[tPoint]) {
					aika += vali[tPoint];
					tPoint++;
					if (tPoint == Nt)
						break;
				}
			}
		}
	}
#ifdef MATLAB
	mexPrintf("End time %f\n", ms);
	mexEvalString("pause(.0001);");
#else
	printf("End time %f\n", ms);
#endif
	fclose(streami);
	return 1;
}
