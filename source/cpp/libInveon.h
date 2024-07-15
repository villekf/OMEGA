#pragma once
#if defined(_MSC_VER)
#ifndef DLL_FUNCTION 
#define DLL_FUNCTION __declspec(dllimport)
#endif
#else
#define DLL_FUNCTION  __attribute__((visibility("default")))
#endif
#include "inveon.h"

extern "C" DLL_FUNCTION
int inveonMain(const double* vali, const double alku, const double loppu, const uint64_t pituus, const uint32_t detectors, const bool randoms_correction, const uint64_t sinoSize,
	const bool saveRawData, const uint32_t Ndist, const uint32_t Nang, const uint32_t ringDifference, const uint32_t span, const uint32_t * seg, const uint64_t NT,
	const int32_t nDistSide, const bool storeCoordinates, const char* argv, uint16_t * LL1, uint16_t * LL2, uint16_t * tpoints, uint16_t * DD1, uint16_t * DD2, uint16_t * Sino, uint16_t * SinoD);