#pragma once
#if defined(_MSC_VER)
#ifndef DLL_FUNCTION 
#define DLL_FUNCTION __declspec(dllimport)
#endif
#else
#define DLL_FUNCTION  __attribute__((visibility("default")))
#endif
#include "rootImport.h"

extern "C" DLL_FUNCTION
int rootMain(const char* rootFile, const double* tPoints, const double alku, const double loppu, bool source, const uint32_t linear_multip, const uint32_t* cryst_per_block, const uint32_t blocks_per_ring,
	const uint32_t* det_per_ring, uint16_t * S, uint16_t * SC, uint16_t * RA, uint16_t* trIndex, uint16_t* axIndex, uint16_t * DtrIndex, uint16_t * DaxIndex, bool obtain_trues, bool store_scatter, bool store_randoms, uint8_t * scatter_components,
	bool randoms_correction, float* coord, float* Dcoord, bool store_coordinates, const uint32_t* cryst_per_block_z,
	const uint32_t transaxial_multip, const uint32_t* rings, const uint64_t* sinoSize, const uint32_t Ndist, const uint32_t* Nang, const uint32_t ringDifference, const uint32_t span,
	const uint32_t * seg, const int64_t Nt, const uint64_t TOFSize, const int32_t nDistSide, uint16_t * Sino, uint16_t * SinoT, uint16_t * SinoC, uint16_t * SinoR, uint16_t * SinoD,
	const uint32_t* detWPseudo, const int32_t nPseudos, const double binSize, const double FWHM, const bool verbose, const int32_t nLayers, const float dx, const float dy, const float dz,
	const float bx, const float by, const float bz, const int64_t Nx, const int64_t Ny, const int64_t Nz, const bool dualLayerSubmodule, const bool indexBased, uint16_t * tIndex);