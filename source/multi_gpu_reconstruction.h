#pragma once

void reconstruction(const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, const mxArray* Sin, const int Nx, const int Ny,
	const int Nz, const int Niter, const mxArray* options, const float dx, const float dy, const float dz, const float bx, const float by, const float bz,
	const float bzb, const float maxxx, const float maxyy, const float zmax, const float NSlices, const uint32_t* pituus, const size_t koko_l,
	const uint32_t* xy_index, const uint16_t* z_index, const int size_x, const int TotSinos, mxArray* cell, const mwSize* dimmi, const bool verbose,
	const int attenuation_correction, const float* atten, const size_t size_atten, const int subsets, const float epps, const uint8_t* rekot,
	const char* k_path, const size_t size_rekot, const int Nt, const int* pseudos, const int det_per_ring, const int prows, const uint16_t* L,
	const uint8_t raw, const size_t size_z, const bool osem_bool, const char* fileName, const bool force_build, const int device, float kerroin);