#pragma once
#include "projector_functions.h"

template <typename T>
void medianFilter3D(const T* grad, T* output, const int Nx, const int Ny, const int Nz, const int NxOrig, const int NyOrig, const int NzOrig, const int32_t search_window_x,
	const int32_t search_window_y, const int32_t search_window_z, const uint8_t* maskBP = nullptr, const uint8_t* fovIndices = nullptr) {

	setThreads();
	int64_t start = 0;
	int64_t end = static_cast<int64_t>(Nx) * static_cast<int64_t>(Ny) * static_cast<int64_t>(Nz);
	int64_t Nxy = static_cast<int64_t>(Nx) * static_cast<int64_t>(Ny);
#ifdef _OPENMP
#if _OPENMP >= 201511
#pragma omp parallel for schedule(monotonic:dynamic, nChunks)
#else
#pragma omp parallel for schedule(dynamic, nChunks)
#endif
#endif
	for (int64_t n = start; n < end; n++) {
		int64_t z = n / Nxy;
		int64_t y = (n - z * Nxy) / static_cast<int64_t>(Nx);
		int64_t x = n - z * Nxy - y * static_cast<int64_t>(Nx);
		x += search_window_x;
		y += search_window_y;
		z += search_window_z;
		//if (fovIndices[xyz.z] == 0)
		//	return;
		//const int maskVal = read_imageui(maskBP, sampler_MASK, (int2)(xyz.x, xyz.y)).w;
		//if (maskVal == 0)
		//	return;
		int koko = (search_window_x * 2 + 1) * (search_window_y * 2 + 1) * (search_window_z * 2 + 1);
		std::vector<T> median(koko);
		std::vector<T> medianF(koko);
		for (int ll = 0; ll < koko; ll++) {
			medianF[ll] = (T)0.;
			median[ll] = (T)0.;
		}
		int uu = 0;
		for (int64_t xx = -search_window_x; xx <= search_window_x; xx++) {
			for (int64_t yy = -search_window_y; yy <= search_window_y; yy++) {
				for (int64_t zz = -search_window_z; zz <= search_window_z; zz++) {
					int64_t pikseli = (x + (xx)) + (y + (yy)) * static_cast<int64_t>((Nx + search_window_x * 2)) + (z + (zz)) * static_cast<int64_t>((Nx + search_window_x * 2)) * static_cast<int64_t>((Ny + search_window_y * 2));
					median[uu] = grad[pikseli];
					uu++;
				}
			}
		}
		for (int hh = 0; hh < koko; hh++) {
			int ind = 0;
			for (int ll = 0; ll < koko; ll++) {
				if (median[hh] > median[ll] || (median[hh] == median[ll] && hh < ll))
					ind++;
			}
			medianF[ind] = median[hh];
			if (ind == koko / 2)
				break;
		}
		output[n] = medianF[koko / 2];
	}
}