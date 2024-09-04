#pragma once
#include "projector_functions.h"


template <typename T>
void hyperbolicKernel(T* grad, const T* u, const T* w, const int Nx, const int Ny, const int Nz, const T sigma, const T beta,
	const int32_t Ndx, const int32_t Ndy, const int32_t Ndz, const uint8_t* maskBP = nullptr, const uint8_t* fovIndices = nullptr) {

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
		const int64_t z = n / Nxy;
		const int64_t y = (n - z * Nxy) / static_cast<int64_t>(Nx);
		const int64_t x = n - z * Nxy - y * static_cast<int64_t>(Nx);
		T output = (T)0.;
		const T uj = u[n];
		int uu = 0;
		for (int k = -Ndz; k <= Ndz; k++) {
			const int z_n = z + k;
			for (int j = -Ndy; j <= Ndy; j++) {
				const int y_n = y + j;
				for (int i = -Ndx; i <= Ndx; i++) {
					const int x_n = x + i;
					const int dim_n = z_n * Nxy + y_n * static_cast<int32_t>(Nx) + x_n;
					T uk;
					if (x_n < 0 || y_n < 0 || z_n < 0 || x_n >= Nx || y_n >= Ny || z_n >= Nz)
						uk = (T)0.;
					else
						uk = u[dim_n];
					const T ux = (uj - uk) / sigma;
					output += (ux / sigma) / std::sqrt((T)1. + ux * ux) * w[uu];
					uu++;
				}
			}
		}
		grad[n] += beta * output;
	}
}