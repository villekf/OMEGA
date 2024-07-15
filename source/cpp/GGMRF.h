#pragma once
#include "projector_functions.h"


template <typename T>
void GGMRFKernel(T* grad, const T* u, T* weight, const int Nx, const int Ny, const int Nz, const T p, const T q, const T c, const T pqc, const T beta, const int32_t Ndx,
	const int32_t Ndy, const int32_t Ndz) {


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
		const T uj = u[n];
		int uu = 0;
		T output = (T)0.;
		for (int i = -Ndx; i <= Ndx; i++) {
			const int x_n = x + i;
			for (int j = -Ndy; j <= Ndy; j++) {
				const int y_n = y + j;
				for (int k = -Ndz; k <= Ndz; k++) {
					if (i == 0 && j == 0 && k == 0)
						continue;
					const int z_n = z + k;
					const int dim_n = z_n * Nxy + y_n * static_cast<int32_t>(Nx) + x_n;
					T uk;
					if (x_n < 0 || y_n < 0 || z_n < 0 || x_n >= Nx || y_n >= Ny || z_n >= Nz)
						uk = (T)0.;
					else
						uk = u[dim_n];
					const T delta = uj - uk;
					const T deltapqc = 1.f + std::pow(std::abs(delta / c), p - q);
					if (delta != (T)0.)
						output += weight[uu++] * (std::pow(std::abs(delta), p - (T)1.) / deltapqc) * (p - pqc * (std::pow(std::abs(delta), p - q) / deltapqc)) * sign(delta);
				}
			}
		}
		grad[n] += beta * output;
	}
}