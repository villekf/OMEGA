#pragma once
#include "projector_functions.h"


// Neighbor read, each coordinate is clamped to the volume (clamp-to-edge, i.e. Neumann boundary), same as the GPU kernels
template <typename T>
inline T ggmrfRead(const T* buf, const int64_t x, const int64_t y, const int64_t z, const int64_t Nx, const int64_t Ny, const int64_t Nz) {
	const int64_t xc = std::max<int64_t>(0, std::min<int64_t>(x, Nx - 1));
	const int64_t yc = std::max<int64_t>(0, std::min<int64_t>(y, Ny - 1));
	const int64_t zc = std::max<int64_t>(0, std::min<int64_t>(z, Nz - 1));
	return buf[xc + yc * Nx + zc * Nx * Ny];
}

template <typename T>
void GGMRFKernel(T* grad, const T* u, T* weight, const int Nx, const int Ny, const int Nz, const T p, const T q, const T c, const T pqc, const T beta, const int32_t Ndx,
	const int32_t Ndy, const int32_t Ndz, const uint8_t* maskBP = nullptr, const uint8_t* fovIndices = nullptr, const uint32_t maskBPZ = 1U) {


	setThreads();
	int64_t start = 0;
	int64_t end = static_cast<int64_t>(Nx) * static_cast<int64_t>(Ny) * static_cast<int64_t>(Nz);
	int64_t Nxy = static_cast<int64_t>(Nx) * static_cast<int64_t>(Ny);
	const int64_t Nx64 = static_cast<int64_t>(Nx);
	const int64_t Ny64 = static_cast<int64_t>(Ny);
	const int64_t Nz64 = static_cast<int64_t>(Nz);
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
		// Voxels outside the extended FOV or the prior mask do not contribute, same as the GPU kernel
		if (fovIndices != nullptr && fovIndices[z] == 0)
			continue;
		if (maskBP != nullptr) {
			const int64_t maskInd = maskBPZ > 1U ? n : x + y * Nx64;
			if (maskBP[maskInd] == 0)
				continue;
		}
		const T uj = u[n];
		int uu = 0;
		T output = (T)0.;
		for (int i = -Ndx; i <= Ndx; i++) {
			const int64_t x_n = x + i;
			for (int j = -Ndy; j <= Ndy; j++) {
				const int64_t y_n = y + j;
				for (int k = -Ndz; k <= Ndz; k++) {
					if (i == 0 && j == 0 && k == 0)
						continue;
					const int64_t z_n = z + k;
					// Neighbors outside the volume are clamped to the edge (Neumann boundary), identical to the GPU kernel
					const T uk = ggmrfRead(u, x_n, y_n, z_n, Nx64, Ny64, Nz64);
					const T delta = uj - uk;
					if (delta != (T)0.) {
						const T deltapqc = (T)1. + std::pow(std::abs(delta / c), p - q);
						output += weight[uu] * (std::pow(std::abs(delta), p - (T)1.) / deltapqc) * (p - pqc * (std::pow(std::abs(delta), p - q) / deltapqc)) * sign(delta);
					}
					uu++;
				}
			}
		}
		grad[n] += beta * output;
	}
}