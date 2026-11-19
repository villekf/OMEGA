#pragma once
#include <cstdint>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include <thread>
#ifdef _OPENMP
#include <omp.h>
#ifndef nChunks
#define nChunks 100
#endif
#endif
#include "projector_functions.h"

// Non-local means
template <typename T>
void NLMFunc(T* grad, const T* u_ref, const T* u, const T* gaussian, const int32_t search_window_x, const int32_t search_window_y,
	const int32_t search_window_z, const int32_t patch_window_x, const int32_t patch_window_y, const int32_t patch_window_z, const uint32_t Nx, 
	const uint32_t Ny, const uint32_t Nz, const int64_t Nxy, const T h, const int32_t type, const T gamma, const T epps, const T p = (T)0., const T q = (T)0., const T c = (T)0., 
	const bool anatomical = false, const uint8_t* maskBP = nullptr, const uint8_t* fovIndices = nullptr, const uint32_t maskBPZ = 1U) {

	setThreads();

	const int window_x = search_window_x + patch_window_x;
	const int window_y = search_window_y + patch_window_y;
	const int window_z = search_window_z + patch_window_z;

	const int min_x = window_x;
	const int max_x = static_cast<int32_t>(Nx) - window_x;
	const int min_y = window_y;
	const int max_y = static_cast<int32_t>(Ny) - window_y;
	const int min_z = window_z;
	const int max_z = static_cast<int32_t>(Nz) - window_z;

	// The volume is padded by the search window + patch window in each direction, the mask and the extended FOV indices are in the original volume
	const int64_t NxOrig = static_cast<int64_t>(Nx) - 2 * window_x;
	const int64_t NyOrig = static_cast<int64_t>(Ny) - 2 * window_y;
	const int64_t Nx64 = static_cast<int64_t>(Nx);

	const int64_t start = static_cast<int64_t>(min_z) * Nxy;
	const int64_t end = static_cast<int64_t>(max_z) * Nxy;

#ifdef _OPENMP
#if _OPENMP >= 201511
#pragma omp parallel for schedule(monotonic:dynamic, nChunks)
#else
#pragma omp parallel for schedule(dynamic, nChunks)
#endif
#endif
	for (int64_t n = start; n < end; n++) {
		const int z = static_cast<int>(n / Nxy);
		const int y = static_cast<int>((n - z * Nxy) / Nx64);
		const int x = static_cast<int>(n - z * Nxy - y * Nx64);
		if (z < min_z || z >= max_z || x < min_x || x >= max_x || y < min_y || y >= max_y)
			continue;
		// Voxels outside the extended FOV or the prior mask do not contribute (the output is left untouched), same as the GPU kernel
		if (fovIndices != nullptr && fovIndices[z - window_z] == 0)
			continue;
		if (maskBP != nullptr) {
			const int64_t xOrig = x - window_x;
			const int64_t yOrig = y - window_y;
			const int64_t zOrig = z - window_z;
			const int64_t maskInd = maskBPZ > 1U ? xOrig + yOrig * NxOrig + zOrig * NxOrig * NyOrig : xOrig + yOrig * NxOrig;
			if (maskBP[maskInd] == 0)
				continue;
		}
		T weight_sum = (T)0.;
		T output = (T)0.;
		// NLTV starts from epps, same as the GPU kernel
		T outputAla = epps;
		const T uj = u[n];
		for (int k = -search_window_z; k <= search_window_z; k++) {
			const int z_n = z + k;
			for (int j = -search_window_y; j <= search_window_y; j++) {
				const int y_n = y + j;
				for (int i = -search_window_x; i <= search_window_x; i++) {
					if (i == 0 && j == 0 && k == 0)
						continue;
					const int x_n = x + i;
					const int64_t dim_n = z_n * Nxy + y_n * Nx64 + x_n;
					const T uk = u[dim_n];
					T distance = 0.;
					T weight = 0.;

					for (int pz = -patch_window_z; pz <= patch_window_z; pz++) {
						const int64_t z_k = (z_n + pz) * Nxy;
						const int64_t z_j = (z + pz) * Nxy;
						for (int py = -patch_window_y; py <= patch_window_y; py++) {
							const int64_t y_k = (y_n + py) * Nx64;
							const int64_t y_j = (y + py) * Nx64;
							int dim_g = (pz + patch_window_z) * (patch_window_x * 2 + 1) * (patch_window_y * 2 + 1) + (py + patch_window_y) * (patch_window_x * 2 + 1);
							for (int px = -patch_window_x; px <= patch_window_x; px++) {
								const T gg = gaussian[dim_g++];
								const int x_k = x_n + px;
								const int64_t dim_k = z_k + y_k + x_k;
								T Pj, Pk;
								if (anatomical)
									Pj = u_ref[dim_k];
								else
									Pj = u[dim_k];
								const int x_j = x + px;
								const int64_t dim_j = z_j + y_j + x_j;
								if (anatomical)
									Pk = u_ref[dim_j];
								else
									Pk = u[dim_j];
								distance += gg * (Pj - Pk) * (Pj - Pk);
							}
						}
					}
					weight = std::exp(-distance / h);
					weight_sum += weight;
					if (type == 2 || type == 5)
						output += weight * uk;
					else if (type == 0) {
						output += (weight * (uj - uk));
					}
					else if (type == 3) {
						// RD
						const T u = (uj - uk);
						const T divPow = (uj + uk + gamma * std::fabs(u) + epps);
						output += weight * u * (gamma * std::fabs(u) + uj + (T)3. * uk + epps * epps) / (divPow * divPow);
					}
					else if (type == 4) {
						// Lange
						const T u = (uj - uk);
						const T uabs = sign(u);
						output += weight * (uabs - uabs / (std::fabs(u) / gamma + (T)1.));
					}
					else if (type == 7) {
						// Geman-McClure, i.e. the derivative of the potential delta^2 / (delta^2 + gamma^2)
						const T delta = uj - uk;
						const T apu = delta * delta + gamma * gamma;
						output += weight * ((T)2. * gamma * gamma * delta) / (apu * apu);
					}
					else if (type == 6) {
						// GGMRF
						// Identical voxels are skipped, same as the GPU kernel
						const T delta = uj - uk;
						if (delta != (T)0.) {
							const T deltapqc = (T)1. + std::pow(std::fabs(delta / c), p - q);
							output += weight * (std::pow(std::fabs(delta), p - (T)1.) / deltapqc) * (p - gamma * (std::pow(std::fabs(delta), p - q) / deltapqc)) * sign(delta);
						}
					}
					else {
						const T apuU = uj - uk;
						output += (weight * apuU);
						outputAla += weight * apuU * apuU;
					}
				}
			}
		}
		if (weight_sum != (T)0.)
			weight_sum = (T)1. / (weight_sum + epps);
		output *= weight_sum;
		if (type == 2)
			output = uj - output;
		else if (type == 5) {
			// Lange with NLMRP
			output = uj - output;
			const T uabs = sign(output);
			output = (uabs - uabs / (std::fabs(output) / gamma + (T)1.));
		}
		else if (type == 1)
			output /= std::sqrt(outputAla * weight_sum + epps);

		grad[n] = output;
	}
}