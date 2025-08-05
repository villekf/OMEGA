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
	const uint32_t Ny, const uint32_t Nz, const int32_t Nxy, const T h, const int32_t type, const T gamma, const T epps, const T p = (T)0., const T q = (T)0., const T c = (T)0., 
	const bool anatomical = false) {

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

	int start = min_z * Nxy + 1;
	int end = max_z * Nxy;

#ifdef _OPENMP
#if _OPENMP >= 201511
#pragma omp parallel for schedule(monotonic:dynamic, nChunks)
#else
#pragma omp parallel for schedule(dynamic, nChunks)
#endif
#endif
	for (int n = start; n < end; n++) {
		const int z = n / Nxy;
		const int y = (n - z * Nxy) / static_cast<int32_t>(Nx);
		const int x = n - z * Nxy - y * static_cast<int32_t>(Nx);
		if (z < min_z || z >= max_z || x < min_x || x >= max_x || y < min_y || y >= max_y)
			continue;
		T weight_sum = (T)0.;
		T output = (T)0.;
		T outputAla = (T)0.;
		const T uj = u[n];
		for (int k = -search_window_z; k <= search_window_z; k++) {
			const int z_n = z + k;
			for (int j = -search_window_y; j <= search_window_y; j++) {
				const int y_n = y + j;
				for (int i = -search_window_x; i <= search_window_x; i++) {
					if (i == 0 && j == 0 && k == 0)
						continue;
					const int x_n = x + i;
					const int dim_n = z_n * Nxy + y_n * static_cast<int32_t>(Nx) + x_n;
					const T uk = u[dim_n];
					T distance = 0.;
					T weight = 0.;

					for (int pz = -patch_window_z; pz <= patch_window_z; pz++) {
						const int z_k = (z_n + pz) * Nxy;
						const int z_j = (z + pz) * Nxy;
						for (int py = -patch_window_y; py <= patch_window_y; py++) {
							const int y_k = (y_n + py) * static_cast<int32_t>(Nx);
							const int y_j = (y + py) * static_cast<int32_t>(Nx);
							int dim_g = (pz + patch_window_z) * (patch_window_x * 2 + 1) * (patch_window_y * 2 + 1) + (py + patch_window_y) * (patch_window_x * 2 + 1);
							for (int px = -patch_window_x; px <= patch_window_x; px++) {
								const T gg = gaussian[dim_g++];
								const int x_k = x_n + px;
								const int dim_k = z_k + y_k + x_k;
								T Pj, Pk;
								if (anatomical)
									Pj = u_ref[dim_k];
								else
									Pj = u[dim_k];
								const int x_j = x + px;
								const int dim_j = z_j + y_j + x_j;
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
					else if (type == 5) {
						// GGMRF
						const T delta = uj - uk;
						const T deltapqc = 1.f + std::pow(std::fabs(delta / c), p - q);
						output += weight * (std::pow(abs(delta), p - 1.f) / deltapqc) * (p - gamma * (std::pow(fabs(delta), p - q) / deltapqc)) * sign(delta);
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
			output /= std::sqrt(outputAla * weight_sum);

		grad[n] = output;
	}
}