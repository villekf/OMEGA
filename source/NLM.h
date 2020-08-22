#pragma once
#include <cstdint>
#include <algorithm>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

void NLM(double* grad, const double* u_ref, const double* u, const double* gaussian, const int32_t search_window_x, const int32_t search_window_y, 
	const int32_t search_window_z, const int32_t patch_window_x, const int32_t patch_window_y, const int32_t patch_window_z, const uint32_t Nx, 
	const uint32_t Ny, const uint32_t Nz, const int32_t Nxy, const double h, const int32_t type, const double epps) {

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

#pragma omp parallel for
	for (int n = start; n < end; n++) {
		const int z = n / Nxy;
		const int y = (n - z * Nxy) / static_cast<int32_t>(Nx);
		const int x = n - z * Nxy - y * static_cast<int32_t>(Nx);
		//const int dim = z * Nxy + y * static_cast<int32_t>(Nx) + x;
		if (z < min_z || z >= max_z || x < min_x || x >= max_x || y < min_y || y >= max_y)
			continue;
		double weight_sum = 0.;
		double output = 0.;
		const double uj = u[n];
		for (int k = -search_window_z; k <= search_window_z; k++) {
			const int z_n = z + k;
			for (int j = -search_window_y; j <= search_window_y; j++) {
				const int y_n = y + j;
				for (int i = -search_window_x; i <= search_window_x; i++) {
					const int x_n = x + i;
					const int dim_n = z_n * Nxy + y_n * static_cast<int32_t>(Nx) + x_n;
					const double uk = u[dim_n];
					double distance = 0.;
					double weight = 0.;

					for (int pz = -patch_window_z; pz <= patch_window_z; pz++) {
						const int z_k = (z_n + pz) * Nxy;
						const int z_j = (z + pz) * Nxy;
						for (int py = -patch_window_y; py <= patch_window_y; py++) {
							const int y_k = (y_n + py) * static_cast<int32_t>(Nx);
							const int y_j = (y + py) * static_cast<int32_t>(Nx);
							int dim_g = (pz + patch_window_z) * (patch_window_x * 2 + 1) * (patch_window_y * 2 + 1) + (py + patch_window_y) * (patch_window_x * 2 + 1);
							for (int px = -patch_window_x; px <= patch_window_x; px++) {
								const double gg = gaussian[dim_g++];
								//const double gg = 1.;
								const int x_k = x_n + px;
								const int dim_k = z_k + y_k + x_k;
								const double Pj = u_ref[dim_k];
								const int x_j = x + px;
								const int dim_j = z_j + y_j + x_j;
								const double Pk = u_ref[dim_j];
								distance += gg * (Pj - Pk) * (Pj - Pk);
							}
						}
					}
					weight = exp(-distance / h);
					weight_sum += weight;
					if (type == 2)
						output += weight * uk;
					else if (type == 0) {
						output += (weight * (uj - uk));
					}
					else {
						output += ((weight * (uj - uk)) / sqrt(weight * (uj - uk) * (uj - uk) + epps));
					}
				}
			}
		}
		//if (output == 1.)
		//	mexPrintf("n = %d\n", n);
		//if (n < start + 10000) {
		//	mexPrintf("weight_sum = %f\n", weight_sum);
		//	mexPrintf("output = %f\n", output);
		//}
		weight_sum = 1. / weight_sum;
		output *= weight_sum;

		grad[n] = output;
	}
}