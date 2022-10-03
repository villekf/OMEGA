/**************************************************************************
* Implements the Orthogonal Siddon's algorithm (Implementation 1).
* This version requires precomputation step; the number of voxels each LOR
* traverses needs to be known in advance.
* This version computes the system matrix column indices and elements for
* the preallocated MATLAB sparse matrix. Due to MATLAB's CSR format, this
* is essentially a transposed version of the system matrix.
*
* Uses OpenMP for parallelization.
*
* Copyright (C) 2020 Ville-Veikko Wettenhovi
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <https://www.gnu.org/licenses/>.
***************************************************************************/
#ifndef CT
#include "projector_functions.h"

// if 0, then determines whether the LOR intercepts the FOV
const static int TYPE = 1;

// Whether to use the OpenMP code or not
const static bool OMP = false;

// Using non-OpenMP with either precomputation or without
const static bool PRECOMPUTE = true;

const static bool DISCARD = false;

const static bool RHS = false, SUMMA = false;

const static bool no_norm = true;

using namespace std;

void orth_siddon_precomputed(const int64_t loop_var_par, const uint32_t size_x, const double zmax, size_t* indices, double* rhs, const double maxyy,
	const double maxxx, const vector<double>& xx_vec, const double dy, const vector<double>& yy_vec, const double* atten, const float* norm_coef, 
	const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const double dx, 
	const double dz, const double bx, const double by, const double bz, const bool attenuation_correction, const bool normalization, const uint16_t* lor1, 
	const uint64_t* lor2, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos, const uint16_t* L, const uint32_t* pseudos, 
	const uint32_t pRows, const uint32_t det_per_ring, const bool raw, const bool attenuation_phase, double* length, const double crystal_size, 
	const double crystal_size_z, double* y_center, double* x_center, const double* z_center, const double global_factor, const bool scatter, 
	const double* scatter_coef, const uint32_t nCores, const uint8_t list_mode) {

#ifdef _OPENMP
	if (nCores == 1U)
		setThreads();
	else
		omp_set_num_threads(nCores);
#endif

	// Precompute
	const double bzb = bz + static_cast<double>(Nz) * dz;
	const uint32_t Nyx = Ny * Nx;

#ifdef _OPENMP
#if _OPENMP >= 201511 && defined(MATLAB)
#pragma omp parallel for schedule(monotonic:dynamic, nChunks)
#else
#pragma omp parallel for schedule(dynamic, nChunks)
#endif
#endif
	for (int64_t lo = 0LL; lo < loop_var_par; lo++) {

		Det detectors;
		double kerroin, jelppi = 0., LL;

		// Load the number of voxels the LOR traverses (precomputed)
		uint32_t Np = static_cast<uint32_t>(lor1[lo]);
		if (Np == 0U)
			continue;

		// Raw list-mode data
		if (raw) {
			get_detector_coordinates_raw(det_per_ring, x, y, z_det, detectors, L, lo, pseudos, pRows, list_mode);
		}
		// Sinogram data
		else {
			get_detector_coordinates(x, y, z_det, size_x, detectors, xy_index, z_index, TotSinos, lo);
		}

		// Calculate the x, y and z distances of the detector pair
		double x_diff = (detectors.xd - detectors.xs);
		double y_diff = (detectors.yd - detectors.ys);
		const double z_diff = (detectors.zd - detectors.zs);
		// The initial index for the sparse matrix elements
		const uint64_t N12 = lor2[lo];
		const uint64_t N22 = lor2[lo + 1];

		std::vector<double> store_elements(N22 - N12, 0.);
		std::vector<uint32_t> store_indices(N22 - N12, 0U);

		size_t idx = 0ULL;
		uint8_t xyz = 0u;
		uint32_t ind = 0u;

		double ax = 0.;
		double* osem_apu = nullptr;
		double* Summ = nullptr;
		vector<double> elements;
		vector<uint32_t> v_indices;
		const double local_sino = 0.;
		const uint32_t tid = 0u;

		uint32_t N0 = Nx;
		uint32_t N1 = Ny;
		uint32_t N2 = 1u;
		uint32_t N3 = Nx;
		uint32_t N4 = Nz;

		double* xcenter = x_center;
		double* ycenter = y_center;

		// Precompute constants
		if (crystal_size_z == 0.) {
			kerroin = norm(x_diff, y_diff, z_diff) * crystal_size;
		}
		else {
			kerroin = norm(x_diff, y_diff, z_diff) * crystal_size_z;
		}
		double local_norm = 0.;
		if (normalization)
			local_norm = static_cast<double>(norm_coef[lo]);

		// If the measurement is on a same ring
		if (fabs(z_diff) < 1e-8 && (fabs(y_diff) < 1e-8 || fabs(x_diff) < 1e-8)) {

			// Z-coordinate (ring number)
			const int32_t tempk = static_cast<int32_t>(fabs(detectors.zs - bz) / dz);

			// If the LOR is perpendicular in the y-direction (Siddon cannot be used)
			if (fabs(y_diff) < 1e-8) {
				if (detectors.yd <= maxyy && detectors.yd >= by) {
					double temppi = detectors.xs;
					detectors.xs = detectors.ys;
					detectors.ys = temppi;
					double temp = 0.;
					if (crystal_size_z == 0.) {
						orth_distance_denominator_perpendicular_mfree(y_center, x_center[0], z_center, kerroin, temp, attenuation_correction, normalization, ax,
							by, detectors.yd, dy, Ny, Nx, tempk, atten, local_norm, local_sino, Ny, 1u, osem_apu, detectors, y_diff, x_diff, z_diff, store_elements, store_indices, tid,
							ind, rhs, indices, lo, PRECOMPUTE, global_factor, scatter, scatter_coef, N12);
					}
					else {
						orth_distance_denominator_perpendicular_mfree_3D(y_center, x_center[0], z_center, temp, attenuation_correction, normalization, ax,
							by, detectors.yd, dy, Ny, Nx, tempk, atten, local_norm, local_sino, Ny, 1u, osem_apu, detectors, y_diff, x_diff, z_diff, kerroin, Nyx, Nz, store_elements, store_indices, tid,
							ind, rhs, indices, lo, PRECOMPUTE, global_factor, scatter, scatter_coef, N12);
					}
				}
			}
			// Same as for the y-case above
			else if (fabs(x_diff) < 1e-8) {
				if (detectors.xd <= maxxx && detectors.xd >= bx) {
					double temp = 0.;
					if (crystal_size_z == 0.) {
						orth_distance_denominator_perpendicular_mfree(x_center, y_center[0], z_center, kerroin, temp, attenuation_correction, normalization, ax,
							bx, detectors.xd, dx, Nx, Ny, tempk, atten, local_norm, local_sino, 1u, Nx, osem_apu, detectors, x_diff, y_diff, z_diff, store_elements, store_indices, tid,
							ind, rhs, indices, lo, PRECOMPUTE, global_factor, scatter, scatter_coef, N12);
					}
					else {
						orth_distance_denominator_perpendicular_mfree_3D(x_center, y_center[0], z_center, temp, attenuation_correction, normalization, ax,
							bx, detectors.xd, dx, Nx, Ny, tempk, atten, local_norm, local_sino, 1u, Nx, osem_apu, detectors, x_diff, y_diff, z_diff, kerroin, Nyx, Nz, store_elements, store_indices, tid,
							ind, rhs, indices, lo, PRECOMPUTE, global_factor, scatter, scatter_coef, N12);
					}
				}
			}
		}
		else {
			int32_t tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;
			double txu = 0., tyu = 0., tzu = 0., tc = 0., tx0 = 1e8, ty0 = 1e8, tz0 = 1e8;
			bool skip = false;

			if (std::fabs(z_diff) < 1e-8) {
				tempk = static_cast<int32_t>(fabs(detectors.zs - bz) / dz);
				skip = siddon_pre_loop_2D(bx, by, x_diff, y_diff, maxxx, maxyy, dx, dy, Nx, Ny, tempi, tempj, txu, tyu, Np, TYPE,
					detectors.ys, detectors.xs, detectors.yd, detectors.xd, tc, iu, ju, tx0, ty0);
			}
			//Detectors on different rings (e.g. oblique sinograms)
			else if (std::fabs(y_diff) < 1e-8) {
				skip = siddon_pre_loop_2D(bx, bz, x_diff, z_diff, maxxx, bzb, dx, dz, Nx, Nz, tempi, tempk, txu, tzu, Np, TYPE,
					detectors.zs, detectors.xs, detectors.zd, detectors.xd, tc, iu, ku, tx0, tz0);
				tempj = perpendicular_start(by, detectors.yd, dy, Ny);
			}
			else if (std::fabs(x_diff) < 1e-8) {
				skip = siddon_pre_loop_2D(by, bz, y_diff, z_diff, maxyy, bzb, dy, dz, Ny, Nz, tempj, tempk, tyu, tzu, Np, TYPE,
					detectors.zs, detectors.ys, detectors.zd, detectors.yd, tc, ju, ku, ty0, tz0);
				tempi = perpendicular_start(bx, detectors.xd, dx, Nx);
				int32_t apu_tempi = tempi;
				double apu_txu = txu;
				double apu_tx0 = tx0;
				double apu_xdiff = x_diff;
				int32_t apu_iu = iu;
				const double temp_x = detectors.xs;
				detectors.xs = detectors.ys;
				detectors.ys = temp_x;
				iu = ju;
				ju = apu_iu;
				tempi = tempj;
				tempj = apu_tempi;
				txu = tyu;
				tyu = apu_txu;
				tx0 = ty0;
				ty0 = apu_tx0;
				x_diff = y_diff;
				y_diff = apu_xdiff;
				N0 = Ny;
				N1 = Nx;
				N2 = Ny;
				N3 = 1u;
				ycenter = x_center;
				xcenter = y_center;
			}
			else {
				skip = siddon_pre_loop_3D(bx, by, bz, x_diff, y_diff, z_diff, maxxx, maxyy, bzb, dx, dy, dz, Nx, Ny, Nz, tempi, tempj, tempk, tyu, txu, tzu,
					Np, TYPE, detectors, tc, iu, ju, ku, tx0, ty0, tz0);
			}
			if (attenuation_correction)
				LL = sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
			double temp = 0.;
			int alku, loppu;
			if (crystal_size_z == 0.) {
				alku = tempk + 1;
				loppu = tempk;
			}
			else {
				alku = Nz;
				loppu = 0;
				if (ku > 0) {
					alku = tempk + 1;
				}
				else if (ku < 0) {
					loppu = tempk;
				}
			}
			orth_distance_3D_full(tempi, N0, N4, y_diff, x_diff, z_diff, ycenter, xcenter, z_center, temp, N2, tempj, tempk, local_sino, ax,
				osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP, PRECOMPUTE, DISCARD, rhs, Summ, indices, elements, v_indices, idx,
				N1, N3, alku, iu, ju, loppu, store_elements, store_indices, tid, ind, N12, N22);
			// Compute the indices and matrix elements
			for (uint32_t ii = 0u; ii < Np; ii++) {
				if (tx0 < ty0 && tx0 < tz0) {
					if (attenuation_correction)
						compute_attenuation(tc, jelppi, LL, tx0, tempi, tempj, tempk, Nx, Nyx, atten);
					tempi += iu;
					tx0 += txu;
					xyz = 1U;
				}
				else if (ty0 < tz0) {
					if (attenuation_correction)
						compute_attenuation(tc, jelppi, LL, ty0, tempi, tempj, tempk, Nx, Nyx, atten);
					tempj += ju;
					ty0 += tyu;
					xyz = 2U;
				}
				else {
					if (attenuation_correction)
						compute_attenuation(tc, jelppi, LL, ty0, tempi, tempj, tempk, Nx, Nyx, atten);
					tempk += ku;
					tz0 += tzu;
					xyz = 3U;
					if (tempk < Nz && tempk >= 0) {
						alku = tempk + 1;
						loppu = tempk;
						orth_distance_3D_full(tempi, N0, N4, y_diff, x_diff, z_diff, ycenter, xcenter, z_center, temp, N2, tempj, tempk, local_sino, ax,
							osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP, PRECOMPUTE, DISCARD, rhs, Summ, indices, elements, v_indices, idx,
							N1, N3, alku, iu, ju, loppu, store_elements, store_indices, tid, ind, N12, N22);
					}
				}
			}
			if (xyz < 3 && crystal_size_z > 0. && std::fabs(z_diff) >= 1e-8) {
				if (xyz == 1)
					tempi -= iu;
				else if (xyz == 2)
					tempj -= ju;
				if ((tempk >= (Nz - 1) && ku > 0) || (tempk <= 0 && ku < 0)) {}
				else {
					tempk += ku;
					alku = Nz;
					loppu = 0;
					if (ku > 0) {
						loppu = tempk;
					}
					else if (ku < 0) {
						alku = tempk + 1;
					}
					orth_distance_3D_full(tempi, N0, N4, y_diff, x_diff, z_diff, ycenter, xcenter, z_center, temp, N2, tempj, tempk, local_sino, ax,
						osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP, PRECOMPUTE, DISCARD, rhs, Summ, indices, elements, v_indices, idx,
						N1, N3, alku, iu, ju, loppu, store_elements, store_indices, tid, ind, N12, N22);
				}
			}
			if (attenuation_phase)
				length[lo] = temp;

			temp = 1. / temp;
			if (attenuation_correction)
				temp *= exp(jelppi);
			if (normalization)
				temp *= local_norm;
			if (scatter)
				temp *= scatter_coef[lo];
			temp *= global_factor;
			for (size_t ii = 0u; ii < idx; ii++) {
				rhs[N12 + ii] *= temp;
			}
		}
	}
}

#endif