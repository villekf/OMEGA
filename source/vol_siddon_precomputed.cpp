/**************************************************************************
* Implements the volume-based Siddon's algorithm (Implementation 1).
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
#include "projector_functions.h"

// if 0, then determines whether the LOR intercepts the FOV
const static int TYPE = 1;

// Whether to use the OpenMP code or not
const static bool OMP = false;

// Using non-OpenMP with either precomputation or without
const static bool PRECOMPUTE = true;

const static bool DISCARD = false;

using namespace std;

void vol_siddon_precomputed(const int64_t loop_var_par, const uint32_t size_x, const double zmax, mwIndex* indices, double* rhs, const double maxyy,
	const double maxxx, const vector<double>& xx_vec, const double dy, const vector<double>& yy_vec, const double* atten, const double* norm_coef,
	const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const double dx,
	const double dz, const double bx, const double by, const double bz, const bool attenuation_correction, const bool normalization, const uint16_t* lor1,
	const uint64_t* lor2, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos, const uint16_t* L, const uint32_t* pseudos,
	const uint32_t pRows, const uint32_t det_per_ring, const bool raw, const bool attenuation_phase, double* length, const double crystal_size,
	const double crystal_size_z, const double* y_center, const double* x_center, const double* z_center, const double global_factor, const double bmin, 
	const double bmax, const double Vmax, const double* V) {

	setThreads();

	// Precompute
	const double bzb = bz + static_cast<double>(Nz) * dz;
	const uint32_t Nyx = Ny * Nx;

	size_t threads = omp_get_max_threads();
	std::vector<double> store_elements;
	std::vector<uint32_t> store_indices;


#pragma omp parallel for schedule(dynamic)
	for (int64_t lo = 0LL; lo < loop_var_par; lo++) {

		Det detectors;
		double kerroin, jelppi = 0., LL;

		// Raw list-mode data
		if (raw) {
			get_detector_coordinates_raw(det_per_ring, x, y, z_det, detectors, L, lo, pseudos, pRows);
		}
		// Sinogram data
		else {
			get_detector_coordinates(x, y, z_det, size_x, detectors, xy_index, z_index, TotSinos, lo);
		}

		// Calculate the x, y and z distances of the detector pair
		const double x_diff = (detectors.xd - detectors.xs);
		const double y_diff = (detectors.yd - detectors.ys);
		const double z_diff = (detectors.zd - detectors.zs);

		// Load the number of voxels the LOR traverses (precomputed)
		uint32_t Np = static_cast<uint32_t>(lor1[lo]);
		// The initial index for the sparse matrix elements
		const uint64_t N2 = lor2[lo];
		const uint64_t N22 = lor2[lo + 1];

		size_t idx = 0ULL;
		int8_t start = 1;
		uint8_t xyz = 0u;
		uint8_t xyz_w = 0u;
		uint32_t ind = 0u;

		double ax = 0.;
		double* osem_apu = nullptr;
		double* Summ = nullptr;
		const bool no_norm = true;
		vector<double> elements;
		vector<uint32_t> v_indices;
		const double local_sino = 0.;
		const bool RHS = false, SUMMA = false;
		const uint32_t tid = 0u;

		// Precompute constants
		kerroin = norm(x_diff, y_diff, z_diff);

		// If the measurement is on a same ring
		if (fabs(z_diff) < 1e-8) {

			// Z-coordinate (ring number)
			const uint32_t tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));

			// If the LOR is perpendicular in the y-direction (Siddon cannot be used)
			if (fabs(y_diff) < 1e-8) {
				if (detectors.yd <= maxyy && detectors.yd >= by) {
					double temppi = detectors.xs;
					detectors.xs = detectors.ys;
					detectors.ys = temppi;
					double temp = 0.;
					volume_distance_denominator_perpendicular_mfree_3D(y_center, x_center[0], z_center, temp, attenuation_correction, normalization, ax,
						by, detectors.yd, dy, Ny, Nx, tempk, atten, norm_coef, local_sino, Ny, 1u, osem_apu, detectors, y_diff, x_diff, z_diff, kerroin, Nyx, Nz, store_elements, store_indices, tid,
						ind, rhs, indices, lo, PRECOMPUTE, global_factor, bmax, bmin, Vmax, V, N2);
				}
			}
			// Same as for the y-case above
			else if (fabs(x_diff) < 1e-8) {
				if (detectors.xd <= maxxx && detectors.xd >= bx) {
					double temp = 0.;
					volume_distance_denominator_perpendicular_mfree_3D(x_center, y_center[0], z_center, temp, attenuation_correction, normalization, ax,
						bx, detectors.xd, dx, Nx, Ny, tempk, atten, norm_coef, local_sino, 1u, Nx, osem_apu, detectors, x_diff, y_diff, z_diff, kerroin, Nyx, Nz, store_elements, store_indices, tid,
						ind, rhs, indices, lo, PRECOMPUTE, global_factor, bmax, bmin, Vmax, V, N2);
				}
			}
			else {
				// If neither x- nor y-directions are perpendicular
				// Correspond to the equations (9) and (10) from reference [2]
				int32_t tempi = 0, tempj = 0, iu = 0, ju = 0;
				double txu = 0., tyu = 0., tc = 0., tx0 = 0., ty0 = 0.;
				const bool skip = siddon_pre_loop_2D(bx, by, x_diff, y_diff, maxxx, maxyy, dx, dy, Nx, Ny, tempi, tempj, txu, tyu, Np, TYPE,
					detectors.ys, detectors.xs, detectors.yd, detectors.xd, tc, iu, ju, tx0, ty0);
				if (attenuation_correction)
					LL = sqrt(x_diff * x_diff + y_diff * y_diff);
				double temp = 0.;
				int alku, loppu;
				alku = Nz;
				loppu = 0;
				volume_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
					tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
					PRECOMPUTE, DISCARD, rhs, Summ, indices, elements, v_indices, idx, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind, 
					bmax, bmin, Vmax, V, N2, N22);
				if (attenuation_correction) {
					// Compute the indices and matrix elements
					for (uint32_t ii = 0u; ii < Np; ii++) {
						// Ray goes along the x-axis
						if (tx0 < ty0) {
							if (attenuation_correction)
								compute_attenuation(tc, jelppi, LL, tx0, tempi, tempj, tempk, Nx, Nyx, atten);
							tempi += iu;
							tx0 += txu;
							xyz = 1u;
						}
						// Ray goes along the y-axis
						else {
							if (attenuation_correction)
								compute_attenuation(tc, jelppi, LL, ty0, tempi, tempj, tempk, Nx, Nyx, atten);
							tempj += ju;
							ty0 += tyu;
							xyz = 2u;
						}
					}
				}
				if (attenuation_phase)
					length[lo] = temp;

				temp = 1. / temp;
				if (attenuation_correction)
					temp *= exp(jelppi);
				if (normalization)
					temp *= norm_coef[lo];
				temp *= global_factor;
				for (size_t ii = 0u; ii < idx; ii++) {
					rhs[N2 + ii] *= temp;
				}
			}
		}
		//If the z-detector coordinates are not on the same ring
		//All computations follow the same logic as above
		else {
			if (fabs(y_diff) < 1e-8) {
				if (detectors.yd <= maxyy && detectors.yd >= by) {
					int32_t tempi = 0, tempk = 0, tempj = 0, iu = 0, ku = 0;
					double txu = 0., tzu = 0., tc = 0., tx0 = 0., tz0 = 0.;
					const bool skip = siddon_pre_loop_2D(bx, bz, x_diff, z_diff, maxxx, bzb, dx, dz, Nx, Nz, tempi, tempk, txu, tzu, Np, TYPE,
						detectors.zs, detectors.xs, detectors.zd, detectors.xd, tc, iu, ku, tx0, tz0);
					double apu1;
					if (attenuation_correction)
						LL = sqrt(x_diff * x_diff + z_diff * z_diff);
					for (size_t ii = 0ULL; ii < static_cast<size_t>(Ny); ii++) {
						apu1 = (yy_vec[ii + 1ULL] - detectors.yd);
						if (apu1 > 0.) {
							tempj = static_cast<int32_t>(ii);
							break;
						}
					}
					double temp = 0.;
					int alku, loppu;
					alku = Nz;
					loppu = 0;
					if (ku > 0) {
						alku = tempk + 1;
					}
					else if (ku < 0) {
						loppu = tempk;
					}
					volume_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
						tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
						PRECOMPUTE, DISCARD, rhs, Summ, indices, elements, v_indices, idx, Ny, Nx, alku, iu, 0, loppu, store_elements, store_indices, tid, ind,
						bmax, bmin, Vmax, V, N2, N22);
					for (uint32_t ii = 0u; ii < Np; ii++) {
						if (tx0 < tz0) {
							if (attenuation_correction)
								compute_attenuation(tc, jelppi, LL, tx0, tempi, tempj, tempk, Nx, Nyx, atten);
							tempi += iu;
							tx0 += txu;
							xyz = 1u;
						}
						else {
							if (attenuation_correction)
								compute_attenuation(tc, jelppi, LL, tz0, tempi, tempj, tempk, Nx, Nyx, atten);
							tempk += ku;
							tz0 += tzu;
							xyz = 3u;
							if (tempk < Nz && tempk >= 0) {
								alku = tempk + 1;
								loppu = tempk;
								volume_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
									tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
									PRECOMPUTE, DISCARD, rhs, Summ, indices, elements, v_indices, idx, Ny, Nx, alku, iu, 0, loppu, store_elements, store_indices, tid, ind,
									bmax, bmin, Vmax, V, N2, N22);
							}
						}
					}
					if (xyz < 3) {
						tempi -= iu;
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
							volume_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
								tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
								PRECOMPUTE, DISCARD, rhs, Summ, indices, elements, v_indices, idx, Ny, Nx, alku, iu, 0, loppu, store_elements, store_indices, tid, ind,
								bmax, bmin, Vmax, V, N2, N22);
						}
					}
					if (attenuation_phase)
						length[lo] = temp;

					temp = 1. / temp;
					if (attenuation_correction)
						temp *= exp(jelppi);
					if (normalization)
						temp *= norm_coef[lo];
					temp *= global_factor;
					for (size_t ii = 0u; ii < idx; ii++) {
						rhs[N2 + ii] *= temp;
					}
				}
			}
			else if (fabs(x_diff) < 1e-8) {
				if (detectors.xd <= maxxx && detectors.xd >= bx) {
					int32_t tempi = 0, tempk = 0, tempj = 0, ju = 0, ku = 0;
					double tyu = 0., tzu = 0., tc = 0., ty0 = 0., tz0 = 0.;
					const bool skip = siddon_pre_loop_2D(by, bz, y_diff, z_diff, maxyy, bzb, dy, dz, Ny, Nz, tempj, tempk, tyu, tzu, Np, TYPE,
						detectors.zs, detectors.ys, detectors.zd, detectors.yd, tc, ju, ku, ty0, tz0);
					if (attenuation_correction)
						LL = sqrt(z_diff * z_diff + y_diff * y_diff);
					double apu1;
					double temp = 0.;
					for (size_t ii = 0ULL; ii < static_cast<size_t>(Nx); ii++) {
						apu1 = (xx_vec[ii + 1ULL] - detectors.xd);
						if (apu1 > 0.) {
							tempi = static_cast<int32_t>(ii);
							break;
						}
					}
					const double temp_x = detectors.xs;
					detectors.xs = detectors.ys;
					detectors.ys = temp_x;
					int alku, loppu;
					alku = Nz;
					loppu = 0;
					if (ku > 0) {
						alku = tempk + 1;
					}
					else if (ku < 0) {
						loppu = tempk;
					}
					volume_distance_3D_full(tempj, Ny, Nz, x_diff, y_diff, z_diff, x_center, y_center, z_center, temp, Nx,
						tempi, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
						PRECOMPUTE, DISCARD, rhs, Summ, indices, elements, v_indices, idx, Nx, 1u, alku, ju, 0, loppu, store_elements, store_indices, tid, ind,
						bmax, bmin, Vmax, V, N2, N22);
					for (uint32_t ii = 0u; ii < Np; ii++) {
						if (ty0 < tz0) {
							if (attenuation_correction)
								compute_attenuation(tc, jelppi, LL, ty0, tempi, tempj, tempk, Nx, Nyx, atten);
							tempj += ju;
							ty0 += tyu;
							xyz = 2u;
						}
						else {
							if (attenuation_correction)
								compute_attenuation(tc, jelppi, LL, tz0, tempi, tempj, tempk, Nx, Nyx, atten);
							tempk += ku;
							tz0 += tzu;
							xyz = 3u;

							if (tempk < Nz && tempk >= 0) {
								alku = tempk + 1;
								loppu = tempk;
								volume_distance_3D_full(tempj, Ny, Nz, x_diff, y_diff, z_diff, x_center, y_center, z_center, temp, Nx,
									tempi, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
									PRECOMPUTE, DISCARD, rhs, Summ, indices, elements, v_indices, idx, Nx, 1u, alku, ju, 0, loppu, store_elements, store_indices, tid, ind,
									bmax, bmin, Vmax, V, N2, N22);
							}
						}
					}
					if (xyz < 3) {
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
							volume_distance_3D_full(tempj, Ny, Nz, x_diff, y_diff, z_diff, x_center, y_center, z_center, temp, Nx,
								tempi, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
								PRECOMPUTE, DISCARD, rhs, Summ, indices, elements, v_indices, idx, Nx, 1u, alku, ju, 0, loppu, store_elements, store_indices, tid, ind,
								bmax, bmin, Vmax, V, N2, N22);
						}
					}
					if (attenuation_phase)
						length[lo] = temp;

					temp = 1. / temp;
					if (attenuation_correction)
						temp *= exp(jelppi);
					if (normalization)
						temp *= norm_coef[lo];
					temp *= global_factor;
					for (size_t ii = 0u; ii < idx; ii++) {
						rhs[N2 + ii] *= temp;
					}
				}
			}
			else {
				int32_t tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 1;
				double txu = 0., tyu = 0., tzu = 0., tc = 0., tx0 = 0., ty0 = 0., tz0 = 0.;
				const bool skip = siddon_pre_loop_3D(bx, by, bz, x_diff, y_diff, z_diff, maxxx, maxyy, bzb, dx, dy, dz, Nx, Ny, Nz, tempi, tempj, tempk,
					tyu, txu, tzu, Np, TYPE, detectors, tc, iu, ju, ku, tx0, ty0, tz0);

				double temp = 0.;

				if (attenuation_correction)
					LL = sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
				int alku, loppu;
				alku = Nz;
				loppu = 0;
				if (ku > 0) {
					alku = tempk + 1;
				}
				else if (ku < 0) {
					loppu = tempk;
				}
				volume_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
					tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
					PRECOMPUTE, DISCARD, rhs, Summ, indices, elements, v_indices, idx, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind,
					bmax, bmin, Vmax, V, N2, N22);

				for (uint32_t ii = 0u; ii < Np; ii++) {
					if (tz0 < ty0 && tz0 < tx0) {
						if (attenuation_correction)
							compute_attenuation(tc, jelppi, LL, tz0, tempi, tempj, tempk, Nx, Nyx, atten);
						tempk += ku;
						tz0 += tzu;
						xyz = 3u;
						if (tempk < Nz && tempk >= 0) {
							alku = tempk + 1;
							loppu = tempk;
							volume_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
								tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
								PRECOMPUTE, DISCARD, rhs, Summ, indices, elements, v_indices, idx, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind,
								bmax, bmin, Vmax, V, N2, N22);
						}
					}
					else if (ty0 < tx0) {
						if (attenuation_correction)
							compute_attenuation(tc, jelppi, LL, ty0, tempi, tempj, tempk, Nx, Nyx, atten);
						tempj += ju;
						ty0 += tyu;
						xyz = 2u;
					}
					else {
						if (attenuation_correction)
							compute_attenuation(tc, jelppi, LL, tx0, tempi, tempj, tempk, Nx, Nyx, atten);
						tempi += iu;
						tx0 += txu;
						xyz = 1u;
					}
				}
				if (xyz < 3) {
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
						volume_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
							tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
							PRECOMPUTE, DISCARD, rhs, Summ, indices, elements, v_indices, idx, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind,
							bmax, bmin, Vmax, V, N2, N22);
					}
				}
				if (attenuation_phase)
					length[lo] = temp;

				temp = 1. / temp;
				if (attenuation_correction)
					temp *= exp(jelppi);
				if (normalization)
					temp *= norm_coef[lo];
				temp *= global_factor;
				for (size_t ii = 0u; ii < idx; ii++) {
					rhs[N2 + ii] *= temp;
				}
			}
		}
	}
}
