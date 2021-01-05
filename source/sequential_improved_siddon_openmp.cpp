/**************************************************************************
* Implements both the improved Siddon's algorithm and Orthogonal Siddon's
* algorithm for OMEGA (Implementation 4).
* Uses the precomputed data (faster).
*
* Uses OpenMP for parallelization. If OpenMP is not available, the code
* is serial with no parallelization.
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

// if 0, then determines whether the LOR intercepts the FOV (i.e. no precomputation phase performed)
const static int TYPE = 1;

// Whether to use the OpenMP code or not
const static bool OMP = true;

// Using non-OpenMP with either precomputation or without
const static bool PRECOMPUTE = false;

const static bool DISCARD = false;

using namespace std;

void sequential_improved_siddon(const int64_t loop_var_par, const uint32_t size_x, const double zmax, double* Summ, double* rhs, const double maxyy,
	const double maxxx,	const vector<double>& xx_vec, const double dy, const vector<double>& yy_vec, const double* atten, const double* norm_coef, 
	const double* randoms, const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, 
	const uint32_t Nz, const double dx, const double dz, const double bx, const double by, const double bz, const bool attenuation_correction, 
	const bool normalization,  const bool randoms_correction, const uint16_t* lor1, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos, 
	const double epps,  const double* Sino, double* osem_apu, const uint16_t* L, const uint32_t* pseudos, const uint32_t pRows, const uint32_t det_per_ring, 
	const bool raw, const bool no_norm, const double global_factor, const uint8_t fp, const bool scatter, const double* scatter_coef, const bool TOF, 
	const int64_t TOFSize, const double sigma_x, const double* TOFCenter, const int64_t nBins, const uint32_t dec_v, const uint32_t nCores) {

	if (nCores == 1U)
		setThreads();
	else
		omp_set_num_threads(nCores);

	const uint32_t Nyx = Ny * Nx;

	const double bzb = bz + static_cast<double>(Nz) * dz;

#ifdef _OPENMP
	size_t threads = omp_get_max_threads();
#else
	size_t threads = 1ULL;
#endif
	vector<double> TOFVal(nBins * dec_v * threads, 0.);

	//mexPrintf("fp = %u\n", fp);

#pragma omp parallel for ordered schedule(dynamic)
	for (int64_t lo = 0LL; lo < loop_var_par; lo++) {

		double local_sino = 0.;
		if (TOF) {
			for (int64_t to = 0LL; to < nBins; to++)
				local_sino += Sino[lo + TOFSize * to];
		}
		else {
			local_sino = Sino[lo];
		}
		if (no_norm && local_sino == 0.)
			continue;
		Det detectors;

		// Raw data
		if (raw) {
			get_detector_coordinates_raw(det_per_ring, x, y, z_det, detectors, L, lo, pseudos, pRows);
		}
		// Sinogram data
		else {
			get_detector_coordinates(x, y, z_det, size_x, detectors, xy_index, z_index, TotSinos, lo);
		}

		// Calculate the x, y and z distances of the detector pair
		const double y_diff = (detectors.yd - detectors.ys);
		const double x_diff = (detectors.xd - detectors.xs);
		const double z_diff = (detectors.zd - detectors.zs);

		// Load the number of voxels the LOR traverses (precomputed)
		uint32_t Np = static_cast<uint32_t>(lor1[lo]);
		double jelppi = 0.;
		double D = 0., DD = 0.;
		double xI = 0., yI = 0., zI = 0.;

		vector<double> ax(nBins, 0.);
		vector<double> yax(nBins, 0.);


		if (fp == 2) {
			for (int64_t to = 0LL; to < nBins; to++)
				yax[to] = osem_apu[lo + to * loop_var_par];
		}

#ifdef _OPENMP
		const int64_t tid = omp_get_thread_num() * dec_v * nBins;
#else
		const int64_t tid = 1LL;
#endif

		if (fabs(z_diff) < 1e-8 && (fabs(y_diff) < 1e-8 || fabs(x_diff) < 1e-8)) {

			const uint32_t tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));

			if (fabs(y_diff) < 1e-8) {

				if (detectors.yd <= maxyy && detectors.yd >= by) {
					uint32_t temp_ijk = 0;

					const double element = perpendicular_elements(Ny, detectors.yd, yy_vec, dx, tempk, Nx, Ny, atten, norm_coef, attenuation_correction,
						normalization, temp_ijk, 1u, lo, global_factor, scatter, scatter_coef);

					if (TOF) {
						xI = (dx * Nx) / 2.;
						if (x_diff > 0.)
							xI = -xI;
						D = xI;
						for (uint32_t k = 0; k < Np; k++) {
							TOFWeightsFP(D, nBins, dx, TOFVal, TOFCenter, sigma_x, xI, osem_apu, temp_ijk + k, ax, epps, tid + static_cast<int64_t>(k) * nBins);
						}
						double temp = element / dx;
						double val_rhs = 0.;
						double val = 0.;
						if (fp == 1) {
							for (int64_t to = 0LL; to < nBins; to++) {
								if (ax[to] < epps)
									ax[to] = epps;
								else
									ax[to] *= temp;
								if (randoms_correction)
									ax[to] += randoms[lo];
								rhs[lo + to * loop_var_par] = ax[to];
							}
							continue;
						}
						if (local_sino > 0.) {
							if (fp != 2) {
								for (int64_t to = 0LL; to < nBins; to++) {
									if (ax[to] < epps)
										ax[to] = epps;
									else
										ax[to] *= temp;
									if (randoms_correction)
										ax[to] += randoms[lo];
									yax[to] = Sino[lo + to * TOFSize] / ax[to];
								}
							}
							for (uint32_t k = 0; k < Np; k++) {
								if (k == 0) {
									xI = D;
								}
								val_rhs = TOFWeightsBP(D, nBins, dx, TOFVal, TOFCenter, sigma_x, xI, yax, epps, temp, val, rhs, tid + static_cast<int64_t>(k) * nBins);
#pragma omp atomic
								rhs[temp_ijk + k] += (val_rhs);
								if (no_norm == 0 && val > 0.) {
#pragma omp atomic
									Summ[temp_ijk + k] += val;
								}
							}
						}
						else {
							for (uint32_t k = 0; k < Np; k++) {
								if (k == 0) {
									xI = D;
								}
								TOFWeightsSumm(D, nBins, dx, TOFVal, TOFCenter, sigma_x, xI, epps, temp, val, tid + static_cast<int64_t>(k) * nBins);
								if (no_norm == 0 && val > 0.) {
#pragma omp atomic
									Summ[temp_ijk + k] += val;
								}
							}
						}
					}
					else {
						if (fp == 1) {
							for (uint32_t k = 0; k < Np; k++) {
								ax[0] += (element * osem_apu[temp_ijk + k]);
							}
							if (ax[0] < epps)
								ax[0] = epps;
							if (randoms_correction)
								ax[0] += randoms[lo];
							rhs[lo] = ax[0];
							continue;
						}
						if (local_sino > 0.) {
							if (fp != 2) {
								for (uint32_t k = 0; k < Np; k++) {
									ax[0] += (element * osem_apu[temp_ijk + k]);
								}
								if (ax[0] < epps)
									ax[0] = epps;
								if (randoms_correction)
									ax[0] += randoms[lo];
								yax[0] = local_sino / ax[0];
							}
							for (uint32_t k = 0; k < Np; k++) {
#pragma omp atomic
								rhs[temp_ijk + k] += (element * yax[0]);
								if (no_norm == 0) {
#pragma omp atomic
									Summ[temp_ijk + k] += element;
								}
							}
						}
						else {
							for (uint32_t k = 0; k < Np; k++) {
								if (no_norm == 0) {
#pragma omp atomic
									Summ[temp_ijk + k] += element;
								}
							}
						}
					}
				}
			}
			else if (fabs(x_diff) < 1e-8) {

				if (detectors.xd <= maxxx && detectors.xd >= bx) {
					uint32_t temp_ijk = 0;

					const double element = perpendicular_elements(1, detectors.xd, xx_vec, dy, tempk, Ny, Nx, atten, norm_coef, attenuation_correction,
						normalization, temp_ijk, Nx, lo, global_factor, scatter, scatter_coef);

					if (TOF) {
						yI = (dy * Ny) / 2.;
						if (y_diff > 0.)
							yI = -yI;
						D = yI;
						for (uint32_t k = 0; k < Np; k++) {
							TOFWeightsFP(D, nBins, dy, TOFVal, TOFCenter, sigma_x, yI, osem_apu, temp_ijk + k * Nx, ax, epps, tid + static_cast<int64_t>(k) * nBins);
						}
						double temp = element / dy;
						double val_rhs = 0.;
						double val = 0.;
						if (fp == 1) {
							for (int64_t to = 0LL; to < nBins; to++) {
								if (ax[to] < epps)
									ax[to] = epps;
								else
									ax[to] *= temp;
								if (randoms_correction)
									ax[to] += randoms[lo];
								rhs[lo + to * loop_var_par] = ax[to];
							}
							continue;
						}
						if (local_sino > 0.) {
							if (fp != 2) {
								for (int64_t to = 0LL; to < nBins; to++) {
									if (ax[to] < epps)
										ax[to] = epps;
									else
										ax[to] *= temp;
									if (randoms_correction)
										ax[to] += randoms[lo];
									yax[to] = Sino[lo + to * TOFSize] / ax[to];
								}
							}
							for (uint32_t k = 0; k < Np; k++) {
								val_rhs = TOFWeightsBP(D, nBins, dy, TOFVal, TOFCenter, sigma_x, yI, yax, epps, temp, val, rhs, tid + k * nBins);
#pragma omp atomic
								rhs[temp_ijk + k * Nx] += (val_rhs);
								if (no_norm == 0 && val > 0.) {
#pragma omp atomic
									Summ[temp_ijk + k * Nx] += val;
								}
							}
						}
						else {
							for (uint32_t k = 0; k < Np; k++) {
								TOFWeightsSumm(D, nBins, dy, TOFVal, TOFCenter, sigma_x, yI, epps, temp, val, tid + k * nBins);
								if (no_norm == 0 && val > 0.) {
#pragma omp atomic
									Summ[temp_ijk + k * Nx] += val;
								}
							}
						}
					}
					else {
						if (fp == 1) {
							for (uint32_t k = 0; k < Np; k++) {
								ax[0] += (element * osem_apu[temp_ijk + k * Nx]);
							}
							if (ax[0] < epps)
								ax[0] = epps;
							if (randoms_correction)
								ax[0] += randoms[lo];
							rhs[lo] = ax[0];
							continue;
						}
						if (local_sino > 0.) {
							if (fp != 2) {
								for (uint32_t k = 0; k < Np; k++) {
									ax[0] += (element * osem_apu[temp_ijk + k * Nx]);
								}
								if (ax[0] < epps)
									ax[0] = epps;
								if (randoms_correction)
									ax[0] += randoms[lo];
								yax[0] = local_sino / ax[0];
							}
							for (uint32_t k = 0; k < Np; k++) {
#pragma omp atomic
								rhs[temp_ijk + k * Nx] += (element * yax[0]);
								if (no_norm == 0) {
#pragma omp atomic
									Summ[temp_ijk + k * Nx] += element;
								}
							}
						}
						else {
							for (uint32_t k = 0; k < Np; k++) {
								if (no_norm == 0) {
#pragma omp atomic
									Summ[temp_ijk + k * Nx] += element;
								}
							}
						}
					}
				}
			}
		}
		else {
			int32_t tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;
			double txu = 0., tyu = 0., tzu = 0., tc = 0., tx0 = 1e8, ty0 = 1e8, tz0 = 1e8;
			bool skip = false;

			if (std::fabs(z_diff) < 1e-8) {
				tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));
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
			}
			else {
				skip = siddon_pre_loop_3D(bx, by, bz, x_diff, y_diff, z_diff, maxxx, maxyy, bzb, dx, dy, dz, Nx, Ny, Nz, tempi, tempj, tempk, tyu, txu, tzu,
					Np, TYPE, detectors, tc, iu, ju, ku, tx0, ty0, tz0);
			}

			const double LL = sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);

			double temp = 0.;
			double tx0_a = tx0, ty0_a = ty0, tz0_a = tz0, tc_a = tc;
			int32_t tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
			uint32_t tempijk = static_cast<uint32_t>(tempk) * Nyx + static_cast<uint32_t>(tempj) * Nx + static_cast<uint32_t>(tempi);

			if (TOF) {
				TOFDis(x_diff, y_diff, z_diff, tc, LL, D, DD);
			}

			for (uint32_t ii = 0; ii < Np; ii++) {
				if (tx0 < ty0 && tx0 < tz0) {
					ForwardProject(tx0, tc, txu, LL, attenuation_correction, jelppi, atten, tempijk, TOF, DD, nBins, TOFVal, TOFCenter,
						sigma_x, D, osem_apu, ax, epps, temp, tempi, iu, 1U, tid, ii);
				}
				else if (ty0 < tz0) {
					ForwardProject(ty0, tc, tyu, LL, attenuation_correction, jelppi, atten, tempijk, TOF, DD, nBins, TOFVal, TOFCenter,
						sigma_x, D, osem_apu, ax, epps, temp, tempj, ju, Nx, tid, ii);
				}
				else {
					ForwardProject(tz0, tc, tzu, LL, attenuation_correction, jelppi, atten, tempijk, TOF, DD, nBins, TOFVal, TOFCenter,
						sigma_x, D, osem_apu, ax, epps, temp, tempk, ku, Nyx, tid, ii);
				}
			}

			temp = 1. / temp;
			tx0 = tx0_a;
			ty0 = ty0_a;
			tz0 = tz0_a;
			tempi = tempi_a, tempj = tempj_a, tempk = tempk_a;
			tempijk = static_cast<uint32_t>(tempk) * Nyx + static_cast<uint32_t>(tempj) * Nx + static_cast<uint32_t>(tempi);
			tc = tc_a;
			if (attenuation_correction)
				temp *= exp(jelppi);
			if (normalization)
				temp *= norm_coef[lo];
			if (scatter)
				temp *= scatter_coef[lo];
			temp *= global_factor;
			D = DD;
			if (fp == 1) {
				if (TOF) {
					for (int64_t to = 0LL; to < nBins; to++) {
						if (ax[to] < epps)
							ax[to] = epps;
						else
							ax[to] *= temp;
						if (randoms_correction)
							ax[to] += randoms[lo];
						rhs[lo + to * loop_var_par] = ax[to];
					}
				}
				else {
					if (ax[0] < epps)
						ax[0] = epps;
					else
						ax[0] *= temp;
					if (randoms_correction)
						ax[0] += randoms[lo];
					rhs[lo] = ax[0];
				}
				continue;
			}

			if (local_sino != 0.) {
				if (fp != 2) {
					if (TOF) {
						for (int64_t to = 0LL; to < nBins; to++) {
							if (ax[to] < epps) {
								ax[to] = epps;
							}
							else {
								ax[to] *= temp;
							}
							if (randoms_correction)
								ax[to] += randoms[lo];
							yax[to] = Sino[lo + to * TOFSize] / ax[to];
						}
					}
					else {
						if (ax[0] < epps) {
							ax[0] = epps;
						}
						else {
							ax[0] *= temp;
						}
						if (randoms_correction)
							ax[0] += randoms[lo];
						yax[0] = local_sino / ax[0];
					}
				}
				for (uint32_t ii = 0; ii < Np; ii++) {
					if (tx0 < ty0 && tx0 < tz0) {
						backwardProjection(tx0, tc, txu, LL, tempijk, TOF, DD, nBins, TOFVal, TOFCenter, sigma_x, D, yax, epps, temp,
							iu, 1U, no_norm, rhs, Summ, tid, ii);
					}
					else if (ty0 < tz0) {
						backwardProjection(ty0, tc, tyu, LL, tempijk, TOF, DD, nBins, TOFVal, TOFCenter, sigma_x, D, yax, epps, temp,
							ju, Nx, no_norm, rhs, Summ, tid, ii);
					}
					else {
						backwardProjection(tz0, tc, tzu, LL, tempijk, TOF, DD, nBins, TOFVal, TOFCenter, sigma_x, D, yax, epps, temp,
							ku, Nyx, no_norm, rhs, Summ, tid, ii);
					}
				}
			}
			else {
				for (uint32_t ii = 0; ii < Np; ii++) {
					if (tx0 < ty0 && tx0 < tz0) {
						sensImage(tx0, tc, txu, LL, tempijk, TOF, DD, nBins, TOFVal, TOFCenter, sigma_x, D, epps, temp,
							iu, 1U, no_norm, Summ, tid, ii);
					}
					else if (ty0 < tz0) {
						sensImage(ty0, tc, tyu, LL, tempijk, TOF, DD, nBins, TOFVal, TOFCenter, sigma_x, D, epps, temp,
							ju, Nx, no_norm, Summ, tid, ii);
					}
					else {
						sensImage(tz0, tc, tzu, LL, tempijk, TOF, DD, nBins, TOFVal, TOFCenter, sigma_x, D, epps, temp,
							ku, Nyx, no_norm, Summ, tid, ii);
					}
				}
			}
		}
	}
}

void sequential_orth_siddon(const int64_t loop_var_par, const uint32_t size_x, const double zmax, double* Summ, double* rhs, const double maxyy,
	const double maxxx, const vector<double>& xx_vec, const double dy, const vector<double>& yy_vec, const double* atten, const double* norm_coef,
	const double* randoms, const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny,
	const uint32_t Nz, const double dx, const double dz, const double bx, const double by, const double bz, const bool attenuation_correction,
	const bool normalization, const bool randoms_correction, const uint16_t* lor1, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos,
	const double epps, const double* Sino, double* osem_apu, const uint16_t* L, const uint32_t* pseudos, const uint32_t pRows, const uint32_t det_per_ring,
	const bool raw, const double crystal_size_xy, double* x_center, double* y_center, const double* z_center, const double crystal_size_z,
	const bool no_norm, const uint32_t dec_v, const double global_factor, const uint8_t fp, const bool scatter, const double* scatter_coef, const uint32_t nCores) {

	if (nCores == 1U)
		setThreads();
	else
		omp_set_num_threads(nCores);

	const uint32_t Nyx = Ny * Nx;

	const double bzb = bz + static_cast<double>(Nz) * dz;

	size_t idx = 0ULL;
	vector<double> elements;
	vector<uint32_t> v_indices;

	size_t* indi = 0ULL;

#ifdef _OPENMP
	size_t threads = omp_get_max_threads();
#else
	size_t threads = 1ULL;
#endif
	std::vector<double> store_elements(threads * dec_v, 0.);
	std::vector<uint32_t> store_indices(threads * dec_v, 0u);

#pragma omp parallel for ordered schedule(dynamic)
	for (int64_t lo = 0LL; lo < loop_var_par; lo++) {


		const double local_sino = Sino[lo];
		if (no_norm && local_sino == 0.)
			continue;
#ifdef _OPENMP
		const uint32_t tid = omp_get_thread_num() * dec_v;
#else
		const uint32_t tid = 0U;
#endif
		Det detectors;
		double kerroin, length_;

		// Raw list-mode data
		if (raw) {
			get_detector_coordinates_raw(det_per_ring, x, y, z_det, detectors, L, lo, pseudos, pRows);
		}
		// Sinogram data
		else {
			get_detector_coordinates(x, y, z_det, size_x, detectors, xy_index, z_index, TotSinos, lo);
		}

		// Calculate the x, y and z distances of the detector pair
		double x_diff = (detectors.xd - detectors.xs);
		double y_diff = (detectors.yd - detectors.ys);
		const double z_diff = (detectors.zd - detectors.zs);

		// Load the number of voxels the LOR traverses (precomputed)
		uint32_t Np = static_cast<uint32_t>(lor1[lo]);
		double ax = 0., jelppi = 0., LL;
		uint8_t xyz = 0u;
		bool RHS = false, SUMMA = false;
		uint32_t ind = 0u;

		uint32_t N0 = Nx;
		uint32_t N1 = Ny;
		uint32_t N2 = 1u;
		uint32_t N3 = Nx;
		uint32_t N4 = Nz;

		double* xcenter = x_center;
		double* ycenter = y_center;

		if (fp == 2) {
			ax = osem_apu[lo];
		}

		if (crystal_size_z == 0.) {
			kerroin = norm(x_diff, y_diff, z_diff) * crystal_size_xy;
		}
		else {
			kerroin = norm(x_diff, y_diff, z_diff) * crystal_size_z;
		}
		
		if (fabs(z_diff) < 1e-8 && (fabs(y_diff) < 1e-8 || fabs(x_diff) < 1e-8)) {

			const uint32_t tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));

			if (fabs(y_diff) < 1e-8) {
				if (detectors.yd <= maxyy && detectors.yd >= by) {
					double temppi = detectors.xs;
					detectors.xs = detectors.ys;
					detectors.ys = temppi;
					if (crystal_size_z == 0.) {
						double temp = 0.;
						orth_distance_denominator_perpendicular_mfree(y_center, x_center[0], z_center, kerroin, temp, attenuation_correction, normalization, ax,
							by, detectors.yd, dy, Ny, Nx, tempk, atten, norm_coef, local_sino, Ny, 1u, osem_apu, detectors, y_diff, x_diff, z_diff, store_elements, store_indices, tid,
							ind, rhs, indi, lo, PRECOMPUTE, global_factor, scatter, scatter_coef);
						if (fp == 1) {
							if (ax == 0.)
								ax = epps;
							else
								ax *= temp;
							if (randoms_correction)
								ax += randoms[lo];
							rhs[lo] = ax;
							continue;
						}
						if (local_sino > 0.) {
							if (fp != 2) {
								nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
							}
							orth_distance_rhs_perpendicular_mfree(y_center, x_center[0], z_center, kerroin, temp, ax, by, detectors.yd, dy, Ny, Nx, tempk, Ny,
								1u, no_norm, rhs, Summ, true, false, detectors, y_diff, x_diff, z_diff, store_elements, store_indices, tid, ind, 0, indi);
						}
						else {
							orth_distance_rhs_perpendicular_mfree(y_center, x_center[0], z_center, kerroin, temp, ax, by, detectors.yd, dy, Ny, Nx, tempk, Ny,
								1u, no_norm, rhs, Summ, false, true, detectors, y_diff, x_diff, z_diff, store_elements, store_indices, tid, ind, 0, indi);
						}
					}
					else {
						double temp = 0.;
						orth_distance_denominator_perpendicular_mfree_3D(y_center, x_center[0], z_center, temp, attenuation_correction, normalization, ax,
							by, detectors.yd, dy, Ny, Nx, tempk, atten, norm_coef, local_sino, Ny, 1u, osem_apu, detectors, y_diff, x_diff, z_diff, kerroin, Nyx, Nz, store_elements, store_indices, tid,
							ind, rhs, indi, lo, PRECOMPUTE, global_factor, scatter, scatter_coef);
						if (local_sino > 0.) {
							if (fp != 2) {
								nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
							}
							orth_distance_rhs_perpendicular_mfree(y_center, x_center[0], z_center, kerroin, temp, ax, by, detectors.yd, dy, Ny, Nx, tempk, Ny,
								1u, no_norm, rhs, Summ, true, false, detectors, y_diff, x_diff, z_diff, store_elements, store_indices, tid, ind, 0, indi);
						}
						else {
							orth_distance_rhs_perpendicular_mfree(y_center, x_center[0], z_center, kerroin, temp, ax, by, detectors.yd, dy, Ny, Nx, tempk, Ny,
								1u, no_norm, rhs, Summ, false, true, detectors, y_diff, x_diff, z_diff, store_elements, store_indices, tid, ind, 0, indi);
						}
					}
				}
			}
			else if (fabs(x_diff) < 1e-8) {

				if (detectors.xd <= maxxx && detectors.xd >= bx) {
					if (crystal_size_z == 0.) {
						double temp = 0.;
						orth_distance_denominator_perpendicular_mfree(x_center, y_center[0], z_center, kerroin, temp, attenuation_correction, normalization, ax,
							bx, detectors.xd, dx, Nx, Ny, tempk, atten, norm_coef, local_sino, 1u, Nx, osem_apu, detectors, x_diff, y_diff, z_diff, store_elements, store_indices, tid,
							ind, rhs, indi, lo, PRECOMPUTE, global_factor, scatter, scatter_coef);
						if (fp == 1) {
							if (ax == 0.)
								ax = epps;
							else
								ax *= temp;
							if (randoms_correction)
								ax += randoms[lo];
							rhs[lo] = ax;
							continue;
						}
						if (local_sino > 0.) {
							if (fp != 2) {
								nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
							}
							orth_distance_rhs_perpendicular_mfree(x_center, y_center[0], z_center, kerroin, temp, ax, bx, detectors.xd, dx, Nx, Ny, tempk,
								1u, Nx, no_norm, rhs, Summ, true, false, detectors, x_diff, y_diff, z_diff, store_elements, store_indices, tid, ind, 0, indi);
						}
						else {
							orth_distance_rhs_perpendicular_mfree(x_center, y_center[0], z_center, kerroin, temp, ax, bx, detectors.xd, dx, Nx, Ny, tempk,
								1u, Nx, no_norm, rhs, Summ, false, true, detectors, x_diff, y_diff, z_diff, store_elements, store_indices, tid, ind, 0, indi);
						}
					}
					else {
						double temp = 0.;
						orth_distance_denominator_perpendicular_mfree_3D(x_center, y_center[0], z_center, temp, attenuation_correction, normalization, ax,
							bx, detectors.xd, dx, Nx, Ny, tempk, atten, norm_coef, local_sino, 1u, Nx, osem_apu, detectors, x_diff, y_diff, z_diff, kerroin, Nyx, Nz, store_elements, store_indices, tid,
							ind, rhs, indi, lo, PRECOMPUTE, global_factor, scatter, scatter_coef);
						if (local_sino > 0.) {
							nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
							orth_distance_rhs_perpendicular_mfree(x_center, y_center[0], z_center, kerroin, temp, ax, bx, detectors.xd, dx, Nx, Ny, tempk,
								1u, Nx, no_norm, rhs, Summ, true, false, detectors, x_diff, y_diff, z_diff, store_elements, store_indices, tid, ind, 0, indi);
						}
						else {
							orth_distance_rhs_perpendicular_mfree(x_center, y_center[0], z_center, kerroin, temp, ax, bx, detectors.xd, dx, Nx, Ny, tempk,
								1u, Nx, no_norm, rhs, Summ, false, true, detectors, x_diff, y_diff, z_diff, store_elements, store_indices, tid, ind, 0, indi);
						}
					}
				}
			}
		}
		else {
			int32_t tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;
			double txu = 0., tyu = 0., tzu = 0., tc = 0., tx0 = 1e8, ty0 = 1e8, tz0 = 1e8;
			bool skip = false;

			if (std::fabs(z_diff) < 1e-8) {
				tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));
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
				osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP, PRECOMPUTE, DISCARD, rhs, Summ, indi, elements, v_indices, idx,
				N1, N3, alku, iu, ju, loppu, store_elements, store_indices, tid, ind);

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
							osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP, PRECOMPUTE, DISCARD, rhs, Summ, indi, elements, v_indices, idx,
							N1, N3, alku, iu, ju, loppu, store_elements, store_indices, tid, ind);
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
						osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP, PRECOMPUTE, DISCARD, rhs, Summ, indi, elements, v_indices, idx,
						N1, N3, alku, iu, ju, loppu, store_elements, store_indices, tid, ind);
				}
			}

			temp = 1. / temp;
			if (attenuation_correction)
				temp *= exp(jelppi);
			if (normalization)
				temp *= norm_coef[lo];
			if (scatter)
				temp *= scatter_coef[lo];
			temp *= global_factor;

			if (fp == 1) {
				if (ax == 0.)
					ax = epps;
				else
					ax *= temp;
				if (randoms_correction)
					ax += randoms[lo];
				rhs[lo] = ax;
				continue;
			}
			if (local_sino > 0.) {
				if (fp != 2) {
					nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
				}
				RHS = true;
			}
			else
				SUMMA = true;

			orth_distance_3D_full(tempi, N0, N4, y_diff, x_diff, z_diff, ycenter, xcenter, z_center, temp, N2, tempj, tempk, local_sino, ax,
				osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP, PRECOMPUTE, DISCARD, rhs, Summ, indi, elements, v_indices, idx,
				N1, N3, alku, iu, ju, loppu, store_elements, store_indices, tid, ind);
		}
	}
}

void sequential_volume_siddon(const int64_t loop_var_par, const uint32_t size_x, const double zmax, double* Summ, double* rhs, const double maxyy,
	const double maxxx, const vector<double>& xx_vec, const double dy, const vector<double>& yy_vec, const double* atten, const double* norm_coef,
	const double* randoms, const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny,
	const uint32_t Nz, const double dx, const double dz, const double bx, const double by, const double bz, const bool attenuation_correction,
	const bool normalization, const bool randoms_correction, const uint16_t* lor1, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos,
	const double epps, const double* Sino, double* osem_apu, const uint16_t* L, const uint32_t* pseudos, const uint32_t pRows, const uint32_t det_per_ring,
	const bool raw, const double Vmax, double* x_center, double* y_center, const double* z_center, const double bmin, const double bmax, const double* V,
	const bool no_norm, const uint32_t dec_v, const double global_factor, const uint8_t fp, const bool scatter, const double* scatter_coef, const uint32_t nCores) {

	if (nCores == 1U)
		setThreads();
	else
		omp_set_num_threads(nCores);

	const uint32_t Nyx = Ny * Nx;

	const double bzb = bz + static_cast<double>(Nz) * dz;

	size_t idx = 0ULL;
	vector<double> elements;
	vector<uint32_t> v_indices;
	size_t* indi = 0ULL;

#ifdef _OPENMP
	size_t threads = omp_get_max_threads();
#else
	size_t threads = 1ULL;
#endif
	std::vector<double> store_elements(threads * dec_v, 0.);
	std::vector<uint32_t> store_indices(threads * dec_v, 0u);

#pragma omp parallel for ordered schedule(dynamic)
	for (int64_t lo = 0LL; lo < loop_var_par; lo++) {


		const double local_sino = Sino[lo];
		if (no_norm && local_sino == 0.)
			continue;

#ifdef _OPENMP
		const uint32_t tid = omp_get_thread_num() * dec_v;
#else
		const uint32_t tid = 0U;
#endif
		Det detectors;
		double kerroin, length_;

		// Raw list-mode data
		if (raw) {
			get_detector_coordinates_raw(det_per_ring, x, y, z_det, detectors, L, lo, pseudos, pRows);
		}
		// Sinogram data
		else {
			get_detector_coordinates(x, y, z_det, size_x, detectors, xy_index, z_index, TotSinos, lo);
		}

		// Calculate the x, y and z distances of the detector pair
		double x_diff = (detectors.xd - detectors.xs);
		double y_diff = (detectors.yd - detectors.ys);
		const double z_diff = (detectors.zd - detectors.zs);

		// Load the number of voxels the LOR traverses (precomputed)
		uint32_t Np = static_cast<uint32_t>(lor1[lo]);
		double ax = 0., jelppi = 0., LL;
		int8_t start = 1;
		uint8_t xyz = 0u;
		uint8_t xyz_w = 0u;
		bool RHS = false, SUMMA = false;
		uint32_t ind = 0u;

		uint32_t N0 = Nx;
		uint32_t N1 = Ny;
		uint32_t N2 = 1u;
		uint32_t N3 = Nx;
		uint32_t N4 = Nz;

		double* xcenter = x_center;
		double* ycenter = y_center;

		if (fp == 2) {
			ax = osem_apu[lo];
		}

		kerroin = norm(x_diff, y_diff, z_diff);

		if (fabs(z_diff) < 1e-8 && (fabs(y_diff) < 1e-8 || fabs(x_diff) < 1e-8)) {

			const uint32_t tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));

			if (fabs(y_diff) < 1e-8) {

				if (detectors.yd <= maxyy && detectors.yd >= by) {
					double temppi = detectors.xs;
					detectors.xs = detectors.ys;
					detectors.ys = temppi;
					double temp = 0.;
					volume_distance_denominator_perpendicular_mfree_3D(y_center, x_center[0], z_center, temp, attenuation_correction, normalization, ax,
						by, detectors.yd, dy, Ny, Nx, tempk, atten, norm_coef, local_sino, Ny, 1u, osem_apu, detectors, y_diff, x_diff, z_diff, kerroin, Nyx, Nz, store_elements, store_indices, tid,
						ind, rhs, indi, lo, PRECOMPUTE, global_factor, bmax, bmin, Vmax, V, scatter, scatter_coef);
					if (fp == 1) {
						if (ax == 0.)
							ax = epps;
						else
							ax *= temp;
						if (randoms_correction)
							ax += randoms[lo];
						rhs[lo] = ax;
						continue;
					}
					if (local_sino > 0.) {
						if (fp != 2) {
							nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
						}
						orth_distance_rhs_perpendicular_mfree(y_center, x_center[0], z_center, kerroin, temp, ax, by, detectors.yd, dy, Ny, Nx, tempk, Ny,
							1u, no_norm, rhs, Summ, true, false, detectors, y_diff, x_diff, z_diff, store_elements, store_indices, tid, ind, 0, indi);
					}
					else {
						orth_distance_rhs_perpendicular_mfree(y_center, x_center[0], z_center, kerroin, temp, ax, by, detectors.yd, dy, Ny, Nx, tempk, Ny,
							1u, no_norm, rhs, Summ, false, true, detectors, y_diff, x_diff, z_diff, store_elements, store_indices, tid, ind, 0, indi);
					}
				}
			}
			else if (fabs(x_diff) < 1e-8) {

				if (detectors.xd <= maxxx && detectors.xd >= bx) {
					double temp = 0.;
					volume_distance_denominator_perpendicular_mfree_3D(x_center, y_center[0], z_center, temp, attenuation_correction, normalization, ax,
						bx, detectors.xd, dx, Nx, Ny, tempk, atten, norm_coef, local_sino, 1u, Nx, osem_apu, detectors, x_diff, y_diff, z_diff, kerroin, Nyx, Nz, store_elements, store_indices, tid,
						ind, rhs, indi, lo, PRECOMPUTE, global_factor, bmax, bmin, Vmax, V, scatter, scatter_coef);
					if (fp == 1) {
						if (ax == 0.)
							ax = epps;
						else
							ax *= temp;
						if (randoms_correction)
							ax += randoms[lo];
						rhs[lo] = ax;
						continue;
					}
					if (local_sino > 0.) {
						if (fp != 2) {
							nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
						}
						orth_distance_rhs_perpendicular_mfree(x_center, y_center[0], z_center, kerroin, temp, ax, bx, detectors.xd, dx, Nx, Ny, tempk,
							1u, Nx, no_norm, rhs, Summ, true, false, detectors, x_diff, y_diff, z_diff, store_elements, store_indices, tid, ind, 0, indi);
					}
					else {
						orth_distance_rhs_perpendicular_mfree(x_center, y_center[0], z_center, kerroin, temp, ax, bx, detectors.xd, dx, Nx, Ny, tempk,
							1u, Nx, no_norm, rhs, Summ, false, true, detectors, x_diff, y_diff, z_diff, store_elements, store_indices, tid, ind, 0, indi);
					}
				}
			}
		}
		else {
			int32_t tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;
			double txu = 0., tyu = 0., tzu = 0., tc = 0., tx0 = 1e8, ty0 = 1e8, tz0 = 1e8;
			bool skip = false;

			if (std::fabs(z_diff) < 1e-8) {
				tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));
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
			alku = Nz;
			loppu = 0;
			if (ku > 0) {
				alku = tempk + 1;
			}
			else if (ku < 0) {
				loppu = tempk;
			}
			volume_distance_3D_full(tempi, N0, N4, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, N2, tempj, tempk, local_sino,
				ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP, PRECOMPUTE, DISCARD, rhs, Summ, indi, elements, v_indices,
				idx, N1, N3, alku, iu, ju, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);

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
						volume_distance_3D_full(tempi, N0, N4, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, N2, tempj, tempk, local_sino,
							ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP, PRECOMPUTE, DISCARD, rhs, Summ, indi, elements, v_indices,
							idx, N1, N3, alku, iu, ju, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);
					}
				}
			}
			if (xyz < 3 && std::fabs(z_diff) >= 1e-8) {
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
					volume_distance_3D_full(tempi, N0, N4, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, N2, tempj, tempk, local_sino,
						ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP, PRECOMPUTE, DISCARD, rhs, Summ, indi, elements, v_indices,
						idx, N1, N3, alku, iu, ju, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);
				}
			}

			temp = 1. / temp;
			if (attenuation_correction)
				temp *= exp(jelppi);
			if (normalization)
				temp *= norm_coef[lo];
			if (scatter)
				temp *= scatter_coef[lo];
			temp *= global_factor;

			if (fp == 1) {
				if (ax == 0.)
					ax = epps;
				else
					ax *= temp;
				if (randoms_correction)
					ax += randoms[lo];
				rhs[lo] = ax;
				continue;
			}
			if (local_sino > 0.) {
				if (fp != 2) {
					nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
				}
				RHS = true;
			}
			else
				SUMMA = true;

			volume_distance_3D_full(tempi, N0, N4, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, N2, tempj, tempk, local_sino,
				ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP, PRECOMPUTE, DISCARD, rhs, Summ, indi, elements, v_indices,
				idx, N1, N3, alku, iu, ju, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);
		}
	}
}