/**************************************************************************
* Implements both the improved Siddon's algorithm and Orthogonal Siddon's 
* algorithm for OMEGA (Implementation 4).
* Determines which LORs intercept the FOV on-the-fly (slower). Improved
* Siddon can use n rays.
*
* Uses OpenMP for parallellization. If OpenMP is not available, the code
* is serial with no parallellization.
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
const static int TYPE = 0;

// Whether to use the OpenMP code or not
const static bool OMP = true;

// Using non-OpenMP with either precomputation or without
const static bool PRECOMPUTE = false;

const static bool DISCARD = false;


using namespace std;

// Improved multi-ray Siddon
void sequential_improved_siddon_no_precompute(const int64_t loop_var_par, const uint32_t size_x, const double zmax, double* Summ, double* rhs, const double maxyy,
	const double maxxx,	const vector<double>& xx_vec, const double dy, const vector<double>& yy_vec, const double* atten, const double* norm_coef, 
	const double* randoms, const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, 
	const uint32_t Nz, const double dx, const double dz, const double bx, const double by, const double bz, const bool attenuation_correction, 
	const bool normalization, const bool randoms_correction, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos,
	const double epps, const double* Sino, double* osem_apu, const uint16_t* L, const uint32_t* pseudos, const size_t pRows, const uint32_t det_per_ring,
	const bool raw, const double cr_pz, const bool no_norm, const uint16_t n_rays, const uint16_t n_rays3D, const double global_factor, const uint8_t fp, 
	const bool list_mode_format, const bool scatter, const double* scatter_coef, const bool TOF, const int64_t TOFSize, const double sigma_x, const double* TOFCenter,
	const int64_t nBins, const uint32_t dec_v) {

	setThreads();

	// Size of single 2D image
	const uint32_t Nyx = Ny * Nx;

	uint32_t bin = 0U;

	// Distance of the last slice from origin
	const double bzb = bz + static_cast<double>(Nz) * dz;

	// Distance between rays in multi-ray Siddon
	const double dc_z = cr_pz / static_cast<double>(n_rays3D + 1);

	int64_t lo = 0LL;
#ifdef _OPENMP
	size_t threads = omp_get_max_threads();
#else
	size_t threads = 1ULL;
#endif

	const size_t nRays = static_cast<size_t>(n_rays) * static_cast<size_t>(n_rays3D);
	vector<double> TOFVal(nRays * nBins * dec_v * threads, 0.);

#pragma omp parallel for schedule(dynamic)
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

		// Form vectors that store the necessary multi-ray information
		vector<int32_t> tempi_a(nRays, 0);
		vector<int32_t> tempj_a(nRays, 0);
		vector<int32_t> tempk_a(nRays, 0);
		vector<int32_t> iu_a(nRays, 0);
		vector<int32_t> ju_a(nRays, 0);
		vector<int32_t> ku_a(nRays, 0);
		vector<double> tx0_a(nRays, 0.);
		vector<double> ty0_a(nRays, 0.);
		vector<double> tz0_a(nRays, 0.);
		vector<double> tc_a(nRays, 0.);
		vector<double> txu_a(nRays, 0.);
		vector<double> tyu_a(nRays, 0.);
		vector<double> tzu_a(nRays, 0.);
		vector<double> x_diff(nRays, 0.);
		vector<double> y_diff(nRays, 0.);
		vector<double> z_diff(nRays, 0.);
		vector<double> LL(nRays, 0.);
		vector<double> D_a(nRays, 0.);
		vector<uint32_t> Np_n(nRays, 0u);


		vector<double> ax(nBins, 0.);
		vector<double> yax(nBins, 0.);
		//vector<double> TOFVal(nRays * nBins, 0.);

		double temp = 0.;
		double jelppi = 0.;
		double D = 0., DD = 0.;
		double xI = 0., yI = 0., zI = 0.;

		vector<bool> pass(nRays, false);

		if (fp == 2) {
			for (int64_t to = 0LL; to < nBins; to++)
				yax[to] = osem_apu[lo + to * loop_var_par];
		}

		// Loop through the rays
		for (uint16_t lor = 0u; lor < nRays; lor++) {

#ifdef _OPENMP
			const int64_t tid = omp_get_thread_num() * dec_v * nRays * nBins + static_cast<int64_t>(lor) * nBins * static_cast<int64_t>(dec_v);
#else
			const int64_t tid = static_cast<int64_t>(lor) * nBins * static_cast<int64_t>(dec_v);
#endif

			Det detectors;

			// Raw list-mode data
			if (raw) {
				// Pure list-mode format (e.g. event-by-event)
				if (list_mode_format)
					get_detector_coordinates_raw(det_per_ring, x, y, z_det, detectors, L, lo, pseudos, pRows, list_mode_format);
				// Raw list-mode format
				else
					get_detector_coordinates_raw_N(det_per_ring, x, y, z_det, detectors, L, lo, pseudos, pRows, lor + 1u, dc_z, n_rays, n_rays3D);
			}
			// Sinogram data
			else {
				get_detector_coordinates_mr(x, y, z_det, size_x, detectors, xy_index, z_index, TotSinos, lo, lor + 1u, dc_z, n_rays, n_rays3D);
			}

			// Calculate the x, y and z distances of the detector pair
			y_diff[lor] = (detectors.yd - detectors.ys);
			x_diff[lor] = (detectors.xd - detectors.xs);
			z_diff[lor] = (detectors.zd - detectors.zs);
			// Skip certain cases (e.g. if the x- and y-coordinates are the same for both detectors, LOR between detector n and n)
			if ((y_diff[lor] == 0. && x_diff[lor] == 0. && z_diff[lor] == 0.) || (y_diff[lor] == 0. && x_diff[lor] == 0.)) {
				continue;
			}

			// Number of voxels the ray traverses
			uint32_t Np = 0U;


			if (fabs(z_diff[lor]) < 1e-8 && (fabs(y_diff[lor]) < 1e-8 || fabs(x_diff[lor]) < 1e-8)) {

				// Ring number
				const uint32_t tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));

				// Detectors are perpendicular
				// Siddon cannot be applied --> trivial to compute
				if (fabs(y_diff[lor]) < 1e-8) {

					if (detectors.yd <= maxyy && detectors.yd >= by) {
						int32_t apu = 0;

						// Determine the starting coordinate, ray length and compute attenuation effects
						double element = perpendicular_elements_multiray(Ny, detectors.yd, yy_vec, dx, tempk, Nx, Ny, atten, attenuation_correction,
							apu, 1u, jelppi);

						if (TOF) {
							xI = (dx * Nx) / 2.;
							if (x_diff[lor] > 0.)
								xI = -xI;
							D = xI;
							for (uint32_t k = 0; k < Nx; k++) {
								TOFWeightsFP(D, nBins, dx, TOFVal, TOFCenter, sigma_x, xI, osem_apu, apu + k, ax, epps, tid + static_cast<int64_t>(k) * nBins);
							}
						}
						else {
							// Forward projection
							for (uint32_t k = 0; k < Nx; k++) {
								ax[0] += (dx * osem_apu[apu + k]);
							}
						}
						// Total length
						temp += element;
						tempk_a[lor] = apu;
						pass[lor] = true;
					}
				}
				else if (fabs(x_diff[lor]) < 1e-8) {

					if (detectors.xd <= maxxx && detectors.xd >= bx) {
						int32_t apu = 0;
						double element = perpendicular_elements_multiray(1u, detectors.xd, xx_vec, dy, tempk, Ny, Nx, atten, attenuation_correction,
							apu, Nx, jelppi);

						if (TOF) {
							yI = (dy * Ny) / 2.;
							if (y_diff[lor] > 0.)
								yI = -yI;
							D = yI;
							for (uint32_t k = 0; k < Ny; k++) {
								TOFWeightsFP(D, nBins, dy, TOFVal, TOFCenter, sigma_x, yI, osem_apu, apu + k * Nx, ax, epps, tid + static_cast<int64_t>(k) * nBins);
							}
						}
						else {
							// Forward projection
							for (uint32_t k = 0; k < Ny; k++) {
								ax[0] += (dy * osem_apu[apu + k * Nx]);
							}
						}
						// Total length
						temp += element;
						tempk_a[lor] = apu;
						pass[lor] = true;
					}
				}
			}
			else {
				int32_t tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;
				double txu = 0., tyu = 0., tzu = 0., tc = 0., tx0 = 1e8, ty0 = 1e8, tz0 = 1e8;
				bool skip = false;

				// Determine the above values and whether the ray intersects the FOV
				// Both detectors are on the same ring, but not perpendicular
				if (std::fabs(z_diff[lor]) < 1e-8) {
					tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));
					skip = siddon_pre_loop_2D(bx, by, x_diff[lor], y_diff[lor], maxxx, maxyy, dx, dy, Nx, Ny, tempi, tempj, txu, tyu, Np, TYPE,
						detectors.ys, detectors.xs, detectors.yd, detectors.xd, tc, iu, ju, tx0, ty0);
				}
				//Detectors on different rings (e.g. oblique sinograms)
				else if (std::fabs(y_diff[lor]) < 1e-8) {
					skip = siddon_pre_loop_2D(bx, bz, x_diff[lor], z_diff[lor], maxxx, bzb, dx, dz, Nx, Nz, tempi, tempk, txu, tzu, Np, TYPE,
						detectors.zs, detectors.xs, detectors.zd, detectors.xd, tc, iu, ku, tx0, tz0);
					if (detectors.yd > maxyy || detectors.yd < by)
						skip = true;
					tempj = perpendicular_start(by, detectors.yd, dy, Ny);
				}
				else if (std::fabs(x_diff[lor]) < 1e-8) {
					skip = siddon_pre_loop_2D(by, bz, y_diff[lor], z_diff[lor], maxyy, bzb, dy, dz, Ny, Nz, tempj, tempk, tyu, tzu, Np, TYPE,
						detectors.zs, detectors.ys, detectors.zd, detectors.yd, tc, ju, ku, ty0, tz0);
					if (detectors.xd > maxxx || detectors.xd < bx)
						skip = true;
					tempi = perpendicular_start(bx, detectors.xd, dx, Nx);
				}
				else {
					skip = siddon_pre_loop_3D(bx, by, bz, x_diff[lor], y_diff[lor], z_diff[lor], maxxx, maxyy, bzb, dx, dy, dz, Nx, Ny, Nz, tempi, tempj, tempk, tyu, txu, tzu,
						Np, TYPE, detectors, tc, iu, ju, ku, tx0, ty0, tz0);
				}

				// Skip if the LOR does not intersect with the FOV
				if (skip) {
					continue;
				}

				// Save the total length of the ray
				LL[lor] = sqrt(x_diff[lor] * x_diff[lor] + y_diff[lor] * y_diff[lor] + z_diff[lor] * z_diff[lor]);

				tempi_a[lor] = tempi, tempj_a[lor] = tempj, tempk_a[lor] = tempk;
				tx0_a[lor] = tx0, ty0_a[lor] = ty0, tz0_a[lor] = tz0, tc_a[lor] = tc;
				txu_a[lor] = txu, tyu_a[lor] = tyu, tzu_a[lor] = tzu;
				iu_a[lor] = iu, ju_a[lor] = ju, ku_a[lor] = ku;
				uint32_t tempijk = static_cast<uint32_t>(tempk) * Nyx + static_cast<uint32_t>(tempj) * Nx + static_cast<uint32_t>(tempi);

				if (TOF) {
					xI = x_diff[lor] * tc;
					yI = y_diff[lor] * tc;
					zI = z_diff[lor] * tc;
					D = std::sqrt(xI * xI + yI * yI + zI * zI) - LL[lor] / 2.;
					D_a[lor] = D;
					DD = D;
				}

				 //Compute the total distance traveled by this ray in the FOV
				for (uint32_t ii = 0; ii < Np; ii++) {
					if (tx0 < ty0 && tx0 < tz0) {
						ForwardProject(tx0, tc, txu, LL[lor], attenuation_correction, jelppi, atten, tempijk, TOF, DD, nBins, TOFVal, TOFCenter,
							sigma_x, D, osem_apu, ax, epps, temp, tempi, iu, 1U, tid, ii);
					}
					else if (ty0 < tz0) {
						ForwardProject(ty0, tc, tyu, LL[lor], attenuation_correction, jelppi, atten, tempijk, TOF, DD, nBins, TOFVal, TOFCenter,
							sigma_x, D, osem_apu, ax, epps, temp, tempj, ju, Nx, tid, ii);
					}
					else {
						ForwardProject(tz0, tc, tzu, LL[lor], attenuation_correction, jelppi, atten, tempijk, TOF, DD, nBins, TOFVal, TOFCenter,
							sigma_x, D, osem_apu, ax, epps, temp, tempk, ku, Nyx, tid, ii);
					}
					// Number of voxels traversed
					Np_n[lor]++;
					// Check if ray has left the FOV
					if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= static_cast<int32_t>(Nx) || tempj >= static_cast<int32_t>(Ny) || tempk >= static_cast<int32_t>(Nz))
						break;
				}
				// This ray passed the FOV
				pass[lor] = true;
			}
		}

		bool alku = true;
		//double yax = 0.;

		// Compute the probabilities for the current LOR
		// Sum all the rays together
		for (uint16_t lor = 0u; lor < static_cast<size_t>(n_rays) * static_cast<size_t>(n_rays3D); lor++) {

			if (pass[lor]) {

				// Compute only before the first ray
				if (alku) {
					// To compute probability of each voxel interaction
					temp = 1. / temp;
					// Average attenuation
					if (attenuation_correction) {
						double n_r_summa = 0.;
						for (uint16_t ln_r = 0u; ln_r < static_cast<size_t>(n_rays) * static_cast<size_t>(n_rays3D); ln_r++)
							n_r_summa += static_cast<double>(pass[ln_r]);
						temp *= exp(jelppi / n_r_summa);
					}
					// Include normalization
					if (normalization)
						temp *= norm_coef[lo];
					if (scatter)
						temp *= scatter_coef[lo];
					// Global correction factor
					temp *= global_factor;
					// Special, only forward projection, case
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
						break;
					}
					// Forward projection when backprojection is also computed
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
					}
					alku = false;
				}
#ifdef _OPENMP
				const int64_t tid = omp_get_thread_num() * dec_v * nRays * nBins + static_cast<int64_t>(lor) * nBins * static_cast<int64_t>(dec_v);
#else
				const int64_t tid = static_cast<int64_t>(lor) * nBins * static_cast<int64_t>(dec_v);
#endif

				if (fabs(z_diff[lor]) < 1e-8 && (fabs(y_diff[lor]) < 1e-8 || fabs(x_diff[lor]) < 1e-8)) {

					if (fabs(y_diff[lor]) < 1e-8) {

						if (local_sino != 0.) {
							for (uint32_t k = 0; k < Nx; k++) {
								// "Right-hand side", backprojection
								double val = dx;
								double val_rhs = 0.;

								if (TOF) {
									if (k == 0) {
										xI = (dx * Nx - x_diff[lor]) / 2.;
										D = xI;
									}
									val_rhs = TOFWeightsBP(D, nBins, dx, TOFVal, TOFCenter, sigma_x, xI, yax, epps, temp, val, rhs, tid + k * nBins);
								}
								else {
									val *= temp;
									val_rhs = val * yax[0];
								}
#pragma omp atomic
								rhs[tempk_a[lor] + k] += (val_rhs);
								if (no_norm == 0 && val > 0.) {
#pragma omp atomic
									Summ[tempk_a[lor] + k] += val;
								}
							}
						}
						else {
							for (uint32_t k = 0; k < Nx; k++) {
								double val = dx;

								if (TOF) {
									if (k == 0) {
										xI = (dx * Nx - x_diff[lor]) / 2.;
										D = xI;
									}
									TOFWeightsSumm(D, nBins, dx, TOFVal, TOFCenter, sigma_x, xI, epps, temp, val, tid + k * nBins);
								}
								else
									val *= temp;

								if (no_norm == 0 && val > 0.) {
#pragma omp atomic
									Summ[tempk_a[lor] + k] += val;
								}
							}
						}
					}
					else if (fabs(x_diff[lor]) < 1e-8) {
						if (local_sino != 0.) {
							for (uint32_t k = 0; k < Ny; k++) {
								double val = dy;
								double val_rhs = 0.;

								if (TOF) {
									if (k == 0) {
										yI = (dy * Ny - y_diff[lor]) / 2.;
										D = yI;
									}
									val_rhs = TOFWeightsBP(D, nBins, dy, TOFVal, TOFCenter, sigma_x, yI, yax, epps, temp, val, rhs, tid + k * nBins);
								}
								else {
									val *= temp;
									val_rhs = val * yax[0];
								}
#pragma omp atomic
								rhs[tempk_a[lor] + k * Nx] += (val_rhs);
								if (no_norm == 0 && val > 0.) {
#pragma omp atomic
									Summ[tempk_a[lor] + k * Nx] += val;
								}
							}
						}
						else {
							for (uint32_t k = 0; k < Ny; k++) {
								double val = dy;

								if (TOF) {
									if (k == 0) {
										yI = (dy * Ny - y_diff[lor]) / 2.;
										D = yI;
									}
									TOFWeightsSumm(D, nBins, dy, TOFVal, TOFCenter, sigma_x, yI, epps, temp, val, tid + k * nBins);
								}
								else
									val *= temp;

								if (no_norm == 0 && val > 0.) {
#pragma omp atomic
									Summ[tempk_a[lor] + k * Nx] += val;
								}
							}
						}
					}
				}
				else {
					double tx0 = tx0_a[lor];
					double ty0 = ty0_a[lor];
					double tz0 = tz0_a[lor];
					const double txu = txu_a[lor];
					const double tyu = tyu_a[lor];
					const double tzu = tzu_a[lor];
					const int32_t tempi = tempi_a[lor];
					const int32_t tempj = tempj_a[lor];
					const int32_t tempk = tempk_a[lor];
					const int32_t iu = iu_a[lor];
					const int32_t ju = ju_a[lor];
					const int32_t ku = ku_a[lor];
					double tc = tc_a[lor];
					uint32_t tempijk = static_cast<uint32_t>(tempk) * Nyx + static_cast<uint32_t>(tempj) * Nx + static_cast<uint32_t>(tempi);
					D = D_a[lor];
					DD = D;

					if (local_sino != 0.) {
						for (uint32_t ii = 0; ii < Np_n[lor]; ii++) {
							if (tx0 < ty0 && tx0 < tz0) {
								backwardProjection(tx0, tc, txu, LL[lor], tempijk, TOF, DD, nBins, TOFVal, TOFCenter, sigma_x, D, yax, epps, temp,
									iu, 1U, no_norm, rhs, Summ, tid, ii);
							}
							else if (ty0 < tz0) {
								backwardProjection(ty0, tc, tyu, LL[lor], tempijk, TOF, DD, nBins, TOFVal, TOFCenter, sigma_x, D, yax, epps, temp,
									ju, Nx, no_norm, rhs, Summ, tid, ii);
							}
							else {
								backwardProjection(tz0, tc, tzu, LL[lor], tempijk, TOF, DD, nBins, TOFVal, TOFCenter, sigma_x, D, yax, epps, temp,
									ku, Nyx, no_norm, rhs, Summ, tid, ii);
							}
						}
					}
					else {
						for (uint32_t ii = 0; ii < Np_n[lor]; ii++) {
							if (tx0 < ty0 && tx0 < tz0) {
								sensImage(tx0, tc, txu, LL[lor], tempijk, TOF, DD, nBins, TOFVal, TOFCenter, sigma_x, D, epps, temp,
									iu, 1U, no_norm, Summ, tid, ii);
							}
							else if (ty0 < tz0) {
								sensImage(ty0, tc, tyu, LL[lor], tempijk, TOF, DD, nBins, TOFVal, TOFCenter, sigma_x, D, epps, temp,
									ju, Nx, no_norm, Summ, tid, ii);
							}
							else {
								sensImage(tz0, tc, tzu, LL[lor], tempijk, TOF, DD, nBins, TOFVal, TOFCenter, sigma_x, D, epps, temp,
									ku, Nyx, no_norm, Summ, tid, ii);
							}
						}
					}
				}
			}
		}
	}
}

// Orthogonal ray tracer
void sequential_orth_siddon_no_precomp(const int64_t loop_var_par, const uint32_t size_x, const double zmax, double* Summ, double* rhs, const double maxyy,
	const double maxxx, const vector<double>& xx_vec, const double dy, const vector<double>& yy_vec, const double* atten, const double* norm_coef,
	const double* randoms, const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny,
	const uint32_t Nz, const double dx, const double dz, const double bx, const double by, const double bz, const bool attenuation_correction,
	const bool normalization, const bool randoms_correction, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos,
	const double epps, const double* Sino, double* osem_apu, const uint16_t* L, const uint32_t* pseudos, const uint32_t pRows, const uint32_t det_per_ring,
	const bool raw, const double crystal_size_xy, double* x_center, double* y_center, const double* z_center, const double crystal_size_z,
	const bool no_norm, const uint32_t dec_v, const double global_factor, const uint8_t fp, const bool list_mode_format, const bool scatter, const double* scatter_coef) {

	setThreads();

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

#pragma omp parallel for schedule(dynamic)
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
		double kerroin;

		// Raw list-mode data
		if (raw) {
			get_detector_coordinates_raw(det_per_ring, x, y, z_det, detectors, L, lo, pseudos, pRows, list_mode_format);
		}
		// Sinogram data
		else {
			get_detector_coordinates(x, y, z_det, size_x, detectors, xy_index, z_index, TotSinos, lo);
		}

		// Calculate the x, y and z distances of the detector pair
		double x_diff = (detectors.xd - detectors.xs);
		double y_diff = (detectors.yd - detectors.ys);
		const double z_diff = (detectors.zd - detectors.zs);
		if ((y_diff == 0. && x_diff == 0. && z_diff == 0.) || (y_diff == 0. && x_diff == 0.))
			continue;

		double ax = 0., jelppi = 0., LL;
		uint32_t Np = 0u;
		uint32_t Np_n = 0u;
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

			// Ray tracing only needed for attenuation
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
				Np_n++;
				if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= static_cast<int32_t>(N0) || tempj >= static_cast<int32_t>(N1) || tempk >= static_cast<int32_t>(N4)) {
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
					break;
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

void sequential_volume_siddon_no_precomp(const int64_t loop_var_par, const uint32_t size_x, const double zmax, double* Summ, double* rhs, const double maxyy,
	const double maxxx, const vector<double>& xx_vec, const double dy, const vector<double>& yy_vec, const double* atten, const double* norm_coef,
	const double* randoms, const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny,
	const uint32_t Nz, const double dx, const double dz, const double bx, const double by, const double bz, const bool attenuation_correction,
	const bool normalization, const bool randoms_correction, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos,
	const double epps, const double* Sino, double* osem_apu, const uint16_t* L, const uint32_t* pseudos, const uint32_t pRows, const uint32_t det_per_ring,
	const bool raw, const double Vmax, double* x_center, double* y_center, const double* z_center, const double bmin, const double bmax, const double* V,
	const bool no_norm, const uint32_t dec_v, const double global_factor, const uint8_t fp, const bool list_mode_format, const bool scatter, const double* scatter_coef) {

	setThreads();

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

#pragma omp parallel for schedule(dynamic)
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
		double kerroin;

		// Raw list-mode data
		if (raw) {
			get_detector_coordinates_raw(det_per_ring, x, y, z_det, detectors, L, lo, pseudos, pRows, list_mode_format);
		}
		// Sinogram data
		else {
			get_detector_coordinates(x, y, z_det, size_x, detectors, xy_index, z_index, TotSinos, lo);
		}

		// Calculate the x, y and z distances of the detector pair
		double x_diff = (detectors.xd - detectors.xs);
		double y_diff = (detectors.yd - detectors.ys);
		const double z_diff = (detectors.zd - detectors.zs);
		if ((y_diff == 0. && x_diff == 0. && z_diff == 0.) || (y_diff == 0. && x_diff == 0.))
			continue;

		double ax = 0., jelppi = 0., LL;
		uint32_t Np = 0u;
		uint32_t Np_n = 0u;
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

			if (skip)
				continue;

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
				Np_n++;
				if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= static_cast<int32_t>(N0) || tempj >= static_cast<int32_t>(N1) || tempk >= static_cast<int32_t>(N4)) {
					if (xyz < 3 && std::fabs(z_diff) >= 1e-8) {
						if (xyz == 1)
							tempi -= iu;
						else if (xyz == 2)
							tempj -= ju;
						if ((tempk >= (Nz - 1) && ku > 0) || (tempk <= 0 && ku < 0))
							break;
						else
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
					break;
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