/**************************************************************************
* Implements both the improved Siddon's algorithm and Orthogonal Siddon's
* algorithm for OMEGA (Implementation 4).
* Uses the precomputed data (faster).
*
* Uses OpenMP for parallelization. If OpenMP is not available, the code
* is serial with no parallelization.
*
* Copyright (C) 2019 Ville-Veikko Wettenhovi
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
constexpr int TYPE = 1;

// Whether to use the OpenMP code or not
constexpr bool OMP = true;

// Using non-OpenMP with either precomputation or without
constexpr bool PRECOMPUTE = false;

constexpr bool DISCARD = false;

using namespace std;

void sequential_improved_siddon(const int64_t loop_var_par, const uint32_t size_x, const double zmax, double* Summ, double* rhs, const double maxyy,
	const double maxxx,	const vector<double>& xx_vec, const double dy, const vector<double>& yy_vec, const double* atten, const double* norm_coef, 
	const double* randoms, const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, 
	const uint32_t Nz, const double dx, const double dz, const double bx, const double by, const double bz, const bool attenuation_correction, 
	const bool normalization,  const bool randoms_correction, const uint16_t* lor1, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos, 
	const double epps,  const double* Sino, double* osem_apu, const uint16_t* L, const uint32_t* pseudos, const uint32_t pRows, const uint32_t det_per_ring, 
	const bool raw, const bool no_norm, const double global_factor, const bool fp, const bool scatter, const double* scatter_coef) {

	setThreads();

	const uint32_t Nyx = Ny * Nx;

	const double bzb = bz + static_cast<double>(Nz) * dz;

#pragma omp parallel for schedule(dynamic)
	for (int64_t lo = 0LL; lo < loop_var_par; lo++) {

		const double local_sino = Sino[lo];
		if (no_norm && local_sino == 0.)
			continue;
		Det detectors;

		// Raw list-mode data
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
		double ax = 0., jelppi = 0.;

		if (fabs(z_diff) < 1e-8) {

			const uint32_t tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));

			if (fabs(y_diff) < 1e-8) {

				if (detectors.yd <= maxyy && detectors.yd >= by) {
					uint32_t temp_ijk = 0;

					const double element = perpendicular_elements(Ny, detectors.yd, yy_vec, dx, tempk, Nx, Ny, atten, norm_coef, attenuation_correction, 
						normalization, temp_ijk, 1u, lo, global_factor, scatter, scatter_coef);

					if (fp) {
						for (uint32_t k = 0; k < Np; k++) {
							ax += (element * osem_apu[temp_ijk + k]);
						}
						if (ax == 0.)
							ax = epps;
						if (randoms_correction)
							ax += randoms[lo];
						rhs[lo] = ax;
						continue;
					}
					if (local_sino > 0.) {
						for (uint32_t k = 0; k < Np; k++) {
							ax += (element * osem_apu[temp_ijk + k]);
						}
						if (ax == 0.)
							ax = epps;
						if (randoms_correction)
							ax += randoms[lo];
						const double yax = local_sino / ax;
						for (uint32_t k = 0; k < Np; k++) {
#pragma omp atomic
							rhs[temp_ijk + k] += (element * yax);
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
			else if (fabs(x_diff) < 1e-8) {

				if (detectors.xd <= maxxx && detectors.xd >= bx) {
					uint32_t temp_ijk = 0;

					const double element = perpendicular_elements(1, detectors.xd, xx_vec, dy, tempk, Ny, Nx, atten, norm_coef, attenuation_correction, 
						normalization, temp_ijk, Nx, lo, global_factor, scatter, scatter_coef);

					if (fp) {
						for (uint32_t k = 0; k < Np; k++) {
							ax += (element * osem_apu[temp_ijk + k * Nx]);
						}
						if (ax == 0.)
							ax = epps;
						if (randoms_correction)
							ax += randoms[lo];
						rhs[lo] = ax;
						continue;
					}
					if (local_sino > 0.) {
						for (uint32_t k = 0; k < Np; k++) {
							ax += (element * osem_apu[temp_ijk + k * Nx]);
						}
						if (ax == 0.)
							ax = epps;
						if (randoms_correction)
							ax += randoms[lo];
						const double yax = local_sino / ax;
						for (uint32_t k = 0; k < Np; k++) {
#pragma omp atomic
							rhs[temp_ijk + k * Nx] += (element * yax);
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
			else {
				int32_t tempi = 0, tempj = 0, iu = 0, ju = 0;
				double txu = 0., tyu = 0., tc = 0., tx0 = 0., ty0 = 0.;

				const bool skip = siddon_pre_loop_2D(bx, by, x_diff, y_diff, maxxx, maxyy, dx, dy, Nx, Ny, tempi, tempj, txu, tyu, Np, TYPE,
					detectors.ys, detectors.xs, detectors.yd, detectors.xd, tc, iu, ju, tx0, ty0);

				const double LL = sqrt(x_diff * x_diff + y_diff * y_diff);

				double temp = 0.;
				double tx0_a = tx0, ty0_a = ty0, tc_a = tc;
				uint32_t tempijk = tempk * Nyx + static_cast<uint32_t>(tempj) * Nx + static_cast<uint32_t>(tempi);

				for (uint32_t ii = 0; ii < Np; ii++) {

					if (tx0 < ty0) {
						const double element = (tx0 - tc) * LL;

						temp += element;
						ax += (element * osem_apu[tempijk]);
						if (attenuation_correction)
							jelppi += (element * -atten[tempijk]);

						if (iu > 0)
							tempijk++;
						else
							tempijk--;
						tc = tx0;
						tx0 += txu;


					}
					else {

						const double element = (ty0 - tc) * LL;

						temp += element;
						ax += (element * osem_apu[tempijk]);
						if (attenuation_correction)
							jelppi += (element * -atten[tempijk]);

						if (ju > 0)
							tempijk += Nx;
						else
							tempijk -= Nx;
						tc = ty0;
						ty0 += tyu;
					}
				}

				temp = 1. / temp;
				tx0 = tx0_a;
				ty0 = ty0_a;
				tempijk = tempk * Nyx + static_cast<uint32_t>(tempj) * Nx + static_cast<uint32_t>(tempi);
				tc = tc_a;
				if (attenuation_correction)
					temp *= exp(jelppi);
				if (normalization)
					temp *= norm_coef[lo];
				if (scatter)
					temp *= scatter_coef[lo];
				temp *= global_factor;
				if (fp) {
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
					ax *= temp;
					if (ax == 0.) {
						ax = epps;
					}
					if (randoms_correction)
						ax += randoms[lo];
					const double yax = local_sino / ax;
					for (uint32_t ii = 0; ii < Np; ii++) {

						if (tx0 < ty0) {
							const double element = (tx0 - tc) * LL * temp;


#pragma omp atomic
							rhs[tempijk] += (element * yax);
							if (no_norm == 0) {
#pragma omp atomic
								Summ[tempijk] += element;
							}

							if (iu > 0)
								tempijk++;
							else
								tempijk--;
							tc = tx0;
							tx0 += txu;


						}
						else {

							const double element = (ty0 - tc) * LL * temp;

#pragma omp atomic
							rhs[tempijk] += (element * yax);
							if (no_norm == 0) {
#pragma omp atomic
								Summ[tempijk] += element;
							}

							if (ju > 0)
								tempijk += Nx;
							else
								tempijk -= Nx;
							tc = ty0;
							ty0 += tyu;
						}
					}
				}
				else {
					for (uint32_t ii = 0; ii < Np; ii++) {

						if (tx0 < ty0) {
							const double element = (tx0 - tc) * LL * temp;

							if (no_norm == 0) {
#pragma omp atomic
								Summ[tempijk] += element;
							}

							if (iu > 0)
								tempijk++;
							else
								tempijk--;
							tc = tx0;
							tx0 += txu;


						}
						else {

							const double element = (ty0 - tc) * LL * temp;

							if (no_norm == 0) {
#pragma omp atomic
								Summ[tempijk] += element;
							}

							if (ju > 0)
								tempijk += Nx;
							else
								tempijk -= Nx;
							tc = ty0;
							ty0 += tyu;
						}
					}
				}
			}
		}
		else {

			if (fabs(y_diff) < 1e-8) {
				if (detectors.yd <= maxyy && detectors.yd >= by) {

					int32_t tempi = 0, tempk = 0, tempj = 0, iu = 0, ku = 0;
					double txu = 0., tzu = 0., tc = 0., tx0 = 0., tz0 = 0.;

					const bool skip = siddon_pre_loop_2D(bx, bz, x_diff, z_diff, maxxx, bzb, dx, dz, Nx, Nz, tempi, tempk, txu, tzu, Np, TYPE,
						detectors.zs, detectors.xs, detectors.zd, detectors.xd, tc, iu, ku, tx0, tz0);

					const double LL = sqrt((x_diff * x_diff + z_diff * z_diff));
					double apu1;

					for (size_t ii = 0ULL; ii < static_cast<size_t>(Ny); ii++) {
						apu1 = (yy_vec[ii + 1ULL] - detectors.yd);
						if (apu1 > 0.) {
							tempj = static_cast<int32_t>(ii);
							break;
						}
					}

					double temp = 0.;
					double tx0_a = tx0, tz0_a = tz0, tc_a = tc;
					uint32_t tempijk = Nyx * static_cast<uint32_t>(tempk) + static_cast<uint32_t>(tempj) * Nx + static_cast<uint32_t>(tempi);

					for (uint32_t ii = 0; ii < Np; ii++) {

						if (tx0 < tz0) {

							const double element = (tx0 - tc) * LL;

							temp += element;
							ax += (element * osem_apu[tempijk]);
							if (attenuation_correction)
								jelppi += (element * -atten[tempijk]);

							if (iu > 0)
								tempijk++;
							else
								tempijk--;
							tc = tx0;
							tx0 += txu;
						}
						else {

							const double element = (tz0 - tc) * LL;

							temp += element;
							ax += (element * osem_apu[tempijk]);
							if (attenuation_correction)
								jelppi += (element * -atten[tempijk]);

							if (ku > 0)
								tempijk += Nyx;
							else
								tempijk -= Nyx;
							tc = tz0;
							tz0 += tzu;
						}

					}

					temp = 1. / temp;
					tx0 = tx0_a;
					tz0 = tz0_a;
					tc = tc_a;
					tempijk = Nyx * static_cast<uint32_t>(tempk) + static_cast<uint32_t>(tempj) * Nx + static_cast<uint32_t>(tempi);
					if (attenuation_correction)
						temp *= exp(jelppi);
					if (normalization)
						temp *= norm_coef[lo];
					if (scatter)
						temp *= scatter_coef[lo];
					temp *= global_factor;

					if (fp) {
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
						ax *= temp;
						if (ax == 0.) {
							ax = epps;
						}
						if (randoms_correction)
							ax += randoms[lo];
						const double yax = local_sino / ax;
						for (uint32_t ii = 0; ii < Np; ii++) {

							if (tx0 < tz0) {

								const double element = (tx0 - tc) * LL * temp;

#pragma omp atomic
								rhs[tempijk] += (element * yax);
								if (no_norm == 0) {
#pragma omp atomic
									Summ[tempijk] += element;
								}

								if (iu > 0)
									tempijk++;
								else
									tempijk--;
								tc = tx0;
								tx0 += txu;

							}
							else {

								const double element = (tz0 - tc) * LL * temp;

#pragma omp atomic
								rhs[tempijk] += (element * yax);
								if (no_norm == 0) {
#pragma omp atomic
									Summ[tempijk] += element;
								}

								if (ku > 0)
									tempijk += Nyx;
								else
									tempijk -= Nyx;
								tc = tz0;
								tz0 += tzu;

							}

						}
					}
					else {
						for (uint32_t ii = 0; ii < Np; ii++) {

							if (tx0 < tz0) {

								const double element = (tx0 - tc) * LL * temp;

								if (no_norm == 0) {
#pragma omp atomic
									Summ[tempijk] += element;
								}

								if (iu > 0)
									tempijk++;
								else
									tempijk--;
								tc = tx0;
								tx0 += txu;

							}
							else {

								const double element = (tz0 - tc) * LL * temp;

								if (no_norm == 0) {
#pragma omp atomic
									Summ[tempijk] += element;
								}

								if (ku > 0)
									tempijk += Nyx;
								else
									tempijk -= Nyx;
								tc = tz0;
								tz0 += tzu;

							}

						}

					}

				}
			}
			else if (fabs(x_diff) < 1e-8) {
				if (detectors.xd <= maxxx && detectors.xd >= bx) {

					int32_t tempi = 0, tempk = 0, tempj = 0, ju = 0, ku = 0;
					double tyu = 0., tzu = 0., tc = 0., ty0 = 0., tz0 = 0.;
					const bool skip = siddon_pre_loop_2D(by, bz, y_diff, z_diff, maxyy, bzb, dy, dz, Ny, Nz, tempj, tempk, tyu, tzu, Np, TYPE,
						detectors.zs, detectors.ys, detectors.zd, detectors.yd, tc, ju, ku, ty0, tz0);

					const double LL = sqrt((y_diff * y_diff + z_diff * z_diff));
					double apu1;

					double temp = 0.;

					for (size_t ii = 0ULL; ii < static_cast<size_t>(Nx); ii++) {
						apu1 = (xx_vec[ii + 1ULL] - detectors.xd);
						if (apu1 > 0.) {
							tempi = static_cast<int32_t>(ii);
							break;
						}
					}

					double ty0_a = ty0, tz0_a = tz0, tc_a = tc;
					uint32_t tempijk = Nyx * static_cast<uint32_t>(tempk) + static_cast<uint32_t>(tempj) * Nx + static_cast<uint32_t>(tempi);

					for (uint32_t ii = 0; ii < Np; ii++) {

						if (ty0 < tz0) {

							const double element = (ty0 - tc) * LL;

							temp += element;
							ax += (element * osem_apu[tempijk]);
							if (attenuation_correction)
								jelppi += (element * -atten[tempijk]);

							if (ju > 0)
								tempijk += Nx;
							else
								tempijk -= Nx;
							tc = ty0;
							ty0 += tyu;
						}
						else {

							const double element = (tz0 - tc) * LL;

							temp += element;
							ax += (element * osem_apu[tempijk]);
							if (attenuation_correction)
								jelppi += (element * -atten[tempijk]);

							if (ku > 0)
								tempijk += Nyx;
							else
								tempijk -= Nyx;
							tc = tz0;
							tz0 += tzu;
						}
					}

					temp = 1. / temp;
					ty0 = ty0_a;
					tz0 = tz0_a;
					tc = tc_a;
					tempijk = Nyx * static_cast<uint32_t>(tempk) + static_cast<uint32_t>(tempj) * Nx + static_cast<uint32_t>(tempi);
					if (attenuation_correction)
						temp *= exp(jelppi);
					if (normalization)
						temp *= norm_coef[lo];
					if (scatter)
						temp *= scatter_coef[lo];
					temp *= global_factor;

					if (fp) {
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
						ax *= temp;
						if (ax == 0.) {
							ax = epps;
						}
						if (randoms_correction)
							ax += randoms[lo];
						const double yax = local_sino / ax;
						for (uint32_t ii = 0; ii < Np; ii++) {

							if (ty0 < tz0) {

								const double element = (ty0 - tc) * LL * temp;

#pragma omp atomic
								rhs[tempijk] += (element * yax);
								if (no_norm == 0) {
#pragma omp atomic
									Summ[tempijk] += element;
								}

								if (ju > 0)
									tempijk += Nx;
								else
									tempijk -= Nx;
								tc = ty0;
								ty0 += tyu;
							}
							else {

								const double element = (tz0 - tc) * LL * temp;

#pragma omp atomic
								rhs[tempijk] += (element * yax);

								if (no_norm == 0) {
#pragma omp atomic
									Summ[tempijk] += element;
								}

								if (ku > 0)
									tempijk += Nyx;
								else
									tempijk -= Nyx;
								tc = tz0;
								tz0 += tzu;
							}
						}
					}
					else {
						for (uint32_t ii = 0; ii < Np; ii++) {

							if (ty0 < tz0) {

								const double element = (ty0 - tc) * LL * temp;
								if (no_norm == 0) {
#pragma omp atomic
									Summ[tempijk] += element;
								}

								if (ju > 0)
									tempijk += Nx;
								else
									tempijk -= Nx;
								tc = ty0;
								ty0 += tyu;

							}
							else {

								const double element = (tz0 - tc) * LL * temp;

								if (no_norm == 0) {
#pragma omp atomic
									Summ[tempijk] += element;
								}

								if (ku > 0)
									tempijk += Nyx;
								else
									tempijk -= Nyx;
								tc = tz0;
								tz0 += tzu;

							}
						}
					}
				}
			}
			else {

				int32_t tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;
				double txu = 0., tyu = 0., tzu = 0., tc = 0., tx0 = 0., ty0 = 0., tz0 = 0.;
				const bool skip = siddon_pre_loop_3D(bx, by, bz, x_diff, y_diff, z_diff, maxxx, maxyy, bzb, dx, dy, dz, Nx, Ny, Nz, tempi, tempj, tempk, tyu, 
					txu, tzu, Np, TYPE, detectors, tc, iu, ju, ku, tx0, ty0, tz0);

				const double LL = sqrt(x_diff * x_diff + z_diff * z_diff + y_diff * y_diff);

				double temp = 0.;

				double ty0_a = ty0, tz0_a = tz0, tc_a = tc, tx0_a = tx0;
				uint32_t tempijk = Nyx * static_cast<uint32_t>(tempk) + static_cast<uint32_t>(tempj) * Nx + static_cast<uint32_t>(tempi);

				for (uint32_t ii = 0; ii < Np; ii++) {
					if (tz0 < ty0 && tz0 < tx0) {

						const double element = (tz0 - tc) * LL;

						temp += element;
						ax += (element * osem_apu[tempijk]);
						if (attenuation_correction)
							jelppi += (element * -atten[tempijk]);

						if (ku > 0)
							tempijk += Nyx;
						else
							tempijk -= Nyx;
						tc = tz0;
						tz0 += tzu;
					}
					else if (ty0 < tx0) {
						const double element = (ty0 - tc) * LL;

						temp += element;
						ax += (element * osem_apu[tempijk]);
						if (attenuation_correction)
							jelppi += (element * -atten[tempijk]);

						if (ju > 0)
							tempijk += Nx;
						else
							tempijk -= Nx;
						tc = ty0;
						ty0 += tyu;
					}
					else {
						const double element = (tx0 - tc) * LL;

						temp += element;
						ax += (element * osem_apu[tempijk]);
						if (attenuation_correction)
							jelppi += (element * -atten[tempijk]);

						if (iu > 0)
							tempijk++;
						else
							tempijk--;
						tc = tx0;
						tx0 += txu;
					}

				}

				temp = 1. / temp;
				ty0 = ty0_a;
				tx0 = tx0_a;
				tz0 = tz0_a;
				tc = tc_a;
				tempijk = Nyx * static_cast<uint32_t>(tempk) + static_cast<uint32_t>(tempj) * Nx + static_cast<uint32_t>(tempi);
				if (attenuation_correction)
					temp *= exp(jelppi);
				if (normalization)
					temp *= norm_coef[lo];
				if (scatter)
					temp *= scatter_coef[lo];
				temp *= global_factor;


				if (fp) {
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
					ax *= temp;
					if (ax == 0.) {
						ax = epps;
					}
					if (randoms_correction)
						ax += randoms[lo];
					const double yax = local_sino / ax;
					for (uint32_t ii = 0; ii < Np; ii++) {
						if (tz0 < ty0 && tz0 < tx0) {

							const double element = (tz0 - tc) * LL * temp;

#pragma omp atomic
							rhs[tempijk] += (element * yax);
							if (no_norm == 0) {
#pragma omp atomic
								Summ[tempijk] += element;
							}

							if (ku > 0)
								tempijk += Nyx;
							else
								tempijk -= Nyx;
							tc = tz0;
							tz0 += tzu;
						}
						else if (ty0 < tx0) {
							const double element = (ty0 - tc) * LL * temp;

#pragma omp atomic
							rhs[tempijk] += (element * yax);
							if (no_norm == 0) {
#pragma omp atomic
								Summ[tempijk] += element;
							}

							if (ju > 0)
								tempijk += Nx;
							else
								tempijk -= Nx;
							tc = ty0;
							ty0 += tyu;
						}
						else {
							const double element = (tx0 - tc) * LL * temp;

#pragma omp atomic
							rhs[tempijk] += (element * yax);
							if (no_norm == 0) {
#pragma omp atomic
								Summ[tempijk] += element;
							}

							if (iu > 0)
								tempijk++;
							else
								tempijk--;
							tc = tx0;
							tx0 += txu;
						}

					}

				}
				else {
					for (uint32_t ii = 0; ii < Np; ii++) {
						if (tz0 < ty0 && tz0 < tx0) {

							const double element = (tz0 - tc) * LL * temp;

							if (no_norm == 0) {
#pragma omp atomic
								Summ[tempijk] += element;
							}

							if (ku > 0)
								tempijk += Nyx;
							else
								tempijk -= Nyx;
							tc = tz0;
							tz0 += tzu;
						}
						else if (ty0 < tx0) {
							const double element = (ty0 - tc) * LL * temp;

							if (no_norm == 0) {
#pragma omp atomic
								Summ[tempijk] += element;
							}

							if (ju > 0)
								tempijk += Nx;
							else
								tempijk -= Nx;
							tc = ty0;
							ty0 += tyu;
						}
						else {
							const double element = (tx0 - tc) * LL * temp;

							if (no_norm == 0) {
#pragma omp atomic
								Summ[tempijk] += element;
							}

							if (iu > 0)
								tempijk++;
							else
								tempijk--;
							tc = tx0;
							tx0 += txu;
						}
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
	const bool raw, const double crystal_size_xy, const double* x_center, const double* y_center, const double* z_center, const double crystal_size_z,
	const bool no_norm, const uint32_t dec_v, const double global_factor, const bool fp, const bool scatter, const double* scatter_coef) {

	setThreads();

	const uint32_t Nyx = Ny * Nx;

	const double bzb = bz + static_cast<double>(Nz) * dz;

	size_t idx = 0ULL;
	vector<double> elements;
	vector<uint32_t> v_indices;

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
		const double x_diff = (detectors.xd - detectors.xs);
		const double y_diff = (detectors.yd - detectors.ys);
		const double z_diff = (detectors.zd - detectors.zs);

		// Load the number of voxels the LOR traverses (precomputed)
		uint32_t Np = static_cast<uint32_t>(lor1[lo]);
		double ax = 0., jelppi = 0., LL;
		uint8_t xyz = 0u;
		bool RHS = false, SUMMA = false;
		uint32_t ind = 0u;

		if (crystal_size_z == 0.) {
			kerroin = norm(x_diff, y_diff, z_diff) * crystal_size_xy;
		}
		else {
			kerroin = norm(x_diff, y_diff, z_diff) * crystal_size_z;
		}

		if (fabs(z_diff) < 1e-8) {

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
							ind, rhs, 0, lo, PRECOMPUTE, global_factor, scatter, scatter_coef);
						if (fp) {
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
							nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
							orth_distance_rhs_perpendicular_mfree(y_center, x_center[0], z_center, kerroin, temp, ax, by, detectors.yd, dy, Ny, Nx, tempk, Ny,
								1u, no_norm, rhs, Summ, true, false, detectors, y_diff, x_diff, z_diff, store_elements, store_indices, tid, ind, 0, 0);
						}
						else {
							orth_distance_rhs_perpendicular_mfree(y_center, x_center[0], z_center, kerroin, temp, ax, by, detectors.yd, dy, Ny, Nx, tempk, Ny,
								1u, no_norm, rhs, Summ, false, true, detectors, y_diff, x_diff, z_diff, store_elements, store_indices, tid, ind, 0, 0);
						}
					}
					else {
						double temp = 0.;
						orth_distance_denominator_perpendicular_mfree_3D(y_center, x_center[0], z_center, temp, attenuation_correction, normalization, ax,
							by, detectors.yd, dy, Ny, Nx, tempk, atten, norm_coef, local_sino, Ny, 1u, osem_apu, detectors, y_diff, x_diff, z_diff, kerroin, Nyx, Nz, store_elements, store_indices, tid,
							ind, rhs, 0, lo, PRECOMPUTE, global_factor, scatter, scatter_coef);
						if (local_sino > 0.) {
							nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
							orth_distance_rhs_perpendicular_mfree(y_center, x_center[0], z_center, kerroin, temp, ax, by, detectors.yd, dy, Ny, Nx, tempk, Ny,
								1u, no_norm, rhs, Summ, true, false, detectors, y_diff, x_diff, z_diff, store_elements, store_indices, tid, ind, 0, 0);
						}
						else {
							orth_distance_rhs_perpendicular_mfree(y_center, x_center[0], z_center, kerroin, temp, ax, by, detectors.yd, dy, Ny, Nx, tempk, Ny,
								1u, no_norm, rhs, Summ, false, true, detectors, y_diff, x_diff, z_diff, store_elements, store_indices, tid, ind, 0, 0);
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
							ind, rhs, 0, lo, PRECOMPUTE, global_factor, scatter, scatter_coef);
						if (fp) {
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
							nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
							orth_distance_rhs_perpendicular_mfree(x_center, y_center[0], z_center, kerroin, temp, ax, bx, detectors.xd, dx, Nx, Ny, tempk,
								1u, Nx, no_norm, rhs, Summ, true, false, detectors, x_diff, y_diff, z_diff, store_elements, store_indices, tid, ind, 0, 0);
						}
						else {
							orth_distance_rhs_perpendicular_mfree(x_center, y_center[0], z_center, kerroin, temp, ax, bx, detectors.xd, dx, Nx, Ny, tempk,
								1u, Nx, no_norm, rhs, Summ, false, true, detectors, x_diff, y_diff, z_diff, store_elements, store_indices, tid, ind, 0, 0);
						}
					}
					else {
						double temp = 0.;
						orth_distance_denominator_perpendicular_mfree_3D(x_center, y_center[0], z_center, temp, attenuation_correction, normalization, ax,
							bx, detectors.xd, dx, Nx, Ny, tempk, atten, norm_coef, local_sino, 1u, Nx, osem_apu, detectors, x_diff, y_diff, z_diff, kerroin, Nyx, Nz, store_elements, store_indices, tid,
							ind, rhs, 0, lo, PRECOMPUTE, global_factor, scatter, scatter_coef);
						if (local_sino > 0.) {
							nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
							orth_distance_rhs_perpendicular_mfree(x_center, y_center[0], z_center, kerroin, temp, ax, bx, detectors.xd, dx, Nx, Ny, tempk,
								1u, Nx, no_norm, rhs, Summ, true, false, detectors, x_diff, y_diff, z_diff, store_elements, store_indices, tid, ind, 0, 0);
						}
						else {
							orth_distance_rhs_perpendicular_mfree(x_center, y_center[0], z_center, kerroin, temp, ax, bx, detectors.xd, dx, Nx, Ny, tempk,
								1u, Nx, no_norm, rhs, Summ, false, true, detectors, x_diff, y_diff, z_diff, store_elements, store_indices, tid, ind, 0, 0);
						}
					}
				}
			}
			else {
				int32_t tempi = 0, tempj = 0, iu = 0, ju = 0;
				double txu = 0., tyu = 0., tc = 0., tx0 = 0., ty0 = 0.;
				const bool skip = siddon_pre_loop_2D(bx, by, x_diff, y_diff, maxxx, maxyy, dx, dy, Nx, Ny, tempi, tempj, txu, tyu, Np, TYPE,
					detectors.ys, detectors.xs, detectors.yd, detectors.xd, tc, iu, ju, tx0, ty0);
				if (attenuation_correction)
					LL = sqrt(x_diff * x_diff + y_diff * y_diff);
				double temp = 0.;
				int alku, loppu;
				if (crystal_size_z == 0.) {
					alku = tempk + 1;
					loppu = tempk;
				}
				else {
					alku = Nz;
					loppu = 0;
				}
				orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
					tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
					PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind);

				if (attenuation_correction) {
					for (uint32_t ii = 0u; ii < Np; ii++) {
						if (tx0 < ty0) {
							if (attenuation_correction)
								compute_attenuation(tc, jelppi, LL, tx0, tempi, tempj, tempk, Nx, Nyx, atten);

							tempi += iu;
							tx0 += txu;
						}
						else {
							if (attenuation_correction)
								compute_attenuation(tc, jelppi, LL, ty0, tempi, tempj, tempk, Nx, Nyx, atten);
							tempj += ju;
							ty0 += tyu;
						}
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

				if (fp) {
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
					nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
					RHS = true;
				}
				else
					SUMMA = true;

				orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
					tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
					PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind);
			}
		}
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
					orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
						tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
						PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Ny, Nx, alku, iu, 0, loppu, store_elements, store_indices, tid, ind);
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
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
									tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
									PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Ny, Nx, alku, iu, 0, loppu, store_elements, store_indices, tid, ind);
							}
						}
					}
					if (xyz < 3 && crystal_size_z > 0.) {
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
							orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
								tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
								PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Ny, Nx, alku, iu, 0, loppu, store_elements, store_indices, tid, ind);
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
					if (fp) {
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
						nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
						RHS = true;
					}
					else
						SUMMA = true;
					orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
						tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
						PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Ny, Nx, alku, iu, 0, loppu, store_elements, store_indices, tid, ind);
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
					orth_distance_3D_full(tempj, Ny, Nz, x_diff, y_diff, z_diff, x_center, y_center, z_center, temp, Nx,
						tempi, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
						PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Nx, 1u, alku, ju, 0, loppu, store_elements, store_indices, tid, ind);

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
								orth_distance_3D_full(tempj, Ny, Nz, x_diff, y_diff, z_diff, x_center, y_center, z_center, temp, Nx,
									tempi, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
									PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Nx, 1u, alku, ju, 0, loppu, store_elements, store_indices, tid, ind);
							}
						}
					}
					if (xyz < 3 && crystal_size_z > 0.) {
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
							orth_distance_3D_full(tempj, Ny, Nz, x_diff, y_diff, z_diff, x_center, y_center, z_center, temp, Nx,
								tempi, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
								PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Nx, 1u, alku, ju, 0, loppu, store_elements, store_indices, tid, ind);
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


					if (fp) {
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
						nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
						RHS = true;
					}
					else
						SUMMA = true;

					orth_distance_3D_full(tempj, Ny, Nz, x_diff, y_diff, z_diff, x_center, y_center, z_center, temp, Nx,
						tempi, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
						PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Nx, 1u, alku, ju, 0, loppu, store_elements, store_indices, tid, ind);
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
				orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
					tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
					PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind);

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
							orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
								tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
								PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind);
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
				if (xyz < 3 && crystal_size_z > 0.) {
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
						orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
							tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
							PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind);
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

				if (fp) {
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
					nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
					RHS = true;
				}
				else
					SUMMA = true;
				orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
					tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
					PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind);

			}
		}
	}
}

void sequential_volume_siddon(const int64_t loop_var_par, const uint32_t size_x, const double zmax, double* Summ, double* rhs, const double maxyy,
	const double maxxx, const vector<double>& xx_vec, const double dy, const vector<double>& yy_vec, const double* atten, const double* norm_coef,
	const double* randoms, const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny,
	const uint32_t Nz, const double dx, const double dz, const double bx, const double by, const double bz, const bool attenuation_correction,
	const bool normalization, const bool randoms_correction, const uint16_t* lor1, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos,
	const double epps, const double* Sino, double* osem_apu, const uint16_t* L, const uint32_t* pseudos, const uint32_t pRows, const uint32_t det_per_ring,
	const bool raw, const double Vmax, const double* x_center, const double* y_center, const double* z_center, const double bmin, const double bmax, const double* V,
	const bool no_norm, const uint32_t dec_v, const double global_factor, const bool fp, const bool scatter, const double* scatter_coef) {

	setThreads();

	const uint32_t Nyx = Ny * Nx;

	const double bzb = bz + static_cast<double>(Nz) * dz;

	size_t idx = 0ULL;
	vector<double> elements;
	vector<uint32_t> v_indices;

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
		const double x_diff = (detectors.xd - detectors.xs);
		const double y_diff = (detectors.yd - detectors.ys);
		const double z_diff = (detectors.zd - detectors.zs);

		// Load the number of voxels the LOR traverses (precomputed)
		uint32_t Np = static_cast<uint32_t>(lor1[lo]);
		double ax = 0., jelppi = 0., LL;
		int8_t start = 1;
		uint8_t xyz = 0u;
		uint8_t xyz_w = 0u;
		bool RHS = false, SUMMA = false;
		uint32_t ind = 0u;

		kerroin = norm(x_diff, y_diff, z_diff);

		if (fabs(z_diff) < 1e-8) {

			const uint32_t tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));

			if (fabs(y_diff) < 1e-8) {

				if (detectors.yd <= maxyy && detectors.yd >= by) {
					double temppi = detectors.xs;
					detectors.xs = detectors.ys;
					detectors.ys = temppi;
					double temp = 0.;
					volume_distance_denominator_perpendicular_mfree_3D(y_center, x_center[0], z_center, temp, attenuation_correction, normalization, ax,
						by, detectors.yd, dy, Ny, Nx, tempk, atten, norm_coef, local_sino, Ny, 1u, osem_apu, detectors, y_diff, x_diff, z_diff, kerroin, Nyx, Nz, store_elements, store_indices, tid,
						ind, rhs, 0, lo, PRECOMPUTE, global_factor, bmax, bmin, Vmax, V, scatter, scatter_coef);
					if (fp) {
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
						nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
						orth_distance_rhs_perpendicular_mfree(y_center, x_center[0], z_center, kerroin, temp, ax, by, detectors.yd, dy, Ny, Nx, tempk, Ny,
							1u, no_norm, rhs, Summ, true, false, detectors, y_diff, x_diff, z_diff, store_elements, store_indices, tid, ind, 0, 0);
					}
					else {
						orth_distance_rhs_perpendicular_mfree(y_center, x_center[0], z_center, kerroin, temp, ax, by, detectors.yd, dy, Ny, Nx, tempk, Ny,
							1u, no_norm, rhs, Summ, false, true, detectors, y_diff, x_diff, z_diff, store_elements, store_indices, tid, ind, 0, 0);
					}
				}
			}
			else if (fabs(x_diff) < 1e-8) {

				if (detectors.xd <= maxxx && detectors.xd >= bx) {
					double temp = 0.;
					volume_distance_denominator_perpendicular_mfree_3D(x_center, y_center[0], z_center, temp, attenuation_correction, normalization, ax,
						bx, detectors.xd, dx, Nx, Ny, tempk, atten, norm_coef, local_sino, 1u, Nx, osem_apu, detectors, x_diff, y_diff, z_diff, kerroin, Nyx, Nz, store_elements, store_indices, tid,
						ind, rhs, 0, lo, PRECOMPUTE, global_factor, bmax, bmin, Vmax, V, scatter, scatter_coef);
					if (fp) {
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
						nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
						orth_distance_rhs_perpendicular_mfree(x_center, y_center[0], z_center, kerroin, temp, ax, bx, detectors.xd, dx, Nx, Ny, tempk,
							1u, Nx, no_norm, rhs, Summ, true, false, detectors, x_diff, y_diff, z_diff, store_elements, store_indices, tid, ind, 0, 0);
					}
					else {
						orth_distance_rhs_perpendicular_mfree(x_center, y_center[0], z_center, kerroin, temp, ax, bx, detectors.xd, dx, Nx, Ny, tempk,
							1u, Nx, no_norm, rhs, Summ, false, true, detectors, x_diff, y_diff, z_diff, store_elements, store_indices, tid, ind, 0, 0);
					}
				}
			}
			else {
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
					PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);

				if (attenuation_correction) {
					for (uint32_t ii = 0u; ii < Np; ii++) {
						if (tx0 < ty0) {
							if (attenuation_correction)
								compute_attenuation(tc, jelppi, LL, tx0, tempi, tempj, tempk, Nx, Nyx, atten);
							tempi += iu;
							tx0 += txu;
							xyz = 1u;
						}
						else {
							if (attenuation_correction)
								compute_attenuation(tc, jelppi, LL, ty0, tempi, tempj, tempk, Nx, Nyx, atten);
							tempj += ju;
							ty0 += tyu;
							xyz = 2u;
						}
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

				if (fp) {
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
					nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
					RHS = true;
				}
				else
					SUMMA = true;

				volume_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
					tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
					PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);
			}
		}
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
						PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Ny, Nx, alku, iu, 0, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);

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
									PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Ny, Nx, alku, iu, 0, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);
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
								PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Ny, Nx, alku, iu, 0, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);
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
					if (fp) {
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
						nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
						RHS = true;
					}
					else
						SUMMA = true;
					volume_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
						tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
						PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Ny, Nx, alku, iu, 0, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);
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
						PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Nx, 1u, alku, ju, 0, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);

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
									PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Nx, 1u, alku, ju, 0, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);
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
								PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Nx, 1u, alku, ju, 0, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);
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


					if (fp) {
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
						nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
						RHS = true;
					}
					else
						SUMMA = true;

					volume_distance_3D_full(tempj, Ny, Nz, x_diff, y_diff, z_diff, x_center, y_center, z_center, temp, Nx,
						tempi, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
						PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Nx, 1u, alku, ju, 0, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);
				}
			}
			else {
				int32_t tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;
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
					PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);

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
								PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);
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
							PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);
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

				if (fp) {
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
					nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
					RHS = true;
				}
				else
					SUMMA = true;
				volume_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
					tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
					PRECOMPUTE, DISCARD, rhs, Summ, 0, elements, v_indices, idx, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);
			}
		}
	}
}