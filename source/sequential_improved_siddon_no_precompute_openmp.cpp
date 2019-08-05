/**************************************************************************
* Implements both the improved Siddon's algorithm and Orthogonal Siddon's 
* algorithm for OMEGA.
* Determines which LORs intercept the FOV on-the-fly (slower).
*
* Uses OpenMP for parallellization. If OpenMP is not available, the code
* is serial with no parallellization.
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
#ifdef _OPENMP
#include <omp.h>
#endif

// if 0, then determines whether the LOR intercepts the FOV
constexpr int TYPE = 0;

// Whether to use the OpenMP code or not
constexpr bool OMP = true;

// Using non-OpenMP with either precomputation or without
constexpr bool PRECOMPUTE = false;

// Normalized distances below this are discarded in orthogonal ray tracer
constexpr auto THR = 0.01;

using namespace std;

void sequential_improved_siddon_no_precompute(const size_t loop_var_par, const uint32_t size_x, const double zmax, double* Summ, double* rhs, const double maxyy, const double maxxx,
	const vector<double>& xx_vec, const double dy, const vector<double>& yy_vec, const double* atten, const double* norm_coef, const double* randoms, const double* x, const double* y, const double* z_det,
	const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const double dx, const double dz, const double bx, const double by, const double bz,
	const bool attenuation_correction, const bool normalization, const bool randoms_correction, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos,
	const double epps, const double* Sino, double* osem_apu, const uint16_t* L, const uint32_t* pseudos, const uint32_t pRows, const uint32_t det_per_ring,
	const bool raw, const bool no_norm) {


	const uint32_t Nyx = Ny * Nx;

	const double bzb = bz + static_cast<double>(Nz) * dz;

#pragma omp parallel for
	for (uint32_t lo = 0; lo < loop_var_par; lo++) {

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
		if ((y_diff == 0. && x_diff == 0. && z_diff == 0.) || (y_diff == 0. && x_diff == 0.))
			continue;

		uint32_t Np = 0u;
		uint32_t Np_n = 0u;
		double ax = 0., jelppi = 0.;


		if (fabs(z_diff) < 1e-8) {

			const uint32_t tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));

			if (fabs(y_diff) < 1e-8) {


//#pragma omp critical
				//mexPrintf("lo = %d\n", lo);

				if (detectors.yd <= maxyy && detectors.yd >= by) {
					uint32_t temp_ijk = 0;

					const double element = perpendicular_elements(Ny, detectors.yd, yy_vec, dx, tempk, Nx, Ny, atten, norm_coef, attenuation_correction, normalization, temp_ijk, 1u, lo);

					if (local_sino != 0.) {
						for (uint32_t k = 0; k < Nx; k++) {
							ax += (element * osem_apu[temp_ijk + k]);
						}
						if (ax == 0.)
							ax = epps;
						if (randoms_correction)
							ax += randoms[lo];
						const double yax = local_sino / ax;
						for (uint32_t k = 0; k < Nx; k++) {
#pragma omp atomic
							rhs[temp_ijk + k] += (element * yax);
							if (no_norm == 0) {
#pragma omp atomic
								Summ[temp_ijk + k] += element;
							}
						}
					}
					else {
						for (uint32_t k = 0; k < Nx; k++) {
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

					const double element = perpendicular_elements(1, detectors.xd, xx_vec, dy, tempk, Ny, Nx, atten, norm_coef, attenuation_correction, normalization, temp_ijk, Nx, lo);

					if (local_sino != 0.) {
						for (uint32_t k = 0; k < Ny; k++) {
							ax += (element * osem_apu[temp_ijk + k * Nx]);
						}
						if (ax == 0.)
							ax = epps;
						if (randoms_correction)
							ax += randoms[lo];
						const double yax = local_sino / ax;
						for (uint32_t k = 0; k < Ny; k++) {
#pragma omp atomic
							rhs[temp_ijk + k * Nx] += (element * yax);
							if (no_norm == 0) {
#pragma omp atomic
								Summ[temp_ijk + k * Nx] += element;
							}
						}
					}
					else {
						for (uint32_t k = 0; k < Ny; k++) {
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

				if (tempi < 0 || tempj < 0 || tempi >= static_cast<int32_t>(Nx) || tempj >= static_cast<int32_t>(Ny))
					continue;
				const double LL = sqrt(x_diff * x_diff + y_diff * y_diff);

				double temp = 0.;
				//int32_t tempi_a = tempi, tempj_a = tempj;
				double tx0_a = tx0, ty0_a = ty0, tc_a = tc;
				uint32_t tempijk = tempk * Nyx + static_cast<uint32_t>(tempj) * Nx + static_cast<uint32_t>(tempi);

				for (uint32_t ii = 0; ii < Np; ii++) {

					if (tx0 < ty0) {
						const double element = (tx0 - tc) * LL;

						temp += element;
						ax += (element * osem_apu[tempijk]);
						if (attenuation_correction)
							jelppi += (element * -atten[tempijk]);

						//tempi += iu;
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

						//tempj += ju;
						if (ju > 0)
							tempijk += Nx;
						else
							tempijk -= Nx;
						tc = ty0;
						ty0 += tyu;
					}
					Np_n++;
					if (tempj < 0 || tempi < 0 || tempi >= static_cast<int32_t>(Nx) || tempj >= static_cast<int32_t>(Ny))
						break;
				}

				temp = 1. / temp;
				tx0 = tx0_a;
				ty0 = ty0_a;
				//tempi = tempi_a;
				//tempj = tempj_a;
				tempijk = tempk * Nyx + static_cast<uint32_t>(tempj) * Nx + static_cast<uint32_t>(tempi);
				tc = tc_a;
				if (attenuation_correction)
					temp *= exp(jelppi);
				if (normalization)
					temp *= norm_coef[lo];

				if (local_sino != 0.) {
					if (ax == 0.) {
						ax = epps;
					}
					else {
						ax *= temp;
					}
					if (randoms_correction)
						ax += randoms[lo];
					const double yax = local_sino / ax;
					for (uint32_t ii = 0; ii < Np_n; ii++) {

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
					for (uint32_t ii = 0; ii < Np_n; ii++) {

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

					if (tempi < 0 || tempk < 0 || tempi >= static_cast<int32_t>(Nx) || tempk >= static_cast<int32_t>(Nz))
						continue;
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
					//uint32_t tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
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
						Np_n++;
						if (tempk < 0 || tempi < 0 || tempi >= static_cast<int32_t>(Nx) || tempk >= static_cast<int32_t>(Nz))
							break;
					}

					temp = 1. / temp;
					tx0 = tx0_a;
					tz0 = tz0_a;
					//tempi = tempi_a;
					//tempk = tempk_a;
					tc = tc_a;
					tempijk = Nyx * static_cast<uint32_t>(tempk) + static_cast<uint32_t>(tempj) * Nx + static_cast<uint32_t>(tempi);
					if (attenuation_correction)
						temp *= exp(jelppi);
					if (normalization)
						temp *= norm_coef[lo];


					if (local_sino != 0.) {
						if (ax == 0.)
							ax = epps;
						else
							ax *= temp;
						if (randoms_correction)
							ax += randoms[lo];
						const double yax = local_sino / ax;
						for (uint32_t ii = 0; ii < Np_n; ii++) {

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
						for (uint32_t ii = 0; ii < Np_n; ii++) {

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

					if (tempk < 0 || tempj < 0 || tempk >= static_cast<int32_t>(Nz) || tempj >= static_cast<int32_t>(Ny))
						continue;
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

					//uint32_t tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
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
						Np_n++;
						if (tempj < 0 || tempk < 0 || tempk >= static_cast<int32_t>(Nz) || tempj >= static_cast<int32_t>(Ny))
							break;
					}

					temp = 1. / temp;
					ty0 = ty0_a;
					tz0 = tz0_a;
					//tempj = tempj_a;
					//tempk = tempk_a;
					tc = tc_a;
					tempijk = Nyx * static_cast<uint32_t>(tempk) + static_cast<uint32_t>(tempj) * Nx + static_cast<uint32_t>(tempi);
					if (attenuation_correction)
						temp *= exp(jelppi);
					if (normalization)
						temp *= norm_coef[lo];

					if (local_sino != 0.) {
						if (ax == 0.)
							ax = epps;
						else
							ax *= temp;
						if (randoms_correction)
							ax += randoms[lo];
						const double yax = local_sino / ax;
						for (uint32_t ii = 0; ii < Np_n; ii++) {

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
						for (uint32_t ii = 0; ii < Np_n; ii++) {

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
				const bool skip = siddon_pre_loop_3D(bx, by, bz, x_diff, y_diff, z_diff, maxxx, maxyy, bzb, dx, dy, dz, Nx, Ny, Nz, tempi, tempj, tempk, tyu, txu, tzu,
					Np, TYPE, detectors, tc, iu, ju, ku, tx0, ty0, tz0);

				if (tempi < 0 || tempj < 0 || tempk < 0 || tempi >= static_cast<int32_t>(Nx) || tempj >= static_cast<int32_t>(Ny) || tempk >= static_cast<int32_t>(Nz))
					continue;
				const double LL = sqrt(x_diff * x_diff + z_diff * z_diff + y_diff * y_diff);

				double temp = 0.;

				//uint32_t tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
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
					Np_n++;
					if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= static_cast<int32_t>(Nx) || tempj >= static_cast<int32_t>(Ny) || tempk >= static_cast<int32_t>(Nz))
						break;
				}

				temp = 1. / temp;
				ty0 = ty0_a;
				tx0 = tx0_a;
				tz0 = tz0_a;
				//tempi = tempi_a;
				//tempj = tempj_a;
				//tempk = tempk_a;
				tc = tc_a;
				tempijk = Nyx * static_cast<uint32_t>(tempk) + static_cast<uint32_t>(tempj) * Nx + static_cast<uint32_t>(tempi);
				if (attenuation_correction)
					temp *= exp(jelppi);
				if (normalization)
					temp *= norm_coef[lo];


				if (local_sino != 0.) {
					if (ax == 0.)
						ax = epps;
					else
						ax *= temp;
					if (randoms_correction)
						ax += randoms[lo];
					const double yax = local_sino / ax;
					for (uint32_t ii = 0; ii < Np_n; ii++) {
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
					for (uint32_t ii = 0; ii < Np_n; ii++) {
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

void sequential_orth_siddon_no_precomp(const size_t loop_var_par, const uint32_t size_x, const double zmax, double* Summ, double* rhs, const double maxyy, const double maxxx,
	const vector<double>& xx_vec, const double dy, const vector<double>& yy_vec, const double* atten, const double* norm_coef, const double* randoms, const double* x, const double* y, const double* z_det,
	const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const double dx, const double dz, const double bx, const double by, const double bz,
	const bool attenuation_correction, const bool normalization, const bool randoms_correction, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos,
	const double epps, const double* Sino, double* osem_apu, const uint16_t* L, const uint32_t* pseudos, const uint32_t pRows, const uint32_t det_per_ring,
	const bool raw, const double crystal_size_xy, const double* x_center, const double* y_center, const double* z_center, const double crystal_size_z, const bool no_norm, const int32_t dec_v) {

	const uint32_t Nyx = Ny * Nx;

	const double bzb = bz + static_cast<double>(Nz) * dz;

	const int32_T dec = static_cast<int32_T>(ceil(crystal_size_z / sqrt(dz * dz * 2.))) * dec_v;

	size_t idx = 0ULL;
	vector<double> elements;
	vector<uint32_t> v_indices;

#pragma omp parallel for
	for (uint32_t lo = 0u; lo < loop_var_par; lo++) {


		const double local_sino = Sino[lo];
		if (no_norm && local_sino == 0.)
			continue;

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
		if ((y_diff == 0. && x_diff == 0. && z_diff == 0.) || (y_diff == 0. && x_diff == 0.))
			continue;

		double ax = 0., jelppi = 0., LL;
		uint32_t Np = 0u;
		uint32_t Np_n = 0u;
		uint8_t xyz = 0u;

		if (crystal_size_z == 0.) {
			kerroin = detectors.xd * detectors.ys - detectors.yd * detectors.xs;
			length_ = sqrt(y_diff * y_diff + x_diff * x_diff) * crystal_size_xy;
		}
		else {
			kerroin = norm(x_diff, y_diff, z_diff) * crystal_size_z;
		}

		if (fabs(z_diff) < 1e-8) {

			const uint32_t tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));

			if (fabs(y_diff) < 1e-8) {

				if (detectors.yd <= maxyy && detectors.yd >= by) {
					if (crystal_size_z == 0.) {
						double temp = 0.;
						orth_distance_denominator_perpendicular_mfree(-x_diff, y_center, kerroin, length_, temp, attenuation_correction, ax,
							by, detectors.yd, dy, Ny, Nx, tempk, atten, local_sino, Ny, 1u, osem_apu);
						if (local_sino != 0.) {
							nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
							orth_distance_rhs_perpendicular_mfree(-x_diff, y_center, kerroin, length_, temp, ax, by, detectors.yd, dy, Ny, Nx, tempk, Ny, 1u, no_norm, 
								rhs, Summ);
						}
						else {
							orth_distance_summ_perpendicular_mfree(-x_diff, y_center, kerroin, length_, temp, ax, by, detectors.yd, dy, Ny, Nx, tempk, Ny, 1u, Summ);
						}
					}
					else {
						double temppi = detectors.xs;
						detectors.xs = detectors.ys;
						detectors.ys = temppi;
						double temp = 0.;
						orth_distance_denominator_perpendicular_mfree_3D(y_center, x_center[0], z_center, temp, attenuation_correction, ax,
							by, detectors.yd, dy, Ny, Nx, tempk, atten, local_sino, Ny, 1u, osem_apu, detectors, y_diff, x_diff, z_diff, kerroin, Nyx, Nz);
						if (local_sino != 0.) {
							nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
							orth_distance_rhs_perpendicular_mfree_3D(y_center, x_center[0], z_center, temp, ax, by, detectors.yd, dy, Ny, Nx, tempk, Ny, 1u, no_norm,
								rhs, Summ, detectors, y_diff, x_diff, z_diff, kerroin, Nyx, Nz);
						}
						else {
							orth_distance_summ_perpendicular_mfree_3D(y_center, x_center[0], z_center, temp, ax, by, detectors.yd, dy, Ny, Nx, tempk, Ny, 1u, 
								Summ, detectors, y_diff, x_diff, z_diff, kerroin, Nyx, Nz);
						}
					}
				}
			}
			else if (fabs(x_diff) < 1e-8) {

				if (detectors.xd <= maxxx && detectors.xd >= bx) {
					if (crystal_size_z == 0.) {
						double temp = 0.;
						orth_distance_denominator_perpendicular_mfree(y_diff, x_center, kerroin, length_, temp, attenuation_correction, ax,
							bx, detectors.xd, dx, Nx, Ny, tempk, atten, local_sino, 1u, Nx, osem_apu);
						if (local_sino != 0.) {
							nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
							orth_distance_rhs_perpendicular_mfree(y_diff, x_center, kerroin, length_, temp, ax, bx, detectors.xd, dx, Nx, Ny, tempk, 1u, Nx, no_norm,
								rhs, Summ);
						}
						else {
							orth_distance_summ_perpendicular_mfree(y_diff, x_center, kerroin, length_, temp, ax, bx, detectors.xd, dx, Nx, Ny, tempk, 1u, Nx, Summ);
						}
					}
					else {
						double temp = 0.;
						orth_distance_denominator_perpendicular_mfree_3D(x_center, y_center[0], z_center, temp, attenuation_correction, ax,
							bx, detectors.xd, dx, Nx, Ny, tempk, atten, local_sino, 1u, Nx, osem_apu, detectors, x_diff, y_diff, z_diff, kerroin, Nyx, Nz);
						if (local_sino != 0.) {
							nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
							orth_distance_rhs_perpendicular_mfree_3D(x_center, y_center[0], z_center, temp, ax, bx, detectors.xd, dx, Nx, Ny, tempk, 1u, Nx, no_norm,
								rhs, Summ, detectors, x_diff, y_diff, z_diff, kerroin, Nyx, Nz);
						}
						else {
							orth_distance_summ_perpendicular_mfree_3D(x_center, y_center[0], z_center, temp, ax, bx, detectors.xd, dx, Nx, Ny, tempk, 1u, Nx, 
								Summ, detectors, x_diff, y_diff, z_diff, kerroin, Nyx, Nz);
						}
					}
				}
			}
			else {
				int32_t tempi = 0, tempj = 0, iu = 0, ju = 0;
				double txu = 0., tyu = 0., tc = 0., tx0 = 0., ty0 = 0.;

				const bool skip = siddon_pre_loop_2D(bx, by, x_diff, y_diff, maxxx, maxyy, dx, dy, Nx, Ny, tempi, tempj, txu, tyu, Np, TYPE,
					detectors.ys, detectors.xs, detectors.yd, detectors.xd, tc, iu, ju, tx0, ty0);

				if (skip)
					continue;

				if (attenuation_correction)
					LL = sqrt(x_diff * x_diff + y_diff * y_diff);
				double temp = 0.;
				int32_t tempi_a = tempi, tempj_a = tempj;
				double tx0_a = tx0, ty0_a = ty0;
				uint32_t tempijk;
				if (crystal_size_z == 0.)
					tempijk = Nyx * tempk + static_cast<uint32_t>(tempj) * Nx;
				else
					tempijk = static_cast<uint32_t>(tempj) * Nx;

				for (uint32_t ii = 0u; ii < Np; ii++) {


					if (tx0 < ty0) {
						if (attenuation_correction)
							compute_attenuation(tc, jelppi, LL, tx0, tempi, tempj, tempk, Nx, Nyx, atten);
						if (ii == Np - 1u) {
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
									tempj, jelppi, local_sino, ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
							else {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
									tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
						}
						else {
							tempi += iu;
							tx0 += txu;
						}
						xyz = 1u;
					}
					else {
						if (attenuation_correction)
							compute_attenuation(tc, jelppi, LL, ty0, tempi, tempj, tempk, Nx, Nyx, atten);
						if (crystal_size_z == 0.) {
							orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
								tempj, jelppi, local_sino, ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
						}
						else {
							orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
								tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
						}

						if (ju > 0) {
							tempijk += Nx;
						}
						else {
							tempijk -= Nx;
						}
						tempj += ju;
						ty0 += tyu;
						xyz = 2u;
					}
					Np_n++;
					if (tempj < 0 || tempi < 0 || tempi >= static_cast<int32_t>(Nx) || tempj >= static_cast<int32_t>(Ny)) {
						if (xyz == 1u && ii != Np - 1u) {
							tempi -= iu;
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
									tempj, jelppi, local_sino, ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
							else {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
									tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
						}
						break;
					}
				}

				temp = 1. / temp;
				tx0 = tx0_a;
				ty0 = ty0_a;
				tempi = tempi_a;
				tempj = tempj_a;
				if (crystal_size_z == 0.)
					tempijk = Nyx * tempk + static_cast<uint32_t>(tempj) * Nx;
				else
					tempijk = static_cast<uint32_t>(tempj)* Nx;
				if (attenuation_correction)
					temp *= exp(jelppi);
				if (normalization)
					temp *= norm_coef[lo];

				if (local_sino != 0.) {
					nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
					for (uint32_t ii = 0u; ii < Np_n; ii++) {
						if (tx0 < ty0) {
							if (ii == Np_n - 1u) {
								if (crystal_size_z == 0.) {
									orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
										tempj, jelppi, local_sino, ax, osem_apu, no_norm, true, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								}
								else {
									orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
										tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, true, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								}
							}
							else {
								tempi += iu;
								tx0 += txu;
							}
						}
						else {
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
									tempj, jelppi, local_sino, ax, osem_apu, no_norm, true, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
							else {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
									tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, true, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
							if (ju > 0) {
								tempijk += Nx;
							}
							else {
								tempijk -= Nx;
							}
							tempj += ju;
							ty0 += tyu;
						}
					}
				}
				else {
					for (uint32_t ii = 0u; ii < Np_n; ii++) {
						if (tx0 < ty0) {
							if (ii == Np_n - 1u) {
								if (crystal_size_z == 0.) {
									orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
										tempj, jelppi, local_sino, ax, osem_apu, no_norm, false, true, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								}
								else {
									orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
										tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, true, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								}
							}
							else {
								tempi += iu;
								tx0 += txu;
							}
						}
						else {
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
									tempj, jelppi, local_sino, ax, osem_apu, no_norm, false, true, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
							else {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
									tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, true, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
							if (ju > 0) {
								tempijk += Nx;
							}
							else {
								tempijk -= Nx;
							}
							tempj += ju;
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

					if (skip)
						continue;

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
					int32_t tempi_a = tempi, tempk_a = tempk;
					double tx0_a = tx0, tz0_a = tz0;
					uint32_t tempijk;
					if (crystal_size_z == 0.)
						tempijk = Nyx * static_cast<uint32_t>(tempk) + static_cast<uint32_t>(tempi);
					else {
						tempijk = static_cast<uint32_t>(tempi);
						const double temp_x = detectors.xs;
						detectors.xs = detectors.ys;
						detectors.ys = temp_x;
					}

					for (uint32_t ii = 0u; ii < Np; ii++) {
						if (tx0 < tz0) {
							if (attenuation_correction)
								compute_attenuation(tc, jelppi, LL, tx0, tempi, tempj, tempk, Nx, Nyx, atten);
							if (crystal_size_z == 0.) {
								orth_distance_full(tempj, Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, temp, tempijk, Nx,
									tempi, jelppi, local_sino, ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
							else {
								orth_distance_3D_full(tempj, Ny, Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, temp, tempijk, Nx,
									tempi, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
							if (iu > 0) {
								tempijk++;
							}
							else {
								tempijk--;
							}
							tempi += iu;
							tx0 += txu;
							xyz = 1u;
						}
						else {
							if (attenuation_correction)
								compute_attenuation(tc, jelppi, LL, tz0, tempi, tempj, tempk, Nx, Nyx, atten);
							if (crystal_size_z == 0.) {
								orth_distance_full(tempj, Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, temp, tempijk, Nx,
									tempi, jelppi, local_sino, ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								if (ku > 0)
									tempijk += Nyx;
								else
									tempijk -= Nyx;
							}
							else if (ii == Np - 1u) {
								orth_distance_3D_full(tempj, Ny, Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, temp, tempijk, Nx,
									tempi, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
							tempk += ku;
							tz0 += tzu;
							xyz = 3u;
						}
						Np_n++;
						if (tempk < 0 || tempi < 0 || tempi >= static_cast<int32_t>(Nx) || tempk >= static_cast<int32_t>(Nz)) {
							if (crystal_size_z != 0.f && xyz == 3u && ii != Np - 1u) {
								tempk -= ku;
								orth_distance_3D_full(tempj, Ny, Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, temp, tempijk, Nx,
									tempi, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
							break;
						}
					}
					temp = 1. / temp;
					tx0 = tx0_a;
					tz0 = tz0_a;
					tempi = tempi_a;
					tempk = tempk_a;
					if (crystal_size_z == 0.)
						tempijk = Nyx * static_cast<uint32_t>(tempk) + static_cast<uint32_t>(tempi);
					else
						tempijk = static_cast<uint32_t>(tempi);
					if (attenuation_correction)
						temp *= exp(jelppi);
					if (normalization)
						temp *= norm_coef[lo];
					if (local_sino != 0.) {
						nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
						for (uint32_t ii = 0u; ii < Np_n; ii++) {

							if (tx0 < tz0) {
								if (crystal_size_z == 0.) {
									orth_distance_full(tempj, Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, temp, tempijk, Nx,
										tempi, jelppi, local_sino, ax, osem_apu, no_norm, true, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								}
								else {
									orth_distance_3D_full(tempj, Ny, Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, temp, tempijk, Nx,
										tempi, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, true, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								}

								if (iu > 0) {
									tempijk++;
								}
								else {
									tempijk--;
								}
								tempi += iu;
								tx0 += txu;
							}
							else {
								if (crystal_size_z == 0.) {
									orth_distance_full(tempj, Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, temp, tempijk, Nx,
										tempi, jelppi, local_sino, ax, osem_apu, no_norm, true, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);

									if (ku > 0)
										tempijk += Nyx;
									else
										tempijk -= Nyx;
								}
								else if (ii == Np_n - 1u) {
									orth_distance_3D_full(tempj, Ny, Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, temp, tempijk, Nx,
										tempi, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, true, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								}
								tempk += ku;
								tz0 += tzu;
							}
						}
					}
					else {
						for (uint32_t ii = 0u; ii < Np_n; ii++) {

							if (tx0 < tz0) {
								if (crystal_size_z == 0.) {
									orth_distance_full(tempj, Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, temp, tempijk, Nx,
										tempi, jelppi, local_sino, ax, osem_apu, no_norm, false, true, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								}
								else {
									orth_distance_3D_full(tempj, Ny, Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, temp, tempijk, Nx,
										tempi, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, true, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								}

								if (iu > 0) {
									tempijk++;
								}
								else {
									tempijk--;
								}
								tempi += iu;
								tx0 += txu;
							}
							else {
								if (crystal_size_z == 0.) {
									orth_distance_full(tempj, Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, temp, tempijk, Nx,
										tempi, jelppi, local_sino, ax, osem_apu, no_norm, false, true, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);

									if (ku > 0)
										tempijk += Nyx;
									else
										tempijk -= Nyx;
								}
								else if (ii == Np_n - 1u) {
									orth_distance_3D_full(tempj, Ny, Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, temp, tempijk, Nx,
										tempi, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, true, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								}
								tempk += ku;
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

					if (skip)
						continue;

					double apu1;
					double temp = 0.;

					if (attenuation_correction)
						LL = sqrt(z_diff * z_diff + y_diff * y_diff);
					for (size_t ii = 0ULL; ii < static_cast<size_t>(Nx); ii++) {
						apu1 = (xx_vec[ii + 1ULL] - detectors.xd);
						if (apu1 > 0.) {
							tempi = static_cast<int32_t>(ii);
							break;
						}
					}

					int32_t tempj_a = tempj, tempk_a = tempk;
					double ty0_a = ty0, tz0_a = tz0;
					uint32_t tempijk;
					if (crystal_size_z == 0.)
						tempijk = Nyx * static_cast<uint32_t>(tempk) + static_cast<uint32_t>(tempj) * Nx;
					else
						tempijk = static_cast<uint32_t>(tempj) * Nx;

					for (uint32_t ii = 0u; ii < Np; ii++) {

						if (ty0 < tz0) {
							if (attenuation_correction)
								compute_attenuation(tc, jelppi, LL, ty0, tempi, tempj, tempk, Nx, Nyx, atten);
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
									tempj, jelppi, local_sino, ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
							else {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
									tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  ju, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}

							if (ju > 0) {
								tempijk += Nx;
							}
							else {
								tempijk -= Nx;
							}
							tempj += ju;
							ty0 += tyu;
							xyz = 2u;
						}
						else {
							if (attenuation_correction)
								compute_attenuation(tc, jelppi, LL, tz0, tempi, tempj, tempk, Nx, Nyx, atten);
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
									tempj, jelppi, local_sino, ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);

								if (ku > 0)
									tempijk += Nyx;
								else
									tempijk -= Nyx;
							}
							else if (ii == Np - 1u) {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
									tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  ju, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
							tempk += ku;
							tz0 += tzu;
							xyz = 3u;
						}
						Np_n++;
						if (tempj < 0 || tempk < 0 || tempk >= static_cast<int32_t>(Nz) || tempj >= static_cast<int32_t>(Ny)) {
							if (xyz == 3u && crystal_size_z != 0.f && ii != Np - 1u) {
								tempk -= ku;
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
									tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  ju, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
							break;
						}
					}

					temp = 1. / temp;
					ty0 = ty0_a;
					tz0 = tz0_a;
					tempj = tempj_a;
					tempk = tempk_a;
					if (crystal_size_z == 0.)
						tempijk = Nyx * static_cast<uint32_t>(tempk) + static_cast<uint32_t>(tempj) * Nx;
					else
						tempijk = static_cast<uint32_t>(tempj) * Nx;
					if (attenuation_correction)
						temp *= exp(jelppi);
					if (normalization)
						temp *= norm_coef[lo];


					if (local_sino != 0.) {
						nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
						for (uint32_t ii = 0u; ii < Np_n; ii++) {

							if (ty0 < tz0) {
								if (crystal_size_z == 0.) {
									orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
										tempj, jelppi, local_sino, ax, osem_apu, no_norm, true, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								}
								else {
									orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
										tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  ju, no_norm, true, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								}

								if (ju > 0) {
									tempijk += Nx;
								}
								else {
									tempijk -= Nx;
								}
								tempj += ju;
								ty0 += tyu;
							}
							else {
								if (crystal_size_z == 0.) {
									orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
										tempj, jelppi, local_sino, ax, osem_apu, no_norm, true, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);

									if (ku > 0)
										tempijk += Nyx;
									else
										tempijk -= Nyx;
								}
								else if (ii == Np_n - 1u) {
									orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
										tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  ju, no_norm, true, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								}
								tempk += ku;
								tz0 += tzu;
							}
						}
					}
					else {
						for (uint32_t ii = 0u; ii < Np_n; ii++) {

							if (ty0 < tz0) {
								if (crystal_size_z == 0.) {
									orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
										tempj, jelppi, local_sino, ax, osem_apu, no_norm, false, true, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								}
								else {
									orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
										tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  ju, no_norm, true, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								}

								if (ju > 0) {
									tempijk += Nx;
								}
								else {
									tempijk -= Nx;
								}
								tempj += ju;
								ty0 += tyu;
							}
							else {
								if (crystal_size_z == 0.) {
									orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
										tempj, jelppi, local_sino, ax, osem_apu, no_norm, false, true, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);

									if (ku > 0)
										tempijk += Nyx;
									else
										tempijk -= Nyx;
								}
								else if (ii == Np_n - 1u) {
									orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
										tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  ju, no_norm, true, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								}
								tempk += ku;
								tz0 += tzu;
							}
						}
					}
				}
			}
			else {

				int32_t tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 1;
				double txu = 0., tyu = 0., tzu = 0., tc = 0., tx0 = 0., ty0 = 0., tz0 = 0.;
				const bool skip = siddon_pre_loop_3D(bx, by, bz, x_diff, y_diff, z_diff, maxxx, maxyy, bzb, dx, dy, dz, Nx, Ny, Nz, tempi, tempj, tempk, tyu, txu, tzu,
					Np, TYPE, detectors, tc, iu, ju, ku, tx0, ty0, tz0);

				if (skip)
					continue;

				double temp = 0.;

				if (attenuation_correction)
					LL = sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
				const uint32_t tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
				const double ty0_a = ty0, tz0_a = tz0, tx0_a = tx0;
				uint32_t tempijk;
				if (crystal_size_z == 0.)
					tempijk = Nyx * static_cast<uint32_t>(tempk) + static_cast<uint32_t>(tempj) * Nx;
				else
					tempijk = static_cast<uint32_t>(tempj) * Nx;

				for (uint32_t ii = 0u; ii < Np; ii++) {
					if (tz0 < ty0 && tz0 < tx0) {
						if (attenuation_correction)
							compute_attenuation(tc, jelppi, LL, tz0, tempi, tempj, tempk, Nx, Nyx, atten);
						if (crystal_size_z == 0.) {
							orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
								tempj, jelppi, local_sino, ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							if (ku > 0)
								tempijk += Nyx;
							else
								tempijk -= Nyx;
						}
						else if (ii == Np - 1u) {
							orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
								tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
						}
						tempk += ku;
						tz0 += tzu;
						xyz = 3u;
					}
					else if (ty0 < tx0) {
						if (attenuation_correction)
							compute_attenuation(tc, jelppi, LL, ty0, tempi, tempj, tempk, Nx, Nyx, atten);
						if (crystal_size_z == 0.) {
							orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
								tempj, jelppi, local_sino, ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
						}
						else {
							orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
								tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
						}

						if (ju > 0) {
							tempijk += Nx;
						}
						else {
							tempijk -= Nx;
						}
						tempj += ju;
						ty0 += tyu;
						xyz = 2u;
					}
					else {
						if (attenuation_correction)
							compute_attenuation(tc, jelppi, LL, tx0, tempi, tempj, tempk, Nx, Nyx, atten);
						if (ii == Np - 1u) {
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
									tempj, jelppi, local_sino, ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
							else {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
									tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
						}
						else {
							tempi += iu;
							tx0 += txu;
						}
						xyz = 1u;
					}
					Np_n++;
					if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= static_cast<int32_t>(Nx) || tempj >= static_cast<int32_t>(Ny) || tempk >= static_cast<int32_t>(Nz)) {
						if (ii != Np - 1u && (xyz == 1u || (xyz == 3u && crystal_size_z != 0.))) {
							if (xyz == 1u)
								tempi -= iu;
							else
								tempk -= ku;
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
									tempj, jelppi, local_sino, ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
							else {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
									tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
						}
						break;
					}
				}

				temp = 1. / temp;
				ty0 = ty0_a;
				tx0 = tx0_a;
				tz0 = tz0_a;
				tempi = tempi_a;
				tempj = tempj_a;
				tempk = tempk_a;
				if (crystal_size_z == 0.)
					tempijk = Nyx * static_cast<uint32_t>(tempk) + static_cast<uint32_t>(tempj) * Nx;
				else
					tempijk = static_cast<uint32_t>(tempj) * Nx;
				if (attenuation_correction)
					temp *= exp(jelppi);
				if (normalization)
					temp *= norm_coef[lo];


				if (local_sino != 0.) {
					nominator_mfree(ax, local_sino, epps, temp, randoms_correction, randoms, lo);
					for (uint32_t ii = 0u; ii < Np_n; ii++) {
						if (tz0 < ty0 && tz0 < tx0) {
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
									tempj, jelppi, local_sino, ax, osem_apu, no_norm, true, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								if (ku > 0)
									tempijk += Nyx;
								else
									tempijk -= Nyx;
							}
							else if (ii == Np_n - 1u) {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
									tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, true, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
							tempk += ku;
							tz0 += tzu;
						}
						else if (ty0 < tx0) {
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
									tempj, jelppi, local_sino, ax, osem_apu, no_norm, true, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
							else {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
									tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, true, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}

							if (ju > 0) {
								tempijk += Nx;
							}
							else {
								tempijk -= Nx;
							}
							tempj += ju;
							ty0 += tyu;
						}
						else {
							if (ii == Np_n - 1u) {
								if (crystal_size_z == 0.) {
									orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
										tempj, jelppi, local_sino, ax, osem_apu, no_norm, true, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								}
								else {
									orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
										tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, true, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								}
							}
							else {
								tempi += iu;
								tx0 += txu;
							}
						}
					}
				}
				else {
					for (uint32_t ii = 0u; ii < Np_n; ii++) {
						if (tz0 < ty0 && tz0 < tx0) {
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
									tempj, jelppi, local_sino, ax, osem_apu, no_norm, false, true, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								if (ku > 0)
									tempijk += Nyx;
								else
									tempijk -= Nyx;
							}
							else if (ii == Np_n - 1u) {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
									tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, true, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
							tempk += ku;
							tz0 += tzu;
						}
						else if (ty0 < tx0) {
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
									tempj, jelppi, local_sino, ax, osem_apu, no_norm, false, true, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}
							else {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
									tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, true, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
							}

							if (ju > 0) {
								tempijk += Nx;
							}
							else {
								tempijk -= Nx;
							}
							tempj += ju;
							ty0 += tyu;
						}
						else {
							if (ii == Np_n - 1u) {
								if (crystal_size_z == 0.) {
									orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
										tempj, jelppi, local_sino, ax, osem_apu, no_norm, false, true, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								}
								else {
									orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
										tempj, tempk, jelppi, local_sino, ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, true, OMP, PRECOMPUTE, rhs, Summ, 0, elements, v_indices, idx);
								}
							}
							else {
								tempi += iu;
								tx0 += txu;
							}
						}
					}
				}
			}
		}
	}
}