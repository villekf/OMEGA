/**************************************************************************
* Implements the original Siddon's algorithm for OMEGA.
* This version also checks whether the LOR/ray intersects the pixel space,
* i.e. it doesn't require any precomputation.
* This version outputs the row and column indices and values that can be
* used to create a sparse matrix.
* This is the slowest Siddon code, but faster projector than orthogonal.
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

using namespace std;

int original_siddon_no_precompute(const int64_t loop_var_par, const uint32_t size_x, const double zmax, const uint32_t TotSinos, vector<uint32_t>& indices,
	vector<double>& elements, uint16_t* lor, const double maxyy, const double maxxx, const vector<double>& xx_vec, const double dy, 
	const vector<double>& yy_vec, const double* atten, const double* norm_coef, const double* x, const double* y, const double* z_det, const uint32_t NSlices, 
	const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const double dx, const double dz, const double bx, const double by, const double bz, 
	const uint32_t *index, const bool attenuation_correction, const bool normalization, const bool raw, const uint32_t det_per_ring, const uint32_t blocks, 
	const uint32_t block1, const uint16_t *L, const uint32_t *pseudos, const uint32_t pRows, const vector<double>& iij_vec, const vector<double>& jjk_vec, 
	const vector<double>& kkj_vec, const double global_factor, const bool scatter, const double* scatter_coef) {

	int ll;
	if (raw)
		ll = 0;
	else
		ll = -1;
	int lz = -1;
	int lj = 0;

	// Precompute variables
	const double bzb = bz + static_cast<double>(Nz) * dz;

	Det detectors;


	for (int64_t lo = 0LL; lo < loop_var_par; lo++) {

		if (raw)
			get_detector_coordinates_raw(det_per_ring, x, y, z_det, detectors, L, ll, pseudos, pRows);
		else
			get_detector_coordinates_noalloc(x, y, z_det, size_x, detectors, ll, index, lz, TotSinos, lo);

		// Calculate the x, y and z distances of the detector pair
		const double y_diff = (detectors.yd - detectors.ys);
		const double x_diff = (detectors.xd - detectors.xs);
		const double z_diff = (detectors.zd - detectors.zs);

		if (abs(z_diff) < 1e-8f) {

			uint32_t tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));

			if (abs(y_diff) < 1e-8) {


				if (detectors.yd <= maxyy && detectors.yd >= by) {
					uint32_t temp_ijk = 0u;

					const double element = perpendicular_elements(Ny, detectors.yd, yy_vec, dx, tempk, Nx, Ny, atten, norm_coef, attenuation_correction, 
						normalization, temp_ijk, 1, lo, global_factor, scatter, scatter_coef);

					// Calculate the next index and store it as well as the probability of emission
					for (uint32_t ii = 0u; ii < Nx; ii++) {
						indices.emplace_back((tempk + ii));
						elements.emplace_back(fabs(element));
					}
					lj++;

					lor[lo] = static_cast<uint16_t>(Nx);
					continue;
				}
			}
			else if (abs(x_diff) < 1e-8f) {

				if (detectors.xd <= maxxx && detectors.xd >= bx) {
					uint32_t temp_ijk = 0u;

					const double element = perpendicular_elements(1, detectors.xd, xx_vec, dy, tempk, Ny, Nx, atten, norm_coef, attenuation_correction, 
						normalization, temp_ijk, Nx, lo, global_factor, scatter, scatter_coef);

					for (uint32_t ii = 0u; ii < Ny; ii++) {
						indices.emplace_back((tempk + ii * Ny));
						elements.emplace_back(fabs(element));
					}
					lj++;

					lor[lo] = static_cast<uint16_t>(Ny);
					continue;
				}
				continue;
			}

			double tx0 = (bx - detectors.xs) / (x_diff);
			double ty0 = (by - detectors.ys) / (y_diff);
			const double txback = (maxxx - detectors.xs) / (x_diff);
			const double tyback = (maxyy - detectors.ys) / (y_diff);

			const double txmin = min(tx0, txback);
			const double txmax = max(tx0, txback);
			const double tymin = min(ty0, tyback);
			const double tymax = max(ty0, tyback);

			const double tmin = max(txmin, tymin);
			const double tmax = min(txmax, tymax);

			if (tmin >= tmax)
				continue;

			uint32_t imin, imax, jmin, jmax;
			int iu, ju;

			if (detectors.xs < detectors.xd) {
				d_g_s_precomp(tmin, txmin, tmax, txmax, imin, imax, tx0, iu, x_diff, bx, dx, detectors.xs, Nx);
			}
			else if (detectors.xs > detectors.xd) {
				s_g_d_precomp(tmin, txmin, tmax, txmax, imin, imax, tx0, iu, x_diff, bx, dx, detectors.xs, Nx);
			}
			if (detectors.ys < detectors.yd) {
				d_g_s_precomp(tmin, tymin, tmax, tymax, jmin, jmax, ty0, ju, y_diff, by, dy, detectors.ys, Ny);
			}
			else if (detectors.ys > detectors.yd) {
				s_g_d_precomp(tmin, tymin, tmax, tymax, jmin, jmax, ty0, ju, y_diff, by, dy, detectors.ys, Ny);
			}


			vector<double> tx_n(static_cast<int64_t>(imax) - static_cast<int64_t>(imin) + 1LL, 0.);
			vector<double> ty_n(static_cast<int64_t>(jmax) - static_cast<int64_t>(jmin) + 1LL, 0.);

			if (imax > imin || imin != 0u) {
				uint32_t tl = 0u;
				for (uint32_t ii = imin; ii <= imax; ii++) {
					tx_n[tl] = (bx + iij_vec[ii] * dx - detectors.xs) / (x_diff);
					tl++;
				}
			}
			else
				tx_n[0] = tx0;

			if (jmax > jmin || jmin != 0u) {
				uint32_t tl = 0u;
				for (uint32_t ii = jmin; ii <= jmax; ii++) {
					ty_n[tl] = (by + jjk_vec[ii] * dy - detectors.ys) / (y_diff);

					tl++;
				}
			}
			else
				ty_n[0] = ty0;

			const double LL = sqrt(x_diff * x_diff + y_diff * y_diff);

			vector<double> tt;
			tt.reserve(tx_n.size() + ty_n.size() + 1);
			tt.emplace_back(tmin);
			tt.insert(tt.end(), tx_n.begin(), tx_n.end());
			tt.insert(tt.end(), ty_n.begin(), ty_n.end());

			sort(tt.begin(), tt.end());
			tt.erase(unique(tt.begin(), tt.end()), tt.end());

			uint32_t koko = static_cast<uint32_t>(tt.size() - 2ULL);

			uint32_t tempi, tempj;
			vector<double> templ_ijk(tt.size() - 1ULL, 0.);
			vector<uint32_t> temp_koko(tt.size() - 1ULL, 0);

			tempk *= (Ny * Nx);

			double jelppi;

			double temp = 0.;

			for (size_t ii = 0ULL; ii <= static_cast<size_t>(koko); ii++) {
				jelppi = tt[ii + 1ULL] + tt[ii];
				double pt = detectors.xs + (jelppi / 2.) * x_diff;
				tempi = static_cast<uint32_t>(floor((pt - bx) / dx));

				pt = detectors.ys + (jelppi / 2.) * y_diff;
				tempj = static_cast<uint32_t>(floor((pt - by) / dy));

				temp_koko[ii] = tempj * Ny + tempi + tempk;

				jelppi = tt[ii + 1ULL] - tt[ii];

				templ_ijk[ii] = jelppi * LL;

				temp += templ_ijk[ii];
			}


			temp = 1. / temp;

			//auto i = sort_indexes(temp_koko);

			if (attenuation_correction == true) {

				jelppi = 0.;

				for (uint32_t iii = 0u; iii <= koko; iii++) {
					jelppi += templ_ijk[iii] * -atten[temp_koko[iii]];
				}
				temp *= exp(jelppi);
			}
			if (normalization)
				temp *= norm_coef[lo];
			temp *= global_factor;

			for (uint32_t ii = 0u; ii <= koko; ii++) {
				//indices.emplace_back(temp_koko[i[ii]]);
				//elements.emplace_back(static_cast<float>(templ_ijk[i[ii]] * temp));
				indices.emplace_back(temp_koko[ii]);
				elements.emplace_back(static_cast<float>(templ_ijk[ii] * temp));

			}

			lor[lo] = static_cast<uint16_t>(templ_ijk.size());
			lj++;
			continue;
		}

		else {
			if (abs(y_diff) < 1e-8f) {
				if (detectors.yd <= maxyy && detectors.yd >= by) {

					double tx0 = (bx - detectors.xs) / (x_diff);
					double tz0 = (bz - detectors.zs) / (z_diff);
					const double txback = (maxxx - detectors.xs) / (x_diff);
					const double tzback = (bzb - detectors.zs) / (z_diff);

					double txmin = min(tx0, txback);
					double txmax = max(tx0, txback);
					double tzmin = min(tz0, tzback);
					double tzmax = max(tz0, tzback);

					double tmin = max(txmin, tzmin);
					double tmax = min(txmax, tzmax);

					if (tmin >= tmax)
						continue;

					uint32_t imin, imax, kmin, kmax;
					int iu, ku;


					if (detectors.xs < detectors.xd) {
						d_g_s_precomp(tmin, txmin, tmax, txmax, imin, imax, tx0, iu, x_diff, bx, dx, detectors.xs, Nx);
					}
					else if (detectors.xs > detectors.xd) {
						s_g_d_precomp(tmin, txmin, tmax, txmax, imin, imax, tx0, iu, x_diff, bx, dx, detectors.xs, Nx);
					}

					if (detectors.zs < detectors.zd) {
						d_g_s_precomp(tmin, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);

					}
					else if (detectors.zs > detectors.zd) {
						s_g_d_precomp(tmin, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);
					}

					vector<double> tx_n(static_cast<int64_t>(imax) - static_cast<int64_t>(imin) + 1LL, 0.);
					vector<double> tz_n(static_cast<int64_t>(kmax) - static_cast<int64_t>(kmin) + 1LL, 0.);

					if (imax > imin || imin != 0u) {
						uint32_t tl = 0u;
						for (uint32_t ii = imin; ii <= imax; ii++) {
							tx_n[tl] = (bx + iij_vec[ii] * dx - detectors.xs) / (x_diff);
							tl++;
						}
					}
					else
						tx_n[0] = tx0;

					if (kmax > kmin || kmin != 0u) {
						uint32_t tl = 0u;
						for (uint32_t ii = kmin; ii <= kmax; ii++) {
							tz_n[tl] = (bz + kkj_vec[ii] * dz - detectors.zs) / (z_diff);
							tl++;
						}
					}
					else
						tz_n[0] = tz0;


					const double LL = sqrt(x_diff * x_diff + z_diff * z_diff);

					vector<double> tt;
					tt.reserve(tx_n.size() + tz_n.size() + 1);
					tt.emplace_back(tmin);
					tt.insert(tt.end(), tx_n.begin(), tx_n.end());
					tt.insert(tt.end(), tz_n.begin(), tz_n.end());

					sort(tt.begin(), tt.end());
					tt.erase(unique(tt.begin(), tt.end()), tt.end());


					uint32_t tempi, tempj, tempk;
					double apu2, apu1;

					for (uint32_t ii = 0u; ii < Ny; ii++) {
						apu1 = abs(yy_vec[ii] - detectors.yd);
						if (ii > 0 && apu1 < apu2 || ii == 0u) {
							tempj = ii;
							apu2 = apu1;
						}

					}

					tempj = tempj * Ny;

					uint32_t koko = static_cast<uint32_t>(tt.size() - 2ULL);

					vector<double> templ_ijk(tt.size() - 1ULL, 0.);
					vector<uint32_t> temp_koko(tt.size() - 1ULL, 0u);

					double jelppi;

					double temp = 0.;

					for (size_t ii = 0ULL; ii < tt.size() - 1u; ii++) {
						jelppi = (tt[ii + 1ULL] + tt[ii]);

						double pt = detectors.xs + (jelppi / 2.) * x_diff;
						tempi = static_cast<uint32_t>(floor((pt - bx) / dx));

						pt = detectors.zs + (jelppi / 2.) * z_diff;
						tempk = static_cast<uint32_t>(floor((pt - bz) / dz));

						temp_koko[ii] = (tempj + tempi + tempk * Ny * Nx);

						jelppi = tt[ii + 1ULL] - tt[ii];

						templ_ijk[ii] = jelppi * LL;

						temp += templ_ijk[ii];
					}

					temp = 1 / temp;

					//auto i = sort_indexes(temp_koko);

					if (attenuation_correction == true) {

						jelppi = 0.;

						for (uint32_t iii = 0u; iii <= koko; iii++) {
							jelppi += templ_ijk[iii] * -atten[temp_koko[iii]];
						}
						temp *= exp(jelppi);
					}
					if (normalization)
						temp *= norm_coef[lo];
					if (scatter)
						temp *= scatter_coef[lo];
					temp *= global_factor;

					for (uint32_t ii = 0u; ii <= koko; ii++) {
						//indices.emplace_back(temp_koko[i[ii]]);
						//elements.emplace_back(static_cast<float>(templ_ijk[i[ii]] * temp));
						indices.emplace_back(temp_koko[ii]);
						elements.emplace_back(static_cast<float>(templ_ijk[ii] * temp));
					}


					lor[lo] = static_cast<uint16_t>(templ_ijk.size());
					lj++;
					continue;

				}
				continue;
			}
			else if (abs(x_diff) < 1e-8f) {
				if (detectors.xd <= maxxx && detectors.xd >= bx) {

					double ty0 = (by - detectors.ys) / (y_diff);
					const double tyback = (maxyy - detectors.ys) / (y_diff);
					double tz0 = (bz - detectors.zs) / (z_diff);
					const double tzback = (bzb - detectors.zs) / (z_diff);

					double tymin = min(ty0, tyback);
					double tymax = max(ty0, tyback);
					double tzmin = min(tz0, tzback);
					double tzmax = max(tz0, tzback);

					double tmin = max(tymin, tzmin);
					double tmax = min(tymax, tzmax);

					if (tmin >= tmax)
						continue;


					uint32_t jmin, jmax, kmin, kmax;

					int ku, ju;

					if (detectors.ys < detectors.yd) {
						d_g_s_precomp(tmin, tymin, tmax, tymax, jmin, jmax, ty0, ju, y_diff, by, dy, detectors.ys, Ny);

					}
					else if (detectors.ys > detectors.yd) {
						s_g_d_precomp(tmin, tymin, tmax, tymax, jmin, jmax, ty0, ju, y_diff, by, dy, detectors.ys, Ny);
					}

					if (detectors.zs < detectors.zd) {
						d_g_s_precomp(tmin, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);

					}
					else if (detectors.zs > detectors.zd) {
						s_g_d_precomp(tmin, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);
					}

					vector<double> ty_n(static_cast<int64_t>(jmax) - static_cast<int64_t>(jmin) + 1LL, 0.);
					vector<double> tz_n(static_cast<int64_t>(kmax) - static_cast<int64_t>(kmin) + 1LL, 0.);

					if (jmax > jmin || jmin != 0u) {
						uint32_t tl = 0u;
						for (uint32_t ii = jmin; ii <= jmax; ii++) {
							ty_n[tl] = (by + jjk_vec[ii] * dy - detectors.ys) / (y_diff);
							tl++;
						}
					}
					else
						ty_n[0] = ty0;

					if (kmax > kmin || kmin != 0u) {
						uint32_t tl = 0u;
						for (uint32_t ii = kmin; ii <= kmax; ii++) {
							tz_n[tl] = (bz + kkj_vec[ii] * dz - detectors.zs) / (z_diff);
							tl++;
						}
					}
					else
						tz_n[0] = tz0;


					const double LL = sqrt(y_diff * y_diff + z_diff * z_diff);

					vector<double> tt;
					tt.reserve(ty_n.size() + tz_n.size() + 1u);
					tt.emplace_back(tmin);
					tt.insert(tt.end(), ty_n.begin(), ty_n.end());
					tt.insert(tt.end(), tz_n.begin(), tz_n.end());

					sort(tt.begin(), tt.end());
					tt.erase(unique(tt.begin(), tt.end()), tt.end());

					uint32_t tempi, tempj, tempk;
					double apu2, apu1;

					for (uint32_t ii = 0u; ii < Nx; ii++) {
						apu1 = abs(xx_vec[ii] - detectors.xd);
						if (ii > 0u && apu1 < apu2 || ii == 0u) {
							tempi = ii;
							apu2 = apu1;
						}

					}

					//tempi = tempi * Ny;

					uint32_t koko = static_cast<uint32_t>(tt.size() - 2ULL);

					vector<double> templ_ijk(tt.size() - 1ULL, 0.);
					vector<uint32_t> temp_koko(tt.size() - 1ULL, 0u);

					double jelppi;

					double temp = 0.;

					for (size_t ii = 0ULL; ii < tt.size() - 1u; ii++) {
						jelppi = (tt[ii + 1ULL] + tt[ii]);

						double pt = detectors.ys + (jelppi / 2.) * y_diff;
						tempj = static_cast<uint32_t>(floor((pt - by) / dy));

						pt = detectors.zs + (jelppi / 2.) * z_diff;
						tempk = static_cast<uint32_t>(floor((pt - bz) / dz));

						temp_koko[ii] = (tempj * Ny + tempi + tempk * Ny * Nx);

						jelppi = tt[ii + 1ULL] - tt[ii];

						templ_ijk[ii] = jelppi * LL;

						temp += templ_ijk[ii];

					}

					temp = 1 / temp;

					//auto i = sort_indexes(temp_koko);

					if (attenuation_correction == true) {

						jelppi = 0.;

						for (uint32_t iii = 0u; iii <= koko; iii++) {
							jelppi += templ_ijk[iii] * -atten[temp_koko[iii]];
						}
						temp *= exp(jelppi);

					}
					if (normalization)
						temp *= norm_coef[lo];
					if (scatter)
						temp *= scatter_coef[lo];
					temp *= global_factor;

					for (uint32_t ii = 0u; ii <= koko; ii++) {
						//indices.emplace_back(temp_koko[i[ii]]);
						//elements.emplace_back(static_cast<float>(templ_ijk[i[ii]] * temp));
						indices.emplace_back(temp_koko[ii]);
						elements.emplace_back(static_cast<float>(templ_ijk[ii] * temp));
					}


					lor[lo] = static_cast<uint16_t>(templ_ijk.size());
					lj++;
					continue;

				}
				continue;
			}



			double tx0 = (bx - detectors.xs) / (x_diff);
			double ty0 = (by - detectors.ys) / (y_diff);
			double tz0 = (bz - detectors.zs) / (z_diff);
			const double txback = (maxxx - detectors.xs) / (x_diff);
			const double tyback = (maxyy - detectors.ys) / (y_diff);
			const double tzback = (bzb - detectors.zs) / (z_diff);

			double txmin = min(tx0, txback);
			double txmax = max(tx0, txback);
			double tymin = min(ty0, tyback);
			double tymax = max(ty0, tyback);
			double tzmin = min(tz0, tzback);
			double tzmax = max(tz0, tzback);

			double tmin = max(max(txmin, tzmin), tymin);
			double tmax = min(min(txmax, tzmax), tymax);

			if (tmin >= tmax)
				continue;

			uint32_t imin, imax, jmin, jmax, kmin, kmax;
			int iu, ju, ku;

			if (detectors.xs < detectors.xd) {
				d_g_s_precomp(tmin, txmin, tmax, txmax, imin, imax, tx0, iu, x_diff, bx, dx, detectors.xs, Nx);
			}
			else {
				s_g_d_precomp(tmin, txmin, tmax, txmax, imin, imax, tx0, iu, x_diff, bx, dx, detectors.xs, Nx);
			}
			if (detectors.ys < detectors.yd) {
				d_g_s_precomp(tmin, tymin, tmax, tymax, jmin, jmax, ty0, ju, y_diff, by, dy, detectors.ys, Ny);
			}
			else {
				s_g_d_precomp(tmin, tymin, tmax, tymax, jmin, jmax, ty0, ju, y_diff, by, dy, detectors.ys, Ny);
			}
			if (detectors.zs < detectors.zd) {
				d_g_s_precomp(tmin, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);
			}
			else {
				s_g_d_precomp(tmin, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);
			}

			vector<double> tx_n(static_cast<int64_t>(imax) - static_cast<int64_t>(imin) + 1LL, 0.);
			vector<double> tz_n(static_cast<int64_t>(kmax) - static_cast<int64_t>(kmin) + 1LL, 0.);
			vector<double> ty_n(static_cast<int64_t>(jmax) - static_cast<int64_t>(jmin) + 1LL, 0.);

			if (jmax > jmin || jmin != 0u) {
				uint32_t tl = 0u;
				for (uint32_t ii = jmin; ii <= jmax; ii++) {
					ty_n[tl] = (by + jjk_vec[ii] * dy - detectors.ys) / (y_diff);
					tl++;
				}
			}
			else
				ty_n[0] = ty0;

			if (imax > imin || imin != 0u) {
				uint32_t tl = 0u;
				for (uint32_t ii = imin; ii <= imax; ii++) {
					tx_n[tl] = (bx + iij_vec[ii] * dx - detectors.xs) / (x_diff);
					tl++;
				}
			}
			else
				tx_n[0] = tx0;

			if (kmax > kmin || kmin != 0u) {
				uint32_t tl = 0u;
				for (uint32_t ii = kmin; ii <= kmax; ii++) {
					tz_n[tl] = (bz + kkj_vec[ii] * dz - detectors.zs) / (z_diff);
					tl++;
				}
			}
			else
				tz_n[0] = tz0;

			const double LL = sqrt(x_diff * x_diff + z_diff * z_diff + y_diff * y_diff);

			vector<double> tt;
			tt.reserve(tx_n.size() + tz_n.size() + ty_n.size() + 1);
			tt.emplace_back(tmin);
			tt.insert(tt.end(), tx_n.begin(), tx_n.end());
			tt.insert(tt.end(), ty_n.begin(), ty_n.end());
			tt.insert(tt.end(), tz_n.begin(), tz_n.end());

			sort(tt.begin(), tt.end());
			tt.erase(unique(tt.begin(), tt.end()), tt.end());

			uint32_t tempi, tempj, tempk;

			uint32_t koko = static_cast<uint32_t>(tt.size() - 2ULL);

			vector<double> templ_ijk(tt.size() - 1ULL, 0.);
			vector<uint32_t> temp_koko(tt.size() - 1ULL, 0u);

			double jelppi;

			double temp = 0.;

			for (size_t ii = 0ULL; ii < tt.size() - 1; ii++) {
				jelppi = (tt[ii + 1ULL] + tt[ii]);

				double pt = detectors.ys + (jelppi / 2.) * y_diff;
				tempj = static_cast<uint32_t>(floor((pt - by) / dy));

				pt = detectors.zs + (jelppi / 2.) * z_diff;
				tempk = static_cast<uint32_t>(floor((pt - bz) / dz));

				pt = detectors.xs + (jelppi / 2.) * x_diff;
				tempi = static_cast<uint32_t>(floor((pt - bx) / dx));

				temp_koko[ii] = (tempj * Ny + tempi + tempk * Ny * Nx);

				jelppi = tt[ii + 1ULL] - tt[ii];

				templ_ijk[ii] = jelppi * LL;

				temp += templ_ijk[ii];
			}

			temp = 1 / temp;

			if (attenuation_correction == true) {

				jelppi = 0.;

				for (uint32_t iii = 0u; iii <= koko; iii++) {
					jelppi = jelppi + templ_ijk[iii] * -atten[temp_koko[iii]];
				}
				temp *= exp(jelppi);

			}
			if (normalization)
				temp *= norm_coef[lo];
			temp *= global_factor;

			//auto i = sort_indexes(temp_koko);

			for (uint32_t ii = 0u; ii <= koko; ii++) {
				//indices.emplace_back(temp_koko[i[ii]]);
				//elements.emplace_back(static_cast<float>(templ_ijk[i[ii]] * temp));
				indices.emplace_back(temp_koko[ii]);
				elements.emplace_back(static_cast<float>(templ_ijk[ii] * temp));
			}

			lor[lo] = static_cast<uint16_t>(templ_ijk.size());
			lj++;
			continue;

		}

	}
	return lj;
}