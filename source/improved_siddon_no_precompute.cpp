/**************************************************************************
* Implements the improved Siddon's algorithm (implementation 1).
* This version also checks whether the LOR/ray intersects the pixel space,
* i.e. it doesn't require any precomputation.
* This version outputs the row and column indices and values that can be
* used to create a sparse matrix.
* This is faster than original Siddon, but slower than precomputed versions.
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

using namespace std;

int improved_siddon_no_precompute(const int64_t loop_var_par, const uint32_t size_x, const double zmax, const uint32_t TotSinos, vector<uint32_t>& indices,
	vector<double>& elements, uint16_t* lor, const double maxyy, const double maxxx, const vector<double>& xx_vec, const double dy,
	const vector<double>& yy_vec, const double* atten, const float* norm_coef, const double* x, const double* y, const double* z_det, const uint32_t NSlices, 
	const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const double dx, const double dz, const double bx, const double by, const double bz, 
	const uint32_t *index, const bool attenuation_correction, const bool normalization, const bool raw, const uint32_t det_per_ring, const uint32_t blocks, 
	const uint32_t block1, const uint16_t *L, const uint32_t *pseudos, const uint32_t pRows, const double global_factor, const bool scatter, const double* scatter_coef, 
	const uint32_t subsets, const double* angles, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t size_y,	const double dPitch, const int64_t nProjections, 
	const uint8_t list_mode) {

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

	const uint32_t Ny_max = Nx * Ny;
	const uint32_t Nz_max = Ny_max * Nz;

	int32_t tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;
	uint32_t Np = 0u;
	double txu = 0., tyu = 0., tzu = 0., tc = 0., tx0 = 0., ty0 = 0., tz0 = 0.;

	for (int64_t lo = 0LL; lo < loop_var_par; lo++) {

#ifndef CT
		if (raw)
			get_detector_coordinates_raw(det_per_ring, x, y, z_det, detectors, L, lo, pseudos, pRows, list_mode);
		else
			get_detector_coordinates_noalloc(x, y, z_det, size_x, detectors, ll, index, lz, TotSinos, lo);
#else
		get_detector_coordinates_CT(x, y, z_det, size_x, detectors, lo, subsets, angles, xy_index, z_index, size_y, dPitch, nProjections, list_mode);
#endif

		// Calculate the x, y and z distances of the detector pair, i.e. the distance between them in the corresponding dimension
		const double y_diff = (detectors.yd - detectors.ys);
		const double x_diff = (detectors.xd - detectors.xs);
		const double z_diff = (detectors.zd - detectors.zs);
		
		if ((y_diff == 0. && x_diff == 0. && z_diff == 0.) || (y_diff == 0. && x_diff == 0.)) {
			continue;
		}
		double local_norm = 0.;
		if (normalization)
			local_norm = static_cast<double>(norm_coef[lo]);

		// If the measurement is on the same ring
		if (fabs(z_diff) < 1e-8 && (fabs(y_diff) < 1e-8 || fabs(x_diff) < 1e-8)) {

			// Z-coordinate
			const uint32_t tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));

			// If the LOR is perpendicular in the y-direction (Siddon cannot be used)
			if (fabs(y_diff) < 1e-8) {

				// Check whether the ray is inside the pixel space
				if (detectors.yd <= maxyy && detectors.yd >= by) {
					uint32_t temp_ijk = 0;

					const double element = perpendicular_elements(Ny, detectors.yd, yy_vec, dx, tempk, Nx, Ny, atten, local_norm, attenuation_correction,
						normalization, temp_ijk, 1u, lo, global_factor, scatter, scatter_coef);

					// Calculate the next index and store it as well as the probability of emission
					for (uint32_t ii = 0u; ii < Nx; ii++) {
						indices.emplace_back((temp_ijk + ii));
						elements.emplace_back(fabs(element));
					}
					lj++;

					lor[lo] = static_cast<uint16_t>(Nx);
					continue;
				}
				continue;
			}
			// Same as for the y-case above
			else if (fabs(x_diff) < 1e-8) {

				if (detectors.xd <= maxxx && detectors.xd >= bx) {
					uint32_t temp_ijk = 0u;

					const double element = perpendicular_elements(1u, detectors.xd, xx_vec, dy, tempk, Ny, Nx, atten, local_norm, attenuation_correction,
						normalization, temp_ijk, Nx, lo, global_factor, scatter, scatter_coef);

					for (uint32_t ii = 0u; ii < Ny_max; ii += Ny) {
						indices.emplace_back((temp_ijk + ii));
						elements.emplace_back(fabs(element));
					}
					lj++;

					lor[lo] = static_cast<uint16_t>(Ny);
					continue;
				}
				continue;
			}
		}
		else {
			int32_t tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;
			double txu = 0., tyu = 0., tzu = 0., tc = 0., tx0 = 1e8, ty0 = 1e8, tz0 = 1e8;
			bool skip = false;

			// Detectors on same ring
			if (std::fabs(z_diff) < 1e-8) {
				tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));
				skip = siddon_pre_loop_2D(bx, by, x_diff, y_diff, maxxx, maxyy, dx, dy, Nx, Ny, tempi, tempj, txu, tyu, Np, TYPE,
					detectors.ys, detectors.xs, detectors.yd, detectors.xd, tc, iu, ju, tx0, ty0);
			}
			// Detectors on different rings (e.g. oblique sinograms)
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

			//if (lo == 651) {
			//	mexPrintf("tx0 = %.10f\n", tx0);
			//	mexPrintf("ty0 = %.10f\n", ty0);
			//	mexPrintf("tz0 = %.10f\n", tz0);
			//	mexPrintf("tzu = %.10f\n", tzu);
			//	mexPrintf("tyu = %.10f\n", tyu);
			//	mexPrintf("txu = %.10f\n", txu);
			//	mexPrintf("tc = %.10f\n", tc);
			//	mexPrintf("tempk = %d\n", tempk);
			//	mexPrintf("tempi = %d\n", tempi);
			//	mexPrintf("tempj = %d\n", tempj);
			//	mexPrintf("detectors.xs = %f\n", detectors.xs);
			//	mexPrintf("detectors.xd = %f\n", detectors.xd);
			//	mexPrintf("detectors.ys = %f\n", detectors.ys);
			//	mexPrintf("detectors.yd = %f\n", detectors.yd);
			//	mexPrintf("detectors.zs = %f\n", detectors.zs);
			//	mexPrintf("detectors.zd = %f\n", detectors.zd);
			//	mexEvalString("pause(.001);");
			//}
			if (skip)
				continue;

			// d_conv
			const double LL = sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);

			// \alpha_c

#ifndef CT
			double temp = 0.;
#endif

			vector<double> templ_ijk;
			vector<uint32_t> temp_koko;
			temp_koko.reserve(Np);
			templ_ijk.reserve(Np);

			uint32_t tempijk = static_cast<uint32_t>(tempk) * Ny_max + static_cast<uint32_t>(tempj) * Nx + static_cast<uint32_t>(tempi);

			// Compute the indices and matrix elements
			for (uint32_t ii = 0u; ii < Np; ii++) {

				temp_koko.emplace_back(tempijk);

				if (tx0 < ty0 && tx0 < tz0) {

					// (30)
					templ_ijk.emplace_back(pixel_value(tx0, tc, LL));

					// (32)
					tempi += iu;
					if (iu > 0)
						tempijk++;
					else
						tempijk--;
					// (33)
					tc = tx0;
					// (34)
					tx0 += txu;

#ifndef CT
					temp += templ_ijk[ii];
#endif

				}
				else if (ty0 < tz0) {

					// (35)
					templ_ijk.emplace_back(pixel_value(ty0, tc, LL));

					// (37)
					tempj += ju;
					if (ju > 0)
						tempijk += Nx;
					else
						tempijk -= Nx;
					// (38)
					tc = ty0;
					// (39)
					ty0 += tyu;

#ifndef CT
					temp += templ_ijk[ii];
#endif
				}
				else {
					templ_ijk.emplace_back(pixel_value(tz0, tc, LL));

					tempk += ku;
					if (ku > 0)
						tempijk += Ny_max;
					else
						tempijk -= Ny_max;
					tc = tz0;
					tz0 += tzu;

#ifndef CT
					temp += templ_ijk[ii];
#endif
				}
				// If the ray/LOR has reached the end of the pixel space
				if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= static_cast<int32_t>(Nx) || tempj >= static_cast<int32_t>(Ny) || tempk >= static_cast<int32_t>(Nz))
					break;
			}

#ifndef CT
			temp = 1. / temp;

			if (attenuation_correction)
				att_corr_vec(templ_ijk, temp_koko, atten, temp, templ_ijk.size());
			if (normalization)
				temp *= local_norm;
			temp *= global_factor;
#endif

			for (uint32_t ii = 0u; ii < templ_ijk.size(); ii++) {
#ifndef CT
				elements.emplace_back((fabs(templ_ijk[ii]) * temp));
#else
				elements.emplace_back(fabs(templ_ijk[ii]));
#endif
				indices.emplace_back(temp_koko[ii]);
			}

			lor[lo] = static_cast<uint16_t>(templ_ijk.size());
			lj++;
			continue;
		}
	}
	return lj;
}