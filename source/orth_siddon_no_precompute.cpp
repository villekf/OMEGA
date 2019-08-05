/**************************************************************************
* Implements the Orthogonal Siddon's algorithm for OMEGA.
* This version also checks whether the LOR/ray intersects the pixel space,
* i.e. it doesn't require any precomputation.
* This version outputs the row and column indices and values that can be
* used to create a sparse matrix.
* This is even slower than the original Siddon code.
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

constexpr int TYPE = 0;

// Whether to use the OpenMP code or not
constexpr bool OMP = false;

// Using non-OpenMP with either precomputation or without
constexpr bool PRECOMPUTE = false;

using namespace std;

int orth_siddon_no_precompute(const uint32_t loop_var_par, const uint32_t size_x, const double zmax, const uint32_t TotSinos, vector<uint32_t>& indices, 
	vector<double>& elements, uint16_t* lor, const double maxyy, const double maxxx, const vector<double>& xx_vec, const double dy, 
	const vector<double>& yy_vec, const double* atten, const double* norm_coef, const double* x, const double* y, const double* z_det, const uint32_t NSlices, 
	const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const double dx, const double dz, const double bx, const double by, const double bz, 
	const uint32_t *index, const bool attenuation_correction, const bool normalization, const bool raw, const uint32_t det_per_ring, const uint32_t blocks, 
	const uint32_t block1, const uint16_t *L, const uint32_t *pseudos, const uint32_t pRows, const double crystal_size, const double crystal_size_z, 
	const double* y_center, const double* x_center, const double* z_center, const int32_t dec_v) {

	int ll;
	if (raw)
		ll = 0;
	else
		ll = -1;
	int lz = -1;
	int lj = 0;
	int lo = -1;

	// Precompute variables
	const double bzb = bz + static_cast<double>(Nz) * dz;
	const uint32_t Nyx = Ny * Nx;

	Det detectors;

	const uint32_t Ny_max = Nx * Ny;
	const uint32_t Nz_max = Ny_max * Nz;


	double ax = 0.;
	double* osem_apu = nullptr;
	double* Summ = nullptr;
	double* rhs = nullptr;
	const bool no_norm = true;

	const int32_T dec = static_cast<int32_T>(ceil(crystal_size_z / sqrt(dz * dz * 2.))) * dec_v;


	for (uint32_t lo = 0; lo < loop_var_par; lo++) {

		if (raw)
			get_detector_coordinates_raw(det_per_ring, x, y, z_det, detectors, L, ll, pseudos, pRows);
		else
			get_detector_coordinates_noalloc(x, y, z_det, size_x, detectors, ll, index, lz, TotSinos, lo);

		if (detectors.zd < detectors.zs) {
			double tempa = detectors.zd;
			detectors.zd = detectors.zs;
			detectors.zs = tempa;
			tempa = detectors.xd;
			detectors.xd = detectors.xs;
			detectors.xs = tempa;
			tempa = detectors.yd;
			detectors.yd = detectors.ys;
			detectors.ys = tempa;
		}

		// Calculate the x, y and z distances of the detector pair
		const double y_diff = (detectors.yd - detectors.ys);
		const double x_diff = (detectors.xd - detectors.xs);
		const double z_diff = (detectors.zd - detectors.zs);

		double kerroin, length_, jelppi = 0., LL;

		if (crystal_size_z == 0.) {
			kerroin = detectors.xd * detectors.ys - detectors.yd * detectors.xs;
			length_ = sqrt(y_diff * y_diff + x_diff * x_diff) * crystal_size;
		}
		else {
			kerroin = norm(x_diff, y_diff, z_diff) * crystal_size_z;
		}
		uint32_t Np = 0u;
		uint32_t xyz = 0u;
		size_t idx = 0ULL;

		// If the measurement is on a same ring
		if (fabs(z_diff) < 1e-8) {

			// Z-coordinate
			const uint32_t tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));

			// If the LOR is perpendicular in the y-direction (Siddon cannot be used)
			if (fabs(y_diff) < 1e-8) {

				// Check whether the ray is inside the pixel space
				if (detectors.yd <= maxyy && detectors.yd >= by) {

					const double length = std::sqrt(x_diff * x_diff) * crystal_size;
					uint32_t temp_ijk = 0;
					double temp = 0.;
					int hpk = 0;
					if (crystal_size_z == 0.) {
						std::vector<double> vec1;
						std::vector<double> vec2;
						
						orth_perpendicular(Ny, detectors.yd, yy_vec, dx, tempk, Nx, Ny, -x_diff, y_center, kerroin, length_, vec1, vec2, temp, temp_ijk);

						if (attenuation_correction)
							att_corr_scalar_orth(temp_ijk, atten, temp, 1, Nx, dx);

						// Calculate the next index and store it as well as the probability of emission
						for (int32_t jj = -static_cast<int32_t>(vec1.size()) + 1; jj <= static_cast<int32_t>(vec2.size()); jj++) {
							double element;
							if (jj <= 0)
								element = vec1[-jj];
							else
								element = vec2[static_cast<size_t>(jj) - 1ULL];
							element *= temp;
							for (uint32_t ii = 0; ii < Ny; ii++) {
								indices.emplace_back((temp_ijk + ii + jj * Ny));
								elements.emplace_back(element);
								hpk++;
							}
						}
						lj++;
						lor[lo] = hpk;
						continue;
					}
					else {
						double temppi = detectors.xs;
						detectors.xs = detectors.ys;
						detectors.ys = temppi;
						orth_perpendicular_np_3D(detectors.yd, yy_vec, tempk, Nx, Ny, Nz, Nyx, Ny, 1u, detectors, y_diff, x_diff, z_diff, y_center, x_center[0],
							z_center, kerroin, hpk, temp, temp_ijk, indices, elements);

						if (attenuation_correction)
							att_corr_scalar_orth(temp_ijk, atten, temp, 1u, Nx, dx);
						if (normalization)
							temp *= norm_coef[lo];

						for (size_t jj = elements.size() - hpk; jj < elements.size(); jj++) {
							elements[jj] *= temp;
						}
						lj++;
						lor[lo] = hpk;
						continue;
					}
				}
				continue;
			}
			// Same as for the y-case above
			else if (fabs(x_diff) < 1e-8) {

				if (detectors.xd <= maxxx && detectors.xd >= bx) {
					const double length = std::sqrt(y_diff * y_diff) * crystal_size;
					uint32_t temp_ijk = 0;
					double temp = 0.;
					int hpk = 0;
					if (crystal_size_z == 0.) {
						std::vector<double> vec1;
						std::vector<double> vec2;

						orth_perpendicular(1u, detectors.xd, xx_vec, dy, tempk, Ny, Nx, y_diff, x_center, kerroin, length_, vec1, vec2, temp, temp_ijk);


						if (attenuation_correction)
							att_corr_scalar_orth(temp_ijk, atten, temp, Ny, 1, dy);
						if (normalization)
							temp *= norm_coef[lo];

						// Calculate the next index and store it as well as the probability of emission
						for (uint32_t ii = 0; ii < Nx; ii++) {
							for (int32_t jj = -static_cast<int32_t>(vec1.size()) + 1; jj <= static_cast<int32_t>(vec2.size()); jj++) {
								double element;
								if (jj <= 0)
									element = vec1[-jj];
								else
									element = vec2[static_cast<size_t>(jj) - 1ULL];
								element *= temp;
								indices.emplace_back((temp_ijk + ii * Ny + jj));
								elements.emplace_back(element);
								hpk++;
							}
						}

						lj++;
						lor[lo] = hpk;
						continue;
					}
					else {
						orth_perpendicular_np_3D(detectors.xd, xx_vec, tempk, Ny, Nx, Nz, Nyx, 1u, Nx, detectors, x_diff, y_diff, z_diff, x_center, y_center[0],
							z_center, kerroin, hpk, temp, temp_ijk, indices, elements);

						if (attenuation_correction)
							att_corr_scalar_orth(temp_ijk, atten, temp, Ny, 1u, dy);
						if (normalization)
							temp *= norm_coef[lo];

						for (size_t jj = elements.size() - hpk; jj < elements.size(); jj++) {
							elements[jj] *= temp;
						}
						lj++;
						lor[lo] = hpk;
						continue;
					}
				}
				continue;
			}
			else {
				// If neither x- nor y-directions are perpendicular
				// Correspond to the equations (9) and (10) from reference [2]
				int32_t tempi = 0, tempj = 0, iu = 0, ju = 0;
				double txu = 0., tyu = 0., tc = 0., tx0 = 0., ty0 = 0.;
				const bool skip = siddon_pre_loop_2D(bx, by, x_diff, y_diff, maxxx, maxyy, dx, dy, Nx, Ny, tempi, tempj, txu, tyu, Np, TYPE,
					detectors.ys, detectors.xs, detectors.yd, detectors.xd, tc, iu, ju, tx0, ty0);

				if (skip)
					continue;

				//vector<double> templ_ijk;
				//vector<uint32_t> indices;
				//indices.reserve(Np * 2ULL);
				//templ_ijk.reserve(Np * 2ULL);
				double d_ort = 0.;
				if (attenuation_correction)
					LL = sqrt(x_diff * x_diff + y_diff * y_diff);
				double temp = 0.;
				uint32_t tempijk;
				if (crystal_size_z == 0.)
					tempijk = Nyx * tempk + static_cast<uint32_t>(tempj) * Nx;
				else
					tempijk = static_cast<uint32_t>(tempj) * Nx;

				size_t koko = elements.size();

				// Compute the indices and matrix elements
				for (uint32_t ii = 0; ii < Np; ii++) {

					if (tx0 < ty0) {
						if (attenuation_correction)
							compute_attenuation(tc, jelppi, LL, tx0, tempi, tempj, tempk, Nx, Nyx, atten);
						if (ii == Np - 1u) {
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u, tempj, jelppi, 
									0., ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, indices, idx);
							}
							else {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u, tempj, 
									tempk, jelppi, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 
									0, elements, indices, idx);
							}
						}
						else {
							// (32)
							tempi += iu;
							// (34)
							tx0 += txu;
						}
						xyz = 1u;

					}
					else {
						if (attenuation_correction)
							compute_attenuation(tc, jelppi, LL, ty0, tempi, tempj, tempk, Nx, Nyx, atten);
						if (crystal_size_z == 0.) {
							orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u, tempj, jelppi, 
								0., ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, indices, idx);
						}
						else {
							orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u, tempj, 
								tempk, jelppi, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 
								0, elements, indices, idx);
						}

						// (37)
						tempj += (ju);
						if (ju > 0)
							tempijk += Nx;
						else
							tempijk -= Nx;
						// (39)
						ty0 += tyu;

						xyz = 2u;
					}
					// If the ray/LOR has reached the end of the pixel space
					if (tempj < 0 || tempi < 0 || tempi >= static_cast<int32_t>(Nx) || tempj >= static_cast<int32_t>(Ny)) {
						if (ii != Np - 1u && xyz == 1) {
							tempi -= iu;
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u, tempj, jelppi, 
									0., ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, indices, idx);
							}
							else {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u, tempj, 
									tempk, jelppi, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 
									0, elements, indices, idx);
							}
						}
						break;
					}
				}

				temp = 1. / temp;

				if (attenuation_correction)
					temp *= jelppi;
				if (normalization)
					temp *= norm_coef[lo];

				for (size_t ii = koko; ii < koko + idx; ii++) {
					elements[ii] *= temp;
				}

				//for (uint32_t ii = 0u; ii < templ_ijk.size(); ii++) {
				//	elements.emplace_back(templ_ijk[ii] * temp);
				//	//indices.emplace_back(indices[ii]);
				//}


				lor[lo] = static_cast<uint16_t>(idx);
				lj++;
				continue;
			}
		}
		// If the z-detector coordinates are not on the same ring
		// All computations follow the same logic as above
		else {
			if (fabs(y_diff) < 1e-8) {
				if (detectors.yd <= maxyy && detectors.yd >= by) {

					int32_t tempi = 0, tempk = 0, tempj = 0, iu = 0, ku = 0;
					double txu = 0., tzu = 0., tc = 0., tx0 = 0., tz0 = 0.;

					const bool skip = siddon_pre_loop_2D(bx, bz, x_diff, z_diff, maxxx, bzb, dx, dz, Nx, Nz, tempi, tempk, txu, tzu, Np, TYPE, detectors.zs, 
						detectors.xs, detectors.zd, detectors.xd, tc, iu, ku, tx0, tz0);

					if (skip)
						continue;

					double apu1;
					if (attenuation_correction)
						LL = sqrt(x_diff * x_diff + z_diff * z_diff);

					//vector<double> templ_ijk;
					//vector<uint32_t> indices;

					//indices.reserve(Np * 2ULL);
					//templ_ijk.reserve(Np * 2ULL);

					for (size_t ii = 0; ii < static_cast<size_t>(Ny); ii++) {
						apu1 = (yy_vec[ii + 1ULL] - detectors.yd);
						if (apu1 > 0.) {
							tempj = static_cast<int32_t>(ii);
							break;
						}
					}

					double temp = 0.;
					uint32_t tempijk;
					if (crystal_size_z == 0.)
						tempijk = Nyx * static_cast<uint32_t>(tempk) + static_cast<uint32_t>(tempi);
					else {
						tempijk = static_cast<uint32_t>(tempi);
						const double temp_x = detectors.xs;
						detectors.xs = detectors.ys;
						detectors.ys = temp_x;
					}
					size_t koko = elements.size();

					for (uint32_t ii = 0; ii < Np; ii++) {

						if (tx0 < tz0) {

							if (attenuation_correction)
								compute_attenuation(tc, jelppi, LL, tx0, tempi, tempj, tempk, Nx, Nyx, atten);
							if (crystal_size_z == 0.) {
								orth_distance_full(tempj, Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, temp, tempijk, Nx, tempi, jelppi, 
									0., ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, indices, idx);
							}
							else {
								orth_distance_3D_full(tempj, Ny, Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, temp, tempijk, Nx, tempi, 
									tempk, jelppi, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 
									0, elements, indices, idx);
							}
							tempi += iu;
							if (iu > 0)
								tempijk++;
							else
								tempijk--;
							tx0 += txu;
							xyz = 1u;
						}
						else {
							if (attenuation_correction)
								compute_attenuation(tc, jelppi, LL, tz0, tempi, tempj, tempk, Nx, Nyx, atten);
							if (crystal_size_z == 0.) {
								orth_distance_full(tempj, Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, temp, tempijk, Nx, tempi, jelppi, 
									0., ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, indices, idx);
								if (ku > 0)
									tempijk += Ny_max;
								else
									tempijk -= Ny_max;
							}
							else if (ii == Np - 1u) {
								orth_distance_3D_full(tempj, Ny, Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, temp, tempijk, Nx, tempi, 
									tempk, jelppi, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 
									0, elements, indices, idx);
							}
							tempk += ku;
							tz0 += tzu;
							xyz = 3u;
						}
						if (tempk < 0 || tempi < 0 || tempi >= static_cast<int32_t>(Nx) || tempk >= static_cast<int32_t>(Nz)) {
							if (ii != Np - 1u && xyz == 3u && crystal_size_z != 0.) {
								tempk -= ku;
								orth_distance_3D_full(tempj, Ny, Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, temp, tempijk, Nx, tempi, 
									tempk, jelppi, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 
									0, elements, indices, idx);
							}
							break;
						}
					}

					temp = 1. / temp;

					if (attenuation_correction)
						temp *= jelppi;
					if (normalization)
						temp *= norm_coef[lo];

					for (size_t ii = koko; ii < koko + idx; ii++) {
						elements[ii] *= temp;
					}

					//for (uint32_t ii = 0; ii < templ_ijk.size(); ii++) {
					//	elements.emplace_back(templ_ijk[ii] * temp);
					//}

					lor[lo] = static_cast<uint16_t>(idx);
					lj++;

					continue;

				}
				continue;
			}
			else if (fabs(x_diff) < 1e-8) {
				if (detectors.xd <= maxxx && detectors.xd >= bx) {

					int32_t tempi = 0, tempk = 0, tempj = 0, ju = 0, ku = 0;
					double tyu = 0., tzu = 0., tc = 0., ty0 = 0., tz0 = 0.;

					const bool skip = siddon_pre_loop_2D(by, bz, y_diff, z_diff, maxyy, bzb, dy, dz, Ny, Nz, tempj, tempk, tyu, tzu, Np, TYPE, detectors.zs, 
						detectors.ys, detectors.zd, detectors.yd, tc, ju, ku, ty0, tz0);

					if (skip)
						continue;

					double apu1;

					//vector<double> templ_ijk;
					//vector<uint32_t> indices;

					//indices.reserve(Np * 2ULL);
					//templ_ijk.reserve(Np * 2ULL);

					for (size_t ii = 0ULL; ii < static_cast<size_t>(Nx); ii++) {
						apu1 = (xx_vec[ii + 1ULL] - detectors.xd);
						if (apu1 > 0.) {
							tempi = static_cast<int32_t>(ii);
							break;
						}
					}
					double temp = 0.;
					uint32_t tempijk;
					if (crystal_size_z == 0.)
						tempijk = Nyx * static_cast<uint32_t>(tempk) + static_cast<uint32_t>(tempj) * Nx;
					else
						tempijk = static_cast<uint32_t>(tempj) * Nx;

					if (attenuation_correction)
						LL = sqrt(y_diff * y_diff + z_diff * z_diff);

					size_t koko = elements.size();

					for (uint32_t ii = 0u; ii < Np; ii++) {

						if (ty0 < tz0) {
							if (attenuation_correction)
								compute_attenuation(tc, jelppi, LL, ty0, tempi, tempj, tempk, Nx, Nyx, atten);
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u, tempj, jelppi, 
									0., ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, indices, idx);
							}
							else {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u, tempj, 
									tempk, jelppi, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  ju, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 
									0, elements, indices, idx);
							}
							tempj += (ju);
							if (ju > 0)
								tempijk += Nx;
							else
								tempijk -= Nx;
							ty0 += tyu;
							xyz = 2u;
						}
						else {
							if (attenuation_correction)
								compute_attenuation(tc, jelppi, LL, tz0, tempi, tempj, tempk, Nx, Nyx, atten);
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u, tempj, jelppi, 
									0., ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, indices, idx);
								if (ku > 0)
									tempijk += Ny_max;
								else
									tempijk -= Ny_max;
							}
							else if (ii == Np - 1u) {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u, tempj, 
									tempk, jelppi, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  ju, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 
									0, elements, indices, idx);
							}
							tempk += ku;
							tz0 += tzu;
							xyz = 3u;
						}
						if (tempj < 0 || tempk < 0 || tempk >= static_cast<int32_t>(Nz) || tempj >= static_cast<int32_t>(Ny)) {
							if (ii != Np - 1u && xyz == 3u && crystal_size_z != 0.) {
								tempk -= ku;
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u, tempj, 
									tempk, jelppi, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  ju, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 
									0, elements, indices, idx);
							}
							break;
						}

					}

					temp = 1. / temp;

					if (attenuation_correction)
						temp *= jelppi;
					if (normalization)
						temp *= norm_coef[lo];

					for (size_t ii = koko; ii < koko + idx; ii++) {
						elements[ii] *= temp;
					}

					//for (uint32_t ii = 0; ii < templ_ijk.size(); ii++) {
					//	elements.emplace_back(templ_ijk[ii] * temp);
					//	//indices.emplace_back(indices[ii]);
					//}

					lor[lo] = static_cast<uint16_t>(idx);
					lj++;

					continue;
				}
				continue;
			}
			else {
				int32_t tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 1;
				double txu = 0., tyu = 0., tzu = 0., tc = 0., tx0 = 0., ty0 = 0., tz0 = 0.;
				const bool skip = siddon_pre_loop_3D(bx, by, bz, x_diff, y_diff, z_diff, maxxx, maxyy, bzb, dx, dy, dz, Nx, Ny, Nz, tempi, tempj, tempk, 
					tyu, txu, tzu, Np, TYPE, detectors, tc, iu, ju, ku, tx0, ty0, tz0);

				if (skip)
					continue;

				double temp = 0.;

				if (attenuation_correction)
					LL = sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
				uint32_t tempijk;
				if (crystal_size_z == 0.)
					tempijk = Nyx * static_cast<uint32_t>(tempk) + static_cast<uint32_t>(tempj) * Nx;
				else
					tempijk = static_cast<uint32_t>(tempj) * Nx;
				//vector<double> templ_ijk;
				//vector<uint32_t> indices;

				//indices.reserve(Np * 2ULL);
				//templ_ijk.reserve(Np * 2ULL);

				size_t koko = elements.size();


				for (uint32_t ii = 0u; ii < Np; ii++) {

					if (tz0 < ty0 && tz0 < tx0) {
						if (attenuation_correction)
							compute_attenuation(tc, jelppi, LL, tz0, tempi, tempj, tempk, Nx, Nyx, atten);
						if (crystal_size_z == 0.) {
							orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u, tempj, jelppi, 
								0., ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, indices, idx);
							if (ku > 0)
								tempijk += Nyx;
							else
								tempijk -= Nyx;
						}
						else if (ii == Np - 1u) {
							orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u, tempj, 
								tempk, jelppi, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 
								0, elements, indices, idx);
						}
						tempk += ku;
						tz0 += tzu;
						xyz = 3u;
					}
					else if (ty0 < tx0) {
						if (attenuation_correction)
							compute_attenuation(tc, jelppi, LL, ty0, tempi, tempj, tempk, Nx, Nyx, atten);
						if (crystal_size_z == 0.) {
							orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u, tempj, jelppi, 
								0., ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, indices, idx);
						}
						else {
							orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u, tempj, 
								tempk, jelppi, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 
								0, elements, indices, idx);
						}
						tempj += (ju);
						if (ju > 0)
							tempijk += Nx;
						else
							tempijk -= Nx;
						ty0 += tyu;
						xyz = 2u;
					}
					else {
						if (attenuation_correction)
							compute_attenuation(tc, jelppi, LL, tx0, tempi, tempj, tempk, Nx, Nyx, atten);
						if (ii == Np - 1u) {
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u, tempj, jelppi, 
									0., ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, indices, idx);
							}
							else {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u, tempj, 
									tempk, jelppi, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 
									0, elements, indices, idx);
							}
						}
						else {
							tempi += iu;
							tx0 += txu;
						}
						xyz = 1u;
					}
					if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= static_cast<int32_t>(Nx) || tempj >= static_cast<int32_t>(Ny) 
						|| tempk >= static_cast<int32_t>(Nz)) {
						if (ii != Np - 1u && (xyz == 1u || (xyz == 3u && crystal_size_z > 0.))) {
							if (xyz == 1u)
								tempi -= iu;
							else
								tempk -= ku;
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u, tempj, jelppi, 
									0., ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 0, elements, indices, idx);
							}
							else {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u, tempj, 
									tempk, jelppi, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, rhs, Summ, 
									0, elements, indices, idx);
							}
						}
						break;
					}
				}

				temp = 1. / temp;

				if (attenuation_correction)
					temp *= jelppi;
				if (normalization)
					temp *= norm_coef[lo];

				for (size_t ii = koko; ii < koko + idx; ii++) {
					elements[ii] *= temp;
				}

				//for (uint32_t ii = 0u; ii < templ_ijk.size(); ii++) {
				//	elements.emplace_back(templ_ijk[ii] * temp);
				//	//indices.emplace_back(indices[ii]);
				//}

				lor[lo] = static_cast<uint16_t>(idx);
				lj++;

				continue;
			}
		}
	}
	return lj;
}