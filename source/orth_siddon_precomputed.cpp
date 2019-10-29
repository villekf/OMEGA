/**************************************************************************
* Implements the Orthogonal Siddon's algorithm (Implementation 1).
* This version requires precomputation step; the number of voxels each LOR
* traverses needs to be known in advance.
* This version computes the system matrix column indices and elements for
* the preallocated MATLAB sparse matrix. Due to MATLAB's CSR format, this
* is essentially a transposed version of the system matrix.
*
* Uses C++11 threads for parallelization.
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

// if 0, then determines whether the LOR intercepts the FOV
constexpr int TYPE = 1;

// Whether to use the OpenMP code or not
constexpr bool OMP = false;

// Using non-OpenMP with either precomputation or without
constexpr bool PRECOMPUTE = true;

using namespace std;

void orth_siddon_precomputed(const size_t loop_var_par, const uint32_t size_x, const double zmax, mwIndex* indices, double* elements, const double maxyy,
	const double maxxx, const vector<double>& xx_vec, const double dy, const vector<double>& yy_vec, const double* atten, const double* norm_coef, 
	const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const double dx, 
	const double dz, const double bx, const double by, const double bz, const bool attenuation_correction, const bool normalization, const uint16_t* lor1, 
	const uint64_t* lor2, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos, const uint16_t* L, const uint32_t* pseudos, 
	const uint32_t pRows, const uint32_t det_per_ring, const bool raw, const bool attenuation_phase, double* length, const double crystal_size, 
	const double crystal_size_z, const double* y_center, const double* x_center, const double* z_center, const int32_t dec_v) {


	// Precompute
	const double bzb = bz + static_cast<double>(Nz) * dz;
	const uint32_t Nyx = Ny * Nx;

	double ax = 0.;
	double* osem_apu = nullptr;
	double* Summ = nullptr;
	const bool no_norm = true;
	vector<double> v_elements;
	vector<uint32_t> v_indices;

	const int32_T dec = static_cast<int32_T>(ceil(crystal_size_z / sqrt(dz * dz * 2.))) * dec_v;


	ThreadPool::ParallelFor(static_cast<size_t>(0), loop_var_par, [&](uint32_t lo) {

		Det detectors;
		double kerroin, length_, jelppi = 0., LL;

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
		// The initial index for the sparse matrix elements
		const uint64_t N2 = lor2[lo];

		size_t idx = 0ULL;

		// Precompute constants
		if (crystal_size_z == 0.) {
			kerroin = detectors.xd * detectors.ys - detectors.yd * detectors.xs;
			length_ = sqrt(y_diff * y_diff + x_diff * x_diff) * crystal_size;
		}
		else {
			kerroin = norm(x_diff, y_diff, z_diff) * crystal_size_z;
		}

		// If the measurement is on a same ring
		if (fabs(z_diff) < 1e-8) {

			// Z-coordinate (ring number)
			const uint32_t tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));

			// If the LOR is perpendicular in the y-direction (Siddon cannot be used)
			if (fabs(y_diff) < 1e-8) {

				if (detectors.yd <= maxyy && detectors.yd >= by) {
					uint32_t temp_ijk = 0u;
					double temp = 0.;
					int hpk = 0;
					if (crystal_size_z == 0.) {
						std::vector<double> vec1;
						std::vector<double> vec2;

						orth_perpendicular(Ny, detectors.yd, yy_vec, dx, tempk, Nx, Ny, -x_diff, y_center, kerroin, length_, vec1, vec2, temp, temp_ijk);

						if (attenuation_correction)
							att_corr_scalar_orth(temp_ijk, atten, temp, 1u, Nx, dx);
						if (normalization)
							temp *= norm_coef[lo];

						// Calculate the next index and store it as well as the probability of emission
						for (int32_t jj = -static_cast<int32_t>(vec1.size()) + 1; jj <= static_cast<int32_t>(vec2.size()); jj++) {
							double element;
							if (jj <= 0) {
								element = vec1[-jj];
							}
							else {
								element = vec2[static_cast<size_t>(jj) - 1ULL];
							}
							element *= temp;
							for (int64_t ii = 0LL; ii < static_cast<int64_t>(Ny); ii++) {
								indices[N2 + hpk] = static_cast<mwIndex>(static_cast<int64_t>(temp_ijk) + ii + static_cast<int64_t>(jj) * 
									static_cast<int64_t>(Ny));
								elements[N2 + hpk] = element;
								hpk++;
							}
						}
					}
					else {
						double temppi = detectors.xs;
						detectors.xs = detectors.ys;
						detectors.ys = temppi;
						orth_perpendicular_3D(detectors.yd, yy_vec, tempk, Nx, Ny, Nz, Nyx, Ny, 1u, detectors, y_diff, x_diff, z_diff, y_center, x_center[0],
							z_center, kerroin, hpk, temp, temp_ijk, indices, elements, N2);

						if (attenuation_correction)
							att_corr_scalar_orth(temp_ijk, atten, temp, 1u, Nx, dx);
						if (normalization)
							temp *= norm_coef[lo];

						for (size_t jj = 0ULL; jj < static_cast<size_t>(hpk); jj++) {
							elements[N2 + jj] *= temp;
						}
					}
				}
			}
			// Same as for the y-case above
			else if (fabs(x_diff) < 1e-8) {

				if (detectors.xd <= maxxx && detectors.xd >= bx) {
					uint32_t temp_ijk = 0u;
					double temp = 0.;
					int hpk = 0;
					if (crystal_size_z == 0.) {
						std::vector<double> vec1;
						std::vector<double> vec2;

						orth_perpendicular(1u, detectors.xd, xx_vec, dy, tempk, Ny, Nx, y_diff, x_center, kerroin, length_, vec1, vec2, temp, temp_ijk);


						if (attenuation_correction)
							att_corr_scalar_orth(temp_ijk, atten, temp, Ny, 1u, dy);
						if (normalization)
							temp *= norm_coef[lo];

						// Calculate the next index and store it as well as the probability of emission
						for (int64_t ii = 0LL; ii < static_cast<int64_t>(Nx); ii++) {
							for (int32_t jj = -static_cast<int32_t>(vec1.size()) + 1; jj <= static_cast<int32_t>(vec2.size()); jj++) {
								double element;
								if (jj <= 0)
									element = vec1[-jj];
								else
									element = vec2[static_cast<size_t>(jj) - 1ULL];
								element *= temp;
								indices[N2 + hpk] = static_cast<mwIndex>(static_cast<int64_t>(temp_ijk) + static_cast<int64_t>(ii) * 
									static_cast<int64_t>(Ny) + jj);
								elements[N2 + hpk] = element;
								hpk++;
							}
						}
					}
					else {
						orth_perpendicular_3D(detectors.xd, xx_vec, tempk, Ny, Nx, Nz, Nyx, 1u, Nx, detectors, x_diff, y_diff, z_diff, x_center, y_center[0],
							z_center, kerroin, hpk, temp, temp_ijk, indices, elements, N2);

						if (attenuation_correction)
							att_corr_scalar_orth(temp_ijk, atten, temp, Ny, 1u, dy);
						if (normalization)
							temp *= norm_coef[lo];

						for (int32_t jj = 0; jj < hpk; jj++) {
							elements[N2 + jj] *= temp;
						}
					}
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
				uint32_t tempijk;
				if (crystal_size_z == 0.)
					tempijk = Nyx * tempk + static_cast<uint32_t>(tempj) * Nx;
				else
					tempijk = static_cast<uint32_t>(tempj) * Nx;

				// Compute the indices and matrix elements
				for (uint32_t ii = 0u; ii < Np; ii++) {

					// Ray goes along the x-axis
					if (tx0 < ty0) {
						if (attenuation_correction)
							compute_attenuation(tc, jelppi, LL, tx0, tempi, tempj, tempk, Nx, Nyx, atten);
						if (ii == Np - 1u) {
							// 2D case
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u, tempj, 
									0., ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, elements, Summ, indices, v_elements, v_indices, idx, N2);
							}
							// 3D case
							else {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u, tempj, 
									tempk, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, elements, Summ, 
									indices, v_elements, v_indices, idx, N2);
							}
						}
						else {
							// (32)
							tempi += iu;
							// (34)
							tx0 += txu;
						}
					}
					// Ray goes along the y-axis
					else {
						if (attenuation_correction)
							compute_attenuation(tc, jelppi, LL, ty0, tempi, tempj, tempk, Nx, Nyx, atten);
						if (crystal_size_z == 0.) {
							orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u, tempj, 
								0., ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, elements, Summ, indices, v_elements, v_indices, idx, N2);
						}
						else {
							orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u, tempj, tempk, 
								0., ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, elements, Summ, indices, 
								v_elements, v_indices, idx, N2);
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

				if (attenuation_phase)
					length[lo] = temp;

				temp = 1. / temp;

				if (attenuation_correction)
					temp *= exp(jelppi);
				if (normalization)
					temp *= norm_coef[lo];
				for (size_t ii = 0u; ii < idx; ii++) {
					elements[N2 + ii] *= temp;
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
									tempi, 0., ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, elements, Summ, indices, v_elements, v_indices, 
									idx, N2);
							}
							else {
								orth_distance_3D_full(tempj, Ny, Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, temp, tempijk, Nx,
									tempi, tempk, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, 
									elements, Summ, indices, v_elements, v_indices, idx, N2);
							}
							tempi += iu;
							if (iu > 0)
								tempijk++;
							else
								tempijk--;
							tx0 += txu;
						}
						else {
							if (attenuation_correction)
								compute_attenuation(tc, jelppi, LL, tz0, tempi, tempj, tempk, Nx, Nyx, atten);
							if (crystal_size_z == 0.) {
								orth_distance_full(tempj, Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, temp, tempijk, Nx,
									tempi, 0., ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, elements, Summ, indices, v_elements, v_indices, 
									idx, N2);
								if (ku > 0)
									tempijk += Nyx;
								else
									tempijk -= Nyx;
							}
							else if (ii == Np - 1u) {
								orth_distance_3D_full(tempj, Ny, Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, temp, tempijk, Nx,
									tempi, tempk, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, 
									elements, Summ, indices, v_elements, v_indices, idx, N2);
							}
							tempk += ku;
							tz0 += tzu;
						}
					}

					if (attenuation_phase)
						length[lo] = temp;

					temp = 1. / temp;

					if (attenuation_correction)
						temp *= exp(jelppi);
					if (normalization)
						temp *= norm_coef[lo];

					for (size_t ii = 0u; ii < idx; ii++) {
						elements[N2 + ii] *= temp;
					}

				}
			}
			else if (fabs(x_diff) < 1e-8) {
				if (detectors.xd <= maxxx && detectors.xd >= bx) {

					int32_t tempi = 0, tempk = 0, tempj = 0, ju = 0, ku = 0;
					double tyu = 0., tzu = 0., tc = 0., ty0 = 0., tz0 = 0.;
					const bool skip = siddon_pre_loop_2D(by, bz, y_diff, z_diff, maxyy, bzb, dy, dz, Ny, Nz, tempj, tempk, tyu, tzu, Np, TYPE,
						detectors.zs, detectors.ys, detectors.zd, detectors.yd, tc, ju, ku, ty0, tz0);

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
									tempj, 0., ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, elements, Summ, indices, v_elements, v_indices, 
									idx, N2);
							}
							else {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
									tempj, tempk, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  ju, no_norm, false, false, OMP, PRECOMPUTE, 
									elements, Summ, indices, v_elements, v_indices, idx, N2);
							}
							tempj += ju;
							if (ju > 0) {
								tempijk += Nx;
							}
							else {
								tempijk -= Nx;
							}
							ty0 += tyu;
						}
						else {
							if (attenuation_correction)
								compute_attenuation(tc, jelppi, LL, tz0, tempi, tempj, tempk, Nx, Nyx, atten);
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
									tempj, 0., ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, elements, Summ, indices, v_elements, v_indices, 
									idx, N2);
								if (ku > 0)
									tempijk += Nyx;
								else
									tempijk -= Nyx;
							}
							else if (ii == Np - 1u) {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
									tempj, tempk, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  ju, no_norm, false, false, OMP, PRECOMPUTE, 
									elements, Summ, indices, v_elements, v_indices, idx, N2);
							}
							tempk += ku;
							tz0 += tzu;
						}
					}

					if (attenuation_phase)
						length[lo] = temp;

					temp = 1. / temp;

					if (attenuation_correction)
						temp *= exp(jelppi);
					if (normalization)
						temp *= norm_coef[lo];

					for (size_t ii = 0u; ii < idx; ii++) {
						elements[N2 + ii] *= temp;
					}

				}
			}
			else {

				int32_t tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 1;
				double txu = 0., tyu = 0., tzu = 0., tc = 0., tx0 = 0., ty0 = 0., tz0 = 0.;
				const bool skip = siddon_pre_loop_3D(bx, by, bz, x_diff, y_diff, z_diff, maxxx, maxyy, bzb, dx, dy, dz, Nx, Ny, Nz, tempi, tempj, tempk, 
					tyu, txu, tzu,
					Np, TYPE, detectors, tc, iu, ju, ku, tx0, ty0, tz0);

				double temp = 0.;

				if (attenuation_correction)
					LL = sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
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
								tempj, 0., ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, elements, Summ, indices, v_elements, v_indices, 
								idx, N2);
							if (ku > 0)
								tempijk += Nyx;
							else
								tempijk -= Nyx;
						}
						else if (ii == Np - 1u) {
							orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
								tempj, tempk, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, elements, 
								Summ, indices, v_elements, v_indices, idx, N2);
						}
						tempk += ku;
						tz0 += tzu;
					}
					else if (ty0 < tx0) {
						if (attenuation_correction)
							compute_attenuation(tc, jelppi, LL, ty0, tempi, tempj, tempk, Nx, Nyx, atten);
						if (crystal_size_z == 0.) {
							orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
								tempj, 0., ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, elements, Summ, indices, v_elements, v_indices, 
								idx, N2);
						}
						else {
							orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
								tempj, tempk, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, 
								elements, Summ, indices, v_elements, v_indices, idx, N2);
						}
						if (ju > 0)
							tempijk += Nx;
						else
							tempijk -= Nx;

						tempj += ju;
						ty0 += tyu;
					}
					else {
						if (attenuation_correction)
							compute_attenuation(tc, jelppi, LL, tx0, tempi, tempj, tempk, Nx, Nyx, atten);
						if (ii == Np - 1u) {
							if (crystal_size_z == 0.) {
								orth_distance_full(tempi, Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, tempijk, 1u,
									tempj, 0., ax, osem_apu, no_norm, false, false, OMP, PRECOMPUTE, elements, Summ, indices, v_elements, v_indices, 
									idx, N2);
							}
							else {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, temp, tempijk, 1u,
									tempj, tempk, 0., ax, osem_apu, detectors, Nyx, kerroin, dec,  iu, no_norm, false, false, OMP, PRECOMPUTE, 
									elements, Summ, indices, v_elements, v_indices, idx, N2);
							}
						}
						else {
							tempi += iu;
							tx0 += txu;
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

				for (size_t ii = 0u; ii < idx; ii++) {
					elements[N2 + ii] *= temp;
				}
			}
		}
	});
}
