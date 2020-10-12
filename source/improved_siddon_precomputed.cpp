/**************************************************************************
* Implements the improved Siddon's algorithm (Implementation 1).
* This version requires precomputation step; the number of voxels each LOR
* traverses needs to be known in advance.
* This version computes the system matrix column indices and elements for
* the preallocated MATLAB sparse matrix. Due to MATLAB's CSC format, this
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

// if 0, then determines whether the LOR intercepts the FOV (i.e. no precomputation phase performed)
const static int TYPE = 1;

using namespace std;


void improved_siddon_precomputed(const int64_t loop_var_par, const uint32_t size_x, const double zmax, size_t* indices, double* elements, const double maxyy,
	const double maxxx, const vector<double>& xx_vec, const double dy, const vector<double>& yy_vec, const double* atten, const double* norm_coef, 
	const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const double dx, 
	const double dz, const double bx, const double by, const double bz, const bool attenuation_correction, const bool normalization, const uint16_t* lor1, 
	const uint64_t* lor2, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos, const uint16_t* L, const uint32_t* pseudos, 
	const uint32_t pRows, const uint32_t det_per_ring, const bool raw, const bool attenuation_phase, double* length, const double global_factor, 
	const bool scatter, const double* scatter_coef) {

	setThreads();

	const uint32_t Nyx = Ny * Nx;

	// Precompute variables
	const double bzb = bz + static_cast<double>(Nz) * dz;

	// Parallel for-loop
#pragma omp parallel for schedule(dynamic)
	for (int64_t lo = 0LL; lo < loop_var_par; lo++) {

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
		// The initial index for the sparse matrix elements
		const uint64_t N2 = lor2[lo];
		const uint64_t N22 = lor2[lo + 1];

		// If the measurement is on a same ring
		if (fabs(z_diff) < 1e-8 && (fabs(y_diff) < 1e-8 || fabs(x_diff) < 1e-8)) {

			// Z-coordinate (ring number)
			const uint32_t tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));

			// If the LOR is perpendicular in the y-direction (Siddon cannot be used)
			if (fabs(y_diff) < 1e-8) {

				if (detectors.yd <= maxyy && detectors.yd >= by) {
					uint32_t temp_ijk = 0u;

					const double element = perpendicular_elements(Ny, detectors.yd, yy_vec, dx, tempk, Nx, Ny, atten, norm_coef, attenuation_correction,
						normalization, temp_ijk, 1u, lo, global_factor, scatter, scatter_coef);

					// Calculate the next index and store it as well as the probability of emission
					for (uint64_t ii = 0ULL; ii < static_cast<uint64_t>(Nx); ii++) {
						indices[N2 + ii] = static_cast<size_t>(temp_ijk) + static_cast<size_t>(ii);
						elements[N2 + ii] = fabs(element);
					}
				}
			}
			// Same as for the y-case above
			else if (fabs(x_diff) < 1e-8) {

				if (detectors.xd <= maxxx && detectors.xd >= bx) {
					uint32_t temp_ijk = 0u;

					const double element = perpendicular_elements(1u, detectors.xd, xx_vec, dy, tempk, Ny, Nx, atten, norm_coef, attenuation_correction,
						normalization, temp_ijk, Nx, lo, global_factor, scatter, scatter_coef);

					for (uint64_t ii = 0u; ii < static_cast<uint64_t>(Ny); ii++) {
						indices[N2 + ii] = static_cast<size_t>(temp_ijk) + static_cast<size_t>(ii) * static_cast<size_t>(Nx);
						elements[N2 + ii] = fabs(element);
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

			// d_conv
			const double LL = sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);

			double temp = 0., apu = 0.;

			uint32_t tempijk = static_cast<uint32_t>(tempk) * Nyx + static_cast<uint32_t>(tempj) * Nx + static_cast<uint32_t>(tempi);

			// Compute the indices and matrix elements
			for (uint32_t ii = 0u; ii < Np; ii++) {

				indices[N2 + ii] = static_cast<size_t>(tempijk);
				if (tx0 < ty0 && tx0 < tz0) {
					apu = pixel_value(tx0, tc, LL);
					tempi += iu;
					if (iu > 0)
						tempijk++;
					else
						tempijk--;
					tc = tx0;
					tx0 += txu;
				}
				else if (ty0 < tz0) {
					apu = pixel_value(ty0, tc, LL);
					tempj += ju;
					if (ju > 0)
						tempijk += Nx;
					else
						tempijk -= Nx;
					tc = ty0;
					ty0 += tyu;
				}
				else {
					apu = pixel_value(tz0, tc, LL);
					tempk += ku;
					if (ku > 0)
						tempijk += Nyx;
					else
						tempijk -= Nyx;
					tc = tz0;
					tz0 += tzu;
				}

				temp += apu;
				elements[N2 + ii] = (apu);
			}

			// If computing the attenuation image (Inveon)
			if (attenuation_phase)
				length[lo] = temp;

			temp = 1. / temp;

			// Apply attenuation correction
			if (attenuation_correction)
				att_corr_vec_precomp(elements, atten, indices, Np, N2, temp);
			if (normalization)
				temp *= norm_coef[lo];
			if (scatter)
				temp *= scatter_coef[lo];
			temp *= global_factor;


			for (uint32_t ii = 0u; ii < Np; ii++) {
				elements[N2 + ii] = fabs(elements[N2 + ii]) * (temp);
			}
		}
	}
}
