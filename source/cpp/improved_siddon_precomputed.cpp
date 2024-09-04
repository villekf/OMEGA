/**************************************************************************
* Implements the improved Siddon's algorithm (Implementation 1).
* This version requires precomputation step; the number of voxels each LOR
* traverses needs to be known in advance.
* This version computes the system matrix column indices and elements for
* the preallocated MATLAB sparse matrix. Due to MATLAB's CSC format, this
* is a transposed version of the system matrix!
*
* Uses OpenMP for parallelization.
*
* Copyright (C) 2020-2024 Ville-Veikko Wettenhovi
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


void improved_siddon_precomputed(paramStruct<double>& param, const int64_t nMeas, const double* x, const double* z, double* elements, size_t* indices, const uint16_t* lor1, const uint64_t* lor2,
	const bool CT, const uint16_t* detIndex, const uint32_t nCores) {

#ifdef _OPENMP
	if (nCores == 0U)
		setThreads();
	else
		omp_set_num_threads(nCores);
#endif

	const uint32_t Nyx = param.Ny * param.Nx;

	// Precompute variables
	const double bmaxz = param.bz + static_cast<double>(param.Nz) * param.dz;
	const double bmaxy = static_cast<double>(param.Ny) * param.dy + param.by;
	const double bmaxx = static_cast<double>(param.Nx) * param.dx + param.bx;

	// Parallel for-loop
#ifdef _OPENMP
#if _OPENMP >= 201511 && defined(MATLAB)
#pragma omp parallel for schedule(monotonic:dynamic, nChunks)
#else
#pragma omp parallel for schedule(dynamic, nChunks)
#endif
#endif
	for (int64_t lo = 0LL; lo < nMeas; lo++) {
		int64_t ix = lo, iy = 0, iz = 0;
		if (param.subsets > 1 && param.subsetType == 1) {
			ix *= param.subsets;
			ix += param.currentSubset;
			iy = ix % param.size_y;
			iz = ix / ((int64_t)param.size_x * (int64_t)param.size_y);
			ix /= param.size_y;
			ix = ix % param.size_x;
		}
		else if (param.subsets > 1 && param.subsetType == 2) {
			ix *= param.subsets;
			ix += param.currentSubset;
			iz = ix / ((int64_t)param.size_x * (int64_t)param.size_y);
			iy = ix / param.size_x;
			iy = iy % param.size_y;
			ix = ix % param.size_x;
		}
		else if (param.subsets > 1 && param.subsetType == 4) {
			ix = (ix / param.size_x) * (int64_t)param.size_x * (int64_t)param.subsets + (int64_t)param.size_x * (int64_t)param.currentSubset + ix % param.size_x;
			iz = ix / ((int64_t)param.size_x * (int64_t)param.size_y);
			iy = (ix / param.size_x) % param.size_y;
			ix = ix % param.size_x;
		}
		else if (param.subsets > 1 && param.subsetType == 5) {
			iy = ix % param.size_y;
			ix = (ix / param.size_y * param.subsets + param.currentSubset);
			iz = ix / param.size_x;
			ix = ix % param.size_x;
		}
		else if (param.subsetType >= 8 || param.subsets == 1) {
			iz = lo / ((int64_t)param.size_x * (int64_t)param.size_y);
			ix = (lo - iz * (int64_t)param.size_x * (int64_t)param.size_y) % param.size_x;
			iy = (lo - iz * (int64_t)param.size_x * (int64_t)param.size_y) / param.size_x;
			iz += param.nMeas;
		}
		else {
			if (!CT) {
				ix = lo;
				iy = 0;
				iz = 0;
			}
			else {
				iz = ix / ((int64_t)param.size_x * (int64_t)param.size_y);
				iy = ix % param.size_y;
				ix = (ix / param.size_y) % param.size_x;
			}
		}
		Det<double> detectors;

		if (!CT) {
			get_detector_coordinates(x, z, param.size_x, param.size_y, detectors, param.xy_index, param.z_index, lo, param.subsetType, param.subsets, ix, iy, iz, param.nRays2D, param.nRays3D, 1, 1, param.dPitchZ, param.dPitchXY);
		}
		else
			get_detector_coordinates_CT(x, z, param.size_x, detectors, lo, param.subsets, param.size_y, ix, iy, iz, param.dPitchZ, param.nProjections, param.listMode, param.pitch);

		// Calculate the x, y and z distances of the detector pair
		double y_diff = (detectors.yd - detectors.ys);
		double x_diff = (detectors.xd - detectors.xs);
		double z_diff = (detectors.zd - detectors.zs);

		// Load the number of voxels the LOR traverses (precomputed)
		uint32_t Np = static_cast<uint32_t>(lor1[lo]);
		if (Np == 0U)
			continue;
		double L = norm(x_diff, y_diff, z_diff);
		double temp = 1.;
		// The initial index for the sparse matrix elements
		const uint64_t N2 = lor2[lo];
		const uint64_t N22 = lor2[lo + 1];
		uint32_t d_N0 = param.Nx;
		uint32_t d_N1 = param.Ny;
		uint32_t d_N2 = 1u;
		uint32_t d_N3 = param.Nx;
		uint32_t local_ind = 0u;
		int localIndX = 0, localIndY = 0, localIndZ = 0;
		bool XY = false;

		double local_norm = 0., local_scat = 0.;
		if (param.normalizationCorrection)
			local_norm = static_cast<double>(param.normCoefs[lo]);
		if (param.scatterCorrectionMult)
			local_scat = static_cast<double>(param.scatterCoefs[lo]);

		// If the measurement is on a same ring
		if (std::fabs(z_diff) < 1e-8 && (std::fabs(y_diff) < 1e-8 || std::fabs(x_diff) < 1e-8)) {

			// Z-coordinate (ring number)
			int32_t tempk = static_cast<int32_t>(std::fabs(detectors.zs - param.bz) / param.dz);
			double d_b, dd, d_db, d_d2;

			// If the LOR is perpendicular in the y-direction (Siddon cannot be used)
			if (std::fabs(y_diff) < 1e-8) {
				d_b = param.by;
				dd = detectors.yd;
				d_db = param.dy;
				d_d2 = param.dx;
				double xs_apu = detectors.xs;
				detectors.xs = detectors.ys;
				detectors.ys = xs_apu;
				double xdiff_apu = x_diff;
				x_diff = y_diff;
				y_diff = xdiff_apu;
				d_N0 = param.Ny;
				d_N1 = param.Nx;
				d_N2 = param.Ny;
				d_N3 = 1u;
			}
			else if (std::fabs(x_diff) < 1e-8) {
				d_b = param.bx;
				dd = detectors.xd;
				d_d2 = param.dy;
				d_db = param.dx;
			}
			temp = perpendicular_elements(d_N2, dd, d_d2, d_b, d_db, d_N0, d_N1, param.atten, local_norm, param.attenuationCorrection, param.normalizationCorrection, param.CTAttenuation, tempk, d_N3, param.globalFactor, param.scatterCorrectionMult,
				local_scat, localIndX, localIndY, localIndZ, L, 1, param.projType, param.Nx, param.Ny, CT, lo);
			local_ind = tempk;

			// Calculate the next index and store it as well as the probability of emission
			// In transmission tomography, store the length of intersection
			for (uint32_t ii = 0u; ii < d_N1; ii++) {
				indices[N2 + ii] = static_cast<size_t>(local_ind);
				elements[N2 + ii] = temp;
				local_ind += d_N3;
			}
		}
		else {
			int32_t tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;
			double txu = 0., tyu = 0., tzu = 0., tc = 0., tx0 = 1e8, ty0 = 1e8, tz0 = 1e8;
			bool skip = false;

			if (std::fabs(z_diff) < 1e-8) {
				tempk = static_cast<uint32_t>(std::fabs(detectors.zs - param.bz) / param.dz);
				skip = siddon_pre_loop_2D(param.bx, param.by, x_diff, y_diff, bmaxx, bmaxy, param.dx, param.dy, param.Nx, param.Ny, tempi, tempj, txu, tyu, Np, TYPE,
					detectors.ys, detectors.xs, detectors.yd, detectors.xd, tc, iu, ju, tx0, ty0, param.projType, XY);
			}
			//Detectors on different rings (e.g. oblique sinograms)
			else if (std::fabs(y_diff) < 1e-8) {
				tempj = perpendicular_start(param.by, detectors.yd, param.dy, param.Ny);
				skip = siddon_pre_loop_2D(param.bx, param.bz, x_diff, z_diff, bmaxx, bmaxz, param.dx, param.dz, param.Nx, param.Nz, tempi, tempk, txu, tzu, Np, TYPE,
					detectors.zs, detectors.xs, detectors.zd, detectors.xd, tc, iu, ku, tx0, tz0, param.projType, XY);
			}
			else if (std::fabs(x_diff) < 1e-8) {
				tempi = perpendicular_start(param.bx, detectors.xd, param.dx, param.Nx);
				skip = siddon_pre_loop_2D(param.by, param.bz, y_diff, z_diff, bmaxy, bmaxz, param.dy, param.dz, param.Ny, param.Nz, tempj, tempk, tyu, tzu, Np, TYPE,
					detectors.zs, detectors.ys, detectors.zd, detectors.yd, tc, ju, ku, ty0, tz0, param.projType, XY);
			}
			else {
				skip = siddon_pre_loop_3D(param.bx, param.by, param.bz, x_diff, y_diff, z_diff, bmaxx, bmaxy, bmaxz, param.dx, param.dy, param.dz, param.Nx, param.Ny, param.Nz, tempi, tempj, tempk, tyu, txu, tzu,
					Np, TYPE, detectors, tc, iu, ju, ku, tx0, ty0, tz0, param.projType, XY);
			}
			if (!CT) {
				temp = 1. / L;
				if (param.normalizationCorrection)
					temp *= local_norm;
				if (param.scatterCorrectionMult)
					temp *= local_scat;
				if (param.attenuationCorrection && !param.CTAttenuation)
					temp *= param.atten[lo];
				temp *= param.globalFactor;
			}
			double local_ele = 0., jelppi = 0.;

			// Compute the indices and matrix elements
			for (uint64_t ii = 0ULL; ii < Np; ii++) {
				local_ind = compute_ind(tempj, tempi, tempk, d_N3, Nyx);

				indices[N2 + ii] = static_cast<size_t>(local_ind);
				// Apply attenuation correction
				if (tz0 < ty0 && tz0 < tx0) {
					local_ele = compute_element(tz0, tc, L, tzu, ku, tempk);
				}
				else if (ty0 < tx0) {
					local_ele = compute_element(ty0, tc, L, tyu, ju, tempj);
				}
				else {
					local_ele = compute_element(tx0, tc, L, txu, iu, tempi);
				}
				if (param.attenuationCorrection && param.CTAttenuation)
					compute_attenuation(local_ele, local_ind, param.atten, jelppi);

				elements[N2 + ii] = local_ele;
			}

			if (!CT) {
				if (param.attenuationCorrection && param.CTAttenuation)
					temp *= std::exp(jelppi);

				for (uint32_t ii = 0u; ii < Np; ii++) {
					elements[N2 + ii] *= temp;
				}
			}
		}
	}
}
