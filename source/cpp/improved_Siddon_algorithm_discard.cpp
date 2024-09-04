/**************************************************************************
* This function is used to check for the number of voxels each LOR/ray
* traverses and also whether the LOR/ray actually intersects with the pixel
* space.
* Output is the number of voxels the ray has traversed (if the LOR does not
* traverse the pixel space, this value will be 0).
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

const static int TYPE = 0;

void improved_siddon_precomputation_phase(paramStruct<double>& param, const int64_t nMeas, const double* x, const double* z, uint16_t* lor, const bool CT, 
	const uint16_t* detIndex, const uint32_t nCores) {

#ifdef _OPENMP
	if (nCores == 0U)
		setThreads();
	else
		omp_set_num_threads(nCores);
#endif

	const double bmaxz = param.bz + static_cast<double>(param.Nz) * param.dz;
	const double bmaxy = static_cast<double>(param.Ny) * param.dy + param.by;
	const double bmaxx = static_cast<double>(param.Nx) * param.dx + param.bx;
	const uint32_t Nyx = param.Nx * param.Ny;

	param.subsets = 1;

#ifdef _OPENMP
#if _OPENMP >= 201511 && defined(MATLAB)
#pragma omp parallel for schedule(monotonic:dynamic, nChunks)
#else
#pragma omp parallel for schedule(dynamic, nChunks)
#endif
#endif
	for (int64_t lo = 0LL; lo < nMeas; lo++) {
		int64_t ix = lo, iy = 0, iz = 0;
		iz = lo / ((int64_t)param.size_x * (int64_t)param.size_y);
		ix = (lo - iz * (int64_t)param.size_x * (int64_t)param.size_y) % (int64_t)param.size_x;
		iy = (lo - iz * (int64_t)param.size_x * (int64_t)param.size_y) / (int64_t)param.size_x;
		Det<double> detectors;

		if (!CT) {
			get_detector_coordinates(x, z, param.size_x, param.size_y, detectors, param.xy_index, param.z_index, lo, param.subsetType, param.subsets, ix, iy, iz, param.nRays2D, param.nRays3D, 1, 1, param.dPitchZ, param.dPitchXY);
		}
		else
			get_detector_coordinates_CT(x, z, param.size_x, detectors, lo, param.subsets, param.size_y, ix, iy, iz, param.dPitchZ, param.nProjections, param.listMode, param.pitch);


		const double x_diff = (detectors.xd - detectors.xs);
		const double y_diff = (detectors.yd - detectors.ys);
		const double z_diff = (detectors.zd - detectors.zs);
		if ((y_diff == 0. && x_diff == 0. && z_diff == 0.) || (y_diff == 0. && x_diff == 0.) || std::isinf(y_diff) || std::isinf(x_diff))
			continue;
		bool XY = false;
		uint16_t temp_koko = 0u;

		uint32_t Np = 0u;

		if (std::fabs(z_diff) < 1e-8 && (std::fabs(y_diff) < 1e-8 || std::fabs(x_diff) < 1e-8)) {

			const int32_t tempk = static_cast<int32_t>(fabs(detectors.zs - param.bz) / param.dz);
			if (tempk >= param.Nz || tempk < 0)
				continue;

			if (std::fabs(y_diff) < 1e-8) {

				if (detectors.yd <= bmaxy && detectors.yd >= param.by) {
					// The number of voxels the LOR/ray traverses
					lor[lo] = static_cast<uint16_t>(param.Nx);
				}
				// LOR/ray doesn't traverse through the pixel space
			}
			else if (std::fabs(x_diff) < 1e-8) {

				if (detectors.xd <= bmaxx && detectors.xd >= param.bx) {
					// The number of voxels the LOR/ray traverses
					lor[lo] = static_cast<uint16_t>(param.Ny);
				}
				// LOR/ray doesn't traverse through the pixel space
			}
		}
		else {
			int32_t tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;
			double txu = 0., tyu = 0., tzu = 0., tc = 0., tx0 = 1e8, ty0 = 1e8, tz0 = 1e8;
			bool skip = false;

			if (std::fabs(z_diff) < 1e-8) {
				tempk = static_cast<int32_t>(std::fabs(detectors.zs - param.bz) / param.dz);
				if (tempk >= param.Nz || tempk < 0)
					skip = true;
				skip = siddon_pre_loop_2D(param.bx, param.by, x_diff, y_diff, bmaxx, bmaxy, param.dx, param.dy, param.Nx, param.Ny, tempi, tempj, txu, tyu, Np, TYPE,
					detectors.ys, detectors.xs, detectors.yd, detectors.xd, tc, iu, ju, tx0, ty0, param.projType, XY);
			}
			//Detectors on different rings (e.g. oblique sinograms)
			else if (std::fabs(y_diff) < 1e-8) {
				if (detectors.yd > bmaxy || detectors.yd < param.by)
					skip = true;
				tempj = perpendicular_start(param.by, detectors.yd, param.dy, param.Ny);
				skip = siddon_pre_loop_2D(param.bx, param.bz, x_diff, z_diff, bmaxx, bmaxz, param.dx, param.dz, param.Nx, param.Nz, tempi, tempk, txu, tzu, Np, TYPE,
					detectors.zs, detectors.xs, detectors.zd, detectors.xd, tc, iu, ku, tx0, tz0, param.projType, XY);
			}
			else if (std::fabs(x_diff) < 1e-8) {
				if (detectors.xd > bmaxx || detectors.xd < param.bx)
					skip = true;
				tempi = perpendicular_start(param.bx, detectors.xd, param.dx, param.Nx);
				skip = siddon_pre_loop_2D(param.by, param.bz, y_diff, z_diff, bmaxy, bmaxz, param.dy, param.dz, param.Ny, param.Nz, tempj, tempk, tyu, tzu, Np, TYPE,
					detectors.zs, detectors.ys, detectors.zd, detectors.yd, tc, ju, ku, ty0, tz0, param.projType, XY);
			}
			else {
				skip = siddon_pre_loop_3D(param.bx, param.by, param.bz, x_diff, y_diff, z_diff, bmaxx, bmaxy, bmaxz, param.dx, param.dy, param.dz, param.Nx, param.Ny, param.Nz, tempi, tempj, tempk, tyu, txu, tzu,
					Np, TYPE, detectors, tc, iu, ju, ku, tx0, ty0, tz0, param.projType, XY);
			}

			if (!skip) {

				for (uint32_t ii = 0u; ii < Np; ii++) {

					temp_koko++;

					if (tx0 < ty0 && tx0 < tz0) {
						tempi += iu;
						tx0 += txu;
					}
					else if (ty0 < tz0) {
						tempj += ju;
						ty0 += tyu;
					}
					else {
						tempk += ku;
						tz0 += tzu;
					}
					if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= static_cast<int32_t>(param.Nx) || tempj >= static_cast<int32_t>(param.Ny) || tempk >= static_cast<int32_t>(param.Nz)) {
						break;
					}
				}
				// The number of voxels the LOR/ray traverses
				lor[lo] = temp_koko;
			}
		}
	}
	
}