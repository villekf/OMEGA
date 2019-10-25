/**************************************************************************
* Compute the number of voxels that each LOR traverses, by using the
* improved Siddon's algorithm. OpenCL version (single, f32, precision).
*
* Used by implementations 2 and 3.
*
* Copyright (C) 2019  Ville-Veikko Wettenhovi
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
#include "opencl_functions.h"
#define TYPE 0

// Matrix free Improved Siddon's algorithm
__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(64, 1, 1)))
void siddon_precomp(const uint d_Nxy, const uint d_N, const uint d_Nx, const uint d_Ny, const uint d_Nz,
	const float d_dz, const float d_dx, const float d_dy, const float d_bz, const float d_bx, const float d_by, 
	const float d_bzb, const float d_maxxx,	const float d_maxyy, const float d_zmax, const float d_NSlices, const uint d_size_x, 
	const ushort d_TotSinos, const uint d_det_per_ring, const uchar d_raw, const uint d_pRows, 
	__constant uint* d_pseudos, const __global float* d_x, 
	const __global float* d_y, const __global float* d_zdet, __global ushort* d_lor, const __global ushort* d_L, const ulong m_size) {
	// Get the current global index
	uint idx = get_global_id(0);
	if (idx >= m_size)
		return;
	float xs, xd, ys, yd, zs, zd;
	// Load the next detector index
	// raw list-mode data
	if (d_raw) {
		get_detector_coordinates_raw(d_x, d_y, d_zdet, d_L, d_det_per_ring, idx, d_pseudos, d_pRows, &xs, &xd, &ys, &yd, &zs, &zd);
	}
	// Sinogram data
	else {
		get_detector_coordinates_precomp(d_size_x, idx, d_TotSinos, &xs, &xd, &ys, &yd, &zs, &zd, d_x, d_y, d_zdet);
	}
	// Calculate the x, y and z distances of the detector pair
	const float y_diff = (yd - ys);
	const float x_diff = (xd - xs);
	const float z_diff = (zd - zs);
	uint Np = 0u;
	if ((y_diff == 0.f && x_diff == 0.f && z_diff == 0.f) || (y_diff == 0.f && x_diff == 0.f))
		return;
	// Load the number of voxels the LOR traverses (precomputed)
	ushort temp_koko = 0u;
	// If the measurement is on a same ring
	if (fabs(z_diff) < 1e-6f) {

		// If the LOR is perpendicular in the y-direction
		if (fabs(y_diff) < 1e-6f) {
			if (yd <= d_maxyy && yd >= d_by) {
				d_lor[idx] = convert_ushort(d_Nx);
			}
		}
		else if (fabs(x_diff) < 1e-6f) {
			if (xd <= d_maxxx && xd >= d_bx) {
				d_lor[idx] = convert_ushort(d_Ny);
			}
		}
		else {
			int tempi = 0, tempj = 0, iu = 0, ju = 0;
			float txu = 0.f, tyu = 0.f, tc = 0.f, tx0 = 0.f, ty0 = 0.f;
			const bool skip = siddon_pre_loop_2D(d_bx, d_by, x_diff, y_diff, d_maxxx, d_maxyy, d_dx, d_dy, d_Nx, d_Ny, &tempi, &tempj, &txu, &tyu, &Np, TYPE,
				ys, xs, yd, xd, &tc, &iu, &ju, &tx0, &ty0);
			if (tempi < 0 || tempj < 0 || tempi >= d_Nx || tempj >= d_Ny)
				return;
			for (uint ii = 0u; ii < Np; ii++) {
				temp_koko++;
				if (tx0 < ty0) {
					tempi += iu;
					tx0 += txu;
				}
				else {
					tempj += ju;
					ty0 += tyu;
				}
				if (tempj < 0 || tempi < 0 || tempi >= d_Nx || tempj >= d_Ny)
					break;
			}
			d_lor[idx] = temp_koko;
		}
	}
	else {
		if (fabs(y_diff) < 1e-6f) {
			if (yd <= d_maxyy && yd >= d_by) {
				int tempi = 0, tempk = 0, tempj = 0, iu = 0, ku = 0;
				float txu = 0.f, tzu = 0.f, tc = 0.f, tx0 = 0.f, tz0 = 0.f;
				const bool skip = siddon_pre_loop_2D(d_bx, d_bz, x_diff, z_diff, d_maxxx, d_bzb, d_dx, d_dz, d_Nx, d_Nz, &tempi, &tempk, &txu, &tzu, &Np, TYPE,
					zs, xs, zd, xd, &tc, &iu, &ku, &tx0, &tz0);
				if (tempi < 0 || tempk < 0 || tempi >= d_Nx || tempk >= d_Nz)
					return;
				for (uint ii = 0u; ii < Np; ii++) {
					temp_koko++;
					if (tx0 < tz0) {
						tempi += iu;
						tx0 += txu;
					}
					else {
						tempk += ku;
						tz0 += tzu;
					}
					if (tempk < 0 || tempi < 0 || tempi >= d_Nx || tempk >= d_Nz)
						break;
				}
				d_lor[idx] = temp_koko;
			}
		}
		else if (fabs(x_diff) < 1e-6f) {
			if (xd <= d_maxxx && xd >= d_bx) {
				int tempi = 0, tempk = 0, tempj = 0, ju = 0, ku = 0;
				float tyu = 0.f, tzu = 0.f, tc = 0.f, ty0 = 0.f, tz0 = 0.f;
				const bool skip = siddon_pre_loop_2D(d_by, d_bz, y_diff, z_diff, d_maxyy, d_bzb, d_dy, d_dz, d_Ny, d_Nz, &tempj, &tempk, &tyu, &tzu, &Np, TYPE,
					zs, ys, zd, yd, &tc, &ju, &ku, &ty0, &tz0);
				if (tempk < 0 || tempj < 0 || tempk >= d_Nz || tempj >= d_Ny)
					return;
				for (uint ii = 0u; ii < Np; ii++) {
					temp_koko++;
					if (tz0 < ty0) {
						tempk += ku;
						tz0 += tzu;
					}
					else {
						tempj += ju;
						ty0 += tyu;
					}
					if (tempj < 0 || tempk < 0 || tempk >= d_Nz || tempj >= d_Ny)
						break;
				}
				d_lor[idx] = temp_koko;
			}
		}
		else {
			int tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;
			float txu = 0.f, tyu = 0.f, tzu = 0.f, tc = 0.f, tx0 = 0.f, ty0 = 0.f, tz0 = 0.f;
			const bool skip = siddon_pre_loop_3D(d_bx, d_by, d_bz, x_diff, y_diff, z_diff, d_maxxx, d_maxyy, d_bzb, d_dx, d_dy, d_dz, d_Nx, d_Ny, d_Nz, &tempi, &tempj, &tempk, &tyu, &txu, &tzu,
				&Np, TYPE, ys, xs, yd, xd, zs, zd, &tc, &iu, &ju, &ku, &tx0, &ty0, &tz0);
			if (tempi < 0 || tempj < 0 || tempk < 0 || tempi >= d_Nx || tempj >= d_Ny || tempk >= d_Nz)
				return;
			for (uint ii = 0u; ii < Np; ii++) {
				temp_koko++;
				if (tz0 < ty0 && tz0 < tx0) {
					tempk += ku;
					tz0 += tzu;
				}
				else if (ty0 < tx0) {
					tempj += ju;
					ty0 += tyu;
				}
				else {
					tempi += iu;
					tx0 += txu;
				}
				if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= d_Nx || tempj >= d_Ny || tempk >= d_Nz)
					break;
			}
			d_lor[idx] = temp_koko;
		}
	}
}



//__kernel void combine_lors(const __global ushort* d_lor, __global float* d_lor1, const uint m_loc) {
//
//	uint gid = get_global_id(0);
//	uint gid2 = gid + m_loc;
//
//	d_lor1[gid2] = d_lor[gid];
//}