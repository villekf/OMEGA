/**************************************************************************
* A matrix free orthogonal distance based ray tracer combined with all the 
* reconstruction functions available in OMEGA. This function calculates 
* Summ = sum(A,1) (sum of every row) and rhs = A*(y./(A'*x)), where A is 
* the system matrix, y the measurements and x the estimate/image.
*
* Used by implementation 2.
*
* This version goes through all the LORs and determines on-the-fly if they
* intersect with the voxel space.
*
* INPUTS:
* MethodList = The type of reconstruction algorithms used (e.g. 2 means
* COSEM)
* d_raw = if 1 then raw list-mode data is used otherwise sinogram
* data
* d_h = power factor for ACOSEM,
* d_Nx/y/z = image size in x/y/z- dimension,
* d_dz/x/y = distance between adjecent voxels in z/x/y-dimension,
* d_bz/x/y = distance from the pixel space to origin (z/x/y-dimension),
* d_bzb = part in parenthesis of equation (9) in [1] precalculated when
* k = Nz,
* d_maxxx/yy = maximum distance of the pixel space from origin in
* x/y-dimension,
* d_zmax = maximum value of d_zdet,
* d_NSlices = the number of image slices,
* d_x/y/z_det = detector x/y/z-coordinates,
* d_size_x = the number of detector elements,
* d_TotSinos = Total number of sinograms,
* d_attenuation_correction = if attenuation is included this is 1 otherwise
* 0,
* d_normalization = if normalization is included this is 1 otherwise 0,
* d_randoms = if randoms/scatter correction is included this is 1
* otherwise 0,
* d_atten = attenuation data (images),
* d_norm = normalization coefficients,
* d_epps = a small constant to prevent division by zero,
* d_N = d_Nx * d_Ny * d_Nz,
* d_pseudos = location of pseudo rings,
* pRows = number of pseudo rings,
* d_Nxy = d_Nx * d_Ny,
* d_det_per_ring = number of detectors per ring,
* n_rekos = number of reconstruction algorithms used,
* tube_width = the width of of the strip used for orthogonal distance based
* projector (2D),
* crystal_size_z = the width of of the tube used for orthogonal distance based
* projector (3D),
* dec = accuracy factor,
* x/y/z_center = Cartesian coordinates for the center of the voxels
* (x/y/z-axis),
* d_Summ = buffer for Summ,
* d_lor = number of pixels that each LOR traverses,
* d_xy/zindex = for sinogram format they determine the detector
* indices corresponding to each sinogram bin (unused with raw data),
* d_L = detector numbers for raw data (unused for sinogram format),
* d_epsilon_mramla = epsilon value for MRAMLA/MBSREM,
* d_Sino = Sinogram/raw data,
* d_sc_ra = Randoms and/or scatter data,
* d_OSEM = buffer for all estimates,
* d_rhs_OSEM = buffer for all RHS elements,
* no_norm = If 1, normalization constant is not computed,
* m_size = Total number of LORs for this subset,
* ax = Local buffer for forward projection data
*
* OUTPUTS:
* d_rhs_OSEM = rhs values for all algorithms/priors X,
* d_OSEM = estimates of all algorithm
* d_Summ = Normalization constant
*
* [1] Jacobs, F., Sundermann, E., De Sutter, B., Christiaens, M. Lemahieu,
* I. (1998). A Fast Algorithm to Calculate the Exact Radiological Path
* through a Pixel or Voxel Space. Journal of computing and information
* technology, 6 (1), 89-94.
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
#include "opencl_AF_functions.h"
#include "orthogonal_kernel_matrixfree.cl"
#define TYPE1 0

// Matrix free Orthogonal distance based ray tracer algorithm
__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(64, 1, 1)))
void orth_noprecomp(__constant uchar* MethodList, const uchar d_raw, const float d_h, const uint d_Nx, const uint d_Ny, const uint d_Nz,
	const float d_dz, const float d_dx, const float d_dy, const float d_bz, const float d_bx, const float d_by, const float d_bzb, const float d_maxxx,
	const float d_maxyy, const float d_zmax, const float d_NSlices, const __global float* d_x, const __global float* d_y, const __global float* d_zdet,
	const uint d_size_x, const uint d_TotSinos, const uint d_attenuation_correction, const uint d_normalization, const uint d_randoms,
	const __global float* d_atten, const __global float* d_norm, const float d_epps, const uint d_N, __constant uint* d_pseudos, const uint d_pRows,
	const uint d_Nxy, const uint d_det_per_ring, const uint n_rekos, const float tube_width_xy, const float crystal_size_z, const int dec,
	__constant float* x_center, __constant float* y_center, __constant float* z_center, __global float* d_Summ, const __global ushort* d_lor,
	const __global uint* d_xyindex, const __global ushort* d_zindex, const __global ushort* d_L, const float d_epsilon_mramla,
	const __global float* d_Sino, const __global float* d_sc_ra,
	const __global float* d_OSEM, __global float* d_rhs_OSEM, const uchar no_norm, const ulong m_size, __local float* ax) {
	// Get the current global index
	uint idx = get_global_id(0);
	if (idx >= m_size)
		return;
	float xs, xd, ys, yd, zs, zd;
	const float local_sino = (d_Sino[idx]);
	if (no_norm == 1u && local_sino == 0.f)
		return;
	uint lid = get_local_id(0);
	// Load the next detector index
	// raw list-mode data
	if (d_raw) {
		get_detector_coordinates_raw(d_x, d_y, d_zdet, d_L, d_det_per_ring, idx, d_pseudos, d_pRows, &xs, &xd, &ys, &yd, &zs, &zd);
	}
	// Sinogram data
	else {
		get_detector_coordinates(d_xyindex, d_zindex, d_size_x, idx, d_TotSinos, &xs, &xd, &ys, &yd, &zs, &zd, d_x, d_y, d_zdet);
	}
	// Calculate the x, y and z distances of the detector pair
	const float y_diff = (yd - ys);
	const float x_diff = (xd - xs);
	const float z_diff = (zd - zs);
	if ((y_diff == 0.f && x_diff == 0.f && z_diff == 0.f) || (y_diff == 0.f && x_diff == 0.f))
		return;
	uint Np = 0u;
	uint Np_n = 0u;
	uint xyz = 0u;
	float jelppi = 0.f, LL;
	for (uint kk = lid * n_rekos; kk < (n_rekos * (lid + 1u)); kk++)
		ax[kk] = 0.f;
	float kerroin, length_;
	if (crystal_size_z == 0.f) {
		kerroin = xd * ys - yd * xs;
		length_ = native_sqrt(y_diff * y_diff + x_diff * x_diff) * tube_width_xy;
	}
	else
		kerroin = e_norm(x_diff, y_diff, z_diff) * crystal_size_z;
	// If the measurement is on a same ring
	if (fabs(z_diff) < 1e-6f) {
		// Z-coordinate (ring)
		const uint z_loop = convert_uint((zs / d_zmax) * (d_NSlices - 1.f));
		// If the LOR is perpendicular in the y-direction (Siddon cannot be used)
		if (fabs(y_diff) < 1e-6f) {
			if (yd <= d_maxyy && yd >= d_by) {
				float temp = 0.f;
				if (crystal_size_z == 0.f) {
					orth_distance_denominator_perpendicular(x_diff, y_center, kerroin, length_, &temp, d_attenuation_correction, d_normalization, ax,
						d_by, yd, d_dy, d_Ny, d_Nx, z_loop, d_atten, d_norm, idx, local_sino, d_Ny, 1u, lid, n_rekos, d_N, d_OSEM);
					if (local_sino > 0.f) {
						nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, temp, d_randoms, d_sc_ra, idx, lid, n_rekos);
						orth_distance_rhs_perpendicular(x_diff, y_center, kerroin, length_, temp, MethodList, ax, d_by, yd, d_dy, d_Ny, d_Nx, z_loop, d_Ny, 1u,
							no_norm, lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
					}
					else {
						orth_distance_summ_perpendicular(x_diff, y_center, kerroin, length_, temp, d_by, yd, d_dy, d_Ny, d_Nx, z_loop, d_Ny, 1u, d_Summ);
					}
				}
				else {
					//const float temppi = xs;
					//xs = ys;
					//ys = temppi;
					orth_distance_denominator_perpendicular_3D(y_center, x_center[0], z_center, &temp, d_attenuation_correction, d_normalization, ax,
						d_by, yd, d_dy, d_Ny, d_Nx, z_loop, d_atten, d_norm, idx, local_sino, d_Ny, 1u, ys, xs, zs, y_diff, x_diff, z_diff, kerroin, d_Nxy, d_Nz,
						lid, n_rekos, d_N, d_OSEM);
					if (local_sino > 0.f) {
						nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, temp, d_randoms, d_sc_ra, idx, lid, n_rekos);
						orth_distance_rhs_perpendicular_3D(y_center, x_center[0], z_center, temp, MethodList, ax, d_by, yd, d_dy, d_Ny, d_Nx, z_loop, d_Ny, 1u,
							no_norm, ys, xs, zs, y_diff, x_diff, z_diff, kerroin, d_Nxy, d_Nz, lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
					}
					else {
						orth_distance_summ_perpendicular_3D(y_center, x_center[0], z_center, temp, d_by, yd, d_dy, d_Ny, d_Nx, z_loop, d_Ny, 1u,
							ys, xs, zs, y_diff, x_diff, z_diff, kerroin, d_Nxy, d_Nz, d_Summ);
					}
				}
			}
		}
		else if (fabs(x_diff) < 1e-6f) {
			if (xd <= d_maxxx && xd >= d_bx) {
				float temp = 0.f;
				if (crystal_size_z == 0.f) {
					orth_distance_denominator_perpendicular(-y_diff, x_center, kerroin, length_, &temp, d_attenuation_correction, d_normalization, ax,
						d_bx, xd, d_dx, d_Nx, d_Ny, z_loop, d_atten, d_norm, idx, local_sino, 1u, d_Nx, lid, n_rekos, d_N, d_OSEM);
					if (local_sino > 0.f) {
						nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, temp, d_randoms, d_sc_ra, idx, lid, n_rekos);
						orth_distance_rhs_perpendicular(-y_diff, x_center, kerroin, length_, temp, MethodList, ax, d_bx, xd, d_dx, d_Nx, d_Ny, z_loop, 1u, d_Nx,
							no_norm, lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
					}
					else {
						orth_distance_summ_perpendicular(-y_diff, x_center, kerroin, length_, temp, d_bx, xd, d_dx, d_Nx, d_Ny, z_loop, 1u, d_Nx, d_Summ);
					}
				}
				else {
					orth_distance_denominator_perpendicular_3D(x_center, y_center[0], z_center, &temp, d_attenuation_correction, d_normalization, ax,
						d_bx, xd, d_dx, d_Nx, d_Ny, z_loop, d_atten, d_norm, idx, local_sino, 1u, d_Nx, xs, ys, zs, x_diff, y_diff, z_diff, kerroin, d_Nxy, d_Nz,
						lid, n_rekos, d_N, d_OSEM);
					if (local_sino > 0.f) {
						nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, temp, d_randoms, d_sc_ra, idx, lid, n_rekos);
						orth_distance_rhs_perpendicular_3D(x_center, y_center[0], z_center, temp, MethodList, ax, d_bx, xd, d_dx, d_Nx, d_Ny, z_loop, 1u, d_Nx,
							no_norm, xs, ys, zs, x_diff, y_diff, z_diff, kerroin, d_Nxy, d_Nz, lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
					}
					else {
						orth_distance_summ_perpendicular_3D(x_center, y_center[0], z_center, temp, d_bx, xd, d_dx, d_Nx, d_Ny, z_loop, 1u, d_Nx,
							xs, ys, zs, x_diff, y_diff, z_diff, kerroin, d_Nxy, d_Nz, d_Summ);
					}
				}
			}
		}
		else {
			int tempi = 0, tempj = 0, iu = 0, ju = 0;
			float txu = 0.f, tyu = 0.f, tc = 0.f, tx0 = 0.f, ty0 = 0.f;
			const bool skip = siddon_pre_loop_2D(d_bx, d_by, x_diff, y_diff, d_maxxx, d_maxyy, d_dx, d_dy, d_Nx, d_Ny, &tempi, &tempj, &txu, &tyu, &Np, TYPE1,
				ys, xs, yd, xd, &tc, &iu, &ju, &tx0, &ty0);
			if (skip)
				return;
			float temp = 0.f;
			float tx0_a = tx0, ty0_a = ty0;
			int tempi_a = tempi, tempj_a = tempj;
			int tempk = z_loop;
			uint temp_ijk;
			if (crystal_size_z == 0.f)
				temp_ijk = z_loop * d_Nxy + convert_uint_sat(tempj) * d_Nx;
			else
				temp_ijk = convert_uint_sat(tempj) * d_Nx;
			if (d_attenuation_correction == 1u)
				LL = native_sqrt(x_diff * x_diff + y_diff * y_diff);
			for (uint ii = 0u; ii < Np; ii++) {
				if (tx0 < ty0) {
					if (d_attenuation_correction == 1u)
						compute_attenuation(&tc, &jelppi, LL, tx0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
					if (ii == Np - 1u) {
						if (crystal_size_z == 0.f) {
							orth_distance_denominator(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, lid, n_rekos, d_N, d_OSEM);
						}
						else {
							orth_distance_denominator_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, lid, n_rekos, d_N, d_OSEM);
						}
					}
					else {
						tempi += iu;
						tx0 += txu;
					}
					xyz = 1u;
				}
				else {
					if (d_attenuation_correction == 1u)
						compute_attenuation(&tc, &jelppi, LL, ty0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
					if (crystal_size_z == 0.f) {
						orth_distance_denominator(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, lid, n_rekos, d_N, d_OSEM);
					}
					else {
						orth_distance_denominator_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, lid, n_rekos, d_N, d_OSEM);
					}
					tempj += ju;
					if (ju > 0)
						temp_ijk += d_Nx;
					else
						temp_ijk -= d_Nx;
					ty0 += tyu;
					xyz = 2u;
				}
				Np_n++;
				if (tempj < 0 || tempi < 0 || tempi >= d_Nx || tempj >= d_Ny) {
					if (xyz == 1u) {
						if (crystal_size_z == 0.f) {
							orth_distance_denominator(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk,
								local_sino, ax, d_Ny, 1u, lid, n_rekos, d_N, d_OSEM);
						}
						else {
							orth_distance_denominator_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, lid, n_rekos, d_N, d_OSEM);
						}
					}
					break;
				}
			}
			temp = 1.f / temp;
			if (d_attenuation_correction == 1u)
				temp *= native_exp(jelppi);
			if (d_normalization == 1u)
				temp *= d_norm[idx];
			tx0 = tx0_a, ty0 = ty0_a;
			tempi = tempi_a, tempj = tempj_a;
			if (crystal_size_z == 0.f)
				temp_ijk = z_loop * d_Nxy + convert_uint_sat(tempj) * d_Nx;
			else
				temp_ijk = convert_uint_sat(tempj) * d_Nx;
			if (local_sino > 0.f) {
				nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, temp, d_randoms, d_sc_ra, idx, lid, n_rekos);
				for (uint ii = 0u; ii < Np_n; ii++) {
					if (tx0 < ty0) {
						if (ii == Np_n - 1u) {
							if (crystal_size_z == 0.f) {
								orth_distance_rhs(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, temp_ijk, MethodList, ax, 
									d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
							}
							else {
								orth_distance_rhs_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, temp, temp_ijk, 
									MethodList, ax, d_Ny, 1u, d_Nxy, xs, ys, zs, tempk, dec, no_norm, lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
							}
						}
						else {
							tempi += iu;
							tx0 += txu;
						}
					}
					else {
						if (crystal_size_z == 0.f) {
							orth_distance_rhs(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, temp_ijk, MethodList, ax, 
									d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
						}
						else {
							orth_distance_rhs_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, temp, temp_ijk, 
									MethodList, ax, d_Ny, 1u, d_Nxy, xs, ys, zs, tempk, dec, no_norm, lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
						}
						tempj += ju;
						if (ju > 0)
							temp_ijk += d_Nx;
						else
							temp_ijk -= d_Nx;
						ty0 += tyu;
					}
				}
			}
			else {
				for (uint ii = 0u; ii < Np_n; ii++) {
					if (tx0 < ty0) {
						if (ii == Np_n - 1u) {
							if (crystal_size_z == 0.f) {
								orth_distance_summ(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, temp_ijk, d_Summ, d_Ny, 1u);
							}
							else {
								orth_distance_summ_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, temp, temp_ijk,
									d_Ny, 1u, d_Nxy, xs, ys, zs, tempk, dec, d_Summ);
							}
						}
						else {
							tempi += iu;
							tx0 += txu;
						}
					}
					else {
						if (crystal_size_z == 0.f) {
							orth_distance_summ(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, temp_ijk, d_Summ, d_Ny, 1u);
						}
						else {
							orth_distance_summ_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, temp, temp_ijk,
								d_Ny, 1u, d_Nxy, xs, ys, zs, tempk, dec, d_Summ);
						}
						tempj += ju;
						if (ju > 0)
							temp_ijk += d_Nx;
						else
							temp_ijk -= d_Nx;
						ty0 += tyu;
					}
				}
			}
		}
	}
	else {
		if (fabs(y_diff) < 1e-6f) {
			if (yd <= d_maxyy && yd >= d_by) {
				int tempi = 0, tempk = 0, tempj = 0, iu = 0, ku = 0;
				float txu = 0.f, tzu = 0.f, tc = 0.f, tx0 = 0.f, tz0 = 0.f;
				const bool skip = siddon_pre_loop_2D(d_bx, d_bz, x_diff, z_diff, d_maxxx, d_bzb, d_dx, d_dz, d_Nx, d_Nz, &tempi, &tempk, &txu, &tzu, &Np, TYPE1,
					zs, xs, zd, xd, &tc, &iu, &ku, &tx0, &tz0);
				if (skip)
					return;
				tempj = perpendicular_start(d_by, yd, d_dy, d_Ny);
				float temp = 0.f;
				float tx0_a = tx0, tz0_a = tz0;
				int tempi_a = tempi, tempk_a = tempk;
				uint temp_ijk;
				if (crystal_size_z == 0.f)
					temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempi);
				else {
					temp_ijk = convert_uint_sat(tempi);
					const float temp_x = xs;
					xs = ys;
					ys = temp_x;
				}
				if (d_attenuation_correction == 1u)
					LL = native_sqrt(x_diff * x_diff + z_diff * z_diff);
				for (uint ii = 0u; ii < Np; ii++) {
					if (tx0 < tz0) {
						if (d_attenuation_correction == 1u)
							compute_attenuation(&tc, &jelppi, LL, tx0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
						if (crystal_size_z == 0.f) {
							orth_distance_denominator(tempj, d_Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Nx, d_Nx, lid, n_rekos, d_N, d_OSEM);
						}
						else {
							orth_distance_denominator_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Nx, d_Nx, tempk, d_Nxy, xs, ys, zs, dec, lid, n_rekos, d_N, d_OSEM);
						}
						if (iu > 0)
							temp_ijk++;
						else
							temp_ijk--;
						tempi += iu;
						tx0 += txu;
						xyz = 1u;
					}
					else {
						if (d_attenuation_correction == 1u)
							compute_attenuation(&tc, &jelppi, LL, tz0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
						if (crystal_size_z == 0.f) {
							orth_distance_denominator(tempj, d_Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Nx, d_Nx, lid, n_rekos, d_N, d_OSEM);
							if (ku > 0)
								temp_ijk += d_Nxy;
							else
								temp_ijk -= d_Nxy;
						}
						else if (ii == Np - 1u) {
							orth_distance_denominator_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Nx, d_Nx, tempk, d_Nxy, xs, ys, zs, dec, lid, n_rekos, d_N, d_OSEM);
						}
						tempk += ku;
						tz0 += tzu;
						xyz = 3u;
					}
					Np_n++;
					if (tempk < 0 || tempi < 0 || tempi >= d_Nx || tempk >= d_Nz) {
						if (crystal_size_z != 0.f && xyz == 3u) {
							orth_distance_denominator_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, &temp, temp_ijk,
								local_sino, ax, d_Nx, d_Nx, tempk, d_Nxy, xs, ys, zs, dec, lid, n_rekos, d_N, d_OSEM);
						}
						break;
					}
				}
				temp = 1.f / temp;
				if (d_attenuation_correction == 1u)
					temp *= native_exp(jelppi);
				if (d_normalization == 1u)
					temp *= d_norm[idx];
				tx0 = tx0_a, tz0 = tz0_a;
				tempi = tempi_a, tempk = tempk_a;
				if (crystal_size_z == 0.f)
					temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempi);
				else
					temp_ijk = convert_uint_sat(tempi);
				if (local_sino > 0.f) {
					nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, temp, d_randoms, d_sc_ra, idx, lid, n_rekos);
					for (uint ii = 0u; ii < Np_n; ii++) {
						if (tx0 < tz0) {
							if (crystal_size_z == 0.f) {
								orth_distance_rhs(tempj, d_Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, temp, temp_ijk, MethodList, ax, d_Nx, d_Nx, no_norm,
									lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
							}
							else {
								orth_distance_rhs_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, temp, temp_ijk, 
									MethodList, ax, d_Nx, d_Nx, d_Nxy, xs, ys, zs, tempk, dec, no_norm, lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
							}
							if (iu > 0)
								temp_ijk++;
							else
								temp_ijk--;
							tempi += iu;
							tx0 += txu;
						}
						else {
							if (crystal_size_z == 0.f) {
								orth_distance_rhs(tempj, d_Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, temp, temp_ijk, MethodList, ax, d_Nx, d_Nx, no_norm,
									lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
								if (ku > 0)
									temp_ijk += d_Nxy;
								else
									temp_ijk -= d_Nxy;
							}
							else if (ii == Np_n - 1u) {
								orth_distance_rhs_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, temp, temp_ijk, 
									MethodList, ax, d_Nx, d_Nx, d_Nxy, xs, ys, zs, tempk, dec, no_norm, lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
							}
							tz0 += tzu;
							tempk += ku;
						}
					}
				}
				else {
					for (uint ii = 0u; ii < Np_n; ii++) {
						if (tx0 < tz0) {
							if (crystal_size_z == 0.f) {
								orth_distance_summ(tempj, d_Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, temp, temp_ijk, d_Summ, d_Nx, d_Nx);
							}
							else {
								orth_distance_summ_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, temp, temp_ijk,
									d_Nx, d_Nx, d_Nxy, xs, ys, zs, tempk, dec, d_Summ);
							}
							if (iu > 0)
								temp_ijk++;
							else
								temp_ijk--;
							tempi += iu;
						}
						else {
							if (crystal_size_z == 0.f) {
								orth_distance_summ(tempj, d_Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, temp, temp_ijk, d_Summ, d_Nx, d_Nx);
								if (ku > 0)
									temp_ijk += d_Nxy;
								else
									temp_ijk -= d_Nxy;
							}
							else if (ii == Np_n - 1u) {
								orth_distance_summ_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, temp, temp_ijk,
									d_Nx, d_Nx, d_Nxy, xs, ys, zs, tempk, dec, d_Summ);
							}
							tempk += ku;
							tz0 += tzu;
						}
					}
				}
			}
		}
		else if (fabs(x_diff) < 1e-6f) {
			if (xd <= d_maxxx && xd >= d_bx) {
				int tempi = 0, tempk = 0, tempj = 0, ju = 0, ku = 0;
				float tyu = 0.f, tzu = 0.f, tc = 0.f, ty0 = 0.f, tz0 = 0.f;
				const bool skip = siddon_pre_loop_2D(d_by, d_bz, y_diff, z_diff, d_maxyy, d_bzb, d_dy, d_dz, d_Ny, d_Nz, &tempj, &tempk, &tyu, &tzu, &Np, TYPE1,
					zs, ys, zd, yd, &tc, &ju, &ku, &ty0, &tz0);
				if (skip)
					return;
				float temp = 0.f;
				tempi = perpendicular_start(d_bx, xd, d_dx, d_Nx);
				const float tz0_a = tz0, ty0_a = ty0;
				const int tempj_a = tempj, tempk_a = tempk;
				uint temp_ijk;
				if (crystal_size_z == 0.f)
					temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempj) * d_Nx;
				else
					temp_ijk = convert_uint_sat(tempj) * d_Nx;
				if (d_attenuation_correction == 1u)
					LL = native_sqrt(y_diff * y_diff + z_diff * z_diff);
				for (uint ii = 0u; ii < Np; ii++) {
					if (tz0 < ty0) {
						if (d_attenuation_correction == 1u)
							compute_attenuation(&tc, &jelppi, LL, tz0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
						if (crystal_size_z == 0.f) {
							orth_distance_denominator(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, lid, n_rekos, d_N, d_OSEM);
							if (ku > 0)
								temp_ijk += d_Nxy;
							else
								temp_ijk -= d_Nxy;
						}
						else if (ii == Np - 1u) {
							orth_distance_denominator_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, lid, n_rekos, d_N, d_OSEM);
						}
						tempk += ku;
						tz0 += tzu;
						xyz = 3u;
					}
					else {
						if (d_attenuation_correction == 1u)
							compute_attenuation(&tc, &jelppi, LL, ty0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
						if (crystal_size_z == 0.f) {
							orth_distance_denominator(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, lid, n_rekos, d_N, d_OSEM);
						}
						else {
							orth_distance_denominator_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, lid, n_rekos, d_N, d_OSEM);
						}
						tempj += (ju);
						if (ju > 0)
							temp_ijk += d_Nx;
						else
							temp_ijk -= d_Nx;
						ty0 += tyu;
						xyz = 2u;
					}
					Np_n++;
					if (tempj < 0 || tempk < 0 || tempk >= d_Nz || tempj >= d_Ny) {
						if (xyz == 3u && crystal_size_z != 0.f) {
							orth_distance_denominator_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, lid, n_rekos, d_N, d_OSEM);
						}
						break;
					}
				}
				temp = 1.f / temp;
				if (d_attenuation_correction == 1u)
					temp *= native_exp(jelppi);
				if (d_normalization == 1u)
					temp *= d_norm[idx];
				tz0 = tz0_a, ty0 = ty0_a;
				tempj = tempj_a, tempk = tempk_a;
				if (crystal_size_z == 0.f)
					temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempj) * d_Nx;
				else
					temp_ijk = convert_uint_sat(tempj) * d_Nx;
				if (local_sino > 0.f) {
					nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, temp, d_randoms, d_sc_ra, idx, lid, n_rekos);
					for (uint ii = 0u; ii < Np_n; ii++) {
						if (tz0 < ty0) {
							if (crystal_size_z == 0.f) {
								orth_distance_rhs(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, temp_ijk, MethodList, ax, 
									d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
								if (ku > 0)
									temp_ijk += d_Nxy;
								else
									temp_ijk -= d_Nxy;
							}
							else if (ii == Np_n - 1u) {
								orth_distance_rhs_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, temp, temp_ijk, 
									MethodList, ax, d_Ny, 1u, d_Nxy, xs, ys, zs, tempk, dec, no_norm, lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
							}
							tempk += ku;
							tz0 += tzu;
						}
						else {
							if (crystal_size_z == 0.f) {
								orth_distance_rhs(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, temp_ijk, MethodList, ax, 
									d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
							}
							else {
								orth_distance_rhs_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, temp, temp_ijk, 
									MethodList, ax, d_Ny, 1u, d_Nxy, xs, ys, zs, tempk, dec, no_norm, lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);									
							}
							if (ju > 0)
								temp_ijk += d_Nx;
							else
								temp_ijk -= d_Nx;
							tempj += (ju);
							ty0 += tyu;
						}
					}
				}
				else {
					for (uint ii = 0u; ii < Np_n; ii++) {
						if (tz0 < ty0) {
							if (crystal_size_z == 0.f) {
								orth_distance_summ(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, temp_ijk, d_Summ, d_Ny, 1u);
								if (ku > 0)
									temp_ijk += d_Nxy;
								else
									temp_ijk -= d_Nxy;
							}
							else if (ii == Np_n - 1u) {
								orth_distance_summ_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, temp, temp_ijk,
									d_Ny, 1u, d_Nxy, xs, ys, zs, tempk, dec, d_Summ);
							}
							tempk += ku;
							tz0 += tzu;
						}
						else {
							if (crystal_size_z == 0.f) {
								orth_distance_summ(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, temp_ijk, d_Summ, d_Ny, 1u);
							}
							else {
								orth_distance_summ_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, temp, temp_ijk,
									d_Ny, 1u, d_Nxy, xs, ys, zs, tempk, dec, d_Summ);
							}
							if (ju > 0)
								temp_ijk += d_Nx;
							else
								temp_ijk -= d_Nx;
							tempj += (ju);
							ty0 += tyu;
						}
					}
				}
			}
		}
		else {
			int tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;
			float txu = 0.f, tyu = 0.f, tzu = 0.f, tc = 0.f, tx0 = 0.f, ty0 = 0.f, tz0 = 0.f;
			const bool skip = siddon_pre_loop_3D(d_bx, d_by, d_bz, x_diff, y_diff, z_diff, d_maxxx, d_maxyy, d_bzb, d_dx, d_dy, d_dz, d_Nx, d_Ny, d_Nz, &tempi,
				&tempj, &tempk, &tyu, &txu, &tzu, &Np, TYPE1, ys, xs, yd, xd, zs, zd, &tc, &iu, &ju, &ku, &tx0, &ty0, &tz0);
			if (skip)
				return;
			float temp = 0.f;
			const float tx0_a = tx0, ty0_a = ty0, tz0_a = tz0;
			const int tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
			uint temp_ijk;
			if (crystal_size_z == 0.f)
				temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempj) * d_Nx;
			else
				temp_ijk = convert_uint_sat(tempj) * d_Nx;
			if (d_attenuation_correction == 1u)
				LL = native_sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
			for (uint ii = 0u; ii < Np; ii++) {
				if (tz0 < ty0 && tz0 < tx0) {
					if (d_attenuation_correction == 1u)
						compute_attenuation(&tc, &jelppi, LL, tz0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
					if (crystal_size_z == 0.f) {
						orth_distance_denominator(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, lid, n_rekos, d_N, d_OSEM);
						if (ku > 0)
							temp_ijk += d_Nxy;
						else
							temp_ijk -= d_Nxy;
					}
					else if (ii == Np - 1u) {
						orth_distance_denominator_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, lid, n_rekos, d_N, d_OSEM);
					}
					tempk += ku;
					tz0 += tzu;
					xyz = 3u;
				}
				else if (ty0 < tx0) {
					if (d_attenuation_correction == 1u)
						compute_attenuation(&tc, &jelppi, LL, ty0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
					if (crystal_size_z == 0.f) {
						orth_distance_denominator(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, lid, n_rekos, d_N, d_OSEM);
					}
					else {
						orth_distance_denominator_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, lid, n_rekos, d_N, d_OSEM);
					}
					if (ju > 0)
						temp_ijk += d_Nx;
					else
						temp_ijk -= d_Nx;
					tempj += (ju);
					ty0 += tyu;
					xyz = 2u;
				}
				else {
					if (d_attenuation_correction == 1u)
						compute_attenuation(&tc, &jelppi, LL, tx0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
					if (ii == Np - 1) {
						if (crystal_size_z == 0.f) {
							orth_distance_denominator(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, lid, n_rekos, d_N, d_OSEM);
						}
						else {
							orth_distance_denominator_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, lid, n_rekos, d_N, d_OSEM);
						}
					}
					else {
						tempi += iu;
						tx0 += txu;
					}
					xyz = 1u;
				}
				Np_n++;
				if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= d_Nx || tempj >= d_Ny || tempk >= d_Nz) {
					if (xyz == 1u || (xyz == 3u && crystal_size_z != 0.f)) {
						if (crystal_size_z == 0.f) {
							orth_distance_denominator(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk,
								local_sino, ax, d_Ny, 1u, lid, n_rekos, d_N, d_OSEM);
						}
						else {
							orth_distance_denominator_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, lid, n_rekos, d_N, d_OSEM);
						}
					}
					break;
				}
			}
			temp = 1.f / temp;
			if (d_attenuation_correction == 1u)
				temp *= native_exp(jelppi);
			if (d_normalization == 1u)
				temp *= d_norm[idx];
			tx0 = tx0_a, ty0 = ty0_a, tz0 = tz0_a;
			tempi = tempi_a, tempj = tempj_a, tempk = tempk_a;
			if (crystal_size_z == 0.f)
				temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempj) * d_Nx;
			else
				temp_ijk = convert_uint_sat(tempj) * d_Nx;
			if (local_sino > 0.f) {
				nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, temp, d_randoms, d_sc_ra, idx, lid, n_rekos);
				for (uint ii = 0u; ii < Np_n; ii++) {
					if (tz0 < ty0 && tz0 < tx0) {
						if (crystal_size_z == 0.f) {
							orth_distance_rhs(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, temp_ijk, MethodList, ax, 
									d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
							if (ku > 0)
								temp_ijk += d_Nxy;
							else
								temp_ijk -= d_Nxy;
						}
						else if (ii == Np_n - 1u) {
							orth_distance_rhs_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, temp, temp_ijk, 
									MethodList, ax, d_Ny, 1u, d_Nxy, xs, ys, zs, tempk, dec, no_norm, lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
						}
						tempk += ku;
						tz0 += tzu;
					}
					else if (ty0 < tx0) {
						if (crystal_size_z == 0.f) {
							orth_distance_rhs(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, temp_ijk, MethodList, ax, 
									d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
						}
						else {
							orth_distance_rhs_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, temp, temp_ijk, 
									MethodList, ax, d_Ny, 1u, d_Nxy, xs, ys, zs, tempk, dec, no_norm, lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
						}
						if (ju > 0)
							temp_ijk += d_Nx;
						else
							temp_ijk -= d_Nx;
						tempj += (ju);
						ty0 += tyu;
					}
					else {
						if (ii == Np_n - 1) {
							if (crystal_size_z == 0.f) {
								orth_distance_rhs(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, temp_ijk, MethodList, ax, 
									d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
							}
							else {
								orth_distance_rhs_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, temp, temp_ijk, 
									MethodList, ax, d_Ny, 1u, d_Nxy, xs, ys, zs, tempk, dec, no_norm, lid, n_rekos, d_N, d_rhs_OSEM, d_h, d_OSEM, d_Summ);
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
				for (uint ii = 0u; ii < Np_n; ii++) {
					if (tz0 < ty0 && tz0 < tx0) {
						if (crystal_size_z == 0.f) {
							orth_distance_summ(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, temp_ijk, d_Summ, d_Ny, 1u);
							if (ku > 0)
								temp_ijk += d_Nxy;
							else
								temp_ijk -= d_Nxy;
						}
						else if (ii == Np_n - 1u) {
							orth_distance_summ_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, temp, temp_ijk,
								d_Ny, 1u, d_Nxy, xs, ys, zs, tempk, dec, d_Summ);
						}
						tempk += ku;
						tz0 += tzu;
					}
					else if (ty0 < tx0) {
						if (crystal_size_z == 0.f) {
							orth_distance_summ(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, temp_ijk, d_Summ, d_Ny, 1u);
						}
						else {
							orth_distance_summ_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, temp, temp_ijk,
								d_Ny, 1u, d_Nxy, xs, ys, zs, tempk, dec, d_Summ);
						}
						if (ju > 0)
							temp_ijk += d_Nx;
						else
							temp_ijk -= d_Nx;
						tempj += (ju);
						ty0 += tyu;
					}
					else {
						if (ii == Np_n - 1) {
							if (crystal_size_z == 0.f) {
								orth_distance_summ(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, temp, temp_ijk, d_Summ, d_Ny, 1u);
							}
							else {
								orth_distance_summ_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, temp, temp_ijk,
									d_Ny, 1u, d_Nxy, xs, ys, zs, tempk, dec, d_Summ);
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






__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(64, 1, 1)))
void orth_ml_noprecomp(const uchar d_raw, const uint d_Nx, const uint d_Ny, const uint d_Nz,
	const float d_dz, const float d_dx, const float d_dy, const float d_bz, const float d_bx, const float d_by, const float d_bzb, const float d_maxxx,
	const float d_maxyy, const float d_zmax, const float d_NSlices, const __global float* d_x, const __global float* d_y, const __global float* d_zdet,
	const uint d_size_x, const uint d_TotSinos, const uint d_attenuation_correction, const uint d_normalization, const uint d_randoms,
	const __global float* d_atten, const __global float* d_norm, const float d_epps, const uint d_N, __constant uint* d_pseudos, const uint d_pRows,
	const uint d_Nxy, const uint d_det_per_ring, const uint n_rekos, const float tube_width_xy, const float crystal_size_z, const int dec, __constant float* x_center,
	__constant float* y_center, __constant float* z_center, __global float* d_Summ, const __global ushort* d_lor, const __global uint* d_xyindex,
	const __global ushort* d_zindex, const __global ushort* d_L, const __global float* d_Sino, const __global float* d_sc_ra,
	const __global float* d_MLEM, __global float* d_rhs_MLEM, const uchar no_norm, const ulong m_size, __local float* ax) {
	uint idx = get_global_id(0);
	if (idx >= m_size)
		return;
	const float local_sino = (d_Sino[idx]);
	if (no_norm == 1u && local_sino == 0.f)
		return;
	float xs, xd, ys, yd, zs, zd;
	uint lid = get_local_id(0);
	if (d_raw) {
		get_detector_coordinates_raw(d_x, d_y, d_zdet, d_L, d_det_per_ring, idx, d_pseudos, d_pRows, &xs, &xd, &ys, &yd, &zs, &zd);
	}
	// Sinogram data
	else {
		get_detector_coordinates(d_xyindex, d_zindex, d_size_x, idx, d_TotSinos, &xs, &xd, &ys, &yd, &zs, &zd, d_x, d_y, d_zdet);
	}
	const float y_diff = (yd - ys);
	const float x_diff = (xd - xs);
	const float z_diff = (zd - zs);
	if ((y_diff == 0.f && x_diff == 0.f && z_diff == 0.f) || (y_diff == 0.f && x_diff == 0.f))
		return;
	uint Np = 0u;
	uint Np_n = 0u;
	uint xyz = 0u;
	float jelppi = 0.f, LL;
	for (uint kk = lid * n_rekos; kk < (n_rekos * (lid + 1u)); kk++)
		ax[kk] = 0.f;
	float kerroin, length_;
	if (crystal_size_z == 0.f) {
		kerroin = xd * ys - yd * xs;
		length_ = native_sqrt(y_diff * y_diff + x_diff * x_diff) * tube_width_xy;
	}
	else
		kerroin = e_norm(x_diff, y_diff, z_diff) * crystal_size_z;
	float temp = 0.f;
	if (fabs(z_diff) < 1e-8f) {
		const uint z_loop = convert_uint((zs / d_zmax) * (d_NSlices - 1.f));
		// If the LOR is perpendicular in the y-direction (Siddon cannot be used)
		if (fabs(y_diff) < 1e-6f) {
			if (yd <= d_maxyy && yd >= d_by) {
				if (crystal_size_z == 0.f) {
					orth_distance_perpendicular_mlem(-x_diff, y_center, kerroin, length_, &temp, d_attenuation_correction, d_normalization, ax,
						d_by, yd, d_dy, d_Ny, d_Nx, z_loop, d_atten, d_norm, idx, local_sino, d_Ny, 1u, no_norm, lid, n_rekos, d_N,
						d_MLEM, d_rhs_MLEM, d_Summ, false, false);
					if (local_sino > 0.f) {
						nominator_mlem(ax, local_sino, d_epps, temp, d_randoms, d_sc_ra, idx, lid, n_rekos);
						orth_distance_perpendicular_mlem(-x_diff, y_center, kerroin, length_, &temp, d_attenuation_correction, d_normalization, ax,
							d_by, yd, d_dy, d_Ny, d_Nx, z_loop, d_atten, d_norm, idx, local_sino, d_Ny, 1u, no_norm, lid, n_rekos, d_N,
							d_MLEM, d_rhs_MLEM, d_Summ, true, false);
					}
					else {
						orth_distance_perpendicular_mlem(-x_diff, y_center, kerroin, length_, &temp, d_attenuation_correction, d_normalization, ax,
							d_by, yd, d_dy, d_Ny, d_Nx, z_loop, d_atten, d_norm, idx, local_sino, d_Ny, 1u, no_norm, lid, n_rekos, d_N,
							d_MLEM, d_rhs_MLEM, d_Summ, false, true);
					}
				}
				else {
					//const float temppi = xs;
					//xs = ys;
					//ys = temppi;
					orth_distance_perpendicular_mlem_3D(y_center, x_center[0], z_center, &temp, d_attenuation_correction, d_normalization, ax,
						d_by, yd, d_dy, d_Ny, d_Nx, z_loop, d_atten, d_norm, idx, local_sino, d_Ny, 1u, ys, xs, zs, y_diff, x_diff, z_diff, kerroin, d_Nxy, d_Nz,
						no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
					if (local_sino > 0.f) {
						nominator_mlem(ax, local_sino, d_epps, temp, d_randoms, d_sc_ra, idx, lid, n_rekos);
						orth_distance_perpendicular_mlem_3D(y_center, x_center[0], z_center, &temp, d_attenuation_correction, d_normalization, ax,
							d_by, yd, d_dy, d_Ny, d_Nx, z_loop, d_atten, d_norm, idx, local_sino, d_Ny, 1u, ys, xs, zs, y_diff, x_diff, z_diff, kerroin, d_Nxy, d_Nz,
							no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, true, false);
					}
					else {
						orth_distance_perpendicular_mlem_3D(y_center, x_center[0], z_center, &temp, d_attenuation_correction, d_normalization, ax,
							d_by, yd, d_dy, d_Ny, d_Nx, z_loop, d_atten, d_norm, idx, local_sino, d_Ny, 1u, ys, xs, zs, y_diff, x_diff, z_diff, kerroin, d_Nxy, d_Nz,
							no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, true);
					}
				}
			}
		}
		else if (fabs(x_diff) < 1e-6f) {
			if (xd <= d_maxxx && xd >= d_bx) {
				if (crystal_size_z == 0.f) {
					orth_distance_perpendicular_mlem(y_diff, x_center, kerroin, length_, &temp, d_attenuation_correction, d_normalization, ax,
						d_bx, xd, d_dx, d_Nx, d_Ny, z_loop, d_atten, d_norm, idx, local_sino, 1u, d_Nx, no_norm, lid, n_rekos, d_N,
						d_MLEM, d_rhs_MLEM, d_Summ, false, false);
					if (local_sino > 0.f) {
						nominator_mlem(ax, local_sino, d_epps, temp, d_randoms, d_sc_ra, idx, lid, n_rekos);
						orth_distance_perpendicular_mlem(y_diff, x_center, kerroin, length_, &temp, d_attenuation_correction, d_normalization, ax,
							d_bx, xd, d_dx, d_Nx, d_Ny, z_loop, d_atten, d_norm, idx, local_sino, 1u, d_Nx, no_norm, lid, n_rekos, d_N,
							d_MLEM, d_rhs_MLEM, d_Summ, true, false);
					}
					else {
						orth_distance_perpendicular_mlem(y_diff, x_center, kerroin, length_, &temp, d_attenuation_correction, d_normalization, ax,
							d_bx, xd, d_dx, d_Nx, d_Ny, z_loop, d_atten, d_norm, idx, local_sino, 1u, d_Nx, no_norm, lid, n_rekos, d_N,
							d_MLEM, d_rhs_MLEM, d_Summ, false, true);
					}
				}
				else {
					orth_distance_perpendicular_mlem_3D(x_center, y_center[0], z_center, &temp, d_attenuation_correction, d_normalization, ax,
						d_bx, xd, d_dx, d_Nx, d_Ny, z_loop, d_atten, d_norm, idx, local_sino, 1u, d_Nx, xs, ys, zs, x_diff, y_diff, z_diff, kerroin, d_Nxy, d_Nz,
						no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
					if (local_sino > 0.f) {
						nominator_mlem(ax, local_sino, d_epps, temp, d_randoms, d_sc_ra, idx, lid, n_rekos);
						orth_distance_perpendicular_mlem_3D(x_center, y_center[0], z_center, &temp, d_attenuation_correction, d_normalization, ax,
							d_bx, xd, d_dx, d_Nx, d_Ny, z_loop, d_atten, d_norm, idx, local_sino, 1u, d_Nx, xs, ys, zs, x_diff, y_diff, z_diff, kerroin, d_Nxy, d_Nz,
							no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, true, false);
					}
					else {
						orth_distance_perpendicular_mlem_3D(x_center, y_center[0], z_center, &temp, d_attenuation_correction, d_normalization, ax,
							d_bx, xd, d_dx, d_Nx, d_Ny, z_loop, d_atten, d_norm, idx, local_sino, 1u, d_Nx, xs, ys, zs, x_diff, y_diff, z_diff, kerroin, d_Nxy, d_Nz,
							no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, true);
					}
				}
			}
		}
		else {
			int tempi = 0, tempj = 0, iu = 0, ju = 0;
			float txu = 0.f, tyu = 0.f, tc = 0.f, tx0 = 0.f, ty0 = 0.f;
			const bool skip = siddon_pre_loop_2D(d_bx, d_by, x_diff, y_diff, d_maxxx, d_maxyy, d_dx, d_dy, d_Nx, d_Ny, &tempi, &tempj, &txu, &tyu, &Np, TYPE1,
				ys, xs, yd, xd, &tc, &iu, &ju, &tx0, &ty0);
			if (skip)
				return;
			float tx0_a = tx0, ty0_a = ty0;
			int tempi_a = tempi, tempj_a = tempj;
			int tempk = z_loop;
			uint temp_ijk;
			if (crystal_size_z == 0.f)
				temp_ijk = z_loop * d_Nxy + convert_uint_sat(tempj) * d_Nx;
			else
				temp_ijk = convert_uint_sat(tempj) * d_Nx;
			if (d_attenuation_correction == 1u)
				LL = native_sqrt(x_diff * x_diff + y_diff * y_diff);
			for (uint ii = 0u; ii < Np; ii++) {
				if (tx0 < ty0) {
					if (d_attenuation_correction == 1u)
						compute_attenuation(&tc, &jelppi, LL, tx0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
					if (ii == Np - 1u) {
						if (crystal_size_z == 0.f) {
							orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
						}
						else {
							orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
						}
					}
					else {
						tempi += iu;
						tx0 += txu;
						xyz = 1u;
					}
				}
				else {
					if (d_attenuation_correction == 1u)
						compute_attenuation(&tc, &jelppi, LL, ty0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
					if (crystal_size_z == 0.f) {
						orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
					}
					else {
						orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
					}
					if (ju > 0)
						temp_ijk += d_Nx;
					else
						temp_ijk -= d_Nx;
					tempj += ju;
					ty0 += tyu;
					xyz = 2u;
				}
				Np_n++;
				if (tempj < 0 || tempi < 0 || tempi >= d_Nx || tempj >= d_Ny) {
					if (xyz == 1u) {
						if (crystal_size_z == 0.f) {
							orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk,
								local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
						}
						else {
							orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
						}
					}
					break;
				}
			}
			temp = 1.f / temp;
			if (d_attenuation_correction == 1u)
				temp *= native_exp(jelppi);
			if (d_normalization == 1u)
				temp *= d_norm[idx];
			tx0 = tx0_a, ty0 = ty0_a;
			tempi = tempi_a, tempj = tempj_a;
			if (crystal_size_z == 0.f)
				temp_ijk = tempk * d_Nxy + convert_uint_sat(tempj) * d_Nx;
			else
				temp_ijk = convert_uint_sat(tempj) * d_Nx;
			if (local_sino > 0.f) {
				nominator_mlem(ax, local_sino, d_epps, temp, d_randoms, d_sc_ra, idx, lid, n_rekos);
				for (uint ii = 0u; ii < Np_n; ii++) {
					if (tx0 < ty0) {
						if (ii == Np_n - 1u) {
							if (crystal_size_z == 0.f) {
								orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
									local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM,  d_rhs_MLEM,  d_Summ, true, false);
							}
							else {
								orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
									local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, true, false);
							}
						}
						else {
							tempi += iu;
							tx0 += txu;
						}
					}
					else {
						if (crystal_size_z == 0.f) {
							orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
									local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM,  d_rhs_MLEM,  d_Summ, true, false);
						}
						else {
							orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
									local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, true, false);
						}
						if (ju > 0)
							temp_ijk += d_Nx;
						else
							temp_ijk -= d_Nx;
						ty0 += tyu;
						tempj += ju;
					}
				}
			}
			else {
				for (uint ii = 0u; ii < Np_n; ii++) {
					if (tx0 < ty0) {
						if (ii == Np_n - 1u) {
							if (crystal_size_z == 0.f) {
								orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM,  d_rhs_MLEM, d_Summ, false, true);
							}
							else {
								orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
							}
						}
						else {
							tempi += iu;
							tx0 += txu;
						}
					}
					else {
						if (crystal_size_z == 0.f) {
							orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM,  d_rhs_MLEM, d_Summ, false, true);
						}
						else {
							orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
						}
						if (ju > 0)
							temp_ijk += d_Nx;
						else
							temp_ijk -= d_Nx;
						ty0 += tyu;
						tempj += ju;
					}
				}
			}
		}
	}
	else {
		if (fabs(y_diff) < 1e-6f) {
			if (yd <= d_maxyy && yd >= d_by) {
				int tempi = 0, tempk = 0, tempj = 0, iu = 0, ku = 0;
				float txu = 0.f, tzu = 0.f, tc = 0.f, tx0 = 0.f, tz0 = 0.f;
				const bool skip = siddon_pre_loop_2D(d_bx, d_bz, x_diff, z_diff, d_maxxx, d_bzb, d_dx, d_dz, d_Nx, d_Nz, &tempi, &tempk, &txu, &tzu, &Np, TYPE1,
					zs, xs, zd, xd, &tc, &iu, &ku, &tx0, &tz0);
				if (skip)
					return;
				tempj = perpendicular_start(d_by, yd, d_dy, d_Ny);
				float tx0_a = tx0, tz0_a = tz0;
				int tempi_a = tempi, tempk_a = tempk;
				uint temp_ijk;
				if (crystal_size_z == 0.f)
					temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempi);
				else {
					temp_ijk = convert_uint_sat(tempi);
					//const float temp_x = xs;
					//xs = ys;
					//ys = temp_x;
				}
				if (d_attenuation_correction == 1u)
					LL = native_sqrt(x_diff * x_diff + z_diff * z_diff);
				for (uint ii = 0u; ii < Np; ii++) {
					if (tx0 < tz0) {
						if (d_attenuation_correction == 1u)
							compute_attenuation(&tc, &jelppi, LL, tx0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
						if (crystal_size_z == 0.f) {
							orth_distance_mlem(tempj, d_Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Nx, d_Nx, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
						}
						else {
							orth_distance_mlem_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Nx, d_Nx, tempk, d_Nxy, ys, xs, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
						}
						if (iu > 0)
							temp_ijk++;
						else
							temp_ijk--;
						tempi += iu;
						tx0 += txu;
						xyz = 1u;
					}
					else {
						if (d_attenuation_correction == 1u)
							compute_attenuation(&tc, &jelppi, LL, tz0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
						if (crystal_size_z == 0.f) {
							orth_distance_mlem(tempj, d_Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Nx, d_Nx, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
							if (ku > 0)
								temp_ijk += d_Nxy;
							else
								temp_ijk -= d_Nxy;
						}
						else if (ii == Np - 1u) {
							orth_distance_mlem_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Nx, d_Nx, tempk, d_Nxy, ys, xs, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
						}
						tempk += ku;
						tz0 += tzu;
						xyz = 3u;
					}
					Np_n++;
					if (tempk < 0 || tempi < 0 || tempi >= d_Nx || tempk >= d_Nz) {
						if (crystal_size_z != 0.f && xyz == 3u) {
							orth_distance_mlem_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, &temp, temp_ijk,
								local_sino, ax, d_Nx, d_Nx, tempk, d_Nxy, ys, xs, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
						}
						break;
					}
				}
				temp = 1.f / temp;
				if (d_attenuation_correction == 1u)
					temp *= native_exp(jelppi);
				if (d_normalization == 1u)
					temp *= d_norm[idx];
				tx0 = tx0_a, tz0 = tz0_a;
				tempi = tempi_a, tempk = tempk_a;
				if (crystal_size_z == 0.f)
					temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempi);
				else
					temp_ijk = convert_uint_sat(tempi);
				if (local_sino > 0.f) {
					nominator_mlem(ax, local_sino, d_epps, temp, d_randoms, d_sc_ra, idx, lid, n_rekos);
					for (uint ii = 0u; ii < Np_n; ii++) {
						if (tx0 < tz0) {
							if (crystal_size_z == 0.f) {
								orth_distance_mlem(tempj, d_Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, &temp, temp_ijk, 
									local_sino, ax, d_Nx, d_Nx, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, true, false);
							}
							else {
								orth_distance_mlem_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Nx, d_Nx, tempk, d_Nxy, ys, xs, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, true, false);
							}
							if (iu > 0)
								temp_ijk++;
							else
								temp_ijk--;
							tempi += iu;
							tx0 += txu;
						}
						else {
							if (crystal_size_z == 0.f) {
								orth_distance_mlem(tempj, d_Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, &temp, temp_ijk, 
									local_sino, ax, d_Nx, d_Nx, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, true, false);
								if (ku > 0)
									temp_ijk += d_Nxy;
								else
									temp_ijk -= d_Nxy;
							}
							else if (ii == Np_n - 1u) {
								orth_distance_mlem_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Nx, d_Nx, tempk, d_Nxy, ys, xs, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, true, false);
							}
							tz0 += tzu;
							tempk += ku;
						}
					}
				}
				else {
					for (uint ii = 0u; ii < Np_n; ii++) {
						if (tx0 < tz0) {
							if (crystal_size_z == 0.f) {
								orth_distance_mlem(tempj, d_Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, &temp, temp_ijk, 
									local_sino, ax, d_Nx, d_Nx, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, true);
							}
							else {
								orth_distance_mlem_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Nx, d_Nx, tempk, d_Nxy, ys, xs, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, true);
							}
							if (iu > 0)
								temp_ijk++;
							else
								temp_ijk--;
							tempi += iu;
						}
						else {
							if (crystal_size_z == 0.f) {
								orth_distance_mlem(tempj, d_Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, &temp, temp_ijk, 
									local_sino, ax, d_Nx, d_Nx, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, true);
								if (ku > 0)
									temp_ijk += d_Nxy;
								else
									temp_ijk -= d_Nxy;
							}
							else if (ii == Np_n - 1u) {
								orth_distance_mlem_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Nx, d_Nx, tempk, d_Nxy, ys, xs, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, true);
							}
							tempk += ku;
							tz0 += tzu;
						}
					}
				}
			}
		}
		else if (fabs(x_diff) < 1e-6f) {
			if (xd <= d_maxxx && xd >= d_bx) {
				int tempi = 0, tempk = 0, tempj = 0, ju = 0, ku = 0;
				float tyu = 0.f, tzu = 0.f, tc = 0.f, ty0 = 0.f, tz0 = 0.f;
				const bool skip = siddon_pre_loop_2D(d_by, d_bz, y_diff, z_diff, d_maxyy, d_bzb, d_dy, d_dz, d_Ny, d_Nz, &tempj, &tempk, &tyu, &tzu, &Np, TYPE1,
					zs, ys, zd, yd, &tc, &ju, &ku, &ty0, &tz0);
				if (skip)
					return;
				tempi = perpendicular_start(d_bx, xd, d_dx, d_Nx);
				const float tz0_a = tz0, ty0_a = ty0;
				const int tempj_a = tempj, tempk_a = tempk;
				uint temp_ijk;
				if (crystal_size_z == 0.f)
					temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempj) * d_Nx;
				else
					temp_ijk = convert_uint_sat(tempj) * d_Nx;
				if (d_attenuation_correction == 1u)
					LL = native_sqrt(y_diff * y_diff + z_diff * z_diff);
				for (uint ii = 0u; ii < Np; ii++) {
					if (tz0 < ty0) {
						if (d_attenuation_correction == 1u)
							compute_attenuation(&tc, &jelppi, LL, tz0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
						if (crystal_size_z == 0.f) {
							orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
							if (ku > 0)
								temp_ijk += d_Nxy;
							else
								temp_ijk -= d_Nxy;
						}
						else if (ii == Np - 1u) {
							orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
						}
						tempk += ku;
						tz0 += tzu;
						xyz = 3u;
					}
					else {
						if (d_attenuation_correction == 1u)
							compute_attenuation(&tc, &jelppi, LL, ty0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
						if (crystal_size_z == 0.f) {
							orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
						}
						else {
							orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
						}
						if (ju > 0)
							temp_ijk += d_Nx;
						else
							temp_ijk -= d_Nx;
						ty0 += tyu;
						tempj += (ju);
						xyz = 2u;
					}
					Np_n++;
					if (tempj < 0 || tempk < 0 || tempk >= d_Nz || tempj >= d_Ny) {
						if (xyz == 3u && crystal_size_z != 0.f) {
							orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
						}
						break;
					}
				}
				temp = 1.f / temp;
				if (d_attenuation_correction == 1u)
					temp *= native_exp(jelppi);
				if (d_normalization == 1u)
					temp *= d_norm[idx];
				tz0 = tz0_a, ty0 = ty0_a;
				tempj = tempj_a, tempk = tempk_a;
				if (crystal_size_z == 0.f)
					temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempj) * d_Nx;
				else
					temp_ijk = convert_uint_sat(tempj) * d_Nx;
				if (local_sino > 0.f) {
					nominator_mlem(ax, local_sino, d_epps, temp, d_randoms, d_sc_ra, idx, lid, n_rekos);
					for (uint ii = 0u; ii < Np_n; ii++) {
						if (tz0 < ty0) {
							if (crystal_size_z == 0.f) {
								orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
									local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM,  d_rhs_MLEM,  d_Summ, true, false);
								if (ku > 0)
									temp_ijk += d_Nxy;
								else
									temp_ijk -= d_Nxy;
							}
							else if (ii == Np_n - 1u) {
								orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
									local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, true, false);
							}
							tempk += ku;
							tz0 += tzu;
						}
						else {
							if (crystal_size_z == 0.f) {
								orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
									local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM,  d_rhs_MLEM,  d_Summ, true, false);
							}
							else {
								orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
									local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, true, false);
							}
							if (ju > 0)
								temp_ijk += d_Nx;
							else
								temp_ijk -= d_Nx;
							ty0 += tyu;
							tempj += (ju);
						}
					}
				}
				else {
					for (uint ii = 0u; ii < Np_n; ii++) {
						if (tz0 < ty0) {
							if (crystal_size_z == 0.f) {
								orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM,  d_rhs_MLEM, d_Summ, false, true);
								if (ku > 0)
									temp_ijk += d_Nxy;
								else
									temp_ijk -= d_Nxy;
							}
							else if (ii == Np_n - 1u) {
								orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
									local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, true);
							}
							tempk += ku;
							tz0 += tzu;
						}
						else {
							if (crystal_size_z == 0.f) {
								orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM,  d_rhs_MLEM, d_Summ, false, true);
							}
							else {
								orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
									local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, true);
							}
							if (ju > 0)
								temp_ijk += d_Nx;
							else
								temp_ijk -= d_Nx;
							ty0 += tyu;
							tempj += (ju);
						}
					}
				}
			}
		}
		else {
			int tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;
			float txu = 0.f, tyu = 0.f, tzu = 0.f, tc = 0.f, tx0 = 0.f, ty0 = 0.f, tz0 = 0.f;
			const bool skip = siddon_pre_loop_3D(d_bx, d_by, d_bz, x_diff, y_diff, z_diff, d_maxxx, d_maxyy, d_bzb, d_dx, d_dy, d_dz, d_Nx, d_Ny, d_Nz, &tempi,
				&tempj, &tempk, &tyu, &txu, &tzu, &Np, TYPE1, ys, xs, yd, xd, zs, zd, &tc, &iu, &ju, &ku, &tx0, &ty0, &tz0);
			if (skip)
				return;
			float temp = 0.f;
			const float tx0_a = tx0, ty0_a = ty0, tz0_a = tz0;
			const int tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
			uint temp_ijk;
			if (crystal_size_z == 0.f)
				temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempj) * d_Nx;
			else
				temp_ijk = convert_uint_sat(tempj) * d_Nx;
			if (d_attenuation_correction == 1u)
				LL = native_sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
			for (uint ii = 0u; ii < Np; ii++) {
				if (tz0 < ty0 && tz0 < tx0) {
					if (d_attenuation_correction == 1u)
						compute_attenuation(&tc, &jelppi, LL, tz0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
					if (crystal_size_z == 0.f) {
						orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
						if (ku > 0)
							temp_ijk += d_Nxy;
						else
							temp_ijk -= d_Nxy;
					}
					else if (ii == Np - 1u) {
						orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
					}
					tempk += ku;
					tz0 += tzu;
					xyz = 3u;
				}
				else if (ty0 < tx0) {
					if (d_attenuation_correction == 1u)
						compute_attenuation(&tc, &jelppi, LL, ty0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
					if (crystal_size_z == 0.f) {
						orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
					}
					else {
						orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
					}
					if (ju > 0)
						temp_ijk += d_Nx;
					else
						temp_ijk -= d_Nx;
					tempj += (ju);
					ty0 += tyu;
					xyz = 2u;
				}
				else {
					if (d_attenuation_correction == 1u)
						compute_attenuation(&tc, &jelppi, LL, tx0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
					if (ii == Np - 1) {
						if (crystal_size_z == 0.f) {
							orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
						}
						else {
							orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
						}
					}
					else {
						tempi += iu;
						tx0 += txu;
					}
					xyz = 1u;
				}
				Np_n++;
				if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= d_Nx || tempj >= d_Ny || tempk >= d_Nz) {
					if (xyz == 1u || (xyz == 3u && crystal_size_z != 0.f)) {
						if (crystal_size_z == 0.f) {
							orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk,
								local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
						}
						else {
							orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
								local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, false);
						}
					}
					break;
				}
			}
			temp = 1.f / temp;
			if (d_attenuation_correction == 1u)
				temp *= native_exp(jelppi);
			if (d_normalization == 1u)
				temp *= d_norm[idx];
			tx0 = tx0_a, ty0 = ty0_a, tz0 = tz0_a;
			tempi = tempi_a, tempj = tempj_a, tempk = tempk_a;
			if (crystal_size_z == 0.f)
				temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempj) * d_Nx;
			else
				temp_ijk = convert_uint_sat(tempj) * d_Nx;
			if (local_sino > 0.f) {
				nominator_mlem(ax, local_sino, d_epps, temp, d_randoms, d_sc_ra, idx, lid, n_rekos);
				for (uint ii = 0u; ii < Np_n; ii++) {
					if (tz0 < ty0 && tz0 < tx0) {
						if (crystal_size_z == 0.f) {
							orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
									local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM,  d_rhs_MLEM,  d_Summ, true, false);
							if (ku > 0)
								temp_ijk += d_Nxy;
							else
								temp_ijk -= d_Nxy;
						}
						else if (ii == Np_n - 1u) {
							orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
									local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, true, false);
						}
						tempk += ku;
						tz0 += tzu;
					}
					else if (ty0 < tx0) {
						if (crystal_size_z == 0.f) {
							orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
									local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM,  d_rhs_MLEM,  d_Summ, true, false);
						}
						else {
							orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
									local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, true, false);
						}
						if (ju > 0)
							temp_ijk += d_Nx;
						else
							temp_ijk -= d_Nx;
						tempj += (ju);
						ty0 += tyu;
					}
					else {
						if (ii == Np_n - 1) {
							if (crystal_size_z == 0.f) {
								orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
									local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM,  d_rhs_MLEM,  d_Summ, true, false);
							}
							else {
								orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
									local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, true, false);
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
				for (uint ii = 0u; ii < Np_n; ii++) {
					if (tz0 < ty0 && tz0 < tx0) {
						if (crystal_size_z == 0.f) {
							orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM,  d_rhs_MLEM, d_Summ, false, true);
							if (ku > 0)
								temp_ijk += d_Nxy;
							else
								temp_ijk -= d_Nxy;
						}
						else if (ii == Np_n - 1u) {
							orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
									local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, true);
						}
						tempk += ku;
						tz0 += tzu;
					}
					else if (ty0 < tx0) {
						if (crystal_size_z == 0.f) {
							orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM,  d_rhs_MLEM, d_Summ, false, true);
						}
						else {
							orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
									local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, true);
						}
						if (ju > 0)
							temp_ijk += d_Nx;
						else
							temp_ijk -= d_Nx;
						tempj += (ju);
						ty0 += tyu;
					}
					else {
						if (ii == Np_n - 1) {
							if (crystal_size_z == 0.f) {
								orth_distance_mlem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								local_sino, ax, d_Ny, 1u, no_norm, lid, n_rekos, d_N, d_MLEM,  d_rhs_MLEM, d_Summ, false, true);
							}
							else {
								orth_distance_mlem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk, 
									local_sino, ax, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, no_norm, lid, n_rekos, d_N, d_MLEM, d_rhs_MLEM, d_Summ, false, true);
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




__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(64, 1, 1)))
void MRAMLA_prepass_orth_noprecomp(const uint d_det_per_ring, const uchar d_raw, __constant uint* d_pseudos, const uint d_pRows, const float d_h,
	const uint d_Nx, const uint d_Ny, const uint d_Nz, const float d_dz, const float d_dx, const float d_dy, const float d_bz, const float d_bx,
	const float d_by, const float d_bzb, const float d_maxxx, const float d_maxyy, const float d_zmax, const float d_NSlices, const __global float* d_x,
	const __global float* d_y, const __global float* d_zdet, const uint d_size_x, const uint d_TotSinos, const uint d_attenuation_correction,
	const uint d_normalization, const uint d_randoms, const __global float* d_atten, const __global float* d_norm, const float d_epps, const uint d_N,
	const RecMethodsOpenCL MethodList, const uint d_Nxy, const float tube_width_xy, const float crystal_size_z, const int dec, __constant float* x_center,
	__constant float* y_center, __constant float* z_center, const __global uint* d_xyindex, const __global ushort* d_zindex, const __global float* d_COSEM,
	const __global float* d_Sino, const __global float* d_sc_ra, const uint d_alku, const __global ushort* d_L, const uchar MBSREM_prepass,
	__global float* d_ACOSEM_lhs, __global float* d_Amin, __global float* d_co, __global float* d_aco, __global float* d_Summ, const __global ushort* d_lor,
	__global float* d_E, const ulong m_size) {
	uint idx = get_global_id(0);
	if (idx >= m_size)
		return;
	float xs, xd, ys, yd, zs, zd;
	if (d_raw) {
		get_detector_coordinates_raw(d_x, d_y, d_zdet, d_L, d_det_per_ring, idx, d_pseudos, d_pRows, &xs, &xd, &ys, &yd, &zs, &zd);
	}
	// Sinogram data
	else {
		get_detector_coordinates(d_xyindex, d_zindex, d_size_x, idx, d_TotSinos, &xs, &xd, &ys, &yd, &zs, &zd, d_x, d_y, d_zdet);
	}
	float y_diff = (yd - ys);
	float x_diff = (xd - xs);
	float z_diff = (zd - zs);
	if ((y_diff == 0.f && x_diff == 0.f && z_diff == 0.f) || (y_diff == 0.f && x_diff == 0.f))
		return;
	uint Np = 0u;
	uint Np_n = 0u;
	uint xyz = 0u;
	const float local_sino = (d_Sino[idx]);
	float jelppi = 0.f, LL;
	float axCOSEM = 0.f;
	float axACOSEM = 0.f;
	float minimi = 1e8f;
	float kerroin, length_;
	if (crystal_size_z == 0.f) {
		kerroin = xd * ys - yd * xs;
		length_ = native_sqrt(y_diff * y_diff + x_diff * x_diff) * tube_width_xy;
	}
	else
		kerroin = e_norm(x_diff, y_diff, z_diff) * crystal_size_z;
	if (fabs(z_diff) < 1e-8f) {
		const uint z_loop = convert_uint((zs / d_zmax) * (d_NSlices - 1.f));
		if (fabs(y_diff) < 1e-8f) {
			if (yd <= d_maxyy && yd >= d_by) {
				float temp = 0.f;
				if (crystal_size_z == 0.f) {
					orth_distance_perpendicular_cosem(-x_diff, y_center, kerroin, length_, &temp, d_attenuation_correction, d_normalization, d_by, yd, d_dy, d_Ny,
						d_Nx, z_loop, d_atten, d_norm, idx, d_Ny, local_sino, 1u, MethodList, d_alku, &axCOSEM, d_COSEM, d_h, d_E,
						d_co, d_aco, &minimi, MBSREM_prepass, d_Summ, &axACOSEM, d_sc_ra, d_randoms, d_Amin, d_ACOSEM_lhs, false);
					if (d_alku == 0 && (MethodList.COSEM == 1 || MethodList.ECOSEM == 1 || MethodList.ACOSEM == 1
						|| MethodList.OSLCOSEM > 0) && local_sino != 0.f) {
						nominator_cosem(&axCOSEM, local_sino, d_epps, temp, d_randoms, d_sc_ra, idx);
					}
					orth_distance_perpendicular_cosem(-x_diff, y_center, kerroin, length_, &temp, d_attenuation_correction, d_normalization, d_by, yd, d_dy, d_Ny,
						d_Nx, z_loop, d_atten, d_norm, idx, d_Ny, local_sino, 1u, MethodList, d_alku, &axCOSEM, d_COSEM, d_h, d_E,
						d_co, d_aco, &minimi, MBSREM_prepass, d_Summ, &axACOSEM, d_sc_ra, d_randoms, d_Amin, d_ACOSEM_lhs, true);
				}
				else {
					//const float temppi = xs;
					//xs = ys;
					//ys = temppi;
					orth_distance_perpendicular_cosem_3D(y_center, x_center[0], z_center, &temp, d_attenuation_correction, d_normalization, d_by, yd, d_dy, d_Ny,
						d_Nx, z_loop, d_atten, d_norm, idx, d_Ny, local_sino, 1u, ys, xs, zs, y_diff, x_diff, z_diff, kerroin, d_Nxy, d_Nz,
						MethodList, d_alku, &axCOSEM, d_COSEM, d_h, d_E, d_co, d_aco, &minimi, MBSREM_prepass, d_Summ, &axACOSEM, d_sc_ra, d_randoms, d_Amin, 
						d_ACOSEM_lhs, false);
					if (d_alku == 0 && (MethodList.COSEM == 1 || MethodList.ECOSEM == 1 || MethodList.ACOSEM == 1
						|| MethodList.OSLCOSEM > 0) && local_sino != 0.f) {
						nominator_cosem(&axCOSEM, local_sino, d_epps, temp, d_randoms, d_sc_ra, idx);
					}
					orth_distance_perpendicular_cosem_3D(y_center, x_center[0], z_center, &temp, d_attenuation_correction, d_normalization, d_by, yd, d_dy, d_Ny,
						d_Nx, z_loop, d_atten, d_norm, idx, d_Ny, local_sino, 1u, ys, xs, zs, y_diff, x_diff, z_diff, kerroin, d_Nxy, d_Nz,
						MethodList, d_alku, &axCOSEM, d_COSEM, d_h, d_E, d_co, d_aco, &minimi, MBSREM_prepass, d_Summ, &axACOSEM, d_sc_ra, d_randoms, d_Amin, 
						d_ACOSEM_lhs, true);
				}
			}
		}
		else if (fabs(x_diff) < 1e-8f) {
			if (xd <= d_maxxx && xd >= d_bx) {
				float temp = 0.f;
				if (crystal_size_z == 0.f) {
					orth_distance_perpendicular_cosem(y_diff, y_center, kerroin, length_, &temp, d_attenuation_correction, d_normalization, d_bx, xd, d_dx, d_Nx,
						d_Ny, z_loop, d_atten, d_norm, idx, 1u, local_sino, d_Nx, MethodList, d_alku, &axCOSEM, d_COSEM, d_h, d_E,
						d_co, d_aco, &minimi, MBSREM_prepass, d_Summ, &axACOSEM, d_sc_ra, d_randoms, d_Amin, d_ACOSEM_lhs, false);
					if (d_alku == 0 && (MethodList.COSEM == 1 || MethodList.ECOSEM == 1 || MethodList.ACOSEM == 1 || MethodList.OSLCOSEM > 0)
						&& local_sino != 0.f) {
						nominator_cosem(&axCOSEM, local_sino, d_epps, temp, d_randoms, d_sc_ra, idx);
					}
					orth_distance_perpendicular_cosem(y_diff, y_center, kerroin, length_, &temp, d_attenuation_correction, d_normalization, d_bx, xd, d_dx, d_Nx,
						d_Ny, z_loop, d_atten, d_norm, idx, 1u, local_sino, d_Nx, MethodList, d_alku, &axCOSEM, d_COSEM, d_h, d_E,
						d_co, d_aco, &minimi, MBSREM_prepass, d_Summ, &axACOSEM, d_sc_ra, d_randoms, d_Amin, d_ACOSEM_lhs, true);
				}
				else {
					orth_distance_perpendicular_cosem_3D(x_center, y_center[0], z_center, &temp, d_attenuation_correction, d_normalization, d_bx, xd, d_dx, d_Nx,
						d_Ny, z_loop, d_atten, d_norm, idx, 1u, local_sino, d_Nx, xs, ys, zs, x_diff, y_diff, z_diff, kerroin, d_Nxy, d_Nz,
						MethodList, d_alku, &axCOSEM, d_COSEM, d_h, d_E, d_co, d_aco, &minimi, MBSREM_prepass, d_Summ, &axACOSEM, d_sc_ra, d_randoms, d_Amin, 
						d_ACOSEM_lhs, false);
					if (d_alku == 0 && (MethodList.COSEM == 1 || MethodList.ECOSEM == 1 || MethodList.ACOSEM == 1 || MethodList.OSLCOSEM > 0)
						&& local_sino != 0.f) {
						nominator_cosem(&axCOSEM, local_sino, d_epps, temp, d_randoms, d_sc_ra, idx);
					}
					orth_distance_perpendicular_cosem_3D(x_center, y_center[0], z_center, &temp, d_attenuation_correction, d_normalization, d_bx, xd, d_dx, d_Nx,
						d_Ny, z_loop, d_atten, d_norm, idx, 1u, local_sino, d_Nx, xs, ys, zs, x_diff, y_diff, z_diff, kerroin, d_Nxy, d_Nz,
						MethodList, d_alku, &axCOSEM, d_COSEM, d_h, d_E, d_co, d_aco, &minimi, MBSREM_prepass, d_Summ, &axACOSEM, d_sc_ra, d_randoms, d_Amin, 
						d_ACOSEM_lhs, true);
				}
			}
		}
		else {
			int tempi = 0, tempj = 0, iu = 0, ju = 0;
			float txu = 0.f, tyu = 0.f, tc = 0.f, tx0 = 0.f, ty0 = 0.f;
			const bool skip = siddon_pre_loop_2D(d_bx, d_by, x_diff, y_diff, d_maxxx, d_maxyy, d_dx, d_dy, d_Nx, d_Ny, &tempi, &tempj, &txu, &tyu, &Np, TYPE1,
				ys, xs, yd, xd, &tc, &iu, &ju, &tx0, &ty0);
			if (skip)
				return;
			float temp = 0.f;
			float tx0_a = tx0, ty0_a = ty0;
			int tempi_a = tempi, tempj_a = tempj;
			int tempk = z_loop;
			uint temp_ijk;
			if (crystal_size_z == 0.f)
				temp_ijk = tempk * d_Nxy + convert_uint_sat(tempj) * d_Nx;
			else
				temp_ijk = convert_uint_sat(tempj) * d_Nx;
			if (d_attenuation_correction == 1u)
				LL = native_sqrt(x_diff * x_diff + y_diff * y_diff);
			for (uint ii = 0u; ii < Np; ii++) {
				if (tx0 < ty0) {
					if (d_attenuation_correction == 1u)
						compute_attenuation(&tc, &jelppi, LL, tx0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
					if (ii == Np - 1u) {
						if (crystal_size_z == 0.f) {
							orth_distance_cosem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, MethodList,
								d_alku, local_sino, d_Ny, 1u, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, false);
						}
						else {
							orth_distance_cosem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
								MethodList, d_alku, local_sino, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, false);
						}
					}
					else {
						tempi += iu;
						tx0 += txu;
					}
					xyz = 1u;
				}
				else {
					if (d_attenuation_correction == 1u)
						compute_attenuation(&tc, &jelppi, LL, ty0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
					if (crystal_size_z == 0.f) {
						orth_distance_cosem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, MethodList,
							d_alku, local_sino, d_Ny, 1u, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
							&minimi, MBSREM_prepass, d_Summ, &axACOSEM, false);
					}
					else {
						orth_distance_cosem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
							MethodList, d_alku, local_sino, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
							&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, false);
					}
					tempj += ju;
					if (ju > 0)
						temp_ijk += d_Nx;
					else
						temp_ijk -= d_Nx;
					ty0 += tyu;
					xyz = 2u;
				}
				Np_n++;
				if (tempj < 0 || tempi < 0 || tempi >= d_Nx || tempj >= d_Ny) {
					if (xyz == 1u) {
						if (crystal_size_z == 0.f) {
							orth_distance_cosem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, MethodList,
								d_alku, local_sino, d_Ny, 1u, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, false);
						}
						else {
							orth_distance_cosem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
								MethodList, d_alku, local_sino, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, false);
						}
					}
					break;
				}
			}
			temp = 1.f / temp;
			if (d_attenuation_correction == 1u)
				temp *= native_exp(jelppi);
			if (d_normalization == 1u)
				temp *= d_norm[idx];
			tx0 = tx0_a, ty0 = ty0_a;
			tempi = tempi_a, tempj = tempj_a;
			if (crystal_size_z == 0.f)
				temp_ijk = tempk * d_Nxy + convert_uint_sat(tempj) * d_Nx;
			else
				temp_ijk = convert_uint_sat(tempj) * d_Nx;
			if (d_alku == 0 && (MethodList.COSEM == 1 || MethodList.ECOSEM == 1 || MethodList.ACOSEM == 1 || MethodList.OSLCOSEM > 0)
				&& local_sino != 0.f) {
				if (axCOSEM == 0.f)
					axCOSEM = d_epps;
				else
					axCOSEM *= temp;
				if (d_randoms == 1u)
					axCOSEM += d_sc_ra[idx];
				axCOSEM = local_sino / axCOSEM;
			}
			for (uint ii = 0u; ii < Np_n; ii++) {
				if (tx0 < ty0) {
					if (ii == Np_n - 1u) {
						if (crystal_size_z == 0.f) {
							orth_distance_cosem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, MethodList,
								d_alku, local_sino, d_Ny, 1u, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, true);
						}
						else {
							orth_distance_cosem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
								MethodList, d_alku, local_sino, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, true);
						}
					}
					else {
						tempi += iu;
						tx0 += txu;
					}
				}
				else {
					if (crystal_size_z == 0.f) {
						orth_distance_cosem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, MethodList,
							d_alku, local_sino, d_Ny, 1u, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
							&minimi, MBSREM_prepass, d_Summ, &axACOSEM, true);
					}
					else {
						orth_distance_cosem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
							MethodList, d_alku, local_sino, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
							&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, true);
					}
					tempj += ju;
					if (ju > 0)
						temp_ijk += d_Nx;
					else
						temp_ijk -= d_Nx;
					ty0 += tyu;
				}
			}
			if ((MethodList.MRAMLA == 1 || MethodList.MBSREM == 1) && MBSREM_prepass == 1)
				d_Amin[idx] = minimi;
			if ((MethodList.ACOSEM == 1 || MethodList.OSLCOSEM == 1) && d_alku > 0) {
				if (d_randoms == 1u)
					axACOSEM += d_sc_ra[idx];
				d_ACOSEM_lhs[idx] = axACOSEM;
			}
		}
	}
	else {
		if (fabs(y_diff) < 1e-8f) {
			if (yd <= d_maxyy && yd >= d_by) {
				int tempi = 0, tempk = 0, tempj = 0, iu = 0, ku = 0;
				float txu = 0.f, tzu = 0.f, tc = 0.f, tx0 = 0.f, tz0 = 0.f;
				const bool skip = siddon_pre_loop_2D(d_bx, d_bz, x_diff, z_diff, d_maxxx, d_bzb, d_dx, d_dz, d_Nx, d_Nz, &tempi, &tempk, &txu, &tzu, &Np, TYPE1,
					zs, xs, zd, xd, &tc, &iu, &ku, &tx0, &tz0);
				if (skip)
					return;
				tempj = perpendicular_start(d_by, yd, d_dy, d_Ny);
				float temp = 0.f;
				float tx0_a = tx0, tz0_a = tz0;
				int tempi_a = tempi, tempk_a = tempk;
				float minimi = 1e8f;
				uint temp_ijk;
				if (crystal_size_z == 0.f)
					temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempi);
				else {
					temp_ijk = convert_uint_sat(tempi);
				}
				if (d_attenuation_correction == 1u)
					LL = native_sqrt(x_diff * x_diff + z_diff * z_diff);
				for (uint ii = 0u; ii < Np; ii++) {
					if (tx0 < tz0) {
						if (d_attenuation_correction == 1u)
							compute_attenuation(&tc, &jelppi, LL, tx0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
						if (crystal_size_z == 0.f) {
							orth_distance_cosem(tempj, d_Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, &temp, temp_ijk, MethodList,
								d_alku, local_sino, d_Nx, d_Nx, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, false);
						}
						else {
							orth_distance_cosem_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, &temp, temp_ijk,
								MethodList, d_alku, local_sino, d_Nx, d_Nx, tempk, d_Nxy, ys, xs, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, false);
						}
						if (iu > 0)
							temp_ijk++;
						else
							temp_ijk--;
						tempi += iu;
						tx0 += txu;
						xyz = 1u;
					}
					else {
						if (d_attenuation_correction == 1u)
							compute_attenuation(&tc, &jelppi, LL, tz0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
						if (crystal_size_z == 0.f) {
							orth_distance_cosem(tempj, d_Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, &temp, temp_ijk, MethodList,
								d_alku, local_sino, d_Nx, d_Nx, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, false);
							if (ku > 0)
								temp_ijk += d_Nxy;
							else
								temp_ijk -= d_Nxy;
						}
						else if (ii == Np - 1u) {
							orth_distance_cosem_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, &temp, temp_ijk,
								MethodList, d_alku, local_sino, d_Nx, d_Nx, tempk, d_Nxy, ys, xs, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, false);
						}
						tempk += ku;
						tz0 += tzu;
						xyz = 3u;
					}
					Np_n++;
					if (tempk < 0 || tempi < 0 || tempi >= d_Nx || tempk >= d_Nz) {
						if (crystal_size_z != 0.f && xyz == 3u) {
							orth_distance_cosem_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, &temp, temp_ijk,
								MethodList, d_alku, local_sino, d_Nx, d_Nx, tempk, d_Nxy, ys, xs, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, false);
						}
						break;
					}
				}
				temp = 1.f / temp;
				if (d_attenuation_correction == 1u)
					temp *= native_exp(jelppi);
				if (d_normalization == 1u)
					temp *= d_norm[idx];
				tx0 = tx0_a, tz0 = tz0_a;
				tempi = tempi_a, tempk = tempk_a;
				if (crystal_size_z == 0.f)
					temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempi);
				else
					temp_ijk = convert_uint_sat(tempi);
				if (d_alku == 0 && (MethodList.COSEM == 1 || MethodList.ECOSEM == 1 || MethodList.ACOSEM == 1 || MethodList.OSLCOSEM > 0)
					&& local_sino != 0.f) {
					if (axCOSEM == 0.f)
						axCOSEM = d_epps;
					else
						axCOSEM *= temp;
					if (d_randoms == 1u)
						axCOSEM += d_sc_ra[idx];
					axCOSEM = local_sino / axCOSEM;
				}
				for (uint ii = 0u; ii < Np_n; ii++) {
					if (tx0 < tz0) {
						if (crystal_size_z == 0.f) {
							orth_distance_cosem(tempj, d_Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, &temp, temp_ijk, MethodList,
								d_alku, local_sino, d_Nx, d_Nx, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, true);
						}
						else {
							orth_distance_cosem_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, &temp, temp_ijk,
								MethodList, d_alku, local_sino, d_Nx, d_Nx, tempk, d_Nxy, ys, xs, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, true);
						}
						if (iu > 0)
							temp_ijk++;
						else
							temp_ijk--;
						tempi += iu;
						tx0 += txu;
					}
					else {
						if (crystal_size_z == 0.f) {
							orth_distance_cosem(tempj, d_Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, &temp, temp_ijk, MethodList,
								d_alku, local_sino, d_Nx, d_Nx, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, true);

							if (ku > 0)
								temp_ijk += d_Nxy;
							else
								temp_ijk -= d_Nxy;
						}
						else if (ii == Np_n - 1u) {
							orth_distance_cosem_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, &temp, temp_ijk,
								MethodList, d_alku, local_sino, d_Nx, d_Nx, tempk, d_Nxy, ys, xs, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, true);
						}
						tempk += ku;
						tz0 += tzu;
					}
				}
				if ((MethodList.MRAMLA == 1 || MethodList.MBSREM == 1) && MBSREM_prepass == 1)
					d_Amin[idx] = minimi;
				if ((MethodList.ACOSEM == 1 || MethodList.OSLCOSEM == 1) && d_alku > 0) {
					if (d_randoms == 1u)
						axACOSEM += d_sc_ra[idx];
					d_ACOSEM_lhs[idx] = axACOSEM;
				}
			}
		}
		else if (fabs(x_diff) < 1e-8f) {
			if (xd <= d_maxxx && xd >= d_bx) {
				int tempi = 0, tempk = 0, tempj = 0, ju = 0, ku = 0;
				float tyu = 0.f, tzu = 0.f, tc = 0.f, ty0 = 0.f, tz0 = 0.f;
				const bool skip = siddon_pre_loop_2D(d_by, d_bz, y_diff, z_diff, d_maxyy, d_bzb, d_dy, d_dz, d_Ny, d_Nz, &tempj, &tempk, &tyu, &tzu, &Np, TYPE1,
					zs, ys, zd, yd, &tc, &ju, &ku, &ty0, &tz0);
				if (skip)
					return;
				float temp = 0.f;
				tempi = perpendicular_start(d_bx, xd, d_dx, d_Nx);
				const float tz0_a = tz0, ty0_a = ty0;
				const int tempj_a = tempj, tempk_a = tempk;
				uint temp_ijk;
				if (crystal_size_z == 0.f)
					temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempj) * d_Nx;
				else
					temp_ijk = convert_uint_sat(tempj) * d_Nx;
				if (d_attenuation_correction == 1u)
					LL = native_sqrt(y_diff * y_diff + z_diff * z_diff);
				for (uint ii = 0u; ii < Np; ii++) {
					if (tz0 < ty0) {
						if (d_attenuation_correction == 1u)
							compute_attenuation(&tc, &jelppi, LL, tz0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
						if (crystal_size_z == 0.f) {
							orth_distance_cosem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, MethodList,
								d_alku, local_sino, d_Ny, 1u, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, false);
							if (ku > 0)
								temp_ijk += d_Nxy;
							else
								temp_ijk -= d_Nxy;
						}
						else if (ii == Np - 1u) {
							orth_distance_cosem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
								MethodList, d_alku, local_sino, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, false);
						}
						tempk += ku;
						tz0 += tzu;
						xyz = 3u;
					}
					else {
						if (d_attenuation_correction == 1u)
							compute_attenuation(&tc, &jelppi, LL, ty0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
						if (crystal_size_z == 0.f) {
							orth_distance_cosem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, MethodList,
								d_alku, local_sino, d_Ny, 1u, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, false);
						}
						else {
							orth_distance_cosem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
								MethodList, d_alku, local_sino, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, false);
						}
						tempj += ju;
						if (ju > 0)
							temp_ijk += d_Nx;
						else
							temp_ijk -= d_Nx;
						ty0 += tyu;
						xyz = 2u;
					}
					Np_n++;
					if (tempj < 0 || tempk < 0 || tempk >= d_Nz || tempj >= d_Ny) {
						if (xyz == 3u && crystal_size_z != 0.f) {
							orth_distance_cosem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
								MethodList, d_alku, local_sino, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, false);
						}
						break;
					}
				}
				temp = 1.f / temp;
				if (d_attenuation_correction == 1u)
					temp *= native_exp(jelppi);
				if (d_normalization == 1u)
					temp *= d_norm[idx];
				ty0 = ty0_a, tz0 = tz0_a;
				tempj = tempj_a, tempk = tempk_a;
				if (crystal_size_z == 0.f)
					temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempj) * d_Nx;
				else
					temp_ijk = convert_uint_sat(tempj) * d_Nx;
				float minimi = 1e8f;
				if (d_alku == 0 && (MethodList.COSEM == 1 || MethodList.ECOSEM == 1 || MethodList.ACOSEM == 1 || MethodList.OSLCOSEM > 0)
					&& local_sino != 0.f) {
					if (axCOSEM == 0.f)
						axCOSEM = d_epps;
					else
						axCOSEM *= temp;
					if (d_randoms == 1u)
						axCOSEM += d_sc_ra[idx];
					axCOSEM = local_sino / axCOSEM;
				}
				for (uint ii = 0u; ii < Np_n; ii++) {
					if (tz0 < ty0) {
						if (crystal_size_z == 0.f) {
							orth_distance_cosem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, MethodList,
								d_alku, local_sino, d_Ny, 1u, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, true);
							if (ku > 0)
								temp_ijk += d_Nxy;
							else
								temp_ijk -= d_Nxy;
						}
						else if (ii == Np_n - 1u) {
							orth_distance_cosem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
								MethodList, d_alku, local_sino, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, true);
						}
						tempk += ku;
						tz0 += tzu;
					}
					else {
						if (crystal_size_z == 0.f) {
							orth_distance_cosem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, MethodList,
								d_alku, local_sino, d_Ny, 1u, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, true);
						}
						else {
							orth_distance_cosem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
								MethodList, d_alku, local_sino, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, true);
						}
						tempj += ju;
						if (ju > 0)
							temp_ijk += d_Nx;
						else
							temp_ijk -= d_Nx;
						ty0 += tyu;
					}
				}
				if ((MethodList.MRAMLA == 1 || MethodList.MBSREM == 1) && MBSREM_prepass == 1)
					d_Amin[idx] = minimi;
				if ((MethodList.ACOSEM == 1 || MethodList.OSLCOSEM == 1) && d_alku > 0) {
					if (d_randoms == 1u)
						axACOSEM += d_sc_ra[idx];
					d_ACOSEM_lhs[idx] = axACOSEM;
				}
			}
		}
		else {
			int tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;
			float txu = 0.f, tyu = 0.f, tzu = 0.f, tc = 0.f, tx0 = 0.f, ty0 = 0.f, tz0 = 0.f;
			const bool skip = siddon_pre_loop_3D(d_bx, d_by, d_bz, x_diff, y_diff, z_diff, d_maxxx, d_maxyy, d_bzb, d_dx, d_dy, d_dz, d_Nx, d_Ny, d_Nz, &tempi,
				&tempj, &tempk, &tyu, &txu, &tzu, &Np, TYPE1, ys, xs, yd, xd, zs, zd, &tc, &iu, &ju, &ku, &tx0, &ty0, &tz0);
			if (skip)
				return;
			float temp = 0.f;
			const float tx0_a = tx0, ty0_a = ty0, tz0_a = tz0;
			const int tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
			float minimi = 1e8f;
			uint temp_ijk;
			if (crystal_size_z == 0.f)
				temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempj) * d_Nx;
			else
				temp_ijk = convert_uint_sat(tempj) * d_Nx;
			if (d_attenuation_correction == 1u)
				LL = native_sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
			for (uint ii = 0u; ii < Np; ii++) {
				if (tz0 < ty0 && tz0 < tx0) {
					if (d_attenuation_correction == 1u)
						compute_attenuation(&tc, &jelppi, LL, tx0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
					if (crystal_size_z == 0.f) {
						orth_distance_cosem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, MethodList,
							d_alku, local_sino, d_Ny, 1u, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
							&minimi, MBSREM_prepass, d_Summ, &axACOSEM, false);
						if (ku > 0)
							temp_ijk += d_Nxy;
						else
							temp_ijk -= d_Nxy;
					}
					else if (ii == Np - 1u) {
						orth_distance_cosem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
							MethodList, d_alku, local_sino, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
							&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, false);
					}
					tempk += ku;
					tz0 += tzu;
					xyz = 3u;
				}
				else if (tx0 < ty0) {
					if (d_attenuation_correction == 1u)
						compute_attenuation(&tc, &jelppi, LL, ty0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
					if (ii == Np - 1u) {
						if (crystal_size_z == 0.f) {
							orth_distance_cosem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, MethodList,
								d_alku, local_sino, d_Ny, 1u, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, false);
						}
						else {
							orth_distance_cosem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
								MethodList, d_alku, local_sino, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, false);
						}
					}
					else {
						tempi += iu;
						tx0 += txu;
					}
					xyz = 1u;
				}
				else {
					if (d_attenuation_correction == 1u)
						compute_attenuation(&tc, &jelppi, LL, tz0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
					if (crystal_size_z == 0.f) {
						orth_distance_cosem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, MethodList,
							d_alku, local_sino, d_Ny, 1u, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
							&minimi, MBSREM_prepass, d_Summ, &axACOSEM, false);
					}
					else {
						orth_distance_cosem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
							MethodList, d_alku, local_sino, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
							&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, false);
					}
					tempj += ju;
					if (ju > 0)
						temp_ijk += d_Nx;
					else
						temp_ijk -= d_Nx;
					ty0 += tyu;
					xyz = 2u;
				}
				Np_n++;
				if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= d_Nx || tempj >= d_Ny || tempk >= d_Nz) {
					if (xyz == 1u || (xyz == 3u && crystal_size_z != 0.f)) {
						if (crystal_size_z == 0.f) {
							orth_distance_cosem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, MethodList,
								d_alku, local_sino, d_Ny, 1u, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, false);
						}
						else {
							orth_distance_cosem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
								MethodList, d_alku, local_sino, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, false);
						}
					}
					break;
				}
			}
			temp = 1.f / temp;
			if (d_attenuation_correction == 1u)
				temp *= native_exp(jelppi);
			if (d_normalization == 1u)
				temp *= d_norm[idx];
			tx0 = tx0_a, ty0 = ty0_a, tz0 = tz0_a;
			tempi = tempi_a, tempj = tempj_a, tempk = tempk_a;
			if (crystal_size_z == 0.f)
				temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempj) * d_Nx;
			else
				temp_ijk = convert_uint_sat(tempj) * d_Nx;
			if (d_alku == 0 && (MethodList.COSEM == 1 || MethodList.ECOSEM == 1 || MethodList.ACOSEM == 1 || MethodList.OSLCOSEM > 0)
				&& local_sino != 0.f) {
				if (axCOSEM == 0.f)
					axCOSEM = d_epps;
				else
					axCOSEM *= temp;
				if (d_randoms == 1u)
					axCOSEM += d_sc_ra[idx];
				axCOSEM = local_sino / axCOSEM;
			}
			for (uint ii = 0u; ii < Np_n; ii++) {
				if (tz0 < ty0 && tz0 < tx0) {
					if (crystal_size_z == 0.f) {
						orth_distance_cosem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, MethodList,
							d_alku, local_sino, d_Ny, 1u, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
							&minimi, MBSREM_prepass, d_Summ, &axACOSEM, true);
						if (ku > 0)
							temp_ijk += d_Nxy;
						else
							temp_ijk -= d_Nxy;
					}
					else if (ii == Np_n - 1u) {
						orth_distance_cosem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
							MethodList, d_alku, local_sino, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
							&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, true);
					}
					tempk += ku;
					tz0 += tzu;
				}
				else if (tx0 < ty0) {
					if (ii == Np_n - 1u) {
						if (crystal_size_z == 0.f) {
							orth_distance_cosem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, MethodList,
								d_alku, local_sino, d_Ny, 1u, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, true);
						}
						else {
							orth_distance_cosem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
								MethodList, d_alku, local_sino, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
								&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, true);
						}
					}
					else {
						tempi += iu;
						tx0 += txu;
					}
				}
				else {
					if (crystal_size_z == 0.f) {
						orth_distance_cosem(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, MethodList,
							d_alku, local_sino, d_Ny, 1u, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
							&minimi, MBSREM_prepass, d_Summ, &axACOSEM, true);
					}
					else {
						orth_distance_cosem_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, temp_ijk,
							MethodList, d_alku, local_sino, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, &axCOSEM, d_COSEM, d_h, d_E, idx, d_co, d_aco,
							&minimi, MBSREM_prepass, d_Summ, &axACOSEM, dec, true);
					}
					tempj += ju;
					if (ju > 0)
						temp_ijk += d_Nx;
					else
						temp_ijk -= d_Nx;
					ty0 += tyu;
				}
			}
			if ((MethodList.MRAMLA == 1 || MethodList.MBSREM == 1) && MBSREM_prepass == 1)
				d_Amin[idx] = minimi;
			if ((MethodList.ACOSEM == 1 || MethodList.OSLCOSEM == 1) && d_alku > 0) {
				if (d_randoms == 1u)
					axACOSEM += d_sc_ra[idx];
				d_ACOSEM_lhs[idx] = axACOSEM;
			}
		}
	}
}