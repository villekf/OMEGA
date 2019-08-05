/**************************************************************************
* A matrix free improved Siddon's combined with all the reconstruction
* functions available in OMEGA. This function calculates Summ = sum(A,1) 
* (sum of every row) and rhs = A*(y./(A'*x)), where A is the system matrix, 
* y the measurements and x the estimate/image.
*
* INPUTS:
* d_rhs_X = buffer for rhs, d_Summ = buffer for Summ, d_lor = row
* indices, d_Nx = image size in x-dimension, d_dz = distance between
* adjecent voxels in z-dimension, d_bz = distance from the pixel space
* to origin (z-dimension), d_bzb = part in parenthesis of equation
* (9) in [1] precalculated when k = Nz, d_maxxx = maximum distance of the
* pixel space from origin in x-dimension, d_zmax = maximum value of d_zdet,
* d_NSlices = the number of image slices, d_x = detector x-coordinates,
* d_size_x = the number of detector elements, d_row = how many voxels have
* been traversed so far, d_xyindex = for sinogram format they determine the
* detector indices corresponding to each sinogram bin (unused with raw data)
* d_TotSinos = Total number of sinograms, d_attenuation_correction = if
* attenuation is included this is 1 otherwise 0, d_atten = attenuation
* data (images), d_L = detector numbers for raw data (unused for sinogram
* format), d_det_per_ring = number of detectors per ring, d_pseudos =
* location of pseudo rings, pRows = number of pseudo rings, d_raw = if
* 1 then raw list-mode data is used otherwise sinogram data, d_rekot = 
* vector containing 1 if the reconstruction method is included (e.g. if
* rekot[1] = 1, then OSEM is calculated) and 0 if not, d_X = buffer for 
* algorithm X, d_X_OSEM = buffer for the OSL estimate of prior X, 
* d_X_BSREM = same as before, but for BSREM
*
* OUTPUTS:
* d_rsh_X = rhs values for algorithm/prior X, d_Summ = sum values, 
* d_X = estimate of algorithm X, d_X_OSEM = OSL estimate of prior X, 
* d_X_BSREM = BSREM estimate of prior X
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
#include "opencl_functions.h"
#define TYPE 0

// Matrix free Improved Siddon's algorithm
__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(64, 1, 1)))
void f_b_project_orth(const uchar fp, const uint d_N, const uint d_Nx, const uint d_Ny, const uint d_Nz, const float d_dz, const float d_dx,
	const float d_dy, const float d_bz, const float d_bx, const float d_by, const float d_bzb, const float d_maxxx, const float d_maxyy,
	const float d_zmax, const float d_NSlices, const uint d_size_x, const ushort d_TotSinos, const uint d_attenuation_correction, const uint d_normalization, const uint d_randoms, const uint d_det_per_ring,
	const uchar d_raw, const uint d_pRows, const uchar no_norm, const uint d_Nxy, const float tube_width_xy, const float crystal_size_z, const int dec,
	const __global float* d_atten, const __global float* d_norm, __global float* d_Summ, const __global ushort* d_lor, __constant uint* d_pseudos, const __global float* d_x,
	const __global float* d_y, const __global float* d_zdet, __constant float* x_center, __constant float* y_center, __constant float* z_center, const __global uint* d_xyindex,
	const __global ushort* d_zindex, const __global ushort* d_L, const __global float* d_sc_ra, const __global float* d_rhs, __global float* d_output, const ulong m_size) {
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
		get_detector_coordinates(d_xyindex, d_zindex, d_size_x, idx, d_TotSinos, &xs, &xd, &ys, &yd, &zs, &zd, d_x, d_y, d_zdet);
	}
	float local_sc_ra = 0.f;
	float local_norm = 0.f;
	if (d_normalization == 1u)
		local_norm = d_norm[idx];
	if (d_randoms == 1u)
		local_sc_ra = d_sc_ra[idx];
	// Calculate the x, y and z distances of the detector pair
	const float y_diff = (yd - ys);
	const float x_diff = (xd - xs);
	const float z_diff = (zd - zs);
	if ((y_diff == 0.f && x_diff == 0.f && z_diff == 0.f) || (y_diff == 0.f && x_diff == 0.f))
		return;
	uint Np = 0u;
	uint Np_n = 0u;
	uint xyz = 0u;
	float d_rhs_local = 0.f, d_output_local = 0.f, LL;
	if (!fp)
		d_rhs_local = d_rhs[idx];
	float jelppi = 0.f;
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
		const uint tempk = convert_uint((zs / d_zmax) * (d_NSlices - 1.f));
		// If the LOR is perpendicular in the y-direction (Siddon cannot be used)
		if (fabs(y_diff) < 1e-6f) {
			if (yd <= d_maxyy && yd >= d_by) {
				float temp = 0.f;
				if (crystal_size_z == 0.f) {
					orth_distance_perpendicular_bpfp(-x_diff, y_center, kerroin, length_, &temp, d_attenuation_correction, d_normalization, d_by, yd, d_dy, d_Ny,
						d_Nx, tempk, d_atten, local_norm, d_Ny, 1u, d_rhs, fp, &d_output_local, d_output, d_Summ, no_norm, d_rhs_local, true);
					if (fp)
						d_output[idx] = d_output_local * temp + local_sc_ra;
					else {
						orth_distance_perpendicular_bpfp(-x_diff, y_center, kerroin, length_, &temp, d_attenuation_correction, d_normalization, d_by, yd, d_dy, d_Ny,
							d_Nx, tempk, d_atten, local_norm, d_Ny, 1u, d_rhs, fp, &d_output_local, d_output, d_Summ, no_norm, d_rhs_local, false);
					}
				}
				else {
					//const float temppi = xs;
					//xs = ys;
					//ys = temppi;
					orth_distance_perpendicular_bpfp_3D(y_center, x_center[0], z_center, &temp, d_attenuation_correction, d_normalization,
						d_by, yd, d_dy, d_Ny, d_Nx, tempk, d_atten, local_norm, d_Ny, 1u, d_rhs, fp, &d_output_local, ys, xs, zs, y_diff, x_diff, z_diff, kerroin, d_Nxy, d_Nz,
						d_output, d_Summ, no_norm, d_rhs_local, true);
					if (fp)
						d_output[idx] = d_output_local * temp + local_sc_ra;
					else {
						orth_distance_perpendicular_bpfp_3D(y_center, x_center[0], z_center, &temp, d_attenuation_correction, d_normalization,
							d_by, yd, d_dy, d_Ny, d_Nx, tempk, d_atten, local_norm, d_Ny, 1u, d_rhs, fp, &d_output_local, ys, xs, zs, y_diff, x_diff, z_diff, kerroin, d_Nxy, d_Nz,
							d_output, d_Summ, no_norm, d_rhs_local, false);
					}
				}
			}
		}
		else if (fabs(x_diff) < 1e-6f) {
			if (xd <= d_maxxx && xd >= d_bx) {
				float temp = 0.f;
				if (crystal_size_z == 0.f) {
					orth_distance_perpendicular_bpfp(y_diff, x_center, kerroin, length_, &temp, d_attenuation_correction, d_normalization, d_bx, xd, d_dx, d_Nx,
						d_Ny, tempk, d_atten, local_norm, 1u, d_Nx, d_rhs, fp, &d_output_local, d_output, d_Summ, no_norm, d_rhs_local, true);
					if (fp)
						d_output[idx] = d_output_local * temp + local_sc_ra;
					else {
						orth_distance_perpendicular_bpfp(y_diff, x_center, kerroin, length_, &temp, d_attenuation_correction, d_normalization, d_bx, xd, d_dx, d_Nx,
							d_Ny, tempk, d_atten, local_norm, 1u, d_Nx, d_rhs, fp, &d_output_local, d_output, d_Summ, no_norm, d_rhs_local, false);
					}
				}
				else {
					orth_distance_perpendicular_bpfp_3D(x_center, y_center[0], z_center, &temp, d_attenuation_correction, d_normalization,
						d_bx, xd, d_dx, d_Nx, d_Ny, tempk, d_atten, local_norm, 1u, d_Nx, d_rhs, fp, &d_output_local, xs, ys, zs, x_diff, y_diff, z_diff,
						kerroin, d_Nxy, d_Nz, d_output, d_Summ, no_norm, d_rhs_local, true);
					if (fp)
						d_output[idx] = d_output_local * temp + local_sc_ra;
					else {
						orth_distance_perpendicular_bpfp_3D(x_center, y_center[0], z_center, &temp, d_attenuation_correction, d_normalization,
							d_bx, xd, d_dx, d_Nx, d_Ny, tempk, d_atten, local_norm, 1u, d_Nx, d_rhs, fp, &d_output_local, xs, ys, zs, x_diff, y_diff, z_diff,
							kerroin, d_Nxy, d_Nz, d_output, d_Summ, no_norm, d_rhs_local, false);
					}
				}
			}
		}
		else {
			int tempi = 0, tempj = 0, iu = 0, ju = 0;
			float txu = 0.f, tyu = 0.f, tc = 0.f, tx0 = 0.f, ty0 = 0.f;
			const bool skip = siddon_pre_loop_2D(d_bx, d_by, x_diff, y_diff, d_maxxx, d_maxyy, d_dx, d_dy, d_Nx, d_Ny, &tempi, &tempj, &txu, &tyu, &Np, TYPE,
				ys, xs, yd, xd, &tc, &iu, &ju, &tx0, &ty0);
			if (skip)
				return;
			float temp = 0.f;
			float tx0_a = tx0, ty0_a = ty0;
			int tempi_a = tempi, tempj_a = tempj;
			uint temp_ijk;
			if (crystal_size_z == 0.f)
				temp_ijk = tempk * d_Nxy + convert_uint_sat(tempj) * d_Nx;
			else
				temp_ijk = convert_uint_sat(tempj) * d_Nx;
			if (d_attenuation_correction)
				LL = native_sqrt(x_diff * x_diff + y_diff * y_diff);
			for (uint ii = 0u; ii < Np; ii++) {
				if (tx0 < ty0) {
					if (d_attenuation_correction)
						compute_attenuation(&tc, &jelppi, LL, tx0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
					if (ii == Np - 1u) {
						if (crystal_size_z == 0.f) {
							orth_distance_bpfp(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk, 
								d_Ny, 1u, tempj, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, 
								true);
						}
						else {
							orth_distance_bpfp_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp, 
								temp_ijk, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
						}
					}
					else {
						tempi += iu;
						tx0 += txu;
					}
					xyz = 1u;
				}
				else {
					if (d_attenuation_correction)
						compute_attenuation(&tc, &jelppi, LL, ty0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
					if (crystal_size_z == 0.f) {
						orth_distance_bpfp(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk,
							d_Ny, 1u, tempj, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
					}
					else {
						orth_distance_bpfp_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp,
							temp_ijk, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
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
							orth_distance_bpfp(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk,
								d_Ny, 1u, tempj, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
						}
						else {
							orth_distance_bpfp_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp,
								temp_ijk, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
						}
					}
					break;
				}
			}
			temp = 1.f / temp;
			if (d_attenuation_correction)
				temp *= native_exp(jelppi);
			if (d_normalization == 1u)
				temp *= local_norm;
			if (fp)
				d_output[idx] = d_output_local * temp + local_sc_ra;
			else {
				tx0 = tx0_a, ty0 = ty0_a;
				tempi = tempi_a, tempj = tempj_a;
				if (crystal_size_z == 0.f)
					temp_ijk = tempk * d_Nxy + convert_uint_sat(tempj) * d_Nx;
				else
					temp_ijk = convert_uint_sat(tempj) * d_Nx;
				for (uint ii = 0u; ii < Np_n; ii++) {
					if (tx0 < ty0) {
						if (ii == Np_n - 1u) {
							if (crystal_size_z == 0.f) {
								orth_distance_bpfp(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk,
									d_Ny, 1u, tempj, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local,
									false);
							}
							else {
								orth_distance_bpfp_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp,
									temp_ijk, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, false);
							}
						}
						else {
							tempi += iu;
							tx0 += txu;
						}
					}
					else {
						if (crystal_size_z == 0.f) {
							orth_distance_bpfp(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk,
								d_Ny, 1u, tempj, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local,
								false);
						}
						else {
							orth_distance_bpfp_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp,
								temp_ijk, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, false);
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
		if (fabs(y_diff) < 1e-6f) {
			if (yd <= d_maxyy && yd >= d_by) {
				int tempi = 0, tempk = 0, tempj = 0, iu = 0, ku = 0;
				float txu = 0.f, tzu = 0.f, tx0 = 0.f, tz0 = 0.f, tc = 0.f;
				const bool skip = siddon_pre_loop_2D(d_bx, d_bz, x_diff, z_diff, d_maxxx, d_bzb, d_dx, d_dz, d_Nx, d_Nz, &tempi, &tempk, &txu, &tzu, &Np, TYPE,
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
					//const float temp_x = xs;
					//xs = ys;
					//ys = temp_x;
				}
				if (d_attenuation_correction)
					LL = native_sqrt(x_diff * x_diff + z_diff * z_diff);
				for (uint ii = 0u; ii < Np; ii++) {
					if (tx0 < tz0) {
						if (d_attenuation_correction)
							compute_attenuation(&tc, &jelppi, LL, tx0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
						if (crystal_size_z == 0.f) {
							orth_distance_bpfp(tempj, d_Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, &temp, temp_ijk, 
								d_Nx, d_Nx, tempj, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
						}
						else {
							orth_distance_bpfp_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, &temp, temp_ijk,
								d_Nx, d_Nx, tempk, d_Nxy, ys, xs, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
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
						if (d_attenuation_correction)
							compute_attenuation(&tc, &jelppi, LL, tz0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
						if (crystal_size_z == 0.f) {
							orth_distance_bpfp(tempj, d_Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, &temp, temp_ijk, 
								d_Nx, d_Nx, tempj, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
							if (ku > 0)
								temp_ijk += d_Nxy;
							else
								temp_ijk -= d_Nxy;
						}
						else if (ii == Np - 1u) {
							orth_distance_bpfp_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, &temp, temp_ijk,
								d_Nx, d_Nx, tempk, d_Nxy, ys, xs, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
						}
						tempk += ku;
						tz0 += tzu;
						xyz = 3u;
					}
					Np_n++;
					if (tempk < 0 || tempi < 0 || tempi >= d_Nx || tempk >= d_Nz) {
						if (crystal_size_z != 0.f && xyz == 3u) {
							orth_distance_bpfp_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, &temp, temp_ijk,
								d_Nx, d_Nx, tempk, d_Nxy, ys, xs, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
						}
						break;
					}
				}
				temp = 1.f / temp;
				if (d_attenuation_correction)
					temp *= native_exp(jelppi);
				if (d_normalization == 1u)
					temp *= local_norm;
				if (fp)
					d_output[idx] = d_output_local * temp + local_sc_ra;
				else {
					tx0 = tx0_a, tz0 = tz0_a;
					tempi = tempi_a, tempk = tempk_a;
					if (crystal_size_z == 0.f)
						temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempi);
					else
						temp_ijk = convert_uint_sat(tempi);
					for (uint ii = 0u; ii < Np_n; ii++) {
						if (tx0 < tz0) {
							if (crystal_size_z == 0.f) {
								orth_distance_bpfp(tempj, d_Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, &temp, temp_ijk, 
									d_Nx, d_Nx, tempj, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, false);
							}
							else {
								orth_distance_bpfp_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, &temp, temp_ijk,
									d_Nx, d_Nx, tempk, d_Nxy, ys, xs, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, false);
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
								orth_distance_bpfp(tempj, d_Ny, -x_diff, -y_diff, x_center[tempi], y_center, kerroin, length_, &temp, temp_ijk, 
									d_Nx, d_Nx, tempj, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, false);
								if (ku > 0)
									temp_ijk += d_Nxy;
								else
									temp_ijk -= d_Nxy;
							}
							else if (ii == Np_n - 1u) {
								orth_distance_bpfp_3D(tempj, d_Ny, d_Nz, x_diff, y_diff, z_diff, x_center[tempi], y_center, z_center, kerroin, &temp, temp_ijk,
									d_Nx, d_Nx, tempk, d_Nxy, ys, xs, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, false);
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
				const bool skip = siddon_pre_loop_2D(d_by, d_bz, y_diff, z_diff, d_maxyy, d_bzb, d_dy, d_dz, d_Ny, d_Nz, &tempj, &tempk, &tyu, &tzu, &Np, TYPE,
					zs, ys, zd, yd, &tc, &ju, &ku, &ty0, &tz0);
				if (skip)
					return;
				float temp = 0.f;
				tempi = perpendicular_start(d_bx, xd, d_dx, d_Nx);
				float tz0_a = tz0, ty0_a = ty0;
				int tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
				uint temp_ijk;
				if (crystal_size_z == 0.f)
					temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempj) * d_Nx;
				else
					temp_ijk = convert_uint_sat(tempj) * d_Nx;
				if (d_attenuation_correction)
					LL = native_sqrt(y_diff * y_diff + z_diff * z_diff);
				for (uint ii = 0u; ii < Np; ii++) {
					if (tz0 < ty0) {
						if (d_attenuation_correction)
							compute_attenuation(&tc, &jelppi, LL, tz0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
						if (crystal_size_z == 0.f) {
							orth_distance_bpfp(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk,
								d_Ny, 1u, tempj, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local,
								true);
							if (ku > 0)
								temp_ijk += d_Nxy;
							else
								temp_ijk -= d_Nxy;
						}
						else if (ii == Np - 1u) {
							orth_distance_bpfp_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp,
								temp_ijk, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
						}
						tempk += ku;
						tz0 += tzu;
						xyz = 3u;
					}
					else {
						if (d_attenuation_correction)
							compute_attenuation(&tc, &jelppi, LL, ty0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
						if (crystal_size_z == 0.f) {
							orth_distance_bpfp(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk,
								d_Ny, 1u, tempj, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local,
								true);
						}
						else {
							orth_distance_bpfp_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp,
								temp_ijk, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
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
							orth_distance_bpfp_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp,
								temp_ijk, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
						}
						break;
					}
				}
				temp = 1.f / temp;
				if (d_attenuation_correction)
					temp *= native_exp(jelppi);
				if (d_normalization == 1u)
					temp *= local_norm;
				if (fp)
					d_output[idx] = d_output_local * temp + local_sc_ra;
				else {
					tz0 = tz0_a, ty0 = ty0_a;
					tempi = tempi_a, tempj = tempj_a, tempk = tempk_a;
					if (crystal_size_z == 0.f)
						temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempj) * d_Nx;
					else
						temp_ijk = convert_uint_sat(tempj) * d_Nx;
					for (uint ii = 0u; ii < Np_n; ii++) {
						if (tz0 < ty0) {
							if (crystal_size_z == 0.f) {
								orth_distance_bpfp(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk,
									d_Ny, 1u, tempj, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, false);
								if (ku > 0)
									temp_ijk += d_Nxy;
								else
									temp_ijk -= d_Nxy;
							}
							else if (ii == Np_n - 1u) {
								orth_distance_bpfp_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp,
									temp_ijk, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, false);
							}
							tempk += ku;
							tz0 += tzu;
						}
						else {
							if (crystal_size_z == 0.f) {
								orth_distance_bpfp(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk,
									d_Ny, 1u, tempj, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, false);
							}
							else {
								orth_distance_bpfp_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp,
									temp_ijk, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, false);
							}
							tempj += (ju);
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
			int tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 1;
			float txu = 0.f, tyu = 0.f, tzu = 0.f, tc = 0.f, tx0 = 0.f, ty0 = 0.f, tz0 = 0.f;
			const bool skip = siddon_pre_loop_3D(d_bx, d_by, d_bz, x_diff, y_diff, z_diff, d_maxxx, d_maxyy, d_bzb, d_dx, d_dy, d_dz, d_Nx, d_Ny, d_Nz, &tempi, &tempj, &tempk, &tyu, &txu, &tzu,
				&Np, TYPE, ys, xs, yd, xd, zs, zd, &tc, &iu, &ju, &ku, &tx0, &ty0, &tz0);
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
			if (d_attenuation_correction)
				LL = native_sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
			for (uint ii = 0u; ii < Np; ii++) {
				if (tz0 < ty0 && tz0 < tx0) {
					if (d_attenuation_correction)
						compute_attenuation(&tc, &jelppi, LL, tz0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
					if (crystal_size_z == 0.f) {
						orth_distance_bpfp(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk,
							d_Ny, 1u, tempj, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
						if (ku > 0)
							temp_ijk += d_Nxy;
						else
							temp_ijk -= d_Nxy;
					}
					else if (ii == Np - 1u) {
						orth_distance_bpfp_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp,
							temp_ijk, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
					}
					tempk += ku;
					tz0 += tzu;
					xyz = 3u;
				}
				else if (ty0 < tx0) {
					if (d_attenuation_correction)
						compute_attenuation(&tc, &jelppi, LL, ty0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
					if (crystal_size_z == 0.f) {
						orth_distance_bpfp(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk,
							d_Ny, 1u, tempj, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
					}
					else {
						orth_distance_bpfp_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp,
							temp_ijk, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
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
					if (d_attenuation_correction)
						compute_attenuation(&tc, &jelppi, LL, tx0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
					if (ii == Np - 1u) {
						if (crystal_size_z == 0.f) {
							orth_distance_bpfp(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk,
								d_Ny, 1u, tempj, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
						}
						else {
							orth_distance_bpfp_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp,
								temp_ijk, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
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
							orth_distance_bpfp(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk,
								d_Ny, 1u, tempj, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
						}
						else {
							orth_distance_bpfp_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp,
								temp_ijk, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, true);
						}
					}
					break;
				}
			}
			temp = 1.f / temp;
			if (d_attenuation_correction)
				temp *= native_exp(jelppi);
			if (d_normalization == 1u)
				temp *= local_norm;
			if (fp)
				d_output[idx] = d_output_local * temp + local_sc_ra;
			else {
				tx0 = tx0_a, ty0 = ty0_a, tz0 = tz0_a;
				tempi = tempi_a, tempj = tempj_a, tempk = tempk_a;
				if (crystal_size_z == 0.f)
					temp_ijk = convert_uint_sat(tempk) * d_Nxy + convert_uint_sat(tempj) * d_Nx;
				else
					temp_ijk = convert_uint_sat(tempj) * d_Nx;
				for (uint ii = 0u; ii < Np_n; ii++) {
					if (tz0 < ty0 && tz0 < tx0) {
						if (crystal_size_z == 0.f) {
							orth_distance_bpfp(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk,
								d_Ny, 1u, tempj, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, false);
							if (ku > 0)
								temp_ijk += d_Nxy;
							else
								temp_ijk -= d_Nxy;
						}
						else if (ii == Np_n - 1u) {
							orth_distance_bpfp_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp,
								temp_ijk, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, false);
						}
						tempk += ku;
						tz0 += tzu;
					}
					else if (ty0 < tx0) {
						if (crystal_size_z == 0.f) {
							orth_distance_bpfp(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk,
								d_Ny, 1u, tempj, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, false);
						}
						else {
							orth_distance_bpfp_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp,
								temp_ijk, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, false);
						}
						tempj += (ju);
						if (ju > 0)
							temp_ijk += d_Nx;
						else
							temp_ijk -= d_Nx;
						ty0 += tyu;
					}
					else {
						if (ii == Np_n - 1u) {
							if (crystal_size_z == 0.f) {
								orth_distance_bpfp(tempi, d_Nx, y_diff, x_diff, y_center[tempj], x_center, kerroin, length_, &temp, temp_ijk,
									d_Ny, 1u, tempj, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, false);
							}
							else {
								orth_distance_bpfp_3D(tempi, d_Nx, d_Nz, y_diff, x_diff, z_diff, y_center[tempj], x_center, z_center, kerroin, &temp,
									temp_ijk, d_Ny, 1u, tempk, d_Nxy, xs, ys, zs, dec, fp, &d_output_local, d_rhs, d_Summ, d_output, no_norm, d_rhs_local, false);
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

__kernel void summa(const __global float* d_Summ_device, __global float* d_Summ_local, const __global float* d_rhs_device, __global float* d_rhs_local, 
	const ulong size_rhs, const uint im_dim, const uchar no_norm) {

	uint gid = get_global_id(0);

	for (uint i = gid; i < size_rhs; i += get_global_size(0)) {
		d_rhs_local[i] += d_rhs_device[i];
		if (!no_norm && i < im_dim)
			d_Summ_local[i] += d_Summ_device[i];
	}
}