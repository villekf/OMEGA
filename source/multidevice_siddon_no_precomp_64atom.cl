/**************************************************************************
* A matrix free improved Siddon's for OSEM or MLEM.
* This function calculates Summ = sum(A,1) (sum of every row) and
* rhs = A*(y./(A'*x)), where A is the system matrix, y the measurements
* and x the estimate/image.
*
* Used by implementation 3.
*
* This version goes through all the LORs and determines on-the-fly if they
* intersect with the voxel space. Uses integer 64-bit atomics.
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
* dc_z = Distance between the rays (z-direction),
* n_rays = number of rays used,
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
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
#include "opencl_functions_a64.h"
#define TYPE 0
#define TH 100000000000.f

// Matrix free Improved Siddon's algorithm
__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(64, 1, 1)))
void siddon_multi(const float d_epps, const uint d_N, const uint d_Nx, const uint d_Ny, const uint d_Nz, const float d_dz, const float d_dx,
	const float d_dy, const float d_bz, const float d_bx, const float d_by, const float d_bzb, const float d_maxxx,	const float d_maxyy, 
	const float d_zmax, const float d_NSlices, const uint d_size_x, const ushort d_TotSinos, const uint d_attenuation_correction, const uint d_normalization, const uint d_randoms, const uint d_det_per_ring,
	const uchar d_raw, const uint d_pRows, const uint d_Nxy, const float dc_z, const ushort n_rays,
	const __global float* d_atten, const __global float* d_norm, __global ulong* d_Summ, const __global ushort* d_lor, __constant uint* d_pseudos, const __global float* d_x,
	const __global float* d_y, const __global float* d_zdet, 
	const __global uint* d_xyindex, const __global ushort* d_zindex, const __global ushort* d_L, const __global float* d_Sino, const __global float* d_sc_ra, const __global float* d_OSEM,
	__global ulong* d_rhs_OSEM, const uchar no_norm, const ulong m_size) {
	// Get the current global index
	uint idx = get_global_id(0);
	if (idx >= m_size)
		return;
	const float local_sino = (d_Sino[idx]);
	if (no_norm == 1u && local_sino == 0.f)
		return;
	float jelppi = 0.f;
	float axOSEM = 0.f;
	float temp = 0.f;
	int tempi_a[5];
	int tempj_a[5];
	int tempk_a[5];
	int iu_a[5];
	int ju_a[5];
	int ku_a[5];
	float tx0_a[5];
	float ty0_a[5];
	float tz0_a[5];
	float tc_a[5];
	float txu_a[5];
	float tyu_a[5];
	float tzu_a[5];
	float LL[5];
	uint Np_n[5];
	bool pass[5];
	// Load the next detector index
	// raw list-mode data
	for (ushort lor = 0u; lor < n_rays; lor++) {
		float xs, xd, ys, yd, zs, zd;
		if (d_raw) {
			get_detector_coordinates_raw_multiray(d_x, d_y, d_zdet, d_L, d_det_per_ring, idx, d_pseudos, d_pRows, &xs, &xd, &ys, &yd, &zs, &zd, lor + 1u, dc_z);
		}
		// Sinogram data
		else {
			get_detector_coordinates_multiray(d_xyindex, d_zindex, d_size_x, idx, d_TotSinos, &xs, &xd, &ys, &yd, &zs, &zd, d_x, d_y, d_zdet, lor + 1u, dc_z);
		}
		//float local_sc_ra = 0.f;
		//float local_norm = 0.f;
		//if (d_normalization == 1u)
		//	local_norm = d_norm[idx];
		//if (d_randoms == 1u)
		//	local_sc_ra = d_sc_ra[idx];
		// Calculate the x, y and z distances of the detector pair
		const float y_diff = (yd - ys);
		const float x_diff = (xd - xs);
		const float z_diff = (zd - zs);
		if ((y_diff == 0.f && x_diff == 0.f && z_diff == 0.f) || (y_diff == 0.f && x_diff == 0.f)) {
			pass[lor] = false;
			continue;
		}
		uint Np = 0u;
		Np_n[lor] = 0u;
		// If the measurement is on a same ring
		if (fabs(z_diff) < 1e-6f) {
			// Z-coordinate (ring)
			const uint z_loop = convert_uint((zs / d_zmax) * (d_NSlices - 1.f));
			// If the LOR is perpendicular in the y-direction (Siddon cannot be used)
			if (fabs(y_diff) < 1e-6f) {
				if (yd <= d_maxyy && yd >= d_by) {
					uint apu = 0u;

					const float element = perpendicular_elements_multiray(d_by, d_dy, d_Ny, yd, d_dx, d_Nx, d_atten, &apu, d_attenuation_correction, z_loop, d_Ny, 1u, &jelppi);
					temp += element;
					tempk_a[lor] = convert_int(apu);
					if (local_sino > 0.f) {
						for (uint k = 0u; k < d_Nx; k++) {
							axOSEM += (d_dx * d_OSEM[apu + k]);
						}
					}
					tx0_a[lor] = 1e7f, ty0_a[lor] = 1e9f;
					pass[lor] = true;
				}
				else
					pass[lor] = false;
			}
			else if (fabs(x_diff) < 1e-6f) {
				if (xd <= d_maxxx && xd >= d_bx) {
					uint apu = 0;

					const float element = perpendicular_elements_multiray(d_bx, d_dx, d_Nx, xd, d_dy, d_Ny, d_atten, &apu, d_attenuation_correction, z_loop, 1u, d_Nx, &jelppi);

					temp += element;
					tempk_a[lor] = convert_int(apu);
					if (local_sino > 0.f) {
						for (uint k = 0u; k < d_Ny; k++) {
							axOSEM += (d_dy * d_OSEM[apu + k * d_Nx]);
						}
					}
					tx0_a[lor] = 1e9f, ty0_a[lor] = 1e7f;
					pass[lor] = true;
				}
				else
					pass[lor] = false;
			}
			else {
				int tempi = 0, tempj = 0, iu = 0, ju = 0;
				float txu = 0.f, tyu = 0.f, tc = 0.f, tx0 = 0.f, ty0 = 0.f;
				const bool skip = siddon_pre_loop_2D(d_bx, d_by, x_diff, y_diff, d_maxxx, d_maxyy, d_dx, d_dy, d_Nx, d_Ny, &tempi, &tempj, &txu, &tyu, &Np, TYPE,
					ys, xs, yd, xd, &tc, &iu, &ju, &tx0, &ty0);
				if (skip || tempi < 0 || tempj < 0 || tempi >= convert_int(d_Nx) || tempj >= convert_int(d_Ny)) {
					pass[lor] = false;
					continue;
				}
				LL[lor] = hypot(x_diff, y_diff); //native_sqrt(x_diff * x_diff + y_diff * y_diff);
				tx0_a[lor] = tx0, ty0_a[lor] = ty0, tz0_a[lor] = 1e8f, tc_a[lor] = tc;
				txu_a[lor] = txu, tyu_a[lor] = tyu, tzu_a[lor] = 1e8f;
				tempi_a[lor] = tempi, tempj_a[lor] = tempj, tempk_a[lor] = convert_int(z_loop);
				iu_a[lor] = iu, ju_a[lor] = ju, ku_a[lor] = 0;
				float local_ele;
				for (uint ii = 0u; ii < Np; ii++) {
					const uint local_ind = compute_ind(tempj, tempi, convert_int(z_loop), d_Nx, d_Ny, d_N, d_Nx, d_Nxy);
					if (tx0 < ty0) {
						local_ele = compute_element(&tx0, &tc, LL[lor], txu, iu, &tempi, &temp);
					}
					else {
						local_ele = compute_element(&ty0, &tc, LL[lor], tyu, ju, &tempj, &temp);
					}
					if (d_attenuation_correction) {
						jelppi += (local_ele * -d_atten[local_ind]);
					}
					if (local_sino > 0.f) {
						denominator_multi(local_ele, &axOSEM, &d_OSEM[local_ind]);
					}
					Np_n[lor]++;
					if (tempj < 0 || tempi < 0 || tempi >= d_Nx || tempj >= d_Ny)
						break;
				}
				pass[lor] = true;
			}
		}
		else {
			if (fabs(y_diff) < 1e-6f) {
				if (yd <= d_maxyy && yd >= d_by) {
					int tempi = 0, tempk = 0, tempj = 0, iu = 0, ku = 0;
					float txu = 0.f, tzu = 0.f, tc = 0.f, tx0 = 0.f, tz0 = 0.f;
					const bool skip = siddon_pre_loop_2D(d_bx, d_bz, x_diff, z_diff, d_maxxx, d_bzb, d_dx, d_dz, d_Nx, d_Nz, &tempi, &tempk, &txu, &tzu, &Np, TYPE,
						zs, xs, zd, xd, &tc, &iu, &ku, &tx0, &tz0);
					if (skip || tempi < 0 || tempk < 0 || tempi >= convert_int(d_Nx) || tempk >= convert_int(d_Nz)) {
						pass[lor] = false;
						continue;
					}
					LL[lor] = hypot(x_diff, z_diff); //native_sqrt((x_diff * x_diff + z_diff * z_diff));
					tempj = perpendicular_start(d_by, yd, d_dy, d_Ny);
					//tempj *= convert_int(d_Nx);
					tx0_a[lor] = tx0, ty0_a[lor] = 1e8f, tz0_a[lor] = tz0, tc_a[lor] = tc;
					txu_a[lor] = txu, tyu_a[lor] = 1e8f, tzu_a[lor] = tzu;
					tempi_a[lor] = tempi, tempj_a[lor] = tempj, tempk_a[lor] = tempk;
					iu_a[lor] = iu, ju_a[lor] = 0, ku_a[lor] = ku;
					float local_ele;
					for (uint ii = 0u; ii < Np; ii++) {
						const uint local_ind = compute_ind(tempj, tempi, tempk, d_Nx, d_Ny, d_N, d_Nx, d_Nxy);
						if (tx0 < tz0) {
							local_ele = compute_element(&tx0, &tc, LL[lor], txu, iu, &tempi, &temp);
						}
						else {
							local_ele = compute_element(&tz0, &tc, LL[lor], tzu, ku, &tempk, &temp);
						}
						if (d_attenuation_correction) {
							jelppi += (local_ele * -d_atten[local_ind]);
						}
						if (local_sino > 0.f) {
							denominator_multi(local_ele, &axOSEM, &d_OSEM[local_ind]);
						}
						Np_n[lor]++;
						if (tempk < 0 || tempi < 0 || tempi >= d_Nx || tempk >= d_Nz)
							break;
					}
					pass[lor] = true;
				}
				else
					pass[lor] = false;
			}
			else if (fabs(x_diff) < 1e-6f) {
				if (xd <= d_maxxx && xd >= d_bx) {
					int tempi = 0, tempk = 0, tempj = 0, ju = 0, ku = 0;
					float tyu = 0.f, tzu = 0.f, tc = 0.f, ty0 = 0.f, tz0 = 0.f;
					const bool skip = siddon_pre_loop_2D(d_by, d_bz, y_diff, z_diff, d_maxyy, d_bzb, d_dy, d_dz, d_Ny, d_Nz, &tempj, &tempk, &tyu, &tzu, &Np, TYPE,
						zs, ys, zd, yd, &tc, &ju, &ku, &ty0, &tz0);
					if (skip || tempk < 0 || tempj < 0 || tempk >= convert_int(d_Nz) || tempj >= convert_int(d_Ny)) {
						pass[lor] = false;
						continue;
					}
					LL[lor] = hypot(y_diff, z_diff); //native_sqrt((y_diff * y_diff + z_diff * z_diff));
					tempi = perpendicular_start(d_bx, xd, d_dx, d_Nx);
					tx0_a[lor] = 1e8f, ty0_a[lor] = ty0, tz0_a[lor] = tz0, tc_a[lor] = tc;
					txu_a[lor] = 1e8f, tyu_a[lor] = tyu, tzu_a[lor] = tzu;
					tempi_a[lor] = tempi, tempj_a[lor] = tempj, tempk_a[lor] = tempk;
					iu_a[lor] = 0, ju_a[lor] = ju, ku_a[lor] = ku;
					float local_ele;
					for (uint ii = 0u; ii < Np; ii++) {
						const uint local_ind = compute_ind(tempj, tempi, tempk, d_Nx, d_Ny, d_N, d_Nx, d_Nxy);
						if (tz0 < ty0) {
							local_ele = compute_element(&tz0, &tc, LL[lor], tzu, ku, &tempk, &temp);
						}
						else {
							local_ele = compute_element(&ty0, &tc, LL[lor], tyu, ju, &tempj, &temp);
						}
						if (d_attenuation_correction) {
							jelppi += (local_ele * -d_atten[local_ind]);
						}
						if (local_sino > 0.f) {
							denominator_multi(local_ele, &axOSEM, &d_OSEM[local_ind]);
						}
						Np_n[lor]++;
						if (tempj < 0 || tempk < 0 || tempk >= d_Nz || tempj >= d_Ny)
							break;
					}
					pass[lor] = true;
				}
				else
					pass[lor] = false;
			}
			else {
				int tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;
				float txu = 0.f, tyu = 0.f, tzu = 0.f, tc = 0.f, tx0 = 0.f, ty0 = 0.f, tz0 = 0.f;
				const bool skip = siddon_pre_loop_3D(d_bx, d_by, d_bz, x_diff, y_diff, z_diff, d_maxxx, d_maxyy, d_bzb, d_dx, d_dy, d_dz, d_Nx, d_Ny, d_Nz, &tempi, &tempj, &tempk, &tyu, &txu, &tzu,
					&Np, TYPE, ys, xs, yd, xd, zs, zd, &tc, &iu, &ju, &ku, &tx0, &ty0, &tz0);
				if (skip || tempi < 0 || tempj < 0 || tempk < 0 || tempi >= convert_int(d_Nx) || tempj >= convert_int(d_Ny)
					|| tempk >= convert_int(d_Nz)) {
					pass[lor] = false;
					continue;
				}
				LL[lor] = native_sqrt(x_diff * x_diff + z_diff * z_diff + y_diff * y_diff);
				tx0_a[lor] = tx0, ty0_a[lor] = ty0, tz0_a[lor] = tz0, tc_a[lor] = tc;
				txu_a[lor] = txu, tyu_a[lor] = tyu, tzu_a[lor] = tzu;
				tempi_a[lor] = tempi, tempj_a[lor] = tempj, tempk_a[lor] = tempk;
				iu_a[lor] = iu, ju_a[lor] = ju, ku_a[lor] = ku;
				float local_ele;
				for (uint ii = 0u; ii < Np; ii++) {
					const uint local_ind = compute_ind(tempj, tempi, tempk, d_Nx, d_Ny, d_N, d_Nx, d_Nxy);
					if (tz0 < ty0 && tz0 < tx0) {
						local_ele = compute_element(&tz0, &tc, LL[lor], tzu, ku, &tempk, &temp);
					}
					else if (ty0 < tx0) {
						local_ele = compute_element(&ty0, &tc, LL[lor], tyu, ju, &tempj, &temp);
					}
					else {
						local_ele = compute_element(&tx0, &tc, LL[lor], txu, iu, &tempi, &temp);
					}
					if (d_attenuation_correction) {
						jelppi += (local_ele * -d_atten[local_ind]);
					}
					if (local_sino > 0.f) {
						denominator_multi(local_ele, &axOSEM, &d_OSEM[local_ind]);
					}
					Np_n[lor]++;
					if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= d_Nx || tempj >= d_Ny || tempk >= d_Nz)
						break;
				}
				pass[lor] = true;
			}
			//pass[lor] = false;
		}
	}
	bool alku = true;
	for (ushort lor = 0u; lor < n_rays; lor++) {

		if (pass[lor]) {

			if (alku) {
				temp = 1.f / temp;
				if (d_attenuation_correction == 1u) {
					float n_r_summa = 0.f;
					for (ushort ln_r = 0u; ln_r < n_rays; ln_r++)
						n_r_summa += convert_float(pass[ln_r]);
					temp *= native_exp(jelppi / n_r_summa);
				}
				if (d_normalization == 1u)
					temp *= d_norm[idx];
				if (local_sino != 0.f) {
					if (axOSEM == 0.f) {
						axOSEM = d_epps;
					}
					else {
						axOSEM *= temp;
					}
					if (d_randoms == 1u)
						axOSEM += d_sc_ra[idx];
					axOSEM = local_sino / axOSEM;
				}
				alku = false;
			}
			if (tx0_a[lor] > 1e6f && ty0_a[lor] > 1e6f) {
				if (ty0_a[lor] > tx0_a[lor]) {

					const int tempk = tempk_a[lor];
					if (local_sino != 0.f) {
						for (uint k = 0; k < d_Nx; k++) {
							if (no_norm == 0u)
								atomicAdd_g_f(&d_Summ[tempk + k], convert_ulong_sat(d_dx * temp* TH));
							atomicAdd_g_f(&d_rhs_OSEM[tempk + k], convert_ulong_sat(d_dx * temp * axOSEM* TH));
						}
					}
					else {
						for (uint k = 0; k < d_Nx; k++) {
							atomicAdd_g_f(&d_Summ[tempk + k], convert_ulong_sat(d_dx * temp* TH));
						}
					}
				}
				else {
					const int tempk = tempk_a[lor];
					if (local_sino != 0.f) {
						for (uint k = 0; k < d_Ny; k++) {
							if (no_norm == 0u)
								atomicAdd_g_f(&d_Summ[tempk + k * d_Nx], convert_ulong_sat(d_dy * temp* TH));
							atomicAdd_g_f(&d_rhs_OSEM[tempk + k * d_Nx], convert_ulong_sat((d_dy * temp) * axOSEM* TH));
						}
					}
					else {
						for (uint k = 0; k < d_Ny; k++) {
							atomicAdd_g_f(&d_Summ[tempk + k * d_Nx], convert_ulong_sat(d_dy * temp* TH));
						}
					}
				}
			}
			else {
				float tx0 = tx0_a[lor];
				float ty0 = ty0_a[lor];
				float tz0 = tz0_a[lor];
				const float txu = txu_a[lor];
				const float tyu = tyu_a[lor];
				const float tzu = tzu_a[lor];
				int tempi = tempi_a[lor];
				int tempj = tempj_a[lor];
				int tempk = tempk_a[lor];
				const int iu = iu_a[lor];
				const int ju = ju_a[lor];
				const int ku = ku_a[lor];
				float tc = tc_a[lor];

				float local_ele;
				if (local_sino != 0.f) {
					for (uint ii = 0u; ii < Np_n[lor]; ii++) {
						const uint local_ind = compute_ind(tempj, tempi, tempk, d_Nx, d_Ny, d_N, d_Nx, d_Nxy);
						if (tz0 < ty0 && tz0 < tx0) {
							local_ele = compute_element_2nd(&tz0, &tc, LL[lor], tzu, ku, &tempk, temp);
							if (no_norm == 0u)
								atom_add(&d_Summ[local_ind], convert_ulong_sat(local_ele* TH));
							atom_add(&d_rhs_OSEM[local_ind], convert_ulong_sat(local_ele* axOSEM* TH));
						}
						else if (ty0 < tx0 && ty0 <= tz0) {
							local_ele = compute_element_2nd(&ty0, &tc, LL[lor], tyu, ju, &tempj, temp);
							if (no_norm == 0u)
								atom_add(&d_Summ[local_ind], convert_ulong_sat(local_ele* TH));
							atom_add(&d_rhs_OSEM[local_ind], convert_ulong_sat(local_ele* axOSEM* TH));
						}
						else if (tx0 <= ty0 && tx0 <= tz0) {
							local_ele = compute_element_2nd(&tx0, &tc, LL[lor], txu, iu, &tempi, temp);
							if (no_norm == 0u)
								atom_add(&d_Summ[local_ind], convert_ulong_sat(local_ele* TH));
							atom_add(&d_rhs_OSEM[local_ind], convert_ulong_sat(local_ele* axOSEM * TH));
						}
					}

				}
				else {
					for (uint ii = 0u; ii < Np_n[lor]; ii++) {
						const uint local_ind = compute_ind(tempj, tempi, tempk, d_Nx, d_Ny, d_N, d_Nx, d_Nxy);
						if (tz0 < ty0 && tz0 < tx0) {
							local_ele = compute_element_2nd(&tz0, &tc, LL[lor], tzu, ku, &tempk, temp);
							atom_add(&d_Summ[local_ind], convert_ulong_sat(local_ele* TH));
						}
						else if (ty0 < tx0 && ty0 <= tz0) {
							local_ele = compute_element_2nd(&ty0, &tc, LL[lor], tyu, ju, &tempj, temp);
							atom_add(&d_Summ[local_ind], convert_ulong_sat(local_ele* TH));
						}
						else if (tx0 <= ty0 && tx0 <= tz0) {
							local_ele = compute_element_2nd(&tx0, &tc, LL[lor], txu, iu, &tempi, temp);
							atom_add(&d_Summ[local_ind], convert_ulong_sat(local_ele* TH));
						}
					}
				}
			}
		}
	}
}



__kernel void summa(const __global ulong* d_Summ_device, __global ulong* d_Summ_local, const __global ulong* d_rhs_device, __global ulong* d_rhs_local,
	const uint im_dim, const uchar no_norm) {

	uint gid = get_global_id(0);

	for (uint i = gid; i < im_dim; i += get_global_size(0)) {
		if (no_norm == 0u)
			d_Summ_local[i] += d_Summ_device[i];
		d_rhs_local[i] += d_rhs_device[i];
	}
}


__kernel void mlem(const __global ulong* d_Summ, const __global ulong* d_rhs, __global float* d_mlem, const uint im_dim, const float d_epps) {

	uint gid = get_global_id(0);

	for (uint i = gid; i < im_dim; i += get_global_size(0)) {
		float rhs = convert_float(d_rhs[i]);
		float Summ = convert_float(d_Summ[i]);
		if (rhs != 0.f) {
			if (Summ == 0.f)
				d_mlem[i] = d_mlem[i] / d_epps * (rhs / TH);
			else
				d_mlem[i] = d_mlem[i] / (Summ / TH) * (rhs / TH);
		}
		else {
			if (Summ != 0.f)
				d_mlem[i] = d_mlem[i] / (Summ / TH) * d_epps;
		}
	}
}