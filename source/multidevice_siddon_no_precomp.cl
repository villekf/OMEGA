/**************************************************************************
* A matrix free improved Siddon's for multiple rays.
* This function calculates Summ = sum(A,1) (sum of every row) and
* rhs = A*(y./(A'*x)), where A is the system matrix, y the measurements
* and x the estimate/image.
*
* Used by implementations 2 and 3.
*
* This version goes through all the LORs and determines on-the-fly if they
* intersect with the voxel space. Uses (optionally) multiple rays.
*
* Compiler preprocessing is utilized heavily, for example all the 
* corrections are implemented as compiler preprocesses. The code for 
* specific correction is thus only applied if it has been selected. The
* kernels are always compiled on-the-fly, though when using same input 
* parameters the kernel should be loaded from cache leading to a slightly
* faster startup time.
*
* INPUTS:
* global_factor = a global correction factor, e.g. dead time
* d_epps = a small constant to prevent division by zero,
* d_N = d_Nx * d_Ny * d_Nz,
* d_Nx/y/z = image size in x/y/z- dimension,
* d_dz/x/y = distance between adjecent voxels in z/x/y-dimension,
* d_bz/x/y = distance from the pixel space to origin (z/x/y-dimension),
* d_bzb = part in parenthesis of equation (9) in [1] precalculated when
* k = Nz,
* d_maxxx/yy = maximum distance of the pixel space from origin in
* x/y-dimension,
* d_zmax = maximum value of d_zdet,
* d_NSlices = the number of image slices,
* d_size_x = the number of detector elements,
* d_TotSinos = Total number of sinograms,
* d_det_per_ring = number of detectors per ring,
* d_raw = if 1 then raw list-mode data is used otherwise sinogram
* data
* pRows = number of pseudo rings,
* d_Nxy = d_Nx * d_Ny,
* fp = if 1, then only forward projection is computed, if 2 only 
* backprojection, if 0 then both
* d_atten = attenuation data (images),
* d_norm = normalization coefficients,
* d_Summ = buffer for d_Summ,
* d_lor = number of pixels that each LOR traverses,
* d_pseudos = location of pseudo rings,
* d_x/y/z_det = detector x/y/z-coordinates,
* d_xy/zindex = for sinogram format they determine the detector
* indices corresponding to each sinogram bin (unused with raw data),
* d_L = detector numbers for raw data (unused for sinogram format),
* d_Sino = Sinogram/raw data,
* d_sc_ra = Randoms and/or scatter data,
* d_OSEM = buffer for OSEM/MLEM estimates,
* d_rhs_OSEM = buffer for OSEM/MLEM RHS elements,
* no_norm = If 1, normalization constant is not computed,
* m_size = Total number of LORs for this subset,
* cumsum = offset for input vector b in backprojection
*
* OUTPUTS:
* d_rhs_OSEM = rhs values for OSEM/MLEM,
* d_OSEM = OSEM/MLEM estimate,
* d_Summ = Sensitivity image
*
* [1] Jacobs, F., Sundermann, E., De Sutter, B., Christiaens, M. Lemahieu,
* I. (1998). A Fast Algorithm to Calculate the Exact Radiological Path
* through a Pixel or Voxel Space. Journal of computing and information
* technology, 6 (1), 89-94.
*
* Copyright (C) 2020  Ville-Veikko Wettenhovi
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
#ifdef ATOMIC
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
//#define TH 100000000000.f
#endif
#include "general_opencl_functions.h"
#define TYPE 0

// Matrix free Improved Siddon's algorithm
__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, 1, 1)))
void siddon_multi(const float global_factor, const float d_epps, const uint d_N, const uint d_Nx, const uint d_Ny, const uint d_Nz, const float d_dz, const float d_dx,
	const float d_dy, const float d_bz, const float d_bx, const float d_by, const float d_bzb, const float d_maxxx, const float d_maxyy,
	const float d_zmax, const float d_NSlices, const uint d_size_x, const ushort d_TotSinos, const uint d_det_per_ring, const uint d_pRows,
	const uint d_Nxy, const uchar fp, const float dc_z, const ushort n_rays, const float d_epsilon_mramla,
	const __global float* d_atten, __constant uint* d_pseudos, const __global float* d_x, const __global float* d_y, const __global float* d_zdet,
	__constant uchar* MethodList, const __global float* d_norm, const __global float* d_scat, __global CAST* d_Summ, const __global ushort* d_lor,
	const __global uint* d_xyindex, const __global ushort* d_zindex, const __global ushort* d_L, const __global float* d_Sino, const __global float* d_sc_ra, const __global float* d_OSEM,
#ifndef MBSREM
	__global CAST* d_rhs_OSEM, const uchar no_norm, const ulong m_size, const ulong cumsum
#else
	const uint d_alku, const uchar MBSREM_prepass, __global float* d_ACOSEM_lhs, __global float* d_Amin, __global CAST* d_co,
	__global CAST* d_aco, __global float* d_E, const ulong m_size, const RecMethodsOpenCL MethodListOpenCL
#endif
) {
	// Get the current global index
	uint idx = get_global_id(0);
	if (idx >= m_size)
		return;
	const float local_sino = (d_Sino[idx]);
#ifndef MBSREM
	if (no_norm == 1u && local_sino == 0.f)
		return;
#else
	const uchar no_norm = 0u;
#endif

#ifdef AF
#ifdef MBSREM
	float axCOSEM = 0.f;
	float axACOSEM = 0.f;
	float minimi = 1e8f;
	bool RHS = true;
#else
	bool RHS = local_sino > 0.f ? true : false;
	float ax[N_REKOS];
#pragma unroll N_REKOS
	for (uint kk = 0; kk < N_REKOS; kk++)
		ax[kk] = 0.f;
#endif
#else
	bool RHS = local_sino > 0.f ? true : false;
	float axOSEM = 0.f;
#endif
#ifndef AF
	if (fp == 2) {
		axOSEM = d_OSEM[idx + cumsum];
	}
#endif
	uint d_N0 = d_Nx;
	uint d_N1 = d_Ny;
	uint d_N2 = 1u;
	uint d_N3 = d_Nx;
	float jelppi = 0.f;
	float temp = 0.f;
	int tempi_a[N_RAYS];
	int tempj_a[N_RAYS];
	int tempk_a[N_RAYS];
	uint d_N0_a[N_RAYS];
	int iu_a[N_RAYS];
	int ju_a[N_RAYS];
	int ku_a[N_RAYS];
	float tx0_a[N_RAYS];
	float ty0_a[N_RAYS];
	float tz0_a[N_RAYS];
	float tc_a[N_RAYS];
	float txu_a[N_RAYS];
	float tyu_a[N_RAYS];
	float tzu_a[N_RAYS];
	float LL[N_RAYS];
	uint Np_n[N_RAYS];
	bool pass[N_RAYS];
	// Load the next detector index
	// raw list-mode data
#pragma unroll N_RAYS
	for (ushort lor = 0u; lor < N_RAYS; lor++) {
		float xs, xd, ys, yd, zs, zd;
#ifdef RAW
		get_detector_coordinates_raw_multiray(d_x, d_y, d_zdet, d_L, d_det_per_ring, idx, d_pseudos, d_pRows, &xs, &xd, &ys, &yd, &zs, &zd, lor + 1u, dc_z);
		// Sinogram data
#else
		get_detector_coordinates_multiray(d_xyindex, d_zindex, d_size_x, idx, d_TotSinos, &xs, &xd, &ys, &yd, &zs, &zd, d_x, d_y, d_zdet, lor + 1u, dc_z);
#endif
		// Calculate the x, y and z distances of the detector pair
		float y_diff = (yd - ys);
		float x_diff = (xd - xs);
		const float z_diff = (zd - zs);
		if ((y_diff == 0.f && x_diff == 0.f && z_diff == 0.f) || (y_diff == 0.f && x_diff == 0.f)) {
			//pass[lor] = false;
			return;
		}
		pass[lor] = false;
		uint Np = 0u;
		Np_n[lor] = 0u;
		// If the measurement is on a same ring
		if (fabs(z_diff) < 1e-6f && (fabs(y_diff) < 1e-6f || fabs(x_diff) < 1e-6f)) {
			float d_b, dd, d_d, d_d2;
			if (fabs(y_diff) < 1e-6f && yd <= d_maxyy && yd >= d_by && ys <= d_maxyy && ys >= d_by) {
				d_b = d_by;
				dd = yd;
				d_d = d_dy;
				d_d2 = d_dx;
				float xs_apu = xs;
				xs = ys;
				ys = xs_apu;
				float xdiff_apu = x_diff;
				x_diff = y_diff;
				y_diff = xdiff_apu;
				d_N0 = d_Ny;
				d_N1 = d_Nx;
				d_N2 = d_Ny;
				d_N3 = 1u;
				tx0_a[lor] = 1e7f, ty0_a[lor] = 1e9f;
				pass[lor] = true;
			}
			else if (fabs(x_diff) < 1e-6f && xd <= d_maxxx && xd >= d_bx && xs <= d_maxxx && xs >= d_bx) {
				d_b = d_bx;
				dd = xd;
				d_d = d_dx;
				d_d2 = d_dy;
				tx0_a[lor] = 1e9f, ty0_a[lor] = 1e7f;
				pass[lor] = true;
			}
			if (pass[lor]) {
				Np_n[lor] = d_N1;
				uint tempk = convert_uint((zs / d_zmax) * (d_NSlices - 1.f));
				uint apu = 0u;
				const float element = perpendicular_elements_multiray(d_b, d_d, d_N0, dd, d_d2, d_N1, d_atten, &apu, tempk, d_N2, d_N3, &jelppi);
				temp += element;
				tempk_a[lor] = apu;
#ifdef FP
#ifdef MBSREM
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && local_sino > 0.f && d_alku == 0u) {
					for (uint k = 0u; k < d_N1; k++)
						axCOSEM += (d_d * d_OSEM[apu + k * d_N3]);
				}
#else
				if (RHS) {
					for (uint k = 0u; k < d_N1; k++) {
#ifdef AF
						denominator(d_d, ax, apu + k * d_N3, d_N, d_OSEM);
#else
						axOSEM += (d_d * d_OSEM[apu + k * d_N3]);
#endif
					}
				}
#endif
#endif
			}
		}
		else {
			int tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;
			float txu = 0.f, tyu = 0.f, tzu = 0.f, tc = 0.f, tx0 = 1e8f, ty0 = 1e8f, tz0 = 1e8f;
			bool skip = false;
			if (fabs(z_diff) < 1e-6f) {
				tempk = convert_int((zs / d_zmax) * (d_NSlices - 1.f));
				skip = siddon_pre_loop_2D(d_bx, d_by, x_diff, y_diff, d_maxxx, d_maxyy, d_dx, d_dy, d_Nx, d_Ny, &tempi, &tempj, &txu, &tyu, &Np, TYPE,
					ys, xs, yd, xd, &tc, &iu, &ju, &tx0, &ty0);
			}
			else if (fabs(y_diff) < 1e-6f) {
				tempj = perpendicular_start(d_by, yd, d_dy, d_Ny);
				skip = siddon_pre_loop_2D(d_bx, d_bz, x_diff, z_diff, d_maxxx, d_bzb, d_dx, d_dz, d_Nx, d_Nz, &tempi, &tempk, &txu, &tzu, &Np, TYPE,
					zs, xs, zd, xd, &tc, &iu, &ku, &tx0, &tz0);
				//int apu_tempi = tempi;
				//float apu_txu = txu;
				//float apu_tx0 = tx0;
				//float apu_xdiff = x_diff;
				//float apu_xs = xs;
				//int apu_iu = iu;
				//iu = ju;
				//ju = apu_iu;
				//tempi = tempj;
				//tempj = apu_tempi;
				//txu = tyu;
				//tyu = apu_txu;
				//tx0 = ty0;
				//ty0 = apu_tx0;
				//x_diff = y_diff;
				//y_diff = apu_xdiff;
				//xs = ys;
				//ys = apu_xs;
				//d_N0 = d_Ny;
				//d_N1 = d_Nx;
				//d_N2 = d_Ny;
				//d_N3 = 1u;
				if (yd > d_maxyy || yd < d_by)
					skip = true;
			}
			else if (fabs(x_diff) < 1e-6f) {
				tempi = perpendicular_start(d_bx, xd, d_dx, d_Nx);
				skip = siddon_pre_loop_2D(d_by, d_bz, y_diff, z_diff, d_maxyy, d_bzb, d_dy, d_dz, d_Ny, d_Nz, &tempj, &tempk, &tyu, &tzu, &Np, TYPE,
					zs, ys, zd, yd, &tc, &ju, &ku, &ty0, &tz0);
				if (xd > d_maxxx || xd < d_bx)
					skip = true;
			}
			else {
				skip = siddon_pre_loop_3D(d_bx, d_by, d_bz, x_diff, y_diff, z_diff, d_maxxx, d_maxyy, d_bzb, d_dx, d_dy, d_dz, d_Nx, d_Ny, d_Nz, &tempi, &tempj, &tempk, &tyu, &txu, &tzu,
					&Np, TYPE, ys, xs, yd, xd, zs, zd, &tc, &iu, &ju, &ku, &tx0, &ty0, &tz0);
			}
			if (!skip) {
				pass[lor] = true;
				//return;
			}
			if (pass[lor]) {
				LL[lor] = native_sqrt(x_diff * x_diff + z_diff * z_diff + y_diff * y_diff);
				tx0_a[lor] = tx0, ty0_a[lor] = ty0, tz0_a[lor] = tz0, tc_a[lor] = tc;
				txu_a[lor] = txu, tyu_a[lor] = tyu, tzu_a[lor] = tzu;
				tempi_a[lor] = tempi, tempj_a[lor] = tempj, tempk_a[lor] = tempk;
				d_N0_a[lor] = d_N0;
				//uint temp_ijk = compute_ind(tempj, tempi, tempk, d_N0, d_N1, d_N, d_N3, d_Nxy);
				//temp_ijk_a[lor] = temp_ijk;
				iu_a[lor] = iu, ju_a[lor] = ju, ku_a[lor] = ku;
				float local_ele;
				for (uint ii = 0u; ii < Np; ii++) {
					const uint local_ind = compute_ind(tempj, tempi, tempk, d_N0, d_N1, d_N, d_N0, d_Nxy);
					if (tz0 < ty0 && tz0 < tx0) {
						local_ele = compute_element(&tz0, &tc, LL[lor], tzu, ku, &tempk, &temp);
					}
					else if (ty0 < tx0) {
						local_ele = compute_element(&ty0, &tc, LL[lor], tyu, ju, &tempj, &temp);
					}
					else {
						local_ele = compute_element(&tx0, &tc, LL[lor], txu, iu, &tempi, &temp);
					}
#ifdef ATN
					jelppi += (local_ele * -d_atten[local_ind]);
#endif
#ifdef FP
#ifdef MBSREM
					if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && local_sino > 0.f && d_alku == 0u)
						axCOSEM += (local_ele * d_OSEM[local_ind]);
#else
					if (RHS) {
#ifdef AF
						denominator(local_ele, ax, local_ind, d_N, d_OSEM);
#else
						denominator_multi(local_ele, &axOSEM, &d_OSEM[local_ind]);
#endif
					}
#endif
#endif
					Np_n[lor]++;
					if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= d_N0 || tempj >= d_N1 || tempk >= d_Nz)
						break;
				}
				//pass[lor] = true;
			}
		}
	}
	bool alku = true;
	for (ushort lor = 0u; lor < N_RAYS; lor++) {
		if (pass[lor]) {
			if (alku) {
				//if (temp == 0.f)
				//	return;
				temp = 1.f / temp;
#ifdef ATN
				float n_r_summa = 0.f;
				for (ushort ln_r = 0u; ln_r < N_RAYS; ln_r++)
					n_r_summa += convert_float(pass[ln_r]);
				temp *= native_exp(jelppi / n_r_summa);
#endif
#ifdef NORM
				temp *= d_norm[idx];
#endif
#ifdef SCATTER
				temp *= d_scat[idx];
#endif
				temp *= global_factor;
#ifdef FP
#ifdef MBSREM
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && local_sino > 0.f && d_alku == 0u) {
					if (axCOSEM == 0.f)
						axCOSEM = d_epps;
					else
						axCOSEM *= temp;
#ifdef RANDOMS
					axCOSEM += d_sc_ra[idx];
#endif
					axCOSEM = local_sino / axCOSEM;
				}
#else
#ifdef AF
				if (RHS) {
					nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx);
				}
#else
				if (RHS) {
					if (axOSEM == 0.f) {
						axOSEM = d_epps;
					}
					else {
						axOSEM *= temp;
					}
#ifdef RANDOMS
					axOSEM += d_sc_ra[idx];
#endif
					if (fp == 1)
						d_rhs_OSEM[idx] = axOSEM;
					else {
						axOSEM = local_sino / axOSEM;
					}
				}
				if (fp == 1)
					return;
#endif
#endif
#endif
				alku = false;
			}
			if (tx0_a[lor] > 1e6f && ty0_a[lor] > 1e6f) {
				const uint tempk = tempk_a[lor];
				//if (tempk >= d_N)
				//	return;
				if (ty0_a[lor] > tx0_a[lor]) {
					if (RHS) {
						for (uint k = 0; k < Np_n[lor]; k++) {
#ifdef MBSREM

							if (d_alku == 0u) {
								if (MBSREM_prepass == 1)
#ifdef ATOMIC
									atom_add(&d_Summ[tempk + k], convert_ulong_sat(d_dx * temp * TH));
#else
									atomicAdd_g_f(&d_Summ[tempk + k], (d_dx * temp));
#endif
								if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
									minimi = d_dx * temp;
									d_E[idx] += d_dx * temp;
								}
								if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
#ifdef ATOMIC
									atom_add(&d_co[tempk + k], convert_ulong_sat(axCOSEM * d_dx * temp * TH));
#else
									atomicAdd_g_f(&d_co[tempk + k], axCOSEM * d_dx * temp);
#endif
								if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
#ifdef ATOMIC
									atom_add(&d_aco[tempk + k], convert_ulong_sat(axCOSEM * d_dx * temp * TH));
#else
									atomicAdd_g_f(&d_aco[tempk + k], axCOSEM * d_dx * temp);
#endif
							}
							if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u)
								axACOSEM += (d_dx * temp * d_OSEM[tempk + k]);
#else
							if (no_norm == 0u)
#ifdef ATOMIC
								atom_add(&d_Summ[tempk + k], convert_ulong_sat(d_dx * temp * TH));
#else
								atomicAdd_g_f(&d_Summ[tempk + k], (d_dx * temp));
#endif
#ifdef AF
							rhs(MethodList, d_dx * temp, ax, tempk + k, d_N, d_rhs_OSEM);
#else
#ifdef ATOMIC
							atom_add(&d_rhs_OSEM[tempk + k], convert_ulong_sat(d_dx * temp * axOSEM * TH));
#else
							atomicAdd_g_f(&d_rhs_OSEM[tempk + k], ((d_dx * temp) * axOSEM));
#endif
#endif
#endif
						}
					}
					else {
						for (uint k = 0; k < Np_n[lor]; k++) {
#ifdef ATOMIC
							atom_add(&d_Summ[tempk + k], convert_ulong_sat(d_dx * temp * TH));
#else
							atomicAdd_g_f(&d_Summ[tempk + k], (d_dx * temp));
#endif
						}
					}
				}
				else {
					if (RHS) {
						for (uint k = 0; k < Np_n[lor]; k++) {
#ifdef MBSREM

							if (d_alku == 0u) {
								if (MBSREM_prepass == 1)
#ifdef ATOMIC
									atom_add(&d_Summ[tempk + k * d_N0], convert_ulong_sat(d_dy * temp * TH));
#else
									atomicAdd_g_f(&d_Summ[tempk + k * d_N0], (d_dy * temp));
#endif
								if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
									minimi = d_dy * temp;
									d_E[idx] += d_dy * temp;
								}
								if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
#ifdef ATOMIC
									atom_add(&d_co[tempk + k * d_N0], convert_ulong_sat(axCOSEM * d_dy * temp * TH));
#else
									atomicAdd_g_f(&d_co[tempk + k * d_N0], axCOSEM * d_dy * temp);
#endif
								if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
#ifdef ATOMIC
									atom_add(&d_aco[tempk + k * d_N0], convert_ulong_sat(axCOSEM * d_dy * temp * TH));
#else
									atomicAdd_g_f(&d_aco[tempk + k * d_N0], axCOSEM * d_dy * temp);
#endif
							}
							if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u)
								axACOSEM += (d_dy * temp * d_OSEM[tempk + k * d_Nx]);
#else
#ifdef AF
							rhs(MethodList, d_dy * temp, ax, tempk + k * d_N0, d_N, d_rhs_OSEM);
#else
#ifdef ATOMIC
							atom_add(&d_rhs_OSEM[tempk + k * d_N0], convert_ulong_sat(d_dy * temp * axOSEM * TH));
#else
							atomicAdd_g_f(&d_rhs_OSEM[tempk + k * d_N0], ((d_dy * temp) * axOSEM));
#endif
#endif
							if (no_norm == 0u)
#ifdef ATOMIC
								atom_add(&d_Summ[tempk + k * d_N0], convert_ulong_sat(d_dy * temp * TH));
#else
								atomicAdd_g_f(&d_Summ[tempk + k * d_N0], (d_dy * temp));
#endif
#endif
						}
					}
					else {
						for (uint k = 0; k < Np_n[lor]; k++) {
#ifdef ATOMIC
							atom_add(&d_Summ[tempk + k * d_N0], convert_ulong_sat(d_dy * temp * TH));
#else
							atomicAdd_g_f(&d_Summ[tempk + k * d_N0], (d_dy * temp));
#endif
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
				if (RHS) {
					for (uint ii = 0u; ii < Np_n[lor]; ii++) {
						const uint local_ind = compute_ind(tempj, tempi, tempk, d_N0, d_N1, d_N, d_N0_a[lor], d_Nxy);
						if (tz0 < ty0 && tz0 < tx0) {
							local_ele = compute_element_2nd(&tz0, &tc, LL[lor], tzu, ku, &tempk, temp);
						}
						else if (ty0 < tx0 && ty0 <= tz0) {
							local_ele = compute_element_2nd(&ty0, &tc, LL[lor], tyu, ju, &tempj, temp);
						}
						else if (tx0 <= ty0 && tx0 <= tz0) {
							local_ele = compute_element_2nd(&tx0, &tc, LL[lor], txu, iu, &tempi, temp);
						}
#ifdef MBSREM
						if (d_alku == 0u) {
							if (MBSREM_prepass == 1)
#ifdef ATOMIC
								atom_add(&d_Summ[local_ind], convert_ulong_sat(local_ele * TH));
#else
								atomicAdd_g_f(&d_Summ[local_ind], local_ele);
#endif
							if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
								if (local_ele < minimi && local_ele > 0.f)
									minimi = local_ele;
								d_E[idx] += local_ele;
							}
							if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
#ifdef ATOMIC
								atom_add(&d_co[local_ind], convert_ulong_sat(axCOSEM * local_ele * TH));
#else
								atomicAdd_g_f(&d_co[local_ind], axCOSEM * local_ele);
#endif
							if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
#ifdef ATOMIC
								atom_add(&d_aco[local_ind], convert_ulong_sat(axCOSEM * local_ele * TH));
#else
								atomicAdd_g_f(&d_aco[local_ind], axCOSEM * local_ele);
#endif
						}
						if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u)
							axACOSEM += (local_ele * d_OSEM[local_ind]);
#else
						if (no_norm == 0u)
#ifdef ATOMIC
							atom_add(&d_Summ[local_ind], convert_ulong_sat(local_ele * TH));
#else
							atomicAdd_g_f(&d_Summ[local_ind], local_ele);
#endif
#ifdef AF
						rhs(MethodList, local_ele, ax, local_ind, d_N, d_rhs_OSEM);
#else
#ifdef ATOMIC
						atom_add(&d_rhs_OSEM[local_ind], convert_ulong_sat(local_ele * axOSEM * TH));
#else
						atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * axOSEM));
#endif
#endif
#endif
					}
				}
				else {
					for (uint ii = 0u; ii < Np_n[lor]; ii++) {
						const uint local_ind = compute_ind(tempj, tempi, tempk, d_N0, d_N1, d_N, d_N0_a[lor], d_Nxy);
						if (tz0 < ty0 && tz0 < tx0) {
							local_ele = compute_element_2nd(&tz0, &tc, LL[lor], tzu, ku, &tempk, temp);
						}
						else if (ty0 < tx0 && ty0 <= tz0) {
							local_ele = compute_element_2nd(&ty0, &tc, LL[lor], tyu, ju, &tempj, temp);
						}
						else if (tx0 <= ty0 && tx0 <= tz0) {
							local_ele = compute_element_2nd(&tx0, &tc, LL[lor], txu, iu, &tempi, temp);
						}
#ifdef ATOMIC
						atom_add(&d_Summ[local_ind], convert_ulong_sat(local_ele * TH));
#else
						atomicAdd_g_f(&d_Summ[local_ind], local_ele);
#endif
					}
				}
			}
		}
	}
#ifdef MBSREM
	if (!alku) {
		if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1 && d_alku == 0u)
			d_Amin[idx] = minimi;
		if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u) {
#ifdef RANDOMS
			axACOSEM += d_sc_ra[idx];
#endif
			d_ACOSEM_lhs[idx] = axACOSEM;
		}
	}
#endif
}


#ifndef AF
__kernel void summa(const __global CAST* d_Summ_device, __global CAST* d_Summ_local, const __global CAST* d_rhs_device, __global CAST* d_rhs_local,
	const uint im_dim, const uchar no_norm) {

	uint gid = get_global_id(0);

	for (uint i = gid; i < im_dim; i += get_global_size(0)) {
		if (no_norm == 0u)
			d_Summ_local[i] += d_Summ_device[i];
		d_rhs_local[i] += d_rhs_device[i];
	}
}

__kernel void mlem(const __global CAST* d_Summ, const __global CAST* d_rhs, __global float* d_mlem, const uint im_dim, const float d_epps) {

	uint gid = get_global_id(0);

	for (uint i = gid; i < im_dim; i += get_global_size(0)) {
#ifdef ATOMIC
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
#else
		if (d_rhs[i] != 0.f) {
			if (d_Summ[i] == 0.f)
				d_mlem[i] = d_mlem[i] / d_epps * d_rhs[i];
			else
				d_mlem[i] = d_mlem[i] / d_Summ[i] * d_rhs[i];
		}
		else {
			if (d_Summ[i] != 0.f)
				d_mlem[i] = d_mlem[i] / d_Summ[i] * d_epps;
		}
#endif
	}
}


#ifdef PSF
__kernel void Convolution3D(const __global CAST* input, __global CAST* output,
	__constant float* convolution_window, int window_size_x, int window_size_y, int window_size_z) {
	int4 ind = (int4)(get_global_id(0), get_global_id(1), get_global_id(2), 0);
	int4 ind_uus = (0, 0, 0, 0);
	float result = 0.f;
	const uint Nyx = get_global_size(0) * get_global_size(1);
	//int radius_x = floor((float)window_size_x / 2.0f);
	//int radius_y = floor((float)window_size_y / 2.0f);
	//int radius_z = floor((float)window_size_z / 2.0f);
	int c = 0;
	for (int k = -window_size_z; k <= window_size_z; k++) {
		for (int j = -window_size_y; j <= window_size_y; j++) {
			for (int i = -window_size_x; i <= window_size_x; i++) {
				ind_uus.x = ind.x + i;
				ind_uus.y = ind.y + j;
				ind_uus.z = ind.z + k;
				if (ind_uus.x >= get_global_size(0))
					ind_uus.x = ind.x - i + 1;
				else if (ind_uus.x < 0)
					ind_uus.x = ind.x - (i + 1);
				if (ind_uus.y >= get_global_size(1))
					ind_uus.y = ind.y - j + 1;
				else if (ind_uus.y < 0)
					ind_uus.y = ind.y - (j + 1);
				if (ind_uus.z >= get_global_size(2))
					ind_uus.z = ind.z - k + 1;
				else if (ind_uus.z < 0)
					ind_uus.z = ind.z - (k + 1);
				uint indeksi = ind_uus.x + ind_uus.y * get_global_size(0) + ind_uus.z * Nyx;
#ifdef ATOMIC
				float p = convert_float(input[indeksi]) / TH;
#else
				float p = input[indeksi];
#endif
				p *= convolution_window[c];
				result += p;
				c++;
			}
		}
	}
#ifdef ATOMIC
	output[ind.x + ind.y * get_global_size(0) + ind.z * Nyx] = convert_ulong_sat(result * TH);
#else
	output[ind.x + ind.y * get_global_size(0) + ind.z * Nyx] = result;
#endif
}


__kernel void Convolution3D_f(const __global float* input, __global float* output,
	__constant float* convolution_window, int window_size_x, int window_size_y, int window_size_z) {
	int4 ind = (int4)(get_global_id(0), get_global_id(1), get_global_id(2), 0);
	int4 ind_uus = (0, 0, 0, 0);
	float result = 0.f;
	const uint Nyx = get_global_size(0) * get_global_size(1);
	//int radius_x = floor((float)window_size_x / 2.0f);
	//int radius_y = floor((float)window_size_y / 2.0f);
	//int radius_z = floor((float)window_size_z / 2.0f);
	int c = 0;
	for (int k = -window_size_z; k <= window_size_z; k++) {
		for (int j = -window_size_y; j <= window_size_y; j++) {
			for (int i = -window_size_x; i <= window_size_x; i++) {
				ind_uus.x = ind.x + i;
				ind_uus.y = ind.y + j;
				ind_uus.z = ind.z + k;
				if (ind_uus.x >= get_global_size(0))
					ind_uus.x = ind.x - i + 1;
				else if (ind_uus.x < 0)
					ind_uus.x = ind.x - (i + 1);
				if (ind_uus.y >= get_global_size(1))
					ind_uus.y = ind.y - j + 1;
				else if (ind_uus.y < 0)
					ind_uus.y = ind.y - (j + 1);
				if (ind_uus.z >= get_global_size(2))
					ind_uus.z = ind.z - k + 1;
				else if (ind_uus.z < 0)
					ind_uus.z = ind.z - (k + 1);
				uint indeksi = ind_uus.x + ind_uus.y * get_global_size(0) + ind_uus.z * Nyx;
				float p = input[indeksi];
				p *= convolution_window[c];
				result += p;
				c++;
			}
		}
	}
	output[ind.x + ind.y * get_global_size(0) + ind.z * Nyx] = result;
}

__kernel void vectorDiv(const __global float* input, __global float* output) {
	uint id = get_global_id(0);
	output[id] = output[id] / input[id];
}

__kernel void vectorMult(const __global float* input, __global float* output) {
	uint id = get_global_id(0);
	output[id] *= input[id];
}
#endif
#endif

#ifdef NLM_
__kernel void NLM(__global float* grad, const __global float* u, const __global float* u_ref, __constant float* gaussian, const int search_window_x, const int search_window_y,
	const int search_window_z, const int patch_window_x, const int patch_window_y, const int patch_window_z, const uint Nx, const uint Ny, const uint Nz,
	const float h, const float epps, const int Nxy, const int min_x, const int max_x, const int min_y, const int max_y, const int min_z,
	const int max_z, const int type) {

	int n = get_global_id(0);
	const int z = n / Nxy;
	const int y = (n - z * Nxy) / convert_int(Nx);
	const int x = n - z * Nxy - y * convert_int(Nx);
	if (z < min_z || z >= max_z || x < min_x || x >= max_x || y < min_y || y >= max_y)
		return;
	float weight_sum = 0.f;
	float output = 0.f;
	const float uj = u[n];
	for (int k = -search_window_z; k <= search_window_z; k++) {
		const int z_n = z + k;
		for (int j = -search_window_y; j <= search_window_y; j++) {
			const int y_n = y + j;
			for (int i = -search_window_x; i <= search_window_x; i++) {
				const int x_n = x + i;
				const int dim_n = z_n * Nxy + y_n * convert_int(Nx) + x_n;
				const float uk = u[dim_n];
				float distance = 0.f;
				float weight = 0.f;

				for (int pz = -patch_window_z; pz <= patch_window_z; pz++) {
					const int z_k = (z_n + pz) * Nxy;
					const int z_j = (z + pz) * Nxy;
					for (int py = -patch_window_y; py <= patch_window_y; py++) {
						const int y_k = (y_n + py) * convert_int(Nx);
						const int y_j = (y + py) * convert_int(Nx);
						int dim_g = (pz + patch_window_z) * (patch_window_x * 2 + 1) * (patch_window_y * 2 + 1) + (py + patch_window_y) * (patch_window_x * 2 + 1);
						for (int px = -patch_window_x; px <= patch_window_x; px++) {
							const float gg = gaussian[dim_g++];
							//const float gg = 1.;
							const int x_k = x_n + px;
							const int dim_k = z_k + y_k + x_k;
							const float Pj = u_ref[dim_k];
							const int x_j = x + px;
							const int dim_j = z_j + y_j + x_j;
							const float Pk = u_ref[dim_j];
							distance += gg * (Pj - Pk) * (Pj - Pk);
						}
					}
				}
				weight = exp(-distance / h);
				weight_sum += weight;
				if (type == 2)
					output += weight * uk;
				else if (type == 0) {
					output += (weight * (uj - uk));
				}
				else {
					output += ((weight * (uj - uk)) / sqrt(weight * (uj - uk) * (uj - uk) + epps));
				}
			}
		}
	}
	weight_sum = 1.f / weight_sum;
	output *= weight_sum;

	grad[n] = output;

}
#endif