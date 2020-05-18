/**************************************************************************
* Special functions for the 3D orthogonal distance-based ray tracer (CUDA).
*
* Copyright (C) 2020 Ville-Veikko Wettenhovi
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) ad_N1 later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT Ad_N1 WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <https://www.gnu.org/licenses/>.
***************************************************************************/
#pragma once
#include "general_orth_cuda_functions.cuh"

// Denominator (forward projection), multi-GPU version, orthogonal distance based ray tracer
__device__ void orth_distance_multi_3D(int tempi, const unsigned int d_N0, const unsigned int d_N4, const float y_diff, const float x_diff, const float z_diff,
	const float* y_center, const float* x_center, const float* z_center, float* temp, const unsigned int d_N2, int tempj, int tempk,
	const float local_sino, const float* d_OSEM, const float xs, const float ys, const float zs, const unsigned int d_Nxy, const float kerroin,
	const bool no_norm, const bool RHS, const bool SUMMA, CAST* Summ, const unsigned int d_N1, const unsigned int d_N3, const char
	start, int iu, int ju, unsigned char loppu, const float bmin, const float bmax, const float Vmax, const float* V, float* d_store_elements,
#ifdef MBSREM
	unsigned int* d_store_indices, unsigned int* ind, float* ax, const RecMethodsOpenCL MethodListOpenCL, float* d_E, const unsigned int idx,
	CAST* d_co, CAST* d_aco, float* minimi, const unsigned char MBSREM_prepass, float* axCOSEM, const unsigned int d_alku
#else
	unsigned int* d_store_indices, unsigned int* ind, CAST* d_rhs_OSEM, const unsigned int im_dim, const unsigned char* MethodList, float* ax
#endif
) {
	int yy1 = 0;
	int uu1 = 0;
	int yy2 = 0;
	int uu2 = 0;
	int alku_y1 = tempj;
	int alku_x1 = tempi;
	int alku_y2 = tempj - 1;
	int alku_x2 = tempi - 1;
	bool breikki1 = false;
	bool breikki2 = false;
	bool breikki3 = false;
	bool breikki4 = false;
	float zs_apu = zs;
	float x_diff_apu = x_diff;
	float ys_apu = ys;
	for (int zz = tempk; zz < start; zz++) {
		zs_apu = (z_center[zz] - zs);
		x_diff_apu = x_diff * zs_apu;
		zs_apu *= y_diff;
		yy1 = 0;
		uu1 = 0;
		yy2 = 0;
		uu2 = 0;
		alku_y1 = tempj;
		alku_x1 = tempi;
		alku_y2 = tempj - 1;
		alku_x2 = tempi - 1;
		breikki1 = false;
		breikki2 = false;
		breikki3 = false;
		breikki4 = false;
		for (yy1 = alku_y1; yy1 < d_N1; yy1++) {
			ys_apu = y_center[yy1] - ys;
			float zs_apu2 = zs_apu - z_diff * ys_apu;
			ys_apu *= x_diff;
			computeOrthVoxelIncreasing(xs, ys_apu, zs_apu2, x_diff_apu, y_diff, z_diff, kerroin, x_center, y_center[yy1], z_center[zz], temp, local_sino, ax,
				&breikki1, alku_x1, zz, yy1, &uu1, no_norm, RHS, SUMMA, d_Nxy, d_N3, d_N2, d_N0, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
#ifdef MBSREM
				d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
				d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
#endif
			computeOrthVoxelDecreasing(xs, ys_apu, zs_apu2, x_diff_apu, y_diff, z_diff, kerroin, x_center, y_center[yy1], z_center[zz], temp, local_sino, ax,
				&breikki2, alku_x2, zz, yy1, &uu2, no_norm, RHS, SUMMA, d_Nxy, d_N3, d_N2, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
#ifdef MBSREM
				d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
				d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
#endif
			if (iu > 0) {
				if (ju > 0) {
					alku_x1 = uu1 - 1;
					alku_x2 = uu1 - 2;
				}
				else {
					alku_x2 = uu2 + 1;
					alku_x1 = uu2 + 2;
				}
			}
			else {
				if (ju > 0) {
					alku_x2 = uu2 + 1;
					alku_x1 = uu2 + 2;
				}
				else {
					alku_x1 = uu1 - 1;
					alku_x2 = uu1 - 2;
				}
			}
			if (breikki1 && breikki2) {
				if (yy1 == alku_y1)
					breikki3 = true;
				break;
			}
		}
		breikki1 = false;
		breikki2 = false;
		alku_x1 = tempi;
		alku_x2 = tempi - 1;
		for (yy2 = alku_y2; yy2 >= 0; yy2--) {
			ys_apu = y_center[yy2] - ys;
			float zs_apu2 = zs_apu - z_diff * ys_apu;
			ys_apu *= x_diff;
			computeOrthVoxelIncreasing(xs, ys_apu, zs_apu2, x_diff_apu, y_diff, z_diff, kerroin, x_center, y_center[yy2], z_center[zz], temp, local_sino, ax,
				&breikki1, alku_x1, zz, yy2, &uu1, no_norm, RHS, SUMMA, d_Nxy, d_N3, d_N2, d_N0, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
#ifdef MBSREM
				d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
				d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
#endif
			computeOrthVoxelDecreasing(xs, ys_apu, zs_apu2, x_diff_apu, y_diff, z_diff, kerroin, x_center, y_center[yy2], z_center[zz], temp, local_sino, ax,
				&breikki2, alku_x2, zz, yy2, &uu2, no_norm, RHS, SUMMA, d_Nxy, d_N3, d_N2, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
#ifdef MBSREM
				d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
				d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
#endif
			if (iu > 0) {
				if (ju > 0) {
					alku_x2 = uu2 + 1;
					alku_x1 = uu2 + 2;
				}
				else {
					alku_x1 = uu1 - 1;
					alku_x2 = uu1 - 2;
				}
			}
			else {
				if (ju > 0) {
					alku_x1 = uu1 - 1;
					alku_x2 = uu1 - 2;
				}
				else {
					alku_x2 = uu2 + 1;
					alku_x1 = uu2 + 2;
				}
			}
			if (breikki1 && breikki2) {
				if (yy2 == alku_y2)
					breikki4 = true;
				break;
			}
		}
		if (breikki3 && breikki4) {
			break;
		}
	}
	for (int zz = tempk - 1; zz >= loppu; zz--) {
		zs_apu = (z_center[zz] - zs);
		x_diff_apu = x_diff * zs_apu;
		zs_apu *= y_diff;
		yy1 = 0;
		uu1 = 0;
		yy2 = 0;
		uu2 = 0;
		alku_y1 = tempj;
		alku_x1 = tempi;
		alku_y2 = tempj - 1;
		alku_x2 = tempi - 1;
		breikki1 = false;
		breikki2 = false;
		breikki3 = false;
		breikki4 = false;
		for (yy1 = alku_y1; yy1 < d_N1; yy1++) {
			ys_apu = y_center[yy1] - ys;
			float zs_apu2 = zs_apu - z_diff * ys_apu;
			ys_apu *= x_diff;
			computeOrthVoxelIncreasing(xs, ys_apu, zs_apu2, x_diff_apu, y_diff, z_diff, kerroin, x_center, y_center[yy1], z_center[zz], temp, local_sino, ax,
				&breikki1, alku_x1, zz, yy1, &uu1, no_norm, RHS, SUMMA, d_Nxy, d_N3, d_N2, d_N0, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
#ifdef MBSREM
				d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
				d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
#endif
			computeOrthVoxelDecreasing(xs, ys_apu, zs_apu2, x_diff_apu, y_diff, z_diff, kerroin, x_center, y_center[yy1], z_center[zz], temp, local_sino, ax,
				&breikki2, alku_x2, zz, yy1, &uu2, no_norm, RHS, SUMMA, d_Nxy, d_N3, d_N2, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
#ifdef MBSREM
				d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
				d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
#endif
			if (iu > 0) {
				if (ju > 0) {
					alku_x1 = uu1 - 1;
					alku_x2 = uu1 - 2;
				}
				else {
					alku_x2 = uu2 + 1;
					alku_x1 = uu2 + 2;
				}
			}
			else {
				if (ju > 0) {
					alku_x2 = uu2 + 1;
					alku_x1 = uu2 + 2;
				}
				else {
					alku_x1 = uu1 - 1;
					alku_x2 = uu1 - 2;
				}
			}
			if (breikki1 && breikki2) {
				if (yy1 == alku_y1)
					breikki3 = true;
				break;
			}
		}
		breikki1 = false;
		breikki2 = false;
		alku_x1 = tempi;
		alku_x2 = tempi - 1;
		for (yy2 = alku_y2; yy2 >= 0; yy2--) {
			ys_apu = y_center[yy2] - ys;
			float zs_apu2 = zs_apu - z_diff * ys_apu;
			ys_apu *= x_diff;
			computeOrthVoxelIncreasing(xs, ys_apu, zs_apu2, x_diff_apu, y_diff, z_diff, kerroin, x_center, y_center[yy2], z_center[zz], temp, local_sino, ax,
				&breikki1, alku_x1, zz, yy2, &uu1, no_norm, RHS, SUMMA, d_Nxy, d_N3, d_N2, d_N0, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
#ifdef MBSREM
				d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
				d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
#endif
			computeOrthVoxelDecreasing(xs, ys_apu, zs_apu2, x_diff_apu, y_diff, z_diff, kerroin, x_center, y_center[yy2], z_center[zz], temp, local_sino, ax,
				&breikki2, alku_x2, zz, yy2, &uu2, no_norm, RHS, SUMMA, d_Nxy, d_N3, d_N2, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
#ifdef MBSREM
				d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
				d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
#endif
			if (iu > 0) {
				if (ju > 0) {
					alku_x2 = uu2 + 1;
					alku_x1 = uu2 + 2;
				}
				else {
					alku_x1 = uu1 - 1;
					alku_x2 = uu1 - 2;
				}
			}
			else {
				if (ju > 0) {
					alku_x1 = uu1 - 1;
					alku_x2 = uu1 - 2;
				}
				else {
					alku_x2 = uu2 + 1;
					alku_x1 = uu2 + 2;
				}
			}
			if (breikki1 && breikki2) {
				if (yy2 == alku_y2)
					breikki4 = true;
				break;
			}
		}
		if (breikki3 && breikki4) {
			break;
		}
	}
}

__device__ void orth_distance_perpendicular_multi_3D(const float* center1, const float center2, const float* z_center,
	float* temp, const unsigned int d_attenuation_correction, const unsigned int d_normalization, float* ax, const float d_b, const float d, const float d_d1,
	const unsigned int d_N1, const unsigned int d_N2, const unsigned int z_loop, const float* d_atten, const float d_norm, const float local_sino, const unsigned int d_N,
	const unsigned int d_NN, const float* d_OSEM, const float xs, const float ys, const float zs, const float xl, const float yl, const float zl,
	const float crystal_size_z, const unsigned int d_N1x, const unsigned int d_N4, const unsigned char no_norm, CAST* Summ,
	const bool FPbool, const bool RHS, const float global_factor, const float bmin, const float bmax, const float Vmax, const float* V,
#ifdef MBSREM
	const RecMethodsOpenCL MethodListOpenCL, const unsigned int d_alku, float* axCOSEM, 
	float* d_E, CAST* d_co, CAST* d_aco, float* minimi, const unsigned char MBSREM_prepass,
	const float* d_sc_ra, const unsigned int d_randoms, float* d_Amin, float* d_ACOSEM_lhs, const unsigned int idx
#else
	CAST* d_rhs_OSEM, const unsigned int im_dim, const unsigned char* MethodList
#endif
) {
	//const unsigned int zz = z_loop * d_N2 * d_N1;
	const unsigned int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	float jelppi = 0.f;
	for (int zz = (int)(z_loop); zz >= 0; zz--) {
		for (int uu = (int)(apu); uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D_per(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
#ifdef VOL
			if (local_ele >= bmax)
				break;
			if (local_ele < bmin)
				local_ele = Vmax;
			else
				local_ele = V[__float2uint_rn((local_ele - bmin) * CC)];
#else
			if (local_ele <= THR)
				break;
#endif
			unsigned int local_ind = uu * d_N + zz * d_N1x;
			if (FPbool) {
				*temp += (local_ele * d_N2);
#if defined(ATN) || defined(MBSREM)
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (d_attenuation_correction && zz == (int)(z_loop) && uu == (int)(apu))
						jelppi += (d_d1 * -d_atten[local_ind]);
#ifdef MBSREM
					if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
						*axCOSEM += (local_ele * d_OSEM[local_ind]);
					}
#else
					if (local_sino > 0.f) {
						denominator(local_ele, ax, local_ind, im_dim, d_OSEM);
					}
#endif
					local_ind += d_NN;
				}
#endif
			}
			else if (RHS) {
				local_ele *= *temp;
#ifdef MBSREM
				if (d_alku == 0 && (MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
					if (local_ele < *minimi && local_ele > 0.f)
						*minimi = local_ele;
					d_E[idx] += local_ele;
				}
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (d_alku == 0) {
						if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
#ifdef ATOMIC
							atomicAdd(&d_co[local_ind], __float2ull_rn(*axCOSEM * local_ele * TH));
#else
							atomicAdd(&d_co[local_ind], *axCOSEM * local_ele);
#endif
						if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
#ifdef ATOMIC
							atomicAdd(&d_aco[local_ind], __float2ull_rn(*axCOSEM * local_ele * TH));
#else
							atomicAdd(&d_aco[local_ind], *axCOSEM * local_ele);
#endif
						if (MBSREM_prepass == 1)
#ifdef ATOMIC
							atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
							atomicAdd(&Summ[local_ind], local_ele);
#endif
					}
					else
						*ax += (local_ele * d_OSEM[local_ind]);
					local_ind += d_NN;
				}
#else
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
#ifdef ATOMIC
						atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
						atomicAdd(&Summ[local_ind], local_ele);
#endif
					rhs(MethodList, local_ele, ax, local_ind, im_dim, d_rhs_OSEM);
					local_ind += d_NN;
				}
#endif
			}
			else {
				local_ele *= *temp;
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
#ifdef ATOMIC
					atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
					atomicAdd(&Summ[local_ind], local_ele);
#endif
					local_ind += d_NN;
				}
			}
		}
		for (unsigned int uu = apu + 1; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D_per(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
#ifdef VOL
			if (local_ele >= bmax)
				break;
			if (local_ele < bmin)
				local_ele = Vmax;
			else
				local_ele = V[__float2uint_rn((local_ele - bmin) * CC)];
#else
			if (local_ele <= THR)
				break;
#endif
			unsigned int local_ind = uu * d_N + zz * d_N1x;
			if (FPbool) {
				*temp += (local_ele * d_N2);
#ifdef MBSREM
				if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
					*axCOSEM += (local_ele * d_OSEM[local_ind]);
				}
#else
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (local_sino > 0.f) {
						denominator(local_ele, ax, local_ind, im_dim, d_OSEM);
					}
					local_ind += d_NN;
				}
#endif
			}
			else if (RHS) {
				local_ele *= *temp;
#ifdef MBSREM
				if (d_alku == 0 && (MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
					if (local_ele < *minimi && local_ele > 0.f)
						*minimi = local_ele;
					d_E[idx] += local_ele;
				}
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (d_alku == 0) {
						if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
#ifdef ATOMIC
							atomicAdd(&d_co[local_ind], __float2ull_rn(*axCOSEM * local_ele * TH));
#else
							atomicAdd(&d_co[local_ind], *axCOSEM * local_ele);
#endif
						if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
#ifdef ATOMIC
							atomicAdd(&d_aco[local_ind], __float2ull_rn(*axCOSEM * local_ele * TH));
#else
							atomicAdd(&d_aco[local_ind], *axCOSEM * local_ele);
#endif
						if (MBSREM_prepass == 1)
#ifdef ATOMIC
							atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
							atomicAdd(&Summ[local_ind], local_ele);
#endif
					}
					else
						*ax += (local_ele * d_OSEM[local_ind]);
					local_ind += d_NN;
				}
#else
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
#ifdef ATOMIC
						atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
						atomicAdd(&Summ[local_ind], local_ele);
#endif
					rhs(MethodList, local_ele, ax, local_ind, im_dim, d_rhs_OSEM);
					local_ind += d_NN;
				}
#endif
			}
			else {
				local_ele *= *temp;
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
#ifdef ATOMIC
					atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
					atomicAdd(&Summ[local_ind], local_ele);
#endif
					local_ind += d_NN;
				}
			}
		}
	}
	for (unsigned int zz = z_loop + 1u; zz < d_N4; zz++) {
		for (int uu = (int)(apu); uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D_per(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
#ifdef VOL
			if (local_ele >= bmax)
				break;
			if (local_ele < bmin)
				local_ele = Vmax;
			else
				local_ele = V[__float2uint_rn((local_ele - bmin) * CC)];
#else
			if (local_ele <= THR)
				break;
#endif
			unsigned int local_ind = uu * d_N + zz * d_N1x;
			if (FPbool) {
				*temp += (local_ele * d_N2);
#ifdef MBSREM
				if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
					*axCOSEM += (local_ele * d_OSEM[local_ind]);
				}
#else
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (local_sino > 0.f) {
						denominator(local_ele, ax, local_ind, im_dim, d_OSEM);
					}
					local_ind += d_NN;
				}
#endif
			}
			else if (RHS) {
				local_ele *= *temp;
#ifdef MBSREM
				if (d_alku == 0 && (MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
					if (local_ele < *minimi && local_ele > 0.f)
						*minimi = local_ele;
					d_E[idx] += local_ele;
				}
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (d_alku == 0) {
						if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
#ifdef ATOMIC
							atomicAdd(&d_co[local_ind], __float2ull_rn(*axCOSEM * local_ele * TH));
#else
							atomicAdd(&d_co[local_ind], *axCOSEM * local_ele);
#endif
						if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
#ifdef ATOMIC
							atomicAdd(&d_aco[local_ind], __float2ull_rn(*axCOSEM * local_ele * TH));
#else
							atomicAdd(&d_aco[local_ind], *axCOSEM * local_ele);
#endif
						if (MBSREM_prepass == 1)
#ifdef ATOMIC
							atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
							atomicAdd(&Summ[local_ind], local_ele);
#endif
					}
					else
						*ax += (local_ele * d_OSEM[local_ind]);
					local_ind += d_NN;
				}
#else
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
#ifdef ATOMIC
						atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
						atomicAdd(&Summ[local_ind], local_ele);
#endif
					rhs(MethodList, local_ele, ax, local_ind, im_dim, d_rhs_OSEM);
					local_ind += d_NN;
				}
#endif
			}
			else {
				local_ele *= *temp;
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
#ifdef ATOMIC
					atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
					atomicAdd(&Summ[local_ind], local_ele);
#endif
					local_ind += d_NN;
				}
			}
		}
		for (unsigned int uu = apu + 1; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D_per(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
#ifdef VOL
			if (local_ele >= bmax)
				break;
			if (local_ele < bmin)
				local_ele = Vmax;
			else
				local_ele = V[__float2uint_rn((local_ele - bmin) * CC)];
#else
			if (local_ele <= THR)
				break;
#endif
			unsigned int local_ind = uu * d_N + zz * d_N1x;
			if (FPbool) {
				*temp += (local_ele * d_N2);
#ifdef MBSREM
				if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
					*axCOSEM += (local_ele * d_OSEM[local_ind]);
				}
#else
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (local_sino > 0.f) {
						denominator(local_ele, ax, local_ind, im_dim, d_OSEM);
					}
					local_ind += d_NN;
				}
#endif
			}
			else if (RHS) {
				local_ele *= *temp;
#ifdef MBSREM
				if (d_alku == 0 && (MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
					if (local_ele < *minimi && local_ele > 0.f)
						*minimi = local_ele;
					d_E[idx] += local_ele;
				}
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (d_alku == 0) {
						if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
#ifdef ATOMIC
							atomicAdd(&d_co[local_ind], __float2ull_rn(*axCOSEM * local_ele * TH));
#else
							atomicAdd(&d_co[local_ind], *axCOSEM * local_ele);
#endif
						if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
#ifdef ATOMIC
							atomicAdd(&d_aco[local_ind], __float2ull_rn(*axCOSEM * local_ele * TH));
#else
							atomicAdd(&d_aco[local_ind], *axCOSEM * local_ele);
#endif
						if (MBSREM_prepass == 1)
#ifdef ATOMIC
							atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
							atomicAdd(&Summ[local_ind], local_ele);
#endif
					}
					else
						*ax += (local_ele * d_OSEM[local_ind]);
					local_ind += d_NN;
				}
#else
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
#ifdef ATOMIC
						atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
						atomicAdd(&Summ[local_ind], local_ele);
#endif
					rhs(MethodList, local_ele, ax, local_ind, im_dim, d_rhs_OSEM);
					local_ind += d_NN;
				}
#endif
			}
			else {
				local_ele *= *temp;
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
#ifdef ATOMIC
					atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
					atomicAdd(&Summ[local_ind], local_ele);
#endif
					local_ind += d_NN;
				}
			}
		}
	}
#ifdef MBSREM
	if (!RHS) {
		*temp = 1.f / *temp;
#ifdef ATN
		* temp *= expf(jelppi);
#endif
#ifdef NORM
		* temp *= d_norm[idx];
#endif
		* temp *= global_factor;
	}
	else {
		if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1)
			d_Amin[idx] = *minimi;
		if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u) {
			if (d_randoms == 1u)
				*ax += d_sc_ra[idx];
			d_ACOSEM_lhs[idx] = *ax;
		}
	}
#else
	if (FPbool) {
		*temp = 1.f / *temp;
#ifdef ATN
		* temp *= expf(jelppi);
#endif
#ifdef NORM
		* temp *= d_norm;
#endif
		* temp *= global_factor;
	}
#endif
}
