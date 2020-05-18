/**************************************************************************
* Special functions for the 2.5D orthogonal distance-based ray tracer
* (CUDA).
*
* Copyright(C) 2020 Ville-Veikko Wettenhovi
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
#pragma once
#include "general_orth_cuda_functions.cuh"

// Denominator (forward projection), multi-GPU version, orthogonal distance based ray tracer
//__device__ void orth_distance_multi(const int tempi, const unsigned int d_N0, const float y_diff, const float x_diff, const float z_diff,
//	const float* y_center, const float* x_center, const float* z_center, const float kerroin, float* temp, const unsigned int d_Nxy,
//	const float local_sino, const unsigned int d_N1, const unsigned int d_N2, const unsigned int d_N3, const int tempj, const int tempk, const float xs,
//	const float ys, const float zs, const char start, int ju, int ku, unsigned char xyz, const unsigned char no_norm, const bool RHS, const bool SUMMA,
//	const float* d_OSEM, float* Summ, const float bmin, const float bmax, const float Vmax,
//	const float* V, float* d_store_elements,
//#ifdef MBSREM
//	unsigned int* d_store_indices, unsigned int* ind, float* ax, const RecMethodsOpenCL MethodListOpenCL, float* d_E, const unsigned int idx,
//	CAST* d_co, CAST* d_aco, float* minimi, const unsigned char MBSREM_prepass, float* axCOSEM, const unsigned int d_alku
//#else
//	unsigned int* d_store_indices, unsigned int* ind, CAST* d_rhs_OSEM, const unsigned int im_dim, const unsigned char* MethodList, float* ax
//#endif
//) {
//	if (start != 0) {
//		bool breikki1 = false;
//		bool breikki2 = false;
//		if (start == 10) {
//			if (xyz == 1u) {
//				for (int uu = tempi; uu >= 0; uu--) {
//					computeOrthVoxelDecreasing(xs, zs, ys, x_diff, z_diff, y_diff, kerroin, x_center[uu], z_center[tempk], y_center, temp, local_sino, ax,
//						&breikki1, tempj, uu, tempk, no_norm, RHS, SUMMA, d_N2, d_Nxy, d_N3, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
//#ifdef MBSREM
//						d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
//#else
//						d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
//#endif
//					computeOrthVoxelIncreasing(xs, zs, ys, x_diff, z_diff, y_diff, kerroin, x_center[uu], z_center[tempk], y_center, temp, local_sino, ax,
//						&breikki2, tempj + 1, uu, tempk, no_norm, RHS, SUMMA, d_N2, d_Nxy, d_N3, d_N1, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
//#ifdef MBSREM
//						d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
//#else
//						d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
//#endif
//					if (breikki1 && breikki2)
//						break;
//				}
//				breikki1 = false;
//				breikki2 = false;
//				for (int uu = tempi + 1; uu < d_N0; uu++) {
//					computeOrthVoxelDecreasing(xs, zs, ys, x_diff, z_diff, y_diff, kerroin, x_center[uu], z_center[tempk], y_center, temp, local_sino, ax,
//						&breikki1, tempj, uu, tempk, no_norm, RHS, SUMMA, d_N2, d_Nxy, d_N3, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
//#ifdef MBSREM
//						d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
//#else
//						d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
//#endif
//					computeOrthVoxelIncreasing(xs, zs, ys, x_diff, z_diff, y_diff, kerroin, x_center[uu], z_center[tempk], y_center, temp, local_sino, ax,
//						&breikki2, tempj + 1, uu, tempk, no_norm, RHS, SUMMA, d_N2, d_Nxy, d_N3, d_N1, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
//#ifdef MBSREM
//						d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
//#else
//						d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
//#endif
//					if (breikki1 && breikki2)
//						break;
//				}
//			}
//			else {
//				for (int yy = tempj; yy >= 0; yy--) {
//					computeOrthVoxelDecreasing(zs, ys, xs, z_diff, y_diff, x_diff, kerroin, z_center[tempk], y_center[yy], x_center, temp, local_sino, ax,
//						&breikki1, tempi, tempk, yy, no_norm, RHS, SUMMA, d_Nxy, d_N3, d_N2, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
//#ifdef MBSREM
//						d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
//#else
//						d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
//#endif
//					computeOrthVoxelIncreasing(zs, ys, xs, z_diff, y_diff, x_diff, kerroin, z_center[tempk], y_center[yy], x_center, temp, local_sino, ax,
//						&breikki2, tempi + 1, tempk, yy, no_norm, RHS, SUMMA, d_Nxy, d_N3, d_N2, d_N0, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
//#ifdef MBSREM
//						d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
//#else
//						d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
//#endif
//					if (breikki1 && breikki2)
//						break;
//				}
//				breikki1 = false;
//				breikki2 = false;
//				for (int yy = tempj + 1; yy < d_N1; yy++) {
//					computeOrthVoxelDecreasing(zs, ys, xs, z_diff, y_diff, x_diff, kerroin, z_center[tempk], y_center[yy], x_center, temp, local_sino, ax,
//						&breikki1, tempi, tempk, yy, no_norm, RHS, SUMMA, d_Nxy, d_N3, d_N2, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
//#ifdef MBSREM
//						d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
//#else
//						d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
//#endif
//					computeOrthVoxelIncreasing(zs, ys, xs, z_diff, y_diff, x_diff, kerroin, z_center[tempk], y_center[yy], x_center, temp, local_sino, ax,
//						&breikki2, tempi + 1, tempk, yy, no_norm, RHS, SUMMA, d_Nxy, d_N3, d_N2, d_N0, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
//#ifdef MBSREM
//						d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
//#else
//						d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
//#endif
//					if (breikki1 && breikki2)
//						break;
//				}
//			}
//		}
//		else {
//			if (start < 0)
//				ju *= -1;
//			if (ju > 0) {
//				if (xyz == 1u) {
//					for (int uu = tempi; uu >= 0; uu--) {
//						computeOrthVoxelDecreasing(xs, zs, ys, x_diff, z_diff, y_diff, kerroin, x_center[uu], z_center[tempk], y_center, temp, local_sino, ax,
//							&breikki1, tempj, uu, tempk, no_norm, RHS, SUMMA, d_N2, d_Nxy, d_N3, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
//#ifdef MBSREM
//							d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
//#else
//							d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
//#endif
//						if (breikki1)
//							break;
//					}
//					for (int uu = tempi + 1; uu < d_N0; uu++) {
//						computeOrthVoxelDecreasing(xs, zs, ys, x_diff, z_diff, y_diff, kerroin, x_center[uu], z_center[tempk], y_center, temp, local_sino, ax,
//							&breikki2, tempj, uu, tempk, no_norm, RHS, SUMMA, d_N2, d_Nxy, d_N3, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
//#ifdef MBSREM
//							d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
//#else
//							d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
//#endif
//						if (breikki2)
//							break;
//					}
//				}
//				else {
//					for (int yy = tempj; yy >= 0; yy--) {
//						computeOrthVoxelDecreasing(zs, ys, xs, z_diff, y_diff, x_diff, kerroin, z_center[tempk], y_center[yy], x_center, temp, local_sino, ax,
//							&breikki1, tempi, tempk, tempj, no_norm, RHS, SUMMA, d_Nxy, d_N3, d_N2, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
//#ifdef MBSREM
//							d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
//#else
//							d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
//#endif
//						computeOrthVoxelIncreasing(zs, ys, xs, z_diff, y_diff, x_diff, kerroin, z_center[tempk], y_center[yy], x_center, temp, local_sino, ax,
//							&breikki2, tempi + 1, tempk, yy, no_norm, RHS, SUMMA, d_Nxy, d_N3, d_N2, d_N0, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
//#ifdef MBSREM
//							d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
//#else
//							d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
//#endif
//						if (breikki1 && breikki2)
//							break;
//					}
//				}
//			}
//			else {
//				if (xyz == 1u) {
//					for (int uu = tempi; uu >= 0; uu--) {
//						computeOrthVoxelIncreasing(xs, zs, ys, x_diff, z_diff, y_diff, kerroin, x_center[uu], z_center[tempk], y_center, temp, local_sino, ax,
//							&breikki1, tempj, uu, tempk, no_norm, RHS, SUMMA, d_N2, d_Nxy, d_N3, d_N1, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
//#ifdef MBSREM
//							d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
//#else
//							d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
//#endif
//						if (breikki1)
//							break;
//					}
//					for (int uu = tempi + 1; uu < d_N0; uu++) {
//						computeOrthVoxelIncreasing(xs, zs, ys, x_diff, z_diff, y_diff, kerroin, x_center[uu], z_center[tempk], y_center, temp, local_sino, ax,
//							&breikki2, tempj, uu, tempk, no_norm, RHS, SUMMA, d_N2, d_Nxy, d_N3, d_N1, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
//#ifdef MBSREM
//							d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
//#else
//							d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
//#endif
//						if (breikki2)
//							break;
//					}
//				}
//				else {
//					for (int yy = tempj; yy < d_N1; yy++) {
//						computeOrthVoxelDecreasing(zs, ys, xs, z_diff, y_diff, x_diff, kerroin, z_center[tempk], y_center[yy], x_center, temp, local_sino, ax,
//							&breikki1, tempi, tempk, tempj, no_norm, RHS, SUMMA, d_Nxy, d_N3, d_N2, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
//#ifdef MBSREM
//							d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
//#else
//							d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
//#endif
//						computeOrthVoxelIncreasing(zs, ys, xs, z_diff, y_diff, x_diff, kerroin, z_center[tempk], y_center[yy], x_center, temp, local_sino, ax,
//							&breikki2, tempi + 1, tempk, yy, no_norm, RHS, SUMMA, d_Nxy, d_N3, d_N2, d_N0, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
//#ifdef MBSREM
//							d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
//#else
//							d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
//#endif
//						if (breikki1 && breikki2)
//							break;
//					}
//				}
//			}
//		}
//	}
//	else {
//		bool breikki1 = false;
//		computeOrthVoxelDecreasing(zs, ys, xs, z_diff, y_diff, x_diff, kerroin, z_center[tempk], y_center[tempj], x_center, temp, local_sino, ax,
//			&breikki1, tempi, tempk, tempj, no_norm, RHS, SUMMA, d_Nxy, d_N3, d_N2, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
//#ifdef MBSREM
//			d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
//#else
//			d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
//#endif
//		computeOrthVoxelIncreasing(zs, ys, xs, z_diff, y_diff, x_diff, kerroin, z_center[tempk], y_center[tempj], x_center, temp, local_sino, ax,
//			&breikki1, tempi + 1, tempk, tempj, no_norm, RHS, SUMMA, d_Nxy, d_N3, d_N2, d_N0, Summ, d_OSEM, bmin, bmax, Vmax, V, d_store_elements,
//#ifdef MBSREM
//			d_store_indices, ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
//#else
//			d_store_indices, ind, d_rhs_OSEM, im_dim, MethodList);
//#endif
//	}
//}

// Denominator (forward projection), orthogonal distance based ray tracer, multi-GPU
__device__ void orth_distance_perpendicular_multi(const float* center1, const float center2, const float* z_center, const float kerroin,
	float* temp, const unsigned int d_attenuation_correction, const unsigned int d_normalization, float* ax, const float d_b, const float d,
	const float d_d1, const unsigned int d_N1, const unsigned int d_N2, const unsigned int z_loop, const float* d_atten, const float d_norm, const float local_sino,
	const float xs, const float ys, const float zs, const float xl, const float yl, const float zl, const unsigned int d_N, const unsigned int d_NN,
	const float* d_OSEM, const unsigned char no_norm, CAST* Summ, const bool FP_bool, const bool RHS, const float global_factor,
#ifdef MBSREM
	const RecMethodsOpenCL MethodListOpenCL, const unsigned int d_alku, float* axCOSEM, 
	float* d_E, CAST* d_co, CAST* d_aco, float* minimi, const unsigned char MBSREM_prepass,
	const float* d_sc_ra, const unsigned int d_randoms, float* d_Amin, float* d_ACOSEM_lhs, const unsigned int idx
#else
	CAST* d_rhs_OSEM, const unsigned int im_dim, const unsigned char* MethodList
#endif
) {
	const unsigned int zz = z_loop * d_N2 * d_N1;
	const int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	float jelppi = 0.f;
	for (int uu = apu; uu >= 0; uu--) {
		float local_ele = compute_element_orth_3D_per(xs, ys, zs, xl, yl, zl, kerroin, center1[uu], center2, z_center[z_loop]);
		if (local_ele <= THR)
			break;
		unsigned int local_ind = (unsigned int)(uu) * d_N + zz;
		if (FP_bool) {
			*temp += (local_ele * d_N2);
#if defined(ATN) || defined(MBSREM)
			for (unsigned int kk = 0u; kk < d_N2; kk++) {
				if (d_attenuation_correction && uu == apu)
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
						atomicAdd(&d_aco[local_ind], *axCOSEM * local_ele));
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
	for (unsigned int uu = (unsigned int)(apu) + 1u; uu < d_N1; uu++) {
		float local_ele = compute_element_orth_3D_per(xs, ys, zs, xl, yl, zl, kerroin, center1[uu], center2, z_center[z_loop]);
		if (local_ele <= THR)
			break;
		unsigned int local_ind = uu * d_N + zz;
		if (FP_bool) {
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
	if (FP_bool) {
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