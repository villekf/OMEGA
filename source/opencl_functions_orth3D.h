
/*******************************************************************************************************************************************
* Special functions for the 3D orthogonal distance-based ray tracer.
*
* Copyright (C) 2022 Ville-Veikko Wettenhovi
*
* This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************************************************************************/
//#pragma once
//#include "general_orth_opencl_functions.h"
#ifdef _MSC_VER
#include "general_orth_opencl_functions.h"
#endif

bool orthogonalHelper3D(const int tempi, const int uu, const uint d_N2, const uint d_N3, const uint d_Nxy, const int zz, const float s2, const float l3, const float l1, const float l2,
	const float diff1, const float diffZ, const float kerroin, const float center2, const float bmin, const float bmax, const float Vmax,
	__constant float* V, const bool RHS, const bool XY,  const bool SUMMA, uint* indeksi, 
	__global CAST* Summ,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
	__read_only image3d_t d_OSEM, 
#endif
	float* temp, const float local_sino, __private float* d_store_elements, const int3 i, 
#ifdef MBSREM
	__private uint* d_store_indices, float* ax, const RecMethodsOpenCL MethodListOpenCL, __global float* d_E, const size_t idx,
	__global CAST* d_co, __global CAST* d_aco, float* minimi, const uchar MBSREM_prepass, float* axCOSEM, const uint d_alku
#else
	__private uint* d_store_indices, __global CAST* d_rhs_OSEM, const uint im_dim, __constant uchar* MethodList, float* ax, const bool no_norm
#endif
) {
	int4 ind;
	if (!RHS) {
		if (XY)
			ind = (int4)(tempi, uu, zz, 0);
		else
			ind = (int4)(uu, tempi, zz, 0);
	}
	float local_ele = compute_element_orth_3D(s2, l3, l1, l2, diff1, diffZ, kerroin, center2);
	//const float x0 = center2 - s2;
	//const float y1 = mad(diffZ, x0, -l2);
	//const float z1 = mad(-diff1, x0, l3);
	//if ((i.x == 89 && i.y == 89 && i.z == 20)) {
	//if (i.x == 30 && i.y == 86 && i.z == 543) {
	////if (i.x == 93 && i.y == 55 && i.z == 167 && tempi == 36 && uu > 88 && uu < 92 && zz > 40 && zz < 47) {
	//		//printf("s2 = %f\n", s2);
	//		//printf("l3 = %f\n", l3);
	//		//printf("l1 = %f\n", l1);
	//		//printf("l2 = %f\n", l2);
	//		//printf("diff1 = %f\n", diff1);
	//		//printf("diffZ = %f\n", diffZ);
	//		//printf("kerroin = %f\n", kerroin);
	//		//printf("center2 = %f\n", center2);
	//		//printf("x0 = %f\n", x0);
	//		//printf("y1 = %f\n", y1);
	//		//printf("z1 = %f\n", z1);
	//		printf("local_ele = %f\n", local_ele);
	//		printf("uu = %d\n", uu);
	//		printf("zz = %d\n", zz);
	//		//printf("tempi = %d\n", tempi);
	//}
#ifdef VOL
	if (local_ele >= bmax) {
		return true;
	}
	if (local_ele < bmin)
		local_ele = Vmax;
	else
		local_ele = V[convert_uint_rte((local_ele - bmin) * CC)];
#else
	if (local_ele <= THR) {
		return true;
	}
#endif
	uint local_ind = 0;
	if (RHS)
		local_ind = compute_ind_orth_3D(convert_uint(tempi), uu * d_N3, (zz), d_N2, d_Nxy);
#ifdef AF
#ifdef MBSREM
	computeIndicesOrth_cosem(RHS, local_ele, temp, ax, Summ, local_sino, d_OSEM, local_ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx, ind);
#else
	computeIndicesOrth_af(RHS, SUMMA, local_ele, temp, ax, no_norm, Summ, d_rhs_OSEM, local_sino, d_OSEM, local_ind, im_dim, MethodList, ind);
#endif
#else
	computeIndicesOrth(RHS, SUMMA, local_ele, temp, ax, no_norm, Summ, d_rhs_OSEM, local_sino,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
		d_OSEM, 
#endif
		local_ind, ind);
#endif
//	if (i.x == 93 && i.y == 55 && i.z == 167 && tempi > 30 && tempi < 80) {
//		printf("ax = %f", *ax);
//		printf("  uu = %d", uu);
//		printf("  zz = %d", zz);
//		printf("  tempi = %d\n", tempi);
//}
#ifdef DEC
	local_ind = compute_ind_orth_3D(convert_uint(tempi), uu * d_N3, (zz), d_N2, d_Nxy);
	d_store_elements[*indeksi] = local_ele;
	d_store_indices[*indeksi] = local_ind;
	*indeksi = *indeksi + 1u;
#endif
	return false;
}
 
// Denominator (forward projection), multi-GPU version, orthogonal distance based ray tracer
//void orth_distance_multi_3D(const int tempi, const uint d_N0, const uint d_N4, const float3 diff, 
//	__constant float* y_center, __constant float* x_center, __constant float* z_center, float* temp, const uint d_N2, const int tempj, const int tempk,
//	const float local_sino, const __global float* d_OSEM, const float3 s, const uint d_Nxy, const float kerroin,
//	const bool no_norm, const bool RHS, const bool SUMMA, __global CAST* Summ, const uint d_N1, const uint d_N3, const int start, const int iu, 
//	const int ju, const int loppu, const float bmin, const float bmax, const float Vmax, __constant float* V, __private float* d_store_elements,
void orthDistance3D(const int tempi, const float diff1, const float diff2, const float diffZ, 
	const float center1, __constant float* center2, __constant float* centerZ, float* temp, int temp2, const int tempk, const float local_sino,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
	__read_only image3d_t d_OSEM, 
#endif
	const float s1, const float s2, const float sZ, const uint d_Nxy, const float kerroin,
	const bool no_norm, const bool RHS, const bool SUMMA, __global CAST* Summ, const uint d_N1, const uint d_N2, const uint d_N3, const uint d_Nz, const float bmin, 
	const float bmax, const float Vmax, __constant float* V, __private float* d_store_elements, const bool XY, const int3 i, 
#ifdef MBSREM
	__private uint* d_store_indices, uint* ind, float* ax, const RecMethodsOpenCL MethodListOpenCL, __global float* d_E, const size_t idx,
	__global CAST* d_co, __global CAST* d_aco, float* minimi, const uchar MBSREM_prepass, float* axCOSEM, const uint d_alku
#else
	__private uint* d_store_indices, uint* ind, __global CAST* d_rhs_OSEM, const uint im_dim, __constant uchar* MethodList, float* ax
#endif
) {
	uint local_ind = 0;
	bool breikki = false;
	//int tempj = temp2;
	//if (tempjapu == temp2)
	//	uy = 0;
	// y0
	const float v0 = center1 - s1;
	// xl * y0
	const float l3 = diff1 * v0;
	// zl * y0
	const float apu1 = diffZ * v0;
	//int maksimiZ, minimiZ, minimiXY, maksimiXY;
	//if (uz > 0) {
	//	//maksimiZ = min(tempk + NSTEPS, convert_int(d_Nz) - 1);
	//	maksimiZ = convert_int(d_Nz);
	//	//minimiZ = max(tempkapu - NSTEPS, 0);
	//	minimiZ = 0;
	//}
	//else {
		//maksimiZ = min(tempkapu + NSTEPS, convert_int(d_Nz) - 1);
		//maksimiZ = convert_int(d_Nz);
		const int maksimiZ = convert_int(d_Nz);
		//minimiZ = max(tempk - NSTEPS, 0);
		//minimiZ = 0;
		const int minimiZ = 0;
	//}
	//if (uy > 0) {
		//maksimiXY = min(temp2 + NSTEPS, convert_int(d_N1) - 1);
		//maksimiXY = convert_int(d_N1);
		const int maksimiXY = convert_int(d_N1);
		//minimiXY = max(tempjapu - NSTEPS, 0);
		//minimiXY = 0;
		const int minimiXY = 0;
	//}
	//else {
	//	//maksimiXY = min(tempjapu + NSTEPS, convert_int(d_N1) - 1);
	//	maksimiXY = convert_int(d_N1);
	//	//minimiXY = max(temp2 - NSTEPS, 0);
	//	minimiXY = 0;
	//}
	int uu1 = 0, uu2 = 0;
	//if (i.x == 126 && i.y == 35 && i.z == 697 && tempi < 30) {
	////	printf("s.x = %f\n", s.x);
	////	printf("s.y = %f\n", s.y);
	//	//printf("s1 = %f\n", s1);
	//	printf("center1 = %f\n", center1);
	//	//printf("diff2 = %f\n", diff2);
	//	//printf("v0 = %f\n", v0);
	//	//printf("maksimiZ = %d\n", maksimiZ);
	//	//printf("minimiZ = %d\n", minimiZ);
	//	//printf("maksimiXY = %d\n", maksimiXY);
	//	//printf("minimiXY = %d\n", minimiXY);
	//	//printf("tempi = %d\n", tempi);
	//	printf("tempk = %d\n", tempk);
	//	printf("temp2 = %d\n", temp2);
	//	//printf("tempkapu = %d\n", tempkapu);
	//	//printf("tempjapu = %d\n", tempjapu);
	//}
#ifdef CRYSTZ
	for (int zz = tempk; zz < maksimiZ; zz++) {
		//int hh1 = 2;
		//int hh2 = 2;
#else
	zz = tempk;
#endif
		// z0
		const float z0 = centerZ[zz] - sZ;
		// x1 = yl * z0 - zl * y0
		const float l1 = diff2 * z0 - apu1;
		// xl * z0
		const float l2 = diff1 * z0;
#ifdef CRYSTXY
		for (uu1 = temp2; uu1 < maksimiXY; uu1++) {
#else
		uu1 = temp2;
#endif
		//if (i.x == 93 && i.y == 55 && i.z == 167 && tempi == 36 && uu1 > 88 && uu1 < 92 && zz > 40 && zz < 47) {
		//	printf("z0 = %f\n", z0);
		//	printf("centerZ[zz] = %f\n", centerZ[zz]);
		//	printf("sZ = %f\n", sZ);
		//	printf("diff2 = %f\n", diff2);
		//	printf("diff1 = %f\n", diff1);
		//	printf("l1 = %f\n", l1);
		//	printf("l2 = %f\n", l2);
		//	printf("apu1 = %f\n", apu1);
		//	printf("center1 = %f\n", center1);
		//	printf("center2[uu1] = %f\n", center2[uu1]);
		//	printf("uu1 = %d\n", uu1);
		//	printf("zz = %d\n", zz);
		//}
			breikki = orthogonalHelper3D(tempi, uu1, d_N2, d_N3, d_Nxy, zz, s2, l3, l1, l2, diff2, diffZ, kerroin, center2[uu1], bmin, bmax, Vmax, V, RHS,
				XY, SUMMA, ind, Summ,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
				d_OSEM,
#endif
				temp, local_sino, d_store_elements, i,
#ifdef MBSREM
				d_store_indices, ax, MethodListOpenCL, d_E, idx, d_co, d_aco, minimi, MBSREM_prepass, axCOSEM, d_alku
#else
				d_store_indices, d_rhs_OSEM, im_dim, MethodList, ax, no_norm
#endif
			);
#ifdef CRYSTXY
			//if (i.x == 93 && i.y == 55 && i.z == 167) {
			//}
			//else {
				//if (breikki && uu1 > temp2) {
				if (breikki) {
					//hh1++;
					break;
				}
			//}
		}
#endif
#ifdef CRYSTXY
		for (uu2 = temp2 - 1; uu2 >= minimiXY; uu2--) {
		//	if (i.x == 93 && i.y == 55 && i.z == 167 && tempi == 36 && uu1 > 88 && uu1 < 92 && zz > 40 && zz < 47) {
		//		printf("z0 = %f\n", z0);
		//		printf("centerZ[zz] = %f\n", centerZ[zz]);
		//		printf("sZ = %f\n", sZ);
		//		printf("diff2 = %f\n", diff2);
		//		printf("diff1 = %f\n", diff1);
		//		printf("l1 = %f\n", l1);
		//		printf("l2 = %f\n", l2);
		//		printf("apu1 = %f\n", apu1);
		//		printf("center1 = %f\n", center1);
		//		printf("center2[uu2] = %f\n", center2[uu2]);
		//		printf("uu2 = %d\n", uu2);
		//		printf("zz = %d\n", zz);
		//}
			breikki = orthogonalHelper3D(tempi, uu2, d_N2, d_N3, d_Nxy, zz, s2, l3, l1, l2, diff2, diffZ, kerroin, center2[uu2], bmin, bmax, Vmax, V, RHS,
				XY, SUMMA, ind, Summ,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
				d_OSEM,
#endif
				temp, local_sino, d_store_elements, i,
#ifdef MBSREM
				d_store_indices, ax, MethodListOpenCL, d_E, idx, d_co, d_aco, minimi, MBSREM_prepass, axCOSEM, d_alku
#else
				d_store_indices, d_rhs_OSEM, im_dim, MethodList, ax, no_norm
#endif
			);
			//if (i.x == 93 && i.y == 55 && i.z == 167) {
			//}
			//else {
				//if (breikki && uu1 > temp2) {
				if (breikki) {
					//hh1++;
					break;
				}
			//}
		}
#else
	uu2 = temp2 - 1;
#endif
#ifdef CRYSTZ
	if (uu1 == temp2 && uu2 == temp2 - 1 && breikki)
		break;
	//temp2 -= uy;
	}
	//temp2 = tempj;
	for (int zz = tempk - 1; zz >= minimiZ; zz--) {
		const float z0 = centerZ[zz] - sZ;
		const float l1 = diff2 * z0 - apu1;
		const float l2 = diff1 * z0;
		//int hh1 = 2;
		//int hh2 = 2;
#ifdef CRYSTXY
		for (uu1 = temp2; uu1 < maksimiXY; uu1++) {
#else
		uu1 = temp2;
#endif
		//if (i.x == 93 && i.y == 55 && i.z == 167 && tempi == 36 && uu1 > 88 && uu1 < 92 && zz > 40 && zz < 47) {
		//	printf("z0 = %f\n", z0);
		//	printf("centerZ[zz] = %f\n", centerZ[zz]);
		//	printf("sZ = %f\n", sZ);
		//	printf("diff2 = %f\n", diff2);
		//	printf("diff1 = %f\n", diff1);
		//	printf("l1 = %f\n", l1);
		//	printf("l2 = %f\n", l2);
		//	printf("apu1 = %f\n", apu1);
		//	printf("center1 = %f\n", center1);
		//	printf("center2[uu1] = %f\n", center2[uu1]);
		//	printf("uu1 = %d\n", uu1);
		//	printf("zz = %d\n", zz);
		//}
			breikki = orthogonalHelper3D(tempi, uu1, d_N2, d_N3, d_Nxy, zz, s2, l3, l1, l2, diff2, diffZ, kerroin, center2[uu1], bmin, bmax, Vmax, V, RHS,
				XY, SUMMA, ind, Summ,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
				d_OSEM,
#endif
				temp, local_sino, d_store_elements, i,
#ifdef MBSREM
				d_store_indices, ax, MethodListOpenCL, d_E, idx, d_co, d_aco, minimi, MBSREM_prepass, axCOSEM, d_alku
#else
				d_store_indices, d_rhs_OSEM, im_dim, MethodList, ax, no_norm
#endif
			);
#ifdef CRYSTXY
			//if (i.x == 93 && i.y == 55 && i.z == 167) {
			//}
			//else {
				//if (breikki && uu1 > temp2) {
				if (breikki) {
					//hh1++;
					break;
				}
			//}
		}
#endif
#ifdef CRYSTXY
		for (uu2 = temp2 - 1; uu2 >= minimiXY; uu2--) {
			//if (i.x == 93 && i.y == 55 && i.z == 167 && tempi == 36 && uu1 > 88 && uu1 < 92 && zz > 40 && zz < 47) {
			//	printf("z0 = %f\n", z0);
			//	printf("centerZ[zz] = %f\n", centerZ[zz]);
			//	printf("sZ = %f\n", sZ);
			//	printf("diff2 = %f\n", diff2);
			//	printf("diff1 = %f\n", diff1);
			//	printf("l1 = %f\n", l1);
			//	printf("l2 = %f\n", l2);
			//	printf("apu1 = %f\n", apu1);
			//	printf("center1 = %f\n", center1);
			//	printf("center2[uu2] = %f\n", center2[uu2]);
			//	printf("uu2 = %d\n", uu2);
			//	printf("zz = %d\n", zz);
			//}
			breikki = orthogonalHelper3D(tempi, uu2, d_N2, d_N3, d_Nxy, zz, s2, l3, l1, l2, diff2, diffZ, kerroin, center2[uu2], bmin, bmax, Vmax, V, RHS,
				XY, SUMMA, ind, Summ,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
				d_OSEM,
#endif
				temp, local_sino, d_store_elements, i,
#ifdef MBSREM
				d_store_indices, ax, MethodListOpenCL, d_E, idx, d_co, d_aco, minimi, MBSREM_prepass, axCOSEM, d_alku
#else
				d_store_indices, d_rhs_OSEM, im_dim, MethodList, ax, no_norm
#endif
			);
			//if (i.x == 93 && i.y == 55 && i.z == 167) {
			//}
			//else {
				//if (breikki && uu1 > temp2) {
				if (breikki) {
					//hh1++;
					break;
				}
			//}
		}
#else
	uu2 = temp2 - 1;
#endif
		if (uu1 == temp2 && uu2 == temp2 - 1 && breikki)
			break;
		//temp2 -= uy;
	}
#endif
	//if (i.x == 30 && i.y == 86 && i.z == 543) {
	//	printf("ax = %f\n", *ax);
	//}
}

void orthHelperPerpendicular(float* temp, float* ax, const float d_d1, const uint d_N2, const int z_loop, const int zz, float* jelppi,  
	const float local_sino, float local_ele, const uint apu, const int uu, 
	const uint d_NN,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
__read_only image3d_t d_OSEM, 
#endif
const int4 local_ind, uint ind, const uchar no_norm, __global CAST* Summ, const bool FPbool, const bool RHS, const float global_factor, 
#ifdef MBSREM
	const RecMethodsOpenCL MethodListOpenCL, const uint d_alku, float* axCOSEM, __global float* d_E, __global CAST* d_co, 
	__global CAST* d_aco, float* minimi, const uchar MBSREM_prepass, const size_t idx
#else
	__global CAST* d_rhs_OSEM, const uint im_dim, __constant uchar* MethodList
#endif
#if !defined(CT)
	, __read_only image3d_t d_atten
#endif
) {
	if (FPbool) {
#ifndef CT
		* temp += (local_ele * d_N2);
#endif
#if defined(ATN) || defined(FP)
		for (uint kk = 0u; kk < d_N2; kk++) {
#ifdef ATN
			if (zz == convert_int_sat(z_loop) && uu == convert_int_sat(apu))
				*jelppi += (d_d1 * -read_imagef(d_atten, samplerIm, local_ind).w);
				//jelppi += (d_d1 * -d_atten[local_ind]);
#endif
#ifdef MBSREM
			if (local_sino != 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || 
				MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
				*axCOSEM += (local_ele * read_imagef(d_OSEM, samplerIm, local_ind).w);
				//*axCOSEM += (local_ele * d_OSEM[local_ind]);
			}
#elif defined(FP)
#ifdef BP
			if (local_sino != 0.f) {
#endif
#ifdef AF
				denominator(local_ele, ax, local_ind, im_dim, d_OSEM);
#else
				//denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
				denominator_multi(local_ele, ax, read_imagef(d_OSEM, samplerIm, local_ind).w);
#endif
#ifdef BP
			}
#endif
#endif
			if (d_NN == 1)
				local_ind.x++;
			else
				local_ind.y++;
			//local_ind += d_NN;
		}
#endif
	}
	else if (RHS) {
#ifndef CT
		local_ele *= *temp;
#endif
#ifdef MBSREM
		if (d_alku == 0 && (MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
			if (local_ele < *minimi && local_ele > 0.f)
				*minimi = local_ele;
			d_E[idx] += local_ele;
		}
		for (uint kk = 0u; kk < d_N2; kk++) {
			if (d_alku == 0) {
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino != 0.f)
#ifdef ATOMIC
					atom_add(&d_co[ind], convert_long(*axCOSEM * local_ele * TH));
				//atom_add(&d_co[local_ind], convert_long(*axCOSEM * local_ele * TH));
#else
					atomicAdd_g_f(&d_co[ind], *axCOSEM * local_ele);
				//atomicAdd_g_f(&d_co[local_ind], *axCOSEM * local_ele);
#endif
				if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino != 0.f)
#ifdef ATOMIC
					atom_add(&d_aco[ind], convert_long(*axCOSEM * local_ele * TH));
				//atom_add(&d_aco[local_ind], convert_long(*axCOSEM * local_ele * TH));
#else
					atomicAdd_g_f(&d_aco[ind], *axCOSEM * local_ele);
				//atomicAdd_g_f(&d_aco[local_ind], *axCOSEM * local_ele);
#endif
				if (MBSREM_prepass == 1)
#ifdef ATOMIC
					atom_add(&Summ[ind], convert_long(local_ele * TH));
				//atom_add(&Summ[local_ind], convert_long(local_ele * TH));
#else
					atomicAdd_g_f(&Summ[ind], local_ele);
				//atomicAdd_g_f(&Summ[local_ind], local_ele);
#endif
			}
			else
				*ax += (local_ele * read_imagef(d_OSEM, samplerIm, local_ind).w);
				//*ax += (local_ele * d_OSEM[local_ind]);
			ind += d_NN;
			if (d_NN == 1)
				local_ind.x++;
			else
				local_ind.y++;
			//local_ind += d_NN;
		}
#else
		for (uint kk = 0u; kk < d_N2; kk++) {
			if (no_norm == 0u)
#ifdef ATOMIC
				atom_add(&Summ[ind], convert_long(local_ele * TH));
			//atom_add(&Summ[local_ind], convert_long(local_ele * TH));
#else
				atomicAdd_g_f(&Summ[ind], local_ele);
			//atomicAdd_g_f(&Summ[local_ind], local_ele);
#endif
#if defined(FP) && defined(BP)
			if (local_sino != 0.f) {
#endif
#ifdef AF
				rhs(MethodList, local_ele, ax, ind, im_dim, d_rhs_OSEM);
				//rhs(MethodList, local_ele, ax, local_ind, im_dim, d_rhs_OSEM);
#else
#ifdef ATOMIC
				atom_add(&d_rhs_OSEM[ind], convert_long(local_ele * *ax * TH));
				//atom_add(&d_rhs_OSEM[local_ind], convert_long(local_ele * *ax * TH));
#else
				atomicAdd_g_f(&d_rhs_OSEM[ind], (local_ele * *ax));
				//atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * *ax));
#endif
#endif
#if defined(FP) && defined(BP)
			}
#endif
			ind += d_NN;
			//local_ind += d_NN;
		}
#endif
	}
//	else {
//#ifndef CT
//		local_ele *= *temp;
//#endif
//		for (uint kk = 0u; kk < d_N2; kk++) {
//#ifdef ATOMIC
//			atom_add(&Summ[ind], convert_long(local_ele * TH));
//			//atom_add(&Summ[local_ind], convert_long(local_ele * TH));
//#else
//			atomicAdd_g_f(&Summ[ind], local_ele);
//			//atomicAdd_g_f(&Summ[local_ind], local_ele);
//#endif
//			ind += d_NN;
//			//local_ind += d_NN;
//		}
//	}
}


void orthDistancePerpendicularMulti3D(__constant float* center1, const float center2, __constant float* z_center, float* temp, float* ax, 
	const float d_b, const float d, const float d_d1, const uint d_N1, const uint d_N2, const uint d_Nz, const int z_loop, 
	const float d_norm, const float local_sino, const uint d_N, const uint d_NN,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
	__read_only image3d_t d_OSEM, 
#endif
	const float3 s, const float3 diff, const float crystal_size_z, const uint d_Nxy, const uchar no_norm, __global CAST* Summ,
	const bool FPbool, const bool RHS, const float global_factor, const float bmin, const float bmax, const float Vmax, __constant float* V,
#ifdef MBSREM
	const RecMethodsOpenCL MethodListOpenCL, const uint d_alku, float* axCOSEM,
	__global float* d_E, __global CAST* d_co, __global CAST* d_aco, float* minimi, const uchar MBSREM_prepass,
	const __global float* d_sc_ra, __global float* d_Amin, __global float* d_ACOSEM_lhs, const size_t idx
#else
	__global CAST* d_rhs_OSEM, const uint im_dim, __constant uchar* MethodList
#endif
#if !defined(CT) && defined(ATN)
	, __read_only image3d_t d_atten
#endif
) {
	int4 local_ind = { 0, 0, 0, 0 };
	//const uint zz = z_loop * d_N2 * d_N1;
	const int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	float jelppi = 0.f;
	const int maksimiXY = min(apu + NSTEPS, convert_int(d_N1));
	const int minimiXY = max(apu - NSTEPS, 0);
#ifdef CRYSTZ
	const int maksimiZ = min(z_loop + NSTEPS, convert_int(d_Nz));
	const int minimiZ = max(z_loop - NSTEPS, 0);
	//if (z_loop == 28 && s.x > 171.651f && s.x < 171.652f) {
	//	printf("s.x = %f\n", s.x);
	//	printf("s.y = %f\n", s.y);
	//	printf("s.z = %f\n", s.z);
	//	//printf("maksimiZ = %d\n", maksimiZ);
	//	//printf("minimiZ = %d\n", minimiZ);
	//	//printf("maksimiXY = %d\n", maksimiXY);
	//	//printf("minimiXY = %d\n", minimiXY);
	//	//printf("apu = %d\n", apu);
	//}
	for (int zz = z_loop; zz >= minimiZ; zz--) {
#else
	int zz = z_loop;
#endif
#ifdef CRYSTXY
	for (int uu = apu; uu >= minimiXY; uu--) {
#else
	int uu = convert_int_sat(apu);
#endif
	//const float3 p0 = (float3)(center1[uu], center2, z_center[zz]) - s;
	//const float3 p1 = cross(diff, p0);
	float local_ele = computeElementOrth3DPer(s, diff, crystal_size_z, center1[uu], center2, z_center[zz]);
	//if (z_loop == 28 && s.x > 171.651f && s.x < 171.652f) {
	//	//printf("p0.x = %f\n", p0.x);
	//	//printf("p0.y = %f\n", p0.y);
	//	//printf("p0.z = %f\n", p0.z);
	//	//printf("p1.x = %f\n", p1.x);
	//	//printf("p1.y = %f\n", p1.y);
	//	//printf("p1.z = %f\n", p1.z);
	//	printf("local_ele = %f\n", local_ele);
	//	printf("center1[uu] = %f\n", center1[uu]);
	//	printf("center2 = %f\n", center2);
	//	printf("z_center[zz] = %f\n", z_center[zz]);
	//}
#ifdef VOL
	if (local_ele >= bmax)
		break;
	if (local_ele < bmin)
		local_ele = Vmax;
	else
		local_ele = V[convert_uint_rte((local_ele - bmin) * CC)];
#else
	if (local_ele <= THR)
		break;
#endif
	//uint local_ind = uu * d_N + zz * d_N1x;
	uint ind = uu * d_N + zz * d_Nxy;
	if (d_N == 1)
		local_ind.x = uu;
	else
		local_ind.y = uu;
	local_ind.z = zz;
	orthHelperPerpendicular(temp, ax, d_d1, d_N2, z_loop, zz, &jelppi, local_sino, local_ele, apu, uu,
		d_NN,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
		d_OSEM,
#endif
		local_ind, ind, no_norm, Summ, FPbool, RHS, global_factor,
#ifdef MBSREM
		MethodListOpenCL, d_alku, &axCOSEM, d_E, d_co, d_aco, &minimi, MBSREM_prepass, idx
#else
		d_rhs_OSEM, im_dim, MethodList
#endif
#if !defined(CT)
		, d_atten
#endif
	);
#ifdef CRYSTXY
	}
for (uint uu = apu + 1; uu < maksimiXY; uu++) {
	float local_ele = computeElementOrth3DPer(s, diff, crystal_size_z, center1[uu], center2, z_center[zz]);
#ifdef VOL
	if (local_ele >= bmax)
		break;
	if (local_ele < bmin)
		local_ele = Vmax;
	else
		local_ele = V[convert_uint_rte((local_ele - bmin) * CC)];
#else
	if (local_ele <= THR)
		break;
#endif
	uint ind = uu * d_N + zz * d_Nxy;
	int4 local_ind = { 0, 0, 0, 0 };
	if (d_N == 1)
		local_ind.x = uu;
	else
		local_ind.y = uu;
	local_ind.z = zz;
	//uint local_ind = uu * d_N + zz * d_N1x;
	orthHelperPerpendicular(temp, ax, d_d1, d_N2, z_loop, zz, &jelppi, local_sino, local_ele, apu, uu,
		d_NN,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
		d_OSEM,
#endif
		local_ind, ind, no_norm, Summ, FPbool, RHS, global_factor,
#ifdef MBSREM
		MethodListOpenCL, d_alku, &axCOSEM, d_E, d_co, d_aco, &minimi, MBSREM_prepass, idx
#else
		d_rhs_OSEM, im_dim, MethodList
#endif
#if !defined(CT)
		, d_atten
#endif
	);
}
	}
#endif
	//if (z_loop == 28 && s.x > 171.651f && s.x < 171.652f) {
	//	printf("ax = %f\n", *ax);
	//}
#ifdef CRYSTZ
	for (uint zz = z_loop + 1u; zz < maksimiZ; zz++) {
#ifdef CRYSTXY
		for (int uu = (apu); uu >= minimiXY; uu--) {
#else
		uu = convert_int_sat(apu);
#endif
		float local_ele = computeElementOrth3DPer(s, diff, crystal_size_z, center1[uu], center2, z_center[zz]);
#ifdef VOL
		if (local_ele >= bmax)
			break;
		if (local_ele < bmin)
			local_ele = Vmax;
		else
			local_ele = V[convert_uint_rte((local_ele - bmin) * CC)];
#else
		if (local_ele <= THR)
			break;
#endif
		uint ind = uu * d_N + zz * d_Nxy;
		int4 local_ind = { 0, 0, 0, 0 };
		if (d_N == 1)
			local_ind.x = uu;
		else
			local_ind.y = uu;
		local_ind.z = zz;
		//uint local_ind = uu * d_N + zz * d_N1x;
		orthHelperPerpendicular(temp, ax, d_d1, d_N2, z_loop, zz, &jelppi, local_sino, local_ele, apu, uu,
			d_NN,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
			d_OSEM,
#endif
			local_ind, ind, no_norm, Summ, FPbool, RHS, global_factor,
#ifdef MBSREM
			MethodListOpenCL, d_alku, &axCOSEM, d_E, d_co, d_aco, &minimi, MBSREM_prepass, idx
#else
			d_rhs_OSEM, im_dim, MethodList
#endif
#if !defined(CT)
			, d_atten
#endif
		);
		}
	for (uint uu = apu + 1; uu < maksimiXY; uu++) {
		float local_ele = computeElementOrth3DPer(s, diff, crystal_size_z, center1[uu], center2, z_center[zz]);
#ifdef VOL
		if (local_ele >= bmax)
			break;
		if (local_ele < bmin)
			local_ele = Vmax;
		else
			local_ele = V[convert_uint_rte((local_ele - bmin) * CC)];
#else
		if (local_ele <= THR)
			break;
#endif
		uint ind = uu * d_N + zz * d_Nxy;
		int4 local_ind = { 0, 0, 0, 0 };
		if (d_N == 1)
			local_ind.x = uu;
		else
			local_ind.y = uu;
		local_ind.z = zz;
		//uint local_ind = uu * d_N + zz * d_N1x;
		orthHelperPerpendicular(temp, ax, d_d1, d_N2, z_loop, zz, &jelppi, local_sino, local_ele, apu, uu,
			d_NN,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
			d_OSEM,
#endif
			local_ind, ind, no_norm, Summ, FPbool, RHS, global_factor,
#ifdef MBSREM
			MethodListOpenCL, d_alku, &axCOSEM, d_E, d_co, d_aco, &minimi, MBSREM_prepass, idx
#else
			d_rhs_OSEM, im_dim, MethodList
#endif
#if !defined(CT)
			, d_atten
#endif
		);
	}
	}
#endif
#ifdef MBSREM
	if (!RHS) {
#ifndef CT
		* temp = 1.f / *temp;
#ifdef ATN
		* temp *= native_exp(jelppi);
#endif
#ifdef NORM
		* temp *= d_norm;
#endif
#ifdef SCATTER
		* temp *= d_scat[idx];
#endif
		* temp *= global_factor;
#endif
	}
	else {
		if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1)
			d_Amin[idx] = *minimi;
		if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u) {
#ifdef RANDOMS
			* ax += d_sc_ra[idx];
#endif
			d_ACOSEM_lhs[idx] = *ax;
		}
	}
#else
#ifndef CT
	if (FPbool) {
		*temp = 1.f / *temp;
#ifdef ATN
		* temp *= native_exp(jelppi);
#endif
#ifdef NORM
		* temp *= d_norm;
#endif
#ifdef SCATTER
		* temp *= d_scat[idx];
#endif
		* temp *= global_factor;
	}
#endif
#endif
}