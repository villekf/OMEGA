/*******************************************************************************************************************************************
* General functions for all the OpenCL kernel files. Contains functions that compute the necessary source and detector coordinates, atomics,
* forward and backward projections. Special functions are available for different cases such as TOF, listmode data, CT data, etc.
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

#ifdef ATOMIC
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
//#define TH 100000000000.f
#endif
#define THR 0.01f
#ifndef N_REKOS
#define N_REKOS 1
#endif
#ifndef NBINS
#define NBINS 1
#endif
#define NROLLS (N_REKOS * NBINS)
#ifdef PRECOMPUTE
#define TYPE 1
#else
#define TYPE 0
#endif
//#include "general_opencl_functions.h"
#ifdef VOL
#define CC 1e3f
#endif
#define TRAPZ_BINS 4.f
#ifdef PITCH
#define NA 6
#else
#define NA 2
#endif
#if defined(PTYPE41)
#define T float4
#else
#define T int4
#endif
__constant sampler_t samplerIm = CLK_NORMALIZED_COORDS_TRUE | CLK_FILTER_LINEAR | CLK_ADDRESS_CLAMP;
__constant sampler_t samplerForw = CLK_NORMALIZED_COORDS_TRUE | CLK_FILTER_LINEAR | CLK_ADDRESS_CLAMP;
//#else
__constant sampler_t samplerSiddon = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP_TO_EDGE;
//#endif

__constant sampler_t sampler_MASK = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_NONE;


// This function was taken from: https://streamhpc.com/blog/2016-02-09/atomic-operations-for-floats-in-opencl-improved/
// Computes the atomic_add for floats
#if !defined(ATOMIC) && !defined(ATOMIC32) && defined(SIDDON)
void atomicAdd_g_f(volatile __global float *addr, float val) {
	union {
		unsigned int u32;
		float        f32;
	} next, expected, current;
	current.f32 = *addr;
	do {
		expected.f32 = current.f32;
		next.f32 = expected.f32 + val;
		current.u32 = atomic_cmpxchg((volatile __global unsigned int *)addr, expected.u32, next.u32);
	} while (current.u32 != expected.u32);
}
#endif

#if defined(AF)//&& defined(SIDDON)
// Computes the forward projection
//void forwardProject(const float local_ele, float* ax, const uint kk, const uint local_ind, const __global float* d_OSEM) {
void forwardProject(const float local_ele, float* ax, const uint kk, const T local_ind, __read_only image3d_t d_OSEM) {
#ifdef PTYPE4
	ax[kk] += (local_ele * read_imagef(d_OSEM, samplerForw, local_ind).x);
#else
	ax[kk] += (local_ele * read_imagef(d_OSEM, samplerSiddon, local_ind).x);
#endif
	//ax[kk] += (local_ele);
	//ax[kk] += (local_ele * d_OSEM[local_ind]);
}

#ifdef BP
// Computes y / (x + r), where x is the forward projection and r randoms and/or scatter
void yDivFP(float* ax, const float d_Sino, const uint kk, const float local_rand) {
#if defined(RANDOMS) && defined(BP) && defined(FP)
	ax[kk] += local_rand;
#endif
#ifdef CT
	ax[kk] = native_exp(-ax[kk]) / d_Sino;
#else
	ax[kk] = d_Sino / ax[kk];
	//float joku = d_Sino;
#endif
}
#endif

// Same as above, but for special MBSREM/MRAMLA case
//#ifdef MRAMLA
//void yDivFPMBSREM(float* ax, const float d_Sino, const uint kk, const float local_rand, const float d_epsilon_mramla) {
//#if defined(RANDOMS) && defined(BP) && defined(FP)
//	ax[kk] += local_rand;
//#endif
//#ifdef CT
//	ax[kk] = native_exp(-d_epsilon_mramla) / d_Sino - (native_exp(-d_epsilon_mramla) / d_Sino) * (ax[kk] - d_epsilon_mramla);
//#else
//	ax[kk] = d_Sino / d_epsilon_mramla - (d_Sino / native_powr(d_epsilon_mramla, 2)) * (ax[kk] - d_epsilon_mramla);
//#endif
//}
//#endif
//#ifdef MBSREM
//// Struct for boolean operators indicating whether a certain method is selected (OpenCL)
//typedef struct _RecMethodsOpenCL {
//	char MLEM, OSEM, MRAMLA_, RAMLA, ROSEM, RBI, DRAMA, COSEM, ECOSEM, ACOSEM;
//	char MRP, Quad, Huber, L, FMH, WeightedMean, TV, AD, APLS, TGV, NLM;
//	char OSLMLEM, OSLOSEM, MBSREM_, BSREM, ROSEMMAP, RBIOSL, OSLCOSEM, PKMA;
//} RecMethodsOpenCL;
//#endif

#ifndef MBSREM
// Denominator (forward projection)
//void denominator(float local_ele, float* ax, uint local_ind, const uint d_N, const __global float* d_OSEM) {
void denominator(float local_ele, float* ax, T localInd, const uint d_N, __read_only image3d_t d_OSEM) {
#ifdef NREKOS1
	forwardProject(local_ele, ax, 0, localInd, d_OSEM);
#elif defined(NREKOS2)
	forwardProject(local_ele, ax, 0, localInd, d_OSEM);
	localInd.z += d_N;
	forwardProject(local_ele, ax, 1, localInd, d_OSEM);
	//forwardProject(local_ele, ax, 1, local_ind + d_N, d_OSEM);
#else
#pragma unroll N_REKOS
	for (uint kk = 0; kk < N_REKOS; kk++) {
		forwardProject(local_ele, ax, kk, localInd, d_OSEM);
		localInd.z += d_N;
		//local_ind += d_N;
	}
#endif
}

// Nominator (backprojection) in MLEM
void nominator(float* ax, const float d_Sino, const float d_epsilon_mramla, const float d_epps, 
	const float temp, const size_t idx
#if defined(BP) && defined(FP)
	, __constant uchar* MethodList
#if defined(RANDOMS)
	, const __global float* restrict d_sc_ra
#endif
#endif
) {
#if defined(BP) && defined(FP)
	float local_rand = 0.f;
#endif
#if defined(RANDOMS) && defined(BP) && defined(FP)
	local_rand = d_sc_ra[idx];
#endif
#ifdef NREKOS1
#ifndef CT
	ax[0] *= temp;
#if defined(BP) && defined(FP)
	if (ax[0] < d_epps)
		ax[0] = d_epps;
#endif
#endif
#if defined(BP) && defined(FP)
#ifdef MRAMLA
	if (MethodList[0] != 1u)
		yDivFP(ax, d_Sino, 0, local_rand);
	else if (MethodList[0] == 1u) { // MRAMLA/MBSREM
		if (ax[0] <= d_epsilon_mramla && local_rand == 0.f && d_Sino > 0.f)
			yDivFPMBSREM(ax, d_Sino, 0, local_rand, d_epsilon_mramla);
		else
			yDivFP(ax, d_Sino, 0, local_rand);
	}
#else
	yDivFP(ax, d_Sino, 0, local_rand);
#endif
#endif
#elif defined(NREKOS2)
#ifndef CT
	ax[0] *= temp;
#if defined(BP) && defined(FP)
	if (ax[0] < d_epps)
		ax[0] = d_epps;
#endif
#endif
#if defined(BP) && defined(FP)
#ifdef MRAMLA
	if (MethodList[0] != 1u)
		yDivFP(ax, d_Sino, 0, local_rand);
	else if (MethodList[0] == 1u) { // MRAMLA/MBSREM
		if (ax[0] < d_epsilon_mramla && local_rand == 0.f && d_Sino > 0.f)
			yDivFPMBSREM(ax, d_Sino, 0, local_rand, d_epsilon_mramla);
		else
			yDivFP(ax, d_Sino, 0, local_rand);
	}
#else
	yDivFP(ax, d_Sino, 0, local_rand);
#endif
#endif
#ifndef CT
	ax[1] *= temp;
#if defined(BP) && defined(FP)
	if (ax[1] < d_epps)
		ax[1] = d_epps;
#endif
#endif
#if defined(BP) && defined(FP)
#ifdef MRAMLA
	if (MethodList[1] != 1u)
		yDivFP(ax, d_Sino, 1, local_rand);
	else if (MethodList[1] == 1u) { // MRAMLA/MBSREM
		if (ax[1] < d_epsilon_mramla && local_rand == 0.f && d_Sino > 0.f)
			yDivFPMBSREM(ax, d_Sino, 1, local_rand, d_epsilon_mramla);
		else
			yDivFP(ax, d_Sino, 1, local_rand);
	}
#else
	yDivFP(ax, d_Sino, 1, local_rand);
#endif
#endif
#else
#if (defined(BP) && defined(FP)) || !defined(CT)
#pragma unroll N_REKOS
	for (uint kk = 0; kk < N_REKOS; kk++) {
#ifndef CT
		ax[kk] *= temp;
#if defined(BP) && defined(FP)
		if (ax[kk] < d_epps)
			ax[kk] = d_epps;
#endif
#endif
#if defined(BP) && defined(FP)
#ifdef MRAMLA
		if (MethodList[kk] != 1u)
			yDivFP(ax, d_Sino, kk, local_rand);
		else if (MethodList[kk] == 1u) { // MRAMLA/MBSREM
			if (ax[kk] < d_epsilon_mramla && local_rand == 0.f && d_Sino > 0.f)
				yDivFPMBSREM(ax, d_Sino, kk, local_rand, d_epsilon_mramla);
			else
				yDivFP(ax, d_Sino, kk, local_rand);
		}
#else
		yDivFP(ax, d_Sino, kk, local_rand);
#endif
#endif
	}
#endif
#endif
}
#endif

//#ifdef MBSREM
//// Nominator (backprojection), COSEM
//void nominator_cosem(float* axCOSEM, const float local_sino, const float d_epps, const float temp, const __global float* d_sc_ra,
//	const size_t idx) {
//#ifndef CT
//	* axCOSEM *= temp;
//	if (*axCOSEM < d_epps)
//		*axCOSEM = d_epps;
//#endif
//#if defined(RANDOMS) && defined(BP) && defined(FP)
//	* axCOSEM += d_sc_ra[idx];
//#endif
//#ifdef CT
//	* axCOSEM = native_exp(-*axCOSEM) / local_sino;
//#else
//	*axCOSEM = local_sino / *axCOSEM;
//#endif
//}
//#endif

#ifndef MBSREM
#ifdef BP
// Compute the backprojection
void rhs(__constant uchar* MethodList, const float local_ele, const float* ax, const uint local_ind,
	const uint d_N, __global CAST* d_rhs_OSEM) {
#ifdef NREKOS1
#ifdef ATOMIC
	atom_add(&d_rhs_OSEM[local_ind], convert_long(local_ele * ax[0] * TH));
	//atom_add(&d_rhs_OSEM[local_ind], convert_long(TH));
	//atom_add(&d_rhs_OSEM[local_ind], convert_long(local_ele * TH));
#elif defined(ATOMIC32)
	atomic_add(&d_rhs_OSEM[local_ind], convert_int(local_ele * ax[0] * TH));
#else
	atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * ax[0]));
#endif
#elif defined(NREKOS2)
#ifdef ATOMIC
	atom_add(&d_rhs_OSEM[local_ind], convert_long(local_ele * ax[0] * TH));
	atom_add(&d_rhs_OSEM[local_ind + d_N], convert_long(local_ele * ax[1] * TH));
#elif defined(ATOMIC32)
	atomic_add(&d_rhs_OSEM[local_ind], convert_int(local_ele * ax[0] * TH));
	atomic_add(&d_rhs_OSEM[local_ind + d_N], convert_int(local_ele * ax[1] * TH));
#else
	atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * ax[0]));
	atomicAdd_g_f(&d_rhs_OSEM[local_ind + d_N], (local_ele * ax[1]));
#endif
#else
	uint yy = local_ind;
#pragma unroll N_REKOS
	for (uint kk = 0; kk < N_REKOS; kk++) {
#ifdef ATOMIC
		atom_add(&d_rhs_OSEM[yy], convert_long(local_ele * ax[kk] * TH));
#elif defined(ATOMIC32)
		atomic_add(&d_rhs_OSEM[yy], convert_int(local_ele * ax[kk] * TH));
#else
		atomicAdd_g_f(&d_rhs_OSEM[yy], (local_ele * ax[kk]));
#endif
		yy += d_N;
	}
#endif
}
#endif
#endif

#else

#if !defined(PTYPE4) && !defined(PROJ5)
// Nominator (backprojection), multi-GPU version
void nominator_multi(float* axOSEM, const float d_Sino, const float d_epps, const float temp, 
	const size_t idx
#if defined(RANDOMS) && defined(BP) && defined(FP)
	, const __global float* d_sc_ra
#endif
) {
#ifndef CT
	*axOSEM *= temp;
#ifdef BP
	if (*axOSEM < d_epps)
		* axOSEM = d_epps;
#endif
#endif
#if defined(RANDOMS) && defined(BP) && defined(FP)
		* axOSEM += d_sc_ra[idx];
#endif
#ifdef BP
#ifdef CT
	* axOSEM = native_exp(-*axOSEM) / d_Sino;
#else
	*axOSEM = d_Sino / *axOSEM;
#endif
#endif
}
#endif
#endif

#if (defined(MBSREM) || !defined(AF)) && defined(SIDDON)
// Denominator (forward projection), multi-GPU version
//void denominator_multi(const float local_ele, float* axOSEM, const __global float* d_OSEM) {
void denominator_multi(const float local_ele, float* axOSEM, const float d_OSEM) {
	*axOSEM += (local_ele * d_OSEM);
	//*axOSEM += (local_ele * *d_OSEM);
}
#endif

// Detector coordinates for listmode data
#ifdef LISTMODE
void getDetectorCoordinatesListmode(const __global float* d_xyz, float3* s, float3* d, const size_t idx
#if defined(N_RAYS)
	, const int lorXY, const int lorZ, const float2 cr
#endif
) {
	const size_t i = idx * 6;
	*s = (float3)(d_xyz[i], d_xyz[i + 1], d_xyz[i + 2]);
	*d = (float3)(d_xyz[i + 3], d_xyz[i + 4], d_xyz[i + 5]);
#if defined(N_RAYS)
	if (N_RAYS3D > 1)
		multirayCoordinateShiftZ(s, d, lorZ, cr.y);
	if (N_RAYS2D > 1)
		multirayCoordinateShiftXY(s, d, lorXY, cr.x);
#endif
}
#endif

// Detector coordinates for CT or SPECT data
#if defined(CT) || (defined(SPECT) && !defined(SPECTMASK))
void getDetectorCoordinatesCT(__constant float* d_xyz, __constant float* d_uv, float3* s, float3* d, const int3 i, const uint d_size_x,
	const uint d_sizey, const float2 d_dPitch
#ifdef PROJ5
	, float3* dR, float3* dL, float3* dU, float3* dD
#endif
) {
	int id = i.z * 6;
	*s = (float3)(d_xyz[id], d_xyz[id + 1], d_xyz[id + 2]);
	*d = (float3)(d_xyz[id + 3], d_xyz[id + 4], d_xyz[id + 5]);
	//*s = vload3(i.z * 2, d_xyz);
	//*d = vload3(i.z * 2, &d_xyz[3]);
	const float2 indeksi = { convert_float(i.x) - convert_float(d_size_x) / 2.f + .5f, convert_float(i.y) - convert_float(d_sizey) / 2.f + .5f };
	id = i.z * NA;
#if defined(PITCH)
	const float3 apuX = (float3)(d_uv[id], d_uv[id + 1], d_uv[id + 2]);
	const float3 apuY = (float3)(d_uv[id + 3], d_uv[id + 4], d_uv[id + 5]);
	*d += apuX * indeksi.x + apuY * indeksi.y;
#if defined(PROJ5) && defined(FP)
	*dR = *d - apuX * 0.5f;
	*dL = *d + apuX * 0.5f;
	*dU = *d + apuY * 0.5f;
	*dD = *d - apuY * 0.5f;
	//*dR = (float3)((*d).x - apuX.x * 0.5f, (*d).y - apuX.y * 0.5f, (*d).z);
	//*dL = (float3)((*d).x + apuX.x * 0.5f, (*d).y + apuX.y * 0.5f, (*d).z);
	//*dU = (float3)((*d).x, (*d).y, (*d).z + d_dPitch.y * 0.5f);
	//*dD = (float3)((*d).x, (*d).y, (*d).z - d_dPitch.y * 0.5f);
#endif
	//*d.x += indeksi.x * d_uv[i.z * NA] + indeksi.y * d_uv[i.z * NA + 1];
	//*d.y += indeksi.x * d_uv[i.z * NA + 2] + indeksi.y * d_uv[i.z * NA + 3];
	//*d.z += indeksi.x * d_uv[i.z * NA + 4] + indeksi.y * d_uv[i.z * NA + 5];
#else
	const float apuX = d_uv[id];
	const float apuY = d_uv[id + 1];
	//if (i.x == 0 && i.y == 0 && i.z == 0) {
	//	printf("apuX = % f\n", apuX);
	//	printf("apuY = % f\n", apuY);
	//	printf("indeksi.x = %f\n", indeksi.x);
	//	printf("indeksi.y = %f\n", indeksi.y);
	//	printf("indeksi.x * apuX = %f\n", indeksi.x * apuX);
	//	printf("indeksi.x * apuY = %f\n", indeksi.x * apuY);
	//}
	(*d).x += indeksi.x * apuX;
	(*d).y += indeksi.x * apuY;
	(*d).z += indeksi.y * d_dPitch.y;
#if defined(PROJ5) && defined(FP)
	*dR = (float3)((*d).x - apuX * 0.5f, (*d).y - apuY * 0.5f, (*d).z);
	*dL = (float3)((*d).x + apuX * 0.5f, (*d).y + apuY * 0.5f, (*d).z);
	*dU = (float3)((*d).x, (*d).y, (*d).z + d_dPitch.y * 0.5f);
	*dD = (float3)((*d).x, (*d).y, (*d).z - d_dPitch.y * 0.5f);
#endif
#endif
}
#else
#if defined(RAW)
// Get the detector coordinates for the current (raw) measurement
void getDetectorCoordinatesRaw(__constant float *d_xy, __constant float *d_z, const __global ushort* d_L,
	const uint d_detPerRing, const size_t idx, float3* s, float3* d
#if defined(N_RAYS)
	, const int lorXY, const int lorZ, const float2 cr
#endif
) {
	// Get the current detector numbers
	uint2 detektorit = { convert_uint(d_L[idx * 2u]) - 1U, convert_uint(d_L[idx * 2u + 1u]) - 1U };
	// Which ring
	const uint2 loop = ((detektorit) / d_detPerRing);
	// same ring
	(*s).z = d_z[loop.x];
	if (loop.x == loop.y) {
		(*d).z = (*s).z;
	}
	else {
		(*d).z = d_z[loop.y];
	}
	detektorit -= loop * d_detPerRing;
	// Get the current x- and y-detector coordinates
	(*s).x = d_xy[detektorit.x * 2];
	(*s).y = d_xy[detektorit.x * 2 + 1];
	(*d).x = d_xy[detektorit.y * 2];
	(*d).y = d_xy[detektorit.y * 2 + 1];
#if defined(N_RAYS)
	if (N_RAYS3D > 1)
		multirayCoordinateShiftZ(s, d, lorZ, cr.y);
	if (N_RAYS2D > 1)
		multirayCoordinateShiftXY(s, d, lorXY, cr.x);
#endif
}
#endif


#if !defined(RAW) && !defined(LISTMODE) && !defined(CT) && !defined(SPECT) && !defined(PET) && defined(SIDDON)
// Get the detector coordinates for the current sinogram bin
void getDetectorCoordinates(const __global uint *d_xyindex, const __global ushort *d_zindex, const size_t idx,
	float3* s, float3* d, __constant float *d_xy, __constant float *d_z
#if defined(N_RAYS)
	, const int lorXY, const int lorZ, const float2 cr
#endif
#if defined(NLAYERS)
	, const uint d_sizey, const uint d_size_x
#endif
) {
	const uint ind = d_xyindex[idx] * 4;
#if defined(NLAYERS)
	const uint indz = d_zindex[idx] * 3;
	const int layer = d_z[indz];
	(*s).x = d_xy[ind + layer * d_size_x * d_sizey];
	(*s).y = d_xy[ind + 1 + layer * d_size_x * d_sizey];
	(*d).x = d_xy[ind + 2 + layer * d_size_x * d_sizey];
	(*d).y = d_xy[ind + 3 + layer * d_size_x * d_sizey];
	(*s).z = d_z[indz + 1];
	(*d).z = d_z[indz + 2];
#else
	const uint indz = d_zindex[idx] * 2;
	(*s).x = d_xy[ind];
	(*s).y = d_xy[ind + 1];
	(*d).x = d_xy[ind + 2];
	(*d).y = d_xy[ind + 3];
	(*s).z = d_z[indz];
	(*d).z = d_z[indz + 1];
#endif
#if defined(N_RAYS)
	if (N_RAYS3D > 1)
		multirayCoordinateShiftZ(s, d, lorZ, cr.y);
	if (N_RAYS2D > 1)
		multirayCoordinateShiftXY(s, d, lorXY, cr.x);
#endif
}
#endif
#endif

#if defined(N_RAYS) && defined(SIDDON)
void multirayCoordinateShiftXY(float3* s, float3* d, const int lor, const float cr) {
	float interval = cr / (convert_float(N_RAYS2D * 2));
	(*s).x += (interval - cr / 2.f);
	(*d).x += (interval - cr / 2.f);
	(*s).y += (interval - cr / 2.f);
	(*d).y += (interval - cr / 2.f);
	interval *= 2.f;
	(*s).x += interval * lor;
	(*d).x += interval * lor;
	(*s).y += interval * lor;
	(*d).y += interval * lor;
}

void multirayCoordinateShiftZ(float3* s, float3* d, const int lor, const float cr) {
	float interval = cr / (convert_float(N_RAYS3D * 2));
	(*s).z += (interval - cr / 2.f);
	(*d).z += (interval - cr / 2.f);
	interval *= 2.f;
	(*s).z += interval * lor;
	(*d).z += interval * lor;
}
#endif

#if (defined(FIND_LORS) || !defined(SUBSETS)) && !defined(CT)
//#ifdef FIND_LORS
// Get the detector coordinates for the current measurement (no subsets or using full sinogram subsets)
void getDetectorCoordinatesFullSinogram(const uint d_size_x, const int3 i, float3* s, float3* d, 
#ifdef CT
	__constant float* d_xy,
#else
	const __global float* restrict d_xy,
#endif
	__constant float* d_z
#if defined(N_RAYS)
	, const int lorXY, const int lorZ, const float2 cr
#endif
#if defined(NLAYERS)
	, const long NSinos, const uint d_sizey
#endif
) {
	const int id = (i.x + i.y * d_size_x) * 4;
#if defined(NLAYERS)
	const int idz = i.z * 3;
	const int layer = d_z[idz];
	*s = (float3)(d_xy[id + layer * d_size_x * d_sizey], d_xy[id + layer * d_size_x * d_sizey + 1], d_z[idz + 1]);
	*d = (float3)(d_xy[id + layer * d_size_x * d_sizey + 2], d_xy[id + layer * d_size_x * d_sizey + 3], d_z[idz + 2]);
#else
	const int idz = i.z * 2;
	*s = (float3)(d_xy[id], d_xy[id + 1], d_z[idz]);
	*d = (float3)(d_xy[id + 2], d_xy[id + 3], d_z[idz + 1]);
#endif
#if defined(N_RAYS)
	if (N_RAYS3D > 1)
		multirayCoordinateShiftZ(s, d, lorZ, cr.y);
	if (N_RAYS2D > 1)
		multirayCoordinateShiftXY(s, d, lorXY, cr.x);
#endif
}
#endif

#if defined(SIDDON)
// Compute the voxel index where the current perpendicular measurement starts
int perpendicular_start(const float d_b, const float d, const float d_d, const uint d_N) {
	int tempi = 0;
	float start = d_b - d + d_d;
	for (uint ii = 0u; ii < d_N; ii++) {
		if (start > 0.f) {
			tempi = convert_int(ii);
			break;
		}
		start += d_d;
	}
	return tempi;
}

// Compute the probability for the perpendicular elements
void perpendicular_elements(const float d_b, const float d_d1, const uint d_N1, const float d, const float d_d2, const uint d_N2, 
	float* templ_ijk, int4* tempi, int* z_loop, const uint d_N, const uint d_NN,
	const size_t idx, const float global_factor
#if defined(SCATTER)
	, const __global float* d_scat
#endif
#if !defined(CT) && defined(ATN)
	, __read_only image3d_t d_atten
#endif
#if defined(NORM)
	, const __global float* d_norm
#endif
) {
	int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	*z_loop = convert_int_sat(apu) * d_N + *z_loop * d_N1 * d_N2;
	if (d_N == 1)
		(*tempi).x = apu;
	else
		(*tempi).y = apu;
#ifdef CT
	* templ_ijk = d_d2;
#else
	float temp = d_d2 * convert_float(d_N2);
	// Probability
#if !defined(N_RAYS)
	temp = 1.f / temp;
#endif
#ifdef ATN
		float jelppi = 0.f;
		int4 atnind = *tempi;
		for (uint iii = 0u; iii < d_N2; iii++) {
			//jelppi += (*templ_ijk * (-d_atten[*tempk + iii * d_NN]));
			if (d_NN == 1)
				atnind.x = iii;
			else
				atnind.y = iii;
			jelppi += (*templ_ijk * (-read_imagef(d_atten, samplerSiddon, atnind).x));
			//jelppi += (*templ_ijk * 0.0095f);
		}
#if !defined(N_RAYS)
		temp *= native_exp(jelppi);
#endif
#endif
#if !defined(N_RAYS)
#ifdef NORM
		temp *= d_norm[idx];
#endif
#ifdef SCATTER
		temp *= d_scat[idx];
#endif
	temp *= global_factor;
	*templ_ijk = temp * d_d2;
#else
	*templ_ijk = temp;
#endif
#endif
}

#ifdef N_RAYS
// Compute the probability for the perpendicular elements
float perpendicular_elements_multiray(const float d_b, const float d_d1, const uint d_N1, const float d, const float d_d2, const uint d_N2, 
	const __global float* d_atten, uint* tempk, const uint z_loop, const uint d_N, const uint d_NN, float* jelppi) {
	int apu = perpendicular_start(d_b, d, d_d1, d_N1);

	*tempk = convert_uint_sat(apu) * d_N + z_loop * d_N1 * d_N2;

	return d_d2 * convert_float(d_N2);
}
#endif

// Compute functions (9) and (29) (detector larger than source)
void d_g_s(const float tmin, const float t_min, uint* v_min, float* t_0, int* v_u, const float diff, const float b, const float d, const float s) {

	if (tmin == t_min)
		// (11)
		*v_min = 1u;
	else {
		// (2) and (19)
		const float p_t = s + tmin * (diff);
		// (12)
		*v_min = convert_uint_sat_rtp((p_t - b) / d);
	}
	// (9)
	*t_0 += ((convert_float(*v_min) * d) / (diff));
	//  (29)
	*v_u = 1;
}

// Compute functions (9) and (29) (source larger than detector)
void s_g_d(const float tmin, const float t_min, uint* v_max, float* t_0, int* v_u, const float diff, const float b, const float d, const float s, 
	const uint N) {

	if (tmin == t_min)
		// (15)
		*v_max = N - 1u;
	else {
		// (2) and (19)
		const float p_t = s + tmin * (diff);
		// (16)
		*v_max = convert_uint_sat_rtz((p_t - b) / d);
	}
	// (9)
	*t_0 += ((convert_float(*v_max) * d) / (diff));
	// (29)
	*v_u = -1;
}

// same as above, but for precomputation phase
void d_g_s_precomp(const float tmin, const float t_min, const float tmax, const float t_max, uint* v_min, uint* v_max, float* t_0, int* v_u, 
	const float diff, const float b, const float d, const float s, const uint N) {

	if (tmin == t_min)
		// (11)
		*v_min = 1u;
	else {
		// (2) and (19)
		const float p_t = s + tmin * (diff);
		// (12)
		*v_min = convert_uint_sat_rtp((p_t - b) / d);
	}
	if (tmax == t_max)
		// (13)
		*v_max = N;
	else {
		// (2) and (19)
		const float p_t = s + tmax * (diff);
		// (14)
		*v_max = convert_uint_sat_rtz((p_t - b) / d);
	}
	// (9)
	*t_0 += ((convert_float(*v_min) * d) / (diff));
	//  (29)
	*v_u = 1;
}

// same as above, but for precomputation phase
void s_g_d_precomp(const float tmin, const float t_min, const float tmax, const float t_max, uint* v_min, uint* v_max, float* t_0, int* v_u, 
	const float diff, const float b, const float d, const float s, const uint N) {

	if (tmin == t_min)
		// (15)
		*v_max = N - 1u;
	else {
		// (2) and (19)
		const float p_t = s + tmin * (diff);
		// (16)
		*v_max = convert_uint_sat_rtz((p_t - b) / d);
	}
	if (tmax == t_max)
		// (17)
		*v_min = 0u;
	else {
		// (2) and (19)
		const float p_t = s + tmax * (diff);
		// (18)
		*v_min = convert_uint_rtp((p_t - b) / d);
	}
	// (9)
	*t_0 += ((convert_float(*v_max) * d) / (diff));
	// (29)
	*v_u = -1;
}

// Compute the index of the current voxel
uint compute_ind(const int tempj, const int tempi, const int tempk, const uint d_N1, const uint d_N2, const uint d_N, const uint d_Nx, 
	const uint d_Nyx) {
	uint local_ind = convert_uint_sat(tempj) * d_Nx + convert_uint_sat(tempi) + convert_uint_sat(tempk) * d_Nyx;
//#ifndef PRECOMPUTE
//	if (local_ind >= d_N) {
//		if (local_ind - d_N1 >= d_N)
//			local_ind -= (d_N1 * d_N2);
//		else if (local_ind - 1u >= d_N)
//			local_ind -= d_N1;
//		else
//			local_ind--;
//	}
//#endif
	return local_ind;
}

#if defined(ATN) && !defined(CT)
float compute_matrix_element(const float t0, const float tc, const float L) {
	return (t0 - tc) * L;
}

void compute_attenuation(float* tc, float* jelppi, const float LL, const float t0, const int tempi, const int tempj, const int tempk, const uint Nx, 
//void compute_attenuation(float* tc, float* jelppi, const float LL, const float t0, const int tempi, const int tempj, const int tempk, const uint Nx, 
	//const uint Nyx, const __global float* d_atten) {
	const uint Nyx, __read_only image3d_t d_atten) {
	*jelppi += (compute_matrix_element(t0, *tc, LL) * -read_imagef(d_atten, samplerSiddon, (int4)(tempi, tempj, tempk, 0)).x);
	*tc = t0;
}
#endif

#ifdef SIDDON
// compute the distance that the ray traverses in the current voxel
float compute_element(float* t0, float* tc, const float L, const float tu, const int u, int* temp_ijk, float* temp) {
	float local_ele = (*t0 - *tc) * L;
	*temp_ijk += u;
	*tc = *t0;
	*t0 += tu;
#ifndef CT
	*temp += local_ele;
#endif
	return local_ele;
}

// compute the probability of emission in the current voxel
float compute_element_2nd(float* t0, float* tc, const float L, const float tu, const int u, int* temp_ijk, const float temp) {
#ifdef CT
	float local_ele = (*t0 - *tc) * L;
#else
	float local_ele = (*t0 - *tc) * L * temp;
#endif
	*temp_ijk += u;
	*tc = *t0;
	*t0 += tu;
	return local_ele;
}
#endif

// compute the initial voxel index (beginning of the ray)
int voxel_index(const float pt, const float diff, const float d, const float apu) {
	return convert_int_rtz((pt * diff - apu) / d);
}

bool siddon_pre_loop_2D(const float b1, const float b2, const float diff1, const float diff2, const float max1, const float max2,
	const float d1, const float d2, const uint N1, const uint N2, int* temp1, int* temp2, float* t1u, float* t2u, uint* Np,
	const int TYYPPI, const float ys, const float xs, const float yd, const float xd, float* tc, int* u1, int* u2, float* t10, float* t20, bool* xy) {
	// If neither x- nor y-directions are perpendicular
// Correspond to the equations (9) and (10) from reference [2]
	//const float2 apuT = { b1 - xs, b2 - ys };
	const float apu_tx = b1 - xs;
	const float apu_ty = b2 - ys;
	*t10 = (apu_tx) / (diff1);
	*t20 = (apu_ty) / (diff2);
	const float txback = (max1 - xs) / (diff1);
	const float tyback = (max2 - ys) / (diff2);

	// Equations (5-8)
	const float txmin = fmin(*t10, txback);
	const float txmax = fmax(*t10, txback);
	const float tymin = fmin(*t20, tyback);
	const float tymax = fmax(*t20, tyback);

	// (3-4)
	*tc = fmax(txmin, tymin);
	const float tmax = fmin(txmax, tymax);
#ifdef ORTH
	if (*tc == *t10 || *tc == txback)
		*xy = true;
	else
		*xy = false;
#endif
	//int3 i = { get_global_id(0), get_global_id(1), get_global_id(2) };
	//if (i.z == 0 && i.y == 100 && i.x == 150) {
	//	printf("tx0 = %f\n", *t10);
	//	printf("ty0 = %f\n", *t20);
	//	printf("tc = %f\n", *tc);
	//	printf("tBack.x = %f\n", txback);
	//	printf("tBack.y = %f\n", tyback);
	//}

	uint imin, imax, jmin, jmax;

	if (TYYPPI == 0) {
		// If true, then the ray/LOR does not intersect the pixel space --> continue to the next LOR
		if (*tc >= tmax) {
			return true;
		}

		// (11-14)
		if (xs < xd)
			d_g_s_precomp(*tc, txmin, tmax, txmax, &imin, &imax, t10, u1, diff1, b1, d1, xs, N1);
		// (15-18)
		else
			s_g_d_precomp(*tc, txmin, tmax, txmax, &imin, &imax, t10, u1, diff1, b1, d1, xs, N1);

		//Same as above
		if (ys < yd)
			d_g_s_precomp(*tc, tymin, tmax, tymax, &jmin, &jmax, t20, u2, diff2, b2, d2, ys, N2);
		else
			s_g_d_precomp(*tc, tymin, tmax, tymax, &jmin, &jmax, t20, u2, diff2, b2, d2, ys, N2);

		*Np = imax + 1u + jmax + 1u - imin - jmin;
	}
	else {
		// (11-14)
		if (xs < xd)
			d_g_s(*tc, txmin, &imin, t10, u1, diff1, b1, d1, xs);
		// (15-18)
		else
			s_g_d(*tc, txmin, &imax, t10, u1, diff1, b1, d1, xs, N1);
		//Same as above
		if (ys < yd)
			d_g_s(*tc, tymin, &jmin, t20, u2, diff2, b2, d2, ys);
		else
			s_g_d(*tc, tymin, &jmax, t20, u2, diff2, b2, d2, ys, N2);
	}

	// (2) and (19)
	const float pt = ((fmin(*t10, *t20) + *tc) / 2.f);

	// (26)
	*temp1 = voxel_index(pt, diff1, d1, apu_tx);
	// (27)
	*temp2 = voxel_index(pt, diff2, d2, apu_ty);

	// (28)
	*t1u = d1 / fabs(diff1);
	*t2u = d2 / fabs(diff2);

	if (TYYPPI == 0) {
		if (*temp1 < 0 || *temp2 < 0 || *temp1 >= N1 || *temp2 >= N2)
			return true;
	}

	//if (*tc == *t10 || *tc == *t20)
	//	*tc -= 1e-7f;

	return false;
}

//bool siddon_pre_loop_3D(const float bx, const float by, const float bz, const float x_diff, const float y_diff, const float z_diff,
//	const float maxxx, const float maxyy, const float bzb, const float dx, const float dy, const float dz,
//	const uint Nx, const uint Ny, const uint Nz, int* tempi, int* tempj, int* tempk, float* tyu, float* txu, float* tzu,
//	uint* Np, const int TYYPPI, const float ys, const float xs, const float yd, const float xd, const float zs, const float zd, float* tc, 
//	int* iu, int* ju, int* ku, float* tx0, float* ty0, float* tz0) {
bool siddon_pre_loop_3D(const float3 b, const float3 diff, const float3 max, const float3 dd, const uint3 N, int* tempi, int* tempj, int* tempk, 
	float* txu, float* tyu, float* tzu, uint* Np, const int TYYPPI, const float3 s, const float3 d, float* tc, int* i, int* j, int* k, float* tx0, 
	float* ty0, float* tz0, bool* xy, const int3 ii) {

	const float3 apuT = b - s;
	//const float3 t0 = native_divide(apuT, diff);
	const float3 t0 = apuT / diff;
	//*tx0 = (apu_tx) / (x_diff);
	//*ty0 = (apu_ty) / (y_diff);
	//*tz0 = (apu_tz) / (z_diff);
	const float3 tBack = native_divide(max - s, diff);
	//const float txback = (maxxx - xs) / (x_diff);
	//const float tyback = (maxyy - ys) / (y_diff);
	//const float tzback = (bzb - zs) / (z_diff);

	const float3 tMin = fmin(t0, tBack);
	const float3 tMax = fmax(t0, tBack);
	//const float txmin = fmin(*tx0, txback);
	//const float txmax = fmax(*tx0, txback);
	//const float tymin = fmin(*ty0, tyback);
	//const float tymax = fmax(*ty0, tyback);
	//const float tzmin = fmin(*tz0, tzback);
	//const float tzmax = fmax(*tz0, tzback);



	*tc = fmax(fmax(tMin.x, tMin.z), tMin.y);
	const float tmax = fmin(fmin(tMax.x, tMax.z), tMax.y);
	*tx0 = t0.x;
	*ty0 = t0.y;
	*tz0 = t0.z;
#ifdef ORTH
	if (*tc == *tx0 || *tc == tBack.x)
		*xy = true;
	else
		*xy = false;
	//int3 ii = { get_global_id(0), get_global_id(1), get_global_id(2) };
	//if (ii.z == 543 && ii.y == 4 && ii.x == 115) {
	//	printf("tx0 = %f\n", *tx0);
	//	printf("ty0 = %f\n", *ty0);
	//	printf("tc = %f\n", *tc);
	//	printf("tBack.x = %f\n", tBack.x);
	//	printf("tBack.y = %f\n", tBack.y);
	//}
#endif

	uint imin, imax, jmin, jmax, kmin, kmax;

	if (TYYPPI == 0) {
		if (*tc >= tmax) {
			return true;
		}
		//float ax = 
		if (s.x < d.x)
			d_g_s_precomp(*tc, tMin.x, tmax, tMax.x, &imin, &imax, tx0, i, diff.x, b.x, dd.x, s.x, N.x);
		else
			s_g_d_precomp(*tc, tMin.x, tmax, tMax.x, &imin, &imax, tx0, i, diff.x, b.x, dd.x, s.x, N.x);

		if (s.y < d.y)
			d_g_s_precomp(*tc, tMin.y, tmax, tMax.y, &jmin, &jmax, ty0, j, diff.y, b.y, dd.y, s.y, N.y);
		else
			s_g_d_precomp(*tc, tMin.y, tmax, tMax.y, &jmin, &jmax, ty0, j, diff.y, b.y, dd.y, s.y, N.y);

		if (s.z < d.z)
			d_g_s_precomp(*tc, tMin.z, tmax, tMax.z, &kmin, &kmax, tz0, k, diff.z, b.z, dd.z, s.z, N.z);
		else
			s_g_d_precomp(*tc, tMin.z, tmax, tMax.z, &kmin, &kmax, tz0, k, diff.z, b.z, dd.z, s.z, N.z);

		*Np = (kmax - kmin + 1) + (jmax - jmin + 1) + (imax - imin + 1);
	}
	else {
		if (s.x < d.x)
			d_g_s(*tc, tMin.x, &imin, tx0, i, diff.x, b.x, dd.x, s.x);
		else
			s_g_d(*tc, tMin.x, &imax, tx0, i, diff.x, b.x, dd.x, s.x, N.x);

		if (s.y < d.y)
			d_g_s(*tc, tMin.y, &jmin, ty0, j, diff.y, b.y, dd.y, s.y);
		else
			s_g_d(*tc, tMin.y, &jmax, ty0, j, diff.y, b.y, dd.y, s.y, N.y);

		if (s.z < d.z)
			d_g_s(*tc, tMin.z, &kmin, tz0, k, diff.z, b.z, dd.z, s.z);
		else
			s_g_d(*tc, tMin.z, &kmax, tz0, k, diff.z, b.z, dd.z, s.z, N.z);
	}

	const float pt = ((fmin(fmin(*tz0, *ty0), *tx0) + *tc) / 2.f);

	const float3 tempijkF = clamp(mad(pt, diff, -apuT) / dd, 0.f, convert_float3(N - 1));
	const int3 tempijk = convert_int3_rtz(tempijkF);
	//const int3 tempijk = convert_int3_rtz(mad(pt, diff, -apuT) / dd);
	*tempi = tempijk.x;
	*tempj = tempijk.y;
	*tempk = tempijk.z;
	//*tempijk = convert_int3_rtz((pt * diff - apuT) / dd);

	//*tempi = voxel_index(pt, diff.x, dd.x, apuT.x);
	//*tempj = voxel_index(pt, diff.y, dd.y, apuT.y);
	//*tempk = voxel_index(pt, diff.z, dd.z, apuT.z);

	//if (ii.x == 526 && ii.y == 345 && ii.z == 0) {
	//	printf("*tc = %f\n", *tc);
	//	printf("t0.x = %f\n", t0.x);
	//	printf("t0.y = %f\n", t0.y);
	//	printf("t0.z = %f\n", t0.z);
	//	printf("tmax = %f\n", tmax);
	//	//printf("((pt * diff - apuT) / dd).y = %f\n", ((pt * diff - apuT) / dd).y);
	//}

	//if (TYYPPI == 0) {
	//	//if (*tempi < 0 || *tempj < 0 || *tempk < 0 || *tempi >= N.x || *tempj >= N.y || *tempk >= N.z) {
	//	//if (*tempijk.x < 0 || *tempijk.y < 0 || *tempijk.z < 0 || *tempijk.x >= N.x || *tempijk.y >= N.y || *tempijk.z >= N.z) {
	//	if (any(tempijk < 0) || any(tempijk >= convert_int3(N))) {
	//		//*tc = (pt * y_diff - apu_ty) / dy;
	//		return true;
	//	}
	//}

	//*tu = dd / fabs(diff);
	*txu = dd.x / fabs(diff.x);
	*tyu = dd.y / fabs(diff.y);
	*tzu = dd.z / fabs(diff.z);

	return false;
}

#ifdef TOF
#define _2PI 0.3989423f

float normPDF(const float x, const float mu, const float sigma) {

	const float a = (x - mu) / sigma;

	return _2PI / sigma * native_exp(-0.5f * a * a);
}

void TOFDis(const float3 diff, const float tc, const float LL, float* D, float* DD) {
	//const float xI = x_diff * tc;
	//const float yI = y_diff * tc;
	//const float zI = z_diff * tc;
	*D = length(diff * tc) - LL / 2.f;
	*DD = *D;
}

float TOFWeight(const float element, const float sigma_x, const float D, const float DD, const float TOFCenter, float dX) {
	float output = normPDF(D, TOFCenter, sigma_x);
	dX *= sign(DD);
	for (long tr = 1L; tr < convert_long(TRAPZ_BINS) - 1; tr++)
		output += (normPDF(D - dX * convert_float(tr), TOFCenter, sigma_x) * 2.f);
	output += normPDF(D - element * sign(DD), TOFCenter, sigma_x);
	//return (element * (normPDF(D, TOFCenter, sigma_x) + normPDF(D - element * sign(DD), TOFCenter, sigma_x)) / 2.f);
	return output;
}


float TOFLoop(const float DD, const float element, __private float* TOFVal, __constant float* TOFCenter,
	const float sigma_x, float* D, const uint tid, const float epps) {
	float TOFSum = 0.f;
#ifdef DEC
	__private float apu[NBINS];
#endif
	const float dX = element / (TRAPZ_BINS - 1.f);
#pragma unroll NBINS
	for (long to = 0L; to < NBINS; to++) {
#ifdef DEC
		apu[to] = TOFWeight(element, sigma_x, *D, DD, TOFCenter[to], dX) * dX;
		TOFSum += apu[to];
#else
		const float apu = TOFWeight(element, sigma_x, *D, DD, TOFCenter[to], dX) * dX;
		TOFSum += apu;
#endif
	}
#ifdef DEC
	*D -= (element * sign(DD));
#endif
	if (TOFSum < epps)
		TOFSum = epps;
#ifdef DEC
#pragma unroll NBINS
	for (long to = 0L; to < NBINS; to++)
		TOFVal[to + tid] = apu[to] / TOFSum;
#endif
	return TOFSum;
}


//void denominatorTOF(float* ax, const float element, const __global float* d_OSEM, uint local_ind, const float TOFSum, __private float* TOFVal,
//	const float DD, __constant float* TOFCenter, const float sigma_x, float* D, const uint tid, const float epps, const uint d_N) {
void denominatorTOF(float* ax, const float element, __read_only image3d_t d_OSEM, int4 ind, const float TOFSum, __private float* TOFVal,
	const float DD, __constant float* TOFCenter, const float sigma_x, float* D, const uint tid, const float epps, const uint d_N) {
#if defined(AF) && !defined(MBSREM)
	uint ll = NBINS;
#pragma unroll N_REKOS
	for (uint kk = 0U; kk < N_REKOS; kk++) {
		uint ii = ll * kk;
#else
	const uint ii = 0U;
#endif
	//float apu = element * d_OSEM[local_ind];
	float apu = element * read_imagef(d_OSEM, samplerSiddon, ind).x;
#ifndef DEC
	const float dX = element / (TRAPZ_BINS - 1.f);
#endif
#pragma unroll NBINS
	for (long to = 0L; to < NBINS; to++) {
#ifdef DEC
		ax[to + ii] += (apu * TOFVal[to + tid]);
#else
		const float jelppi = (TOFWeight(element, sigma_x, *D, DD, TOFCenter[to], dX) * dX) / TOFSum;
		ax[to + ii] += (apu * jelppi);
#endif
	}
#if defined(AF) && !defined(MBSREM)
	ind.z += d_N;
		//local_ind += d_N;
	}
#endif
#ifndef DEC
*D -= (element * sign(DD));
#endif
}



// Nominator (y for backprojection)
void nominatorTOF(__constant uchar* MethodList, float* ax, const __global float* d_Sino, const float d_epsilon_mramla, const float d_epps,
	const float temp, const size_t idx, const long TOFSize, const float local_sino
#if defined(RANDOMS) && defined(BP) && defined(FP)
	, const __global float* restrict d_sc_ra
#endif
	{
	float local_rand = 0.f;
#if defined(RANDOMS) && defined(BP) && defined(FP)
	local_rand = d_sc_ra[idx];
#endif
#if defined(AF) && !defined(MBSREM)
	uint ll = NBINS;
#pragma unroll N_REKOS
	for (uint kk = 0; kk < N_REKOS; kk++) {
		uint ii = ll * kk;
#else
	const uint ii = 0U;
#endif
#pragma unroll NBINS
	for (long to = 0L; to < NBINS; to++) {
		ax[to + ii] *= temp;
		if (ax[to + ii] < d_epps)
			ax[to + ii] = d_epps;
#if defined(RANDOMS) && defined(BP) && defined(FP)
		ax[to + ii] += local_rand;
#endif
#if defined(AF) && defined(MRAMLA) && !defined(MBSREM)
		if (MethodList[kk] != 1u)
			ax[to + ii] = d_Sino[idx + to * TOFSize] / ax[to + ii];
		else if (MethodList[kk] == 1u) { // MRAMLA/MBSREM
			if (ax[to + ii] <= d_epsilon_mramla && local_rand == 0.f && local_sino > 0.f)
				ax[to + ii] = d_Sino[idx + to * TOFSize] / d_epsilon_mramla - (d_Sino[idx + to * TOFSize] / native_powr(d_epsilon_mramla, 2))
				* (ax[to + ii] - d_epsilon_mramla);
			else
				ax[to + ii] = d_Sino[idx + to * TOFSize] / ax[to + ii];
		}
#else
#if defined(AF) || defined(BP)
		ax[to + ii] = d_Sino[idx + to * TOFSize] / ax[to + ii];
#endif
#endif
	}
#if defined(AF) && !defined(MBSREM)
	}
#endif
}


void backprojectTOF(const uint local_ind, const int4 localInd, const float local_ele, const uint tid, const __private float* TOFVal, 
	const float* yax, __global CAST* d_Summ, const float local_sino
#ifndef DEC
	, const float temp, const float sigma_x, float* D, const float DD, __constant float* TOFCenter, const float epps, const float TOFSum
#endif
#ifdef MBSREM
	, const RecMethodsOpenCL MethodListOpenCL, const uint d_alku, const uchar MBSREM_prepass, float* minimi, float* axACOSEM, 
	const __global float* d_OSEM, __global float* d_E, __global CAST* d_co, __global CAST* d_aco, const size_t idx, const long TOFSize
#else
	, __global CAST* d_rhs, const uchar no_norm, const uint d_N
#endif
	) {
	uint yy = local_ind;
	float val = 0.f;
#if defined(AF) && !defined(MBSREM)
	uint ll = NBINS;
#pragma unroll N_REKOS
	for (uint kk = 0; kk < N_REKOS; kk++) {
		uint ii = ll * kk;
#else
	const uint ii = 0U;
#endif
#ifndef DEC
	const float dX = element / (TRAPZ_BINS - 1.f);
#endif

	float yaxTOF = 0.f;
#pragma unroll NBINS
	for (long to = 0L; to < NBINS; to++) {

#ifdef DEC
		const float apu = local_ele * TOFVal[to + tid];
#else
		const float apu = local_ele * ((TOFWeight(local_ele / temp, sigma_x, *D, DD, TOFCenter[to], dX) * dX) / TOFSum);
#endif

#ifdef MBSREM
		if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1 && d_alku == 0u) {
			if (apu < minimi[to] && apu > 0.f)
				minimi[to] = apu;
			d_E[idx + to * TOFSize] += apu;
		}
		if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u)
			axACOSEM[to] += (apu * read_imagef(d_OSEM, samplerSiddon, localInd).x);
			//axACOSEM[to] += (apu * d_OSEM[local_ind]);
#endif

		val += apu;
		yaxTOF += apu * yax[to + ii];
	}

#ifdef MBSREM
	if (d_alku == 0u) {
		if (MBSREM_prepass == 1)
#ifdef ATOMIC
			atom_add(&d_Summ[local_ind], convert_long(val * TH));
#elif defined(ATOMIC32)
			atomic_add(&d_Summ[local_ind], convert_int(val * TH));
#else
			atomicAdd_g_f(&d_Summ[local_ind], val);
#endif
		if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino != 0.f)
#ifdef ATOMIC
			atom_add(&d_co[local_ind], convert_long(yaxTOF * TH));
#elif defined(ATOMIC32)
			atomic_add(&d_co[local_ind], convert_int(yaxTOF * TH));
#else
			atomicAdd_g_f(&d_co[local_ind], yaxTOF);
#endif
		if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino != 0.f)
#ifdef ATOMIC
			atom_add(&d_aco[local_ind], convert_long(yaxTOF * TH));
#elif defined(ATOMIC32)
			atomic_add(&d_aco[local_ind], convert_int(yaxTOF * TH));
#else
			atomicAdd_g_f(&d_aco[local_ind], yaxTOF);
#endif
	}
#else
	if (no_norm == 0u && ii == 0)
#ifdef ATOMIC
		atom_add(&d_Summ[local_ind], convert_long(val * TH));
#elif defined(ATOMIC32)
		atomic_add(&d_Summ[local_ind], convert_int(val * TH));
#else
		atomicAdd_g_f(&d_Summ[local_ind], val);
#endif
#if defined(FP) && defined(BP)
	if (local_sino != 0.f) {
#endif
#ifdef ATOMIC
		atom_add(&d_rhs[yy], convert_long(yaxTOF * TH));
#elif defined(ATOMIC32)
		atomic_add(&d_rhs[yy], convert_int(yaxTOF * TH));
#else
		atomicAdd_g_f(&d_rhs[yy], (yaxTOF));
#endif
#if defined(FP) && defined(BP)
	}
#endif
#endif

#if defined(AF) && !defined(MBSREM)
		yy += d_N;
	}
#endif

#ifndef DEC
*D -= (local_ele / temp * sign(DD));
#endif
}


void sensTOF(const uint local_ind, const int4 localInd, const float local_ele, const uint tid, const __private float* TOFVal, 
	__global CAST * d_Summ, 
#ifndef DEC
	const float temp, const float sigma_x, float* D, const float DD, __constant float* TOFCenter, const float epps, const float TOFSum, 
#endif
#ifdef MBSREM
	const RecMethodsOpenCL MethodListOpenCL, const uint d_alku, const uchar MBSREM_prepass, float* minimi, float* axACOSEM, 
	const __global float* d_OSEM, __global float* d_E, const size_t idx, const long TOFSize,
#endif
	const uchar no_norm) {
	float val = 0.f;
#ifndef DEC
	const float dX = element / (TRAPZ_BINS - 1.f);
#endif
#pragma unroll NBINS
	for (long to = 0L; to < NBINS; to++) {

#ifdef DEC
		const float apu = local_ele * TOFVal[to + tid];
#else
		const float apu = local_ele * ((TOFWeight(local_ele / temp, sigma_x, *D, DD, TOFCenter[to], dX) * dX) / TOFSum);
#endif
#ifdef MBSREM
		if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1 && d_alku == 0u) {
			if (apu < minimi[to] && apu > 0.f)
				minimi[to] = apu;
			d_E[idx + to * TOFSize] += apu;
		}
		if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u)
			axACOSEM[to] += (apu * read_imagef(d_OSEM, samplerSiddon, localInd).x);
			//axACOSEM[to] += (apu * d_OSEM[local_ind]);
#endif
		val += apu;
	}
#ifdef MBSREM
	if (no_norm == 0u && d_alku == 0u)
#ifdef ATOMIC
		atom_add(&d_Summ[local_ind], convert_long(val * TH));
#elif defined(ATOMIC32)
		atomic_add(&d_Summ[local_ind], convert_int(val * TH));
#else
		atomicAdd_g_f(&d_Summ[local_ind], val);
#endif
#else
	if (no_norm == 0u)
#ifdef ATOMIC
		atom_add(&d_Summ[local_ind], convert_long(val * TH));
#elif defined(ATOMIC32)
		atomic_add(&d_Summ[local_ind], convert_int(val * TH));
#else
		atomicAdd_g_f(&d_Summ[local_ind], val);
#endif
#endif

#ifndef DEC
*D -= (local_ele / temp * sign(DD));
#endif
}
#endif

#endif

#if defined(AF) && defined(FP) && !defined(BP)
void forwardProjectAF(__global float* output, float* ax, size_t idx, const uint N) {

#ifdef NREKOS1
	output[idx] = ax[0];
#elif defined(NREKOS2)
	output[idx] = ax[0];
	output[idx + N] = ax[1];
#else
#pragma unroll N_REKOS
	for (uint kk = 0; kk < N_REKOS; kk++) {
		output[idx + kk * N] = ax[kk];
	}
#endif
}
#endif
