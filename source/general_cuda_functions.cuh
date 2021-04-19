/**************************************************************************
* General functions for all the CUDA kernel files.
*
* Copyright (C) 2020 Ville-Veikko Wettenhovi
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
#define TRAPZ_BINS 5.f

#ifdef MBSREM
// Struct for boolean operators indicating whether a certain method is selected (OpenCL)
typedef struct _RecMethodsOpenCL {
	char MLEM, OSEM, MRAMLA_, RAMLA, ROSEM, RBI, DRAMA, COSEM, ECOSEM, ACOSEM;
	char MRP, Quad, Huber, L, FMH, WeightedMean, TV, AD, APLS, TGV, NLM;
	char OSLMLEM, OSLOSEM, MBSREM_, BSREM, ROSEMMAP, RBIOSL, OSLCOSEM, PKMA;
} RecMethodsOpenCL;
#endif

//#ifdef AF
#ifndef MBSREM
// Denominator (forward projection)
__device__ void denominator(float local_ele, float* ax, unsigned int local_ind, const unsigned int d_N, const float* d_OSEM) {
#ifdef NREKOS1
	ax[0] += (local_ele * d_OSEM[local_ind]);
#elif defined(NREKOS2)
	ax[0] += (local_ele * d_OSEM[local_ind]);
	ax[1] += (local_ele * d_OSEM[local_ind + d_N]);
#else
#pragma unroll N_REKOS
	for (unsigned int kk = 0; kk < N_REKOS; kk++) {
		ax[kk] += (local_ele * d_OSEM[local_ind]);
		local_ind += d_N;
	}
#endif
}

// Nominator (backprojection) in MLEM
__device__ void nominator(const unsigned char* MethodList, float* ax, const float d_Sino, const float d_epsilon_mramla, const float d_epps,
	const float temp, const float* d_sc_ra, const unsigned int idx) {
	float local_rand = 0.f;
#ifdef RANDOMS
		local_rand = d_sc_ra[idx];
#endif
#ifdef NREKOS1
#ifndef CT
		ax[0] *= temp;
	if (ax[0] < d_epps)
		ax[0] = d_epps;
#endif
#ifdef RANDOMS
		ax[0] += local_rand;
#endif
#ifdef MRAMLA
#ifdef CT
	if (MethodList[0] != 1u)
		ax[0] = expf(-ax[0]) / d_Sino;
	else if (MethodList[0] == 1u) { // MRAMLA/MBSREM
		if (ax[0] <= d_epsilon_mramla && local_rand == 0.f && d_Sino > 0.f)
			ax[0] = expf(-d_epsilon_mramla) / d_Sino - (expf(-d_epsilon_mramla) / d_Sino) * (ax[0] - d_epsilon_mramla);
		else
			ax[0] = expf(-ax[0]) / d_Sino;
	}
#else
	if (MethodList[0] != 1u)
		ax[0] = d_Sino / ax[0];
	else if (MethodList[0] == 1u) { // MRAMLA/MBSREM
		if (ax[0] < d_epsilon_mramla && local_rand == 0.f && d_Sino > 0.f)
			ax[0] = d_Sino / d_epsilon_mramla - (d_Sino / powf(d_epsilon_mramla, 2)) * (ax[0] - d_epsilon_mramla);
		else
			ax[0] = d_Sino / ax[0];
		ax[0];
#endif
#else
#ifdef CT
	ax[0] = expf(-ax[0]) / d_Sino;
#else
	ax[0] = d_Sino / ax[0];
#endif
#endif
#elif defined(NREKOS2)
#ifndef CT
		ax[0] *= temp;
	if (ax[0] < d_epps)
		ax[0] = d_epps;
#endif
#ifdef RANDOMS
		ax[0] += local_rand;
#endif
#ifdef MRAMLA
#ifdef CT
		if (MethodList[0] != 1u)
			ax[0] = expf(-ax[0]) / d_Sino;
		else if (MethodList[0] == 1u) { // MRAMLA/MBSREM
			if (ax[0] <= d_epsilon_mramla && local_rand == 0.f && d_Sino > 0.f)
				ax[0] = expf(-d_epsilon_mramla) / d_Sino - (expf(-d_epsilon_mramla) / d_Sino) * (ax[0] - d_epsilon_mramla);
			else
				ax[0] = expf(-ax[0]) / d_Sino;
		}
#else
	if (MethodList[0] != 1u)
		ax[0] = d_Sino / ax[0];
	else if (MethodList[0] == 1u) { // MRAMLA/MBSREM
		if (ax[0] < d_epsilon_mramla && local_rand == 0.f && d_Sino > 0.f)
			ax[0] = d_Sino / d_epsilon_mramla - (d_Sino / powf(d_epsilon_mramla, 2)) * (ax[0] - d_epsilon_mramla);
		else
			ax[0] = d_Sino / ax[0];
	}
#endif
#else
#ifdef CT
		ax[0] = expf(-ax[0]) / d_Sino;
#else
		ax[0] = d_Sino / ax[0];
#endif
#endif
#ifndef CT
	ax[1] *= temp;
	if (ax[1] < d_epps)
		ax[1] = d_epps;
#endif
#ifdef RANDOMS
		ax[1] += local_rand;
#endif
#ifdef MRAMLA
#ifdef CT
		if (MethodList[1] != 1u)
			ax[1] = expf(-ax[1]) / d_Sino;
		else if (MethodList[1] == 1u) { // MRAMLA/MBSREM
			if (ax[1] <= d_epsilon_mramla && local_rand == 0.f && d_Sino > 0.f)
				ax[1] = expf(-d_epsilon_mramla) / d_Sino - (expf(-d_epsilon_mramla) / d_Sino) * (ax[1] - d_epsilon_mramla);
			else
				ax[1] = expf(-ax[1]) / d_Sino;
		}
#else
	if (MethodList[1] != 1u)
		ax[1] = d_Sino / ax[1];
	else if (MethodList[1] == 1u) { // MRAMLA/MBSREM
		if (ax[1] < d_epsilon_mramla && local_rand == 0.f && d_Sino > 0.f)
			ax[1] = d_Sino / d_epsilon_mramla - (d_Sino / powf(d_epsilon_mramla, 2)) * (ax[1] - d_epsilon_mramla);
		else
			ax[1] = d_Sino / ax[1];
	}
#endif
#else
#ifdef CT
		ax[1] = expf(-ax[1]) / d_Sino;
#else
		ax[1] = d_Sino / ax[1];
#endif
#endif
#else
#pragma unroll N_REKOS
	for (unsigned int kk = 0; kk < N_REKOS; kk++) {
#ifndef CT
		ax[kk] *= temp;
		if (ax[kk] < d_epps)
			ax[kk] = d_epps;
#endif
#ifdef RANDOMS
			ax[kk] += local_rand;
#endif
#ifdef MRAMLA
#ifdef CT
			if (MethodList[kk] != 1u)
				ax[kk] = expf(-ax[kk]) / d_Sino;
			else if (MethodList[kk] == 1u) { // MRAMLA/MBSREM
				if (ax[kk] <= d_epsilon_mramla && local_rand == 0.f && d_Sino > 0.f)
					ax[kk] = expf(-d_epsilon_mramla) / d_Sino - (expf(-d_epsilon_mramla) / d_Sino) * (ax[kk] - d_epsilon_mramla);
				else
					ax[kk] = expf(-ax[kk]) / d_Sino;
	}
#else
		if (MethodList[kk] != 1u)
			ax[kk] = d_Sino / ax[kk];
		else if (MethodList[kk] == 1u) { // MRAMLA/MBSREM
			if (ax[kk] < d_epsilon_mramla && local_rand == 0.f && d_Sino > 0.f)
				ax[kk] = d_Sino / d_epsilon_mramla - (d_Sino / powf(d_epsilon_mramla, 2)) * (ax[kk] - d_epsilon_mramla);
			else
				ax[kk] = d_Sino / ax[kk];
		}
#endif
#else
#ifdef CT
		ax[kk] = expf(-ax[kk]) / d_Sino;
#else
		ax[kk] = d_Sino / ax[kk];
#endif
#endif
	}
#endif
}
#endif
#ifdef MBSREM
// Nominator (backprojection), COSEM
__device__ void nominator_cosem(float* axCOSEM, const float local_sino, const float d_epps, const float temp, const float* d_sc_ra,
	const unsigned int idx) {
#ifndef CT
	* axCOSEM *= temp;
	if (*axCOSEM < d_epps)
		*axCOSEM = d_epps;
#endif
#ifdef RANDOMS
	* axCOSEM += d_sc_ra[idx];
#endif
#ifdef CT
	* axCOSEM = expf(-*axCOSEM) / local_sino;
#else
	* axCOSEM = local_sino / *axCOSEM;
#endif
}
#endif
#ifndef MBSREM
__device__ void rhs(const unsigned char* MethodList, const float local_ele, const float* ax, const unsigned int local_ind,
	const unsigned int d_N, CAST* d_rhs_OSEM) {
#ifdef NREKOS1
#ifdef ATOMIC
	atomicAdd(&d_rhs_OSEM[local_ind], __float2ull_rn(local_ele * ax[0] * TH));
#else
	atomicAdd(&d_rhs_OSEM[local_ind], (local_ele * ax[0]));
#endif
#elif defined(NREKOS2)
#ifdef ATOMIC
	atomicAdd(&d_rhs_OSEM[local_ind], __float2ull_rn(local_ele * ax[0] * TH));
	atomicAdd(&d_rhs_OSEM[local_ind + d_N], __float2ull_rn(local_ele * ax[1] * TH));
#else
		atomicAdd(&d_rhs_OSEM[local_ind], (local_ele * ax[0]));
		atomicAdd(&d_rhs_OSEM[local_ind + d_N], (local_ele * ax[1]));
#endif
#else
	unsigned int yy = local_ind;
#pragma unroll N_REKOS
	for (unsigned int kk = 0; kk < N_REKOS; kk++) {
#ifdef ATOMIC
		atomicAdd(&d_rhs_OSEM[yy], __float2ull_rn(local_ele * ax[kk] * TH));
#else
		atomicAdd(&d_rhs_OSEM[yy], (local_ele * ax[kk]));
#endif
		yy += d_N;
	}
#endif
}
#endif
//#else
//
//// Nominator (backprojection), multi-GPU version
//__device__ void nominator_multi(float* axOSEM, const float d_Sino, const float d_epps, const float temp, const float* d_sc_ra,
//	const size_t idx) {
//#ifndef CT
//	* axOSEM *= temp;
//#ifdef BP
//	if (*axOSEM < d_epps)
//		*axOSEM = d_epps;
//#endif
//#endif
//#ifdef RANDOMS
//	* axOSEM += d_sc_ra[idx];
//#endif
//#ifdef BP
//#ifdef CT
//	* axOSEM = expf(-*axOSEM) / d_Sino;
//#else
//	* axOSEM = d_Sino / *axOSEM;
//#endif
//#endif
//}
//#endif
//
//#if defined(MBSREM) || !defined(AF)
//// Denominator (forward projection), multi-GPU version
//__device__ void denominator_multi(const float local_ele, float* axOSEM, const float* d_OSEM) {
//	*axOSEM += (local_ele * *d_OSEM);
//}
//#endif


#ifdef CT
__device__ void get_detector_coordinates_CT(const float* x, const float* y, const float* z, const unsigned int size_x, const size_t idx, const unsigned int subsets,
	const float* angles, const unsigned int* d_xyindex, const unsigned short* d_zindex, const unsigned int size_z, const float dPitch, const long long int nProjections,
	float* xs, float* xd, float* ys, float* yd, float* zs, float* zd) {
#ifdef LISTMODE
	*xs = x[idx];
	*xd = x[idx + size_x];
	*ys = y[idx];
	*yd = y[idx + size_x];
	*zs = z[idx];
	*zd = z[idx + size_x];
#else
	if (subsets > 1U) {
		const unsigned int lu = d_xyindex[idx];
		*xs = x[d_zindex[idx] + nProjections];
		*xd = x[d_zindex[idx]] - dPitch * float(lu % size_x) * cosf(angles[d_zindex[idx]]);
		*ys = y[d_zindex[idx] + nProjections];
		*yd = y[d_zindex[idx]] - dPitch * float(lu % size_x) * sinf(angles[d_zindex[idx]]);
		*zs = z[d_zindex[idx] + nProjections];
		*zd = z[d_zindex[idx]] + dPitch * float(lu / size_x);
	}
	else {
		const long long int ll = idx / ((size_x) * (size_z));
		const long long int lu = idx % ((size_x) * (size_z));
		*xs = x[ll + nProjections];
		*xd = x[ll] - dPitch * float(lu % size_x) * cosf(angles[ll]);
		*ys = y[ll + nProjections];
		*yd = y[ll] - dPitch * float(lu % size_x) * sinf(angles[ll]);
		*zs = z[ll + nProjections];
		*zd = z[ll] + dPitch * float(lu / size_x);
		}
#endif
}
#else
#if defined(RAW) && !defined(N_RAYS)
// Get the detector coordinates for the current (raw) measurement
__device__ void get_detector_coordinates_raw(const float* d_x, const float* d_y, const float* d_zdet, const unsigned short* d_L,
	const unsigned int d_det_per_ring, const unsigned int idx, const unsigned int* d_pseudos, const unsigned int d_pRows, float* xs, float* xd, float* ys, float* yd, float* zs,
	float* zd) {
	unsigned int ps = 0;
	// Get the current detector numbers
	const unsigned int detektorit1 = (unsigned int)(d_L[idx * 2u]) - 1u;
	const unsigned int detektorit2 = (unsigned int)(d_L[idx * 2u + 1u]) - 1u;
	// Which ring
	const unsigned int loop1 = ((detektorit1) / d_det_per_ring);
	const unsigned int loop2 = ((detektorit2) / d_det_per_ring);
	// same ring
	if (loop1 == loop2) {
		*zs = d_zdet[loop1];
		*zd = *zs;
	}
	else {
		*zs = d_zdet[loop1];
		*zd = d_zdet[loop2];
	}
	// Get the current x- and y-detector coordinates
	*xs = d_x[detektorit1 - d_det_per_ring * (loop1)];
	*xd = d_x[detektorit2 - d_det_per_ring * (loop2)];
	*ys = d_y[detektorit1 - d_det_per_ring * (loop1)];
	*yd = d_y[detektorit2 - d_det_per_ring * (loop2)];
}
#endif



#if defined(N_RAYS) && (defined(RAW) && !defined(LISTMODE))
// Get the detector coordinates for the current (raw) measurement
__device__ void get_detector_coordinates_raw_multiray(const float* d_x, const float* d_y, const float* d_zdet, const unsigned short* d_L,
	const unsigned int d_det_per_ring, const unsigned int idx, const unsigned int* d_pseudos, const unsigned int d_pRows, float* xs, float* xd, float* ys, float* yd, float* zs, float* zd,
	const unsigned short lor, const float cr_pz) {
	unsigned int ps = 0;
	// Get the current detector numbers
	const unsigned int detektorit1 = (unsigned int)(d_L[idx * 2u]) - 1u;
	const unsigned int detektorit2 = (unsigned int)(d_L[idx * 2u + 1u]) - 1u;
	// Which ring
	const unsigned int loop1 = ((detektorit1) / d_det_per_ring);
	const unsigned int loop2 = ((detektorit2) / d_det_per_ring);
	// same ring
	if (loop1 == loop2) {
		*zs = d_zdet[loop1];
		*zd = *zs;
	}
	else {
		*zs = d_zdet[loop1];
		*zd = d_zdet[loop2];
	}
	// Get the current x- and y-detector coordinates
	if (N_RAYS3D > 1) {
		unsigned short rays3D = N_RAYS3D;
		if (N_RAYS3D % 2 == 0)
			rays3D++;
		int r[N_RAYS3D + 1];
		int hh = 0;
		for (int kk = -N_RAYS3D / 2; kk <= N_RAYS3D / 2; kk++) {
			if (kk == 0 && (N_RAYS3D % 2 == 0)) {}
			else {
				r[hh] = kk;
				hh++;
			}
		}
		const float rr = (float)(r[(lor - 1) % N_RAYS3D]);
		*zs += (cr_pz * rr);
		*zd += (cr_pz * rr);
	}
#ifdef N_RAYS2D
		const unsigned short ll = (lor - 1) / N_RAYS3D;
		*xs = d_x[detektorit1 - d_det_per_ring * loop1 + d_det_per_ring * ll];
		*xd = d_x[detektorit1 - d_det_per_ring * loop1 + d_det_per_ring * ll];
		*ys = d_y[detektorit1 - d_det_per_ring * loop1 + d_det_per_ring * ll];
		*yd = d_y[detektorit1 - d_det_per_ring * loop1 + d_det_per_ring * ll];
#else
		*xs = d_x[detektorit1 - d_det_per_ring * loop1];
		*xd = d_x[detektorit2 - d_det_per_ring * loop2];
		*ys = d_y[detektorit1 - d_det_per_ring * loop1];
		*yd = d_y[detektorit2 - d_det_per_ring * loop2];
#endif
}
#endif

#if !defined(N_RAYS) && (!defined(RAW) || defined(LISTMODE))
// Get the detector coordinates for the current sinogram bin
__device__ void get_detector_coordinates(const unsigned int* d_xyindex, const unsigned short* d_zindex, const unsigned int d_size_x, const unsigned int idx,
	const unsigned short d_TotSinos, float* xs, float* xd, float* ys, float* yd, float* zs, float* zd, const float* d_x, const float* d_y,
	const float* d_zdet, const unsigned long long int cumsum) {

#ifdef LISTMODE
	* xs = d_x[idx + cumsum];
	*xd = d_x[idx + d_size_x + cumsum];
	*ys = d_y[idx + cumsum];
	*yd = d_y[idx + d_size_x + cumsum];
	*zs = d_zdet[idx + cumsum];
	*zd = d_zdet[idx + d_size_x + cumsum];
#else
	if (d_xyindex[idx] >= d_size_x) {
		*xs = d_x[d_xyindex[idx]];
		*xd = d_x[d_xyindex[idx] - d_size_x];
		*ys = d_y[d_xyindex[idx]];
		*yd = d_y[d_xyindex[idx] - d_size_x];
	}
	else {
		*xs = d_x[d_xyindex[idx]];
		*xd = d_x[d_xyindex[idx] + d_size_x];
		*ys = d_y[d_xyindex[idx]];
		*yd = d_y[d_xyindex[idx] + d_size_x];
	}
	if (d_zindex[idx] >= d_TotSinos) {
		*zs = d_zdet[d_zindex[idx]];
		*zd = d_zdet[d_zindex[idx] - d_TotSinos];
	}
	else {
		*zs = d_zdet[d_zindex[idx]];
		*zd = d_zdet[d_zindex[idx] + d_TotSinos];
	}
#endif
}
#endif

#if defined(N_RAYS) && !defined(RAW)
// Get the detector coordinates for the current sinogram bin
__device__ void get_detector_coordinates_multiray(const unsigned int* d_xyindex, const unsigned short* d_zindex, const unsigned int d_size_x, const unsigned int idx,
	const unsigned short d_TotSinos, float* xs, float* xd, float* ys, float* yd, float* zs, float* zd, const float* d_x, const float* d_y,
	const float* d_zdet, const unsigned short lor, const float cr_pz) {

	if (N_RAYS3D > 1) {
		unsigned short rays3D = N_RAYS3D;
		if (N_RAYS3D % 2 == 0)
			rays3D++;
		int r[N_RAYS3D + 1];
		int hh = 0;
		for (int kk = -N_RAYS3D / 2; kk <= N_RAYS3D / 2; kk++) {
			if (kk == 0 && (N_RAYS3D % 2 == 0)) {}
			else {
				r[hh] = kk;
				hh++;
			}
		}
		const float rr = (float)(r[(lor - 1) % N_RAYS3D]);
		*zs = d_zdet[d_zindex[idx]] + cr_pz * rr;
		if (d_zindex[idx] >= d_TotSinos) {
			*zd = d_zdet[d_zindex[idx] - d_TotSinos] + cr_pz * rr;
		}
		else {
			*zd = d_zdet[d_zindex[idx] + d_TotSinos] + cr_pz * rr;
		}
	}
	else {
		*zs = d_zdet[d_zindex[idx]];
		if (d_zindex[idx] >= d_TotSinos) {
			*zd = d_zdet[d_zindex[idx] - d_TotSinos];
		}
		else {
			*zd = d_zdet[d_zindex[idx] + d_TotSinos];
		}
	}
#ifdef N_RAYS2D
		const unsigned short ll = (lor - 1) / N_RAYS3D * 2;
		*xs = d_x[d_xyindex[idx] + d_size_x * ll];
		*ys = d_y[d_xyindex[idx] + d_size_x * ll];
		if (d_xyindex[idx] >= d_size_x) {
			*xd = d_x[d_xyindex[idx] + d_size_x * (ll - 1)];
			*yd = d_y[d_xyindex[idx] + d_size_x * (ll - 1)];
		}
		else {
			*xd = d_x[d_xyindex[idx] + d_size_x * (ll + 1)];
			*yd = d_y[d_xyindex[idx] + d_size_x * (ll + 1)];
		}
#else
		*xs = d_x[d_xyindex[idx]];
		*ys = d_y[d_xyindex[idx]];
		if (d_xyindex[idx] >= d_size_x) {
			*xd = d_x[d_xyindex[idx] - d_size_x];
			*yd = d_y[d_xyindex[idx] - d_size_x];
		}
		else {
			*xd = d_x[d_xyindex[idx] + d_size_x];
			*yd = d_y[d_xyindex[idx] + d_size_x];
		}
#endif
}
#endif

#ifdef FIND_LORS
// Get the detector coordinates for the current measurement (precomputation phase)
__device__ void get_detector_coordinates_precomp(const unsigned int d_size_x, const unsigned int idx, const unsigned short d_TotSinos, float* xs, float* xd, float* ys, float* yd, float* zs,
	float* zd, const float* d_x, const float* d_y, const float* d_zdet) {

	const unsigned int id = idx % d_size_x;
	const unsigned int idz = idx / d_size_x;
	*xs = d_x[id];
	*xd = d_x[id + d_size_x];
	*ys = d_y[id];
	*yd = d_y[id + d_size_x];
	*zs = d_zdet[idz];
	*zd = d_zdet[idz + d_TotSinos];
}
#endif
#endif

// Compute the voxel index where the current perpendicular measurement starts
__device__ int perpendicular_start(const float d_b, const float d, const float d_d, const unsigned int d_N) {
	int tempi = 0;
	float start = d_b - d + d_d;
	for (unsigned int ii = 0u; ii < d_N; ii++) {
		if (start > 0.f) {
			tempi = (int)(ii);
			break;
		}
		start += d_d;
	}
	return tempi;
}

// Compute the probability for the perpendicular elements
__device__ void perpendicular_elements(const float d_b, const float d_d1, const unsigned int d_N1, const float d, const float d_d2, const unsigned int d_N2,
	const float* d_atten, float* templ_ijk, unsigned int* tempk, const unsigned int z_loop, const unsigned int d_N, const unsigned int d_NN,
	const float* d_norm, const unsigned int idx, const float global_factor) {
	int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	//float start = d_b - d + d_d1;
	//// Find the closest y-index value by finding the smallest y-distance between detector 2 and all the y-pixel coordinates
	//for (int ii = 0; ii < d_N1; ii++) {
	//	if (start > 0.f) {
	//		apu = ii;
	//		break;
	//	}
	//	start += d_d1;
	//}
	//*templ_ijk = d_d2;
	* tempk = (unsigned int)(apu) * d_N + z_loop * d_N1 * d_N2;
#ifdef CT
	* templ_ijk = d_d2;
#else
	float temp = d_d2 * (float)(d_N2);
	// Probability
	temp = 1.f / temp;
#ifdef ATN
		float jelppi = 0.f;
		for (unsigned int iii = 0u; iii < d_N2; iii++) {
			jelppi += (*templ_ijk * (-d_atten[*tempk + iii * d_NN]));
		}
		temp *= expf(jelppi);
#endif
#ifdef NORM
		temp *= d_norm[idx];
#endif
	temp *= global_factor;
	*templ_ijk = temp * d_d2;
#endif
}

#ifdef N_RAYS
// Compute the probability for the perpendicular elements
__device__ float perpendicular_elements_multiray(const float d_b, const float d_d1, const unsigned int d_N1, const float d, const float d_d2, const unsigned int d_N2,
	const float* d_atten, unsigned int* tempk, const unsigned int z_loop, const unsigned int d_N, const unsigned int d_NN, float* jelppi) {
	int apu = perpendicular_start(d_b, d, d_d1, d_N1);

	* tempk = (unsigned int)(apu) * d_N + z_loop * d_N1 * d_N2;

	return d_d2 * (float)(d_N2);
}
#endif

// Compute functions (9) and (29) (detector larger than source)
__device__ void d_g_s(const float tmin, const float t_min, unsigned int* v_min, float* t_0, int* v_u, const float diff, const float b, const float d, const float s) {

	if (tmin == t_min)
		// (11)
		*v_min = 1u;
	else {
		// (2) and (19)
		const float p_t = s + tmin * (diff);
		// (12)
		*v_min = __float2uint_ru((p_t - b) / d);
	}
	//if (tmax == t_max)
	//	// (13)
	//	*v_max = N;
	//else {
	//	// (2) and (19)
	//	p_t = s + tmax * (diff);
	//	// (14)
	//	*v_max = __float2int_rz((p_t - b) / d);
	//}
	// (9)
	* t_0 += (((float)(*v_min) * d) / (diff));
	//  (29)
	*v_u = 1;
}

// Compute functions (9) and (29) (source larger than detector)
__device__ void s_g_d(const float tmin, const float t_min, unsigned int* v_max, float* t_0, int* v_u, const float diff, const float b, const float d, const float s,
	const unsigned int N) {

	if (tmin == t_min)
		// (15)
		*v_max = N - 1u;
	else {
		// (2) and (19)
		const float p_t = s + tmin * (diff);
		// (16)
		*v_max = __float2uint_rd((p_t - b) / d);
	}
	//if (tmax == t_max)
	//	// (17)
	//	*v_min = 0;
	//else {
	//	// (2) and (19)
	//	p_t = s + tmax * (diff);
	//	// (18)
	//	*v_min = (int)_rtp((p_t - b) / d);
	//}
	// (9)
	* t_0 += (((float)(*v_max) * d) / (diff));
	// (29)
	*v_u = -1;
}

// same as above, but for precomputation phase
__device__ void d_g_s_precomp(const float tmin, const float t_min, const float tmax, const float t_max, unsigned int* v_min, unsigned int* v_max, float* t_0, int* v_u,
	const float diff, const float b, const float d, const float s, const unsigned int N) {

	if (tmin == t_min)
		// (11)
		*v_min = 1u;
	else {
		// (2) and (19)
		const float p_t = s + tmin * (diff);
		// (12)
		*v_min = __float2uint_ru((p_t - b) / d);
	}
	if (tmax == t_max)
		// (13)
		*v_max = N;
	else {
		// (2) and (19)
		const float p_t = s + tmax * (diff);
		// (14)
		*v_max = __float2uint_rd((p_t - b) / d);
	}
	// (9)
	*t_0 += (((float)(*v_min) * d) / (diff));
	//  (29)
	*v_u = 1;
}

// same as above, but for precomputation phase
__device__ void s_g_d_precomp(const float tmin, const float t_min, const float tmax, const float t_max, unsigned int* v_min, unsigned int* v_max, float* t_0, int* v_u,
	const float diff, const float b, const float d, const float s, const unsigned int N) {

	if (tmin == t_min)
		// (15)
		*v_max = N - 1u;
	else {
		// (2) and (19)
		const float p_t = s + tmin * (diff);
		// (16)
		*v_max = __float2uint_rd((p_t - b) / d);
	}
	if (tmax == t_max)
		// (17)
		*v_min = 0u;
	else {
		// (2) and (19)
		const float p_t = s + tmax * (diff);
		// (18)
		*v_min = __float2uint_ru((p_t - b) / d);
	}
	// (9)
	*t_0 += (((float)(*v_max) * d) / (diff));
	// (29)
	*v_u = -1;
}

// Compute the index of the current voxel
__device__ unsigned int compute_ind(const int tempj, const int tempi, const int tempk, const unsigned int d_N1, const unsigned int d_N2, const unsigned int d_N, const unsigned int d_Nx,
	const unsigned int d_Nyx) {
	unsigned int local_ind = (unsigned int)(tempj) * d_Nx + (unsigned int)(tempi) + (unsigned int)(tempk) * d_Nyx;
	//if (local_ind >= d_N) {
	//	if (local_ind - d_N1 >= d_N)
	//		local_ind -= (d_N1 * d_N2);
	//	else if (local_ind - 1u >= d_N)
	//		local_ind -= d_N1;
	//	else
	//		local_ind--;
	//}
	//else if (local_ind < 0) {
	//	if (local_ind + d_N1 < 0)
	//		local_ind += (d_N1 * d_N2);
	//	else if (local_ind + 1 < 0)
	//		local_ind += d_N1;
	//	else
	//		local_ind++;
	//}
	return local_ind;
}

#ifdef ATN
__device__ float compute_matrix_element(const float t0, const float tc, const float L) {
	return (t0 - tc) * L;
}

__device__ void compute_attenuation(float* tc, float* jelppi, const float LL, const float t0, const int tempi, const int tempj, const int tempk, const unsigned int Nx,
	const unsigned int Nyx, const float* d_atten) {
	*jelppi += (compute_matrix_element(t0, *tc, LL) * -d_atten[tempi + tempj * Nx + Nyx * tempk]);
	*tc = t0;
}
#endif

#ifdef SIDDON
// compute the distance that the ray traverses in the current voxel
__device__ float compute_element(float* t0, float* tc, const float L, const float tu, const int u, int* temp_ijk, float* temp) {
	float local_ele = (*t0 - *tc) * L;
	*temp_ijk += u;
	*tc = *t0;
	*t0 += tu;
#ifndef CT
	* temp += local_ele;
#endif
	return local_ele;
}

// compute the probability of emission in the current voxel
__device__ float compute_element_2nd(float* t0, float* tc, const float L, const float tu, const int u, int* temp_ijk, const float temp) {
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
__device__ int voxel_index(const float pt, const float diff, const float d, const float apu) {
	return __float2int_rz((pt * diff - apu) / d);
}

__device__ bool siddon_pre_loop_2D(const float b1, const float b2, const float diff1, const float diff2, const float max1, const float max2,
	const float d1, const float d2, const unsigned int N1, const unsigned int N2, int* temp1, int* temp2, float* t1u, float* t2u, unsigned int* Np,
	const int TYYPPI, const float ys, const float xs, const float yd, const float xd, float* tc, int* u1, int* u2, float* t10, float* t20) {
	// If neither x- nor y-directions are perpendicular
// Correspond to the equations (9) and (10) from reference [2]
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

	unsigned int imin, imax, jmin, jmax;

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

	if (*tc == *t10 || *tc == *t20)
		*tc -= 1e-7f;

	return false;
}

__device__ bool siddon_pre_loop_3D(const float bx, const float by, const float bz, const float x_diff, const float y_diff, const float z_diff,
	const float maxxx, const float maxyy, const float bzb, const float dx, const float dy, const float dz,
	const unsigned int Nx, const unsigned int Ny, const unsigned int Nz, int* tempi, int* tempj, int* tempk, float* tyu, float* txu, float* tzu,
	unsigned int* Np, const int TYYPPI, const float ys, const float xs, const float yd, const float xd, const float zs, const float zd, float* tc,
	int* iu, int* ju, int* ku, float* tx0, float* ty0, float* tz0) {

	const float apu_tx = bx - xs;
	const float apu_ty = by - ys;
	const float apu_tz = bz - zs;
	*tx0 = (apu_tx) / (x_diff);
	*ty0 = (apu_ty) / (y_diff);
	*tz0 = (apu_tz) / (z_diff);
	const float txback = (maxxx - xs) / (x_diff);
	const float tyback = (maxyy - ys) / (y_diff);
	const float tzback = (bzb - zs) / (z_diff);

	const float txmin = fmin(*tx0, txback);
	const float txmax = fmax(*tx0, txback);
	const float tymin = fmin(*ty0, tyback);
	const float tymax = fmax(*ty0, tyback);
	const float tzmin = fmin(*tz0, tzback);
	const float tzmax = fmax(*tz0, tzback);

	*tc = fmax(fmax(txmin, tzmin), tymin);
	const float tmax = fmin(fmin(txmax, tzmax), tymax);

	unsigned int imin, imax, jmin, jmax, kmin, kmax;

	if (TYYPPI == 0) {
		if (*tc >= tmax) {
			return true;
		}

		if (xs < xd)
			d_g_s_precomp(*tc, txmin, tmax, txmax, &imin, &imax, tx0, iu, x_diff, bx, dx, xs, Nx);
		else
			s_g_d_precomp(*tc, txmin, tmax, txmax, &imin, &imax, tx0, iu, x_diff, bx, dx, xs, Nx);

		if (ys < yd)
			d_g_s_precomp(*tc, tymin, tmax, tymax, &jmin, &jmax, ty0, ju, y_diff, by, dy, ys, Ny);
		else
			s_g_d_precomp(*tc, tymin, tmax, tymax, &jmin, &jmax, ty0, ju, y_diff, by, dy, ys, Ny);

		if (zs < zd)
			d_g_s_precomp(*tc, tzmin, tmax, tzmax, &kmin, &kmax, tz0, ku, z_diff, bz, dz, zs, Nz);
		else
			s_g_d_precomp(*tc, tzmin, tmax, tzmax, &kmin, &kmax, tz0, ku, z_diff, bz, dz, zs, Nz);

		*Np = (kmax - kmin + 1) + (jmax - jmin + 1) + (imax - imin + 1);
	}
	else {
		if (xs < xd)
			d_g_s(*tc, txmin, &imin, tx0, iu, x_diff, bx, dx, xs);
		else
			s_g_d(*tc, txmin, &imax, tx0, iu, x_diff, bx, dx, xs, Nx);

		if (ys < yd)
			d_g_s(*tc, tymin, &jmin, ty0, ju, y_diff, by, dy, ys);
		else
			s_g_d(*tc, tymin, &jmax, ty0, ju, y_diff, by, dy, ys, Ny);

		if (zs < zd)
			d_g_s(*tc, tzmin, &kmin, tz0, ku, z_diff, bz, dz, zs);
		else
			s_g_d(*tc, tzmin, &kmax, tz0, ku, z_diff, bz, dz, zs, Nz);
	}

	const float pt = ((fmin(fmin(*tz0, *ty0), *tx0) + *tc) / 2.f);

	*tempi = voxel_index(pt, x_diff, dx, apu_tx);
	*tempj = voxel_index(pt, y_diff, dy, apu_ty);
	*tempk = voxel_index(pt, z_diff, dz, apu_tz);

	if (TYYPPI == 0) {
		if (*tempi < 0 || *tempj < 0 || *tempk < 0 || *tempi >= Nx || *tempj >= Ny || *tempk >= Nz)
			return true;
	}

	*txu = dx / fabs(x_diff);
	*tyu = dy / fabs(y_diff);
	*tzu = dz / fabs(z_diff);

	return false;
}

#ifdef TOF
#define _2PI 0.3989423f

__device__ float normPDF(const float x, const float mu, const float sigma) {

	const float a = (x - mu) / sigma;

	return _2PI / sigma * expf(-0.5f * a * a);
}

__device__ void TOFDis(const float x_diff, const float y_diff, const float z_diff, const float tc, const float LL, float* D, float* DD) {
	const float xI = x_diff * tc;
	const float yI = y_diff * tc;
	const float zI = z_diff * tc;
	*D = sqrtf(xI * xI + yI * yI + zI * zI) - LL / 2.f;
	*DD = *D;
}

__device__ float TOFWeight(const float element, const float sigma_x, const float D, const float DD, const float TOFCenter, float dX) {
	float output = normPDF(D, TOFCenter, sigma_x);
	dX = copysignf(dX, DD);
	for (long long int tr = 1LL; tr < (long long int)(TRAPZ_BINS) - 1LL; tr++)
		output += (normPDF(D - dX * (float)tr, TOFCenter, sigma_x) * 2.f);
	output += normPDF(D - copysignf(element, DD), TOFCenter, sigma_x);
	//return (element * (normPDF(D, TOFCenter, sigma_x) + normPDF(D - element * sign(DD), TOFCenter, sigma_x)) / 2.f);
	return output;
}


__device__ float TOFLoop(const float DD, const float element, float* TOFVal, const float* TOFCenter,
	const float sigma_x, float* D, const unsigned int tid, const float epps) {
	float TOFSum = 0.f;
#ifdef DEC
	float apu[NBINS];
#endif
	const float dX = element / (TRAPZ_BINS - 1.f);
#pragma unroll NBINS
	for (long long int to = 0L; to < NBINS; to++) {
#ifdef DEC
		apu[to] = TOFWeight(element, sigma_x, *D, DD, TOFCenter[to], dX) * dX;
		TOFSum += apu[to];
#else
		const float apu = TOFWeight(element, sigma_x, *D, DD, TOFCenter[to], dX) * dX;
		TOFSum += apu;
#endif
	}
#ifdef DEC
	* D -= (copysignf(element, DD));
#endif
	if (TOFSum < epps)
		TOFSum = epps;
#ifdef DEC
#pragma unroll NBINS
	for (long long int to = 0L; to < NBINS; to++)
		TOFVal[to + tid] = apu[to] / TOFSum;
#endif
	return TOFSum;
}


__device__ void denominatorTOF(float* ax, const float element, const float* d_OSEM, unsigned int local_ind, const float TOFSum, float* TOFVal,
	const float DD, const float* TOFCenter, const float sigma_x, float* D, const unsigned int tid, const float epps, const unsigned int d_N) {
#if !defined(MBSREM)
	unsigned int ll = NBINS;
#pragma unroll N_REKOS
	for (unsigned int kk = 0U; kk < N_REKOS; kk++) {
		unsigned int ii = ll * kk;
#else
	const unsigned int ii = 0U;
#endif
	float apu = element * d_OSEM[local_ind];
#ifndef DEC
	const float dX = element / (TRAPZ_BINS - 1.f);
#endif
#pragma unroll NBINS
	for (long long int to = 0L; to < NBINS; to++) {
#ifdef DEC
		ax[to + ii] += (apu * TOFVal[to + tid]);
#else
		const float jelppi = (TOFWeight(element, sigma_x, *D, DD, TOFCenter[to], dX) * dX) / TOFSum;
		ax[to + ii] += (apu * jelppi);
#endif
	}
#if !defined(MBSREM)
	local_ind += d_N;
	}
#endif
#ifndef DEC
* D -= (copysignf(element, DD));
#endif
}



// Nominator (y for backprojection)
__device__ void nominatorTOF(const unsigned char* MethodList, float* ax, const float* d_Sino, const float d_epsilon_mramla, const float d_epps,
	const float temp, const float* d_sc_ra, const unsigned int idx, const long long int TOFSize, const float local_sino) {
	float local_rand = 0.f;
#ifdef RANDOMS
	local_rand = d_sc_ra[idx];
#endif
#if !defined(MBSREM)
	unsigned int ll = NBINS;
#pragma unroll N_REKOS
	for (unsigned int kk = 0; kk < N_REKOS; kk++) {
		unsigned int ii = ll * kk;
#else
	const unsigned int ii = 0U;
#endif
#pragma unroll NBINS
	for (long long int to = 0L; to < NBINS; to++) {
		ax[to + ii] *= temp;
		if (ax[to + ii] < d_epps)
			ax[to + ii] = d_epps;
#ifdef RANDOMS
		ax[to + ii] += local_rand;
#endif
#if defined(MRAMLA) && !defined(MBSREM)
		if (MethodList[kk] != 1u)
			ax[to + ii] = d_Sino[idx + to * TOFSize] / ax[to + ii];
		else if (MethodList[kk] == 1u) { // MRAMLA/MBSREM
			if (ax[to + ii] <= d_epsilon_mramla && local_rand == 0.f && local_sino > 0.f)
				ax[to + ii] = d_Sino[idx + to * TOFSize] / d_epsilon_mramla - (d_Sino[idx + to * TOFSize] / (d_epsilon_mramla * d_epsilon_mramla)) * (ax[to + ii] - d_epsilon_mramla);
			else
				ax[to + ii] = d_Sino[idx + to * TOFSize] / ax[to + ii];
		}
#else
		ax[to + ii] = d_Sino[idx + to * TOFSize] / ax[to + ii];
#endif
	}
#if !defined(MBSREM)
	}
#endif
}


__device__ void backprojectTOF(const unsigned int local_ind, const float local_ele, const unsigned int tid, const float* TOFVal, const float* yax,
	CAST* d_Summ
#ifndef DEC
	, const float temp, const float sigma_x, float* D, const float DD, const float* TOFCenter, const float epps, const float TOFSum
#endif
#ifdef MBSREM
	, const RecMethodsOpenCL MethodListOpenCL, const unsigned int d_alku, const unsigned char MBSREM_prepass, float* minimi, float* axACOSEM, const float* d_OSEM,
	float* d_E, CAST* d_co, CAST* d_aco, const float local_sino, const unsigned int idx, const long long int TOFSize
#else
	, CAST* d_rhs, const unsigned char no_norm, const unsigned int d_N
#endif
) {
	unsigned int yy = local_ind;
	float val = 0.f;
#if !defined(MBSREM)
	unsigned int ll = NBINS;
#pragma unroll N_REKOS
	for (unsigned int kk = 0; kk < N_REKOS; kk++) {
		unsigned int ii = ll * kk;
#else
	const unsigned int ii = 0U;
#endif
#ifndef DEC
	const float dX = element / (TRAPZ_BINS - 1.f);
#endif

	float yaxTOF = 0.f;
#pragma unroll NBINS
	for (long long int to = 0L; to < NBINS; to++) {

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
			axACOSEM[to] += (apu * d_OSEM[local_ind]);
#endif

		val += apu;
		yaxTOF += apu * yax[to + ii];
	}

#ifdef MBSREM
	if (d_alku == 0u) {
		if (MBSREM_prepass == 1)
#ifdef ATOMIC
			atomicAdd(&d_Summ[local_ind], __float2ull_rn(val * TH));
#else
			atomicAdd(&d_Summ[local_ind], val);
#endif
		if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino != 0.f)
#ifdef ATOMIC
			atomicAdd(&d_co[local_ind], __float2ull_rn(yaxTOF * TH));
#else
			atomicAdd(&d_co[local_ind], yaxTOF);
#endif
		if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino != 0.f)
#ifdef ATOMIC
			atomicAdd(&d_aco[local_ind], __float2ull_rn(yaxTOF * TH));
#else
			atomicAdd(&d_aco[local_ind], yaxTOF);
#endif
	}
#else
	if (no_norm == 0u && ii == 0)
#ifdef ATOMIC
		atomicAdd(&d_Summ[local_ind], __float2ull_rn(val * TH));
#else
		atomicAdd(&d_Summ[local_ind], val);
#endif
#ifdef ATOMIC
	atomicAdd(&d_rhs[yy], __float2ull_rn(yaxTOF * TH));
#else
	atomicAdd(&d_rhs[yy], (yaxTOF));
#endif
#endif

#if !defined(MBSREM)
	yy += d_N;
	}
#endif

#ifndef DEC
* D -= (copysignf(local_ele / temp, DD));
#endif
}


__device__ void sensTOF(const unsigned int local_ind, const float local_ele, const unsigned int tid, const float* TOFVal, CAST* d_Summ,
#ifndef DEC
	const float temp, const float sigma_x, float* D, const float DD, const float* TOFCenter, const float epps, const float TOFSum,
#endif
#ifdef MBSREM
	const RecMethodsOpenCL MethodListOpenCL, const unsigned int d_alku, const unsigned char MBSREM_prepass, float* minimi, float* axACOSEM, const float* d_OSEM,
	float* d_E, const unsigned int idx, const long long int TOFSize,
#endif
	const unsigned char no_norm) {
	float val = 0.f;
#ifndef DEC
	const float dX = element / (TRAPZ_BINS - 1.f);
#endif
#pragma unroll NBINS
	for (long long int to = 0L; to < NBINS; to++) {

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
			axACOSEM[to] += (apu * d_OSEM[local_ind]);
#endif
		val += apu;
	}
#ifdef MBSREM
	if (no_norm == 0u && d_alku == 0u)
#ifdef ATOMIC
		atomicAdd(&d_Summ[local_ind], __float2ull_rn(val * TH));
#else
		atomicAdd(&d_Summ[local_ind], val);
#endif
#else
	if (no_norm == 0u)
#ifdef ATOMIC
		atomicAdd(&d_Summ[local_ind], __float2ull_rn(val * TH));
#else
		atomicAdd(&d_Summ[local_ind], val);
#endif
#endif

#ifndef DEC
	* D -= (copysignf(local_ele / temp, DD));
#endif
}
#endif