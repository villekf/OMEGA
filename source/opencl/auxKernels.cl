
/*******************************************************************************************************************************************
* This file contains "auxliary" kernels. What this means is that this file contains kernels for many of the regularization techniques/priors
* as well as kernels for some algorithm computations. The latter contains also convolution, element-wise multiplication and division, 
* derivatives and other functions. This file uses the preprocessor definitions and functions from general_opencl_functions.h. Note that
* the inclusion is not done here but rather during the compilation.
*
* Copyright (C) 2019-2026 Ville-Veikko Wettenhovi
*
* This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************************************************************************/

#ifndef LOCAL_SIZE3
#define LOCAL_SIZE3 1
#endif
#ifndef LTYPE
#define LTYPE int
#endif
#ifndef LTYPE3
#define LTYPE3 int3
#endif

// The total size of the local/shared memory region
#if defined(NLM_) || defined(PROXNLM) || defined(NLM4_)
#define SIZEX LOCAL_SIZE + SWINDOWX * 2 + PWINDOWX * 2
#define SIZEY LOCAL_SIZE2 + SWINDOWY * 2 + PWINDOWY * 2
#define SIZEZ LOCAL_SIZE3 + SWINDOWZ * 2 + PWINDOWZ * 2
#ifndef NLTYPE
#define NLTYPE 0
#endif
#endif
#if defined(GGMRF) || defined(HYPER) || defined(RDPCORNERS)
#define SIZEX LOCAL_SIZE + SWINDOWX * 2
#define SIZEY LOCAL_SIZE2 + SWINDOWY * 2
#define SIZEZ LOCAL_SIZE3 + SWINDOWZ * 2
#endif
#ifdef MEDIAN
#define KOKO (SEARCH_WINDOW_X * 2 + 1) * (SEARCH_WINDOW_Y * 2 + 1) * (SEARCH_WINDOW_Z * 2 + 1)
#endif
#ifndef AF // START NOTAF


// Combine multi-GPU backprojections and sensitivity images
// No longer used as of v2.0
KERN
void summa(const CLGLOBAL CAST* d_Summ_device, CLGLOBAL CAST* d_Summ_local, const CLGLOBAL CAST* d_rhs_device, CLGLOBAL CAST* d_rhs_local,
	const uint im_dim, const uchar no_norm) {

	uint gid = GID0;

	for (uint i = gid; i < im_dim; i += GSIZE0) {
		if (no_norm == 0u)
			d_Summ_local[i] += d_Summ_device[i];
		d_rhs_local[i] += d_rhs_device[i];
	}
}

// MLEM/OSEM steps before backprojection
// Used only by implementation 3
KERN
void forward(CLGLOBAL float* d_outputFP, const CLGLOBAL float* CLRESTRICT meas
#ifdef CT
	, CLGLOBAL float* d_outputCT
#endif
#ifdef RANDOMS
	, const CLGLOBAL float* CLRESTRICT rand
#endif
	) {
	uint gid = GID0;
#ifdef CT
	const float apu = EXP(-d_outputFP[gid]);
	d_outputFP[gid] = apu;
#ifdef RANDOMS
	d_outputCT[gid] = (meas[gid] * apu) / (apu + rand[gid]);
#endif
#else
#ifdef RANDOMS
	d_outputFP[gid] = meas[gid] / (d_outputFP[gid] + rand[gid] + 1e-6f);
#else
	d_outputFP[gid] = meas[gid] / (d_outputFP[gid] + 1e-6f);
#endif
#endif
}


// Compute MLEM/OSEM image estimate
// Used only by implementation 3
KERNEL3
void computeEstimate(const CLGLOBAL CAST* CLRESTRICT d_Summ, const CLGLOBAL CAST* CLRESTRICT d_rhs, CLGLOBAL float* d_im, const float d_epps, const int3 d_N, const uchar no_norm
#ifdef CT
	, const float flat
#endif
	) {

	int3 ind = CMINT3(GID0, GID1, GID2);
	size_t idx = GID0 + GID1 * GSIZE0 + GID2 * GSIZE0 * GSIZE1;
#if defined(CUDA) || defined(HIP)
	if (ind.x >= d_N.x || ind.y >= d_N.y || ind.z >= d_N.z)
#else
	if (any(ind >= d_N))
#endif
		return;
	float apu = d_im[idx];
#ifdef CT
	apu *= flat;
#endif
#if defined(ATOMIC) || defined(ATOMIC32) // START ATOMIC/ATOMIC32
	d_im[idx] = apu / (convert_float(d_Summ[idx]) / TH + d_epps) * (convert_float(d_rhs[idx]) / TH + d_epps);
#else
	d_im[idx] = apu / (d_Summ[idx] + d_epps) * (d_rhs[idx] + d_epps);
#endif
}

// PSF blurring
// This is mainly for non-float inputs
#ifdef PSF // START PSF
KERN
void Convolution3D(const CLGLOBAL CAST* input, CLGLOBAL CAST* output,
	CONSTANT float* convolution_window, int window_size_x, int window_size_y, int window_size_z) {
	int4 ind = CMINT4(GID0, GID1, GID2, 0);
	int4 ind_uus = CMINT4(0, 0, 0, 0);
	const uint Nyx = GSIZE0 * GSIZE1;
	float result = FLOAT_ZERO;
	int c = 0;
	for (int k = -window_size_z; k <= window_size_z; k++) {
		if (ind.z < window_size_z) {
			if (k < -ind.z)
				ind_uus.z = abs(k) - 1 - ind.z;
			else
				ind_uus.z = k + ind.z;
		}
		else {
			ind_uus.z = ind.z + k;
			if (ind_uus.z >= GSIZE2)
				ind_uus.z = GSIZE2 - 1 - (ind_uus.z - GSIZE2);
		}
		ind_uus.z *= Nyx;
		for (int j = -window_size_y; j <= window_size_y; j++) {
			if (ind.y < window_size_y) {
				if (j < -ind.y)
					ind_uus.y = abs(j) - 1 - ind.y;
				else
					ind_uus.y = j + ind.y;
			}
			else {
				ind_uus.y = ind.y + j;
				if (ind_uus.y >= GSIZE1)
					ind_uus.y = GSIZE1 - 1 - (ind_uus.y - GSIZE1);
			}
			ind_uus.y *= GSIZE0;
			for (int i = (-window_size_x); i <= window_size_x; i++) {
				if (ind.x < window_size_x) {
					if (i < -ind.x)
						ind_uus.x = abs(i) - 1 - ind.x;
					else
						ind_uus.x = i + ind.x;
				}
				else {
					ind_uus.x = ind.x + i;
					if (ind_uus.x >= GSIZE0)
						ind_uus.x = GSIZE0 - 1 - (ind_uus.x - GSIZE0);
				}
				int indeksi = ind_uus.x + ind_uus.y + ind_uus.z;
#if defined(ATOMIC) || defined(ATOMIC32) // START ATOMIC/ATOMIC32
				float p = convert_float(input[indeksi]) / TH;
#else
				float p = input[indeksi];
#endif // END ATOMIC/ATOMIC32
				p *= convolution_window[c];
				result += p;
				c++;
			}
		}
	}
#ifdef ATOMIC // START ATOMIC
	output[ind.x + ind.y * GSIZE0 + ind.z * Nyx] = convert_long(result * TH);
#elif defined(ATOMIC32)
	output[ind.x + ind.y * GSIZE0 + ind.z * Nyx] = convert_int(result * TH);
#else
	output[ind.x + ind.y * GSIZE0 + ind.z * Nyx] = result;
#endif // END ATOMIC
}

// PSF blurring, floats
KERNEL3
void Convolution3D_f(const CLGLOBAL float* input, CLGLOBAL float* output,
	CONSTANT float* convolution_window, int window_size_x, int window_size_y, int window_size_z) {
	int4 ind = CMINT4(GID0, GID1, GID2, 0);
	int4 ind_uus = CMINT4(0, 0, 0, 0);
	float result = FLOAT_ZERO;
	const uint Nyx = GSIZE0 * GSIZE1;
	int c = 0;
	for (int k = -window_size_z; k <= window_size_z; k++) {
		if (ind.z < window_size_z) {
			if (k < -ind.z)
				ind_uus.z = abs(k) - 1 - ind.z;
			else
				ind_uus.z = k + ind.z;
		}
		else {
			ind_uus.z = ind.z + k;
			if (ind_uus.z >= GSIZE2)
				ind_uus.z = GSIZE2 - 1 - (ind_uus.z - GSIZE2);
		}
		ind_uus.z *= Nyx;
		for (int j = -window_size_y; j <= window_size_y; j++) {
			if (ind.y < window_size_y) {
				if (j < -ind.y)
					ind_uus.y = abs(j) - 1 - ind.y;
				else
					ind_uus.y = j + ind.y;
			}
			else {
				ind_uus.y = ind.y + j;
				if (ind_uus.y >= GSIZE1)
					ind_uus.y = GSIZE1 - 1 - (ind_uus.y - GSIZE1);
			}
			ind_uus.y *= GSIZE0;
			for (int i = (-window_size_x); i <= window_size_x; i++) {
				if (ind.x < window_size_x) {
					if (i < -ind.x)
						ind_uus.x = abs(i) - 1 - ind.x;
					else
						ind_uus.x = i + ind.x;
				}
				else {
					ind_uus.x = ind.x + i;
					if (ind_uus.x >= GSIZE0)
						ind_uus.x = GSIZE0 - 1 - (ind_uus.x - GSIZE0);
				}
				int indeksi = ind_uus.x + ind_uus.y + ind_uus.z;
				float p = input[indeksi];
				p *= convolution_window[c];
				result += p;
				c++;
			}
		}
	}
	output[ind.x + ind.y * GSIZE0 + ind.z * Nyx] = result;
}

// Division by the sensitivity image
// Used only by implementation 3
KERN
void vectorDiv(const CLGLOBAL float* input, CLGLOBAL float* output, const float epps) {
	uint id = GID0;
	output[id] = output[id] / (input[id] + epps);
}

// Elementwise multiplication
KERN
void vectorMult(const CLGLOBAL float* input, CLGLOBAL float* output) {
	uint id = GID0;
	output[id] *= input[id];
}
#endif // END PSF
#endif // END NOTAF

// Complex elementwise multiplication
// Used by the filtering
// This kernel assumes that the imaginary element is right after the real element, i.e. [real,imaginary,real,imaginary,...]
#if !defined(METAL)
KERN
void vectorElementMultiply(const CLGLOBAL float* CLRESTRICT input, CLGLOBAL float* output, const uchar D2) {
	const LTYPE3 xyz = MINT3(GID0, GID1, GID2);
	const LTYPE n = xyz.x + xyz.y * GSIZE0 + xyz.z * GSIZE0 * GSIZE1;
	float mult;
	if (D2)
		mult = input[xyz.x + xyz.y * GSIZE0];
	else
		mult = input[xyz.x];
	output[2 * n] *= mult;
	output[2 * n + 1] *= mult;
}

// Complex elementwise division
// Used by the filtering
// This kernel assumes that the imaginary element is right after the real element, i.e. [real,imaginary,real,imaginary,...]
KERN
void vectorElementDivision(const CLGLOBAL float* CLRESTRICT input, CLGLOBAL float* output) {
	const LTYPE3 xyz = MINT3(GID0, GID1, GID2);
	const LTYPE n = xyz.x + xyz.y * GSIZE0 + xyz.z * GSIZE0 * GSIZE1;
	float div = input[xyz.x];
	// Make sure there is no division by zero
	if (fabs(div) < 1e-12f)
		div = (div < FLOAT_ZERO) ? -1e-12f : 1e-12f;
	output[2 * n] /= div;
	output[2 * n + 1] /= div;
}
#endif

// Non-local means
#ifdef NLM_ // START NLM
#if defined(USEIMAGES) && defined(OPENCL)
CONSTANT sampler_t samplerNLM = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP_TO_EDGE;
#endif

KERNEL3
#ifdef USEIMAGES
void NLM(CLGLOBAL float* CLRESTRICT grad, IMAGE3D u, CONSTANT float* gaussian, 
#else
void NLM(CLGLOBAL float* CLRESTRICT grad, const CLGLOBAL float* CLRESTRICT u, CONSTANT float* gaussian, 
#endif
#ifdef PYTHON
	const int Nx, const int Ny, const int Nz, const int NOrigx, const int NOrigy, const int NOrigz, 
#else
	const int3 N, const int3 NOrig, 
#endif
	const float h, const float epps, const float beta
#if NLTYPE >= 3
	, const float gamma
#endif
#if NLTYPE == 6 // NLGGMRF
	, const float p, const float q, const float c
#endif
#if defined(NLMADAPTIVE)
	, const float s
#endif
// Reference image
#ifdef NLMREF // START NLMREF
#ifdef USEIMAGES
	, IMAGE3D u_ref
#else
	, const CLGLOBAL float* CLRESTRICT u_ref
#endif
#endif // END NLMREF
#ifdef MASKPRIOR
	, MASKBPTYPE maskBP
#endif
#ifdef EFOVZ // Compute only in the voxels of the actual FOV (when using extended FOV)
	, CONSTANT uchar* fovIndices
#endif
#ifdef LARGEDIM
	, const uint2 nOffset
#endif
) {
#ifdef PYTHON
	const int3 N = MINT3(Nx, Ny, Nz);
#endif
	LTYPE3 ii = MINT3(GID0, GID1, GID2);
	const LTYPE n = (ii.x) + (ii.y) * (N.x) + (ii.z) * (N.x * N.y);
	LTYPE startX = GRID0 * LSIZE0 - SWINDOWX - PWINDOWX + LID0;
	LTYPE startY = GRID1 * LSIZE1 - SWINDOWY - PWINDOWY + LID1;
	LTYPE startZ = GRID2 * LSIZE2 - SWINDOWZ - PWINDOWZ + LID2;
	LTYPE endX = (GRID0 + 1) * LSIZE0 + SWINDOWX + PWINDOWX;
	LTYPE endY = (GRID1 + 1) * LSIZE1 + SWINDOWY + PWINDOWY;
	LTYPE endZ = (GRID2 + 1) * LSIZE2 + SWINDOWZ + PWINDOWZ;
	LOCAL float lCache[NLM_TILEX * NLM_TILEY * NLM_TILEZ];
#ifdef NLMREF
	LOCAL float lCacheRef[NLM_TILEX * NLM_TILEY * NLM_TILEZ];
#endif
	LTYPE indZ = LID2;
	for (LTYPE zz = startZ; zz < endZ; zz += LSIZE2) {
		LTYPE indY = LID1;
		for (LTYPE yy = startY; yy < endY; yy += LSIZE1) {
			LTYPE indX = LID0;
			for (LTYPE xx = startX; xx < endX; xx += LSIZE0) {
#if defined(NLMREF) // START NLMREF
#ifdef USEIMAGES
#if defined(CUDA) || defined(HIP)
				lCacheRef[indX + indY * NLM_TILEX + indZ * NLM_TILEX * NLM_TILEY] = tex3D<float>(u_ref, xx, yy, zz);
#else
				lCacheRef[indX + indY * NLM_TILEX + indZ * NLM_TILEX * NLM_TILEY] = read_imagef(u_ref, samplerNLM, (int4)(xx, yy, zz, 0)).w;
#endif
#else
				if (xx < 0 || yy < 0 || zz < 0 || xx >= N.x || yy >= N.y || zz >= N.z)
					lCacheRef[indX + indY * NLM_TILEX + indZ * NLM_TILEX * NLM_TILEY] = FLOAT_ZERO;
				else
					lCacheRef[indX + indY * NLM_TILEX + indZ * NLM_TILEX * NLM_TILEY] = u_ref[(xx) + (yy) * N.x + (zz) * N.x * N.y];
#endif
#endif // END NLMREF
#ifdef USEIMAGES
#if defined(CUDA) || defined(HIP)
				lCache[indX + indY * NLM_TILEX + indZ * NLM_TILEX * NLM_TILEY] = tex3D<float>(u, xx, yy, zz);
#else
				lCache[indX + indY * NLM_TILEX + indZ * NLM_TILEX * NLM_TILEY] = read_imagef(u, samplerNLM, (int4)(xx, yy, zz, 0)).w;
#endif
#else
				if (xx < 0 || yy < 0 || zz < 0 || xx >= N.x || yy >= N.y || zz >= N.z)
					lCache[indX + indY * NLM_TILEX + indZ * NLM_TILEX * NLM_TILEY] = FLOAT_ZERO;
				else
					lCache[indX + indY * NLM_TILEX + indZ * NLM_TILEX * NLM_TILEY] = u[(xx) + (yy) * N.x + (zz) * N.x * N.y];
#endif
				indX += LSIZE0;
			}
			indY += LSIZE1;
		}
		indZ += LSIZE2;
	}
	BARRIER
#if defined(CUDA) || defined(HIP)
	if (ii.x >= N.x || ii.y >= N.y || ii.z >= N.z)
#else
	if (any(ii >= N))
#endif
		return;
#ifdef MASKPRIOR
	const int maskVal = readMaskBP(maskBP, CINT3(ii), CUINT3(N));
#ifndef MASKSCALE
    if (maskVal == 0)
        return;
#endif
#endif
	const int3 xxyyzz = CMINT3(LID0 + SWINDOWX + PWINDOWX, LID1 + SWINDOWY + PWINDOWY, LID2 + SWINDOWZ + PWINDOWZ);
#ifdef MASKSCALE
	const float uj = NLMFETCH(lCache, xxyyzz.x, xxyyzz.y, xxyyzz.z);
#endif
#if NLTYPE == 6
	// Precompute for NLGGMRF
	const float cpq = POWR(c, p - q);
#endif
	float output;
#ifdef MASKSCALE
	if (maskVal == 0) {
		float weight_sum = epps;
		output = FLOAT_ZERO;
#if NLTYPE == 1
		float outputAla = epps;
#endif
#if defined(NLMADAPTIVE)
		float hh = FLOAT_ZERO;
		const float pSize = CFLOAT((PWINDOWX * 2 + 1) * (PWINDOWY * 2 + 1) * (PWINDOWZ * 2 + 1));
#endif
#pragma unroll
		for (int i = -1; i <= 1; i++) {
#pragma unroll
			for (int j = -1; j <= 1; j++) {
				int k = 0;
				if (i == 0 && j == 0)
					continue;
				float weight = FLOAT_ZERO;
				float distance = FLOAT_ZERO;
				int pz = 0;
#pragma unroll
					for (int py = -1; py <= 1; py++) {
						int dim_g = (pz + PWINDOWZ) * (PWINDOWX * 2 + 1) * (PWINDOWY * 2 + 1) + (py + PWINDOWY) * (PWINDOWX * 2 + 1) + (PWINDOWX - 1);
#pragma unroll
						for (int px = -1; px <= 1; px++) {
							const float gg = gaussian[dim_g++];
#ifdef NLMREF
							const float Pk = NLMFETCH(lCacheRef, xxyyzz.x + i + px, xxyyzz.y + j + py, xxyyzz.z);
							const float Pj = NLMFETCH(lCacheRef, xxyyzz.x + px, xxyyzz.y + py, xxyyzz.z);
#else
							const float Pk = NLMFETCH(lCache, xxyyzz.x + i + px, xxyyzz.y + j + py, xxyyzz.z);
							const float Pj = NLMFETCH(lCache, xxyyzz.x + px, xxyyzz.y + py, xxyyzz.z);
#endif
							const float PP = Pj - Pk;
							distance += gg * PP * PP;
						}
					}
#if defined(NLMADAPTIVE)
				hh = distance / pSize;
				weight = EXP(-distance / (hh * h + s));
#else
 				weight = EXP(-distance / h);
#endif
 				weight_sum += weight;
				const float uk = NLMFETCH(lCache, xxyyzz.x + i, xxyyzz.y + j, xxyyzz.z);
 				// Different NLM regularization methods
				// NLTYPE 0 = MRF NLM
				// NLTYPE 1 = NLTV
				// NLTYPE 2 = NLM filtered (i.e. similar to MRP)
				// NLTYPE 3 = NLRD
				// NLTYPE 4 = NL Lange
				// NLTYPE 5 = NLM filtered with Lange
				// NLTYPE 6 = NLGGMRF
				// NLTYPE 7 = ?
#if NLTYPE == 2 || NLTYPE == 5 // START NLM NLTYPE
				// NLMRP
 				output += weight * uk;
#elif NLTYPE == 0
 				output += (weight * (uj - uk));
#elif NLTYPE == 3
				// NLRD
				const float u = (uj - uk);
#ifndef USEMAD // START FMAD
				const float divPow = (uj + uk + gamma * fabs(u) + epps);
				output += weight * u * (gamma * fabs(u) + uj + 3.f * uk + epps * epps) / (divPow * divPow); 
#else
				const float divPow = FMAD(gamma, fabs(u), uj + uk + epps);
				output += weight * u * (FMAD(gamma, fabs(u), uj + 3.f * uk + epps * epps)) / (divPow * divPow); 
#endif // END FMAD
#elif NLTYPE == 4
				// Lange
				const float u = (uj - uk);
				const float uabs = sign(u);
				output += weight * (uabs - uabs / (fabs(u) / gamma + FLOAT_ONE));
#elif NLTYPE == 6
				// NLGGMRF
				const float delta = uj - uk;
				const float dcpq = POWR(fabs(delta / c), p - q);
				const float deltapqc = FLOAT_ONE + dcpq;
				output += weight * (POWR(fabs(delta), p - FLOAT_ONE) / deltapqc) * (p - gamma * ((dcpq * cpq) / deltapqc)) * sign(delta);
#elif NLTYPE == 7
				const float u = (uk - uj);
				const float apu = (u * u + gamma * gamma);
// #ifndef USEMAD // START FMAD
				output += ((FLOAT_TWO * u * u * u) / (apu * apu) - FLOAT_TWO * (u / apu));
// #else
// 				output += ((FLOAT_TWO * u * u * u) / FMAD(apu, apu, -FLOAT_TWO * (u / apu)));
// #endif // END FMAD
#else
 				//NLTV
				const float apuU = uj - uk;
 				output += (weight * apuU);
 				outputAla += weight * apuU * apuU;
#endif // END NLM NLTYPE
				}
			}
		weight_sum = FLOAT_ONE / weight_sum;
		output *= weight_sum;
#if NLTYPE == 2 // START NLM NLTYPE
		output = uj - output;
#elif NLTYPE == 5
		// Lange with NLMRP
		output = uj - output;
		const float uabs = sign(output);
		output = (uabs - uabs / (fabs(output) / gamma + FLOAT_ONE));
#elif NLTYPE == 1
#ifndef USEMAD // START FMAD
		output /= SQRT(outputAla * weight_sum + epps);
#else
		output /= SQRT(FMAD(outputAla, weight_sum, epps));
#endif // END FMAD
#endif // END NLM NLTYPE
		}
	else {
#endif
		output = NLMGradient(lCache, gaussian, xxyyzz.x, xxyyzz.y, xxyyzz.z, h, epps
#if NLTYPE >= 3
			, gamma
#endif
#if NLTYPE == 6
			, p, q, c
#endif
#if defined(NLMADAPTIVE)
			, s
#endif
#ifdef NLMREF
			, lCacheRef
#endif
			);
#ifdef MASKSCALE
	}
#endif
#ifdef LARGEDIM
	if (ii.z >= nOffset.x && ii.z < nOffset.y)
		grad[n - N.x * N.y * nOffset.x] += beta * output;
#else
	grad[n] += beta * output;
#endif
}
#endif // END NLM

// Relative difference prior
#ifdef RDP // START RDP
#ifdef OPENCL
#ifdef USEIMAGES
CONSTANT sampler_t samplerRDP = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP_TO_EDGE;
#endif
#ifdef RDPCORNERS
__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, LOCAL_SIZE3)))
#else
__kernel __attribute__((vec_type_hint(float2))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, LOCAL_SIZE3)))
#endif
#else
extern "C" __global__
#endif
#ifdef USEIMAGES
void RDPKernel(CLGLOBAL float* CLRESTRICT grad, IMAGE3D u, 
#else
void RDPKernel(CLGLOBAL float* CLRESTRICT grad, const CLGLOBAL float* CLRESTRICT u, 
#endif
#ifdef PYTHON
	const int Nx, const int Ny, const int Nz, const int NOrigx, const int NOrigy, const int NOrigz, 
#else
	const int3 N, const int3 NOrig, 
#endif
	const float gamma, const float epps, const float beta
#ifdef MASKPRIOR
	, MASKBPTYPE maskBP
#endif
#ifdef EFOVZ
	, CONSTANT uchar* fovIndices
#endif
#ifdef RDPCORNERS
	, CONSTANT float* weight
#endif
#ifdef RDPREF
#ifdef USEIMAGES
	, IMAGE3D u_ref
#else
	, const CLGLOBAL float* CLRESTRICT u_ref
#endif
#endif
#ifdef LARGEDIM
	, const uint2 nOffset
#endif
) {
#ifdef PYTHON
	const int3 N = MINT3(Nx, Ny, Nz);
#endif
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#ifdef RDPCORNERS // START RDPCORNERS
	float output = FLOAT_ZERO;
	LTYPE startX = GRID0 * LSIZE0 - SWINDOWX + LID0;
	LTYPE startY = GRID1 * LSIZE1 - SWINDOWY + LID1;
	LTYPE startZ = GRID2 * LSIZE2 - SWINDOWZ + LID2;
	LTYPE endX = (GRID0 + 1) * LSIZE0 + SWINDOWX;
	LTYPE endY = (GRID1 + 1) * LSIZE1 + SWINDOWY;
	LTYPE endZ = (GRID2 + 1) * LSIZE2 + SWINDOWZ;
	LOCAL float lCache[SW_TILEX * SW_TILEY * SW_TILEZ];
#ifdef RDPREF
	LOCAL float lCacheRef[SW_TILEX * SW_TILEY * SW_TILEZ];
#endif
	LTYPE indZ = LID2;
	for (LTYPE zz = startZ; zz < endZ; zz += LSIZE2) {
		LTYPE indY = LID1;
		for (LTYPE yy = startY; yy < endY; yy += LSIZE1) {
			LTYPE indX = LID0;
			for (LTYPE xx = startX; xx < endX; xx += LSIZE0) {
#if defined(RDPREF) // START RDPREF
#ifdef USEIMAGES
#if defined(CUDA) || defined(HIP)
				lCacheRef[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = tex3D<float>(u_ref, xx, yy, zz);
#else
				lCacheRef[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = read_imagef(u_ref, samplerRDP, (int4)(xx, yy, zz, 0)).w;
#endif
#else
				if (xx < 0 || yy < 0 || zz < 0 || xx >= N.x || yy >= N.y || zz >= N.z)
					lCacheRef[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = FLOAT_ZERO;
				else
					lCacheRef[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = u_ref[(xx) + (yy) * N.x + (zz) * N.x * N.y];
#endif
#endif // END RDPREF
#ifdef USEIMAGES
#if defined(CUDA) || defined(HIP)
				lCache[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = tex3D<float>(u, xx, yy, zz);
#else
				lCache[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = read_imagef(u, samplerRDP, (int4)(xx, yy, zz, 0)).w;
#endif
#else
				if (xx < 0 || yy < 0 || zz < 0 || xx >= N.x || yy >= N.y || zz >= N.z)
					lCache[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = FLOAT_ZERO;
				else
					lCache[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = u[(xx) + (yy) * N.x + (zz) * N.x * N.y];
#endif
				indX += LSIZE0;
			}
			indY += LSIZE1;
		}
		indZ += LSIZE2;
	}
	BARRIER
#endif // END RDPCORNERS
#if defined(CUDA) || defined(HIP)
	if (xyz.x >= N.x || xyz.y >= N.y || xyz.z >= N.z)
#else
	if (any(xyz >= N))
#endif
		return;
#ifdef EFOVZ
	if (fovIndices[xyz.z] == 0)
        return;
#endif
#ifdef MASKPRIOR
	const int maskVal = readMaskBP(maskBP, CINT3(xyz), CUINT3(N));
    if (maskVal == 0)
        return;
#endif
	const LTYPE n = (xyz.x) + (xyz.y) * (N.x) + (xyz.z) * (N.x * N.y);
// #endif
#ifdef RDPCORNERS // START RDPCORNERS
	const int3 xxyyzz = CMINT3(LID0 + SWINDOWX, LID1 + SWINDOWY, LID2 + SWINDOWZ);
	// Shared search-window gradient (general_opencl_functions.h)
	output = priorGradientSW(lCache, weight, xxyyzz.x, xxyyzz.y, xxyyzz.z, epps, gamma
#ifdef RDPREF
		, lCacheRef
#endif
	);
#ifdef LARGEDIM
	if (xyz.z >= nOffset.x && xyz.z < nOffset.y)
		grad[n - N.x * N.y * nOffset.x] += beta * output;
#else
	grad[n] += beta * output;
#endif
#else
	const float output = RDPGradientNorm(u, xyz.x, xyz.y, xyz.z, gamma, epps
#ifndef USEIMAGES
		, N
#endif
	);
#ifdef LARGEDIM
	if (xyz.z >= nOffset.x && xyz.z < nOffset.y)
		grad[n - N.x * N.y * nOffset.x] += beta * output;
#else
	grad[n] += beta * output;
#endif
#endif // END RDPCORNERS
}
#endif // END RDP


// Generalized Gaussian Markov random field
#ifdef GGMRF // START GGMRF
#if defined(USEIMAGES) && defined(OPENCL)
CONSTANT sampler_t samplerNLM = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP_TO_EDGE;
#endif

KERNEL3
#ifdef USEIMAGES
void GGMRFKernel(CLGLOBAL float* CLRESTRICT grad, IMAGE3D u, 
#else
void GGMRFKernel(CLGLOBAL float* CLRESTRICT grad, const CLGLOBAL float* CLRESTRICT u, 
#endif
	CONSTANT float* weight, const int3 N, const float p, const float q, const float c, const float pqc, const float beta
#ifdef MASKPRIOR
	, MASKBPTYPE maskBP
#endif
#ifdef EFOVZ
	, CONSTANT uchar* fovIndices
#endif
#ifdef LARGEDIM
	, const uint2 nOffset
#endif
) {

	LTYPE3 ii = MINT3(GID0, GID1, GID2);
	const LTYPE n = (ii.x) + (ii.y) * (N.x) + (ii.z) * (N.x * N.y);
	float output = FLOAT_ZERO;
	LTYPE startX = GRID0 * LSIZE0 - SWINDOWX + LID0;
	LTYPE startY = GRID1 * LSIZE1 - SWINDOWY + LID1;
	LTYPE startZ = GRID2 * LSIZE2 - SWINDOWZ + LID2;
	LTYPE endX = (GRID0 + 1) * LSIZE0 + SWINDOWX;
	LTYPE endY = (GRID1 + 1) * LSIZE1 + SWINDOWY;
	LTYPE endZ = (GRID2 + 1) * LSIZE2 + SWINDOWZ;
	LOCAL float lCache[SW_TILEX * SW_TILEY * SW_TILEZ];
	LTYPE indZ = LID2;
	for (LTYPE zz = startZ; zz < endZ; zz += LSIZE2) {
		LTYPE indY = LID1;
		for (LTYPE yy = startY; yy < endY; yy += LSIZE1) {
			LTYPE indX = LID0;
			for (LTYPE xx = startX; xx < endX; xx += LSIZE0) {
#ifdef USEIMAGES
#if defined(CUDA) || defined(HIP)
				lCache[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = tex3D<float>(u, xx, yy, zz);
#else
				lCache[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = read_imagef(u, samplerNLM, (int4)(xx, yy, zz, 0)).w;
#endif
#else
				if (xx < 0 || yy < 0 || zz < 0 || xx >= N.x || yy >= N.y || zz >= N.z)
					lCache[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = FLOAT_ZERO;
				else
					lCache[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = u[(xx) + (yy) * N.x + (zz) * N.x * N.y];
#endif
				indX += LSIZE0;
			}
			indY += LSIZE1;
		}
		indZ += LSIZE2;
	}
	BARRIER
#if defined(CUDA) || defined(HIP)
	if (ii.x >= N.x || ii.y >= N.y || ii.z >= N.z)
#else
	if (any(ii >= N))
#endif
		return;
#ifdef EFOVZ
	if (fovIndices[ii.z] == 0)
        return;
#endif
#ifdef MASKPRIOR
	const int maskVal = readMaskBP(maskBP, CINT3(ii), CUINT3(N));
    if (maskVal == 0)
        return;
#endif
	const int3 xxyyzz = CMINT3(LID0 + SWINDOWX, LID1 + SWINDOWY, LID2 + SWINDOWZ);
	// Shared search-window gradient (general_opencl_functions.h)
	const float epps = FLOAT_ZERO;
	output = priorGradientSW(lCache, weight, xxyyzz.x, xxyyzz.y, xxyyzz.z, epps, p, q, c, pqc);
#ifdef LARGEDIM
	if (ii.z >= nOffset.x && ii.z < nOffset.y)
		grad[n - N.x * N.y * nOffset.x] += beta * output;
#else
	grad[n] += beta * output;
#endif
}
#endif

// Median root prior
#ifdef MEDIAN // START MEDIAN
KERNEL3
void medianFilter3D(const CLGLOBAL float* grad, CLGLOBAL float* output, const int3 N, const int3 NOrig
#ifdef MASKPRIOR
	, MASKBPTYPE maskBP
#endif
#ifdef EFOVZ
	, CONSTANT uchar* fovIndices
#endif
) {
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
	if (xyz.x >= N.x + SEARCH_WINDOW_X || xyz.y >= N.y + SEARCH_WINDOW_Y || xyz.z >= N.z + SEARCH_WINDOW_Z || xyz.x < SEARCH_WINDOW_X || xyz.y < SEARCH_WINDOW_Y || xyz.z < SEARCH_WINDOW_Z)
		return;
	const LTYPE n = (xyz.x - SEARCH_WINDOW_X) + (xyz.y - SEARCH_WINDOW_Y) * (N.x) + (xyz.z - SEARCH_WINDOW_Z) * (N.x * N.y);
#if defined(EFOVZ) || defined(MASKPRIOR)
	// MRP is computed as (im - grad) / (grad + epps), which can cause issues if grad = 0 as in the case of mask use
	// This guarantees that the voxel value is returned instead (grad == im --> im - grad = 0)
	const LTYPE nCenter = xyz.x + xyz.y * (N.x + SEARCH_WINDOW_X * 2) + xyz.z * (N.x + SEARCH_WINDOW_X * 2) * (N.y + SEARCH_WINDOW_Y * 2);
#endif
#ifdef EFOVZ
	if (fovIndices[xyz.z - SEARCH_WINDOW_Z] == 0) {
		output[n] = grad[nCenter];
		return;
	}
#endif
#ifdef MASKPRIOR
	const LTYPE3 xyzOrig = MINT3(xyz.x - SEARCH_WINDOW_X, xyz.y - SEARCH_WINDOW_Y, xyz.z - SEARCH_WINDOW_Z);
	const int maskVal = readMaskBP(maskBP, CINT3(xyzOrig), CUINT3(N));
	if (maskVal == 0) {
		output[n] = grad[nCenter];
		return;
	}
#endif
	float median[KOKO];
	int uu = 0;
	for (LTYPE x = -SEARCH_WINDOW_X; x <= SEARCH_WINDOW_X; x++) {
		for (LTYPE y = -SEARCH_WINDOW_Y; y <= SEARCH_WINDOW_Y; y++) {
			for (LTYPE z = -SEARCH_WINDOW_Z; z <= SEARCH_WINDOW_Z; z++) {
				LTYPE pikseli = (xyz.x + (x)) + (xyz.y + (y)) * (N.x + SEARCH_WINDOW_X * 2) + (xyz.z + (z)) * (N.x + SEARCH_WINDOW_X * 2) * (N.y + SEARCH_WINDOW_Y * 2);
				median[uu] = grad[pikseli];
				uu++;
			}
		}
	}
	for (int hh = 0; hh < KOKO; hh++) {
		int ind = 0;
		for (int ll = 0; ll < KOKO; ll++) {
			if (median[hh] > median[ll] || (median[hh] == median[ll] && hh < ll))
				ind++;
		}
		if (ind == KOKO / 2) {
			output[n] = median[hh];
			return;
		}
	}
}
#endif // END MEDIAN

// TODO: Add actual support for backward difference
#if defined(PROXTV) || defined(TVGRAD) // START PROXTV || TVGRAD
// Backward difference for X-axis
DEVICE void backwardDiffX(float* apuVal, const LTYPE3 xyz, const int3 N, const LTYPE xx, const CLGLOBAL float* im) {
		const LTYPE xh = ((xyz.x - 1) + xyz.y * N.x + xyz.z * N.x * N.y);
		if (xyz.x == 0)
			*apuVal += im[xx];
		else
			*apuVal += (im[xx] - im[xh]);
}

// Backward difference for Y-axis
DEVICE void backwardDiffY(float* apuVal, const LTYPE3 xyz, const int3 N, const LTYPE xx, const CLGLOBAL float* im) {
		const LTYPE xh = ((xyz.x) + (xyz.y - 1) * N.x + xyz.z * N.x * N.y);
		if (xyz.y == 0)
			*apuVal += im[xx];
		else
			*apuVal += (im[xx] - im[xh]);
}

// Backward difference for Z-axis
DEVICE void backwardDiffZ(float* apuVal, const LTYPE3 xyz, const int3 N, const LTYPE xx, const CLGLOBAL float* im) {
		const LTYPE xh = ((xyz.x) + (xyz.y) * N.x + (xyz.z - 1) * N.x * N.y);
		if (xyz.z == 0)
			*apuVal += im[xx];
		else
			*apuVal += (im[xx] - im[xh]);
}

// Backward difference for X-axis (not the current voxel)
DEVICE void backwardDiffX2(float* apuVal, const LTYPE3 xyz, const int3 N, const CLGLOBAL float* im, const float imApu) {
		const LTYPE xh = ((xyz.x - 1) + xyz.y * N.x + xyz.z * N.x * N.y);
		if (xyz.x == 0)
			*apuVal += imApu;
		else
			*apuVal += (imApu - im[xh]);
}

// Backward difference for Y-axis (not the current voxel)
DEVICE void backwardDiffY2(float* apuVal, const LTYPE3 xyz, const int3 N, const CLGLOBAL float* im, const float imApu) {
		const LTYPE xh = ((xyz.x) + (xyz.y - 1) * N.x + (xyz.z) * N.x * N.y);
		if (xyz.y == 0)
			*apuVal += imApu;
		else
			*apuVal += (imApu - im[xh]);
}

// Backward difference for Z-axis (not the current voxel)
DEVICE void backwardDiffZ2(float* apuVal, const LTYPE3 xyz, const int3 N, const CLGLOBAL float* im, const float imApu) {
		const LTYPE xh = ((xyz.x) + (xyz.y) * N.x + (xyz.z - 1) * N.x * N.y);
		if (xyz.z == 0)
			*apuVal += imApu;
		else
			*apuVal += (imApu - im[xh]);
}

// Forward difference for X-axis
DEVICE void forwardDiffX(float* apuVal, const LTYPE3 xyz, const int3 N, const LTYPE xx, const CLGLOBAL float* im) {
		const LTYPE xh = ((xyz.x + 1) + xyz.y * N.x + xyz.z * N.x * N.y);
		if (xyz.x == N.x - 1)
			*apuVal -= im[xx];
		else
			*apuVal += (im[xh] - im[xx]);
}

// Forward difference for Y-axis
DEVICE void forwardDiffY(float* apuVal, const LTYPE3 xyz, const int3 N, const LTYPE xx, const CLGLOBAL float* im) {
		const LTYPE xh = ((xyz.x) + (xyz.y + 1) * N.x + xyz.z * N.x * N.y);
		if (xyz.y == N.y - 1)
			*apuVal -= im[xx];
		else
			*apuVal += (im[xh] - im[xx]);
}

// Forward difference for Z-axis
DEVICE void forwardDiffZ(float* apuVal, const LTYPE3 xyz, const int3 N, const LTYPE xx, const CLGLOBAL float* im) {
		const LTYPE xh = ((xyz.x) + (xyz.y) * N.x + (xyz.z + 1) * N.x * N.y);
		if (xyz.z == N.z - 1)
			*apuVal -= im[xx];
		else
			*apuVal += (im[xh] - im[xx]);
}

// Forward difference for X-axis (not the current voxel)
DEVICE void forwardDiffX2(float* apuVal, const LTYPE3 xyz, const int3 N, const CLGLOBAL float* im, const float imApu) {
		const LTYPE xh = ((xyz.x + 1) + xyz.y * N.x + xyz.z * N.x * N.y);
		if (xyz.x == N.x - 1)
			*apuVal -= imApu;
		else
			*apuVal += (im[xh] - imApu);
}

// Forward difference for Y-axis (not the current voxel)
DEVICE void forwardDiffY2(float* apuVal, const LTYPE3 xyz, const int3 N, const CLGLOBAL float* im, const float imApu) {
		const LTYPE xh = ((xyz.x) + (xyz.y + 1) * N.x + (xyz.z) * N.x * N.y);
		if (xyz.y == N.y - 1)
			*apuVal -= imApu;
		else
			*apuVal += (im[xh] - imApu);
}

// Forward difference for Z-axis (not the current voxel)
DEVICE void forwardDiffZ2(float* apuVal, const LTYPE3 xyz, const int3 N, const CLGLOBAL float* im, const float imApu) {
		const LTYPE xh = ((xyz.x) + (xyz.y) * N.x + (xyz.z + 1) * N.x * N.y);
		if (xyz.z == N.z - 1)
			*apuVal -= imApu;
		else
			*apuVal += (im[xh] - imApu);
}
#endif // END PROXTV || TVGRAD


#ifdef PROXTV // START PROXTV
#ifndef DIFFTYPE 
#define DIFFTYPE 0
#endif
// Computing q of proximal TV (see http://dx.doi.org/10.1088/0031-9155/57/10/3065)
KERN
void ProxTVq(CLGLOBAL float* inputX, CLGLOBAL float* inputY, CLGLOBAL float* inputZ, const float alpha) {
	LTYPE idx = GID0;
	const float3 apu = MFLOAT3(inputX[idx], inputY[idx], inputZ[idx]);
#ifdef L2 // START L2
// L2 norm
	const float scale = fmax(FLOAT_ONE, length(apu) / alpha);
#else
// L1 norm
	const float scale = fmax(fmax(fabs(apu.z), fmax(fabs(apu.x), fabs(apu.y))) / alpha, FLOAT_ONE);
#endif // END L2
	inputX[idx] = apu.x / scale;
	inputY[idx] = apu.y / scale;
	inputZ[idx] = apu.z / scale;
}

#ifdef PROXTGV // START PROXTGV
KERN
// Same as above, but for TGV
// TGVZ refers to full 3D TGV, i.e. it takes into account x-, y- and z-axis voxels
#ifdef TGVZ
void ProxTGVq(CLGLOBAL float* inputX, CLGLOBAL float* inputY, CLGLOBAL float* inputZ, CLGLOBAL float* input2XY, CLGLOBAL float* input2XZ, CLGLOBAL float* input2YZ, const float alpha) {
#else
void ProxTGVq(CLGLOBAL float* inputX, CLGLOBAL float* inputY, CLGLOBAL float* input2XY, const float alpha) {
#endif
	LTYPE idx = GID0;
#ifdef TGVZ
	const float3 apu = MFLOAT3(inputX[idx], inputY[idx], inputZ[idx]);
	const float3 apu2 = MFLOAT3(input2XY[idx], input2XZ[idx], input2YZ[idx]);
#else
	const float3 apu = MFLOAT3(inputX[idx], inputY[idx], FLOAT_ZERO);
	const float3 apu2 = MFLOAT3(input2XY[idx], FLOAT_ZERO, FLOAT_ZERO);
#endif
#ifdef L2 // START L2
	const float scale = fmax(FLOAT_ONE, SQRT(apu.x * apu.x + apu.y * apu.y + apu.z * apu.z + (apu2.x * apu2.x) * FLOAT_TWO + (apu2.y * apu2.y) * FLOAT_TWO + (apu2.z * apu2.z) * FLOAT_TWO) / alpha);
#else
	const float scale = fmax(fmax(fabs(apu2.z),fmax(fabs(apu2.y), fmax(fabs(apu2.x), fmax(fabs(apu.z), fmax(fabs(apu.x), fabs(apu.y)))))) / alpha, FLOAT_ONE);
#endif // END L2
	inputX[idx] = apu.x / scale;
	inputY[idx] = apu.y / scale;
	input2XY[idx] = apu2.x / scale;
#ifdef TGVZ
	inputZ[idx] = apu.z / scale;
	input2XZ[idx] = apu2.y / scale;
	input2YZ[idx] = apu2.z / scale;
#endif
}
#endif // END PROXTGV

// Proximal TV divergence
KERNEL3
void ProxTVDivergence(const int3 N, const int3 NOrig, const CLGLOBAL float* CLRESTRICT gradX, const CLGLOBAL float* CLRESTRICT gradY, const CLGLOBAL float* CLRESTRICT gradZ, CLGLOBAL float* output
// The (optional) logical mask should be zero in regions where the prior is not needed
#ifdef MASKPRIOR
	, MASKBPTYPE maskBP
#endif
// The (optional) logical vector should be zero in axial slices where the extended FOV is
#ifdef EFOVZ
	, CONSTANT uchar* fovIndices
#endif
) {
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#if defined(CUDA) || defined(HIP)
	if (xyz.x >= N.x || xyz.y >= N.y || xyz.z >= N.z)
#else
	if (any(xyz >= N))
#endif
		return;
#ifdef EFOVZ
	if (fovIndices[xyz.z] == 0)
        return;
#endif
#ifdef MASKPRIOR
	const int maskVal = readMaskBP(maskBP, CINT3(xyz), CUINT3(N));
    if (maskVal == 0)
        return;
#endif
	const LTYPE x = xyz.x + xyz.y * N.x + xyz.z * N.x * N.y;
#if defined(EFOVZ)
	const LTYPE3 NDiff = (N - NOrig) / 2;
	xyz = xyz - NDiff;
	const LTYPE y = (xyz.x) + (xyz.y) * NOrig.x + (xyz.z) * NOrig.x * NOrig.y;
#else
	const LTYPE y = x;
#endif
	float apuVal = FLOAT_ZERO;
// Transpose of forward difference (backward difference)
#if DIFFTYPE == 0 // START DIFFTYPE == 0
		backwardDiffX(&apuVal, xyz, NOrig, y, gradX);
		backwardDiffY(&apuVal, xyz, NOrig, y, gradY);
		backwardDiffZ(&apuVal, xyz, NOrig, y, gradZ);
// Transpose of backward difference (forward difference)
#elif DIFFTYPE == 1  // START DIFFTYPE == 1
		forwardDiffX(&apuVal, xyz, NOrig, y, gradX);
		forwardDiffY(&apuVal, xyz, NOrig, y, gradY);
		forwardDiffZ(&apuVal, xyz, NOrig, y, gradZ);
#else
#endif // END DIFFTYPE
	output[x] -= apuVal;
}

// Proximal TV or TGV gradient computation
KERNEL3
void ProxTVGradient(const int3 N, const int3 NOrig, const CLGLOBAL float* CLRESTRICT im, CLGLOBAL float* outputX, CLGLOBAL float* outputY, CLGLOBAL float* outputZ, const float sigma2
#ifdef PROXTGV
#ifdef TGVZ
	, const CLGLOBAL float* CLRESTRICT vX, const CLGLOBAL float* CLRESTRICT vY, const CLGLOBAL float* CLRESTRICT vZ
#else
	, const CLGLOBAL float* CLRESTRICT vX, const CLGLOBAL float* CLRESTRICT vY
#endif
#endif
#ifdef MASKPRIOR
	, MASKBPTYPE maskBP
#endif
#ifdef EFOVZ
	, CONSTANT uchar* fovIndices
#endif
) {
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#if defined(CUDA) || defined(HIP)
	if (xyz.x >= N.x || xyz.y >= N.y || xyz.z >= N.z)
#else
	if (any(xyz >= N))
#endif
		return;
#ifdef EFOVZ
	if (fovIndices[xyz.z] == 0)
        return;
#endif
#ifdef MASKPRIOR
	const int maskVal = readMaskBP(maskBP, CINT3(xyz), CUINT3(N));
    if (maskVal == 0)
        return;
#endif
		
	const LTYPE x = xyz.x + xyz.y * N.x + xyz.z * N.x * N.y;
	float apuVal = FLOAT_ZERO;
	float imApu = im[x];
#if defined(EFOVZ)
	const LTYPE3 NDiff = (N - NOrig) / 2;
	const LTYPE y = (xyz.x - NDiff.x) + (xyz.y - NDiff.y) * NOrig.x + (xyz.z - NDiff.z) * NOrig.x * NOrig.y;
#else
	const LTYPE y = x;
#endif
////////////////////////
// Forward difference //
////////////////////////
#if DIFFTYPE == 0 // START DIFFTYPE == 0
	forwardDiffX2(&apuVal, xyz, N, im, imApu);
#ifdef PROXTGV
	apuVal -= vX[y];
#endif
	apuVal *=  sigma2;
	outputX[y] += apuVal;
	apuVal = FLOAT_ZERO;
	forwardDiffY2(&apuVal, xyz, N, im, imApu);
#ifdef PROXTGV
	apuVal -= vY[y];
#endif
	apuVal *=  sigma2;
	outputY[y] += apuVal;
	apuVal = FLOAT_ZERO;
	forwardDiffZ2(&apuVal, xyz, N, im, imApu);
#if defined(PROXTGV) && defined(TGVZ)
	apuVal -= vZ[y];
#endif
	apuVal *=  sigma2;
	outputZ[y] += apuVal;
/////////////////////////
// Backward difference //
/////////////////////////
#elif DIFFTYPE == 1 // START DIFFTYPE == 1
	backwardDiffX2(&apuVal, xyz, N, im, imApu);
#ifdef PROXTGV
	apuVal -= vX[y];
#endif
	apuVal *=  sigma2;
	outputX[y] += apuVal;
	apuVal = FLOAT_ZERO;
	backwardDiffY2(&apuVal, xyz, N, im, imApu);
#ifdef PROXTGV
	apuVal -= vY[y];
#endif
	apuVal *=  sigma2;
	outputY[y] += apuVal;
	apuVal = FLOAT_ZERO;
	backwardDiffZ2(&apuVal, xyz, N, im, imApu);
#if defined(PROXTGV) && defined(TGVZ)
	apuVal -= vZ[y];
#endif
	apuVal *=  sigma2;
	outputZ[y] += apuVal;
#else
#endif // END DIFFTYPE
}
#endif // END PROXTV

#ifdef PROXTGV // START PROXTGV
KERNEL3
// Symmetric derivative for TGV
void ProxTGVSymmDeriv(const int3 N, const int3 NOrig, const CLGLOBAL float* CLRESTRICT vX, const CLGLOBAL float* CLRESTRICT vY, 
#ifdef TGVZ
	const CLGLOBAL float* CLRESTRICT vZ, CLGLOBAL float* qX, CLGLOBAL float* qY, CLGLOBAL float* qZ, CLGLOBAL float* q2XY, CLGLOBAL float* q2XZ, CLGLOBAL float* q2YZ, 
#else
	CLGLOBAL float* qX, CLGLOBAL float* qY, CLGLOBAL float* q2XY, 
#endif
	const float sigma2
#ifdef MASKPRIOR
	, MASKBPTYPE maskBP
#endif
) {
	const LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#if defined(CUDA) || defined(HIP)
	if (xyz.x >= NOrig.x || xyz.y >= NOrig.y || xyz.z >= NOrig.z)
#else
	if (any(xyz >= NOrig))
#endif
		return;
#ifdef MASKPRIOR
	const LTYPE3 NDiff = (N - NOrig) / 2;
	const int maskVal = readMaskBP(maskBP, CINT3(xyz + NDiff), CUINT3(N));
    if (maskVal == 0)
        return;
#endif
		
	const LTYPE x = xyz.x + xyz.y * NOrig.x + xyz.z * NOrig.x * NOrig.y;
/////////////// X ///////////////
	float apuVal = FLOAT_ZERO;
	float imApuX = vX[x];
// Forward difference
#if DIFFTYPE == 0 // START DIFFTYPE == 0
	forwardDiffX2(&apuVal, xyz, NOrig, vX, imApuX);
	apuVal *= sigma2;
	qX[x] += apuVal;
/////////////// Y ///////////////
	apuVal = FLOAT_ZERO;
	float imApuY = vY[x];
	forwardDiffY2(&apuVal, xyz, NOrig, vY, imApuY);
	apuVal *= sigma2;
	qY[x] += apuVal;
#ifdef TGVZ
/////////////// Z ///////////////
	apuVal = FLOAT_ZERO;
	float imApuZ = vZ[x];
	forwardDiffZ2(&apuVal, xyz, NOrig, vZ, imApuZ);
	apuVal *= sigma2;
	qZ[x] += apuVal;
#endif
/////////////// XY/YX ///////////////
	apuVal = FLOAT_ZERO;
	forwardDiffY2(&apuVal, xyz, NOrig, vX, imApuX);
	forwardDiffX2(&apuVal, xyz, NOrig, vY, imApuY);
	apuVal *= sigma2 * FLOAT_HALF;
	q2XY[x] += apuVal;
#ifdef TGVZ
/////////////// XZ/ZX ///////////////
	apuVal = FLOAT_ZERO;
	forwardDiffZ2(&apuVal, xyz, NOrig, vX, imApuX);
	forwardDiffX2(&apuVal, xyz, NOrig, vZ, imApuZ);
	apuVal *= sigma2 * FLOAT_HALF;
	q2XZ[x] += apuVal;
/////////////// YZ/ZY ///////////////
	apuVal = FLOAT_ZERO;
	forwardDiffZ2(&apuVal, xyz, NOrig, vY, imApuY);
	forwardDiffY2(&apuVal, xyz, NOrig, vZ, imApuZ);
	apuVal *= sigma2 * FLOAT_HALF;
	q2YZ[x] += apuVal;
#endif
/////////////////////////
// Backward difference //
/////////////////////////
#elif DIFFTYPE == 1 // START DIFFTYPE == 1
	backwardDiffX2(&apuValX, xyz, NOrig, vX, imApuX);
	apuValX *=  sigma2;
	qX[x] += apuValX;
/////////////// Y ///////////////
	apuVal = FLOAT_ZERO;
	float imApuY = vY[x];
	backwardDiffY2(&apuValY, xyz, NOrig, vY, imApuY);
	apuValY *=  sigma2;
	qY[x] += apuValY;
#ifdef TGVZ
/////////////// Z ///////////////
	apuVal = FLOAT_ZERO;
	float imApuZ = vZ[x];
	backwardDiffZ2(&apuValZ, xyz, NOrig, vZ, imApuZ);
	apuValZ *=  sigma2;
	qZ[x] += apuValZ;
#endif
/////////////// XY/YX ///////////////
	apuVal = FLOAT_ZERO;
	backwardDiffY2(&apuVal, xyz, NOrig, vX, imApuX);
	backwardDiffX2(&apuVal, xyz, NOrig, vY, imApuY);
	apuVal *=  sigma2 * FLOAT_HALF;
	q2XY[x] += apuVal;
#ifdef TGVZ
/////////////// XZ/ZX ///////////////
	apuVal = FLOAT_ZERO;
	backwardDiffZ2(&apuVal, xyz, NOrig, vX, imApuX);
	backwardDiffX2(&apuVal, xyz, NOrig, vZ, imApuZ);
	apuVal *=  sigma2 * FLOAT_HALF;
	q2XZ[x] += apuVal;
/////////////// YZ/ZY ///////////////
	apuVal = FLOAT_ZERO;
	backwardDiffZ2(&apuVal, xyz, NOrig, vY, imApuY);
	backwardDiffY2(&apuVal, xyz, NOrig, vZ, imApuZ);
	apuVal *=  sigma2 * FLOAT_HALF;
	q2YZ[x] += apuVal;
#endif
// Central difference?
#else
#endif // END DIFFTYPE
}

// TGV divergence
KERNEL3
void ProxTGVDivergence(const int3 N, const int3 NOrig, const CLGLOBAL float* CLRESTRICT qX, const CLGLOBAL float* CLRESTRICT qY, 
#ifdef TGVZ
	const CLGLOBAL float* CLRESTRICT qZ, const CLGLOBAL float* CLRESTRICT q2XY, const CLGLOBAL float* CLRESTRICT q2XZ, const CLGLOBAL float* CLRESTRICT q2YZ, 
	CLGLOBAL float* vX, CLGLOBAL float* vY, CLGLOBAL float* vZ, 
#else
	const CLGLOBAL float* CLRESTRICT q2XY, CLGLOBAL float* vX, CLGLOBAL float* vY, 
#endif
	const CLGLOBAL float* CLRESTRICT pX, const CLGLOBAL float* CLRESTRICT pY, const CLGLOBAL float* CLRESTRICT pZ, const float theta, const float tau
#ifdef MASKPRIOR
	, MASKBPTYPE maskBP
#endif
) {
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#if defined(CUDA) || defined(HIP)
	if (xyz.x >= NOrig.x || xyz.y >= NOrig.y || xyz.z >= NOrig.z)
#else
	if (any(xyz >= NOrig))
#endif
		return;
#ifdef MASKPRIOR
	const LTYPE3 NDiff = (N - NOrig) / 2;
	const int maskVal = readMaskBP(maskBP, CINT3(xyz + NDiff), CUINT3(N));
    if (maskVal == 0)
        return;
#endif
	const LTYPE x = xyz.x + xyz.y * NOrig.x + xyz.z * NOrig.x * NOrig.y;
	float apuVal = FLOAT_ZERO;
// Transpose of forward difference
#if DIFFTYPE == 0 // START DIFFTYPE == 0
/////////////// X ///////////////
		backwardDiffX(&apuVal, xyz, NOrig, x, qX);
		backwardDiffY(&apuVal, xyz, NOrig, x, q2XY);
#ifdef TGVZ
		backwardDiffZ(&apuVal, xyz, NOrig, x, q2XZ);
#endif
		float vApu = vX[x];
#ifdef USEMAD
		float v1 = FMAD(tau, pX[x] + apuVal, vApu);
		vX[x] = FMAD(theta, v1 - vApu, v1);
#else
		float v1 = vApu + tau * FLOAT_ONE * (pX[x] + apuVal);
		vX[x] = v1 + theta * (v1 - vApu);
#endif
/////////////// Y ///////////////
		apuVal = FLOAT_ZERO;
		backwardDiffY(&apuVal, xyz, NOrig, x, qY);
		backwardDiffX(&apuVal, xyz, NOrig, x, q2XY);
#ifdef TGVZ
		backwardDiffZ(&apuVal, xyz, NOrig, x, q2YZ);
#endif
		vApu = vY[x];
#ifdef USEMAD
		v1 = FMAD(tau, pY[x] + apuVal, vApu);
		vY[x] = FMAD(theta, v1 - vApu, v1);
#else
		v1 = vApu + tau * FLOAT_ONE * (pY[x] + apuVal);
		vY[x] = v1 + theta * (v1 - vApu);
#endif
#ifdef TGVZ
/////////////// Z ///////////////
		apuVal = FLOAT_ZERO;
		backwardDiffZ(&apuVal, xyz, NOrig, x, qZ);
		backwardDiffX(&apuVal, xyz, NOrig, x, q2XZ);
		backwardDiffY(&apuVal, xyz, NOrig, x, q2YZ);
		vApu = vZ[x];
#ifdef USEMAD
		v1 = FMAD(tau, pZ[x] + apuVal, vApu);
		vZ[x] = FMAD(theta, v1 - vApu, v1);
#else
		v1 = vApu + tau * FLOAT_ONE * (pZ[x] + apuVal);
		vZ[x] = v1 + theta * (v1 - vApu);
#endif
#endif
// Transpose of backward difference
#elif DIFFTYPE == 1 // START DIFFTYPE == 1
/////////////// X ///////////////
		forwardDiffX(&apuVal, xyz, NOrig, x, qX);
		forwardDiffY(&apuVal, xyz, NOrig, x, q2XY);
#ifdef TGVZ
		forwardDiffZ(&apuVal, xyz, NOrig, x, q2XZ);
#endif
		float vApu = vX[x];
#ifdef USEMAD
		float v1 = FMAD(tau, pX[x] + apuVal, vApu);
		vX[x] = FMAD(theta, v1 - vApu, v1);
#else
		float v1 = vApu + tau * (pX[x] + apuVal);
		vX[x] = v1 + theta * (v1 - vApu);
#endif
/////////////// Y ///////////////
		apuVal = FLOAT_ZERO;
		forwardDiffY(&apuVal, xyz, NOrig, x, qY);
		forwardDiffX(&apuVal, xyz, NOrig, x, q2XY);
#ifdef TGVZ
		forwardDiffZ(&apuVal, xyz, NOrig, x, q2YZ);
#endif
		vApu = vY[x];
#ifdef USEMAD
		v1 = FMAD(tau, pY[x] + apuVal, vApu);
		vY[x] = FMAD(theta, v1 - vApu, v1);
#else
		v1 = vApu + tau * (pY[x] + apuVal);
		vY[x] = v1 + theta * (v1 - vApu);
#endif
#ifdef TGVZ
/////////////// Z ///////////////
		apuVal = FLOAT_ZERO;
		forwardDiffZ(&apuVal, xyz, NOrig, x, qZ);
		forwardDiffX(&apuVal, xyz, NOrig, x, q2XZ);
		forwardDiffY(&apuVal, xyz, NOrig, x, q2YZ);
		vApu = vZ[x];
#ifdef USEMAD
		v1 = FMAD(tau, pZ[x] + apuVal, vApu);
		vZ[x] = FMAD(theta, v1 - vApu, v1);
#else
		v1 = vApu + tau * (pZ[x] + apuVal);
		vZ[x] = v1 + theta * (v1 - vApu);
#endif
#endif
#else
#endif // END DIFFTYPE
}
#endif // END PROXTGV


// Gradient of hyperbolic prior
#if defined(HYPER) // START HYPER
#ifdef OPENCL
CONSTANT sampler_t samplerTV = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP_TO_EDGE;
#endif

#ifdef OPENCL
__kernel __attribute__((vec_type_hint(float2))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, LOCAL_SIZE3)))
#else
extern "C" __global__
#endif
#ifdef USEIMAGES
void hyperbolicKernel(CLGLOBAL float* CLRESTRICT grad, IMAGE3D u, 
#else
void hyperbolicKernel(CLGLOBAL float* CLRESTRICT grad, const CLGLOBAL float* CLRESTRICT u,
#endif
	const int3 N, const int3 NOrig, const float sigma, const float epps, const float beta, CONSTANT float* w
#ifdef MASKPRIOR
	, MASKBPTYPE maskBP
#endif
#ifdef EFOVZ
	, CONSTANT uchar* fovIndices
#endif
#ifdef LARGEDIM
	, const uint2 nOffset
#endif
) {
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#if defined(CUDA) || defined(HIP)
	if (xyz.x >= N.x || xyz.y >= N.y || xyz.z >= N.z)
#else
	if (any(xyz >= N))
#endif
		return;
#ifdef EFOVZ
	if (fovIndices[xyz.z] == 0)
        return;
#endif
#ifdef MASKPRIOR
	const int maskVal = readMaskBP(maskBP, CINT3(xyz), CUINT3(N));
    if (maskVal == 0)
        return;
#endif
	const LTYPE n = (xyz.x) + (xyz.y) * (N.x) + (xyz.z) * (N.x * N.y);
	float output = FLOAT_ZERO;
	LTYPE startX = GRID0 * LSIZE0 - SWINDOWX + LID0;
	LTYPE startY = GRID1 * LSIZE1 - SWINDOWY + LID1;
	LTYPE startZ = GRID2 * LSIZE2 - SWINDOWZ + LID2;
	LTYPE endX = (GRID0 + 1) * LSIZE0 + SWINDOWX;
	LTYPE endY = (GRID1 + 1) * LSIZE1 + SWINDOWY;
	LTYPE endZ = (GRID2 + 1) * LSIZE2 + SWINDOWZ;
	LOCAL float lCache[SW_TILEX * SW_TILEY * SW_TILEZ];
	LTYPE indZ = LID2;
	for (LTYPE zz = startZ; zz < endZ; zz += LSIZE2) {
		LTYPE indY = LID1;
		for (LTYPE yy = startY; yy < endY; yy += LSIZE1) {
			LTYPE indX = LID0;
			for (LTYPE xx = startX; xx < endX; xx += LSIZE0) {
#ifdef USEIMAGES
#if defined(CUDA) || defined(HIP)
				lCache[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = tex3D<float>(u, xx, yy, zz);
#else
				lCache[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = read_imagef(u, samplerTV, (int4)(xx, yy, zz, 0)).w;
#endif
#else
				if (xx < 0 || yy < 0 || zz < 0 || xx >= N.x || yy >= N.y || zz >= N.z)
					lCache[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = FLOAT_ZERO;
				else
					lCache[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = u[(xx) + (yy) * N.x + (zz) * N.x * N.y];
#endif
				indX += LSIZE0;
			}
			indY += LSIZE1;
		}
		indZ += LSIZE2;
	}
	BARRIER
	const int3 xxyyzz = CMINT3(LID0 + SWINDOWX, LID1 + SWINDOWY, LID2 + SWINDOWZ);
	// Shared search-window gradient (general_opencl_functions.h)
	output = priorGradientSW(lCache, w, xxyyzz.x, xxyyzz.y, xxyyzz.z, epps, sigma);
#ifdef LARGEDIM
	if (xyz.z >= nOffset.x && xyz.z < nOffset.y)
		grad[n - N.x * N.y * nOffset.x] += beta * output;
#else
	grad[n] += beta * output;
#endif
}
#endif

// Gradient of TV prior
// This is different from the proximal TV!
// Note that this contains several different TV methods
// SATV = Modified Lange prior
// JPTV = TV type 2
// ANATOMICAL1 = TV type 1, with anatomical reference image
// ANATOMICAL2 = TV type 2, with anatomical reference image
// ANATOMICAL3 = APLS
// TVW1 = Weighted TV
// Non-reference image TVs are identical (not counting Lange or weighted)
#if defined(TVGRAD) // START TVGRAD
// samplerTV and sqrtVal used to live here --> moved to general_opencl_functions.h

#ifdef OPENCL
__kernel __attribute__((vec_type_hint(float3))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, LOCAL_SIZE3)))
#else
extern "C" __global__
#endif
#ifdef USEIMAGES
void TVKernel(CLGLOBAL float* CLRESTRICT grad, IMAGE3D u, 
#else
void TVKernel(CLGLOBAL float* CLRESTRICT grad, const CLGLOBAL float* CLRESTRICT u,
#endif
#ifdef PYTHON
	const int Nx, const int Ny, const int Nz, const int NOrigx, const int NOrigy, const int NOrigz, 
#else
	const int3 N, const int3 NOrig, 
#endif
	const float sigma, const float epps, const float beta
#ifdef MASKPRIOR
	, MASKBPTYPE maskBP
#endif
#ifdef EFOVZ
	, CONSTANT uchar* fovIndices
#endif
#if defined(ANATOMICAL2) || defined(ANATOMICAL3)
	, const float C
#endif
#if defined(ANATOMICAL1) || defined(ANATOMICAL2) || defined(ANATOMICAL3)
	, CLGLOBAL float* CLRESTRICT S
#endif
#ifdef LARGEDIM
	, const uint2 nOffset
#endif
) {
#ifdef PYTHON
	const int3 N = MINT3(Nx, Ny, Nz);
#endif
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#if defined(CUDA) || defined(HIP)
	if (xyz.x >= N.x || xyz.y >= N.y || xyz.z >= N.z)
#else
	if (any(xyz >= N))
#endif
		return;
#ifdef EFOVZ
	if (fovIndices[xyz.z] == 0)
        return;
#endif
#ifdef MASKPRIOR
	const int maskVal = readMaskBP(maskBP, CINT3(xyz), CUINT3(N));
    if (maskVal == 0)
        return;
#endif
	const LTYPE n = (xyz.x) + (xyz.y) * (N.x) + (xyz.z) * (N.x * N.y);
	const float output = TVGradient(u, xyz.x, xyz.y, xyz.z, sigma, epps
#ifdef TVNEEDN
		, N
#endif
#if defined(ANATOMICAL2) || defined(ANATOMICAL3)
		, C
#endif
#if defined(ANATOMICAL1) || defined(ANATOMICAL2) || defined(ANATOMICAL3)
		, S
#endif
	);
#ifdef LARGEDIM
	if (xyz.z >= nOffset.x && xyz.z < nOffset.y)
		grad[n - N.x * N.y * nOffset.x] += beta * output;
#else
	grad[n] += beta * output;
#endif
}
#endif // END TVGRAD

// This kernel computes the image estimate for PKMA, MBSREM and BSREM
// I.e. the step after backprojection has been computed
#if defined(PKMA) || defined(MBSREM) || defined(BSREM)
KERNEL3
void PoissonUpdate(CLGLOBAL float* CLRESTRICT im, const CLGLOBAL float* CLRESTRICT rhs,
	const int3 N, const float lambda, const float epps, const float alpha, const uchar enforcePositivity) {
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#if defined(CUDA) || defined(HIP)
	if (xyz.x >= N.x || xyz.y >= N.y || xyz.z >= N.z)
#else
	if (any(xyz >= N))
#endif
		return;
	const LTYPE n = (xyz.x) + (xyz.y) * (N.x) + (xyz.z) * (N.x * N.y);
	// Use common function with fastPDHG, see general_opencl_functions.h
	im[n] = PoissonUpdateVoxel(im[n], rhs[n], lambda, epps, alpha, enforcePositivity);
}
#endif

// PDHG image update
// Similar to above, after backprojection, but for PDHG and its variants
// Different variations for subset and non-subset versions
#if defined(PDHG)
KERNEL3
void PDHGUpdate(
	CLGLOBAL float* CLRESTRICT im BUF0,
	const CLGLOBAL float* CLRESTRICT rhs BUF1,
	CLGLOBAL float* CLRESTRICT u BUF2,
#ifdef METAL
	SCALAR_PARAMS(scalarParams) BUF3,
	uint3 metalGlobalId [[thread_position_in_grid]]
#else
	const int3 N,
	const float epps,
	const float theta,
	const float tau,
	const uchar enforcePositivity
#endif
) {
#ifdef METAL
	UNPACK_SCALAR_PARAMS_PDHG(scalarParams)
	LTYPE3 xyz = MINT3(metalGlobalId.x, metalGlobalId.y, metalGlobalId.z);
#else
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#endif
#if defined(CUDA) || defined(HIP)
	if (xyz.x >= N.x || xyz.y >= N.y || xyz.z >= N.z)
#else
	if (ANY(xyz >= N))
#endif
		return;
	const LTYPE n = (xyz.x) + (xyz.y) * (N.x) + (xyz.z) * (N.x * N.y);
#ifdef SUBSETS
	// Use common function with fastPDHG, see general_opencl_functions.h
	im[n] = PDHGSubsetPrimal(im[n], rhs[n], tau, epps, enforcePositivity);
#else
	const float uPrev = u[n];
	float uNew = uPrev;
	uNew -= tau * rhs[n];
	if (enforcePositivity)
		uNew = FMAX(epps, uNew);
	u[n] = uNew;
	im[n] = uNew + theta * (uNew - uPrev);
#endif
}
#endif


#ifdef ROTATE
#if defined(USEIMAGES) && defined(OPENCL)
CONSTANT sampler_t samplerRotate = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_LINEAR | CLK_ADDRESS_CLAMP_TO_EDGE;
#elif defined(USEIMAGES) && defined(METAL)
constexpr metal::sampler samplerRotate(metal::coord::pixel, metal::filter::linear, metal::address::clamp_to_edge);
#endif
KERNEL void rotate(
	CLGLOBAL float* CLRESTRICT rotim BUF0,
	IMTYPE im TEX1,
#ifdef METAL
	SCALAR_PARAMS(scalarParams) BUF2,
	uint3 metalGlobalId [[thread_position_in_grid]]
#else
	const int Nx, const int Ny, const int Nz, const float cosa, const float sina
#endif
) {
	// Initial version from: https://stackoverflow.com/questions/9833316/cuda-image-rotation/10008412#10008412
#ifdef METAL
	UNPACK_SCALAR_PARAMS_ROTATE(scalarParams)
	LTYPE3 xyz = MINT3(metalGlobalId.x, metalGlobalId.y, metalGlobalId.z);
#else
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#endif
	if (xyz.x >= Nx || xyz.y >= Ny || xyz.z >= Nz)
		return;
	const LTYPE n = (xyz.x) + (xyz.y) * (Nx) + (xyz.z) * (Nx * Ny);

    const float xA = (float)(xyz.x - Nx/2) + FLOAT_HALF;
    const float yA = (float)(xyz.y - Ny/2) + FLOAT_HALF;

    const float src_x = (xA * cosa - yA * sina + Nx/2) - FLOAT_HALF;
    const float src_y = (xA * sina + yA * cosa + Ny/2) - FLOAT_HALF;

    if (src_x >= 0.0f && src_x < Nx && src_y >= 0.0f && src_y < Ny) {
		float val = FLOAT_ZERO;
#ifdef USEIMAGES
#if defined(CUDA) || defined(HIP)
        val = tex3D<float>(im, src_x + FLOAT_HALF, src_y + FLOAT_HALF, CFLOAT(xyz.z) + FLOAT_HALF);
#elif defined(OPENCL)
        val = read_imagef(im, samplerRotate, (float4)(src_x + FLOAT_HALF, src_y + FLOAT_HALF, CFLOAT(xyz.z) + FLOAT_HALF, FLOAT_ZERO)).w;
#elif defined(METAL)
		val = im.sample(samplerRotate, float3(src_x + FLOAT_HALF, src_y + FLOAT_HALF, CFLOAT(xyz.z) + FLOAT_HALF)).r;
#endif
#else
        // BILINEAR INTERPOLATION
        const int src_x0 = (int)(src_x);
        const int src_x1 = (src_x0 + 1);
        const int src_y0 = (int)(src_y);
        const int src_y1 = (src_y0 + 1);

        const float sx = (src_x - src_x0);
        const float sy = (src_y - src_y0);

        const int idx_src00 = MIN(MAX(0, src_x0 + src_y0 * Nx), (Nx * Ny) - 1);
        const int idx_src10 = MIN(MAX(0, src_x1 + src_y0 * Nx), (Nx * Ny) - 1);
        const int idx_src01 = MIN(MAX(0, src_x0 + src_y1 * Nx), (Nx * Ny) - 1);
        const int idx_src11 = MIN(MAX(0, src_x1 + src_y1 * Nx), (Nx * Ny) - 1);

        val  = (FLOAT_ONE - sx) * (FLOAT_ONE - sy) * im[idx_src00 + xyz.z * Nx * Ny];
        val += (       sx) * (FLOAT_ONE - sy) * im[idx_src10 + xyz.z * Nx * Ny];
        val += (FLOAT_ONE - sx) * (       sy) * im[idx_src01 + xyz.z * Nx * Ny];
        val += (       sx) * (       sy) * im[idx_src11 + xyz.z * Nx * Ny];
#endif
		rotim[n] = val;
    } 
	else {
        rotim[n] = 0.0f;
    }

}
#endif
