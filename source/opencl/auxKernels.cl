
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
#if defined(CUDA) || defined(HIP)
	if (ind.x >= d_N.x || ind.y >= d_N.y || ind.z >= d_N.z)
#else
	if (any(ind >= d_N))
#endif
		return;
	const size_t idx = (size_t)(ind.x) + (size_t)(ind.y) * (size_t)(d_N.x) + (size_t)(ind.z) * (size_t)(d_N.x) * (size_t)(d_N.y);
	float apu = d_im[idx];
#ifdef CT
	apu *= flat;
#endif
#if defined(ATOMIC) || defined(ATOMIC32) // START ATOMIC/ATOMIC32
#ifdef CT
	d_im[idx] = apu / (convert_float(d_Summ[idx]) / TH + d_epps) * (convert_float(d_rhs[idx]) / TH + d_epps);
#else
	// Sensitivity and rhs are floored at epps (matches subiterStep.h): voxels that the subset does not see stay unchanged
	d_im[idx] = apu / (FMAX(convert_float(d_Summ[idx]) / TH, d_epps) + d_epps) * (FMAX(convert_float(d_rhs[idx]) / TH, d_epps) + d_epps);
#endif
#else
#ifdef CT
	d_im[idx] = apu / (d_Summ[idx] + d_epps) * (d_rhs[idx] + d_epps);
#else
	d_im[idx] = apu / (FMAX(d_Summ[idx], d_epps) + d_epps) * (FMAX(d_rhs[idx], d_epps) + d_epps);
#endif
#endif
}

// PSF blurring
// This is mainly for non-float inputs
#ifdef PSF // START PSF
// A separate mirroring function rather than retyping the same operations for every dimension
// Slightly reduces code usage
DEVICE int mirrorIndex(int i, const int n) {
	while (i < 0 || i >= n) {
		if (i < 0)
			i = -i - 1;
		if (i >= n)
			i = 2 * n - 1 - i;
	}
	return i;
}
// Indexing was previously a bit messed up here, but this function hasn't been, and isn't, supported
KERN
void Convolution3D(const CLGLOBAL CAST* input, CLGLOBAL CAST* output,
	CONSTANT float* convolution_window, int window_size_x, int window_size_y, int window_size_z, const int3 N) {
	const int4 ind = CMINT4(GID0, GID1, GID2, 0);
	if (ind.x >= N.x || ind.y >= N.y || ind.z >= N.z)
		return;
	const LONG Nyx = CLONG_rtz(N.x) * CLONG_rtz(N.y);
	float result = FLOAT_ZERO;
	int c = 0;
	for (int k = -window_size_z; k <= window_size_z; k++) {
		const LONG indZ = CLONG_rtz(mirrorIndex(ind.z + k, N.z)) * Nyx;
		for (int j = -window_size_y; j <= window_size_y; j++) {
			const LONG indY = CLONG_rtz(mirrorIndex(ind.y + j, N.y)) * CLONG_rtz(N.x);
			for (int i = (-window_size_x); i <= window_size_x; i++) {
				const LONG indeksi = CLONG_rtz(mirrorIndex(ind.x + i, N.x)) + indY + indZ;
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
	const LONG outInd = CLONG_rtz(ind.x) + CLONG_rtz(ind.y) * CLONG_rtz(N.x) + CLONG_rtz(ind.z) * Nyx;
#ifdef ATOMIC // START ATOMIC
	output[outInd] = convert_long(result * TH);
#elif defined(ATOMIC32)
	output[outInd] = convert_int(result * TH);
#else
	output[outInd] = result;
#endif // END ATOMIC
}

// PSF blurring, floats
KERNEL3
void Convolution3D_f(const CLGLOBAL float* input, CLGLOBAL float* output,
	CONSTANT float* convolution_window, int window_size_x, int window_size_y, int window_size_z,
#ifdef PYTHON
	const int Nx, const int Ny, const int Nz) {
	const int3 N = MINT3(Nx, Ny, Nz);
#else
	const int3 N) {
#endif
	const int4 ind = CMINT4(GID0, GID1, GID2, 0);
	// See Convolution3D above
	if (ind.x >= N.x || ind.y >= N.y || ind.z >= N.z)
		return;
	float result = FLOAT_ZERO;
	const LONG Nyx = CLONG_rtz(N.x) * CLONG_rtz(N.y);
	int c = 0;
	for (int k = -window_size_z; k <= window_size_z; k++) {
		const LONG indZ = CLONG_rtz(mirrorIndex(ind.z + k, N.z)) * Nyx;
		for (int j = -window_size_y; j <= window_size_y; j++) {
			const LONG indY = CLONG_rtz(mirrorIndex(ind.y + j, N.y)) * CLONG_rtz(N.x);
			for (int i = (-window_size_x); i <= window_size_x; i++) {
				const LONG indeksi = CLONG_rtz(mirrorIndex(ind.x + i, N.x)) + indY + indZ;
				float p = input[indeksi];
				p *= convolution_window[c];
				result += p;
				c++;
			}
		}
	}
	output[CLONG_rtz(ind.x) + CLONG_rtz(ind.y) * CLONG_rtz(N.x) + CLONG_rtz(ind.z) * Nyx] = result;
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

#if defined(METAL)
KERNEL
void vectorElementMultiply(const CLGLOBAL float* CLRESTRICT input BUF0, CLGLOBAL float* output BUF1, constant uchar& D2 BUF2,
	uint3 metalGlobalId [[thread_position_in_grid]], uint3 metalGridSize [[threads_per_grid]]) {
	const LTYPE3 xyz = MINT3(metalGlobalId.x, metalGlobalId.y, metalGlobalId.z);
	const LTYPE gridX = metalGridSize.x;
	const LTYPE gridY = metalGridSize.y;
#else
KERN
void vectorElementMultiply(const CLGLOBAL float* CLRESTRICT input, CLGLOBAL float* output, const uchar D2) {
	const LTYPE3 xyz = MINT3(GID0, GID1, GID2);
	const LTYPE gridX = GSIZE0;
	const LTYPE gridY = GSIZE1;
#endif
	const LTYPE n = xyz.x + xyz.y * gridX + xyz.z * gridX * gridY;
	float mult;
	if (D2)
		mult = input[xyz.x + xyz.y * gridX];
	else
		mult = input[xyz.x];
	output[2 * n] *= mult;
	output[2 * n + 1] *= mult;
}

// Complex elementwise division
// Used by the filtering
// This kernel assumes that the imaginary element is right after the real element, i.e. [real,imaginary,real,imaginary,...]

#if defined(METAL)
KERNEL
void vectorElementDivision(const CLGLOBAL float* CLRESTRICT input BUF0, CLGLOBAL float* output BUF1,
	uint3 metalGlobalId [[thread_position_in_grid]], uint3 metalGridSize [[threads_per_grid]]) {
	const LTYPE3 xyz = MINT3(metalGlobalId.x, metalGlobalId.y, metalGlobalId.z);
	const LTYPE gridX = metalGridSize.x;
	const LTYPE gridY = metalGridSize.y;
#else
KERN
void vectorElementDivision(const CLGLOBAL float* CLRESTRICT input, CLGLOBAL float* output) {
	const LTYPE3 xyz = MINT3(GID0, GID1, GID2);
	const LTYPE gridX = GSIZE0;
	const LTYPE gridY = GSIZE1;
#endif
	const LTYPE n = xyz.x + xyz.y * gridX + xyz.z * gridX * gridY;
	float div = input[xyz.x];
	// Make sure there is no division by zero
	if (FABS(div) < 1e-12f)
		div = (div < FLOAT_ZERO) ? -1e-12f : 1e-12f;
	output[2 * n] /= div;
	output[2 * n + 1] /= div;
}

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
#ifdef NLMAXWEIGHT
		float maxWeight = FLOAT_ZERO;
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
#ifdef NLMAXWEIGHT
				maxWeight = FMAX(maxWeight, weight);
#endif
				const float uk = NLMFETCH(lCache, xxyyzz.x + i, xxyyzz.y + j, xxyyzz.z);
 				// Different NLM regularization methods
				// NLTYPE 0 = MRF NLM
				// NLTYPE 1 = NLTV
				// NLTYPE 2 = NLM filtered (i.e. similar to MRP)
				// NLTYPE 3 = NLRD
				// NLTYPE 4 = NL Lange
				// NLTYPE 5 = NLM filtered with Lange
				// NLTYPE 6 = NLGGMRF
				// NLTYPE 7 = NLGM
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
				const float delta = uj - uk;
				const float apu = (delta * delta + gamma * gamma);
				output += weight * (FLOAT_TWO * gamma * gamma * delta) / (apu * apu);
#else
 				//NLTV
				const float apuU = uj - uk;
 				output += (weight * apuU);
 				outputAla += weight * apuU * apuU;
#endif // END NLM NLTYPE
				}
			}
#ifdef NLMAXWEIGHT
#if NLTYPE == 2 || NLTYPE == 5
		output += maxWeight * uj;
#endif
		weight_sum += maxWeight;
#endif
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
void ProxTVq(CLGLOBAL float* inputX, CLGLOBAL float* inputY, CLGLOBAL float* inputZ, const float alpha, const LTYPE N) {
	LTYPE idx = GID0;
	if (idx >= N)
		return;
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
void ProxTGVq(CLGLOBAL float* inputX, CLGLOBAL float* inputY, CLGLOBAL float* inputZ, CLGLOBAL float* input2XY, CLGLOBAL float* input2XZ, CLGLOBAL float* input2YZ, const float alpha, const LTYPE N) {
#else
void ProxTGVq(CLGLOBAL float* inputX, CLGLOBAL float* inputY, CLGLOBAL float* input2XY, const float alpha, const LTYPE N) {
#endif
	LTYPE idx = GID0;
	if (idx >= N)
		return;
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

#ifndef METAL
// Convolves a volume of N0 x N1 x N2 with a separate 1D filter for each slice along dimension 2,
// along dimension 0 when axis is 0 and along dimension 1 when axis is 1. filt stores 2 * fRadius + 1
// taps for each of the N2 slices. The filters are symmetric, so correlation and convolution are the same.
// This is used by the SPECT rotation-based projector, where the collimator-detector response is a
// separable Gaussian whose width is different for every image plane. Splitting the response into two
// 1D passes, instead of one 2D convolution, is what makes a depth-dependent response affordable.
KERNEL void convolveCDRF(CLGLOBAL float* CLRESTRICT out, const CLGLOBAL float* CLRESTRICT in, const CLGLOBAL float* CLRESTRICT filt,
	const int N0, const int N1, const int N2, const int fRadius, const int axis) {
	const int x = GID0;
	const int y = GID1;
	const int z = GID2;
	if (x >= N0 || y >= N1 || z >= N2)
		return;
	const LTYPE ind = (LTYPE)x + (LTYPE)y * (LTYPE)N0 + (LTYPE)z * (LTYPE)N0 * (LTYPE)N1;
	const LTYPE fInd = (LTYPE)z * (LTYPE)(2 * fRadius + 1) + (LTYPE)fRadius;
	float summa = FLOAT_ZERO;
	if (axis == 0) {
		// The taps outside the volume are simply left out, which is the same as zero padding
		const int lo = MAX(-fRadius, -x);
		const int hi = MIN(fRadius, N0 - 1 - x);
		for (int kk = lo; kk <= hi; kk++)
			summa += filt[fInd + (LTYPE)kk] * in[ind + (LTYPE)kk];
	}
	else {
		const int lo = MAX(-fRadius, -y);
		const int hi = MIN(fRadius, N1 - 1 - y);
		for (int kk = lo; kk <= hi; kk++)
			summa += filt[fInd + (LTYPE)kk] * in[ind + (LTYPE)kk * (LTYPE)N0];
	}
	out[ind] = summa;
}

// The fused type-6 kernels below replace the per-view chain of ArrayFire operations (filter
// construction, reorders, slicing, convolutions, depth sum, detector resampling) by a few launches. They
// work directly on the layout of the rotation output, x (the depth, i.e. the image plane) fastest, then y
// (transaxial) and z (axial), so no reordering is needed. The collimator-detector response is built
// inside the kernels that use it, one group of planes at a time in shared memory (see cdrfAxial), or for the
// radii that do not fit there once per view and volume by cdrfWeights, with exactly the formulas of psfHelperType6:
//   depth of plane p (mm)   zmm = ((p - firstPlane) + (0.5 + firstPlane) - dist) * dxv
//   FWHM_axial              = sqrt(max(Az * zmm + Bz, 0)^2 + C^2) / dzv,  sigma = FWHM * sigmaCoeff
//   FWHM_trans              = sqrt(max(Axy * zmm + Bxy, 0)^2 + C^2) / dyv
//   sigma_trans             = sqrt(max(FWHM_trans^2 * sigmaCoeff^2 - comp, 1e-16 * sigmaCoeff^2)), comp the variance (voxels^2) that the
//                             interpolation of the rotation of the view adds transaxially (rotVarType6 in functions.hpp, per view and volume;
//                             a negative comp is the constant sigmaCoeff^2 of before, the formula below is then sqrt(max(FWHM_trans^2 - 1,
//                             1e-16)) * sigmaCoeff)
//   tap kk in [-r_p, r_p]   w(kk) = exp(-kk^2 / (2 sigma_p^2)), normalized by the sum over those 2 r_p + 1 taps,
//                           r_p = min(ceil(5.2565 sigma_p), r), the radius of plane p (per plane and direction)
// The planes behind the collimator (p < firstPlane, zmm <= 0) get all-zero weights. The weights of cdrfWeights
// are stored tap-major, tap kk of plane x is at (kk + r) * nPlanes + x (r the global radius, the stride is 2 r + 1), so
// that the threads of the convolution kernels, which are adjacent in x, read adjacent elements. Only the
// taps |kk| <= r_p of a plane are written and read, the radius r_p of every plane is stored in rAx / rTr. The
// arrays are zero padded at their edges, i.e. the taps outside of the array are left out of the convolution
// but not out of the normalization, which is what the filter bank convolution of convolveCDRF does as well.
// xOff is the absolute plane index of the first plane of the arrays
DEVICE float cdrfDepth(const int p, const int firstPlane, const float dist, const float dxv) {
	return ((CFLOAT(p - firstPlane) + (FLOAT_HALF + CFLOAT(firstPlane))) - dist) * dxv;
}

// The model of the geometric response of the collimator holes, CDRF_MODEL (a build option, the host passes -DCDRF_MODEL=n from
// the option CDRFModel). Every model gives, for a plane, the 1D response of a tap offset k (in voxels of the direction), which
// cdrfTap makes from the two numbers that cdrfSigma and cdrfHalfBase return for the plane: sigma and hb.
//   0 (default, 'gaussian'): the geometric response and the intrinsic resolution are one Gaussian of
//        FWHM = sqrt(Rg^2 + C^2): sigma is its sigma and hb is not used. w(k) = exp(-k^2 / (2 sigma^2))
//   1 ('square'): square holes. The response averaged over the hole positions factorises into x and y parts, each of them a
//        triangle with FWHM (the half base) Rg. The response is that triangle convolved with the intrinsic Gaussian
//        (FWHM C): sigma is the sigma of that Gaussian and hb the half base of the triangle (both in voxels), see cdrfSquareTap
// To add a model (e.g. hexagonal holes): give it a number here and in the host (the option CDRFModel in setMissingValues.m,
// SPECTParameters.m, proj.py, the cdrfModel6 scalar and the extent of the response in psfHelperType6 / computePadZType6),
// return its two parameters from cdrfSigma and cdrfHalfBase, its weight from cdrfTap, and its radius from cdrfRadius.
// The weights are point samples at the tap centres, normalized over the taps of the plane, whatever the model
#ifndef CDRF_MODEL
#define CDRF_MODEL 0
#endif
#if defined(CUDA) || defined(HIP)
#define CDRF_EXPF expf
#define CDRF_ERFCF erfcf
#else
#define CDRF_EXPF exp
#define CDRF_ERFCF erfc
#endif

#if CDRF_MODEL == 1
// P(x) = E[(S - x)+] for S ~ N(0, s^2), i.e. s * phi(x / s) - x * Q(x / s) with the density phi and the upper tail Q of the standard
// normal distribution. Its second derivative is the density of S
DEVICE float cdrfSquareP(const float x, const float s) {
	const float z = x / s;
	return s * 0.3989422804014327f * CDRF_EXPF(-0.5f * z * z) - x * 0.5f * CDRF_ERFCF(z * 0.7071067811865476f);
}

// The triangle t(u) = max(1 - |u| / R, 0) / R (FWHM R, the half base) convolved with the Gaussian of sigma s, at the offset k.
// t is the second difference of the ramp, t(u) = [(u + R)+ - 2 u+ + (u - R)+] / R^2, and the ramp convolved with the Gaussian is
// Q(x) = x + P(x), so  h(u) = [P(u + R) - 2 P(u) + P(u - R)] / R^2  with P(-x) = P(x) + x (the linear terms cancel in the
// difference, which keeps the tails of the result accurate): for u >= R all arguments are non-negative, for 0 <= u < R the
// last term P(u - R) = P(R - u) + (R - u) is the triangle itself plus the (small) smoothing. When R is a small fraction of s the
// difference of the closed form loses digits in float, so then the response is the Gaussian times the series
// sum_n m_2n / (2n)! He_2n(u / s) / s^2n, with the even moments m_2n = 2 R^2n / ((2n + 1) (2n + 2)) of the triangle and the
// Hermite polynomials He_2n, up to n = 4 (the error is below 2e-7 of the value for R < s / 2, where the series is used)
DEVICE float cdrfSquareTap(const int k, const float R, const float s) {
	const float u = CFLOAT(k < 0 ? -k : k);
	const float a = R / s;
	if (a < 0.5f) {
		const float z = u / s;
		const float z2 = z * z;
		if (z2 > 100.f)
			return FLOAT_ZERO;
		const float a2 = a * a;
		const float he2 = z2 - 1.f;
		const float he4 = (z2 - 6.f) * z2 + 3.f;
		const float he6 = ((z2 - 15.f) * z2 + 45.f) * z2 - 15.f;
		const float he8 = (((z2 - 28.f) * z2 + 210.f) * z2 - 420.f) * z2 + 105.f;
		const float f = 1.f + a2 * (he2 / 12.f + a2 * (he4 / 360.f + a2 * (he6 / 20160.f + a2 * he8 / 1814400.f)));
		return 0.3989422804014327f / s * CDRF_EXPF(-0.5f * z2) * f;
	}
	const float pp = cdrfSquareP(u + R, s);
	const float pu = cdrfSquareP(u, s);
	if (u >= R)
		return (pp - 2.f * pu + cdrfSquareP(u - R, s)) / (R * R);
	return (pp - 2.f * pu + cdrfSquareP(R - u, s) + (R - u)) / (R * R);
}

// The parameters of the square model of a plane, from the geometric half base Rg and the intrinsic FWHM C (mm) of the direction, dd
// the voxel size (mm): which = false the sigma of the intrinsic Gaussian and true the half base of the triangle (voxels). The
// transaxial pass is also applied to an image that the rotation has blurred by comp voxels^2 (variance, per view and volume, see
// rotVarType6 in functions.hpp; negative: the constant sigmaCoeff^2), which, as in the
// Gaussian model, is taken out of the response: out of the variance of the Gaussian (if it has that much, it is not allowed to
// be narrower than sMin, a delta), and what is missing out of the variance Rg^2 / 6 of the triangle
DEVICE float cdrfSquareParam(const bool trans, const float Rg, const float C, const float dd, const float sigmaCoeff, const float compIn, const bool which) {
	const float sMin = 1e-3f;
	const float Rv = Rg / dd;
	const float si = C / dd * sigmaCoeff;
	const float comp = trans ? (compIn < FLOAT_ZERO ? sigmaCoeff * sigmaCoeff : compIn) : FLOAT_ZERO;
	const float v = si * si - comp;
	if (v >= sMin * sMin)
		return which ? Rv : SQRT(v);
	if (!which)
		return sMin;
	return SQRT(FMAX(Rv * Rv - 6.f * (sMin * sMin - v), FLOAT_ZERO));
}
#endif

// The sigma (in voxels of the direction) of the response of plane p, trans = false for the axial direction (A = Az,
// B = Bz, dd = dzv) and true for the transaxial one (A = Axy, B = Bxy, dd = dyv). -1 for the planes that the
// detector does not see (p < firstPlane, or zmm <= 0). This is the single place where the formulas above are
// evaluated, cdrfWeights and the kernels that make their weights on their own use it. comp is the variance (voxels^2 of dd) that the
// rotation of the view blurs the image by, only used by the transaxial direction (the axial passes give zero); negative means the
// constant sigmaCoeff^2 that was used before the compensation depended on the view
DEVICE float cdrfSigma(const bool trans, const int p, const int firstPlane, const float dist, const float dxv, const float dd,
	const float A, const float B, const float C, const float sigmaCoeff, const float comp) {
	if (p < firstPlane)
		return -FLOAT_ONE;
	const float zmm = cdrfDepth(p, firstPlane, dist, dxv);
	if (!(zmm > FLOAT_ZERO))
		return -FLOAT_ONE;
	const float Rg = FMAX(A * zmm + B, FLOAT_ZERO);
#if CDRF_MODEL == 1
	return cdrfSquareParam(trans, Rg, C, dd, sigmaCoeff, comp, false);
#endif
	const float FWHM = SQRT(Rg * Rg + C * C) / dd;
	if (trans) {
		// The interpolation of the rotation blurs the image transaxially by comp voxels^2 (variance), which is compensated for here
		if (comp < FLOAT_ZERO)
			return SQRT(FMAX(FWHM * FWHM - FLOAT_ONE, 1e-16f)) * sigmaCoeff;
		// (the square of the sigma of the FWHM, not FWHM * FWHM: the product that the constant compensation above uses is then not shared,
		// and the compiler fuses it in the same way as it did before the compensation depended on the view)
		const float sg = FWHM * sigmaCoeff;
		return SQRT(FMAX(sg * sg - comp, 1e-16f * (sigmaCoeff * sigmaCoeff)));
	}
	return FWHM * sigmaCoeff;
}

// The same truncation as the global radius r (where the response is 1e-6 of its peak), but per plane
// The half base (voxels) of the triangle of the square model of plane p, see cdrfTap. Zero for the planes that the detector does not
// see, and for the models that do not use it (the Gaussian one)
DEVICE float cdrfHalfBase(const bool trans, const int p, const int firstPlane, const float dist, const float dxv, const float dd,
	const float A, const float B, const float C, const float sigmaCoeff, const float comp) {
#if CDRF_MODEL == 1
	if (p < firstPlane)
		return FLOAT_ZERO;
	const float zmm = cdrfDepth(p, firstPlane, dist, dxv);
	if (!(zmm > FLOAT_ZERO))
		return FLOAT_ZERO;
	return cdrfSquareParam(trans, FMAX(A * zmm + B, FLOAT_ZERO), C, dd, sigmaCoeff, comp, true);
#else
	return FLOAT_ZERO;
#endif
}

// The weight (before the normalization) of tap k of a plane with the sigma and the half base that cdrfSigma and cdrfHalfBase
// returned, denom = 2 sigma^2
DEVICE float cdrfTap(const int k, const float sigma, const float hb, const float denom) {
#if CDRF_MODEL == 1
	return cdrfSquareTap(k, hb, sigma);
#else
	return EXP(-CFLOAT(k * k) / denom);
#endif
}

// The radius of the plane: where the response has fallen to 1e-6 of its peak, in the same way for every model (the Gaussian model
// by 5.2565 sigma, the square model by the half base plus 5.2565 sigma of the Gaussian that it is convolved with), at most r
DEVICE int cdrfRadius(const float sigma, const float hb, const int r) {
#if CDRF_MODEL == 1
	return MIN(MAX((int)CEIL(hb + 5.2565f * sigma), 0), r);
#else
	return MIN(MAX((int)CEIL(5.2565f * sigma), 0), r);
#endif
}

// The normalized weights of the axial (wAx, GID1 = 0) and the transaxial (wTr, GID1 = 1) response of
// every plane, both of them (2 * r + 1) x nPlanes, and the radius of every plane, rAx and rTr (nPlanes
// each, zero for the planes that are not visible). One thread per plane and direction. The host makes them
// per view and volume when the region is large (the kernels that read them are cdrfAxialGW, cdrfTransFPDirect and
// cdrfTransBPDirect, which are also the "global weights" variants of the transaxial passes), and always for the radii
// that the kernels below cannot handle on their own (r > CDRF_RMAX, the *Direct kernels). The other kernels
// (cdrfAxial, cdrfTransFP, cdrfTransBP) make the weights of their planes in shared memory, with the same formulas
KERNEL void cdrfWeights(CLGLOBAL float* CLRESTRICT wAx, CLGLOBAL float* CLRESTRICT wTr,
	CLGLOBAL int* CLRESTRICT rAx, CLGLOBAL int* CLRESTRICT rTr,
	const int nPlanes, const int firstPlane, const int xOff, const int r,
	const float dist, const float dxv, const float dyv, const float dzv,
	const float Az, const float Bz, const float Axy, const float Bxy, const float C, const float sigmaCoeff, const float comp) {
	const int x = GID0;
	const int dir = GID1;
	if (x >= nPlanes || dir >= 2)
		return;
	CLGLOBAL float* w = (dir == 0) ? wAx : wTr;
	CLGLOBAL int* rp = (dir == 0) ? rAx : rTr;
	const int p = x + xOff;
	const float sigma = (dir == 0) ? cdrfSigma(false, p, firstPlane, dist, dxv, dzv, Az, Bz, C, sigmaCoeff, FLOAT_ZERO)
		: cdrfSigma(true, p, firstPlane, dist, dxv, dyv, Axy, Bxy, C, sigmaCoeff, comp);
	const float hb = (dir == 0) ? cdrfHalfBase(false, p, firstPlane, dist, dxv, dzv, Az, Bz, C, sigmaCoeff, FLOAT_ZERO)
		: cdrfHalfBase(true, p, firstPlane, dist, dxv, dyv, Axy, Bxy, C, sigmaCoeff, comp);
	const float denom = 2.f * sigma * sigma;
	const bool visible = sigma >= FLOAT_ZERO;
	if (!visible) {
		rp[x] = 0;
		for (int kk = -r; kk <= r; kk++)
			w[(LTYPE)(kk + r) * (LTYPE)nPlanes + (LTYPE)x] = FLOAT_ZERO;
		return;
	}
	const int radP = cdrfRadius(sigma, hb, r);
	rp[x] = radP;
	float wsum = FLOAT_ZERO;
	for (int kk = -radP; kk <= radP; kk++)
		wsum += cdrfTap(kk, sigma, hb, denom);
	for (int kk = -radP; kk <= radP; kk++)
		w[(LTYPE)(kk + r) * (LTYPE)nPlanes + (LTYPE)x] = cdrfTap(kk, sigma, hb, denom) / wsum;
}

// The tiled kernels below (cdrfAxial, cdrfTransFP, cdrfTransBP) use work-groups with a flat thread index
// tid = LID0 + LOCAL_SIZE * LID1 (LOCAL_SIZE x LOCAL_SIZE2 threads, the local size of the auxiliary
// kernels), which is split as tx = tid % CDRF_TX, the thread along the x run (the depth, the fastest array
// dimension, so that a warp / wavefront reads and writes contiguous memory), and l = tid / CDRF_TX, one of
// CDRF_NL lanes. CDRF_TX is 64 (a wavefront on AMD, two warps on NVIDIA) when the group size allows it.
// CDRF_TZPT is the number of axial outputs of a thread of cdrfAxial (a group makes CDRF_TZ = CDRF_TZPT *
// CDRF_NL rows), CDRF_YPT the number of transaxial outputs of a thread of cdrfTransBP. CDRF_RMAX is the
// largest radius that cdrfAxial can tile, see its shared memory array. The host passes the same values
// (they are defined with -D there) and launches cdrfAxialDirect instead when the radius is larger than CDRF_RMAX
#define CDRF_NT (LOCAL_SIZE * LOCAL_SIZE2)
#ifndef CDRF_TX
#if (CDRF_NT % 64) == 0
#define CDRF_TX 64
#elif (CDRF_NT % 32) == 0
#define CDRF_TX 32
#else
#define CDRF_TX 1
#endif
#endif
#define CDRF_NL (CDRF_NT / CDRF_TX)
#ifndef CDRF_TZPT
#define CDRF_TZPT 8
#endif
#ifndef CDRF_YPT
#define CDRF_YPT 4
#endif
#ifndef CDRF_RMAX
#define CDRF_RMAX 48
#endif
#define CDRF_TZ (CDRF_TZPT * CDRF_NL)
#define CDRF_ROWS (CDRF_TZ + 2 * CDRF_RMAX)

// The attenuation factor of element (x, y, z) of the region that cdrfAxial processes. The attenuation is read from its own array,
// of the row and slice strides attSx and attSy (the dimensions the array is allocated with) at the offset (attX0, attY0, attZ0)
// of that array, so a region that is a part of a larger attenuation array reads it in place. With attDx != 0 the array holds the
// accumulated attenuation (the sum along the depth, without the voxel size, the sum S_k including the plane k itself) and the factor
// is made of it as attConv says (TYPE6_ATTCONV, type6AttConv in precomp.h): 0 = exp(-attDx * S_k) (through the whole voxel),
// 1 = exp(-attDx * (S_(k-1) + S_k) / 2) (to the centre of the voxel), 2 = exp(-attDx * S_(k-1)) * g(t), t = attDx * (S_k - S_(k-1)),
// g(t) = (1 - exp(-t)) / t (the average over the voxel). S_(-1) = 0: the array starts at the detector side, so the plane before
// the first one of the array does not exist. Otherwise (attDx = 0) the array holds the factors
DEVICE float cdrfAttAt(const CLGLOBAL float* CLRESTRICT att, const int x, const int y, const int z,
	const int attX0, const int attY0, const int attZ0, const int attSx, const int attSy, const float attDx, const int attConv) {
	const LTYPE idx = (LTYPE)(x + attX0) + (LTYPE)(y + attY0) * (LTYPE)attSx + (LTYPE)(z + attZ0) * (LTYPE)attSx * (LTYPE)attSy;
	const float a = att[idx];
	if (attDx == FLOAT_ZERO)
		return a;
	if (attConv == 0)
		return EXP(-attDx * a);
	const float aPrev = (x + attX0) > 0 ? att[idx - (LTYPE)1] : FLOAT_ZERO;
	if (attConv == 1)
		return EXP(-attDx * (0.5f * (aPrev + a)));
	const float t = attDx * (a - aPrev);
	// The series of g for the small t, where 1 - exp(-t) loses the digits
	const float g = (t > -0.1f && t < 0.1f) ? (1.f - t * (0.5f - t * (1.f / 6.f - t * (1.f / 24.f - t * (1.f / 120.f))))) : ((1.f - EXP(-t)) / t);
	return EXP(-attDx * aPrev) * g;
}
#define CDRF_ATT(att, x, y, z) cdrfAttAt(att, x, y, z, attX0, attY0, attZ0, attSx, attSy, attDx, attConv)

// Axial pass. The input and the output (and the attenuation, when used, see cdrfAttAt) are logically N0 x N1 x N2 (x, y, z) and
// N0 x N1 x nZOut, the output holding the axial rows zStart ... zStart + nZOut - 1 of the input frame. The arrays may be allocated
// larger than that in the first two dimensions: inSx x inSy (the input) and outSx x outSy (the output) are the dimensions they are
// allocated with, i.e. their row and slice strides (N0 x N1 when they are not larger)
// attMode: 0 = no attenuation, 1 = the attenuation multiplies the input before the convolution (forward
// projection), 2 = the attenuation multiplies the output after the convolution (backprojection, the
// adjoint order)
// Tiled: a work-group is a run of CDRF_TX planes (x) of one y and a block of CDRF_TZ output rows, the group
// index is (x run, y, z block), so the launch is (ceil(N0 / CDRF_TX), N1, ceil(nZOut / CDRF_TZ)) groups. The
// group makes the weights of the axial response of its planes itself (the formulas of cdrfWeights, see
// cdrfSigma) into shared memory: the Gaussian is symmetric, so only the taps 0 ... r_p of a plane are kept,
// the exponentials are spread over the lanes of the group, every thread sums the same 2 r_p + 1 taps in
// the same order as cdrfWeights does and normalizes its share of the taps. The group then reads the rows of
// the block plus the radius on both sides, the largest radius of its planes, into shared memory (multiplied by
// the attenuation with attMode 1) once, and every thread, which has its own plane and so its own weights, makes
// CDRF_TZPT consecutive outputs from there. The rows outside of the input are zero, which is the zero
// padding. Needs r <= CDRF_RMAX. The shared memory is (CDRF_TZ + 3 * CDRF_RMAX + 1) * CDRF_TX floats (the
// rows, the weights) and the radii, 44.5 KB with 256 threads
KERNEL void cdrfAxial(CLGLOBAL float* CLRESTRICT out, const CLGLOBAL float* CLRESTRICT in, const CLGLOBAL float* CLRESTRICT att,
	const int N0, const int N1, const int N2, const int zStart, const int nZOut, const int attMode,
	const int inSx, const int inSy, const int outSx, const int outSy,
	const int attX0, const int attY0, const int attZ0, const int attSx, const int attSy, const float attDx, const int attConv,
	const int xOff, const int firstPlane, const int r,
	const float dist, const float dxv, const float dzv, const float Az, const float Bz, const float C, const float sigmaCoeff) {
	LOCAL float lIn[CDRF_ROWS * CDRF_TX];
	LOCAL float lW[(CDRF_RMAX + 1) * CDRF_TX];
	LOCAL int lR[CDRF_TX];
	const int tid = (int)LID0 + LOCAL_SIZE * (int)LID1;
	const int tx = tid % CDRF_TX;
	const int zl = tid / CDRF_TX;
	const int x = (int)GRID0 * CDRF_TX + tx;
	const int y = (int)GRID1;
	const int zo0 = (int)GRID2 * CDRF_TZ;
	const LTYPE base = (LTYPE)x + (LTYPE)y * (LTYPE)inSx;
	const LTYPE sN = (LTYPE)inSx * (LTYPE)inSy;
	const LTYPE oBase = (LTYPE)x + (LTYPE)y * (LTYPE)outSx;
	const LTYPE oSn = (LTYPE)outSx * (LTYPE)outSy;
	// The planes that the detector sees have a positive sigma, the rest have all-zero weights
	const float sigma = (x < N0) ? cdrfSigma(false, x + xOff, firstPlane, dist, dxv, dzv, Az, Bz, C, sigmaCoeff, FLOAT_ZERO) : -FLOAT_ONE;
	const float hb = (x < N0) ? cdrfHalfBase(false, x + xOff, firstPlane, dist, dxv, dzv, Az, Bz, C, sigmaCoeff, FLOAT_ZERO) : FLOAT_ZERO;
	const float denom = 2.f * sigma * sigma;
	const bool vis = x < N0 && sigma >= FLOAT_ZERO && denom > FLOAT_ZERO;
	const int rp = vis ? cdrfRadius(sigma, hb, r) : 0;
	if (zl == 0)
		lR[tx] = rp;
	if (vis) {
		for (int k = zl; k <= rp; k += CDRF_NL)
			lW[k * CDRF_TX + tx] = cdrfTap(k, sigma, hb, denom);
	}
	BARRIER
	int rg = 0;
	for (int i = 0; i < CDRF_TX; i++)
		rg = MAX(rg, lR[i]);
	// The normalization: the sum over the taps -rp ... rp in this order
	float wsum = FLOAT_ONE;
	if (vis) {
		wsum = FLOAT_ZERO;
		for (int k = rp; k >= 1; k--)
			wsum += lW[k * CDRF_TX + tx];
		for (int k = 0; k <= rp; k++)
			wsum += lW[k * CDRF_TX + tx];
	}
	// Shared rows 0 ... CDRF_TZ + 2 rg - 1 are the input rows z0 - rg ... z0 + CDRF_TZ + rg - 1
	const int z0 = zStart + zo0;
	const int nRows = CDRF_TZ + 2 * rg;
	if (vis) {
		for (int row = zl; row < nRows; row += CDRF_NL) {
			const int zz = z0 - rg + row;
			float v = FLOAT_ZERO;
			if (zz >= 0 && zz < N2) {
				const LTYPE idx = base + (LTYPE)zz * sN;
				v = in[idx];
				if (attMode == 1)
					v *= CDRF_ATT(att, x, y, zz);
			}
			lIn[row * CDRF_TX + tx] = v;
		}
	}
	// Everybody has summed the taps, so they can be normalized
	BARRIER
	if (vis) {
		for (int k = zl; k <= rp; k += CDRF_NL)
			lW[k * CDRF_TX + tx] = lW[k * CDRF_TX + tx] / wsum;
	}
	BARRIER
	if (x >= N0)
		return;
	float acc[CDRF_TZPT];
	for (int j = 0; j < CDRF_TZPT; j++)
		acc[j] = FLOAT_ZERO;
	if (vis) {
		const int t0 = zl * CDRF_TZPT;
		for (int kk = -rp; kk <= rp; kk++) {
			const float w = lW[(kk < 0 ? -kk : kk) * CDRF_TX + tx];
			const int rowB = t0 + kk + rg;
			for (int j = 0; j < CDRF_TZPT; j++)
				acc[j] += w * lIn[(rowB + j) * CDRF_TX + tx];
		}
	}
	for (int j = 0; j < CDRF_TZPT; j++) {
		const int zo = zo0 + zl * CDRF_TZPT + j;
		if (zo < nZOut) {
			float summa = acc[j];
			if (vis && attMode == 2)
				summa *= CDRF_ATT(att, x, y, zStart + zo);
			out[oBase + (LTYPE)zo * oSn] = summa;
		}
	}
}

// The same pass, tiled in the same way, but reading the weights and the radii that cdrfWeights made (wAx are
// the weights of cdrfWeights with nPlanes = N0, rAx the radii) instead of making them in shared memory: the
// variant that the host chooses for the large regions, where the weights are made once per view and volume and
// read by every group, which is cheaper than making them again in every group (see type6GlobalWeights in
// functions.hpp). The arguments and the launch are those of cdrfAxialDirect, the groups those of cdrfAxial.
// The shared memory is (CDRF_TZ + 2 * CDRF_RMAX) * CDRF_TX floats (the rows) and the radii, 32.3 KB with 256
// threads. Needs r <= CDRF_RMAX. The result is that of cdrfAxial (to the order of the sums)
KERNEL void cdrfAxialGW(CLGLOBAL float* CLRESTRICT out, const CLGLOBAL float* CLRESTRICT in, const CLGLOBAL float* CLRESTRICT att,
	const CLGLOBAL float* CLRESTRICT wAx, const CLGLOBAL int* CLRESTRICT rAx,
	const int N0, const int N1, const int N2, const int zStart, const int nZOut, const int attMode,
	const int inSx, const int inSy, const int outSx, const int outSy,
	const int attX0, const int attY0, const int attZ0, const int attSx, const int attSy, const float attDx, const int attConv,
	const int xOff, const int firstPlane, const int r) {
	LOCAL float lIn[CDRF_ROWS * CDRF_TX];
	LOCAL int lR[CDRF_TX];
	const int tid = (int)LID0 + LOCAL_SIZE * (int)LID1;
	const int tx = tid % CDRF_TX;
	const int zl = tid / CDRF_TX;
	const int x = (int)GRID0 * CDRF_TX + tx;
	const int y = (int)GRID1;
	const int zo0 = (int)GRID2 * CDRF_TZ;
	const LTYPE base = (LTYPE)x + (LTYPE)y * (LTYPE)inSx;
	const LTYPE sN = (LTYPE)inSx * (LTYPE)inSy;
	const LTYPE oBase = (LTYPE)x + (LTYPE)y * (LTYPE)outSx;
	const LTYPE oSn = (LTYPE)outSx * (LTYPE)outSy;
	// The weight of the center tap is positive for the planes that the detector sees and zero for the rest
	const bool vis = x < N0 && (x + xOff) >= firstPlane && wAx[(LTYPE)r * (LTYPE)N0 + (LTYPE)x] > FLOAT_ZERO;
	const int rp = vis ? rAx[x] : 0;
	if (zl == 0)
		lR[tx] = rp;
	BARRIER
	int rg = 0;
	for (int i = 0; i < CDRF_TX; i++)
		rg = MAX(rg, lR[i]);
	// Shared rows 0 ... CDRF_TZ + 2 rg - 1 are the input rows z0 - rg ... z0 + CDRF_TZ + rg - 1
	const int z0 = zStart + zo0;
	const int nRows = CDRF_TZ + 2 * rg;
	if (vis) {
		for (int row = zl; row < nRows; row += CDRF_NL) {
			const int zz = z0 - rg + row;
			float v = FLOAT_ZERO;
			if (zz >= 0 && zz < N2) {
				const LTYPE idx = base + (LTYPE)zz * sN;
				v = in[idx];
				if (attMode == 1)
					v *= CDRF_ATT(att, x, y, zz);
			}
			lIn[row * CDRF_TX + tx] = v;
		}
	}
	BARRIER
	if (x >= N0)
		return;
	float acc[CDRF_TZPT];
	for (int j = 0; j < CDRF_TZPT; j++)
		acc[j] = FLOAT_ZERO;
	if (vis) {
		const int t0 = zl * CDRF_TZPT;
		for (int kk = -rp; kk <= rp; kk++) {
			const float w = wAx[(LTYPE)(kk + r) * (LTYPE)N0 + (LTYPE)x];
			const int rowB = t0 + kk + rg;
			for (int j = 0; j < CDRF_TZPT; j++)
				acc[j] += w * lIn[(rowB + j) * CDRF_TX + tx];
		}
	}
	for (int j = 0; j < CDRF_TZPT; j++) {
		const int zo = zo0 + zl * CDRF_TZPT + j;
		if (zo < nZOut) {
			float summa = acc[j];
			if (vis && attMode == 2)
				summa *= CDRF_ATT(att, x, y, zStart + zo);
			out[oBase + (LTYPE)zo * oSn] = summa;
		}
	}
}

// The same pass without the tiling, one thread per output, for the radii that cdrfAxial cannot tile
// (r > CDRF_RMAX, where cdrfAxialGW cannot either). It reads the weights and the radii that cdrfWeights made (wAx are the weights of
// cdrfWeights with nPlanes = N0). Launched with (N0, N1, nZOut) threads rounded up to the local size
KERNEL void cdrfAxialDirect(CLGLOBAL float* CLRESTRICT out, const CLGLOBAL float* CLRESTRICT in, const CLGLOBAL float* CLRESTRICT att,
	const CLGLOBAL float* CLRESTRICT wAx, const CLGLOBAL int* CLRESTRICT rAx,
	const int N0, const int N1, const int N2, const int zStart, const int nZOut, const int attMode,
	const int inSx, const int inSy, const int outSx, const int outSy,
	const int attX0, const int attY0, const int attZ0, const int attSx, const int attSy, const float attDx, const int attConv,
	const int xOff, const int firstPlane, const int r) {
	const int x = GID0;
	const int y = GID1;
	const int zo = GID2;
	if (x >= N0 || y >= N1 || zo >= nZOut)
		return;
	const LTYPE ind = (LTYPE)x + (LTYPE)y * (LTYPE)outSx + (LTYPE)zo * (LTYPE)outSx * (LTYPE)outSy;
	float summa = FLOAT_ZERO;
	const int p = x + xOff;
	// The weight of the center tap is positive for the planes that the detector sees and zero for the rest
	if (p >= firstPlane && wAx[(LTYPE)r * (LTYPE)N0 + (LTYPE)x] > FLOAT_ZERO) {
		const int z = zStart + zo;
		const LTYPE base = (LTYPE)x + (LTYPE)y * (LTYPE)inSx;
		const LTYPE sN = (LTYPE)inSx * (LTYPE)inSy;
		const int rp = rAx[x];
		for (int kk = -rp; kk <= rp; kk++) {
			const int zz = z + kk;
			if (zz >= 0 && zz < N2) {
				const LTYPE idx = base + (LTYPE)zz * sN;
				float v = in[idx];
				if (attMode == 1)
					v *= CDRF_ATT(att, x, y, zz);
				summa += wAx[(LTYPE)(kk + r) * (LTYPE)N0 + (LTYPE)x] * v;
			}
		}
		if (attMode == 2)
			summa *= CDRF_ATT(att, x, y, z);
	}
	out[ind] = summa;
}

// Transaxial pass of the forward projection together with the sum over the depth, which is the projection.
// The input is N0 x N1 x N2 (x, y, z), allocated inSx x inSy x (at least N2) (the row and slice strides, N0 x N1 when it is not
// allocated larger), the output is nOut x N2 (trans, axial) of rows yStart ... yStart +
// nOut - 1 of the input, in the view'th slot of the array out (slots of nOut x N2 elements), multiplied by
// norm. Only the planes that the detector sees contribute (the others are zero). With accumulate the sum is
// added to what is there, otherwise it replaces it
// A work-group makes CDRF_NL * CDRF_YPT outputs (CDRF_YPT per lane, CDRF_NL apart) of one z, the group index
// is (ceil(nOut / (CDRF_NL * CDRF_YPT)), N2): the CDRF_TX threads of a lane step along the planes in chunks of
// CDRF_TX, each adds up its own planes and the partial sums are reduced in shared memory. The weights of the
// transaxial response of the planes of a chunk are made by the group itself into shared memory, in the
// way that cdrfAxial does, and are used for all of the CDRF_YPT outputs of a thread. Needs r <= CDRF_RMAX
KERNEL void cdrfTransFP(CLGLOBAL float* CLRESTRICT out, const CLGLOBAL float* CLRESTRICT in,
	const int N0, const int N1, const int N2, const int inSx, const int inSy, const int yStart, const int nOut,
	const int xOff, const int firstPlane, const int r, const float norm, const int view, const int accumulate,
	const float dist, const float dxv, const float dyv, const float Axy, const float Bxy, const float C, const float sigmaCoeff, const float comp) {
	LOCAL float lRed[CDRF_YPT * CDRF_NT];
	LOCAL float lW[(CDRF_RMAX + 1) * CDRF_TX];
	const int tid = (int)LID0 + LOCAL_SIZE * (int)LID1;
	const int tx = tid % CDRF_TX;
	const int yl = tid / CDRF_TX;
	const int yo0 = (int)GRID0 * (CDRF_NL * CDRF_YPT) + yl;
	const int z = (int)GRID1;
	const LTYPE zBase = (LTYPE)z * (LTYPE)inSx * (LTYPE)inSy;
	float part[CDRF_YPT];
	for (int q = 0; q < CDRF_YPT; q++)
		part[q] = FLOAT_ZERO;
	for (int x0 = 0; x0 < N0; x0 += CDRF_TX) {
		const int x = x0 + tx;
		const float sigma = (x < N0) ? cdrfSigma(true, x + xOff, firstPlane, dist, dxv, dyv, Axy, Bxy, C, sigmaCoeff, comp) : -FLOAT_ONE;
		const float hb = (x < N0) ? cdrfHalfBase(true, x + xOff, firstPlane, dist, dxv, dyv, Axy, Bxy, C, sigmaCoeff, comp) : FLOAT_ZERO;
		const float denom = 2.f * sigma * sigma;
		const bool vis = x < N0 && sigma >= FLOAT_ZERO && denom > FLOAT_ZERO;
		const int rp = vis ? cdrfRadius(sigma, hb, r) : 0;
		if (vis) {
			for (int k = yl; k <= rp; k += CDRF_NL)
				lW[k * CDRF_TX + tx] = cdrfTap(k, sigma, hb, denom);
		}
		BARRIER
		float wsum = FLOAT_ONE;
		if (vis) {
			wsum = FLOAT_ZERO;
			for (int k = rp; k >= 1; k--)
				wsum += lW[k * CDRF_TX + tx];
			for (int k = 0; k <= rp; k++)
				wsum += lW[k * CDRF_TX + tx];
		}
		BARRIER
		if (vis) {
			for (int k = yl; k <= rp; k += CDRF_NL)
				lW[k * CDRF_TX + tx] = lW[k * CDRF_TX + tx] / wsum;
		}
		BARRIER
		if (vis) {
			const LTYPE base = (LTYPE)x + zBase;
			float summa[CDRF_YPT];
			for (int q = 0; q < CDRF_YPT; q++)
				summa[q] = FLOAT_ZERO;
			for (int kk = -rp; kk <= rp; kk++) {
				const float w = lW[(kk < 0 ? -kk : kk) * CDRF_TX + tx];
				for (int q = 0; q < CDRF_YPT; q++) {
					const int yo = yo0 + q * CDRF_NL;
					const int yy = yStart + yo + kk;
					if (yo < nOut && yy >= 0 && yy < N1)
						summa[q] += w * in[base + (LTYPE)yy * (LTYPE)inSx];
				}
			}
			for (int q = 0; q < CDRF_YPT; q++)
				part[q] += summa[q];
		}
		// The next chunk overwrites the weights
		BARRIER
	}
	for (int q = 0; q < CDRF_YPT; q++)
		lRed[q * CDRF_NT + tid] = part[q];
	BARRIER
	for (int s = CDRF_TX / 2; s > 0; s >>= 1) {
		if (tx < s) {
			for (int q = 0; q < CDRF_YPT; q++)
				lRed[q * CDRF_NT + tid] += lRed[q * CDRF_NT + tid + s];
		}
		BARRIER
	}
	if (tx == 0) {
		for (int q = 0; q < CDRF_YPT; q++) {
			const int yo = yo0 + q * CDRF_NL;
			if (yo < nOut) {
				const LTYPE o = (LTYPE)view * (LTYPE)nOut * (LTYPE)N2 + (LTYPE)yo + (LTYPE)z * (LTYPE)nOut;
				const float summa = lRed[q * CDRF_NT + tid] * norm;
				if (accumulate != 0)
					out[o] += summa;
				else
					out[o] = summa;
			}
		}
	}
}

// The same pass with the weights and the radii that cdrfWeights made (wTr are the weights of cdrfWeights with
// nPlanes = N0). It is used for the radii that cdrfTransFP cannot handle (r > CDRF_RMAX) and as the global
// weights variant of it for the large regions (see type6GlobalWeights in functions.hpp), where reading the weights
// is cheaper than making them in every group. A work-group makes CDRF_NL outputs (one per lane) of one z, the
// group index is (ceil(nOut / CDRF_NL), N2)
KERNEL void cdrfTransFPDirect(CLGLOBAL float* CLRESTRICT out, const CLGLOBAL float* CLRESTRICT in, const CLGLOBAL float* CLRESTRICT wTr,
	const CLGLOBAL int* CLRESTRICT rTr,
	const int N0, const int N1, const int N2, const int inSx, const int inSy, const int yStart, const int nOut,
	const int xOff, const int firstPlane, const int r, const float norm, const int view, const int accumulate) {
	LOCAL float lRed[CDRF_NT];
	const int tid = (int)LID0 + LOCAL_SIZE * (int)LID1;
	const int tx = tid % CDRF_TX;
	const int yl = tid / CDRF_TX;
	const int yo = (int)GRID0 * CDRF_NL + yl;
	const int z = (int)GRID1;
	float part = FLOAT_ZERO;
	if (yo < nOut) {
		const int y = yStart + yo;
		const LTYPE zBase = (LTYPE)z * (LTYPE)inSx * (LTYPE)inSy;
		for (int x = tx; x < N0; x += CDRF_TX) {
			if ((x + xOff) >= firstPlane && wTr[(LTYPE)r * (LTYPE)N0 + (LTYPE)x] > FLOAT_ZERO) {
				const LTYPE base = (LTYPE)x + zBase;
				const int rp = rTr[x];
				float summa = FLOAT_ZERO;
				for (int kk = -rp; kk <= rp; kk++) {
					const int yy = y + kk;
					if (yy >= 0 && yy < N1)
						summa += wTr[(LTYPE)(kk + r) * (LTYPE)N0 + (LTYPE)x] * in[base + (LTYPE)yy * (LTYPE)inSx];
				}
				part += summa;
			}
		}
	}
	lRed[tid] = part;
	BARRIER
	for (int s = CDRF_TX / 2; s > 0; s >>= 1) {
		if (tx < s)
			lRed[tid] += lRed[tid + s];
		BARRIER
	}
	if (tx == 0 && yo < nOut) {
		const LTYPE o = (LTYPE)view * (LTYPE)nOut * (LTYPE)N2 + (LTYPE)yo + (LTYPE)z * (LTYPE)nOut;
		const float summa = lRed[tid] * norm;
		if (accumulate != 0)
			out[o] += summa;
		else
			out[o] = summa;
	}
}

// One element of the lattice of the rotated grid, P(i, j), from the view'th measurement (mBase is its offset) M
// (nRowsD x nColsD): the adjoint of resampleToDetector, see resampleFromDetector
DEVICE float cdrfFromDetectorAt(const CLGLOBAL float* CLRESTRICT M, const LTYPE mBase, const int i, const int j,
	const int nRowsD, const int nColsD, const float trA, const float trB, const float axA, const float axB) {
	const float fi = CFLOAT(i);
	const float fj = CFLOAT(j);
	const int r0 = (int)FMAX(FLOOR((fi - FLOAT_ONE - trB) / trA), FLOAT_ZERO);
	const int r1 = (int)FMIN(CEIL((fi + FLOAT_ONE - trB) / trA), CFLOAT(nRowsD - 1));
	const int c0 = (int)FMAX(FLOOR((fj - FLOAT_ONE - axB) / axA), FLOAT_ZERO);
	const int c1 = (int)FMIN(CEIL((fj + FLOAT_ONE - axB) / axA), CFLOAT(nColsD - 1));
	float summa = FLOAT_ZERO;
	for (int c = c0; c <= c1; c++) {
		const float wc = FMAX(FLOAT_ONE - FABS(axA * CFLOAT(c) + axB - fj), FLOAT_ZERO);
		if (wc <= FLOAT_ZERO)
			continue;
		float inner = FLOAT_ZERO;
		for (int rr = r0; rr <= r1; rr++) {
			const float wr = FMAX(FLOAT_ONE - FABS(trA * CFLOAT(rr) + trB - fi), FLOAT_ZERO);
			inner += wr * M[mBase + (LTYPE)rr + (LTYPE)c * (LTYPE)nRowsD];
		}
		summa += wc * inner;
	}
	return summa;
}

// Transaxial pass of the backprojection. The input is the projection P, pNy x pNz (trans, axial), which
// sits at (pY0, pZ0) inside the N1 x N2 region and is zero everywhere else, so the spreading of the
// projection over the depth and the multiplication by norm are done here as well. P is the view'th one of
// the projections that are stored one after another (view = 0 for a single one). The output is
// N0 x N1 x N2 (x, y, z), allocated outSx x outSy x (at least N2) (the row and slice strides, N0 x N1 when it is not allocated
// larger). With resample set P is not an array of projections but the measurements, nRowsD x
// nColsD each, and the projection is made here: its element (i, j) (pNy = N1, pNz = N2, pY0 = pZ0 = 0) is what
// resampleFromDetector makes from the view'th measurement (cdrfFromDetectorAt)
// A work-group is a run of CDRF_TX planes of one z and CDRF_NL * CDRF_YPT rows y, the group index is
// (ceil(N0 / CDRF_TX), ceil(N1 / (CDRF_NL * CDRF_YPT)), N2). The group makes the weights of its planes in
// shared memory as cdrfAxial does. Its rows of the projection, the rows of the group plus the largest radius
// on both sides, are read (or made) once into shared memory too (the projection does not depend on the
// plane, so all the threads of a warp read the same element of it). Every thread makes CDRF_YPT consecutive
// rows to read each of its weights once for all of them. Needs r <= CDRF_RMAX
KERNEL void cdrfTransBP(CLGLOBAL float* CLRESTRICT out, const CLGLOBAL float* CLRESTRICT P,
	const int N0, const int N1, const int N2, const int outSx, const int outSy, const int pNy, const int pNz, const int pY0, const int pZ0, const int view,
	const int xOff, const int firstPlane, const int r, const float norm,
	const float dist, const float dxv, const float dyv, const float Axy, const float Bxy, const float C, const float sigmaCoeff, const float comp,
	const int resample, const int nRowsD, const int nColsD, const float trA, const float trB, const float axA, const float axB) {
	LOCAL float lW[(CDRF_RMAX + 1) * CDRF_TX];
	LOCAL float lP[CDRF_NL * CDRF_YPT + 2 * CDRF_RMAX];
	LOCAL int lR[CDRF_TX];
	const int tid = (int)LID0 + LOCAL_SIZE * (int)LID1;
	const int tx = tid % CDRF_TX;
	const int yl = tid / CDRF_TX;
	const int x = (int)GRID0 * CDRF_TX + tx;
	const int z = (int)GRID2;
	const int yg0 = (int)GRID1 * (CDRF_NL * CDRF_YPT);
	const int yb = yg0 + yl * CDRF_YPT;
	const int pz = z - pZ0;
	const LTYPE zOut = (LTYPE)z * (LTYPE)outSx * (LTYPE)outSy;
	if (pz < 0 || pz >= pNz) {
		// The projection does not reach this z, which holds for the whole group
		if (x < N0) {
			for (int j = 0; j < CDRF_YPT; j++) {
				const int y = yb + j;
				if (y < N1)
					out[(LTYPE)x + (LTYPE)y * (LTYPE)outSx + zOut] = FLOAT_ZERO;
			}
		}
		return;
	}
	const float sigma = (x < N0) ? cdrfSigma(true, x + xOff, firstPlane, dist, dxv, dyv, Axy, Bxy, C, sigmaCoeff, comp) : -FLOAT_ONE;
	const float hb = (x < N0) ? cdrfHalfBase(true, x + xOff, firstPlane, dist, dxv, dyv, Axy, Bxy, C, sigmaCoeff, comp) : FLOAT_ZERO;
	const float denom = 2.f * sigma * sigma;
	const bool vis = x < N0 && sigma >= FLOAT_ZERO && denom > FLOAT_ZERO;
	const int rp = vis ? cdrfRadius(sigma, hb, r) : 0;
	if (yl == 0)
		lR[tx] = rp;
	if (vis) {
		for (int k = yl; k <= rp; k += CDRF_NL)
			lW[k * CDRF_TX + tx] = cdrfTap(k, sigma, hb, denom);
	}
	BARRIER
	int rg = 0;
	for (int i = 0; i < CDRF_TX; i++)
		rg = MAX(rg, lR[i]);
	float wsum = FLOAT_ONE;
	if (vis) {
		wsum = FLOAT_ZERO;
		for (int k = rp; k >= 1; k--)
			wsum += lW[k * CDRF_TX + tx];
		for (int k = 0; k <= rp; k++)
			wsum += lW[k * CDRF_TX + tx];
	}
	// Shared elements 0 ... CDRF_NL * CDRF_YPT + 2 rg - 1 are the rows yg0 - rg ... of the projection, zero
	// where the projection is not
	const int nW = CDRF_NL * CDRF_YPT + 2 * rg;
	const LTYPE pBase = (LTYPE)view * (LTYPE)pNy * (LTYPE)pNz + (LTYPE)pz * (LTYPE)pNy;
	const LTYPE mBase = (LTYPE)view * (LTYPE)nRowsD * (LTYPE)nColsD;
	for (int i = tid; i < nW; i += CDRF_NT) {
		const int yy = yg0 - rg + i;
		float v = FLOAT_ZERO;
		if (yy >= 0 && yy < N1) {
			if (resample != 0) {
				v = cdrfFromDetectorAt(P, mBase, yy, z, nRowsD, nColsD, trA, trB, axA, axB);
			}
			else {
				const int py = yy - pY0;
				if (py >= 0 && py < pNy)
					v = P[pBase + (LTYPE)py];
			}
		}
		lP[i] = v;
	}
	BARRIER
	if (vis) {
		for (int k = yl; k <= rp; k += CDRF_NL)
			lW[k * CDRF_TX + tx] = lW[k * CDRF_TX + tx] / wsum;
	}
	BARRIER
	if (x >= N0 || yb >= N1)
		return;
	float acc[CDRF_YPT];
	for (int j = 0; j < CDRF_YPT; j++)
		acc[j] = FLOAT_ZERO;
	if (vis) {
		// Row yb + j + kk is the shared element yb + j + kk - (yg0 - rg), which the zeros of the window make
		// the same as leaving out the rows that are not in the projection
		const int e0 = yl * CDRF_YPT + rg;
		for (int kk = -rp; kk <= rp; kk++) {
			const float w = lW[(kk < 0 ? -kk : kk) * CDRF_TX + tx];
			for (int j = 0; j < CDRF_YPT; j++)
				acc[j] += w * lP[e0 + j + kk];
		}
	}
	for (int j = 0; j < CDRF_YPT; j++) {
		const int y = yb + j;
		if (y < N1)
			out[(LTYPE)x + (LTYPE)y * (LTYPE)outSx + zOut] = vis ? acc[j] * norm : FLOAT_ZERO;
	}
}

// The same pass with the weights and the radii that cdrfWeights made (wTr are the weights of cdrfWeights with
// nPlanes = N0) and the projection P, which is an array (no resampling here, resampleFromDetector makes it). It is
// used for the radii that cdrfTransBP cannot handle (r > CDRF_RMAX) and as the global weights variant of it for the
// large regions (see type6GlobalWeights in functions.hpp). A work-group is a run of CDRF_TX planes of one z
// and CDRF_NL * CDRF_YPT rows y, the group index is (ceil(N0 / CDRF_TX), ceil(N1 / (CDRF_NL * CDRF_YPT)), N2)
KERNEL void cdrfTransBPDirect(CLGLOBAL float* CLRESTRICT out, const CLGLOBAL float* CLRESTRICT P, const CLGLOBAL float* CLRESTRICT wTr,
	const CLGLOBAL int* CLRESTRICT rTr,
	const int N0, const int N1, const int N2, const int outSx, const int outSy, const int pNy, const int pNz, const int pY0, const int pZ0, const int view,
	const int xOff, const int firstPlane, const int r, const float norm) {
	const int tid = (int)LID0 + LOCAL_SIZE * (int)LID1;
	const int tx = tid % CDRF_TX;
	const int yl = tid / CDRF_TX;
	const int x = (int)GRID0 * CDRF_TX + tx;
	const int z = (int)GRID2;
	const int yb = ((int)GRID1 * CDRF_NL + yl) * CDRF_YPT;
	if (x >= N0 || yb >= N1)
		return;
	float acc[CDRF_YPT];
	for (int j = 0; j < CDRF_YPT; j++)
		acc[j] = FLOAT_ZERO;
	const int pz = z - pZ0;
	const bool act = (x + xOff) >= firstPlane && pz >= 0 && pz < pNz && wTr[(LTYPE)r * (LTYPE)N0 + (LTYPE)x] > FLOAT_ZERO;
	if (act) {
		const LTYPE base = (LTYPE)view * (LTYPE)pNy * (LTYPE)pNz + (LTYPE)pz * (LTYPE)pNy;
		const int rp = rTr[x];
		for (int kk = -rp; kk <= rp; kk++) {
			const float w = wTr[(LTYPE)(kk + r) * (LTYPE)N0 + (LTYPE)x];
			for (int j = 0; j < CDRF_YPT; j++) {
				const int yy = yb + j + kk;
				const int py = yy - pY0;
				if (yy >= 0 && yy < N1 && py >= 0 && py < pNy)
					acc[j] += w * P[base + (LTYPE)py];
			}
		}
	}
	for (int j = 0; j < CDRF_YPT; j++) {
		const int y = yb + j;
		if (y < N1)
			out[(LTYPE)x + (LTYPE)y * (LTYPE)outSx + (LTYPE)z * (LTYPE)outSx * (LTYPE)outSy] = act ? acc[j] * norm : FLOAT_ZERO;
	}
}

// Resamples the projection P, nTr x nAx (trans, axial) on the lattice of the rotated grid, onto the
// detector and adds it to the view'th projection of fProj (nRowsD x nColsD each):
//   fProj(r, c, view) += sum_i sum_j tent(trA * r + trB - i) * tent(axA * c + axB - j) * P(i, j)
// with tent(t) = max(0, 1 - |t|), so only the 2 x 2 neighbours of the sampling position have a weight.
// Samples that are outside of P have no weight, i.e. read as zero. One thread per detector bin
KERNEL void resampleToDetector(CLGLOBAL float* CLRESTRICT fProj, const CLGLOBAL float* CLRESTRICT P,
	const int nTr, const int nAx, const int nRowsD, const int nColsD,
	const float trA, const float trB, const float axA, const float axB, const int view) {
	const int r = GID0;
	const int c = GID1;
	if (r >= nRowsD || c >= nColsD)
		return;
	const float tt = trA * CFLOAT(r) + trB;
	const float ta = axA * CFLOAT(c) + axB;
	const int i0 = (int)FLOOR(tt);
	const int j0 = (int)FLOOR(ta);
	float summa = FLOAT_ZERO;
	for (int jj = 0; jj < 2; jj++) {
		const int j = j0 + jj;
		if (j < 0 || j >= nAx)
			continue;
		const float wj = FMAX(FLOAT_ONE - FABS(ta - CFLOAT(j)), FLOAT_ZERO);
		for (int ii = 0; ii < 2; ii++) {
			const int i = i0 + ii;
			if (i < 0 || i >= nTr)
				continue;
			const float wi = FMAX(FLOAT_ONE - FABS(tt - CFLOAT(i)), FLOAT_ZERO);
			summa += wi * wj * P[(LTYPE)i + (LTYPE)j * (LTYPE)nTr];
		}
	}
	fProj[(LTYPE)view * (LTYPE)nRowsD * (LTYPE)nColsD + (LTYPE)r + (LTYPE)c * (LTYPE)nRowsD] += summa;
}

// The exact adjoint of resampleToDetector, from the view'th measurement M (nRowsD x nColsD) of an array
// of them to the lattice of the rotated grid P, nTr x nAx (trans, axial):
//   P(i, j) = sum_r sum_c tent(trA * r + trB - i) * tent(axA * c + axB - j) * M(r, c)
// The detector bins with a weight are those with |trA * r + trB - i| < 1, i.e. r between
// (i - 1 - trB) / trA and (i + 1 - trB) / trA (trA > 0), and likewise for c, so only those are visited.
// The ends of the ranges are rounded outwards, a bin that is there but has no weight adds exactly zero.
// One thread per element of P, every one of them is written
KERNEL void resampleFromDetector(CLGLOBAL float* CLRESTRICT P, const CLGLOBAL float* CLRESTRICT M,
	const int nTr, const int nAx, const int nRowsD, const int nColsD,
	const float trA, const float trB, const float axA, const float axB, const int view) {
	const int i = GID0;
	const int j = GID1;
	if (i >= nTr || j >= nAx)
		return;
	const LTYPE mBase = (LTYPE)view * (LTYPE)nRowsD * (LTYPE)nColsD;
	P[(LTYPE)i + (LTYPE)j * (LTYPE)nTr] = cdrfFromDetectorAt(M, mBase, i, j, nRowsD, nColsD, trA, trB, axA, axB);
}

// ---------------------------------------------------------------------------------------------------------------
// Batched variants of the fused type-6 kernels (TYPE6_BATCH, see functions.hpp). One launch handles nb
// consecutive views of one volume: the third dimension of the launch is the batch index (times the number of z
// blocks / z slices of a view), the per-view arrays are stacked with a batch stride (a view's slice b is at b * the
// stride), and everything that differs from view to view is read from the parameter buffer bp, nb blocks of
// T6P_STRIDE floats, one per view, starting at the block bpOff (the buffer holds the blocks of all the launches of a projection). Every integer in a block is stored exactly in a float (they are all below 2^24). The
// indices below are mirrored by type6BatchParam in precomp.h. The arithmetic of every output element is that of the
// unbatched kernel (cdrfAxial, cdrfTransFP, cdrfTransBP, resampleToDetector, rotate); only the tiles that the view
// does not need (its extents are smaller than the largest of the batch) return at once, which they can, since the
// tests only use the group index. The batched kernels make the weights of the response in shared memory like the
// unbatched ones do, so only the views that do not use the global weights (type6GlobalWeights) are batched
// The rotation parameters (source and output windows of the rotation of a view)
// The coordinates of the rotation kernels (rotate, rotateB, rotateBackB) are made with explicit fused multiply-adds, in the same order in all
// of them: left to the compiler, the choice of the product that is fused in a*b + c*d differed between the kernels, which made the batched and
// the unbatched rotations (and so everything that is interpolated at the rotated positions) differ in the last bit
#if defined(CUDA) || defined(HIP)
#define T6FMA(a, b, c) fmaf((a), (b), (c))
#elif defined(OPENCL)
#define T6FMA(a, b, c) fma((a), (b), (c))
#endif
// The number of slices (z) of the volume that a thread of the rotation back of the backprojection handles (rotate with accumulate,
// rotateBackB), ROT_ZPT, a build option of the host that also makes the launch ranges of those kernels ROT_ZPT times shorter in z.
// The weights of the exact transpose of the rotation do not depend on z, so the thread makes them once for its slices
#ifndef ROT_ZPT
#define ROT_ZPT 1
#endif
#define T6P_STRIDE 53
#define T6P_SRCNX 0
#define T6P_SRCNY 1
#define T6P_SRCNZ 2
#define T6P_OUTNX 3
#define T6P_OUTNYF 4
#define T6P_COSA 5
#define T6P_SINA 6
#define T6P_OFFX 7
#define T6P_OFFY 8
#define T6P_OUTY 9
#define T6P_OUTZ 10
#define T6P_OUTNY 11
#define T6P_OUTNZ 12
#define T6P_SRCY 13
#define T6P_SRCZ 14
#define T6P_SRCNYB 15
#define T6P_SRCNZB 16
#define T6P_OUTX 17
#define T6P_OUTNXW 18
#define T6P_SRCX 19
#define T6P_SRCNXB 20
// The region of the axial pass (and of the transaxial ones) and the response of the view
#define T6P_N0 21
#define T6P_N1 22
#define T6P_N2 23
#define T6P_ZSTART 24
#define T6P_NZOUT 25
#define T6P_XOFF 26
#define T6P_FIRSTPLANE 27
#define T6P_RADIUS 28
// The attenuation of the view: the offset of its region in its slice of the attenuation array (the offset along z includes
// the slice), the voxel size, and the distance from the first plane to the collimator surface
#define T6P_ATTX0 29
#define T6P_ATTY0 30
#define T6P_ATTZ0 31
#define T6P_ATTDX 32
#define T6P_DIST 33
// The transaxial passes: the rows that the forward projection keeps, the slot of the output that the view goes to, the
// view of the projections array and the projection that the backprojection spreads (see cdrfTransBP), the map onto the
// detector and the size of the projection that resampleToDetectorB resamples
#define T6P_YSTART 34
#define T6P_NOUT 35
#define T6P_OUTSLOT 36
#define T6P_VIEW 37
#define T6P_PNY 38
#define T6P_PNZ 39
#define T6P_PY0 40
#define T6P_PZ0 41
#define T6P_TRA 42
#define T6P_TRB 43
#define T6P_AXA 44
#define T6P_AXB 45
#define T6P_NTR 46
#define T6P_NAX 47
// The window of the stacked arrays of a view (region coordinates, see type6WindowsType6 in functions.hpp): the arrays hold the rows
// T6P_Y0 ... T6P_Y0 + T6P_NY - 1 (index y - T6P_Y0) of the region along y, and the input of the axial pass (the rotated region in the forward
// projection, the output of cdrfTransBPB in the backprojection) the axial rows T6P_IZ0 ... T6P_IZ0 + T6P_INZ - 1 (index z - T6P_IZ0); the
// elements outside of the windows are not read (they are zero). The default (the whole region) is the behaviour without the windows
#define T6P_Y0 48
#define T6P_NY 49
#define T6P_IZ0 50
#define T6P_INZ 51
// The variance (voxels^2 of the transaxial voxel size) that the interpolation of the rotation of the view adds transaxially, which the
// transaxial response is narrowed by (rotVarType6 in functions.hpp), or negative for the constant sigmaCoeff^2
#define T6P_COMP 52
#define T6PI(p, i) ((int)(p)[i])

// Batched rotation of the forward projection, see rotate. The output is nb slices, slice b (the one of view b, outSx x outSy x
// outSz elements, its window being outNxW x outNy x outNz of them) at b * outSx * outSy * outSz. The source is the same for
// every view: the texture or buffer of the volume (or of the attenuation source). The launch is (max outNxW, max outNy, nzB * nb)
// threads, nzB the largest outNz of the views: thread z is the slice z % nzB of the view z / nzB. Always a plain rotation
// (the elements that the rotation does not reach are zero), never accumulates
KERNEL void rotateB(
	CLGLOBAL float* CLRESTRICT rotim BUF0,
	IMTYPE im TEX1,
	const CLGLOBAL float* CLRESTRICT bp, const int bpOff,
	const int outSx, const int outSy, const int outSz, const int nzB) {
	const int ix = (int)GID0;
	const int iy = (int)GID1;
	const int gz = (int)GID2;
	const int b = gz / nzB;
	const int iz = gz - b * nzB;
	const CLGLOBAL float* p = bp + ((LTYPE)b + (LTYPE)bpOff) * (LTYPE)T6P_STRIDE;
	const int outNxW = T6PI(p, T6P_OUTNXW);
	const int outNy = T6PI(p, T6P_OUTNY);
	const int outNz = T6PI(p, T6P_OUTNZ);
	if (ix >= outNxW || iy >= outNy || iz >= outNz)
		return;
	const int Nx = T6PI(p, T6P_SRCNX);
	const int Ny = T6PI(p, T6P_SRCNY);
	const int outNx = T6PI(p, T6P_OUTNX);
	const int outNyF = T6PI(p, T6P_OUTNYF);
	const float cosa = p[T6P_COSA];
	const float sina = p[T6P_SINA];
	const float offX = p[T6P_OFFX];
	const float offY = p[T6P_OFFY];
	const int outX = T6PI(p, T6P_OUTX);
	const int outY = T6PI(p, T6P_OUTY);
	const int outZ = T6PI(p, T6P_OUTZ);
	const int srcX = T6PI(p, T6P_SRCX);
	const int srcY = T6PI(p, T6P_SRCY);
	const int srcZ = T6PI(p, T6P_SRCZ);
	const int srcNxB = T6PI(p, T6P_SRCNXB);
	const int srcNy = T6PI(p, T6P_SRCNYB);
	const int srcNz = T6PI(p, T6P_SRCNZB);
	const LTYPE n = (LTYPE)ix + (LTYPE)iy * (LTYPE)outSx + (LTYPE)iz * ((LTYPE)outSx * (LTYPE)outSy)
		+ (LTYPE)b * ((LTYPE)outSx * (LTYPE)outSy * (LTYPE)outSz);
	// Position in the full volume
	const int x = outX + ix;
	const int y = outY + iy;
	const int z = outZ + iz;

	// The rotation centre of a grid of n elements is its exact centre, n / 2 (real division): element x is at x + 0.5 - n / 2 from it, i.e. (2 x + 1 - n) / 2,
	// and the centre of the source in the source index is Nx / 2 - 0.5. For the odd n it is the middle of the middle element (not the boundary n / 2 (integer
	// division) of before, which was half an element off the axis that the host assumes, see offsetsFPType6, depthWindowType6, type6RotCentre in precomp.h)
	const float xA = (float)(2 * x + 1 - outNx) * FLOAT_HALF;
	const float yA = (float)(2 * y + 1 - outNyF) * FLOAT_HALF;

	const float src_x = (T6FMA(xA, cosa, -(yA * sina)) + offX + (float)Nx * FLOAT_HALF) - FLOAT_HALF;
	const float src_y = (T6FMA(yA, cosa, xA * sina) + offY + (float)Ny * FLOAT_HALF) - FLOAT_HALF;

	// The source is bilinearly interpolated with zero padding, exactly as in rotate
	const int xLo = srcX;
	const int xHi = srcX + srcNxB - 1;
	const int yLo = srcY;
	const int yHi = srcY + srcNy - 1;

	if (src_x > CFLOAT(xLo) - FLOAT_ONE && src_x < CFLOAT(xHi) + FLOAT_ONE && src_y > CFLOAT(yLo) - FLOAT_ONE && src_y < CFLOAT(yHi) + FLOAT_ONE
		&& z >= srcZ && z < srcZ + srcNz) {
		float val = FLOAT_ZERO;
#ifdef USEIMAGES
		if (src_x >= CFLOAT(xLo) && src_x <= CFLOAT(xHi) && src_y >= CFLOAT(yLo) && src_y <= CFLOAT(yHi)) {
#if defined(CUDA) || defined(HIP)
			val = tex3D<float>(im, src_x + FLOAT_HALF, src_y + FLOAT_HALF, CFLOAT(z) + FLOAT_HALF);
#elif defined(OPENCL)
			val = read_imagef(im, samplerRotate, (float4)(src_x + FLOAT_HALF, src_y + FLOAT_HALF, CFLOAT(z) + FLOAT_HALF, FLOAT_ZERO)).w;
#endif
		}
		else {
			const int x0 = (int)FLOOR(src_x);
			const int y0 = (int)FLOOR(src_y);
			const float fx = src_x - CFLOAT(x0);
			const float fy = src_y - CFLOAT(y0);
			const float zc = CFLOAT(z) + FLOAT_HALF;
			for (int jj = 0; jj < 2; jj++) {
				const int yy = y0 + jj;
				if (yy < yLo || yy > yHi)
					continue;
				const float wy = (jj == 0) ? (FLOAT_ONE - fy) : fy;
				for (int ii = 0; ii < 2; ii++) {
					const int xx = x0 + ii;
					if (xx < xLo || xx > xHi)
						continue;
					const float w = ((ii == 0) ? (FLOAT_ONE - fx) : fx) * wy;
#if defined(CUDA) || defined(HIP)
					val += w * tex3D<float>(im, CFLOAT(xx) + FLOAT_HALF, CFLOAT(yy) + FLOAT_HALF, zc);
#elif defined(OPENCL)
					val += w * read_imagef(im, samplerRotate, (float4)(CFLOAT(xx) + FLOAT_HALF, CFLOAT(yy) + FLOAT_HALF, zc, FLOAT_ZERO)).w;
#endif
				}
			}
		}
#else
		const int x0 = (int)FLOOR(src_x);
		const int y0 = (int)FLOOR(src_y);
		const float fx = src_x - CFLOAT(x0);
		const float fy = src_y - CFLOAT(y0);
		const LTYPE plane = (LTYPE)(z - srcZ) * (LTYPE)srcNxB * (LTYPE)srcNy;
		for (int jj = 0; jj < 2; jj++) {
			const int yy = y0 + jj;
			if (yy < yLo || yy > yHi)
				continue;
			const float wy = (jj == 0) ? (FLOAT_ONE - fy) : fy;
			for (int ii = 0; ii < 2; ii++) {
				const int xx = x0 + ii;
				if (xx < xLo || xx > xHi)
					continue;
				const float w = ((ii == 0) ? (FLOAT_ONE - fx) : fx) * wy;
				val += w * im[plane + (LTYPE)((xx - srcX) + (yy - srcY) * srcNxB)];
			}
		}
#endif
		rotim[n] = val;
	}
	else {
		rotim[n] = 0.0f;
	}
}

// Batched rotation back of the backprojection into the image out (outNxD x outNyD x outNzD, the whole volume), from the nb
// regions of src, slice b (srcSx x srcSy x any elements, the window of the view in it, see rotate with the source and output
// swapped) at b * srcBS. Every output voxel is handled by one thread, which adds the views in their order to the value of the
// voxel and writes it once, so there is no conflict between the views and the result is that of rotating the views one after the
// other with accumulate (the elements that no view reaches are not written). The source is always read from a buffer
KERNEL void rotateBackB(
	CLGLOBAL float* CLRESTRICT out,
	const CLGLOBAL float* CLRESTRICT src,
	const CLGLOBAL float* CLRESTRICT bp, const int bpOff,
	const int nb, const int outNxD, const int outNyD, const int outNzD,
	const int srcSx, const int srcSy, const int srcBS) {
#ifdef TYPE6_BP_FAST
	const int ix = (int)GID0;
	const int iy = (int)GID1;
	const int iz = (int)GID2;
	if (ix >= outNxD || iy >= outNyD || iz >= outNzD)
		return;
	const LTYPE n = (LTYPE)ix + (LTYPE)iy * (LTYPE)outNxD + (LTYPE)iz * ((LTYPE)outNxD * (LTYPE)outNyD);
	float acc = out[n];
	bool any = false;
	for (int b = 0; b < nb; b++) {
		const CLGLOBAL float* p = bp + ((LTYPE)b + (LTYPE)bpOff) * (LTYPE)T6P_STRIDE;
		const int Nx = T6PI(p, T6P_SRCNX);
		const int Ny = T6PI(p, T6P_SRCNY);
		const int outNx = T6PI(p, T6P_OUTNX);
		const int outNyF = T6PI(p, T6P_OUTNYF);
		const float cosa = p[T6P_COSA];
		const float sina = p[T6P_SINA];
		const float offX = p[T6P_OFFX];
		const float offY = p[T6P_OFFY];
		const int x = T6PI(p, T6P_OUTX) + ix;
		const int y = T6PI(p, T6P_OUTY) + iy;
		const int z = T6PI(p, T6P_OUTZ) + iz;
		const float xA = (float)(2 * x + 1 - outNx) * FLOAT_HALF;
		const float yA = (float)(2 * y + 1 - outNyF) * FLOAT_HALF;
		const float src_x = (T6FMA(xA, cosa, -(yA * sina)) + offX + (float)Nx * FLOAT_HALF) - FLOAT_HALF;
		const float src_y = (T6FMA(yA, cosa, xA * sina) + offY + (float)Ny * FLOAT_HALF) - FLOAT_HALF;
		const int srcX = T6PI(p, T6P_SRCX);
		const int srcY = T6PI(p, T6P_SRCY);
		const int srcZ = T6PI(p, T6P_SRCZ);
		const int xLo = srcX;
		const int xHi = srcX + T6PI(p, T6P_SRCNXB) - 1;
		const int yLo = srcY;
		const int yHi = srcY + T6PI(p, T6P_SRCNYB) - 1;
		if (src_x > CFLOAT(xLo) - FLOAT_ONE && src_x < CFLOAT(xHi) + FLOAT_ONE && src_y > CFLOAT(yLo) - FLOAT_ONE && src_y < CFLOAT(yHi) + FLOAT_ONE
			&& z >= srcZ && z < srcZ + T6PI(p, T6P_SRCNZB)) {
			float val = FLOAT_ZERO;
			const int x0 = (int)FLOOR(src_x);
			const int y0 = (int)FLOOR(src_y);
			const float fx = src_x - CFLOAT(x0);
			const float fy = src_y - CFLOAT(y0);
			const LTYPE plane = (LTYPE)b * (LTYPE)srcBS + (LTYPE)(z - srcZ) * (LTYPE)srcSx * (LTYPE)srcSy;
			for (int jj = 0; jj < 2; jj++) {
				const int yy = y0 + jj;
				if (yy < yLo || yy > yHi)
					continue;
				const float wy = (jj == 0) ? (FLOAT_ONE - fy) : fy;
				for (int ii = 0; ii < 2; ii++) {
					const int xx = x0 + ii;
					if (xx < xLo || xx > xHi)
						continue;
					const float w = ((ii == 0) ? (FLOAT_ONE - fx) : fx) * wy;
					val += w * src[plane + (LTYPE)((xx - srcX) + (yy - srcY) * srcSx)];
				}
			}
			acc += val;
			any = true;
		}
	}
	if (any)
		out[n] = acc;
}
#else
	const int ix = (int)GID0;
	const int iy = (int)GID1;
	// Every thread handles ROT_ZPT consecutive slices of its (x, y) column, so that the weights of the elements of the views, which do not
	// depend on z, are made once for them
	const int iz0 = (int)GID2 * ROT_ZPT;
	if (ix >= outNxD || iy >= outNyD || iz0 >= outNzD)
		return;
	const LTYPE nxy = (LTYPE)outNxD * (LTYPE)outNyD;
	const LTYPE n0 = (LTYPE)ix + (LTYPE)iy * (LTYPE)outNxD + (LTYPE)iz0 * nxy;
	// The slices are written only if some view reaches them with an element (as in the gather of the fast rotation), so the part of the
	// volume that no element of the regions reaches is not written
	float acc[ROT_ZPT];
	for (int j = 0; j < ROT_ZPT; j++)
		acc[j] = (iz0 + j < outNzD) ? out[n0 + (LTYPE)j * nxy] : FLOAT_ZERO;
	int anyMask = 0;
#pragma unroll 4
	for (int b = 0; b < nb; b++) {
		const CLGLOBAL float* p = bp + ((LTYPE)b + (LTYPE)bpOff) * (LTYPE)T6P_STRIDE;
		// The slices of the column that the region of the view reaches: slice iz of the volume is slice iz - zlo of the region
		const int zlo = T6PI(p, T6P_SRCZ) - T6PI(p, T6P_OUTZ);
		const int zhi = zlo + T6PI(p, T6P_SRCNZB);
		int act = 0;
		for (int j = 0; j < ROT_ZPT; j++) {
			const int izj = iz0 + j;
			if (izj < outNzD && izj >= zlo && izj < zhi)
				act |= 1 << j;
		}
		if (act == 0)
			continue;
		const int Nx = T6PI(p, T6P_SRCNX);
		const int Ny = T6PI(p, T6P_SRCNY);
		const int outNx = T6PI(p, T6P_OUTNX);
		const int outNyF = T6PI(p, T6P_OUTNYF);
		const float cosa = p[T6P_COSA];
		const float sina = p[T6P_SINA];
		const float offX = p[T6P_OFFX];
		const float offY = p[T6P_OFFY];
		const int x = T6PI(p, T6P_OUTX) + ix;
		const int y = T6PI(p, T6P_OUTY) + iy;
		const float xA = (float)(2 * x + 1 - outNx) * FLOAT_HALF;
		const float yA = (float)(2 * y + 1 - outNyF) * FLOAT_HALF;
		const float src_x = (T6FMA(xA, cosa, -(yA * sina)) + offX + (float)Nx * FLOAT_HALF) - FLOAT_HALF;
		const float src_y = (T6FMA(yA, cosa, xA * sina) + offY + (float)Ny * FLOAT_HALF) - FLOAT_HALF;
		const int srcX = T6PI(p, T6P_SRCX);
		const int srcY = T6PI(p, T6P_SRCY);
		const int xLo = srcX;
		const int xHi = srcX + T6PI(p, T6P_SRCNXB) - 1;
		const int yLo = srcY;
		const int yHi = srcY + T6PI(p, T6P_SRCNYB) - 1;
		// The exact transpose of the bilinear rotation of the forward projection: that rotation samples the image at the position
		// P(q) = R^T (q - off) (q the element of the region, both centered on their own grids) with the weight tri(P(q) - voxel), so the
		// voxel receives sum_q tri(P(q) - voxel) * region(q) over the elements q whose sample position is within one voxel of it. Those
		// are the q strictly within |cosa| + |sina| of src_x, src_y (the position that the inverse rotation would sample), the others have
		// no weight. The inverse rotation itself (the interpolation at src_x, src_y) is not the transpose of the interpolation of the
		// forward projection, which made the projector pair mismatched for everything that is not smooth
		const float hh = FABS(cosa) + FABS(sina);
		const int qx0 = MAX((int)FLOOR(src_x - hh) + 1, xLo);
		const int qx1 = MIN((int)CEIL(src_x + hh) - 1, xHi);
		const int qy0 = MAX((int)FLOOR(src_y - hh) + 1, yLo);
		const int qy1 = MIN((int)CEIL(src_y + hh) - 1, yHi);
		const LTYPE sxy = (LTYPE)srcSx * (LTYPE)srcSy;
		const LTYPE plane = (LTYPE)b * (LTYPE)srcBS + (LTYPE)(iz0 - zlo) * sxy;
		float val[ROT_ZPT];
		for (int j = 0; j < ROT_ZPT; j++)
			val[j] = FLOAT_ZERO;
		bool reached = false;
		for (int qy = qy0; qy <= qy1; qy++) {
			const float uy = (float)(2 * qy + 1 - Ny) * FLOAT_HALF - offY;
			for (int qx = qx0; qx <= qx1; qx++) {
				const float ux = (float)(2 * qx + 1 - Nx) * FLOAT_HALF - offX;
				const float vx = (T6FMA(ux, cosa, uy * sina) + ((float)outNx * FLOAT_HALF)) - FLOAT_HALF;
				const float vy = (T6FMA(uy, cosa, -(ux * sina)) + ((float)outNyF * FLOAT_HALF)) - FLOAT_HALF;
				const float wx = FMAX(FLOAT_ONE - FABS(vx - (float)x), FLOAT_ZERO);
				const float wy = FMAX(FLOAT_ONE - FABS(vy - (float)y), FLOAT_ZERO);
				if (wx > FLOAT_ZERO && wy > FLOAT_ZERO) {
					const float w = wx * wy;
					reached = true;
					const LTYPE o = plane + (LTYPE)((qx - srcX) + (qy - srcY) * srcSx);
					for (int j = 0; j < ROT_ZPT; j++) {
						if ((act >> j) & 1)
							val[j] += w * src[o + (LTYPE)j * sxy];
					}
				}
			}
		}
		if (!reached)
			continue;
		for (int j = 0; j < ROT_ZPT; j++) {
			if ((act >> j) & 1)
				acc[j] += val[j];
		}
		anyMask |= act;
	}
	for (int j = 0; j < ROT_ZPT; j++) {
		if ((anyMask >> j) & 1)
			out[n0 + (LTYPE)j * nxy] = acc[j];
	}
}
#endif

// Batched axial pass, see cdrfAxial (the weights are made in shared memory, the radius of the view is at most CDRF_RMAX). The input
// is nb slices, slice b at b * inBS (inSx x inSy x any elements), the output slices at b * outBS (outSx x outSy x any), the
// attenuation, which is the one array for the batch, is read at the offsets of the view (its offset in z includes the slice of the
// view, see T6P_ATTZ0). The groups are (ceil(max N0 / CDRF_TX), max N1, nzb * nb), nzb the number of the z blocks of a view that the
// largest view needs: group z is the z block z % nzb of the view z / nzb. The groups that are outside of the region of their view
// return before the first barrier, which they can since the test is on the group index only
KERNEL void cdrfAxialB(CLGLOBAL float* CLRESTRICT out, const CLGLOBAL float* CLRESTRICT in, const CLGLOBAL float* CLRESTRICT att,
	const CLGLOBAL float* CLRESTRICT bp, const int bpOff, const int nzb, const int attMode, const int attConv,
	const int inSx, const int inSy, const int inBS, const int outSx, const int outSy, const int outBS, const int attSx, const int attSy,
	const float dxv, const float dzv, const float Az, const float Bz, const float C, const float sigmaCoeff) {
	LOCAL float lIn[CDRF_ROWS * CDRF_TX];
	LOCAL float lW[(CDRF_RMAX + 1) * CDRF_TX];
	LOCAL int lR[CDRF_TX];
	const int gz = (int)GRID2;
	const int b = gz / nzb;
	const int zb = gz - b * nzb;
	const CLGLOBAL float* p = bp + ((LTYPE)b + (LTYPE)bpOff) * (LTYPE)T6P_STRIDE;
	const int N0 = T6PI(p, T6P_N0);
	const int NY = T6PI(p, T6P_NY);
	const int y0W = T6PI(p, T6P_Y0);
	const int iz0 = T6PI(p, T6P_IZ0);
	const int iz1 = iz0 + T6PI(p, T6P_INZ);
	const int zStart = T6PI(p, T6P_ZSTART);
	const int nZOut = T6PI(p, T6P_NZOUT);
	const int zo0 = zb * CDRF_TZ;
	if ((int)GRID0 * CDRF_TX >= N0 || (int)GRID1 >= NY || zo0 >= nZOut)
		return;
	const int xOff = T6PI(p, T6P_XOFF);
	const int firstPlane = T6PI(p, T6P_FIRSTPLANE);
	const int r = T6PI(p, T6P_RADIUS);
	const int attX0 = T6PI(p, T6P_ATTX0);
	const int attY0 = T6PI(p, T6P_ATTY0);
	const int attZ0 = T6PI(p, T6P_ATTZ0);
	const float attDx = p[T6P_ATTDX];
	const float dist = p[T6P_DIST];
	const int tid = (int)LID0 + LOCAL_SIZE * (int)LID1;
	const int tx = tid % CDRF_TX;
	const int zl = tid / CDRF_TX;
	const int x = (int)GRID0 * CDRF_TX + tx;
	const int y = (int)GRID1;
	const LTYPE base = (LTYPE)x + (LTYPE)y * (LTYPE)inSx + (LTYPE)b * (LTYPE)inBS;
	const LTYPE sN = (LTYPE)inSx * (LTYPE)inSy;
	const LTYPE oBase = (LTYPE)x + (LTYPE)y * (LTYPE)outSx + (LTYPE)b * (LTYPE)outBS;
	const LTYPE oSn = (LTYPE)outSx * (LTYPE)outSy;
	// The planes that the detector sees have a positive sigma, the rest have all-zero weights
	const float sigma = (x < N0) ? cdrfSigma(false, x + xOff, firstPlane, dist, dxv, dzv, Az, Bz, C, sigmaCoeff, FLOAT_ZERO) : -FLOAT_ONE;
	const float hb = (x < N0) ? cdrfHalfBase(false, x + xOff, firstPlane, dist, dxv, dzv, Az, Bz, C, sigmaCoeff, FLOAT_ZERO) : FLOAT_ZERO;
	const float denom = 2.f * sigma * sigma;
	const bool vis = x < N0 && sigma >= FLOAT_ZERO && denom > FLOAT_ZERO;
	const int rp = vis ? cdrfRadius(sigma, hb, r) : 0;
	if (zl == 0)
		lR[tx] = rp;
	if (vis) {
		for (int k = zl; k <= rp; k += CDRF_NL)
			lW[k * CDRF_TX + tx] = cdrfTap(k, sigma, hb, denom);
	}
	BARRIER
	int rg = 0;
	for (int i = 0; i < CDRF_TX; i++)
		rg = MAX(rg, lR[i]);
	// The normalization: the sum over the taps -rp ... rp in this order
	float wsum = FLOAT_ONE;
	if (vis) {
		wsum = FLOAT_ZERO;
		for (int k = rp; k >= 1; k--)
			wsum += lW[k * CDRF_TX + tx];
		for (int k = 0; k <= rp; k++)
			wsum += lW[k * CDRF_TX + tx];
	}
	// Shared rows 0 ... CDRF_TZ + 2 rg - 1 are the input rows z0 - rg ... z0 + CDRF_TZ + rg - 1
	const int z0 = zStart + zo0;
	const int nRows = CDRF_TZ + 2 * rg;
	if (vis) {
		for (int row = zl; row < nRows; row += CDRF_NL) {
			const int zz = z0 - rg + row;
			float v = FLOAT_ZERO;
			if (zz >= iz0 && zz < iz1) {
				const LTYPE idx = base + (LTYPE)(zz - iz0) * sN;
				v = in[idx];
				if (attMode == 1)
					v *= CDRF_ATT(att, x, y + y0W, zz);
			}
			lIn[row * CDRF_TX + tx] = v;
		}
	}
	// Everybody has summed the taps, so they can be normalized
	BARRIER
	if (vis) {
		for (int k = zl; k <= rp; k += CDRF_NL)
			lW[k * CDRF_TX + tx] = lW[k * CDRF_TX + tx] / wsum;
	}
	BARRIER
	if (x >= N0)
		return;
	float acc[CDRF_TZPT];
	for (int j = 0; j < CDRF_TZPT; j++)
		acc[j] = FLOAT_ZERO;
	if (vis) {
		const int t0 = zl * CDRF_TZPT;
		for (int kk = -rp; kk <= rp; kk++) {
			const float w = lW[(kk < 0 ? -kk : kk) * CDRF_TX + tx];
			const int rowB = t0 + kk + rg;
			for (int j = 0; j < CDRF_TZPT; j++)
				acc[j] += w * lIn[(rowB + j) * CDRF_TX + tx];
		}
	}
	for (int j = 0; j < CDRF_TZPT; j++) {
		const int zo = zo0 + zl * CDRF_TZPT + j;
		if (zo < nZOut) {
			float summa = acc[j];
			if (vis && attMode == 2)
				summa *= CDRF_ATT(att, x, y + y0W, zStart + zo);
			out[oBase + (LTYPE)zo * oSn] = summa;
		}
	}
}

// Batched transaxial pass of the forward projection, see cdrfTransFP. The input is nb slices, slice b at b * inBS (inSx x inSy x any).
// The output is the nOut x N2 projection of the view in the slot T6P_OUTSLOT (view slots of the projections, or the slices of the array
// of the projections that resampleToDetectorB resamples), slot s at s * outSlotStride, the rows of the projection outLd apart
// (outLd = nOut for the projections of the measurements). The groups are (ceil(max nOut / (CDRF_NL * CDRF_YPT)), max N2, nb)
KERNEL void cdrfTransFPB(CLGLOBAL float* CLRESTRICT out, const CLGLOBAL float* CLRESTRICT in, const CLGLOBAL float* CLRESTRICT bp, const int bpOff,
	const int inSx, const int inSy, const int inBS, const int outSlotStride, const int outLd, const int accumulate, const float norm,
	const float dxv, const float dyv, const float Axy, const float Bxy, const float C, const float sigmaCoeff) {
	LOCAL float lRed[CDRF_YPT * CDRF_NT];
	LOCAL float lW[(CDRF_RMAX + 1) * CDRF_TX];
	const int b = (int)GRID2;
	const CLGLOBAL float* p = bp + ((LTYPE)b + (LTYPE)bpOff) * (LTYPE)T6P_STRIDE;
	const float comp = p[T6P_COMP];
	const int N0 = T6PI(p, T6P_N0);
	const int NY = T6PI(p, T6P_NY);
	const int y0W = T6PI(p, T6P_Y0);
	const int N2 = T6PI(p, T6P_NZOUT);
	const int nOut = T6PI(p, T6P_NOUT);
	const int z = (int)GRID1;
	if ((int)GRID0 * (CDRF_NL * CDRF_YPT) >= nOut || z >= N2)
		return;
	const int yStart = T6PI(p, T6P_YSTART);
	const int xOff = T6PI(p, T6P_XOFF);
	const int firstPlane = T6PI(p, T6P_FIRSTPLANE);
	const int r = T6PI(p, T6P_RADIUS);
	const float dist = p[T6P_DIST];
	const LTYPE slot = (LTYPE)T6PI(p, T6P_OUTSLOT);
	const int tid = (int)LID0 + LOCAL_SIZE * (int)LID1;
	const int tx = tid % CDRF_TX;
	const int yl = tid / CDRF_TX;
	const int yo0 = (int)GRID0 * (CDRF_NL * CDRF_YPT) + yl;
	const LTYPE zBase = (LTYPE)z * (LTYPE)inSx * (LTYPE)inSy + (LTYPE)b * (LTYPE)inBS;
	float part[CDRF_YPT];
	for (int q = 0; q < CDRF_YPT; q++)
		part[q] = FLOAT_ZERO;
	for (int x0 = 0; x0 < N0; x0 += CDRF_TX) {
		const int x = x0 + tx;
		const float sigma = (x < N0) ? cdrfSigma(true, x + xOff, firstPlane, dist, dxv, dyv, Axy, Bxy, C, sigmaCoeff, comp) : -FLOAT_ONE;
		const float hb = (x < N0) ? cdrfHalfBase(true, x + xOff, firstPlane, dist, dxv, dyv, Axy, Bxy, C, sigmaCoeff, comp) : FLOAT_ZERO;
		const float denom = 2.f * sigma * sigma;
		const bool vis = x < N0 && sigma >= FLOAT_ZERO && denom > FLOAT_ZERO;
		const int rp = vis ? cdrfRadius(sigma, hb, r) : 0;
		if (vis) {
			for (int k = yl; k <= rp; k += CDRF_NL)
				lW[k * CDRF_TX + tx] = cdrfTap(k, sigma, hb, denom);
		}
		BARRIER
		float wsum = FLOAT_ONE;
		if (vis) {
			wsum = FLOAT_ZERO;
			for (int k = rp; k >= 1; k--)
				wsum += lW[k * CDRF_TX + tx];
			for (int k = 0; k <= rp; k++)
				wsum += lW[k * CDRF_TX + tx];
		}
		BARRIER
		if (vis) {
			for (int k = yl; k <= rp; k += CDRF_NL)
				lW[k * CDRF_TX + tx] = lW[k * CDRF_TX + tx] / wsum;
		}
		BARRIER
		if (vis) {
			const LTYPE base = (LTYPE)x + zBase;
			float summa[CDRF_YPT];
			for (int q = 0; q < CDRF_YPT; q++)
				summa[q] = FLOAT_ZERO;
			for (int kk = -rp; kk <= rp; kk++) {
				const float w = lW[(kk < 0 ? -kk : kk) * CDRF_TX + tx];
				for (int q = 0; q < CDRF_YPT; q++) {
					const int yo = yo0 + q * CDRF_NL;
					const int yy = yStart + yo + kk - y0W;
					if (yo < nOut && yy >= 0 && yy < NY)
						summa[q] += w * in[base + (LTYPE)yy * (LTYPE)inSx];
				}
			}
			for (int q = 0; q < CDRF_YPT; q++)
				part[q] += summa[q];
		}
		// The next chunk overwrites the weights
		BARRIER
	}
	for (int q = 0; q < CDRF_YPT; q++)
		lRed[q * CDRF_NT + tid] = part[q];
	BARRIER
	for (int s = CDRF_TX / 2; s > 0; s >>= 1) {
		if (tx < s) {
			for (int q = 0; q < CDRF_YPT; q++)
				lRed[q * CDRF_NT + tid] += lRed[q * CDRF_NT + tid + s];
		}
		BARRIER
	}
	if (tx == 0) {
		for (int q = 0; q < CDRF_YPT; q++) {
			const int yo = yo0 + q * CDRF_NL;
			if (yo < nOut) {
				const LTYPE o = slot * (LTYPE)outSlotStride + (LTYPE)yo + (LTYPE)z * (LTYPE)outLd;
				const float summa = lRed[q * CDRF_NT + tid] * norm;
				if (accumulate != 0)
					out[o] += summa;
				else
					out[o] = summa;
			}
		}
	}
}

// Batched transaxial pass of the backprojection, see cdrfTransBP (the weights are made in shared memory). The projection of the
// view is the T6P_VIEW'th one of P (pNy x pNz each, at (pY0, pZ0) of the region), or, with resample, it is made from the T6P_VIEW'th of
// the measurements, nRowsD x nColsD each, with the map of the view. The output is nb slices, slice b at b * outBS (outSx x outSy x
// any). The groups are (ceil(max N0 / CDRF_TX), ceil(max N1 / (CDRF_NL * CDRF_YPT)), n2b * nb), n2b the largest N2 of the views: group
// z is the slice z % n2b of the view z / n2b
KERNEL void cdrfTransBPB(CLGLOBAL float* CLRESTRICT out, const CLGLOBAL float* CLRESTRICT P, const CLGLOBAL float* CLRESTRICT bp, const int bpOff,
	const int n2b, const int outSx, const int outSy, const int outBS, const float norm,
	const float dxv, const float dyv, const float Axy, const float Bxy, const float C, const float sigmaCoeff,
	const int resample, const int nRowsD, const int nColsD) {
	LOCAL float lW[(CDRF_RMAX + 1) * CDRF_TX];
	LOCAL float lP[CDRF_NL * CDRF_YPT + 2 * CDRF_RMAX];
	LOCAL int lR[CDRF_TX];
	const int gz = (int)GRID2;
	const int b = gz / n2b;
	const int z = gz - b * n2b;
	const CLGLOBAL float* p = bp + ((LTYPE)b + (LTYPE)bpOff) * (LTYPE)T6P_STRIDE;
	const float comp = p[T6P_COMP];
	const int N0 = T6PI(p, T6P_N0);
	const int N1 = T6PI(p, T6P_N1);
	const int NY = T6PI(p, T6P_NY);
	const int y0W = T6PI(p, T6P_Y0);
	const int INZ = T6PI(p, T6P_INZ);
	const int zAbs = T6PI(p, T6P_IZ0) + z;
	const int yg0 = (int)GRID1 * (CDRF_NL * CDRF_YPT);
	if ((int)GRID0 * CDRF_TX >= N0 || yg0 >= NY || z >= INZ)
		return;
	const int pNy = T6PI(p, T6P_PNY);
	const int pNz = T6PI(p, T6P_PNZ);
	const int pY0 = T6PI(p, T6P_PY0);
	const int pZ0 = T6PI(p, T6P_PZ0);
	const int view = T6PI(p, T6P_VIEW);
	const int xOff = T6PI(p, T6P_XOFF);
	const int firstPlane = T6PI(p, T6P_FIRSTPLANE);
	const int r = T6PI(p, T6P_RADIUS);
	const float dist = p[T6P_DIST];
	const float trA = p[T6P_TRA];
	const float trB = p[T6P_TRB];
	const float axA = p[T6P_AXA];
	const float axB = p[T6P_AXB];
	const int tid = (int)LID0 + LOCAL_SIZE * (int)LID1;
	const int tx = tid % CDRF_TX;
	const int yl = tid / CDRF_TX;
	const int x = (int)GRID0 * CDRF_TX + tx;
	const int yb = yg0 + yl * CDRF_YPT;
	const int pz = zAbs - pZ0;
	const LTYPE zOut = (LTYPE)z * (LTYPE)outSx * (LTYPE)outSy + (LTYPE)b * (LTYPE)outBS;
	if (pz < 0 || pz >= pNz) {
		// The projection does not reach this z, which holds for the whole group
		if (x < N0) {
			for (int j = 0; j < CDRF_YPT; j++) {
				const int y = yb + j;
				if (y < NY)
					out[(LTYPE)x + (LTYPE)y * (LTYPE)outSx + zOut] = FLOAT_ZERO;
			}
		}
		return;
	}
	const float sigma = (x < N0) ? cdrfSigma(true, x + xOff, firstPlane, dist, dxv, dyv, Axy, Bxy, C, sigmaCoeff, comp) : -FLOAT_ONE;
	const float hb = (x < N0) ? cdrfHalfBase(true, x + xOff, firstPlane, dist, dxv, dyv, Axy, Bxy, C, sigmaCoeff, comp) : FLOAT_ZERO;
	const float denom = 2.f * sigma * sigma;
	const bool vis = x < N0 && sigma >= FLOAT_ZERO && denom > FLOAT_ZERO;
	const int rp = vis ? cdrfRadius(sigma, hb, r) : 0;
	if (yl == 0)
		lR[tx] = rp;
	if (vis) {
		for (int k = yl; k <= rp; k += CDRF_NL)
			lW[k * CDRF_TX + tx] = cdrfTap(k, sigma, hb, denom);
	}
	BARRIER
	int rg = 0;
	for (int i = 0; i < CDRF_TX; i++)
		rg = MAX(rg, lR[i]);
	float wsum = FLOAT_ONE;
	if (vis) {
		wsum = FLOAT_ZERO;
		for (int k = rp; k >= 1; k--)
			wsum += lW[k * CDRF_TX + tx];
		for (int k = 0; k <= rp; k++)
			wsum += lW[k * CDRF_TX + tx];
	}
	// Shared elements 0 ... CDRF_NL * CDRF_YPT + 2 rg - 1 are the rows yg0 - rg ... of the projection, zero
	// where the projection is not
	const int nW = CDRF_NL * CDRF_YPT + 2 * rg;
	const LTYPE pBase = (LTYPE)view * (LTYPE)pNy * (LTYPE)pNz + (LTYPE)pz * (LTYPE)pNy;
	const LTYPE mBase = (LTYPE)view * (LTYPE)nRowsD * (LTYPE)nColsD;
	for (int i = tid; i < nW; i += CDRF_NT) {
		const int yy = y0W + yg0 - rg + i;
		float v = FLOAT_ZERO;
		if (yy >= 0 && yy < N1) {
			if (resample != 0) {
				v = cdrfFromDetectorAt(P, mBase, yy, zAbs, nRowsD, nColsD, trA, trB, axA, axB);
			}
			else {
				const int py = yy - pY0;
				if (py >= 0 && py < pNy)
					v = P[pBase + (LTYPE)py];
			}
		}
		lP[i] = v;
	}
	BARRIER
	if (vis) {
		for (int k = yl; k <= rp; k += CDRF_NL)
			lW[k * CDRF_TX + tx] = lW[k * CDRF_TX + tx] / wsum;
	}
	BARRIER
	if (x >= N0 || yb >= NY)
		return;
	float acc[CDRF_YPT];
	for (int j = 0; j < CDRF_YPT; j++)
		acc[j] = FLOAT_ZERO;
	if (vis) {
		// Row yb + j + kk is the shared element yb + j + kk - (yg0 - rg), which the zeros of the window make
		// the same as leaving out the rows that are not in the projection
		const int e0 = yl * CDRF_YPT + rg;
		for (int kk = -rp; kk <= rp; kk++) {
			const float w = lW[(kk < 0 ? -kk : kk) * CDRF_TX + tx];
			for (int j = 0; j < CDRF_YPT; j++)
				acc[j] += w * lP[e0 + j + kk];
		}
	}
	for (int j = 0; j < CDRF_YPT; j++) {
		const int y = yb + j;
		if (y < NY)
			out[(LTYPE)x + (LTYPE)y * (LTYPE)outSx + zOut] = vis ? acc[j] * norm : FLOAT_ZERO;
	}
}

// Batched resampling of the projections onto the detector, see resampleToDetector. P is nb slices, slice b (nTr x nAx of them used,
// the rows pLd apart) at b * pBS; the result of the view is added to its T6P_VIEW'th projection of fProj (slots of fSlotStride
// elements, nRowsD x nColsD). The launch is (nRowsD, nColsD, nb) threads
KERNEL void resampleToDetectorB(CLGLOBAL float* CLRESTRICT fProj, const CLGLOBAL float* CLRESTRICT P, const CLGLOBAL float* CLRESTRICT bp, const int bpOff,
	const int nRowsD, const int nColsD, const int pLd, const int pBS, const int fSlotStride) {
	const int r = GID0;
	const int c = GID1;
	const int b = (int)GID2;
	if (r >= nRowsD || c >= nColsD)
		return;
	const CLGLOBAL float* p = bp + ((LTYPE)b + (LTYPE)bpOff) * (LTYPE)T6P_STRIDE;
	const int nTr = T6PI(p, T6P_NTR);
	const int nAx = T6PI(p, T6P_NAX);
	// The rows and the axial rows of the projection that P holds start at the T6P_YSTART'th and the T6P_ZSTART'th of the lattice
	const int yS = T6PI(p, T6P_YSTART);
	const int zS = T6PI(p, T6P_ZSTART);
	const float trA = p[T6P_TRA];
	const float trB = p[T6P_TRB];
	const float axA = p[T6P_AXA];
	const float axB = p[T6P_AXB];
	const LTYPE view = (LTYPE)T6PI(p, T6P_VIEW);
	const float tt = trA * CFLOAT(r) + trB;
	const float ta = axA * CFLOAT(c) + axB;
	const int i0 = (int)FLOOR(tt);
	const int j0 = (int)FLOOR(ta);
	float summa = FLOAT_ZERO;
	for (int jj = 0; jj < 2; jj++) {
		const int j = j0 + jj;
		if (j < zS || j >= zS + nAx)
			continue;
		const float wj = FMAX(FLOAT_ONE - FABS(ta - CFLOAT(j)), FLOAT_ZERO);
		for (int ii = 0; ii < 2; ii++) {
			const int i = i0 + ii;
			if (i < yS || i >= yS + nTr)
				continue;
			const float wi = FMAX(FLOAT_ONE - FABS(tt - CFLOAT(i)), FLOAT_ZERO);
			summa += wi * wj * P[(LTYPE)b * (LTYPE)pBS + (LTYPE)(i - yS) + (LTYPE)(j - zS) * (LTYPE)pLd];
		}
	}
	fProj[view * (LTYPE)fSlotStride + (LTYPE)r + (LTYPE)c * (LTYPE)nRowsD] += summa;
}
#endif

// Initial version from: https://stackoverflow.com/questions/9833316/cuda-image-rotation/10008412#10008412
// offX and offY shift the image, in voxels, at the same time as it is rotated, i.e. the output is
// rotim(q) = im(R * q + (offX, offY)). The shifts may be non-integer. This is required by the SPECT
// rotation-based projector with detector heads that can swivel, in which case the detector panel is
// not centered with respect to the rotation axis.
//
// The output does not have to be sampled on the input grid. Nx and Ny are the dimensions of the input
// while outNx and outNyF are those of the output, both of them centered on their own center, and the
// two grids are assumed to have the same voxel size. This is what the multi-resolution volumes of the
// projector need, since they are rotated into a grid that is centered on the rotation axis of the
// gantry while the volume itself is not. With outNx == Nx and outNyF == Ny the input grid is also the
// output grid, which is the case without multi-resolution
//
// Only a band of the output is written and only a band of the input is read, which is what the
// projector needs when the projections cover just the detector instead of the whole FOV:
//   outY, outZ, outNy, outNz  the part of dimensions 1 and 2 that is written, the output array has the
//                             dimensions outNx x outNy x outNz
//   srcY, srcZ, srcNy, srcNz  the part of dimensions 1 and 2 that holds data, everything else reads as
//                             zero. With images the source sits at its natural place inside a full
//                             sized image, with buffers the source array holds only the band
// The x dimension (the depth of the detector) has a window of its own, like the other two:
//   outX, outNxW  the part of dimension 0 that is written: the output array is outNxW x outNy x outNz and
//                 its x index 0 is x index outX of the output grid, whose full width outNx (like outNyF for y)
//                 is the one the output is centered on
//   srcX, srcNxB  the part of dimension 0 of the source that holds data, everything else reads as zero.
//                 With images the source sits at its natural place inside a full sized image, with buffers
//                 the source array holds only the band (srcNxB wide, indexed relative to srcX)
// With accumulate set the result is added to rotim instead of replacing it, and the elements that the
// rotation does not reach (outside of the source) are neither read nor written, i.e. they keep their value.
// This is what the backprojection uses to add a rotated view directly into the image. It is not available
// on METAL. Every element is handled by a single thread, so no atomics are needed
// With accumulate set the rotation is, by default, not a gather at the position of the inverse rotation but the exact transpose (adjoint)
// of the rotation that the forward projection applies (the same bilinear weights, summed over the elements that sample the voxel), so
// that the forward and the backprojection of the projector are an adjoint pair (the dot test <A x, y> = <x, A^T y> holds to the float
// precision). Each voxel then sums at most 9 elements of the source (3 on average, which is all that has a weight), but it is still handled by
// a single thread, and only the voxels that an element reaches are written. A thread handles ROT_ZPT slices of a column of the volume (the
// launch is ROT_ZPT times shorter in z), since the weights do not depend on z. With TYPE6_BP_FAST (the host option type6BPRotation =
// 'fast') the rotation back is the cheaper gather at the position of the inverse rotation, which is only an approximate transpose, and
// every thread handles one voxel
// outSx and outSy are the dimensions that the output array is allocated with (its row and slice strides), outNxW and outNy when it is
// not allocated larger than the window; the elements outside of the window are not written
// Setting both windows to the whole volume gives a plain rotation of the entire image. On METAL, which
// keeps the previous, plain (non-multires-aware) rotation, the extra parameters are simply given their
// values for that case: outNx=Nx, outNyF=Ny, offX=offY=0, outY=outZ=0, outNy=Ny, outNz=Nz, srcY=srcZ=0,
// srcNy=Ny, srcNz=Nz, outX=0, outNxW=Nx, srcX=0, srcNxB=Nx
KERNEL void rotate(
	CLGLOBAL float* CLRESTRICT rotim BUF0,
	IMTYPE im TEX1,
#ifdef METAL
	SCALAR_PARAMS(scalarParams) BUF2,
	uint3 metalGlobalId [[thread_position_in_grid]]
#else
	const int Nx, const int Ny, const int Nz, const int outNx, const int outNyF,
	const float cosa, const float sina, const float offX, const float offY,
	const int outY, const int outZ, const int outNy, const int outNz,
	const int srcY, const int srcZ, const int srcNy, const int srcNz, const int accumulate,
	const int outX, const int outNxW, const int srcX, const int srcNxB,
	const int outSx, const int outSy
#endif
) {
#ifdef METAL
	UNPACK_SCALAR_PARAMS_ROTATE(scalarParams)
	const int outNx = Nx;
	const int outNyF = Ny;
	const float offX = FLOAT_ZERO;
	const float offY = FLOAT_ZERO;
	const int outY = 0;
	const int outZ = 0;
	const int outNy = Ny;
	const int outNz = Nz;
	const int srcY = 0;
	const int srcZ = 0;
	const int srcNy = Ny;
	const int srcNz = Nz;
	const int outX = 0;
	const int outNxW = Nx;
	const int srcX = 0;
	const int srcNxB = Nx;
	const int outSx = Nx;
	const int outSy = Ny;
	LTYPE3 xyz = MINT3(metalGlobalId.x, metalGlobalId.y, metalGlobalId.z);
#else
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#endif
#if defined(METAL) || defined(TYPE6_BP_FAST)
	if (xyz.x >= outNxW || xyz.y >= outNy || xyz.z >= outNz)
		return;
#else
	// With accumulate the z range of the launch is in blocks of ROT_ZPT slices
	if (xyz.x >= outNxW || xyz.y >= outNy || (accumulate != 0 ? xyz.z * (LTYPE)ROT_ZPT : xyz.z) >= outNz)
		return;
#endif
	const LTYPE n = (xyz.x) + (xyz.y) * (LTYPE)(outSx) + (xyz.z) * ((LTYPE)outSx * (LTYPE)outSy);
	// Position in the full volume
	const int x = outX + xyz.x;
	const int y = outY + xyz.y;
	const int z = outZ + xyz.z;

    // The rotation centre of a grid of n elements is its exact centre, n / 2 (real division): element x is at x + 0.5 - n / 2 from it, i.e. (2 x + 1 - n) / 2,
    // and the centre of the source in the source index is Nx / 2 - 0.5. For the odd n it is the middle of the middle element (not the boundary n / 2 (integer
    // division) of before, which was half an element off the axis that the host assumes, see offsetsFPType6, depthWindowType6, type6RotCentre in precomp.h)
    const float xA = (float)(2 * x + 1 - outNx) * FLOAT_HALF;
    const float yA = (float)(2 * y + 1 - outNyF) * FLOAT_HALF;

#ifdef METAL
    const float src_x = (xA * cosa - yA * sina + offX + (float)Nx * FLOAT_HALF) - FLOAT_HALF;
    const float src_y = (xA * sina + yA * cosa + offY + (float)Ny * FLOAT_HALF) - FLOAT_HALF;
#else
    const float src_x = (T6FMA(xA, cosa, -(yA * sina)) + offX + (float)Nx * FLOAT_HALF) - FLOAT_HALF;
    const float src_y = (T6FMA(yA, cosa, xA * sina) + offY + (float)Ny * FLOAT_HALF) - FLOAT_HALF;
#endif

#ifdef METAL
	// Metal keeps the previous, plain rotation
	// The interpolation reads one row past src_y, so with a partial source the last row is left out.
	// There the values are below the truncation level of the collimator response anyway
	const float srcYlo = CFLOAT(srcY);
	const float srcYhi = (srcY + srcNy >= Ny) ? CFLOAT(Ny) : CFLOAT(srcY + srcNy - 1);

    if (src_x >= 0.0f && src_x < Nx && src_y >= srcYlo && src_y < srcYhi && z >= srcZ && z < srcZ + srcNz) {
		float val = FLOAT_ZERO;
		// The interpolation reads one element past the sampling position, which in the outermost row
		// and column of the source is outside it. Clamping the position there repeats the outermost
		// value, which is what an image with the dimensions of the source would do on its own. The
		// image is shared by every volume and is therefore larger than the source, so without this the
		// read would return whatever the previous, larger, upload left in it. With buffers the read
		// would wrap around to the next row
		const float sampX = (src_x > CFLOAT(Nx - 1)) ? CFLOAT(Nx - 1) : src_x;
		const float sampY = (src_y > CFLOAT(srcY + srcNy - 1)) ? CFLOAT(srcY + srcNy - 1) : src_y;
#ifdef USEIMAGES
		val = im.sample(samplerRotate, float3(sampX + FLOAT_HALF, sampY + FLOAT_HALF, CFLOAT(z) + FLOAT_HALF)).r;
#else
        // BILINEAR INTERPOLATION
        // The buffer holds only the band, so the source is indexed relative to its origin
        const float src_yL = sampY - CFLOAT(srcY);
        const int src_x0 = (int)(sampX);
        const int src_x1 = (src_x0 + 1);
        const int src_y0 = (int)(src_yL);
        const int src_y1 = (src_y0 + 1);

        const float sx = (sampX - src_x0);
        const float sy = (src_yL - src_y0);

        const int idx_src00 = MIN(MAX(0, src_x0 + src_y0 * Nx), (Nx * srcNy) - 1);
        const int idx_src10 = MIN(MAX(0, src_x1 + src_y0 * Nx), (Nx * srcNy) - 1);
        const int idx_src01 = MIN(MAX(0, src_x0 + src_y1 * Nx), (Nx * srcNy) - 1);
        const int idx_src11 = MIN(MAX(0, src_x1 + src_y1 * Nx), (Nx * srcNy) - 1);
        const LTYPE plane = (LTYPE)(z - srcZ) * (LTYPE)Nx * (LTYPE)srcNy;

        val  = (FLOAT_ONE - sx) * (FLOAT_ONE - sy) * im[plane + idx_src00];
        val += (       sx) * (FLOAT_ONE - sy) * im[plane + idx_src10];
        val += (FLOAT_ONE - sx) * (       sy) * im[plane + idx_src01];
        val += (       sx) * (       sy) * im[plane + idx_src11];
#endif
		rotim[n] = val;
    }
	else {
        rotim[n] = 0.0f;
    }
#else
#ifndef TYPE6_BP_FAST
	if (accumulate != 0) {
		// The backprojection (the rotation back of a view into the image): the exact transpose of the bilinear rotation of the forward
		// projection, which samples the image at the position P(q) = R^T (q - off) (q the element of the source here, the region of the
		// view, both centered on their own grids) with the weight tri(P(q) - voxel). The voxel therefore receives sum_q tri(P(q) - voxel)
		// * source(q) over the elements q whose sample position is within one voxel of it. Those are the q strictly within |cosa| + |sina|
		// of src_x, src_y (the position that the inverse rotation would sample), the others have no weight. That inverse rotation (the
		// interpolation at src_x, src_y) is not the transpose of the interpolation of the forward projection, which made the two
		// mismatched for everything that is not smooth. Elements outside of the source have no weight. The weights do not depend on z, so
		// they are made once for the ROT_ZPT slices z0 ... z0 + ROT_ZPT - 1 that the thread handles
		const int tLoX = srcX;
		const int tHiX = srcX + srcNxB - 1;
		const int tLoY = srcY;
		const int tHiY = srcY + srcNy - 1;
		const int z0 = outZ + (int)xyz.z * ROT_ZPT;
		int act = 0;
		for (int j = 0; j < ROT_ZPT; j++) {
			const int zj = z0 + j;
			if (zj < outZ + outNz && zj >= srcZ && zj < srcZ + srcNz)
				act |= 1 << j;
		}
		if (act == 0)
			return;
		const float hh = FABS(cosa) + FABS(sina);
		const int qx0 = MAX((int)FLOOR(src_x - hh) + 1, tLoX);
		const int qx1 = MIN((int)CEIL(src_x + hh) - 1, tHiX);
		const int qy0 = MAX((int)FLOOR(src_y - hh) + 1, tLoY);
		const int qy1 = MIN((int)CEIL(src_y + hh) - 1, tHiY);
		float val[ROT_ZPT];
		for (int j = 0; j < ROT_ZPT; j++)
			val[j] = FLOAT_ZERO;
		// Only the voxels that an element of the source reaches are written (the others would only get a zero added, which costs a
		// read and a write of the whole volume for every view, while the band that the source covers is a fraction of it)
		bool reached = false;
#ifndef USEIMAGES
		const LTYPE pl = (LTYPE)srcNxB * (LTYPE)srcNy;
#endif
		for (int qy = qy0; qy <= qy1; qy++) {
			const float uy = (float)(2 * qy + 1 - Ny) * FLOAT_HALF - offY;
			for (int qx = qx0; qx <= qx1; qx++) {
				const float ux = (float)(2 * qx + 1 - Nx) * FLOAT_HALF - offX;
				const float vx = (T6FMA(ux, cosa, uy * sina) + ((float)outNx * FLOAT_HALF)) - FLOAT_HALF;
				const float vy = (T6FMA(uy, cosa, -(ux * sina)) + ((float)outNyF * FLOAT_HALF)) - FLOAT_HALF;
				const float wx = FMAX(FLOAT_ONE - FABS(vx - (float)x), FLOAT_ZERO);
				const float wy = FMAX(FLOAT_ONE - FABS(vy - (float)y), FLOAT_ZERO);
				if (wx > FLOAT_ZERO && wy > FLOAT_ZERO) {
					const float w = wx * wy;
					reached = true;
#ifdef USEIMAGES
					const float xc = CFLOAT(qx) + FLOAT_HALF;
					const float yc = CFLOAT(qy) + FLOAT_HALF;
					for (int j = 0; j < ROT_ZPT; j++) {
						if ((act >> j) & 1) {
							const float zcT = CFLOAT(z0 + j) + FLOAT_HALF;
#if defined(CUDA) || defined(HIP)
							val[j] += w * tex3D<float>(im, xc, yc, zcT);
#elif defined(OPENCL)
							val[j] += w * read_imagef(im, samplerRotate, (float4)(xc, yc, zcT, FLOAT_ZERO)).w;
#endif
						}
					}
#else
					const LTYPE o = (LTYPE)(z0 - srcZ) * pl + (LTYPE)((qx - srcX) + (qy - srcY) * srcNxB);
					for (int j = 0; j < ROT_ZPT; j++) {
						if ((act >> j) & 1)
							val[j] += w * im[o + (LTYPE)j * pl];
					}
#endif
				}
			}
		}
		if (!reached)
			return;
		const LTYPE nBase = (xyz.x) + (xyz.y) * (LTYPE)(outSx) + (LTYPE)(z0 - outZ) * ((LTYPE)outSx * (LTYPE)outSy);
		for (int j = 0; j < ROT_ZPT; j++) {
			if ((act >> j) & 1)
				rotim[nBase + (LTYPE)j * ((LTYPE)outSx * (LTYPE)outSy)] += val[j];
		}
		return;
	}
#endif
	// The source is bilinearly interpolated with zero padding: the valid voxels are x in [xLo, xHi] and
	// y in [yLo, yHi], which is the band for a partial source, and everything outside of them is zero.
	// This way a sample that falls between the outermost voxel and the first voxel outside the source
	// fades linearly to zero instead of repeating the outermost value, which keeps the footprint of a
	// voxel at the edge of the volume, or of the band, identical to that of an interior voxel. The
	// z dimension is not interpolated, so it is just a range check
	const int xLo = srcX;
	const int xHi = srcX + srcNxB - 1;
	const int yLo = srcY;
	const int yHi = srcY + srcNy - 1;

	if (src_x > CFLOAT(xLo) - FLOAT_ONE && src_x < CFLOAT(xHi) + FLOAT_ONE && src_y > CFLOAT(yLo) - FLOAT_ONE && src_y < CFLOAT(yHi) + FLOAT_ONE
		&& z >= srcZ && z < srcZ + srcNz) {
		float val = FLOAT_ZERO;
#ifdef USEIMAGES
		if (src_x >= CFLOAT(xLo) && src_x <= CFLOAT(xHi) && src_y >= CFLOAT(yLo) && src_y <= CFLOAT(yHi)) {
			// Interior, all the voxels the filter reads are inside the source (or have zero weight)
#if defined(CUDA) || defined(HIP)
			val = tex3D<float>(im, src_x + FLOAT_HALF, src_y + FLOAT_HALF, CFLOAT(z) + FLOAT_HALF);
#elif defined(OPENCL)
			val = read_imagef(im, samplerRotate, (float4)(src_x + FLOAT_HALF, src_y + FLOAT_HALF, CFLOAT(z) + FLOAT_HALF, FLOAT_ZERO)).w;
#endif
		}
		else {
			// Edge, the neighbours outside the source contribute nothing. Sampling at a texel center returns
			// the value of that texel as such
			const int x0 = (int)FLOOR(src_x);
			const int y0 = (int)FLOOR(src_y);
			const float fx = src_x - CFLOAT(x0);
			const float fy = src_y - CFLOAT(y0);
			const float zc = CFLOAT(z) + FLOAT_HALF;
			for (int jj = 0; jj < 2; jj++) {
				const int yy = y0 + jj;
				if (yy < yLo || yy > yHi)
					continue;
				const float wy = (jj == 0) ? (FLOAT_ONE - fy) : fy;
				for (int ii = 0; ii < 2; ii++) {
					const int xx = x0 + ii;
					if (xx < xLo || xx > xHi)
						continue;
					const float w = ((ii == 0) ? (FLOAT_ONE - fx) : fx) * wy;
#if defined(CUDA) || defined(HIP)
					val += w * tex3D<float>(im, CFLOAT(xx) + FLOAT_HALF, CFLOAT(yy) + FLOAT_HALF, zc);
#elif defined(OPENCL)
					val += w * read_imagef(im, samplerRotate, (float4)(CFLOAT(xx) + FLOAT_HALF, CFLOAT(yy) + FLOAT_HALF, zc, FLOAT_ZERO)).w;
#endif
				}
			}
		}
#else
        // BILINEAR INTERPOLATION
        // The buffer holds only the band, so the source is indexed relative to its origin. Neighbours
        // outside of the source have zero weight and are not read, which for interior samples is the
        // same as reading all four of them
		const int x0 = (int)FLOOR(src_x);
		const int y0 = (int)FLOOR(src_y);
		const float fx = src_x - CFLOAT(x0);
		const float fy = src_y - CFLOAT(y0);
		const LTYPE plane = (LTYPE)(z - srcZ) * (LTYPE)srcNxB * (LTYPE)srcNy;
		for (int jj = 0; jj < 2; jj++) {
			const int yy = y0 + jj;
			if (yy < yLo || yy > yHi)
				continue;
			const float wy = (jj == 0) ? (FLOAT_ONE - fy) : fy;
			for (int ii = 0; ii < 2; ii++) {
				const int xx = x0 + ii;
				if (xx < xLo || xx > xHi)
					continue;
				const float w = ((ii == 0) ? (FLOAT_ONE - fx) : fx) * wy;
				val += w * im[plane + (LTYPE)((xx - srcX) + (yy - srcY) * srcNxB)];
			}
		}
#endif
		if (accumulate != 0)
			rotim[n] += val;
		else
			rotim[n] = val;
    }
	else if (accumulate == 0) {
        rotim[n] = 0.0f;
    }
#endif

}
#endif
