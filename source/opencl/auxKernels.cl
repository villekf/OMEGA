
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
#if defined(NLM_) || defined(NLM4_)
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
#define KOKO ((SEARCH_WINDOW_X * 2 + 1) * (SEARCH_WINDOW_Y * 2 + 1) * (SEARCH_WINDOW_Z * 2 + 1))
// Median selection. Windows of at most MED_NET_LIMIT voxels use a fully unrolled min/max pruning network (the window stays in registers),
// larger windows use a bisection on the order-preserving bit pattern of the values (32 read-only counting passes).
// Both give exactly the same median as a ranking search. The best limit depends on the compiler:
// NVIDIA RTX 4090 (CUDA and OpenCL, 512^3): the network is about 10x faster than a ranking search at 25-27 voxels, but slower than it
// from 45 voxels on, where the bisection is 3-90x faster.
// AMD Radeon RX 7900 XT (HIP, 256^3): the network is the fastest method at every tested size up to 147 voxels (2-6x faster than the bisection).
#ifndef MED_NET_LIMIT
#ifdef HIP
#define MED_NET_LIMIT 147
#else
#define MED_NET_LIMIT 27
#endif
#endif
#if KOKO <= MED_NET_LIMIT
#define MED_SWAP(a, b) { const float medTmp = fmin(a, b); b = fmax(a, b); a = medTmp; }
#else
#ifdef OPENCL
#define MED_F2U(x) as_uint(x)
#define MED_U2F(x) as_float(x)
#else
#define MED_F2U(x) __float_as_uint(x)
#define MED_U2F(x) __uint_as_float(x)
#endif
// Order-preserving key: non-negative floats get the sign bit set, negative floats are bit-inverted
#define MED_KEY(u) (((u) & 0x80000000u) ? ~(u) : ((u) | 0x80000000u))
#define MED_UNKEY(k) (((k) & 0x80000000u) ? MED_U2F((k) ^ 0x80000000u) : MED_U2F(~(k)))
#endif
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
	d_im[idx] = apu / (convert_float(d_Summ[idx]) / TH + d_epps) * (convert_float(d_rhs[idx]) / TH + d_epps);
#else
	d_im[idx] = apu / (d_Summ[idx] + d_epps) * (d_rhs[idx] + d_epps);
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
				{
					const LTYPE xc = MAX((LTYPE)0, MIN(xx, (LTYPE)(N.x - 1)));
					const LTYPE yc = MAX((LTYPE)0, MIN(yy, (LTYPE)(N.y - 1)));
					const LTYPE zc = MAX((LTYPE)0, MIN(zz, (LTYPE)(N.z - 1)));
					lCacheRef[indX + indY * NLM_TILEX + indZ * NLM_TILEX * NLM_TILEY] = u_ref[xc + yc * N.x + zc * N.x * N.y];
				}
#endif
#endif // END NLMREF
#ifdef USEIMAGES
#if defined(CUDA) || defined(HIP)
				lCache[indX + indY * NLM_TILEX + indZ * NLM_TILEX * NLM_TILEY] = tex3D<float>(u, xx, yy, zz);
#else
				lCache[indX + indY * NLM_TILEX + indZ * NLM_TILEX * NLM_TILEY] = read_imagef(u, samplerNLM, (int4)(xx, yy, zz, 0)).w;
#endif
#else
				{
					const LTYPE xc = MAX((LTYPE)0, MIN(xx, (LTYPE)(N.x - 1)));
					const LTYPE yc = MAX((LTYPE)0, MIN(yy, (LTYPE)(N.y - 1)));
					const LTYPE zc = MAX((LTYPE)0, MIN(zz, (LTYPE)(N.z - 1)));
					lCache[indX + indY * NLM_TILEX + indZ * NLM_TILEX * NLM_TILEY] = u[xc + yc * N.x + zc * N.x * N.y];
				}
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
	const int3 xxyyzz = CMINT3(LID0 + SWINDOWX + PWINDOWX,LID1 + SWINDOWY + PWINDOWY, LID2 + SWINDOWZ + PWINDOWZ);
	const float output = NLMGradient(lCache, gaussian, xxyyzz.x, xxyyzz.y, xxyyzz.z, h, epps
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
				lCacheRef[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = SWBufRead(u_ref, (int)xx, (int)yy, (int)zz, N);
#endif
#endif // END RDPREF
#ifdef USEIMAGES
#if defined(CUDA) || defined(HIP)
				lCache[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = tex3D<float>(u, xx, yy, zz);
#else
				lCache[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = read_imagef(u, samplerRDP, (int4)(xx, yy, zz, 0)).w;
#endif
#else
				lCache[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = SWBufRead(u, (int)xx, (int)yy, (int)zz, N);
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
				lCache[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = SWBufRead(u, (int)xx, (int)yy, (int)zz, N);
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
#if KOKO <= MED_NET_LIMIT
	float median[KOKO];
#else
	uint median[KOKO];
#endif
	int uu = 0;
	for (LTYPE x = -SEARCH_WINDOW_X; x <= SEARCH_WINDOW_X; x++) {
		for (LTYPE y = -SEARCH_WINDOW_Y; y <= SEARCH_WINDOW_Y; y++) {
			for (LTYPE z = -SEARCH_WINDOW_Z; z <= SEARCH_WINDOW_Z; z++) {
				LTYPE pikseli = (xyz.x + (x)) + (xyz.y + (y)) * (N.x + SEARCH_WINDOW_X * 2) + (xyz.z + (z)) * (N.x + SEARCH_WINDOW_X * 2) * (N.y + SEARCH_WINDOW_Y * 2);
#if KOKO <= MED_NET_LIMIT
				median[uu] = grad[pikseli];
#else
				const uint medBits = MED_F2U(grad[pikseli]);
				median[uu] = MED_KEY(medBits);
#endif
				uu++;
			}
		}
	}
#if KOKO <= MED_NET_LIMIT
	// Repeatedly find the smallest and largest of KOKO / 2 + 2 values, drop both and take in one new value; the last middle value is the median
	const int medE = KOKO / 2 + 1;
#pragma unroll
	for (int s = 0; s < KOKO / 2; s++) {
		const int nAct = medE - s + 1;
#pragma unroll
		for (int i = 0; i < nAct / 2; i++)
			MED_SWAP(median[s + i], median[medE - i]);
#pragma unroll
		for (int i = 1; i < (nAct + 1) / 2; i++)
			MED_SWAP(median[s], median[s + i]);
#pragma unroll
		for (int i = nAct / 2; i < nAct - 1; i++)
			MED_SWAP(median[s + i], median[medE]);
		if (s < KOKO / 2 - 1)
			median[medE] = median[medE + 1 + s];
	}
	output[n] = median[KOKO / 2];
#else
	// Build the key of rank KOKO / 2 bit by bit: a bit is kept when at most KOKO / 2 keys are smaller than the candidate
	uint medAns = 0u;
	for (int bb = 31; bb >= 0; bb--) {
		const uint cand = medAns | (1u << bb);
		int cnt = 0;
		for (int ll = 0; ll < KOKO; ll++)
			cnt += (median[ll] < cand) ? 1 : 0;
		medAns = (cnt <= KOKO / 2) ? cand : medAns;
	}
	output[n] = MED_UNKEY(medAns);
#endif
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
// Reads the dual variable of the proximal TV/TGV. The dual variable exists only in the volume where the prior is computed
// (the original FOV) and is zero elsewhere
DEVICE float proxQRead(const CLGLOBAL float* CLRESTRICT q, const LTYPE x, const LTYPE y, const LTYPE z, const int3 N) {
	if (x < 0 || y < 0 || z < 0 || x >= N.x || y >= N.y || z >= N.z)
		return FLOAT_ZERO;
	return q[x + y * N.x + z * N.x * N.y];
}

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
	const LTYPE x = xyz.x + xyz.y * N.x + xyz.z * N.x * N.y;
#if defined(EFOVZ)
	// Coordinates in the original FOV, i.e. in the volume where the prior is computed
	const LTYPE3 NDiff = (N - NOrig) / 2;
	const LTYPE3 xyzF = xyz - NDiff;
#else
	const LTYPE3 xyzF = xyz;
#endif
	// The prior is computed only in the original FOV: voxels of the extended FOV are skipped
	// Note that the gradient of the voxels at the edge of the original FOV uses the neighboring voxel of the extended FOV, but
	// the divergence does not include that voxel. An exact adjoint, which would include it, was tested and removes the anchoring
	// effect of the extended FOV on the edge of the original FOV, degrading the edge slices with a strong regularization
	if (xyzF.x < 0 || xyzF.y < 0 || xyzF.z < 0 || xyzF.x >= NOrig.x || xyzF.y >= NOrig.y || xyzF.z >= NOrig.z)
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
	float apuVal = FLOAT_ZERO;
// Transpose of forward difference (backward difference)
#if DIFFTYPE == 0 // START DIFFTYPE == 0
	apuVal += (proxQRead(gradX, xyzF.x, xyzF.y, xyzF.z, NOrig) - proxQRead(gradX, xyzF.x - 1, xyzF.y, xyzF.z, NOrig));
	apuVal += (proxQRead(gradY, xyzF.x, xyzF.y, xyzF.z, NOrig) - proxQRead(gradY, xyzF.x, xyzF.y - 1, xyzF.z, NOrig));
	apuVal += (proxQRead(gradZ, xyzF.x, xyzF.y, xyzF.z, NOrig) - proxQRead(gradZ, xyzF.x, xyzF.y, xyzF.z - 1, NOrig));
// Transpose of backward difference (forward difference)
#elif DIFFTYPE == 1  // START DIFFTYPE == 1
	apuVal += (proxQRead(gradX, xyzF.x + 1, xyzF.y, xyzF.z, NOrig) - proxQRead(gradX, xyzF.x, xyzF.y, xyzF.z, NOrig));
	apuVal += (proxQRead(gradY, xyzF.x, xyzF.y + 1, xyzF.z, NOrig) - proxQRead(gradY, xyzF.x, xyzF.y, xyzF.z, NOrig));
	apuVal += (proxQRead(gradZ, xyzF.x, xyzF.y, xyzF.z + 1, NOrig) - proxQRead(gradZ, xyzF.x, xyzF.y, xyzF.z, NOrig));
#else
#endif // END DIFFTYPE
	output[x] -= apuVal;
}

// Proximal TV or TGV gradient computation, i.e. the dual update with the projection
KERNEL3
void ProxTVGradient(const int3 N, const int3 NOrig, const CLGLOBAL float* CLRESTRICT im, CLGLOBAL float* outputX, CLGLOBAL float* outputY, CLGLOBAL float* outputZ, const float sigma2, const float alpha
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
	float imApu = im[x];
#if defined(EFOVZ)
	const LTYPE3 NDiff = (N - NOrig) / 2;
	const LTYPE y = (xyz.x - NDiff.x) + (xyz.y - NDiff.y) * NOrig.x + (xyz.z - NDiff.z) * NOrig.x * NOrig.y;
#else
	const LTYPE y = x;
#endif
	// Differences of the image
	float dX = FLOAT_ZERO;
	float dY = FLOAT_ZERO;
	float dZ = FLOAT_ZERO;
#if DIFFTYPE == 0 // START DIFFTYPE == 0
	// Forward difference
	forwardDiffX2(&dX, xyz, N, im, imApu);
	forwardDiffY2(&dY, xyz, N, im, imApu);
	forwardDiffZ2(&dZ, xyz, N, im, imApu);
#elif DIFFTYPE == 1 // START DIFFTYPE == 1
	// Backward difference
	backwardDiffX2(&dX, xyz, N, im, imApu);
	backwardDiffY2(&dY, xyz, N, im, imApu);
	backwardDiffZ2(&dZ, xyz, N, im, imApu);
#else
#endif // END DIFFTYPE
#ifdef PROXTGV
	dX -= vX[y];
	dY -= vY[y];
#ifdef TGVZ
	dZ -= vZ[y];
#endif
#endif
	// Dual update
	const float3 apu = MFLOAT3(outputX[y] + dX * sigma2, outputY[y] + dY * sigma2, outputZ[y] + dZ * sigma2);
	// Projection to the L2 or L1 ball (see http://dx.doi.org/10.1088/0031-9155/57/10/3065)
	// Voxels skipped above keep their previous, already projected, values
#ifdef L2 // START L2
// L2 norm
	const float scale = fmax(FLOAT_ONE, length(apu) / alpha);
	outputX[y] = apu.x / scale;
	outputY[y] = apu.y / scale;
	outputZ[y] = apu.z / scale;
#else
// L1 norm (anisotropic): the projection is a separate clamp of each component to [-alpha, alpha]
	outputX[y] = apu.x / fmax(FLOAT_ONE, fabs(apu.x) / alpha);
	outputY[y] = apu.y / fmax(FLOAT_ONE, fabs(apu.y) / alpha);
	outputZ[y] = apu.z / fmax(FLOAT_ONE, fabs(apu.z) / alpha);
#endif // END L2
}
#endif // END PROXTV

#ifdef PROXTGV // START PROXTGV
KERNEL3
// Symmetric derivative for TGV, i.e. the dual update with the projection
void ProxTGVSymmDeriv(const int3 N, const int3 NOrig, const CLGLOBAL float* CLRESTRICT vX, const CLGLOBAL float* CLRESTRICT vY, 
#ifdef TGVZ
	const CLGLOBAL float* CLRESTRICT vZ, CLGLOBAL float* qX, CLGLOBAL float* qY, CLGLOBAL float* qZ, CLGLOBAL float* q2XY, CLGLOBAL float* q2XZ, CLGLOBAL float* q2YZ, 
#else
	CLGLOBAL float* qX, CLGLOBAL float* qY, CLGLOBAL float* q2XY, 
#endif
	const float sigma2, const float alpha
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
	const float imApuX = vX[x];
	const float imApuY = vY[x];
	float dXX = FLOAT_ZERO;
	float dYY = FLOAT_ZERO;
	float dXY = FLOAT_ZERO;
#ifdef TGVZ
	const float imApuZ = vZ[x];
	float dZZ = FLOAT_ZERO;
	float dXZ = FLOAT_ZERO;
	float dYZ = FLOAT_ZERO;
#endif
#if DIFFTYPE == 0 // START DIFFTYPE == 0
	// Forward difference
	forwardDiffX2(&dXX, xyz, NOrig, vX, imApuX);
	forwardDiffY2(&dYY, xyz, NOrig, vY, imApuY);
	// XY/YX
	forwardDiffY2(&dXY, xyz, NOrig, vX, imApuX);
	forwardDiffX2(&dXY, xyz, NOrig, vY, imApuY);
#ifdef TGVZ
	forwardDiffZ2(&dZZ, xyz, NOrig, vZ, imApuZ);
	// XZ/ZX
	forwardDiffZ2(&dXZ, xyz, NOrig, vX, imApuX);
	forwardDiffX2(&dXZ, xyz, NOrig, vZ, imApuZ);
	// YZ/ZY
	forwardDiffZ2(&dYZ, xyz, NOrig, vY, imApuY);
	forwardDiffY2(&dYZ, xyz, NOrig, vZ, imApuZ);
#endif
#elif DIFFTYPE == 1 // START DIFFTYPE == 1
	// Backward difference
	backwardDiffX2(&dXX, xyz, NOrig, vX, imApuX);
	backwardDiffY2(&dYY, xyz, NOrig, vY, imApuY);
	// XY/YX
	backwardDiffY2(&dXY, xyz, NOrig, vX, imApuX);
	backwardDiffX2(&dXY, xyz, NOrig, vY, imApuY);
#ifdef TGVZ
	backwardDiffZ2(&dZZ, xyz, NOrig, vZ, imApuZ);
	// XZ/ZX
	backwardDiffZ2(&dXZ, xyz, NOrig, vX, imApuX);
	backwardDiffX2(&dXZ, xyz, NOrig, vZ, imApuZ);
	// YZ/ZY
	backwardDiffZ2(&dYZ, xyz, NOrig, vY, imApuY);
	backwardDiffY2(&dYZ, xyz, NOrig, vZ, imApuZ);
#endif
#else
#endif // END DIFFTYPE
	// Dual update
	const float sigmaHalf = sigma2 * FLOAT_HALF;
#ifdef TGVZ
	const float3 apu = MFLOAT3(qX[x] + dXX * sigma2, qY[x] + dYY * sigma2, qZ[x] + dZZ * sigma2);
	const float3 apu2 = MFLOAT3(q2XY[x] + dXY * sigmaHalf, q2XZ[x] + dXZ * sigmaHalf, q2YZ[x] + dYZ * sigmaHalf);
#else
	const float3 apu = MFLOAT3(qX[x] + dXX * sigma2, qY[x] + dYY * sigma2, FLOAT_ZERO);
	const float3 apu2 = MFLOAT3(q2XY[x] + dXY * sigmaHalf, FLOAT_ZERO, FLOAT_ZERO);
#endif
	// Projection to the L2 or L1 ball
	// Voxels skipped above keep their previous, already projected, values
#ifdef L2 // START L2
	const float scale = fmax(FLOAT_ONE, SQRT(apu.x * apu.x + apu.y * apu.y + apu.z * apu.z + (apu2.x * apu2.x) * FLOAT_TWO + (apu2.y * apu2.y) * FLOAT_TWO + (apu2.z * apu2.z) * FLOAT_TWO) / alpha);
	qX[x] = apu.x / scale;
	qY[x] = apu.y / scale;
	q2XY[x] = apu2.x / scale;
#ifdef TGVZ
	qZ[x] = apu.z / scale;
	q2XZ[x] = apu2.y / scale;
	q2YZ[x] = apu2.z / scale;
#endif
#else
// L1 norm (anisotropic): the projection is a separate clamp of each component to [-alpha, alpha]
	qX[x] = apu.x / fmax(FLOAT_ONE, fabs(apu.x) / alpha);
	qY[x] = apu.y / fmax(FLOAT_ONE, fabs(apu.y) / alpha);
	q2XY[x] = apu2.x / fmax(FLOAT_ONE, fabs(apu2.x) / alpha);
#ifdef TGVZ
	qZ[x] = apu.z / fmax(FLOAT_ONE, fabs(apu.z) / alpha);
	q2XZ[x] = apu2.y / fmax(FLOAT_ONE, fabs(apu2.y) / alpha);
	q2YZ[x] = apu2.z / fmax(FLOAT_ONE, fabs(apu2.z) / alpha);
#endif
#endif // END L2
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
				lCache[indX + indY * SW_TILEX + indZ * SW_TILEX * SW_TILEY] = SWBufRead(u, (int)xx, (int)yy, (int)zz, N);
#endif
				indX += LSIZE0;
			}
			indY += LSIZE1;
		}
		indZ += LSIZE2;
	}
	BARRIER
	// The bounds/EFOV/mask early returns must come after the cooperative tile load and barrier,
	// otherwise edge work-groups leave tile entries uninitialised
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
