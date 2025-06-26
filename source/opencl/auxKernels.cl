
/*******************************************************************************************************************************************
* This file contains "auxliary" kernels. What this means is that this file contains kernels for many of the regularization techniques/priors
* as well as kernels for some algorithm computations. The latter contains also convolution, element-wise multiplication and division, 
* derivatives and other functions. This file uses the preprocessor definitions and functions from general_opencl_functions.h. Note that
* the inclusion is not done here but rather during the compilation.
*
* Copyright (C) 2019-2024 Ville-Veikko Wettenhovi
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
#if defined(NLM_) || defined(PROXNLM)
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
void forward(CLGLOBAL FLOAT* d_outputFP, const CLGLOBAL FLOAT* CLRESTRICT meas
#ifdef CT
	, CLGLOBAL FLOAT* d_outputCT
#endif
#ifdef RANDOMS
	, const CLGLOBAL FLOAT* CLRESTRICT rand
#endif
	) {
	uint gid = GID0;
#ifdef CT
	const FLOAT apu = EXP(-d_outputFP[gid]);
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
void computeEstimate(const CLGLOBAL CAST* CLRESTRICT d_Summ, const CLGLOBAL CAST* CLRESTRICT d_rhs, CLGLOBAL FLOAT* d_im, const FLOAT d_epps, const int3 d_N, const uchar no_norm
#ifdef CT
	, const FLOAT flat
#endif
	) {

	int3 ind = CMINT3(GID0, GID1, GID2);
	size_t idx = GID0 + GID1 * GSIZE0 + GID2 * GSIZE0 * GSIZE1;
#ifdef CUDA
	if (ind.x >= d_N.x || ind.y >= d_N.y || ind.z >= d_N.z)
#else
	if (any(ind >= d_N))
#endif
		return;
	FLOAT apu = d_im[idx];
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
// This is mainly for non-FLOAT inputs
#ifdef PSF // START PSF
KERN
void Convolution3D(const CLGLOBAL CAST* input, CLGLOBAL CAST* output,
	CONSTANT FLOAT* convolution_window, int window_size_x, int window_size_y, int window_size_z) {
	int4 ind = CMINT4(GID0, GID1, GID2, 0);
	int4 ind_uus = CMINT4(0, 0, 0, 0);
	const uint Nyx = GSIZE0 * GSIZE1;
	FLOAT result = 0.f;
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
				FLOAT p = convert_float(input[indeksi]) / TH;
#else
				FLOAT p = input[indeksi];
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
void Convolution3D_f(const CLGLOBAL FLOAT* input, CLGLOBAL FLOAT* output,
	CONSTANT FLOAT* convolution_window, int window_size_x, int window_size_y, int window_size_z) {
	int4 ind = CMINT4(GID0, GID1, GID2, 0);
	int4 ind_uus = CMINT4(0, 0, 0, 0);
	FLOAT result = 0.f;
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
				FLOAT p = input[indeksi];
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
void vectorDiv(const CLGLOBAL FLOAT* input, CLGLOBAL FLOAT* output, const FLOAT epps) {
	uint id = GID0;
	output[id] = output[id] / (input[id] + epps);
}

// Elementwise multiplication
KERN
void vectorMult(const CLGLOBAL FLOAT* input, CLGLOBAL FLOAT* output) {
	uint id = GID0;
	output[id] *= input[id];
}
#endif // END PSF
#endif // END NOTAF

// Complex elementwise multiplication
// Used by the filtering
// This kernel assumes that the imaginary element is right after the real element, i.e. [real,imaginary,real,imaginary,...]
KERN
void vectorElementMultiply(const CLGLOBAL float* CLRESTRICT input, CLGLOBAL float* output, const uchar D2) {
	const LTYPE3 xyz = MINT3(GID0, GID1, GID2);
	const LTYPE n = xyz.x + xyz.y * GSIZE0 + xyz.z * GSIZE0 * GSIZE1;
	FLOAT mult;
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
	const FLOAT div = input[xyz.x];
	output[2 * n] /= div;
	output[2 * n + 1] /= div;
}

// Non-local means
#ifdef NLM_ // START NLM
#if defined(USEIMAGES) && defined(OPENCL)
CONSTANT sampler_t samplerNLM = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP;
#endif

KERNEL3
#ifdef USEIMAGES
void NLM(CLGLOBAL FLOAT* CLRESTRICT grad, IMAGE3D u, CONSTANT float* gaussian, 
#else
void NLM(CLGLOBAL FLOAT* CLRESTRICT grad, const CLGLOBAL FLOAT* CLRESTRICT u, CONSTANT float* gaussian, 
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
	, const CLGLOBAL FLOAT* CLRESTRICT u_ref
#endif
#endif // END NLMREF
#ifdef MASKPRIOR
#ifdef USEIMAGES
#ifdef MASKBP3D
	, IMAGE3D maskBP
#else
	, IMAGE2D maskBP
#endif
#else
	, const CLGLOBAL uchar* CLRESTRICT maskBP
#endif
#endif
#ifdef EFOVZ // Compute only in the voxels of the actual FOV (when using extended FOV)
	, CONSTANT uchar* fovIndices
#endif
) {
#ifdef PYTHON
	const int3 N = MINT3(Nx, Ny, Nz);
#endif
	LTYPE3 ii = MINT3(GID0, GID1, GID2);
	const LTYPE n = (ii.x) + (ii.y) * (N.x) + (ii.z) * (N.x * N.y);
	float weight_sum = (float)epps;
	float output = (float)0.f;
#if NLTYPE == 1
	float outputAla = epps;
#endif
	LTYPE startX = GRID0 * LSIZE0 - SWINDOWX - PWINDOWX + LID0;
	LTYPE startY = GRID1 * LSIZE1 - SWINDOWY - PWINDOWY + LID1;
	LTYPE startZ = GRID2 * LSIZE2 - SWINDOWZ - PWINDOWZ + LID2;
	LTYPE endX = (GRID0 + 1) * LSIZE0 + SWINDOWX + PWINDOWX;
	LTYPE endY = (GRID1 + 1) * LSIZE1 + SWINDOWY + PWINDOWY;
	LTYPE endZ = (GRID2 + 1) * LSIZE2 + SWINDOWZ + PWINDOWZ;
	LOCAL FLOAT lCache[SIZEX][SIZEY][SIZEZ];
#ifdef NLMREF
	LOCAL float lCacheRef[SIZEX][SIZEY][SIZEZ];
#endif
	LTYPE indZ = LID2;
	for (LTYPE zz = startZ; zz < endZ; zz += LSIZE2) {
		LTYPE indY = LID1;
		for (LTYPE yy = startY; yy < endY; yy += LSIZE1) {
			LTYPE indX = LID0;
			for (LTYPE xx = startX; xx < endX; xx += LSIZE0) {
#if defined(NLMREF) // START NLMREF
#ifdef USEIMAGES
#ifdef CUDA
				lCacheRef[indX][indY][indZ] = tex3D<float>(u_ref, xx, yy, zz);
#else
#ifdef HALF
				lCacheRef[indX][indY][indZ] = read_imageh(u_ref, samplerNLM, (int4)(xx, yy, zz, 0)).w;
#else
				lCacheRef[indX][indY][indZ] = read_imagef(u_ref, samplerNLM, (int4)(xx, yy, zz, 0)).w;
#endif
#endif
#else
				if (xx < 0 || yy < 0 || zz < 0 || xx >= N.x || yy >= N.y || zz >= N.z)
					lCacheRef[indX][indY][indZ] = (float)0.f;
				else
					lCacheRef[indX][indY][indZ] = u_ref[(xx) + (yy) * N.x + (zz) * N.x * N.y];
#endif
#endif // END NLMREF
#ifdef USEIMAGES
#ifdef CUDA
				lCache[indX][indY][indZ] = tex3D<float>(u, xx, yy, zz);
#else
#ifdef HALF
				lCache[indX][indY][indZ] = read_imageh(u, samplerNLM, (int4)(xx, yy, zz, 0)).w;
#else
				lCache[indX][indY][indZ] = read_imagef(u, samplerNLM, (int4)(xx, yy, zz, 0)).w;
#endif
#endif
#else
				if (xx < 0 || yy < 0 || zz < 0 || xx >= N.x || yy >= N.y || zz >= N.z)
					lCache[indX][indY][indZ] = (float)0.f;
				else
					lCache[indX][indY][indZ] = (float)u[(xx) + (yy) * N.x + (zz) * N.x * N.y];
#endif
				indX += LSIZE0;
			}
			indY += LSIZE1;
		}
		indZ += LSIZE2;
	}
	BARRIER
#ifdef CUDA
	if (ii.x >= N.x || ii.y >= N.y || ii.z >= N.z)
#else
	if (any(ii >= N))
#endif
		return;
#ifdef MASKPRIOR
#ifdef USEIMAGES
#ifdef CUDA
#ifdef MASKBP3D
    const int maskVal = tex3D<unsigned char>(maskBP, ii.x, ii.y, ii.z);
#else
    const int maskVal = tex2D<unsigned char>(maskBP, ii.x, ii.y);
#endif
#else
#ifdef MASKBP3D
    const int maskVal = read_imageui(maskBP, sampler_MASK, (int4)(ii.x, ii.y, ii.z, 0)).w;
#else
    const int maskVal = read_imageui(maskBP, sampler_MASK, (int2)(ii.x, ii.y)).w;
#endif
#endif
#else
	const int maskVal = maskBP[n];
#endif
#ifndef MASKSCALE
    if (maskVal == 0)
        return;
#endif
#endif
#if defined(NLMADAPTIVE)
	float hh = 0.f;
	const float pSize = CFLOAT((PWINDOWX * 2 + 1) * (PWINDOWY * 2 + 1) * (PWINDOWZ * 2 + 1));
#endif
	const int3 xxyyzz = CMINT3(LID0 + SWINDOWX + PWINDOWX, LID1 + SWINDOWY + PWINDOWY, LID2 + SWINDOWZ + PWINDOWZ);
	const float uj = (float)lCache[xxyyzz.x][xxyyzz.y][xxyyzz.z];
#ifdef MASKSCALE
	if (maskVal == 0) {
#pragma unroll
		for (int i = -1; i <= 1; i++) {
#pragma unroll
			for (int j = -1; j <= 1; j++) {
				int k = 0;
				if (i == 0 && j == 0)
					continue;
				float weight = 0.f;
				float distance = 0.f;
				int pz = 0;
#pragma unroll
					for (int py = -1; py <= 1; py++) {
						int dim_g = (pz + 1) * (3) * (3) + (py + 1) * (3);
#pragma unroll
						for (int px = -1; px <= 1; px++) {
							const float gg = (float)gaussian[dim_g++];
#ifdef NLMREF
							const float Pk = lCacheRef[xxyyzz.x + i + px][xxyyzz.y + j + py][xxyyzz.z];
							const float Pj = lCacheRef[xxyyzz.x + px][xxyyzz.y + py][xxyyzz.z];
#else
							const float Pk = lCache[xxyyzz.x + i + px][xxyyzz.y + j + py][xxyyzz.z];
							const float Pj = lCache[xxyyzz.x + px][xxyyzz.y + py][xxyyzz.z];
#endif
							const float PP = Pj - Pk;
							distance += gg * PP * PP;
						}
					}
#if defined(NLMADAPTIVE)
				hh = distance / pSize;
				weight = EXP(-distance / (float)(hh * h + s));
#else
 				weight = EXP(-distance / (float)h);
#endif
 				weight_sum += weight;
				const float uk = lCache[xxyyzz.x + i][xxyyzz.y + j][xxyyzz.z];
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
				output += weight * u * ((float)gamma * fabs(u) + uj + (float)3.f * uk + (float)(epps * epps)) / (divPow * divPow); 
#else
				const float divPow = FMAD((float)gamma, fabs(u), uj + uk + epps);
				output += weight * u * (FMAD((float)gamma, fabs(u), uj + (float)3.f * uk + (float)(epps * epps))) / (divPow * divPow); 
#endif // END FMAD
#elif NLTYPE == 4
				// Lange
				const float u = (uj - uk);
				const float uabs = sign(u);
				output += weight * (uabs - uabs / (fabs(u) / (float)gamma + (float)1.f));
#elif NLTYPE == 6
				// NLGGMRF
				const float delta = uj - uk;
				const float deltapqc = 1.f + POWR((fabs(delta / c)), (float)(p - q));
				output += weight * (float)((POWR(fabs(delta), p - 1.f) / deltapqc) * (p - gamma * (POWR(fabs(delta), p - q) / deltapqc))) * sign(delta);
#elif NLTYPE == 7
				const float u = (uk - uj);
				const float apu = (u * u + (float)(gamma * gamma));
// #ifndef USEMAD // START FMAD
				output += (((float)2.f * u * u * u) / (apu * apu) - (float)2.f * (u / apu));
// #else
// 				output += ((2.f * u * u * u) / FMAD(apu, apu, -2.f * (u / apu)));
// #endif // END FMAD
#else
 				//NLTV
				const float apuU = uj - uk;
 				output += (weight * apuU);
 				outputAla += weight * apuU * apuU;
#endif // END NLM NLTYPE
				}
			}
		}
	else {
#endif
		FLOAT pCache[PWINDOWX * 2 + 1][PWINDOWY * 2 + 1][PWINDOWZ * 2 + 1];
#pragma unroll
				for (int pz = -PWINDOWZ; pz <= PWINDOWZ; pz++) {
#pragma unroll
					for (int py = -PWINDOWY; py <= PWINDOWY; py++) {
#pragma unroll
						for (int px = -PWINDOWX; px <= PWINDOWX; px++) {
							pCache[px + PWINDOWX][py + PWINDOWY][pz + PWINDOWZ] = lCache[xxyyzz.x + px][xxyyzz.y + py][xxyyzz.z + pz];
						}
					}
				}
#if PWINDOWZ > 0
#pragma unroll
#endif
	for (int k = -SWINDOWZ; k <= SWINDOWZ; k++) {
#if PWINDOWZ > 0
#pragma unroll
#endif
		for (int j = -SWINDOWY; j <= SWINDOWY; j++) {
#if PWINDOWZ > 0
#pragma unroll
#endif
			for (int i = -SWINDOWX; i <= SWINDOWX; i++) {
				if (i == 0 && j == 0 && k == 0)
					continue;
				float weight = (float)0.f;
				float distance = (float)0.f;
				// nz, ny, nx = 0;
				// float uk = 0.f;
				// const float gg = EXP(-(Cfloat((i)*(i)) / (2.f * 25.f) + Cfloat((j)*(j)) / (2.f * 25.f) + Cfloat((k)*(k)) / (2.f * 25.f)));
#pragma unroll
				for (int pz = -PWINDOWZ; pz <= PWINDOWZ; pz++) {
#pragma unroll
					for (int py = -PWINDOWY; py <= PWINDOWY; py++) {
						int dim_g = (pz + PWINDOWZ) * (PWINDOWX * 2 + 1) * (PWINDOWY * 2 + 1) + (py + PWINDOWY) * (PWINDOWX * 2 + 1);
#pragma unroll
						for (int px = -PWINDOWX; px <= PWINDOWX; px++) {
							const float gg = (float)gaussian[dim_g++];
#ifdef NLMREF
							const float Pk = lCacheRef[xxyyzz.x + i + px][xxyyzz.y + j + py][xxyyzz.z + k + pz];
							const float Pj = lCacheRef[xxyyzz.x + px][xxyyzz.y + py][xxyyzz.z + pz];
#else
							// float Pk = 0.f;
							// if (i + px >= -PWINDOWX && i + px <= PWINDOWX && j + py >= -PWINDOWY && j + py <= PWINDOWY && k + pz >= -PWINDOWZ && k + pz <= PWINDOWZ)
							// 	Pk = pCache[px + i + PWINDOWX][py + j + PWINDOWY][pz + k + PWINDOWZ];
							// else
								// Pk = lCache[xxyyzz.x + i + px][xxyyzz.y + j + py][xxyyzz.z + k + pz];
							const float Pk = (float)lCache[xxyyzz.x + i + px][xxyyzz.y + j + py][xxyyzz.z + k + pz];
							// if (i == -SWINDOWX && j == -SWINDOWY && k == -SWINDOWZ)
								// testi[nx][ny][nz] = lCache[xxyyzz.x + px][xxyyzz.y + py][xxyyzz.z + pz];
							// testi[px][py][pz] = lCache[xxyyzz.x + px][xxyyzz.y + py][xxyyzz.z + pz];
							// const float Pj = lCache[xxyyzz.x + px][xxyyzz.y + py][xxyyzz.z + pz];
							const float Pj = (float)pCache[px + PWINDOWX][py + PWINDOWY][pz + PWINDOWZ];
							// if (pz == 0 && py == 0 && pz == 0)
							// 	uk = Pk;
#endif
							// const float PP = testi[nx][ny][nz] - Pk;
							const float PP = Pj - Pk;
							distance += gg * PP * PP;
							// distance += PP * PP;
							// nx++;
						}
						// ny++;
					}
					// nz++;
				}
#if defined(NLMADAPTIVE)
				hh = distance / pSize;
				weight = EXP(-distance / (float)(hh * h + s));
#else
 				weight = EXP(-distance / (float)h);
#endif
 				weight_sum += weight;
				const float uk = (float)lCache[xxyyzz.x + i][xxyyzz.y + j][xxyyzz.z + k];
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
				const float divPow = (uj + uk + (float)gamma * fabs(u) + (float)epps);
				output += weight * u * ((float)gamma * fabs(u) + uj + (float)3.f * uk + (float)(epps * epps)) / (divPow * divPow); 
#else
				const float divPow = FMAD((float)gamma, fabs(u), uj + uk + (float)epps);
				output += weight * u * (FMAD((float)gamma, fabs(u), uj + (float)3.f * uk + (float)(epps * epps))) / (divPow * divPow); 
#endif // END FMAD
#elif NLTYPE == 4
				// Lange
				const float u = (uj - uk);
				const float uabs = sign(u);
				output += weight * (uabs - uabs / (fabs(u) / (float)gamma + (float)1.f));
#elif NLTYPE == 6
				// NLGGMRF
				const float delta = uj - uk;
				const float deltapqc = 1.f + POWR(fabs(delta / c), p - q);
				output += weight * (float)((POWR(fabs(delta), p - 1.f) / deltapqc) * (p - gamma * (POWR(fabs(delta), p - q) / deltapqc))) * sign(delta);
#elif NLTYPE == 7
				const float u = (uk - uj);
				const float apu = (u * u + (float)(gamma * gamma));
// #ifndef USEMAD // START FMAD
				output += (((float)2.f * u * u * u) / (apu * apu) - (float)2.f * (u / apu));
// #else
// 				output += ((2.f * u * u * u) / FMAD(apu, apu, -2.f * (u / apu)));
// #endif // END FMAD
#else
 				//NLTV
				const float apuU = uj - uk;
 				output += (weight * apuU);
 				outputAla += weight * apuU * apuU;
#endif // END NLM NLTYPE
			}
		}
	}
#ifdef MASKSCALE
	}
#endif
	weight_sum = (float)1.f / weight_sum;
	output *= weight_sum;
#if NLTYPE == 2 // START NLM NLTYPE
	output = uj - output;
#elif NLTYPE == 5
	// Lange with NLMRP
	output = uj - output;
	const float uabs = sign(output);
	output = (uabs - uabs / (fabs(output) / (float)gamma + (float)1.f));
#elif NLTYPE == 1
#ifndef USEMAD // START FMAD
	output /= SQRT(outputAla * weight_sum + (float)epps);
#else
	output /= SQRT(FMAD(outputAla, weight_sum, (float)epps));
#endif // END FMAD
#endif // END NLM NLTYPE
	grad[n] += (FLOAT)(beta * output);
}
#endif // END NLM

// Relative difference prior
#ifdef RDP // START RDP
#ifdef OPENCL
#ifdef USEIMAGES
CONSTANT sampler_t samplerRDP = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP;
#endif
#ifdef RDPCORNERS
__kernel __attribute__((vec_type_hint(FLOAT))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, LOCAL_SIZE3)))
#else
__kernel __attribute__((vec_type_hint(float2))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, LOCAL_SIZE3)))
#endif
#else
extern "C" __global__
#endif
#ifdef USEIMAGES
void RDPKernel(CLGLOBAL FLOAT* CLRESTRICT grad, IMAGE3D u, 
#else
void RDPKernel(CLGLOBAL FLOAT* CLRESTRICT grad, const CLGLOBAL FLOAT* CLRESTRICT u, 
#endif
#ifdef PYTHON
	const int Nx, const int Ny, const int Nz, const int NOrigx, const int NOrigy, const int NOrigz, 
#else
	const int3 N, const int3 NOrig, 
#endif
	const float gamma, const float epps, const float beta
#ifdef MASKPRIOR
#ifdef MASKBP3D
	, IMAGE3D maskBP
#else
	, IMAGE2D maskBP
#endif
#endif
#ifdef EFOVZ
	, CONSTANT uchar* fovIndices
#endif
#ifdef RDPCORNERS
	, CONSTANT FLOAT* weight
#endif
#ifdef RDPREF
#ifdef USEIMAGES
	, IMAGE3D u_ref
#else
	, const CLGLOBAL FLOAT* CLRESTRICT u_ref
#endif
#endif
) {
#ifdef PYTHON
	const int3 N = MINT3(Nx, Ny, Nz);
#endif
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#ifdef RDPCORNERS // START RDPCORNERS
	FLOAT output = (FLOAT)0.f;
	LTYPE startX = GRID0 * LSIZE0 - SWINDOWX + LID0;
	LTYPE startY = GRID1 * LSIZE1 - SWINDOWY + LID1;
	LTYPE startZ = GRID2 * LSIZE2 - SWINDOWZ + LID2;
	LTYPE endX = (GRID0 + 1) * LSIZE0 + SWINDOWX;
	LTYPE endY = (GRID1 + 1) * LSIZE1 + SWINDOWY;
	LTYPE endZ = (GRID2 + 1) * LSIZE2 + SWINDOWZ;
	LOCAL FLOAT lCache[SIZEX][SIZEY][SIZEZ];
#ifdef RDPREF
	LOCAL FLOAT lCacheRef[SIZEX][SIZEY][SIZEZ];
#endif
	LTYPE indZ = LID2;
	for (LTYPE zz = startZ; zz < endZ; zz += LSIZE2) {
		LTYPE indY = LID1;
		for (LTYPE yy = startY; yy < endY; yy += LSIZE1) {
			LTYPE indX = LID0;
			for (LTYPE xx = startX; xx < endX; xx += LSIZE0) {
#if defined(RDPREF) // START RDPREF
#ifdef USEIMAGES
#ifdef CUDA
				lCacheRef[indX][indY][indZ] = tex3D<FLOAT>(u_ref, xx, yy, zz);
#else
				lCacheRef[indX][indY][indZ] = read_imagef(u_ref, samplerRDP, (int4)(xx, yy, zz, 0)).w;
#endif
#else
				if (xx < 0 || yy < 0 || zz < 0 || xx >= N.x || yy >= N.y || zz >= N.z)
					lCacheRef[indX][indY][indZ] = (FLOAT)0.f;
				else
					lCacheRef[indX][indY][indZ] = u_ref[(xx) + (yy) * N.x + (zz) * N.x * N.y];
#endif
#endif // END RDPREF
#ifdef USEIMAGES
#ifdef CUDA
				lCache[indX][indY][indZ] = tex3D<FLOAT>(u, xx, yy, zz);
#else
				lCache[indX][indY][indZ] = read_imagef(u, samplerRDP, (int4)(xx, yy, zz, 0)).w;
#endif
#else
				if (xx < 0 || yy < 0 || zz < 0 || xx >= N.x || yy >= N.y || zz >= N.z)
					lCache[indX][indY][indZ] = (FLOAT)0.f;
				else
					lCache[indX][indY][indZ] = u[(xx) + (yy) * N.x + (zz) * N.x * N.y];
#endif
				indX += LSIZE0;
			}
			indY += LSIZE1;
		}
		indZ += LSIZE2;
	}
	BARRIER
#endif // END RDPCORNERS
#ifdef CUDA
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
#ifdef CUDA
#ifdef MASKBP3D
    const int maskVal = tex3D<unsigned char>(maskBP, xyz.x, xyz.y, xyz.z);
#else
    const int maskVal = tex2D<unsigned char>(maskBP, xyz.x, xyz.y);
#endif
#else
#ifdef MASKBP3D
    const int maskVal = read_imageui(maskBP, sampler_MASK, (int4)(xyz.x, xyz.y, xyz.z, 0)).w;
#else
    const int maskVal = read_imageui(maskBP, sampler_MASK, (int2)(xyz.x, xyz.y)).w;
#endif
#endif
    if (maskVal == 0)
        return;
#endif
	const LTYPE n = (xyz.x) + (xyz.y) * (N.x) + (xyz.z) * (N.x * N.y);
// #endif
#ifdef RDPCORNERS // START RDPCORNERS
	const int3 xxyyzz = CMINT3(LID0 + SWINDOWX, LID1 + SWINDOWY, LID2 + SWINDOWZ);
#if defined(RDPREF) // START RDPREF
	const FLOAT kj = lCacheRef[xxyyzz.x][xxyyzz.y][xxyyzz.z];
#endif // END RDPREF
	const FLOAT uj = lCache[xxyyzz.x][xxyyzz.y][xxyyzz.z];
	int uu = 0;
	for (int i = -SWINDOWX; i <= SWINDOWX; i++) {
		for (int j = -SWINDOWY; j <= SWINDOWY; j++) {
			for (int k = -SWINDOWZ; k <= SWINDOWZ; k++) {
				if (i == 0 && j == 0 && k == 0)
					continue;
				const FLOAT uk = lCache[xxyyzz.x + i][xxyyzz.y + j][xxyyzz.z + k];
#if defined(RDPREF) // START RDPREF
				const FLOAT kk = lCacheRef[xxyyzz.x + i][xxyyzz.y + j][xxyyzz.z + k];
#endif // END RDPREF
				const FLOAT delta = uj - uk;
				const FLOAT divPow2 = FMAD((FLOAT)gamma, fabs(delta), uj + uk);
#if defined(RDPREF) // START RDPREF
				output += weight[uu++] * SQRT(kk * kj) * delta * ((FLOAT)gamma * fabs(delta) + uj + (FLOAT)3.f * uk + (FLOAT)(epps * epps)) / (divPow2 * divPow2 + (FLOAT)epps);
#else
				output += weight[uu++] * delta * ((FLOAT)gamma * fabs(delta) + uj + (FLOAT)3.f * uk + (FLOAT)(epps * epps)) / (divPow2 * divPow2 + (FLOAT)epps);
#endif // END RDPREF
			}
		}
	}
	if (isnan(output))
		output = (FLOAT)0.f;
	grad[n] += (FLOAT)(beta * (float)output);
#else
#ifdef USEIMAGES
#ifdef CUDA
	// Current voxel
	const FLOAT uj = tex3D<FLOAT>(u, xyz.x, xyz.y, xyz.z);
	// Left-right
	const float2 ux = make_float2(tex3D<FLOAT>(u, xyz.x + 1, xyz.y, xyz.z), tex3D<FLOAT>(u, xyz.x - 1, xyz.y, xyz.z));
	// Top-bottom
	const float2 uy = make_float2(tex3D<FLOAT>(u, xyz.x, xyz.y + 1, xyz.z), tex3D<FLOAT>(u, xyz.x, xyz.y - 1, xyz.z));
	// Front-back
	const float2 uz = make_float2(tex3D<FLOAT>(u, xyz.x, xyz.y, xyz.z + 1), tex3D<FLOAT>(u, xyz.x, xyz.y, xyz.z - 1));
#else
	// Current voxel
	const FLOAT uj = read_imagef(u, samplerRDP, (int4)(xyz.x, xyz.y, xyz.z, 0)).w;
	// Left-right
	const float2 ux = { read_imagef(u, samplerRDP, (int4)(xyz.x + 1, xyz.y, xyz.z, 0)).w, read_imagef(u, samplerRDP, (int4)(xyz.x - 1, xyz.y, xyz.z, 0)).w };
	// Top-bottom
	const float2 uy = { read_imagef(u, samplerRDP, (int4)(xyz.x, xyz.y + 1, xyz.z, 0)).w, read_imagef(u, samplerRDP, (int4)(xyz.x, xyz.y - 1, xyz.z, 0)).w };
	// Front-back
	const float2 uz = { read_imagef(u, samplerRDP, (int4)(xyz.x, xyz.y, xyz.z + 1, 0)).w, read_imagef(u, samplerRDP, (int4)(xyz.x, xyz.y, xyz.z - 1, 0)).w };
#endif
#else
	// Current voxel
	const FLOAT uj = u[xyz.x + xyz.y * N.x + xyz.z * N.x * N.y];
	float2 ux = MFLOAT2((FLOAT)0.f, (FLOAT)0.f);
	float2 uy = MFLOAT2((FLOAT)0.f, (FLOAT)0.f);
	float2 uz = MFLOAT2((FLOAT)0.f, (FLOAT)0.f);
	// Left-right
	if (xyz.x < N.x - 1)
		ux.x = u[(xyz.x + 1) + xyz.y * N.x + xyz.z * N.x * N.y];
	if (xyz.x > 0)
		ux.y = u[(xyz.x - 1) + xyz.y * N.x + xyz.z * N.x * N.y];
	// Top-bottom
	if (xyz.y < N.y - 1)
		uy.x = u[(xyz.x) + (xyz.y + 1) * N.x + (xyz.z) * N.x * N.y];
	if (xyz.y > 0)
		uy.y = u[(xyz.x) + (xyz.y - 1) * N.x + (xyz.z) * N.x * N.y];
	// Front-back
	if (xyz.z < N.z - 1)
		uz.x = u[(xyz.x) + (xyz.y) * N.x + (xyz.z + 1) * N.x * N.y];
	if (xyz.z > 0)
		uz.y = u[(xyz.x) + (xyz.y) * N.x + (xyz.z - 1) * N.x * N.y];
#endif
	const float2 uj_ux = uj - ux;
	const float2 uj_uy = uj - uy;
	const float2 uj_uz = uj - uz;
#ifndef USEMAD // START FMAD
	const float2 divPow2X = (uj + ux + gamma * fabs(uj_ux));
	const float2 divPow2Y = (uj + uy + gamma * fabs(uj_uy));
	const float2 divPow2Z = (uj + uz + gamma * fabs(uj_uz));
	float2 output = uj_ux * (gamma * fabs(uj_ux) + uj + 3.f * ux + epps * epps) / (divPow2X * divPow2X + epps) 
		+ uj_uy * (gamma * fabs(uj_uy) + uj + 3.f * uy + epps * epps) / (divPow2Y * divPow2Y + epps)
		+ uj_uz * (gamma * fabs(uj_uz) + uj + 3.f * uz + epps * epps) / (divPow2Z * divPow2Z + epps);
#else
	const float2 divPow2X = FMAD2(gamma, fabs(uj_ux), uj + ux);
	const float2 divPow2Y = FMAD2(gamma, fabs(uj_uy), uj + uy);
	const float2 divPow2Z = FMAD2(gamma, fabs(uj_uz), uj + uz);
	float2 output = uj_ux * FMAD2(gamma, fabs(uj_ux), uj + 3.f * ux + epps * epps) / (divPow2X * divPow2X + epps) 
		+ uj_uy * FMAD2(gamma, fabs(uj_uy), uj + 3.f * uy + epps * epps) / (divPow2Y * divPow2Y + epps)
		+ uj_uz * FMAD2(gamma, fabs(uj_uz), uj + 3.f * uz + epps * epps) / (divPow2Z * divPow2Z + epps);
#endif // END FMAD
	if (isnan(output.x))
		output.x = 0.f;
	if (isnan(output.y))
		output.y = 0.f;
	grad[n] += beta * (output.x + output.y);
#endif // END RDPCORNERS
}
#endif // END RDP


// Generalized Gaussian Markov random field
#ifdef GGMRF // START GGMRF
#if defined(USEIMAGES) && defined(OPENCL)
CONSTANT sampler_t samplerNLM = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP;
#endif

KERNEL3
#ifdef USEIMAGES
void GGMRFKernel(CLGLOBAL FLOAT* CLRESTRICT grad, IMAGE3D u, 
#else
void GGMRFKernel(CLGLOBAL FLOAT* CLRESTRICT grad, const CLGLOBAL FLOAT* CLRESTRICT u, 
#endif
	CONSTANT FLOAT* weight, const int3 N, const FLOAT p, const FLOAT q, const FLOAT c, const FLOAT pqc, const FLOAT beta
) {

	LTYPE3 ii = MINT3(GID0, GID1, GID2);
	const LTYPE n = (ii.x) + (ii.y) * (N.x) + (ii.z) * (N.x * N.y);
	FLOAT output = 0.f;
	LTYPE startX = GRID0 * LSIZE0 - SWINDOWX + LID0;
	LTYPE startY = GRID1 * LSIZE1 - SWINDOWY + LID1;
	LTYPE startZ = GRID2 * LSIZE2 - SWINDOWZ + LID2;
	LTYPE endX = (GRID0 + 1) * LSIZE0 + SWINDOWX;
	LTYPE endY = (GRID1 + 1) * LSIZE1 + SWINDOWY;
	LTYPE endZ = (GRID2 + 1) * LSIZE2 + SWINDOWZ;
	LOCAL FLOAT lCache[SIZEX][SIZEY][SIZEZ];
	LTYPE indZ = LID2;
	for (LTYPE zz = startZ; zz < endZ; zz += LSIZE2) {
		LTYPE indY = LID1;
		for (LTYPE yy = startY; yy < endY; yy += LSIZE1) {
			LTYPE indX = LID0;
			for (LTYPE xx = startX; xx < endX; xx += LSIZE0) {
#ifdef USEIMAGES
#ifdef CUDA
				lCache[indX][indY][indZ] = tex3D<FLOAT>(u, xx, yy, zz);
#else
				lCache[indX][indY][indZ] = read_imagef(u, samplerNLM, (int4)(xx, yy, zz, 0)).w;
#endif
#else
				if (xx < 0 || yy < 0 || zz < 0 || xx >= N.x || yy >= N.y || zz >= N.z)
					lCache[indX][indY][indZ] = 0.f;
				else
					lCache[indX][indY][indZ] = u[(xx) + (yy) * N.x + (zz) * N.x * N.y];
#endif
				indX += LSIZE0;
			}
			indY += LSIZE1;
		}
		indZ += LSIZE2;
	}
	BARRIER
#ifdef CUDA
	if (ii.x >= N.x || ii.y >= N.y || ii.z >= N.z)
#else
	if (any(ii >= N))
#endif
		return;
	const int3 xxyyzz = CMINT3(LID0 + SWINDOWX, LID1 + SWINDOWY, LID2 + SWINDOWZ);
	const FLOAT uj = lCache[xxyyzz.x][xxyyzz.y][xxyyzz.z];
	int uu = 0;
	for (int i = -SWINDOWX; i <= SWINDOWX; i++) {
		for (int j = -SWINDOWY; j <= SWINDOWY; j++) {
			for (int k = -SWINDOWZ; k <= SWINDOWZ; k++) {
				if (i == 0 && j == 0 && k == 0)
					continue;
				const FLOAT uk = lCache[xxyyzz.x + i][xxyyzz.y + j][xxyyzz.z + k];
				const FLOAT delta = uj - uk;
				const FLOAT deltapqc = 1.f + POWR(fabs(delta / c), p - q);
				if (delta != 0.f)
					output += weight[uu] * (POWR(fabs(delta), p - 1.f) / deltapqc) * (p - pqc * (POWR(fabs(delta), p - q) / deltapqc)) * sign(delta);
				uu++;
			}
		}
	}
	grad[n] += beta * output;
}
#endif

// Median root prior
#ifdef MEDIAN // START MEDIAN
KERNEL3
void medianFilter3D(const CLGLOBAL FLOAT* grad, CLGLOBAL FLOAT* output, const int3 N, const int3 NOrig
#ifdef MASKPRIOR
#ifdef MASKBP3D
	, IMAGE3D maskBP
#else
	, IMAGE2D maskBP
#endif
#endif
#ifdef EFOVZ
	, CONSTANT uchar* fovIndices
#endif
) {
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
	if (xyz.x >= N.x + SEARCH_WINDOW_X || xyz.y >= N.y + SEARCH_WINDOW_Y || xyz.z >= N.z + SEARCH_WINDOW_Z || xyz.x < SEARCH_WINDOW_X || xyz.y < SEARCH_WINDOW_Y || xyz.z < SEARCH_WINDOW_Z)
		return;
#ifdef EFOVZ
	if (fovIndices[xyz.z] == 0)
        return;
#endif
#ifdef MASKPRIOR
#ifdef CUDA
    const int maskVal = tex2D<unsigned char>(maskBP, xyz.x, xyz.y);
#else
#ifdef MASKBP3D
    const int maskVal = read_imageui(maskBP, sampler_MASK, (int4)(xyz.x, xyz.y, xyz.z, 0)).w;
#else
    const int maskVal = read_imageui(maskBP, sampler_MASK, (int2)(xyz.x, xyz.y)).w;
#endif
#endif
    if (maskVal == 0)
        return;
#endif
	const LTYPE n = (xyz.x - SEARCH_WINDOW_X) + (xyz.y - SEARCH_WINDOW_Y) * (N.x) + (xyz.z - SEARCH_WINDOW_Z) * (N.x * N.y);
	FLOAT median[KOKO];
	FLOAT medianF[KOKO];
	for (int ll = 0; ll < KOKO; ll++) {
		medianF[ll] = 0.f;
		median[ll] = 0.f;
	}
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
		medianF[ind] = median[hh];
		if (ind == KOKO / 2)
			break;
	}
	output[n] = medianF[KOKO / 2];
}
#endif // END MEDIAN

#if defined(PROXTV) || defined(TVGRAD) // START PROXTV || TVGRAD
// Backward difference for X-axis
DEVICE void backwardDiffX(FLOAT* apuVal, const LTYPE3 xyz, const int3 N, const LTYPE xx, const CLGLOBAL FLOAT* im) {
		const LTYPE xh = ((xyz.x - 1) + xyz.y * N.x + xyz.z * N.x * N.y);
		if (xyz.x == 0)
			*apuVal += im[xx];
		else
			*apuVal += (im[xx] - im[xh]);
}

// Backward difference for Y-axis
DEVICE void backwardDiffY(FLOAT* apuVal, const LTYPE3 xyz, const int3 N, const LTYPE xx, const CLGLOBAL FLOAT* im) {
		const LTYPE xh = ((xyz.x) + (xyz.y - 1) * N.x + xyz.z * N.x * N.y);
		if (xyz.y == 0)
			*apuVal += im[xx];
		else
			*apuVal += (im[xx] - im[xh]);
}

// Backward difference for Z-axis
DEVICE void backwardDiffZ(FLOAT* apuVal, const LTYPE3 xyz, const int3 N, const LTYPE xx, const CLGLOBAL FLOAT* im) {
		const LTYPE xh = ((xyz.x) + (xyz.y) * N.x + (xyz.z - 1) * N.x * N.y);
		if (xyz.z == 0)
			*apuVal += im[xx];
		else
			*apuVal += (im[xx] - im[xh]);
}

// Backward difference for X-axis (not the current voxel)
DEVICE void backwardDiffX2(FLOAT* apuVal, const LTYPE3 xyz, const int3 N, const CLGLOBAL FLOAT* im, const FLOAT imApu) {
		const LTYPE xh = ((xyz.x - 1) + xyz.y * N.x + xyz.z * N.x * N.y);
		if (xyz.x == 0)
			*apuVal += imApu;
		else
			*apuVal += (imApu - im[xh]);
}

// Backward difference for Y-axis (not the current voxel)
DEVICE void backwardDiffY2(FLOAT* apuVal, const LTYPE3 xyz, const int3 N, const CLGLOBAL FLOAT* im, const FLOAT imApu) {
		const LTYPE xh = ((xyz.x) + (xyz.y - 1) * N.x + (xyz.z) * N.x * N.y);
		if (xyz.y == 0)
			*apuVal += imApu;
		else
			*apuVal += (imApu - im[xh]);
}

// Backward difference for Z-axis (not the current voxel)
DEVICE void backwardDiffZ2(FLOAT* apuVal, const LTYPE3 xyz, const int3 N, const CLGLOBAL FLOAT* im, const FLOAT imApu) {
		const LTYPE xh = ((xyz.x) + (xyz.y) * N.x + (xyz.z - 1) * N.x * N.y);
		if (xyz.z == 0)
			*apuVal += imApu;
		else
			*apuVal += (imApu - im[xh]);
}

// Forward difference for X-axis
DEVICE void forwardDiffX(FLOAT* apuVal, const LTYPE3 xyz, const int3 N, const LTYPE xx, const CLGLOBAL FLOAT* im) {
		const LTYPE xh = ((xyz.x + 1) + xyz.y * N.x + xyz.z * N.x * N.y);
		if (xyz.x == N.x - 1)
			*apuVal -= im[xx];
		else
			*apuVal += (im[xh] - im[xx]);
}

// Forward difference for Y-axis
DEVICE void forwardDiffY(FLOAT* apuVal, const LTYPE3 xyz, const int3 N, const LTYPE xx, const CLGLOBAL FLOAT* im) {
		const LTYPE xh = ((xyz.x) + (xyz.y + 1) * N.x + xyz.z * N.x * N.y);
		if (xyz.y == N.y - 1)
			*apuVal -= im[xx];
		else
			*apuVal += (im[xh] - im[xx]);
}

// Forward difference for Z-axis
DEVICE void forwardDiffZ(FLOAT* apuVal, const LTYPE3 xyz, const int3 N, const LTYPE xx, const CLGLOBAL FLOAT* im) {
		const LTYPE xh = ((xyz.x) + (xyz.y) * N.x + (xyz.z + 1) * N.x * N.y);
		if (xyz.z == N.z - 1)
			*apuVal -= im[xx];
		else
			*apuVal += (im[xh] - im[xx]);
}

// Forward difference for X-axis (not the current voxel)
DEVICE void forwardDiffX2(FLOAT* apuVal, const LTYPE3 xyz, const int3 N, const CLGLOBAL FLOAT* im, const FLOAT imApu) {
		const LTYPE xh = ((xyz.x + 1) + xyz.y * N.x + xyz.z * N.x * N.y);
		if (xyz.x == N.x - 1)
			*apuVal -= imApu;
		else
			*apuVal += (im[xh] - imApu);
}

// Forward difference for Y-axis (not the current voxel)
DEVICE void forwardDiffY2(FLOAT* apuVal, const LTYPE3 xyz, const int3 N, const CLGLOBAL FLOAT* im, const FLOAT imApu) {
		const LTYPE xh = ((xyz.x) + (xyz.y + 1) * N.x + (xyz.z) * N.x * N.y);
		if (xyz.y == N.y - 1)
			*apuVal -= imApu;
		else
			*apuVal += (im[xh] - imApu);
}

// Forward difference for Z-axis (not the current voxel)
DEVICE void forwardDiffZ2(FLOAT* apuVal, const LTYPE3 xyz, const int3 N, const CLGLOBAL FLOAT* im, const FLOAT imApu) {
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
void ProxTVq(CLGLOBAL FLOAT* inputX, CLGLOBAL FLOAT* inputY, CLGLOBAL FLOAT* inputZ, const FLOAT alpha) {
	LTYPE idx = GID0;
	const float3 apu = MFLOAT3(inputX[idx], inputY[idx], inputZ[idx]);
#ifdef L2 // START L2
// L2 norm
	const FLOAT scale = fmax(1.f, length(apu) / alpha);
#else
// L1 norm
	const FLOAT scale = fmax(fmax(fabs(apu.z), fmax(fabs(apu.x), fabs(apu.y))) / alpha, 1.f);
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
void ProxTGVq(CLGLOBAL FLOAT* inputX, CLGLOBAL FLOAT* inputY, CLGLOBAL FLOAT* inputZ, CLGLOBAL FLOAT* input2XY, CLGLOBAL FLOAT* input2XZ, CLGLOBAL FLOAT* input2YZ, const FLOAT alpha) {
#else
void ProxTGVq(CLGLOBAL FLOAT* inputX, CLGLOBAL FLOAT* inputY, CLGLOBAL FLOAT* input2XY, const FLOAT alpha) {
#endif
	LTYPE idx = GID0;
#ifdef TGVZ
	const float3 apu = MFLOAT3(inputX[idx], inputY[idx], inputZ[idx]);
	const float3 apu2 = MFLOAT3(input2XY[idx], input2XZ[idx], input2YZ[idx]);
#else
	const float3 apu = MFLOAT3(inputX[idx], inputY[idx], 0.f);
	const float3 apu2 = MFLOAT3(input2XY[idx], 0.f, 0.f);
#endif
#ifdef L2 // START L2
	const FLOAT scale = fmax(1.f, SQRT(apu.x * apu.x + apu.y * apu.y + apu.z * apu.z + (apu2.x * apu2.x) * 2.f + (apu2.y * apu2.y) * 2.f + (apu2.z * apu2.z) * 2.f) / alpha);
#else
	const FLOAT scale = fmax(fmax(fabs(apu2.z),fmax(fabs(apu2.y), fmax(fabs(apu2.x), fmax(fabs(apu.z), fmax(fabs(apu.x), fabs(apu.y)))))) / alpha, 1.f);
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
void ProxTVDivergence(const int3 N, const int3 NOrig, const CLGLOBAL FLOAT* CLRESTRICT gradX, const CLGLOBAL FLOAT* CLRESTRICT gradY, const CLGLOBAL FLOAT* CLRESTRICT gradZ, CLGLOBAL FLOAT* output
// The (optional) logical mask should be zero in regions where the prior is not needed
#ifdef MASKPRIOR
#ifdef MASKBP3D
	, IMAGE3D maskBP
#else
	, IMAGE2D maskBP
#endif
#endif
// The (optional) logical vector should be zero in axial slices where the extended FOV is
#ifdef EFOVZ
	, CONSTANT uchar* fovIndices
#endif
) {
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#ifdef CUDA
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
#ifdef CUDA
    const int maskVal = tex2D<unsigned char>(maskBP, xyz.x, xyz.y);
#else
#ifdef MASKBP3D
    const int maskVal = read_imageui(maskBP, sampler_MASK, (int4)(xyz.x, xyz.y, xyz.z, 0)).w;
#else
    const int maskVal = read_imageui(maskBP, sampler_MASK, (int2)(xyz.x, xyz.y)).w;
#endif
#endif
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
	FLOAT apuVal = 0.f;
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
void ProxTVGradient(const int3 N, const int3 NOrig, const CLGLOBAL FLOAT* CLRESTRICT im, CLGLOBAL FLOAT* outputX, CLGLOBAL FLOAT* outputY, CLGLOBAL FLOAT* outputZ, const FLOAT sigma2
#ifdef PROXTGV
#ifdef TGVZ
	, const CLGLOBAL FLOAT* CLRESTRICT vX, const CLGLOBAL FLOAT* CLRESTRICT vY, const CLGLOBAL FLOAT* CLRESTRICT vZ
#else
	, const CLGLOBAL FLOAT* CLRESTRICT vX, const CLGLOBAL FLOAT* CLRESTRICT vY
#endif
#endif
#ifdef MASKPRIOR
#ifdef MASKBP3D
	, IMAGE3D maskBP
#else
	, IMAGE2D maskBP
#endif
#endif
#ifdef EFOVZ
	, CONSTANT uchar* fovIndices
#endif
) {
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#ifdef CUDA
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
#ifdef CUDA
    const int maskVal = tex2D<unsigned char>(maskBP, xyz.x, xyz.y);
#else
#ifdef MASKBP3D
    const int maskVal = read_imageui(maskBP, sampler_MASK, (int4)(xyz.x, xyz.y, xyz.z, 0)).w;
#else
    const int maskVal = read_imageui(maskBP, sampler_MASK, (int2)(xyz.x, xyz.y)).w;
#endif
#endif
    if (maskVal == 0)
        return;
#endif
		
	const LTYPE x = xyz.x + xyz.y * N.x + xyz.z * N.x * N.y;
	FLOAT apuVal = 0.f;
	FLOAT imApu = im[x];
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
	apuVal = 0.f;
	forwardDiffY2(&apuVal, xyz, N, im, imApu);
#ifdef PROXTGV
	apuVal -= vY[y];
#endif
	apuVal *=  sigma2;
	outputY[y] += apuVal;
	apuVal = 0.f;
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
	apuVal = 0.f;
	backwardDiffY2(&apuVal, xyz, N, im, imApu);
#ifdef PROXTGV
	apuVal -= vY[y];
#endif
	apuVal *=  sigma2;
	outputY[y] += apuVal;
	apuVal = 0.f;
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
void ProxTGVSymmDeriv(const int3 N, const int3 NOrig, const CLGLOBAL FLOAT* CLRESTRICT vX, const CLGLOBAL FLOAT* CLRESTRICT vY, 
#ifdef TGVZ
	const CLGLOBAL FLOAT* CLRESTRICT vZ, CLGLOBAL FLOAT* qX, CLGLOBAL FLOAT* qY, CLGLOBAL FLOAT* qZ, CLGLOBAL FLOAT* q2XY, CLGLOBAL FLOAT* q2XZ, CLGLOBAL FLOAT* q2YZ, 
#else
	CLGLOBAL FLOAT* qX, CLGLOBAL FLOAT* qY, CLGLOBAL FLOAT* q2XY, 
#endif
	const FLOAT sigma2
#ifdef MASKPRIOR
#ifdef MASKBP3D
	, IMAGE3D maskBP
#else
	, IMAGE2D maskBP
#endif
#endif
) {
	const LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#ifdef CUDA
	if (xyz.x >= NOrig.x || xyz.y >= NOrig.y || xyz.z >= NOrig.z)
#else
	if (any(xyz >= NOrig))
#endif
		return;
#ifdef MASKPRIOR
	const LTYPE3 NDiff = (N - NOrig) / 2;
#ifdef CUDA
    const int maskVal = tex2D<unsigned char>(maskBP, xyz.x + NDiff.x, xyz.y + NDiff.y);
#else
    const int maskVal = read_imageui(maskBP, sampler_MASK, (int2)(xyz.x + NDiff.x, xyz.y + NDiff.y)).w;
#endif
    if (maskVal == 0)
        return;
#endif
		
	const LTYPE x = xyz.x + xyz.y * NOrig.x + xyz.z * NOrig.x * NOrig.y;
/////////////// X ///////////////
	FLOAT apuVal = 0.f;
	FLOAT imApuX = vX[x];
// Forward difference
#if DIFFTYPE == 0 // START DIFFTYPE == 0
	forwardDiffX2(&apuVal, xyz, NOrig, vX, imApuX);
	apuVal *= sigma2;
	qX[x] += apuVal;
/////////////// Y ///////////////
	apuVal = 0.f;
	FLOAT imApuY = vY[x];
	forwardDiffY2(&apuVal, xyz, NOrig, vY, imApuY);
	apuVal *= sigma2;
	qY[x] += apuVal;
#ifdef TGVZ
/////////////// Z ///////////////
	apuVal = 0.f;
	FLOAT imApuZ = vZ[x];
	forwardDiffZ2(&apuVal, xyz, NOrig, vZ, imApuZ);
	apuVal *= sigma2;
	qZ[x] += apuVal;
#endif
/////////////// XY/YX ///////////////
	apuVal = 0.f;
	forwardDiffY2(&apuVal, xyz, NOrig, vX, imApuX);
	forwardDiffX2(&apuVal, xyz, NOrig, vY, imApuY);
	apuVal *= sigma2 * .5f;
	q2XY[x] += apuVal;
#ifdef TGVZ
/////////////// XZ/ZX ///////////////
	apuVal = 0.f;
	forwardDiffZ2(&apuVal, xyz, NOrig, vX, imApuX);
	forwardDiffX2(&apuVal, xyz, NOrig, vZ, imApuZ);
	apuVal *= sigma2 * .5f;
	q2XZ[x] += apuVal;
/////////////// YZ/ZY ///////////////
	apuVal = 0.f;
	forwardDiffZ2(&apuVal, xyz, NOrig, vY, imApuY);
	forwardDiffY2(&apuVal, xyz, NOrig, vZ, imApuZ);
	apuVal *= sigma2 * .5f;
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
	apuVal = 0.f;
	FLOAT imApuY = vY[x];
	backwardDiffY2(&apuValY, xyz, NOrig, vY, imApuY);
	apuValY *=  sigma2;
	qY[x] += apuValY;
#ifdef TGVZ
/////////////// Z ///////////////
	apuVal = 0.f;
	FLOAT imApuZ = vZ[x];
	backwardDiffZ2(&apuValZ, xyz, NOrig, vZ, imApuZ);
	apuValZ *=  sigma2;
	qZ[x] += apuValZ;
#endif
/////////////// XY/YX ///////////////
	apuVal = 0.f;
	backwardDiffY2(&apuVal, xyz, NOrig, vX, imApuX);
	backwardDiffX2(&apuVal, xyz, NOrig, vY, imApuY);
	apuVal *=  sigma2 * .5f;
	q2XY[x] += apuVal;
#ifdef TGVZ
/////////////// XZ/ZX ///////////////
	apuVal = 0.f;
	backwardDiffZ2(&apuVal, xyz, NOrig, vX, imApuX);
	backwardDiffX2(&apuVal, xyz, NOrig, vZ, imApuZ);
	apuVal *=  sigma2 * .5f;
	q2XZ[x] += apuVal;
/////////////// YZ/ZY ///////////////
	apuVal = 0.f;
	backwardDiffZ2(&apuVal, xyz, NOrig, vY, imApuY);
	backwardDiffY2(&apuVal, xyz, NOrig, vZ, imApuZ);
	apuVal *=  sigma2 * .5f;
	q2YZ[x] += apuVal;
#endif
// Central difference?
#else
#endif // END DIFFTYPE
}

// TGV divergence
KERNEL3
void ProxTGVDivergence(const int3 N, const int3 NOrig, const CLGLOBAL FLOAT* CLRESTRICT qX, const CLGLOBAL FLOAT* CLRESTRICT qY, 
#ifdef TGVZ
	const CLGLOBAL FLOAT* CLRESTRICT qZ, const CLGLOBAL FLOAT* CLRESTRICT q2XY, const CLGLOBAL FLOAT* CLRESTRICT q2XZ, const CLGLOBAL FLOAT* CLRESTRICT q2YZ, 
	CLGLOBAL FLOAT* vX, CLGLOBAL FLOAT* vY, CLGLOBAL FLOAT* vZ, 
#else
	const CLGLOBAL FLOAT* CLRESTRICT q2XY, CLGLOBAL FLOAT* vX, CLGLOBAL FLOAT* vY, 
#endif
	const CLGLOBAL FLOAT* CLRESTRICT pX, const CLGLOBAL FLOAT* CLRESTRICT pY, const CLGLOBAL FLOAT* CLRESTRICT pZ, const FLOAT theta, const FLOAT tau
#ifdef MASKPRIOR
#ifdef MASKBP3D
	, IMAGE3D maskBP
#else
	, IMAGE2D maskBP
#endif
#endif
) {
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#ifdef CUDA
	if (xyz.x >= NOrig.x || xyz.y >= NOrig.y || xyz.z >= NOrig.z)
#else
	if (any(xyz >= NOrig))
#endif
		return;
#ifdef MASKPRIOR
	const LTYPE3 NDiff = (N - NOrig) / 2;
#ifdef CUDA
    const int maskVal = tex2D<unsigned char>(maskBP, xyz.x + NDiff.x, xyz.y + NDiff.y);
#else
    const int maskVal = read_imageui(maskBP, sampler_MASK, (int2)(xyz.x + NDiff.x, xyz.y + NDiff.y)).w;
#endif
    if (maskVal == 0)
        return;
#endif
	const LTYPE x = xyz.x + xyz.y * NOrig.x + xyz.z * NOrig.x * NOrig.y;
	FLOAT apuVal = 0.f;
// Transpose of forward difference
#if DIFFTYPE == 0 // START DIFFTYPE == 0
/////////////// X ///////////////
		backwardDiffX(&apuVal, xyz, NOrig, x, qX);
		backwardDiffY(&apuVal, xyz, NOrig, x, q2XY);
#ifdef TGVZ
		backwardDiffZ(&apuVal, xyz, NOrig, x, q2XZ);
#endif
		FLOAT vApu = vX[x];
#ifdef USEMAD
		FLOAT v1 = FMAD(tau, pX[x] + apuVal, vApu);
		vX[x] = FMAD(theta, v1 - vApu, v1);
#else
		FLOAT v1 = vApu + tau * 1.f * (pX[x] + apuVal);
		vX[x] = v1 + theta * (v1 - vApu);
#endif
/////////////// Y ///////////////
		apuVal = 0.f;
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
		v1 = vApu + tau * 1.f * (pY[x] + apuVal);
		vY[x] = v1 + theta * (v1 - vApu);
#endif
#ifdef TGVZ
/////////////// Z ///////////////
		apuVal = 0.f;
		backwardDiffZ(&apuVal, xyz, NOrig, x, qZ);
		backwardDiffX(&apuVal, xyz, NOrig, x, q2XZ);
		backwardDiffY(&apuVal, xyz, NOrig, x, q2YZ);
		vApu = vZ[x];
#ifdef USEMAD
		v1 = FMAD(tau, pZ[x] + apuVal, vApu);
		vZ[x] = FMAD(theta, v1 - vApu, v1);
#else
		v1 = vApu + tau * 1.f * (pZ[x] + apuVal);
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
		FLOAT vApu = vX[x];
#ifdef USEMAD
		FLOAT v1 = FMAD(tau, pX[x] + apuVal, vApu);
		vX[x] = FMAD(theta, v1 - vApu, v1);
#else
		FLOAT v1 = vApu + tau * (pX[x] + apuVal);
		vX[x] = v1 + theta * (v1 - vApu);
#endif
/////////////// Y ///////////////
		apuVal = 0.f;
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
		apuVal = 0.f;
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
CONSTANT sampler_t samplerTV = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP;
#endif

#ifdef OPENCL
__kernel __attribute__((vec_type_hint(float2))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, LOCAL_SIZE3)))
#else
extern "C" __global__
#endif
#ifdef USEIMAGES
void hyperbolicKernel(CLGLOBAL FLOAT* CLRESTRICT grad, IMAGE3D u, 
#else
void hyperbolicKernel(CLGLOBAL FLOAT* CLRESTRICT grad, const CLGLOBAL FLOAT* CLRESTRICT u,
#endif
	const int3 N, const int3 NOrig, const FLOAT sigma, const FLOAT epps, const FLOAT beta, CONSTANT FLOAT* w
#ifdef MASKPRIOR
#ifdef MASKBP3D
	, IMAGE3D maskBP
#else
	, IMAGE2D maskBP
#endif
#endif
#ifdef EFOVZ
	, CONSTANT uchar* fovIndices
#endif
) {
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#ifdef CUDA
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
#ifdef CUDA
    const int maskVal = tex2D<unsigned char>(maskBP, xyz.x, xyz.y);
#else
#ifdef MASKBP3D
    const int maskVal = read_imageui(maskBP, sampler_MASK, (int4)(xyz.x, xyz.y, xyz.z, 0)).w;
#else
    const int maskVal = read_imageui(maskBP, sampler_MASK, (int2)(xyz.x, xyz.y)).w;
#endif
#endif
    if (maskVal == 0)
        return;
#endif
	const LTYPE n = (xyz.x) + (xyz.y) * (N.x) + (xyz.z) * (N.x * N.y);
	FLOAT output = 0.f;
	LTYPE startX = GRID0 * LSIZE0 - SWINDOWX + LID0;
	LTYPE startY = GRID1 * LSIZE1 - SWINDOWY + LID1;
	LTYPE startZ = GRID2 * LSIZE2 - SWINDOWZ + LID2;
	LTYPE endX = (GRID0 + 1) * LSIZE0 + SWINDOWX;
	LTYPE endY = (GRID1 + 1) * LSIZE1 + SWINDOWY;
	LTYPE endZ = (GRID2 + 1) * LSIZE2 + SWINDOWZ;
	LOCAL FLOAT lCache[SIZEX][SIZEY][SIZEZ];
	LTYPE indZ = LID2;
	for (LTYPE zz = startZ; zz < endZ; zz += LSIZE2) {
		LTYPE indY = LID1;
		for (LTYPE yy = startY; yy < endY; yy += LSIZE1) {
			LTYPE indX = LID0;
			for (LTYPE xx = startX; xx < endX; xx += LSIZE0) {
#ifdef USEIMAGES
#ifdef CUDA
				lCache[indX][indY][indZ] = tex3D<FLOAT>(u, xx, yy, zz);
#else
				lCache[indX][indY][indZ] = read_imagef(u, samplerTV, (int4)(xx, yy, zz, 0)).w;
#endif
#else
				if (xx < 0 || yy < 0 || zz < 0 || xx >= N.x || yy >= N.y || zz >= N.z)
					lCache[indX][indY][indZ] = 0.f;
				else
					lCache[indX][indY][indZ] = u[(xx) + (yy) * N.x + (zz) * N.x * N.y];
#endif
				indX += LSIZE0;
			}
			indY += LSIZE1;
		}
		indZ += LSIZE2;
	}
	BARRIER
	const int3 xxyyzz = CMINT3(LID0 + SWINDOWX, LID1 + SWINDOWY, LID2 + SWINDOWZ);
	const FLOAT uj = lCache[xxyyzz.x][xxyyzz.y][xxyyzz.z];
	int uu = 0;
	for (int k = -SWINDOWZ; k <= SWINDOWZ; k++) {
		for (int j = -SWINDOWY; j <= SWINDOWY; j++) {
			for (int i = -SWINDOWX; i <= SWINDOWX; i++) {
				if (i == 0 && j == 0 && k == 0)
					continue;
				const FLOAT u = lCache[xxyyzz.x + i][xxyyzz.y + j][xxyyzz.z + k];
				const FLOAT ux = (uj - u) / sigma;
				output += (ux / sigma) / SQRT(1.f + ux * ux) * w[uu++];
			}
		}
	}
	grad[n] += beta * output;
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
#ifdef OPENCL
CONSTANT sampler_t samplerTV = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP;
#endif

DEVICE FLOAT sqrtVal(const float3 input, const FLOAT epps
#ifdef TVW1
	, const float3 w
#endif
) {
#ifdef TVW1
#ifdef USEMAD
	return SQRT(FMAD(w.x, input.x * input.x, FMAD(w.y, input.y * input.y, FMAD(w.z, input.z * input.z, epps))));
#else
	return SQRT(w.x * input.x * input.x + w.y * input.y * input.y + w.z * input.z * input.z + epps);
#endif
#else
#ifdef USEMAD
	return SQRT(FMAD(input.x, input.x, FMAD(input.y, input.y, FMAD(input.z, input.z, epps))));
#else
	return SQRT(input.x * input.x + input.y * input.y + input.z * input.z + epps);
#endif
#endif
}

#ifdef OPENCL
__kernel __attribute__((vec_type_hint(float3))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, LOCAL_SIZE3)))
#else
extern "C" __global__
#endif
#ifdef USEIMAGES
void TVKernel(CLGLOBAL FLOAT* CLRESTRICT grad, IMAGE3D u, 
#else
void TVKernel(CLGLOBAL FLOAT* CLRESTRICT grad, const CLGLOBAL FLOAT* CLRESTRICT u,
#endif
#ifdef PYTHON
	const int Nx, const int Ny, const int Nz, const int NOrigx, const int NOrigy, const int NOrigz, 
#else
	const int3 N, const int3 NOrig, 
#endif
	const FLOAT sigma, const FLOAT epps, const FLOAT beta
#ifdef MASKPRIOR
#ifdef MASKBP3D
	, IMAGE3D maskBP
#else
	, IMAGE2D maskBP
#endif
#endif
#ifdef EFOVZ
	, CONSTANT uchar* fovIndices
#endif
#if defined(ANATOMICAL2) || defined(ANATOMICAL3)
	, const FLOAT C
#endif
#if defined(ANATOMICAL1) || defined(ANATOMICAL2) || defined(ANATOMICAL3)
	, CLGLOBAL FLOAT* CLRESTRICT S
#endif
) {
#ifdef PYTHON
	const int3 N = MINT3(Nx, Ny, Nz);
#endif
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#ifdef CUDA
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
#ifdef CUDA
    const int maskVal = tex2D<unsigned char>(maskBP, xyz.x, xyz.y);
#else
#ifdef MASKBP3D
    const int maskVal = read_imageui(maskBP, sampler_MASK, (int4)(xyz.x, xyz.y, xyz.z, 0)).w;
#else
    const int maskVal = read_imageui(maskBP, sampler_MASK, (int2)(xyz.x, xyz.y)).w;
#endif
#endif
    if (maskVal == 0)
        return;
#endif
	const LTYPE n = (xyz.x) + (xyz.y) * (N.x) + (xyz.z) * (N.x * N.y);
#ifdef USEIMAGES
#ifdef CUDA
	const FLOAT uijk = tex3D<FLOAT>(u, xyz.x, xyz.y, xyz.z);
#else
	const FLOAT uijk = read_imagef(u, samplerTV, (int4)(xyz.x, xyz.y, xyz.z, 0)).w;
#endif
#else
	const FLOAT uijk = u[n];
#endif
#if defined(SATV) // START JPTV || SATV
#ifdef USEIMAGES
#ifdef CUDA
	float2 ux = make_float2(tex3D<FLOAT>(u, xyz.x + 1, xyz.y, xyz.z), tex3D<FLOAT>(u, xyz.x - 1, xyz.y, xyz.z));
	float2 uy = make_float2(tex3D<FLOAT>(u, xyz.x, xyz.y + 1, xyz.z), tex3D<FLOAT>(u, xyz.x, xyz.y - 1, xyz.z));
	float2 uz = make_float2(tex3D<FLOAT>(u, xyz.x, xyz.y, xyz.z + 1), tex3D<FLOAT>(u, xyz.x, xyz.y, xyz.z - 1));
#else
	float2 ux = { read_imagef(u, samplerTV, (int4)(xyz.x + 1, xyz.y, xyz.z, 0)).w, read_imagef(u, samplerTV, (int4)(xyz.x - 1, xyz.y, xyz.z, 0)).w };
	float2 uy = { read_imagef(u, samplerTV, (int4)(xyz.x, xyz.y + 1, xyz.z, 0)).w, read_imagef(u, samplerTV, (int4)(xyz.x, xyz.y - 1, xyz.z, 0)).w };
	float2 uz = { read_imagef(u, samplerTV, (int4)(xyz.x, xyz.y, xyz.z + 1, 0)).w, read_imagef(u, samplerTV, (int4)(xyz.x, xyz.y, xyz.z - 1, 0)).w };
#endif
#else
	float2 ux = MFLOAT2(0.f, 0.f);
	float2 uy = MFLOAT2(0.f, 0.f);
	float2 uz = MFLOAT2(0.f, 0.f);
	if (xyz.x < N.x - 1)
		ux.x = u[(xyz.x + 1) + (xyz.y) * N.x + (xyz.z) * N.x * N.y];
	if (xyz.x > 0)
		ux.y = u[(xyz.x - 1) + (xyz.y) * N.x + (xyz.z) * N.x * N.y];
	if (xyz.y < N.y - 1)
		uy.x = u[(xyz.x) + (xyz.y + 1) * N.x + (xyz.z) * N.x * N.y];
	if (xyz.y > 0)
		uy.y = u[(xyz.x) + (xyz.y - 1) * N.x + (xyz.z) * N.x * N.y];
	if (xyz.z < N.z - 1)
		uz.x = u[(xyz.x) + (xyz.y) * N.x + (xyz.z + 1) * N.x * N.y];
	if (xyz.z > 0)
		uz.y = u[(xyz.x) + (xyz.y) * N.x + (xyz.z - 1) * N.x * N.y];
#endif
    ux = uijk - ux;
    uy = uijk - uy;
    uz = uijk - uz;
	const float2 uabsx = ux / (fabs(ux) + epps);
	const float2 uabsy = uy / (fabs(uy) + epps);
	const float2 uabsz = uz / (fabs(uz) + epps);
	float2 output = uabsx - uabsx / (fabs(ux) / sigma + 1.f) + uabsy - uabsy / (fabs(uy) / sigma + 1.f) + uabsz - uabsz / (fabs(uz) / sigma + 1.f);
	grad[n] += beta * (output.x + output.y);
#else
#ifdef USEIMAGES
#ifdef CUDA
	const float3 uijkP = make_float3(tex3D<FLOAT>(u, xyz.x + 1, xyz.y, xyz.z), tex3D<FLOAT>(u, xyz.x, xyz.y + 1, xyz.z), tex3D<FLOAT>(u, xyz.x, xyz.y, xyz.z + 1));
	const float3 uijkM = make_float3(tex3D<FLOAT>(u, xyz.x - 1, xyz.y, xyz.z), tex3D<FLOAT>(u, xyz.x, xyz.y - 1, xyz.z), tex3D<FLOAT>(u, xyz.x, xyz.y, xyz.z - 1));
	const float2 ui = make_float2(tex3D<FLOAT>(u, xyz.x - 1, xyz.y + 1, xyz.z), tex3D<FLOAT>(u, xyz.x - 1, xyz.y, xyz.z + 1));
	const float2 uj = make_float2(tex3D<FLOAT>(u, xyz.x + 1, xyz.y - 1, xyz.z), tex3D<FLOAT>(u, xyz.x, xyz.y - 1, xyz.z + 1));
	const float2 uk = make_float2(tex3D<FLOAT>(u, xyz.x + 1, xyz.y, xyz.z - 1), tex3D<FLOAT>(u, xyz.x, xyz.y + 1, xyz.z - 1));
#else
	const float3 uijkP = {read_imagef(u, samplerTV, (int4)(xyz.x + 1, xyz.y, xyz.z, 0)).w, read_imagef(u, samplerTV, (int4)(xyz.x, xyz.y + 1, xyz.z, 0)).w, read_imagef(u, samplerTV, (int4)(xyz.x, xyz.y, xyz.z + 1, 0)).w};
	const float3 uijkM = {read_imagef(u, samplerTV, (int4)(xyz.x - 1, xyz.y, xyz.z, 0)).w, read_imagef(u, samplerTV, (int4)(xyz.x, xyz.y - 1, xyz.z, 0)).w, read_imagef(u, samplerTV, (int4)(xyz.x, xyz.y, xyz.z - 1, 0)).w};
	const float2 ui = {read_imagef(u, samplerTV, (int4)(xyz.x - 1, xyz.y + 1, xyz.z, 0)).w, read_imagef(u, samplerTV, (int4)(xyz.x - 1, xyz.y, xyz.z + 1, 0)).w};
	const float2 uj = {read_imagef(u, samplerTV, (int4)(xyz.x + 1, xyz.y - 1, xyz.z, 0)).w, read_imagef(u, samplerTV, (int4)(xyz.x, xyz.y - 1, xyz.z + 1, 0)).w};
	const float2 uk = {read_imagef(u, samplerTV, (int4)(xyz.x + 1, xyz.y, xyz.z - 1, 0)).w, read_imagef(u, samplerTV, (int4)(xyz.x, xyz.y + 1, xyz.z - 1, 0)).w};
#endif
#else
	float3 uijkP = MFLOAT3(0.f, 0.f, 0.f);
	float3 uijkM = MFLOAT3(0.f, 0.f, 0.f);
	float2 ui = MFLOAT2(0.f, 0.f);
	float2 uj = MFLOAT2(0.f, 0.f);
	float2 uk = MFLOAT2(0.f, 0.f);
	if (xyz.x < N.x - 1)
		uijkP.x = u[(xyz.x + 1) + (xyz.y) * N.x + (xyz.z) * N.x * N.y];
	if (xyz.y < N.y - 1)
		uijkP.y = u[(xyz.x) + (xyz.y + 1) * N.x + (xyz.z) * N.x * N.y];
	if (xyz.z < N.z - 1)
		uijkP.z = u[(xyz.x) + (xyz.y) * N.x + (xyz.z + 1) * N.x * N.y];
	if (xyz.x > 0)
		uijkM.x = u[(xyz.x - 1) + (xyz.y) * N.x + (xyz.z) * N.x * N.y];
	if (xyz.y > 0)
		uijkM.y = u[(xyz.x) + (xyz.y - 1) * N.x + (xyz.z) * N.x * N.y];
	if (xyz.z > 0)
		uijkM.z = u[(xyz.x) + (xyz.y) * N.x + (xyz.z - 1) * N.x * N.y];
		
	if (xyz.x > 0 && xyz.y < N.y - 1)
		ui.x = u[(xyz.x - 1) + (xyz.y + 1) * N.x + (xyz.z) * N.x * N.y];
	if (xyz.x > 0 && xyz.z < N.z - 1)
		ui.y = u[(xyz.x - 1) + (xyz.y) * N.x + (xyz.z + 1) * N.x * N.y];
	if (xyz.y > 0 && xyz.x < N.x - 1)
		uj.x = u[(xyz.x + 1) + (xyz.y - 1) * N.x + (xyz.z) * N.x * N.y];
	if (xyz.y > 0 && xyz.z < N.z - 1)
		uj.y = u[(xyz.x) + (xyz.y - 1) * N.x + (xyz.z + 1) * N.x * N.y];
	if (xyz.z > 0 && xyz.x < N.x - 1)
		uk.x = u[(xyz.x + 1) + (xyz.y) * N.x + (xyz.z - 1) * N.x * N.y];
	if (xyz.z > 0 && xyz.y < N.y - 1)
		uk.y = u[(xyz.x) + (xyz.y + 1) * N.x + (xyz.z - 1) * N.x * N.y];
#endif
	const float3 u1 = MFLOAT3(uijk - uijkM.x, ui.x - uijkM.x, ui.y - uijkM.x);
	const float3 u2 = MFLOAT3(uj.x - uijkM.y, uijk - uijkM.y, uj.y - uijkM.y);
	const float3 u3 = MFLOAT3(uk.x - uijkM.z, uk.y - uijkM.z, uijk - uijkM.z);
#ifdef TVW1 // START TVW1
	const float3 u4 = uijkP - uijk;
	float3 w4 = (u4) / sigma;
	w4 = EXP3(-w4 * w4);
	const FLOAT pvalijk = sqrtVal(u4, epps, w4);
	float3 w1 = (u1) / sigma;
	w1 = EXP3(-w1 * w1);
	float3 w2 = (u2) / sigma;
	w2 = EXP3(-w2 * w2);
	float3 w3 = (u3) / sigma;
	w3 = EXP3(-w3 * w3);
#ifdef USEMAD
	grad[n] += beta * (-(FMAD(w4.x, u4.x, FMAD(w4.y, u4.y, w4.z * u4.z))) / pvalijk + (w1.x * (uijk - uijkM.x)) / sqrtVal(u1, epps, w1) + (w2.y * (uijk - uijkM.y)) / sqrtVal(u2, epps, w2) + (w3.y * (uijk - uijkM.z)) / sqrtVal(u3, epps, w3));
#else
	grad[n] += beta * (-(w4.x * u4.x + w4.y * u4.y + w4.z * u4.z) / pvalijk + (w1.x * (uijk - uijkM.x)) / sqrtVal(u1, epps, w1) + (w2.y * (uijk - uijkM.y)) / sqrtVal(u2, epps, w2) + (w3.y * (uijk - uijkM.z)) / sqrtVal(u3, epps, w3));
#endif
#else
#ifdef ANATOMICAL1 // TV type 1
	const LTYPE NN = N.x * N.y * N.z;
	__private FLOAT s[9];
	for (int kk = 0; kk < 9; kk++)
		s[kk] = S[n + NN * kk];
	const float3 val = uijkP - uijk;
	const FLOAT pvalijk = SQRT(val.x * val.x * s[0] + val.y * val.y * s[4] + val.z * val.z * s[8] + s[1] * (val.x) * (val.y) + s[3] * (val.x) * (val.y) + s[2] * (val.x) * (val.z) + s[6] * (val.x) * (val.z) + 
		s[5] * (val.y) * (val.z) + s[7] * (val.y) * (val.z) + epps);
	const FLOAT pvalijkX = SQRT(u1.x * u1.x * s[0] + u1.y * u1.y * s[4] + u1.z * u1.z * s[8] + s[1] * (u1.x) * (u1.y) + s[3] * (u1.x) * (u1.y) + s[2] * (u1.x) * (u1.z) + s[6] * (u1.x) * (u1.z) + 
		s[5] * (u1.y) * (u1.z) + s[7] * (u1.y) * (u1.z) + epps);
	const FLOAT pvalijkY = SQRT(u2.x * u2.x * s[0] + u2.y * u2.y * s[4] + u2.z * u2.z * s[8] + s[1] * (u2.x) * (u2.y) + s[3] * (u2.x) * (u2.y) + s[2] * (u2.x) * (u2.z) + s[6] * (u2.x) * (u2.z) + 
		s[5] * (u2.y) * (u2.z) + s[7] * (u2.y) * (u2.z) + epps);
	const FLOAT pvalijkZ = SQRT(u3.x * u3.x * s[0] + u3.y * u3.y * s[4] + u3.z * u3.z * s[8] + s[1] * (u3.x) * (u3.y) + s[3] * (u3.x) * (u3.y) + s[2] * (u3.x) * (u3.z) + s[6] * (u3.x) * (u3.z) + 
		s[5] * (u3.y) * (u3.z) + s[7] * (u3.y) * (u3.z) + epps);
	const FLOAT dx = s[0] * (2.f * (uijk - uijkM.x)) + s[3] * u1.y + s[2] * u1.z + s[6] * u1.z + s[1] * u1.y;
	const FLOAT dy = s[4] * (2.f * (uijk - uijkM.y)) + s[5] * u2.z + s[3] * u2.x + s[1] * u2.x + s[7] * u2.z;
	const FLOAT dz = s[8] * (2.f * (uijk - uijkM.z)) + s[6] * u3.x + s[5] * u3.y + s[7] * u3.y + s[2] * u3.x;
	const FLOAT d = s[1] * val.x + s[2] * val.x + s[3] * val.x + s[6] * val.x + s[1] * val.y + s[3] * val.y + s[5] * val.y + s[7] * val.y + s[2] * val.z + s[5] * val.z + s[6] * val.z + s[7] * val.z + s[0] * 2.f * val.x + s[4] * 2.f * val.y + s[8] * 2.f * val.z;
	grad[n] += beta * .5f * (d / pvalijk + dx / pvalijkX + dy / pvalijkY + dz / pvalijkZ);
#elif defined(ANATOMICAL2) // TV type 2
	float3 uijkR = MFLOAT3(0.f, 0.f, 0.f);
	if (xyz.x < N.x - 1)
		uijkR.x = S[(xyz.x + 1) + (xyz.y) * N.x + (xyz.z) * N.x * N.y];
	if (xyz.y < N.y - 1)
		uijkR.y = S[(xyz.x) + (xyz.y + 1) * N.x + (xyz.z) * N.x * N.y];
	if (xyz.z < N.z - 1)
		uijkR.z = S[(xyz.x) + (xyz.y) * N.x + (xyz.z + 1) * N.x * N.y];
	const float3 apuS = (uijkR - S[n]);
	const float3 apu = uijkP - uijk;
	const FLOAT pvalijk = SQRT(apu.x * apu.x + apu.y * apu.y + apu.z * apu.z + C * (apuS.x * apuS.x + apuS.y * apuS.y + apuS.z * apuS.z) + epps);
	grad[n] += beta * ((3.f * uijk - uijkP.x - uijkP.y - uijkP.z) / pvalijk + (uijk - uijkM.x) / sqrtVal(u1, epps) + (uijk - uijkM.y) / sqrtVal(u2, epps) + (uijk - uijkM.z) / sqrtVal(u3, epps) + 1e-7f);
#elif defined(ANATOMICAL3) // APLS
	float3 uijkR = MFLOAT3(0.f, 0.f, 0.f);
	if (xyz.x < N.x - 1)
		uijkR.x = S[(xyz.x + 1) + (xyz.y) * N.x + (xyz.z) * N.x * N.y];
	if (xyz.y < N.y - 1)
		uijkR.y = S[(xyz.x) + (xyz.y + 1) * N.x + (xyz.z) * N.x * N.y];
	if (xyz.z < N.z - 1)
		uijkR.z = S[(xyz.x) + (xyz.y) * N.x + (xyz.z + 1) * N.x * N.y];
	float3 epsilon = (uijkR - S[n]);
	epsilon = epsilon / SQRT(epsilon.x * epsilon.x + epsilon.y * epsilon.y + epsilon.z * epsilon.z + C * C);
	const float3 apu = uijkP - uijk;
	const FLOAT apuR = uijkR.x * apu.x + uijkR.y * apu.y + uijkR.z * apu.z;
	const FLOAT pvalijk = SQRT(apu.x * apu.x + apu.y * apu.y + apu.z * apu.z - apuR * apuR + epps);
	FLOAT apuRXYZ = uijkR.x * u1.x + uijkR.y * u1.y + uijkR.z * u1.z;
	const FLOAT pvalijkX = SQRT(u1.x * u1.x + u1.y * u1.y + u1.z * u1.z + apuRXYZ * apuRXYZ + epps);
	apuRXYZ = uijkR.x * u2.x + uijkR.y * u2.y + uijkR.z * u2.z;
	const FLOAT pvalijkY = SQRT(u2.x * u2.x + u2.y * u2.y + u2.z * u2.z + apuRXYZ * apuRXYZ + epps);
	apuRXYZ = uijkR.x * u3.x + uijkR.y * u3.y + uijkR.z * u3.z;
	const FLOAT pvalijkZ = SQRT(u3.x * u3.x + u3.y * u3.y + u3.z * u3.z + apuRXYZ * apuRXYZ + epps);
	grad[n] += beta * .5f * ((6.f * uijk - 2.f * uijkP.x - 2.f * uijkP.y - 2.f * uijkP.z + 2.f * (epsilon.x*(uijk - uijkP.x) + epsilon.y*(uijk - uijkP.y) + epsilon.z*(uijk - uijkP.z)) * (epsilon.x + epsilon.y + epsilon.z)) / pvalijk + 
		2.f * (u1.x - epsilon.x * (epsilon.x * u1.x + epsilon.y * u1.y + epsilon.z * u1.z)) / pvalijkX + 2.f * (u2.y - epsilon.y * (epsilon.x * u2.x + epsilon.y * u2.y + epsilon.z * u2.z)) / pvalijkY + 
		2.f * (u3.z - epsilon.z * (epsilon.x * u3.x + epsilon.y * u3.y + epsilon.z * u3.z))/ pvalijkZ + 1e-7f);
#else // Non-reference image TV
	const FLOAT pvalijk = sqrtVal(uijkP - uijk, epps);
	grad[n] += beta * ((3.f * uijk - uijkP.x - uijkP.y - uijkP.z) / pvalijk + (uijk - uijkM.x) / sqrtVal(u1, epps) + (uijk - uijkM.y) / sqrtVal(u2, epps) + (uijk - uijkM.z) / sqrtVal(u3, epps) + 1e-7f);
#endif
#endif // END TVW1
#endif // END JPTV || SATV
}
#endif // END TVGRAD

// This kernel computes the image estimate for PKMA, MBSREM and BSREM
// I.e. the step after backprojection has been computed
#if defined(PKMA) || defined(MBSREM) || defined(BSREM)
KERNEL3
void PoissonUpdate(CLGLOBAL FLOAT* CLRESTRICT im, const CLGLOBAL FLOAT* CLRESTRICT rhs,
	const int3 N, const FLOAT lambda, const FLOAT epps, const FLOAT alpha, const uchar enforcePositivity) {
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#ifdef CUDA
	if (xyz.x >= N.x || xyz.y >= N.y || xyz.z >= N.z)
#else
	if (any(xyz >= N))
#endif
		return;
	const LTYPE n = (xyz.x) + (xyz.y) * (N.x) + (xyz.z) * (N.x * N.y);
	const FLOAT imOld = im[n];
#ifdef PKMA
	FLOAT imApu = imOld - lambda * rhs[n];
#elif defined(MBSREM)
	FLOAT imApu = imOld + lambda * rhs[n];
#elif defined(BSREM)
	FLOAT imApu = imOld + lambda * rhs[n] * imOld;
#endif
	if (enforcePositivity)
		imApu = fmax((FLOAT)epps, imApu);
#ifdef PKMA
	im[n] = ((FLOAT)1.f - alpha) * imOld  + alpha * imApu;
#elif defined(MBSREM)
	if (imApu >= alpha)
		imApu = alpha - (FLOAT)epps;
	im[n] = imApu;
#elif defined(BSREM)
	im[n] = imApu;
#endif
}
#endif

// PDHG image update
// Similar to above, after backprojection, but for PDHG and its variants
// Different variations for subset and non-subset versions
#if defined(PDHG)
KERNEL3
void PDHGUpdate(CLGLOBAL FLOAT* CLRESTRICT im, const CLGLOBAL FLOAT* CLRESTRICT rhs, CLGLOBAL FLOAT* CLRESTRICT u,
	const int3 N, const float epps, const float theta, const float tau, const uchar enforcePositivity) {
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#ifdef CUDA
	if (xyz.x >= N.x || xyz.y >= N.y || xyz.z >= N.z)
#else
	if (any(xyz >= N))
#endif
		return;
	const LTYPE n = (xyz.x) + (xyz.y) * (N.x) + (xyz.z) * (N.x * N.y);
#ifdef SUBSETS
	FLOAT imApu = im[n];
	imApu -= tau * rhs[n];
	if (enforcePositivity)
		imApu = (FLOAT)(fmax((epps), (float)imApu));
	im[n] = imApu;
#else
	const FLOAT uPrev = u[n];
	FLOAT uNew = uPrev;
	uNew -= tau * rhs[n];
	if (enforcePositivity)
		uNew = (FLOAT)(fmax((epps), (float)uNew));
	u[n] = uNew;
	im[n] = uNew + theta * (uNew - uPrev);
#endif
}
#endif


#ifdef ROTATE
#if defined(USEIMAGES) && defined(OPENCL)
CONSTANT sampler_t samplerRotate = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_LINEAR | CLK_ADDRESS_CLAMP;
#endif
KERNEL
// Initial version from: https://stackoverflow.com/questions/9833316/cuda-image-rotation/10008412#10008412
void rotate(
    CLGLOBAL FLOAT* CLRESTRICT rotim
#ifdef USEIMAGES
    , IMAGE3D im
#else
    , const CLGLOBAL FLOAT* CLRESTRICT im
#endif
    , const int Nx, const int Ny, const int Nz, const float cosa, const float sina
) {
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
	if (xyz.x >= Nx || xyz.y >= Ny || xyz.z >= Nz)
		return;
	const LTYPE n = (xyz.x) + (xyz.y) * (Nx) + (xyz.z) * (Nx * Ny);

    const FLOAT xA = (FLOAT)(xyz.x - Nx/2) + (FLOAT)0.5f;
    const FLOAT yA = (FLOAT)(xyz.y - Ny/2) + (FLOAT)0.5f;

    const FLOAT src_x = (xA * (FLOAT)cosa - yA * (FLOAT)sina + Nx/2) - (FLOAT)0.5f;
    const FLOAT src_y = (xA * (FLOAT)sina + yA * (FLOAT)cosa + Ny/2) - (FLOAT)0.5f;

    if (src_x >= (FLOAT)0.0f && src_x < Nx && src_y >= (FLOAT)0.0f && src_y < Ny) {
		FLOAT val = (FLOAT)0.f;
#ifdef USEIMAGES
#if defined(CUDA)
        val = tex3D<FLOAT>(im, src_x, src_y, xyz.z);
#elif defined(OPENCL)
        val = read_imagef(im, samplerRotate, (int4)(src_x, src_y, xyz.z, 0)).w;
#endif
#else
        // BILINEAR INTERPOLATION
        const int src_x0 = (int)(src_x);
        const int src_x1 = (src_x0 + 1);
        const int src_y0 = (int)(src_y);
        const int src_y1 = (src_y0 + 1);

        const FLOAT sx = (src_x - src_x0);
        const FLOAT sy = (src_y - src_y0);

        const int idx_src00 = min(max(0, src_x0 + src_y0 * Nx), (Nx * Ny) - 1);
        const int idx_src10 = min(max(0, src_x1 + src_y0 * Nx), (Nx * Ny) - 1);
        const int idx_src01 = min(max(0, src_x0 + src_y1 * Nx), (Nx * Ny) - 1);
        const int idx_src11 = min(max(0, src_x1 + src_y1 * Nx), (Nx * Ny) - 1);

        val  = ((FLOAT)1.0f - sx) * ((FLOAT)1.0f - sy) * im[idx_src00 + xyz.z * Nx * Ny];
        val += (       sx) * ((FLOAT)1.0f - sy) * im[idx_src10 + xyz.z * Nx * Ny];
        val += ((FLOAT)1.0f - sx) * (       sy) * im[idx_src01 + xyz.z * Nx * Ny];
        val += (       sx) * (       sy) * im[idx_src11 + xyz.z * Nx * Ny];
#endif
		rotim[n] = val;
    } 
	else {
        rotim[n] = (FLOAT)0.0f;
    }

}
#endif