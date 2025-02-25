
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
#ifdef CUDA
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
	float result = 0.f;
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
	float result = 0.f;
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
	const float div = input[xyz.x];
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
void NLM(CLGLOBAL float* CLRESTRICT grad, IMAGE3D u, CONSTANT float* gaussian, 
#else
void NLM(CLGLOBAL float* CLRESTRICT grad, const CLGLOBAL float* CLRESTRICT u, CONSTANT float* gaussian, 
#endif
	const int3 N, const int3 NOrig, const float h, const float epps, const float beta
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
#ifdef MASKBP3D
	, IMAGE3D maskBP
#else
	, IMAGE2D maskBP
#endif
#endif
#ifdef EFOVZ // Compute only in the voxels of the actual FOV (when using extended FOV)
	, CONSTANT uchar* fovIndices
#endif
) {

	LTYPE3 ii = MINT3(GID0, GID1, GID2);
	const LTYPE n = (ii.x) + (ii.y) * (N.x) + (ii.z) * (N.x * N.y);
	float weight_sum = epps;
	float output = 0.f;
#if NLTYPE == 1
	float outputAla = epps;
#endif
	LTYPE startX = GRID0 * LSIZE0 - SWINDOWX - PWINDOWX + LID0;
	LTYPE startY = GRID1 * LSIZE1 - SWINDOWY - PWINDOWY + LID1;
	LTYPE startZ = GRID2 * LSIZE2 - SWINDOWZ - PWINDOWZ + LID2;
	LTYPE endX = (GRID0 + 1) * LSIZE0 + SWINDOWX + PWINDOWX;
	LTYPE endY = (GRID1 + 1) * LSIZE1 + SWINDOWY + PWINDOWY;
	LTYPE endZ = (GRID2 + 1) * LSIZE2 + SWINDOWZ + PWINDOWZ;
	LOCAL float lCache[SIZEX][SIZEY][SIZEZ];
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
				lCacheRef[indX][indY][indZ] = read_imagef(u_ref, samplerNLM, (int4)(xx, yy, zz, 0)).w;
#endif
#else
				if (xx < 0 || yy < 0 || zz < 0 || xx >= N.x || yy >= N.y || zz >= N.z)
					lCacheRef[indX][indY][indZ] = 0.f;
				else
					lCacheRef[indX][indY][indZ] = u_ref[(xx) + (yy) * N.x + (zz) * N.x * N.y];
#endif
#endif // END NLMREF
#ifdef USEIMAGES
#ifdef CUDA
				lCache[indX][indY][indZ] = tex3D<float>(u, xx, yy, zz);
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
#if defined(NLMADAPTIVE)
	float hh = 0.f;
	const float pSize = CFLOAT((PWINDOWX * 2 + 1) * (PWINDOWY * 2 + 1) * (PWINDOWZ * 2 + 1));
#endif
	const int3 xxyyzz = CMINT3(LID0 + SWINDOWX + PWINDOWX, LID1 + SWINDOWY + PWINDOWY, LID2 + SWINDOWZ + PWINDOWZ);
	const float uj = lCache[xxyyzz.x][xxyyzz.y][xxyyzz.z];
	for (int i = -SWINDOWX; i <= SWINDOWX; i++) {
		for (int j = -SWINDOWY; j <= SWINDOWY; j++) {
			for (int k = -SWINDOWZ; k <= SWINDOWZ; k++) {
				if (i == 0 && j == 0 && k == 0)
					continue;
				float weight = 0.f;
				float distance = 0.f;
// #pragma unroll
				for (int pz = -PWINDOWZ; pz <= PWINDOWZ; pz++) {
// #pragma unroll
					for (int py = -PWINDOWY; py <= PWINDOWY; py++) {
						int dim_g = (pz + PWINDOWZ) * (PWINDOWX * 2 + 1) * (PWINDOWY * 2 + 1) + (py + PWINDOWY) * (PWINDOWX * 2 + 1);
// #pragma unroll
						for (int px = -PWINDOWX; px <= PWINDOWX; px++) {
							const float gg = gaussian[dim_g++];
#ifdef NLMREF
							const float Pk = lCacheRef[xxyyzz.x + i + px][xxyyzz.y + j + py][xxyyzz.z + k + pz];
							const float Pj = lCacheRef[xxyyzz.x + px][xxyyzz.y + py][xxyyzz.z + pz];
#else
							const float Pk = lCache[xxyyzz.x + i + px][xxyyzz.y + j + py][xxyyzz.z + k + pz];
							const float Pj = lCache[xxyyzz.x + px][xxyyzz.y + py][xxyyzz.z + pz];
#endif
							const float PP = Pj - Pk;
							distance += gg * PP * PP;
						}
					}
				}
#if defined(NLMADAPTIVE)
				hh = distance / pSize;
				weight = EXP(-distance / (hh * h + s));
#else
 				weight = EXP(-distance / h);
#endif
 				weight_sum += weight;
				const float uk = lCache[xxyyzz.x + i][xxyyzz.y + j][xxyyzz.z + k];
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
				output += weight * (uabs - uabs / (fabs(u) / gamma + 1.f));
#elif NLTYPE == 6
				// NLGGMRF
				const float delta = uj - uk;
				const float deltapqc = 1.f + POWR(fabs(delta / c), p - q);
				output += weight * (POWR(fabs(delta), p - 1.f) / deltapqc) * (p - gamma * (POWR(fabs(delta), p - q) / deltapqc)) * sign(delta);
#elif NLTYPE == 7
				const float u = (uk - uj);
				const float apu = (u * u + gamma * gamma);
// #ifndef USEMAD // START FMAD
				output += ((2.f * u * u * u) / (apu * apu) - 2.f * (u / apu));
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
	weight_sum = 1.f / weight_sum;
	output *= weight_sum;
#if NLTYPE == 2 // START NLM NLTYPE
	output = uj - output;
#elif NLTYPE == 5
	// Lange with NLMRP
	output = uj - output;
	const float uabs = sign(output);
	output = (uabs - uabs / (fabs(output) / gamma + 1.f));
#elif NLTYPE == 1
#ifndef USEMAD // START FMAD
	output /= SQRT(outputAla * weight_sum + epps);
#else
	output /= SQRT(FMAD(outputAla, weight_sum, epps));
#endif // END FMAD
#endif // END NLM NLTYPE
	grad[n] += beta * output;
}
#endif // END NLM

// Relative difference prior
#ifdef RDP // START RDP
#ifdef OPENCL
#ifdef USEIMAGES
CONSTANT sampler_t samplerRDP = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP;
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
	const int3 N, const int3 NOrig, const float gamma, const float epps, const float beta
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
	, CONSTANT float* weight
#endif
#ifdef RDPREF
#ifdef USEIMAGES
	, IMAGE3D u_ref
#else
	, const CLGLOBAL float* CLRESTRICT u_ref
#endif
#endif
) {

	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#ifdef RDPCORNERS // START RDPCORNERS
	float output = 0.f;
	LTYPE startX = GRID0 * LSIZE0 - SWINDOWX + LID0;
	LTYPE startY = GRID1 * LSIZE1 - SWINDOWY + LID1;
	LTYPE startZ = GRID2 * LSIZE2 - SWINDOWZ + LID2;
	LTYPE endX = (GRID0 + 1) * LSIZE0 + SWINDOWX;
	LTYPE endY = (GRID1 + 1) * LSIZE1 + SWINDOWY;
	LTYPE endZ = (GRID2 + 1) * LSIZE2 + SWINDOWZ;
	LOCAL float lCache[SIZEX][SIZEY][SIZEZ];
#ifdef RDPREF
	LOCAL float lCacheRef[SIZEX][SIZEY][SIZEZ];
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
				lCacheRef[indX][indY][indZ] = tex3D<float>(u_ref, xx, yy, zz);
#else
				lCacheRef[indX][indY][indZ] = read_imagef(u_ref, samplerRDP, (int4)(xx, yy, zz, 0)).w;
#endif
#else
				if (xx < 0 || yy < 0 || zz < 0 || xx >= N.x || yy >= N.y || zz >= N.z)
					lCacheRef[indX][indY][indZ] = 0.f;
				else
					lCacheRef[indX][indY][indZ] = u_ref[(xx) + (yy) * N.x + (zz) * N.x * N.y];
#endif
#endif // END RDPREF
#ifdef USEIMAGES
#ifdef CUDA
				lCache[indX][indY][indZ] = tex3D<float>(u, xx, yy, zz);
#else
				lCache[indX][indY][indZ] = read_imagef(u, samplerRDP, (int4)(xx, yy, zz, 0)).w;
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
// #endif
#ifdef RDPCORNERS // START RDPCORNERS
	const int3 xxyyzz = CMINT3(LID0 + SWINDOWX, LID1 + SWINDOWY, LID2 + SWINDOWZ);
#if defined(RDPREF) // START RDPREF
	const float kj = lCacheRef[xxyyzz.x][xxyyzz.y][xxyyzz.z];
#endif // END RDPREF
	const float uj = lCache[xxyyzz.x][xxyyzz.y][xxyyzz.z];
	int uu = 0;
	for (int i = -SWINDOWX; i <= SWINDOWX; i++) {
		for (int j = -SWINDOWY; j <= SWINDOWY; j++) {
			for (int k = -SWINDOWZ; k <= SWINDOWZ; k++) {
				if (i == 0 && j == 0 && k == 0)
					continue;
				const float uk = lCache[xxyyzz.x + i][xxyyzz.y + j][xxyyzz.z + k];
#if defined(RDPREF) // START RDPREF
				const float kk = lCacheRef[xxyyzz.x + i][xxyyzz.y + j][xxyyzz.z + k];
#endif // END RDPREF
				const float delta = uj - uk;
				const float divPow2 = FMAD(gamma, fabs(delta), uj + uk);
#if defined(RDPREF) // START RDPREF
				output += weight[uu++] * SQRT(kk * kj) * delta * (gamma * fabs(delta) + uj + 3.f * uk + epps * epps) / (divPow2 * divPow2 + epps);
#else
				output += weight[uu++] * delta * (gamma * fabs(delta) + uj + 3.f * uk + epps * epps) / (divPow2 * divPow2 + epps);
#endif // END RDPREF
			}
		}
	}
	if (isnan(output))
		output = 0.f;
	grad[n] += beta * output;
#else
#ifdef USEIMAGES
#ifdef CUDA
	// Current voxel
	const float uj = tex3D<float>(u, xyz.x, xyz.y, xyz.z);
	// Left-right
	const float2 ux = make_float2(tex3D<float>(u, xyz.x + 1, xyz.y, xyz.z), tex3D<float>(u, xyz.x - 1, xyz.y, xyz.z));
	// Top-bottom
	const float2 uy = make_float2(tex3D<float>(u, xyz.x, xyz.y + 1, xyz.z), tex3D<float>(u, xyz.x, xyz.y - 1, xyz.z));
	// Front-back
	const float2 uz = make_float2(tex3D<float>(u, xyz.x, xyz.y, xyz.z + 1), tex3D<float>(u, xyz.x, xyz.y, xyz.z - 1));
#else
	// Current voxel
	const float uj = read_imagef(u, samplerRDP, (int4)(xyz.x, xyz.y, xyz.z, 0)).w;
	// Left-right
	const float2 ux = { read_imagef(u, samplerRDP, (int4)(xyz.x + 1, xyz.y, xyz.z, 0)).w, read_imagef(u, samplerRDP, (int4)(xyz.x - 1, xyz.y, xyz.z, 0)).w };
	// Top-bottom
	const float2 uy = { read_imagef(u, samplerRDP, (int4)(xyz.x, xyz.y + 1, xyz.z, 0)).w, read_imagef(u, samplerRDP, (int4)(xyz.x, xyz.y - 1, xyz.z, 0)).w };
	// Front-back
	const float2 uz = { read_imagef(u, samplerRDP, (int4)(xyz.x, xyz.y, xyz.z + 1, 0)).w, read_imagef(u, samplerRDP, (int4)(xyz.x, xyz.y, xyz.z - 1, 0)).w };
#endif
#else
	// Current voxel
	const float uj = u[xyz.x + xyz.y * N.x + xyz.z * N.x * N.y];
	float2 ux = MFLOAT2(0.f, 0.f);
	float2 uy = MFLOAT2(0.f, 0.f);
	float2 uz = MFLOAT2(0.f, 0.f);
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
void GGMRFKernel(CLGLOBAL float* CLRESTRICT grad, IMAGE3D u, 
#else
void GGMRFKernel(CLGLOBAL float* CLRESTRICT grad, const CLGLOBAL float* CLRESTRICT u, 
#endif
	CONSTANT float* weight, const int3 N, const float p, const float q, const float c, const float pqc, const float beta
) {

	LTYPE3 ii = MINT3(GID0, GID1, GID2);
	const LTYPE n = (ii.x) + (ii.y) * (N.x) + (ii.z) * (N.x * N.y);
	float output = 0.f;
	LTYPE startX = GRID0 * LSIZE0 - SWINDOWX + LID0;
	LTYPE startY = GRID1 * LSIZE1 - SWINDOWY + LID1;
	LTYPE startZ = GRID2 * LSIZE2 - SWINDOWZ + LID2;
	LTYPE endX = (GRID0 + 1) * LSIZE0 + SWINDOWX;
	LTYPE endY = (GRID1 + 1) * LSIZE1 + SWINDOWY;
	LTYPE endZ = (GRID2 + 1) * LSIZE2 + SWINDOWZ;
	LOCAL float lCache[SIZEX][SIZEY][SIZEZ];
	LTYPE indZ = LID2;
	for (LTYPE zz = startZ; zz < endZ; zz += LSIZE2) {
		LTYPE indY = LID1;
		for (LTYPE yy = startY; yy < endY; yy += LSIZE1) {
			LTYPE indX = LID0;
			for (LTYPE xx = startX; xx < endX; xx += LSIZE0) {
#ifdef USEIMAGES
#ifdef CUDA
				lCache[indX][indY][indZ] = tex3D<float>(u, xx, yy, zz);
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
	const float uj = lCache[xxyyzz.x][xxyyzz.y][xxyyzz.z];
	int uu = 0;
	for (int i = -SWINDOWX; i <= SWINDOWX; i++) {
		for (int j = -SWINDOWY; j <= SWINDOWY; j++) {
			for (int k = -SWINDOWZ; k <= SWINDOWZ; k++) {
				if (i == 0 && j == 0 && k == 0)
					continue;
				const float uk = lCache[xxyyzz.x + i][xxyyzz.y + j][xxyyzz.z + k];
				const float delta = uj - uk;
				const float deltapqc = 1.f + POWR(fabs(delta / c), p - q);
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
void medianFilter3D(const CLGLOBAL float* grad, CLGLOBAL float* output, const int3 N, const int3 NOrig
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
	float median[KOKO];
	float medianF[KOKO];
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
	const float scale = fmax(1.f, length(apu) / alpha);
#else
// L1 norm
	const float scale = fmax(fmax(fabs(apu.z), fmax(fabs(apu.x), fabs(apu.y))) / alpha, 1.f);
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
	const float3 apu = MFLOAT3(inputX[idx], inputY[idx], 0.f);
	const float3 apu2 = MFLOAT3(input2XY[idx], 0.f, 0.f);
#endif
#ifdef L2 // START L2
	const float scale = fmax(1.f, SQRT(apu.x * apu.x + apu.y * apu.y + apu.z * apu.z + (apu2.x * apu2.x) * 2.f + (apu2.y * apu2.y) * 2.f + (apu2.z * apu2.z) * 2.f) / alpha);
#else
	const float scale = fmax(fmax(fabs(apu2.z),fmax(fabs(apu2.y), fmax(fabs(apu2.x), fmax(fabs(apu.z), fmax(fabs(apu.x), fabs(apu.y)))))) / alpha, 1.f);
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
	float apuVal = 0.f;
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
	float apuVal = 0.f;
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
void ProxTGVSymmDeriv(const int3 N, const int3 NOrig, const CLGLOBAL float* CLRESTRICT vX, const CLGLOBAL float* CLRESTRICT vY, 
#ifdef TGVZ
	const CLGLOBAL float* CLRESTRICT vZ, CLGLOBAL float* qX, CLGLOBAL float* qY, CLGLOBAL float* qZ, CLGLOBAL float* q2XY, CLGLOBAL float* q2XZ, CLGLOBAL float* q2YZ, 
#else
	CLGLOBAL float* qX, CLGLOBAL float* qY, CLGLOBAL float* q2XY, 
#endif
	const float sigma2
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
	float apuVal = 0.f;
	float imApuX = vX[x];
// Forward difference
#if DIFFTYPE == 0 // START DIFFTYPE == 0
	forwardDiffX2(&apuVal, xyz, NOrig, vX, imApuX);
	apuVal *= sigma2;
	qX[x] += apuVal;
/////////////// Y ///////////////
	apuVal = 0.f;
	float imApuY = vY[x];
	forwardDiffY2(&apuVal, xyz, NOrig, vY, imApuY);
	apuVal *= sigma2;
	qY[x] += apuVal;
#ifdef TGVZ
/////////////// Z ///////////////
	apuVal = 0.f;
	float imApuZ = vZ[x];
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
	float imApuY = vY[x];
	backwardDiffY2(&apuValY, xyz, NOrig, vY, imApuY);
	apuValY *=  sigma2;
	qY[x] += apuValY;
#ifdef TGVZ
/////////////// Z ///////////////
	apuVal = 0.f;
	float imApuZ = vZ[x];
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
void ProxTGVDivergence(const int3 N, const int3 NOrig, const CLGLOBAL float* CLRESTRICT qX, const CLGLOBAL float* CLRESTRICT qY, 
#ifdef TGVZ
	const CLGLOBAL float* CLRESTRICT qZ, const CLGLOBAL float* CLRESTRICT q2XY, const CLGLOBAL float* CLRESTRICT q2XZ, const CLGLOBAL float* CLRESTRICT q2YZ, 
	CLGLOBAL float* vX, CLGLOBAL float* vY, CLGLOBAL float* vZ, 
#else
	const CLGLOBAL float* CLRESTRICT q2XY, CLGLOBAL float* vX, CLGLOBAL float* vY, 
#endif
	const CLGLOBAL float* CLRESTRICT pX, const CLGLOBAL float* CLRESTRICT pY, const CLGLOBAL float* CLRESTRICT pZ, const float theta, const float tau
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
	float apuVal = 0.f;
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
		float v1 = vApu + tau * 1.f * (pX[x] + apuVal);
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
		float vApu = vX[x];
#ifdef USEMAD
		float v1 = FMAD(tau, pX[x] + apuVal, vApu);
		vX[x] = FMAD(theta, v1 - vApu, v1);
#else
		float v1 = vApu + tau * (pX[x] + apuVal);
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
void hyperbolicKernel(CLGLOBAL float* CLRESTRICT grad, IMAGE3D u, 
#else
void hyperbolicKernel(CLGLOBAL float* CLRESTRICT grad, const CLGLOBAL float* CLRESTRICT u,
#endif
	const int3 N, const int3 NOrig, const float sigma, const float epps, const float beta, CONSTANT float* w
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
	float output = 0.f;
	LTYPE startX = GRID0 * LSIZE0 - SWINDOWX + LID0;
	LTYPE startY = GRID1 * LSIZE1 - SWINDOWY + LID1;
	LTYPE startZ = GRID2 * LSIZE2 - SWINDOWZ + LID2;
	LTYPE endX = (GRID0 + 1) * LSIZE0 + SWINDOWX;
	LTYPE endY = (GRID1 + 1) * LSIZE1 + SWINDOWY;
	LTYPE endZ = (GRID2 + 1) * LSIZE2 + SWINDOWZ;
	LOCAL float lCache[SIZEX][SIZEY][SIZEZ];
	LTYPE indZ = LID2;
	for (LTYPE zz = startZ; zz < endZ; zz += LSIZE2) {
		LTYPE indY = LID1;
		for (LTYPE yy = startY; yy < endY; yy += LSIZE1) {
			LTYPE indX = LID0;
			for (LTYPE xx = startX; xx < endX; xx += LSIZE0) {
#ifdef USEIMAGES
#ifdef CUDA
				lCache[indX][indY][indZ] = tex3D<float>(u, xx, yy, zz);
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
	if (any(xyz >= N))
		return;
	const int3 xxyyzz = CMINT3(LID0 + SWINDOWX, LID1 + SWINDOWY, LID2 + SWINDOWZ);
	const float uj = lCache[xxyyzz.x][xxyyzz.y][xxyyzz.z];
	int uu = 0;
	for (int k = -SWINDOWZ; k <= SWINDOWZ; k++) {
		for (int j = -SWINDOWY; j <= SWINDOWY; j++) {
			for (int i = -SWINDOWX; i <= SWINDOWX; i++) {
				if (i == 0 && j == 0 && k == 0)
					continue;
				const float u = lCache[xxyyzz.x + i][xxyyzz.y + j][xxyyzz.z + k];
				const float ux = (uj - u) / sigma;
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

inline float sqrtVal(const float3 input, const float epps
#ifdef TVW1
	, const float3 w
#endif
) {
#ifdef TVW1
#ifdef USEMAD
	return sqrt(FMAD(w.x, input.x * input.x, FMAD(w.y, input.y * input.y, FMAD(w.z, input.z * input.z, epps))));
#else
	return sqrt(w.x * input.x * input.x + w.y * input.y * input.y + w.z * input.z * input.z + epps);
#endif
#else
#ifdef USEMAD
	return sqrt(FMAD(input.x, input.x, FMAD(input.y, input.y, FMAD(input.z, input.z, epps))));
#else
	return sqrt(input.x * input.x + input.y * input.y + input.z * input.z + epps);
#endif
#endif
}

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
	const int3 N, const int3 NOrig, const float sigma, const float epps, const float beta
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
	, const float C
#endif
#if defined(ANATOMICAL1) || defined(ANATOMICAL2) || defined(ANATOMICAL3)
	, CLGLOBAL float* CLRESTRICT S
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
#ifdef USEIMAGES
	const float uijk = read_imagef(u, samplerTV, (int4)(xyz.x, xyz.y, xyz.z, 0)).w;
#else
	const float uijk = u[n];
#endif
#if defined(SATV) // START JPTV || SATV
#ifdef USEIMAGES
#ifdef CUDA
	float2 ux = make_float2(tex3D<float>(u, xyz.x + 1, xyz.y, xyz.z), tex3D<float>(u, xyz.x - 1, xyz.y, xyz.z));
	float2 uy = make_float2(tex3D<float>(u, xyz.x, xyz.y + 1, xyz.z), tex3D<float>(u, xyz.x, xyz.y - 1, xyz.z));
	float2 uz = make_float2(tex3D<float>(u, xyz.x, xyz.y, xyz.z + 1), tex3D<float>(u, xyz.x, xyz.y, xyz.z - 1));
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
	const float3 uijkP = make_float3(tex3D<float>(u, xyz.x + 1, xyz.y, xyz.z), tex3D<float>(u, xyz.x, xyz.y + 1, xyz.z), tex3D<float>(u, xyz.x, xyz.y, xyz.z + 1));
	const float3 uijkM = make_float3(tex3D<float>(u, xyz.x - 1, xyz.y, xyz.z), tex3D<float>(u, xyz.x, xyz.y - 1, xyz.z), tex3D<float>(u, xyz.x, xyz.y, xyz.z - 1));
	const float2 ui = make_float2(tex3D<float>(u, xyz.x - 1, xyz.y + 1, xyz.z), tex3D<float>(u, xyz.x - 1, xyz.y, xyz.z + 1));
	const float2 uj = make_float2(tex3D<float>(u, xyz.x + 1, xyz.y - 1, xyz.z), tex3D<float>(u, xyz.x, xyz.y - 1, xyz.z + 1));
	const float2 uk = make_float2(tex3D<float>(u, xyz.x + 1, xyz.y, xyz.z - 1), tex3D<float>(u, xyz.x, xyz.y + 1, xyz.z - 1));
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
	w4 = EXP(-w4 * w4);
	const float pvalijk = sqrtVal(u4, epps, w4);
	float3 w1 = (u1) / sigma;
	w1 = EXP(-w1 * w1);
	float3 w2 = (u2) / sigma;
	w2 = EXP(-w2 * w2);
	float3 w3 = (u3) / sigma;
	w3 = EXP(-w3 * w3);
#ifdef USEMAD
	grad[n] += beta * (-(FMAD(w4.x, u4.x, FMAD(w4.y, u4.y, w4.z * u4.z))) / pvalijk + (w1.x * (uijk - uijkM.x)) / sqrtVal(u1, epps, w1) + (w2.y * (uijk - uijkM.y)) / sqrtVal(u2, epps, w2) + (w3.y * (uijk - uijkM.z)) / sqrtVal(u3, epps, w3));
#else
	grad[n] += beta * (-(w4.x * u4.x + w4.y * u4.y + w4.z * u4.z) / pvalijk + (w1.x * (uijk - uijkM.x)) / sqrtVal(u1, epps, w1) + (w2.y * (uijk - uijkM.y)) / sqrtVal(u2, epps, w2) + (w3.y * (uijk - uijkM.z)) / sqrtVal(u3, epps, w3));
#endif
#else
#ifdef ANATOMICAL1 // TV type 1
	const LTYPE NN = N.x * N.y * N.z;
	__private float s[9];
	for (int kk = 0; kk < 9; kk++)
		s[kk] = S[n + NN * kk];
	const float3 val = uijkP - uijk;
	const float pvalijk = sqrt(val.x * val.x * s[0] + val.y * val.y * s[4] + val.z * val.z * s[8] + s[1] * (val.x) * (val.y) + s[3] * (val.x) * (val.y) + s[2] * (val.x) * (val.z) + s[6] * (val.x) * (val.z) + 
		s[5] * (val.y) * (val.z) + s[7] * (val.y) * (val.z) + epps);
	const float pvalijkX = sqrt(u1.x * u1.x * s[0] + u1.y * u1.y * s[4] + u1.z * u1.z * s[8] + s[1] * (u1.x) * (u1.y) + s[3] * (u1.x) * (u1.y) + s[2] * (u1.x) * (u1.z) + s[6] * (u1.x) * (u1.z) + 
		s[5] * (u1.y) * (u1.z) + s[7] * (u1.y) * (u1.z) + epps);
	const float pvalijkY = sqrt(u2.x * u2.x * s[0] + u2.y * u2.y * s[4] + u2.z * u2.z * s[8] + s[1] * (u2.x) * (u2.y) + s[3] * (u2.x) * (u2.y) + s[2] * (u2.x) * (u2.z) + s[6] * (u2.x) * (u2.z) + 
		s[5] * (u2.y) * (u2.z) + s[7] * (u2.y) * (u2.z) + epps);
	const float pvalijkZ = sqrt(u3.x * u3.x * s[0] + u3.y * u3.y * s[4] + u3.z * u3.z * s[8] + s[1] * (u3.x) * (u3.y) + s[3] * (u3.x) * (u3.y) + s[2] * (u3.x) * (u3.z) + s[6] * (u3.x) * (u3.z) + 
		s[5] * (u3.y) * (u3.z) + s[7] * (u3.y) * (u3.z) + epps);
	const float dx = s[0] * (2.f * (uijk - uijkM.x)) + s[3] * u1.y + s[2] * u1.z + s[6] * u1.z + s[1] * u1.y;
	const float dy = s[4] * (2.f * (uijk - uijkM.y)) + s[5] * u2.z + s[3] * u2.x + s[1] * u2.x + s[7] * u2.z;
	const float dz = s[8] * (2.f * (uijk - uijkM.z)) + s[6] * u3.x + s[5] * u3.y + s[7] * u3.y + s[2] * u3.x;
	const float d = s[1] * val.x + s[2] * val.x + s[3] * val.x + s[6] * val.x + s[1] * val.y + s[3] * val.y + s[5] * val.y + s[7] * val.y + s[2] * val.z + s[5] * val.z + s[6] * val.z + s[7] * val.z + s[0] * 2.f * val.x + s[4] * 2.f * val.y + s[8] * 2.f * val.z;
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
	const float pvalijk = native_sqrt(apu.x * apu.x + apu.y * apu.y + apu.z * apu.z + C * (apuS.x * apuS.x + apuS.y * apuS.y + apuS.z * apuS.z) + epps);
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
	epsilon = epsilon / native_sqrt(epsilon.x * epsilon.x + epsilon.y * epsilon.y + epsilon.z * epsilon.z + C * C);
	const float3 apu = uijkP - uijk;
	const float apuR = uijkR.x * apu.x + uijkR.y * apu.y + uijkR.z * apu.z;
	const float pvalijk = native_sqrt(apu.x * apu.x + apu.y * apu.y + apu.z * apu.z - apuR * apuR + epps);
	float apuRXYZ = uijkR.x * u1.x + uijkR.y * u1.y + uijkR.z * u1.z;
	const float pvalijkX = native_sqrt(u1.x * u1.x + u1.y * u1.y + u1.z * u1.z + apuRXYZ * apuRXYZ + epps);
	apuRXYZ = uijkR.x * u2.x + uijkR.y * u2.y + uijkR.z * u2.z;
	const float pvalijkY = native_sqrt(u2.x * u2.x + u2.y * u2.y + u2.z * u2.z + apuRXYZ * apuRXYZ + epps);
	apuRXYZ = uijkR.x * u3.x + uijkR.y * u3.y + uijkR.z * u3.z;
	const float pvalijkZ = native_sqrt(u3.x * u3.x + u3.y * u3.y + u3.z * u3.z + apuRXYZ * apuRXYZ + epps);
	grad[n] += beta * .5f * ((6.f * uijk - 2.f * uijkP.x - 2.f * uijkP.y - 2.f * uijkP.z + 2.f * (epsilon.x*(uijk - uijkP.x) + epsilon.y*(uijk - uijkP.y) + epsilon.z*(uijk - uijkP.z)) * (epsilon.x + epsilon.y + epsilon.z)) / pvalijk + 
		2.f * (u1.x - epsilon.x * (epsilon.x * u1.x + epsilon.y * u1.y + epsilon.z * u1.z)) / pvalijkX + 2.f * (u2.y - epsilon.y * (epsilon.x * u2.x + epsilon.y * u2.y + epsilon.z * u2.z)) / pvalijkY + 
		2.f * (u3.z - epsilon.z * (epsilon.x * u3.x + epsilon.y * u3.y + epsilon.z * u3.z))/ pvalijkZ + 1e-7f);
#else // Non-reference image TV
	const float pvalijk = sqrtVal(uijkP - uijk, epps);
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
void PoissonUpdate(CLGLOBAL float* CLRESTRICT im, const CLGLOBAL float* CLRESTRICT rhs,
	const int3 N, const float lambda, const float epps, const float alpha, const uchar enforcePositivity) {
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
#ifdef CUDA
	if (xyz.x >= N.x || xyz.y >= N.y || xyz.z >= N.z)
#else
	if (any(xyz >= N))
#endif
		return;
	const LTYPE n = (xyz.x) + (xyz.y) * (N.x) + (xyz.z) * (N.x * N.y);
	const float imOld = im[n];
#ifdef PKMA
	float imApu = imOld - lambda * rhs[n];
#elif defined(MBSREM)
	float imApu = imOld + lambda * rhs[n];
#elif defined(BSREM)
	float imApu = imOld + lambda * rhs[n] * imOld;
#endif
	if (enforcePositivity)
		imApu = fmax(epps, imApu);
#ifdef PKMA
	im[n] = (1.f - alpha) * imOld  + alpha * imApu;
#elif defined(MBSREM)
	if (imApu >= alpha)
		imApu = alpha - epps;
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
void PDHGUpdate(CLGLOBAL float* CLRESTRICT im, const CLGLOBAL float* CLRESTRICT rhs, CLGLOBAL float* CLRESTRICT u,
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
	float imApu = im[n];
	imApu -= tau * rhs[n];
	if (enforcePositivity)
		imApu = fmax(epps, imApu);
	im[n] = imApu;
#else
	const float uPrev = u[n];
	float uNew = uPrev;
	uNew -= tau * rhs[n];
	if (enforcePositivity)
		uNew = fmax(epps, uNew);
	u[n] = uNew;
	im[n] = uNew + theta * (uNew - uPrev);
#endif
}
#endif


#ifdef ROTATE
KERNEL
// Initial version from: https://stackoverflow.com/questions/9833316/cuda-image-rotation/10008412#10008412
void rotate(CLGLOBAL float* CLRESTRICT rotim, const CLGLOBAL float* CLRESTRICT im, const int Nx, const int Ny, const int Nz, const float cosa, const float sina) {
	LTYPE3 xyz = MINT3(GID0, GID1, GID2);
	if (xyz.x >= Nx || xyz.y >= Ny || xyz.z >= Nz)
		return;
	const LTYPE n = (xyz.x) + (xyz.y) * (Nx) + (xyz.z) * (Nx * Ny);

    const float xA = (float)(xyz.x - Nx/2) + 0.5f;
    const float yA = (float)(xyz.y - Ny/2) + 0.5f;

    const float src_x = (xA * cosa - yA * sina + Nx/2) - 0.5f;
    const float src_y = (xA * sina + yA * cosa + Ny/2) - 0.5f;

    if (src_x >= 0.0f && src_x < Nx && src_y >= 0.0f && src_y < Ny) {
        // BILINEAR INTERPOLATION
        const int src_x0 = (int)(src_x);
        const int src_x1 = (src_x0 + 1);
        const int src_y0 = (int)(src_y);
        const int src_y1 = (src_y0 + 1);

        const float sx = (src_x - src_x0);
        const float sy = (src_y - src_y0);

        const int idx_src00 = min(max(0, src_x0 + src_y0 * Nx), (Nx * Ny) - 1);
        const int idx_src10 = min(max(0, src_x1 + src_y0 * Nx), (Nx * Ny) - 1);
        const int idx_src01 = min(max(0, src_x0 + src_y1 * Nx), (Nx * Ny) - 1);
        const int idx_src11 = min(max(0, src_x1 + src_y1 * Nx), (Nx * Ny) - 1);

		float val = 0.f;

        val  = (1.0f - sx) * (1.0f - sy) * im[idx_src00 + xyz.z * Nx * Ny];
        val += (       sx) * (1.0f - sy) * im[idx_src10 + xyz.z * Nx * Ny];
        val += (1.0f - sx) * (       sy) * im[idx_src01 + xyz.z * Nx * Ny];
        val += (       sx) * (       sy) * im[idx_src11 + xyz.z * Nx * Ny];
		rotim[n] = val;
    } 
	else {
        rotim[n] = 0.0f;
    }

}
#endif