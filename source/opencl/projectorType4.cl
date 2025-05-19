
/*******************************************************************************************************************************************
* Matrix free projectors for forward and backward projections. This is the interpolation-based projector and contains a ray-based forward
* projector and a voxel-based backprojector. The ray-based forward projector can be used with any data while the backprojector should be
* used for CT-like data only (i.e. projection-type data). The ray-based forward projector does support backprojection as well but as with
* projector types 1-3, it is not efficient backprojector, though it is completely adjoint. The voxel-based backprojector is not exactly
* adjoint with the forward projection, but the difference is less than 1%.
*
* Used by implementations 2, 3 and 5.
*
* The forward and backward projections are separate kernels. 
*
* Compiler preprocessing is utilized heavily, for example all the corrections are implemented as compiler preprocesses. The code for 
* specific correction is thus only applied if it has been selected. The kernels are always compiled on-the-fly, though when using same input 
* parameters the kernel should be loaded from cache leading to a slightly faster startup time.
* 
* The ray-based forward projection projects a ray from the source to each detector pixel. Once the ray enters the FOV, the values are
* linearly interpolated at each d_L step, weighted by the fixed length of d_L. The ray-based forward projector works for all data.
*
* INPUTS:
* d_size_x = the number of detector elements (rows),
* d_sizey = the number of detector elements (columns),
* d_dPitch = Either a vector of float2 or two floats if PYTHON is defined, the detector size/pitch in both "row" and "column" directions
* d_L = Interpolation length, i.e. the length that is moved everytime the interpolation is done, forward projection only
* global_factor = a global correction factor, e.g. dead time, can be simply 1.f
* maskFP = 2D/3D Forward projection mask, i.e. LORs/measurements with 0 will be skipped
* maskBP = 2D/3D backward projection mask, i.e. voxels with 0 will be skipped
* d_TOFCenter = Offset of the TOF center from the first center of the FOV,
* sigma_x = TOF STD, 
* d_atten = attenuation data (images, if USEIMAGES is defined, or buffer),
* d_N = image size in x/y/z- dimension, uint3 or three uints (if PYTHON is defined),
* d_b = distance from the pixel space to origin (z/x/y-dimension), float3 or three floats (if PYTHON is defined),
* d_bmax = part in parenthesis of equation (9) in [1] precalculated when k = Nz, float3 or three floats (if PYTHON is defined),
* d_scale = precomputed scaling value, float3 or three floats (if PYTHON is defined), see computeProjectorScalingValues.m
* d_OSEM = image for current estimates or the input buffer for backward projection,
* d_output = forward or backward projection,
* d_xy/z = detector x/y/z-coordinates,
* rings = Number of detector rings, PET only and only when computing listmode sensitivity image,
* d_det_per_ring = number of detectors per ring, (only for listmode data sensitivity image computation, can be any value otherwise)
* d_norm = normalization coefficients,
* d_scat = scatter coefficients when using the system matrix method (multiplication), 
* d_Summ = buffer for d_Summ (sensitivity image),
* no_norm = If 1, sensitivity image is not computed,
* m_size = Total number of LORs/measurements for this subset,
* currentSubset = current subset
* aa = The current volume, for multi-resolution reconstruction, 0 can be used when not using multi-resolution
*
* OUTPUTS:
* d_output = forward projection,
*
* Copyright (C) 2022-2025 Ville-Veikko Wettenhovi
*
* This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it wiL be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************************************************************************/

#ifndef NVOXELS
#define NVOXELS 1
#endif
//#define NVOXELS 8
//#define PITCH

#if defined(FP) || (defined(BP) && !defined(CT))// START FP
#define NSTEPS 10000
// #undef PITCH
// #undef NA
// #define NA 6

// Forward projection
KERNEL
void projectorType4Forward(const uint d_size_x, const uint d_sizey, 
#ifdef PYTHON
	const float d_dPitchX, const float d_dPitchY, 
#else
	const float2 d_dPitch, 
#endif
    const float dL, const float global_factor, 
	///////////////////////// TOF BINS /////////////////////////
#ifdef TOF
	CONSTANT float* TOFCenter, const float sigma_x, 
#endif
	///////////////////////// END TOF BINS /////////////////////////
    ////////////////////////////////////////////////////////////////////////
#if !defined(CT) && defined(ATN) && !defined(ATNM)
	IMAGE3D d_atten,
#elif !defined(CT) && !defined(ATN) && defined(ATNM)
	const CLGLOBAL float* CLRESTRICT d_atten,
#endif
    ////////////////////////////////////////////////////////////////////////
#ifdef PYTHON
	const uint d_Nx, const uint d_Ny, const uint d_Nz, const float bx, const float by, const float bz, 
    const float d_bmaxx, const float d_bmaxy, const float d_bmaxz, const float d_scalex, const float d_scaley, const float d_scalez,
#else
	const uint3 d_N, const float3 b, const float3 bmax, const float3 d_scale, 
#endif
#ifdef FP
    IMAGE3D d_OSEM, CLGLOBAL float* CLRESTRICT d_output,
#else
    const CLGLOBAL float* CLRESTRICT d_OSEM, CLGLOBAL CAST* CLRESTRICT d_output,
#endif
#if !defined(USEGLOBAL)
	CONSTANT float* d_xyz,
#else
	const CLGLOBAL float* CLRESTRICT d_xyz,
#endif
#if (defined(LISTMODE) && !defined(SENS) && !defined(INDEXBASED))
    const CLGLOBAL float* CLRESTRICT d_uv,
#else
    CONSTANT float* d_uv,
#endif
#ifdef SENS
	const int rings, const uint d_det_per_ring,
#endif
#ifdef MASKFP
#ifdef MASKFP3D
    IMAGE3D maskFP,
#else
    IMAGE2D maskFP,
#endif
#endif
#if defined(MASKBP) && (defined(BP) || defined(SENS)) && !defined(CT)
#ifdef MASKBP3D
    IMAGE3D maskBP,
#else
    IMAGE2D maskBP,
#endif
#endif
    const LONG d_nProjections,
#if defined(SUBSETS) && !defined(LISTMODE)
    const CLGLOBAL uint* CLRESTRICT d_xyindex, const CLGLOBAL ushort* CLRESTRICT d_zindex,
#endif
#if defined(INDEXBASED) && defined(LISTMODE) && !defined(SENS)
	const CLGLOBAL ushort* CLRESTRICT trIndex, const CLGLOBAL ushort* CLRESTRICT axIndex,
#endif
#if defined(LISTMODE) && defined(TOF)
	const CLGLOBAL uchar* CLRESTRICT TOFIndex, 
#endif
#ifdef RAW
    const CLGLOBAL ushort* CLRESTRICT d_L, const uint d_det_per_ring,
#endif
	///////////////////////// PET NORMALIZATION DATA /////////////////////////
#ifdef NORM
    const CLGLOBAL float* CLRESTRICT d_norm,
#endif
	///////////////////////// END PET NORMALIZATION DATA /////////////////////////
	///////////////////////// EXTRA CORRECTION DATA /////////////////////////
#ifdef SCATTER
    const CLGLOBAL float* CLRESTRICT d_scat,
#endif
	///////////////////////// END EXTRA CORRECTION DATA /////////////////////////
#if defined(BP) && !defined(CT)
    CLGLOBAL CAST* CLRESTRICT d_Summ,
#endif
	const uchar no_norm, const ULONG m_size, const uint currentSubset, const int aa)
{
    //
    int3 i = MINT3(GID0, GID1, GID2);
#if STYPE == 1 || STYPE == 2 || STYPE == 4 || STYPE == 5
    getIndex(&i, d_size_x, d_sizey, currentSubset);
#endif

#if (defined(CT) || defined(SPECT) || defined(PET)) && !defined(LISTMODE)
    size_t idx = GID0 + GID1 * d_size_x + GID2 * d_sizey * d_size_x;
    if (i.x >= d_size_x || i.y >= d_sizey || i.z >= d_nProjections)
#elif defined(SENS)
	int2 indz;
	indz.x = i.z / rings;
	indz.y = i.z % rings;
	size_t idx = i.x + i.y * d_det_per_ring + i.z * d_det_per_ring * d_det_per_ring;
	if (i.x >= d_det_per_ring || i.y >= d_det_per_ring || i.z >= rings * rings)
#elif defined(LISTMODE)
	size_t idx = GID0;
	if (idx >= m_size)
#else
    size_t idx = GID0 + GID1 * GSIZE0 + GID2 * GSIZE1 * GSIZE0;
    if (i.x >= d_size_x || idx >= m_size)
#endif
        return; 
#ifdef MASKFP
#ifdef CUDA
#ifdef MASKFP3D
	const int maskVal = tex3D<unsigned char>(maskFP, i.x, i.y, i.z);
#else
    const int maskVal = tex2D<unsigned char>(maskFP, i.x, i.y);
#endif
#else
#ifdef MASKFP3D
    const int maskVal = read_imageui(maskFP, sampler_MASK, (int4)(i.x, i.y, i.z, 0)).w;
#else
    const int maskVal = read_imageui(maskFP, sampler_MASK, (int2)(i.x, i.y)).w;
#endif
#endif
    if (maskVal == 0)
        return;
#endif
#ifdef PYTHON
	const float2 d_dPitch = make_float2(d_dPitchX, d_dPitchY);
	const uint3 d_N = make_uint3(d_Nx, d_Ny, d_Nz);
	const float3 d_scale = make_float3(d_scalex, d_scaley, d_scalez);
	const float3 b = make_float3(bx, by, bz);
	const float3 bmax = make_float3(d_bmaxx, d_bmaxy, d_bmaxz);
#endif
#if defined(LISTMODE) && defined(TOF)
	const int TOFid = TOFIndex[idx];
#endif
#if defined(N_RAYS) && defined(FP)
	float ax[NBINS * N_RAYS];
#else
	float ax[NBINS];
#endif
#if defined(BP) && !defined(CT)
#if defined(LISTMODE) && defined(TOF) && !defined(SENS)
	for (int to = 0; to < NBINS; to++)
		ax[to] = 0.f;
	ax[TOFid] = d_OSEM[idx];
#else
#ifndef __CUDACC__ 
#pragma unroll NBINS
#endif
	for (int to = 0; to < NBINS; to++)
#ifdef SENS
		ax[to] = 1.f;
#else
		ax[to] = d_OSEM[idx + to * m_size];
#endif
#endif
#else
#if defined(N_RAYS) && defined(FP)
#ifndef __CUDACC__ 
#pragma unroll NBINS * N_RAYS
#endif
	for (int to = 0; to < NBINS * N_RAYS; to++)
		ax[to] = 0.f;
#else
#ifndef __CUDACC__ 
#pragma unroll NBINS
#endif
	for (int to = 0; to < NBINS; to++)
		ax[to] = 0.f;
#endif
#endif

#ifndef CT

#ifdef NORM // Normalization included
	float local_norm = 0.f;
	local_norm = d_norm[idx];
#endif
#ifdef SCATTER // Scatter data included
	float local_scat = 0.f;
	local_scat = d_scat[idx];
#endif
#endif

#ifdef N_RAYS //////////////// MULTIRAY ////////////////
	int lor = -1;
	// Load the next detector index
	// raw list-mode data
	for (int lorZ = 0u; lorZ < N_RAYS3D; lorZ++) {
		for (int lorXY = 0u; lorXY < N_RAYS2D; lorXY++) {
			lor++;
#endif  //////////////// END MULTIRAY ////////////////
	float3 s, d;
#if (defined(CT) || defined(SPECT)) && !defined(LISTMODE) && !defined(PET) // CT data
	getDetectorCoordinatesCT(d_xyz, d_uv, &s, &d, i, d_size_x, d_sizey, d_dPitch);
#elif defined(LISTMODE) && !defined(SENS) // Listmode data
#if defined(INDEXBASED)
	getDetectorCoordinatesListmode(d_xyz, d_uv, trIndex, axIndex, &s, &d, idx
#else
	getDetectorCoordinatesListmode(d_xyz, &s, &d, idx
#endif
#if defined(N_RAYS)
		, lorXY, lorZ, d_dPitch
#endif
	);
#elif defined(RAW) || defined(SENS) // raw data
	getDetectorCoordinatesRaw(d_xyz, d_uv, i, &s, &d, indz
#if defined(N_RAYS)
		, lorXY, lorZ, d_dPitch
#endif
	);
#elif !defined(SUBSETS) // Precomputation phase, subset types 1, 2, 4, 5, 8, 9, 10, 11
	getDetectorCoordinatesFullSinogram(d_size_x, i, &s, &d, d_xyz, d_uv
#if defined(N_RAYS)
		, lorXY, lorZ, d_dPitch
#endif
#if defined(NLAYERS)
		, d_sizey
#endif
	);
#else // Subset types 3, 6, 7
	getDetectorCoordinates(d_xyindex, d_zindex, idx, &s, &d, d_xyz, d_uv
#if defined(N_RAYS)
		, lorXY, lorZ, d_dPitch
#endif
#if defined(NLAYERS)
		, d_sizey, d_size_x
#endif
	);
#endif
    float3 v = d - s;
    const float3 bmin = b;
    const float3 tBack = (bmin - s) / v;
    const float3 tFront = (bmax - s) / v;

    const float3 tMin = fmin(tFront, tBack);
    const float3 tMax = fmax(tFront, tBack);

    const float tStart = fmax(fmax(tMin.x, tMin.y), tMin.z);
    const float tEnd = fmin(fmin(tMax.x, tMax.y), tMax.z);
#ifdef TOF //////////////// TOF ////////////////
    float TOFSum = 0.f;
#endif //////////////// END TOF ////////////////
    if (tStart >= tEnd)
#ifdef N_RAYS //////////////// MULTIRAY ////////////////
		continue;
#else
		return;
#endif  //////////////// END MULTIRAY ////////////////
	float L = length(v);
#ifndef CT
#if !defined(TOTLENGTH)
	float LL = 0.f;
#endif
#ifdef ATN // Attenuation included
	float jelppi = 0.f;
#endif
#ifdef TOF //////////////// TOF ////////////////
	float D = 0.f;
	float DD = 0.f;
	TOFDis(v, tStart, L, &D, &DD);
	float TOFWeights[NBINS];
#endif //////////////// END TOF ////////////////
#endif
    const float tStep = DIVIDE(dL, L);

    s = (s - bmin) * d_scale;
    v *= d_scale;
    float temp = 0.f;

    float t = tStart;
#if (defined(ATN) && defined(BP)) || (defined(BP) && !defined(TOTLENGTH) && !defined(CT))
    for (uint ii = 0; ii < NSTEPS; ii++) {
#if (defined(ATN) && defined(BP))
#ifndef USEMAD
        const float3 p = t * v + s;
#else
        const float3 p = FMAD3(t, v, s);
#endif
		compute_attenuation(dL, p, d_atten, &jelppi, aa);
#else
        LL += dL;
#endif
        t += tStep;
        if (t >= tEnd)
            break;
    }
    t = tStart;
#endif
#if !defined(CT) //////////////// PET ////////////////
#if defined(TOTLENGTH) //////////////// TOTLENGTH ////////////////
#ifdef N_RAYS //////////////// MULTIRAY ////////////////
		temp = 1.f / (L * CFLOAT(N_RAYS));
#else
		temp = 1.f / L;
#endif //////////////// END MULTIRAY ////////////////
#elif !defined(TOTLENGTH) && defined(BP) //////////////// NOTTOTLENGTH+BP ////////////////
		if (LL == 0.f)
			LL = L;
#ifdef N_RAYS  //////////////// MULTIRAY ////////////////
		temp = 1.f / (LL * CFLOAT(N_RAYS));
#else //////////////// SINGLERAY ////////////////
		temp = 1.f / LL;
#endif //////////////// END MULTIRAY ////////////////
#endif //////////////// END TOTLENGTH ////////////////
#if defined(ATN) && defined(BP)
		temp *= EXP(jelppi);
#endif
#if defined(TOTLENGTH) || defined(BP)
#ifdef NORM
		temp *= local_norm;
#endif
#ifdef SCATTER
		temp *= local_scat;
#endif
#ifdef ATNM
        temp *= d_atten[idx];
#endif
		temp *= global_factor;
#endif
#endif //////////////// END PET ////////////////
    for (uint ii = 0; ii < NSTEPS; ii++) {
#ifndef USEMAD
        const float3 p = v * t + s;
#else
        const float3 p = FMAD3(t, v, s);
#endif
#if defined(ATN) && defined(FP)
		compute_attenuation(dL, p, d_atten, &jelppi, aa);
#endif
#ifdef TOF //////////////// TOF ////////////////
			TOFSum = TOFLoop(DD, dL, TOFCenter, sigma_x, &D, 1e-6f, TOFWeights);
#endif //////////////// END TOF ////////////////
#if defined(FP) //////////////// FORWARD PROJECTION ////////////////
			denominator(ax, p, dL, d_OSEM
#ifdef TOF //////////////// TOF ////////////////
			, dL, TOFSum, DD, sigma_x, &D, TOFWeights
#ifdef LISTMODE
			, TOFid
#endif
#endif //////////////// END TOF ////////////////
#ifdef N_RAYS
			, lor
#endif
			);
    // if (idx == 0) {
    //     printf("ax[0] = %f\n", ax[0]);
    //     printf("dL = %f\n", dL);
    //     printf("d_OSEM = %f\n", tex3D<float>(d_OSEM, p.x, p.y, p.z));
    //     // printf("d_OSEM = %f\n", read_imagef(d_OSEM, samplerForw, (T4)(p, (typeTT)0)).w);
    // }
#endif  //////////////// END FORWARD PROJECTION ////////////////
#if defined(BP) && !defined(CT) //////////////// BACKWARD PROJECTION ////////////////
            const uint local_ind = CUINT_rtz(p.x * CFLOAT(d_N.x)) + CUINT_rtz(p.y * CFLOAT(d_N.y)) * d_N.x + CUINT_rtz(p.z * CFLOAT(d_N.z)) * d_N.x * d_N.y;
            if (local_ind <= 0 || local_ind >= d_N.x * d_N.y * d_N.z) {
                t += tStep;
                if (t > tEnd)
                    break;
                else
                    continue;
            }
#if defined(MASKBP)
            int maskVal = 1;
            if (aa == 0) {
#ifdef CUDA
#ifdef MASKBP3D
		        maskVal = tex3D<unsigned char>(maskBP, p.x, p.y. p.z);
#else
		        maskVal = tex2D<unsigned char>(maskBP, p.x, p.y);
#endif
#else
#ifdef MASKBP3D
		        maskVal = read_imageui(maskBP, sampler_MASK, (int4)(p.x, p.y, p.z, 0)).w;
#else
		        maskVal = read_imageui(maskBP, sampler_MASK, (int2)(p.x, p.y)).w;
#endif
#endif
            }
            if (maskVal > 0)
#endif
			rhs(dL * temp, ax, local_ind, d_output, no_norm, d_Summ
#ifdef TOF
			, dL, sigma_x, &D, DD, TOFSum, TOFWeights
#ifdef LISTMODE
			, TOFid
#endif
#endif
			);
#endif //////////////// END BACKWARD PROJECTION ////////////////
#if defined(TOF)
			D -= (dL * sign(DD));
#endif
#if !defined(TOTLENGTH) && defined(FP) && !defined(CT)
			LL += dL;
#endif
        t += tStep;
        if (t > tEnd)
            break;
    }
#if !defined(TOTLENGTH) && !defined(CT) && defined(FP)
        if (LL == 0.f)
            LL = L;
#if defined(N_RAYS) && !defined(ORTH) //////////////// MULTIRAY ////////////////
        temp = 1.f / (LL * CFLOAT(N_RAYS));
#else //////////////// SINGLERAY ////////////////
		temp = 1.f / LL;
#endif //////////////// END MULTIRAY ////////////////
#ifdef NORM
        temp *= local_norm;
#endif
#ifdef SCATTER
        temp *= local_scat;
#endif
        temp *= global_factor;
#ifdef ATNM
        temp *= d_atten[idx];
#endif
#endif
#if defined(ATN) && defined(FP)
		temp *= EXP(jelppi);
#endif
#if defined(FP) && !defined(N_RAYS) //////////////// FORWARD PROJECTION ////////////////
#if defined(TOF) && defined(LISTMODE)
		size_t to = TOFid;
#else
#ifndef __CUDACC__ 
#pragma unroll NBINS
#endif
		for (size_t to = 0; to < NBINS; to++) {
#endif
			forwardProjectAF(d_output, ax, idx, temp, to);
#ifdef TOF
			idx += m_size;
#endif
#if defined(TOF) && defined(LISTMODE)
#else
		}
#endif
#elif defined(FP) && defined(N_RAYS)
#if defined(TOF) && defined(LISTMODE)
	int to = TOFid;
#else
#ifndef __CUDACC__ 
#pragma unroll NBINS
#endif
	for (int to = 0; to < NBINS; to++)
#endif
		ax[to + NBINS * lor] *= temp;
#endif //////////////// END FORWARD PROJECTION ////////////////
#ifdef N_RAYS //////////////// MULTIRAY ////////////////
		}
	}
#if defined(FP) //////////////// FORWARD PROJECTION ////////////////
#if defined(TOF) && defined(LISTMODE)
		size_t to = TOFid;
#else
#ifndef __CUDACC__ 
#pragma unroll NBINS
#endif
    for (size_t to = 0; to < NBINS; to++) {
#endif
        float apu = 0.f;
#ifndef __CUDACC__ 
#pragma unroll N_RAYS
#endif
        for (size_t kk = 0; kk < N_RAYS; kk++) {
            apu += ax[to + NBINS * kk];
        }
        ax[to] = apu;
#if defined(TOF) && defined(LISTMODE)
#else
	}
#endif
#if defined(TOF) && defined(LISTMODE)
		to = TOFid;
#else
#ifndef __CUDACC__ 
#pragma unroll NBINS
#endif
    for (size_t to = 0; to < NBINS; to++) {
#endif
        forwardProjectAF(d_output, ax, idx, 1.f, to);
#ifdef TOF
        idx += m_size;
#endif
#if defined(TOF) && defined(LISTMODE)
#else
	}
#endif
#endif //////////////// END FORWARD PROJECTION ////////////////
#endif //////////////// END MULTIRAY ////////////////
}
#endif // END FP

/*******************************************************************************************************************************************
* Voxel-based backprojection. A ray is projected from the source through the center of each voxel and onto the detector panel. This is
* repeated for each projection. Linear or nearest neighbor interpolation is used at the detector panel (linear is default). CT-type data
* only!
*
* INPUTS:
* d_size_x = the number of detector elements (rows),
* d_sizey = the number of detector elements (columns),
* d_dPitch = Either a vector of float2 or two floats if PYTHON is defined, the detector size/pitch in both "row" and "column" directions
* maskBP = 2D/3D backward projection mask, i.e. voxels with 0 will be skipped
* T = redundancy weights for offset imaging,
* d_N = image size in x/y/z- dimension, float3 or three floats (if PYTHON is defined),
* d_b = distance from the pixel space to origin (z/x/y-dimension), float3 or three floats (if PYTHON is defined),
* d_d = distance between adjecent voxels in z/x/y-dimension, float3 or three floats (if PYTHON is defined),
* kerroin = precomputed scaling value, see computeProjectorScalingValues.m
* d_forw = forward projection
* angle = Projection angles for FDK reconstruction
* DSC = Detector to center of rotation distance
* d_OSEM = backward projection,
* d_xy/z = detector x/y/z-coordinates,
* d_uv = Direction coordinates for the detector panel,
* d_Summ = buffer for d_Summ (sensitivity image),
* no_norm = If 1, sensitivity image is not computed,
* d_nProjections = Number of projections/sinograms,
* ii = The current volume, for multi-resolution reconstruction, 0 can be used when not using multi-resolution
*
* OUTPUTS:
* d_OSEM = backprojection,
*******************************************************************************************************************************************/

#if defined(BP) && defined(CT)// START BP
KERNEL2
void projectorType4Backward(const uint d_size_x, const uint d_sizey, 
#ifdef PYTHON
	const float d_dPitchX, const float d_dPitchY, 
#else
	const float2 d_dPitch, 
#endif
#ifdef OFFSET
    CONSTANT float* T, 
#endif
#ifdef PYTHON
	const uint d_Nx, const uint d_Ny, const uint d_Nz, const float bx, const float by, const float bz, 
    const float d_dx, const float d_dy, const float d_dz,
#else
	const uint3 d_N, const float3 b, const float3 d_d, 
#endif
    const float kerroin, 
#ifdef USEIMAGES
    IMAGE3D d_forw, 
#else
    const CLGLOBAL float* CLRESTRICT d_forw, 
#endif
#ifdef FDK
    CONSTANT float* angle, const float DSC, 
#endif
    CLGLOBAL float* CLRESTRICT d_OSEM, 
#if defined(USEGLOBAL)
	const CLGLOBAL float* CLRESTRICT d_xyz,
#else
    CONSTANT float* d_xyz, 
#endif
    CONSTANT float* d_uv, CLGLOBAL float* CLRESTRICT d_Summ, 
    const uchar no_norm, 
#ifdef USEIMAGES
#ifdef MASKBP
#ifdef MASKBP3D
    IMAGE3D maskBP,
#else
    IMAGE2D maskBP,
#endif
#endif
#else
#ifdef MASKBP
    const CLGLOBAL uchar* CLRESTRICT maskBP,
#endif
#endif
    const LONG d_nProjections, const int ii) {

    const int3 i = MINT3(GID0, GID1, GID2 * NVOXELS);

#ifdef PYTHON
	const uint3 d_N = make_uint3(d_Nx, d_Ny, d_Nz);
#endif
    if (i.x >= d_N.x || i.y >= d_N.y || i.z >= d_N.z)
        return;
    size_t idx = GID0 + GID1 * d_N.x + GID2 * NVOXELS * d_N.y * d_N.x;
#ifdef MASKBP
    if (ii == 0) {
#ifdef USEIMAGES
#ifdef CUDA
#ifdef MASKBP3D
        const int maskVal = tex3D<unsigned char>(maskBP, i.x, i.y, i.z);
#else
        const int maskVal = tex2D<unsigned char>(maskBP, i.x, i.y);
#endif
#else
#ifdef MASKBP3D
        const int maskVal = read_imageui(maskBP, sampler_MASK, (int4)(i.x, i.y, i.z, 0)).w;
#else
        const int maskVal = read_imageui(maskBP, sampler_MASK, (int2)(i.x, i.y)).w;
#endif
#endif
#else
#ifdef MASKBP3D
		const int maskVal = maskBP[i.x + i.y * d_N.x + i.z * d_N.x * d_N.y];
#else
        const int maskVal = maskBP[i.x + i.y * d_N.x];
#endif
#endif
        if (maskVal == 0)
            return;
    }
#endif
#ifdef PYTHON
	const float2 d_dPitch = make_float2(d_dPitchX, d_dPitchY);
	const float3 d_d = make_float3(d_dx, d_dy, d_dz);
	const float3 b = make_float3(bx, by, bz);
#endif
    float temp[NVOXELS];
    float wSum[NVOXELS];
    for (int zz = 0; zz < NVOXELS; zz++) {
        temp[zz] = 0.f;
        if (no_norm == 0u)
            wSum[zz] = 0.f;
    }
    const float3 dV = CFLOAT3(i) * d_d + d_d / 2.f + b;
    const float2 koko = MFLOAT2(CFLOAT(d_size_x) * d_dPitch.x, CFLOAT(d_sizey) * d_dPitch.y );
    const float2 indeksi = MFLOAT2(CFLOAT(d_size_x) / 2.f, CFLOAT(d_sizey) / 2.f );
    for (int kk = 0; kk < d_nProjections; kk++) {
        float3 d1, d2, d3;
        float3 s;
#ifndef PARALLEL
        s = CMFLOAT3(d_xyz[kk * 6], d_xyz[kk * 6 + 1], d_xyz[kk * 6 + 2]);
#endif
        d1 = CMFLOAT3(d_xyz[kk * 6 + 3], d_xyz[kk * 6 + 4], d_xyz[kk * 6 + 5]);
#if defined(PITCH)
        const float3 apuX = CMFLOAT3(d_uv[kk * NA], d_uv[kk * NA + 1], d_uv[kk * NA + 2]) * indeksi.x;
        const float3 apuY = CMFLOAT3(d_uv[kk * NA + 3], d_uv[kk * NA + 4], d_uv[kk * NA + 5]) * indeksi.y;
#else
        const float3 apuX = MFLOAT3(indeksi.x * d_uv[kk * NA], indeksi.x * d_uv[kk * NA + 1], 0.f);
        const float3 apuY = MFLOAT3(0.f, 0.f, indeksi.y * d_dPitch.y);
#endif
        d2 = apuX - apuY;
        d3 = d1 - apuX - apuY;
        const float3 normX = normalize(apuX) / koko.x;
        const float3 normY = normalize(apuY) / koko.y;
        const float3 cP = cross(d2, d3 - d1);
        const float pz = (CFLOAT(kk) + 0.5f) / CFLOAT(d_nProjections);
        const float dApu = d_d.z * cP.z;
#ifdef PARALLEL
        const float apuXP = d_uv[kk * NA];
        const float apuYP = d_uv[kk * NA + 1];
        const float3 ss = CMFLOAT3(d_xyz[kk * 6], d_xyz[kk * 6 + 1], d_xyz[kk * 6 + 2]);
        for (int xx = 0; xx < d_size_x; xx++) {
            for (int yy = 0; yy < d_sizey; yy++) {
                s = ss;
	            const float2 indeksiP = MFLOAT2(CFLOAT(xx) - CFLOAT(d_size_x) / 2.f + .5f, CFLOAT(yy) - CFLOAT(d_sizey) / 2.f + .5f);
                (s).x += indeksiP.x * apuXP;
                (s).y += indeksiP.x * apuYP;
                (s).z += indeksiP.y * d_dPitch.y;
                float3 d4 = d1, d5 = d1;
                (d4).x += (indeksiP.x + .5f) * apuXP;
                (d4).y += (indeksiP.x + .5f) * apuYP;
                (d4).z += (indeksiP.y + .5f) * d_dPitch.y;
                (d5).x += (indeksiP.x - .5f) * apuXP;
                (d5).y += (indeksiP.x - .5f) * apuYP;
                (d5).z += (indeksiP.y - .5f) * d_dPitch.y;
                const float3 d44 = d4 - d3;
                const float3 d55 = d5 - d3;
                const float d4x = dot(d44, normX);
                const float d4y = dot(d44, normY);
                const float d5x = dot(d55, normX);
                const float d5y = dot(d55, normY);
                float2 scale = MFLOAT2(1.f / CFLOAT(d_size_x), 1.f / CFLOAT(d_sizey));
                scale /= 2.f;
#endif
        const float upperPart = dot(cP, s - d1);
        float3 v = dV - s;
        float lowerPart = -dot(v, cP);
#ifndef USEMAD
        const float vApu = v.x * v.x + v.y * v.y;
#else
        const float vApu = FMAD(v.x, v.x, v.y * v.y);
#endif
#ifndef __CUDACC__ 
#pragma unroll NVOXELS
#endif
        for (int zz = 0; zz < NVOXELS; zz++) {
            const uint ind = i.z + zz;
            if (ind >= d_N.z)
                break;
            const float t = DIVIDE(upperPart, lowerPart);
#ifndef USEMAD
            float3 p = s + v * t;
            const float l1 = vApu + v.z * v.z;
#else
            float3 p = FMAD3(t, v, s);
            const float l1 = FMAD(v.z, v.z, vApu);
#endif
// #ifdef PARALLEL
//             if (p.x > d4.x || p.x < d5.x || p.y > d4.y || p.y < d5.y || p.z > d4.z || p.z < d5.z) {
//                 v.z += d_d.z;
//                 lowerPart -= dApu;
//                 continue;
//             }
// #endif
#ifdef FDK
            float weight = (DSC + dV.x * COSF(angle[kk]) - dV.y * SINF(angle[kk]));
            weight = (DSC * DSC) / (weight * weight) * (M_PI_F / (CFLOAT(d_nProjections) * d_dPitch.x));
#else
            const float L = distance(p, s);
            const float weight = (L * L * L) / (l1)*kerroin;
#endif
            p -= d3;
            float px = dot(p, normX);
            float py = dot(p, normY);
#ifdef PARALLEL
            if (px > d4x + scale.x || px < d5x - scale.x || py > d4y + scale.y || py < d5y - scale.y) {
                v.z += d_d.z;
                lowerPart -= dApu;
                continue;
            }
#endif
            float yVar = 0.f;
#ifdef USEIMAGES
#ifdef CUDA
            if (px <= 1.f && py <= 1.f && pz <= 1.f && px >= 0.f && py >= 0.f && pz >= 0.f)
                yVar = tex3D<float>(d_forw, px, py, pz);
#else
            if (px <= 1.f && py <= 1.f && pz <= 1.f && px >= 0.f && py >= 0.f && pz >= 0.f)
                yVar = read_imagef(d_forw, samplerIm, CFLOAT4(px, py, pz, 0.f)).w;
#endif
#else
            if (px < 1.f && py < 1.f && pz < 1.f && px >= 0.f && py >= 0.f && pz >= 0.f) {
                const LONG indX = CLONG_rtz(px * CFLOAT(d_size_x));
                const LONG indY = CLONG_rtz(py * CFLOAT(d_sizey)) * CLONG_rtz(d_size_x);
                const LONG indZ = CLONG_rtz(pz * CFLOAT(d_nProjections)) * CLONG_rtz(d_sizey) * CLONG_rtz(d_size_x);
                yVar = d_forw[indX + indY + indZ];
            }
#endif
#ifdef OFFSET
            float TT;
            float Tloc = T[kk];
            if (Tloc > koko.x / 2.f) {
                px = fabs(px - 1.f);
                TT = koko.x - Tloc;
            }
            else
                TT = Tloc;
            px *= koko.x;
            px -= TT;
#endif
            if (yVar != 0.f) {
#ifdef OFFSET
                if (px <= TT && px >= -TT) {
                    float w = .5f * (1.f + SINF(M_PI_F * px / (TT * 2.f)));
                    temp[zz] += w * yVar * weight;
                    if (no_norm == 0u)
                        wSum[zz] += w * weight;
                }
                else if (px < -TT) {
                }
                else {
#endif
                    temp[zz] += yVar * weight;
                    if (no_norm == 0u)
                        wSum[zz] += weight;
#ifdef OFFSET
                }
#endif
            }
            v.z += d_d.z;
            lowerPart -= dApu;
        }
#ifdef PARALLEL
        }
    }
#endif
    }
    for (int zz = 0; zz < NVOXELS; zz++) {
        const uint ind = i.z + zz;
        if (ind >= d_N.z)
            break;
            d_OSEM[idx] += temp[zz];

            if (no_norm == 0u)
                d_Summ[idx] = wSum[zz];
        idx += d_N.y * d_N.x;
    }
}
#endif // END BP