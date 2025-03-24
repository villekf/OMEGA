
/*******************************************************************************************************************************************
* Matrix free projectors for forward and backward projections. This function calculates the sensitivity image d_Summ = sum(A,1) (sum of every 
* row) and A*x and A^Ty (forward and backward projections), where A is the system matrix, y data in the measurement space and x data in the 
* image space.
*
* Used by implementations 2, 3 and 5.
*
* This file contains three different projectors (Siddon, orthogonal, volume-based). Preprocessor commands are used to separate
* different areas of code for the different projectors. Furthermore the forward-backward projection example uses this same file. 
* 64-bit atomics are also currently included in the same file and used if supported. Also nowadays includes the multi-ray case as well for 
* Siddon.
*
* Compiler preprocessing is utilized heavily, for example all the corrections are implemented as compiler preprocesses. The code for 
* specific correction is thus only applied if it has been selected. The kernels are always compiled on-the-fly, though when using same input 
* parameters the kernel should be loaded from cache leading to a slightly faster startup time.
*
* INPUTS:
* global_factor = a global correction factor, e.g. dead time, can be 1.f
* d_epps = a small constant to prevent division by zero, 
* d_size_x = the number of detector elements (rows),
* d_det_per_ring = number of detectors per ring, (only for listmode data sensitivity image computation, can be any value otherwise)
* sigma_x = TOF STD, can be any value otherwise
* crystalSize = Either a vector of float2 or two floats if PYTHON is defined, the detector size/pitch in both "row" and "column" directions
* orthWidth = the width of of the tube used for orthogonal distance based projector (3D), or the radius of the tube with volume projector
* bmin = smaller orthogonal distances than this are fully inside the TOR, volume projector only,
* bmax = Distances greater than this do not touch the TOR, volume projector only,
* Vmax = Full volume of the spherical "voxel", volume projector only,
* maskFP = 2D/3D Forward projection mask, i.e. LORs/measurements with 0 will be skipped
* maskBP = 2D/3D backward projection mask, i.e. voxels with 0 will be skipped
* d_TOFCenter = Offset of the TOF center from the first center of the FOV,
* x/y/z_center = Cartesian coordinates for the center of the voxels (x/y/z-axis),
* V = Precomputed volumes for specific orthogonal distances, volume projector only
* d_sizey = the number of detector elements (columns),
* d_atten = attenuation data (images, if USEIMAGES is defined, or buffer),
* d_nProjections = Number of projections/sinograms,
* d_xy/z = detector x/y/z-coordinates,
* rings = Number of detector rings, PET only and only when computing listmode sensitivity image,
* d_norm = normalization coefficients,
* d_scat = scatter coefficients when using the system matrix method (multiplication), 
* d_Summ = buffer for d_Summ (sensitivity image),
* d_Nxyz = image size in x/y/z- dimension, float3 or three floats (if PYTHON is defined),
* d_d = distance between adjecent voxels in z/x/y-dimension, float3 or three floats (if PYTHON is defined),
* d_b = distance from the pixel space to origin (z/x/y-dimension), float3 or three floats (if PYTHON is defined),
* d_bmax = part in parenthesis of equation (9) in [1] precalculated when k = Nz, float3 or three floats (if PYTHON is defined),
* d_xy/zindex = for sinogram format they determine the detector indices corresponding to each sinogram bin (subset type 3/6/7 only),
* d_L = detector numbers for raw data (unused for sinogram/projection or listmode format),
* d_OSEM = image for current estimates or the input buffer for backward projection,
* d_output = forward or backward projection,
* no_norm = If 1, sensitivity image is not computed,
* m_size = Total number of LORs/measurements for this subset,
* currentSubset = current subset number, zero-based
*
* OUTPUTS:
* d_output = forward or backward projection,
* d_Summ = Sensitivity image (optional, though must be present)
*
* [1] Jacobs, F., Sundermann, E., De Sutter, B., Christiaens, M. Lemahieu, I. (1998). A Fast Algorithm to Calculate the Exact Radiological 
* Path through a Pixel or Voxel Space. Journal of computing and information technology, 6 (1), 89-94.
*
* Copyright (C) 2019-2024 Ville-Veikko Wettenhovi, Niilo Saarlemo
*
* This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it wiL be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************************************************************************/

#ifdef OPENCL
#if (defined(CT) || defined(PET) || defined(SPECT)) && !defined(LISTMODE)
__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, 1)))
#else
__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, 1, 1)))
#endif
#else
extern "C" __global__
#endif
void projectorType123(const float global_factor, const float d_epps, const uint d_size_x, const uint d_det_per_ring,
	const float sigma_x,
#if defined(SPECT)
	const CLGLOBAL float* CLRESTRICT d_rayShiftsDetector,
	const CLGLOBAL float* CLRESTRICT d_rayShiftsSource,
#endif
#ifdef PYTHON
	const float crystalSizeX, const float crystalSizeY, 
#else
	const float2 crystalSize, 
#endif
	///////////////////////// ORTHOGONAL-BASED RAY TRACER /////////////////////////
#ifdef ORTH
	const float orthWidth, const float bmin, const float bmax, const float Vmax,	
#endif
	///////////////////////// END ORTHOGONAL-BASED RAY TRACER /////////////////////////
	///////////////////////// TOF BINS /////////////////////////
#ifdef TOF
	CONSTANT float* TOFCenter, 
#endif
	///////////////////////// END TOF BINS /////////////////////////
	///////////////////////// ORTHOGONAL-BASED RAY TRACER /////////////////////////
#ifdef ORTH
	CONSTANT float* V, 
#endif
	///////////////////////// END ORTHOGONAL-BASED RAY TRACER /////////////////////////
	const uint d_sizey, 
	///////////////////////// PET ATTENUATION CORRECTION /////////////////////////
#if !defined(CT) && defined(ATN) && !defined(ATNM)
#ifdef USEIMAGES
	IMAGE3D d_atten,
#else
	const CLGLOBAL float* CLRESTRICT d_atten,
#endif
#elif !defined(CT) && !defined(ATN) && defined(ATNM)
	const CLGLOBAL float* CLRESTRICT d_atten,
#endif
	///////////////////////// END PET ATTENUATION CORRECTION /////////////////////////
	///////////////////////// FORWARD PROJECTION MASK /////////////////////////
#ifdef USEIMAGES
#ifdef MASKFP
#ifdef MASKFP3D
    IMAGE3D maskFP,
#else
    IMAGE2D maskFP,
#endif
#endif
#if defined(MASKBP) && defined(BP) && !defined(FP)
#ifdef MASKBP3D
    IMAGE3D maskBP,
#else
    IMAGE2D maskBP,
#endif
#endif
#else
#ifdef MASKFP
	const CLGLOBAL uchar* CLRESTRICT maskFP,
#endif
#if defined(MASKBP) && defined(BP) && !defined(FP)
    const CLGLOBAL uchar* CLRESTRICT maskBP,
#endif
#endif
	///////////////////////// END FORWARD PROJECTION MASK /////////////////////////
	///////////////////////// FULL PROJECTIONS/SINOGRAMS /////////////////////////
#if (defined(CT) || defined(SPECT) || defined(PET)) && !defined(LISTMODE)
	const LONG d_nProjections,
#endif
	///////////////////////// END FULL PROJECTIONS/SINOGRAMS /////////////////////////
#if !defined(USEGLOBAL)
	CONSTANT float* d_xy,
#else
	const CLGLOBAL float* CLRESTRICT d_xy,
#endif
	///////////////////////// LISTMODE DATA /////////////////////////
#if (defined(LISTMODE) && !defined(SENS) && !defined(INDEXBASED))
	const CLGLOBAL float* CLRESTRICT d_z,
#else
	CONSTANT float* d_z,
#endif
#ifdef SENS
	const uint rings,
#endif
	///////////////////////// END LISTMODE DATA /////////////////////////
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
	CLGLOBAL CAST* CLRESTRICT d_Summ, 
#ifdef PYTHON
	const uint d_Nx, const uint d_Ny, const uint d_Nz, const float d_dx, const float d_dy, const float d_dz,
	const float bx, const float by, const float bz, const float d_bmaxx, const float d_bmaxy, const float d_bmaxz, 
#else
	const uint3 d_Nxyz, const float3 d_d, const float3 b, const float3 d_bmax,
#endif
	///////////////////////// SUBSET TYPES 3, 5 and 6 /////////////////////////
#if defined(SUBSETS) && !defined(LISTMODE)
	const CLGLOBAL uint* CLRESTRICT d_xyindex, const CLGLOBAL ushort* CLRESTRICT d_zindex, 
#endif
	///////////////////////// END SUBSET TYPES 3, 5 and 6 /////////////////////////
#if defined(INDEXBASED) && defined(LISTMODE) && !defined(SENS)
	const CLGLOBAL ushort* CLRESTRICT trIndex, const CLGLOBAL ushort* CLRESTRICT axIndex,
#endif
#if defined(LISTMODE) && defined(TOF)
	const CLGLOBAL uchar* CLRESTRICT TOFIndex, 
#endif
	///////////////////////// RAW PET DATA /////////////////////////
#ifdef RAW
	const CLGLOBAL ushort* CLRESTRICT d_L, 
#endif
	///////////////////////// END RAW PET DATA /////////////////////////
	///////////////////////// FORWARD OR BACKWARD PROJECTIONS /////////////////////////
#if defined(FP) //&& defined(USEIMAGES)
#ifdef USEIMAGES
	IMAGE3D d_OSEM, 
#else
	const CLGLOBAL float* CLRESTRICT d_OSEM, 
#endif
#else
	const CLGLOBAL float* CLRESTRICT d_OSEM, 
#endif
	///////////////////////// END FORWARD OR BACKWARD PROJECTIONS /////////////////////////
	/* Not yet implemented //////////////
#ifdef MATRIX
	CLGLOBAL int* CLRESTRICT indices, CLGLOBAL int* CLRESTRICT rowInd, CLGLOBAL float* CLRESTRICT values, const int maxLOR, 
#endif
	////////////// Not yet implemented */
	///////////////////////// FORWARD OR BACKWARD PROJECTIONS /////////////////////////
#if defined(BP)
	CLGLOBAL CAST* d_output,
#else
	CLGLOBAL float* d_output,
#endif
	///////////////////////// END FORWARD OR BACKWARD PROJECTIONS /////////////////////////
	const uchar no_norm, const ULONG m_size, const uint currentSubset, const int aa) {
	///*
	// Get the current global index
	int3 i = MINT3(GID0, GID1, GID2);
#if STYPE == 1 || STYPE == 2 || STYPE == 4 || STYPE == 5
	getIndex(&i, d_size_x, d_sizey, currentSubset);
#endif
#if (defined(CT) || defined(SPECT) || defined(PET)) && !defined(LISTMODE)
	size_t idx = i.x + i.y * d_size_x + i.z * d_sizey * d_size_x;
	if (i.x >= d_size_x || i.y >= d_sizey || i.z >= d_nProjections)
#elif defined(SENS)
	int2 indz;
	indz.x = i.z / rings;
	indz.y = i.z % rings;
	size_t idx = i.x + i.y * d_det_per_ring + i.z * d_det_per_ring * d_det_per_ring;
#if defined(NLAYERS)
	if (i.z >= rings * rings) {
		// return;
		i.z -= rings * rings;
		i.x += d_det_per_ring;
		i.y += d_det_per_ring;
		indz.x = i.z / rings + rings;
		indz.y = i.z % rings + rings;
	}
	if (i.x >= d_det_per_ring * NLAYERS || i.y >= d_det_per_ring * NLAYERS || i.z >= rings * rings * NLAYERS)
#else
	if (i.x >= d_det_per_ring || i.y >= d_det_per_ring || i.z >= rings * rings)
#endif
#else
	size_t idx = GID0;
	if (idx >= m_size)
#endif
		return;
#ifdef PYTHON
	const float2 crystalSize = make_float2(crystalSizeX, crystalSizeY);
	const uint3 d_Nxyz = make_uint3(d_Nx, d_Ny, d_Nz);
	const float3 d_d = make_float3(d_dx, d_dy, d_dz);
	const float3 b = make_float3(bx, by, bz);
	const float3 d_bmax = make_float3(d_bmaxx, d_bmaxy, d_bmaxz);
#endif
#ifdef MASKFP // Mask image
#ifdef USEIMAGES
#ifdef CUDA
#ifdef MASKFP3D
	const int maskVal = tex3D<unsigned char>(maskFP, i.x, i.y. i.z);
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
#else
#ifdef MASKFP3D
	const int maskVal = maskFP[i.x + i.y * d_size_x + i.z * d_size_x * d_sizey];
#else
	const int maskVal = maskFP[i.x + i.y * d_size_x];
#endif
#endif
	if (maskVal == 0)
		return;
#endif
#if defined(LISTMODE) && defined(TOF)
	const int TOFid = TOFIndex[idx];
#endif
#if defined(N_RAYS) && defined(FP)
	float ax[NBINS * N_RAYS];
#else
	float ax[NBINS];
#endif
#ifdef BP
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
	const uint d_Nxy = d_Nxyz.x * d_Nxyz.y;
	float local_norm = 0.f;
	float local_scat = 0.f;

#ifdef NORM // Normalization included
	local_norm = d_norm[idx];
#endif
#ifdef SCATTER // Scatter data included
	local_scat = d_scat[idx];
#endif
#ifdef ORTH // Orthogonal or volume-based
	float b1, b2, d1, d2;
	float _bz = b.z, dz = d_d.z;
#endif

#if defined(N_RAYS) //////////////// MULTIRAY ////////////////
	int lor = -1;
	// Load the next detector index
	// raw data
	for (int lorZ = 0u; lorZ < N_RAYS3D; lorZ++) {
		for (int lorXY = 0u; lorXY < N_RAYS2D; lorXY++) {
			lor++;
#elif defined(SPECT)
	const int lorXY = 0;
#endif  //////////////// END MULTIRAY ////////////////
	float3 s, d;
#if defined(NLAYERS) && !defined(LISTMODE)
	const uint layer = i.z / NLAYERS;
	// if (layer != 0)
	// if (layer == 1 || layer == 2)
		// return;
#endif
	// Load the next detector index
#if defined(CT) && !defined(LISTMODE) && !defined(PET) // CT data
	getDetectorCoordinatesCT(d_xy, d_z, &s, &d, i, d_size_x, d_sizey, crystalSize);
#elif defined(SPECT) && !defined(LISTMODE) && !defined(PET) // SPECT data
	getDetectorCoordinatesSPECT(d_xy, d_z, &s, &d, i, d_size_x, d_sizey, crystalSize, d_rayShiftsDetector, d_rayShiftsSource, lorXY);
#elif defined(LISTMODE) && !defined(SENS) // Listmode data
#if defined(INDEXBASED)
	getDetectorCoordinatesListmode(d_xy, d_z, trIndex, axIndex, &s, &d, idx
#else
	getDetectorCoordinatesListmode(d_xy, &s, &d, idx
#endif
#if defined(N_RAYS)
		, lorXY, lorZ, crystalSize
#endif
	);
#elif defined(RAW) || defined(SENS) // raw data
	getDetectorCoordinatesRaw(d_xy, d_z, i, &s, &d, indz
#if defined(N_RAYS)
		, lorXY, lorZ, crystalSize
#endif
	);
#elif !defined(SUBSETS) // Subset types 1, 2, 4, 5, 8, 9, 10, 11
	getDetectorCoordinatesFullSinogram(d_size_x, i, &s, &d, d_xy, d_z
#if defined(N_RAYS)
		, lorXY, lorZ, crystalSize
#endif
#if defined(NLAYERS)
		, d_sizey, layer
#endif
	);
#else // Subset types 3, 6, 7
	getDetectorCoordinates(d_xyindex, d_zindex, idx, &s, &d, d_xy, d_z
#if defined(N_RAYS)
		, lorXY, lorZ, crystalSize
#endif
#if defined(NLAYERS)
		, d_sizey, d_size_x
#endif
	);
#endif
	// Calculate the x, y and z distances of the detector pair
	float3 diff = d - s;

#ifdef CUDA
	if ((diff.x == 0.f && diff.y == 0.f && diff.z == 0.f) || (diff.x == 0.f && diff.y == 0.f) || isinf(diff.x) || isinf(diff.y) || isinf(diff.z) || isnan(diff.x) || isnan(diff.y) || isnan(diff.z))
#else
	if (all(diff == 0.f) || (diff.x == 0.f && diff.y == 0.f) || any(isinf(diff)) || any(isnan(diff)))
#endif
#ifdef N_RAYS //////////////// MULTIRAY ////////////////
		continue;
#else
		return;
#endif  //////////////// END MULTIRAY ////////////////
	uint Np = 0u;

#ifdef ATN // Attenuation included
	float jelppi = 0.f;
#endif
	float temp = 1.f;
	uint d_N0 = d_Nxyz.x;
	uint d_N1 = d_Nxyz.y;
	uint d_N2 = 1u;
	uint d_N3 = d_Nxyz.x;

	int tempi = 0, tempj = 0, tempk = 0, ux = 0, uy = 0, uz = 0;

	float L = length(diff);
#ifndef TOTLENGTH
	float LL = 0.f;
#endif
	LONG local_ind = 0u;
	int3 localInd = MINT3(0, 0, 0);
#ifdef TOF //////////////// TOF ////////////////
	float D = 0.f;
	float DD = 0.f;
#endif //////////////// END TOF ////////////////
	float local_ele = 0.f;

#ifdef ORTH //////////////// ORTHOGONAL OR VOLUME-BASED RAY TRACER ////////////////
	bool XY = false;
	float kerroin = 0.f;
#if !defined(VOL) // Orthogonal
	kerroin = L * orthWidth;
#elif defined(VOL) // Volume-based
	kerroin = L;
#if defined(TOTLENGTH)
	float TotV = L * M_1_PI_F * orthWidth * orthWidth;
#else
	float TotV = M_1_PI_F * orthWidth * orthWidth;
#endif
#endif
#endif //////////////// END ORTHOGONAL OR VOLUME-BASED RAY TRACER OR SIDDON ////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//If the LOR is perpendicular in the y-direction (Siddon cannot be used)
	if (fabs(diff.z) < 1e-6f && (fabs(diff.y) < 1e-6f || fabs(diff.x) < 1e-6f)) {

		//return;
		tempk = CINT(fabs(s.z - b.z) / d_d.z);
		if (tempk < 0 || tempk >= d_Nxyz.z)
#ifdef N_RAYS //////////////// MULTIRAY ////////////////
			continue;
#else
			return;
#endif  //////////////// END MULTIRAY ////////////////

#if defined(ORTH) //////////////// ORTHOGONAL OR VOLUME-BASED RAY TRACER ////////////////
		int indO = 0;
#endif //////////////// END ORTHOGONAL OR VOLUME-BASED RAY TRACER ////////////////
		float d_b, dd, d_db, d_d2;
		int apuX1, apuX2;
// #if defined(LISTMODE)
		float dT1, dT2;
// #endif
		if (fabs(diff.y) < 1e-6f && d.y <= d_bmax.y && d.y >= b.y && s.y <= d_bmax.y && s.y >= b.y) {
			apuX1 = 0;
			apuX2 = d_Nxyz.x - 1;
// #if defined(LISTMODE)
			float dist1, dist2 = 0.f;
			if (s.x > d.x) {
				dist1 = (b.x - d.x);
				dist2 = (b.x + CFLOAT(d_Nxyz.x) * d_d.x - s.x);
			}
			else {
				dist1 = (b.x - s.x);
				dist2 = (b.x + CFLOAT(d_Nxyz.x) * d_d.x - d.x);
			}
			for (int kk = 0; kk < d_Nxyz.x; kk++) {
				if (dist1 >= 0.f) {
					apuX1 = kk;
					if (kk == 0)
						dT1 = d_d.x;
					else
						dT1 = min(dist1, d_d.x);
					break;
				}
				dist1 += d_d.x;
			}
			for (int kk = d_Nxyz.x - 1; kk >= apuX1; kk--) {
				if (dist2 <= 0.f) {
					apuX2 = kk;
					if (kk == d_Nxyz.x - 1)
						dT2 = d_d.x;
					else
						dT2 = min(-dist2, d_d.x);
					break;
				}
				dist2 -= d_d.x;
			}
// #endif
			d_b = b.y;
			dd = d.y;
			d_db = d_d.y;
			d_d2 = d_d.x;
			d_N0 = d_Nxyz.y;
			d_N1 = d_Nxyz.x;
			d_N2 = d_Nxyz.y;
			d_N3 = 1u;
#if defined(ORTH) //////////////// ORTHOGONAL OR VOLUME-BASED RAY TRACER ////////////////
			b1 = b.x;
			b2 = b.y;
			d1 = d_d.x;
			d2 = d_d.y;
			XY = true;
#endif //////////////// END ORTHOGONAL OR VOLUME-BASED RAY TRACER OR SIDDON ////////////////
		}
		else if (fabs(diff.x) < 1e-6f && d.x <= d_bmax.x && d.x >= b.x && s.x <= d_bmax.x && s.x >= b.x) {
			apuX1 = 0;
			apuX2 = d_Nxyz.y - 1;
// #if defined(LISTMODE)
			float dist1, dist2 = 0.f;
			if (s.y > d.y) {
				dist1 = (b.y - d.y);
				dist2 = (b.y + CFLOAT(d_Nxyz.y) * d_d.y - s.y);
			}
			else {
				dist1 = (b.y - s.y);
				dist2 = (b.y + CFLOAT(d_Nxyz.y) * d_d.y - d.y);
			}
			for (int kk = 0; kk < d_Nxyz.y; kk++) {
				if (dist1 >= 0.f) {
					apuX1 = kk;
					if (kk == 0)
						dT1 = d_d.y;
					else
						dT1 = min(dist1, d_d.y);
					break;
				}
				dist1 += d_d.y;
			}
			for (int kk = d_Nxyz.y - 1; kk >= apuX1; kk--) {
				if (dist2 <= 0.f) {
					apuX2 = kk;
					if (kk == d_Nxyz.y - 1)
						dT2 = d_d.y;
					else
						dT2 = min(-dist2, d_d.y);
					break;
				}
				dist2 -= d_d.y;
			}
// #endif
			d_b = b.x;
			dd = d.x;
			d_db = d_d.x;
			d_d2 = d_d.y;
#if defined(ORTH) //////////////// ORTHOGONAL OR VOLUME-BASED RAY TRACER ////////////////
			b1 = b.y;
			b2 = b.x;
			d1 = d_d.y;
			d2 = d_d.x;
			float s_b = s.x;
			s.x = s.y;
			s.y = s_b;
			float diff_b = diff.x;
			diff.x = diff.y;
			diff.y = diff_b;
#endif //////////////// END ORTHOGONAL OR VOLUME-BASED RAY TRACER OR SIDDON ////////////////

		}
		else
#ifdef N_RAYS //////////////// MULTIRAY ////////////////
		    continue;
#else
		    return;
#endif  //////////////// END MULTIRAY ////////////////
		localInd.z = tempk;
		local_ind = CLONG_rtz(tempk);
		perpendicular_elements(d_b, d_db, d_N0, dd, d_d2, d_N1, &temp, &localInd, &local_ind, d_N2, d_N3, idx, global_factor, local_scat, 
#if !defined(CT) && defined(ATN) && !defined(ATNM)
			d_atten, aa, 
#elif !defined(CT) && !defined(ATN) && defined(ATNM)
			d_atten, 
#endif
		    local_norm, L);
#if defined(ORTH)
#if defined(VOL)
#if defined(TOTLENGTH)
		temp *= (1.f / TotV);
#else
		temp *= (1.f / (TotV * d_d2 * CFLOAT(apuX2 - apuX1)));
#endif
#endif
		if (d_N2 == 1)
			indO = localInd.x;
		else
			indO = localInd.y;
#endif
#ifdef TOF //////////////// TOF ////////////////
			float dI = (d_d2 * d_N1) / 2.f * sign(diff.y);
			D = dI;
			DD = D;
#endif //////////////// END TOF ////////////////
			for (uint ii = apuX1; ii <= apuX2; ii++) {
				float d_in = d_d2;
// #if defined(LISTMODE)
				if (ii == apuX1) {
					local_ind += CLONG_rtz(d_N3 * ii);
					if (d_N3 == 1)
						localInd.x += ii;
					else
						localInd.y += ii;
				}
				if (apuX1 > 0 && ii == apuX1)
					d_in = dT1;
				else if (apuX2 < d_N1 - 1 && ii == apuX2)
					d_in = dT2;
// #endif
#if defined(ATN) && defined(SPECT)
#ifdef USEIMAGES
				compute_attenuation(d_in, localInd, d_atten, &jelppi, aa);
#else
				compute_attenuation(d_in, local_ind, d_atten, &jelppi, aa);
#endif
#if defined(SPECT) && !defined(ORTH)
				d_in *= EXP(jelppi);
#endif
#endif
#ifdef TOF //////////////// TOF ////////////////
				const float TOFSum = TOFLoop(DD, d_in, TOFCenter, sigma_x, &D, d_epps);
#endif //////////////// END TOF ////////////////
#ifdef ORTH //////////////// ORTH/VOL ////////////////
				const float xcenter = b1 + d1 * CFLOAT(ii) + d1 / 2.f;
				orthDistance3D(ii, diff.y, diff.x, diff.z, xcenter, b2, d2, _bz, dz, temp, indO, localInd.z, s.x, s.y, s.z, d_Nxy, kerroin, d_N1, d_N3, d_N2, d_Nxyz.z, bmin, bmax, Vmax, V, XY, ax, 
#if defined(FP)
					d_OSEM
#else
					no_norm, d_Summ, d_output 
#endif
#ifdef TOF
				, d_d2, sigma_x, &D, DD, TOFCenter, TOFSum
#if defined(LISTMODE)
				, TOFid
#endif
#endif
#if defined(SPECT) && defined(ATN)
				, EXP(jelppi)
#endif
#if defined(MASKBP) && defined(BP)
				, aa, maskBP
#endif
				);
#else //////////////// SIDDON ////////////////
#if defined(FP) //////////////// FORWARD PROJECTION ////////////////
#ifdef USEIMAGES
				denominator(ax, localInd, d_in, d_OSEM
#else
				denominator(ax, local_ind, d_in, d_OSEM
#endif
#ifdef TOF //////////////// TOF ////////////////
				, d_in, TOFSum, DD, TOFCenter, sigma_x, &D
#ifdef LISTMODE
				, TOFid
#endif
#endif //////////////// END TOF ////////////////
#ifdef N_RAYS
				, lor
#endif
				);
#endif  //////////////// END FORWARD PROJECTION ////////////////
#if defined(BP) //////////////// BACKWARD PROJECTION ////////////////
#if defined(MASKBP) //////////////// MASKBP ////////////////
				int maskVal = 1;
				if (aa == 0) {
#ifdef USEIMAGES
#ifdef CUDA
#ifdef MASKBP3D
					maskVal = tex3D<unsigned char>(maskBP, localInd.x, localInd.y, localInd.z);
#else
					maskVal = tex2D<unsigned char>(maskBP, localInd.x, localInd.y);
#endif
#else
#ifdef MASKBP3D
					maskVal = read_imageui(maskBP, sampler_MASK, (int4)(localInd.x, localInd.y,  localInd.z, 0)).w;
#else
					maskVal = read_imageui(maskBP, sampler_MASK, (int2)(localInd.x, localInd.y)).w;
#endif
#endif
#else
#ifdef MASKBP3D
					maskVal = maskBP[localInd.x + localInd.y * d_Nxyz.x + localInd.z * d_Nxyz.x * d_Nxyz.y];
#else
					maskVal = maskBP[localInd.x + localInd.y * d_Nxyz.x];
#endif
#endif
				}
				if (maskVal > 0)
#endif //////////////// END MASKBP ////////////////
				rhs(temp * d_in, ax, local_ind, d_output, no_norm, d_Summ
#ifdef TOF //////////////// TOF ////////////////
				, d_in, sigma_x, &D, DD, TOFCenter, TOFSum
#ifdef LISTMODE
				, TOFid
#endif
#endif //////////////// END TOF ////////////////
				);
#endif //////////////// END BACKWARD PROJECTION ////////////////
#endif //////////////// END SIDDON/ORTH/VOL ////////////////
#if (!defined(USEIMAGES) && defined(FP)) || defined(BP) || (defined(ATN) && defined(SPECT) && !defined(USEIMAGES))
				local_ind += CLONG_rtz(d_N3);
#endif 
#if defined(FP) || (defined(MASKBP) && defined(BP)) || (defined(ATN) && defined(SPECT) && defined(USEIMAGES))
				if (d_N3 == 1)
					localInd.x++;
				else
					localInd.y++;
#endif
#if defined(TOF)
				D -= (d_d2 * sign(DD));
#endif
			}
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
		size_t to = TOFid;
#else
#ifndef __CUDACC__ 
#pragma unroll NBINS
#endif
	for (int to = 0; to < NBINS; to++)
#endif
		ax[to + NBINS * lor] *= temp;
#endif //////////////// END FORWARD PROJECTION ////////////////

	}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// /*
	else {
		// return;
		float txu = 0.f, tyu = 0.f, tzu = 0.f, tc = 0.f, tx0 = 1e8f, ty0 = 1e8f, tz0 = 1e8f;
		bool skip = false, XY = true;
#ifdef TOF //////////////// TOF ////////////////
        float TOFSum = 0.f;
#endif //////////////// END TOF ////////////////

		// If the measurement is on a same ring
			// Z-coordinate (ring)
		if (fabs(diff.z) < 1e-6f) {
		    	// return;
			tempk = CINT(fabs(s.z - b.z) / d_d.z);
			if (tempk < 0 || tempk >= d_Nxyz.z)
#ifdef N_RAYS //////////////// MULTIRAY ////////////////
		    	continue;
#else
		    	return;
#endif  //////////////// END MULTIRAY ////////////////
			skip = siddon_pre_loop_2D(b.x, b.y, diff.x, diff.y, d_bmax.x, d_bmax.y, d_d.x, d_d.y, d_Nxyz.x, d_Nxyz.y, &tempi, &tempj, &txu, &tyu, &Np, TYPE,
				s.y, s.x, d.y, d.x, &tc, &ux, &uy, &tx0, &ty0, &XY);
		}
		else if (fabs(diff.y) < 1e-6f) {
			// return;
			tempj = perpendicular_start(b.y, d.y, d_d.y, d_Nxyz.y);
			skip = siddon_pre_loop_2D(b.x, b.z, diff.x, diff.z, d_bmax.x, d_bmax.z, d_d.x, d_d.z, d_Nxyz.x, d_Nxyz.z, &tempi, &tempk, &txu, &tzu, &Np, TYPE,
				s.z, s.x, d.z, d.x, &tc, &ux, &uz, &tx0, &tz0, &XY);
			XY = true;
			if (d.y > d_bmax.y || d.y < b.y)
				skip = true;
		}
		else if (fabs(diff.x) < 1e-6f) {
			// return;
			tempi = perpendicular_start(b.x, d.x, d_d.x, d_Nxyz.x);
			skip = siddon_pre_loop_2D(b.y, b.z, diff.y, diff.z, d_bmax.y, d_bmax.z, d_d.y, d_d.z, d_Nxyz.y, d_Nxyz.z, &tempj, &tempk, &tyu, &tzu, &Np, TYPE,
				s.z, s.y, d.z, d.y, &tc, &uy, &uz, &ty0, &tz0, &XY);
			XY = false;
			if (d.x > d_bmax.x || d.x < b.x)
				skip = true;
		}
		else {
			skip = siddon_pre_loop_3D(b, diff, d_bmax, d_d, d_Nxyz, &tempi, &tempj, &tempk, &txu, &tyu, &tzu, &Np, TYPE, s, d, &tc, &ux, &uy, &uz, &tx0, &ty0, &tz0, &XY, i);
		}
		if (skip)
#ifdef N_RAYS //////////////// MULTIRAY ////////////////
		    continue;
#else
		    return;
#endif  //////////////// END MULTIRAY ////////////////
		// L = length(s + (d - s) * tc);
#ifdef TOF
		TOFDis(diff, tc, L, &D, &DD);
#endif
#if (defined(ATN) && defined(BP)) || defined(ORTH) || (defined(BP) && !defined(TOTLENGTH) && !defined(CT))
		float tx0_a = tx0, ty0_a = ty0, tz0_a = tz0;
		int tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
#endif
#if ((defined(ATN)) || defined(TOF)) && defined(ORTH)
		float tx0_c = tx0, ty0_c = ty0, tz0_c = tz0, txu_c = txu, tyu_c = tyu, tzu_c = tzu, tc_c = tc;
		int tempi_c = tempi, tempj_c = tempj, tempk_c = tempk, ux_c = ux, uy_c = uy, uz_c = uz;
#endif
#if (defined(ATN) && defined(BP) && !defined(SPECT)) || (defined(BP) && !defined(TOTLENGTH) && !defined(CT))
		float tc_a = tc;
		for (uint ii = 0u; ii < Np; ii++) {
			local_ele = 0.f;
#if defined(ATN) && !defined(SPECT)
#ifdef USEIMAGES
			localInd = CMINT3(tempi, tempj, tempk);
#else
			local_ind = compute_ind(tempj, tempi * d_N2, tempk, d_N3, d_Nxy);
#endif
#endif
			if (tz0 < ty0 && tz0 < tx0) {
// #if defined(LISTMODE)
				if (tz0 >= 0.f && tz0 <= 1.f) {
					if (tc < 0.f) {
						local_ele = tz0 * L;
						compute_element(&tz0, &tc, L, tzu, uz, &tempk);
					}
					else
						local_ele = compute_element(&tz0, &tc, L, tzu, uz, &tempk);
				}
				else if (tc >= 0.f && tc <= 1.f) {
					if (tz0 > 1.f) {
						local_ele = (1.f - tc) * L;
						compute_element(&tz0, &tc, L, tzu, uz, &tempk);
					}
					else
					local_ele = compute_element(&tz0, &tc, L, tzu, uz, &tempk);
				}
				else
					compute_element(&tz0, &tc, L, tzu, uz, &tempk);
// #else
// 				local_ele = compute_element(&tz0, &tc, L, tzu, uz, &tempk);
// #endif
			}
			else if (ty0 < tx0) {
// #if defined(LISTMODE)
				if (ty0 >= 0.f && ty0 <= 1.f) {
					if (tc < 0.f) {
						local_ele = ty0 * L;
						compute_element(&ty0, &tc, L, tyu, uy, &tempj);
					}
					else
						local_ele = compute_element(&ty0, &tc, L, tyu, uy, &tempj);
				}
				else if (tc >= 0.f && tc <= 1.f) {
					if (ty0 > 1.f) {
						local_ele = (1.f - tc) * L;
						compute_element(&ty0, &tc, L, tyu, uy, &tempj);
					}
					else
						local_ele = compute_element(&ty0, &tc, L, tyu, uy, &tempj);
				}
				else
					compute_element(&ty0, &tc, L, tyu, uy, &tempj);
// #else
// 				local_ele = compute_element(&ty0, &tc, L, tyu, uy, &tempj);
// #endif
			}
			else {
// #if defined(LISTMODE)
				if (tx0 >= 0.f && tx0 <= 1.f) {
					if (tc < 0.f) {
						local_ele = tx0 * L;
						compute_element(&tx0, &tc, L, txu, ux, &tempi);
					}
					else
						local_ele = compute_element(&tx0, &tc, L, txu, ux, &tempi);
				}
				else if (tc >= 0.f && tc <= 1.f) {
					if (tx0 > 1.f) {
						local_ele = (1.f - tc) * L;
						compute_element(&tx0, &tc, L, txu, ux, &tempi);
					}
					else
						local_ele = compute_element(&tx0, &tc, L, txu, ux, &tempi);
				}
				else
					compute_element(&tx0, &tc, L, txu, ux, &tempi);
// #else
// 				local_ele = compute_element(&tx0, &tc, L, txu, ux, &tempi);
// #endif
			}
#if defined(ATN) && !defined(SPECT)
#ifdef USEIMAGES
			compute_attenuation(local_ele, localInd, d_atten, &jelppi, aa);
#else
			compute_attenuation(local_ele, local_ind, d_atten, &jelppi, aa);
#endif
#endif
#if !defined(TOTLENGTH)
			LL += local_ele;
#endif
			if (tempi < 0 || tempi >= d_Nxyz.x || tempj < 0 || tempj >= d_Nxyz.y || tempk < 0 || tempk >= d_Nxyz.z) {
				break;
			}
		}
		tempi = tempi_a;
		tempj = tempj_a;
		tempk = tempk_a;
		tx0 = tx0_a;
		ty0 = ty0_a;
		tz0 = tz0_a;
		tc = tc_a;
#endif
#ifdef ORTH // Orthogonal or volume-based
		int tempi_b, u_b;
		float t0_b, tu_b, diff_b, s_b, b1, b2, d1, d2;
		float _bz = b.z, dz = d_d.z;
		uint3 d_NN = MUINT3(d_Nxyz.x, d_Nxyz.y, d_Nxyz.z);
		if (!XY) {
			b1 = b.y;
			b2 = b.x;
			d1 = d_d.y;
			d2 = d_d.x;
			tempi_b = tempi;
			tempi = tempj;
			tempj = tempi_b;
			s_b = s.x;
			s.x = s.y;
			s.y = s_b;
			u_b = ux;
			ux = uy;
			uy = u_b;
			t0_b = tx0;
			tx0 = ty0;
			ty0 = t0_b;
			tu_b = txu;
			txu = tyu;
			tyu = tu_b;
			diff_b = diff.x;
			diff.x = diff.y;
			diff.y = diff_b;
			d_NN.x = d_Nxyz.y;
			d_NN.y = d_Nxyz.x;
			d_N0 = d_Nxyz.y;
			d_N1 = d_Nxyz.x;
			d_N2 = d_Nxyz.x;
			d_N3 = 1;
		}
		else {
			b1 = b.x;
			b2 = b.y;
			d1 = d_d.x;
			d2 = d_d.y;
			d_N0 = d_Nxyz.x;
			d_N1 = d_Nxyz.y;
			d_N2 = 1;
			d_N3 = d_Nxyz.x;
		}
		tx0_a = tx0, ty0_a = ty0, tz0_a = tz0;
		tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
#endif

#if !defined(CT) //////////////// PET/SPECT ////////////////
#ifdef ORTH //////////////// ORTH/VOL ////////////////
#if defined(VOL) && defined(TOTLENGTH) //////////////// VOL+TOTLENGTH ////////////////
		temp = 1.f / TotV;
#elif defined(VOL) && !defined(TOTLENGTH) && defined(BP)
		if (LL == 0.f)
			LL = L;
		temp = 1.f / (TotV * LL);
#endif //////////////// END VOL+TOTLENGTH ////////////////
#else //////////////// SIDDON ////////////////
#if defined(TOTLENGTH) //////////////// TOTLENGTH ////////////////
#ifdef N_RAYS  //////////////// MULTIRAY ////////////////
		temp = 1.f / (L * CFLOAT(N_RAYS));
#else //////////////// SINGLERAY ////////////////
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
#endif //////////////// END ORTH/VOL/SIDDON ////////////////
#if defined(ATN) && defined(BP) && !defined(SPECT)
		temp *= EXP(jelppi);
#endif
#if defined(TOTLENGTH) || defined(BP)
#ifdef NORM //////////////// NORM ////////////////
		temp *= local_norm;
#endif //////////////// END NORM ////////////////
#ifdef SCATTER //////////////// SCATTER ////////////////
		temp *= local_scat;
#endif //////////////// END SCATTER ////////////////
		temp *= global_factor;
#ifdef ATNM //////////////// ATTENUATIONLOR ////////////////
		temp *= d_atten[idx];
#endif //////////////// END ATTENUATIONLOR ////////////////
#endif
#endif //////////////// END PET/SPECT ////////////////
// #if (defined(SPECT) && !defined(ORTH)) // Ray length inside BP mask
// 		float L_SPECT = 0.f;
// //#ifdef N_RAYS
// 		temp = temp * (L * CFLOAT(N_RAYS));
// //#else
// //		temp = temp * L;
// //#endif
// #endif

		for (uint ii = 0u; ii < Np; ii++) {
// #if defined(LISTMODE)
			local_ele = 0.f;
// #endif
			local_ind = compute_ind(tempj, tempi * d_N2, tempk, d_N3, d_Nxy);
			localInd = CMINT3(tempi, tempj, tempk);
#if defined(ATN) && (defined(FP) || defined(SPECT)) //////////////// ATTENUATION ////////////////
#ifdef USEIMAGES
			int3 localInd2 = CMINT3(tempi, tempj, tempk);
#else
			LONG localInd2 = local_ind;
#endif
#endif //////////////// END ATTENUATION ////////////////
#ifdef ORTH //////////////// ORTH/VOL ////////////////
			tx0_a = tx0;
			ty0_a = ty0;
			tz0_a = tz0;
			tempi_a = tempi;
			tempj_a = tempj;
			tempk_a = tempk;
#endif //////////////// END ORTH/VOL ////////////////
			bool pass = false;
			if (tz0 < ty0 && tz0 < tx0) {
// #if defined(LISTMODE)
				if (tz0 >= 0.f && tz0 <= 1.f) {
					if (tc < 0.f) {
						local_ele = tz0 * L;
						compute_element(&tz0, &tc, L, tzu, uz, &tempk);
					}
					else
						local_ele = compute_element(&tz0, &tc, L, tzu, uz, &tempk);
					pass = true;
				}
				else if (tc >= 0.f && tc <= 1.f) {
					if (tz0 > 1.f) {
						local_ele = (1.f - tc) * L;
						compute_element(&tz0, &tc, L, tzu, uz, &tempk);
					}
					else
						local_ele = compute_element(&tz0, &tc, L, tzu, uz, &tempk);
					pass = true;
				}
				else
					compute_element(&tz0, &tc, L, tzu, uz, &tempk);
// #else
// 				local_ele = compute_element(&tz0, &tc, L, tzu, uz, &tempk);
// #endif
			}
			else if (ty0 < tx0) {
// #if defined(LISTMODE)
				if (ty0 >= 0.f && ty0 <= 1.f) {
					if (tc < 0.f) {
						local_ele = ty0 * L;
						compute_element(&ty0, &tc, L, tyu, uy, &tempj);
					}
					else
						local_ele = compute_element(&ty0, &tc, L, tyu, uy, &tempj);
					pass = true;
				}
				else if (tc >= 0.f && tc <= 1.f) {
					if (ty0 > 1.f) {
						local_ele = (1.f - tc) * L;
						compute_element(&ty0, &tc, L, tyu, uy, &tempj);
					}
					else
						local_ele = compute_element(&ty0, &tc, L, tyu, uy, &tempj);
					pass = true;
				}
				else
					compute_element(&ty0, &tc, L, tyu, uy, &tempj);
// #else
// 				local_ele = compute_element(&ty0, &tc, L, tyu, uy, &tempj);
// #endif
			}
			else {
// #if defined(LISTMODE)
				if (tx0 >= 0.f && tx0 <= 1.f) {
					if (tc < 0.f) {
						local_ele = tx0 * L;
						compute_element(&tx0, &tc, L, txu, ux, &tempi);
					}
					else
						local_ele = compute_element(&tx0, &tc, L, txu, ux, &tempi);
					pass = true;
				}
				else if (tc >= 0.f && tc <= 1.f) {
					if (tx0 > 1.f) {
						local_ele = (1.f - tc) * L;
						compute_element(&tx0, &tc, L, txu, ux, &tempi);
					}
					else
						local_ele = compute_element(&tx0, &tc, L, txu, ux, &tempi);
					pass = true;
				}
				else
					compute_element(&tx0, &tc, L, txu, ux, &tempi);
// #else
// 				local_ele = compute_element(&tx0, &tc, L, txu, ux, &tempi);
// #endif
			}
#if !defined(TOTLENGTH) && defined(FP) && !defined(CT)
			LL += local_ele;
#endif
#if (defined(ATN) && (defined(FP) || defined(SPECT))) || defined(TOF)
			float local_ele2 = local_ele;
#endif
#if ((defined(ATN) && defined(FP)) || defined(TOF)) && defined(ORTH)
#if defined(ATN)
#ifdef USEIMAGES
			localInd2 = CMINT3(tempi_c, tempj_c, tempk_c);
#else
			localInd2 = compute_ind(tempj_c, tempi_c, tempk_c, d_Nxyz.x, d_Nxy);
#endif
#endif
			if (tz0_c < ty0_c && tz0_c < tx0_c) {
				local_ele2 = compute_element(&tz0_c, &tc_c, L, tzu_c, uz_c, &tempk_c);
			}
			else if (ty0_c < tx0_c) {
				local_ele2 = compute_element(&ty0_c, &tc_c, L, tyu_c, uy_c, &tempj_c);
			}
			else {
				local_ele2 = compute_element(&tx0_c, &tc_c, L, txu_c, ux_c, &tempi_c);
			}
#endif
#if defined(ATN) && (defined(FP) || defined(SPECT))
			if (pass) {
				compute_attenuation(local_ele2, localInd2, d_atten, &jelppi, aa);
#if defined(SPECT) && !defined(ORTH)
				local_ele *= EXP(jelppi);
#endif
			}
#endif
#ifdef TOF //////////////// TOF ////////////////
			TOFSum = TOFLoop(DD, local_ele2, TOFCenter, sigma_x, &D, d_epps);
#endif //////////////// END TOF ////////////////
#ifdef ORTH //////////////// ORTH/VOL ////////////////
            if (ii == 0) {
                if (ux >= 0) {
                    for (int kk = tempi_a - 1; kk >= 0; kk--) {
						const float xcenter = b1 + d1 * CFLOAT(kk) + d1 / 2.f;
                        int uu = orthDistance3D(kk, diff.y, diff.x, diff.z, xcenter, b2, d2, _bz, dz, temp, tempj_a, tempk_a, s.x, s.y, s.z, d_Nxy, kerroin, d_N1, d_N2, d_N3, d_Nxyz.z, bmin, bmax, Vmax, V, XY, ax, 
#if defined(FP) //////////////// FP ////////////////
                            d_OSEM
#else //////////////// BP ////////////////
                            no_norm, d_Summ, d_output 
#endif //////////////// END FP/BP ////////////////
#ifdef TOF //////////////// TOF ////////////////
                        , local_ele2, sigma_x, &D, DD, TOFCenter, TOFSum
#if defined(LISTMODE)
						, TOFid
#endif
#endif //////////////// END TOF ////////////////
#if defined(SPECT) && defined(ATN)
						, EXP(jelppi)
#endif
#if defined(MASKBP) && defined(BP) //////////////// MASKBP ////////////////
						, aa, maskBP
#endif //////////////// END MASKBP ////////////////
                        );
                        if (uu == 0)
                            break;
                    }
                }
                else {
                    for (int kk = tempi_a + 1; kk < d_NN.x; kk++) {
						const float xcenter = b1 + d1 * CFLOAT(kk) + d1 / 2.f;
                        int uu = orthDistance3D(kk, diff.y, diff.x, diff.z, xcenter, b2, d2, _bz, dz, temp, tempj_a, tempk_a, s.x, s.y, s.z, d_Nxy, kerroin, d_N1, d_N2, d_N3, d_Nxyz.z, bmin, bmax, Vmax, V, XY, ax, 
#if defined(FP) //////////////// FP ////////////////
                            d_OSEM
#else //////////////// BP ////////////////
                            no_norm, d_Summ, d_output 
#endif //////////////// END FP/BP ////////////////
#ifdef TOF //////////////// TOF ////////////////
                        , local_ele2, sigma_x, &D, DD, TOFCenter, TOFSum
#if defined(LISTMODE)
						, TOFid
#endif
#endif //////////////// END TOF ////////////////
#if defined(SPECT) && defined(ATN)
						, EXP(jelppi)
#endif
#if defined(MASKBP) && defined(BP) //////////////// MASKBP ////////////////
						, aa, maskBP
#endif //////////////// END MASKBP ////////////////
                        );
                        if (uu == 0)
                            break;
                    }
                }
            }
			if (tz0_a >= tx0_a && ty0_a >= tx0_a) {
				const float xcenter = b1 + d1 * CFLOAT(localInd.x) + d1 / 2.f;
				orthDistance3D(localInd.x, diff.y, diff.x, diff.z, xcenter, b2, d2, _bz, dz, temp, localInd.y, localInd.z, s.x, s.y, s.z, d_Nxy, kerroin, d_N1, d_N2, d_N3, d_Nxyz.z, bmin, bmax, Vmax, V, XY, ax, 
#if defined(FP) //////////////// FP ////////////////
					d_OSEM
#else //////////////// BP ////////////////
					no_norm, d_Summ, d_output 
#endif //////////////// END FP/BP ////////////////
#ifdef TOF //////////////// TOF ////////////////
				, local_ele, sigma_x, &D, DD, TOFCenter, TOFSum
#if defined(LISTMODE)
				, TOFid
#endif
#endif //////////////// END TOF ////////////////
#if defined(SPECT) && defined(ATN)
				, EXP(jelppi)
#endif
#if defined(MASKBP) && defined(BP) //////////////// MASKBP ////////////////
				, aa, maskBP
#endif //////////////// END MASKBP ////////////////
				);
			}
#else //////////////// SIDDON ////////////////
#if defined(FP) //////////////// FORWARD PROJECTION ////////////////
#ifdef USEIMAGES
			denominator(ax, localInd, local_ele, d_OSEM
#else
			denominator(ax, local_ind, local_ele, d_OSEM
#endif
#ifdef TOF //////////////// TOF ////////////////
			, local_ele, TOFSum, DD, TOFCenter, sigma_x, &D
#ifdef LISTMODE
			, TOFid
#endif
#endif //////////////// END TOF ////////////////
#ifdef N_RAYS
			, lor
#endif
			);
#endif  //////////////// END FORWARD PROJECTION ////////////////
#if defined(BP) //////////////// BACKWARD PROJECTION ////////////////
#if defined(MASKBP) //////////////// MASKBP ////////////////
			int maskVal = 1;
			if (aa == 0) {
#ifdef USEIMAGES
#ifdef CUDA
#ifdef MASKBP3D
				maskVal = tex3D<unsigned char>(maskBP, localInd.x, localInd.y, localInd.z);
#else
				maskVal = tex2D<unsigned char>(maskBP, localInd.x, localInd.y);
#endif
#else
#ifdef MASKBP3D
				maskVal = read_imageui(maskBP, sampler_MASK, (int4)(localInd.x, localInd.y,  localInd.z, 0)).w;
#else
				maskVal = read_imageui(maskBP, sampler_MASK, (int2)(localInd.x, localInd.y)).w;
#endif
#endif
#else
#ifdef MASKBP3D
				maskVal = maskBP[localInd.x * d_N2 + localInd.y * d_N3 + localInd.z * d_Nxyz.x * d_Nxyz.y];
#else
				maskVal = maskBP[localInd.x * d_N2 + localInd.y * d_N3];
#endif
#endif
			}
			if (maskVal > 0)
#endif //////////////// END MASKBP ////////////////
			rhs(local_ele * temp, ax, local_ind, d_output, no_norm, d_Summ
#ifdef TOF
			, local_ele, sigma_x, &D, DD, TOFCenter, TOFSum
#ifdef LISTMODE
			, TOFid
#endif
#endif
			);
#endif //////////////// END BACKWARD PROJECTION ////////////////
#endif //////////////// END ORTH/VOL ////////////////
#if defined(TOF)
			D -= (local_ele2 * sign(DD));
#endif
			if (tempi < 0 || tempi >= d_Nxyz.x || tempj < 0 || tempj >= d_Nxyz.y || tempk < 0 || tempk >= d_Nxyz.z) {
				break;
			}
// #if (defined(SPECT) && !defined(ORTH)) // Ray length inside BP mask
// 		L_SPECT += local_ele;
// #endif
		}

// #if (defined(SPECT) && !defined(ORTH)) // Ray length inside BP mask
// 		temp /= L_SPECT;
// #endif
#if !defined(TOTLENGTH) && !defined(CT) && defined(FP)
			if (LL == 0.f)
				LL = L;
#if defined(N_RAYS) && !defined(ORTH)
			temp = 1.f / (LL * CFLOAT(N_RAYS));
#elif !defined(ORTH)
			temp = 1.f / LL;
#endif
#if defined(VOL)
			temp = 1.f / (TotV * LL);
#endif
#if defined(ATN) && defined(BP) && !defined(SPECT)
			temp *= EXP(jelppi);
#endif
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

#ifdef ORTH
		if (ux < 0) {
			for (int ii = tempi_a; ii >= 0; ii--) {
				const float xcenter = b1 + d1 * CFLOAT(ii) + d1 / 2.f;
				int uu = orthDistance3D(ii, diff.y, diff.x, diff.z, xcenter, b2, d2, _bz, dz, temp, tempj_a, tempk_a, s.x, s.y, s.z, d_Nxy, kerroin, d_N1, d_N2, d_N3, d_Nxyz.z, bmin, bmax, Vmax, V, XY, ax, 
#if defined(FP)
					d_OSEM
#else
					no_norm, d_Summ, d_output 
#endif
#ifdef TOF
				, local_ele, sigma_x, &D, DD, TOFCenter, TOFSum
#ifdef LISTMODE
				, TOFid
#endif
#endif
#if defined(SPECT) && defined(ATN)
				, EXP(jelppi)
#endif
#if defined(MASKBP) && defined(BP)
				, aa, maskBP
#endif
				);
				if (uu == 0)
					break;
			}
		}
		else {
			for (int ii = tempi_a; ii < d_NN.x; ii++) {
				const float xcenter = b1 + d1 * CFLOAT(ii) + d1 / 2.f;
				int uu = orthDistance3D(ii, diff.y, diff.x, diff.z, xcenter, b2, d2, _bz, dz, temp, tempj_a, tempk_a, s.x, s.y, s.z, d_Nxy, kerroin, d_N1, d_N2, d_N3, d_Nxyz.z, bmin, bmax, Vmax, V, XY, ax, 
#if defined(FP)
					d_OSEM
#else
					no_norm, d_Summ, d_output 
#endif
#ifdef TOF
				, local_ele, sigma_x, &D, DD, TOFCenter, TOFSum
#ifdef LISTMODE
				, TOFid
#endif
#endif
#if defined(SPECT) && defined(ATN)
				, EXP(jelppi)
#endif
#if defined(MASKBP) && defined(BP)
				, aa, maskBP
#endif
				);
				if (uu == 0)
					break;
			}
		}
#endif
#if defined(ATN) && defined(FP) && !defined(SPECT)
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
	}
//  */
#ifdef N_RAYS //////////////// MULTIRAY ////////////////
		}
	}


#if defined(FP) //////////////// FORWARD PROJECTION ////////////////
#ifndef __CUDACC__ 
#pragma unroll NBINS
#endif
    for (size_t to = 0; to < NBINS; to++) {
        float apu = 0.f;
#pragma unroll N_RAYS
		for (size_t kk = 0; kk < N_RAYS; kk++) {
            apu += ax[to + NBINS * kk];
        }
        ax[to] = apu;
    }
#ifndef __CUDACC__ 
#pragma unroll NBINS
#endif
    for (size_t to = 0; to < NBINS; to++) {
        forwardProjectAF(d_output, ax, idx, 1.f, to);
#ifdef TOF
        idx += m_size;
#endif
    }
#endif //////////////// END FORWARD PROJECTION ////////////////
#endif //////////////// END MULTIRAY ////////////////
}

